module derivatives
use mpiinfo
use declaration
implicit none


contains
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0

selectcase (nderivative)
case(1)  
                    
s=1.0d0

case(4)  
                    
s=2.0d0*xd1
case(5)  
                    
s=yd1

case(7)  
                    
s=zd1


case(10)  
                    
s=3.d0*xd1**2
case(11)  
                    
s=2.d0*xd1*yd1
case(12)  
                    
s=yd1**2

case(14)  
                    
s=2.d0*xd1*zd1
case(15)  
                    
s=yd1*zd1

case(17)  
                    
s=zd1**2


case(20)  
                    
s=4.d0*xd1**3
case(21)  
                    
s=3.d0*xd1**2*yd1
case(22)  
                    
s=2.d0*xd1*yd1**2
case(23)  
                    
s=yd1**3

case(25)  
                    
s=3.d0*xd1**2*zd1
case(26)  
                    
s=2.d0*xd1*yd1*zd1
case(27)  
                    
s=yd1**2*zd1

case(29)  
                    
s=2.d0*xd1*zd1**2
case(30)  
                    
s=yd1*zd1**2

case(32)  
                    
s=zd1**3


case(35)  
                    
s=5.d0*xd1**4
case(36)  
                    
s=4.d0*xd1**3*yd1
case(37)  
                    
s=3.d0*xd1**2*yd1**2
case(38)  
                    
s=2.d0*xd1*yd1**3
case(39)  
                    
s=yd1**4

case(41)  
                    
s=4.d0*xd1**3*zd1
case(42)  
                    
s=3.d0*xd1**2*yd1*zd1
case(43)  
                    
s=2.d0*xd1*yd1**2*zd1
case(44)  
                    
s=yd1**3*zd1

case(46)  
                    
s=3.d0*xd1**2*zd1**2
case(47)  
                    
s=2.d0*xd1*yd1*zd1**2
case(48)  
                    
s=yd1**2*zd1**2

case(50)  
                    
s=2.d0*xd1*zd1**3
case(51)  
                    
s=yd1*zd1**3

case(53)  
                    
s=zd1**4

case(56)  
                    
s=6.d0*xd1**5
case(57)  
                    
s=5.d0*xd1**4*yd1
case(58)  
                    
s=4.d0*xd1**3*yd1**2
case(59)  
                    
s=3.d0*xd1**2*yd1**3
case(60)  
                    
s=2.d0*xd1*yd1**4
case(61)  
                    
s=yd1**5

case(63)  
                    
s=5.d0*xd1**4*zd1
case(64)  
                    
s=4.d0*xd1**3*yd1*zd1
case(65)  
                    
s=3.d0*xd1**2*yd1**2*zd1
case(66)  
                    
s=2.d0*xd1*yd1**3*zd1
case(67)  
                    
s=yd1**4*zd1

case(69)  
                    
s=4.d0*xd1**3*zd1**2
case(70)  
                    
s=3.d0*xd1**2*yd1*zd1**2
case(71)  
                    
s=2.d0*xd1*yd1**2*zd1**2
case(72)  
                    
s=yd1**3*zd1**2

case(74)  
                    
s=3.d0*xd1**2*zd1**3
case(75)  
                    
s=2.d0*xd1*yd1*zd1**3
case(76)  
                    
s=yd1**2*zd1**3

case(78)  
                    
s=2.d0*xd1*zd1**4
case(79)  
                    
s=yd1*zd1**4

case(81)  
                    
s=zd1**5

case default

s=0.0

end select
dfx=s
end function dfx
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)

case(2) 
                    
s=1d0

case(5) 
                    
s=xd1
case(6) 
                    
s=2.d0*yd1

case(8) 
                    
s=zd1


case(11) 
                    
s=xd1**2
case(12) 
                    
s=2.d0*xd1*yd1

case(13)
s=3.0d0*yd1**2

case(15) 
                    
s=xd1*zd1
case(16) 
                    
s=2.d0*yd1*zd1

case(18) 
                    
s=zd1**2

case(21) 
                    
s=xd1**3
case(22) 
                    
s=2.d0*xd1**2*yd1
case(23) 
                    
s=3.d0*xd1*yd1**2
case(24) 
                    
s=4.d0*yd1**3

case(26) 
                    
s=xd1**2*zd1
case(27) 
                    
s=2.d0*xd1*yd1*zd1
case(28) 
                    
s=3.d0*yd1**2*zd1

case(30) 
                    
s=xd1*zd1**2
case(31) 
                    
s=2.d0*yd1*zd1**2

case(33) 
                    
s=zd1**3

case(36) 
                    
s=xd1**4
case(37) 
                    
s=2.d0*xd1**3*yd1
case(38) 
                    
s=3.d0*xd1**2*yd1**2
case(39) 
                    
s=4.d0*xd1*yd1**3
case(40) 
                    
s=5.d0*yd1**4

case(42) 
                    
s=xd1**3*zd1
case(43) 
                    
s=2.d0*xd1**2*yd1*zd1
case(44) 
                    
s=3.d0*xd1*yd1**2*zd1
case(45) 
                    
s=4.d0*yd1**3*zd1

case(47) 
                    
s=xd1**2*zd1**2
case(48) 
                    
s=2.d0*xd1*yd1*zd1**2
case(49) 
                    
s=3.d0*yd1**2*zd1**2

case(51) 
                    
s=xd1*zd1**3
case(52) 
                    
s=2.d0*yd1*zd1**3

case(54) 
                    
s=zd1**4

case(57) 
                    
s=xd1**5
case(58) 
                    
s=2.d0*xd1**4*yd1
case(59) 
                    
s=3.d0*xd1**3*yd1**2
case(60) 
                    
s=4.d0*xd1**2*yd1**3
case(61) 
                    
s=5.d0*xd1*yd1**4
case(62) 
                    
s=6.d0*yd1**5

case(64) 
                    
s=xd1**4*zd1
case(65) 
                    
s=2.d0*xd1**3*yd1*zd1
case(66) 
                    
s=3.d0*xd1**2*yd1**2*zd1
case(67) 
                    
s=4.d0*xd1*yd1**3*zd1
case(68) 
                    
s=5.d0*yd1**4*zd1

case(70) 
                    
s=xd1**3*zd1**2
case(71) 
                    
s=2.d0*xd1**2*yd1*zd1**2
case(72) 
                    
s=3.d0*xd1*yd1**2*zd1**2
case(73) 
                    
s=4.d0*yd1**3*zd1**2

case(75) 
                    
s=xd1**2*zd1**3
case(76) 
                    
s=2.d0*xd1*yd1*zd1**3
case(77) 
                    
s=3.d0*yd1**2*zd1**3

case(79) 
                    
s=xd1*zd1**4
case(80) 
                    
s=2.d0*yd1*zd1**4

case(82) 
                    
s=zd1**5


case default
s=0.0d0



end select
dfy=s
end function dfy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)

case (3)  
                    
s=1

case (7)  
                    
s=xd1
case (8)  
                    
s=yd1
case (9)  
                    
s=2.d0*zd1


case (14)  
                    
s=xd1**2
case (15)  
                    
s=xd1*yd1
case (16)  
                    
s=yd1**2
case (17)  
                    
s=2.d0*xd1*zd1
case (18)  
                    
s=2.d0*yd1*zd1
case (19)  
                    
s=3.d0*zd1**2

case (25)  
                    
s=xd1**3
case (26)  
                    
s=xd1**2*yd1
case (27)  
                    
s=xd1*yd1**2
case (28)  
                    
s=yd1**3
case (29)  
                    
s=2.d0*xd1**2*zd1
case (30)  
                    
s=2.d0*xd1*yd1*zd1
case (31)  
                    
s=2.d0*yd1**2*zd1
case (32)  
                    
s=3.d0*xd1*zd1**2
case (33)  
                    
s=3.d0*yd1*zd1**2
case (34)  
                    
s=4.d0*zd1**3

                    
s=xd1**4
case (42)  
                    
s=xd1**3*yd1
case (43)  
                    
s=xd1**2*yd1**2
case (44)  
                    
s=xd1*yd1**3
case (45)  
                    
s=yd1**4
case (46)  
                    
s=2.d0*xd1**3*zd1
case (47)  
                    
s=2.d0*xd1**2*yd1*zd1
case (48)  
                    
s=2.d0*xd1*yd1**2*zd1
case (49)  
                    
s=2.d0*yd1**3*zd1
case (50)  
                    
s=3.d0*xd1**2*zd1**2
case (51)  
                    
s=3.d0*xd1*yd1*zd1**2
case (52)  
                    
s=3.d0*yd1**2*zd1**2
case (53)  
                    
s=4.d0*xd1*zd1**3
case (54)  
                    
s=4.d0*yd1*zd1**3
case (55)  
                    
s=5.d0*zd1**4

                    
s=xd1**5
case (64)  
                    
s=xd1**4*yd1
case (65)  
                    
s=xd1**3*yd1**2
case (66)  
                    
s=xd1**2*yd1**3
case (67)  
                    
s=xd1*yd1**4
case (68)  
                    
s=yd1**5
case (69)  
                    
s=2.d0*xd1**4*zd1
case (70)  
                    
s=2.d0*xd1**3*yd1*zd1
case (71)  
                    
s=2.d0*xd1**2*yd1**2*zd1
case (72)  
                    
s=2.d0*xd1*yd1**3*zd1
case (73)  
                    
s=2.d0*yd1**4*zd1
case (74)  
                    
s=3.d0*xd1**3*zd1**2
case (75)  
                    
s=3.d0*xd1**2*yd1*zd1**2
case (76)  
                    
s=3.d0*xd1*yd1**2*zd1**2
case (77)  
                    
s=3.d0*yd1**3*zd1**2
case (78)  
                    
s=4.d0*xd1**2*zd1**3
case (79)  
                    
s=4.d0*xd1*yd1*zd1**3
case (80)  
                    
s=4.d0*yd1**2*zd1**3
case (81)  
                    
s=5.d0*xd1*zd1**4
case (82)  
                    
s=5.d0*yd1*zd1**4
case (83)  
                    
s=6.d0*zd1**5

case default
s=0.0d0

end select
dfz=s
end function dfz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=2
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=6.d0*xd1
case(11) 
                    
s=2.d0*yd1
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=2.d0*zd1
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=12.d0*xd1**2
case(21) 
                    
s=6.d0*xd1*yd1
case(22) 
                    
s=2.d0*yd1**2
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=6.d0*xd1*zd1
case(26) 
                    
s=2.d0*yd1*zd1
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=2.d0*zd1**2
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=20*xd1**3
case(36) 
                    
s=12.d0*xd1**2*yd1
case(37) 
                    
s=6.d0*xd1*yd1**2
case(38) 
                    
s=2.d0*yd1**3
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=12.d0*xd1**2*zd1
case(42) 
                    
s=6.d0*xd1*yd1*zd1
case(43) 
                    
s=2.d0*yd1**2*zd1
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=6.d0*xd1*zd1**2
case(47) 
                    
s=2.d0*yd1*zd1**2
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=2.d0*zd1**3
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=30*xd1**4
case(57) 
                    
s=20*xd1**3*yd1
case(58) 
                    
s=12.d0*xd1**2*yd1**2
case(59) 
                    
s=6.d0*xd1*yd1**3
case(60) 
                    
s=2.d0*yd1**4
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=20*xd1**3*zd1
case(64) 
                    
s=12.d0*xd1**2*yd1*zd1
case(65) 
                    
s=6.d0*xd1*yd1**2*zd1
case(66) 
                    
s=2.d0*yd1**3*zd1
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=12.d0*xd1**2*zd1**2
case(70) 
                    
s=6.d0*xd1*yd1*zd1**2
case(71) 
                    
s=2.d0*yd1**2*zd1**2
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=6.d0*xd1*zd1**3
case(75) 
                    
s=2.d0*yd1*zd1**3
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=2.d0*zd1**4
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2=s
end function dfx2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=2
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=2.d0*xd1
case(13) 
                    
s=6.d0*yd1
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=2.d0*zd1
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=2.d0*xd1**2
case(23) 
                    
s=6.d0*xd1*yd1
case(24) 
                    
s=12.d0*yd1**2
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=2.d0*xd1*zd1
case(28) 
                    
s=6.d0*yd1*zd1
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=2.d0*zd1**2
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=2.d0*xd1**3
case(38) 
                    
s=6.d0*xd1**2*yd1
case(39) 
                    
s=12.d0*xd1*yd1**2
case(40) 
                    
s=20*yd1**3
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=2.d0*xd1**2*zd1
case(44) 
                    
s=6.d0*xd1*yd1*zd1
case(45) 
                    
s=12.d0*yd1**2*zd1
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=2.d0*xd1*zd1**2
case(49) 
                    
s=6.d0*yd1*zd1**2
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=2.d0*zd1**3
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=2.d0*xd1**4
case(59) 
                    
s=6.d0*xd1**3*yd1
case(60) 
                    
s=12.d0*xd1**2*yd1**2
case(61) 
                    
s=20*xd1*yd1**3
case(62) 
                    
s=30*yd1**4
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=2.d0*xd1**3*zd1
case(66) 
                    
s=6.d0*xd1**2*yd1*zd1
case(67) 
                    
s=12.d0*xd1*yd1**2*zd1
case(68) 
                    
s=20*yd1**3*zd1
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=2.d0*xd1**2*zd1**2
case(72) 
                    
s=6.d0*xd1*yd1*zd1**2
case(73) 
                    
s=12.d0*yd1**2*zd1**2
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=2.d0*xd1*zd1**3
case(77) 
                    
s=6.d0*yd1*zd1**3
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=2.d0*zd1**4
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy2=s
end function dfy2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1)  
                    
s=0.0d0
case (2)  
                    
s=0.0d0
case (3)  
                    
s=0.0d0
case (4)  
                    
s=0.0d0
case (5)  
                    
s=0.0d0
case (6)  
                    
s=0.0d0
case (7)  
                    
s=0.0d0
case (8)  
                    
s=0.0d0
case (9)  
                    
s=2
case (10)  
                    
s=0.0d0
case (11)  
                    
s=0.0d0
case (12)  
                    
s=0.0d0
case (13)  
                    
s=0.0d0
case (14)  
                    
s=0.0d0
case (15)  
                    
s=0.0d0
case (16)  
                    
s=0.0d0
case (17)  
                    
s=2.d0*xd1
case (18)  
                    
s=2.d0*yd1
case (19)  
                    
s=6.d0*zd1
case (20)  
                    
s=0.0d0
case (21)  
                    
s=0.0d0
case (22)  
                    
s=0.0d0
case (23)  
                    
s=0.0d0
case (24)  
                    
s=0.0d0
case (25)  
                    
s=0.0d0
case (26)  
                    
s=0.0d0
case (27)  
                    
s=0.0d0
case (28)  
                    
s=0.0d0
case (29)  
                    
s=2.d0*xd1**2
case (30)  
                    
s=2.d0*xd1*yd1
case (31)  
                    
s=2.d0*yd1**2
case (32)  
                    
s=6.d0*xd1*zd1
case (33)  
                    
s=6.d0*yd1*zd1
case (34)  
                    
s=12.d0*zd1**2
case (35)  
                    
s=0.0d0
case (36)  
                    
s=0.0d0
case (37)  
                    
s=0.0d0
case (38)  
                    
s=0.0d0
case (39)  
                    
s=0.0d0
case (40)  
                    
s=0.0d0
case (41)  
                    
s=0.0d0
case (42)  
                    
s=0.0d0
case (43)  
                    
s=0.0d0
case (44)  
                    
s=0.0d0
case (45)  
                    
s=0.0d0
case (46)  
                    
s=2.d0*xd1**3
case (47)  
                    
s=2.d0*xd1**2*yd1
case (48)  
                    
s=2.d0*xd1*yd1**2
case (49)  
                    
s=2.d0*yd1**3
case (50)  
                    
s=6.d0*xd1**2*zd1
case (51)  
                    
s=6.d0*xd1*yd1*zd1
case (52)  
                    
s=6.d0*yd1**2*zd1
case (53)  
                    
s=12.d0*xd1*zd1**2
case (54)  
                    
s=12.d0*yd1*zd1**2
case (55)  
                    
s=20*zd1**3
case (56)  
                    
s=0.0d0
case (57)  
                    
s=0.0d0
case (58)  
                    
s=0.0d0
case (59)  
                    
s=0.0d0
case (60)  
                    
s=0.0d0
case (61)  
                    
s=0.0d0
case (62)  
                    
s=0.0d0
case (63)  
                    
s=0.0d0
case (64)  
                    
s=0.0d0
case (65)  
                    
s=0.0d0
case (66)  
                    
s=0.0d0
case (67)  
                    
s=0.0d0
case (68)  
                    
s=0.0d0
case (69)  
                    
s=2.d0*xd1**4
case (70)  
                    
s=2.d0*xd1**3*yd1
case (71)  
                    
s=2.d0*xd1**2*yd1**2
case (72)  
                    
s=2.d0*xd1*yd1**3
case (73)  
                    
s=2.d0*yd1**4
case (74)  
                    
s=6.d0*xd1**3*zd1
case (75)  
                    
s=6.d0*xd1**2*yd1*zd1
case (76)  
                    
s=6.d0*xd1*yd1**2*zd1
case (77)  
                    
s=6.d0*yd1**3*zd1
case (78)  
                    
s=12.d0*xd1**2*zd1**2
case (79)  
                    
s=12.d0*xd1*yd1*zd1**2
case (80)  
                    
s=12.d0*yd1**2*zd1**2
case (81)  
                    
s=20*xd1*zd1**3
case (82)  
                    
s=20*yd1*zd1**3
case (83)  
                    
s=30*zd1**4
end select
dfz2=s
end function dfz2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=1
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=2.d0*xd1
case(12) 
                    
s=2.d0*yd1
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=zd1
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=3.d0*xd1**2
case(22) 
                    
s=4.d0*xd1*yd1
case(23) 
                    
s=3.d0*yd1**2
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=2.d0*xd1*zd1
case(27) 
                    
s=2.d0*yd1*zd1
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=zd1**2
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=4.d0*xd1**3
case(37) 
                    
s=6.d0*xd1**2*yd1
case(38) 
                    
s=6.d0*xd1*yd1**2
case(39) 
                    
s=4.d0*yd1**3
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=3.d0*xd1**2*zd1
case(43) 
                    
s=4.d0*xd1*yd1*zd1
case(44) 
                    
s=3.d0*yd1**2*zd1
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=2.d0*xd1*zd1**2
case(48) 
                    
s=2.d0*yd1*zd1**2
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=zd1**3
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=5.d0*xd1**4
case(58) 
                    
s=8.d0*xd1**3*yd1
case(59) 
                    
s=9*xd1**2*yd1**2
case(60) 
                    
s=8.d0*xd1*yd1**3
case(61) 
                    
s=5.d0*yd1**4
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=4.d0*xd1**3*zd1
case(65) 
                    
s=6.d0*xd1**2*yd1*zd1
case(66) 
                    
s=6.d0*xd1*yd1**2*zd1
case(67) 
                    
s=4.d0*yd1**3*zd1
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=3.d0*xd1**2*zd1**2
case(71) 
                    
s=4.d0*xd1*yd1*zd1**2
case(72) 
                    
s=3.d0*yd1**2*zd1**2
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=2.d0*xd1*zd1**3
case(76) 
                    
s=2.d0*yd1*zd1**3
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=zd1**4
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy=s
end function dfxy
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfyz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1)  
                    
s=0.0d0
case(2)  
                    
s=0.0d0
case(3)  
                    
s=0.0d0
case(4)  
                    
s=0.0d0
case(5)  
                    
s=0.0d0
case(6)  
                    
s=0.0d0
case(7)  
                    
s=0.0d0
case(8)  
                    
s=1
case(9)  
                    
s=0.0d0
case(10)  
                    
s=0.0d0
case(11)  
                    
s=0.0d0
case(12)  
                    
s=0.0d0
case(13)  
                    
s=0.0d0
case(14)  
                    
s=0.0d0
case(15)  
                    
s=xd1
case(16)  
                    
s=2.d0*yd1
case(17)  
                    
s=0.0d0
case(18)  
                    
s=2.d0*zd1
case(19)  
                    
s=0.0d0
case(20)  
                    
s=0.0d0
case(21)  
                    
s=0.0d0
case(22)  
                    
s=0.0d0
case(23)  
                    
s=0.0d0
case(24)  
                    
s=0.0d0
case(25)  
                    
s=0.0d0
case(26)  
                    
s=xd1**2
case(27)  
                    
s=2.d0*xd1*yd1
case(28)  
                    
s=3.d0*yd1**2
case(29)  
                    
s=0.0d0
case(30)  
                    
s=2.d0*xd1*zd1
case(31)  
                    
s=4.d0*yd1*zd1
case(32)  
                    
s=0.0d0
case(33)  
                    
s=3.d0*zd1**2
case(34)  
                    
s=0.0d0
case(35)  
                    
s=0.0d0
case(36)  
                    
s=0.0d0
case(37)  
                    
s=0.0d0
case(38)  
                    
s=0.0d0
case(39)  
                    
s=0.0d0
case(40)  
                    
s=0.0d0
case(41)  
                    
s=0.0d0
case(42)  
                    
s=xd1**3
case(43)  
                    
s=2.d0*xd1**2*yd1
case(44)  
                    
s=3.d0*xd1*yd1**2
case(45)  
                    
s=4.d0*yd1**3
case(46)  
                    
s=0.0d0
case(47)  
                    
s=2.d0*xd1**2*zd1
case(48)  
                    
s=4.d0*xd1*yd1*zd1
case(49)  
                    
s=6.d0*yd1**2*zd1
case(50)  
                    
s=0.0d0
case(51)  
                    
s=3.d0*xd1*zd1**2
case(52)  
                    
s=6.d0*yd1*zd1**2
case(53)  
                    
s=0.0d0
case(54)  
                    
s=4.d0*zd1**3
case(55)  
                    
s=0.0d0
case(56)  
                    
s=0.0d0
case(57)  
                    
s=0.0d0
case(58)  
                    
s=0.0d0
case(59)  
                    
s=0.0d0
case(60)  
                    
s=0.0d0
case(61)  
                    
s=0.0d0
case(62)  
                    
s=0.0d0
case(63)  
                    
s=0.0d0
case(64)  
                    
s=xd1**4
case(65)  
                    
s=2.d0*xd1**3*yd1
case(66)  
                    
s=3.d0*xd1**2*yd1**2
case(67)  
                    
s=4.d0*xd1*yd1**3
case(68)  
                    
s=5.d0*yd1**4
case(69)  
                    
s=0.0d0
case(70)  
                    
s=2.d0*xd1**3*zd1
case(71)  
                    
s=4.d0*xd1**2*yd1*zd1
case(72)  
                    
s=6.d0*xd1*yd1**2*zd1
case(73)  
                    
s=8.d0*yd1**3*zd1
case(74)  
                    
s=0.0d0
case(75)  
                    
s=3.d0*xd1**2*zd1**2
case(76)  
                    
s=6.d0*xd1*yd1*zd1**2
case(77)  
                    
s=9*yd1**2*zd1**2
case(78)  
                    
s=0.0d0
case(79)  
                    
s=4.d0*xd1*zd1**3
case(80)  
                    
s=8.d0*yd1*zd1**3
case(81)  
                    
s=0.0d0
case(82)  
                    
s=5.d0*zd1**4
case(83)  
                    
s=0.0d0
end select
dfyz=s
end function dfyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=1
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=2.d0*xd1
case(15) 
                    
s=yd1
case(16) 
                    
s=0.0d0
case(17) 
                    
s=2.d0*zd1
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=3.d0*xd1**2
case(26) 
                    
s=2.d0*xd1*yd1
case(27) 
                    
s=yd1**2
case(28) 
                    
s=0.0d0
case(29) 
                    
s=4.d0*xd1*zd1
case(30) 
                    
s=2.d0*yd1*zd1
case(31) 
                    
s=0.0d0
case(32) 
                    
s=3.d0*zd1**2
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=4.d0*xd1**3
case(42) 
                    
s=3.d0*xd1**2*yd1
case(43) 
                    
s=2.d0*xd1*yd1**2
case(44) 
                    
s=yd1**3
case(45) 
                    
s=0.0d0
case(46) 
                    
s=6.d0*xd1**2*zd1
case(47) 
                    
s=4.d0*xd1*yd1*zd1
case(48) 
                    
s=2.d0*yd1**2*zd1
case(49) 
                    
s=0.0d0
case(50) 
                    
s=6.d0*xd1*zd1**2
case(51) 
                    
s=3.d0*yd1*zd1**2
case(52) 
                    
s=0.0d0
case(53) 
                    
s=4.d0*zd1**3
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=5.d0*xd1**4
case(64) 
                    
s=4.d0*xd1**3*yd1
case(65) 
                    
s=3.d0*xd1**2*yd1**2
case(66) 
                    
s=2.d0*xd1*yd1**3
case(67) 
                    
s=yd1**4
case(68) 
                    
s=0.0d0
case(69) 
                    
s=8.d0*xd1**3*zd1
case(70) 
                    
s=6.d0*xd1**2*yd1*zd1
case(71) 
                    
s=4.d0*xd1*yd1**2*zd1
case(72) 
                    
s=2.d0*yd1**3*zd1
case(73) 
                    
s=0.0d0
case(74) 
                    
s=9*xd1**2*zd1**2
case(75) 
                    
s=6.d0*xd1*yd1*zd1**2
case(76) 
                    
s=3.d0*yd1**2*zd1**2
case(77) 
                    
s=0.0d0
case(78) 
                    
s=8.d0*xd1*zd1**3
case(79) 
                    
s=4.d0*yd1*zd1**3
case(80) 
                    
s=0.0d0
case(81) 
                    
s=5.d0*zd1**4
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxz=s
end function dfxz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=6
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=24.d0*xd1
case(21) 
                    
s=6.d0*yd1
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=6.d0*zd1
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=60*xd1**2
case(36) 
                    
s=24.d0*xd1*yd1
case(37) 
                    
s=6.d0*yd1**2
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=24.d0*xd1*zd1
case(42) 
                    
s=6.d0*yd1*zd1
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=6.d0*zd1**2
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=120*xd1**3
case(57) 
                    
s=60*xd1**2*yd1
case(58) 
                    
s=24.d0*xd1*yd1**2
case(59) 
                    
s=6.d0*yd1**3
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=60*xd1**2*zd1
case(64) 
                    
s=24.d0*xd1*yd1*zd1
case(65) 
                    
s=6.d0*yd1**2*zd1
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=24.d0*xd1*zd1**2
case(70) 
                    
s=6.d0*yd1*zd1**2
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=6.d0*zd1**3
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3=s
end function dfx3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=2
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=6.d0*xd1
case(22) 
                    
s=4.d0*yd1
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=2.d0*zd1
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=12.d0*xd1**2
case(37) 
                    
s=12.d0*xd1*yd1
case(38) 
                    
s=6.d0*yd1**2
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=6.d0*xd1*zd1
case(43) 
                    
s=4.d0*yd1*zd1
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=2.d0*zd1**2
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=20*xd1**3
case(58) 
                    
s=24.d0*xd1**2*yd1
case(59) 
                    
s=18.d0*xd1*yd1**2
case(60) 
                    
s=8.d0*yd1**3
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=12.d0*xd1**2*zd1
case(65) 
                    
s=12.d0*xd1*yd1*zd1
case(66) 
                    
s=6.d0*yd1**2*zd1
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=6.d0*xd1*zd1**2
case(71) 
                    
s=4.d0*yd1*zd1**2
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=2.d0*zd1**3
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2y=s
end function dfx2y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=2
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=4.d0*xd1
case(23) 
                    
s=6.d0*yd1
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=2.d0*zd1
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=6.d0*xd1**2
case(38) 
                    
s=12.d0*xd1*yd1
case(39) 
                    
s=12.d0*yd1**2
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=4.d0*xd1*zd1
case(44) 
                    
s=6.d0*yd1*zd1
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=2.d0*zd1**2
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=8.d0*xd1**3
case(59) 
                    
s=18.d0*xd1**2*yd1
case(60) 
                    
s=24.d0*xd1*yd1**2
case(61) 
                    
s=20*yd1**3
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=6.d0*xd1**2*zd1
case(66) 
                    
s=12.d0*xd1*yd1*zd1
case(67) 
                    
s=12.d0*yd1**2*zd1
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=4.d0*xd1*zd1**2
case(72) 
                    
s=6.d0*yd1*zd1**2
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=2.d0*zd1**3
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy2=s
end function dfxy2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=6
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=6.d0*xd1
case(24) 
                    
s=24.d0*yd1
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=6.d0*zd1
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=6.d0*xd1**2
case(39) 
                    
s=24.d0*xd1*yd1
case(40) 
                    
s=60*yd1**2
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=6.d0*xd1*zd1
case(45) 
                    
s=24.d0*yd1*zd1
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=6.d0*zd1**2
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=6.d0*xd1**3
case(60) 
                    
s=24.d0*xd1**2*yd1
case(61) 
                    
s=60*xd1*yd1**2
case(62) 
                    
s=120*yd1**3
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=6.d0*xd1**2*zd1
case(67) 
                    
s=24.d0*xd1*yd1*zd1
case(68) 
                    
s=60*yd1**2*zd1
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=6.d0*xd1*zd1**2
case(73) 
                    
s=24.d0*yd1*zd1**2
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=6.d0*zd1**3
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy3=s
end function dfy3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=2
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=6.d0*xd1
case(26) 
                    
s=2.d0*yd1
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=4.d0*zd1
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=12.d0*xd1**2
case(42) 
                    
s=6.d0*xd1*yd1
case(43) 
                    
s=2.d0*yd1**2
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=12.d0*xd1*zd1
case(47) 
                    
s=4.d0*yd1*zd1
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=6.d0*zd1**2
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=20*xd1**3
case(64) 
                    
s=12.d0*xd1**2*yd1
case(65) 
                    
s=6.d0*xd1*yd1**2
case(66) 
                    
s=2.d0*yd1**3
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=24.d0*xd1**2*zd1
case(70) 
                    
s=12.d0*xd1*yd1*zd1
case(71) 
                    
s=4.d0*yd1**2*zd1
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=18.d0*xd1*zd1**2
case(75) 
                    
s=6.d0*yd1*zd1**2
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=8.d0*zd1**3
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2z=s
end function dfx2z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxyz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=1
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=2.d0*xd1
case(27) 
                    
s=2.d0*yd1
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=2.d0*zd1
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=3.d0*xd1**2
case(43) 
                    
s=4.d0*xd1*yd1
case(44) 
                    
s=3.d0*yd1**2
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=4.d0*xd1*zd1
case(48) 
                    
s=4.d0*yd1*zd1
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=3.d0*zd1**2
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=4.d0*xd1**3
case(65) 
                    
s=6.d0*xd1**2*yd1
case(66) 
                    
s=6.d0*xd1*yd1**2
case(67) 
                    
s=4.d0*yd1**3
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=6.d0*xd1**2*zd1
case(71) 
                    
s=8.d0*xd1*yd1*zd1
case(72) 
                    
s=6.d0*yd1**2*zd1
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=6.d0*xd1*zd1**2
case(76) 
                    
s=6.d0*yd1*zd1**2
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=4.d0*zd1**3
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxyz=s
end function dfxyz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=2
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=2.d0*xd1
case(28) 
                    
s=6.d0*yd1
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=4.d0*zd1
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=2.d0*xd1**2
case(44) 
                    
s=6.d0*xd1*yd1
case(45) 
                    
s=12.d0*yd1**2
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=4.d0*xd1*zd1
case(49) 
                    
s=12.d0*yd1*zd1
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=6.d0*zd1**2
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=2.d0*xd1**3
case(66) 
                    
s=6.d0*xd1**2*yd1
case(67) 
                    
s=12.d0*xd1*yd1**2
case(68) 
                    
s=20*yd1**3
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=4.d0*xd1**2*zd1
case(72) 
                    
s=12.d0*xd1*yd1*zd1
case(73) 
                    
s=24.d0*yd1**2*zd1
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=6.d0*xd1*zd1**2
case(77) 
                    
s=18.d0*yd1*zd1**2
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=8.d0*zd1**3
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy2z=s
end function dfy2z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=2
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=4.d0*xd1
case(30) 
                    
s=2.d0*yd1
case(31) 
                    
s=0.0d0
case(32) 
                    
s=6.d0*zd1
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=6.d0*xd1**2
case(47) 
                    
s=4.d0*xd1*yd1
case(48) 
                    
s=2.d0*yd1**2
case(49) 
                    
s=0.0d0
case(50) 
                    
s=12.d0*xd1*zd1
case(51) 
                    
s=6.d0*yd1*zd1
case(52) 
                    
s=0.0d0
case(53) 
                    
s=12.d0*zd1**2
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=8.d0*xd1**3
case(70) 
                    
s=6.d0*xd1**2*yd1
case(71) 
                    
s=4.d0*xd1*yd1**2
case(72) 
                    
s=2.d0*yd1**3
case(73) 
                    
s=0.0d0
case(74) 
                    
s=18.d0*xd1**2*zd1
case(75) 
                    
s=12.d0*xd1*yd1*zd1
case(76) 
                    
s=6.d0*yd1**2*zd1
case(77) 
                    
s=0.0d0
case(78) 
                    
s=24.d0*xd1*zd1**2
case(79) 
                    
s=12.d0*yd1*zd1**2
case(80) 
                    
s=0.0d0
case(81) 
                    
s=20*zd1**3
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxz2=s
end function dfxz2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfyz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1) 
                    
s=0.0d0
case (2) 
                    
s=0.0d0
case (3) 
                    
s=0.0d0
case (4) 
                    
s=0.0d0
case (5) 
                    
s=0.0d0
case (6) 
                    
s=0.0d0
case (7) 
                    
s=0.0d0
case (8) 
                    
s=0.0d0
case (9) 
                    
s=0.0d0
case (10) 
                    
s=0.0d0
case (11) 
                    
s=0.0d0
case (12) 
                    
s=0.0d0
case (13) 
                    
s=0.0d0
case (14) 
                    
s=0.0d0
case (15) 
                    
s=0.0d0
case (16) 
                    
s=0.0d0
case (17) 
                    
s=0.0d0
case (18) 
                    
s=2
case (19) 
                    
s=0.0d0
case (20) 
                    
s=0.0d0
case (21) 
                    
s=0.0d0
case (22) 
                    
s=0.0d0
case (23) 
                    
s=0.0d0
case (24) 
                    
s=0.0d0
case (25) 
                    
s=0.0d0
case (26) 
                    
s=0.0d0
case (27) 
                    
s=0.0d0
case (28) 
                    
s=0.0d0
case (29) 
                    
s=0.0d0
case (30) 
                    
s=2.d0*xd1
case (31) 
                    
s=4.d0*yd1
case (32) 
                    
s=0.0d0
case (33) 
                    
s=6.d0*zd1
case (34) 
                    
s=0.0d0
case (35) 
                    
s=0.0d0
case (36) 
                    
s=0.0d0
case (37) 
                    
s=0.0d0
case (38) 
                    
s=0.0d0
case (39) 
                    
s=0.0d0
case (40) 
                    
s=0.0d0
case (41) 
                    
s=0.0d0
case (42) 
                    
s=0.0d0
case (43) 
                    
s=0.0d0
case (44) 
                    
s=0.0d0
case (45) 
                    
s=0.0d0
case (46) 
                    
s=0.0d0
case (47) 
                    
s=2.d0*xd1**2
case (48) 
                    
s=4.d0*xd1*yd1
case (49) 
                    
s=6.d0*yd1**2
case (50) 
                    
s=0.0d0
case (51) 
                    
s=6.d0*xd1*zd1
case (52) 
                    
s=12.d0*yd1*zd1
case (53) 
                    
s=0.0d0
case (54) 
                    
s=12.d0*zd1**2
case (55) 
                    
s=0.0d0
case (56) 
                    
s=0.0d0
case (57) 
                    
s=0.0d0
case (58) 
                    
s=0.0d0
case (59) 
                    
s=0.0d0
case (60) 
                    
s=0.0d0
case (61) 
                    
s=0.0d0
case (62) 
                    
s=0.0d0
case (63) 
                    
s=0.0d0
case (64) 
                    
s=0.0d0
case (65) 
                    
s=0.0d0
case (66) 
                    
s=0.0d0
case (67) 
                    
s=0.0d0
case (68) 
                    
s=0.0d0
case (69) 
                    
s=0.0d0
case (70) 
                    
s=2.d0*xd1**3
case (71) 
                    
s=4.d0*xd1**2*yd1
case (72) 
                    
s=6.d0*xd1*yd1**2
case (73) 
                    
s=8.d0*yd1**3
case (74) 
                    
s=0.0d0
case (75) 
                    
s=6.d0*xd1**2*zd1
case (76) 
                    
s=12.d0*xd1*yd1*zd1
case (77) 
                    
s=18.d0*yd1**2*zd1
case (78) 
                    
s=0.0d0
case (79) 
                    
s=12.d0*xd1*zd1**2
case (80) 
                    
s=24.d0*yd1*zd1**2
case (81) 
                    
s=0.0d0
case (82) 
                    
s=20*zd1**3
case (83) 
                    
s=0.0d0
end select
dfyz2=s
end function dfyz2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1)  
                    
s=0.0d0
case (2)  
                    
s=0.0d0
case (3)  
                    
s=0.0d0
case (4)  
                    
s=0.0d0
case (5)  
                    
s=0.0d0
case (6)  
                    
s=0.0d0
case (7)  
                    
s=0.0d0
case (8)  
                    
s=0.0d0
case (9)  
                    
s=0.0d0
case (10)  
                    
s=0.0d0
case (11)  
                    
s=0.0d0
case (12)  
                    
s=0.0d0
case (13)  
                    
s=0.0d0
case (14)  
                    
s=0.0d0
case (15)  
                    
s=0.0d0
case (16)  
                    
s=0.0d0
case (17)  
                    
s=0.0d0
case (18)  
                    
s=0.0d0
case (19)  
                    
s=6
case (20)  
                    
s=0.0d0
case (21)  
                    
s=0.0d0
case (22)  
                    
s=0.0d0
case (23)  
                    
s=0.0d0
case (24)  
                    
s=0.0d0
case (25)  
                    
s=0.0d0
case (26)  
                    
s=0.0d0
case (27)  
                    
s=0.0d0
case (28)  
                    
s=0.0d0
case (29)  
                    
s=0.0d0
case (30)  
                    
s=0.0d0
case (31)  
                    
s=0.0d0
case (32)  
                    
s=6.d0*xd1
case (33)  
                    
s=6.d0*yd1
case (34)  
                    
s=24.d0*zd1
case (35)  
                    
s=0.0d0
case (36)  
                    
s=0.0d0
case (37)  
                    
s=0.0d0
case (38)  
                    
s=0.0d0
case (39)  
                    
s=0.0d0
case (40)  
                    
s=0.0d0
case (41)  
                    
s=0.0d0
case (42)  
                    
s=0.0d0
case (43)  
                    
s=0.0d0
case (44)  
                    
s=0.0d0
case (45)  
                    
s=0.0d0
case (46)  
                    
s=0.0d0
case (47)  
                    
s=0.0d0
case (48)  
                    
s=0.0d0
case (49)  
                    
s=0.0d0
case (50)  
                    
s=6.d0*xd1**2
case (51)  
                    
s=6.d0*xd1*yd1
case (52)  
                    
s=6.d0*yd1**2
case (53)  
                    
s=24.d0*xd1*zd1
case (54)  
                    
s=24.d0*yd1*zd1
case (55)  
                    
s=60*zd1**2
case (56)  
                    
s=0.0d0
case (57)  
                    
s=0.0d0
case (58)  
                    
s=0.0d0
case (59)  
                    
s=0.0d0
case (60)  
                    
s=0.0d0
case (61)  
                    
s=0.0d0
case (62)  
                    
s=0.0d0
case (63)  
                    
s=0.0d0
case (64)  
                    
s=0.0d0
case (65)  
                    
s=0.0d0
case (66)  
                    
s=0.0d0
case (67)  
                    
s=0.0d0
case (68)  
                    
s=0.0d0
case (69)  
                    
s=0.0d0
case (70)  
                    
s=0.0d0
case (71)  
                    
s=0.0d0
case (72)  
                    
s=0.0d0
case (73)  
                    
s=0.0d0
case (74)  
                    
s=6.d0*xd1**3
case (75)  
                    
s=6.d0*xd1**2*yd1
case (76)  
                    
s=6.d0*xd1*yd1**2
case (77)  
                    
s=6.d0*yd1**3
case (78)  
                    
s=24.d0*xd1**2*zd1
case (79)  
                    
s=24.d0*xd1*yd1*zd1
case (80)  
                    
s=24.d0*yd1**2*zd1
case (81)  
                    
s=60*xd1*zd1**2
case (82)  
                    
s=60*yd1*zd1**2
case (83)  
                    
s=120*zd1**3
end select
dfz3=s
end function dfz3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=24
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=120*xd1
case(36) 
                    
s=24.d0*yd1
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=24.d0*zd1
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=360*xd1**2
case(57) 
                    
s=120*xd1*yd1
case(58) 
                    
s=24.d0*yd1**2
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=120*xd1*zd1
case(64) 
                    
s=24.d0*yd1*zd1
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=24.d0*zd1**2
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx4=s
end function dfx4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=6
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=24.d0*xd1
case(37) 
                    
s=12.d0*yd1
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=6.d0*zd1
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=60*xd1**2
case(58) 
                    
s=48.d0*xd1*yd1
case(59) 
                    
s=18.d0*yd1**2
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=24.d0*xd1*zd1
case(65) 
                    
s=12.d0*yd1*zd1
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=6.d0*zd1**2
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3y=s
end function dfx3y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=4
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=12.d0*xd1
case(38) 
                    
s=12.d0*yd1
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=4.d0*zd1
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=24.d0*xd1**2
case(59) 
                    
s=36.d0*xd1*yd1
case(60) 
                    
s=24.d0*yd1**2
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=12.d0*xd1*zd1
case(66) 
                    
s=12.d0*yd1*zd1
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=4.d0*zd1**2
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2y2=s
end function dfx2y2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=6
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=12.d0*xd1
case(39) 
                    
s=24.d0*yd1
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=6.d0*zd1
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=18.d0*xd1**2
case(60) 
                    
s=48.d0*xd1*yd1
case(61) 
                    
s=60*yd1**2
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=12.d0*xd1*zd1
case(67) 
                    
s=24.d0*yd1*zd1
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=6.d0*zd1**2
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy3=s
end function dfxy3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=24
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=24.d0*xd1
case(40) 
                    
s=120*yd1
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=24.d0*zd1
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=24.d0*xd1**2
case(61) 
                    
s=120*xd1*yd1
case(62) 
                    
s=360*yd1**2
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=24.d0*xd1*zd1
case(68) 
                    
s=120*yd1*zd1
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=24.d0*zd1**2
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy4=s
end function dfy4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=6
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=24.d0*xd1
case(42) 
                    
s=6.d0*yd1
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=12.d0*zd1
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=60*xd1**2
case(64) 
                    
s=24.d0*xd1*yd1
case(65) 
                    
s=6.d0*yd1**2
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=48.d0*xd1*zd1
case(70) 
                    
s=12.d0*yd1*zd1
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=18.d0*zd1**2
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3z=s
end function dfx3z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=2
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=6.d0*xd1
case(43) 
                    
s=4.d0*yd1
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=4.d0*zd1
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=12.d0*xd1**2
case(65) 
                    
s=12.d0*xd1*yd1
case(66) 
                    
s=6.d0*yd1**2
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=12.d0*xd1*zd1
case(71) 
                    
s=8.d0*yd1*zd1
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=6.d0*zd1**2
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2yz=s
end function dfx2yz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=2
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=4.d0*xd1
case(44) 
                    
s=6.d0*yd1
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=4.d0*zd1
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=6.d0*xd1**2
case(66) 
                    
s=12.d0*xd1*yd1
case(67) 
                    
s=12.d0*yd1**2
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=8.d0*xd1*zd1
case(72) 
                    
s=12.d0*yd1*zd1
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=6.d0*zd1**2
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy2z=s
end function dfxy2z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=6
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=6.d0*xd1
case(45) 
                    
s=24.d0*yd1
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=12.d0*zd1
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=6.d0*xd1**2
case(67) 
                    
s=24.d0*xd1*yd1
case(68) 
                    
s=60*yd1**2
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=12.d0*xd1*zd1
case(73) 
                    
s=48.d0*yd1*zd1
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=18.d0*zd1**2
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy3z=s
end function dfy3z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=4
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=12.d0*xd1
case(47) 
                    
s=4.d0*yd1
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=12.d0*zd1
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=24.d0*xd1**2
case(70) 
                    
s=12.d0*xd1*yd1
case(71) 
                    
s=4.d0*yd1**2
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=36.d0*xd1*zd1
case(75) 
                    
s=12.d0*yd1*zd1
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=24.d0*zd1**2
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2z2=s
end function dfx2z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxyz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=2
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=4.d0*xd1
case(48) 
                    
s=4.d0*yd1
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=6.d0*zd1
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=6.d0*xd1**2
case(71) 
                    
s=8.d0*xd1*yd1
case(72) 
                    
s=6.d0*yd1**2
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=12.d0*xd1*zd1
case(76) 
                    
s=12.d0*yd1*zd1
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=12.d0*zd1**2
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxyz2=s
end function dfxyz2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=4
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=4.d0*xd1
case(49) 
                    
s=12.d0*yd1
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=12.d0*zd1
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=4.d0*xd1**2
case(72) 
                    
s=12.d0*xd1*yd1
case(73) 
                    
s=24.d0*yd1**2
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=12.d0*xd1*zd1
case(77) 
                    
s=36.d0*yd1*zd1
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=24.d0*zd1**2
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy2z2=s
end function dfy2z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=6
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=12.d0*xd1
case(51) 
                    
s=6.d0*yd1
case(52) 
                    
s=0.0d0
case(53) 
                    
s=24.d0*zd1
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=18.d0*xd1**2
case(75) 
                    
s=12.d0*xd1*yd1
case(76) 
                    
s=6.d0*yd1**2
case(77) 
                    
s=0.0d0
case(78) 
                    
s=48.d0*xd1*zd1
case(79) 
                    
s=24.d0*yd1*zd1
case(80) 
                    
s=0.0d0
case(81) 
                    
s=60*zd1**2
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxz3=s
end function dfxz3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfyz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1)  
                    
s=0.0d0
case (2)  
                    
s=0.0d0
case (3)  
                    
s=0.0d0
case (4)  
                    
s=0.0d0
case (5)  
                    
s=0.0d0
case (6)  
                    
s=0.0d0
case (7)  
                    
s=0.0d0
case (8)  
                    
s=0.0d0
case (9)  
                    
s=0.0d0
case (10)  
                    
s=0.0d0
case (11)  
                    
s=0.0d0
case (12)  
                    
s=0.0d0
case (13)  
                    
s=0.0d0
case (14)  
                    
s=0.0d0
case (15)  
                    
s=0.0d0
case (16)  
                    
s=0.0d0
case (17)  
                    
s=0.0d0
case (18)  
                    
s=0.0d0
case (19)  
                    
s=0.0d0
case (20)  
                    
s=0.0d0
case (21)  
                    
s=0.0d0
case (22)  
                    
s=0.0d0
case (23)  
                    
s=0.0d0
case (24)  
                    
s=0.0d0
case (25)  
                    
s=0.0d0
case (26)  
                    
s=0.0d0
case (27)  
                    
s=0.0d0
case (28)  
                    
s=0.0d0
case (29)  
                    
s=0.0d0
case (30)  
                    
s=0.0d0
case (31)  
                    
s=0.0d0
case (32)  
                    
s=0.0d0
case (33)  
                    
s=6
case (34)  
                    
s=0.0d0
case (35)  
                    
s=0.0d0
case (36)  
                    
s=0.0d0
case (37)  
                    
s=0.0d0
case (38)  
                    
s=0.0d0
case (39)  
                    
s=0.0d0
case (40)  
                    
s=0.0d0
case (41)  
                    
s=0.0d0
case (42)  
                    
s=0.0d0
case (43)  
                    
s=0.0d0
case (44)  
                    
s=0.0d0
case (45)  
                    
s=0.0d0
case (46)  
                    
s=0.0d0
case (47)  
                    
s=0.0d0
case (48)  
                    
s=0.0d0
case (49)  
                    
s=0.0d0
case (50)  
                    
s=0.0d0
case (51)  
                    
s=6.d0*xd1
case (52)  
                    
s=12.d0*yd1
case (53)  
                    
s=0.0d0
case (54)  
                    
s=24.d0*zd1
case (55)  
                    
s=0.0d0
case (56)  
                    
s=0.0d0
case (57)  
                    
s=0.0d0
case (58)  
                    
s=0.0d0
case (59)  
                    
s=0.0d0
case (60)  
                    
s=0.0d0
case (61)  
                    
s=0.0d0
case (62)  
                    
s=0.0d0
case (63)  
                    
s=0.0d0
case (64)  
                    
s=0.0d0
case (65)  
                    
s=0.0d0
case (66)  
                    
s=0.0d0
case (67)  
                    
s=0.0d0
case (68)  
                    
s=0.0d0
case (69)  
                    
s=0.0d0
case (70)  
                    
s=0.0d0
case (71)  
                    
s=0.0d0
case (72)  
                    
s=0.0d0
case (73)  
                    
s=0.0d0
case (74)  
                    
s=0.0d0
case (75)  
                    
s=6.d0*xd1**2
case (76)  
                    
s=12.d0*xd1*yd1
case (77)  
                    
s=18.d0*yd1**2
case (78)  
                    
s=0.0d0
case (79)  
                    
s=24.d0*xd1*zd1
case (80)  
                    
s=48.d0*yd1*zd1
case (81)  
                    
s=0.0d0
case (82)  
                    
s=60*zd1**2
case (83)  
                    
s=0.0d0
end select
dfyz3=s
end function dfyz3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1)  
                    
s=0.0d0
case (2)  
                    
s=0.0d0
case (3)  
                    
s=0.0d0
case (4)  
                    
s=0.0d0
case (5)  
                    
s=0.0d0
case (6)  
                    
s=0.0d0
case (7)  
                    
s=0.0d0
case (8)  
                    
s=0.0d0
case (9)  
                    
s=0.0d0
case (10)  
                    
s=0.0d0
case (11)  
                    
s=0.0d0
case (12)  
                    
s=0.0d0
case (13)  
                    
s=0.0d0
case (14)  
                    
s=0.0d0
case (15)  
                    
s=0.0d0
case (16)  
                    
s=0.0d0
case (17)  
                    
s=0.0d0
case (18)  
                    
s=0.0d0
case (19)  
                    
s=0.0d0
case (20)  
                    
s=0.0d0
case (21)  
                    
s=0.0d0
case (22)  
                    
s=0.0d0
case (23)  
                    
s=0.0d0
case (24)  
                    
s=0.0d0
case (25)  
                    
s=0.0d0
case (26)  
                    
s=0.0d0
case (27)  
                    
s=0.0d0
case (28)  
                    
s=0.0d0
case (29)  
                    
s=0.0d0
case (30)  
                    
s=0.0d0
case (31)  
                    
s=0.0d0
case (32)  
                    
s=0.0d0
case (33)  
                    
s=0.0d0
case (34)  
                    
s=24
case (35)  
                    
s=0.0d0
case (36)  
                    
s=0.0d0
case (37)  
                    
s=0.0d0
case (38)  
                    
s=0.0d0
case (39)  
                    
s=0.0d0
case (40)  
                    
s=0.0d0
case (41)  
                    
s=0.0d0
case (42)  
                    
s=0.0d0
case (43)  
                    
s=0.0d0
case (44)  
                    
s=0.0d0
case (45)  
                    
s=0.0d0
case (46)  
                    
s=0.0d0
case (47)  
                    
s=0.0d0
case (48)  
                    
s=0.0d0
case (49)  
                    
s=0.0d0
case (50)  
                    
s=0.0d0
case (51)  
                    
s=0.0d0
case (52)  
                    
s=0.0d0
case (53)  
                    
s=24.d0*xd1
case (54)  
                    
s=24.d0*yd1
case (55)  
                    
s=120*zd1
case (56)  
                    
s=0.0d0
case (57)  
                    
s=0.0d0
case (58)  
                    
s=0.0d0
case (59)  
                    
s=0.0d0
case (60)  
                    
s=0.0d0
case (61)  
                    
s=0.0d0
case (62)  
                    
s=0.0d0
case (63)  
                    
s=0.0d0
case (64)  
                    
s=0.0d0
case (65)  
                    
s=0.0d0
case (66)  
                    
s=0.0d0
case (67)  
                    
s=0.0d0
case (68)  
                    
s=0.0d0
case (69)  
                    
s=0.0d0
case (70)  
                    
s=0.0d0
case (71)  
                    
s=0.0d0
case (72)  
                    
s=0.0d0
case (73)  
                    
s=0.0d0
case (74)  
                    
s=0.0d0
case (75)  
                    
s=0.0d0
case (76)  
                    
s=0.0d0
case (77)  
                    
s=0.0d0
case (78)  
                    
s=24.d0*xd1**2
case (79)  
                    
s=24.d0*xd1*yd1
case (80)  
                    
s=24.d0*yd1**2
case (81)  
                    
s=120*xd1*zd1
case (82)  
                    
s=120*yd1*zd1
case (83)  
                    
s=360*zd1**2
end select
dfz4=s
end function dfz4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=120
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=720*xd1
case(57) 
                    
s=120*yd1
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=120*zd1
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx5=s
end function dfx5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx4y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=24
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=120*xd1
case(58) 
                    
s=48.d0*yd1
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=24.d0*zd1
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx4y=s
end function dfx4y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=12
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=48.d0*xd1
case(59) 
                    
s=36.d0*yd1
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=12.d0*zd1
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3y2=s
end function dfx3y2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2y3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=12
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=36.d0*xd1
case(60) 
                    
s=48.d0*yd1
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=12.d0*zd1
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2y3=s
end function dfx2y3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=24
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=48.d0*xd1
case(61) 
                    
s=120*yd1
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=24.d0*zd1
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy4=s
end function dfxy4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=120
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=120*xd1
case(62) 
                    
s=720*yd1
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=120*zd1
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy5=s
end function dfy5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=24
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=120*xd1
case(64) 
                    
s=24.d0*yd1
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=48.d0*zd1
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx4z=s
end function dfx4z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=6
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=24.d0*xd1
case(65) 
                    
s=12.d0*yd1
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=12.d0*zd1
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3yz=s
end function dfx3yz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2y2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=4
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=12.d0*xd1
case(66) 
                    
s=12.d0*yd1
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=8.d0*zd1
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2y2z=s
end function dfx2y2z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=6
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=12.d0*xd1
case(67) 
                    
s=24.d0*yd1
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=12.d0*zd1
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy3z=s
end function dfxy3z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=24
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=24.d0*xd1
case(68) 
                    
s=120*yd1
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=48.d0*zd1
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy4z=s
end function dfy4z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=12
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=48.d0*xd1
case(70) 
                    
s=12.d0*yd1
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=36.d0*zd1
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3z2=s
end function dfx3z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2yz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=4
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=12.d0*xd1
case(71) 
                    
s=8.d0*yd1
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=12.d0*zd1
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2yz2=s
end function dfx2yz2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=4
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=8.d0*xd1
case(72) 
                    
s=12.d0*yd1
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=12.d0*zd1
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy2z2=s
end function dfxy2z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=12
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=12.d0*xd1
case(73) 
                    
s=48.d0*yd1
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=36.d0*zd1
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy3z2=s
end function dfy3z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=12
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=36.d0*xd1
case(75) 
                    
s=12.d0*yd1
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=48.d0*zd1
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2z3=s
end function dfx2z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxyz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=6
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=12.d0*xd1
case(76) 
                    
s=12.d0*yd1
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=24.d0*zd1
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxyz3=s
end function dfxyz3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=12
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=12.d0*xd1
case(77) 
                    
s=36.d0*yd1
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=48.d0*zd1
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy2z3=s
end function dfy2z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=24
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=48.d0*xd1
case(79) 
                    
s=24.d0*yd1
case(80) 
                    
s=0.0d0
case(81) 
                    
s=120*zd1
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxz4=s
end function dfxz4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfyz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1)  
                    
s=0.0d0
case (2)  
                    
s=0.0d0
case (3)  
                    
s=0.0d0
case (4)  
                    
s=0.0d0
case (5)  
                    
s=0.0d0
case (6)  
                    
s=0.0d0
case (7)  
                    
s=0.0d0
case (8)  
                    
s=0.0d0
case (9)  
                    
s=0.0d0
case (10)  
                    
s=0.0d0
case (11)  
                    
s=0.0d0
case (12)  
                    
s=0.0d0
case (13)  
                    
s=0.0d0
case (14)  
                    
s=0.0d0
case (15)  
                    
s=0.0d0
case (16)  
                    
s=0.0d0
case (17)  
                    
s=0.0d0
case (18)  
                    
s=0.0d0
case (19)  
                    
s=0.0d0
case (20)  
                    
s=0.0d0
case (21)  
                    
s=0.0d0
case (22)  
                    
s=0.0d0
case (23)  
                    
s=0.0d0
case (24)  
                    
s=0.0d0
case (25)  
                    
s=0.0d0
case (26)  
                    
s=0.0d0
case (27)  
                    
s=0.0d0
case (28)  
                    
s=0.0d0
case (29)  
                    
s=0.0d0
case (30)  
                    
s=0.0d0
case (31)  
                    
s=0.0d0
case (32)  
                    
s=0.0d0
case (33)  
                    
s=0.0d0
case (34)  
                    
s=0.0d0
case (35)  
                    
s=0.0d0
case (36)  
                    
s=0.0d0
case (37)  
                    
s=0.0d0
case (38)  
                    
s=0.0d0
case (39)  
                    
s=0.0d0
case (40)  
                    
s=0.0d0
case (41)  
                    
s=0.0d0
case (42)  
                    
s=0.0d0
case (43)  
                    
s=0.0d0
case (44)  
                    
s=0.0d0
case (45)  
                    
s=0.0d0
case (46)  
                    
s=0.0d0
case (47)  
                    
s=0.0d0
case (48)  
                    
s=0.0d0
case (49)  
                    
s=0.0d0
case (50)  
                    
s=0.0d0
case (51)  
                    
s=0.0d0
case (52)  
                    
s=0.0d0
case (53)  
                    
s=0.0d0
case (54)  
                    
s=24
case (55)  
                    
s=0.0d0
case (56)  
                    
s=0.0d0
case (57)  
                    
s=0.0d0
case (58)  
                    
s=0.0d0
case (59)  
                    
s=0.0d0
case (60)  
                    
s=0.0d0
case (61)  
                    
s=0.0d0
case (62)  
                    
s=0.0d0
case (63)  
                    
s=0.0d0
case (64)  
                    
s=0.0d0
case (65)  
                    
s=0.0d0
case (66)  
                    
s=0.0d0
case (67)  
                    
s=0.0d0
case (68)  
                    
s=0.0d0
case (69)  
                    
s=0.0d0
case (70)  
                    
s=0.0d0
case (71)  
                    
s=0.0d0
case (72)  
                    
s=0.0d0
case (73)  
                    
s=0.0d0
case (74)  
                    
s=0.0d0
case (75)  
                    
s=0.0d0
case (76)  
                    
s=0.0d0
case (77)  
                    
s=0.0d0
case (78)  
                    
s=0.0d0
case (79)  
                    
s=24.d0*xd1
case (80)  
                    
s=48.d0*yd1
case (81)  
                    
s=0.0d0
case (82)  
                    
s=120*zd1
case (83)  
                    
s=0.0d0
end select
dfyz4=s
end function dfyz4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfz5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1)  
                    
s=0.0d0
case(2)  
                    
s=0.0d0
case(3)  
                    
s=0.0d0
case(4)  
                    
s=0.0d0
case(5)  
                    
s=0.0d0
case(6)  
                    
s=0.0d0
case(7)  
                    
s=0.0d0
case(8)  
                    
s=0.0d0
case(9)  
                    
s=0.0d0
case(10)  
                    
s=0.0d0
case(11)  
                    
s=0.0d0
case(12)  
                    
s=0.0d0
case(13)  
                    
s=0.0d0
case(14)  
                    
s=0.0d0
case(15)  
                    
s=0.0d0
case(16)  
                    
s=0.0d0
case(17)  
                    
s=0.0d0
case(18)  
                    
s=0.0d0
case(19)  
                    
s=0.0d0
case(20)  
                    
s=0.0d0
case(21)  
                    
s=0.0d0
case(22)  
                    
s=0.0d0
case(23)  
                    
s=0.0d0
case(24)  
                    
s=0.0d0
case(25)  
                    
s=0.0d0
case(26)  
                    
s=0.0d0
case(27)  
                    
s=0.0d0
case(28)  
                    
s=0.0d0
case(29)  
                    
s=0.0d0
case(30)  
                    
s=0.0d0
case(31)  
                    
s=0.0d0
case(32)  
                    
s=0.0d0
case(33)  
                    
s=0.0d0
case(34)  
                    
s=0.0d0
case(35)  
                    
s=0.0d0
case(36)  
                    
s=0.0d0
case(37)  
                    
s=0.0d0
case(38)  
                    
s=0.0d0
case(39)  
                    
s=0.0d0
case(40)  
                    
s=0.0d0
case(41)  
                    
s=0.0d0
case(42)  
                    
s=0.0d0
case(43)  
                    
s=0.0d0
case(44)  
                    
s=0.0d0
case(45)  
                    
s=0.0d0
case(46)  
                    
s=0.0d0
case(47)  
                    
s=0.0d0
case(48)  
                    
s=0.0d0
case(49)  
                    
s=0.0d0
case(50)  
                    
s=0.0d0
case(51)  
                    
s=0.0d0
case(52)  
                    
s=0.0d0
case(53)  
                    
s=0.0d0
case(54)  
                    
s=0.0d0
case(55)  
                    
s=120
case(56)  
                    
s=0.0d0
case(57)  
                    
s=0.0d0
case(58)  
                    
s=0.0d0
case(59)  
                    
s=0.0d0
case(60)  
                    
s=0.0d0
case(61)  
                    
s=0.0d0
case(62)  
                    
s=0.0d0
case(63)  
                    
s=0.0d0
case(64)  
                    
s=0.0d0
case(65)  
                    
s=0.0d0
case(66)  
                    
s=0.0d0
case(67)  
                    
s=0.0d0
case(68)  
                    
s=0.0d0
case(69)  
                    
s=0.0d0
case(70)  
                    
s=0.0d0
case(71)  
                    
s=0.0d0
case(72)  
                    
s=0.0d0
case(73)  
                    
s=0.0d0
case(74)  
                    
s=0.0d0
case(75)  
                    
s=0.0d0
case(76)  
                    
s=0.0d0
case(77)  
                    
s=0.0d0
case(78)  
                    
s=0.0d0
case(79)  
                    
s=0.0d0
case(80)  
                    
s=0.0d0
case(81)  
                    
s=120*xd1
case(82)  
                    
s=120*yd1
case(83)  
                    
s=720*zd1
end select
dfz5=s
end function dfz5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=720
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx6=s
end function dfx6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx5y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=120
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx5y=s
end function dfx5y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx4y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=48
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx4y2=s
end function dfx4y2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3y3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=36
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3y3=s
end function dfx3y3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2y4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=48
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2y4=s
end function dfx2y4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1) 
                    
s=0.0d0
case (2) 
                    
s=0.0d0
case (3) 
                    
s=0.0d0
case (4) 
                    
s=0.0d0
case (5) 
                    
s=0.0d0
case (6) 
                    
s=0.0d0
case (7) 
                    
s=0.0d0
case (8) 
                    
s=0.0d0
case (9) 
                    
s=0.0d0
case (10) 
                    
s=0.0d0
case (11) 
                    
s=0.0d0
case (12) 
                    
s=0.0d0
case (13) 
                    
s=0.0d0
case (14) 
                    
s=0.0d0
case (15) 
                    
s=0.0d0
case (16) 
                    
s=0.0d0
case (17) 
                    
s=0.0d0
case (18) 
                    
s=0.0d0
case (19) 
                    
s=0.0d0
case (20) 
                    
s=0.0d0
case (21) 
                    
s=0.0d0
case (22) 
                    
s=0.0d0
case (23) 
                    
s=0.0d0
case (24) 
                    
s=0.0d0
case (25) 
                    
s=0.0d0
case (26) 
                    
s=0.0d0
case (27) 
                    
s=0.0d0
case (28) 
                    
s=0.0d0
case (29) 
                    
s=0.0d0
case (30) 
                    
s=0.0d0
case (31) 
                    
s=0.0d0
case (32) 
                    
s=0.0d0
case (33) 
                    
s=0.0d0
case (34) 
                    
s=0.0d0
case (35) 
                    
s=0.0d0
case (36) 
                    
s=0.0d0
case (37) 
                    
s=0.0d0
case (38) 
                    
s=0.0d0
case (39) 
                    
s=0.0d0
case (40) 
                    
s=0.0d0
case (41) 
                    
s=0.0d0
case (42) 
                    
s=0.0d0
case (43) 
                    
s=0.0d0
case (44) 
                    
s=0.0d0
case (45) 
                    
s=0.0d0
case (46) 
                    
s=0.0d0
case (47) 
                    
s=0.0d0
case (48) 
                    
s=0.0d0
case (49) 
                    
s=0.0d0
case (50) 
                    
s=0.0d0
case (51) 
                    
s=0.0d0
case (52) 
                    
s=0.0d0
case (53) 
                    
s=0.0d0
case (54) 
                    
s=0.0d0
case (55) 
                    
s=0.0d0
case (56) 
                    
s=0.0d0
case (57) 
                    
s=0.0d0
case (58) 
                    
s=0.0d0
case (59) 
                    
s=0.0d0
case (60) 
                    
s=0.0d0
case (61) 
                    
s=120
case (62) 
                    
s=0.0d0
case (63) 
                    
s=0.0d0
case (64) 
                    
s=0.0d0
case (65) 
                    
s=0.0d0
case (66) 
                    
s=0.0d0
case (67) 
                    
s=0.0d0
case (68) 
                    
s=0.0d0
case (69) 
                    
s=0.0d0
case (70) 
                    
s=0.0d0
case (71) 
                    
s=0.0d0
case (72) 
                    
s=0.0d0
case (73) 
                    
s=0.0d0
case (74) 
                    
s=0.0d0
case (75) 
                    
s=0.0d0
case (76) 
                    
s=0.0d0
case (77) 
                    
s=0.0d0
case (78) 
                    
s=0.0d0
case (79) 
                    
s=0.0d0
case (80) 
                    
s=0.0d0
case (81) 
                    
s=0.0d0
case (82) 
                    
s=0.0d0
case (83) 
                    
s=0.0d0
end select
dfxy5=s
end function dfxy5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=720
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy6=s
end function dfy6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx5z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=120
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx5z=s
end function dfx5z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx4yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=24
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx4yz=s
end function dfx4yz
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3y2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=12
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3y2z=s
end function dfx3y2z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2y3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1)
                    
s=0.0d0
case(2)
                    
s=0.0d0
case(3)
                    
s=0.0d0
case(4)
                    
s=0.0d0
case(5)
                    
s=0.0d0
case(6)
                    
s=0.0d0
case(7)
                    
s=0.0d0
case(8)
                    
s=0.0d0
case(9)
                    
s=0.0d0
case(10)
                    
s=0.0d0
case(11)
                    
s=0.0d0
case(12)
                    
s=0.0d0
case(13)
                    
s=0.0d0
case(14)
                    
s=0.0d0
case(15)
                    
s=0.0d0
case(16)
                    
s=0.0d0
case(17)
                    
s=0.0d0
case(18)
                    
s=0.0d0
case(19)
                    
s=0.0d0
case(20)
                    
s=0.0d0
case(21)
                    
s=0.0d0
case(22)
                    
s=0.0d0
case(23)
                    
s=0.0d0
case(24)
                    
s=0.0d0
case(25)
                    
s=0.0d0
case(26)
                    
s=0.0d0
case(27)
                    
s=0.0d0
case(28)
                    
s=0.0d0
case(29)
                    
s=0.0d0
case(30)
                    
s=0.0d0
case(31)
                    
s=0.0d0
case(32)
                    
s=0.0d0
case(33)
                    
s=0.0d0
case(34)
                    
s=0.0d0
case(35)
                    
s=0.0d0
case(36)
                    
s=0.0d0
case(37)
                    
s=0.0d0
case(38)
                    
s=0.0d0
case(39)
                    
s=0.0d0
case(40)
                    
s=0.0d0
case(41)
                    
s=0.0d0
case(42)
                    
s=0.0d0
case(43)
                    
s=0.0d0
case(44)
                    
s=0.0d0
case(45)
                    
s=0.0d0
case(46)
                    
s=0.0d0
case(47)
                    
s=0.0d0
case(48)
                    
s=0.0d0
case(49)
                    
s=0.0d0
case(50)
                    
s=0.0d0
case(51)
                    
s=0.0d0
case(52)
                    
s=0.0d0
case(53)
                    
s=0.0d0
case(54)
                    
s=0.0d0
case(55)
                    
s=0.0d0
case(56)
                    
s=0.0d0
case(57)
                    
s=0.0d0
case(58)
                    
s=0.0d0
case(59)
                    
s=0.0d0
case(60)
                    
s=0.0d0
case(61)
                    
s=0.0d0
case(62)
                    
s=0.0d0
case(63)
                    
s=0.0d0
case(64)
                    
s=0.0d0
case(65)
                    
s=0.0d0
case(66)
                    
s=12
case(67)
                    
s=0.0d0
case(68)
                    
s=0.0d0
case(69)
                    
s=0.0d0
case(70)
                    
s=0.0d0
case(71)
                    
s=0.0d0
case(72)
                    
s=0.0d0
case(73)
                    
s=0.0d0
case(74)
                    
s=0.0d0
case(75)
                    
s=0.0d0
case(76)
                    
s=0.0d0
case(77)
                    
s=0.0d0
case(78)
                    
s=0.0d0
case(79)
                    
s=0.0d0
case(80)
                    
s=0.0d0
case(81)
                    
s=0.0d0
case(82)
                    
s=0.0d0
case(83)
                    
s=0.0d0
end select
dfx2y3z=s
end function dfx2y3z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=24
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy4z=s
end function dfxy4z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy5z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=24
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy5z=s
end function dfy5z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx4z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=48
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx4z2=s
end function dfx4z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3yz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=12
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3yz2=s
end function dfx3yz2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2y2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=8
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2y2z2=s
end function dfx2y2z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1) 
                    
s=0.0d0
case (2) 
                    
s=0.0d0
case (3) 
                    
s=0.0d0
case (4) 
                    
s=0.0d0
case (5) 
                    
s=0.0d0
case (6) 
                    
s=0.0d0
case (7) 
                    
s=0.0d0
case (8) 
                    
s=0.0d0
case (9) 
                    
s=0.0d0
case (10) 
                    
s=0.0d0
case (11) 
                    
s=0.0d0
case (12) 
                    
s=0.0d0
case (13) 
                    
s=0.0d0
case (14) 
                    
s=0.0d0
case (15) 
                    
s=0.0d0
case (16) 
                    
s=0.0d0
case (17) 
                    
s=0.0d0
case (18) 
                    
s=0.0d0
case (19) 
                    
s=0.0d0
case (20) 
                    
s=0.0d0
case (21) 
                    
s=0.0d0
case (22) 
                    
s=0.0d0
case (23) 
                    
s=0.0d0
case (24) 
                    
s=0.0d0
case (25) 
                    
s=0.0d0
case (26) 
                    
s=0.0d0
case (27) 
                    
s=0.0d0
case (28) 
                    
s=0.0d0
case (29) 
                    
s=0.0d0
case (30) 
                    
s=0.0d0
case (31) 
                    
s=0.0d0
case (32) 
                    
s=0.0d0
case (33) 
                    
s=0.0d0
case (34) 
                    
s=0.0d0
case (35) 
                    
s=0.0d0
case (36) 
                    
s=0.0d0
case (37) 
                    
s=0.0d0
case (38) 
                    
s=0.0d0
case (39) 
                    
s=0.0d0
case (40) 
                    
s=0.0d0
case (41) 
                    
s=0.0d0
case (42) 
                    
s=0.0d0
case (43) 
                    
s=0.0d0
case (44) 
                    
s=0.0d0
case (45) 
                    
s=0.0d0
case (46) 
                    
s=0.0d0
case (47) 
                    
s=0.0d0
case (48) 
                    
s=0.0d0
case (49) 
                    
s=0.0d0
case (50) 
                    
s=0.0d0
case (51) 
                    
s=0.0d0
case (52) 
                    
s=0.0d0
case (53) 
                    
s=0.0d0
case (54) 
                    
s=0.0d0
case (55) 
                    
s=0.0d0
case (56) 
                    
s=0.0d0
case (57) 
                    
s=0.0d0
case (58) 
                    
s=0.0d0
case (59) 
                    
s=0.0d0
case (60) 
                    
s=0.0d0
case (61) 
                    
s=0.0d0
case (62) 
                    
s=0.0d0
case (63) 
                    
s=0.0d0
case (64) 
                    
s=0.0d0
case (65) 
                    
s=0.0d0
case (66) 
                    
s=0.0d0
case (67) 
                    
s=0.0d0
case (68) 
                    
s=0.0d0
case (69) 
                    
s=0.0d0
case (70) 
                    
s=0.0d0
case (71) 
                    
s=0.0d0
case (72) 
                    
s=12
case (73) 
                    
s=0.0d0
case (74) 
                    
s=0.0d0
case (75) 
                    
s=0.0d0
case (76) 
                    
s=0.0d0
case (77) 
                    
s=0.0d0
case (78) 
                    
s=0.0d0
case (79) 
                    
s=0.0d0
case (80) 
                    
s=0.0d0
case (81) 
                    
s=0.0d0
case (82) 
                    
s=0.0d0
case (83) 
                    
s=0.0d0
end select
dfxy3z2=s
end function dfxy3z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy4z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=48
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy4z2=s
end function dfy4z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx3z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=36
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx3z3=s
end function dfx3z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2yz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=12
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2yz3=s
end function dfx2yz3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxy2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=12
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxy2z3=s
end function dfxy2z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy3z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=36
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy3z3=s
end function dfy3z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfx2z4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=48
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfx2z4=s
end function dfx2z4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxyz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=24
case(80) 
                    
s=0.0d0
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxyz4=s
end function dfxyz4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfy2z4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=48
case(81) 
                    
s=0.0d0
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfy2z4=s
end function dfy2z4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfxz5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1) 
                    
s=0.0d0
case(2) 
                    
s=0.0d0
case(3) 
                    
s=0.0d0
case(4) 
                    
s=0.0d0
case(5) 
                    
s=0.0d0
case(6) 
                    
s=0.0d0
case(7) 
                    
s=0.0d0
case(8) 
                    
s=0.0d0
case(9) 
                    
s=0.0d0
case(10) 
                    
s=0.0d0
case(11) 
                    
s=0.0d0
case(12) 
                    
s=0.0d0
case(13) 
                    
s=0.0d0
case(14) 
                    
s=0.0d0
case(15) 
                    
s=0.0d0
case(16) 
                    
s=0.0d0
case(17) 
                    
s=0.0d0
case(18) 
                    
s=0.0d0
case(19) 
                    
s=0.0d0
case(20) 
                    
s=0.0d0
case(21) 
                    
s=0.0d0
case(22) 
                    
s=0.0d0
case(23) 
                    
s=0.0d0
case(24) 
                    
s=0.0d0
case(25) 
                    
s=0.0d0
case(26) 
                    
s=0.0d0
case(27) 
                    
s=0.0d0
case(28) 
                    
s=0.0d0
case(29) 
                    
s=0.0d0
case(30) 
                    
s=0.0d0
case(31) 
                    
s=0.0d0
case(32) 
                    
s=0.0d0
case(33) 
                    
s=0.0d0
case(34) 
                    
s=0.0d0
case(35) 
                    
s=0.0d0
case(36) 
                    
s=0.0d0
case(37) 
                    
s=0.0d0
case(38) 
                    
s=0.0d0
case(39) 
                    
s=0.0d0
case(40) 
                    
s=0.0d0
case(41) 
                    
s=0.0d0
case(42) 
                    
s=0.0d0
case(43) 
                    
s=0.0d0
case(44) 
                    
s=0.0d0
case(45) 
                    
s=0.0d0
case(46) 
                    
s=0.0d0
case(47) 
                    
s=0.0d0
case(48) 
                    
s=0.0d0
case(49) 
                    
s=0.0d0
case(50) 
                    
s=0.0d0
case(51) 
                    
s=0.0d0
case(52) 
                    
s=0.0d0
case(53) 
                    
s=0.0d0
case(54) 
                    
s=0.0d0
case(55) 
                    
s=0.0d0
case(56) 
                    
s=0.0d0
case(57) 
                    
s=0.0d0
case(58) 
                    
s=0.0d0
case(59) 
                    
s=0.0d0
case(60) 
                    
s=0.0d0
case(61) 
                    
s=0.0d0
case(62) 
                    
s=0.0d0
case(63) 
                    
s=0.0d0
case(64) 
                    
s=0.0d0
case(65) 
                    
s=0.0d0
case(66) 
                    
s=0.0d0
case(67) 
                    
s=0.0d0
case(68) 
                    
s=0.0d0
case(69) 
                    
s=0.0d0
case(70) 
                    
s=0.0d0
case(71) 
                    
s=0.0d0
case(72) 
                    
s=0.0d0
case(73) 
                    
s=0.0d0
case(74) 
                    
s=0.0d0
case(75) 
                    
s=0.0d0
case(76) 
                    
s=0.0d0
case(77) 
                    
s=0.0d0
case(78) 
                    
s=0.0d0
case(79) 
                    
s=0.0d0
case(80) 
                    
s=0.0d0
case(81) 
                    
s=120
case(82) 
                    
s=0.0d0
case(83) 
                    
s=0.0d0
end select
dfxz5=s
end function dfxz5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfyz5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case (1) 
                    
s=0.0d0
case (2) 
                    
s=0.0d0
case (3) 
                    
s=0.0d0
case (4) 
                    
s=0.0d0
case (5) 
                    
s=0.0d0
case (6) 
                    
s=0.0d0
case (7) 
                    
s=0.0d0
case (8) 
                    
s=0.0d0
case (9) 
                    
s=0.0d0
case (10) 
                    
s=0.0d0
case (11) 
                    
s=0.0d0
case (12) 
                    
s=0.0d0
case (13) 
                    
s=0.0d0
case (14) 
                    
s=0.0d0
case (15) 
                    
s=0.0d0
case (16) 
                    
s=0.0d0
case (17) 
                    
s=0.0d0
case (18) 
                    
s=0.0d0
case (19) 
                    
s=0.0d0
case (20) 
                    
s=0.0d0
case (21) 
                    
s=0.0d0
case (22) 
                    
s=0.0d0
case (23) 
                    
s=0.0d0
case (24) 
                    
s=0.0d0
case (25) 
                    
s=0.0d0
case (26) 
                    
s=0.0d0
case (27) 
                    
s=0.0d0
case (28) 
                    
s=0.0d0
case (29) 
                    
s=0.0d0
case (30) 
                    
s=0.0d0
case (31) 
                    
s=0.0d0
case (32) 
                    
s=0.0d0
case (33) 
                    
s=0.0d0
case (34) 
                    
s=0.0d0
case (35) 
                    
s=0.0d0
case (36) 
                    
s=0.0d0
case (37) 
                    
s=0.0d0
case (38) 
                    
s=0.0d0
case (39) 
                    
s=0.0d0
case (40) 
                    
s=0.0d0
case (41) 
                    
s=0.0d0
case (42) 
                    
s=0.0d0
case (43) 
                    
s=0.0d0
case (44) 
                    
s=0.0d0
case (45) 
                    
s=0.0d0
case (46) 
                    
s=0.0d0
case (47) 
                    
s=0.0d0
case (48) 
                    
s=0.0d0
case (49) 
                    
s=0.0d0
case (50) 
                    
s=0.0d0
case (51) 
                    
s=0.0d0
case (52) 
                    
s=0.0d0
case (53) 
                    
s=0.0d0
case (54) 
                    
s=0.0d0
case (55) 
                    
s=0.0d0
case (56) 
                    
s=0.0d0
case (57) 
                    
s=0.0d0
case (58) 
                    
s=0.0d0
case (59) 
                    
s=0.0d0
case (60) 
                    
s=0.0d0
case (61) 
                    
s=0.0d0
case (62) 
                    
s=0.0d0
case (63) 
                    
s=0.0d0
case (64) 
                    
s=0.0d0
case (65) 
                    
s=0.0d0
case (66) 
                    
s=0.0d0
case (67) 
                    
s=0.0d0
case (68) 
                    
s=0.0d0
case (69) 
                    
s=0.0d0
case (70) 
                    
s=0.0d0
case (71) 
                    
s=0.0d0
case (72) 
                    
s=0.0d0
case (73) 
                    
s=0.0d0
case (74) 
                    
s=0.0d0
case (75) 
                    
s=0.0d0
case (76) 
                    
s=0.0d0
case (77) 
                    
s=0.0d0
case (78) 
                    
s=0.0d0
case (79) 
                    
s=0.0d0
case (80) 
                    
s=0.0d0
case (81) 
                    
s=0.0d0
case (82) 
                    
s=120
case (83) 
                    
s=0.0d0
end select
dfyz5=s
end function dfyz5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dfz6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(1)  
                    
s=0.0d0
case(2)  
                    
s=0.0d0
case(3)  
                    
s=0.0d0
case(4)  
                    
s=0.0d0
case(5)  
                    
s=0.0d0
case(6)  
                    
s=0.0d0
case(7)  
                    
s=0.0d0
case(8)  
                    
s=0.0d0
case(9)  
                    
s=0.0d0
case(10)  
                    
s=0.0d0
case(11)  
                    
s=0.0d0
case(12)  
                    
s=0.0d0
case(13)  
                    
s=0.0d0
case(14)  
                    
s=0.0d0
case(15)  
                    
s=0.0d0
case(16)  
                    
s=0.0d0
case(17)  
                    
s=0.0d0
case(18)  
                    
s=0.0d0
case(19)  
                    
s=0.0d0
case(20)  
                    
s=0.0d0
case(21)  
                    
s=0.0d0
case(22)  
                    
s=0.0d0
case(23)  
                    
s=0.0d0
case(24)  
                    
s=0.0d0
case(25)  
                    
s=0.0d0
case(26)  
                    
s=0.0d0
case(27)  
                    
s=0.0d0
case(28)  
                    
s=0.0d0
case(29)  
                    
s=0.0d0
case(30)  
                    
s=0.0d0
case(31)  
                    
s=0.0d0
case(32)  
                    
s=0.0d0
case(33)  
                    
s=0.0d0
case(34)  
                    
s=0.0d0
case(35)  
                    
s=0.0d0
case(36)  
                    
s=0.0d0
case(37)  
                    
s=0.0d0
case(38)  
                    
s=0.0d0
case(39)  
                    
s=0.0d0
case(40)  
                    
s=0.0d0
case(41)  
                    
s=0.0d0
case(42)  
                    
s=0.0d0
case(43)  
                    
s=0.0d0
case(44)  
                    
s=0.0d0
case(45)  
                    
s=0.0d0
case(46)  
                    
s=0.0d0
case(47)  
                    
s=0.0d0
case(48)  
                    
s=0.0d0
case(49)  
                    
s=0.0d0
case(50)  
                    
s=0.0d0
case(51)  
                    
s=0.0d0
case(52)  
                    
s=0.0d0
case(53)  
                    
s=0.0d0
case(54)  
                    
s=0.0d0
case(55)  
                    
s=0.0d0
case(56)  
                    
s=0.0d0
case(57)  
                    
s=0.0d0
case(58)  
                    
s=0.0d0
case(59)  
                    
s=0.0d0
case(60)  
                    
s=0.0d0
case(61)  
                    
s=0.0d0
case(62)  
                    
s=0.0d0
case(63)  
                    
s=0.0d0
case(64)  
                    
s=0.0d0
case(65)  
                    
s=0.0d0
case(66)  
                    
s=0.0d0
case(67)  
                    
s=0.0d0
case(68)  
                    
s=0.0d0
case(69)  
                    
s=0.0d0
case(70)  
                    
s=0.0d0
case(71)  
                    
s=0.0d0
case(72)  
                    
s=0.0d0
case(73)  
                    
s=0.0d0
case(74)  
                    
s=0.0d0
case(75)  
                    
s=0.0d0
case(76)  
                    
s=0.0d0
case(77)  
                    
s=0.0d0
case(78)  
                    
s=0.0d0
case(79)  
                    
s=0.0d0
case(80)  
                    
s=0.0d0
case(81)  
                    
s=0.0d0
case(82)  
                    
s=0.0d0
case(83)  
                    
s=720
end select
dfz6=s
end function dfz6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real function dlx(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)

case(1) 
s= 2
case(6) 
s= -2 + 4.d0*zd1
case(19) 
s= 2 - 12.d0*zd1 + 12.d0*zd1**2
case(29) 
s= -2 + 24.d0*zd1 - 60*zd1**2 + 40*zd1**3
case(55) 
s= 2 - 40*zd1 + 180*zd1**2 - 280*zd1**3 + 140*zd1**4
case(76) 
s= -2 + 60*zd1 - 420*zd1**2 + 1120*zd1**3 - 1260*zd1**4 + 504.d0*zd1**5
case(5) 
s= -2 + 4.d0*yd1
case(17) 
s= 2 - 4.d0*yd1 - 4.d0*zd1 + 8.d0*yd1*zd1
case(28) 
s= -2 + 4.d0*yd1 + 12.d0*zd1 - 24.d0*yd1*zd1 - 12.d0*zd1**2 + 24.d0*yd1*zd1**2
case(47) 
s= 2 - 4.d0*yd1 - 24.d0*zd1 + 48.d0*yd1*zd1 + 60*zd1**2 - 120*yd1*zd1**2 - 40*zd1**3 + 80*yd1*zd1**3
case(75) 
s= -2 + 4.d0*yd1 + 40*zd1 - 80*yd1*zd1 - 180*zd1**2 + 360*yd1*zd1**2 + 280*zd1**3 - 560*yd1*zd1**3 - 140*zd1**4 + 280*yd1*zd1**4
case(12) 
s= 2 - 12.d0*yd1 + 12.d0*yd1**2
case(27) 
s= -2 + 12.d0*yd1 - 12.d0*yd1**2 + 4.d0*zd1 - 24.d0*yd1*zd1 + 24.d0*yd1**2*zd1
case(48) 
s= 2 - 12.d0*yd1 + 12.d0*yd1**2 - 12.d0*zd1 + 72.d0*yd1*zd1 - 72.d0*yd1**2*zd1 + 12.d0*zd1**2 - 72.d0*yd1*zd1**2 + 72.d0*yd1**2*zd1**2
case(74) 
s= -2 + 12.d0*yd1 - 12.d0*yd1**2 + 24.d0*zd1 - 144.d0*yd1*zd1 + 144.d0*yd1**2*zd1 - 60*zd1**2 + 360*yd1*zd1**2 - 360*yd1**2*zd1**2 &
+ 40*zd1**3 - 240*yd1*zd1**3 + 240*yd1**2*zd1**3
case(26) 
s= -2 + 24.d0*yd1 - 60*yd1**2 + 40*yd1**3
case(46) 
s= 2 - 24.d0*yd1 + 60*yd1**2 - 40*yd1**3 - 4.d0*zd1 + 48.d0*yd1*zd1 - 120*yd1**2*zd1 + 80*yd1**3*zd1
case(73) 
s= -2 + 24.d0*yd1 - 60*yd1**2 + 40*yd1**3 + 12.d0*zd1 - 144.d0*yd1*zd1 + 360*yd1**2*zd1 - 240*yd1**3*zd1 - 12.d0*zd1**2 &
+ 144.d0*yd1*zd1**2 - 360*yd1**2*zd1**2 + 240*yd1**3*zd1**2
case(45) 
s= 2 - 40*yd1 + 180*yd1**2 - 280*yd1**3 + 140*yd1**4
case(72) 
s= -2 + 40*yd1 - 180*yd1**2 + 280*yd1**3 - 140*yd1**4 + 4.d0*zd1 - 80*yd1*zd1 + 360*yd1**2*zd1 - 560*yd1**3*zd1 + 280*yd1**4*zd1
case(71) 
s= -2 + 60*yd1 - 420*yd1**2 + 1120*yd1**3 - 1260*yd1**4 + 504.d0*yd1**5
case(4) 
s= -6 + 12.d0*xd1
case(18) 
s= 6 - 12.d0*xd1 - 12.d0*zd1 + 24.d0*xd1*zd1
case(25) 
s= -6 + 12.d0*xd1 + 36.d0*zd1 - 72.d0*xd1*zd1 - 36.d0*zd1**2 + 72.d0*xd1*zd1**2
case(44) 
s= 6 - 12.d0*xd1 - 72.d0*zd1 + 144.d0*xd1*zd1 + 180*zd1**2 - 360*xd1*zd1**2 - 120*zd1**3 + 240*xd1*zd1**3
case(70) 
s= -6 + 12.d0*xd1 + 120*zd1 - 240*xd1*zd1 - 540*zd1**2 + 1080*xd1*zd1**2 + 840*zd1**3 - 1680*xd1*zd1**3 - 420*zd1**4 &
+ 840*xd1*zd1**4
case(11) 
s= 6 - 12.d0*xd1 - 12.d0*yd1 + 24.d0*xd1*yd1
case(24) 
s= -6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 12.d0*zd1 - 24.d0*xd1*zd1 - 24.d0*yd1*zd1 + 48.d0*xd1*yd1*zd1
case(43) 
s= 6 - 12.d0*xd1 - 12.d0*yd1 + 24.d0*xd1*yd1 - 36.d0*zd1 + 72.d0*xd1*zd1 + 72.d0*yd1*zd1 - 144.d0*xd1*yd1*zd1 + 36.d0*zd1**2 - 72.d0*xd1*zd1**2 &
- 72.d0*yd1*zd1**2 + 144.d0*xd1*yd1*zd1**2
case(69) 
s= -6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 - 180*zd1**2 + 360*xd1*zd1**2 &
+ 360*yd1*zd1**2 - 720*xd1*yd1*zd1**2 + 120*zd1**3 - 240*xd1*zd1**3 - 240*yd1*zd1**3 + 480*xd1*yd1*zd1**3
case(23) 
s= -6 + 12.d0*xd1 + 36.d0*yd1 - 72.d0*xd1*yd1 - 36.d0*yd1**2 + 72.d0*xd1*yd1**2
case(42) 
s= 6 - 12.d0*xd1 - 36.d0*yd1 + 72.d0*xd1*yd1 + 36.d0*yd1**2 - 72.d0*xd1*yd1**2 - 12.d0*zd1 + 24.d0*xd1*zd1 + 72.d0*yd1*zd1 - 144.d0*xd1*yd1*zd1 &
- 72.d0*yd1**2*zd1 + 144.d0*xd1*yd1**2*zd1
case(68) 
s= -6 + 12.d0*xd1 + 36.d0*yd1 - 72.d0*xd1*yd1 - 36.d0*yd1**2 + 72.d0*xd1*yd1**2 + 36.d0*zd1 - 72.d0*xd1*zd1 - 216.d0*yd1*zd1 + 432.d0*xd1*yd1*zd1 &
+ 216.d0*yd1**2*zd1 - 432.d0*xd1*yd1**2*zd1 - 36.d0*zd1**2 + 72.d0*xd1*zd1**2 + 216.d0*yd1*zd1**2 - 432.d0*xd1*yd1*zd1**2 - 216.d0*yd1**2*zd1**2&
 + 432.d0*xd1*yd1**2*zd1**2
case(41) 
s= 6 - 12.d0*xd1 - 72.d0*yd1 + 144.d0*xd1*yd1 + 180*yd1**2 - 360*xd1*yd1**2 - 120*yd1**3 + 240*xd1*yd1**3
case(67) 
s= -6 + 12.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 180*yd1**2 + 360*xd1*yd1**2 + 120*yd1**3 - 240*xd1*yd1**3 + 12.d0*zd1 - 24.d0*xd1*zd1 &
- 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 + 360*yd1**2*zd1 - 720*xd1*yd1**2*zd1 - 240*yd1**3*zd1 + 480*xd1*yd1**3*zd1
case(66) 
s= -6 + 12.d0*xd1 + 120*yd1 - 240*xd1*yd1 - 540*yd1**2 + 1080*xd1*yd1**2 + 840*yd1**3 - 1680*xd1*yd1**3 - 420*yd1**4&
 + 840*xd1*yd1**4
case(10) 
s= 12 - 60*xd1 + 60*xd1**2
case(22) 
s= -12 + 60*xd1 - 60*xd1**2 + 24.d0*zd1 - 120*xd1*zd1 + 120*xd1**2*zd1
case(39) 
s= 12 - 60*xd1 + 60*xd1**2 - 72.d0*zd1 + 360*xd1*zd1 - 360*xd1**2*zd1 + 72.d0*zd1**2 - 360*xd1*zd1**2 + 360*xd1**2*zd1**2
case(65) 
s= -12 + 60*xd1 - 60*xd1**2 + 144.d0*zd1 - 720*xd1*zd1 + 720*xd1**2*zd1 - 360*zd1**2 + 1800*xd1*zd1**2 - 1800*xd1**2*zd1**2&
 + 240*zd1**3 - 1200*xd1*zd1**3 + 1200*xd1**2*zd1**3
case(21) 
s= -12 + 60*xd1 - 60*xd1**2 + 24.d0*yd1 - 120*xd1*yd1 + 120*xd1**2*yd1
case(40) 
s= 12 - 60*xd1 + 60*xd1**2 - 24.d0*yd1 + 120*xd1*yd1 - 120*xd1**2*yd1 - 24.d0*zd1 + 120*xd1*zd1 - 120*xd1**2*zd1 + 48.d0*yd1*zd1 &
- 240*xd1*yd1*zd1 + 240*xd1**2*yd1*zd1
case(64) 
s= -12 + 60*xd1 - 60*xd1**2 + 24.d0*yd1 - 120*xd1*yd1 + 120*xd1**2*yd1 + 72.d0*zd1 - 360*xd1*zd1 + 360*xd1**2*zd1 - 144.d0*yd1*zd1 &
+ 720*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 - 72.d0*zd1**2 + 360*xd1*zd1**2 - 360*xd1**2*zd1**2 + 144.d0*yd1*zd1**2 - 720*xd1*yd1*zd1**2&
 + 720*xd1**2*yd1*zd1**2
case(38) 
s= 12 - 60*xd1 + 60*xd1**2 - 72.d0*yd1 + 360*xd1*yd1 - 360*xd1**2*yd1 + 72.d0*yd1**2 - 360*xd1*yd1**2 + 360*xd1**2*yd1**2
case(63) 
s= -12 + 60*xd1 - 60*xd1**2 + 72.d0*yd1 - 360*xd1*yd1 + 360*xd1**2*yd1 - 72.d0*yd1**2 + 360*xd1*yd1**2 - 360*xd1**2*yd1**2 + 24.d0*zd1&
 - 120*xd1*zd1 + 120*xd1**2*zd1 - 144.d0*yd1*zd1 + 720*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 + 144.d0*yd1**2*zd1 - 720*xd1*yd1**2*zd1 &
+ 720*xd1**2*yd1**2*zd1
case(62) 
s= -12 + 60*xd1 - 60*xd1**2 + 144.d0*yd1 - 720*xd1*yd1 + 720*xd1**2*yd1 - 360*yd1**2 + 1800*xd1*yd1**2 - 1800*xd1**2*yd1**2 &
+ 240*yd1**3 - 1200*xd1*yd1**3 + 1200*xd1**2*yd1**3
case(20) 
s= -20.0d0+ 180*xd1 - 420*xd1**2 + 280*xd1**3
case(37) 
s= 20.0d0- 180*xd1 + 420*xd1**2 - 280*xd1**3 - 40*zd1 + 360*xd1*zd1 - 840*xd1**2*zd1 + 560*xd1**3*zd1
case(61) 
s= -20.0d0+ 180*xd1 - 420*xd1**2 + 280*xd1**3 + 120*zd1 - 1080*xd1*zd1 + 2520*xd1**2*zd1 - 1680*xd1**3*zd1 - 120*zd1**2 &
+ 1080*xd1*zd1**2 - 2520*xd1**2*zd1**2 + 1680*xd1**3*zd1**2
case(36) 
s= 20.0d0- 180*xd1 + 420*xd1**2 - 280*xd1**3 - 40*yd1 + 360*xd1*yd1 - 840*xd1**2*yd1 + 560*xd1**3*yd1
case(60) 
s= -20.0d0+ 180*xd1 - 420*xd1**2 + 280*xd1**3 + 40*yd1 - 360*xd1*yd1 + 840*xd1**2*yd1 - 560*xd1**3*yd1 + 40*zd1 - 360*xd1*zd1&
 + 840*xd1**2*zd1 - 560*xd1**3*zd1 - 80*yd1*zd1 + 720*xd1*yd1*zd1 - 1680*xd1**2*yd1*zd1 + 1120*xd1**3*yd1*zd1
case(59) 
s= -20.0d0+ 180*xd1 - 420*xd1**2 + 280*xd1**3 + 120*yd1 - 1080*xd1*yd1 + 2520*xd1**2*yd1 - 1680*xd1**3*yd1 - 120*yd1**2 &
+ 1080*xd1*yd1**2 - 2520*xd1**2*yd1**2 + 1680*xd1**3*yd1**2
case(35) 
s= 30.0d0- 420*xd1 + 1680*xd1**2 - 2520*xd1**3 + 1260*xd1**4
case(58) 
s= -30.0d0+ 420*xd1 - 1680*xd1**2 + 2520*xd1**3 - 1260*xd1**4 + 60*zd1 - 840*xd1*zd1 + 3360*xd1**2*zd1 - 5040*xd1**3*zd1 &
+ 2520*xd1**4*zd1
case(57)  
s=-30.0d0+ 420*xd1 - 1680*xd1**2 + 2520*xd1**3 - 1260*xd1**4 + 60*yd1 - 840*xd1*yd1 + 3360*xd1**2*yd1 - 5040*xd1**3*yd1 &
+ 2520*xd1**4*yd1
case(56) 
s= -42 + 840*xd1 - 5040*xd1**2 + 12600*xd1**3 - 13860*xd1**4 + 5544.d0*xd1**5

case default
s=0.0d0


end select  
	dlx=s    
	end function dlx 


real function dly(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)

case(2)
 s=2
case(8)
 s=-2 + 4.d0*zd1
case(15)
 s=2 - 12.d0*zd1 + 12.d0*zd1**2
case(32)
 s=-2 + 24.d0*zd1 - 60*zd1**2 + 40*zd1**3
case(53)
 s=2 - 40*zd1 + 180*zd1**2 - 280*zd1**3 + 140*zd1**4
case(82)
 s=-2 + 60*zd1 - 420*zd1**2 + 1120*zd1**3 - 1260*zd1**4 + 504.d0*zd1**5
case(7)
 s=-6 + 12.d0*yd1
case(14)
 s=6 - 12.d0*yd1 - 12.d0*zd1 + 24.d0*yd1*zd1
case(31)
 s=-6 + 12.d0*yd1 + 36.d0*zd1 - 72.d0*yd1*zd1 - 36.d0*zd1**2 + 72.d0*yd1*zd1**2
case(52)
 s=6 - 12.d0*yd1 - 72.d0*zd1 + 144.d0*yd1*zd1 + 180*zd1**2 - 360*yd1*zd1**2 - 120*zd1**3 + 240*yd1*zd1**3
case(81)
 s=-6 + 12.d0*yd1 + 120*zd1 - 240*yd1*zd1 - 540*zd1**2 + 1080*yd1*zd1**2 + 840*zd1**3 - 1680*yd1*zd1**3 - 420*zd1**4 &
+ 840*yd1*zd1**4
case(13)
 s=12 - 60*yd1 + 60*yd1**2
case(30)
 s=-12 + 60*yd1 - 60*yd1**2 + 24.d0*zd1 - 120*yd1*zd1 + 120*yd1**2*zd1
case(51)
 s=12 - 60*yd1 + 60*yd1**2 - 72.d0*zd1 + 360*yd1*zd1 - 360*yd1**2*zd1 + 72.d0*zd1**2 - 360*yd1*zd1**2 + 360*yd1**2*zd1**2
case(80)
 s=-12 + 60*yd1 - 60*yd1**2 + 144.d0*zd1 - 720*yd1*zd1 + 720*yd1**2*zd1 - 360*zd1**2 + 1800*yd1*zd1**2 - 1800*yd1**2*zd1**2 &
+ 240*zd1**3 - 1200*yd1*zd1**3 + 1200*yd1**2*zd1**3
case(34)
 s=-20.0d0+ 180*yd1 - 420*yd1**2 + 280*yd1**3
case(50)
 s=20.0d0- 180*yd1 + 420*yd1**2 - 280*yd1**3 - 40*zd1 + 360*yd1*zd1 - 840*yd1**2*zd1 + 560*yd1**3*zd1
case(79)
 s=-20.0d0+ 180*yd1 - 420*yd1**2 + 280*yd1**3 + 120*zd1 - 1080*yd1*zd1 + 2520*yd1**2*zd1 - 1680*yd1**3*zd1 - 120*zd1**2&
 + 1080*yd1*zd1**2 - 2520*yd1**2*zd1**2 + 1680*yd1**3*zd1**2
case(49)
 s=30.0d0- 420*yd1 + 1680*yd1**2 - 2520*yd1**3 + 1260*yd1**4
case(78)
 s=-30.0d0+ 420*yd1 - 1680*yd1**2 + 2520*yd1**3 - 1260*yd1**4 + 60*zd1 - 840*yd1*zd1 + 3360*yd1**2*zd1 - 5040*yd1**3*zd1&
 + 2520*yd1**4*zd1
case(77)
 s=-42 + 840*yd1 - 5040*yd1**2 + 12600*yd1**3 - 13860*yd1**4 + 5544.d0*yd1**5

case(5)
 s=-2 + 4.d0*xd1
case(17)
 s=2 - 4.d0*xd1 - 4.d0*zd1 + 8.d0*xd1*zd1
case(28)
 s=-2 + 4.d0*xd1 + 12.d0*zd1 - 24.d0*xd1*zd1 - 12.d0*zd1**2 + 24.d0*xd1*zd1**2
case(47)
 s=2 - 4.d0*xd1 - 24.d0*zd1 + 48.d0*xd1*zd1 + 60*zd1**2 - 120*xd1*zd1**2 - 40*zd1**3 + 80*xd1*zd1**3
case(75)
 s=-2 + 4.d0*xd1 + 40*zd1 - 80*xd1*zd1 - 180*zd1**2 + 360*xd1*zd1**2 + 280*zd1**3 - 560*xd1*zd1**3 - 140*zd1**4 + 280*xd1*zd1**4
case(12)
 s=6 - 12.d0*xd1 - 12.d0*yd1 + 24.d0*xd1*yd1
case(27)
 s=-6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 12.d0*zd1 - 24.d0*xd1*zd1 - 24.d0*yd1*zd1 + 48.d0*xd1*yd1*zd1
case(48)
 s=6 - 12.d0*xd1 - 12.d0*yd1 + 24.d0*xd1*yd1 - 36.d0*zd1 + 72.d0*xd1*zd1 + 72.d0*yd1*zd1 - 144.d0*xd1*yd1*zd1 + 36.d0*zd1**2 - 72.d0*xd1*zd1**2 &
- 72.d0*yd1*zd1**2 + 144.d0*xd1*yd1*zd1**2
case(74)
 s=-6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 - 180*zd1**2 + 360*xd1*zd1**2 &
+ 360*yd1*zd1**2 - 720*xd1*yd1*zd1**2 + 120*zd1**3 - 240*xd1*zd1**3 - 240*yd1*zd1**3 + 480*xd1*yd1*zd1**3
case(26)
 s=-12 + 24.d0*xd1 + 60*yd1 - 120*xd1*yd1 - 60*yd1**2 + 120*xd1*yd1**2
case(46)
 s=12 - 24.d0*xd1 - 60*yd1 + 120*xd1*yd1 + 60*yd1**2 - 120*xd1*yd1**2 - 24.d0*zd1 + 48.d0*xd1*zd1 + 120*yd1*zd1 - 240*xd1*yd1*zd1 &
- 120*yd1**2*zd1 + 240*xd1*yd1**2*zd1
case(73)
 s=-12 + 24.d0*xd1 + 60*yd1 - 120*xd1*yd1 - 60*yd1**2 + 120*xd1*yd1**2 + 72.d0*zd1 - 144.d0*xd1*zd1 - 360*yd1*zd1 + 720*xd1*yd1*zd1 &
+ 360*yd1**2*zd1 - 720*xd1*yd1**2*zd1 - 72.d0*zd1**2 + 144.d0*xd1*zd1**2 + 360*yd1*zd1**2 - 720*xd1*yd1*zd1**2 - 360*yd1**2*zd1**2 &
+ 720*xd1*yd1**2*zd1**2
case(45)
 s=20.0d0- 40*xd1 - 180*yd1 + 360*xd1*yd1 + 420*yd1**2 - 840*xd1*yd1**2 - 280*yd1**3 + 560*xd1*yd1**3
case(72)
 s=-20.0d0+ 40*xd1 + 180*yd1 - 360*xd1*yd1 - 420*yd1**2 + 840*xd1*yd1**2 + 280*yd1**3 - 560*xd1*yd1**3 + 40*zd1 - 80*xd1*zd1 &
- 360*yd1*zd1 + 720*xd1*yd1*zd1 + 840*yd1**2*zd1 - 1680*xd1*yd1**2*zd1 - 560*yd1**3*zd1 + 1120*xd1*yd1**3*zd1
case(71)
 s=-30.0d0+ 60*xd1 + 420*yd1 - 840*xd1*yd1 - 1680*yd1**2 + 3360*xd1*yd1**2 + 2520*yd1**3 - 5040*xd1*yd1**3 - 1260*yd1**4 &
+ 2520*xd1*yd1**4
case(11)
 s=2 - 12.d0*xd1 + 12.d0*xd1**2
case(24)
 s=-2 + 12.d0*xd1 - 12.d0*xd1**2 + 4.d0*zd1 - 24.d0*xd1*zd1 + 24.d0*xd1**2*zd1
case(43)
 s=2 - 12.d0*xd1 + 12.d0*xd1**2 - 12.d0*zd1 + 72.d0*xd1*zd1 - 72.d0*xd1**2*zd1 + 12.d0*zd1**2 - 72.d0*xd1*zd1**2 + 72.d0*xd1**2*zd1**2
case(69)
 s=-2 + 12.d0*xd1 - 12.d0*xd1**2 + 24.d0*zd1 - 144.d0*xd1*zd1 + 144.d0*xd1**2*zd1 - 60*zd1**2 + 360*xd1*zd1**2 - 360*xd1**2*zd1**2 &
+ 40*zd1**3 - 240*xd1*zd1**3 + 240*xd1**2*zd1**3
case(23)
 s=-6 + 36.d0*xd1 - 36.d0*xd1**2 + 12.d0*yd1 - 72.d0*xd1*yd1 + 72.d0*xd1**2*yd1
case(42)
 s=6 - 36.d0*xd1 + 36.d0*xd1**2 - 12.d0*yd1 + 72.d0*xd1*yd1 - 72.d0*xd1**2*yd1 - 12.d0*zd1 + 72.d0*xd1*zd1 - 72.d0*xd1**2*zd1 + 24.d0*yd1*zd1 &
- 144.d0*xd1*yd1*zd1 + 144.d0*xd1**2*yd1*zd1
case(68)
 s=-6 + 36.d0*xd1 - 36.d0*xd1**2 + 12.d0*yd1 - 72.d0*xd1*yd1 + 72.d0*xd1**2*yd1 + 36.d0*zd1 - 216.d0*xd1*zd1 + 216.d0*xd1**2*zd1 - 72.d0*yd1*zd1 &
+ 432.d0*xd1*yd1*zd1 - 432.d0*xd1**2*yd1*zd1 - 36.d0*zd1**2 + 216.d0*xd1*zd1**2 - 216.d0*xd1**2*zd1**2 + 72.d0*yd1*zd1**2 - 432.d0*xd1*yd1*zd1**2&
 + 432.d0*xd1**2*yd1*zd1**2
case(41)
 s=12 - 72.d0*xd1 + 72.d0*xd1**2 - 60*yd1 + 360*xd1*yd1 - 360*xd1**2*yd1 + 60*yd1**2 - 360*xd1*yd1**2 + 360*xd1**2*yd1**2
case(67)
 s=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 60*yd1 - 360*xd1*yd1 + 360*xd1**2*yd1 - 60*yd1**2 + 360*xd1*yd1**2 - 360*xd1**2*yd1**2 &
+ 24.d0*zd1 - 144.d0*xd1*zd1 + 144.d0*xd1**2*zd1 - 120*yd1*zd1 + 720*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 + 120*yd1**2*zd1 &
- 720*xd1*yd1**2*zd1 + 720*xd1**2*yd1**2*zd1
case(66)
 s=-20.0d0+ 120*xd1 - 120*xd1**2 + 180*yd1 - 1080*xd1*yd1 + 1080*xd1**2*yd1 - 420*yd1**2 + 2520*xd1*yd1**2 - 2520*xd1**2*yd1**2&
 + 280*yd1**3 - 1680*xd1*yd1**3 + 1680*xd1**2*yd1**3

case(21)
 s=-2 + 24.d0*xd1 - 60*xd1**2 + 40*xd1**3
case(40)
 s=2 - 24.d0*xd1 + 60*xd1**2 - 40*xd1**3 - 4.d0*zd1 + 48.d0*xd1*zd1 - 120*xd1**2*zd1 + 80*xd1**3*zd1
case(64)
 s=-2 + 24.d0*xd1 - 60*xd1**2 + 40*xd1**3 + 12.d0*zd1 - 144.d0*xd1*zd1 + 360*xd1**2*zd1 - 240*xd1**3*zd1 - 12.d0*zd1**2 + 144.d0*xd1*zd1**2&
 - 360*xd1**2*zd1**2 + 240*xd1**3*zd1**2
case(38)
 s=6 - 72.d0*xd1 + 180*xd1**2 - 120*xd1**3 - 12.d0*yd1 + 144.d0*xd1*yd1 - 360*xd1**2*yd1 + 240*xd1**3*yd1
case(63)
 s=-6 + 72.d0*xd1 - 180*xd1**2 + 120*xd1**3 + 12.d0*yd1 - 144.d0*xd1*yd1 + 360*xd1**2*yd1 - 240*xd1**3*yd1 + 12.d0*zd1 - 144.d0*xd1*zd1 &
+ 360*xd1**2*zd1 - 240*xd1**3*zd1 - 24.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 + 480*xd1**3*yd1*zd1
case(62)
 s=-12 + 144.d0*xd1 - 360*xd1**2 + 240*xd1**3 + 60*yd1 - 720*xd1*yd1 + 1800*xd1**2*yd1 - 1200*xd1**3*yd1 - 60*yd1**2 &
+ 720*xd1*yd1**2 - 1800*xd1**2*yd1**2 + 1200*xd1**3*yd1**2

case(36)
 s=2 - 40*xd1 + 180*xd1**2 - 280*xd1**3 + 140*xd1**4
case(60)
 s=-2 + 40*xd1 - 180*xd1**2 + 280*xd1**3 - 140*xd1**4 + 4.d0*zd1 - 80*xd1*zd1 + 360*xd1**2*zd1 - 560*xd1**3*zd1 + 280*xd1**4*zd1
case(59)
 s=-6 + 120*xd1 - 540*xd1**2 + 840*xd1**3 - 420*xd1**4 + 12.d0*yd1 - 240*xd1*yd1 + 1080*xd1**2*yd1 - 1680*xd1**3*yd1 + 840*xd1**4*yd1

case(57)
 s=-2 + 60*xd1 - 420*xd1**2 + 1120*xd1**3 - 1260*xd1**4 + 504.d0*xd1**5
case default
 s=0.0d0
  end select
 dly=s    
	end function dly
	
real function dlz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=2
case(9)
 s=-6 + 12.d0*zd1
case(16)
 s=12 - 60*zd1 + 60*zd1**2
case(33)
 s=-20.0d0+ 180*zd1 - 420*zd1**2 + 280*zd1**3
case(54)
 s=30.0d0- 420*zd1 + 1680*zd1**2 - 2520*zd1**3 + 1260*zd1**4
case(83)
 s=-42 + 840*zd1 - 5040*zd1**2 + 12600*zd1**3 - 13860*zd1**4 + 5544.d0*zd1**5

case(8)
 s=-2 + 4.d0*yd1
case(15)
 s=6 - 12.d0*yd1 - 12.d0*zd1 + 24.d0*yd1*zd1
case(32)
 s=-12 + 24.d0*yd1 + 60*zd1 - 120*yd1*zd1 - 60*zd1**2 + 120*yd1*zd1**2
case(53)
 s=20.0d0- 40*yd1 - 180*zd1 + 360*yd1*zd1 + 420*zd1**2 - 840*yd1*zd1**2 - 280*zd1**3 + 560*yd1*zd1**3
case(82)
 s=-30.0d0+ 60*yd1 + 420*zd1 - 840*yd1*zd1 - 1680*zd1**2 + 3360*yd1*zd1**2 + 2520*zd1**3 - 5040*yd1*zd1**3 - 1260*zd1**4 &
+ 2520*yd1*zd1**4

case(14)
 s=2 - 12.d0*yd1 + 12.d0*yd1**2
case(31)
 s=-6 + 36.d0*yd1 - 36.d0*yd1**2 + 12.d0*zd1 - 72.d0*yd1*zd1 + 72.d0*yd1**2*zd1
case(52)
 s=12 - 72.d0*yd1 + 72.d0*yd1**2 - 60*zd1 + 360*yd1*zd1 - 360*yd1**2*zd1 + 60*zd1**2 - 360*yd1*zd1**2 + 360*yd1**2*zd1**2
case(81)
 s=-20.0d0+ 120*yd1 - 120*yd1**2 + 180*zd1 - 1080*yd1*zd1 + 1080*yd1**2*zd1 - 420*zd1**2 + 2520*yd1*zd1**2 - 2520*yd1**2*zd1**2&
 + 280*zd1**3 - 1680*yd1*zd1**3 + 1680*yd1**2*zd1**3
case(13)
 
case(30)
 s=-2 + 24.d0*yd1 - 60*yd1**2 + 40*yd1**3
case(51)
 s=6 - 72.d0*yd1 + 180*yd1**2 - 120*yd1**3 - 12.d0*zd1 + 144.d0*yd1*zd1 - 360*yd1**2*zd1 + 240*yd1**3*zd1
case(80)
 s=-12 + 144.d0*yd1 - 360*yd1**2 + 240*yd1**3 + 60*zd1 - 720*yd1*zd1 + 1800*yd1**2*zd1 - 1200*yd1**3*zd1 - 60*zd1**2 &
+ 720*yd1*zd1**2 - 1800*yd1**2*zd1**2 + 1200*yd1**3*zd1**2

case(50)
 s=2 - 40*yd1 + 180*yd1**2 - 280*yd1**3 + 140*yd1**4
case(79)
 s=-6 + 120*yd1 - 540*yd1**2 + 840*yd1**3 - 420*yd1**4 + 12.d0*zd1 - 240*yd1*zd1 + 1080*yd1**2*zd1 - 1680*yd1**3*zd1 &
+ 840*yd1**4*zd1

case(78)
 s=-2 + 60*yd1 - 420*yd1**2 + 1120*yd1**3 - 1260*yd1**4 + 504.d0*yd1**5


case(6)
 s=-2 + 4.d0*xd1
case(19)
 s=6 - 12.d0*xd1 - 12.d0*zd1 + 24.d0*xd1*zd1
case(29)
 s=-12 + 24.d0*xd1 + 60*zd1 - 120*xd1*zd1 - 60*zd1**2 + 120*xd1*zd1**2
case(55)
 s=20.0d0- 40*xd1 - 180*zd1 + 360*xd1*zd1 + 420*zd1**2 - 840*xd1*zd1**2 - 280*zd1**3 + 560*xd1*zd1**3
case(76)
 s=-30.0d0+ 60*xd1 + 420*zd1 - 840*xd1*zd1 - 1680*zd1**2 + 3360*xd1*zd1**2 + 2520*zd1**3 - 5040*xd1*zd1**3 &
- 1260*zd1**4 + 2520*xd1*zd1**4

case(17)
 s=2 - 4.d0*xd1 - 4.d0*yd1 + 8.d0*xd1*yd1
case(28)
 s=-6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 12.d0*zd1 - 24.d0*xd1*zd1 - 24.d0*yd1*zd1 + 48.d0*xd1*yd1*zd1
case(47)
 s=12 - 24.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1 - 60*zd1 + 120*xd1*zd1 + 120*yd1*zd1 - 240*xd1*yd1*zd1 + 60*zd1**2 &
- 120*xd1*zd1**2 - 120*yd1*zd1**2 + 240*xd1*yd1*zd1**2
case(75)
 s=-20.0d0+ 40*xd1 + 40*yd1 - 80*xd1*yd1 + 180*zd1 - 360*xd1*zd1 - 360*yd1*zd1 + 720*xd1*yd1*zd1 - 420*zd1**2&
 + 840*xd1*zd1**2 + 840*yd1*zd1**2 - 1680*xd1*yd1*zd1**2 + 280*zd1**3 - 560*xd1*zd1**3 - 560*yd1*zd1**3&
 + 1120*xd1*yd1*zd1**3

case(27)
 s=-2 + 4.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 - 12.d0*yd1**2 + 24.d0*xd1*yd1**2
case(48)
 s=6 - 12.d0*xd1 - 36.d0*yd1 + 72.d0*xd1*yd1 + 36.d0*yd1**2 - 72.d0*xd1*yd1**2 - 12.d0*zd1 + 24.d0*xd1*zd1 + 72.d0*yd1*zd1 - 144.d0*xd1*yd1*zd1 &
- 72.d0*yd1**2*zd1 + 144.d0*xd1*yd1**2*zd1
case(74)
 s=-12 + 24.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 72.d0*yd1**2 + 144.d0*xd1*yd1**2 + 60*zd1 - 120*xd1*zd1 - 360*yd1*zd1 &
+ 720*xd1*yd1*zd1 + 360*yd1**2*zd1 - 720*xd1*yd1**2*zd1 - 60*zd1**2 + 120*xd1*zd1**2 + 360*yd1*zd1**2 &
- 720*xd1*yd1*zd1**2 - 360*yd1**2*zd1**2 + 720*xd1*yd1**2*zd1**2

case(46)
 s=2 - 4.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1 + 60*yd1**2 - 120*xd1*yd1**2 - 40*yd1**3 + 80*xd1*yd1**3
case(73)
 s=-6 + 12.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 180*yd1**2 + 360*xd1*yd1**2 + 120*yd1**3 - 240*xd1*yd1**3 + 12.d0*zd1 - 24.d0*xd1*zd1 &
- 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 + 360*yd1**2*zd1 - 720*xd1*yd1**2*zd1 - 240*yd1**3*zd1 + 480*xd1*yd1**3*zd1

case(72)
 s=-2 + 4.d0*xd1 + 40*yd1 - 80*xd1*yd1 - 180*yd1**2 + 360*xd1*yd1**2 + 280*yd1**3 - 560*xd1*yd1**3 - 140*yd1**4 &
+ 280*xd1*yd1**4


case(18)
 s=2 - 12.d0*xd1 + 12.d0*xd1**2
case(25)
 s=-6 + 36.d0*xd1 - 36.d0*xd1**2 + 12.d0*zd1 - 72.d0*xd1*zd1 + 72.d0*xd1**2*zd1
case(44)
 s=12 - 72.d0*xd1 + 72.d0*xd1**2 - 60*zd1 + 360*xd1*zd1 - 360*xd1**2*zd1 + 60*zd1**2 - 360*xd1*zd1**2 + 360*xd1**2*zd1**2
case(70)
 s=-20.0d0+ 120*xd1 - 120*xd1**2 + 180*zd1 - 1080*xd1*zd1 + 1080*xd1**2*zd1 - 420*zd1**2 + 2520*xd1*zd1**2 &
- 2520*xd1**2*zd1**2 + 280*zd1**3 - 1680*xd1*zd1**3 + 1680*xd1**2*zd1**3

case(24)
 s=-2 + 12.d0*xd1 - 12.d0*xd1**2 + 4.d0*yd1 - 24.d0*xd1*yd1 + 24.d0*xd1**2*yd1
case(43)
 s=6 - 36.d0*xd1 + 36.d0*xd1**2 - 12.d0*yd1 + 72.d0*xd1*yd1 - 72.d0*xd1**2*yd1 - 12.d0*zd1 + 72.d0*xd1*zd1 - 72.d0*xd1**2*zd1 + 24.d0*yd1*zd1 &
- 144.d0*xd1*yd1*zd1 + 144.d0*xd1**2*yd1*zd1
case(69)
 s=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*yd1 - 144.d0*xd1*yd1 + 144.d0*xd1**2*yd1 + 60*zd1 - 360*xd1*zd1 + 360*xd1**2*zd1 &
- 120*yd1*zd1 + 720*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 - 60*zd1**2 + 360*xd1*zd1**2 - 360*xd1**2*zd1**2 &
+ 120*yd1*zd1**2 - 720*xd1*yd1*zd1**2 + 720*xd1**2*yd1*zd1**2

case(42)
 s=2 - 12.d0*xd1 + 12.d0*xd1**2 - 12.d0*yd1 + 72.d0*xd1*yd1 - 72.d0*xd1**2*yd1 + 12.d0*yd1**2 - 72.d0*xd1*yd1**2 + 72.d0*xd1**2*yd1**2
case(68)
 s=-6 + 36.d0*xd1 - 36.d0*xd1**2 + 36.d0*yd1 - 216.d0*xd1*yd1 + 216.d0*xd1**2*yd1 - 36.d0*yd1**2 + 216.d0*xd1*yd1**2 &
- 216.d0*xd1**2*yd1**2 + 12.d0*zd1 - 72.d0*xd1*zd1 + 72.d0*xd1**2*zd1 - 72.d0*yd1*zd1 + 432.d0*xd1*yd1*zd1 - 432.d0*xd1**2*yd1*zd1 &
+ 72.d0*yd1**2*zd1 - 432.d0*xd1*yd1**2*zd1 + 432.d0*xd1**2*yd1**2*zd1

case(67)
 s=-2 + 12.d0*xd1 - 12.d0*xd1**2 + 24.d0*yd1 - 144.d0*xd1*yd1 + 144.d0*xd1**2*yd1 - 60*yd1**2 + 360*xd1*yd1**2 &
- 360*xd1**2*yd1**2 + 40*yd1**3 - 240*xd1*yd1**3 + 240*xd1**2*yd1**3

case(22)
 s=-2 + 24.d0*xd1 - 60*xd1**2 + 40*xd1**3
case(39)
 s=6 - 72.d0*xd1 + 180*xd1**2 - 120*xd1**3 - 12.d0*zd1 + 144.d0*xd1*zd1 - 360*xd1**2*zd1 + 240*xd1**3*zd1
case(65)
 s=-12 + 144.d0*xd1 - 360*xd1**2 + 240*xd1**3 + 60*zd1 - 720*xd1*zd1 + 1800*xd1**2*zd1 - 1200*xd1**3*zd1 &
- 60*zd1**2 + 720*xd1*zd1**2 - 1800*xd1**2*zd1**2 + 1200*xd1**3*zd1**2

case(40)
 s=2 - 24.d0*xd1 + 60*xd1**2 - 40*xd1**3 - 4.d0*yd1 + 48.d0*xd1*yd1 - 120*xd1**2*yd1 + 80*xd1**3*yd1
case(64)
 s=-6 + 72.d0*xd1 - 180*xd1**2 + 120*xd1**3 + 12.d0*yd1 - 144.d0*xd1*yd1 + 360*xd1**2*yd1 - 240*xd1**3*yd1 + 12.d0*zd1 &
- 144.d0*xd1*zd1 + 360*xd1**2*zd1 - 240*xd1**3*zd1 - 24.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 + 480*xd1**3*yd1*zd1

case(63)
 s=-2 + 24.d0*xd1 - 60*xd1**2 + 40*xd1**3 + 12.d0*yd1 - 144.d0*xd1*yd1 + 360*xd1**2*yd1 - 240*xd1**3*yd1 - 12.d0*yd1**2 &
+ 144.d0*xd1*yd1**2 - 360*xd1**2*yd1**2 + 240*xd1**3*yd1**2

case(37)
 s=2 - 40*xd1 + 180*xd1**2 - 280*xd1**3 + 140*xd1**4
case(61)
 s=-6 + 120*xd1 - 540*xd1**2 + 840*xd1**3 - 420*xd1**4 + 12.d0*zd1 - 240*xd1*zd1 + 1080*xd1**2*zd1 &
- 1680*xd1**3*zd1 + 840*xd1**4*zd1

case(60)
 s=-2 + 40*xd1 - 180*xd1**2 + 280*xd1**3 - 140*xd1**4 + 4.d0*yd1 - 80*xd1*yd1 + 360*xd1**2*yd1 &
- 560*xd1**3*yd1 + 280*xd1**4*yd1

case(58)
 s=-2 + 60*xd1 - 420*xd1**2 + 1120*xd1**3 - 1260*xd1**4 + 504.d0*xd1**5

case default
 s=0.0d0
  end select
 dlz=s    
	end function dlz
 
real function dlx2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=12
case(18)
 s=-12 + 24.d0*zd1
case(25)
 s=12 - 72.d0*zd1 + 72.d0*zd1**2
case(44)
 s=-12 + 144.d0*zd1 - 360*zd1**2 + 240*zd1**3
case(70)
 s=12 - 240*zd1 + 1080*zd1**2 - 1680*zd1**3 + 840*zd1**4
case(11)
 s=-12 + 24.d0*yd1
case(24)
 s=12 - 24.d0*yd1 - 24.d0*zd1 + 48.d0*yd1*zd1
case(43)
 s=-12 + 24.d0*yd1 + 72.d0*zd1 - 144.d0*yd1*zd1 - 72.d0*zd1**2 + 144.d0*yd1*zd1**2
case(69)
 s=12 - 24.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1 + 360*zd1**2 - 720*yd1*zd1**2 - 240*zd1**3 + 480*yd1*zd1**3
case(23)
 s=12 - 72.d0*yd1 + 72.d0*yd1**2
case(42)
 s=-12 + 72.d0*yd1 - 72.d0*yd1**2 + 24.d0*zd1 - 144.d0*yd1*zd1 + 144.d0*yd1**2*zd1
case(68)
 s=12 - 72.d0*yd1 + 72.d0*yd1**2 - 72.d0*zd1 + 432.d0*yd1*zd1 - 432.d0*yd1**2*zd1 + 72.d0*zd1**2 - 432.d0*yd1*zd1**2 + 432.d0*yd1**2*zd1**2
case(41)
 s=-12 + 144.d0*yd1 - 360*yd1**2 + 240*yd1**3
case(67)
 s=12 - 144.d0*yd1 + 360*yd1**2 - 240*yd1**3 - 24.d0*zd1 + 288.d0*yd1*zd1 - 720*yd1**2*zd1 + 480*yd1**3*zd1
case(66)
 s=12 - 240*yd1 + 1080*yd1**2 - 1680*yd1**3 + 840*yd1**4
case(10)
 s=-60.0d0+ 120*xd1
case(22)
 s=60.0d0- 120*xd1 - 120*zd1 + 240*xd1*zd1
case(39)
 s=-60.0d0+ 120*xd1 + 360*zd1 - 720*xd1*zd1 - 360*zd1**2 + 720*xd1*zd1**2
case(65)
 s=60.0d0- 120*xd1 - 720*zd1 + 1440*xd1*zd1 + 1800*zd1**2 - 3600*xd1*zd1**2 - 1200*zd1**3 + 2400*xd1*zd1**3
case(21)
 s=60.0d0- 120*xd1 - 120*yd1 + 240*xd1*yd1
case(40)
 s=-60.0d0+ 120*xd1 + 120*yd1 - 240*xd1*yd1 + 120*zd1 - 240*xd1*zd1 - 240*yd1*zd1 + 480*xd1*yd1*zd1
case(64)
 s=60.0d0- 120*xd1 - 120*yd1 + 240*xd1*yd1 - 360*zd1 + 720*xd1*zd1 + 720*yd1*zd1 - 1440*xd1*yd1*zd1 + 360*zd1**2 &
- 720*xd1*zd1**2 - 720*yd1*zd1**2 + 1440*xd1*yd1*zd1**2
case(38)
 s=-60.0d0+ 120*xd1 + 360*yd1 - 720*xd1*yd1 - 360*yd1**2 + 720*xd1*yd1**2
case(63)
 s=60.0d0- 120*xd1 - 360*yd1 + 720*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 120*zd1 + 240*xd1*zd1 + 720*yd1*zd1 &
- 1440*xd1*yd1*zd1 - 720*yd1**2*zd1 + 1440*xd1*yd1**2*zd1
case(62)
 s=60.0d0- 120*xd1 - 720*yd1 + 1440*xd1*yd1 + 1800*yd1**2 - 3600*xd1*yd1**2 - 1200*yd1**3 + 2400*xd1*yd1**3
case(20)
 s=180.0d0- 840*xd1 + 840*xd1**2
case(37)
 s=-180.0d0+ 840*xd1 - 840*xd1**2 + 360*zd1 - 1680*xd1*zd1 + 1680*xd1**2*zd1
case(61)
 s=180.0d0- 840*xd1 + 840*xd1**2 - 1080*zd1 + 5040*xd1*zd1 - 5040*xd1**2*zd1 + 1080*zd1**2 - 5040*xd1*zd1**2 + 5040*xd1**2*zd1**2
case(36)
 s=-180.0d0+ 840*xd1 - 840*xd1**2 + 360*yd1 - 1680*xd1*yd1 + 1680*xd1**2*yd1
case(60)
 s=180.0d0- 840*xd1 + 840*xd1**2 - 360*yd1 + 1680*xd1*yd1 - 1680*xd1**2*yd1 - 360*zd1 + 1680*xd1*zd1 - 1680*xd1**2*zd1 &
+ 720*yd1*zd1 - 3360*xd1*yd1*zd1 + 3360*xd1**2*yd1*zd1
case(59)
 s=180.0d0- 840*xd1 + 840*xd1**2 - 1080*yd1 + 5040*xd1*yd1 - 5040*xd1**2*yd1 + 1080*yd1**2 - 5040*xd1*yd1**2 + 5040*xd1**2*yd1**2
case(35)
 s=-420.0d0+ 3360*xd1 - 7560*xd1**2 + 5040*xd1**3
case(58)
 s=420.0d0- 3360*xd1 + 7560*xd1**2 - 5040*xd1**3 - 840*zd1 + 6720*xd1*zd1 - 15120*xd1**2*zd1 + 10080*xd1**3*zd1
case(57)
 s=420.0d0- 3360*xd1 + 7560*xd1**2 - 5040*xd1**3 - 840*yd1 + 6720*xd1*yd1 - 15120*xd1**2*yd1 + 10080*xd1**3*yd1
case(56)
 s=840.0d0- 10080*xd1 + 37800*xd1**2 - 55440*xd1**3 + 27720*xd1**4
  end select
 dlx2=s    
	end function dlx2
 
 
real function dly2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=12
case(14)
 s=-12 + 24.d0*zd1
case(31)
 s=12 - 72.d0*zd1 + 72.d0*zd1**2
case(52)
 s=-12 + 144.d0*zd1 - 360*zd1**2 + 240*zd1**3
case(81)
 s=12 - 240*zd1 + 1080*zd1**2 - 1680*zd1**3 + 840*zd1**4
case(13)
 s=-60.0d0+ 120*yd1
case(30)
 s=60.0d0- 120*yd1 - 120*zd1 + 240*yd1*zd1
case(51)
 s=-60.0d0+ 120*yd1 + 360*zd1 - 720*yd1*zd1 - 360*zd1**2 + 720*yd1*zd1**2
case(80)
 s=60.0d0- 120*yd1 - 720*zd1 + 1440*yd1*zd1 + 1800*zd1**2 - 3600*yd1*zd1**2 - 1200*zd1**3 + 2400*yd1*zd1**3
case(34)
 s=180.0d0- 840*yd1 + 840*yd1**2
case(50)
 s=-180.0d0+ 840*yd1 - 840*yd1**2 + 360*zd1 - 1680*yd1*zd1 + 1680*yd1**2*zd1
case(79)
 s=180.0d0- 840*yd1 + 840*yd1**2 - 1080*zd1 + 5040*yd1*zd1 - 5040*yd1**2*zd1 + 1080*zd1**2 - 5040*yd1*zd1**2 &
+ 5040*yd1**2*zd1**2
case(49)
 s=-420.0d0+ 3360*yd1 - 7560*yd1**2 + 5040*yd1**3
case(78)
 s=420.0d0- 3360*yd1 + 7560*yd1**2 - 5040*yd1**3 - 840*zd1 + 6720*yd1*zd1 - 15120*yd1**2*zd1 + 10080*yd1**3*zd1
case(77)
 s=840.0d0- 10080*yd1 + 37800*yd1**2 - 55440*yd1**3 + 27720*yd1**4
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=-12 + 24.d0*xd1
case(27)
 s=12 - 24.d0*xd1 - 24.d0*zd1 + 48.d0*xd1*zd1
case(48)
 s=-12 + 24.d0*xd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 72.d0*zd1**2 + 144.d0*xd1*zd1**2
case(74)
 s=12 - 24.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1 + 360*zd1**2 - 720*xd1*zd1**2 - 240*zd1**3 + 480*xd1*zd1**3
case(26)
 s=60.0d0- 120*xd1 - 120*yd1 + 240*xd1*yd1
case(46)
 s=-60.0d0+ 120*xd1 + 120*yd1 - 240*xd1*yd1 + 120*zd1 - 240*xd1*zd1 - 240*yd1*zd1 + 480*xd1*yd1*zd1
case(73)
 s=60.0d0- 120*xd1 - 120*yd1 + 240*xd1*yd1 - 360*zd1 + 720*xd1*zd1 + 720*yd1*zd1 - 1440*xd1*yd1*zd1 + 360*zd1**2 &
- 720*xd1*zd1**2 - 720*yd1*zd1**2 + 1440*xd1*yd1*zd1**2
case(45)
 s=-180.0d0+ 360*xd1 + 840*yd1 - 1680*xd1*yd1 - 840*yd1**2 + 1680*xd1*yd1**2
case(72)
 s=180.0d0- 360*xd1 - 840*yd1 + 1680*xd1*yd1 + 840*yd1**2 - 1680*xd1*yd1**2 - 360*zd1 + 720*xd1*zd1 + 1680*yd1*zd1 &
- 3360*xd1*yd1*zd1 - 1680*yd1**2*zd1 + 3360*xd1*yd1**2*zd1
case(71)
 s=420.0d0- 840*xd1 - 3360*yd1 + 6720*xd1*yd1 + 7560*yd1**2 - 15120*xd1*yd1**2 - 5040*yd1**3 + 10080*xd1*yd1**3
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=12 - 72.d0*xd1 + 72.d0*xd1**2
case(42)
 s=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*zd1 - 144.d0*xd1*zd1 + 144.d0*xd1**2*zd1
case(68)
 s=12 - 72.d0*xd1 + 72.d0*xd1**2 - 72.d0*zd1 + 432.d0*xd1*zd1 - 432.d0*xd1**2*zd1 + 72.d0*zd1**2 - 432.d0*xd1*zd1**2 &
+ 432.d0*xd1**2*zd1**2
case(41)
 s=-60.0d0+ 360*xd1 - 360*xd1**2 + 120*yd1 - 720*xd1*yd1 + 720*xd1**2*yd1
case(67)
 s=60.0d0- 360*xd1 + 360*xd1**2 - 120*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 - 120*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1&
 + 240*yd1*zd1 - 1440*xd1*yd1*zd1 + 1440*xd1**2*yd1*zd1
case(66)
 s=180.0d0- 1080*xd1 + 1080*xd1**2 - 840*yd1 + 5040*xd1*yd1 - 5040*xd1**2*yd1 + 840*yd1**2 - 5040*xd1*yd1**2 &
+ 5040*xd1**2*yd1**2
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=-12 + 144.d0*xd1 - 360*xd1**2 + 240*xd1**3
case(63)
 s=12 - 144.d0*xd1 + 360*xd1**2 - 240*xd1**3 - 24.d0*zd1 + 288.d0*xd1*zd1 - 720*xd1**2*zd1 + 480*xd1**3*zd1
case(62)
 s=60.0d0- 720*xd1 + 1800*xd1**2 - 1200*xd1**3 - 120*yd1 + 1440*xd1*yd1 - 3600*xd1**2*yd1 + 2400*xd1**3*yd1
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=12 - 240*xd1 + 1080*xd1**2 - 1680*xd1**3 + 840*xd1**4
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly2=s    
	end function dly2
 
real function dlz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=12
case(16)
 s=-60.0d0+ 120*zd1
case(33)
 s=180.0d0- 840*zd1 + 840*zd1**2
case(54)
 s=-420.0d0+ 3360*zd1 - 7560*zd1**2 + 5040*zd1**3
case(83)
 s=840.0d0- 10080*zd1 + 37800*zd1**2 - 55440*zd1**3 + 27720*zd1**4
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=-12 + 24.d0*yd1
case(32)
 s=60.0d0- 120*yd1 - 120*zd1 + 240*yd1*zd1
case(53)
 s=-180.0d0+ 360*yd1 + 840*zd1 - 1680*yd1*zd1 - 840*zd1**2 + 1680*yd1*zd1**2
case(82)
 s=420.0d0- 840*yd1 - 3360*zd1 + 6720*yd1*zd1 + 7560*zd1**2 - 15120*yd1*zd1**2 - 5040*zd1**3 + 10080*yd1*zd1**3
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=12 - 72.d0*yd1 + 72.d0*yd1**2
case(52)
 s=-60.0d0+ 360*yd1 - 360*yd1**2 + 120*zd1 - 720*yd1*zd1 + 720*yd1**2*zd1
case(81)
 s=180.0d0- 1080*yd1 + 1080*yd1**2 - 840*zd1 + 5040*yd1*zd1 - 5040*yd1**2*zd1 + 840*zd1**2 - 5040*yd1*zd1**2&
 + 5040*yd1**2*zd1**2
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=-12 + 144.d0*yd1 - 360*yd1**2 + 240*yd1**3
case(80)
 s=60.0d0- 720*yd1 + 1800*yd1**2 - 1200*yd1**3 - 120*zd1 + 1440*yd1*zd1 - 3600*yd1**2*zd1 + 2400*yd1**3*zd1
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=12 - 240*yd1 + 1080*yd1**2 - 1680*yd1**3 + 840*yd1**4
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=-12 + 24.d0*xd1
case(29)
 s=60.0d0- 120*xd1 - 120*zd1 + 240*xd1*zd1
case(55)
 s=-180.0d0+ 360*xd1 + 840*zd1 - 1680*xd1*zd1 - 840*zd1**2 + 1680*xd1*zd1**2
case(76)
 s=420.0d0- 840*xd1 - 3360*zd1 + 6720*xd1*zd1 + 7560*zd1**2 - 15120*xd1*zd1**2 - 5040*zd1**3 + 10080*xd1*zd1**3
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=12 - 24.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1
case(47)
 s=-60.0d0+ 120*xd1 + 120*yd1 - 240*xd1*yd1 + 120*zd1 - 240*xd1*zd1 - 240*yd1*zd1 + 480*xd1*yd1*zd1
case(75)
 s=180.0d0- 360*xd1 - 360*yd1 + 720*xd1*yd1 - 840*zd1 + 1680*xd1*zd1 + 1680*yd1*zd1 - 3360*xd1*yd1*zd1 + 840*zd1**2&
 - 1680*xd1*zd1**2 - 1680*yd1*zd1**2 + 3360*xd1*yd1*zd1**2
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=-12 + 24.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 72.d0*yd1**2 + 144.d0*xd1*yd1**2
case(74)
 s=60.0d0- 120*xd1 - 360*yd1 + 720*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 120*zd1 + 240*xd1*zd1 + 720*yd1*zd1 &
- 1440*xd1*yd1*zd1 - 720*yd1**2*zd1 + 1440*xd1*yd1**2*zd1
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=12 - 24.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 240*yd1**3 + 480*xd1*yd1**3
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=12 - 72.d0*xd1 + 72.d0*xd1**2
case(44)
 s=-60.0d0+ 360*xd1 - 360*xd1**2 + 120*zd1 - 720*xd1*zd1 + 720*xd1**2*zd1
case(70)
 s=180.0d0- 1080*xd1 + 1080*xd1**2 - 840*zd1 + 5040*xd1*zd1 - 5040*xd1**2*zd1 + 840*zd1**2 - 5040*xd1*zd1**2 &
+ 5040*xd1**2*zd1**2
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*yd1 - 144.d0*xd1*yd1 + 144.d0*xd1**2*yd1
case(69)
 s=60.0d0- 360*xd1 + 360*xd1**2 - 120*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 - 120*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1 &
+ 240*yd1*zd1 - 1440*xd1*yd1*zd1 + 1440*xd1**2*yd1*zd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=12 - 72.d0*xd1 + 72.d0*xd1**2 - 72.d0*yd1 + 432.d0*xd1*yd1 - 432.d0*xd1**2*yd1 + 72.d0*yd1**2 - 432.d0*xd1*yd1**2 + 432.d0*xd1**2*yd1**2
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=-12 + 144.d0*xd1 - 360*xd1**2 + 240*xd1**3
case(65)
 s=60.0d0- 720*xd1 + 1800*xd1**2 - 1200*xd1**3 - 120*zd1 + 1440*xd1*zd1 - 3600*xd1**2*zd1 + 2400*xd1**3*zd1
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=12 - 144.d0*xd1 + 360*xd1**2 - 240*xd1**3 - 24.d0*yd1 + 288.d0*xd1*yd1 - 720*xd1**2*yd1 + 480*xd1**3*yd1
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=12 - 240*xd1 + 1080*xd1**2 - 1680*xd1**3 + 840*xd1**4
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlz2=s    
	end function dlz2
 
real function dlxy(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=4
case(17)
 s=-4 + 8.d0*zd1
case(28)
 s=4 - 24.d0*zd1 + 24.d0*zd1**2
case(47)
 s=-4 + 48.d0*zd1 - 120*zd1**2 + 80*zd1**3
case(75)
 s=4 - 80*zd1 + 360*zd1**2 - 560*zd1**3 + 280*zd1**4
case(12)
 s=-12 + 24.d0*yd1
case(27)
 s=12 - 24.d0*yd1 - 24.d0*zd1 + 48.d0*yd1*zd1
case(48)
 s=-12 + 24.d0*yd1 + 72.d0*zd1 - 144.d0*yd1*zd1 - 72.d0*zd1**2 + 144.d0*yd1*zd1**2
case(74)
 s=12 - 24.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1 + 360*zd1**2 - 720*yd1*zd1**2 - 240*zd1**3 + 480*yd1*zd1**3
case(26)
 s=24 - 120*yd1 + 120*yd1**2
case(46)
 s=-24 + 120*yd1 - 120*yd1**2 + 48.d0*zd1 - 240*yd1*zd1 + 240*yd1**2*zd1
case(73)
 s=24 - 120*yd1 + 120*yd1**2 - 144.d0*zd1 + 720*yd1*zd1 - 720*yd1**2*zd1 + 144.d0*zd1**2 - 720*yd1*zd1**2 + 720*yd1**2*zd1**2
case(45)
 s=-40.0d0+ 360*yd1 - 840*yd1**2 + 560*yd1**3
case(72)
 s=40.0d0- 360*yd1 + 840*yd1**2 - 560*yd1**3 - 80*zd1 + 720*yd1*zd1 - 1680*yd1**2*zd1 + 1120*yd1**3*zd1
case(71)
 s=60.0d0- 840*yd1 + 3360*yd1**2 - 5040*yd1**3 + 2520*yd1**4
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=-12 + 24.d0*xd1
case(24)
 s=12 - 24.d0*xd1 - 24.d0*zd1 + 48.d0*xd1*zd1
case(43)
 s=-12 + 24.d0*xd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 72.d0*zd1**2 + 144.d0*xd1*zd1**2
case(69)
 s=12 - 24.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1 + 360*zd1**2 - 720*xd1*zd1**2 - 240*zd1**3 + 480*xd1*zd1**3
case(23)
 s=36 - 72.d0*xd1 - 72.d0*yd1 + 144.d0*xd1*yd1
case(42)
 s=-36 + 72.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1
case(68)
 s=36 - 72.d0*xd1 - 72.d0*yd1 + 144.d0*xd1*yd1 - 216.d0*zd1 + 432.d0*xd1*zd1 + 432.d0*yd1*zd1 - 864.d0*xd1*yd1*zd1 + 216.d0*zd1**2 &
- 432.d0*xd1*zd1**2 - 432.d0*yd1*zd1**2 + 864.d0*xd1*yd1*zd1**2
case(41)
 s=-72 + 144.d0*xd1 + 360*yd1 - 720*xd1*yd1 - 360*yd1**2 + 720*xd1*yd1**2
case(67)
 s=72 - 144.d0*xd1 - 360*yd1 + 720*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 144.d0*zd1 + 288.d0*xd1*zd1 + 720*yd1*zd1 &
- 1440*xd1*yd1*zd1 - 720*yd1**2*zd1 + 1440*xd1*yd1**2*zd1
case(66)
 s=120.0d0- 240*xd1 - 1080*yd1 + 2160*xd1*yd1 + 2520*yd1**2 - 5040*xd1*yd1**2 - 1680*yd1**3 + 3360*xd1*yd1**3
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=24 - 120*xd1 + 120*xd1**2
case(40)
 s=-24 + 120*xd1 - 120*xd1**2 + 48.d0*zd1 - 240*xd1*zd1 + 240*xd1**2*zd1
case(64)
 s=24 - 120*xd1 + 120*xd1**2 - 144.d0*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1 + 144.d0*zd1**2 - 720*xd1*zd1**2 &
+ 720*xd1**2*zd1**2
case(38)
 s=-72 + 360*xd1 - 360*xd1**2 + 144.d0*yd1 - 720*xd1*yd1 + 720*xd1**2*yd1
case(63)
 s=72 - 360*xd1 + 360*xd1**2 - 144.d0*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 - 144.d0*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1 &
+ 288.d0*yd1*zd1 - 1440*xd1*yd1*zd1 + 1440*xd1**2*yd1*zd1
case(62)
 s=144 - 720*xd1 + 720*xd1**2 - 720*yd1 + 3600*xd1*yd1 - 3600*xd1**2*yd1 + 720*yd1**2 - 3600*xd1*yd1**2 &
+ 3600*xd1**2*yd1**2
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=-40.0d0+ 360*xd1 - 840*xd1**2 + 560*xd1**3
case(60)
 s=40.0d0- 360*xd1 + 840*xd1**2 - 560*xd1**3 - 80*zd1 + 720*xd1*zd1 - 1680*xd1**2*zd1 + 1120*xd1**3*zd1
case(59)
 s=120.0d0- 1080*xd1 + 2520*xd1**2 - 1680*xd1**3 - 240*yd1 + 2160*xd1*yd1 - 5040*xd1**2*yd1 + 3360*xd1**3*yd1
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=60.0d0- 840*xd1 + 3360*xd1**2 - 5040*xd1**3 + 2520*xd1**4
case(56)
 s=0.0d0
 
  end select
 dlxy=s    
	end function dlxy
 
real function dlyz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=4
case(15)
 s=-12 + 24.d0*zd1
case(32)
 s=24 - 120*zd1 + 120*zd1**2
case(53)
 s=-40.0d0+ 360*zd1 - 840*zd1**2 + 560*zd1**3
case(82)
 s=60.0d0- 840*zd1 + 3360*zd1**2 - 5040*zd1**3 + 2520*zd1**4
case(7)
 s=0.0d0
case(14)
 s=-12 + 24.d0*yd1
case(31)
 s=36 - 72.d0*yd1 - 72.d0*zd1 + 144.d0*yd1*zd1
case(52)
 s=-72 + 144.d0*yd1 + 360*zd1 - 720*yd1*zd1 - 360*zd1**2 + 720*yd1*zd1**2
case(81)
 s=120.0d0- 240*yd1 - 1080*zd1 + 2160*yd1*zd1 + 2520*zd1**2 - 5040*yd1*zd1**2 - 1680*zd1**3 + 3360*yd1*zd1**3
case(13)
 s=0.0d0
case(30)
 s=24 - 120*yd1 + 120*yd1**2
case(51)
 s=-72 + 360*yd1 - 360*yd1**2 + 144.d0*zd1 - 720*yd1*zd1 + 720*yd1**2*zd1
case(80)
 s=144 - 720*yd1 + 720*yd1**2 - 720*zd1 + 3600*yd1*zd1 - 3600*yd1**2*zd1 + 720*zd1**2 - 3600*yd1*zd1**2 + 3600*yd1**2*zd1**2
case(34)
 s=0.0d0
case(50)
 s=-40.0d0+ 360*yd1 - 840*yd1**2 + 560*yd1**3
case(79)
 s=120.0d0- 1080*yd1 + 2520*yd1**2 - 1680*yd1**3 - 240*zd1 + 2160*yd1*zd1 - 5040*yd1**2*zd1 + 3360*yd1**3*zd1
case(49)
 s=0.0d0
case(78)
 s=60.0d0- 840*yd1 + 3360*yd1**2 - 5040*yd1**3 + 2520*yd1**4
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=-4 + 8.d0*xd1
case(28)
 s=12 - 24.d0*xd1 - 24.d0*zd1 + 48.d0*xd1*zd1
case(47)
 s=-24 + 48.d0*xd1 + 120*zd1 - 240*xd1*zd1 - 120*zd1**2 + 240*xd1*zd1**2
case(75)
 s=40.0d0- 80*xd1 - 360*zd1 + 720*xd1*zd1 + 840*zd1**2 - 1680*xd1*zd1**2 - 560*zd1**3 + 1120*xd1*zd1**3
case(12)
 s=0.0d0
case(27)
 s=12 - 24.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1
case(48)
 s=-36 + 72.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1
case(74)
 s=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1 - 360*zd1 + 720*xd1*zd1 + 720*yd1*zd1 - 1440*xd1*yd1*zd1 + 360*zd1**2 &
- 720*xd1*zd1**2 - 720*yd1*zd1**2 + 1440*xd1*yd1*zd1**2
case(26)
 s=0.0d0
case(46)
 s=-24 + 48.d0*xd1 + 120*yd1 - 240*xd1*yd1 - 120*yd1**2 + 240*xd1*yd1**2
case(73)
 s=72 - 144.d0*xd1 - 360*yd1 + 720*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 144.d0*zd1 + 288.d0*xd1*zd1 + 720*yd1*zd1&
 - 1440*xd1*yd1*zd1 - 720*yd1**2*zd1 + 1440*xd1*yd1**2*zd1
case(45)
 s=0.0d0
case(72)
 s=40.0d0- 80*xd1 - 360*yd1 + 720*xd1*yd1 + 840*yd1**2 - 1680*xd1*yd1**2 - 560*yd1**3 + 1120*xd1*yd1**3
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=4 - 24.d0*xd1 + 24.d0*xd1**2
case(43)
 s=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*zd1 - 144.d0*xd1*zd1 + 144.d0*xd1**2*zd1
case(69)
 s=24 - 144.d0*xd1 + 144.d0*xd1**2 - 120*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1 + 120*zd1**2 - 720*xd1*zd1**2 + 720*xd1**2*zd1**2
case(23)
 s=0.0d0
case(42)
 s=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*yd1 - 144.d0*xd1*yd1 + 144.d0*xd1**2*yd1
case(68)
 s=36 - 216.d0*xd1 + 216.d0*xd1**2 - 72.d0*yd1 + 432.d0*xd1*yd1 - 432.d0*xd1**2*yd1 - 72.d0*zd1 + 432.d0*xd1*zd1 - 432.d0*xd1**2*zd1 &
+ 144.d0*yd1*zd1 - 864.d0*xd1*yd1*zd1 + 864.d0*xd1**2*yd1*zd1
case(41)
 s=0.0d0
case(67)
 s=24 - 144.d0*xd1 + 144.d0*xd1**2 - 120*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 + 120*yd1**2 - 720*xd1*yd1**2 + 720*xd1**2*yd1**2
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=-4 + 48.d0*xd1 - 120*xd1**2 + 80*xd1**3
case(64)
 s=12 - 144.d0*xd1 + 360*xd1**2 - 240*xd1**3 - 24.d0*zd1 + 288.d0*xd1*zd1 - 720*xd1**2*zd1 + 480*xd1**3*zd1
case(38)
 s=0.0d0
case(63)
 s=12 - 144.d0*xd1 + 360*xd1**2 - 240*xd1**3 - 24.d0*yd1 + 288.d0*xd1*yd1 - 720*xd1**2*yd1 + 480*xd1**3*yd1
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=4 - 80*xd1 + 360*xd1**2 - 560*xd1**3 + 280*xd1**4
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlyz=s    
	end function dlyz
 
real function dlxz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=4
case(19)
 s=-12 + 24.d0*zd1
case(29)
 s=24 - 120*zd1 + 120*zd1**2
case(55)
 s=-40.0d0+ 360*zd1 - 840*zd1**2 + 560*zd1**3
case(76)
 s=60.0d0- 840*zd1 + 3360*zd1**2 - 5040*zd1**3 + 2520*zd1**4
case(5)
 s=0.0d0
case(17)
 s=-4 + 8.d0*yd1
case(28)
 s=12 - 24.d0*yd1 - 24.d0*zd1 + 48.d0*yd1*zd1
case(47)
 s=-24 + 48.d0*yd1 + 120*zd1 - 240*yd1*zd1 - 120*zd1**2 + 240*yd1*zd1**2
case(75)
 s=40.0d0- 80*yd1 - 360*zd1 + 720*yd1*zd1 + 840*zd1**2 - 1680*yd1*zd1**2 - 560*zd1**3 + 1120*yd1*zd1**3
case(12)
 s=0.0d0
case(27)
 s=4 - 24.d0*yd1 + 24.d0*yd1**2
case(48)
 s=-12 + 72.d0*yd1 - 72.d0*yd1**2 + 24.d0*zd1 - 144.d0*yd1*zd1 + 144.d0*yd1**2*zd1
case(74)
 s=24 - 144.d0*yd1 + 144.d0*yd1**2 - 120*zd1 + 720*yd1*zd1 - 720*yd1**2*zd1 + 120*zd1**2 - 720*yd1*zd1**2 + 720*yd1**2*zd1**2
case(26)
 s=0.0d0
case(46)
 s=-4 + 48.d0*yd1 - 120*yd1**2 + 80*yd1**3
case(73)
 s=12 - 144.d0*yd1 + 360*yd1**2 - 240*yd1**3 - 24.d0*zd1 + 288.d0*yd1*zd1 - 720*yd1**2*zd1 + 480*yd1**3*zd1
case(45)
 s=0.0d0
case(72)
 s=4 - 80*yd1 + 360*yd1**2 - 560*yd1**3 + 280*yd1**4
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=-12 + 24.d0*xd1
case(25)
 s=36 - 72.d0*xd1 - 72.d0*zd1 + 144.d0*xd1*zd1
case(44)
 s=-72 + 144.d0*xd1 + 360*zd1 - 720*xd1*zd1 - 360*zd1**2 + 720*xd1*zd1**2
case(70)
 s=120.0d0- 240*xd1 - 1080*zd1 + 2160*xd1*zd1 + 2520*zd1**2 - 5040*xd1*zd1**2 - 1680*zd1**3 + 3360*xd1*zd1**3
case(11)
 s=0.0d0
case(24)
 s=12 - 24.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1
case(43)
 s=-36 + 72.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1
case(69)
 s=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1 - 360*zd1 + 720*xd1*zd1 + 720*yd1*zd1 - 1440*xd1*yd1*zd1 + 360*zd1**2 &
- 720*xd1*zd1**2 - 720*yd1*zd1**2 + 1440*xd1*yd1*zd1**2
case(23)
 s=0.0d0
case(42)
 s=-12 + 24.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 72.d0*yd1**2 + 144.d0*xd1*yd1**2
case(68)
 s=36 - 72.d0*xd1 - 216.d0*yd1 + 432.d0*xd1*yd1 + 216.d0*yd1**2 - 432.d0*xd1*yd1**2 - 72.d0*zd1 + 144.d0*xd1*zd1 + 432.d0*yd1*zd1 &
- 864.d0*xd1*yd1*zd1 - 432.d0*yd1**2*zd1 + 864.d0*xd1*yd1**2*zd1
case(41)
 s=0.0d0
case(67)
 s=12 - 24.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 240*yd1**3 + 480*xd1*yd1**3
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=24 - 120*xd1 + 120*xd1**2
case(39)
 s=-72 + 360*xd1 - 360*xd1**2 + 144.d0*zd1 - 720*xd1*zd1 + 720*xd1**2*zd1
case(65)
 s=144 - 720*xd1 + 720*xd1**2 - 720*zd1 + 3600*xd1*zd1 - 3600*xd1**2*zd1 + 720*zd1**2 - 3600*xd1*zd1**2 &
+ 3600*xd1**2*zd1**2
case(21)
 s=0.0d0
case(40)
 s=-24 + 120*xd1 - 120*xd1**2 + 48.d0*yd1 - 240*xd1*yd1 + 240*xd1**2*yd1
case(64)
 s=72 - 360*xd1 + 360*xd1**2 - 144.d0*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 - 144.d0*zd1 + 720*xd1*zd1 &
- 720*xd1**2*zd1 + 288.d0*yd1*zd1 - 1440*xd1*yd1*zd1 + 1440*xd1**2*yd1*zd1
case(38)
 s=0.0d0
case(63)
 s=24 - 120*xd1 + 120*xd1**2 - 144.d0*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 + 144.d0*yd1**2 - 720*xd1*yd1**2 &
+ 720*xd1**2*yd1**2
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=-40.0d0+ 360*xd1 - 840*xd1**2 + 560*xd1**3
case(61)
 s=120.0d0- 1080*xd1 + 2520*xd1**2 - 1680*xd1**3 - 240*zd1 + 2160*xd1*zd1 - 5040*xd1**2*zd1 + 3360*xd1**3*zd1
case(36)
 s=0.0d0
case(60)
 s=40.0d0- 360*xd1 + 840*xd1**2 - 560*xd1**3 - 80*yd1 + 720*xd1*yd1 - 1680*xd1**2*yd1 + 1120*xd1**3*yd1
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=60.0d0- 840*xd1 + 3360*xd1**2 - 5040*xd1**3 + 2520*xd1**4
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxz=s    
	end function dlxz
 
 
real function dlx3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=120
case(22)
 s=-120.0d0+ 240*zd1
case(39)
 s=120.0d0- 720*zd1 + 720*zd1**2
case(65)
 s=-120.0d0+ 1440*zd1 - 3600*zd1**2 + 2400*zd1**3
case(21)
 s=-120.0d0+ 240*yd1
case(40)
 s=120.0d0- 240*yd1 - 240*zd1 + 480*yd1*zd1
case(64)
 s=-120.0d0+ 240*yd1 + 720*zd1 - 1440*yd1*zd1 - 720*zd1**2 + 1440*yd1*zd1**2
case(38)
 s=120.0d0- 720*yd1 + 720*yd1**2
case(63)
 s=-120.0d0+ 720*yd1 - 720*yd1**2 + 240*zd1 - 1440*yd1*zd1 + 1440*yd1**2*zd1
case(62)
 s=-120.0d0+ 1440*yd1 - 3600*yd1**2 + 2400*yd1**3
case(20)
 s=-840.0d0+ 1680*xd1
case(37)
 s=840.0d0- 1680*xd1 - 1680*zd1 + 3360*xd1*zd1
case(61)
 s=-840.0d0+ 1680*xd1 + 5040*zd1 - 10080*xd1*zd1 - 5040*zd1**2 + 10080*xd1*zd1**2
case(36)
 s=840.0d0- 1680*xd1 - 1680*yd1 + 3360*xd1*yd1
case(60)
 s=-840.0d0+ 1680*xd1 + 1680*yd1 - 3360*xd1*yd1 + 1680*zd1 - 3360*xd1*zd1 - 3360*yd1*zd1 + 6720*xd1*yd1*zd1
case(59)
 s=-840.0d0+ 1680*xd1 + 5040*yd1 - 10080*xd1*yd1 - 5040*yd1**2 + 10080*xd1*yd1**2
case(35)
 s=3360.0d0- 15120*xd1 + 15120*xd1**2
case(58)
 s=-3360.0d0+ 15120*xd1 - 15120*xd1**2 + 6720*zd1 - 30240*xd1*zd1 + 30240*xd1**2*zd1
case(57)
 s=-3360.0d0+ 15120*xd1 - 15120*xd1**2 + 6720*yd1 - 30240*xd1*yd1 + 30240*xd1**2*yd1
case(56)
 s=-10080.0d0+ 75600*xd1 - 166320*xd1**2 + 110880*xd1**3
  end select
 dlx3=s    
	end function dlx3
 
real function dlx2y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0

selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=24
case(24)
 s=-24 + 48.d0*zd1
case(43)
 s=24 - 144.d0*zd1 + 144.d0*zd1**2
case(69)
 s=-24 + 288.d0*zd1 - 720*zd1**2 + 480*zd1**3
case(23)
 s=-72 + 144.d0*yd1
case(42)
 s=72 - 144.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1
case(68)
 s=-72 + 144.d0*yd1 + 432.d0*zd1 - 864.d0*yd1*zd1 - 432.d0*zd1**2 + 864.d0*yd1*zd1**2
case(41)
 s=144 - 720*yd1 + 720*yd1**2
case(67)
 s=-144 + 720*yd1 - 720*yd1**2 + 288.d0*zd1 - 1440*yd1*zd1 + 1440*yd1**2*zd1
case(66)
 s=-240.0d0+ 2160*yd1 - 5040*yd1**2 + 3360*yd1**3
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=-120.0d0+ 240*xd1
case(40)
 s=120.0d0- 240*xd1 - 240*zd1 + 480*xd1*zd1
case(64)
 s=-120.0d0+ 240*xd1 + 720*zd1 - 1440*xd1*zd1 - 720*zd1**2 + 1440*xd1*zd1**2
case(38)
 s=360.0d0- 720*xd1 - 720*yd1 + 1440*xd1*yd1
case(63)
 s=-360.0d0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
case(62)
 s=-720.0d0+ 1440*xd1 + 3600*yd1 - 7200*xd1*yd1 - 3600*yd1**2 + 7200*xd1*yd1**2
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=360.0d0- 1680*xd1 + 1680*xd1**2
case(60)
 s=-360.0d0+ 1680*xd1 - 1680*xd1**2 + 720*zd1 - 3360*xd1*zd1 + 3360*xd1**2*zd1
case(59)
 s=-1080.0d0+ 5040*xd1 - 5040*xd1**2 + 2160*yd1 - 10080*xd1*yd1 + 10080*xd1**2*yd1
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=-840.0d0+ 6720*xd1 - 15120*xd1**2 + 10080*xd1**3
case(56)
 s=0.0d0
  end select
 dlx2y=s    
	end function dlx2y
 
real function dlxy2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=24
case(27)
 s=-24 + 48.d0*zd1
case(48)
 s=24 - 144.d0*zd1 + 144.d0*zd1**2
case(74)
 s=-24 + 288.d0*zd1 - 720*zd1**2 + 480*zd1**3
case(26)
 s=-120.0d0+ 240*yd1
case(46)
 s=120.0d0- 240*yd1 - 240*zd1 + 480*yd1*zd1
case(73)
 s=-120.0d0+ 240*yd1 + 720*zd1 - 1440*yd1*zd1 - 720*zd1**2 + 1440*yd1*zd1**2
case(45)
 s=360.0d0- 1680*yd1 + 1680*yd1**2
case(72)
 s=-360.0d0+ 1680*yd1 - 1680*yd1**2 + 720*zd1 - 3360*yd1*zd1 + 3360*yd1**2*zd1
case(71)
 s=-840.0d0+ 6720*yd1 - 15120*yd1**2 + 10080*yd1**3
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=-72 + 144.d0*xd1
case(42)
 s=72 - 144.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1
case(68)
 s=-72 + 144.d0*xd1 + 432.d0*zd1 - 864.d0*xd1*zd1 - 432.d0*zd1**2 + 864.d0*xd1*zd1**2
case(41)
 s=360.0d0- 720*xd1 - 720*yd1 + 1440*xd1*yd1
case(67)
 s=-360.0d0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
case(66)
 s=-1080.0d0+ 2160*xd1 + 5040*yd1 - 10080*xd1*yd1 - 5040*yd1**2 + 10080*xd1*yd1**2
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=144 - 720*xd1 + 720*xd1**2
case(63)
 s=-144 + 720*xd1 - 720*xd1**2 + 288.d0*zd1 - 1440*xd1*zd1 + 1440*xd1**2*zd1
case(62)
 s=-720.0d0+ 3600*xd1 - 3600*xd1**2 + 1440*yd1 - 7200*xd1*yd1 + 7200*xd1**2*yd1
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=-240.0d0+ 2160*xd1 - 5040*xd1**2 + 3360*xd1**3
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy2=s    
	end function dlxy2
 
 
real function dly3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=120
case(30)
 s=-120.0d0+ 240*zd1
case(51)
 s=120.0d0- 720*zd1 + 720*zd1**2
case(80)
 s=-120.0d0+ 1440*zd1 - 3600*zd1**2 + 2400*zd1**3
case(34)
 s=-840.0d0+ 1680*yd1
case(50)
 s=840.0d0- 1680*yd1 - 1680*zd1 + 3360*yd1*zd1
case(79)
 s=-840.0d0+ 1680*yd1 + 5040*zd1 - 10080*yd1*zd1 - 5040*zd1**2 + 10080*yd1*zd1**2
case(49)
 s=3360.0d0- 15120*yd1 + 15120*yd1**2
case(78)
 s=-3360.0d0+ 15120*yd1 - 15120*yd1**2 + 6720*zd1 - 30240*yd1*zd1 + 30240*yd1**2*zd1
case(77)
 s=-10080.0d0+ 75600*yd1 - 166320*yd1**2 + 110880*yd1**3
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=-120.0d0+ 240*xd1
case(46)
 s=120.0d0- 240*xd1 - 240*zd1 + 480*xd1*zd1
case(73)
 s=-120.0d0+ 240*xd1 + 720*zd1 - 1440*xd1*zd1 - 720*zd1**2 + 1440*xd1*zd1**2
case(45)
 s=840.0d0- 1680*xd1 - 1680*yd1 + 3360*xd1*yd1
case(72)
 s=-840.0d0+ 1680*xd1 + 1680*yd1 - 3360*xd1*yd1 + 1680*zd1 - 3360*xd1*zd1 - 3360*yd1*zd1 + 6720*xd1*yd1*zd1
case(71)
 s=-3360.0d0+ 6720*xd1 + 15120*yd1 - 30240*xd1*yd1 - 15120*yd1**2 + 30240*xd1*yd1**2
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=120.0d0- 720*xd1 + 720*xd1**2
case(67)
 s=-120.0d0+ 720*xd1 - 720*xd1**2 + 240*zd1 - 1440*xd1*zd1 + 1440*xd1**2*zd1
case(66)
 s=-840.0d0+ 5040*xd1 - 5040*xd1**2 + 1680*yd1 - 10080*xd1*yd1 + 10080*xd1**2*yd1
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=-120.0d0+ 1440*xd1 - 3600*xd1**2 + 2400*xd1**3
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly3=s    
	end function dly3
 
 
real function dlx2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=24
case(25)
 s=-72 + 144.d0*zd1
case(44)
 s=144 - 720*zd1 + 720*zd1**2
case(70)
 s=-240.0d0+ 2160*zd1 - 5040*zd1**2 + 3360*zd1**3
case(11)
 s=0.0d0
case(24)
 s=-24 + 48.d0*yd1
case(43)
 s=72 - 144.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1
case(69)
 s=-144 + 288.d0*yd1 + 720*zd1 - 1440*yd1*zd1 - 720*zd1**2 + 1440*yd1*zd1**2
case(23)
 s=0.0d0
case(42)
 s=24 - 144.d0*yd1 + 144.d0*yd1**2
case(68)
 s=-72 + 432.d0*yd1 - 432.d0*yd1**2 + 144.d0*zd1 - 864.d0*yd1*zd1 + 864.d0*yd1**2*zd1
case(41)
 s=0.0d0
case(67)
 s=-24 + 288.d0*yd1 - 720*yd1**2 + 480*yd1**3
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=-120.0d0+ 240*xd1
case(39)
 s=360.0d0- 720*xd1 - 720*zd1 + 1440*xd1*zd1
case(65)
 s=-720.0d0+ 1440*xd1 + 3600*zd1 - 7200*xd1*zd1 - 3600*zd1**2 + 7200*xd1*zd1**2
case(21)
 s=0.0d0
case(40)
 s=120.0d0- 240*xd1 - 240*yd1 + 480*xd1*yd1
case(64)
 s=-360.0d0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
case(38)
 s=0.0d0
case(63)
 s=-120.0d0+ 240*xd1 + 720*yd1 - 1440*xd1*yd1 - 720*yd1**2 + 1440*xd1*yd1**2
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=360.0d0- 1680*xd1 + 1680*xd1**2
case(61)
 s=-1080.0d0+ 5040*xd1 - 5040*xd1**2 + 2160*zd1 - 10080*xd1*zd1 + 10080*xd1**2*zd1
case(36)
 s=0.0d0
case(60)
 s=-360.0d0+ 1680*xd1 - 1680*xd1**2 + 720*yd1 - 3360*xd1*yd1 + 3360*xd1**2*yd1
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=-840.0d0+ 6720*xd1 - 15120*xd1**2 + 10080*xd1**3
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2z=s    
	end function dlx2z
 
real function dlxyz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=8
case(28)
 s=-24 + 48.d0*zd1
case(47)
 s=48 - 240*zd1 + 240*zd1**2
case(75)
 s=-80.0d0+ 720*zd1 - 1680*zd1**2 + 1120*zd1**3
case(12)
 s=0.0d0
case(27)
 s=-24 + 48.d0*yd1
case(48)
 s=72 - 144.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1
case(74)
 s=-144 + 288.d0*yd1 + 720*zd1 - 1440*yd1*zd1 - 720*zd1**2 + 1440*yd1*zd1**2
case(26)
 s=0.0d0
case(46)
 s=48 - 240*yd1 + 240*yd1**2
case(73)
 s=-144 + 720*yd1 - 720*yd1**2 + 288.d0*zd1 - 1440*yd1*zd1 + 1440*yd1**2*zd1
case(45)
 s=0.0d0
case(72)
 s=-80.0d0+ 720*yd1 - 1680*yd1**2 + 1120*yd1**3
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=-24 + 48.d0*xd1
case(43)
 s=72 - 144.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1
case(69)
 s=-144 + 288.d0*xd1 + 720*zd1 - 1440*xd1*zd1 - 720*zd1**2 + 1440*xd1*zd1**2
case(23)
 s=0.0d0
case(42)
 s=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1
case(68)
 s=-216 + 432.d0*xd1 + 432.d0*yd1 - 864.d0*xd1*yd1 + 432.d0*zd1 - 864.d0*xd1*zd1 - 864.d0*yd1*zd1 + 1728.d0*xd1*yd1*zd1
case(41)
 s=0.0d0
case(67)
 s=-144 + 288.d0*xd1 + 720*yd1 - 1440*xd1*yd1 - 720*yd1**2 + 1440*xd1*yd1**2
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=48 - 240*xd1 + 240*xd1**2
case(64)
 s=-144 + 720*xd1 - 720*xd1**2 + 288.d0*zd1 - 1440*xd1*zd1 + 1440*xd1**2*zd1
case(38)
 s=0.0d0
case(63)
 s=-144 + 720*xd1 - 720*xd1**2 + 288.d0*yd1 - 1440*xd1*yd1 + 1440*xd1**2*yd1
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=-80.0d0+ 720*xd1 - 1680*xd1**2 + 1120*xd1**3
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxyz=s    
	end function dlxyz
 
 
real function dly2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=24
case(31)
 s=-72 + 144.d0*zd1
case(52)
 s=144 - 720*zd1 + 720*zd1**2
case(81)
 s=-240.0d0+ 2160*zd1 - 5040*zd1**2 + 3360*zd1**3
case(13)
 s=0.0d0
case(30)
 s=-120.0d0+ 240*yd1
case(51)
 s=360.0d0- 720*yd1 - 720*zd1 + 1440*yd1*zd1
case(80)
 s=-720.0d0+ 1440*yd1 + 3600*zd1 - 7200*yd1*zd1 - 3600*zd1**2 + 7200*yd1*zd1**2
case(34)
 s=0.0d0
case(50)
 s=360.0d0- 1680*yd1 + 1680*yd1**2
case(79)
 s=-1080.0d0+ 5040*yd1 - 5040*yd1**2 + 2160*zd1 - 10080*yd1*zd1 + 10080*yd1**2*zd1
case(49)
 s=0.0d0
case(78)
 s=-840.0d0+ 6720*yd1 - 15120*yd1**2 + 10080*yd1**3
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=-24 + 48.d0*xd1
case(48)
 s=72 - 144.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1
case(74)
 s=-144 + 288.d0*xd1 + 720*zd1 - 1440*xd1*zd1 - 720*zd1**2 + 1440*xd1*zd1**2
case(26)
 s=0.0d0
case(46)
 s=120.0d0- 240*xd1 - 240*yd1 + 480*xd1*yd1
case(73)
 s=-360.0d0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
case(45)
 s=0.0d0
case(72)
 s=-360.0d0+ 720*xd1 + 1680*yd1 - 3360*xd1*yd1 - 1680*yd1**2 + 3360*xd1*yd1**2
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=24 - 144.d0*xd1 + 144.d0*xd1**2
case(68)
 s=-72 + 432.d0*xd1 - 432.d0*xd1**2 + 144.d0*zd1 - 864.d0*xd1*zd1 + 864.d0*xd1**2*zd1
case(41)
 s=0.0d0
case(67)
 s=-120.0d0+ 720*xd1 - 720*xd1**2 + 240*yd1 - 1440*xd1*yd1 + 1440*xd1**2*yd1
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=-24 + 288.d0*xd1 - 720*xd1**2 + 480*xd1**3
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly2z=s    
	end function dly2z
 
 
real function dlxz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=24
case(29)
 s=-120.0d0+ 240*zd1
case(55)
 s=360.0d0- 1680*zd1 + 1680*zd1**2
case(76)
 s=-840.0d0+ 6720*zd1 - 15120*zd1**2 + 10080*zd1**3
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=-24 + 48.d0*yd1
case(47)
 s=120.0d0- 240*yd1 - 240*zd1 + 480*yd1*zd1
case(75)
 s=-360.0d0+ 720*yd1 + 1680*zd1 - 3360*yd1*zd1 - 1680*zd1**2 + 3360*yd1*zd1**2
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=24 - 144.d0*yd1 + 144.d0*yd1**2
case(74)
 s=-120.0d0+ 720*yd1 - 720*yd1**2 + 240*zd1 - 1440*yd1*zd1 + 1440*yd1**2*zd1
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=-24 + 288.d0*yd1 - 720*yd1**2 + 480*yd1**3
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=-72 + 144.d0*xd1
case(44)
 s=360.0d0- 720*xd1 - 720*zd1 + 1440*xd1*zd1
case(70)
 s=-1080.0d0+ 2160*xd1 + 5040*zd1 - 10080*xd1*zd1 - 5040*zd1**2 + 10080*xd1*zd1**2
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1
case(69)
 s=-360.0d0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=-72 + 144.d0*xd1 + 432.d0*yd1 - 864.d0*xd1*yd1 - 432.d0*yd1**2 + 864.d0*xd1*yd1**2
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=144 - 720*xd1 + 720*xd1**2
case(65)
 s=-720.0d0+ 3600*xd1 - 3600*xd1**2 + 1440*zd1 - 7200*xd1*zd1 + 7200*xd1**2*zd1
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=-144 + 720*xd1 - 720*xd1**2 + 288.d0*yd1 - 1440*xd1*yd1 + 1440*xd1**2*yd1
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=-240.0d0+ 2160*xd1 - 5040*xd1**2 + 3360*xd1**3
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxz2=s    
	end function dlxz2
 
 
real function dlyz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=24
case(32)
 s=-120.0d0+ 240*zd1
case(53)
 s=360.0d0- 1680*zd1 + 1680*zd1**2
case(82)
 s=-840.0d0+ 6720*zd1 - 15120*zd1**2 + 10080*zd1**3
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=-72 + 144.d0*yd1
case(52)
 s=360.0d0- 720*yd1 - 720*zd1 + 1440*yd1*zd1
case(81)
 s=-1080.0d0+ 2160*yd1 + 5040*zd1 - 10080*yd1*zd1 - 5040*zd1**2 + 10080*yd1*zd1**2
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=144 - 720*yd1 + 720*yd1**2
case(80)
 s=-720.0d0+ 3600*yd1 - 3600*yd1**2 + 1440*zd1 - 7200*yd1*zd1 + 7200*yd1**2*zd1
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=-240.0d0+ 2160*yd1 - 5040*yd1**2 + 3360*yd1**3
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=-24 + 48.d0*xd1
case(47)
 s=120.0d0- 240*xd1 - 240*zd1 + 480*xd1*zd1
case(75)
 s=-360.0d0+ 720*xd1 + 1680*zd1 - 3360*xd1*zd1 - 1680*zd1**2 + 3360*xd1*zd1**2
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1
case(74)
 s=-360.0d0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=-144 + 288.d0*xd1 + 720*yd1 - 1440*xd1*yd1 - 720*yd1**2 + 1440*xd1*yd1**2
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=24 - 144.d0*xd1 + 144.d0*xd1**2
case(69)
 s=-120.0d0+ 720*xd1 - 720*xd1**2 + 240*zd1 - 1440*xd1*zd1 + 1440*xd1**2*zd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=-72 + 432.d0*xd1 - 432.d0*xd1**2 + 144.d0*yd1 - 864.d0*xd1*yd1 + 864.d0*xd1**2*yd1
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=-24 + 288.d0*xd1 - 720*xd1**2 + 480*xd1**3
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlyz2=s    
	end function dlyz2
 
 
real function dlz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=120
case(33)
 s=-840.0d0+ 1680*zd1
case(54)
 s=3360.0d0- 15120*zd1 + 15120*zd1**2
case(83)
 s=-10080.0d0+ 75600*zd1 - 166320*zd1**2 + 110880*zd1**3
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=-120.0d0+ 240*yd1
case(53)
 s=840.0d0- 1680*yd1 - 1680*zd1 + 3360*yd1*zd1
case(82)
 s=-3360.0d0+ 6720*yd1 + 15120*zd1 - 30240*yd1*zd1 - 15120*zd1**2 + 30240*yd1*zd1**2
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=120.0d0- 720*yd1 + 720*yd1**2
case(81)
 s=-840.0d0+ 5040*yd1 - 5040*yd1**2 + 1680*zd1 - 10080*yd1*zd1 + 10080*yd1**2*zd1
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=-120.0d0+ 1440*yd1 - 3600*yd1**2 + 2400*yd1**3
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=-120.0d0+ 240*xd1
case(55)
 s=840.0d0- 1680*xd1 - 1680*zd1 + 3360*xd1*zd1
case(76)
 s=-3360.0d0+ 6720*xd1 + 15120*zd1 - 30240*xd1*zd1 - 15120*zd1**2 + 30240*xd1*zd1**2
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=120.0d0- 240*xd1 - 240*yd1 + 480*xd1*yd1
case(75)
 s=-840.0d0+ 1680*xd1 + 1680*yd1 - 3360*xd1*yd1 + 1680*zd1 - 3360*xd1*zd1 - 3360*yd1*zd1 + 6720*xd1*yd1*zd1
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=-120.0d0+ 240*xd1 + 720*yd1 - 1440*xd1*yd1 - 720*yd1**2 + 1440*xd1*yd1**2
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=120.0d0- 720*xd1 + 720*xd1**2
case(70)
 s=-840.0d0+ 5040*xd1 - 5040*xd1**2 + 1680*zd1 - 10080*xd1*zd1 + 10080*xd1**2*zd1
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=-120.0d0+ 720*xd1 - 720*xd1**2 + 240*yd1 - 1440*xd1*yd1 + 1440*xd1**2*yd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=-120.0d0+ 1440*xd1 - 3600*xd1**2 + 2400*xd1**3
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlz3=s    
	end function dlz3
 
 
 
real function dlx4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=1680
case(37)
 s=-1680.0d0+ 3360*zd1
case(61)
 s=1680.0d0- 10080*zd1 + 10080*zd1**2
case(36)
 s=-1680.0d0+ 3360*yd1
case(60)
 s=1680.0d0- 3360*yd1 - 3360*zd1 + 6720*yd1*zd1
case(59)
 s=1680.0d0- 10080*yd1 + 10080*yd1**2
case(35)
 s=-15120.0d0+ 30240*xd1
case(58)
 s=15120.0d0- 30240*xd1 - 30240*zd1 + 60480*xd1*zd1
case(57)
 s=15120.0d0- 30240*xd1 - 30240*yd1 + 60480*xd1*yd1
case(56)
 s=75600.0d0- 332640*xd1 + 332640*xd1**2
  end select
 dlx4=s    
	end function dlx4
 
 
real function dlx3y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=240
case(40)
 s=-240.0d0+ 480*zd1
case(64)
 s=240.0d0- 1440*zd1 + 1440*zd1**2
case(38)
 s=-720.0d0+ 1440*yd1
case(63)
 s=720.0d0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
case(62)
 s=1440.0d0- 7200*yd1 + 7200*yd1**2
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=-1680.0d0+ 3360*xd1
case(60)
 s=1680.0d0- 3360*xd1 - 3360*zd1 + 6720*xd1*zd1
case(59)
 s=5040.0d0- 10080*xd1 - 10080*yd1 + 20160*xd1*yd1
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=6720.0d0- 30240*xd1 + 30240*xd1**2
case(56)
 s=0.0d0
  end select
 dlx3y=s    
	end function dlx3y
 
 
 
 
real function dlx2y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=144
case(42)
 s=-144 + 288.d0*zd1
case(68)
 s=144 - 864.d0*zd1 + 864.d0*zd1**2
case(41)
 s=-720.0d0+ 1440*yd1
case(67)
 s=720.0d0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
case(66)
 s=2160.0d0- 10080*yd1 + 10080*yd1**2
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=-720.0d0+ 1440*xd1
case(63)
 s=720.0d0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
case(62)
 s=3600.0d0- 7200*xd1 - 7200*yd1 + 14400*xd1*yd1
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=2160.0d0- 10080*xd1 + 10080*xd1**2
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2y2=s    
	end function dlx2y2
 
 
 
real function dlxy3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=240
case(46)
 s=-240.0d0+ 480*zd1
case(73)
 s=240.0d0- 1440*zd1 + 1440*zd1**2
case(45)
 s=-1680.0d0+ 3360*yd1
case(72)
 s=1680.0d0- 3360*yd1 - 3360*zd1 + 6720*yd1*zd1
case(71)
 s=6720.0d0- 30240*yd1 + 30240*yd1**2
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=-720.0d0+ 1440*xd1
case(67)
 s=720.0d0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
case(66)
 s=5040.0d0- 10080*xd1 - 10080*yd1 + 20160*xd1*yd1
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=1440.0d0- 7200*xd1 + 7200*xd1**2
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy3=s    
	end function dlxy3
 
 
 
real function dly4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=1680
case(50)
 s=-1680.0d0+ 3360*zd1
case(79)
 s=1680.0d0- 10080*zd1 + 10080*zd1**2
case(49)
 s=-15120.0d0+ 30240*yd1
case(78)
 s=15120.0d0- 30240*yd1 - 30240*zd1 + 60480*yd1*zd1
case(77)
 s=75600.0d0- 332640*yd1 + 332640*yd1**2
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=-1680.0d0+ 3360*xd1
case(72)
 s=1680.0d0- 3360*xd1 - 3360*zd1 + 6720*xd1*zd1
case(71)
 s=15120.0d0- 30240*xd1 - 30240*yd1 + 60480*xd1*yd1
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=1680.0d0- 10080*xd1 + 10080*xd1**2
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly4=s    
	end function dly4
 
 
real function dlx3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=240
case(39)
 s=-720.0d0+ 1440*zd1
case(65)
 s=1440.0d0- 7200*zd1 + 7200*zd1**2
case(21)
 s=0.0d0
case(40)
 s=-240.0d0+ 480*yd1
case(64)
 s=720.0d0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
case(38)
 s=0.0d0
case(63)
 s=240.0d0- 1440*yd1 + 1440*yd1**2
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=-1680.0d0+ 3360*xd1
case(61)
 s=5040.0d0- 10080*xd1 - 10080*zd1 + 20160*xd1*zd1
case(36)
 s=0.0d0
case(60)
 s=1680.0d0- 3360*xd1 - 3360*yd1 + 6720*xd1*yd1
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=6720.0d0- 30240*xd1 + 30240*xd1**2
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx3z=s    
	end function dlx3z
 
 
real function dlx2yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=48
case(43)
 s=-144 + 288.d0*zd1
case(69)
 s=288 - 1440*zd1 + 1440*zd1**2
case(23)
 s=0.0d0
case(42)
 s=-144 + 288.d0*yd1
case(68)
 s=432 - 864.d0*yd1 - 864.d0*zd1 + 1728.d0*yd1*zd1
case(41)
 s=0.0d0
case(67)
 s=288 - 1440*yd1 + 1440*yd1**2
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=-240.0d0+ 480*xd1
case(64)
 s=720.0d0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
case(38)
 s=0.0d0
case(63)
 s=720.0d0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=720.0d0- 3360*xd1 + 3360*xd1**2
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2yz=s    
	end function dlx2yz
 
 
real function dlxy2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=48
case(48)
 s=-144 + 288.d0*zd1
case(74)
 s=288 - 1440*zd1 + 1440*zd1**2
case(26)
 s=0.0d0
case(46)
 s=-240.0d0+ 480*yd1
case(73)
 s=720.0d0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
case(45)
 s=0.0d0
case(72)
 s=720.0d0- 3360*yd1 + 3360*yd1**2
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=-144 + 288.d0*xd1
case(68)
 s=432 - 864.d0*xd1 - 864.d0*zd1 + 1728.d0*xd1*zd1
case(41)
 s=0.0d0
case(67)
 s=720.0d0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=288 - 1440*xd1 + 1440*xd1**2
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy2z=s    
	end function dlxy2z
 
 
 
 
real function dly3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=240
case(51)
 s=-720.0d0+ 1440*zd1
case(80)
 s=1440.0d0- 7200*zd1 + 7200*zd1**2
case(34)
 s=0.0d0
case(50)
 s=-1680.0d0+ 3360*yd1
case(79)
 s=5040.0d0- 10080*yd1 - 10080*zd1 + 20160*yd1*zd1
case(49)
 s=0.0d0
case(78)
 s=6720.0d0- 30240*yd1 + 30240*yd1**2
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=-240.0d0+ 480*xd1
case(73)
 s=720.0d0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
case(45)
 s=0.0d0
case(72)
 s=1680.0d0- 3360*xd1 - 3360*yd1 + 6720*xd1*yd1
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=240.0d0- 1440*xd1 + 1440*xd1**2
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly3z=s    
	end function dly3z
 
 
 
 
 
 
real function dlx2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=144
case(44)
 s=-720.0d0+ 1440*zd1
case(70)
 s=2160.0d0- 10080*zd1 + 10080*zd1**2
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=-144 + 288.d0*yd1
case(69)
 s=720.0d0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=144 - 864.d0*yd1 + 864.d0*yd1**2
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=-720.0d0+ 1440*xd1
case(65)
 s=3600.0d0- 7200*xd1 - 7200*zd1 + 14400*xd1*zd1
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=720.0d0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=2160.0d0- 10080*xd1 + 10080*xd1**2
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2z2=s    
	end function dlx2z2
 
 
 
real function dlxyz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=48
case(47)
 s=-240.0d0+ 480*zd1
case(75)
 s=720.0d0- 3360*zd1 + 3360*zd1**2
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=-144 + 288.d0*yd1
case(74)
 s=720.0d0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=288 - 1440*yd1 + 1440*yd1**2
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=-144 + 288.d0*xd1
case(69)
 s=720.0d0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=432 - 864.d0*xd1 - 864.d0*yd1 + 1728.d0*xd1*yd1
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=288 - 1440*xd1 + 1440*xd1**2
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxyz2=s    
	end function dlxyz2
 
 
real function dly2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=144
case(52)
 s=-720.0d0+ 1440*zd1
case(81)
 s=2160.0d0- 10080*zd1 + 10080*zd1**2
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=-720.0d0+ 1440*yd1
case(80)
 s=3600.0d0- 7200*yd1 - 7200*zd1 + 14400*yd1*zd1
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=2160.0d0- 10080*yd1 + 10080*yd1**2
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=-144 + 288.d0*xd1
case(74)
 s=720.0d0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=720.0d0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=144 - 864.d0*xd1 + 864.d0*xd1**2
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly2z2=s    
	end function dly2z2
 
 
 
real function dlxz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=240
case(55)
 s=-1680.0d0+ 3360*zd1
case(76)
 s=6720.0d0- 30240*zd1 + 30240*zd1**2
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=-240.0d0+ 480*yd1
case(75)
 s=1680.0d0- 3360*yd1 - 3360*zd1 + 6720*yd1*zd1
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=240.0d0- 1440*yd1 + 1440*yd1**2
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=-720.0d0+ 1440*xd1
case(70)
 s=5040.0d0- 10080*xd1 - 10080*zd1 + 20160*xd1*zd1
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=720.0d0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=1440.0d0- 7200*xd1 + 7200*xd1**2
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxz3=s    
	end function dlxz3
 
 
 
real function dlyz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=240
case(53)
 s=-1680.0d0+ 3360*zd1
case(82)
 s=6720.0d0- 30240*zd1 + 30240*zd1**2
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=-720.0d0+ 1440*yd1
case(81)
 s=5040.0d0- 10080*yd1 - 10080*zd1 + 20160*yd1*zd1
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=1440.0d0- 7200*yd1 + 7200*yd1**2
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=-240.0d0+ 480*xd1
case(75)
 s=1680.0d0- 3360*xd1 - 3360*zd1 + 6720*xd1*zd1
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=720.0d0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=240.0d0- 1440*xd1 + 1440*xd1**2
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlyz3=s    
	end function dlyz3
 
 
 
 
real function dlz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=1680
case(54)
 s=-15120.0d0+ 30240*zd1
case(83)
 s=75600.0d0- 332640*zd1 + 332640*zd1**2
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=-1680.0d0+ 3360*yd1
case(82)
 s=15120.0d0- 30240*yd1 - 30240*zd1 + 60480*yd1*zd1
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=1680.0d0- 10080*yd1 + 10080*yd1**2
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=-1680.0d0+ 3360*xd1
case(76)
 s=15120.0d0- 30240*xd1 - 30240*zd1 + 60480*xd1*zd1
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=1680.0d0- 3360*xd1 - 3360*yd1 + 6720*xd1*yd1
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=1680.0d0- 10080*xd1 + 10080*xd1**2
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlz4=s    
	end function dlz4
 
 
 
real function dlx5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=30240
case(58)
 s=-30240.0d0+ 60480*zd1
case(57)
 s=-30240.0d0+ 60480*yd1
case(56)
 s=-332640.0d0+ 665280*xd1
  end select
 dlx5=s    
	end function dlx5
 
real function dlx4y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=3360
case(60)
 s=-3360.0d0+ 6720*zd1
case(59)
 s=-10080.0d0+ 20160*yd1
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=-30240.0d0+ 60480*xd1
case(56)
 s=0.0d0
  end select
 dlx4y=s    
	end function dlx4y
 
 
 
 
real function dlx3y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=1440
case(63)
 s=-1440.0d0+ 2880*zd1
case(62)
 s=-7200.0d0+ 14400*yd1
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=-10080.0d0+ 20160*xd1
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx3y2=s    
	end function dlx3y2
 
 
 
real function dlx2y3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=1440
case(67)
 s=-1440.0d0+ 2880*zd1
case(66)
 s=-10080.0d0+ 20160*yd1
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=-7200.0d0+ 14400*xd1
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2y3=s    
	end function dlx2y3
 
 
 
 
real function dlxy4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=3360
case(72)
 s=-3360.0d0+ 6720*zd1
case(71)
 s=-30240.0d0+ 60480*yd1
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=-10080.0d0+ 20160*xd1
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy4=s    
	end function dlxy4
 
 
real function dly5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=30240
case(78)
 s=-30240.0d0+ 60480*zd1
case(77)
 s=-332640.0d0+ 665280*yd1
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=-30240.0d0+ 60480*xd1
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly5=s    
	end function dly5
 
real function dlx4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=3360
case(61)
 s=-10080.0d0+ 20160*zd1
case(36)
 s=0.0d0
case(60)
 s=-3360.0d0+ 6720*yd1
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=-30240.0d0+ 60480*xd1
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx4z=s    
	end function dlx4z
 
 
 
real function dlx3yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=480
case(64)
 s=-1440.0d0+ 2880*zd1
case(38)
 s=0.0d0
case(63)
 s=-1440.0d0+ 2880*yd1
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=-3360.0d0+ 6720*xd1
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx3yz=s    
	end function dlx3yz
 
 
 
real function dlx2y2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=288
case(68)
 s=-864 + 1728.d0*zd1
case(41)
 s=0.0d0
case(67)
 s=-1440.0d0+ 2880*yd1
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=-1440.0d0+ 2880*xd1
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2y2z=s    
	end function dlx2y2z
 
real function dlxy3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=480
case(73)
 s=-1440.0d0+ 2880*zd1
case(45)
 s=0.0d0
case(72)
 s=-3360.0d0+ 6720*yd1
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=-1440.0d0+ 2880*xd1
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy3z=s    
	end function dlxy3z
 
real function dly4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=3360
case(79)
 s=-10080.0d0+ 20160*zd1
case(49)
 s=0.0d0
case(78)
 s=-30240.0d0+ 60480*yd1
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=-3360.0d0+ 6720*xd1
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly4z=s    
	end function dly4z
 
 
real function dlx3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=1440
case(65)
 s=-7200.0d0+ 14400*zd1
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=-1440.0d0+ 2880*yd1
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=-10080.0d0+ 20160*xd1
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx3z2=s    
	end function dlx3z2
 
 
real function dlx2yz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=288
case(69)
 s=-1440.0d0+ 2880*zd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=-864 + 1728.d0*yd1
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=-1440.0d0+ 2880*xd1
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2yz2=s    
	end function dlx2yz2
	
 real function dlxy2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
 
 
 case(3)
 s=0.0d0

case(9)
 s=0.0d0

case(16)
 s=0.0d0

case(33)
 s=0.0d0

case(54)
 s=0.0d0

case(83)
 s=0.0d0

case(2)
 s=0.0d0

case(8)
 s=0.0d0

case(15)
 s=0.0d0

case(32)
 s=0.0d0

case(53)
 s=0.0d0

case(82)
 s=0.0d0

case(7)
 s=0.0d0

case(14)
 s=0.0d0

case(31)
 s=0.0d0

case(52)
 s=0.0d0

case(81)
 s=0.0d0

case(13)
 s=0.0d0

case(30)
 s=0.0d0

case(51)
 s=0.0d0

case(80)
 s=0.0d0

case(34)
 s=0.0d0

case(50)
 s=0.0d0

case(79)
 s=0.0d0

case(49)
 s=0.0d0

case(78)
 s=0.0d0

case(77)
 s=0.0d0

case(1)
 s=0.0d0

case(6)
 s=0.0d0

case(19)
 s=0.0d0

case(29)
 s=0.0d0

case(55)
 s=0.0d0

case(76)
 s=0.0d0

case(5)
 s=0.0d0

case(17)
 s=0.0d0

case(28)
 s=0.0d0

case(47)
 s=0.0d0

case(75)
 s=0.0d0

case(12)
 s=0.0d0

case(27)
 s=0.0d0

case(48)
 s=288

case(74)
 s=-1440.0d0+ 2880*zd1

case(26)
 s=0.0d0

case(46)
 s=0.0d0

case(73)
 s=-1440.0d0+ 2880*yd1

case(45)
 s=0.0d0

case(72)
 s=0.0d0

case(71)
 s=0.0d0

case(4)
 s=0.0d0

case(18)
 s=0.0d0

case(25)
 s=0.0d0

case(44)
 s=0.0d0

case(70)
 s=0.0d0

case(11)
 s=0.0d0

case(24)
 s=0.0d0

case(43)
 s=0.0d0

case(69)
 s=0.0d0

case(23)
 s=0.0d0

case(42)
 s=0.0d0

case(68)
 s=-864 + 1728.d0*xd1

case(41)
 s=0.0d0

case(67)
 s=0.0d0

case(66)
 s=0.0d0

case(10)
 s=0.0d0

case(22)
 s=0.0d0

case(39)
 s=0.0d0

case(65)
 s=0.0d0

case(21)
 s=0.0d0

case(40)
 s=0.0d0

case(64)
 s=0.0d0

case(38)
 s=0.0d0

case(63)
 s=0.0d0

case(62)
 s=0.0d0

case(20)
 s=0.0d0

case(37)
 s=0.0d0

case(61)
 s=0.0d0

case(36)
 s=0.0d0

case(60)
 s=0.0d0

case(59)
 s=0.0d0

case(35)
 s=0.0d0

case(58)
 s=0.0d0

case(57)
 s=0.0d0

case(56)
 s=0.0d0
 
   end select
 dlxy2z2=s    
	end function dlxy2z2
 
 
 
 
 
 
 
 
 
 
 
 
 
real function dly3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=1440
case(80)
 s=-7200.0d0+ 14400*zd1
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=-10080.0d0+ 20160*yd1
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=-1440.0d0+ 2880*xd1
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly3z2=s    
	end function dly3z2
 
 
 
real function dlx2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=1440
case(70)
 s=-10080.0d0+ 20160*zd1
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=-1440.0d0+ 2880*yd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=-7200.0d0+ 14400*xd1
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2z3=s    
	end function dlx2z3
 
 
real function dlxyz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=480
case(75)
 s=-3360.0d0+ 6720*zd1
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=-1440.0d0+ 2880*yd1
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=-1440.0d0+ 2880*xd1
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxyz3=s    
	end function dlxyz3
 
 
 
real function dly2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=1440
case(81)
 s=-10080.0d0+ 20160*zd1
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=-7200.0d0+ 14400*yd1
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=-1440.0d0+ 2880*xd1
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly2z3=s    
	end function dly2z3
 
 
real function dlxz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=3360
case(76)
 s=-30240.0d0+ 60480*zd1
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=-3360.0d0+ 6720*yd1
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=-10080.0d0+ 20160*xd1
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxz4=s    
	end function dlxz4
 
 
real function dlyz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=3360
case(82)
 s=-30240.0d0+ 60480*zd1
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=-10080.0d0+ 20160*yd1
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=-3360.0d0+ 6720*xd1
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlyz4=s    
	end function dlyz4
 
 
real function dlz5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=30240
case(83)
 s=-332640.0d0+ 665280*zd1
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=-30240.0d0+ 60480*yd1
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=-30240.0d0+ 60480*xd1
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlz5=s    
	end function dlz5
 
 
real function dlx6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=665280
  end select
 dlx6=s    
	end function dlx6
 
 
 
real function dlx5y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=60480
case(56)
 s=0.0d0
  end select
 dlx5y=s    
	end function dlx5y
 
 
real function dlx4y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=20160
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx4y2=s    
	end function dlx4y2
 
 
real function dlx3y3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=14400
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx3y3=s    
	end function dlx3y3
 
 
real function dlx2y4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=20160
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2y4=s    
	end function dlx2y4
 
 
 
real function dlxy5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=60480
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy5=s    
	end function dlxy5
 
 
real function dly6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=665280
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly6=s    
	end function dly6
 
 
real function dlx5z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=60480
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx5z=s    
	end function dlx5z
 
 
 
 
real function dlx4yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=6720
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx4yz=s    
	end function dlx4yz
 
 
 
real function dlx3y2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=2880
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx3y2z=s    
	end function dlx3y2z
 
 
 
real function dlx2y3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=2880
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2y3z=s    
	end function dlx2y3z
 
 
 
 
 
 
real function dlxy4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=6720
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy4z=s    
	end function dlxy4z
 
 
 
real function dly5z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=60480
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly5z=s    
	end function dly5z
 
 
 
real function dlx4z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=20160
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx4z2=s    
	end function dlx4z2
 
 
real function dlx3yz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=2880
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx3yz2=s    
	end function dlx3yz2
 
 
 
real function dlx2y2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=1728
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2y2z2=s    
	end function dlx2y2z2
 
 
real function dlxy3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=2880
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy3z2=s    
	end function dlxy3z2
 
 
 
 
real function dly4z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=20160
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly4z2=s    
	end function dly4z2
 
 
real function dlx3z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=14400
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx3z3=s    
	end function dlx3z3
 
 
real function dlx2yz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=2880
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2yz3=s    
	end function dlx2yz3
 
 
 
 
 
real function dlxy2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=2880
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlxy2z3=s    
	end function dlxy2z3
 
 
 
 
real function dly3z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=14400
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dly3z3=s    
	end function dly3z3
 
 
 
real function dlx2z4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=20160
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
  end select
 dlx2z4=s    
	end function dlx2z4
 
 
real function dlxyz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=6720
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
 end select
 dlxyz4=s    
	end function dlxyz4
 
 
real function dly2z4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=20160
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
 end select  
	dly2z4=s    
	end function dly2z4 
 
 
 
real function dlxz5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=60480
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
 end select  
	dlxz5=s    
	end function dlxz5 
 
 
 
real function dlyz5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=0.0d0
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=60480
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
 end select  
	dlyz5=s    
	end function dlyz5
 
 
real function dlz6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s
                    
s=0.0d0
selectcase (nderivative)
case(3)
 s=0.0d0
case(9)
 s=0.0d0
case(16)
 s=0.0d0
case(33)
 s=0.0d0
case(54)
 s=0.0d0
case(83)
 s=665280
case(2)
 s=0.0d0
case(8)
 s=0.0d0
case(15)
 s=0.0d0
case(32)
 s=0.0d0
case(53)
 s=0.0d0
case(82)
 s=0.0d0
case(7)
 s=0.0d0
case(14)
 s=0.0d0
case(31)
 s=0.0d0
case(52)
 s=0.0d0
case(81)
 s=0.0d0
case(13)
 s=0.0d0
case(30)
 s=0.0d0
case(51)
 s=0.0d0
case(80)
 s=0.0d0
case(34)
 s=0.0d0
case(50)
 s=0.0d0
case(79)
 s=0.0d0
case(49)
 s=0.0d0
case(78)
 s=0.0d0
case(77)
 s=0.0d0
case(1)
 s=0.0d0
case(6)
 s=0.0d0
case(19)
 s=0.0d0
case(29)
 s=0.0d0
case(55)
 s=0.0d0
case(76)
 s=0.0d0
case(5)
 s=0.0d0
case(17)
 s=0.0d0
case(28)
 s=0.0d0
case(47)
 s=0.0d0
case(75)
 s=0.0d0
case(12)
 s=0.0d0
case(27)
 s=0.0d0
case(48)
 s=0.0d0
case(74)
 s=0.0d0
case(26)
 s=0.0d0
case(46)
 s=0.0d0
case(73)
 s=0.0d0
case(45)
 s=0.0d0
case(72)
 s=0.0d0
case(71)
 s=0.0d0
case(4)
 s=0.0d0
case(18)
 s=0.0d0
case(25)
 s=0.0d0
case(44)
 s=0.0d0
case(70)
 s=0.0d0
case(11)
 s=0.0d0
case(24)
 s=0.0d0
case(43)
 s=0.0d0
case(69)
 s=0.0d0
case(23)
 s=0.0d0
case(42)
 s=0.0d0
case(68)
 s=0.0d0
case(41)
 s=0.0d0
case(67)
 s=0.0d0
case(66)
 s=0.0d0
case(10)
 s=0.0d0
case(22)
 s=0.0d0
case(39)
 s=0.0d0
case(65)
 s=0.0d0
case(21)
 s=0.0d0
case(40)
 s=0.0d0
case(64)
 s=0.0d0
case(38)
 s=0.0d0
case(63)
 s=0.0d0
case(62)
 s=0.0d0
case(20)
 s=0.0d0
case(37)
 s=0.0d0
case(61)
 s=0.0d0
case(36)
 s=0.0d0
case(60)
 s=0.0d0
case(59)
 s=0.0d0
case(35)
 s=0.0d0
case(58)
 s=0.0d0
case(57)
 s=0.0d0
case(56)
 s=0.0d0
 end select  
	dlz6=s    
	end function dlz6
	
	
real function tl3dz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
  case(3) ; s= 1/hxc
  case(9) ; s= zd1/hxc**2
  case(19) ; s= zd1**2/(2.0d0*hxc**3)
  case(34) ; s= zd1**3/(6.0d0*hxc**4)
  case(55) ; s= zd1**4/(24.0d0*hxc**5)
  case(83) ; s= zd1**5/(120.0d0*hxc**6)
  case(2) ; s= 0.0d0
  case(8) ; s= yd1/hxc**2
  case(18) ; s= (yd1*zd1)/hxc**3
  case(33) ; s= (yd1*zd1**2)/(2.0d0*hxc**4)
  case(54) ; s= (yd1*zd1**3)/(6.0d0*hxc**5)
  case(82) ; s= (yd1*zd1**4)/(24.0d0*hxc**6)
  case(7) ; s= 0.0d0
  case(17) ; s= yd1**2/(2.0d0*hxc**3)
  case(32) ; s= (yd1**2*zd1)/(2.0d0*hxc**4)
  case(53) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**5)
  case(81) ; s= (yd1**2*zd1**3)/(12.0d0*hxc**6)
  case(16) ; s= 0.0d0
  case(31) ; s= yd1**3/(6.0d0*hxc**4)
  case(52) ; s= (yd1**3*zd1)/(6.0d0*hxc**5)
  case(80) ; s= (yd1**3*zd1**2)/(12.0d0*hxc**6)
  case(30) ; s= 0.0d0
  case(51) ; s= yd1**4/(24.0d0*hxc**5)
  case(79) ; s= (yd1**4*zd1)/(24.0d0*hxc**6)
  case(50) ; s= 0.0d0
  case(78) ; s= yd1**5/(120.0d0*hxc**6)
  case(77) ; s= 0.0d0
  case(1) ; s= 0.0d0
  case(6) ; s= xd1/hxc**2
  case(14) ; s= (xd1*zd1)/hxc**3
  case(27) ; s= (xd1*zd1**2)/(2.0d0*hxc**4)
  case(49) ; s= (xd1*zd1**3)/(6.0d0*hxc**5)
  case(76) ; s= (xd1*zd1**4)/(24.0d0*hxc**6)
  case(5) ; s= 0.0d0
  case(15) ; s= (xd1*yd1)/hxc**3
  case(29) ; s= (xd1*yd1*zd1)/hxc**4
  case(48) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**5)
  case(75) ; s= (xd1*yd1*zd1**3)/(6.0d0*hxc**6)
  case(13) ; s= 0.0d0
  case(28) ; s= (xd1*yd1**2)/(2.0d0*hxc**4)
  case(47) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**5)
  case(74) ; s= (xd1*yd1**2*zd1**2)/(4.0d0*hxc**6)
  case(26) ; s= 0.0d0
  case(46) ; s= (xd1*yd1**3)/(6.0d0*hxc**5)
  case(73) ; s= (xd1*yd1**3*zd1)/(6.0d0*hxc**6)
  case(45) ; s= 0.0d0
  case(72) ; s= (xd1*yd1**4)/(24.0d0*hxc**6)
  case(71) ; s= 0.0d0
  case(4) ; s= 0.0d0
  case(12) ; s= xd1**2/(2.0d0*hxc**3)
  case(24) ; s= (xd1**2*zd1)/(2.0d0*hxc**4)
  case(44) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**5)
  case(70) ; s= (xd1**2*zd1**3)/(12.0d0*hxc**6)
  case(11) ; s= 0.0d0
  case(25) ; s= (xd1**2*yd1)/(2.0d0*hxc**4)
  case(43) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**5)
  case(69) ; s= (xd1**2*yd1*zd1**2)/(4.0d0*hxc**6)
  case(23) ; s= 0.0d0
  case(42) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**5)
  case(68) ; s= (xd1**2*yd1**2*zd1)/(4.0d0*hxc**6)
  case(41) ; s= 0.0d0
  case(67) ; s= (xd1**2*yd1**3)/(12.0d0*hxc**6)
  case(66) ; s= 0.0d0
  case(10) ; s= 0.0d0
  case(22) ; s= xd1**3/(6.0d0*hxc**4)
  case(39) ; s= (xd1**3*zd1)/(6.0d0*hxc**5)
  case(65) ; s= (xd1**3*zd1**2)/(12.0d0*hxc**6)
  case(21) ; s= 0.0d0
  case(40) ; s= (xd1**3*yd1)/(6.0d0*hxc**5)
  case(64) ; s= (xd1**3*yd1*zd1)/(6.0d0*hxc**6)
  case(38) ; s= 0.0d0
  case(63) ; s= (xd1**3*yd1**2)/(12.0d0*hxc**6)
  case(62) ; s= 0.0d0
  case(20) ; s= 0.0d0
  case(37) ; s= xd1**4/(24.0d0*hxc**5)
  case(61) ; s= (xd1**4*zd1)/(24.0d0*hxc**6)
  case(36) ; s= 0.0d0
  case(60) ; s= (xd1**4*yd1)/(24.0d0*hxc**6)
  case(59) ; s= 0.0d0
  case(35) ; s= 0.0d0
  case(58) ; s= xd1**5/(120.0d0*hxc**6)
  case(57) ; s= 0.0d0
  case(56) ; s= 0.0d0
 
 
 end select

    tl3dz = s

  end function
 
  
 
 
  real function tl3dz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
  case(3) ; s= 0.0d0
  case(9) ; s= hxc**(-2)
  case(19) ; s= zd1/hxc**3
  case(34) ; s= zd1**2/(2.0d0*hxc**4)
  case(55) ; s= zd1**3/(6.0d0*hxc**5)
  case(83) ; s= zd1**4/(24.0d0*hxc**6)
  case(2) ; s= 0.0d0
  case(8) ; s= 0.0d0
  case(18) ; s= yd1/hxc**3
  case(33) ; s= (yd1*zd1)/hxc**4
  case(54) ; s= (yd1*zd1**2)/(2.0d0*hxc**5)
  case(82) ; s= (yd1*zd1**3)/(6.0d0*hxc**6)
  case(7) ; s= 0.0d0
  case(17) ; s= 0.0d0
  case(32) ; s= yd1**2/(2.0d0*hxc**4)
  case(53) ; s= (yd1**2*zd1)/(2.0d0*hxc**5)
  case(81) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**6)
  case(16) ; s= 0.0d0
  case(31) ; s= 0.0d0
  case(52) ; s= yd1**3/(6.0d0*hxc**5)
  case(80) ; s= (yd1**3*zd1)/(6.0d0*hxc**6)
  case(30) ; s= 0.0d0
  case(51) ; s= 0.0d0
  case(79) ; s= yd1**4/(24.0d0*hxc**6)
  case(50) ; s= 0.0d0
  case(78) ; s= 0.0d0
  case(77) ; s= 0.0d0
  case(1) ; s= 0.0d0
  case(6) ; s= 0.0d0
  case(14) ; s= xd1/hxc**3
  case(27) ; s= (xd1*zd1)/hxc**4
  case(49) ; s= (xd1*zd1**2)/(2.0d0*hxc**5)
  case(76) ; s= (xd1*zd1**3)/(6.0d0*hxc**6)
  case(5) ; s= 0.0d0
  case(15) ; s= 0.0d0
  case(29) ; s= (xd1*yd1)/hxc**4
  case(48) ; s= (xd1*yd1*zd1)/hxc**5
  case(75) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**6)
  case(13) ; s= 0.0d0
  case(28) ; s= 0.0d0
  case(47) ; s= (xd1*yd1**2)/(2.0d0*hxc**5)
  case(74) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**6)
  case(26) ; s= 0.0d0
  case(46) ; s= 0.0d0
  case(73) ; s= (xd1*yd1**3)/(6.0d0*hxc**6)
  case(45) ; s= 0.0d0
  case(72) ; s= 0.0d0
  case(71) ; s= 0.0d0
  case(4) ; s= 0.0d0
  case(12) ; s= 0.0d0
  case(24) ; s= xd1**2/(2.0d0*hxc**4)
  case(44) ; s= (xd1**2*zd1)/(2.0d0*hxc**5)
  case(70) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**6)
  case(11) ; s= 0.0d0
  case(25) ; s= 0.0d0
  case(43) ; s= (xd1**2*yd1)/(2.0d0*hxc**5)
  case(69) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**6)
  case(23) ; s= 0.0d0
  case(42) ; s= 0.0d0
  case(68) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**6)
  case(41) ; s= 0.0d0
  case(67) ; s= 0.0d0
  case(66) ; s= 0.0d0
  case(10) ; s= 0.0d0
  case(22) ; s= 0.0d0
  case(39) ; s= xd1**3/(6.0d0*hxc**5)
  case(65) ; s= (xd1**3*zd1)/(6.0d0*hxc**6)
  case(21) ; s= 0.0d0
  case(40) ; s= 0.0d0
  case(64) ; s= (xd1**3*yd1)/(6.0d0*hxc**6)
  case(38) ; s= 0.0d0
  case(63) ; s= 0.0d0
  case(62) ; s= 0.0d0
  case(20) ; s= 0.0d0
  case(37) ; s= 0.0d0
  case(61) ; s= xd1**4/(24.0d0*hxc**6)
  case(36) ; s= 0.0d0
  case(60) ; s= 0.0d0
  case(59) ; s= 0.0d0
  case(35) ; s= 0.0d0
  case(58) ; s= 0.0d0
  case(57) ; s= 0.0d0
  case(56) ; s= 0.0d0
  
  
  
  
   end select

    tl3dz2 = s

  end function
 
 
 
  real function tl3dz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
  case(3) ; s= 0.0d0
  case(9) ; s= 0.0d0
  case(19) ; s= hxc**(-3)
  case(34) ; s= zd1/hxc**4
  case(55) ; s= zd1**2/(2.0d0*hxc**5)
  case(83) ; s= zd1**3/(6.0d0*hxc**6)
  case(2) ; s= 0.0d0
  case(8) ; s= 0.0d0
  case(18) ; s= 0.0d0
  case(33) ; s= yd1/hxc**4
  case(54) ; s= (yd1*zd1)/hxc**5
  case(82) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
  case(7) ; s= 0.0d0
  case(17) ; s= 0.0d0
  case(32) ; s= 0.0d0
  case(53) ; s= yd1**2/(2.0d0*hxc**5)
  case(81) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
  case(16) ; s= 0.0d0
  case(31) ; s= 0.0d0
  case(52) ; s= 0.0d0
  case(80) ; s= yd1**3/(6.0d0*hxc**6)
  case(30) ; s= 0.0d0
  case(51) ; s= 0.0d0
  case(79) ; s= 0.0d0
  case(50) ; s= 0.0d0
  case(78) ; s= 0.0d0
  case(77) ; s= 0.0d0
  case(1) ; s= 0.0d0
  case(6) ; s= 0.0d0
  case(14) ; s= 0.0d0
  case(27) ; s= xd1/hxc**4
  case(49) ; s= (xd1*zd1)/hxc**5
  case(76) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
  case(5) ; s= 0.0d0
  case(15) ; s= 0.0d0
  case(29) ; s= 0.0d0
  case(48) ; s= (xd1*yd1)/hxc**5
  case(75) ; s= (xd1*yd1*zd1)/hxc**6
  case(13) ; s= 0.0d0
  case(28) ; s= 0.0d0
  case(47) ; s= 0.0d0
  case(74) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
  case(26) ; s= 0.0d0
  case(46) ; s= 0.0d0
  case(73) ; s= 0.0d0
  case(45) ; s= 0.0d0
  case(72) ; s= 0.0d0
  case(71) ; s= 0.0d0
  case(4) ; s= 0.0d0
  case(12) ; s= 0.0d0
  case(24) ; s= 0.0d0
  case(44) ; s= xd1**2/(2.0d0*hxc**5)
  case(70) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
  case(11) ; s= 0.0d0
  case(25) ; s= 0.0d0
  case(43) ; s= 0.0d0
  case(69) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
  case(23) ; s= 0.0d0
  case(42) ; s= 0.0d0
  case(68) ; s= 0.0d0
  case(41) ; s= 0.0d0
  case(67) ; s= 0.0d0
  case(66) ; s= 0.0d0
  case(10) ; s= 0.0d0
  case(22) ; s= 0.0d0
  case(39) ; s= 0.0d0
  case(65) ; s= xd1**3/(6.0d0*hxc**6)
  case(21) ; s= 0.0d0
  case(40) ; s= 0.0d0
  case(64) ; s= 0.0d0
  case(38) ; s= 0.0d0
  case(63) ; s= 0.0d0
  case(62) ; s= 0.0d0
  case(20) ; s= 0.0d0
  case(37) ; s= 0.0d0
  case(61) ; s= 0.0d0
  case(36) ; s= 0.0d0
  case(60) ; s= 0.0d0
  case(59) ; s= 0.0d0
  case(35) ; s= 0.0d0
  case(58) ; s= 0.0d0
  case(57) ; s= 0.0d0
  case(56) ; s= 0.0d0
  
  
  
  
   end select

    tl3dz3 = s

  end function
 
 
 
  real function tl3dz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
  case(3) ; s= 0.0d0
  case(9) ; s= 0.0d0
  case(19) ; s= 0.0d0
  case(34) ; s= hxc**(-4)
  case(55) ; s= zd1/hxc**5
  case(83) ; s= zd1**2/(2.0d0*hxc**6)
  case(2) ; s= 0.0d0
  case(8) ; s= 0.0d0
  case(18) ; s= 0.0d0
  case(33) ; s= 0.0d0
  case(54) ; s= yd1/hxc**5
  case(82) ; s= (yd1*zd1)/hxc**6
  case(7) ; s= 0.0d0
  case(17) ; s= 0.0d0
  case(32) ; s= 0.0d0
  case(53) ; s= 0.0d0
  case(81) ; s= yd1**2/(2.0d0*hxc**6)
  case(16) ; s= 0.0d0
  case(31) ; s= 0.0d0
  case(52) ; s= 0.0d0
  case(80) ; s= 0.0d0
  case(30) ; s= 0.0d0
  case(51) ; s= 0.0d0
  case(79) ; s= 0.0d0
  case(50) ; s= 0.0d0
  case(78) ; s= 0.0d0
  case(77) ; s= 0.0d0
  case(1) ; s= 0.0d0
  case(6) ; s= 0.0d0
  case(14) ; s= 0.0d0
  case(27) ; s= 0.0d0
  case(49) ; s= xd1/hxc**5
  case(76) ; s= (xd1*zd1)/hxc**6
  case(5) ; s= 0.0d0
  case(15) ; s= 0.0d0
  case(29) ; s= 0.0d0
  case(48) ; s= 0.0d0
  case(75) ; s= (xd1*yd1)/hxc**6
  case(13) ; s= 0.0d0
  case(28) ; s= 0.0d0
  case(47) ; s= 0.0d0
  case(74) ; s= 0.0d0
  case(26) ; s= 0.0d0
  case(46) ; s= 0.0d0
  case(73) ; s= 0.0d0
  case(45) ; s= 0.0d0
  case(72) ; s= 0.0d0
  case(71) ; s= 0.0d0
  case(4) ; s= 0.0d0
  case(12) ; s= 0.0d0
  case(24) ; s= 0.0d0
  case(44) ; s= 0.0d0
  case(70) ; s= xd1**2/(2.0d0*hxc**6)
  case(11) ; s= 0.0d0
  case(25) ; s= 0.0d0
  case(43) ; s= 0.0d0
  case(69) ; s= 0.0d0
  case(23) ; s= 0.0d0
  case(42) ; s= 0.0d0
  case(68) ; s= 0.0d0
  case(41) ; s= 0.0d0
  case(67) ; s= 0.0d0
  case(66) ; s= 0.0d0
  case(10) ; s= 0.0d0
  case(22) ; s= 0.0d0
  case(39) ; s= 0.0d0
  case(65) ; s= 0.0d0
  case(21) ; s= 0.0d0
  case(40) ; s= 0.0d0
  case(64) ; s= 0.0d0
  case(38) ; s= 0.0d0
  case(63) ; s= 0.0d0
  case(62) ; s= 0.0d0
  case(20) ; s= 0.0d0
  case(37) ; s= 0.0d0
  case(61) ; s= 0.0d0
  case(36) ; s= 0.0d0
  case(60) ; s= 0.0d0
  case(59) ; s= 0.0d0
  case(35) ; s= 0.0d0
  case(58) ; s= 0.0d0
  case(57) ; s= 0.0d0
  case(56) ; s= 0.0d0
  
  
  
  end select

    tl3dz4 = s

  end function
 
 
 
  real function  tl3dz5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= hxc**(-5)
 case(83) ; s= zd1/hxc**6
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= yd1/hxc**6
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= xd1/hxc**6
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
   end select

   tl3dz5 = s

  end function
 
 
 
  real function tl3dz6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
  case(3) ; s= 0.0d0
  case(9) ; s= 0.0d0
  case(19) ; s= 0.0d0
  case(34) ; s= 0.0d0
  case(55) ; s= 0.0d0
  case(83) ; s= hxc**(-6)
  case(2) ; s= 0.0d0
  case(8) ; s= 0.0d0
  case(18) ; s= 0.0d0
  case(33) ; s= 0.0d0
  case(54) ; s= 0.0d0
  case(82) ; s= 0.0d0
  case(7) ; s= 0.0d0
  case(17) ; s= 0.0d0
  case(32) ; s= 0.0d0
  case(53) ; s= 0.0d0
  case(81) ; s= 0.0d0
  case(16) ; s= 0.0d0
  case(31) ; s= 0.0d0
  case(52) ; s= 0.0d0
  case(80) ; s= 0.0d0
  case(30) ; s= 0.0d0
  case(51) ; s= 0.0d0
  case(79) ; s= 0.0d0
  case(50) ; s= 0.0d0
  case(78) ; s= 0.0d0
  case(77) ; s= 0.0d0
  case(1) ; s= 0.0d0
  case(6) ; s= 0.0d0
  case(14) ; s= 0.0d0
  case(27) ; s= 0.0d0
  case(49) ; s= 0.0d0
  case(76) ; s= 0.0d0
  case(5) ; s= 0.0d0
  case(15) ; s= 0.0d0
  case(29) ; s= 0.0d0
  case(48) ; s= 0.0d0
  case(75) ; s= 0.0d0
  case(13) ; s= 0.0d0
  case(28) ; s= 0.0d0
  case(47) ; s= 0.0d0
  case(74) ; s= 0.0d0
  case(26) ; s= 0.0d0
  case(46) ; s= 0.0d0
  case(73) ; s= 0.0d0
  case(45) ; s= 0.0d0
  case(72) ; s= 0.0d0
  case(71) ; s= 0.0d0
  case(4) ; s= 0.0d0
  case(12) ; s= 0.0d0
  case(24) ; s= 0.0d0
  case(44) ; s= 0.0d0
  case(70) ; s= 0.0d0
  case(11) ; s= 0.0d0
  case(25) ; s= 0.0d0
  case(43) ; s= 0.0d0
  case(69) ; s= 0.0d0
  case(23) ; s= 0.0d0
  case(42) ; s= 0.0d0
  case(68) ; s= 0.0d0
  case(41) ; s= 0.0d0
  case(67) ; s= 0.0d0
  case(66) ; s= 0.0d0
  case(10) ; s= 0.0d0
  case(22) ; s= 0.0d0
  case(39) ; s= 0.0d0
  case(65) ; s= 0.0d0
  case(21) ; s= 0.0d0
  case(40) ; s= 0.0d0
  case(64) ; s= 0.0d0
  case(38) ; s= 0.0d0
  case(63) ; s= 0.0d0
  case(62) ; s= 0.0d0
  case(20) ; s= 0.0d0
  case(37) ; s= 0.0d0
  case(61) ; s= 0.0d0
  case(36) ; s= 0.0d0
  case(60) ; s= 0.0d0
  case(59) ; s= 0.0d0
  case(35) ; s= 0.0d0
  case(58) ; s= 0.0d0
  case(57) ; s= 0.0d0
  case(56) ; s= 0.0d0
  
  
  
  end select

    tl3dz6= s

  end function
 
 
 
  real function tl3dy(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 1/hxc
 case(8) ; s= zd1/hxc**2
 case(18) ; s= zd1**2/(2.0d0*hxc**3)
 case(33) ; s= zd1**3/(6.0d0*hxc**4)
 case(54) ; s= zd1**4/(24.0d0*hxc**5)
 case(82) ; s= zd1**5/(120.0d0*hxc**6)
 case(7) ; s= yd1/hxc**2
 case(17) ; s= (yd1*zd1)/hxc**3
 case(32) ; s= (yd1*zd1**2)/(2.0d0*hxc**4)
 case(53) ; s= (yd1*zd1**3)/(6.0d0*hxc**5)
 case(81) ; s= (yd1*zd1**4)/(24.0d0*hxc**6)
 case(16) ; s= yd1**2/(2.0d0*hxc**3)
 case(31) ; s= (yd1**2*zd1)/(2.0d0*hxc**4)
 case(52) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**5)
 case(80) ; s= (yd1**2*zd1**3)/(12.0d0*hxc**6)
 case(30) ; s= yd1**3/(6.0d0*hxc**4)
 case(51) ; s= (yd1**3*zd1)/(6.0d0*hxc**5)
 case(79) ; s= (yd1**3*zd1**2)/(12.0d0*hxc**6)
 case(50) ; s= yd1**4/(24.0d0*hxc**5)
 case(78) ; s= (yd1**4*zd1)/(24.0d0*hxc**6)
 case(77) ; s= yd1**5/(120.0d0*hxc**6)
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= xd1/hxc**2
 case(15) ; s= (xd1*zd1)/hxc**3
 case(29) ; s= (xd1*zd1**2)/(2.0d0*hxc**4)
 case(48) ; s= (xd1*zd1**3)/(6.0d0*hxc**5)
 case(75) ; s= (xd1*zd1**4)/(24.0d0*hxc**6)
 case(13) ; s= (xd1*yd1)/hxc**3
 case(28) ; s= (xd1*yd1*zd1)/hxc**4
 case(47) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**5)
 case(74) ; s= (xd1*yd1*zd1**3)/(6.0d0*hxc**6)
 case(26) ; s= (xd1*yd1**2)/(2.0d0*hxc**4)
 case(46) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**5)
 case(73) ; s= (xd1*yd1**2*zd1**2)/(4.0d0*hxc**6)
 case(45) ; s= (xd1*yd1**3)/(6.0d0*hxc**5)
 case(72) ; s= (xd1*yd1**3*zd1)/(6.0d0*hxc**6)
 case(71) ; s= (xd1*yd1**4)/(24.0d0*hxc**6)
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= xd1**2/(2.0d0*hxc**3)
 case(25) ; s= (xd1**2*zd1)/(2.0d0*hxc**4)
 case(43) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**5)
 case(69) ; s= (xd1**2*zd1**3)/(12.0d0*hxc**6)
 case(23) ; s= (xd1**2*yd1)/(2.0d0*hxc**4)
 case(42) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**5)
 case(68) ; s= (xd1**2*yd1*zd1**2)/(4.0d0*hxc**6)
 case(41) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**5)
 case(67) ; s= (xd1**2*yd1**2*zd1)/(4.0d0*hxc**6)
 case(66) ; s= (xd1**2*yd1**3)/(12.0d0*hxc**6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= xd1**3/(6.0d0*hxc**4)
 case(40) ; s= (xd1**3*zd1)/(6.0d0*hxc**5)
 case(64) ; s= (xd1**3*zd1**2)/(12.0d0*hxc**6)
 case(38) ; s= (xd1**3*yd1)/(6.0d0*hxc**5)
 case(63) ; s= (xd1**3*yd1*zd1)/(6.0d0*hxc**6)
 case(62) ; s= (xd1**3*yd1**2)/(12.0d0*hxc**6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= xd1**4/(24.0d0*hxc**5)
 case(60) ; s= (xd1**4*zd1)/(24.0d0*hxc**6)
 case(59) ; s= (xd1**4*yd1)/(24.0d0*hxc**6)
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= xd1**5/(120.0d0*hxc**6)
 case(56) ; s= 0.0d0
 
 
 
   end select

    tl3dy = s

  end function
 
 
 
  real function tl3dyz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
  case(3) ; s= 0.0d0
  case(9) ; s= 0.0d0
  case(19) ; s= 0.0d0
  case(34) ; s= 0.0d0
  case(55) ; s= 0.0d0
  case(83) ; s= 0.0d0
  case(2) ; s= 0.0d0
  case(8) ; s= hxc**(-2)
  case(18) ; s= zd1/hxc**3
  case(33) ; s= zd1**2/(2.0d0*hxc**4)
  case(54) ; s= zd1**3/(6.0d0*hxc**5)
  case(82) ; s= zd1**4/(24.0d0*hxc**6)
  case(7) ; s= 0.0d0
  case(17) ; s= yd1/hxc**3
  case(32) ; s= (yd1*zd1)/hxc**4
  case(53) ; s= (yd1*zd1**2)/(2.0d0*hxc**5)
  case(81) ; s= (yd1*zd1**3)/(6.0d0*hxc**6)
  case(16) ; s= 0.0d0
  case(31) ; s= yd1**2/(2.0d0*hxc**4)
  case(52) ; s= (yd1**2*zd1)/(2.0d0*hxc**5)
  case(80) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**6)
  case(30) ; s= 0.0d0
  case(51) ; s= yd1**3/(6.0d0*hxc**5)
  case(79) ; s= (yd1**3*zd1)/(6.0d0*hxc**6)
  case(50) ; s= 0.0d0
  case(78) ; s= yd1**4/(24.0d0*hxc**6)
  case(77) ; s= 0.0d0
  case(1) ; s= 0.0d0
  case(6) ; s= 0.0d0
  case(14) ; s= 0.0d0
  case(27) ; s= 0.0d0
  case(49) ; s= 0.0d0
  case(76) ; s= 0.0d0
  case(5) ; s= 0.0d0
  case(15) ; s= xd1/hxc**3
  case(29) ; s= (xd1*zd1)/hxc**4
  case(48) ; s= (xd1*zd1**2)/(2.0d0*hxc**5)
  case(75) ; s= (xd1*zd1**3)/(6.0d0*hxc**6)
  case(13) ; s= 0.0d0
  case(28) ; s= (xd1*yd1)/hxc**4
  case(47) ; s= (xd1*yd1*zd1)/hxc**5
  case(74) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**6)
  case(26) ; s= 0.0d0
  case(46) ; s= (xd1*yd1**2)/(2.0d0*hxc**5)
  case(73) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**6)
  case(45) ; s= 0.0d0
  case(72) ; s= (xd1*yd1**3)/(6.0d0*hxc**6)
  case(71) ; s= 0.0d0
  case(4) ; s= 0.0d0
  case(12) ; s= 0.0d0
  case(24) ; s= 0.0d0
  case(44) ; s= 0.0d0
  case(70) ; s= 0.0d0
  case(11) ; s= 0.0d0
  case(25) ; s= xd1**2/(2.0d0*hxc**4)
  case(43) ; s= (xd1**2*zd1)/(2.0d0*hxc**5)
  case(69) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**6)
  case(23) ; s= 0.0d0
  case(42) ; s= (xd1**2*yd1)/(2.0d0*hxc**5)
  case(68) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**6)
  case(41) ; s= 0.0d0
  case(67) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**6)
  case(66) ; s= 0.0d0
  case(10) ; s= 0.0d0
  case(22) ; s= 0.0d0
  case(39) ; s= 0.0d0
  case(65) ; s= 0.0d0
  case(21) ; s= 0.0d0
  case(40) ; s= xd1**3/(6.0d0*hxc**5)
  case(64) ; s= (xd1**3*zd1)/(6.0d0*hxc**6)
  case(38) ; s= 0.0d0
  case(63) ; s= (xd1**3*yd1)/(6.0d0*hxc**6)
  case(62) ; s= 0.0d0
  case(20) ; s= 0.0d0
  case(37) ; s= 0.0d0
  case(61) ; s= 0.0d0
  case(36) ; s= 0.0d0
  case(60) ; s= xd1**4/(24.0d0*hxc**6)
  case(59) ; s= 0.0d0
  case(35) ; s= 0.0d0
  case(58) ; s= 0.0d0
  case(57) ; s= 0.0d0
  case(56) ; s= 0.0d0
  
  
  
  
   end select

   tl3dyz = s

  end function
 
 
 
  real function tl3dyz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
  case(3) ; s= 0.0d0
  case(9) ; s= 0.0d0
  case(19) ; s= 0.0d0
  case(34) ; s= 0.0d0
  case(55) ; s= 0.0d0
  case(83) ; s= 0.0d0
  case(2) ; s= 0.0d0
  case(8) ; s= 0.0d0
  case(18) ; s= hxc**(-3)
  case(33) ; s= zd1/hxc**4
  case(54) ; s= zd1**2/(2.0d0*hxc**5)
  case(82) ; s= zd1**3/(6.0d0*hxc**6)
  case(7) ; s= 0.0d0
  case(17) ; s= 0.0d0
  case(32) ; s= yd1/hxc**4
  case(53) ; s= (yd1*zd1)/hxc**5
  case(81) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
  case(16) ; s= 0.0d0
  case(31) ; s= 0.0d0
  case(52) ; s= yd1**2/(2.0d0*hxc**5)
  case(80) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
  case(30) ; s= 0.0d0
  case(51) ; s= 0.0d0
  case(79) ; s= yd1**3/(6.0d0*hxc**6)
  case(50) ; s= 0.0d0
  case(78) ; s= 0.0d0
  case(77) ; s= 0.0d0
  case(1) ; s= 0.0d0
  case(6) ; s= 0.0d0
  case(14) ; s= 0.0d0
  case(27) ; s= 0.0d0
  case(49) ; s= 0.0d0
  case(76) ; s= 0.0d0
  case(5) ; s= 0.0d0
  case(15) ; s= 0.0d0
  case(29) ; s= xd1/hxc**4
  case(48) ; s= (xd1*zd1)/hxc**5
  case(75) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
  case(13) ; s= 0.0d0
  case(28) ; s= 0.0d0
  case(47) ; s= (xd1*yd1)/hxc**5
  case(74) ; s= (xd1*yd1*zd1)/hxc**6
  case(26) ; s= 0.0d0
  case(46) ; s= 0.0d0
  case(73) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
  case(45) ; s= 0.0d0
  case(72) ; s= 0.0d0
  case(71) ; s= 0.0d0
  case(4) ; s= 0.0d0
  case(12) ; s= 0.0d0
  case(24) ; s= 0.0d0
  case(44) ; s= 0.0d0
  case(70) ; s= 0.0d0
  case(11) ; s= 0.0d0
  case(25) ; s= 0.0d0
  case(43) ; s= xd1**2/(2.0d0*hxc**5)
  case(69) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
  case(23) ; s= 0.0d0
  case(42) ; s= 0.0d0
  case(68) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
  case(41) ; s= 0.0d0
  case(67) ; s= 0.0d0
  case(66) ; s= 0.0d0
  case(10) ; s= 0.0d0
  case(22) ; s= 0.0d0
  case(39) ; s= 0.0d0
  case(65) ; s= 0.0d0
  case(21) ; s= 0.0d0
  case(40) ; s= 0.0d0
  case(64) ; s= xd1**3/(6.0d0*hxc**6)
  case(38) ; s= 0.0d0
  case(63) ; s= 0.0d0
  case(62) ; s= 0.0d0
  case(20) ; s= 0.0d0
  case(37) ; s= 0.0d0
  case(61) ; s= 0.0d0
  case(36) ; s= 0.0d0
  case(60) ; s= 0.0d0
  case(59) ; s= 0.0d0
  case(35) ; s= 0.0d0
  case(58) ; s= 0.0d0
  case(57) ; s= 0.0d0
  case(56) ; s= 0.0d0
  
  
  
  end select

    tl3dyz2 = s

  end function
 
 
 
  real function tl3dyz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= hxc**(-4)
 case(54) ; s= zd1/hxc**5
 case(82) ; s= zd1**2/(2.0d0*hxc**6)
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= yd1/hxc**5
 case(81) ; s= (yd1*zd1)/hxc**6
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= yd1**2/(2.0d0*hxc**6)
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= xd1/hxc**5
 case(75) ; s= (xd1*zd1)/hxc**6
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= (xd1*yd1)/hxc**6
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= xd1**2/(2.0d0*hxc**6)
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

   tl3dyz3 = s

  end function
 
 
 
  real function tl3dyz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= hxc**(-5)
 case(82) ; s= zd1/hxc**6
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= yd1/hxc**6
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= xd1/hxc**6
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dyz4 = s

  end function
 
 
 
  real function tl3dyz5 (xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= hxc**(-6)
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dyz5 = s

  end function
 
 
 
  real function tl3dy2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= hxc**(-2)
 case(17) ; s= zd1/hxc**3
 case(32) ; s= zd1**2/(2.0d0*hxc**4)
 case(53) ; s= zd1**3/(6.0d0*hxc**5)
 case(81) ; s= zd1**4/(24.0d0*hxc**6)
 case(16) ; s= yd1/hxc**3
 case(31) ; s= (yd1*zd1)/hxc**4
 case(52) ; s= (yd1*zd1**2)/(2.0d0*hxc**5)
 case(80) ; s= (yd1*zd1**3)/(6.0d0*hxc**6)
 case(30) ; s= yd1**2/(2.0d0*hxc**4)
 case(51) ; s= (yd1**2*zd1)/(2.0d0*hxc**5)
 case(79) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**6)
 case(50) ; s= yd1**3/(6.0d0*hxc**5)
 case(78) ; s= (yd1**3*zd1)/(6.0d0*hxc**6)
 case(77) ; s= yd1**4/(24.0d0*hxc**6)
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= xd1/hxc**3
 case(28) ; s= (xd1*zd1)/hxc**4
 case(47) ; s= (xd1*zd1**2)/(2.0d0*hxc**5)
 case(74) ; s= (xd1*zd1**3)/(6.0d0*hxc**6)
 case(26) ; s= (xd1*yd1)/hxc**4
 case(46) ; s= (xd1*yd1*zd1)/hxc**5
 case(73) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**6)
 case(45) ; s= (xd1*yd1**2)/(2.0d0*hxc**5)
 case(72) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**6)
 case(71) ; s= (xd1*yd1**3)/(6.0d0*hxc**6)
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= xd1**2/(2.0d0*hxc**4)
 case(42) ; s= (xd1**2*zd1)/(2.0d0*hxc**5)
 case(68) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**6)
 case(41) ; s= (xd1**2*yd1)/(2.0d0*hxc**5)
 case(67) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**6)
 case(66) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= xd1**3/(6.0d0*hxc**5)
 case(63) ; s= (xd1**3*zd1)/(6.0d0*hxc**6)
 case(62) ; s= (xd1**3*yd1)/(6.0d0*hxc**6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= xd1**4/(24.0d0*hxc**6)
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
 
  end select

    tl3dy2 = s

  end function
 
 
 
  real function tl3dy2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= hxc**(-3)
 case(32) ; s= zd1/hxc**4
 case(53) ; s= zd1**2/(2.0d0*hxc**5)
 case(81) ; s= zd1**3/(6.0d0*hxc**6)
 case(16) ; s= 0.0d0
 case(31) ; s= yd1/hxc**4
 case(52) ; s= (yd1*zd1)/hxc**5
 case(80) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
 case(30) ; s= 0.0d0
 case(51) ; s= yd1**2/(2.0d0*hxc**5)
 case(79) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
 case(50) ; s= 0.0d0
 case(78) ; s= yd1**3/(6.0d0*hxc**6)
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= xd1/hxc**4
 case(47) ; s= (xd1*zd1)/hxc**5
 case(74) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
 case(26) ; s= 0.0d0
 case(46) ; s= (xd1*yd1)/hxc**5
 case(73) ; s= (xd1*yd1*zd1)/hxc**6
 case(45) ; s= 0.0d0
 case(72) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= xd1**2/(2.0d0*hxc**5)
 case(68) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
 case(41) ; s= 0.0d0
 case(67) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= xd1**3/(6.0d0*hxc**6)
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dy2z = s

  end function
 
 
 
  real function tl3dy2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= hxc**(-4)
 case(53) ; s= zd1/hxc**5
 case(81) ; s= zd1**2/(2.0d0*hxc**6)
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= yd1/hxc**5
 case(80) ; s= (yd1*zd1)/hxc**6
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= yd1**2/(2.0d0*hxc**6)
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= xd1/hxc**5
 case(74) ; s= (xd1*zd1)/hxc**6
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= (xd1*yd1)/hxc**6
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= xd1**2/(2.0d0*hxc**6)
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dy2z2 = s

  end function
 
 
 
  real function tl3dy2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= hxc**(-5)
 case(81) ; s= zd1/hxc**6
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= yd1/hxc**6
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= xd1/hxc**6
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dy2z3 = s

  end function
 
 
 
  real function tl3dy2z4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= hxc**(-6)
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dy2z4 = s

  end function
 
 
 
  real function tl3dy3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= hxc**(-3)
 case(31) ; s= zd1/hxc**4
 case(52) ; s= zd1**2/(2.0d0*hxc**5)
 case(80) ; s= zd1**3/(6.0d0*hxc**6)
 case(30) ; s= yd1/hxc**4
 case(51) ; s= (yd1*zd1)/hxc**5
 case(79) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
 case(50) ; s= yd1**2/(2.0d0*hxc**5)
 case(78) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
 case(77) ; s= yd1**3/(6.0d0*hxc**6)
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= xd1/hxc**4
 case(46) ; s= (xd1*zd1)/hxc**5
 case(73) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
 case(45) ; s= (xd1*yd1)/hxc**5
 case(72) ; s= (xd1*yd1*zd1)/hxc**6
 case(71) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= xd1**2/(2.0d0*hxc**5)
 case(67) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
 case(66) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= xd1**3/(6.0d0*hxc**6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dy3 = s

  end function
 
 
 
  real function tl3dy3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= hxc**(-4)
 case(52) ; s= zd1/hxc**5
 case(80) ; s= zd1**2/(2.0d0*hxc**6)
 case(30) ; s= 0.0d0
 case(51) ; s= yd1/hxc**5
 case(79) ; s= (yd1*zd1)/hxc**6
 case(50) ; s= 0.0d0
 case(78) ; s= yd1**2/(2.0d0*hxc**6)
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= xd1/hxc**5
 case(73) ; s= (xd1*zd1)/hxc**6
 case(45) ; s= 0.0d0
 case(72) ; s= (xd1*yd1)/hxc**6
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= xd1**2/(2.0d0*hxc**6)
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dy3z = s

  end function
 
 
 
  real function tl3dy3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= hxc**(-5)
 case(80) ; s= zd1/hxc**6
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= yd1/hxc**6
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= xd1/hxc**6
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dy3z2 = s

  end function
 
 
 
  real function tl3dy3z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= hxc**(-6)
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dy3z3 = s

  end function
 
 
 
  real function tl3dy4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= hxc**(-4)
 case(51) ; s= zd1/hxc**5
 case(79) ; s= zd1**2/(2.0d0*hxc**6)
 case(50) ; s= yd1/hxc**5
 case(78) ; s= (yd1*zd1)/hxc**6
 case(77) ; s= yd1**2/(2.0d0*hxc**6)
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= xd1/hxc**5
 case(72) ; s= (xd1*zd1)/hxc**6
 case(71) ; s= (xd1*yd1)/hxc**6
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= xd1**2/(2.0d0*hxc**6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dy4 = s

  end function
 
 
 
  real function tl3dy4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= hxc**(-5)
 case(79) ; s= zd1/hxc**6
 case(50) ; s= 0.0d0
 case(78) ; s= yd1/hxc**6
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= xd1/hxc**6
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dy4z = s

  end function
 
 
 
  real function tl3dy4z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= hxc**(-6)
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dy4z2 = s

  end function
 
 
 
  real function tl3dy5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= hxc**(-5)
 case(78) ; s= zd1/hxc**6
 case(77) ; s= yd1/hxc**6
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= xd1/hxc**6
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dy5 = s

  end function
 
 
 
  real function tl3dy5z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= hxc**(-6)
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dy5z = s

  end function
 
 
 
  real function tl3dy6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= hxc**(-6)
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dy6 = s

  end function
 
 
 
  real function tl3dx(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 1/hxc
 case(6) ; s= zd1/hxc**2
 case(14) ; s= zd1**2/(2.0d0*hxc**3)
 case(27) ; s= zd1**3/(6.0d0*hxc**4)
 case(49) ; s= zd1**4/(24.0d0*hxc**5)
 case(76) ; s= zd1**5/(120.0d0*hxc**6)
 case(5) ; s= yd1/hxc**2
 case(15) ; s= (yd1*zd1)/hxc**3
 case(29) ; s= (yd1*zd1**2)/(2.0d0*hxc**4)
 case(48) ; s= (yd1*zd1**3)/(6.0d0*hxc**5)
 case(75) ; s= (yd1*zd1**4)/(24.0d0*hxc**6)
 case(13) ; s= yd1**2/(2.0d0*hxc**3)
 case(28) ; s= (yd1**2*zd1)/(2.0d0*hxc**4)
 case(47) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**5)
 case(74) ; s= (yd1**2*zd1**3)/(12.0d0*hxc**6)
 case(26) ; s= yd1**3/(6.0d0*hxc**4)
 case(46) ; s= (yd1**3*zd1)/(6.0d0*hxc**5)
 case(73) ; s= (yd1**3*zd1**2)/(12.0d0*hxc**6)
 case(45) ; s= yd1**4/(24.0d0*hxc**5)
 case(72) ; s= (yd1**4*zd1)/(24.0d0*hxc**6)
 case(71) ; s= yd1**5/(120.0d0*hxc**6)
 case(4) ; s= xd1/hxc**2
 case(12) ; s= (xd1*zd1)/hxc**3
 case(24) ; s= (xd1*zd1**2)/(2.0d0*hxc**4)
 case(44) ; s= (xd1*zd1**3)/(6.0d0*hxc**5)
 case(70) ; s= (xd1*zd1**4)/(24.0d0*hxc**6)
 case(11) ; s= (xd1*yd1)/hxc**3
 case(25) ; s= (xd1*yd1*zd1)/hxc**4
 case(43) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**5)
 case(69) ; s= (xd1*yd1*zd1**3)/(6.0d0*hxc**6)
 case(23) ; s= (xd1*yd1**2)/(2.0d0*hxc**4)
 case(42) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**5)
 case(68) ; s= (xd1*yd1**2*zd1**2)/(4.0d0*hxc**6)
 case(41) ; s= (xd1*yd1**3)/(6.0d0*hxc**5)
 case(67) ; s= (xd1*yd1**3*zd1)/(6.0d0*hxc**6)
 case(66) ; s= (xd1*yd1**4)/(24.0d0*hxc**6)
 case(10) ; s= xd1**2/(2.0d0*hxc**3)
 case(22) ; s= (xd1**2*zd1)/(2.0d0*hxc**4)
 case(39) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**5)
 case(65) ; s= (xd1**2*zd1**3)/(12.0d0*hxc**6)
 case(21) ; s= (xd1**2*yd1)/(2.0d0*hxc**4)
 case(40) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**5)
 case(64) ; s= (xd1**2*yd1*zd1**2)/(4.0d0*hxc**6)
 case(38) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**5)
 case(63) ; s= (xd1**2*yd1**2*zd1)/(4.0d0*hxc**6)
 case(62) ; s= (xd1**2*yd1**3)/(12.0d0*hxc**6)
 case(20) ; s= xd1**3/(6.0d0*hxc**4)
 case(37) ; s= (xd1**3*zd1)/(6.0d0*hxc**5)
 case(61) ; s= (xd1**3*zd1**2)/(12.0d0*hxc**6)
 case(36) ; s= (xd1**3*yd1)/(6.0d0*hxc**5)
 case(60) ; s= (xd1**3*yd1*zd1)/(6.0d0*hxc**6)
 case(59) ; s= (xd1**3*yd1**2)/(12.0d0*hxc**6)
 case(35) ; s= xd1**4/(24.0d0*hxc**5)
 case(58) ; s= (xd1**4*zd1)/(24.0d0*hxc**6)
 case(57) ; s= (xd1**4*yd1)/(24.0d0*hxc**6)
 case(56) ; s= xd1**5/(120.0d0*hxc**6)
 
  end select

   tl3dx = s

  end function
 
 
 
  real function tl3dxz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= hxc**(-2)
 case(14) ; s= zd1/hxc**3
 case(27) ; s= zd1**2/(2.0d0*hxc**4)
 case(49) ; s= zd1**3/(6.0d0*hxc**5)
 case(76) ; s= zd1**4/(24.0d0*hxc**6)
 case(5) ; s= 0.0d0
 case(15) ; s= yd1/hxc**3
 case(29) ; s= (yd1*zd1)/hxc**4
 case(48) ; s= (yd1*zd1**2)/(2.0d0*hxc**5)
 case(75) ; s= (yd1*zd1**3)/(6.0d0*hxc**6)
 case(13) ; s= 0.0d0
 case(28) ; s= yd1**2/(2.0d0*hxc**4)
 case(47) ; s= (yd1**2*zd1)/(2.0d0*hxc**5)
 case(74) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**6)
 case(26) ; s= 0.0d0
 case(46) ; s= yd1**3/(6.0d0*hxc**5)
 case(73) ; s= (yd1**3*zd1)/(6.0d0*hxc**6)
 case(45) ; s= 0.0d0
 case(72) ; s= yd1**4/(24.0d0*hxc**6)
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= xd1/hxc**3
 case(24) ; s= (xd1*zd1)/hxc**4
 case(44) ; s= (xd1*zd1**2)/(2.0d0*hxc**5)
 case(70) ; s= (xd1*zd1**3)/(6.0d0*hxc**6)
 case(11) ; s= 0.0d0
 case(25) ; s= (xd1*yd1)/hxc**4
 case(43) ; s= (xd1*yd1*zd1)/hxc**5
 case(69) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**6)
 case(23) ; s= 0.0d0
 case(42) ; s= (xd1*yd1**2)/(2.0d0*hxc**5)
 case(68) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**6)
 case(41) ; s= 0.0d0
 case(67) ; s= (xd1*yd1**3)/(6.0d0*hxc**6)
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= xd1**2/(2.0d0*hxc**4)
 case(39) ; s= (xd1**2*zd1)/(2.0d0*hxc**5)
 case(65) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**6)
 case(21) ; s= 0.0d0
 case(40) ; s= (xd1**2*yd1)/(2.0d0*hxc**5)
 case(64) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**6)
 case(38) ; s= 0.0d0
 case(63) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**6)
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= xd1**3/(6.0d0*hxc**5)
 case(61) ; s= (xd1**3*zd1)/(6.0d0*hxc**6)
 case(36) ; s= 0.0d0
 case(60) ; s= (xd1**3*yd1)/(6.0d0*hxc**6)
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= xd1**4/(24.0d0*hxc**6)
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxz = s

  end function
 
 
 
  real function tl3dxz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= hxc**(-3)
 case(27) ; s= zd1/hxc**4
 case(49) ; s= zd1**2/(2.0d0*hxc**5)
 case(76) ; s= zd1**3/(6.0d0*hxc**6)
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= yd1/hxc**4
 case(48) ; s= (yd1*zd1)/hxc**5
 case(75) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= yd1**2/(2.0d0*hxc**5)
 case(74) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= yd1**3/(6.0d0*hxc**6)
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= xd1/hxc**4
 case(44) ; s= (xd1*zd1)/hxc**5
 case(70) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= (xd1*yd1)/hxc**5
 case(69) ; s= (xd1*yd1*zd1)/hxc**6
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= xd1**2/(2.0d0*hxc**5)
 case(65) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= xd1**3/(6.0d0*hxc**6)
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxz2 = s

  end function
 
 
 
  real function tl3dxz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= hxc**(-4)
 case(49) ; s= zd1/hxc**5
 case(76) ; s= zd1**2/(2.0d0*hxc**6)
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= yd1/hxc**5
 case(75) ; s= (yd1*zd1)/hxc**6
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= yd1**2/(2.0d0*hxc**6)
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= xd1/hxc**5
 case(70) ; s= (xd1*zd1)/hxc**6
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= (xd1*yd1)/hxc**6
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= xd1**2/(2.0d0*hxc**6)
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxz3 = s

  end function
 
 
 
  real function tl3dxz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= hxc**(-5)
 case(76) ; s= zd1/hxc**6
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= yd1/hxc**6
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= xd1/hxc**6
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxz4 = s

  end function
 
 
 
  real function tl3dxz5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= hxc**(-6)
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dxz5 = s

  end function
 
 
 
  real function tl3dxy(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= hxc**(-2)
 case(15) ; s= zd1/hxc**3
 case(29) ; s= zd1**2/(2.0d0*hxc**4)
 case(48) ; s= zd1**3/(6.0d0*hxc**5)
 case(75) ; s= zd1**4/(24.0d0*hxc**6)
 case(13) ; s= yd1/hxc**3
 case(28) ; s= (yd1*zd1)/hxc**4
 case(47) ; s= (yd1*zd1**2)/(2.0d0*hxc**5)
 case(74) ; s= (yd1*zd1**3)/(6.0d0*hxc**6)
 case(26) ; s= yd1**2/(2.0d0*hxc**4)
 case(46) ; s= (yd1**2*zd1)/(2.0d0*hxc**5)
 case(73) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**6)
 case(45) ; s= yd1**3/(6.0d0*hxc**5)
 case(72) ; s= (yd1**3*zd1)/(6.0d0*hxc**6)
 case(71) ; s= yd1**4/(24.0d0*hxc**6)
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= xd1/hxc**3
 case(25) ; s= (xd1*zd1)/hxc**4
 case(43) ; s= (xd1*zd1**2)/(2.0d0*hxc**5)
 case(69) ; s= (xd1*zd1**3)/(6.0d0*hxc**6)
 case(23) ; s= (xd1*yd1)/hxc**4
 case(42) ; s= (xd1*yd1*zd1)/hxc**5
 case(68) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**6)
 case(41) ; s= (xd1*yd1**2)/(2.0d0*hxc**5)
 case(67) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**6)
 case(66) ; s= (xd1*yd1**3)/(6.0d0*hxc**6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= xd1**2/(2.0d0*hxc**4)
 case(40) ; s= (xd1**2*zd1)/(2.0d0*hxc**5)
 case(64) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**6)
 case(38) ; s= (xd1**2*yd1)/(2.0d0*hxc**5)
 case(63) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**6)
 case(62) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= xd1**3/(6.0d0*hxc**5)
 case(60) ; s= (xd1**3*zd1)/(6.0d0*hxc**6)
 case(59) ; s= (xd1**3*yd1)/(6.0d0*hxc**6)
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= xd1**4/(24.0d0*hxc**6)
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dxy = s

  end function
 
 
 
  real function tl3dxyz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= hxc**(-3)
 case(29) ; s= zd1/hxc**4
 case(48) ; s= zd1**2/(2.0d0*hxc**5)
 case(75) ; s= zd1**3/(6.0d0*hxc**6)
 case(13) ; s= 0.0d0
 case(28) ; s= yd1/hxc**4
 case(47) ; s= (yd1*zd1)/hxc**5
 case(74) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
 case(26) ; s= 0.0d0
 case(46) ; s= yd1**2/(2.0d0*hxc**5)
 case(73) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
 case(45) ; s= 0.0d0
 case(72) ; s= yd1**3/(6.0d0*hxc**6)
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= xd1/hxc**4
 case(43) ; s= (xd1*zd1)/hxc**5
 case(69) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
 case(23) ; s= 0.0d0
 case(42) ; s= (xd1*yd1)/hxc**5
 case(68) ; s= (xd1*yd1*zd1)/hxc**6
 case(41) ; s= 0.0d0
 case(67) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= xd1**2/(2.0d0*hxc**5)
 case(64) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
 case(38) ; s= 0.0d0
 case(63) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= xd1**3/(6.0d0*hxc**6)
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dxyz = s

  end function
 
 
 
  real function tl3dxyz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= hxc**(-4)
 case(48) ; s= zd1/hxc**5
 case(75) ; s= zd1**2/(2.0d0*hxc**6)
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= yd1/hxc**5
 case(74) ; s= (yd1*zd1)/hxc**6
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= yd1**2/(2.0d0*hxc**6)
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= xd1/hxc**5
 case(69) ; s= (xd1*zd1)/hxc**6
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= (xd1*yd1)/hxc**6
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= xd1**2/(2.0d0*hxc**6)
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxyz2 = s

  end function
 
 
 
  real function tl3dxyz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= hxc**(-5)
 case(75) ; s= zd1/hxc**6
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= yd1/hxc**6
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= xd1/hxc**6
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxyz3 = s

  end function
 
 
 
  real function tl3dxyz4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= hxc**(-6)
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxyz4 = s

  end function
 
 
 
  real function tl3dxy2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= hxc**(-3)
 case(28) ; s= zd1/hxc**4
 case(47) ; s= zd1**2/(2.0d0*hxc**5)
 case(74) ; s= zd1**3/(6.0d0*hxc**6)
 case(26) ; s= yd1/hxc**4
 case(46) ; s= (yd1*zd1)/hxc**5
 case(73) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
 case(45) ; s= yd1**2/(2.0d0*hxc**5)
 case(72) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
 case(71) ; s= yd1**3/(6.0d0*hxc**6)
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= xd1/hxc**4
 case(42) ; s= (xd1*zd1)/hxc**5
 case(68) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
 case(41) ; s= (xd1*yd1)/hxc**5
 case(67) ; s= (xd1*yd1*zd1)/hxc**6
 case(66) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= xd1**2/(2.0d0*hxc**5)
 case(63) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
 case(62) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= xd1**3/(6.0d0*hxc**6)
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxy2 = s

  end function
 
 
 
  real function tl3dxy2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= hxc**(-4)
 case(47) ; s= zd1/hxc**5
 case(74) ; s= zd1**2/(2.0d0*hxc**6)
 case(26) ; s= 0.0d0
 case(46) ; s= yd1/hxc**5
 case(73) ; s= (yd1*zd1)/hxc**6
 case(45) ; s= 0.0d0
 case(72) ; s= yd1**2/(2.0d0*hxc**6)
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= xd1/hxc**5
 case(68) ; s= (xd1*zd1)/hxc**6
 case(41) ; s= 0.0d0
 case(67) ; s= (xd1*yd1)/hxc**6
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= xd1**2/(2.0d0*hxc**6)
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dxy2z = s

  end function
 
 
 
  real function tl3dxy2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= hxc**(-5)
 case(74) ; s= zd1/hxc**6
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= yd1/hxc**6
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= xd1/hxc**6
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dxy2z2 = s

  end function
 
 
 
  real function tl3dxy2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= hxc**(-6)
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxy2z3 = s

  end function
 
 
 
  real function tl3dxy3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= hxc**(-4)
 case(46) ; s= zd1/hxc**5
 case(73) ; s= zd1**2/(2.0d0*hxc**6)
 case(45) ; s= yd1/hxc**5
 case(72) ; s= (yd1*zd1)/hxc**6
 case(71) ; s= yd1**2/(2.0d0*hxc**6)
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= xd1/hxc**5
 case(67) ; s= (xd1*zd1)/hxc**6
 case(66) ; s= (xd1*yd1)/hxc**6
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= xd1**2/(2.0d0*hxc**6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dxy3 = s

  end function
 
 
 
  real function tl3dxy3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= hxc**(-5)
 case(73) ; s= zd1/hxc**6
 case(45) ; s= 0.0d0
 case(72) ; s= yd1/hxc**6
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= xd1/hxc**6
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

   tl3dxy3z = s

  end function
 
 
 
  real function tl3dxy3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= hxc**(-6)
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dxy3z2 = s

  end function
 
 
 
  real function tl3dxy4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= hxc**(-5)
 case(72) ; s= zd1/hxc**6
 case(71) ; s= yd1/hxc**6
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= xd1/hxc**6
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dxy4 = s

  end function
 
 
 
  real function tl3dxy4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= hxc**(-6)
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxy4z = s

  end function
 
 
 
  real function tl3dxy5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= hxc**(-6)
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dxy5 = s

  end function
 
 
 
  real function tl3dx2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= hxc**(-2)
 case(12) ; s= zd1/hxc**3
 case(24) ; s= zd1**2/(2.0d0*hxc**4)
 case(44) ; s= zd1**3/(6.0d0*hxc**5)
 case(70) ; s= zd1**4/(24.0d0*hxc**6)
 case(11) ; s= yd1/hxc**3
 case(25) ; s= (yd1*zd1)/hxc**4
 case(43) ; s= (yd1*zd1**2)/(2.0d0*hxc**5)
 case(69) ; s= (yd1*zd1**3)/(6.0d0*hxc**6)
 case(23) ; s= yd1**2/(2.0d0*hxc**4)
 case(42) ; s= (yd1**2*zd1)/(2.0d0*hxc**5)
 case(68) ; s= (yd1**2*zd1**2)/(4.0d0*hxc**6)
 case(41) ; s= yd1**3/(6.0d0*hxc**5)
 case(67) ; s= (yd1**3*zd1)/(6.0d0*hxc**6)
 case(66) ; s= yd1**4/(24.0d0*hxc**6)
 case(10) ; s= xd1/hxc**3
 case(22) ; s= (xd1*zd1)/hxc**4
 case(39) ; s= (xd1*zd1**2)/(2.0d0*hxc**5)
 case(65) ; s= (xd1*zd1**3)/(6.0d0*hxc**6)
 case(21) ; s= (xd1*yd1)/hxc**4
 case(40) ; s= (xd1*yd1*zd1)/hxc**5
 case(64) ; s= (xd1*yd1*zd1**2)/(2.0d0*hxc**6)
 case(38) ; s= (xd1*yd1**2)/(2.0d0*hxc**5)
 case(63) ; s= (xd1*yd1**2*zd1)/(2.0d0*hxc**6)
 case(62) ; s= (xd1*yd1**3)/(6.0d0*hxc**6)
 case(20) ; s= xd1**2/(2.0d0*hxc**4)
 case(37) ; s= (xd1**2*zd1)/(2.0d0*hxc**5)
 case(61) ; s= (xd1**2*zd1**2)/(4.0d0*hxc**6)
 case(36) ; s= (xd1**2*yd1)/(2.0d0*hxc**5)
 case(60) ; s= (xd1**2*yd1*zd1)/(2.0d0*hxc**6)
 case(59) ; s= (xd1**2*yd1**2)/(4.0d0*hxc**6)
 case(35) ; s= xd1**3/(6.0d0*hxc**5)
 case(58) ; s= (xd1**3*zd1)/(6.0d0*hxc**6)
 case(57) ; s= (xd1**3*yd1)/(6.0d0*hxc**6)
 case(56) ; s= xd1**4/(24.0d0*hxc**6)
 
 
  end select

    tl3dx2 = s

  end function
 
 
 
  real function tl3dx2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= hxc**(-3)
 case(24) ; s= zd1/hxc**4
 case(44) ; s= zd1**2/(2.0d0*hxc**5)
 case(70) ; s= zd1**3/(6.0d0*hxc**6)
 case(11) ; s= 0.0d0
 case(25) ; s= yd1/hxc**4
 case(43) ; s= (yd1*zd1)/hxc**5
 case(69) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
 case(23) ; s= 0.0d0
 case(42) ; s= yd1**2/(2.0d0*hxc**5)
 case(68) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
 case(41) ; s= 0.0d0
 case(67) ; s= yd1**3/(6.0d0*hxc**6)
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= xd1/hxc**4
 case(39) ; s= (xd1*zd1)/hxc**5
 case(65) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
 case(21) ; s= 0.0d0
 case(40) ; s= (xd1*yd1)/hxc**5
 case(64) ; s= (xd1*yd1*zd1)/hxc**6
 case(38) ; s= 0.0d0
 case(63) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= xd1**2/(2.0d0*hxc**5)
 case(61) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
 case(36) ; s= 0.0d0
 case(60) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= xd1**3/(6.0d0*hxc**6)
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dx2z = s

  end function
 
 
 
  real function tl3dx2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= hxc**(-4)
 case(44) ; s= zd1/hxc**5
 case(70) ; s= zd1**2/(2.0d0*hxc**6)
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= yd1/hxc**5
 case(69) ; s= (yd1*zd1)/hxc**6
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= yd1**2/(2.0d0*hxc**6)
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= xd1/hxc**5
 case(65) ; s= (xd1*zd1)/hxc**6
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= (xd1*yd1)/hxc**6
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= xd1**2/(2.0d0*hxc**6)
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
  end select

    tl3dx2z2 = s

  end function
 
 
 
  real function tl3dx2z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= hxc**(-5)
 case(70) ; s= zd1/hxc**6
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= yd1/hxc**6
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= xd1/hxc**6
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx2z3 = s

  end function
 
 
 
  real function tl3dx2z4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= hxc**(-6)
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dx2z4 = s

  end function
 
 
 
  real function tl3dx2y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= hxc**(-3)
 case(25) ; s= zd1/hxc**4
 case(43) ; s= zd1**2/(2.0d0*hxc**5)
 case(69) ; s= zd1**3/(6.0d0*hxc**6)
 case(23) ; s= yd1/hxc**4
 case(42) ; s= (yd1*zd1)/hxc**5
 case(68) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
 case(41) ; s= yd1**2/(2.0d0*hxc**5)
 case(67) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
 case(66) ; s= yd1**3/(6.0d0*hxc**6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= xd1/hxc**4
 case(40) ; s= (xd1*zd1)/hxc**5
 case(64) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
 case(38) ; s= (xd1*yd1)/hxc**5
 case(63) ; s= (xd1*yd1*zd1)/hxc**6
 case(62) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= xd1**2/(2.0d0*hxc**5)
 case(60) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
 case(59) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= xd1**3/(6.0d0*hxc**6)
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx2y = s

  end function
 
 
 
  real function tl3dx2yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= hxc**(-4)
 case(43) ; s= zd1/hxc**5
 case(69) ; s= zd1**2/(2.0d0*hxc**6)
 case(23) ; s= 0.0d0
 case(42) ; s= yd1/hxc**5
 case(68) ; s= (yd1*zd1)/hxc**6
 case(41) ; s= 0.0d0
 case(67) ; s= yd1**2/(2.0d0*hxc**6)
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= xd1/hxc**5
 case(64) ; s= (xd1*zd1)/hxc**6
 case(38) ; s= 0.0d0
 case(63) ; s= (xd1*yd1)/hxc**6
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= xd1**2/(2.0d0*hxc**6)
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dx2yz = s

  end function
 
 
 
  real function tl3dx2yz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= hxc**(-5)
 case(69) ; s= zd1/hxc**6
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= yd1/hxc**6
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= xd1/hxc**6
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dx2yz2 = s

  end function
 
 
 
  real function tl3dx2yz3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= hxc**(-6)
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx2yz3 = s

  end function
 
 
 
  real function tl3dx2y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= hxc**(-4)
 case(42) ; s= zd1/hxc**5
 case(68) ; s= zd1**2/(2.0d0*hxc**6)
 case(41) ; s= yd1/hxc**5
 case(67) ; s= (yd1*zd1)/hxc**6
 case(66) ; s= yd1**2/(2.0d0*hxc**6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= xd1/hxc**5
 case(63) ; s= (xd1*zd1)/hxc**6
 case(62) ; s= (xd1*yd1)/hxc**6
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= xd1**2/(2.0d0*hxc**6)
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dx2y2 = s

  end function
 
 
 
  real function tl3dx2y2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= hxc**(-5)
 case(68) ; s= zd1/hxc**6
 case(41) ; s= 0.0d0
 case(67) ; s= yd1/hxc**6
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= xd1/hxc**6
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx2y2z = s

  end function
 
 
 
  real function tl3dx2y2z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= hxc**(-6)
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx2y2z2 = s

  end function
 
 
 
  real function tl3dx2y3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= hxc**(-5)
 case(67) ; s= zd1/hxc**6
 case(66) ; s= yd1/hxc**6
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= xd1/hxc**6
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx2y3 = s

  end function
 
 
 
  real function tl3dx2y3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= hxc**(-6)
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx2y3z = s

  end function
 
 
 
  real function tl3dx2y4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= hxc**(-6)
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx2y4 = s

  end function
 
 
 
  real function tl3dx3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= hxc**(-3)
 case(22) ; s= zd1/hxc**4
 case(39) ; s= zd1**2/(2.0d0*hxc**5)
 case(65) ; s= zd1**3/(6.0d0*hxc**6)
 case(21) ; s= yd1/hxc**4
 case(40) ; s= (yd1*zd1)/hxc**5
 case(64) ; s= (yd1*zd1**2)/(2.0d0*hxc**6)
 case(38) ; s= yd1**2/(2.0d0*hxc**5)
 case(63) ; s= (yd1**2*zd1)/(2.0d0*hxc**6)
 case(62) ; s= yd1**3/(6.0d0*hxc**6)
 case(20) ; s= xd1/hxc**4
 case(37) ; s= (xd1*zd1)/hxc**5
 case(61) ; s= (xd1*zd1**2)/(2.0d0*hxc**6)
 case(36) ; s= (xd1*yd1)/hxc**5
 case(60) ; s= (xd1*yd1*zd1)/hxc**6
 case(59) ; s= (xd1*yd1**2)/(2.0d0*hxc**6)
 case(35) ; s= xd1**2/(2.0d0*hxc**5)
 case(58) ; s= (xd1**2*zd1)/(2.0d0*hxc**6)
 case(57) ; s= (xd1**2*yd1)/(2.0d0*hxc**6)
 case(56) ; s= xd1**3/(6.0d0*hxc**6)
 
 
 
  end select

    tl3dx3 = s

  end function
 
 
 
  real function tl3dx3z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= hxc**(-4)
 case(39) ; s= zd1/hxc**5
 case(65) ; s= zd1**2/(2.0d0*hxc**6)
 case(21) ; s= 0.0d0
 case(40) ; s= yd1/hxc**5
 case(64) ; s= (yd1*zd1)/hxc**6
 case(38) ; s= 0.0d0
 case(63) ; s= yd1**2/(2.0d0*hxc**6)
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= xd1/hxc**5
 case(61) ; s= (xd1*zd1)/hxc**6
 case(36) ; s= 0.0d0
 case(60) ; s= (xd1*yd1)/hxc**6
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= xd1**2/(2.0d0*hxc**6)
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx3z = s

  end function
 
 
 
  real function tl3dx3z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= hxc**(-5)
 case(65) ; s= zd1/hxc**6
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= yd1/hxc**6
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= xd1/hxc**6
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx3z2 = s

  end function
 
 
 
  real function tl3dx3z3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= hxc**(-6)
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx3z3 = s

  end function
 
 
 
  real function tl3dx3y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= hxc**(-4)
 case(40) ; s= zd1/hxc**5
 case(64) ; s= zd1**2/(2.0d0*hxc**6)
 case(38) ; s= yd1/hxc**5
 case(63) ; s= (yd1*zd1)/hxc**6
 case(62) ; s= yd1**2/(2.0d0*hxc**6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= xd1/hxc**5
 case(60) ; s= (xd1*zd1)/hxc**6
 case(59) ; s= (xd1*yd1)/hxc**6
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= xd1**2/(2.0d0*hxc**6)
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx3y = s

  end function
 
 
 
  real function tl3dx3yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= hxc**(-5)
 case(64) ; s= zd1/hxc**6
 case(38) ; s= 0.0d0
 case(63) ; s= yd1/hxc**6
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= xd1/hxc**6
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dx3yz = s

  end function
 
 
 
  real function tl3dx3yz2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= hxc**(-6)
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dx3yz2 = s

  end function
 
 
 
  real function tl3dx3y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= hxc**(-5)
 case(63) ; s= zd1/hxc**6
 case(62) ; s= yd1/hxc**6
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= xd1/hxc**6
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx3y2 = s

  end function
 
 
 
  real function tl3dx3y2z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= hxc**(-6)
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx3y2z = s

  end function
 
 
 
  real function tl3dx3y3(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= hxc**(-6)
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx3y3 = s

  end function
 
 
 
  real function tl3dx4(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= hxc**(-4)
 case(37) ; s= zd1/hxc**5
 case(61) ; s= zd1**2/(2.0d0*hxc**6)
 case(36) ; s= yd1/hxc**5
 case(60) ; s= (yd1*zd1)/hxc**6
 case(59) ; s= yd1**2/(2.0d0*hxc**6)
 case(35) ; s= xd1/hxc**5
 case(58) ; s= (xd1*zd1)/hxc**6
 case(57) ; s= (xd1*yd1)/hxc**6
 case(56) ; s= xd1**2/(2.0d0*hxc**6)
 
 
  end select

    tl3dx4 = s

  end function
 
 
 
  real function tl3dx4z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= hxc**(-5)
 case(61) ; s= zd1/hxc**6
 case(36) ; s= 0.0d0
 case(60) ; s= yd1/hxc**6
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= xd1/hxc**6
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx4z = s

  end function
 
 
 
  real function tl3dx4z2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= hxc**(-6)
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx4z2 = s

  end function
 
 
 
  real function tl3dx4y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= hxc**(-5)
 case(60) ; s= zd1/hxc**6
 case(59) ; s= yd1/hxc**6
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= xd1/hxc**6
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx4y = s

  end function
 
 
 
  real function tl3dx4yz(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= hxc**(-6)
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
 
  end select

    tl3dx4yz = s

  end function
 
 
 
  real function tl3dx4y2(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= hxc**(-6)
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx4y2 = s

  end function
 
 
 
  real function tl3dx5(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= hxc**(-5)
 case(58) ; s= zd1/hxc**6
 case(57) ; s= yd1/hxc**6
 case(56) ; s= xd1/hxc**6
 
 
  end select

    tl3dx5= s

  end function
 
 
 
  real function tl3dx5z(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= hxc**(-6)
 case(57) ; s= 0.0d0
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx5z= s

  end function
 
 
 
  real function tl3dx5y(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= hxc**(-6)
 case(56) ; s= 0.0d0
 
 
  end select

    tl3dx5y = s

  end function
 
 
 
  real function tl3dx6(xd1,yd1,zd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1,zd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(3) ; s= 0.0d0
 case(9) ; s= 0.0d0
 case(19) ; s= 0.0d0
 case(34) ; s= 0.0d0
 case(55) ; s= 0.0d0
 case(83) ; s= 0.0d0
 case(2) ; s= 0.0d0
 case(8) ; s= 0.0d0
 case(18) ; s= 0.0d0
 case(33) ; s= 0.0d0
 case(54) ; s= 0.0d0
 case(82) ; s= 0.0d0
 case(7) ; s= 0.0d0
 case(17) ; s= 0.0d0
 case(32) ; s= 0.0d0
 case(53) ; s= 0.0d0
 case(81) ; s= 0.0d0
 case(16) ; s= 0.0d0
 case(31) ; s= 0.0d0
 case(52) ; s= 0.0d0
 case(80) ; s= 0.0d0
 case(30) ; s= 0.0d0
 case(51) ; s= 0.0d0
 case(79) ; s= 0.0d0
 case(50) ; s= 0.0d0
 case(78) ; s= 0.0d0
 case(77) ; s= 0.0d0
 case(1) ; s= 0.0d0
 case(6) ; s= 0.0d0
 case(14) ; s= 0.0d0
 case(27) ; s= 0.0d0
 case(49) ; s= 0.0d0
 case(76) ; s= 0.0d0
 case(5) ; s= 0.0d0
 case(15) ; s= 0.0d0
 case(29) ; s= 0.0d0
 case(48) ; s= 0.0d0
 case(75) ; s= 0.0d0
 case(13) ; s= 0.0d0
 case(28) ; s= 0.0d0
 case(47) ; s= 0.0d0
 case(74) ; s= 0.0d0
 case(26) ; s= 0.0d0
 case(46) ; s= 0.0d0
 case(73) ; s= 0.0d0
 case(45) ; s= 0.0d0
 case(72) ; s= 0.0d0
 case(71) ; s= 0.0d0
 case(4) ; s= 0.0d0
 case(12) ; s= 0.0d0
 case(24) ; s= 0.0d0
 case(44) ; s= 0.0d0
 case(70) ; s= 0.0d0
 case(11) ; s= 0.0d0
 case(25) ; s= 0.0d0
 case(43) ; s= 0.0d0
 case(69) ; s= 0.0d0
 case(23) ; s= 0.0d0
 case(42) ; s= 0.0d0
 case(68) ; s= 0.0d0
 case(41) ; s= 0.0d0
 case(67) ; s= 0.0d0
 case(66) ; s= 0.0d0
 case(10) ; s= 0.0d0
 case(22) ; s= 0.0d0
 case(39) ; s= 0.0d0
 case(65) ; s= 0.0d0
 case(21) ; s= 0.0d0
 case(40) ; s= 0.0d0
 case(64) ; s= 0.0d0
 case(38) ; s= 0.0d0
 case(63) ; s= 0.0d0
 case(62) ; s= 0.0d0
 case(20) ; s= 0.0d0
 case(37) ; s= 0.0d0
 case(61) ; s= 0.0d0
 case(36) ; s= 0.0d0
 case(60) ; s= 0.0d0
 case(59) ; s= 0.0d0
 case(35) ; s= 0.0d0
 case(58) ; s= 0.0d0
 case(57) ; s= 0.0d0
 case(56) ; s= hxc**(-6)

 end select

    tl3dx6 = s

  end function	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
real function df2dx(xd1,yd1,nderivative,iconsidered)
   implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s






   s=0.0d0
   select case(nderivative)
     case(1)
      s = 1.d0
     case(3)
      s = 2.d0*xd1
     case(4)
      s = yd1
     case(6)
      s = 3.d0*xd1**2
     case(7)
      s = 2.d0*xd1*yd1
     case(8)
      s = yd1**2
     case(10)
      s = 4.d0*xd1**3
     case(11)
      s = 3.d0*xd1**2*yd1
     case(12)
      s = 2.d0*xd1*yd1**2
     case(13)
      s = yd1**3
     case(15) 
      s = 5.d0*xd1**4
     case(16) 
      s = 4.d0*xd1**3*yd1
     case(17) 
      s = 3.d0*xd1**2*yd1**2
     case(18) 
      s = 2.d0*xd1*yd1**3
     case(19) 
      s = yd1**4
     case(21) 
      s = 6.d0*xd1**5
     case(22) 
      s = 5.d0*xd1**4*yd1
     case(23) 
      s = 4.d0*xd1**3*yd1**2
     case(24) 
      s = 3.d0*xd1**2*yd1**3
     case(25)
      s = 2.d0*xd1*yd1**4
     case(26) 
      s = yd1**5
    end select

    df2dx = s

  end function


 ! *****************************************************************************
  real function df2dy(xd1,yd1,nderivative,iconsidered)
     implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s



   s=0.0d0
   select case(nderivative)

     case(2) 
      s = 1.d0
     case(4) 
      s = xd1
     case(5) 
      s = 2.d0*yd1
     case(7) 
      s = xd1**2
     case(8) 
      s = 2.d0*xd1*yd1
     case(9) 
      s = 3.d0*yd1**2
     case(11) 
      s = xd1**3
     case(12)
      s = 2.d0*xd1**2*yd1
     case(13)
      s = 3.d0*xd1*yd1**2
     case(14) 
      s = 4.d0*yd1**3
     case(16) 
      s = xd1**4
     case(17) 
      s = 2.d0*xd1**3*yd1
     case(18) 
      s = 3.d0*xd1**2*yd1**2
     case(19) 
      s = 4.d0*xd1*yd1**3
     case(20) 
      s = 5.d0*yd1**4
     case(22) 
      s = xd1**5
     case(23) 
      s = 2.d0*xd1**4*yd1
     case(24) 
      s = 3.d0*xd1**3*yd1**2
     case(25) 
      s = 4.d0*xd1**2*yd1**3
     case(26) 
      s = 5.d0*xd1*yd1**4
     case(27) 
      s = 6.d0*yd1**5
    end select

    df2dy= s

  end function


 ! *****************************************************************************
  real function df2dx2(xd1,yd1,nderivative,iconsidered)
   implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(  3  ) 
      s = 2.d0
     case(  6  ) 
      s = 6.d0*xd1
     case(  7  ) 
      s = 2.d0*yd1
     case(  10  ) 
      s = 12.d0*xd1**2
     case(  11  ) 
      s = 6.d0*xd1*yd1
     case(12)
      s = 2.d0*yd1**2
     case(  15  ) 
      s = 20.d0*xd1**3
     case(  16  ) 
      s = 12.d0*xd1**2*yd1
     case(  17  ) 
      s = 6.d0*xd1*yd1**2
     case(  18  ) 
      s = 2.d0*yd1**3
     case(  21  ) 
      s = 30.d0*xd1**4
     case(  22  ) 
      s = 20.d0*xd1**3*yd1
     case(  23  ) 
      s = 12.d0*xd1**2*yd1**2
     case(  24  ) 
      s = 6.d0*xd1*yd1**3
     case(  25  ) 
      s = 2.d0*yd1**4
   end select

    df2dx2 = s

  end function



 ! *****************************************************************************
  real function df2dy2(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s = 0.0d0

   select case(nderivative)

    case(5) 
      s = 2.d0
    case(8) 
      s = 2.d0*xd1
    case(9) 
      s = 6.d0*yd1
    case(12)
      s = 2.d0*xd1**2
    case(13)
      s = 6.d0*xd1*yd1
    case(14)
      s = 12.d0*yd1**2
    case(17) 
      s = 2.d0*xd1**3
    case(18) 
      s = 6.d0*xd1**2*yd1
    case(19) 
      s = 12.d0*xd1*yd1**2
    case(20) 
      s = 20.d0*yd1**3
    case(23) 
      s = 2.d0*xd1**4
    case(24) 
      s = 6.d0*xd1**3*yd1
    case(25) 
      s = 12.d0*xd1**2*yd1**2
    case(26) 
      s = 20.d0*xd1*yd1**3
    case(27) 
      s = 30.d0*yd1**4

   end select
   df2dy2 = s

  end function


 ! *****************************************************************************
  real function df2dxy(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s


   s=0.0d0

   select case(nderivative)

    case( 4 ) 
      s = 1.d0
    case( 7 ) 
      s = 2.d0*xd1
    case( 8 ) 
      s = 2.d0*yd1
    case( 11 ) 
      s = 3.d0*xd1**2
    case(12)
      s = 4.d0*xd1*yd1
    case(13)
      s = 3.d0*yd1**2
    case( 16 ) 
      s = 4.d0*xd1**3
    case( 17 ) 
      s = 6.d0*xd1**2*yd1
    case( 18 ) 
      s = 6.d0*xd1*yd1**2
    case( 19 ) 
      s = 4.d0*yd1**3
    case( 22 ) 
      s = 5.d0*xd1**4
    case( 23 ) 
      s = 8.d0*xd1**3*yd1
    case( 24 ) 
      s = 9.d0*xd1**2*yd1**2
    case( 25 ) 
      s = 8.d0*xd1*yd1**3
    case( 26 ) 
      s = 5.d0*yd1**4
    end select

   df2dxy = s

  end function

! %%%%%%%%%%%%%%%%
 ! *****************************************************************************
  real function df2dx3(xd1,yd1,nderivative,iconsidered)
   implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s


   s=0.0d0

   select case(nderivative)
    case( 6 ) 
      s = 6.d0
    case( 10 ) 
      s = 24.d0*xd1
    case( 11 ) 
      s = 6.d0*yd1
    case( 15 ) 
      s = 60.d0*xd1**2
    case( 16 ) 
      s = 24.d0*xd1*yd1
    case( 17 ) 
      s = 6.d0*yd1**2
    case( 21 ) 
      s = 120.d0*xd1**3
    case( 22 ) 
      s = 60.d0*xd1**2*yd1
    case( 23 ) 
      s = 24.d0*xd1*yd1**2
    case( 24 ) 
      s = 6.d0*yd1**3
   end select

    df2dx3 = s
  end function


 ! *****************************************************************************
  real function df2dx2y(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s


   s=0.0d0
   select case(nderivative)
     case(7) 
      s = 2.d0
     case(11) 
      s = 6.d0*xd1
     case(12)
      s = 4.d0*yd1
     case(16) 
      s = 12.d0*xd1**2
     case(17) 
      s = 12.d0*xd1*yd1
     case(18) 
      s = 6.d0*yd1**2
     case(22) 
      s = 20.d0*xd1**3
     case(23) 
      s = 24.d0*xd1**2*yd1
     case(24) 
      s = 18.d0*xd1*yd1**2
     case(25) 
      s = 8.d0*yd1**3
   end select
   df2dx2y = s

  end function



 ! *****************************************************************************
  real function df2dxy2(xd1,yd1,nderivative,iconsidered)
   implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s  = 0.0d0
   select case(nderivative)

     case( 8  ) 
      s = 2.d0
     case(12)
      s = 4.d0*xd1
    case(13)
      s = 6.d0*yd1
     case( 17  ) 
      s = 6.d0*xd1**2
     case( 18  ) 
      s = 12*xd1*yd1
     case( 19  ) 
      s = 12.d0*yd1**2
     case( 23  ) 
      s = 8.d0*xd1**3
     case( 24  ) 
      s = 18.d0*xd1**2*yd1
     case( 25  ) 
       s = 24.d0*xd1*yd1**2
     case( 26  ) 
      s = 20.d0*yd1**3
   end select

   df2dxy2 = s

  end function


 ! *****************************************************************************
  real function df2dy3(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(  9  ) 
      s = 6.d0
     case(13)
      s = 6.d0*xd1
     case(  14  ) 
      s = 24.d0*yd1
     case(  18  ) 
      s = 6.d0*xd1**2
     case(  19  ) 
      s = 24.d0*xd1*yd1
     case(  20  ) 
      s = 60.d0*yd1**2
     case(  24  ) 
      s = 6.d0*xd1**3
     case(  25  ) 
      s = 24.d0*xd1**2*yd1
     case(  26  ) 
      s = 60.d0*xd1*yd1**2
     case(  27  ) 
      s = 120.d0*yd1**3
   end select

   df2dy3 = s

  end function



 ! *****************************************************************************
  real function df2dx4(xd1,yd1,nderivative,iconsidered)
     implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

    s=0.0d0
    select case(nderivative)
      case(10)
       s = 24.d0
      case(15)
       s = 120.d0*xd1
      case(16)
       s = 24.d0*yd1
      case(21)
       s = 360.d0*xd1**2
      case(22)
       s = 120.d0*xd1*yd1
      case(23)
       s = 24.d0*yd1**2
    end select

    df2dx4 = s
  end function




 ! *****************************************************************************
  real function df2dx3y(xd1,yd1,nderivative,iconsidered)
   implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0

   select case(nderivative)
     case(11)
      s = 6.d0
     case(16)
      s = 24.d0*xd1
     case(17)
      s = 12.d0*yd1
     case(22)
      s = 60.d0*xd1**2
     case(23)
      s = 48.d0*xd1*yd1
     case(24)
      s = 18.d0*yd1**2
   end select

   df2dx3y = s

  end function


 ! *****************************************************************************
  real function df2dx2y2(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(12)
      s = 4.d0
     case(17)
      s = 12.d0*xd1
     case(18)
      s = 12.d0*yd1
     case(23)
      s = 24.d0*xd1**2
     case(24)
      s = 36.d0*xd1*yd1
     case(25)
      s = 24.d0*yd1**2
   end select

   df2dx2y2 = s

  end function

 ! *****************************************************************************
  real function df2dxy3(xd1,yd1,nderivative,iconsidered)
     implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
    case(12)
      s = 6.d0
    case(18)
      s = 12.d0*xd1
    case(19)
      s = 24.d0*yd1
    case(24)
      s = 18.d0*xd1**2
    case(25)
      s = 48.d0*xd1*yd1
    case(26)
      s = 60.d0*yd1**2
   end select

   df2dxy3 = s

  end function

 ! *****************************************************************************
  real function df2dy4(xd1,yd1,nderivative,iconsidered)
    implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
    case(14)
      s = 24.d0
    case(19)
      s = 24.d0*xd1
    case(20)
      s = 120.d0*yd1
    case(25)
      s = 24.d0*xd1**2
    case(26)
      s = 120.d0*xd1*yd1
    case(27)
      s = 360.d0*yd1**2
  end select

   df2dy4 = s

  end function


! for  6th order weno

 ! *****************************************************************************
  real function df2dx5(xd1,yd1,nderivative,iconsidered)
    implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
    case(15)
      s = 120.d0
    case(21)
      s = 720.d0*xd1
    case(22)
      s = 120.d0*yd1
   end select

   df2dx5 = s

  end function


 ! *****************************************************************************
  real function df2dx4y(xd1,yd1,nderivative,iconsidered)
   implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(16)
      s = 24.d0
     case(22)
      s = 120.d0*xd1
     case(23)
      s = 48.d0*yd1
   end select

   df2dx4y = s

  end function


 ! *****************************************************************************
  real function df2dx3y2(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(17)
      s = 12.d0
     case(23)
      s = 48.d0*xd1
     case(24)
      s = 36.d0*yd1
    end select

   df2dx3y2 = s

  end function


 ! *****************************************************************************
  real function df2dx2y3(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
    case(18)
      s = 12.d0
    case(24)
      s = 36.d0*xd1
    case(25)
      s = 48.d0*yd1
   end select

   df2dx2y3 = s

  end function


 ! *****************************************************************************
  real function df2dxy4(xd1,yd1,nderivative,iconsidered)
    implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
    case(19)
      s = 24.d0
    case(25)
      s = 48.d0*xd1
    case(26)
      s = 120.d0*yd1
    end select
   df2dxy4 = s

  end function


 ! *****************************************************************************
  real function df2dy5(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
    case(20)
      s = 120.d0
    case(26)
      s = 120.d0*xd1
   case(27)
      s = 720.d0*yd1
   end select

   df2dy5 = s

  end function


 ! *****************************************************************************
  real function df2dx6(xd1,yd1,nderivative,iconsidered)
    implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(21)
      s = 720.d0
   end select

   df2dx6 = s

  end function


 ! *****************************************************************************
  real function df2dx5y(xd1,yd1,nderivative,iconsidered)
    implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(22)
      s = 120.d0
   end select

   df2dx5y = s

  end function


 ! *****************************************************************************
  real function df2dx4y2(xd1,yd1,nderivative,iconsidered)
   implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(23)
      s = 48.d0
   end select

   df2dx4y2 = s

  end function


 ! *****************************************************************************
  real function df2dx3y3(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
    case(24)
      s = 36.d0
   end select

   df2dx3y3 = s

  end function



 ! *****************************************************************************
  real function df2dx2y4(xd1,yd1,nderivative,iconsidered)
    implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
    case(25)
      s = 48.d0

   end select

   df2dx2y4 = s

  end function


 ! *****************************************************************************
  real function df2dxy5(xd1,yd1,nderivative,iconsidered)
    implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(26)
      s = 120.d0
   end select

  df2dxy5 = s

  end function


 ! *****************************************************************************
  real function df2dy6(xd1,yd1,nderivative,iconsidered)
  implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s

   s=0.0d0
   select case(nderivative)
     case(27)
      s = 720.d0
   end select

   df2dy6 = s

  end function	


  
  
  real function tl2dy(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=1/hxc
 case(5)  ; s=yd1/hxc**2
 case(9)  ; s=yd1**2/(2.0d0*hxc**3)
 case(14)  ; s=yd1**3/(6.0d0*hxc**4)
 case(20)  ; s=yd1**4/(24.0d0*hxc**5)
 case(27)  ; s=yd1**5/(120.0d0*hxc**6)
 case(1)  ; s=0.0d0
 case(4)  ; s=xd1/hxc**2
 case(8)  ; s=(xd1*yd1)/hxc**3
 case(13)  ; s=(xd1*yd1**2)/(2.0d0*hxc**4)
 case(19)  ; s=(xd1*yd1**3)/(6.0d0*hxc**5)
 case(26)  ; s=(xd1*yd1**4)/(24.0d0*hxc**6)
 case(3)  ; s=0.0d0
 case(7)  ; s=xd1**2/(2.0d0*hxc**3)
 case(12)  ; s=(xd1**2*yd1)/(2.0d0*hxc**4)
 case(18)  ; s=(xd1**2*yd1**2)/(4.0d0*hxc**5)
 case(25)  ; s=(xd1**2*yd1**3)/(12.0d0*hxc**6)
 case(6)  ; s=0.0d0
 case(11)  ; s=xd1**3/(6.0d0*hxc**4)
 case(17)  ; s=(xd1**3*yd1)/(6.0d0*hxc**5)
 case(24)  ; s=(xd1**3*yd1**2)/(12.0d0*hxc**6)
 case(10)  ; s=0.0d0
 case(16)  ; s=xd1**4/(24.0d0*hxc**5)
 case(23)  ; s=(xd1**4*yd1)/(24.0d0*hxc**6)
 case(15)  ; s=0.0d0
 case(22)  ; s=xd1**5/(120.0d0*hxc**6)
 case(21)  ; s=0.0d0
 
 
 
 
  end select

    tl2dy = s

  end function

real function  tl2dy2(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=hxc**(-2)
 case(9)  ; s=yd1/hxc**3
 case(14)  ; s=yd1**2/(2.0d0*hxc**4)
 case(20)  ; s=yd1**3/(6.0d0*hxc**5)
 case(27)  ; s=yd1**4/(24.0d0*hxc**6)
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=xd1/hxc**3
 case(13)  ; s=(xd1*yd1)/hxc**4
 case(19)  ; s=(xd1*yd1**2)/(2.0d0*hxc**5)
 case(26)  ; s=(xd1*yd1**3)/(6.0d0*hxc**6)
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=xd1**2/(2.0d0*hxc**4)
 case(18)  ; s=(xd1**2*yd1)/(2.0d0*hxc**5)
 case(25)  ; s=(xd1**2*yd1**2)/(4.0d0*hxc**6)
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=xd1**3/(6.0d0*hxc**5)
 case(24)  ; s=(xd1**3*yd1)/(6.0d0*hxc**6)
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=xd1**4/(24.0d0*hxc**6)
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
 
  end select

   tl2dy2 = s

  end function

real function  tl2dy3(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=hxc**(-3)
 case(14)  ; s=yd1/hxc**4
 case(20)  ; s=yd1**2/(2.0d0*hxc**5)
 case(27)  ; s=yd1**3/(6.0d0*hxc**6)
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=xd1/hxc**4
 case(19)  ; s=(xd1*yd1)/hxc**5
 case(26)  ; s=(xd1*yd1**2)/(2.0d0*hxc**6)
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=xd1**2/(2.0d0*hxc**5)
 case(25)  ; s=(xd1**2*yd1)/(2.0d0*hxc**6)
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=xd1**3/(6.0d0*hxc**6)
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dy3 = s

  end function

real function  tl2dy4(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=hxc**(-4)
 case(20)  ; s=yd1/hxc**5
 case(27)  ; s=yd1**2/(2.0d0*hxc**6)
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=xd1/hxc**5
 case(26)  ; s=(xd1*yd1)/hxc**6
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=xd1**2/(2.0d0*hxc**6)
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
  end select

    tl2dy4 = s

  end function

real function  tl2dy5(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=hxc**(-5)
 case(27)  ; s=yd1/hxc**6
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=xd1/hxc**6
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dy5 = s

  end function

real function  tl2dy6(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=hxc**(-6)
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dy6 = s

  end function

real function  tl2dx(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=1/hxc
 case(4)  ; s=yd1/hxc**2
 case(8)  ; s=yd1**2/(2.0d0*hxc**3)
 case(13)  ; s=yd1**3/(6.0d0*hxc**4)
 case(19)  ; s=yd1**4/(24.0d0*hxc**5)
 case(26)  ; s=yd1**5/(120.0d0*hxc**6)
 case(3)  ; s=xd1/hxc**2
 case(7)  ; s=(xd1*yd1)/hxc**3
 case(12)  ; s=(xd1*yd1**2)/(2.0d0*hxc**4)
 case(18)  ; s=(xd1*yd1**3)/(6.0d0*hxc**5)
 case(25)  ; s=(xd1*yd1**4)/(24.0d0*hxc**6)
 case(6)  ; s=xd1**2/(2.0d0*hxc**3)
 case(11)  ; s=(xd1**2*yd1)/(2.0d0*hxc**4)
 case(17)  ; s=(xd1**2*yd1**2)/(4.0d0*hxc**5)
 case(24)  ; s=(xd1**2*yd1**3)/(12.0d0*hxc**6)
 case(10)  ; s=xd1**3/(6.0d0*hxc**4)
 case(16)  ; s=(xd1**3*yd1)/(6.0d0*hxc**5)
 case(23)  ; s=(xd1**3*yd1**2)/(12.0d0*hxc**6)
 case(15)  ; s=xd1**4/(24.0d0*hxc**5)
 case(22)  ; s=(xd1**4*yd1)/(24.0d0*hxc**6)
 case(21)  ; s=xd1**5/(120.0d0*hxc**6)
 
 
 
  end select

   tl2dx = s

  end function

real function  tl2dxy(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=hxc**(-2)
 case(8)  ; s=yd1/hxc**3
 case(13)  ; s=yd1**2/(2.0d0*hxc**4)
 case(19)  ; s=yd1**3/(6.0d0*hxc**5)
 case(26)  ; s=yd1**4/(24.0d0*hxc**6)
 case(3)  ; s=0.0d0
 case(7)  ; s=xd1/hxc**3
 case(12)  ; s=(xd1*yd1)/hxc**4
 case(18)  ; s=(xd1*yd1**2)/(2.0d0*hxc**5)
 case(25)  ; s=(xd1*yd1**3)/(6.0d0*hxc**6)
 case(6)  ; s=0.0d0
 case(11)  ; s=xd1**2/(2.0d0*hxc**4)
 case(17)  ; s=(xd1**2*yd1)/(2.0d0*hxc**5)
 case(24)  ; s=(xd1**2*yd1**2)/(4.0d0*hxc**6)
 case(10)  ; s=0.0d0
 case(16)  ; s=xd1**3/(6.0d0*hxc**5)
 case(23)  ; s=(xd1**3*yd1)/(6.0d0*hxc**6)
 case(15)  ; s=0.0d0
 case(22)  ; s=xd1**4/(24.0d0*hxc**6)
 case(21)  ; s=0.0d0
 
 
  end select

    tl2dxy = s

  end function

real function  tl2dxy2(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=hxc**(-3)
 case(13)  ; s=yd1/hxc**4
 case(19)  ; s=yd1**2/(2.0d0*hxc**5)
 case(26)  ; s=yd1**3/(6.0d0*hxc**6)
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=xd1/hxc**4
 case(18)  ; s=(xd1*yd1)/hxc**5
 case(25)  ; s=(xd1*yd1**2)/(2.0d0*hxc**6)
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=xd1**2/(2.0d0*hxc**5)
 case(24)  ; s=(xd1**2*yd1)/(2.0d0*hxc**6)
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=xd1**3/(6.0d0*hxc**6)
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dxy2 = s

  end function

real function  tl2dxy3(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=hxc**(-4)
 case(19)  ; s=yd1/hxc**5
 case(26)  ; s=yd1**2/(2.0d0*hxc**6)
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=xd1/hxc**5
 case(25)  ; s=(xd1*yd1)/hxc**6
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=xd1**2/(2.0d0*hxc**6)
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
  end select

    tl2dxy3 = s

  end function

real function  tl2dxy4(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=hxc**(-5)
 case(26)  ; s=yd1/hxc**6
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=xd1/hxc**6
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
  end select

    tl2dxy4 = s

  end function

real function  tl2dxy5(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=hxc**(-6)
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
  end select

   tl2dxy5 = s

  end function

real function  tl2dx2(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=hxc**(-2)
 case(7)  ; s=yd1/hxc**3
 case(12)  ; s=yd1**2/(2.0d0*hxc**4)
 case(18)  ; s=yd1**3/(6.0d0*hxc**5)
 case(25)  ; s=yd1**4/(24.0d0*hxc**6)
 case(6)  ; s=xd1/hxc**3
 case(11)  ; s=(xd1*yd1)/hxc**4
 case(17)  ; s=(xd1*yd1**2)/(2.0d0*hxc**5)
 case(24)  ; s=(xd1*yd1**3)/(6.0d0*hxc**6)
 case(10)  ; s=xd1**2/(2.0d0*hxc**4)
 case(16)  ; s=(xd1**2*yd1)/(2.0d0*hxc**5)
 case(23)  ; s=(xd1**2*yd1**2)/(4.0d0*hxc**6)
 case(15)  ; s=xd1**3/(6.0d0*hxc**5)
 case(22)  ; s=(xd1**3*yd1)/(6.0d0*hxc**6)
 case(21)  ; s=xd1**4/(24.0d0*hxc**6)
 
 
 
  end select

    tl2dx2 = s

  end function

real function  tl2dx2y(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=hxc**(-3)
 case(12)  ; s=yd1/hxc**4
 case(18)  ; s=yd1**2/(2.0d0*hxc**5)
 case(25)  ; s=yd1**3/(6.0d0*hxc**6)
 case(6)  ; s=0.0d0
 case(11)  ; s=xd1/hxc**4
 case(17)  ; s=(xd1*yd1)/hxc**5
 case(24)  ; s=(xd1*yd1**2)/(2.0d0*hxc**6)
 case(10)  ; s=0.0d0
 case(16)  ; s=xd1**2/(2.0d0*hxc**5)
 case(23)  ; s=(xd1**2*yd1)/(2.0d0*hxc**6)
 case(15)  ; s=0.0d0
 case(22)  ; s=xd1**3/(6.0d0*hxc**6)
 case(21)  ; s=0.0d0
 
 
  end select

    tl2dx2y = s

  end function

real function  tl2dx2y2(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=hxc**(-4)
 case(18)  ; s=yd1/hxc**5
 case(25)  ; s=yd1**2/(2.0d0*hxc**6)
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=xd1/hxc**5
 case(24)  ; s=(xd1*yd1)/hxc**6
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=xd1**2/(2.0d0*hxc**6)
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dx2y2= s

  end function

real function  tl2dx2y3(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=hxc**(-5)
 case(25)  ; s=yd1/hxc**6
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=xd1/hxc**6
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
  end select

    tl2dx2y3 = s

  end function

real function  tl2dx2y4(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=hxc**(-6)
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
  end select

    tl2dx2y4 = s

  end function

real function  tl2dx3(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=hxc**(-3)
 case(11)  ; s=yd1/hxc**4
 case(17)  ; s=yd1**2/(2.0d0*hxc**5)
 case(24)  ; s=yd1**3/(6.0d0*hxc**6)
 case(10)  ; s=xd1/hxc**4
 case(16)  ; s=(xd1*yd1)/hxc**5
 case(23)  ; s=(xd1*yd1**2)/(2.0d0*hxc**6)
 case(15)  ; s=xd1**2/(2.0d0*hxc**5)
 case(22)  ; s=(xd1**2*yd1)/(2.0d0*hxc**6)
 case(21)  ; s=xd1**3/(6.0d0*hxc**6)
 
 
  end select

    tl2dx3 = s

  end function

real function  tl2dx3y(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=hxc**(-4)
 case(17)  ; s=yd1/hxc**5
 case(24)  ; s=yd1**2/(2.0d0*hxc**6)
 case(10)  ; s=0.0d0
 case(16)  ; s=xd1/hxc**5
 case(23)  ; s=(xd1*yd1)/hxc**6
 case(15)  ; s=0.0d0
 case(22)  ; s=xd1**2/(2.0d0*hxc**6)
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dx3y = s

  end function

real function  tl2dx3y2(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=hxc**(-5)
 case(24)  ; s=yd1/hxc**6
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=xd1/hxc**6
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dx3y2 = s

  end function

real function  tl2dx3y3(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=hxc**(-6)
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
  end select

    tl2dx3y3 = s

  end function

real function  tl2dx4(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=hxc**(-4)
 case(16)  ; s=yd1/hxc**5
 case(23)  ; s=yd1**2/(2.0d0*hxc**6)
 case(15)  ; s=xd1/hxc**5
 case(22)  ; s=(xd1*yd1)/hxc**6
 case(21)  ; s=xd1**2/(2.0d0*hxc**6)
 
 
 
  end select

    tl2dx4 = s

  end function

real function  tl2dx4y(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=hxc**(-5)
 case(23)  ; s=yd1/hxc**6
 case(15)  ; s=0.0d0
 case(22)  ; s=xd1/hxc**6
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dx4y= s

  end function

real function  tl2dx4y2(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=hxc**(-6)
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dx4y2 = s

  end function

real function  tl2dx5(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=hxc**(-5)
 case(22)  ; s=yd1/hxc**6
 case(21)  ; s=xd1/hxc**6
 
 
  end select

    tl2dx5 = s

  end function

real function  tl2dx5y(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=hxc**(-6)
 case(21)  ; s=0.0d0
 
 
 
  end select

    tl2dx5y = s

  end function

real function  tl2dx6(xd1,yd1,nderivative,iconsidered)
implicit none
#ifdef gpu
!$omp declare target
#endif
integer,intent(in)::nderivative,iconsidered
real,intent(in)::xd1,yd1
real::s,hxc   
hxc=sqrt(ielem_totvolume(iconsidered))




   s=0.0d0
   select case(nderivative)
 case(2)  ; s=0.0d0
 case(5)  ; s=0.0d0
 case(9)  ; s=0.0d0
 case(14)  ; s=0.0d0
 case(20)  ; s=0.0d0
 case(27)  ; s=0.0d0
 case(1)  ; s=0.0d0
 case(4)  ; s=0.0d0
 case(8)  ; s=0.0d0
 case(13)  ; s=0.0d0
 case(19)  ; s=0.0d0
 case(26)  ; s=0.0d0
 case(3)  ; s=0.0d0
 case(7)  ; s=0.0d0
 case(12)  ; s=0.0d0
 case(18)  ; s=0.0d0
 case(25)  ; s=0.0d0
 case(6)  ; s=0.0d0
 case(11)  ; s=0.0d0
 case(17)  ; s=0.0d0
 case(24)  ; s=0.0d0
 case(10)  ; s=0.0d0
 case(16)  ; s=0.0d0
 case(23)  ; s=0.0d0
 case(15)  ; s=0.0d0
 case(22)  ; s=0.0d0
 case(21)  ; s=hxc**(-6)

end select

    tl2dx6 = s

  end function
  
  
  
  
  
  
  

end module derivatives