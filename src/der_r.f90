MODULE DERIVATIVES
USE MPIINFO
USE DECLARATION
IMPLICIT NONE

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0

SELECTCASE (NDERIVATIVE)
CASE(1)  
                    
S=1.0d0

CASE(4)  
                    
S=2.0d0*XD1
CASE(5)  
                    
S=YD1

CASE(7)  
                    
S=ZD1


CASE(10)  
                    
S=3.d0*XD1**2
CASE(11)  
                    
S=2.d0*XD1*YD1
CASE(12)  
                    
S=YD1**2

CASE(14)  
                    
S=2.d0*XD1*ZD1
CASE(15)  
                    
S=YD1*ZD1

CASE(17)  
                    
S=ZD1**2


CASE(20)  
                    
S=4.d0*XD1**3
CASE(21)  
                    
S=3.d0*XD1**2*YD1
CASE(22)  
                    
S=2.d0*XD1*YD1**2
CASE(23)  
                    
S=YD1**3

CASE(25)  
                    
S=3.d0*XD1**2*ZD1
CASE(26)  
                    
S=2.d0*XD1*YD1*ZD1
CASE(27)  
                    
S=YD1**2*ZD1

CASE(29)  
                    
S=2.d0*XD1*ZD1**2
CASE(30)  
                    
S=YD1*ZD1**2

CASE(32)  
                    
S=ZD1**3


CASE(35)  
                    
S=5.d0*XD1**4
CASE(36)  
                    
S=4.d0*XD1**3*YD1
CASE(37)  
                    
S=3.d0*XD1**2*YD1**2
CASE(38)  
                    
S=2.d0*XD1*YD1**3
CASE(39)  
                    
S=YD1**4

CASE(41)  
                    
S=4.d0*XD1**3*ZD1
CASE(42)  
                    
S=3.d0*XD1**2*YD1*ZD1
CASE(43)  
                    
S=2.d0*XD1*YD1**2*ZD1
CASE(44)  
                    
S=YD1**3*ZD1

CASE(46)  
                    
S=3.d0*XD1**2*ZD1**2
CASE(47)  
                    
S=2.d0*XD1*YD1*ZD1**2
CASE(48)  
                    
S=YD1**2*ZD1**2

CASE(50)  
                    
S=2.d0*XD1*ZD1**3
CASE(51)  
                    
S=YD1*ZD1**3

CASE(53)  
                    
S=ZD1**4

CASE(56)  
                    
S=6.d0*XD1**5
CASE(57)  
                    
S=5.d0*XD1**4*YD1
CASE(58)  
                    
S=4.d0*XD1**3*YD1**2
CASE(59)  
                    
S=3.d0*XD1**2*YD1**3
CASE(60)  
                    
S=2.d0*XD1*YD1**4
CASE(61)  
                    
S=YD1**5

CASE(63)  
                    
S=5.d0*XD1**4*ZD1
CASE(64)  
                    
S=4.d0*XD1**3*YD1*ZD1
CASE(65)  
                    
S=3.d0*XD1**2*YD1**2*ZD1
CASE(66)  
                    
S=2.d0*XD1*YD1**3*ZD1
CASE(67)  
                    
S=YD1**4*ZD1

CASE(69)  
                    
S=4.d0*XD1**3*ZD1**2
CASE(70)  
                    
S=3.d0*XD1**2*YD1*ZD1**2
CASE(71)  
                    
S=2.d0*XD1*YD1**2*ZD1**2
CASE(72)  
                    
S=YD1**3*ZD1**2

CASE(74)  
                    
S=3.d0*XD1**2*ZD1**3
CASE(75)  
                    
S=2.d0*XD1*YD1*ZD1**3
CASE(76)  
                    
S=YD1**2*ZD1**3

CASE(78)  
                    
S=2.d0*XD1*ZD1**4
CASE(79)  
                    
S=YD1*ZD1**4

CASE(81)  
                    
S=ZD1**5

CASE DEFAULT

S=0.0

END SELECT
DFX=S
END FUNCTION DFX
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)

CASE(2) 
                    
S=1d0

CASE(5) 
                    
S=XD1
CASE(6) 
                    
S=2.d0*YD1

CASE(8) 
                    
S=ZD1


CASE(11) 
                    
S=XD1**2
CASE(12) 
                    
S=2.d0*XD1*YD1

case(13)
S=3.0D0*YD1**2

CASE(15) 
                    
S=XD1*ZD1
CASE(16) 
                    
S=2.d0*YD1*ZD1

CASE(18) 
                    
S=ZD1**2

CASE(21) 
                    
S=XD1**3
CASE(22) 
                    
S=2.d0*XD1**2*YD1
CASE(23) 
                    
S=3.d0*XD1*YD1**2
CASE(24) 
                    
S=4.d0*YD1**3

CASE(26) 
                    
S=XD1**2*ZD1
CASE(27) 
                    
S=2.d0*XD1*YD1*ZD1
CASE(28) 
                    
S=3.d0*YD1**2*ZD1

CASE(30) 
                    
S=XD1*ZD1**2
CASE(31) 
                    
S=2.d0*YD1*ZD1**2

CASE(33) 
                    
S=ZD1**3

CASE(36) 
                    
S=XD1**4
CASE(37) 
                    
S=2.d0*XD1**3*YD1
CASE(38) 
                    
S=3.d0*XD1**2*YD1**2
CASE(39) 
                    
S=4.d0*XD1*YD1**3
CASE(40) 
                    
S=5.d0*YD1**4

CASE(42) 
                    
S=XD1**3*ZD1
CASE(43) 
                    
S=2.d0*XD1**2*YD1*ZD1
CASE(44) 
                    
S=3.d0*XD1*YD1**2*ZD1
CASE(45) 
                    
S=4.d0*YD1**3*ZD1

CASE(47) 
                    
S=XD1**2*ZD1**2
CASE(48) 
                    
S=2.d0*XD1*YD1*ZD1**2
CASE(49) 
                    
S=3.d0*YD1**2*ZD1**2

CASE(51) 
                    
S=XD1*ZD1**3
CASE(52) 
                    
S=2.d0*YD1*ZD1**3

CASE(54) 
                    
S=ZD1**4

CASE(57) 
                    
S=XD1**5
CASE(58) 
                    
S=2.d0*XD1**4*YD1
CASE(59) 
                    
S=3.d0*XD1**3*YD1**2
CASE(60) 
                    
S=4.d0*XD1**2*YD1**3
CASE(61) 
                    
S=5.d0*XD1*YD1**4
CASE(62) 
                    
S=6.d0*YD1**5

CASE(64) 
                    
S=XD1**4*ZD1
CASE(65) 
                    
S=2.d0*XD1**3*YD1*ZD1
CASE(66) 
                    
S=3.d0*XD1**2*YD1**2*ZD1
CASE(67) 
                    
S=4.d0*XD1*YD1**3*ZD1
CASE(68) 
                    
S=5.d0*YD1**4*ZD1

CASE(70) 
                    
S=XD1**3*ZD1**2
CASE(71) 
                    
S=2.d0*XD1**2*YD1*ZD1**2
CASE(72) 
                    
S=3.d0*XD1*YD1**2*ZD1**2
CASE(73) 
                    
S=4.d0*YD1**3*ZD1**2

CASE(75) 
                    
S=XD1**2*ZD1**3
CASE(76) 
                    
S=2.d0*XD1*YD1*ZD1**3
CASE(77) 
                    
S=3.d0*YD1**2*ZD1**3

CASE(79) 
                    
S=XD1*ZD1**4
CASE(80) 
                    
S=2.d0*YD1*ZD1**4

CASE(82) 
                    
S=ZD1**5


CASE DEFAULT
S=0.0D0



END SELECT
DFY=S
END FUNCTION DFY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)

CASE (3)  
                    
S=1

CASE (7)  
                    
S=XD1
CASE (8)  
                    
S=YD1
CASE (9)  
                    
S=2.d0*ZD1


CASE (14)  
                    
S=XD1**2
CASE (15)  
                    
S=XD1*YD1
CASE (16)  
                    
S=YD1**2
CASE (17)  
                    
S=2.d0*XD1*ZD1
CASE (18)  
                    
S=2.d0*YD1*ZD1
CASE (19)  
                    
S=3.d0*ZD1**2

CASE (25)  
                    
S=XD1**3
CASE (26)  
                    
S=XD1**2*YD1
CASE (27)  
                    
S=XD1*YD1**2
CASE (28)  
                    
S=YD1**3
CASE (29)  
                    
S=2.d0*XD1**2*ZD1
CASE (30)  
                    
S=2.d0*XD1*YD1*ZD1
CASE (31)  
                    
S=2.d0*YD1**2*ZD1
CASE (32)  
                    
S=3.d0*XD1*ZD1**2
CASE (33)  
                    
S=3.d0*YD1*ZD1**2
CASE (34)  
                    
S=4.d0*ZD1**3

                    
S=XD1**4
CASE (42)  
                    
S=XD1**3*YD1
CASE (43)  
                    
S=XD1**2*YD1**2
CASE (44)  
                    
S=XD1*YD1**3
CASE (45)  
                    
S=YD1**4
CASE (46)  
                    
S=2.d0*XD1**3*ZD1
CASE (47)  
                    
S=2.d0*XD1**2*YD1*ZD1
CASE (48)  
                    
S=2.d0*XD1*YD1**2*ZD1
CASE (49)  
                    
S=2.d0*YD1**3*ZD1
CASE (50)  
                    
S=3.d0*XD1**2*ZD1**2
CASE (51)  
                    
S=3.d0*XD1*YD1*ZD1**2
CASE (52)  
                    
S=3.d0*YD1**2*ZD1**2
CASE (53)  
                    
S=4.d0*XD1*ZD1**3
CASE (54)  
                    
S=4.d0*YD1*ZD1**3
CASE (55)  
                    
S=5.d0*ZD1**4

                    
S=XD1**5
CASE (64)  
                    
S=XD1**4*YD1
CASE (65)  
                    
S=XD1**3*YD1**2
CASE (66)  
                    
S=XD1**2*YD1**3
CASE (67)  
                    
S=XD1*YD1**4
CASE (68)  
                    
S=YD1**5
CASE (69)  
                    
S=2.d0*XD1**4*ZD1
CASE (70)  
                    
S=2.d0*XD1**3*YD1*ZD1
CASE (71)  
                    
S=2.d0*XD1**2*YD1**2*ZD1
CASE (72)  
                    
S=2.d0*XD1*YD1**3*ZD1
CASE (73)  
                    
S=2.d0*YD1**4*ZD1
CASE (74)  
                    
S=3.d0*XD1**3*ZD1**2
CASE (75)  
                    
S=3.d0*XD1**2*YD1*ZD1**2
CASE (76)  
                    
S=3.d0*XD1*YD1**2*ZD1**2
CASE (77)  
                    
S=3.d0*YD1**3*ZD1**2
CASE (78)  
                    
S=4.d0*XD1**2*ZD1**3
CASE (79)  
                    
S=4.d0*XD1*YD1*ZD1**3
CASE (80)  
                    
S=4.d0*YD1**2*ZD1**3
CASE (81)  
                    
S=5.d0*XD1*ZD1**4
CASE (82)  
                    
S=5.d0*YD1*ZD1**4
CASE (83)  
                    
S=6.d0*ZD1**5

case default
s=0.0d0

END SELECT
DFZ=S
END FUNCTION DFZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=2
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=6.d0*XD1
CASE(11) 
                    
S=2.d0*YD1
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=2.d0*ZD1
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=12.d0*XD1**2
CASE(21) 
                    
S=6.d0*XD1*YD1
CASE(22) 
                    
S=2.d0*YD1**2
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=6.d0*XD1*ZD1
CASE(26) 
                    
S=2.d0*YD1*ZD1
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=2.d0*ZD1**2
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=20*XD1**3
CASE(36) 
                    
S=12.d0*XD1**2*YD1
CASE(37) 
                    
S=6.d0*XD1*YD1**2
CASE(38) 
                    
S=2.d0*YD1**3
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=12.d0*XD1**2*ZD1
CASE(42) 
                    
S=6.d0*XD1*YD1*ZD1
CASE(43) 
                    
S=2.d0*YD1**2*ZD1
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=6.d0*XD1*ZD1**2
CASE(47) 
                    
S=2.d0*YD1*ZD1**2
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=2.d0*ZD1**3
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=30*XD1**4
CASE(57) 
                    
S=20*XD1**3*YD1
CASE(58) 
                    
S=12.d0*XD1**2*YD1**2
CASE(59) 
                    
S=6.d0*XD1*YD1**3
CASE(60) 
                    
S=2.d0*YD1**4
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=20*XD1**3*ZD1
CASE(64) 
                    
S=12.d0*XD1**2*YD1*ZD1
CASE(65) 
                    
S=6.d0*XD1*YD1**2*ZD1
CASE(66) 
                    
S=2.d0*YD1**3*ZD1
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=12.d0*XD1**2*ZD1**2
CASE(70) 
                    
S=6.d0*XD1*YD1*ZD1**2
CASE(71) 
                    
S=2.d0*YD1**2*ZD1**2
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=6.d0*XD1*ZD1**3
CASE(75) 
                    
S=2.d0*YD1*ZD1**3
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=2.d0*ZD1**4
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2=S
END FUNCTION DFX2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=2
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=2.d0*XD1
CASE(13) 
                    
S=6.d0*YD1
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=2.d0*ZD1
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=2.d0*XD1**2
CASE(23) 
                    
S=6.d0*XD1*YD1
CASE(24) 
                    
S=12.d0*YD1**2
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=2.d0*XD1*ZD1
CASE(28) 
                    
S=6.d0*YD1*ZD1
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=2.d0*ZD1**2
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=2.d0*XD1**3
CASE(38) 
                    
S=6.d0*XD1**2*YD1
CASE(39) 
                    
S=12.d0*XD1*YD1**2
CASE(40) 
                    
S=20*YD1**3
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=2.d0*XD1**2*ZD1
CASE(44) 
                    
S=6.d0*XD1*YD1*ZD1
CASE(45) 
                    
S=12.d0*YD1**2*ZD1
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=2.d0*XD1*ZD1**2
CASE(49) 
                    
S=6.d0*YD1*ZD1**2
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=2.d0*ZD1**3
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=2.d0*XD1**4
CASE(59) 
                    
S=6.d0*XD1**3*YD1
CASE(60) 
                    
S=12.d0*XD1**2*YD1**2
CASE(61) 
                    
S=20*XD1*YD1**3
CASE(62) 
                    
S=30*YD1**4
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=2.d0*XD1**3*ZD1
CASE(66) 
                    
S=6.d0*XD1**2*YD1*ZD1
CASE(67) 
                    
S=12.d0*XD1*YD1**2*ZD1
CASE(68) 
                    
S=20*YD1**3*ZD1
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=2.d0*XD1**2*ZD1**2
CASE(72) 
                    
S=6.d0*XD1*YD1*ZD1**2
CASE(73) 
                    
S=12.d0*YD1**2*ZD1**2
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=2.d0*XD1*ZD1**3
CASE(77) 
                    
S=6.d0*YD1*ZD1**3
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=2.d0*ZD1**4
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY2=S
END FUNCTION DFY2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1)  
                    
S=0.0D0
CASE (2)  
                    
S=0.0D0
CASE (3)  
                    
S=0.0D0
CASE (4)  
                    
S=0.0D0
CASE (5)  
                    
S=0.0D0
CASE (6)  
                    
S=0.0D0
CASE (7)  
                    
S=0.0D0
CASE (8)  
                    
S=0.0D0
CASE (9)  
                    
S=2
CASE (10)  
                    
S=0.0D0
CASE (11)  
                    
S=0.0D0
CASE (12)  
                    
S=0.0D0
CASE (13)  
                    
S=0.0D0
CASE (14)  
                    
S=0.0D0
CASE (15)  
                    
S=0.0D0
CASE (16)  
                    
S=0.0D0
CASE (17)  
                    
S=2.d0*XD1
CASE (18)  
                    
S=2.d0*YD1
CASE (19)  
                    
S=6.d0*ZD1
CASE (20)  
                    
S=0.0D0
CASE (21)  
                    
S=0.0D0
CASE (22)  
                    
S=0.0D0
CASE (23)  
                    
S=0.0D0
CASE (24)  
                    
S=0.0D0
CASE (25)  
                    
S=0.0D0
CASE (26)  
                    
S=0.0D0
CASE (27)  
                    
S=0.0D0
CASE (28)  
                    
S=0.0D0
CASE (29)  
                    
S=2.d0*XD1**2
CASE (30)  
                    
S=2.d0*XD1*YD1
CASE (31)  
                    
S=2.d0*YD1**2
CASE (32)  
                    
S=6.d0*XD1*ZD1
CASE (33)  
                    
S=6.d0*YD1*ZD1
CASE (34)  
                    
S=12.d0*ZD1**2
CASE (35)  
                    
S=0.0D0
CASE (36)  
                    
S=0.0D0
CASE (37)  
                    
S=0.0D0
CASE (38)  
                    
S=0.0D0
CASE (39)  
                    
S=0.0D0
CASE (40)  
                    
S=0.0D0
CASE (41)  
                    
S=0.0D0
CASE (42)  
                    
S=0.0D0
CASE (43)  
                    
S=0.0D0
CASE (44)  
                    
S=0.0D0
CASE (45)  
                    
S=0.0D0
CASE (46)  
                    
S=2.d0*XD1**3
CASE (47)  
                    
S=2.d0*XD1**2*YD1
CASE (48)  
                    
S=2.d0*XD1*YD1**2
CASE (49)  
                    
S=2.d0*YD1**3
CASE (50)  
                    
S=6.d0*XD1**2*ZD1
CASE (51)  
                    
S=6.d0*XD1*YD1*ZD1
CASE (52)  
                    
S=6.d0*YD1**2*ZD1
CASE (53)  
                    
S=12.d0*XD1*ZD1**2
CASE (54)  
                    
S=12.d0*YD1*ZD1**2
CASE (55)  
                    
S=20*ZD1**3
CASE (56)  
                    
S=0.0D0
CASE (57)  
                    
S=0.0D0
CASE (58)  
                    
S=0.0D0
CASE (59)  
                    
S=0.0D0
CASE (60)  
                    
S=0.0D0
CASE (61)  
                    
S=0.0D0
CASE (62)  
                    
S=0.0D0
CASE (63)  
                    
S=0.0D0
CASE (64)  
                    
S=0.0D0
CASE (65)  
                    
S=0.0D0
CASE (66)  
                    
S=0.0D0
CASE (67)  
                    
S=0.0D0
CASE (68)  
                    
S=0.0D0
CASE (69)  
                    
S=2.d0*XD1**4
CASE (70)  
                    
S=2.d0*XD1**3*YD1
CASE (71)  
                    
S=2.d0*XD1**2*YD1**2
CASE (72)  
                    
S=2.d0*XD1*YD1**3
CASE (73)  
                    
S=2.d0*YD1**4
CASE (74)  
                    
S=6.d0*XD1**3*ZD1
CASE (75)  
                    
S=6.d0*XD1**2*YD1*ZD1
CASE (76)  
                    
S=6.d0*XD1*YD1**2*ZD1
CASE (77)  
                    
S=6.d0*YD1**3*ZD1
CASE (78)  
                    
S=12.d0*XD1**2*ZD1**2
CASE (79)  
                    
S=12.d0*XD1*YD1*ZD1**2
CASE (80)  
                    
S=12.d0*YD1**2*ZD1**2
CASE (81)  
                    
S=20*XD1*ZD1**3
CASE (82)  
                    
S=20*YD1*ZD1**3
CASE (83)  
                    
S=30*ZD1**4
END SELECT
DFZ2=S
END FUNCTION DFZ2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=1
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=2.d0*XD1
CASE(12) 
                    
S=2.d0*YD1
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=ZD1
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=3.d0*XD1**2
CASE(22) 
                    
S=4.d0*XD1*YD1
CASE(23) 
                    
S=3.d0*YD1**2
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=2.d0*XD1*ZD1
CASE(27) 
                    
S=2.d0*YD1*ZD1
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=ZD1**2
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=4.d0*XD1**3
CASE(37) 
                    
S=6.d0*XD1**2*YD1
CASE(38) 
                    
S=6.d0*XD1*YD1**2
CASE(39) 
                    
S=4.d0*YD1**3
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=3.d0*XD1**2*ZD1
CASE(43) 
                    
S=4.d0*XD1*YD1*ZD1
CASE(44) 
                    
S=3.d0*YD1**2*ZD1
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=2.d0*XD1*ZD1**2
CASE(48) 
                    
S=2.d0*YD1*ZD1**2
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=ZD1**3
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=5.d0*XD1**4
CASE(58) 
                    
S=8.d0*XD1**3*YD1
CASE(59) 
                    
S=9*XD1**2*YD1**2
CASE(60) 
                    
S=8.d0*XD1*YD1**3
CASE(61) 
                    
S=5.d0*YD1**4
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=4.d0*XD1**3*ZD1
CASE(65) 
                    
S=6.d0*XD1**2*YD1*ZD1
CASE(66) 
                    
S=6.d0*XD1*YD1**2*ZD1
CASE(67) 
                    
S=4.d0*YD1**3*ZD1
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=3.d0*XD1**2*ZD1**2
CASE(71) 
                    
S=4.d0*XD1*YD1*ZD1**2
CASE(72) 
                    
S=3.d0*YD1**2*ZD1**2
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=2.d0*XD1*ZD1**3
CASE(76) 
                    
S=2.d0*YD1*ZD1**3
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=ZD1**4
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY=S
END FUNCTION DFXY
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFYZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1)  
                    
S=0.0D0
CASE(2)  
                    
S=0.0D0
CASE(3)  
                    
S=0.0D0
CASE(4)  
                    
S=0.0D0
CASE(5)  
                    
S=0.0D0
CASE(6)  
                    
S=0.0D0
CASE(7)  
                    
S=0.0D0
CASE(8)  
                    
S=1
CASE(9)  
                    
S=0.0D0
CASE(10)  
                    
S=0.0D0
CASE(11)  
                    
S=0.0D0
CASE(12)  
                    
S=0.0D0
CASE(13)  
                    
S=0.0D0
CASE(14)  
                    
S=0.0D0
CASE(15)  
                    
S=XD1
CASE(16)  
                    
S=2.d0*YD1
CASE(17)  
                    
S=0.0D0
CASE(18)  
                    
S=2.d0*ZD1
CASE(19)  
                    
S=0.0D0
CASE(20)  
                    
S=0.0D0
CASE(21)  
                    
S=0.0D0
CASE(22)  
                    
S=0.0D0
CASE(23)  
                    
S=0.0D0
CASE(24)  
                    
S=0.0D0
CASE(25)  
                    
S=0.0D0
CASE(26)  
                    
S=XD1**2
CASE(27)  
                    
S=2.d0*XD1*YD1
CASE(28)  
                    
S=3.d0*YD1**2
CASE(29)  
                    
S=0.0D0
CASE(30)  
                    
S=2.d0*XD1*ZD1
CASE(31)  
                    
S=4.d0*YD1*ZD1
CASE(32)  
                    
S=0.0D0
CASE(33)  
                    
S=3.d0*ZD1**2
CASE(34)  
                    
S=0.0D0
CASE(35)  
                    
S=0.0D0
CASE(36)  
                    
S=0.0D0
CASE(37)  
                    
S=0.0D0
CASE(38)  
                    
S=0.0D0
CASE(39)  
                    
S=0.0D0
CASE(40)  
                    
S=0.0D0
CASE(41)  
                    
S=0.0D0
CASE(42)  
                    
S=XD1**3
CASE(43)  
                    
S=2.d0*XD1**2*YD1
CASE(44)  
                    
S=3.d0*XD1*YD1**2
CASE(45)  
                    
S=4.d0*YD1**3
CASE(46)  
                    
S=0.0D0
CASE(47)  
                    
S=2.d0*XD1**2*ZD1
CASE(48)  
                    
S=4.d0*XD1*YD1*ZD1
CASE(49)  
                    
S=6.d0*YD1**2*ZD1
CASE(50)  
                    
S=0.0D0
CASE(51)  
                    
S=3.d0*XD1*ZD1**2
CASE(52)  
                    
S=6.d0*YD1*ZD1**2
CASE(53)  
                    
S=0.0D0
CASE(54)  
                    
S=4.d0*ZD1**3
CASE(55)  
                    
S=0.0D0
CASE(56)  
                    
S=0.0D0
CASE(57)  
                    
S=0.0D0
CASE(58)  
                    
S=0.0D0
CASE(59)  
                    
S=0.0D0
CASE(60)  
                    
S=0.0D0
CASE(61)  
                    
S=0.0D0
CASE(62)  
                    
S=0.0D0
CASE(63)  
                    
S=0.0D0
CASE(64)  
                    
S=XD1**4
CASE(65)  
                    
S=2.d0*XD1**3*YD1
CASE(66)  
                    
S=3.d0*XD1**2*YD1**2
CASE(67)  
                    
S=4.d0*XD1*YD1**3
CASE(68)  
                    
S=5.d0*YD1**4
CASE(69)  
                    
S=0.0D0
CASE(70)  
                    
S=2.d0*XD1**3*ZD1
CASE(71)  
                    
S=4.d0*XD1**2*YD1*ZD1
CASE(72)  
                    
S=6.d0*XD1*YD1**2*ZD1
CASE(73)  
                    
S=8.d0*YD1**3*ZD1
CASE(74)  
                    
S=0.0D0
CASE(75)  
                    
S=3.d0*XD1**2*ZD1**2
CASE(76)  
                    
S=6.d0*XD1*YD1*ZD1**2
CASE(77)  
                    
S=9*YD1**2*ZD1**2
CASE(78)  
                    
S=0.0D0
CASE(79)  
                    
S=4.d0*XD1*ZD1**3
CASE(80)  
                    
S=8.d0*YD1*ZD1**3
CASE(81)  
                    
S=0.0D0
CASE(82)  
                    
S=5.d0*ZD1**4
CASE(83)  
                    
S=0.0D0
END SELECT
DFYZ=S
END FUNCTION DFYZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=1
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=2.d0*XD1
CASE(15) 
                    
S=YD1
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=2.d0*ZD1
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=3.d0*XD1**2
CASE(26) 
                    
S=2.d0*XD1*YD1
CASE(27) 
                    
S=YD1**2
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=4.d0*XD1*ZD1
CASE(30) 
                    
S=2.d0*YD1*ZD1
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=3.d0*ZD1**2
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=4.d0*XD1**3
CASE(42) 
                    
S=3.d0*XD1**2*YD1
CASE(43) 
                    
S=2.d0*XD1*YD1**2
CASE(44) 
                    
S=YD1**3
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=6.d0*XD1**2*ZD1
CASE(47) 
                    
S=4.d0*XD1*YD1*ZD1
CASE(48) 
                    
S=2.d0*YD1**2*ZD1
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=6.d0*XD1*ZD1**2
CASE(51) 
                    
S=3.d0*YD1*ZD1**2
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=4.d0*ZD1**3
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=5.d0*XD1**4
CASE(64) 
                    
S=4.d0*XD1**3*YD1
CASE(65) 
                    
S=3.d0*XD1**2*YD1**2
CASE(66) 
                    
S=2.d0*XD1*YD1**3
CASE(67) 
                    
S=YD1**4
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=8.d0*XD1**3*ZD1
CASE(70) 
                    
S=6.d0*XD1**2*YD1*ZD1
CASE(71) 
                    
S=4.d0*XD1*YD1**2*ZD1
CASE(72) 
                    
S=2.d0*YD1**3*ZD1
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=9*XD1**2*ZD1**2
CASE(75) 
                    
S=6.d0*XD1*YD1*ZD1**2
CASE(76) 
                    
S=3.d0*YD1**2*ZD1**2
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=8.d0*XD1*ZD1**3
CASE(79) 
                    
S=4.d0*YD1*ZD1**3
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=5.d0*ZD1**4
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXZ=S
END FUNCTION DFXZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=6
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=24.d0*XD1
CASE(21) 
                    
S=6.d0*YD1
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=6.d0*ZD1
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=60*XD1**2
CASE(36) 
                    
S=24.d0*XD1*YD1
CASE(37) 
                    
S=6.d0*YD1**2
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=24.d0*XD1*ZD1
CASE(42) 
                    
S=6.d0*YD1*ZD1
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=6.d0*ZD1**2
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=120*XD1**3
CASE(57) 
                    
S=60*XD1**2*YD1
CASE(58) 
                    
S=24.d0*XD1*YD1**2
CASE(59) 
                    
S=6.d0*YD1**3
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=60*XD1**2*ZD1
CASE(64) 
                    
S=24.d0*XD1*YD1*ZD1
CASE(65) 
                    
S=6.d0*YD1**2*ZD1
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=24.d0*XD1*ZD1**2
CASE(70) 
                    
S=6.d0*YD1*ZD1**2
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=6.d0*ZD1**3
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3=S
END FUNCTION DFX3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=2
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=6.d0*XD1
CASE(22) 
                    
S=4.d0*YD1
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=2.d0*ZD1
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=12.d0*XD1**2
CASE(37) 
                    
S=12.d0*XD1*YD1
CASE(38) 
                    
S=6.d0*YD1**2
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=6.d0*XD1*ZD1
CASE(43) 
                    
S=4.d0*YD1*ZD1
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=2.d0*ZD1**2
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=20*XD1**3
CASE(58) 
                    
S=24.d0*XD1**2*YD1
CASE(59) 
                    
S=18.d0*XD1*YD1**2
CASE(60) 
                    
S=8.d0*YD1**3
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=12.d0*XD1**2*ZD1
CASE(65) 
                    
S=12.d0*XD1*YD1*ZD1
CASE(66) 
                    
S=6.d0*YD1**2*ZD1
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=6.d0*XD1*ZD1**2
CASE(71) 
                    
S=4.d0*YD1*ZD1**2
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=2.d0*ZD1**3
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Y=S
END FUNCTION DFX2Y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=2
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=4.d0*XD1
CASE(23) 
                    
S=6.d0*YD1
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=2.d0*ZD1
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=6.d0*XD1**2
CASE(38) 
                    
S=12.d0*XD1*YD1
CASE(39) 
                    
S=12.d0*YD1**2
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=4.d0*XD1*ZD1
CASE(44) 
                    
S=6.d0*YD1*ZD1
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=2.d0*ZD1**2
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=8.d0*XD1**3
CASE(59) 
                    
S=18.d0*XD1**2*YD1
CASE(60) 
                    
S=24.d0*XD1*YD1**2
CASE(61) 
                    
S=20*YD1**3
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=6.d0*XD1**2*ZD1
CASE(66) 
                    
S=12.d0*XD1*YD1*ZD1
CASE(67) 
                    
S=12.d0*YD1**2*ZD1
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=4.d0*XD1*ZD1**2
CASE(72) 
                    
S=6.d0*YD1*ZD1**2
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=2.d0*ZD1**3
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY2=S
END FUNCTION DFXY2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=6
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=6.d0*XD1
CASE(24) 
                    
S=24.d0*YD1
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=6.d0*ZD1
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=6.d0*XD1**2
CASE(39) 
                    
S=24.d0*XD1*YD1
CASE(40) 
                    
S=60*YD1**2
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=6.d0*XD1*ZD1
CASE(45) 
                    
S=24.d0*YD1*ZD1
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=6.d0*ZD1**2
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=6.d0*XD1**3
CASE(60) 
                    
S=24.d0*XD1**2*YD1
CASE(61) 
                    
S=60*XD1*YD1**2
CASE(62) 
                    
S=120*YD1**3
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=6.d0*XD1**2*ZD1
CASE(67) 
                    
S=24.d0*XD1*YD1*ZD1
CASE(68) 
                    
S=60*YD1**2*ZD1
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=6.d0*XD1*ZD1**2
CASE(73) 
                    
S=24.d0*YD1*ZD1**2
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=6.d0*ZD1**3
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY3=S
END FUNCTION DFY3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=2
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=6.d0*XD1
CASE(26) 
                    
S=2.d0*YD1
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=4.d0*ZD1
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=12.d0*XD1**2
CASE(42) 
                    
S=6.d0*XD1*YD1
CASE(43) 
                    
S=2.d0*YD1**2
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=12.d0*XD1*ZD1
CASE(47) 
                    
S=4.d0*YD1*ZD1
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=6.d0*ZD1**2
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=20*XD1**3
CASE(64) 
                    
S=12.d0*XD1**2*YD1
CASE(65) 
                    
S=6.d0*XD1*YD1**2
CASE(66) 
                    
S=2.d0*YD1**3
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=24.d0*XD1**2*ZD1
CASE(70) 
                    
S=12.d0*XD1*YD1*ZD1
CASE(71) 
                    
S=4.d0*YD1**2*ZD1
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=18.d0*XD1*ZD1**2
CASE(75) 
                    
S=6.d0*YD1*ZD1**2
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=8.d0*ZD1**3
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Z=S
END FUNCTION DFX2Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXYZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=1
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=2.d0*XD1
CASE(27) 
                    
S=2.d0*YD1
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=2.d0*ZD1
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=3.d0*XD1**2
CASE(43) 
                    
S=4.d0*XD1*YD1
CASE(44) 
                    
S=3.d0*YD1**2
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=4.d0*XD1*ZD1
CASE(48) 
                    
S=4.d0*YD1*ZD1
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=3.d0*ZD1**2
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=4.d0*XD1**3
CASE(65) 
                    
S=6.d0*XD1**2*YD1
CASE(66) 
                    
S=6.d0*XD1*YD1**2
CASE(67) 
                    
S=4.d0*YD1**3
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=6.d0*XD1**2*ZD1
CASE(71) 
                    
S=8.d0*XD1*YD1*ZD1
CASE(72) 
                    
S=6.d0*YD1**2*ZD1
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=6.d0*XD1*ZD1**2
CASE(76) 
                    
S=6.d0*YD1*ZD1**2
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=4.d0*ZD1**3
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXYZ=S
END FUNCTION DFXYZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=2
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=2.d0*XD1
CASE(28) 
                    
S=6.d0*YD1
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=4.d0*ZD1
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=2.d0*XD1**2
CASE(44) 
                    
S=6.d0*XD1*YD1
CASE(45) 
                    
S=12.d0*YD1**2
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=4.d0*XD1*ZD1
CASE(49) 
                    
S=12.d0*YD1*ZD1
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=6.d0*ZD1**2
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=2.d0*XD1**3
CASE(66) 
                    
S=6.d0*XD1**2*YD1
CASE(67) 
                    
S=12.d0*XD1*YD1**2
CASE(68) 
                    
S=20*YD1**3
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=4.d0*XD1**2*ZD1
CASE(72) 
                    
S=12.d0*XD1*YD1*ZD1
CASE(73) 
                    
S=24.d0*YD1**2*ZD1
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=6.d0*XD1*ZD1**2
CASE(77) 
                    
S=18.d0*YD1*ZD1**2
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=8.d0*ZD1**3
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY2Z=S
END FUNCTION DFY2Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=2
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=4.d0*XD1
CASE(30) 
                    
S=2.d0*YD1
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=6.d0*ZD1
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=6.d0*XD1**2
CASE(47) 
                    
S=4.d0*XD1*YD1
CASE(48) 
                    
S=2.d0*YD1**2
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=12.d0*XD1*ZD1
CASE(51) 
                    
S=6.d0*YD1*ZD1
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=12.d0*ZD1**2
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=8.d0*XD1**3
CASE(70) 
                    
S=6.d0*XD1**2*YD1
CASE(71) 
                    
S=4.d0*XD1*YD1**2
CASE(72) 
                    
S=2.d0*YD1**3
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=18.d0*XD1**2*ZD1
CASE(75) 
                    
S=12.d0*XD1*YD1*ZD1
CASE(76) 
                    
S=6.d0*YD1**2*ZD1
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=24.d0*XD1*ZD1**2
CASE(79) 
                    
S=12.d0*YD1*ZD1**2
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=20*ZD1**3
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXZ2=S
END FUNCTION DFXZ2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFYZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1) 
                    
S=0.0D0
CASE (2) 
                    
S=0.0D0
CASE (3) 
                    
S=0.0D0
CASE (4) 
                    
S=0.0D0
CASE (5) 
                    
S=0.0D0
CASE (6) 
                    
S=0.0D0
CASE (7) 
                    
S=0.0D0
CASE (8) 
                    
S=0.0D0
CASE (9) 
                    
S=0.0D0
CASE (10) 
                    
S=0.0D0
CASE (11) 
                    
S=0.0D0
CASE (12) 
                    
S=0.0D0
CASE (13) 
                    
S=0.0D0
CASE (14) 
                    
S=0.0D0
CASE (15) 
                    
S=0.0D0
CASE (16) 
                    
S=0.0D0
CASE (17) 
                    
S=0.0D0
CASE (18) 
                    
S=2
CASE (19) 
                    
S=0.0D0
CASE (20) 
                    
S=0.0D0
CASE (21) 
                    
S=0.0D0
CASE (22) 
                    
S=0.0D0
CASE (23) 
                    
S=0.0D0
CASE (24) 
                    
S=0.0D0
CASE (25) 
                    
S=0.0D0
CASE (26) 
                    
S=0.0D0
CASE (27) 
                    
S=0.0D0
CASE (28) 
                    
S=0.0D0
CASE (29) 
                    
S=0.0D0
CASE (30) 
                    
S=2.d0*XD1
CASE (31) 
                    
S=4.d0*YD1
CASE (32) 
                    
S=0.0D0
CASE (33) 
                    
S=6.d0*ZD1
CASE (34) 
                    
S=0.0D0
CASE (35) 
                    
S=0.0D0
CASE (36) 
                    
S=0.0D0
CASE (37) 
                    
S=0.0D0
CASE (38) 
                    
S=0.0D0
CASE (39) 
                    
S=0.0D0
CASE (40) 
                    
S=0.0D0
CASE (41) 
                    
S=0.0D0
CASE (42) 
                    
S=0.0D0
CASE (43) 
                    
S=0.0D0
CASE (44) 
                    
S=0.0D0
CASE (45) 
                    
S=0.0D0
CASE (46) 
                    
S=0.0D0
CASE (47) 
                    
S=2.d0*XD1**2
CASE (48) 
                    
S=4.d0*XD1*YD1
CASE (49) 
                    
S=6.d0*YD1**2
CASE (50) 
                    
S=0.0D0
CASE (51) 
                    
S=6.d0*XD1*ZD1
CASE (52) 
                    
S=12.d0*YD1*ZD1
CASE (53) 
                    
S=0.0D0
CASE (54) 
                    
S=12.d0*ZD1**2
CASE (55) 
                    
S=0.0D0
CASE (56) 
                    
S=0.0D0
CASE (57) 
                    
S=0.0D0
CASE (58) 
                    
S=0.0D0
CASE (59) 
                    
S=0.0D0
CASE (60) 
                    
S=0.0D0
CASE (61) 
                    
S=0.0D0
CASE (62) 
                    
S=0.0D0
CASE (63) 
                    
S=0.0D0
CASE (64) 
                    
S=0.0D0
CASE (65) 
                    
S=0.0D0
CASE (66) 
                    
S=0.0D0
CASE (67) 
                    
S=0.0D0
CASE (68) 
                    
S=0.0D0
CASE (69) 
                    
S=0.0D0
CASE (70) 
                    
S=2.d0*XD1**3
CASE (71) 
                    
S=4.d0*XD1**2*YD1
CASE (72) 
                    
S=6.d0*XD1*YD1**2
CASE (73) 
                    
S=8.d0*YD1**3
CASE (74) 
                    
S=0.0D0
CASE (75) 
                    
S=6.d0*XD1**2*ZD1
CASE (76) 
                    
S=12.d0*XD1*YD1*ZD1
CASE (77) 
                    
S=18.d0*YD1**2*ZD1
CASE (78) 
                    
S=0.0D0
CASE (79) 
                    
S=12.d0*XD1*ZD1**2
CASE (80) 
                    
S=24.d0*YD1*ZD1**2
CASE (81) 
                    
S=0.0D0
CASE (82) 
                    
S=20*ZD1**3
CASE (83) 
                    
S=0.0D0
END SELECT
DFYZ2=S
END FUNCTION DFYZ2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1)  
                    
S=0.0D0
CASE (2)  
                    
S=0.0D0
CASE (3)  
                    
S=0.0D0
CASE (4)  
                    
S=0.0D0
CASE (5)  
                    
S=0.0D0
CASE (6)  
                    
S=0.0D0
CASE (7)  
                    
S=0.0D0
CASE (8)  
                    
S=0.0D0
CASE (9)  
                    
S=0.0D0
CASE (10)  
                    
S=0.0D0
CASE (11)  
                    
S=0.0D0
CASE (12)  
                    
S=0.0D0
CASE (13)  
                    
S=0.0D0
CASE (14)  
                    
S=0.0D0
CASE (15)  
                    
S=0.0D0
CASE (16)  
                    
S=0.0D0
CASE (17)  
                    
S=0.0D0
CASE (18)  
                    
S=0.0D0
CASE (19)  
                    
S=6
CASE (20)  
                    
S=0.0D0
CASE (21)  
                    
S=0.0D0
CASE (22)  
                    
S=0.0D0
CASE (23)  
                    
S=0.0D0
CASE (24)  
                    
S=0.0D0
CASE (25)  
                    
S=0.0D0
CASE (26)  
                    
S=0.0D0
CASE (27)  
                    
S=0.0D0
CASE (28)  
                    
S=0.0D0
CASE (29)  
                    
S=0.0D0
CASE (30)  
                    
S=0.0D0
CASE (31)  
                    
S=0.0D0
CASE (32)  
                    
S=6.d0*XD1
CASE (33)  
                    
S=6.d0*YD1
CASE (34)  
                    
S=24.d0*ZD1
CASE (35)  
                    
S=0.0D0
CASE (36)  
                    
S=0.0D0
CASE (37)  
                    
S=0.0D0
CASE (38)  
                    
S=0.0D0
CASE (39)  
                    
S=0.0D0
CASE (40)  
                    
S=0.0D0
CASE (41)  
                    
S=0.0D0
CASE (42)  
                    
S=0.0D0
CASE (43)  
                    
S=0.0D0
CASE (44)  
                    
S=0.0D0
CASE (45)  
                    
S=0.0D0
CASE (46)  
                    
S=0.0D0
CASE (47)  
                    
S=0.0D0
CASE (48)  
                    
S=0.0D0
CASE (49)  
                    
S=0.0D0
CASE (50)  
                    
S=6.d0*XD1**2
CASE (51)  
                    
S=6.d0*XD1*YD1
CASE (52)  
                    
S=6.d0*YD1**2
CASE (53)  
                    
S=24.d0*XD1*ZD1
CASE (54)  
                    
S=24.d0*YD1*ZD1
CASE (55)  
                    
S=60*ZD1**2
CASE (56)  
                    
S=0.0D0
CASE (57)  
                    
S=0.0D0
CASE (58)  
                    
S=0.0D0
CASE (59)  
                    
S=0.0D0
CASE (60)  
                    
S=0.0D0
CASE (61)  
                    
S=0.0D0
CASE (62)  
                    
S=0.0D0
CASE (63)  
                    
S=0.0D0
CASE (64)  
                    
S=0.0D0
CASE (65)  
                    
S=0.0D0
CASE (66)  
                    
S=0.0D0
CASE (67)  
                    
S=0.0D0
CASE (68)  
                    
S=0.0D0
CASE (69)  
                    
S=0.0D0
CASE (70)  
                    
S=0.0D0
CASE (71)  
                    
S=0.0D0
CASE (72)  
                    
S=0.0D0
CASE (73)  
                    
S=0.0D0
CASE (74)  
                    
S=6.d0*XD1**3
CASE (75)  
                    
S=6.d0*XD1**2*YD1
CASE (76)  
                    
S=6.d0*XD1*YD1**2
CASE (77)  
                    
S=6.d0*YD1**3
CASE (78)  
                    
S=24.d0*XD1**2*ZD1
CASE (79)  
                    
S=24.d0*XD1*YD1*ZD1
CASE (80)  
                    
S=24.d0*YD1**2*ZD1
CASE (81)  
                    
S=60*XD1*ZD1**2
CASE (82)  
                    
S=60*YD1*ZD1**2
CASE (83)  
                    
S=120*ZD1**3
END SELECT
DFZ3=S
END FUNCTION DFZ3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=24
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=120*XD1
CASE(36) 
                    
S=24.d0*YD1
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=24.d0*ZD1
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=360*XD1**2
CASE(57) 
                    
S=120*XD1*YD1
CASE(58) 
                    
S=24.d0*YD1**2
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=120*XD1*ZD1
CASE(64) 
                    
S=24.d0*YD1*ZD1
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=24.d0*ZD1**2
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX4=S
END FUNCTION DFX4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=6
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=24.d0*XD1
CASE(37) 
                    
S=12.d0*YD1
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=6.d0*ZD1
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=60*XD1**2
CASE(58) 
                    
S=48.d0*XD1*YD1
CASE(59) 
                    
S=18.d0*YD1**2
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=24.d0*XD1*ZD1
CASE(65) 
                    
S=12.d0*YD1*ZD1
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=6.d0*ZD1**2
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3Y=S
END FUNCTION DFX3Y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=4
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=12.d0*XD1
CASE(38) 
                    
S=12.d0*YD1
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=4.d0*ZD1
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=24.d0*XD1**2
CASE(59) 
                    
S=36.d0*XD1*YD1
CASE(60) 
                    
S=24.d0*YD1**2
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=12.d0*XD1*ZD1
CASE(66) 
                    
S=12.d0*YD1*ZD1
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=4.d0*ZD1**2
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Y2=S
END FUNCTION DFX2Y2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=6
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=12.d0*XD1
CASE(39) 
                    
S=24.d0*YD1
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=6.d0*ZD1
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=18.d0*XD1**2
CASE(60) 
                    
S=48.d0*XD1*YD1
CASE(61) 
                    
S=60*YD1**2
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=12.d0*XD1*ZD1
CASE(67) 
                    
S=24.d0*YD1*ZD1
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=6.d0*ZD1**2
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY3=S
END FUNCTION DFXY3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=24
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=24.d0*XD1
CASE(40) 
                    
S=120*YD1
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=24.d0*ZD1
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=24.d0*XD1**2
CASE(61) 
                    
S=120*XD1*YD1
CASE(62) 
                    
S=360*YD1**2
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=24.d0*XD1*ZD1
CASE(68) 
                    
S=120*YD1*ZD1
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=24.d0*ZD1**2
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY4=S
END FUNCTION DFY4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=6
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=24.d0*XD1
CASE(42) 
                    
S=6.d0*YD1
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=12.d0*ZD1
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=60*XD1**2
CASE(64) 
                    
S=24.d0*XD1*YD1
CASE(65) 
                    
S=6.d0*YD1**2
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=48.d0*XD1*ZD1
CASE(70) 
                    
S=12.d0*YD1*ZD1
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=18.d0*ZD1**2
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3Z=S
END FUNCTION DFX3Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=2
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=6.d0*XD1
CASE(43) 
                    
S=4.d0*YD1
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=4.d0*ZD1
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=12.d0*XD1**2
CASE(65) 
                    
S=12.d0*XD1*YD1
CASE(66) 
                    
S=6.d0*YD1**2
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=12.d0*XD1*ZD1
CASE(71) 
                    
S=8.d0*YD1*ZD1
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=6.d0*ZD1**2
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2YZ=S
END FUNCTION DFX2YZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=2
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=4.d0*XD1
CASE(44) 
                    
S=6.d0*YD1
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=4.d0*ZD1
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=6.d0*XD1**2
CASE(66) 
                    
S=12.d0*XD1*YD1
CASE(67) 
                    
S=12.d0*YD1**2
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=8.d0*XD1*ZD1
CASE(72) 
                    
S=12.d0*YD1*ZD1
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=6.d0*ZD1**2
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY2Z=S
END FUNCTION DFXY2Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=6
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=6.d0*XD1
CASE(45) 
                    
S=24.d0*YD1
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=12.d0*ZD1
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=6.d0*XD1**2
CASE(67) 
                    
S=24.d0*XD1*YD1
CASE(68) 
                    
S=60*YD1**2
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=12.d0*XD1*ZD1
CASE(73) 
                    
S=48.d0*YD1*ZD1
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=18.d0*ZD1**2
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY3Z=S
END FUNCTION DFY3Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=4
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=12.d0*XD1
CASE(47) 
                    
S=4.d0*YD1
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=12.d0*ZD1
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=24.d0*XD1**2
CASE(70) 
                    
S=12.d0*XD1*YD1
CASE(71) 
                    
S=4.d0*YD1**2
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=36.d0*XD1*ZD1
CASE(75) 
                    
S=12.d0*YD1*ZD1
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=24.d0*ZD1**2
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Z2=S
END FUNCTION DFX2Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXYZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=2
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=4.d0*XD1
CASE(48) 
                    
S=4.d0*YD1
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=6.d0*ZD1
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=6.d0*XD1**2
CASE(71) 
                    
S=8.d0*XD1*YD1
CASE(72) 
                    
S=6.d0*YD1**2
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=12.d0*XD1*ZD1
CASE(76) 
                    
S=12.d0*YD1*ZD1
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=12.d0*ZD1**2
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXYZ2=S
END FUNCTION DFXYZ2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=4
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=4.d0*XD1
CASE(49) 
                    
S=12.d0*YD1
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=12.d0*ZD1
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=4.d0*XD1**2
CASE(72) 
                    
S=12.d0*XD1*YD1
CASE(73) 
                    
S=24.d0*YD1**2
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=12.d0*XD1*ZD1
CASE(77) 
                    
S=36.d0*YD1*ZD1
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=24.d0*ZD1**2
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY2Z2=S
END FUNCTION DFY2Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=6
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=12.d0*XD1
CASE(51) 
                    
S=6.d0*YD1
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=24.d0*ZD1
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=18.d0*XD1**2
CASE(75) 
                    
S=12.d0*XD1*YD1
CASE(76) 
                    
S=6.d0*YD1**2
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=48.d0*XD1*ZD1
CASE(79) 
                    
S=24.d0*YD1*ZD1
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=60*ZD1**2
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXZ3=S
END FUNCTION DFXZ3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFYZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1)  
                    
S=0.0D0
CASE (2)  
                    
S=0.0D0
CASE (3)  
                    
S=0.0D0
CASE (4)  
                    
S=0.0D0
CASE (5)  
                    
S=0.0D0
CASE (6)  
                    
S=0.0D0
CASE (7)  
                    
S=0.0D0
CASE (8)  
                    
S=0.0D0
CASE (9)  
                    
S=0.0D0
CASE (10)  
                    
S=0.0D0
CASE (11)  
                    
S=0.0D0
CASE (12)  
                    
S=0.0D0
CASE (13)  
                    
S=0.0D0
CASE (14)  
                    
S=0.0D0
CASE (15)  
                    
S=0.0D0
CASE (16)  
                    
S=0.0D0
CASE (17)  
                    
S=0.0D0
CASE (18)  
                    
S=0.0D0
CASE (19)  
                    
S=0.0D0
CASE (20)  
                    
S=0.0D0
CASE (21)  
                    
S=0.0D0
CASE (22)  
                    
S=0.0D0
CASE (23)  
                    
S=0.0D0
CASE (24)  
                    
S=0.0D0
CASE (25)  
                    
S=0.0D0
CASE (26)  
                    
S=0.0D0
CASE (27)  
                    
S=0.0D0
CASE (28)  
                    
S=0.0D0
CASE (29)  
                    
S=0.0D0
CASE (30)  
                    
S=0.0D0
CASE (31)  
                    
S=0.0D0
CASE (32)  
                    
S=0.0D0
CASE (33)  
                    
S=6
CASE (34)  
                    
S=0.0D0
CASE (35)  
                    
S=0.0D0
CASE (36)  
                    
S=0.0D0
CASE (37)  
                    
S=0.0D0
CASE (38)  
                    
S=0.0D0
CASE (39)  
                    
S=0.0D0
CASE (40)  
                    
S=0.0D0
CASE (41)  
                    
S=0.0D0
CASE (42)  
                    
S=0.0D0
CASE (43)  
                    
S=0.0D0
CASE (44)  
                    
S=0.0D0
CASE (45)  
                    
S=0.0D0
CASE (46)  
                    
S=0.0D0
CASE (47)  
                    
S=0.0D0
CASE (48)  
                    
S=0.0D0
CASE (49)  
                    
S=0.0D0
CASE (50)  
                    
S=0.0D0
CASE (51)  
                    
S=6.d0*XD1
CASE (52)  
                    
S=12.d0*YD1
CASE (53)  
                    
S=0.0D0
CASE (54)  
                    
S=24.d0*ZD1
CASE (55)  
                    
S=0.0D0
CASE (56)  
                    
S=0.0D0
CASE (57)  
                    
S=0.0D0
CASE (58)  
                    
S=0.0D0
CASE (59)  
                    
S=0.0D0
CASE (60)  
                    
S=0.0D0
CASE (61)  
                    
S=0.0D0
CASE (62)  
                    
S=0.0D0
CASE (63)  
                    
S=0.0D0
CASE (64)  
                    
S=0.0D0
CASE (65)  
                    
S=0.0D0
CASE (66)  
                    
S=0.0D0
CASE (67)  
                    
S=0.0D0
CASE (68)  
                    
S=0.0D0
CASE (69)  
                    
S=0.0D0
CASE (70)  
                    
S=0.0D0
CASE (71)  
                    
S=0.0D0
CASE (72)  
                    
S=0.0D0
CASE (73)  
                    
S=0.0D0
CASE (74)  
                    
S=0.0D0
CASE (75)  
                    
S=6.d0*XD1**2
CASE (76)  
                    
S=12.d0*XD1*YD1
CASE (77)  
                    
S=18.d0*YD1**2
CASE (78)  
                    
S=0.0D0
CASE (79)  
                    
S=24.d0*XD1*ZD1
CASE (80)  
                    
S=48.d0*YD1*ZD1
CASE (81)  
                    
S=0.0D0
CASE (82)  
                    
S=60*ZD1**2
CASE (83)  
                    
S=0.0D0
END SELECT
DFYZ3=S
END FUNCTION DFYZ3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1)  
                    
S=0.0D0
CASE (2)  
                    
S=0.0D0
CASE (3)  
                    
S=0.0D0
CASE (4)  
                    
S=0.0D0
CASE (5)  
                    
S=0.0D0
CASE (6)  
                    
S=0.0D0
CASE (7)  
                    
S=0.0D0
CASE (8)  
                    
S=0.0D0
CASE (9)  
                    
S=0.0D0
CASE (10)  
                    
S=0.0D0
CASE (11)  
                    
S=0.0D0
CASE (12)  
                    
S=0.0D0
CASE (13)  
                    
S=0.0D0
CASE (14)  
                    
S=0.0D0
CASE (15)  
                    
S=0.0D0
CASE (16)  
                    
S=0.0D0
CASE (17)  
                    
S=0.0D0
CASE (18)  
                    
S=0.0D0
CASE (19)  
                    
S=0.0D0
CASE (20)  
                    
S=0.0D0
CASE (21)  
                    
S=0.0D0
CASE (22)  
                    
S=0.0D0
CASE (23)  
                    
S=0.0D0
CASE (24)  
                    
S=0.0D0
CASE (25)  
                    
S=0.0D0
CASE (26)  
                    
S=0.0D0
CASE (27)  
                    
S=0.0D0
CASE (28)  
                    
S=0.0D0
CASE (29)  
                    
S=0.0D0
CASE (30)  
                    
S=0.0D0
CASE (31)  
                    
S=0.0D0
CASE (32)  
                    
S=0.0D0
CASE (33)  
                    
S=0.0D0
CASE (34)  
                    
S=24
CASE (35)  
                    
S=0.0D0
CASE (36)  
                    
S=0.0D0
CASE (37)  
                    
S=0.0D0
CASE (38)  
                    
S=0.0D0
CASE (39)  
                    
S=0.0D0
CASE (40)  
                    
S=0.0D0
CASE (41)  
                    
S=0.0D0
CASE (42)  
                    
S=0.0D0
CASE (43)  
                    
S=0.0D0
CASE (44)  
                    
S=0.0D0
CASE (45)  
                    
S=0.0D0
CASE (46)  
                    
S=0.0D0
CASE (47)  
                    
S=0.0D0
CASE (48)  
                    
S=0.0D0
CASE (49)  
                    
S=0.0D0
CASE (50)  
                    
S=0.0D0
CASE (51)  
                    
S=0.0D0
CASE (52)  
                    
S=0.0D0
CASE (53)  
                    
S=24.d0*XD1
CASE (54)  
                    
S=24.d0*YD1
CASE (55)  
                    
S=120*ZD1
CASE (56)  
                    
S=0.0D0
CASE (57)  
                    
S=0.0D0
CASE (58)  
                    
S=0.0D0
CASE (59)  
                    
S=0.0D0
CASE (60)  
                    
S=0.0D0
CASE (61)  
                    
S=0.0D0
CASE (62)  
                    
S=0.0D0
CASE (63)  
                    
S=0.0D0
CASE (64)  
                    
S=0.0D0
CASE (65)  
                    
S=0.0D0
CASE (66)  
                    
S=0.0D0
CASE (67)  
                    
S=0.0D0
CASE (68)  
                    
S=0.0D0
CASE (69)  
                    
S=0.0D0
CASE (70)  
                    
S=0.0D0
CASE (71)  
                    
S=0.0D0
CASE (72)  
                    
S=0.0D0
CASE (73)  
                    
S=0.0D0
CASE (74)  
                    
S=0.0D0
CASE (75)  
                    
S=0.0D0
CASE (76)  
                    
S=0.0D0
CASE (77)  
                    
S=0.0D0
CASE (78)  
                    
S=24.d0*XD1**2
CASE (79)  
                    
S=24.d0*XD1*YD1
CASE (80)  
                    
S=24.d0*YD1**2
CASE (81)  
                    
S=120*XD1*ZD1
CASE (82)  
                    
S=120*YD1*ZD1
CASE (83)  
                    
S=360*ZD1**2
END SELECT
DFZ4=S
END FUNCTION DFZ4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=120
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=720*XD1
CASE(57) 
                    
S=120*YD1
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=120*ZD1
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX5=S
END FUNCTION DFX5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX4Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=24
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=120*XD1
CASE(58) 
                    
S=48.d0*YD1
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=24.d0*ZD1
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX4Y=S
END FUNCTION DFX4Y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=12
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=48.d0*XD1
CASE(59) 
                    
S=36.d0*YD1
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=12.d0*ZD1
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3Y2=S
END FUNCTION DFX3Y2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Y3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=12
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=36.d0*XD1
CASE(60) 
                    
S=48.d0*YD1
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=12.d0*ZD1
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Y3=S
END FUNCTION DFX2Y3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=24
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=48.d0*XD1
CASE(61) 
                    
S=120*YD1
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=24.d0*ZD1
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY4=S
END FUNCTION DFXY4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=120
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=120*XD1
CASE(62) 
                    
S=720*YD1
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=120*ZD1
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY5=S
END FUNCTION DFY5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=24
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=120*XD1
CASE(64) 
                    
S=24.d0*YD1
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=48.d0*ZD1
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX4Z=S
END FUNCTION DFX4Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=6
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=24.d0*XD1
CASE(65) 
                    
S=12.d0*YD1
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=12.d0*ZD1
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3YZ=S
END FUNCTION DFX3YZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Y2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=4
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=12.d0*XD1
CASE(66) 
                    
S=12.d0*YD1
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=8.d0*ZD1
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Y2Z=S
END FUNCTION DFX2Y2Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=6
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=12.d0*XD1
CASE(67) 
                    
S=24.d0*YD1
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=12.d0*ZD1
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY3Z=S
END FUNCTION DFXY3Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=24
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=24.d0*XD1
CASE(68) 
                    
S=120*YD1
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=48.d0*ZD1
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY4Z=S
END FUNCTION DFY4Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=12
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=48.d0*XD1
CASE(70) 
                    
S=12.d0*YD1
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=36.d0*ZD1
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3Z2=S
END FUNCTION DFX3Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2YZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=4
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=12.d0*XD1
CASE(71) 
                    
S=8.d0*YD1
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=12.d0*ZD1
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2YZ2=S
END FUNCTION DFX2YZ2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=4
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=8.d0*XD1
CASE(72) 
                    
S=12.d0*YD1
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=12.d0*ZD1
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY2Z2=S
END FUNCTION DFXY2Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=12
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=12.d0*XD1
CASE(73) 
                    
S=48.d0*YD1
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=36.d0*ZD1
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY3Z2=S
END FUNCTION DFY3Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=12
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=36.d0*XD1
CASE(75) 
                    
S=12.d0*YD1
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=48.d0*ZD1
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Z3=S
END FUNCTION DFX2Z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXYZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=6
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=12.d0*XD1
CASE(76) 
                    
S=12.d0*YD1
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=24.d0*ZD1
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXYZ3=S
END FUNCTION DFXYZ3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=12
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=12.d0*XD1
CASE(77) 
                    
S=36.d0*YD1
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=48.d0*ZD1
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY2Z3=S
END FUNCTION DFY2Z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=24
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=48.d0*XD1
CASE(79) 
                    
S=24.d0*YD1
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=120*ZD1
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXZ4=S
END FUNCTION DFXZ4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFYZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1)  
                    
S=0.0D0
CASE (2)  
                    
S=0.0D0
CASE (3)  
                    
S=0.0D0
CASE (4)  
                    
S=0.0D0
CASE (5)  
                    
S=0.0D0
CASE (6)  
                    
S=0.0D0
CASE (7)  
                    
S=0.0D0
CASE (8)  
                    
S=0.0D0
CASE (9)  
                    
S=0.0D0
CASE (10)  
                    
S=0.0D0
CASE (11)  
                    
S=0.0D0
CASE (12)  
                    
S=0.0D0
CASE (13)  
                    
S=0.0D0
CASE (14)  
                    
S=0.0D0
CASE (15)  
                    
S=0.0D0
CASE (16)  
                    
S=0.0D0
CASE (17)  
                    
S=0.0D0
CASE (18)  
                    
S=0.0D0
CASE (19)  
                    
S=0.0D0
CASE (20)  
                    
S=0.0D0
CASE (21)  
                    
S=0.0D0
CASE (22)  
                    
S=0.0D0
CASE (23)  
                    
S=0.0D0
CASE (24)  
                    
S=0.0D0
CASE (25)  
                    
S=0.0D0
CASE (26)  
                    
S=0.0D0
CASE (27)  
                    
S=0.0D0
CASE (28)  
                    
S=0.0D0
CASE (29)  
                    
S=0.0D0
CASE (30)  
                    
S=0.0D0
CASE (31)  
                    
S=0.0D0
CASE (32)  
                    
S=0.0D0
CASE (33)  
                    
S=0.0D0
CASE (34)  
                    
S=0.0D0
CASE (35)  
                    
S=0.0D0
CASE (36)  
                    
S=0.0D0
CASE (37)  
                    
S=0.0D0
CASE (38)  
                    
S=0.0D0
CASE (39)  
                    
S=0.0D0
CASE (40)  
                    
S=0.0D0
CASE (41)  
                    
S=0.0D0
CASE (42)  
                    
S=0.0D0
CASE (43)  
                    
S=0.0D0
CASE (44)  
                    
S=0.0D0
CASE (45)  
                    
S=0.0D0
CASE (46)  
                    
S=0.0D0
CASE (47)  
                    
S=0.0D0
CASE (48)  
                    
S=0.0D0
CASE (49)  
                    
S=0.0D0
CASE (50)  
                    
S=0.0D0
CASE (51)  
                    
S=0.0D0
CASE (52)  
                    
S=0.0D0
CASE (53)  
                    
S=0.0D0
CASE (54)  
                    
S=24
CASE (55)  
                    
S=0.0D0
CASE (56)  
                    
S=0.0D0
CASE (57)  
                    
S=0.0D0
CASE (58)  
                    
S=0.0D0
CASE (59)  
                    
S=0.0D0
CASE (60)  
                    
S=0.0D0
CASE (61)  
                    
S=0.0D0
CASE (62)  
                    
S=0.0D0
CASE (63)  
                    
S=0.0D0
CASE (64)  
                    
S=0.0D0
CASE (65)  
                    
S=0.0D0
CASE (66)  
                    
S=0.0D0
CASE (67)  
                    
S=0.0D0
CASE (68)  
                    
S=0.0D0
CASE (69)  
                    
S=0.0D0
CASE (70)  
                    
S=0.0D0
CASE (71)  
                    
S=0.0D0
CASE (72)  
                    
S=0.0D0
CASE (73)  
                    
S=0.0D0
CASE (74)  
                    
S=0.0D0
CASE (75)  
                    
S=0.0D0
CASE (76)  
                    
S=0.0D0
CASE (77)  
                    
S=0.0D0
CASE (78)  
                    
S=0.0D0
CASE (79)  
                    
S=24.d0*XD1
CASE (80)  
                    
S=48.d0*YD1
CASE (81)  
                    
S=0.0D0
CASE (82)  
                    
S=120*ZD1
CASE (83)  
                    
S=0.0D0
END SELECT
DFYZ4=S
END FUNCTION DFYZ4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFZ5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1)  
                    
S=0.0D0
CASE(2)  
                    
S=0.0D0
CASE(3)  
                    
S=0.0D0
CASE(4)  
                    
S=0.0D0
CASE(5)  
                    
S=0.0D0
CASE(6)  
                    
S=0.0D0
CASE(7)  
                    
S=0.0D0
CASE(8)  
                    
S=0.0D0
CASE(9)  
                    
S=0.0D0
CASE(10)  
                    
S=0.0D0
CASE(11)  
                    
S=0.0D0
CASE(12)  
                    
S=0.0D0
CASE(13)  
                    
S=0.0D0
CASE(14)  
                    
S=0.0D0
CASE(15)  
                    
S=0.0D0
CASE(16)  
                    
S=0.0D0
CASE(17)  
                    
S=0.0D0
CASE(18)  
                    
S=0.0D0
CASE(19)  
                    
S=0.0D0
CASE(20)  
                    
S=0.0D0
CASE(21)  
                    
S=0.0D0
CASE(22)  
                    
S=0.0D0
CASE(23)  
                    
S=0.0D0
CASE(24)  
                    
S=0.0D0
CASE(25)  
                    
S=0.0D0
CASE(26)  
                    
S=0.0D0
CASE(27)  
                    
S=0.0D0
CASE(28)  
                    
S=0.0D0
CASE(29)  
                    
S=0.0D0
CASE(30)  
                    
S=0.0D0
CASE(31)  
                    
S=0.0D0
CASE(32)  
                    
S=0.0D0
CASE(33)  
                    
S=0.0D0
CASE(34)  
                    
S=0.0D0
CASE(35)  
                    
S=0.0D0
CASE(36)  
                    
S=0.0D0
CASE(37)  
                    
S=0.0D0
CASE(38)  
                    
S=0.0D0
CASE(39)  
                    
S=0.0D0
CASE(40)  
                    
S=0.0D0
CASE(41)  
                    
S=0.0D0
CASE(42)  
                    
S=0.0D0
CASE(43)  
                    
S=0.0D0
CASE(44)  
                    
S=0.0D0
CASE(45)  
                    
S=0.0D0
CASE(46)  
                    
S=0.0D0
CASE(47)  
                    
S=0.0D0
CASE(48)  
                    
S=0.0D0
CASE(49)  
                    
S=0.0D0
CASE(50)  
                    
S=0.0D0
CASE(51)  
                    
S=0.0D0
CASE(52)  
                    
S=0.0D0
CASE(53)  
                    
S=0.0D0
CASE(54)  
                    
S=0.0D0
CASE(55)  
                    
S=120
CASE(56)  
                    
S=0.0D0
CASE(57)  
                    
S=0.0D0
CASE(58)  
                    
S=0.0D0
CASE(59)  
                    
S=0.0D0
CASE(60)  
                    
S=0.0D0
CASE(61)  
                    
S=0.0D0
CASE(62)  
                    
S=0.0D0
CASE(63)  
                    
S=0.0D0
CASE(64)  
                    
S=0.0D0
CASE(65)  
                    
S=0.0D0
CASE(66)  
                    
S=0.0D0
CASE(67)  
                    
S=0.0D0
CASE(68)  
                    
S=0.0D0
CASE(69)  
                    
S=0.0D0
CASE(70)  
                    
S=0.0D0
CASE(71)  
                    
S=0.0D0
CASE(72)  
                    
S=0.0D0
CASE(73)  
                    
S=0.0D0
CASE(74)  
                    
S=0.0D0
CASE(75)  
                    
S=0.0D0
CASE(76)  
                    
S=0.0D0
CASE(77)  
                    
S=0.0D0
CASE(78)  
                    
S=0.0D0
CASE(79)  
                    
S=0.0D0
CASE(80)  
                    
S=0.0D0
CASE(81)  
                    
S=120*XD1
CASE(82)  
                    
S=120*YD1
CASE(83)  
                    
S=720*ZD1
END SELECT
DFZ5=S
END FUNCTION DFZ5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=720
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX6=S
END FUNCTION DFX6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX5Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=120
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX5Y=S
END FUNCTION DFX5Y
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX4Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=48
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX4Y2=S
END FUNCTION DFX4Y2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3Y3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=36
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3Y3=S
END FUNCTION DFX3Y3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Y4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=48
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Y4=S
END FUNCTION DFX2Y4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1) 
                    
S=0.0D0
CASE (2) 
                    
S=0.0D0
CASE (3) 
                    
S=0.0D0
CASE (4) 
                    
S=0.0D0
CASE (5) 
                    
S=0.0D0
CASE (6) 
                    
S=0.0D0
CASE (7) 
                    
S=0.0D0
CASE (8) 
                    
S=0.0D0
CASE (9) 
                    
S=0.0D0
CASE (10) 
                    
S=0.0D0
CASE (11) 
                    
S=0.0D0
CASE (12) 
                    
S=0.0D0
CASE (13) 
                    
S=0.0D0
CASE (14) 
                    
S=0.0D0
CASE (15) 
                    
S=0.0D0
CASE (16) 
                    
S=0.0D0
CASE (17) 
                    
S=0.0D0
CASE (18) 
                    
S=0.0D0
CASE (19) 
                    
S=0.0D0
CASE (20) 
                    
S=0.0D0
CASE (21) 
                    
S=0.0D0
CASE (22) 
                    
S=0.0D0
CASE (23) 
                    
S=0.0D0
CASE (24) 
                    
S=0.0D0
CASE (25) 
                    
S=0.0D0
CASE (26) 
                    
S=0.0D0
CASE (27) 
                    
S=0.0D0
CASE (28) 
                    
S=0.0D0
CASE (29) 
                    
S=0.0D0
CASE (30) 
                    
S=0.0D0
CASE (31) 
                    
S=0.0D0
CASE (32) 
                    
S=0.0D0
CASE (33) 
                    
S=0.0D0
CASE (34) 
                    
S=0.0D0
CASE (35) 
                    
S=0.0D0
CASE (36) 
                    
S=0.0D0
CASE (37) 
                    
S=0.0D0
CASE (38) 
                    
S=0.0D0
CASE (39) 
                    
S=0.0D0
CASE (40) 
                    
S=0.0D0
CASE (41) 
                    
S=0.0D0
CASE (42) 
                    
S=0.0D0
CASE (43) 
                    
S=0.0D0
CASE (44) 
                    
S=0.0D0
CASE (45) 
                    
S=0.0D0
CASE (46) 
                    
S=0.0D0
CASE (47) 
                    
S=0.0D0
CASE (48) 
                    
S=0.0D0
CASE (49) 
                    
S=0.0D0
CASE (50) 
                    
S=0.0D0
CASE (51) 
                    
S=0.0D0
CASE (52) 
                    
S=0.0D0
CASE (53) 
                    
S=0.0D0
CASE (54) 
                    
S=0.0D0
CASE (55) 
                    
S=0.0D0
CASE (56) 
                    
S=0.0D0
CASE (57) 
                    
S=0.0D0
CASE (58) 
                    
S=0.0D0
CASE (59) 
                    
S=0.0D0
CASE (60) 
                    
S=0.0D0
CASE (61) 
                    
S=120
CASE (62) 
                    
S=0.0D0
CASE (63) 
                    
S=0.0D0
CASE (64) 
                    
S=0.0D0
CASE (65) 
                    
S=0.0D0
CASE (66) 
                    
S=0.0D0
CASE (67) 
                    
S=0.0D0
CASE (68) 
                    
S=0.0D0
CASE (69) 
                    
S=0.0D0
CASE (70) 
                    
S=0.0D0
CASE (71) 
                    
S=0.0D0
CASE (72) 
                    
S=0.0D0
CASE (73) 
                    
S=0.0D0
CASE (74) 
                    
S=0.0D0
CASE (75) 
                    
S=0.0D0
CASE (76) 
                    
S=0.0D0
CASE (77) 
                    
S=0.0D0
CASE (78) 
                    
S=0.0D0
CASE (79) 
                    
S=0.0D0
CASE (80) 
                    
S=0.0D0
CASE (81) 
                    
S=0.0D0
CASE (82) 
                    
S=0.0D0
CASE (83) 
                    
S=0.0D0
END SELECT
DFXY5=S
END FUNCTION DFXY5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=720
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY6=S
END FUNCTION DFY6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX5Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=120
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX5Z=S
END FUNCTION DFX5Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX4YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=24
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX4YZ=S
END FUNCTION DFX4YZ
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3Y2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=12
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3Y2Z=S
END FUNCTION DFX3Y2Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Y3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1)
                    
S=0.0D0
CASE(2)
                    
S=0.0D0
CASE(3)
                    
S=0.0D0
CASE(4)
                    
S=0.0D0
CASE(5)
                    
S=0.0D0
CASE(6)
                    
S=0.0D0
CASE(7)
                    
S=0.0D0
CASE(8)
                    
S=0.0D0
CASE(9)
                    
S=0.0D0
CASE(10)
                    
S=0.0D0
CASE(11)
                    
S=0.0D0
CASE(12)
                    
S=0.0D0
CASE(13)
                    
S=0.0D0
CASE(14)
                    
S=0.0D0
CASE(15)
                    
S=0.0D0
CASE(16)
                    
S=0.0D0
CASE(17)
                    
S=0.0D0
CASE(18)
                    
S=0.0D0
CASE(19)
                    
S=0.0D0
CASE(20)
                    
S=0.0D0
CASE(21)
                    
S=0.0D0
CASE(22)
                    
S=0.0D0
CASE(23)
                    
S=0.0D0
CASE(24)
                    
S=0.0D0
CASE(25)
                    
S=0.0D0
CASE(26)
                    
S=0.0D0
CASE(27)
                    
S=0.0D0
CASE(28)
                    
S=0.0D0
CASE(29)
                    
S=0.0D0
CASE(30)
                    
S=0.0D0
CASE(31)
                    
S=0.0D0
CASE(32)
                    
S=0.0D0
CASE(33)
                    
S=0.0D0
CASE(34)
                    
S=0.0D0
CASE(35)
                    
S=0.0D0
CASE(36)
                    
S=0.0D0
CASE(37)
                    
S=0.0D0
CASE(38)
                    
S=0.0D0
CASE(39)
                    
S=0.0D0
CASE(40)
                    
S=0.0D0
CASE(41)
                    
S=0.0D0
CASE(42)
                    
S=0.0D0
CASE(43)
                    
S=0.0D0
CASE(44)
                    
S=0.0D0
CASE(45)
                    
S=0.0D0
CASE(46)
                    
S=0.0D0
CASE(47)
                    
S=0.0D0
CASE(48)
                    
S=0.0D0
CASE(49)
                    
S=0.0D0
CASE(50)
                    
S=0.0D0
CASE(51)
                    
S=0.0D0
CASE(52)
                    
S=0.0D0
CASE(53)
                    
S=0.0D0
CASE(54)
                    
S=0.0D0
CASE(55)
                    
S=0.0D0
CASE(56)
                    
S=0.0D0
CASE(57)
                    
S=0.0D0
CASE(58)
                    
S=0.0D0
CASE(59)
                    
S=0.0D0
CASE(60)
                    
S=0.0D0
CASE(61)
                    
S=0.0D0
CASE(62)
                    
S=0.0D0
CASE(63)
                    
S=0.0D0
CASE(64)
                    
S=0.0D0
CASE(65)
                    
S=0.0D0
CASE(66)
                    
S=12
CASE(67)
                    
S=0.0D0
CASE(68)
                    
S=0.0D0
CASE(69)
                    
S=0.0D0
CASE(70)
                    
S=0.0D0
CASE(71)
                    
S=0.0D0
CASE(72)
                    
S=0.0D0
CASE(73)
                    
S=0.0D0
CASE(74)
                    
S=0.0D0
CASE(75)
                    
S=0.0D0
CASE(76)
                    
S=0.0D0
CASE(77)
                    
S=0.0D0
CASE(78)
                    
S=0.0D0
CASE(79)
                    
S=0.0D0
CASE(80)
                    
S=0.0D0
CASE(81)
                    
S=0.0D0
CASE(82)
                    
S=0.0D0
CASE(83)
                    
S=0.0D0
END SELECT
DFX2Y3Z=S
END FUNCTION DFX2Y3Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=24
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY4Z=S
END FUNCTION DFXY4Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY5Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=24
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY5Z=S
END FUNCTION DFY5Z
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX4Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=48
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX4Z2=S
END FUNCTION DFX4Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3YZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=12
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3YZ2=S
END FUNCTION DFX3YZ2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Y2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=8
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Y2Z2=S
END FUNCTION DFX2Y2Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1) 
                    
S=0.0D0
CASE (2) 
                    
S=0.0D0
CASE (3) 
                    
S=0.0D0
CASE (4) 
                    
S=0.0D0
CASE (5) 
                    
S=0.0D0
CASE (6) 
                    
S=0.0D0
CASE (7) 
                    
S=0.0D0
CASE (8) 
                    
S=0.0D0
CASE (9) 
                    
S=0.0D0
CASE (10) 
                    
S=0.0D0
CASE (11) 
                    
S=0.0D0
CASE (12) 
                    
S=0.0D0
CASE (13) 
                    
S=0.0D0
CASE (14) 
                    
S=0.0D0
CASE (15) 
                    
S=0.0D0
CASE (16) 
                    
S=0.0D0
CASE (17) 
                    
S=0.0D0
CASE (18) 
                    
S=0.0D0
CASE (19) 
                    
S=0.0D0
CASE (20) 
                    
S=0.0D0
CASE (21) 
                    
S=0.0D0
CASE (22) 
                    
S=0.0D0
CASE (23) 
                    
S=0.0D0
CASE (24) 
                    
S=0.0D0
CASE (25) 
                    
S=0.0D0
CASE (26) 
                    
S=0.0D0
CASE (27) 
                    
S=0.0D0
CASE (28) 
                    
S=0.0D0
CASE (29) 
                    
S=0.0D0
CASE (30) 
                    
S=0.0D0
CASE (31) 
                    
S=0.0D0
CASE (32) 
                    
S=0.0D0
CASE (33) 
                    
S=0.0D0
CASE (34) 
                    
S=0.0D0
CASE (35) 
                    
S=0.0D0
CASE (36) 
                    
S=0.0D0
CASE (37) 
                    
S=0.0D0
CASE (38) 
                    
S=0.0D0
CASE (39) 
                    
S=0.0D0
CASE (40) 
                    
S=0.0D0
CASE (41) 
                    
S=0.0D0
CASE (42) 
                    
S=0.0D0
CASE (43) 
                    
S=0.0D0
CASE (44) 
                    
S=0.0D0
CASE (45) 
                    
S=0.0D0
CASE (46) 
                    
S=0.0D0
CASE (47) 
                    
S=0.0D0
CASE (48) 
                    
S=0.0D0
CASE (49) 
                    
S=0.0D0
CASE (50) 
                    
S=0.0D0
CASE (51) 
                    
S=0.0D0
CASE (52) 
                    
S=0.0D0
CASE (53) 
                    
S=0.0D0
CASE (54) 
                    
S=0.0D0
CASE (55) 
                    
S=0.0D0
CASE (56) 
                    
S=0.0D0
CASE (57) 
                    
S=0.0D0
CASE (58) 
                    
S=0.0D0
CASE (59) 
                    
S=0.0D0
CASE (60) 
                    
S=0.0D0
CASE (61) 
                    
S=0.0D0
CASE (62) 
                    
S=0.0D0
CASE (63) 
                    
S=0.0D0
CASE (64) 
                    
S=0.0D0
CASE (65) 
                    
S=0.0D0
CASE (66) 
                    
S=0.0D0
CASE (67) 
                    
S=0.0D0
CASE (68) 
                    
S=0.0D0
CASE (69) 
                    
S=0.0D0
CASE (70) 
                    
S=0.0D0
CASE (71) 
                    
S=0.0D0
CASE (72) 
                    
S=12
CASE (73) 
                    
S=0.0D0
CASE (74) 
                    
S=0.0D0
CASE (75) 
                    
S=0.0D0
CASE (76) 
                    
S=0.0D0
CASE (77) 
                    
S=0.0D0
CASE (78) 
                    
S=0.0D0
CASE (79) 
                    
S=0.0D0
CASE (80) 
                    
S=0.0D0
CASE (81) 
                    
S=0.0D0
CASE (82) 
                    
S=0.0D0
CASE (83) 
                    
S=0.0D0
END SELECT
DFXY3Z2=S
END FUNCTION DFXY3Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY4Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=48
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY4Z2=S
END FUNCTION DFY4Z2
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX3Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=36
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX3Z3=S
END FUNCTION DFX3Z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2YZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=12
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2YZ3=S
END FUNCTION DFX2YZ3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXY2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=12
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXY2Z3=S
END FUNCTION DFXY2Z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY3Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=36
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY3Z3=S
END FUNCTION DFY3Z3
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFX2Z4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=48
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFX2Z4=S
END FUNCTION DFX2Z4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXYZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=24
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXYZ4=S
END FUNCTION DFXYZ4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFY2Z4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=48
CASE(81) 
                    
S=0.0D0
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFY2Z4=S
END FUNCTION DFY2Z4
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFXZ5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1) 
                    
S=0.0D0
CASE(2) 
                    
S=0.0D0
CASE(3) 
                    
S=0.0D0
CASE(4) 
                    
S=0.0D0
CASE(5) 
                    
S=0.0D0
CASE(6) 
                    
S=0.0D0
CASE(7) 
                    
S=0.0D0
CASE(8) 
                    
S=0.0D0
CASE(9) 
                    
S=0.0D0
CASE(10) 
                    
S=0.0D0
CASE(11) 
                    
S=0.0D0
CASE(12) 
                    
S=0.0D0
CASE(13) 
                    
S=0.0D0
CASE(14) 
                    
S=0.0D0
CASE(15) 
                    
S=0.0D0
CASE(16) 
                    
S=0.0D0
CASE(17) 
                    
S=0.0D0
CASE(18) 
                    
S=0.0D0
CASE(19) 
                    
S=0.0D0
CASE(20) 
                    
S=0.0D0
CASE(21) 
                    
S=0.0D0
CASE(22) 
                    
S=0.0D0
CASE(23) 
                    
S=0.0D0
CASE(24) 
                    
S=0.0D0
CASE(25) 
                    
S=0.0D0
CASE(26) 
                    
S=0.0D0
CASE(27) 
                    
S=0.0D0
CASE(28) 
                    
S=0.0D0
CASE(29) 
                    
S=0.0D0
CASE(30) 
                    
S=0.0D0
CASE(31) 
                    
S=0.0D0
CASE(32) 
                    
S=0.0D0
CASE(33) 
                    
S=0.0D0
CASE(34) 
                    
S=0.0D0
CASE(35) 
                    
S=0.0D0
CASE(36) 
                    
S=0.0D0
CASE(37) 
                    
S=0.0D0
CASE(38) 
                    
S=0.0D0
CASE(39) 
                    
S=0.0D0
CASE(40) 
                    
S=0.0D0
CASE(41) 
                    
S=0.0D0
CASE(42) 
                    
S=0.0D0
CASE(43) 
                    
S=0.0D0
CASE(44) 
                    
S=0.0D0
CASE(45) 
                    
S=0.0D0
CASE(46) 
                    
S=0.0D0
CASE(47) 
                    
S=0.0D0
CASE(48) 
                    
S=0.0D0
CASE(49) 
                    
S=0.0D0
CASE(50) 
                    
S=0.0D0
CASE(51) 
                    
S=0.0D0
CASE(52) 
                    
S=0.0D0
CASE(53) 
                    
S=0.0D0
CASE(54) 
                    
S=0.0D0
CASE(55) 
                    
S=0.0D0
CASE(56) 
                    
S=0.0D0
CASE(57) 
                    
S=0.0D0
CASE(58) 
                    
S=0.0D0
CASE(59) 
                    
S=0.0D0
CASE(60) 
                    
S=0.0D0
CASE(61) 
                    
S=0.0D0
CASE(62) 
                    
S=0.0D0
CASE(63) 
                    
S=0.0D0
CASE(64) 
                    
S=0.0D0
CASE(65) 
                    
S=0.0D0
CASE(66) 
                    
S=0.0D0
CASE(67) 
                    
S=0.0D0
CASE(68) 
                    
S=0.0D0
CASE(69) 
                    
S=0.0D0
CASE(70) 
                    
S=0.0D0
CASE(71) 
                    
S=0.0D0
CASE(72) 
                    
S=0.0D0
CASE(73) 
                    
S=0.0D0
CASE(74) 
                    
S=0.0D0
CASE(75) 
                    
S=0.0D0
CASE(76) 
                    
S=0.0D0
CASE(77) 
                    
S=0.0D0
CASE(78) 
                    
S=0.0D0
CASE(79) 
                    
S=0.0D0
CASE(80) 
                    
S=0.0D0
CASE(81) 
                    
S=120
CASE(82) 
                    
S=0.0D0
CASE(83) 
                    
S=0.0D0
END SELECT
DFXZ5=S
END FUNCTION DFXZ5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFYZ5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE (1) 
                    
S=0.0D0
CASE (2) 
                    
S=0.0D0
CASE (3) 
                    
S=0.0D0
CASE (4) 
                    
S=0.0D0
CASE (5) 
                    
S=0.0D0
CASE (6) 
                    
S=0.0D0
CASE (7) 
                    
S=0.0D0
CASE (8) 
                    
S=0.0D0
CASE (9) 
                    
S=0.0D0
CASE (10) 
                    
S=0.0D0
CASE (11) 
                    
S=0.0D0
CASE (12) 
                    
S=0.0D0
CASE (13) 
                    
S=0.0D0
CASE (14) 
                    
S=0.0D0
CASE (15) 
                    
S=0.0D0
CASE (16) 
                    
S=0.0D0
CASE (17) 
                    
S=0.0D0
CASE (18) 
                    
S=0.0D0
CASE (19) 
                    
S=0.0D0
CASE (20) 
                    
S=0.0D0
CASE (21) 
                    
S=0.0D0
CASE (22) 
                    
S=0.0D0
CASE (23) 
                    
S=0.0D0
CASE (24) 
                    
S=0.0D0
CASE (25) 
                    
S=0.0D0
CASE (26) 
                    
S=0.0D0
CASE (27) 
                    
S=0.0D0
CASE (28) 
                    
S=0.0D0
CASE (29) 
                    
S=0.0D0
CASE (30) 
                    
S=0.0D0
CASE (31) 
                    
S=0.0D0
CASE (32) 
                    
S=0.0D0
CASE (33) 
                    
S=0.0D0
CASE (34) 
                    
S=0.0D0
CASE (35) 
                    
S=0.0D0
CASE (36) 
                    
S=0.0D0
CASE (37) 
                    
S=0.0D0
CASE (38) 
                    
S=0.0D0
CASE (39) 
                    
S=0.0D0
CASE (40) 
                    
S=0.0D0
CASE (41) 
                    
S=0.0D0
CASE (42) 
                    
S=0.0D0
CASE (43) 
                    
S=0.0D0
CASE (44) 
                    
S=0.0D0
CASE (45) 
                    
S=0.0D0
CASE (46) 
                    
S=0.0D0
CASE (47) 
                    
S=0.0D0
CASE (48) 
                    
S=0.0D0
CASE (49) 
                    
S=0.0D0
CASE (50) 
                    
S=0.0D0
CASE (51) 
                    
S=0.0D0
CASE (52) 
                    
S=0.0D0
CASE (53) 
                    
S=0.0D0
CASE (54) 
                    
S=0.0D0
CASE (55) 
                    
S=0.0D0
CASE (56) 
                    
S=0.0D0
CASE (57) 
                    
S=0.0D0
CASE (58) 
                    
S=0.0D0
CASE (59) 
                    
S=0.0D0
CASE (60) 
                    
S=0.0D0
CASE (61) 
                    
S=0.0D0
CASE (62) 
                    
S=0.0D0
CASE (63) 
                    
S=0.0D0
CASE (64) 
                    
S=0.0D0
CASE (65) 
                    
S=0.0D0
CASE (66) 
                    
S=0.0D0
CASE (67) 
                    
S=0.0D0
CASE (68) 
                    
S=0.0D0
CASE (69) 
                    
S=0.0D0
CASE (70) 
                    
S=0.0D0
CASE (71) 
                    
S=0.0D0
CASE (72) 
                    
S=0.0D0
CASE (73) 
                    
S=0.0D0
CASE (74) 
                    
S=0.0D0
CASE (75) 
                    
S=0.0D0
CASE (76) 
                    
S=0.0D0
CASE (77) 
                    
S=0.0D0
CASE (78) 
                    
S=0.0D0
CASE (79) 
                    
S=0.0D0
CASE (80) 
                    
S=0.0D0
CASE (81) 
                    
S=0.0D0
CASE (82) 
                    
S=120
CASE (83) 
                    
S=0.0D0
END SELECT
DFYZ5=S
END FUNCTION DFYZ5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DFZ6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(1)  
                    
S=0.0D0
CASE(2)  
                    
S=0.0D0
CASE(3)  
                    
S=0.0D0
CASE(4)  
                    
S=0.0D0
CASE(5)  
                    
S=0.0D0
CASE(6)  
                    
S=0.0D0
CASE(7)  
                    
S=0.0D0
CASE(8)  
                    
S=0.0D0
CASE(9)  
                    
S=0.0D0
CASE(10)  
                    
S=0.0D0
CASE(11)  
                    
S=0.0D0
CASE(12)  
                    
S=0.0D0
CASE(13)  
                    
S=0.0D0
CASE(14)  
                    
S=0.0D0
CASE(15)  
                    
S=0.0D0
CASE(16)  
                    
S=0.0D0
CASE(17)  
                    
S=0.0D0
CASE(18)  
                    
S=0.0D0
CASE(19)  
                    
S=0.0D0
CASE(20)  
                    
S=0.0D0
CASE(21)  
                    
S=0.0D0
CASE(22)  
                    
S=0.0D0
CASE(23)  
                    
S=0.0D0
CASE(24)  
                    
S=0.0D0
CASE(25)  
                    
S=0.0D0
CASE(26)  
                    
S=0.0D0
CASE(27)  
                    
S=0.0D0
CASE(28)  
                    
S=0.0D0
CASE(29)  
                    
S=0.0D0
CASE(30)  
                    
S=0.0D0
CASE(31)  
                    
S=0.0D0
CASE(32)  
                    
S=0.0D0
CASE(33)  
                    
S=0.0D0
CASE(34)  
                    
S=0.0D0
CASE(35)  
                    
S=0.0D0
CASE(36)  
                    
S=0.0D0
CASE(37)  
                    
S=0.0D0
CASE(38)  
                    
S=0.0D0
CASE(39)  
                    
S=0.0D0
CASE(40)  
                    
S=0.0D0
CASE(41)  
                    
S=0.0D0
CASE(42)  
                    
S=0.0D0
CASE(43)  
                    
S=0.0D0
CASE(44)  
                    
S=0.0D0
CASE(45)  
                    
S=0.0D0
CASE(46)  
                    
S=0.0D0
CASE(47)  
                    
S=0.0D0
CASE(48)  
                    
S=0.0D0
CASE(49)  
                    
S=0.0D0
CASE(50)  
                    
S=0.0D0
CASE(51)  
                    
S=0.0D0
CASE(52)  
                    
S=0.0D0
CASE(53)  
                    
S=0.0D0
CASE(54)  
                    
S=0.0D0
CASE(55)  
                    
S=0.0D0
CASE(56)  
                    
S=0.0D0
CASE(57)  
                    
S=0.0D0
CASE(58)  
                    
S=0.0D0
CASE(59)  
                    
S=0.0D0
CASE(60)  
                    
S=0.0D0
CASE(61)  
                    
S=0.0D0
CASE(62)  
                    
S=0.0D0
CASE(63)  
                    
S=0.0D0
CASE(64)  
                    
S=0.0D0
CASE(65)  
                    
S=0.0D0
CASE(66)  
                    
S=0.0D0
CASE(67)  
                    
S=0.0D0
CASE(68)  
                    
S=0.0D0
CASE(69)  
                    
S=0.0D0
CASE(70)  
                    
S=0.0D0
CASE(71)  
                    
S=0.0D0
CASE(72)  
                    
S=0.0D0
CASE(73)  
                    
S=0.0D0
CASE(74)  
                    
S=0.0D0
CASE(75)  
                    
S=0.0D0
CASE(76)  
                    
S=0.0D0
CASE(77)  
                    
S=0.0D0
CASE(78)  
                    
S=0.0D0
CASE(79)  
                    
S=0.0D0
CASE(80)  
                    
S=0.0D0
CASE(81)  
                    
S=0.0D0
CASE(82)  
                    
S=0.0D0
CASE(83)  
                    
S=720
END SELECT
DFZ6=S
END FUNCTION DFZ6
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL FUNCTION DLX(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)

CASE(1) 
S= 2
CASE(6) 
S= -2 + 4.d0*zd1
CASE(19) 
S= 2 - 12.d0*zd1 + 12.d0*zd1**2
CASE(29) 
S= -2 + 24.d0*zd1 - 60*zd1**2 + 40*zd1**3
CASE(55) 
S= 2 - 40*zd1 + 180*zd1**2 - 280*zd1**3 + 140*zd1**4
CASE(76) 
S= -2 + 60*zd1 - 420*zd1**2 + 1120*zd1**3 - 1260*zd1**4 + 504.d0*zd1**5
CASE(5) 
S= -2 + 4.d0*yd1
CASE(17) 
S= 2 - 4.d0*yd1 - 4.d0*zd1 + 8.d0*yd1*zd1
CASE(28) 
S= -2 + 4.d0*yd1 + 12.d0*zd1 - 24.d0*yd1*zd1 - 12.d0*zd1**2 + 24.d0*yd1*zd1**2
CASE(47) 
S= 2 - 4.d0*yd1 - 24.d0*zd1 + 48.d0*yd1*zd1 + 60*zd1**2 - 120*yd1*zd1**2 - 40*zd1**3 + 80*yd1*zd1**3
CASE(75) 
S= -2 + 4.d0*yd1 + 40*zd1 - 80*yd1*zd1 - 180*zd1**2 + 360*yd1*zd1**2 + 280*zd1**3 - 560*yd1*zd1**3 - 140*zd1**4 + 280*yd1*zd1**4
CASE(12) 
S= 2 - 12.d0*yd1 + 12.d0*yd1**2
CASE(27) 
S= -2 + 12.d0*yd1 - 12.d0*yd1**2 + 4.d0*zd1 - 24.d0*yd1*zd1 + 24.d0*yd1**2*zd1
CASE(48) 
S= 2 - 12.d0*yd1 + 12.d0*yd1**2 - 12.d0*zd1 + 72.d0*yd1*zd1 - 72.d0*yd1**2*zd1 + 12.d0*zd1**2 - 72.d0*yd1*zd1**2 + 72.d0*yd1**2*zd1**2
CASE(74) 
S= -2 + 12.d0*yd1 - 12.d0*yd1**2 + 24.d0*zd1 - 144.d0*yd1*zd1 + 144.d0*yd1**2*zd1 - 60*zd1**2 + 360*yd1*zd1**2 - 360*yd1**2*zd1**2 &
+ 40*zd1**3 - 240*yd1*zd1**3 + 240*yd1**2*zd1**3
CASE(26) 
S= -2 + 24.d0*yd1 - 60*yd1**2 + 40*yd1**3
CASE(46) 
S= 2 - 24.d0*yd1 + 60*yd1**2 - 40*yd1**3 - 4.d0*zd1 + 48.d0*yd1*zd1 - 120*yd1**2*zd1 + 80*yd1**3*zd1
CASE(73) 
S= -2 + 24.d0*yd1 - 60*yd1**2 + 40*yd1**3 + 12.d0*zd1 - 144.d0*yd1*zd1 + 360*yd1**2*zd1 - 240*yd1**3*zd1 - 12.d0*zd1**2 &
+ 144.d0*yd1*zd1**2 - 360*yd1**2*zd1**2 + 240*yd1**3*zd1**2
CASE(45) 
S= 2 - 40*yd1 + 180*yd1**2 - 280*yd1**3 + 140*yd1**4
CASE(72) 
S= -2 + 40*yd1 - 180*yd1**2 + 280*yd1**3 - 140*yd1**4 + 4.d0*zd1 - 80*yd1*zd1 + 360*yd1**2*zd1 - 560*yd1**3*zd1 + 280*yd1**4*zd1
CASE(71) 
S= -2 + 60*yd1 - 420*yd1**2 + 1120*yd1**3 - 1260*yd1**4 + 504.d0*yd1**5
CASE(4) 
S= -6 + 12.d0*xd1
CASE(18) 
S= 6 - 12.d0*xd1 - 12.d0*zd1 + 24.d0*xd1*zd1
CASE(25) 
S= -6 + 12.d0*xd1 + 36.d0*zd1 - 72.d0*xd1*zd1 - 36.d0*zd1**2 + 72.d0*xd1*zd1**2
CASE(44) 
S= 6 - 12.d0*xd1 - 72.d0*zd1 + 144.d0*xd1*zd1 + 180*zd1**2 - 360*xd1*zd1**2 - 120*zd1**3 + 240*xd1*zd1**3
CASE(70) 
S= -6 + 12.d0*xd1 + 120*zd1 - 240*xd1*zd1 - 540*zd1**2 + 1080*xd1*zd1**2 + 840*zd1**3 - 1680*xd1*zd1**3 - 420*zd1**4 &
+ 840*xd1*zd1**4
CASE(11) 
S= 6 - 12.d0*xd1 - 12.d0*yd1 + 24.d0*xd1*yd1
CASE(24) 
S= -6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 12.d0*zd1 - 24.d0*xd1*zd1 - 24.d0*yd1*zd1 + 48.d0*xd1*yd1*zd1
CASE(43) 
S= 6 - 12.d0*xd1 - 12.d0*yd1 + 24.d0*xd1*yd1 - 36.d0*zd1 + 72.d0*xd1*zd1 + 72.d0*yd1*zd1 - 144.d0*xd1*yd1*zd1 + 36.d0*zd1**2 - 72.d0*xd1*zd1**2 &
- 72.d0*yd1*zd1**2 + 144.d0*xd1*yd1*zd1**2
CASE(69) 
S= -6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 - 180*zd1**2 + 360*xd1*zd1**2 &
+ 360*yd1*zd1**2 - 720*xd1*yd1*zd1**2 + 120*zd1**3 - 240*xd1*zd1**3 - 240*yd1*zd1**3 + 480*xd1*yd1*zd1**3
CASE(23) 
S= -6 + 12.d0*xd1 + 36.d0*yd1 - 72.d0*xd1*yd1 - 36.d0*yd1**2 + 72.d0*xd1*yd1**2
CASE(42) 
S= 6 - 12.d0*xd1 - 36.d0*yd1 + 72.d0*xd1*yd1 + 36.d0*yd1**2 - 72.d0*xd1*yd1**2 - 12.d0*zd1 + 24.d0*xd1*zd1 + 72.d0*yd1*zd1 - 144.d0*xd1*yd1*zd1 &
- 72.d0*yd1**2*zd1 + 144.d0*xd1*yd1**2*zd1
CASE(68) 
S= -6 + 12.d0*xd1 + 36.d0*yd1 - 72.d0*xd1*yd1 - 36.d0*yd1**2 + 72.d0*xd1*yd1**2 + 36.d0*zd1 - 72.d0*xd1*zd1 - 216.d0*yd1*zd1 + 432.d0*xd1*yd1*zd1 &
+ 216.d0*yd1**2*zd1 - 432.d0*xd1*yd1**2*zd1 - 36.d0*zd1**2 + 72.d0*xd1*zd1**2 + 216.d0*yd1*zd1**2 - 432.d0*xd1*yd1*zd1**2 - 216.d0*yd1**2*zd1**2&
 + 432.d0*xd1*yd1**2*zd1**2
CASE(41) 
S= 6 - 12.d0*xd1 - 72.d0*yd1 + 144.d0*xd1*yd1 + 180*yd1**2 - 360*xd1*yd1**2 - 120*yd1**3 + 240*xd1*yd1**3
CASE(67) 
S= -6 + 12.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 180*yd1**2 + 360*xd1*yd1**2 + 120*yd1**3 - 240*xd1*yd1**3 + 12.d0*zd1 - 24.d0*xd1*zd1 &
- 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 + 360*yd1**2*zd1 - 720*xd1*yd1**2*zd1 - 240*yd1**3*zd1 + 480*xd1*yd1**3*zd1
CASE(66) 
S= -6 + 12.d0*xd1 + 120*yd1 - 240*xd1*yd1 - 540*yd1**2 + 1080*xd1*yd1**2 + 840*yd1**3 - 1680*xd1*yd1**3 - 420*yd1**4&
 + 840*xd1*yd1**4
CASE(10) 
S= 12 - 60*xd1 + 60*xd1**2
CASE(22) 
S= -12 + 60*xd1 - 60*xd1**2 + 24.d0*zd1 - 120*xd1*zd1 + 120*xd1**2*zd1
CASE(39) 
S= 12 - 60*xd1 + 60*xd1**2 - 72.d0*zd1 + 360*xd1*zd1 - 360*xd1**2*zd1 + 72.d0*zd1**2 - 360*xd1*zd1**2 + 360*xd1**2*zd1**2
CASE(65) 
S= -12 + 60*xd1 - 60*xd1**2 + 144.d0*zd1 - 720*xd1*zd1 + 720*xd1**2*zd1 - 360*zd1**2 + 1800*xd1*zd1**2 - 1800*xd1**2*zd1**2&
 + 240*zd1**3 - 1200*xd1*zd1**3 + 1200*xd1**2*zd1**3
CASE(21) 
S= -12 + 60*xd1 - 60*xd1**2 + 24.d0*yd1 - 120*xd1*yd1 + 120*xd1**2*yd1
CASE(40) 
S= 12 - 60*xd1 + 60*xd1**2 - 24.d0*yd1 + 120*xd1*yd1 - 120*xd1**2*yd1 - 24.d0*zd1 + 120*xd1*zd1 - 120*xd1**2*zd1 + 48.d0*yd1*zd1 &
- 240*xd1*yd1*zd1 + 240*xd1**2*yd1*zd1
CASE(64) 
S= -12 + 60*xd1 - 60*xd1**2 + 24.d0*yd1 - 120*xd1*yd1 + 120*xd1**2*yd1 + 72.d0*zd1 - 360*xd1*zd1 + 360*xd1**2*zd1 - 144.d0*yd1*zd1 &
+ 720*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 - 72.d0*zd1**2 + 360*xd1*zd1**2 - 360*xd1**2*zd1**2 + 144.d0*yd1*zd1**2 - 720*xd1*yd1*zd1**2&
 + 720*xd1**2*yd1*zd1**2
CASE(38) 
S= 12 - 60*xd1 + 60*xd1**2 - 72.d0*yd1 + 360*xd1*yd1 - 360*xd1**2*yd1 + 72.d0*yd1**2 - 360*xd1*yd1**2 + 360*xd1**2*yd1**2
CASE(63) 
S= -12 + 60*xd1 - 60*xd1**2 + 72.d0*yd1 - 360*xd1*yd1 + 360*xd1**2*yd1 - 72.d0*yd1**2 + 360*xd1*yd1**2 - 360*xd1**2*yd1**2 + 24.d0*zd1&
 - 120*xd1*zd1 + 120*xd1**2*zd1 - 144.d0*yd1*zd1 + 720*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 + 144.d0*yd1**2*zd1 - 720*xd1*yd1**2*zd1 &
+ 720*xd1**2*yd1**2*zd1
CASE(62) 
S= -12 + 60*xd1 - 60*xd1**2 + 144.d0*yd1 - 720*xd1*yd1 + 720*xd1**2*yd1 - 360*yd1**2 + 1800*xd1*yd1**2 - 1800*xd1**2*yd1**2 &
+ 240*yd1**3 - 1200*xd1*yd1**3 + 1200*xd1**2*yd1**3
CASE(20) 
S= -20.0D0+ 180*xd1 - 420*xd1**2 + 280*xd1**3
CASE(37) 
S= 20.0D0- 180*xd1 + 420*xd1**2 - 280*xd1**3 - 40*zd1 + 360*xd1*zd1 - 840*xd1**2*zd1 + 560*xd1**3*zd1
CASE(61) 
S= -20.0D0+ 180*xd1 - 420*xd1**2 + 280*xd1**3 + 120*zd1 - 1080*xd1*zd1 + 2520*xd1**2*zd1 - 1680*xd1**3*zd1 - 120*zd1**2 &
+ 1080*xd1*zd1**2 - 2520*xd1**2*zd1**2 + 1680*xd1**3*zd1**2
CASE(36) 
S= 20.0D0- 180*xd1 + 420*xd1**2 - 280*xd1**3 - 40*yd1 + 360*xd1*yd1 - 840*xd1**2*yd1 + 560*xd1**3*yd1
CASE(60) 
S= -20.0D0+ 180*xd1 - 420*xd1**2 + 280*xd1**3 + 40*yd1 - 360*xd1*yd1 + 840*xd1**2*yd1 - 560*xd1**3*yd1 + 40*zd1 - 360*xd1*zd1&
 + 840*xd1**2*zd1 - 560*xd1**3*zd1 - 80*yd1*zd1 + 720*xd1*yd1*zd1 - 1680*xd1**2*yd1*zd1 + 1120*xd1**3*yd1*zd1
CASE(59) 
S= -20.0D0+ 180*xd1 - 420*xd1**2 + 280*xd1**3 + 120*yd1 - 1080*xd1*yd1 + 2520*xd1**2*yd1 - 1680*xd1**3*yd1 - 120*yd1**2 &
+ 1080*xd1*yd1**2 - 2520*xd1**2*yd1**2 + 1680*xd1**3*yd1**2
CASE(35) 
S= 30.0D0- 420*xd1 + 1680*xd1**2 - 2520*xd1**3 + 1260*xd1**4
CASE(58) 
S= -30.0D0+ 420*xd1 - 1680*xd1**2 + 2520*xd1**3 - 1260*xd1**4 + 60*zd1 - 840*xd1*zd1 + 3360*xd1**2*zd1 - 5040*xd1**3*zd1 &
+ 2520*xd1**4*zd1
CASE(57)  
S=-30.0D0+ 420*xd1 - 1680*xd1**2 + 2520*xd1**3 - 1260*xd1**4 + 60*yd1 - 840*xd1*yd1 + 3360*xd1**2*yd1 - 5040*xd1**3*yd1 &
+ 2520*xd1**4*yd1
CASE(56) 
S= -42 + 840*xd1 - 5040*xd1**2 + 12600*xd1**3 - 13860*xd1**4 + 5544.d0*xd1**5

case default
s=0.0d0


END SELECT  
	DLX=S    
	END FUNCTION DLX 


REAL FUNCTION DLY(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)

CASE(2)
 S=2
CASE(8)
 S=-2 + 4.d0*zd1
CASE(15)
 S=2 - 12.d0*zd1 + 12.d0*zd1**2
CASE(32)
 S=-2 + 24.d0*zd1 - 60*zd1**2 + 40*zd1**3
CASE(53)
 S=2 - 40*zd1 + 180*zd1**2 - 280*zd1**3 + 140*zd1**4
CASE(82)
 S=-2 + 60*zd1 - 420*zd1**2 + 1120*zd1**3 - 1260*zd1**4 + 504.d0*zd1**5
CASE(7)
 S=-6 + 12.d0*yd1
CASE(14)
 S=6 - 12.d0*yd1 - 12.d0*zd1 + 24.d0*yd1*zd1
CASE(31)
 S=-6 + 12.d0*yd1 + 36.d0*zd1 - 72.d0*yd1*zd1 - 36.d0*zd1**2 + 72.d0*yd1*zd1**2
CASE(52)
 S=6 - 12.d0*yd1 - 72.d0*zd1 + 144.d0*yd1*zd1 + 180*zd1**2 - 360*yd1*zd1**2 - 120*zd1**3 + 240*yd1*zd1**3
CASE(81)
 S=-6 + 12.d0*yd1 + 120*zd1 - 240*yd1*zd1 - 540*zd1**2 + 1080*yd1*zd1**2 + 840*zd1**3 - 1680*yd1*zd1**3 - 420*zd1**4 &
+ 840*yd1*zd1**4
CASE(13)
 S=12 - 60*yd1 + 60*yd1**2
CASE(30)
 S=-12 + 60*yd1 - 60*yd1**2 + 24.d0*zd1 - 120*yd1*zd1 + 120*yd1**2*zd1
CASE(51)
 S=12 - 60*yd1 + 60*yd1**2 - 72.d0*zd1 + 360*yd1*zd1 - 360*yd1**2*zd1 + 72.d0*zd1**2 - 360*yd1*zd1**2 + 360*yd1**2*zd1**2
CASE(80)
 S=-12 + 60*yd1 - 60*yd1**2 + 144.d0*zd1 - 720*yd1*zd1 + 720*yd1**2*zd1 - 360*zd1**2 + 1800*yd1*zd1**2 - 1800*yd1**2*zd1**2 &
+ 240*zd1**3 - 1200*yd1*zd1**3 + 1200*yd1**2*zd1**3
CASE(34)
 S=-20.0D0+ 180*yd1 - 420*yd1**2 + 280*yd1**3
CASE(50)
 S=20.0D0- 180*yd1 + 420*yd1**2 - 280*yd1**3 - 40*zd1 + 360*yd1*zd1 - 840*yd1**2*zd1 + 560*yd1**3*zd1
CASE(79)
 S=-20.0D0+ 180*yd1 - 420*yd1**2 + 280*yd1**3 + 120*zd1 - 1080*yd1*zd1 + 2520*yd1**2*zd1 - 1680*yd1**3*zd1 - 120*zd1**2&
 + 1080*yd1*zd1**2 - 2520*yd1**2*zd1**2 + 1680*yd1**3*zd1**2
CASE(49)
 S=30.0D0- 420*yd1 + 1680*yd1**2 - 2520*yd1**3 + 1260*yd1**4
CASE(78)
 S=-30.0D0+ 420*yd1 - 1680*yd1**2 + 2520*yd1**3 - 1260*yd1**4 + 60*zd1 - 840*yd1*zd1 + 3360*yd1**2*zd1 - 5040*yd1**3*zd1&
 + 2520*yd1**4*zd1
CASE(77)
 S=-42 + 840*yd1 - 5040*yd1**2 + 12600*yd1**3 - 13860*yd1**4 + 5544.d0*yd1**5

CASE(5)
 S=-2 + 4.d0*xd1
CASE(17)
 S=2 - 4.d0*xd1 - 4.d0*zd1 + 8.d0*xd1*zd1
CASE(28)
 S=-2 + 4.d0*xd1 + 12.d0*zd1 - 24.d0*xd1*zd1 - 12.d0*zd1**2 + 24.d0*xd1*zd1**2
CASE(47)
 S=2 - 4.d0*xd1 - 24.d0*zd1 + 48.d0*xd1*zd1 + 60*zd1**2 - 120*xd1*zd1**2 - 40*zd1**3 + 80*xd1*zd1**3
CASE(75)
 S=-2 + 4.d0*xd1 + 40*zd1 - 80*xd1*zd1 - 180*zd1**2 + 360*xd1*zd1**2 + 280*zd1**3 - 560*xd1*zd1**3 - 140*zd1**4 + 280*xd1*zd1**4
CASE(12)
 S=6 - 12.d0*xd1 - 12.d0*yd1 + 24.d0*xd1*yd1
CASE(27)
 S=-6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 12.d0*zd1 - 24.d0*xd1*zd1 - 24.d0*yd1*zd1 + 48.d0*xd1*yd1*zd1
CASE(48)
 S=6 - 12.d0*xd1 - 12.d0*yd1 + 24.d0*xd1*yd1 - 36.d0*zd1 + 72.d0*xd1*zd1 + 72.d0*yd1*zd1 - 144.d0*xd1*yd1*zd1 + 36.d0*zd1**2 - 72.d0*xd1*zd1**2 &
- 72.d0*yd1*zd1**2 + 144.d0*xd1*yd1*zd1**2
CASE(74)
 S=-6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 - 180*zd1**2 + 360*xd1*zd1**2 &
+ 360*yd1*zd1**2 - 720*xd1*yd1*zd1**2 + 120*zd1**3 - 240*xd1*zd1**3 - 240*yd1*zd1**3 + 480*xd1*yd1*zd1**3
CASE(26)
 S=-12 + 24.d0*xd1 + 60*yd1 - 120*xd1*yd1 - 60*yd1**2 + 120*xd1*yd1**2
CASE(46)
 S=12 - 24.d0*xd1 - 60*yd1 + 120*xd1*yd1 + 60*yd1**2 - 120*xd1*yd1**2 - 24.d0*zd1 + 48.d0*xd1*zd1 + 120*yd1*zd1 - 240*xd1*yd1*zd1 &
- 120*yd1**2*zd1 + 240*xd1*yd1**2*zd1
CASE(73)
 S=-12 + 24.d0*xd1 + 60*yd1 - 120*xd1*yd1 - 60*yd1**2 + 120*xd1*yd1**2 + 72.d0*zd1 - 144.d0*xd1*zd1 - 360*yd1*zd1 + 720*xd1*yd1*zd1 &
+ 360*yd1**2*zd1 - 720*xd1*yd1**2*zd1 - 72.d0*zd1**2 + 144.d0*xd1*zd1**2 + 360*yd1*zd1**2 - 720*xd1*yd1*zd1**2 - 360*yd1**2*zd1**2 &
+ 720*xd1*yd1**2*zd1**2
CASE(45)
 S=20.0D0- 40*xd1 - 180*yd1 + 360*xd1*yd1 + 420*yd1**2 - 840*xd1*yd1**2 - 280*yd1**3 + 560*xd1*yd1**3
CASE(72)
 S=-20.0D0+ 40*xd1 + 180*yd1 - 360*xd1*yd1 - 420*yd1**2 + 840*xd1*yd1**2 + 280*yd1**3 - 560*xd1*yd1**3 + 40*zd1 - 80*xd1*zd1 &
- 360*yd1*zd1 + 720*xd1*yd1*zd1 + 840*yd1**2*zd1 - 1680*xd1*yd1**2*zd1 - 560*yd1**3*zd1 + 1120*xd1*yd1**3*zd1
CASE(71)
 S=-30.0D0+ 60*xd1 + 420*yd1 - 840*xd1*yd1 - 1680*yd1**2 + 3360*xd1*yd1**2 + 2520*yd1**3 - 5040*xd1*yd1**3 - 1260*yd1**4 &
+ 2520*xd1*yd1**4
CASE(11)
 S=2 - 12.d0*xd1 + 12.d0*xd1**2
CASE(24)
 S=-2 + 12.d0*xd1 - 12.d0*xd1**2 + 4.d0*zd1 - 24.d0*xd1*zd1 + 24.d0*xd1**2*zd1
CASE(43)
 S=2 - 12.d0*xd1 + 12.d0*xd1**2 - 12.d0*zd1 + 72.d0*xd1*zd1 - 72.d0*xd1**2*zd1 + 12.d0*zd1**2 - 72.d0*xd1*zd1**2 + 72.d0*xd1**2*zd1**2
CASE(69)
 S=-2 + 12.d0*xd1 - 12.d0*xd1**2 + 24.d0*zd1 - 144.d0*xd1*zd1 + 144.d0*xd1**2*zd1 - 60*zd1**2 + 360*xd1*zd1**2 - 360*xd1**2*zd1**2 &
+ 40*zd1**3 - 240*xd1*zd1**3 + 240*xd1**2*zd1**3
CASE(23)
 S=-6 + 36.d0*xd1 - 36.d0*xd1**2 + 12.d0*yd1 - 72.d0*xd1*yd1 + 72.d0*xd1**2*yd1
CASE(42)
 S=6 - 36.d0*xd1 + 36.d0*xd1**2 - 12.d0*yd1 + 72.d0*xd1*yd1 - 72.d0*xd1**2*yd1 - 12.d0*zd1 + 72.d0*xd1*zd1 - 72.d0*xd1**2*zd1 + 24.d0*yd1*zd1 &
- 144.d0*xd1*yd1*zd1 + 144.d0*xd1**2*yd1*zd1
CASE(68)
 S=-6 + 36.d0*xd1 - 36.d0*xd1**2 + 12.d0*yd1 - 72.d0*xd1*yd1 + 72.d0*xd1**2*yd1 + 36.d0*zd1 - 216.d0*xd1*zd1 + 216.d0*xd1**2*zd1 - 72.d0*yd1*zd1 &
+ 432.d0*xd1*yd1*zd1 - 432.d0*xd1**2*yd1*zd1 - 36.d0*zd1**2 + 216.d0*xd1*zd1**2 - 216.d0*xd1**2*zd1**2 + 72.d0*yd1*zd1**2 - 432.d0*xd1*yd1*zd1**2&
 + 432.d0*xd1**2*yd1*zd1**2
CASE(41)
 S=12 - 72.d0*xd1 + 72.d0*xd1**2 - 60*yd1 + 360*xd1*yd1 - 360*xd1**2*yd1 + 60*yd1**2 - 360*xd1*yd1**2 + 360*xd1**2*yd1**2
CASE(67)
 S=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 60*yd1 - 360*xd1*yd1 + 360*xd1**2*yd1 - 60*yd1**2 + 360*xd1*yd1**2 - 360*xd1**2*yd1**2 &
+ 24.d0*zd1 - 144.d0*xd1*zd1 + 144.d0*xd1**2*zd1 - 120*yd1*zd1 + 720*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 + 120*yd1**2*zd1 &
- 720*xd1*yd1**2*zd1 + 720*xd1**2*yd1**2*zd1
CASE(66)
 S=-20.0D0+ 120*xd1 - 120*xd1**2 + 180*yd1 - 1080*xd1*yd1 + 1080*xd1**2*yd1 - 420*yd1**2 + 2520*xd1*yd1**2 - 2520*xd1**2*yd1**2&
 + 280*yd1**3 - 1680*xd1*yd1**3 + 1680*xd1**2*yd1**3

CASE(21)
 S=-2 + 24.d0*xd1 - 60*xd1**2 + 40*xd1**3
CASE(40)
 S=2 - 24.d0*xd1 + 60*xd1**2 - 40*xd1**3 - 4.d0*zd1 + 48.d0*xd1*zd1 - 120*xd1**2*zd1 + 80*xd1**3*zd1
CASE(64)
 S=-2 + 24.d0*xd1 - 60*xd1**2 + 40*xd1**3 + 12.d0*zd1 - 144.d0*xd1*zd1 + 360*xd1**2*zd1 - 240*xd1**3*zd1 - 12.d0*zd1**2 + 144.d0*xd1*zd1**2&
 - 360*xd1**2*zd1**2 + 240*xd1**3*zd1**2
CASE(38)
 S=6 - 72.d0*xd1 + 180*xd1**2 - 120*xd1**3 - 12.d0*yd1 + 144.d0*xd1*yd1 - 360*xd1**2*yd1 + 240*xd1**3*yd1
CASE(63)
 S=-6 + 72.d0*xd1 - 180*xd1**2 + 120*xd1**3 + 12.d0*yd1 - 144.d0*xd1*yd1 + 360*xd1**2*yd1 - 240*xd1**3*yd1 + 12.d0*zd1 - 144.d0*xd1*zd1 &
+ 360*xd1**2*zd1 - 240*xd1**3*zd1 - 24.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 + 480*xd1**3*yd1*zd1
CASE(62)
 S=-12 + 144.d0*xd1 - 360*xd1**2 + 240*xd1**3 + 60*yd1 - 720*xd1*yd1 + 1800*xd1**2*yd1 - 1200*xd1**3*yd1 - 60*yd1**2 &
+ 720*xd1*yd1**2 - 1800*xd1**2*yd1**2 + 1200*xd1**3*yd1**2

CASE(36)
 S=2 - 40*xd1 + 180*xd1**2 - 280*xd1**3 + 140*xd1**4
CASE(60)
 S=-2 + 40*xd1 - 180*xd1**2 + 280*xd1**3 - 140*xd1**4 + 4.d0*zd1 - 80*xd1*zd1 + 360*xd1**2*zd1 - 560*xd1**3*zd1 + 280*xd1**4*zd1
CASE(59)
 S=-6 + 120*xd1 - 540*xd1**2 + 840*xd1**3 - 420*xd1**4 + 12.d0*yd1 - 240*xd1*yd1 + 1080*xd1**2*yd1 - 1680*xd1**3*yd1 + 840*xd1**4*yd1

CASE(57)
 S=-2 + 60*xd1 - 420*xd1**2 + 1120*xd1**3 - 1260*xd1**4 + 504.d0*xd1**5
CASE default
 S=0.0D0
  END SELECT
 DLY=S    
	END FUNCTION DLY
	
REAL FUNCTION DLZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=2
CASE(9)
 S=-6 + 12.d0*zd1
CASE(16)
 S=12 - 60*zd1 + 60*zd1**2
CASE(33)
 S=-20.0D0+ 180*zd1 - 420*zd1**2 + 280*zd1**3
CASE(54)
 S=30.0D0- 420*zd1 + 1680*zd1**2 - 2520*zd1**3 + 1260*zd1**4
CASE(83)
 S=-42 + 840*zd1 - 5040*zd1**2 + 12600*zd1**3 - 13860*zd1**4 + 5544.d0*zd1**5

CASE(8)
 S=-2 + 4.d0*yd1
CASE(15)
 S=6 - 12.d0*yd1 - 12.d0*zd1 + 24.d0*yd1*zd1
CASE(32)
 S=-12 + 24.d0*yd1 + 60*zd1 - 120*yd1*zd1 - 60*zd1**2 + 120*yd1*zd1**2
CASE(53)
 S=20.0D0- 40*yd1 - 180*zd1 + 360*yd1*zd1 + 420*zd1**2 - 840*yd1*zd1**2 - 280*zd1**3 + 560*yd1*zd1**3
CASE(82)
 S=-30.0D0+ 60*yd1 + 420*zd1 - 840*yd1*zd1 - 1680*zd1**2 + 3360*yd1*zd1**2 + 2520*zd1**3 - 5040*yd1*zd1**3 - 1260*zd1**4 &
+ 2520*yd1*zd1**4

CASE(14)
 S=2 - 12.d0*yd1 + 12.d0*yd1**2
CASE(31)
 S=-6 + 36.d0*yd1 - 36.d0*yd1**2 + 12.d0*zd1 - 72.d0*yd1*zd1 + 72.d0*yd1**2*zd1
CASE(52)
 S=12 - 72.d0*yd1 + 72.d0*yd1**2 - 60*zd1 + 360*yd1*zd1 - 360*yd1**2*zd1 + 60*zd1**2 - 360*yd1*zd1**2 + 360*yd1**2*zd1**2
CASE(81)
 S=-20.0D0+ 120*yd1 - 120*yd1**2 + 180*zd1 - 1080*yd1*zd1 + 1080*yd1**2*zd1 - 420*zd1**2 + 2520*yd1*zd1**2 - 2520*yd1**2*zd1**2&
 + 280*zd1**3 - 1680*yd1*zd1**3 + 1680*yd1**2*zd1**3
CASE(13)
 
CASE(30)
 S=-2 + 24.d0*yd1 - 60*yd1**2 + 40*yd1**3
CASE(51)
 S=6 - 72.d0*yd1 + 180*yd1**2 - 120*yd1**3 - 12.d0*zd1 + 144.d0*yd1*zd1 - 360*yd1**2*zd1 + 240*yd1**3*zd1
CASE(80)
 S=-12 + 144.d0*yd1 - 360*yd1**2 + 240*yd1**3 + 60*zd1 - 720*yd1*zd1 + 1800*yd1**2*zd1 - 1200*yd1**3*zd1 - 60*zd1**2 &
+ 720*yd1*zd1**2 - 1800*yd1**2*zd1**2 + 1200*yd1**3*zd1**2

CASE(50)
 S=2 - 40*yd1 + 180*yd1**2 - 280*yd1**3 + 140*yd1**4
CASE(79)
 S=-6 + 120*yd1 - 540*yd1**2 + 840*yd1**3 - 420*yd1**4 + 12.d0*zd1 - 240*yd1*zd1 + 1080*yd1**2*zd1 - 1680*yd1**3*zd1 &
+ 840*yd1**4*zd1

CASE(78)
 S=-2 + 60*yd1 - 420*yd1**2 + 1120*yd1**3 - 1260*yd1**4 + 504.d0*yd1**5


CASE(6)
 S=-2 + 4.d0*xd1
CASE(19)
 S=6 - 12.d0*xd1 - 12.d0*zd1 + 24.d0*xd1*zd1
CASE(29)
 S=-12 + 24.d0*xd1 + 60*zd1 - 120*xd1*zd1 - 60*zd1**2 + 120*xd1*zd1**2
CASE(55)
 S=20.0D0- 40*xd1 - 180*zd1 + 360*xd1*zd1 + 420*zd1**2 - 840*xd1*zd1**2 - 280*zd1**3 + 560*xd1*zd1**3
CASE(76)
 S=-30.0D0+ 60*xd1 + 420*zd1 - 840*xd1*zd1 - 1680*zd1**2 + 3360*xd1*zd1**2 + 2520*zd1**3 - 5040*xd1*zd1**3 &
- 1260*zd1**4 + 2520*xd1*zd1**4

CASE(17)
 S=2 - 4.d0*xd1 - 4.d0*yd1 + 8.d0*xd1*yd1
CASE(28)
 S=-6 + 12.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 + 12.d0*zd1 - 24.d0*xd1*zd1 - 24.d0*yd1*zd1 + 48.d0*xd1*yd1*zd1
CASE(47)
 S=12 - 24.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1 - 60*zd1 + 120*xd1*zd1 + 120*yd1*zd1 - 240*xd1*yd1*zd1 + 60*zd1**2 &
- 120*xd1*zd1**2 - 120*yd1*zd1**2 + 240*xd1*yd1*zd1**2
CASE(75)
 S=-20.0D0+ 40*xd1 + 40*yd1 - 80*xd1*yd1 + 180*zd1 - 360*xd1*zd1 - 360*yd1*zd1 + 720*xd1*yd1*zd1 - 420*zd1**2&
 + 840*xd1*zd1**2 + 840*yd1*zd1**2 - 1680*xd1*yd1*zd1**2 + 280*zd1**3 - 560*xd1*zd1**3 - 560*yd1*zd1**3&
 + 1120*xd1*yd1*zd1**3

CASE(27)
 S=-2 + 4.d0*xd1 + 12.d0*yd1 - 24.d0*xd1*yd1 - 12.d0*yd1**2 + 24.d0*xd1*yd1**2
CASE(48)
 S=6 - 12.d0*xd1 - 36.d0*yd1 + 72.d0*xd1*yd1 + 36.d0*yd1**2 - 72.d0*xd1*yd1**2 - 12.d0*zd1 + 24.d0*xd1*zd1 + 72.d0*yd1*zd1 - 144.d0*xd1*yd1*zd1 &
- 72.d0*yd1**2*zd1 + 144.d0*xd1*yd1**2*zd1
CASE(74)
 S=-12 + 24.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 72.d0*yd1**2 + 144.d0*xd1*yd1**2 + 60*zd1 - 120*xd1*zd1 - 360*yd1*zd1 &
+ 720*xd1*yd1*zd1 + 360*yd1**2*zd1 - 720*xd1*yd1**2*zd1 - 60*zd1**2 + 120*xd1*zd1**2 + 360*yd1*zd1**2 &
- 720*xd1*yd1*zd1**2 - 360*yd1**2*zd1**2 + 720*xd1*yd1**2*zd1**2

CASE(46)
 S=2 - 4.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1 + 60*yd1**2 - 120*xd1*yd1**2 - 40*yd1**3 + 80*xd1*yd1**3
CASE(73)
 S=-6 + 12.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 180*yd1**2 + 360*xd1*yd1**2 + 120*yd1**3 - 240*xd1*yd1**3 + 12.d0*zd1 - 24.d0*xd1*zd1 &
- 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 + 360*yd1**2*zd1 - 720*xd1*yd1**2*zd1 - 240*yd1**3*zd1 + 480*xd1*yd1**3*zd1

CASE(72)
 S=-2 + 4.d0*xd1 + 40*yd1 - 80*xd1*yd1 - 180*yd1**2 + 360*xd1*yd1**2 + 280*yd1**3 - 560*xd1*yd1**3 - 140*yd1**4 &
+ 280*xd1*yd1**4


CASE(18)
 S=2 - 12.d0*xd1 + 12.d0*xd1**2
CASE(25)
 S=-6 + 36.d0*xd1 - 36.d0*xd1**2 + 12.d0*zd1 - 72.d0*xd1*zd1 + 72.d0*xd1**2*zd1
CASE(44)
 S=12 - 72.d0*xd1 + 72.d0*xd1**2 - 60*zd1 + 360*xd1*zd1 - 360*xd1**2*zd1 + 60*zd1**2 - 360*xd1*zd1**2 + 360*xd1**2*zd1**2
CASE(70)
 S=-20.0D0+ 120*xd1 - 120*xd1**2 + 180*zd1 - 1080*xd1*zd1 + 1080*xd1**2*zd1 - 420*zd1**2 + 2520*xd1*zd1**2 &
- 2520*xd1**2*zd1**2 + 280*zd1**3 - 1680*xd1*zd1**3 + 1680*xd1**2*zd1**3

CASE(24)
 S=-2 + 12.d0*xd1 - 12.d0*xd1**2 + 4.d0*yd1 - 24.d0*xd1*yd1 + 24.d0*xd1**2*yd1
CASE(43)
 S=6 - 36.d0*xd1 + 36.d0*xd1**2 - 12.d0*yd1 + 72.d0*xd1*yd1 - 72.d0*xd1**2*yd1 - 12.d0*zd1 + 72.d0*xd1*zd1 - 72.d0*xd1**2*zd1 + 24.d0*yd1*zd1 &
- 144.d0*xd1*yd1*zd1 + 144.d0*xd1**2*yd1*zd1
CASE(69)
 S=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*yd1 - 144.d0*xd1*yd1 + 144.d0*xd1**2*yd1 + 60*zd1 - 360*xd1*zd1 + 360*xd1**2*zd1 &
- 120*yd1*zd1 + 720*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 - 60*zd1**2 + 360*xd1*zd1**2 - 360*xd1**2*zd1**2 &
+ 120*yd1*zd1**2 - 720*xd1*yd1*zd1**2 + 720*xd1**2*yd1*zd1**2

CASE(42)
 S=2 - 12.d0*xd1 + 12.d0*xd1**2 - 12.d0*yd1 + 72.d0*xd1*yd1 - 72.d0*xd1**2*yd1 + 12.d0*yd1**2 - 72.d0*xd1*yd1**2 + 72.d0*xd1**2*yd1**2
CASE(68)
 S=-6 + 36.d0*xd1 - 36.d0*xd1**2 + 36.d0*yd1 - 216.d0*xd1*yd1 + 216.d0*xd1**2*yd1 - 36.d0*yd1**2 + 216.d0*xd1*yd1**2 &
- 216.d0*xd1**2*yd1**2 + 12.d0*zd1 - 72.d0*xd1*zd1 + 72.d0*xd1**2*zd1 - 72.d0*yd1*zd1 + 432.d0*xd1*yd1*zd1 - 432.d0*xd1**2*yd1*zd1 &
+ 72.d0*yd1**2*zd1 - 432.d0*xd1*yd1**2*zd1 + 432.d0*xd1**2*yd1**2*zd1

CASE(67)
 S=-2 + 12.d0*xd1 - 12.d0*xd1**2 + 24.d0*yd1 - 144.d0*xd1*yd1 + 144.d0*xd1**2*yd1 - 60*yd1**2 + 360*xd1*yd1**2 &
- 360*xd1**2*yd1**2 + 40*yd1**3 - 240*xd1*yd1**3 + 240*xd1**2*yd1**3

CASE(22)
 S=-2 + 24.d0*xd1 - 60*xd1**2 + 40*xd1**3
CASE(39)
 S=6 - 72.d0*xd1 + 180*xd1**2 - 120*xd1**3 - 12.d0*zd1 + 144.d0*xd1*zd1 - 360*xd1**2*zd1 + 240*xd1**3*zd1
CASE(65)
 S=-12 + 144.d0*xd1 - 360*xd1**2 + 240*xd1**3 + 60*zd1 - 720*xd1*zd1 + 1800*xd1**2*zd1 - 1200*xd1**3*zd1 &
- 60*zd1**2 + 720*xd1*zd1**2 - 1800*xd1**2*zd1**2 + 1200*xd1**3*zd1**2

CASE(40)
 S=2 - 24.d0*xd1 + 60*xd1**2 - 40*xd1**3 - 4.d0*yd1 + 48.d0*xd1*yd1 - 120*xd1**2*yd1 + 80*xd1**3*yd1
CASE(64)
 S=-6 + 72.d0*xd1 - 180*xd1**2 + 120*xd1**3 + 12.d0*yd1 - 144.d0*xd1*yd1 + 360*xd1**2*yd1 - 240*xd1**3*yd1 + 12.d0*zd1 &
- 144.d0*xd1*zd1 + 360*xd1**2*zd1 - 240*xd1**3*zd1 - 24.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1 - 720*xd1**2*yd1*zd1 + 480*xd1**3*yd1*zd1

CASE(63)
 S=-2 + 24.d0*xd1 - 60*xd1**2 + 40*xd1**3 + 12.d0*yd1 - 144.d0*xd1*yd1 + 360*xd1**2*yd1 - 240*xd1**3*yd1 - 12.d0*yd1**2 &
+ 144.d0*xd1*yd1**2 - 360*xd1**2*yd1**2 + 240*xd1**3*yd1**2

CASE(37)
 S=2 - 40*xd1 + 180*xd1**2 - 280*xd1**3 + 140*xd1**4
CASE(61)
 S=-6 + 120*xd1 - 540*xd1**2 + 840*xd1**3 - 420*xd1**4 + 12.d0*zd1 - 240*xd1*zd1 + 1080*xd1**2*zd1 &
- 1680*xd1**3*zd1 + 840*xd1**4*zd1

CASE(60)
 S=-2 + 40*xd1 - 180*xd1**2 + 280*xd1**3 - 140*xd1**4 + 4.d0*yd1 - 80*xd1*yd1 + 360*xd1**2*yd1 &
- 560*xd1**3*yd1 + 280*xd1**4*yd1

CASE(58)
 S=-2 + 60*xd1 - 420*xd1**2 + 1120*xd1**3 - 1260*xd1**4 + 504.d0*xd1**5

CASE default
 S=0.0D0
  END SELECT
 DLZ=S    
	END FUNCTION DLZ
 
REAL FUNCTION DLX2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=12
CASE(18)
 S=-12 + 24.d0*zd1
CASE(25)
 S=12 - 72.d0*zd1 + 72.d0*zd1**2
CASE(44)
 S=-12 + 144.d0*zd1 - 360*zd1**2 + 240*zd1**3
CASE(70)
 S=12 - 240*zd1 + 1080*zd1**2 - 1680*zd1**3 + 840*zd1**4
CASE(11)
 S=-12 + 24.d0*yd1
CASE(24)
 S=12 - 24.d0*yd1 - 24.d0*zd1 + 48.d0*yd1*zd1
CASE(43)
 S=-12 + 24.d0*yd1 + 72.d0*zd1 - 144.d0*yd1*zd1 - 72.d0*zd1**2 + 144.d0*yd1*zd1**2
CASE(69)
 S=12 - 24.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1 + 360*zd1**2 - 720*yd1*zd1**2 - 240*zd1**3 + 480*yd1*zd1**3
CASE(23)
 S=12 - 72.d0*yd1 + 72.d0*yd1**2
CASE(42)
 S=-12 + 72.d0*yd1 - 72.d0*yd1**2 + 24.d0*zd1 - 144.d0*yd1*zd1 + 144.d0*yd1**2*zd1
CASE(68)
 S=12 - 72.d0*yd1 + 72.d0*yd1**2 - 72.d0*zd1 + 432.d0*yd1*zd1 - 432.d0*yd1**2*zd1 + 72.d0*zd1**2 - 432.d0*yd1*zd1**2 + 432.d0*yd1**2*zd1**2
CASE(41)
 S=-12 + 144.d0*yd1 - 360*yd1**2 + 240*yd1**3
CASE(67)
 S=12 - 144.d0*yd1 + 360*yd1**2 - 240*yd1**3 - 24.d0*zd1 + 288.d0*yd1*zd1 - 720*yd1**2*zd1 + 480*yd1**3*zd1
CASE(66)
 S=12 - 240*yd1 + 1080*yd1**2 - 1680*yd1**3 + 840*yd1**4
CASE(10)
 S=-60.0D0+ 120*xd1
CASE(22)
 S=60.0D0- 120*xd1 - 120*zd1 + 240*xd1*zd1
CASE(39)
 S=-60.0D0+ 120*xd1 + 360*zd1 - 720*xd1*zd1 - 360*zd1**2 + 720*xd1*zd1**2
CASE(65)
 S=60.0D0- 120*xd1 - 720*zd1 + 1440*xd1*zd1 + 1800*zd1**2 - 3600*xd1*zd1**2 - 1200*zd1**3 + 2400*xd1*zd1**3
CASE(21)
 S=60.0D0- 120*xd1 - 120*yd1 + 240*xd1*yd1
CASE(40)
 S=-60.0D0+ 120*xd1 + 120*yd1 - 240*xd1*yd1 + 120*zd1 - 240*xd1*zd1 - 240*yd1*zd1 + 480*xd1*yd1*zd1
CASE(64)
 S=60.0D0- 120*xd1 - 120*yd1 + 240*xd1*yd1 - 360*zd1 + 720*xd1*zd1 + 720*yd1*zd1 - 1440*xd1*yd1*zd1 + 360*zd1**2 &
- 720*xd1*zd1**2 - 720*yd1*zd1**2 + 1440*xd1*yd1*zd1**2
CASE(38)
 S=-60.0D0+ 120*xd1 + 360*yd1 - 720*xd1*yd1 - 360*yd1**2 + 720*xd1*yd1**2
CASE(63)
 S=60.0D0- 120*xd1 - 360*yd1 + 720*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 120*zd1 + 240*xd1*zd1 + 720*yd1*zd1 &
- 1440*xd1*yd1*zd1 - 720*yd1**2*zd1 + 1440*xd1*yd1**2*zd1
CASE(62)
 S=60.0D0- 120*xd1 - 720*yd1 + 1440*xd1*yd1 + 1800*yd1**2 - 3600*xd1*yd1**2 - 1200*yd1**3 + 2400*xd1*yd1**3
CASE(20)
 S=180.0D0- 840*xd1 + 840*xd1**2
CASE(37)
 S=-180.0D0+ 840*xd1 - 840*xd1**2 + 360*zd1 - 1680*xd1*zd1 + 1680*xd1**2*zd1
CASE(61)
 S=180.0D0- 840*xd1 + 840*xd1**2 - 1080*zd1 + 5040*xd1*zd1 - 5040*xd1**2*zd1 + 1080*zd1**2 - 5040*xd1*zd1**2 + 5040*xd1**2*zd1**2
CASE(36)
 S=-180.0D0+ 840*xd1 - 840*xd1**2 + 360*yd1 - 1680*xd1*yd1 + 1680*xd1**2*yd1
CASE(60)
 S=180.0D0- 840*xd1 + 840*xd1**2 - 360*yd1 + 1680*xd1*yd1 - 1680*xd1**2*yd1 - 360*zd1 + 1680*xd1*zd1 - 1680*xd1**2*zd1 &
+ 720*yd1*zd1 - 3360*xd1*yd1*zd1 + 3360*xd1**2*yd1*zd1
CASE(59)
 S=180.0D0- 840*xd1 + 840*xd1**2 - 1080*yd1 + 5040*xd1*yd1 - 5040*xd1**2*yd1 + 1080*yd1**2 - 5040*xd1*yd1**2 + 5040*xd1**2*yd1**2
CASE(35)
 S=-420.0D0+ 3360*xd1 - 7560*xd1**2 + 5040*xd1**3
CASE(58)
 S=420.0D0- 3360*xd1 + 7560*xd1**2 - 5040*xd1**3 - 840*zd1 + 6720*xd1*zd1 - 15120*xd1**2*zd1 + 10080*xd1**3*zd1
CASE(57)
 S=420.0D0- 3360*xd1 + 7560*xd1**2 - 5040*xd1**3 - 840*yd1 + 6720*xd1*yd1 - 15120*xd1**2*yd1 + 10080*xd1**3*yd1
CASE(56)
 S=840.0D0- 10080*xd1 + 37800*xd1**2 - 55440*xd1**3 + 27720*xd1**4
  END SELECT
 DLX2=S    
	END FUNCTION DLX2
 
 
REAL FUNCTION DLY2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=12
CASE(14)
 S=-12 + 24.d0*zd1
CASE(31)
 S=12 - 72.d0*zd1 + 72.d0*zd1**2
CASE(52)
 S=-12 + 144.d0*zd1 - 360*zd1**2 + 240*zd1**3
CASE(81)
 S=12 - 240*zd1 + 1080*zd1**2 - 1680*zd1**3 + 840*zd1**4
CASE(13)
 S=-60.0D0+ 120*yd1
CASE(30)
 S=60.0D0- 120*yd1 - 120*zd1 + 240*yd1*zd1
CASE(51)
 S=-60.0D0+ 120*yd1 + 360*zd1 - 720*yd1*zd1 - 360*zd1**2 + 720*yd1*zd1**2
CASE(80)
 S=60.0D0- 120*yd1 - 720*zd1 + 1440*yd1*zd1 + 1800*zd1**2 - 3600*yd1*zd1**2 - 1200*zd1**3 + 2400*yd1*zd1**3
CASE(34)
 S=180.0D0- 840*yd1 + 840*yd1**2
CASE(50)
 S=-180.0D0+ 840*yd1 - 840*yd1**2 + 360*zd1 - 1680*yd1*zd1 + 1680*yd1**2*zd1
CASE(79)
 S=180.0D0- 840*yd1 + 840*yd1**2 - 1080*zd1 + 5040*yd1*zd1 - 5040*yd1**2*zd1 + 1080*zd1**2 - 5040*yd1*zd1**2 &
+ 5040*yd1**2*zd1**2
CASE(49)
 S=-420.0D0+ 3360*yd1 - 7560*yd1**2 + 5040*yd1**3
CASE(78)
 S=420.0D0- 3360*yd1 + 7560*yd1**2 - 5040*yd1**3 - 840*zd1 + 6720*yd1*zd1 - 15120*yd1**2*zd1 + 10080*yd1**3*zd1
CASE(77)
 S=840.0D0- 10080*yd1 + 37800*yd1**2 - 55440*yd1**3 + 27720*yd1**4
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=-12 + 24.d0*xd1
CASE(27)
 S=12 - 24.d0*xd1 - 24.d0*zd1 + 48.d0*xd1*zd1
CASE(48)
 S=-12 + 24.d0*xd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 72.d0*zd1**2 + 144.d0*xd1*zd1**2
CASE(74)
 S=12 - 24.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1 + 360*zd1**2 - 720*xd1*zd1**2 - 240*zd1**3 + 480*xd1*zd1**3
CASE(26)
 S=60.0D0- 120*xd1 - 120*yd1 + 240*xd1*yd1
CASE(46)
 S=-60.0D0+ 120*xd1 + 120*yd1 - 240*xd1*yd1 + 120*zd1 - 240*xd1*zd1 - 240*yd1*zd1 + 480*xd1*yd1*zd1
CASE(73)
 S=60.0D0- 120*xd1 - 120*yd1 + 240*xd1*yd1 - 360*zd1 + 720*xd1*zd1 + 720*yd1*zd1 - 1440*xd1*yd1*zd1 + 360*zd1**2 &
- 720*xd1*zd1**2 - 720*yd1*zd1**2 + 1440*xd1*yd1*zd1**2
CASE(45)
 S=-180.0D0+ 360*xd1 + 840*yd1 - 1680*xd1*yd1 - 840*yd1**2 + 1680*xd1*yd1**2
CASE(72)
 S=180.0D0- 360*xd1 - 840*yd1 + 1680*xd1*yd1 + 840*yd1**2 - 1680*xd1*yd1**2 - 360*zd1 + 720*xd1*zd1 + 1680*yd1*zd1 &
- 3360*xd1*yd1*zd1 - 1680*yd1**2*zd1 + 3360*xd1*yd1**2*zd1
CASE(71)
 S=420.0D0- 840*xd1 - 3360*yd1 + 6720*xd1*yd1 + 7560*yd1**2 - 15120*xd1*yd1**2 - 5040*yd1**3 + 10080*xd1*yd1**3
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=12 - 72.d0*xd1 + 72.d0*xd1**2
CASE(42)
 S=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*zd1 - 144.d0*xd1*zd1 + 144.d0*xd1**2*zd1
CASE(68)
 S=12 - 72.d0*xd1 + 72.d0*xd1**2 - 72.d0*zd1 + 432.d0*xd1*zd1 - 432.d0*xd1**2*zd1 + 72.d0*zd1**2 - 432.d0*xd1*zd1**2 &
+ 432.d0*xd1**2*zd1**2
CASE(41)
 S=-60.0D0+ 360*xd1 - 360*xd1**2 + 120*yd1 - 720*xd1*yd1 + 720*xd1**2*yd1
CASE(67)
 S=60.0D0- 360*xd1 + 360*xd1**2 - 120*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 - 120*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1&
 + 240*yd1*zd1 - 1440*xd1*yd1*zd1 + 1440*xd1**2*yd1*zd1
CASE(66)
 S=180.0D0- 1080*xd1 + 1080*xd1**2 - 840*yd1 + 5040*xd1*yd1 - 5040*xd1**2*yd1 + 840*yd1**2 - 5040*xd1*yd1**2 &
+ 5040*xd1**2*yd1**2
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=-12 + 144.d0*xd1 - 360*xd1**2 + 240*xd1**3
CASE(63)
 S=12 - 144.d0*xd1 + 360*xd1**2 - 240*xd1**3 - 24.d0*zd1 + 288.d0*xd1*zd1 - 720*xd1**2*zd1 + 480*xd1**3*zd1
CASE(62)
 S=60.0D0- 720*xd1 + 1800*xd1**2 - 1200*xd1**3 - 120*yd1 + 1440*xd1*yd1 - 3600*xd1**2*yd1 + 2400*xd1**3*yd1
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=12 - 240*xd1 + 1080*xd1**2 - 1680*xd1**3 + 840*xd1**4
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY2=S    
	END FUNCTION DLY2
 
REAL FUNCTION DLZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=12
CASE(16)
 S=-60.0D0+ 120*zd1
CASE(33)
 S=180.0D0- 840*zd1 + 840*zd1**2
CASE(54)
 S=-420.0D0+ 3360*zd1 - 7560*zd1**2 + 5040*zd1**3
CASE(83)
 S=840.0D0- 10080*zd1 + 37800*zd1**2 - 55440*zd1**3 + 27720*zd1**4
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=-12 + 24.d0*yd1
CASE(32)
 S=60.0D0- 120*yd1 - 120*zd1 + 240*yd1*zd1
CASE(53)
 S=-180.0D0+ 360*yd1 + 840*zd1 - 1680*yd1*zd1 - 840*zd1**2 + 1680*yd1*zd1**2
CASE(82)
 S=420.0D0- 840*yd1 - 3360*zd1 + 6720*yd1*zd1 + 7560*zd1**2 - 15120*yd1*zd1**2 - 5040*zd1**3 + 10080*yd1*zd1**3
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=12 - 72.d0*yd1 + 72.d0*yd1**2
CASE(52)
 S=-60.0D0+ 360*yd1 - 360*yd1**2 + 120*zd1 - 720*yd1*zd1 + 720*yd1**2*zd1
CASE(81)
 S=180.0D0- 1080*yd1 + 1080*yd1**2 - 840*zd1 + 5040*yd1*zd1 - 5040*yd1**2*zd1 + 840*zd1**2 - 5040*yd1*zd1**2&
 + 5040*yd1**2*zd1**2
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=-12 + 144.d0*yd1 - 360*yd1**2 + 240*yd1**3
CASE(80)
 S=60.0D0- 720*yd1 + 1800*yd1**2 - 1200*yd1**3 - 120*zd1 + 1440*yd1*zd1 - 3600*yd1**2*zd1 + 2400*yd1**3*zd1
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=12 - 240*yd1 + 1080*yd1**2 - 1680*yd1**3 + 840*yd1**4
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=-12 + 24.d0*xd1
CASE(29)
 S=60.0D0- 120*xd1 - 120*zd1 + 240*xd1*zd1
CASE(55)
 S=-180.0D0+ 360*xd1 + 840*zd1 - 1680*xd1*zd1 - 840*zd1**2 + 1680*xd1*zd1**2
CASE(76)
 S=420.0D0- 840*xd1 - 3360*zd1 + 6720*xd1*zd1 + 7560*zd1**2 - 15120*xd1*zd1**2 - 5040*zd1**3 + 10080*xd1*zd1**3
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=12 - 24.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1
CASE(47)
 S=-60.0D0+ 120*xd1 + 120*yd1 - 240*xd1*yd1 + 120*zd1 - 240*xd1*zd1 - 240*yd1*zd1 + 480*xd1*yd1*zd1
CASE(75)
 S=180.0D0- 360*xd1 - 360*yd1 + 720*xd1*yd1 - 840*zd1 + 1680*xd1*zd1 + 1680*yd1*zd1 - 3360*xd1*yd1*zd1 + 840*zd1**2&
 - 1680*xd1*zd1**2 - 1680*yd1*zd1**2 + 3360*xd1*yd1*zd1**2
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=-12 + 24.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 72.d0*yd1**2 + 144.d0*xd1*yd1**2
CASE(74)
 S=60.0D0- 120*xd1 - 360*yd1 + 720*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 120*zd1 + 240*xd1*zd1 + 720*yd1*zd1 &
- 1440*xd1*yd1*zd1 - 720*yd1**2*zd1 + 1440*xd1*yd1**2*zd1
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=12 - 24.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 240*yd1**3 + 480*xd1*yd1**3
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=12 - 72.d0*xd1 + 72.d0*xd1**2
CASE(44)
 S=-60.0D0+ 360*xd1 - 360*xd1**2 + 120*zd1 - 720*xd1*zd1 + 720*xd1**2*zd1
CASE(70)
 S=180.0D0- 1080*xd1 + 1080*xd1**2 - 840*zd1 + 5040*xd1*zd1 - 5040*xd1**2*zd1 + 840*zd1**2 - 5040*xd1*zd1**2 &
+ 5040*xd1**2*zd1**2
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*yd1 - 144.d0*xd1*yd1 + 144.d0*xd1**2*yd1
CASE(69)
 S=60.0D0- 360*xd1 + 360*xd1**2 - 120*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 - 120*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1 &
+ 240*yd1*zd1 - 1440*xd1*yd1*zd1 + 1440*xd1**2*yd1*zd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=12 - 72.d0*xd1 + 72.d0*xd1**2 - 72.d0*yd1 + 432.d0*xd1*yd1 - 432.d0*xd1**2*yd1 + 72.d0*yd1**2 - 432.d0*xd1*yd1**2 + 432.d0*xd1**2*yd1**2
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=-12 + 144.d0*xd1 - 360*xd1**2 + 240*xd1**3
CASE(65)
 S=60.0D0- 720*xd1 + 1800*xd1**2 - 1200*xd1**3 - 120*zd1 + 1440*xd1*zd1 - 3600*xd1**2*zd1 + 2400*xd1**3*zd1
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=12 - 144.d0*xd1 + 360*xd1**2 - 240*xd1**3 - 24.d0*yd1 + 288.d0*xd1*yd1 - 720*xd1**2*yd1 + 480*xd1**3*yd1
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=12 - 240*xd1 + 1080*xd1**2 - 1680*xd1**3 + 840*xd1**4
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLZ2=S    
	END FUNCTION DLZ2
 
REAL FUNCTION DLXY(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=4
CASE(17)
 S=-4 + 8.d0*zd1
CASE(28)
 S=4 - 24.d0*zd1 + 24.d0*zd1**2
CASE(47)
 S=-4 + 48.d0*zd1 - 120*zd1**2 + 80*zd1**3
CASE(75)
 S=4 - 80*zd1 + 360*zd1**2 - 560*zd1**3 + 280*zd1**4
CASE(12)
 S=-12 + 24.d0*yd1
CASE(27)
 S=12 - 24.d0*yd1 - 24.d0*zd1 + 48.d0*yd1*zd1
CASE(48)
 S=-12 + 24.d0*yd1 + 72.d0*zd1 - 144.d0*yd1*zd1 - 72.d0*zd1**2 + 144.d0*yd1*zd1**2
CASE(74)
 S=12 - 24.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1 + 360*zd1**2 - 720*yd1*zd1**2 - 240*zd1**3 + 480*yd1*zd1**3
CASE(26)
 S=24 - 120*yd1 + 120*yd1**2
CASE(46)
 S=-24 + 120*yd1 - 120*yd1**2 + 48.d0*zd1 - 240*yd1*zd1 + 240*yd1**2*zd1
CASE(73)
 S=24 - 120*yd1 + 120*yd1**2 - 144.d0*zd1 + 720*yd1*zd1 - 720*yd1**2*zd1 + 144.d0*zd1**2 - 720*yd1*zd1**2 + 720*yd1**2*zd1**2
CASE(45)
 S=-40.0D0+ 360*yd1 - 840*yd1**2 + 560*yd1**3
CASE(72)
 S=40.0D0- 360*yd1 + 840*yd1**2 - 560*yd1**3 - 80*zd1 + 720*yd1*zd1 - 1680*yd1**2*zd1 + 1120*yd1**3*zd1
CASE(71)
 S=60.0D0- 840*yd1 + 3360*yd1**2 - 5040*yd1**3 + 2520*yd1**4
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=-12 + 24.d0*xd1
CASE(24)
 S=12 - 24.d0*xd1 - 24.d0*zd1 + 48.d0*xd1*zd1
CASE(43)
 S=-12 + 24.d0*xd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 72.d0*zd1**2 + 144.d0*xd1*zd1**2
CASE(69)
 S=12 - 24.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1 + 360*zd1**2 - 720*xd1*zd1**2 - 240*zd1**3 + 480*xd1*zd1**3
CASE(23)
 S=36 - 72.d0*xd1 - 72.d0*yd1 + 144.d0*xd1*yd1
CASE(42)
 S=-36 + 72.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1
CASE(68)
 S=36 - 72.d0*xd1 - 72.d0*yd1 + 144.d0*xd1*yd1 - 216.d0*zd1 + 432.d0*xd1*zd1 + 432.d0*yd1*zd1 - 864.d0*xd1*yd1*zd1 + 216.d0*zd1**2 &
- 432.d0*xd1*zd1**2 - 432.d0*yd1*zd1**2 + 864.d0*xd1*yd1*zd1**2
CASE(41)
 S=-72 + 144.d0*xd1 + 360*yd1 - 720*xd1*yd1 - 360*yd1**2 + 720*xd1*yd1**2
CASE(67)
 S=72 - 144.d0*xd1 - 360*yd1 + 720*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 144.d0*zd1 + 288.d0*xd1*zd1 + 720*yd1*zd1 &
- 1440*xd1*yd1*zd1 - 720*yd1**2*zd1 + 1440*xd1*yd1**2*zd1
CASE(66)
 S=120.0D0- 240*xd1 - 1080*yd1 + 2160*xd1*yd1 + 2520*yd1**2 - 5040*xd1*yd1**2 - 1680*yd1**3 + 3360*xd1*yd1**3
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=24 - 120*xd1 + 120*xd1**2
CASE(40)
 S=-24 + 120*xd1 - 120*xd1**2 + 48.d0*zd1 - 240*xd1*zd1 + 240*xd1**2*zd1
CASE(64)
 S=24 - 120*xd1 + 120*xd1**2 - 144.d0*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1 + 144.d0*zd1**2 - 720*xd1*zd1**2 &
+ 720*xd1**2*zd1**2
CASE(38)
 S=-72 + 360*xd1 - 360*xd1**2 + 144.d0*yd1 - 720*xd1*yd1 + 720*xd1**2*yd1
CASE(63)
 S=72 - 360*xd1 + 360*xd1**2 - 144.d0*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 - 144.d0*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1 &
+ 288.d0*yd1*zd1 - 1440*xd1*yd1*zd1 + 1440*xd1**2*yd1*zd1
CASE(62)
 S=144 - 720*xd1 + 720*xd1**2 - 720*yd1 + 3600*xd1*yd1 - 3600*xd1**2*yd1 + 720*yd1**2 - 3600*xd1*yd1**2 &
+ 3600*xd1**2*yd1**2
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=-40.0D0+ 360*xd1 - 840*xd1**2 + 560*xd1**3
CASE(60)
 S=40.0D0- 360*xd1 + 840*xd1**2 - 560*xd1**3 - 80*zd1 + 720*xd1*zd1 - 1680*xd1**2*zd1 + 1120*xd1**3*zd1
CASE(59)
 S=120.0D0- 1080*xd1 + 2520*xd1**2 - 1680*xd1**3 - 240*yd1 + 2160*xd1*yd1 - 5040*xd1**2*yd1 + 3360*xd1**3*yd1
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=60.0D0- 840*xd1 + 3360*xd1**2 - 5040*xd1**3 + 2520*xd1**4
CASE(56)
 S=0.0D0
 
  END SELECT
 DLXY=S    
	END FUNCTION DLXY
 
REAL FUNCTION DLYZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=4
CASE(15)
 S=-12 + 24.d0*zd1
CASE(32)
 S=24 - 120*zd1 + 120*zd1**2
CASE(53)
 S=-40.0D0+ 360*zd1 - 840*zd1**2 + 560*zd1**3
CASE(82)
 S=60.0D0- 840*zd1 + 3360*zd1**2 - 5040*zd1**3 + 2520*zd1**4
CASE(7)
 S=0.0D0
CASE(14)
 S=-12 + 24.d0*yd1
CASE(31)
 S=36 - 72.d0*yd1 - 72.d0*zd1 + 144.d0*yd1*zd1
CASE(52)
 S=-72 + 144.d0*yd1 + 360*zd1 - 720*yd1*zd1 - 360*zd1**2 + 720*yd1*zd1**2
CASE(81)
 S=120.0D0- 240*yd1 - 1080*zd1 + 2160*yd1*zd1 + 2520*zd1**2 - 5040*yd1*zd1**2 - 1680*zd1**3 + 3360*yd1*zd1**3
CASE(13)
 S=0.0D0
CASE(30)
 S=24 - 120*yd1 + 120*yd1**2
CASE(51)
 S=-72 + 360*yd1 - 360*yd1**2 + 144.d0*zd1 - 720*yd1*zd1 + 720*yd1**2*zd1
CASE(80)
 S=144 - 720*yd1 + 720*yd1**2 - 720*zd1 + 3600*yd1*zd1 - 3600*yd1**2*zd1 + 720*zd1**2 - 3600*yd1*zd1**2 + 3600*yd1**2*zd1**2
CASE(34)
 S=0.0D0
CASE(50)
 S=-40.0D0+ 360*yd1 - 840*yd1**2 + 560*yd1**3
CASE(79)
 S=120.0D0- 1080*yd1 + 2520*yd1**2 - 1680*yd1**3 - 240*zd1 + 2160*yd1*zd1 - 5040*yd1**2*zd1 + 3360*yd1**3*zd1
CASE(49)
 S=0.0D0
CASE(78)
 S=60.0D0- 840*yd1 + 3360*yd1**2 - 5040*yd1**3 + 2520*yd1**4
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=-4 + 8.d0*xd1
CASE(28)
 S=12 - 24.d0*xd1 - 24.d0*zd1 + 48.d0*xd1*zd1
CASE(47)
 S=-24 + 48.d0*xd1 + 120*zd1 - 240*xd1*zd1 - 120*zd1**2 + 240*xd1*zd1**2
CASE(75)
 S=40.0D0- 80*xd1 - 360*zd1 + 720*xd1*zd1 + 840*zd1**2 - 1680*xd1*zd1**2 - 560*zd1**3 + 1120*xd1*zd1**3
CASE(12)
 S=0.0D0
CASE(27)
 S=12 - 24.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1
CASE(48)
 S=-36 + 72.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1
CASE(74)
 S=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1 - 360*zd1 + 720*xd1*zd1 + 720*yd1*zd1 - 1440*xd1*yd1*zd1 + 360*zd1**2 &
- 720*xd1*zd1**2 - 720*yd1*zd1**2 + 1440*xd1*yd1*zd1**2
CASE(26)
 S=0.0D0
CASE(46)
 S=-24 + 48.d0*xd1 + 120*yd1 - 240*xd1*yd1 - 120*yd1**2 + 240*xd1*yd1**2
CASE(73)
 S=72 - 144.d0*xd1 - 360*yd1 + 720*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 144.d0*zd1 + 288.d0*xd1*zd1 + 720*yd1*zd1&
 - 1440*xd1*yd1*zd1 - 720*yd1**2*zd1 + 1440*xd1*yd1**2*zd1
CASE(45)
 S=0.0D0
CASE(72)
 S=40.0D0- 80*xd1 - 360*yd1 + 720*xd1*yd1 + 840*yd1**2 - 1680*xd1*yd1**2 - 560*yd1**3 + 1120*xd1*yd1**3
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=4 - 24.d0*xd1 + 24.d0*xd1**2
CASE(43)
 S=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*zd1 - 144.d0*xd1*zd1 + 144.d0*xd1**2*zd1
CASE(69)
 S=24 - 144.d0*xd1 + 144.d0*xd1**2 - 120*zd1 + 720*xd1*zd1 - 720*xd1**2*zd1 + 120*zd1**2 - 720*xd1*zd1**2 + 720*xd1**2*zd1**2
CASE(23)
 S=0.0D0
CASE(42)
 S=-12 + 72.d0*xd1 - 72.d0*xd1**2 + 24.d0*yd1 - 144.d0*xd1*yd1 + 144.d0*xd1**2*yd1
CASE(68)
 S=36 - 216.d0*xd1 + 216.d0*xd1**2 - 72.d0*yd1 + 432.d0*xd1*yd1 - 432.d0*xd1**2*yd1 - 72.d0*zd1 + 432.d0*xd1*zd1 - 432.d0*xd1**2*zd1 &
+ 144.d0*yd1*zd1 - 864.d0*xd1*yd1*zd1 + 864.d0*xd1**2*yd1*zd1
CASE(41)
 S=0.0D0
CASE(67)
 S=24 - 144.d0*xd1 + 144.d0*xd1**2 - 120*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 + 120*yd1**2 - 720*xd1*yd1**2 + 720*xd1**2*yd1**2
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=-4 + 48.d0*xd1 - 120*xd1**2 + 80*xd1**3
CASE(64)
 S=12 - 144.d0*xd1 + 360*xd1**2 - 240*xd1**3 - 24.d0*zd1 + 288.d0*xd1*zd1 - 720*xd1**2*zd1 + 480*xd1**3*zd1
CASE(38)
 S=0.0D0
CASE(63)
 S=12 - 144.d0*xd1 + 360*xd1**2 - 240*xd1**3 - 24.d0*yd1 + 288.d0*xd1*yd1 - 720*xd1**2*yd1 + 480*xd1**3*yd1
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=4 - 80*xd1 + 360*xd1**2 - 560*xd1**3 + 280*xd1**4
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLYZ=S    
	END FUNCTION DLYZ
 
REAL FUNCTION DLXZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=4
CASE(19)
 S=-12 + 24.d0*zd1
CASE(29)
 S=24 - 120*zd1 + 120*zd1**2
CASE(55)
 S=-40.0D0+ 360*zd1 - 840*zd1**2 + 560*zd1**3
CASE(76)
 S=60.0D0- 840*zd1 + 3360*zd1**2 - 5040*zd1**3 + 2520*zd1**4
CASE(5)
 S=0.0D0
CASE(17)
 S=-4 + 8.d0*yd1
CASE(28)
 S=12 - 24.d0*yd1 - 24.d0*zd1 + 48.d0*yd1*zd1
CASE(47)
 S=-24 + 48.d0*yd1 + 120*zd1 - 240*yd1*zd1 - 120*zd1**2 + 240*yd1*zd1**2
CASE(75)
 S=40.0D0- 80*yd1 - 360*zd1 + 720*yd1*zd1 + 840*zd1**2 - 1680*yd1*zd1**2 - 560*zd1**3 + 1120*yd1*zd1**3
CASE(12)
 S=0.0D0
CASE(27)
 S=4 - 24.d0*yd1 + 24.d0*yd1**2
CASE(48)
 S=-12 + 72.d0*yd1 - 72.d0*yd1**2 + 24.d0*zd1 - 144.d0*yd1*zd1 + 144.d0*yd1**2*zd1
CASE(74)
 S=24 - 144.d0*yd1 + 144.d0*yd1**2 - 120*zd1 + 720*yd1*zd1 - 720*yd1**2*zd1 + 120*zd1**2 - 720*yd1*zd1**2 + 720*yd1**2*zd1**2
CASE(26)
 S=0.0D0
CASE(46)
 S=-4 + 48.d0*yd1 - 120*yd1**2 + 80*yd1**3
CASE(73)
 S=12 - 144.d0*yd1 + 360*yd1**2 - 240*yd1**3 - 24.d0*zd1 + 288.d0*yd1*zd1 - 720*yd1**2*zd1 + 480*yd1**3*zd1
CASE(45)
 S=0.0D0
CASE(72)
 S=4 - 80*yd1 + 360*yd1**2 - 560*yd1**3 + 280*yd1**4
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=-12 + 24.d0*xd1
CASE(25)
 S=36 - 72.d0*xd1 - 72.d0*zd1 + 144.d0*xd1*zd1
CASE(44)
 S=-72 + 144.d0*xd1 + 360*zd1 - 720*xd1*zd1 - 360*zd1**2 + 720*xd1*zd1**2
CASE(70)
 S=120.0D0- 240*xd1 - 1080*zd1 + 2160*xd1*zd1 + 2520*zd1**2 - 5040*xd1*zd1**2 - 1680*zd1**3 + 3360*xd1*zd1**3
CASE(11)
 S=0.0D0
CASE(24)
 S=12 - 24.d0*xd1 - 24.d0*yd1 + 48.d0*xd1*yd1
CASE(43)
 S=-36 + 72.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 + 72.d0*zd1 - 144.d0*xd1*zd1 - 144.d0*yd1*zd1 + 288.d0*xd1*yd1*zd1
CASE(69)
 S=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1 - 360*zd1 + 720*xd1*zd1 + 720*yd1*zd1 - 1440*xd1*yd1*zd1 + 360*zd1**2 &
- 720*xd1*zd1**2 - 720*yd1*zd1**2 + 1440*xd1*yd1*zd1**2
CASE(23)
 S=0.0D0
CASE(42)
 S=-12 + 24.d0*xd1 + 72.d0*yd1 - 144.d0*xd1*yd1 - 72.d0*yd1**2 + 144.d0*xd1*yd1**2
CASE(68)
 S=36 - 72.d0*xd1 - 216.d0*yd1 + 432.d0*xd1*yd1 + 216.d0*yd1**2 - 432.d0*xd1*yd1**2 - 72.d0*zd1 + 144.d0*xd1*zd1 + 432.d0*yd1*zd1 &
- 864.d0*xd1*yd1*zd1 - 432.d0*yd1**2*zd1 + 864.d0*xd1*yd1**2*zd1
CASE(41)
 S=0.0D0
CASE(67)
 S=12 - 24.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1 + 360*yd1**2 - 720*xd1*yd1**2 - 240*yd1**3 + 480*xd1*yd1**3
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=24 - 120*xd1 + 120*xd1**2
CASE(39)
 S=-72 + 360*xd1 - 360*xd1**2 + 144.d0*zd1 - 720*xd1*zd1 + 720*xd1**2*zd1
CASE(65)
 S=144 - 720*xd1 + 720*xd1**2 - 720*zd1 + 3600*xd1*zd1 - 3600*xd1**2*zd1 + 720*zd1**2 - 3600*xd1*zd1**2 &
+ 3600*xd1**2*zd1**2
CASE(21)
 S=0.0D0
CASE(40)
 S=-24 + 120*xd1 - 120*xd1**2 + 48.d0*yd1 - 240*xd1*yd1 + 240*xd1**2*yd1
CASE(64)
 S=72 - 360*xd1 + 360*xd1**2 - 144.d0*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 - 144.d0*zd1 + 720*xd1*zd1 &
- 720*xd1**2*zd1 + 288.d0*yd1*zd1 - 1440*xd1*yd1*zd1 + 1440*xd1**2*yd1*zd1
CASE(38)
 S=0.0D0
CASE(63)
 S=24 - 120*xd1 + 120*xd1**2 - 144.d0*yd1 + 720*xd1*yd1 - 720*xd1**2*yd1 + 144.d0*yd1**2 - 720*xd1*yd1**2 &
+ 720*xd1**2*yd1**2
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=-40.0D0+ 360*xd1 - 840*xd1**2 + 560*xd1**3
CASE(61)
 S=120.0D0- 1080*xd1 + 2520*xd1**2 - 1680*xd1**3 - 240*zd1 + 2160*xd1*zd1 - 5040*xd1**2*zd1 + 3360*xd1**3*zd1
CASE(36)
 S=0.0D0
CASE(60)
 S=40.0D0- 360*xd1 + 840*xd1**2 - 560*xd1**3 - 80*yd1 + 720*xd1*yd1 - 1680*xd1**2*yd1 + 1120*xd1**3*yd1
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=60.0D0- 840*xd1 + 3360*xd1**2 - 5040*xd1**3 + 2520*xd1**4
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXZ=S    
	END FUNCTION DLXZ
 
 
REAL FUNCTION DLX3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=120
CASE(22)
 S=-120.0D0+ 240*zd1
CASE(39)
 S=120.0D0- 720*zd1 + 720*zd1**2
CASE(65)
 S=-120.0D0+ 1440*zd1 - 3600*zd1**2 + 2400*zd1**3
CASE(21)
 S=-120.0D0+ 240*yd1
CASE(40)
 S=120.0D0- 240*yd1 - 240*zd1 + 480*yd1*zd1
CASE(64)
 S=-120.0D0+ 240*yd1 + 720*zd1 - 1440*yd1*zd1 - 720*zd1**2 + 1440*yd1*zd1**2
CASE(38)
 S=120.0D0- 720*yd1 + 720*yd1**2
CASE(63)
 S=-120.0D0+ 720*yd1 - 720*yd1**2 + 240*zd1 - 1440*yd1*zd1 + 1440*yd1**2*zd1
CASE(62)
 S=-120.0D0+ 1440*yd1 - 3600*yd1**2 + 2400*yd1**3
CASE(20)
 S=-840.0D0+ 1680*xd1
CASE(37)
 S=840.0D0- 1680*xd1 - 1680*zd1 + 3360*xd1*zd1
CASE(61)
 S=-840.0D0+ 1680*xd1 + 5040*zd1 - 10080*xd1*zd1 - 5040*zd1**2 + 10080*xd1*zd1**2
CASE(36)
 S=840.0D0- 1680*xd1 - 1680*yd1 + 3360*xd1*yd1
CASE(60)
 S=-840.0D0+ 1680*xd1 + 1680*yd1 - 3360*xd1*yd1 + 1680*zd1 - 3360*xd1*zd1 - 3360*yd1*zd1 + 6720*xd1*yd1*zd1
CASE(59)
 S=-840.0D0+ 1680*xd1 + 5040*yd1 - 10080*xd1*yd1 - 5040*yd1**2 + 10080*xd1*yd1**2
CASE(35)
 S=3360.0D0- 15120*xd1 + 15120*xd1**2
CASE(58)
 S=-3360.0D0+ 15120*xd1 - 15120*xd1**2 + 6720*zd1 - 30240*xd1*zd1 + 30240*xd1**2*zd1
CASE(57)
 S=-3360.0D0+ 15120*xd1 - 15120*xd1**2 + 6720*yd1 - 30240*xd1*yd1 + 30240*xd1**2*yd1
CASE(56)
 S=-10080.0D0+ 75600*xd1 - 166320*xd1**2 + 110880*xd1**3
  END SELECT
 DLX3=S    
	END FUNCTION DLX3
 
REAL FUNCTION DLX2Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0

SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=24
CASE(24)
 S=-24 + 48.d0*zd1
CASE(43)
 S=24 - 144.d0*zd1 + 144.d0*zd1**2
CASE(69)
 S=-24 + 288.d0*zd1 - 720*zd1**2 + 480*zd1**3
CASE(23)
 S=-72 + 144.d0*yd1
CASE(42)
 S=72 - 144.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1
CASE(68)
 S=-72 + 144.d0*yd1 + 432.d0*zd1 - 864.d0*yd1*zd1 - 432.d0*zd1**2 + 864.d0*yd1*zd1**2
CASE(41)
 S=144 - 720*yd1 + 720*yd1**2
CASE(67)
 S=-144 + 720*yd1 - 720*yd1**2 + 288.d0*zd1 - 1440*yd1*zd1 + 1440*yd1**2*zd1
CASE(66)
 S=-240.0D0+ 2160*yd1 - 5040*yd1**2 + 3360*yd1**3
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=-120.0D0+ 240*xd1
CASE(40)
 S=120.0D0- 240*xd1 - 240*zd1 + 480*xd1*zd1
CASE(64)
 S=-120.0D0+ 240*xd1 + 720*zd1 - 1440*xd1*zd1 - 720*zd1**2 + 1440*xd1*zd1**2
CASE(38)
 S=360.0D0- 720*xd1 - 720*yd1 + 1440*xd1*yd1
CASE(63)
 S=-360.0D0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
CASE(62)
 S=-720.0D0+ 1440*xd1 + 3600*yd1 - 7200*xd1*yd1 - 3600*yd1**2 + 7200*xd1*yd1**2
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=360.0D0- 1680*xd1 + 1680*xd1**2
CASE(60)
 S=-360.0D0+ 1680*xd1 - 1680*xd1**2 + 720*zd1 - 3360*xd1*zd1 + 3360*xd1**2*zd1
CASE(59)
 S=-1080.0D0+ 5040*xd1 - 5040*xd1**2 + 2160*yd1 - 10080*xd1*yd1 + 10080*xd1**2*yd1
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=-840.0D0+ 6720*xd1 - 15120*xd1**2 + 10080*xd1**3
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Y=S    
	END FUNCTION DLX2Y
 
REAL FUNCTION DLXY2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=24
CASE(27)
 S=-24 + 48.d0*zd1
CASE(48)
 S=24 - 144.d0*zd1 + 144.d0*zd1**2
CASE(74)
 S=-24 + 288.d0*zd1 - 720*zd1**2 + 480*zd1**3
CASE(26)
 S=-120.0D0+ 240*yd1
CASE(46)
 S=120.0D0- 240*yd1 - 240*zd1 + 480*yd1*zd1
CASE(73)
 S=-120.0D0+ 240*yd1 + 720*zd1 - 1440*yd1*zd1 - 720*zd1**2 + 1440*yd1*zd1**2
CASE(45)
 S=360.0D0- 1680*yd1 + 1680*yd1**2
CASE(72)
 S=-360.0D0+ 1680*yd1 - 1680*yd1**2 + 720*zd1 - 3360*yd1*zd1 + 3360*yd1**2*zd1
CASE(71)
 S=-840.0D0+ 6720*yd1 - 15120*yd1**2 + 10080*yd1**3
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=-72 + 144.d0*xd1
CASE(42)
 S=72 - 144.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1
CASE(68)
 S=-72 + 144.d0*xd1 + 432.d0*zd1 - 864.d0*xd1*zd1 - 432.d0*zd1**2 + 864.d0*xd1*zd1**2
CASE(41)
 S=360.0D0- 720*xd1 - 720*yd1 + 1440*xd1*yd1
CASE(67)
 S=-360.0D0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
CASE(66)
 S=-1080.0D0+ 2160*xd1 + 5040*yd1 - 10080*xd1*yd1 - 5040*yd1**2 + 10080*xd1*yd1**2
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=144 - 720*xd1 + 720*xd1**2
CASE(63)
 S=-144 + 720*xd1 - 720*xd1**2 + 288.d0*zd1 - 1440*xd1*zd1 + 1440*xd1**2*zd1
CASE(62)
 S=-720.0D0+ 3600*xd1 - 3600*xd1**2 + 1440*yd1 - 7200*xd1*yd1 + 7200*xd1**2*yd1
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=-240.0D0+ 2160*xd1 - 5040*xd1**2 + 3360*xd1**3
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY2=S    
	END FUNCTION DLXY2
 
 
REAL FUNCTION DLY3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=120
CASE(30)
 S=-120.0D0+ 240*zd1
CASE(51)
 S=120.0D0- 720*zd1 + 720*zd1**2
CASE(80)
 S=-120.0D0+ 1440*zd1 - 3600*zd1**2 + 2400*zd1**3
CASE(34)
 S=-840.0D0+ 1680*yd1
CASE(50)
 S=840.0D0- 1680*yd1 - 1680*zd1 + 3360*yd1*zd1
CASE(79)
 S=-840.0D0+ 1680*yd1 + 5040*zd1 - 10080*yd1*zd1 - 5040*zd1**2 + 10080*yd1*zd1**2
CASE(49)
 S=3360.0D0- 15120*yd1 + 15120*yd1**2
CASE(78)
 S=-3360.0D0+ 15120*yd1 - 15120*yd1**2 + 6720*zd1 - 30240*yd1*zd1 + 30240*yd1**2*zd1
CASE(77)
 S=-10080.0D0+ 75600*yd1 - 166320*yd1**2 + 110880*yd1**3
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=-120.0D0+ 240*xd1
CASE(46)
 S=120.0D0- 240*xd1 - 240*zd1 + 480*xd1*zd1
CASE(73)
 S=-120.0D0+ 240*xd1 + 720*zd1 - 1440*xd1*zd1 - 720*zd1**2 + 1440*xd1*zd1**2
CASE(45)
 S=840.0D0- 1680*xd1 - 1680*yd1 + 3360*xd1*yd1
CASE(72)
 S=-840.0D0+ 1680*xd1 + 1680*yd1 - 3360*xd1*yd1 + 1680*zd1 - 3360*xd1*zd1 - 3360*yd1*zd1 + 6720*xd1*yd1*zd1
CASE(71)
 S=-3360.0D0+ 6720*xd1 + 15120*yd1 - 30240*xd1*yd1 - 15120*yd1**2 + 30240*xd1*yd1**2
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=120.0D0- 720*xd1 + 720*xd1**2
CASE(67)
 S=-120.0D0+ 720*xd1 - 720*xd1**2 + 240*zd1 - 1440*xd1*zd1 + 1440*xd1**2*zd1
CASE(66)
 S=-840.0D0+ 5040*xd1 - 5040*xd1**2 + 1680*yd1 - 10080*xd1*yd1 + 10080*xd1**2*yd1
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=-120.0D0+ 1440*xd1 - 3600*xd1**2 + 2400*xd1**3
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY3=S    
	END FUNCTION DLY3
 
 
REAL FUNCTION DLX2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=24
CASE(25)
 S=-72 + 144.d0*zd1
CASE(44)
 S=144 - 720*zd1 + 720*zd1**2
CASE(70)
 S=-240.0D0+ 2160*zd1 - 5040*zd1**2 + 3360*zd1**3
CASE(11)
 S=0.0D0
CASE(24)
 S=-24 + 48.d0*yd1
CASE(43)
 S=72 - 144.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1
CASE(69)
 S=-144 + 288.d0*yd1 + 720*zd1 - 1440*yd1*zd1 - 720*zd1**2 + 1440*yd1*zd1**2
CASE(23)
 S=0.0D0
CASE(42)
 S=24 - 144.d0*yd1 + 144.d0*yd1**2
CASE(68)
 S=-72 + 432.d0*yd1 - 432.d0*yd1**2 + 144.d0*zd1 - 864.d0*yd1*zd1 + 864.d0*yd1**2*zd1
CASE(41)
 S=0.0D0
CASE(67)
 S=-24 + 288.d0*yd1 - 720*yd1**2 + 480*yd1**3
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=-120.0D0+ 240*xd1
CASE(39)
 S=360.0D0- 720*xd1 - 720*zd1 + 1440*xd1*zd1
CASE(65)
 S=-720.0D0+ 1440*xd1 + 3600*zd1 - 7200*xd1*zd1 - 3600*zd1**2 + 7200*xd1*zd1**2
CASE(21)
 S=0.0D0
CASE(40)
 S=120.0D0- 240*xd1 - 240*yd1 + 480*xd1*yd1
CASE(64)
 S=-360.0D0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
CASE(38)
 S=0.0D0
CASE(63)
 S=-120.0D0+ 240*xd1 + 720*yd1 - 1440*xd1*yd1 - 720*yd1**2 + 1440*xd1*yd1**2
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=360.0D0- 1680*xd1 + 1680*xd1**2
CASE(61)
 S=-1080.0D0+ 5040*xd1 - 5040*xd1**2 + 2160*zd1 - 10080*xd1*zd1 + 10080*xd1**2*zd1
CASE(36)
 S=0.0D0
CASE(60)
 S=-360.0D0+ 1680*xd1 - 1680*xd1**2 + 720*yd1 - 3360*xd1*yd1 + 3360*xd1**2*yd1
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=-840.0D0+ 6720*xd1 - 15120*xd1**2 + 10080*xd1**3
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Z=S    
	END FUNCTION DLX2Z
 
REAL FUNCTION DLXYZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=8
CASE(28)
 S=-24 + 48.d0*zd1
CASE(47)
 S=48 - 240*zd1 + 240*zd1**2
CASE(75)
 S=-80.0D0+ 720*zd1 - 1680*zd1**2 + 1120*zd1**3
CASE(12)
 S=0.0D0
CASE(27)
 S=-24 + 48.d0*yd1
CASE(48)
 S=72 - 144.d0*yd1 - 144.d0*zd1 + 288.d0*yd1*zd1
CASE(74)
 S=-144 + 288.d0*yd1 + 720*zd1 - 1440*yd1*zd1 - 720*zd1**2 + 1440*yd1*zd1**2
CASE(26)
 S=0.0D0
CASE(46)
 S=48 - 240*yd1 + 240*yd1**2
CASE(73)
 S=-144 + 720*yd1 - 720*yd1**2 + 288.d0*zd1 - 1440*yd1*zd1 + 1440*yd1**2*zd1
CASE(45)
 S=0.0D0
CASE(72)
 S=-80.0D0+ 720*yd1 - 1680*yd1**2 + 1120*yd1**3
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=-24 + 48.d0*xd1
CASE(43)
 S=72 - 144.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1
CASE(69)
 S=-144 + 288.d0*xd1 + 720*zd1 - 1440*xd1*zd1 - 720*zd1**2 + 1440*xd1*zd1**2
CASE(23)
 S=0.0D0
CASE(42)
 S=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1
CASE(68)
 S=-216 + 432.d0*xd1 + 432.d0*yd1 - 864.d0*xd1*yd1 + 432.d0*zd1 - 864.d0*xd1*zd1 - 864.d0*yd1*zd1 + 1728.d0*xd1*yd1*zd1
CASE(41)
 S=0.0D0
CASE(67)
 S=-144 + 288.d0*xd1 + 720*yd1 - 1440*xd1*yd1 - 720*yd1**2 + 1440*xd1*yd1**2
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=48 - 240*xd1 + 240*xd1**2
CASE(64)
 S=-144 + 720*xd1 - 720*xd1**2 + 288.d0*zd1 - 1440*xd1*zd1 + 1440*xd1**2*zd1
CASE(38)
 S=0.0D0
CASE(63)
 S=-144 + 720*xd1 - 720*xd1**2 + 288.d0*yd1 - 1440*xd1*yd1 + 1440*xd1**2*yd1
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=-80.0D0+ 720*xd1 - 1680*xd1**2 + 1120*xd1**3
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXYZ=S    
	END FUNCTION DLXYZ
 
 
REAL FUNCTION DLY2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=24
CASE(31)
 S=-72 + 144.d0*zd1
CASE(52)
 S=144 - 720*zd1 + 720*zd1**2
CASE(81)
 S=-240.0D0+ 2160*zd1 - 5040*zd1**2 + 3360*zd1**3
CASE(13)
 S=0.0D0
CASE(30)
 S=-120.0D0+ 240*yd1
CASE(51)
 S=360.0D0- 720*yd1 - 720*zd1 + 1440*yd1*zd1
CASE(80)
 S=-720.0D0+ 1440*yd1 + 3600*zd1 - 7200*yd1*zd1 - 3600*zd1**2 + 7200*yd1*zd1**2
CASE(34)
 S=0.0D0
CASE(50)
 S=360.0D0- 1680*yd1 + 1680*yd1**2
CASE(79)
 S=-1080.0D0+ 5040*yd1 - 5040*yd1**2 + 2160*zd1 - 10080*yd1*zd1 + 10080*yd1**2*zd1
CASE(49)
 S=0.0D0
CASE(78)
 S=-840.0D0+ 6720*yd1 - 15120*yd1**2 + 10080*yd1**3
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=-24 + 48.d0*xd1
CASE(48)
 S=72 - 144.d0*xd1 - 144.d0*zd1 + 288.d0*xd1*zd1
CASE(74)
 S=-144 + 288.d0*xd1 + 720*zd1 - 1440*xd1*zd1 - 720*zd1**2 + 1440*xd1*zd1**2
CASE(26)
 S=0.0D0
CASE(46)
 S=120.0D0- 240*xd1 - 240*yd1 + 480*xd1*yd1
CASE(73)
 S=-360.0D0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
CASE(45)
 S=0.0D0
CASE(72)
 S=-360.0D0+ 720*xd1 + 1680*yd1 - 3360*xd1*yd1 - 1680*yd1**2 + 3360*xd1*yd1**2
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=24 - 144.d0*xd1 + 144.d0*xd1**2
CASE(68)
 S=-72 + 432.d0*xd1 - 432.d0*xd1**2 + 144.d0*zd1 - 864.d0*xd1*zd1 + 864.d0*xd1**2*zd1
CASE(41)
 S=0.0D0
CASE(67)
 S=-120.0D0+ 720*xd1 - 720*xd1**2 + 240*yd1 - 1440*xd1*yd1 + 1440*xd1**2*yd1
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=-24 + 288.d0*xd1 - 720*xd1**2 + 480*xd1**3
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY2Z=S    
	END FUNCTION DLY2Z
 
 
REAL FUNCTION DLXZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=24
CASE(29)
 S=-120.0D0+ 240*zd1
CASE(55)
 S=360.0D0- 1680*zd1 + 1680*zd1**2
CASE(76)
 S=-840.0D0+ 6720*zd1 - 15120*zd1**2 + 10080*zd1**3
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=-24 + 48.d0*yd1
CASE(47)
 S=120.0D0- 240*yd1 - 240*zd1 + 480*yd1*zd1
CASE(75)
 S=-360.0D0+ 720*yd1 + 1680*zd1 - 3360*yd1*zd1 - 1680*zd1**2 + 3360*yd1*zd1**2
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=24 - 144.d0*yd1 + 144.d0*yd1**2
CASE(74)
 S=-120.0D0+ 720*yd1 - 720*yd1**2 + 240*zd1 - 1440*yd1*zd1 + 1440*yd1**2*zd1
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=-24 + 288.d0*yd1 - 720*yd1**2 + 480*yd1**3
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=-72 + 144.d0*xd1
CASE(44)
 S=360.0D0- 720*xd1 - 720*zd1 + 1440*xd1*zd1
CASE(70)
 S=-1080.0D0+ 2160*xd1 + 5040*zd1 - 10080*xd1*zd1 - 5040*zd1**2 + 10080*xd1*zd1**2
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1
CASE(69)
 S=-360.0D0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=-72 + 144.d0*xd1 + 432.d0*yd1 - 864.d0*xd1*yd1 - 432.d0*yd1**2 + 864.d0*xd1*yd1**2
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=144 - 720*xd1 + 720*xd1**2
CASE(65)
 S=-720.0D0+ 3600*xd1 - 3600*xd1**2 + 1440*zd1 - 7200*xd1*zd1 + 7200*xd1**2*zd1
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=-144 + 720*xd1 - 720*xd1**2 + 288.d0*yd1 - 1440*xd1*yd1 + 1440*xd1**2*yd1
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=-240.0D0+ 2160*xd1 - 5040*xd1**2 + 3360*xd1**3
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXZ2=S    
	END FUNCTION DLXZ2
 
 
REAL FUNCTION DLYZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=24
CASE(32)
 S=-120.0D0+ 240*zd1
CASE(53)
 S=360.0D0- 1680*zd1 + 1680*zd1**2
CASE(82)
 S=-840.0D0+ 6720*zd1 - 15120*zd1**2 + 10080*zd1**3
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=-72 + 144.d0*yd1
CASE(52)
 S=360.0D0- 720*yd1 - 720*zd1 + 1440*yd1*zd1
CASE(81)
 S=-1080.0D0+ 2160*yd1 + 5040*zd1 - 10080*yd1*zd1 - 5040*zd1**2 + 10080*yd1*zd1**2
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=144 - 720*yd1 + 720*yd1**2
CASE(80)
 S=-720.0D0+ 3600*yd1 - 3600*yd1**2 + 1440*zd1 - 7200*yd1*zd1 + 7200*yd1**2*zd1
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=-240.0D0+ 2160*yd1 - 5040*yd1**2 + 3360*yd1**3
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=-24 + 48.d0*xd1
CASE(47)
 S=120.0D0- 240*xd1 - 240*zd1 + 480*xd1*zd1
CASE(75)
 S=-360.0D0+ 720*xd1 + 1680*zd1 - 3360*xd1*zd1 - 1680*zd1**2 + 3360*xd1*zd1**2
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=72 - 144.d0*xd1 - 144.d0*yd1 + 288.d0*xd1*yd1
CASE(74)
 S=-360.0D0+ 720*xd1 + 720*yd1 - 1440*xd1*yd1 + 720*zd1 - 1440*xd1*zd1 - 1440*yd1*zd1 + 2880*xd1*yd1*zd1
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=-144 + 288.d0*xd1 + 720*yd1 - 1440*xd1*yd1 - 720*yd1**2 + 1440*xd1*yd1**2
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=24 - 144.d0*xd1 + 144.d0*xd1**2
CASE(69)
 S=-120.0D0+ 720*xd1 - 720*xd1**2 + 240*zd1 - 1440*xd1*zd1 + 1440*xd1**2*zd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=-72 + 432.d0*xd1 - 432.d0*xd1**2 + 144.d0*yd1 - 864.d0*xd1*yd1 + 864.d0*xd1**2*yd1
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=-24 + 288.d0*xd1 - 720*xd1**2 + 480*xd1**3
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLYZ2=S    
	END FUNCTION DLYZ2
 
 
REAL FUNCTION DLZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=120
CASE(33)
 S=-840.0D0+ 1680*zd1
CASE(54)
 S=3360.0D0- 15120*zd1 + 15120*zd1**2
CASE(83)
 S=-10080.0D0+ 75600*zd1 - 166320*zd1**2 + 110880*zd1**3
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=-120.0D0+ 240*yd1
CASE(53)
 S=840.0D0- 1680*yd1 - 1680*zd1 + 3360*yd1*zd1
CASE(82)
 S=-3360.0D0+ 6720*yd1 + 15120*zd1 - 30240*yd1*zd1 - 15120*zd1**2 + 30240*yd1*zd1**2
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=120.0D0- 720*yd1 + 720*yd1**2
CASE(81)
 S=-840.0D0+ 5040*yd1 - 5040*yd1**2 + 1680*zd1 - 10080*yd1*zd1 + 10080*yd1**2*zd1
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=-120.0D0+ 1440*yd1 - 3600*yd1**2 + 2400*yd1**3
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=-120.0D0+ 240*xd1
CASE(55)
 S=840.0D0- 1680*xd1 - 1680*zd1 + 3360*xd1*zd1
CASE(76)
 S=-3360.0D0+ 6720*xd1 + 15120*zd1 - 30240*xd1*zd1 - 15120*zd1**2 + 30240*xd1*zd1**2
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=120.0D0- 240*xd1 - 240*yd1 + 480*xd1*yd1
CASE(75)
 S=-840.0D0+ 1680*xd1 + 1680*yd1 - 3360*xd1*yd1 + 1680*zd1 - 3360*xd1*zd1 - 3360*yd1*zd1 + 6720*xd1*yd1*zd1
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=-120.0D0+ 240*xd1 + 720*yd1 - 1440*xd1*yd1 - 720*yd1**2 + 1440*xd1*yd1**2
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=120.0D0- 720*xd1 + 720*xd1**2
CASE(70)
 S=-840.0D0+ 5040*xd1 - 5040*xd1**2 + 1680*zd1 - 10080*xd1*zd1 + 10080*xd1**2*zd1
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=-120.0D0+ 720*xd1 - 720*xd1**2 + 240*yd1 - 1440*xd1*yd1 + 1440*xd1**2*yd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=-120.0D0+ 1440*xd1 - 3600*xd1**2 + 2400*xd1**3
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLZ3=S    
	END FUNCTION DLZ3
 
 
 
REAL FUNCTION DLX4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=1680
CASE(37)
 S=-1680.0D0+ 3360*zd1
CASE(61)
 S=1680.0D0- 10080*zd1 + 10080*zd1**2
CASE(36)
 S=-1680.0D0+ 3360*yd1
CASE(60)
 S=1680.0D0- 3360*yd1 - 3360*zd1 + 6720*yd1*zd1
CASE(59)
 S=1680.0D0- 10080*yd1 + 10080*yd1**2
CASE(35)
 S=-15120.0D0+ 30240*xd1
CASE(58)
 S=15120.0D0- 30240*xd1 - 30240*zd1 + 60480*xd1*zd1
CASE(57)
 S=15120.0D0- 30240*xd1 - 30240*yd1 + 60480*xd1*yd1
CASE(56)
 S=75600.0D0- 332640*xd1 + 332640*xd1**2
  END SELECT
 DLX4=S    
	END FUNCTION DLX4
 
 
REAL FUNCTION DLX3Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=240
CASE(40)
 S=-240.0D0+ 480*zd1
CASE(64)
 S=240.0D0- 1440*zd1 + 1440*zd1**2
CASE(38)
 S=-720.0D0+ 1440*yd1
CASE(63)
 S=720.0D0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
CASE(62)
 S=1440.0D0- 7200*yd1 + 7200*yd1**2
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=-1680.0D0+ 3360*xd1
CASE(60)
 S=1680.0D0- 3360*xd1 - 3360*zd1 + 6720*xd1*zd1
CASE(59)
 S=5040.0D0- 10080*xd1 - 10080*yd1 + 20160*xd1*yd1
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=6720.0D0- 30240*xd1 + 30240*xd1**2
CASE(56)
 S=0.0D0
  END SELECT
 DLX3Y=S    
	END FUNCTION DLX3Y
 
 
 
 
REAL FUNCTION DLX2Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=144
CASE(42)
 S=-144 + 288.d0*zd1
CASE(68)
 S=144 - 864.d0*zd1 + 864.d0*zd1**2
CASE(41)
 S=-720.0D0+ 1440*yd1
CASE(67)
 S=720.0D0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
CASE(66)
 S=2160.0D0- 10080*yd1 + 10080*yd1**2
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=-720.0D0+ 1440*xd1
CASE(63)
 S=720.0D0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
CASE(62)
 S=3600.0D0- 7200*xd1 - 7200*yd1 + 14400*xd1*yd1
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=2160.0D0- 10080*xd1 + 10080*xd1**2
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Y2=S    
	END FUNCTION DLX2Y2
 
 
 
REAL FUNCTION DLXY3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=240
CASE(46)
 S=-240.0D0+ 480*zd1
CASE(73)
 S=240.0D0- 1440*zd1 + 1440*zd1**2
CASE(45)
 S=-1680.0D0+ 3360*yd1
CASE(72)
 S=1680.0D0- 3360*yd1 - 3360*zd1 + 6720*yd1*zd1
CASE(71)
 S=6720.0D0- 30240*yd1 + 30240*yd1**2
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=-720.0D0+ 1440*xd1
CASE(67)
 S=720.0D0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
CASE(66)
 S=5040.0D0- 10080*xd1 - 10080*yd1 + 20160*xd1*yd1
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=1440.0D0- 7200*xd1 + 7200*xd1**2
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY3=S    
	END FUNCTION DLXY3
 
 
 
REAL FUNCTION DLY4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=1680
CASE(50)
 S=-1680.0D0+ 3360*zd1
CASE(79)
 S=1680.0D0- 10080*zd1 + 10080*zd1**2
CASE(49)
 S=-15120.0D0+ 30240*yd1
CASE(78)
 S=15120.0D0- 30240*yd1 - 30240*zd1 + 60480*yd1*zd1
CASE(77)
 S=75600.0D0- 332640*yd1 + 332640*yd1**2
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=-1680.0D0+ 3360*xd1
CASE(72)
 S=1680.0D0- 3360*xd1 - 3360*zd1 + 6720*xd1*zd1
CASE(71)
 S=15120.0D0- 30240*xd1 - 30240*yd1 + 60480*xd1*yd1
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=1680.0D0- 10080*xd1 + 10080*xd1**2
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY4=S    
	END FUNCTION DLY4
 
 
REAL FUNCTION DLX3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=240
CASE(39)
 S=-720.0D0+ 1440*zd1
CASE(65)
 S=1440.0D0- 7200*zd1 + 7200*zd1**2
CASE(21)
 S=0.0D0
CASE(40)
 S=-240.0D0+ 480*yd1
CASE(64)
 S=720.0D0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
CASE(38)
 S=0.0D0
CASE(63)
 S=240.0D0- 1440*yd1 + 1440*yd1**2
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=-1680.0D0+ 3360*xd1
CASE(61)
 S=5040.0D0- 10080*xd1 - 10080*zd1 + 20160*xd1*zd1
CASE(36)
 S=0.0D0
CASE(60)
 S=1680.0D0- 3360*xd1 - 3360*yd1 + 6720*xd1*yd1
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=6720.0D0- 30240*xd1 + 30240*xd1**2
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX3Z=S    
	END FUNCTION DLX3Z
 
 
REAL FUNCTION DLX2YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=48
CASE(43)
 S=-144 + 288.d0*zd1
CASE(69)
 S=288 - 1440*zd1 + 1440*zd1**2
CASE(23)
 S=0.0D0
CASE(42)
 S=-144 + 288.d0*yd1
CASE(68)
 S=432 - 864.d0*yd1 - 864.d0*zd1 + 1728.d0*yd1*zd1
CASE(41)
 S=0.0D0
CASE(67)
 S=288 - 1440*yd1 + 1440*yd1**2
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=-240.0D0+ 480*xd1
CASE(64)
 S=720.0D0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
CASE(38)
 S=0.0D0
CASE(63)
 S=720.0D0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=720.0D0- 3360*xd1 + 3360*xd1**2
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2YZ=S    
	END FUNCTION DLX2YZ
 
 
REAL FUNCTION DLXY2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=48
CASE(48)
 S=-144 + 288.d0*zd1
CASE(74)
 S=288 - 1440*zd1 + 1440*zd1**2
CASE(26)
 S=0.0D0
CASE(46)
 S=-240.0D0+ 480*yd1
CASE(73)
 S=720.0D0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
CASE(45)
 S=0.0D0
CASE(72)
 S=720.0D0- 3360*yd1 + 3360*yd1**2
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=-144 + 288.d0*xd1
CASE(68)
 S=432 - 864.d0*xd1 - 864.d0*zd1 + 1728.d0*xd1*zd1
CASE(41)
 S=0.0D0
CASE(67)
 S=720.0D0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=288 - 1440*xd1 + 1440*xd1**2
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY2Z=S    
	END FUNCTION DLXY2Z
 
 
 
 
REAL FUNCTION DLY3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=240
CASE(51)
 S=-720.0D0+ 1440*zd1
CASE(80)
 S=1440.0D0- 7200*zd1 + 7200*zd1**2
CASE(34)
 S=0.0D0
CASE(50)
 S=-1680.0D0+ 3360*yd1
CASE(79)
 S=5040.0D0- 10080*yd1 - 10080*zd1 + 20160*yd1*zd1
CASE(49)
 S=0.0D0
CASE(78)
 S=6720.0D0- 30240*yd1 + 30240*yd1**2
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=-240.0D0+ 480*xd1
CASE(73)
 S=720.0D0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
CASE(45)
 S=0.0D0
CASE(72)
 S=1680.0D0- 3360*xd1 - 3360*yd1 + 6720*xd1*yd1
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=240.0D0- 1440*xd1 + 1440*xd1**2
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY3Z=S    
	END FUNCTION DLY3Z
 
 
 
 
 
 
REAL FUNCTION DLX2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=144
CASE(44)
 S=-720.0D0+ 1440*zd1
CASE(70)
 S=2160.0D0- 10080*zd1 + 10080*zd1**2
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=-144 + 288.d0*yd1
CASE(69)
 S=720.0D0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=144 - 864.d0*yd1 + 864.d0*yd1**2
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=-720.0D0+ 1440*xd1
CASE(65)
 S=3600.0D0- 7200*xd1 - 7200*zd1 + 14400*xd1*zd1
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=720.0D0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=2160.0D0- 10080*xd1 + 10080*xd1**2
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Z2=S    
	END FUNCTION DLX2Z2
 
 
 
REAL FUNCTION DLXYZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=48
CASE(47)
 S=-240.0D0+ 480*zd1
CASE(75)
 S=720.0D0- 3360*zd1 + 3360*zd1**2
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=-144 + 288.d0*yd1
CASE(74)
 S=720.0D0- 1440*yd1 - 1440*zd1 + 2880*yd1*zd1
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=288 - 1440*yd1 + 1440*yd1**2
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=-144 + 288.d0*xd1
CASE(69)
 S=720.0D0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=432 - 864.d0*xd1 - 864.d0*yd1 + 1728.d0*xd1*yd1
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=288 - 1440*xd1 + 1440*xd1**2
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXYZ2=S    
	END FUNCTION DLXYZ2
 
 
REAL FUNCTION DLY2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=144
CASE(52)
 S=-720.0D0+ 1440*zd1
CASE(81)
 S=2160.0D0- 10080*zd1 + 10080*zd1**2
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=-720.0D0+ 1440*yd1
CASE(80)
 S=3600.0D0- 7200*yd1 - 7200*zd1 + 14400*yd1*zd1
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=2160.0D0- 10080*yd1 + 10080*yd1**2
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=-144 + 288.d0*xd1
CASE(74)
 S=720.0D0- 1440*xd1 - 1440*zd1 + 2880*xd1*zd1
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=720.0D0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=144 - 864.d0*xd1 + 864.d0*xd1**2
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY2Z2=S    
	END FUNCTION DLY2Z2
 
 
 
REAL FUNCTION DLXZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=240
CASE(55)
 S=-1680.0D0+ 3360*zd1
CASE(76)
 S=6720.0D0- 30240*zd1 + 30240*zd1**2
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=-240.0D0+ 480*yd1
CASE(75)
 S=1680.0D0- 3360*yd1 - 3360*zd1 + 6720*yd1*zd1
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=240.0D0- 1440*yd1 + 1440*yd1**2
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=-720.0D0+ 1440*xd1
CASE(70)
 S=5040.0D0- 10080*xd1 - 10080*zd1 + 20160*xd1*zd1
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=720.0D0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=1440.0D0- 7200*xd1 + 7200*xd1**2
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXZ3=S    
	END FUNCTION DLXZ3
 
 
 
REAL FUNCTION DLYZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=240
CASE(53)
 S=-1680.0D0+ 3360*zd1
CASE(82)
 S=6720.0D0- 30240*zd1 + 30240*zd1**2
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=-720.0D0+ 1440*yd1
CASE(81)
 S=5040.0D0- 10080*yd1 - 10080*zd1 + 20160*yd1*zd1
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=1440.0D0- 7200*yd1 + 7200*yd1**2
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=-240.0D0+ 480*xd1
CASE(75)
 S=1680.0D0- 3360*xd1 - 3360*zd1 + 6720*xd1*zd1
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=720.0D0- 1440*xd1 - 1440*yd1 + 2880*xd1*yd1
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=240.0D0- 1440*xd1 + 1440*xd1**2
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLYZ3=S    
	END FUNCTION DLYZ3
 
 
 
 
REAL FUNCTION DLZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=1680
CASE(54)
 S=-15120.0D0+ 30240*zd1
CASE(83)
 S=75600.0D0- 332640*zd1 + 332640*zd1**2
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=-1680.0D0+ 3360*yd1
CASE(82)
 S=15120.0D0- 30240*yd1 - 30240*zd1 + 60480*yd1*zd1
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=1680.0D0- 10080*yd1 + 10080*yd1**2
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=-1680.0D0+ 3360*xd1
CASE(76)
 S=15120.0D0- 30240*xd1 - 30240*zd1 + 60480*xd1*zd1
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=1680.0D0- 3360*xd1 - 3360*yd1 + 6720*xd1*yd1
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=1680.0D0- 10080*xd1 + 10080*xd1**2
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLZ4=S    
	END FUNCTION DLZ4
 
 
 
REAL FUNCTION DLX5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=30240
CASE(58)
 S=-30240.0D0+ 60480*zd1
CASE(57)
 S=-30240.0D0+ 60480*yd1
CASE(56)
 S=-332640.0D0+ 665280*xd1
  END SELECT
 DLX5=S    
	END FUNCTION DLX5
 
REAL FUNCTION DLX4Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=3360
CASE(60)
 S=-3360.0D0+ 6720*zd1
CASE(59)
 S=-10080.0D0+ 20160*yd1
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=-30240.0D0+ 60480*xd1
CASE(56)
 S=0.0D0
  END SELECT
 DLX4Y=S    
	END FUNCTION DLX4Y
 
 
 
 
REAL FUNCTION DLX3Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=1440
CASE(63)
 S=-1440.0D0+ 2880*zd1
CASE(62)
 S=-7200.0D0+ 14400*yd1
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=-10080.0D0+ 20160*xd1
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX3Y2=S    
	END FUNCTION DLX3Y2
 
 
 
REAL FUNCTION DLX2Y3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=1440
CASE(67)
 S=-1440.0D0+ 2880*zd1
CASE(66)
 S=-10080.0D0+ 20160*yd1
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=-7200.0D0+ 14400*xd1
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Y3=S    
	END FUNCTION DLX2Y3
 
 
 
 
REAL FUNCTION DLXY4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=3360
CASE(72)
 S=-3360.0D0+ 6720*zd1
CASE(71)
 S=-30240.0D0+ 60480*yd1
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=-10080.0D0+ 20160*xd1
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY4=S    
	END FUNCTION DLXY4
 
 
REAL FUNCTION DLY5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=30240
CASE(78)
 S=-30240.0D0+ 60480*zd1
CASE(77)
 S=-332640.0D0+ 665280*yd1
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=-30240.0D0+ 60480*xd1
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY5=S    
	END FUNCTION DLY5
 
REAL FUNCTION DLX4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=3360
CASE(61)
 S=-10080.0D0+ 20160*zd1
CASE(36)
 S=0.0D0
CASE(60)
 S=-3360.0D0+ 6720*yd1
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=-30240.0D0+ 60480*xd1
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX4Z=S    
	END FUNCTION DLX4Z
 
 
 
REAL FUNCTION DLX3YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=480
CASE(64)
 S=-1440.0D0+ 2880*zd1
CASE(38)
 S=0.0D0
CASE(63)
 S=-1440.0D0+ 2880*yd1
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=-3360.0D0+ 6720*xd1
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX3YZ=S    
	END FUNCTION DLX3YZ
 
 
 
REAL FUNCTION DLX2Y2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=288
CASE(68)
 S=-864 + 1728.d0*zd1
CASE(41)
 S=0.0D0
CASE(67)
 S=-1440.0D0+ 2880*yd1
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=-1440.0D0+ 2880*xd1
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Y2Z=S    
	END FUNCTION DLX2Y2Z
 
REAL FUNCTION DLXY3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=480
CASE(73)
 S=-1440.0D0+ 2880*zd1
CASE(45)
 S=0.0D0
CASE(72)
 S=-3360.0D0+ 6720*yd1
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=-1440.0D0+ 2880*xd1
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY3Z=S    
	END FUNCTION DLXY3Z
 
REAL FUNCTION DLY4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=3360
CASE(79)
 S=-10080.0D0+ 20160*zd1
CASE(49)
 S=0.0D0
CASE(78)
 S=-30240.0D0+ 60480*yd1
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=-3360.0D0+ 6720*xd1
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY4Z=S    
	END FUNCTION DLY4Z
 
 
REAL FUNCTION DLX3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=1440
CASE(65)
 S=-7200.0D0+ 14400*zd1
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=-1440.0D0+ 2880*yd1
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=-10080.0D0+ 20160*xd1
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX3Z2=S    
	END FUNCTION DLX3Z2
 
 
REAL FUNCTION DLX2YZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=288
CASE(69)
 S=-1440.0D0+ 2880*zd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=-864 + 1728.d0*yd1
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=-1440.0D0+ 2880*xd1
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2YZ2=S    
	END FUNCTION DLX2YZ2
	
 REAL FUNCTION DLXY2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
 
 
 CASE(3)
 S=0.0D0

CASE(9)
 S=0.0D0

CASE(16)
 S=0.0D0

CASE(33)
 S=0.0D0

CASE(54)
 S=0.0D0

CASE(83)
 S=0.0D0

CASE(2)
 S=0.0D0

CASE(8)
 S=0.0D0

CASE(15)
 S=0.0D0

CASE(32)
 S=0.0D0

CASE(53)
 S=0.0D0

CASE(82)
 S=0.0D0

CASE(7)
 S=0.0D0

CASE(14)
 S=0.0D0

CASE(31)
 S=0.0D0

CASE(52)
 S=0.0D0

CASE(81)
 S=0.0D0

CASE(13)
 S=0.0D0

CASE(30)
 S=0.0D0

CASE(51)
 S=0.0D0

CASE(80)
 S=0.0D0

CASE(34)
 S=0.0D0

CASE(50)
 S=0.0D0

CASE(79)
 S=0.0D0

CASE(49)
 S=0.0D0

CASE(78)
 S=0.0D0

CASE(77)
 S=0.0D0

CASE(1)
 S=0.0D0

CASE(6)
 S=0.0D0

CASE(19)
 S=0.0D0

CASE(29)
 S=0.0D0

CASE(55)
 S=0.0D0

CASE(76)
 S=0.0D0

CASE(5)
 S=0.0D0

CASE(17)
 S=0.0D0

CASE(28)
 S=0.0D0

CASE(47)
 S=0.0D0

CASE(75)
 S=0.0D0

CASE(12)
 S=0.0D0

CASE(27)
 S=0.0D0

CASE(48)
 S=288

CASE(74)
 S=-1440.0D0+ 2880*zd1

CASE(26)
 S=0.0D0

CASE(46)
 S=0.0D0

CASE(73)
 S=-1440.0D0+ 2880*yd1

CASE(45)
 S=0.0D0

CASE(72)
 S=0.0D0

CASE(71)
 S=0.0D0

CASE(4)
 S=0.0D0

CASE(18)
 S=0.0D0

CASE(25)
 S=0.0D0

CASE(44)
 S=0.0D0

CASE(70)
 S=0.0D0

CASE(11)
 S=0.0D0

CASE(24)
 S=0.0D0

CASE(43)
 S=0.0D0

CASE(69)
 S=0.0D0

CASE(23)
 S=0.0D0

CASE(42)
 S=0.0D0

CASE(68)
 S=-864 + 1728.d0*xd1

CASE(41)
 S=0.0D0

CASE(67)
 S=0.0D0

CASE(66)
 S=0.0D0

CASE(10)
 S=0.0D0

CASE(22)
 S=0.0D0

CASE(39)
 S=0.0D0

CASE(65)
 S=0.0D0

CASE(21)
 S=0.0D0

CASE(40)
 S=0.0D0

CASE(64)
 S=0.0D0

CASE(38)
 S=0.0D0

CASE(63)
 S=0.0D0

CASE(62)
 S=0.0D0

CASE(20)
 S=0.0D0

CASE(37)
 S=0.0D0

CASE(61)
 S=0.0D0

CASE(36)
 S=0.0D0

CASE(60)
 S=0.0D0

CASE(59)
 S=0.0D0

CASE(35)
 S=0.0D0

CASE(58)
 S=0.0D0

CASE(57)
 S=0.0D0

CASE(56)
 S=0.0D0
 
   END SELECT
 DLXY2Z2=S    
	END FUNCTION DLXY2Z2
 
 
 
 
 
 
 
 
 
 
 
 
 
REAL FUNCTION DLY3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=1440
CASE(80)
 S=-7200.0D0+ 14400*zd1
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=-10080.0D0+ 20160*yd1
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=-1440.0D0+ 2880*xd1
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY3Z2=S    
	END FUNCTION DLY3Z2
 
 
 
REAL FUNCTION DLX2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=1440
CASE(70)
 S=-10080.0D0+ 20160*zd1
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=-1440.0D0+ 2880*yd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=-7200.0D0+ 14400*xd1
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Z3=S    
	END FUNCTION DLX2Z3
 
 
REAL FUNCTION DLXYZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=480
CASE(75)
 S=-3360.0D0+ 6720*zd1
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=-1440.0D0+ 2880*yd1
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=-1440.0D0+ 2880*xd1
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXYZ3=S    
	END FUNCTION DLXYZ3
 
 
 
REAL FUNCTION DLY2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=1440
CASE(81)
 S=-10080.0D0+ 20160*zd1
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=-7200.0D0+ 14400*yd1
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=-1440.0D0+ 2880*xd1
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY2Z3=S    
	END FUNCTION DLY2Z3
 
 
REAL FUNCTION DLXZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=3360
CASE(76)
 S=-30240.0D0+ 60480*zd1
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=-3360.0D0+ 6720*yd1
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=-10080.0D0+ 20160*xd1
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXZ4=S    
	END FUNCTION DLXZ4
 
 
REAL FUNCTION DLYZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=3360
CASE(82)
 S=-30240.0D0+ 60480*zd1
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=-10080.0D0+ 20160*yd1
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=-3360.0D0+ 6720*xd1
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLYZ4=S    
	END FUNCTION DLYZ4
 
 
REAL FUNCTION DLZ5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=30240
CASE(83)
 S=-332640.0D0+ 665280*zd1
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=-30240.0D0+ 60480*yd1
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=-30240.0D0+ 60480*xd1
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLZ5=S    
	END FUNCTION DLZ5
 
 
REAL FUNCTION DLX6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=665280
  END SELECT
 DLX6=S    
	END FUNCTION DLX6
 
 
 
REAL FUNCTION DLX5Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=60480
CASE(56)
 S=0.0D0
  END SELECT
 DLX5Y=S    
	END FUNCTION DLX5Y
 
 
REAL FUNCTION DLX4Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=20160
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX4Y2=S    
	END FUNCTION DLX4Y2
 
 
REAL FUNCTION DLX3Y3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=14400
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX3Y3=S    
	END FUNCTION DLX3Y3
 
 
REAL FUNCTION DLX2Y4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=20160
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Y4=S    
	END FUNCTION DLX2Y4
 
 
 
REAL FUNCTION DLXY5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=60480
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY5=S    
	END FUNCTION DLXY5
 
 
REAL FUNCTION DLY6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=665280
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY6=S    
	END FUNCTION DLY6
 
 
REAL FUNCTION DLX5Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=60480
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX5Z=S    
	END FUNCTION DLX5Z
 
 
 
 
REAL FUNCTION DLX4YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=6720
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX4YZ=S    
	END FUNCTION DLX4YZ
 
 
 
REAL FUNCTION DLX3Y2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=2880
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX3Y2Z=S    
	END FUNCTION DLX3Y2Z
 
 
 
REAL FUNCTION DLX2Y3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=2880
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Y3Z=S    
	END FUNCTION DLX2Y3Z
 
 
 
 
 
 
REAL FUNCTION DLXY4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=6720
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY4Z=S    
	END FUNCTION DLXY4Z
 
 
 
REAL FUNCTION DLY5Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=60480
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY5Z=S    
	END FUNCTION DLY5Z
 
 
 
REAL FUNCTION DLX4Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=20160
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX4Z2=S    
	END FUNCTION DLX4Z2
 
 
REAL FUNCTION DLX3YZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=2880
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX3YZ2=S    
	END FUNCTION DLX3YZ2
 
 
 
REAL FUNCTION DLX2Y2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=1728
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Y2Z2=S    
	END FUNCTION DLX2Y2Z2
 
 
REAL FUNCTION DLXY3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=2880
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY3Z2=S    
	END FUNCTION DLXY3Z2
 
 
 
 
REAL FUNCTION DLY4Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=20160
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY4Z2=S    
	END FUNCTION DLY4Z2
 
 
REAL FUNCTION DLX3Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=14400
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX3Z3=S    
	END FUNCTION DLX3Z3
 
 
REAL FUNCTION DLX2YZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=2880
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2YZ3=S    
	END FUNCTION DLX2YZ3
 
 
 
 
 
REAL FUNCTION DLXY2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=2880
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLXY2Z3=S    
	END FUNCTION DLXY2Z3
 
 
 
 
REAL FUNCTION DLY3Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=14400
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLY3Z3=S    
	END FUNCTION DLY3Z3
 
 
 
REAL FUNCTION DLX2Z4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=20160
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
  END SELECT
 DLX2Z4=S    
	END FUNCTION DLX2Z4
 
 
REAL FUNCTION DLXYZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=6720
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
 END SELECT
 DLXYZ4=S    
	END FUNCTION DLXYZ4
 
 
REAL FUNCTION DLY2Z4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=20160
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
 END SELECT  
	DLY2Z4=S    
	END FUNCTION DLY2Z4 
 
 
 
REAL FUNCTION DLXZ5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=60480
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
 END SELECT  
	DLXZ5=S    
	END FUNCTION DLXZ5 
 
 
 
REAL FUNCTION DLYZ5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=0.0D0
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=60480
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
 END SELECT  
	DLYZ5=S    
	END FUNCTION DLYZ5
 
 
REAL FUNCTION DLZ6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S
                    
S=0.0D0
SELECTCASE (NDERIVATIVE)
CASE(3)
 S=0.0D0
CASE(9)
 S=0.0D0
CASE(16)
 S=0.0D0
CASE(33)
 S=0.0D0
CASE(54)
 S=0.0D0
CASE(83)
 S=665280
CASE(2)
 S=0.0D0
CASE(8)
 S=0.0D0
CASE(15)
 S=0.0D0
CASE(32)
 S=0.0D0
CASE(53)
 S=0.0D0
CASE(82)
 S=0.0D0
CASE(7)
 S=0.0D0
CASE(14)
 S=0.0D0
CASE(31)
 S=0.0D0
CASE(52)
 S=0.0D0
CASE(81)
 S=0.0D0
CASE(13)
 S=0.0D0
CASE(30)
 S=0.0D0
CASE(51)
 S=0.0D0
CASE(80)
 S=0.0D0
CASE(34)
 S=0.0D0
CASE(50)
 S=0.0D0
CASE(79)
 S=0.0D0
CASE(49)
 S=0.0D0
CASE(78)
 S=0.0D0
CASE(77)
 S=0.0D0
CASE(1)
 S=0.0D0
CASE(6)
 S=0.0D0
CASE(19)
 S=0.0D0
CASE(29)
 S=0.0D0
CASE(55)
 S=0.0D0
CASE(76)
 S=0.0D0
CASE(5)
 S=0.0D0
CASE(17)
 S=0.0D0
CASE(28)
 S=0.0D0
CASE(47)
 S=0.0D0
CASE(75)
 S=0.0D0
CASE(12)
 S=0.0D0
CASE(27)
 S=0.0D0
CASE(48)
 S=0.0D0
CASE(74)
 S=0.0D0
CASE(26)
 S=0.0D0
CASE(46)
 S=0.0D0
CASE(73)
 S=0.0D0
CASE(45)
 S=0.0D0
CASE(72)
 S=0.0D0
CASE(71)
 S=0.0D0
CASE(4)
 S=0.0D0
CASE(18)
 S=0.0D0
CASE(25)
 S=0.0D0
CASE(44)
 S=0.0D0
CASE(70)
 S=0.0D0
CASE(11)
 S=0.0D0
CASE(24)
 S=0.0D0
CASE(43)
 S=0.0D0
CASE(69)
 S=0.0D0
CASE(23)
 S=0.0D0
CASE(42)
 S=0.0D0
CASE(68)
 S=0.0D0
CASE(41)
 S=0.0D0
CASE(67)
 S=0.0D0
CASE(66)
 S=0.0D0
CASE(10)
 S=0.0D0
CASE(22)
 S=0.0D0
CASE(39)
 S=0.0D0
CASE(65)
 S=0.0D0
CASE(21)
 S=0.0D0
CASE(40)
 S=0.0D0
CASE(64)
 S=0.0D0
CASE(38)
 S=0.0D0
CASE(63)
 S=0.0D0
CASE(62)
 S=0.0D0
CASE(20)
 S=0.0D0
CASE(37)
 S=0.0D0
CASE(61)
 S=0.0D0
CASE(36)
 S=0.0D0
CASE(60)
 S=0.0D0
CASE(59)
 S=0.0D0
CASE(35)
 S=0.0D0
CASE(58)
 S=0.0D0
CASE(57)
 S=0.0D0
CASE(56)
 S=0.0D0
 END SELECT  
	DLZ6=S    
	END FUNCTION DLZ6
	
	
REAL FUNCTION TL3DZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
  CASE(3) ; S= 1/hxc
  CASE(9) ; S= zd1/hxc**2
  CASE(19) ; S= zd1**2/(2.0D0*hxc**3)
  CASE(34) ; S= zd1**3/(6.0D0*hxc**4)
  CASE(55) ; S= zd1**4/(24.0D0*hxc**5)
  CASE(83) ; S= zd1**5/(120.0D0*hxc**6)
  CASE(2) ; S= 0.0D0
  CASE(8) ; S= yd1/hxc**2
  CASE(18) ; S= (yd1*zd1)/hxc**3
  CASE(33) ; S= (yd1*zd1**2)/(2.0D0*hxc**4)
  CASE(54) ; S= (yd1*zd1**3)/(6.0D0*hxc**5)
  CASE(82) ; S= (yd1*zd1**4)/(24.0D0*hxc**6)
  CASE(7) ; S= 0.0D0
  CASE(17) ; S= yd1**2/(2.0D0*hxc**3)
  CASE(32) ; S= (yd1**2*zd1)/(2.0D0*hxc**4)
  CASE(53) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**5)
  CASE(81) ; S= (yd1**2*zd1**3)/(12.0D0*hxc**6)
  CASE(16) ; S= 0.0D0
  CASE(31) ; S= yd1**3/(6.0D0*hxc**4)
  CASE(52) ; S= (yd1**3*zd1)/(6.0D0*hxc**5)
  CASE(80) ; S= (yd1**3*zd1**2)/(12.0D0*hxc**6)
  CASE(30) ; S= 0.0D0
  CASE(51) ; S= yd1**4/(24.0D0*hxc**5)
  CASE(79) ; S= (yd1**4*zd1)/(24.0D0*hxc**6)
  CASE(50) ; S= 0.0D0
  CASE(78) ; S= yd1**5/(120.0D0*hxc**6)
  CASE(77) ; S= 0.0D0
  CASE(1) ; S= 0.0D0
  CASE(6) ; S= xd1/hxc**2
  CASE(14) ; S= (xd1*zd1)/hxc**3
  CASE(27) ; S= (xd1*zd1**2)/(2.0D0*hxc**4)
  CASE(49) ; S= (xd1*zd1**3)/(6.0D0*hxc**5)
  CASE(76) ; S= (xd1*zd1**4)/(24.0D0*hxc**6)
  CASE(5) ; S= 0.0D0
  CASE(15) ; S= (xd1*yd1)/hxc**3
  CASE(29) ; S= (xd1*yd1*zd1)/hxc**4
  CASE(48) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**5)
  CASE(75) ; S= (xd1*yd1*zd1**3)/(6.0D0*hxc**6)
  CASE(13) ; S= 0.0D0
  CASE(28) ; S= (xd1*yd1**2)/(2.0D0*hxc**4)
  CASE(47) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**5)
  CASE(74) ; S= (xd1*yd1**2*zd1**2)/(4.0D0*hxc**6)
  CASE(26) ; S= 0.0D0
  CASE(46) ; S= (xd1*yd1**3)/(6.0D0*hxc**5)
  CASE(73) ; S= (xd1*yd1**3*zd1)/(6.0D0*hxc**6)
  CASE(45) ; S= 0.0D0
  CASE(72) ; S= (xd1*yd1**4)/(24.0D0*hxc**6)
  CASE(71) ; S= 0.0D0
  CASE(4) ; S= 0.0D0
  CASE(12) ; S= xd1**2/(2.0D0*hxc**3)
  CASE(24) ; S= (xd1**2*zd1)/(2.0D0*hxc**4)
  CASE(44) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**5)
  CASE(70) ; S= (xd1**2*zd1**3)/(12.0D0*hxc**6)
  CASE(11) ; S= 0.0D0
  CASE(25) ; S= (xd1**2*yd1)/(2.0D0*hxc**4)
  CASE(43) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**5)
  CASE(69) ; S= (xd1**2*yd1*zd1**2)/(4.0D0*hxc**6)
  CASE(23) ; S= 0.0D0
  CASE(42) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**5)
  CASE(68) ; S= (xd1**2*yd1**2*zd1)/(4.0D0*hxc**6)
  CASE(41) ; S= 0.0D0
  CASE(67) ; S= (xd1**2*yd1**3)/(12.0D0*hxc**6)
  CASE(66) ; S= 0.0D0
  CASE(10) ; S= 0.0D0
  CASE(22) ; S= xd1**3/(6.0D0*hxc**4)
  CASE(39) ; S= (xd1**3*zd1)/(6.0D0*hxc**5)
  CASE(65) ; S= (xd1**3*zd1**2)/(12.0D0*hxc**6)
  CASE(21) ; S= 0.0D0
  CASE(40) ; S= (xd1**3*yd1)/(6.0D0*hxc**5)
  CASE(64) ; S= (xd1**3*yd1*zd1)/(6.0D0*hxc**6)
  CASE(38) ; S= 0.0D0
  CASE(63) ; S= (xd1**3*yd1**2)/(12.0D0*hxc**6)
  CASE(62) ; S= 0.0D0
  CASE(20) ; S= 0.0D0
  CASE(37) ; S= xd1**4/(24.0D0*hxc**5)
  CASE(61) ; S= (xd1**4*zd1)/(24.0D0*hxc**6)
  CASE(36) ; S= 0.0D0
  CASE(60) ; S= (xd1**4*yd1)/(24.0D0*hxc**6)
  CASE(59) ; S= 0.0D0
  CASE(35) ; S= 0.0D0
  CASE(58) ; S= xd1**5/(120.0D0*hxc**6)
  CASE(57) ; S= 0.0D0
  CASE(56) ; S= 0.0D0
 
 
 End select

    TL3DZ = s

  end function
 
  
 
 
  REAL FUNCTION TL3DZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
  CASE(3) ; S= 0.0D0
  CASE(9) ; S= hxc**(-2)
  CASE(19) ; S= zd1/hxc**3
  CASE(34) ; S= zd1**2/(2.0D0*hxc**4)
  CASE(55) ; S= zd1**3/(6.0D0*hxc**5)
  CASE(83) ; S= zd1**4/(24.0D0*hxc**6)
  CASE(2) ; S= 0.0D0
  CASE(8) ; S= 0.0D0
  CASE(18) ; S= yd1/hxc**3
  CASE(33) ; S= (yd1*zd1)/hxc**4
  CASE(54) ; S= (yd1*zd1**2)/(2.0D0*hxc**5)
  CASE(82) ; S= (yd1*zd1**3)/(6.0D0*hxc**6)
  CASE(7) ; S= 0.0D0
  CASE(17) ; S= 0.0D0
  CASE(32) ; S= yd1**2/(2.0D0*hxc**4)
  CASE(53) ; S= (yd1**2*zd1)/(2.0D0*hxc**5)
  CASE(81) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**6)
  CASE(16) ; S= 0.0D0
  CASE(31) ; S= 0.0D0
  CASE(52) ; S= yd1**3/(6.0D0*hxc**5)
  CASE(80) ; S= (yd1**3*zd1)/(6.0D0*hxc**6)
  CASE(30) ; S= 0.0D0
  CASE(51) ; S= 0.0D0
  CASE(79) ; S= yd1**4/(24.0D0*hxc**6)
  CASE(50) ; S= 0.0D0
  CASE(78) ; S= 0.0D0
  CASE(77) ; S= 0.0D0
  CASE(1) ; S= 0.0D0
  CASE(6) ; S= 0.0D0
  CASE(14) ; S= xd1/hxc**3
  CASE(27) ; S= (xd1*zd1)/hxc**4
  CASE(49) ; S= (xd1*zd1**2)/(2.0D0*hxc**5)
  CASE(76) ; S= (xd1*zd1**3)/(6.0D0*hxc**6)
  CASE(5) ; S= 0.0D0
  CASE(15) ; S= 0.0D0
  CASE(29) ; S= (xd1*yd1)/hxc**4
  CASE(48) ; S= (xd1*yd1*zd1)/hxc**5
  CASE(75) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**6)
  CASE(13) ; S= 0.0D0
  CASE(28) ; S= 0.0D0
  CASE(47) ; S= (xd1*yd1**2)/(2.0D0*hxc**5)
  CASE(74) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**6)
  CASE(26) ; S= 0.0D0
  CASE(46) ; S= 0.0D0
  CASE(73) ; S= (xd1*yd1**3)/(6.0D0*hxc**6)
  CASE(45) ; S= 0.0D0
  CASE(72) ; S= 0.0D0
  CASE(71) ; S= 0.0D0
  CASE(4) ; S= 0.0D0
  CASE(12) ; S= 0.0D0
  CASE(24) ; S= xd1**2/(2.0D0*hxc**4)
  CASE(44) ; S= (xd1**2*zd1)/(2.0D0*hxc**5)
  CASE(70) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**6)
  CASE(11) ; S= 0.0D0
  CASE(25) ; S= 0.0D0
  CASE(43) ; S= (xd1**2*yd1)/(2.0D0*hxc**5)
  CASE(69) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**6)
  CASE(23) ; S= 0.0D0
  CASE(42) ; S= 0.0D0
  CASE(68) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**6)
  CASE(41) ; S= 0.0D0
  CASE(67) ; S= 0.0D0
  CASE(66) ; S= 0.0D0
  CASE(10) ; S= 0.0D0
  CASE(22) ; S= 0.0D0
  CASE(39) ; S= xd1**3/(6.0D0*hxc**5)
  CASE(65) ; S= (xd1**3*zd1)/(6.0D0*hxc**6)
  CASE(21) ; S= 0.0D0
  CASE(40) ; S= 0.0D0
  CASE(64) ; S= (xd1**3*yd1)/(6.0D0*hxc**6)
  CASE(38) ; S= 0.0D0
  CASE(63) ; S= 0.0D0
  CASE(62) ; S= 0.0D0
  CASE(20) ; S= 0.0D0
  CASE(37) ; S= 0.0D0
  CASE(61) ; S= xd1**4/(24.0D0*hxc**6)
  CASE(36) ; S= 0.0D0
  CASE(60) ; S= 0.0D0
  CASE(59) ; S= 0.0D0
  CASE(35) ; S= 0.0D0
  CASE(58) ; S= 0.0D0
  CASE(57) ; S= 0.0D0
  CASE(56) ; S= 0.0D0
  
  
  
  
   End select

    TL3DZ2 = s

  end function
 
 
 
  REAL FUNCTION TL3DZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
  CASE(3) ; S= 0.0D0
  CASE(9) ; S= 0.0D0
  CASE(19) ; S= hxc**(-3)
  CASE(34) ; S= zd1/hxc**4
  CASE(55) ; S= zd1**2/(2.0D0*hxc**5)
  CASE(83) ; S= zd1**3/(6.0D0*hxc**6)
  CASE(2) ; S= 0.0D0
  CASE(8) ; S= 0.0D0
  CASE(18) ; S= 0.0D0
  CASE(33) ; S= yd1/hxc**4
  CASE(54) ; S= (yd1*zd1)/hxc**5
  CASE(82) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
  CASE(7) ; S= 0.0D0
  CASE(17) ; S= 0.0D0
  CASE(32) ; S= 0.0D0
  CASE(53) ; S= yd1**2/(2.0D0*hxc**5)
  CASE(81) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
  CASE(16) ; S= 0.0D0
  CASE(31) ; S= 0.0D0
  CASE(52) ; S= 0.0D0
  CASE(80) ; S= yd1**3/(6.0D0*hxc**6)
  CASE(30) ; S= 0.0D0
  CASE(51) ; S= 0.0D0
  CASE(79) ; S= 0.0D0
  CASE(50) ; S= 0.0D0
  CASE(78) ; S= 0.0D0
  CASE(77) ; S= 0.0D0
  CASE(1) ; S= 0.0D0
  CASE(6) ; S= 0.0D0
  CASE(14) ; S= 0.0D0
  CASE(27) ; S= xd1/hxc**4
  CASE(49) ; S= (xd1*zd1)/hxc**5
  CASE(76) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
  CASE(5) ; S= 0.0D0
  CASE(15) ; S= 0.0D0
  CASE(29) ; S= 0.0D0
  CASE(48) ; S= (xd1*yd1)/hxc**5
  CASE(75) ; S= (xd1*yd1*zd1)/hxc**6
  CASE(13) ; S= 0.0D0
  CASE(28) ; S= 0.0D0
  CASE(47) ; S= 0.0D0
  CASE(74) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
  CASE(26) ; S= 0.0D0
  CASE(46) ; S= 0.0D0
  CASE(73) ; S= 0.0D0
  CASE(45) ; S= 0.0D0
  CASE(72) ; S= 0.0D0
  CASE(71) ; S= 0.0D0
  CASE(4) ; S= 0.0D0
  CASE(12) ; S= 0.0D0
  CASE(24) ; S= 0.0D0
  CASE(44) ; S= xd1**2/(2.0D0*hxc**5)
  CASE(70) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
  CASE(11) ; S= 0.0D0
  CASE(25) ; S= 0.0D0
  CASE(43) ; S= 0.0D0
  CASE(69) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
  CASE(23) ; S= 0.0D0
  CASE(42) ; S= 0.0D0
  CASE(68) ; S= 0.0D0
  CASE(41) ; S= 0.0D0
  CASE(67) ; S= 0.0D0
  CASE(66) ; S= 0.0D0
  CASE(10) ; S= 0.0D0
  CASE(22) ; S= 0.0D0
  CASE(39) ; S= 0.0D0
  CASE(65) ; S= xd1**3/(6.0D0*hxc**6)
  CASE(21) ; S= 0.0D0
  CASE(40) ; S= 0.0D0
  CASE(64) ; S= 0.0D0
  CASE(38) ; S= 0.0D0
  CASE(63) ; S= 0.0D0
  CASE(62) ; S= 0.0D0
  CASE(20) ; S= 0.0D0
  CASE(37) ; S= 0.0D0
  CASE(61) ; S= 0.0D0
  CASE(36) ; S= 0.0D0
  CASE(60) ; S= 0.0D0
  CASE(59) ; S= 0.0D0
  CASE(35) ; S= 0.0D0
  CASE(58) ; S= 0.0D0
  CASE(57) ; S= 0.0D0
  CASE(56) ; S= 0.0D0
  
  
  
  
   End select

    TL3DZ3 = s

  end function
 
 
 
  REAL FUNCTION TL3DZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
  CASE(3) ; S= 0.0D0
  CASE(9) ; S= 0.0D0
  CASE(19) ; S= 0.0D0
  CASE(34) ; S= hxc**(-4)
  CASE(55) ; S= zd1/hxc**5
  CASE(83) ; S= zd1**2/(2.0D0*hxc**6)
  CASE(2) ; S= 0.0D0
  CASE(8) ; S= 0.0D0
  CASE(18) ; S= 0.0D0
  CASE(33) ; S= 0.0D0
  CASE(54) ; S= yd1/hxc**5
  CASE(82) ; S= (yd1*zd1)/hxc**6
  CASE(7) ; S= 0.0D0
  CASE(17) ; S= 0.0D0
  CASE(32) ; S= 0.0D0
  CASE(53) ; S= 0.0D0
  CASE(81) ; S= yd1**2/(2.0D0*hxc**6)
  CASE(16) ; S= 0.0D0
  CASE(31) ; S= 0.0D0
  CASE(52) ; S= 0.0D0
  CASE(80) ; S= 0.0D0
  CASE(30) ; S= 0.0D0
  CASE(51) ; S= 0.0D0
  CASE(79) ; S= 0.0D0
  CASE(50) ; S= 0.0D0
  CASE(78) ; S= 0.0D0
  CASE(77) ; S= 0.0D0
  CASE(1) ; S= 0.0D0
  CASE(6) ; S= 0.0D0
  CASE(14) ; S= 0.0D0
  CASE(27) ; S= 0.0D0
  CASE(49) ; S= xd1/hxc**5
  CASE(76) ; S= (xd1*zd1)/hxc**6
  CASE(5) ; S= 0.0D0
  CASE(15) ; S= 0.0D0
  CASE(29) ; S= 0.0D0
  CASE(48) ; S= 0.0D0
  CASE(75) ; S= (xd1*yd1)/hxc**6
  CASE(13) ; S= 0.0D0
  CASE(28) ; S= 0.0D0
  CASE(47) ; S= 0.0D0
  CASE(74) ; S= 0.0D0
  CASE(26) ; S= 0.0D0
  CASE(46) ; S= 0.0D0
  CASE(73) ; S= 0.0D0
  CASE(45) ; S= 0.0D0
  CASE(72) ; S= 0.0D0
  CASE(71) ; S= 0.0D0
  CASE(4) ; S= 0.0D0
  CASE(12) ; S= 0.0D0
  CASE(24) ; S= 0.0D0
  CASE(44) ; S= 0.0D0
  CASE(70) ; S= xd1**2/(2.0D0*hxc**6)
  CASE(11) ; S= 0.0D0
  CASE(25) ; S= 0.0D0
  CASE(43) ; S= 0.0D0
  CASE(69) ; S= 0.0D0
  CASE(23) ; S= 0.0D0
  CASE(42) ; S= 0.0D0
  CASE(68) ; S= 0.0D0
  CASE(41) ; S= 0.0D0
  CASE(67) ; S= 0.0D0
  CASE(66) ; S= 0.0D0
  CASE(10) ; S= 0.0D0
  CASE(22) ; S= 0.0D0
  CASE(39) ; S= 0.0D0
  CASE(65) ; S= 0.0D0
  CASE(21) ; S= 0.0D0
  CASE(40) ; S= 0.0D0
  CASE(64) ; S= 0.0D0
  CASE(38) ; S= 0.0D0
  CASE(63) ; S= 0.0D0
  CASE(62) ; S= 0.0D0
  CASE(20) ; S= 0.0D0
  CASE(37) ; S= 0.0D0
  CASE(61) ; S= 0.0D0
  CASE(36) ; S= 0.0D0
  CASE(60) ; S= 0.0D0
  CASE(59) ; S= 0.0D0
  CASE(35) ; S= 0.0D0
  CASE(58) ; S= 0.0D0
  CASE(57) ; S= 0.0D0
  CASE(56) ; S= 0.0D0
  
  
  
  End select

    TL3DZ4 = s

  end function
 
 
 
  REAL FUNCTION  TL3DZ5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= hxc**(-5)
 CASE(83) ; S= zd1/hxc**6
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= yd1/hxc**6
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= xd1/hxc**6
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
   End select

   TL3DZ5 = s

  end function
 
 
 
  REAL FUNCTION TL3DZ6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
  CASE(3) ; S= 0.0D0
  CASE(9) ; S= 0.0D0
  CASE(19) ; S= 0.0D0
  CASE(34) ; S= 0.0D0
  CASE(55) ; S= 0.0D0
  CASE(83) ; S= hxc**(-6)
  CASE(2) ; S= 0.0D0
  CASE(8) ; S= 0.0D0
  CASE(18) ; S= 0.0D0
  CASE(33) ; S= 0.0D0
  CASE(54) ; S= 0.0D0
  CASE(82) ; S= 0.0D0
  CASE(7) ; S= 0.0D0
  CASE(17) ; S= 0.0D0
  CASE(32) ; S= 0.0D0
  CASE(53) ; S= 0.0D0
  CASE(81) ; S= 0.0D0
  CASE(16) ; S= 0.0D0
  CASE(31) ; S= 0.0D0
  CASE(52) ; S= 0.0D0
  CASE(80) ; S= 0.0D0
  CASE(30) ; S= 0.0D0
  CASE(51) ; S= 0.0D0
  CASE(79) ; S= 0.0D0
  CASE(50) ; S= 0.0D0
  CASE(78) ; S= 0.0D0
  CASE(77) ; S= 0.0D0
  CASE(1) ; S= 0.0D0
  CASE(6) ; S= 0.0D0
  CASE(14) ; S= 0.0D0
  CASE(27) ; S= 0.0D0
  CASE(49) ; S= 0.0D0
  CASE(76) ; S= 0.0D0
  CASE(5) ; S= 0.0D0
  CASE(15) ; S= 0.0D0
  CASE(29) ; S= 0.0D0
  CASE(48) ; S= 0.0D0
  CASE(75) ; S= 0.0D0
  CASE(13) ; S= 0.0D0
  CASE(28) ; S= 0.0D0
  CASE(47) ; S= 0.0D0
  CASE(74) ; S= 0.0D0
  CASE(26) ; S= 0.0D0
  CASE(46) ; S= 0.0D0
  CASE(73) ; S= 0.0D0
  CASE(45) ; S= 0.0D0
  CASE(72) ; S= 0.0D0
  CASE(71) ; S= 0.0D0
  CASE(4) ; S= 0.0D0
  CASE(12) ; S= 0.0D0
  CASE(24) ; S= 0.0D0
  CASE(44) ; S= 0.0D0
  CASE(70) ; S= 0.0D0
  CASE(11) ; S= 0.0D0
  CASE(25) ; S= 0.0D0
  CASE(43) ; S= 0.0D0
  CASE(69) ; S= 0.0D0
  CASE(23) ; S= 0.0D0
  CASE(42) ; S= 0.0D0
  CASE(68) ; S= 0.0D0
  CASE(41) ; S= 0.0D0
  CASE(67) ; S= 0.0D0
  CASE(66) ; S= 0.0D0
  CASE(10) ; S= 0.0D0
  CASE(22) ; S= 0.0D0
  CASE(39) ; S= 0.0D0
  CASE(65) ; S= 0.0D0
  CASE(21) ; S= 0.0D0
  CASE(40) ; S= 0.0D0
  CASE(64) ; S= 0.0D0
  CASE(38) ; S= 0.0D0
  CASE(63) ; S= 0.0D0
  CASE(62) ; S= 0.0D0
  CASE(20) ; S= 0.0D0
  CASE(37) ; S= 0.0D0
  CASE(61) ; S= 0.0D0
  CASE(36) ; S= 0.0D0
  CASE(60) ; S= 0.0D0
  CASE(59) ; S= 0.0D0
  CASE(35) ; S= 0.0D0
  CASE(58) ; S= 0.0D0
  CASE(57) ; S= 0.0D0
  CASE(56) ; S= 0.0D0
  
  
  
  End select

    TL3DZ6= s

  end function
 
 
 
  REAL FUNCTION TL3DY(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 1/hxc
 CASE(8) ; S= zd1/hxc**2
 CASE(18) ; S= zd1**2/(2.0D0*hxc**3)
 CASE(33) ; S= zd1**3/(6.0D0*hxc**4)
 CASE(54) ; S= zd1**4/(24.0D0*hxc**5)
 CASE(82) ; S= zd1**5/(120.0D0*hxc**6)
 CASE(7) ; S= yd1/hxc**2
 CASE(17) ; S= (yd1*zd1)/hxc**3
 CASE(32) ; S= (yd1*zd1**2)/(2.0D0*hxc**4)
 CASE(53) ; S= (yd1*zd1**3)/(6.0D0*hxc**5)
 CASE(81) ; S= (yd1*zd1**4)/(24.0D0*hxc**6)
 CASE(16) ; S= yd1**2/(2.0D0*hxc**3)
 CASE(31) ; S= (yd1**2*zd1)/(2.0D0*hxc**4)
 CASE(52) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**5)
 CASE(80) ; S= (yd1**2*zd1**3)/(12.0D0*hxc**6)
 CASE(30) ; S= yd1**3/(6.0D0*hxc**4)
 CASE(51) ; S= (yd1**3*zd1)/(6.0D0*hxc**5)
 CASE(79) ; S= (yd1**3*zd1**2)/(12.0D0*hxc**6)
 CASE(50) ; S= yd1**4/(24.0D0*hxc**5)
 CASE(78) ; S= (yd1**4*zd1)/(24.0D0*hxc**6)
 CASE(77) ; S= yd1**5/(120.0D0*hxc**6)
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= xd1/hxc**2
 CASE(15) ; S= (xd1*zd1)/hxc**3
 CASE(29) ; S= (xd1*zd1**2)/(2.0D0*hxc**4)
 CASE(48) ; S= (xd1*zd1**3)/(6.0D0*hxc**5)
 CASE(75) ; S= (xd1*zd1**4)/(24.0D0*hxc**6)
 CASE(13) ; S= (xd1*yd1)/hxc**3
 CASE(28) ; S= (xd1*yd1*zd1)/hxc**4
 CASE(47) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**5)
 CASE(74) ; S= (xd1*yd1*zd1**3)/(6.0D0*hxc**6)
 CASE(26) ; S= (xd1*yd1**2)/(2.0D0*hxc**4)
 CASE(46) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**5)
 CASE(73) ; S= (xd1*yd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(45) ; S= (xd1*yd1**3)/(6.0D0*hxc**5)
 CASE(72) ; S= (xd1*yd1**3*zd1)/(6.0D0*hxc**6)
 CASE(71) ; S= (xd1*yd1**4)/(24.0D0*hxc**6)
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= xd1**2/(2.0D0*hxc**3)
 CASE(25) ; S= (xd1**2*zd1)/(2.0D0*hxc**4)
 CASE(43) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**5)
 CASE(69) ; S= (xd1**2*zd1**3)/(12.0D0*hxc**6)
 CASE(23) ; S= (xd1**2*yd1)/(2.0D0*hxc**4)
 CASE(42) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**5)
 CASE(68) ; S= (xd1**2*yd1*zd1**2)/(4.0D0*hxc**6)
 CASE(41) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**5)
 CASE(67) ; S= (xd1**2*yd1**2*zd1)/(4.0D0*hxc**6)
 CASE(66) ; S= (xd1**2*yd1**3)/(12.0D0*hxc**6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= xd1**3/(6.0D0*hxc**4)
 CASE(40) ; S= (xd1**3*zd1)/(6.0D0*hxc**5)
 CASE(64) ; S= (xd1**3*zd1**2)/(12.0D0*hxc**6)
 CASE(38) ; S= (xd1**3*yd1)/(6.0D0*hxc**5)
 CASE(63) ; S= (xd1**3*yd1*zd1)/(6.0D0*hxc**6)
 CASE(62) ; S= (xd1**3*yd1**2)/(12.0D0*hxc**6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= xd1**4/(24.0D0*hxc**5)
 CASE(60) ; S= (xd1**4*zd1)/(24.0D0*hxc**6)
 CASE(59) ; S= (xd1**4*yd1)/(24.0D0*hxc**6)
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= xd1**5/(120.0D0*hxc**6)
 CASE(56) ; S= 0.0D0
 
 
 
   End select

    TL3DY = s

  end function
 
 
 
  REAL FUNCTION TL3DYZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
  CASE(3) ; S= 0.0D0
  CASE(9) ; S= 0.0D0
  CASE(19) ; S= 0.0D0
  CASE(34) ; S= 0.0D0
  CASE(55) ; S= 0.0D0
  CASE(83) ; S= 0.0D0
  CASE(2) ; S= 0.0D0
  CASE(8) ; S= hxc**(-2)
  CASE(18) ; S= zd1/hxc**3
  CASE(33) ; S= zd1**2/(2.0D0*hxc**4)
  CASE(54) ; S= zd1**3/(6.0D0*hxc**5)
  CASE(82) ; S= zd1**4/(24.0D0*hxc**6)
  CASE(7) ; S= 0.0D0
  CASE(17) ; S= yd1/hxc**3
  CASE(32) ; S= (yd1*zd1)/hxc**4
  CASE(53) ; S= (yd1*zd1**2)/(2.0D0*hxc**5)
  CASE(81) ; S= (yd1*zd1**3)/(6.0D0*hxc**6)
  CASE(16) ; S= 0.0D0
  CASE(31) ; S= yd1**2/(2.0D0*hxc**4)
  CASE(52) ; S= (yd1**2*zd1)/(2.0D0*hxc**5)
  CASE(80) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**6)
  CASE(30) ; S= 0.0D0
  CASE(51) ; S= yd1**3/(6.0D0*hxc**5)
  CASE(79) ; S= (yd1**3*zd1)/(6.0D0*hxc**6)
  CASE(50) ; S= 0.0D0
  CASE(78) ; S= yd1**4/(24.0D0*hxc**6)
  CASE(77) ; S= 0.0D0
  CASE(1) ; S= 0.0D0
  CASE(6) ; S= 0.0D0
  CASE(14) ; S= 0.0D0
  CASE(27) ; S= 0.0D0
  CASE(49) ; S= 0.0D0
  CASE(76) ; S= 0.0D0
  CASE(5) ; S= 0.0D0
  CASE(15) ; S= xd1/hxc**3
  CASE(29) ; S= (xd1*zd1)/hxc**4
  CASE(48) ; S= (xd1*zd1**2)/(2.0D0*hxc**5)
  CASE(75) ; S= (xd1*zd1**3)/(6.0D0*hxc**6)
  CASE(13) ; S= 0.0D0
  CASE(28) ; S= (xd1*yd1)/hxc**4
  CASE(47) ; S= (xd1*yd1*zd1)/hxc**5
  CASE(74) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**6)
  CASE(26) ; S= 0.0D0
  CASE(46) ; S= (xd1*yd1**2)/(2.0D0*hxc**5)
  CASE(73) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**6)
  CASE(45) ; S= 0.0D0
  CASE(72) ; S= (xd1*yd1**3)/(6.0D0*hxc**6)
  CASE(71) ; S= 0.0D0
  CASE(4) ; S= 0.0D0
  CASE(12) ; S= 0.0D0
  CASE(24) ; S= 0.0D0
  CASE(44) ; S= 0.0D0
  CASE(70) ; S= 0.0D0
  CASE(11) ; S= 0.0D0
  CASE(25) ; S= xd1**2/(2.0D0*hxc**4)
  CASE(43) ; S= (xd1**2*zd1)/(2.0D0*hxc**5)
  CASE(69) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**6)
  CASE(23) ; S= 0.0D0
  CASE(42) ; S= (xd1**2*yd1)/(2.0D0*hxc**5)
  CASE(68) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**6)
  CASE(41) ; S= 0.0D0
  CASE(67) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**6)
  CASE(66) ; S= 0.0D0
  CASE(10) ; S= 0.0D0
  CASE(22) ; S= 0.0D0
  CASE(39) ; S= 0.0D0
  CASE(65) ; S= 0.0D0
  CASE(21) ; S= 0.0D0
  CASE(40) ; S= xd1**3/(6.0D0*hxc**5)
  CASE(64) ; S= (xd1**3*zd1)/(6.0D0*hxc**6)
  CASE(38) ; S= 0.0D0
  CASE(63) ; S= (xd1**3*yd1)/(6.0D0*hxc**6)
  CASE(62) ; S= 0.0D0
  CASE(20) ; S= 0.0D0
  CASE(37) ; S= 0.0D0
  CASE(61) ; S= 0.0D0
  CASE(36) ; S= 0.0D0
  CASE(60) ; S= xd1**4/(24.0D0*hxc**6)
  CASE(59) ; S= 0.0D0
  CASE(35) ; S= 0.0D0
  CASE(58) ; S= 0.0D0
  CASE(57) ; S= 0.0D0
  CASE(56) ; S= 0.0D0
  
  
  
  
   End select

   TL3DYZ = s

  end function
 
 
 
  REAL FUNCTION TL3DYZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
  CASE(3) ; S= 0.0D0
  CASE(9) ; S= 0.0D0
  CASE(19) ; S= 0.0D0
  CASE(34) ; S= 0.0D0
  CASE(55) ; S= 0.0D0
  CASE(83) ; S= 0.0D0
  CASE(2) ; S= 0.0D0
  CASE(8) ; S= 0.0D0
  CASE(18) ; S= hxc**(-3)
  CASE(33) ; S= zd1/hxc**4
  CASE(54) ; S= zd1**2/(2.0D0*hxc**5)
  CASE(82) ; S= zd1**3/(6.0D0*hxc**6)
  CASE(7) ; S= 0.0D0
  CASE(17) ; S= 0.0D0
  CASE(32) ; S= yd1/hxc**4
  CASE(53) ; S= (yd1*zd1)/hxc**5
  CASE(81) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
  CASE(16) ; S= 0.0D0
  CASE(31) ; S= 0.0D0
  CASE(52) ; S= yd1**2/(2.0D0*hxc**5)
  CASE(80) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
  CASE(30) ; S= 0.0D0
  CASE(51) ; S= 0.0D0
  CASE(79) ; S= yd1**3/(6.0D0*hxc**6)
  CASE(50) ; S= 0.0D0
  CASE(78) ; S= 0.0D0
  CASE(77) ; S= 0.0D0
  CASE(1) ; S= 0.0D0
  CASE(6) ; S= 0.0D0
  CASE(14) ; S= 0.0D0
  CASE(27) ; S= 0.0D0
  CASE(49) ; S= 0.0D0
  CASE(76) ; S= 0.0D0
  CASE(5) ; S= 0.0D0
  CASE(15) ; S= 0.0D0
  CASE(29) ; S= xd1/hxc**4
  CASE(48) ; S= (xd1*zd1)/hxc**5
  CASE(75) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
  CASE(13) ; S= 0.0D0
  CASE(28) ; S= 0.0D0
  CASE(47) ; S= (xd1*yd1)/hxc**5
  CASE(74) ; S= (xd1*yd1*zd1)/hxc**6
  CASE(26) ; S= 0.0D0
  CASE(46) ; S= 0.0D0
  CASE(73) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
  CASE(45) ; S= 0.0D0
  CASE(72) ; S= 0.0D0
  CASE(71) ; S= 0.0D0
  CASE(4) ; S= 0.0D0
  CASE(12) ; S= 0.0D0
  CASE(24) ; S= 0.0D0
  CASE(44) ; S= 0.0D0
  CASE(70) ; S= 0.0D0
  CASE(11) ; S= 0.0D0
  CASE(25) ; S= 0.0D0
  CASE(43) ; S= xd1**2/(2.0D0*hxc**5)
  CASE(69) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
  CASE(23) ; S= 0.0D0
  CASE(42) ; S= 0.0D0
  CASE(68) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
  CASE(41) ; S= 0.0D0
  CASE(67) ; S= 0.0D0
  CASE(66) ; S= 0.0D0
  CASE(10) ; S= 0.0D0
  CASE(22) ; S= 0.0D0
  CASE(39) ; S= 0.0D0
  CASE(65) ; S= 0.0D0
  CASE(21) ; S= 0.0D0
  CASE(40) ; S= 0.0D0
  CASE(64) ; S= xd1**3/(6.0D0*hxc**6)
  CASE(38) ; S= 0.0D0
  CASE(63) ; S= 0.0D0
  CASE(62) ; S= 0.0D0
  CASE(20) ; S= 0.0D0
  CASE(37) ; S= 0.0D0
  CASE(61) ; S= 0.0D0
  CASE(36) ; S= 0.0D0
  CASE(60) ; S= 0.0D0
  CASE(59) ; S= 0.0D0
  CASE(35) ; S= 0.0D0
  CASE(58) ; S= 0.0D0
  CASE(57) ; S= 0.0D0
  CASE(56) ; S= 0.0D0
  
  
  
  End select

    TL3DYZ2 = s

  end function
 
 
 
  REAL FUNCTION TL3DYZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= hxc**(-4)
 CASE(54) ; S= zd1/hxc**5
 CASE(82) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= yd1/hxc**5
 CASE(81) ; S= (yd1*zd1)/hxc**6
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= xd1/hxc**5
 CASE(75) ; S= (xd1*zd1)/hxc**6
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= (xd1*yd1)/hxc**6
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

   TL3DYZ3 = s

  end function
 
 
 
  REAL FUNCTION TL3DYZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= hxc**(-5)
 CASE(82) ; S= zd1/hxc**6
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= yd1/hxc**6
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= xd1/hxc**6
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DYZ4 = s

  end function
 
 
 
  REAL FUNCTION TL3DYZ5 (XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= hxc**(-6)
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DYZ5 = s

  end function
 
 
 
  REAL FUNCTION TL3DY2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= hxc**(-2)
 CASE(17) ; S= zd1/hxc**3
 CASE(32) ; S= zd1**2/(2.0D0*hxc**4)
 CASE(53) ; S= zd1**3/(6.0D0*hxc**5)
 CASE(81) ; S= zd1**4/(24.0D0*hxc**6)
 CASE(16) ; S= yd1/hxc**3
 CASE(31) ; S= (yd1*zd1)/hxc**4
 CASE(52) ; S= (yd1*zd1**2)/(2.0D0*hxc**5)
 CASE(80) ; S= (yd1*zd1**3)/(6.0D0*hxc**6)
 CASE(30) ; S= yd1**2/(2.0D0*hxc**4)
 CASE(51) ; S= (yd1**2*zd1)/(2.0D0*hxc**5)
 CASE(79) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(50) ; S= yd1**3/(6.0D0*hxc**5)
 CASE(78) ; S= (yd1**3*zd1)/(6.0D0*hxc**6)
 CASE(77) ; S= yd1**4/(24.0D0*hxc**6)
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= xd1/hxc**3
 CASE(28) ; S= (xd1*zd1)/hxc**4
 CASE(47) ; S= (xd1*zd1**2)/(2.0D0*hxc**5)
 CASE(74) ; S= (xd1*zd1**3)/(6.0D0*hxc**6)
 CASE(26) ; S= (xd1*yd1)/hxc**4
 CASE(46) ; S= (xd1*yd1*zd1)/hxc**5
 CASE(73) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(45) ; S= (xd1*yd1**2)/(2.0D0*hxc**5)
 CASE(72) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(71) ; S= (xd1*yd1**3)/(6.0D0*hxc**6)
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= xd1**2/(2.0D0*hxc**4)
 CASE(42) ; S= (xd1**2*zd1)/(2.0D0*hxc**5)
 CASE(68) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(41) ; S= (xd1**2*yd1)/(2.0D0*hxc**5)
 CASE(67) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**6)
 CASE(66) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= xd1**3/(6.0D0*hxc**5)
 CASE(63) ; S= (xd1**3*zd1)/(6.0D0*hxc**6)
 CASE(62) ; S= (xd1**3*yd1)/(6.0D0*hxc**6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= xd1**4/(24.0D0*hxc**6)
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
 
  End select

    TL3DY2 = s

  end function
 
 
 
  REAL FUNCTION TL3DY2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= hxc**(-3)
 CASE(32) ; S= zd1/hxc**4
 CASE(53) ; S= zd1**2/(2.0D0*hxc**5)
 CASE(81) ; S= zd1**3/(6.0D0*hxc**6)
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= yd1/hxc**4
 CASE(52) ; S= (yd1*zd1)/hxc**5
 CASE(80) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= yd1**2/(2.0D0*hxc**5)
 CASE(79) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= yd1**3/(6.0D0*hxc**6)
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= xd1/hxc**4
 CASE(47) ; S= (xd1*zd1)/hxc**5
 CASE(74) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= (xd1*yd1)/hxc**5
 CASE(73) ; S= (xd1*yd1*zd1)/hxc**6
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= xd1**2/(2.0D0*hxc**5)
 CASE(68) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= xd1**3/(6.0D0*hxc**6)
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DY2Z = s

  end function
 
 
 
  REAL FUNCTION TL3DY2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= hxc**(-4)
 CASE(53) ; S= zd1/hxc**5
 CASE(81) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= yd1/hxc**5
 CASE(80) ; S= (yd1*zd1)/hxc**6
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= xd1/hxc**5
 CASE(74) ; S= (xd1*zd1)/hxc**6
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= (xd1*yd1)/hxc**6
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DY2Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DY2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= hxc**(-5)
 CASE(81) ; S= zd1/hxc**6
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= yd1/hxc**6
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= xd1/hxc**6
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DY2Z3 = s

  end function
 
 
 
  REAL FUNCTION TL3DY2Z4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= hxc**(-6)
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DY2Z4 = s

  end function
 
 
 
  REAL FUNCTION TL3DY3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= hxc**(-3)
 CASE(31) ; S= zd1/hxc**4
 CASE(52) ; S= zd1**2/(2.0D0*hxc**5)
 CASE(80) ; S= zd1**3/(6.0D0*hxc**6)
 CASE(30) ; S= yd1/hxc**4
 CASE(51) ; S= (yd1*zd1)/hxc**5
 CASE(79) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(50) ; S= yd1**2/(2.0D0*hxc**5)
 CASE(78) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(77) ; S= yd1**3/(6.0D0*hxc**6)
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= xd1/hxc**4
 CASE(46) ; S= (xd1*zd1)/hxc**5
 CASE(73) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
 CASE(45) ; S= (xd1*yd1)/hxc**5
 CASE(72) ; S= (xd1*yd1*zd1)/hxc**6
 CASE(71) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= xd1**2/(2.0D0*hxc**5)
 CASE(67) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
 CASE(66) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= xd1**3/(6.0D0*hxc**6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DY3 = s

  end function
 
 
 
  REAL FUNCTION TL3DY3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= hxc**(-4)
 CASE(52) ; S= zd1/hxc**5
 CASE(80) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= yd1/hxc**5
 CASE(79) ; S= (yd1*zd1)/hxc**6
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= xd1/hxc**5
 CASE(73) ; S= (xd1*zd1)/hxc**6
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= (xd1*yd1)/hxc**6
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DY3Z = s

  end function
 
 
 
  REAL FUNCTION TL3DY3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= hxc**(-5)
 CASE(80) ; S= zd1/hxc**6
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= yd1/hxc**6
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= xd1/hxc**6
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DY3Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DY3Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= hxc**(-6)
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DY3Z3 = s

  end function
 
 
 
  REAL FUNCTION TL3DY4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= hxc**(-4)
 CASE(51) ; S= zd1/hxc**5
 CASE(79) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(50) ; S= yd1/hxc**5
 CASE(78) ; S= (yd1*zd1)/hxc**6
 CASE(77) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= xd1/hxc**5
 CASE(72) ; S= (xd1*zd1)/hxc**6
 CASE(71) ; S= (xd1*yd1)/hxc**6
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DY4 = s

  end function
 
 
 
  REAL FUNCTION TL3DY4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= hxc**(-5)
 CASE(79) ; S= zd1/hxc**6
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= yd1/hxc**6
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= xd1/hxc**6
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DY4Z = s

  end function
 
 
 
  REAL FUNCTION TL3DY4Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= hxc**(-6)
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DY4Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DY5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= hxc**(-5)
 CASE(78) ; S= zd1/hxc**6
 CASE(77) ; S= yd1/hxc**6
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= xd1/hxc**6
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DY5 = s

  end function
 
 
 
  REAL FUNCTION TL3DY5Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= hxc**(-6)
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DY5Z = s

  end function
 
 
 
  REAL FUNCTION TL3DY6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= hxc**(-6)
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DY6 = s

  end function
 
 
 
  REAL FUNCTION TL3DX(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 1/hxc
 CASE(6) ; S= zd1/hxc**2
 CASE(14) ; S= zd1**2/(2.0D0*hxc**3)
 CASE(27) ; S= zd1**3/(6.0D0*hxc**4)
 CASE(49) ; S= zd1**4/(24.0D0*hxc**5)
 CASE(76) ; S= zd1**5/(120.0D0*hxc**6)
 CASE(5) ; S= yd1/hxc**2
 CASE(15) ; S= (yd1*zd1)/hxc**3
 CASE(29) ; S= (yd1*zd1**2)/(2.0D0*hxc**4)
 CASE(48) ; S= (yd1*zd1**3)/(6.0D0*hxc**5)
 CASE(75) ; S= (yd1*zd1**4)/(24.0D0*hxc**6)
 CASE(13) ; S= yd1**2/(2.0D0*hxc**3)
 CASE(28) ; S= (yd1**2*zd1)/(2.0D0*hxc**4)
 CASE(47) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**5)
 CASE(74) ; S= (yd1**2*zd1**3)/(12.0D0*hxc**6)
 CASE(26) ; S= yd1**3/(6.0D0*hxc**4)
 CASE(46) ; S= (yd1**3*zd1)/(6.0D0*hxc**5)
 CASE(73) ; S= (yd1**3*zd1**2)/(12.0D0*hxc**6)
 CASE(45) ; S= yd1**4/(24.0D0*hxc**5)
 CASE(72) ; S= (yd1**4*zd1)/(24.0D0*hxc**6)
 CASE(71) ; S= yd1**5/(120.0D0*hxc**6)
 CASE(4) ; S= xd1/hxc**2
 CASE(12) ; S= (xd1*zd1)/hxc**3
 CASE(24) ; S= (xd1*zd1**2)/(2.0D0*hxc**4)
 CASE(44) ; S= (xd1*zd1**3)/(6.0D0*hxc**5)
 CASE(70) ; S= (xd1*zd1**4)/(24.0D0*hxc**6)
 CASE(11) ; S= (xd1*yd1)/hxc**3
 CASE(25) ; S= (xd1*yd1*zd1)/hxc**4
 CASE(43) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**5)
 CASE(69) ; S= (xd1*yd1*zd1**3)/(6.0D0*hxc**6)
 CASE(23) ; S= (xd1*yd1**2)/(2.0D0*hxc**4)
 CASE(42) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**5)
 CASE(68) ; S= (xd1*yd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(41) ; S= (xd1*yd1**3)/(6.0D0*hxc**5)
 CASE(67) ; S= (xd1*yd1**3*zd1)/(6.0D0*hxc**6)
 CASE(66) ; S= (xd1*yd1**4)/(24.0D0*hxc**6)
 CASE(10) ; S= xd1**2/(2.0D0*hxc**3)
 CASE(22) ; S= (xd1**2*zd1)/(2.0D0*hxc**4)
 CASE(39) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**5)
 CASE(65) ; S= (xd1**2*zd1**3)/(12.0D0*hxc**6)
 CASE(21) ; S= (xd1**2*yd1)/(2.0D0*hxc**4)
 CASE(40) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**5)
 CASE(64) ; S= (xd1**2*yd1*zd1**2)/(4.0D0*hxc**6)
 CASE(38) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**5)
 CASE(63) ; S= (xd1**2*yd1**2*zd1)/(4.0D0*hxc**6)
 CASE(62) ; S= (xd1**2*yd1**3)/(12.0D0*hxc**6)
 CASE(20) ; S= xd1**3/(6.0D0*hxc**4)
 CASE(37) ; S= (xd1**3*zd1)/(6.0D0*hxc**5)
 CASE(61) ; S= (xd1**3*zd1**2)/(12.0D0*hxc**6)
 CASE(36) ; S= (xd1**3*yd1)/(6.0D0*hxc**5)
 CASE(60) ; S= (xd1**3*yd1*zd1)/(6.0D0*hxc**6)
 CASE(59) ; S= (xd1**3*yd1**2)/(12.0D0*hxc**6)
 CASE(35) ; S= xd1**4/(24.0D0*hxc**5)
 CASE(58) ; S= (xd1**4*zd1)/(24.0D0*hxc**6)
 CASE(57) ; S= (xd1**4*yd1)/(24.0D0*hxc**6)
 CASE(56) ; S= xd1**5/(120.0D0*hxc**6)
 
  End select

   TL3DX = s

  end function
 
 
 
  REAL FUNCTION TL3DXZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= hxc**(-2)
 CASE(14) ; S= zd1/hxc**3
 CASE(27) ; S= zd1**2/(2.0D0*hxc**4)
 CASE(49) ; S= zd1**3/(6.0D0*hxc**5)
 CASE(76) ; S= zd1**4/(24.0D0*hxc**6)
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= yd1/hxc**3
 CASE(29) ; S= (yd1*zd1)/hxc**4
 CASE(48) ; S= (yd1*zd1**2)/(2.0D0*hxc**5)
 CASE(75) ; S= (yd1*zd1**3)/(6.0D0*hxc**6)
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= yd1**2/(2.0D0*hxc**4)
 CASE(47) ; S= (yd1**2*zd1)/(2.0D0*hxc**5)
 CASE(74) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= yd1**3/(6.0D0*hxc**5)
 CASE(73) ; S= (yd1**3*zd1)/(6.0D0*hxc**6)
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= yd1**4/(24.0D0*hxc**6)
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= xd1/hxc**3
 CASE(24) ; S= (xd1*zd1)/hxc**4
 CASE(44) ; S= (xd1*zd1**2)/(2.0D0*hxc**5)
 CASE(70) ; S= (xd1*zd1**3)/(6.0D0*hxc**6)
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= (xd1*yd1)/hxc**4
 CASE(43) ; S= (xd1*yd1*zd1)/hxc**5
 CASE(69) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= (xd1*yd1**2)/(2.0D0*hxc**5)
 CASE(68) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= (xd1*yd1**3)/(6.0D0*hxc**6)
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= xd1**2/(2.0D0*hxc**4)
 CASE(39) ; S= (xd1**2*zd1)/(2.0D0*hxc**5)
 CASE(65) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= (xd1**2*yd1)/(2.0D0*hxc**5)
 CASE(64) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**6)
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**6)
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= xd1**3/(6.0D0*hxc**5)
 CASE(61) ; S= (xd1**3*zd1)/(6.0D0*hxc**6)
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= (xd1**3*yd1)/(6.0D0*hxc**6)
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= xd1**4/(24.0D0*hxc**6)
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXZ = s

  end function
 
 
 
  REAL FUNCTION TL3DXZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= hxc**(-3)
 CASE(27) ; S= zd1/hxc**4
 CASE(49) ; S= zd1**2/(2.0D0*hxc**5)
 CASE(76) ; S= zd1**3/(6.0D0*hxc**6)
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= yd1/hxc**4
 CASE(48) ; S= (yd1*zd1)/hxc**5
 CASE(75) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= yd1**2/(2.0D0*hxc**5)
 CASE(74) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= yd1**3/(6.0D0*hxc**6)
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= xd1/hxc**4
 CASE(44) ; S= (xd1*zd1)/hxc**5
 CASE(70) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= (xd1*yd1)/hxc**5
 CASE(69) ; S= (xd1*yd1*zd1)/hxc**6
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= xd1**2/(2.0D0*hxc**5)
 CASE(65) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= xd1**3/(6.0D0*hxc**6)
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXZ2 = s

  end function
 
 
 
  REAL FUNCTION TL3DXZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= hxc**(-4)
 CASE(49) ; S= zd1/hxc**5
 CASE(76) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= yd1/hxc**5
 CASE(75) ; S= (yd1*zd1)/hxc**6
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= xd1/hxc**5
 CASE(70) ; S= (xd1*zd1)/hxc**6
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= (xd1*yd1)/hxc**6
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXZ3 = s

  end function
 
 
 
  REAL FUNCTION TL3DXZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= hxc**(-5)
 CASE(76) ; S= zd1/hxc**6
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= yd1/hxc**6
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= xd1/hxc**6
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXZ4 = s

  end function
 
 
 
  REAL FUNCTION TL3DXZ5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= hxc**(-6)
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DXZ5 = s

  end function
 
 
 
  REAL FUNCTION TL3DXY(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= hxc**(-2)
 CASE(15) ; S= zd1/hxc**3
 CASE(29) ; S= zd1**2/(2.0D0*hxc**4)
 CASE(48) ; S= zd1**3/(6.0D0*hxc**5)
 CASE(75) ; S= zd1**4/(24.0D0*hxc**6)
 CASE(13) ; S= yd1/hxc**3
 CASE(28) ; S= (yd1*zd1)/hxc**4
 CASE(47) ; S= (yd1*zd1**2)/(2.0D0*hxc**5)
 CASE(74) ; S= (yd1*zd1**3)/(6.0D0*hxc**6)
 CASE(26) ; S= yd1**2/(2.0D0*hxc**4)
 CASE(46) ; S= (yd1**2*zd1)/(2.0D0*hxc**5)
 CASE(73) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(45) ; S= yd1**3/(6.0D0*hxc**5)
 CASE(72) ; S= (yd1**3*zd1)/(6.0D0*hxc**6)
 CASE(71) ; S= yd1**4/(24.0D0*hxc**6)
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= xd1/hxc**3
 CASE(25) ; S= (xd1*zd1)/hxc**4
 CASE(43) ; S= (xd1*zd1**2)/(2.0D0*hxc**5)
 CASE(69) ; S= (xd1*zd1**3)/(6.0D0*hxc**6)
 CASE(23) ; S= (xd1*yd1)/hxc**4
 CASE(42) ; S= (xd1*yd1*zd1)/hxc**5
 CASE(68) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(41) ; S= (xd1*yd1**2)/(2.0D0*hxc**5)
 CASE(67) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(66) ; S= (xd1*yd1**3)/(6.0D0*hxc**6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= xd1**2/(2.0D0*hxc**4)
 CASE(40) ; S= (xd1**2*zd1)/(2.0D0*hxc**5)
 CASE(64) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(38) ; S= (xd1**2*yd1)/(2.0D0*hxc**5)
 CASE(63) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**6)
 CASE(62) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= xd1**3/(6.0D0*hxc**5)
 CASE(60) ; S= (xd1**3*zd1)/(6.0D0*hxc**6)
 CASE(59) ; S= (xd1**3*yd1)/(6.0D0*hxc**6)
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= xd1**4/(24.0D0*hxc**6)
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DXY = s

  end function
 
 
 
  REAL FUNCTION TL3DXYZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= hxc**(-3)
 CASE(29) ; S= zd1/hxc**4
 CASE(48) ; S= zd1**2/(2.0D0*hxc**5)
 CASE(75) ; S= zd1**3/(6.0D0*hxc**6)
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= yd1/hxc**4
 CASE(47) ; S= (yd1*zd1)/hxc**5
 CASE(74) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= yd1**2/(2.0D0*hxc**5)
 CASE(73) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= yd1**3/(6.0D0*hxc**6)
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= xd1/hxc**4
 CASE(43) ; S= (xd1*zd1)/hxc**5
 CASE(69) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= (xd1*yd1)/hxc**5
 CASE(68) ; S= (xd1*yd1*zd1)/hxc**6
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= xd1**2/(2.0D0*hxc**5)
 CASE(64) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= xd1**3/(6.0D0*hxc**6)
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DXYZ = s

  end function
 
 
 
  REAL FUNCTION TL3DXYZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= hxc**(-4)
 CASE(48) ; S= zd1/hxc**5
 CASE(75) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= yd1/hxc**5
 CASE(74) ; S= (yd1*zd1)/hxc**6
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= xd1/hxc**5
 CASE(69) ; S= (xd1*zd1)/hxc**6
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= (xd1*yd1)/hxc**6
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXYZ2 = s

  end function
 
 
 
  REAL FUNCTION TL3DXYZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= hxc**(-5)
 CASE(75) ; S= zd1/hxc**6
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= yd1/hxc**6
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= xd1/hxc**6
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXYZ3 = s

  end function
 
 
 
  REAL FUNCTION TL3DXYZ4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= hxc**(-6)
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXYZ4 = s

  end function
 
 
 
  REAL FUNCTION TL3DXY2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= hxc**(-3)
 CASE(28) ; S= zd1/hxc**4
 CASE(47) ; S= zd1**2/(2.0D0*hxc**5)
 CASE(74) ; S= zd1**3/(6.0D0*hxc**6)
 CASE(26) ; S= yd1/hxc**4
 CASE(46) ; S= (yd1*zd1)/hxc**5
 CASE(73) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(45) ; S= yd1**2/(2.0D0*hxc**5)
 CASE(72) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(71) ; S= yd1**3/(6.0D0*hxc**6)
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= xd1/hxc**4
 CASE(42) ; S= (xd1*zd1)/hxc**5
 CASE(68) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
 CASE(41) ; S= (xd1*yd1)/hxc**5
 CASE(67) ; S= (xd1*yd1*zd1)/hxc**6
 CASE(66) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= xd1**2/(2.0D0*hxc**5)
 CASE(63) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
 CASE(62) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= xd1**3/(6.0D0*hxc**6)
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXY2 = s

  end function
 
 
 
  REAL FUNCTION TL3DXY2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= hxc**(-4)
 CASE(47) ; S= zd1/hxc**5
 CASE(74) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= yd1/hxc**5
 CASE(73) ; S= (yd1*zd1)/hxc**6
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= xd1/hxc**5
 CASE(68) ; S= (xd1*zd1)/hxc**6
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= (xd1*yd1)/hxc**6
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DXY2Z = s

  end function
 
 
 
  REAL FUNCTION TL3DXY2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= hxc**(-5)
 CASE(74) ; S= zd1/hxc**6
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= yd1/hxc**6
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= xd1/hxc**6
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DXY2Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DXY2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= hxc**(-6)
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXY2Z3 = s

  end function
 
 
 
  REAL FUNCTION TL3DXY3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= hxc**(-4)
 CASE(46) ; S= zd1/hxc**5
 CASE(73) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(45) ; S= yd1/hxc**5
 CASE(72) ; S= (yd1*zd1)/hxc**6
 CASE(71) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= xd1/hxc**5
 CASE(67) ; S= (xd1*zd1)/hxc**6
 CASE(66) ; S= (xd1*yd1)/hxc**6
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DXY3 = s

  end function
 
 
 
  REAL FUNCTION TL3DXY3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= hxc**(-5)
 CASE(73) ; S= zd1/hxc**6
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= yd1/hxc**6
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= xd1/hxc**6
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

   TL3DXY3Z = s

  end function
 
 
 
  REAL FUNCTION TL3DXY3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= hxc**(-6)
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DXY3Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DXY4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= hxc**(-5)
 CASE(72) ; S= zd1/hxc**6
 CASE(71) ; S= yd1/hxc**6
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= xd1/hxc**6
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DXY4 = s

  end function
 
 
 
  REAL FUNCTION TL3DXY4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= hxc**(-6)
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXY4Z = s

  end function
 
 
 
  REAL FUNCTION TL3DXY5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= hxc**(-6)
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DXY5 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= hxc**(-2)
 CASE(12) ; S= zd1/hxc**3
 CASE(24) ; S= zd1**2/(2.0D0*hxc**4)
 CASE(44) ; S= zd1**3/(6.0D0*hxc**5)
 CASE(70) ; S= zd1**4/(24.0D0*hxc**6)
 CASE(11) ; S= yd1/hxc**3
 CASE(25) ; S= (yd1*zd1)/hxc**4
 CASE(43) ; S= (yd1*zd1**2)/(2.0D0*hxc**5)
 CASE(69) ; S= (yd1*zd1**3)/(6.0D0*hxc**6)
 CASE(23) ; S= yd1**2/(2.0D0*hxc**4)
 CASE(42) ; S= (yd1**2*zd1)/(2.0D0*hxc**5)
 CASE(68) ; S= (yd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(41) ; S= yd1**3/(6.0D0*hxc**5)
 CASE(67) ; S= (yd1**3*zd1)/(6.0D0*hxc**6)
 CASE(66) ; S= yd1**4/(24.0D0*hxc**6)
 CASE(10) ; S= xd1/hxc**3
 CASE(22) ; S= (xd1*zd1)/hxc**4
 CASE(39) ; S= (xd1*zd1**2)/(2.0D0*hxc**5)
 CASE(65) ; S= (xd1*zd1**3)/(6.0D0*hxc**6)
 CASE(21) ; S= (xd1*yd1)/hxc**4
 CASE(40) ; S= (xd1*yd1*zd1)/hxc**5
 CASE(64) ; S= (xd1*yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(38) ; S= (xd1*yd1**2)/(2.0D0*hxc**5)
 CASE(63) ; S= (xd1*yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(62) ; S= (xd1*yd1**3)/(6.0D0*hxc**6)
 CASE(20) ; S= xd1**2/(2.0D0*hxc**4)
 CASE(37) ; S= (xd1**2*zd1)/(2.0D0*hxc**5)
 CASE(61) ; S= (xd1**2*zd1**2)/(4.0D0*hxc**6)
 CASE(36) ; S= (xd1**2*yd1)/(2.0D0*hxc**5)
 CASE(60) ; S= (xd1**2*yd1*zd1)/(2.0D0*hxc**6)
 CASE(59) ; S= (xd1**2*yd1**2)/(4.0D0*hxc**6)
 CASE(35) ; S= xd1**3/(6.0D0*hxc**5)
 CASE(58) ; S= (xd1**3*zd1)/(6.0D0*hxc**6)
 CASE(57) ; S= (xd1**3*yd1)/(6.0D0*hxc**6)
 CASE(56) ; S= xd1**4/(24.0D0*hxc**6)
 
 
  End select

    TL3DX2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= hxc**(-3)
 CASE(24) ; S= zd1/hxc**4
 CASE(44) ; S= zd1**2/(2.0D0*hxc**5)
 CASE(70) ; S= zd1**3/(6.0D0*hxc**6)
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= yd1/hxc**4
 CASE(43) ; S= (yd1*zd1)/hxc**5
 CASE(69) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= yd1**2/(2.0D0*hxc**5)
 CASE(68) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= yd1**3/(6.0D0*hxc**6)
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= xd1/hxc**4
 CASE(39) ; S= (xd1*zd1)/hxc**5
 CASE(65) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= (xd1*yd1)/hxc**5
 CASE(64) ; S= (xd1*yd1*zd1)/hxc**6
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= xd1**2/(2.0D0*hxc**5)
 CASE(61) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= xd1**3/(6.0D0*hxc**6)
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DX2Z = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= hxc**(-4)
 CASE(44) ; S= zd1/hxc**5
 CASE(70) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= yd1/hxc**5
 CASE(69) ; S= (yd1*zd1)/hxc**6
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= xd1/hxc**5
 CASE(65) ; S= (xd1*zd1)/hxc**6
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= (xd1*yd1)/hxc**6
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
  End select

    TL3DX2Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= hxc**(-5)
 CASE(70) ; S= zd1/hxc**6
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= yd1/hxc**6
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= xd1/hxc**6
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX2Z3 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Z4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= hxc**(-6)
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DX2Z4 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= hxc**(-3)
 CASE(25) ; S= zd1/hxc**4
 CASE(43) ; S= zd1**2/(2.0D0*hxc**5)
 CASE(69) ; S= zd1**3/(6.0D0*hxc**6)
 CASE(23) ; S= yd1/hxc**4
 CASE(42) ; S= (yd1*zd1)/hxc**5
 CASE(68) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(41) ; S= yd1**2/(2.0D0*hxc**5)
 CASE(67) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(66) ; S= yd1**3/(6.0D0*hxc**6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= xd1/hxc**4
 CASE(40) ; S= (xd1*zd1)/hxc**5
 CASE(64) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
 CASE(38) ; S= (xd1*yd1)/hxc**5
 CASE(63) ; S= (xd1*yd1*zd1)/hxc**6
 CASE(62) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= xd1**2/(2.0D0*hxc**5)
 CASE(60) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
 CASE(59) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= xd1**3/(6.0D0*hxc**6)
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX2Y = s

  end function
 
 
 
  REAL FUNCTION TL3DX2YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= hxc**(-4)
 CASE(43) ; S= zd1/hxc**5
 CASE(69) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= yd1/hxc**5
 CASE(68) ; S= (yd1*zd1)/hxc**6
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= xd1/hxc**5
 CASE(64) ; S= (xd1*zd1)/hxc**6
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= (xd1*yd1)/hxc**6
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DX2YZ = s

  end function
 
 
 
  REAL FUNCTION TL3DX2YZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= hxc**(-5)
 CASE(69) ; S= zd1/hxc**6
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= yd1/hxc**6
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= xd1/hxc**6
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DX2YZ2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2YZ3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= hxc**(-6)
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX2YZ3 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= hxc**(-4)
 CASE(42) ; S= zd1/hxc**5
 CASE(68) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(41) ; S= yd1/hxc**5
 CASE(67) ; S= (yd1*zd1)/hxc**6
 CASE(66) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= xd1/hxc**5
 CASE(63) ; S= (xd1*zd1)/hxc**6
 CASE(62) ; S= (xd1*yd1)/hxc**6
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DX2Y2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Y2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= hxc**(-5)
 CASE(68) ; S= zd1/hxc**6
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= yd1/hxc**6
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= xd1/hxc**6
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX2Y2Z = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Y2Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= hxc**(-6)
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX2Y2Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Y3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= hxc**(-5)
 CASE(67) ; S= zd1/hxc**6
 CASE(66) ; S= yd1/hxc**6
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= xd1/hxc**6
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX2Y3 = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Y3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= hxc**(-6)
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX2Y3Z = s

  end function
 
 
 
  REAL FUNCTION TL3DX2Y4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= hxc**(-6)
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX2Y4 = s

  end function
 
 
 
  REAL FUNCTION TL3DX3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= hxc**(-3)
 CASE(22) ; S= zd1/hxc**4
 CASE(39) ; S= zd1**2/(2.0D0*hxc**5)
 CASE(65) ; S= zd1**3/(6.0D0*hxc**6)
 CASE(21) ; S= yd1/hxc**4
 CASE(40) ; S= (yd1*zd1)/hxc**5
 CASE(64) ; S= (yd1*zd1**2)/(2.0D0*hxc**6)
 CASE(38) ; S= yd1**2/(2.0D0*hxc**5)
 CASE(63) ; S= (yd1**2*zd1)/(2.0D0*hxc**6)
 CASE(62) ; S= yd1**3/(6.0D0*hxc**6)
 CASE(20) ; S= xd1/hxc**4
 CASE(37) ; S= (xd1*zd1)/hxc**5
 CASE(61) ; S= (xd1*zd1**2)/(2.0D0*hxc**6)
 CASE(36) ; S= (xd1*yd1)/hxc**5
 CASE(60) ; S= (xd1*yd1*zd1)/hxc**6
 CASE(59) ; S= (xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(35) ; S= xd1**2/(2.0D0*hxc**5)
 CASE(58) ; S= (xd1**2*zd1)/(2.0D0*hxc**6)
 CASE(57) ; S= (xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(56) ; S= xd1**3/(6.0D0*hxc**6)
 
 
 
  End select

    TL3DX3 = s

  end function
 
 
 
  REAL FUNCTION TL3DX3Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= hxc**(-4)
 CASE(39) ; S= zd1/hxc**5
 CASE(65) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= yd1/hxc**5
 CASE(64) ; S= (yd1*zd1)/hxc**6
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= xd1/hxc**5
 CASE(61) ; S= (xd1*zd1)/hxc**6
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= (xd1*yd1)/hxc**6
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX3Z = s

  end function
 
 
 
  REAL FUNCTION TL3DX3Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= hxc**(-5)
 CASE(65) ; S= zd1/hxc**6
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= yd1/hxc**6
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= xd1/hxc**6
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX3Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX3Z3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= hxc**(-6)
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX3Z3 = s

  end function
 
 
 
  REAL FUNCTION TL3DX3Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= hxc**(-4)
 CASE(40) ; S= zd1/hxc**5
 CASE(64) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(38) ; S= yd1/hxc**5
 CASE(63) ; S= (yd1*zd1)/hxc**6
 CASE(62) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= xd1/hxc**5
 CASE(60) ; S= (xd1*zd1)/hxc**6
 CASE(59) ; S= (xd1*yd1)/hxc**6
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= xd1**2/(2.0D0*hxc**6)
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX3Y = s

  end function
 
 
 
  REAL FUNCTION TL3DX3YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= hxc**(-5)
 CASE(64) ; S= zd1/hxc**6
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= yd1/hxc**6
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= xd1/hxc**6
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DX3YZ = s

  end function
 
 
 
  REAL FUNCTION TL3DX3YZ2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= hxc**(-6)
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DX3YZ2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX3Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= hxc**(-5)
 CASE(63) ; S= zd1/hxc**6
 CASE(62) ; S= yd1/hxc**6
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= xd1/hxc**6
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX3Y2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX3Y2Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= hxc**(-6)
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX3Y2Z = s

  end function
 
 
 
  REAL FUNCTION TL3DX3Y3(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= hxc**(-6)
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX3Y3 = s

  end function
 
 
 
  REAL FUNCTION TL3DX4(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= hxc**(-4)
 CASE(37) ; S= zd1/hxc**5
 CASE(61) ; S= zd1**2/(2.0D0*hxc**6)
 CASE(36) ; S= yd1/hxc**5
 CASE(60) ; S= (yd1*zd1)/hxc**6
 CASE(59) ; S= yd1**2/(2.0D0*hxc**6)
 CASE(35) ; S= xd1/hxc**5
 CASE(58) ; S= (xd1*zd1)/hxc**6
 CASE(57) ; S= (xd1*yd1)/hxc**6
 CASE(56) ; S= xd1**2/(2.0D0*hxc**6)
 
 
  End select

    TL3DX4 = s

  end function
 
 
 
  REAL FUNCTION TL3DX4Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= hxc**(-5)
 CASE(61) ; S= zd1/hxc**6
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= yd1/hxc**6
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= xd1/hxc**6
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX4Z = s

  end function
 
 
 
  REAL FUNCTION TL3DX4Z2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= hxc**(-6)
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX4Z2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX4Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= hxc**(-5)
 CASE(60) ; S= zd1/hxc**6
 CASE(59) ; S= yd1/hxc**6
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= xd1/hxc**6
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX4Y = s

  end function
 
 
 
  REAL FUNCTION TL3DX4YZ(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= hxc**(-6)
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
 
  End select

    TL3DX4YZ = s

  end function
 
 
 
  REAL FUNCTION TL3DX4Y2(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= hxc**(-6)
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX4Y2 = s

  end function
 
 
 
  REAL FUNCTION TL3DX5(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= hxc**(-5)
 CASE(58) ; S= zd1/hxc**6
 CASE(57) ; S= yd1/hxc**6
 CASE(56) ; S= xd1/hxc**6
 
 
  End select

    TL3DX5= s

  end function
 
 
 
  REAL FUNCTION TL3DX5Z(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= hxc**(-6)
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX5Z= s

  end function
 
 
 
  REAL FUNCTION TL3DX5Y(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= hxc**(-6)
 CASE(56) ; S= 0.0D0
 
 
  End select

    TL3DX5Y = s

  end function
 
 
 
  REAL FUNCTION TL3DX6(XD1,YD1,ZD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1,ZD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(3) ; S= 0.0D0
 CASE(9) ; S= 0.0D0
 CASE(19) ; S= 0.0D0
 CASE(34) ; S= 0.0D0
 CASE(55) ; S= 0.0D0
 CASE(83) ; S= 0.0D0
 CASE(2) ; S= 0.0D0
 CASE(8) ; S= 0.0D0
 CASE(18) ; S= 0.0D0
 CASE(33) ; S= 0.0D0
 CASE(54) ; S= 0.0D0
 CASE(82) ; S= 0.0D0
 CASE(7) ; S= 0.0D0
 CASE(17) ; S= 0.0D0
 CASE(32) ; S= 0.0D0
 CASE(53) ; S= 0.0D0
 CASE(81) ; S= 0.0D0
 CASE(16) ; S= 0.0D0
 CASE(31) ; S= 0.0D0
 CASE(52) ; S= 0.0D0
 CASE(80) ; S= 0.0D0
 CASE(30) ; S= 0.0D0
 CASE(51) ; S= 0.0D0
 CASE(79) ; S= 0.0D0
 CASE(50) ; S= 0.0D0
 CASE(78) ; S= 0.0D0
 CASE(77) ; S= 0.0D0
 CASE(1) ; S= 0.0D0
 CASE(6) ; S= 0.0D0
 CASE(14) ; S= 0.0D0
 CASE(27) ; S= 0.0D0
 CASE(49) ; S= 0.0D0
 CASE(76) ; S= 0.0D0
 CASE(5) ; S= 0.0D0
 CASE(15) ; S= 0.0D0
 CASE(29) ; S= 0.0D0
 CASE(48) ; S= 0.0D0
 CASE(75) ; S= 0.0D0
 CASE(13) ; S= 0.0D0
 CASE(28) ; S= 0.0D0
 CASE(47) ; S= 0.0D0
 CASE(74) ; S= 0.0D0
 CASE(26) ; S= 0.0D0
 CASE(46) ; S= 0.0D0
 CASE(73) ; S= 0.0D0
 CASE(45) ; S= 0.0D0
 CASE(72) ; S= 0.0D0
 CASE(71) ; S= 0.0D0
 CASE(4) ; S= 0.0D0
 CASE(12) ; S= 0.0D0
 CASE(24) ; S= 0.0D0
 CASE(44) ; S= 0.0D0
 CASE(70) ; S= 0.0D0
 CASE(11) ; S= 0.0D0
 CASE(25) ; S= 0.0D0
 CASE(43) ; S= 0.0D0
 CASE(69) ; S= 0.0D0
 CASE(23) ; S= 0.0D0
 CASE(42) ; S= 0.0D0
 CASE(68) ; S= 0.0D0
 CASE(41) ; S= 0.0D0
 CASE(67) ; S= 0.0D0
 CASE(66) ; S= 0.0D0
 CASE(10) ; S= 0.0D0
 CASE(22) ; S= 0.0D0
 CASE(39) ; S= 0.0D0
 CASE(65) ; S= 0.0D0
 CASE(21) ; S= 0.0D0
 CASE(40) ; S= 0.0D0
 CASE(64) ; S= 0.0D0
 CASE(38) ; S= 0.0D0
 CASE(63) ; S= 0.0D0
 CASE(62) ; S= 0.0D0
 CASE(20) ; S= 0.0D0
 CASE(37) ; S= 0.0D0
 CASE(61) ; S= 0.0D0
 CASE(36) ; S= 0.0D0
 CASE(60) ; S= 0.0D0
 CASE(59) ; S= 0.0D0
 CASE(35) ; S= 0.0D0
 CASE(58) ; S= 0.0D0
 CASE(57) ; S= 0.0D0
 CASE(56) ; S= hxc**(-6)

 End select

    TL3DX6 = s

  end function	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
real function DF2DX(xd1,yd1,NDERIVATIVE)
   IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s






   S=0.0D0
   select case(NDERIVATIVE)
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
    End select

    DF2DX = s

  end function


 ! *****************************************************************************
  real function DF2DY(xd1,yd1,NDERIVATIVE)
     IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s



   S=0.0D0
   select case(NDERIVATIVE)

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
    End select

    DF2DY= s

  end function


 ! *****************************************************************************
  real function DF2DX2(xd1,yd1,NDERIVATIVE)
   IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   Select case(NDERIVATIVE)
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
   End select

    DF2DX2 = s

  end function



 ! *****************************************************************************
  real function DF2DY2(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   s = 0.0d0

   select case(NDERIVATIVE)

    case(5) 
      s = 2.d0
    case(8) 
      s = 2.d0*xd1
    case(9) 
      s = 6.d0*yd1
    case(12) 
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
   DF2DY2 = s

  end function


 ! *****************************************************************************
  real function DF2DXY(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s


   S=0.0D0

   select case(NDERIVATIVE)

    case( 4 ) 
      s = 1.d0
    case( 7 ) 
      s = 2.d0*xd1
    case( 8 ) 
      s = 2.d0*yd1
    case( 11 ) 
      s = 3.d0*xd1**2
    case( 12 ) 
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

   DF2DXY = s

  end function

! %%%%%%%%%%%%%%%%
 ! *****************************************************************************
  real function DF2DX3(xd1,yd1,NDERIVATIVE)
   IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s


   S=0.0D0

   select case(NDERIVATIVE)
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

    DF2DX3 = s
  end function


 ! *****************************************************************************
  real function DF2DX2Y(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s


   S=0.0D0
   select case(NDERIVATIVE)
     case(7) 
      s = 2.d0
     case(11) 
      s = 6.d0*xd1
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
   DF2DX2Y = s

  end function



 ! *****************************************************************************
  real function DF2DXY2(xd1,yd1,NDERIVATIVE)
   IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   s  = 0.0d0
   select case(NDERIVATIVE)

     case( 8  ) 
      s = 2.d0
     case( 12  ) 
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

   DF2DXY2 = s

  end function


 ! *****************************************************************************
  real function DF2DY3(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
     case(  9  ) 
      s = 6.d0
     case(  12  ) 
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

   DF2DY3 = s

  end function



 ! *****************************************************************************
  real function DF2DX4(xd1,yd1,NDERIVATIVE)
     IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

    S=0.0D0
    select case(NDERIVATIVE)
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
    End select

    DF2DX4 = s
  end function




 ! *****************************************************************************
  real function DF2DX3Y(xd1,yd1,NDERIVATIVE)
   IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0

   select case(NDERIVATIVE)
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

   DF2DX3Y = s

  end function


 ! *****************************************************************************
  real function DF2DX2Y2(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
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

   DF2DX2Y2 = s

  end function

 ! *****************************************************************************
  real function DF2DXY3(xd1,yd1,NDERIVATIVE)
     IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
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

   DF2DXY3 = s

  end function

 ! *****************************************************************************
  real function DF2DY4(xd1,yd1,NDERIVATIVE)
    IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
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

   DF2DY4 = s

  end function


! for  6th order WENO

 ! *****************************************************************************
  real function DF2DX5(xd1,yd1,NDERIVATIVE)
    IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
    case(15)
      s = 120.d0
    case(21)
      s = 720.d0*xd1
    case(22)
      s = 120.d0*yd1
   end select

   DF2DX5 = s

  end function


 ! *****************************************************************************
  real function DF2DX4Y(xd1,yd1,NDERIVATIVE)
   IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
     case(16)
      s = 24.d0
     case(22)
      s = 120.d0*xd1
     case(23)
      s = 48.d0*yd1
   end select

   DF2DX4Y = s

  end function


 ! *****************************************************************************
  real function DF2DX3Y2(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
     case(17)
      s = 12.d0
     case(23)
      s = 48.d0*xd1
     case(24)
      s = 36.d0*yd1
    end select

   DF2DX3Y2 = s

  end function


 ! *****************************************************************************
  real function DF2DX2Y3(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
    case(18)
      s = 12.d0
    case(24)
      s = 36.d0*xd1
    case(25)
      s = 48.d0*yd1
   end select

   DF2DX2Y3 = s

  end function


 ! *****************************************************************************
  real function DF2DXY4(xd1,yd1,NDERIVATIVE)
    IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
    case(19)
      s = 24.d0
    case(25)
      s = 48.d0*xd1
    case(26)
      s = 120.d0*yd1
    end select
   DF2DXY4 = s

  end function


 ! *****************************************************************************
  real function DF2DY5(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
    case(20)
      s = 120.d0
    case(26)
      s = 120.d0*xd1
   case(27)
      s = 720.d0*yd1
   end select

   DF2DY5 = s

  end function


 ! *****************************************************************************
  real function DF2DX6(xd1,yd1,NDERIVATIVE)
    IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
     case(21)
      s = 720.d0
   end select

   DF2DX6 = s

  end function


 ! *****************************************************************************
  real function DF2DX5Y(xd1,yd1,NDERIVATIVE)
    IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
     case(22)
      s = 120.d0
   end select

   DF2DX5Y = s

  end function


 ! *****************************************************************************
  real function DF2DX4Y2(xd1,yd1,NDERIVATIVE)
   IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
     case(23)
      s = 48.d0
   end select

   DF2DX4Y2 = s

  end function


 ! *****************************************************************************
  real function DF2DX3Y3(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
    case(24)
      s = 36.d0
   end select

   DF2DX3Y3 = s

  end function



 ! *****************************************************************************
  real function DF2DX2Y4(xd1,yd1,NDERIVATIVE)
    IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
    case(25)
      s = 48.d0

   end select

   DF2DX2Y4 = s

  end function


 ! *****************************************************************************
  real function DF2DXY5(xd1,yd1,NDERIVATIVE)
    IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
     case(26)
      s = 120.d0
   end select

  DF2DXY5 = s

  end function


 ! *****************************************************************************
  real function DF2DY6(xd1,yd1,NDERIVATIVE)
  IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
real::s

   S=0.0D0
   select case(NDERIVATIVE)
     case(27)
      s = 720.d0
   end select

   DF2DY6 = s

  end function	


  
  
  REAL FUNCTION TL2DY(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=1/hxc
 CASE(5)  ; S=yd1/hxc**2
 CASE(9)  ; S=yd1**2/(2.0D0*hxc**3)
 CASE(14)  ; S=yd1**3/(6.0D0*hxc**4)
 CASE(20)  ; S=yd1**4/(24.0D0*hxc**5)
 CASE(27)  ; S=yd1**5/(120.0D0*hxc**6)
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=xd1/hxc**2
 CASE(8)  ; S=(xd1*yd1)/hxc**3
 CASE(13)  ; S=(xd1*yd1**2)/(2.0D0*hxc**4)
 CASE(19)  ; S=(xd1*yd1**3)/(6.0D0*hxc**5)
 CASE(26)  ; S=(xd1*yd1**4)/(24.0D0*hxc**6)
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=xd1**2/(2.0D0*hxc**3)
 CASE(12)  ; S=(xd1**2*yd1)/(2.0D0*hxc**4)
 CASE(18)  ; S=(xd1**2*yd1**2)/(4.0D0*hxc**5)
 CASE(25)  ; S=(xd1**2*yd1**3)/(12.0D0*hxc**6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=xd1**3/(6.0D0*hxc**4)
 CASE(17)  ; S=(xd1**3*yd1)/(6.0D0*hxc**5)
 CASE(24)  ; S=(xd1**3*yd1**2)/(12.0D0*hxc**6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=xd1**4/(24.0D0*hxc**5)
 CASE(23)  ; S=(xd1**4*yd1)/(24.0D0*hxc**6)
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=xd1**5/(120.0D0*hxc**6)
 CASE(21)  ; S=0.0D0
 
 
 
 
  End select

    TL2DY = S

  end function

REAL FUNCTION  TL2DY2(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=hxc**(-2)
 CASE(9)  ; S=yd1/hxc**3
 CASE(14)  ; S=yd1**2/(2.0D0*hxc**4)
 CASE(20)  ; S=yd1**3/(6.0D0*hxc**5)
 CASE(27)  ; S=yd1**4/(24.0D0*hxc**6)
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=xd1/hxc**3
 CASE(13)  ; S=(xd1*yd1)/hxc**4
 CASE(19)  ; S=(xd1*yd1**2)/(2.0D0*hxc**5)
 CASE(26)  ; S=(xd1*yd1**3)/(6.0D0*hxc**6)
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=xd1**2/(2.0D0*hxc**4)
 CASE(18)  ; S=(xd1**2*yd1)/(2.0D0*hxc**5)
 CASE(25)  ; S=(xd1**2*yd1**2)/(4.0D0*hxc**6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=xd1**3/(6.0D0*hxc**5)
 CASE(24)  ; S=(xd1**3*yd1)/(6.0D0*hxc**6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=xd1**4/(24.0D0*hxc**6)
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
 
  End select

   TL2DY2 = s

  end function

REAL FUNCTION  TL2DY3(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=hxc**(-3)
 CASE(14)  ; S=yd1/hxc**4
 CASE(20)  ; S=yd1**2/(2.0D0*hxc**5)
 CASE(27)  ; S=yd1**3/(6.0D0*hxc**6)
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=xd1/hxc**4
 CASE(19)  ; S=(xd1*yd1)/hxc**5
 CASE(26)  ; S=(xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=xd1**2/(2.0D0*hxc**5)
 CASE(25)  ; S=(xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=xd1**3/(6.0D0*hxc**6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DY3 = s

  end function

REAL FUNCTION  TL2DY4(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=hxc**(-4)
 CASE(20)  ; S=yd1/hxc**5
 CASE(27)  ; S=yd1**2/(2.0D0*hxc**6)
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=xd1/hxc**5
 CASE(26)  ; S=(xd1*yd1)/hxc**6
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=xd1**2/(2.0D0*hxc**6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
  End select

    TL2DY4 = s

  end function

REAL FUNCTION  TL2DY5(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=hxc**(-5)
 CASE(27)  ; S=yd1/hxc**6
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=xd1/hxc**6
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DY5 = s

  end function

REAL FUNCTION  TL2DY6(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=hxc**(-6)
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DY6 = s

  end function

REAL FUNCTION  TL2DX(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=1/hxc
 CASE(4)  ; S=yd1/hxc**2
 CASE(8)  ; S=yd1**2/(2.0D0*hxc**3)
 CASE(13)  ; S=yd1**3/(6.0D0*hxc**4)
 CASE(19)  ; S=yd1**4/(24.0D0*hxc**5)
 CASE(26)  ; S=yd1**5/(120.0D0*hxc**6)
 CASE(3)  ; S=xd1/hxc**2
 CASE(7)  ; S=(xd1*yd1)/hxc**3
 CASE(12)  ; S=(xd1*yd1**2)/(2.0D0*hxc**4)
 CASE(18)  ; S=(xd1*yd1**3)/(6.0D0*hxc**5)
 CASE(25)  ; S=(xd1*yd1**4)/(24.0D0*hxc**6)
 CASE(6)  ; S=xd1**2/(2.0D0*hxc**3)
 CASE(11)  ; S=(xd1**2*yd1)/(2.0D0*hxc**4)
 CASE(17)  ; S=(xd1**2*yd1**2)/(4.0D0*hxc**5)
 CASE(24)  ; S=(xd1**2*yd1**3)/(12.0D0*hxc**6)
 CASE(10)  ; S=xd1**3/(6.0D0*hxc**4)
 CASE(16)  ; S=(xd1**3*yd1)/(6.0D0*hxc**5)
 CASE(23)  ; S=(xd1**3*yd1**2)/(12.0D0*hxc**6)
 CASE(15)  ; S=xd1**4/(24.0D0*hxc**5)
 CASE(22)  ; S=(xd1**4*yd1)/(24.0D0*hxc**6)
 CASE(21)  ; S=xd1**5/(120.0D0*hxc**6)
 
 
 
  End select

   TL2DX = s

  end function

REAL FUNCTION  TL2DXY(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=hxc**(-2)
 CASE(8)  ; S=yd1/hxc**3
 CASE(13)  ; S=yd1**2/(2.0D0*hxc**4)
 CASE(19)  ; S=yd1**3/(6.0D0*hxc**5)
 CASE(26)  ; S=yd1**4/(24.0D0*hxc**6)
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=xd1/hxc**3
 CASE(12)  ; S=(xd1*yd1)/hxc**4
 CASE(18)  ; S=(xd1*yd1**2)/(2.0D0*hxc**5)
 CASE(25)  ; S=(xd1*yd1**3)/(6.0D0*hxc**6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=xd1**2/(2.0D0*hxc**4)
 CASE(17)  ; S=(xd1**2*yd1)/(2.0D0*hxc**5)
 CASE(24)  ; S=(xd1**2*yd1**2)/(4.0D0*hxc**6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=xd1**3/(6.0D0*hxc**5)
 CASE(23)  ; S=(xd1**3*yd1)/(6.0D0*hxc**6)
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=xd1**4/(24.0D0*hxc**6)
 CASE(21)  ; S=0.0D0
 
 
  End select

    TL2DXY = s

  end function

REAL FUNCTION  TL2DXY2(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=hxc**(-3)
 CASE(13)  ; S=yd1/hxc**4
 CASE(19)  ; S=yd1**2/(2.0D0*hxc**5)
 CASE(26)  ; S=yd1**3/(6.0D0*hxc**6)
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=xd1/hxc**4
 CASE(18)  ; S=(xd1*yd1)/hxc**5
 CASE(25)  ; S=(xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=xd1**2/(2.0D0*hxc**5)
 CASE(24)  ; S=(xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=xd1**3/(6.0D0*hxc**6)
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DXY2 = s

  end function

REAL FUNCTION  TL2DXY3(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=hxc**(-4)
 CASE(19)  ; S=yd1/hxc**5
 CASE(26)  ; S=yd1**2/(2.0D0*hxc**6)
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=xd1/hxc**5
 CASE(25)  ; S=(xd1*yd1)/hxc**6
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=xd1**2/(2.0D0*hxc**6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
  End select

    TL2DXY3 = s

  end function

REAL FUNCTION  TL2DXY4(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=hxc**(-5)
 CASE(26)  ; S=yd1/hxc**6
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=xd1/hxc**6
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
  End select

    TL2DXY4 = s

  end function

REAL FUNCTION  TL2DXY5(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=hxc**(-6)
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
  End select

   TL2DXY5 = s

  end function

REAL FUNCTION  TL2DX2(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=hxc**(-2)
 CASE(7)  ; S=yd1/hxc**3
 CASE(12)  ; S=yd1**2/(2.0D0*hxc**4)
 CASE(18)  ; S=yd1**3/(6.0D0*hxc**5)
 CASE(25)  ; S=yd1**4/(24.0D0*hxc**6)
 CASE(6)  ; S=xd1/hxc**3
 CASE(11)  ; S=(xd1*yd1)/hxc**4
 CASE(17)  ; S=(xd1*yd1**2)/(2.0D0*hxc**5)
 CASE(24)  ; S=(xd1*yd1**3)/(6.0D0*hxc**6)
 CASE(10)  ; S=xd1**2/(2.0D0*hxc**4)
 CASE(16)  ; S=(xd1**2*yd1)/(2.0D0*hxc**5)
 CASE(23)  ; S=(xd1**2*yd1**2)/(4.0D0*hxc**6)
 CASE(15)  ; S=xd1**3/(6.0D0*hxc**5)
 CASE(22)  ; S=(xd1**3*yd1)/(6.0D0*hxc**6)
 CASE(21)  ; S=xd1**4/(24.0D0*hxc**6)
 
 
 
  End select

    TL2DX2 = s

  end function

REAL FUNCTION  TL2DX2Y(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=hxc**(-3)
 CASE(12)  ; S=yd1/hxc**4
 CASE(18)  ; S=yd1**2/(2.0D0*hxc**5)
 CASE(25)  ; S=yd1**3/(6.0D0*hxc**6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=xd1/hxc**4
 CASE(17)  ; S=(xd1*yd1)/hxc**5
 CASE(24)  ; S=(xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=xd1**2/(2.0D0*hxc**5)
 CASE(23)  ; S=(xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=xd1**3/(6.0D0*hxc**6)
 CASE(21)  ; S=0.0D0
 
 
  End select

    TL2DX2Y = s

  end function

REAL FUNCTION  TL2DX2Y2(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=hxc**(-4)
 CASE(18)  ; S=yd1/hxc**5
 CASE(25)  ; S=yd1**2/(2.0D0*hxc**6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=xd1/hxc**5
 CASE(24)  ; S=(xd1*yd1)/hxc**6
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=xd1**2/(2.0D0*hxc**6)
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DX2Y2= s

  end function

REAL FUNCTION  TL2DX2Y3(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=hxc**(-5)
 CASE(25)  ; S=yd1/hxc**6
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=xd1/hxc**6
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
  End select

    TL2DX2Y3 = s

  end function

REAL FUNCTION  TL2DX2Y4(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=hxc**(-6)
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
  End select

    TL2DX2Y4 = s

  end function

REAL FUNCTION  TL2DX3(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=hxc**(-3)
 CASE(11)  ; S=yd1/hxc**4
 CASE(17)  ; S=yd1**2/(2.0D0*hxc**5)
 CASE(24)  ; S=yd1**3/(6.0D0*hxc**6)
 CASE(10)  ; S=xd1/hxc**4
 CASE(16)  ; S=(xd1*yd1)/hxc**5
 CASE(23)  ; S=(xd1*yd1**2)/(2.0D0*hxc**6)
 CASE(15)  ; S=xd1**2/(2.0D0*hxc**5)
 CASE(22)  ; S=(xd1**2*yd1)/(2.0D0*hxc**6)
 CASE(21)  ; S=xd1**3/(6.0D0*hxc**6)
 
 
  End select

    TL2DX3 = s

  end function

REAL FUNCTION  TL2DX3Y(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=hxc**(-4)
 CASE(17)  ; S=yd1/hxc**5
 CASE(24)  ; S=yd1**2/(2.0D0*hxc**6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=xd1/hxc**5
 CASE(23)  ; S=(xd1*yd1)/hxc**6
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=xd1**2/(2.0D0*hxc**6)
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DX3Y = s

  end function

REAL FUNCTION  TL2DX3Y2(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=hxc**(-5)
 CASE(24)  ; S=yd1/hxc**6
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=xd1/hxc**6
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DX3Y2 = s

  end function

REAL FUNCTION  TL2DX3Y3(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=hxc**(-6)
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
  End select

    TL2DX3Y3 = s

  end function

REAL FUNCTION  TL2DX4(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=hxc**(-4)
 CASE(16)  ; S=yd1/hxc**5
 CASE(23)  ; S=yd1**2/(2.0D0*hxc**6)
 CASE(15)  ; S=xd1/hxc**5
 CASE(22)  ; S=(xd1*yd1)/hxc**6
 CASE(21)  ; S=xd1**2/(2.0D0*hxc**6)
 
 
 
  End select

    TL2DX4 = s

  end function

REAL FUNCTION  TL2DX4Y(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=hxc**(-5)
 CASE(23)  ; S=yd1/hxc**6
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=xd1/hxc**6
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DX4Y= s

  end function

REAL FUNCTION  TL2DX4Y2(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=hxc**(-6)
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DX4Y2 = s

  end function

REAL FUNCTION  TL2DX5(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=hxc**(-5)
 CASE(22)  ; S=yd1/hxc**6
 CASE(21)  ; S=xd1/hxc**6
 
 
  End select

    TL2DX5 = s

  end function

REAL FUNCTION  TL2DX5Y(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=hxc**(-6)
 CASE(21)  ; S=0.0D0
 
 
 
  End select

    TL2DX5Y = s

  end function

REAL FUNCTION  TL2DX6(XD1,YD1,NDERIVATIVE)
IMPLICIT NONE
INTEGER,INTENT(IN)::NDERIVATIVE
REAL,INTENT(IN)::XD1,YD1
REAL::S,HXC   
hxc=sqrt(IELEM(N,Iconsidered)%totvolume)




   S=0.0D0
   select case(NDERIVATIVE)
 CASE(2)  ; S=0.0D0
 CASE(5)  ; S=0.0D0
 CASE(9)  ; S=0.0D0
 CASE(14)  ; S=0.0D0
 CASE(20)  ; S=0.0D0
 CASE(27)  ; S=0.0D0
 CASE(1)  ; S=0.0D0
 CASE(4)  ; S=0.0D0
 CASE(8)  ; S=0.0D0
 CASE(13)  ; S=0.0D0
 CASE(19)  ; S=0.0D0
 CASE(26)  ; S=0.0D0
 CASE(3)  ; S=0.0D0
 CASE(7)  ; S=0.0D0
 CASE(12)  ; S=0.0D0
 CASE(18)  ; S=0.0D0
 CASE(25)  ; S=0.0D0
 CASE(6)  ; S=0.0D0
 CASE(11)  ; S=0.0D0
 CASE(17)  ; S=0.0D0
 CASE(24)  ; S=0.0D0
 CASE(10)  ; S=0.0D0
 CASE(16)  ; S=0.0D0
 CASE(23)  ; S=0.0D0
 CASE(15)  ; S=0.0D0
 CASE(22)  ; S=0.0D0
 CASE(21)  ; S=hxc**(-6)

End select

    TL2DX6 = s

  end function
  
  
  
  
  
  
  

END MODULE DERIVATIVES
