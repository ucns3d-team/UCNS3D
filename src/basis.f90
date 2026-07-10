module basis
!> @brief
!> this module includes all the basis functions for the polynomials used in 2d and 3d
use declaration
use library
implicit none

 contains

real function basis_real_power(x,p)
implicit none
!> @brief
!> small integer power helper for scalar basis evaluation
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
real,intent(in)::x
integer,intent(in)::p
integer::i

basis_real_power=1.0d0
do i=1,p
   basis_real_power=basis_real_power*x
end do

end function basis_real_power


real function basis_rec_value(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,kbasis)
implicit none
!> @brief
!> scalar 3D generic monomial basis value for linear-scheme GPU isolation
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
integer,intent(in)::number,iconsidered,number_of_dog,icompwrt,kbasis
real,intent(in)::x1,y1,z1
integer::degree,zpow,ypow,xpow,k
real::raw,oov

basis_rec_value=zero
if ((kbasis.lt.1).or.(kbasis.gt.number_of_dog)) return

k=0
raw=zero
degree_loop: do degree=1,number
   do zpow=0,degree
      do ypow=0,degree-zpow
         xpow=degree-zpow-ypow
         k=k+1
         if (k.eq.kbasis)then
            raw=basis_real_power(x1,xpow)*basis_real_power(y1,ypow)*basis_real_power(z1,zpow)
            exit degree_loop
         end if
      end do
   end do
end do degree_loop

if (icompwrt.eq.-2)then
   oov=1.0d0/(ielem_totvolume(iconsidered))
   basis_rec_value=raw-(integ_basis_dg_value(kbasis,iconsidered))*oov
else if (icompwrt.eq.0)then
   oov=1.0d0/(rec_volume(1,1,iconsidered))
   basis_rec_value=raw-(integ_basis_value(kbasis,iconsidered))*oov
else if (icompwrt.eq.1)then
   oov=1.0d0/(rec_volume(1,1,iconsidered))
   basis_rec_value=raw-(integ_basis_valuec(kbasis,iconsidered))*oov
else
   basis_rec_value=raw
end if

end function basis_rec_value


real function basis_rec2d_value(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,kbasis)
implicit none
!> @brief
!> scalar 2D generic monomial basis value for linear-scheme GPU isolation
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
integer,intent(in)::number,iconsidered,number_of_dog,icompwrt,kbasis
real,intent(in)::x1,y1
integer::degree,ypow,xpow,k
real::raw,oov

basis_rec2d_value=zero
if ((kbasis.lt.1).or.(kbasis.gt.number_of_dog)) return

k=0
raw=zero
degree_loop_2d: do degree=1,number
   do ypow=0,degree
      xpow=degree-ypow
      k=k+1
      if (k.eq.kbasis)then
         raw=basis_real_power(x1,xpow)*basis_real_power(y1,ypow)
         exit degree_loop_2d
      end if
   end do
end do degree_loop_2d

if (icompwrt.eq.-2)then
   oov=1.0d0/(ielem_totvolume(iconsidered))
   basis_rec2d_value=raw-(integ_basis_dg_value(kbasis,iconsidered))*oov
else if (icompwrt.eq.0)then
   oov=1.0d0/(rec_volume(1,1,iconsidered))
   basis_rec2d_value=raw-(integ_basis_value(kbasis,iconsidered))*oov
else if (icompwrt.eq.1)then
   oov=1.0d0/(rec_volume(1,1,iconsidered))
   basis_rec2d_value=raw-(integ_basis_valuec(kbasis,iconsidered))*oov
else
   basis_rec2d_value=raw
end if

end function basis_rec2d_value


subroutine basis_rec_fill(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt,basis_out)
implicit none
!> @brief
!> fixed-size subroutine variant of basis_rec for accelerator kernels
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
integer,intent(in)::number,iconsidered,number_of_dog,icompwrt
real,intent(in)::x1,y1,z1
real,intent(out)::basis_out(1:gpu_max_dof)
integer::degree,zpow,ypow,xpow,k,ip,kstart
real::raw,oov,xterm,yterm,zterm

basis_out=zero
oov=1.0d0/(rec_volume(1,1,iconsidered))
if (icompwrt.eq.-2) oov=1.0d0/(ielem_totvolume(iconsidered))

k=0
do degree=1,number
   do zpow=0,degree
      zterm=1.0d0
      do ip=1,zpow
         zterm=zterm*z1
      end do
      do ypow=0,degree-zpow
         xpow=degree-zpow-ypow
         k=k+1
         if (k.le.number_of_dog)then
            xterm=1.0d0
            do ip=1,xpow
               xterm=xterm*x1
            end do
            yterm=1.0d0
            do ip=1,ypow
               yterm=yterm*y1
            end do
            raw=xterm*yterm*zterm
            select case (icompwrt)
            case(-2)
               basis_out(k)=raw-(integ_basis_dg_value(k,iconsidered))*oov
            case(0)
               basis_out(k)=raw-(integ_basis_value(k,iconsidered))*oov
            case(1)
               basis_out(k)=raw-(integ_basis_valuec(k,iconsidered))*oov
            case default
               basis_out(k)=raw
            end select
         end if
      end do
   end do
end do

kstart=k+1
do k=kstart,number_of_dog
   select case (icompwrt)
   case(-2)
      basis_out(k)=-(integ_basis_dg_value(k,iconsidered))*oov
   case(0)
      basis_out(k)=-(integ_basis_value(k,iconsidered))*oov
   case(1)
      basis_out(k)=-(integ_basis_valuec(k,iconsidered))*oov
   case default
      basis_out(k)=zero
   end select
end do

end subroutine basis_rec_fill


subroutine basis_rec2d_fill(n,x1,y1,number,iconsidered,number_of_dog,icompwrt,basis_out)
implicit none
!> @brief
!> fixed-size subroutine variant of basis_rec2d for accelerator kernels
#if defined(gpu) || defined(xpu)
!$omp declare target
#endif
integer,intent(in)::n
integer,intent(in)::number,iconsidered,number_of_dog,icompwrt
real,intent(in)::x1,y1
real,intent(out)::basis_out(1:gpu_max_dof)
integer::degree,ypow,xpow,k,ip,kstart
real::raw,oov,xterm,yterm

basis_out=zero
oov=1.0d0/(rec_volume(1,1,iconsidered))
if (icompwrt.eq.-2) oov=1.0d0/(ielem_totvolume(iconsidered))

k=0
do degree=1,number
   do ypow=0,degree
      xpow=degree-ypow
      k=k+1
      if (k.le.number_of_dog)then
         xterm=1.0d0
         do ip=1,xpow
            xterm=xterm*x1
         end do
         yterm=1.0d0
         do ip=1,ypow
            yterm=yterm*y1
         end do
         raw=xterm*yterm
         select case (icompwrt)
         case(-2)
            basis_out(k)=raw-(integ_basis_dg_value(k,iconsidered))*oov
         case(0)
            basis_out(k)=raw-(integ_basis_value(k,iconsidered))*oov
         case(1)
            basis_out(k)=raw-(integ_basis_valuec(k,iconsidered))*oov
         case default
            basis_out(k)=raw
         end select
      end if
   end do
end do

kstart=k+1
do k=kstart,number_of_dog
   select case (icompwrt)
   case(-2)
      basis_out(k)=-(integ_basis_dg_value(k,iconsidered))*oov
   case(0)
      basis_out(k)=-(integ_basis_value(k,iconsidered))*oov
   case(1)
      basis_out(k)=-(integ_basis_valuec(k,iconsidered))*oov
   case default
      basis_out(k)=zero
   end select
end do

end subroutine basis_rec2d_fill
 
 
function basis_rec(n,x1,y1,z1,number,iconsidered,number_of_dog,icompwrt)
implicit none
!> @brief
!> this function returns the value of the basis function for a specific polynomial order and coordinates
#ifdef gpu
!!$omp declare target disabled: array-return basis routine replaced by basis_rec_fill
#endif
integer,intent(in)::n
integer,intent(in)::number,iconsidered,number_of_dog,icompwrt
real,intent(in)::x1,y1,z1
real::oov,hxc
real,dimension(number_of_dog)::basis_rec
real,dimension(number_of_dog)::sb
sb=zero
oov=1.0d0/(rec_volume(1,1,iconsidered))


hxc=(sqrt(ielem_totvolume(iconsidered)))

 if (poly.eq.1)then
 select case(number)
   case(1)
    	!first order functions (2nd-order of accuracy 3)
    sb(1)= x1
    sb(2)= y1
    sb(3)= z1

case(2)
	! second order functions (3rd-order of accuracy 4-9)
    sb(1)= x1
    sb(2)= y1
    sb(3)= z1	
    sb(4)= x1**2
    sb(5)= x1*y1
    sb(6)= y1**2
    sb(7)= x1*z1
    sb(8)= y1*z1
    sb(9)= z1**2
    case(3)
	! third order functions (4th-order of accuracy  10-19)
    sb(1)= x1
    sb(2)= y1
    sb(3)= z1	
    sb(4)= x1**2
    sb(5)= x1*y1
    sb(6)= y1**2
    sb(7)= x1*z1
    sb(8)= y1*z1
    sb(9)= z1**2	
    sb(10)= (x1**3)
    sb(11)= (sb(4))*y1
    sb(12)= (x1)*(sb(6))
    sb(13)= (y1**3)
    sb(14)= (sb(4))*z1
    sb(15)= (x1)*(y1)*(z1)
    sb(16)= (sb(6))*z1
    sb(17)= (x1)*(sb(9))
    sb(18)= (y1)*(sb(9))
    sb(19)= (z1**3)
    case(4)
	! fourth order functions (5th-order of accuracy 20-34)
     sb(1)= x1
    sb(2)= y1
    sb(3)= z1	
    sb(4)= x1**2
    sb(5)= x1*y1
    sb(6)= y1**2
    sb(7)= x1*z1
    sb(8)= y1*z1
    sb(9)= z1**2	
    sb(10)= (x1**3)
    sb(11)= (sb(4))*y1
    sb(12)= (x1)*(sb(6))
    sb(13)= (y1**3)
    sb(14)= (sb(4))*z1
    sb(15)= (x1)*sb(8)
    sb(16)= (sb(6))*z1
    sb(17)= (x1)*(sb(9))
    sb(18)= (y1)*(sb(9))
    sb(19)= (z1**3)	
    sb(20)= (x1**4)
    sb(21)= (sb(10))*y1
    sb(22)= (sb(4))*(sb(6))
    sb(23)= (x1)*(sb(13))
    sb(24)= (y1**4)
    sb(25)= (sb(10))*(z1)
    sb(26)= (sb(4))*sb(8)
    sb(27)= sb(7)*(sb(6))
    sb(28)= (sb(13))*(z1)
    sb(29)= (sb(4))*(sb(9))
    sb(30)= sb(5)*(sb(9))
    sb(31)= (sb(6))*(sb(9))
    sb(32)= (x1)*(sb(19))
    sb(33)= (y1)*(sb(19))
    sb(34)= (z1**4)
    case(5)
	! fifth order functions (6th-order of accuracy 35-55)
      sb(1)= x1
    sb(2)= y1
    sb(3)= z1	
    sb(4)= x1**2
    sb(5)= x1*y1
    sb(6)= y1**2
    sb(7)= x1*z1
    sb(8)= y1*z1
    sb(9)= z1**2	
    sb(10)= (x1**3)
    sb(11)= (x1**2)*y1
    sb(12)= (x1)*(y1**2)
    sb(13)= (y1**3)
    sb(14)= (x1**2)*z1
    sb(15)= (x1)*(y1)*(z1)
    sb(16)= (y1**2)*z1
    sb(17)= (x1)*(z1**2)
    sb(18)= (y1)*(z1**2)
    sb(19)= (z1**3)	
    sb(20)= (x1**4)
    sb(21)= (x1**3)*y1
    sb(22)= (x1**2)*(y1**2)
    sb(23)= (x1)*(y1**3)
    sb(24)= (y1**4)
    sb(25)= (x1**3)*(z1)
    sb(26)= (x1**2)*(y1)*(z1)
    sb(27)= (x1)*(y1**2)*(z1)
    sb(28)= (y1**3)*(z1)
    sb(29)= (x1**2)*(z1**2)
    sb(30)= (x1)*(y1)*(z1**2)
    sb(31)= (y1**2)*(z1**2)
    sb(32)= (x1)*(z1**3)
    sb(33)= (y1)*(z1**3)
    sb(34)= (z1**4)	
    sb(35)= (x1**5)
    sb(36)= (x1**4)*(y1)
    sb(37)= (x1**3)*(y1**2)
    sb(38)= (x1**2)*(y1**3)
    sb(39)= (x1)*(y1**4)
    sb(40)= (y1**5)
    sb(41)= (x1**4)*(z1)
    sb(42)= (x1**3)*(y1)*(z1)
    sb(43)= (x1**2)*(y1**2)*(z1)
    sb(44)= (x1)*(y1**3)*(z1)
    sb(45)= (y1**4)*(z1)
    sb(46)= (x1**3)*(z1**2)
    sb(47)= (x1**2)*(y1)*(z1**2)
    sb(48)= (x1)*(y1**2)*(z1**2)
    sb(49)= (y1**3)*(z1**2)
    sb(50)= (x1**2)*(z1**3)
    sb(51)= (x1)*(y1)*(z1**3)
    sb(52)= (y1**2)*(z1**3)
    sb(53)= (x1)*(z1**4)
    sb(54)= (y1)*(z1**4)
    sb(55)= (z1**5)
    
   ! sixth order functions (7th-order of accuracy 56-83)
     case(6)
      sb(1)= x1
    sb(2)= y1
    sb(3)= z1	
    sb(4)= x1**2
    sb(5)= x1*y1
    sb(6)= y1**2
    sb(7)= x1*z1
    sb(8)= y1*z1
    sb(9)= z1**2	
    sb(10)= (x1**3)
    sb(11)= (x1**2)*y1
    sb(12)= (x1)*(y1**2)
    sb(13)= (y1**3)
    sb(14)= (x1**2)*z1
    sb(15)= (x1)*(y1)*(z1)
    sb(16)= (y1**2)*z1
    sb(17)= (x1)*(z1**2)
    sb(18)= (y1)*(z1**2)
    sb(19)= (z1**3)	
    sb(20)= (x1**4)
    sb(21)= (x1**3)*y1
    sb(22)= (x1**2)*(y1**2)
    sb(23)= (x1)*(y1**3)
    sb(24)= (y1**4)
    sb(25)= (x1**3)*(z1)
    sb(26)= (x1**2)*(y1)*(z1)
    sb(27)= (x1)*(y1**2)*(z1)
    sb(28)= (y1**3)*(z1)
    sb(29)= (x1**2)*(z1**2)
    sb(30)= (x1)*(y1)*(z1**2)
    sb(31)= (y1**2)*(z1**2)
    sb(32)= (x1)*(z1**3)
    sb(33)= (y1)*(z1**3)
    sb(34)= (z1**4)	
    sb(35)= (x1**5)
    sb(36)= (x1**4)*(y1)
    sb(37)= (x1**3)*(y1**2)
    sb(38)= (x1**2)*(y1**3)
    sb(39)= (x1)*(y1**4)
    sb(40)= (y1**5)
    sb(41)= (x1**4)*(z1)
    sb(42)= (x1**3)*(y1)*(z1)
    sb(43)= (x1**2)*(y1**2)*(z1)
    sb(44)= (x1)*(y1**3)*(z1)
    sb(45)= (y1**4)*(z1)
    sb(46)= (x1**3)*(z1**2)
    sb(47)= (x1**2)*(y1)*(z1**2)
    sb(48)= (x1)*(y1**2)*(z1**2)
    sb(49)= (y1**3)*(z1**2)
    sb(50)= (x1**2)*(z1**3)
    sb(51)= (x1)*(y1)*(z1**3)
    sb(52)= (y1**2)*(z1**3)
    sb(53)= (x1)*(z1**4)
    sb(54)= (y1)*(z1**4)
    sb(55)= (z1**5)
    sb(56)= (x1**6)
    sb(57)= (x1**5)*(y1)
    sb(58)= (x1**4)*(y1**2)
    sb(59)= (x1**3)*(y1**3)
    sb(60)= (x1**2)*(y1**4)
    sb(61)= (x1)*(y1**5)
    sb(62)= (y1**6)
    sb(63)= (x1**5)*(z1)
    sb(64)= (x1**4)*(y1)*(z1)
    sb(65)= (x1**3)*(y1**2)*(z1)
    sb(66)= (x1**2)*(y1**3)*(z1)
    sb(67)=(x1)*(y1**4)*(z1)
    sb(68)= (y1**5)*(z1)
    sb(69)= (x1**4)*(z1**2)
    sb(70)= (x1**3)*(y1)*(z1**2)
    sb(71)= (x1**2)*(y1**2)*(z1**2)
    sb(72)= (x1)*(y1**3)*(z1**2)
    sb(73)= (y1**4)*(z1**2)
    sb(74)= (x1**3)*(z1**3)
    sb(75)= (x1**2)*(y1)*(z1**3)
    sb(76)= (x1)*(y1**2)*(z1**3)
    sb(77)= (y1**3)*(z1**3)
    sb(78)= (x1**2)*(z1**4)
    sb(79)= (x1)*(y1)*(z1**4)
    sb(80)= (y1**2)*(z1**4)
    sb(81)= (x1)*(z1**5)
    sb(82)= (y1)*(z1**5)
    sb(83)= (z1**6)
    
    case(7)
    ! seventh order functions (8th-order of accuracy 84-119)
      sb(1)= x1
    sb(2)= y1
    sb(3)= z1	
    sb(4)= x1**2
    sb(5)= x1*y1
    sb(6)= y1**2
    sb(7)= x1*z1
    sb(8)= y1*z1
    sb(9)= z1**2	
    sb(10)= (x1**3)
    sb(11)= (x1**2)*y1
    sb(12)= (x1)*(y1**2)
    sb(13)= (y1**3)
    sb(14)= (x1**2)*z1
    sb(15)= (x1)*(y1)*(z1)
    sb(16)= (y1**2)*z1
    sb(17)= (x1)*(z1**2)
    sb(18)= (y1)*(z1**2)
    sb(19)= (z1**3)	
    sb(20)= (x1**4)
    sb(21)= (x1**3)*y1
    sb(22)= (x1**2)*(y1**2)
    sb(23)= (x1)*(y1**3)
    sb(24)= (y1**4)
    sb(25)= (x1**3)*(z1)
    sb(26)= (x1**2)*(y1)*(z1)
    sb(27)= (x1)*(y1**2)*(z1)
    sb(28)= (y1**3)*(z1)
    sb(29)= (x1**2)*(z1**2)
    sb(30)= (x1)*(y1)*(z1**2)
    sb(31)= (y1**2)*(z1**2)
    sb(32)= (x1)*(z1**3)
    sb(33)= (y1)*(z1**3)
    sb(34)= (z1**4)	
    sb(35)= (x1**5)
    sb(36)= (x1**4)*(y1)
    sb(37)= (x1**3)*(y1**2)
    sb(38)= (x1**2)*(y1**3)
    sb(39)= (x1)*(y1**4)
    sb(40)= (y1**5)
    sb(41)= (x1**4)*(z1)
    sb(42)= (x1**3)*(y1)*(z1)
    sb(43)= (x1**2)*(y1**2)*(z1)
    sb(44)= (x1)*(y1**3)*(z1)
    sb(45)= (y1**4)*(z1)
    sb(46)= (x1**3)*(z1**2)
    sb(47)= (x1**2)*(y1)*(z1**2)
    sb(48)= (x1)*(y1**2)*(z1**2)
    sb(49)= (y1**3)*(z1**2)
    sb(50)= (x1**2)*(z1**3)
    sb(51)= (x1)*(y1)*(z1**3)
    sb(52)= (y1**2)*(z1**3)
    sb(53)= (x1)*(z1**4)
    sb(54)= (y1)*(z1**4)
    sb(55)= (z1**5)
    sb(56)= (x1**6)
    sb(57)= (x1**5)*(y1)
    sb(58)= (x1**4)*(y1**2)
    sb(59)= (x1**3)*(y1**3)
    sb(60)= (x1**2)*(y1**4)
    sb(61)= (x1)*(y1**5)
    sb(62)= (y1**6)
    sb(63)= (x1**5)*(z1)
    sb(64)= (x1**4)*(y1)*(z1)
    sb(65)= (x1**3)*(y1**2)*(z1)
    sb(66)= (x1**2)*(y1**3)*(z1)
    sb(67)=(x1)*(y1**4)*(z1)
    sb(68)= (y1**5)*(z1)
    sb(69)= (x1**4)*(z1**2)
    sb(70)= (x1**3)*(y1)*(z1**2)
    sb(71)= (x1**2)*(y1**2)*(z1**2)
    sb(72)= (x1)*(y1**3)*(z1**2)
    sb(73)= (y1**4)*(z1**2)
    sb(74)= (x1**3)*(z1**3)
    sb(75)= (x1**2)*(y1)*(z1**3)
    sb(76)= (x1)*(y1**2)*(z1**3)
    sb(77)= (y1**3)*(z1**3)
    sb(78)= (x1**2)*(z1**4)
    sb(79)= (x1)*(y1)*(z1**4)
    sb(80)= (y1**2)*(z1**4)
    sb(81)= (x1)*(z1**5)
    sb(82)= (y1)*(z1**5)
    sb(83)= (z1**6)
    sb(84)= (x1**7)
    sb(85)= (x1**6)*(y1)
    sb(86)= (x1**5)*(y1**2)
    sb(87)= (x1**4)*(y1**3)
    sb(88)= (x1**3)*(y1**4)
    sb(89)= (x1**2)*(y1**5)
    sb(90)= (x1)*(y1**6)
    sb(91)= (y1**7)
    sb(92)= (x1**6)*(z1)
    sb(93)= (x1**5)*(y1)*(z1)
    sb(94)= (x1**4)*(y1**2)*(z1)
    sb(95)= (x1**3)*(y1**3)*(z1)
    sb(96)= (x1**2)*(y1**4)*(z1)
    sb(97)= (x1)*(y1**5)*(z1)
    sb(98)= (y1**6)*(z1)
    sb(99)= (x1**5)*(z1**2)
    sb(100)= (x1**4)*(y1)*(z1**2)
    sb(101)= (x1**3)*(y1**2)*(z1**2)
    sb(102)= (x1**2)*(y1**3)*(z1**2)
    sb(103)= (x1)*(y1**4)*(z1**2)
    sb(104)= (y1**5)*(z1**2)
    sb(105)= (x1**4)*(z1**3)
    sb(106)= (x1**3)*(y1)*(z1**3)
    sb(107)= (x1**2)*(y1**2)*(z1**3)
    sb(108)= (x1)*(y1**3)*(z1**3)
    sb(109)= (y1**4)*(z1**3)
    sb(110)= (x1**3)*(z1**4)
    sb(111)= (x1**2)*(y1)*(z1**4)
    sb(112)= (x1)*(y1**2)*(z1**4)
    sb(113)= (y1**3)*(z1**4)
    sb(114)= (x1**2)*(z1**5)
    sb(115)= (x1)*(y1)*(z1**5)
    sb(116)= (y1**2)*(z1**5)
    sb(117)= (x1)*(z1**6)
    sb(118)= (y1)*(z1**6)
    sb(119)= (z1**7)
    case(8)
	!eighth order functions (9th-order of accuracy 85-164)
      sb(1)= x1
    sb(2)= y1
    sb(3)= z1	
    sb(4)= x1**2
    sb(5)= x1*y1
    sb(6)= y1**2
    sb(7)= x1*z1
    sb(8)= y1*z1
    sb(9)= z1**2	
    sb(10)= (x1**3)
    sb(11)= (x1**2)*y1
    sb(12)= (x1)*(y1**2)
    sb(13)= (y1**3)
    sb(14)= (x1**2)*z1
    sb(15)= (x1)*(y1)*(z1)
    sb(16)= (y1**2)*z1
    sb(17)= (x1)*(z1**2)
    sb(18)= (y1)*(z1**2)
    sb(19)= (z1**3)	
    sb(20)= (x1**4)
    sb(21)= (x1**3)*y1
    sb(22)= (x1**2)*(y1**2)
    sb(23)= (x1)*(y1**3)
    sb(24)= (y1**4)
    sb(25)= (x1**3)*(z1)
    sb(26)= (x1**2)*(y1)*(z1)
    sb(27)= (x1)*(y1**2)*(z1)
    sb(28)= (y1**3)*(z1)
    sb(29)= (x1**2)*(z1**2)
    sb(30)= (x1)*(y1)*(z1**2)
    sb(31)= (y1**2)*(z1**2)
    sb(32)= (x1)*(z1**3)
    sb(33)= (y1)*(z1**3)
    sb(34)= (z1**4)	
    sb(35)= (x1**5)
    sb(36)= (x1**4)*(y1)
    sb(37)= (x1**3)*(y1**2)
    sb(38)= (x1**2)*(y1**3)
    sb(39)= (x1)*(y1**4)
    sb(40)= (y1**5)
    sb(41)= (x1**4)*(z1)
    sb(42)= (x1**3)*(y1)*(z1)
    sb(43)= (x1**2)*(y1**2)*(z1)
    sb(44)= (x1)*(y1**3)*(z1)
    sb(45)= (y1**4)*(z1)
    sb(46)= (x1**3)*(z1**2)
    sb(47)= (x1**2)*(y1)*(z1**2)
    sb(48)= (x1)*(y1**2)*(z1**2)
    sb(49)= (y1**3)*(z1**2)
    sb(50)= (x1**2)*(z1**3)
    sb(51)= (x1)*(y1)*(z1**3)
    sb(52)= (y1**2)*(z1**3)
    sb(53)= (x1)*(z1**4)
    sb(54)= (y1)*(z1**4)
    sb(55)= (z1**5)
    sb(56)= (x1**6)
    sb(57)= (x1**5)*(y1)
    sb(58)= (x1**4)*(y1**2)
    sb(59)= (x1**3)*(y1**3)
    sb(60)= (x1**2)*(y1**4)
    sb(61)= (x1)*(y1**5)
    sb(62)= (y1**6)
    sb(63)= (x1**5)*(z1)
    sb(64)= (x1**4)*(y1)*(z1)
    sb(65)= (x1**3)*(y1**2)*(z1)
    sb(66)= (x1**2)*(y1**3)*(z1)
    sb(67)=(x1)*(y1**4)*(z1)
    sb(68)= (y1**5)*(z1)
    sb(69)= (x1**4)*(z1**2)
    sb(70)= (x1**3)*(y1)*(z1**2)
    sb(71)= (x1**2)*(y1**2)*(z1**2)
    sb(72)= (x1)*(y1**3)*(z1**2)
    sb(73)= (y1**4)*(z1**2)
    sb(74)= (x1**3)*(z1**3)
    sb(75)= (x1**2)*(y1)*(z1**3)
    sb(76)= (x1)*(y1**2)*(z1**3)
    sb(77)= (y1**3)*(z1**3)
    sb(78)= (x1**2)*(z1**4)
    sb(79)= (x1)*(y1)*(z1**4)
    sb(80)= (y1**2)*(z1**4)
    sb(81)= (x1)*(z1**5)
    sb(82)= (y1)*(z1**5)
    sb(83)= (z1**6)
    sb(84)= (x1**7)
    sb(85)= (x1**6)*(y1)
    sb(86)= (x1**5)*(y1**2)
    sb(87)= (x1**4)*(y1**3)
    sb(88)= (x1**3)*(y1**4)
    sb(89)= (x1**2)*(y1**5)
    sb(90)= (x1)*(y1**6)
    sb(91)= (y1**7)
    sb(92)= (x1**6)*(z1)
    sb(93)= (x1**5)*(y1)*(z1)
    sb(94)= (x1**4)*(y1**2)*(z1)
    sb(95)= (x1**3)*(y1**3)*(z1)
    sb(96)= (x1**2)*(y1**4)*(z1)
    sb(97)= (x1)*(y1**5)*(z1)
    sb(98)= (y1**6)*(z1)
    sb(99)= (x1**5)*(z1**2)
    sb(100)= (x1**4)*(y1)*(z1**2)
    sb(101)= (x1**3)*(y1**2)*(z1**2)
    sb(102)= (x1**2)*(y1**3)*(z1**2)
    sb(103)= (x1)*(y1**4)*(z1**2)
    sb(104)= (y1**5)*(z1**2)
    sb(105)= (x1**4)*(z1**3)
    sb(106)= (x1**3)*(y1)*(z1**3)
    sb(107)= (x1**2)*(y1**2)*(z1**3)
    sb(108)= (x1)*(y1**3)*(z1**3)
    sb(109)= (y1**4)*(z1**3)
    sb(110)= (x1**3)*(z1**4)
    sb(111)= (x1**2)*(y1)*(z1**4)
    sb(112)= (x1)*(y1**2)*(z1**4)
    sb(113)= (y1**3)*(z1**4)
    sb(114)= (x1**2)*(z1**5)
    sb(115)= (x1)*(y1)*(z1**5)
    sb(116)= (y1**2)*(z1**5)
    sb(117)= (x1)*(z1**6)
    sb(118)= (y1)*(z1**6)
    sb(119)= (z1**7)	
    sb(120)= (x1**8)
    sb(121)= (x1**7)*(y1)
    sb(122)= (x1**6)*(y1**2)
    sb(123)= (x1**5)*(y1**3)
    sb(124)= (x1**4)*(y1**4)
    sb(125)= (x1**3)*(y1**5)
    sb(126)= (x1**2)*(y1**6)
    sb(127)= (x1)*(y1**7)
    sb(128)= (y1**8)
    sb(129)= (x1**7)*(z1)
    sb(130)= (x1**6)*(y1)*(z1)
    sb(131)= (x1**5)*(y1**2)*(z1)
    sb(132)= (x1**4)*(y1**3)*(z1)
    sb(133)= (x1**3)*(y1**4)*(z1)
    sb(134)= (x1**2)*(y1**5)*(z1)
    sb(135)= (x1)*(y1**6)*(z1)
    sb(136)= (y1**7)*(z1)
    sb(137)= (x1**6)*(z1**2)
    sb(138)= (x1**5)*(y1)*(z1**2)
    sb(139)= (x1**4)*(y1**2)*(z1**2)
    sb(140)= (x1**3)*(y1**3)*(z1**2)
    sb(141)= (x1**2)*(y1**4)*(z1**2)
    sb(142)= (x1)*(y1**5)*(z1**2)
    sb(143)= (y1**6)*(z1**2)
    sb(144)= (x1**5)*(z1**3)
    sb(145)= (x1**4)*(y1)*(z1**3)
    sb(146)= (x1**3)*(y1**2)*(z1**3)
    sb(147)= (x1**2)*(y1**3)*(z1**3)
    sb(148)= (x1)*(y1**4)*(z1**3)
    sb(149)= (y1**5)*(z1**3)
    sb(150)= (x1**4)*(z1**4)
    sb(151)= (x1**3)*(y1)*(z1**4)
    sb(152)= (x1**2)*(y1**2)*(z1**4)
    sb(153)= (x1)*(y1**3)*(z1**4)
    sb(154)= (y1**4)*(z1**4)
    sb(155)= (x1**3)*(z1**5)
    sb(156)= (x1**2)*(y1)*(z1**5)
    sb(157)= (x1)*(y1**2)*(z1**5)
    sb(158)= (y1**3)*(z1**5)
    sb(159)= (x1**2)*(z1**6)
    sb(160)= (x1)*(y1)*(z1**6)
    sb(161)= (y1**2)*(z1**6)
    sb(162)= (x1)*(z1**7)
    sb(163)= (y1)*(z1**7)
    sb(164)= (z1**8)
    case(9)
    !ninth order functions (10th-order of accuracy 165-219)
      sb(1)= x1
    sb(2)= y1
    sb(3)= z1	
    sb(4)= x1**2
    sb(5)= x1*y1
    sb(6)= y1**2
    sb(7)= x1*z1
    sb(8)= y1*z1
    sb(9)= z1**2	
    sb(10)= (x1**3)
    sb(11)= (x1**2)*y1
    sb(12)= (x1)*(y1**2)
    sb(13)= (y1**3)
    sb(14)= (x1**2)*z1
    sb(15)= (x1)*(y1)*(z1)
    sb(16)= (y1**2)*z1
    sb(17)= (x1)*(z1**2)
    sb(18)= (y1)*(z1**2)
    sb(19)= (z1**3)	
    sb(20)= (x1**4)
    sb(21)= (x1**3)*y1
    sb(22)= (x1**2)*(y1**2)
    sb(23)= (x1)*(y1**3)
    sb(24)= (y1**4)
    sb(25)= (x1**3)*(z1)
    sb(26)= (x1**2)*(y1)*(z1)
    sb(27)= (x1)*(y1**2)*(z1)
    sb(28)= (y1**3)*(z1)
    sb(29)= (x1**2)*(z1**2)
    sb(30)= (x1)*(y1)*(z1**2)
    sb(31)= (y1**2)*(z1**2)
    sb(32)= (x1)*(z1**3)
    sb(33)= (y1)*(z1**3)
    sb(34)= (z1**4)	
    sb(35)= (x1**5)
    sb(36)= (x1**4)*(y1)
    sb(37)= (x1**3)*(y1**2)
    sb(38)= (x1**2)*(y1**3)
    sb(39)= (x1)*(y1**4)
    sb(40)= (y1**5)
    sb(41)= (x1**4)*(z1)
    sb(42)= (x1**3)*(y1)*(z1)
    sb(43)= (x1**2)*(y1**2)*(z1)
    sb(44)= (x1)*(y1**3)*(z1)
    sb(45)= (y1**4)*(z1)
    sb(46)= (x1**3)*(z1**2)
    sb(47)= (x1**2)*(y1)*(z1**2)
    sb(48)= (x1)*(y1**2)*(z1**2)
    sb(49)= (y1**3)*(z1**2)
    sb(50)= (x1**2)*(z1**3)
    sb(51)= (x1)*(y1)*(z1**3)
    sb(52)= (y1**2)*(z1**3)
    sb(53)= (x1)*(z1**4)
    sb(54)= (y1)*(z1**4)
    sb(55)= (z1**5)
    sb(56)= (x1**6)
    sb(57)= (x1**5)*(y1)
    sb(58)= (x1**4)*(y1**2)
    sb(59)= (x1**3)*(y1**3)
    sb(60)= (x1**2)*(y1**4)
    sb(61)= (x1)*(y1**5)
    sb(62)= (y1**6)
    sb(63)= (x1**5)*(z1)
    sb(64)= (x1**4)*(y1)*(z1)
    sb(65)= (x1**3)*(y1**2)*(z1)
    sb(66)= (x1**2)*(y1**3)*(z1)
    sb(67)=(x1)*(y1**4)*(z1)
    sb(68)= (y1**5)*(z1)
    sb(69)= (x1**4)*(z1**2)
    sb(70)= (x1**3)*(y1)*(z1**2)
    sb(71)= (x1**2)*(y1**2)*(z1**2)
    sb(72)= (x1)*(y1**3)*(z1**2)
    sb(73)= (y1**4)*(z1**2)
    sb(74)= (x1**3)*(z1**3)
    sb(75)= (x1**2)*(y1)*(z1**3)
    sb(76)= (x1)*(y1**2)*(z1**3)
    sb(77)= (y1**3)*(z1**3)
    sb(78)= (x1**2)*(z1**4)
    sb(79)= (x1)*(y1)*(z1**4)
    sb(80)= (y1**2)*(z1**4)
    sb(81)= (x1)*(z1**5)
    sb(82)= (y1)*(z1**5)
    sb(83)= (z1**6)
    sb(84)= (x1**7)
    sb(85)= (x1**6)*(y1)
    sb(86)= (x1**5)*(y1**2)
    sb(87)= (x1**4)*(y1**3)
    sb(88)= (x1**3)*(y1**4)
    sb(89)= (x1**2)*(y1**5)
    sb(90)= (x1)*(y1**6)
    sb(91)= (y1**7)
    sb(92)= (x1**6)*(z1)
    sb(93)= (x1**5)*(y1)*(z1)
    sb(94)= (x1**4)*(y1**2)*(z1)
    sb(95)= (x1**3)*(y1**3)*(z1)
    sb(96)= (x1**2)*(y1**4)*(z1)
    sb(97)= (x1)*(y1**5)*(z1)
    sb(98)= (y1**6)*(z1)
    sb(99)= (x1**5)*(z1**2)
    sb(100)= (x1**4)*(y1)*(z1**2)
    sb(101)= (x1**3)*(y1**2)*(z1**2)
    sb(102)= (x1**2)*(y1**3)*(z1**2)
    sb(103)= (x1)*(y1**4)*(z1**2)
    sb(104)= (y1**5)*(z1**2)
    sb(105)= (x1**4)*(z1**3)
    sb(106)= (x1**3)*(y1)*(z1**3)
    sb(107)= (x1**2)*(y1**2)*(z1**3)
    sb(108)= (x1)*(y1**3)*(z1**3)
    sb(109)= (y1**4)*(z1**3)
    sb(110)= (x1**3)*(z1**4)
    sb(111)= (x1**2)*(y1)*(z1**4)
    sb(112)= (x1)*(y1**2)*(z1**4)
    sb(113)= (y1**3)*(z1**4)
    sb(114)= (x1**2)*(z1**5)
    sb(115)= (x1)*(y1)*(z1**5)
    sb(116)= (y1**2)*(z1**5)
    sb(117)= (x1)*(z1**6)
    sb(118)= (y1)*(z1**6)
    sb(119)= (z1**7)	
    sb(120)= (x1**8)
    sb(121)= (x1**7)*(y1)
    sb(122)= (x1**6)*(y1**2)
    sb(123)= (x1**5)*(y1**3)
    sb(124)= (x1**4)*(y1**4)
    sb(125)= (x1**3)*(y1**5)
    sb(126)= (x1**2)*(y1**6)
    sb(127)= (x1)*(y1**7)
    sb(128)= (y1**8)
    sb(129)= (x1**7)*(z1)
    sb(130)= (x1**6)*(y1)*(z1)
    sb(131)= (x1**5)*(y1**2)*(z1)
    sb(132)= (x1**4)*(y1**3)*(z1)
    sb(133)= (x1**3)*(y1**4)*(z1)
    sb(134)= (x1**2)*(y1**5)*(z1)
    sb(135)= (x1)*(y1**6)*(z1)
    sb(136)= (y1**7)*(z1)
    sb(137)= (x1**6)*(z1**2)
    sb(138)= (x1**5)*(y1)*(z1**2)
    sb(139)= (x1**4)*(y1**2)*(z1**2)
    sb(140)= (x1**3)*(y1**3)*(z1**2)
    sb(141)= (x1**2)*(y1**4)*(z1**2)
    sb(142)= (x1)*(y1**5)*(z1**2)
    sb(143)= (y1**6)*(z1**2)
    sb(144)= (x1**5)*(z1**3)
    sb(145)= (x1**4)*(y1)*(z1**3)
    sb(146)= (x1**3)*(y1**2)*(z1**3)
    sb(147)= (x1**2)*(y1**3)*(z1**3)
    sb(148)= (x1)*(y1**4)*(z1**3)
    sb(149)= (y1**5)*(z1**3)
    sb(150)= (x1**4)*(z1**4)
    sb(151)= (x1**3)*(y1)*(z1**4)
    sb(152)= (x1**2)*(y1**2)*(z1**4)
    sb(153)= (x1)*(y1**3)*(z1**4)
    sb(154)= (y1**4)*(z1**4)
    sb(155)= (x1**3)*(z1**5)
    sb(156)= (x1**2)*(y1)*(z1**5)
    sb(157)= (x1)*(y1**2)*(z1**5)
    sb(158)= (y1**3)*(z1**5)
    sb(159)= (x1**2)*(z1**6)
    sb(160)= (x1)*(y1)*(z1**6)
    sb(161)= (y1**2)*(z1**6)
    sb(162)= (x1)*(z1**7)
    sb(163)= (y1)*(z1**7)
    sb(164)= (z1**8)
    sb(165)= (x1**9)
    sb(166)= (x1**8)*(y1)
    sb(167)= (x1**7)*(y1**2)
    sb(168)= (x1**6)*(y1**3)
    sb(169)= (x1**5)*(y1**4)
    sb(170)= (x1**4)*(y1**5)
    sb(171)= (x1**3)*(y1**6)
    sb(172)= (x1**2)*(y1**7)
    sb(173)= (x1)*(y1**8)
    sb(174)= (y1**9)
    sb(175)= (x1**8)*(z1)
    sb(176)= (x1**7)*(y1)*(z1)
    sb(177)= (x1**6)*(y1**2)*(z1)
    sb(178)= (x1**5)*(y1**3)*(z1)
    sb(179)= (x1**4)*(y1**4)*(z1)
    sb(180)= (x1**3)*(y1**5)*(z1)
    sb(181)= (x1**2)*(y1**6)*(z1)
    sb(182)= (x1)*(y1**7)*(z1)
    sb(183)= (y1**8)*(z1)
    sb(184)= (x1**7)*(z1**2)
    sb(185)= (x1**6)*(y1)*(z1**2)
    sb(186)= (x1**5)*(y1**2)*(z1**2)
    sb(187)= (x1**4)*(y1**3)*(z1**2)
    sb(188)= (x1**3)*(y1**4)*(z1**2)
    sb(189)= (x1**2)*(y1**5)*(z1**2)
    sb(190)= (x1)*(y1**6)*(z1**2)
    sb(191)= (y1**7)*(z1**2)
    sb(192)= (x1**6)*(z1**3)
    sb(193)= (x1**5)*(y1)*(z1**3)
    sb(194)= (x1**4)*(y1**2)*(z1**3)
    sb(195)= (x1**3)*(y1**3)*(z1**3)
    sb(196)= (x1**2)*(y1**4)*(z1**3)
    sb(197)= (x1)*(y1**5)*(z1**3)
    sb(198)= (y1**6)*(z1**3)
    sb(199)= (x1**5)*(z1**4)
    sb(200)= (x1**4)*(y1)*(z1**4)
    sb(201)= (x1**3)*(y1**2)*(z1**4)
    sb(202)= (x1**2)*(y1**3)*(z1**4)
    sb(203)= (x1)*(y1**4)*(z1**4)
    sb(204)= (y1**5)*(z1**4)
    sb(205)= (x1**4)*(z1**5)
    sb(206)= (x1**3)*(y1)*(z1**5)
    sb(207)= (x1**2)*(y1**2)*(z1**5)
    sb(208)= (x1)*(y1**3)*(z1**5)
    sb(209)= (y1**4)*(z1**5)
    sb(210)= (x1**3)*(z1**6)
    sb(211)= (x1**2)*(y1)*(z1**6)
    sb(212)= (x1)*(y1**2)*(z1**6)
    sb(213)= (y1**3)*(z1**6)
    sb(214)= (x1**2)*(z1**7)
    sb(215)= (x1)*(y1)*(z1**7)
    sb(216)= (y1**2)*(z1**7)
    sb(217)= (x1)*(z1**8)
    sb(218)= (y1)*(z1**8)
    sb(219)= (z1**9)
    
    end select
   end if
   
   if (poly.eq.2)then
   
   
   
   select case(number)
   case(1)
    !first order functions (2nd-order of accuracy 3)
    sb(1)=-1.0d0 + 2.0d0*x1
    sb(2)=-1.0d0 + 2.0d0*y1
    sb(3)=-1.0d0 + 2.0d0*z1 
    case(2)
! second order functions (3rd-order of accuracy 4-9)
  sb(1)=-1.0d0 + 2.0d0*x1
    sb(2)=-1.0d0 + 2.0d0*y1
    sb(3)=-1.0d0 + 2.0d0*z1 
 sb(4)=1.0d0 - 6.0d0*x1 + 6.0d0*x1**2
 sb(5)=sb(1)*sb(2)
 sb(6)=sb(1)*sb(3)
 sb(7)=1.0d0 - 6.0d0*y1 + 6.0d0*y1**2 
 sb(8)=sb(2)*sb(3)
 sb(9)=1.0d0 - 6.0d0*z1 + 6.0d0*z1**2
 case(3)
! third order functions (4th-order of accuracy  10-19)
 sb(1)=-1.0d0 + 2.0d0*x1
 sb(2)=-1.0d0 + 2.0d0*y1
 sb(3)=-1.0d0 + 2.0d0*z1 
 sb(4)=1.0d0 - 6.0d0*x1 + 6.0d0*x1**2
 sb(5)=sb(1)*sb(2)
 sb(6)=sb(1)*sb(3)
 sb(7)=1.0d0 - 6.0d0*y1 + 6.0d0*y1**2 
 sb(8)=sb(2)*sb(3)
 sb(9)=1.0d0 - 6.0d0*z1 + 6.0d0*z1**2
 sb(10)=-1.0d0 + 12.0d0*x1 - 30.0d0*x1**2 + 20.0d0*x1**3
 sb(11)=sb(4)*sb(2)
 sb(12)=sb(1)*sb(7)
 sb(13)=-1.0d0 + 12.0d0*y1 - 30.0d0*y1**2 + 20.0d0*y1**3
 sb(14)=sb(7)*sb(3)
 sb(15)=sb(2)*sb(9)
 sb(16)=-1.0d0 + 12.0d0*z1 - 30.0d0*z1**2 + 20.0d0*z1**3
 sb(17)=sb(1)*sb(8)
 sb(18)=sb(4)*sb(3)
 sb(19)=sb(1)*sb(9)
 
 case(4)
! fourth order functions (5th-order of accuracy 20-34)
  sb(1)=-1.0d0 + 2.0d0*x1
 sb(2)=-1.0d0 + 2.0d0*y1
 sb(3)=-1.0d0 + 2.0d0*z1 
 sb(4)=1.0d0 - 6.0d0*x1 + 6.0d0*x1**2
 sb(5)=sb(1)*sb(2)
 sb(6)=sb(1)*sb(3)
 sb(7)=1.0d0 - 6.0d0*y1 + 6.0d0*y1**2 
 sb(8)=sb(2)*sb(3)
 sb(9)=1.0d0 - 6.0d0*z1 + 6.0d0*z1**2
 sb(10)=-1.0d0 + 12.0d0*x1 - 30.0d0*x1**2 + 20.0d0*x1**3
 sb(11)=sb(4)*sb(2)
 sb(12)=sb(1)*sb(7)
 sb(13)=-1.0d0 + 12.0d0*y1 - 30.0d0*y1**2 + 20.0d0*y1**3
 sb(14)=sb(7)*sb(3)
 sb(15)=sb(2)*sb(9)
 sb(16)=-1.0d0 + 12.0d0*z1 - 30.0d0*z1**2 + 20.0d0*z1**3
 sb(17)=sb(1)*sb(8)
 sb(18)=sb(4)*sb(3)
 sb(19)=sb(1)*sb(9)
 sb(20)=1.0d0 - 20.0d0*x1 + 90.0d0*x1**2 - 140.0d0*x1**3 + 70.0d0*x1**4
 sb(21)=sb(10)*sb(2)
 sb(22)=sb(10)*sb(3)
 sb(23)=sb(4)*sb(7)
 sb(24)=sb(11)*sb(3)
 sb(25)=sb(4)*sb(9)
 sb(26)=sb(1)*sb(13)
 sb(27)=sb(12)*sb(3)
 sb(28)=sb(1)*sb(15)
 sb(29)=sb(1)*sb(16)
 sb(30)=sb(13)*sb(3)
 sb(31)=sb(7)*sb(9)
 sb(32)=sb(2)*sb(16)
 sb(33)=1.0d0 - 20.0d0*z1 + 90.0d0*z1**2 - 140.0d0*z1**3 + 70.0d0*z1**4
 sb(34)=1.0d0 - 20.0d0*y1 + 90.0d0*y1**2 - 140.0d0*y1**3 + 70.0d0*y1**4
 case(5)
! fifth order functions (6th-order of accuracy 35-55)
  sb(1)=-1.0d0 + 2.0d0*x1
 sb(2)=-1.0d0 + 2.0d0*y1
 sb(3)=-1.0d0 + 2.0d0*z1 
 sb(4)=1.0d0 - 6.0d0*x1 + 6.0d0*x1**2
 sb(5)=sb(1)*sb(2)
 sb(6)=sb(1)*sb(3)
 sb(7)=1.0d0 - 6.0d0*y1 + 6.0d0*y1**2 
 sb(8)=sb(2)*sb(3)
 sb(9)=1.0d0 - 6.0d0*z1 + 6.0d0*z1**2
 sb(10)=-1.0d0 + 12.0d0*x1 - 30.0d0*x1**2 + 20.0d0*x1**3
 sb(11)=sb(4)*sb(2)
 sb(12)=sb(1)*sb(7)
 sb(13)=-1.0d0 + 12.0d0*y1 - 30.0d0*y1**2 + 20.0d0*y1**3
 sb(14)=sb(7)*sb(3)
 sb(15)=sb(2)*sb(9)
 sb(16)=-1.0d0 + 12.0d0*z1 - 30.0d0*z1**2 + 20.0d0*z1**3
 sb(17)=sb(1)*sb(8)
 sb(18)=sb(4)*sb(3)
 sb(19)=sb(1)*sb(9)
 sb(20)=1.0d0 - 20.0d0*x1 + 90.0d0*x1**2 - 140.0d0*x1**3 + 70.0d0*x1**4
 sb(21)=sb(10)*sb(2)
 sb(22)=sb(10)*sb(3)
 sb(23)=sb(4)*sb(7)
 sb(24)=sb(11)*sb(3)
 sb(25)=sb(4)*sb(9)
 sb(26)=sb(1)*sb(13)
 sb(27)=sb(12)*sb(3)
 sb(28)=sb(1)*sb(15)
 sb(29)=sb(1)*sb(16)
 sb(30)=sb(13)*sb(3)
 sb(31)=sb(7)*sb(9)
 sb(32)=sb(2)*sb(16)
 sb(33)=1.0d0 - 20.0d0*z1 + 90.0d0*z1**2 - 140.0d0*z1**3 + 70.0d0*z1**4
 sb(34)=1.0d0 - 20.0d0*y1 + 90.0d0*y1**2 - 140.0d0*y1**3 + 70.0d0*y1**4
 sb(35)=-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5
 sb(36)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)
 sb(37)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*z1)
 sb(38)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)
 sb(39)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*z1 + 6*z1**2)
 sb(40)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(41)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(42)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(43)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(44)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(45)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(46)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(47)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(48)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(49)=-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5
 sb(50)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(51)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(52)=(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(53)=(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(54)=-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5
 sb(55)=(-1 + 2*x1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 case(6)
   ! sixth order functions (7th-order of accuracy 56-83)
   sb(1)=-1 + 2*x1
    sb(2)=-1 + 2*y1
    sb(3)=-1 + 2*z1 
 sb(4)=1 - 6*x1 + 6*x1**2
 sb(5)=(-1 + 2*x1)*(-1 + 2*y1)
 sb(6)=(-1 + 2*x1)*(-1 + 2*z1)
 sb(7)=1 - 6*y1 + 6*y1**2 
 sb(8)=(-1 + 2*y1)*(-1 + 2*z1)
 sb(9)=1 - 6*z1 + 6*z1**2
 sb(10)=-1 + 12*x1 - 30*x1**2 + 20*x1**3
 sb(11)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)
 sb(12)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)
 sb(13)=-1 + 12*y1 - 30*y1**2 + 20*y1**3
 sb(14)=(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(15)=(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(16)=-1 + 12*z1 - 30*z1**2 + 20*z1**3
 sb(17)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(18)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*z1)
 sb(19)=(-1 + 2*x1)*(1 - 6*z1 + 6*z1**2)
 sb(20)=1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4
 sb(21)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)
 sb(22)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*z1)
 sb(23)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)
 sb(24)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(25)=(1 - 6*x1 + 6*x1**2)*(1 - 6*z1 + 6*z1**2)
 sb(26)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(27)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(28)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(29)=(-1 + 2*x1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(30)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(31)=(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(32)=(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(33)=1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4
 sb(34)=1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4
 sb(35)=-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5
 sb(36)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)
 sb(37)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*z1)
 sb(38)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)
 sb(39)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*z1 + 6*z1**2)
 sb(40)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(41)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(42)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(43)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(44)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(45)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(46)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(47)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(48)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(49)=-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5
 sb(50)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(51)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(52)=(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(53)=(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(54)=-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5
 sb(55)=(-1 + 2*x1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(56)=1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6
   sb(57)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)
   sb(58)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*z1)
   sb(59)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)
   sb(60)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(-1 + 2*z1)
   sb(61)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*z1 + 6*z1**2)
   sb(62)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
   sb(63)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
   sb(64)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
   sb(65)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(66)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
   sb(67)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
   sb(68)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
   sb(69)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(70)=(1 - 6*x1 + 6*x1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(71)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
   sb(72)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
   sb(73)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
   sb(74)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(75)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(76)=(-1 + 2*x1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
   sb(77)=1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6
   sb(78)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
   sb(79)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
   sb(80)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(81)=(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(82)=(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
   sb(83)=1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6
   case(7)
   ! seventh order functions (8th-order of accuracy 84-119)
     sb(1)=-1 + 2*x1
    sb(2)=-1 + 2*y1
    sb(3)=-1 + 2*z1 
 sb(4)=1 - 6*x1 + 6*x1**2
 sb(5)=(-1 + 2*x1)*(-1 + 2*y1)
 sb(6)=(-1 + 2*x1)*(-1 + 2*z1)
 sb(7)=1 - 6*y1 + 6*y1**2 
 sb(8)=(-1 + 2*y1)*(-1 + 2*z1)
 sb(9)=1 - 6*z1 + 6*z1**2
 sb(10)=-1 + 12*x1 - 30*x1**2 + 20*x1**3
 sb(11)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)
 sb(12)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)
 sb(13)=-1 + 12*y1 - 30*y1**2 + 20*y1**3
 sb(14)=(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(15)=(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(16)=-1 + 12*z1 - 30*z1**2 + 20*z1**3
 sb(17)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(18)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*z1)
 sb(19)=(-1 + 2*x1)*(1 - 6*z1 + 6*z1**2)
 sb(20)=1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4
 sb(21)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)
 sb(22)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*z1)
 sb(23)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)
 sb(24)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(25)=(1 - 6*x1 + 6*x1**2)*(1 - 6*z1 + 6*z1**2)
 sb(26)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(27)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(28)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(29)=(-1 + 2*x1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(30)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(31)=(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(32)=(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(33)=1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4
 sb(34)=1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4
 sb(35)=-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5
 sb(36)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)
 sb(37)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*z1)
 sb(38)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)
 sb(39)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*z1 + 6*z1**2)
 sb(40)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(41)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(42)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(43)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(44)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(45)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(46)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(47)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(48)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(49)=-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5
 sb(50)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(51)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(52)=(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(53)=(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(54)=-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5
 sb(55)=(-1 + 2*x1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(56)=1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6
   sb(57)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)
   sb(58)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*z1)
   sb(59)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)
   sb(60)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(-1 + 2*z1)
   sb(61)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*z1 + 6*z1**2)
   sb(62)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
   sb(63)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
   sb(64)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
   sb(65)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(66)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
   sb(67)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
   sb(68)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
   sb(69)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(70)=(1 - 6*x1 + 6*x1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(71)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
   sb(72)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
   sb(73)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
   sb(74)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(75)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(76)=(-1 + 2*x1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
   sb(77)=1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6
   sb(78)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
   sb(79)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
   sb(80)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(81)=(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(82)=(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
   sb(83)=1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6
 sb(84)=-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7
 sb(85)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*y1)
 sb(86)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*z1)
 sb(87)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*y1 + 6*y1**2)
 sb(88)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(89)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*z1 + 6*z1**2)
 sb(90)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(91)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(92)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(93)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(94)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(95)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(96)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(97)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(98)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(99)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
 sb(100)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(101)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(102)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(103)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(104)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(105)=(-1 + 2*x1)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)
 sb(106)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
 sb(107)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
 sb(108)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(109)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(110)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(111)=(-1 + 2*x1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(112)=-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7
 sb(113)=(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(-1 + 2*z1)
 sb(114)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(1 - 6*z1 + 6*z1**2)
 sb(115)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(116)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(117)=(1 - 6*y1 + 6*y1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(118)=(-1 + 2*y1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(119)=-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7 
 case(8)
 !eighth order functions (9th-order of accuracy 120-164)
  sb(1)=-1 + 2*x1
    sb(2)=-1 + 2*y1
    sb(3)=-1 + 2*z1 
 sb(4)=1 - 6*x1 + 6*x1**2
 sb(5)=(-1 + 2*x1)*(-1 + 2*y1)
 sb(6)=(-1 + 2*x1)*(-1 + 2*z1)
 sb(7)=1 - 6*y1 + 6*y1**2 
 sb(8)=(-1 + 2*y1)*(-1 + 2*z1)
 sb(9)=1 - 6*z1 + 6*z1**2
 sb(10)=-1 + 12*x1 - 30*x1**2 + 20*x1**3
 sb(11)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)
 sb(12)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)
 sb(13)=-1 + 12*y1 - 30*y1**2 + 20*y1**3
 sb(14)=(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(15)=(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(16)=-1 + 12*z1 - 30*z1**2 + 20*z1**3
 sb(17)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(18)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*z1)
 sb(19)=(-1 + 2*x1)*(1 - 6*z1 + 6*z1**2)
 sb(20)=1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4
 sb(21)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)
 sb(22)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*z1)
 sb(23)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)
 sb(24)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(25)=(1 - 6*x1 + 6*x1**2)*(1 - 6*z1 + 6*z1**2)
 sb(26)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(27)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(28)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(29)=(-1 + 2*x1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(30)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(31)=(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(32)=(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(33)=1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4
 sb(34)=1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4
 sb(35)=-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5
 sb(36)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)
 sb(37)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*z1)
 sb(38)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)
 sb(39)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*z1 + 6*z1**2)
 sb(40)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(41)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(42)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(43)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(44)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(45)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(46)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(47)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(48)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(49)=-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5
 sb(50)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(51)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(52)=(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(53)=(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(54)=-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5
 sb(55)=(-1 + 2*x1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(56)=1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6
   sb(57)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)
   sb(58)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*z1)
   sb(59)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)
   sb(60)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(-1 + 2*z1)
   sb(61)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*z1 + 6*z1**2)
   sb(62)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
   sb(63)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
   sb(64)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
   sb(65)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(66)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
   sb(67)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
   sb(68)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
   sb(69)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(70)=(1 - 6*x1 + 6*x1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(71)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
   sb(72)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
   sb(73)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
   sb(74)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(75)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(76)=(-1 + 2*x1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
   sb(77)=1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6
   sb(78)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
   sb(79)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
   sb(80)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(81)=(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(82)=(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
   sb(83)=1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6
 sb(84)=-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7
 sb(85)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*y1)
 sb(86)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*z1)
 sb(87)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*y1 + 6*y1**2)
 sb(88)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(89)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*z1 + 6*z1**2)
 sb(90)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(91)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(92)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(93)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(94)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(95)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(96)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(97)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(98)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(99)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
 sb(100)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(101)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(102)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(103)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(104)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(105)=(-1 + 2*x1)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)
 sb(106)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
 sb(107)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
 sb(108)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(109)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(110)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(111)=(-1 + 2*x1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(112)=-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7
 sb(113)=(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(-1 + 2*z1)
 sb(114)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(1 - 6*z1 + 6*z1**2)
 sb(115)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(116)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(117)=(1 - 6*y1 + 6*y1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(118)=(-1 + 2*y1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(119)=-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7 
 sb(164)=1 - 72*z1 + 1260*z1**2 - 9240*z1**3 + 34650*z1**4 - 72072*z1**5 + 84084*z1**6 - 51480*z1**7 + 12870*z1**8
 sb(163)=(-1 + 2*y1)*(-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7)
 sb(162)=(1 - 6*y1 + 6*y1**2)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(161)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(160)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(159)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(158)=(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(1 - 6*z1 + 6*z1**2)
 sb(157)=(-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7)*(-1 + 2*z1)
 sb(156)=1 - 72*y1 + 1260*y1**2 - 9240*y1**3 + 34650*y1**4 - 72072*y1**5 + 84084*y1**6 - 51480*y1**7 + 12870*y1**8
 sb(155)=(-1 + 2*x1)*(-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7)  
 sb(154)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(153)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(152)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(151)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(150)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(1 - 6*z1 + 6*z1**2)
 sb(149)=(-1 + 2*x1)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(-1 + 2*z1)
 sb(148)=(-1 + 2*x1)*(-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7)
 sb(147)=(1 - 6*x1 + 6*x1**2)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(146)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(145)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(144)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(143)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
 sb(142)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
 sb(141)=(1 - 6*x1 + 6*x1**2)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)
 sb(140)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(139)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(138)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(137)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(136)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(135)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
 sb(134)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(133)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(132)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(131)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(130)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(129)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(128)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(127)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(126)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(125)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(1 - 6*z1 + 6*z1**2)
 sb(124)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(123)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(1 - 6*y1 + 6*y1**2)
 sb(122)=(-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7)*(-1 + 2*z1)
 sb(121)=(-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7)*(-1 + 2*y1)
 sb(120)=1 - 72*x1 + 1260*x1**2 - 9240*x1**3 + 34650*x1**4 - 72072*x1**5 + 84084*x1**6 - 51480*x1**7 + 12870*x1**8
 case(9)
   !ninth order functions (10th-order of accuracy 165-219)
 sb(1)=-1 + 2*x1
    sb(2)=-1 + 2*y1
    sb(3)=-1 + 2*z1 
 sb(4)=1 - 6*x1 + 6*x1**2
 sb(5)=(-1 + 2*x1)*(-1 + 2*y1)
 sb(6)=(-1 + 2*x1)*(-1 + 2*z1)
 sb(7)=1 - 6*y1 + 6*y1**2 
 sb(8)=(-1 + 2*y1)*(-1 + 2*z1)
 sb(9)=1 - 6*z1 + 6*z1**2
 sb(10)=-1 + 12*x1 - 30*x1**2 + 20*x1**3
 sb(11)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)
 sb(12)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)
 sb(13)=-1 + 12*y1 - 30*y1**2 + 20*y1**3
 sb(14)=(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(15)=(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(16)=-1 + 12*z1 - 30*z1**2 + 20*z1**3
 sb(17)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(18)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*z1)
 sb(19)=(-1 + 2*x1)*(1 - 6*z1 + 6*z1**2)
 sb(20)=1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4
 sb(21)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)
 sb(22)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*z1)
 sb(23)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)
 sb(24)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(25)=(1 - 6*x1 + 6*x1**2)*(1 - 6*z1 + 6*z1**2)
 sb(26)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(27)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(28)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(29)=(-1 + 2*x1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(30)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(31)=(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(32)=(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(33)=1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4
 sb(34)=1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4
 sb(35)=-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5
 sb(36)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)
 sb(37)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*z1)
 sb(38)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)
 sb(39)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*z1 + 6*z1**2)
 sb(40)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(41)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(42)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(43)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(44)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(45)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(46)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(47)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(48)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(49)=-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5
 sb(50)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(51)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(52)=(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(53)=(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(54)=-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5
 sb(55)=(-1 + 2*x1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(56)=1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6
   sb(57)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)
   sb(58)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*z1)
   sb(59)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)
   sb(60)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(-1 + 2*z1)
   sb(61)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*z1 + 6*z1**2)
   sb(62)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
   sb(63)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
   sb(64)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
   sb(65)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(66)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
   sb(67)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
   sb(68)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
   sb(69)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(70)=(1 - 6*x1 + 6*x1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(71)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
   sb(72)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
   sb(73)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
   sb(74)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(75)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(76)=(-1 + 2*x1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
   sb(77)=1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6
   sb(78)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
   sb(79)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
   sb(80)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
   sb(81)=(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
   sb(82)=(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
   sb(83)=1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6
 sb(84)=-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7
 sb(85)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*y1)
 sb(86)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*z1)
 sb(87)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*y1 + 6*y1**2)
 sb(88)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(89)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*z1 + 6*z1**2)
 sb(90)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(91)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(92)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(93)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(94)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(95)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(96)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(97)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(98)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(99)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
 sb(100)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(101)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(102)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(103)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(104)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(105)=(-1 + 2*x1)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)
 sb(106)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
 sb(107)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
 sb(108)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(109)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(110)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(111)=(-1 + 2*x1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(112)=-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7
 sb(113)=(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(-1 + 2*z1)
 sb(114)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(1 - 6*z1 + 6*z1**2)
 sb(115)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(116)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(117)=(1 - 6*y1 + 6*y1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(118)=(-1 + 2*y1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(119)=-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7 
 sb(164)=1 - 72*z1 + 1260*z1**2 - 9240*z1**3 + 34650*z1**4 - 72072*z1**5 + 84084*z1**6 - 51480*z1**7 + 12870*z1**8
 sb(163)=(-1 + 2*y1)*(-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7)
 sb(162)=(1 - 6*y1 + 6*y1**2)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(161)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(160)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(159)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(158)=(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(1 - 6*z1 + 6*z1**2)
 sb(157)=(-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7)*(-1 + 2*z1)
 sb(156)=1 - 72*y1 + 1260*y1**2 - 9240*y1**3 + 34650*y1**4 - 72072*y1**5 + 84084*y1**6 - 51480*y1**7 + 12870*y1**8
 sb(155)=(-1 + 2*x1)*(-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7)  
 sb(154)=(-1 + 2*x1)*(-1 + 2*y1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(153)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(152)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(151)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(150)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(1 - 6*z1 + 6*z1**2)
 sb(149)=(-1 + 2*x1)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(-1 + 2*z1)
 sb(148)=(-1 + 2*x1)*(-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7)
 sb(147)=(1 - 6*x1 + 6*x1**2)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(146)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(145)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(144)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(143)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
 sb(142)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
 sb(141)=(1 - 6*x1 + 6*x1**2)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)
 sb(140)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(139)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(138)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(137)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(136)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(135)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
 sb(134)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(133)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(132)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(131)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(130)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(129)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(128)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(127)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(126)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(125)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(1 - 6*z1 + 6*z1**2)
 sb(124)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(123)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(1 - 6*y1 + 6*y1**2)
 sb(122)=(-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7)*(-1 + 2*z1)
 sb(121)=(-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7)*(-1 + 2*y1)
 sb(120)=1 - 72*x1 + 1260*x1**2 - 9240*x1**3 + 34650*x1**4 - 72072*x1**5 + 84084*x1**6 - 51480*x1**7 + 12870*x1**8  
 sb(219)=-1 + 90*z1 - 1980*z1**2 + 18480*z1**3 - 90090*z1**4 + 252252*z1**5 - 420420*z1**6 + 411840*z1**7 - 218790*z1**8 + 48620*z1**9
 sb(218)=(-1 + 2*y1)*(1 - 72*z1 + 1260*z1**2 - 9240*z1**3 + 34650*z1**4 - 72072*z1**5 + 84084*z1**6 - 51480*z1**7 + 12870*z1**8)
 sb(217)=(1 - 6*y1 + 6*y1**2)*(-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7)
 sb(216)=(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(215)=(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(214)=(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(213)=(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(212)=(-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7)*(1 - 6*z1 + 6*z1**2)
 sb(211)=(1 - 72*y1 + 1260*y1**2 - 9240*y1**3 + 34650*y1**4 - 72072*y1**5 + 84084*y1**6 - 51480*y1**7 + 12870*y1**8)*(-1 + 2*z1)
 sb(210)=-1 + 90*y1 - 1980*y1**2 + 18480*y1**3 - 90090*y1**4 + 252252*y1**5 - 420420*y1**6 + 411840*y1**7 - 218790*y1**8 + 48620*y1**9
 sb(209)=(-1 + 2*x1)*(1 - 72*z1 + 1260*z1**2 - 9240*z1**3 + 34650*z1**4 - 72072*z1**5 + 84084*z1**6 - 51480*z1**7 + 12870*z1**8)
 sb(208)=(-1 + 2*x1)*(-1 + 2*y1)*(-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7)
 sb(207)=(-1 + 2*x1)*(1 - 6*y1 + 6*y1**2)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(206)=(-1 + 2*x1)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(204)=(-1 + 2*x1)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(203)=(-1 + 2*x1)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(1 - 6*z1 + 6*z1**2)
 sb(202)=(-1 + 2*x1)*(-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7)*(-1 + 2*z1)
 sb(201)=(-1 + 2*x1)*(1 - 72*y1 + 1260*y1**2 - 9240*y1**3 + 34650*y1**4 - 72072*y1**5 + 84084*y1**6 - 51480*y1**7 + 12870*y1**8)
 sb(200)=(1 - 6*x1 + 6*x1**2)*(-1 + 56*z1 - 756*z1**2 + 4200*z1**3 - 11550*z1**4 + 16632*z1**5 - 12012*z1**6 + 3432*z1**7)
 sb(205)=(-1 + 2*x1)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(198)=(1 - 6*x1 + 6*x1**2)*(1 - 6*y1 + 6*y1**2)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(199)=(1 - 6*x1 + 6*x1**2)*(-1 + 2*y1)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(197)=(1 - 6*x1 + 6*x1**2)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(196)=(1 - 6*x1 + 6*x1**2)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(195)=(1 - 6*x1 + 6*x1**2)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(1 - 6*z1 + 6*z1**2)
 sb(194)=(1 - 6*x1 + 6*x1**2)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)*(-1 + 2*z1)
 sb(193)=(1 - 6*x1 + 6*x1**2)*(-1 + 56*y1 - 756*y1**2 + 4200*y1**3 - 11550*y1**4 + 16632*y1**5 - 12012*y1**6 + 3432*y1**7)
 sb(192)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 42*z1 + 420*z1**2 - 1680*z1**3 + 3150*z1**4 - 2772*z1**5 + 924*z1**6)
 sb(191)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 2*y1)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(190)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 6*y1 + 6*y1**2)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(189)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(188)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(1 - 6*z1 + 6*z1**2)
 sb(187)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)*(-1 + 2*z1)
 sb(186)=(-1 + 12*x1 - 30*x1**2 + 20*x1**3)*(1 - 42*y1 + 420*y1**2 - 1680*y1**3 + 3150*y1**4 - 2772*y1**5 + 924*y1**6)
 sb(185)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 30*z1 - 210*z1**2 + 560*z1**3 - 630*z1**4 + 252*z1**5)
 sb(184)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 2*y1)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(183)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 6*y1 + 6*y1**2)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(182)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(1 - 6*z1 + 6*z1**2)
 sb(181)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)*(-1 + 2*z1)
 sb(180)=(1 - 20*x1 + 90*x1**2 - 140*x1**3 + 70*x1**4)*(-1 + 30*y1 - 210*y1**2 + 560*y1**3 - 630*y1**4 + 252*y1**5)
 sb(179)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 20*z1 + 90*z1**2 - 140*z1**3 + 70*z1**4)
 sb(178)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 2*y1)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(177)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 6*y1 + 6*y1**2)*(1 - 6*z1 + 6*z1**2)
 sb(176)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)*(-1 + 2*z1)
 sb(175)=(-1 + 30*x1 - 210*x1**2 + 560*x1**3 - 630*x1**4 + 252*x1**5)*(1 - 20*y1 + 90*y1**2 - 140*y1**3 + 70*y1**4)
 sb(173)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 12*z1 - 30*z1**2 + 20*z1**3)
 sb(174)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 2*y1)*(1 - 6*z1 + 6*z1**2)
 sb(172)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(1 - 6*y1 + 6*y1**2)*(-1 + 2*z1)
 sb(171)=(1 - 42*x1 + 420*x1**2 - 1680*x1**3 + 3150*x1**4 - 2772*x1**5 + 924*x1**6)*(-1 + 12*y1 - 30*y1**2 + 20*y1**3)
 sb(170)=(-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7)*(1 - 6*z1 + 6*z1**2)
 sb(169)=(-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7)*(-1 + 2*y1)*(-1 + 2*z1)
 sb(168)=(-1 + 56*x1 - 756*x1**2 + 4200*x1**3 - 11550*x1**4 + 16632*x1**5 - 12012*x1**6 + 3432*x1**7)*(1 - 6*y1 + 6*y1**2)
 sb(167)=(1 - 72*x1 + 1260*x1**2 - 9240*x1**3 + 34650*x1**4 - 72072*x1**5 + 84084*x1**6 - 51480*x1**7 + 12870*x1**8)*(-1 + 2*z1)
 sb(166)=(1 - 72*x1 + 1260*x1**2 - 9240*x1**3 + 34650*x1**4 - 72072*x1**5 + 84084*x1**6 - 51480*x1**7 + 12870*x1**8)*(-1 + 2*y1)
 sb(165)=-1 + 90*x1 - 1980*x1**2 + 18480*x1**3 - 90090*x1**4 + 252252*x1**5 - 420420*x1**6 + 411840*x1**7 - 218790*x1**8 + 48620*x1**9
 
     end select
   end if
   
   
   
    if (poly.eq.4)then
   
   
   
   select case(number)
   case(1)
    !first order functions (2nd-order of accuracy 3)
    sb(1)=x1/hxc
    sb(2)=y1/hxc
    sb(3)=z1/hxc
    case(2)
! second order functions (3rd-order of accuracy 4-9)
        sb(1)=x1/hxc
        sb(2)=y1/hxc
        sb(3)=z1/hxc
        sb(4)=x1**2/(2.0d0*hxc**2)
        sb(5)=(x1*y1)/hxc**2
        sb(6)=(x1*z1)/hxc**2
        sb(7)=y1**2/(2.0d0*hxc**2)
        sb(8)=(y1*z1)/hxc**2
        sb(9)=z1**2/(2.0d0*hxc**2)
 case(3)
! third order functions (4th-order of accuracy  10-19)
        sb(1)=x1/hxc
        sb(2)=y1/hxc
        sb(3)=z1/hxc
        sb(4)=x1**2/(2.0d0*hxc**2)
        sb(5)=(x1*y1)/hxc**2
        sb(6)=(x1*z1)/hxc**2
        sb(7)=y1**2/(2.0d0*hxc**2)
        sb(8)=(y1*z1)/hxc**2
        sb(9)=z1**2/(2.0d0*hxc**2)
        sb(10)=x1**3/(6.0d0*hxc**3)
        sb(11)=(x1**2*y1)/(2.0d0*hxc**3)
        sb(12)=(x1**2*z1)/(2.0d0*hxc**3)
        sb(13)=(x1*y1**2)/(2.0d0*hxc**3)
        sb(14)=(x1*z1**2)/(2.0d0*hxc**3)
        sb(15)=(x1*y1*z1)/hxc**3
        sb(16)=y1**3/(6.0d0*hxc**3)
        sb(17)=(y1**2*z1)/(2.0d0*hxc**3)
        sb(18)=(y1*z1**2)/(2.0d0*hxc**3)
        sb(19)=z1**3/(6.0d0*hxc**3)
 
 case(4)
! fourth order functions (5th-order of accuracy 20-34)
        sb(1)=x1/hxc
        sb(2)=y1/hxc
        sb(3)=z1/hxc
        sb(4)=x1**2/(2.0d0*hxc**2)
        sb(5)=(x1*y1)/hxc**2
        sb(6)=(x1*z1)/hxc**2
        sb(7)=y1**2/(2.0d0*hxc**2)
        sb(8)=(y1*z1)/hxc**2
        sb(9)=z1**2/(2.0d0*hxc**2)
        sb(10)=x1**3/(6.0d0*hxc**3)
        sb(11)=(x1**2*y1)/(2.0d0*hxc**3)
        sb(12)=(x1**2*z1)/(2.0d0*hxc**3)
        sb(13)=(x1*y1**2)/(2.0d0*hxc**3)
        sb(14)=(x1*z1**2)/(2.0d0*hxc**3)
        sb(15)=(x1*y1*z1)/hxc**3
        sb(16)=y1**3/(6.0d0*hxc**3)
        sb(17)=(y1**2*z1)/(2.0d0*hxc**3)
        sb(18)=(y1*z1**2)/(2.0d0*hxc**3)
        sb(19)=z1**3/(6.0d0*hxc**3)
        sb(20)=x1**4/(24.0d0*hxc**4)
        sb(21)=(x1**3*y1)/(6.0d0*hxc**4)
        sb(22)=(x1**3*z1)/(6.0d0*hxc**4)
        sb(23)=(x1**2*y1**2)/(4.0d0*hxc**4)
        sb(24)=(x1**2*z1**2)/(4.0d0*hxc**4)
        sb(25)=(x1**2*y1*z1)/(2.0d0*hxc**4)
        sb(26)=(x1*y1**3)/(6.0d0*hxc**4)
        sb(27)=(x1*z1**3)/(6.0d0*hxc**4)
        sb(28)=(x1*y1**2*z1)/(2.0d0*hxc**4)
        sb(29)=(x1*y1*z1**2)/(2.0d0*hxc**4)
        sb(30)=y1**4/(24.0d0*hxc**4)
        sb(31)=(y1**3*z1)/(6.0d0*hxc**4)
        sb(32)=(y1**2*z1**2)/(4.0d0*hxc**4)
        sb(33)=(y1*z1**3)/(6.0d0*hxc**4)
        sb(34)=z1**4/(24.0d0*hxc**4)
 case(5)
! fifth order functions (6th-order of accuracy 35-55)
  sb(1)=x1/hxc
        sb(2)=y1/hxc
        sb(3)=z1/hxc
        sb(4)=x1**2/(2.0d0*hxc**2)
        sb(5)=(x1*y1)/hxc**2
        sb(6)=(x1*z1)/hxc**2
        sb(7)=y1**2/(2.0d0*hxc**2)
        sb(8)=(y1*z1)/hxc**2
        sb(9)=z1**2/(2.0d0*hxc**2)
        sb(10)=x1**3/(6.0d0*hxc**3)
        sb(11)=(x1**2*y1)/(2.0d0*hxc**3)
        sb(12)=(x1**2*z1)/(2.0d0*hxc**3)
        sb(13)=(x1*y1**2)/(2.0d0*hxc**3)
        sb(14)=(x1*z1**2)/(2.0d0*hxc**3)
        sb(15)=(x1*y1*z1)/hxc**3
        sb(16)=y1**3/(6.0d0*hxc**3)
        sb(17)=(y1**2*z1)/(2.0d0*hxc**3)
        sb(18)=(y1*z1**2)/(2.0d0*hxc**3)
        sb(19)=z1**3/(6.0d0*hxc**3)
        sb(20)=x1**4/(24.0d0*hxc**4)
        sb(21)=(x1**3*y1)/(6.0d0*hxc**4)
        sb(22)=(x1**3*z1)/(6.0d0*hxc**4)
        sb(23)=(x1**2*y1**2)/(4.0d0*hxc**4)
        sb(24)=(x1**2*z1**2)/(4.0d0*hxc**4)
        sb(25)=(x1**2*y1*z1)/(2.0d0*hxc**4)
        sb(26)=(x1*y1**3)/(6.0d0*hxc**4)
        sb(27)=(x1*z1**3)/(6.0d0*hxc**4)
        sb(28)=(x1*y1**2*z1)/(2.0d0*hxc**4)
        sb(29)=(x1*y1*z1**2)/(2.0d0*hxc**4)
        sb(30)=y1**4/(24.0d0*hxc**4)
        sb(31)=(y1**3*z1)/(6.0d0*hxc**4)
        sb(32)=(y1**2*z1**2)/(4.0d0*hxc**4)
        sb(33)=(y1*z1**3)/(6.0d0*hxc**4)
        sb(34)=z1**4/(24.0d0*hxc**4)
        sb(35)=x1**5/(120.0d0*hxc**5)
        sb(36)=(x1**4*y1)/(24.0d0*hxc**5)
        sb(37)=(x1**4*z1)/(24.0d0*hxc**5)
        sb(38)=(x1**3*y1**2)/(12.0d0*hxc**5)
        sb(39)=(x1**3*z1**2)/(12.0d0*hxc**5)
        sb(40)=(x1**3*y1*z1)/(6.0d0*hxc**5)
        sb(41)=(x1**2*y1**3)/(12.0d0*hxc**5)
        sb(42)=(x1**2*y1**2*z1)/(4.0d0*hxc**5)
        sb(43)=(x1**2*y1*z1**2)/(4.0d0*hxc**5)
        sb(44)=(x1**2*z1**3)/(12.0d0*hxc**5)
        sb(45)=(x1*y1**4)/(24.0d0*hxc**5)
        sb(46)=(x1*y1**3*z1)/(6.0d0*hxc**5)
        sb(47)=(x1*y1**2*z1**2)/(4.0d0*hxc**5)
        sb(48)=(x1*y1*z1**3)/(6.0d0*hxc**5)
        sb(49)=(x1*z1**4)/(24.0d0*hxc**5)
        sb(50)=y1**5/(120.0d0*hxc**5)
        sb(51)=(y1**4*z1)/(24.0d0*hxc**5)
        sb(52)=(y1**3*z1**2)/(12.0d0*hxc**5)
        sb(53)=(y1**2*z1**3)/(12.0d0*hxc**5)
        sb(54)=(y1*z1**4)/(24.0d0*hxc**5)
        sb(55)=z1**5/(120.0d0*hxc**5)
 case(6)
   ! sixth order functions (7th-order of accuracy 56-83)
   sb(1)=x1/hxc
        sb(2)=y1/hxc
        sb(3)=z1/hxc
        sb(4)=x1**2/(2.0d0*hxc**2)
        sb(5)=(x1*y1)/hxc**2
        sb(6)=(x1*z1)/hxc**2
        sb(7)=y1**2/(2.0d0*hxc**2)
        sb(8)=(y1*z1)/hxc**2
        sb(9)=z1**2/(2.0d0*hxc**2)
        sb(10)=x1**3/(6.0d0*hxc**3)
        sb(11)=(x1**2*y1)/(2.0d0*hxc**3)
        sb(12)=(x1**2*z1)/(2.0d0*hxc**3)
        sb(13)=(x1*y1**2)/(2.0d0*hxc**3)
        sb(14)=(x1*z1**2)/(2.0d0*hxc**3)
        sb(15)=(x1*y1*z1)/hxc**3
        sb(16)=y1**3/(6.0d0*hxc**3)
        sb(17)=(y1**2*z1)/(2.0d0*hxc**3)
        sb(18)=(y1*z1**2)/(2.0d0*hxc**3)
        sb(19)=z1**3/(6.0d0*hxc**3)
        sb(20)=x1**4/(24.0d0*hxc**4)
        sb(21)=(x1**3*y1)/(6.0d0*hxc**4)
        sb(22)=(x1**3*z1)/(6.0d0*hxc**4)
        sb(23)=(x1**2*y1**2)/(4.0d0*hxc**4)
        sb(24)=(x1**2*z1**2)/(4.0d0*hxc**4)
        sb(25)=(x1**2*y1*z1)/(2.0d0*hxc**4)
        sb(26)=(x1*y1**3)/(6.0d0*hxc**4)
        sb(27)=(x1*z1**3)/(6.0d0*hxc**4)
        sb(28)=(x1*y1**2*z1)/(2.0d0*hxc**4)
        sb(29)=(x1*y1*z1**2)/(2.0d0*hxc**4)
        sb(30)=y1**4/(24.0d0*hxc**4)
        sb(31)=(y1**3*z1)/(6.0d0*hxc**4)
        sb(32)=(y1**2*z1**2)/(4.0d0*hxc**4)
        sb(33)=(y1*z1**3)/(6.0d0*hxc**4)
        sb(34)=z1**4/(24.0d0*hxc**4)
        sb(35)=x1**5/(120.0d0*hxc**5)
        sb(36)=(x1**4*y1)/(24.0d0*hxc**5)
        sb(37)=(x1**4*z1)/(24.0d0*hxc**5)
        sb(38)=(x1**3*y1**2)/(12.0d0*hxc**5)
        sb(39)=(x1**3*z1**2)/(12.0d0*hxc**5)
        sb(40)=(x1**3*y1*z1)/(6.0d0*hxc**5)
        sb(41)=(x1**2*y1**3)/(12.0d0*hxc**5)
        sb(42)=(x1**2*y1**2*z1)/(4.0d0*hxc**5)
        sb(43)=(x1**2*y1*z1**2)/(4.0d0*hxc**5)
        sb(44)=(x1**2*z1**3)/(12.0d0*hxc**5)
        sb(45)=(x1*y1**4)/(24.0d0*hxc**5)
        sb(46)=(x1*y1**3*z1)/(6.0d0*hxc**5)
        sb(47)=(x1*y1**2*z1**2)/(4.0d0*hxc**5)
        sb(48)=(x1*y1*z1**3)/(6.0d0*hxc**5)
        sb(49)=(x1*z1**4)/(24.0d0*hxc**5)
        sb(50)=y1**5/(120.0d0*hxc**5)
        sb(51)=(y1**4*z1)/(24.0d0*hxc**5)
        sb(52)=(y1**3*z1**2)/(12.0d0*hxc**5)
        sb(53)=(y1**2*z1**3)/(12.0d0*hxc**5)
        sb(54)=(y1*z1**4)/(24.0d0*hxc**5)
        sb(55)=z1**5/(120.0d0*hxc**5)
        sb(56)=x1**6/(720.0d0*hxc**6)
        sb(57)=(x1**5*y1)/(120.0d0*hxc**6)
        sb(58)=(x1**5*z1)/(120.0d0*hxc**6)
        sb(59)=(x1**4*y1**2)/(48.0d0*hxc**6)
        sb(60)=(x1**4*y1*z1)/(24.0d0*hxc**6)
        sb(61)=(x1**4*z1**2)/(48.0d0*hxc**6)
        sb(62)=(x1**3*y1**3)/(36.0d0*hxc**6)
        sb(63)=(x1**3*y1**2*z1)/(12.0d0*hxc**6)
        sb(64)=(x1**3*y1*z1**2)/(12.0d0*hxc**6)
        sb(65)=(x1**3*z1**3)/(36.0d0*hxc**6)
        sb(66)=(x1**2*y1**4)/(48.0d0*hxc**6)
        sb(67)=(x1**2*y1**3*z1)/(12.0d0*hxc**6)
        sb(68)=(x1**2*y1**2*z1**2)/(8.0d0*hxc**6)
        sb(69)=(x1**2*y1*z1**3)/(12.0d0*hxc**6)
        sb(70)=(x1**2*z1**4)/(48.0d0*hxc**6)
        sb(71)=(x1*y1**5)/(120.0d0*hxc**6)
        sb(72)=(x1*y1**4*z1)/(24.0d0*hxc**6)
        sb(73)=(x1*y1**3*z1**2)/(12.0d0*hxc**6)
        sb(74)=(x1*y1**2*z1**3)/(12.0d0*hxc**6)
        sb(75)=(x1*y1*z1**4)/(24.0d0*hxc**6)
        sb(76)=(x1*z1**5)/(120.0d0*hxc**6)
        sb(77)=y1**6/(720.0d0*hxc**6)
        sb(78)=(y1**5*z1)/(120.0d0*hxc**6)
        sb(79)=(y1**4*z1**2)/(48.0d0*hxc**6)
        sb(80)=(y1**3*z1**3)/(36.0d0*hxc**6)
        sb(81)=(y1**2*z1**4)/(48.0d0*hxc**6)
        sb(82)=(y1*z1**5)/(120.0d0*hxc**6)
        sb(83)=z1**6/(720.0d0*hxc**6)
   
   
    end select
    end if

   select case (icompwrt)
    
    
    case(-2)

        oov=1.0d0/(ielem_totvolume(iconsidered))
        basis_rec(1:number_of_dog) = sb(1:number_of_dog)-((integ_basis_dg_value(1:number_of_dog,iconsidered))*oov)
    
    case(0)
        oov=1.0d0/(rec_volume(1,1,iconsidered))
        basis_rec(1:number_of_dog)=sb(1:number_of_dog)-((integ_basis_value(1:number_of_dog,iconsidered))*oov)
    
    case(1)
        oov=1.0d0/(rec_volume(1,1,iconsidered))
        basis_rec(1:number_of_dog)=sb(1:number_of_dog)-((integ_basis_valuec(1:number_of_dog,iconsidered))*oov)
    
    end select





    

end function basis_rec



function basis_rec2d(n,x1,y1,number,iconsidered,number_of_dog,icompwrt)
implicit none
!> @brief
!> this function returns the value of the basis function for a specific polynomial order and coordinates in 2d \n
!> requires: x1, y1: coordinates of basis evaluation wrt ?; number: order of basis; iconsidered: considered cell?; number_of_dog: number of degrees of freedom
! number and number_of_dog redundant?
#ifdef gpu
!!$omp declare target disabled: array-return basis routine replaced by basis_rec2d_fill
#endif
integer,intent(in)::n
integer,intent(in)::number,iconsidered,number_of_dog,icompwrt
real,intent(in)::x1,y1
integer::i_node,i_qp,n_qp
real::oov,hxc
real,dimension(number_of_dog)::basis_rec2d
real,dimension(number_of_dog)::sb
sb=zero


oov=1.0d0/(rec_volume(1,1,iconsidered))

hxc=(sqrt(ielem_totvolume(iconsidered)))




select case(poly)
case(1) ! generic
    select case(number) ! order of basis
    
    case(1)
    !first order
    sb(1)=x1
    sb(2)=y1
    
    case(2)
    !second order
    sb(1)=x1
    sb(2)=y1
    sb(3)=x1*x1
    sb(4) = x1*y1
    sb(5)= y1*y1
    
    case(3)
    !third order
    sb(1)=x1
    sb(2)=y1
    sb(3)=x1*x1
    sb(4) = x1*y1
    sb(5)= y1*y1
    sb(6)= x1*x1*x1
    sb(7) = x1*x1*y1
    sb(8)= x1*y1*y1
    sb(9)= y1*y1*y1
    
    case(4)
    !fourth order
    sb(1)=x1
    sb(2)=y1
    sb(3)=x1*x1
    sb(4) = x1*y1
    sb(5)= y1*y1
    sb(6)= x1*x1*x1
    sb(7) = x1*x1*y1
    sb(8)= x1*y1*y1
    sb(9)= y1*y1*y1
    sb(10)= x1*x1*x1*x1
    sb(11)= x1*x1*x1*y1
    sb(12)= x1*x1*y1*y1
    sb(13)= x1*y1*y1*y1
    sb(14)= y1*y1*y1*y1
    
    case(5)
    !fifth order
    sb(1)=x1
    sb(2)=y1
    sb(3)=x1*x1
    sb(4) = x1*y1
    sb(5)= y1*y1
    sb(6)= x1*x1*x1
    sb(7) = x1*x1*y1
    sb(8)= x1*y1*y1
    sb(9)= y1*y1*y1
    sb(10)= x1*x1*x1*x1
    sb(11)= x1*x1*x1*y1
    sb(12)= x1*x1*y1*y1
    sb(13)= x1*y1*y1*y1
    sb(14)= y1*y1*y1*y1
    sb(15) = x1*x1*x1*x1*x1
    sb(16)= x1*x1*x1*x1*y1
    sb(17)= x1*x1*x1*y1*y1
    sb(18)= x1*x1*y1*y1*y1
    sb(19)= x1*y1*y1*y1*y1
    sb(20)= y1*y1*y1*y1*y1
    
    case(6)
    !sixth order
     sb(1)=x1
    sb(2)=y1
    sb(3)=x1*x1
    sb(4) = x1*y1
    sb(5)= y1*y1
    sb(6)= x1*x1*x1
    sb(7) = x1*x1*y1
    sb(8)= x1*y1*y1
    sb(9)= y1*y1*y1
    sb(10)= x1*x1*x1*x1
    sb(11)= x1*x1*x1*y1
    sb(12)= x1*x1*y1*y1
    sb(13)= x1*y1*y1*y1
    sb(14)= y1*y1*y1*y1
    sb(15) = x1*x1*x1*x1*x1
    sb(16)= x1*x1*x1*x1*y1
    sb(17)= x1*x1*x1*y1*y1
    sb(18)= x1*x1*y1*y1*y1
    sb(19)= x1*y1*y1*y1*y1
    sb(20)= y1*y1*y1*y1*y1
    sb(21) = x1**6
    sb(22) = x1**5*y1
    sb(23)= x1**4*y1**2
    sb(24)= x1**3*y1**3
    sb(25)= x1**2*y1**4
    sb(26)= x1*y1**5
    sb(27)= y1**6

    
    end select
    
case(2) ! legendre
    select case(number)
    case(1)
        !first order functions (2nd-order of accuracy 3)
        sb(1)=-1.0d0 + 2.0d0*x1
        sb(2)=-1.0d0 + 2.0d0*y1

    case(2)
        ! second order functions (3rd-order of accuracy 4-9)
        sb(1)=-1.0d0 + 2.0d0*x1
        sb(2)=-1.0d0 + 2.0d0*y1
        sb(3)=1.0d0 - 6.0d0*x1 + 6.0d0*x1**2
        sb(4)=sb(1)*sb(2)
        sb(5)=1.0d0 - 6.0d0*y1 + 6.0d0*y1**2

    case(3)
        ! third order functions (4th-order of accuracy  10-19)
        sb(1)=-1.0d0 + 2.0d0*x1
        sb(2)=-1.0d0 + 2.0d0*y1
        sb(3)=1.0d0 - 6.0d0*x1 + 6.0d0*x1**2
        sb(4)=sb(1)*sb(2)
        sb(5)=1.0d0 - 6.0d0*y1 + 6.0d0*y1**2 
        sb(6)=-1.0d0 + 12.0d0*x1 - 30.0d0*x1**2 + 20.0d0*x1**3
        sb(7)=sb(3)*sb(2)
        sb(8)=sb(1)*sb(5)
        sb(9)=-1.0d0 + 12.0d0*y1 - 30.0d0*y1**2 + 20.0d0*y1**3  
    end select
    



case(4) ! !taylor
    select case(number)


    case (1)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
           
    
    case (2)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
            sb(3)=x1**2/(2.*hxc**2)
            sb(4)=(x1*y1)/hxc**2
            sb(5)=y1**2/(2.*hxc**2)
            
            
    case (3)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
            sb(3)=x1**2/(2.0d0*hxc**2)
            sb(4)=(x1*y1)/hxc**2
            sb(5)=y1**2/(2.0d0*hxc**2)
            sb(6)=x1**3/(6.0d0*hxc**3)
            sb(7)=(x1**2*y1)/(2.0d0*hxc**3)
            sb(8)=(x1*y1**2)/(2.0d0*hxc**3)
            sb(9)=y1**3/(6.0d0*hxc**3)
            
            
    case (4)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
            sb(3)=x1**2/(2.0d0*hxc**2)
            sb(4)=(x1*y1)/hxc**2
            sb(5)=y1**2/(2.0d0*hxc**2)
            sb(6)=x1**3/(6.0d0*hxc**3)
            sb(7)=(x1**2*y1)/(2.0d0*hxc**3)
            sb(8)=(x1*y1**2)/(2.0d0*hxc**3)
            sb(9)=y1**3/(6.0d0*hxc**3)
            sb(10)=x1**4/(24.0d0*hxc**4)
            sb(11)=(x1**3*y1)/(6.0d0*hxc**4)
            sb(12)=(x1**2*y1**2)/(4.0d0*hxc**4)
            sb(13)=(x1*y1**3)/(6.0d0*hxc**4)
            sb(14)=y1**4/(24.0d0*hxc**4)
           
            
    case (5)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
            sb(3)=x1**2/(2.*hxc**2)
            sb(4)=(x1*y1)/hxc**2
            sb(5)=y1**2/(2.*hxc**2)
            sb(6)=x1**3/(6.*hxc**3)
            sb(7)=(x1**2*y1)/(2.*hxc**3)
            sb(8)=(x1*y1**2)/(2.*hxc**3)
            sb(9)=y1**3/(6.*hxc**3)
            sb(10)=x1**4/(24.*hxc**4)
            sb(11)=(x1**3*y1)/(6.*hxc**4)
            sb(12)=(x1**2*y1**2)/(4.*hxc**4)
            sb(13)=(x1*y1**3)/(6.*hxc**4)
            sb(14)=y1**4/(24.*hxc**4)
            sb(15)=x1**5/(120.*hxc**5)
            sb(16)=(x1**4*y1)/(24.*hxc**5)
            sb(17)=(x1**3*y1**2)/(12.*hxc**5)
            sb(18)=(x1**2*y1**3)/(12.*hxc**5)
            sb(19)=(x1*y1**4)/(24.*hxc**5)
            sb(20)=y1**5/(120.*hxc**5)
            
            
    case (6)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
            sb(3)=x1**2/(2.*hxc**2)
            sb(4)=(x1*y1)/hxc**2
            sb(5)=y1**2/(2.*hxc**2)
            sb(6)=x1**3/(6.*hxc**3)
            sb(7)=(x1**2*y1)/(2.*hxc**3)
            sb(8)=(x1*y1**2)/(2.*hxc**3)
            sb(9)=y1**3/(6.*hxc**3)
            sb(10)=x1**4/(24.*hxc**4)
            sb(11)=(x1**3*y1)/(6.*hxc**4)
            sb(12)=(x1**2*y1**2)/(4.*hxc**4)
            sb(13)=(x1*y1**3)/(6.*hxc**4)
            sb(14)=y1**4/(24.*hxc**4)
            sb(15)=x1**5/(120.*hxc**5)
            sb(16)=(x1**4*y1)/(24.*hxc**5)
            sb(17)=(x1**3*y1**2)/(12.*hxc**5)
            sb(18)=(x1**2*y1**3)/(12.*hxc**5)
            sb(19)=(x1*y1**4)/(24.*hxc**5)
            sb(20)=y1**5/(120.*hxc**5)
            sb(21)=x1**6/(720.*hxc**6)
            sb(22)=(x1**5*y1)/(120.*hxc**6)
            sb(23)=(x1**4*y1**2)/(48.*hxc**6)
            sb(24)=(x1**3*y1**3)/(36.*hxc**6)
            sb(25)=(x1**2*y1**4)/(48.*hxc**6)
            sb(26)=(x1*y1**5)/(120.*hxc**6)
            sb(27)=y1**6/(720.*hxc**6)
           
            
    case (7)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
            sb(3)=x1**2/(2.*hxc**2)
            sb(4)=(x1*y1)/hxc**2
            sb(5)=y1**2/(2.*hxc**2)
            sb(6)=x1**3/(6.*hxc**3)
            sb(7)=(x1**2*y1)/(2.*hxc**3)
            sb(8)=(x1*y1**2)/(2.*hxc**3)
            sb(9)=y1**3/(6.*hxc**3)
            sb(10)=x1**4/(24.*hxc**4)
            sb(11)=(x1**3*y1)/(6.*hxc**4)
            sb(12)=(x1**2*y1**2)/(4.*hxc**4)
            sb(13)=(x1*y1**3)/(6.*hxc**4)
            sb(14)=y1**4/(24.*hxc**4)
            sb(15)=x1**5/(120.*hxc**5)
            sb(16)=(x1**4*y1)/(24.*hxc**5)
            sb(17)=(x1**3*y1**2)/(12.*hxc**5)
            sb(18)=(x1**2*y1**3)/(12.*hxc**5)
            sb(19)=(x1*y1**4)/(24.*hxc**5)
            sb(20)=y1**5/(120.*hxc**5)
            sb(21)=x1**6/(720.*hxc**6)
            sb(22)=(x1**5*y1)/(120.*hxc**6)
            sb(23)=(x1**4*y1**2)/(48.*hxc**6)
            sb(24)=(x1**3*y1**3)/(36.*hxc**6)
            sb(25)=(x1**2*y1**4)/(48.*hxc**6)
            sb(26)=(x1*y1**5)/(120.*hxc**6)
            sb(27)=y1**6/(720.*hxc**6)
            sb(28)=x1**7/(5040.*hxc**7)
            sb(29)=(x1**6*y1)/(720.*hxc**7)
            sb(30)=(x1**5*y1**2)/(240.*hxc**7)
            sb(31)=(x1**4*y1**3)/(144.*hxc**7)
            sb(32)=(x1**3*y1**4)/(144.*hxc**7)
            sb(33)=(x1**2*y1**5)/(240.*hxc**7)
            sb(34)=(x1*y1**6)/(720.*hxc**7)
            sb(35)=y1**7/(5040.*hxc**7)
            
            
    case (8)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
            sb(3)=x1**2/(2.*hxc**2)
            sb(4)=(x1*y1)/hxc**2
            sb(5)=y1**2/(2.*hxc**2)
            sb(6)=x1**3/(6.*hxc**3)
            sb(7)=(x1**2*y1)/(2.*hxc**3)
            sb(8)=(x1*y1**2)/(2.*hxc**3)
            sb(9)=y1**3/(6.*hxc**3)
            sb(10)=x1**4/(24.*hxc**4)
            sb(11)=(x1**3*y1)/(6.*hxc**4)
            sb(12)=(x1**2*y1**2)/(4.*hxc**4)
            sb(13)=(x1*y1**3)/(6.*hxc**4)
            sb(14)=y1**4/(24.*hxc**4)
            sb(15)=x1**5/(120.*hxc**5)
            sb(16)=(x1**4*y1)/(24.*hxc**5)
            sb(17)=(x1**3*y1**2)/(12.*hxc**5)
            sb(18)=(x1**2*y1**3)/(12.*hxc**5)
            sb(19)=(x1*y1**4)/(24.*hxc**5)
            sb(20)=y1**5/(120.*hxc**5)
            sb(21)=x1**6/(720.*hxc**6)
            sb(22)=(x1**5*y1)/(120.*hxc**6)
            sb(23)=(x1**4*y1**2)/(48.*hxc**6)
            sb(24)=(x1**3*y1**3)/(36.*hxc**6)
            sb(25)=(x1**2*y1**4)/(48.*hxc**6)
            sb(26)=(x1*y1**5)/(120.*hxc**6)
            sb(27)=y1**6/(720.*hxc**6)
            sb(28)=x1**7/(5040.*hxc**7)
            sb(29)=(x1**6*y1)/(720.*hxc**7)
            sb(30)=(x1**5*y1**2)/(240.*hxc**7)
            sb(31)=(x1**4*y1**3)/(144.*hxc**7)
            sb(32)=(x1**3*y1**4)/(144.*hxc**7)
            sb(33)=(x1**2*y1**5)/(240.*hxc**7)
            sb(34)=(x1*y1**6)/(720.*hxc**7)
            sb(35)=y1**7/(5040.*hxc**7)
            sb(36)=x1**8/(40320.*hxc**8)
            sb(37)=(x1**7*y1)/(5040.*hxc**8)
            sb(38)=(x1**6*y1**2)/(1440.*hxc**8)
            sb(39)=(x1**5*y1**3)/(720.*hxc**8)
            sb(40)=(x1**4*y1**4)/(576.*hxc**8)
            sb(41)=(x1**3*y1**5)/(720.*hxc**8)
            sb(42)=(x1**2*y1**6)/(1440.*hxc**8)
            sb(43)=(x1*y1**7)/(5040.*hxc**8)
            sb(44)=y1**8/(40320.*hxc**8)
            
            
    case (9)
    
            sb(1)=x1/hxc
            sb(2)=y1/hxc
            sb(3)=x1**2/(2.*hxc**2)
            sb(4)=(x1*y1)/hxc**2
            sb(5)=y1**2/(2.*hxc**2)
            sb(6)=x1**3/(6.*hxc**3)
            sb(7)=(x1**2*y1)/(2.*hxc**3)
            sb(8)=(x1*y1**2)/(2.*hxc**3)
            sb(9)=y1**3/(6.*hxc**3)
            sb(10)=x1**4/(24.*hxc**4)
            sb(11)=(x1**3*y1)/(6.*hxc**4)
            sb(12)=(x1**2*y1**2)/(4.*hxc**4)
            sb(13)=(x1*y1**3)/(6.*hxc**4)
            sb(14)=y1**4/(24.*hxc**4)
            sb(15)=x1**5/(120.*hxc**5)
            sb(16)=(x1**4*y1)/(24.*hxc**5)
            sb(17)=(x1**3*y1**2)/(12.*hxc**5)
            sb(18)=(x1**2*y1**3)/(12.*hxc**5)
            sb(19)=(x1*y1**4)/(24.*hxc**5)
            sb(20)=y1**5/(120.*hxc**5)
            sb(21)=x1**6/(720.*hxc**6)
            sb(22)=(x1**5*y1)/(120.*hxc**6)
            sb(23)=(x1**4*y1**2)/(48.*hxc**6)
            sb(24)=(x1**3*y1**3)/(36.*hxc**6)
            sb(25)=(x1**2*y1**4)/(48.*hxc**6)
            sb(26)=(x1*y1**5)/(120.*hxc**6)
            sb(27)=y1**6/(720.*hxc**6)
            sb(28)=x1**7/(5040.*hxc**7)
            sb(29)=(x1**6*y1)/(720.*hxc**7)
            sb(30)=(x1**5*y1**2)/(240.*hxc**7)
            sb(31)=(x1**4*y1**3)/(144.*hxc**7)
            sb(32)=(x1**3*y1**4)/(144.*hxc**7)
            sb(33)=(x1**2*y1**5)/(240.*hxc**7)
            sb(34)=(x1*y1**6)/(720.*hxc**7)
            sb(35)=y1**7/(5040.*hxc**7)
            sb(36)=x1**8/(40320.*hxc**8)
            sb(37)=(x1**7*y1)/(5040.*hxc**8)
            sb(38)=(x1**6*y1**2)/(1440.*hxc**8)
            sb(39)=(x1**5*y1**3)/(720.*hxc**8)
            sb(40)=(x1**4*y1**4)/(576.*hxc**8)
            sb(41)=(x1**3*y1**5)/(720.*hxc**8)
            sb(42)=(x1**2*y1**6)/(1440.*hxc**8)
            sb(43)=(x1*y1**7)/(5040.*hxc**8)
            sb(44)=y1**8/(40320.*hxc**8)
            sb(45)=x1**9/(362880.*hxc**9)
            sb(46)=(x1**8*y1)/(40320.*hxc**9)
            sb(47)=(x1**7*y1**2)/(10080.*hxc**9)
            sb(48)=(x1**6*y1**3)/(4320.*hxc**9)
            sb(49)=(x1**5*y1**4)/(2880.*hxc**9)
            sb(50)=(x1**4*y1**5)/(2880.*hxc**9)
            sb(51)=(x1**3*y1**6)/(4320.*hxc**9)
            sb(52)=(x1**2*y1**7)/(10080.*hxc**9)
            sb(53)=(x1*y1**8)/(40320.*hxc**9)
            sb(54)=y1**9/(362880.*hxc**9)
            
    
    
    




    end select




end select
   
   
    select case (icompwrt)
    
    
    case(-2)

        oov=1.0d0/(ielem_totvolume(iconsidered))
        basis_rec2d(1:number_of_dog) = sb(1:number_of_dog)-(integ_basis_dg_value(1:number_of_dog,iconsidered)*oov)
        
        
    
        
    
    case(0)
        oov=1.0d0/(rec_volume(1,1,iconsidered))
        basis_rec2d(1:number_of_dog)=sb(1:number_of_dog)-(integ_basis_value(1:number_of_dog,iconsidered)*oov)
    
    case(1)
        oov=1.0d0/(rec_volume(1,1,iconsidered))
        basis_rec2d(1:number_of_dog)=sb(1:number_of_dog)-(integ_basis_valuec(1:number_of_dog,iconsidered)*oov)
    
    end select
    

end function basis_rec2d




end module basis
