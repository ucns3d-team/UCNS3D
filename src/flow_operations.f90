MODULE FLOW_OPERATIONS
USE DECLARATION
USE MPIINFO
USE TRANSFORM


IMPLICIT NONE


 CONTAINS
 
 

 subroutine MRFSWITCH(N,ICONSIDERED,FACEX,POINTX,pox,poy)

  !> @brief
!> This subroutine  check if the element is on the rotational/stationary reference frame and update the MRF_ORIGIN and SRF_VELOCITY SRF accordingly

implicit none
INTEGER, INTENT(IN)::N,ICONSIDERED,FACEX,POINTX
! TYPE(LOCAL_RECON3),ALLOCATABLE,DIMENSION(:),INTENT(INOUT)::ILOCAL_RECON3
!output:
!ILOCAL_RECON3%MRF_ORIGIN; ILOCAL_RECON3%MRF_VELOCITY; ILOCAL_RECON3%ROTVEL, ILOCAL_RECON3%MRF
real, dimension(3) ::MRF_ORIGIN, MRF_VELOCITY,ROTVEL
integer:: ROTFRAME_ON
INTEGER::NINV
real,dimension(1:dimensiona),intent(inout)::POX,POY
!internal variables
real, dimension(3) :: P1P2, PC, POPC,PO,PGP,POZ !element coordinates, roation_axys, Cylinder_center_coordinates, vector_element_center, rotational velocity at gaussian points, Gausian points coordinates
real :: d1, d2, r1, theta, dPOPC

PO=(POX(1:3))
PGP=(POY(1:3))
!body
DO NINV=1,NROTORS

PC(1:3)= (point1_GL(NINV,1:3)+point2_GL(NINV,1:3))/2  !center of cYlinder
P1P2(1:3)=point2_GL(NINV,1:3)-point1_GL(NINV,1:3)          !axysvector
POPC(1:3)=PO(1:3)-PC(1:3)              ! vector elelement-centre
dPOPC=((PO(1)-PC(1))**2+(PO(2)-PC(2))**2+(PO(3)-PC(3))**2)**0.5 !distance between element and center

theta= ACOS((dot_product(POPC,P1P2))/(sqrt(POPC(1)**2+POPC(2)**2+POPC(3)**2)*sqrt(P1P2(1)**2+P1P2(2)**2+P1P2(3)**2))) !angle between element vector and axys
d2=  dPOPC*abs(cos(theta))
r1=dPOPC*abs(sin(theta))
d1=((point1_GL(NINV,1)-PC(1))**2+(point1_GL(NINV,2)-PC(2))**2+(point1_GL(NINV,3)-PC(3))**2)**0.5

if ((d1.ge.d2).and.(r1.le.Radius_GL(NINV))) then
   ROTFRAME_ON=1
    MRF_ORIGIN(1:3)=PC(1:3)
    POX(1:3)=PGP(1:3)-MRF_ORIGIN(1:3)
    MRF_VELOCITY(1:3)=MRF_ROT_GL(NINV)*(P1P2)/(P1P2(1)**2+P1P2(2)**2+P1P2(3)**2)**0.5
!     SRF_VELOCITY(1)=0.0
!     SRF_VELOCITY(2)=MRF_ROT_GL
!     SRF_VELOCITY(3)=0.0
    POY(1:3)=MRF_VELOCITY(1:3)
    ROTVEL(1:3)=VECT_FUNCTION(POX,POY)



    GO TO 606
else
    MRF_ORIGIN(1:3)=0.0
    ROTFRAME_ON=0
    MRF_VELOCITY(1:3)=0.0
    ROTVEL(1:3)=0.0
end if

END DO

606 CONTINUE

ILOCAL_RECON3(ICONSIDERED)%MRF_ORIGIN=MRF_ORIGIN
ILOCAL_RECON3(ICONSIDERED)%MRF_VELOCITY=MRF_VELOCITY
ILOCAL_RECON3(ICONSIDERED)%ROTVEL(FACEX,POINTX,1:3)=ROTVEL
ILOCAL_RECON3(ICONSIDERED)%MRF=ROTFRAME_ON






end subroutine MRFSWITCH





SUBROUTINE MULTISPECIES_MIXTURES_RG(LEFTV,MP_mu_mix,MP_ktr_mix,MP_kve,MP_D_eff,MP_HTR,MP_HVIB,GAMMAL)
implicit none
real,dimension(1:nof_Variables),intent(in)::leftv
real,intent(inout)::MP_mu_mix,MP_ktr_mix,MP_kve,gammal
real,dimension(1:nof_species),intent(inout)::MP_HTR,MP_HVIB,MP_D_eff
real,dimension(1:nof_species)::RGS_htr_i, RGS_hvib_i
real:: RGS_Ttr, RGS_Tv, RGS_P_pa, RGS_rho,RGS_mu_mix,MP_PINFl
real, dimension(1:nof_SPECIES):: RGS_x, RGS_Y, RGS_Deff, RGS_mu_i, RGS_ktr_i,tmp
real, dimension(1:nof_SPECIES,1:nof_SPECIES):: RGS_Dij
real ::RGS_ktr_mix, RGS_kve,sumx,sumY
real,dimension(1:nof_variables)::temp_vect1,temp_vect2,temp_vect3,temp_vect4
integer::rg_i,k
real, parameter :: epsY = 1.0d-20   ! floor for tiny/negative noise
real, parameter :: epsx = 1.0d-20   ! floor for tiny/negative noise
real, parameter :: tiny = 1.0d-300  ! protect against division by zero

temp_vect1=leftv  !copy left vector
temp_vect2=leftv

CALL CONS2DIV(N,temp_vect1,MP_PINFl,gammal)

!now this vector holds the following variables rho, u,v,w,Ttt,tvib,y_1,y_2,y_3,y_4,y_5 etc

RGS_Ttr=temp_vect1(dimensiona+2)
RGS_Tv=temp_vect1(dimensiona+3)

RGS_Y(1:nof_SPECIES)=temp_vect1(dimensiona+4:nof_Variables)



! convert to mole-fraction *proportional* values
do k = 1, nof_species
    tmp(k) = RGS_Y(k) / RG_molm(k)
end do

! normalise to get actual mole fractions
sumx = 0.d0
do k = 1, nof_species
    sumx = sumx + tmp(k)
end do

if (sumx > tiny) then
    do k = 1, nof_species
        RGS_x(k) = tmp(k) / sumx
    end do
else
    ! fallback if something is seriously wrong
    do k = 1, nof_species
        RGS_x(k) = 1.d0 / nof_species
    end do
end if






call cons2prim(n,temp_Vect2,mp_pinfl,gammal)

RGS_P_pa=temp_Vect2(dimensiona+2)
RGS_rho=temp_Vect2(1)



! RGS_Ttr=RGS_Ttr
! RGS_Tv=RGS_Tv


    CALL COMPUTE_REAL_GAS_DIFFUSION(RGS_Ttr, RGS_Tv, RGS_P_pa, RGS_rho, RGS_x, RGS_Y, &
                                     RGS_Dij, RGS_Deff, RGS_mu_i, RGS_mu_mix, RGS_ktr_i, RGS_ktr_mix, RGS_kve, &
                                     RGS_htr_i, RGS_hvib_i)

       MP_mu_mix=RGS_mu_mix
       MP_ktr_mix=RGS_ktr_mix
       MP_kve=RGS_kve
       MP_D_eff=RGS_Deff
       MP_HTR(1:nof_SPECIES)=RGS_htr_i(1:nof_species)
       MP_Hvib(1:nof_SPECIES)=RGS_hvib_i(1:nof_species)



END SUBROUTINE MULTISPECIES_MIXTURES_RG






  subroutine COMPUTE_REAL_GAS_DIFFUSION( RGS_Ttr, RGS_Tv, RGS_P_pa, RGS_rho, RGS_x, RGS_Y, &
                                       RGS_Dij, RGS_Deff, RGS_mu_i, RGS_mu_mix, &
                                       RGS_ktr_i, RGS_ktr_mix, RGS_kve, &
                                       RGS_htr_i, RGS_hvib_i )

  implicit none
  real, intent(in)  :: RGS_Ttr, RGS_Tv, RGS_P_pa, RGS_rho
  real, intent(in)  :: RGS_x(5), RGS_Y(5)
  real, intent(out) :: RGS_Dij(5,5), RGS_Deff(5)
  real, intent(out) :: RGS_mu_i(5), RGS_mu_mix
  real, intent(out) :: RGS_ktr_i(5), RGS_ktr_mix, RGS_kve
  real, intent(out) :: RGS_htr_i(5), RGS_hvib_i(5)

  ! ---- Binary diffusion ----
  call compute_binary_diffusion( RGS_Ttr, RGS_P_pa, RGS_Dij )

  ! ---- Mixture-average diffusion coefficient ----
  call compute_effective_diffusion( RGS_x, RGS_Dij, RGS_Deff )

  ! ---- Species viscosity (Blottner) ----
  call blottner_mu_species( RGS_Ttr, RGS_mu_i )

  ! ---- Mixture viscosity (Wilke) ----
  call wilke_mixture_viscosity( RGS_x, RGS_mu_i, RGS_mu_mix )

  ! ---- Species thermal conductivity (Eucken) ----
  call eucken_ktr_species( RGS_mu_i, RGS_ktr_i )

  ! ---- Mixture translational conductivity ----
  call mason_saxena_ktr_mixture( RGS_x, RGS_ktr_i, RGS_ktr_mix )

  ! ---- Vibrational conductivity ----
  call vibrational_conductivity( RGS_rho, RGS_Y, RGS_Deff, RGS_Tv, RGS_kve )

  ! ---- Species translational and vibrational enthalpy ----
  call htr_air5 ( RGS_Ttr, RGS_htr_i )
  call hvib_air5( RGS_Tv , RGS_hvib_i )

end subroutine COMPUTE_REAL_GAS_DIFFUSION





  function omega11_neufeld(RGS_Tstar) result(RGS_omega11)
  implicit none
  real, intent(in) :: RGS_Tstar
  real             :: RGS_omega11
  real, parameter  :: RGS_A=1.06036d0, RGS_B=0.15610d0, RGS_C=0.19300d0, RGS_D=0.47635d0, &
                      RGS_E=1.03587d0, RGS_F=1.52996d0, RGS_G=1.76474d0, RGS_H=3.89411d0
  real :: Tstar_eff

  Tstar_eff = max(RGS_Tstar, 1.0d-6)   ! avoid 0^(-B) and tiny T*

  RGS_omega11 = RGS_A / (Tstar_eff**RGS_B)                              &
              + RGS_C*exp(-RGS_D*Tstar_eff)                             &
              + RGS_E*exp(-RGS_F*Tstar_eff)                             &
              + RGS_G*exp(-RGS_H*Tstar_eff)

end function omega11_neufeld




subroutine compute_binary_diffusion(RGS_T, RGS_P_pa, RGS_Dij)
  implicit none
  real, intent(in)  :: RGS_T, RGS_P_pa
  real, intent(out) :: RGS_Dij(5,5)

  integer :: i, j
  real :: sig_ij, eps_ij, Tstar, omega, P_atm, denom
  real, parameter :: tiny = 1d-30

  RGS_Dij = 0.0d0

  ! Pressure in atm
  P_atm = max(RGS_P_pa / RGS_Pa_per_atm, 1d-12)

  do i = 1, 5
     do j = i+1, 5

        sig_ij = 0.5d0 * (RGS_sigmaA(i) + RGS_sigmaA(j))
        eps_ij = sqrt( RGS_eps_over_k(i) * RGS_eps_over_k(j) )

        ! Non-dimensional temperature
        Tstar  = max(RGS_T / eps_ij, 1d-6)
        omega  = omega11_neufeld(Tstar)

        denom = P_atm * sig_ij**2 * omega * sqrt(1.d0/RGS_Mg(i) + 1.d0/RGS_Mg(j))
        denom = max(denom, tiny)

        RGS_Dij(i,j) = (0.001858d0 * RGS_T**1.5d0) / denom * RGS_cm2s_to_m2s
        RGS_Dij(j,i) = RGS_Dij(i,j)

     end do
  end do

end subroutine compute_binary_diffusion










 subroutine compute_effective_diffusion(RGS_x, RGS_Dij, RGS_Deff)
  implicit none
  real, intent(in)  :: RGS_x(5)
  real, intent(in)  :: RGS_Dij(5,5)
  real, intent(out) :: RGS_Deff(5)

  integer :: i, j
  real :: sumj, one_minus_xi, sumx
  real :: xloc(5)
  real, parameter :: tiny = 1d-20
  real, parameter :: Dmin = 1d-10   ! lower bound on diffusivity [m^2/s]

  ! Work on a local copy to enforce positivity and normalisation


xloc(:) = RGS_x(:)

! enforce positivity
do i = 1, 5
    if (xloc(i) < 0.d0) xloc(i) = 0.d0
end do

! enforce normalisation
sumx = sum(xloc)
if (sumx > tiny) then
    xloc(:) = xloc(:) / sumx
else
    xloc(:) = 1.d0 / 5.d0
end if

  do i = 1, 5
     sumj = 0.0d0

     do j = 1, 5
        if (j /= i) then
           sumj = sumj + xloc(j) / max(RGS_Dij(i,j), tiny)
        end if
     end do

     one_minus_xi = max(1.0d0 - xloc(i), 0.0d0)

     if (one_minus_xi <= tiny .or. sumj <= tiny) then
        ! Pure or nearly pure gas: use small floor (essentially no diffusive transport)
        RGS_Deff(i) = Dmin
     else
        RGS_Deff(i) = one_minus_xi / sumj
        RGS_Deff(i) = max(RGS_Deff(i), Dmin)
     end if
  end do

end subroutine compute_effective_diffusion



   subroutine blottner_mu_species(RGS_T, RGS_mu)
  implicit none
  real, intent(in)  :: RGS_T
  real, intent(out) :: RGS_mu(5)
  integer :: i
  real :: Tlog

  Tlog = log10(max(RGS_T, 50.d0))

  do i = 1, 5
     RGS_mu(i) = 1.d-7 * 10.d0**( RGS_aB(i)*Tlog*Tlog + RGS_bB(i)*Tlog + RGS_cB(i) )
     RGS_mu(i) = max(RGS_mu(i), 1d-12)
  end do
end subroutine blottner_mu_species

   subroutine wilke_mixture_viscosity(RGS_x, RGS_mu_i, RGS_mu_mix)
  implicit none
  real, intent(in)  :: RGS_x(5), RGS_mu_i(5)
  real, intent(out) :: RGS_mu_mix
  integer :: i, j
  real :: phi_ij, denom

  RGS_mu_mix = 0.d0

  do i = 1, 5
     denom = 0.d0

     do j = 1, 5
        if (i == j) then
            denom = denom + RGS_x(j)
        else
            phi_ij = ( 1.d0 + sqrt(RGS_mu_i(i)/RGS_mu_i(j)) * (RG_MOLM(j)/RG_MOLM(i))**0.25 )**2 &
                     / ( sqrt(8.d0) * sqrt(1.d0 + RG_MOLM(i)/RG_MOLM(j)) )

            denom = denom + RGS_x(j) * phi_ij
        end if
     end do

     RGS_mu_mix = RGS_mu_mix + RGS_x(i) * RGS_mu_i(i) / max(denom, 1d-20)
  end do

end subroutine wilke_mixture_viscosity

   subroutine cp_tr_species(RGS_Cp_tr)
  implicit none
    real, intent(out) :: RGS_Cp_tr(5)  ! J/kg-K
    integer :: RGS_i
    real :: RGS_Rspec
    do RGS_i = 1, 5
      RGS_Rspec = RGS_Ru / RG_MOLM(RGS_i)
      if (RGS_i <= 3) then
        RGS_Cp_tr(RGS_i) = 3.5d0 * RGS_Rspec   ! diatomics: 7/2 R
      else
        RGS_Cp_tr(RGS_i) = 2.5d0 * RGS_Rspec   ! atoms:     5/2 R
      end if
    end do
  end subroutine cp_tr_species

   subroutine eucken_ktr_species(RGS_mu_i, RGS_ktr_i)
  implicit none
  real, intent(in)  :: RGS_mu_i(5)
  real, intent(out) :: RGS_ktr_i(5)

  integer :: i
  real :: Rspec, Cp

  do i = 1, 5
     Rspec = RGS_Ru / RG_MOLM(i)

     if (i <= 3) then
         Cp = 3.5d0 * Rspec          ! diatomic
         RGS_ktr_i(i) = RGS_mu_i(i) * (Cp + 1.25d0*Rspec)
     else
         Cp = 2.5d0 * Rspec          ! atomic
         RGS_ktr_i(i) = RGS_mu_i(i) * (Cp + 1.50d0*Rspec)
     end if

  end do

end subroutine eucken_ktr_species

  subroutine mason_saxena_ktr_mixture(RGS_x, RGS_ktr_i, RGS_ktr_mix)
  implicit none
  real, intent(in)  :: RGS_x(5), RGS_ktr_i(5)
  real, intent(out) :: RGS_ktr_mix

  integer :: i, j
  real :: psi_ij, denom

  RGS_ktr_mix = 0.d0

  do i = 1, 5
     denom = 0.d0

     do j = 1, 5
        if (i == j) then
            denom = denom + RGS_x(j)
        else
            psi_ij = (1.d0 + sqrt(RGS_ktr_i(i)/RGS_ktr_i(j)) * (RG_MOLM(j)/RG_MOLM(i))**0.25 )**2 &
                      / ( sqrt(8.d0) * sqrt(1.d0 + RG_MOLM(i)/RG_MOLM(j)) )

            denom = denom + RGS_x(j) * psi_ij
        end if
     end do

     RGS_ktr_mix = RGS_ktr_mix + RGS_x(i) * RGS_ktr_i(i) / max(denom, 1d-20)

  end do

end subroutine mason_saxena_ktr_mixture

function RGS_cv_vibrational_diatomic(RGS_Tv, RGS_theta) result(RGS_cv)
  implicit none
  real, intent(in) :: RGS_Tv, RGS_theta
  real             :: RGS_cv, RGS_x, ex

  if (RGS_Tv <= 1.0d0 .or. RGS_theta <= 0.0d0) then
    RGS_cv = 0.0d0
  else
    RGS_x = RGS_theta / RGS_Tv
    if (RGS_x > 60.0d0) then
      ! Asymptotic form for large x: cv ∝ x^2 e^(-x)
      RGS_cv = RGS_Ru * (RGS_x*RGS_x) * exp(-RGS_x)
    else
      ex     = exp(RGS_x)
      RGS_cv = RGS_Ru * (RGS_x*RGS_x) * ex / ( (ex - 1.0d0)**2 )
    end if
  end if
end function RGS_cv_vibrational_diatomic


function RGS_hvib_species(RGS_Tv, RGS_theta, RGS_Mi) result(RGS_hv)
  implicit none
  real, intent(in) :: RGS_Tv, RGS_theta, RGS_Mi
  real             :: RGS_hv, RGS_x, RGS_Ri, ex

  if (RGS_theta <= 0.0d0 .or. RGS_Tv <= 1.0d0) then
    RGS_hv = 0.0d0
    return
  end if

  RGS_Ri = RGS_Ru / RGS_Mi
  RGS_x  = RGS_theta / RGS_Tv

  if (RGS_x > 60.0d0) then
     ! Asymptotic form: h_v ≈ R_i θ e^(-x)
     RGS_hv = RGS_Ri * RGS_theta * exp(-RGS_x)
  else
     ex     = exp(RGS_x)
     RGS_hv = RGS_Ri * RGS_theta / (ex - 1.0d0)
  end if

end function RGS_hvib_species





 subroutine hvib_air5(RGS_Tv, RGS_hvib_i)
  implicit none
  real, intent(in)  :: RGS_Tv
  real, intent(out) :: RGS_hvib_i(5)     ! J/kg
  integer :: RGS_i

  do RGS_i = 1, 5
    RGS_hvib_i(RGS_i) = RGS_hvib_species(RGS_Tv, RG_thetaG(RGS_i), RG_MOLM(RGS_i))
  end do

end subroutine hvib_air5






subroutine htr_air5(RGS_Ttr, RGS_htr_i)
  implicit none
  real, intent(in)  :: RGS_Ttr
  real, intent(out) :: RGS_htr_i(5)
  real :: Cp(5), h0_mass(5)
  integer :: i
  real, parameter :: Tref = 298.15d0

  call cp_tr_species(Cp)

  do i = 1, 5
     h0_mass(i)   = RG_hzero(i) / RG_MOLM(i)
     RGS_htr_i(i) = h0_mass(i) + Cp(i) * (RGS_Ttr - Tref)
  end do

end subroutine htr_air5







subroutine vibrational_conductivity(RGS_rho, RGS_Y, RGS_Deff, RGS_Tv, RGS_kve)
  implicit none

  real, intent(in)  :: RGS_rho, RGS_Y(5), RGS_Deff(5), RGS_Tv
  real, intent(out) :: RGS_kve

  real :: Cv_vib_mass(5)
  integer :: i

  do i = 1, 5
     if (RG_thetaG(i) > 0.0d0) then
        Cv_vib_mass(i) = RGS_cv_vibrational_diatomic(RGS_Tv, RG_thetaG(i)) / RG_MOLM(i)
     else
        Cv_vib_mass(i) = 0.0d0
     end if
  end do

  RGS_kve = 0.0d0
  do i = 1, 5
     RGS_kve = RGS_kve + RGS_rho * RGS_Y(i) * RGS_Deff(i) * Cv_vib_mass(i)
  end do

  RGS_kve = max(RGS_kve, 0.0d0)

end subroutine vibrational_conductivity












SUBROUTINE MULTISPECIES_MIXTURES(LEFTV,MP_Temp,MP_mu_mix,MP_k_mix,MP_Cp_mix,gammal)
  IMPLICIT NONE
  REAL, INTENT(INOUT)  :: MP_TEMP    ! Temperature in Kelvin
  REAL, INTENT(IN),DIMENSION(1:NOF_variables)  :: LEFTV      ! Vector of conserved variables
  REAL, INTENT(OUT) :: MP_mu_mix  ! Mixture viscosity [Pa.s]
  REAL, INTENT(OUT) :: MP_k_mix   ! Mixture thermal conductivity [W/m.K]
  REAL, INTENT(OUT) :: MP_Cp_mix  ! Mixture Cp [J/mol.K]
  REAL, INTENT(OUT) :: GAMMAL     ! total
  REAL::MOLAR_SUM,MASS_SUM
  integer::i,j,k
  ! Constants
  REAL,DIMENSION(1:NOF_SPECIEs) :: MP_VISCL,MP_DENOM
  REAL,DIMENSION(1:NOF_SPECIES) ::MP_MOLE_FRACTION,MP_MASS_FRACTION
  REAL,DIMENSION(1:NOF_SPECIEs) :: MP_Cp
  REAL,DIMENSION(1:NOF_SPECIEs) :: ML_LAML,MP_DENOL
  REAL,DIMENSION(1:NOF_SPECIES,1:NOF_SPECIES) :: MP_phi





  ! Molar masses
  !MP_M(1)  = 2.016E-3    ! kg/mol   !defined in file
  !MP_M(2)  = 28.97E-3    ! kg/mol   !defined in file

  !MASS SUM
  MASS_SUM=ZERO
  DO I=1,NOF_SPECIES
      MP_MASS_FRACTION(I)=LEFTV(NOF_VARIABLES-NOF_SPECIES-1+I)
      MASS_SUM=MASS_SUM+MP_MASS_FRACTION(I)
  END DO
  DO I=1,NOF_SPECIES
        MP_MASS_FRACTION(I)=MP_MASS_FRACTION(I)/MASS_SUM
  END DO



  !Mole fraction of EACH SPECIES
  MOLAR_SUM=ZERO
  DO I=1,NOF_SPECIES
      MP_MOLE_FRACTION(I)=MP_MASS_FRACTION(i)/MP_M(I)
      MOLAR_SUM=MOLAR_SUM+MP_MOLE_FRACTION(I)
  END DO
  DO I=1,NOF_SPECIES
        MP_MOLE_FRACTION(I)=MP_MOLE_FRACTION(I)/MOLAR_SUM
  END DO



  !NOW GET THE TEMPERATURE

  CALL GET_TEMP_JANAF(LEFTV,MP_TEMP,MP_MASS_FRACTION,MP_MOLE_FRACTION)




  ! --- Compute individual Cp values [J/mol.K]
  DO I=1,NOF_SPECIES
   MP_Cp(I)=RGS_Ru*(MP_JANAF(I,j,1) + MP_JANAF(I,j,2)*mp_Temp + MP_JANAF(I,j,3)*mp_Temp**2 + MP_JANAF(I,j,4)*mp_Temp**3 + MP_JANAF(I,j,5)*mp_Temp**4)
  END DO
  ! --- Compute individual viscosities [Pa.s]
  DO I=1,NOF_SPECIES
    MP_VISCL(I)=MP_BROK_A(I) * mp_Temp**MP_BROK_B(I) + MP_BROK_C(I)
  END DO

  ! --- Compute interaction terms
  DO i = 1, NOF_SPECIES
    DO j = 1, NOF_SPECIES
      IF (i == j) THEN
        MP_phi(i,j) = 1.0D0
      ELSE
        MP_phi(i,j) = (1.0D0 + SQRT(MP_VISCL(i)/MP_VISCL(j)) * (MP_M(j)/MP_M(i))**0.25D0)**2 / &
                   SQRT(8.0D0 * (1.0D0 + MP_M(i)/MP_M(j)))
      END IF
    END DO
  END DO

  ! Wilke’s mixture viscosity
  MP_mu_mix = 0.0D0
  DO i = 1, NOF_SPECIES
    MP_DENOM(i) = 0.0D0
    DO j = 1, NOF_SPECIES
      MP_DENOM(i) = MP_DENOM(i) + MP_MOLE_FRACTION(j) * MP_phi(i,j)
    END DO
    MP_mu_mix = MP_mu_mix + MP_MOLE_FRACTION(i) * MP_VISCL(i) / MP_DENOM(i)
  END DO

  ! Compute thermal conductivity via Eucken
  DO i = 1, NOF_SPECIES
    ML_LAML(i) = (MP_Cp(i) + 1.25D0 * RGS_Ru) * MP_VISCL(i)
  END DO


  MP_k_mix = 0.0D0
  DO i = 1, nof_species
    MP_DENOL(i) = 0.0D0
    DO j = 1, nof_species
      MP_DENOL(i) = MP_DENOL(i) + MP_MOLE_FRACTION(j) * MP_phi(i,j)
    END DO
    MP_k_mix = MP_k_mix + MP_MOLE_FRACTION(i) * ML_LAML(i) / MP_DENOL(i)
  END DO

  ! Mixture Cp
  MP_Cp_mix = 0.0D0
  DO i = 1, nof_species
    MP_Cp_mix = MP_Cp_mix + MP_MOLE_FRACTION(i) * MP_Cp(i)
  END DO

  !  Mixture gamma
  gammaL = MP_Cp_mix / (MP_Cp_mix - RGS_Ru)




END SUBROUTINE MULTISPECIES_MIXTURES













SUBROUTINE MULTISPECIES_temp(LEFTV,MP_Temp)
  IMPLICIT NONE
  REAL, INTENT(INOUT)  :: MP_TEMP    ! Temperature in Kelvin
  REAL, INTENT(IN),DIMENSION(1:NOF_variables)  :: LEFTV      ! Vector of conserved variables
  REAL::MOLAR_SUM,MASS_SUM
  integer::i,j,k
  ! Constants
  REAL,DIMENSION(1:NOF_SPECIEs) :: MP_VISCL,MP_DENOM
  REAL,DIMENSION(1:NOF_SPECIES) ::MP_MOLE_FRACTION,MP_MASS_FRACTION
  REAL,DIMENSION(1:NOF_SPECIEs) :: MP_Cp
  REAL,DIMENSION(1:NOF_SPECIEs) :: ML_LAML,MP_DENOL
  REAL,DIMENSION(1:NOF_SPECIES,1:NOF_SPECIES) :: MP_phi





  ! Molar masses
  !MP_M(1)  = 2.016E-3    ! kg/mol   !PLEASE DEFINE IN FILE
  !MP_M(2)  = 28.97E-3    ! kg/mol   !PLEASE DEFINE IN FILE

  !MASS SUM
  MASS_SUM=ZERO
  DO I=1,NOF_SPECIES
      MP_MASS_FRACTION(I)=LEFTV(NOF_VARIABLES-NOF_SPECIES-1+I)
      MASS_SUM=MASS_SUM+MP_MASS_FRACTION(I)
  END DO
  DO I=1,NOF_SPECIES
        MP_MASS_FRACTION(I)=MP_MASS_FRACTION(I)/MASS_SUM
  END DO



  !Mole fraction of EACH SPECIES
  MOLAR_SUM=ZERO
  DO I=1,NOF_SPECIES
      MP_MOLE_FRACTION(I)=MP_MASS_FRACTION(i)/MP_M(I)
      MOLAR_SUM=MOLAR_SUM+MP_MOLE_FRACTION(I)
  END DO
  DO I=1,NOF_SPECIES
        MP_MOLE_FRACTION(I)=MP_MOLE_FRACTION(I)/MOLAR_SUM
  END DO



  !NOW GET THE TEMPERATURE

  CALL GET_TEMP_JANAF(LEFTV,MP_TEMP,MP_MASS_FRACTION,MP_MOLE_FRACTION)









END SUBROUTINE MULTISPECIES_temp




  subroutine normalize_MF(MPC_MP)
    implicit none
    !! Ensure mole fractions sum to 1 (robust to roundoff)
    real,DIMENSION(1:NOF_SPECIES),intent(inout) :: MPC_MP
    real:: s
    integer :: k
    s = ZERO
    do k=1,NOF_SPECIES
    s = s + max(ZERO, MPC_MP(k))
    end do
    if (s > ZERO) then
      do k=1,NOF_SPECIES
      MPC_MP(k) = max(ZERO, MPC_MP(k))/s
      end do
    else
      do k=1,NOF_SPECIES
      MPC_MP(k) = ZERO
      end do
    end if
  end subroutine normalize_MF

  subroutine MOLEF_to_MASSF(MPC_MP,MPC_MASS)
    implicit none
    !! Convert mole fractions to mass fractions (for mass-specific energy)
    real,DIMENSION(1:NOF_SPECIES), intent(in)  :: MPC_MP
    real,DIMENSION(1:NOF_SPECIES), intent(out) :: MPC_MASS
    real :: denom
    integer :: s
    denom = ZERO
    do s=1,NOF_SPECIES
      denom = denom + MPC_MP(s)*MP_M(s)
    end do
    if (denom <= zero) then
      MPC_MASS = zero
    else
      do s=1,NOF_SPECIES
        MPC_MASS(s) = MPC_MP(s)*MP_M(s) / denom
      end do
    end if
  end subroutine MOLEF_to_MASSF

   subroutine species_cp_h_sensible(mp_Temp, is, mp_cpl, mp_hl)
    implicit none
    !! JANAF 5-coeffs, sensible enthalpy with h(298.15K)=0 for each species
    !!   cp/R = a1 + a2*T + a3*T^2 + a4*T^3 + a5*T^-2
    !!    h/R = a1*T + a2*T^2/2 + a3*T^3/3 + a4*T^4/4 - a5/T  (then shifted so h(298.15)=0)
    real, intent(in)  :: mp_Temp
    integer,  intent(in)  :: is
    real, intent(out) :: mp_cpl, mp_hl
    integer :: r
    real:: Rs, mp_tt,mp_t2,mp_t3, mp_invT, mp_invT2, mp_a1,mp_a2,mp_a3,mp_a4,mp_a5, mp_href, mp_Tref

    Rs = RGS_Ru / MP_M(is)
    r  = merge(1,2, mp_Temp <= MP_Tmid_in(is))
    mp_a1=MP_JANAF(1,r,is); mp_a2=MP_JANAF(2,r,is); mp_a3=MP_JANAF(3,r,is); mp_a4=MP_JANAF(4,r,is); mp_a5=MP_JANAF(5,r,is)

     mp_tt=mp_temp; mp_t2=mp_tt*mp_tt; mp_t3=mp_t2*mp_tt; mp_invT=1.0d0/mp_tt; mp_invT2=mp_invT*mp_invT

    mp_cpl = Rs * ( mp_a1 + mp_a2*mp_tt + mp_a3*mp_t2 + mp_a4*mp_t3 + mp_a5*mp_invT2 )

    ! raw enthalpy
    mp_hl  = Rs * ( mp_a1*mp_tt + 0.5d0*mp_a2*mp_t2 + (mp_a3/3.0d0)*mp_t3 + 0.25d0*mp_a4*mp_t2*mp_t2 - mp_a5*mp_invT )

    !only for sensible enthalpies
    ! subtract reference to exclude formation/constant offsets: h(298.15 K) = 0
    mp_Tref = 298.15d0
    mp_href = Rs * ( mp_a1*mp_Tref + 0.5d0*mp_a2*mp_Tref**2 + (mp_a3/3.0d0)*mp_Tref**3 + 0.25d0*mp_a4*mp_Tref**4 - mp_a5/mp_Tref )
    mp_hl = mp_hl - mp_href
  end subroutine species_cp_h_sensible

  subroutine mixture_props(MP_TEMP, mp_mol_x, MP_RMIX, MP_cp_mix, MP_cv_mix, MP_e_mix)
  implicit none
    !! Mixture properties with mass-weighted energies; energies are sensible
    real, intent(in)  :: MP_TEMP
    real,DIMENSION(1:NOF_SPECIES), intent(in)  :: mp_mol_x
    real, intent(out) :: MP_Rmix, MP_cp_mix, MP_cv_mix, MP_e_mix
    real, DIMENSION(1:NOF_SPECIES) :: MPC_MASS
    real::MP_cp_s, MP_h_s
    integer :: s
    call MOLEF_to_MASSF(mp_mol_x,MPC_MASS)
    MP_Rmix  = ZERO
    MP_cp_mix= ZERO
    MP_e_mix = ZERO
    do s=1,NOF_SPECIES
      if (MPC_MASS(s) <= ZERO) cycle
      MP_Rmix = MP_Rmix + MPC_MASS(s) * (RGS_Ru/mp_M(s))
      call species_cp_h_sensible(MP_TEMP, s, MP_cp_s, MP_h_s)
      MP_e_mix  = MP_e_mix  + MPC_MASS(s) * (MP_h_s - (RGS_Ru/mp_M(s))*MP_TEMP)  ! e = h - Rs*T (sensible e)
      MP_cp_mix = MP_cp_mix + MPC_MASS(s) * MP_cp_s
    end do
    MP_cv_mix = MP_cp_mix - MP_Rmix
  end subroutine mixture_props


  subroutine normalize_X(X)
  implicit none
    !! Ensure mole fractions sum to 1 (robust to roundoff)
    real, intent(inout) :: X(nof_species)
    real :: s
    integer :: k
    s = 0.0d0
    do k=1,Nof_species; s = s + max(0.0d0, X(k)); end do
    if (s > 0.0d0) then
      do k=1,Nof_species; X(k) = max(0.0d0, X(k))/s; end do
    else
      do k=1,Nof_species; X(k) = 0.0d0; end do
    end if
  end subroutine normalize_X



  subroutine GET_TEMP_JANAF(LEFTV, MP_TEMP,MP_MASS_FRACTION,MP_MOLE_FRACTION)
   implicit none
    !! Recover temperature from conserved variables using mole fractions X
    !!
    !! Inputs:
    !!   rho   [kg/m^3], u,v,w [m/s], rhoE [J/m^3], X_in(NS) mole fractions (sum≈1)
    !! Optional:
    !!   T_low,T_high [K] bracket (default 150..20000 clipped to species ranges)
    !!   tol residual on |e_mix - e_target| [J/kg] (default ~1e-6*|e_target|+1e-3)
    !!   max_iter (default 60)
    !!
    !! Outputs:
    !!   T [K], info: 0 ok; 1 bracket width tol; 2 iter limit; <0 bracketing issue
    real, dimension(1:nof_Variables),intent(in):: leftv
    real, intent(out) :: mp_temp
    REAL,DIMENSION(1:NOF_SPECIES),intent(in) ::MP_MOLE_FRACTION,MP_MASS_FRACTION
    integer:: mp_info,s
    real:: mp_tol,mp_SUM,u,v,w
    integer:: mp_max_iter,i,j,k
    real,dimension(1:nof_SPECIES)::mp_mol_x
    real:: mp_Rmix, mp_cp_mix, mp_cv_mix, mp_e_mix
    real:: mp_kin, mp_e_tgt, mp_Tlo, mp_Thi, mp_fa, mp_fb, mp_f, mp_df, mp_aT, mp_bT, mp_Tn, mp_fn, mp_Tguess, mp_atol
    integer  :: mp_it, mp_itmax, mp_s
    logical::mp_bracket


    ! copy/normalize X
    mp_mol_x = zero





  mp_mol_x(1:NOF_SPECIES)=MP_MOLE_FRACTION(1:nof_SPECIES)


    call normalize_X(mp_mol_x)

    U=LEFTV(2)/LEFTV(1)
    V=LEFTV(3)/LEFTV(1)

    IF (DIMENSIONA.EQ.3)THEN
    W=LEFTV(4)/LEFTV(1)
    ELSE
    W=ZERO
    END IF

    mp_itmax = 50
    mp_kin   = 0.5D0*(u*u + v*v + w*w)
    mp_e_tgt = (LEFTV(NOF_VARIABLES-NOF_SPECIES-1)/LEFTV(1)) - MP_kin        ! target sensible internal energy [J/kg]

    ! temperature bracket honoring species ranges
    mp_Tlo = 150.0D0
    mp_Thi = 20000.0D0
    do s=1,NOF_SPECIES
      if (mp_mol_x(s) <= zero) cycle
      mp_Tlo = max(mp_Tlo, MP_Tlow_in(s))
      mp_Thi = min(mp_Thi, MP_Thigh_in(s))
    end do
    if (mp_Tlo >= mp_Thi) then
      MP_info = -3; MP_TEMP = 300.0d0; return
    end if

    ! initial guess from cv around 300 K
    call mixture_props(300.0D0, mp_mol_x, mp_Rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
    mp_Tguess = 300.0d0 + (MP_e_tgt - MP_e_mix)/max(1.0e-8, MP_cv_mix)
    mp_Temp = min(max(mp_Tguess, mp_Tlo), mp_Thi)

    ! establish/verify bracket
    call mixture_props(mp_Tlo, mp_mol_x, mp_Rmix, mp_cp_mix, mp_cv_mix, mp_e_mix); mp_fa = mp_e_mix - mp_e_tgt
    call mixture_props(mp_Thi, mp_mol_x, mp_Rmix, mp_cp_mix, mp_cv_mix, mp_e_mix); mp_fb = mp_e_mix - mp_e_tgt
    mp_bracket = (mp_fa*mp_fb <= 0.0)
    mp_aT = mp_Tlo; mp_bT = mp_Thi

    do mp_it = 1, mp_itmax
      call mixture_props(mp_Temp, mp_mol_x, mp_Rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
      mp_f  = mp_e_mix - mp_e_tgt
      mp_df = max(mp_cv_mix, 1.0e-20)

      mp_atol = 1.0e-6*max(1.0d0,abs(mp_e_tgt)) + 1.0e-3
      if (abs(mp_f) <= mp_atol) then
        mp_info = 0; return
      end if

      ! safeguarded Newton step
      mp_Tn = mp_Temp - mp_f/mp_df
      if ((.not. mp_bracket) .or. (mp_Tn <= mp_aT) .or. (mp_Tn >= mp_bT)) mp_Tn = 0.5d0*(mp_aT + mp_bT)

      call mixture_props(mp_Tn, mp_mol_x, mp_Rmix, mp_cp_mix, mp_cv_mix, mp_e_mix)
      mp_fn = mp_e_mix - mp_e_tgt

      if (mp_fa*mp_fn <= 0.0) then
        mp_bT = mp_Tn; mp_fb = mp_fn; mp_bracket = .true.
      else
        mp_aT = mp_Tn; mp_fa = mp_fn
      end if
      mp_Temp = mp_Tn

      if (abs(mp_bT - mp_aT) <= 1.0e-8*max(1.0,MP_TEMP)) then
        mp_info = 1; return
      end if
    end do

    mp_info = 2  ! iteration limit
  end subroutine GET_TEMP_JANAF
























 
 
 FUNCTION FLUXEVAL2D(LEFTV)
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::FLUXEVAL2D
REAL,DIMENSION(1:nof_Variables),INTENT(IN)::LEFTV
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI
R=LEFTV(1)
U=LEFTV(2)
V=LEFTV(3)
P=LEFTV(4)
GM=GAMMA
!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)

FLUXEVAL2D(1)=R*U
FLUXEVAL2D(2)=(R*(U**2))+P
FLUXEVAL2D(3)=R*U*V
FLUXEVAL2D(4)=U*(E+P)

END FUNCTION

FUNCTION FLUXEVAL3D(LEFTV)
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::FLUXEVAL3D
REAL,DIMENSION(1:nof_Variables),INTENT(IN)::LEFTV
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI
R=LEFTV(1)
U=LEFTV(2)
V=LEFTV(3)
w=LEFTV(4)
P=LEFTV(5)
GM=GAMMA
!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)

FLUXEVAL3D(1)=R*U
FLUXEVAL3D(2)=(R*(U**2))+P
FLUXEVAL3D(3)=R*U*V
FLUXEVAL3D(4)=R*U*W
FLUXEVAL3D(5)=U*(E+P)

END FUNCTION
 
 
SUBROUTINE CONS2PRIM(N,leftv,MP_PINFl,gammal)
!> @brief
!> This subroutine transforms one vector of conservative variables to primitive variables
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:NOF_VARIABLES)::TEMPS
real,dimension(1:nof_Variables),INTENT(INOUT)::leftv
real,INTENT(INOUT)::MP_PINFL,gammal
REAL::OODENSITY,MP_DENSITY,MP_STIFF,skinx
REAL::P_SAT,P_TOL, RHO_G,RHO_L, SS_G, SS_L, PP, P_GL, VOID_FRAC,p_temp,etr
REAL,DIMENSION(1:NOF_SPECIES)::MP_AR,MP_IE,mp_vft
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:NOF_SPECIES)::RG_VFTEMP
REAL,DIMENSION(1:NOF_SPECIES)::RG_CVS,rg_cps
REAL,DIMENSION(1:NOF_SPECIES)::RG_TVSL,RG_EV
REAL::RG_VE,RG_TR,RG_CHEM,RG_DENSITY,RG_EV_TOTAL,rg_rmix,rg_kin,rg_cvt,rg_cpt
integer::rg_iter,rg_maxiter
real::rgf_vib,rgdf_vib,rg_tol,tve,TTR,rg_f,rg_df,sum1,sum2,sum3,U,V,W
REAL:: rhoE, rhoEv
REAL,dimension(1:NOF_SPECIES):: rho_i,y
REAL:: KE, etot, echem, rmix,evib
REAL:: Cv_mix, Cp_mix,rho
INTEGER :: i, idxE, idxEv






rg_maxiter=20

if (nof_Variables.gt.1)then     !only then

P_SAT =2000
P_TOL =10E-5

      IF (MULTISPECIES.EQ.1) then   !multispecies


            if (mp_modelc.eq.0)then
                  sum1=zero;
                  sum2=zero;
                  sum3=zero
                  do rg_i=1,nof_species
                      sum1=sum1+LEFTV(nof_Variables+rg_i)
                  end do
                  MP_DENSITY=sum1

                  do rg_i=1,nof_SPECIES-1
                  sum2=sum2+LEFTV(dimensiona+2+nof_species+rg_i)
                  MP_AR(rg_i)=LEFTV(dimensiona+2+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
                  mp_vft(rg_i)=LEFTV(dimensiona+2+nof_species+rg_i)
                  end do

                  MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
                  mp_vft(nof_SPECIES)=(1.0d0-sum2)
                  sum3=zero
                  do rg_i=1,nof_SPECIES
                    sum3=sum3+MP_AR(rg_i)
                    end do
                    GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
                    OODENSITY=1.0D0/MP_DENSITY

                  TEMPS(1)=MP_DENSITY
                  TEMPS(2)=LEFTV(2)*OODENSITY; u=temps(2)
                  TEMPS(3)=LEFTV(3)*OODENSITY;v=temps(3)
                  w=zero
                    if (dimensiona.eq.3)then
                    TEMPS(4)=LEFTV(4)*OODENSITY;w=TEMPS(4)
                    end if
                  sum3=zero
                  do rg_i=1,nof_SPECIES
                    sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
                    end do

                  MP_STIFF=sum3*(GAMMAL-1.0D0)
                  sum2=zero
                  do rg_i=1,nof_SPECIES
                    sum2=sum2+(mp_vft(rg_i)*MP_PINF(rg_i))
                    end do

                    MP_PINFL=sum2

                    skinx=(u*u)+(v*v)+(w*w)

                  TEMPS(dimensiona+2)=(((GAMMAL-1.0D0))*((LEFTV(dimensiona+2))-OO2*TEMPS(1)*skinx))-MP_STIFF
                  TEMPS(dimensiona+3:NOF_VARIABLES)=LEFTV(dimensiona+3:NOF_VARIABLES)
                  LEFTV(1:nof_Variables)=TEMPS(1:nof_Variables)
            end if  !mp modelc



      ELSE

                        IF (REALGAS.EQ.1)then

                                TEMPS(:)=0.0d0

                                !first get total density-correct

!
!                                 DO RG_I = 1, nof_species
!                                   IF (leftv(dimensiona+3+RG_I).LT.0.0D0)THEN
!                                       leftv(dimensiona+3+RG_I)=0.0D0
!                                   END IF
!                                 END DO







                                DO RG_I=1,nof_species
                                TEMPS(1)=temps(1)+leftv(dimensiona+3+RG_I)
                                END DO






                                RHO=temps(1)

                                U=LEFTV(2)/rho
                                v=LEFTV(3)/rho
                                IF (DIMENSIONA.EQ.3)THEN
                                w=LEFTV(4)/rho
                                Else
                                w=zero
                                end if

                               idxE  = DIMENSIONA + 2          ! index of ρE
                                idxEv = DIMENSIONA + 3          ! index of ρe_vib

                                rhoE  = LEFTV(idxE)
                                rhoEv = LEFTV(idxEv)

                                DO i = 1, NOF_SPECIES
                                  rho_i(i) = LEFTV(idxEv + i)   ! species densities ρ_i
                                END DO


                                ! -------------------------------
                                ! 2. Mass fractions
                                ! -------------------------------
                                DO i = 1, NOF_SPECIES
                                  Y(i) = rho_i(i) / RHO
                                END DO

                                ! -------------------------------
                                ! 3. Energies
                                ! -------------------------------
                                KE   = 0.5D0 * (U*U + V*V + W*W)
                                etot = rhoE / RHO
                                EVIB = rhoEv / RHO

                                ! Chemical energy: e_chem = - Σ Y_i * h°_i / M_i  [J/kg]
                                  echem = 0.0D0
                                  DO i = 1, NOF_SPECIES
                                    IF (RG_HZERO(i) > 0.0D0) THEN
                                      echem = echem - Y(i) * (RG_HZERO(i) / RG_MOLM(i))
                                    END IF
                                  END DO

                                  ! Translational–rotational internal energy
                                  etr = etot - EVIB - echem - KE
                                  IF (etr < 0.0D0) etr = 1.0D-12   ! safety

                                  ! -------------------------------
                                  ! 4. Mixture gas constant Rmix
                                  ! -------------------------------
                                  rmix = 0.0D0
                                  DO i = 1, NOF_SPECIES
                                    rmix = rmix + Y(i) / RG_MOLM(i)
                                  END DO
                                  rmix = RGS_Ru * rmix    ! [J/kg/K]

                                  ! -------------------------------
                                  ! 5. Mixture Cv, Cp, gamma
                                  !    Here Cv_i = 5/2 R for i<=3 (diatomic),
                                  !               3/2 R for i>3  (monatomic)
                                  ! -------------------------------
                                  Cv_mix = 0.0D0
                                  DO i = 1, NOF_SPECIES
                                    IF (i <= 3) THEN
                                      Cv_mix = Cv_mix + Y(i) * (5.0D0/2.0D0) * (RGS_Ru / RG_MOLM(i))
                                    ELSE
                                      Cv_mix = Cv_mix + Y(i) * (3.0D0/2.0D0) * (RGS_Ru / RG_MOLM(i))
                                    END IF
                                  END DO

                                  Cp_mix = Cv_mix + rmix
                                  GAMMAl  = Cp_mix / Cv_mix

                                  ! -------------------------------
                                  ! 6. Translational temperature Ttr and pressure
                                  !     etr = Cv_mix * Ttr  → Ttr = etr/Cv_mix
                                  !     P   = ρ Rmix Ttr
                                  ! -------------------------------
                                  IF (Cv_mix > 1.0D-20) THEN
                                    TEMPS(DIMENSIONA+2) = RHO * rmix * (etr / Cv_mix)
                                  ELSE
                                    TEMPS(DIMENSIONA+2) = 0.0D0
                                  END IF

                                mp_pinfl=zero


                                LEFTV(1)=RHO
                                LEFTV(2)=U
                                LEFTV(3)=V
                                IF (DIMENSIONA.EQ.3)THEN
                                LEFTV(4)=W
                                END IF
                                LEFTV(DIMENSIONA+2)=TEMPS(DIMENSIONA+2) !pressure
                                LEFTV(DIMENSIONA+3)=EVIB
                                DO i = 1, NOF_SPECIES
                                  LEFTV(idxEv + i)=y(i)
                                END DO


                                else      !NO REAL GAS


                                OODENSITY=1.0D0/LEFTV(1)

                                TEMPS(1)=LEFTV(1)


                                TEMPS(2)=LEFTV(2)*OODENSITY; u=temps(2)
                                TEMPS(3)=LEFTV(3)*OODENSITY; V=temps(3)
                                W=ZERO
                                IF (DIMENSIONA.EQ.3)THEN
                                TEMPS(4)=LEFTV(4)*OODENSITY; W=temps(4)
                                END IF
                                SKINX=(U**2 + V**2+ W**2)

                                TEMPS(DIMENSIONA+2)=((GAMMA-1.0D0))*((LEFTV(DIMENSIONA+2))-OO2*LEFTV(1)*SKINX)

                                LEFTV(1:nof_Variables)=TEMPS(1:nof_Variables)
                                END IF                    !

        END IF              !MULTISPECIES IF

END IF


END SUBROUTINE CONS2PRIM




SUBROUTINE CONS2DIV(N,leftv,MP_PINFl,gammal)
!> @brief
!> This subroutine transforms one vector of conservative variables to the variables required for the diffusion (r,u,v,w,T,y1,y2,y3..) or (r,u,v,w,Ttr,tve,y1,y2,etc)
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:NOF_VARIABLES)::TEMPS
real,dimension(1:nof_Variables),INTENT(INOUT)::leftv
real,INTENT(INOUT)::MP_PINFL,gammal
REAL::OODENSITY,MP_DENSITY,MP_STIFF,mp_temp
REAL::P_SAT,P_TOL, RHO_G,RHO_L, SS_G, SS_L, PP, P_GL, VOID_FRAC,p_temp,etr
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE,mp_vft
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:NOF_SPECIES)::RG_VFTEMP
REAL,DIMENSION(1:NOF_SPECIES)::RG_CVS
REAL,DIMENSION(1:NOF_SPECIES)::RG_TVSL,RG_EV
REAL::RG_VE,RG_TR,RG_CHEM,RG_DENSITY,RG_EV_TOTAL,rg_rmix,rg_kin
integer::rg_iter,rg_maxiter
real::rgf_vib,rgdf_vib,rg_tol,tve,TTR,rg_f,rg_df,sum1,sum2,sum3,Cv_i
REAL:: rhoE, rhoev, KE, etot, ev, u,v,w,Cv_mix
REAL:: g, f, df, theta, ei, dEi
REAL, DIMENSION(1:Nof_species) :: rho_i,  Cv_s
REAL :: echem, sumY
REAL:: RG_Temp, rg_T_old,denom,tmpExp
INTEGER :: i, iter
REAL, PARAMETER :: rgtol = 1.0d-10
REAL, PARAMETER :: rg_T_lo = 50.0d0, rg_T_hi = 20000.0d0
REAL:: RHO
REAL:: RG_TTR2, rg_TV
REAL,dimension(1:nof_species):: Y







if (nof_Variables.gt.1)then

IF (DIMENSIONA.EQ.3)THEN

P_SAT =2000
P_TOL =10E-5

      IF (MULTISPECIES.EQ.1) then
            if (mp_modelc.eq.0)then
            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+LEFTV(nof_Variables+rg_i)
            end do
            MP_DENSITY=sum1

            do rg_i=1,nof_SPECIES-1
            sum2=sum2+LEFTV(5+nof_species+rg_i)
            MP_AR(rg_i)=LEFTV(5+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
            mp_vft(rg_i)=LEFTV(5+nof_species+rg_i)
            end do

            MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
            mp_vft(nof_SPECIES)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_SPECIES
              sum3=sum3+MP_AR(rg_i)
              end do
              GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
              OODENSITY=1.0D0/MP_DENSITY

            TEMPS(1)=MP_DENSITY
            TEMPS(2)=LEFTV(2)*OODENSITY
            TEMPS(3)=LEFTV(3)*OODENSITY
            TEMPS(4)=LEFTV(4)*OODENSITY
            sum3=zero
            do rg_i=1,nof_SPECIES
              sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
              end do

            MP_STIFF=sum3*(GAMMAL-1.0D0)
             sum2=zero
            do rg_i=1,nof_SPECIES
              sum2=sum2+(mp_vft(rg_i)*MP_PINF(rg_i))
              end do

              MP_PINFL=sum2

            TEMPS(5)=(((GAMMAL-1.0D0))*((LEFTV(5))-OO2*TEMPS(1)*(((TEMPS(2))**2)+((TEMPS(3))**2)+((TEMPS(4))**2))))-MP_STIFF

            TEMPS(6:NOF_VARIABLES)=LEFTV(6:NOF_VARIABLES)



            call MULTISPECIES_temp(LEFTV,MP_Temp)


            TEMPS(5)=MP_Temp


!             IF(CAVITATION.EQ.1)THEN
!             RHO_G = LEFTV(6)/LEFTV(8)
!             RHO_L = LEFTV(7)/LEFTV(8)
!             SS_G = sqrt(GAMMA_IN(1)*(TEMPS(5)+MP_PINF(1))/RHO_G)
!             SS_L = sqrt(GAMMA_IN(2)*(TEMPS(5)+MP_PINF(2))/RHO_L)
!
!             P_GL=RHO_G*SS_G*SS_G*RHO_L*SS_L*SS_L*(RHO_G-RHO_l)/((RHO_G*RHO_G*SS_G*SS_G)-(RHO_l*RHO_l*SS_l*SS_l))
!             VOID_FRAC=(RHO_G*SS_G*SS_G*RHO_L*SS_L*SS_L*(RHO_L+(LEFTV(8)*(RHO_G-RHO_l))))/(RHO_l*((RHO_G*SS_G*SS_G)-LEFTV(8)*((RHO_G*SS_G*SS_G)-(RHO_L*SS_L*SS_L))))
!
!             P_tEMP=TEMPS(5)
!
!               if ((TEMPS(5).GT.P_TOL).AND.(TEMPS(5).LT.P_SAT))THEN
!               p_temp=P_SAT+P_GL*LOG(VOID_FRAC)
!               END IF
!
!               IF (TEMPS(5).LT.P_TOL)THEN
!               p_temp=P_TOL
!               end if
!             TEMPS(5)=p_temp
!
!             end if




      LEFTV(1:nof_Variables)=TEMPS(1:nof_Variables)

            end if !(mp_modelc.eq.0)then

      ELSE



      IF (REALGAS.EQ.1)then

              TEMPS(:)=0.0d0

              !first get total density-correct

!                                 DO RG_I = 1, nof_species
!                                   IF (leftv(dimensiona+3+RG_I).LT.0.0D0)THEN
!                                       leftv(dimensiona+3+RG_I)=0.0D0
!                                   END IF
!                                 END DO











              DO RG_I=1,nof_species
              TEMPS(1)=temps(1)+leftv(dimensiona+3+RG_I)
              END DO



              ! ============================================================
              ! 1. Unpack conservative variables
              ! ============================================================

              RHO  = TEMPS(1)
              U    = LEFTV(2) / RHO
              V    = LEFTV(3) / RHO
              if (dimensiona.eq.3)then
              w   =LEFTV(4) / RHO
              else
              w=zero
              end if

              rhoE = LEFTV(dimensiona+2)
              rhoev = LEFTV(dimensiona+3)

              DO i=1,Nof_species
                rho_i(i) = LEFTV(dimensiona+3+i)
              END DO

              ! ============================================================
              ! 2. Recover mass fractions
              ! ============================================================
              sumY = 0.0d0
              DO i=1,Nof_species
                Y(i) = rho_i(i) / RHO
                sumY = sumY + Y(i)
              END DO

!               ! Normalize (safety)
!               DO i=1,Nof_species
!                 Y(i) = Y(i) / sumY
!               END DO

              ! ============================================================
              ! 3. Compute kinetic energy and total specific energy
              ! ============================================================
              KE   = 0.5d0*(U*U + V*V + W*W)
              etot = rhoE / RHO

              ! Vibrational specific energy (already known from conservative vars)
              ev = rhoev / RHO

              ! ============================================================
              ! 4. Compute chemical energy
              ! ============================================================
              echem = 0.0d0
              DO i=1,Nof_species
                IF (rg_hzero(i) > 0.0d0) THEN
                  ! hzero is J/mol → convert to J/kg
                  echem = echem - Y(i) * (rg_hzero(i) / RG_MOLM(i))
                END IF
              END DO

              ! ============================================================
              ! 5. Compute remaining translational energy
              !     etr = etot - ev - echem - KE
              ! ============================================================
              etr = etot - ev - echem - KE
              IF (etr < 0.0d0) etr = 1d-12    ! safety

              ! ============================================================
              ! 6. Compute Ttr (no Newton needed: etr = Cv_mix * Ttr)
              ! ============================================================
              Cv_mix = 0.0d0
              DO i = 1, nof_species
                IF (i <= 3) THEN
                  Cv_i = (5.0d0/2.0d0)*(RGS_Ru/RG_MOLM(i))
                ELSE
                  Cv_i = (3.0d0/2.0d0)*(RGS_Ru/RG_MOLM(i))
                END IF
                Cv_mix = Cv_mix + Y(i)*Cv_i
              END DO

              IF (Cv_mix > rgs_tiny) THEN
                rg_TTR2 = etr / Cv_mix
              ELSE
                rg_TTR2 = rg_T_lo
              END IF

              rg_TTR2 = MAX(rg_T_lo, MIN(rg_T_hi, rg_TTR2))

              ! ============================================================
              ! 7. Solve for Tv using Newton iteration
              !     ev = Σ Y_i * (R/M_i) * θ/(exp(θ/Tv)-1)
              !     --> use Ttr as initial guess
              ! ============================================================

              rg_TV = rg_TTR2   ! <-- use Ttr as initial guess for Tv
              rg_TV = MAX(rg_T_lo, MIN(rg_T_hi, rg_TV))

              DO iter=1,200
                g  = 0.0d0
                df = 0.0d0

                DO i=1,3   ! only N2,O2,NO vibrate
                  IF (Y(i) < 1d-16) CYCLE
                  theta = rg_thetag(i)
                  IF (theta <= 0.0d0) CYCLE

                  ei = theta / rg_TV

                  IF (ei > 60.0d0) THEN
                    ! overflow-safe asymptotic form
                    tmpExp = EXP(-ei)
                    g  = g  + Y(i)*(RGS_Ru/RG_MOLM(i))*theta*tmpExp
                    df = df - Y(i)*(RGS_Ru/RG_MOLM(i))*theta*(ei/rg_TV)*tmpExp
                  ELSE
                    tmpExp = EXP(ei)
                    denom  = tmpExp - 1.0d0
                    g  = g  + Y(i)*(RGS_Ru/RG_MOLM(i))*(theta/denom)
                    df = df + Y(i)*(RGS_Ru/RG_MOLM(i))*theta*ei*tmpExp / (denom*denom*rg_TV)
                  END IF
                END DO

                g = g - ev        ! target equation ev(Tv) - ev_known = 0

                IF (ABS(g) < rg_tol) EXIT
                IF (ABS(df) < 1d-20) EXIT

                rg_TV = rg_TV - g/df

                ! enforce bounds
                IF (rg_TV < rg_T_lo) rg_TV = rg_T_lo
                IF (rg_TV > rg_T_hi) rg_TV = rg_T_hi
              END DO



              leftv(1)=RHO
              leftv(2)=U
              leftv(3)=V

              IF (DIMENSIONA.EQ.3)THEN
              leftv(4)=W
              END IF

              LEFTV(DIMENSIONA+2)=rg_TTR2
              LEFTV(DIMENSIONA+3)=rg_TV



              DO i=1,Nof_species
                 LEFTV(dimensiona+3+i)=Y(i)
              END DO



              else















      OODENSITY=1.0D0/LEFTV(1)

      TEMPS(1)=LEFTV(1)
      TEMPS(2)=LEFTV(2)*OODENSITY
      TEMPS(3)=LEFTV(3)*OODENSITY
      TEMPS(4)=LEFTV(4)*OODENSITY
      TEMPS(5)=((GAMMA-1.0D0))*((LEFTV(5))-OO2*LEFTV(1)*(((TEMPS(2))**2)+((TEMPS(3))**2)+((TEMPS(4))**2)))

      TEMPS(5)=  leftv(5)/(leftv(1)*R_gas)  !temperature

      LEFTV(1:nof_Variables)=TEMPS(1:nof_Variables)

      END IF

END IF




ELSE

P_SAT =2000
P_TOL =10E-5


IF ((MULTISPECIES.EQ.1)) then
            if (mp_modelc.eq.0)then
            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+LEFTV(nof_Variables+rg_i)
            end do
            MP_DENSITY=sum1

            do rg_i=1,nof_SPECIES-1
            sum2=sum2+LEFTV(4+nof_species+rg_i)
            MP_AR(rg_i)=LEFTV(4+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
            mp_vft(rg_i)=LEFTV(4+nof_species+rg_i)
            end do

             MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
            mp_vft(nof_SPECIES)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_SPECIES
              sum3=sum3+MP_AR(rg_i)
              end do
              GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
              OODENSITY=1.0D0/MP_DENSITY

              TEMPS(1)=MP_DENSITY
              TEMPS(2)=LEFTV(2)*OODENSITY
              TEMPS(3)=LEFTV(3)*OODENSITY

            sum3=zero
            do rg_i=1,nof_SPECIES
              sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
              end do

            MP_STIFF=sum3*(GAMMAL-1.0D0)

             sum2=zero
            do rg_i=1,nof_SPECIES
              sum2=sum2+(mp_vft(rg_i)*MP_PINF(rg_i))
              end do

              MP_PINFL=sum2

            TEMPS(4)=(((GAMMAL-1.0D0))*((LEFTV(4))-OO2*TEMPS(1)*(((TEMPS(2))**2)+((TEMPS(3))**2))))-MP_STIFF

            TEMPS(5:NOF_VARIABLES)=LEFTV(5:NOF_VARIABLES)

            call MULTISPECIES_temp(LEFTV,MP_Temp)


            TEMPS(4)=MP_Temp








 LEFTV(1:nof_Variables)=TEMPS(1:nof_Variables)

           end  if !(mp_modelc.eq.0)then

ELSE


IF (REALGAS.EQ.1)then

TEMPS(:)=0.0d0

              !first get total density-correct

!                DO RG_I = 1, nof_species
!                                   IF (leftv(dimensiona+3+RG_I).LT.0.0D0)THEN
!                                       leftv(dimensiona+3+RG_I)=0.0D0
!                                   END IF
!                                 END DO






              DO RG_I=1,nof_species
              TEMPS(1)=temps(1)+leftv(dimensiona+3+RG_I)
              END DO



              ! ============================================================
              ! 1. Unpack conservative variables
              ! ============================================================

              RHO  = TEMPS(1)
              U    = LEFTV(2) / RHO
              V    = LEFTV(3) / RHO
              if (dimensiona.eq.3)then
              w   =LEFTV(4) / RHO
              else
              w=zero
              end if

              rhoE = LEFTV(dimensiona+2)
              rhoev = LEFTV(dimensiona+3)

              DO i=1,Nof_species
                rho_i(i) = LEFTV(dimensiona+3+i)
              END DO

              ! ============================================================
              ! 2. Recover mass fractions
              ! ============================================================
              sumY = 0.0d0
              DO i=1,Nof_species
                Y(i) = rho_i(i) / RHO
                sumY = sumY + Y(i)
              END DO

!               ! Normalize (safety)
!               DO i=1,Nof_species
!                 Y(i) = Y(i) / sumY
!               END DO

              ! ============================================================
              ! 3. Compute kinetic energy and total specific energy
              ! ============================================================
              KE   = 0.5d0*(U*U + V*V + W*W)
              etot = rhoE / RHO

              ! Vibrational specific energy (already known from conservative vars)
              ev = rhoev / RHO

              ! ============================================================
              ! 4. Compute chemical energy
              ! ============================================================
              echem = 0.0d0
              DO i=1,Nof_species
                IF (rg_hzero(i) > 0.0d0) THEN
                  ! hzero is J/mol → convert to J/kg
                  echem = echem - Y(i) * (rg_hzero(i) / RG_MOLM(i))
                END IF
              END DO

              ! ============================================================
              ! 5. Compute remaining translational energy
              !     etr = etot - ev - echem - KE
              ! ============================================================
              etr = etot - ev - echem - KE
              IF (etr < 0.0d0) etr = 1d-12    ! safety

              ! ============================================================
              ! 6. Compute Ttr (no Newton needed: etr = Cv_mix * Ttr)
              ! ============================================================
              Cv_mix = 0.0d0
              DO i = 1, nof_species
                IF (i <= 3) THEN
                  Cv_i = (5.0d0/2.0d0)*(RGS_Ru/RG_MOLM(i))
                ELSE
                  Cv_i = (3.0d0/2.0d0)*(RGS_Ru/RG_MOLM(i))
                END IF
                Cv_mix = Cv_mix + Y(i)*Cv_i
              END DO

              IF (Cv_mix > rgs_tiny) THEN
                rg_TTR2 = etr / Cv_mix
              ELSE
                rg_TTR2 = rg_T_lo
              END IF

              rg_TTR2 = MAX(rg_T_lo, MIN(rg_T_hi, rg_TTR2))

              ! ============================================================
              ! 7. Solve for Tv using Newton iteration
              !     ev = Σ Y_i * (R/M_i) * θ/(exp(θ/Tv)-1)
              !     --> use Ttr as initial guess
              ! ============================================================

              rg_TV = rg_TTR2   ! <-- use Ttr as initial guess for Tv
              rg_TV = MAX(rg_T_lo, MIN(rg_T_hi, rg_TV))

              DO iter=1,200
                g  = 0.0d0
                df = 0.0d0

                DO i=1,3   ! only N2,O2,NO vibrate
                  IF (Y(i) < 1d-16) CYCLE
                  theta = rg_thetag(i)
                  IF (theta <= 0.0d0) CYCLE

                  ei = theta / rg_TV

                  IF (ei > 60.0d0) THEN
                    ! overflow-safe asymptotic form
                    tmpExp = EXP(-ei)
                    g  = g  + Y(i)*(RGS_Ru/RG_MOLM(i))*theta*tmpExp
                    df = df - Y(i)*(RGS_Ru/RG_MOLM(i))*theta*(ei/rg_TV)*tmpExp
                  ELSE
                    tmpExp = EXP(ei)
                    denom  = tmpExp - 1.0d0
                    g  = g  + Y(i)*(RGS_Ru/RG_MOLM(i))*(theta/denom)
                    df = df + Y(i)*(RGS_Ru/RG_MOLM(i))*theta*ei*tmpExp / (denom*denom*rg_TV)
                  END IF
                END DO

                g = g - ev        ! target equation ev(Tv) - ev_known = 0

                IF (ABS(g) < rg_tol) EXIT
                IF (ABS(df) < 1d-20) EXIT

                rg_TV = rg_TV - g/df

                ! enforce bounds
                IF (rg_TV < rg_T_lo) rg_TV = rg_T_lo
                IF (rg_TV > rg_T_hi) rg_TV = rg_T_hi
              END DO



              leftv(1)=RHO
              leftv(2)=U
              leftv(3)=V

              IF (DIMENSIONA.EQ.3)THEN
              leftv(4)=W
              END IF

              LEFTV(DIMENSIONA+2)=rg_TTR2
              LEFTV(DIMENSIONA+3)=rg_TV



              DO i=1,Nof_species
                 LEFTV(dimensiona+3+i)=Y(i)
              END DO









else















OODENSITY=1.0D0/LEFTV(1)

TEMPS(1)=LEFTV(1)
TEMPS(2)=LEFTV(2)*OODENSITY
TEMPS(3)=LEFTV(3)*OODENSITY
TEMPS(4)=((GAMMA-1.0D0))*((LEFTV(4))-OO2*LEFTV(1)*(((TEMPS(2))**2)+((TEMPS(3))**2)))

TEMPS(4)=  leftv(4)/(leftv(1)*R_gas)  !temperature

LEFTV(1:nof_Variables)=TEMPS(1:nof_Variables)



END IF





END IF




end if


end if


END SUBROUTINE CONS2DIV














SUBROUTINE CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
!> @brief
!> This subroutine transforms two vector of conservative variables to primitive variables
IMPLICIT NONE
INTEGER,INTENT(IN)::N
real,dimension(1:nof_Variables),INTENT(INOUT)::leftv
real,INTENT(INOUT)::MP_PINFL,gammal
real,dimension(1:nof_Variables),INTENT(INOUT)::RIGHTv
real,INTENT(INOUT)::MP_PINFR,gammaR

CALL cons2prim(N,leftv,MP_PINFl,gammal)
CALL cons2prim(N,RIGHTV,MP_PINFR,gammaR)


END SUBROUTINE CONS2PRIM2


SUBROUTINE LMACHT(N,LEFTV,RIGHTV)
!> @brief
!> This subroutine applies the low-Mach number correction to two vectors of conserved variables
IMPLICIT NONE
INTEGER,INTENT(IN)::N
real,dimension(1:nof_Variables),INTENT(INOUT)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables),INTENT(INOUT)::RIGHTv
real::MP_PINFR,gammaR
REAL::Q2L,Q2R,UUL,UUR,VVL,VVR,WWR,WWL,RHOL,RHOR,ETAL,ETAR,DUU,DVV,DWW
REAL::MACH2,MACH,CMA,DUS,DVS,DWS,DIFF,C1o2,SSL,SSR,ppl,ppr,eel,eer,TOLE,MLM

TOLE=ZERO

 EEL=LEFTV(5)
 EER=RIGHTV(5)

 CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)


      C1o2=0.50D0
      RHOL=LEFTV(1)
      UUL=LEFTV(2)
      VVL=LEFTV(3)
      WWL=LEFTV(4)
       PPL=LEFTV(5)
     RHOR=RIGHTV(1)
      UUR=RIGHTV(2)
      VVR=RIGHTV(3)
      WWR=RIGHTV(4)
       PPR=RIGHTV(5)
		
      Q2L=(UUL*uul)+(vvl*vvl)+(wwl*wwl)
      Q2R=(uur*uur)+(vvr*vvr)+(wwr*wwr)

!	SSL=((GAMMA*PPL)/(RHOL)); SSR=((GAMMA*PPR)/(RHOR))
                                  
!	MLM=((ABS(PPR-PPL))/(0.5*RRES*((GAMMA*PRES/RRES))))




      if ((multispecies.eq.1).OR.(REALGAS.EQ.1))then
		 SSL=SQRT((LEFTV(5)+MP_PINFL)*GAMMAl/LEFTV(1))
		  Ssr=SQRT((rightV(5)+MP_PINFr)*GAMMAr/rightV(1))
		else

      SSL=((GAMMA*PPL)/(RHOL))
      SSR=((GAMMA*PPR)/(RHOR))
      end if


      

      CMA=1.0D0

      DUU=UUR-UUL
      DVV=VVR-VVL
      DWW=WWR-WWL

!       IF(LMACH.EQ.1) THEN !Standard proportional to du^2
         MACH2=MAX(Q2L/SSL,Q2R/SSR)
         MACH=sqrt(MACH2)
         MACH=MIN(CMA*MACH,1.0D0)
!       END IF

      DUS=UUR+UUL
      DVS=VVR+VVL
      DWS=WWR+WWL

      DUU=MACH*DUU
      DVV=MACH*DVV
      DWW=MACH*DWW


      DIFF=C1O2*DUU

      UUL=(DUS*C1O2-DIFF)
      UUR=(DUS*C1O2+DIFF)

      if (lmach_style.eq.1)then
       DIFF=C1O2*DVV

      VVL=(DVS*C1O2-DIFF)
      VVR=(DVS*C1O2+DIFF)

      DIFF=C1O2*DWW

      WWL=(DWS*C1O2-DIFF)
      WWR=(DWS*C1O2+DIFF)
      end if
	
    	
       LEFTV(1)=RHOL
       LEFTV(2)=UUL*RHOL
       LEFTV(3)=VVL*RHOL
       LEFTV(4)=WWL*RHOL
       LEFTV(5)=EEL
	
      RIGHTV(1)=RHOR 	
       RIGHTV(2)=UUR*RHOR
       RIGHTV(3)=VVR*RHOR
       RIGHTV(4)=WWR*RHOR
       RIGHTV(5)=EER




END SUBROUTINE LMACHT



SUBROUTINE LMACHT2D(N,LEFTV,RIGHTV)
!> @brief
!> This subroutine applies the low-Mach number correction to two vectors of conserved variables 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
real,dimension(1:nof_Variables),INTENT(INOUT)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables),INTENT(INOUT)::RIGHTv
real::MP_PINFR,gammaR
REAL::Q2L,Q2R,UUL,UUR,VVL,VVR,WWR,WWL,RHOL,RHOR,ETAL,ETAR,DUU,DVV,DWW
REAL::MACH2,MACH,CMA,DUS,DVS,DWS,DIFF,C1o2,SSL,SSR,ppl,ppr,eel,eer,TOLE,MLM

TOLE=tolsmall



 EEL=LEFTV(4)
 EER=RIGHTV(4)


 
 CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)



      C1o2=0.5D0
      RHOL=LEFTV(1)
      UUL=LEFTV(2)
      VVL=LEFTV(3)
       PPL=LEFTV(4)
       
     RHOR=RIGHTV(1)
      UUR=RIGHTV(2)
      VVR=RIGHTV(3)
       PPR=RIGHTV(4)
		
       Q2L=(UUL*uul)+(vvl*vvl)
       Q2R=(uur*uur)+(vvr*vvr)
      
      
      
    if ((multispecies.eq.1).OR.(REALGAS.EQ.1))then
		 SSL=SQRT((LEFTV(4)+MP_PINFL)*GAMMAl/LEFTV(1))
		  Ssr=SQRT((rightV(4)+MP_PINFr)*GAMMAr/rightV(1))
		else

      SSL=((GAMMA*PPL)/(RHOL))
      SSR=((GAMMA*PPR)/(RHOR))
      end if

      CMA=1.0D0

      DUU=UUR-UUL
      DVV=VVR-VVL
      
      

!       IF(LMACH.EQ.1) THEN !Standard proportional to du^2
         MACH2=MAX(Q2L/SSL,Q2R/SSR)
         MACH=sqrt(MACH2)
         MACH=MIN(CMA*MACH,1.0D0)
!       END IF

      DUS=UUR+UUL
      DVS=VVR+VVL
     

       !UL+ZUL+UR-ZUR=(UL+UR)-Z(UR-UL))
      DUU=MACH*DUU
      DVV=MACH*DVV
      


      DIFF=C1O2*DUU
	!if (uul*uur.gt.tole)then
	
      UUL=(DUS*C1O2)-DIFF
      UUR=(DUS*C1O2)+DIFF
      
if (lmach_style.eq.1)then
      DIFF=C1O2*DVV
      VVL=(DVS*C1O2-DIFF)
      VVR=(DVS*C1O2+DIFF)
end if
	
       LEFTV(1)=RHOL
       LEFTV(2)=UUL*RHOL
       LEFTV(3)=VVL*RHOL
       LEFTV(4)=EEL
	
       RIGHTV(1)=RHOR 	
       RIGHTV(2)=UUR*RHOR
       RIGHTV(3)=VVR*RHOR
       RIGHTV(4)=EER




END SUBROUTINE LMACHT2D































SUBROUTINE PRIM2CONS(N,leftv)
!> @brief
! !> This subroutine transforms one vector of primitive variables to conservative variables
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:nof_Variables)::TEMPS
REAL::OODENSITY,skin1,ie1,MP_DENSITY,mp_stiff,RG_VE,RG_TR,RG_CHEM
real,dimension(1:nof_Variables),INTENT(INOUT)::leftv
real::MP_PINFL,gammal,p
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE,mp_vft
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:nof_species)::RG_VFTEMP
REAL,DIMENSION(1:nof_species)::RG_CVS,rg_theta
REAL,DIMENSION(1:nof_species)::RG_TVSL,RG_EV,y
REAL::RG_DENSITY,RG_EV_TOTAL,rg_rmix,rg_kin
integer::rg_iter,rg_maxiter
real::rgf_vib,rgdf_vib,rg_tol,tve,TTR,rg_f,rg_df,sum1,sum2,sum3,u,v,w
REAL:: rmix
REAL:: Cv_mix, e_tr, e_chem, KE, e_tot,evib,rho
INTEGER :: i, idxE, idxEv

rg_maxiter=20

if (nof_Variables.gt.1)then

IF (MULTISPECIES.EQ.1) then

            if (mp_modelc.eq.0)then
            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+LEFTV(nof_Variables+rg_i)
            end do
            MP_DENSITY=sum1

            do rg_i=1,nof_SPECIES-1
            sum2=sum2+LEFTV(dimensiona+2+nof_species+rg_i)
            MP_AR(rg_i)=LEFTV(dimensiona+2+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
            mp_vft(rg_i)=LEFTV(dimensiona+2+nof_species+rg_i)
            end do
            MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
            mp_vft(nof_SPECIES)=(1.0d0-sum2)

             sum3=zero
            do rg_i=1,nof_SPECIES
              sum3=sum3+MP_AR(rg_i)
              end do
              GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
 

 
          TEMPS(1)=MP_DENSITY
          TEMPS(2)=LEFTV(2)*TEMPS(1);u=LEFTV(2)
          TEMPS(3)=LEFTV(3)*TEMPS(1);v=LEFTV(3)
          w=zero
          if (dimensiona.eq.3)then
                    TEMPS(4)=LEFTV(4)*TEMPS(1);w=LEFTV(4)
                    end if
          skin1=(oo2)*((u*u)+(v*v)+(w*w))

            sum3=zero
            do rg_i=1,nof_SPECIES
              sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
              end do

            MP_STIFF=sum3*(GAMMAL-1.0D0)


! MP_STIFF=((LEFTV(8)*(GAMMA_IN(1)/(GAMMA_IN(1)-1.0D0))*MP_PINF(1))+((1.0D0-LEFTV(8))*(GAMMA_IN(2)/(GAMMA_IN(2)-1.0D0))*MP_PINF(2)))*(GAMMAL-1.0D0)

      ie1=((leftv(dimensiona+2)+mp_stiff)/((GAMMAL-1.0D0)*TEMPS(1)))
      TEMPS(dimensiona+2)=TEMPS(1)*(ie1+skin1)
      TEMPS(dimensiona+3:NOF_VARIABLES)=LEFTV(dimensiona+3:NOF_VARIABLES)
      LEFTV(1:nof_Variables)=TEMPS(1:nof_Variables)
 
            end if !(mp_modelc.eq.0)then
 
 else


            IF (REALGAS.EQ.1)THEN


            RHO=leftv(1)
            U=leftv(2)
            v=leftv(3)

            if (dimensiona.eq.3)then
            w=leftv(4)
            ELSE
            w=zero
            end if

            idxE  = DIMENSIONA + 2
            idxEv = DIMENSIONA + 3

            p=leftv(idxE)
            evib=leftv(idxEv)

                              DO i = 1, NOF_SPECIES
                                  Y(i) = LEFTV(idxEv + i)   ! species mass fractions
                                END DO




             DO i = 1, NOF_SPECIES
              rmix = rmix + Y(i) / RG_MOLM(i)
            END DO




            ! -------------------------------
            ! 1. Mixture gas constant Rmix
            ! P = ρ Rmix Ttr → Ttr = P / (ρ Rmix)
            ! -------------------------------
            rmix = 0.0D0
            DO i = 1, NOF_SPECIES
              rmix = rmix + Y(i) / RG_MOLM(i)
            END DO
            rmix = RGS_Ru * rmix

            ! -------------------------------
            ! 2. Mixture Cv and translational energy
            !    Cv_i = 5/2 R for i<=3 (diatomic),
            !           3/2 R for i>3  (monatomic)
            ! -------------------------------
            Cv_mix = 0.0D0
            DO i = 1, NOF_SPECIES
              IF (i <= 3) THEN
                Cv_mix = Cv_mix + Y(i) * (5.0D0/2.0D0) * (RGS_Ru / RG_MOLM(i))
              ELSE
                Cv_mix = Cv_mix + Y(i) * (3.0D0/2.0D0) * (RGS_Ru / RG_MOLM(i))
              END IF
            END DO

            ! Translational temperature Ttr from P = ρ Rmix Ttr
            ! Then e_tr = Cv_mix * Ttr
            IF (rmix > 1.0D-20 .AND. Cv_mix > 1.0D-20) THEN
              e_tr = Cv_mix * (P / (RHO * rmix))
            ELSE
              e_tr = 0.0D0
            END IF

            ! -------------------------------
            ! 3. Chemical energy: e_chem = - Σ Y_i * h°_i/M_i
            ! -------------------------------
            e_chem = 0.0D0
!             DO i = 1, NOF_SPECIES
!               IF (RG_HZERO(i) > 0.0D0) THEN
!                 e_chem = e_chem - Y(i) * (RG_HZERO(i) / RG_MOLM(i))
!               END IF
!             END DO

            ! -------------------------------
            ! 4. Kinetic and total specific energy
            ! -------------------------------
            KE    = 0.5D0 * (U*U + V*V + W*W)
            e_tot = e_tr + EVIB + e_chem + KE

            ! -------------------------------
            ! 5. Fill conservative vector
            ! -------------------------------
            LEFTV(1) = RHO
            LEFTV(2) = RHO * U
            LEFTV(3) = RHO * V
            IF (DIMENSIONA == 3) THEN
              LEFTV(4) = RHO * W
            ELSE
              LEFTV(4) = 0.0D0
            END IF



            LEFTV(idxE)  = RHO * e_tot      ! ρE_total
            LEFTV(idxEv) = RHO * EVIB       ! ρ e_vib

            DO i = 1, NOF_SPECIES
              LEFTV(idxEv + i) = RHO * Y(i) ! ρ_i
            END DO




            else



                      u=LEFTV(2)
                      v=LEFTV(3)
                      w=zero
                            if (dimensiona.eq.3)then
                                w=LEFTV(4)
                                end if
                      skin1=(oo2)*((u*u)+(v*v)+(w*w))
            ie1=((leftv(dimensiona+2))/((GAMMA-1.0D0)*leftv(1)))

            OODENSITY=1.0D0/LEFTV(1)

            TEMPS(1)=LEFTV(1)
            TEMPS(2)=LEFTV(2)*LEFTV(1)
            TEMPS(3)=LEFTV(3)*LEFTV(1)
            if (dimensiona.eq.3)then
            TEMPS(4)=LEFTV(4)*LEFTV(1)
            end if
            TEMPS(dimensiona+2)=leftv(1)*(ie1+skin1)

            LEFTV(1:nof_Variables)=TEMPS(1:nof_Variables)

            end if

end if

end if






END SUBROUTINE PRIM2CONS







SUBROUTINE PRIM2CONS2(N,LEFTV,RIGHTV)
!> @brief
!> This subroutine transforms two vectors of primitive variables to conservative variables
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:nof_Variables)::TEMPS
REAL::OODENSITY,skin1,ie1,MP_DENSITY,mp_stiff
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE,mp_vft
real,dimension(1:nof_Variables),INTENT(INOUT)::leftv,RIGHTV
real::MP_PINFL,gammal,MP_PINFR,gammaR
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:nof_species)::RG_VFTEMP
REAL,DIMENSION(1:nof_species)::RG_CVS,rg_theta
REAL,DIMENSION(1:nof_species)::RG_TVSL,RG_EV
REAL::RG_DENSITY,RG_EV_TOTAL,rg_rmix,rg_kin
integer::rg_iter,rg_maxiter
real::rgf_vib,rgdf_vib,rg_tol,tve,TTR,rg_f,rg_df,sum1,sum2,sum3

if (nof_Variables.gt.1)then

call PRIM2CONS(N,leftv)
call PRIM2CONS(N,rightv)


end if


END SUBROUTINE PRIM2CONS2












FUNCTION INFLOW(INITCOND,POX,POY,POZ)
!> @brief
!> This function applies a prescribed boundary condition to  the inflow in 3D
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::INFLOW
INTEGER,INTENT(IN)::INITCOND
REAL,DIMENSION(1:dimensiona),INTENT(IN)::POX,POY,POZ
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI
REAL::XF,YF,ZF
REAL:: Theta_0,vtang, vradial,GAMMAR
REAL::MP_DENSITY,MP_STIFF,SUM1,SUM2,SUM3
REAL,DIMENSION(1:NOF_VARIABLES)::VECT_IN
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE
INTEGER::RG_I,RG_J
REAL::KHX,T1L,VHX,AMP,DVEL,rgg,tt1,khi_slope,khi_b,theeta,reeta,RG_VE,RG_TR,RG_CHEM,RG_DENSITY,RG_EV_TOTAL,rg_rmix
REAL,DIMENSION(1:NOF_SPECIES)::RG_CVS
REAL,DIMENSION(1:NOF_SPECIES)::RG_TVSL,RG_EV
real::rg_tv,rg_Ttr0,rg_Tve0

IF (multispecies.EQ.1) then



P=PRES
U=uvel
V=vvel
w=wvel


VECT_IN(1)=R
VECT_IN(2)=U
VECT_IN(3)=V
VECT_IN(4)=w
VECT_IN(5)=p

DO RG_I=1,NOF_SPECIES
VECT_IN(6+RG_I)=MP_R_IN(RG_I)*MP_A_IN(RG_I)
END DO
DO RG_I=1,NOF_SPECIES-1
VECT_IN(6+NOF_SPECIES+RG_I)=MP_A_IN(RG_I)
END DO

   CALL PRIM2CONS(N,VECT_IN)


INFLOW(1:NOF_VARIABLES)=VECT_IN(1:NOF_VARIABLES)


!
!
!
!
!
!
! MP_AR(1)=MP_A_IN(1)/(GAMMA_IN(1)-1.0D0)
! MP_AR(2)=MP_A_IN(2)/(GAMMA_IN(2)-1.0D0)
! GAMMAR=(1.0D0/(MP_AR(1)+MP_AR(2)))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
!
! GM=GAMMAR
!
! R=(MP_R_IN(1)*MP_A_IN(1))+(MP_R_IN(2)*MP_A_IN(2))
! MP_IE(1)=((P+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
! MP_IE(2)=((P+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
!
! IEn=(MP_IE(1)*MP_A_IN(1))+(MP_IE(2)*MP_A_IN(2))
! ! !KINETIC ENERGY FIRST!
! SKIN=(OO2)*((U**2)+(V**2)+(w**2))
! ! !TOTAL ENERGY
! E=(R*SKIN)+IEN
!
! !VECTOR OF CONSERVED VARIABLES NOW
! INFLOW(1)=R
! INFLOW(2)=R*U
! INFLOW(3)=R*V
! INFLOW(4)=R*w
! INFLOW(5)=E
! INFLOW(6)=MP_R_IN(1)*MP_A_IN(1)
! INFLOW(7)=MP_R_IN(2)*MP_A_IN(2)
! INFLOW(8)=MP_A_IN(1)



ELSE

R=RRES
GM=GAMMA
P=PRES
U=uvel
V=vvel
W=wvel



if (initcond.eq.10000)then
if (sqrt(((poy(1)-0.0)**2)+((poz(1)-0.5)**2)).le.0.05)then
!if (((poy(1).ge.-0.05).and.(poy(1).le.0.05)).and.((poz(1).ge.-0.05).and.(poz(1).le.0.05)))then
p=0.4127
	R=5
	u=30.0
	v=0.0d0
	w=0.0d0

end if
end if







!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)
!VECTOR OF CONSERVED VARIABLES NOW
INFLOW(1)=R
INFLOW(2)=R*U
INFLOW(3)=R*V
INFLOW(4)=R*W
INFLOW(5)=E

end if


IF (SWIRL.EQ.1)THEN

IF (POX(1).LT.-0.03)THEN
XF=POX(1)
YF=POY(1)
ZF=POZ(1)
Theta_0=atan2(ZF,YF)
Vtang=18.0375D0
Vradial=-12.63D0
U=0.0D0

V=-Vtang*sin(Theta_0)+Vradial*cos(Theta_0)
W=Vtang*cos(Theta_0)+Vradial*sin(Theta_0)



ELSE

 U=70.06D0
  V=0.0D0
  W=0.0D0


END IF

R=RRES
GM=GAMMA
P=PRES
S=SQRT((GM*P)/(R))
  
  
!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)
!VECTOR OF CONSERVED VARIABLES NOW
INFLOW(1)=R
INFLOW(2)=R*U
INFLOW(3)=R*V
INFLOW(4)=R*W
INFLOW(5)=E  
  
  
END IF


IF (REALGAS.EQ.1)THEN
U=UVEL
V=VVEL
w=wVEL
P=PRES
R=RRES



!First build mixture gas constant
rg_rmix=zero
do rg_i=1,nof_species
  rg_rmix=rg_rmix+RG_VF(rg_i)/RG_MOLM(rg_i)
end do

  rg_rmix=RGS_Ru*rg_rmix

  rg_ttr0=rg_ttr
  rg_tve0=rg_tve


! Translational-rotational internal energy
rg_tr = 0.0D0
    do rg_i = 1, nof_species
      if (rg_i <= 3) then
        RG_CVS(rg_i) = (5.0D0 / 2.0D0) * RGS_Ru / RG_MOLM(rg_i)
      else
        RG_CVS(rg_i) = (3.0D0 / 2.0D0) * RGS_Ru / RG_MOLM(rg_i)
      end if
      rg_tr = rg_tr + RG_VF(rg_i) * RG_CVS(rg_i) * rg_Ttr0
    end do


! Vibrational energy
    RG_EV_TOTAL= 0.0D0
    do rg_i = 1, 3
      RG_EV_TOTAL = RG_EV_TOTAL + RG_VF(rg_i)  * (RGS_Ru / RG_MOLM(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i)/rg_Tve0) - 1.0D0))
    end do

RG_CHEM=zero

  ! Chemical energy
 do rg_i=1,nof_species
        if (rg_hzero(RG_i).gt.1.0e-12)then
        RG_CHEM=RG_CHEM-(RG_VF(rg_i)*rg_hzero(RG_i)/RG_MOLM(rg_i))
        end if
END DO
!RG_CHEM=zero





! Kinetic energy

SKIN=(oo2)*((U**2)+(V**2)+(W**2))


inflow(1)=R
inflow(2)=R*U
inflow(3)=R*V
inflow(4)=R*w
inflow(5)=R*(RG_EV_TOTAL+RG_TR+RG_CHEM+skin)
inflow(6)=R*RG_EV_TOTAL

do rg_i=1,nof_species
inflow(6+rg_i)=r * RG_VF(rg_i)
end do



end if






END FUNCTION INFLOW




FUNCTION VECT_FUNCTION(POX,POY)
!> @brief
!> This makes a multipliciation between two vectors
IMPLICIT NONE
REAL,DIMENSION(3)::VECT_FUNCTION
REAL,DIMENSION(1:dimensiona),INTENT(IN)::POX,POY

VECT_FUNCTION(1)=(POY(2)*POX(3))-(POY(3)*POX(2))
VECT_FUNCTION(2)=(POY(3)*POX(1))-(POY(1)*POX(3))
VECT_FUNCTION(3)=(POY(1)*POX(2))-(POY(2)*POX(1))





END FUNCTION VECT_FUNCTION


FUNCTION INFLOW2d(INITCOND,POX,POY)
!> @brief
!> This function applies a prescribed boundary condition to  the inflow in 2D
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::INFLOW2d
INTEGER,INTENT(IN)::INITCOND
REAL,DIMENSION(1:2),INTENT(IN)::POX,POY
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI,PS
REAL::XF,YF,ZF,LIT_A,LIT_O
REAL:: Theta_0,vtang, vradial,GAMMAR
REAL::MP_DENSITY,MP_STIFF,SUM1,SUM2,SUM3
REAL,DIMENSION(1:NOF_VARIABLES)::VECT_IN
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE
INTEGER::RG_I,RG_J
REAL::KHX,T1L,VHX,AMP,DVEL,rgg,tt1,khi_slope,khi_b,theeta,reeta,RG_VE,RG_TR,RG_CHEM,RG_DENSITY,RG_EV_TOTAL,rg_rmix
REAL,DIMENSION(1:NOF_SPECIES)::RG_CVS
REAL,DIMENSION(1:NOF_SPECIES)::RG_TVSL,RG_EV
real::rg_tv,rg_Ttr0,rg_Tve0





IF (multispecies.EQ.1) then



P=PRES
U=uvel
V=vvel

!time variable boundary condition
if (initcond.eq.430)then
v=(179299.375638680*(t**5)) - (82455.0868677361*(t**4)) + (14299.8472891299*(t**3)) - (1281.65548492021*(t**2)) + (62.2260666329356*t) - (0.0419033282181554)
end if


! IF (INITCOND.EQ.157)THEN
! PS=35*10e6
! LIT_A=1.48*10E8
! LIT_o=1.21*10E8
! p=PRES+2.0D0*PS*EXP(-LIT_A*T)*COS((LIT_O*T)+(PI/3.0))
! UVEL=0.0
! vVEL=0.0
! end if



VECT_IN(1)=R
VECT_IN(2)=U
VECT_IN(3)=V
VECT_IN(4)=P
DO RG_I=1,NOF_SPECIES
VECT_IN(5+RG_I)=MP_R_IN(RG_I)*MP_A_IN(RG_I)
END DO
DO RG_I=1,NOF_SPECIES-1
VECT_IN(5+NOF_SPECIES+RG_I)=MP_A_IN(RG_I)
END DO

   CALL PRIM2CONS(N,VECT_IN)


INFLOW2D(1:NOF_VARIABLES)=VECT_IN(1:NOF_VARIABLES)

! !VECTOR OF CONSERVED VARIABLES NOW
! INFLOW2d(1)=R
! INFLOW2d(2)=R*U
! INFLOW2d(3)=R*V
! INFLOW2d(4)=E
! INFLOW2d(5)=MP_R_IN(1)*MP_A_IN(1)
! INFLOW2d(6)=MP_R_IN(2)*MP_A_IN(2)
! INFLOW2d(7)=MP_A_IN(1)



ELSE




R=RRES
GM=GAMMA
P=PRES
U=uvel
V=vvel

if (initcond.eq.133)then
p=195557.25
	R=p/(350.5d0*287.058d0)
	u=168.62
	v=0.0d0

end if

if (initcond.eq.10000)then
if ((poy(1).ge.-0.05).and.(poy(1).le.0.05))then
p=0.4127
	R=5
	u=30.0
	v=0.0d0

end if
end if




IF (INITCOND.EQ.790)THEN
R=(2.4D0*6**2)/((0.4*6**2)+2)
U=(6*SQRT(1.4))*(70/(2.4*36))
V=0.0D0
P=(2.8*36-0.4)/(2.4)


END IF




!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)
!VECTOR OF CONSERVED VARIABLES NOW
INFLOW2d(1)=R
INFLOW2d(2)=R*U
INFLOW2d(3)=R*V
INFLOW2d(4)=E


ENDIF

IF (REALGAS.EQ.1)THEN
U=UVEL
V=VVEL
P=PRES
R=RRES



!First build mixture gas constant
rg_rmix=zero
do rg_i=1,nof_species
  rg_rmix=rg_rmix+RG_VF(rg_i)/RG_MOLM(rg_i)
end do

  rg_rmix=RGS_Ru*rg_rmix

  rg_ttr0=rg_ttr
  rg_tve0=rg_tve


! Translational-rotational internal energy
rg_tr = 0.0D0
    do rg_i = 1, nof_species
      if (rg_i <= 3) then
        RG_CVS(rg_i) = (5.0D0 / 2.0D0) * RGS_Ru / RG_MOLM(rg_i)
      else
        RG_CVS(rg_i) = (3.0D0 / 2.0D0) * RGS_Ru / RG_MOLM(rg_i)
      end if
      rg_tr = rg_tr + RG_VF(rg_i) * RG_CVS(rg_i) * rg_Ttr0
    end do


! Vibrational energy
    RG_EV_TOTAL= 0.0D0
    do rg_i = 1, 3
      RG_EV_TOTAL = RG_EV_TOTAL + RG_VF(rg_i)  * (RGS_Ru / RG_MOLM(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i)/rg_Tve0) - 1.0D0))
    end do

RG_CHEM=zero

  ! Chemical energy
 do rg_i=1,nof_species
        if (rg_hzero(RG_i).gt.1.0e-12)then
        RG_CHEM=RG_CHEM-(RG_VF(rg_i)*rg_hzero(RG_i)/RG_MOLM(rg_i))
        end if
END DO
!RG_CHEM=zero





! Kinetic energy

SKIN=(oo2)*((U**2)+(V**2))


inflow2d(1)=R
inflow2d(2)=R*U
inflow2d(3)=R*V
inflow2d(4)=R*(RG_EV_TOTAL+RG_TR+RG_CHEM+skin)
inflow2d(5)=R*RG_EV_TOTAL

do rg_i=1,nof_species
inflow2d(5+rg_i)=r * RG_VF(rg_i)
end do



end if






END FUNCTION INFLOW2d


FUNCTION OUTFLOW2d(INITCOND,POX,POY)
!> @brief
!> This function applies a prescribed boundary condition to  the outflow in 2D
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::OUTFLOW2d
INTEGER,INTENT(IN)::INITCOND
REAL,DIMENSION(1:2),INTENT(IN)::POX,POY
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI
REAL::XF,YF,ZF,GAMMAR
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE


IF (multispecies.EQ.1) then



P=PRES
U=uvel
V=vvel
MP_AR(1)=MP_A_IN(1)/(GAMMA_IN(1)-1.0D0)  
MP_AR(2)=MP_A_IN(2)/(GAMMA_IN(2)-1.0D0)
GAMMAR=(1.0D0/(MP_AR(1)+MP_AR(2)))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION

GM=GAMMAR

R=(MP_R_IN(1)*MP_A_IN(1))+(MP_R_IN(2)*MP_A_IN(2))
MP_IE(1)=((P+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IEn=(MP_IE(1)*MP_A_IN(1))+(MP_IE(2)*MP_A_IN(2))
! !KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2))
! !TOTAL ENERGY
E=(R*SKIN)+IEn

!VECTOR OF CONSERVED VARIABLES NOW
OUTFLOW2d(1)=R
OUTFLOW2d(2)=R*U
OUTFLOW2d(3)=R*V
OUTFLOW2d(4)=E
OUTFLOW2d(5)=MP_R_IN(1)*MP_A_IN(1)
OUTFLOW2d(6)=MP_R_IN(2)*MP_A_IN(2)
OUTFLOW2d(7)=MP_A_IN(1)

ELSE

R=RRES
GM=GAMMA
P=PRES
U=uvel
V=vvel


!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)
!VECTOR OF CONSERVED VARIABLES NOW
OUTFLOW2d(1)=R
OUTFLOW2d(2)=R*U
OUTFLOW2d(3)=R*V
OUTFLOW2d(4)=E
















END IF

END FUNCTION OUTFLOW2d

FUNCTION OUTFLOW(INITCOND,POX,POY,POZ)
!> @brief
!> This function applies a prescribed boundary condition to  the outflow in 3D
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::OUTFLOW
INTEGER,INTENT(IN)::INITCOND
REAL,DIMENSION(1:dimensiona),INTENT(IN)::POX,POY,POZ
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI
REAL::XF,YF,ZF,GAMMAR
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE

IF (multispecies.EQ.1) then



P=PRES
U=uvel
V=vvel
w=wvel
MP_AR(1)=MP_A_IN(1)/(GAMMA_IN(1)-1.0D0)  
MP_AR(2)=MP_A_IN(2)/(GAMMA_IN(2)-1.0D0)
GAMMAR=(1.0D0/(MP_AR(1)+MP_AR(2)))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION

GM=GAMMAR

R=(MP_R_IN(1)*MP_A_IN(1))+(MP_R_IN(2)*MP_A_IN(2))
MP_IE(1)=((P+(GAMMA_IN(1)*MP_PINF(1)))/((GAMMA_IN(1)-1.0D0)))
MP_IE(2)=((P+(GAMMA_IN(2)*MP_PINF(2)))/((GAMMA_IN(2)-1.0D0)))
IEn=(MP_IE(1)*MP_A_IN(1))+(MP_IE(2)*MP_A_IN(2))
! !KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(w**2))
! !TOTAL ENERGY
E=(R*SKIN)+IEn

!VECTOR OF CONSERVED VARIABLES NOW
OUTFLOW(1)=R
OUTFLOW(2)=R*U
OUTFLOW(3)=R*V
OUTFLOW(4)=R*w
OUTFLOW(5)=E
OUTFLOW(6)=MP_R_IN(1)*MP_A_IN(1)
OUTFLOW(7)=MP_R_IN(2)*MP_A_IN(2)
OUTFLOW(8)=MP_A_IN(1)

ELSE

R=RRES
GM=GAMMA
P=PRES

if (initcond.eq.977)then
 p=101325
 end if 

U=uvel
V=vvel
W=wvel

!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)
!VECTOR OF CONSERVED VARIABLES NOW
OUTFLOW(1)=R
OUTFLOW(2)=R*U
OUTFLOW(3)=R*V
OUTFLOW(4)=R*W
OUTFLOW(5)=E

end if

END FUNCTION OUTFLOW


FUNCTION OUTFLOW2(INITCOND,POX,POY,POZ)
!> @brief
!> This function applies a prescribed boundary condition to  the outflow in 3D
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::OUTFLOW2
INTEGER,INTENT(IN)::INITCOND
REAL,DIMENSION(1:dimensiona),INTENT(IN)::POX,POY,POZ
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI
REAL::XF,YF,ZF
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE




R=RRES
GM=GAMMA
P=PRESS_OUTLET
U=uvel
V=vvel
W=wvel






!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)
!VECTOR OF CONSERVED VARIABLES NOW
OUTFLOW2(1)=R
OUTFLOW2(2)=R*U
OUTFLOW2(3)=R*V
OUTFLOW2(4)=R*W
OUTFLOW2(5)=E



END FUNCTION OUTFLOW2



FUNCTION BLEED2D(Iconsidered,facex,POX,POY)
!> @brief
!> This function applies a prescribed boundary condition to  the outflow in 3D
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::BLEED2D
INTEGER,INTENT(IN)::iconsidered, facex
REAL,dimension(1:dimensiona),INTENT(IN)::POX,POY
INTEGER::IBLEEDn
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI
REAL::XF,YF,ZF
REAL::CELL_AREA
!------COEFFICIENTS FOR BLEED-----!
REAL::BLEED_QSONIC_S,BLEED_MDOTSONIC_S,BLEED_AREA,BLEED_MDOTSONIC,bleed_region
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL::ANGLE1,ANGLE2
REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::CLEFT,CRIGHT,CLEFT_ROT,CRIGHT_ROT






CELL_AREA=IELEM(N,Iconsidered)%SURF(facex)




!LEFTV IS THE INPUT IN CONSERVATIVE VARIABLES THAT WE TRANSFORM TO PRIMITIVE VARIABLES
CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
R=LEFTV(1)
U=LEFTV(2)
V=LEFTV(3)
P=LEFTV(4)







!NOW FIND THE BLEED NUMBER CONDITIONS FOR THIS FACE
IBLEEDn=IELEM(N,ICONSIDERED)%BLEEDN(FACEX)



bleed_region=sqrt(((bleed_start(IBLEEDn,1)-bleed_end(IBLEEDn,1))**2)+((bleed_start(IBLEEDn,2)-bleed_end(IBLEEDn,2))**2))

!NOW COMPUTE THE BLEED AREA

BLEED_AREA=BLEED_POROSITY(IBLEEDN)*BLEED_REGION

!NOW COMPUTE THE BLEED M DOT SONIC-S MASS FLOW RATE EQ.17 , https://doi.org/10.2514/1.B37474

BLEED_MDOTSONIC_S=BLEED_AREA*P*(SQRT(GAMMA*R/P))*(((GAMMA+1.0D0)/(2))**((GAMMA+1)/(2*(1-GAMMA))))

!NOW COMPUTE QSONIC EQ. 22
BLEED_QSONIC_S=0.598+0.0307*(BLEED_PLENUM(IBLEEDN)/P)-0.5936*((BLEED_PLENUM(IBLEEDN)/P)**2)

!equation 16
BLEED_MDOTSONIC=BLEED_QSONIC_S*BLEED_MDOTSONIC_S

!NOW THAT I HAVE COMPUTED THE BLEED MASS FLOW RATE I HAVE TO DISTRIBUTE TO THE SURFACE AREA OF THIS TAGGED CELL IN THE NORMAL DIRECTION




CALL PRIM2CONS2(N,LEFTV,RIGHTV)

RIGHTV(1:NOF_VARIABLES)=LEFTV(1:NOF_VARIABLES)


CALL ROTATEF2d(N,CRIGHT_ROT,RIGHTV,ANGLE1,ANGLE2)



CRIGHT_ROT(2)=BLEED_MDOTSONIC/CELL_AREA



CALL ROTATEb2d(N,rightv,Cright_ROT,ANGLE1,ANGLE2)




BLEED2D(1:NOF_VARIABLES)=RIGHTV(1:NOF_VARIABLES)


END FUNCTION BLEED2D


FUNCTION BLEED3D(Iconsidered,facex,POX,POY,poz)
!> @brief
!> This function applies a prescribed boundary condition to  the outflow in 3D
IMPLICIT NONE
REAL,DIMENSION(1:nof_Variables)::BLEED3D
INTEGER,INTENT(IN)::iconsidered, facex
REAL,dimension(1:dimensiona),INTENT(IN)::POX,POY,poz
INTEGER::IBLEEDn
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI
REAL::XF,YF,ZF
REAL::CELL_AREA
!------COEFFICIENTS FOR BLEED-----!
REAL::BLEED_QSONIC_S,BLEED_MDOTSONIC_S,BLEED_AREA,BLEED_MDOTSONIC,bleed_region
real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
REAL::ANGLE1,ANGLE2
REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR)::CLEFT,CRIGHT,CLEFT_ROT,CRIGHT_ROT



PRINT*,"NOT READY YET"

STOP
BLEED3D(1:NOF_VARIABLES)=RIGHTV(1:NOF_VARIABLES)


END FUNCTION BLEED3D




FUNCTION PASS_INLET(INITCOND,POX,POY,POZ)
!> @brief
!> This function applies a prescribed boundary condition to  the inlet for a passive scalar
IMPLICIT NONE
REAL,DIMENSION(1:PASSIVESCALAR)::PASS_INLET
INTEGER,INTENT(IN)::INITCOND
REAL,DIMENSION(1:dimensiona),INTENT(IN)::POX,POY,POZ


PASS_INLET(1:PASSIVESCALAR)=1.0D0*RRES

END FUNCTION PASS_INLET


FUNCTION PASS_INLET2d(INITCOND,POX,POY)
!> @brief
!> This function applies a prescribed boundary condition to  the inlet for a passive scalar in 2d
IMPLICIT NONE
REAL,DIMENSION(1:PASSIVESCALAR)::PASS_INLET2d
INTEGER,INTENT(IN)::INITCOND
REAL,DIMENSION(1:dimensiona),INTENT(IN)::POX,POY



PASS_INLET2d(1:PASSIVESCALAR)=1.0D0

END FUNCTION PASS_INLET2d









SUBROUTINE SHEAR_X(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the shear stresses in x-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml

 

 
VORTET1(1:3,1:3) = ILOCAL_RECON3(ICONSIDERED)%GRADS(1:3,1:3)


 ux = Vortet1(1,1);uy = Vortet1(1,2);uz = Vortet1(1,3)
 vx = Vortet1(2,1);vy = Vortet1(2,2);vz = Vortet1(2,3)
 wx = Vortet1(3,1);wy = Vortet1(3,2);wz = Vortet1(3,3)


 ANGLE1=IELEM(N,ICONSIDERED)%FACEANGLEX(FACEX)
 ANGLE2=IELEM(N,ICONSIDERED)%FACEANGLEY(FACEX)
 NX=(COS(ANGLE1)*SIN(ANGLE2))
 NY=(SIN(ANGLE1)*SIN(ANGLE2))
 NZ=(COS(ANGLE2))
 LEFTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)
 RIGHTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)

 CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

SSX=ZERO; SSP=ZERO; SSY=ZERO; SSZ=ZERO

TAUXX=(4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
TAUYY=(4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
TAUZZ=(4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
TAUYX=(UY + VX)
TAUZX=(WX + UZ)
TAUZY=(VZ + WY)
SSX=(VISCL(1)*((NX*TAUXX)+(NY*TAUYX)+(NZ*TAUZX)))
SSY=(VISCL(1)*((NX*TAUYX)+(NY*TAUYY)+(NZ*TAUZY)))
SSZ=(VISCL(1)*((NX*TAUZX)+(NY*TAUZY)+(NZ*TAUZZ)))




IF(RFRAME.EQ.0)THEN
SHEAR_TEMP=-SSX/(0.5*rres*ufreestream*ufreestream)
ELSE
SHEAR_TEMP=-SSX/(0.5*rres*V_REF*V_REF)
END IF

END SUBROUTINE SHEAR_X




SUBROUTINE SHEAR_Y(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the shear stresses in Y-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 
VORTET1(1:3,1:3) = ILOCAL_RECON3(ICONSIDERED)%GRADS(1:3,1:3)


 ux = Vortet1(1,1);uy = Vortet1(1,2);uz = Vortet1(1,3)
 vx = Vortet1(2,1);vy = Vortet1(2,2);vz = Vortet1(2,3)
 wx = Vortet1(3,1);wy = Vortet1(3,2);wz = Vortet1(3,3)


 ANGLE1=IELEM(N,ICONSIDERED)%FACEANGLEX(FACEX)
 ANGLE2=IELEM(N,ICONSIDERED)%FACEANGLEY(FACEX)
 NX=(COS(ANGLE1)*SIN(ANGLE2))
 NY=(SIN(ANGLE1)*SIN(ANGLE2))
 NZ=(COS(ANGLE2))
 LEFTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)
 RIGHTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)

 CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

SSX=ZERO; SSP=ZERO; SSY=ZERO; SSZ=ZERO

TAUXX=(4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
TAUYY=(4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
TAUZZ=(4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
TAUYX=(UY + VX)
TAUZX=(WX + UZ)
TAUZY=(VZ + WY)
SSX=(VISCL(1)*((NX*TAUXX)+(NY*TAUYX)+(NZ*TAUZX)))
SSY=(VISCL(1)*((NX*TAUYX)+(NY*TAUYY)+(NZ*TAUZY)))
SSZ=(VISCL(1)*((NX*TAUZX)+(NY*TAUZY)+(NZ*TAUZZ)))


IF(RFRAME.EQ.0)THEN
SHEAR_TEMP=-SSY/(0.5*rres*ufreestream*ufreestream)
ELSE
SHEAR_TEMP=-SSY/(0.5*rres*V_REF*V_REF)
END IF

END SUBROUTINE SHEAR_Y


SUBROUTINE SHEAR_Z(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the shear stresses in Z-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 
VORTET1(1:3,1:3) = ILOCAL_RECON3(ICONSIDERED)%GRADS(1:3,1:3)


 ux = Vortet1(1,1);uy = Vortet1(1,2);uz = Vortet1(1,3)
 vx = Vortet1(2,1);vy = Vortet1(2,2);vz = Vortet1(2,3)
 wx = Vortet1(3,1);wy = Vortet1(3,2);wz = Vortet1(3,3)


 ANGLE1=IELEM(N,ICONSIDERED)%FACEANGLEX(FACEX)
 ANGLE2=IELEM(N,ICONSIDERED)%FACEANGLEY(FACEX)
 NX=(COS(ANGLE1)*SIN(ANGLE2))
 NY=(SIN(ANGLE1)*SIN(ANGLE2))
 NZ=(COS(ANGLE2))
 LEFTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)
 RIGHTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)

 CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

SSX=ZERO; SSP=ZERO; SSY=ZERO; SSZ=ZERO

TAUXX=(4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
TAUYY=(4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
TAUZZ=(4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
TAUYX=(UY + VX)
TAUZX=(WX + UZ)
TAUZY=(VZ + WY)
SSX=(VISCL(1)*((NX*TAUXX)+(NY*TAUYX)+(NZ*TAUZX)))
SSY=(VISCL(1)*((NX*TAUYX)+(NY*TAUYY)+(NZ*TAUZY)))
SSZ=(VISCL(1)*((NX*TAUZX)+(NY*TAUZY)+(NZ*TAUZZ)))


IF(RFRAME.EQ.0)THEN
SHEAR_TEMP=-SSZ/(0.5*rres*ufreestream*ufreestream)
ELSE
SHEAR_TEMP=-SSZ/(0.5*rres*V_REF*V_REF)
END IF

END SUBROUTINE SHEAR_Z







SUBROUTINE HEAT_X(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the shear stresses in x-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP,lam_qflux
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 INTEGER::I,K,J,KMAXE,gqi_points,nnd,IM
 real,dimension(1:nof_Variables)::leftv
 real,dimension(1:3)::TEMP_gRAD
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz,surface_temp
real,dimension(1:4)::viscl,laml
REAL,DIMENSION(1:2)::TURBMV
REAL,DIMENSION(1)::ETVM
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
REAL,DIMENSION(1:DIMENSIONA,1:NUMBEROFPOINTS2)::QPOINTS2D
REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D


I=ICONSIDERED
SSX=ZERO;SSy=ZERO;SSz=ZERO
J=FACEX
			      ANGLE1=IELEM(N,I)%FACEANGLEX(j)
			      ANGLE2=IELEM(N,I)%FACEANGLEY(j)
			      NX=(COS(ANGLE1)*SIN(ANGLE2))
			      NY=(SIN(ANGLE1)*SIN(ANGLE2))
			      NZ=(COS(ANGLE2))

                select case(ielem(n,i)%types_faces(j))
				case (5)
					  gqi_points=qp_quad_n



					  if(reduce_comp.eq.1)then
					  WEqua2d=1.0d0;
					  else
					    NND=4
				      do K=1,nnd
					VEXT(k,1:dims)=inoder4(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
				      END DO
					  call  QUADRATUREQUAD3D(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					  end if
					  surface_temp=IELEM(N,I)%SURF(J)


				case(6)
					gqi_points=qp_triangle_n


					if(reduce_comp.eq.1)then
					  WEqua2d=1.0d0;
					  else
					  NND=3
					do K=1,nnd
					  VEXT(k,1:dims)=inoder4(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
					END DO
					call QUADRATURETRIANG(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					end if
 					    surface_temp=IELEM(N,I)%SURF(J)



				end select




				do im=1,gqi_points
				TEMP_gRAD(1:3)=ILOCAL_RECON3(i)%ULEFTV(1:3,dimensiona+1,J,IM)
				IF (DG.EQ.1)THEN
				  LEFTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT_DG(1:nof_Variables, J,IM)
				  RIGHTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT_DG(1:nof_Variables, J,IM)


				  ELSE
				  LEFTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  RIGHTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  END IF


					CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

                              IF (TURBULENCEMODEL.EQ.1)THEN

                              TURBMV(1)=ILOCAL_RECON3(I)%ULEFTTURB(1,j,im)

							  TURBMV(2)=ILOCAL_RECON3(I)%ULEFTTURB(1,j,im)
							  eddyfl(2)=turbmv(1);
							  eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF

					if (turbulence .eq. 1) then
					lam_qflux=LAML(3)
					else
					lam_qflux=LAML(1)
					end if



				  SSX=SSX+lam_qflux*TEMP_gRAD(1)*WEQUA2D(im)*nx!*surface_temp
				  SSy=SSy+lam_qflux*TEMP_gRAD(2)*WEQUA2D(im)*ny!*surface_temp
				  SSz=SSz+lam_qflux*TEMP_gRAD(3)*WEQUA2D(im)*nz!*surface_temp
               END DO








SHEAR_TEMP=-(SSX+ssy+ssz)



END SUBROUTINE HEAT_X








SUBROUTINE HEAT_X2D(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the shear stresses in x-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 INTEGER::I,K,J,KMAXE,gqi_points,nnd,IM
 real,dimension(1:nof_Variables)::leftv
 real,dimension(1:2)::TEMP_gRAD
real::MP_PINFL,gammal,lam_qflux
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz,surface_temp
real,dimension(1:4)::viscl,laml
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
REAL,DIMENSION(1:DIMENSIONA,1:NUMBEROFPOINTS2)::QPOINTS2D
REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
REAL,DIMENSION(1:2)::TURBMV
REAL,DIMENSION(1)::ETVM
REAL,DIMENSION(1:20)::EDDYFL,EDDYFR


I=ICONSIDERED
SSX=ZERO;ssy=zero
J=FACEX
			     nx=IELEM(N,I)%FACEANGLEX(j)
			      ny=IELEM(N,I)%FACEANGLEY(j)

                gqi_points=qp_line_n
					   if(reduce_comp.eq.1)then
					  WEqua2d=1.0d0;
					  else
					  NND=2
				      do K=1,nnd
					VEXT(k,1:dims)=inoder4(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
				      END DO

					  call  QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					  end if
					  surface_temp=IELEM(N,I)%SURF(J)




				do im=1,gqi_points
				TEMP_gRAD(1:2)=ILOCAL_RECON3(i)%ULEFTV(1:2,dimensiona+1,J,IM)

				IF (DG.EQ.1)THEN
				  LEFTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT_DG(1:nof_Variables, J,IM)
				  RIGHTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT_DG(1:nof_Variables, J,IM)


				  ELSE
				  LEFTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  RIGHTV(1:nof_Variables)=ILOCAL_RECON3(I)%ULEFT(:,j,im)
				  END IF

				  CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

                              IF (TURBULENCEMODEL.EQ.1)THEN

                              TURBMV(1)=ILOCAL_RECON3(I)%ULEFTTURB(1,j,im)

							  TURBMV(2)=ILOCAL_RECON3(I)%ULEFTTURB(1,j,im)
							  eddyfl(2)=turbmv(1);
							  eddyfr(2)=turbmv(2)
							  Call EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
						      END IF

					if (turbulence .eq. 1) then
					lam_qflux=LAML(3)
					else
					lam_qflux=LAML(1)
					end if

                      SSX=SSX+lam_qflux*TEMP_gRAD(1)*WEQUA2D(im)*nx!*surface_temp
                      SSy=SSy+lam_qflux*TEMP_gRAD(2)*WEQUA2D(im)*ny!*surface_temp





               END DO



SHEAR_TEMP=-(SSX+ssy)



END SUBROUTINE HEAT_X2D


SUBROUTINE HEAT_Y2D(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the shear stresses in x-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 INTEGER::I,K,J,KMAXE,gqi_points,nnd,IM
 real,dimension(1:nof_Variables)::leftv
 real,dimension(1:3)::TEMP_gRAD
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz,surface_temp
real,dimension(1:4)::viscl,laml
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
REAL,DIMENSION(1:DIMENSIONA,1:NUMBEROFPOINTS2)::QPOINTS2D
REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D


I=ICONSIDERED
SSY=ZERO
J=FACEX
			     nx=IELEM(N,I)%FACEANGLEX(j)
			      ny=IELEM(N,I)%FACEANGLEY(j)

                gqi_points=qp_line_n
					   if(reduce_comp.eq.1)then
					  WEqua2d=1.0d0;
					  else
					  NND=2
				      do K=1,nnd
					VEXT(k,1:dims)=inoder4(IELEM(N,I)%NODES_FACES(J,K))%CORD(1:dims)
				      END DO

					  call  QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
					  end if
					  surface_temp=IELEM(N,I)%SURF(J)




				do im=1,gqi_points
				TEMP_gRAD(1:2)=ILOCAL_RECON3(i)%ULEFTV(1:2,dimensiona+1,J,IM)
				  SSY=SSY-0.026*TEMP_gRAD(2)*WEQUA2D(im)*surface_temp
               END DO



SHEAR_TEMP=SSY



END SUBROUTINE HEAT_Y2D




















SUBROUTINE SHEAR_X_av(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the AVERAGE shear stresses in x-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 integer::ind1
 
 
 if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if

VORTET1(1:3,1:3) = ILOCAL_RECON3(ICONSIDERED)%GRADSAV(1:3,1:3)


 ux = Vortet1(1,1);uy = Vortet1(1,2);uz = Vortet1(1,3)
 vx = Vortet1(2,1);vy = Vortet1(2,2);vz = Vortet1(2,3)
 wx = Vortet1(3,1);wy = Vortet1(3,2);wz = Vortet1(3,3)


 ANGLE1=IELEM(N,ICONSIDERED)%FACEANGLEX(FACEX)
 ANGLE2=IELEM(N,ICONSIDERED)%FACEANGLEY(FACEX)
 NX=(COS(ANGLE1)*SIN(ANGLE2))
 NY=(SIN(ANGLE1)*SIN(ANGLE2))
 NZ=(COS(ANGLE2))
 LEFTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(IND1,1:nof_Variables)
 RIGHTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(IND1,1:nof_Variables)

 CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

SSX=ZERO; SSP=ZERO; SSY=ZERO; SSZ=ZERO

TAUXX=(4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
TAUYY=(4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
TAUZZ=(4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
TAUYX=(UY + VX)
TAUZX=(WX + UZ)
TAUZY=(VZ + WY)
SSX=(VISCL(1)*((NX*TAUXX)+(NY*TAUYX)+(NZ*TAUZX)))
SSY=(VISCL(1)*((NX*TAUYX)+(NY*TAUYY)+(NZ*TAUZY)))
SSZ=(VISCL(1)*((NX*TAUZX)+(NY*TAUZY)+(NZ*TAUZZ)))


SHEAR_TEMP=-SSX/(0.5*rres*ufreestream*ufreestream)



END SUBROUTINE SHEAR_X_av



SUBROUTINE SHEAR_Y_av(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the AVERAGE shear stresses in Y-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 integer::ind1

 
 
 if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if
 
VORTET1(1:3,1:3) = ILOCAL_RECON3(ICONSIDERED)%GRADSAV(1:3,1:3)


 ux = Vortet1(1,1);uy = Vortet1(1,2);uz = Vortet1(1,3)
 vx = Vortet1(2,1);vy = Vortet1(2,2);vz = Vortet1(2,3)
 wx = Vortet1(3,1);wy = Vortet1(3,2);wz = Vortet1(3,3)


 ANGLE1=IELEM(N,ICONSIDERED)%FACEANGLEX(FACEX)
 ANGLE2=IELEM(N,ICONSIDERED)%FACEANGLEY(FACEX)
 NX=(COS(ANGLE1)*SIN(ANGLE2))
 NY=(SIN(ANGLE1)*SIN(ANGLE2))
 NZ=(COS(ANGLE2))
 LEFTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(IND1,1:nof_Variables)
 RIGHTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(IND1,1:nof_Variables)

 CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

SSX=ZERO; SSP=ZERO; SSY=ZERO; SSZ=ZERO

TAUXX=(4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
TAUYY=(4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
TAUZZ=(4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
TAUYX=(UY + VX)
TAUZX=(WX + UZ)
TAUZY=(VZ + WY)
SSX=(VISCL(1)*((NX*TAUXX)+(NY*TAUYX)+(NZ*TAUZX)))
SSY=(VISCL(1)*((NX*TAUYX)+(NY*TAUYY)+(NZ*TAUZY)))
SSZ=(VISCL(1)*((NX*TAUZX)+(NY*TAUZY)+(NZ*TAUZZ)))


SHEAR_TEMP=-SSY/(0.5*rres*ufreestream*ufreestream)



END SUBROUTINE SHEAR_Y_av


SUBROUTINE SHEAR_Z_av(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the AVERGAGE shear stresses in Z-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 integer::ind1
 
 
 if (rungekutta.eq.4)then
	      ind1=7
	      else
	      ind1=5
	      end if
 
VORTET1(1:3,1:3) = ILOCAL_RECON3(ICONSIDERED)%GRADSAV(1:3,1:3)


 ux = Vortet1(1,1);uy = Vortet1(1,2);uz = Vortet1(1,3)
 vx = Vortet1(2,1);vy = Vortet1(2,2);vz = Vortet1(2,3)
 wx = Vortet1(3,1);wy = Vortet1(3,2);wz = Vortet1(3,3)


 ANGLE1=IELEM(N,ICONSIDERED)%FACEANGLEX(FACEX)
 ANGLE2=IELEM(N,ICONSIDERED)%FACEANGLEY(FACEX)
 NX=(COS(ANGLE1)*SIN(ANGLE2))
 NY=(SIN(ANGLE1)*SIN(ANGLE2))
 NZ=(COS(ANGLE2))
 LEFTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(IND1,1:nof_Variables)
 RIGHTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(IND1,1:nof_Variables)

 CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

SSX=ZERO; SSP=ZERO; SSY=ZERO; SSZ=ZERO

TAUXX=(4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
TAUYY=(4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
TAUZZ=(4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY
TAUYX=(UY + VX)
TAUZX=(WX + UZ)
TAUZY=(VZ + WY)
SSX=(VISCL(1)*((NX*TAUXX)+(NY*TAUYX)+(NZ*TAUZX)))
SSY=(VISCL(1)*((NX*TAUYX)+(NY*TAUYY)+(NZ*TAUZY)))
SSZ=(VISCL(1)*((NX*TAUZX)+(NY*TAUZY)+(NZ*TAUZZ)))


SHEAR_TEMP=-SSZ/(0.5*rres*ufreestream*ufreestream)


END SUBROUTINE SHEAR_Z_av

SUBROUTINE SHEAR_X2d(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the shear stresses in x-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 integer::ind1
integer::gqi_points,im
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
 SSX=zero; SSP=zero; SSY=zero; 
 
 gqi_points=qp_line_n
CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
 
 do im=1,gqi_points
 if (ielem(n,iconsiDERED)%ggs.eq.1)then
VORTET1(1:2,1:2) = ILOCAL_RECON3(ICONSIDERED)%GRADS(1:2,1:2)
else
 vortet1(1,1:2)=ILOCAL_RECON3(iconsiDERED)%ULEFTV(1:2,1,FACEX,IM)
vortet1(2,1:2)=ILOCAL_RECON3(iconsiDERED)%ULEFTV(1:2,2,FACEX,IM)

end if

 ux = Vortet1(1,1);uy = Vortet1(1,2)
 vx = Vortet1(2,1);vy = Vortet1(2,2)
 


 ANGLE1=IELEM(N,ICONSIDERED)%FACEANGLEX(FACEX)
 ANGLE2=IELEM(N,ICONSIDERED)%FACEANGLEY(FACEX)
 NX=ANGLE1
 NY=ANGLE2
 
 LEFTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)
 RIGHTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)

 CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

SSX=ZERO; SSP=ZERO; SSY=ZERO; SSZ=ZERO

TAUXX=2.0D0*UX
TAUYY=2.0D0*VY

TAUYX=(UY + VX)

SSX=SSX+((VISCL(1)*((NY*TAUYX)))*WEQUA2D(im))
SSY=SSY+((VISCL(1)*((NX*TAUYX)))*WEQUA2D(im))
end do

SHEAR_TEMP=-SSX/(0.5*rres*ufreestream*ufreestream)




END SUBROUTINE SHEAR_X2d



SUBROUTINE SHEAR_Y2d(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the shear stresses in Y-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,TAUXX,TAUYY,TAUZZ,TAUYX,TAUZX,TAUZY
 REAL::SSX,SSY,SSZ,SSP
 REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml
 integer::ind1
integer::gqi_points,im
REAL,DIMENSION(1:8,1:DIMENSIONA)::VEXT
	REAL,DIMENSION(1:dimensiona,1:NUMBEROFPOINTS2)::QPOINTS2D
	REAL,DIMENSION(1:NUMBEROFPOINTS2)::WEQUA2D
 SSX=zero; SSP=zero; SSY=zero; 
 
 gqi_points=qp_line_n
CALL QUADRATURELINE(N,IGQRULES,VEXT,QPOINTS2D,WEQUA2D)
 
 do im=1,gqi_points
 if (ielem(n,iconsiDERED)%ggs.eq.1)then
VORTET1(1:2,1:2) = ILOCAL_RECON3(ICONSIDERED)%GRADS(1:2,1:2)
else
 vortet1(1,1:2)=ILOCAL_RECON3(iconsiDERED)%ULEFTV(1:2,1,FACEX,IM)
vortet1(2,1:2)=ILOCAL_RECON3(iconsiDERED)%ULEFTV(1:2,2,FACEX,IM)

end if

 ux = Vortet1(1,1);uy = Vortet1(1,2)
 vx = Vortet1(2,1);vy = Vortet1(2,2)
 


 ANGLE1=IELEM(N,ICONSIDERED)%FACEANGLEX(FACEX)
 ANGLE2=IELEM(N,ICONSIDERED)%FACEANGLEY(FACEX)
 NX=ANGLE1
 NY=ANGLE2
 
 LEFTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)
 RIGHTV(1:nof_Variables)=U_C(ICONSIDERED)%VAL(1,1:nof_Variables)

 CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)

SSX=ZERO; SSP=ZERO; SSY=ZERO; SSZ=ZERO

TAUXX=2.0D0*UX
TAUYY=2.0D0*VY

TAUYX=(UY + VX)

SSX=SSX+((VISCL(1)*((NY*TAUYX)))*WEQUA2D(im))
SSY=SSY+((VISCL(1)*((NX*TAUYX)))*WEQUA2D(im))
end do

SHEAR_TEMP=-SSX/(0.5*rres*ufreestream*ufreestream)


END SUBROUTINE SHEAR_Y2d




SUBROUTINE SHEAR_X2d_av(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the AVERAGE shear stresses in x-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP





SHEAR_TEMP=0.0D0


END SUBROUTINE SHEAR_X2d_av



SUBROUTINE SHEAR_Y2d_av(ICONSIDERED,FACEX,SHEAR_TEMP)
!> @brief
!> This subroutine computes the AVERAGE shear stresses in Y-axis
IMPLICIT NONE
INTEGER,INTENT(IN)::ICONSIDERED,FACEX
REAL,INTENT(INOUT)::SHEAR_TEMP

SHEAR_TEMP=0.0D0


END SUBROUTINE SHEAR_Y2d_av



SUBROUTINE GET_visc_conduct(N,leftv,rightv,VISCL,LAML)
!> @brief
!> This subroutine computes the viscosity and thermal conductivity according to sutherland's law or others
	IMPLICIT NONE
        REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::LEFTV,RIGHTV
        REAL,DIMENSION(1:4),INTENT(INOUT)::VISCL,LAML
	INTEGER,INTENT(IN)::N
	REAL::KINETIC,U,V,W,T0L,T1L,T0R,T1R
	REAL::MP_Temp,MP_mu_mix,MP_k_mix,MP_Cp_mix,gammal,GAMMAR
	REAL,DIMENSION(1:NOF_SPECIES)::MP_D_eff,mp_mu_i,MP_ktr_i
	REAL::MP_ktr_mix,MP_kve,mp_pinfl,mp_pinfr
	REAL,DIMENSION(1:NOF_VARIABLES)::LEFT_Temp,RIGHT_temp


          LEFT_Temp=leftv
          right_temp=rightv
	if ((multispecies.eq.1).or.(realgas.eq.1))then

          if (multispecies.eq.1)then

              CALL cons2prim(N,LEFT_Temp,MP_PINFl,gammal)
              CALL PRIM2CONS(N,LEFT_Temp)

              CALL MULTISPECIES_MIXTURES(LEFT_Temp,MP_Temp,MP_mu_mix,MP_k_mix,MP_Cp_mix,gammal)

              VISCL(1)=MP_mu_mix
              LAML(1)=MP_k_mix

              CALL cons2prim(N,right_temp,MP_PINFR,gammaR)
              CALL PRIM2CONS(N,right_temp)
              CALL MULTISPECIES_MIXTURES(right_temp,MP_Temp,MP_mu_mix,MP_k_mix,MP_Cp_mix,gammaR)

              VISCL(2)=MP_mu_mix
              LAML(2)=MP_k_mix

          end if


          if (realgas.eq.1)then
            call MULTISPECIES_MIXTURES_RG(LEFT_Temp,MP_mu_mix,MP_ktr_mix,MP_kve,MP_D_eff,mp_mu_i,MP_ktr_i,GAMMAL)
								LAML(1)=MP_ktr_mix
								VISCL(1)=MP_mu_mix

            call MULTISPECIES_MIXTURES_RG(right_temp,MP_mu_mix,MP_ktr_mix,MP_kve,MP_D_eff,mp_mu_i,MP_ktr_i,GAMMAr)
								LAML(2)=MP_ktr_mix
								VISCL(2)=MP_mu_mix


          end if



	else


        CALL CONS2PRIM2(N,LEFT_Temp,right_temp,MP_PINFl,MP_PINFr,gammal,gammar)
		
		T1L=LEFT_Temp(dimensiona+2)/(LEFT_Temp(1)*R_gas)
		T0L=PRES/(RRES*R_gas)
		
		
		T1R=right_temp(dimensiona+2)/(right_temp(1)*R_gas)
		T0R=PRES/(RRES*R_gas)

             
	      	
              VISCL(1)=VISC*((T1L/T0L)**BETAAS)*((T0L+(SUTHER*T0L))/(T1L+(SUTHER*T0L)))
              VISCL(2)=VISC*((T1R/T0R)**BETAAS)*((T0R+(SUTHER*T0R))/(T1R+(SUTHER*T0R)))

	      
	      LAML(1)=VISCL(1)*R_gas*GAMMA/(PRANDTL*(GAMMA-1.d0))
	      LAML(2)=VISCL(2)*R_gas*GAMMA/(PRANDTL*(GAMMA-1.d0))
	  


	  end if
	
	     
	   

      

  END SUBROUTINE GET_visc_conduct
  



  
SUBROUTINE VORTEXCALC(N)
!> @brief
!> This subroutine computes the q-criterion
  
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I,IHGT,IHGJ
REAL::SNORM,ONORM
REAL,DIMENSION(3,3)::TVORT,SVORT,OVORT
REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 	 KMAXE=XMPIELRANK(N)
!$OMP DO
DO I=1,KMAXE     
	    
                VORTET1(1:3,1:3)=ILOCAL_RECON3(I)%GRADS(1:3,1:3)    
	    
	    DO IHGT=1,3; DO IHGJ=1,3
	    TVORT(IHGT,IHGJ)=VORTET1(IHGJ,IHGT)
	      END DO; END DO
	      SVORT=0.5D0*(VORTET1+TVORT)
	      OVORT=0.5D0*(VORTET1-TVORT)
	      SNORM=SQRT((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+&
(SVORT(1,3)*SVORT(1,3))+(SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))+(SVORT(2,3)*SVORT(2,3))&
+(SVORT(3,1)*SVORT(3,1))+(SVORT(3,2)*SVORT(3,2))+(SVORT(3,3)*SVORT(3,3)))
	      ONORM=SQRT((OVORT(1,1)*OVORT(1,1))+(OVORT(1,2)*OVORT(1,2))+(OVORT(1,3)*OVORT(1,3))+&
(OVORT(2,1)*OVORT(2,1))+(OVORT(2,2)*OVORT(2,2))+(OVORT(2,3)*OVORT(2,3))+(OVORT(3,1)*OVORT(3,1))+&
(OVORT(3,2)*OVORT(3,2))+(OVORT(3,3)*OVORT(3,3)))
	      
	      IELEM(N,I)%VORTEX(1)=(0.5D0*((ONORM**2)-(SNORM**2)))
		
END DO
!$OMP END DO



END SUBROUTINE VORTEXCALC



SUBROUTINE ENSTROPHY_CALC(N)
!> @brief
!> This subroutine computes the q-criterion
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I,IHGT,IHGJ
REAL::SNORM,ONORM
REAL,DIMENSION(3,3)::TVORT,SVORT,OVORT
real,dimension(3,3)::taul,taur,TAU
REAL,DIMENSION(3)::Q,NNN,nall
REAL::UX,UY,UZ,VX,VY,VZ,WX,WY,WZ,RHO12,U12,V12,W12 ,damp,vdamp,TEMPXX
REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1
 real,dimension(1:nof_Variables)::leftv
real::MP_PINFL,gammal
real,dimension(1:nof_Variables)::RIGHTv
real::MP_PINFR,gammaR
real::angle1,angle2,nx,ny,nz
real,dimension(1:4)::viscl,laml


 	 KMAXE=XMPIELRANK(N)
!$OMP DO
DO I=1,KMAXE

                VORTET1(1:3,1:3)=ILOCAL_RECON3(I)%GRADS(1:3,1:3)

	    DO IHGT=1,3; DO IHGJ=1,3
	    TVORT(IHGT,IHGJ)=VORTET1(IHGJ,IHGT)
	      END DO; END DO

	      OVORT=(VORTET1-TVORT)
	      ONORM=((OVORT(1,1)*OVORT(1,1))+(OVORT(1,2)*OVORT(1,2))+(OVORT(1,3)*OVORT(1,3))+&
(OVORT(2,1)*OVORT(2,1))+(OVORT(2,2)*OVORT(2,2))+(OVORT(2,3)*OVORT(2,3))+(OVORT(3,1)*OVORT(3,1))+&
(OVORT(3,2)*OVORT(3,2))+(OVORT(3,3)*OVORT(3,3)))

           if(boundtype.eq.1)then

	      IELEM(N,I)%VORTEX(2)=(0.5D0*(ONORM*u_c(i)%val(1,1)))

	      else


	      LEFTV(1:NOF_vARIABLES)=U_C(I)%VAL(1,1:NOF_vARIABLES)
		RIGHTV(1:NOF_vARIABLES)=LEFTV(1:NOF_vARIABLES)
		CALL GET_visc_conduct(N,LEFTV,RIGHTV,VISCL,LAML)




                      UX = ILOCAL_RECON3(I)%GRADS(1,1); UY = ILOCAL_RECON3(I)%GRADS(1,2); UZ = ILOCAL_RECON3(I)%GRADS(1,3);
					  VX = ILOCAL_RECON3(I)%GRADS(2,1); VY = ILOCAL_RECON3(I)%GRADS(2,2); VZ = ILOCAL_RECON3(I)%GRADS(2,3);
					  WX = ILOCAL_RECON3(I)%GRADS(3,1); WY = ILOCAL_RECON3(I)%GRADS(3,2); WZ = ILOCAL_RECON3(I)%GRADS(3,3);






					  ! TAU_XX
					  TAUL(1,1) = (4.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY - (2.0D0/3.0D0)*WZ
					  ! TAU_YY
					  TAUL(2,2) = (4.0D0/3.0D0)*VY - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*WZ
					  ! TAU_ZZ
					  TAUL(3,3) = (4.0D0/3.0D0)*WZ - (2.0D0/3.0D0)*UX - (2.0D0/3.0D0)*VY

					  ! tau_xy
					  TAUL(1,2) = (UY + VX);TAUL(2,1) = TAUL(1,2)

					  ! TAU_XZ
					  TAUL(1,3) = (WX + UZ);TAUL(3,1) = TAUL(1,3)

					  ! TAU_YZ
					  TAUL(2,3) = (VZ + WY);TAUL(3,2) = TAUL(2,3)







          SNORM=((WY-VZ)**2)+((UZ-WX)**2)+((VX-UY)**2)
          ielem(n,i)%vortex(2)=ielem(n,i)%TOTVOLUME*SNORM*viscl(1)
          IELEM(N,I)%VORTEX(3)=(4.0/3.0)*VISCL(1)*((ux+vy+wz)**2)*ielem(n,i)%TOTVOLUME





	      end if

END DO

!$OMP END DO



END SUBROUTINE ENSTROPHY_CALC


SUBROUTINE VORTEXCALC2D(N)
!> @brief
!> This subroutine computes the q criterion for 2D
  
IMPLICIT NONE
INTEGER,INTENT(IN)::N
INTEGER::KMAXE,I,IHGT,IHGJ
REAL::SNORM,ONORM
REAL,DIMENSION(2,2)::TVORT,SVORT,OVORT
REAL,DIMENSION(1:DIMS,1:DIMS)::VORTET1

 	 KMAXE=XMPIELRANK(N)
 	 
 	
!$OMP DO
DO I=1,KMAXE     
	    
                VORTET1(1:2,1:2)=ILOCAL_RECON3(I)%GRADS(1:2,1:2)    
	    
	    DO IHGT=1,2; DO IHGJ=1,2
	    TVORT(IHGT,IHGJ)=VORTET1(IHGJ,IHGT)
	      END DO; END DO
	      SVORT=0.5D0*(VORTET1+TVORT)
	      OVORT=0.5D0*(VORTET1-TVORT)
	      SNORM=SQRT((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+&
 (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2)))
	      ONORM=SQRT((OVORT(1,1)*OVORT(1,1))+(OVORT(1,2)*OVORT(1,2))+&
(OVORT(2,1)*OVORT(2,1))+(OVORT(2,2)*OVORT(2,2)))
	      
	      IELEM(N,I)%VORTEX(1)=(0.5D0*((ONORM**2)-(SNORM**2)))
		
END DO
!$OMP END DO



END SUBROUTINE VORTEXCALC2D





SUBROUTINE BOUNDARYS(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
!> @brief
!> This subroutine applies the boundary condition to each bounded cell
implicit none
integer,intent(in)::n,b_code,ICONSIDERED,facex
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::LEFTV,RIGHTV
INTEGER,INTENT(INOUT)::IBFC
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::SRF_SPEEDROT,SRF_SPEED
REAL,DIMENSION(1:dimensiona),INTENT(IN)::POX,POY,POZ
REAL,INTENT(IN)::ANGLE1,ANGLE2,NX,NY,NZ
REAL,DIMENSION(TURBULENCEEQUATIONS),INTENT(INOUT)::CTURBL,CTURBR
REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR),INTENT(INOUT)::CRIGHT_ROT,CLEFT_ROT
REAL,DIMENSION(1:NOF_VARIABLES)::SUBSON1,SUBSON2,SUBSON3,tempxv,TEMPVECT1,TEMPVECT2
REAl::MP_PINFL,MP_PINFR,GAMMAL,GAMMAR,u,v,w
REAL::SPS,SKINS,IKINS,VEL,vnb,theeta,reeta
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX,VHX,AMP,DVEL,rgg,tt1
REAL::MN,UN,AGRT,rho_g,rg_rmix,Rmix,T_g,Tv_g,rg_tr,RG_EV_TOTAL,RG_CHEM
REAL,DIMENSION(1:NOF_SPECIES)::rho_s_g1,Y_s_g,RG_CVS
integer::rg_i





SELECT CASE(B_CODE)

    
    CASE(1)!INFLOW SUBSONIC OR SUPERSONIC WILL BE CHOSEN BASED ON MACH NUMBER
    if (boundtype.eq.0)then	!SUPERSONIC
    
    RIGHTV(1:nof_Variables)=INFLOW(INITCOND,POX,POY,POZ)
    
    
    
        
    
    
    ELSE		!SUBSONIC
    RIGHTV(1:nof_Variables)=INFLOW(INITCOND,POX,POY,POZ)
    
    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
    
    SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
    SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
    SPS=SQRT((GAMMA*SUBSON2(5))/(SUBSON2(1)))
    VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2+SUBSON2(4)**2)
    CALL PRIM2CONS2(N,LEFTV,RIGHTV)
    
      IF (VEL/(SPS+TOLSMALL).GT.1.0D0)THEN	!SUPERSONIC
      
      
      RIGHTV(1:nof_Variables)=INFLOW(INITCOND,POX,POY,POZ)
      
      
      
      ELSE		!SUBSONIC
      
      SUBSON3(5)=0.5*((SUBSON1(5))+(SUBSON2(5))-(SUBSON2(1)*SPS*((NX*(SUBSON1(2)-SUBSON2(2)))+(NY*(SUBSON1(3)-SUBSON2(3)))&
+(NZ*(SUBSON1(4)-SUBSON2(4))))))
      SUBSON3(1)=SUBSON1(1)+(SUBSON3(5)-SUBSON1(5))/(SPS**2)
      SUBSON3(2)=SUBSON1(2)-(NX*(SUBSON1(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
      SUBSON3(3)=SUBSON1(3)-(NY*(SUBSON1(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
      SUBSON3(4)=SUBSON1(4)-(NZ*(SUBSON1(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
      
     
      
       rightv(1)=SUBSON3(1)
    rightv(2)=SUBSON3(2)*SUBSON3(1)
    rightv(3)=SUBSON3(3)*SUBSON3(1)
    rightv(4)=SUBSON3(4)*SUBSON3(1)
    SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2)+(SUBSON3(4)**2))
    IKINS=SUBSON3(5)/((GAMMA-1.0d0)*(SUBSON3(1)))
    rightv(5)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)
      
      
      
      END IF

      
      
      
      
      
      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
	      
      
      
      
	      IF (TURBULENCEMODEL.EQ.1)THEN
		  CTURBR(1)=VISC*TURBINIT*rightv(1)
	      END IF
	      IF (TURBULENCEMODEL.EQ.2)THEN	 
		CTURBR(1)=(1.5D0*I_turb_inlet*(ufreestream**2))*RIGHTV(1)!K INITIALIZATION
		CTURBR(2)=RIGHTV(1)*CTURBR(1)/(10.0e-5*visc)!OMEGA INITIALIZATION
        IF (ILOCAL_RECON3(ICONSIDERED)%MRF.EQ.1)THEN
        CTURBR(1)=(1.5D0*I_turb_inlet*(KINIT_SRF**2))*RIGHTV(1)!K INITIALIZATION
		CTURBR(2)=RIGHTV(1)*CTURBR(1)/(10.0e-5*visc)!OMEGA INITIALIZATION
        END IF
	      END IF
 
	      IF (PASSIVESCALAR.GT.0)THEN
	      CTURBR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=PASS_INLET(INITCOND,POX,POY,POZ)*RIGHTV(1)
	      END IF
      END IF
      
      
      
      
      
      
      
      
    
    END IF
    
    
    CASE(2)!OUTFLOW SUBSONIC OR SUPERSONIC WILL BE CHOSEN BASED ON MACH NUMBER
     if (boundtype.eq.0)then
      rightv(1:nof_Variables)=leftv(1:nof_Variables)
      
     else

      if (Realgas.eq.1)then



                TEMPVECT1=LEFTV
                 TEMPVECT2=LEFTV

                call CONS2PRIM(N,TEMPVECT1,MP_PINFl,gammal)
                AGRT=SQRT((TEMPVECT1(5)+MP_PINFL)*GAMMAl/TEMPVECT1(1))
                U=TEMPVECT1(2)
                V=TEMPVECT1(3)
                w=TEMPVECT1(4)
                un  = U*nx + V*ny +w*nz


                call cons2div(N,TEMPVECT2,MP_PINFl,gammal)


                if (un.lt.0.0d0)then

                rightv=leftv


                else



                Mn  = un / AGRT

                if (Mn >= 1.0D0) then
                    ! Supersonic outflow: copy state
                    rightv(1:nof_Variables)=leftv(1:nof_Variables)
                else
                    ! Subsonic pressure outflow
                    rho_g   = TEMPVECT1(1)
                    !rho_s_g1(1:nof_SPECIES) = TEMPVECT1(dimensiona+4:nof_Variables)
                    Y_s_g(1:nof_SPECIES)   = TEMPVECT1(dimensiona+4:nof_Variables)

                  rg_rmix=zero
                    do rg_i=1,nof_species
                      rg_rmix=rg_rmix+Y_s_g(rg_i)/RG_MOLM(rg_i)
                    end do



              rg_rmix=RGS_Ru*rg_rmix




                    Rmix    = rg_rmix
                    T_g     = PRES / (rho_g * Rmix)

                    Tv_g    = TEMPVECT2(6)              ! simple extrapolation

                      ! Translational-rotational internal energy
                    rg_tr = 0.0D0
                        do rg_i = 1, nof_species
                          if (rg_i <= 3) then
                            RG_CVS(rg_i) = (5.0D0 / 2.0D0) *  RGS_Ru / RG_MOLM(rg_i)
                          else
                            RG_CVS(rg_i) = (3.0D0 / 2.0D0) *  RGS_Ru / RG_MOLM(rg_i)
                          end if
                          rg_tr = rg_tr + Y_s_g(rg_i) * RG_CVS(rg_i) * T_g
                        end do

                    ! Vibrational energy
                        RG_EV_TOTAL= 0.0D0
                        do rg_i = 1, 3
                          RG_EV_TOTAL = RG_EV_TOTAL + Y_s_g(rg_i)  * (RGS_Ru / RG_MOLM(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i)/Tv_g) - 1.0D0))
                        end do

                    RG_CHEM=zero

                        ! Chemical energy
                    do rg_i=1,nof_species
                            if (rg_hzero(RG_i).gt.1.0e-12)then
                            RG_CHEM=RG_CHEM-(Y_s_g(rg_i)*rg_hzero(RG_i)/RG_MOLM(rg_i))
                            end if
                    END DO

                      !RG_CHEM=zero  !set it to zero for testing
                    ! Kinetic energy


                    SKIN1=(oo2)*((U**2)+(V**2))

                      RIGHTV(1)=rho_g
                      RIGHTV(2)=rho_g*U
                      RIGHTV(3)=rho_g*V
                      RIGHTV(4)=rho_g*w
                      RIGHTV(5)=rho_g*(RG_EV_TOTAL+RG_TR+RG_CHEM+skin1)
                      RIGHTV(6)=rho_g*RG_EV_TOTAL


                      do rg_i=1,nof_species
                      RIGHTV(6+rg_i)=rho_g * Y_s_g(rg_i)
                      end do

                  end if
                  end if


                  else





     
     rightv(1:nof_Variables)=OUTFLOW(INITCOND,pox,poy,poz)
     CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
    
    SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
    SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
     
     SPS=SQRT((GAMMA*SUBSON2(5))/(SUBSON2(1)))
    VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2+SUBSON2(4)**2)
     
    SPS=SQRT((GAMMA*SUBSON2(5))/(SUBSON2(1)))
    
    CALL PRIM2CONS2(N,LEFTV,RIGHTV)
    
    IF (VEL/(SPS+TOLSMALL).GT.1.0D0)THEN	!SUPERSONIC
    CALL PRIM2CONS2(N,LEFTV,RIGHTV)
    rightv(1:nof_Variables)=leftv(1:nof_Variables)
    
    
      Else
    SUBSON3(5)=SUBSON1(5)
    SUBSON3(1)=SUBSON2(1)+(SUBSON3(5)-SUBSON2(5))/(SPS**2)
    SUBSON3(2)=SUBSON2(2)+(NX*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
    SUBSON3(3)=SUBSON2(3)+(NY*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
    SUBSON3(4)=SUBSON2(4)+(NZ*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
! 							
    rightv(1)=SUBSON3(1)
    rightv(2)=SUBSON3(2)*SUBSON3(1)
    rightv(3)=SUBSON3(3)*SUBSON3(1)
    rightv(4)=SUBSON3(4)*SUBSON3(1)
    SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2)+(SUBSON3(4)**2))
    IKINS=SUBSON3(5)/((GAMMA-1.0d0)*(SUBSON3(1)))
    rightv(5)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)
     
     
    end if
    end if
  end if
    
    
      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
	     
		  CTURBR(:)=CTURBL(:)
	
      END IF
    
    
    
    
    
    
    
    
    
    CASE(9)!OUTLETS SUBSONIC OR SUPERSONIC WILL BE CHOSEN BASED ON MACH NUMBER
    ! if (boundtype.eq.0)then
     ! rightv(1:nof_Variables)=leftv(1:nof_Variables)

     !else

     rightv(1:nof_Variables)=OUTFLOW2(INITCOND,pox,poy,poz)
     CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)

    SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)


    SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)

     SPS=SQRT((GAMMA*SUBSON2(5))/(SUBSON2(1)))
    VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2+SUBSON2(4)**2)

    SPS=SQRT((GAMMA*SUBSON2(5))/(SUBSON2(1)))

    CALL PRIM2CONS2(N,LEFTV,RIGHTV)

    IF (VEL/(SPS+TOLSMALL).GT.1.0D0)THEN	!SUPERSONIC
    CALL PRIM2CONS2(N,LEFTV,RIGHTV)
    rightv(1:nof_Variables)=leftv(1:nof_Variables)


      Else
    rightv(1:nof_Variables)=leftv(1:nof_Variables)
    SUBSON3(5)=subson1(5)
    SUBSON3(1)=SUBSON2(1)+(SUBSON3(5)-SUBSON2(5))/(SPS**2)
    SUBSON3(2)=SUBSON2(2)+(NX*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
    SUBSON3(3)=SUBSON2(3)+(NY*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
    SUBSON3(4)=SUBSON2(4)+(NZ*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
!
    rightv(1)=SUBSON3(1)
    rightv(2)=SUBSON3(2)*SUBSON3(1)
    rightv(3)=SUBSON3(3)*SUBSON3(1)
    rightv(4)=SUBSON3(4)*SUBSON3(1)
    SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2)+(SUBSON3(4)**2))
    IKINS=SUBSON3(5)/((GAMMA-1.0d0)*(SUBSON3(1)))
    rightv(5)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)



    !end if
    end if



      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN

		  CTURBR(:)=CTURBL(:)

      END IF
    
    

     CASE(99)    !BLEED BOUNDARY



     rightv(1:nof_Variables)=BLEED3D(Iconsidered,facex,pox,poy,poz)

      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN

		  CTURBR(:)=CTURBL(:)

      END IF



    
    CASE(3)!SYMMETRY
    
			      CALL ROTATEF(N,Cleft_ROT,leftV,ANGLE1,ANGLE2)
				IF (ILOCAL_RECON3(ICONSIDERED)%MRF.EQ.1)THEN
                    CRIGHT_ROT(1)=CLEFT_ROT(1)
                    CRIGHT_ROT(2)=-(CLEFT_ROT(2))+2.0D0*CLEFT_ROT(1)*SRF_SPEEDROT(2)
                    CRIGHT_ROT(3)=CLEFT_ROT(3)
                    CRIGHT_ROT(4)=CLEFT_ROT(4)
                    CRIGHT_ROT(5)=CLEFT_ROT(5)+2.0D0*CLEFT_ROT(1)*(SRF_SPEEDROT(2)**2)-2.0D0*CLEFT_ROT(2)*SRF_SPEEDROT(2)
                ELSE
                    CRIGHT_ROT(1)=CLEFT_ROT(1)
                    CRIGHT_ROT(2)=-CLEFT_ROT(2)
                    CRIGHT_ROT(3)=CLEFT_ROT(3)
                    CRIGHT_ROT(4)=CLEFT_ROT(4)
                    CRIGHT_ROT(5)=CLEFT_ROT(5)

                IF((MULTISPECIES.EQ.1).or.(realgas.eq.1))THEN
                      CRIGHT_ROT(6:nof_Variables)=CLEFT_ROT(6:nof_Variables)
                    END IF

				END IF
				     IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					    CTURBR(:)=CTURBL(:)

					  if (passivescalar.gt.0)then
					  cturbR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=&
						    ctURBL(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)

					  end if
				      END IF
				
			      CALL ROTATEb(N,rightv,Cright_ROT,ANGLE1,ANGLE2)
    
    
    
    CASE(4)!WALL
    
    
			      IF (ITESTCASE.EQ.3)THEN
			      
			       CALL ROTATEF(N,Cleft_ROT,leftV,ANGLE1,ANGLE2)
			      IF (ILOCAL_RECON3(ICONSIDERED)%MRF.EQ.1)THEN
                        CRIGHT_ROT(1)=CLEFT_ROT(1)
                        CRIGHT_ROT(2)=-(CLEFT_ROT(2))+2.0D0*CLEFT_ROT(1)*SRF_SPEEDROT(2)
                        CRIGHT_ROT(3)=CLEFT_ROT(3)
                        CRIGHT_ROT(4)=CLEFT_ROT(4)
                        CRIGHT_ROT(5)=CLEFT_ROT(5)+CLEFT_ROT(1)*(SRF_SPEEDROT(2)**2)*2.0D0-2.0D0*CLEFT_ROT(2)*SRF_SPEEDROT(2)
			      ELSE
         		      CRIGHT_ROT(:)=CLEFT_ROT(:)
			      CRIGHT_ROT(2)=-CLEFT_ROT(2)
                  END IF
				     IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					    CTURBR(:)=CTURBL(:)

					  if (passivescalar.gt.0)then
					  cturbR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=&
						    ctURBL(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)

					  end if
				      END IF
				
				
				
			      CALL ROTATEB(N,rightv,Cright_ROT,ANGLE1,ANGLE2)
			      
			      
			      
			      ELSE
                  IF (ILOCAL_RECON3(ICONSIDERED)%MRF.EQ.1)THEN
                    rightv(1)=leftv(1)
                    rightv(2)=-leftv(2)+2.0D0*leftv(1)*SRF_SPEED(2)
                    rightv(3)=-leftv(3)+2.0D0*leftv(1)*SRF_SPEED(3)
                    rightv(4)=-leftv(4)+2.0D0*leftv(1)*SRF_SPEED(4)
                    rightv(5)=leftv(5)+2.0D0*leftv(1)*(SRF_SPEED(2)**2+SRF_SPEED(3)**2+SRF_SPEED(4)**2)&
                                            -2.0D0*(leftv(2)*SRF_SPEED(2)+leftv(3)*SRF_SPEED(3)+leftv(4)*SRF_SPEED(4))
    
                  ELSE
                    rightv(1:nof_variables)=leftv(1:nof_variables)
                    rightv(2)=-leftv(2)
                    rightv(3)=-leftv(3)
                    rightv(4)=-leftv(4)










                  END IF
    
    
    
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					IF (TURBULENCEMODEL.NE.2)THEN
					    CTURBR(:)=-CTURBL(:)

					  if (passivescalar.gt.0)then
					  cturbR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=&
						    -ctURBL(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)

					  end if
					ELSE
					     CTURBR(1)=-CTURBL(1)
					     CTURBR(2)=60.0D0*VISC/(BETA_I1*(IELEM(N,ICONSIDERED)%WallDist**2))

					  if (passivescalar.gt.0)then
					  cturbR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=&
						    -ctURBL(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)

					  end if
					  
					
					
					
					
					END IF
				      END IF
    
    
				END IF
    
    
    
    
    
    CASE(6)!FARFIELD INFLOW OR OUTFLOW, SUBSONIC OR SUPERSONIC WILL BE CHOSEN BASED ON MACH NUMBER
	    CALL ROTATEF(N,Cleft_ROT,leftV,ANGLE1,ANGLE2)
	    vnb=cleft_rot(2)
	  
	    CALL CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
    
	    SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
	    SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
	    SPS=SQRT((GAMMA*SUBSON2(5))/(SUBSON2(1)))
	    VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2+SUBSON2(4)**2)
	  
	    CALL PRIM2CONS2(N,LEFTV,RIGHTV)

	  if (vnb.le.0.0d0)then		!inflow
			ibfc=-1
	
		  if ((abs(vnb)).ge.sps)then
				!supersonic
				rightv=INFLOW(INITCOND,POX,POY,POZ)
					
		  else
				!subsonic
			
			
			rightv=INFLOW(INITCOND,POX,POY,POZ)
	  	        
			  call CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
			  
			
			SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
			SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
			SPS=SQRT((GAMMA*SUBSON2(5))/(SUBSON2(1)))
			VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2+SUBSON2(4)**2)
	             CALL PRIM2CONS2(N,LEFTV,RIGHTV)
				    
		    SUBSON3(5)=0.5*((SUBSON1(5))+(SUBSON2(5))-(SUBSON2(1)*SPS*((NX*(SUBSON1(2)-SUBSON2(2)))+(NY*(SUBSON1(3)-SUBSON2(3)))&
	      +(NZ*(SUBSON1(4)-SUBSON2(4))))))
		    SUBSON3(1)=SUBSON1(1)+(SUBSON3(5)-SUBSON1(5))/(SPS**2)
		    SUBSON3(2)=SUBSON1(2)-(NX*(SUBSON1(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
		    SUBSON3(3)=SUBSON1(3)-(NY*(SUBSON1(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
		    SUBSON3(4)=SUBSON1(4)-(NZ*(SUBSON1(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
		    
		    
		      rightv(1)=SUBSON3(1)
    rightv(2)=SUBSON3(2)*SUBSON3(1)
    rightv(3)=SUBSON3(3)*SUBSON3(1)
    rightv(4)=SUBSON3(4)*SUBSON3(1)
    SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2)+(SUBSON3(4)**2))
    IKINS=SUBSON3(5)/((GAMMA-1.0d0)*(SUBSON3(1)))
    rightv(5)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)
		  
		  END IF
		  
		  
		  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
	      
      
      
      
	       IF (TURBULENCEMODEL.EQ.1)THEN
		  CTURBR(1)=VISC*TURBINIT*rightv(1)
	      END IF
	      IF (TURBULENCEMODEL.EQ.2)THEN	 
		CTURBR(1)=(1.5D0*I_turb_inlet*(ufreestream**2))*RIGHTV(1)!K INITIALIZATION
		CTURBR(2)=RIGHTV(1)*CTURBR(1)/(10.0e-5*visc)!OMEGA INITIALIZATION
		IF (ILOCAL_RECON3(ICONSIDERED)%MRF.EQ.1)THEN
            CTURBR(1)=(1.5D0*I_turb_inlet*(KINIT_SRF**2))*RIGHTV(1)!K INITIALIZATION
            CTURBR(2)=RIGHTV(1)*CTURBR(1)/(10.0e-5*visc)!OMEGA INITIALIZATION
		END IF
	      END IF
 
	      IF (PASSIVESCALAR.GT.0)THEN
	      CTURBR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=PASS_INLET(INITCOND,POX,POY,POZ)*RIGHTV(1)
	      END IF
		    END IF
		  
		  
      
	else
      
	  !outflow
	
	    ibfc=-2
		if ((abs(vnb)).ge.sps)then
		      
		      rightv(1:nof_Variables)=leftv(1:nof_Variables)
		
		else
		
		
		
		rightv(1:nof_Variables)=OUTFLOW(INITCOND,pox,poy,poz)
		
		call CONS2PRIM2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
		
    
    SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
    SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
     
     
     CALL PRIM2CONS2(N,LEFTV,RIGHTV)
   
    
    
    SUBSON3(5)=SUBSON1(5)
    SUBSON3(1)=SUBSON2(1)+(SUBSON3(5)-SUBSON2(5))/(SPS**2)
    SUBSON3(2)=SUBSON2(2)+(NX*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
    SUBSON3(3)=SUBSON2(3)+(NY*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
    SUBSON3(4)=SUBSON2(4)+(NZ*(SUBSON2(5)-SUBSON3(5)))/(SPS*SUBSON2(1))
! 							
    rightv(1)=SUBSON3(1)
    rightv(2)=SUBSON3(2)*SUBSON3(1)
    rightv(3)=SUBSON3(3)*SUBSON3(1)
    rightv(4)=SUBSON3(4)*SUBSON3(1)
    SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2)+(SUBSON3(4)**2))
    IKINS=SUBSON3(5)/((GAMMA-1.0d0)*(SUBSON3(1)))
    rightv(5)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)
	
	      
	      
	       IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
	     
		  CTURBR(:)=CTURBL(:)
	
		END IF
	      end if
	      

	END IF

	
! 	CALL ROTATEF(N,Cright_ROT,RIGHTV,ANGLE1,ANGLE2)
    
    

			      


END SELECT



END SUBROUTINE BOUNDARYS


SUBROUTINE BOUNDARYS2d(N,B_CODE,ICONSIDERED,facex,LEFTV,RIGHTV,POX,POY,POZ,ANGLE1,ANGLE2,NX,NY,NZ,CTURBL,CTURBR,CRIGHT_ROT,CLEFT_ROT,SRF_SPEED,SRF_SPEEDROT,IBFC)
!> @brief
!> This subroutine applies the boundary condition to each bounded cell
implicit none
integer,intent(in)::n,b_code,ICONSIDERED,facex
INTEGER,INTENT(INOUT)::IBFC
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::LEFTV,RIGHTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::SRF_SPEEDROT,SRF_SPEED
REAL,DIMENSION(1:dimensiona),INTENT(IN)::POX,POY,POZ
REAL,INTENT(IN)::ANGLE1,ANGLE2,NX,NY,NZ
REAL,DIMENSION(TURBULENCEEQUATIONS),INTENT(INOUT)::CTURBL,CTURBR
REAL,DIMENSION(1:NOF_VARIABLES+TURBULENCEEQUATIONS+PASSIVESCALAR),INTENT(INOUT)::CRIGHT_ROT,CLEFT_ROT
REAL,DIMENSION(1:NOF_VARIABLES)::SUBSON1,SUBSON2,SUBSON3,tempxv,TEMPVECT1,TEMPVECT2
REAl::MP_PINFL,MP_PINFR,GAMMAL,GAMMAR,u,v
REAL::SPS,SKINS,IKINS,VEL,vnb,theeta,reeta
REAL::INTENERGY,R1,U1,V1,W1,ET1,S1,IE1,P1,SKIN1,E1,RS,US,VS,WS,KHX,VHX,AMP,DVEL,rgg,tt1
REAL::MN,UN,AGRT,rho_g,rg_rmix,Rmix,T_g,Tv_g,rg_tr,RG_EV_TOTAL,RG_CHEM
REAL,DIMENSION(1:NOF_SPECIES)::rho_s_g1,Y_s_g,RG_CVS
integer::rg_i





SELECT CASE(B_CODE)

    
    CASE(1)!INFLOW SUBSONIC OR SUPERSONIC WILL BE CHOSEN BASED ON MACH NUMBER
    if (boundtype.eq.0)then	!SUPERSONIC
    
    RIGHTV(1:nof_Variables)=INFLOW2d(INITCOND,POX,POY)
    
    
    
        
    
    
    ELSE		!SUBSONIC
    RIGHTV(1:nof_Variables)=INFLOW2d(INITCOND,POX,POY)

    CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
    
    SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
    SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
    SPS=SQRT((GAMMA*SUBSON2(4))/(SUBSON2(1)))
    VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2)
    

    CALL PRIM2CONS2(N,LEFTV,RIGHTV)
    

      IF (VEL/(SPS+TOLSMALL).GT.1.0D0)THEN	!SUPERSONIC
      
      
      RIGHTV(1:nof_Variables)=INFLOW2d(INITCOND,POX,POY)
      
      
      
      ELSE		!SUBSONIC
      
      
      SUBSON3(4)=0.5*((SUBSON1(4))+(SUBSON2(4))-(SUBSON2(1)*SPS*((NX*(SUBSON1(2)-SUBSON2(2)))+(NY*(SUBSON1(3)-SUBSON2(3)))&
)))
      SUBSON3(1)=SUBSON1(1)+(SUBSON3(4)-SUBSON1(4))/(SPS**2)
      SUBSON3(2)=SUBSON1(2)-(NX*(SUBSON1(4)-SUBSON3(4)))/(SPS*SUBSON2(1))
      SUBSON3(3)=SUBSON1(3)-(NY*(SUBSON1(4)-SUBSON3(4)))/(SPS*SUBSON2(1))
      
      
      		    
      
      
    rightv(1)=SUBSON3(1)
    rightv(2)=SUBSON3(2)*SUBSON3(1)
    rightv(3)=SUBSON3(3)*SUBSON3(1)
    
    SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2))
    IKINS=SUBSON3(4)/((GAMMA-1.0d0)*(SUBSON3(1)))
    rightv(4)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)
      
      
       
      
      END IF

      
       END IF
      
      
      
      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
	      
      
      
      
	      IF (TURBULENCEMODEL.EQ.1)THEN
		  CTURBR(1)=VISC*TURBINIT
	      END IF
	     IF (TURBULENCEMODEL.EQ.2)THEN	 
		CTURBR(1)=(1.5D0*I_turb_inlet*(ufreestream**2))*RIGHTV(1)!K INITIALIZATION
		CTURBR(2)=RIGHTV(1)*CTURBR(1)/(10.0e-5*visc)!OMEGA INITIALIZATION
		
		
	      END IF
 
	      IF (PASSIVESCALAR.GT.0)THEN
	      CTURBR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=PASS_INLET2d(INITCOND,POX,POY)*RIGHTV(1)
	      END IF
      END IF
      
      
      
      
      
      
      
      
    
   
    
    
    CASE(2)!OUTFLOW SUBSONIC OR SUPERSONIC WILL BE CHOSEN BASED ON MACH NUMBER
     if (boundtype.eq.0)then
      rightv(1:nof_Variables)=leftv(1:nof_Variables)
      
     else
     




                if (Realgas.eq.1)then



                TEMPVECT1=LEFTV
                 TEMPVECT2=LEFTV

                call CONS2PRIM(N,TEMPVECT1,MP_PINFl,gammal)
                AGRT=SQRT((TEMPVECT1(4)+MP_PINFL)*GAMMAl/TEMPVECT1(1))
                U=TEMPVECT1(2)
                V=TEMPVECT1(3)
                un  = U*nx + V*ny


                call cons2div(N,TEMPVECT2,MP_PINFl,gammal)


                if (un.lt.0.0d0)then

                rightv=leftv


                else



                Mn  = un / AGRT

                if (Mn >= 1.0D0) then
                    ! Supersonic outflow: copy state
                    rightv(1:nof_Variables)=leftv(1:nof_Variables)
                else
                    ! Subsonic pressure outflow
                    rho_g   = TEMPVECT1(1)
                    !rho_s_g1(1:nof_SPECIES) = TEMPVECT1(dimensiona+4:nof_Variables)
                    Y_s_g(1:nof_SPECIES)   = TEMPVECT1(dimensiona+4:nof_Variables)

                  rg_rmix=zero
                    do rg_i=1,nof_species
                      rg_rmix=rg_rmix+Y_s_g(rg_i)/RG_MOLM(rg_i)
                    end do



              rg_rmix=RGS_Ru*rg_rmix




                    Rmix    = rg_rmix
                    T_g     = PRES / (rho_g * Rmix)

                    Tv_g    = TEMPVECT2(5)              ! simple extrapolation

                      ! Translational-rotational internal energy
                    rg_tr = 0.0D0
                        do rg_i = 1, nof_species
                          if (rg_i <= 3) then
                            RG_CVS(rg_i) = (5.0D0 / 2.0D0) *  RGS_Ru / RG_MOLM(rg_i)
                          else
                            RG_CVS(rg_i) = (3.0D0 / 2.0D0) *  RGS_Ru / RG_MOLM(rg_i)
                          end if
                          rg_tr = rg_tr + Y_s_g(rg_i) * RG_CVS(rg_i) * T_g
                        end do

                    ! Vibrational energy
                        RG_EV_TOTAL= 0.0D0
                        do rg_i = 1, 3
                          RG_EV_TOTAL = RG_EV_TOTAL + Y_s_g(rg_i)  * (RGS_Ru / RG_MOLM(rg_i)) * (rg_thetag(rg_i) / (exp(rg_thetag(rg_i)/Tv_g) - 1.0D0))
                        end do

                    RG_CHEM=zero

                        ! Chemical energy
                    do rg_i=1,nof_species
                            if (rg_hzero(RG_i).gt.1.0e-12)then
                            RG_CHEM=RG_CHEM-(Y_s_g(rg_i)*rg_hzero(RG_i)/RG_MOLM(rg_i))
                            end if
                    END DO

                      !RG_CHEM=zero  !set it to zero for testing
                    ! Kinetic energy


                    SKIN1=(oo2)*((U**2)+(V**2))

                      RIGHTV(1)=rho_g
                      RIGHTV(2)=rho_g*U
                      RIGHTV(3)=rho_g*V
                      RIGHTV(4)=rho_g*(RG_EV_TOTAL+RG_TR+RG_CHEM+skin1)
                      RIGHTV(5)=rho_g*RG_EV_TOTAL


                      do rg_i=1,nof_species
                      RIGHTV(5+rg_i)=rho_g * Y_s_g(rg_i)
                      end do

                  end if
                  end if

                  else
                        rightv(1:nof_Variables)=OUTFLOW2d(INITCOND,pox,poy)
                        CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)

                        SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
                        SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)

                        SPS=SQRT((GAMMA*SUBSON2(4))/(SUBSON2(1)))
                        VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2)

                        SPS=SQRT((GAMMA*SUBSON2(4))/(SUBSON2(1)))

                        CALL PRIM2CONS2(N,LEFTV,RIGHTV)

                        IF (VEL/(SPS+TOLSMALL).GT.1.0D0)THEN	!SUPERSONIC
                        rightv(1:nof_Variables)=leftv(1:nof_Variables)


                          Else
                        SUBSON3(4)=SUBSON1(4)
                        SUBSON3(1)=SUBSON2(1)+(SUBSON3(4)-SUBSON2(4))/(SPS**2)
                        SUBSON3(2)=SUBSON2(2)+(NX*(SUBSON2(4)-SUBSON3(4)))/(SPS*SUBSON2(1))
                        SUBSON3(3)=SUBSON2(3)+(NY*(SUBSON2(4)-SUBSON3(4)))/(SPS*SUBSON2(1))

                    !
                        rightv(1)=SUBSON3(1)
                        rightv(2)=SUBSON3(2)*SUBSON3(1)
                        rightv(3)=SUBSON3(3)*SUBSON3(1)

                        SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2))
                        IKINS=SUBSON3(4)/((GAMMA-1.0d0)*(SUBSON3(1)))
                        rightv(4)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)


                        end if

                  end if
    end if

    
    
      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
	     
		  CTURBR(:)=CTURBL(:)
	
      END IF
    
    
    
    
    
    
     CASE(99)  !BLEED


     rightv(1:nof_Variables)=BLEED2d(Iconsidered,facex,pox,poy)



      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN

		  CTURBR(:)=CTURBL(:)

      END IF
    
    
    
    
    
    
    CASE(3)!SYMMETRY
    
    
    
                            IF ((INITCOND.EQ.102).or.(INITCOND.EQ.30).or.(initcond.eq.222))THEN	!shock density interaction
                                   IF ((INITCOND.EQ.102))THEN	!shock density interaction
                            if (pox(1).lt.((1.0d0/6.0d0)+((1.0d0+20.0d0*T)/(sqrt(3.0d0)))))then
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
                            rightv(1)=R1
                            rightv(2)=R1*U1
                            rightv(3)=R1*V1
                            rightv(4)=E1
                            end if
    
                            
                            
                            
                            
                            
                             IF ((INITCOND.EQ.222))THEN
                             if (sqrt((pox(1)**2)+(poy(1)**2)).lt.(T/3.0d0))then
                            r1=16.0d0
                            u1=0.0
                            v1=0.0
                            p1=16.0d0/3.0d0
                            else
                            r1=1.0d0+(t/sqrt((pox(1)**2)+(poy(1)**2)))
                            reeta=-1
                            THEETA=ATAN(POY(1)/POX(1))
                            U1=REETA*COS(THEETA)
                            V1=REETA*SIN(THEETA)
                            P1=1.0E-6
                            end if
                            SKIN1=(OO2)*((U1**2)+(V1**2))
                            !INTERNAL ENERGY 
                            IE1=((P1)/((GAMMA-1.0D0)*R1))
                            !TOTAL ENERGY
                            E1=(P1/(GAMMA-1))+(R1*SKIN1)
                            !VECTOR OF CONSERVED VARIABLES NOW
                            rightv(1)=R1
                            rightv(2)=R1*U1
                            rightv(3)=R1*V1
                            rightv(4)=E1
                             
                             
                             
                             end if
			     
			       IF ((INITCOND.EQ.30))THEN
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
                           ! rightv(1)=leftv(1)
                            !rightv(2)=0.0
                            !rightv(3)=0.0
                            !rightv(4)=leftv(4)
                            
                             rightv(1)=R1
                            rightv(2)=R1*U1
                            rightv(3)=R1*V1
                            rightv(4)=E1
			      
			      end if
			      
			      
			      
			      
			      else
			      
			      
			       CALL ROTATEF2d(N,Cleft_ROT,leftV,ANGLE1,ANGLE2)
			      
                  CRIGHT_ROT(1)=CLEFT_ROT(1)
			      CRIGHT_ROT(2)=-CLEFT_ROT(2)
			      CRIGHT_ROT(3)=CLEFT_ROT(3)
			      CRIGHT_ROT(4)=CLEFT_ROT(4)
			      
			       IF((MULTISPECIES.EQ.1).or.(realgas.eq.1))THEN
                      CRIGHT_ROT(5:nof_Variables)=CLEFT_ROT(5:nof_Variables)
                    
                    END IF
			     
					 
				     IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					    CTURBR(:)=CTURBL(:)

					  if (passivescalar.gt.0)then
					  cturbR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=&
						    ctURBL(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)

					  end if
				      END IF
				
			    	
				
			      CALL ROTATEb2d(N,rightv,Cright_ROT,ANGLE1,ANGLE2)
                            end if
    
    
    CASE(4)!WALL
    
			     IF (ITESTCASE.EQ.3)THEN
			      
			       CALL ROTATEF2D(N,Cleft_ROT,leftV,ANGLE1,ANGLE2)
			      
			      
                      IF ((multispecies.EQ.1).or.(realgas.eq.1))then
                          CRIGHT_ROT(:)=CLEFT_ROT(:)
                      CRIGHT_ROT(2)=-CLEFT_ROT(2)


                      else
                    CRIGHT_ROT(:)=CLEFT_ROT(:)
                          CRIGHT_ROT(1)=CLEFT_ROT(1)
                      CRIGHT_ROT(2)=-CLEFT_ROT(2)
                      CRIGHT_ROT(3)=CLEFT_ROT(3)
                      CRIGHT_ROT(4)=CLEFT_ROT(4)
                      end if
			     
			      
					 
				     IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					    CTURBR(:)=CTURBL(:)

					  if (passivescalar.gt.0)then
					  cturbR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=&
						    ctURBL(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)

					  end if
				      END IF
				
				
				
			      CALL ROTATEb2D(N,rightv,Cright_ROT,ANGLE1,ANGLE2)
			      
			    
			      
			      ELSE


                  rightv(1:nof_Variables)=leftv(1:nof_Variables)


			      rightv(2)=-leftv(2)
			      rightv(3)=-leftv(3)
			      
! 			      rightv(4)=leftv(4)
    

    
    
				      IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
					IF (TURBULENCEMODEL.NE.2)THEN
					    CTURBR(:)=-CTURBL(:)

					  if (passivescalar.gt.0)then
					  cturbR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=&
						    -ctURBL(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)

					  end if
					ELSE
					     CTURBR(1)=-CTURBL(1)
					     CTURBR(2)=60.0D0*VISC/(BETA_I1*(IELEM(N,ICONSIDERED)%WallDist**2))

					  if (passivescalar.gt.0)then
					  cturbR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=&
						    -ctURBL(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)

					  end if
					  
					
					
					
					
					END IF
				      END IF
    
    
				END IF
    
    
    
    
    
    
    
    
    
    CASE(6)!FARFIELD INFLOW OR OUTFLOW, SUBSONIC OR SUPERSONIC WILL BE CHOSEN BASED ON MACH NUMBER
	 CALL ROTATEF2d(N,Cleft_ROT,leftV,ANGLE1,ANGLE2)
	    vnb=cleft_rot(2)/CLEFT_ROT(1)
	    
	    
	  
	    
	    CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
    
	    SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
	    SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
	    SPS=SQRT((GAMMA*SUBSON2(4))/(SUBSON2(1)))
	    VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2)
	  
	  CALL PRIM2CONS2(N,LEFTV,RIGHTV)

	  if (vnb.le.0.0d0)then		!inflow
			ibfc=-1
	
		  if ((abs(vnb)).ge.sps)then
				!supersonic
				rightv(1:nof_Variables)=INFLOW2d(INITCOND,POX,POY)
					
		  else
				!subsonic
				
			rightv(1:nof_Variables)=INFLOW2d(INITCOND,POX,POY)
	  	  
			  
			  CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
			
			SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
			SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
			SPS=SQRT((GAMMA*SUBSON2(4))/(SUBSON2(1)))
			VEL=sqrt(SUBSON2(2)**2+SUBSON2(3)**2)
			  CALL PRIM2CONS2(N,LEFTV,RIGHTV)
				    
		    SUBSON3(4)=0.5d0*((SUBSON1(4))+(SUBSON2(4))-(SUBSON2(1)*SPS*((NX*(SUBSON1(2)-SUBSON2(2)))+(NY*(SUBSON1(3)-SUBSON2(3))))))
		    SUBSON3(1)=SUBSON1(1)+(SUBSON3(4)-SUBSON1(4))/(SPS**2)
		    SUBSON3(2)=SUBSON1(2)-(NX*(SUBSON1(4)-SUBSON3(4)))/(SPS*SUBSON2(1))
		    SUBSON3(3)=SUBSON1(3)-(NY*(SUBSON1(4)-SUBSON3(4)))/(SPS*SUBSON2(1))
		    
		    

		    rightv(1)=SUBSON3(1)
    rightv(2)=SUBSON3(2)*SUBSON3(1)
    rightv(3)=SUBSON3(3)*SUBSON3(1)
    
    SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2))
    IKINS=SUBSON3(4)/((GAMMA-1.0d0)*(SUBSON3(1)))
    rightv(4)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)
		  
		  END IF
		  
		  
		  IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
	      
      
      
      
	        IF (TURBULENCEMODEL.EQ.1)THEN
		  CTURBR(1)=VISC*TURBINIT*rightv(1)
	      END IF
	     IF (TURBULENCEMODEL.EQ.2)THEN	 
		CTURBR(1)=(1.5D0*I_turb_inlet*(ufreestream**2))*RIGHTV(1)!K INITIALIZATION
		CTURBR(2)=RIGHTV(1)*CTURBR(1)/(10.0e-5*visc)!OMEGA INITIALIZATION
	      END IF
 
	      IF (PASSIVESCALAR.GT.0)THEN
	      CTURBR(TURBULENCEEQUATIONS+1:TURBULENCEEQUATIONS+PASSIVESCALAR)=PASS_INLET2d(INITCOND,POX,POY)*RIGHTV(1)
	      END IF
		    END IF
		  
		  
      
	else
      
	  !outflow
	
	    ibfc=-2
		if ((abs(vnb)).ge.sps)then
		
		      rightv(1:nof_Variables)=leftv(1:nof_Variables)
		
		else
		
		
		rightv(1:nof_Variables)=OUTFLOW2d(INITCOND,pox,poy)
		 CALL cons2prim2(N,LEFTV,RIGHTV,MP_PINFL,MP_PINFR,GAMMAL,GAMMAR)
    
    SUBSON1(1:nof_Variables)=RIGHTV(1:nof_Variables)
    SUBSON2(1:nof_Variables)=LEFTV(1:nof_Variables)
     
     CALL PRIM2CONS2(N,LEFTV,RIGHTV)
     
   
    
    
    SUBSON3(4)=SUBSON1(4)
    SUBSON3(1)=SUBSON2(1)+(SUBSON3(4)-SUBSON2(4))/(SPS**2)
    SUBSON3(2)=SUBSON2(2)+(NX*(SUBSON2(4)-SUBSON3(4)))/(SPS*SUBSON2(1))
    SUBSON3(3)=SUBSON2(3)+(NY*(SUBSON2(4)-SUBSON3(4)))/(SPS*SUBSON2(1))
    
! 							
    rightv(1)=SUBSON3(1)
    rightv(2)=SUBSON3(2)*SUBSON3(1)
    rightv(3)=SUBSON3(3)*SUBSON3(1)
    
    SKINS=oo2*((SUBSON3(2)**2)+(SUBSON3(3)**2))
    IKINS=SUBSON3(4)/((GAMMA-1.0d0)*(SUBSON3(1)))
    rightv(4)=(SUBSON3(1)*(IKINS))+(SUBSON3(1)*SKINS)
	
	      
	      
	       IF ((TURBULENCE.EQ.1).OR.(PASSIVESCALAR.GT.0))THEN
	     
		  CTURBR(:)=CTURBL(:)
	
		END IF
	        end if
	      

	END IF

	

    
    

			      


END SELECT



END SUBROUTINE BOUNDARYS2d


SUBROUTINE COMPUTE_EIGENVECTORS(N,RVEIGL,RVEIGR,EIGVL,EIGVR,GAMMA)
!> @brief
!> This subroutine computes the left and right eigenvectors 
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::RVEIGL,RVEIGR
REAL,INTENT(IN)::GAMMA
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES),INTENT(INOUT)::EIGVL,EIGVR
REAL::RS,US,VS,WS,ES,PS,VVS,AS,HS,GAMMAM1,vsd,OOR1,OOR2
INTEGER::IVGT

EIGVR=ZERO
GAMMAM1=GAMMA-1.0D0
OOR1=1.0D0/RVEIGL(1)
OOR2=1.0D0/RVEIGR(1)

RS=OO2*(RVEIGL(1)+RVEIGR(1))
US=OO2*((RVEIGL(2)*OOR1)+(RVEIGR(2)*OOR2))
VS=OO2*((RVEIGL(3)*OOR1)+(RVEIGR(3)*OOR2))
WS=OO2*((RVEIGL(4)*OOR1)+(RVEIGR(4)*OOR2))
ES=OO2*(RVEIGL(5)+RVEIGR(5))
 
VVS=(US**2)+(VS**2)+(WS**2)
VSD=OO2*VVS
PS=(GAMMA-1.0D0)*(ES - OO2*RS*VVS)
AS=SQRT(GAMMA*PS/RS)
HS=(OO2*VVS) + ((AS**2)/GAMMAM1)
EIGVR(1,1)=1.0D0		; EIGVR(1,2)=1.0D0	; EIGVR(1,3)=0.0D0	; EIGVR(1,4)=0.0D0	; EIGVR(1,5)=1.0D0
EIGVR(2,1)=US-AS	; EIGVR(2,2)=US		; EIGVR(2,3)=0.0D0	; EIGVR(2,4)=0.0D0	; EIGVR(2,5)=US+AS
EIGVR(3,1)=VS		; EIGVR(3,2)=VS		; EIGVR(3,3)=1.0D0	; EIGVR(3,4)=0.0D0	; EIGVR(3,5)=VS
EIGVR(4,1)=WS		; EIGVR(4,2)=WS		; EIGVR(4,3)=0.0D0	; EIGVR(4,4)=1.0D0	; EIGVR(4,5)=WS
EIGVR(5,1)=HS-(US*AS)	; EIGVR(5,2)=OO2*VVS	; EIGVR(5,3)=VS		; EIGVR(5,4)=WS		; EIGVR(5,5)=HS+(US*AS)

EIGVL(1,1)=(HS*us + as*us**2 + as*vs**2 - as*VSD - us*VSD + as*ws**2)/ (2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(1,2)=(-HS - as*us + VSD)/(2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(1,3)= -((as*vs)/(2.0D0*as*HS - 2.0D0*as*VSD))
EIGVL(1,4)=  -((as*ws)/(2.0D0*as*HS - 2.0D0*as*VSD))
EIGVL(1,5)=as/(2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(2,1)=(2.0*as*HS - 2.0*as*us**2 - 2.0D0*as*vs**2 -2.0D0*as*ws**2)/(2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(2,2)=(2.0D0*as*us)/(2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(2,3)=(2.0D0*as*vs)/(2.0D0*as*HS - 2.0D0*as*VSD)
 EIGVL(2,4)=(2.0D0*as*ws)/(2.0D0*as*HS - 2.0D0*as*VSD) 
 EIGVL(2,5)=(-2.0D0*as)/(2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(3,1)=(-2.0*as*HS*vs + 2.0*as*vs*VSD)/(2.0*as*HS - 2.0*as*VSD)
EIGVL(3,2)=0.0D0
EIGVL(3,3)=1.0D0
EIGVL(3,4)=0.0D0
 EIGVL(3,5)=0.0D0
EIGVL(4,1)=(-2.0D0*as*HS*ws + 2.0D0*as*VSD*ws)/(2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(4,2)=0.0D0
 EIGVL(4,3)=0.0D0
EIGVL(4,4)=1.0D0
EIGVL(4,5)=0.0D0
EIGVL(5,1)=(-(HS*us) + as*us**2 + as*vs**2 - as*VSD + us*VSD + as*ws**2)/ (2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(5,2)=(HS - as*us - VSD)/(2.0D0*as*HS - 2.0D0*as*VSD)
EIGVL(5,3)=-((as*vs)/(2.0D0*as*HS - 2.0D0*as*VSD))
 EIGVL(5,4)=-((as*ws)/(2.0D0*as*HS - 2.0D0*as*VSD))
EIGVL(5,5)= as/(2.0D0*as*HS - 2.0D0*as*VSD)
 
 

END SUBROUTINE COMPUTE_EIGENVECTORS








SUBROUTINE COMPUTE_JACOBIANSE(N,ICONSIDERED,EIGVL,RVEIGL,GAMMA,ANGLE1,ANGLE2,SRF_SPEEDROT,nx,ny,nz)
!> @brief
!> This subroutine computes the Jacobians for the implicit time stepping
IMPLICIT NONE
INTEGER,INTENT(IN)::N,ICONSIDERED
REAL,INTENT(IN)::ANGLE1,ANGLE2
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::RVEIGL,SRF_SPEEDROT
REAL,INTENT(IN)::GAMMA
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES),INTENT(INOUT)::EIGVL
REAL::RS,US,VS,WS,ES,PS,VVS,AS,HS,GAMMAM1,vsd,PHI,A1,A2,A3,OORS,NX,NY,NZ
INTEGER::IVGT,i


A2=GAMMA-1.0D0
A3=GAMMA-2.0D0
OORS=1.0D0/RVEIGL(1)
RS=(RVEIGL(1))
US=(RVEIGL(2)*OORS)
VS=(RVEIGL(3)*OORS)
WS=(RVEIGL(4)*OORS)
ES=(RVEIGL(5)*OORS)
PHI=OO2*(A2)*((US*US)+(VS*VS)+(WS*WS))
 A1=GAMMA*ES-PHI

 
 
VVS=NX*US+NY*VS+NZ*WS

IF (ILOCAL_RECON3(ICONSIDERED)%MRF.EQ.1)THEN
    EIGVL(1,1)=0.0D0-SRF_SPEEDROT(2);EIGVL(1,2)=NX	; 		EIGVL(1,3)=NY	; 		EIGVL(1,4)=NZ	; 		EIGVL(1,5)=0.0D0
    EIGVL(2,1)=NX*PHI-US*VVS	;EIGVL(2,2)=VVS-A3*NX*US-SRF_SPEEDROT(2)	;EIGVL(2,3)=NY*US-A2*NX*VS	; EIGVL(2,4)=NZ*US-A2*NX*WS; EIGVL(2,5)=A2*NX
    EIGVL(3,1)=NY*PHI-VS*VVS	;EIGVL(3,2)=NX*VS-A2*NY*US	; EIGVL(3,3)=VVS-A3*NY*VS-SRF_SPEEDROT(2)	; EIGVL(3,4)=NZ*VS-A2*NY*WS; EIGVL(3,5)=A2*NY
    EIGVL(4,1)=NZ*PHI-WS*VVS	;EIGVL(4,2)=NX*WS-A2*NZ*US	; EIGVL(4,3)=NY*WS-A2*NZ*VS	; EIGVL(4,4)=VVS-A3*NZ*WS-SRF_SPEEDROT(2); EIGVL(4,5)=A2*NZ
    EIGVL(5,1)=VVS*(PHI-A1)         ; EIGVL(5,2)=NX*A1-A2*US*VVS	; EIGVL(5,3)=NY*A1-A3*VS*VVS; EIGVL(5,4)=NZ*A1-A2*WS*VVS; EIGVL(5,5)=GAMMA*VVS-SRF_SPEEDROT(2)
ELSE
    EIGVL(1,1)=0.0D0		; 		EIGVL(1,2)=NX	; 		EIGVL(1,3)=NY	; 		EIGVL(1,4)=NZ	; 		EIGVL(1,5)=0.0D0
    EIGVL(2,1)=NX*PHI-US*VVS	; 	EIGVL(2,2)=VVS-A3*NX*US	; 	EIGVL(2,3)=NY*US-A2*NX*VS	; EIGVL(2,4)=NZ*US-A2*NX*WS	; EIGVL(2,5)=A2*NX
    EIGVL(3,1)=NY*PHI-VS*VVS		; EIGVL(3,2)=NX*VS-A2*NY*US	; EIGVL(3,3)=VVS-A3*NY*VS	; EIGVL(3,4)=NZ*VS-A2*NY*WS	; EIGVL(3,5)=A2*NY
    EIGVL(4,1)=NZ*PHI-WS*VVS		; EIGVL(4,2)=NX*WS-A2*NZ*US	; EIGVL(4,3)=NY*WS-A2*NZ*VS	; EIGVL(4,4)=VVS-A3*NZ*WS	; EIGVL(4,5)=A2*NZ
    EIGVL(5,1)=VVS*(PHI-A1)	;		 EIGVL(5,2)=NX*A1-A2*US*VVS	; EIGVL(5,3)=NY*A1-A3*VS*VVS	; EIGVL(5,4)=NZ*A1-A2*WS*VVS	; EIGVL(5,5)=GAMMA*VVS
END IF

 
 if (realgas.eq.1)then
 ! --- Extra equations: 6 = vibrational energy, 7..11 = species densities ---
! Approximate scalar advection: only diagonal = VVS

  ! zero all extra rows first
  DO i = 6, NOF_VARIABLES
    EIGVL(i,1:NOF_VARIABLES) = 0.0D0
  END DO

  ! vibrational energy equation
  EIGVL(6,6) = VVS

  ! species equations (5 species -> rows 7..11)
  DO i = 7, 11
    EIGVL(i,i) = VVS
  END DO


 end if




 

END SUBROUTINE COMPUTE_JACOBIANSE




SUBROUTINE COMPUTE_EIGENVECTORS2D(N,RVEIGL,RVEIGR,EIGVL,EIGVR,GAMMA)
!> @brief
!> This subroutine computes the left and right eigenvectors  in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::RVEIGL,RVEIGR
REAL,INTENT(IN)::GAMMA
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES),INTENT(INOUT)::EIGVL,EIGVR
REAL::RS,US,VS,WS,ES,PS,VVS,AS,HS,GAMMAM1,vsd,G8,s1,s2,vtots,OOR1,OOR2
INTEGER::IVGT,J,K

EIGVR=ZERO
GAMMAM1=GAMMA-1.0D0
OOR1=1.0D0/RVEIGL(1)
OOR2=1.0D0/RVEIGR(1)


EIGVR=0.0D0
GAMMAM1=GAMMA-1.0D0
RS=0.5*(RVEIGL(1)+RVEIGR(1))
US=0.5*((RVEIGL(2)*OOR1)+(RVEIGR(2)*OOR2))
VS=0.5*((RVEIGL(3)*OOR1)+(RVEIGR(3)*OOR2))
ES=0.5*(RVEIGL(4)+RVEIGR(4))
G8 = gamma - 1.0D0
 
VTOTS = US**2 + VS**2
   PS = (GAMMA-1)*(ES - OO2*RS*VTOTS)
   AS = SQRT(GAMMA*PS/RS)
   HS = OO2*VTOTS + (AS**2)/(G8)
  

   EIGVR(1,1) =1.D0;        EIGVR(1,2) = 1.D0;                EIGVR(1,3) =0.D0;  EIGVR(1,4) = 1.D0  
   EIGVR(2,1) =Us-As;     EIGVR(2,2) = Us;                EIGVR(2,3) =0.D0;  EIGVR(2,4) = Us+As 
   EIGVR(3,1) =Vs;        EIGVR(3,2) = Vs;                EIGVR(3,3) =1.D0;  EIGVR(3,4) = Vs   
   EIGVR(4,1) =Hs-Us*As;  EIGVR(4,2) = OO2*VTOTS; EIGVR(4,3) =Vs ; EIGVR(4,4) = Hs+Us*As 

   S1 = AS/(Gamma-1D0)
   S2 = AS**2/(Gamma-1D0)

   EIGVl(1,1) =   HS + S1*(US-AS);   EIGVl(1,2) = -(US+S1); EIGVl(1,3) = -VS;            EIGVl(1,4) = 1.D0
   EIGVl(2,1) =-2D0*HS+4D0*S2;           EIGVl(2,2) = 2D0*US;     EIGVl(2,3) = 2D0*VS;           EIGVl(2,4) = -2.D0
   EIGVl(3,1) =-2D0*VS*S2;             EIGVl(3,2) =   0.D0;     EIGVl(3,3) = 2D0*S2;           EIGVl(3,4) = 0.D0
   EIGVl(4,1) = HS- S1*(US+AS);      EIGVl(4,2) = -US+S1;   EIGVl(4,3) = - VS;           EIGVl(4,4) = 1.D0 

DO J=1,4
	DO K=1,4

   EIGVl(J,K) = EIGVl(J,K)/(2D0*S2)
	END DO
END DO
 
 

END SUBROUTINE COMPUTE_EIGENVECTORS2D


 subroutine COMPUTE_JACOBIANSE2D(N,EIGVL,rveigl,GAMMA,ANGLE1,ANGLE2,nx,ny,nz)
 !> @brief
!> This subroutine computes the Jacobians for the implicit time stepping in 2D
IMPLICIT NONE
INTEGER,INTENT(IN)::N
REAL,INTENT(IN)::ANGLE1,ANGLE2,nx,ny,nz
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::RVEIGL
REAL,INTENT(IN)::GAMMA
REAL,DIMENSION(1:NOF_VARIABLES,1:NOF_VARIABLES),INTENT(INOUT)::EIGVL
REAL::RS,US,VS,ES,PS,VVS,AS,HS,GAMMAM1,vsd,PHI,A1,A2,A3,OORS
INTEGER::IVGT


A2=GAMMA-1.0D0
A3=GAMMA-2.0D0
OORS=1.0D0/RVEIGL(1)
RS=(RVEIGL(1))
US=(RVEIGL(2)*OORS)
VS=(RVEIGL(3)*OORS)
ES=(RVEIGL(4)*OORS)
PHI=OO2*(A2)*((US*US)+(VS*VS))
 A1=GAMMA*ES-PHI

 
VVS=NX*US+NY*VS


EIGVL(1,1)=0.0D0		; 		EIGVL(1,2)=NX	; 		EIGVL(1,3)=NY	; 		EIGVL(1,4)=0.0D0
EIGVL(2,1)=NX*PHI-US*VVS	; 	EIGVL(2,2)=VVS-A3*NX*US	; 	EIGVL(2,3)=NY*US-A2*NX*VS	; EIGVL(2,4)=A2*NX
EIGVL(3,1)=NY*PHI-VS*VVS		; EIGVL(3,2)=NX*VS-A2*NY*US	; EIGVL(3,3)=VVS-A3*NY*VS	; EIGVL(3,4)=A2*NY

EIGVL(4,1)=VVS*(PHI-A1)	;		 EIGVL(4,2)=NX*A1-A2*US*VVS	; EIGVL(4,3)=NY*A1-A2*VS*VVS	; EIGVL(4,4)=GAMMA*VVS


 
 

END SUBROUTINE COMPUTE_JACOBIANSE2D








SUBROUTINE EDDYVISCO(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
!> @brief
!> This subroutine computes the tubulent eddy viscosity for turbulence models
	IMPLICIT NONE
    REAL,DIMENSION(1:2),INTENT(INOUT)::TURBMV
    real,dimension(1:nof_Variables),INTENT(IN)::leftv,rightv
    REAL,DIMENSION(1),INTENT(INOUT)::ETVM
    REAL,DIMENSION(1:4),INTENT(INOUT)::VISCL,LAML
    REAL,DIMENSION(1:20),INTENT(INOUT)::EDDYFL,EDDYFR
	INTEGER,INTENT(IN)::N
	REAL::ML,SMMM,MTT,chi,tolepsma,chipow3,fv1
	INTEGER:: IHGT, IHGJ
	REAL,DIMENSION(3,3)::VORTET,TVORT,SVORT
	REAL:: ux,uy,uz,vx,vy,vz,wx,wy,wz
	REAL:: wally, D_omplus, phi_2, phi_1, F_1, F_2, k_0, om_0,alpha_inf, alpha_star, Mu_turb
	REAL:: sigma_k_L, sigma_om_L, sigma_k_r, sigma_om_r
	REAL::SNORM,dervk_dervom,Re_t_SST,RHO_0,beta_i




	tolepsma = 10e-16
	
	IF (TURBULENCE.EQ.0)THEN
	VISCL(4)=ZERO;VISCL(3)=ZERO
	LAML(4)=ZERO;LAML(3)=ZERO
	ELSE
	
	
	
!Modified on 19/6/2013
  SELECT CASE(TURBULENCEMODEL)
  
   CASE(1)
	  TURBMV(1)=EDDYFL(2)
	  TURBMV(2)=EDDYFR(2)  
         chi     = abs ( max(TURBMV(1),tolepsma) / max(VISCL(1),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
	VISCL(3) = TURBMV(1)*fv1
	chi     = abs ( max(TURBMV(2),tolepsma) / max(VISCL(2),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
		VISCL(4) = TURBMV(2)*fv1


  CASE(2)

  
  !FOR LEFT CELL
  
 VORTET(1,1:3) = EDDYFL(4:6)
 VORTET(2,1:3) = EDDYFL(7:9) 
 VORTET(3,1:3) = EDDYFL(10:12)

  ux = Vortet(1,1);uy = Vortet(1,2);uz = Vortet(1,3)
  vx = Vortet(2,1);vy = Vortet(2,2);vz = Vortet(2,3)
  wx = Vortet(3,1);wy = Vortet(3,2);wz = Vortet(3,3)

  DO IHGT=1,3
  DO IHGJ=1,3
  TVORT(IHGT,IHGJ)=VORTET(IHGJ,IHGT)
  END DO
  END DO

sVORT=0.5*(VORTET+TVORT)
SNORM=SQRT(2.0*((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+(SVORT(1,3)*SVORT(1,3))+&
	       (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))+(SVORT(2,3)*SVORT(2,3))+& 
	       (SVORT(3,1)*SVORT(3,1))+(SVORT(3,2)*SVORT(3,2))+(SVORT(3,3)*SVORT(3,3))))
		   
 wally=EDDYFL(1)
 rho_0=LEFTV(1)
 k_0=MAX(tolepsma,EDDYFL(2)/LEFTV(1))
 om_0=max(EDDYFL(3)/LEFTV(1),ufreestream/charlength/10.0)


 dervk_dervom=(EDDYFL(13)*EDDYFL(16))+(EDDYFL(14)*EDDYFL(17))+(EDDYFL(15)*EDDYFL(18))
		   
D_omplus=max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    
 phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*VISCL(1)/(rho_0*wally*wally*om_0))
 phi_1=min(phi_2, 4.0*rho_0*k_0/(sigma_om2*D_omplus*wally*wally)) 

F_1=tanh(phi_1**4)
F_2=tanh(phi_2**2)
RE_T_SST=RHO_0*K_0/(VISCL(1)*OM_0)
beta_i=F_1*beta_i1+(1.0-F_1)*beta_i2
alpha_star0=beta_i/3.0    
alpha_inf=F_1*alpha_inf1+(1.0-F_1)*alpha_inf2
alpha_star=alpha_starinf*(alpha_star0+Re_t_SST/R_k_SST)/(1.0+Re_t_SST/R_k_SST)


VISCL(3)=rho_0*k_0/om_0/max(1.0/alpha_star,SNORM*F_2/(aa_1*om_0))


!Added 20/6/2013
sigma_k_l=sigma_k1/F_1+sigma_k2/F_2
sigma_om_l=sigma_om1/F_1+sigma_om2/F_2



  IF (EDDYFR(1).GT.0.0)THEN
!FOR RIGHT
 VORTET(1,1:3) = EDDYFR(4:6)
 VORTET(2,1:3) = EDDYFR(7:9) 
 VORTET(3,1:3) = EDDYFR(10:12)

  ux = Vortet(1,1);uy = Vortet(1,2);uz = Vortet(1,3)
  vx = Vortet(2,1);vy = Vortet(2,2);vz = Vortet(2,3)
  wx = Vortet(3,1);wy = Vortet(3,2);wz = Vortet(3,3)

  DO IHGT=1,3
  DO IHGJ=1,3
  TVORT(IHGT,IHGJ)=VORTET(IHGJ,IHGT)
  END DO
  END DO

sVORT=0.5*(VORTET+TVORT)
SNORM=SQRT(2.0*((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+(SVORT(1,3)*SVORT(1,3))+&
	       (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))+(SVORT(2,3)*SVORT(2,3))+& 
	       (SVORT(3,1)*SVORT(3,1))+(SVORT(3,2)*SVORT(3,2))+(SVORT(3,3)*SVORT(3,3))))
		   
 wally=EDDYFR(1)
 rho_0=RIGHTV(1)
 k_0=EDDYFR(2)/RIGHTV(1)
 om_0=max(EDDYFR(3)/RIGHTV(1),1.0e-6)

 ! EDDYFL(13:15)=ILOCAL_RECON3(K)%GRADS(4,1:3)
!EDDYFL(16:18)=ILOCAL_RECON3(K)%GRADS(5,1:3)

 dervk_dervom=(EDDYFR(13)*EDDYFR(16))+(EDDYFR(14)*EDDYFR(17))+(EDDYFR(15)*EDDYFR(18))
		   
D_omplus=max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    !I need derivative of k
 phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*VISCL(2)/(rho_0*wally*wally*om_0))
 phi_1=min(phi_2, 4.0*rho_0*k_0/(sigma_om2*D_omplus*wally*wally)) 

F_1=tanh(phi_1**4)
F_2=tanh(phi_2**2)
RE_T_SST=RHO_0*K_0/(VISCL(2)*OM_0)
beta_i=F_1*beta_i1+(1.0-F_1)*beta_i2
alpha_star0=beta_i/3.0    
alpha_inf=F_1*alpha_inf1+(1.0-F_1)*alpha_inf2
alpha_star=alpha_starinf*(alpha_star0+Re_t_SST/R_k_SST)/(1.0+Re_t_SST/R_k_SST)


VISCL(4)=rho_0*k_0/om_0/max(1.0/alpha_star,SNORM*F_2/(aa_1*om_0))

!Added 20/6/2013
sigma_k_r=sigma_k1/F_1+sigma_k2/F_2
sigma_om_r=sigma_om1/F_1+sigma_om2/F_2


ELSE
VISCL(4)=-VISCL(3)
SIGMA_K_R=SIGMA_K_L
SIGMA_OM_R=SIGMA_OM_L



END IF


END SELECT
		  
		  
    
		  Viscl(3) = MIN(10000000*visc,VISCL(3))  
		  Viscl(4) = MIN(10000000*visc,VISCL(4))		  
  

	 LAML(3)=( VISCL(3)*R_GAS*GAMMA/(PRTU*(GAMMA-1)) ) + ( VISCL(1)*R_GAS*GAMMA/(PRANDTL*(GAMMA-1)) )
	 LAML(4)=( VISCL(4)*R_GAS*GAMMA/(PRTU*(GAMMA-1)) ) + ( VISCL(2)*R_GAS*GAMMA/(PRANDTL*(GAMMA-1)) )
	 VISCL(3)=MAX(0.0D0,VISCL(3))
	 VISCL(4)=MAX(0.0D0,VISCL(4))
	 
	 IF ((TURBMV(1).LT.ZERO).OR.(TURBMV(2).LT.ZERO))THEN
	 VISCL(3)=0.0D0
	 VISCL(4)=0.0D0
	 END IF
	 
	 ETVM(1) = ( 0.5*(VISCL(1)+VISCL(2)) ) +  ( 0.5*(VISCL(3)+VISCL(4)) )




	  
!Added on 20/6/2013---------------------------------------------------------------
!After limiting these variables, we compute the diffusion for the turbulent variables
if (TURBULENCEMODEL .eq. 2) then
!--------EDDYFL/R(19)=GAMMA_k_L/R
!--------EDDYFL/R(20)=GAMMA_om_L/R  

EDDYFL(19)=VISCL(1)+VISCL(3)/sigma_k_l
EDDYFR(19)=VISCL(2)+VISCL(4)/sigma_k_r

EDDYFL(20)=VISCL(1)+VISCL(3)/sigma_om_l
EDDYFR(20)=VISCL(2)+VISCL(4)/sigma_om_r
end if
END IF




  END SUBROUTINE EDDYVISCO






SUBROUTINE EDDYVISCO2D(N,VISCL,LAML,TURBMV,ETVM,EDDYFL,EDDYFR,LEFTV,RIGHTV)
!> @brief
!> This subroutine computes the tubulent eddy viscosity for turbulence models
	IMPLICIT NONE
    REAL,DIMENSION(1:2),INTENT(INOUT)::TURBMV
    real,dimension(1:nof_Variables),INTENT(IN)::leftv,rightv
    REAL,DIMENSION(1),INTENT(INOUT)::ETVM
    REAL,DIMENSION(1:4),INTENT(INOUT)::VISCL,LAML
    REAL,DIMENSION(1:20),INTENT(INOUT)::EDDYFL,EDDYFR
	INTEGER,INTENT(IN)::N
	REAL::ML,SMMM,MTT,chi,tolepsma,chipow3,fv1
	INTEGER:: IHGT, IHGJ
	REAL,DIMENSION(2,2)::VORTET,TVORT,SVORT
	REAL:: ux,uy,uz,vx,vy,vz,wx,wy,wz
	REAL:: wally, D_omplus, phi_2, phi_1, F_1, F_2, k_0, om_0,alpha_inf, alpha_star, Mu_turb
	REAL:: sigma_k_L, sigma_om_L, sigma_k_r, sigma_om_r
	REAL::SNORM,dervk_dervom,Re_t_SST,RHO_0,beta_i




	tolepsma = TOLSMALL
	
	IF (TURBULENCE.EQ.0)THEN
	VISCL(4)=ZERO;VISCL(3)=ZERO
	LAML(4)=ZERO;LAML(3)=ZERO
	ELSE
	
	
	
!Modified on 19/6/2013
  SELECT CASE(TURBULENCEMODEL)
  
   CASE(1)
            if (ispal.eq.2)then
	  TURBMV(1)=EDDYFL(2)
	  TURBMV(2)=EDDYFR(2)  
	  
	  chi     = abs ((TURBMV(1)) / (VISCL(1)))
         !chi     = abs ( max(TURBMV(1),tolepsma) / max(VISCL(1),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
	VISCL(3) = TURBMV(1)*fv1
	 chi     = abs ((TURBMV(2)) / (VISCL(2)))
! 	chi     = abs ( max(TURBMV(2),tolepsma) / max(VISCL(2),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
		VISCL(4) = TURBMV(2)*fv1
        else
        TURBMV(1)=EDDYFL(2)
	  TURBMV(2)=EDDYFR(2)  
	  
	 
         chi     = abs ( max(TURBMV(1),tolepsma) / max(VISCL(1),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
	VISCL(3) = TURBMV(1)*fv1
! 	 chi     = abs ((TURBMV(2)) / (VISCL(2)))
 	chi     = abs ( max(TURBMV(2),tolepsma) / max(VISCL(2),tolepsma))
         chipow3 = chi * chi * chi
         fv1 = chipow3 / (chipow3 + (cv1*cv1*cv1))
		VISCL(4) = TURBMV(2)*fv1
        
        
        end if

  CASE(2)

  
  !FOR LEFT CELL
  
 VORTET(1,1:2) = EDDYFL(4:5)
 VORTET(2,1:2) = EDDYFL(6:7) 
 

  ux = Vortet(1,1);uy = Vortet(1,2)
  vx = Vortet(2,1);vy = Vortet(2,2)
 

  DO IHGT=1,2
  DO IHGJ=1,2
  TVORT(IHGT,IHGJ)=VORTET(IHGJ,IHGT)
  END DO
  END DO

sVORT=0.5*(VORTET+TVORT)
SNORM=SQRT(2.0*((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+&
	       (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))))
		   
 wally=EDDYFL(1)
 rho_0=LEFTV(1)
 k_0=MAX(tolepsma,EDDYFL(2)/LEFTV(1))
 om_0=max(EDDYFL(3)/LEFTV(1),ufreestream/charlength/10.0)


 dervk_dervom=(EDDYFL(8)*EDDYFL(10))+(EDDYFL(9)*EDDYFL(11))
		   
D_omplus=max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    
 phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*VISCL(1)/(rho_0*wally*wally*om_0))
 phi_1=min(phi_2, 4.0*rho_0*k_0/(sigma_om2*D_omplus*wally*wally)) 

F_1=tanh(phi_1**4)
F_2=tanh(phi_2**2)
RE_T_SST=RHO_0*K_0/(VISCL(1)*OM_0)
beta_i=F_1*beta_i1+(1.0-F_1)*beta_i2
alpha_star0=beta_i/3.0    
alpha_inf=F_1*alpha_inf1+(1.0-F_1)*alpha_inf2
alpha_star=alpha_starinf*(alpha_star0+Re_t_SST/R_k_SST)/(1.0+Re_t_SST/R_k_SST)


VISCL(3)=rho_0*k_0/om_0/max(1.0/alpha_star,SNORM*F_2/(aa_1*om_0))


!Added 20/6/2013
sigma_k_l=sigma_k1/F_1+sigma_k2/F_2
sigma_om_l=sigma_om1/F_1+sigma_om2/F_2



  IF (EDDYFR(1).GT.0.0)THEN
VORTET(1,1:2) = EDDYFL(4:5)
 VORTET(2,1:2) = EDDYFL(6:7) 
 

  ux = Vortet(1,1);uy = Vortet(1,2)
  vx = Vortet(2,1);vy = Vortet(2,2)
 

  DO IHGT=1,2
  DO IHGJ=1,2
  TVORT(IHGT,IHGJ)=VORTET(IHGJ,IHGT)
  END DO
  END DO

sVORT=0.5*(VORTET+TVORT)
SNORM=SQRT(2.0*((SVORT(1,1)*SVORT(1,1))+(SVORT(1,2)*SVORT(1,2))+&
	       (SVORT(2,1)*SVORT(2,1))+(SVORT(2,2)*SVORT(2,2))))
		   
 wally=EDDYFR(1)
 rho_0=RIGHTV(1)
 k_0=EDDYFR(2)/RIGHTV(1)
 om_0=max(EDDYFR(3)/RIGHTV(1),1.0e-6)

 ! EDDYFL(13:15)=ILOCAL_RECON3(K)%GRADS(4,1:3)
!EDDYFL(16:18)=ILOCAL_RECON3(K)%GRADS(5,1:3)

 dervk_dervom=(EDDYFL(8)*EDDYFL(10))+(EDDYFL(9)*EDDYFL(11))
		   
D_omplus=max(2*rho_0/sigma_om2/om_0*dervk_dervom, 1.0e-10)    !I need derivative of k
 phi_2=max(sqrt(k_0)/(0.09*om_0*wally),500.0*VISCL(2)/(rho_0*wally*wally*om_0))
 phi_1=min(phi_2, 4.0*rho_0*k_0/(sigma_om2*D_omplus*wally*wally)) 

F_1=tanh(phi_1**4)
F_2=tanh(phi_2**2)
RE_T_SST=RHO_0*K_0/(VISCL(2)*OM_0)
beta_i=F_1*beta_i1+(1.0-F_1)*beta_i2
alpha_star0=beta_i/3.0    
alpha_inf=F_1*alpha_inf1+(1.0-F_1)*alpha_inf2
alpha_star=alpha_starinf*(alpha_star0+Re_t_SST/R_k_SST)/(1.0+Re_t_SST/R_k_SST)


VISCL(4)=rho_0*k_0/om_0/max(1.0/alpha_star,SNORM*F_2/(aa_1*om_0))

!Added 20/6/2013
sigma_k_r=sigma_k1/F_1+sigma_k2/F_2
sigma_om_r=sigma_om1/F_1+sigma_om2/F_2


ELSE
VISCL(4)=-VISCL(3)
SIGMA_K_R=SIGMA_K_L
SIGMA_OM_R=SIGMA_OM_L



END IF


END SELECT
		  
		  
    
		  Viscl(3) = MIN(10000000*visc,VISCL(3))  
		  Viscl(4) = MIN(10000000*visc,VISCL(4))		  
  

	 LAML(3)=( VISCL(3)*R_GAS*GAMMA/(PRTU*(GAMMA-1)) ) + ( VISCL(1)*R_GAS*GAMMA/(PRANDTL*(GAMMA-1)) )
	 LAML(4)=( VISCL(4)*R_GAS*GAMMA/(PRTU*(GAMMA-1)) ) + ( VISCL(2)*R_GAS*GAMMA/(PRANDTL*(GAMMA-1)) )
	 VISCL(3)=MAX(0.0D0,VISCL(3))
	 VISCL(4)=MAX(0.0D0,VISCL(4))
	 
	 IF ((TURBMV(1).LT.ZERO).OR.(TURBMV(2).LT.ZERO))THEN
	 VISCL(3)=0.0D0
	 VISCL(4)=0.0D0
	 END IF
	 
	 ETVM(1) = ( 0.5*(VISCL(1)+VISCL(2)) ) +  ( 0.5*(VISCL(3)+VISCL(4)) )




	  
!Added on 20/6/2013---------------------------------------------------------------
!After limiting these variables, we compute the diffusion for the turbulent variables
if (TURBULENCEMODEL .eq. 2) then
!--------EDDYFL/R(19)=GAMMA_k_L/R
!--------EDDYFL/R(20)=GAMMA_om_L/R  

EDDYFL(12)=VISCL(1)+VISCL(3)/sigma_k_l
EDDYFR(13)=VISCL(2)+VISCL(4)/sigma_k_r

EDDYFL(12)=VISCL(1)+VISCL(3)/sigma_om_l
EDDYFR(13)=VISCL(2)+VISCL(4)/sigma_om_r
end if
END IF




  END SUBROUTINE EDDYVISCO2d



  
  
  

SUBROUTINE TRAJECTORIES
IMPLICIT NONE
INTEGER::I,J,K,TRAJ1,TRAJ2,TRAJ3,TRAJ4,kmaxe,writeid,writeconf,num_Vg
REAL::WIN1,WIN2,WIN3,WIN4,POST,POST1,POST2,POST3,POST4
real,dimension(1:4)::pos_l,pos_g
REAL,DIMENSION(1:NOF_VARIABLES)::LEFTV
REAL::MP_PINFl,GAMMAL
KMAXE=XMPIELRANK(N)
POST1=TOLBIG
POST2=TOLBIG
POST3=-TOLBIG
traj1=0
TRAJ2=0
TRAJ3=0

if (dimensiona.eq.3)then
num_Vg=5
else
num_Vg=4

end if


IF (INITCOND.EQ.157)THEN
pos_l(1:2)=zero
pos_g(1:2)=zero
DO I=1,KMAXE
      LEFTV(1:NOF_VARIABLES)=U_C(I)%VAL(1,1:NOF_VARIABLES)
      if (dimensiona.eq.3)then
      CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)
      else
      CALL cons2prim(N,leftv,MP_PINFl,gammal)
      end if

     IF ((LEFTV(NOF_VARIABLES)).GT.0.0D0)THEN  !total volume of gas evolution
            pos_l(2)=pos_l(2)+(LEFTV(NOF_VARIABLES))*ielem(n,i)%totvolume
     end if
     IF (LEFTV(num_Vg).GT.pos_l(1))THEN  !total volume of gas evolution
            pos_l(1)=MAX(pos_l(1),LEFTV(num_Vg))
     end if



END Do
     !find position globally
CALL MPI_ALLREDUCE(pos_l(1),pos_g(1),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)
CALL MPI_ALLREDUCE(pos_l(2),pos_g(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)





IF (n.eq.0)THEN

OPEN(70,FILE='Volumex.DAT',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7)')T,POS_G(1),POS_G(2)
close(70)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)






END IF




IF (INITCOND.EQ.430)THEN
pos_l(1:2)=zero
pos_g(1:2)=zero
DO I=1,KMAXE
    IF (U_C(I)%VAL(1,7).GT.0.1D0)THEN   !VOLUME FRACTION OF GAS to be used for lowest location tracking
        IF (IELEM(N,I)%YYC.Le.POST3)THEN
            POST3=IELEM(N,I)%YYC
            TRAJ1=I
        END IF
     END IF   
     IF (U_C(I)%VAL(1,7).GT.0.0D0)THEN  !total volume of gas evolution
            pos_l(2)=pos_l(2)+U_C(I)%VAL(1,7)*ielem(n,i)%totvolume
     end if
END Do

pos_l(1)=post3  !lowest position in my local  cpu


!find position globally
CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
CALL MPI_ALLREDUCE(pos_l(1),pos_g(1),1,MPI_DOUBLE_PRECISION,MPI_min,MPI_COMM_WORLD,IERROR)

!NOW SELECT WHICH CPU HAS THE SMALLEST
writeid=100000000
IF (ABS(POS_G(1)-POS_L(1)).LE.TOLSMALL/1000)THEN
IF (TRAJ1.GT.0)THEN
writeid=n
END IF
end if
!AND IF MORE THAN ONE, SELECT THE SMALLEST ID ONE
call mpi_barrier(mpi_comm_world,ierror)
CALL MPI_ALLREDUCE(writeid,writeconf,1,MPI_INTEGER,MPI_min,MPI_COMM_WORLD,IERROR)



if (writeid.EQ.writeconf)then

leftv(1:nof_Variables)=U_C(TRAJ1)%VAL(1,1:nof_Variables)

CALL cons2prim(N,leftv,MP_PINFl,gammal)
OPEN(70,FILE='position.DAT',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POS_G(1),LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5),LEFTV(6),LEFTV(7)
close(70)
!the CPU that holds this cell will write its position
end if



CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  
CALL MPI_ALLREDUCE(pos_l(2:2),pos_g(2:2),1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

IF (n.eq.0)THEN

OPEN(70,FILE='volume.DAT',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7)')T,POS_G(2)
close(70)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)




END IF


if (initcond.eq.405)then


DO I=1,KMAXE
    IF (DIMENSIONA.EQ.3)THEN
    IF (U_C(I)%VAL(1,8).GT.0.1D0)THEN
        IF (IELEM(N,I)%XXC.le.POST1)THEN
            POST1=IELEM(N,I)%XXC
            TRAJ1=I
        END IF
        IF (((IELEM(N,I)%YYC.LE.0.0515).AND.(IELEM(N,I)%YYC.GE.0.0485)).and.((IELEM(N,I)%zzC.LE.0.0515).AND.(IELEM(N,I)%zzC.GE.0.0485)))THEN
            IF(IELEM(N,I)%XXC.lE.POST2)THEN
                POST2=IELEM(N,I)%XXC
                TRAJ2=i
            end if
            if (ielem(n,i)%xxc.ge.post3)then
                post3=ielem(n,i)%xxc
                traj3=i
            end if
        END IF
    END IF
    ELSE
    IF (U_C(I)%VAL(1,7).Ge.0.1D0)THEN
        IF (IELEM(N,I)%XXC.lE.POST1)THEN
            POST1=IELEM(N,I)%XXC
            TRAJ1=I
        END IF
        IF (((IELEM(N,I)%YYC.LE.0.0515).AND.(IELEM(N,I)%YYC.GE.0.0485)))THEN
            IF(IELEM(N,I)%XXC.lE.POST2)THEN
                POST2=IELEM(N,I)%XXC
                TRAJ2=i
            end if
            if (ielem(n,i)%xxc.GE.post3)then
                post3=ielem(n,i)%xxc
                traj3=i
            end if
        END IF
    END IF


    END IF
END DO

pos_l(1)=post1
pos_l(2)=post2
pos_l(3)=post3

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  
CALL MPI_ALLREDUCE(pos_l(1:2),pos_g(1:2),2,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,IERROR)

!NOW SELECT WHICH CPU HAS THE SMALLEST
writeid=100000000
writeconf=-1
IF (ABS(POS_G(1)-POS_L(1)).LE.TOLSMALL/1000)THEN
IF (TRAJ1.GT.0)THEN
writeid=n
END IF
end if
!AND IF MORE THAN ONE, SELECT THE SMALLEST ID ONE

!IF (TRAJ1.GT.0)THEN
CALL MPI_ALLREDUCE(writeid,writeconf,1,MPI_INTEGER,MPI_min,MPI_COMM_WORLD,IERROR)
!END IF



if (writeid.EQ.writeconf)then

        IF (traj1.gt.0)then
leftv(1:nof_Variables)=U_C(TRAJ1)%VAL(1,1:nof_Variables)


IF (DIMENSIONA.EQ.3)THEN
CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)

OPEN(70,FILE='pos1.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POS_G(1),LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5),LEFTV(6),LEFTV(7),LEFTV(8)
close(70)
ELSE
CALL cons2prim(N,leftv,MP_PINFl,gammal)

OPEN(70,FILE='pos1.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POS_G(1),LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5),LEFTV(6),LEFTV(7)
close(70)

END IF




end if




END IF



writeid=100000000
writeconf=-1
IF (ABS(POS_G(2)-POS_L(2)).LE.TOLSMALL/1000)THEN
IF (TRAJ2.GT.0)THEN
writeid=n
END IF
end if

!IF (TRAJ2.GT.0)THEN
CALL MPI_ALLREDUCE(writeid,writeconf,1,MPI_INTEGER,MPI_min,MPI_COMM_WORLD,IERROR)
!END IF


if (writeid.EQ.writeconf)then

        if (traj2.gt.0)then
leftv(1:nof_Variables)=U_C(TRAJ2)%VAL(1,1:nof_Variables)


IF (DIMENSIONA.EQ.3)THEN

CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)

OPEN(71,FILE='pos2.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(71,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POS_G(2),LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5),LEFTV(6),LEFTV(7),LEFTV(8)
close(71)

Else
CALL cons2prim(N,leftv,MP_PINFl,gammal)

OPEN(70,FILE='pos2.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POS_G(2),LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5),LEFTV(6),LEFTV(7)
close(70)


END IF



end if
END IF


CALL MPI_ALLREDUCE(pos_l(3),pos_g(3),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)



writeid=100000000
writeconf=-1
IF (ABS(POS_G(3)-POS_L(3)).LE.TOLSMALL/1000)THEN
IF (TRAJ3.GT.0)THEN
writeid=n
END IF
end if


!IF (TRAJ3.GT.0)THEN
CALL MPI_ALLREDUCE(writeid,writeconf,1,MPI_INTEGER,MPI_min,MPI_COMM_WORLD,IERROR)
!END IF



if (writeid.EQ.writeconf)then
        if (traj3.gt.0)then
leftv(1:nof_Variables)=U_C(TRAJ3)%VAL(1,1:nof_Variables)


IF (DIMENSIONA.EQ.3)THEN



CALL CONS2PRIM(N,leftv,MP_PINFl,gammal)

OPEN(71,FILE='pos3.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(71,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POS_G(3),LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5),LEFTV(6),LEFTV(7),LEFTV(8)
close(71)


Else
CALL cons2prim(N,leftv,MP_PINFl,gammal)

OPEN(70,FILE='pos3.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POS_G(3),LEFTV(1),LEFTV(2),LEFTV(3),LEFTV(4),LEFTV(5),LEFTV(6),LEFTV(7)
close(70)

END IF




end if






END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


end if

if ((initcond.eq.408).or.(initcond.eq.422))then

pos_l(1:2)=zero
pos_g(1:2)=zero

DO I=1,KMAXE
    IF (U_C(I)%VAL(1,8).GT.0.0D0)THEN
            pos_l(1)=pos_l(1)+ielem(n,i)%totvolume
            pos_l(2)=pos_l(2)+U_C(I)%VAL(1,8)*ielem(n,i)%totvolume
    end if
END DO


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  
CALL MPI_ALLREDUCE(pos_l(1:2),pos_g(1:2),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

IF (n.eq.0)THEN

OPEN(70,FILE='volume.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7)')T,POS_G(1),pos_g(2)
close(70)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


end if


if ((initcond.eq.411).or.(initcond.eq.444))then

pos_l(1:3)=zero
pos_g(1:3)=zero

DO I=1,KMAXE
    leftv(1:nof_variables)=U_C(I)%VAL(1,1:nof_Variables)
    IF (DIMENSIONA.EQ.2)THEN
    call cons2prim(N,leftv,MP_PINFl,gammal)
    pos_l(3)=max(abs(leftv(4)),pos_l(3))

    IF (U_C(I)%VAL(1,7).GT.0.0D0)THEN
            pos_l(1)=pos_l(1)+ielem(n,i)%totvolume
            pos_l(2)=pos_l(2)+U_C(I)%VAL(1,7)*ielem(n,i)%totvolume
    end if
    ELSE
    call CONS2PRIM(N,leftv,MP_PINFl,gammal)
    pos_l(3)=max(abs(leftv(5)),pos_l(3))

    IF (U_C(I)%VAL(1,8).GT.0.0D0)THEN
            pos_l(1)=pos_l(1)+ielem(n,i)%totvolume
            pos_l(2)=pos_l(2)+U_C(I)%VAL(1,7)*ielem(n,i)%totvolume
    end if


    END IF
    
END DO


CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)
  
CALL MPI_ALLREDUCE(pos_l(1:2),pos_g(1:2),2,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,IERROR)

CALL MPI_ALLREDUCE(pos_l(3),pos_g(3),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,IERROR)

IF (n.eq.0)THEN

OPEN(70,FILE='volume.dat',FORM='FORMATTED',ACTION='WRITE',POSITION='APPEND')
WRITE(70,'(E14.7,1X,E14.7,1X,E14.7,1X,E14.7)')T,POS_G(1),pos_g(2),POS_G(3)
close(70)

END IF

CALL MPI_BARRIER(MPI_COMM_WORLD,IERROR)


end if


END SUBROUTINE

SUBROUTINE FLUX2DX(FLUX_TERM_X,LEFTV)
IMPLICIT NONE
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI, IE1, MP_STIFF, MP_DENSITY,GAMMAL,GAMMAR,sum1,sum2,sum3
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::LEFTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::FLUX_TERM_X
REAL,DIMENSION(NOF_SPECIES)::mp_vft

IF(MULTISPECIES.EQ.1)THEN


            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+LEFTV(nof_Variables+rg_i)
            end do
            MP_DENSITY=sum1
            do rg_i=1,nof_SPECIES-1
            sum2=sum2+LEFTV(DIMENSIONA+2+nof_species+rg_i)
            MP_AR(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
            mp_vft(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)
            end do
            MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
            mp_vft(nof_SPECIES)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+MP_AR(rg_i)
            end do
            GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
             sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
            end do

            MP_STIFF=sum3*(GAMMAL-1.0D0)
             sum2=zero
            do rg_i=1,nof_SPECIES
            sum2=sum2+(mp_vft(rg_i)*MP_PINF(rg_i))
            end do

 
R=MP_DENSITY
U=LEFTV(2)
V=LEFTV(3)
P=LEFTV(4)

!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2))
!INTERNAL ENERGY 
IE1=((P+MP_STIFF)/((GAMMAL-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IE1)
FLUX_TERM_X(1)=R*U
FLUX_TERM_X(2)=(R*(U**2))+P
FLUX_TERM_X(3)=R*U*V
FLUX_TERM_X(4)=U*(E+P)
FLUX_TERM_X(5:NOF_VARIABLES)=LEFTV(5:NOF_VARIABLES)*U
ELSE

R=LEFTV(1)
U=LEFTV(2)
V=LEFTV(3)
P=LEFTV(4)
GM=GAMMA
!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)
FLUX_TERM_X(1)=R*U
FLUX_TERM_X(2)=(R*(U**2))+P
FLUX_TERM_X(3)=R*U*V
FLUX_TERM_X(4)=U*(E+P)

END IF

END SUBROUTINE


SUBROUTINE FLUX2DY(FLUX_TERM_Y,LEFTV)
IMPLICIT NONE
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI, IE1, MP_STIFF, MP_DENSITY,GAMMAL,GAMMAR,sum1,sum2,sum3
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE,mp_vft
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::LEFTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::FLUX_TERM_Y
IF(MULTISPECIES.EQ.1)THEN


            sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+LEFTV(nof_Variables+rg_i)
            end do
            MP_DENSITY=sum1
            do rg_i=1,nof_SPECIES-1
            sum2=sum2+LEFTV(DIMENSIONA+2+nof_species+rg_i)
            MP_AR(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
            mp_vft(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)
            end do
            MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
            mp_vft(nof_SPECIES)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+MP_AR(rg_i)
            end do
            GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
             sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
            end do

            MP_STIFF=sum3*(GAMMAL-1.0D0)
             sum2=zero
            do rg_i=1,nof_SPECIES
            sum2=sum2+(mp_vft(rg_i)*MP_PINF(rg_i))
            end do










! MP_AR(1)=LEFTV(7)/(GAMMA_IN(1)-1.0D0)
! MP_AR(2)=(1.0D0-LEFTV(7))/(GAMMA_IN(2)-1.0D0)
! GAMMAL=(1.0D0/(MP_AR(1)+MP_AR(2)))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTIO
! MP_STIFF=((LEFTV(7)*(GAMMA_IN(1)/(GAMMA_IN(1)-1.0D0))*MP_PINF(1))+((1.0D0-LEFTV(7))*(GAMMA_IN(2)/(GAMMA_IN(2)-1.0D0))*MP_PINF(2)))*(GAMMAL-1.0D0)
! MP_DENSITY = LEFTV(5)+LEFTV(6)
 
R=MP_DENSITY
U=LEFTV(2)
V=LEFTV(3)
P=LEFTV(4)

!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2))
!INTERNAL ENERGY 
IE1=((P+MP_STIFF)/((GAMMAL-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IE1)
FLUX_TERM_Y(1)=R*V
FLUX_TERM_Y(2)=R*U*V
FLUX_TERM_Y(3)=(R*(V**2))+P
FLUX_TERM_Y(4)=V*(E+P)
FLUX_TERM_Y(5:7)=LEFTV(5:7)*V
ELSE


R=LEFTV(1)
U=LEFTV(2)
V=LEFTV(3)
P=LEFTV(4)
GM=GAMMA
!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2))
!INTERNAL ENERGY
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)

FLUX_TERM_Y(1)=R*v
FLUX_TERM_Y(2)=R*u*v
FLUX_TERM_Y(3)=(R*v*v)+P
FLUX_TERM_Y(4)=v*(E+P)

ENDIF

END subroutine







SUBROUTINE FLUX3Dx(FLUX_TERM_X,LEFTV)
IMPLICIT NONE
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI,IE1,MP_STIFF,MP_DENSITY,GAMMAL,GAMMAR,sum1,sum2,sum3
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE,mp_vft
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::LEFTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::FLUX_TERM_X

IF(MULTISPECIES.EQ.1)THEN

 sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+LEFTV(nof_Variables+rg_i)
            end do
            MP_DENSITY=sum1
            do rg_i=1,nof_SPECIES-1
            sum2=sum2+LEFTV(DIMENSIONA+2+nof_species+rg_i)
            MP_AR(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
            mp_vft(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)
            end do
            MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
            mp_vft(nof_SPECIES)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+MP_AR(rg_i)
            end do
            GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
             sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
            end do

            MP_STIFF=sum3*(GAMMAL-1.0D0)
             sum2=zero
            do rg_i=1,nof_SPECIES
            sum2=sum2+(mp_vft(rg_i)*MP_PINF(rg_i))
            end do



! MP_AR(1)=LEFTV(8)/(GAMMA_IN(1)-1.0D0)
! MP_AR(2)=(1.0D0-LEFTV(8))/(GAMMA_IN(2)-1.0D0)
! GAMMAL=(1.0D0/(MP_AR(1)+MP_AR(2)))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTIO
! MP_STIFF=((LEFTV(8)*(GAMMA_IN(1)/(GAMMA_IN(1)-1.0D0))*MP_PINF(1))+((1.0D0-LEFTV(8))*(GAMMA_IN(2)/(GAMMA_IN(2)-1.0D0))*MP_PINF(2)))*(GAMMAL-1.0D0)


! MP_DENSITY = LEFTV(6)+LEFTV(7)

R=MP_DENSITY

U=LEFTV(2)
V=LEFTV(3)
W=LEFTV(4)
P=LEFTV(5)

!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY
IE1=((P+MP_stiff)/((GAMMAL-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IE1)
FLUX_TERM_X(1)=R*U
FLUX_TERM_X(2)=(R*(U**2))+P
FLUX_TERM_X(3)=R*U*V
FLUX_TERM_X(4)=R*U*W
FLUX_TERM_X(5)=U*(E+P)
FLUX_TERM_X(6:8)=LEFTV(6:8)*U
ELSE
R=LEFTV(1)
U=LEFTV(2)
V=LEFTV(3)
w=LEFTV(4)
P=LEFTV(5)
GM=GAMMA
!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)

FLUX_TERM_X(1)=R*U
FLUX_TERM_X(2)=(R*(U**2))+P
FLUX_TERM_X(3)=R*U*V
FLUX_TERM_X(4)=R*U*W
FLUX_TERM_X(5)=U*(E+P)

end if

END SUBROUTINE



SUBROUTINE FLUX3DY(FLUX_TERM_Y,LEFTV)
IMPLICIT NONE
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI,IE1, MP_STIFF,MP_DENSITY,GAMMAL,GAMMAR,sum1,sum2,sum3
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE,mp_vft
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::LEFTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::FLUX_TERM_Y

IF(MULTISPECIES.EQ.1)THEN



 sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+LEFTV(nof_Variables+rg_i)
            end do
            MP_DENSITY=sum1
            do rg_i=1,nof_SPECIES-1
            sum2=sum2+LEFTV(DIMENSIONA+2+nof_species+rg_i)
            MP_AR(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
            mp_vft(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)
            end do
            MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
            mp_vft(nof_SPECIES)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+MP_AR(rg_i)
            end do
            GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
             sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
            end do

            MP_STIFF=sum3*(GAMMAL-1.0D0)
             sum2=zero
            do rg_i=1,nof_SPECIES
            sum2=sum2+(mp_vft(rg_i)*MP_PINF(rg_i))
            end do







! MP_AR(1)=LEFTV(8)/(GAMMA_IN(1)-1.0D0)
! MP_AR(2)=(1.0D0-LEFTV(8))/(GAMMA_IN(2)-1.0D0)
! GAMMAL=(1.0D0/(MP_AR(1)+MP_AR(2)))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTIO
! MP_STIFF=((LEFTV(8)*(GAMMA_IN(1)/(GAMMA_IN(1)-1.0D0))*MP_PINF(1))+((1.0D0-LEFTV(8))*(GAMMA_IN(2)/(GAMMA_IN(2)-1.0D0))*MP_PINF(2)))*(GAMMAL-1.0D0)
!
!
! MP_DENSITY = LEFTV(6)+LEFTV(7)

R=MP_DENSITY


U=LEFTV(2)
V=LEFTV(3)
W=LEFTV(4)
P=LEFTV(5)

!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY
IE1=((P+MP_stiff)/((GAMMAL-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IE1)
FLUX_TERM_Y(1)=R*V
FLUX_TERM_Y(2)=R*U*V
FLUX_TERM_Y(3)=(R*(V**2))+P
FLUX_TERM_Y(4)=R*v*W
FLUX_TERM_Y(5)=V*(E+P)
FLUX_TERM_Y(6:8)=LEFTV(6:8)*V
ELSE
R=LEFTV(1)
U=LEFTV(2)
V=LEFTV(3)
w=LEFTV(4)
P=LEFTV(5)
GM=GAMMA
!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)

FLUX_TERM_Y(1)=R*v
FLUX_TERM_Y(2)=R*u*v
FLUX_TERM_Y(3)=(R*(v*v))+p
FLUX_TERM_Y(4)=R*v*W
FLUX_TERM_Y(5)=v*(E+P)

end if
END SUBROUTINE

SUBROUTINE FLUX3DZ(FLUX_TERM_Z,LEFTV)
IMPLICIT NONE
REAL::P,U,V,W,E,R,S,GM,SKIN,IEN,PI,IE1, MP_STIFF,MP_DENSITY,GAMMAL,GAMMAR,sum1,sum2,sum3
REAL,DIMENSION(NOF_SPECIES)::MP_AR,MP_IE,mp_vft
INTEGER::rg_i,rg_j
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::LEFTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::FLUX_TERM_Z

IF(MULTISPECIES.EQ.1)THEN



 sum1=zero;
            sum2=zero;
            sum3=zero
            do rg_i=1,nof_species
                sum1=sum1+LEFTV(nof_Variables+rg_i)
            end do
            MP_DENSITY=sum1
            do rg_i=1,nof_SPECIES-1
            sum2=sum2+LEFTV(DIMENSIONA+2+nof_species+rg_i)
            MP_AR(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)/(GAMMA_IN(rg_i)-1.0D0)
            mp_vft(rg_i)=LEFTV(DIMENSIONA+2+nof_species+rg_i)
            end do
            MP_AR(nof_SPECIES)=(1.0d0-sum2)/(GAMMA_IN(nof_SPECIES)-1.0D0)
            mp_vft(nof_SPECIES)=(1.0d0-sum2)
            sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+MP_AR(rg_i)
            end do
            GAMMAL=(1.0D0/(sum3))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTION
             sum3=zero
            do rg_i=1,nof_SPECIES
            sum3=sum3+(mp_vft(rg_i)*(GAMMA_IN(rg_i)/(GAMMA_IN(rg_i)-1.0D0))*MP_PINF(rg_i))
            end do

            MP_STIFF=sum3*(GAMMAL-1.0D0)
             sum2=zero
            do rg_i=1,nof_SPECIES
            sum2=sum2+(mp_vft(rg_i)*MP_PINF(rg_i))
            end do










! MP_AR(1)=LEFTV(8)/(GAMMA_IN(1)-1.0D0)
! MP_AR(2)=(1.0D0-LEFTV(8))/(GAMMA_IN(2)-1.0D0)
! GAMMAL=(1.0D0/(MP_AR(1)+MP_AR(2)))+1.0D0    !MIXTURE GAMMA ISOBARIC ASSUMPTIO
! MP_STIFF=((LEFTV(8)*(GAMMA_IN(1)/(GAMMA_IN(1)-1.0D0))*MP_PINF(1))+((1.0D0-LEFTV(8))*(GAMMA_IN(2)/(GAMMA_IN(2)-1.0D0))*MP_PINF(2)))*(GAMMAL-1.0D0)
!
! MP_DENSITY = LEFTV(6)+LEFTV(7)
!
 R=MP_DENSITY
U=LEFTV(2)
V=LEFTV(3)
W=LEFTV(4)
P=LEFTV(5)

!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY
IE1=((P+MP_stiff)/((GAMMAL-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IE1)
FLUX_TERM_Z(1)=R*W
FLUX_TERM_Z(2)=R*U*W
FLUX_TERM_Z(3)=R*V*W
FLUX_TERM_Z(4)=(R*(W*W))+p
FLUX_TERM_Z(5)=W*(E+P)
FLUX_TERM_Z(6:8)=LEFTV(6:8)*W
ELSE
R=LEFTV(1)
U=LEFTV(2)
V=LEFTV(3)
w=LEFTV(4)
P=LEFTV(5)
GM=GAMMA
!KINETIC ENERGY FIRST!
SKIN=(OO2)*((U**2)+(V**2)+(W**2))
!INTERNAL ENERGY 
IEN=((P)/((GM-1.0D0)*R))
!TOTAL ENERGY
E=R*(SKIN+IEN)

FLUX_TERM_Z(1)=R*w
FLUX_TERM_Z(2)=r*u*w
FLUX_TERM_Z(3)=R*v*w
FLUX_TERM_Z(4)=(R*(w**2))+P
FLUX_TERM_Z(5)=w*(E+P)
END IF
END SUBROUTINE






SUBROUTINE FLUX_VISC2D(FLUX_TERM_X,FLUX_TERM_Y,LEFTV,LEFTV_DER)
IMPLICIT NONE
REAL:: U, V, UX, UY, VX, VY, TX, TY, TAUXX, TAUXY, TAUYY
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::LEFTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::FLUX_TERM_X,FLUX_TERM_Y
REAL,DIMENSION(1:4)::VISCL,LAML
REAL,DIMENSION(1:NOF_VARIABLES,1:DIMENSIONA),INTENT(IN)::LEFTV_DER

    CALL GET_visc_conduct(N, LEFTV, LEFTV,VISCL,LAML)

       U = LEFTV(2)
       V = LEFTV(3)
      UX = LEFTV_DER(2,1)
      UY = LEFTV_DER(2,2)
      VX = LEFTV_DER(3,1)
      VY = LEFTV_DER(3,2)
      TX = LEFTV_DER(4,1)
      TY = LEFTV_DER(4,2)

    TAUXX = 2.0D0 / 3.0D0 * VISCL(1) * (2 * UX - VY)
    TAUYY = 2.0D0 / 3.0D0 * VISCL(1) * (2 * VY - UX)
    TAUXY = VISCL(1) * (UY + VX)

    FLUX_TERM_X(2) = FLUX_TERM_X(2) - TAUXX
    FLUX_TERM_X(3) = FLUX_TERM_X(3) - TAUXY
    FLUX_TERM_X(4) = FLUX_TERM_X(4) - (U * TAUXX + V * TAUXY + LAML(1) * TX)
    FLUX_TERM_Y(2) = FLUX_TERM_Y(2) - TAUXY
    FLUX_TERM_Y(3) = FLUX_TERM_Y(3) - TAUYY
    FLUX_TERM_Y(4) = FLUX_TERM_Y(4) - (U * TAUXY + V * TAUYY + LAML(1) * TY)


END SUBROUTINE


SUBROUTINE FLUX_VISC3D(FLUX_TERM_X,FLUX_TERM_Y,FLUX_TERM_Z,LEFTV,LEFTV_DER)
IMPLICIT NONE
REAL:: U, V, UX, UY, VX, VY, TX, TY, TAUXX, TAUXY, TAUYY, W, UZ, VZ, WX, WY, WZ, TZ, TAUXZ, TAUYZ, TAUZZ
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::LEFTV
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(INOUT)::FLUX_TERM_X,FLUX_TERM_Y,FLUX_TERM_Z
REAL,DIMENSION(1:NOF_VARIABLES,1:DIMENSIONA),INTENT(IN)::LEFTV_DER
REAL,DIMENSION(1:4)::VISCL,LAML


    CALL GET_visc_conduct(N, LEFTV, LEFTV,VISCL,LAML)

! Variables extrapolated at boundary
       U = LEFTV(2)
       V = LEFTV(3)
       W = LEFTV(4)

      UX = LEFTV_DER(2,1)
      UY = LEFTV_DER(2,2)
      UZ = LEFTV_DER(2,3)

      VX = LEFTV_DER(3,1)
      VY = LEFTV_DER(3,2)
      VZ = LEFTV_DER(3,3)

      WX = LEFTV_DER(4,1)
      WY = LEFTV_DER(4,2)
      WZ = LEFTV_DER(4,3)

      TX = LEFTV_DER(5,1)
      TY = LEFTV_DER(5,2)
      TZ = LEFTV_DER(5,3)

    TAUXX = 2.0D0 / 3.0D0 * VISCL(1) * (2 * UX - VY - WZ)
    TAUYY = 2.0D0 / 3.0D0 * VISCL(1) * (2 * VY - UX - WZ)
    TAUXY = VISCL(1) * (UY + VX)

    TAUXZ = VISCL(1) * (UZ + WX)
    TAUYZ = VISCL(1) * (WY + VZ)
    TAUZZ = 2.0D0 / 3.0D0 * VISCL(1) * (2 * WZ - UX - VY)

    FLUX_TERM_X(2) = FLUX_TERM_X(2) - TAUXX
    FLUX_TERM_X(3) = FLUX_TERM_X(3) - TAUXY
    FLUX_TERM_X(4) = FLUX_TERM_X(4) - TAUXZ
    FLUX_TERM_X(5) = FLUX_TERM_X(5) - (U * TAUXX + V * TAUXY + W * TAUXZ + LAML(1) * TX)

    FLUX_TERM_Y(2) = FLUX_TERM_Y(2) - TAUXY
    FLUX_TERM_Y(3) = FLUX_TERM_Y(3) - TAUYY
    FLUX_TERM_Y(4) = FLUX_TERM_Y(4) - TAUYZ
    FLUX_TERM_Y(5) = FLUX_TERM_Y(5) - (U * TAUXY + V * TAUYY + W * TAUYZ + LAML(1) * TY)

    FLUX_TERM_Z(2) = FLUX_TERM_Z(2) - TAUXZ
    FLUX_TERM_Z(3) = FLUX_TERM_Z(3) - TAUYZ
    FLUX_TERM_Z(4) = FLUX_TERM_Z(4) - TAUZZ
    FLUX_TERM_Z(5) = FLUX_TERM_Z(5) - (U * TAUXZ + V * TAUYZ + W * TAUZZ + LAML(1) * TZ)

END SUBROUTINE








SUBROUTINE DCONS2DPRIM(LEFTV_DER,LEFTV)
IMPLICIT NONE
REAL,DIMENSION(1:NOF_VARIABLES,1:DIMENSIONA),INTENT(INOUT)::LEFTV_DER
REAL,DIMENSION(1:NOF_VARIABLES),INTENT(IN)::LEFTV
INTEGER:: I_DIM, I_VAR


    DO I_DIM = 1, DIMENSIONA
        DO I_VAR = 2, NOF_VARIABLES-1
            LEFTV_DER(I_VAR,I_DIM) = (LEFTV_DER(I_VAR,I_DIM) - LEFTV_DER(1,I_DIM) * LEFTV(I_VAR)) / LEFTV(1) ! UX = (RHOUX - RHOX * U) / RHO
        END DO


        LEFTV_DER(NOF_VARIABLES,I_DIM) = (GAMMA - 1.0D0) * ((LEFTV(1) * LEFTV_DER(NOF_VARIABLES,I_DIM) - LEFTV(NOF_VARIABLES) * LEFTV_DER(1,I_DIM)) / LEFTV(1) ** 2 - DOT_PRODUCT(LEFTV(2:NOF_VARIABLES-1), LEFTV_DER(2:NOF_VARIABLES-1,I_DIM)))
        ! TX = (GAMMA-1) * ((RHO * EX - E * RHOX) / RHO ** 2 - (U * UX + V * VX + W * WX))
    END DO

END SUBROUTINE DCONS2DPRIM















END MODULE FLOW_OPERATIONS
