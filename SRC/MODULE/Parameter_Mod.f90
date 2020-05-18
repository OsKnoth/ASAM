MODULE Parameter_Mod

  USE Kind_Mod

  IMPLICIT NONE

! ---Numerical Parameter-----------------------------------------------------------------------------------------
!  REAL(RealKind), PARAMETER :: Eps        = 1.d-40
!  REAL(RealKind), PARAMETER :: Zero       = 0.0d0
!  REAL(RealKind), PARAMETER :: Half       = 0.5d0
!  REAL(RealKind), PARAMETER :: One        = 1.0d0
!  REAL(RealKind), PARAMETER :: Two        = 2.0d0
!  REAL(RealKind), PARAMETER :: Three      = 3.0d0
!  REAL(RealKind), PARAMETER :: Four       = 4.0d0
  REAL(RealKind) :: Pi        ! 4.0d0*ATAN(1.0d0)    
  REAL(RealKind) :: PiHalf    ! 0.5d0*Pi

! ---Physical Parameter-------------------------------------------------------------------------------------------
  
! ---Thermodynamic parameters---
  REAL(RealKind), PARAMETER :: Cpd        = 1004.0d0   ! specific heat of dry air at constant pressure
  REAL(RealKind), PARAMETER :: Cvd        = 717.0d0    ! specific heat of dry air at constant volume (=Cpd-Rd)
  REAL(RealKind), PARAMETER :: Rd=Cpd-Cvd              ! gas constant for dry air 

  REAL(RealKind), PARAMETER :: Cpv        = 1885.0d0   ! specific heat of vapor at constant pressure
  REAL(RealKind), PARAMETER :: Rv         = 461.00d0   ! gas constant for vapor 
  REAL(RealKind), PARAMETER :: Cvv        = 1424.0d0   ! specific heat of vapor at constant volume Cvv=Cpv-Rv

  REAL(RealKind), PARAMETER :: Cpl        = 4186.0d0   ! specific heat of liquid water at constant pressure 
  REAL(RealKind), PARAMETER :: Cpi        = 2110.0d0   ! specific heat of ice at constant pressure 

  REAL(RealKind), PARAMETER :: L0         = 2.5000d6   ! latent heat of vaporisation at 0 Celsius
  REAL(RealKind)            :: L00                     ! latent heat of vaporisation at 0 Kelvin (=L0+(Cpv-Ct)*T_freeze
  REAL(RealKind), PARAMETER :: Lhs        = 2.83458d6  ! latent heat of sublimation
  REAL(RealKind), PARAMETER :: es0        = 611.d0     ! saturation vapor pressure at tes0
  REAL(RealKind), PARAMETER :: p0Star     = 610.7d0    ! saturation vapor pressure at T_freeze
  REAL(RealKind), PARAMETER :: p0StarIce  = 610.6d0    ! saturation vapor pressure at T_freeze
  REAL(RealKind), PARAMETER :: tes0       = 273.15d0   ! temperature value for es0
  REAL(RealKind), PARAMETER :: p0         = 1.0d5      ! reference pressure
  REAL(RealKind), PARAMETER :: T_freeze   = 273.15d0   ! freezing point temperature, reference temperature

  REAL(RealKind), PARAMETER :: e0         = 6.1078d2   ! saturation water vapor pressure at 0 degree celcius
  REAL(RealKind)            :: gamma      = 1.4d0      ! =Cp/Cv
  REAL(RealKind)            :: kappa      != 0.286d0    ! =Rd/Cp=(Cp-Cv)/Cp

  REAL(RealKind)            :: kappa_liq  = 4.02d-10   ! kappa for water
  REAL(RealKind)            :: alpha_liq  = 4.02d-4    ! alpha for water
  REAL(RealKind)            :: t0_liq     = 293.15d0   ! t0 for water
  REAL(RealKind)            :: rho0_liq   = 998.2d0    ! Rho0 for water
  REAL(RealKind)            :: p0_liq     = 1013.25d2  ! Rho0 for water
  REAL(RealKind)            :: p0_tam     = 3.0d8      ! p0 for water Tammann
  REAL(RealKind)            :: gamma_tam  = 7.15d0     ! gamma for water Tammann
  REAL(RealKind)            :: Cpl_tam  = 1.1626016260162602d0       ! Cpl for water Tammann
  REAL(RealKind)            :: alphaBous  = 0.294d-3    ! Expansion coefficient

  REAL(RealKind), PARAMETER :: Cpair      = 1004.67d0  ! specific heat for air at constant pressure (J/kg*K)
  REAL(RealKind), PARAMETER :: Cvair      = 717.63d0   ! specific heat for air at constant volume (J/kg*K)
  REAL(RealKind), PARAMETER :: Cpwv       = 1865.1d0   ! specific heat for water vapour at constant pressure (J/kg*K)
  REAL(RealKind), PARAMETER :: Cvwv       = 1403.2d0   ! specific heat for water vapour at constant volume (J/kg*K)
  REAL(RealKind), PARAMETER :: CWater     = 4185.5d0   ! specific heat for liquid water (J/kg*K)
  REAL(RealKind), PARAMETER :: CIce       = 2110.d0    ! specific heat for ice (J/kg*K)
  REAL(RealKind), PARAMETER :: LatwaEvap  = 2.501d6    ! latent heat for evaporation and sublimation of water (J/kg)
  REAL(RealKind), PARAMETER :: LatwaSUBL  = 2.834d6    ! latent heat for evaporation and sublimation of water (J/kg)
  REAL(RealKind), PARAMETER :: LatwaFREEZE= 0.334d6    ! latent heat for freezing of water (J/kg)
  REAL(RealKind), PARAMETER :: lv         = 2.50078d6  ! latente W�rme
 
  REAL(RealKind), PARAMETER :: Grav       = 9.81d0     ! Gravitational parameter
  REAL(RealKind)            :: GravComp   = 9.81d0     ! Gravitational parameter
  REAL(RealKind), PARAMETER :: Karm       = 0.4d0      ! von Karman Konstante
  REAL(RealKind), PARAMETER :: Cmy0       = 0.09d0     ! Wirbelviskositaetskoeffizient
  REAL(RealKind), PARAMETER :: Cmy1       = 1.44d0     ! Prod.-term-Konstante in Diss.-Gleichung
  REAL(RealKind), PARAMETER :: Cmy2       = 1.92d0     ! Diss.-term-Konstante in Diss.-Gleichung
  REAL(RealKind), PARAMETER :: Cmy3       = 0.80d0     ! Rich.-term-Konstante in Dis-Gleichung
  REAL(RealKind), PARAMETER :: Cmy4       = 0.93d0     ! Diss.-term-Konstante in TkeSGS-Gleichung
  REAL(RealKind), PARAMETER :: Come1      = 0.52d0     ! Constants from Wilcox (2002)
  REAL(RealKind), PARAMETER :: Come2      = 0.833d0    ! Constants from Wilcox (2002)
  REAL(RealKind), PARAMETER :: Cmiu       = 0.08d0     ! Constants from Wilcox (2002)
  REAL(RealKind), PARAMETER :: SigT       = 0.74d0     ! Rich.-term-Konstante in Tke/Dis-Gleichung
  REAL(RealKind), PARAMETER :: PrandtlDis = 1.3d0      ! Prandtl Zahl for dissipation
  REAL(RealKind), PARAMETER :: CsmagK      = 0.18d0     ! Smagorinsky-Konstante (0.1 - 0.2)

! ---Parameter Turbulence_Mod---
! Constants for tkeHVLen
  REAL(RealKind), PARAMETER :: c0         = 0.32d0
  REAL(RealKind), PARAMETER :: c1         = 0.80d0
  REAL(RealKind), PARAMETER :: c2         = 0.43d0
! Constants for tkeSGS
  REAL(RealKind), PARAMETER :: RichC      = 0.25 ! OSSI 0.25d0     ! Critical Richardson Number
  REAL(RealKind), PARAMETER :: Prn        = 0.7d0      ! Prandtl Number for Neutral Condition
  REAL(RealKind), PARAMETER :: Ceps       = 0.93d0
! REAL(RealKind), PARAMETER :: Cs         = 0.25d0
  REAL(RealKind) :: Cs         = 0.20 ! 0.125d0  !OSSI

! Factor for tke-Production
  REAL(RealKind), PARAMETER :: ftke=0.1d0

  REAL(RealKind)            :: RadEarth   = 6371.220d3 ! Earth radius
  REAL(RealKind)            :: Omega                   ! Winkelgeschwindigkeit
  REAL(RealKind)            :: Curv                    ! Prameter f�r Rotation

  REAL(RealKind)            :: TempAbs    = 293.15d0
  REAL(RealKind)            :: ConvAir                 !ConvAir
  REAL(RealKind)            :: MassSalMin              !MinimalSoluteMass
  REAL(RealKind)            :: VISC_air                !air dynamic viscosity
  REAL(RealKind)            :: KinViscAir              !kinematic viscosity
  REAL(RealKind)            :: TERM_velair             !air thermal velocity
  REAL(RealKind)            :: LAMDA_air               !air mean free path
  REAL(RealKind)            :: BolzMol                 !Boltzmanns constant(kg*m�/(s�*K*Mol))

  REAL(RealKind), PARAMETER :: SBsigma    = 5.670373d-8    ! [W m-2 K-4], Stefan-Boltzmann constant
! REAL(RealKind), PARAMETER :: BOLZ       = 1.380662d-20   !Boltzmanns constant(g*m�/(s�*K*moecule))
  REAL(RealKind), PARAMETER :: BOLZ       = 1.380662d-23   !Boltzmanns constant(kg*m�/(s�*K*moecule))
! REAL(RealKind), PARAMETER :: BOLZ       = 1.380662d-16   !Boltzmanns constant(g*cm�/s�*K*moecule)
! REAL(RealKind), PARAMETER :: BOLZ       = 1.380662d-14   !Boltzmanns constant(\mu g*m^2/s�*K*moecule)
  REAL(RealKind), PARAMETER :: AVGA       = 6.02295D+23    !Avogadros constant (molecule/mole)
  REAL(RealKind), PARAMETER :: ELCH       = 1.6021892d-19  !electron charge (C)
  REAL(RealKind), PARAMETER :: ELMA       = 9.109534d-31   !electron mass (kg)
  REAL(RealKind), PARAMETER :: PRMA       = 1.6605655d-27  !proton mass (kg)
  REAL(RealKind), PARAMETER :: CLIGHT     = 2.99792458d+08 !vacuum velocity of light (m/s)
  REAL(RealKind), PARAMETER :: PLANCK     = 6.6260755d-34  !Plancks constant (Js)
  REAL(RealKind), PARAMETER :: GRAVCONST  = 6.672d-11      !gravitation constant (Nm�/kg�)
  REAL(RealKind), PARAMETER :: CELCIUS_0  = 273.15d0       !0 degree celcius (K)
  REAL(RealKind), PARAMETER :: GASCONST   = 8.3144087d0    !gas constant (J/K*mole)
  REAL(RealKind), PARAMETER :: ATM        = 1.01325d0      !standard pressure or 1 atm (Pa)
  REAL(RealKind), PARAMETER :: Hg2Pa      = 133.322d0      !conversion factor from mm Hg -> Pa (Pa/mm)
  REAL(RealKind), PARAMETER :: RHOLUFT25  = 1.23d-3        !air density at 25 deg celsius (g/cm�)
! REAL(RealKind), PARAMETER :: RHOLUFT25  = 1.185d+3       !air density at 25 deg celsius (\mu g/m�)
  REAL(RealKind), PARAMETER :: MOLUFT     = 28.9644d0      !molecular weight of air (g/mole)
  REAL(RealKind), PARAMETER :: MOLECUL    = 4.8096d-23     !mass of an air molecule (moluft/avga :g)
! REAL(RealKind), PARAMETER :: MOLECUL    = 4.8096d-17     !mass of an air molecule (moluft/avga :\mu g)
  REAL(RealKind), PARAMETER :: Tenswater  = 71.97d-3       !water surface tension at 25 deg celcius (kg/s�)
  REAL(RealKind), PARAMETER :: CAL        = 4.1868d0       !conversion factor cal -> J
  REAL(RealKind), PARAMETER :: StandardN  = 2.5471d25      !molecule density at standard conditions (1/m�)
  REAL(RealKind), PARAMETER :: Sigma0     = 7.28D-2        !Oberfl�chenspannung
  REAL(RealKind), PARAMETER :: Rw         = 461.5D0        !Gaskonstante f�r Wasserdampf (R/Mw) (universelle Gaskonstante/molare Masse von Wasser)
  REAL(RealKind), PARAMETER :: Mol2Part   = 6.023E17

  REAL(RealKind), ALLOCATABLE :: LatentHeat(:)
  REAL(RealKind) :: TotalPrecip=0.0d0

! Pseudo
  REAL(RealKind)            :: cS0=340.0d0
  REAL(RealKind)            :: Rho0=1.0d0

! Chemie
  REAL(RealKind), PARAMETER :: mO2=0.2095d0  ! proportion of O2
  REAL(RealKind), PARAMETER :: mN2=0.7809e0  ! proportion of N2
  REAL(RealKind), PARAMETER :: mH2=0.e0      ! proportion of H2 (not used as 



CONTAINS


SUBROUTINE ComputeParameter

  PI          = 2.0d0*ASIN(1.0d0)
  PiHalf      = 0.5d0*Pi
  ConvAir     = One/Mol2Part
  TERM_velair = SQRT(8.0d0*BOLZ*TempAbs/(PI*MOLECUL*1d-3))
  VISC_air    = 1.8325d-5*(416.16d0/(TempAbs+120.d0))*(TempAbs/296.16d0)**1.5d0
  KinViscAir  = VISC_air/(RHOLUFT25*1d3)
  LAMDA_air   = 2.0d0*VISC_air/(RHOLUFT25*1d3*TERM_velair)
  BolzMol     = Bolz*Avga
  Omega       = 2.0d0*Pi/(24.0d0*3600.0d0)
! Omega       = 2*Pi*5.0d0/60.0d0
  Curv        = 1.0d0
! Rd=Cpd-Cvd
  kappa       = Rd/Cpd
  gamma       = Cpd/Cvd
  L00         = L0+(Cpl-Cpv)*Tes0
  Cs=0.15d0
! WRITE(*,*) 'Omega',Omega
 
END SUBROUTINE ComputeParameter

FUNCTION TERM_vel(MOL)

  REAL(RealKind) :: TERM_vel
  REAL(RealKind) :: MOL

 TERM_vel = SQRT(8.0d0*BOLZ*TempAbs/(PI*MOL))

END FUNCTION TERM_vel

FUNCTION SLIP(r)

  REAL(RealKind) :: SLIP
  REAL(RealKind) :: KNUD
  REAL(RealKind) :: r

  KNUD = LAMDA_air/r
  SLIP = 1.0d0+KNUD*(1.249d0+0.42d0*EXP(-0.87d0/KNUD))

END FUNCTION SLIP

FUNCTION DIFF_coef(r)
 
  REAL(RealKind) :: DIFF_coef
  REAL(RealKind) :: r

  DIFF_coef = BOLZ*TempAbs*SLIP(r)/(6.0d0*PI*r*VISC_air)

END FUNCTION DIFF_coef 


END MODULE Parameter_Mod
