MODULE DCMIP_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE

! GlobalOro
  REAL(RealKind) :: VelMax=35.0d0
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: ProfFac=0.0d0
! REAL(RealKind), PARAMETER :: Rho0=1.0d0
  REAL(RealKind), PARAMETER :: th0=293.16d0
  REAL(RealKind), PARAMETER :: D0=1.d-1
  CHARACTER*40 :: Problem

! Advection 1  
  REAL(RealKind) :: tau, &
                    u0Adv, &
                    k0, &
                    omega0, &
                    T0Adv, &
                    ScaleHeight, &
                    RR, &
                    ZZ, &
                    z0, &
                    lam0Adv, &
                    lam1Adv, &
                    phi0Adv, &
                    phi1Adv, &
                    alphaAdv, &
                    lamp,     &
                    phip,     &
                    zp1,      &
                    zp2,      &
                    zp3,      &
                    dzp1,      &
                    dzp2,      &
                    dzp3,      &
                    Rp,       &
                    zTop

  REAL(RealKind) :: z1, &
                    z2, &
                    K,  &
                    w0Adv
  REAL(RealKind) :: Teq, &
                    LapsRate, &
                    ueq, &
                    peq, &
                    c,   &
                    cSchar,  &
                    XS
  REAL(RealKind) :: Lz,  &
                    delta_theta, &
                    lambdac,     &
                    phic,        &
                    dd
  LOGICAL :: Shear                  
                    
! GloablAcoustic
  REAL(RealKind), PARAMETER :: T0=300.0d0
  REAL(RealKind), PARAMETER :: DeltaP=100.0d0 !Pa
  REAL(RealKind) :: zT
  REAL(RealKind) :: R
  REAL(RealKind) :: phi0=0.0d0
  REAL(RealKind) :: lam0=0.0d0

! GlobalGravity
  REAL(RealKind), PARAMETER :: DeltaTh=10.0d0
  REAL(RealKind) :: OmegaLoc=0.0d0
  REAL(RealKind) :: CurvLoc=0.0d0
  REAL(RealKind) :: nv=1.0d0

! BaroIn / TestCases jw
  REAL(RealKind)            :: exponent1
  REAL(RealKind)            :: a_omega 
  REAL(RealKind)            :: GeoPot0 
  REAL(RealKind)            :: deg2rad
  REAL(RealKind)            :: rad2deg
  REAL(RealKind)            :: rotation_angle
  REAL(RealKind)            :: OmegaFac=1.0d0
  REAL(RealKind), PARAMETER :: delta_T     = 480000.d0  ! in K, for T mean calculation
  REAL(RealKind), PARAMETER :: eta_tropo   = 0.2d0      ! tropopause level
  REAL(RealKind), PARAMETER :: Temp0       = 288.d0     ! horizontal mean T at surface
  REAL(RealKind), PARAMETER :: eta0        = 0.252d0    ! center of jets (hybrid)
  REAL(RealKind), PARAMETER :: perturbation_amplitude =  1.0d0 ! amplitude of u perturbation 1 m/s
  REAL(RealKind), PARAMETER :: radius      = 10.0d0     ! reciprocal radius of the perturbation without 'a'
  REAL(RealKind), PARAMETER :: LapRate     = 0.005d0  ! lapse rate

! DCMIP_4_2
  REAL(RealKind), PARAMETER :: q0=21.0d-3      ! Kg/kg
  REAL(RealKind), PARAMETER :: pw=340.d2     ! Pa
  REAL(RealKind)            :: phiw

! DCMIP_4_1
  REAL(RealKind), PARAMETER :: GammaEPV=0.005d0   ! K/m

! TestCases 3 4 5 6^
  REAL(RealKind), PARAMETER :: MoistLapRate= 0.0065d0  ! lapse rate for moist air
  REAL(RealKind), PARAMETER :: p_ref   = 95500.d0      ! reference pressure
  REAL(RealKind)            :: eta_top
  LOGICAL                   :: Perturbation
  INTEGER                   :: ScalingFactor=1

! Shallow
  REAL(RealKind), PARAMETER :: H06=8000.0d0    ! m
  REAL(RealKind)            :: omega6         ! s^-1 
  REAL(RealKind)            :: K6             ! s^-1
  REAL(RealKind), PARAMETER :: R6=4.0d0       ! wave number 

! MoistVortex
  REAL(RealKind), PARAMETER :: zTP=15000.0d0     ! m
  REAL(RealKind), PARAMETER :: zq1Vo=3000.0d0    ! m
  REAL(RealKind), PARAMETER :: zq2Vo=8000.0d0    ! m
  REAL(RealKind), PARAMETER :: q0Vo=21.0d-3      ! Kg/kg
  REAL(RealKind), PARAMETER :: qTP=1.d-11        ! kg/kg
  REAL(RealKind), PARAMETER :: GammaVo=0.007d0   ! K/m
  REAL(RealKind), PARAMETER :: T0Vo=302.15d0     ! K
  REAL(realKind), PARAMETER :: zpVo=7000.0d0     ! m
  REAL(RealKind)            :: rpVo=282.0d3      ! m
  REAL(RealKind), PARAMETER :: p0Vo=101500.0d0   ! Pascal
  REAL(RealKind)            :: phicVo=10.0d0     ! Grad
  REAL(RealKind)            :: lamcVo=180.0d0    ! Grad
  REAL(RealKind)            :: DeltapVo=1115.0d0 ! Pascal   
  NAMELIST /Example/    &
 
                    VelMax , &
                    ProfFac , &
                    Problem , &
                    Perturbation, &
                    OmegaFac, &
                    N, &
                    phi0, &
                    lam0, &
                    nv, &
                    rpVo, &
                    DeltapVo, &
                    ScalingFactor

CONTAINS

FUNCTION eta_n1(lam,phi,z,U0)
  REAL(RealKind) :: eta_n1
  REAL(RealKind) :: lam,phi,z,U0

  REAL(RealKind) :: eta
  REAL(RealKind) :: Accuracy
  Accuracy=1.d-8
  eta=exp(-Grav*z/(Rd*Temp0))

  Newtonverfahren: DO
    eta_n1=eta-F_lam_phi_eta(lam,phi,eta,z,U0)/dFdeta(lam,phi,eta,U0)
    IF (Accuracy.GE.ABS(eta_n1-eta)) THEN
      EXIT Newtonverfahren 
    ELSE 
      eta=eta_n1
      CYCLE Newtonverfahren 
    END IF
  END DO Newtonverfahren
END FUNCTION eta_n1

FUNCTION F_lam_phi_eta(lam,phi,eta,z,U0)
  REAL(RealKind) :: F_lam_phi_eta
  REAL(RealKind) :: lam,phi,eta,U0,z
  F_lam_phi_eta=-Grav*z+geopotential(lam,phi,eta,U0)-GeoPot0
END FUNCTION F_lam_phi_eta

FUNCTION dFdeta(lam,phi,eta,U0)
  REAL(RealKind) :: dFdeta
  REAL(RealKind) :: lam,phi,eta,U0
  dFdeta=-Rd/eta*temperature(lam,phi,eta,U0)
END FUNCTION dFdeta

FUNCTION pressure(lam,phi,z,U0) 
  REAL(RealKind) :: lam,phi,z
  REAL(RealKind) :: U0
  REAL(RealKind) :: pressure
  REAL(RealKind) :: surface_pressure,phi_perturb
  REAL(RealKind) :: geo_deviation
  REAL(RealKind) :: A,B,C
  A = Half*omega6*(Two*Omega+omega6)*COS(phi)*COS(phi) &
      +0.25d0*K6*K6*(COS(phi))**(Two*R6)*((R6+One)*COS(phi)*COS(phi) &
      +(Two*R6*R6-R6-Two)-Two*R6*R6*(COS(phi))**(-Two))
  B = Two*(Omega+omega6)*K6/(R6+One)/(R6+Two) &
      *(COS(phi))**R6*((R6*R6+Two*R6+Two) &
      -(R6+One)*(R6+One)*COS(phi)*COS(phi))
  C = 0.25d0*K6*K6*(COS(phi))**(Two*R6)*((R6+One)*COS(phi)*COS(phi) &
      -(R6+Two))
  phi_perturb      = RadEarth**2.0d0*(A+B*COS(R6*lam)+C*COS(2.0d0*R6*lam))
  surface_pressure = p_ref*(1.0d0+MoistLapRate/(Grav*Temp0)*phi_perturb)  &
                            **(Grav/(MoistLapRate*Rd))   ! surface pressure
  pressure         = p_ref*((surface_pressure/p_ref)**((MoistLapRate*Rd)/Grav) & 
                               -MoistLapRate/Temp0*z)**(Grav/(MoistLapRate*Rd))
END FUNCTION pressure

!Moist Vortex
FUNCTION PressurePrime(r,z)

  REAL(RealKind) :: PressurePrime
  REAL(RealKind) :: r,z

  REAL(RealKind) :: Tv,qv

  IF(z<=zT) THEN 
    Tv=T0Vo*(1.0d0+0.608d0*q0Vo)-GammaVo*z
    PressurePrime=-DeltapVo &
                 *EXP(-(r/rpVo)**(3.0d0/2.0d0))*EXP(-(z/zpVo)**2.0d0) &
                 *(Tv/(T0Vo*(1.0d0+0.608d0*q0Vo)))**(Grav/(Rd*GammaVo))
  ELSE               
    PressurePrime=Zero
  END IF               
END FUNCTION PressurePrime

FUNCTION TemperaturePrime(r,z)

  REAL(RealKind) :: TemperaturePrime
  REAL(RealKind) :: r,z

  REAL(RealKind) :: Tv,qv

  IF(z<=zT) THEN 
    Tv=T0Vo*(1.0d0+0.608d0*q0Vo)-GammaVo*z
    qv=q0Vo*EXP(-z/zq1Vo)*EXP(-(z/zq2Vo)**2.0d0)
    TemperaturePrime=Tv/(One+qv)*(-Two*Rd*Tv*z)/(Two*Rd*Tv*z &
                    +Grav*zpVo**Two*(One-p0Vo/DeltapVo &
                    *EXP((r/rpVo)**(3.0d0/2.0d0))*EXP((z/zpVo)**2.0d0)))
  ELSE               
    TemperaturePrime=Zero
  END IF               
END FUNCTION TemperaturePrime

FUNCTION temperature(lam,phi,eta,U0) 
  REAL(RealKind) :: temperature
  REAL(RealKind) :: eta,U0
  REAL(RealKind) :: lam,phi

  temperature = t_mean(eta) + t_deviation(lam,phi,eta,U0)
END FUNCTION temperature
  !
  ! Horizontally averaged temperature (equation (4) and (5) in Jablonowski and Williamson (2006))
  !
FUNCTION t_mean(eta)
  REAL(RealKind) :: t_mean
  REAL(RealKind) :: eta

  IF ((eta).GE.(eta_tropo)) THEN
     t_mean = Temp0*eta**exponent1                                ! mean temperature at each level (troposphere)
  ELSE
     t_mean = Temp0*eta**exponent1 + delta_T*(eta_tropo-eta)**5  ! mean temperature at each level (stratosphere)
  ENDIF
END FUNCTION t_mean
  !
  ! Temperature deviation from the horizontal mean 
  ! (equation (6) minus horizontally averaged temperature)
  !
FUNCTION t_deviation(lam,phi,eta,U0)
  REAL(RealKind) :: t_deviation
  REAL(RealKind) :: eta,U0
  REAL(RealKind) :: lam,phi
  REAL(RealKind) :: factor,eta_vertical

  factor       = eta*pi*U0/Rd             
  eta_vertical = (eta - eta0) * 0.5d0*pi
  t_deviation  = factor * 1.5d0 * SIN(eta_vertical) * (cos(eta_vertical))**0.5d0 *                        &
                 ((-2.d0*(SIN(phi))**6 * ((COS(phi))**2 + 1.d0/3.d0) + 10.d0/63.d0)*              &
                 U0 * (COS(eta_vertical))**1.5d0  +                                                      &
                 (8.d0/5.d0*(COS(phi))**3 * ((SIN(phi))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega*0.5d0 )
END FUNCTION t_deviation

FUNCTION geo_deviation(lam,phi,eta,U0)
  REAL(RealKind) :: geo_deviation
  REAL(RealKind) :: eta,U0
  REAL(RealKind) :: lam,phi
  REAL(RealKind) :: A,B,C
  REAL(RealKind) :: cos_tmp

  cos_tmp    = U0 * (COS((eta-eta0)*pi*0.5d0))**1.5d0
  geo_deviation = ( & 
                    ( -2.d0*(SIN(phi))**6.0d0      &
                     * ((COS(phi))**2.0d0 + 1.d0/3.d0) + 10.d0/63.d0 )*cos_tmp     &
                   + ( 8.d0/5.d0*(COS(phi))**3.0d0   &
                       * ((SIN(phi))**2.0d0 + 2.d0/3.d0) - pi/4.d0 )*a_omega      & 
                  )*cos_tmp
END FUNCTION geo_deviation

FUNCTION geo_mean(eta)
  REAL(RealKind) :: geo_mean
  REAL(RealKind) :: factor,eta

      factor=Temp0*Grav/LapRate
      IF (eta.GE.eta_tropo) THEN
        geo_mean = factor*(1.d0-eta**exponent1)
      ELSE
        geo_mean = factor*(1.d0-eta**exponent1)    &
                  - Rd*delta_T*( (LOG(eta/eta_tropo)+137.d0/60.d0)*eta_tropo**5.d0   &
                                - 5.d0*eta_tropo**4.d0*eta                &
                                + 5.d0*eta_tropo**3.d0*eta**2.d0          &
                                - 10.d0/3.d0*eta_tropo**2.d0*eta**3.d0    &
                                + 5.d0/4.d0*eta_tropo*eta**4.d0           &
                                - 1.d0/5.d0*eta**5.d0      )
      END IF
END FUNCTION geo_mean

FUNCTION geopotential(lam,phi,eta,U0)
  REAL(RealKind) :: geopotential
  REAL(RealKind) :: lam,phi,eta,U0,z
  geopotential=geo_mean(eta) + geo_deviation(lam,phi,eta,U0)

END FUNCTION geopotential

FUNCTION u_pert(lam,phi)
  REAL(RealKind) :: u_pert
  REAL(RealKind) :: lam,phi
  REAL(RealKind) :: phi_vertical,sin_tmp,cos_tmp,circdis
  REAL(RealKind) :: perturb_lam,perturb_phi
  REAL(RealKind) :: rot_perturb_lam,rot_perturb_phi
  REAL(RealKind) :: phicVo,lamcVo
 

  perturb_lam = Pi/9.0d0
  perturb_phi = 2.0d0*Pi/9.0d0
  CALL Rotate(perturb_lam,perturb_phi,rot_perturb_lam,rot_perturb_phi)
  sin_tmp = SIN(rot_perturb_phi)*SIN(phi)
  cos_tmp = COS(rot_perturb_phi)*COS(phi)
  circdis = ACOS( sin_tmp + cos_tmp*COS(lam-rot_perturb_lam) )    ! great circle distance without radius 'a'
  u_pert = perturbation_amplitude*EXP(- (circdis*radius)**2.0d0 )
  IF (u_pert <= 1.0d-6) u_pert = 0.0d0
END FUNCTION u_pert

!-----------------------------------------------------------------------
! Ertel's potential vorticity
!-----------------------------------------------------------------------
FUNCTION EPV(lam,phi,eta)
  IMPLICIT NONE
  REAL(RealKind) :: EPV
  REAL(RealKind), INTENT(IN) :: eta,lam,phi
  REAL(RealKind) :: perturb_lam,perturb_phi
  REAL(RealKind) :: rot_perturb_lam,rot_perturb_phi
  REAL(RealKind) :: eta_nu,cos_tmp,Y,circdist,zeta,tmp,Dtmp,F
  REAL(RealKind) :: dudeta,dmeanthetadeta,dthetadeta
  REAL(RealKind) :: dthetadphi1,dthetadphi2,dthetadphi
  REAL(RealKind) :: u0
  REAL(8) :: epsilon

  u0=VelMax
  eta_nu  = (eta-eta0)*pi*0.5d0
  cos_tmp = u0 * (COS(eta_nu))**1.5d0
! initialize epsilon, used to judge the distance of a point to prescribed points
  epsilon = 1.d-6


  Y = (-2.d0*(SIN(phi))**6 * ((COS(phi))**2 + 1.d0/3.d0) + 10.d0/63.d0)*2.d0*COS_tmp   &
      + (8.d0/5.d0*(COS(phi))**3 * ((SIN(phi))**2 + 2.d0/3.d0) - pi/4.d0)*a_omega

  perturb_lam = Pi/9.0d0
  perturb_phi = 2.0d0*Pi/9.0d0
  CALL Rotate(perturb_lam,perturb_phi,rot_perturb_lam,rot_perturb_phi)
  ! great circle distance without radius 'a'
!  K  = SIN(perturb_lat)*SIN(phi) + COS(perturb_lat)*COS(phi)*COS(lam-perturb_lon)
!  DK = SIN(perturb_lat)*COS(phi) - COS(perturb_lat)*SIN(phi)*COS(lam-perturb_lon)
  tmp  = SIN(rot_perturb_phi)*SIN(phi) + COS(rot_perturb_phi)*COS(phi)*COS(lam-rot_perturb_lam)
  Dtmp = SIN(rot_perturb_phi)*COS(phi) - COS(rot_perturb_phi)*SIN(phi)*COS(lam-rot_perturb_lam)
  circdist  = ACOS(tmp)         

  ! relative vorticity
  ! if construct prevents DIVISION BY ZERO in zeta calculation

  IF ((ABS(phi-perturb_phi).LE.epsilon).AND.(ABS(lam-perturb_lam).LE.epsilon)) THEN  ! grid point is at center position
    zeta = -4.d0/RadEarth*cos_tmp*SIN(phi)*COS(phi)*(2.d0-5.d0*(SIN(phi))**2) &
                + perturbation_amplitude/RadEarth*TAN(phi)
  ELSE IF ( ((ABS(phi+perturb_phi).LE.epsilon).AND.(abs(lam-(perturb_lam+Pi)).LE.epsilon)) & ! antipode
            .OR.(ABS(phi-Pi*0.5d0).LE.epsilon)                                             & ! north pole
            .OR.(ABS(phi+Pi*0.5d0).LE.epsilon) ) THEN                                        ! south pole
    zeta = -4.d0/RadEarth*cos_tmp*SIN(phi)*COS(phi)*(2.d0-5.d0*(SIN(phi))**2)
  ELSE                                                                                     ! all other positions
    zeta = -4.d0/RadEarth*cos_tmp*SIN(phi)*COS(phi)*(2.d0-5.d0*(SIN(phi))**2)    &
           + perturbation_amplitude/RadEarth*EXP(- (circdist*radius)**2 )        &
           * (TAN(phi) - 2.d0*radius**2*ACOS(tmp)*Dtmp/SQRT(1.d0-tmp**2))
  END IF
! derivative of u with respect to eta
  dudeta = -u0*(SIN(2.d0*phi))**2*3.d0*pi/4.d0*SQRT(COS(eta_nu))*SIN(eta_nu)

! derivative of mean theta with respect to eta
  dmeanthetadeta = T0*(Rd*GammaEPV/Grav - kappa)*eta**(Rd*GammaEPV/Grav-kappa-1.d0)

  IF (eta < eta_tropo) THEN
    dmeanthetadeta = dmeanthetadeta                          &
                    - delta_T * (                                          &
                      5.d0*(eta_tropo - eta)**4*eta**(-kappa)              &
                      + kappa * (eta_tropo - eta)**5 * eta**(-kappa-1.d0))
  END IF 

  ! derivative of theta with respect to eta
  dthetadeta = dmeanthetadeta                                                             &
               + 3.d0/4.d0*pi*u0/Rd*(1.d0-kappa)*eta**(-kappa)*SIN(eta_nu)*SQRT(COS(eta_nu))*Y       &
               + 3.d0/8.d0*pi*pi/Rd*eta**(1.d0-kappa)*cos_tmp*Y                                      &
               - 3.d0/16.d0*pi*pi*u0/Rd*eta**(1.d0-kappa)*(SIN(eta_nu))**2*(COS(eta_nu))**(-0.5d0)*Y &
               - 9.d0/8.d0*pi*pi*u0*u0/Rd*eta**(1.d0-kappa)*(SIN(eta_nu))**2*COS(eta_nu)             &
                 * (-2.d0*(SIN(phi))**6*((COS(phi))**2+1.d0/3.d0) + 10.d0/63.d0)

  ! derivative of theta with respect to phi
  dthetadphi1 = 2.d0*cos_tmp                                    &
               * (- 12.d0*COS(phi)*(SIN(phi))**5*((COS(phi))**2+1.d0/3.d0) &
                   + 4.d0*COS(phi)*(SIN(phi))**7)

  dthetadphi2 = a_omega                                             &
               * (-24.d0/5.d0*SIN(phi)*(COS(phi))**2*((SIN(phi))**2+2.d0/3.d0) &
                 + 16.d0/5.d0*(COS(phi))**4*SIN(phi))

! corrected on June/6/2012
  dthetadphi = 3.d0/4.d0*pi*u0/Rd*eta**(1.d0-kappa)*SIN(eta_nu)*SQRT(COS(eta_nu)) &
               * (dthetadphi1 + dthetadphi2)

  ! Coriolis parameter
  f = 2.d0 * omega * SIN(phi)

  ! EPV, unscaled since real-Earth 'a' and 'Omega' were used
  EPV = Grav/p0*(-dudeta*dthetadphi/RadEarth - (zeta+f) * dthetadeta)

END FUNCTION EPV

END MODULE DCMIP_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE DCMIP_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

  SELECT CASE(Problem)
    CASE('MoistVortex','DCMIP_5_2')
      BoundCellLoc%TeS=T0Vo
    CASE DEFAULT
      BoundCellLoc%TeS=0.0d0
  END SELECT

END SUBROUTINE SetBoundCells

SUBROUTINE InputExample(FileName)
  USE DCMIP_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&Example')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=Example)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
  Omega=OmegaFac*Omega
  zT=Domain%z1
  R=RadEarth/3.0d0
  deg2rad     = Pi/180.0d0
  rad2deg     = 180.0d0/Pi
  a_omega     = RadEarth*Omega 
  exponent1   = Rd*Laprate/Grav 
  rotation_angle = RotAngle
  lam0         = lam0*deg2rad
  phi0         = phi0*deg2rad
  phicVo         =(phicVo/180.0d0)*pi
  lamcVo         =(lamcVo/180.0d0)*pi
  eta_top=EXP(-Grav*H06/(Rd*Temp0))
  omega6=50/(R6*RadEarth)
  K6=50/(R6*RadEarth)     !K6=omega6=U0/(R6*RadEarth), R6=4=wave number, U0=50m/s
  GeoPot0=geopotential(0.0d0,Half*Pi,1.0d0,VelMax)
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1')
      tau=12.0d0*86400.0d0        ! period of motion 12 days
      u0Adv=(2.0d0*pi*RadEarth)/tau         ! 2 pi a / 12 days
      k0= (10.0d0*RadEarth)/tau          ! Velocity Magnitude
      omega0=(23000.0d0*pi)/tau   ! Velocity Magnitude
      T0Adv=300.0d0                  ! temperature
      ScaleHeight=Rd*T0/Grav                ! scale height
      RR=1.0d0/2.0d0              ! horizontal half width divided by 'a'
      ZZ=1000.0d0                 ! vertical half width
      z0=5000.0d0                 ! center point in z 
      lam0Adv=5.0d0*pi/6.0d0      ! center point in longitudes
      lam1Adv=7.0d0*pi/6.0d0      ! center point in longitudes
      phi0Adv=0.0d0                  ! center point in latitudes
      phi1Adv=0.0d0                  ! center point in latitudes
      ztop=12000.0d0              ! model top         
      z1=domain%z0
      z2=domain%z1
    CASE ('DCMIP_1_2')
      tau=1.0d0*86400.0d0        ! period of motion 1 day
      u0Adv=40.0d0               ! 2 pi a / 12 days
      w0Adv=0.25d0               ! 2 pi a / 12 days
      ScaleHeight=Rd*T0/Grav                ! scale height
      T0Adv=300.0d0                  ! temperature
      K= 5.0d0                    ! number of hadley cells
      z1=3500.0d0                 ! lower tracer bound         
      z2=6500.0d0                 ! upper tracer bound         
      z0=0.5d0*(z1+z2)              ! midpoint        
      ztop=12000.0d0              ! model top         
    CASE ('DCMIP_1_3')
      tau=12.0d0*86400.0d0        ! period of motion 1 day
      u0Adv=2.0d0*Pi*RadEarth/tau
      ScaleHeight=Rd*T0/Grav                ! scale height
      T0Adv=300.0d0                  ! temperature
      ztop=12000.0d0              ! model top         
      alphaAdv=Pi/6.0d0
      lamP=Pi/2.0d0
      phiP=0.0d0
      zp1=3050.0d0
      zp2=5050.0d0
      zp3=8200.0d0
      dzp1=1000.0d0
      dzp2=1000.0d0
      dzp3=400.0d0
      Rp=Pi/4.0d0
    CASE ('DCMIP_2_0')  
      LapsRate=0.0065d0           ! Lapse rate in K/m
      Teq=300.0d0 
    CASE ('DCMIP_2_1')
      XS=500.0d0                  ! Reduced Earth reduction factor
!     Om=0.0d0                    ! Rotation Rate of Earth
      RadEarth=RadEarth/XS        ! New Radius of small Earth     
      ueq=20.0d0                  ! Reference Velocity 
      Teq=300.0d0                 ! Temperature at Equator    
      Peq=100000.0d0              ! Reference PS at Equator
      ztop=30000.0d0              ! Model Top       
!     lambdac=pi/4.0d0            ! Lon of Schar Mountain Center
!     phic=0.0d0                  ! Lat of Schar Mountain Center
!     h0=250.0d0                  ! Height of Mountain
!     d=5000.0d0                  ! Mountain Half-Width
!     xi=4000.0d0                 ! Mountain Wavelength
      c=0.0d0
    CASE ('DCMIP_2_2')
!      Shear=.TRUE.
      XS=500.0d0                  ! Reduced Earth reduction factor
!     Om=0.0d0                    ! Rotation Rate of Earth
      RadEarth=RadEarth/XS        ! New Radius of small Earth     
      ueq=20.0d0                  ! Reference Velocity 
      Teq=300.0d0                 ! Temperature at Equator    
      Peq=100000.0d0              ! Reference PS at Equator
      ztop=30000.0d0              ! Model Top       
!     lambdac=pi/4.0d0            ! Lon of Schar Mountain Center
!     phic=0.0d0                  ! Lat of Schar Mountain Center
!     h0=250.0d0                  ! Height of Mountain
!     d=5000.0d0                  ! Mountain Half-Width
!     xi=4000.0d0                 ! Mountain Wavelength
      cSchar=0.00025d0
      c=cSchar
!      IF (Shear) THEN
!        c=cSchar
!      ELSE  
!        c=0.0d0
!      END IF  
    CASE ('DCMIP_3')
      Omega=0.0d0                 ! Rotation Rate of Earth
      XS=125.0d0
      RadEarth=RadEarth/XS        ! New Radius of small Earth     
      VelMax=20.0d0
      delta_theta=1.d0            ! Max Amplitude of Pert
      Lz=20000.d0          ! Vertical Wavelength of Pert
      N=0.01d0
      lambdac=Pi/4.0d0
      phic=0.0d0
      dd=5000.d0
      Peq=100000.0d0              ! Reference PS at Equator
      Teq=300.0d0                 ! Temperature at Equator    
    CASE ('DCMIP_4_1')
      VelMax=35.0d0
      Omega=Omega*ScalingFactor
      RadEarth=RadEarth/ScalingFactor
      a_omega=RadEarth*Omega
    CASE ('DCMIP_4_2')
      phiw=2.d0*Pi/9.d0
      VelMax=35.0d0
  END SELECT  

END SUBROUTINE InputExample


FUNCTION RhoFun(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G,F,U0
  REAL(RealKind) :: pLoc,rad
  REAL(RealKind) :: eta
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: T,Tv0,Tv,TvWirbel,Tges,Tprime
  REAL(RealKind) :: TvLoc,ThVLoc
  REAL(RealKind) :: p,pt,pprime,pges
  REAL(RealKind) :: qV,qvLoc
  REAL(RealKind) :: RhoD,RhoV
  REAL(RealKind) :: TLoc
  REAL(RealKind) :: ps,TS,N2,bigG

  U0=VelMax
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1','DCMIP_1_2','DCMIP_1_3')
      pLoc =p0*EXP(-z/ScaleHeight)
      tLoc=t0Adv
      RhoFun =pLoc/(Rd*tLoc)
    CASE ('DCMIP_2_0')
      tLoc=Teq-LapsRate*z
      pLoc=p0*(1.0d0-LapsRate/Teq*z)**(Grav/(Rd*LapsRate))
      RhoFun =pLoc/(Rd*tLoc)
    CASE ('DCMIP_2_1','DCMIP_2_2')
      tLoc=Teq*(1.0d0-(c*ueq*ueq/(Grav))*(SIN(phi)**2) )
      pLoc=peq*EXP(-(ueq*ueq/(2.0d0*Rd*Teq))*(SIN(phi)**2)-Grav*z/(Rd*tLoc))
      RhoFun =pLoc/(Rd*tLoc)
    CASE ('DCMIP_3')
      N2=N*N                  
      bigG=(Grav*Grav)/(N2*cpD)
      TS=bigG+(Teq-bigG)*EXP(-(U0*N2/(4.0d0*Grav*Grav))*(U0+2.0d0*Omega*RadEarth)*(COS(2.d0*phi)-1.0d0))
      ps=peq*EXP((U0/(4.0d0*bigG*Rd))*(U0+2.0*Omega*RadEarth)*(COS(2.0d0*phi)-1.0d0)) &
        *(TS/Teq)**(cpD/Rd)
      pLoc=ps*((bigG/TS)*EXP(-N2*z/Grav)+1.0-(bigG/TS))**(cpD/Rd)
      TLoc=bigG*(1.0d0-EXP(N2*z/Grav))+TS*EXP(N2*z/Grav)
      RhoFun =pLoc/(Rd*tLoc)
    CASE('GlobalOro')
      G=Grav*Grav/(N*N*cpD)
      S=N*N/Grav
      F=-U0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*U0)
      IF (N>Zero) THEN
        ExPre=(One-G/th0)+G/th0*EXP(-S*(z+F/Grav))
      ELSE
        ExPre=One-Grav/(th0*cpD)*(z+F/Grav)
      END IF
      pLoc=p0*ExPre**(1.0d0/kappa)
      ThLoc=th0*exp(S*(z+F/Grav))
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
    CASE('GlobalAcoustic')
      ploc=p0*EXP(-z/(Rd*t0))
      rad=RadEarth*ACOS(SIN(phi0)*SIN(phi)+COS(phi0)*COS(phi)*COS(lam-lam0))
      F=Zero
      IF (rad<R) THEN
        F=Half*(One+COS(Pi*rad/R))
      END IF
      g=SIN(Pi*z/zT)
      pLoc=ploc+DeltaP*F*g
      RhoFun=pLoc/(Rd*t0)
    CASE('GlobalGravity')
      G=Grav*Grav/(N*N*Cpd)
      S=N*N/Grav
      F=-U0*COS(phi)**2*(RadEarth*OmegaLoc+CurvLoc*0.5d0*U0)
      IF (N>Zero) THEN
        ExPre=(One-G/t0)+G/t0*EXP(-S*(z+F/Grav))
      ELSE
        ExPre=One-Grav/(t0*Cpd)*(z+F/Grav)
      END IF
      pLoc=p0*ExPre**(1.0d0/kappa)
      ThLoc=t0*exp(S*(z+F/Grav))
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

    CASE('MoistVortex','DCMIP_5_2')
      r=RadEarth*(ACOS(SIN(phicVo)*SIN(phi)+COS(phicVo)*COS(phi)*COS(lam-lamcVo)))
      Tv0=T0Vo*(1.0d0+0.608d0*q0Vo)
      IF(z<=zT) THEN 
        qV=q0Vo*EXP(-z/zq1Vo)*EXP(-(z/zq2Vo)**2.0d0)
        Tv=Tv0-GammaVo*z
        T=Tv/(1.0d0+0.608d0*qV)
        p=p0Vo*(Tv/Tv0)**(Grav/(Rd*GammaVo))
      ELSE  IF(z>zT) THEN
        qV=qTP
        Tv=Tv0-GammaVo*zTP
        T=Tv
        pt=p0Vo*((Tv/Tv0)**(Grav/(Rd*GammaVo))) 
        p=pt*EXP(-((Grav*(z-zTP))/(Rd*Tv)))
      END IF
      Tges=T+TemperaturePrime(r,z)
      pges=p+PressurePrime(r,z)
      RhoFun=pGes/(TGes*((Rv-Rd)*qV+Rd))
    CASE('BaroIn','DCMIP_4_1')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      eta=eta_n1(rot_lam,rot_phi,z,ABS(U0))
      ThLoc=temperature(rot_lam,rot_phi,eta,ABS(U0)) 
      RhoFun=p0*eta/(Rd*ThLoc)
    CASE ('DCMIP_4_2')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      eta=eta_n1(rot_lam,rot_phi,z,ABS(U0))
      qvLoc=q0*EXP(-(rot_phi/phiw)**4.d0)*EXP(-((eta-1.d0)*p0/pw)**2.d0)
      ThVLoc=temperature(rot_lam,rot_phi,eta,ABS(U0))
      pLoc=p0*eta
      ThLoc=ThVLoc !/(1.d0+0.608d0*qvLoc)
      RhoFun=pLoc/(Rd*ThVLoc) !(ThLoc*((Rv-Rd)*qVLoc+Rd))
    CASE ('William6')
      pLoc=pressure(lam,phi,z,U0)
      ThLoc=Temp0*(pLoc/p_ref)**((MoistLapRate*Rd)/Grav)
      RhoFun=pLoc/(Rd*ThLoc)
  END SELECT
END FUNCTION RhoFun

FUNCTION RhoProf(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G,F,U0,pLoc
  REAL(RealKind) :: T,Tv0,Tv,TvWirbel,Tges,Tprime
  REAL(RealKind) :: p,pt,pprime,pges
  REAL(RealKind) :: qV
  REAL(RealKind) :: RhoD,RhoV
  U0=VelMax
  SELECT CASE(Problem)
    CASE('GlobalOro')
      G=Grav*Grav/(N*N*cpD)
      S=N*N/Grav
      F=-U0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*U0)
      IF (N>Zero) THEN
        ExPre=(One-G/th0)+G/th0*EXP(-S*(z+F/Grav))
      ELSE
        ExPre=One-Grav/(th0*cpD)*(z+F/Grav)
      END IF
      pLoc=p0*ExPre**(1.0d0/kappa)
      ThLoc=th0*exp(S*(z+F/Grav))
      RhoProf=ProfFac*pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
    CASE('GlobalAcoustic')
      ploc=p0*EXP(-z/(Rd*t0))
      RhoProf=pLoc/(Rd*t0)
    CASE('GlobalGravity')
      G=Grav*Grav/(N*N*Cpd)
      S=N*N/Grav
      IF (N>Zero) THEN
        ExPre=(One-G/t0)+G/t0*EXP(-S*z)
      ELSE
        ExPre=One-Grav/(t0*Cpd)*z
      END IF
      pLoc=p0*ExPre**(1.0d0/kappa)
      ThLoc=t0*exp(S*z)
      RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
    CASE('MoistVortex','DCMIP_5_2')
      r=RadEarth*(ACOS(SIN(phicVo)*SIN(phi)+COS(phicVo)*COS(phi)*COS(lam-lamcVo)))
      Tv0=T0Vo*(1.0d0+0.608d0*q0Vo)
      IF(z<=zT) THEN 
        qV=q0Vo*EXP(-z/zq1Vo)*EXP(-(z/zq2Vo)**2.0d0)
        Tv=Tv0-GammaVo*z
        T=Tv/(1.0d0+0.608d0*qV)
        p=p0Vo*(Tv/Tv0)**(Grav/(Rd*GammaVo))
      ELSE  IF(z>zT) THEN
        qV=qTP
        Tv=Tv0-GammaVo*zTP
        T=Tv
        pt=p0Vo*((Tv/Tv0)**(Grav/(Rd*GammaVo))) 
        p=pt*EXP(-((Grav*(z-zTP))/(Rd*Tv)))
      END IF
      Tges=T
      pges=p
      RhoProf=pGes/(TGes*((Rv-Rd)*qV+Rd))
    CASE DEFAULT
      RhoProf=0.0d0
  END SELECT
END FUNCTION RhoProf

FUNCTION PreStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G,F,U0
  U0=VelMax
  SELECT CASE(Problem)
    CASE('GlobalOro')
      G=Grav*Grav/(N*N*cpD)
      S=N*N/Grav
      F=-U0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*U0)
      IF (N>Zero) THEN
        ExPre=(One-G/th0)+G/th0*EXP(-S*(z+F/Grav))
      ELSE
        ExPre=One-Grav/(th0*cpD)*(z+F/Grav)
      END IF
      PreStart=p0*ExPre**(1.0d0/kappa)
    CASE('GlobalAcoustic')
      PreStart=p0*EXP(-z/(Rd*t0))
    CASE('GlobalGravity')
      G=Grav*Grav/(N*N*Cpd)
      S=N*N/Grav
      F=-U0*COS(phi)**2*(RadEarth*OmegaLoc+CurvLoc*0.5d0*U0)
      IF (N>Zero) THEN
        ExPre=(One-G/t0)+G/t0*EXP(-S*(z+F/Grav))
      ELSE
        ExPre=One-Grav/(t0*Cpd)*(z+F/Grav)
      END IF
      PreStart=p0*ExPre**(1.0d0/kappa)
    CASE DEFAULT
      PreStart=0.0d0
  END SELECT
END FUNCTION PreStart

FUNCTION PreProf(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G
  SELECT CASE(Problem)
    CASE('GlobalOro')
      PreProf=0.0d0
    CASE('GlobalAcoustic')
      PreProf=p0*EXP(-z/(Rd*t0))
    CASE('GlobalGravity')
      S=N*N/Grav
      G=Grav*Grav/(N*N*Cpd)
      IF (N>Zero) THEN
        ExPre=(One-G/t0)+G/t0*EXP(-S*z)
      ELSE
        ExPre=One-Grav/(t0*Cpd)*z
      END IF
      PreProf=p0*ExPre**(1.0d0/kappa)
    CASE DEFAULT
      PreProf=0.0d0
  END SELECT
END FUNCTION PreProf

FUNCTION ThProfFun(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,F,U0
  REAL(RealKind) :: TV0,TV,qV
  REAL(RealKind) :: Rho,RhoV,RhoD
  REAL(RealKind) :: KappaLoc
  REAL(RealKind) :: T,p,pt
  REAL(RealKind) :: d1,d2,d
  SELECT CASE(Problem)
    CASE('GlobalOro')
      U0=VelMax
      F=-U0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*U0)
      S=N*N/Grav
      ThProfFun=ProfFac*th0*exp(S*(z+F/Grav))
    CASE('GlobalAcoustic')
      ThProfFun=t0*EXP(kappa*z/(Rd*t0))
    CASE('GlobalGravity')
      S=N*N/Grav
      ThProfFun=t0*EXP(S*z)
    CASE ('MoistVortex','DCMIP_5_2')
      Tv0=T0Vo*(1+0.608d0*q0Vo)
      IF(z<=zT) THEN
        qV=q0Vo*EXP(-z/zq1Vo)*EXP(-(z/zq2Vo)**2.0d0)
        Tv=Tv0-GammaVo*z
        T=Tv/(1.0d0+0.608d0*qV)
        p=p0Vo*(((Tv0-GammaVo*z)/Tv0)**(Grav/(Rd*GammaVo)))
      ELSE  IF(z>zT) THEN
        qV=qTP
        Tv=Tv0-GammaVo*zTP
        T=Tv
        pt=p0Vo*((Tv/Tv0)**(Grav/(Rd*GammaVo)))
        p=pt*EXP(-((Grav*(z-zTP))/(Rd*Tv)))
      END IF
      Rho=p/(T*((Rv-Rd)*qV+Rd))
      RhoV=Rho*qV
      RhoD=Rho-RhoV
      KappaLoc=(Rd*RhoD+Rv*RhoV)/(cpD*RhoD+cpV*RhoV)
      ThProfFun=(p0Vo/p)**KappaLoc*(RhoD+Rv/Rd*RhoV)/(RhoD+RhoV)*T
    CASE DEFAULT
      ThProfFun=0.0d0
  END SELECT
END FUNCTION ThProfFun

FUNCTION QvProfFun(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: RhoV, RhoFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  SELECT CASE(Problem)
    CASE ('MoistVortex','DCMIP_5_2')
      r=RadEarth*(ACOS(SIN(phicVo)*SIN(phi)+COS(phicVo)*COS(phi)*COS(lam-lamcVo)))
      Tv0=T0Vo*(1.0d0+0.608d0*q0Vo)
      IF(z<=zT) THEN 
        qV=q0Vo*EXP(-z/zq1Vo)*EXP(-(z/zq2Vo)**2.0d0)
        Tv=Tv0-GammaVo*z

        TvWirbel=(Tv0-GammaVo*z)/((1.0d0+((2*Rd*(Tv0-GammaVo*z)*z)/(Grav*(zpVo**2)&
                *(1.0d0-(p0Vo/DeltapVo)*EXP((r/rpVo)**(3.0d0/2.0d0))*EXP((z/zpVo)**2.0d0))))))
!       T'=Temperatureprime
!       p'=pressureprime
        T=Tv/(1.0d0+0.608d0*qV)
        p=p0Vo*(((Tv0-GammaVo*z)/Tv0)**(Grav/Rd*GammaVo))
        Tges=T
        pges=p
      ELSE  IF(z>zT) THEN
        qV=qTP
        Tv=Tv0-GammaVo*zTP
        TvWirbel=Tv
        T=Tv
        Tprime=0 
        pt=p0Vo*((Tv/Tv0)**(Grav/(Rd*GammaVo))) 
        p=pt*EXP(-((Grav*(z-zTP))/(Rd*Tv)))
        pprime=0
        Tges=T
        pges=p
      END IF
      RhoD=(p*(1.0d0-qV))/(T*(Rd*(1.0d0-qV)+qV*Rv))
      RhoV=(p*qV)/(T*(1.0d0-qV)*Rd+Rv*qV)
      RhoFun=RhoD+RhoV
      QvProfFun=RhoV/RhoFun
    CASE DEFAULT
      QvProf=0.d0
  END SELECT
END FUNCTION QvProfFun

FUNCTION UStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  USE Physics_Mod
  USE Rho_Mod
  REAL(RealKind) :: UStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: phiL,phiR
  REAL(RealKind) :: U0,u_lat,v_lat,v_tmp
  REAL(RealKind) :: Tv0
  REAL(RealKind) :: vT
  REAL(RealKind) :: d1,d2,d
  REAL(RealKind) :: fc
  REAL(RealKind) :: lonp,tLoc
  U0=VelMax
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1')
      lonp=lam-2.0d0*Pi*Time/tau
      Ustart=k0*SIN(lonp)*SIN(lonp)*SIN(2.0*phi)*COS(Pi*Time/tau)+u0Adv*COS(phi)
      Ustart=RhoFun(lam,phi,z,zHeight,Time)*UStart
    CASE ('DCMIP_2_0')
      UStart=0.0d0
    CASE ('DCMIP_1_2')
      UStart=RhoFun(lam,phi,z,zHeight,Time)*u0Adv*COS(phi)
    CASE ('DCMIP_1_3')
      UStart=RhoFun(lam,phi,z,zHeight,Time)*u0Adv*(COS(phi)*COS(alphaAdv)+SIN(phi)*COS(lam)*SIN(alphaAdv))
    CASE ('DCMIP_2_1','DCMIP_2_2')
      tLoc=Teq*(1.0d0-(c*ueq*ueq/(Grav))*(SIN(phi)**2) )
      uStart=ueq*COS(phi)*SQRT((2.0d0*Teq/(tLoc))*c*z+tLoc/(Teq))
      Ustart=UStart
    CASE ('DCMIP_3')
      Ustart=U0*COS(phi)
    CASE('GlobalOro')
      UStart=U0*COS(phi)
    CASE('GlobalAcoustic')
      UStart=U0*COS(phi)
    CASE('GlobalGravity')
      UStart=U0*COS(phi)
    CASE('BaroIn','DCMIP_4_2','DCMIP_4_1')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      eta=eta_n1(rot_lam,rot_phi,z,ABS(U0))
      eta_v=(eta-eta0)*0.5d0*Pi
      u_lat=U0*COS(eta_v)**1.5d0*SIN(2.0d0*rot_phi)**2.0d0
      IF (Perturbation) THEN
        u_lat=u_lat+u_pert(rot_lam,rot_phi)
      END IF
      IF (ABS(rotation_angle)<1.0E-8) THEN
        UStart = u_lat
      ELSE
        v_lat = 0.0d0
        ! rotate wind components
        CALL turnwi(u_lat,v_lat,UStart,v_tmp,lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*Pi+rotation_angle,-1)
        IF (ABS(UStart)<1.0E-10) UStart=0.0d0
      ENDIF
    CASE ('William6')
      UStart=RadEarth*omega6*COS(phi) &
            +RadEarth*K6*(COS(phi))**(R6-One) &
            *(R6*SIN(phi)*SIN(phi)-COS(phi)*COS(phi))*COS(R6*lam)

    CASE ('MoistVortex','DCMIP_5_2')
      r=RadEarth*(ACOS(SIN(phicVo)*SIN(phi)+COS(phicVo)*COS(phi)*COS(lam-lamcVo)))
      fc=Two*Omega*SIN(phicVo)
      Tv0=T0Vo*(1.0d0+0.608d0*q0Vo)
      IF (z<=zTP) THEN  
!       WRITE(*,*) (fc*r)/2.0d0,-3.0d0/2.0d0 &
!         *(r/rpVo)**(3.0d0/2.0d0)*(Tv0-GammaVo*z)*Rd/(1.0d0+2.0d0*Rd*(Tv0-GammaVo*z)&
!         *z/(Grav*zpVo**2)-p0Vo/DeltapVo*EXP((r/rpVo)**(3.0d0/2.0d0))*EXP((z/zpVo)**2.0d0))
        vT=-(fc*r)/2.0d0+SQRT((fc*r)**2.0d0/4.0d0-3.0d0/2.0d0 &
          *(r/rpVo)**(3.0d0/2.0d0)*(Tv0-GammaVo*z)*Rd/(1.0d0+2.0d0*Rd*(Tv0-GammaVo*z)&
          *z/(Grav*zpVo**2)-p0Vo/DeltapVo*EXP((r/rpVo)**(3.0d0/2.0d0))*EXP((z/zpVo)**2.0d0)))
      ELSE IF (z>zTP) THEN 
        vT=0
      END IF
      d1=SIN(phicVo)*COS(phi)-COS(phicVo)*SIN(phi)*COS(lam-lamcVo)
      d2=COS(phicVo)*SIN(lam-lamcVo)
!     d=MAX(Epsilon, SQRT(d1**2.0d0+d2**2.0d0))
      d=MAX(1.d-40, SQRT(d1**2.0d0+d2**2.0d0))
      UStart=(vT*d1)/d

  END SELECT
END FUNCTION UStart

FUNCTION UStartE(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  USE UVW_Mod
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  UStartE=UStart(lam,phi,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  USE Physics_Mod
  USE Rho_Mod
  REAL(RealKind) :: VStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: phiL,phiR
  REAL(RealKind) :: U0,u_lat,v_lat,u_tmp
  REAL(RealKind) :: Tv0
  REAL(RealKind) :: vT
  REAL(RealKind) :: fc
  REAL(RealKind) :: d1,d2,d
  REAL(RealKind) :: lonp
  U0=VelMax
  VStart=0.0d0
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1')
      lonp=lam-2.0d0*Pi*Time/tau
      VStart=k0*SIN(2.0*lonp)*COS(phi)*COS(Pi*Time/tau)
      Vstart=RhoFun(lam,phi,z,zHeight,Time)*VStart
    CASE ('DCMIP_2_0')
      VStart=0.0d0
    CASE ('DCMIP_1_2')
      VStart=-RadEarth*w0Adv*Pi/(K*zTop)*COS(phi)*SIN(K*phi)*COS(Pi*z/ztop)*COS(Pi*Time/tau)
      Vstart=RhoFun(lam,phi,z,zHeight,Time)*VStart
    CASE ('DCMIP_2_1','DCMIP_2_2')
      VStart=0.0d0
    CASE ('DCMIP_1_3')
      VStart=RhoFun(lam,phi,z,zHeight,Time)*u0Adv*(-SIN(lam)*SIN(alphaAdv))
    CASE ('DCMIP_3')
      VStart=0.0d0
    CASE ('William6')
      VStart=-RadEarth*K6*R6*(COS(phi))**(R6-One)*SIN(phi)*SIN(R6*lam)
    CASE ('BaroIn','DCMIP_4_2','DCMIP_4_1')
      IF (ABS(rotation_angle)<1.0E-8) THEN
        VStart = 0.0d0
      ELSE
        CALL Rotate(lam,phi,rot_lam,rot_phi)
        eta=eta_n1(rot_lam,rot_phi,z,ABS(U0))
        eta_v=(eta-eta0)*0.5d0*Pi
        u_lat=U0*COS(eta_v)**1.5d0*SIN(2.0d0*rot_phi)**2.0d0
        IF (Perturbation) THEN
          u_lat=u_lat+u_pert(rot_lam,rot_phi)
        END IF
        v_lat = 0.0d0
        ! pole point velocities are not well-defined
        IF (ABS(Pi*0.5d0-phi)<1.0E-8.OR.ABS(Pi*0.5d0+phi)<1.0E-8) THEN
          VStart = 0.0d0
        ELSE
          ! rotate wind components
          CALL turnwi(u_lat,v_lat,u_tmp,VStart,lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*Pi+rotation_angle,-1)
        ENDIF
      ENDIF
    CASE ('MoistVortex','DCMIP_5_2')
      r=RADEARTH*(ACOS(SIN(phicVo)*SIN(phi)+COS(phicVo)*COS(phi)*COS(lam-lamcVo)))
      fc=Two*Omega*SIN(phicVo)
      Tv0=T0Vo*(1.0d0+0.608d0*q0Vo)
      IF(z<=zT) THEN 
        vT=-(fc*r)/2.0d0+SQRT((fc*r)**2.0d0/4.0d0-3.0d0/2.0d0 &
          *(r/rpVo)**(3.0d0/2.0d0)*(Tv0-GammaVo*z)*Rd/(1.0d0+2.0d0*Rd*(Tv0-GammaVo*z)&
          *z/(Grav*zpVo**2)-p0Vo/DeltapVo*EXP((r/rpVo)**(3.0d0/2.0d0))*EXP((z/zpVo)**2.0d0)))
      ELSE  IF(z>zT) THEN 
        vT=0.0d0
      END IF
      d1=SIN(phicVo)*COS(phi)-COS(phicVo)*SIN(phi)*COS(lam-lamcVo)
      d2=COS(phicVo)*SIN(lam-lamcVo)
!     d=max(Epsilon, SQRT(d1**2.0d0+d2**2.0d0))
      d=MAX(1.d-40, SQRT(d1**2.0d0+d2**2.0d0))
      VStart=(vT*d2)/d
    CASE DEFAULT
      VStart=0.0d0
  END SELECT
                            
END FUNCTION VStart

FUNCTION VStartE(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  USE UVW_Mod
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  VStartE=VStart(lam,phi,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  USE Rho_Mod
  REAL(RealKind) :: WStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1')
      WStart=0.0d0
    CASE ('DCMIP_1_2')
      WStart=w0Adv/K*(-2.0d0*SIN(K*phi)*SIN(phi)+K*COS(phi)*COS(K*phi))*SIN(Pi*z/ztop)*COS(Pi*Time/tau)
      Wstart=RhoFun(lam,phi,z,zHeight,Time)*WStart
    CASE ('DCMIP_1_3')
      WStart=0.0d0
    CASE ('DCMIP_2_0')
      WStart=0.0d0
    CASE ('DCMIP_2_1','DCMIP_2_2')
      WStart=0.0d0
    CASE ('DCMIP_3')
      WStart=0.0d0
    CASE DEFAULT  
      WStart=0.0d0
  END SELECT    
END FUNCTION WStart

FUNCTION EnStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  USE ThProf_Mod
  REAL(RealKind) :: EnStart
  EnStart=Zero
END FUNCTION EnStart

FUNCTION ThStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,F,ff,gg,rad
  REAL(RealKind) :: U0,eta,T,ThLoc,pLoc,rhoLoc,tLoc
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: TV0,TV,ThVLoc,qV,qvLoc
  REAL(RealKind) :: Rho,RhoV,RhoD
  REAL(RealKind) :: KappaLoc
  REAL(RealKind) :: p,pt
  REAL(RealKind) :: d1,d2,d
  REAL(RealKind) :: N2,bigG,TS,ps
  REAL(RealKind) :: sin_tmp,cos_tmp,ss,rLoc,theta_pert
  U0=VelMax
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1','DCMIP_1_2','DCMIP_1_3')
      pLoc =p0*exp(-z/ScaleHeight)
      tLoc=T0Adv
      rhoLoc =pLoc/(Rd*tLoc)
      ThStart=(p0/pLoc)**(Rd/Cpd)*TLoc
    CASE ('DCMIP_2_0')
      tLoc=Teq-LapsRate*z
      pLoc=p0*(1.0d0-LapsRate/Teq*z)**(Grav/(Rd*LapsRate))
      ThStart=(p0/pLoc)**(Rd/Cpd)*TLoc
    CASE ('DCMIP_2_1','DCMIP_2_2')
      tLoc=Teq*(1.0d0-(c*ueq*ueq/(Grav))*(SIN(phi)**2) )
      pLoc=peq*EXP(-(ueq*ueq/(2.0d0*Rd*Teq))*(SIN(phi)**2)-Grav*z/(Rd*tLoc))
      ThStart=(p0/pLoc)**(Rd/Cpd)*TLoc
    CASE ('DCMIP_3')
      N2=N*N                  
      bigG=(Grav*Grav)/(N2*cpD)
      TS=bigG+(Teq-bigG)*EXP(-(U0*N2/(4.0d0*Grav*Grav))*(U0+2.0d0*Omega*RadEarth)*(COS(2.d0*phi)-1.0d0))
      ps=peq*EXP((U0/(4.0d0*bigG*Rd))*(U0+2.0*Omega*RadEarth)*(COS(2.0d0*phi)-1.0d0)) &
        *(TS/Teq)**(cpD/Rd)
      pLoc=ps*((bigG/TS)*EXP(-N2*z/Grav)+1.0-(bigG/TS))**(cpD/Rd)
      TLoc=bigG*(1.0d0-EXP(N2*z/Grav))+TS*EXP(N2*z/Grav)
      ThStart=(p0/pLoc)**(Rd/Cpd)*TLoc
      sin_tmp=SIN(phi)*SIN(phic)
      cos_tmp=COS(phi)*COS(phic)
! great circle distance with 'a/X' 
      rLoc=RadEarth*ACOS(sin_tmp+cos_tmp*COS(lam-lambdac))
      ss=(dd**2)/(dd**2+rLoc**2)
      theta_pert=delta_theta*ss*SIN(2.d0*Pi*z/Lz)
      ThStart=ThStart+theta_pert
    CASE('GlobalOro')
      F=-U0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*U0)
      S=N*N/Grav
      ThStart=th0*exp(S*(z+F/Grav))
    CASE('GlobalAcoustic')
      ThStart=t0*EXP(kappa*z/(Rd*t0))
    CASE('GlobalGravity')
      F=-U0*COS(phi)**2*(RadEarth*OmegaLoc+CurvLoc*0.5d0*U0)
      S=N*N/Grav
      ThStart=t0*exp(S*(z+F/Grav))
      rad=RadEarth*ACOS(SIN(phi0)*SIN(phi)+COS(phi0)*COS(phi)*COS(lam-lam0))
      ff=Zero
      IF (rad<R) THEN
        ff=Half*(One+COS(Pi*rad/R))
      END IF
      gg=SIN(nv*Pi*z/zT)
      ThStart=ThStart+DeltaTh*ff*gg
    CASE('DCMIP_4_2')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      eta=eta_n1(rot_lam,rot_phi,z,ABS(U0))
      qvLoc=q0*EXP(-(rot_phi/phiw)**4.d0)*EXP(-((eta-1.d0)*p0/pw)**2.d0)
      ThVLoc=temperature(rot_lam,rot_phi,eta,ABS(U0))
      pLoc=p0*eta
      ThLoc=ThVLoc !/(1.d0+0.608d0*qvLoc)
      Rho=pLoc/(Rd*ThVLoc) !(ThLoc*((Rv-Rd)*qVLoc+Rd))
      RhoV=Rho*qVLoc
      RhoD=Rho-RhoV
      KappaLoc=(Rd*RhoD+Rv*RhoV)/(cpD*RhoD+cpV*RhoV)
      ThStart=ThLoc*(p0/pLoc)**KappaLoc     
    CASE('BaroIn','DCMIP_4_1')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      eta=eta_n1(rot_lam,rot_phi,z,ABS(U0))
      T=temperature(rot_lam,rot_phi,eta,ABS(U0)) 
      ThStart=T/eta**kappa
    CASE ('William6')
      pLoc=pressure(lam,phi,z,U0)
      ThLoc=Temp0*(pressure(lam,phi,z,U0)/p_ref)**((MoistLapRate*Rd)/Grav)
      ThStart=ThLoc*(p0/pLoc)**kappa
    CASE ('MoistVortex','DCMIP_5_2')
      Tv0=T0Vo*(1+0.608d0*q0Vo)
      IF(z<=zT) THEN
        qV=q0Vo*EXP(-z/zq1Vo)*EXP(-(z/zq2Vo)**2.0d0)
        Tv=Tv0-GammaVo*z
        T=Tv/(1.0d0+0.608d0*qV)
        p=p0Vo*(((Tv0-GammaVo*z)/Tv0)**(Grav/(Rd*GammaVo)))
      ELSE  IF(z>zT) THEN
        qV=qTP
        Tv=Tv0-GammaVo*zTP
        T=Tv
        pt=p0Vo*((Tv/Tv0)**(Grav/(Rd*GammaVo)))
        p=pt*EXP(-((Grav*(z-zTP))/(Rd*Tv)))
      END IF
      T=T+TemperaturePrime(r,z)
      p=p+PressurePrime(r,z)
      Rho=p/(T*((Rv-Rd)*qV+Rd))
      RhoV=Rho*qV
      RhoD=Rho-RhoV
      KappaLoc=(Rd*RhoD+Rv*RhoV)/(cpD*RhoD+cpV*RhoV)
      ThStart=(p0Vo/p)**KappaLoc*(RhoD+Rv/Rd*RhoV)/(RhoD+RhoV)*T
  END SELECT
END FUNCTION ThStart

FUNCTION TkeStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  TkeStart=0.0d0
END FUNCTION TkeStart

FUNCTION DisStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DisStart=0.0d0
END FUNCTION DisStart

FUNCTION QvStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: RhoV, RhoFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: U0
  U0=VelMax
  SELECT CASE(Problem)
    CASE ('MoistVortex','DCMIP_5_2')
      r=RadEarth*(ACOS(SIN(phicVo)*SIN(phi)+COS(phicVo)*COS(phi)*COS(lam-lamcVo)))
      Tv0=T0Vo*(1.0d0+0.608d0*q0Vo)
      IF(z<=zT) THEN 
        qV=q0Vo*EXP(-z/zq1Vo)*EXP(-(z/zq2Vo)**2.0d0)
        Tv=Tv0-GammaVo*z

        TvWirbel=(Tv0-GammaVo*z)/((1.0d0+((2*Rd*(Tv0-GammaVo*z)*z)/(Grav*(zpVo**2)&
                *(1.0d0-(p0Vo/DeltapVo)*EXP((r/rpVo)**(3.0d0/2.0d0))*EXP((z/zpVo)**2.0d0))))))
!       T'=Temperatureprime
!       p'=pressureprime
        T=Tv/(1.0d0+0.608d0*qV)
        p=p0Vo*(((Tv0-GammaVo*z)/Tv0)**(Grav/Rd*GammaVo))
        Tges=T+TemperaturePrime(r,z)
        pges=p+PressurePrime(r,z)
      ELSE  IF(z>zT) THEN
        qV=qTP
        Tv=Tv0-GammaVo*zTP
        TvWirbel=Tv
        T=Tv
        Tprime=0 
        pt=p0Vo*((Tv/Tv0)**(Grav/(Rd*GammaVo))) 
        p=pt*EXP(-((Grav*(z-zTP))/(Rd*Tv)))
        pprime=0
        Tges=T+Tprime
        pges=p+pprime
      END IF
      RhoD=(p*(1.0d0-qV))/(T*(Rd*(1.0d0-qV)+qV*Rv))
      RhoV=(p*qV)/(T*(1.0d0-qV)*Rd+Rv*qV)
      RhoFun=RhoD+RhoV
      QvStart=RhoV/RhoFun
    CASE('DCMIP_4_2')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      eta=eta_n1(rot_lam,rot_phi,z,ABS(U0))
      QvStart=q0*EXP(-(rot_phi/phiw)**4.d0)*EXP(-((eta-1.d0)*p0/pw)**2.d0)
    CASE DEFAULT
      QvStart=0.d0
  END SELECT
END FUNCTION QvStart

FUNCTION QcStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: DStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: U0,eta_c,eta
  REAL(RealKind) :: sin_tmp,cos_tmp,circdis
  REAL(RealKind) :: sin_tmp2,cos_tmp2
  REAL(RealKind) :: r2,d1,d2
  REAL(RealKind) :: perturb_lam,perturb_phi,tmp
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1')
      sin_tmp=SIN(phi)*SIN(phi0Adv)
      cos_tmp=COS(phi)*COS(phi0Adv)
      sin_tmp2=SIN(phi)*SIN(phi1Adv)
      cos_tmp2=COS(phi)*COS(phi1Adv)
! great circle distance without 'a'
      r= ACOS(sin_tmp + cos_tmp*COS(lam-lam0Adv)) 
      r2=ACOS(sin_tmp2 + cos_tmp2*COS(lam-lam1Adv)) 
      d1=MIN(1.0d0, (r/RR)**2 + ((z-z0)/ZZ)**2 )
      d2=MIN(1.0d0, (r2/RR)**2 + ((z-z0)/ZZ)**2 )
      DummyStart1=0.5d0*(1.0d0+COS(pi*d1))+0.5d0*(1.0d0+COS(pi*d2))
    CASE ('DCMIP_1_2')
      IF (z<z2 .AND.z>z1) then
        DummyStart1=0.5d0*(1.0d0+COS(2.0d0*Pi*(z-z0)/(z2-z1)))
      ELSE
        DummyStart1=0.0d0
      END IF
    CASE ('DCMIP_1_3')
      r=ACOS(SIN(phiP)*SIN(phi)+COS(phiP)*COS(phi)*COS(lam-lamP))
      d1=ABS(z-zp1)
      DummyStart1=0.0d0
      IF (d1<0.5d0*dzp1.AND.r<Rp) THEN
        DummyStart1=0.25d0*(1.0d0+COS(2.0d0*Pi*d1/dzp1)) &
                          *(1.0d0+COS(Pi*r/Rp)) 
      END IF  
    CASE ('DCMIP_4_1')
      U0=VelMax
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      eta=eta_n1(rot_lam,rot_phi,z,U0)
      DummyStart1=ABS(EPV(rot_lam,rot_phi,eta))*ScalingFactor ! the value of |EPV| scales with X
    CASE ('BaroIn')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      perturb_lam = Pi/9.d0
      perturb_phi = 11.d0*Pi/36.d0
      eta_c=0.6d0
      U0=VelMax
      eta=eta_n1(rot_lam,rot_phi,z,U0)
      sin_tmp = SIN(perturb_phi)*SIN(rot_phi)
      cos_tmp = COS(perturb_phi)*COS(rot_phi)
      circdis = ACOS( sin_tmp + cos_tmp*COS(rot_lam-perturb_lam) )    ! great circle distance
      tmp = EXP(- ((circdis*radius)**2.0d0 + ((eta-eta_c)/0.1d0)**2.0d0))
      IF (ABS(tmp)<1.0d-8) tmp = 0.0d0
      DummyStart1 = tmp
    CASE DEFAULT
      DummyStart1=0.d0
  END SELECT   
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: U0,eta_c,eta
  REAL(RealKind) :: sin_tmp,cos_tmp,circdis
  REAL(RealKind) :: sin_tmp2,cos_tmp2
  REAL(RealKind) :: r2,d1,d2
  REAL(RealKind) :: perturb_lam,perturb_phi,tmp
  REAL(RealKind) :: T
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1')
      sin_tmp=SIN(phi)*SIN(phi0Adv)
      cos_tmp=COS(phi)*COS(phi0Adv)
      sin_tmp2=SIN(phi)*SIN(phi1Adv)
      cos_tmp2=COS(phi)*COS(phi1Adv)
! great circle distance without 'a'
      r= ACOS(sin_tmp + cos_tmp*COS(lam-lam0Adv)) 
      r2=ACOS(sin_tmp2 + cos_tmp2*COS(lam-lam1Adv)) 
      d1=MIN(1.0d0, (r/RR)**2 + ((z-z0)/ZZ)**2 )
      d2=MIN(1.0d0, (r2/RR)**2 + ((z-z0)/ZZ)**2 )
      DummyStart2=0.9d0-0.8d0*(0.5d0*(1.0d0+COS(pi*d1))+0.5d0*(1.0d0+COS(pi*d2)))**2
    CASE ('DCMIP_1_2')
      DummyStart2=1.0d0
    CASE ('DCMIP_1_3')
      r=ACOS(SIN(phiP)*SIN(phi)+COS(phiP)*COS(phi)*COS(lam-lamP))
      d1=ABS(z-zp2)
      DummyStart2=0.0d0
      IF (d1<0.5d0*dzp2.AND.r<Rp) THEN
        DummyStart2=0.25d0*(1.0d0+COS(2.0d0*Pi*d1/dzp2)) &
                          *(1.0d0+COS(Pi*r/Rp)) 
      END IF  
    CASE ('BaroIn')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      perturb_lam = Pi/9.d0
      perturb_phi = 11.d0*Pi/36.d0
      eta_c=1.d0
      U0=VelMax
      eta=eta_n1(rot_lam,rot_phi,z,U0)
      sin_tmp = SIN(perturb_phi)*SIN(rot_phi)
      cos_tmp = COS(perturb_phi)*COS(rot_phi)
      circdis = ACOS( sin_tmp + cos_tmp*COS(rot_lam-perturb_lam) )    ! great circle distance
      tmp = EXP(- ((circdis*radius)**2.0d0 + ((eta-eta_c)/0.1d0)**2.0d0))
      IF (ABS(tmp)<1.0d-8) tmp = 0.0d0
      DummyStart2 = tmp
    CASE ('DCMIP_4_1')
      U0=VelMax
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      eta=eta_n1(rot_lam,rot_phi,z,ABS(U0))
      T=temperature(rot_lam,rot_phi,eta,ABS(U0)) 
      DummyStart2=T/eta**kappa
    CASE DEFAULT
      DummyStart2 = 0.d0
  END SELECT   
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: rot_lam,rot_phi
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: sin_tmp,cos_tmp
  REAL(RealKind) :: sin_tmp2,cos_tmp2
  REAL(RealKind) :: r2,d1,d2
  REAL(RealKind) :: eta
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1')
      sin_tmp=SIN(phi)*SIN(phi0Adv)
      cos_tmp=COS(phi)*COS(phi0Adv)
      sin_tmp2=SIN(phi)*SIN(phi1Adv)
      cos_tmp2=COS(phi)*COS(phi1Adv)
! great circle distance without 'a'
      r= ACOS(sin_tmp + cos_tmp*COS(lam-lam0Adv)) 
      r2=ACOS(sin_tmp2 + cos_tmp2*COS(lam-lam1Adv)) 
      d1=MIN(1.0d0, (r/RR)**2 + ((z-z0)/ZZ)**2 )
      d2=MIN(1.0d0, (r2/RR)**2 + ((z-z0)/ZZ)**2 )
      IF (d1<RR) THEN
        DummyStart3=1.0d0
      ELSE IF (d2<RR) THEN
        DummyStart3=1.0d0
      ELSE  
        DummyStart3=0.1d0
      END IF  
      IF (ABS(phi)<0.125d0) THEN
        DummyStart3=0.1d0 
      END IF  
    CASE ('DCMIP_1_3')
      r=ACOS(SIN(phiP)*SIN(phi)+COS(phiP)*COS(phi)*COS(lam-lamP))
      d1=ABS(z-zp3)
      DummyStart3=0.0d0
      IF (d1<0.5d0*dzp3.AND.r<Rp) THEN
        DummyStart3=1.0d0
      END IF  
    CASE ('BaroIn')
      CALL Rotate(lam,phi,rot_lam,rot_phi)
      DummyStart3 = 0.5d0*( TANH( 3.0d0*ABS(rot_phi)-pi ) + 1.0d0)
    CASE DEFAULT
      DummyStart3 = 0.0d0
  END SELECT   
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: sin_tmp,cos_tmp
  REAL(RealKind) :: sin_tmp2,cos_tmp2
  REAL(RealKind) :: r2,d1,d2
  REAL(RealKind) :: Dummy1,Dummy2,Dummy3
  SELECT CASE(Problem)
    CASE ('DCMIP_1_1')
      sin_tmp=SIN(phi)*SIN(phi0Adv)
      cos_tmp=COS(phi)*COS(phi0Adv)
      sin_tmp2=SIN(phi)*SIN(phi1Adv)
      cos_tmp2=COS(phi)*COS(phi1Adv)
! great circle distance without 'a'
      r= ACOS(sin_tmp + cos_tmp*COS(lam-lam0Adv)) 
      r2=ACOS(sin_tmp2 + cos_tmp2*COS(lam-lam1Adv)) 
      d1=MIN(1.0d0, (r/RR)**2 + ((z-z0)/ZZ)**2 )
      d2=MIN(1.0d0, (r2/RR)**2 + ((z-z0)/ZZ)**2 )
      Dummy1=0.5d0*(1.0d0+COS(pi*d1))+0.5d0*(1.0d0+COS(pi*d2))
      Dummy2=0.9d0-0.8d0*Dummy1**2
      IF (d1<RR) THEN
        Dummy3=1.0d0
      ELSE IF (d2<RR) THEN
        Dummy3=1.0d0
      ELSE  
        Dummy3=0.1d0
      END IF  
      IF (ABS(phi)<0.125d0) THEN
        Dummy3=0.1d0 
      END IF  
      DummyStart4=1.0d0-0.3d0*(Dummy1+Dummy2+Dummy3)
    CASE ('DCMIP_1_3')
      r=ACOS(SIN(phiP)*SIN(phi)+COS(phiP)*COS(phi)*COS(lam-lamP))
      DummyStart4=0.0d0
      d1=ABS(z-zp1)
      IF (d1<0.5d0*dzp1.AND.r<Rp) THEN
        DummyStart4=DummyStart4+0.25d0*(1.0d0+COS(2.0d0*Pi*d1/dzp1)) &
                                      *(1.0d0+COS(Pi*r/Rp)) 
      END IF  
      d1=ABS(z-zp2)
      IF (d1<0.5d0*dzp2.AND.r<Rp) THEN
        DummyStart4=DummyStart4+0.25d0*(1.0d0+COS(2.0d0*Pi*d1/dzp2)) &
                                      *(1.0d0+COS(Pi*r/Rp)) 
      END IF  
      d1=ABS(z-zp3)
      IF (d1<0.5d0*dzp3.AND.r<Rp) THEN
        DummyStart4=DummyStart4+1.0d0
      END IF  
    CASE ('BaroIn')
      DummyStart4=1.0d0
    CASE DEFAULT
      DummyStart4=0.0d0
  END SELECT   
END FUNCTION DummyStart4

FUNCTION QiStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION RhoStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  RhoStart=One
END FUNCTION RhoStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =Zero
END FUNCTION TStart


FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=Zero
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=Zero
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=Zero
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=Zero
END FUNCTION ForceRho


FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=Zero
END FUNCTION HeightFun

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE DCMIP_Mod
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  ThStartSoil=T0Vo
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  QvStartSoil=0.0d0
END FUNCTION QvStartSoil


FUNCTION QsStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QsStart = 0.0
END FUNCTION QsStart

FUNCTION NvStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NvStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NvStart = 0.0
END FUNCTION NvStart

FUNCTION NcStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NcStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NcStart = 0.0
END FUNCTION NcStart

FUNCTION NrStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NrStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NrStart = 0.0
END FUNCTION NrStart

FUNCTION NiStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NiStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NiStart = 0.0
END FUNCTION NiStart

FUNCTION NsStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NsStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NsStart = 0.0
END FUNCTION NsStart

FUNCTION OmeStart(lam,phi,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  OmeStart = 0.0
END FUNCTION OmeStart


SUBROUTINE Rotate(lam,phi,rot_lam,rot_phi)
  USE DCMIP_Mod
  IMPLICIT NONE!
  REAL(RealKind) :: lam,phi,rot_lam,rot_phi
  IF (ABS(rotation_angle)<1.0E-8) THEN
    rot_lam = lam
    rot_phi = phi
  ELSE
     CALL regrot(lam,phi,rot_lam,rot_phi,0.0d0,-0.5d0*pi+rotation_angle,1)
  ENDIF
END SUBROUTINE Rotate

SUBROUTINE regrot(pxreg,pyreg,pxrot,pyrot,pxcen,pycen,kcall)
  USE DCMIP_Mod
  IMPLICIT NONE!
!----------------------------------------------------------------------
!
!*    conversion between regular and rotated spherical coordinates.
!*
!*    pxreg     longitudes of the regular coordinates
!*    pyreg     latitudes of the regular coordinates
!*    pxrot     longitudes of the rotated coordinates
!*    pyrot     latitudes of the rotated coordinates
!*              all coordinates given in degrees n (negative for s)
!*              and degrees e (negative values for w)
!*    pxcen     regular longitude of the south pole of the rotated grid
!*    pycen     regular latitude of the south pole of the rotated grid
!*
!*    kcall=-1: find regular as functions of rotated coordinates.
!*    kcall= 1: find rotated as functions of regular coordinates.
!
!-----------------------------------------------------------------------
!
  INTEGER :: kxdim,kydim,kx,ky,kcall
  REAL(RealKind) :: pxreg,pyreg,&
                    pxrot,pyrot,&
                    pxcen,pycen
  REAL(RealKind) :: zsycen,zcycen,zxmxc,zsxmxc,zcxmxc,zsyreg,zcyreg, &
                    zsyrot,zcyrot,zcxrot,zsxrot,zpi,zpih
  INTEGER  :: jy,jx

  zpih = Pi*0.5d0
  zsycen = SIN((pycen+zpih))
  zcycen = COS((pycen+zpih))
!
  IF (kcall.EQ.1) THEN
     zxmxc  = pxreg - pxcen
     zsxmxc = SIN(zxmxc)
     zcxmxc = COS(zxmxc)
     zsyreg = SIN(pyreg)
     zcyreg = COS(pyreg)
     zsyrot = zcycen*zsyreg - zsycen*zcyreg*zcxmxc
     zsyrot = max(zsyrot,-1.d0)
     zsyrot = min(zsyrot,+1.d0)
     
     pyrot = ASIN(zsyrot)
    
     zcyrot = COS(pyrot)
     zcxrot = (zcycen*zcyreg*zcxmxc +zsycen*zsyreg)/zcyrot
     zcxrot = max(zcxrot,-1.d0)
     zcxrot = min(zcxrot,+1.d0)
     zsxrot = zcyreg*zsxmxc/zcyrot
   
     pxrot = ACOS(zcxrot)
  
     IF (zsxrot<0.0) pxrot = -pxrot
  ELSE IF (kcall.EQ.-1) THEN
     zsxrot = SIN(pxrot)
     zcxrot = COS(pxrot)
     zsyrot = SIN(pyrot)
     zcyrot = COS(pyrot)
     zsyreg = zcycen*zsyrot + zsycen*zcyrot*zcxrot
     zsyreg = max(zsyreg,-1.d0)
     zsyreg = min(zsyreg,+1.d0)
       
     pyreg = ASIN(zsyreg)
         
     zcyreg = COS(pyreg)
     zcxmxc = (zcycen*zcyrot*zcxrot -&
              zsycen*zsyrot)/zcyreg
     zcxmxc = max(zcxmxc,-1.d0)
     zcxmxc = min(zcxmxc,+1.d0)
     zsxmxc = zcyrot*zsxrot/zcyreg
     zxmxc  = ACOS(zcxmxc)
     IF (zsxmxc<0.0) zxmxc = -zxmxc
         
     pxreg = zxmxc + pxcen
          
  ELSE
     WRITE(6,'(1x,''invalid kcall in regrot'')')
     STOP
  ENDIF
END SUBROUTINE regrot

SUBROUTINE turnwi(puarg,pvarg,pures,pvres,         &
                  pxreg,pyreg,pxrot,pyrot,   &
                  pxcen,pycen,kcall)
!
  USE DCMIP_Mod
  IMPLICIT NONE
!
!-----------------------------------------------------------------------
!
!*    turn horizontal velocity components between regular and
!*    rotated spherical coordinates.
!
!*    puarg : input u components
!*    pvarg : input v components
!*    pures : output u components
!*    pvres : output v components
!*    p_a    : transformation coefficients
!*    p_b    :    -"-
!*    p_c    :    -"-
!*    p_d    :    -"-
!*    pxreg : regular longitudes
!*    pyreg : regular latitudes
!*    pxrot : rotated longitudes
!*    pyrot : rotated latitudes
!*    kxdim              : dimension in the x (longitude) direction
!*    kydim              : dimension in the y (latitude) direction
!*    kx                 : number of gridpoints in the x direction
!*    ky                 : number of gridpoints in the y direction
!*    pxcen              : regular longitude of the south pole of the
!*                         transformed grid
!*    pycen              : regular latitude of the south pole of the
!*                         transformed grid
!*
!*    kcall < 0          : find wind components in regular coordinates
!*                         from wind components in rotated coordinates
!*    kcall > 0          : find wind components in rotated coordinates
!*                         from wind components in regular coordinates
!*    note that all coordinates are given in degrees n and degrees e.
!*       (negative values for s and w)
!
!-----------------------------------------------------------------------

  INTEGER kxdim,kydim,kx,ky,kcall
  REAL(RealKind) puarg,pvarg,    &
                 pures,pvres,    &
                 p_a,   p_b,       &
                 p_c,   p_d,       &
                 pxreg,pyreg,    &
                 pxrot,pyrot
  REAL(RealKind) pxcen,pycen
!-----------------------------------------------------------------------
  INTEGER jy,jx
  REAL(RealKind) zpih,zsyc,zcyc,zsxreg,zcxreg,zsyreg,zcyreg,zxmxc,&
                 zsxmxc,zcxmxc,zsxrot,zcxrot,zsyrot,zcyrot
!-----------------------------------------------------------------------
  IF (kcall.EQ.1) THEN
     zpih = Pi*0.5d0
     zsyc = SIN(pycen+zpih)
     zcyc = COS(pycen+zpih)
     !
     zsxreg = SIN(pxreg)
     zcxreg = COS(pxreg)
     zsyreg = SIN(pyreg)
     zcyreg = COS(pyreg)
     !
     zxmxc  = pxreg - pxcen
     zsxmxc = SIN(zxmxc)
     zcxmxc = COS(zxmxc)
     !
     zsxrot = SIN(pxrot)
     zcxrot = COS(pxrot)
     zsyrot = SIN(pyrot)
     zcyrot = COS(pyrot)
     !
     p_a = zcyc*zsxmxc*zsxrot + zcxmxc*zcxrot
     p_b = zcyc*zcxmxc*zsyreg*zsxrot - zsyc*zcyreg*zsxrot - &
          zsxmxc*zsyreg*zcxrot
     p_c = zsyc*zsxmxc/zcyrot
     p_d = (zsyc*zcxmxc*zsyreg + zcyc*zcyreg)/zcyrot
     !
     pures = p_a*puarg + p_b*pvarg
     pvres = p_c*puarg + p_d*pvarg
  ELSE IF (kcall.EQ.-1) THEN
     zpih = Pi*0.5d0
     zsyc = SIN(pycen+zpih)
     zcyc = COS(pycen+zpih)
     !
     zsxreg = SIN(pxreg)
     zcxreg = COS(pxreg)
     zsyreg = SIN(pyreg)
     zcyreg = COS(pyreg)
     !
     zxmxc  = pxreg - pxcen
     zsxmxc = SIN(zxmxc)
     zcxmxc = COS(zxmxc)
     !
     zsxrot = SIN(pxrot)
     zcxrot = COS(pxrot)
     zsyrot = SIN(pyrot)
     zcyrot = COS(pyrot)
     !
     p_a = zcxmxc*zcxrot + zcyc*zsxmxc*zsxrot
     p_b = zcyc*zsxmxc*zcxrot*zsyrot + zsyc*zsxmxc*zcyrot -&
          zcxmxc*zsxrot*zsyrot
     p_c =-zsyc*zsxrot/zcyreg
     p_d = (zcyc*zcyrot - zsyc*zcxrot*zsyrot)/zcyreg
     !
     pures = p_a*puarg + p_b*pvarg
     pvres = p_c*puarg + p_d*pvarg
  ELSE
     WRITE(6,'(1x,''invalid kcall in turnwi'')')
     STOP
  ENDIF
END SUBROUTINE turnwi


FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start
FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE DCMIP_Mod
  USE Kind_Mod
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=Zero
END FUNCTION ForceTh
