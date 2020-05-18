MODULE Thermodynamic_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Control_Mod
  USE DataType_Mod
  USE Transport_Mod
  USE Physics_Mod

  IMPLICIT NONE

  INTERFACE Pressure
    MODULE PROCEDURE PressureT,PressureTheta
  END INTERFACE Pressure


CONTAINS

FUNCTION gammln(xx)
  REAL(RealKind) :: gammln,xx
  INTEGER :: j
  REAL(RealKind) :: ser,tmp,x,y
  REAL(RealKind), PARAMETER ::  &
     cof(1:6)=(/76.18009172947146d0,-86.50532032941677d0 &
               ,24.01409824083091d0,-1.231739572450155d0 &
               ,.1208650973866179d-2,-.5395239384953d-5/) &
    ,stp=2.5066282746310005d0
! REAL(RealKind) :: cof,stp/76.18009172947146d0,-86.50532032941677d0,
!*24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,
!*-.5395239384953d-5,2.5066282746310005d0/
  x=xx
  y=x
  tmp=x+5.5d0
  tmp=(x+0.5d0)*log(tmp)-tmp
  ser=1.000000000190015d0
  DO j=1,6
    y=y+1.d0
    ser=ser+cof(j)/y
  END DO
  gammln=tmp+log(stp*ser/x)
END FUNCTION gammln

FUNCTION gfct(x)
  REAL(RealKind) :: gfct
  REAL(RealKind) :: x,xx,tmp,ser,gamma
  INTEGER :: j
  REAL(RealKind), PARAMETER:: &
       cof(1:6)=(/76.18009173d0,-86.50532033d0,24.01409822d0,  &
          &     -1.231739516d0,.120858003d-2,-.536382d-5/)&
           ,stp=2.50662827465d0
     ! DATA half,one,fpf/0.5d0,1.0d0,5.5d0/

    xx  = x  - 1.0d0
    tmp = xx + 5.5d0
    tmp = (xx + 0.5d0) * LOG(tmp) - tmp
    ser = 1.0d0
    DO j = 1,6
      xx  = xx  + 1.0d0
      ser = ser + cof(j) / xx
    END DO
    gamma = tmp + LOG(stp*ser)
    gamma = EXP(gamma)

    gfct = gamma
END FUNCTION gfct

FUNCTION FallF(RhoL)

  REAL(RealKind) :: FallF 
  REAL(RealKind) :: RhoL

  REAL(RealKind) :: Lambda

  REAL(RealKind), PARAMETER :: a_qc=842.0d0     ! m^(1-b_qc) s^-1
  REAL(RealKind), PARAMETER :: b_qc=0.8d0       !
  REAL(RealKind), PARAMETER :: rho_qc=1000.0d0  ! kg m^-3
  REAL(RealKind), PARAMETER :: N_qc=8.0d6       ! m^-4
  REAL(RealKind), PARAMETER :: c_r=130.0d0      ! m^(1/2) s^-1
  
  IF (RhoL>1d-10.AND.PrecipRain) THEN
!   Lambda=(Pi*Rho_qc*N_qc/(RhoL+Eps))**(1.0d0/3.0d0)
!   FallF=-a_qc*exp(gammln(4.0d0+b_qc))/(6.0d0*lambda**b_qc) ! Reisner
    FallF=-c_r*EXP(gammln(4.5e0_RealKind))/6.0d0*(RhoL/(Pi*Rho_qc*N_qc))**(1.0d0/8.0d0) ! Wacker
  ELSE
    FallF=Zero
  END IF
END FUNCTION FallF

FUNCTION LatHeat(T)

  REAL(RealKind) ::  LatHeat
  REAL(RealKind) ::  T 

! T Absolute temperature
! L00=3.14762e6 J kg^-1 Latent heat at 0 K 
! Cpv=1850.0 J kg^-1 K^-1 specific heat of vapor at constant pressure
! Cpl=4218.0 J kg^-1 K^-1 specific heat of liquid water

  LatHeat=L00+(Cpv-Cpl)*T

END FUNCTION LatHeat

FUNCTION LatHeatFreeze(T)

  REAL(RealKind) :: LatHeatFreeze
  REAL(RealKind) :: T
! T Absolute Temperatur
! LatwaFREEZE= 0.334d6 latent heat at T=273.16K
! CIce = 2110 J kg^-1 K^-1 latent heat of ice
! CWater = 4186 J kg^-1 K^-1 latent heat of liquid water

  LatHeatFreeze=LatwaFREEZE+(Cpl-Cpi)*T

END FUNCTION LatHeatFreeze

FUNCTION SaturVapor(T)

  REAL(RealKind) ::  SaturVapor
  REAL(RealKind) ::  T
  REAL(RealKind) ::  T_C

! Guide to Meteorological Instruments and Methods of Observation (CIMO Guide) (WMO, 2008)
  T_C=T-273.15d0
  SaturVapor=6.112d2*EXP(17.62d0*T_C/(243.12d0+T_C))  

END FUNCTION SaturVapor

FUNCTION SaturVaporIce(T)

  REAL(RealKind) ::  SaturVaporIce
  REAL(RealKind) ::  T
  REAL(RealKind) ::  T_C

! Guide to Meteorological Instruments and Methods of Observation (CIMO Guide) (WMO, 2008)
  T_C=T-273.15d0
  SaturVaporIce=6.112d2*EXP(22.46d0*T_C/(272.62d0+T_C))

END FUNCTION SaturVaporIce


FUNCTION AbsTemp(RhoD,RhoV,p)

  REAL(RealKind) ::  AbsTemp
  REAL(RealKind) ::  RhoD,RhoV,p

! Rho density
! RhoV density of water vapor
! T Absolute temperature
! Rd=287.0 J kg^-1 K^-1 gas constant of dry air
! Rv=461.50 J kg^-1 K^-1 gas constant of vapor

  AbsTemp=p/(Rd*RhoD+Rv*RhoV+Eps)

END FUNCTION AbsTemp

FUNCTION CpmlFun(RhoD,RhoV,RhoL,RhoI)

  REAL(RealKind) :: CpmlFun
  REAL(RealKind) :: RhoD,RhoV,RhoL,RhoI

  CpmlFun=RhoD*Cpd+RhoV*Cpv+RhoL*Cpl+RhoI*Cpi+Eps

END FUNCTION CpmlFun

FUNCTION CvmlFun(RhoD,RhoV,RhoL,RhoI)

  REAL(RealKind) :: CvmlFun
  REAL(RealKind) :: RhoD,RhoV,RhoL,RhoI

  CvmlFun=RhoD*Cvd+RhoV*Cvv+RhoL*Cpl+RhoI*Cpi+Eps

END FUNCTION CvmlFun

FUNCTION RmFun(RhoD,RhoV)

  REAL(RealKind) :: RmFun
  REAL(RealKind) :: RhoD,RhoV

  RmFun=RhoD*Rd+RhoV*Rv+Eps

END FUNCTION RmFun

FUNCTION PressureT(RhoD,RhoV,T)

  REAL(RealKind) ::  PressureT
  REAL(RealKind) ::  RhoD,RhoV,T

! Rho density 
! RhoV density of water vapor 
! T Absolute temperature
! Rd=287.0 J kg^-1 K^-1 gas constant of dry air
! Rv=461.50 J kg^-1 K^-1 gas constant of vapor

  PressureT=(Rd*RhoD+Rv*RhoV)*T

END FUNCTION PressureT

FUNCTION PressureV(RhoV,T)

  REAL(RealKind) ::  PressureV
  REAL(RealKind) ::  RhoV,T

! Rho density
! RhoV density of water vapor
! T Absolute temperature
! Rd=287.0 J kg^-1 K^-1 gas constant of dry air
! Rv=461.50 J kg^-1 K^-1 gas constant of vapor

  PressureV=Rv*RhoV*T

END FUNCTION PressureV


FUNCTION PressureTheta(RhoD,RhoV,RhoL,RhoI,RhoPot)

  REAL(RealKind) ::  PressureTheta
  REAL(RealKind) ::  RhoD,RhoV,RhoL,RhoI,RhoPot

  REAL(RealKind) :: Rm,Cpml,KappaLoc

! RhoPotT generalized potential temperature times density
! Rd=287.0 J kg^-1 K^-1 gas constant of dry air
! Rv=461.50 J kg^-1 K^-1 gas constant of vapor

! Rm=Rd*RhoD+Rv*RhoV
! Cpml=Cpd*RhoD+Cpv*RhoV+Cpl*RhoL+Eps
  Rm=RmFun(RhoD,RhoV)
  Cpml=CpmlFun(RhoD,RhoV,RhoL,RhoI)
  KappaLoc=Rm/(Cpml+Eps)
  PressureTheta=p0*(Rd*RhoPot/p0)**(One/(One-KappaLoc))

END FUNCTION PressureTheta


FUNCTION qvs(RhoDry,PotM,RhoVap)

  REAL(RealKind) :: qvs
  REAL(RealKind) :: RhoDry,PotM,RhoVap

  REAL(RealKind) :: p,tAbs,e

  p=p0*(Rd*PotM/p0)**(gamma)
  tAbs=p/(Rd*RhoDry+Rv*RhoVap)
  IF (tAbs>=273.15d0) THEN
    e=611.2d0*EXP(17.67d0*(tAbs-273.15d0)/(tAbs-29.65d0))
  ELSE
    e=611.2d0*EXP(22.514d0-6150.d0/tAbs)
  END IF
  qvs=Rd/Rv*(e/(p-e))

END FUNCTION qvs

FUNCTION qvsP(RhoDry,p,qv)

  REAL(RealKind) :: qvsP
  REAL(RealKind) :: RhoDry,p,qv

  REAL(RealKind) :: tAbs,e

  tAbs=p/(RhoDry*(Rd+Rv*qv))
  IF (tAbs>=273.15d0) THEN
    e=611.2d0*EXP(17.67d0*(tAbs-273.15d0)/(tAbs-29.65d0))
  ELSE
!   e=611.2d0*EXP(22.514d0-6150.d0/tAbs)
    e=611.2d0*EXP(22.515101592d0-6150.d0/tAbs)
  END IF
  qvsP=Rd/Rv*(e/(p-e))

END FUNCTION qvsP



SUBROUTINE qvsJac(RhoDry,PotM,RhoVap,qvs,dqvsdRhoDry,dqvsdPotM,dqvsdRhoVap  &
		   ,dpdPotM,dqvsdp)

  REAL(RealKind) :: RhoDry,PotM,RhoVap
  REAL(RealKind) :: dpdPotM,dpdRhoDry,dpdRhoVap
  REAL(RealKind) :: dtAbsdp,dtAbsdPotM,dtAbsdRhoDry,dtAbsdRhoVap
  REAL(RealKind) :: dedtAbs,dedPotM,dedRhoDry,dedRhoVap
  REAL(RealKind) :: dqvsdp,dqvsde,dqvsdRhoDry,dqvsdPotM,dqvsdRhoVap

  REAL(RealKind) :: p,tAbs,e,qvs

  p=p0*(Rd*PotM/p0)**(gamma)
  dpdPotM=p0*(Rd/p0)**(gamma)*gamma*PotM**(gamma-1)
  dpdRhoDry=Zero
  dpdRhoVap=Zero
  tAbs=p/(Rd*RhoDry+Rv*RhoVap)
  dtAbsdp=1.0d0/(Rd*RhoDry+Rv*RhoVap)
  dtAbsdPotM=dtAbsdp*dpdPotM
  dtAbsdRhoDry=-p*Rd/(Rd*RhoDry+Rv*RhoVap)**2
  dtAbsdRhoVap=-p*Rv/(Rd*RhoDry+Rv*RhoVap)**2
  IF (tAbs>=273.15d0) THEN
    e=611.2d0*EXP(17.67d0*(tAbs-273.15d0)/(tAbs-29.65d0))
    dedtAbs=(17.67d0/(tAbs-29.65d0) &
           *EXP(17.67d0*(tAbs-273.15d0)/(tAbs-29.65d0))) &
           -(17.67d0*(tAbs-273.15d0)/(tAbs-29.65d0)**2  &
           *EXP(17.67d0*(tAbs-273.15d0)/(tAbs-29.65d0)))
    dedPotM=dedtAbs*dtAbsdPotM
    dedRhoDry=dedtAbs*dtAbsdRhoDry
    dedRhoVap=dedtAbs*dtAbsdRhoVap
  ELSE
!   e=611.2d0*EXP(22.514d0-6150.d0/tAbs)
    e=611.2d0*EXP(22.515101592d0-6150.d0/tAbs)
!   dedtAbs=-611.2d0*6150.d0/tAbs**2*EXP(22.514d0-6150.d0/tAbs)
    dedtAbs=-611.2d0*6150.d0/tAbs**2*EXP(22.515101592d0-6150.d0/tAbs)
    dedPotM=dedtAbs*dtAbsdPotM
    dedRhoDry=dedtAbs*dtAbsdRhoDry
    dedRhoVap=dedtAbs*dtAbsdRhoVap
  END IF
  qvs=Rd/Rv*(e/(p-e))
  dqvsdp=-Rd/Rv*e/((p-e)**2)
  dqvsde=Rd/Rv*((1/(p-e))+(e/((p-e)**2)))
  dqvsdPotM=dqvsdp*dpdPotM+dqvsde*dedPotM
  dqvsdRhoDry=dqvsdp*dpdRhoDry+dqvsde*dedRhoDry
  dqvsdRhoVap=dqvsdp*dpdRhoVap+dqvsde*dedRhoVap

END SUBROUTINE qvsJac


FUNCTION RelHumidity(T,RhoV)

  REAL(RealKind) :: RelHumidity
  REAL(RealKind) :: T,RhoV

  RelHumidity=RhoV*Rv*T/SaturVapor(T)

END FUNCTION RelHumidity

FUNCTION ThEquiv(T,RhoD,RhoV,RhoL)

  REAL(RealKind) :: ThEquiv
  REAL(RealKind) :: T,RhoD,RhoV,RhoL

  REAL(RealKind) :: Cp_eff

  Cp_eff=Cpd+(RhoV+RhoL)/RhoD*Cpl
  ThEquiv=T*(p0/(RhoD*Rd*T))**(Rd/Cp_eff) &
           *RelHumidity(T,RhoV)**(-RhoV*Rv/(RhoD*Cp_eff)) &
           *EXP(LatHeat(T)*RhoV/(RhoD*Cp_eff*T))

END FUNCTION ThEquiv

FUNCTION DThEquivDT(T,RhoD,RhoV,RhoL)

  REAL(RealKind) :: DThEquivDT
  REAL(RealKind) :: T,RhoD,RhoV,RhoL

  REAL(RealKind) :: Cpml,Rm,Cp_eff

  Rm=Rd+Rv*RhoV/RhoD
  Cpml=Cpd+Cpv*RhoV/RhoD+Cpl*RhoL/RhoD
  Cp_eff=Cpd+(RhoV+RhoL)/RhoD*Cpl
  DThEquivDT=ThEquiv(T,RhoD,RhoV,RhoL)/(Cp_eff*T)*(Cpml-Rm)

END FUNCTION DThEquivDT

SUBROUTINE AbsTNewton(TAbs,Th,RhoD,RhoV,RhoL)
  REAL(RealKind) :: TAbs,Th,RhoD,RhoV,RhoL
  REAL(RealKind) :: Acc,NewtonError,NewtonFunc
  Acc=1.0d-12
  NewtonError=One
! WRITE(*,*) 'TAbs Newton'
  DO WHILE (NewtonError>Acc)
      NewtonFunc=TAbs-(ThEquiv(TAbs,RhoD,RhoV,RhoL)-Th) &
                      /DThEquivDT(TAbs,RhoD,RhoV,RhoL)
      NewtonError=ABS(NewtonFunc-TAbs)
      TAbs=NewtonFunc
!     WRITE(*,*) 'TAbs',TAbs,ThEquiv(TAbs,RhoD,RhoV,RhoL)-Th
  END DO
END SUBROUTINE AbsTNewton

FUNCTION AbsTEnergy(Energy,Rho,RhoV,RhoL,Kinetic,z)
  REAL(RealKind) :: AbsTEnergy
  REAL(RealKind) :: Energy,Rho,RhoV,RhoL,Kinetic,z

  REAL(RealKind) :: RhoD

   RhoD=Rho-RhoV-RhoL
   AbsTEnergy=(Energy &
               -Rho*Kinetic &
               -Rho*Half*z*Grav &
               -RhoV*L00) &
              /(RhoD*Cvd+RhoV*Cvv+RhoL*Cpl)
END FUNCTION AbsTEnergy

FUNCTION PressureEnergy(Energy,Rho,RhoV,RhoL,Kinetic,z)

  REAL(RealKind) :: PressureEnergy
  REAL(RealKind) :: Energy,Rho,RhoV,RhoL,Kinetic,z

  REAL(RealKind) :: RhoD

   RhoD=Rho-RhoV-RhoL
   PressureEnergy=(RhoD*Rd+RhoV*Rv) &
               *(Energy &
               -Rho*Kinetic &
               -Rho*Half*z*Grav &
               -RhoV*L00) &
              /(RhoD*Cvd+RhoV*Cvv+RhoL*Cpl)

END FUNCTION PressureEnergy

END MODULE Thermodynamic_Mod

