MODULE GlobalAcoustic_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: VelMax=0.0d0
  REAL(RealKind), PARAMETER :: t0=300.0d0
  REAL(RealKind), PARAMETER :: D0=1.d-1
  REAL(RealKind), PARAMETER :: DeltaP=100.0d0 !Pa
  REAL(RealKind) :: zT
  REAL(RealKind) :: R
  REAL(RealKind) :: phi0=0.0d0
  REAL(RealKind) :: lam0=0.0d0


  NAMELIST /Example/ phi0 &
                    ,lam0   

END MODULE GlobalAcoustic_Mod

SUBROUTINE InputExample(FileName)
  USE GlobalAcoustic_Mod
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
  zT=Domain%z1
  R=RadEarth/3.0d0
END SUBROUTINE InputExample

FUNCTION RhoProf(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: pLoc

  ploc=p0*EXP(-z/(Rd*t0))
  RhoProf=pLoc/(Rd*t0)
END FUNCTION RhoProf

FUNCTION ThProfFun(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time

  ThProfFun=t0*EXP(kappa*z/(Rd*t0))

END FUNCTION ThProfFun

FUNCTION PreProf(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: lam,phi,z,zHeight,Time

  PreProf=p0*EXP(-z/(Rd*t0))

END FUNCTION PreProf

FUNCTION ThStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  USE ThProf_Mod
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: pLoc

  ThStart=t0*EXP(kappa*z/(Rd*t0))
  
END FUNCTION ThStart

FUNCTION RhoFun(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: pLoc,rad,f,g

  ploc=p0*EXP(-z/(Rd*t0))
  rad=RadEarth*ACOS(SIN(phi0)*SIN(phi)+COS(phi0)*COS(phi)*COS(lam-lam0))
  f=Zero
  IF (rad<R) THEN
    f=Half*(One+COS(Pi*rad/R))
  END IF
  g=SIN(Pi*z/zT)
  pLoc=ploc+DeltaP*f*g
  RhoFun=pLoc/(Rd*t0)
END FUNCTION RhoFun

FUNCTION PreStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time

  PreStart=p0*EXP(-z/(Rd*t0))

END FUNCTION PreStart

FUNCTION QvProfFun(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  USE ThProf_Mod
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QvProfFun=0.0d0
END FUNCTION QvProfFun

FUNCTION UStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: UStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: phiL,phiR
  REAL(RealKind) :: U0
  U0=VelMax
  UStart=U0*COS(phi)
END FUNCTION UStart

FUNCTION UStartE(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  USE UVW_Mod
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  UStartE=UStart(lam,phi,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: VStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  VStart=0.0d0
END FUNCTION VStart

FUNCTION VStartE(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  VStartE=VStart(lam,phi,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: WStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  TkeStart=0.0d0
END FUNCTION TkeStart

FUNCTION DisStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DisStart=0.0d0
END FUNCTION DisStart

FUNCTION QvStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QvStart=0.d0
END FUNCTION QvStart

FUNCTION QcStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: DStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DummyStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION QiStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION RhoStart(lam,phi,z,zHeight,Time)
  USE GlobalAcoustic_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  RhoStart=RhoFun(lam,phi,z,zHeight,Time)  &
          -RhoProf(lam,phi,z,zHeight,Time)
END FUNCTION RhoStart


