MODULE Klein_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind), PARAMETER :: th0=300.0d0
  REAL(RealKind), PARAMETER :: D0=1.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind) :: uMax=10.0d0,vMax=0.0d0
  REAL(RealKind) :: TkeMax=0.0d0,DisMax=0.0d0
  REAL(RealKind) :: zDiff
  REAL(RealKind) :: Rho0=1.0d0
  REAL(RealKind) :: c
  REAL(RealKind) :: gamma1
  REAL(RealKind) :: Delta=1.0d-4
  REAL(RealKind) :: Lenx

  NAMELIST /Example/    &
                    uMax , &
                    vMax , &
                    zDiff, &
                    Delta



END MODULE Klein_Mod

SUBROUTINE InputExample(FileName)
  USE Klein_Mod
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
  gamma1=Cpd/Cvd
  Rho0=p0/(Rd*Th0)
  c=SQRT(gamma1*p0/rho0)
  Lenx=domain%x1-domain%x0
  GravComp=Zero
  WRITE(*,*) 'Gamma1',Gamma
  WRITE(*,*) 'Rho0',Rho0
  WRITE(*,*) 'c',c
END SUBROUTINE InputExample

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time

  PreProf=p0+Delta*Rho0*c*SIN(Two*Pi*x/Lenx)
END FUNCTION PreProf


FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc
  pLoc=p0+Delta*Rho0*c*SIN(Two*Pi*x/Lenx)
  ThLoc=Th0
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc
  pLoc=p0
  ThLoc=Th0
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S
  ThProfFun=th0
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE Rho_Mod 
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: RhoLoc,ThLoc
  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: U0
  UStart=Delta*SIN(Two*Pi*x/Lenx)
END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=vMax
END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThStart=ThProfFun(x,y,z,zHeight,Time)
END FUNCTION ThStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE Rho_Mod
  USE RhoProf_Mod

  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)  &
          -RhoProf(x,y,z,zHeight,Time)

END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  IF (z>=zDiff) THEN
    DStart=D0
  ELSE
    DStart=Zero
  END IF
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0+Delta*Rho0*c*SIN(Two*Pi*x/Lenx)
END FUNCTION PreStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Klein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart



