MODULE Profile_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: th0=293.16d0
  REAL(RealKind) :: D0=1.0d0
  REAL(RealKind) :: uMax=10.0d0
  REAL(RealKind) :: vMax=0.0d0
  REAL(RealKind) :: TkeMax=0.0d0
  REAL(RealKind) :: DisMax=0.0d0
  REAL(RealKind) :: qvr=0.95d0

  NAMELIST /Example/        &
                    uMax ,  &
                    vMax ,  &
                    th0  ,  &
                    TkeMax, &
                    DisMax, &
                    D0,     &
                    qvr,    &
                    N

END MODULE Profile_Mod

SUBROUTINE InputExample(FileName)
  USE Profile_Mod
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
END SUBROUTINE InputExample

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE Profile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,pLoc
  S=N*N/Grav
  ThLoc=th0*exp(z*S)
  IF (N>Zero) THEN
    pLoc=p0*(One-Grav/(cp*th0*S)*(One-EXP(-S*z)))**(cP/Rd)
  ELSE
    pLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(cP/Rd)
  END IF
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE Profile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,pLoc
  S=N*N/Grav
  ThLoc=th0*exp(z*S)
  IF (N>Zero) THEN
    pLoc=p0*(One-Grav/(cp*th0*S)*(One-EXP(-S*z)))**(cP/Rd)
  ELSE
    pLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(cP/Rd)
  END IF
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoProf


FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S
  S=N*N/Grav
  ThProfFun=th0*exp(z*S)*RhoFun(x,y,z,zHeight,Time)
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod 
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: RhoLoc,ThLoc
  RhoLoc=RhoFun(x,y,z,zHeight,Time)
  ThLoc=thProfFun(x,y,z,zHeight,Time)
  QvProfFun=qvr*qvs(RhoLoc,ThLoc,ThLoc)
END FUNCTION QvProfFun

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: U0
  UStart=uMax*RhoFun(x,y,z,zHeight,Time)
END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=vMax*RhoFun(x,y,z,zHeight,Time)
END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0*RhoFun(x,y,z,zHeight,Time)
END FUNCTION WStart

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThStart=ThProfFun(x,y,z,zHeight,Time)
END FUNCTION ThStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=Zero
END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax*RhoFun(x,y,z,zHeight,Time)

END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax*RhoFun(x,y,z,zHeight,Time)
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=D0*RhoFun(x,y,z,zHeight,Time)
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,RadLoc,ThBack
  S=N*N/Grav
  IF (N>Zero) THEN
    PreStart=p0*(One-Grav/(cp*th0*S)*(One-EXP(-S*z)))**(cP/Rd)
  ELSE
    PreStart=p0*(One-kappa*Grav*z/(Rd*th0))**(cP/Rd)
  END IF
END FUNCTION PreStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Profile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart



