MODULE Turm_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind), PARAMETER :: VelMax=1.0d0
  REAL(RealKind), PARAMETER :: N=1.0d-2
  REAL(RealKind), PARAMETER :: Rho0=1.0d0
  REAL(RealKind) :: th0=293.16d0
  REAL(RealKind), PARAMETER :: D0=10.0d0
  REAL(RealKind) :: qvr=0.60d0
  REAL(RealKind) :: uMax=0.0d0,vMax=0.0d0
  REAL(RealKind) :: xQ0,xQ1,yQ0,yQ1,zQ0,zQ1
  REAL(RealKind) :: qcValue
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-6

  NAMELIST /Example/ uMax &
                    ,vMax &
                    ,xQ0  &
                    ,xQ1  &
                    ,yQ0  &
                    ,yQ1  &
                    ,zQ0  &
                    ,zQ1  &
                    ,qcValue &
                    ,TkeMax &
                    ,DisMax &
                    ,qvr &
                    ,th0

END MODULE Turm_Mod

SUBROUTINE InputExample(FileName)
  USE Turm_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: InputUnit,Pos
  CHARACTER(300) :: Line

! Find line
  InputUnit=1
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
  WRITE(*,*) 'uMax=',uMax
END SUBROUTINE InputExample

FUNCTION QcEmission(x,y,z,Time)
  USE Turm_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcEmission
  REAL(RealKind) :: x,y,z,Time
  QcEmission=Zero
  IF (x>=xQ0.AND.x<=xQ1.AND. &
      y>=yQ0.AND.y<=yQ1.AND. &
      z>=zQ0.AND.z<=zQ1) THEN
     QcEmission=QcValue
  END IF
END FUNCTION QcEmission


FUNCTION RhoFun(x,y,z)
  USE Turm_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z
  REAL(RealKind) :: S
  S=N*N/Grav
  RhoFun=Rho0*EXP(-S*z)*(One-Grav/(cp*th0*S)*(One-EXP(-S*z)))**(cv/Rd)
END FUNCTION RhoFun

FUNCTION ThProfFun(x,y,z)
  USE Turm_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z
  REAL(RealKind) :: S
  S=N*N/Grav
  ThProfFun=th0*exp(z*S)*RhoFun(x,y,z)
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z)
  USE Turm_Mod
  USE Rho_Mod 
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z
  REAL(RealKind) :: RhoLoc,ThLoc
  RhoLoc=RhoFun(x,y,z)
  ThLoc=thProfFun(x,y,z)
  QvProfFun=qvr*qvs(RhoLoc,ThLoc,ThLoc)
END FUNCTION QvProfFun

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Turm_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: U0
  UStart=uMax*RhoFun(x,y,z)
END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Turm_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Turm_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=vMax*RhoFun(x,y,z)
END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Turm_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,Time)
  USE Turm_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,Time
  WStart=0.0d0*RhoFun(x,y,z)
END FUNCTION WStart

FUNCTION ThStart(x,y,z,Time)
  USE Turm_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,Time
  ThStart=ThProfFun(x,y,z)
END FUNCTION ThStart

FUNCTION TkeStart(x,y,z,Time)
  USE Turm_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,Time
  TkeStart=TkeMax*RhoFun(x,y,z)
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,Time)
  USE Turm_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,Time
  DisStart=DisMax*RhoFun(x,y,z)
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,Time)
  USE Turm_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,Time
  QvStart=QvProfFun(x,y,z)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,Time)
  USE Turm_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QiStart(x,y,z,Time)
  USE Turm_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,Time
  QiStart=0.0d0
END FUNCTION QiStart


FUNCTION QrStart(x,y,z,Time)
  USE Turm_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE Turm_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=D0*RhoFun(x,y,z)
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,Time)
  USE Turm_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,Time
  DummyStart=0.0d0
END FUNCTION DummyStart


