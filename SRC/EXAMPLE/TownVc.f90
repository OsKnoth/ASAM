MODULE TownVc_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind), PARAMETER :: VelMax=1.0d0
  REAL(RealKind), PARAMETER :: N=1.0d-2
  REAL(RealKind), PARAMETER :: Rho0=1.0d0
  REAL(RealKind), PARAMETER :: th0=290d0
  REAL(RealKind), PARAMETER :: D0=1.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind), PARAMETER :: L=10000.0d0
  REAL(RealKind), PARAMETER :: h=5000.0d0
  REAL(RealKind), PARAMETER :: x0=0.0e0
  REAL(RealKind), PARAMETER :: y0=L
  REAL(RealKind), PARAMETER :: z0=0.5d0*h
  REAL(RealKind), PARAMETER :: rx=L
  REAL(RealKind), PARAMETER :: ry=L
  REAL(RealKind), PARAMETER :: rz=h
  REAL(RealKind) :: uMax=0.0d0,vMax=0.0d0
  REAL(RealKind) :: DInit=1.0d0
  REAL(RealKind) :: TkeMax=0.0d0,DisMax=0.0d0
  REAL(RealKind) :: Froude=.5d0
  REAL(RealKind) :: Height=5000.0d0
  NAMELIST /Example/ uMax &
                    ,vMax &
                    ,TkeMax &
                    ,DisMax &
                    ,DInit


END MODULE TownVc_Mod

SUBROUTINE InputExample(FileName)
  USE TownVc_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  CHARACTER(300) :: Line

! Find line
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
  if(MyID==0)WRITE(*,*) 'uMax=',uMax
  if(MyID==0)WRITE(*,*) 'vMax=',vMax
END SUBROUTINE InputExample

FUNCTION UStart(x,y,z,zHeight,Time)

  USE ReadProfile_Mod
  USE Rho_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,a,zHeight,Time,one=1.0
  REAL(RealKind), POINTER, SAVE :: uInt(:,:)
  CHARACTER*10, SAVE :: uType
  LOGICAL, SAVE :: Load=.TRUE.

  a=1 !MIN(Time/50.0,one)
  
  IF (Load) THEN
    CALL ReadProfile(uInt,uType,'uProf')
    Load=.FALSE.
  END IF
  IF (uType=='Stretch') THEN
    UStart=a*RhoFun(x,y,z,zHeight,Time)*ProfileStretch(uInt,z,zHeight)
  ELSE IF (uType=='Equal') THEN
    UStart=a*RhoFun(x,y,z,zHeight,Time)*ProfileEqual(uInt,z) 
  END IF

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE ReadProfile_Mod
  USE Rho_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: vInt(:,:)
  CHARACTER*10, SAVE :: vType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Load) THEN
    CALL ReadProfile(vInt,vType,'vProf')
    Load=.FALSE.
  END IF
  IF (vType=='Stretch') THEN
    VStart=RhoFun(x,y,z,zHeight,Time)*ProfileStretch(vInt,z,zHeight)
  ELSE IF (vType=='Equal') THEN
    VStart=RhoFun(x,y,z,zHeight,Time)*ProfileEqual(vInt,z) 
  END IF

END FUNCTION VStart

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThStart=ThProfFun(x,y,z,zHeight,Time)
END FUNCTION ThStart


FUNCTION DStart(x,y,z,zHeight,Time)

  USE ReadProfile_Mod
  USE Rho_Mod
  USE TownVc_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: DType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Dinit<0.0d0) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,DType,'DProf')
      Load=.FALSE.
    END IF
!   DStart=RhoFun(x,y,z,zHeight,Time)*ProfileEqual(cInt,z)
    DStart=RhoFun(x,y,z,zHeight,Time)*ProfileStretch(cInt,z,zHeight)
  ELSE
    DStart=RhoFun(x,y,z,zHeight,Time)*DInit
  END IF

END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE TownVc_Mod
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

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE TownVc_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S
  S=N*N/Grav
  ThProfFun=th0*exp(z*S)*RhoFun(x,y,z,zHeight,Time)
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE TownVc_Mod
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

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE TownVc_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE TownVc_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0*RhoFun(x,y,z,zHeight,Time)
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=Zero
END FUNCTION RhoStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE TownVc_Mod
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
