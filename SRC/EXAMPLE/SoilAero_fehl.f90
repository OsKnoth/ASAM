MODULE SauleBoden_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: NL=1.0d-2
  REAL(RealKind) :: NU=1.0d-2
  REAL(RealKind) :: H=1.0d-2
  REAL(RealKind) :: RhoL=1.15d0
  REAL(RealKind) :: RhoU=0.0d0
  REAL(RealKind) :: D0=10.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind) :: xL1,xL2,xL3,xL4

  REAL(RealKind) :: th0=293.16d0
  REAL(RealKind) :: N=0.0d0
  REAL(RealKind) :: ThSoil=287.0d0
  REAL(RealKind) :: uMax=8.5d0,vMax=0.0d0
  REAL(RealKind) :: QvSoil=0.2 !1.0d0
  
  REAL(RealKind) :: rx1,ry1,rz1
  REAL(RealKind) :: x0Qc,y0Qc,z0Qc
  REAL(RealKind) :: ThInit,DMax
  REAL(RealKind) :: DeltaQc
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-4
  REAL(RealKind) :: TkeHMax=1.0d-2,TkeVMax=1.0d-2,LenMax=1.0d3
  LOGICAL :: ProfIn=.FALSE.

  NAMELIST /Example/ th0, &
                    NL,    &
                    NU,    &
                    D0,    &
                    H,     &
                    xL1,   &
                    xL2,   &
                    xL3,   &
                    xL4,   &
                    N,   &
                    ThSoil,   &
                    uMax, &
                    vMax, &
                    QvSoil,   &
                    x0Qc, &
                    y0Qc, &
                    z0Qc, &
                    rx1, &
                    ry1, &
                    rz1, &
                    DMax, &
                    ThInit, &
                    ProfIn


CONTAINS


FUNCTION ThetaLoc(z)

  REAL(RealKind) :: ThetaLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: SL,SU

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    ThetaLoc=th0*EXP(SL)
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
    ThetaLoc=th0*EXP(SL+SU)
  END IF

END FUNCTION ThetaLoc

FUNCTION PressLoc(z)

  REAL(RealKind) :: PressLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: SL,SU
  REAL(RealKind) :: t1,t2,pH
  REAL(RealKind) :: s1,s2

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    t1=Grav*Grav*kappa/(Rd*th0*NL*NL)
    t2=Grav*Grav*kappa/(Rd*th0*NL*NL)*EXP(-SL)
    PressLoc=p0 &
            *(1.0d0+t2-t1)**(1.0d0/kappa)
           
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
    t1=Grav*Grav*kappa/(Rd*th0*NL*NL)
    t2=Grav*Grav*kappa/(Rd*th0*NL*NL)*EXP(-SL)
    pH=p0 &
       *(1.0d0+t2-t1)**(1.0d0/kappa)
    s1=Grav*Grav*kappa/(Rd*th0*NU*NU)*EXP(-SL)
    s2=Grav*Grav*kappa/(Rd*th0*NU*NU)*EXP(-SL-SU)
    PressLoc=p0 &
       *(1.0d0+t2-t1+s2-s1)**(1.0d0/kappa)
  END IF

END FUNCTION PressLoc


END MODULE SauleBoden_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)
  USE SauleBoden_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc
END SUBROUTINE SetBoundCells

SUBROUTINE PerturbProfile(VecC)
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)
END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE SauleBoden_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
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
END SUBROUTINE InputExample

FUNCTION UStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStart=UMax
END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
!  VStart=Zero
  VStart=vMax
END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThProfFun=Zero
END FUNCTION ThProfFun

!========================================
FUNCTION ThStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThProfFun
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,thType,'ThProf')
      Load=.FALSE.
    END IF
    ThStart=ProfileEqual(cInt,z)
  ELSE
    IF (z<=700) THEN
      ThStart=thInit
    ELSE
      ThStart=thInit+(z-700)*0.006
    END IF
  END IF
END FUNCTION ThStart

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil
  INTEGER :: LandClass,SoilType
  REAL(RealKind) :: Time
!  ThStartSoil=ThSoil
  IF (zSoil<=0.1d0) THEN
!    ThStartSoil=276.50d0
!    ThStartSoil=290.00d0
    ThStartSoil=289.55d0
  ELSE IF (zSoil<=0.5d0) THEN
!    ThStartSoil=277.70d0
!    ThStartSoil=290.20d0
    ThStartSoil=289.75d0
  ELSE IF (zSoil<=1.0d0) THEN
!    ThStartSoil=280.90d0
!    ThStartSoil=289.90d0
    ThStartSoil=289.45d0
  ELSE
!    ThStartSoil=286.50d0
!    ThStartSoil=288.50d0
    ThStartSoil=288.50d0
  END IF
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil
  INTEGER :: LandClass,SoilType
  REAL(RealKind) :: Time
  IF (zSoil<=0.1d0) THEN
    QvStartSoil=2.d-1
  ELSE IF (zSoil<=1.0d0) THEN
    QvStartSoil=1.8d-1
  ELSE
    QvStartSoil=1.2d-1
  END IF
END FUNCTION QvStartSoil

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
!  USE QvProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.
  IF (ProfIn) THEN
    IF (Load) THEN
!      CALL ReadProfile(cInt,qvType,'QvProf')
      CALL ReadProfile(cInt,qvType,'RhoVProf')
      Load=.FALSE.
    END IF
    qvStart=ProfileEqual(cInt,z)
  ELSE
!    IF (z<=700) THEN
!      qvStart=8d-3
!    ELSE
!      qvStart=4d-3
!    END IF  
     qvStart=0.0d0
  END IF
END FUNCTION QvStart
!========================================

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: Rad
  QcStart=Zero
END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=Zero
END FUNCTION QiStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=DMax
END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE Parameter_Mod
  USE SauleBoden_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: RhoType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S,pLoc,ThLoc,L
  
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,RhoType,'RhoProf')
      Load=.FALSE.
    END IF
    RhoFun=ProfileEqual(cInt,z)
  ELSE
    S=N*N/Grav
    ThLoc=ThInit*exp(z*S)
    IF (N>Zero) THEN
      pLoc=p0*(One-Grav/(Cpd*ThInit*S)*(One-EXP(-S*z)))**(Cpd/Rd)
    ELSE
      pLoc=p0*(One-kappa*Grav*z/(Rd*ThInit))**(Cpd/Rd)
    END IF
    RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END IF
END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE Parameter_Mod
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc
  RhoProf=Zero
END FUNCTION RhoProf

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE SauleBoden_Mod
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE SauleBoden_Mod
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE SauleBoden_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE SauleBoden_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc,Fpot,RhoLoc
  ThLoc  =ThetaLoc(z)
  pLoc   =PressLoc(z)
  RhoLoc =pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  Fpot   = (Rd*RhoLoc*ThLoc/p0)**(-kappa/(One-kappa))
  TStart =ThLoc/Fpot
!  TStart =0.0d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=Zero
END FUNCTION HeightFun

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=0.0d0
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=0.0d0
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=0.0d0
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=0.0d0
END FUNCTION ForceRho


