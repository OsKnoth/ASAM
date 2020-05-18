MODULE SauleBoden_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE SoilData_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: th0=300.0d0
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: ThSoil=287.0d0
  REAL(RealKind) :: uMax=5.0d0
  REAL(RealKind) :: vMax=5.0d0
  REAL(RealKind) :: QvSoil=0.2 !1.0d0
  
  REAL(RealKind) :: rx1,ry1,rz1
  REAL(RealKind) :: x0Qc,y0Qc,z0Qc
  REAL(RealKind) :: ThInit,DMax
  REAL(RealKind) :: DeltaQc
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0
  LOGICAL :: ProfIn=.FALSE.
  LOGICAL :: WangProfil=.FALSE.
  CHARACTER*40 :: SauleCase=''

  NAMELIST /Example/ th0 &
                    ,N   &
                    ,ThSoil   &
                    ,uMax &
                    ,vMax &
                    ,QvSoil   &
                    ,x0Qc &
                    ,y0Qc &
                    ,z0Qc &
                    ,rx1 &
                    ,ry1 &
                    ,rz1 &
                    ,DMax &
                    ,ThInit &
                    ,ProfIn & 
                    ,WangProfil &
                    ,SauleCase &
                    ,SoilParam

END MODULE SauleBoden_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)
  USE SauleBoden_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

  BoundCell%qv=0.0d0
  BoundCell%TeS=ThInit
END SUBROUTINE SetBoundCells

SUBROUTINE PerturbProfile(VecC)
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)
END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE SauleBoden_Mod
  USE SoilData_Mod
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
  SELECT CASE(SauleCase)
    CASE('3Stripes')
      zrough(0:9)=0.05d0
  END SELECT
END SUBROUTINE InputExample

FUNCTION UStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: uType
  LOGICAL, SAVE :: Load=.TRUE.
 
  SELECT CASE(SauleCase)
  CASE('Valley')
    UStart=uMax
  CASE DEFAULT
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,uType,'uProf')
      Load=.FALSE.
    END IF
    UStart=ProfileEqual(cInt,z)
  ELSE IF (WangProfil) THEN
    UStart=MIN(uMax,uMax*LOG((z/0.01d0))/LOG((1000.0d0/0.01d0)))
  ELSE
    UStart=UMax
  END IF
  END SELECT
END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: vType
  LOGICAL, SAVE :: Load=.TRUE.
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,vType,'vProf')
      Load=.FALSE.
    END IF
    VStart=ProfileEqual(cInt,z)
  ELSE IF (WangProfil) THEN
    VStart=MIN(vMax,vMax*LOG((z/0.01d0))/LOG((1000.0d0/0.01d0)))
  ELSE
    VStart=vMax
  END IF
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
  ELSE IF (WangProfil) THEN
    IF (z<=420.0d0) THEN
      ThStart=289.2d0
    ELSE IF (z<=550.0d0) THEN
      ThStart=289.2d0+1.d-2*(z-420.0d0)
    ELSE
      ThStart=290.5+4.21d-3*(z-550.0d0)
    END IF
  ELSE
    S=N*N/Grav
    ThStart=thInit*exp(z*S)
  END IF
END FUNCTION ThStart

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE SauleBoden_Mod
  USE Start_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil,ThAir
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  REAL(RealKind) :: Rm,Cpml,KappaLoc
  REAL(RealKind) :: RhoDLoc,RhoLoc,RhoVLoc,RhoLLoc
  RhoLoc=RhoFun(x,y,z,zHeight,Time)
  RhoVLoc=QvStart(x,y,z,zHeight,Time)
  RhoLLoc=QcStart(x,y,z,zHeight,Time)
  RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
  Rm=Rd*RhoDLoc+Rv*RhoVLoc
  Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc
  KappaLoc=Rm/Cpml
  ThAir=ThStart(x,y,z,zHeight,Time)
  ThAir=(Rd*RhoLoc*ThAir/p0**KappaLoc)**(One/(One-KappaLoc))/(Rd*RhoDLoc+Rv*RhoVLoc)
  IF (zSoil<=0.05d0) THEN
!    ThStartSoil=280.50d0
    ThStartSoil=ThAir
  ELSE IF (zSoil<=0.5d0) THEN
    ThStartSoil=ThAir-0.5d0
!    ThStartSoil=280.00d0
  ELSE IF (zSoil<=1.0d0) THEN
    ThStartSoil=ThAir-2.0d0
!    ThStartSoil=279.50d0
  ELSE
    ThStartSoil=ThAir-3.5d0
!    ThStartSoil=278.00d0
  END IF
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  IF (zSoil==0.0d0) THEN
 !  QvStartSoil=0.d0  ! interception reservoir
 !  QvStartSoil=MIN(QvStartSoil,5.d-4)
    QvStartSoil=5.d-4
  ELSE 
 !  QvStartSoil=0.0d0 ! 0.2
 !  QvStartSoil=MIN(QvStartSoil,cporv(SoilType))
    QvStartSoil = (cfcap(SoilType) - cpwp(SoilType))*50.0/100.0  ! = 50 % plant-available water
  END IF
  IF (SoilType<=2) QvStartSoil=0.0d0
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
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qvType,'RhoVProf')
      Load=.FALSE.
    END IF
    qvStart=ProfileEqual(cInt,z)
  ELSE IF (WangProfil) THEN
    IF (z<=450.0d0) THEN
      QvStart=MAX(6.d-3,6.0d-3+0.002d0*EXP(-z/100.d0))
    ELSE IF (z<=550.0d0) THEN
      QvStart=MAX(0.0d0,6.d-3-6.d-5*(z-450.d0))
    ELSE
      QvStart=0.0d0
    END IF
  ELSE
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

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=Zero
END FUNCTION QiStart

FUNCTION QsStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  QsStart=Zero
END FUNCTION QsStart

FUNCTION NvStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NvStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  NvStart=Zero
END FUNCTION NvStart

FUNCTION NcStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NcStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  NcStart=Zero
END FUNCTION NcStart

FUNCTION NrStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NrStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  NrStart=Zero
END FUNCTION NrStart

FUNCTION NiStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NiStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  NiStart=Zero
END FUNCTION NiStart

FUNCTION NsStart(x,y,z,zHeight,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NsStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  NsStart=Zero
END FUNCTION NsStart

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
  WStart=Zero
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ! TkeStart=Zero
  TkeStart=0.01d0
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION OmeStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ! OmeStart=Zero
  OmeStart=0.01d0
END FUNCTION OmeStart

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


FUNCTION DummyStart(lam,phi,z,zHeight,Time)
  USE SauleBoden_Mod
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

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
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: preType
  LOGICAL, SAVE :: Load=.TRUE.
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,preType,'PreProf')
      Load=.FALSE.
    END IF
    PreStart=ProfileEqual(cInt,z)
  ELSE
    PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
  END IF
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
  TStart =0.0d0
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

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=0.0d0
END FUNCTION ForceTh

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  EnStart=Zero
END FUNCTION EnStart

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE SauleBoden_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start

FUNCTION DampFun(z,Name)
  USE SauleBoden_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
END FUNCTION DampFun
