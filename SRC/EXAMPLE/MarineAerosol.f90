MODULE MarineAerosol_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE SoilData_Mod

  IMPLICIT NONE 

! NAMELIST-Variablen
  REAL(RealKind) :: uMax=10.0d0
  REAL(RealKind) :: vMax=0.0d0
  REAL(RealKind) :: DMax=0.0d0
  REAL(RealKind) :: th0=293.16d0
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d3
  REAL(RealKind) :: H=5.d2 ! Altitude of Inversionlayer
  REAL(RealKind) :: z0Water=1.0d-4
  LOGICAL :: ProfIn=.FALSE.
  LOGICAL :: Perturb=.FALSE.
  CHARACTER(12) :: ProfileWind='Logarithmic'

! Sonstige Variablen
  REAL(RealKind) :: LapseRateTh=0.01d0

  NAMELIST /Example/ uMax,    &
                     vMax,    &
                     DMax,    &
                     th0,     &
                     TkeHMax, &
                     TkeVMax, &
                     LenMax,  &
                     H,       &
                     z0Water, &
                     ProfIn,  &
                     Perturb, &
                     ProfileWind

END MODULE MarineAerosol_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)
  USE MarineAerosol_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc
  BoundCellLoc%TeS=298.d0
  BoundCellLoc%ThetaS=298.d0
  BoundCellLoc%LandClass=9
  BoundCellLoc%zRauh=z0Water
  BoundCellLoc%zRauhT=z0Water
  BoundCellLoc%qv=0.0d0
  BoundCellLoc%qv=0.75*SaturVapor(BoundCellLoc%TeS)/(Rv*BoundCellLoc%TeS)
END SUBROUTINE SetBoundCells

SUBROUTINE PerturbProfile(VecC)
  USE MarineAerosol_Mod
  USE DataType_Mod
  USE Floor_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: rRand
  IF (Perturb) THEN
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        IF (zP(iz)<80.0d0) THEN
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1
              CALL Random_Number(rRand)
              VecC(ibLoc)%Vec(ThPos)%c(ix,iy,iz,1)=VecC(ibLoc)%Vec(ThPos)%c(ix,iy,iz,1)+1.0d-1*Two*(rRand-0.5d0)
            END DO
          END DO
        END IF
      END DO
    END DO
  END IF
END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
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

FUNCTION UStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: ustar
  REAL(RealKind) :: x,y,z,zHeight,Time
  SELECT CASE(ProfileWind)
    CASE('Const')
      UStart=uMax
    CASE('Logarithmic')
      uStar=uMax*Karm/LOG(H/z0Water+One)
      UStart=uStar/Karm*LOG(z/z0Water+One)
  END SELECT
END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: ustar
  REAL(RealKind) :: x,y,z,zHeight,Time
  SELECT CASE(ProfileWind)
    CASE('Const')
      VStart=uMax
    CASE('Logarithmic')
      uStar=VMax*Karm/LOG(H/z0Water+One)
      VStart=uStar/Karm*LOG(z/z0Water+One)
  END SELECT
END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER(10), SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  ThProfFun=Zero
END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE Rho_Mod
  USE QvProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER(10), SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S
  IF (ProfIn) THEN
    IF (Load) THEN
      SELECT CASE(ThetaKind)
      CASE('Density')
        CALL ReadProfile(cInt,thType,'ThDensProf')
      CASE('Equiv')
        CALL ReadProfile(cInt,thType,'ThetaEProf')
      CASE('Energy')
        CALL ReadProfile(cInt,thType,'EnergyProf')
      CASE('PreEn')
        CALL ReadProfile(cInt,thType,'PreProf')
      CASE DEFAULT
        WRITE(*,*) 'ThetaKind pruefen.'
      END SELECT
      Load=.FALSE.
    END IF
    ThStart=ProfileEqual(cInt,z)
  ELSE
    ThStart=Th0+LapseRateTh*z
  END IF
END FUNCTION ThStart


FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  !USE Parameter_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER(10), SAVE :: RhoType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: pLoc,ThLoc
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,RhoType,'RhoProf')
      Load=.FALSE.
    END IF
    RhoFun=ProfileEqual(cInt,z)
  ELSE
    ThLoc=Th0+LapseRateTh*z
    pLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
    RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END IF
END FUNCTION RhoFun

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER(10), SAVE :: tType
  LOGICAL, SAVE :: Load=.TRUE.
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,tType,'PreProf')
      Load=.FALSE.
    END IF
    PreStart=ProfileEqual(cInt,z)
  ELSE
    PreStart=p0*(One-kappa*Grav*z/(Rd*Th0))**(Cpd/Rd)
  END IF
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=Zero
END FUNCTION PreProf

FUNCTION TFun(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE Rho_Mod
  USE PreProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,RhoLoc,GammaDry
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER(10), SAVE :: tType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,tType,'TempProf')
      Load=.FALSE.
    END IF
    TFun=ProfileEqual(cInt,z)
  ELSE
    TFun=300.0d0
  END IF
END FUNCTION TFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: TLoc,RhoLoc
  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE QvProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER(10), SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.
  !IF (ProfIn) THEN
    !IF (Load) THEN
      !CALL ReadProfile(cInt,qvType,'RhoVProf')
      !Load=.FALSE.
    !END IF
    !QvStart=ProfileEqual(cInt,z)
  !ELSE
    QvStart=1.2d-2
  !END IF
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER(10), SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.
  !IF (ProfIn) THEN
    !IF (Load) THEN
      !CALL ReadProfile(cInt,qcType,'QcProf')
      !Load=.FALSE.
    !END IF
    !qcStart=ProfileEqual(cInt,z)
  !ELSE
    QcStart=Zero
  !END IF
END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=Zero
END FUNCTION QiStart

FUNCTION QsStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QsStart=Zero
END FUNCTION QsStart

FUNCTION NvStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NvStart=Zero
END FUNCTION NvStart

FUNCTION NcStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NcStart=Zero
END FUNCTION NcStart

FUNCTION NrStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NrStart=Zero
END FUNCTION NrStart

FUNCTION NiStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NiStart=Zero
END FUNCTION NiStart

FUNCTION NsStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NsStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NsStart=Zero
END FUNCTION NsStart

FUNCTION OmeStart(lam,phi,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  OmeStart=Zero
END FUNCTION OmeStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=DMax
END FUNCTION DStart

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoProf=Zero
END FUNCTION RhoProf

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=Zero
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=1.0d-2
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=Zero
END FUNCTION QrStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=Zero
END FUNCTION TStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=Zero
END FUNCTION DummyStart

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=Zero
END FUNCTION HeightFun

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=Zero
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=Zero
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=Zero
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE Rho_Mod
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: fTime,fSpace
  ForceRho=Zero
END FUNCTION ForceRho

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE Rho_Mod
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=Zero
END FUNCTION ForceTh

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE MarineAerosol_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  ThStartSoil=300.0d0
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE MarineAerosol_Mod
  !USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  QvStartSoil = (cfcap(SoilType) - cpwp(SoilType))*50.0/100.0
END FUNCTION QvStartSoil

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE MarineAerosol_Mod
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1=1.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE MarineAerosol_Mod
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE MarineAerosol_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE MarineAerosol_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE MarineAerosol_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  EnStart=Zero
END FUNCTION EnStart

FUNCTION DampFun(z,Name)
  USE MarineAerosol_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
  IF (Name=='VelProf') THEN
    IF (z<=825.0d0) THEN
      DampFun=1.0d0/7200.0d0*(1.0d0-COS(Pi*z/825.0d0))/2.0d0
    ELSE
      DampFun=1.0d0/7200.0d0
    END IF
  ELSE IF (Name=='ThProf'.OR.Name=='RhoVProf') THEN
    IF (z<=1200.0d0) THEN
      DampFun=0.0d0
    ELSE IF (z<=1500.0d0) THEN
      DampFun=1.0d0/3600.0d0*(1.0d0-COS(Pi*(z-1200.0d0)/(1500.0d0-1200.0d0)))/2.0d0
    ELSE
      DampFun=1.0d0/3600.0d0
    END IF
  END IF
END FUNCTION DampFun
