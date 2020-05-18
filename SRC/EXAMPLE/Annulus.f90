MODULE Annulus_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: x0,y0,z0
  REAL(RealKind) :: rOut
  REAL(RealKind) :: rIn
  REAL(RealKind) :: z1
  REAL(RealKind) :: UMax=0.4d0
  REAL(RealKind) :: OmegaLoc=5.0d0
  REAL(RealKind) :: cS0Loc=1.0d0
  REAL(RealKind) :: rMean
  REAL(RealKind) :: ThetaIn
  REAL(RealKind) :: ThetaOut
  REAL(RealKind) :: Theta0
  REAL(RealKind) :: DMax
  REAL(RealKind) :: DragH
  REAL(RealKind) :: Ta

  CHARACTER*40 :: Problem 
  
  NAMELIST /Example/ x0        &
                    ,y0        &
                    ,rOut      &
                    ,rIn       &
                    ,z1        &
                    ,UMax      &
                    ,OmegaLoc  & 
                    ,cS0Loc  & 
                    ,Problem   &
                    ,rMean     &
                    ,ThetaIn   &
                    ,ThetaOut  &
                    ,Theta0    &
                    ,DMax      &
                    ,DragH     &
                    ,Ta



END MODULE Annulus_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Annulus_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

  REAL(RealKind) :: x,y,rLoc

  SELECT CASE (Problem)
    CASE ('AnnulusZW')
      IF (BoundCellLoc%yS>rMean) THEN
        BoundCellLoc%ThetaS=ThetaOut
      ELSE
        BoundCellLoc%ThetaS=ThetaIn
      END IF
      BoundCellLoc%DragH=DragH
      BoundCellLoc%DragM=Zero
    CASE ('AnnulusAir')
      x=BoundCellLoc%xS 
      y=BoundCellLoc%yS 
      rLoc=SQRT(x*x+y*y)
      IF (rLoc>rMean) THEN
        BoundCellLoc%ThetaS=ThetaOut
      ELSE
        BoundCellLoc%ThetaS=ThetaIn
      END IF
      BoundCellLoc%DragH=DragH
      BoundCellLoc%DragM=Zero
  END SELECT


END SUBROUTINE SetBoundCells


SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile


SUBROUTINE InputExample(FileName)
  USE Annulus_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

  REAL(RealKind) :: L

  rMean=0.05d0
  ThetaIn=285.0d0
  ThetaOut=291.0d0
  Theta0=288.0d0 
  z1=0.1d0
! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&Example')>0) THEN
      BACKSPACE(InputUnit)
      WRITE(*,*) 'read Exam'
      READ(InputUnit,NML=Example)
      WRITE(*,*) 'Ta',Ta
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
  Omega=OmegaLoc
  DiffMin=DMax
  cS0=cS0Loc
  Rho0=1.0d3
END SUBROUTINE InputExample

FUNCTION UStart(x,y,z,zHeight,Time)

  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: r
  SELECT CASE (Problem)
    CASE ('Annulus')
      UStart=-UMax*Pi*y
    CASE ('AnnulusZ')
      r=y
      UStart=-UMax*r
    CASE ('AnnulusZW')
      UStart=Zero
    CASE DEFAULT
      UStart=Zero
  END SELECT

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  SELECT CASE (Problem)
    CASE ('Annulus')
      VStart=UMax*Pi*x
    CASE ('AnnulusZ')
      VStart=Zero
    CASE ('AnnulusZW')
      VStart=Zero
    CASE DEFAULT
      VStart=Zero
  END SELECT

END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: rLoc

  SELECT CASE (Problem)
    CASE ('Annulus')
      ThProfFun=One
    CASE ('AnnulusZ')
      ThProfFun=One
    CASE ('AnnulusZW')
      rLoc=y
      rLoc=MAX(MIN(rLoc,rOut),rIn)
      ThProfFun=(ThetaIn*(rOut-rLoc)+ThetaOut*(rLoc-rIn))/(rOut-rIn)
      ThProfFun=Theta0
    CASE ('AnnulusAir')
      ThProfFun=Zero
  END SELECT
END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: rLoc

  SELECT CASE (Problem)
    CASE ('Annulus')
      ThStart=One
    CASE ('AnnulusZ')
      ThStart=One
    CASE ('AnnulusZW')
      rLoc=y
      rLoc=MAX(MIN(rLoc,rOut),rIn)
      IF (Time>0.0d0) THEN
        ThStart=(ThetaIn*(rOut-rLoc)+ThetaOut*(rLoc-rIn))/(rOut-rIn)
      ELSE  
        ThStart=Theta0
      END IF  
    CASE ('AnnulusAir')
      rLoc=SQRT(x*x+y*y)
      rLoc=MAX(MIN(rLoc,rOut),rIn)
      ThStart=(ThetaIn*(rOut-rLoc)+ThetaOut*(rLoc-rIn))/(rOut-rIn)
  END SELECT

END FUNCTION ThStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvStart=Zero

END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Annulus_Mod
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

  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=DMax

END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: pLoc,ThLoc,rLoc

  SELECT CASE (Problem)
    CASE ('Annulus')
      RhoFun=One
    CASE ('AnnulusZ')
      RhoFun=One
    CASE ('AnnulusZW')
      rLoc=y
      ThLoc=(ThetaIn*(rOut-rLoc)+ThetaOut*(rLoc-rIn))/(rOut-rIn)
      pLoc=p0*(One-kappa*Grav*z/(Rd*ThLoc))**(cpD/Rd)
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
      RhoFun=Rho0
    CASE ('AnnulusAir')
      rLoc=SQRT(x*x+y*y)
      ThLoc=(ThetaIn*(rOut-rLoc)+ThetaOut*(rLoc-rIn))/(rOut-rIn)
      pLoc=p0*(One-kappa*Grav*z/(Rd*ThLoc))**(cpD/Rd)
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END SELECT

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Annulus_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: pLoc,ThLoc

  SELECT CASE (Problem)
    CASE ('Annulus')
      RhoProf=One
    CASE ('AnnulusZ')
      RhoProf=One
    CASE ('AnnulusZW')
      pLoc=p0*(One-kappa*Grav*z/(Rd*Theta0))**(cpD/Rd)
      ThLoc=Theta0
      RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
      RhoProf=Rho0
    CASE ('AnnulusAir')
      ThLoc=Theta0
      pLoc=p0*(One-kappa*Grav*z/(Rd*ThLoc))**(cpD/Rd)
      RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
      RhoProf=Zero
  END SELECT

END FUNCTION RhoProf

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Annulus_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Annulus_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart1(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: phi

  SELECT CASE (Problem)
    CASE ('Annulus')
      phi=ATAN2(y,x)-0.5d0*Pi+Pi
      DummyStart1=0.5d0*(ERF(5.0d0*(Pi/6.0d0-Phi)) &
                       +ERF(5.0d0*(Pi/6.0d0+Phi)))
    CASE ('AnnulusZ')
      phi=x-Pi
      DummyStart1=0.5d0*(ERF(5.0d0*(Pi/6.0d0-Phi)) &
                       +ERF(5.0d0*(Pi/6.0d0+Phi)))
    CASE DEFAULT
      DummyStart1=Zero
  END SELECT
END FUNCTION DummyStart1

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart


FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=Zero
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=Zero
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  HeightFun=Zero

END FUNCTION HeightFun
FUNCTION DummyStart(lam,phi,z,zHeight,Time)
  USE Annulus_Mod
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart
FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE Annulus_Mod
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2
FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE Annulus_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3
FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE Annulus_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4
FUNCTION ForceU(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=Zero
END FUNCTION ForceU
FUNCTION ForceV(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=Zero
END FUNCTION ForceV
FUNCTION ForceW(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=Zero
END FUNCTION ForceW
FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=0.0d0
END FUNCTION ForceRho

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  EnStart=Zero
END FUNCTION EnStart

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass,SoilType
  REAL(RealKind) :: Time
  ThStartSoil=Zero
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass,SoilType
  REAL(RealKind) :: Time
  QvStartSoil=Zero
END FUNCTION QvStartSoil

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

FUNCTION OmeStart(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  OmeStart=0.0d0
END FUNCTION OmeStart


FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE Annulus_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=Zero
END FUNCTION ForceTh
