MODULE Fedkiw_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  LOGICAL :: ProfIn=.FALSE.
 
  REAL(RealKind) :: xSod
  REAL(RealKind) :: RhoLSod,RhoRSod
  REAL(RealKind) :: pLSod,pRSod
  REAL(RealKind) :: uLSod,uRSod
  CHARACTER*40 :: Problem
  
  NAMELIST /Example/ xSod &
                    ,RhoLSod &
                    ,RhoRSod &
                    ,pLSod &
                    ,pRSod &
                    ,uLSod &
                    ,uRSod &
                    ,Problem


END MODULE Fedkiw_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Fedkiw_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

END SUBROUTINE SetBoundCells

SUBROUTINE PerturbProfile(VecC)

  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Parameter_Mod
  USE Floor_Mod
  USE Fedkiw_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: ThVirt(:,:,:,:)
  REAL(RealKind) :: RhoLoc,RhoNew,RhoVLoc,RhoLLoc,RhoDry,RhoDryNew
  REAL(RealKind) :: KappaLoc,pLoc,pvs,rt,rvs
  REAL(RealKind) :: ThDensLoc,ThDensNew,ThLoc,ThNew,TLoc,ThEquivLoc
  REAL(RealKind) :: xPLoc,zPloc,L,Delta

END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE Fedkiw_Mod
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

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Fedkiw_Mod
  USE Rho_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time

  RhoProf=0.0d0

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  USE Start_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time


  ThProfFun=0.0d0

END FUNCTION ThProfFun

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Fedkiw_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  SELECT CASE (Problem)
    CASE ('Sod')
      IF (x<=xSod) THEN
        RhoFun=RhoLSod
      ELSE
        RhoFun=RhoRSod
      END IF
    CASE ('Advection')
      RhoFun=One  
  END SELECT 

END FUNCTION RhoFun


FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  HeightFun=Zero

END FUNCTION HeightFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThProfFun

  REAL(RealKind) :: RhoLoc,pLoc,uLoc,TLoc

  SELECT CASE (Problem)
    CASE ('Sod')
      IF (x<=xSod) THEN
        RhoLoc=RhoLSod
        pLoc=pLSod
        uLoc=uLSod
      ELSE
        RhoLoc=RhoRSod
        pLoc=pRSod
        uLoc=uRSod
      END IF
      SELECT CASE (ThetaKind)
        CASE ('Density')
          TLoc=pLoc/(Rd*RhoLoc)
          ThStart=(p0/pLoc)**Kappa*TLoc
        CASE ('Energy')
          TLoc=pLoc/(Rd*RhoLoc)
          ThStart=Cvd*TLoc+Half*uLoc*uLoc+z*Grav
        CASE('PreEn')
          ThStart=pLoc/RhoLoc
        END SELECT
    CASE ('Advection')
      ThStart=One
  END SELECT 
END FUNCTION ThStart

FUNCTION EnStart(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThProfFun

  REAL(RealKind) :: RhoLoc,pLoc,uLoc,TLoc

  EnStart=Zero
  SELECT CASE (Problem)
    CASE ('Sod')
      IF (x<=xSod) THEN
        RhoLoc=RhoLSod
        pLoc=pLSod
        uLoc=uLSod
      ELSE
        RhoLoc=RhoRSod
        pLoc=pRSod
        uLoc=uRSod
      END IF
      SELECT CASE (ThetaKind)
        CASE('PreEn')
          TLoc=pLoc/(Rd*RhoLoc)
          EnStart=Cvd*TLoc+Half*uLoc*uLoc+z*Grav
        END SELECT
    CASE ('Advection')
      EnStart=One
  END SELECT 
END FUNCTION EnStart

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  qvStart=0.0d0

END FUNCTION QvStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  qvProfFun=0.0d0

END FUNCTION QvProfFun


FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time

 qcStart=0.0d0

END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=Zero

END FUNCTION QiStart

FUNCTION UStart(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: lam,phi

  SELECT CASE (Problem)
    CASE ('Advection')
      UStart=One
    CASE ('Sod')
      IF (x<=xSod) THEN
        Ustart=uLSod
      ELSE
        Ustart=uRSod
      END IF
    CASE DEFAULT
      UStart=Zero
  END SELECT 

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: lam,phi

  SELECT CASE (Problem)
    CASE DEFAULT
      VStart=Zero
  END SELECT 

END FUNCTION VStart



FUNCTION DStart(x,y,z,zHeight,Time)

  USE Fedkiw_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=Zero

END FUNCTION DStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  SELECT CASE (Problem)
    CASE ('Sod')
      DummyStart=Zero
    CASE ('Advection')
      IF (x>=.3d0.AND.x<=.5d0) THEN
        DummyStart=One
      ELSE 
        DummyStart=Zero
      END IF
  END SELECT
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=Zero
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=Zero
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=0.0d0
END FUNCTION TStart


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

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Kind_Mod
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  ForceRho=0.0d0
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Kind_Mod
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  ForceRho=0.0d0
END FUNCTION QvStartSoil


FUNCTION DummyStart1(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(x,y,z,zHeight,Time)
  USE Fedkiw_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

