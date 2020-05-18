MODULE TestAdv1D_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: x0,y0,z0
  REAL(RealKind) :: UMax
  
  NAMELIST /Example/ x0 &
                    ,y0 &
                    ,UMax


END MODULE TestAdv1D_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE TestAdv1D_Mod
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
  USE TestAdv1D_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

! Find line
  x0=0.5d0
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

  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  UStart=UMax

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  VStart=Zero

END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  ThProfFun=Zero
END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad,pLoc

  ThStart=One

END FUNCTION ThStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvStart=Zero

END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE TestAdv1D_Mod
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

  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=Zero

END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  RhoFun=One

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE TestAdv1D_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time

  RhoProf=Zero

END FUNCTION RhoProf


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: phi

  DummyStart=ABS(COS(Pi*(x-x0)))**10
  IF (ABS(x-x0)<=2500.0d0) THEN
    DummyStart=1.0d0
  ELSE
    DummyStart=0.0d0
  END IF
 
END FUNCTION DummyStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart


FUNCTION PreStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=Zero
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=Zero
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE TestAdv1D_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun =0.0d0
END FUNCTION HeightFun





