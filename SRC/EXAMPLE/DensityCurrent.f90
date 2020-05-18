MODULE DensityCurrent_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: x0,y0,z0
  REAL(RealKind) :: x0Qc,y0Qc,z0Qc
  REAL(RealKind) :: xr,yr,zr
  REAL(RealKind) :: DeltaTh,uMax,vMax,ThInit,DMax
  REAL(RealKind) :: DeltaQc
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0
  CHARACTER*20 :: TypePert='Temp2D'
  
  NAMELIST /Example/ x0 &
                    ,y0 &
                    ,z0 &
                    ,x0Qc &
                    ,y0Qc &
                    ,z0Qc &
                    ,xr &
                    ,yr &
                    ,zr &
                    ,uMax &
                    ,vMax &
                    ,DMax &
                    ,ThInit &
                    ,DeltaTh &
                    ,DeltaQc &
                    ,TypePert


END MODULE DensityCurrent_Mod

SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile


SUBROUTINE InputExample(FileName)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

! Find line
  x0=Zero
  x0Qc=Zero
  xr=4000.0d0
  y0=Zero
  y0Qc=Zero
  yr=4000.0d0
  z0=3000.0d0
  z0Qc=3000.0d0
  zr=2000.0d0
  DMax=75.0d0
  uMax=0.0d0
  vMax=0.0d0
  ThInit=300.0d0
  DeltaTh=-15.0d0
  DeltaQc=2.0d0
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
  WRITE(*,*) x0,z0
  WRITE(*,*) DeltaTh
END SUBROUTINE InputExample

FUNCTION UStart(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  UStart=uMax

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  VStart=VMax

END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  ThProfFun=0.0d0
END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad,pLoc

  SELECT CASE (TypePert)
    CASE ('Temp2D')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      Rad=SQRT(((x-x0)/xr)**2 &
              +((z-z0)/zr)**2 &
                ) 
      ThStart=ThInit
      IF (Rad<1.0d0) THEN
        ThStart=ThStart+DeltaTh*(COS(Pi*Rad)+1.0d0)/2.0d0*(pLoc/p0)**(-kappa) 
      END IF
    CASE ('Temp3D')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      Rad=SQRT(((x-x0)/xr)**2 &
              +((y-y0)/yr)**2 &
              +((z-z0)/zr)**2 &
            ) 
      ThStart=ThInit
      IF (Rad<1.0d0) THEN
        ThStart=ThStart+DeltaTh*(COS(Pi*Rad)+1.0d0)/2.0d0*(pLoc/p0)**(-kappa) 
      END IF
    CASE ('PotTemp2D')
      Rad=SQRT(((x-x0)/xr)**2 &
              +((z-z0)/zr)**2 &
                ) 
      ThStart=ThInit
      IF (Rad<1.0d0) THEN
        ThStart=ThStart+DeltaTh*(COS(Pi*Rad)+1.0d0)/2.0d0
      END IF
   END SELECT   

END FUNCTION ThStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvStart=Zero

END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad

  Rad=SQRT(((x-x0Qc)/xr)**2 &
          +((y-y0Qc)/yr)**2 &
          +((z-z0Qc)/zr)**2 &
            ) 
  QcStart=Zero
  IF (Rad<1.0d0) THEN
    QcStart=QcStart+DeltaQc*(COS(Pi*Rad)+1.0d0)/2.0d0 
  END IF


END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=Zero

END FUNCTION QiStart


FUNCTION DStart(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=DMax

END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,ThLoc,Rad

  SELECT CASE (TypePert)
    CASE ('Temp2D')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      Rad=SQRT(((x-x0)/xr)**2 &
              +((z-z0)/zr)**2 &
                )
      ThLoc=ThInit
      IF (Rad<1.0d0) THEN
        ThLoc=ThLoc+DeltaTh*(COS(Pi*Rad)+1.0d0)/2.0d0*(pLoc/p0)**(-kappa)
      END IF
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
    CASE ('Temp3D')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      Rad=SQRT(((x-x0)/xr)**2 &
              +((y-y0)/yr)**2 &
              +((z-z0)/zr)**2 &
                )
      ThLoc=ThInit
      IF (Rad<1.0d0) THEN
        ThLoc=ThLoc+DeltaTh*(COS(Pi*Rad)+1.0d0)/2.0d0*(pLoc/p0)**(-kappa)
      END IF
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
    CASE ('PotTemp2D')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      Rad=SQRT(((x-x0)/xr)**2 &
              +((z-z0)/zr)**2 &
                )
      ThLoc=ThInit
      IF (Rad<1.0d0) THEN
        ThLoc=ThLoc+DeltaTh*(COS(Pi*Rad)+1.0d0)/2.0d0
      END IF
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END SELECT   


END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE DensityCurrent_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc

  RhoProf=Zero

END FUNCTION RhoProf


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart


FUNCTION qsStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QsStart=0.d0
END FUNCTION QsStart

FUNCTION NcStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NcStart=0.d0
END FUNCTION NcStart

FUNCTION NrStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NrStart=0.d0
END FUNCTION NrStart

FUNCTION NiStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NiStart=0.d0
END FUNCTION NiStart


FUNCTION OmeStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  OmeStart=0.d0
END FUNCTION OmeStart

FUNCTION NsStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NsStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NsStart=0.d0
END FUNCTION NsStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart


FUNCTION PreStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
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

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=0.0d0
END FUNCTION ForceTh

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE DensityCurrent_Mod
  USE Start_Mod
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  ThAir=ThStart(x,y,z,zHeight,Time)
  IF (zSoil<=0.05d0) THEN
    ThStartSoil=ThAir-0.5d0
  ELSE IF (zSoil<=0.5d0) THEN
    ThStartSoil=ThAir-1.0d0
  ELSE IF (zSoil<=1.0d0) THEN
    ThStartSoil=ThAir-2.0d0
  ELSE
    ThStartSoil=ThAir-3.5d0
  END IF
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE DensityCurrent_Mod
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  QvStartSoil=0.d0 
END FUNCTION QvStartSoil

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=0.0d0
END FUNCTION Tracer1Start
FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=0.0d0
END FUNCTION Tracer2Start

FUNCTION DummyStart1(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  EnStart=0.0d0
END FUNCTION EnStart


SUBROUTINE SetBoundCells(BoundCellLoc)

  USE DensityCurrent_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

END SUBROUTINE SetBoundCells

FUNCTION HeightFun(x,y,z,zHeight,Time)

  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  HeightFun=Zero

END FUNCTION HeightFun

FUNCTION NvStart(lam,phi,z,zHeight,Time)
  USE DensityCurrent_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NvStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  NvStart = 0.0
END FUNCTION NvStart

FUNCTION DampFun(z,Name)
  USE DensityCurrent_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
END FUNCTION DampFun
