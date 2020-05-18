MODULE DensityBryan_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: H0=31000.0d0
  REAL(RealKind) :: x1=-5000.0d0
  REAL(RealKind) :: x2=5000.0d0
  REAL(RealKind) :: z1=5000.0d0
  REAL(RealKind) :: z2=10000.0d0
  REAL(RealKind) :: uMax=.0d0
  REAL(RealKind) :: ThInit=302.9d0
  REAL(RealKind) :: Ratio=0.9d0
  REAL(RealKind) :: B=-0.12d0
  INTEGER :: NumX=4000
  REAL(RealKind) :: H
  
  NAMELIST /Example/ uMax &
                    ,H0   &
                    ,ThInit & 
                    ,B &
                    ,x1 &
                    ,x2 & 
                    ,z1 &
                    ,z2 

CONTAINS

FUNCTION DeltaThF(x,z)
  REAL(RealKind) :: DeltaThF
  REAL(RealKind) :: x,z

  IF (z>=z1.AND.z<=z2.AND. &
      x>=x1.AND.x<=x2) THEN
    DeltaThF=B*ThInit/Grav
  ELSE
     DeltaThF=0.0d0
  END IF

END FUNCTION DeltaThF


END MODULE DensityBryan_Mod

SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile


SUBROUTINE InputExample(FileName)
  USE DensityBryan_Mod
  USE Control_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

  REAL(RealKind) :: Delta,Time

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
  H=Ratio*H0
  Delta=H/100.0d0
  Time=SQRT(H/ABS(B))
  EndTime=60.0d0*Time
  WRITE(*,*) 'H   =',H
  WRITE(*,*) 'Time=',Time
END SUBROUTINE InputExample

FUNCTION UStart(x,y,z,zHeight,Time)

  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  UStart=uMax

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  VStart=0.0d0

END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  ThProfFun=ThInit
  ThProfFun=0.0d0
END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad,pLoc
  REAL(RealKind) :: DeltaTh

  DeltaTh=DeltaThF(x,z)
  ThStart=ThInit+DeltaTh

END FUNCTION ThStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvStart=Zero

END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad

  QcStart=0.0d0

END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=Zero

END FUNCTION QiStart


FUNCTION DStart(x,y,z,zHeight,Time)

  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=0.0d0

END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,ThLoc,Rad
  REAL(RealKind) :: DeltaTh

  DeltaTh=DeltaThF(x,z)
  pLoc=p0*(One-kappa*Grav*z/(Rd*ThInit))**(cpD/Rd)
  thLoc=ThInit+DeltaTh
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE DensityBryan_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc

  pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThInit)
  RhoProf=0.0d0

END FUNCTION RhoProf


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=0.0d0
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=0.0d0
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=0.0d0
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE DensityBryan_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart





