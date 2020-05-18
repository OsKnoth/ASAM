MODULE GravityWaveGenerator_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: k1,k2,l1,l2
  REAL(RealKind) :: Lk1=10000.0d0
  REAL(RealKind) :: Lk2=11000.0d0
  REAL(RealKind) :: Ll1=2500.0d0
  REAL(RealKind) :: Ll2=1500.0d0
  REAL(RealKind) :: OmegaG=0.002d0
  REAL(RealKind) :: alpha=0.2d0
  REAL(RealKind) :: zShift=-5000.0d0
  REAL(RealKind) :: H=10000.0d0
  REAL(RealKind) :: UBot=5.0d0
  REAL(RealKind) :: Utop=10.0d0
  REAL(RealKind) :: th0=293.16d0
  REAL(RealKind) :: D0=10.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind) :: uMax=10.0d0,vMax=0.0d0
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-4
  REAL(RealKind) :: TkeHMax=1.0d-2,TkeVMax=1.0d-2,LenMax=1.0d3

  NAMELIST /Example/    &
                    uBot , &
                    uTop , &
                    zShift, &
                    N,    &
                    th0

CONTAINS

FUNCTION E(x,z)

  REAL(RealKind) :: E
  REAL(RealKind) :: x,z

  IF (ABS(x)<=Pi/k2.AND.ABS(z-zShift)<=Pi/l2) THEN
    E=alpha*(One+COS(k2*x))*(One+COS(l2*(z-zShift)))
  ELSE
    E=Zero
  END IF

END FUNCTION E

FUNCTION DxE(x,z)

  REAL(RealKind) :: DxE
  REAL(RealKind) :: x,z

  IF (ABS(x)<=Pi/k2.AND.ABS(z-zShift)<=Pi/l2) THEN
    DxE=-alpha*k2*SIN(k2*x)*(One+COS(l2*(z-zShift)))
  ELSE
    DxE=Zero
  END IF

END FUNCTION DxE

FUNCTION DzE(x,z)

  REAL(RealKind) :: DzE
  REAL(RealKind) :: x,z

  IF (ABS(x)<=Pi/k2.AND.ABS(z-zShift)<=Pi/l2) THEN
    DzE=-alpha*l2*(One+COS(k2*x))*SIN(l2*(z-zShift))
  ELSE
    DzE=Zero
  END IF

END FUNCTION DzE

FUNCTION DxPsi(x,z,t)

  REAL(RealKind) :: DxPsi
  REAL(RealKind) :: x,z,t

  DxPsi=DxE(x,z)*SIN(OmegaG*t)*SIN(k1*x)*COS(l1*z-zShift) &
       +E(x,z)*SIN(OmegaG*t)*k1*COS(k1*x)*COS(l1*z-zShift) 

END FUNCTION DxPsi

FUNCTION DzPsi(x,z,t)

  REAL(RealKind) :: DzPsi
  REAL(RealKind) :: x,z,t

  DzPsi=DzE(x,z)*SIN(OmegaG*t)*SIN(k1*x)*COS(l1*z-zShift) &
       -E(x,z)*SIN(OmegaG*t)*SIN(k1*x)*l1*SIN(l1*z-zShift)

END FUNCTION DzPsi

FUNCTION ThetaLoc(z)

  REAL(RealKind) :: ThetaLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: S

  S=N*N*z/Grav
  ThetaLoc=th0*EXP(S)

END FUNCTION ThetaLoc

FUNCTION PressLoc(z)

  REAL(RealKind) :: PressLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: S

  S=N*N/Grav
  PressLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
           
END FUNCTION PressLoc

END MODULE GravityWaveGenerator_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE GravityWaveGenerator_Mod
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
  USE GravityWaveGenerator_Mod
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
  k1=Two*Pi/Lk1
  k2=Two*Pi/Lk2
  l1=Two*Pi/Ll1
  l2=Two*Pi/Ll2
END SUBROUTINE InputExample

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time

  ForceU=-RhoFun(x,y,z,zHeight,Time)*DzPsi(x,z,Time)

END FUNCTION ForceU


FUNCTION ForceV(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time

  ForceV=Zero
  
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time

  ForceW=RhoFun(x,y,z,zHeight,Time)*DxPsi(x,z,Time)

END FUNCTION ForceW

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  pLoc=PressLoc(z)
  ThLoc=ThetaLoc(z)
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  pLoc=PressLoc(z)
  ThLoc=ThetaLoc(z)
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  RhoProf=Zero

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc
  ThLoc=ThetaLoc(z)
  ThProfFun=Zero
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION UStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: U0
  UStart=(z*UTop+(H-z)*UBot)/H
END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=vMax
END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThStart=ThetaLoc(z)
END FUNCTION ThStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,RadLoc,ThBack
  PreStart=PressLoc(z)
END FUNCTION PreStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart


FUNCTION TStart(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc,Fpot,RhoLoc
  ThLoc  =ThetaLoc(z)
  pLoc   =PressLoc(z)
  RhoLoc =pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  Fpot   = (Rd*RhoLoc*ThLoc/p0)**(-kappa/(One-kappa))
  TStart =ThLoc/Fpot
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE GravityWaveGenerator_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=Zero
END FUNCTION HeightFun



