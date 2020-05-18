MODULE GlobalOro_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: VelMax=20.0d0
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: ProfFac=0.0d0
  REAL(RealKind), PARAMETER :: Rho0=1.0d0
  REAL(RealKind), PARAMETER :: th0=293.16d0
  REAL(RealKind), PARAMETER :: D0=1.d-1


  NAMELIST /Example/    &
                    VelMax , &
                    ProfFac , &
                    N

END MODULE GlobalOro_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE GlobalOro_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

END SUBROUTINE SetBoundCells


SUBROUTINE InputExample(FileName)
  USE GlobalOro_Mod
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


FUNCTION RhoFun(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G,F,u0,pLoc
  u0=VelMax
  G=Grav*Grav/(N*N*cpD)
  S=N*N/Grav
  F=-u0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*u0)
  IF (N>Zero) THEN
    ExPre=(One-G/th0)+G/th0*EXP(-S*(z+F/Grav))
  ELSE
    ExPre=One-Grav/(th0*cpD)*(z+F/Grav)
  END IF
  pLoc=p0*ExPre**(1.0d0/kappa)
  ThLoc=th0*exp(S*(z+F/Grav))
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
END FUNCTION RhoFun

FUNCTION RhoProf(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G,F,u0,pLoc
  u0=VelMax
  G=Grav*Grav/(N*N*cpD)
  S=N*N/Grav
  F=-u0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*u0)
  IF (N>Zero) THEN
    ExPre=(One-G/th0)+G/th0*EXP(-S*(z+F/Grav))
  ELSE
    ExPre=One-Grav/(th0*cpD)*(z+F/Grav)
  END IF
  pLoc=p0*ExPre**(1.0d0/kappa)
  ThLoc=th0*exp(S*(z+F/Grav))
  RhoProf=ProfFac*pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoProf

FUNCTION PreStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G,F,u0
  u0=VelMax
  G=Grav*Grav/(N*N*cpD)
  S=N*N/Grav
  F=-u0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*u0)
  IF (N>Zero) THEN
    ExPre=(One-G/th0)+G/th0*EXP(-S*(z+F/Grav))
  ELSE
    ExPre=One-Grav/(th0*cpD)*(z+F/Grav)
  END IF
  PreStart=p0*ExPre**(1.0d0/kappa)
END FUNCTION PreStart

FUNCTION PreProf(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G
  S=N*N/Grav
  G=Grav*Grav/(N*N*cpD)
  IF (N>Zero) THEN
    ExPre=(One-G/th0)+G/th0*EXP(-S*z)
  ELSE
    ExPre=One-Grav/(th0*cpD)*z
  END IF
  PreProf=p0*ExPre**(1.0d0/kappa)
END FUNCTION PreProf
FUNCTION ThProfFun(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,F,u0
  u0=VelMax
  F=-u0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*u0)
  S=N*N/Grav
  ThProfFun=ProfFac*th0*exp(S*(z+F/Grav))
END FUNCTION ThProfFun

FUNCTION QvProfFun(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  USE ThProf_Mod
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QvProfFun=0.0d0
END FUNCTION QvProfFun

FUNCTION UStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: UStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: phiL,phiR
  REAL(RealKind) :: U0
  U0=VelMax
  UStart=U0*COS(phi)
END FUNCTION UStart

FUNCTION UStartE(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  USE UVW_Mod
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  UStartE=UStart(lam,phi,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: VStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  VStart=0.0d0
END FUNCTION VStart

FUNCTION VStartE(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  VStartE=VStart(lam,phi,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: WStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION ThStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  USE ThProf_Mod
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,F,u0
  u0=VelMax
  F=-u0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*u0)
  S=N*N/Grav
  ThStart=th0*exp(S*(z+F/Grav))
END FUNCTION ThStart

FUNCTION TkeStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  TkeStart=0.0d0
END FUNCTION TkeStart

FUNCTION DisStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DisStart=0.0d0
END FUNCTION DisStart

FUNCTION QvStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QvStart=0.d0
END FUNCTION QvStart

FUNCTION QcStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: DStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DummyStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION QiStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

!FUNCTION RhoStart(lam,phi,z,zHeight,Time)
!  USE GlobalOro_Mod
!  USE Rho_Mod
!  USE RhoProf_Mod
!  IMPLICIT NONE
!  REAL(RealKind) :: RhoStart
!  REAL(RealKind) :: lam,phi,z,zHeight,Time
!  RhoStart=RhoFun(lam,phi,z,zHeight,Time)  &
!          -RhoProf(lam,phi,z,zHeight,Time)
!END FUNCTION RhoStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =Zero
END FUNCTION TStart


FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=Zero
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=Zero
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=Zero
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=Zero
END FUNCTION ForceRho


FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=Zero
END FUNCTION HeightFun






