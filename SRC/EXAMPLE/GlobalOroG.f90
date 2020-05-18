MODULE GlobalOro_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: VelMax=20.0d0
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind), PARAMETER :: Rho0=1.0d0
  REAL(RealKind), PARAMETER :: th0=293.16d0
  REAL(RealKind), PARAMETER :: D0=1.d-1


  NAMELIST /Example/    &
                    VelMax , &
                    N

END MODULE GlobalOro_Mod

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
  S=N*N/Grav
  F=-u0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*u0)
  IF (N>Zero) THEN
    G=Grav*Grav/(N*N*cp)
    ExPre=(One-G/th0)+G/th0*EXP(-S*(z+F/Grav))
  ELSE
    ExPre=One-Grav/(th0*cp)*(z+F/Grav)
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
  REAL(RealKind) :: S,ThLoc,ExPre,G,pLoc
  S=N*N/Grav
  IF (N>Zero) THEN
    G=Grav*Grav/(N*N*cp)
    ExPre=(One-G/th0)+G/th0*EXP(-S*z)
  ELSE
    ExPre=One-Grav/(th0*cp)*z
  END IF
  pLoc=p0*ExPre**(1.0d0/kappa)
  ThLoc=th0*exp(S*z)
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoProf

FUNCTION PreStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,ExPre,G,F,u0
  u0=VelMax
  S=N*N/Grav
  F=-u0*COS(phi)**2*(RadEarth*Omega+Curv*0.5d0*u0)
  IF (N>Zero) THEN
    G=Grav*Grav/(N*N*cp)
    ExPre=(One-G/th0)+G/th0*EXP(-S*(z+F/Grav))
  ELSE
    ExPre=One-Grav/(th0*cp)*(z+F/Grav)
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
  IF (N>Zero) THEN
    G=Grav*Grav/(N*N*cp)
    ExPre=(One-G/th0)+G/th0*EXP(-S*z)
  ELSE
    ExPre=One-Grav/(th0*cp)*z
  END IF
  PreProf=p0*ExPre**(1.0d0/kappa)
END FUNCTION PreProf
FUNCTION ThProfFun(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  REAL(RealKind) :: S
  S=N*N/Grav
  ThProfFun=th0*EXP(S*z)
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

FUNCTION RhoStart(lam,phi,z,zHeight,Time)
  USE GlobalOro_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  RhoStart=RhoFun(lam,phi,z,zHeight,Time)  &
          -RhoProf(lam,phi,z,zHeight,Time)
END FUNCTION RhoStart


