MODULE Town_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Thermodynamic_Mod
  USE Physics_Mod

  IMPLICIT NONE 
  REAL(RealKind) :: N=0.0d-2
  REAL(RealKind), PARAMETER :: th0=293.16d0
  REAL(RealKind), PARAMETER :: D0=1.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind) :: uMax=10.0d0,vMax=0.0d0
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-4
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0  

  NAMELIST /Example/    &
                    uMax , &
                    vMax , &
                    TkeMax ,&
                    DisMax ,&
                    N



END MODULE Town_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Town_Mod
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
  USE Town_Mod
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

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,pLoc
  S=N*N/Grav
  ThLoc=th0*exp(z*S)
  IF (N>Zero) THEN
    pLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
  ELSE
    pLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
  END IF
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,ThLoc,pLoc
  S=N*N/Grav
  ThLoc=th0*exp(z*S)
  IF (N>Zero) THEN
    pLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
  ELSE
    pLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
  END IF
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S
  S=N*N/Grav
  ThProfFun=th0*exp(z*S)
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod 
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: RhoLoc,ThLoc
  RhoLoc=RhoFun(x,y,z,zHeight,Time)
  ThLoc=thProfFun(x,y,z,zHeight,Time)
  QvProfFun=qvr*qvs(RhoLoc,ThLoc,ThLoc)
END FUNCTION QvProfFun

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: U0
  UStart=uMax
END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Town_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod 
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=vMax
END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Town_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThStart=ThProfFun(x,y,z,zHeight,Time)
END FUNCTION ThStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=Zero
END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DHStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DHStart=Zero
END FUNCTION DHStart

FUNCTION DVStart(x,y,z,zHeight,Time)
  USE Town_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DVStart=Zero
END FUNCTION DVStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,RadLoc,ThBack
  S=N*N/Grav
  IF (N>Zero) THEN
    PreStart=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
  ELSE
    PreStart=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
  END IF
END FUNCTION PreStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=0.d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  HeightFun=Zero

END FUNCTION HeightFun


