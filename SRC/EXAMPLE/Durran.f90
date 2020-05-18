MODULE Durran_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 
  REAL(RealKind) :: NL=1.0d-2
  REAL(RealKind) :: NU=1.0d-2
  REAL(RealKind) :: H=1.0d-2
  REAL(RealKind) :: RhoL=1.15d0
  REAL(RealKind) :: RhoU=0.0d0
  REAL(RealKind) :: th0=293.16d0
  REAL(RealKind) :: D0=10.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind) :: uMax=10.0d0,vMax=0.0d0
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-4
  REAL(RealKind) :: TkeHMax=1.0d-2,TkeVMax=1.0d-2,LenMax=1.0d3
  REAL(RealKind) :: xL1,xL2,xL3,xL4

  NAMELIST /Example/    &
                    uMax , &
                    vMax , &
                    NL,    &
                    NU,    &
                    th0,   &
                    D0,    &
                    H,     &
                    xL1,   &
                    xL2,   &
                    xL3,   &
                    xL4

CONTAINS


FUNCTION ThetaLoc(z)

  REAL(RealKind) :: ThetaLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: SL,SU

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    ThetaLoc=th0*EXP(SL)
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
    ThetaLoc=th0*EXP(SL+SU)
  END IF

END FUNCTION ThetaLoc

FUNCTION PressLoc(z)

  REAL(RealKind) :: PressLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: SL,SU
  REAL(RealKind) :: t1,t2,pH
  REAL(RealKind) :: s1,s2

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    t1=Grav*Grav*kappa/(Rd*th0*NL*NL)
    t2=Grav*Grav*kappa/(Rd*th0*NL*NL)*EXP(-SL)
    PressLoc=p0 &
            *(1.0d0+t2-t1)**(1.0d0/kappa)
           
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
    t1=Grav*Grav*kappa/(Rd*th0*NL*NL)
    t2=Grav*Grav*kappa/(Rd*th0*NL*NL)*EXP(-SL)
    pH=p0 &
       *(1.0d0+t2-t1)**(1.0d0/kappa)
    s1=Grav*Grav*kappa/(Rd*th0*NU*NU)*EXP(-SL)
    s2=Grav*Grav*kappa/(Rd*th0*NU*NU)*EXP(-SL-SU)
    PressLoc=p0 &
       *(1.0d0+t2-t1+s2-s1)**(1.0d0/kappa)
  END IF

END FUNCTION PressLoc

FUNCTION TotalEnergy(u,v,w,p,Rho,RhoV,RhoL,z)
  REAL(RealKind) :: TotalEnergy
  REAL(RealKind) :: u,v,w,p,Rho,RhoV,RhoL,z

  REAL(RealKind) :: RhoDLoc,TLoc,Cvml

  RhoDLoc=Rho-RhoV-RhoL
  Cvml=Cvd*RhoDLoc+Cvv*RhoV+Cpl*RhoL
  TLoc=p/(RhoV*Rv+RhoDLoc*Rd)
  TotalEnergy=(Cvml*TLoc & 
             +Half*Rho*(u*u+v*v+w*w) &
             +RhoV*L00 &
             +Rho*z*Grav)/(Rho+Eps)
END FUNCTION TotalEnergy

END MODULE Durran_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Durran_Mod
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
  USE Durran_Mod
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
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  pLoc=PressLoc(z)
  ThLoc=ThetaLoc(z)
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  RhoProf=0.0d0

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc,RhoLoc
  ThProfFun=0.0d0
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: U0
  UStart=uMax
!  UStart=0.0d0
!  UStart=5.0
!  UStart=5*sin(x/5000*Pi)*cos(z/3000*Pi)
END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=vMax
END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
!  WStart=5*cos(x/5000*Pi+Pi)*sin(z/3000*Pi)
  WStart=0.0d0
END FUNCTION WStart

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,ThLoc,RhoLoc
  SELECT CASE (ThetaKind)
    CASE ('Density')
      ThStart=ThetaLoc(z)
    CASE('PreEn')
      pLoc=PressLoc(z)
      ThLoc=ThetaLoc(z)
      RhoLoc=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
      ThStart=pLoc/RhoLoc
  END SELECT
END FUNCTION ThStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=(RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time))/RhoFun(x,y,z,zHeight,Time)
END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,RadLoc,ThBack
  PreStart=PressLoc(z)
END FUNCTION PreStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart


FUNCTION TStart(x,y,z,zHeight,Time)
  USE Durran_Mod
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
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=Zero
END FUNCTION HeightFun

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
  INTEGER :: LandClass, SoilType
  ThStartSoil=0.0d0
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Kind_Mod
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass, SoilType
  QvStartSoil=0.0d0
END FUNCTION QvStartSoil


FUNCTION DummyStart1(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE UVW_Mod
  USE Rho_Mod
  USE Start_Mod, ONLY: QvStart,QcStart
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: u,v,w,Rho,RhoV,RhoC,p

  u=UStart(x,y,z,zHeight,Time)
  v=VStart(x,y,z,zHeight,Time)
  w=WStart(x,y,z,zHeight,Time)
  Rho=RhoFun(x,y,z,zHeight,Time)
  RhoV=QvStart(x,y,z,zHeight,Time)*Rho
  RhoC=QcStart(x,y,z,zHeight,Time)*Rho
  p=PressLoc(z)
  EnStart=TotalEnergy(u,v,w,p,Rho,RhoV,RhoC,z)
END FUNCTION EnStart


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
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  OmeStart=0.0d0
END FUNCTION OmeStart

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=Zero
END FUNCTION ForceTh
FUNCTION DampFun(z,Name)
  USE Durran_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
  DampFun=0.0d0
END FUNCTION DampFun
