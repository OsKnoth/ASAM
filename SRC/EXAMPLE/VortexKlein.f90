MODULE VortexKlein_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: r0=0.125d4 !Core radius of vortex
  REAL(RealKind) :: UIn=-1024.0d0 !Core velocity
  REAL(RealKind) :: ThInit=300.0d0 !Core velocity
  REAL(RealKind) :: RhoSurf=0.5d0
  REAL(RealKind) :: pSurf=1.d3
  REAL(RealKind) :: UTrans=0.d0
  
  NAMELIST /Example/ r0 &
                    ,RhoSurf &
                    ,pSurf &
                    ,UTrans &
                    ,UIn 


CONTAINS

FUNCTION VRad(r)
  REAL(RealKind) :: VRad
  REAL(RealKind) :: r

  REAL(RealKind) :: rLoc

  rLoc=r/r0
  IF (rLoc<=1.0d0) THEN
    VRad=UIn*(1.d0-rLoc)**6.d0*rLoc**6.d0
  ELSE
    VRad=0.0d0
  END IF  
END FUNCTION VRad

FUNCTION RhoRad(r)
  REAL(RealKind) :: RhoRad
  REAL(RealKind) :: r

  REAL(RealKind) :: rLoc

  rLoc=r/r0
  IF (rLoc<=1.0d0) THEN
    RhoRad=RhoSurf*(1.0d0+(1.d0-rLoc**2.d0)**6.d0)
  ELSE
    RhoRad=RhoSurf
  END IF  
END FUNCTION RhoRad

FUNCTION PresRad(r)
  REAL(RealKind) :: PresRad
  REAL(RealKind) :: r

  REAL(RealKind) :: rLoc

  rLoc=MIN(r/r0,1.0d0)
!   PresRad=p0+ABS(UIn)*27.0d0/4.0d0*(1.0d0/30.0d0+(1.0d0-rLoc)**5*(-1.0d0/30.0d0-rLoc/6.0d0))
    PresRad=pSurf+UIn**2.0d0*(rLoc**(36.0d0)/72.0d0 - (6.0d0*rLoc**(35.0d0))/35.0d0            &
                      + (15.0d0*rLoc**(34.0d0))/17.0d0 - (74.0d0*rLoc**(33.0d0))/33.0d0     &
                      + (57.0d0*rLoc**(32.0d0))/32.0d0 + (174.0d0*rLoc**(31.0d0))/31.0d0    &
                      - (269.0d0*rLoc**(30.0d0))/15.0d0 + (450.0d0*rLoc**(29.0d0))/29.0d0   &
                      + (153.0d0*rLoc**(28.0d0))/8.0d0 - (1564.0d0*rLoc**(27.0d0))/27.0d0   &
                      + (510.0d0*rLoc**(26.0d0))/13.0d0 + (204.0d0*rLoc**(25.0d0))/5.0d0    &
                      - (1473.0d0*rLoc**(24.0d0))/16.0d0 + (1014.0d0*rLoc**(23.0d0))/23.0d0 &
                      + (1053.0d0*rLoc**(22.0d0))/22.0d0 - (558.0d0*rLoc**(21.0d0))/7.0d0   &
                      + (783.0d0*rLoc**(20.0d0))/20.0d0 + (54.0d0*rLoc**(19.0d0))/19.0d0    &
                      - (38.0d0*rLoc**(18.0d0))/9.0d0 - (222.0d0*rLoc**(17.0d0))/17.0d0     &
                      + (609.0d0*rLoc**(16.0d0))/32.0d0 - (184.0d0*rLoc**(15.0d0))/15.0d0   &
                      + (9.0d0*rLoc**(14.0d0))/2.0d0 - (12.0d0*rLoc**(13.0d0))/13.0d0       &
                      + rLoc**(12.d0)/12.d0 - 34373.d0/1805044411170.d0)

!r^36/36 - (12*r^35)/35 + (30*r^34)/17 - (148*r^33)/33 + (57*r^32)/16 
!+ (348*r^31)/31 - (538*r^30)/15 + (900*r^29)/29 + (153*r^28)/4 - (3128*r^27)/27 
!+ (1020*r^26)/13 + (408*r^25)/5 - (1473*r^24)/8 + (2028*r^23)/23 
!+ (1053*r^22)/11 - (1116*r^21)/7 + (783*r^20)/10 + (108*r^19)/19 
!- (76*r^18)/9 - (444*r^17)/17 + (609*r^16)/16 
!- (368*r^15)/15 + 9*r^14 - (24*r^13)/13 + r^12/6

!r^36/72 - (6*r^35)/35 + (15*r^34)/17 - (74*r^33)/33 + (57*r^32)/32 + (174*r^31)/31 
!- (269*r^30)/15 + (450*r^29)/29 + (153*r^28)/8 - (1564*r^27)/27 + (510*r^26)/13 
!+ (204*r^25)/5 - (1473*r^24)/16 + (1014*r^23)/23 + (1053*r^22)/22 
!- (558*r^21)/7 + (783*r^20)/20 + (54*r^19)/19 - (38*r^18)/9 - (222*r^17)/17 
!+ (609*r^16)/32 - (184*r^15)/15 + (9*r^14)/2 - (12*r^13)/13 + r^12/12
 
END FUNCTION PresRad

END MODULE VortexKlein_Mod


SUBROUTINE SetBoundCells(BoundCellLoc)

  USE VortexKlein_Mod
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
  USE VortexKlein_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

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
END SUBROUTINE InputExample

FUNCTION UStart(x,y,z,zHeight,Time)

  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: r,U
  r=SQRT((x*x)+(y*y))
  U=VRad(r)
  UStart=-U*y/r+UTrans

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: r,U
  r=SQRT((x*x)+(y*y))
  U=VRad(r)
  VStart=U*x/r

END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  ThProfFun=ThInit
END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: r,pLoc,RhoLoc

! INTEGER :: i
! REAL(RealKind) :: p1,p2,V,DpDr,DpDrA
!
! r=1.0d0
! DO i=1,100
!   p1=PresRad(r)
!   p2=PresRad((1.0d0+1.d-8)*r)
!   DpDrA=(p2-p1)/((1.0d0+1.d-8)*r-r)
!   V=VRad(r)
!   RhoLoc=RhoRad(r)
!   DpDr=V*V/r*RhoLoc
!   WRITE(*,*) r,DpDrA,DpDr,DpDrA/DpDr
!   r=r+20.0d0
! END DO  
! STOP

  r=SQRT((x*x)+(y*y))
  pLoc=PresRad(r)
  RhoLoc=RhoRad(r)
  ThStart=pLoc/((pLoc/p0)**kappa*RhoLoc*Rd)
END FUNCTION ThStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvStart=Zero

END FUNCTION QvStart

FUNCTION DStart(x,y,z,zHeight,Time)

  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=Zero

END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: r

  r=SQRT((x*x)+(y*y))
  RhoFun=RhoRad(r)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE VortexKlein_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc

  RhoProf=Zero

END FUNCTION RhoProf


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart


FUNCTION PreStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE VortexKlein_Mod
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
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil,ThAir
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  ThStartSoil=0.0d0
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  QvStartSoil=0.0d0
END FUNCTION QvStartSoil

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  EnStart=Zero
END FUNCTION EnStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.0d0
END FUNCTION QiStart

FUNCTION QsStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QsStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QsStart=0.0d0
END FUNCTION QsStart

FUNCTION NvStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NvStart=0.0d0
END FUNCTION NvStart

FUNCTION NcStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NcStart=0.0d0
END FUNCTION NcStart

FUNCTION NrStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NrStart=0.0d0
END FUNCTION NrStart

FUNCTION NiStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NiStart=0.0d0
END FUNCTION NiStart

FUNCTION NsStart(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: NsStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  NsStart=0.0d0
END FUNCTION NsStart

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=0.0d0
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=0.0d0
END FUNCTION Tracer2Start

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=0.0d0
END FUNCTION HeightFun

FUNCTION OmeStart(lam,phi,z,zHeight,Time)
  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  OmeStart = 0.0
END FUNCTION OmeStart
