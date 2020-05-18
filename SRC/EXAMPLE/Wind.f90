MODULE Wind_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE SoilData_Mod

  IMPLICIT NONE 

  
  REAL(RealKind) :: x01,y01,z01
  REAL(RealKind) :: rx1,ry1,rz1
  REAL(RealKind) :: r1,s1
  REAL(RealKind) :: DeltaTh1
  REAL(RealKind) :: x02,y02,z02
  REAL(RealKind) :: rx2,ry2,rz2
  REAL(RealKind) :: r2,s2
  REAL(RealKind) :: DeltaTh2
  REAL(RealKind) :: x0Qc,y0Qc,z0Qc
  REAL(RealKind) :: DeltaTh,uMax,vMax,ThInit,DMax
  REAL(RealKind) :: DeltaQc
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0
  REAL(RealKind) :: FacProf=0.0d0
  REAL(RealKind) :: zRauh=1.0d-1
  REAL(RealKind) :: uStar=0.5d0
  CHARACTER*8 :: BubbleType
  LOGICAL :: ProfIn
  
  NAMELIST /Example/ x01 &
                    ,y01 &
                    ,z01 &
                    ,rx1 &
                    ,ry1 &
                    ,rz1 &
                    ,r1  &
                    ,s1  &
                    ,DeltaTh1 &
                    ,x02 &

                    ,y02 &
                    ,z02 &
                    ,rx2 &
                    ,ry2 &
                    ,rz2 &
                    ,DeltaTh2 &
                    ,r2  &
                    ,s2  &
                    ,x0Qc &
                    ,y0Qc &
                    ,z0Qc &
                    ,uMax &
                    ,vMax &
                    ,DMax &
                    ,ThInit &
                    ,FacProf &
                    ,uStar &
                    ,zRauh &
                    ,ProfIn &
                    ,BubbleType

CONTAINS

FUNCTION ThPert(x,y,z)

  REAL(RealKind) :: ThPert
  REAL(RealKind) :: x,y,z
  REAL(RealKind) :: Rad
  SELECT CASE(BubbleType)
   CASE('Giraldo')
     Rad=SQRT(((x-x01)/rx1)**2 &
             +((z-z01)/rz1)**2 &
            ) 
     ThPert=Zero
     IF (Rad<One) THEN
       ThPert=0.5d0*DeltaTh1*(1.0d0+COS(Pi*Rad)) 
     END IF
   CASE('Cos2B')
     Rad=SQRT(((x-x01)/rx1)**2 &
             +((z-z01)/rz1)**2 &
            ) 
     ThPert=Zero
     IF (Rad<One) THEN
       ThPert=DeltaTh1*COS(Pi*Rad/2.0d0)**2 
     END IF
  
   CASE('Cos2BY')
     Rad=SQRT(((y-y01)/ry1)**2 &
             +((z-z01)/rz1)**2 &
            ) 
     ThPert=Zero
     IF (Rad<One) THEN
       ThPert=DeltaTh1*COS(Pi*Rad/2.0d0)**2 
     END IF
  
   CASE ('CB')
     Rad=SQRT(((x-x01)/rx1)**2 &
             +((z-z01)/rz1)**2 &
            ) 
     ThPert=Zero
     IF (Rad<One) THEN
       ThPert=DeltaTh1
     END IF
   CASE ('ExpB')
     Rad=SQRT((x-x01)**2+(z-z01)**2) 
     IF (Rad<r1) THEN
       ThPert=DeltaTh1
     ELSE
       ThPert=DeltaTh1*EXP(-((Rad-r1)/s1)**2)
     END IF
   CASE ('ExpB2')
     ThPert=Zero
     Rad=SQRT((x-x01)**2+(z-z01)**2) 
     IF (Rad<r1) THEN
       ThPert=ThPert+DeltaTh1
     ELSE
       ThPert=ThPert+DeltaTh1*EXP(-((Rad-r1)/s1)**2)
     END IF
     Rad=SQRT((x-x02)**2+(z-z02)**2) 
     IF (Rad<r2) THEN
       ThPert=ThPert+DeltaTh2
     ELSE
       ThPert=ThPert+DeltaTh2*EXP(-((Rad-r2)/s2)**2)
     END IF
   CASE DEFAULT
    STOP
  END SELECT

END FUNCTION ThPert


END MODULE Wind_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Wind_Mod
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
  USE Wind_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

! Find line
  x01=Zero
  x0Qc=Zero
  rx1=2000.0d0
  y01=Zero
  y0Qc=Zero
  ry1=2000.0d0
  z01=2000.0d0
  z0Qc=3000.0d0
  rz1=2000.0d0
  DMax=75.0d0
  uMax=0.0d0
  vMax=0.0d0
  ThInit=300.0d0
  DeltaTh=2.0d0
  DeltaQc=2.0d0
  BubbleType='Cos2B'
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


!FUNCTION RhoFun(x,y,z,zHeight,Time)
!
!  USE Parameter_Mod
!  USE Wind_Mod
!  USE ReadProfile_Mod
!  IMPLICIT NONE
!
!  REAL(RealKind) :: RhoFun
!
!  REAL(RealKind) :: x,y,z,zHeight,Time
!  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
!  CHARACTER*10, SAVE :: RhoType
!  LOGICAL, SAVE :: Load=.TRUE.
!  REAL(RealKind) :: S,pLoc,ThLoc,L


!  IF (ProfIn) THEN
!    IF (Load) THEN
!      CALL ReadProfile(cInt,RhoType,'RhoProf')
!      Load=.FALSE.
!    END IF
!    RhoFun=ProfileEqual(cInt,z)
!  ELSE
!    S=N*N/Grav
!    ThLoc=ThInit*exp(z*S)
!    IF (N>Zero) THEN
!      pLoc=p0*(One-Grav/(Cpd*ThInit*S)*(One-EXP(-S*z)))**(Cpd/Rd)
!    ELSE
!      pLoc=p0*(One-kappa*Grav*z/(Rd*ThInit))**(Cpd/Rd)
!    END IF
!    RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
!  END IF
!

FUNCTION UStart(x,y,z,zHeight,Time)

  USE Wind_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: UStartType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN !reads U(z) profile, file specified in grid file
    IF (Load) THEN
      CALL ReadProfile(cInt,UStartType,'UProf')
      Load=.FALSE.
    END IF
    UStart=ProfileEqual(cInt,z)
  ELSE
   uStart=uStar/Karm*LOG(z/zRauh+One)
  END IF  


END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  
  IF (ProfIn) THEN
    VStart=VMax
  ELSE
    VStart=0.0d0
!  VStart=uStar/Karm*LOG(z/zRauh+One)
  !VStart=vMax
  !VStart=VStart*Time
  END IF
END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  SELECT CASE (ThetaKind)
    CASE ('Density','PreEn','Exner','Energy')
      ThProfFun=0.0d0
    CASE ('Pseudo')
      ThProfFun=ThInit
    CASE ('Anelastic')  
      ThProfFun=ThInit
  END SELECT    
  
END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: pLoc,RhoLoc,ThLoc,TLoc

  SELECT CASE (ThetaKind)
    CASE ('Density','Pseudo','Anelastic')
      ThStart=ThInit+ThPert(x,y,z)
    CASE ('Energy')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      ThLoc=ThInit+ThPert(x,y,z)
      RhoLoc=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
      TLoc=pLoc/(Rd*RhoLoc)
      ThStart=Cvd*TLoc+Half*uMax*uMax+z*Grav
    CASE('PreEn')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      ThLoc=ThInit+ThPert(x,y,z)
      RhoLoc=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
      ThStart=pLoc/RhoLoc
    CASE('Exner')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      ThLoc=ThInit+ThPert(x,y,z)
      RhoLoc=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
      ThStart=(pLoc/p0)**kappa/RhoLoc
  END SELECT

END FUNCTION ThStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvStart=Zero

END FUNCTION QvStart

FUNCTION QrStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QrStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QrStart=Zero

END FUNCTION QrStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad

  Rad=SQRT(((x-x0Qc)/rx1)**2 &
          +((y-y0Qc)/ry1)**2 &
          +((z-z0Qc)/rz1)**2 &
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

FUNCTION DStart(x,y,z,zHeight,Time)

  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=DMax

END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,ThLoc

  SELECT CASE (ThetaKind)
    CASE ('Density','PreEn','Energy','Exner')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      ThLoc=ThInit+ThPert(x,y,z)
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
    CASE ('Anelastic')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      ThLoc=ThInit
      RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
    CASE ('Pseudo')
      RhoFun=Rho0
  END SELECT    

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Wind_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc

  SELECT CASE (ThetaKind)
    CASE ('Density','PreEn','Exner','Energy')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThInit)
      RhoProf=FacProf*RhoProf
    CASE ('Pseudo')
      RhoProf=Rho0
    CASE ('Anelastic')
      pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
      RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThInit)
  END SELECT    

END FUNCTION RhoProf


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Wind_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Wind_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=(RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time))/RhoFun(x,y,z,zHeight,Time)
END FUNCTION RhoStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreStart

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: pLoc,RhoLoc,ThLoc,TLoc
  pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
  ThLoc=ThInit+ThPert(x,y,z)
  RhoLoc=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  TLoc=pLoc/(Rd*RhoLoc)
  EnStart=Cvd*TLoc+Half*uMax*uMax+z*Grav
END FUNCTION EnStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Wind_Mod
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
  USE Wind_Mod
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
  USE Wind_Mod
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  IF (zSoil==0.0d0) THEN ! interception reservoir
    QvStartSoil=0.d0 
    QvStartSoil=MIN(QvStartSoil,5.d-4)
  ELSE
    QvStartSoil=0.0d0
    QvStartSoil=MIN(QvStartSoil,cporv(SoilType))
  END IF
  IF (SoilType<=2) QvStartSoil=0.0d0
END FUNCTION QvStartSoil

FUNCTION OmeStart(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  OmeStart=0.0d0
END FUNCTION OmeStart

FUNCTION DummyStart1(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start
FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE Wind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=Zero
END FUNCTION ForceTh
FUNCTION DampFun(z,Name)
  USE Wind_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
  DampFun=0.0d0
END FUNCTION DampFun
