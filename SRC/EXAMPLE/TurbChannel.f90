MODULE TurbChannel_Mod

  USE Kind_Mod
  USE Thermodynamic_Mod !added
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE Random_Mod ! added

  IMPLICIT NONE 

  REAL(RealKind) :: uMax=0.0d0
  REAL(RealKind) :: vMax=0.0d0
  REAL(RealKind) :: wMax=1.0d0
  REAL(RealKind) :: thL=290.0d0
  REAL(RealKind) :: thR=300.0d0
  REAL(RealKind) :: uStar=0.4d0
  REAL(RealKind) :: alpha=0.00d0
  REAL(RealKind) :: N=0.0d-2

! REAL(RealKind), PARAMETER :: th0=289.2d0 !293.16d0
  REAL(RealKind), PARAMETER :: D0=1.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
! REAL(RealKind) :: uMax=10.0d0,vMax=0.0d0
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-4
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0
  REAL(RealKind) :: Inflow_lenx=1.0d0
  REAL(RealKind) :: Inflow_leny=1.0d0
  REAL(RealKind) :: Inflow_lenz=1.0d0
  REAL(RealKind) :: offset_x=1.0d0
  REAL(RealKind) :: offset_x1=1.0d0
  REAL(RealKind) :: offset_y=1.0d0
  REAL(RealKind) :: offset_y1=1.0d0
  REAL(RealKind) :: offset_z=1.0d0
  REAL(RealKind) :: offset_z1=1.0d0
  REAL(RealKind) :: intenz=0.1d0

  NAMELIST /Example/    &
           uMax , &
           vMax , &
           wMax , &
           uStar , &
           thL , &
           thR , &
           alpha, &
           TkeMax ,&
           DisMax ,&
           intenz ,&
           N
                                                                
CONTAINS
FUNCTION uInFlow(x,y,z,Time)
  REAL(RealKind) :: uInFlow
  REAL(RealKind) :: x,y,z,Time

  INTEGER        :: k_end
  COMPLEX(REALKIND) :: UStart1,UStart2
  REAL(REALKIND) :: UStart3,UStart4,UStart5
  complex, parameter    :: i  = (0.0,1.0)
  REAL(RealKind) :: yL,yR,r
  REAL(RealKind) :: U0,w,wx,wy,wz,k,knull,AMP=1.0
  INTEGER        :: ss,s
!WRITE(*,*) 'x',offset_x,offset_y,inflow_lenx,inflow_leny
!  //////turb = ran1()
  IF (Time>1.D-1) THEN
  CALL Random_Number(r)
  uInFlow=0.0d0
  knull=18.0
  k_end = int(Inflow_leny/4.0)
  DO s = 1,10 !k_end
     k = ((10.0)*(0.9**s)) !*(1.0 + ((- 0.5d-1 + 1.0d-1*r)/((6.0*(0.8**1.0)))))
!40.0*(0.97)**s
!     k = 5.0
      DO ss = 1,2
!     w = ((-1.0)**ss)*k*(SQRT(uMax))/Inflow_leny 
     wx = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_x1-offset_x)
     wy = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_y1-offset_y)
     wz = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_z1-offset_z)
!IF(y.GT.offset_y-1.and.x.gt.1) THEN
!UStart2 = 1.0*AMP*SIN(k*((x-offset_x)/Inflow_lenx)*2.0*3.1415 + Time*(w)*2.0*3.1415) !- 0.5d-1 + 1.0d-1*r  
!UStart = UStart + UStart2
!ELSE
     UStart3 = AMP*SIN(k*((x-offset_x)/(offset_x1-offset_x))*2.0*3.1415 + &
               Time*(wx)*2.0*3.1415 + k**(2.0))  &!- 0.5d-1 + 1.0d-1*r  
             + AMP*SIN(k*((y-offset_y)/(offset_y1-offset_y))*2.0*3.1415 + &
               Time*(wy)*2.0*3.1415 + k**(2.0))  &
             + AMP*SIN(k*((z-offset_z)/(offset_z1-offset_z))*2.0*3.1415 + &
               Time*(wz)*2.0*3.1415 + k**(2.0))
  uInFlow=uInFlow + UStart3
!ENDIF
      ENDDO
  ENDDO
    uInFlow=(intenz*(SQRT((uMax**2.0)+(vMax**2.0)))*uInFlow/((s-1)*(ss-1))) + uMax
  ELSE
    uInFlow=UMax
  END IF  
 
!  WRITE(*,*)'offset_z',offset_z,'offset_z1',offset_z1,'Inflow_lenz',Inflow_lenz
!  IF(Time.LT.1.0d-4) UStart = uMax
!  IF(x.GT.(offset_x)) UStart = uMax
!  UStart=
!  (((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(16.0/0.1d0))))*UStart -
!  0.5d-1 + 1.0d-1*r
!  0.5d-1 + 1.0d-1*r
!  UStart=
!  (((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(16.0/0.1d0))))*uMax -
!  0.5d-1 + 1.0d-1*r
!   UStart=uMax +0.01*(y-offset_y)
!   UStart = uMax !- 0.5d-1 +1.0d-1*r
!   IF(y.lt.0.51) UStart =uMax -0.5d-1 +1.0d-1*r
!  UStart=UStart - 0.5d0 + 1.0d0*r
!!   IF(z.lt.(3.0/186.0)) UStart = 0.0d0 -0.5d-1 +1.0d-1*r


END FUNCTION uInFlow



END MODULE TurbChannel_Mod


SUBROUTINE SetBoundCells(BoundCellLoc)

  USE TurbChannel_Mod
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
  USE TurbChannel_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

  Inflow_lenx=domain%nx ! added
  Inflow_leny=domain%ny
  Inflow_lenz=domain%nz
  offset_x=domain%x0
  offset_x1=domain%x1
  offset_y=domain%y0
  offset_y1=domain%y1
  offset_z=domain%z0
  offset_z1=domain%z1

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
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  IF (y<0.0d0) THEN
    ThLoc=thL
  ELSE
    ThLoc=thR
  END IF
  pLoc=p0*(One-kappa*Grav*z/(Rd*thLoc))**(Cpd/Rd)
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  RhoProf=0.0d0

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc,RhoLoc
  ThProfFun=0.0d0
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvProfFun=Zero
END FUNCTION QvProfFun


  


FUNCTION UStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  UStart=0.0d0

END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE



FUNCTION VStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart,UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=0.0d0

END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  IF (y<0.0d0) THEN 
    WStart=-wMax
  ELSE
    WStart=-1.0d0*wMax
  END IF  
END FUNCTION WStart

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,ThLoc,RhoLoc

  IF (y<0.0d0) THEN
    ThStart=thL
  ELSE  
    ThStart=thR
  END IF
END FUNCTION ThStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=(RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time))/RhoFun(x,y,z,zHeight,Time)
END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=0.0d0
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,RadLoc,ThBack
  PreStart=0.0d0
END FUNCTION PreStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=0.0d0
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=0.0d0
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=0.0d0
END FUNCTION LenStart


FUNCTION TStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc,Fpot,RhoLoc
  TStart =0.0d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
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
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  USE UVW_Mod
  USE Rho_Mod
  USE Start_Mod, ONLY: QvStart,QcStart
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: u,v,w,Rho,RhoV,RhoC,p

  EnStart=0.0d0
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
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE TurbChannel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=Zero
END FUNCTION ForceTh
FUNCTION DampFun(z,Name)
  USE TurbChannel_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
  DampFun=0.0d0
END FUNCTION DampFun
