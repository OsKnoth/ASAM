MODULE Town_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Thermodynamic_Mod
  USE Physics_Mod
  USE Random_Mod

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
  REAL(RealKind) :: Inflow_len=1.0d0  
  REAL(RealKind) :: offset_y=1.0d0  

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
  Inflow_len=domain%ny
  offset_y=domain%y0
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
  COMPLEX(REALKIND) :: UStart1,UStart2
  REAL(REALKIND) :: UStart3,UStart4,UStart5
  complex, parameter    :: i  = (0.0,1.0)
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: U0,w,k,knull,intenz
  INTEGER        :: ss,s

!  turb = ran1()

  UStart=0.0
  knull=5.0
  intenz=0.50 !Turbulenzintensitaet 0 - 1
  DO s = 1,2 !10 
     k = s**3 !7.0*(0.9)**s 
      DO ss = 1,2
     w = ((-1.0)**ss)*k*uMax/Inflow_len
!  UStart3=((k**4.0)*(exp(-2.0*(k/knull)**2.0)))/((knull**4.0)*(exp(-2.0*(knull/knull)**2.0))) &
!  UStart3=((((8.0/3.14)*(1.0**2.0)*(((k/knull)**4.0)/(knull*(1.0+(k/knull)**2.0)**3.0))/(4.0*3.14*k**2.0)))**(1.0/2.0)) &
!              *SIN(k*(y/Inflow_len)*2.0*3.1415+Time*(w)*2.0*3.1415 + k/10.3*2.0*3.1415)
     UStart3 =SIN(k*((y-offset_y)/Inflow_len)*2.0*3.1415 + Time*(w)*2.0*3.1415 + k*2.0*3.1415)
  UStart=UStart + UStart3
      ENDDO
  ENDDO
  UStart=(intenz*(SQRT((uMax**2.0)+(vMax**2.0)))*UStart/((s-1)*(ss-1))) + uMax
!  UStart=(((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(32.0/0.1d0))))*UStart
!   UStart=uMax! +0.01*(y-offset_y)
END FUNCTION UStart

!FUNCTION UStart(x,y,z,zHeight,Time)
!  USE Town_Mod
!  USE Rho_Mod 
!  IMPLICIT NONE
!  REAL(RealKind) :: UStart
!  REAL(RealKind) :: x,y,z,zHeight,Time
!  REAL(RealKind) :: yL,yR
!  REAL(RealKind) :: U0
!!  IF(z.LT.500) THEN
!!  turb = ran1()
!!  UStart=uMax-0.005d0+1.0d-2*turb
!  UStart=uMax
!!  ELSE
!!  USTART=uMax+Four
!!  ENDIF
!END FUNCTION UStart

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
  COMPLEX(REALKIND) :: VStart1,VStart2
  REAL(REALKIND) :: VStart3,VStart4,VStart5
  complex, parameter    :: i  = (0.0,1.0)
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: V0,w,k,knull, intenz
  INTEGER        :: ss,s

  VStart=0.0
  knull=5.0
  intenz=0.50 ! Turbulenzintensitaet 0 - 1 Anteil von uMax
  DO s = 1,2 !10 
     k = s**3 !7.0*(0.9)**s 
    DO ss = 1,2
     w = ((-1.0)**ss)*k*uMax/Inflow_len
!  VStart3=((k**4.0)*(exp(-2.0*(k/knull)**2.0)))/((knull**4.0)*(exp(-2.0*(knull/knull)**2.0))) &
!  VStart3=((((8.0/3.14)*(1.0**2.0)*(((k/knull)**4.0)/(knull*(1.0+(k/knull)**2.0)**3.0))/(4.0*3.14*k**2.0)))**(1.0/2.0)) &
!           *COS(k*(y/Inflow_len)*2.0*3.1415+Time*(w)*2.0*3.1415 + k/10.3*2.0*3.1415)
     VStart3 =COS(k*((y-offset_y)/Inflow_len)*2.0*3.1415  + Time*(w)*2.0*3.1415 + k*2.0*3.1415)
  VStart=VStart + ((VStart3))
     ENDDO
  ENDDO
  VStart=(intenz*(SQRT((uMax**2.0)+(vMax**2.0)))*VStart/((s-1)*(ss-1))) + vMax
!  VStart=(((0.44d0/0.4d0)*(log(z/0.1d0)))/((0.44d0/0.4d0)*(log(32.0/0.1d0))))*VStart
!  VStart = 0.0d0 !vMax
END FUNCTION VStart

!FUNCTION VStart(x,y,z,zHeight,Time)
!  USE Town_Mod
!  USE Rho_Mod 
!  IMPLICIT NONE
!  REAL(RealKind) :: VStart,turb
!  REAL(RealKind) :: x,y,z,zHeight,Time
!!  IF(z.LT.900) THEN
!!  turb = ran1()
!!  VStart=-0.5d0+ 1.0d0*turb
!!  ELSE
!  VSTART=0.0d0
!!  ENDIF
!END FUNCTION VStart

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
  REAL(RealKind) :: WStart,turb
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
!  IF(z.LT.900) THEN
!  turb = ran1()
!  WStart=-0.5d0+ 1.0d0*turb
!  ELSE
!  WSTART=0.0d0
!  ENDIF
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

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=Zero
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=Zero
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=Zero
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=Zero
END FUNCTION ForceRho

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1 = 0.0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2 = 0.0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE Town_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3 = 0.0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE Town_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil
  INTEGER :: LandClass,SoilType
  REAL(RealKind) :: Time
!  ThStartSoil=ThSoil
  IF (zSoil<=0.1d0) THEN
!    ThStartSoil=276.50d0
    ThStartSoil=290.50d0
  ELSE IF (zSoil<=0.5d0) THEN
!    ThStartSoil=277.70d0
    ThStartSoil=290.00d0
  ELSE IF (zSoil<=1.0d0) THEN
!    ThStartSoil=280.90d0
    ThStartSoil=289.50d0
  ELSE
!    ThStartSoil=286.50d0
    ThStartSoil=288.00d0
  END IF
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Town_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil
  INTEGER :: LandClass,SoilType
  REAL(RealKind) :: Time
  IF (zSoil<=0.1d0) THEN
    QvStartSoil=2.d-1
  ELSE IF (zSoil<=1.0d0) THEN
    QvStartSoil=1.8d-1
  ELSE
    QvStartSoil=1.2d-1
  END IF
END FUNCTION QvStartSoil

