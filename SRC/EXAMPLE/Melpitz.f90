MODULE Melpitz_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE SoilData_Mod

  IMPLICIT NONE 

! Wind
  REAL(RealKind) :: uMax=5.0d0
  REAL(RealKind) :: vMax=5.0d0
! Soil temperature
  REAL(RealKind) :: TS05=292.15d0
  REAL(RealKind) :: TS15=296.05d0
  REAL(RealKind) :: TS25=296.75d0
  REAL(RealKind) :: TS35=295.85d0
  REAL(RealKind) :: TS45=294.95d0
  REAL(RealKind) :: TS55=293.85d0
  REAL(RealKind) :: TS80=290.35d0
  REAL(RealKind) :: TS150=288.95d0
! Soil moisture (% plant-available water)
  REAL(RealKind) :: BF10=9.0d0
  REAL(RealKind) :: BF20=8.0d0
  REAL(RealKind) :: BF30=11.0d0
  REAL(RealKind) :: BF40=37.5d0
  REAL(RealKind) :: BF50=60.0d0
  REAL(RealKind) :: BF60=70.0d0
! Perturbations
  LOGICAL        :: Perturb=.FALSE.
  REAL(RealKind) :: DeltaT=0.02d0
  REAL(RealKind) :: HeightDeltaT=1000.0d0
  REAL(RealKind) :: DeltaQv=0.05d-3
  REAL(RealKind) :: HeightDeltaQv=1000.0d0
! Others  
  REAL(RealKind) :: ThInit=285.15d0 ! 13 degC
  REAL(RealKind) :: N=0.01d0
  LOGICAL :: ProfIn=.FALSE.
  LOGICAL :: ProfInVel=.FALSE.

  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0

  NAMELIST /Example/ uMax &
                    ,vMax &
                    ,TS05 &
                    ,TS15 &
                    ,TS25 &
                    ,TS35 &
                    ,TS45 &
                    ,TS55 &
                    ,TS80 &
                    ,TS150 &
                    ,BF10 &
                    ,BF20 &
                    ,BF30 &
                    ,BF40 &
                    ,BF50 &
                    ,BF60 &
                    ,ThInit &
                    ,N &
                    ,Perturb &
                    ,DeltaT &
                    ,HeightDeltaQv &
                    ,DeltaQv &
                    ,HeightDeltaT &
                    ,ProfIn & 
                    ,ProfInVel  

END MODULE Melpitz_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)
  USE Melpitz_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc
END SUBROUTINE SetBoundCells

SUBROUTINE PerturbProfile(VecC)

  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Parameter_Mod
  USE Floor_Mod
  USE Melpitz_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: i,ix,iy,iz
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: ThVirt(:,:,:,:)
  REAL(RealKind), POINTER :: En(:,:,:,:)
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDry
  REAL(RealKind) :: KappaLoc,pLoc,rvLoc,rlLoc,rtLoc,qvLoc,qlLoc,qtLoc
  REAL(RealKind) :: ThDensLoc,ThLoc,TLoc,ThEquivLoc
  REAL(RealKind) :: zPloc,r

  IF (Perturb) THEN
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Rho=>RhoCell(ibLoc)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    ThVirt=>VecC(ibLoc)%Vec(thPos)%c
    IF (EnPos>0) THEN
      En=>VecC(ibLoc)%Vec(EnPos)%c
    END IF
    DO iz=iz0+1,iz1
      zpLoc=zP(iz-1)+0.5d0*dz(iz)
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          RhoDry=RhoLoc-RhoVLoc-RhoLLoc
          rvLoc=RhoVLoc/(RhoDry+Eps)
          rlLoc=RhoLLoc/(RhoDry+Eps)
          rtLoc=rvLoc+rlLoc
          qvLoc=RhoVLoc/(RhoLoc+Eps)
          qlLoc=RhoLLoc/(RhoLoc+Eps)
          qtLoc=qvLoc+qlLoc
          KappaLoc=(Rd*(RhoLoc-RhoLLoc) &
                   +(Rv-Rd)*RhoVLoc) &
                   /(Cpd*RhoLoc+(Cpv-Cpd)*RhoVLoc &
                   +(Cpl-Cpd)*RhoLLoc)
          pLoc=p0*(Rd*ThVirt(ix,iy,iz,1)/p0)**(One/(One-KappaLoc))
          TLoc=pLoc/(Rd*RhoDry+Rv*RhoVLoc)
          IF (zpLoc<HeightDeltaT) THEN
            CALL Random_number(r)
            TLoc=TLoc+(r-0.5d0)*Two*DeltaT
          END IF
          IF (zpLoc<HeightDeltaQv) THEN
            CALL Random_number(r)
            qvLoc=qvLoc+(r-0.5d0)*Two*DeltaQv
            RhoV(ix,iy,iz,1)=qvLoc*Rho(ix,iy,iz,1)
          END IF
          RhoVLoc=RhoV(ix,iy,iz,1)
          rvLoc=qvLoc*RhoLoc/(RhoDry+Eps)
          rlLoc=qlLoc*RhoLoc/(RhoDry+Eps)
          rtLoc=rvLoc+rlLoc
          KappaLoc=(Rd*(RhoLoc-RhoLLoc) &
                   +(Rv-Rd)*RhoVLoc) &
                   /(Cpd*RhoLoc+(Cpv-Cpd)*RhoVLoc &
                   +(Cpl-Cpd)*RhoLLoc)
          SELECT CASE(ThetaKind)
          CASE('Density')
            ThVirt(ix,iy,iz,1)=RhoLoc*TLoc*(p0/pLoc)**KappaLoc*(One+rvLoc*Rv/Rd)/(One+rtLoc)
          CASE('PotTemp')
            ThVirt(ix,iy,iz,1)=TLoc*(p0/pLoc)**(Rd/Cpd)*RhoLoc
          CASE('Equiv')
            ThEquivLoc=ThEquiv(TLoc,RhoDry,RhoVLoc,RhoLLoc)
            ThVirt(ix,iy,iz,1)=ThEquivLoc*RhoLoc
          CASE('Energy')
            ThVirt(ix,iy,iz,1)= &
              (RhoDry*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl)*TLoc+RhoVLoc*L00 &
              +Half*RhoLoc*(zP(iz-1)+zP(iz))*Grav
          CASE('PreEn')
            ThVirt(ix,iy,iz,1)=pLoc
            En(ix,iy,iz,1)= &
              (RhoDry*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl)*TLoc+RhoVLoc*L00 &
              +Half*RhoLoc*(zP(iz-1)+zP(iz))*Grav+Half*RhoLoc*uMax*uMax
          END SELECT
        END DO
      END DO
    END DO
  END DO
  END IF
END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE Melpitz_Mod
  USE SoilData_Mod
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
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: uType
  LOGICAL, SAVE :: Load=.TRUE.
 
  IF (ProfInVel) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,uType,'uProf')
      Load=.FALSE.
    END IF
    UStart=ProfileEqual(cInt,z)
  ELSE
    UStart=UMax
  END IF
END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: vType
  LOGICAL, SAVE :: Load=.TRUE.
  IF (ProfInVel) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,vType,'vProf')
      Load=.FALSE.
    END IF
    VStart=ProfileEqual(cInt,z)
  ELSE
    VStart=vMax
  END IF
END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThProfFun=Zero
END FUNCTION ThProfFun

!========================================
FUNCTION ThStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThProfFun
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,thType,'ThDensProf')
      Load=.FALSE.
    END IF
    ThStart=ProfileEqual(cInt,z)
  ELSE
    S=N*N/Grav
    ThStart=thInit*exp(z*S)
  END IF
END FUNCTION ThStart

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Melpitz_Mod
  USE Start_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil,ThAir
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  IF (zSoil<=0.10d0) THEN
    ThStartSoil=TS05
  ELSE IF (zSoil<=0.20d0) THEN
    ThStartSoil=TS15
  ELSE IF (zSoil<=0.30d0) THEN
    ThStartSoil=TS25
  ELSE IF (zSoil<=0.40d0) THEN
    ThStartSoil=TS35
  ELSE IF (zSoil<=0.50d0) THEN
    ThStartSoil=TS45
  ELSE IF (zSoil<=0.60d0) THEN
    ThStartSoil=TS55
  ELSE IF (zSoil<=1.00d0) THEN
    ThStartSoil=TS80
  ELSE IF (zSoil<=2.00d0) THEN
    ThStartSoil=TS150
  ELSE 
    ThStartSoil=TS150
  END IF
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  IF (zSoil<=0.10d0) THEN
    QvStartSoil=cpwp(SoilType)+(cfcap(SoilType)-cpwp(SoilType))*BF10/100.0  
  ELSE IF (zSoil<=0.20d0) THEN
    QvStartSoil=cpwp(SoilType)+(cfcap(SoilType)-cpwp(SoilType))*BF20/100.0  
  ELSE IF (zSoil<=0.30d0) THEN
    QvStartSoil=cpwp(SoilType)+(cfcap(SoilType)-cpwp(SoilType))*BF30/100.0  
  ELSE IF (zSoil<=0.40d0) THEN
    QvStartSoil=cpwp(SoilType)+(cfcap(SoilType)-cpwp(SoilType))*BF40/100.0  
  ELSE IF (zSoil<=0.50d0) THEN
    QvStartSoil=cpwp(SoilType)+(cfcap(SoilType)-cpwp(SoilType))*BF50/100.0  
  ELSE 
    QvStartSoil=cpwp(SoilType)+(cfcap(SoilType)-cpwp(SoilType))*BF60/100.0  
  END IF
  IF (SoilType<=2) QvStartSoil=0.0d0
END FUNCTION QvStartSoil

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qvType,'QvProf')
      Load=.FALSE.
    END IF
    qvStart=ProfileEqual(cInt,z)
  ELSE
    qvStart=0.0d0
  END IF
END FUNCTION QvStart
!========================================


FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: Rad
  QcStart=Zero
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

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

  NvStart=300d6
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
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=0.0d0
END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE Parameter_Mod
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: RhoType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S,pLoc,ThLoc,L
  
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,RhoType,'RhoProf')
      Load=.FALSE.
    END IF
    RhoFun=ProfileEqual(cInt,z)
  ELSE
    S=N*N/Grav
    ThLoc=ThInit*exp(z*S)
    IF (N>Zero) THEN
      pLoc=p0*(One-Grav/(Cpd*ThInit*S)*(One-EXP(-S*z)))**(Cpd/Rd)
    ELSE
      pLoc=p0*(One-kappa*Grav*z/(Rd*ThInit))**(Cpd/Rd)
    END IF
    RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END IF
END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE Parameter_Mod
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc
  RhoProf=Zero
END FUNCTION RhoProf

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=Zero
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
 !TkeStart=0.01d0
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION OmeStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  OmeStart=Zero
! OmeStart=0.01d0
END FUNCTION OmeStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart


FUNCTION DummyStart(lam,phi,z,zHeight,Time)
  USE Melpitz_Mod
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE Melpitz_Mod
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE Melpitz_Mod
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE Melpitz_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE Melpitz_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: preType
  LOGICAL, SAVE :: Load=.TRUE.
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,preType,'PreProf')
      Load=.FALSE.
    END IF
    PreStart=ProfileEqual(cInt,z)
  ELSE
    PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
  END IF
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Melpitz_Mod
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

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=0.0d0
END FUNCTION ForceTh

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=0.0d0
END FUNCTION ForceRho

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  EnStart=Zero
END FUNCTION EnStart

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE Melpitz_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start

FUNCTION DampFun(z,Name)
  USE Melpitz_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
  DampFun=6.0d0
END FUNCTION DampFun
