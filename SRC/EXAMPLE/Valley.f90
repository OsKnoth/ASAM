MODULE Valley_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE SoilData_Mod

  IMPLICIT NONE 


  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: ThSoil=283.0d0
!  REAL(RealKind) :: UMax=5.0d0
  REAL(RealKind) :: QvSoil=0.2 !1.0d0
  REAL(RealKind) :: DeltaT=0.0d0
  REAL(RealKind) :: HeightDeltaT=0.0d0
  REAL(RealKind) :: DeltaQv=0.0d0
  REAL(RealKind) :: HeightDeltaQv=0.0d0

  REAL(RealKind) :: xC,yC,zC
  REAL(RealKind) :: xr,yr,zr
  REAL(RealKind) :: DeltaTh,uMax,vMax,ThInit
  REAL(RealKind) :: qc0,qv0,qt0
  REAL(RealKind) :: RelHum !,N=0.0d0
  LOGICAL :: ProfIn=.FALSE.
  LOGICAL :: Perturb=.TRUE.
  
  NAMELIST /Example/ N   &
                    ,uMax &
                    ,vMax &
                    ,ThSoil   &
                    ,QvSoil   &
                    ,xC &
                    ,yC &
                    ,zC &
                    ,xr &
                    ,yr &
                    ,zr &
                    ,qv0 &
                    ,qc0 &
                    ,qt0 &
                    ,RelHum &
                    ,ThInit &
                    ,DeltaTh &
                    ,ProfIn &
                    ,Perturb


END MODULE Valley_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Valley_Mod
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
  USE Valley_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
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
          IF (VolC(ix,iy,iz)>Zero) THEN
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
            rvLoc=qvLoc*RhoLoc/(RhoDry+Eps)
            rlLoc=qlLoc*RhoLoc/(RhoDry+Eps)
            rtLoc=rvLoc+rlLoc
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
          END IF
        END DO
      END DO
    END DO
  END DO
  END IF
END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE Valley_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos 
  CHARACTER(300) :: Line
  INTEGER :: i
  REAL(RealKind) :: r

! Find line
  uMax=0.0d0
  vMax=0.0d0
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

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Valley_Mod
  USE Rho_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time

!  RhoProf=RhoFun(x,y,z,zHeight,Time)
  RhoProf=0.0d0

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE Valley_Mod
  USE Start_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  ThProfFun=0.0d0

END FUNCTION ThProfFun

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE Parameter_Mod
  USE Valley_Mod
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

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  USE Rho_Mod
  USE QvProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad,pLoc,TLoc
  REAL(RealKind) :: RhoLoc,RhoDLoc,RhoVLoc,RhoLLoc
  REAL(RealKind) :: rvLoc,rlLoc 
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      SELECT CASE(ThetaKind)
      CASE('Density')
        CALL ReadProfile(cInt,thType,'ThDensProf')
      CASE('Equiv')
        CALL ReadProfile(cInt,thType,'ThetaEProf')
      CASE('Energy')
        CALL ReadProfile(cInt,thType,'EnergyProf')
      CASE('PreEn')
        CALL ReadProfile(cInt,thType,'PreProf')
      CASE DEFAULT
        WRITE(*,*) 'ThetaKind pruefen.' 
      END SELECT
      Load=.FALSE.
    END IF
    ThStart=ProfileEqual(cInt,z)
  ELSE
    ThStart=0.0d0
  END IF

END FUNCTION ThStart

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Valley_Mod
  USE Start_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil,ThAir
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  REAL(RealKind) :: Rm,Cpml,KappaLoc
  REAL(RealKind) :: RhoDLoc,RhoLoc,RhoVLoc,RhoLLoc
  RhoLoc=RhoFun(x,y,z,zHeight,Time)
  RhoVLoc=QvStart(x,y,z,zHeight,Time)
  RhoLLoc=QcStart(x,y,z,zHeight,Time)
  RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
  Rm=Rd*RhoDLoc+Rv*RhoVLoc
  Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc
  KappaLoc=Rm/Cpml
  ThAir=ThStart(x,y,z,zHeight,Time)
  ThAir=(Rd*RhoLoc*ThAir/p0**KappaLoc)**(One/(One-KappaLoc))/(Rd*RhoDLoc+Rv*RhoVLoc)
  IF (zSoil<=0.05d0) THEN
!    ThStartSoil=280.50d0
    ThStartSoil=ThAir
  ELSE IF (zSoil<=0.5d0) THEN
    ThStartSoil=ThAir-0.5d0
!    ThStartSoil=280.00d0
  ELSE IF (zSoil<=1.0d0) THEN
    ThStartSoil=ThAir-2.0d0
!    ThStartSoil=279.50d0
  ELSE
    ThStartSoil=ThAir-3.5d0
!    ThStartSoil=278.00d0
  END IF
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  IF (zSoil==0.0d0) THEN
    QvStartSoil=0.d0  ! interception reservoir
    QvStartSoil=MIN(QvStartSoil,5.d-4)
  ELSE
    QvStartSoil=0.0d0
    QvStartSoil=MIN(QvStartSoil,cporv(SoilType))
  END IF
  IF (SoilType<=2) QvStartSoil=0.0d0
END FUNCTION QvStartSoil

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE Valley_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qvType,'QvProf')
      Load=.FALSE.
    END IF
    qvProfFun=ProfileEqual(cInt,z)
  ELSE
   qvProfFun=0.0d0
  END IF
! qvProfFun=0.0d0

END FUNCTION QvProfFun


FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

  qcStart=0.0d0

END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=Zero

END FUNCTION QiStart

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: uType
  LOGICAL, SAVE :: Load=.TRUE.
  UStart=uMax
END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: vType
  LOGICAL, SAVE :: Load=.TRUE.
  VStart=VMax
END FUNCTION VStart

FUNCTION DStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=Zero

END FUNCTION DStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Valley_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Valley_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(Cpd/Rd)
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(Cpd/Rd)
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=0.0d0
END FUNCTION TStart

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=Zero
END FUNCTION

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=Zero
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=Zero
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=Zero
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=Zero
END FUNCTION ForceRho

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=Zero
END FUNCTION ForceTh

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1 = 0.0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2 = 0.0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE Valley_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3 = 0.0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE Valley_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  EnStart=Zero
END FUNCTION EnStart

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  USE QvProf_Mod
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
    QvStart=0.0d0
  END IF

END FUNCTION QvStart

FUNCTION QsStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QsStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QsStart=Zero

END FUNCTION QsStart


FUNCTION NvStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NvStart
  REAL(RealKind) :: QcLoc

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

  QcLoc=Zero
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qcType,'QcProf')
      Load=.FALSE.
    END IF
    QcLoc=ProfileEqual(cInt,z)
  END IF
  NvStart=100.0d6

END FUNCTION NvStart

FUNCTION NcStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NcStart
  REAL(RealKind) :: QcLoc

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

  QcLoc=Zero
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qcType,'QcProf')
      Load=.FALSE.
    END IF
    QcLoc=ProfileEqual(cInt,z)
  END IF
  NcStart=Zero

END FUNCTION NcStart

FUNCTION NrStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NrStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NrStart=Zero

END FUNCTION NrStart

FUNCTION NiStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NiStart=Zero

END FUNCTION NiStart

FUNCTION NsStart(x,y,z,zHeight,Time)

  USE Valley_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NsStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NsStart=Zero

END FUNCTION NsStart

FUNCTION OmeStart(lam,phi,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  OmeStart = 0.0
END FUNCTION OmeStart

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE Valley_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start
