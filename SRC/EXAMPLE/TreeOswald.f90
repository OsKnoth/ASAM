MODULE Tree_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: alpha=0.00d0
  REAL(RealKind) :: uMax=5.0d0
  REAL(RealKind) :: vMax=0.0d0 ,DMax
  REAL(RealKind) :: N=1.0d-2
  REAL(RealKind) :: th0=293.15d0
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0
  REAL(RealKind) :: AeroStart
  REAL(RealKind) :: DeltaT=0.00d0
  REAL(RealKind) :: QcProfMax=6.91457757446715515d-4
  REAL(RealKind) :: TkeMax=1.0d-2
  REAL(RealKind) :: DisMax=1.0d-4
  REAL(RealKind) :: zRauh=1.0d-1
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
  LOGICAL :: ProfIn=.FALSE.
  LOGICAL :: Perturb=.FALSE.
  LOGICAL :: LogWind=.FALSE.
  LOGICAL :: SineWind=.FALSE.
 
  NAMELIST /Example/ uMax &
                    ,vMax &
                    ,DMax &
                    ,AeroStart &
                    ,DeltaT &
                    ,QcProfMax &
                    ,ProfIn &
                    ,Perturb &
                    ,N &
                    ,th0 &
                    ,TkeMax &
                    ,DisMax &
                    ,zRauh &
                    ,LogWind &
                    ,SineWind

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
  
!  WRITE(*,*) 'x',offset_x,offset_y,inflow_lenx,inflow_leny

  IF (Time>1.0d0) THEN
    CALL Random_Number(r)
    uInFlow=0.0d0
    knull=18.0
    k_end = int(Inflow_leny/4.0)
      DO s = 1,10 !k_end
        k = ((10.0)*(0.9**s)) !*(1.0 + ((- 0.5d-1 + 1.0d-1*r)/((6.0*(0.8**1.0)))))
          DO ss = 1,2
            wx = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_x1-offset_x)
            wy = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_y1-offset_y)
            wz = ((-1.0)**ss)*k*(SQRT(uMax**Two+vMax**Two))/(offset_z1-offset_z)

            UStart3 = AMP*SIN(k*((x-offset_x)/(offset_x1-offset_x))*2.0*3.1415 + &
                  Time*(wx)*2.0*3.1415 + k**(2.0))  &
                + AMP*SIN(k*((y-offset_y)/(offset_y1-offset_y))*2.0*3.1415 + &
                  Time*(wy)*2.0*3.1415 + k**(2.0))  &
                + AMP*SIN(k*((z-offset_z)/(offset_z1-offset_z))*2.0*3.1415 + &
                  Time*(wz)*2.0*3.1415 + k**(2.0))
            uInFlow=uInFlow + UStart3
          ENDDO
      ENDDO
    uInFlow=(intenz*(SQRT((uMax**2.0)+(vMax**2.0)))*uInFlow/((s-1)*(ss-1))) + uMax
  ELSE
    uInFlow=UMax
  END IF

END FUNCTION uInFlow

END MODULE Tree_Mod



SUBROUTINE SetBoundCells(BoundCellLoc)

  USE Tree_Mod
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
  USE Tree_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: i,ix,iy,iz !,ib,ibloc
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: ThVirt(:,:,:,:)
  REAL(RealKind), POINTER :: En(:,:,:,:)
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDry
  REAL(RealKind) :: KappaLoc,pLoc,rvLoc,rlLoc,rtLoc
  REAL(RealKind) :: ThDensLoc,ThLoc,TLoc,ThEquivLoc
  REAL(RealKind) :: zPloc,r

  IF (Perturb) THEN
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Rho=>RhoCell(ib)%c
    RhoV=>VecC(ib)%Vec(RhoVPos)%c
    RhoL=>VecC(ib)%Vec(RhoCPos)%c
    RhoR=>VecC(ib)%Vec(RhoRPos)%c
    ThVirt=>VecC(ib)%Vec(thPos)%c
    IF (EnPos>0) THEN
      En=>VecC(ib)%Vec(EnPos)%c
    END IF 
    DO iz=iz0+1,iz1
      zpLoc=zP(iz-1)+0.5d0*dz(iz)
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          RhoDry=RhoLoc-RhoVLoc-RhoLLoc
          rvLoc=RhoVLoc/RhoDry+Eps
          rlLoc=RhoLLoc/RhoDry+Eps
          rtLoc=rvLoc+rlLoc
          KappaLoc=(Rd*(RhoLoc-RhoLLoc) &
                   +(Rv-Rd)*RhoVLoc) &
                   /(Cpd*RhoLoc+(Cpv-Cpd)*RhoVLoc &
                   +(Cpl-Cpd)*RhoLLoc)
          pLoc=p0*(Rd*ThVirt(ix,iy,iz,1)/p0)**(One/(One-KappaLoc))
          TLoc=pLoc/(Rd*RhoDry+Rv*RhoVLoc)
          IF (zpLoc<600.d0) THEN
            CALL Random_number(r)
            TLoc=TLoc+(r-0.5d0)*DeltaT
          END IF
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
  USE Tree_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

! Find line
  DMax=75.0d0
  uMax=0.0d0
  vMax=0.0d0

  Inflow_lenx=domain%nx
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

FUNCTION UStart(x,y,z,zHeight,Time)

  USE Tree_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart
  REAL(RealKind) :: uStar=0.4

  REAL(RealKind) :: x,y,z,zHeight,Time

! CM: Distinction into either a logarithmic wind profile or a continuous wind field over the whole vertical column
  IF (LogWind) THEN
    UStart=((log((z+1.d-1)/zRauh))/(log(100.0/1d-1)))*uMax !*uMax
!    UStart=uStar/0.4*LOG((z+1.0d-1)/zRauh)
  ELSE IF (SineWind) THEN
    UStart=ABS(COS(alpha*PI/180.0d0))*uInFlow(x,y,z,Time)
  ELSE
    UStart=uMax
  END IF

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE Tree_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  VStart=VMax

END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE Tree_Mod
  USE ReadProfile_Mod

  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.

  REAL(RealKind) :: pLoc,TLoc

  ThProfFun=Zero

END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Tree_Mod
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
  REAL(RealKind) :: S
  S=N*N/Grav
  ThStart=th0*exp(z*S)

END FUNCTION ThStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Tree_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: RhoType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S,pLoc,ThLoc,L
  S=N*N/Grav
  ThLoc=th0*exp(z*S)
  IF (N>Zero) THEN
    pLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
  ELSE
    pLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
  END IF
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
END FUNCTION RhoFun


FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: tType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,tType,'PreProf')
      Load=.FALSE.
    END IF
    PreStart=ProfileEqual(cInt,z)
  ELSE
    PreStart=p0*(One-kappa*Grav*z/(Rd*300.0d0))**(Cpd/Rd)
  END IF

END FUNCTION PreStart


FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=Zero
END FUNCTION PreProf


FUNCTION TFun(x,y,z,zHeight,Time)
  USE Tree_Mod
  USE Rho_Mod
  USE PreProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: TFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: pLoc,RhoLoc,GammaDry
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: tType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,tType,'TProf')
      Load=.FALSE.
    END IF
    TFun=ProfileEqual(cInt,z)
  ELSE
    TFun=300.0d0
  END IF

END FUNCTION TFun


FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Tree_Mod
  USE Rho_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: TLoc,RhoLoc

  QvProfFun=0.0d0

END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Tree_Mod
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

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Tree_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qcType,'QcProf')
      Load=.FALSE.
    END IF
    qcStart=ProfileEqual(cInt,z)
  ELSE
    QcStart=0.0d0
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

  USE Tree_Mod
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
  NvStart=AeroStart*(One-QcLoc/QcProfMax)

END FUNCTION NvStart

FUNCTION NcStart(x,y,z,zHeight,Time)

  USE Tree_Mod
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
  NcStart=QcLoc/QcProfMax*AeroStart

END FUNCTION NcStart

FUNCTION NrStart(x,y,z,zHeight,Time)

  USE Tree_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NrStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NrStart=Zero

END FUNCTION NrStart

FUNCTION NiStart(x,y,z,zHeight,Time)

  USE Tree_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NiStart=Zero

END FUNCTION NiStart

FUNCTION NsStart(x,y,z,zHeight,Time)

  USE Tree_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NsStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NsStart=Zero

END FUNCTION NsStart

FUNCTION OmeStart(lam,phi,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE 
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  OmeStart = 0.0
END FUNCTION OmeStart

FUNCTION DStart(x,y,z,zHeight,Time)

  USE Tree_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=DMax

END FUNCTION DStart


FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Tree_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,TLoc,ThLoc

  RhoProf=Zero

END FUNCTION RhoProf


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Tree_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Tree_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=(RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time))/RhoFun(x,y,z,zHeight,Time)
END FUNCTION RhoStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=Zero
END FUNCTION DummyStart

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer1Start=Zero
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE Tree_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  Tracer2Start=Zero
END FUNCTION Tracer2Start


FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE Tree_Mod
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
  USE Tree_Mod
  USE Kind_Mod
  USE Rho_Mod
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: fTime,fSpace 
  ForceRho=0.0d0
END FUNCTION ForceRho

FUNCTION ForceTh(x,y,z,zHeight,Time)
  USE Tree_Mod
  USE Kind_Mod
  USE Rho_Mod
  REAL(RealKind) :: ForceTh
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceTh=0.0d0
END FUNCTION ForceTh

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  ThStartSoil=0.0d0
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE Parameter_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,z,zHeight,zSoil,Time
  INTEGER :: LandClass,SoilType
  QvStartSoil=0.0d0
END FUNCTION QvStartSoil

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE Tree_Mod
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE Tree_Mod
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE Tree_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE Tree_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE Tree_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  EnStart=Zero

END FUNCTION EnStart

FUNCTION DampFun(z,Name)
  USE Tree_Mod
  REAL(RealKind) :: DampFun
  REAL(RealKind) :: z
  CHARACTER(*) :: Name
  DampFun=0.0d0
END FUNCTION DampFun
