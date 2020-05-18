MODULE Physics_Mod

  USE Kind_Mod
  USE Control_Mod
  USE DataType_Mod
  USE Transport_Mod

  IMPLICIT NONE

  TYPE RealBoundary_T
    REAL(RealKind) :: R_W
    REAL(RealKind) :: R_E
    REAL(RealKind) :: R_S
    REAL(RealKind) :: R_N
    REAL(RealKind) :: R_T
    REAL(RealKind) :: R_B
  END TYPE RealBoundary_T

  TYPE (VelocityFace_T), POINTER :: PGrad(:)
  TYPE (Vector4Cell_T), POINTER :: VecEnv1(:)
  TYPE (Vector4Cell_T), POINTER :: VecAmb(:)
  REAL(RealKind) :: Time1
  TYPE (Vector4Cell_T), POINTER :: VecEnv2(:)
  REAL(RealKind) :: Time2
  TYPE (ScalarCell_T), POINTER :: uCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: vCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: wCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: PreCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: SoundCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: RhoLCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: LWPCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: CpmlCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: HeightG(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: thProfG(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: PreKin(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DivCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: RhoProfG(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: PreProfG(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: RhoVProfG(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: RhoCProfG(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: PProfG(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: RhoCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: ThetaCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: TCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DummyCell(:)=>NULL()
  TYPE (Vector4Cell_T), POINTER :: DummyVec(:)=>NULL()
  TYPE (Vector4Cell_T), POINTER :: TAbsCell(:)=>NULL()
  TYPE(Vector4Cell_T), POINTER :: Act(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: KinEnCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: ECell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: SoSCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: EStartCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: PStartCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: qvsCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DampKoeffCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: ShadowCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: RaddirCell(:)=>NULL() ! direct solar radiation
  TYPE (ScalarCell_T), POINTER :: RaddifCell(:)=>NULL() ! diffuse solar radiation
  TYPE (ScalarCell_T), POINTER :: RadinfCell(:)=>NULL() ! Infrared radiation, namely the long wave radiation from the atmosphere
  TYPE (ScalarCell_T), POINTER :: SensFluxCell(:)=>NULL() ! Sensible heat flux
  TYPE (ScalarCell_T), POINTER :: LatFluxCell(:)=>NULL() ! Latent heat flux 
  TYPE (ScalarCell_T), POINTER :: HeatRateCell(:)=>NULL() ! heating rate from radiation parameterization
  TYPE (ScalarCell_T), POINTER :: ForceThCell(:)=>NULL() ! forcing term from nudging (theta)
  TYPE (ScalarCell_T), POINTER :: ForceRhoVCell(:)=>NULL() ! forcing term from nudging (qv)
  TYPE (ScalarCell_T), POINTER :: BulkCoeffDragCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: BulkCoeffHeatCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: BulkCoeffMoistCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: LandClassCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: RoughnessLengthCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: AlbedoCell(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffHKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffVKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffPotKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffPotHKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffPotVKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffMomKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffMomHKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: DiffMomVKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: LenKoeff(:)=>NULL()
  TYPE (ScalarCell_T), POINTER :: SediCell(:)=>NULL()
  REAL(RealKind), ALLOCATABLE :: PrandtlNumber(:)
  REAL(RealKind), ALLOCATABLE :: FallVelocity(:)
  REAL(RealKind), ALLOCATABLE :: MeanProfile(:,:) 
  TYPE (VelocityFace_T), POINTER :: DTUG(:)=>NULL()
  TYPE (VelocityFace_T), POINTER :: DUUG(:)=>NULL()
  TYPE (VelocityFace_T), POINTER :: DUJac(:)=>NULL()
  TYPE (Pressure_T), POINTER :: DUTG(:)=>NULL()
  TYPE (Pressure_T), POINTER :: DURG(:)=>NULL()
  TYPE (Pressure_T), POINTER :: DTTG(:)=>NULL()
  INTEGER :: NumMet 
  REAL(RealKind), ALLOCATABLE :: TimeEnvi(:) 

  LOGICAL :: Parcel
  LOGICAL :: TrajIn
  LOGICAL :: DilutParcel
  LOGICAL :: Shallow
  LOGICAL :: Liquid
  LOGICAL :: LiquidTam
  LOGICAL :: Height
  LOGICAL :: Advection
  LOGICAL :: Diffusion
  LOGICAL :: CrossDiff
  REAL(RealKind) :: DiffMin
  REAL(RealKind) :: DiffMax
  LOGICAL :: Transport
  LOGICAL :: PGradient
  LOGICAL :: Anelastic
  LOGICAL :: PseudoIn
  LOGICAL :: PressureUpdate
  CHARACTER*20 :: BousinesqType
  CHARACTER*20 :: ThetaKind
  CHARACTER*20 :: PreAdv
  CHARACTER*20 :: MicroScheme
  REAL(RealKind):: StartTimeIce
  LOGICAL :: JacMicro
  REAL(RealKind) :: FacAnela=One
  LOGICAL :: NoTke
  LOGICAL :: TkeSGS
  LOGICAL :: TkeSmag
  LOGICAL :: DynSmag
  LOGICAL :: TkeLen
  LOGICAL :: TkeDis
  LOGICAL :: TkeOme
  LOGICAL :: DTkeDis
  LOGICAL :: TkeDisRich
  LOGICAL :: DTkeDisRich
  LOGICAL :: TkeHVLen  
  LOGICAL :: Buoyancy
  LOGICAL :: Radiation
  LOGICAL :: SunMove
  LOGICAL :: RadiationProfile
  LOGICAL :: RadiationValues
  CHARACTER(100) :: RadiationValuesFile
  INTEGER :: RadiationValuesTime
  LOGICAL :: Damping
  LOGICAL :: DampAverage
  LOGICAL :: Forcing
  LOGICAL :: ForcingExtern
  LOGICAL :: ForcingExternNudging=.FALSE.
  LOGICAL :: ForcingExternTendency=.FALSE.
  LOGICAL :: ForcingExternSurface=.FALSE.
  CHARACTER*80 :: FileForceNudging=''
  CHARACTER*80 :: FileForceTendency=''
  CHARACTER*80 :: FileForceSurface=''
  LOGICAL :: ForcingCellPert
  INTEGER :: ForcingCellPertTime
  REAL(RealKind) :: ForcingCellPertTh
  REAL(RealKind) :: ForcingCellPertHeight
  CHARACTER(64) :: ForcingCellPertFile
  REAL(RealKind) :: ForcingStartTime
  INTEGER :: ForcingShift
  LOGICAL :: Subsidence
  INTEGER :: NumDamp
  TYPE(RealBoundary_T) :: Damp
  TYPE(RealBoundary_T) :: StrideDamp
  LOGICAL :: Cloud
  LOGICAL :: PrecipIce
  LOGICAL :: PrecipSnow
  LOGICAL :: PrecipRain
  LOGICAL :: RadiationISDAC
  LOGICAL :: Spectral
  REAL(RealKind) :: RelCloud
  REAL(RealKind) :: RelCloud2
  REAL(RealKind) :: RelCloudIce
  LOGICAL :: N_IN_Fix
  REAL(RealKind) :: N_IN
  REAL(RealKind) :: S_i_min
  REAL(RealKind) :: IceNucQc
  LOGICAL :: Coriolis
  LOGICAL :: CoriolisFree
  LOGICAL :: CoriolisFix
  LOGICAL :: CoriolisProfile
  REAL(RealKind) :: PhiCor
  REAL(RealKind) :: RotAngle
  LOGICAL :: Curvature
  LOGICAL :: Centrifugal
  LOGICAL :: Baum
  LOGICAL :: BaumFast
  LOGICAL :: DragSurf
  LOGICAL :: SensFluxFix
  REAL(RealKind) :: SensFluxConst
  LOGICAL :: LatFluxFix
  REAL(RealKind) :: LatFluxConst
  LOGICAL :: MomFlux
  LOGICAL :: ScalarFlux
  LOGICAL :: FluxDistribute
  LOGICAL :: DragScalar
  LOGICAL :: FreeSlip
  LOGICAL :: FluxSurf
  CHARACTER*10 :: FluxType
  REAL(RealKind) :: SensFluxAmp
  REAL(RealKind) :: SensFluxMin
  REAL(RealKind) :: LatFluxAmp
  REAL(RealKind) :: LatFluxMin
  REAL(RealKind) :: FluxShift
  REAL(RealKind) :: FluxClearness
  REAL(RealKind) :: FluxRandomSens
  REAL(RealKind) :: FluxRandomLat
  INTEGER :: FluxIxL
  INTEGER :: FluxIxR
  INTEGER :: FluxIyL
  INTEGER :: FluxIyR
  CHARACTER(32) :: IslandFile
  REAL(RealKind) :: SunRise
  REAL(RealKind) :: SunSet
  LOGICAL :: FireEmiss
  REAL(RealKind) :: FireStartTime
  REAL(RealKind) :: FireSensFlux
  INTEGER :: FireX0
  INTEGER :: FireX1
  INTEGER :: FireY0
  INTEGER :: FireY1
  REAL(RealKind) :: vSubConst
  CHARACTER(32) :: vSubProfile
  LOGICAL :: RainSurf
  LOGICAL :: IceSurf
  LOGICAL :: SnowSurf
  LOGICAL :: Wall
  LOGICAL :: Louis
  LOGICAL :: LouisMom
  LOGICAL :: DynamicSoil
  INTEGER :: LandClassDef
  LOGICAL :: Canopy
  LOGICAL :: CanopyEmission
  LOGICAL :: Wind
  LOGICAL :: WindTurbineTurn
  LOGICAL :: Traffic
  LOGICAL :: RealIsland
  CHARACTER*20 :: WindDragScheme
  CHARACTER*20 :: WindFarmFileName
  CHARACTER*20 :: DragScheme
  LOGICAL :: NoSlip
  LOGICAL :: SourceP
  LOGICAL :: Source
  LOGICAL :: Sphere
  LOGICAL :: Cylinder
  LOGICAL :: BetaPlane
  REAL(Realkind) :: OmegaPlane=2.0d-11
  CHARACTER(80) :: Profile    ! Name of Profile Input File


!!!!!!!! Photo
  REAL(RealKind) :: sss,cosz
  REAL(RealKind) :: cosSun
  REAL(RealKind) :: radn1,radn2,radn3
!  REAL(RealKind) :: raddiffus,raddirekt

  REAL(RealKind), ALLOCATABLE :: aphot(:,:)
  INTEGER :: nmadro

  NAMELIST /ModelPhysics/ Parcel,  &
                          TrajIn,  &
                          Profile, &
                          DilutParcel,  &
                          Shallow, &
                          Liquid,  &
                          LiquidTam,&
                          Height, &
                          Advection, &
                          Diffusion, &
                          CrossDiff, &
                          DiffMin,   &
                          DiffMax,   &
                          PGradient, &
                          Anelastic, &
                          PseudoIn, &
                          PressureUpdate, &
                          BousinesqType, &
                          ThetaKind,     &
                          PreAdv,   &
                          MicroScheme, &
                          StartTimeIce, &
                          JacMicro, &
                          NoTke,    &
                          TkeSGS,    &
                          TkeSmag,    &
                          DynSmag,    &
                          TkeLen,    &
                          TkeDis,    &
                          TkeOme,    &
                          TkeDisRich,    &
                          DTkeDis,   &
                          TkeHVLen,  &
                          Buoyancy,  &
                          Radiation, &
                          SunMove,   &
                          RadiationProfile, &
                          RadiationValues, &
                          RadiationValuesFile, &
                          RadiationValuesTime, &
                          Damping,   &
                          DampAverage,   &
                          Forcing,   &
                          ForcingExtern,   &
                          ForcingExternNudging,   &
                          ForcingExternTendency,   &
                          ForcingExternSurface,   &
                          FileForceNudging, &
                          FileForceTendency, &
                          FileForceSurface, &
                          ForcingCellPert,   &
                          ForcingCellPertTime,   &
                          ForcingCellPertFile,   &
                          ForcingCellPertTh, &
                          ForcingCellPertHeight, &
                          ForcingShift,   &
                          Subsidence,   &
                          NumDamp,   &
                          Damp,      &
                          StrideDamp,& 
                          Cloud,     &
                          PrecipIce, &
                          PrecipSnow, &
                          PrecipRain, &
                          RadiationISDAC, &
                          Spectral,  &
                          RelCloud,    & 
                          RelCloud2,    & 
                          RelCloudIce, &
                          N_IN_Fix, &
                          N_IN, &
                          S_i_min, &
                          IceNucQc, &
                          Coriolis,  &
                          CoriolisFree, &
                          CoriolisFix,&
                          CoriolisProfile,&
                          PhiCor,&
                          RotAngle, &
                          Sphere,&
                          Cylinder,&
                          BetaPlane,&
                          OmegaPlane,&
                          Curvature, &
                          Centrifugal, &
                          Baum, & 
                          BaumFast, & 
                          DragSurf,  &
                          LatFluxFix,  &
                          LatFluxConst, &
                          SensFluxFix,  &
                          SensFluxConst, &
                          MomFlux,  &
                          ScalarFlux,  &
                          FluxDistribute,  &
                          DragScalar,  &
                          FreeSlip,  &
                          FluxSurf,  &
                          FluxType,  &
                          FluxShift,  &
                          FluxClearness, &
                          FluxRandomSens,  &
                          FluxRandomLat,  &
                          FluxIxL,  &
                          FluxIxR,  &
                          FluxIyL,  &
                          FluxIyR,  &
                          ForcingStartTime, &
                          IslandFile, &
                          SunRise,  &
                          SunSet,   &
                          FireEmiss, &
                          FireStartTime, &
                          FireSensFlux, &
                          FireX0, &
                          FireX1, &
                          FireY0, &
                          FireY1, &
                          RainSurf,  &
                          IceSurf,  &
                          SnowSurf,  &
                          vSubConst, &
                          vSubProfile, &
                          Wall,      &
                          Louis,      &
                          LouisMom,      &
                          DynamicSoil,      &
                          LandClassDef, &
                          Canopy,      &
                          CanopyEmission,      &
                          Wind,         &
                          WindDragScheme,         &
                          WindFarmFileName,         &
                          WindTurbineTurn,         &
                          Traffic,   &
                          RealIsland,         &
                          DragScheme,      &
                          NoSlip    


CONTAINS

SUBROUTINE InitPrandtlNumber(FileName)

  CHARACTER(*) :: FileName

  REAL(RealKind) :: PrandtlTh
  INTEGER :: Pos
  CHARACTER(300) :: Line
  NAMELIST /Prandtl/ PrandtlTh
  PrandtlTh=0.7d0
! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&Prandtl')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=Prandtl)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)

! Prandtl number
  ALLOCATE(PrandtlNumber(VectorComponentsT))
  PrandtlNumber=One
  IF (disPos>0) THEN
    PrandtlNumber(disPos)=Karm*Karm/((Cmy2-Cmy1)*SQRT(Cmy0))
  END IF
  IF (ThPos>0) THEN
    PrandtlNumber(thPos)=PrandtlTh
  END IF
! FallVelocity
END SUBROUTINE InitPrandtlNumber

SUBROUTINE InputModelPhysics(FileName)

  CHARACTER(*) :: FileName
  
  INTEGER :: Pos
  CHARACTER(300) :: Line
  REAL(RealKind) :: deg2rad

  Parcel=.FALSE.
  TrajIn=.FALSE.
  Profile='Profile'
  DilutParcel=.FALSE.
  Shallow=.FALSE.
  Liquid=.FALSE.
  LiquidTam=.FALSE.
  Height=.FALSE.
  Advection=.FALSE.
  Diffusion=.FALSE.
  CrossDiff=.FALSE.
  DiffMin=1.d-4
  DiffMax=1.d+4
  Radiation=.FALSE.
  SunMove=.FALSE.
  RadiationProfile=.FALSE.
  RadiationValues=.FALSE.
  RadiationValuesFile=''
  RadiationValuesTime=600
  PGradient=.FALSE.
  Anelastic=.FALSE.
  PseudoIn=.FALSE.
  PressureUpdate=.FALSE.
  BousinesqType='Scaled'
  ThetaKind='Density'
  PreAdv='Outer'
  MicroScheme='Bulk2'
  StartTimeIce=0.0d0
  JacMicro=.FALSE.
  NoTke=.FALSE.
  TkeSGS=.FALSE.
  TkeSmag=.FALSE.
  DynSmag=.FALSE.
  TkeLen=.FALSE.
  TkeDis=.FALSE.
  DTkeDis=.FALSE.
  TkeOme=.FALSE.
  TkeDisRich=.FALSE.
  DTkeDisRich=.FALSE.
  TkeHVLen=.FALSE.
  Cloud=.FALSE.
  PrecipIce=.FALSE.
  PrecipSnow=.FALSE.
  PrecipRain=.FALSE.
  RadiationISDAC=.FALSE.
  Spectral=.FALSE.
  Damp%R_W=One
  Damp%R_E=One
  Damp%R_S=One
  Damp%R_N=One
  Damp%R_B=One
  Damp%R_T=One
  StrideDamp%R_W=Zero
  StrideDamp%R_E=Zero
  StrideDamp%R_S=Zero
  StrideDamp%R_N=Zero
  StrideDamp%R_B=Zero
  StrideDamp%R_T=Zero
  RelCloud=1.0d-1
  RelCloud2=1.0d-1
  RelCloudIce=1.0d-1
  N_IN_Fix=.FALSE.
  N_IN=4.0d0
  S_i_min=0.0d0
  IceNucQc=0.0d0
  Buoyancy=.FALSE.
  Coriolis=.FALSE.
  CoriolisFree=.FALSE.
  CoriolisFix=.TRUE.
  CoriolisProfile=.FALSE.
  Sphere=.FALSE.
  Cylinder=.FALSE.
  BetaPlane=.FALSE.
  PhiCor=1.0d-4
  RotAngle=Zero
  Curvature=.FALSE.
  Centrifugal=.FALSE.
  Baum=.FALSE. 
  BaumFast=.FALSE. 
  DragSurf=.FALSE.
  LatFluxFix=.FALSE.
  LatFluxConst=Zero
  SensFluxFix=.FALSE.
  SensFluxConst=Zero
  MomFlux=.TRUE.
  ScalarFlux=.TRUE.
  FluxDistribute=.FALSE.
  DragScalar=.FALSE.
  FreeSlip=.TRUE.
  FluxSurf=.FALSE.
  FluxType='CONST'
  SensFluxAmp=Zero
  SensFluxMin=Zero
  LatFluxAmp=Zero
  LatFluxMin=Zero
  FluxShift=0.0d0
  FluxClearness=One
  FluxRandomSens=Zero
  FluxRandomLat=Zero
  FluxIxL=0
  FluxIxR=1
  FluxIyL=0
  FluxIyR=1
  ForcingStartTime=0.0d0
  ForcingShift=0
  IslandFile='Barbados100x100.csv'
  SunRise=28800.0d0
  SunSet=66960.0d0
  FireEmiss=.FALSE.
  FireStartTime=3600.0d0
  FireSensFlux=75000.0d0
  FireX0=0
  FireX1=0
  FireY0=0
  FireY1=0
  RainSurf=.FALSE.
  IceSurf=.FALSE.
  SnowSurf=.FALSE.
  vSubConst=Zero
  vSubProfile='BOMEX'
  Wall=.FALSE.
  Louis=.FALSE.
  LouisMom=.FALSE.
  DynamicSoil=.FALSE.
  LandClassDef=5
  Canopy=.FALSE.
  CanopyEmission=.FALSE.
  DragScheme='revLouis'
  NoSlip=.FALSE.
  Damping=.FALSE.
  DampAverage=.FALSE.
  Forcing=.FALSE.
  ForcingExtern=.FALSE.
  ForcingCellPert=.FALSE.       ! turn cell perturbation method on/off
  ForcingCellPertTime=200       ! time betweed cell perturbation are applied
  ForcingCellPertFile=''        ! cell perturbation file with block information
  ForcingCellPertTh=0.5d0       ! pot. temperature amplitude
  ForcingCellPertHeight=565.0d0 ! height until perturbations are applied (2/3 z_i recommended)
  Subsidence=.FALSE.
  Wind=.FALSE.
  WindTurbineTurn=.FALSE.
  RealIsland=.FALSE.
  WindDragScheme='Disk'
  WindFarmFileName=''
  Traffic=.FALSE.

! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ModelPhysics')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelPhysics)
      EXIT
    END IF
  END DO
  CLOSE(UNIT=InputUnit)
1 CONTINUE

  CLOSE(UNIT=InputUnit)
  
  IF (Baum.AND.BaumFast) THEN
    Baum=.FALSE.
  END IF  
  RotAngle=RotAngle/180.0d0*Pi
  IF (.NOT.Coriolis) THEN
!   Omega=Zero  OSSI
  END IF
  IF (.NOT.Curvature) THEN
    Curv=Zero
  END IF
  IF (Anelastic.OR.PseudoIn) THEN
    FacAnela=Zero
  END IF
  IF (Sphere) THEN
    deg2rad=Pi/180.0d0
    StrideDamp%R_W=StrideDamp%R_W*deg2rad
    StrideDamp%R_E=StrideDamp%R_E*deg2rad
    StrideDamp%R_S=StrideDamp%R_S*deg2rad
    StrideDamp%R_N=StrideDamp%R_N*deg2rad
  END IF
  Transport=Advection.OR.Diffusion
  SourceP   =Buoyancy.OR.NoTke.OR.TkeSGS.OR.TkeSmag.OR.DynSmag.OR.TkeLen.OR.TkeDis.OR.TkeHVLen.OR. &
             Cloud.OR.Coriolis.OR.Curvature.OR.Centrifugal.OR.Baum.OR.DragSurf.OR.Wall.OR.      &
             Louis.OR.LouisMom.OR.NoSlip.OR.Damping
  Source=.TRUE.

  IF (MyId==0.AND.PrintNameLists) THEN
    WRITE(TermUnit,*) 'Parcel      = ',Parcel 
    WRITE(TermUnit,*) 'TrajIn      = ',TrajIn 
    WRITE(TermUnit,*) 'DilutParcel = ',DilutParcel 
    WRITE(TermUnit,*) 'Shallow     = ',Shallow
    WRITE(TermUnit,*) 'Liquid      = ',Liquid
    WRITE(TermUnit,*) 'LiquidTam   = ',LiquidTam
    WRITE(TermUnit,*) 'Height      = ',Height 
    WRITE(TermUnit,*) 'Advection   = ',Advection
    WRITE(TermUnit,*) 'Diffusion   = ',Diffusion
    WRITE(TermUnit,*) 'CrossDiff   = ',CrossDiff
    WRITE(TermUnit,*) 'Buoyancy    = ',Buoyancy
    WRITE(TermUnit,*) 'NoTke       = ',NoTke
    WRITE(TermUnit,*) 'TkeSGS      = ',TkeSGS
    WRITE(TermUnit,*) 'TkeSmag     = ',TkeSmag
    WRITE(TermUnit,*) 'DynSmag     = ',DynSmag
    WRITE(TermUnit,*) 'TkeLen      = ',TkeLen
    WRITE(TermUnit,*) 'TkeDis      = ',TkeDis
    WRITE(TermUnit,*) 'TkeDisRich  = ',TkeDisRich
    WRITE(TermUnit,*) 'TkeHVLen    = ',TkeHVLen  
    WRITE(TermUnit,*) 'Cloud       = ',Cloud
    WRITE(TermUnit,*) 'PrecipIce   = ',PrecipIce
    WRITE(TermUnit,*) 'PrecipSnow  = ',PrecipSnow
    WRITE(TermUnit,*) 'PrecipRain  = ',PrecipRain
    WRITE(TermUnit,*) 'Coriolis    = ',Coriolis
    WRITE(TermUnit,*) 'Sphere      = ',Sphere
    WRITE(TermUnit,*) 'Cylinder    = ',Cylinder
    WRITE(TermUnit,*) 'BetaPlane   = ',BetaPlane
    WRITE(TermUnit,*) 'Curvature   = ',Curvature
    WRITE(TermUnit,*) 'Centrifugal = ',Centrifugal
    WRITE(TermUnit,*) 'Baum        = ',Baum 
    WRITE(TermUnit,*) 'BaumFast    = ',BaumFast
    WRITE(TermUnit,*)
    WRITE(TermUnit,*) 'Radiation   = ',Radiation
    WRITE(TermUnit,*) 'SunMove     = ',SunMove  
    WRITE(TermUnit,*) 'RadiationPro= ',RadiationProfile
    WRITE(TermUnit,*) 'RadiationVal= ',RadiationValues
    WRITE(TermUnit,*) 'RadiationValFile= ',RadiationValuesFile
    WRITE(TermUnit,*) 'RadiationValTime= ',RadiationValuesTime
    WRITE(TermUnit,*) 'DragSurf    = ',DragSurf
    WRITE(TermUnit,*) 'SensFluxFix = ',SensFluxFix
    WRITE(TermUnit,*) 'LatFluxFix  = ',LatFluxFix
    WRITE(TermUnit,*) 'MomFlux     = ',MomFlux   
    WRITE(TermUnit,*) 'ScalarFlux  = ',ScalarFlux   
    WRITE(TermUnit,*) 'DragSacalar = ',DragScalar
    WRITE(TermUnit,*) 'FreeSlip    = ',FreeSlip
    WRITE(TermUnit,*) 'FluxSurf    = ',FluxSurf
    WRITE(TermUnit,*) 'FluxType    = ',FluxType
    WRITE(TermUnit,*) 'RainSurf    = ',RainSurf
    WRITE(TermUnit,*) 'IceSurf     = ',IceSurf
    WRITE(TermUnit,*) 'SnowSurf    = ',SnowSurf
    WRITE(TermUnit,*) 'Wall        = ',Wall
    WRITE(TermUnit,*) 'Louis       = ',Louis
    WRITE(TermUnit,*) 'LouisMom    = ',LouisMom
    WRITE(TermUnit,*) 'DynamicSoil = ',DynamicSoil
    WRITE(TermUnit,*) 'LandClassDef= ',LandClassDef
    WRITE(TermUnit,*) 'Canopy      = ',Canopy
    WRITE(TermUnit,*) 'CanopyEmission = ',CanopyEmission
    WRITE(TermUnit,*) 'Wind        = ',Wind
    WRITE(TermUnit,*) 'NoSlip      = ',NoSlip
    WRITE(TermUnit,*)
    WRITE(TermUnit,*) 'DTkeDis     = ',DTkeDis
    WRITE(TermUnit,*)
  END IF

!
! Parameter Check
!
! IF (Buoyancy.AND.thPos.EQ.0) THEN
!    STOP 'Error: Buoyancy with thPos=0'
! END IF
! IF (TkeDis.AND.tkePos*disPos.EQ.0) THEN
!    STOP 'Error: TkeDis with tke/disPos=0'
! END IF
! IF (Cloud.AND.RhoVPos*RhoCPos*RhoRPos.EQ.0) THEN
!    STOP 'Error: Cloud with qv/RhoC/RhoRPos=0'
! END IF
! IF (DTkeDis.AND.tkePos*disPos.EQ.0) THEN
!    STOP 'Error: DTkeDis with tke/disPos=0'
! END IF


END SUBROUTINE InputModelPhysics

END MODULE Physics_Mod
