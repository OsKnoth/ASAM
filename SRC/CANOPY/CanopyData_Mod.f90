MODULE CanopyData_Mod
  USE Kind_Mod
  USE Control_Mod


IMPLICIT NONE
  INTEGER, PARAMETER ::                     &
    ireals    = SELECTED_REAL_KIND (12,200),&
                     ! number of desired significant digits for
                     ! real variables
                     ! corresponds to 8 byte real variables

    iintegers = KIND  (1)
                     ! kind-type parameter of the integer values
                     ! corresponds to the default integers

  REAL(KIND=ireals) ::  &
!   a) parameters describing the soil water budget
    cporv (1:8) , &   ! pore volume (fraction of volume)
    cfcap (1:8) , &   ! field capacity (fraction of volume)
    cpwp  (1:8) , &   ! plant wilting point (fraction of volume)
    cadp  (1:8) , &   ! air dryness point (fraction of volume)
    cik2  (1:8) , &   ! minimum infiltration rate (kg/s*m**2)
    ckw0  (1:8) , &   ! parameter for determination of hydr. conductivity (m/s)
    ckw1  (1:8) , &   ! parameter for determination of hydr. conductivity (1)
    cdw0  (1:8) , &   ! parameter for determination of hydr. diffusivity (m**2/s)
    cdw1  (1:8) , &   ! parameter for determination of hydr. diffusivity (1)
    smPot (1:8) , &   ! parameter for determination of adhesive and cohesive (capillary) forces between water and soil

!   b) parameters describing the soil heat budget
    crhoc (1:8) , &   ! soil heat capacity  (J/K*m**3)
    cala0 (1:8) , &   ! parameters for the determination of
    cala1 (1:8) , &   ! the soil heat conductivity (W/(K*m))

!   c) additional parameters for the BATS scheme (Dickinson)
    cbedi (1:8) , &   ! (1)
    cres  (1:8) , &   
    alb0  (1:8) , &   ! albedo for dry soil 
    dalb  (1:8) , &   ! difference between albedo of dry soil and of wet soil
    lamb  (1:8) , &   ! Werte nach Rawls
    emiss0(1:8) , &   ! emissivity of dry soil 
    demiss(1:8) , &   ! difference between emissivity of dry soil and of wet soil 

!   d) Parameters depending on land use 
    PlantCover_Max (0:9), & ! Plant cover during summer
    PlantCover_Min (0:9), & ! Plant cover during winter
    LAI_Max (0:9), & ! Leaf Area Index during summer
    LAI_Min (0:9), & ! Leaf Area Index during winter
    r_depth (0:9), & ! root depth [m]
    aroot   (0:9), &
    broot   (0:9), &
    zrough  (0:9), &
    AlbMax  (0:9), &
    AlbMin  (0:9), &
    zCanopy (0:9), &
    zStem   (0:9)
    

!   soil type:     ice    rock    sand    sandy   loam    clay     clay    peat   
!   (by index)                            loam            loam                     
  DATA  cporv   / 1.E-10, 1.E-10,  .364 ,  .445 ,  .455 ,  .475 ,  .507 ,  .863  /  ! w_pv volume of voids [1] [m³/m³]
  DATA  cfcap   / 1.E-10, 1.E-10,  .196 ,  .260 ,  .340 ,  .370 ,  .463 ,  .763  /  ! w_FV field capacity [1] [m³/m³]
  DATA  cpwp    / 0.0   , 0.0   ,  .042 ,  .100 ,  .110 ,  .185 ,  .257 ,  .265  /  ! w_PWP permanent wilting point[1][m³/m³]
  DATA  cadp    / 0.0   , 0.0   ,  .012 ,  .030 ,  .035 ,  .060 ,  .065 ,  .098  /  ! w_ADP air dryness point [1]
  DATA  cik2    / 0.0   , 0.0   , 0.0035, 0.0023, 0.0010, 0.0006, 0.0001, 0.0002 /  ! I_K2 min. infiltration rate [kg/(m²s)]
  DATA  cdw0    / 0.0   , 0.0   , 184d-7, 346d-8, 357d-8, 118d-8, 442d-9, 106d-9 /  ! D_0 hydraulic diffus.param.[m²/s]
  DATA  cdw1    / 0.0   , 0.0   , -8.45 ,  -9.47, -7.44 , -7.76 , -6.74 , -5.97  /  ! D_1 hydraulic diffus.param.[1]
  DATA  ckw0    / 0.0   , 0.0   , 479d-7, 943d-8, 531d-8, 764d-9,  17d-9,  58d-9 /  ! K_0 hydraulic conduct.param.[m/s]
  DATA  ckw1    / 0.0   , 0.0   , -19.27, -20.86, -19.66, -18.52, -16.32, -16.48 /  ! K_1 hydraulic conduct.param.[1]
  DATA  crhoc   / 1.92E6, 2.10E6, 1.28d6, 1.35d6, 1.42d6, 1.50d6, 1.63d6, 0.58d6 /  ! rho_0 c_0 heat capacity [J/(m³K)]
  DATA  cala0   / 2.26  , 2.41  , 0.30  ,  0.28 ,  0.25 ,  0.21 ,  0.18 ,  0.06  /  ! lambda_0 heat conductivity [W/(K m)]
  DATA  cala1   / 2.26  , 2.41  , 2.40  ,  2.40 ,  1.58 ,  1.55 ,  1.50 ,  0.50  /  ! delta lambda heat conductivity[W/(K m)]
  DATA  cbedi   / 1.00  , 1.00  ,  3.5  ,  4.8  , 6.1   , 8.6   , 10.0  , 9.0    /  ! Exponent B [1]
  DATA  lamb    / 0.00  , 0.00  , 0.694 , 0.379 , 0.252 , 0.242 , 0.165 , 0.129  /  ! Exponent B [1], Rawls/Brakensiek
  DATA  smPot   / 1.0d99, 1.0d99, -0.121, -0.218, -0.478, -0.630, -0.450,-0.356  /  ! Soil Saturated Matric Potential[m], Tjernström
  DATA  cres    / 0.0   , 0.0   , 0.044 , 0.031 , 0.052 , 0.092 , 0.075 , 0.110  /  ! Residual water content, Rawls (Schaap Leij 1998)
  DATA  alb0    / 0.50  , 0.22  , 0.35  , 0.32  , 0.30  , 0.28  , 0.25  , 0.15   /  ! Albedo des trockenen Bodens
  DATA  dalb    / 0.10  , 0.12  , 0.15  , 0.14  , 0.13  , 0.11  , 0.10  , 0.10   /  ! Albedodifferenz zw. trockenem + feuchten Boden
  DATA  emiss0  / 0.0   , 0.0   , 0.9   , 0.9   , 0.9   , 0.9   , 0.9   , 0.98   /  ! Emissionsvermögen des trockenen Bodens
  DATA  demiss  / 0.0   , 0.0   , 0.08  , 0.08  , 0.08  , 0.08  , 0.08  , 0.0    /  ! Emissionsdifferenz zw. trockenem + feuchtem Boden

!  DATA  cporv   / 1.E-10, 1.E-10,  .415 ,  .389 ,  .354 ,  .441 ,  .451 ,  .863  /  ! w_pv volume of voids [1] [m³/m³] Rawls (Schaap Leij 1998)
!  DATA  cporv   / 1.E-10, 1.E-10,  .437 ,  .453 ,  .463 ,  .464 ,  .475 ,  .863  /  ! w_pv volume of voids [1] [m³/m³] Rawls/Brakensiek
!  DATA  ckw0    / 0.0   , 0.0   ,5.94E-5,3.92E-6,1.13E-6,5.41E-7,1.01E-6,  8E-6  /  ! K_0 hydraulic conduct.param.[m/s] Rawls (Schaap Leij 1998)
!  DATA  ckw0    / 0.0   , 0.0   ,1.76E-4, 3.4E-3,   7E-6, 2.5E-6, 1.3E-6,  8E-6  /  ! K_0 hydraulic conduct.param.[m/s] Tjernström
!  DATA  ckw0    / 0.0   , 0.0   , 583E-7, 720E-8, 190E-8, 639E-9, 167E-9, 800E-8 /  ! K_0 hydraulic conduct.param.[m/s] Rawls/Brakensiek
!  DATA  lamb    / 1.00  , 1.00  ,  1.88 , 0.413 , 0.549 , 0.549 , 0.288 , 0.129  /  ! Exponent B [1], Rawls et al
!  DATA  cbedi   / 1.00  , 1.00  ,  4.05 ,  4.9  , 5.39  , 8.52  , 11.4  , 9.0    /  ! Exponent B [1]
!  DATA  smPot   / 0.0   , 0.0   ,-0.3715,-0.3715,-1.3183,-0.8913,-0.8511,-2.2909 /  ! Soil Saturated Matric Potential[m], Rawls (Schaap Leij 1998)
!  DATA  smPot   / 0.0   , 0.0   ,-0.1598,-0.3020,-0.4012,-0.5643,-0.8560,-0.3256 /  ! Soil Saturated Matric Potential[m], Rawls/Brakensiek
!  DATA  cres    / 0.0   , 0.0   , 0.020 , 0.041 , 0.027 , 0.075 , 0.090 , 0.089  /  ! Residual water content, Rawls/Brakensiek

!  Land use:	         bare   urban  savannah  deciduous  coniferous mixed              annual grass   sea
!             	         soil   area             forest     forest     forest  shrubland  crops  land  
  DATA LAI_Max         /  0.0  , 4.7 ,  3.0    ,   6.0   ,  8.0      ,  7.0   ,  4.0    ,  5.0   , 4.5  , 0.0  /
  DATA LAI_Min         /  0.0  , 0.1 ,  1.0    ,   0.0   ,  8.0      ,  4.0   ,  0.1    ,  0.2   , 0.5  , 0.0  /
  DATA PlantCover_Max  /  0.0  ,0.1  ,  0.8    ,   1.0   ,  0.9      ,  1.0   ,  0.4    ,  1.0   , 0.95 , 0.0  /
  DATA PlantCover_Min  /  0.0  ,0.05 ,  0.2    ,   0.0   ,  0.9      ,  0.5   ,  0.1    ,  0.5   , 0.9  , 0.0  /
!  DATA r_depth         / 0.0  ,0.6  ,  1.0    ,   1.0   ,  0.6      ,  0.8   ,  0.4    ,  0.8   , 0.6   /
  DATA r_depth         /  0.0  ,1.5  ,  2.4    ,   2.0   ,  1.8      ,  2.4   ,  3.1    ,  1.5   , 1.5  , 0.0  /
  DATA aroot           /  0.0  ,5.558,  8.235  ,   5.990 ,  6.706    ,  4.453 ,  7.718  ,  5.558 ,10.74 , 0.0  /    
  DATA broot           /  0.0  ,2.614,  1.627  ,   1.955 ,  2.175    ,  1.631 ,  1.262  ,  2.614 ,2.608 , 0.0  / 
  DATA zrough          /  0.05 ,1.0  ,  0.15   ,   1.0   ,  1.0      ,  1.0   ,  0.2    ,  0.5   , 0.03 , 0.01 /
  DATA AlbMax          /  0.0  ,0.16 ,  0.15   ,   0.15  ,  0.16     ,  0.16  ,  0.20   ,  0.20  , 0.20 , 0.08 /
  DATA AlbMin          /  0.0  ,0.14 ,  0.15   ,   0.11  ,  0.14     ,  0.12  ,  0.17   ,  0.16  , 0.16 , 0.08 /
  DATA zCanopy         /  0.0  ,5.0  ,  1.00   ,   15.0  ,  30.0     ,  23.0  ,  2.00   ,  1.50  , 0.50 , 0.00 /
  DATA zStem           /  0.0  ,2.5  ,  0.00   ,   4.00  ,  20.0     ,  10.0  ,  0.50   ,  0.00  , 0.00 , 0.00 /

  REAL(RealKind), PARAMETER :: zS(0:8)= (/ 0.0d0, 0.01d0, 0.03d0, 0.09d0, 0.27d0, 0.81d0, 2.43d0, 7.29d0, 21.87d0 /)

  REAL(RealKind)            :: Cdash      = 0.05d0  !  BATS ((m/s)**1/2)
  REAL(RealKind), PARAMETER :: RhoW       = 1.d3       ! [kg/m³] Density of liquid Water
  REAL(RealKind), PARAMETER :: CRhoW      = 4.18d6     ! J/(m^3K)
  REAL(RealKind), PARAMETER :: TSClim     = 283.0d0    ! [K] climatological layer temperature
  REAL(RealKind), PARAMETER :: capsoil    = 1.260d6    ! J m^¯3 K^-1 volumetric heat capacity of soil
  REAL(RealKind), PARAMETER :: capwater   = 4.18d6     ! J m^¯3 K^-1 volumetric heat capacity of water
  REAL(RealKind), PARAMETER :: capair     = 1.298d3    ! J m^¯3 K^-1 volumetric heat capacity of air 

  REAL(RealKind), PARAMETER :: ResMin     = 150.0d0    ! [s/m]
  REAL(RealKind), PARAMETER :: ResMax     = 4000.0d0   ! [s/m]
  REAL(RealKind), PARAMETER :: RadPARcrit = 100.0d0    ! [W/m²] photosynthetically active radiation
  REAL(RealKind), PARAMETER :: TempEnd    = 313.15d0   ! [K] 
  REAL(RealKind), PARAMETER :: Temp0      = 273.15d0   ! [K] Temperature at freezing point
  REAL(RealKind), PARAMETER :: Epotnorm   = 4.75       ! mm/d

  REAL(RealKind), PARAMETER :: sigma      = 5.671d-8   ! Wm^(-2)K^(-4) Stefan-Boltzmann-Konstante
  REAL(RealKind), PARAMETER :: Kg         = 1.818d-5   ! s^(-1) Relaxation Constant
  REAL(RealKind), PARAMETER :: tau_wg     = 86400.0d0   !86400.d0	! restoring timescale (24h)
  REAL(RealKind), PARAMETER :: tau_perc   = 1000.0d0   ! time constant > 2*dtMax 
  REAL(RealKind), PARAMETER :: Cw         = 2.0d-2     ! [m], moist capacity of the ground
  REAL(RealKind), PARAMETER :: bb         = 5.0d0      ! bb=cc=dd 
  REAL(RealKind), PARAMETER :: cc         = 5.0d0      ! free parameters of the stability functions 
  REAL(RealKind), PARAMETER :: dd         = 5.0d0      

  REAL(RealKind), PARAMETER :: Prec       = 0.d0       ! m -> 10^-6m=1mm Precipitaion

  REAL(RealKind)            :: VegBegin
  REAL(RealKind)            :: VegPeriod
  REAL(RealKind)            :: VegFactor
  REAL(RealKind)            :: HeightReduct 
  REAL(RealKind), PARAMETER :: z_rmin     = 0.12d0
  REAL(RealKind), PARAMETER :: z_rmax     = 0.70d0

  CHARACTER*40              :: SoilParam


CONTAINS

SUBROUTINE DataCompute
  REAL(RealKind) :: latitude

  latitude     = lat*180.d0/Pi
  VegBegin     = MAX(1.0d0,3.0d0*(ABS(latitude)-20.0d0))
  VegPeriod    = MIN(365.0d0,345.0d0-4.5d0*(ABS(latitude)-20.0d0))
  VegFactor    = MAX(0.0d0,MIN(1.0d0,1.12d0*SIN(Pi*MAX(0.0d0,(StartDay-VegBegin)/VegPeriod))))
  HeightReduct = 1.0d0 ! EXP(-5.0d-9*GeopotentialHeightOfSurface**2.0d0)
END SUBROUTINE DataCompute

FUNCTION z_root(LandUse)
  REAL(RealKind) :: z_root,rootdepth
  INTEGER        :: LandUse
  CALL DataCompute
  rootdepth = r_depth(LandUse)
  z_root    = MIN(rootdepth,z_rmin+(z_rmax-z_rmin)*VegFactor**2.0d0)
END FUNCTION

FUNCTION xroot(LandUse,rootdepth)
  REAL(RealKind) :: xroot
  REAL(RealKind) :: rootdepth
  INTEGER        :: LandUse
  xroot = 1.0d0-0.5d0*(EXP(-aroot(LandUse)*rootdepth)+EXP(-broot(LandUse)*rootdepth))
END FUNCTION xroot

FUNCTION LeafAreaIndex(LandUse)
  REAL(RealKind) :: LeafAreaIndex
  INTEGER        :: LandUse
  CALL DataCompute
  LeafAreaIndex = LAI_Min(LandUse)+(LAI_Max(LandUse)-LAI_Min(LandUse))  &
                                   *VegFactor*HeightReduct
END FUNCTION

FUNCTION LeafAreaDensity(LandUse,z1,z2,dz)
  REAL(RealKind) :: LeafAreaDensity
  REAL(RealKind) :: z1,z2,dz
  INTEGER        :: LandUse
  IF (LandUse/=0) THEN
    LeafAreaDensity = LeafAreaIndex(LandUse)/(zCanopy(LandUse)-zStem(LandUse))* &
                      (1.0d0-( MAX(MIN(zStem(LandUse)-z1,dz)/(z2-z1),0.0d0) + &
                               MAX(MIN(z2-zCanopy(LandUse),dz)/(z2-z1),0.0d0) )) 
  ELSE
    LeafAreaDensity = 0.0d0
  END IF
END FUNCTION

FUNCTION PlantCover(LandUse)
  REAL(RealKind) :: PlantCover
  INTEGER        :: LandUse
  CALL DataCompute
  PlantCover = PlantCover_Min(LandUse)+(PlantCover_Max(LandUse)-PlantCover_Min(LandUse))  &
                                     *VegFactor*HeightReduct
END FUNCTION

FUNCTION AlbedoSoil(wcontent,SoilClass)
  INTEGER        :: SoilClass
  REAL(RealKind) :: AlbedoSoil 
  REAL(RealKind) :: wcontent
  AlbedoSoil = alb0(SoilClass)-MAX(0.0d0,(wcontent-cres(SoilClass)))/ &
                              (cporv(SoilClass)-cres(SoilClass))*dalb(SoilClass)
END FUNCTION AlbedoSoil

FUNCTION AlbedoVeg(PlantArea,LandClass)
  INTEGER        :: LandClass
  REAL(RealKind) :: AlbedoVeg 
  REAL(RealKind) :: PlantArea
  AlbedoVeg = AlbMin(LandClass)+(AlbMax(LandClass)-AlbMin(LandClass)) &
                                            *VegFactor*HeightReduct 
END FUNCTION AlbedoVeg

FUNCTION Emissivity(wcontent,SoilClass)
  INTEGER        :: SoilClass
  REAL(RealKind) :: Emissivity
  REAL(RealKind) :: wcontent
  Emissivity = emiss0(SoilClass)-MAX(0.0d0,(wcontent-cres(SoilClass)))/ &
                                (cporv(SoilClass)-cres(SoilClass))*demiss(SoilClass)
END FUNCTION Emissivity

END MODULE CanopyData_Mod
