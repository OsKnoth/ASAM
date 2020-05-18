! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Data module for all global meteorological fields
!------------------------------------------------------------------------------

MODULE data_fields

!------------------------------------------------------------------------------
!
! Description:
!  This module declares all meteorological fields that have to reside in 
!  the long term storage, i.e. that are used in more than one module.
!  Fields included are
!    - constant fields defining the reference atmosphere
!    - external parameter fields
!    - prognostic variables
!    - tendency fields for the prognostic variables
!    - fields for surface values
!    - fields that are computed in the parametrization packages 
!      or in the dynamics
!    - fields for model-output and diagnostics
!    - fields for the boundary values
!
!  All fields are declared as allocatable arrays. They are allocated in the
!  setup of the model and deallocated in the cleanup at the end of the
!  program.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8062 3721
!  email:  ulrich.schaettler@dwd.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.3        1998/04/15 Guenther Doms
!  Definition of new tendency arrays for convection
! 1.7        1998/07/16 Guenther Doms
!  Removal of global array 'rrssk'.
! 1.20       1999/01/07 Guenhter Doms
!  Renaming of some global variables.
! 1.24       1999/03/01 Guenther Doms
!  Declaration of a new 3-D prognostic variable (qi).
! 1.30       1999/06/24 Matthias Raschendorfer
!  Declaration of a 3-D prognostic var. (tke) and its tendenci field (tketens).
!  Declaration of a 3-D array (rcld).
!  Declaration of 5 3-D arrays for canopy layers (c_big,c_sml, r_air, t_e,qv_e).
!  Declaration of a 2-D arrays (tfh,tfm) and (h_can,d_pat).
!  Declaration of 4 2-D variables for new convection closures
! 1.33       1999/10/14 Reinhold Hess
!  Declaration of new global 2-D arrays idiv_hum and aevap_s for diagnosis
!  of the model water budget
!  Declaration of a new 2-D array 'sai' for surface area index (M.Raschendorfer)
! 1.34       1999/12/10 Ulrich Schaettler
!  Named all boundary fields with "_bd" (for consistency)
! 1.39       2000/05/03 Ulrich Schaettler
!  Add declaration of a boundary field for w (used for interactive nesting)
! 2.2        2000/08/18 Matthias Raschendorfer
!  Declaration of the 2-D arrays 'eai' and 'tai'.
! 2.4        2001/01/29 Christoph Schraff
!  Declaration of the 2-D array 'prne_con', for humidity balancing at T-nudging.
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new fields for multi-layer soil model and surface fluxes
! 2.11       2001/09/28 Ulrich Schaettler
!  Added new fields for lateral values of cloud ice
! 2.17       2002/05/08 Ulrich Schaettler
!  New fields for Kain-Fritsch convection scheme
! 2.18       2002/07/16 Ulrich Schaettler
!  New fields for specific rain and snow content;
!  included declaration of a1t, a2t from src_leapfrog
! 3.5        2003/09/02 Ulrich Schaettler
!  New fields phi_tot, rla_tot to avoid global communication in the radiation
! 3.6        2003/12/11 Reinhold Schrodin
!  New field freshsnow for new multi-layer soil model added
! 3.7        2004/02/18 Ulrich Schaettler
!  New fields for computing synthetic satellite images (synme5-7, synmsg)
!  New field for storing convective cloud water (clw_con)
!  Renamed alb (alb_rad), idiv_hum (tdiv_hum), phi (rlat), rla (rlon),
!      cphi (crlat), acphir (acrlat), tgphi (tgrlat) (for consistency with GME
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

! VG this routine is modified by Verena Gruetzun, gruetzun@tropos.de

!
! Declarations:
!
USE Kind_Mod
!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. constant fields for the reference atmosphere                     (unit)
! -----------------------------------------------

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
    rho0(:,:,:),    & ! reference density at the full model levels    (kg/m3)
    dp0 (:,:,:),    & ! reference pressure thickness of layers        ( Pa  )
    p0  (:,:,:),    & ! reference pressure at full levels             ( Pa  )
    hhl (:,:,:)       ! geometrical height of half levels             ( m   )

! 2. external parameter fields                                        (unit)
! ----------------------------

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
    hsurf  (:,:),   & ! height of surface topography                  ( m   )
    gz0    (:,:),   & ! surface roughness * g                         (m2/s2)
    fr_land(:,:),   & ! fraction of land in a grid element              --
    soiltyp(:,:),   & ! type of the soil (keys 0-9)                     --
    vio3   (:,:),   & ! vertical integrated ozone contents            (Pa O3)
    hmo3   (:,:),   & ! ozone maximum                                 ( Pa  )
    rlat   (:,:),   & ! geographical latitude                         ( rad )
    rlon   (:,:),   & ! geographical longitude                        ( rad )
    rlattot(:,:),   & ! geographical latitude                         ( rad )
    rlontot(:,:),   & ! geographical longitude                        ( rad )
    fc     (:,:),   & ! coriolis-parameter                            ( 1/s )
    rmy    (:,:,:), & ! Davis-parameter for boundary relaxation         --
    crlat  (:,:),   & ! cosine of transformed latitude                  --
    acrlat (:,:),   & ! 1 / ( crlat * radius of the earth )           ( 1/m )
    tgrlat (:)  ,   & ! tangens of transformed latitude                 --
    aerlan (:,:),   & ! aerosol-distribution on rural areas             --
    aerurb (:,:),   & ! aerosol-distribution on urban areas             --
    aerdes (:,:),   & ! aerosol-distribution on desert areas            --
    aersea (:,:),   & ! aerosol-distribution on the sea                 --
    plcov  (:,:),   & ! fraction of plant cover                         --
    lai    (:,:),   & ! leaf area index of plants                       --
    tai    (:,:),   & ! transpiration area index                        --
    sai    (:,:),   & ! surface area index                              --
    eai    (:,:),   & ! (evaporative) earth area index                  --
    rootdp (:,:),   & ! depth of the roots                            ( m  )

    h_can (:,:),    & ! hight of the vertically resolved canopy       ( m )
    d_pat (:,:),    & ! horizontal pattern length scale               ( m )
    c_big (:,:,:),  & ! effective drag coefficient of canopy elements
                      ! larger than or equal to the tubulent length
                      ! scale                                         (1/m)
    c_sml (:,:,:),  & ! effective drag coefficient of canopy elements
                      ! smaller than the tubulent length scale        (1/m)
    r_air (:,:,:)     ! air containing fraction of a gridbox inside
                      ! the canopy                                    ( 1 )


  LOGICAL, ALLOCATABLE ::           &
    llandmask(:,:)    ! landpoint mask

! 3. prognostic variables                                             (unit)
! -----------------------

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
    u (:,:,:,:),    & ! zonal wind speed                              ( m/s )
    v (:,:,:,:),    & ! meridional wind speed                         ( m/s )
    w (:,:,:,:),    & ! vertical wind speed (defined on half levels)  ( m/s )
    t (:,:,:,:),    & ! temperature                                   (  k  )
    qv(:,:,:,:),    & ! specific water vapor content                  (kg/kg)
    qc(:,:,:,:),    & ! specific cloud water content                  (kg/kg)
    qi(:,:,:,:),    & ! specific cloud ice content                    (kg/kg)
    qr(:,:,:,:),    & ! specific rain content                         (kg/kg)
    qs(:,:,:,:),    & ! specific snow content                         (kg/kg)
    pp(:,:,:,:),    & ! deviation from the reference pressure         ( pa  )
    tke(:,:,:,:)      ! SQRT(2 * turbulent kinetik energy)           ( m/s )
                      ! (defined on half levels)


! 4. tendency fields for the prognostic variables                     (unit )
! -----------------------------------------------
!    time tendencies  by diabatic and adiabatic processes
!    without sound-wave terms

  REAL  (RealKind), ALLOCATABLE ::           &
    utens (:,:,:),  & ! u-tendency without sound-wave terms           ( m/s2)
    vtens (:,:,:),  & ! v-tendency without sound-wave terms           ( m/s2)
    wtens (:,:,:),  & ! w-tendency without sound-wave terms           ( m/s2)
                      ! (defined on half levels)
    ttens (:,:,:),  & ! t-tendency without sound-wave terms           ( K/s )
    qvtens(:,:,:),  & ! qv-tendency                                   ( 1/s )
    qctens(:,:,:),  & ! qc-tendency                                   ( 1/s )
    pptens(:,:,:),  & ! pp-tendency without sound-wave terms          (Pa/s )

    tketens(:,:,:), & ! tke-tendency (defined on half-levels)          ( m/s )
! VG
! VG start
    qctens_c(:,:,:,:)
! VG
 
! 5. fields for surface values and soil/canopy model variables        (unit )
! -----------------------------------------------------

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
    ps (:,:,:),     & ! surface pressure                              ( pa  )
    t_snow(:,:,:),  & ! temperature of the snow-surface               (  k  )
    t_s (:,:,:),    & ! temperature of the ground surface (soil)      (  k  )
    t_g (:,:,:),    & ! weighted surface temperature                  (  k  )
    qv_s(:,:,:),    & ! specific water vapor content at the surface   (kg/kg)
    t_m (:,:,:),    & ! temperature between upper and medium 
                      ! soil layer                                    (  k  )
    t_cl(:,:),      & ! temperature between medium and lower
                      ! soil layer (climatology)                      (  k  )
    t_so(:,:,:,:),  & ! multi-layer soil temperature                  (  k  )
    w_snow(:,:,:),  & ! water content of snow                         (m H2O)
    w_i (:,:,:),    & ! water content of interception water           (m H2O)
    w_g1(:,:,:),    & ! water content of the upper soil layer         (m H2O)
    w_g2(:,:,:),    & ! water content of the medium soil layer        (m H2O)
    w_g3(:,:,:),    & ! water content of the lower soil layer         (m H2O)
                      ! (if nlgw=3, unused otherwise)
    w_so(:,:,:,:),  & ! multi-layer soil moisture                     (m H2O)
    w_ice(:,:,:,:), & ! multi-layer soil ice                          (m H2O)
    w_cl(:,:),      & ! climatological water content                  (m H2O) 
    freshsnow(:,:), & ! weighting function indicating 'freshness' of snow
    t_e(:,:,:),     & ! surface temperature of the canopy elements    (  k  )
    qv_e(:,:,:)       ! surface value of qv of the canopy elements    (Kg/Kg)
 

! 6. fields that are computed in the parametrization and dynamics     (unit )
! ---------------------------------------------------------------

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &

!   air density at present time level (main levels)       
    rho (:,:,:),    & ! total density of air                          (kg/m3)

!   coefficients for turbulent diffusion in the atmosphere
!   (defined on half levels)
    tkvm(:,:,:),    & ! turbulent diffusion coefficient for momentum  (m2/s)
    tkvh(:,:,:),    & ! turbulent diffusion coefficient for heat      (m2/s)
                      ! and moisture

!   vertical varying implicitness of vertical diffusion
!   (def. at half levels)
    a1t(:),         & !                                               ( -- )
    a2t(:),         & !                                               ( -- )

!   turbulence statistics in the atmosphere
!   (defined on full levels)
    rcld(:,:,:),    & ! standard deviation of the saturation deficit    --

!   turbulent coefficients at the surface
    tcm (:,:),      & ! transfer coefficient for momentum             ( -- )
    tch (:,:),      & ! transfer coefficient for heat and moisture    ( -- )

    tfm (:,:),      & ! factor of laminar transfer of momentum           --
    tfh (:,:),      & ! factor of laminar transfer of scalars            --
 
!   fields from the radiation scheme
    sohr (:,:,:),   & ! rate of solar heating                         ( k/s )
    thhr (:,:,:),   & ! rate of thermal heating                       ( k/s )
    clc_sgs(:,:,:), & ! subgrid-scale stratiform cloud cover            --
    alb_rad(:,:),   & ! albedo of the ground                            --
    sobs (:,:),     & ! solar radiation at the ground                 ( w/m2)
    thbs (:,:),     & ! thermal radiation at the ground               ( w/m2)
    pabs (:,:),     & ! photosynthetic active radiation at the ground ( w/m2)
    sobt (:,:),     & ! solar radiation at the upper boundary         ( w/m2)
                      ! of the atmosphere
    thbt (:,:),     & ! thermal radiation at the upper boundary       ( w/m2)
                      ! of the atmosphere
    clch (:,:),     & ! cloud cover with high clouds                    --   
    clcm (:,:),     & ! cloud cover with medium clouds                  --   
    clcl (:,:),     & ! cloud cover with low clouds                     --   
    clct (:,:)        ! total cloud cover                               --   

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
!   fields from the convection scheme
    clc_con(:,:,:), & ! cloud cover due to convection                   --     
    clw_con(:,:,:), & ! convective cloud liquid water
    prr_con(:,:),   & ! precipitation rate of rain, convective        (kg/m2*s)
    prs_con(:,:),   & ! precipitation rate of snow, convective        (kg/m2*s)
    prne_con(:,:),  & ! precipitation rate, no evaporat., convective  (kg/m2*s)
    bas_con(:,:),   & ! level index of convective cloud base            --
    top_con(:,:),   & ! level index of convective cloud top             --
    tt_conv (:,:,:),& ! temperature tendency due to convection        ( K/s  )
    qvt_conv(:,:,:),& ! humidity    tendency due to convection        ( 1/s  )
    qct_conv(:,:,:),& ! qc-tendency tendency due to convection        ( 1/s  )
    qit_conv(:,:,:),& ! qi-tendency tendency due to convection        ( 1/s  )
    qrt_conv(:,:,:),& ! qr-tendency tendency due to convection        ( 1/s  )
    qst_conv(:,:,:),& ! qs-tendency tendency due to convection        ( 1/s  )
    ut_conv (:,:,:),& ! u-tendency due to convection                  ( m/s^2)
    vt_conv (:,:,:),& ! v-tendency due to convection                  ( m/s^2)
    mflx_con(:,:),  & ! cloud base massflux                           (kg/m2*s)
    cape_con(:,:),  & ! convective available energy                   (   J/kg)
    tke_con (:,:),  & ! convective turbulent energy                   (   J/kg)
    qcvg_con(:,:),  & ! moisture convergence for Kuo-type closure     (    1/s)
    w0avg   (:,:,:),& ! running average of w
! VG
! VG start
    qct_conv_c(:,:,:,:)
! VG end

  INTEGER , TARGET, ALLOCATABLE ::           &
    nca     (:,:)     !

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
!   fields from the grid-scale precipitation scheme
    qrs    (:,:,:), & ! precipitation water content (water loading)   (kg/kg)
    prr_gsp(:,:),   & ! precipitation rate of rain, grid-scale        (kg/m2*s)
    prs_gsp(:,:),   & ! precipitation rate of snow, grid-scale        (kg/m2*s)

!   fields that are computed in the dynamics
    dqvdt  (:,:,:), & ! threedimensional moisture convergence         ( 1/s )
    qvsflx (:,:),   & ! surface flux of water vapour                  ( 1/m2s)
    dpsdt  (:,:),   & ! tendency of the surface pressure              ( pa/s)
    umfl_s (:,:),   & ! u-momentum flux (surface)                     ( N/m2)
    vmfl_s (:,:),   & ! v-momentum flux (surface)                     ( N/m2)
    shfl_s (:,:),   & ! sensible heat flux (surface)                  ( W/m2)
    lhfl_s (:,:),   & ! latent heat flux (surface)                    ( W/m2)
    aumfl_s(:,:),   & ! average u-momentum flux (surface)             ( N/m2)
    avmfl_s(:,:),   & ! average v-momentum flux (surface)             ( N/m2)
    ashfl_s(:,:),   & ! average sensible heat flux (surface)          ( W/m2)
    alhfl_s(:,:)      ! average latent heat flux (surface)            ( W/m2)

! 7. fields for model output and diagnostics                          (unit )
! ------------------------------------------

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
    t_2m    (:,:),  & ! temperature in 2m                             (  K  )
    qv_2m   (:,:),  & ! specific water vapor content in 2m            (kg/kg)
    td_2m   (:,:),  & ! dew-point in 2m                               (  K  )
    u_10m   (:,:),  & ! zonal wind in 10m                             ( m/s )
    v_10m   (:,:),  & ! meridional wind in 10m                        ( m/s )
    tmin_2m (:,:),  & ! minimum temperature in 2m                     (  K  )
    tmax_2m (:,:),  & ! maximum temperature in 2m                     (  K  )
    vmax_10m(:,:),  & ! maximal windspeed in 10m                      ( m/s )
    asob_s  (:,:),  & ! average solar radiation budget (surface)      ( W/m2)
    athb_s  (:,:),  & ! average thermal radiation budget (surface)    ( W/m2)
    apab_s  (:,:),  & ! average photosynthetic active radiation (sfc) ( W/m2)
    asob_t  (:,:),  & ! average solar radiation budget (model top)    ( W/m2) 
    athb_t  (:,:),  & ! average thermal radiation budget (model top)  ( W/m2)
    rain_gsp(:,:),  & ! amount of rain from grid-scale precip. (sum)  (kg/m2)
    snow_gsp(:,:),  & ! amount of snow from grid-scale precip. (sum)  (kg/m2)
    rain_con(:,:),  & ! amount of rain from convective precip. (sum)  (kg/m2)
    snow_con(:,:),  & ! amount of snow from convective precip. (sum)  (kg/m2)
    runoff_s(:,:),  & ! surface water runoff; sum over forecast       (kg/m2)
    runoff_g(:,:),  & ! soil water runoff; sum over forecast          (kg/m2)
    tdiv_hum(:,:),  & ! vertical sum for  divergence of humidity      (kg/m2)
    aevap_s (:,:)     ! accumulated surface moisture flux             (kg/m2)

! 8. fields for the boundary values                                   (unit )
! ---------------------------------

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
    u_bd   (:,:,:,:), & ! boundary field for u                        ( m/s )
    v_bd   (:,:,:,:), & ! boundary field for v                        ( m/s )
    w_bd   (:,:,:,:), & ! boundary field for w                        ( m/s )
    t_bd   (:,:,:,:), & ! boundary field for t                        (  k  )
    qv_bd  (:,:,:,:), & ! boundary field for qv                       (kg/kg)
    qc_bd  (:,:,:,:), & ! boundary field for qc                       (kg/kg)
    qi_bd  (:,:,:,:), & ! boundary field for qi                       (kg/kg)
    pp_bd  (:,:,:,:), & ! boundary field for pp                       (  pa )
    qv_s_bd  (:,:,:), & ! boundary field for qv_s                     (kg/kg)
    t_snow_bd(:,:,:), & ! boundary field for t_snow                   (  k  )
    t_s_bd   (:,:,:), & ! boundary field for t_s                      (  k  )
    t_m_bd   (:,:,:), & ! boundary field for t_m                      (  k  )
    w_snow_bd(:,:,:), & ! boundary field for w_snow                   (m H2O)
    w_g1_bd  (:,:,:), & ! boundary field for w_g1                     (m H2O)
    w_g2_bd  (:,:,:), & ! boundary field for wb2                      (m H2O)
    w_g3_bd  (:,:,:), & ! boundary field for wb3                      (m H2O)
! VG
! VG start
    qc_bd_c(:,:,:,:,:)

! 9. fields for the synthetic satellite images
! --------------------------------------------

  REAL  (RealKind), TARGET, ALLOCATABLE ::           &
    synme5 (:,:,:)  , & ! Meteosat 5
    synme6 (:,:,:)  , & ! Meteosat 6
    synme7 (:,:,:)  , & ! Meteosat 7
    synmsg (:,:,:)      ! Meteosat Second Generation

!==============================================================================

END MODULE data_fields
