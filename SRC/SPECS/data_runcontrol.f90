! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Data module for variables controlling the run of the model
!------------------------------------------------------------------------------

MODULE data_runcontrol

!------------------------------------------------------------------------------
!
! Description:
!  This module contains variables for running and controlling the forecast.
!  Concerned are the organization of the forecast and the grib I/O.
!  The variables are divided into several groups:
!    
!    - start and end of the forecast
!    - boundary definition and update
!    - controlling the physics
! VG - controlling the microphysics
!    - controlling the dynamics
!    - controlling the nudging
!    - controlling the upper boundary condition
!    - additional control variables
!    - controlling the grib I/O
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
! 1.2        1998/03/30 Ulrich Schaettler
!  Introduction of Namelist variable lcond for switching on/off condensation
! 1.3        1998/04/15 Guenther Doms
!  Introduction of Namelist variabel nincconv (timestep increment for convection)
! 1.4        1998/05/22 Guenther Doms
!  Introduction of Namelist variable l2tls for switching on/off 
!  the two timelevel RK-integration scheme
! 1.5        1998/06/29 Guenther Doms
!  Introduction of Namelist variable rdheight and nrdtaur to control
!  the Rayleigh damping layer at the upper boundary.
! 1.9        1998/09/16 Guenther Doms
!  The namelist variable 'nincmxn' specifying the time averaging intervall
!  for certain output fields has been replaces by 'nincmxt' and 'nincmxu'
!  to allow for different averaging periods.
! 1.10       1998/09/29 Ulrich Schaettler
!  Added new control variables for semi-imp. time stepping, nudging, llm.
! 1.11       1998/10/13 Michael Buchhold
!  Specification of initial fields to be checked for the time range indicator.
! 1.17       1998/11/17 Ulrich Schaettler
!  New control variables for reading and writing ready files
! 1.19       1998/12/11 Christoph Schraff
!  Initial field to be checked corrected from VMO3 to HMO3.
! 1.21       1999/01/25 Guenhter Doms
!  Character variables in list 'yunaman' corrected to new names
! 1.29       1999/05/11 Reinhold Hess
!  Check soil water levels for additional element number
! 1.30       1999/06/24 Matthias Raschendofer
!  Introduction of the time index ntke belonging to the TKE-field.
!  Introduction of 9 INTEGER-namelist-parameters controlling the physics:
!  These are: itype_(wcld, tran, turb, synd), imode_(tran, turb), icldm_(rad, tran, turb).
!  Introduction of 4 LOGICAL-namelist-parameters controlling the physics:
!  These are: lturhor, lexpcor, lnonloc, lcpfluc.
!  Introduction of 3 REAL-namelist-parameters controlling the physics:
!  These are: lam_h, lam_m, pat_len.
!  Introduction of a LOGICAL namelist-parameter controlling the output (lgpspec)
!  Introduction of 3 LOGICAL Namelist Parameters for the convection closure and
!  the soil model :lcape, lctke, lbats
! 1.32       1999/08/24 Guenther Doms
!  New logical control variable 'l2dim' for 2D-runs added.
! 1.33       1999/10/14 Matthias Raschendorfer
!  Introduction of 2 LOGICAL namelist-parameter controlling the physics 
!  (ltmpcor,lprfcor).
!  Removal of a LOGICAL namelist-parameter (lbats).
!  Introduction of 2 INTEGER-namelist-parameters controlling the evaporation:
!  These are: itype_(trvg, evsl).
!  Introduction of a REAL-namelist-parameter (crsmin) to control transpiration.
!  Introduction of 7 REAL-namelist-parameter controlling the turbulence:
!  (tur_len, a_heat, d_heat, a_mom, d_mom, c_diff, rat_lam, rat_can, 
!   c_lnd, c_see)
! 1.34       1999/12/10 Ulrich Schaettler
!  Added variables for unit numbers for output files
! 1.39       2000/05/03 Ulrich Schaettler
!  Removed some organizational data (--> now in organize_diagnosis)
!  Included switch ldfi             (in data_filter before)
!  Introduced switch lw_freeslip    (for treatment of w in nesting)
! 2.2        2000/08/18 Guenther Doms
!  Introduction of the logical switch 'lconf_avg' on Namelist input to
!  enable (default) or disable a horizontal averaging of the convective
!  forcing functions. Also, two new REAL namelist input parameters
!  'c_soil' and 'e_surf' to control surface fluxes have been introduced. 
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new variables for multi-layer soil model and moved some variables to 
!  data_io.f90. Added variable nvers for documenting purposes in Grib-Code
! 2.9        2001/07/16 Guenther Doms
!  Introduction of new global contol parameters (for Namelist input in
!  group INPUT_DYN to control horizontal diffusion: itype_hdiff, hd_corr_t,
!  hd_corr_q and hd_dhmax.
! 2.11       2001/09/28 Ulrich Schaettler
!  Added another variable lmelt_var for multi-layer soil model
! 2.15       2002/03/20 Matthias Raschendorfer
!  Introduction of a REAL namelist-parameter controlling the physics (z0m_dia)
! 2.17       2002/05/08 Ulrich Schaettler
!  New Namelist parameters lkainfri, ltiedtke for choosing a convection scheme
! 3.7        2004/02/18 Ulrich Schaettler
!  New Namelist parameters for computing synthetic satellite images,
!    for turning on/off prognostic precipitation
!    for control variables for 2 tl scheme         and
!    for the ratio of laminar scaling factors for heat over sea and land
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
!==============================================================================
USE Kind_Mod

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. start and end of the forecast
! --------------------------------

  INTEGER          ::           &
    nstart,       & ! first time step of the forecast
    nstop,        & ! last time step of the forecast
    ntstep,       & ! actual time step
                    ! indices for permutation of three time levels
    nold,         & ! corresponds to ntstep - 1
    nnow,         & ! corresponds to ntstep
    nnew,         & ! corresponds to ntstep + 1
                    ! indices for permutation of two time levels
    ntke,         & ! corresponds to ntstep

! VG Time levels in microphysics loop
    nnowMP,       &
    nnewMP

! 2. boundary definition and update
! ---------------------------------

  INTEGER          ::           &
    nlastbound,   & ! time step of the last boundary update
    nincbound,    & ! time step increment of boundary update
    ndiff_ini_bd, & ! difference between start date and date of boundary data
    nbd1,         & ! indices for permutation of the 
    nbd2,         & ! two boundary time levels
    newbc,        & ! number of times that boundary update is analysis after 1 h
    newbcdt,      & ! time step increment of boundary update being derived from
                    ! the (latest) analysis (rather than forecast) fields
    nincboufac      ! factor to 'nincbound' when new boundary update is analysis

! 3. controlling the physics
! --------------------------

  INTEGER          ::           &
    nincrad,      & ! time step increment for running the radiation
    ninctura,     & ! time step increment for running the vertical diffusion
    nincconv,     & ! time step increment for running the convection scheme 

    itype_trvg,   & ! type of vegetation transpiration parameterization
    itype_evsl,   & ! type of parameterization of bare soil evaporation

    itype_gscp,   & ! type of grid-scale precipitaiton physics

    itype_wcld,   & ! type of water cloud diagnosis
    itype_tran,   & ! type of surface-atmosphere transfer
    itype_turb,   & ! type of turbulent diffusion parametrization
    itype_synd,   & ! type of diagnosis of synop. station values
  
    imode_tran,   & ! mode of surface-atmosphere transfer
    imode_turb,   & ! mode of turbulent diffusion parametrization
  
    icldm_rad,    & ! mode of cloud representation in radiation  parametr.
    icldm_tran,   & ! mode of cloud representation in transfer parametr.
    icldm_turb,   & ! mode of cloud representation in turbulence parametr.

    nlgw_ini,     & ! number of prognostic soil water levels in initial data
    nlgw_bd,      & ! number of prognostic soil water levels in boundary data
    nlgw            ! number of prognostic soil water levels



  LOGICAL                          ::           &
    lphys,        & ! forecast with physical parametrizations
    lrad,         & ! forecast with radiation
    ltur,         & ! forecast with vertical diffusion
    lconv,        & ! forecast with convection
    ltiedtke,     & ! Tiedtke mass-flux scheme
    lkainfri,     & ! Kain-Fritsch mass-flux scheme
    lgsp,         & ! forecast with grid scale precipitation
    lprogprec,    & ! forecast with prognostic rain and snow (qr, qs)
    ltrans_prec,  & ! forecast with transport of rain and snow (qr, qs)
    lsoil,        & ! forecast with soil model
    lmelt,        & ! soil model with melting process
    lmelt_var,    & ! freezing temperature dependent on water content
    lmulti_layer, & ! run multi-layer soil model

    lconf_avg,    & ! average convective forcings in case of massflux closure
    lcape,        & ! convection with CAPE closure
    lctke,        & ! convection with turbulent convective energy closure
                    ! warning: lctke not yet fully implemented
    lturhor,      & ! additional horizontal turbulent diffusion
    lexpcor,      & ! explicit corrections of the implicit calculated
                    ! turbulent diffusion (only if itype_turb=3)
    ltmpcor,      & ! consideration of thermal TKE-sources in the 
                    ! enthalpy budget
    lprfcor,      & ! using the profile values of the lowest main level instead
                    ! of the mean value of the lowest layer for surface flux
                    ! calulations
    lnonloc,      & ! nonlocal calculation of vertical gradients used
                    ! for turbulent diffusion (only if itype_turb=3)
    lcpfluc         ! consideration of fluctuations of the heat capacity of air


  REAL (RealKind) ::           &
    crsmin,       & ! minimum value of stomatal resistance 
                    ! (used by the Pen.-Mont. method for vegetation 
                    !  transpiration, itype_trvg=2)
    rlam_mom,     & ! scaling factor of the laminar boudary layer for momentum
    rlam_heat,    & ! scaling factor of the laminar boudary layer for heat
    rat_lam,      & ! ratio of laminar scaling factors for vapour and heat
    rat_can,      & ! ratio of canop[y height over z0m
    rat_sea,      & ! ratio of laminar scaling factors for heat over sea and land
    c_lnd,        & ! surface area density of the roughness elements over land 
                    ! [1/m]
    c_sea,        & ! surface area density of the waves over sea [1/m]
    c_soil,       & ! surface area density of the (evaporative) soil surface
    e_surf,       & ! exponent to get the effictive surface area

    pat_len,      & ! lenth scale of subscale surface patterns over land
                    ! (should be an external parameter field, which is not 
                    !  available yet)
    tur_len,      & ! maximal turbulent length scale
    z0m_dia,      & ! roughness length of a typical synoptic station

    a_heat,       & ! factor for turbulent heat transport
    a_mom,        & ! factor for turbulent momentum transport
    d_heat,       & ! factor for turbulent heat dissipation
    d_mom,        & ! factor for turbulent momentum dissipation
    c_diff          ! factor for turbulent diffusion of TKE


! VG
! 3. controlling the 2D Surface
! --------------------------
    REAL (RealKind)  ::      &
         hill_halfwidth,        &  ! [m]
         hill_height               ! [m]


    INTEGER  :: &
         topo_type                 ! only value 2 at the moment! see src_artifdata
 
! 3.1 controlling the microphysics
! --------------------------------

  INTEGER          ::           &
       n_class,       &  ! Number of mass classes
       iap,           &  ! Initial spectrum of aerosols, chose 1,2,3 (default: 3)
       itimemax          ! number of microphysical time step per LM time step
                         ! take care that dt/itimemax = INTEGER!

  LOGICAL                          ::           &

       lmp,           &  ! switch on (.TRUE.) / switch off (.FALSE.) microphysics, 
                         ! default is .FALSE.
       lice              ! switch on/off microphysical ice scheme
  
  REAL (RealKind) ::           &
       dnap_init         ! Initial number of Aerosols in unimodal distribution, 1/cm3



! 4. controlling the dynamics
! ---------------------------

  LOGICAL                          ::           &
    l2tls,        & ! time integration by two timelevel RK-scheme (.TRUE.)
                    ! or by default three-time level KW-scheme (.FALSE.)
    lvertad_impl, & ! if =.TRUE.:  implicit vertical advection in RK-scheme
                    !    =.FALSE.: explicit vertical advection ...
    lsemi_imp       ! if .TRUE.,  running with semi-implicit scheme,
                    ! else with split-explicit scheme (only for l2tls=FALSE!)

  INTEGER          ::           &
    irunge_kutta, & ! =0: use scheme from module src_2timelevel,
                    ! =1: use new RK scheme from module src_runge_kutta,
                    ! =2: use new TVD-RK scheme
    irk_order,    & ! order of the Runge-Kutta scheme
    iadv_order,   & ! order of the advection scheme
    ikrylow_si,   & ! dimension of the Krylow space used in the elliptic
                    ! solver for the semi-implicit scheme
    maxit_si,     & ! maximum number of iterations for the elliptic solver
    iprint_si       ! to control whether statistics of the solver are printed

  REAL      (RealKind)          ::           &
    xkd,          & ! coefficient for divergence damping
    eps_si          ! precision limit for the elliptic solver

! 5. controlling the observation processing
! -----------------------------------------

  LOGICAL                          ::           &
    luseobs         ! on - off switch for using observational data for:
                    ! - nudging (of conventional data)
                    ! - latent heat nudging (not implemented yet)
                    ! - 2-dim. analyses (2m-Temperature, 2m-Humidity, precipit.)
                    ! - verification of model data against observations

! 6. controlling the upper boundary condition
! -------------------------------------------

  LOGICAL                          ::           &
    lspubc,       & ! with Rayleigh damping in the upper levels
    lrubc           ! radiative upper boundary condition

  REAL      (RealKind)          ::           &
    rdheight        ! bottom height of Rayleigh damping layer

  INTEGER          ::           &
    nrdtau          ! number of time steps in Rayleigh damping time scale     

! 7. additional control variables
! -------------------------------

  LOGICAL                          ::           &
    llm,          & ! if .TRUE., running with a lowered upper boundary
    lprog_qi,     & ! if .TRUE., running with cloud ice
                    !   (this is set internally by the program,
                    !    depending on itype_gscp)
    lcond,        & ! forecast with condensation/evaporation
    lclock,       & ! system clock is present
    ltime,        & ! detailled timings of the program are given
    ltime_mean,   & ! mean values of the times of all processors are given
    ltime_proc,   & ! the timings are given for every single processors
    lreproduce,   & ! the results are reproducible in parallel mode
    lhordiff,     & ! running with horizontal diffusion
    ldebug,       & ! with control output on a special file 
    lrerun,       & ! not implemented
    lrout,        & ! routine-forecast of the model
    lgen,         & ! forecast with self-defined artificial data
    lperi,        & ! if lgen=.TRUE.: periodic boundary conditions (.TRUE.)
                    !                 or with Davies conditions (.FALSE.)
    l2dim,        & ! if lgen=.TRUE.: 2dimensional model version (.TRUE) or
                    !                 full 3dimensional version (.FALSE)
    lcori,        & ! if lgen=.TRUE.: with Coriolis force (.TRUE.)
                    !                 or without Coriolis force (.FALSE.)
    lmetr,        & ! if lgen=.TRUE.: with metric terms (.TRUE.)
                    !                 or without metric terms (.FALSE.)
    lw_freeslip,  & ! if .TRUE.: with free slip lateral boundary condition and
                    ! if .FALSE. specified lateral boundary values for w
    ldfi,         & ! switch for initialization by digital filtering
                    ! (if .TRUE. do apply dfi)
    lcon_clw,     & ! if .TRUE.: convective liquid water used in rttov
    luse_rttov      ! if rttov-library is used

  INTEGER          ::           &
    itothours,    & ! number of forecast hours
    nincmxt,      & ! time step increment for deleting tmin, tmax
    nincmxu,      & ! time step increment for deleting vbmax
    nvers,        & ! version number of experiment for documentation
    itype_hdiff     ! type of horizontal diffusion (=1: 4th order linear),
                    ! =2: 4th order linear monotonic with orographic limiter)

  REAL      (RealKind)          ::           &
    hd_corr_u,    & ! correction factor for horizontal diffusion flux of u,v,w
    hd_corr_t,    & ! correction factor for horizontal diffusion flux of t,p
    hd_corr_q,    & ! correction factor for horizontal diffusion fluc of qv,qc
    hd_dhmax        ! maximum gridpoint height difference for applying
                    ! horizontal diffusion fluxes between them

! 8. diagnostic calculations
! --------------------------

  LOGICAL                          ::           &
    ldiagnos        ! perform diagnostic calculations

! 9. Variables for Ascii file handling, time measuring, ...
! ---------------------------------------------------------

  INTEGER          ::           &
    ihours          ! actual hour of the forecast

  REAL (RealKind), ALLOCATABLE  ::           &
    timings (:,:)   ! for storing the times for different parts of the program

  CHARACTER (LEN=10)               ::           &
    yakdat1   ! actual date (ydate_ini+ntstep/dt) in the form 
              ! ddmmyyyyhh (day, month, year, hour)
  CHARACTER (LEN=22)               ::           &
    yakdat2   ! actual date (ydate_ini+ntstep/dt) in the form 
              ! wd dd.mm.yyyy  hh UTC  (weekday, ...)
                                                
  INTEGER          ::           &
    nusolver,     & ! unit number for file YUSOLVER
    nudebug,      & ! unit number for file YUDEBUG
    nuspecif        ! unit number for file YUSPECIF

  CHARACTER (LEN= 8) :: yusolver='YUSOLVER'
  CHARACTER (LEN= 7) :: yudebug ='YUDEBUG'
  CHARACTER (LEN= 8) :: yuspecif='YUSPECIF'

!==============================================================================

END MODULE data_runcontrol
