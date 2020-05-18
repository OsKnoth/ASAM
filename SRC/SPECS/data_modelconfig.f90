! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Data module for the configuration of spatial and time discretization
!------------------------------------------------------------------------------

MODULE data_modelconfig

!------------------------------------------------------------------------------
!
! Description:
!  This module contains all variables (scalars and arrays) concerned with the
!  configuration of the horizontal and vertical grid, the time discretization
!  and the reference atmosphere. These are
!    - the vertical coordinate parameters and related variables
!    - horizontal and vertical sizes of the arrays and related variables
!    - start- and end-indices for the computations in the horizontal layers
!    - constants for the horizontal rotated grid and related variables
!    - variables for the time discretization and related variables
!
!  The arrays are declared as allocatable arrays and are allocated in the
!  setup and deallocated in the cleanup of the model. The scalar variables
!  are either part of a NAMELIST (read in the setup) or are initialized
!  in the routine calculate_consts.
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
! 1.5        1998/06/29 Guenther Doms
!  Definition of a new global array hhlr (height of half levels ref. to z=0)
! 1.30       1999/06/24 Matthias Raschendofer
!  Definition of a new INTEGER-parameter kcm used for allocation of canopy fields.
! 1.34       1999/12/10 Ulrich Schaettler
!  Added variables vhmx_vol and vhmx_cfl
! 1.39       2000/05/03 Ulrich Schaettler
!  Renamed variables concerned with latitude and longitude from phi, rla to
!  lat and lon. Included variables for specifying layer index corresponding
!  to certain pressures (has been in data_diagnostics before).
! 2.8        2001/07/06 Ulrich Schaettler
!  Added new variables for multi-layer soil model
! 2.11       2001/09/28 Ulrich Schaettler
!  Renamed a variable for multi-layer soil model
! 2.14       2002/02/15 Ulrich Schaettler
!  New variables for the smooth level vertical coordinate (SLEVE)
! 3.6        2003/12/11 Reinhold Schrodin
!  Added (and changed) some variables for multi-layer soil model
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================
!
! Declarations:
!
USE Kind_Mod

!==============================================================================

IMPLICIT NONE

!==============================================================================

! Global (i.e. public) Declarations:

! 1. vertical coordinate parameters and related variables
! -------------------------------------------------------

  REAL  (RealKind), ALLOCATABLE ::           &
    vcoord(:),    & ! vertical coordinate of LM                    half level
    sigmr (:),    & ! sigma-coordinate referring to PMSL           half level
    hhlr  (:)       ! height-coordinate referring to MSL (z=0)     half level

  REAL  (RealKind)              ::           &
    vcflat,       & ! vertical coordinate where the terrain following system
                    ! changes back to the z-system
    svc1,         & ! vertical decay rate for large-scale topo part of SLEVE
                    ! coordinate (in meter)
    svc2            ! vertical decay rate for small-scale topo part of SLEVE
                    ! coordinate (in meter)

  INTEGER          ::           &
    ivctype,      & ! type of vertical coordinate
                    ! =1 pressure based hybrid coordinate
                    ! =2 height   based hybrid coordinate
                    ! =3 non-specified vertical coordinate
    kflat,        & ! level index where (half) levels become flat
    nfltvc          ! number of filter applications in topo-splitting of
                    ! SLEVE coordinate

! 2. horizontal and vertical sizes of the fields and related variables
! --------------------------------------------------------------------

  INTEGER          ::           &
    ! number of grid points for the total domain
    ie_tot,       & ! number of grid points in zonal direction
    je_tot,       & ! number of grid points in meridional direction
    ke_tot,       & ! number of grid points in vertical direction
    nlandpoints_tot,  & ! number of land points in the grid

    ! number of grid points for this domain
    ie,           & ! number of grid points in zonal direction
    je,           & ! number of grid points in meridional direction
    ke,           & ! number of grid points in vertical direction
    ke_soil,      & ! number of layers in the multi-layer soil model
    nlandpoints,  & ! number of land points in the grid

    ke1,          & ! KE+1
    ieje,         & ! IE*JE
    iejeke,       & ! IE*JE*KE
    ieke,         & ! IE*KE
    ie_max,       & ! Max. of ie on all processors
    je_max,       & ! Max. of je on all processors

    ! start- and end-indices for vertically bounded regions
    kcm             ! index of the lowest model layer, higher than the canopy
 
! 2a. Variables for the new multi-layer soil model
! --------------------------------------------------------------------

  REAL      (RealKind),    ALLOCATABLE       ::           &
    czmls(:)        ! depth of the main soil layers in meters

  INTEGER   , ALLOCATABLE       ::           &
    msoilgrib(:)    ! grib coded depth of main soil levels in centimeters
                    ! (careful: the first level will be coded with 1,
                    !           but is in the depth of 0.5 cm!)

! 3. start- and end-indices for the computations in the horizontal layers
! -----------------------------------------------------------------------
!    These variables give the start- and the end-indices of the 
!    forecast for the prognostic variables in a horizontal layer.
!    Note, that the indices for the wind-speeds u and v differ from 
!    the other ones because of the use of the staggered Arakawa-C-grid.
!    
  INTEGER          ::           &
!   zonal direction
    istart,       & ! start index for the forecast of w, t, qd, qw and pp
    iend,         & ! end index for the forecast of w, t, qd, qw and pp
    istartu,      & ! start index for the forecast of u
    iendu,        & ! end index for the forecast of u
    istartv,      & ! start index for the forecast of v
    iendv,        & ! end index for the forecast of v
    istartpar,    & ! start index for computations in the parallel program
    iendpar,      & ! end index for computations in the parallel program

!   meridional direction
    jstart,       & ! start index for the forecast of w, t, qd, qw and pp
    jend,         & ! end index for the forecast of w, t, qd, qw and pp
    jstartu,      & ! start index for the forecast of u
    jendu,        & ! end index for the forecast of u
    jstartv,      & ! start index for the forecast of v
    jendv,        & ! end index for the forecast of v
    jstartpar,    & ! start index for computations in the parallel program
    jendpar         ! end index for computations in the parallel program

! 4. constants for the horizontal rotated grid and related variables
! ------------------------------------------------------------------

  REAL  (RealKind)              ::           &
    pollon,       & ! longitude of the rotated north pole (in degrees, E>0)
    pollat,       & ! latitude of the rotated north pole (in degrees, N>0)
    dlon,         & ! grid point distance in zonal direction (in degrees)
    dlat,         & ! grid point distance in meridional direction (in degrees)
    startlon_tot, & ! transformed longitude of the lower left grid point
                    ! of the total domain (in degrees, E>0)
    startlat_tot, & ! transformed latitude of the lower left grid point
                    ! of the total domain (in degrees, N>0)
    startlon,     & ! transformed longitude of the lower left grid point
                    ! of this subdomain (in degrees, E>0)
    startlat,     & ! transformed latitude of the lower left grid point
                    ! of this subdomain (in degrees, N>0)
    eddlon,       & ! 1 / dlon
    eddlat,       & ! 1 / dlat
    edadlat,      & ! 1 / (radius of the earth * dlat)
    dlonddlat,    & ! dlon / dlat
    dlatddlon,    & ! dlat / dlon
    degrad,       & ! factor for transforming degree to rad
    raddeg          ! factor for transforming rad to degree

! 5. variables for the time discretization and related variables
! --------------------------------------------------------------

  REAL  (RealKind)              ::           &
    dt,           & ! long time-step
    ed2dt,        & ! 1 / (2 * dt)
    dt2,          & ! 2 * dt
    dtdeh,        & ! dt / 3600 seconds
    epsass,       & ! eps for the Asselin-filter
    betasw,       & ! beta-variable for treatment of soundwaves
    vhmx_vol,     & ! maximum absolute horizontal wind in total model domain
    vhmx_cfl        ! maximum absolute horizontal wind velocity from CFL

  INTEGER          ::           &
    nehddt          ! 3600 seconds / dt

! 6. variables for the reference atmosphere
! -----------------------------------------

  REAL  (RealKind)              ::           &
    p0sl,         & ! reference pressure at sea level
    t0sl,         & ! reference temperature at sea level
    dt0lp           ! d (t0) / d (ln p0)

! 7. Layer index corresponding to a specified pressure
! ----------------------------------------------------

  INTEGER                  ::             &
    klv850,       & ! k index of the LM-mainlevel, on 850 HPa
    klv800,       & ! k index of the LM-mainlevel, on 800 HPa
    klv500,       & ! k index of the LM-mainlevel, on 500 HPa
    klv400,       & ! k index of the LM-mainlevel, on 400 HPa
    klv300          ! k index of the LM-mainlevel, on 300 HPa

!=======================================================================
END MODULE data_modelconfig
