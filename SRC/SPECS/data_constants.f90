! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Data module for mathematical, physical and parametrizational constants
!-------------------------------------------------------------------------------

MODULE data_constants

!-------------------------------------------------------------------------------
!
! Description:
!  This module contains mathematical and physical constants as well as 
!  constants for the parametrizations and related variables.
!
! Current Code Owner: DWD, Ulrich Schaettler
!  phone:  +49  69  8062 2739
!  fax:    +49  69  8236 1493
!  email:  uschaettler@dwd.d400.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.1        1998/03/11 Ulrich Schaettler
!  Initial release
! 1.20       1999/01/07 Guenther Doms
!  Renaming of some global variables.
! 1.30       1999/06/24 Matthias Raschendorfer
!  Declaration of 5 new constants (lhocp, rcpv, rcpl, con_m, con_h)
! 2.18       2002/07/16 Reinhold Schrodin
!  Eliminated variable rhde (will be replaced by cf_snow from data_soil)
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!=======================================================================
!
! Declarations:
!

USE Kind_Mod

!=======================================================================

IMPLICIT NONE

!=======================================================================
!
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
! Global (i.e. public) Declarations:

! 1. mathematical constants
! -------------------------

  REAL  (RealKind)              ::           &
    pi              ! circle constant

! 2. physical constants and related variables
! -------------------------------------------

  REAL  (RealKind)              ::           &
    t0,           & ! melting temperature of ice
    r_d,          & ! gas constant for dry air
    r_v,          & ! gas constant for water vapor
    rdv,          & ! r_d / r_v
    o_m_rdv,      & ! 1 - r_d/r_v
    rvd_m_o,      & ! r_v/r_d - 1
    cp_d,         & ! specific heat of dry air at constant pressure
    cpdr,         & ! 1 / cp_d
    rdocp,        & ! r_d / cp_d
    gamma,        & ! 1 / (1 - rdocp)   ( = cp_d/cv_d)
    lh_v,         & ! latent heat of vapourization
    lh_f,         & ! latent heat of fusion
    lh_s,         & ! latent heat of sublimation
    lhocp,        & ! lh_v / cp_d
    rcpv,         & ! cp_d / cp_v - 1
    rcpl,         & ! cp_d / cp_l - 1
    con_m,        & ! kinematic viscosity (m2/s)
    con_h,        & ! scalar conductivity (m2/s)
    g,            & ! acceleration due to gravity
    gq,           & ! g * g
    gh,           & ! g / 2
    gr,           & ! 1 / g
    r_earth,      & ! mean radius of the earth (m)
    day_len,      & ! mean length of the day (s)
    rho_w,        & ! density of liquid water (kg/m^3)
    sigma,        & ! Boltzmann-constant
    solc            ! solar constant

! 3. constants for parametrizations
! ---------------------------------

  REAL  (RealKind)              ::           &
    b1,           & ! variables for computing the saturation vapour pressure
    b2w,          & ! over water (w) and ice (i) according to Teten's formula
    b2i,          & !               -- " --
    b3,           & !               -- " --
    b4w,          & !               -- " --
    b4i,          & !               -- " --
    b234w,        & ! b2w * (b3 - b4w)
    uc1,          & ! variable for computing the rate of cloud cover in 
    uc2,          & ! the unsaturated case
    ucl,          & !               -- " --
    aks2,         & ! variable for horizontal diffusion of second order
    aks4,         & ! variable for horizontal diffusion of fourth order
    akt             ! von Karman-constant

!=======================================================================
END MODULE data_constants
