! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Variables for the computation of synthetic satellite images
!------------------------------------------------------------------------------

MODULE data_mp_ice

!------------------------------------------------------------------------------
!
! Description:
!  This data module contains all data necessary for the computation of 
!  the microphysical processes. It replaces the common blocks in the 
!  original microphysics routines written in fixed format
!
! Current Code Owner: IfT, Verena Gruetzun
!  phone:  +49  341 235 2228
!  fax:    +49  341 235 2139
!  email:  gruetzun@tropos.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----
! 1.0   DECEMBER 2005    Verena Gruetzun
!  Initial release
! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".

!
! Modules used:
!               
USE Kind_Mod

USE data_modelconfig,   ONLY:    &
     ie,je,ke        ! number of gridpoints in zonal/meridional direction




!-------------------------------------------------------------------------------
!---------------
!Modellparameter
!---------------
!jmax: Anzahl der Tropfenklassen
!nmax: Aufloesung des Tropfenspektrums
!smax: Anzahl der Aerosolklassen
!simax: Anzahl der unloeslichen Aerosolklassen
!tmax: Aufloesung des unloeslichen Aerosolspektrums
!itmax: Anzahl der extern gemischten Aerosoltypen
!ipmax: horizontale Aufloesung
!-----------------------------------------------------
!| jmax | nmax | smax | simax | tmax | itmax | ipmax |
!-----------------------------------------------------
!   66      2      1      66      2      1       1

!==============================================================================

IMPLICIT NONE

!==============================================================================
!   
! Declarations must be of the form:
! <type> <VariableName> ! Description/ purpose of variable
! Declarations:

!==============================================================================


! From AKM

REAL  (RealKind), ALLOCATABLE     ::     &
     ini    (:,:),  & ! Anfangswerte, werden ueberall wieder mit reingegeben
     vt     (:,:),  & ! wahrscheinlich nutzls, frueher vbest, 
                      ! jetzt nicht mehr verwendet (lt. Harald)
     miv    (:  ),  & ! Kontaktwinkel, wird aus datei eingelesen
     rrnw   (: ),   & ! Niederschlag sammeln
     rrqw   (: ),   & ! "
     rrnf   (: ),   & ! "
     rrqf   (: ),   &   ! "
     mquer  (:,:,:), &
     rquer  (:,:,:), &
     squer  (:,:,:), &
     mfquer (:,:,:), &
     rfquer (:,:,:), &
     sfquer (:,:,:)


! COMMON PHYSICS

REAL  (RealKind)                  ::     &
     CP,         &  ! cp_d from LM
     GRAV,       &  ! g from LM
     KTH,        &  !
     LS,         &  !
     LV,         &  ! lh_v from LM
     LV_ICE,     &  !
     mol_w,      &  ! molality of water (not used?)
     PI,         &  ! pi from LM
     RDRY,       &  !
     RW,         &  ! gas constant dry air 
     RHOW,       &  ! rho_w from LM
     RHOI,       &  !
     SIGMA0,     &  ! surface tension of water (not used?)
     mol_AS,     &  ! moles of soluble aerosol mass
     mol_AI,     &  ! moles of insoluble aerosol mass
     rho_AS,     &  ! density of soluble aerosol mass
     rho_AI,     &  ! density of insoluble aerosol mass
     vantS,      &  ! (Stoff konstante)
     E0,         &  ! saturation ratio at ... (not used?)
     C1,         &  ! constant used to calculate saturation ratio 
     C2,         &  ! "
     mol_g(7),   &  ! piecewise linear function used for osmotic coefficient
     phi_g(7),   &  ! "
     a_phi(7),   &  ! "
     b_phi(7),   &  ! "
     CAPPA,      &  ! RDRY/cp
     FAC34,      &  ! 3/(4*pi)
     FACT,       &  ! 4/3 * pi rho_w
     FACTI,      &  !
     fact_AS,    &  ! 4/3 * pi rho_AS
     fact_AI,    &  ! 4/3 * pi rho_AI
     GAMMA,      &  !
     RAD_BUB0,   &  !
     k1_vel,     &  !
     k2_vel,     &  !
     rho_0          ! 

INTEGER , ALLOCATABLE  ::     &
     ifreeze(: )      ! freezing par, Bakterien, Russ etc. pp



! COMMON MIKROPH

INTEGER              ::     &
     ITERMAX        ! maximal number of iteration

REAL  (RealKind)                  ::     &
     R0,         &  ! variates class boundaries (see sconst) 
     RAmin,      &  !
     RAmax,      &  !
     PR,         &  ! constant in apvert (??!)
     TC,         &  ! used for temperature dependance of water surface tension
     mu,         &  !
     b1,         &  ! used for temperature dependance of water surface tension
     b2,         &  ! "
     bb,         &  ! "
     ALPHA_C,    &  ! constant in diffusion and heat conduction
     ALPHA_T,    &  ! "
     DEL_V,      &  ! ? size dependant pars, DIFFSTAR, KTHSTAR
     DEL_T,      &  ! "
     WEI_Q,      &  ! weighing factor for bin shift
     WEI_N,      &  ! "
     SC_cond,    &  ! ? boundaries for saturation ratio
     SC_ggw,     &  !
     SATTabs,    &  ! (not used?)
     size_fact,  &  ! stop criteria for iteration loop
     accu_fact,  &  ! "
     step_fact      ! "


! COMMON EMPIRIE


! modify eventually to be defined locally in break - do a grep
REAL  (RealKind)                  ::     &
     PC1,        &  !
     PC2,        &  !
     QC1,        &  !
     QC2            !


! COMMON NUMERIK

INTEGER              ::     &
     J_D01,      &  !
     J_D05,      &  !
     J_D10          !

REAL  (RealKind)                  ::     &
     DELTAT,     &  ! internal time strp
     TIMMAX,     &  !
     INCOUT,     &  !
     INCHEI,     &  !
     SMALL1,     &  ! small constant
     QU1D3,      &  ! 1/3
     QU4D3,      &  !
     QU5D3          !


! COMMON SPEKGRID

REAL  (RealKind), ALLOCATABLE      ::     &               
     DIFF21(:),  & ! difference of two boundaries
     MGRENZ(:),  & ! boundaries of classes - mass 
     MMITTE(:),  & ! center of classes     - mass 
     RGRENZ(:),  & ! boundaries of classes - radius
     RMITTE(:),  & ! center of classes     - radius
     SDIFF21(:), & !
     SGRENZ(:),  & !
     SMITTE(:),  & !
     RSGRENZ(:), & !
     RSMITTE(:)    !


! COMMON KOLLISION

INTEGER  , ALLOCATABLE ::     &
     KZIEL(:,:)    ! target bin for collision

REAL  (RealKind)      , ALLOCATABLE ::     &
     KOLK(:,:), &  ! collision kernel
     KOLK2D(:,:)   ! " 2D


! COMMON SKOLL

INTEGER  , ALLOCATABLE ::     &
     SZIEL(:,:)    !


! COMMON MET_INI

REAL  (RealKind)                  ::     &
     DELT0,     &  ! 
     TABS0,     &  !
     PTOT0,     &  !
     RHOTOT0,   &  !
     SATT0,     &  !
     RHUM0,     &  !
     T0


! COMMON INIT_R

INTEGER              ::     &
     iread,               &  !
     ip112,               &  !
     iefeu,               &  !
     iprof,               &  !
     iqv,                 &  !
     iap_p,               &  !
     iap_e,               &  !
     iextern,             &  !
     iload,               &  !
     ientr,               &  !
     icond,               &  !
     ikoll,               &  !
     ibrea,               &  !
     iice,                &  !
     iimfr,               &  !
     ikofr,               &  !
     isolutens,           &  !
     iideal,              &  !
     inacl,               &  !
     idepo,               &  !
     imelt,               &  !
     ifrier,              &  !
     ikeis,               &  !
     iinsol                  !

! COMMON EIS_KOLLISION

INTEGER   , ALLOCATABLE  ::     &
     CONTACT(:,:)    !

REAL     (RealKind),     ALLOCATABLE  ::     &
     KOLKI(:,:), &   !
     KOLKI2D(:,:)    !

!-------------------------------------------------------------------------------
! MP MODEL CONFIG PARAMETERS
!-------------------------------------------------------------------------------
! So, Parameter bringen wir erst mal alle hier unter, dann haben wir sie beisammen.
! was mal geaendert wird, wandert in die namelists
! Parameters

INTEGER                   ::     &
     anz,        & !
     itmax,      & ! number of externally mixed aerosols 
     SMAX,       & ! number of           aerosol classes
     SIMAX,      & ! number of insoluble aerosol classes
     ipmax,      & ! horizontal resolution in the cylinder model
     jmax,       & !
     nmax,       &
     tmax,       &
! DEBUG
     icecounter, &
     filecounter,filecounter2,printcounter

REAL (RealKind)                        ::     &
     epsi

CHARACTER :: dateik*100

!!$PARAMETER (         &
!!$     itmax     = 1, &
!!$     SMAX      = 1, &
!!$     SIMAX     = 66,&
!!$     anz       = 30,&
!!$     ipmax     = 1, &
!!$     nmax      = 2, &
!!$     tmax      = 2  &
!!$           )
!!$PARAMETER (         &
!!$!     iap       = 1, &  ! chose extern via runcontrol
!!$     epsi      = 0.5_ireals,  &
!!$     icond     = 1, &
!!$     ikoll     = 2, & 
!!$     ibrea     = 0, &
!!$!
!!$     iice      = 1, &  ! Eis an/aus
!!$     idepo     = 0, &  ! Depositionsgefrieren
!!$     iimfr     = 1, &  ! Immersionsgefrieren
!!$     ikofr     = 1, &  ! Eis - Tropfen Kollision
!!$     iinsol    = 1, &  ! Kollisionsgefrieren 
!!$     ikeis     = 0, &  ! Eis - Eis Kollision
!!$     imelt     = 1, &  ! Schmelzen
!!$     ifrier    = 1, &  ! 
!!$     isolutens = 0, &  ! Oberflaechenspannung im Kelvinterm: loesliche Substanzen beruecksichtigen?
!!$     iideal    = 1  &  ! Ideale Loesung (Tropfen)
!!$           )
!!$
!!$!     iice      = 1, &  ! Eis an/aus
!!$     iice      = 0, &  ! Eis an/aus
!!$     idepo     = 0, &  ! Depositionsgefrieren
!!$     iimfr     = 1, &  ! Immersionsgefrieren
!!$!     iimfr     = 0, &  ! Immersionsgefrieren
!!$!     ikofr     = 1, &  ! Eis - Tropfen Kollision
!!$     ikofr     = 0, &  ! Eis - Tropfen Kollision
!!$     iinsol    = 0, &  ! Kontaktgefrieren 
!!$     ikeis     = 0, &  ! Eis - Eis Kollision
!!$     imelt     = 1, &  ! Schmelzen
!!$!     imelt     = 0, &  ! Schmelzen
!!$     ifrier    = 1, &  ! Gefrieren des Wassermantels
!!$!     ifrier    = 0, &  ! Gefrieren des Wassermantels
!!$     isolutens = 0, &  ! Oberflaechenspannung im Kelvinterm: 
!!$                       ! loesliche Substanzen beruecksichtigen?
!!$     iideal    = 1  &  ! Ideale Loesung (Tropfen)
!!$           )
!!$
!!$
!!$PARAMETER (dateik="DATEN/kontaktwinkel.dat")

REAL  (RealKind),       TARGET, ALLOCATABLE        ::           &
     epsilon_mp (:,:,:,:), &   ! Fraction of ...
     r_dry (:),      &   ! dry aerosol radius
     rho_old(:,:,:), &   ! for linear approximation of forcing, old value is 
                         ! stored here
     pp_old (:,:,:), &   ! as above
     p_tot_old (:,:,:), &! as above
     vmp (:,:),      &   ! droplet velocity per bin and height
     ttend (:,:,:)       ! dynamical temperature tendency

TYPE meteors
   REAL (RealKind),POINTER        :: qc (:,:,:,:,:) ! drops / particles
   REAL (RealKind),POINTER        :: bd (:,:,:,:,:) ! boundary value
   REAL (RealKind),POINTER        :: v  (:,:)       ! falling velocity
   REAL (RealKind),POINTER        :: qr (:,:)     ! accumulated precipitation at ground level
   INTEGER               :: nbins          ! number of bins
   INTEGER               :: nvars          ! number of variables per bin
   CHARACTER (LEN=10)                    :: varname        ! descriptor of meteor
ENDTYPE meteors

TYPE (meteors), POINTER   :: tmeteors (:)

!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
END MODULE data_mp_ice
!-------------------------------------------------------------------------------












