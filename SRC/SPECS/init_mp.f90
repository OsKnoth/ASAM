! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ External procedure for organizing the calls to the physics packages
!------------------------------------------------------------------------------

SUBROUTINE init_mp (                  &
                    nw_lm,      & ! LM   
                    qw_lm,      & ! LM   
                    qws_lm,     & ! LM   
                    qwa_lm,     & ! LM   
                    vtw,        & ! LM
                    nf_lm,      & ! LM   
                    qf_lm,      & ! LM   
                    qfs_lm,     & ! LM   
                    qfa_lm,     & ! LM   
                    qfw_lm,     & ! LM   
                    vtf,        & ! LM
                    rho,        & ! LM   
                    TT,         & ! LM   
                    pp,         & ! LM   
                    ni_lm,      & ! LM   
                    qia_lm,     & ! LM   
                    kz,         & ! LM
                    dz_lm,      & ! LM     
                    jmax,       & ! LM/MP
                    qv          & ! LM
                                       )


!------------------------------------------------------------------------------
!
! Description:
! This procedure initializes microphysical variables
!
! Method:
!
! Current Code Owner: IfT Verena Gruetzun
!  phone:  +49  341 235 2826
!  fax:    +49  341 235 2139
!  email:  gruetzun@tropos.de
!
! History:
! Version    Date       Name
! ---------- ---------- ----

! !VERSION!  !DATE!     <Your name>
!  <Modification comments>
!
! Code Description:
! Language: Fortran 90.
! Software Standards: "European Standards for Writing and
! Documenting Exchangeable Fortran 90 Code".
!==============================================================================

USE Kind_Mod

USE data_mp_ice,        ONLY:     T0

USE data_mp_ice,        ONLY:     &
     ini,     & ! Anfangswerte, werden ueberall wieder mit reingegeben
     vt ,     & ! wahrscheinlich nutzls, frueher vbest, jetzt nicht mehr verwendet (lt. Harald)
     miv,     & ! Kontaktwinkel, wird aus datei eingelesen
     rrnw,    & ! Niederschlag sammeln
     rrqw,    & ! "
     rrnf,    & ! "
     rrqf,    & ! "
     ifreeze    ! freezing par, Bakterien, Russ etc. pp

USE data_runcontrol,    ONLY: &
     iap

USE data_mp_ice,        ONLY: &
     itmax,             &
     smax,              &
     simax,             &
     ipmax,             &
     dateik,            &
     epsi               





!==============================================================================
IMPLICIT NONE
!==============================================================================


! define parameters
INTEGER   ::   &
     anz

PARAMETER (anz=30)


INTEGER    ::   &
     J12,          & ! wird wahrscheinlich nichts mit gemacht
     jmax,         & ! LM
     kz,           & ! LM
     k               ! loop index

REAL    (RealKind)           ::   &
! macroscopic quantities
     TT (kz),        & ! LM
     pp (kz),        & ! LM
     rho(kz),        & ! LM
     dz_lm (kz+1),   & ! LM
     dz (kz),        & ! LM
     qv (kz),        & ! LM

! microphysical quantities
     R0W,     &      ! wird laut Harald nichts mit gemacht
!! water
     nw_lm  (kz,jmax),      & ! LM   
     qw_lm  (kz,jmax),      & ! LM   
     qws_lm (kz,jmax),      & ! LM   
     qwa_lm (kz,jmax),      & ! LM   
     nw(kz,jmax,smax,ipmax),      & ! MP
     qw(kz,jmax,smax,ipmax),      & ! MP
     qws(kz,jmax,smax,ipmax),     & ! MP
     qwa(kz,jmax,smax,ipmax),     & ! MP
     vtw(kz,jmax),              & ! MP
!! ice
     nf(kz,jmax,smax,ipmax),      & ! MP
     qf(kz,jmax,smax,ipmax),      & ! MP
     qfs(kz,jmax,smax,ipmax),     & ! MP
     qfa(kz,jmax,smax,ipmax),     & ! MP
     qfw(kz,jmax,smax,ipmax),     & ! MP
     vtf(kz,jmax),              & ! MP
     nf_lm  (kz,jmax),      & ! LM   
     qf_lm  (kz,jmax),      & ! LM   
     qfs_lm (kz,jmax),      & ! LM   
     qfa_lm (kz,jmax),      & ! LM   
     qfw_lm (kz,jmax),      & ! LM   
!!unsoluble particles
     ni (kz,simax,itmax,ipmax),   & ! MP
     qia(kz,simax,itmax,ipmax),   & ! MP
     ni_lm  (kz,jmax),      & ! LM   
     qia_lm (kz,jmax)         ! LM   

CHARACTER :: blah


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Which data is actually needed globally? E.g. constants, velocities and such
! generally all from common blocks
! mp values are already known via LM
! this routine will be called with the needed values out of LM

!------------------------------------------------------------------------------
!- End of header
!------------------------------------------------------------------------------

!------------------------------------------------------------------------------
!- Begin Subroutine organize_microphysics
!------------------------------------------------------------------------------

! This routine is called once to initialize all necessary 
! variables for the microphysics
! the calls are taken from akmxd.f90

! modify
!!! Hier noch mit den Dimensionen arbeiten


! LM 2 Microphysics routine

! initialize jmax for subroutines

IF (smax  .NE. 1) THEN
   print *,'WARNING! DIMENSION WRONG, SMAX .NE. 1' 
   STOP
ENDIF
IF (ipmax .NE. 1) THEN
   print *,'WARNING! DIMENSION WRONG, IPMAX .NE. 1'
   STOP
ENDIF

! water 
nw (:,:,smax,ipmax) = nw_lm (:,:)
qw (:,:,smax,ipmax) = qw_lm (:,:)
qws(:,:,smax,ipmax) = qws_lm(:,:)
qwa(:,:,smax,ipmax) = qwa_lm(:,:)

! ice
nf (:,:,smax,ipmax) = nf_lm (:,:)
qf (:,:,smax,ipmax) = qf_lm (:,:)   
qfs(:,:,smax,ipmax) = qfs_lm(:,:)   
qfa(:,:,smax,ipmax) = qfa_lm(:,:)   
qfw(:,:,smax,ipmax) = qfw_lm(:,:)

! insoluble aerosol
ni (:,:,smax,ipmax) = ni_lm (:,:)   
qia(:,:,smax,ipmax) = qia_lm(:,:)  

! layer thickness
DO k=1,kz
   dz (k) = ABS(dz_lm (k+1) - dz_lm(k))
ENDDO





CALL konstxd

J12 = 1
R0W = 1.d0
CALL sconstxd(J12,R0W)
CALL kernelxd
CALL contact_freezing_matrixxd

!!$
!!$print *,'nw',nw(1,1,1,1)
!!$print *,'qw',qw(1,1,1,1)
!!$print *,'qws',qws(1,1,1,1)
!!$print *,'qwa',qwa(1,1,1,1)
!!$print *,'nf',nf(1,1,1,1)
!!$print *,'qf',qf(1,1,1,1)
!!$print *,'qfs',qfs(1,1,1,1)
!!$print *,'qfa',qfa(1,1,1,1)
!!$print *,'qfw',qfw(1,1,1,1)
!!$print *,'ni',ni(1,1,1,1)
!!$print *,'qia',qia(1,1,1,1)
!!$print *,'ifreeze',ifreeze(1)
!!$print *,'ini',ini(1,1)
!!$print *,'rho',rho(1)
!!$print *,'TT',TT(1)
!!$print *,'kz',kz
!!$print *,'anz',anz
!!$print *,'dz',dz(1)
!!$print *,'iap',iap
!!$print *,'rrnw',rrnw(1)
!!$print *,'rrqw',rrqw(1)
!!$print *,'rrnf',rrnf(1)
!!$print *,'rrqw',rrqf(1)
!!$print *,'epsi',epsi
!!$print *,'pp',pp(1)


CALL initialxd(nw,      & ! LM
               qw,      & ! LM
               qws,     & ! LM
               qwa,     & ! LM
               nf,      & ! LM
               qf,      & ! LM
               qfs,     & ! LM
               qfa,     & ! LM
               qfw,     & ! LM
               ni,      & ! LM
               qia,     & ! LM
               ifreeze, & ! mp_data
               ini,     & ! mp_data
               rho,     & ! LM
               TT,      & ! LM
               kz,      & ! LM
               anz,     & ! Parameter from MP
               dz,      & ! LM
               dz_lm,   & ! height in LM
               iap,     & ! data_runcontrol
               rrnw,    & ! mp_data
               rrqw,    & ! mp_data
               rrnf,    & ! mp_data
               rrqf,    & ! mp_data
               epsi,    & ! mp_data parameter
               pp,      & ! LM
               qv       & ! LM
                              )

CALL vbestxd  (vt,      & ! useless
               kz,      & ! LM
               dz)        ! LM

CALL vstokesxd(vtw,     & ! LM
               vtf,     & ! LM
               kz,      & ! LM
               dz,      & ! LM
               pp,      & ! LM
               rho,     & ! LM
               TT,      & ! LM
               T0)        ! Parameter from MP

CALL kontaktxd(dateik,  & ! Parameter from MP
               miv,     & ! Kontaktwinkel from MP
               itmax)     ! Parameter from MP



nw_lm (:,:) = nw (:,:,smax,ipmax) 
qw_lm (:,:) = qw (:,:,smax,ipmax) 
qws_lm(:,:) = qws(:,:,smax,ipmax) 
qwa_lm(:,:) = qwa(:,:,smax,ipmax) 

! ice
nf_lm (:,:) = nf (:,:,smax,ipmax) 
qf_lm (:,:) = qf (:,:,smax,ipmax) 
qfs_lm(:,:) = qfs(:,:,smax,ipmax) 
qfa_lm(:,:) = qfa(:,:,smax,ipmax) 
qfw_lm(:,:) = qfw(:,:,smax,ipmax) 

! insoluble aerosol
ni_lm (:,:) = ni (:,:,smax,ipmax) 
qia_lm(:,:) = qia(:,:,smax,ipmax) 



!------------------------------------------------------------------------------
! End of module procedure src_microphysics
!------------------------------------------------------------------------------
END SUBROUTINE init_mp



















