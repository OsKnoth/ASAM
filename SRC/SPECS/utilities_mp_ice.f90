! File %M% from Library %Q%
! Version %I% from %G% extracted: %H%
!+ Variables for the computation of synthetic satellite images
!------------------------------------------------------------------------------

MODULE utilities_mp_ice

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
USE Transport_Mod

USE data_constants ,  ONLY :   &
     b11=>b1,           &
     b2w,               &
     b3,                &
     b4w,               &
     rdv,               &  ! r_d/r_v
     o_m_rdv               ! 1 - r_d/r_v

USE data_mp_ice,      ONLY :   &
     ini    ,  & ! Anfangswerte, werden ueberall wieder mit reingegeben
     vt     ,  & ! wahrscheinlich nutzls, frueher vbest, jetzt nicht mehr verwendet (lt. Harald)
     miv    ,  & ! Kontaktwinkel, wird aus datei eingelesen
     rrnw   ,   & ! Niederschlag sammeln
     rrqw   ,   & ! "
     rrnf   ,   & ! "
     rrqf   ,   & ! "
     ifreeze,   & ! freezing par, Bakterien, Russ etc. pp
     DIFF21,  & ! difference of two boundaries
     MGRENZ,  & ! boundaries of classes - mass 
     MMITTE,  & ! center of classes     - mass 
     RGRENZ,  & ! boundaries of classes - radius
     RMITTE,  & ! center of classes     - radius
     SDIFF21, & !
     SGRENZ,  & !
     SMITTE,  & !
     RSGRENZ, & !
     RSMITTE, & !
     mquer,   &
     rquer,   &
     squer,   &
     mfquer,  &
     rfquer,  &
     sfquer,  &
     anz,        & !
     jmax,       &
     itmax,      & ! number of externally mixed aerosols 
     SMAX,       & ! number of           aerosol classes
     SIMAX,      & ! number of insoluble aerosol classes
     ipmax,      & ! horizontal resolution in the cylinder model
     KOLKI     , & !
     KOLKI2D,    &!
     SZIEL,      & !
     KOLK,       &  ! collision kernel
     KOLK2D,     & ! " 2D
     KZIEL,      & ! target bin for collision
     CONTACT    !

CONTAINS 

!-------------------------------------------------------------------------------
SUBROUTINE saturation(satt,t,p,qv)
!-------------------------------------------------------------------------------
  IMPLICIT NONE
  REAL (Real8Kind) :: &
       t ,        &
       p ,        &
       qv

  REAL (Real8Kind) :: &
       zpvs,  zqvs

  REAL (Real8Kind),INTENT(OUT) :: &
       satt 

      zpvs = b11*EXP( b2w*(t-b3) / (t-b4w) )
      zqvs = rdv*zpvs / (p - o_m_rdv*zpvs)
      satt = qv/zqvs - 1.0d0

END SUBROUTINE saturation


!-------------------------------------------------------------------------------
SUBROUTINE mp_allocate (kz,n_class)
!-------------------------------------------------------------------------------
  IMPLICIT NONE  
  INTEGER                  ::     &
     kz,n_class

! End of header
!-------------------------------------------------------------------------------

  ALLOCATE ( ini(kz,anz)       ); ini = 0.0d0
  ALLOCATE ( vt (kz,n_class)   ); vt  = 0.0d0
  ALLOCATE ( miv(itmax)        ); miv = 0.0d0
  ALLOCATE ( rrnw(n_class)     ); rrnw = 0.0d0
  ALLOCATE ( rrqw(n_class)     ); rrqw = 0.0d0
  ALLOCATE ( rrnf(n_class)     ); rrnf = 0.0d0
  ALLOCATE ( rrqf(n_class)     ); rrqf = 0.0d0
  ALLOCATE ( mquer(1,jmax,smax) ); mquer  = 0.0d0 
  ALLOCATE ( rquer(1,jmax,smax) ); rquer  = 0.0d0 
  ALLOCATE ( squer(1,jmax,smax) ); squer  = 0.0d0 
  ALLOCATE ( mfquer(1,jmax,smax)); mfquer = 0.0d0
  ALLOCATE ( rfquer(1,jmax,smax)); rfquer = 0.0d0
  ALLOCATE ( sfquer(1,jmax,smax)); sfquer = 0.0d0

  ALLOCATE ( ifreeze(itmax)    ); ifreeze = 0
  !
  ALLOCATE ( DIFF21(n_class)   ); diff21 = 0.0d0
  ALLOCATE ( MGRENZ(n_class+1) ); mgrenz = 0.0d0
  ALLOCATE ( MMITTE(n_class)   ); mmitte = 0.0d0
  ALLOCATE ( RGRENZ(n_class+1) ); rgrenz = 0.0d0
  ALLOCATE ( RMITTE(n_class)   ); rmitte = 0.0d0
  ALLOCATE ( SDIFF21(SIMAX)    ); sdiff21 = 0.0d0
  ALLOCATE ( SGRENZ(SIMAX+1)   ); sgrenz = 0.0d0
  ALLOCATE ( SMITTE(SIMAX)     ); smitte = 0.0d0
  ALLOCATE ( RSGRENZ(SIMAX+1)  ); rsgrenz = 0.0d0
  ALLOCATE ( RSMITTE(SIMAX)    ); rsmitte = 0.0d0
  !
  ALLOCATE ( KZIEL(n_class,n_class)      ); kziel  = 0
  ALLOCATE ( KOLK(n_class,n_class)       ); kolk   = 0.0d0
  ALLOCATE ( KOLK2D(n_class+1,n_class+1) ); kolk2d = 0.0d0
  ALLOCATE ( SZIEL(SIMAX,SIMAX)          ); sziel  = 0
  ALLOCATE ( CONTACT(n_class,SMAX)       ); contact= 0
  ALLOCATE ( KOLKI(n_class,n_class)      ); kolki  = 0.0d0
  ALLOCATE ( KOLKI2D(n_class+1,n_class+1)); kolki2d = 0.0d0
  
!-------------------------------------------------------------------------------
END SUBROUTINE mp_allocate
!-------------------------------------------------------------------------------


!-------------------------------------------------------------------------------
SUBROUTINE mp_deallocate
!-------------------------------------------------------------------------------
  DEALLOCATE ( ini  )
  DEALLOCATE ( vt   )
  DEALLOCATE ( miv  )
  DEALLOCATE ( rrnw )
  DEALLOCATE ( rrqw )
  DEALLOCATE ( rrnf )
  DEALLOCATE ( rrqf )
  DEALLOCATE ( mquer )  
  DEALLOCATE ( rquer )
  DEALLOCATE ( squer )
  DEALLOCATE ( mfquer )
  DEALLOCATE ( rfquer )
  DEALLOCATE ( sfquer )

  DEALLOCATE ( ifreeze )
  !
  DEALLOCATE ( DIFF21 )
  DEALLOCATE ( MGRENZ )
  DEALLOCATE ( MMITTE ) 
  DEALLOCATE ( RGRENZ )
  DEALLOCATE ( RMITTE )
  DEALLOCATE ( SDIFF21 )
  DEALLOCATE ( SGRENZ )
  DEALLOCATE ( SMITTE ) 
  DEALLOCATE ( RSGRENZ )
  DEALLOCATE ( RSMITTE )
  !
  DEALLOCATE ( KZIEL )
  DEALLOCATE ( KOLK )
  DEALLOCATE ( KOLK2D )
  DEALLOCATE ( SZIEL )
  DEALLOCATE ( CONTACT )
  DEALLOCATE ( KOLKI )
  DEALLOCATE ( KOLKI2D )

END SUBROUTINE mp_deallocate




END MODULE utilities_mp_ice
