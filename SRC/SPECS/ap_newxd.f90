SUBROUTINE ap_newxd(nw,qw,qs,qa,nwins,qains,ifreeze,satt,tabs,satt_e,   &
                    tabs_e,rho,rho0,hhz,iap,eps)

  ! control module for the soluble and insoluble particle distributions 
  ! and initialisation of QW, QS, QA, NW, QAINS, NWINS
  ! called by wmain.f
  ! calling apvert.f

  USE data_runcontrol, ONLY : dnap_init
  
  INCLUDE 'HEADER2'
  
  IMPLICIT NONE 
  
  INTEGER :: I,J,imode,IP,IT,ifreeze(ITMAX)
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: epsi
  DOUBLE PRECISION :: NWINS(SIMAX,ITMAX,IPMAX),QAINS(SIMAX,ITMAX,IPMAX)
  DOUBLE PRECISION :: TABS,SATT,TABS_E,SATT_E
  DOUBLE PRECISION :: NA0,RA0,SIGMA
  
  ! maximum number of modes for the AP distributions
  INTEGER :: mmax
  PARAMETER(mmax = 4)
  DOUBLE PRECISION :: dp_ap(mmax),dn_ap(mmax),sig_ap(mmax)
  
  INTEGER :: iap,si
  DOUBLE PRECISION :: rho,rho0,hhz,hhmax,eps
    
  ! parcel 
  IP = 1
  ! soluble particles
  ! Jaenicke (1988) siehe PK 97 S.272
  IF(iap.EQ.1) THEN
     it    = 0
! VG trimodale verteilung
     imode = 1
!      imode = 3
    hhmax = 2000.0D0
     !  hhmax = 0.0D0
!!$  WRITE(*,*) "maritim aerosol            ",iap
!!$  WRITE(*,*) "mode(s)                    ",imode
     ! dn_ap in [1/kg] ^= 1D-3 [1/l]
     dn_ap(1) = 133.0D6
     dn_ap(2) = 66.6D6
     dn_ap(3) = 3.06D6
     ! dp_ap in [m]
     dp_ap(1) = 3.9D-9
     dp_ap(2) = 133.0D-9
     dp_ap(3) = 290.0D-9
     ! Standardabweichung
     sig_ap(1) = 4.54D0
     sig_ap(2) = 1.62D0
     sig_ap(3) = 2.49D0
     epsi = eps
     
     ! VG Kreidenweis
     IF (imode == 1) THEN 

        dn_ap(1)  = dnap_init*1.D6 ! Default: 566.0D6 
        dp_ap(1)  = 30.0D-9
        sig_ap(1) = 2.0D0

!print *,dp_ap(1)
     ENDIF

     CALL apvertxd(nw,qw,qs,qa,nwins,qains,satt,tabs,ip,it,imode,         &
                   iap,dn_ap,dp_ap,sig_ap,epsi,hhmax,hhz,rho,rho0)
  END IF
  ! insoluble particles
  IT = 1
  ! imode=2
  imode = 1
  hhmax = 2000.0D0
  !hhmax = 0.0D0
!!$WRITE(*,*) "insoluble                    "
!!$WRITE(*,*) "Young BB aerosol           "
!!$WRITE(*,*) "mode(s)                    ",imode
  !! dn_ap in [1/kg]
  ! dn_ap(1) = 45.0D6
  ! dn_ap(2) = 5.0D6 
  ! dn_ap(1) = 0.0D0
  ! dn_ap(2) = 0.0D0
! "default" is 1e6
  dn_ap(1) = 1.0D6
!print *,'dn_ap for insolubles',dn_ap(1)
  
  !! dp_ap in [m]         
  ! dp_ap(1) = 15.0D-9
  ! dp_ap(2) = 60.0D-9
  ! dp_ap(1) = 0.0D0
  ! dp_ap(2) = 0.0D0
  dp_ap(1) = 0.5D-6
  !! Standardabweichung
  ! sig_ap(1) = 1.5D0
  ! sig_ap(2) = 1.8D0
  ! sig_ap(1) = 0.0D0
  ! sig_ap(2) = 0.0D0
  sig_ap(1) = 1.5D0
  epsi = 0.0D0
  CALL apvertxd(nw,qw,qs,qa,nwins,qains,satt,tabs,ip,it,imode,           &
                iap,dn_ap,dp_ap,sig_ap,epsi,hhmax,hhz,rho,rho0)

  ifreeze(IT) = iinsol

  
!WRITE(*,*) "IN type                    ",IT,ifreeze(IT),iinsol
  
  ! environment
  DO i=1,smax
     DO j=1,jmax
        nw(j,i,ipmax) = nw(j,i,1)
        qw(j,i,ipmax) = qw(j,i,1)
        qs(j,i,ipmax) = qs(j,i,1)
        qa(j,i,ipmax) = qa(j,i,1)
     END DO
  END DO
  DO si=1,simax
     DO it=1,itmax
        nwins(si,it,ipmax) = nwins(si,it,1)
        qains(si,it,ipmax) = qains(si,it,1)
     END DO
  END DO
  
  RETURN
END SUBROUTINE ap_newxd








