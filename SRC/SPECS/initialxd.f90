!*****************************************************************
!*                                                               *
!*    subroutine  initialxd                                      *
!*                                                               *
!*                                                               *
!*****************************************************************
SUBROUTINE initialxd(nw,qw,qws,qwa,nf,qf,qfs,qfa,qfw,ni,qi,ifreeze,&
                     ini,rho,TT,kz,anz,dz,hz_lm,iap,rrnw,rrqw,     &
                     rrnf,rrqf,epsi,pp,qv)

  INCLUDE 'HEADER3'
  
  USE utilities_mp_ice, ONLY : &
       saturation
  
  IMPLICIT NONE
      
  INTEGER :: kz,anz,iap
  
  INTEGER :: ifreeze(itmax)
  
  DOUBLE PRECISION :: nw(kz,jmax,smax,ipmax),qw(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: qws(kz,jmax,smax,ipmax),qwa(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: nf(kz,jmax,smax,ipmax),qf(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: qfs(kz,jmax,smax,ipmax),qfa(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: qfw(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: ni(kz,simax,itmax,ipmax),qi(kz,simax,itmax,ipmax)
  DOUBLE PRECISION :: rho(kz),TT(kz),ini(kz,anz),dz(kz),epsi
  DOUBLE PRECISION :: rrnw(jmax),rrqw(jmax),rrnf(jmax),rrqf(jmax)
  ! VG
  DOUBLE PRECISION :: hz_lm (kz+1),pp(kz),qv(kz)
  
  INTEGER :: k,i,j,hh,ip,si,it
  DOUBLE PRECISION :: nnw(jmax,smax,ipmax),qqw(jmax,smax,ipmax)
  DOUBLE PRECISION :: qqs(jmax,smax,ipmax),qqa(jmax,smax,ipmax)
  DOUBLE PRECISION :: nni(simax,itmax,ipmax),qqi(simax,itmax,ipmax)
  DOUBLE PRECISION :: satt,h,zentrum,hz
  
  DO k=1,kz
     ! Nullen
     ! Wassertropfen
     DO i=1,smax
        DO j=1,jmax
           DO ip = 1,ipmax
              nnw(j,i,ip) = 0.0D0
              qqw(j,i,ip) = 0.0D0
              qqs(j,i,ip) = 0.0D0
              qqa(j,i,ip) = 0.0D0
           END DO
        END DO
     END DO
     ! unloesliche AP
     DO si=1,simax
        DO it=1,itmax
           DO ip=1,ipmax
              nni(si,it,ip) = 0.0D0
              qqi(si,it,ip) = 0.0D0
           END DO
        END DO
     END DO
     ! Spektren initialisieren. iap waehlt die Parameter der Verteilungen.
     
     ! VG
     !  satt = ini(k,15)/ini(k,9) - 1.0D0
     CALL saturation (satt,tt(k),pp(k),qv(k))
     
     ! VG modify? +1?
     !  hz = DBLE(k - 1)*dz
     hz = ABS ( (hz_lm (k+1)+hz_lm(k)) / 2)
     
     ! VG 
     !  WRITE(*,*) 'Hoehe (m)             ', hz
     !  WRITE(*,*) 'Temperatur (°C)       ', TT(k) - 273.15D0
     !  WRITE(*,*) 'Uebersaettigung (%)   ', satt
     

     CALL ap_newxd(nnw,qqw,qqs,qqa,nni,qqi,ifreeze,satt,TT(k),satt,TT(k),    &
                   rho(k),rho(40),hz,iap,epsi)

     DO i=1,smax
        DO j=1,jmax
           DO ip=1,ipmax
              nw(k,j,i,ip)  = nnw(j,i,ip)
              qw(k,j,i,ip)  = qqw(j,i,ip)
              qws(k,j,i,ip) = qqs(j,i,ip)
              qwa(k,j,i,ip) = qqa(j,i,ip)
              nf(k,j,i,ip)  = 0.0D0
              qf(k,j,i,ip)  = 0.0D0
              qfs(k,j,i,ip) = 0.0D0
              qfa(k,j,i,ip) = 0.0D0
              qfw(k,j,i,ip) = 0.0D0
           END DO
        END DO
     END DO
     DO si=1,simax
        DO it=1,itmax
           DO ip=1,ipmax
              ni(k,si,it,ip) = nni(si,it,ip)
              qi(k,si,it,ip) = qqi(si,it,ip)
           END DO
        END DO
     END DO
     ! k-loop
  END DO
  
  !OPEN(10,FILE="RESULTATE/avertu")
  !OPEN(20,FILE="RESULTATE/averto")
  DO j=1,jmax
     rrnw(j) = 0.0D0
     rrqw(j) = 0.0D0
     rrnf(j) = 0.0D0
     rrqf(j) = 0.0D0
     !  DO i=1,smax
     !    WRITE(10,100) rmitte(j),rsmitte(i),MAX(nw(1,j,i,1),1.0D-30),         &
     !                  MAX(qw(1,j,i,1),1.0D-30),MAX(qws(1,j,i,1),1.0D-30),    &
     !                  MAX(qwa(1,j,i,1),1.0D-30)
     !    WRITE(20,100) rmitte(j),rsmitte(i),MAX(nw(kz,j,i,1),1.0D-30),        &
     !                  MAX(qw(kz,j,i,1),1.0D-30),MAX(qws(kz,j,i,1),1.0D-30),  &
     !                  MAX(qwa(kz,j,i,1),1.0D-30)
     !  END DO
     !  IF(smax.NE.1) THEN
     !    WRITE(10,*)
     !    WRITE(20,*)
     !  END IF
  END DO
  !CLOSE(10)
  !CLOSE(20)
  
  !OPEN(30,FILE="RESULTATE/aivertu")
  !OPEN(40,FILE="RESULTATE/aiverto")
  !DO si=1,simax
  !  WRITE(30,200) rsmitte(si),MAX(ni(1,si,1,1),1.0D-30),                    &
  !                MAX(qi(1,si,1,1),1.0D-30)
  !  WRITE(40,200) rsmitte(si),MAX(ni(kz,si,1,1),1.0D-30),                   &
  !                MAX(qi(kz,si,1,1),1.0D-30)
  !END DO
  !CLOSE(30)
  !CLOSE(40)
  
100 FORMAT(6(E13.6,2X))
200 FORMAT(3(E13.6,2X))
  
  RETURN
END SUBROUTINE initialxd




