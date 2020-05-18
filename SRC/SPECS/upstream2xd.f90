!*****************************************************************
!*                                                               *
!*    subroutine upstream2xd                                     *
!*                                                               *
!*    Berechnet aus x(t) den Wert zum Zeitpunkt t+1 gemaess      *
!*    der Advektionsgleichung:                                   *
!*                                                               *
!*         d(x)/dt = -d((w+vt)*x)/dz                             *
!*                                                               *
!*    wm ist die Geschwindigkeit von w an den Zellkanten,        *
!*    vt ist die Endfallgeschwindigkeit, beide muessen im        *
!*    aufrufenden Programm bestueckt werden.                     *
!*                                                               *
!*                                                               *
!*****************************************************************
SUBROUTINE upstream2xd(x,dx,wm,vt,dt,dz,kz,jmax,smax,ipmax,xyz)

IMPLICIT NONE

INTEGER :: kz,jmax,smax,ipmax
DOUBLE PRECISION :: dt,dz,xyz
DOUBLE PRECISION :: x(kz,jmax,smax,ipmax),dx(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: wm(kz+1,ipmax),vt(kz,jmax)

INTEGER :: k,j,i,ip
DOUBLE PRECISION :: wmt,fluss(kz+1,ipmax)
      

DO ip=1,ipmax
  DO j=1,jmax
    DO i=1,smax
! Bestueckung des Flussvektors
      DO k=2,kz
        wmt = wm(k,ip) + xyz*vt(k,j)
        fluss(k,ip) = (wmt*(x(k-1,j,i,ip) + x(k,j,i,ip))              &
                    + DABS(wmt)*(x(k-1,j,i,ip) - x(k,j,i,ip)))/2.0D0
      END DO
      fluss(1,ip) = 0.0D0
      fluss(kz+1,ip) = 0.0D0      
! Berechnug der Vorwaertsdifferenz
      DO k=1,kz
            dx(k,j,i,ip) = -(fluss(k+1,ip) - fluss(k,ip))/dz
      END DO
    END DO
  END DO
END DO

           
RETURN
END
