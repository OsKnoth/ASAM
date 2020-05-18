!*****************************************************************
!*                                                               *
!*    subroutine upstreamxd                                      *
!*                                                               *
!*    Berechnet aus x(t) den Wert zum Zeitpunkt t+1 gemaess      *
!*    der Advektionsgleichung:                                   *
!*                                                               *
!*         d(rho*x)/dt = -d(w*rho*x)/dz                          *
!*                                                               *
!*    wm ist die Geschwindigkeit von w an den Zellkanten, sie    *
!*    muss im aufrufenden Programm bestueckt werden.             *
!*                                                               *
!*                                                               *
!*****************************************************************
SUBROUTINE upstreamxd(x,dx,rho,wm,dt,dz,kz,ipmax)

IMPLICIT NONE

INTEGER :: kz,ipmax
DOUBLE PRECISION :: dt,dz
DOUBLE PRECISION :: x(kz,ipmax),dx(kz,ipmax),rho(kz),wm(kz+1,ipmax)

INTEGER :: k,ip
DOUBLE PRECISION :: xx(kz,ipmax),fluss(kz+1,ipmax)
      

! Feld wird mit der Dichte multipliziert.
DO ip=1,ipmax
  DO k=1,kz
    xx(k,ip) = x(k,ip)*rho(k)
  END DO
END DO

! Bestueckung des Flussvektors
DO ip=1,ipmax
  DO k=2,kz
    fluss(k,ip) = (wm(k,ip)*(xx(k-1,ip) + xx(k,ip))                    &
                + DABS(wm(k,ip))*(xx(k-1,ip) - xx(k,ip)))/2.0D0
  END DO
  fluss(1,ip) = 0.0D0
  fluss(kz+1,ip) = 0.0D0
END DO
      
! Berechnug der Vorwaertsdifferenz und durch Dichte teilen
DO ip=1,ipmax
  DO k=1,kz  
    dx(k,ip) = -(fluss(k+1,ip) - fluss(k,ip))/dz
    IF(DABS(dx(k,ip)).LT.1.0D-8) dx(k,ip) = 0.0D0
    xx(k,ip) = (xx(k,ip) + dx(k,ip)*dt)/rho(k)
    dx(k,ip) = dx(k,ip)/rho(k)
  END DO
END DO

           
RETURN
END
