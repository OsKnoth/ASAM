!*****************************************************************
!*                                                               *
!*    subroutine  teta2txd                                       *
!*                                                               *
!*****************************************************************
SUBROUTINE teta2txd(teta,kz,ipmax,t,p,kappa,frage)

IMPLICIT NONE

INTEGER :: kz,ipmax,frage
DOUBLE PRECISION :: teta(kz,ipmax),t(kz,ipmax),p(kz),kappa

INTEGER :: k,ip


DO ip=1,ipmax        
  DO k=1,kz
    IF(frage.EQ.1) t(k,ip) = teta(k,ip)*(p(k)/p(1))**kappa
    IF(frage.EQ.-1) teta(k,ip) = t(k,ip)*(p(1)/p(k))**kappa
  END DO
END DO

           
RETURN
END
