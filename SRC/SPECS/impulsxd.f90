SUBROUTINE impulsxd(TT,dTT,teta,qq,pp,RR,kz,ipmax,hdT,dz,kappa,T0)
  
  IMPLICIT NONE
  
  INTEGER :: kz,ipmax
  DOUBLE PRECISION :: dTT,dz,hdT,kappa,T0
  DOUBLE PRECISION :: pp(kz),TT(kz),RR(kz),qq(kz,ipmax),teta(kz,ipmax)
  
  INTEGER :: k,hhdT
  DOUBLE PRECISION :: qs(kz)
  
  
  hhdT = hdT/dz
  DO k=1,hhdT
     TT(k) = TT(k) + dTT
     teta(k,1) = TT(k)*(pp(1)/pp(k))**kappa  
     !  qs(k) = (3.8D0*DEXP(17.27D0*(TT(k) - T0)/(TT(k) - 35.86D0)))/(pp(k)/100.0D0)
     !  qq(k,1) = RR(k)*qs(k)
  END DO
  
  
  RETURN
END SUBROUTINE impulsxd
