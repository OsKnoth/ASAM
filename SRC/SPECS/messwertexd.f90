SUBROUTINE messwertexd(dateim,kz,ipmax,dz,Rd,g,T0,kappa,TT,RR,pp,teta,   &
                       qs,qq,ww,rho,ini,nanz,dTT,hdT,kelvin)

IMPLICIT NONE

INTEGER :: kz,nanz,ipmax
DOUBLE PRECISION :: dz,Rd,g,T0,kappa,dTT,hdT
DOUBLE PRECISION :: pp(kz),TT(kz),RR(kz),qs(kz),qq(kz,ipmax)
DOUBLE PRECISION :: ww(kz,ipmax),teta(kz,ipmax),rho(kz),ini(kz,nanz)
CHARACTER :: dateim*100,kelvin*1

INTEGER :: i,ii,k,anz
DOUBLE PRECISION, ALLOCATABLE :: druck(:),Temp(:),tau(:),zh(:),Rf(:)
DOUBLE PRECISION :: TTd(kz),hoehe(kz),a,b,dp,xlnp,dtdp,dtddp,drfdp,Tx


OPEN(10,FILE=dateim)
READ(10,*) anz
ALLOCATE(druck(anz),Temp(anz),tau(anz),zh(anz),Rf(anz))
lesen: DO i=1,anz
  READ(10,*) druck(i),Temp(i),tau(i)
  IF(kelvin.EQ.'K') Tx = T0
  IF(kelvin.EQ.'C') Tx = 0.0D0
  Temp(i) = Temp(i) - Tx
  tau(i) = tau(i) - Tx
  WRITE(*,*) druck(i),Temp(i),tau(i)
  a = 7.45D0*tau(i)/(235.0D0 + tau(i))
  b = 7.45D0*Temp(i)/(235.0D0 + Temp(i))
  Rf(i) = 10.0D0**(a - b)
  Temp(i) = Temp(i) + T0
  tau(i) = tau(i) + T0
END DO lesen
CLOSE(10)

zh(1) = 0.0D0
DO i=2,anz
  zh(i) = zh(i-1) + Rd*(Temp(i) + Temp(i-1))/(g*2.0D0)*    &
          DLOG(druck(i-1)/druck(i))
END DO

hoehe(1) = zh(1)
pp(1) = druck(1)
TT(1) = Temp(1)
TTd(1) = tau(1)
RR(1) = Rf(1)
teta(1,1) = TT(1)*(pp(1)/pp(1))**kappa
teta(1,2) = teta(1,1)
qs(1) = (3.8D0*DEXP(17.27D0*(TT(1) - T0)/(TT(1) - 35.86D0)))/(pp(1)/100.0D0)
qq(1,1) = RR(1)*qs(1)
qq(1,2) = RR(1)*qs(1)
ww(1,1) = 0.0D0
ww(1,2) = 0.0D0
rho(1) = pp(1)/(Rd*TT(1))
DO k=2,kz
  hoehe(k) = DBLE(k-1)*dz
  DO i=2,anz
    ii = i
    IF(hoehe(k).LT.zh(i)) GOTO 111
  END DO
111 CONTINUE
  xlnp = DLOG(druck(ii)/druck(ii-1))
  dtdp = (Temp(ii) - Temp(ii-1))/xlnp
  dtddp = (tau(ii) - tau(ii-1))/xlnp
  drfdp = (Rf(ii) - Rf(ii-1))/xlnp
  pp(k) = druck(ii-1)*EXP(-g*(hoehe(k) - zh(ii-1))/(Rd*(Temp(ii-1)   &
        + Temp(ii))/2.0D0))
  dp = DLOG(pp(k)/druck(ii-1))
  TT(k) = Temp(ii-1) + dtdp*dp
  TTd(k) = tau(ii-1) + dtddp*dp
  RR(k) = Rf(ii-1) + drfdp*dp
  teta(k,1) = TT(k)*(pp(1)/pp(k))**kappa
  teta(k,2) = teta(k,1)
! Fuer die Berechnung von qs wird der Druck in hPa verwendet,
! d.h.: p(z)/100. Es gilt qs = 0.622*es/p, der Saettigungsdampfdruck (es)
! wird nach der Magnusformel bestimmt.
  qs(k) = (3.8D0*DEXP(17.27D0*(TT(k) - T0)/(TT(k) - 35.86D0)))/(pp(k)/100.0D0)
  qq(k,1) = RR(k)*qs(k)
  qq(k,2) = RR(k)*qs(k)
  ww(k,1) = 0.0D0
  ww(k,2) = 0.0D0
! Dichte trockener Luft
  rho(k) = pp(k)/(Rd*TT(k))
END DO

CALL impulsxd(TT,dTT,teta,qq,pp,RR,kz,ipmax,hdT,dz,kappa,T0)

OPEN(20,FILE='DATEN/startfelder')
DO i=1,kz
  ini(i,1) = ww(i,1)
  ini(i,2) = teta(i,1)
  ini(i,3) = qq(i,1)
  ini(i,4) = 0.0D0
  ini(i,5) = ww(i,2)
  ini(i,6) = teta(i,2)
  ini(i,7) = qq(i,2)
  ini(i,8) = 0.0D0
  ini(i,9) = qs(i)
  ini(i,10) = TT(i)
  ini(i,15) = qq(i,1)
  WRITE(20,100) (i-1)*INT(dz),TT(i),TTd(i),RR(i),pp(i),rho(i),ww(i,1),    &
                ww(i,2),teta(i,1),teta(i,2),qs(i),qq(i,1),qq(i,2)
END DO
CLOSE(20)

100  FORMAT(I6,14(2X,D13.6))

RETURN
END
