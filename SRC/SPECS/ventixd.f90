! Ventilationskoeffizient nach Pruppacher und Klett (1997) S.541

DOUBLE PRECISION FUNCTION ventixd(TT,pp,rho,radius,vt)

USE data_mp_ice, ONLY : T0

IMPLICIT NONE
DOUBLE PRECISION :: TT,pp,rho,radius,vt

DOUBLE PRECISION :: p0,Nsc,Nre,Nscre
PARAMETER(p0 = 101325.0D0)

! Functions
DOUBLE PRECISION :: etaxd,dvxd


! Nsc ist die Schmidt-Zahl fuer Wasserdampf (siehe PK 97, S.538)
Nsc = etaxd(TT)/(rho*dvxd(TT,pp))

! Nre ist die Reynolds Zahl (siehe PK 97, S.364)
Nre = 2.0D0*DABS(vt)*radius*rho/etaxd(TT)

Nscre = Nsc**(1.0D0/3.0D0)*Nre**(1.0D0/2.0D0)
IF(Nscre.LT.1.4D0) ventixd = 1.0D0 + 0.108D0*Nscre**2.0D0
IF(Nscre.GE.1.4D0.AND.Nscre.LE.51.4D0) THEN
  ventixd = 0.78D0 + 0.308D0*Nscre
END IF
IF(Nscre.GT.51.4D0) THEN
!  WRITE(*,*) 'Nscre = ',Nscre
!VG  ventixd = 1.0D0
  ventixd = 0.78D0 + 0.308D0*51.4D0
END IF


RETURN
END
