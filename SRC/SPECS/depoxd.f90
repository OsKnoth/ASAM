! Berechnet die heteorogene Nukleationsrate Js fuer die Eisentstehung
! aus der Wasserdampfphase (Depositionsgefrieren).

SUBROUTINE depoxd(TT,usatt,dqq,ni,dni,qia,dqia,dnf,qf,dqf,qfa,dqfa,miv,ip)
  
  INCLUDE 'HEADER2'   
  
  IMPLICIT NONE

  INTEGER :: ip
  DOUBLE PRECISION :: TT,usatt,dqq,miv(itmax)
  ! Felder fuer unloesliche Partikel
  DOUBLE PRECISION :: ni(simax,itmax,ipmax),qia(simax,itmax,ipmax)
  DOUBLE PRECISION :: dni(simax,itmax),dqia(simax,itmax)
  ! Felder fuer Eisteilchen
  DOUBLE PRECISION :: dnf(jmax,smax)
  DOUBLE PRECISION :: qf(jmax,smax,ipmax),qfa(jmax,smax,ipmax)
  DOUBLE PRECISION :: dqf(jmax,smax),dqfa(jmax,smax)
  
  INTEGER :: si,it,jj
  DOUBLE PRECISION :: R,Mw,sigma,konst,kb
  DOUBLE PRECISION :: ag,Sw,Seis,xx,phi,ff,dF,Js(simax,itmax),dd,psi,chi
  DOUBLE PRECISION :: Vg(simax,itmax),meis(simax,itmax)
  
  ! Konstanten belegen
  R = 8.3144D0
  Mw = 0.018015D0
  sigma = 0.1D0
  konst = 1.0D+30
  kb = 1.381D-23
  
  IF(TT.GE.T0) RETURN
  

  
  ! heterogene Nukleationsrate Js berechnen
  Sw   = usatt + 1.0D0
  Seis = MIN(1.0D0,(TT/T0)**2.66D0)
  Seis = Sw/Seis
  IF(Seis.GT.1.0D0) THEN
     ag = 2.0D0*Mw*sigma/(R*TT*rhoi*DLOG(Seis))
     DO it=1,itmax
        DO si=1,simax
           xx  = rsmitte(si)/ag
           phi = DSQRT(1.0D0 - 2.0D0*miv(it)*xx + xx**2.0D0)
           ff  = 1.0D0 + ((1.0D0 - miv(it)*xx)/phi)**3.0D0 + xx**3.0D0*        &
                 (2.0D0 - 3.0D0*(xx - miv(it))/phi +                           &
                 ((xx - miv(it))/phi)**3.0D0)*                                 &
                 3.0D0*miv(it)*xx**2.0D0*((xx - miv(it))/phi - 1)
           ff  = ff/2.0D0
           dF  = 4.0D0*pi*ag**2.0D0*sigma*ff/3.0D0
           IF(dF.LT.0.0D0) dF = 1.0D-12
           !      IF(dF.GT.1.0D-16) Js(si,it) = 0.0D0
           !      IF(dF.LE.1.0D-16.AND.dF.GT.1.0D-21)                               &
           Js(si,it) = konst*rsmitte(si)**2.0D0*DEXP(-dF/(kb*TT))
           !      IF(dF.LE.1.0D-21) Js(si,it) = 100.0D0

           ! Volumen Vg des Eiskeimes berechnen
           dd  = DSQRT(rsmitte(si)**2.0D0 + ag**2.0D0 - 2.0D0*rsmitte(si)*     &
                 ag*miv(it))
           psi = -(ag - rsmitte(si)*miv(it))/dd
           chi = (rsmitte(si) - ag*miv(it))/dd
           Vg(si,it) = pi*ag**3.0D0*(2.0D0 - 3.0D0*psi + psi**3.0D0)          &
                     - pi*rsmitte(si)**3.0D0*(2.0D0 - 3.0D0*chi + chi**3.0D0)
           Vg(si,it) = Vg(si,it)/3.0D0

           ! Eismasse pro Teilchen berechnen
           meis(si,it) = rhoi*Vg(si,it)
           IF(Js(si,it).GE.1.0D0) meis(si,it) = rhoi*Vg(si,it)*Js(si,it)
        END DO
     END DO
  ELSE IF(Seis.LE.1.0D0) THEN
     Js   = 0.0D0
     Vg   = 0.0D0
     meis = 0.0D0
  END IF
  
  ! Aenderungen in den Verteilungsfunktionen der unloeslichen Partikel
  ! berechnen.
  DO it=1,itmax
     DO si=1,simax
        ! Anzahl der neu entstandenen Eisteilchen = -dni
        dni(si,it)  = -ni(si,it,ip)*Js(si,it)
        IF(Js(si,it).GE.1.0D0) dni(si,it) = -ni(si,it,ip)
        dqia(si,it) = dni(si,it)/ni(si,it,ip)*qia(si,it,ip)
        IF(ni(si,it,ip).LE.0.0D0) dqia(si,it) = 0.0D0    
     END DO
  END DO
  
  ! Aenderungen in den Verteilungsfunktionen der Eisteilchen berechnen.
  dqq = 0.0D0
  DO it=1,itmax
     DO si=1,simax
        DO jj=1,jmax
           IF(meis(si,it).GT.0.0D0.AND.meis(si,it).LT.mgrenz(1)) THEN
              ! Fuer die 1d-Variante, smax = 1.
              IF(smax.EQ.1) THEN
                 dnf(1,1) = -dni(si,it)
                 dqf(1,1) = -rhoi*Vg(si,it)*dni(si,it)
                 IF(Js(si,it).GE.1.0D0) dqf(1,1) =                         &
                                        -rhoi*Vg(si,it)*dni(si,it)*Js(si,it)
                 dqq = dqq - dqf(1,1)
                 dqfa(1,1) = -dqia(si,it)
              END IF
              ! Fuer die 2d-Variante mit smax = simax.
              IF(smax.EQ.simax) THEN
                 dnf(1,si) = -dni(si,it)
                 dqf(1,si) = -rhoi*Vg(si,it)*dni(si,it)
                 IF(Js(si,it).GE.1.0D0) dqf(1,si) =                        &
                                        -rhoi*Vg(si,it)*dni(si,it)*Js(si,it)
                 dqq = dqq - dqf(1,si)
                 dqfa(1,si) = -dqia(si,it)
              END IF
           END IF
           IF(meis(si,it).GE.mgrenz(jj).AND.meis(si,it).LT.mgrenz(jj+1)) THEN

              ! Fuer die 1d-Variante, smax = 1.
              IF(smax.EQ.1) THEN
                 dnf(jj,1) = -dni(si,it)
                 dqf(jj,1) = -rhoi*Vg(si,it)*dni(si,it)
                 IF(Js(si,it).GE.1.0D0) dqf(jj,1) =                        &
                                        -rhoi*Vg(si,it)*dni(si,it)*Js(si,it)
                 dqq = dqq - dqf(jj,1)
                 dqfa(jj,1) = -dqia(si,it)
              END IF

              ! Fuer die 2d-Variante mit smax = simax.
              IF(smax.EQ.simax) THEN
                 dnf(jj,si) = -dni(si,it)
                 dqf(jj,si) = -rhoi*Vg(si,it)*dni(si,it)
                 IF(Js(si,it).GE.1.0D0) dqf(jj,si) =                       &
                                        -rhoi*Vg(si,it)*dni(si,it)*Js(si,it)
                 dqq = dqq - dqf(jj,si)
                 dqfa(jj,si) = -dqia(si,it)
              END IF
           END IF
        END DO
     END DO
  END DO
  

  
  RETURN
END SUBROUTINE depoxd
