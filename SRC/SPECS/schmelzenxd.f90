! Schmelzen von Eisteilchen nach PK 97, S. 692 ff.
! Nur Schmelzen, keine Wechselwirkung mit der Gasphase!
! wird nur fuer T > T0, also oberhalb des Gefrierpunkts ausgefuehrt

SUBROUTINE schmelzenxd(TT,pp,rho,vt,nf,dnf,qf,dqf,qfw,dqfw,    &
                       qfs,dqfs,qfa,dqfa,dnw,dqw,dqws,dqwa,ip,rfquer)

INCLUDE 'HEADER2'

IMPLICIT NONE

INTEGER :: ip
DOUBLE PRECISION :: TT,pp,rho,vt(jmax)
! Felder fuer Eisteilchen
DOUBLE PRECISION :: nf(jmax,smax,ipmax),dnf(jmax,smax)
DOUBLE PRECISION :: qf(jmax,smax,ipmax),qfw(jmax,smax,ipmax)
DOUBLE PRECISION :: dqf(jmax,smax),dqfw(jmax,smax),qfs(jmax,smax,ipmax)
DOUBLE PRECISION :: dqfs(jmax,smax),qfa(jmax,smax,ipmax),dqfa(jmax,smax)
DOUBLE PRECISION :: rfquer(jmax,smax)
! Felder fuer Wassertropfen
DOUBLE PRECISION :: dnw(jmax,smax),dqw(jmax,smax),dqws(jmax,smax)
DOUBLE PRECISION :: dqwa(jmax,smax) 

INTEGER :: j,i
DOUBLE PRECISION :: Mw,Rgas,dmi,meis,mwasser,deis
PARAMETER(Mw = 0.018015D0, Rgas = 8.3144D0)

! Functions
DOUBLE PRECISION :: kluftxd,ventixd,dvxd


IF(TT.GT.T0) THEN
  DO j=1,jmax
    DO i=1,smax
      dnf(j,i) = 0.0D0
      dqf(j,i) = 0.0D0
      dqfw(j,i) = 0.0D0
      dqfs(j,i) = 0.0D0
      dqfa(j,i) = 0.0D0
      dnw(j,i) = 0.0D0
      dqw(j,i) = 0.0D0
      dqws(j,i) = 0.0D0
      dqwa(j,i) = 0.0D0
      dmi = 0.0D0
      deis = 0.0D0
! Schmelzen des Eiskerns falls vorhanden
      IF(nf(j,i,ip).GT.1.0D-20.AND.qf(j,i,ip).GT.0.0D0) THEN
        mwasser = (qfw(j,i,ip) - qf(j,i,ip))/nf(j,i,ip)
        meis = qf(j,i,ip)/nf(j,i,ip)

! Was ist dmi? Warum tritt rmitte hier auf (anstatt rquer oder riquer)?
! Die Masse dessen was schmilzt
       dmi = -4.0D0*pi*rmitte(j)*kluftxd(TT)*(TT - T0)                &
              *ventixd(TT,pp,rho,rmitte(j),vt(j))/LV_ICE

! Verdunsten der Wasserhuelle ist ausgeschaltet, nach cond_mixxd.f90 verlagert
        deis = dmi * deltat

! Falls Eis schmilzt 
        IF(deis.LT.0.0D0) THEN

! Falls alles schmilzt und als Wasser uebrig bleibt
          IF(DABS(deis).GE.meis) THEN
! Verlust in der Eisklasse
            dnf(j,i) = -nf(j,i,ip)/deltat
            dqf(j,i) = -qf(j,i,ip)/deltat
            dqfw(j,i) = -qfw(j,i,ip)/deltat
            dqfs(j,i) = -qfs(j,i,ip)/deltat
            dqfa(j,i) = - qfa(j,i,ip)/deltat
! Zugewinn in der entsprechenden Fluessigwasserklasse
            dnw(j,i) = nf(j,i,ip)/deltat
            dqw(j,i) = qfw(j,i,ip)/deltat
            dqws(j,i) = qfs(j,i,ip)/deltat
            dqwa(j,i) = qfa(j,i,ip)/deltat

! Falls nicht alles schmilzt: Umschichtung von Eis zu Wasserhuelle, d.h.
! Gesamtmasse bleibt gleich, Eis wird weniger
          ELSE
            dqf(j,i) = nf(j,i,ip)*deis/deltat
!!!            dqfw(j,i) = -nf(j,i,ip)*deis/deltat
          END IF

        END IF

! keine Wasserdampfbilanz noetig!

     END IF
      
    END DO
  END DO
END IF

           
RETURN
END
