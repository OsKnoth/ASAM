! Schmelzen von Eisteilchen nach PK 97, S. 692 ff.

SUBROUTINE schmelzenxd2(TT,pp,rho,usatt,dqq,vt,nf,dnf,qf,dqf,qfw,dqfw,    &
                        qfs,dqfs,qfa,dqfa,dnw,dqw,dqws,dqwa,ip,rfquer)

INCLUDE 'HEADER3'


IMPLICIT NONE

INTEGER :: ip
DOUBLE PRECISION :: TT,pp,rho,usatt,dqq,vt(jmax)
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
DOUBLE PRECISION :: Mw,Rgas,dmi,meis,dmw,mwasser,satt
PARAMETER(Mw = 0.018015D0, Rgas = 8.3144D0)

! Functions
DOUBLE PRECISION :: kluftxd,ventixd,esattxd,dvxd


IF(TT.GT.T0) THEN
  satt = usatt + 1.0D0
  dqq = 0.0D0
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
      dmw = 0.0D0
! Schmelzen des Eiskerns falls vorhanden
      IF(nf(j,i,ip).GT.1.0D-20.AND.qf(j,i,ip).GT.0.0D0) THEN
        mwasser = (qfw(j,i,ip) - qf(j,i,ip))/nf(j,i,ip)
        meis = qf(j,i,ip)/nf(j,i,ip)
        dmi = -4.0D0*pi*rmitte(j)*kluftxd(TT)*(TT - T0)                &
              *ventixd(TT,pp,rho,rmitte(j),vt(j))/LV_ICE
        IF(DABS(dmi).GE.meis) dmi = -meis

! Verdunsten der Wasserhuelle falls vorhanden
        IF(mwasser.GT.0.0D0) THEN
          dmw = -4.0D0*pi*rmitte(j)*dvxd(TT,pp)*LV*Mw                      &
                *ventixd(TT,pp,rho,rmitte(j),vt(j))*(satt*esattxd(TT)/TT   &
              - esattxd(T0)/T0)/(Rgas*LV_ICE)
          IF(dmw.LT.0.0D0.AND.DABS(dmw).GE.mwasser) dmw = -mwasser
!!          IF(dmw.GT.0.0D0) dmw = 0.0D0
        ELSE
          dmw = 0.0D0
        END IF

! Falls alles schmilzt und verdunstet
        IF(DABS(dmi).GE.meis.AND.DABS(dmw).GE.mwasser) THEN
          dnf(j,i) = -nf(j,i,ip)
          dqf(j,i) = -qf(j,i,ip)
          dqfw(j,i) = -qfw(j,i,ip)
          dqfs(j,i) = -qfs(j,i,ip)
          dqfa(j,i) = - qfa(j,i,ip)
        END IF

! Falls alles schmilzt und noch Wasser uebrig bleibt
        IF(DABS(dmi).GE.meis.AND.DABS(dmw).LT.mwasser) THEN
          dnf(j,i) = -nf(j,i,ip)
          dqf(j,i) = -qf(j,i,ip)
          dqfw(j,i) = -qfw(j,i,ip)
          dqfs(j,i) = -qfs(j,i,ip)
          dqfa(j,i) = - qfa(j,i,ip)
          dnw(j,i) = nf(j,i,ip)
          dqw(j,i) = (mwasser + dmw - dmi)*nf(j,i,ip)
          dqws(j,i) = dnw(j,i)/nf(j,i,ip)*qfs(j,i,ip)
          dqwa(j,i) = dnw(j,i)/nf(j,i,ip)*qfa(j,i,ip)
        END IF

! Falls nicht alles schmilzt
        IF(DABS(dmi).LT.meis.AND.dmi.NE.0.0D0) THEN
          dqf(j,i) = (dmi + dmw)*nf(j,i,ip)
          dqfw(j,i) = dmw*nf(j,i,ip)
        END IF

        dqq = dqq - (dqfw(j,i) + dqw(j,i))

      END IF
      
    END DO
  END DO
END IF

           
RETURN
END
