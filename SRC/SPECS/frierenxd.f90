! Gefrieren von Eisteilchen mit Wassermantel nach PK 97, S. 668 ff.

SUBROUTINE frierenxd(TT,nf,qf,dqf,qfw,rfquer,ip,i_lm,j_lm,k_lm,t_lm)
  
  INCLUDE 'HEADER2' 
  
USE data_parallel,      ONLY: my_cart_id, nboundlines, isubpos

  IMPLICIT NONE
  
  INTEGER :: ip, i_lm, j_lm, k_lm ,t_lm
  DOUBLE PRECISION :: TT
  ! Felder fuer Eisteilchen
  DOUBLE PRECISION :: nf(jmax,smax,ipmax)
  DOUBLE PRECISION :: qf(jmax,smax,ipmax),qfw(jmax,smax,ipmax)
  DOUBLE PRECISION :: dqf(jmax,smax)
  DOUBLE PRECISION :: rfquer(jmax,smax)
  
  INTEGER :: j,i
  DOUBLE PRECISION :: meis,dmeis,mwasser,Gm,rwasser,reis,da,Lm
  PARAMETER(Lm = 0.33D+6)
  
  ! Funktionen
  DOUBLE PRECISION :: ceisxd,cwasserxd,keisxd
  
  
!!$IF(TT.LE.T0) THEN
!!$! Gefrierrate nach PK 97, S. 668
!!$  Gm = 3.5D-04*(T0 - TT)**2.22
!!$  DO j=1,jmax
!!$    DO i=1,smax
!!$      dqf(j,i) = 0.0D0
!!$      mwasser = (qfw(j,i,ip) - qf(j,i,ip))/nf(j,i,ip)
!!$      IF(mwasser.GT.0.0D0) THEN
!!$        meis = qf(j,i,ip)/nf(j,i,ip)
!!$        reis = (3.0D0*meis/(4.0D0*pi*rhoi))**(1.0D0/3.0D0)
!!$        rwasser = rfquer(j,i) - reis
!!$        IF(rwasser.LE.0.0D0) rwasser = 0.0D0
!!$        da = Gm*deltat
!!$        IF(da.GE.rwasser) dqf(j,i) = qfw(j,i,ip) - qf(j,i,ip)
!!$        IF(da.LT.rwasser) dqf(j,i) = 4.0D0*pi*da**3.0D0*rhoi*nf(j,i,ip)/3.0D0
!!$      END IF
!!$    END DO
!!$  END DO
!!$END IF

!!$! Gefrierrate nach PK 97, S. 676, Gl. 16.26 umgestellt
!!$IF(TT.LT.T0) THEN
!!$  DO j=1,jmax
!!$    DO i=1,smax
!!$      dqf(j,i) = 0.0D0
!!$      mwasser = 0.0D0
!!$      meis = 0.0D0
!!$      IF(nf(j,i,ip).GT.1.0D-15) THEN
!!$        mwasser = (qfw(j,i,ip) - qf(j,i,ip))/nf(j,i,ip)
!!$        meis = qf(j,i,ip)/nf(j,i,ip)
!!$      END IF
!!$      IF(mwasser.GT.0.0D0.AND.meis.GT.0.0D0) THEN
!!$        reis = (3.0D0*meis/(4.0D0*pi*rhoi))**(1.0D0/3.0D0)
!!$        rwasser = rfquer(j,i) - reis
!!$        IF(rwasser.LE.0.0D0) rwasser = 0.0D0
!!$! VG Formel korrigiert (rhow quadriert). Ist sie ueberhaupt verwendbar?
!!$        da = SQRT(4.0D0*deltat*keisxd(TT)*rhoi*ceisxd(TT)*(T0 - TT)**2.0D0   &
!!$            /(pi*rhow*rhow*(Lm - cwasserxd(TT)*(T0 - TT))**2.0D0))
!!$        IF(da.GE.rwasser) dqf(j,i) = qfw(j,i,ip) - qf(j,i,ip)
!!$        IF(da.LT.rwasser) dqf(j,i) = 4.0D0*pi*da**3.0D0*rhoi*nf(j,i,ip)/3.0D0
!!$      END IF
!!$    END DO
!!$  END DO
!!$END IF
  
  
  IF(TT.LT.T0) THEN
     DO j=1,jmax
        DO i=1,smax
           dqf(j,i) = 0.0D0
           mwasser = 0.0D0
           meis = 0.0D0
!           IF(nf(j,i,ip).GT.1.0D-15) THEN
           IF(nf(j,i,ip).GT.SMALL1) THEN
              mwasser = (qfw(j,i,ip) - qf(j,i,ip))/nf(j,i,ip)
              meis    =  qf(j,i,ip)/nf(j,i,ip)
           END IF
           IF(mwasser.GT.0.0D0.AND.meis.GT.0.0D0) THEN
              reis  = (3.0D0*meis/(4.0D0*pi*rhoi))**(1.0D0/3.0D0)
              da    = SQRT(4.0D0*deltat*keisxd(TT)*rhoi*ceisxd(TT)*(T0 - TT)**2.0D0   &
                           /(pi*rhow**2*(Lm - cwasserxd(TT)*(T0 - TT))**2.0D0))
              dmeis = 4.D0*pi *rhoi * (reis + da)**3.D0 / 3.D0 - meis
              IF (dmeis .LT. 0.D0) print *,'frierenxd.f: dmeis < 0',t_lm, i_lm, j_lm, k_lm
              IF(dmeis.GE.mwasser) THEN 
                 dqf(j,i) = (qfw(j,i,ip) - qf(j,i,ip)) / deltat
              ELSE
                 dqf(j,i) = dmeis * nf(j,i,ip) / deltat
                 IF (dqf(j,i) .gt. (qfw(j,i,ip) - qf(j,i,ip))/deltat) &
                      print *,'frierenxd.f: dqf > than possible!',t_lm, i_lm, j_lm, k_lm, dqf(j,i),qfw(j,i,ip) - qf(j,i,ip)
              END IF
           ENDIF
        END DO
     END DO
  END IF
  
  RETURN
END SUBROUTINE frierenxd
