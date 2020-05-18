!*****************************************************************
!*                                                               *
!*    subroutine  aktuellxd                                      *
!*                                                               *
!*****************************************************************
SUBROUTINE aktuellxd(nw,qw,qws,qwa,dnw,dqw,dqws,dqwa,dmnw,dmqw,dmqws,       &
                     dmqwa,nf,qf,qfs,qfa,qfw,dnf,dqf,dqfs,dqfa,dqfw,dmnf,   &
                     dmqf,dmqfs,dmqfa,dmqfw,ni,qia,dni,dqia,dmni,dmqia,qq,  &
                     dqq,dmqq,teta,dteta,dmteta,kz,dt)

  INCLUDE 'HEADER2'
  
  IMPLICIT NONE
  
  INTEGER :: kz
  DOUBLE PRECISION :: dt
  DOUBLE PRECISION :: nw(kz,jmax,smax,ipmax),qw(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: qws(kz,jmax,smax,ipmax),qwa(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: nf(kz,jmax,smax,ipmax),qf(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: qfs(kz,jmax,smax,ipmax),qfa(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: qfw(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: ni(kz,simax,itmax,ipmax),qia(kz,simax,itmax,ipmax)
  DOUBLE PRECISION :: qq(kz,ipmax),dqq(kz,ipmax),dmqq(kz,ipmax)
  DOUBLE PRECISION :: teta(kz,ipmax),dteta(kz,ipmax),dmteta(kz,ipmax)
  
  DOUBLE PRECISION :: dnw(kz,jmax,smax,ipmax),dqw(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dqws(kz,jmax,smax,ipmax),dqwa(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dnf(kz,jmax,smax,ipmax),dqf(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dqfs(kz,jmax,smax,ipmax),dqfa(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dqfw(kz,jmax,smax,ipmax)      
  DOUBLE PRECISION :: dmnw(kz,jmax,smax,ipmax),dmqw(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dmqws(kz,jmax,smax,ipmax),dmqwa(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dmnf(kz,jmax,smax,ipmax),dmqf(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dmqfs(kz,jmax,smax,ipmax),dmqfa(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dmqfw(kz,jmax,smax,ipmax)
  DOUBLE PRECISION :: dni(kz,simax,itmax,ipmax),dqia(kz,simax,itmax,ipmax)
  DOUBLE PRECISION :: dmni(kz,simax,itmax,ipmax),dmqia(kz,simax,itmax,ipmax)
  
  INTEGER :: i,j,k,si,it,ip
  DOUBLE PRECISION :: nwalt,qwalt,qwsalt,qwaalt
  DOUBLE PRECISION :: nfalt,qfalt,qfsalt,qfaalt,qfwalt
  DOUBLE PRECISION :: nialt,qiaalt
  
  DOUBLE PRECISION :: small_number
      
  small_number = 1.e-10

DO ip=1,ipmax
  DO k=1,kz
! Aktualisierung der Wasser- und Eis-Verteilungen
    DO j=1,jmax
      DO i=1,smax
        nwalt = nw(k,j,i,ip)
        qwalt = qw(k,j,i,ip)
        qwsalt = qws(k,j,i,ip)
        qwaalt = qwa(k,j,i,ip)

        nfalt = nf(k,j,i,ip)
        qfalt = qf(k,j,i,ip)
        qfsalt = qfs(k,j,i,ip)
        qfaalt = qfa(k,j,i,ip)
        qfwalt = qfw(k,j,i,ip)

        nw(k,j,i,ip) = nw(k,j,i,ip) + (dnw(k,j,i,ip) + dmnw(k,j,i,ip))*dt
        qw(k,j,i,ip) = qw(k,j,i,ip) + (dqw(k,j,i,ip) + dmqw(k,j,i,ip))*dt
        qws(k,j,i,ip) = qws(k,j,i,ip) + (dqws(k,j,i,ip) + dmqws(k,j,i,ip))*dt
        qwa(k,j,i,ip) = qwa(k,j,i,ip) + (dqwa(k,j,i,ip) + dmqwa(k,j,i,ip))*dt

        nf(k,j,i,ip) = nf(k,j,i,ip) + (dnf(k,j,i,ip) + dmnf(k,j,i,ip))*dt
        qf(k,j,i,ip) = qf(k,j,i,ip) + (dqf(k,j,i,ip) + dmqf(k,j,i,ip))*dt
        qfs(k,j,i,ip) = qfs(k,j,i,ip) + (dqfs(k,j,i,ip) + dmqfs(k,j,i,ip))*dt
        qfa(k,j,i,ip) = qfa(k,j,i,ip) + (dqfa(k,j,i,ip) + dmqfa(k,j,i,ip))*dt
        qfw(k,j,i,ip) = qfw(k,j,i,ip) + (dqfw(k,j,i,ip) + dmqfw(k,j,i,ip))*dt


        IF (nwalt.gt.0.D0) THEN 
           IF(ABS(nw(k,j,i,ip)/nwalt).LT.1.0D-15.OR.                   &
                ABS(qw(k,j,i,ip)/qwalt).LT.1.0D-15.OR.                   &
                ABS(qws(k,j,i,ip)/qwsalt).LT.1.0D-15.OR.                 &
                ABS(qwa(k,j,i,ip)/qwaalt).LT.1.0D-15) THEN
              nw(k,j,i,ip) = 0.0D0
              qw(k,j,i,ip) = 0.0D0
              qws(k,j,i,ip) = 0.0D0
              qwa(k,j,i,ip) = 0.0D0
           END IF
        ENDIF

        IF (nfalt.gt.0.D0) THEN 
           IF(ABS(nf(k,j,i,ip)/nfalt).LT.1.0D-15.OR.                   &
                ABS(qf(k,j,i,ip)/qfalt).LT.1.0D-15.OR.                   &
                ABS(qfw(k,j,i,ip)/qfsalt).LT.1.0D-15.OR.                 &
                ABS(qfs(k,j,i,ip)/qfsalt).LT.1.0D-15.OR.                 &
                ABS(qfa(k,j,i,ip)/qfaalt).LT.1.0D-15) THEN
              nf(k,j,i,ip) = 0.0D0
              qf(k,j,i,ip) = 0.0D0
              qfw(k,j,i,ip) = 0.0D0
              qfs(k,j,i,ip) = 0.0D0
              qfa(k,j,i,ip) = 0.0D0
           END IF
        ENDIF
        IF(nw(k,j,i,ip).LT.small_number.OR.qw(k,j,i,ip).LE.0.0D0.OR.     &
           qws(k,j,i,ip).LE.0.0D0.OR.qwa(k,j,i,ip).LE.0.0D0) THEN
          nw(k,j,i,ip) = 0.0D0
          qw(k,j,i,ip) = 0.0D0
          qws(k,j,i,ip) = 0.0D0
          qwa(k,j,i,ip) = 0.0D0
        END IF

        IF(nf(k,j,i,ip).LT.small_number.OR.qf(k,j,i,ip).LE.0.0D0.OR.     &
           qfw(k,j,i,ip).LE.0.0D0.OR.                               &
           qfs(k,j,i,ip).LE.0.0D0.OR.qfa(k,j,i,ip).LE.0.0D0) THEN
          nf(k,j,i,ip) = 0.0D0
          qf(k,j,i,ip) = 0.0D0
          qfw(k,j,i,ip) = 0.0D0
          qfs(k,j,i,ip) = 0.0D0
          qfa(k,j,i,ip) = 0.0D0
        END IF

        IF(qws(k,j,i,ip).GT.qwa(k,j,i,ip)) qws(k,j,i,ip) = qwa(k,j,i,ip)  !MS

!!        IF(nw(k,j,i,ip).LT.1.0D-20.OR.qw(k,j,i,ip).LE.0.0D0) THEN
!!          nw(k,j,i,ip) = 0.0D0
!!          qw(k,j,i,ip) = 0.0D0
!!          qws(k,j,i,ip) = 0.0D0
!!          qwa(k,j,i,ip) = 0.0D0
!!        END IF
!!        IF(qws(k,j,i,ip).LT.0.0D0) qws(k,j,i,ip) = 0.0D0
!!        IF(qwa(k,j,i,ip).LT.0.0D0) qwa(k,j,i,ip) = 0.0D0

!!        IF(nf(k,j,i,ip).LT.1.0D-20.OR.qfw(k,j,i,ip).LE.0.0D0) THEN
!!          nf(k,j,i,ip) = 0.0D0
!!          qf(k,j,i,ip) = 0.0D0
!!          qfs(k,j,i,ip) = 0.0D0
!!          qfa(k,j,i,ip) = 0.0D0
!!          qfw(k,j,i,ip) = 0.0D0
!!        END IF
! MS: 'Numerical melting' of mixed particles
        IF(qf(k,j,i,ip).LE.0.0D0.AND.qfw(k,j,i,ip).GT.0.0D0) THEN
          nw(k,j,i,ip) = nw(k,j,i,ip) + nf(k,j,i,ip)
          qw(k,j,i,ip) = qw(k,j,i,ip) + qfw(k,j,i,ip)
          qws(k,j,i,ip) = qws(k,j,i,ip) + qfs(k,j,i,ip)
          qwa(k,j,i,ip) = qwa(k,j,i,ip) + qfa(k,j,i,ip)
          nf(k,j,i,ip) = 0.0D0
          qf(k,j,i,ip) = 0.0D0
          qfw(k,j,i,ip) = 0.0D0
          qfs(k,j,i,ip) = 0.0D0
          qfa(k,j,i,ip) = 0.0D0
        END IF 
! MS: limitation of frozen mass compared to total mass
        IF(qf(k,j,i,ip).GT.qfw(k,j,i,ip)) THEN
          qf(k,j,i,ip) = qfw(k,j,i,ip)
        END IF 
!!        IF(qfs(k,j,i,ip).LT.0.0D0) qfs(k,j,i,ip) = 0.0D0
        IF(qfs(k,j,i,ip).GT.qfa(k,j,i,ip)) qfs(k,j,i,ip) = qfa(k,j,i,ip)  !MS
!!        IF(qfa(k,j,i,ip).LT.0.0D0) qfa(k,j,i,ip) = 0.0D0

      END DO
    END DO
! Aktualisierung der unloeslichen Verteilungen
    DO si=1,simax
      DO it=1,itmax
        nialt = ni(k,si,it,ip)
        qiaalt = qia(k,si,it,ip)
        ni(k,si,it,ip) = ni(k,si,it,ip) + (dni(k,si,it,ip)          &
                       + dmni(k,si,it,ip))*dt
        qia(k,si,it,ip) = qia(k,si,it,ip) + (dqia(k,si,it,ip)          &
                       + dmqia(k,si,it,ip))*dt
        
!!        IF(ABS(ni(k,si,it,ip)/nialt).LT.1.0D-15.OR.                   &
!!          ABS(qia(k,si,it,ip)/qiaalt).LT.1.0D-15) THEN
!!          ni(k,si,it,ip) = 0.0D0
!!          qia(k,si,it,ip) = 0.0D0
!!        END IF

        IF(ni(k,si,it,ip).LT.1.0D-20) THEN
          ni(k,si,it,ip) = 0.0D0
          qia(k,si,it,ip) = 0.0D0
        END IF
        IF(qia(k,si,it,ip).LT.0.0D0) qia(k,si,it,ip) = 0.0D0

      END DO
    END DO
! Aktualisierung der thermodynamischen Groessen
    teta(k,ip) = teta(k,ip) + (dteta(k,ip) + dmteta(k,ip))*dt
    qq(k,ip) = qq(k,ip) + (dqq(k,ip) + dmqq(k,ip))*dt
    IF(qq(k,ip).LT.0.0D0) qq(k,ip) = 0.0D0
  END DO
END DO

           
RETURN
END
