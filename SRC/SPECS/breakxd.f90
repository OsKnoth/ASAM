SUBROUTINE breakxd(BREAN,BREAQ,BREAS,BREAA,MQUER,NW,QW,QS,QA,ip)
  
  ! break-up process 
  ! new version
  ! called by cloud.f
  ! calling none
  
  INCLUDE 'HEADER2'
  
  IMPLICIT NONE
  
  INTEGER :: I1,I2,J1,J2,JBRE,JTAR,ip
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: DQQW(JMAX,SMAX),DQNW(JMAX,SMAX)
  DOUBLE PRECISION :: DQQS(JMAX,SMAX),DQQA(JMAX,SMAX)
  DOUBLE PRECISION :: DSQW(JMAX,SMAX),DSNW(JMAX,SMAX)
  DOUBLE PRECISION :: DSQS(JMAX,SMAX),DSQA(JMAX,SMAX)
  DOUBLE PRECISION :: BREAN(JMAX,SMAX),BREAQ(JMAX,SMAX),BREAS(JMAX,SMAX)
  DOUBLE PRECISION :: BREAA(JMAX,SMAX)
  DOUBLE PRECISION :: BREA1,BREA2,BREA3,BREA4,QSbreak
  DOUBLE PRECISION :: MQUER(JMAX,SMAX),RQUER
  DOUBLE PRECISION :: SQUER(JMAX,SMAX)
  DOUBLE PRECISION :: PBRE,QBRE,RQUOT
  DOUBLE PRECISION :: WWNW(JMAX,JMAX,SMAX),WWQW(JMAX,JMAX,SMAX)
  DOUBLE PRECISION :: WWQS(JMAX,JMAX,SMAX),WWQA(JMAX,JMAX,SMAX)
  DOUBLE PRECISION :: NORM,NORM1
  
  ! write(*,*) "Hier break.f"
  ! set to zero
  DO I1=1,SMAX
     DO J1=1,JMAX
        DQNW(J1,I1) = 0.D0
        DQQW(J1,I1) = 0.D0
        DQQS(J1,I1) = 0.D0
        DQQA(J1,I1) = 0.D0
        DSNW(J1,I1) = 0.D0
        DSQW(J1,I1) = 0.D0
        DSQS(J1,I1) = 0.D0
        DSQA(J1,I1) = 0.D0
     ENDDO
  ENDDO
  
  ! Calculation of mass distribution due to break-up
  ! Probability distribution, spontaneous break-up after Srivastava 1971
  ! Choice of the radii, masses used:
  ! Breaking drops:       mquer, rquer
  ! New drops:            mmitte, rmitte
  ! seems to be sensitive to the choice of JBRE and JTAR only if both are 
  ! chosen relatively small, e.g., 42/03
  ! For JBRE close to JMAX (>55) or JTAR gt JMAX/2 results do not differ much
  ! Limiting bin for breakup
  !      JBRE=63*JMAX/66
  !      JTAR=54*JMAX/66
  JBRE  =  60*JMAX/66
  JTAR  =  42*JMAX/66
  ! Interaction terms between all bins
  DO I1=1,SMAX
     DO J2=1,JMAX
        DO J1=1,JMAX
           WWNW(J1,J2,I1) = 0.D0
           WWQW(J1,J2,I1) = 0.D0
           WWQA(J1,J2,I1) = 0.D0
           WWQS(J1,J2,I1) = 0.D0
        ENDDO
     ENDDO
  ENDDO
  ! write(*,*) "Hier break.f 0"
  ! J1,I1: breaking drop; J2: target drop
  DO I1=1,SMAX
     DO J1=JBRE,JMAX
        IF(MQUER(J1,I1).GT.0.D0) THEN
!!!     write(*,*) "break.f: bin ",J1,I1
           SQUER(J1,I1) =  QA(J1,I1,ip)/NW(J1,I1,ip)
           RQUER        = (MQUER(J1,I1)/FACT)**(1.D0/3.D0)
           PBRE         =  PC1*EXP(PC2*RQUER)
           NORM         =  0.D0
           DO J2=JTAR,J1-1
              RQUOT     = RMITTE(J2)/RQUER
              QBRE      = QC1/MMITTE(J2) * RQUOT * EXP(-1.D0*QC2*RQUOT)
              NORM1     = PBRE*QBRE*DIFF21(J2)
              WWQW(J1,J2,I1) = NORM1*QW(J1,I1,ip)
              WWQA(J1,J2,I1) = NORM1*QA(J1,I1,ip)
              WWQS(J1,J2,I1) = NORM1*QS(J1,I1,ip)
              NORM           = NORM+NORM1
           ENDDO
           IF(NORM.GT.1.D0) THEN
              DO J2             = 1,J1-1
                 WWQW(J1,J2,I1) = WWQW(J1,J2,I1)/NORM
                 WWQA(J1,J2,I1) = WWQA(J1,J2,I1)/NORM
                 WWQS(J1,J2,I1) = WWQS(J1,J2,I1)/NORM
              ENDDO
           ENDIF
        ENDIF
     ENDDO
  ENDDO
  ! write(*,*) "Hier break.f WW-Ende"
  ! Source terms for J2,I2 from J1,I1
  DO I1=1,SMAX
     DO J1=JBRE,JMAX
        DO J2=JTAR,J1-1
           if(WWQW(J1,J2,I1).le.0.D0) goto 222
           QSbreak = WWQA(J1,J2,I1)/WWQW(J1,J2,I1)*MMITTE(J2)
           if(SMAX.gt.1) then
              DO I2=1,SMAX
                 if(QSbreak.ge.SGRENZ(I2).and.QSbreak.lt.SGRENZ(I2+1)) then
                    DQQW(J2,I2) = DQQW(J2,I2)+WWQW(J1,J2,I1)
                    DQQA(J2,I2) = DQQA(J2,I2)+WWQA(J1,J2,I1)
                    DQQS(J2,I2) = DQQS(J2,I2)+WWQS(J1,J2,I1)
                    DQNW(J2,I2) = DQQW(J2,I2)/MMITTE(J2)
                    !!write(*,*) "Quell",J2,I2,DQQA(J2,I2)/DQNW(J2,I2)/SGRENZ(I2)
                 endif
              ENDDO
           else
              I2          = 1
              DQQW(J2,I2) = DQQW(J2,I2)+WWQW(J1,J2,I1)
              DQQA(J2,I2) = DQQA(J2,I2)+WWQA(J1,J2,I1)
              DQQS(J2,I2) = DQQS(J2,I2)+WWQS(J1,J2,I1)
              DQNW(J2,I2) = DQQW(J2,I2)/MMITTE(J2)
           endif
222        continue
        ENDDO
     ENDDO
  ENDDO
  ! write(*,*) "Hier break.f Quelle-Ende"
  ! Sink terms for J1,I1
  DO I1=1,SMAX   
     DO J1=JBRE,JMAX
        IF(MQUER(J1,I1).GT.0.D0) THEN
           DO J2=JTAR,J1-1
              DSQW(J1,I1) = DSQW(J1,I1)+WWQW(J1,J2,I1) 
              DSQA(J1,I1) = DSQA(J1,I1)+WWQA(J1,J2,I1) 
              DSQS(J1,I1) = DSQS(J1,I1)+WWQS(J1,J2,I1) 
           ENDDO
           DSNW(J1,I1)    = DSQW(J1,I1)/MQUER(J1,I1)
           !      write(*,*) "Senke:",J1,I1,SQUER(J1,I1)/DSQS(J1,I1)*DSNW(J1,I1)
           !      DSQS(J1,I1)=DSNW(J1,I1)*SQUER(J1,I1)
        ENDIF
     ENDDO
  ENDDO
  ! write(*,*) "Hier break.f Senke-Ende"
  ! Total balance
  DO I1=1,SMAX
     DO J1=1,JMAX
        BREAN(J1,I1) = DQNW(J1,I1)-DSNW(J1,I1)
        BREAQ(J1,I1) = DQQW(J1,I1)-DSQW(J1,I1)
        BREAA(J1,I1) = DQQA(J1,I1)-DSQA(J1,I1)
        BREAS(J1,I1) = DQQS(J1,I1)-DSQS(J1,I1)
        IF(ABS(BREAQ(J1,I1)).GT.1.D0)                              &
             WRITE(*,*) "break.f6",J1,I1,BREAQ(J1,I1),DQQW(J1,I1),DSQW(J1,I1)
        IF(BREAQ(J1,I1).GT.1.D0.OR.BREAQ(J1,I1).LE.1.D0) THEN
        ELSE
           WRITE(*,*) "break.f7",J1,I1,BREAQ(J1,I1),DQQW(J1,I1),DSQW(J1,I1)
        ENDIF
        IF(ABS(BREAS(J1,I1)).GT.1.D0) THEN 
           WRITE(*,*) "break.f8",J1,I1,BREAS(J1,I1),DQQS(J1,I1),DSQS(J1,I1)
        ENDIF
     ENDDO
  ENDDO
  ! Mass balance
  BREA1 = 0.D0
  BREA2 = 0.D0
  BREA3 = 0.D0
  BREA4 = 0.D0
  DO I1=1,SMAX
     DO J1=1,JMAX
        BREA1 = BREA1+BREAQ(J1,I1)
        BREA2 = BREA2+BREAN(J1,I1)
        BREA3 = BREA3+BREAA(J1,I1)
        BREA4 = BREA4+BREAS(J1,I1)
     ENDDO
  ENDDO
  IF(BREA1.GT.1.D-10) WRITE(*,*) "break: Massenbilanz: ",    BREA1
  IF(BREA2.LT.0.D0) WRITE(*,*) "break: Anzahlbilanz: ",      BREA2
  IF(BREA3.GT.1.D-10) WRITE(*,*) "break: Massenbilanz Lsg: ",BREA3
  IF(BREA4.GT.1.D-10) WRITE(*,*) "break: Massenbilanz Lsg: ",BREA4
  !      WRITE(*,*) "break: Massenbilanz: ",BREA1
  !      WRITE(*,*) "break: Anzahlbilanz: ",BREA2
  !      WRITE(*,*) "break: Massenbilanz Lsg: ",BREA3
  !      write(*,*) "Hier break.f Bilanz-Ende"
  
  RETURN
END SUBROUTINE breakxd
