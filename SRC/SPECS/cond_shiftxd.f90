SUBROUTINE COND_SHIFTxd(QW,QWF,NW,QS,QA,CONDQ,CONDQF,CONDN,CONDS,CONDA,   &
                         MQUER,SQUER,IP)
  
  ! calculating the shift due to condensation/evaporation/sublimation
  ! for liquid or frozen drops/particles
  ! formerly included in cond_mix.f
  
  ! called by cond_mix.f
  !           cond_e.f 
  ! calling none
  
  INCLUDE 'HEADER2' 
  
  IMPLICIT NONE
  
  INTEGER :: II,JJ,IIM,IIg,JJg,ITIM,ISTEP,IP
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QWF(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: MQUER(JMAX,SMAX),MQneu,MQalt,RQalt
  DOUBLE PRECISION :: SQUER(JMAX,SMAX),SQneu,SQalt
  DOUBLE PRECISION :: CONDQ(JMAX,SMAX),CONDN(JMAX,SMAX)
  DOUBLE PRECISION :: CONDS(JMAX,SMAX),CONDA(JMAX,SMAX)
  DOUBLE PRECISION :: CONDQF(JMAX,SMAX)
  DOUBLE PRECISION :: QWver(JMAX,SMAX),NWver(JMAX,SMAX)
  DOUBLE PRECISION :: QSver(JMAX,SMAX),QAver(JMAX,SMAX)
  DOUBLE PRECISION :: QWFver(JMAX,SMAX)
  DOUBLE PRECISION :: QWneu,QWFneu,NWneu,QSneu,QAneu,MNEU,MQUERNEU,QWver1,NWver1
  DOUBLE PRECISION :: DELM,DELM0,DELM1,DELM2,DELM21,STEIG,STEIGSTAR
  DOUBLE PRECISION :: N0,N1,N2,N0KX0,OG,UG,MSTAR
  DOUBLE PRECISION :: UGQS,QSver1,QAver1
  DOUBLE PRECISION :: N0QS,N1QS,N2QS,PP,STEIGQS,X0,STEIGSTARQS,XSTAR,OGABS
  DOUBLE PRECISION :: BIL,Mug,Mog,Sug,Sog,DELM_MIT
  DOUBLE PRECISION :: CCC1,CCC2,NNN1
  DOUBLE PRECISION :: QSver_fact
  
  INTEGER :: IMAX,IMIN,ifact
    
  !      write(*,*) "cond_shift.f: start"
  
  ! now Linear Discrete Method for the correct bin classification
  
  ! shifts are set to zero for LDM
  DO II=1,SMAX
     DO JJ=1,JMAX
        QWver(JJ,II)  = 0.D0
        QWFver(JJ,II) = 0.D0
        NWver(JJ,II)  = 0.D0
        QSver(JJ,II)  = 0.D0
        QAver(JJ,II)  = 0.D0
     ENDDO
  ENDDO
  
  ! correct mass classification
  DO II=1,SMAX 
     DO JJ=1,JMAX 
        ifact=0
        ! LDM
        IF(CONDQ(JJ,II).ne.0.D0) THEN
           QWneu  = QW(JJ,II,IP)+CONDQ(JJ,II)*DELTAT
           QWFneu = QWF(JJ,II,IP)+CONDQF(JJ,II)*DELTAT
           NWneu  = NW(JJ,II,IP)
           QSneu  = QS(JJ,II,IP)
           QAneu  = QA(JJ,II,IP)
           ! shift of the mean value only if NWNEU and QWNEU GT 0
           MNEU     = QWneu
           MQUERNEU = MNEU/NWneu
           MQalt    = MQUER(JJ,II)
           SQalt    = SQUER(JJ,II)
           DELM     = MQUERNEU-MQalt
           ! calculation of the shifted boundaries
           ! new version using interpolated growth rates
           ! growth rate of left boundary
           !          goto 777
           if(JJ.gt.1) then
              if(NW(JJ-1,II,IP).gt.1.D-30) then
                 CCC1  = QW(JJ-1,II,IP)/NW(JJ-1,II,IP)/MGRENZ(JJ)
                 CCC2  = MGRENZ(JJ)/MQUER(JJ,II)
                 NNN1  = CCC1+CCC2
                 CCC1  = CCC1/NNN1
                 CCC2  = CCC2/NNN1
                 DELM1 = MGRENZ(JJ)+(CCC1*CONDQ(JJ-1,II)/NW(JJ-1,II,IP)     &
                                   + CCC2*CONDQ(JJ,II)/NW(JJ,II,IP))*DELTAT
              else
                 DELM1 = MGRENZ(JJ)+DELM
              endif
           else
              DELM1 = MGRENZ(JJ)+DELM
           endif
!           DELM1 = MGRENZ(JJ)+DELM
           ! growth rate of right boundary
           !          goto 777
           if(JJ.lt.JMAX) then
              if(NW(JJ+1,II,IP).gt.1.D-30) then
                 CCC1  = MQUER(JJ,II)/MGRENZ(JJ+1)
                 CCC2  = MGRENZ(JJ+1)/QW(JJ+1,II,IP)*NW(JJ+1,II,IP)
                 NNN1  = CCC1+CCC2
                 CCC1  = CCC1/NNN1
                 CCC2  = CCC2/NNN1
                 DELM2 = MGRENZ(JJ+1)+(CCC1*CONDQ(JJ,II)/NW(JJ,II,IP)       &
                      + CCC2*CONDQ(JJ+1,II)/NW(JJ+1,II,IP))*DELTAT
              else
                 DELM2 = MGRENZ(JJ+1)+DELM
              endif
           else
              DELM2 = MGRENZ(JJ+1)+DELM
           endif
!           DELM2 = MGRENZ(JJ+1)+DELM
777        continue
           ! test of shifted values
           goto 321
           !          write(*,*) "cond_shift.f: DELM1/2",DELM1,MQUERNEU,DELM2
           if(MQUERNEU.LT.DELM1) write(*,*) JJ,II," cond_shift.f: too small"
           if(MQUERNEU.GT.DELM2) write(*,*) JJ,II," cond_shift.f: too large"
           if(MQUERNEU.LT.DELM1.or.MQUERNEU.GT.DELM2) then
              write(*,*) "M ",DELM1,MQUER(JJ,II),DELM2,MQUERNEU
           endif
           !          write(*,*) "cond_shift.f: shift"
           !     &        ,DELM1-Mug,DELM,DELM2-Mog
321        continue
           
           ! correction
           !          DELM1=MIN(DELM1,.999D0*MQUERNEU)
           DELM1  = MIN(DELM1,.9D0*MQUERNEU)
           DELM1  = MAX(DELM1,.0D0)
           !          DELM2=MAX(DELM2,1.001D0*MQUERNEU)
           DELM2  = MAX(DELM2,1.1D0*MQUERNEU)
           DELM0  = (DELM1+DELM2)/2.D0
           DELM21 = DELM2-DELM1
!!!          if(DELM21.le.0.D0) 
!!!     &      write(*,*) "Alarm DELM21",DELM21,DELM0,MQUERNEU
           ! Calculation of the parameters of the linear function
           ! after the shift
           N0     = NWneu/DELM21
           STEIG  = 12.D0*(MNEU-DELM0*NWneu)/DELM21**3
           N0KX0  = N0-STEIG*DELM0
           N1     = N0+STEIG*(DELM1-DELM0)
           N2     = N0+STEIG*(DELM2-DELM0)
           ! where is the shifted interval located?
           ! QWver, NWver, QSver, QAver of the target intervals are calculated
           ! N1 > 0 und N2 > 0
           !            write(*,*) "cond_shift.f: parameters",DELM21
           IF(N1.GE.0.D0.AND.N2.GE.0.D0) THEN
              IMIN = 1
              IMAX = JMAX
              DO IIM=1,JMAX
                 IF(DELM2.LE.MGRENZ(IIM+1).AND.DELM2.GT.MGRENZ(IIM)) IMAX=IIM
                 IF(DELM1.LE.MGRENZ(IIM+1).AND.DELM1.GT.MGRENZ(IIM)) IMIN=IIM
              ENDDO
              IF(DELM1.LT.MGRENZ(1)) IMIN=1
              IF(DELM2.LT.MGRENZ(1)) IMAX=1
              IF(DELM1.GT.MGRENZ(JMAX+1)) IMIN=JMAX
              IF(DELM2.GT.MGRENZ(JMAX+1)) IMAX=JMAX
              IF(IMIN.EQ.IMAX.OR.DELM.EQ.0.D0) THEN
                 QWver(IMIN,II)  = QWver(IMIN,II)+QWneu
                 QWFver(IMIN,II) = QWFver(IMIN,II)+QWFneu
                 NWver(IMIN,II)  = NWver(IMIN,II)+NWneu
                 QSver(IMIN,II)  = QSver(IMIN,II)+QSneu
                 QAver(IMIN,II)  = QAver(IMIN,II)+QAneu
                 ! Test of QSver/QAver
!!!            if(QSver(IMIN,II)/QAver(IMIN,II).lt.0.D0.or.
!!!     &         QSver(IMIN,II)/QAver(IMIN,II).gt.0.9D0)
!!!     &      print *,'cond.f LDM1a: J,I,QS/Aver(J,I)'
!!!     &      ,IMIN,II,QSver(IMIN,II),QAver(IMIN,II),NWneu
                 BIL = 1.D2
              ELSE
                 DO IIM=IMIN,IMAX
                    IF(IIM.EQ.IMIN) THEN
                       UG = DELM1
                       OG = MGRENZ(IIM+1)
                    ELSE
                       IF(IIM.EQ.IMAX) THEN
                          UG = MGRENZ(IIM)
                          OG = DELM2
                       ELSE
                          UG = MGRENZ(IIM)
                          OG = MGRENZ(IIM+1)
                       ENDIF
                    ENDIF
                    QWver1         = (0.5D0*N0KX0*(OG**2-UG**2)+STEIG*(OG**3-UG**3)/3.D0)
                    NWver1         = N0KX0*(OG-UG)+0.5D0*STEIG*(OG**2-UG**2)
                    QWver(IIM,II)  = QWver(IIM,II)+QWver1
                    NWver(IIM,II)  = NWver(IIM,II)+NWver1
                    QWFver(IIM,II) = QWFver(IIM,II)+QWver1/QWneu*QWFneu
                    if(SMAX.eq.1) then
                       QSver_fact  = NWver1/NWneu*WEI_N + QWver1/QWneu*WEI_Q
                    else
                       QSver_fact  = NWver1/NWneu
                    endif
                    QSver1         = QSneu*QSver_fact
                    QAver1         = QAneu*QSver_fact
                    QSver(IIM,II)  = QSver(IIM,II)+QSver1
                    QAver(IIM,II)  = QAver(IIM,II)+QAver1
                    ! Test of QSver/QAver
!!!            if(QSver1/QAver1.lt.0.D0.or.
!!!     &         QSver1/QAver1.gt.0.9D0)
!!!     &      print *,'cond.f LDM1: J,IM,I,QS/Aver1'
!!!     &      ,JJ,IIM,II,QSver1,QAver1,NWver1
!!!            if(QWver1.le.0.D0.or.
!!!     &         NWver1.le.0.0D0)
!!!     &      print *,'cond_shift.f LDM1: J,IM,I,Q/NWver1'
!!!     &      ,JJ,IIM,II,QWver1,NWver1
                 ENDDO
              ENDIF
           ENDIF
           ! N1 < 0
           IF(N1.LT.0.D0) THEN
              MSTAR     = 3.D0*MQUERNEU-2.D0*DELM2
              STEIGSTAR = 2.D0*NWneu/(DELM2-MSTAR)**2
              IMIN  = 1
              IMAX  = JMAX
              DO IIM=1,JMAX
                 IF(DELM2.LE.MGRENZ(IIM+1).AND.DELM2.GT.MGRENZ(IIM)) IMAX=IIM
                 IF(MSTAR.LE.MGRENZ(IIM+1).AND.MSTAR.GT.MGRENZ(IIM)) IMIN=IIM
              ENDDO
              IF(MSTAR.LT.MGRENZ(1)) IMIN=1
              IF(DELM2.LT.MGRENZ(1)) IMAX=1
              IF(MSTAR.GT.MGRENZ(JMAX+1)) IMIN=JMAX
              IF(DELM2.GT.MGRENZ(JMAX+1)) IMAX=JMAX
              IF(IMIN.EQ.IMAX) THEN
                 QWver(IMIN,II)  = QWver(IMIN,II)+QWneu
                 QWFver(IMIN,II) = QWFver(IMIN,II)+QWFneu
                 NWver(IMIN,II)  = NWver(IMIN,II)+NWneu
                 QSver(IMIN,II)  = QSver(IMIN,II)+QSneu
                 QAver(IMIN,II)  = QAver(IMIN,II)+QAneu
                 ! Test of QSver/QAver
!!!            if(QSver(IMIN,II)/QAver(IMIN,II).lt.0.D0.or.
!!!     &         QSver(IMIN,II)/QAver(IMIN,II).gt.0.9D0)
!!!     &      print *,'cond.f LDM2a: J,I,QS/Aver(J,I)'
!!!     &      ,IMIN,II,QSver(IMIN,II),QAver(IMIN,II),NWneu
                 BIL=1.D2
              ELSE
                 DO IIM=IMIN,IMAX
                    IF(IIM.EQ.IMIN) THEN
                       UG = MSTAR
                       OG = MGRENZ(IIM+1)
                    ELSE
                       IF(IIM.EQ.IMAX) THEN
                          UG = MGRENZ(IIM)
                          OG = DELM2
                       ELSE
                          UG = MGRENZ(IIM)
                          OG = MGRENZ(IIM+1)
                       ENDIF
                    ENDIF
                    QWver1        = STEIGSTAR*((OG**3-UG**3)/3.D0 -0.5D0*MSTAR*(OG**2-UG**2))
                    NWver1        = STEIGSTAR*((OG**2-UG**2)/2.D0-MSTAR*(OG-UG))
                    QWver(IIM,II) = QWver(IIM,II)+QWver1
                    NWver(IIM,II) = NWver(IIM,II)+NWver1
                    QWFver(IIM,II)= QWFver(IIM,II)+QWver1/QWneu*QWFneu
                    if(SMAX.eq.1) then
                       QSver_fact = NWver1/NWneu*WEI_N+QWver1/QWneu*WEI_Q
                    else 
                       QSver_fact = NWver1/NWneu
                    endif
                    QSver1        = QSneu*QSver_fact
                    QAver1        = QAneu*QSver_fact
                    QSver(IIM,II) = QSver(IIM,II)+QSver1
                    QAver(IIM,II) = QAver(IIM,II)+QAver1
                    ! Test of QSver/QAver
!!!            if(QSver1/QAver1.lt.0.D0.or.
!!!     &         QSver1/QAver1.gt.0.9D0)
!!!     &      print *,'cond.f LDM2: J,IM,I,QS/Aver1'
!!!     &      ,JJ,IIM,II,QSver1,QAver1,NWver1
                    
!!!            if(QWver1.le.0.D0.or.
!!!     &         NWver1.le.0.0D0)
!!!     &      print *,'cond_shift.f LDM1: J,IM,I,Q/NWver1'
!!!     &      ,JJ,IIM,II,QWver1,NWver1
                 ENDDO
              ENDIF
           ENDIF
           ! N2 < 0
           IF(N2.LT.0.D0) THEN
              MSTAR     =  3.D0*MQUERNEU-2.D0*DELM1
              STEIGSTAR = -2.D0*NWneu/(DELM1-MSTAR)**2
              IMIN      = 1
              IMAX      = JMAX
              DO IIM=1,JMAX
                 IF(DELM1.LE.MGRENZ(IIM+1).AND.DELM1.GT.MGRENZ(IIM)) IMIN=IIM
                 IF(MSTAR.LE.MGRENZ(IIM+1).AND.MSTAR.GT.MGRENZ(IIM)) IMAX=IIM
              ENDDO
              IF(DELM1.LT.MGRENZ(1)) IMIN=1
              IF(MSTAR.LT.MGRENZ(1)) IMAX=1
              IF(DELM1.GT.MGRENZ(JMAX+1)) IMIN=JMAX
              IF(MSTAR.GT.MGRENZ(JMAX+1)) IMAX=JMAX
              IF(IMIN.EQ.IMAX) THEN
                 QWver(IMIN,II)  = QWver(IMIN,II)+QWneu
                 QWFver(IMIN,II) = QWFver(IMIN,II)+QWFneu
                 NWver(IMIN,II)  = NWver(IMIN,II)+NWneu
                 QSver(IMIN,II)  = QSver(IMIN,II)+QSneu
                 QAver(IMIN,II)  = QAver(IMIN,II)+QAneu
                 ! Test of QSver/QAver
!!!            if(QSver(IMIN,II)/QAver(IMIN,II).lt.0.D0.or.
!!!     &         QSver(IMIN,II)/QAver(IMIN,II).gt.0.9D0)
!!!     &      print *,'cond.f LDM3a: J,I,QS/Aver(J,I)'
!!!     &      ,IMIN,II,QSver(IMIN,II),QAver(IMIN,II),NWneu
                 BIL = 1.D2
              ELSE
                 DO IIM=IMIN,IMAX
                    IF(IIM.EQ.IMIN) THEN
                       UG = DELM1
                       OG = MGRENZ(IIM+1)
                    ELSE
                       IF(IIM.EQ.IMAX) THEN
                          UG = MGRENZ(IIM)
                          OG = MSTAR
                       ELSE
                          UG = MGRENZ(IIM)
                          OG = MGRENZ(IIM+1)
                       ENDIF
                    ENDIF
                    QWver1         = STEIGSTAR*((OG**3-UG**3)/3.D0-0.5D0*MSTAR*(OG**2-UG**2))
                    NWver1         = STEIGSTAR*((OG**2-UG**2)/2.D0-MSTAR*(OG-UG))
                    QWver(IIM,II)  = QWver(IIM,II)+QWver1
                    NWver(IIM,II)  = NWver(IIM,II)+NWver1
                    QWFver(IIM,II) = QWFver(IIM,II)+QWver1/QWneu*QWFneu
                    if(SMAX.eq.1) then
                       QSver_fact  = NWver1/NWneu*WEI_N+QWver1/QWneu*WEI_Q
                    else
                       QSver_fact  = NWver1/NWneu
                    endif
                    QSver1        = QSneu*QSver_fact
                    QAver1        = QAneu*QSver_fact
                    QSver(IIM,II) = QSver(IIM,II)+QSver1
                    QAver(IIM,II) = QAver(IIM,II)+QAver1
                    ! Test of QSver/QAver
!!!            if(QSver1/QAver1.lt.0.D0.or.
!!!     &         QSver1/QAver1.gt.0.9D0)
!!!     &      print *,'cond.f LDM3: J,IM,I,QS/Aver1'
!!!     &      ,JJ,IIM,II,QSver1,QAver1,NWver1
                    
!!!            if(QWver1.le.0.D0.or.
!!!     &         NWver1.le.0.0D0)
!!!     &      print *,'cond_shift.f LDM1: J,IM,I,Q/NWver1'
!!!     &      ,JJ,IIM,II,QWver1,NWver1
                 ENDDO
              ENDIF
              ! End DELM1<0
           ENDIF
           
           !          write(*,*) "QS-Bilanz",JJ,II,BIL
        ELSE
           ! if no shift is calculated
           QWver(JJ,II)  = QWver(JJ,II)+QW(JJ,II,IP)+CONDQ(JJ,II)*DELTAT
           QWFver(JJ,II) = QWFver(JJ,II)+QWF(JJ,II,IP)+CONDQF(JJ,II)*DELTAT
           NWver(JJ,II)  = NWver(JJ,II)+NW(JJ,II,IP)+CONDN(JJ,II)*DELTAT
           QSver(JJ,II)  = QSver(JJ,II)+QS(JJ,II,IP)+CONDS(JJ,II)*DELTAT
           QAver(JJ,II)  = QAver(JJ,II)+QA(JJ,II,IP)+CONDA(JJ,II)*DELTAT
        ENDIF
     ENDDO
  ENDDO
  
  ! calculation of growth rates out of shifted intervals
  DO II=1,SMAX
     DO JJ=1,JMAX
        CONDQ(JJ,II)  = (QWver(JJ,II)-QW(JJ,II,IP))/DELTAT
        CONDQF(JJ,II) = (QWFver(JJ,II)-QWF(JJ,II,IP))/DELTAT
        CONDN(JJ,II)  = (NWver(JJ,II)-NW(JJ,II,IP))/DELTAT
        CONDS(JJ,II)  = (QSver(JJ,II)-QS(JJ,II,IP))/DELTAT
        CONDA(JJ,II)  = (QAver(JJ,II)-QA(JJ,II,IP))/DELTAT
        ! Test of CONDS/CONDA
        if(CONDA(JJ,II).eq.0.D0) goto 79
        goto 79
        if(CONDS(JJ,II)/CONDA(JJ,II).lt.0.D0.or.                  &
             CONDS(JJ,II)/CONDA(JJ,II).gt.0.9D0)                     &
             print *,'cond_shift.f Ende: J,I,CONDS/A(J,I)'           &
             ,JJ,II,CONDS(JJ,II),CONDA(JJ,II),CONDN(JJ,II)           &
             ,QSver(JJ,II),QS(JJ,II,IP),QAver(JJ,II),QA(JJ,II,IP)    &
             ,QSver(JJ,II)/QAver(JJ,II),QS(JJ,II,IP)/QA(JJ,II,IP)
79      continue
        goto 80
        if(NWver(JJ,II).gt.0.D0) then
           MQalt = QWver(JJ,II)/NWver(JJ,II)
           SQalt = QAver(JJ,II)/NWver(JJ,II)
           if(MQalt.lt.MGRENZ(JJ).or.MQalt.gt.MGRENZ(JJ+1)) then
              write(*,*) "cond_shift.f: MQ ",JJ,II,MQalt/MGRENZ(JJ)
              write(*,*) "NWver ",NWver(JJ,II)," QWver ",QWver(JJ,II)
              write(*,*) "QAver ",QAver(JJ,II)," QSver ",QSver(JJ,II)
           endif
           if((SQalt.lt.SGRENZ(II).or.SQalt.gt.SGRENZ(II+1))        &
        .and.SMAX.gt.1)                                        &
        write(*,*) "cond_shift.f: SQ ",JJ,II,SQalt/SGRENZ(II)  &
                                      ,NWver(JJ,II)
        endif
        if(CONDN(JJ,II).ne.0.D0) then
           MQalt = CONDQ(JJ,II)/CONDN(JJ,II)
           SQalt = CONDA(JJ,II)/CONDN(JJ,II)
           if(MQalt.lt.MGRENZ(JJ).or.MQalt.gt.MGRENZ(JJ+1)) then
              write(*,*) "cond_shift.f: COQ",JJ,II,MQalt/MGRENZ(JJ)
              write(*,*) "CONDN ",CONDN(JJ,II)," CONDQ ",CONDQ(JJ,II)
              write(*,*) "CONDA ",CONDA(JJ,II)," CONDS ",CONDS(JJ,II)
           endif
           if((SQalt.lt.SGRENZ(II).or.SQalt.gt.SGRENZ(II+1))        &
                .and.SMAX.gt.1)                                        &
                write(*,*) "cond_shift.f: COS",JJ,II,SQalt/SGRENZ(II)  &
                ,CONDN(JJ,II)
        endif
80      continue
     ENDDO
  ENDDO
300 CONTINUE
  
  
  RETURN
END SUBROUTINE COND_SHIFTxd
