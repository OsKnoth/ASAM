! with water vapor deposition on ice!!!

SUBROUTINE COND_MIXXD(CONDN,CONDQ,CONDS,CONDA,r_qsa,r_sol,r_dry,MQUER,         &
                      RQUER,SQUER,NW,QW,QS,QA,HOEHE,ITIM,TIM,QV,DQV,           &
                      DQVcond,DQVdep,RHOdyn,ES,Pdyn,SATT,TABS,                 &
                      CONDNFROD,CONDQFROD,CONDSFROD,CONDAFROD,CONDQWFROD,      &
                      NFROD,QFROD,QSFROD,QAFROD,qwfrod,MIQUER,RIQUER,SIQUER    &
                      ,ip,i_lm,j_lm,k_lm)

  ! calculation of the condensation/evaporation and deposition
  ! iterative calculation of the saturation ratio
  !
  ! SATT:    saturation ratio due to the dynamics
  ! SATTini: estimated value for the iteration
  ! SATTneu: value after the calculation of the microphysics
  !   => SATTneu = SATTini ???
  !     - yes: everything okay, no more iteration
  !     - no:  new SATTini, one more iteration
  !
  ! called by cloud.f
  !
  ! calling cond_ini.f
  !         krit.f
  !         cond_shift.f
  
  
  INCLUDE 'HEADER2'   
  
  USE utilities_mp_ice,       ONLY:    &
       saturation
  
  IMPLICIT NONE
  
  INTEGER :: II,JJ,ITIM,IP,IPP,i_lm,j_lm,k_lm
  DOUBLE PRECISION :: SATT,TABS,Pdyn,RHOdyn
  DOUBLE PRECISION :: QV,DQV,DQVcond,DQVdep,HOEHE,TIM
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: NFROD(JMAX,SMAX,IPMAX),QFROD(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QSFROD(JMAX,SMAX,IPMAX),QAFROD(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: qwfrod(jmax,smax,ipmax)
  DOUBLE PRECISION :: MQUER(JMAX,SMAX),MQneu,MQalt,RQalt
  DOUBLE PRECISION :: SQUER(JMAX,SMAX),SQneu,SQalt
  DOUBLE PRECISION :: RQUER(JMAX,SMAX)
  DOUBLE PRECISION :: CONDQ(JMAX,SMAX),CONDN(JMAX,SMAX)
  DOUBLE PRECISION :: CONDS(JMAX,SMAX),CONDA(JMAX,SMAX)
  DOUBLE PRECISION :: CONDNFROD(JMAX,SMAX),CONDQFROD(JMAX,SMAX)
  DOUBLE PRECISION :: CONDQWFROD(JMAX,SMAX)
  DOUBLE PRECISION :: CONDSFROD(JMAX,SMAX),CONDAFROD(JMAX,SMAX)
  DOUBLE PRECISION :: MIQUER(JMAX,SMAX),SIQUER(JMAX,SMAX)
  DOUBLE PRECISION :: RIQUER(JMAX,SMAX)
  DOUBLE PRECISION :: S_krit(JMAX,SMAX),r_krit(JMAX,SMAX),m_krit(JMAX,SMAX)
  DOUBLE PRECISION :: a_krit,b_krit
  DOUBLE PRECISION :: m_ggw,r_ggw
  DOUBLE PRECISION :: r_dry(JMAX,SMAX),m_dry(JMAX,SMAX),RHOS
  DOUBLE PRECISION :: r_wet(JMAX,SMAX),m_sol(JMAX,SMAX)
  DOUBLE PRECISION :: r_qsa(JMAX,SMAX),r_sol(JMAX,SMAX)
  DOUBLE PRECISION :: DIFF,ES,FD,FK,DMFAK(JMAX,SMAX),DMW
  DOUBLE PRECISION :: DMFAKice(JMAX,SMAX)
  DOUBLE PRECISION :: AAFAK
  DOUBLE PRECISION :: bhilf,AAA(JMAX,SMAX),BBB(JMAX,SMAX),SATTeq(JMAX,SMAX),a_sig
  DOUBLE PRECISION :: TABSneu,QVneu,SATTneu,SATTini,SATTice,EE
  DOUBLE PRECISION :: SATTtop,SATTbot
  DOUBLE PRECISION :: SATTmax,SATTmin
  DOUBLE PRECISION :: MQ
  DOUBLE PRECISION :: DUMMY(JMAX,SMAX,IPMAX),DDUMMY(JMAX,SMAX)
  DOUBLE PRECISION :: esattxd
  
  INTEGER :: ITER,IBED,irueck,iakt
  DOUBLE PRECISION :: SATTini90(ITERMAX),SATTneu90(ITERMAX)
  DOUBLE PRECISION :: SATTmax90(ITERMAX),SATTmin90(ITERMAX)
  
  
  iakt=2
  ! switches
  IBED=1           ! type of exit condition for iteration
  irueck=1         ! feedback on thermodynamics yes/no

  ! condensation/deposition
  ! with Kelvin, Raoult; without ventilation, kinetics
  ! calculation of new Kelvin and Raoult terms, DMFAK, and critical values
  ! for each size bin
  ! before the iteration and mass loops using a subroutine
  
  call cond_inixd(TABS,Pdyn,RHOdyn,ES,QW,QS,QA,NW,MQUER,RQUER,DMFAK,      &
                  DMFAKice,AAA,BBB,SATTeq,r_wet,r_qsa,r_sol,r_dry,        &
                  m_dry,m_sol,a_sig,S_krit,r_krit,m_krit,NFROD,           &
                  QFROD,QWFROD,RIQUER,ip, &
                  i_lm,j_lm,k_lm)
  
  ! without feedback
  if(irueck.le.0) then
     SATTini=SATT
     TABSneu=TABS
     goto 444
  endif

  ! Iteration to find the correct saturation ratio
  ! initial value: SATTini
  !      write(99,*) "cond.f: SATT",SATT,SATTini
  
  TABSneu=TABS
  SATTini=SATT
  SATTmax=5.D-2
  SATTmin=-7.D-1
  ITER=0
  ! Begin of the iteration loop for SATTini
222 CONTINUE
  ITER=ITER+1
  IF(ITER.GT.ITERMAX) write(*,*) "too many iterations",ITER
  ! IF(ITER.GT.ITERMAX) STOP
  IF(ITER.GT.ITERMAX) GOTO 333
444 continue
  DQVcond=0.D0
  DQVdep=0.D0
  ! condensation for all bins
  ! liquid phase
  
  ! VG modify??! corrected missing IP loop
  
  
  DO II=1,SMAX
     DO JJ=1,JMAX
        IF(NW(JJ,II,IP).GT.SMALL1) THEN
           ! mquer > m_krit: activated drop => growth rate out of DGE
           IF(MQUER(JJ,II).GT.m_krit(JJ,II)) THEN
              DMW=DMFAK(JJ,II) * (SATTini-SATTeq(JJ,II))
              ! if the growth rate is highly negative, then the equilibrium size is used
              ! (possible only for SATTini < S_krit)
              if((MQUER(JJ,II)+DMW*DELTAT).LT.m_krit(JJ,II)) THEN
                 ! call of the equilibrium routine
                 call ggw_newxd(m_ggw,r_ggw,SATTini,TABSneu,r_qsa(JJ,II),       &
                                r_sol(JJ,II),r_dry(JJ,II),m_sol(JJ,II),a_sig,   &
                                r_krit(JJ,II),S_krit(JJ,II))
                 ! evaporation: particle has to be .ge. m_ggw
                 IF((MQUER(JJ,II)+DMW*DELTAT).LT.m_ggw) DMW=(m_ggw-MQUER(JJ,II))/DELTAT
              endif
              !   mquer < m_krit: non-activated   
           ELSE
              ! SATTini < S_krit: particles remain non-activated
              IF(SATTini.LT.S_krit(JJ,II)) THEN
                 ! call of the equilibrium routine
                 call ggw_newxd(m_ggw,r_ggw,SATTini,TABSneu,r_qsa(JJ,II),       &
                                r_sol(JJ,II),r_dry(JJ,II),m_sol(JJ,II),a_sig,   &
                                r_krit(JJ,II),S_krit(JJ,II))
              ELSE
                 ! SATTini > S_krit: particles become activated
                 ! critical size is used
                 m_ggw=m_krit(JJ,II)
              ENDIF
              ! growth rate
              DMW=DMFAK(JJ,II)* (SATTini-SATTeq(JJ,II))
              ! Comparison of growth rate and equilibrium size
              IF(((MQUER(JJ,II)-m_ggw)*DMW.GT.0.D0)  .OR.(m_ggw.GT.MQUER(JJ,II) .and.   &
                 m_ggw.LT.(MQUER(JJ,II)+DMW*DELTAT)) .OR.(m_ggw.LT.MQUER(JJ,II) .and.   &
                 m_ggw.GT.(MQUER(JJ,II)+DMW*DELTAT)))                                   &
              THEN
                 DMW=(m_ggw-MQUER(JJ,II))/DELTAT
              ENDIF
           ENDIF
           CONDQ(JJ,II)=NW(JJ,II,IP)*DMW
           DQVcond=DQVcond-CONDQ(JJ,II)
        ELSE
           CONDQ(JJ,II)=0.D0
        ENDIF
        CONDN(JJ,II)=0.D0
        CONDS(JJ,II)=0.D0
        CONDA(JJ,II)=0.D0
     ENDDO
  ENDDO
  
  ! ice phase 
  ! without form factor, density, Raoult & Kelvin term, ...
  !CMS why only for negative temperatures???
  !Evaporation possible for all temp.!
  IF(TABSneu.LT.273.16) THEN
     SATTice = MIN(1.D0,(TABSneu/273.16D0)**2.66D0)
     SATTice = (SATTini + 1.D0)/SATTice - 1.D0
     DO II=1,SMAX
        DO JJ=1,JMAX
           !!      IF(NFROD(JJ,II,IP).GT.0.0D0) THEN
           !! CMS: Ice particle has to consist out of .ge. 90 % frozen water
           IF(NFROD(JJ,II,IP).GT.SMALL1) THEN
              IF(DABS(qfrod(jj,ii,ip)/qwfrod(jj,ii,ip)).GE.0.9D0) THEN
                 ! growth rate out of DGE
                 DMW = DMFAKice(JJ,II)*SATTice
                 CONDQFROD(JJ,II) = NFROD(JJ,II,IP)*DMW
                 CONDQWFROD(JJ,II) = CONDQFROD(JJ,II)   ! VG
                 !CMS growth for CONDQWFROD missing???
              ELSE
                 ! CMS 
                 ! Ice particle with liquid shell due to freezing/melting
                 ! Growth rate is similar as for supercooled droplet?
                 ! definition of CONDQWFROD(JJ,II) with
                 ! different influence on latent heat than CONDQFROD
                 DMW = DMFAKice(JJ,II)*SATTini
                 CONDQWFROD(JJ,II) = NFROD(JJ,II,IP)*DMW
                 CONDQFROD(JJ,II) = 0.D0   ! VG
              ENDIF
              ! ACHTUNG! Vermeidet bis jetzt nur, dass mehr Eis verdampft als
              ! da ist. Vorlaeufige Variante!
              ! CMS 
              ! Introduction of lower size limit of frozen particles comparable to
              ! interstitial aerosol particles?
              ! Melting due to solved/soluble material within the ice particle?
              IF((qfrod(jj,ii,ip)+condqfrod(jj,ii)*deltat).LT.0.0D0) THEN
                 condqfrod(jj,ii) = -qfrod(jj,ii,ip)/deltat
              END IF
              DQVdep = DQVdep - CONDQFROD(JJ,II)
              DQVcond= DQVcond- (CONDQWFROD(JJ,II)-CONDQFROD(JJ,II))
           ELSE
              CONDQFROD(JJ,II) = 0.D0
              CONDQWFROD(JJ,II) = 0.D0
           END IF
           CONDNFROD(JJ,II) = 0.D0
           CONDSFROD(JJ,II) = 0.D0
           CONDAFROD(JJ,II) = 0.D0
        END DO
     END DO
  END IF

  ! end condensation

  ! new RHOAdyn is not calculated
  QVneu=QV+(DQVcond+DQVdep)*DELTAT
  TABSneu=TABS-(DQVcond*LV+DQVdep*LS)/CP*DELTAT
  ES=esattxd(TABSneu)
  EE=QVneu*RHOdyn*RW*TABSneu

  CALL saturation (SATTneu,TABSneu, Pdyn, QVneu)  ! SATT scalar
  
  if(irueck.eq.0) DQVcond=0.D0
  if(irueck.eq.0) DQVdep=0.D0
  if(irueck.le.0) goto 333
  if(IBED.eq.0) goto 333
  ! boundaries for saturation ratio, new estimated value
  SATTtop=SATTini*(1.D0+SC_cond)
  SATTbot=SATTini*(1.D0-SC_cond)
  SATTmax90(ITER)=SATTmax
  SATTmin90(ITER)=SATTmin
  SATTini90(ITER)=SATTini
  SATTneu90(ITER)=SATTneu
  ! too much QV consumed 
  if(MIN(SATTtop,SATTbot).gt.SATTneu) then
     SATTmax=SATTini
     SATTmin=MAX(SATTneu,SATTmin)
     if((SATTmax-SATTmin).lt.1.D-10) then
!!!          write(*,*) "too much QV consumed",ITER
        goto 322
     endif
     SATTini=(SATTini+SATTmin)*0.5D0
     GOTO 222
  endif
  ! too less QV consumed
  if(MAX(SATTtop,SATTbot).lt.SATTneu) then
     SATTmin=SATTini
     SATTmax=MIN(SATTneu,SATTmax)
     if((SATTmax-SATTmin).lt.1.D-10) then
!!!          write(*,*) "too less QV consumed",ITER
        goto 322
     endif
     SATTini=(SATTini+SATTmax)*0.5D0
     GOTO 222
  endif
  goto 333

  ! correction when SATTini and SATTneu are too different     
322 continue
!!!      do II=1,ITER
!!!        write(90,*) II,SATTini90(II),SATTneu90(II)
!!!     &                ,SATTmax90(II),SATTmin90(II)
!!!      enddo
!!!      write(90,*)
  
  ! End calculation of the condensation/evaporation rates
  
333 continue
  ! only for control!!!
  do JJ=1,JMAX/2
     II=32
     !        write(*,*) "r_krit ",JJ,II,r_krit(JJ,II)*1.D6
     !     &                            ,r_krit(JJ,II+1)*1.D6
  enddo
  
  
!!!      write(75,111) SATTini*100.D0,HOEHE,SATTneu*100.D0
!!!     &             ,SATTice*100.D0,TIM
111 FORMAT(5F14.4)
  !      write(76,*) TIM/TIMMAX,ITER
  
  DQV=DQV+DQVcond+DQVdep
  !      write(*,*) DQV*1.D3,DQVcond*1.D3,DQVdep*1.D3
  SATT=SATTini
  
  
  ! TEST for growth rates
  do JJ=1,JMAX
     do II=1,SMAX
        if(NW(JJ,II,IP).gt.0.D0) then
           MQ=(QW(JJ,II,IP)+CONDQ(JJ,II)*DELTAT)/NW(JJ,II,IP)
           !!           if(MQ.lt.0.D0) then
           !!             write(*,*) JJ,II," cond_mix.f: MQ ",MQ/MGRENZ(JJ)
           !!             write(*,*) QW(JJ,II,IP),CONDQ(JJ,II)*DELTAT,NW(JJ,II,IP)
           !!           endif
        endif
        if(NFROD(JJ,II,IP).gt.0.D0) then
           MQ=(QFROD(JJ,II,IP)+CONDQFROD(JJ,II)*DELTAT)/NFROD(JJ,II,IP)
!!!           if(MQ.lt.0.D0) then
!!!             write(*,*) JJ,II," cond_mix.f: MQI ",MQ/MGRENZ(JJ)
!!!             write(*,*) QFROD(JJ,II,IP),CONDQFROD(JJ,II)*DELTAT
!!!     &                 ,NFROD(JJ,II,IP)
!!!           endif
        endif
     enddo
  enddo
  ! Initialize DUMMY and DDUMMY fields
  do JJ=1,JMAX
     do II=1,SMAX
        do IPP=1,IPMAX
           DUMMY(JJ,II,IPP)=0.D0
        enddo
        DDUMMY(JJ,II)=0.D0
     enddo
  enddo
  
  
  ! now Linear Discrete Method for the correct bin classification
  ! for the liquid water
!!!      IP=1
  !      write(*,*) "cond_mix.f: before cond_shift liquid"
  call COND_SHIFTxd(QW,DUMMY,NW,QS,QA,CONDQ,DDUMMY,CONDN               &
                    ,CONDS,CONDA,MQUER,SQUER,IP)
  !      write(*,*) "cond_mix.f: after cond_shift liquid"
  
  ! for the frozen drops
  !!      call COND_SHIFTxd(QFROD,NFROD,QSFROD,QAFROD
  !!     &      ,CONDQFROD,CONDNFROD,CONDSFROD,CONDAFROD
  !!     &      ,MIQUER,SIQUER,IP)
  call COND_SHIFTxd(QWFROD,QFROD,NFROD,QSFROD,QAFROD,                  &
                    CONDQWFROD,CONDQFROD,CONDNFROD,                     &
! VG 12/04/07 corrected condqfrod to condqwfrod above
                    CONDSFROD,CONDAFROD,MIQUER,SIQUER,IP)
  
  RETURN
END SUBROUTINE COND_MIXXD
