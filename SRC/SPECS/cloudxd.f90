! control module for the microphysical processes
! this one includes deposition of water vapor on ice !!!
!
! calling:
!   cond_mixxd.f90
!   break.f90
!   koll.f90  not active, collision drop-drop, drop-AP, no freezing
!   koll_contact.f90  contact freezing AP-drop and drop-drop+ -AP collision
!                   works for T>T_freeze like koll.f
!                   for T<T_freeze like koll.f plus contact freezing
!   koll_ice_drop.f90  collisions ice particles-drops, contact freezing
!   koll_insol.f90 collisions drops-insoluble particles, contact freezing
!                 if critical temperature is reached

SUBROUTINE cloudxd(                                       &
     DELT,deltdyn,DP,DRHO,DNW,DQW,DQS,DQA,DQV,dqvdyn,     &
     SATT,NW,QW,QS,QA,QV,PTOT,RHOTOT,TABS,                &
     MQUER,RQUER,SQUER,MIQUER,RIQUER,SIQUER,              &
     DNFROD,DQFROD,DQSFROD,DQAFROD,dqwfrod,               &
     NFROD,QFROD,QSFROD,QAFROD,qwfrod,                    &
     DNWINS,DQAINS,NWINS,QAINS,                           &
     TABSold,ifreeze,miv,ip,kk,vvtw,vvtf,                 &
     ! VG debug
     i_lm,j_lm,k_lm,t_lm                                       )

  INCLUDE 'HEADER3'
  
  USE data_parallel, ONLY: my_cart_id
  
  USE data_mp_ice,   ONLY:    &
       ! DEBUG
       icecounter,filecounter,filecounter2,printcounter
  
  
  
  IMPLICIT NONE
  
  INTEGER :: II,JJ,IP,IT,ifreeze(ITMAX),itim,iruf,kk,jjj,i_lm,j_lm,k_lm,l_lm, idummy,t_lm,j
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QWSUM,NWSUM,QSSUM,QASUM
  DOUBLE PRECISION :: DQW(JMAX,SMAX),DNW(JMAX,SMAX)
  DOUBLE PRECISION :: DQS(JMAX,SMAX),DQA(JMAX,SMAX)
  DOUBLE PRECISION :: BREAQ(JMAX,SMAX),BREAN(JMAX,SMAX),BREAS(JMAX,SMAX)
  DOUBLE PRECISION :: KOLLQ(JMAX,SMAX),KOLLN(JMAX,SMAX),KOLLS(JMAX,SMAX)
  DOUBLE PRECISION :: CONDQ(JMAX,SMAX),CONDN(JMAX,SMAX),CONDS(JMAX,SMAX)
  DOUBLE PRECISION :: CONDA(JMAX,SMAX),BREAA(JMAX,SMAX),KOLLA(JMAX,SMAX)
  DOUBLE PRECISION :: MQUER(JMAX,SMAX),RQUER(JMAX,SMAX),SQUER(JMAX,SMAX)
  DOUBLE PRECISION :: MIQUER(JMAX,SMAX),RIQUER(JMAX,SMAX),SIQUER(JMAX,SMAX)
  DOUBLE PRECISION :: r_dry(JMAX,SMAX)
  DOUBLE PRECISION :: r_qsa(JMAX,SMAX),r_sol(JMAX,SMAX)
  DOUBLE PRECISION :: TABS,DELT,DQV,DQVcond,DQVdep
  DOUBLE PRECISION :: EE,ES,QV,PTOT,SATT,RHOTOT,esattxd
  DOUBLE PRECISION :: TABSdyn,QVdyn,Pdyn,RHOdyn,DP,DRHO
  DOUBLE PRECISION :: TABSold,deltdyn,dqvdyn,hoehe,tim
  DOUBLE PRECISION :: vvtw(jmax),vvtf(jmax)
  
  DOUBLE PRECISION :: NFROD(JMAX,SMAX,IPMAX),QFROD(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QSFROD(JMAX,SMAX,IPMAX),QAFROD(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: qwfrod(jmax,smax,ipmax)
  DOUBLE PRECISION :: DNFROD(JMAX,SMAX),DQFROD(JMAX,SMAX),NFRODSUM,QFRODSUM
  DOUBLE PRECISION :: DQSFROD(JMAX,SMAX),DQAFROD(JMAX,SMAX),dqwfrod(jmax,smax)
  DOUBLE PRECISION :: dnfmelt(jmax,smax),dqfmelt(jmax,smax),dqfwmelt(jmax,smax)
  DOUBLE PRECISION :: dqfsmelt(jmax,smax),dqfamelt(jmax,smax),dnwmelt(jmax,smax)
  DOUBLE PRECISION :: dqwmelt(jmax,smax),dqwsmelt(jmax,smax),dqwamelt(jmax,smax)
  DOUBLE PRECISION :: dqffrier(jmax,smax)
  DOUBLE PRECISION :: CONDNFROD(JMAX,SMAX),CONDQFROD(JMAX,SMAX)
  DOUBLE PRECISION :: CONDQWFROD(JMAX,SMAX)
  DOUBLE PRECISION :: CONDSFROD(JMAX,SMAX),CONDAFROD(JMAX,SMAX)
  DOUBLE PRECISION :: KOLLNFROD(JMAX,SMAX),KOLLQFROD(JMAX,SMAX) 
  DOUBLE PRECISION :: KOLLAFROD(JMAX,SMAX),KOLLSFROD(JMAX,SMAX)
  DOUBLE PRECISION :: NWINS(SIMAX,ITMAX,IPMAX),QAINS(SIMAX,ITMAX,IPMAX)
  DOUBLE PRECISION :: DNWINS(SIMAX,ITMAX),DQAINS(SIMAX,ITMAX)
  DOUBLE PRECISION :: NWINSSUM(ITMAX),QAINSSUM(ITMAX)
  DOUBLE PRECISION :: NWINSTOT,QAINSTOT
  DOUBLE PRECISION :: KOLLNINS(SIMAX,ITMAX),KOLLAINS(SIMAX,ITMAX)
  DOUBLE PRECISION :: KOLLNINS_FR(SIMAX,ITMAX),KOLLAINS_FR(SIMAX,ITMAX)
  DOUBLE PRECISION :: KOLLN_INS(JMAX,SMAX),KOLLQ_INS(JMAX,SMAX) 
  DOUBLE PRECISION :: KOLLA_INS(JMAX,SMAX),KOLLS_INS(JMAX,SMAX) 
  DOUBLE PRECISION :: KOLLNFROD_INS(JMAX,SMAX),KOLLQFROD_INS(JMAX,SMAX) 
  DOUBLE PRECISION :: KOLLAFROD_INS(JMAX,SMAX),KOLLSFROD_INS(JMAX,SMAX) 
  ! changes due to contact freezing by collision drop-frozen drop
  DOUBLE PRECISION :: KOLLQI(JMAX,SMAX),KOLLNI(JMAX,SMAX),KOLLSI(JMAX,SMAX)
  DOUBLE PRECISION :: KOLLAI(JMAX,SMAX),kollqwf(jmax,smax)
  DOUBLE PRECISION :: KOLLNFRODI(JMAX,SMAX),KOLLQFRODI(JMAX,SMAX) 
  DOUBLE PRECISION :: KOLLAFRODI(JMAX,SMAX),KOLLSFRODI(JMAX,SMAX)
  ! immersion freezing
  DOUBLE PRECISION :: IMMERQ(JMAX,SMAX),IMMERN(JMAX,SMAX)
  DOUBLE PRECISION :: IMMERA(JMAX,SMAX),IMMERS(JMAX,SMAX) 
  ! temperature change owing to phase transition liquid -> ice
  DOUBLE PRECISION :: DELT_ICE 
  ! Depositionsgefieren
  DOUBLE PRECISION :: deponi(simax,itmax),depoqia(simax,itmax)
  DOUBLE PRECISION :: deponf(jmax,smax),depoqf(jmax,smax)
  DOUBLE PRECISION :: depoqfa(jmax,smax),dqq,miv(itmax)
  ! Kollision Wasser/Eis - Wasser/Eis
  DOUBLE PRECISION :: knf(jmax,smax),kqf(jmax,smax),kqwf(jmax,smax)
  DOUBLE PRECISION :: kqsf(jmax,smax),kqaf(jmax,smax)
  

  ! Initialization
  ! Increments are set to zero
  DO II=1,SMAX
     DO JJ=1,JMAX
        ! Break-up of drops, no latent heat release
        BREAN(JJ,II) = 0.D0
        BREAQ(JJ,II) = 0.D0
        BREAS(JJ,II) = 0.D0
        BREAA(JJ,II) = 0.D0
        
        ! water vapor -> liquid: latent heat release
        CONDN(JJ,II) = 0.D0
        CONDQ(JJ,II) = 0.D0
        CONDS(JJ,II) = 0.D0
        CONDA(JJ,II) = 0.D0
        
        ! water vapor -> solid: latent heat release 
        CONDNFROD(JJ,II) = 0.D0
        CONDQFROD(JJ,II) = 0.D0
        CONDQWFROD(JJ,II)= 0.D0
        CONDSFROD(JJ,II) = 0.D0
        CONDAFROD(JJ,II) = 0.D0
        
        ! Coalescence of drops, no latent heat release
        KOLLN(JJ,II) = 0.D0
        KOLLQ(JJ,II) = 0.D0
        KOLLS(JJ,II) = 0.D0
        KOLLA(JJ,II) = 0.D0
        
        ! "Riming", no latent heat release
        KOLLNFROD(JJ,II) = 0.D0
        KOLLQFROD(JJ,II) = 0.D0
        KOLLAFROD(JJ,II) = 0.D0
        KOLLSFROD(JJ,II) = 0.D0
        
        ! Coalescence of ?????, no latent heat release
        KOLLNI(JJ,II) = 0.D0
        KOLLQI(JJ,II) = 0.D0
        KOLLSI(JJ,II) = 0.D0
        KOLLAI(JJ,II) = 0.D0
        
        ! Coalescence of ?????, no latent heat release
        KOLLNFRODI(JJ,II) = 0.D0
        KOLLQFRODI(JJ,II) = 0.D0
        KOLLAFRODI(JJ,II) = 0.D0
        KOLLSFRODI(JJ,II) = 0.D0
        kollqwf(jj,ii) = 0.0D0
        
        ! immersion freezing, liquid -> solid
        IMMERQ(JJ,II) = 0.D0
        IMMERN(JJ,II) = 0.D0
        IMMERA(JJ,II) = 0.D0 
        IMMERS(JJ,II) = 0.D0
        
        ! Deposition freezing, vapor -> solid?
        deponf(jj,ii) = 0.0D0
        depoqf(jj,ii) = 0.0D0
        depoqfa(jj,ii) = 0.0D0
        
        ! Melting, solid -> liquid
        dnfmelt(jj,ii) = 0.0D0
        dqfmelt(jj,ii) = 0.0D0
        dqfwmelt(jj,ii) = 0.0D0
        dqfsmelt(jj,ii) = 0.0D0
        dqfamelt(jj,ii) = 0.0D0
        dnwmelt(jj,ii) = 0.0D0
        dqwmelt(jj,ii) = 0.0D0
        dqwsmelt(jj,ii) = 0.0D0
        dqwamelt(jj,ii) = 0.0D0
        
        ! Freezing, liquid -> solid
        dqffrier(jj,ii) = 0.0D0
        
        ! collision with insoluble particles, latent heat release?
        KOLLN_INS(JJ,II) = 0.D0
        KOLLQ_INS(JJ,II) = 0.D0
        KOLLS_INS(JJ,II) = 0.D0
        KOLLA_INS(JJ,II) = 0.D0
        
        KOLLNFROD_INS(JJ,II) = 0.D0
        KOLLQFROD_INS(JJ,II) = 0.D0
        KOLLAFROD_INS(JJ,II) = 0.D0
        KOLLSFROD_INS(JJ,II) = 0.D0
        
        ! Kollison Wasser/Eis - Wasser/Eis, no latent heat release
        knf(jj,ii) = 0.0D0
        kqf(jj,ii) = 0.0D0
        kqwf(jj,ii) = 0.0D0
        kqsf(jj,ii) = 0.0D0
        kqaf(jj,ii) = 0.0D0
     END DO
  END DO
  
  DO IT=1,ITMAX
     DO II=1,SIMAX
        ! collision insoluble particles
        KOLLNINS(II,IT)    = 0.D0
        KOLLAINS(II,IT)    = 0.D0
        KOLLNINS_FR(II,IT) = 0.D0
        KOLLAINS_FR(II,IT) = 0.D0
        ! Deposition
        deponi(ii,it)      = 0.0D0
        depoqia(ii,it)     = 0.0D0
     END DO
  END DO
  
  ! Calculation of the saturation ratio due to the "dynamics"
  TABSdyn    = TABS + deltdyn*DELTAT
  QVdyn      = QV + dqvdyn*DELTAT
  Pdyn       = PTOT
  RHOdyn     = RHOTOT
  
  ES = esattxd(TABSdyn)
  EE = QVdyn*RHOdyn*RW*TABSdyn
  ! VG satt kommt von aussen rein!
  !SATT = EE/ES - 1.D0
  
  ! Start of the microphysics
  ! call of the condensation routine
  ! iterative calculation of the saturation ratio, Kelvin and Raoult effects
  ! non-ideal solution with partly insoluble AP
  ! with deposition of water vapor on ice
  
  IF(icond.GE.1) THEN                                             
     CALL COND_MIXxd(                                           &
          CONDN,CONDQ,CONDS,CONDA,r_qsa,r_sol,r_dry,MQUER,      &
          RQUER,SQUER,NW,QW,QS,QA,HOEHE,ITIM,TIM,QVdyn,DQV,     &
          DQVcond,DQVdep,RHOdyn,ES,Pdyn,SATT,TABSdyn,           &
          CONDNFROD,CONDQFROD,CONDSFROD,CONDAFROD,CONDQWFROD,   &
          NFROD,QFROD,QSFROD,QAFROD,qwfrod,                     &
          MIQUER,RIQUER,SIQUER,ip,                              &
          ! VG debug
          i_lm,j_lm,k_lm)
  ENDIF
  ! Depositionsgefrieren
  IF(idepo.GE.1) THEN
     dqq = 0.0D0
     CALL depoxd(TABS,SATT,dqq,NWINS,deponi,QAINS,depoqia,deponf,QFROD,    &
                 depoqf,QAFROD,depoqfa,miv,ip)
     DQV = DQV + dqq
  END IF
  
  ! call the immersion freezing routine
  IF(iimfr.GE.1) THEN 
     CALL immersion_koopxd(TABS,TABSold,NW,QW,QS,QA,IMMERQ,   &
                           IMMERN,IMMERA,IMMERS,ifreeze,ip)
  ENDIF
  
  ! Melting of ice
  IF(imelt.GE.1) THEN
     CALL schmelzenxd(TABS,PTOT,RHOTOT,vvtf,NFROD,dnfmelt,                 &
                      QFROD,dqfmelt,qwfrod,dqfwmelt,QSFROD,dqfsmelt,       &
                      QAFROD,dqfamelt,dnwmelt,dqwmelt,dqwsmelt,dqwamelt,   &
                      ip,RIQUER)
  END IF
  
  ! Gefrieren von Eisteilchen mit Wassermantel
  ! Freezing of ice particles with liquid water shell
  IF(ifrier.GE.1) CALL frierenxd(TABS,NFROD,QFROD,dqffrier,qwfrod,RIQUER,ip &
                                 ,i_lm,j_lm,k_lm,t_lm)
  
  ! call of the break-up routine, no ice included here
  ! no influence on temperature
  IF(ibrea.EQ.1) CALL breakxd(BREAN,BREAQ,BREAS,BREAA,MQUER,NW,QW,QS,QA,ip)
  
  ! call of the collision routines
  ! collision of liquid drops with liquid drops
  ! no influence on temperature
  IF(ikoll.GE.1) CALL koll_contactxd(KOLLN,KOLLQ,KOLLS,KOLLA,NW,QW,QS,    &
                                     QA,RHOTOT,TABS,MQUER,SQUER,ip)
  
  ! collision of liquid drops with frozen drops (riming/drop freezing)
  ! no influence on temperature
  IF(ikofr.GE.1) THEN
     CALL koll_ice_dropsxd(KOLLNI,KOLLQI,KOLLSI,KOLLAI,NW,QW,QS,QA,        &
                           RHOTOT,KOLLNFRODI,KOLLQFRODI,KOLLAFRODI,        &
                           KOLLSFRODI,kollqwf,NFROD,QFROD,QSFROD,QAFROD,   &
                           qwfrod,TABS,MQUER,SQUER,MIQUER,SIQUER,ip,kk     &
                                 ,i_lm,j_lm,k_lm,t_lm)
  ENDIF
  
  ! collision of insoluble particles and liquid drops, contact freezing
  ! if critical temperature is reached
  IF(iinsol.GE.1) THEN 
     CALL koll_insolxd(KOLLN_INS,KOLLQ_INS,KOLLS_INS,KOLLA_INS,NW,QW,QS,   &
                       QA,RHOTOT,KOLLNFROD_INS,KOLLQFROD_INS,              &
                       KOLLAFROD_INS,KOLLSFROD_INS,NWINS,QAINS,KOLLNINS,   &
                       KOLLAINS,TABS,MQUER,SQUER,ifreeze,ip)
  ENDIF
  
  ! Kollision von Eis/Wasser Tropfen mit Eis/Wasser Tropfen
  ! no influence on temperature
  IF(ikeis.GE.1) THEN 
     CALL koll_eis_eisxd(knf,kqf,kqwf,kqsf,kqaf,NFROD,QFROD,qwfrod,        &
                         QSFROD,QAFROD,RHOTOT,TABS,MIQUER,SIQUER,ip)
  ENDIF

  ! calculation of the temperature increment DELT owing to condensation
  DELT = DELT - (DQVcond*LV + DQVdep*LS)/CP 
  
  ! calculation of temperature change owing to phase transitions
  ! liquid -> ice and vapour -> ice
  ! CALL delt_freezexd(KOLLQFROD,KOLLQFRODI,kollqwf,KOLLQI,KOLLQFROD_INS,   &
  !                    IMMERQ,depoqf,dqfwmelt,dqwmelt,dqfmelt,dqffrier,     &
  !                    DELT_ICE)
  CALL delt_freezexd(KOLLQFROD_INS,   &
                     IMMERQ,depoqf,dqfwmelt,dqwmelt,dqfmelt,dqffrier,     &
                     DELT_ICE, i_lm, j_lm, k_lm)
  DELT = DELT + DELT_ICE
  
  ! calculation of DNW,DQW,DQS,DQA for liquid and frozen drops
  DO II=1,SMAX
     DO JJ=1,JMAX
        DNW(JJ,II) = (CONDN(JJ,II) + KOLLN(JJ,II) + KOLLNI(JJ,II)           &
                   + KOLLN_INS(JJ,II) + BREAN(JJ,II) - IMMERN(JJ,II)        &
                   + dnwmelt(jj,ii))*DELTAT
        
        DQW(JJ,II) = (CONDQ(JJ,II) + KOLLQ(JJ,II) + KOLLQI(JJ,II)           &
                   + KOLLQ_INS(JJ,II) + BREAQ(JJ,II) - IMMERQ(JJ,II)        &
                   + dqwmelt(jj,ii))*DELTAT

#ifdef _VG_DEBUG_
        IF   ( &
             (((qw(jj,ii,1)+dqw(jj,ii))/(nw(jj,ii,1)+dnw(jj,ii)))/mgrenz(jj+1).ge.1.D0) .OR.   &
             (((qw(jj,ii,1)+dqw(jj,ii))/(nw(jj,ii,1)+dnw(jj,ii)))/mgrenz(jj  ).le.1.D0)     ) THEN 
           print *,'cloud,wrong classification of qw',jj,ii,nfrod(jj,ii,1),qwfrod(jj,ii,1),i_lm,j_lm,k_lm,t_lm
           print *,'nw,qw          ', nw(jj,ii,1),           qw(jj,ii,1)
           print *,'nw,qw new      ', nw(jj,ii,1)+dnw(jj,ii),qw(jj,ii,1)+dqw(jj,ii)
           print *,'qw/mgrenz j,j+1',((qw(jj,ii,1)+dqw(jj,ii))/(nw(jj,ii,1)+dnw(jj,ii)))/mgrenz(jj  ), &
                                         ((qw(jj,ii,1)+dqw(jj,ii))/(nw(jj,ii,1)+dnw(jj,ii)))/mgrenz(jj+1)

           print *,'processes dnw: ',CONDN(JJ,II),KOLLN(JJ,II),KOLLNI(JJ,II),KOLLN_INS(JJ,II),BREAN(JJ,II)
           print *,'               ',IMMERN(JJ,II),dnwmelt(jj,ii)

           print *,'processes dqw: ',CONDQ(JJ,II),KOLLQ(JJ,II),KOLLQI(JJ,II),KOLLQ_INS(JJ,II),BREAQ(JJ,II) 
           print *,'               ',IMMERQ(JJ,II),dqwmelt(jj,ii)

        ENDIF
        ! VG
        IF (ABS(DQW (JJ,II)).GE.1.D0) THEN 
           print *,'cloud.f DQW',i_lm,j_lm,k_lm,jj,ii, &
                CONDQ(JJ,II),KOLLQ(JJ,II),KOLLQI(JJ,II),   &
                KOLLQ_INS(JJ,II),BREAQ(JJ,II),IMMERQ(JJ,II),   &
                dqwmelt(jj,ii)
        ENDIF
#endif        

        DQS(JJ,II) = (CONDS(JJ,II) + KOLLS(JJ,II) + KOLLSI(JJ,II)           &
                   + KOLLS_INS(JJ,II) + BREAS(JJ,II) - IMMERS(JJ,II)        &
                   + dqwsmelt(jj,ii))*DELTAT

        DQA(JJ,II) = (CONDA(JJ,II) + KOLLA(JJ,II) + KOLLAI(JJ,II)           &
                   + KOLLA_INS(JJ,II) + BREAA(JJ,II) - IMMERA(JJ,II)        &
                   + dqwamelt(jj,ii))*DELTAT

        DNFROD(JJ,II) = (CONDNFROD(JJ,II) + KOLLNFROD(JJ,II)                &
                      + KOLLNFRODI(JJ,II) + KOLLNFROD_INS(JJ,II)            &
                      + knf(jj,ii) + IMMERN(JJ,II) + deponf(jj,ii)          &
                      + dnfmelt(jj,ii))*DELTAT
        
        DQFROD(JJ,II) = (CONDQFROD(JJ,II) + KOLLQFROD(JJ,II)                &
                      + KOLLQFRODI(JJ,II) + KOLLQFROD_INS(JJ,II)            &
                      + kqf(jj,ii) + IMMERQ(JJ,II) + depoqf(jj,ii)          &
                      + dqfmelt(jj,ii) + dqffrier(jj,ii))*DELTAT

#ifdef _VG_DEBUG        
        IF (ABS(DQFROD (JJ,II)).GE.1.D0) THEN 
           print *,'cloud.f DQFROD',i_lm,j_lm,k_lm,jj,ii, &
                CONDQFROD(JJ,II),KOLLQFROD(JJ,II),        &
                KOLLQFRODI(JJ,II),KOLLQFROD_INS(JJ,II),   &
                kqf(jj,ii),IMMERQ(JJ,II),depoqf(jj,ii),   &
                dqfmelt(jj,ii),dqffrier(jj,ii)
        ENDIF
#endif

        DQSFROD(JJ,II) = (CONDSFROD(JJ,II) + KOLLSFROD(JJ,II)               &
                       + KOLLSFRODI(JJ,II) + KOLLSFROD_INS(JJ,II)           &
                       + kqsf(jj,ii) + IMMERS(JJ,II)                        &
                       + dqfsmelt(jj,ii))*DELTAT         
        
        DQAFROD(JJ,II) = (CONDAFROD(JJ,II) + KOLLAFROD(JJ,II)               &
                       + KOLLAFRODI(JJ,II) + KOLLAFROD_INS(JJ,II)           &
                       + kqaf(jj,ii) + IMMERA(JJ,II) + depoqfa(jj,ii)       &
                       + dqfamelt(jj,ii))*DELTAT

        dqwfrod(jj,ii) = (CONDQWFROD(JJ,II) + KOLLQFROD(JJ,II)              &
                       + kollqwf(jj,ii) + KOLLQFROD_INS(JJ,II)              &
                       + kqwf(jj,ii) + IMMERQ(JJ,II) + depoqf(jj,ii)        &
                       + dqfwmelt(jj,ii))*DELTAT

#ifdef _VG_DEBUG
        IF   ( &
             (((qwfrod(jj,ii,1)+dqwfrod(jj,ii))/(nfrod(jj,ii,1)+dnfrod(jj,ii)))/mgrenz(jj+1).ge.1.D0) .OR.   &
             (((qwfrod(jj,ii,1)+dqwfrod(jj,ii))/(nfrod(jj,ii,1)+dnfrod(jj,ii)))/mgrenz(jj  ).le.1.D0)     ) THEN 
           print *,'cloud,wrong classification of qwfrod',jj,ii,nfrod(jj,ii,1),qwfrod(jj,ii,1),i_lm,j_lm,k_lm,t_lm
           print *,'nfrod,qwfrod       ', nfrod(jj,ii,1),              qwfrod(jj,ii,1)
           print *,'nfrod,qwfrod new   ', nfrod(jj,ii,1)+dnfrod(jj,ii),qwfrod(jj,ii,1)+dqwfrod(jj,ii)
           print *,'qwfrod/mgrenz j,j+1',((qwfrod(jj,ii,1)+dqwfrod(jj,ii))/(nfrod(jj,ii,1)+dnfrod(jj,ii)))/mgrenz(jj  ), &
                                         ((qwfrod(jj,ii,1)+dqwfrod(jj,ii))/(nfrod(jj,ii,1)+dnfrod(jj,ii)))/mgrenz(jj+1)

           print *,'processes dnfrod:  ',CONDNFROD(JJ,II),KOLLNFROD(JJ,II),KOLLNFRODI(JJ,II),KOLLNFROD_INS(JJ,II)    
           print *,'                   ',knf(jj,ii),IMMERN(JJ,II),deponf(jj,ii),dnfmelt(jj,ii)

           print *,'processes dqwfrod: ',CONDQWFROD(JJ,II),KOLLQFROD(JJ,II),kollqwf(jj,ii),KOLLQFROD_INS(JJ,II)      
           print *,'                   ',kqwf(jj,ii),IMMERQ(JJ,II),depoqf(jj,ii),dqfwmelt(jj,ii)

        ENDIF

        IF (ABS(DQwFROD (JJ,II)).GE.1.D0) THEN 
           print *,'cloud.f DQWFROD',i_lm,j_lm,k_lm,jj,ii, &
                CONDQWFROD(JJ,II),KOLLQFROD(JJ,II), &
                kollqwf(jj,ii),KOLLQFROD_INS(JJ,II), &
                kqwf(jj,ii),IMMERQ(JJ,II),depoqf(jj,ii), &
                dqfwmelt(jj,ii)
        ENDIF
#endif
        
        !VG zu viel Wasser wird umgesetzt?        
        IF(ABS(DQW(JJ,II)).GT.1.D-2) THEN
           WRITE(*,*)"cloud.f 2:",i_lm, j_lm, k_lm, kk,ip,JJ,II,DQW(JJ,II),QW(JJ,II,IP)
           WRITE(*,*) CONDQ(JJ,II),KOLLQ(JJ,II),BREAQ(JJ,II)
           WRITE(*,*) kollqi(jj,ii),kollq_ins(jj,ii),immerq(jj,ii)
           !!      WRITE(*,'(4(2X,E13.6))') (qw(jjj,ii,ip), qfrod(jjj,ii,ip),        &
           !!        nw(jjj,ii,ip), nfrod(jjj,ii,ip), jjj=1,66)
        ENDIF
     END DO
  END DO
  
  DO IT=1,ITMAX
     DO II=1,SIMAX
        DNWINS(II,IT) = (KOLLNINS(II,IT) + deponi(ii,it))*DELTAT
        DQAINS(II,IT) = (KOLLAINS(II,IT) + depoqia(ii,it))*DELTAT
     END DO
  END DO
  
  ! VG Zu hohe Temperaturaenderung?
  IF(ABS(DELT).GT.1.0D0) THEN
     WRITE(*,*) "cloud.f T",i_lm, j_lm, k_lm, kk,ip,DELT,DELT-DELT_ICE,DELT_ICE,             &
                            -DQVcond*LV/CP,-DQVdep*LS/CP
     DO l_lm=1,66
        print *,'cloud.f, spectra' , &
                 NW(l_lm,1,1) , QW(l_lm,1,1), NFROD(l_lm,1,1), QFROD(l_lm,1,1), QWFROD(l_lm,1,1)
     ENDDO
  END IF




IF (my_cart_id == 2) THEN 
II = 1
IF (t_lm.ge.90 .and. t_lm.le.100) THEN 
   IF (i_lm==3 .and. j_lm==9)    THEN 

!print *,'*********************** test ***************************',my_cart_id,i_lm,j_lm,k_lm

   OPEN (999,FILE='spectra_t97.dat'  ,POSITION='APPEND')
   OPEN (888,FILE='changes_w_t97.dat',POSITION='APPEND')
   OPEN (777,FILE='changes_c1_t97.dat',POSITION='APPEND')
   OPEN (666,FILE='changes_c2_t97.dat',POSITION='APPEND')
   OPEN (555,FILE='changes_cn_t97.dat',POSITION='APPEND')
   OPEN (444,FILE='changes_wn_t97.dat',POSITION='APPEND')
   OPEN (333,FILE='meteorology_t97.dat',POSITION='APPEND')

   DO j=1,66
   write (999,*) my_cart_id,t_lm, i_lm, j_lm, k_lm, j, NW(j,1,1) , QW(j,1,1), NFROD(j,1,1), QFROD(j,1,1), QWFROD(j,1,1)
   ENDDO

   DO j=1,66
   write (888,*) t_lm, i_lm, j_lm, k_lm, j, CONDQ(J,II),  KOLLQ(J,II), KOLLQI(J,II), &
                   KOLLQ_INS(J,II), BREAQ(J,II), IMMERQ(J,II) ,dqwmelt(j,ii)
   ENDDO

   DO j=1,66
   write (777,*) t_lm, i_lm, j_lm, k_lm, j, CONDQWFROD(J,II),KOLLQFROD(J,II),kollqwf(j,ii),KOLLQFROD_INS(J,II), &
                       kqwf(j,ii),IMMERQ(J,II),depoqf(j,ii),dqfwmelt(j,ii)
   ENDDO

   DO j=1,66
   write (666,*) t_lm, i_lm, j_lm, k_lm, j, CONDQFROD(J,II),KOLLQFROD(J,II) , &
                       KOLLQFRODI(J,II) , KOLLQFROD_INS(J,II) , &
                       kqf(j,ii) , IMMERQ(J,II) , depoqf(j,ii) , &
                       dqfmelt(j,ii) , dqffrier(j,ii)
   ENDDO

   DO j=1,66
   write (555,*) t_lm, i_lm, j_lm, k_lm, j, CONDNFROD(J,II) , KOLLNFROD(J,II) , &
                       KOLLNFRODI(J,II) , KOLLNFROD_INS(J,II)            , &
                       knf(j,ii) , IMMERN(J,II) , deponf(j,ii)          , &
                       dnfmelt(j,ii)
   ENDDO


   DO j=1,66
   write (444,*) t_lm, i_lm, j_lm, k_lm, j, CONDN(J,II) , KOLLN(J,II) , KOLLNI(J,II)  , & 
                    KOLLN_INS(J,II) , BREAN(J,II) , IMMERN(J,II) , &
                    dnwmelt(j,ii)
   ENDDO

   write (333,*) t_lm, i_lm, j_lm, k_lm, DELT,deltdyn,DP,DRHO,dqvdyn,SATT,QV,PTOT,RHOTOT,TABS

   CLOSE (999) 
   CLOSE (888) 
   CLOSE (777) 
   CLOSE (666) 
   CLOSE (555) 
   CLOSE (444) 
   CLOSE (333)
   ENDIF
ENDIF
ENDIF



  RETURN
END SUBROUTINE cloudxd
