subroutine delt_freezexd(KOLLQFROD_INS,IMMERQ,depoqf,dqfwmelt,dqwmelt,     &
                         dqfmelt,dqffrier,DELT_ICE, i_lm, j_lm, k_lm)
  
  ! temperature change due to freezing of liquid water caused by collisions 
  ! drop-ap (contact freezing) or drop-ice (riming)
  
  INCLUDE 'HEADER3'
  
  IMPLICIT NONE
  
  INTEGER :: II,JJ, i_lm, j_lm, k_lm  ! VG
  DOUBLE PRECISION :: IMMERQ(jmax,smax)
  DOUBLE PRECISION :: KOLLQFROD_INS(JMAX,SMAX),depoqf(jmax,smax)
  DOUBLE PRECISION :: dqfwmelt(jmax,smax)
  DOUBLE PRECISION :: dqwmelt(jmax,smax)
  DOUBLE PRECISION :: dqfmelt(jmax,smax),dqffrier(jmax,smax)
  
  ! temperature change due to phase transitions:
  !   IMMERQ: liquid -> solid
  !   DQFMELT: solid -> liquid (larger liquid water shell)
  !   DQFWMELT: no additional phase transition compared to DQFMELT
  !   DQWMELT: solid -> liquid 
  !            (total melting of an ice particle, transfer to drops)
  !   DQFFRIER: liquid -> solid
  !   KOLLQFROD_INS: liquid -> solid (contact freezing instantaneously)
  DOUBLE PRECISION :: QFRODneu,qfneu,qfwneu,DELT_ICE  & 
       ! Debug
       ,QFROD_delta                   
  
  QFRODneu    = 0.0D0
  QFROD_delta = 0.0D0
  qfneu = 0.0D0
  qfwneu = 0.0D0
  DO II=1,SMAX
     DO JJ=1,JMAX
        QFRODneu = QFRODneu                     &
             + IMMERQ(jj,ii)                &
             !CMS             + KOLLQFROD(JJ,II)             &
             !CMS             + KOLLQFRODI(JJ,II)            &
             !CMS             + kollqwf(jj,ii)               &
             !CMS             + kollqw(jj,ii)                &
             + dqfmelt(jj,ii)               &
             + dqffrier(jj,ii)              &
             + KOLLQFROD_INS(JJ,II)
        qfneu       = qfneu + depoqf(jj,ii)
!#ifdef _VG_DEBUG_
        ! VG debug

        QFROD_delta = QFRODneu - QFROD_delta
        IF (ABS(QFROD_delta * (LV_ICE/CP)) .GE. 0.5D0) THEN 
           print *,'delt_freezexd:', i_lm, j_lm, k_lm, jj,ii,&
                QFRODneu, IMMERQ(jj,ii),dqfmelt(jj,ii), dqwmelt(jj,ii), &
                dqffrier(jj,ii),KOLLQFROD_INS(JJ,II)
        ENDIF
!#endif
        QFROD_delta = QFRODneu
        !    qfwneu = qfwneu + dqfwmelt(jj,ii) + dqwmelt(jj,ii)
     ENDDO
  ENDDO
  
  ! DELT_ICE = QFRODneu*LV_ICE/CP + qfneu*LS/CP + qfwneu*LV/CP
  DELT_ICE = QFRODneu*LV_ICE/CP + qfneu*LS/CP 
  
  RETURN
END subroutine delt_freezexd
