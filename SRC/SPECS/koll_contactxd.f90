! Calculation of the collisions and the changes in QW, NW, QS, and QA
! using the two-dimensional Linear Discrete Method LDM
! calling koll_s.f

SUBROUTINE koll_contactxd(KOLLN,KOLLQ,KOLLS,KOLLA,NW,QW,QS,QA,RHOTOT,   &
                          TABS,MQUER,SQUER,ip)
  
  INCLUDE 'HEADER3'
  
  IMPLICIT NONE
  
  INTEGER :: I,I1,I2,J,J1,J2,S,K,ip,ikolk
  DOUBLE PRECISION :: RHOTOT,TABS
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: KOLLQ(JMAX,SMAX),KOLLN(JMAX,SMAX)
  DOUBLE PRECISION :: KOLLS(JMAX,SMAX),KOLLA(JMAX,SMAX)
  DOUBLE PRECISION :: MQUER(JMAX,SMAX),SQUER(JMAX,SMAX)
  
  DOUBLE PRECISION :: DDQW12,DDQW21,DDNW12,DDNW21,DDQS12,DDQS21
  DOUBLE PRECISION :: DDQA12,DDQA21,NHILF1,NHILF2,MNEU,SNEU,ANEU
  DOUBLE PRECISION :: XSUM1,XSUM2,NNEU,XGRE,XMIT,XDIF,F1N,F2N,F1M,F2M
  DOUBLE PRECISION :: F1S,F2S,SMALL,KOLK_eff
  
  ! small constant  
! VG
!  SMALL=1.D-05
  SMALL=1.D0
  
  ! checking the collisions
  DO I1=1,SMAX
     DO J1=1,JMAX
        IF(NW(J1,I1,IP).LE.SMALL) GOTO 200
        !!    IF(SQUER(J1,I1).LT.SGRENZ(I1).AND.SMAX.GT.1)               &
        !!      WRITE(*,*) "koll_contact.f: SQ1 wrong",J1,I1,            &
        !!                 SQUER(J1,I1)/SGRENZ(I1),NW(J1,I1,IP)
        DO I2=1,I1
           DO J2=1,J1
              IF(NW(J2,I2,IP).LE.SMALL) GOTO 100
              ikolk = 0
              CALL kolk_avexd(KOLK_eff,J1,I1,J2,I2,MQUER(J1,I1),        &
                              SQUER(J1,I1),MQUER(J2,I2),SQUER(J2,I2),   &
                              ikolk,1)
              IF(KOLK_eff.EQ.0.D0) GOTO 100
              
              ! target bins
              K = KZIEL(J1,J2)
              S = SZIEL(I1,I2)
              IF(SMAX.EQ.1) S = 1
              IF(S.GT.SMAX) WRITE(*,*) 'koll_contactxd S = ', S, SMAX
              
              ! interaction terms
              NHILF1 = NW(J1,I1,IP)*KOLK_eff
              NHILF2 = NW(J2,I2,IP)*KOLK_eff
              ! Unit volume instead of unit mass => change of air density, kg^-1 => m^-3
              ! NW, QW, QS in m^-3 and DD?? also in m^-3, therefore, only one RHOTOT
              NHILF1 = NHILF1*(RHOTOT)
              NHILF2 = NHILF2*(RHOTOT)
              ! correction due to increasing fall velocity with decreasing density
              NHILF1 = NHILF1*SQRT(RHOTOT0/RHOTOT)
              NHILF2 = NHILF2*SQRT(RHOTOT0/RHOTOT)
              
              DDQW12 = QW(J1,I1,IP)*NHILF2
              DDQW21 = QW(J2,I2,IP)*NHILF1
              DDQS12 = QS(J1,I1,IP)*NHILF2
              DDQS21 = QS(J2,I2,IP)*NHILF1
              DDQA12 = QA(J1,I1,IP)*NHILF2
              DDQA21 = QA(J2,I2,IP)*NHILF1
              DDNW12 = NW(J1,I1,IP)*NHILF2
              DDNW21 = DDNW12

              ! collision with drops out of the same size and aerosol bin J1=J2 and I1=I2
              IF((J1.EQ.J2.OR.K.EQ.JMAX).AND.(I1.EQ.I2.OR.S.EQ.SMAX)) THEN
                 
                 ! source for (K,S)
                 KOLLN(K,S) = KOLLN(K,S) + DDNW12/2.D0
                 KOLLQ(K,S) = KOLLQ(K,S) + DDQW12
                 KOLLS(K,S) = KOLLS(K,S) + DDQS12
                 KOLLA(K,S) = KOLLA(K,S) + DDQA12
                 
#ifdef _VG_DEBUG_
                 IF ( DDQW12 / DDNW12/2.D0 / MGRENZ(K+1) .gt.1.D0) THEN 
                    print *,'koll_contact, wrong class - a',J1,J2,K,DDQW12 / DDNW12/2.D0 / MGRENZ(K+1)
                    print *,'source nw,qw                 ',k, nw(K,1,1),qw(k,1,1)
                    print *,'sink nw,qw                   ',j1,nw(j1,1,1),   qw(j1,1,1),j2, nw(j2,1,1),qw(j2,1,1)
                    print *,'updated koll*                ',kolln(k,s),kollq(k,s)
                    print *,'dd*                          ',ddnw12/2.D0,ddqw12,nhilf2,kolk_eff
                 ENDIF
#endif

                 ! sink for (J1,I1)=(J2,I2)
                 KOLLN(J1,I1) = KOLLN(J1,I1)-DDNW12
                 KOLLQ(J1,I1) = KOLLQ(J1,I1)-DDQW12
                 KOLLS(J1,I1) = KOLLS(J1,I1)-DDQS12
                 KOLLA(J1,I1) = KOLLA(J1,I1)-DDQA12
                 
                 ! two drops form one new drop.
                 ! Therefore, QW, QS, QA are conserved, NW is divided by 2
              END IF
                            
              
              ! collision of two drops with the same size J1=J2, but different aerosol
              ! content I1>I2
              IF((J1.EQ.J2.OR.K.EQ.JMAX).AND.(I1.GT.I2.AND.S.LT.SMAX)) THEN
                 
                 MNEU = DDQW12 + DDQW21
                 ! parameters for linear approximation (aerosol)
                 XSUM1 = SGRENZ(I1) + SGRENZ(I2)
                 XSUM2 = SGRENZ(I1+1) + SGRENZ(I2+1)
                 SNEU = DDQS12 + DDQS21
                 ANEU = DDQA12 + DDQA21
                 NNEU = DDNW12
                 XDIF = SDIFF21(I1) + SDIFF21(I2)
                 XMIT = SMITTE(I1) + SMITTE(I2)
                 XGRE = SGRENZ(S+1)
                 
                 CALL koll_sxd(XSUM1,XSUM2,XGRE,ANEU,NNEU,XDIF,XMIT,F2S,F2N,10)
                 
                 F2M = F2N
                 
                 ! source for (K,S) and (K,S+1)
                 KOLLN(K,S) = KOLLN(K,S) + F2N*NNEU
                 KOLLQ(K,S) = KOLLQ(K,S) + F2M*MNEU
                 KOLLS(K,S) = KOLLS(K,S) + F2S*SNEU
                 KOLLA(K,S) = KOLLA(K,S) + F2S*ANEU
                 
                 KOLLN(K,S+1) = KOLLN(K,S+1) + (1.D0 - F2N)*NNEU
                 KOLLQ(K,S+1) = KOLLQ(K,S+1) + (1.D0 - F2M)*MNEU
                 KOLLS(K,S+1) = KOLLS(K,S+1) + (1.D0 - F2S)*SNEU
                 KOLLA(K,S+1) = KOLLA(K,S+1) + (1.D0 - F2S)*ANEU
                 
                 ! sinks for (J1,I1) and (J2,I2)
                 KOLLN(J1,I1) = KOLLN(J1,I1) - DDNW12
                 KOLLQ(J1,I1) = KOLLQ(J1,I1) - DDQW12
                 KOLLS(J1,I1) = KOLLS(J1,I1) - DDQS12
                 KOLLA(J1,I1) = KOLLA(J1,I1) - DDQA12
                 
                 KOLLN(J2,I2) = KOLLN(J2,I2) - DDNW21
                 KOLLQ(J2,I2) = KOLLQ(J2,I2) - DDQW21
                 KOLLS(J2,I2) = KOLLS(J2,I2) - DDQS21
                 KOLLA(J2,I2) = KOLLA(J2,I2) - DDQA21
                 
              END IF
              
              ! collision of two drops with same aerosol content I1=I2, but of
              ! different size J1>J2
              IF((J1.GT.J2.AND.K.LT.JMAX).AND.(I1.EQ.I2.OR.S.EQ.SMAX)) THEN
                 
                 SNEU = DDQS12 + DDQS21
                 ANEU = DDQA12 + DDQA21
                 ! parameters for linear approximation (size)
                 XSUM1 = MGRENZ(J1) + MGRENZ(J2)
                 XSUM2 = MGRENZ(J1+1) + MGRENZ(J2+1)
                 MNEU = DDQW12 + DDQW21
                 NNEU = DDNW12
                 XDIF = DIFF21(J1) + DIFF21(J2)
                 XMIT = MMITTE(J1) + MMITTE(J2)
                 XGRE = MGRENZ(K+1)
                 
                 CALL koll_sxd(XSUM1,XSUM2,XGRE,MNEU,NNEU,XDIF,XMIT,F1M,F1N,11)
                 
                 F1S = F1N
                 
                 ! source for (K,S) and (K+1,S)
                 KOLLN(K,S) = KOLLN(K,S) + F1N*NNEU
                 KOLLQ(K,S) = KOLLQ(K,S) + F1M*MNEU
                 KOLLS(K,S) = KOLLS(K,S) + F1S*SNEU
                 KOLLA(K,S) = KOLLA(K,S) + F1S*ANEU

#ifdef _VG_DEBUG_               
                 IF ((F1M*MNEU)/(F1N*NNEU)/mgrenz(k+1) .gt. 1.D0) THEN 
                    print *,'koll_contact, wrong class - b1',J1,J2,K,(F1M*MNEU)/(F1N*NNEU)/mgrenz(k+1)
                    print *,'source nw,qw                  ',k, nw(k,1,1), qw(k,1,1)
                    print *,'sink nw,qw                    ',j1,nw(j1,1,1),   qw(j1,1,1),j2, nw(j2,1,1),qw(j2,1,1)
                    print *,'updated koll*                 ',kolln(k,s),kollq(k,s)
                    print *,'neu - values n                ',nneu,f1n
                    print *,'neu - values m                ',mneu,f1m
                    print *,'dd*12*                        ',ddnw12,ddqw12,nhilf2,kolk_eff
                    print *,'dd*12*                        ',ddnw21,ddqw21,nhilf1,kolk_eff
                 ENDIF
#endif

                 KOLLN(K+1,S) = KOLLN(K+1,S) + (1.D0 - F1N)*NNEU
                 KOLLQ(K+1,S) = KOLLQ(K+1,S) + (1.D0 - F1M)*MNEU
                 KOLLS(K+1,S) = KOLLS(K+1,S) + (1.D0 - F1S)*SNEU
                 KOLLA(K+1,S) = KOLLA(K+1,S) + (1.D0 - F1S)*ANEU

#ifdef _VG_DEBUG_               
                 IF (((1.D0-F1M)*MNEU)/((1.D0-F1N)*NNEU)/mgrenz(k+2) .gt. 1.D0) THEN 
                    print *,'koll_contact, wrong class - b2',J1,J2,K,'+1',((1-F1M)*MNEU)/((1.D0-F1N)*NNEU)/mgrenz(k+2)
                    print *,'source nw,qw                  ',k, nw(k+1,1,1), qw(k+1,1,1)
                    print *,'sink nw,qw                    ',j1,nw(j1,1,1),      qw(j1,1,1),j2, nw(j2,1,1),qw(j2,1,1)
                    print *,'updated koll*                 ',kolln(k+1,s),kollq(k+1,s)
                    print *,'neu - values n                ',nneu,(1.D0-f1n),nneu*(1.D0-f1n)
                    print *,'neu - values m                ',mneu,(1.D0-f1m)
                    print *,'dd*12*                        ',ddnw12,ddqw12,nhilf2,kolk_eff
                    print *,'dd*12*                        ',ddnw21,ddqw21,nhilf1,kolk_eff
                 ENDIF
#endif

                 
                 ! sinks for(J1,I1) and (J2,I2)
                 KOLLN(J1,I1) = KOLLN(J1,I1) - DDNW12
                 KOLLQ(J1,I1) = KOLLQ(J1,I1) - DDQW12
                 KOLLS(J1,I1) = KOLLS(J1,I1) - DDQS12
                 KOLLA(J1,I1) = KOLLA(J1,I1) - DDQA12
                 
                 KOLLN(J2,I2) = KOLLN(J2,I2) - DDNW21
                 KOLLQ(J2,I2) = KOLLQ(J2,I2) - DDQW21
                 KOLLS(J2,I2) = KOLLS(J2,I2) - DDQS21
                 KOLLA(J2,I2) = KOLLA(J2,I2) - DDQA21
                 
              END IF
              
              
              ! collision of two drops with different size I1>I2 and different aerosol
              ! content J1>J2
              IF((J1.GT.J2.AND.K.LT.JMAX).AND.(I1.GT.I2.AND.S.LT.SMAX)) THEN
                 
                 ! parameters for linear approximation (size)
                 XSUM1 = MGRENZ(J1) + MGRENZ(J2)
                 XSUM2 = MGRENZ(J1+1) + MGRENZ(J2+1)
                 MNEU = DDQW12 + DDQW21
                 NNEU = DDNW12
                 XDIF = DIFF21(J1) + DIFF21(J2)
                 XMIT = MMITTE(J1) + MMITTE(J2)
                 XGRE = MGRENZ(K+1)
                 
                 CALL koll_sxd(XSUM1,XSUM2,XGRE,MNEU,NNEU,XDIF,XMIT,F1M,F1N,12)
                 
                 F1S = F1N
                 
! parameters for linear approximation (aerosol)
                 XSUM1 = SGRENZ(I1) + SGRENZ(I2)
                 XSUM2 = SGRENZ(I1+1) + SGRENZ(I2+1)
                 SNEU = DDQS12 + DDQS21
                 ANEU = DDQA12 + DDQA21
                 NNEU = DDNW12
                 XDIF = SDIFF21(I1) + SDIFF21(I2)
                 XMIT = SMITTE(I1) + SMITTE(I2)
                 XGRE = SGRENZ(S+1)
                 
                 CALL koll_sxd(XSUM1,XSUM2,XGRE,ANEU,NNEU,XDIF,XMIT,F2S,F2N,13)
                 
                 F2M = F2N
                 
                 ! sources for (K,S), (K+1,S), (K,S+1) and (K+1,S+1)
                 KOLLN(K,S) = KOLLN(K,S) + F1N*F2N*NNEU
                 KOLLQ(K,S) = KOLLQ(K,S) + F1M*F2M*MNEU
                 KOLLS(K,S) = KOLLS(K,S) + F1S*F2S*SNEU
                 KOLLA(K,S) = KOLLA(K,S) + F1S*F2S*ANEU
                 
                 KOLLN(K+1,S) = KOLLN(K+1,S) +(1.D0 - F1N)*F2N*NNEU
                 KOLLQ(K+1,S) = KOLLQ(K+1,S) +(1.D0 - F1M)*F2M*MNEU
                 KOLLS(K+1,S) = KOLLS(K+1,S) +(1.D0 - F1S)*F2S*SNEU
                 KOLLA(K+1,S) = KOLLA(K+1,S) +(1.D0 - F1S)*F2S*ANEU
                 
                 KOLLN(K,S+1) = KOLLN(K,S+1) + F1N*(1.D0 - F2N)*NNEU
                 KOLLQ(K,S+1) = KOLLQ(K,S+1) + F1M*(1.D0 - F2M)*MNEU
                 KOLLS(K,S+1) = KOLLS(K,S+1) + F1S*(1.D0 - F2S)*SNEU
                 KOLLA(K,S+1) = KOLLA(K,S+1) + F1S*(1.D0 - F2S)*ANEU
                 
                 KOLLN(K+1,S+1) = KOLLN(K+1,S+1) + (1.D0 - F1N)*(1.D0 - F2N)*NNEU
                 KOLLQ(K+1,S+1) = KOLLQ(K+1,S+1) + (1.D0 - F1M)*(1.D0 - F2M)*MNEU
                 KOLLS(K+1,S+1) = KOLLS(K+1,S+1) + (1.D0 - F1S)*(1.D0 - F2S)*SNEU
                 KOLLA(K+1,S+1) = KOLLA(K+1,S+1) + (1.D0 - F1S)*(1.D0 - F2S)*ANEU
                 
                 ! sinks for (J1,I1) and (J2,I2)
                 KOLLN(J1,I1) = KOLLN(J1,I1) - DDNW12
                 KOLLQ(J1,I1) = KOLLQ(J1,I1) - DDQW12
                 KOLLS(J1,I1) = KOLLS(J1,I1) - DDQS12
                 KOLLA(J1,I1) = KOLLA(J1,I1) - DDQA12
                 
                 KOLLN(J2,I2) = KOLLN(J2,I2) - DDNW21
                 KOLLQ(J2,I2) = KOLLQ(J2,I2) - DDQW21
                 KOLLS(J2,I2) = KOLLS(J2,I2) - DDQS21
                 KOLLA(J2,I2) = KOLLA(J2,I2) - DDQA21
                 
              END IF
              
100           CONTINUE
           END DO
        END DO
        
200     CONTINUE
     END DO
  END DO
  
  
  RETURN
END SUBROUTINE koll_contactxd

