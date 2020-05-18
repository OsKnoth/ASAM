! Calculation of the collisions and the changes in QW, NW, QS, and QA,
! NFROD, QFROD, QWFROD, QSFROD, QAFROD using the two-dimensional Linear
! Discrete Method LDM
! called by cloud.f
! calling koll_s.f
! definition: I1,J1 are the liquid drops
!             I2,J2 are the ice particles potentially with liquid water shell
!
! WARNING!!! THIS WORKS ONLY FOR IDENTICAL GRID RESOLUTIONS FOR DROPS
! AND ICE PARTICLES

SUBROUTINE koll_ice_dropsxd(KOLLN,KOLLQ,KOLLS,KOLLA,NW,QW,QS,QA,RHOTOT,   &
                            KOLLNFROD,KOLLQFROD,KOLLAFROD,KOLLSFROD,      &
                            KOLLQWF,NFROD,QFROD,QSFROD,QAFROD,QWFROD,     &
                            TABS,MQUER,SQUER,MIQUER,SIQUER,IP,kk,         &
                            i_lm,j_lm,k_lm,t_lm)

  INCLUDE 'HEADER3'

USE data_parallel,      ONLY: my_cart_id, nboundlines, isubpos
  
  IMPLICIT NONE

  INTEGER :: IP,kk, i_lm, j_lm, k_lm, t_lm
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: KOLLQ(JMAX,SMAX),KOLLN(JMAX,SMAX)
  DOUBLE PRECISION :: KOLLS(JMAX,SMAX),KOLLA(JMAX,SMAX)
  DOUBLE PRECISION :: MQUER(JMAX,SMAX),SQUER(JMAX,SMAX)
  DOUBLE PRECISION :: KOLLNFROD(JMAX,SMAX),KOLLQFROD(JMAX,SMAX)
  DOUBLE PRECISION :: KOLLSFROD(JMAX,SMAX),KOLLAFROD(JMAX,SMAX)
  DOUBLE PRECISION :: NFROD(JMAX,SMAX,IPMAX),QFROD(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QSFROD(JMAX,SMAX,IPMAX),QAFROD(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: MIQUER(JMAX,SMAX),SIQUER(JMAX,SMAX)
  DOUBLE PRECISION :: QWFROD(JMAX,SMAX,IPMAX),KOLLQWF(JMAX,SMAX)
  DOUBLE PRECISION :: RHOTOT,TABS
  
  INTEGER :: I,I1,I2,J,J1,J2,S,K,ikolk
  DOUBLE PRECISION :: DDQW12,DDQW21,DDNW12,DDNW21,DDQS12,DDQS21,DDQA12
  DOUBLE PRECISION :: DDQA21,DDQWF21,NHILF1,NHILF2,MNEU,SNEU,ANEU
  DOUBLE PRECISION :: XSUM1,XSUM2,NNEU,XGRE,XMIT,XDIF,F1N,F2N,F1M,F2M
  DOUBLE PRECISION :: F1S,F2S,SMALL,klein,KOLK_eff
  
  ! small constant  
! VG
!  SMALL = 1.0D-05
  SMALL = 1.0D0
  klein = 1.0D-15 
  
  ! choice of the drop-ice kernel
  ikolk = 1
  
  ! checking the collisions
  DO I1=1,SMAX             ! drops
     DO J1=1,JMAX            ! drops
! VG
!        IF(NW(J1,I1,IP).LT.SMALL.AND.QW(j1,i1,IP).LT.klein) GOTO 200
        IF(NW(J1,I1,IP).LT.SMALL) GOTO 200
        DO I2=1,SMAX          ! ice
           DO J2=1,JMAX        ! ice
! VG
!              IF(NFROD(J2,I2,IP).LT.SMALL.AND.QWFROD(J2,I2,IP).LT.klein) GOTO 100
              IF(NFROD(J2,I2,IP).LT.SMALL) GOTO 100
              call kolk_avexd(KOLK_eff,J1,I1,J2,I2,MQUER(J1,I1),SQUER(J1,I1),   &
                              MIQUER(J2,I2),SIQUER(J2,I2),ikolk,2)
              IF(KOLK_eff.EQ.0.D0) GOTO 100
              
              ! target bins
              K = KZIEL(J1,J2)
              S = SZIEL(I1,I2)
              IF(SMAX.EQ.1) S = 1
              
              ! interaction terms
              NHILF1  = NW(J1,I1,IP)*KOLK_eff
              NHILF2  = NFROD(J2,I2,IP)*KOLK_eff
              
              ! Unit volume instead of unit mass => change of air density, kg^-1 => m^-3
              ! NW, QW, QS in m^-3 and DD?? also in m^-3, therefore, only one RHOTOT
              NHILF1  = NHILF1*(RHOTOT)
              NHILF2  = NHILF2*(RHOTOT)

              ! correction due to increasing fall velocity with decreasing density
              NHILF1 = NHILF1*SQRT(RHOTOT0/RHOTOT)
              NHILF2 = NHILF2*SQRT(RHOTOT0/RHOTOT)     

              DDNW12  = NW(J1,I1,IP)*NHILF2
              !!        DDNW21=NFROD(J2,I2,IP)*NHILF1
              DDNW21  = DDNW12
              DDQW12  = QW(J1,I1,IP)*NHILF2
              DDQW21  = QFROD(J2,I2,IP)*NHILF1
              DDQS12  = QS(J1,I1,IP)*NHILF2
              DDQS21  = QSFROD(J2,I2,IP)*NHILF1
              DDQA12  = QA(J1,I1,IP)*NHILF2
              DDQA21  = QAFROD(J2,I2,IP)*NHILF1          
              DDQWF21 = QWFROD(J2,I2,IP)*nhilf1
              
              ! collision with drops out of the same size and aerosol bin
              ! J1=J2 and I1=I2
              IF((J1.EQ.J2.OR.K.EQ.JMAX).AND.(I1.EQ.I2.OR.S.EQ.SMAX)) THEN
                 
                 ! Aenderungen in den Zielklassen (K,S)              
                 KOLLNFROD(K,S) = KOLLNFROD(K,S) + DDNW21
                 KOLLQFROD(K,S) = KOLLQFROD(K,S) + DDQW21             
                 KOLLQWF(K,S)   = KOLLQWF(K,S)   + DDQW12 + DDQWF21
                 KOLLSFROD(K,S) = KOLLSFROD(K,S) + DDQS12 + DDQS21
                 KOLLAFROD(K,S) = KOLLAFROD(K,S) + DDQA12 + DDQA21
                 
                 ! Aenderungen in den Ausgangsklassen
                 KOLLN(J1,I1) = KOLLN(J1,I1) - DDNW12
                 KOLLQ(J1,I1) = KOLLQ(J1,I1) - DDQW12
                 KOLLS(J1,I1) = KOLLS(J1,I1) - DDQS12
                 KOLLA(J1,I1) = KOLLA(J1,I1) - DDQA12
#ifdef _VG_DEBUG_
         ! VG Some debugging, it's in every if
                 IF  (KOLLN(J1,I1).le. -1.D0 ) THEN 
                    print  *,'kolln 1',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,KOLLN(J1,I1),&
                    DDNW12,NW(J1,I1,IP),NHILF2,NFROD(J2,I2,IP),KOLK_eff,K,S
                 ENDIF
                 IF (KOLLQ(J1,I1).le.-1.D-10) THEN 
                    print *,'kollq 1',t_lm,i_lm, j_lm, k_lm,  J1,J2,I1,I2,KOLLQ(J1,I1),&
                    DDQW12,QW(J1,I1,IP),NHILF2,NFROD(J2,I2,IP),KOLK_eff,K,S
                 ENDIF        
#endif

                 KOLLNFROD(J2,I2) = KOLLNFROD(J2,I2) - DDNW21
                 KOLLQFROD(J2,I2) = KOLLQFROD(J2,I2) - DDQW21
                 KOLLQWF(J2,I2)   = KOLLQWF(J2,I2)   - DDQWF21

#ifdef _VG_DEBUG_
                 IF ((KOLLNFROD(J2,I2).le.-1.D-1).OR.KOLLQFROD(J2,I2).le.-1.D-20.OR.KOLLQWF(J2,I2).le.-1.D-20) THEN 
                    print *,'koll*frod 1',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,K,S
                    print *,'           ',NFROD(J2,I2,IP) ,QFROD(J2,I2,IP) ,QWFROD(J2,I2,IP), NW(J1,I1,IP),    QW (j1,i1,ip)
                    print *,'           ',KOLLNFROD(J2,I2),KOLLQFROD(J2,I2),KOLLQWF(J2,I2),   KOLLN(j1,i1), KOLLQ(j1,i1)
                    print *,'           ',DDNW21,          DDQW21,          DDQWF21,          DDNW12,          DDQW12
                    print *,'           ',nneu,            mneu,  f2n,  f2m
                    print *,'           ',qwfrod(j2,i2,ip)/nfrod(j2,i2,ip)/mmitte(j2),mmitte(j2),j2
                    print *,'           ',NHILF1,NHILF2,KOLK_eff,RHOtot
                 ENDIF
#endif

                 KOLLSFROD(J2,I2) = KOLLSFROD(J2,I2) - DDQS21
                 KOLLAFROD(J2,I2) = KOLLAFROD(J2,I2) - DDQA21              
              END IF
              
              ! collision of a drop and an ice particle with the same size J1=J2, 
              ! but different aerosol content I1.ne.I2!
              IF((J1.EQ.J2.OR.K.EQ.JMAX).AND.(I1.NE.I2.AND.S.LT.SMAX)) THEN
                 
                 !          MNEU = DDQW12 + DDQWF21
                 ! parameters for linear approximation (aerosol)
                 XSUM1 = SGRENZ(I1)   + SGRENZ(I2)
                 XSUM2 = SGRENZ(I1+1) + SGRENZ(I2+1)
                 SNEU  = DDQS12       + DDQS21
                 ANEU  = DDQA12       + DDQA21
                 NNEU  = DDNW12
                 XDIF  = SDIFF21(I1)  + SDIFF21(I2)
                 XMIT  = SMITTE(I1)   + SMITTE(I2)
                 XGRE  = SGRENZ(S+1)

                 CALL koll_sxd(XSUM1,XSUM2,XGRE,ANEU,NNEU,XDIF,XMIT,F2S,F2N,20)
                 F2M   = F2N
                 
                 ! Aenderungen in den Zielklassen (K,S)
                 KOLLNFROD(K,S) = KOLLNFROD(K,S) + F2N*DDNW21
                 KOLLQFROD(K,S) = KOLLQFROD(K,S) + F2M*DDQW21             
                 KOLLQWF(K,S)   = KOLLQWF(K,S)   + F2M*(DDQW12 + DDQWF21)
                 KOLLSFROD(K,S) = KOLLSFROD(K,S) + F2S*SNEU       
                 KOLLAFROD(K,S) = KOLLAFROD(K,S) + F2S*ANEU          

                 
                 KOLLNFROD(K,S+1) = KOLLNFROD(K,S+1) + (1.D0 - F2N)*DDNW21
                 KOLLQFROD(K,S+1) = KOLLQFROD(K,S+1) + (1.D0 - F2M)*DDQW21
                 KOLLQWF(K,S+1)   = KOLLQWF(K,S+1)   + (1.D0 - F2M)*(DDQW12 + DDQWF21)
                 KOLLSFROD(K,S+1) = KOLLSFROD(K,S+1) + (1.D0 - F2S)*SNEU
                 KOLLAFROD(K,S+1) = KOLLAFROD(K,S+1) + (1.D0 - F2S)*ANEU

                 ! Aenderungen in den Ausgangsklassen
                 KOLLN(J1,I1) = KOLLN(J1,I1) - DDNW12
                 KOLLQ(J1,I1) = KOLLQ(J1,I1) - DDQW12
                 KOLLS(J1,I1) = KOLLS(J1,I1) - DDQS12
                 KOLLA(J1,I1) = KOLLA(J1,I1) - DDQA12

#ifdef _VG_DEBUG_
                 IF  (KOLLN(J1,I1).le. -1.D0 ) THEN 
                    print  *,'kolln 2',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,KOLLN(J1,I1),&
                    & DDNW12,NW(J1,I1,IP),NHILF2,NFROD(J2,I2,IP),KOLK_eff,K,S,f2n
                 ENDIF
                 IF (KOLLQ(J1,I1).le.-1.D-10) THEN 
                    print *, 'kollq 2',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,KOLLQ(J1,I1),&
                    & DDQW12,QW(J1,I1,IP),NHILF2,NFROD(J2,I2,IP),KOLK_eff,K,S,f2m
                 ENDIF
#endif                 
                 KOLLNFROD(J2,I2) = KOLLNFROD(J2,I2) - DDNW21
                 KOLLQFROD(J2,I2) = KOLLQFROD(J2,I2) - DDQW21
                 KOLLQWF(J2,I2)   = KOLLQWF(J2,I2)   - DDQWF21
#ifdef _VG_DEBUG_
                 IF ((KOLLNFROD(J2,I2).le.-1.D-1).OR.KOLLQFROD(J2,I2).le.-1.D-20.OR.KOLLQWF(J2,I2).le.-1.D-20) THEN 
                    print *,'koll*frod 2',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,K,S
                    print *,'           ',NFROD(J2,I2,IP) ,QFROD(J2,I2,IP) ,QWFROD(J2,I2,IP), NW(J1,I1,IP),    QW (j1,i1,ip)
                    print *,'           ',KOLLNFROD(J2,I2),KOLLQFROD(J2,I2),KOLLQWF(J2,I2),   KOLLN(j1,i1), KOLLQ(j1,i1)
                    print *,'           ',DDNW21,          DDQW21,          DDQWF21,          DDNW12,          DDQW12
                    print *,'           ',nneu,            mneu,  f2n,  f2m
                    print *,'           ',qwfrod(j2,i2,ip)/nfrod(j2,i2,ip)/mmitte(j2),mmitte(j2),j2
                    print *,'           ',NHILF1,NHILF2,KOLK_eff,RHOtot
                 ENDIF
#endif

                 KOLLSFROD(J2,I2) = KOLLSFROD(J2,I2) - DDQS21
                 KOLLAFROD(J2,I2) = KOLLAFROD(J2,I2) - DDQA21              
              END IF
              
              ! collision of a drop and an ice particle with same aerosol content I1=I2, 
              ! but of different size J1.ne.J2!
              IF((J1.NE.J2.AND.K.LT.JMAX).AND.(I1.EQ.I2.OR.S.EQ.SMAX)) THEN
                 
                 SNEU  = DDQS12 + DDQS21
                 ANEU  = DDQA12 + DDQA21
                 ! parameters for linear approximation (size)
                 XSUM1 = MGRENZ(J1)   + MGRENZ(J2)
                 XSUM2 = MGRENZ(J1+1) + MGRENZ(J2+1)
                 MNEU  = DDQW12       + DDQWF21
                 NNEU  = DDNW12
                 XDIF  = DIFF21(J1)   + DIFF21(J2)
                 XMIT  = MMITTE(J1)   + MMITTE(J2)
                 XGRE  = MGRENZ(K+1)
                 
                 CALL koll_sxd(XSUM1,XSUM2,XGRE,MNEU,NNEU,XDIF,XMIT,F1M,F1N,21)
                 F1S   = F1N
                 
                 ! Aenderungen in den Zielklassen (K,S)
                 KOLLNFROD(K,S) = KOLLNFROD(K,S) + F1N*DDNW21
                 KOLLQFROD(K,S) = KOLLQFROD(K,S) + F1M*DDQW21
                 KOLLQWF(K,S)   = KOLLQWF(K,S)   + F1M*(DDQW12 + DDQWF21)
                 KOLLSFROD(K,S) = KOLLSFROD(K,S) + F1S*SNEU         
                 KOLLAFROD(K,S) = KOLLAFROD(K,S) + F1S*ANEU              
                 
                 KOLLNFROD(K+1,S) = KOLLNFROD(K+1,S) + (1.D0 - F1N)*DDNW21
                 KOLLQFROD(K+1,S) = KOLLQFROD(K+1,S) + (1.D0 - F1M)*DDQW21
                 KOLLQWF(k+1,s)   = KOLLQWF(k+1,s)   + (1.D0 - F1M)*(DDQW12 + DDQWF21)
                 KOLLSFROD(K+1,S) = KOLLSFROD(K+1,S) + (1.D0 - F1S)*SNEU         
                 KOLLAFROD(K+1,S) = KOLLAFROD(K+1,S) + (1.D0 - F1S)*ANEU

                 ! Aenderungen in den Ausgangsklassen
                 KOLLN(J1,I1) = KOLLN(J1,I1) - DDNW12
                 KOLLQ(J1,I1) = KOLLQ(J1,I1) - DDQW12
                 KOLLS(J1,I1) = KOLLS(J1,I1) - DDQS12
                 KOLLA(J1,I1) = KOLLA(J1,I1) - DDQA12
#ifdef _VG_DEBUG_                 
                IF  (KOLLN(J1,I1).le. -1.D0 ) THEN 
                    print  *,'kolln 3',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,KOLLN(J1,I1),&
                    & DDNW12,NW(J1,I1,IP),NHILF2,NFROD(J2,I2,IP),KOLK_eff,K,S
                 ENDIF
                 IF (KOLLQ(J1,I1).le.-1.D-10) THEN 
                    print *,'kollq 3',t_lm, i_lm, j_lm, k_lm,  J1,J2,I1,I2,KOLLQ(J1,I1),&
                    & DDQW12,QW(J1,I1,IP),NHILF2,NFROD(J2,I2,IP),KOLK_eff,K,S
                 ENDIF
#endif
                 KOLLNFROD(J2,I2) = KOLLNFROD(J2,I2) - DDNW21
                 KOLLQFROD(J2,I2) = KOLLQFROD(J2,I2) - DDQW21
                 KOLLQWF(J2,I2)   = KOLLQWF(J2,I2)   - DDQWF21

#ifdef _VG_DEBUG_
                 IF ((KOLLNFROD(J2,I2).le.-1.D0-1).OR.KOLLQFROD(J2,I2).le.-1.D-20.OR.KOLLQWF(J2,I2).le.-1.D-20) THEN 
                    print *,'koll*frod 3',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,K,S
                    print *,'           ',NFROD(J2,I2,IP) ,QFROD(J2,I2,IP) ,QWFROD(J2,I2,IP), NW(J1,I1,IP),    QW (j1,i1,ip)
                    print *,'           ',KOLLNFROD(J2,I2),KOLLQFROD(J2,I2),KOLLQWF(J2,I2),   KOLLN(j1,i1), KOLLQ(j1,i1)
                    print *,'           ',DDNW21,          DDQW21,          DDQWF21,          DDNW12,          DDQW12
                    print *,'           ',nneu,            mneu,  f1n,  f1m
                    print *,'           ',qwfrod(j2,i2,ip)/nfrod(j2,i2,ip)/mmitte(j2),j2
                    print *,'           ',NHILF1,NHILF2,KOLK_eff,RHOtot
                    print *,'           ',mgrenz(j1),mgrenz(j2),mgrenz(j1+1),mgrenz(j2+1),&
                    & mgrenz(k+1),diff21(j1),diff21(j2),mmitte(j1),mmitte(j2)
                 ENDIF
#endif
                 KOLLSFROD(J2,I2) = KOLLSFROD(J2,I2) - DDQS21
                 KOLLAFROD(J2,I2) = KOLLAFROD(J2,I2) - DDQA21              
              END IF
              
              ! collision of two drops with different size I1>I2; new: I1.ne.I2 
              ! and different aerosol content J1>J2; new: J1.ne.J2
              IF((J1.NE.J2.AND.K.LT.JMAX).AND.(I1.NE.I2.AND.S.LT.SMAX)) THEN
                 
                 ! parameters for linear approximation (size)
                 XSUM1 = MGRENZ(J1) + MGRENZ(J2)
                 XSUM2 = MGRENZ(J1+1) + MGRENZ(J2+1)
                 MNEU  = DDQW12 + DDQWF21
                 NNEU  = DDNW12
                 XDIF  = DIFF21(J1) + DIFF21(J2)
                 XMIT  = MMITTE(J1) + MMITTE(J2)
                 XGRE  = MGRENZ(K+1)
                 
                 CALL koll_sxd(XSUM1,XSUM2,XGRE,MNEU,NNEU,XDIF,XMIT,F1M,F1N,22)
                 F1S   = F1N
                 
                 ! parameters for linear approximation (aerosol)
                 XSUM1 = SGRENZ(I1) + SGRENZ(I2)
                 XSUM2 = SGRENZ(I1+1) + SGRENZ(I2+1)
                 SNEU  = DDQS12 + DDQS21
                 ANEU  = DDQA12 + DDQA21
                 NNEU  = DDNW12
                 XDIF  = SDIFF21(I1) + SDIFF21(I2)
                 XMIT  = SMITTE(I1) + SMITTE(I2)
                 XGRE  = SGRENZ(S+1)
                 
                 CALL koll_sxd(XSUM1,XSUM2,XGRE,ANEU,NNEU,XDIF,XMIT,F2S,F2N,23)
                 F2M   = F2N
                 
                 ! Aenderungen in den Zielklassen (K,S)               
                 KOLLNFROD(K,S) = KOLLNFROD(K,S) + F1N*F2N*DDNW21
                 KOLLQFROD(K,S) = KOLLQFROD(K,S) + F1M*F2M*DDQW21             
                 KOLLQWF(K,S)   = KOLLQWF(K,S)   + F1M*F2M*(DDQW12 + DDQWF21)
                 KOLLSFROD(K,S) = KOLLSFROD(K,S) + F1S*F2S*SNEU         
                 KOLLAFROD(K,S) = KOLLAFROD(K,S) + F1S*F2S*ANEU          
                 
                 KOLLNFROD(K+1,S) = KOLLNFROD(K+1,S) + (1.D0 - F1N)*F2N*DDNW21
                 KOLLQFROD(K+1,S) = KOLLQFROD(K+1,S) + (1.D0 - F1M)*F2M*DDQW21
                 KOLLQWF(k+1,s)   = KOLLQWF(k+1,s)   + (1.D0 - F1M)*F2M*(DDQW12 + DDQWF21)
                 KOLLSFROD(K+1,S) = KOLLSFROD(K+1,S) + (1.D0 - F1S)*F2S*SNEU         
                 KOLLAFROD(K+1,S) = KOLLAFROD(K+1,S) + (1.D0 - F1S)*F2S*ANEU
                 
                 KOLLNFROD(K,S+1) = KOLLNFROD(K,S+1) + F1N*(1.D0 - F2N)*DDNW21
                 KOLLQFROD(K,S+1) = KOLLQFROD(K,S+1) + F1M*(1.D0 - F2M)*DDQW21
                 KOLLQWF(K,S+1)   = KOLLQWF(K,S+1)   + F1M*(1.D0 - F2M)*(DDQW12 + DDQWF21)
                 KOLLSFROD(K,S+1) = KOLLSFROD(K,S+1) + F1S*(1.D0 - F2S)*SNEU         
                 KOLLAFROD(K,S+1) = KOLLAFROD(K,S+1) + F1S*(1.D0 - F2S)*ANEU
                 
                 KOLLNFROD(K+1,S+1) = KOLLNFROD(K+1,S+1) + (1.D0 - F1N)*(1.D0 - F2N)*DDNW21
                 KOLLQFROD(K+1,S+1) = KOLLQFROD(K+1,S+1) + (1.D0 - F1M)*(1.D0 - F2M)*DDQW21             
                 KOLLQWF(k+1,s+1)   = KOLLQWF(k+1,s+1)   + (1.D0 - F1M)*(1.D0 - F2M)*(DDQW12 + DDQWF21)
                 KOLLSFROD(K+1,S+1) = KOLLSFROD(K+1,S+1) + (1.D0 - F1S)*(1.D0 - F2S)*SNEU         
                 KOLLAFROD(K+1,S+1) = KOLLAFROD(K+1,S+1) + (1.D0 - F1S)*(1.D0 - F2S)*ANEU 
                 
                 ! Aenderungen in den Ausgangsklassen
                 KOLLN(J1,I1) = KOLLN(J1,I1) - DDNW12
                 KOLLQ(J1,I1) = KOLLQ(J1,I1) - DDQW12
                 KOLLS(J1,I1) = KOLLS(J1,I1) - DDQS12
                 KOLLA(J1,I1) = KOLLA(J1,I1) - DDQA12

#ifdef _VG_DEBUG_                 
                 IF  (KOLLN(J1,I1).le. -1.D0 ) THEN 
                    print  *,'kolln 4',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,KOLLN(J1,I1),&
                    & DDNW12,NW(J1,I1,IP),NHILF2,NFROD(J2,I2,IP),KOLK_eff,K,S
                 ENDIF
                 IF (KOLLQ(J1,I1).le.-1.D-10) THEN 
                    print *,'kollq 4',t_lm, i_lm, j_lm, k_lm,  J1,J2,I1,I2,KOLLQ(J1,I1),& 
                    & DDQW12,QW(J1,I1,IP),NHILF2,NFROD(J2,I2,IP),KOLK_eff,K,S
                 ENDIF
#endif

                 KOLLNFROD(J2,I2) = KOLLNFROD(J2,I2) - DDNW21
                 KOLLQFROD(J2,I2) = KOLLQFROD(J2,I2) - DDQW21
                 KOLLQWF(J2,I2)   = KOLLQWF(J2,I2) - DDQWF21

#ifdef _VG_DEBUG_
                 IF ((KOLLNFROD(J2,I2).le.-1.D-1).OR.KOLLQFROD(J2,I2).le.-1.D-20.OR.KOLLQWF(J2,I2).le.-1.D-20) THEN 
                    print *,'koll*frod 4',t_lm, i_lm, j_lm, k_lm, J1,J2,I1,I2,K,S
                    print *,'           ',NFROD(J2,I2,IP) ,QFROD(J2,I2,IP) ,QWFROD(J2,I2,IP), NW(J1,I1,IP),    QW (j1,i1,ip)
                    print *,'           ',KOLLNFROD(J2,I2),KOLLQFROD(J2,I2),KOLLQWF(J2,I2),   KOLLN(j1,i1), KOLLQ(j1,i1)
                    print *,'           ',DDNW21,          DDQW21,          DDQWF21,          DDNW12,          DDQW12
                    print *,'           ',nneu,            mneu,  f2n,  f2m
                    print *,'           ',qwfrod(j2,i2,ip)/nfrod(j2,i2,ip)/mmitte(j2),mmitte(j2),j2
                    print *,'           ',NHILF1,NHILF2,KOLK_eff,RHOtot
                 ENDIF
#endif

                 KOLLSFROD(J2,I2) = KOLLSFROD(J2,I2) - DDQS21
                 KOLLAFROD(J2,I2) = KOLLAFROD(J2,I2) - DDQA21              
              END IF
              
100           CONTINUE
           END DO
        END DO
200     CONTINUE
     END DO
  END DO
  
  RETURN
END SUBROUTINE koll_ice_dropsxd

