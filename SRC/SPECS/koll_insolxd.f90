SUBROUTINE koll_insolxd(KOLLN_INS,KOLLQ_INS,KOLLS_INS,KOLLA_INS,NW,QW,   &
                        QS,QA,RHOTOT,KOLLNFROD_INS,KOLLQFROD_INS,        &
                        KOLLAFROD_INS,KOLLSFROD_INS,NWINS,QAINS,         &
                        KOLLNINS,KOLLAINS,TABS,MQUER,SQUER,ifreeze,ip)

! Calculation of the collisions and the
! changes in QW, NW, QS, and QA,
! NFROD, QFROD, QSFROD, QAFROD
! using the two-dimensional Linear Discrete Method LDM
!
! called by cloud.f
!
! calling koll_s.f

INCLUDE 'HEADER3'

implicit none

integer :: I,I1,I2,J,J1,J2,S,K,IP,IT,icontrol,ikolk,ifreeze(ITMAX)
DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
DOUBLE PRECISION :: KOLLQ_INS(JMAX,SMAX),KOLLN_INS(JMAX,SMAX)
DOUBLE PRECISION :: KOLLS_INS(JMAX,SMAX),KOLLA_INS(JMAX,SMAX)
DOUBLE PRECISION :: MQUER(JMAX,SMAX),SQUER(JMAX,SMAX)
DOUBLE PRECISION :: DDQW12,DDQW21,DDNW12,DDNW21
DOUBLE PRECISION :: DDQS12,DDQS21,DDQA12,DDQA21
DOUBLE PRECISION :: NWINS(SIMAX,ITMAX,IPMAX),QAINS(SIMAX,ITMAX,IPMAX)
DOUBLE PRECISION :: KOLLNINS(SIMAX,ITMAX),KOLLAINS(SIMAX,ITMAX)
DOUBLE PRECISION :: RHOTOT
DOUBLE PRECISION :: NHILF1,NHILF2
DOUBLE PRECISION :: MNEU,SNEU,ANEU,XSUM1,XSUM2,NNEU,XGRE,XMIT,XDIF
DOUBLE PRECISION :: F1N,F2N,F1M,F2M,F1S,F2S
DOUBLE PRECISION :: SMALL,small_facfreeze
DOUBLE PRECISION :: AVE,KOLK_eff,AQ,MQ
DOUBLE PRECISION :: KOLLNFROD_INS(JMAX,SMAX),KOLLQFROD_INS(JMAX,SMAX)
DOUBLE PRECISION :: KOLLSFROD_INS(JMAX,SMAX),KOLLAFROD_INS(JMAX,SMAX)

DOUBLE PRECISION :: T_CONTACT
DOUBLE PRECISION :: TABS
DOUBLE PRECISION :: fac_contact,fac_freeze,fac_liquid

 
if(iice.eq.0) T_CONTACT = 100.D0
if(iice.eq.1) T_CONTACT = 273.15D0
if(iice.eq.2) T_CONTACT = 263.15D0

! small constant  
! VG
!SMALL=1.D-05
SMALL_facfreeze=1.D-05
SMALL=1.D0

! control on/off  
icontrol=0

! checking the collisions
do I1=1,SMAX !drops
  do J1=1,JMAX !drops
    if(NW(J1,I1,IP).lt.SMALL) goto 200
      do IT=1,ITMAX !insol. AP, type
        fac_contact=0.D0
        if(TABS.lt.T_CONTACT) call contact_tempxd(TABS,fac_contact,ifreeze,IT)
        do I2=1,SIMAX !insol. AP, size
          if(NWINS(I2,IT,IP).lt.SMALL) goto 100
          ikolk=0
          MQ=0.D0
          AQ=QAINS(I2,IT,IP)/NWINS(I2,IT,IP)
          J2=I2/2       ! only for different resolutions for drops and AP
          call kolk_avexd(KOLK_eff,J1,I1,J2,I2,MQUER(J1,I1),SQUER(J1,I1),   &
                          MQ,AQ,ikolk,3)
          if(KOLK_eff.eq.0.D0) goto 100

! target bins
          K=J1
          S=SZIEL(I1,I2)
          if(SMAX.eq.1) S=1
! interaction terms
          NHILF1=NW(J1,I1,IP)*KOLK_eff
          NHILF2=NWINS(I2,IT,IP)*KOLK_eff
! Unit volume instead of unit mass => change of air density, kg^-1 => m^-3
! NW, QW, QS in m^-3 and DD?? also in m^-3, therefore, only one RHOTOT
          NHILF1=NHILF1*(RHOTOT)
          NHILF2=NHILF2*(RHOTOT)
! correction due to increasing fall velocity with decreasing density
          NHILF1=NHILF1*SQRT(RHOTOT0/RHOTOT)
          NHILF2=NHILF2*SQRT(RHOTOT0/RHOTOT)

          DDQW12=QW(J1,I1,IP)*NHILF2
          DDQS12=QS(J1,I1,IP)*NHILF2
          DDQA12=QA(J1,I1,IP)*NHILF2
          DDNW12=NW(J1,I1,IP)*NHILF2
          DDNW21=DDNW12
          DDQA21=QAINS(I2,IT,IP)*NHILF1

!          fac_freeze=CONTACT(J2,I2)*fac_contact
!          fac_liquid=(1-CONTACT(J2,I2))*fac_contact
          fac_freeze=fac_contact
          fac_liquid=1.D0-fac_contact

! collision with the same aerosol bin
! I1=I2
          if(I1.eq.I2.or.S.eq.SMAX) then
! source for (K,S)

!             if(TABS.GT.T_CONTACT.OR.fac_contact.le.SMALL) then
            if(fac_freeze.le.SMALL_facfreeze) then
              KOLLN_INS(K,S)=KOLLN_INS(K,S)+DDNW12
              KOLLQ_INS(K,S)=KOLLQ_INS(K,S)+DDQW12
              KOLLS_INS(K,S)=KOLLS_INS(K,S)+DDQS12
              KOLLA_INS(K,S)=KOLLA_INS(K,S)+DDQA12+DDQA21
            else
              KOLLN_INS(K,S)=KOLLN_INS(K,S)+DDNW12 * fac_liquid
              KOLLQ_INS(K,S)=KOLLQ_INS(K,S)+DDQW12 * fac_liquid
              KOLLS_INS(K,S)=KOLLS_INS(K,S)+DDQS12 * fac_liquid
              KOLLA_INS(K,S)=KOLLA_INS(K,S)+(DDQA12+DDQA21) * fac_liquid
              KOLLNFROD_INS(K,S)= KOLLNFROD_INS(K,S)+DDNW12 * fac_freeze
              KOLLQFROD_INS(K,S)= KOLLQFROD_INS(K,S)+DDQW12 * fac_freeze
              KOLLSFROD_INS(K,S)= KOLLSFROD_INS(K,S)+DDQS12 * fac_freeze
              KOLLAFROD_INS(K,S)= KOLLAFROD_INS(K,S)+(DDQA12+DDQA21)     &
                                * fac_freeze
            endif
  
! sink for (J1,I1)=(J2,I2)
            KOLLN_INS(J1,I1)=KOLLN_INS(J1,I1)-DDNW12
            KOLLQ_INS(J1,I1)=KOLLQ_INS(J1,I1)-DDQW12
            KOLLS_INS(J1,I1)=KOLLS_INS(J1,I1)-DDQS12
            KOLLA_INS(J1,I1)=KOLLA_INS(J1,I1)-DDQA12
            KOLLNINS(I2,IT)=KOLLNINS(I2,IT)-DDNW21
            KOLLAINS(I2,IT)=KOLLAINS(I2,IT)-DDQA21

          endif




! collision with different aerosol content I1.ne.I2
!          goto 222
          if(I1.ne.I2.and.S.lt.SMAX) then
            MNEU=DDQW12
! parameters for linear approximation (aerosol)
            XSUM1=SGRENZ(I1)+SGRENZ(I2)
            XSUM2=SGRENZ(I1+1)+SGRENZ(I2+1)
            SNEU=DDQS12
            ANEU=DDQA12+DDQA21
            NNEU=DDNW12
            XDIF=SDIFF21(I1)+SDIFF21(I2)
            XMIT=SMITTE(I1)+SMITTE(I2)
            XGRE=SGRENZ(S+1)

            CALL koll_sxd(XSUM1,XSUM2,XGRE,ANEU,NNEU,XDIF,XMIT,F2S,F2N,30)
            F2M=F2N

! source for (K,S) and (K,S+1)
!            if(TABS.GT.T_CONTACT.OR.fac_contact.le.SMALL) then
            if(fac_freeze.lt.SMALL) then
              KOLLN_INS(K,S)=KOLLN_INS(K,S)+F2N*NNEU
              KOLLQ_INS(K,S)=KOLLQ_INS(K,S)+F2M*MNEU
              KOLLS_INS(K,S)=KOLLS_INS(K,S)+F2S*SNEU
              KOLLA_INS(K,S)=KOLLA_INS(K,S)+F2S*ANEU
              KOLLN_INS(K,S+1)=KOLLN_INS(K,S+1)+(1.D0-F2N)*NNEU
              KOLLQ_INS(K,S+1)=KOLLQ_INS(K,S+1)+(1.D0-F2M)*MNEU
              KOLLS_INS(K,S+1)=KOLLS_INS(K,S+1)+(1.D0-F2S)*SNEU
              KOLLA_INS(K,S+1)=KOLLA_INS(K,S+1)+(1.D0-F2S)*ANEU
            else
              KOLLN_INS(K,S)=KOLLN_INS(K,S)+F2N*NNEU * fac_liquid
              KOLLQ_INS(K,S)=KOLLQ_INS(K,S)+F2M*MNEU * fac_liquid
              KOLLS_INS(K,S)=KOLLS_INS(K,S)+F2S*SNEU * fac_liquid
              KOLLA_INS(K,S)=KOLLA_INS(K,S)+F2S*ANEU * fac_liquid

              KOLLN_INS(K,S+1)=KOLLN_INS(K,S+1)+(1.D0-F2N)*NNEU * fac_liquid
              KOLLQ_INS(K,S+1)=KOLLQ_INS(K,S+1)+(1.D0-F2M)*MNEU * fac_liquid
              KOLLS_INS(K,S+1)=KOLLS_INS(K,S+1)+(1.D0-F2S)*SNEU * fac_liquid
              KOLLA_INS(K,S+1)=KOLLA_INS(K,S+1)+(1.D0-F2S)*ANEU * fac_liquid
               
              KOLLNFROD_INS(K,S)=KOLLNFROD_INS(K,S)+F2N*NNEU * fac_freeze
              KOLLQFROD_INS(K,S)=KOLLQFROD_INS(K,S)+F2M*MNEU * fac_freeze
              KOLLSFROD_INS(K,S)=KOLLSFROD_INS(K,S)+F2S*SNEU * fac_freeze
              KOLLAFROD_INS(K,S)=KOLLAFROD_INS(K,S)+F2S*ANEU * fac_freeze

              KOLLNFROD_INS(K,S+1)=KOLLNFROD_INS(K,S+1)+(1.D0-F2N)*NNEU   &
                                  * fac_freeze
              KOLLQFROD_INS(K,S+1)=KOLLQFROD_INS(K,S+1)+(1.D0-F2M)*MNEU   &
                                  * fac_freeze
              KOLLSFROD_INS(K,S+1)=KOLLSFROD_INS(K,S+1)+(1.D0-F2S)*SNEU   &
                                  * fac_freeze
              KOLLAFROD_INS(K,S+1)=KOLLAFROD_INS(K,S+1)+(1.D0-F2S)*ANEU   &
                                  * fac_freeze
            endif

! sinks for (J1,I1) and (I2)
            KOLLN_INS(J1,I1)=KOLLN_INS(J1,I1)-DDNW12
            KOLLQ_INS(J1,I1)=KOLLQ_INS(J1,I1)-DDQW12
            KOLLS_INS(J1,I1)=KOLLS_INS(J1,I1)-DDQS12
            KOLLA_INS(J1,I1)=KOLLA_INS(J1,I1)-DDQA12
            KOLLNINS(I2,IT)=KOLLNINS(I2,IT)-DDNW21
            KOLLAINS(I2,IT)=KOLLAINS(I2,IT)-DDQA21

          endif
222       continue

100       continue

        enddo
      enddo

200   continue

    enddo
  enddo


RETURN
END
