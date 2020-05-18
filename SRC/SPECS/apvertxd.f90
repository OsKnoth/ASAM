SUBROUTINE apvertxd(nw,qw,qs,qa,nwins,qains,sattap,tabsap,ip,it,imode,    &
                    iap,dn_ap,dp_ap,sig_ap,epsi,hhmax,hhz,rho,rho0)
  
  ! Calculation of the dry and wet aerosol particle distributions at the
  ! beginning and initialisation of QW, QS, QA, NW
  ! called by wmain.f
  ! calling krit.f
  !         impact.f
  !         ggw_new.f
  
  INCLUDE 'HEADER2'  
  
  IMPLICIT NONE
  
  INTEGER :: I,J,jmode,iap,IP,IT,imode,ista,stamax
  PARAMETER(stamax = 5)
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QWwet2d(JMAX,SMAX),NWwet2d(JMAX,SMAX),QSwet2d(JMAX,SMAX)
  DOUBLE PRECISION :: QAwet2d(JMAX,SMAX)
  DOUBLE PRECISION :: epsilon_wet2d(JMAX,SMAX)
  DOUBLE PRECISION :: QWwet1d(JMAX),NWwet1d(JMAX),QSwet1d(JMAX)
  DOUBLE PRECISION :: QAwet1d(JMAX)
  DOUBLE PRECISION :: epsilon_wet1d(JMAX)
  DOUBLE PRECISION :: Mswe1d(JMAX),Nswe1d(JMAX),Rswe
  DOUBLE PRECISION :: epsi,epsilon(JMAX,SMAX),m_ggw,r_ggw
  DOUBLE PRECISION :: NWINS(SIMAX,ITMAX,IPMAX),QAINS(SIMAX,ITMAX,IPMAX)
  DOUBLE PRECISION :: TABSap,SATTap
  DOUBLE PRECISION :: epsi_sta(stamax),mass_sta(stamax)
  DOUBLE PRECISION :: S_krit,m_krit,r_krit,a_sig
  ! logNormal distribution
  DOUBLE PRECISION :: NA0,RA0,SIGMA
  DOUBLE PRECISION :: LNSIG,FAC_SIG,EXP_SIG
  DOUBLE PRECISION :: NA,RA,MA,MAS,FAC_PR,R_EXP
  DOUBLE PRECISION :: a_epsi,b_epsi,epsi_l,epsi_h
  ! EFEU
  ! maximum number of modes for the AP distributions
  INTEGER :: mmax
  PARAMETER(mmax = 4)
  DOUBLE PRECISION :: dp_ap(mmax),dn_ap(mmax),sig_ap(mmax)
  DOUBLE PRECISION :: Ndry1d(JMAX),Mdry1d(JMAX)
  DOUBLE PRECISION :: RHOAIR,Nsmall,Nlarge,Nall
  DOUBLE PRECISION :: nnn,r_dry,r_sol
  
  DOUBLE PRECISION :: hhmax,hhz,rho,rho0
  
  ! dry logN-distribution for the AP with subsequent swelling of the AP
  ! to the given relative humidity
  DO J=1,JMAX
     Ndry1d(J) =0.D0
     Mdry1d(J) =0.D0
     Nswe1d(J) =0.D0
     Mswe1d(J) =0.D0
     NWwet1d(J)=0.D0
     QWwet1d(J)=0.D0
     QSwet1d(J)=0.D0
     QAwet1d(J)=0.D0
     DO I=1,SMAX
        NWwet2d(J,I)=0.D0
        QWwet2d(J,I)=0.D0
        QSwet2d(J,I)=0.D0
        QAwet2d(J,I)=0.D0
     END DO
  END DO
  
  DO jmode=1,imode
     NA0=dn_ap(jmode)
     ! VG die Berechnung produziert sinnlos grosse Anzahlen in der Hoehe
     !    --> Unterscheidung gt / le faellt raus
     !    IF(hhz.LE.hhmax.AND.hhmax.NE.0.0D0) NA0 = NA0/EXP(hhz/hhmax)
     !    IF(hhz.GT.hhmax.AND.hhmax.NE.0.0D0) NA0 = NA0/EXP(hhmax/hhmax)
     !    NA0 = NA0/rho0
     IF (hhmax .ne. 0.D0)  NA0 = NA0/EXP(hhz/hhmax)
     NA0    = NA0/rho
     RA0    = dp_ap(jmode)
     SIGMA  = sig_ap(jmode)
     RHOAIR = rho
!!$  WRITE(*,*) "parameters number (mg^-1) ",NA0*1.D-6
!!$  WRITE(*,*) "           radius (nm)    ",RA0*1.D9
!!$  WRITE(*,*) "           sigma          ",SIGMA
!!$  WRITE(*,*)
     LNSIG   = LOG(SIGMA)
     FAC_SIG = 1.D0/(SQRT(2.D0*PI)*LNSIG)
     EXP_SIG = .5D0/LNSIG**2
!!!  RA=RAmin
     RA      = RA0/SIGMA**3.0D0
     RA      = MAX(RA,RAmin)
     ! VG changed to former initialization MODIFY
     RA      = 1.D-9
     FAC_PR  = 2.D0*(PR-1.D0)/(PR+1.D0)
400  CONTINUE
     R_EXP   = LOG(RA/RA0)*LOG(RA/RA0)
     NA      = NA0*FAC_SIG * FAC_PR * EXP(-R_EXP*EXP_SIG)
     
     r_dry   = RA*(1.D0-epsi)**QU1D3
     r_sol   = RA*epsi**QU1D3
     MAS     = r_sol**3.D0 * fact_AS
     MA      = MAS+r_dry**3.D0 * fact_AI
     
     IF(epsi.EQ.0.D0) THEN
        m_ggw = 0.D0
        r_ggw = 0.D0
     ELSE
        CALL kritxd(RA,r_sol,r_dry,MAS,S_krit,r_krit,m_krit,a_sig,TABSap,J,I, &
                    ! VG debug 
                    0,0,0)
        CALL ggw_newxd(m_ggw,r_ggw,SATTap,TABSap,RA,r_sol,r_dry,MAS,a_sig,    &
                       r_krit,S_krit)
     END IF
     ! classification according to water/aerosol mass
     
     IF(epsi.EQ.0.D0) THEN
        DO I=1,SIMAX
           IF(MA.GE.SGRENZ(I).AND.MA.LT.SGRENZ(I+1)) THEN
              NWINS(I,IT,IP) = NWINS(I,IT,IP)+NA
              QAINS(I,IT,IP) = QAINS(I,IT,IP)+NA*MA
              GOTO 500
           END IF
        END DO
     ELSE
        DO J=1,JMAX
           if(m_ggw.ge.MGRENZ(J).and.m_ggw.lt.MGRENZ(J+1)) then
              NWwet1d(J) = NWwet1d(J)+NA
              QWwet1d(J) = QWwet1d(J)+NA*m_ggw
              QSwet1d(J) = QSwet1d(J)+NA*MAS
              QAwet1d(J) = QAwet1d(J)+NA*MA
              epsilon_wet1d(J) = epsilon_wet1d(J) + NA*MAS
              do I=1,SMAX
                 if((MA.ge.SGRENZ(I).and.MA.lt.SGRENZ(I+1)).or.SMAX.eq.1) then
                    NWwet2d(J,I) = NWwet2d(J,I)+NA
                    QWwet2d(J,I) = QWwet2d(J,I)+NA*m_ggw
                    QSwet2d(J,I) = QSwet2d(J,I)+NA*MAS
                    QAwet2d(J,I) = QAwet2d(J,I)+NA*MA
                    epsilon_wet2d(J,I)=epsilon_wet2d(J,I) + NA*MAS
                    goto 500
                 endif
              enddo
           endif
        enddo
     endif
500  continue
     ! classification according to aerosol particle size
     do J=1,JMAX
        if(RA.gt.RGRENZ(J).and.RA.le.RGRENZ(J+1)) then
           Ndry1d(J) = Ndry1d(J)+NA
           Mdry1d(J) = Mdry1d(J)+NA*MA
           goto 600
        endif
     enddo
600  continue
     ! classification according to wet aerosol particle size
     do J=1,JMAX
        Rswe = (r_ggw**3.D0 + RA**3.D0 * (1.D0-epsi))**QU1D3
        if(Rswe.gt.RGRENZ(J).and.Rswe.le.RGRENZ(J+1)) then
           Nswe1d(J) = Nswe1d(J)+NA
           Mswe1d(J) = Mswe1d(J)+NA*(MA+m_ggw)
           goto 700
        endif
     enddo
700  continue
     RA = RA*PR
!!!        if(RA.le.1.D-5) goto 400
     if(RA.LE.MIN(1.0D-5,RA0*SIGMA**3.0D0)) goto 400
777  continue
  enddo
  if(epsi.eq.0.D0) goto 887
  do I=1,SMAX
     do J=1,JMAX
        if(NWwet2d(J,I).lt.1.D-10) then
           NW(J,I,IP)   = 0.D0
           QW(J,I,IP)   = 0.D0
           QS(J,I,IP)   = 0.D0
           QA(J,I,IP)   = 0.D0
           epsilon(J,I) = 0.D0
        else
           NW(J,I,IP)   = NWwet2d(J,I)
           QW(J,I,IP)   = QWwet2d(J,I)
           QS(J,I,IP)   = QSwet2d(J,I)
           QA(J,I,IP)   = QAwet2d(J,I)
           !      epsilon(J,I)=epsilon_wet2d(J,I)/QAwet2d(J,I)
           epsilon(J,I) = QS(J,I,IP)/QA(J,I,IP)
           if(QS(J,I,IP).gt.QA(J,I,IP))                                   &
                write(*,*) "apvert.f: soluble>total AP mass",J,I,QS(J,I,IP)/QA(J,I,IP)
        endif
     enddo
  enddo
  
  
887 continue
  
  if(epsi.eq.0.D0)  return
  
  do I=1,SMAX
     do J=1,JMAX
        !    write(*,*) "ap1: NW(J,I)",J,I,NW(J,I),NWwet2d(J,I)
        !    if(NW(J,I,IP).ge.1.D-10) then
        if(NW(J,I,IP).ge.1.D-30) then
           MA    = (QA(J,I,IP)-QS(J,I,IP))/NW(J,I,IP)
           r_dry = (MA/fact_AI)**QU1D3
           MAS   = QS(J,I,IP)/NW(J,I,IP)
           r_sol = (MAS/fact_AS)**QU1D3
           call kritxd(RA,r_sol,r_dry,MAS,S_krit,r_krit,m_krit,a_sig,TABSap,J,I,0,0,0)
           call ggw_newxd(m_ggw,r_ggw,SATTap,TABSap,RA,r_sol,r_dry,MAS    &
                          ,a_sig,r_krit,S_krit)
           QW(J,I,IP) = NW(J,I,IP)*m_ggw
        else
           NW(J,I,IP) = 0.D0
           QW(J,I,IP) = 0.D0
           QA(J,I,IP) = 0.D0
           QS(J,I,IP) = 0.D0
        endif
     enddo
  enddo
    
  RETURN
END SUBROUTINE apvertxd









