subroutine cond_inixd(TABS,Pdyn,RHOdyn,ES,QW,QS,QA,NW,MQUER,RQUER,DMFAK,   &
                      DMFAKice,AAA,BBB,SATTeq,r_wet,r_qsa,r_sol,r_dry,     &
                      m_dry,m_sol,a_sig,S_krit,r_krit,m_krit,NFROD,        &
                      QFROD,QWFROD,RIQUER,ip,i_lm,j_lm,k_lm)
  
  ! calculation of new Kelvin and Raoult terms, DMFAK, and critical values
  ! for each size bin
  ! called by cond_mixxd.f90
  ! calling molal.f90
  !         krit.f90
  
  INCLUDE 'HEADER3'
  
  IMPLICIT NONE

  INTEGER :: II,JJ,ITER,IP,i_lm,j_lm,k_lm
  DOUBLE PRECISION :: TABS,Pdyn,RHOdyn,ES
  DOUBLE PRECISION :: MQUER(JMAX,SMAX),RQUER(JMAX,SMAX)
  DOUBLE PRECISION :: RIQUER(JMAX,SMAX),NFROD(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QFROD(JMAX,SMAX,IPMAX),QWFROD(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QW(JMAX,SMAX,IPMAX),NW(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: QS(JMAX,SMAX,IPMAX),QA(JMAX,SMAX,IPMAX)
  DOUBLE PRECISION :: DIFF,DIFFSTAR,KTHSTAR
  DOUBLE PRECISION :: DFAC,KFAC
  DOUBLE PRECISION :: FD,FK,DMFAK(JMAX,SMAX),DMFAKice(JMAX,SMAX)
  DOUBLE PRECISION :: tau,a_sig
  DOUBLE PRECISION :: molality,phi_S,AAA(JMAX,SMAX),BBB(JMAX,SMAX)
  DOUBLE PRECISION :: SATTeq(JMAX,SMAX)
  DOUBLE PRECISION :: r_wet(JMAX,SMAX),r_dry(JMAX,SMAX)
  DOUBLE PRECISION :: r_qsa(JMAX,SMAX),r_sol(JMAX,SMAX)
  DOUBLE PRECISION :: m_qsa(JMAX,SMAX),m_dry(JMAX,SMAX),m_sol(JMAX,SMAX)
  DOUBLE PRECISION :: r_krit(JMAX,SMAX),m_krit(JMAX,SMAX),S_krit(JMAX,SMAX)
  DOUBLE PRECISION :: f_vent,vel_part
  
  ! write(*,*) "cond_ini.f"
!!! IP=1
  ! temperature dependence of surface tension
  tau   = 1.D0-TABS/TC
  a_sig = b1*tau**mu*(1.D0+b2*tau)
  ! Diffusion and heat conduction
  ! old
  ! DIFF = 8.79D-5*(TABS**1.81D0)/Pdyn
  ! KTH   = 2.43D-2
  ! new
  DIFF  =  4.0122D-5 * TABS**1.94D0 / Pdyn
  KTH   = (5.69D0+0.017D0*(TABS-273.16D0)) * 418.5D-5
  DFAC  =  DIFF/ALPHA_C * SQRT(2.D0*PI/(RW*TABS))
  KFAC  =  KTH /(ALPHA_T*CP*RHOdyn) * SQRT(2.D0*PI/(RDRY*TABS))
  ! calculation of the size-dependent parameters needed once per time-step
  ! must be done seperately for ice and liquid!
  ! liquid
  DO II=1,SMAX
     DO JJ=1,JMAX
        !    IF(MQUER(JJ,II).GT.0.D0) THEN
        IF(NW(JJ,II,IP).GT.SMALL1) THEN
           DIFFSTAR = DFAC / RQUER(JJ,II)
           DIFFSTAR = DIFFSTAR + RQUER(JJ,II)/(RQUER(JJ,II)+DEL_V)
           DIFFSTAR = DIFF/DIFFSTAR
           !      DIFFSTAR=DIFF                   ! old version
           KTHSTAR  = KFAC / RQUER(JJ,II)
           KTHSTAR  = KTHSTAR + RQUER(JJ,II)/(RQUER(JJ,II)+DEL_T)
           KTHSTAR  = KTH/KTHSTAR
           !      KTHSTAR=KTH                   ! old version
           ! prefactor for growth equation
           ! liquid
           FD           = RW*TABS/(DIFFSTAR*ES)
           FK           = (LV/(RW*TABS)-1.D0)*LV/(KTHSTAR*TABS)
           DMFAK(JJ,II) = 4.D0*PI/(FK+FD)
           ! new Kelvin and Raoult terms are included
           m_sol(JJ,II) = QS(JJ,II,IP)/NW(JJ,II,IP)
           
           m_dry(JJ,II) = (QA(JJ,II,IP)-QS(JJ,II,IP))/NW(JJ,II,IP)
           m_qsa(JJ,II) = QA(JJ,II,IP)/NW(JJ,II,IP)
           r_dry(JJ,II) = (m_dry(JJ,II)/fact_AI)**QU1D3
           r_sol(JJ,II) = (m_sol(JJ,II)/fact_AS)**QU1D3
           r_qsa(JJ,II) = (r_dry(JJ,II)**3+r_sol(JJ,II)**3)**QU1D3
           r_wet(JJ,II) = (r_dry(JJ,II)**3+RQUER(JJ,II)**3)**QU1D3
           r_wet(JJ,II) = MAX(r_qsa(JJ,II),r_wet(JJ,II))
           molality     = m_sol(JJ,II)/(mol_AS*MQUER(JJ,II))
           
           call molalxd(molality,Phi_S)
           ! Kelvin
           if(icond.eq.3) then
              AAA(JJ,II) = 2.D0*a_sig/(RHOW*RW*TABS*RQUER(JJ,II))
           else
              AAA(JJ,II) = 2.D0*(a_sig+bb*molality)/(RHOW*RW*TABS*r_wet(JJ,II))
           endif
           
           ! Raoult
           if(icond.eq.3) then
              BBB(JJ,II) = molality*vantS*mol_w
           else
              BBB(JJ,II) = molality*vantS*mol_w*Phi_S
           endif
           ! equilibrium saturation pressure at surface
           if(icond.eq.1) then
              SATTeq(JJ,II) = exp(AAA(JJ,II)-BBB(JJ,II)) - 1.D0
           else
              SATTeq(JJ,II) = AAA(JJ,II)-BBB(JJ,II)
           endif
           !!          if(SATTeq(JJ,II).gt.10.D0) 
           !!     &      write(*,*) II,JJ,"cond_ini.f: SATTeq",SATTeq(JJ,II),
           !!     &                 RQUER(JJ,II)
           ! terminal velocity and ventilation coefficient
           call term_velxd(vel_part,r_wet(JJ,II),RHOdyn)
           
           call vent_coeffxd(f_vent,DIFFSTAR,vel_part,r_wet(JJ,II),TABS,RHOdyn)
           DMFAK(JJ,II) = DMFAK(JJ,II) * f_vent
           ! prefactor
           DMFAK(JJ,II) = DMFAK(JJ,II) * r_wet(JJ,II)
           
           ! calculation of S_krit and r_krit
           call kritxd(r_qsa(JJ,II),r_sol(JJ,II),r_dry(JJ,II),m_sol(JJ,II),   &
                       S_krit(JJ,II),r_krit(JJ,II),m_krit(JJ,II),a_sig,TABS,  &
                       JJ,II,i_lm,j_lm,k_lm)
        ENDIF
     ENDDO
  ENDDO
  
  ! ice
  DO II=1,SMAX
     DO JJ=1,JMAX
        IF(NFROD(JJ,II,IP).GT.SMALL1) THEN
           DIFFSTAR = DFAC / RIQUER(JJ,II)
           DIFFSTAR = DIFFSTAR + RIQUER(JJ,II)/(RIQUER(JJ,II)+DEL_V)
           DIFFSTAR = DIFF/DIFFSTAR
           !      DIFFSTAR=DIFF                   ! old version
           KTHSTAR  = KFAC / RIQUER(JJ,II)
           KTHSTAR  = KTHSTAR + RIQUER(JJ,II)/(RIQUER(JJ,II)+DEL_T)
           KTHSTAR  = KTH/KTHSTAR
           !      KTHSTAR=KTH                   ! old version
           ! growth like solid ice
           IF(QFROD(JJ,II,IP)/QWFROD(JJ,II,IP).GT.0.9D0) THEN
              FD  = RW*TABS/(DIFFSTAR*ES*MIN(1.D0,(TABS/273.16D0)**2.66D0))
              FK  = (LS/(RW*TABS)-1.D0)*LS/(KTHSTAR*TABS)
           ELSE
              ! growth like liquid water
              FD  =  RW*TABS/(DIFFSTAR*ES)
              FK  = (LV/(RW*TABS)-1.D0)*LV/(KTHSTAR*TABS)
           ENDIF
           DMFAKice(JJ,II) = 4.D0*PI/(FK+FD)
           ! terminal velocity and ventilation coefficient
           ! should be replaced by those for ice somewhen
           ! changes could balance each other (velocity of ice < that of liquid,
           ! ventilation coefficient >?)
           ! use Stokes velocities!!!
           ! Check RIQUER!!!
           call term_velxd(vel_part,RIQUER(JJ,II),RHOdyn)
           
           call vent_coeffxd(f_vent,DIFFSTAR,vel_part,RIQUER(JJ,II),TABS,RHOdyn)
           DMFAKice(JJ,II) = DMFAKice(JJ,II) * f_vent
           DMFAKice(JJ,II) = DMFAKice(JJ,II) * RIQUER(JJ,II)
           ! size contribution of particulate mass is neglected
        ENDIF
     ENDDO
  ENDDO
  
  RETURN
END subroutine cond_inixd

