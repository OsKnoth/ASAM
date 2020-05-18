MODULE TwoMomentBulkMicroPhysics_Mod

  USE Control_Mod
  USE Parameter_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Chemie_Mod
  USE ParameterMicrophys_Mod

  IMPLICIT NONE

CONTAINS 

SUBROUTINE ComputeTp(T,p,RhoD,RhoV,RhoC,RhoI,PotM,c,Task,nSpecies)

  REAL(RealKind) :: T,p,RhoD,RhoV,RhoC,RhoI,PotM
  REAL(RealKind) :: c(:)
  CHARACTER*8    :: Task
  INTEGER        :: nSpecies
  REAL(RealKind) :: Rho,Kinetic,zh

  IF (Task=='Function') THEN
    T=c(nSpecies+1)+Eps
    p=c(nSpecies+2)+Eps
  ELSE
    SELECT CASE(ThetaKind)
    CASE('Density')
      p=PressureTheta(RhoD,RhoV,RhoC,RhoI,PotM)+Eps
      T=AbsTemp(RhoD,RhoV,p)+Eps
    CASE('Equiv')
      T=c(nSpecies+1)+Eps
      Rho=RhoD+RhoV+RhoC+RhoI
      CALL AbsTNewton(T,PotM/(Rho+Eps),RhoD,RhoV,RhoC)
      p=(RhoD*Rd+RhoV*Rv)*T+Eps
    CASE('Energy')
      Kinetic=c(nSpecies+3)
      zh=c(nSpecies+4)
      Rho=RhoD+RhoV+RhoC+RhoI
      T=AbsTEnergy(PotM+Eps,Rho+Eps,RhoV,RhoC,Kinetic,zh)
      p=(RhoD*Rd+RhoV*Rv)*T+Eps
    CASE('PreEn')
      p=PotM+Eps 
      T=AbsTemp(RhoD,RhoV,p)+Eps
    END SELECT
  END IF

END SUBROUTINE ComputeTp

SUBROUTINE ComputeRhsPotM(c,RhsPotM,PotM,QRhs,RhoD,RhoV,RhoC,RhoI,Rm,Cpml,Cvml,T,p,Process)

  REAL(RealKind) :: RhsPotM,PotM,QRhs,RhoD,RhoV,RhoC,RhoI,Rm,Cpml,Cvml,T,p
  REAL(RealKind) :: qD,qV,qL,qR,Rho,RhoL,RhoR
  CHARACTER*8    :: Process
  REAL(RealKind) :: c(:)
  REAL(RealKind) :: Lv,Lil,Cp_eff,eLoc,DpDRhoV,DpDRhoC

  SELECT CASE(ThetaKind)
  CASE('Density')
    IF (Process=='COND'.OR.Process=='NUCL') THEN
!     Cond --> QRhs>0
      Lv=LatHeat(T)
      RhsPotM=PotM*(                                             &
                    (-Lv/(Cpml*T+Eps)                            &
                     -LOG((p+Eps)/P0)*(Rm/Cpml)*(Rv/Rm-Cpv/Cpml) &
                     +Rv/Rm                                      &   
                    )*(-QRhs)                                    &   
                    +(LOG((p+Eps)/P0)*(Rm/Cpml)*(Cpl/Cpml)       &   
                     )*(QRhs)                                    &   
                   )   
    ELSE IF (Process=='FREEZEC'.OR.Process=='FREEZER'.OR.Process=='MELTI') THEN
      Lil=LatHeatFreeze(T)
      RhsPotM=PotM*(Lil/(Cpml*T+Eps)+LOG((p+Eps)/P0)*(Rm/Cpml)*((Cpi-Cpl)/Cpml))*QRhs
    ELSE IF (Process=='NUCLI'.OR.Process=='ICEDEP') THEN
      Lv=LatHeat(T)
      Lil=LatHeatFreeze(T)
      RhsPotM=PotM*((Lv+Lil)/(Cpml*T+Eps)-Rv/Rm+LOG((p+Eps)/P0)*(Rm/Cpml)*((Rv/Rm)+((Cpi-Cpv)/Cpml)))*QRhs
    ELSE IF (Process=='AUTOI')THEN
      Lil=LatHeatFreeze(T)
      RhsPotM=PotM*(-Lil/(Cpml*T+Eps)-2*RhoD/Rho0+LOG((p+Eps)/P0)*(Rm/Cpml)*((Rv/Rm))*((-Cpl-Cpi)/Cpml))*QRhs
    ELSE IF (Process=='AUTO'.OR.Process=='ACC') THEN
      RhsPotM=PotM*(-LOG((p+Eps)/p0)*Cpl/Cpml)*QRhs
    ELSE IF (Process=='EVAP') THEN
!     Evap --> RhsPotM < 0
      RhsPotM=PotM*(Rv/Rm-LOG((p+Eps)/p0)*(Rv/Rm-Cpv/Cpml))*QRhs
    ELSE IF (Process=='DEP') THEN

      Rho=c(1)
      RhoV=c(2)
      RhoL=c(3)
      PotM=c(4)
      RhoR=c(5)
      IF (RhoL<Zero) THEN
        RhoL=Zero
      END IF
      IF (RhoR<Zero) THEN
      RhoR=Zero
      END IF
    RhoL=RhoL+RhoR+Eps
    RhoI=Zero
    RhoD=Rho-RhoV-RhoL+Eps
    qD=RhoD/(Rho+Eps)
    qV=RhoV/(Rho+Eps)
    qL=RhoL/(Rho+Eps)
    Cpml=Cpd*RhoD+Cpv*RhoV+Cpl*RhoL+Eps
    Cvml=Cvd*RhoD+Cvv*RhoV+Cpl*RhoL+Eps
    Rm=Rd*RhoD+Rv*RhoV+Eps

    END IF
  CASE('Equiv')
    Cp_eff=Cpd+(RhoV+RhoC)/RhoD*Cpl
    RhsPotM=PotM/(Cp_eff*RhoD)*Rv*LOG(RelHumidity(T,RhoV))*QRhs
  CASE('Energy')
    RhsPotM=Zero
  CASE('PreEn')
    eLoc=Cvml*T+RhoV*L00
    DpDRhoV=(Rv-Rd)*(eLoc-RhoV*L00)/Cvml &
            -L00*(RhoD*Rd+RhoV*Rv)/Cvml &
            -(Cvv-Cvd)*(eLoc-RhoV*L00)*(RhoD*Rd+RhoV*Rv)/(Cvml**2)
    DpDRhoC=-Rd*(eLoc-RhoV*L00)/Cvml &
            -(Cpl-Cvd)*(eLoc-RhoV*L00)*(RhoD*Rd+RhoV*Rv)/(Cvml**2)
    RhsPotM=-(DpDRhoV-DpDRhoC)*QRhs
  END SELECT

END SUBROUTINE ComputeRhsPotM

SUBROUTINE MicroProcess2(c,f,Process,Task,dtAct,Time,lprint) 

  REAL(RealKind) :: c(:),f(:)
  CHARACTER*8 :: Process,Task
  REAL(RealKind) :: dtAct,Time
  LOGICAL, OPTIONAL :: lprint

  REAL(RealKind) :: Rm,Cpml,Cvml,Cp_eff,RDens
  REAL(RealKind) :: Qcond,QcondNc,Qauto,Qacc,Qdep,Qself,Qevap,Qnucl
  REAL(RealKind) :: PotM,Kinetic,zh,eLoc
  REAL(RealKind) :: RhsRhoV,RhsRhoC,RhsPotM
  REAL(RealKind) :: Forcing,Limiter

  REAL(RealKind) :: Nd,Dd,Psi,Dia,Mass,SedVel,Rey
  REAL(RealKind) :: Lambda,Nmean,Mmean,VentFactor,ThermoG
  REAL(RealKind) :: DmDtAcc,DmDtDep 
  REAL(RealKind) :: temptest 

  REAL(RealKind) :: x_s,phi,tau,D_r,kccn,sc,x_c,x_r,C_ccn,kappa_nuc,T_c,facg,fr_n,fr_q,J_hom,J_het,J_tot,T_a,x_i,D_i,S_i,v_i,G_iv
  REAL(RealKind) :: F_vi,e_v,e_coll,v_s,D_s,theta_n,theta_q,delta_n,delta_q
  REAL(RealKind) :: q_l,N_IN_diag,F_vi_0,a_ven,b_ven,F_vi_1,F_vii,D_T,N_IN_M92,N_IN_Reis
  REAL(RealKind) :: mytheta,rho,rhod,rhov,rhoc,rhor,rhoi,rhos,nv,nc,nr,ni,ns                                         ! inner cell values
  REAL(RealKind) :: mythetazr,rhozr,rhodzr,rhovzr,rhoczr,rhorzr,rhoizr,rhoszr,nvzr,nczr,nrzr,nizr,nszr,cpmlzr,rmlzr  ! upper cell values
  REAL(RealKind) :: mythetazl,rhozl,rhodzl,rhovzl,rhoczl,rhorzl,rhoizl,rhoszl,nvzl,nczl,nrzl,nizl,nszl,cpmlzl,rmlzl  ! lower cell values
  REAL(RealKind) :: drhoTh,drho,drhod,drhov,drhoc,drhor,drhoi,drhos,dnv,dnc,dnr,dni,dns,dcpml,drml                   ! inner cell rates

  REAL(RealKind) :: theta,T,Lv,p,pv,pvs,Ts,ps,pVi 
  REAL(RealKind) :: thetazr,Tzr,Lvzr,pzr,pvzr,pvszr,Tszr,pszr,PotMzr
  REAL(RealKind) :: cond,moistflux,S2,S,dSdzw,dzLoc,wLoc,wLoczr 
  REAL(RealKind) :: QcLoc

  REAL(RealKind) :: g_d,v_r,N_re,f_v,f_q,eva,eva_q,eva_n 
  REAL(RealKind) :: D_vtp
  REAL(RealKind) :: lambda_r,N_0
  INTEGER        :: nSpecies
! redefine Rho0
  Rho0=1.225d0

  SELECT CASE(Process)
  CASE('COND') ! Condensation/Evaporation cloud water (v <-> c)
    Rho =c(1)
    RhoV=c(2)
    RhoC=MAX(c(3),Zero)
    RhoR=MAX(c(4),Zero)
    NV  =c(5)
    RhoI=c(29)
    NI  =c(30)
    RhoS=c(31)
    NS  =c(32)
    PotM=c(6)
    nSpecies=6
    Nc  =c(26)
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    pVs=SaturVapor(T)
    Forcing= RhoV-pVs/(Rv*T+Eps)
    Limiter= c(3)
    Qcond  = RelCloud*(Forcing-Limiter+SQRT(Forcing*Forcing+Limiter*Limiter))
    CALL ComputeRhsPotM(c,RhsPotM,PotM,Qcond,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
    dnc =Zero
    dnc =dnc+0.01*min(Zero,RhoC/p_cloud_x_min-Nc)
!   dnc =dnc+0.01*max(Zero,RhoC/p_cloud_x_max-Nc)
!   Forcing=dnc
!   Limiter=Nv
!   dnc=RelCloud*(Forcing+Limiter-SQRT(Forcing*Forcing+Limiter*Limiter))
    f(1) =-Qcond ! QV
    f(2) = Qcond ! QC
    f(3) =-dnc   ! NV
    f(4) = dnc   ! NC 
    f(5) = RhsPotM
  CASE('NUCL') ! Nucleation (v -> c if S>0)
!   -- .chem info  ---
!   --   2 1 2 1   ---
!   --   3 6 2 4   ---
!   --   3 4 6 7   ---
!   ------------------
    Rho =c(21)
    RhoV=c(22)
    RhoC=c(23)
    RhoR=c(24)
    RhoI=c(29)
    NI  =c(30)
    RhoS=c(31)
    NS  =c(32)
    Nv  =c(25)
    Nc  =c(26)
    Nr  =c(27)
    PotM=c(28)
!   reaction variables
    RhoV=c(1)
    NV  =c(2)
    nSpecies=2 ! RhoV, NV
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    f(5)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    Rhozr =c(nSpecies+5)
    RhoVzr=c(nSpecies+6)
    RhoCzr=c(nSpecies+7)
    RhoRzr=c(nSpecies+8)
    RhoIzr=Zero
    RhoIzr=Zero
    PotMzr=c(nSpecies+9)
    wLoc  =c(nSpecies+10)
    dzLoc =c(nSpecies+11)+Eps
    RhoDzr=Rhozr-RhoVzr-RhoCzr-RhoRzr-RhoIzr-RhoSzr 
    cpmlzr   =CpmlFun(rhodzr,rhovzr,rhoczr+rhorzr,rhoizr+rhoszr)
    rmlzr    =RmFun(rhodzr,rhovzr)
    CALL ComputeTp(Tzr,pzr,RhoDzr,RhoVzr,RhoCzr+RhoRzr,RhoIzr+RhoSzr,PotMzr,c,Task,nSpecies)
    pvzr =(RhoDzr*Rd+RhoVzr*Rv)*Tzr 
    pVszr=SaturVapor(Tzr)
    pVs  =SaturVapor(T)
    S    =100.0*((RhoV  *Rv*T  /pVs  )-1.0)
    S2   =100.0*((RhoVzr*Rv*Tzr/pVszr)-1.0)
    dSdzw=(S2-S)/dzLoc*wLoc 
!   calulate nucleation rate
    Qnucl=Zero
    IF (S>Zero.AND.dSdzw>Zero.AND.T>233.15d0.AND.wLoc>Zero) THEN
      kappa_nuc = 0.462d0 
      C_ccn     = Nv
!     aerosol activation rate
      Qnucl=C_ccn*kappa_nuc*MIN(S,1.1d0)**(kappa_nuc-One)*dSdzw
    ! forcing=Qnucl 
    ! limiter=Nv 
    ! dnc    =RelCloud2*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter)) 
      dnc    =Qnucl
      drhoc  =dnc*1d-12 ! 1d-12 is smallest drop mass
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhoc,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1) =-drhoc   ! QV
      f(2) = drhoc   ! QC
      f(3) =-dnc     ! NV
      f(4) = dnc     ! NC
      f(5) = RhsPotM ! TH
    END IF
  CASE('AUTO') ! Autoconversion cloudwater to rain (c + c -> r)
    Rho =c(1)
    RhoV=MAX(c(2),Zero)
    RhoC=MAX(c(3),Zero)
    RhoR=MAX(c(4),Zero)
    RhoI=c(29)
    NI  =c(30)
    RhoS=c(31)
    NS  =c(32)
    NC  =MAX(c(5),Zero)
    PotM=c(6)
    nSpecies=6
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    f(5)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    x_s=p_cloud_x_max
    Qauto=Zero
    IF (RhoC>1.d-9) THEN
      x_c=min(max(RhoC/(NC+Eps),p_cloud_x_min),p_cloud_x_max)     ! mean mass in SI
      Qauto = p_cloud_k_au*RhoC*RhoC*x_c*x_c*Rho0/(Rho+Eps)       ! autoconversion rate (SB2000)
      IF (RhoC>1.d-6) THEN
        tau=min(max(One-RhoC/(RhoC+RhoR+Eps),Eps),One-1.0d-1) 
        phi=400.d0*tau**0.7d0*(1.-tau**0.7d0)**3.d0 
        Qauto =Qauto*(One+phi/(One-tau)**2.0d0) 
      END IF
    ! Forcing= Qauto  
    ! Limiter= c(3)
    ! drhor  = RelCloud2*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter))
      drhor  = Qauto
    ! Forcing= Qauto/x_s 
    ! Limiter= c(5)
    ! dnr    = RelCloud2*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter))
      dnr    = Qauto/x_s
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhor,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1) =-drhor   ! QC
      f(2) = drhor   ! QR
      f(3) =-dnr     ! NC
      f(4) = dnr     ! NR
      f(5) = RhsPotM ! TH
    END IF
  CASE('SELFC') ! Selfcollection cloud droplets (c + c -> c) 
    Rho    = c(1)
    RhoC   = c(2)
    NC     = c(3)
    f(1)   = Zero
    sc=Zero
    IF (RhoC>1.0d-10) THEN
      sc     = p_cloud_k_sc*rhoc*rhoc*rho0/(rho+Eps) 
    ! forcing=-sc 
    ! limiter= nc  
    ! dnc    = RelCloud2*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter))
      dnc    = -sc
      f(1)   = dnc ! NC 
    END IF
  CASE('SELFR') ! Selfcollection rain (r + r -> r)
    Rho    = c(1)
    RhoR   = c(2)
    NR     = c(3)
    f(1)   = Zero
    IF (RhoR>1.0d-10) THEN
      x_r=MAX(xr_min,MIN(xr_min,RhoR/(NR+Eps)))
      N_0=MAX(N0_min,MIN(N0_max,NR*(Pi*rho0_liq/x_r)**(1./3.)))
      lambda_r=MAX(lambda_min,MIN(lambda_max,(Pi*rho0_liq*N_0/RhoR)**(1./4.)))
      sc     = k_rr*Nr*RhoR*(One+60.7d0/lambda_r)**(-9.0d0)*SQRT(Rho0/(Rho+Eps)) 
    ! forcing=-sc 
    ! limiter= nr 
    ! dnr    = RelCloud2*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter)) 
      dnr    = -sc
      f(1)   = dnr ! NR
    END IF
  CASE('ACC') ! Accretion (r + c -> r)
    Rho    = c(1)
    RhoV   = MAX(c(2),Zero)
    RhoC   = MAX(c(3),Zero)
    RhoR   = MAX(c(4),Zero)
    NC     = MAX(c(5),Zero)
    PotM   = c(6)
    nSpecies=6
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    f(5)=Zero
    RhoI=c(29)
    NI  =c(30)
    RhoS=c(31)
    NS  =c(32)
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    Qacc=Zero
    IF (RhoC>1.0d-10.AND.RhoR>1.0d-10) THEN
      Tau=min(max(One-RhoC/(RhoC+RhoR+Eps),Eps),One-1.0d-1)
      Phi=(Tau/(Tau+5.0d-5))**4.0d0
      Qacc=k_cr*RhoC*RhoR*Phi*SQRT(Rho0/(Rho+Eps)) 
    ! Forcing= Qacc
    ! Limiter= RhoC
    ! drhor  = RelCloud2*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter)) 
      drhor  = Qacc
    ! Forcing= Qacc/x_s
    ! Limiter= Nc 
    ! dnr    = RelCloud2*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter))
      x_s    = p_cloud_x_max
      dnr    = Qacc/x_s
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhor,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1)   =-drhor   ! QC
      f(2)   = drhor   ! QR
      f(3)   =-dnr     ! NC
      f(4)   = dnr     ! NR
      f(5)   = RhsPotM ! TH
    END IF
  CASE('EVAP') ! Evaporation of rain (r + r -> v)
    Rho =c(1)
    RhoV=MAX(c(2),Zero)
    RhoC=MAX(c(3),Zero)
    RhoR=MAX(c(4),Zero)
    Nr  =MAX(c(5),Zero)
    PotM=c(6)
    nSpecies=6
    RhoI=c(29)
    NI  =c(30)
    RhoS=c(31)
    NS  =c(32)
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    f(5)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    pVs=SaturVapor(T)
    pV =RhoV*Rv*T
    S  =pV/(pVs+Eps)-One
    IF (S<-1.0d-3.AND.RhoR>1.0e-8.AND.Nr>Zero) THEN
      D_vtp = (8.7602d-5*T)**1.81d0/p
      Lv    = LatHeat(T)
      g_d   = 4.0d0*3.14d0/(Lv*Lv/(K_T*287.0d0*T*T)+287.0d0*T/(D_vtp*pVs)) 
      x_r   = min(max(RhoR/(Nr+1d-10),p_rain_x_min),p_rain_x_max) 
      d_r   = p_rain_a_geo*x_r**p_rain_b_geo 
      v_r   = p_rain_a_vel*x_r**p_rain_b_vel*(Rho0/(Rho+Eps))**0.5d0  ! corr-factor missing
      N_re  = v_r*d_r/nu_l 
      f_v   = N_sc**n_f*N_re**m_f
      f_q   = a_q+b_q*f_v 
      eva   = g_d*max(Nr,1.d3)*c_r*d_r*S 
      eva_q = min(-5d-9,100.d0*f_q*eva) ! <0
      eva_n = eva_q/x_r 
    ! forcing= eva_q 
    ! limiter= c(4) 
    ! drhor  = RelCloud2*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter)) 
      drhor  = eva_q
      drhov  = -drhor  
    ! forcing= eva_n 
    ! limiter= c(5)
    ! dnr    = RelCloud2*(forcing-limiter+sqrt(forcing*forcing+limiter*limiter)) 
      dnr    = eva_n
      dnv    = -dnr 
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhor,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1)   = drhov   ! QV >0
      f(2)   =-drhov   ! QR <0
      f(3)   = dnv     ! NV
      f(4)   =-dnv     ! NR
      f(5)   = RhsPotM ! TH
    END IF
! ICE MICROPHYSICS
  CASE('FREEZEC')!homogeneous and heterogeneous freezing of cloud droplets
    Rho =c(1)           !density
    RhoV=MAX(c(2),Zero) !water vapour content kg m^-3
    RhoC=MAX(c(3),Zero) !cloud water content kg m^-3
    RhoR=MAX(c(4),Zero) !rain water content kg m^-3
    RhoI=MAX(c(5),Zero) !ice content kg m^-3
    NC  =MAX(c(6),Zero)           !number density of cloud droplets
    NI  =MAX(c(7),Zero)           !number density of ice particle 
    PotM=c(8)           !Thermodynamic variable
    nSpecies=8 
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    f(5)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS                  !density of dry air 
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    IF (T<273.15d0.AND.RhoC>Zero.AND.NC>Zero.AND.(Time.GE.StartTimeIce))THEN 
      T_c=T-273.15d0
      IF(T_c<-50.0d0)THEN
        fr_q=RhoC                                             !..Komplettes hom. Gefrieren
        fr_n=NC                                               !..unterhalb -50 C
      ELSE
        x_c=MIN(MAX(RhoC/(NC+Eps),p_cloud_x_min),p_cloud_x_max) !Mittlere Masse
        !homgeneous freezing of cloud droplets after Cotton and Field (2002) equ. 12 in SI [1/m3 s]
        IF(T_c<=-30.0d0) THEN 
          J_hom=1.0d6*EXP(-243.4d0-14.75d0*T_c-0.307d0*T_c**2.0d0-0.00287d0*T_c**3.0d0-0.0000102d0*T_c**4.0d0)
        ELSE
          J_hom=1.0d6*EXP(-7.63d0-2.996d0*(T_c+30.0d0))                             
        END IF
        J_het = MAX(A_het*(EXP(B_het*(T_3-T)-One)),Zero) ! [1/kg]
        !J_het=Zero                               !bei Wolkentropfengefrieren nicht berücksichigt 
        J_tot=((J_hom/rho0_liq)+J_het)            ![kg^-1]
        facg=2.173d0                              !factor from gamma distribution!
        fr_q=J_tot*RhoC*x_c*facg                  !rate of change of mass densities (liquid water content)
        fr_n=J_tot*RhoC                           !rate of change of number distribution
      END IF
      fr_q=MIN(fr_q,RhoC)
      fr_n=MIN(fr_n,NC)

    ! forcing= fr_q
    ! limiter=RhoC
    ! drhoi  = RelCloudIce*(forcing+limiter-SQRT(forcing*forcing+limiter*limiter)) ! zeitliche Ableitung der Dichte von Eis 
      drhoi  = fr_q

    ! forcing= fr_n
    ! limiter= NC
    ! dni    = RelCloudIce*(forcing+limiter-SQRT(forcing*forcing+limiter*limiter)) ! zeitlich Ableitung der number density rate of ice particles
      dni    = fr_n
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhoi,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1) =-drhoi !RhoC 
      f(2) = drhoi !RhoI
      f(3) =-dni   !NC
      f(4) = dni   !NI 
      f(5) = RhsPotM
    END IF
  CASE('FREEZER')!heterogeneous freezing of rain droplets
    Rho =c(1)           !density
    RhoV=MAX(c(2),Zero) !water vapour content kg m^-3
    RhoC=MAX(c(3),Zero) !cloud water content kg m^-3
    RhoR=MAX(c(4),Zero) !rain water content kg m^-3
    RhoI=MAX(c(5),Zero)           !ice content kg m^-3
    NI  =c(6)           !number density of ice particle
    NR  =c(7)           !number density of rain droplets 
    RhoS=c(31)
    NS  =c(32)
    PotM=c(8)          !Temperature
    nSpecies=8
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    f(5)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    T_a=T
    IF(T_a<T_freeze.AND.RhoR>q_krit_fr.AND.(Time.GE.StartTimeIce))THEN !T_freeze = 273.15 K
    ! WRITE(*,*) 'Cloudice'
      x_r = MIN(MAX(RhoR/(NR+Eps),p_rain_x_min),p_rain_x_max)
      IF(T_a<T_f)THEN                       !T_f = 233.16 K, ist der absolute Gefrierpunkt von -40°C
        fr_q=RhoR
        fr_n=NR
      ELSE
        J_het = MAX(A_het*(EXP(B_het*(T_3-T_a)-One)),Zero) !*dtAct ! [1/kg]
        facg=20.0d0
        fr_q=facg*RhoR*x_r*J_het
        fr_n=NR*x_r*J_het
      END IF
      fr_q=MIN(fr_q,RhoR)
      fr_n=MIN(fr_n,NR)
      forcing=fr_q
      limiter=c(4) !RhoR
    ! drhoi  = RelCloudIce*(forcing+limiter-SQRT(forcing*forcing+limiter*limiter)) ! zeitliche Ableitung der Dichte von Eis
      drhoi  = fr_q

      forcing=fr_n
      limiter=c(7) !NR
    ! dni    = RelCloudIce*(forcing+limiter-SQRT(forcing*forcing+limiter*limiter)) ! zeitlich Ableitung der number density rate of ice particles
      dni    = fr_n
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhoi,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1) =-drhoi !RhoR 
      f(2) = drhoi !RhoI
      f(3) =-dni   !NR
      f(4) = dni   !NI 
      f(5) = RhsPotM
    END IF
  CASE('NUCLI') !Nucleation of ice crystals
    Rho =c(1)           !density
    RhoV=MAX(c(2),Zero) !water vapour content kg m^-3
    RhoC=MAX(c(3),Zero)
    RhoR=MAX(c(4),Zero)
    RhoI=MAX(c(5),Zero)           !ice content kg m^-3
    NI  =MAX(c(6),Zero)           !number density of ice particle
    RhoS=c(31)
    NS  =c(32)
    PotM=c(7)
    nSpecies=7
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    pVi=SaturVaporIce(T)                               !saturation vapour pressure over ice [Pa]
    S_i=(RhoV*Rv*T/(pVi+Eps))-One                    !supersaturation wrt. ice
    QcLoc=RhoC/(Rho+Eps)
    IF (N_IN_Fix) THEN !can be set to a constant value, like 0/1 or 4 L^-1 INP
      N_IN_diag=N_IN
    ELSE
    ! N_IN_diag=N_M92*EXP(a_M92+b_M92*S_i)                  ! Meyers et al. (1992)
      N_IN_diag=0.01d0*EXP(0.6d0*(273.15d0-MAX(T,246.0d0))) ! Reisner (1998)
    ! N_IN_diag=0.0000594d0*(273.16d0-T)**3.33d0*(N_tot_05*1.0d-6)**(0.0264d0*(273.16d0-T)+0.0033d0) ! DeMott et al. (2010), n_{a,0.5} in cm^-3 
    END IF
    IF((T<T_freeze).AND.(S_i.GE.S_i_min).AND.(NI<N_IN_diag).AND.(QcLoc.GE.IceNucQc).AND.(Time.GE.StartTimeIce)) THEN !.AND.(N_IN_M92>0.1*N_IN_Reis).AND.(N_IN_M92<10*N_IN_Reis)
      fr_n=(N_IN_diag-NI)*RelCloudIce
    ! forcing=fr_n
    ! limiter=N_IN_diag
    ! dni=RelCloudIce*(forcing+limiter-SQRT(forcing*forcing+limiter*limiter))
      dni=fr_n
      drhoi=dni*p_ice_x_min
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhoi,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1)=-drhoi !RhoV
      f(2)=drhoi  !RhoI
      f(3)=dni
      f(4)=RhsPotM
     !WRITE (*,*) 'NUCLI: ',T,S_i,NI,N_IN_diag,dni,drhoi,RhsPotM
    END IF
  CASE('ICEDEP') !growth of ice particles by water vapour deposition 
    Rho =c(1)           !density
    RhoV=MAX(c(2),Zero) !water vapour content kg m^-3
    RhoC=MAX(c(3),Zero) !cloud water content kg m^-3
    RhoR=MAX(c(4),Zero) !rain water content kg m^-3
    RhoI=MAX(c(6),Zero) !ice content kg m^-3
    NI  =MAX(c(5),Zero) !number density of ice particle
    PotM=c(7)           !thermodynamic variable
    nSpecies=7
    RhoS=c(31)
    NS  =c(32)
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    pVi=SaturVaporIce(T)                                                 !saturation vapour pressure over ice [Pa]
    S_i=((RhoV*Rv*T/(pVi+Eps))-(One))                                  !supersaturation wrt. ice
    IF(RhoI>Zero.AND.NI>Zero.AND.T<T_freeze.AND.(Time.GE.StartTimeIce)) THEN
      x_i=MIN(MAX(RhoI/(NI+Eps),p_ice_x_min),p_ice_x_max)         !mass of ice [kg]
      D_i=p_ice_a_geo*x_i**p_ice_b_geo                                     !diameter-mass relation [m]
      v_i=p_ice_a_vel_SB*x_i**p_ice_b_vel_SB*SQRT(Rho0/Rho+Eps)            !velocity-mass relation [m s^-1]
      G_iv=(((Rv*T+Eps)/(pVi*D_v))+(LatwaSUBL/(K_T*T+Eps))*((LatwaSUBL/(Rv*T+Eps))-One))**(-One)
      N_Re=(v_i*D_i)/nu_l                                                   !Reynoldsnumber
      F_vi=a_vent(One)+b_vent(One)*N_sc**n_f*N_Re**m_f
      fr_q=4.0d0*G_iv*D_i*F_vi*S_i*NI                                      !mass density of ice 

      drhoi=fr_q
    ! WRITE(*,*) 'ICEDEP',drhoi
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhoi,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1)=-drhoi !RhoV
      f(2)= drhoi !RhoI
      f(3)= RhsPotM
    ! WRITE (*,*) 'ICEDEP: ',T,S_i,RhoI,NI,N_IN_diag,drhoi,RhsPotM
    END IF
  CASE('MELTI') !This process describes the melting process of ice crystals
    Rho =c(1)           !density
    RhoV=MAX(c(2),Zero) !water vapour content kg m^-3
    RhoC=MAX(c(3),Zero) !cloud water content kg m^-3
    RhoR=MAX(c(4),Zero) !rain water content kg m^-3
    RhoI=MAX(c(7),Zero) !ice content kg m^-3
    NC  =MAX(c(5),Zero)
    NI  =MAX(c(6),Zero) !number density of ice particle
    PotM=c(8)           !thermodynamic variable
    nSpecies=8
    RhoS=c(31)
    NS  =c(32)
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    f(5)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    IF(T>T_freeze.AND.RhoI>Zero.AND.NI>Zero.AND.(Time.GE.StartTimeIce))THEN
      pVs=SaturVapor(T)
      x_i=MIN(MAX(RhoI/(NI+Eps),p_ice_x_min),p_ice_x_max)         !mass of ice [kg]
      D_i=p_ice_a_geo*x_i**p_ice_b_geo                                     !diameter-mass relation [m]
      v_i=p_ice_a_vel_SB*x_i**p_ice_b_vel_SB*SQRT(Rho0/Rho+Eps)            !velocity-mass relation [m s^-1]
      N_Re=(v_i*D_i)/nu_l                                                  !Reynoldsnumber
      D_T=K_T/(Cpair*Rho0)                                                 !diffusivity of heat 
      F_vi_0=a_vent(Zero)+b_vent(Zero)*N_sc**n_f*N_Re**m_f                              !averaged ventilation coefficient for n=0
      F_vi_1=a_vent(One)+b_vent(One)*N_sc**n_f*N_Re**m_f                              !averaged ventilation coefficient for n=1
      G_iv=(2*pi/LatwaFREEZE)*(((K_T*D_T*(T-T_freeze))/D_v)+((D_v*LatwaSUBL)/Rv)*((pVs/T)-(p0Star/T_freeze)))
      fr_q=G_iv*NI*D_i*F_vi_1  !*dtAct !-sign neglected
      fr_n=G_iv*NI*D_i*x_i**(-1)*F_vi_0  !*dtAct !-sign neglected
      ! IF(T>283.15)THEN
      !   fr_q=RhoI
      !   fr_n=NI
      ! END IF
    ! forcing=fr_q
    ! limiter=RhoI
    ! drhoc=RelCloudIce*(forcing+limiter-SQRT(forcing*forcing+limiter*limiter))
      drhoc=fr_q

    ! forcing=fr_n
    ! limiter=NI
    ! dnc=RelCloudIce*(forcing+limiter-SQRT(forcing*forcing+limiter*limiter))
      dnc=fr_n
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhoc,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1)=-drhoc !RhoI
      f(2)= drhoc !RhoC
      f(3)=-dnc   !NI 
      f(4)= dnc   !NC
      f(5)= RhsPotM
    END IF
  CASE('AUTOI')
    Rho =c(1)           !density
    RhoV=MAX(c(2),Zero) !water vapour content kg m^-3
    RhoC=MAX(c(3),Zero) !cloud water content kg m^-3
    RhoR=MAX(c(4),Zero) !rain water content kg m^-3
    RhoI=MAX(c(5),Zero) !ice content kg m^-3
    RhoS=MAX(c(7),Zero)  !snow content kg m^-3
    NI  =MAX(c(6),Zero) !number density of ice particle
    NS  =MAX(c(8),Zero)  !number density of snow crystals 
    PotM=c(9)           !thermodynamic variable
    nSpecies=9
    f(1)=Zero
    f(2)=Zero
    f(3)=Zero
    f(4)=Zero
    f(5)=Zero
    RhoD=Rho-RhoV-RhoC-RhoR-RhoI-RhoS
    Cpml=CpmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Cvml=CvmlFun(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS)
    Rm  =RmFun(RhoD,RhoV)
    CALL ComputeTp(T,p,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM,c,Task,nSpecies)
    x_i=MIN(MAX((RhoI/(NI+Eps)),p_ice_x_min),p_ice_x_max)         !mass of ice [kg]
    D_i=p_ice_a_geo*x_i**p_ice_b_geo
    IF(NI>Zero.AND.RhoI>q_krit_ii.AND.D_i>D_krit_ii.AND.(Time.GE.StartTimeIce))THEN
      IF(T>T_freeze)THEN
        e_coll=1.0
      ELSE
        !.. Temperaturabhaengige Efficiency nach Cotton et al. (1986) 
        !   (siehe auch Straka, 1989; S. 53)
        !  e_coll = MIN(10**(0.035*(T_a-T_3)-0.7),0.2d0)
        !.. Temperaturabhaengige Efficiency nach Lin et al. (1983)
        !e_coll = MIN(EXP(0.09*(T_a-T_freeze)),One)
        !e_coll = MAX(e_ii,MIN(exp(0.09*(T_a-T_3)),One))
        e_coll = MAX(0.1d0,MIN(EXP(0.09*(T-273.15)),One))
      END IF
      v_i=p_ice_a_vel_SB*x_i**p_ice_b_vel_SB*SQRT(Rho0/Rho+Eps)            !velocity-mass relation [m s^-1] 
      
      theta_n      = 2.0d0*theta_i(Zero)-theta_ii(Zero)
      theta_q      = theta_i(Zero)-theta_ii(One)+theta_i(One)
      delta_n      = 2.0d0*delta_i(Zero)+delta_ii(Zero)
      delta_q      = delta_i(Zero)+delta_ii(One)+delta_i(One)
      fr_n=-pi*0.25*e_coll*NI*NI  *D_i*D_i*delta_n*SQRT(v_i*v_i*theta_n+2*ice_s_vel*ice_s_vel) !*dtAct
      fr_q=-pi*0.25*e_coll*NI*RhoI*D_i*D_i*delta_q*SQRT(v_i*v_i*theta_q+2*ice_s_vel*ice_s_vel) !*dtAct

    ! forcing=fr_n
    ! limiter=c(6) !NI
    ! dni=RelCloudIce*(forcing-limiter+SQRT(forcing*forcing+limiter*limiter))
      dni=fr_n

    ! forcing=fr_q
    ! limiter=c(5) !RhoI
    ! drhoi=RelCloudIce*(forcing-limiter+SQRT(forcing*forcing+limiter*limiter))
      drhoi=fr_q
      CALL ComputeRhsPotM(c,RhsPotM,PotM,drhoi,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1)=drhoi      !RhoI
      f(2)=-drhoi     !RhoS
      f(3)=dni        !NI 
      f(4)=-dni   !NS
      f(5)= RhsPotM
    END IF
  CASE DEFAULT
    IF (MyId==0) THEN
      WRITE(*,*) '  Error in .chem file - wrong process: ThetaKind :: ', TRIM(ThetaKind)
      WRITE(*,*) '  Stopped in TwoMomentBulkMicroPhysics_Mod.'
      STOP  
    END IF
  END SELECT

END SUBROUTINE MicroProcess2

SUBROUTINE BulkMicroJac2(Gas,Jac,VelFace,dtAct) 
  TYPE(Vec4_T) :: Gas(:)
  TYPE(Vec4_T) :: Jac(:)
  TYPE(VelocityFace_T) :: VelFace(:)
  REAL(RealKind) :: dtAct,Time

  INTEGER :: ix,iy,iz
  INTEGER :: istr,k,l,idf,iSpecies,iReak
  REAL(RealKind) :: Temp,Df
  REAL(RealKind) :: c(64),f(30),fC(30)
  TYPE(Reaction_T), POINTER :: Reaction

  REAL(RealKind) :: Rho,RhoD,RhoV,RhoC,RhoR,RhoI,RhoS,T,p,PotM,RelHum

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Reaction=>MicroFirst
        DO iReak=1,nReakMicro 
          istr=0
          DO k=1,Reaction%NumSpeciesLeftAktiv
            iSpecies=Reaction%SpeciesLeft(k)
            c(k)=Gas(iSpecies)%c(ix,iy,iz,1)
            IF (iSpecies==RhoPos) THEN
              c(k)=c(k)+RhoProfG(ibLoc)%c(ix,iy,iz,1) 
            END IF
          END DO
          c(Reaction%NumSpeciesLeftAktiv+1)=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+2)=PreCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+3)=KinEnCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+4)=zP(iz-1)+zP(iz)
          IF (Reaction%TypeR=='NUCL'.OR.Reaction%TypeR=='NUCL2') THEN
            IF (iz/=iz1) THEN
              c(Reaction%NumSpeciesLeftAktiv+5) =Gas(RhoPos)%c(ix,iy,iz+1,1) ! Rhozr
              c(Reaction%NumSpeciesLeftAktiv+6) =Gas(RhoVPos)%c(ix,iy,iz+1,1)  ! RhoVzr
              c(Reaction%NumSpeciesLeftAktiv+7) =Gas(RhoCPos)%c(ix,iy,iz+1,1)  ! RhoCzr
              c(Reaction%NumSpeciesLeftAktiv+8) =Gas(RhoRPos)%c(ix,iy,iz+1,1)  ! RhoRzr
              c(Reaction%NumSpeciesLeftAktiv+9) =Gas(thPos)%c(ix,iy,iz+1,1)  ! RhoRzr
            ELSE
              c(Reaction%NumSpeciesLeftAktiv+5) =Gas(RhoPos)%c(ix,iy,iz,1)   ! Rhozr
              c(Reaction%NumSpeciesLeftAktiv+6) =Gas(RhoVPos)%c(ix,iy,iz,1)    ! RhoVzr
              c(Reaction%NumSpeciesLeftAktiv+7) =Gas(RhoCPos)%c(ix,iy,iz,1)    ! RhoCzr
              c(Reaction%NumSpeciesLeftAktiv+8) =Gas(RhoRPos)%c(ix,iy,iz,1)    ! RhoRzr
              c(Reaction%NumSpeciesLeftAktiv+9) =Gas(thPos)%c(ix,iy,iz,1)    ! RhoRzr
            END IF
            c(Reaction%NumSpeciesLeftAktiv+10)=VelFace(ibLoc)%wF(ix,iy,iz)      ! wLoc
            c(Reaction%NumSpeciesLeftAktiv+11)=zP(iz)-zP(iz-1)               ! dz
          END IF
          c(21) = Gas(RhoPos)%c(ix,iy,iz,1)
          c(22) = Gas(RhoVPos)%c(ix,iy,iz,1)
          c(23) = Gas(RhoCPos)%c(ix,iy,iz,1)
          c(24) = Gas(RhoRPos)%c(ix,iy,iz,1)
          c(25) = Gas(nvPos)%c(ix,iy,iz,1)
          c(26) = Gas(ncPos)%c(ix,iy,iz,1)
          c(27) = Gas(nrPos)%c(ix,iy,iz,1)
          c(28) = Gas(thPos)%c(ix,iy,iz,1)
          c(29) = Gas(RhoIPos)%c(ix,iy,iz,1)
          c(30) = Gas(niPos)%c(ix,iy,iz,1)
          c(31) = Gas(RhoSPos)%c(ix,iy,iz,1)
          c(32) = Gas(nsPos)%c(ix,iy,iz,1)
          Rho = Gas(RhoPos)%c(ix,iy,iz,1)
          RhoV= Gas(RhoVPos)%c(ix,iy,iz,1)
          RhoC= Gas(RhoCPos)%c(ix,iy,iz,1)
          RhoR= Gas(RhoRPos)%c(ix,iy,iz,1)
          RhoI= Gas(RhoIPos)%c(ix,iy,iz,1)
          RhoS= Gas(RhoSPos)%c(ix,iy,iz,1)
          PotM= Gas(thPos)%c(ix,iy,iz,1)
          RhoD=Rho-RhoV-RhoC-RhoR
          p=PressureTheta(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM)+Eps
          T=AbsTemp(RhoD,RhoV,p)+Eps
          RelHum=RhoV*Rv*T/(SaturVapor(T)+Eps)
          IF (RhoC+RhoR>1.d-10.OR.RelHum>=One.OR.RhoI+RhoS>1.d-10) THEN
            CALL MicroProcess2(c,fC,Reaction%TypeR,'Function',dtAct,Time)
            DO k=1,Reaction%NumSpeciesLeftAktiv
              IF (ABS(c(k))<Eps) THEN
                Temp=c(k)
                c(k)=c(k)*(One+1.e-8_RealKind)+SIGN(1.e-8_RealKind,c(k))
                CALL MicroProcess2(c,f,Reaction%TypeR,'Jacobi  ',dtAct,Time)
                DO l=1,Reaction%NumSpecies
                  istr=istr+1
                  idf=Reaction%Struct(istr)
                  Df=(f(l)-fC(l))/(c(k)-Temp+Eps)
                  IF (ABS(fC(l))<=1.0d-10) THEN
                    Df=Zero 
                  END IF
                  Jac(idf)%c(ix,iy,iz,1)=Jac(idf)%c(ix,iy,iz,1)+Df
                END DO
                c(k)=Temp 
              END IF
            END DO
          END IF
          Reaction=>Reaction%Next
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE BulkMicroJac2

SUBROUTINE BulkMicro2(Gas,fGas,VelFace,dtAct,Time) 

  TYPE(Vec4_T) :: Gas(:),fGas(:)
  TYPE(VelocityFace_T) :: VelFace(:)
  REAL(RealKind) ::dtAct,Time

  INTEGER :: ix,iy,iz
  INTEGER :: k,iSpecies,iReak
  TYPE(Reaction_T), POINTER :: Reaction
  REAL(RealKind) :: c(64),f(30),Rho
  LOGICAL :: printf
  REAL(RealKind) :: RhoD,RhoV,RhoC,RhoR,RhoI,RhoS,T,p,PotM,RelHum

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        Reaction=>MicroFirst
        DO iReak=1,nReakMicro 
          DO k=1,Reaction%NumSpeciesLeftAktiv
            iSpecies=Reaction%SpeciesLeft(k)
            c(k)=Gas(iSpecies)%c(ix,iy,iz,1)
            IF (iSpecies==RhoPos) THEN
              c(k)=c(k)+RhoProfG(ibLoc)%c(ix,iy,iz,1) 
            END IF
          END DO
          c(Reaction%NumSpeciesLeftAktiv+1)=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+2)=PreCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+3)=KinEnCell(ibLoc)%c(ix,iy,iz,1)
          c(Reaction%NumSpeciesLeftAktiv+4)=zP(iz-1)+zP(iz)
          IF (Reaction%TypeR=='NUCL'.OR.Reaction%TypeR=='NUCL2') THEN
            IF (iz/=iz1) THEN
              c(Reaction%NumSpeciesLeftAktiv+5) =Gas(RhoPos)%c(ix,iy,iz+1,1) ! Rhozr
              c(Reaction%NumSpeciesLeftAktiv+6) =Gas(RhoVPos)%c(ix,iy,iz+1,1)  ! RhoVzr
              c(Reaction%NumSpeciesLeftAktiv+7) =Gas(RhoCPos)%c(ix,iy,iz+1,1)  ! RhoCzr
              c(Reaction%NumSpeciesLeftAktiv+8) =Gas(RhoRPos)%c(ix,iy,iz+1,1)  ! RhoRzr
              c(Reaction%NumSpeciesLeftAktiv+9) =Gas(thPos)%c(ix,iy,iz+1,1)  ! Thzr
            ELSE
              c(Reaction%NumSpeciesLeftAktiv+5) =Gas(RhoPos)%c(ix,iy,iz,1)   ! Rho
              c(Reaction%NumSpeciesLeftAktiv+6) =Gas(RhoVPos)%c(ix,iy,iz,1)    ! RhoV
              c(Reaction%NumSpeciesLeftAktiv+7) =Gas(RhoCPos)%c(ix,iy,iz,1)    ! RhoC
              c(Reaction%NumSpeciesLeftAktiv+8) =Gas(RhoRPos)%c(ix,iy,iz,1)    ! RhoR
              c(Reaction%NumSpeciesLeftAktiv+9) =Gas(thPos)%c(ix,iy,iz,1)    ! Th
            END IF
            c(Reaction%NumSpeciesLeftAktiv+10)=VelFace(ibLoc)%wF(ix,iy,iz)      ! wLoc
            c(Reaction%NumSpeciesLeftAktiv+11)=zP(iz)-zP(iz-1)               ! dz
          END IF
          c(21) = Gas(RhoPos)%c(ix,iy,iz,1)
          c(22) = Gas(RhoVPos)%c(ix,iy,iz,1)
          c(23) = Gas(RhoCPos)%c(ix,iy,iz,1)
          c(24) = Gas(RhoRPos)%c(ix,iy,iz,1)
          c(25) = Gas(nvPos)%c(ix,iy,iz,1)
          c(26) = Gas(ncPos)%c(ix,iy,iz,1)
          c(27) = Gas(nrPos)%c(ix,iy,iz,1)
          c(28) = Gas(thPos)%c(ix,iy,iz,1)
          c(29) = Gas(RhoIPos)%c(ix,iy,iz,1)
          c(30) = Gas(niPos)%c(ix,iy,iz,1)
          c(31) = Gas(RhoSPos)%c(ix,iy,iz,1)
          c(32) = Gas(nsPos)%c(ix,iy,iz,1)
          Rho = Gas(RhoPos)%c(ix,iy,iz,1)
          RhoV= Gas(RhoVPos)%c(ix,iy,iz,1)
          RhoC= Gas(RhoCPos)%c(ix,iy,iz,1)
          RhoR= Gas(RhoRPos)%c(ix,iy,iz,1)
          RhoI= Gas(RhoIPos)%c(ix,iy,iz,1)
          RhoS= Gas(RhoSPos)%c(ix,iy,iz,1)
          PotM= Gas(thPos)%c(ix,iy,iz,1)
          RhoD=Rho-RhoV-RhoC-RhoR
          p=PressureTheta(RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,PotM)+Eps
          T=AbsTemp(RhoD,RhoV,p)+Eps
          RelHum=RhoV*Rv*T/(SaturVapor(T)+Eps)
          IF (RhoC+RhoR>1.d-10.OR.RelHum>=One.OR.RhoI+RhoS>1.d-10) THEN
            CALL MicroProcess2(c,f,Reaction%TypeR,'Function',dtAct,Time)
            DO k=1,Reaction%NumSpecies
              iSpecies=Reaction%Species(k)%Species
              fGas(iSpecies)%c(ix,iy,iz,1)=fGas(iSpecies)%c(ix,iy,iz,1)+f(k)
            END DO
          END IF
          Reaction=>Reaction%Next 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE BulkMicro2

END MODULE TwoMomentBulkMicroPhysics_Mod
