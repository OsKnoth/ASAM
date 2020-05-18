MODULE ISDACMicroPhysics_Mod

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

SUBROUTINE ComputeRhsPotM(RhsPotM,PotM,QRhs,RhoD,RhoV,RhoC,RhoI,Rm,Cpml,Cvml,T,p,Process)

  REAL(RealKind) :: RhsPotM,PotM,QRhs,RhoD,RhoV,RhoC,RhoI,Rm,Cpml,Cvml,T,p
  CHARACTER*8    :: Process
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

SUBROUTINE MicroProcessISDAC(c,f,Process,Task,dtAct,Time,lprint) 

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

! Ice deposition parameters
  REAL(RealKind) :: Dep_B,Dep_C,Dep_Dv,Dep_K,Dep_ac,Dep_bc,dmi_dt

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
    CALL ComputeRhsPotM(RhsPotM,PotM,Qcond,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
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
  !  IF (S>0.0d0) WRITE(*,*) 'Nucleation ',S,dSdzW,T,wLoc
    IF (S>0.0d0.AND.dSdzw>0.0d0.AND.T>233.15d0.AND.wLoc>0.0d0) THEN
      kappa_nuc = 0.462d0 
      C_ccn     = Nv
!     aerosol activation rate
      Qnucl=C_ccn*kappa_nuc*MIN(S,1.1d0)**(kappa_nuc-One)*dSdzw
    ! forcing=Qnucl 
    ! limiter=Nv 
    ! dnc    =RelCloud2*(forcing+limiter-sqrt(forcing*forcing+limiter*limiter)) 
      dnc    =Qnucl
      drhoc  =dnc*1d-12 ! 1d-12 is smallest drop mass
      CALL ComputeRhsPotM(RhsPotM,PotM,drhoc,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1) =-drhoc   ! QV
      f(2) = drhoc   ! QC
      f(3) =-dnc     ! NV
      f(4) = dnc     ! NC
      f(5) = RhsPotM ! TH
      IF (ABS(drhoc)>0.0d0) WRITE(*,*) 'ISDAC',drhoc,dnc
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
    S_i=(RhoV*Rv*T/(pVi+Eps))-1.0d0                    !supersaturation wrt. ice
    QcLoc=RhoC/(Rho+Eps)
    IF (N_IN_Fix) THEN !can be set to a constant value, like 0/1 or 4 L^-1 INP
      N_IN_diag=N_IN
    ELSE
      N_IN_M92=N_M92*EXP(a_M92+b_M92*S_i)                !deposition freezing after Meyers et al. (1992) 
      N_IN_diag=N_IN_M92                                 
    ! N_IN_Reis=0.01d0*EXP(-MIN(T,246.15)-273.15d0)      !that's the Reisner parameterization and after Seifert and Beheng the N_IN_M92 is bounded by 0.1 and 10 times N_IN_Reis
    ! N_IN_Reis=0.01d0*EXP(0.6d0*(273.15d0-MAX(T,246.0d0)))  ! Reisner (1998)
    END IF
    IF((T<T_freeze).AND.(S_i.GE.S_i_min).AND.(NI<N_IN_diag).AND.(QcLoc.GE.IceNucQc).AND.(Time.GE.StartTimeIce)) THEN !.AND.(N_IN_M92>0.1*N_IN_Reis).AND.(N_IN_M92<10*N_IN_Reis)
      fr_n=(N_IN_diag-NI)*RelCloudIce
    ! forcing=fr_n
    ! limiter=N_IN_diag
    ! dni=RelCloudIce*(forcing+limiter-SQRT(forcing*forcing+limiter*limiter))
      dni=fr_n
      drhoi=dni*p_ice_x_min
      CALL ComputeRhsPotM(RhsPotM,PotM,drhoi,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
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
    S_i=((RhoV*Rv*T/(pVi+Eps))-(1.0d0))                                  !supersaturation wrt. ice
    IF(RhoI>Zero.AND.NI>Zero.AND.T<T_freeze.AND.(Time.GE.StartTimeIce)) THEN
      x_i=MIN(MAX(RhoI/(NI+Eps),p_ice_x_min),p_ice_x_max)         !mass of ice [kg]
   !  D_i=p_ice_a_geo*x_i**p_ice_b_geo
      D_i=0.09d0*x_i**(1.d0/3.d0)
      Dep_Dv=3.0d-5 ! [m2 s-1] - Diffusivity coefficient of water vapor in air
      Dep_K=2.5d-2  ! [J m s-1 K-1] - coefficient of air heat conductivity
      Dep_B=One/(Rv*T/(pVi*Dep_Dv)+Lhs/(Dep_K*T)*(Lhs/(Rv*T)-One)) 
   !  Dep_ac=0.09d0
   !  Dep_C=Dep_ac*RhoI**(1.0d0/3.0d0)
      drhoi=Two*Pi*Dep_B*S_i*D_i*NI
   !  WRITE(*,*) 'ICEDEP',drhoi,Dep_B
      CALL ComputeRhsPotM(RhsPotM,PotM,drhoi,RhoD,RhoV,RhoC+RhoR,RhoI+RhoS,Rm,Cpml,Cvml,T,p,Process)
      f(1)=-drhoi !RhoV
      f(2)= drhoi !RhoI
      f(3)= RhsPotM
    ! WRITE (*,*) 'ICEDEP: ',T,S_i,RhoI,NI,N_IN_diag,drhoi,RhsPotM
    END IF
  CASE DEFAULT
  IF (MyId==0) THEN
    WRITE(*,*) 'Error in .chem file - wrong process. Stopped.'
    STOP  
  END IF
  END SELECT

END SUBROUTINE MicroProcessISDAC

SUBROUTINE BulkMicroJacISDAC(Gas,Jac,VelFace,dtAct) 
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
          Rho = Gas(rhoPos)%c(ix,iy,iz,1)
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
            CALL MicroProcessISDAC(c,fC,Reaction%TypeR,'Function',dtAct,Time)
            DO k=1,Reaction%NumSpeciesLeftAktiv
              IF (ABS(c(k))<Eps) THEN
                Temp=c(k)
                c(k)=c(k)*(1.0d0+1.e-8_RealKind)+SIGN(1.e-8_RealKind,c(k))
                CALL MicroProcessISDAC(c,f,Reaction%TypeR,'Jacobi  ',dtAct,Time)
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
END SUBROUTINE BulkMicroJacISDAC

SUBROUTINE BulkMicroISDAC(Gas,fGas,VelFace,dtAct,Time) 

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
            CALL MicroProcessISDAC(c,f,Reaction%TypeR,'Function',dtAct,Time)
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

END SUBROUTINE BulkMicroISDAC

END MODULE ISDACMicroPhysics_Mod
