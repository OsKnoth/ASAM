MODULE SpecialOutput_Mod

  !USE Control_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Output_Mod
  USE Operator_Mod
  USE WindFarm_Mod
  USE MicroPhysics_Mod
  USE Names_Mod, only: ql,qd

  IMPLICIT NONE

  INTEGER :: PosOutSh=0
  INTEGER :: PosOutRi=0
  INTEGER :: PosOutThMProf=0
  INTEGER :: PosOutRhoMProf=0
  INTEGER :: PosOutThPseudo=0
  INTEGER :: PosOutAbsTemp=0           ! absolute temperature [K]
  INTEGER :: PosOutIntEne=0            ! internal energy [J]  
  INTEGER :: PosOutTotalEne=0          ! total energy [J] 
  INTEGER :: PosOutEnergyDiff=0        ! difference between current and initial total energy [J] 
  INTEGER :: PosOutPreDiff=0           ! difference between current and initial pressure [Pa] 
  INTEGER :: PosOutThetaDry=0          ! dry potential temperature [K] 
  INTEGER :: PosOutDivergence=0        ! Divergence
  INTEGER :: PosOutxVorticity=0        ! x1 component of vorticity
  INTEGER :: PosOutyVorticity=0        ! x2 component of vorticity
  INTEGER :: PosOutzVorticity=0        ! x3 component of vorticity
  INTEGER :: PosOutVorticity=0         ! vorticity
  INTEGER :: PosOutNonDimStreamVel=0   ! nondimensional StreamVel
  INTEGER :: PosOutxzStress=0          ! Stress in 1,3
  INTEGER :: PosOutLen_turb=0          ! Csmag in dynamic Smagorinsky
  INTEGER :: PosOutAnelaP=0            ! Anelastic
  INTEGER :: PosOutRhoVMProf=0
  INTEGER :: PosOutThDens=0            ! density potential temperature [K] 
  INTEGER :: PosOutThEquiv=0           ! equivalent potential temperature [K] 
  INTEGER :: PosOutThetaL=0            ! liquid-water potential temperature [K] 
  INTEGER :: PosOutRelHum=0            ! relative humidity wrt. water [%] 
  INTEGER :: PosOutRelHumIce=0         ! relative humidity wrt. ice [%] 
  INTEGER :: PosOutRhoVS=0             ! density of water vapor at saturation [kg m^-3] 
  INTEGER :: PosOutRhoVOut=0           ! density of water vapor [kg m^-3] 
  INTEGER :: PosOutVPD=0               ! vapor pressure deficit [Pa] 
  INTEGER :: PosOutSuperSat=0          ! Supersaturation wrt. water [%]
  INTEGER :: PosOutSuperSatIce=0       ! Supersaturation wrt. ice [%]
  INTEGER :: PosOutQt=0                ! total water content:  qt = RhoV + ql + RhoR [kg kg^(-1)]
  INTEGER :: PosOutQLR=0               ! liquid water content: qlr = ql + RhoR [kg kg^(-1)]
  INTEGER :: PosOutLWC=0               ! Liquid water content [kg m^(-3)]
  INTEGER :: PosOutLWP=0               ! Liquid water path [kg m^(-2)]                
  INTEGER :: PosOutIWC=0               ! ice water content [kg m^(-3)]
  INTEGER :: PosOutIWP=0               ! ice water path [kg m^(-2)]
  INTEGER :: PosOutHeatRate=0          ! Heating Rate [K s^(-1)]
  INTEGER :: PosOutThForcing=0         ! Large-scale forcing term on density potential temperature [K s^(-1)]
  INTEGER :: PosOutRhoVForcing=0         ! Large-scale forcing term on specific humidity [kg kg^(-1) s^(-1)]
  INTEGER :: PosOutCpuNumb=0           ! Cpu_Number
  INTEGER :: PosOutBlkNumb=0           ! Cpu_Number
  INTEGER :: PosOutRhoPHeight=0
  INTEGER :: PosOutCdWind=0
  INTEGER :: PosOutCdTree=0            ! drag force of trees on the mean flow, CM
  INTEGER :: PosOutTreePor=0            ! tree porosity, CM
  INTEGER :: PosOutPM1=0               ! Aerosol conentration in PM1 at ambient conditions
  INTEGER :: PosOutPM2_5=0             ! Aerosol conentration in PM2.5 at ambient conditions
  INTEGER :: PosOutPM10=0              ! Aerosol conentration in PM10 at ambient conditions
  INTEGER :: PosOutVol=0               ! Volume per grid cell

CONTAINS


SUBROUTINE InitSpecialOutput

  INTEGER :: i

  DO i=1, LenOutSpecial
    IF (TRIM(NameOutSpecial(i))=='ShProd')          PosOutSh=i
    IF (TRIM(NameOutSpecial(i))=='RichNum')         PosOutRi=i
    IF (TRIM(NameOutSpecial(i))=='ThMProf')         PosOutThMProf=i
    IF (TRIM(NameOutSpecial(i))=='RhoMProf')        PosOutRhoMProf=i
    IF (TRIM(NameOutSpecial(i))=='ThPseudo')        PosOutThPseudo=i
    IF (TRIM(NameOutSpecial(i))=='AbsTemp')         PosOutAbsTemp=i          ! absolute temperature [K]
    IF (TRIM(NameOutSpecial(i))=='IntEne')          PosOutIntEne=i           ! internal energy [J]  
    IF (TRIM(NameOutSpecial(i))=='TotalEne')        PosOutTotalEne=i         ! total energy [J] 
    IF (TRIM(NameOutSpecial(i))=='EnergyDiff')      PosOutEnergyDiff=i       ! difference between current and initial total energy [J] 
    IF (TRIM(NameOutSpecial(i))=='PreDiff')         PosOutPreDiff=i          ! difference between current and initial pressure [Pa] 
    IF (TRIM(NameOutSpecial(i))=='ThetaDry')        PosOutThetaDry=i         ! dry potential temperature [K] 
    IF (TRIM(NameOutSpecial(i))=='Divergence')      PosOutDivergence=i       ! Divergence
    IF (TRIM(NameOutSpecial(i))=='xVorticity')      PosOutxVorticity=i       ! x1 component of vorticity
    IF (TRIM(NameOutSpecial(i))=='yVorticity')      PosOutyVorticity=i       ! x2 component of vorticity
    IF (TRIM(NameOutSpecial(i))=='zVorticity')      PosOutzVorticity=i       ! x3 component of vorticity
    IF (TRIM(NameOutSpecial(i))=='Vorticity')       PosOutVorticity=i        ! vorticity
    IF (TRIM(NameOutSpecial(i))=='NonDimStreamVel') PosOutNonDimStreamVel=i  ! nondimensional StreamVel
    IF (TRIM(NameOutSpecial(i))=='xzStress')        PosOutxzStress=i         ! Stress in 1,3
    IF (TRIM(NameOutSpecial(i))=='Len_turb')        PosOutLen_turb=i         ! Csmag in dynamic Smagorinsky
    IF (TRIM(NameOutSpecial(i))=='AnelaP')          PosOutAnelaP=i           ! Anelastic
    IF (TRIM(NameOutSpecial(i))=='RhoVMProf')         PosOutRhoVMProf=i
    IF (TRIM(NameOutSpecial(i))=='ThDens')          PosOutThDens=i           ! density potential temperature [K] 
    IF (TRIM(NameOutSpecial(i))=='ThEquiv')         PosOutThEquiv=i          ! equivalent potential temperature [K] 
    IF (TRIM(NameOutSpecial(i))=='ThetaL')          PosOutThetaL=i           ! liquid-water potential temperature [K] 
    IF (TRIM(NameOutSpecial(i))=='RelHum')          PosOutRelHum=i           ! relative humidity wrt. water [%] 
    IF (TRIM(NameOutSpecial(i))=='RelHumIce')       PosOutRelHumIce=i        ! relative humidity wrt. ice [%] 
    IF (TRIM(NameOutSpecial(i))=='RhoVS')           PosOutRhoVS=i            ! density of water vapor at saturation [kg m^-3] 
    IF (TRIM(NameOutSpecial(i))=='RhoV')            PosOutRhoVOut=i          ! density of water vapor [kg m^-3] 
    IF (TRIM(NameOutSpecial(i))=='VPD')             PosOutVPD=i              ! vapor pressure deficit [Pa] 
    IF (TRIM(NameOutSpecial(i))=='SuperSat')        PosOutSuperSat=i         ! Supersaturation wrt. water [%]
    IF (TRIM(NameOutSpecial(i))=='SuperSatIce')     PosOutSuperSatIce=i      ! Supersaturation wrt. ice [%]
    IF (TRIM(NameOutSpecial(i))=='Qt')              PosOutQt=i               ! total water content:  qt = RhoV + ql + RhoR [kg kg^(-1)]
    IF (TRIM(NameOutSpecial(i))=='QLR')             PosOutQLR=i              ! liquid water content: qlr = ql + RhoR [kg kg^(-1)]
    IF (TRIM(NameOutSpecial(i))=='LWC')             PosOutLWC=i              ! Liquid water content [kg m^(-3)]
    IF (TRIM(NameOutSpecial(i))=='LWP')             PosOutLWP=i              ! Liquid water path [kg m^(-2)]
    IF (TRIM(NameOutSpecial(i))=='IWC')             PosOutIWC=i              ! ice water content [kg m^(-3)]
    IF (TRIM(NameOutSpecial(i))=='IWP')             PosOutIWP=i              ! ice water path [kg m^(-2)]
    IF (TRIM(NameOutSpecial(i))=='HeatingRate')     PosOutHeatRate=i         ! Heating Rate [K s^(-1)]
    IF (TRIM(NameOutSpecial(i))=='dThdt')           PosOutThForcing=i        ! Large-scale forcing term on density potential temperature [K s^(-1)]
    IF (TRIM(NameOutSpecial(i))=='dRhoVdt')           PosOutRhoVForcing=i        ! Large-scale forcing term on specific humidity [kg kg^(-1) s^(-1)]
    IF (TRIM(NameOutSpecial(i))=='CpuNumb')         PosOutCpuNumb=i          ! Cpu_Number
    IF (TRIM(NameOutSpecial(i))=='BlkNumb')         PosOutBlkNumb=i          ! Cpu_Number
    IF (TRIM(NameOutSpecial(i))=='RhoPHeight')      PosOutRhoPHeight=i
    IF (TRIM(NameOutSpecial(i))=='CdWind')          PosOutCdWind=i
    IF (TRIM(NameOutSpecial(i))=='CdTree')          PosOutCdTree=i
    IF (TRIM(NameOutSpecial(i))=='TreePor')         PosOutTreePor=i
    IF (TRIM(NameOutSpecial(i))=='PM1_aNUMBER')     PosOutPM1=i              ! Aerosol conentration in PM1 at ambient conditions   ! Write "PM1" in OutputNamefile!!
    IF (TRIM(NameOutSpecial(i))=='PM2.5_aNUMBER')   PosOutPM2_5=i            ! Aerosol conentration in PM2.5 at ambient conditions ! Write "PM2.5" in OutputNamefile!!
    IF (TRIM(NameOutSpecial(i))=='PM10_aNUMBER')    PosOutPM10=i             ! Aerosol conentration in PM10 at ambient conditions  ! Write "PM10" in OutputNamefile!!
    IF (TRIM(NameOutSpecial(i))=='VolCell')         PosOutVol=i              ! Volume per grid cell
  END DO

END SUBROUTINE InitSpecialOutput

FUNCTION OutputTimeCheck(Time)

  LOGICAL :: OutputTimeCheck
  REAL(RealKind) :: Time

  OutputTimeCheck=.FALSE.
  IF (Step==OutputStep.OR.Step==OutputEnd) THEN
    IF (Time>=OutputTimeStart) THEN
      OutputTimeCheck=.TRUE.
    END IF
  END IF
  IF (Time>=OutputTime) THEN
    OutputTimeCheck=.TRUE.
  END IF

END FUNCTION OutputTimeCheck

SUBROUTINE SpecialOutput(VecT,VecG,VelF,ActTime)

  TYPE(Vector4Cell_T) :: VecT(:)
  TYPE(Vector4Cell_T) :: VecG(:)
  TYPE(VelocityFace_T),POINTER :: VelF(:)
  REAL(RealKind) :: ActTime

  IF (Step==OutputStart) THEN
    CALL InitSpecialOutput
  ENDIF
  IF (OutputTimeCheck(ActTime)) THEN
    IF (PosOutSh>0 .OR. PosOutRi>0) CALL ShearProdCompute(VecT,OutSpecial)
    IF (PosOutThMProf>0)            CALL ComputeThMinusThProf(VecT,OutSpecial)
    IF (PosOutRhoMProf>0)           CALL ComputeRhoMinusRhoProf(VecT,OutSpecial)
    IF (PosOutThPseudo>0)           CALL ComputeThPseudo(VecT,OutSpecial)
    IF (PosOutAbsTemp>0)            CALL ComputeAbsTemp(VecT,OutSpecial)            ! absolute temperature [K]
    IF (PosOutIntEne>0)             CALL ComputeInterEnergy(VecT,OutSpecial)        ! internal energy [J]  
    IF (PosOutTotalEne>0)           CALL ComputeTotalEnergy(VecT,OutSpecial)        ! total energy [J] 
    IF (PosOutEnergyDiff>0)         CALL ComputeEnergyDiff(VecT,OutSpecial)         ! difference between current and initial total energy [J] 
    IF (PosOutPreDiff>0)            CALL ComputePreDiff(VecT,OutSpecial)            ! difference between current and initial pressure [Pa] 
    IF (PosOutThetaDry>0)           CALL ComputeThetaDry(VecT,OutSpecial)           ! dry potential temperature [K] 
    IF (PosOutDivergence>0)         CALL ComputeDivergence(VecT,VelF,OutSpecial)    ! Divergence
    IF (PosOutxVorticity>0)         CALL ComputexVorticity(VecT,OutSpecial)         ! x1 component of vorticity
    IF (PosOutyVorticity>0)         CALL ComputeyVorticity(VecT,OutSpecial)         ! x2 component of vorticity
    IF (PosOutzVorticity>0)         CALL ComputezVorticity(VecT,OutSpecial)         ! x3 component of vorticity
    IF (PosOutVorticity>0)          CALL ComputeVorticity(VecT,OutSpecial)          ! vorticity
    IF (PosOutNonDimStreamVel>0)    CALL ComputeNonDimStreamVel(VecT,OutSpecial)    ! nondimensional StreamVel
    IF (PosOutxzStress>0)           CALL ComputexzStress(VecT,OutSpecial)           ! Stress in 1,3
    IF (PosOutLen_turb>0)           CALL ComputeLen_turb(VecT,OutSpecial)           ! Csmag in dynamic Smagorinsky
    IF (PosOutAnelaP>0)             CALL ComputeAnelaP(VecG,OutSpecial)             ! Anelastic
    IF (PosOutRhoVMProf>0)            CALL ComputeRhoVMinusRhoVProf(VecT,OutSpecial)
    IF (PosOutThDens>0)             CALL ComputeThDens(VecT,OutSpecial)             ! density potential temperature [K] 
    IF (PosOutThEquiv>0)            CALL ComputeThEquiv(VecT,OutSpecial)            ! equivalent potential temperature [K] 
    IF (PosOutThetaL>0)             CALL ComputeThetaL(VecT,OutSpecial)             ! liquid-water potential temperature [K] 
    IF (PosOutRelHum>0)             CALL ComputeRelHum(VecT,OutSpecial,ActTime)     ! relative humidity wrt. water [%] 
    IF (PosOutRelHumIce>0)          CALL ComputeRelHumIce(VecT,OutSpecial,ActTime)  ! relative humidity wrt. ice [%] 
    IF (PosOutRhoVS>0)              CALL ComputeRhoVS(VecT,OutSpecial)              ! density of water vapor at saturation [kg m^-3] 
    IF (PosOutRhoVOut>0)            CALL ComputeRhoV(VecT,OutSpecial)               ! density of water vapor [kg m^-3] 
    IF (PosOutVPD>0)                CALL ComputeVPD(VecT,OutSpecial)                ! vapor pressure deficit [Pa] 
    IF (PosOutSuperSat>0)           CALL ComputeSuperSat(VecT,OutSpecial)           ! Supersaturation wrt. water [%]
    IF (PosOutSuperSatIce>0)        CALL ComputeSuperSatIce(VecT,OutSpecial)        ! Supersaturation wrt. ice [%]
    IF (PosOutQt>0)                 CALL ComputeQt(VecT,OutSpecial)                 ! total water content:  qt = RhoV + ql + RhoR [kg kg^(-1)]
    IF (PosOutQLR>0)                CALL ComputeQLR(VecT,OutSpecial)                ! liquid water content: qlr = ql + RhoR [kg kg^(-1)]
    IF (PosOutLWC>0)                CALL ComputeLWC(VecT,OutSpecial)                ! Liquid water content [kg m^(-3)]
    IF (PosOutLWP>0)                CALL ComputeLWP(VecT,OutSpecial)                ! Liquid water path [kg m^(-2)]
    IF (PosOutIWC>0)                CALL ComputeIWC(VecT,OutSpecial)                ! ice water content [kg m^(-3)]
    IF (PosOutIWP>0)                CALL ComputeIWP(VecT,OutSpecial)                ! ice water path [kg m^(-2)]
    IF (PosOutHeatRate>0)           CALL ComputeHeatRate(VecT,OutSpecial)           ! Heating Rate [K s^(-1)]
    IF (PosOutThForcing>0)          CALL ComputeThForcing(VecT,OutSpecial)          ! Large-scale forcing term on density potential temperature [K s^(-1)]
    IF (PosOutRhoVForcing>0)          CALL ComputeRhoVForcing(VecT,OutSpecial)          ! Large-scale forcing term on specific humidity [kg kg^(-1) s^(-1)]
    IF (PosOutCpuNumb>0)            CALL ComputeCpuNumb(VecT,OutSpecial)            ! Cpu_Number
    IF (PosOutBlkNumb>0)            CALL ComputeBlkNumb(VecT,OutSpecial)            ! Cpu_Number
    IF (PosOutVol>0)                CALL ComputeVolCell(VecT,OutSpecial)            ! Volume per grid cell
    IF (Shallow .AND. Height .AND. PosOutRhoPHeight>0) THEN
      CALL ComputeRhoPlusHeight(VecT,OutSpecial)
    END IF
    IF (Wind .AND. PosOutCdWind>0) THEN
      CALL WindCdCompute(VecT,OutSpecial)        ! Drag Coefficient of wind turbine
      !CALL WindPowerCompute(VecT,OutSpecial)     ! writes power of wind turbines at each time step in ASCI-file
    END IF  
    IF (PosOutCdTree>0) CALL TreeCdCompute(VecT,OutSpecial)                               ! Drag force of tree strips on the mean wind
    IF (PosOutTreePor>0) CALL TreePorCompute(VecT,OutSpecial)                             ! tree porosity
    IF (PosOutPM1>0) CALL ComputeParticulateMatter(VecT,OutSpecial,PosOutPM1,1.d-6)       ! Aerosol conentration in PM1 at ambient conditions
    IF (PosOutPM2_5>0) CALL ComputeParticulateMatter(VecT,OutSpecial,PosOutPM2_5,2.5d-6)  ! Aerosol conentration in PM2.5 at ambient conditions
    IF (PosOutPM10>0) CALL ComputeParticulateMatter(VecT,OutSpecial,PosOutPM10,1.d-5)     ! Aerosol conentration in PM10 at ambient conditions
    !CALL ProfileVariablesCompute(VecT,OutSpecial,ActTime)
    !CALL RainRateCompute(VecT,ActTime)
    !IF (DynamicSoil.AND.Radiation) THEN
    !  CALL SoilFluxCompute(ActTime)
    !END IF
    !IF (Canopy.AND.Radiation) THEN
       !CALL CanopyFluxCompute(ActTime)
    !END IF
  END IF
    
END SUBROUTINE SpecialOutput


SUBROUTINE WindCdCompute(Vector,Out)

  TYPE(Vector4Cell_T) :: Vector(:),Out(:)
  
  INTEGER :: ibloc
  INTEGER :: i,ib,ix,iy,iz,iW,n
  REAL(RealKind) :: TempU,TempV,TempW,TempVM ! velocity components and mean velocity
  REAL(RealKind) :: Cd ! drag coefficient of the wind turbine 
  
  REAL(RealKind) :: R,Rad,d,Theta,Lambda,RXY
  
  DO ibloc=1,nbLoc
    ib=LocGlob(ibloc)
    Out(ibLoc)%Vec(PosOutCdWind)%c=0.0d0
    CALL Set(Floor(ib))
    DO iW=1,NumberW 
      n=WindTurbineInfluence(iW,ibloc)%Number
      DO i=1,n
        ix=WindTurbineInfluence(iW,ibloc)%index(1,i) 
        iy=WindTurbineInfluence(iW,ibloc)%index(2,i) 
        iz=WindTurbineInfluence(iW,ibloc)%index(3,i) 
         
        R=SQRT(WindTurbineInfluence(iW,ibloc)%RadiusVec(i,1)**2+WindTurbineInfluence(iW,ibloc)%RadiusVec(i,2)**2+ &
&         WindTurbineInfluence(iW,ibloc)%RadiusVec(i,3)**2)    !distance cell to wind turbine center, norm RVec
        Theta=WindTurbineInfluence(iW,ibloc)%Theta(i) !angle R to z-axes in yz plane
        Lambda=WindTurbineInfluence(iW,ibloc)%Lambda(i) !angle R to x axes in xy plane
        RXY=SQRT(WindTurbineInfluence(iW,ibloc)%RadiusVec(i,1)**2+WindTurbineInfluence(iW,ibloc)%RadiusVec(i,2)**2) !RVec proj xy plane
        d=SIN(WindFarm(iW)%Alpha-Lambda)*RXY !distance to WindTurbine plane 
     

        Rad=Windfarm(iW)%Rad !Rotor radius
    
        Cd=WindDragCoefficient(R,Rad,d,ABS(Theta-WindFarm(iW)%RotorPosition))
        Out(ibLoc)%Vec(PosOutCdWind)%c(ix,iy,iz,1)=Cd
      END DO    
    END DO
  END DO

END SUBROUTINE WindCdCompute

SUBROUTINE TreeCdCompute(Vector,Out)

  TYPE(Vector4Cell_T) :: Vector(:),Out(:)
  INTEGER :: NumPoints
  INTEGER :: i,ix,iy,iz
  TYPE(TreePoint_T), POINTER :: TreePoint(:)

  REAL(RealKind) :: uLocL,uLocR
  REAL(RealKind) :: vLocL,vLocR
  REAL(RealKind) :: wLocL,wLocR
  REAL(RealKind), POINTER :: RhoLoc
  REAL(RealKind) :: c
  REAL(RealKind) :: VAbsL,VAbsR,VAbs
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: uCL(:,:,:,:)
  REAL(RealKind), POINTER :: vCL(:,:,:,:)
  REAL(RealKind), POINTER :: wCL(:,:,:,:)
  REAL(RealKind), POINTER :: uCR(:,:,:,:)
  REAL(RealKind), POINTER :: vCR(:,:,:,:)
  REAL(RealKind), POINTER :: wCR(:,:,:,:)

  REAL(RealKind) :: RhsC

!  Outc=>Out(ibLoc)%Vec(PosOutCdTree)%c

  NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
  TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint
  DO i=1,NumPoints
    ix=TreePoint(i)%ix
    iy=TreePoint(i)%iy
    iz=TreePoint(i)%iz
    c=TreePoint(i)%c
!   Evaluate velocity tangential and normal  
    RhoLoc=Rho(ix,iy,iz,1)
    uLocL=uCL(ix,iy,iz,1)/(RhoLoc+Eps)
    uLocR=uCR(ix,iy,iz,1)/(RhoLoc+Eps)
    vLocL=vCL(ix,iy,iz,1)/(RhoLoc+Eps)
    vLocR=vCR(ix,iy,iz,1)/(RhoLoc+Eps)
    wLocL=wCL(ix,iy,iz,1)/(RhoLoc+Eps)
    wLocR=wCR(ix,iy,iz,1)/(RhoLoc+Eps)

    VAbsL=SQRT(uLocL*uLocL+vLocL*vLocL+wLocL*wLocL)
    VAbsR=SQRT(uLocR*uLocR+vLocR*vLocL+wLocR*wLocR)
    VAbs=SQRT(VAbsL*VAbsL+VAbsR*VAbsR)

    RhsC=VolC(ix,iy,iz)*c*VAbs
   ! vRhsL(ix,iy,iz,1)=VolC(ix,iy,iz)*c*vCL(ix,iy,iz,1)*VAbsL
   ! wRhsL(ix,iy,iz,1)=VolC(ix,iy,iz)*c*wCL(ix,iy,iz,1)*VAbsL
   ! uRhsR(ix,iy,iz,1)=VolC(ix,iy,iz)*c*uCR(ix,iy,iz,1)*VAbsR
   ! vRhsR(ix,iy,iz,1)=VolC(ix,iy,iz)*c*vCR(ix,iy,iz,1)*VAbsR
   ! wRhsR(ix,iy,iz,1)=VolC(ix,iy,iz)*c*wCR(ix,iy,iz,1)*VAbsR

    Out(ibLoc)%Vec(PosOutCdTree)%c(ix,iy,iz,1)=RhsC    

  END DO

END SUBROUTINE TreeCdCompute

SUBROUTINE TreePorCompute(Vector,Out)

  TYPE(Vector4Cell_T) :: Vector(:),Out(:)
  INTEGER :: NumPoints
  INTEGER :: i,ix,iy,iz,ib,ibloc
  TYPE(TreePoint_T), POINTER :: TreePoint(:)

  REAL(RealKind) :: c

  DO ibloc=1,nbLoc
    ib=LocGlob(ibloc)
    Out(ibLoc)%Vec(PosOutTreePor)%c=0.0d0
    CALL Set(Floor(ib))

    NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
    TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint

    DO i=1,NumPoints
      ix=TreePoint(i)%ix
      iy=TreePoint(i)%iy
      iz=TreePoint(i)%iz
      c=TreePoint(i)%c

      Out(ibLoc)%Vec(PosOutTreePor)%c(ix,iy,iz,1)=c
    END DO

  END DO

END SUBROUTINE TreePorCompute

SUBROUTINE WindPowerCompute(Vector,Out)

  TYPE(Vector4Cell_T) :: Vector(:),Out(:)
  
  REAL(RealKind) :: Power(NumberW)
  REAL(RealKind), PARAMETER :: CP=0.350
  INTEGER :: iW
  
      DO iW=1,NumberW 
        Power(iW)=CP*MeanRhoGlob(iW)*((MeanWindGlob(iW,1)*SIN(WindFarm(iW)%Alpha)-MeanWindGlob(iW,2)* &
&                 COS(WindFarm(iW)%Alpha))**3)*PI*(WindFarm(iW)%Rad**2)/2
      END DO 
 ! Write Time,WindPower in File
      
        OPEN(99,File='WindPowerOut.dat',status='unknown',action='write',position='append')
        WRITE(99,*) NumberW,TimeAct,Power
        CLOSE(99)
END SUBROUTINE WindPowerCompute


SUBROUTINE TimeProfileCompute(VecOut,PosOut,ibLoc,TimeProfileLabel,Unit1,Unit2,Time)

  TYPE(Vector4Cell_T) :: VecOut(:)
  INTEGER :: PosOut
  INTEGER :: ibLoc
  CHARACTER(8) :: TimeProfileLabel
  INTEGER :: Unit1,Unit2
  REAL(RealKind) :: Time

  ! Local variables
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  INTEGER :: iTime
  INTEGER :: iz
  INTEGER :: Hours,Minutes,Seconds
  CHARACTER(2) :: CHours,CMinutes,CSeconds
  CHARACTER(8) :: TimeOfDay

  Outc=>VecOut(ibLoc)%Vec(PosOut)%c

  iTime=INT(Time)
  Hours=INT(AINT(Time/3600.))
  Minutes=INT(AINT(MOD(iTime,3600)/60.))
  Seconds=(MOD(iTime,60))

  IF (Hours<10) THEN
    WRITE (CHours,"(A1,I1)") '0',Hours 
  ELSE
    WRITE (CHours,"(I2)") Hours 
  ENDIF
  IF (Minutes<10) THEN
    WRITE (CMinutes,"(A1,I1)") '0',Minutes 
  ELSE
    WRITE (CMinutes,"(I2)") Minutes 
  ENDIF
  IF (Seconds<10) THEN
    WRITE (CSeconds,"(A1,I1)") '0',Seconds 
  ELSE
    WRITE (CSeconds,"(I2)") Seconds 
  ENDIF

  TimeOfDay=CHours//':'//CMinutes//':'//CSeconds

  IF (Time==StartTime) THEN
    IF (TimeProfile1.AND.(ix0<TimeProfileIx1).AND.(ix1>=TimeProfileIx1).AND.(iy0<TimeProfileIy1).AND.(iy1>=TimeProfileIy1)) THEN 
      OPEN(UNIT=Unit1,FILE=TRIM(OutputFileName)//'.'//TRIM(TimeProfileLabel)//'.prof1',STATUS='REPLACE')
      WRITE(Unit1,*) '  Time          Height                      ',TRIM(TimeProfileLabel)
    END IF
    IF (TimeProfile2.AND.(ix0<TimeProfileIx2).AND.(ix1>=TimeProfileIx2).AND.(iy0<TimeProfileIy2).AND.(iy1>=TimeProfileIy2)) THEN 
      OPEN(UNIT=Unit2,FILE=TRIM(OutputFileName)//'.'//TRIM(TimeProfileLabel)//'.prof2',STATUS='REPLACE')
      WRITE(Unit2,*) '  Time          Height                      ',TRIM(TimeProfileLabel)
    END IF
  ELSE
    IF (TimeProfile1.AND.(ix0<TimeProfileIx1).AND.(ix1>=TimeProfileIx1).AND.(iy0<TimeProfileIy1).AND.(iy1>=TimeProfileIy1)) THEN 
      OPEN(UNIT=Unit1,FILE=TRIM(OutputFileName)//'.'//TRIM(TimeProfileLabel)//&
      &'.prof1',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
    END IF
    IF (TimeProfile2.AND.(ix0<TimeProfileIx2).AND.(ix1>=TimeProfileIx2).AND.(iy0<TimeProfileIy2).AND.(iy1>=TimeProfileIy2)) THEN 
      OPEN(UNIT=Unit2,FILE=TRIM(OutputFileName)//'.'//TRIM(TimeProfileLabel)//&
      &'.prof2',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
    END IF
  END IF
  DO iz=iz0+1,iz1
    IF (TimeProfile1.AND.(ix0<TimeProfileIx1).AND.(ix1>=TimeProfileIx1).AND.(iy0<TimeProfileIy1).AND.(iy1>=TimeProfileIy1)) THEN
      WRITE(Unit1,*) TimeOfDay,Half*(zP(iz-1)+zP(iz)),OutC(TimeProfileIx1,TimeProfileIy1,iz,1)
    END IF
    IF (TimeProfile2.AND.(ix0<TimeProfileIx2).AND.(ix1>=TimeProfileIx2).AND.(iy0<TimeProfileIy2).AND.(iy1>=TimeProfileIy2)) THEN
      WRITE(Unit2,*) TimeOfDay,Half*(zP(iz-1)+zP(iz)),OutC(TimeProfileIx2,TimeProfileIy2,iz,1)
    END IF
  END DO
  IF (TimeProfile1) CLOSE(Unit1)
  IF (TimeProfile2) CLOSE(Unit2)

END SUBROUTINE TimeProfileCompute


!SUBROUTINE SoilFluxCompute(Time)
!
!  IMPLICIT NONE
!
!  REAL(RealKind) :: Time

!  INTEGER :: i,ix,iy,iz
!  REAL(RealKind) :: VolFull
!  REAL(RealKind) :: SumSHFOcean,SumLHFOcean,SumSHFIsland,SumLHFIsland
!  INTEGER        :: CountSHFOcean,CountLHFOcean,CountSHFIsland,CountLHFIsland
!  INTEGER        :: time_int
!
!  time_int=AINT(Time)
!
!  IF (MOD(time_int,INT(OutputTimeStep))==0) THEN
!    SumSHFOcean    = Zero
!    SumLHFOcean    = Zero
!    CountSHFOcean  = 0
!    CountLHFOcean  = 0
!    SumSHFIsland   = Zero
!    SumLHFIsland   = Zero
!    CountSHFIsland = 0
!    CountLHFIsland = 0
!
!    DO i=1,NumBoundCell
!      ix      = BoundCell(i)%ix
!      iy      = BoundCell(i)%iy
!      iz      = BoundCell(i)%iz
!      VolFull = dx(ix)*dy(iy)*dz(iz) ! full cell volume
!      IF (VolC(ix,iy,iz)==VolFull) THEN ! ocean
!        SumSHFOcean=SumSHFOcean+BoundCell(i)%FluxSens
!        SumLHFOcean=SumLHFOcean+BoundCell(i)%FluxLat
!        CountSHFOcean=CountSHFOcean+1
!        CountLHFOcean=CountLHFOcean+1
!      ELSE
!        SumSHFIsland=SumSHFIsland+BoundCell(i)%FluxSens
!        SumLHFIsland=SumLHFIsland+BoundCell(i)%FluxLat
!        CountSHFIsland=CountSHFIsland+1
!        CountLHFIsland=CountLHFIsland+1
!      END IF
!    END DO
!
!  END IF  ! time is output time
!
!END SUBROUTINE SoilFluxCompute


SUBROUTINE CanopyFluxCompute(Time)

  IMPLICIT NONE

  REAL(RealKind) :: Time

  INTEGER :: i,ix,iy,iz
  INTEGER :: time_int

  REAL(RealKind) :: raddir, raddif, radinf ! direct solar radiation, diffuse solar radiation, long wave radiation
  REAL(RealKind) :: Hs, LEs ! sensible heat flux and latent heat flux on the surface
  REAL(RealKind) :: Tsurf ! surface temperature
  REAL(RealKind) :: albs ! albedo of canopy and soil
  REAL(RealKind), ALLOCATABLE :: SunlitLeafT(:) ! sunlit leaf temperature
  REAL(RealKind), ALLOCATABLE :: ShadedLeafT(:) ! shaded leaf temperature
  REAL(RealKind) :: PRc, GFc ! penetration rate, gap function
  REAL(RealKind), ALLOCATABLE :: Hc(:) ! sensible heat flux from canopy to air
  REAL(RealKind), ALLOCATABLE :: LEc(:) ! latent heat flux from canopy to air

  time_int=AINT(Time)

  ! write(*,*) 'NumBoundCell, NumProcs: ', NumBoundCell, MyID
  ! write(*,*) 'BoundCell size: ', size(BoundCell), allocated(BoundCell)
  IF (MOD(time_int,INT(OutputTimeStep))==0) THEN
    IF (NumBoundCell>0) THEN  ! if there exist boundary cells
      ALLOCATE(SunlitLeafT(1:BoundCell(1)%CanopyCell%NrCanopyLayers))
      ALLOCATE(ShadedLeafT(1:BoundCell(1)%CanopyCell%NrCanopyLayers))
      ALLOCATE(Hc(1:BoundCell(1)%CanopyCell%NrCanopyLayers))
      ALLOCATE(LEc(1:BoundCell(1)%CanopyCell%NrCanopyLayers))

      raddir = 0.0d0
      raddif = 0.0d0
      radinf = 0.0d0
      Hs = 0.0d0
      LEs = 0.0d0
      Tsurf = 0.0d0
      albs = 0.0d0
      SunlitLeafT = 0.0d0
      ShadedLeafT = 0.0d0
      PRc = 0.0d0
      GFc = 0.0d0
      Hc = 0.0d0
      LEc = 0.0d0

      DO i=1,NumBoundCell
        ix = BoundCell(i)%ix
        iy = BoundCell(i)%iy
        iz = BoundCell(i)%iz

        !===== Save the radiation value for the first boundary cell since it is homogeneous.
        write(*,*) 'CanopyFluxCompute 1'
        raddir = raddir + BoundCell(i)%raddirekt
        raddif = raddif + BoundCell(i)%raddiffus
        radinf = radinf + BoundCell(i)%radinfred
        Hs  = Hs + BoundCell(i)%FluxSens
        LEs = LEs + BoundCell(i)%FluxLat
        Tsurf = Tsurf + BoundCell(i)%TeS
        albs = albs + BoundCell(i)%CanopyCell%AlbSoil
        SunlitLeafT = SunlitLeafT + BoundCell(i)%CanopyCell%SunlitLeafT
        ShadedLeafT = ShadedLeafT + BoundCell(i)%CanopyCell%ShadedLeafT
        PRc = PRc + BoundCell(i)%CanopyCell%RadPtrRate
        GFc = GFc + BoundCell(i)%CanopyCell%GapFunc
        Hc = Hc + BoundCell(i)%CanopyCell%Hc
        LEc = LEc + BoundCell(i)%CanopyCell%LEc
      END DO  ! i=1, NumBoundCell

      write(*,*) 'CanopyFluxCompute 2'
      !===== Output radiation and heat fluxes =====!
      IF (time_int == StartTime) then
        OPEN(205,FILE='Radiation.txt',STATUS='REPLACE',ACTION='WRITE')
        WRITE(205, *) NumProcs, '! NumProcs    time', &
          '    raddirekt    raddiffus    radinfred', &
          '    Hs    LEs    TSurf    albs', &
          '    NrCanopyLayers    SunlitLeafT    ShadedLeafT    PRc    GFc', &
          '    Hc    LEc'
        CLOSE(205)
      END IF

      OPEN(205,FILE='Radiation.txt',STATUS='OLD',ACTION='WRITE',POSITION='APPEND')
      WRITE(205,*) Time, raddir/NumBoundCell, raddif/NumBoundCell, radinf/NumBoundCell, &
        Hs/NumBoundCell, LEs/NumBoundCell, TSurf/NumBoundCell, albs/NumBoundCell, &
        BoundCell(1)%CanopyCell%NrCanopyLayers, SunlitLeafT/NumBoundCell, &
        &ShadedLeafT/NumBoundCell, PRc/NumBoundCell, GFc/NumBoundCell, &
        Hc/NumBoundCell, LEc/NumBoundCell
      CLOSE(205)

      write(*,*) 'CanopyFluxCompute 3'
      DEALLOCATE(SunlitLeafT)
      DEALLOCATE(ShadedLeafT)
      DEALLOCATE(Hc)
      DEALLOCATE(LEc)
      write(*,*) 'CanopyFluxCompute 4'
    END IF  ! NumBoundCell>0
  END IF  ! Time is output time

END SUBROUTINE CanopyFluxCompute


SUBROUTINE ProfileVariablesCompute(VecC,Out,Time) 

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)
  REAL(RealKind) :: Time

! Local variables
  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: OutC(:,:,:,:)
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: pre(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoC(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: NCCN(:,:,:,:)
  REAL(RealKind), POINTER :: Nc(:,:,:,:)
  REAL(RealKind), POINTER :: Nr(:,:,:,:)
  REAL(RealKind), POINTER :: Tracer1(:,:,:,:)
  REAL(RealKind), POINTER :: Tracer2(:,:,:,:)
  REAL(RealKind), POINTER :: AbsT(:,:,:,:)

  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,T
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv,RmLoc,CpmlLoc
  REAL(RealKind) :: xPLoc,zPloc,dudz,dvdz,dudy,dvdx,dwdx,dwdy
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB

  CHARACTER(8)   :: Label

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    pre=>PreCell(ibLoc)%c
    Rho=>RhoCell(ibLoc)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoC=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    NCCN=>VecC(ibLoc)%Vec(NvPos)%c
    Nc=>VecC(ibLoc)%Vec(NcPos)%c
    Nr=>VecC(ibLoc)%Vec(NrPos)%c
    Tracer1=>VecC(ibLoc)%Vec(tracer1Pos)%c
    Tracer2=>VecC(ibLoc)%Vec(tracer2Pos)%c
    AbsT=>TAbsCell(ibLoc)%Vec(1)%c
    Outc=>Out(ibLoc)%Vec(LBOUND(NameOutSpecial,1))%c

    ! Relative humidity [%]
    IF (RhoVPos>0.AND.TimeProfileRH.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            T=AbsT(ix,iy,iz,1)
            Outc(ix,iy,iz,1)=RhoV(ix,iy,iz,1)*Rv*T/(SaturVapor(T)+Eps)*100.0d0
          END DO
        END DO
      END DO
      Label='RH'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Absolute humidity [g kg^(-1)]
    IF (RhoVPos>0.AND.TimeProfileRhoV.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=MAX((RhoV(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps))*1.d3,Zero)
          END DO
        END DO
      END DO
      Label='RhoV'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Mass ratio of cloud water [g kg^(-1)]
    IF (RhoCPos>0.AND.TimeProfileRhoC.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=MAX((RhoC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps))*1.d3,Zero)
          END DO
        END DO
      END DO
      Label='RhoC'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Mass ratio of rain water [g kg^(-1)]
    IF (RhoRPos>0.AND.TimeProfileRhoR.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=MAX((RhoR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps))*1.d3,Zero)
          END DO
        END DO
      END DO
      Label='RhoR'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Liquid water content [g m^(-3)]
    IF (RhoCPos>0.AND.RhoRPos>0.AND.TimeProfileLWC.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=MAX((RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1))*1.d3,Zero)
          END DO
        END DO
      END DO
      Label='LWC'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! CCN number concentration [cm^(-3)]
    IF (nvPos>0.AND.TimeProfileNCCN.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=MAX(NCCN(ix,iy,iz,1)*1.d6,Zero)
          END DO
        END DO
      END DO
      Label='NCCN'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Cloud droplets number concentration [cm^(-3)]
    IF (ncPos>0.AND.TimeProfileNc.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=MAX(Nc(ix,iy,iz,1)*1.d6,Zero)
          END DO
        END DO
      END DO
      Label='Nc'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Rain drops number concentration [cm^(-3)]
    IF (nrPos>0.AND.TimeProfileNr.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=MAX(Nr(ix,iy,iz,1)*1.d6,Zero)
          END DO
        END DO
      END DO
      Label='Nr'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Total number concentration [cm^(-3)]
    IF (nvPos>0.AND.ncPos>0.AND.nrPos>0.AND.TimeProfileNTotal.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=MAX((NCCN(ix,iy,iz,1)+Nc(ix,iy,iz,1)+Nr(ix,iy,iz,1))*1.d6,Zero)
          END DO
        END DO
      END DO
      Label='NTotal'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Passive tracer 1 [-]
    IF (tracer1Pos>0.AND.TimeProfileTracer1.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=Tracer1(ix,iy,iz,1)
          END DO
        END DO
      END DO
      Label='Tracer1'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Passive tracer 1 [-]
    IF (tracer2Pos>0.AND.TimeProfileTracer2.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=Tracer2(ix,iy,iz,1)
          END DO
        END DO
      END DO
      Label='Tracer2'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Vertical wind velocity [m s^(-1)]
    IF (TimeProfileW.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=wC(ix,iy,iz,1)
          IF (.NOT.(Outc(ix,iy,iz,1)<9999.9.OR.Outc(ix,iy,iz,1)>-9999.9)) Outc(ix,iy,iz,1)=Zero
        END DO
      END DO
    END DO
      Label='W'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Absolute wind velocity [m s^(-1)]
    IF (TimeProfileAbsVel.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=SQRT( uC(ix,iy,iz,1)**2 + &
                                 vC(ix,iy,iz,1)**2 + & 
                                 wC(ix,iy,iz,1)**2   & 
                               )
          IF (.NOT.(Outc(ix,iy,iz,1)<9999.9.OR.Outc(ix,iy,iz,1)>-9999.9)) Outc(ix,iy,iz,1)=Zero
        END DO
      END DO
    END DO
      Label='AbsVel'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Temperature [K]
    IF (TimeProfileT.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=AbsT(ix,iy,iz,1)
          END DO
        END DO
      END DO
      Label='T'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Temperature [°C]
    IF (TimeProfileTC.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            Outc(ix,iy,iz,1)=AbsT(ix,iy,iz,1)-273.15d0
          END DO
        END DO
      END DO
      Label='TC'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Dry potential temperature [K]
    IF (TimeProfilePotTemp.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            pLoc=pre(ix,iy,iz,1)
            T=AbsT(ix,iy,iz,1)
            Outc(ix,iy,iz,1)=T*(p0/(pLoc+Eps))**(Rd/Cpd)
          END DO
        END DO
      END DO
      Label='PotTemp'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Density potential temperature [K]
    IF (RhoVPos>0.AND.TimeProfileThDens.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            RhoLoc=Rho(ix,iy,iz,1)+Eps
            RhoVLoc=RhoV(ix,iy,iz,1)+Eps
            RhoLLoc=RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1)+Eps
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Three*Eps
            rLLoc=RhoLLoc/(RhoDLoc+Eps)
            rVLoc=RhoVLoc/(RhoDLoc+Eps)
            rt=rLLoc+rVLoc
            pLoc=pre(ix,iy,iz,1)
            T=AbsT(ix,iy,iz,1)
            RmLoc=Rd+rVLoc*Rv
            CpmlLoc=Cpd+rVLoc*Cpv+rLLoc*Cpl
            Outc(ix,iy,iz,1)=T*(p0/(pLoc+Eps))**(RmLoc/CpmlLoc)*(1.+Rv/Rd*rVLoc)/(1.+rt)
          END DO
        END DO
      END DO
      Label='ThDens'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF

    ! Equivalent potential temperature [K]
    IF (RhoVPos>0.AND.TimeProfileThEquiv.AND.MOD(INT(Time),INT(OutputTimeStep))==0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            T=AbsT(ix,iy,iz,1)+Eps
            RhoLoc=Rho(ix,iy,iz,1)+Eps
            RhoVLoc=RhoV(ix,iy,iz,1)+Eps
            RhoLLoc=RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1)+Eps
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Three*Eps
            rLLoc=RhoLLoc/(RhoDLoc+Eps)
            rVLoc=RhoVLoc/(RhoDLoc+Eps)
            rt=rLLoc+rVLoc
            Lv=LatHeat(T)
            Outc(ix,iy,iz,1)=T*((RhoDLoc*Rd*T)/p0)**(-Rd/(Cpd+Cpl*rt)) &
!                             *((RhoVLoc*Rv*T)/(SaturVapor(T)+Eps))**(-rvLoc*Rv/(Cpd+Cpl*rt)) &
                              *EXP(Lv*rvLoc/((Cpd+Cpl*rt)*T+Eps))
          END DO
        END DO
      END DO
      Label='ThEquiv'
      CALL TimeProfileCompute(Out,LBOUND(NameOutSpecial,1),ib,Label,101,102,Time)
    END IF
  END DO

END SUBROUTINE ProfileVariablesCompute


SUBROUTINE RainRateCompute(Vector,Time) 

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: Vector(:)
  REAL(RealKind) :: Time

! Local variables
  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: RR(:,:)
  REAL(RealKind) :: RainRate

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
!   WRITE(*,*) 'ib,ibLoc,nbLoc',ib,ibLoc,nbLoc
    RR=>Vector(ibLoc)%Vec(RhoRPos)%cB
    DO i=1,NumBoundCell
      ix     = BoundCell(i)%ix
      iy     = BoundCell(i)%iy
      RainRate=RR(i,1)
!     IF (ix==50.AND.iy==50) WRITE (*,*) i,'/',NumBoundCell,ix,iy,RainRate,ib
!     WRITE (*,*) i,'/',NumBoundCell,ix,iy,RainRate
    END DO
  END DO

END SUBROUTINE RainRateCompute


SUBROUTINE ShearProdCompute(Vector,Out)

  TYPE(Vector4Cell_T) :: Vector(:),Out(:)
  
  INTEGER :: ibLoc,ib
  INTEGER :: i,ix,iy,iz
  REAL(RealKind) :: S,O
  REAL(RealKind) :: dudz
  REAL(RealKind) :: N2
  REAL(RealKind) :: FM,FH,PhiM
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC
  REAL(RealKind) :: tkeC,rhoC
  REAL(RealKind), POINTER :: uC(:,:,:,:),vC(:,:,:,:),wC(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    Rho=>RhoCell(ibLoc)%c
    th=>Vector(ibLoc)%Vec(thPos)%c
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          VCe  = VolC(ix,iy,iz)
          VB   = VolC(ix,iy,iz-1)
          VT   = VolC(ix,iy,iz+1)
          FB   = FW(ix,iy,iz-1)
          FT   = FW(ix,iy,iz)
          rhoC = Rho(ix,iy,iz,1)
          thC  = th(ix,iy,iz,1)/(rhoC+Eps)
          CALL DeformVortN(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                           vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                           wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                           Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                           VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                           FU(ix-1:ix,iy:iy,iz:iz),              &
                           FV(ix:ix,iy-1:iy,iz:iz),              &
                           FW(ix:ix,iy:iy,iz-1:iz),              &
                           S,O,dudz)
          dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                           ,thC                                      &
                           ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                           ,FB,FT,VB,VCe,VT)
!         Brunt-Vasala-Frequency
          N2 = Grav/(thC+Eps)*dthdz

          IF (PosOutSh > 0) THEN
            Out(ibLoc)%Vec(PosOutSh)%c(ix,iy,iz,1) = S*O
          END IF
          IF (PosOutRi > 0) THEN
            Out(ibLoc)%Vec(PosOutRi)%c(ix,iy,iz,1) = dudz
          END IF
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ShearProdCompute


SUBROUTINE ComputeThetaL(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoC(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: AbsTe(:,:,:,:)
  REAL(RealKind), POINTER :: p(:,:,:,:)
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,Te,pLoc
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: Const1,Const2 

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    AbsTe=>TAbsCell(ibLoc)%Vec(1)%c
    Rho=>RhoCell(ibLoc)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoC=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    p=>PreCell(ibLoc)%c
    Outc=>Out(ibLoc)%Vec(PosOutThetaL)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Te=AbsTe(ix,iy,iz,1)+Eps
          pLoc=p(ix,iy,iz,1)+Eps
          RhoLoc=Rho(ix,iy,iz,1)+Eps
          RhoVLoc=RhoV(ix,iy,iz,1)+Eps
          RhoLLoc=RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1)+Eps
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Three*Eps
          rLLoc=RhoLLoc/(RhoDLoc+Eps)
          rVLoc=RhoVLoc/(RhoDLoc+Eps)
          rt=rLLoc+rVLoc
           Lv=LatHeat(Te)
          Const1=(Rd+rt*Rv)/(Cpd+rt*Cpv)
          Const2=rt*Rv/(Cpd+rt*Cpv)
          Outc(ix,iy,iz,1)=Te*(p0/pLoc)**Const1              &
                            *(One-rLLoc/(Rd/Rv+rt))**Const1 &
                            *(One-rLLoc/rt)**(-Const2)      &
                            *EXP(Lv*rLLoc/((Cpd+rt*Cpv)*Te))
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeThetaL


SUBROUTINE ComputeThEquiv(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoC(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: AbsTe(:,:,:,:)
  REAL(RealKind), POINTER :: p(:,:,:,:)
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,Te,pLoc
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    AbsTe=>TAbsCell(ibLoc)%Vec(1)%c
    Rho=>RhoCell(ibLoc)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoC=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    Outc=>Out(ibLoc)%Vec(PosOutThEquiv)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Te=AbsTe(ix,iy,iz,1)+Eps
          RhoLoc=Rho(ix,iy,iz,1)+Eps
          RhoVLoc=RhoV(ix,iy,iz,1)+Eps
          RhoLLoc=RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1)+Eps
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Three*Eps
          rLLoc=RhoLLoc/(RhoDLoc+Eps)
          rVLoc=RhoVLoc/(RhoDLoc+Eps)
          rt=rLLoc+rVLoc
          Lv=LatHeat(Te)
          Outc(ix,iy,iz,1)=Te*((RhoDLoc*Rd*Te)/p0)**(-Rd/(Cpd+Cpl*rt)) &
                           *((RhoVLoc*Rv*Te)/(SaturVapor(Te)+Eps))**(-rvLoc*Rv/(Cpd+Cpl*rt)) &
                           *EXP(Lv*rvLoc/((Cpd+Cpl*rt)*Te+Eps))
                           WRITE(*,*) 'ThE',ix,iz,Outc(ix,iy,iz,1),rt
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeThEquiv


SUBROUTINE ComputeLen_turb(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: xPLoc,zPloc,dudz,dvdz,dudy,dvdx,dwdx,dwdy
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutLen_turb)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    Rho=>RhoCell(ibLoc)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)= LenKoeff(ibLoc)%c(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeLen_turb


SUBROUTINE ComputeAnelaP(VecG,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecG(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutAnelaP)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)= VecG(ibLoc)%Vec(1)%c(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeAnelaP


SUBROUTINE ComputeDivergence(VecC,VelF,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)
  TYPE(VelocityFace_T), POINTER   :: VelF(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: uF(:,:,:)
  REAL(RealKind), POINTER :: vF(:,:,:)
  REAL(RealKind), POINTER :: wF(:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutDivergence)%c
    uF=>VelF(ibLoc)%uF
    vF=>VelF(ibLoc)%vF
    wF=>VelF(ibLoc)%wF
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          Outc(ix,iy,iz,1)= &
                  (uF(ix,iy,iz)*FU(ix,iy,iz)-uF(ix-1,iy,iz)*FU(ix-1,iy,iz) &
                  +vF(ix,iy,iz)*FV(ix,iy,iz)-vF(ix,iy-1,iz)*FV(ix,iy-1,iz) &
                  +wF(ix,iy,iz)*FW(ix,iy,iz)-wF(ix,iy,iz-1)*FW(ix,iy,iz-1))/Volc(ix,iy,iz)
!         IF (ABS(OutC(ix,iy,iz,1))>1.0d0) THEN
!           WRITE(*,*) 'Div',OutC(ix,iy,iz,1)
!           WRITE(*,*) ix,iy,iz
!           WRITE(*,*) uF(ix,iy,iz),FU(ix,iy,iz)
!           WRITE(*,*) uF(ix-1,iy,iz),FU(ix-1,iy,iz)
!           WRITE(*,*) vF(ix,iy,iz),FV(ix,iy,iz)
!           WRITE(*,*) vF(ix,iy-1,iz),FV(ix,iy-1,iz)
!           WRITE(*,*) wF(ix,iy,iz),FW(ix,iy,iz)
!           WRITE(*,*) wF(ix,iy,iz-1),FW(ix,iy,iz-1)
!           WRITE(*,*) VolC(ix,iy,iz)
!         END IF  
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeDivergence


SUBROUTINE ComputexVorticity(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: xPLoc,zPloc,dudz,dvdz,dudy,dvdx,dwdx,dwdy
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutxVorticity)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    Rho=>RhoCell(ibLoc)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          dwdy = GradCentr(wC(ix,iy-1,iz,1)/(Rho(ix,iy-1,iz,1)+Eps) &
                      ,wC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,wC(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)    &
                       ,FV(ix,iy-1,iz),FV(ix,iy,iz),VolC(ix,iy-1,iz) &
                       ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
          dvdz = GradCentr(vC(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                       ,vC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,vC(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)    &
                       ,FW(ix,iy,iz-1),FW(ix,iy,iz),VolC(ix,iy,iz-1) &
                       ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
          Outc(ix,iy,iz,1)= dwdy-dvdz 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputexVorticity


SUBROUTINE ComputeyVorticity(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: xPLoc,zPloc,dudz,dvdz,dudy,dvdx,dwdx,dwdy
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutyVorticity)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    Rho=>RhoCell(ibLoc)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          dudz = GradCentr(uC(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                       ,uC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,uC(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)    &
                       ,FW(ix,iy,iz-1),FW(ix,iy,iz),VolC(ix,iy,iz-1) &
                       ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
          dwdx = GradCentr(wC(ix-1,iy,iz,1)/(Rho(ix-1,iy,iz,1)+Eps) &
                       ,wC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,wC(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)    &
                       ,FU(ix-1,iy,iz),FU(ix,iy,iz),VolC(ix-1,iy,iz) &
                       ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
          Outc(ix,iy,iz,1)= dudz-dwdx
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeyVorticity


SUBROUTINE ComputeCpuNumb(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: xPLoc,zPloc,dudz,dvdz,dudy,dvdx,dwdx,dwdy
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutCpuNumb)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)= myID 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeCpuNumb


SUBROUTINE ComputeBlkNumb(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutBlkNumb)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)= ib 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeBlkNumb

SUBROUTINE ComputeVolCell(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutVol)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=VolC(ix,iy,iz) 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeVolCell


SUBROUTINE ComputezVorticity(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: xPLoc,zPloc,dudz,dvdz,dudy,dvdx,dwdx,dwdy
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutzVorticity)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    Rho=>RhoCell(ibLoc)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          dudy = GradCentr(uC(ix,iy-1,iz,1)/(Rho(ix,iy-1,iz,1)+Eps) &
                       ,uC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,uC(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)    &
                       ,FV(ix,iy-1,iz),FV(ix,iy,iz),VolC(ix,iy-1,iz) &
                       ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
          dvdx = GradCentr(vC(ix-1,iy,iz,1)/(Rho(ix-1,iy,iz,1)+Eps) &
                       ,vC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,vC(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)    &
                       ,FU(ix-1,iy,iz),FU(ix,iy,iz),VolC(ix-1,iy,iz) &
                       ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
          Outc(ix,iy,iz,1)= dvdx-dudy
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputezVorticity


SUBROUTINE ComputeVorticity(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: xPLoc,zPloc,dudz,dvdz,dudy,dvdx,dwdx,dwdy
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutVorticity)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    Rho=>RhoCell(ibLoc)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          dudy = GradCentr(uC(ix,iy-1,iz,1)/(Rho(ix,iy-1,iz,1)+Eps) &
                       ,uC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,uC(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)    &
                       ,FV(ix,iy-1,iz),FV(ix,iy,iz),VolC(ix,iy-1,iz) &
                       ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
          dwdy = GradCentr(wC(ix,iy-1,iz,1)/(Rho(ix,iy-1,iz,1)+Eps) &
                      ,wC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,wC(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)    &
                       ,FV(ix,iy-1,iz),FV(ix,iy,iz),VolC(ix,iy-1,iz) &
                       ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
          dudz = GradCentr(uC(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                       ,uC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,uC(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)    &
                       ,FW(ix,iy,iz-1),FW(ix,iy,iz),VolC(ix,iy,iz-1) &
                       ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
          dvdx = GradCentr(vC(ix-1,iy,iz,1)/(Rho(ix-1,iy,iz,1)+Eps) &
                       ,vC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,vC(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)    &
                       ,FU(ix-1,iy,iz),FU(ix,iy,iz),VolC(ix-1,iy,iz) &
                       ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
          dwdx = GradCentr(wC(ix-1,iy,iz,1)/(Rho(ix-1,iy,iz,1)+Eps) &
                       ,wC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,wC(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)    &
                       ,FU(ix-1,iy,iz),FU(ix,iy,iz),VolC(ix-1,iy,iz) &
                       ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
          dvdz = GradCentr(vC(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                       ,vC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,vC(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)    &
                       ,FW(ix,iy,iz-1),FW(ix,iy,iz),VolC(ix,iy,iz-1) &
                       ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
          Outc(ix,iy,iz,1)= dwdy-dvdz+dudz-dwdx+dvdx-dudy
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeVorticity


SUBROUTINE ComputeNonDimStreamVel(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: xPLoc,zPloc,dudz
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB,z

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutNonDimStreamVel)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    Rho=>RhoCell(ibLoc)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          dudz = GradCentr(uC(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                        ,uC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,uC(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)    &
                       ,FW(ix,iy,iz-1),FW(ix,iy,iz),VolC(ix,iy,iz-1) &
                       ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
          Outc(ix,iy,iz,1)= (Karm*(iz*(1000/iz1))/0.45d0)*dudz
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeNonDimStreamVel


SUBROUTINE ComputexzStress(VecC,Out) !marcelk

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,Lv
  REAL(RealKind) :: xPLoc,zPloc,dudx,dvdy,dwdz,dudz,dvdz,dudy,dvdx,dwdx,dwdy
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind) :: VT,FB,U,FT,VCE,VB
  REAL(RealKind),DIMENSION(1:6)  :: Sij
  REAL(RealKind) :: S

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutxzStress)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    Rho=>RhoCell(ibLoc)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          dudy = GradCentr(uC(ix,iy-1,iz,1)/(Rho(ix,iy-1,iz,1)+Eps) &
                       ,uC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,uC(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)    &
                       ,FV(ix,iy-1,iz),FV(ix,iy,iz),VolC(ix,iy-1,iz) &
                       ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
          dwdy = GradCentr(wC(ix,iy-1,iz,1)/(Rho(ix,iy-1,iz,1)+Eps) &
                      ,wC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,wC(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)    &
                       ,FV(ix,iy-1,iz),FV(ix,iy,iz),VolC(ix,iy-1,iz) &
                       ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
          dudz = GradCentr(uC(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                       ,uC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,uC(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)    &
                       ,FW(ix,iy,iz-1),FW(ix,iy,iz),VolC(ix,iy,iz-1) &
                       ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
          dvdx = GradCentr(vC(ix-1,iy,iz,1)/(Rho(ix-1,iy,iz,1)+Eps) &
                       ,vC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,vC(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)    &
                       ,FU(ix-1,iy,iz),FU(ix,iy,iz),VolC(ix-1,iy,iz) &
                       ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
          dwdx = GradCentr(wC(ix-1,iy,iz,1)/(Rho(ix-1,iy,iz,1)+Eps) &
                       ,wC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,wC(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)    &
                       ,FU(ix-1,iy,iz),FU(ix,iy,iz),VolC(ix-1,iy,iz) &
                       ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
          dvdz = GradCentr(vC(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                       ,vC(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                       ,vC(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)    &
                       ,FW(ix,iy,iz-1),FW(ix,iy,iz),VolC(ix,iy,iz-1) &
                       ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))

          Sij(1) = Half*(dudx+dudx)
          Sij(2) = Half*(dvdy+dvdy)
          Sij(3) = Half*(dwdz+dwdz)
          Sij(4) = Half*(dudy+dvdx)
          Sij(5) = Half*(dudz+dwdx)
          Sij(6) = Half*(dvdz+dwdy)
          S = SQRT(Two*((Sij(1)*Sij(1))+(Sij(2)*Sij(2))+(Sij(3)*Sij(3))+ &
                      Two*(Sij(4)*Sij(4))+Two*(Sij(5)*Sij(5))+Two*(Sij(6)*Sij(6))))

          Outc(ix,iy,iz,1)= -(Two*(LenKoeff(ibLoc)%c(ix,iy,iz,1)*((dx(ix)+dz(iz))/Two))**Two)* &
                             Sij(5)*ABS(S)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputexzStress


SUBROUTINE ComputeThMinusThProf(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutThMProf)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    ThProf=>ThProfG(ibLoc)%c
    Rho=>RhoCell(ibLoc)%c
    RhoProf=>RhoProfG(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=th(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                         -thProf(ix,iy,iz,1)/(RhoProf(ix,iy,iz,1)+Eps)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeThMinusThProf


SUBROUTINE ComputeRhoMinusRhoProf(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutRhoMProf)%c
    Rho=>RhoCell(ibLoc)%c
    RhoProf=>RhoProfG(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=Rho(ix,iy,iz,1)-RhoProf(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeRhoMinusRhoProf


SUBROUTINE ComputeRhoVMinusRhoVProf(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoVProf(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutRhoVMProf)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoVProf=>RhoVProfG(ibLoc)%c
    Rho=>RhoCell(ibLoc)%c
    RhoProf=>RhoProfG(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=RhoV(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps) &
                         -RhoVProf(ix,iy,iz,1)/(RhoProf(ix,iy,iz,1)+Eps)
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE ComputeRhoVMinusRhoVProf


SUBROUTINE ComputeRhoPlusHeight(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: HeightC(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutRhoPHeight)%c
    Rho=>RhoCell(ibLoc)%c
    HeightC=>HeightG(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=Rho(ix,iy,iz,1)+HeightC(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeRhoPlusHeight


SUBROUTINE ComputeAbsTemp(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: AbsTemp(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind) :: Rm,Cpml,KappaLoc,TLoc
  REAL(RealKind) :: RhoDLoc,RhoLoc,RhoVLoc,RhoLLoc,ThLoc

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    AbsTemp=>TAbsCell(ibLoc)%Vec(1)%c
    th=>ThProfG(ibLoc)%c
    Outc=>Out(ibLoc)%Vec(PosOutAbsTemp)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoLoc=RhoProfG(ibLoc)%c(ix,iy,iz,1)
          IF (RhoVPos>0) THEN
            RhoVLoc=RhoVProfG(ibLoc)%c(ix,iy,iz,1)
          ELSE  
            RhoVLoc=Zero
          END IF  
          IF (RhoCPos>0) THEN
            RhoLLoc=RhoCProfG(ibLoc)%c(ix,iy,iz,1)
          ELSE  
            RhoLLoc=Zero
          END IF  
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
          ThLoc=Th(ix,iy,iz,1)
          Rm=Rd*RhoDLoc+Rv*RhoVLoc
          Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc
          KappaLoc=Rm/Cpml
          TLoc=(Rd*ThLoc/p0**KappaLoc)**(One/(One-KappaLoc))/(Rd*RhoDLoc+Rv*RhoVLoc)     
          Outc(ix,iy,iz,1)=AbsTemp(ix,iy,iz,1) !-TLoc
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeAbsTemp


SUBROUTINE ComputeInterEnergy(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)
  
  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: AbsTemp(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)

  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,TLoc
  REAL(RealKind) :: qd,ql,qv
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    AbsTemp=>TAbsCell(ibLoc)%Vec(1)%c
    Rho=>RhoCell(ibLoc)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    Outc=>Out(ibLoc)%Vec(PosOutIntEne)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc
          TLoc=AbsTemp(ix,iy,iz,1)
          qd=RhoDLoc/(RhoLoc+Eps)
          qv=RhoVLoc/(RhoLoc+Eps)
          ql=RhoLLoc/(RhoLoc+Eps)
!         moist internal energy
          Outc(ix,iy,iz,1)=RhoLoc*((qd*Cvd+qv*Cvv+ql*Cpl)*TLoc+qv*L00)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeInterEnergy


SUBROUTINE ComputeTotalEnergy(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: AbsTemp(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: KinEn(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)

  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,TLoc
  REAL(RealKind) :: Internal,Potential,Kinetic
  REAL(RealKind) :: ql,qd,qv
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    AbsTemp=>TAbsCell(ibLoc)%Vec(1)%c
    KinEn=>KinEnCell(ibLoc)%c
    Rho=>RhoCell(ibLoc)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    Outc=>Out(ibLoc)%Vec(PosOutTotalEne)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc
          TLoc=AbsTemp(ix,iy,iz,1)
          qd=RhoDLoc/(RhoLoc+Eps)
          qv=RhoVLoc/(RhoLoc+Eps)
          ql=RhoLLoc/(RhoLoc+Eps)
          Internal=(qd*Cvd+qv*Cvv+ql*Cpl)*TLoc+qv*L00
          Potential=Grav*Half*(zP(iz-1)+zP(iz))
          Kinetic=KinEn(ix,iy,iz,1)
!         total energy
          Outc(ix,iy,iz,1)=RhoLoc*(Internal+Potential+Kinetic)
        END DO
      END DO
    END DO
  WRITE(70,*) TimeAct,SUM(VolB*Out(ibLoc)%Vec(PosOutTotalEne)%cInt(:,:,:,1))/SUM(VolB) &
                     ,SUM(VolB*VecC(ibLoc)%Vec(RhoPos)%cInt(:,:,:,1)*KinEnCell(ibLoc)%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1))/SUM(VolB)
  WRITE(*,*) 'Energy',TimeAct,SUM(VolB*Out(ibLoc)%Vec(PosOutTotalEne)%cInt(:,:,:,1))/SUM(VolB)
  END DO

END SUBROUTINE ComputeTotalEnergy


SUBROUTINE ComputeEnergyDiff(VecC,Out)

  IMPLICIT NONE 

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: E(:,:,:,:)
  REAL(RealKind), POINTER :: EnergyStart(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)

  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,TLoc
  REAL(RealKind) :: Internal,Potential,Kinetic
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    E=>ECell(ibLoc)%c
    EnergyStart=>EStartCell(ibLoc)%c
    Outc=>Out(ibLoc)%Vec(PosOutEnergyDiff)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=E(ix,iy,iz,1)-EnergyStart(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeEnergyDiff


SUBROUTINE ComputePreDiff(VecC,Out)

  IMPLICIT NONE 

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: p(:,:,:,:)
  REAL(RealKind), POINTER :: PressureStart(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)

 DO ibLoc=1,nbLoc
   ib=LocGlob(ibLoc)
   CALL Set(Floor(ib))
   p=>PreCell(ibLoc)%c
   PressureStart=>PStartCell(ibLoc)%c
   Outc=>Out(ibLoc)%Vec(PosOutPreDiff)%c
   DO iz=iz0+1,iz1
     DO iy=iy0+1,iy1
       DO ix=ix0+1,ix1
         Outc(ix,iy,iz,1)=p(ix,iy,iz,1)-PressureStart(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputePreDiff


SUBROUTINE ComputeThDens(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)


  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoC(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: AbsTe(:,:,:,:)
  REAL(RealKind), POINTER :: pre(:,:,:,:)

  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te
  REAL(RealKind) :: rVLoc,rLLoc,rt,RmLoc,CpmlLoc

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoC=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    Outc=>Out(ibLoc)%Vec(PosOutThDens)%c
    Rho=>RhoCell(ibLoc)%c
    pre=>PreCell(ibLoc)%c
    AbsTe=>TAbsCell(ibLoc)%Vec(1)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoLoc=Rho(ix,iy,iz,1)+Eps
          RhoVLoc=RhoV(ix,iy,iz,1)+Eps
          RhoLLoc=RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1)+Eps
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Three*Eps
          rLLoc=RhoLLoc/(RhoDLoc+Eps)
          rVLoc=RhoVLoc/(RhoDLoc+Eps)
          rt=rLLoc+rVLoc
          pLoc=pre(ix,iy,iz,1)
          Te=AbsTe(ix,iy,iz,1)
          RmLoc=Rd+rVLoc*Rv
          CpmlLoc=Cpd+rVLoc*Cpv+rLLoc*Cpl
          Outc(ix,iy,iz,1)=Te*(p0/(pLoc+Eps))**(RmLoc/CpmlLoc)*(1.+Rv/Rd*rVLoc)/(1.+rt)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeThDens


SUBROUTINE ComputeThetaDry(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: AbsTe(:,:,:,:)
  REAL(RealKind), POINTER :: pre(:,:,:,:)

  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDloc,pLoc,Te

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutThetaDry)%c
    pre=>PreCell(ibLoc)%c
    AbsTe=>TAbsCell(ibLoc)%Vec(1)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          pLoc=pre(ix,iy,iz,1)
          Te=AbsTe(ix,iy,iz,1)
          Outc(ix,iy,iz,1)=Te*(p0/(pLoc+Eps))**(Rd/Cpd)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeThetaDry


SUBROUTINE ComputeThPseudo(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: RhoD(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutThPseudo)%c
    ThProf=>ThProfG(ibLoc)%c
    RhoD=>VecC(ibLoc)%Vec(RhoPos)%c
    RhoProf=>RhoProfG(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=ThProf(ix,iy,iz,1)/(RhoProf(ix,iy,iz,1)+RhoD(ix,iy,iz,1)+Eps)
!                         -ThProf(ix,iy,iz,1)/(RhoProf(ix,iy,iz,1)+Eps) 
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeThPseudo


SUBROUTINE ComputeRelHum(VecC,Out,Time)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)
  REAL(RealKind) :: Time

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)
  REAL(RealKind), POINTER :: AbsT(:,:,:,:)
  REAL(RealKind) ::T
  CHARACTER(8) :: Label
 
  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutRelHum)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    Rho=>RhoCell(ibLoc)%c
    AbsT=>TAbsCell(ibLoc)%Vec(1)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          T=AbsT(ix,iy,iz,1)
          Outc(ix,iy,iz,1)=RhoV(ix,iy,iz,1)*Rv*T/(SaturVapor(T)+Eps)*100.0d0
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeRelHum


SUBROUTINE ComputeRelHumIce(VecC,Out,Time)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)
  REAL(RealKind) :: Time

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)
  REAL(RealKind), POINTER :: AbsT(:,:,:,:)
  REAL(RealKind) ::T
  CHARACTER(8) :: Label
 
  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutRelHumIce)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    Rho=>RhoCell(ibLoc)%c
    AbsT=>TAbsCell(ibLoc)%Vec(1)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          T=AbsT(ix,iy,iz,1)
          Outc(ix,iy,iz,1)=RhoV(ix,iy,iz,1)*Rv*T/(SaturVaporIce(T)+Eps)*100.0d0
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeRelHumIce


SUBROUTINE ComputeRhoV(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
    
  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutRhoVOut)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=RhoV(ix,iy,iz,1)
        END DO
      END DO
    END DO 
  END DO  

END SUBROUTINE ComputeRhoV


SUBROUTINE ComputeRhoVS(VecC,Out)
  
  IMPLICIT NONE
  
  TYPE(Vector4Cell_T) :: VecC(:),Out(:)
  
  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: Pre(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)
  REAL(RealKind) :: RhoD,p,T,PotM
  
  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutRhoVS)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VecC(ibLoc)%Vec(RhoCPos)%c
    ThProf=>DiffKoeff(ibLoc)%c
    Rho=>RhoCell(ibLoc)%c
    Pre=>PreCell(ibLoc)%c
    RhoProf=>RhoProfG(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoD=Rho(ix,iy,iz,1)-RhoV(ix,iy,iz,1)-RhoL(ix,iy,iz,1)+Eps
          p=Pre(ix,iy,iz,1)
          T=AbsTemp(RhoD,RhoV(ix,iy,iz,1),p)
          Outc(ix,iy,iz,1)=SaturVapor(T)/(Rv*T+Eps)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeRhoVS


SUBROUTINE ComputeVPD(VecC,Out) ! Vapor Pressure Deficit

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: Pre(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)
  REAL(RealKind) :: RhoD,p,T,PotM

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutVPD)%c
    th=>VecC(ibLoc)%Vec(thPos)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VecC(ibLoc)%Vec(RhoCPos)%c
    ThProf=>DiffKoeff(ibLoc)%c
    Rho=>RhoCell(ibLoc)%c
    Pre=>PreCell(ibLoc)%c
    RhoProf=>RhoProfG(ibLoc)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoD=Rho(ix,iy,iz,1)-RhoV(ix,iy,iz,1)-RhoL(ix,iy,iz,1)+Three*Eps
          p=Pre(ix,iy,iz,1)
          T=AbsTemp(RhoD,RhoV(ix,iy,iz,1),p)
          Outc(ix,iy,iz,1)=SaturVapor(T)-RhoV(ix,iy,iz,1)*Rv*T
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeVPD


SUBROUTINE ComputeSuperSat(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)
  REAL(RealKind), POINTER :: AbsT(:,:,:,:)
  REAL(RealKind) ::T

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutSuperSat)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    Rho=>RhoCell(ibLoc)%c
    AbsT=>TAbsCell(ibLoc)%Vec(1)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          T=AbsT(ix,iy,iz,1)
          IF (T>Zero) THEN
            Outc(ix,iy,iz,1)=((RhoV(ix,iy,iz,1)*Rv*T/(SaturVapor(T)+eps))-One)*100.0d0  ! [%]
          ELSE
            Outc(ix,iy,iz,1)=Zero
          END IF
          IF (Outc(ix,iy,iz,1)<Zero) Outc(ix,iy,iz,1)=Zero
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeSuperSat


SUBROUTINE ComputeSuperSatIce(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoProf(:,:,:,:)
  REAL(RealKind), POINTER :: AbsT(:,:,:,:)
  REAL(RealKind) ::T

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutSuperSatIce)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    Rho=>RhoCell(ibLoc)%c
    AbsT=>TAbsCell(ibLoc)%Vec(1)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          T=AbsT(ix,iy,iz,1)
          IF (T>Zero) THEN
            Outc(ix,iy,iz,1)=((RhoV(ix,iy,iz,1)*Rv*T/(SaturVaporIce(T)+eps))-One)*100.0d0  ! [%]
          ELSE
            Outc(ix,iy,iz,1)=Zero
          END IF
          IF (Outc(ix,iy,iz,1)<Zero) Outc(ix,iy,iz,1)=Zero
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeSuperSatIce


SUBROUTINE ComputeQLR(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoC(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutQLR)%c
    Rho=>RhoCell(ibLoc)%c
    RhoC=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=(RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1))/(Rho(ix,iy,iz,1)+Eps)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeQLR


SUBROUTINE ComputeQt(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoC(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutQt)%c
    Rho=>RhoCell(ibLoc)%c
    RhoV=>VecC(ibLoc)%Vec(RhoVPos)%c
    RhoC=>VecC(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=(RhoV(ix,iy,iz,1)+RhoC(ix,iy,iz,1)+RhoR(ix,iy,iz,1))/(Rho(ix,iy,iz,1)+Eps)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeQt


SUBROUTINE ComputeLWC(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutLWC)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    RhoL=>VecC(ibLoc)%Vec(RhoCPos)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=RhoR(ix,iy,iz,1)+RhoL(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeLWC


SUBROUTINE ComputeLWP(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind) :: LWP

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutLWP)%c
    RhoR=>VecC(ibLoc)%Vec(RhoRPos)%c
    RhoL=>VecC(ibLoc)%Vec(RhoCPos)%c
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        LWP=Zero
        DO iz=iz0+1,iz1
          LWP=LWP+(zP(iz)-zP(iz-1))*(RhoR(ix,iy,iz,1)+RhoL(ix,iy,iz,1))
          Outc(ix,iy,iz,1)=LWP
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeLWP


SUBROUTINE ComputeIWC(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoI(:,:,:,:)
  REAL(RealKind), POINTER :: RhoS(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutIWC)%c
    RhoI=>VecC(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VecC(ibLoc)%Vec(RhoSPos)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Outc(ix,iy,iz,1)=RhoI(ix,iy,iz,1)+RhoS(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeIWC


SUBROUTINE ComputeIWP(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoI(:,:,:,:)
  REAL(RealKind), POINTER :: RhoS(:,:,:,:)
  REAL(RealKind) :: IWP

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutIWP)%c
    RhoI=>VecC(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VecC(ibLoc)%Vec(RhoSPos)%c
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        IWP=Zero
        DO iz=iz0+1,iz1
          IWP=IWP+(zP(iz)-zP(iz-1))*(RhoI(ix,iy,iz,1)+RhoS(ix,iy,iz,1))
          Outc(ix,iy,iz,1)=IWP
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeIWP


SUBROUTINE ComputeThForcing(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: ThForcing(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutThForcing)%c
    ThForcing=>ForceThCell(ibLoc)%c
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        DO iz=iz0+1,iz1
          Outc(ix,iy,iz,1)=ThForcing(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeThForcing


SUBROUTINE ComputeRhoVForcing(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoVForcing(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutRhoVForcing)%c
    RhoVForcing=>ForceRhoVCell(ibLoc)%c
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        DO iz=iz0+1,iz1
          Outc(ix,iy,iz,1)=RhoVForcing(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeRhoVForcing


SUBROUTINE ComputeHeatRate(VecC,Out)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Outc(:,:,:,:)
  REAL(RealKind), POINTER :: HeatRate(:,:,:,:)

  CALL Set(Floor(1))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Outc=>Out(ibLoc)%Vec(PosOutHeatRate)%c
    HeatRate=>HeatRateCell(ibLoc)%c
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        DO iz=iz0+1,iz1
          Outc(ix,iy,iz,1)=HeatRate(ix,iy,iz,1)
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeHeatRate


SUBROUTINE ComputeParticulateMatter(VecC,Out,PM_StartPos,Max_diameter)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:),Out(:)
  INTEGER :: PM_StartPos
  REAL(RealKind) :: Max_diameter

  INTEGER :: i,ib,ibLoc,ix,iy,iz,is
  REAL(RealKind) :: nC,Rad
  REAL(RealKind) :: MM(0:nFrac+1,nAqua)
  REAL(RealKind) :: mFrac(nAqua)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO i=PM_StartPos,PM_StartPos+NumAeroOut-1
      Out(ibLoc)%Vec(i)%c(:,:,:,:)=0.d0
    END DO
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DO i=1,nFrac
            nC=VecC(ibLoc)%Vec(iNC)%c(ix,iy,iz,i)
            IF (nC.LE.0.d0) CYCLE
            DO is=1,nAqua
              MM(i,is)=VecC(ibLoc)%Vec(is)%c(ix,iy,iz,i)
            END DO
            mFrac=MM(i,:)/nC
            Rad=Radius(mFrac)

            IF (2.d0*Rad.LT.Max_diameter) THEN
              DO is=1,NumAeroOut
                Out(ibLoc)%Vec(PM_StartPos+is-1)%c(ix,iy,iz,1)=Out(ibLoc)%Vec(PM_StartPos+is-1)%c(ix,iy,iz,1)+MM(i,AeroOut(is))
              END DO
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeParticulateMatter


END MODULE SpecialOutput_Mod
