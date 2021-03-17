MODULE Function_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Physics_Mod
  USE Emission_Mod
  USE Deposition_Mod
  USE Koagulation_Mod
  USE CondenAdvect_Mod
  USE Ansatz_Mod 
  USE Adiabatic_Mod
  USE Parameter_Mod
  USE BoundaryCond_Mod
  USE Example_Mod
  USE Activity_Mod
  USE Radiation_Mod
  USE Shadow_Mod
  USE DiffKoeff_Mod
  USE Turbulence_Mod
  USE Operator_Mod
  USE Wall_Mod
  USE BulkMicroPhysics_Mod
  USE TwoMomentBulkMicroPhysics_Mod
  USE ISDACMicroPhysics_Mod
  USE LSCMicroPhysics_Mod
! USE SpectralMicro_Mod
  USE Canopy_Mod
  USE Soil_Mod
  USE WindFarm_Mod
  USE Advection_Mod
  USE Diffusion_Mod
  USE Buoyancy_Mod
  USE VelocityCellFace_Mod
  USE BoundaryCondition_Mod
  USE Diagnostic_Mod
  USE Coriolis_Mod
  USE PressureGrad_Mod
  USE Forcing_Mod
  USE Damping_Mod
  Use Control_Mod

  USE ReadProfile_Mod

  IMPLICIT NONE

  INTERFACE Fcn
    MODULE PROCEDURE FcnMet
  END INTERFACE
  INTERFACE FcnG
    MODULE PROCEDURE FcnMetG
  END INTERFACE
  INTERFACE Update
    MODULE PROCEDURE UpdateScalarCell,UpdateVectorCell
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE ScaleScalarCell,ScaleVectorCell
  END INTERFACE
  INTERFACE Jac
    MODULE PROCEDURE JacMeteo
  END INTERFACE

CONTAINS

SUBROUTINE PressureG(Press,VectorCell)

  TYPE(ScalarCell_T), POINTER :: Press(:)
  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    p=>Press(ibLoc)%c
    Th=>VectorCell(ibLoc)%Vec(ThPos)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
    RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
    IF (RhoCPos>0) THEN 
      RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
    ELSE
      RhoL=>RhoLCell(ibLoc)%c
    END IF
    RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
    CALL AbsTCompute
    IF (ThetaKind=='Energy') THEN
      KinEn=>KinEnCell(ibLoc)%c
      CALL KinEnCompute
    END IF
    CALL PreCompute
  END DO
  CALL ExchangeCell(Press)
END SUBROUTINE PressureG

SUBROUTINE PrepareF(VectorCell,VelocityFace,Time)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)
  REAL(RealKind) :: Time,Time1
  REAL(RealKind) :: RadOutTime
  REAL(RealKind) :: yP0,yP1,yP2,dyP

  INTEGER :: i,j
  INTEGER :: iW

  logical, save :: load = .TRUE.
  logical :: exist
  integer :: iostat
  real(realkind), allocatable, save :: rad(:, :)
  real(realkind) :: glob
  real(realkind) :: rat_dir=0.5d0, lw_rad=262.75433349609375d0
  integer :: nrad

  RadOutTime=OutputTimeStep
  VelocityFaceAct=>VelocityFace
  IF (RadiationValues) THEN
    CALL GetRadiationValues(Time)
  END IF  
  IF (Radiation.AND..NOT.RadiationProfile.AND..NOT.RadiationValues) THEN
    CALL ShadowCompute
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO i=1,NumBoundCell
        ShadowCell(ibLoc)%cB(i,1)=BoundCell(i)%shad
      END DO
    END DO
    IF (Time==StartTime.OR.MOD(INT(Time-StartTime),INT(RadOutTime))==0) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DO i=1, NumBoundCell
          IF (Sphere) THEN
            Time1=Time+86400.0d0/(2.0d0*Pi)*xP(BoundCell(i)%ix)
            CALL SunAngleCompute(Time1,yP(BoundCell(i)%iy))
          ELSE
            CALL SunAngleCompute(Time)
          END IF

          INQUIRE(FILE='rad.prof',EXIST=exist) 
          IF (exist) then
            IF (load) then
              OPEN(UNIT=900,FILE='rad.prof',STATUS='OLD',IOSTAT=iostat)
              READ(900,*) nrad
              ALLOCATE(rad(nrad,2))
              DO j=1,nrad
                READ(900,*) rad(j,:)
              END DO
              CLOSE(900)
              load=.FALSE.
            ELSE
              glob=ProfileEqual(rad,Time)
              BoundCell(i)%raddirekt=glob*rat_dir
              BoundCell(i)%raddiffus=glob-BoundCell(i)%raddirekt
              BoundCell(i)%radinfred=lw_rad
              RaddirCell(ibLoc)%cB(i,1)=BoundCell(i)%raddirekt*BoundCell(i)%shad
              RaddifCell(ibLoc)%cB(i,1)=BoundCell(i)%raddiffus
              RadinfCell(ibLoc)%cB(i,1)=BoundCell(i)%radinfred
            END IF
          ELSE
!           CALL RadiationFuLiouCompute(i)
            BoundCell(i)%raddirekt=400.0d0
            BoundCell(i)%raddiffus=0.0d0
            BoundCell(i)%radinfred=lw_rad
            RaddirCell(ibLoc)%cB(i,1)=BoundCell(i)%raddirekt*BoundCell(i)%shad
            RaddifCell(ibLoc)%cB(i,1)=BoundCell(i)%raddiffus
            RadinfCell(ibLoc)%cB(i,1)=BoundCell(i)%radinfred
          END IF
        END DO
      END DO
    END IF
  ELSE   
    IF (Sphere) THEN
      CALL SunAngleCompute(Time,PhiCor)
    ELSE
      CALL SunAngleCompute(Time)
    END IF  
  END IF
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL SetVelocityFace(ibLoc)
    CALL DomainSet(ib)
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
    RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
    IF (RhoCPos==0) THEN 
      RhoL=>RhoLCell(ibLoc)%c
      IF (iWater>0) THEN
        WaterLiq=>VectorCell(ibLoc)%Vec(iWater)%c
        CALL RhoLCompute
      ELSE
        RhoLCell(ibLoc)%c=0d0
      END IF
    END IF
    IF (ASSOCIATED(KinEnCell)) THEN
      KinEn=>KinEnCell(ibLoc)%c
      CALL KinEnCompute
    END IF  
    IF (ThPos>0) THEN
      Th=>VectorCell(ibLoc)%Vec(ThPos)%c
      T=>TAbsCell(ibLoc)%Vec(1)%c
      CALL AbsTCompute
    END IF  
    IF (ASSOCIATED(ECell)) THEN
      E=>ECell(ibLoc)%c
      CALL ECompute
    END IF  
    IF (PGradient.AND..NOT.(Anelastic.OR.PseudoIn).OR.Parcel) THEN
      p=>PreCell(ibLoc)%c
      CALL PreCompute
      IF (Time==StartTime) THEN
        EnergyStart=>EStartCell(ibLoc)%c
        EnergyStart=E
        PressureStart=>PStartCell(ibLoc)%c
        PressureStart=p
      END IF
      IF (ThetaKind=='PreEn'.OR.ThetaKind=='Exner'.OR.Parcel) THEN
        Sound=>SoundCell(ibLoc)%c
        CALL SoundCompute
      END IF  
    END IF
    IF (DragSurf.OR.DynamicSoil.OR.Canopy.OR.SeaEmiss.OR.FireEmiss) THEN
      CALL DragCoeff(VelocityFace,VectorCell)
    END IF
    IF (uPosL>0) THEN
      uC=>uCell(ibLoc)%c
      vC=>vCell(ibLoc)%c
      wC=>wCell(ibLoc)%c
      CALL VelFaceToCellCompute
      CALL BoundaryVelocityCellCompute(Time)
    END IF
    IF (Aerosol.AND.Depos) THEN
      cVec=>VectorCell(ib)%Vec
      CALL SedimentVelocity(RhoCell(ibLoc),SediCell(ibLoc))
    END IF
    IF (Chemie.AND.DataFile/='') THEN 
      IF (Position('H2O')>0) THEN
        VectorCell(ibLoc)%Vec(Position('H2O'))%c=RhoV/Molmass(Position('H2O')) ! water vapor concentration
      END IF
      IF (Position('O2')>0) THEN
        VectorCell(ibLoc)%Vec(Position('O2'))%c=(Rho-RhoV-RhoL)/Molmass(Position('O2'))*0.20942
      END IF
      IF (Position('N2')>0) THEN
        VectorCell(ibLoc)%Vec(Position('N2'))%c=(Rho-RhoV-RhoL)/Molmass(Position('N2'))*0.78084
      END IF
    END IF
  END DO
  IF (uPosL>0) THEN
    CALL ExchangeCell(uCell) 
    CALL ExchangeCell(vCell) 
    CALL ExchangeCell(wCell) 
  END IF
  IF (ASSOCIATED(KinEnCell)) THEN
    CALL ExchangeCell(KinEnCell) 
  END IF
  IF (ASSOCIATED(PreCell)) THEN
    CALL ExchangeCell(PreCell) 
  END IF
  IF (ASSOCIATED(DivCell)) THEN
    CALL ExchangeCell(DivCell) 
  END IF
  IF (ASSOCIATED(TAbsCell)) THEN
    CALL ExchangeCell(TAbsCell) 
  END IF
  IF (ASSOCIATED(SoundCell)) THEN
    CALL ExchangeCell(SoundCell) 
  END IF
  IF (Parcel.AND.TrajIn) THEN
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL DomainSet(ib)
      CALL MetAmbientComput(VectorCell(ibLoc)%Vec,Time)
    END DO
  END IF
  CALL DiffKoeffSelect(VectorCell)
  IF (ForcingExtern) THEN
    CALL PrepareForcing(VectorCell,Time)
    CALL MeanVerticalProfile(VectorCell)
  END IF  

END SUBROUTINE PrepareF

SUBROUTINE PrepareAbsTemp(VectorCell)
  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
    RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
    RhoEn=>VectorCell(ibLoc)%Vec(EnPos)%c
    p=>PreCell(ibLoc)%c
    Th=>VectorCell(ibLoc)%Vec(ThPos)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    IF (ThPos>0) THEN
      CALL AbsTPreCompute
    END IF  
  END DO  
END SUBROUTINE PrepareAbsTemp
  
SUBROUTINE PrepareFEx(VectorCell,VelocityFace,UVec,Time)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)
  TYPE(Vector4Cell_T), POINTER :: UVec(:)
  REAL(RealKind) :: Time

  INTEGER :: iW

  VelocityFaceAct=>VelocityFace
  IF (RadiationValues) THEN    
    CALL GetRadiationValues(Time)
  END IF 
  IF (Radiation.AND..NOT.RadiationProfile.AND..NOT.RadiationValues) THEN
   CALL SunAngleCompute(Time)
   CALL RadiationFuLiouCompute()
   CALL ShadowCompute
  END IF
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
    RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
    RhoEn=>VectorCell(ibLoc)%Vec(EnPos)%c
    p=>PreCell(ibLoc)%c
    Th=>VectorCell(ibLoc)%Vec(ThPos)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    IF (EnPos>0) THEN
      KinEn=>KinEnCell(ibLoc)%c
      CALL KinEnCompute
      E=>ECell(ibLoc)%c
      CALL ECompute
    END IF  
    IF (ThPos>0) THEN
      WRITE(*,*) 'AbsTPreCompute in PrepareFEX',EnPos
      CALL AbsTPreCompute
    END IF  
    IF (DragSurf.OR.DynamicSoil.OR.Canopy.OR.SeaEmiss.OR.FireEmiss) THEN
      CALL DragCoeff(VelocityFace,VectorCell,UVec)
    END IF
    IF (uPosL>0) THEN
      uC=>uCell(ibLoc)%c
      vC=>vCell(ibLoc)%c
      wC=>wCell(ibLoc)%c
      CALL VelFaceToCellCompute
      CALL BoundaryVelocityCellCompute(Time)
    END IF
  END DO   
  CALL ExchangeCell(TAbsCell)
  CALL ExchangeCell(PreCell)
  IF (EnPos>0) THEN
    CALL ExchangeCell(KinEnCell)
    CALL ExchangeCell(ECell)
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      VectorCell(ibLoc)%Vec(ThPos)%c=PreCell(ibLoc)%c
    END DO  
  END IF  
  IF (uPosL>0) THEN
    CALL ExchangeCell(uCell) 
    CALL ExchangeCell(vCell) 
    CALL ExchangeCell(wCell) 
  END IF
  CALL DiffKoeffSelect(VectorCell,UVec)
  IF (ForcingExtern) THEN
    CALL PrepareForcing(VectorCell,Time)
    CALL MeanVerticalProfile(VectorCell,UVec)
  END IF

END SUBROUTINE PrepareFEx

SUBROUTINE FcnMetSlow(VecUCell,VectorMetCell,VectorChemCell,VelocityFace,RhsMet,RhsChem,ThetaF,PreFacF,SoundFac,Time,dt)

  TYPE(Vector4Cell_T), POINTER :: VecUCell(:)
  TYPE(Vector4Cell_T), POINTER :: VectorMetCell(:)
  TYPE(Vector4Cell_T), POINTER :: VectorChemCell(:)
  TYPE(VelocityFace_T), POINTER :: VelocityFace(:)
  TYPE(Vector4Cell_T), POINTER :: RhsMet(:)
  TYPE(Vector4Cell_T), POINTER :: RhsChem(:)
  TYPE(VectorSFace_T), POINTER :: ThetaF(:)
  TYPE(VelocityFace_T), POINTER :: PreFacF(:)
  TYPE(ScalarCell_T), POINTER :: SoundFac(:)
  REAL(RealKind) :: Time,dt
  CHARACTER*10 :: MethAdvOri

  VelocityFaceAct=>VelocityFace
  IF (ThetaKind=='EnergyBryan') THEN
    ThetaKind='EnergyBryanSlow'
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL DomainSet(ib)
      Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
      RhoV=>VectorMetCell(ibLoc)%Vec(RhoVPos)%c
      RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
      RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
      p=>PreKin(ibLoc)%c
      KinEn=>KinEnCell(ibLoc)%c
      CALL PreCompute
    END DO  
    ThetaKind='EnergyBryan'
    CALL ExchangeCell(PreKin) 
  END IF  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
    IF (ForcingExternTendency) THEN
      wSubs=>Subs(iBloc)%Vec(1)%c
    END IF  
    p=>PreCell(ibLoc)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    DO ic=uPosL,wPosR
      IF (ic>0) THEN
        c=>VecUCell(ibLoc)%Vec(ic)%c
        f=>RhsMet(ibLoc)%Vec(ic)%c
        IF (Advection) THEN
          CALL AdvectionCompute(PhiLim)
        END IF 
      END IF  
    END DO
    IF (Diffusion) THEN
      uCL=>VecUCell(ibLoc)%Vec(uPosL)%c
      vCL=>VecUCell(ibLoc)%Vec(vPosL)%c
      wCL=>VecUCell(ibLoc)%Vec(wPosL)%c
      uCR=>VecUCell(ibLoc)%Vec(uPosR)%c
      vCR=>VecUCell(ibLoc)%Vec(vPosR)%c
      wCR=>VecUCell(ibLoc)%Vec(wPosR)%c
      IF (TkeDis.OR.TkeOme) THEN 
        DH=>DiffKoeff(ibLoc)%c
        DV=>DiffKoeff(ibLoc)%c
      ELSE IF (TkeHVLen) THEN
        DH=>DiffHKoeff(ibLoc)%c
        DV=>DiffVKoeff(ibLoc)%c
      ELSE IF (TkeSGS.OR.NoTke) THEN
        DH=>DiffMomKoeff(ibLoc)%c
        DV=>DiffMomKoeff(ibLoc)%c
      ELSE IF (TkeSmag) THEN
        DH=>DiffMomHKoeff(ibLoc)%c
        DV=>DiffMomVKoeff(ibLoc)%c
      ELSE IF (DynSmag) THEN
        DH=>DiffMomHKoeff(ibLoc)%c
        DV=>DiffMomVKoeff(ibLoc)%c
      ELSE
        DH=>DiffKoeff(ibLoc)%c
        DV=>DiffKoeff(ibLoc)%c
      END IF
      DO ic=uPosL,wPosR
        c=>VecUCell(ibLoc)%Vec(ic)%c
        f=>RhsMet(ibLoc)%Vec(ic)%c
        IF (ic==uPosL.AND.CrossDiff) THEN
          uRhsL=>RhsMet(ibLoc)%Vec(uPosL)%c
          uRhsR=>RhsMet(ibLoc)%Vec(uPosR)%c
          CALL DiffusionComputeU
        ELSE IF (ic==uPosR.AND.CrossDiff) THEN
        ELSE IF (ic==vPosL.AND.CrossDiff) THEN
          vRhsL=>RhsMet(ibLoc)%Vec(vPosL)%c
          vRhsR=>RhsMet(ibLoc)%Vec(vPosR)%c
          CALL DiffusionComputeV
        ELSE IF (ic==vPosR.AND.CrossDiff) THEN
        ELSE IF (ic==wPosL.AND.CrossDiff) THEN
          wRhsL=>RhsMet(ibLoc)%Vec(wPosL)%c
          wRhsR=>RhsMet(ibLoc)%Vec(wPosR)%c
          CALL DiffusionComputeW
        ELSE IF (ic==wPosR.AND.CrossDiff) THEN
        ELSE 
          CALL DiffusionCompute
        END IF
      END DO
!     Diffusion for scalars      
      IF (TkeDis.OR.TkeOme) THEN 
        DH=>DiffKoeff(ibLoc)%c
        DV=>DiffKoeff(ibLoc)%c
      ELSE IF (TkeHVLen) THEN
        DH=>DiffHKoeff(ibLoc)%c
        DV=>DiffVKoeff(ibLoc)%c
      ELSE IF (TkeSGS.OR.NoTke) THEN
        DH=>DiffPotKoeff(ibLoc)%c
        DV=>DiffPotKoeff(ibLoc)%c
      ELSE IF (TkeSmag) THEN
        DH=>DiffPotHKoeff(ibLoc)%c
        DV=>DiffPotVKoeff(ibLoc)%c
      ELSE IF (DynSmag) THEN
        DH=>DiffPotHKoeff(ibLoc)%c
        DV=>DiffPotVKoeff(ibLoc)%c
      ELSE
        DH=>DiffKoeff(ibLoc)%c
        DV=>DiffKoeff(ibLoc)%c
      END IF
      DO ic=1,UBOUND(VectorMetCell(ibLoc)%Vec,1)
        IF (ic<uPosL.OR.ic>wPosR) THEN
          c=>VectorMetCell(ibLoc)%Vec(ic)%c
          f=>RhsMet(ibLoc)%Vec(ic)%c
          IF (TkeDis) THEN 
            DV=DV/PrandtlNumber(ic)
          ELSE IF (TkeHVLen) THEN
            DV=DV/PrandtlNumber(ic)
            DH=DH/PrandtlNumber(ic)
          ELSE
            DV=DV/PrandtlNumber(ic)
          END IF
          CALL DiffusionCompute
          IF (TkeDis) THEN 
            DV=DV*PrandtlNumber(ic)
          ELSE IF (TkeHVLen) THEN
            DV=DV*PrandtlNumber(ic)
            DH=DH*PrandtlNumber(ic)
          ELSE
            DV=DV*PrandtlNumber(ic)
          END IF
        END IF
      END DO
    END IF
    IF (ThetaKind=='PreEn'.OR.ThetaKind=='Exner') THEN 
      Sound=>SoundFac(ibLoc)%c
      Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
      RhoV=>VectorMetCell(ibLoc)%Vec(RhoVPos)%c
      RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
      RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
      Th=>VectorMetCell(ibLoc)%Vec(thPos)%c
      CALL SoundCompute
      WRITE(*,*) 'Sound ',Sound(1,1,1,1),SQRT(Sound(1,1,1,1))
    END IF
    IF (RainSurf) THEN
      CALL RainSurfCompute(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),Time)
    END IF
    IF (IceSurf) THEN
      CALL IceSurfCompute(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),Time)
    END IF
    IF (SnowSurf) THEN
      CALL SnowSurfCompute(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),Time)
    END IF
    ! Surface drag and fluxes of momentum, heat and mass
    IF (DynamicSoil) THEN
      CALL Soil(VectorMetCell(ibLoc),RhsMet(ibLoc),VecUCell(ibLoc),Time=Time)
    ELSE IF (Canopy) THEN
      CALL CanopyCompute(VectorMetCell(ibLoc),RhsMet(ibLoc),Time)
    ELSE IF (DragSurf) THEN
      CALL Drag(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),VecUCell(ibLoc),Time=Time)
    END IF
    IF (Baum) THEN 
      CALL BaumDrag(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),VecUCell(ibLoc)) 
    END IF
  END DO
  IF (ASSOCIATED(SoundFac)) THEN
    CALL ExchangeCell(SoundFac) 
  END IF  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    IF (ThetaKind=='PreEn'.OR.ThetaKind=='Exner') THEN 
      Sound=>SoundFac(ibLoc)%c
      Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
      f=>RhsMet(ibLoc)%Vec(thPos)%c
      c=>VectorMetCell(ibLoc)%Vec(thPos)%c
      IF (PreAdv=='Inner') THEN
        WRITE(*,*) 'PreAdv Inner ',PreAdv
        c=c-Sound
        CALL AdvectionPreCompute(PhiLim)  
        c=c+Sound
      ELSE
        WRITE(*,*) 'PreAdv Outer ',PreAdv
        CALL AdvectionPreCompute(PhiLim)  
      END IF
    END IF  
    IF (PGradient.AND.ThetaKind=='EnergyBryan') THEN
      uRhsL=>RhsMet(ibLoc)%Vec(uPosL)%c
      vRhsL=>RhsMet(ibLoc)%Vec(vPosL)%c
      wRhsL=>RhsMet(ibLoc)%Vec(wPosL)%c
      uRhsR=>RhsMet(ibLoc)%Vec(uPosR)%c
      vRhsR=>RhsMet(ibLoc)%Vec(vPosR)%c
      wRhsR=>RhsMet(ibLoc)%Vec(wPosR)%c
      DUU=>DUUG(ibLoc)%uF
      DUV=>DUUG(ibLoc)%vF
      DUW=>DUUG(ibLoc)%wF
      p=>PreKin(ibLoc)%c
      Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
      CALL PGradComputeC
    END IF
  END DO

  CALL BoundaryFluxCondition(RhsMet,RhsChem)
  CALL ExchangeFlux(RhsMet)
  CALL ExchangeFlux(RhsChem)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    uCL=>VecUCell(ibLoc)%Vec(uPosL)%c
    vCL=>VecUCell(ibLoc)%Vec(vPosL)%c
    wCL=>VecUCell(ibLoc)%Vec(wPosL)%c
    uCR=>VecUCell(ibLoc)%Vec(uPosR)%c
    vCR=>VecUCell(ibLoc)%Vec(vPosR)%c
    wCR=>VecUCell(ibLoc)%Vec(wPosR)%c
    uRhsL=>RhsMet(ibLoc)%Vec(uPosL)%c
    vRhsL=>RhsMet(ibLoc)%Vec(vPosL)%c
    wRhsL=>RhsMet(ibLoc)%Vec(wPosL)%c
    uRhsR=>RhsMet(ibLoc)%Vec(uPosR)%c
    vRhsR=>RhsMet(ibLoc)%Vec(vPosR)%c
    wRhsR=>RhsMet(ibLoc)%Vec(wPosR)%c
    Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
    Th=>VectorMetCell(ibLoc)%Vec(thPos)%c
    RhoV=>VectorMetCell(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
    RhoI=>VectorMetCell(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VectorMetCell(ibLoc)%Vec(RhoSPos)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    p=>PreCell(ibLoc)%c
    DO ic=1,UBOUND(RhsMet(ibLoc)%Vec,1)
      f=>RhsMet(ibLoc)%Vec(ic)%c
      CALL AdvectionScale
    END DO
    IF (Coriolis) THEN
      IF (CoriolisFree) THEN
        Th=>VectorMetCell(ibLoc)%Vec(thPos)%c
        ThProf=>ThProfG(ibLoc)%c
        IF (Sphere) THEN
          CALL CoriolisFreeComputeLR
        ELSE IF (Cylinder) THEN
          CALL CoriolisCylFreeComputeLR
        ELSE
          CALL CoriolisFreeComputeLR
        END IF
      ELSE
        IF (CoriolisProfile) THEN
          uEProf=>uG(ibLoc)%Vec(1)%c
          vEProf=>uG(ibLoc)%Vec(2)%c
          Call CoriolisProfileComputeLR
        ELSE
          ThProf=>ThProfG(ibLoc)%c
          uE=>VecEnv1(ibLoc)%Vec(uPosEnv)%c
          vE=>VecEnv1(ibLoc)%Vec(vPosEnv)%c
          Call CoriolisComputeLR
        END IF  
      END IF
    END IF
    IF (Forcing) THEN
      HeatRate=>HeatRateCell(ibLoc)%c
      ThForcing=>ForceThCell(ibLoc)%c
      RhoVForcing=>ForceRhoVCell(ibLoc)%c
      CALL ForceVelCompute(Time)
      IF (ForcingExtern.AND.ForcingExternTendency) THEN
        CALL ForceScalarCompute(RhsMet(ibLoc)%Vec,Time,TendAdv(ibLoc)%Vec)
      ELSE  
        CALL ForceScalarCompute(RhsMet(ibLoc)%Vec,Time)
      END IF  
    END IF
    IF (ForcingCellPert.AND.(Time==StartTime.OR.(Time>=ForcingCellPertTime.AND.&
    &MOD(INT(Time*1000+0.0001d0),ForcingCellPertTime*1000)==0))) THEN
      CALL CellPert(RhsMet(ibLoc)%Vec,Time)
    END IF
    IF (Subsidence) THEN
      CALL SubsidenceScalar(VectorMetCell(ibLoc)%Vec,RhsMet(ibLoc)%Vec)
    END IF
    IF (Damping) THEN
      CALL ForceDamp(VectorMetCell(ibLoc),RhsMet(ibLoc),Time,VecUCell(ibLoc))
    END IF
    IF (Canopy) THEN
      CALL CanopyDrag(VectorMetCell(ibLoc),RhsMet(ibLoc),VecUCell(ibLoc))
    END IF
    IF (Wind) THEN
      CALL WindDrag(VectorMetCell(ibLoc),RhsMet(ibLoc),VecUCell(ibLoc))
    END IF
    CALL Turbulence(VectorMetCell(ibLoc),RhsMet(ibLoc),VecUCell(ibLoc))
    IF (Chemie) THEN
      TAbs=>TAbsCell(ibLoc)%Vec(1)
      cVec=>VectorMetCell(ibLoc)%Vec
      fVec=>RhsMet(ibLoc)%Vec
      IF (Emiss) THEN
        CALL EmissionCompute(Time,VelocityFace)
        CALL EmissionPointCompute(Time)
      END IF
      IF (EmissStreet) THEN 
        CALL EmissionStreetCompute(Time) 
      END IF
    END IF
    IF (Cloud.AND..NOT.CloudFast) THEN
      IF (MicroScheme=='Bulk') THEN
        CALL BulkMicro(VectorMetCell(ibLoc)%Vec(1:),RhsMet(ibLoc)%Vec(1:))
      ELSE IF (MicroScheme=='Bulk2') THEN
        CALL BulkMicro2(VectorMetCell(ibLoc)%Vec(1:),RhsMet(ibLoc)%Vec(1:),VelocityFace(:),dt,Time)
      ELSE IF (MicroScheme=='ISDAC') THEN
        CALL BulkMicroISDAC(VectorMetCell(ibLoc)%Vec(1:),RhsMet(ibLoc)%Vec(1:),VelocityFace(:),dt,Time)
      ELSE IF (MicroScheme=='LSC') THEN
        CALL LSCMicro(VectorMetCell(ibLoc)%Vec(1:),RhsMet(ibLoc)%Vec(1:))
      END IF    
    END IF    

    DO ic=1,UBOUND(ThetaF(ibLoc)%VecF,1)
      c=>VectorMetCell(ibLoc)%Vec(ic)%c
      cFU=>ThetaF(ibLoc)%VecF(ic)%uF
      cFV=>ThetaF(ibLoc)%VecF(ic)%vF
      cFW=>ThetaF(ibLoc)%VecF(ic)%wF
      cFU=0.0d0
      cFV=0.0d0
      cFW=0.0d0
      IF (ForcingExternTendency) THEN
        wSubs=>Subs(iBloc)%Vec(1)%c
      END IF  
      IF (ic==enPos.AND.ThetaKind=='PreEn') THEN
!       Separation of E (slow) and p (fast)
        WRITE(*,*) 'Face Energy En'
        c=>VectorMetCell(ibLoc)%Vec(EnPos)%c
        p=>VectorMetCell(ibLoc)%Vec(thPos)%c
        c=c+p
        CALL AdvectionFaceCompute(PhiLim)
        c=c-p
      ELSE IF (ic==thPos.AND.(ThetaKind=='PreEn'.OR.ThetaKind=='Exner')) THEN 
        c=>SoundFac(ibLoc)%c
        CALL AdvectionFaceCompute(PhiLim)
      ELSE IF (ic==thPos.AND.(ThetaKind=='Energy'.OR.ThetaKind=='EnergyBryan')) THEN 
        c=>VectorMetCell(ibLoc)%Vec(thPos)%c
        p=>PreCell(ibLoc)%c
        c=c+p
        CALL AdvectionFaceCompute(PhiLim)
        c=c-p
      ELSE IF(ic==RhoRPos) THEN
        RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
        Nrain=>VectorMetCell(ibLoc)%Vec(nrPos)%c
        CALL AdvectionFaceCompute(PhiLim,Species='RhoR')
      ELSE IF (ic==nrPos) THEN
        RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
        Nrain=>VectorMetCell(ibLoc)%Vec(nrPos)%c
        CALL AdvectionFaceCompute(PhiLim,Species='NR')
      ELSE IF (ic==RhoIPos) THEN
        RhoI=>VectorMetCell(ibLoc)%Vec(RhoIPos)%c
        Nice=>VectorMetCell(ibLoc)%Vec(niPos)%c
        CALL AdvectionFaceCompute(PhiLim,Species='RhoI')
      ELSE IF (ic==niPos) THEN
        RhoI=>VectorMetCell(ibLoc)%Vec(RhoIPos)%c
        Nice=>VectorMetCell(ibLoc)%Vec(niPos)%c
        CALL AdvectionFaceCompute(PhiLim,Species='NI')
      ELSE IF (ic==RhoSPos) THEN
        RhoS=>VectorMetCell(ibLoc)%Vec(RhoSPos)%c
        Nsnow=>VectorMetCell(ibLoc)%Vec(nsPos)%c
        CALL AdvectionFaceCompute(PhiLim,Species='RhoS')
      ELSE IF (ic==nsPos) THEN
        RhoS=>VectorMetCell(ibLoc)%Vec(RhoSPos)%c
        Nsnow=>VectorMetCell(ibLoc)%Vec(nsPos)%c
        CALL AdvectionFaceCompute(PhiLim,Species='NS')
      ELSE
        CALL AdvectionFaceCompute(PhiLim)
      END IF
    END DO
!   IF (PrecipRain) THEN
!     CALL DomainSet(ib)
!     CALL SetVelocityFace(ibLoc)
!     Rho=>RhoCell(ibLoc)%c
!     RhoV=>VectorMetCell(ibLoc)%Vec(RhoVPos)%c
!     RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
!     RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
!     p=>PreCell(ibLoc)%c
!     IF (ASSOCIATED(TAbsCell)) THEN
!       T=>TAbsCell(ibLoc)%Vec(1)%c
!     END IF  
!     RhoRhsMet=>RhsMet(ibLoc)%Vec(RhoPos)%c  
!     ThRhsMet=>RhsMet(ibLoc)%Vec(ThPos)%c
!     DO ic=1,UBOUND(VectorMetCell(ibLoc)%Vec,1)
!       c=>VectorMetCell(ibLoc)%Vec(ic)%c
!       IF (ic==RhoPos.OR.ic==ThPos) THEN
!         CALL AdvectionQFallCompute(PhiLim1,FallF)
!       END IF
!     END DO    
!   END IF !Precip
  END DO
  CALL ExchangeScalarFace(ThetaF) 
END SUBROUTINE FcnMetSlow

SUBROUTINE FcnMetFastU(VectorCell,VelF,PreFacF,RhsF,Time)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(VelocityFace_T), POINTER :: VelF(:)
  TYPE(VelocityFace_T), POINTER :: PreFacF(:)
  TYPE(VelocityFace_T), POINTER :: RhsF(:)
  REAL(RealKind) :: Time

  INTEGER :: iz

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    th=>VectorCell(ibLoc)%Vec(thPos)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    Rho=>VectorCell(ibLoc)%Vec(rhoPos)%c
    IF (PGradient) THEN
      IF (RhoCPos==0) THEN 
        RhoL=>RhoLCell(ibLoc)%c
      END IF  
      p=>PreCell(ibLoc)%c
      RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
      RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
      RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
      RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
      RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
      CALL PreCompute
      DummyVec(ibLoc)%Vec(1)%c=>PreCell(ibloc)%c
      DummyVec(ibLoc)%Vec(2)%c=>VectorCell(ibloc)%Vec(RhoPos)%c
    END IF
  END DO
  IF (PGradient) THEN
    CALL ExchangeCell(DummyVec)
  END IF  
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    uFRhs=>RhsF(ibLoc)%uF
    vFRhs=>RhsF(ibLoc)%vF
    wFRhs=>RhsF(ibLoc)%wF
    IF (PGradient) THEN
      p=>PreCell(ibLoc)%c
      CALL PGradComputeF
    END IF
    IF (Buoyancy) THEN
      Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
      RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
      RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
      RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
      RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
      RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
      Th=>VectorCell(ibLoc)%Vec(thPos)%c
      ThProf=>ThProfG(ibLoc)%c
      p=>PreCell(ibLoc)%c
      Call AuftriebComputeF
    END IF
    IF (BaumFast) THEN 
      CALL BaumDragF(VectorCell(ibLoc),VelF(ibLoc),RhsF(ibLoc))
    END IF
  END DO
  
END SUBROUTINE FcnMetFastU

SUBROUTINE FcnMetFastU1(VectorCell,VelF,PreFacF,RhsF,Time)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(VelocityFace_T), POINTER :: VelF(:)
  TYPE(VelocityFace_T), POINTER :: PreFacF(:)
  TYPE(VelocityFace_T), POINTER :: RhsF(:)
  REAL(RealKind) :: Time

  INTEGER :: iz

  ! Only with nonlinear pressure
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    p=>PreCell(ibLoc)%c
    th=>VectorCell(ibLoc)%Vec(thPos)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    Rho=>VectorCell(ibLoc)%Vec(rhoPos)%c
    RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
    RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
    CALL AbsTPreCompute
  END DO
  IF (PGradient) THEN
    CALL ExchangeCell(PreCell)
  END IF  
  CALL Assign(DummyCell,VectorCell,RhoPos)
  CALL ExchangeCell(DummyCell)
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    uFRhs=>RhsF(ibLoc)%uF
    vFRhs=>RhsF(ibLoc)%vF
    wFRhs=>RhsF(ibLoc)%wF
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    IF (PGradient) THEN
      CALL PGradComputeF
    END IF
    IF (Buoyancy) THEN
      Call AuftriebComputeF
    END IF
  END DO
  
END SUBROUTINE FcnMetFastU1

SUBROUTINE FcnMetFastScalar(VectorCell,Velocityface,ThetaF,SoundFac,Rhs,Time)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(VelocityFace_T), POINTER :: VelocityFace(:)
  TYPE(VectorSFace_T), POINTER :: ThetaF(:)
  TYPE(ScalarCell_T), POINTER :: SoundFac(:)
  TYPE(Vector4Cell_T), POINTER :: Rhs(:)
  REAL(RealKind) :: Time

  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    DO ic=1,UBOUND(ThetaF(ibLoc)%VecF,1)
      f=>Rhs(ibLoc)%Vec(ic)%c
      cFU=>ThetaF(ibLoc)%VecF(ic)%uF
      cFV=>ThetaF(ibLoc)%VecF(ic)%vF
      cFW=>ThetaF(ibLoc)%VecF(ic)%wF
      IF (ForcingExternTendency) THEN
        wSubs=>Subs(iBloc)%Vec(1)%c
      END IF  
      IF (Advection) THEN
        IF(ic==RhoRPos) THEN
          RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
          Nrain=>VectorCell(ibLoc)%Vec(nrPos)%c
          CALL DivScalarCompute(Species='RhoR')
        ELSE IF (ic==nrPos) THEN
          RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
          Nrain=>VectorCell(ibLoc)%Vec(nrPos)%c
          CALL DivScalarCompute(Species='NR')
        ELSE IF (ic==RhoIPos) THEN
          RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
          Nice=>VectorCell(ibLoc)%Vec(niPos)%c
          CALL DivScalarCompute(Species='RhoI')
        ELSE IF (ic==niPos) THEN
          RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
          Nice=>VectorCell(ibLoc)%Vec(niPos)%c
          CALL DivScalarCompute(Species='NI')
        ELSE IF (ic==RhoSPos) THEN
          RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
          Nsnow=>VectorCell(ibLoc)%Vec(nsPos)%c
          CALL DivScalarCompute(Species='RhoS')
        ELSE IF (ic==nsPos) THEN
          RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
          Nsnow=>VectorCell(ibLoc)%Vec(nsPos)%c
          CALL DivScalarCompute(Species='NS')
        ELSE
          CALL DivScalarCompute
        END IF
!       IF (ic==RhoRPos.OR.ic==nrPos) THEN ! precipitation
!         vFall=FallVelocity(ic)
!         RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
!         CALL DivScalarCompute(FallF)
!       ELSE
!         CALL DivScalarCompute
!       END IF
      END IF ! Advection
    END DO
!   IF (ThetaKind=='PreEn'.AND.EnPos>0) THEN
!     p=>VectorCell(ibLoc)%Vec(ThPos)%c
!     Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
!     f=>Rhs(ibLoc)%Vec(EnPos)%c
!     CALL DivEnCompute
!   END IF  
  END DO

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    DO ic=1,UBOUND(Rhs(ibLoc)%Vec,1)
      f=>Rhs(ibLoc)%Vec(ic)%c
      CALL AdvectionScale
    END DO
    IF (ThetaKind=='PreEn'.OR.ThetaKind=='Exner') THEN
      Sound=>SoundFac(ibLoc)%c
      thRhs=>Rhs(ibLoc)%Vec(thPos)%c
      Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    END IF  
    IF (Cloud.AND.CloudFast) THEN
      Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
      IF (MicroScheme=='Bulk') THEN
        CALL BulkMicro(VectorCell(ibLoc)%Vec(1:),Rhs(ibLoc)%Vec(1:))
      ELSE IF (MicroScheme=='Bulk2') THEN
        CALL BulkMicro2(VectorCell(ibLoc)%Vec(1:),Rhs(ibLoc)%Vec(1:),Velocityface(:),dt,Time)
      ELSE IF (MicroScheme=='ISDAC') THEN
        CALL BulkMicroISDAC(VectorCell(ibLoc)%Vec(1:),Rhs(ibLoc)%Vec(1:),Velocityface(:),dt,Time)
      ELSE IF (MicroScheme=='LSC') THEN
        CALL LSCMicro(VectorCell(ibLoc)%Vec(1:),Rhs(ibLoc)%Vec(1:))
      END IF 
    END IF
  END DO

END SUBROUTINE FcnMetFastScalar

SUBROUTINE FcnMetG(VectorCell,Velocityface,Rhs,Time,dt)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(VelocityFace_T), POINTER :: VelocityFace(:)
  TYPE(Vector4Cell_T), POINTER :: VectorCellG(:)
  TYPE(Vector4Cell_T), POINTER :: Rhs(:)
  REAL(RealKind) :: Time,dt

  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    Rho=>RhoCell(ibLoc)%c
    IF (Anelastic) THEN
      c=>RhoCell(ibLoc)%c
    ELSE IF (PseudoIn) THEN  
      c=>ThProfG(ibLoc)%c
    END IF  
    f=>Rhs(ibLoc)%Vec(1)%c
    IF (ForcingExternTendency) THEN
      wSubs=>Subs(iBloc)%Vec(1)%c
    END IF  
    CALL AdvectionCompute(PhiLim)
  END DO
  CALL BoundaryFluxCondition(Rhs,Rhs)
  CALL ExchangeFlux(Rhs)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    DO ic=1,SIZE(Rhs(ibLoc)%Vec)
      f=>Rhs(ibLoc)%Vec(ic)%c
      CALL AdvectionScale
    END DO
  END DO
END SUBROUTINE FcnMetG

SUBROUTINE FcnMet(VectorMetCell,VectorChemCell,Velocityface,RhsMet,RhsChem,Time,dt,VectorMetCellG)

  TYPE(Vector4Cell_T), POINTER :: VectorMetCell(:)
  TYPE(Vector4Cell_T), POINTER :: VectorChemCell(:)
  TYPE(VelocityFace_T), POINTER :: VelocityFace(:)
  TYPE(Vector4Cell_T), POINTER :: RhsMet(:)
  TYPE(Vector4Cell_T), POINTER :: RhsChem(:)
  REAL(RealKind) :: Time,dt
  TYPE(Vector4Cell_T), OPTIONAL, POINTER :: VectorMetCellG(:)

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: wLoc
  INTEGER :: l
  REAL(RealKind) :: Temp

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
    IF (Diffusion.AND.CrossDiff.AND.TkePos>0) THEN
      uRhsL=>RhsMet(ibLoc)%Vec(uPosL)%c
      vRhsL=>RhsMet(ibLoc)%Vec(vPosL)%c
      wRhsL=>RhsMet(ibLoc)%Vec(wPosL)%c
      uRhsR=>RhsMet(ibLoc)%Vec(uPosR)%c
      vRhsR=>RhsMet(ibLoc)%Vec(vPosR)%c
      wRhsR=>RhsMet(ibLoc)%Vec(wPosR)%c
      tke=>VectorMetCell(ibLoc)%Vec(tkePos)%c
      CALL TurbTke  
    END IF
  END DO  

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
    IF (ForcingExternTendency) THEN
      wSubs=>Subs(iBloc)%Vec(1)%c
    END IF  
    IF (ASSOCIATED(TAbsCell)) THEN
      T=>TAbsCell(ibLoc)%Vec(1)%c
    END IF  
    IF (Advection) THEN
      DO ic=1,UBOUND(VectorMetCell(ibLoc)%Vec,1)
        c=>VectorMetCell(ibLoc)%Vec(ic)%c
        RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
        f=>RhsMet(ibLoc)%Vec(ic)%c
        vFall=FallVelocity(ic)
        BC=>BCMetVec(ic)
        IF (ic==RhoPos) THEN
          c=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
          CALL AdvectionCompute(PhiLim)
        ELSE IF (ic==ThPos.AND.ThetaKind=='Energy') THEN
          p=>PreCell(ibLoc)%c
          c=c+p
          CALL AdvectionCompute(PhiLim)
          c=c-p
        ELSE IF (ic==enPos.AND.ThetaKind=='PreEn') THEN
          p=>VectorMetCell(ibLoc)%Vec(thPos)%c
          c=c+p
          CALL AdvectionCompute(PhiLim)
          c=c-p
        ELSE IF (ic==thPos.AND.ThetaKind=='PreEn') THEN 
          p=>SoundCell(ibLoc)%c
          IF (PreAdv=='Inner') THEN
            c=c-p
            CALL AdvectionPreCompute(PhiLim) 
            c=c+p
          ELSE
            CALL AdvectionPreCompute(PhiLim)
          END IF
        ELSE IF(ic==RhoRPos) THEN
          RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
          Nrain=>VectorMetCell(ibLoc)%Vec(nrPos)%c
          CALL AdvectionCompute(PhiLim,Species='RhoR')
        ELSE IF (ic==nrPos) THEN
          RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
          Nrain=>VectorMetCell(ibLoc)%Vec(nrPos)%c
          CALL AdvectionCompute(PhiLim,Species='NR')
        ELSE IF (ic==RhoIPos) THEN
          RhoI=>VectorMetCell(ibLoc)%Vec(RhoIPos)%c
          Nice=>VectorMetCell(ibLoc)%Vec(niPos)%c
          CALL AdvectionCompute(PhiLim,Species='RhoI')
        ELSE IF (ic==niPos) THEN
          RhoI=>VectorMetCell(ibLoc)%Vec(RhoIPos)%c
          Nice=>VectorMetCell(ibLoc)%Vec(niPos)%c
          CALL AdvectionCompute(PhiLim,Species='NI')
        ELSE IF (ic==RhoSPos) THEN
          RhoS=>VectorMetCell(ibLoc)%Vec(RhoSPos)%c
          Nsnow=>VectorMetCell(ibLoc)%Vec(nsPos)%c
          CALL AdvectionCompute(PhiLim,Species='RhoS')
        ELSE IF (ic==nsPos) THEN
          RhoS=>VectorMetCell(ibLoc)%Vec(RhoSPos)%c
          Nsnow=>VectorMetCell(ibLoc)%Vec(nsPos)%c
          CALL AdvectionCompute(PhiLim,Species='NS')
        ELSE IF (ic<=nAqua.AND.Depos) THEN
          VelSedi=>SediCell(ibLoc)%c
          CALL AdvectionCompute(PhiLim,Species='AE') ! Aerosol
        ELSE
          CALL AdvectionCompute(PhiLim)
        END IF
      END DO     
      DO ic=1,UBOUND(VectorChemCell(ibLoc)%Vec,1)
        c=>VectorChemCell(ibLoc)%Vec(ic)%c
        f=>RhsChem(ibLoc)%Vec(ic)%c
        BC=>BCChemVec(ic)
        vFall=0.0d0
        CALL AdvectionCompute(PhiLim)
      END DO  
    END IF
!   IF (ThetaKind=='PreEn'.AND.EnPos>0) THEN
!     p=>VectorMetCell(ibLoc)%Vec(thPos)%c
!     Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
!     f=>RhsMet(ibLoc)%Vec(EnPos)%c
!     CALL DivEnCompute
!   END IF
    IF (ThetaKind=='PreEn'.AND.PreAdv=='Inner'.AND.ThPos>0) THEN
      Weight=>SoundCell(ibLoc)%c
      Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
      f=>RhsMet(ibLoc)%Vec(thPos)%c
      CALL DivWeightCompute 
    END IF  
    IF (Diffusion) THEN
!     Diffusion for momentum      
      IF (uPosL>0) THEN
        IF (TkeDis.OR.TkeOme) THEN 
          DH=>DiffKoeff(ibLoc)%c
          DV=>DiffKoeff(ibLoc)%c
        ELSE IF (TkeHVLen) THEN
          DH=>DiffHKoeff(ibLoc)%c
          DV=>DiffVKoeff(ibLoc)%c
        ELSE IF (TkeSGS.OR.NoTke) THEN
          DH=>DiffMomKoeff(ibLoc)%c
          DV=>DiffMomKoeff(ibLoc)%c
        ELSE IF (TkeSmag) THEN
          DH=>DiffMomHKoeff(ibLoc)%c
          DV=>DiffMomVKoeff(ibLoc)%c
        ELSE IF (DynSmag) THEN
          DH=>DiffMomHKoeff(ibLoc)%c
          DV=>DiffMomVKoeff(ibLoc)%c
        ELSE
          DH=>DiffKoeff(ibLoc)%c
          DV=>DiffKoeff(ibLoc)%c
        END IF
        DO ic=uPosL,wPosR
          c=>VectorMetCell(ibLoc)%Vec(ic)%c
          f=>RhsMet(ibLoc)%Vec(ic)%c
          IF (ic==uPosL.AND.CrossDiff) THEN
            uCL=>VectorMetCell(ibLoc)%Vec(uPosL)%c
            uCR=>VectorMetCell(ibLoc)%Vec(uPosR)%c
            vCL=>VectorMetCell(ibLoc)%Vec(vPosL)%c
            vCR=>VectorMetCell(ibLoc)%Vec(vPosR)%c
            wCL=>VectorMetCell(ibLoc)%Vec(wPosL)%c
            wCR=>VectorMetCell(ibLoc)%Vec(wPosR)%c
            uRhsL=>RhsMet(ibLoc)%Vec(uPosL)%c
            uRhsR=>RhsMet(ibLoc)%Vec(uPosR)%c
            CALL DiffusionComputeU
          ELSE IF (ic==uPosR.AND.CrossDiff) THEN
          ELSE IF (ic==vPosL.AND.CrossDiff) THEN
            vRhsL=>RhsMet(ibLoc)%Vec(vPosL)%c
            vRhsR=>RhsMet(ibLoc)%Vec(vPosR)%c
            uCL=>VectorMetCell(ibLoc)%Vec(uPosL)%c
            uCR=>VectorMetCell(ibLoc)%Vec(uPosR)%c
            vCL=>VectorMetCell(ibLoc)%Vec(vPosL)%c
            vCR=>VectorMetCell(ibLoc)%Vec(vPosR)%c
            wCL=>VectorMetCell(ibLoc)%Vec(wPosL)%c
            wCR=>VectorMetCell(ibLoc)%Vec(wPosR)%c
            CALL DiffusionComputeV
          ELSE IF (ic==vPosR.AND.CrossDiff) THEN
          ELSE IF (ic==wPosL.AND.CrossDiff) THEN
            wRhsL=>RhsMet(ibLoc)%Vec(wPosL)%c
            wRhsR=>RhsMet(ibLoc)%Vec(wPosR)%c
            uCL=>VectorMetCell(ibLoc)%Vec(uPosL)%c
            uCR=>VectorMetCell(ibLoc)%Vec(uPosR)%c
            vCL=>VectorMetCell(ibLoc)%Vec(vPosL)%c
            vCR=>VectorMetCell(ibLoc)%Vec(vPosR)%c
            wCL=>VectorMetCell(ibLoc)%Vec(wPosL)%c
            wCR=>VectorMetCell(ibLoc)%Vec(wPosR)%c
            CALL DiffusionComputeW
          ELSE IF (ic==wPosR.AND.CrossDiff) THEN
          ELSE 
            CALL DiffusionCompute
          END IF
        END DO
      END IF
!     Diffusion for scalars      
      IF (TkeDis.OR.TkeOme) THEN 
        DH=>DiffKoeff(ibLoc)%c
        DV=>DiffKoeff(ibLoc)%c
      ELSE IF (TkeHVLen) THEN
        DH=>DiffHKoeff(ibLoc)%c
        DV=>DiffVKoeff(ibLoc)%c
      ELSE IF (TkeSGS.OR.NoTke) THEN
        DH=>DiffPotKoeff(ibLoc)%c
        DV=>DiffPotKoeff(ibLoc)%c
      ELSE IF (TkeSmag) THEN
        DH=>DiffPotHKoeff(ibLoc)%c
        DV=>DiffPotVKoeff(ibLoc)%c
      ELSE IF (DynSmag) THEN
        DH=>DiffPotHKoeff(ibLoc)%c
        DV=>DiffPotVKoeff(ibLoc)%c
      ELSE
        DH=>DiffKoeff(ibLoc)%c
        DV=>DiffKoeff(ibLoc)%c
      END IF
      DO ic=1,UBOUND(VectorMetCell(ibLoc)%Vec,1)
        IF (ic<uPosL.OR.ic>wPosR) THEN
          c=>VectorMetCell(ibLoc)%Vec(ic)%c
          f=>RhsMet(ibLoc)%Vec(ic)%c
          IF (TkeDis) THEN 
            DV=DV/PrandtlNumber(ic)
          ELSE IF (TkeHVLen) THEN
            DV=DV/PrandtlNumber(ic)
            DH=DH/PrandtlNumber(ic)
          ELSE
            DV=DV/PrandtlNumber(ic)
          END IF
          CALL DiffusionCompute
          IF (TkeDis) THEN 
            DV=DV*PrandtlNumber(ic)
          ELSE IF (TkeHVLen) THEN
            DV=DV*PrandtlNumber(ic)
            DH=DH*PrandtlNumber(ic)
          ELSE
            DV=DV*PrandtlNumber(ic)
          END IF
        END IF
      END DO
      DO ic=1,UBOUND(VectorChemCell(ibLoc)%Vec,1)
        c=>VectorChemCell(ibLoc)%Vec(ic)%c
        f=>RhsChem(ibLoc)%Vec(ic)%c
        IF (TkeDis) THEN 
          DV=DV
        ELSE IF (TkeHVLen) THEN
          DV=DV
          DH=DH
        ELSE
          DV=DV
        END IF
        CALL DiffusionCompute
        IF (TkeDis) THEN 
          DV=DV
        ELSE IF (TkeHVLen) THEN
          DV=DV
          DH=DH
        ELSE
          DV=DV
        END IF
      END DO  
    END IF
    IF (RainSurf) THEN
      CALL RainSurfCompute(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),Time)
    END IF
    IF (IceSurf) THEN
      CALL IceSurfCompute(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),Time)
    END IF
    IF (SnowSurf) THEN
      CALL SnowSurfCompute(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),Time)
    END IF
    IF (RadiationProfile) THEN
      CALL RadiationCompute(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc))
    END IF
    ! Surface drag and fluxes of momentum, heat and mass
    IF (DynamicSoil) THEN
      CALL Soil(VectorMetCell(ibLoc),RhsMet(ibLoc),Time=Time)
    ELSE IF (Canopy) THEN
      CALL CanopyCompute(VectorMetCell(ibLoc),RhsMet(ibLoc),Time)
    ELSE IF (DragSurf) THEN
      CALL Drag(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc),Time=Time)
    END IF
    IF (Baum) THEN 
      CALL BaumDrag(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc)) 
    END IF
    IF (Traffic) THEN 
      CALL TrafficDrag(VectorMetCell(ibLoc),Velocityface(ibLoc),RhsMet(ibLoc)) 
    END IF

    IF (PGradient) THEN
      uRhsL=>RhsMet(ibLoc)%Vec(uPosL)%c
      vRhsL=>RhsMet(ibLoc)%Vec(vPosL)%c
      wRhsL=>RhsMet(ibLoc)%Vec(wPosL)%c
      uRhsR=>RhsMet(ibLoc)%Vec(uPosR)%c
      vRhsR=>RhsMet(ibLoc)%Vec(vPosR)%c
      wRhsR=>RhsMet(ibLoc)%Vec(wPosR)%c
      DUU=>DUUG(ibLoc)%uF
      DUV=>DUUG(ibLoc)%vF
      DUW=>DUUG(ibLoc)%wF
      IF (Anelastic.OR.PseudoIn) THEN
        p=>VectorMetCellG(ibLoc)%Vec(1)%c
      ELSE
        p=>PreCell(ibLoc)%c
      END IF
      th=>VectorMetCell(ibLoc)%Vec(thPos)%c
      Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
      CALL PGradComputeC
    END IF
  END DO
  
! IF (PrecipRain) THEN
!   DO ibLoc=1,nbLoc
!     ib=LocGlob(ibLoc)
!     CALL DomainSet(ib)
!     CALL SetVelocityFace(ibLoc)
!     Rho=>RhoCell(ibLoc)%c
!     RhoV=>VectorMetCell(ibLoc)%Vec(RhoVPos)%c
!     RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
!     RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
!     p=>PreCell(ibLoc)%c
!     IF (ASSOCIATED(TAbsCell)) THEN
!       T=>TAbsCell(ibLoc)%Vec(1)%c
!     END IF  
!     RhoRhsMet=>RhsMet(ibLoc)%Vec(RhoPos)%c  
!     ThRhsMet=>RhsMet(ibLoc)%Vec(ThPos)%c
!     DO ic=1,UBOUND(VectorMetCell(ibLoc)%Vec,1)
!       c=>VectorMetCell(ibLoc)%Vec(ic)%c
!       IF (ic==RhoPos.OR.ic==ThPos) THEN
!         CALL AdvectionQFallCompute(PhiLim1,FallF)
!       END IF
!     END DO    
!   END DO
! END IF !Precip

  CALL BoundaryFluxCondition(RhsMet,RhsChem)
  CALL ExchangeFlux(RhsMet)
  CALL ExchangeFlux(RhsChem)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    uCL=>VectorMetCell(ibLoc)%Vec(uPosL)%c
    vCL=>VectorMetCell(ibLoc)%Vec(vPosL)%c
    wCL=>VectorMetCell(ibLoc)%Vec(wPosL)%c
    uCR=>VectorMetCell(ibLoc)%Vec(uPosR)%c
    vCR=>VectorMetCell(ibLoc)%Vec(vPosR)%c
    wCR=>VectorMetCell(ibLoc)%Vec(wPosR)%c
    uRhsL=>RhsMet(ibLoc)%Vec(uPosL)%c
    vRhsL=>RhsMet(ibLoc)%Vec(vPosL)%c
    wRhsL=>RhsMet(ibLoc)%Vec(wPosL)%c
    uRhsR=>RhsMet(ibLoc)%Vec(uPosR)%c
    vRhsR=>RhsMet(ibLoc)%Vec(vPosR)%c
    wRhsR=>RhsMet(ibLoc)%Vec(wPosR)%c
    Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
    RhoV=>VectorMetCell(ibLoc)%Vec(RhoVPos)%c
    RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
    RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
    Th=>VectorMetCell(ibLoc)%Vec(ThPos)%c
    p=>PreCell(ibLoc)%c
    T=>TAbsCell(ibLoc)%Vec(1)%c
    DO ic=1,UBOUND(VectorMetCell(ibLoc)%Vec,1)
      f=>RhsMet(ibLoc)%Vec(ic)%c
      CALL AdvectionScale
    END DO
    DO ic=1,UBOUND(VectorChemCell(ibLoc)%Vec,1)
      f=>RhsChem(ibLoc)%Vec(ic)%c
      CALL AdvectionScale
    END DO
    IF (ThetaKind=='PreEn'.AND.PreAdv=='Outer') THEN
      thRhs=>RhsMet(ibLoc)%Vec(thPos)%c
      Sound=>SoundCell(ibLoc)%c
      CALL DivPreCompute
    END IF
    CALL Turbulence(VectorMetCell(ibLoc),RhsMet(ibLoc))
    IF (Forcing) THEN
      HeatRate=>HeatRateCell(ibLoc)%c
      ThForcing=>ForceThCell(ibLoc)%c
      RhoVForcing=>ForceRhoVCell(ibLoc)%c
      CALL ForceVelCompute(Time)
      IF (ForcingExternTendency) THEN
        CALL ForceScalarCompute(RhsMet(ibLoc)%Vec,Time,TendAdv(ibLoc)%Vec)
      ELSE  
        CALL ForceScalarCompute(RhsMet(ibLoc)%Vec,Time)
      END IF  
    END IF
    IF (ForcingCellPert.AND.(Time==StartTime.OR.(Time>=ForcingCellPertTime&
    & .AND.MOD(INT(Time*1000+0.0001d0),ForcingCellPertTime*1000)==0))) THEN
      CALL CellPert(RhsMet(ibLoc)%Vec,Time)
    END IF
    IF (Subsidence) THEN
      CALL SubsidenceScalar(VectorMetCell(ibLoc)%Vec,RhsMet(ibLoc)%Vec)
    END IF
    IF (Damping) THEN
      CALL ForceDamp(VectorMetCell(ibLoc),RhsMet(ibLoc),Time)
    END IF
    IF (Cloud) THEN
      IF (MicroScheme=='Bulk') THEN
        CALL BulkMicro(VectorMetCell(ibLoc)%Vec(1:),RhsMet(ibLoc)%Vec(1:))
      ELSE IF (MicroScheme=='Bulk2') THEN
        CALL BulkMicro2(VectorMetCell(ibLoc)%Vec(1:),RhsMet(ibLoc)%Vec(1:),Velocityface(:),dt,Time)
      ELSE IF (MicroScheme=='ISDAC') THEN
        CALL BulkMicroISDAC(VectorMetCell(ibLoc)%Vec(1:),RhsMet(ibLoc)%Vec(1:),Velocityface(:),dt,Time)
      ELSE IF (MicroScheme=='LSC') THEN
        CALL LSCMicro(VectorMetCell(ibLoc)%Vec(1:),RhsMet(ibLoc)%Vec(1:))
      END IF 
    END IF
    IF (Buoyancy) THEN
      fL=>RhsMet(ibLoc)%Vec(wPosL)%c
      fR=>RhsMet(ibLoc)%Vec(wPosR)%c
      Call AuftriebDryComputeLR
    END IF
    IF (Parcel) THEN
      cVec=>VectorMetCell(ibLoc)%Vec
      fVec=>RhsMet(ibLoc)%Vec
      IF (TrajIn) THEN 
        CALL RelaxationParcel(cVec,fVec,Time)
      ELSE  
        CALL Adiabatic(cVec,fVec,Time)
      END IF  
      IF (DilutParcel) THEN
        CALL DilutionParcel(cVec,fVec,Time)
      END IF
    END IF
    IF (Shallow.AND.Height) THEN
      HeightC=>HeightG(ibLoc)%c
      CALL HeightComputeLR 
    END IF
    IF (Coriolis) THEN 
      IF (CoriolisFree) THEN
        IF (Sphere) THEN
          CALL CoriolisFreeComputeLR 
        ELSE IF (Cylinder) THEN
          IF (Anelastic) THEN
            Th=>VectorMetCell(ibLoc)%Vec(ThPos)%c
            ThProf=>ThProfG(ibLoc)%c
          END IF   
          CALL CoriolisCylFreeComputeLR 
        ELSE
          CALL CoriolisFreeComputeLR 
        END IF
      ELSE
        uE=>VecEnv1(ibLoc)%Vec(uPosEnv)%c
        vE=>VecEnv1(ibLoc)%Vec(vPosEnv)%c
        Call CoriolisComputeLR
      END IF
    END IF
    ! Add canopy drag terms to the momentum equations
    IF (Canopy) THEN
      CALL CanopyDrag(VectorMetCell(ibLoc),RhsMet(ibLoc))
    END IF
    IF (Wind) THEN
      CALL WindDrag(VectorMetCell(ibLoc),RhsMet(ibLoc))
    END IF
    IF (Curvature) THEN 
      IF (Sphere) THEN
        Call CurvatureComputeLR
      ELSE IF (Cylinder) THEN
        Call CurvatureCylComputeLR
      END IF
    END IF
    IF (Centrifugal) THEN 
      IF (Anelastic) THEN
        ThProf=>ThProfG(ibLoc)%c
      END IF  
      IF (Cylinder) THEN
        Call CentrifugalCylComputeLR
      ELSE
        Call CentrifugalComputeLR
      END IF
    END IF
    IF (Chemie) THEN
      IF (ChemieGas) Then
        TAbs=>TAbsCell(ibLoc)%Vec(1)
        IF (RhoCPos==0) THEN
          RhoLC=>RhoLCell(ibLoc)
        END IF
        cVec=>VectorChemCell(ibLoc)%Vec
        fVec=>RhsChem(ibLoc)%Vec
        CALL GasChemie(cVec,fVec,TAbs)
      END IF
      IF (Aerosol) THEN
        TAbs=>TAbsCell(ibLoc)%Vec(1)
        IF (RhoCPos==0) THEN
          RhoLC=>RhoLCell(ibLoc)
        END IF
        cVec=>VectorChemCell(ibLoc)%Vec
        fVec=>RhsChem(ibLoc)%Vec
        AVec=>Act(ibLoc)%Vec
        IF (ChemieAqua) THEN
          CALL ComputeActivity(cVec,AVec,TAbsCell(ibLoc)%Vec(1))  
          CALL AquaChemie(cVec,fVec,TAbs,AVec)
          CALL DissChemie(cVec,fVec,TAbs,AVec)
          CALL SolidChemie(cVec,fVec,TAbs,AVec)
        END IF
        IF (Condens) THEN
          CALL Condensation(cVec,fVec,HenryTrans,AVec)
        END IF
        IF (Koagul) THEN
          CALL Koagulation(cVec,fVec)
        END IF
      END IF
      IF (Emiss) THEN
        TAbs=>TAbsCell(ibLoc)%Vec(1)
        IF (RhoCPos==0) THEN
          RhoLC=>RhoLCell(ibLoc)
        END IF
        cVec=>VectorChemCell(ibLoc)%Vec
        fVec=>RhsChem(ibLoc)%Vec
        CALL EmissionCompute(Time,VelocityFace)
        CALL EmissionPointCompute(Time)
      END IF
      IF (EmissStreet) THEN 
        cVec=>VectorChemCell(ibLoc)%Vec
        fVec=>RhsChem(ibLoc)%Vec
        CALL EmissionStreetCompute(Time) 
      END IF
      IF (Depos) THEN
        IF (RhoCPos==0) THEN
          RhoLC=>RhoLCell(ibLoc)
        END IF
        IF (Aerosol) SediVelC=>SediCell(ibLoc)
        cVec=>VectorChemCell(ibLoc)%Vec
        fVec=>RhsChem(ibLoc)%Vec
        CALL Deposition(VelocityFace)
      END IF
    END IF
  END DO  

END SUBROUTINE FcnMet

SUBROUTINE JacMeteo(VectorMetCell,VectorChemCell,VelocityFace,JacMatrix,Time)

  TYPE(Vector4Cell_T) :: VectorMetCell(:)
  TYPE(Vector4Cell_T) :: VectorChemCell(:)
  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)
  TYPE(JacSpMatrix4_T) :: JacMatrix(:)
  REAL(RealKind) :: Time

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    Rho=>VectorMetCell(ibLoc)%Vec(RhoPos)%c
    IF (ForcingExternTendency) THEN
      wSubs=>Subs(iBloc)%Vec(1)%c
    END IF  
    AT=>JacMatrix(ibLoc)%JacTMom
    AFall=>JacMatrix(ibLoc)%JacFall
    AFallRhoL=>JacMatrix(ibLoc)%JacFallRhoL
    CALL SetVelocityFace(ibLoc)
    IF (JacAdvection.AND.Advection) THEN
      IF (JacPartial) THEN
        FU=>DUJac(ibLoc)%uF
        FV=>DUJac(ibLoc)%vF
        FW=>DUJac(ibLoc)%wF
      END IF  
      IF (Advection) THEN
        RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
        RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
        IF (RhoCPos==0.AND.iWater.NE.0) THEN
          RhoL=>RhoLCell(ibLoc)%c
        END IF
        IF (RhoRPos>0) THEN
          CALL JacAdvectionCompute(FallF) 
        ELSE
          CALL JacAdvectionCompute
        END IF
      END IF
      IF (JacPartial) THEN
        FU=>Floor(ib)%WeiFU
        FV=>Floor(ib)%WeiFV
        FW=>Floor(ib)%WeiFW
      END IF  
    END IF 
    IF (JacDiffusion.AND.Diffusion) THEN
      IF (JacPartial) THEN
        FU=>DUJac(ibLoc)%uF
        FV=>DUJac(ibLoc)%vF
        FW=>DUJac(ibLoc)%wF
      END IF  
      IF (TkeHVLen) THEN
        DH=>DiffHKoeff(ibLoc)%c
        DV=>DiffVKoeff(ibLoc)%c
        CALL JacDiffusionCompute
      ELSE IF (TkeSGS.OR.NoTke) THEN
        JacMatrix(ibLoc)%JacTPot=JacMatrix(ibLoc)%JacTMom
        DH=>DiffMomKoeff(ibLoc)%c
        DV=>DiffMomKoeff(ibLoc)%c
        CALL JacDiffusionCompute
        DH=>DiffPotKoeff(ibLoc)%c
        DV=>DiffPotKoeff(ibLoc)%c
        AT=>JacMatrix(ibLoc)%JacTPot
        CALL JacDiffusionCompute
      ELSE IF (TkeSmag) THEN
        JacMatrix(ibLoc)%JacTPot=JacMatrix(ibLoc)%JacTMom
        DH=>DiffMomHKoeff(ibLoc)%c
        DV=>DiffMomVKoeff(ibLoc)%c
        CALL JacDiffusionCompute
        DH=>DiffPotHKoeff(ibLoc)%c
        DV=>DiffPotVKoeff(ibLoc)%c
        AT=>JacMatrix(ibLoc)%JacTPot
        CALL JacDiffusionCompute
      ELSE IF (DynSmag) THEN
        JacMatrix(ibLoc)%JacTPot=JacMatrix(ibLoc)%JacTMom
        DH=>DiffMomHKoeff(ibLoc)%c
        DV=>DiffMomVKoeff(ibLoc)%c
        CALL JacDiffusionCompute
        DH=>DiffPotHKoeff(ibLoc)%c
        DV=>DiffPotVKoeff(ibLoc)%c
        AT=>JacMatrix(ibLoc)%JacTPot
        CALL JacDiffusionCompute
      ELSE
        DH=>DiffKoeff(ibLoc)%c
        DV=>DiffKoeff(ibLoc)%c
        CALL JacDiffusionCompute
      END IF
      IF (JacPartial) THEN
        FU=>Floor(ib)%WeiFU
        FV=>Floor(ib)%WeiFV
        FW=>Floor(ib)%WeiFW
      END IF  
    END IF
    ATMom=>JacMatrix(ibLoc)%JacTMom
    ATPot=>JacMatrix(ibLoc)%JacTPot
    CALL JacScaleCompute

    AS=>JacMatrix(ibLoc)%JacSMetLU%Mat
    DiagP=>JacMatrix(ibLoc)%JacSMetLU%Struct%DiagPtr(:)
    Permu=>JacMatrix(ibLoc)%JacSMetLU%Struct%Permu(:)
    IF (DynamicSoil) THEN
      CALL JacSoil(VectorMetCell(ibLoc),AS)
    ELSE IF (Canopy) THEN
      CALL JacCanopyCompute(JacMatrix(ibLoc)%JacSMetLU, VectorMetCell(ibLoc),AS)
    ELSE IF (DragSurf) THEN
      CALL JacDrag(VectorMetCell(ibLoc),Velocityface(ibLoc),AS)
    END IF
    IF (RainSurf) THEN
      CALL JacRainSurf(VectorMetCell(ibLoc),AS)
    END IF
    IF (IceSurf) THEN
      CALL JacIceSurf(VectorMetCell(ibLoc),AS)
    END IF
    IF (SnowSurf) THEN
      CALL JacSnowSurf(VectorMetCell(ibLoc),AS)
    END IF
    IF (Baum) THEN 
      CALL JacBaumDrag(VectorMetCell(ibLoc),AS) 
    END IF
    IF (Traffic) THEN 
      CALL JacTrafficDrag(VectorMetCell(ibLoc),AS) 
    END IF
    CALL JacTurbulence(VectorMetCell(ibLoc),AS)
    IF (Damping) THEN
      CALL JacForceDamp(JacMatrix(ibLoc)%JacSMetLU)
    END IF
    IF (Parcel) THEN
      cVec=>VectorMetCell(ibLoc)%Vec
      dfVec=>JacMatrix(ibLoc)%JacSMetLU%Mat
      IF (TrajIn) THEN 
!       CALL RelaxationParcelJac(cVec,fVec,Time,ib)
      ELSE  
        CALL AdiabaticJac(cVec,dfVec,Time)
      END IF  
      IF (DilutParcel) THEN
        CALL DilutionParcelJac(cVec,dfVec,JacMatrix(ibLoc)%JacSMetLU,Time)
      END IF
    END IF
    IF (Cloud) THEN
      th=>VectorMetCell(ibLoc)%Vec(thPos)%c
      RhoV=>VectorMetCell(ibLoc)%Vec(RhoVPos)%c
      RhoV=>VectorMetCell(ibLoc)%Vec(RhoVPos)%c
      RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
      RhoR=>VectorMetCell(ibLoc)%Vec(RhoRPos)%c
      IF (RhoCPos>0) THEN
        RhoL=>VectorMetCell(ibLoc)%Vec(RhoCPos)%c
      ELSE
        RhoL=>RhoLCell(ibLoc)%c
      END IF  
      IF (JacMicro) THEN
        IF (MicroScheme=='Bulk') THEN
          CALL BulkMicroJac(VectorMetCell(ibLoc)%Vec(1:),JacMatrix(ibLoc)%JacSMetLU%Mat)
        ELSE IF (MicroScheme=='Bulk2') THEN
          !CALL BulkMicroJac2(VectorMetCell(ibLoc)%Vec(1:),JacMatrix(ibLoc)%JacSLU%Mat,ibLoc,Velocityface(:),dt)
        ELSE IF (MicroScheme=='ISDAC') THEN
          !CALL BulkMicroJacISDAC(VectorMetCell(ibLoc)%Vec(1:),JacMatrix(ibLoc)%JacSLU%Mat,ibLoc,Velocityface(:),dt)
        ELSE IF (MicroScheme=='LSC') THEN
          CALL LSCMicroJac(VectorMetCell(ibLoc)%Vec(1:),JacMatrix(ibLoc)%JacSMetLU%Mat)
        END IF
      END IF
      IF (Buoyancy) THEN
        wCL=>VectorMetCell(ibLoc)%Vec(wPosL)%c
        wCR=>VectorMetCell(ibLoc)%Vec(wPosR)%c
        CALL JacAuftriebDryComputeLR
      END IF
    ELSE
      IF (Buoyancy) THEN
        wCL=>VectorMetCell(ibLoc)%Vec(wPosL)%c
        wCR=>VectorMetCell(ibLoc)%Vec(wPosR)%c
        IF (Anelastic) THEN
          ThProf=>ThProfG(ibLoc)%c
          Rho=>RhoCell(ibLoc)%c
        END IF
!       CALL JacAuftriebDryComputeLR 
      END IF
    END IF
    IF (Coriolis) THEN
      Th=>VectorMetCell(ibLoc)%Vec(ThPos)%c
      ThProf=>ThProfG(ibLoc)%c
      IF (Sphere) THEN
        CALL JacCoriolisComputeLR
      ELSE IF (Cylinder) THEN
        CALL JacCoriolisCylCompute
      ELSE
        CALL JacCoriolisComputeLR
      END IF
    END IF
    ! Add canopy drag terms to the momentum equations
    IF (Canopy .AND. (.NOT. DynamicSoil) .AND. (.NOT. DragSurf)) THEN
      CALL JacCanopyDrag(VectorMetCell(ibLoc),AS)
    END IF
    IF (Curvature) THEN
      uCL=>VectorMetCell(ibLoc)%Vec(uPosL)%c
      uCR=>VectorMetCell(ibLoc)%Vec(uPosL)%c
      vCL=>VectorMetCell(ibLoc)%Vec(vPosL)%c
      vCR=>VectorMetCell(ibLoc)%Vec(vPosL)%c
      IF (Sphere) THEN
        CALL JacCurvatureCompute 
      ELSE IF (Cylinder) THEN
        CALL JacCurvatureCylCompute
      END IF
    END IF
    IF (Centrifugal) THEN 
      uPosJac=uPosLJac
      vPosJac=vPosLJac
      wPosJac=wPosLJac
      IF (Anelastic) THEN
        Th=>VectorMetCell(ibLoc)%Vec(ThPos)%c
        ThProf=>ThProfG(ibLoc)%c
      END IF
      IF (Cylinder) THEN
        Call JacCentrifugalCylComputeLR
      ELSE
        Call JacCentrifugalComputeLR
      END IF
      uPosJac=uPosRJac
      vPosJac=vPosRJac
      wPosJac=wPosRJac
      IF (Cylinder) THEN
        Call JacCentrifugalCylComputeLR
      ELSE
        Call JacCentrifugalComputeLR
      END IF
    END IF
    IF (Chemie) THEN
      IF (ChemieGas) THEN
        cVec=>VectorMetCell(ibLoc)%Vec
        dfVec=>JacMatrix(ibLoc)%JacSChemLU%Mat
        TAbs=>TAbsCell(ibLoc)%Vec(1)
        IF (RhoCPos==0) THEN
          RhoLC=>RhoLCell(ibLoc)
        END IF
        CALL GasChemieJac(cVec,dfVec,TAbs)
      END IF
      IF (Aerosol) THEN
        cVec=>VectorMetCell(ibLoc)%Vec
        dfVec=>JacMatrix(ibLoc)%JacSChemLU%Mat
        TAbs=>TAbsCell(ibLoc)%Vec(1)
        IF (RhoCPos==0) THEN
          RhoLC=>RhoLCell(ibLoc)
        END IF
        AVec=>Act(ibLoc)%Vec
        IF (ChemieAqua) THEN
          CALL ComputeActivity(cVec,AVec,TAbsCell(ibLoc)%Vec(1))  
          CALL AquaChemieJac(cVec,dfVec,TAbs,AVec)
          CALL DissChemieJac(cVec,dfVec,TAbs,AVec)
          CALL SolidChemieJac(cVec,dfVec,TAbs,AVec)
        END IF
        IF (Condens) THEN
          CALL CondensationJac(cVec,dfVec,HenryTrans,AVec)
        END IF  
      END IF
    END IF
    IF (Depos) THEN
      IF (Aerosol) SediVelC=>SediCell(ibLoc)
      CALL DepositionJac(JacMatrix(ibLoc)%JacSChemLU,VelocityFace,cVec,dfVec)
    END IF
    CALL Gefa(JacMatrix(ibLoc)%JacSChemLU)
    CALL Gefa(JacMatrix(ibLoc)%JacSMetLU)
  END DO
END SUBROUTINE JacMeteo

SUBROUTINE Limiter(VectorT,VectorT1)

  TYPE (Vector4Cell_T) :: VectorT(:)
  TYPE (Vector4Cell_T) :: VectorT1(:)

  DO ibLoc=1,nbLoc
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL LimiterCompute(VectorT(ibLoc)%Vec(1:nAqua),VectorT1(ibLoc)%Vec(1:nAqua))
  END DO
  CALL ExchangeCell(VectorT)

END SUBROUTINE Limiter

SUBROUTINE JacScaleCompute

  INTEGER :: i
  INTEGER :: ix,iy,iz,iDiag
  INTEGER :: nxP2,nyP2,nzP2

  nxP2=ix1-ix0+2
  nyP2=iy1-iy0+2
  nzP2=iz1-iz0+2

  IF (TypeW(1:1)/='o') THEN
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        i=Index(ix,iy,iz)
        ATMom%Val(i,4)=ATMom%Val(i,4)-ATMom%Val(i,3)
      END DO
    END DO
  END IF
  IF (TypeE(1:1)/='o') THEN
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        i=Index(ix+1,iy,iz)
        ATMom%Val(i,4)=ATMom%Val(i,4)-ATMom%Val(i,5)
      END DO
    END DO
  END IF
  IF (TypeS(1:1)/='o') THEN
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        i=Index(ix,iy,iz)
        ATMom%Val(i,4)=ATMom%Val(i,4)-ATMom%Val(i,2)
      END DO
    END DO
  END IF
  IF (TypeN(1:1)/='o') THEN
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        i=Index(ix,iy+1,iz)
        ATMom%Val(i,4)=ATMom%Val(i,4)-ATMom%Val(i,6)
      END DO
    END DO
  END IF

  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        i=Index(ix,iy,iz)
        ATMom%Val(i,4)=ATMom%Val(i,4)-ATMom%Val(i,1)
      END DO
    END DO
  END IF
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        i=Index(ix,iy,iz+1)
        ATMom%Val(i,4)=ATMom%Val(i,4)-ATMom%Val(i,7)
        AFall%Val(i,1)=AFall%Val(i,1)-AFall%Val(i,2)
      END DO
    END DO
    IF (ASSOCIATED(AFallRhoL%Val)) THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          i=Index(ix,iy,iz+1)
          AFallRhoL%Val(i,1)=AFallRhoL%Val(i,1)-AFallRhol%Val(i,2)
        END DO
      END DO
    END IF
  END IF

  ATMom%Val(1:ATMom%n,ATMom%Diag)=One-(beta0*dt)*ATMom%Val(1:ATMom%n,ATMom%Diag)
  DO iDiag=1,ATMom%NumDiag
    IF (ATMom%DiagPtr(iDiag)/=0) THEN
      ATMom%Val(1:ATMom%n,iDiag)=-(beta0*dt)*ATMom%Val(1:ATMom%n,iDiag) 
    END IF
  END DO
  DO iDiag=1,AFall%NumDiag
    AFall%Val(1:AFall%n,iDiag)=-(beta0*dt)*AFall%Val(1:AFall%n,iDiag) 
  END DO
  IF (ASSOCIATED(AFallRhoL%Val)) THEN
    DO iDiag=1,AFallRhoL%NumDiag
      AFallRhoL%Val(1:AFall%n,iDiag)=-(beta0*dt)*AFallRhoL%Val(1:AFall%n,iDiag) 
    END DO
  END IF
  IF (TkeSGS.OR.NoTke.OR.TkeSmag.OR.DynSmag) THEN
    IF (TypeW(1:1)/='o') THEN
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          i=Index(ix,iy,iz)
          ATPot%Val(i,4)=ATPot%Val(i,4)-ATPot%Val(i,3)
        END DO
      END DO
    END IF
    IF (TypeE(1:1)/='o') THEN
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          i=Index(ix+1,iy,iz)
          ATPot%Val(i,4)=ATPot%Val(i,4)-ATPot%Val(i,5)
        END DO
      END DO
    END IF
    IF (TypeS(1:1)/='o') THEN
      iy=iy0
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          i=Index(ix,iy,iz)
          ATPot%Val(i,4)=ATPot%Val(i,4)-ATPot%Val(i,2)
        END DO
      END DO
    END IF
    IF (TypeN(1:1)/='o') THEN
      iy=iy1
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          i=Index(ix,iy+1,iz)
          ATPot%Val(i,4)=ATPot%Val(i,4)-ATPot%Val(i,6)
        END DO
      END DO
    END IF
    IF (TypeB(1:1)/='o') THEN
      iz=iz0
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          i=Index(ix,iy,iz)
          ATPot%Val(i,4)=ATPot%Val(i,4)-ATPot%Val(i,1)
        END DO
      END DO
    END IF
    IF (TypeT(1:1)/='o') THEN
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          i=Index(ix,iy,iz+1)
          ATPot%Val(i,4)=ATPot%Val(i,4)-ATPot%Val(i,7)
        END DO
      END DO
    END IF
    ATPot%Val(1:ATPot%n,ATPot%Diag)=One-(beta0*dt)*ATPot%Val(1:ATPot%n,ATPot%Diag)
    DO iDiag=1,ATPot%NumDiag
      IF (ATPot%DiagPtr(iDiag)/=0) THEN
        ATPot%Val(1:ATPot%n,iDiag)=-(beta0*dt)*ATPot%Val(1:ATPot%n,iDiag) 
      END IF
    END DO
  END IF
    
CONTAINS
FUNCTION Index(ix,iy,iz)
  INTEGER :: Index,ix,iy,iz
  Index=ix-ix0+1+nxP2*(iy-iy0)+nxP2*nyP2*(iz-iz0)  
END FUNCTION Index
END SUBROUTINE JacScaleCompute

SUBROUTINE Cond(Rho,RhoV,RhoL,PotM,RhsRhoV,RhsRhoL,RhsPotM)

  REAL(RealKind) :: Rho,RhoV,RhoL,PotM,RhsRhoV,RhsRhoL,RhsPotM

  REAL(RealKind) :: RhoD,RhoI,p,qL,T
  REAL(RealKind) :: pVs,DrvDt,DrlDt,Rm,Cpml,Cvml,Lv

  qL=RhoL/(Rho+Eps)
  RhoI=Zero
  RhoD=Rho-RhoV-RhoL-RhoI+Eps
  p=PressureTheta(RhoD,RhoV,RhoL,RhoI,PotM)+Eps
  T=AbsTemp(RhoD,RhoV,p)+Eps
  pVs=SaturVapor(T)
  DrvDt=RelCloud*((pVs/(Rv*T)-RhoV)+RhoL- &
        SQRT((pVs/(Rv*T)-RhoV)*(pVs/(Rv*T)-RhoV)+RhoL*RhoL))
  DrlDt=-DrvDt
  Rm=Rd*RhoD+Rv*RhoV+Eps
  Cpml=Cpd*RhoD+Cpv*RhoV+Cpl*RhoL+Eps
  Cvml=Cvd*RhoD+Cvv*RhoV+Cpl*RhoL+Eps
  Lv=LatHeat(T)
  RhsPotM=PotM*( &
                (-Lv/(Cpml*T) &  
                 -LOG(p/P0)*(Rm/Cpml)*(Rv/Rm-Cpv/Cpml) &
                 +Rv/Rm                                   &
                 )*DrvDt                                   &
                +(LOG(p/P0)*(Rm/Cpml)*(Cpl/Cpml)        &
                 )*DrlDt                                   &
                )
  RhsRhoV=DrvDt
  RhsRhoL=DrlDt

END SUBROUTINE Cond

SUBROUTINE CloudCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc
  REAL(RealKind) :: RhsRhoV,RhsRhoL,RhsPotM

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        PotMLoc=th(ix,iy,iz,1)
        RhoLoc=Rho(ix,iy,iz,1)
        RhoVLoc=RhoV(ix,iy,iz,1)
        RhoLLoc=RhoL(ix,iy,iz,1)
        CALL Cond(RhoLoc,RhoVLoc,RhoLLoc,PotMLoc,RhsRhoV,RhsRhoL,RhsPotM) 
        thRhs(ix,iy,iz,1)=thRhs(ix,iy,iz,1)+RhsPotM
        RhoVRhs(ix,iy,iz,1)=RhoVRhs(ix,iy,iz,1)+RhsRhoV 
        RhoLRhs(ix,iy,iz,1)=RhoLRhs(ix,iy,iz,1)+RhsRhoL 
      END DO
    END DO
  END DO
        
END SUBROUTINE CloudCompute

SUBROUTINE JacCloudCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc
  REAL(RealKind) :: RhsRhoV,RhsRhoL,RhsPotM
  REAL(RealKind) :: RhsRhoV1,RhsRhoL1,RhsPotM1
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        PotMLoc=th(ix,iy,iz,1)
        RhoLoc=Rho(ix,iy,iz,1)
        RhoVLoc=RhoV(ix,iy,iz,1)
        RhoLLoc=RhoL(ix,iy,iz,1)
        CALL Cond(RhoLoc,RhoVLoc,RhoLLoc,PotMLoc,RhsRhoV,RhsRhoL,RhsPotM)
        Temp=PotMLoc
        PotMLoc=Temp*(1.0d0+1.d-8)
        CALL Cond(RhoLoc,RhoVLoc,RhoLLoc,PotMLoc,RhsRhoV1,RhsRhoL1,RhsPotM1)
        AS(IndexMet(thPosJac,thPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(thPosJac,thPosJac))%c(ix,iy,iz,1)+(RhsPotM1-RhsPotM)/(PotMloc-Temp+Eps)
        AS(IndexMet(RhoVPosJac,thPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(RhoVPosJac,thPosJac))%c(ix,iy,iz,1)+(RhsRhoV1-RhsRhoV)/(PotMloc-Temp+Eps)
        AS(IndexMet(RhoCPosJac,thPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(RhoCPosJac,thPosJac))%c(ix,iy,iz,1)+(RhsRhoL1-RhsRhoL)/(PotMloc-Temp+Eps)
        PotMLoc=Temp

        Temp=RhoVLoc
        RhoVLoc=Temp*(1.0d0+1.d-8)
        CALL Cond(RhoLoc,RhoVLoc,RhoLLoc,PotMLoc,RhsRhoV1,RhsRhoL1,RhsPotM1)
        AS(IndexMet(thPosJac,RhoVPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(thPosJac,RhoVPosJac))%c(ix,iy,iz,1)+(RhsPotM1-RhsPotM)/(RhoVLoc-Temp+Eps)
        AS(IndexMet(RhoVPosJac,RhoVPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(RhoVPosJac,RhoVPosJac))%c(ix,iy,iz,1)+(RhsRhoV1-RhsRhoV)/(RhoVLoc-Temp+Eps)
        AS(IndexMet(RhoCPosJac,RhoVPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(RhoCPosJac,RhoVPosJac))%c(ix,iy,iz,1)+(RhsRhoL1-RhsRhoL)/(RhoVLoc-Temp+Eps)
        RhoVLoc=Temp

        Temp=RhoLLoc
        RhoLLoc=Temp*(1.0d0+1.d-8)
        CALL Cond(RhoLoc,RhoVLoc,RhoLLoc,PotMLoc,RhsRhoV1,RhsRhoL1,RhsPotM1)
        AS(IndexMet(thPosJac,RhoCPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(thPosJac,RhoCPosJac))%c(ix,iy,iz,1)+(RhsPotM1-RhsPotM)/(RhoLLoc-Temp+Eps)
        AS(IndexMet(RhoVPosJac,RhoCPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(RhoVPosJac,RhoCPosJac))%c(ix,iy,iz,1)+(RhsRhoV1-RhsRhoV)/(RhoLLoc-Temp+Eps)
        AS(IndexMet(RhoCPosJac,RhoCPosJac))%c(ix,iy,iz,1)= &
        AS(IndexMet(RhoCPosJac,RhoCPosJac))%c(ix,iy,iz,1)+(RhsRhoL1-RhsRhoL)/(RhoLLoc-Temp+Eps)
        RhoLLoc=Temp
      END DO
    END DO
  END DO

END SUBROUTINE JacCloudCompute

SUBROUTINE TurbTke

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FluxL,FluxR

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        FluxL=2.0d0/3.0d0*FU(ix-1,iy,iz)*tke(ix,iy,iz,1)
        uRhsR(ix-1,iy,iz,1)=uRhsR(ix-1,iy,iz,1)-FluxL
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)-FluxL
        FluxR=2.0d0/3.0d0*FU(ix,iy,iz)*tke(ix,iy,iz,1)
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+FluxR
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)+FluxR

        FluxL=2.0d0/3.0d0*FV(ix,iy-1,iz)*tke(ix,iy,iz,1)
        vRhsR(ix,iy-1,iz,1)=vRhsR(ix,iy-1,iz,1)-FluxL
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)-FluxL
        FluxR=2.0d0/3.0d0*FV(ix,iy,iz)*tke(ix,iy,iz,1)
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)+FluxR
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)+FluxR

        FluxL=2.0d0/3.0d0*FW(ix,iy,iz-1)*tke(ix,iy,iz,1)
        wRhsR(ix,iy,iz-1,1)=wRhsR(ix,iy,iz-1,1)-FluxL
        wRhsL(ix,iy,iz,1)=wRhsL(ix,iy,iz,1)-FluxL
        FluxR=2.0d0/3.0d0*FW(ix,iy,iz)*tke(ix,iy,iz,1)
        wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)+FluxR
        wRhsL(ix,iy,iz+1,1)=wRhsL(ix,iy,iz+1,1)+FluxR
      END DO
    END DO
  END DO
  IF (TypeW(1:1)=='o') THEN
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        FluxR=2.0d0/3.0d0*FU(ix,iy,iz)*tke(ix,iy,iz,1)
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)+FluxR
        uRhsR(ix,iy,iz,1)=Zero
      END DO
    END DO
  END IF
  IF (TypeE(1:1)=='o') THEN
    ix=ix1+1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        FluxL=2.0d0/3.0d0*FU(ix-1,iy,iz)*tke(ix,iy,iz,1)
        uRhsR(ix-1,iy,iz,1)=uRhsR(ix-1,iy,iz,1)-FluxL
        uRhsL(ix,iy,iz,1)=Zero
      END DO
    END DO
  END IF
  IF (TypeS(1:1)=='o') THEN
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        FluxR=2.0d0/3.0d0*FV(ix,iy,iz)*tke(ix,iy,iz,1)
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)+FluxR
        vRhsR(ix,iy,iz,1)=Zero
      END DO
    END DO
  END IF
  IF (TypeN(1:1)=='o') THEN
    iy=iy1+1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        FluxL=2.0d0/3.0d0*FV(ix,iy-1,iz)*tke(ix,iy,iz,1)
        vRhsR(ix,iy-1,iz,1)=vRhsR(ix,iy-1,iz,1)-FluxL
        vRhsL(ix,iy,iz,1)=Zero
      END DO
    END DO
  END IF

  IF (TypeB(1:1)=='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        FluxR=2.0d0/3.0d0*FW(ix,iy,iz)*tke(ix,iy,iz,1)
        wRhsL(ix,iy,iz+1,1)=wRhsL(ix,iy,iz+1,1)+FluxR
        wRhsR(ix,iy,iz,1)=Zero
      END DO
    END DO
  END IF
  IF (TypeT(1:1)=='o') THEN
    iz=iz1+1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        FluxL=2.0d0/3.0d0*FW(ix,iy,iz-1)*tke(ix,iy,iz,1)
        wRhsR(ix,iy,iz-1,1)=wRhsR(ix,iy,iz-1,1)-FluxL
        wRhsL(ix,iy,iz,1)=Zero
      END DO
    END DO
  END IF

END SUBROUTINE TurbTke

SUBROUTINE ForceVelCompute(Time)

  REAL(RealKind) :: Time
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: xPL,yPL,zPL

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        xPL=xP(ix-1)+0.5e0*dx(ix)
        yPL=yP(iy-1)+0.5e0*dy(iy)
        zPL=zP(iz-1)+0.5e0*dz(iz)
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)+ForceU(xPL,yPL,zPL,zH(ix,iy),Time)
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)+ForceU(xPL,yPL,zPL,zH(ix,iy),Time)
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)+ForceV(xPL,yPL,zPL,zH(ix,iy),Time)
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)+ForceV(xPL,yPL,zPL,zH(ix,iy),Time)
        wRhsL(ix,iy,iz,1)=wRhsL(ix,iy,iz,1)+ForceW(xPL,yPL,zPL,zH(ix,iy),Time)
        wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)+ForceW(xPL,yPL,zPL,zH(ix,iy),Time)
      END DO
    END DO
  END DO

END SUBROUTINE ForceVelCompute

SUBROUTINE ForceScalarCompute(Rhs,Time,Tend)

  TYPE(Vec4_T), POINTER :: Rhs(:)
  REAL(RealKind) :: Time
  TYPE(Vec1_T), POINTER, OPTIONAL :: Tend(:)
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: xPL,yPL,zPL
  INTEGER :: x0p,x1p,y0p,y1p
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoILoc,RhoDLoc,TLoc,pLoc,eLoc
  REAL(RealKind) :: Rm,Cpml,Cvml,pvs,Lv,Cp_eff,LWPLoc
  REAL(RealKind) :: LWP(iz0+1:iz1)
  REAL(RealKind) :: drvdt,rrv,rrl,dPotdt,PotM,ThDens,ThEquiv,Kinetic,Energy 
  REAL(RealKind) :: DpDRho,DpDe,DpDRhoV,dFdzTop,dFdzBottom,dFdz 
  REAL(RealKind) :: Theta,LWFlux,dTdt,LWTop,LWBottom,DeltaLW 
  REAL(RealKind) :: PreFacRhoV,PreFacRhoVEn
  REAL(RealKind) :: F_0=72.0d0
  Real(RealKind) :: F_1=15.0d0
  REAL(RealKind) :: k=170.0d0 
  REAL(RealKind) :: RhoVForce(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: ThForce(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)

  
  IF (RadiationISDAC) THEN
!   ALLOCATE(LWP(iz1))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        LWPLoc=0.0d0
        LWP=0.0d0
        DO iz=iz0+1,iz1
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          LWPLoc=LWPLoc+RhoLLoc*dz(iz) !LWP in height z (local)
          LWP(iz)=LWPLoc
        END DO
        DO iz=iz0+1,iz1
          RhoLoc=Rho(ix,iy,iz,1)
          IF(iz==iz0+1)THEN
            dFdz=((F_0*EXP(-k*(LWP(iz1)-LWP(iz+1)))+F_1*EXP(-k*LWP(iz+1)))-&
                 (F_0*EXP(-k*(LWP(iz1)-LWP(iz)))+F_1*EXP(-k*LWP(iz))))/dz(iz+1)
          ELSE IF(iz==iz1)THEN
            dFdz=((F_0*EXP(-k*(LWP(iz1)-LWP(iz)))+F_1*EXP(-k*LWP(iz)))-&
                 (F_0*EXP(-k*(LWP(iz1)-LWP(iz-1)))+F_1*EXP(-k*LWP(iz-1))))/dz(iz)
          ELSE  
            dFdzTop   =((F_0*EXP(-k*(LWP(iz1)-LWP(iz+1)))+F_1*EXP(-k*LWP(iz+1)))-&
                        (F_0*EXP(-k*(LWP(iz1)-LWP(iz)))+F_1*EXP(-k*LWP(iz))))/dz(iz+1)
            dFdzBottom=((F_0*EXP(-k*(LWP(iz1)-LWP(iz)))+F_1*EXP(-k*LWP(iz)))-&
                        (F_0*EXP(-k*(LWP(iz1)-LWP(iz-1)))+F_1*EXP(-k*LWP(iz-1))))/dz(iz)
            dFdz=(dFdzTop+dFdzBottom)/2.0d0
          END IF
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoILoc=RhoI(ix,iy,iz,1)+RhoS(ix,iy,iz,1)
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc-RhoILoc
          TLoc=T(ix,iy,iz,1)
          Theta=Th(ix,iy,iz,1)/(RhoLoc+Eps)
          dTdt=-(One/(CpmlFun(RhoDLoc,RhoVLoc,RhoLLoc,RhoILoc)))*dFdz  
          HeatRate(ix,iy,iz,1)=dTdt*Theta/(TLoc+Eps)
          Rhs(ThPos)%c(ix,iy,iz,1)=Rhs(ThPos)%c(ix,iy,iz,1) &
                                   +dTdt*Theta/(TLoc+Eps)*RhoLoc
        END DO
      END DO
    END DO
  END IF


  IF (PRESENT(Tend)) THEN
    DO iz=iz0+1,iz1
      RhoVForce(:,:,iz)=Tend(2)%c(iz)*Rho(ix0+1:ix1,iy0+1:iy1,iz,1)
      ThForce(:,:,iz)=Tend(1)%c(iz)*Rho(ix0+1:ix1,iy0+1:iy1,iz,1)
    ! RhoVForcing(ix,iy,iz,1)=Tend(2)%c(iz)
    ! ThForcing(ix,iy,iz,1)=Tend(1)%c(iz)
    END DO  
  ELSE
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          RhoVForce(ix,iy,iz)=ForceRho(xPL,yPL,zPL,zH(ix,iy),Time)*Rho(ix,iy,iz,1)
          ThForce(ix,iy,iz)=ForceTh(xPL,yPL,zPL,zH(ix,iy),Time)*Rho(ix,iy,iz,1)
        END DO      
      END DO      
    END DO      
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        IF (RhoPos>0) THEN
          Rhs(RhoPos)%c(ix,iy,iz,1)=Rhs(RhoPos)%c(ix,iy,iz,1) &
                       +RhoVForce(ix,iy,iz)
        END IF
        IF (RhoVPos>0) THEN
          Rhs(RhoVPos)%c(ix,iy,iz,1)=Rhs(RhoVPos)%c(ix,iy,iz,1) &
                       +RhoVForce(ix,iy,iz)
        END IF
        IF (ThPos>0.AND.PRESENT(Tend)) THEN
          Rhs(ThPos)%c(ix,iy,iz,1)=Rhs(ThPos)%c(ix,iy,iz,1) &
                       +ThForce(ix,iy,iz)
        ELSE IF (ThPos>0) THEN
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc
          rrv=RhoVLoc/RhoDLoc
          rrl=RhoLLoc/RhoDLoc
          pLoc=p(ix,iy,iz,1)
          TLoc=T(ix,iy,iz,1)
          Theta=Th(ix,iy,iz,1)
          pVs=SaturVapor(TLoc)
          Lv=LatHeat(TLoc)
          Rm=RhoDLoc*Rd+RhoVLoc*Rv
          Cpml=RhoDLoc*Cpd+RhoVLoc*Cpv+RhoLLoc*Cpl
          Cp_eff=Cpd+(rrv+rrl)*Cpl
          PotM=TLoc*(p0/pLoc)**(Rm/Cpml)
          SELECT CASE(ThetaKind)
          CASE('Density')
            PreFacRhoV=Theta*(Rv/Rm &
                            +LOG(p0/((RhoDLoc*Rd+RhoVLoc*Rv)*TLoc)) &
                             *(-cpv*Rd*RhoDLoc+cpd*Rv*RhoDLoc+cpl*Rv*RholLoc)/Cpml**2.0d0  &
                           )
            Rhs(ThPos)%c(ix,iy,iz,1)=Rhs(ThPos)%c(ix,iy,iz,1) &
                                     +ThForce(ix,iy,iz)*RhoLoc
          CASE('Equiv')
            PreFacRhoV=Theta/(RhoDLoc*Cpd+(RhoVLoc+RhoLLoc)*Cpl)**2   &
                      *(-RhoDLoc*Rd*Cpl*LOG(p0/(RhoDLoc*Rd*TLoc))   &
                        +(RhoDLoc*Cpd+RhoLLoc*Cpl)*(Lv/TLoc-Rv*LOG(RelHumidity(TLoc,RhoVLoc))) &
                       ) &
                     +Theta/RhoLoc
          CASE('Energy')
            Kinetic=KinEn(ix,iy,iz,1)
            PreFacRhoV=(Cpv*TLoc+L00)+Grav*zPL+Kinetic
          CASE('PreEn')
            Cvml=RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl
            eLoc=Cvml*TLoc+RhoVLoc*L00
            DpDRho=Rd*(eLoc-RhoVLoc*L00)/Cvml &
                   -Cvd*(eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2 &
                   +eLoc/(RhoLoc+Eps)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml
            DpDe=RhoLoc*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml
            DpDRhoV=(Rv-Rd)*(eLoc-RhoVLoc*L00)/Cvml &
                    -L00*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml &
                    -(Cvv-Cvd)*(eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2 
            PreFacRhoV=Cpv*TLoc+L00
            Rhs(enPos)%c(ix,iy,iz,1)=Rhs(enPos)%c(ix,iy,iz,1) &
                       +RhoVForce(ix,iy,iz)*PreFacRhoV
            PreFacRhoV=DpDRho+(Cpv*TLoc+L00)/(RhoLoc+Eps)*DpDe+DpDRhoV
          END SELECT  
          Rhs(ThPos)%c(ix,iy,iz,1)=Rhs(ThPos)%c(ix,iy,iz,1) &
                     +RhoVForce(ix,iy,iz)*PreFacRhoV
        END IF
        IF (ThPos>0.AND.Shallow) THEN
          Rhs(ThPos)%c(ix,iy,iz,1)=Rhs(ThPos)%c(ix,iy,iz,1) &
                       +RhoVForce(ix,iy,iz)
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE ForceScalarCompute

! Subroutines for cell perturbation method
! to generate turbulent inflow conditions
! when using non-cyclic boundary conditions

SUBROUTINE CellPert(Rhs,TimeMain)

  TYPE(Vec4_T), POINTER :: Rhs(:)
  REAL(RealKind) :: TimeMain
  INTEGER :: ix,iy,iz
  INTEGER :: x0p,x1p,y0p,y1p,BlockSize!,BlockNum
  INTEGER,PARAMETER :: BlockNum=1024
  INTEGER :: BlockX0(1:BlockNum),BlockX1(1:BlockNum),BlockY0(1:BlockNum),BlockY1(1:BlockNum)
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDLoc,TLoc,pLoc,eLoc
  REAL(RealKind) :: Rm,Cpml,Cvml,pvs,Lv,Cp_eff
  REAL(RealKind) :: drvdt,rrv,rrl,dPotdt,PotM,ThDens,ThEquiv,Kinetic,Energy
  REAL(RealKind) :: DpDRho,DpDe,DpDRhoV
  REAL(RealKind) :: Theta
  REAL(RealKind) :: PreFacRhoV,PreFacRhoVEn
  INTEGER :: iostat
  LOGICAL :: exist

  ! Check if perturbation file exists
  INQUIRE(FILE=ForcingCellPertFile, EXIST=exist)
  IF (exist) then
    OPEN(99,FILE=ForcingCellPertFile,STATUS='OLD',ACTION='READ',IOSTAT=iostat)
    READ(99,*) x0p,y0p
    READ(99,*) x1p,y1p
    READ(99,*) BlockSize
    CLOSE(99)
  ELSE
    WRITE(*,*) 'File ',ForcingCellPertFile,' does not exist. Stopped.'
    STOP
  END IF
  CALL CellPertDefineBlocks(BlockSize,BlockNum,BlockX0,BlockY0,BlockX1,BlockY1,x0P,y0P,x1P,y1P)
  CALL CellPertCompute(Rhs,BlockSize,BlockNum,BlockX0,BlockY0,BlockX1,BlockY1,x0P,y0P,x1P,y1P,TimeMain)

END SUBROUTINE CellPert


SUBROUTINE CellPertDefineBlocks(BlockSize,BlockNum,BlockX0,BlockY0,BlockX1,BlockY1,x0P,y0P,x1P,y1P)

  INTEGER :: BlockSize,BlockNum
  INTEGER :: BlockX0(1:BlockNum),BlockX1(1:BlockNum),BlockY0(1:BlockNum),BlockY1(1:BlockNum)
  INTEGER :: x0P,y0P,x1P,y1P

  INTEGER :: nBlock,iBlock,jBlock

  nBlock=0
 !DO iBlock=1,3
    iblock=1
    DO jBlock=1,(y1p-y0p)/BlockSize
      nBlock=nBlock+1
      BlockX0(nBlock)=x0p+BlockSize*(iBlock-1)
      BlockY0(nBlock)=y0p+BlockSize*(jBlock-1)
      BlockX1(nBlock)=x0p+BlockSize*(iBlock)
      BlockY1(nBlock)=y0p+BlockSize*(jBlock)
    END DO
 !END DO

END SUBROUTINE CellPertDefineBlocks


SUBROUTINE CellPertCompute(Rhs,BlockSize,BlockNum,BlockX0,BlockY0,BlockX1,BlockY1,x0P,y0P,x1P,y1P,TimeMain)

  TYPE(Vec4_T), POINTER :: Rhs(:)
  INTEGER :: BlockSize,BlockNum
  INTEGER :: BlockX0(1:BlockNum),BlockX1(1:BlockNum),BlockY0(1:BlockNum),BlockY1(1:BlockNum)
  INTEGER :: x0P,y0P,x1P,y1P
  REAL(RealKind) :: TimeMain

  INTEGER :: ix,iy,iz,i
  REAL(RealKind) :: zPL
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: Perturbation
  REAL(RealKind) :: Ran(1:200000)
  INTEGER :: nBlock,iBlock,jBlock
  REAL(RealKind) :: r(BlockNum)
  INTEGER :: NumTime,RanPos
  INTEGER :: iostat
  LOGICAL :: exist

  ! Check if random number file exists
  INQUIRE(FILE='Random_number.txt', EXIST=exist)
  IF (exist) then
    OPEN(99,FILE='Random_number.txt',STATUS='OLD',ACTION='READ',IOSTAT=iostat)
    DO i=1,200000  ! 200000 random numbers (uniform distribution), one per line
      READ(99,*) Ran(i)
    END DO
    CLOSE(99)
  ELSE
    WRITE(*,*) 'File Random_number.txt does not exist. Stopped.'
    STOP
  END IF
  ! Counting current number of cell perturbation
  NumTime=INT(TimeMain)/INT(ForcingCellPertTime)
  DO iz=iz0+1,iz1
    zPL=zP(iz-1)+0.5e0*dz(iz)
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        IF (ThPos>0) THEN
          nBlock=0
          iBlock=1
          DO jBlock=1,(y1p-y0p)/BlockSize
            nBlock=nBlock+1
            IF (ix>=BlockX0(nBlock).AND.ix<BlockX1(nBlock).AND.iy>=BlockY0(nBlock).AND.iy<BlockY1(nBlock) &
                .AND.zPL<=ForcingCellPertHeight) THEN
              RanPos=MOD((NumTime+1)*nBlock,200000)
              r(nBlock)=Ran(RanPos)
              RhoLoc=Rho(ix,iy,iz,1)
              Perturbation=ForcingCellPertTh*Two*(r(nBlock)-0.5d0)/dtMax
              Rhs(ThPos)%c(ix,iy,iz,1)=Rhs(ThPos)%c(ix,iy,iz,1) &
                                       +RhoLoc*Perturbation
            END IF
          END DO
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE CellPertCompute

SUBROUTINE GetRadiationValues(Time)

  REAL(RealKind) :: Time
  
  REAL(RealKind) :: eko_dir_sw,eko_dif_sw,kip_dir_sw,kip_dif_sw,kip_glo_lw
  CHARACTER(8) :: time_char
  INTEGER :: time_sec,time_h,time_m,time_s,rows
  INTEGER :: i,j,r
  CHARACTER(124) :: tmp
  LOGICAL :: exist

  IF ((Time==StartTime.OR.MOD(INT(Time-StartTime),INT(RadiationValuesTime))==0).AND.Time<EndTime-RadiationValuesTime) THEN
    INQUIRE(FILE=RadiationValuesFile,EXIST=exist)
    IF (exist) THEN
      OPEN(UNIT=12,FILE=RadiationValuesFile,form="formatted")
      DO j=1,3
        READ(12,*) tmp ! read header (3 lines)
      END DO
      IF (RadiationValuesTime==60) THEN
        rows=1440
      END IF
      IF (RadiationValuesTime==600) THEN
        rows=144
      END IF
      DO r=1,rows
        READ(12,*) time_char,eko_dir_sw,eko_dif_sw, &
                   kip_dir_sw,kip_dif_sw,kip_glo_lw
        READ(time_char(1:2),'(i2)') time_h
        READ(time_char(4:5),'(i2)') time_m
        READ(time_char(7:8),'(i2)') time_s
        time_sec = 3600 * time_h & ! hours
                  +  60 * time_m & ! minutes
                  +       time_s   ! seconds
        IF (time_sec==INT(Time)) THEN
          IF (MyID==0) WRITE (*,*) 'Time, radiation values (dir,dif,lw):',time_sec,kip_dir_sw,kip_dif_sw,kip_glo_lw
          DO ibLoc=1,nbLoc 
            ib=LocGlob(ibLoc)
            CALL Set(Floor(ib))
            DO i=1,NumBoundCell 
              BoundCell(i)%raddirekt=kip_dir_sw
              BoundCell(i)%raddiffus=kip_dif_sw
              BoundCell(i)%radinfred=kip_glo_lw
              RaddirCell(ibLoc)%cB(i,1)=BoundCell(i)%raddirekt
              RaddifCell(ibLoc)%cB(i,1)=BoundCell(i)%raddiffus
              RadinfCell(ibLoc)%cB(i,1)=BoundCell(i)%radinfred
            END DO
          END DO
          EXIT
        END IF
      END DO
      CLOSE(12)
    ELSE
      WRITE(*,*) 'RadiationValuesFile ',RadiationValuesFile,' does not exist!'
      STOP
    END IF
  END IF

END SUBROUTINE GetRadiationValues

SUBROUTINE ScaleVectorCell(x)

  TYPE(Vector4Cell_T) :: x(:)

  INTEGER :: ic

  CALL ExchangeFlux(x)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    DO ic=1,SIZE(x(ibLoc)%Vec)
      f=>x(ibLoc)%Vec(ic)%c
      CALL AdvectionScale
    END DO
  END DO
END SUBROUTINE ScaleVectorCell

SUBROUTINE ScaleScalarCell(x)

  TYPE(ScalarCell_T) :: x(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    f=>x(ibLoc)%c
    CALL AdvectionScale
  END DO

END SUBROUTINE ScaleScalarCell

SUBROUTINE UpdateScalarCell(alpha,x,y)
  
  REAL(RealKind) :: alpha
  TYPE(ScalarCell_T) :: x(:)
  TYPE(ScalarCell_T) :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    c=>y(ibLoc)%c
    f=>x(ibLoc)%c
    CALL AdvectionUpdate(alpha)
  END DO

END SUBROUTINE UpdateScalarCell
SUBROUTINE UpdateVectorCell(alpha,x,y)

  REAL(RealKind) :: alpha
  TYPE(Vector4Cell_T) :: x(:)
  TYPE(Vector4Cell_T) :: y(:)

  INTEGER :: ic

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    DO ic=1,SIZE(x(ibLoc)%Vec)
      c=>y(ibLoc)%Vec(ic)%c
      f=>x(ibLoc)%Vec(ic)%c
      CALL AdvectionUpdate(alpha)
    END DO
  END DO

END SUBROUTINE UpdateVectorCell
END MODULE Function_Mod
