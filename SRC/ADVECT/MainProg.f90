PROGRAM MainProg

  USE Domain_Mod
  USE Init_Mod
  USE JacAccGrav_Mod
  USE Int_Mod
! USE IntPeer_Mod
! USE IntLinPeer_Mod
  USE ReadOutput_Mod
  USE Emission_Mod
  USE Koagulation_Mod
  USE ScalarVectorCellPar_Mod
  USE Microphysics_Mod
  USE Activity_Mod
  USE RhoProf_Mod
  USE IntSub_Mod
  USE Turbulence_Mod
  USE Operator_Mod
  USE Parameter_Mod
  USE Physics_Mod
  USE Tools_Mod
  USE SpecialOutput_Mod
  USE Soil_Mod
  USE Canopy_Mod
  USE ReadWeights_Mod
  USE WindFarm_Mod
  USE Forcing_Mod
  USE Restart_Mod
  USE Grid_Mod
  USE Example_Mod
  USE BoundaryCondition_Mod
  USE OutputMeanProfile_Mod
  USE SingleColumnXProfile_Mod
  USE SingleColumnZProfile_Mod
  USE MeanSurface_Mod
  USE PointSurface_Mod 
  USE TimeStep_Mod
  USE ReadWRF_Mod
  USE ReadWRFnc_Mod

  IMPLICIT NONE

  TYPE(VelocityFace_T), POINTER :: VelF1(:)
  TYPE(VecVelocityFace_T), POINTER :: VecVelF(:)
  TYPE(Vector4Cell_T), POINTER :: VecMet(:)
  TYPE(Vector4Cell_T), POINTER :: VecChem(:)
  TYPE(Vector4Cell_T), POINTER :: VecG(:)
  TYPE(VecVector4Cell_T), POINTER :: VecVecMet(:)
  TYPE(Vector4Cell_T), POINTER :: VecMetP(:)
  REAL(RealKind) :: Temp,RhoLoc,rRand
  INTEGER :: Iter,iForce
  INTEGER :: iInt,ix,iy,iz
  INTEGER :: ixVol,iyVol,izVol
  Real(RealKind) :: VolMin
  INTEGER :: i,j,k,iShift
  INTEGER :: iFrac
  INTEGER :: stat
  CHARACTER(80) :: ProfileSTART='ProfileSTART'
  CHARACTER(80) :: ProfileEND='ProfileEND'
  CHARACTER(80) :: InputFileName
  CHARACTER(8) :: ScalarName
  INTEGER :: iFort
  INTEGER :: unr
  INTEGER :: PosDummy
  REAL(RealKind) :: LocTime,MaxTime
  !REAL(RealKind) :: OutputTime !!
  REAL(RealKind) :: dtAct,Time
  REAL(RealKind) :: dtCFL=0.0d0
  INTEGER :: nsCFL

  LOGICAL :: FileExist

  REAL(RealKind) :: TotalStart(1:5),TotalCurrent(1:5)
  INTEGER :: SizeSeed=12
  !INTEGER,Dimension(12) :: seed
  INTEGER, ALLOCATABLE :: seed(:)
  CHARACTER :: ArgRestart

  REAL(RealKind) :: ErrorL1,SumVol
  REAL(RealKind) :: Start,Finish

! Outputs
  LOGICAL :: CheckMean
  LOGICAL :: CheckColumnX
  LOGICAL :: CheckColumnZ
  LOGICAL :: CheckSurfMean
  LOGICAL :: CheckSurfPoint
 
  LOGICAL :: FinishNudging
  LOGICAL :: FinishTendency

  REAL(RealKind) :: tt

  CALL start_MPI

  CALL Random_seed(size=SizeSeed)
! seed = 1 !myID
  ALLOCATE (seed(SizeSeed))
  seed = myID 
  CALL Random_seed(put=seed)

! -- Lesen der Gitterdatei --
  CALL get_command_argument(1,InputFileName)
  CALL get_command_argument(2,ArgRestart)
  IF (ArgRestart=='T') THEN
    Restart=.TRUE.
  ELSE  
    Restart=.FALSE.
  END IF  
  IF (MyId==0) THEN
    IF (InputFileName == ' ') THEN
!      STOP 'Dateiname fehlt'
    END IF
  END IF
  OutputNameSoil=TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))//'.Soil'
  OutputNameCut=TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))//'.Cut'

  CALL MPI_Bcast(InputFileName,80,MPI_CHARACTER,0,MPI_COMM_WORLD,MPIErr)
  IF (MyId==0) THEN
    WRITE(*,*)
    WRITE(*,*) InputFileName
  END IF

! -- Compute Parameter
  CALL ComputeParameter 
! -- Einlesen des Gitters (Block-Struktur) --
  CALL Allocate(Floor)
  CALL inp_part(InputFileName)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL InputExample(InputFileName)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
! CALL Hydrostatic OSSI
  IF (MultiTriTB.OR.MultiTriTR.OR.MultiMuTR) THEN
    CALL SetDim1P(2)
  END IF

  CALL Allocate(Floor)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL ReadWeights(InputFileName)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL MinVol

  CALL InputModelTransport(InputFileName)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)

  IF (ChemieFile/='') THEN 
    WRITE(*,*) 'InputSystem(ChemieFile)',TRIM(ChemieFile)
    CALL InputSystem(ChemieFile)
  END IF

  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL SetIndices
  CALL Allocate(JacTrans)
  JacTrans=Zero
  CALL AllocateVec4Chemie(VecMet,VectorComponentsMet)
  CALL AllocateVec4Chemie(VecChem,VectorComponentsChem)
  WRITE(*,*) 'Nach CALL AllocateVec4Chemie(VecChem,VectorComponentsChem)'

  CALL InitPrandtlNumber(InputFileName)

! FallVelocity
  ALLOCATE(FallVelocity(0:VectorComponentsT))
  FallVelocity=Zero
  IF (RhoRPos>0) THEN
    FallVelocity(RhoRPos)=0.0d0
  END IF
  ALLOCATE(ScaleMat2(0:VectorComponentsT))
  ScaleMat2=FallVelocity

  VectorComponentsM=VectorComponentsT
  CALL InputModelBC(InputFileName)

  IF (BCVel%West  =='MeanFlow'.OR.BCVel%East =='MeanFlow'.OR. &
      BCVel%South =='MeanFlow'.OR.BCVel%North=='MeanFlow'.OR. &
      BCVel%Bottom=='MeanFlow'.OR.BCVel%Top  =='MeanFlow') THEN 
    ALLOCATE(MeanProfile(0:VectorComponentsT,Domain%ix0+1:Domain%ix1))
  END IF

  VecMet=Zero
  VecChem=Zero
!   -- Eingabe Chemie --
  IF (Chemie.AND.DataFile/='') THEN 
  
    CALL InputChemicalData(DataFile)
  END IF

! -- WindFarm--
  IF (Wind.AND.WindFarmFileName/='') THEN 
! -- Input WindFarm --
    CALL InitWind(WindFarmFileName)
  END IF  
   
! -- Fields for 3D output --
  IF (Forcing) THEN
    CALL Allocate(HeatRateCell)
    CALL Allocate(ForceThCell)
    CALL Allocate(ForceRhoVCell)
  END IF
  IF (Radiation.AND..NOT.RadiationProfile) THEN
    CALL Allocate(RaddirCell)
    CALL Allocate(RaddifCell)
    CALL Allocate(RadinfCell)
    CALL Allocate(ShadowCell)
  END IF

! -- Fields for 2D output --
  IF (DynamicSoil) THEN
    CALL Allocate(SensFluxCell)
    CALL Allocate(LatFluxCell)
    CALL Allocate(LandClassCell)
    CALL Allocate(RoughnessLengthCell)
    CALL Allocate(AlbedoCell)
    CALL Allocate(BulkCoeffDragCell)
    CALL Allocate(BulkCoeffHeatCell)
    CALL Allocate(BulkCoeffMoistCell)
  END IF


! -- Meteorologie --
  WRITE(*,*) 'Meteorologie '
  CALL Allocate(VelF1)
  VelF1=Zero
  VecMet=Zero

  Time=StartTime
  IF (Shallow) THEN
    CALL Allocate(HeightG)
    CALL ScalarInit(HeightG,HeightFun,Time)
  END IF
! CALL Allocate(RhoCell)
  ALLOCATE(DummyCell(nbLoc))
  ALLOCATE(DummyVec(nbLoc))
  DO ibLoc=1,nbLoc
    ALLOCATE(DummyVec(ibLoc)%Vec(2))
  END DO
  CALL Allocate(RhoProfG)
  IF (RhoPos>0) THEN
    ALLOCATE(RhoCell(nbLoc))
    CALL Assign(RhoCell,VecMet,RhoPos)
    CALL ScalarInit(RhoProfG,RhoProf,Time)
    CALL ExchangeCell(RhoProfG)
    CALL VectorInit(RhoPos,VecMet,RhoFun,Time)
    CALL Allocate(PreCell)
    CALL Allocate(DivCell)
    CALL Allocate(EStartCell)
    CALL Allocate(PStartCell)
    IF (ThetaKind=='PreEn'.OR.ThetaKind=='Exner'.OR.Parcel) THEN
      CALL Allocate(SoundCell)
    END IF  
    IF (ThetaKind=='Energy') THEN
      CALL Allocate(PreKin)
    END IF  
  ELSE IF (Anelastic) THEN
    ALLOCATE(RhoCell(nbLoc))
    CALL ScalarInit(RhoCell,RhoFun,Time)
    CALL ExchangeCell(RhoCell)
    CALL Allocate(VecG,1)
    CALL ScalarInit(RhoProfG,RhoFun,Time)
    CALL ExchangeCell(RhoProfG)
    CALL Allocate(RhoCell)
    CALL ScalarInit(RhoCell,RhoFun,Time)
    CALL ExchangeCell(RhoCell)
    CALL Allocate(PreCell)
  ELSE
    CALL Allocate(RhoCell)
    RhoCell=1.0d0
  END IF

  IF ( uPosL>0)  THEN
    CALL VectorInit(uPosL,VecMet,UStart,Time)
    CALL Mult(RhoCell,VecMet,uPosL)
    CALL Allocate(uCell)
  END IF
  IF ( uPosR>0)  THEN
    CALL VectorInit(uPosR,VecMet,UStart,Time)
    CALL Mult(RhoCell,VecMet,uPosR)
  END IF
  IF ( vPosL>0) THEN
    CALL VectorInit(vPosL,VecMet,vStart,Time)
    CALL Mult(RhoCell,VecMet,vPosL)
    CALL Allocate(vCell)
  END IF
  IF ( vPosR>0) THEN
    CALL VectorInit(vPosR,VecMet,vStart,Time)
    CALL Mult(RhoCell,VecMet,vPosR)
  END IF
  IF ( wPosL>0) THEN
    CALL VectorInit(wPosL,VecMet,wStart,Time)
    CALL Mult(RhoCell,VecMet,wPosL)
    CALL Allocate(wCell)
  END IF
  IF ( wPosR>0) THEN
    CALL VectorInit(wPosR,VecMet,wStart,Time)
    CALL Mult(RhoCell,VecMet,wPosR)
  END IF
  IF ( thPos>0) THEN
    CALL VectorInit(thPos,VecMet,ThStart,Time)
!   CALL BoundaryInit('Skin',TStart,Time)
    CALL Mult(RhoCell,VecMet,thPos)
    CALL Allocate(thProfG)
    CALL ScalarInit(thProfG,thProfFun,Time)
    CALL Mult(RhoProfG,thProfG)
    CALL ExchangeCell(thProfG)
    CALL Allocate(thProfG)
    CALL Allocate(TAbsCell,1,1)
    CALL Allocate(TCell)
    CALL Allocate(CpmlCell)
    TAbsCell=300.0e0_RealKind
    CALL Allocate(KinEnCell)
    CALL Allocate(ECell)
    KinEnCell=Zero
    IF (DynamicSoil) THEN
      CALL InitSoil
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DEALLOCATE(VecMet(ibLoc)%Vec(thPos)%cB)
        ALLOCATE(VecMet(ibLoc)%Vec(thPos)%cB(1:NumBoundCell,1:Domain%nrsoillayers)) 
        DEALLOCATE(VecMet(ibLoc)%Vec(RhoVPos)%cB)
        ALLOCATE(VecMet(ibLoc)%Vec(RhoVPos)%cB(1:NumBoundCell,1:Domain%nrsoillayers+1))
      END DO
      CALL VectorInitSoilFunction(thPos,VecMet,ThStartSoil,Time)
    ELSE IF (Canopy) THEN
      CALL InitCanopy
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DEALLOCATE(VecMet(ibLoc)%Vec(thPos)%cB)
        ALLOCATE(VecMet(ibLoc)%Vec(thPos)%cB(1:NumBoundCell,1:Domain%nrsoillayers)) 
        DEALLOCATE(VecMet(ibLoc)%Vec(RhoVPos)%cB)
        ALLOCATE(VecMet(ibLoc)%Vec(RhoVPos)%cB(1:NumBoundCell,1:Domain%nrsoillayers+1))
      END DO
      CALL VectorInitCanopyFunction(thPos,VecMet,ThStartSoil,Time)
    END IF
  END IF
  IF (ThPos==0.AND.Chemie) THEN
    CALL Allocate(TAbsCell,1,1)
    TAbsCell=300.0d0
    CALL Allocate(PreCell)
    PreCell=1.0d5
  END IF  
    
  IF (EnPos>0) THEN
    CALL VectorInit(EnPos,VecMet,EnStart,Time)
    CALL Mult(RhoCell,VecMet,EnPos)
  END IF  
  IF ( RhoVPos>0) THEN
    CALL Allocate(RhoVProfG)
    CALL ScalarInit(RhoVProfG,QvStart,Time)
    CALL ExchangeCell(RhoVProfG)
    CALL Mult(RhoProfG,RhoVProfG)
    CALL VectorInit(RhoVPos,VecMet,QvStart,Time)
    IF (DynamicSoil) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
      CALL VectorInitSoilFunction(RhoVPos,VecMet,QvStartSoil,Time)
    ELSE IF (Canopy) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
      CALL VectorInitCanopyFunction(RhoVPos,VecMet,QvStartSoil,Time)
    END IF
  END IF
  IF ( RhoCPos>0) THEN
    CALL Allocate(RhoCProfG)
    RhoCProfG=Zero
    CALL ExchangeCell(RhoCProfG)
    CALL Mult(RhoProfG,RhoCProfG)
    CALL VectorInit(RhoCPos,VecMet,QcStart,Time)
    CALL Mult(RhoCell,VecMet,RhoCPos)
  ELSE 
    CALL Allocate(RhoLCell)
  END IF
  IF ( RhoRPos>0) THEN
    CALL VectorInit(RhoRPos,VecMet,QrStart,Time)
    CALL Mult(RhoCell,VecMet,RhoRPos)
    IF (PrecipRain.AND.RainSurf) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
    END IF
  END IF
  IF ( RhoIPos>0) THEN
    CALL VectorInit(RhoIPos,VecMet,QiStart,Time)
    CALL Mult(RhoCell,VecMet,RhoIPos)
    IF (PrecipIce.AND.IceSurf) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
    END IF
  END IF
  IF ( RhoSPos>0) THEN
    CALL VectorInit(RhoSPos,VecMet,QsStart,Time)
    CALL Mult(RhoCell,VecMet,RhoSPos)
    IF (PrecipSnow.AND.SnowSurf) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
    END IF
  END IF
  IF ( nvPos>0) THEN
    CALL VectorInit(nvPos,VecMet,nvStart,Time)
  END IF
  IF ( ncPos>0) THEN
    CALL VectorInit(ncPos,VecMet,ncStart,Time)
  END IF
  IF ( nrPos>0) THEN
    CALL VectorInit(nrPos,VecMet,nrStart,Time)
  END IF
  IF ( niPos>0) THEN
    CALL VectorInit(niPos,VecMet,niStart,Time)
  END IF
  IF ( nsPos>0) THEN
    CALL VectorInit(nsPos,VecMet,nsStart,Time)
  END IF
  IF (tkePos>0) THEN
    CALL VectorInit(tkePos,VecMet,TkeStart,Time)
    CALL Mult(RhoCell,VecMet,tkePos)
  END IF
  IF (tkeHPos>0) THEN
    CALL VectorInit(tkeHPos,VecMet,TkeHStart,Time)
    CALL Mult(RhoCell,VecMet,tkeHPos)
  END IF
  IF (tkeVPos>0) THEN
    CALL VectorInit(tkeVPos,VecMet,TkeVStart,Time)
    CALL Mult(RhoCell,VecMet,tkeVPos)
  END IF
  IF (disPos>0) THEN
    CALL VectorInit(disPos,VecMet,DisStart,Time)
    CALL Mult(RhoCell,VecMet,disPos)
  END IF
  IF (omePos>0) THEN
    CALL VectorInit(omePos,VecMet,OmeStart,Time)
    CALL Mult(RhoCell,VecMet,omePos)
  END IF
  IF (LenPos>0) THEN
    CALL VectorInit(LenPos,VecMet,LenStart,Time)
    CALL Mult(RhoCell,VecMet,LenPos)
  END IF
  IF (Tracer1Pos>0) THEN
    CALL VectorInit(Tracer1Pos,VecMet,Tracer1Start,Time)
  END IF
  IF (Tracer2Pos>0) THEN
    CALL VectorInit(Tracer2Pos,VecMet,Tracer2Start,Time)
  END IF
  CALL PerturbProfile(VecMet)
  CALL ExchangeCell(VecMet)
  IF (uPosL*uPosR>0) THEN
    CALL VelocityCellToFaceLR(VecMet,VelF1,VelF1,Time)
    CALL BoundaryVelocity(VelF1,Time)
  ELSE
    CALL VelocityInit(VelF1,UStart,VStart,WStart,Time)
  END IF

  IF (Parcel) THEN
    CALL InputAdiabatic(InputFileName,VecMet)
  END IF
  IF (Baum.OR.BaumFast) THEN 
    CALL ReadBaum(InputFileName) 
  END IF
  IF (Aerosol.AND.Depos) THEN
    CALL Allocate(SediCell,nFrac)
  END IF

! Environment
  VectorComponentsME=MIN(VectorComponentsME,VectorComponentsM)
  CALL Allocate(VecEnv1,VectorComponentsME)
  CALL Allocate(VecEnv2,VectorComponentsME)
  VecEnv1=Zero
  VecEnv2=Zero
  CALL SetPosE2Pos
  Time= StartTime
  IF (.NOT.ALLOCATED(TimeEnvi)) THEN
    ALLOCATE(TimeEnvi(2))
    TimeEnvi(1)=StartTime
    TimeEnvi(2)=EndTime
  END IF
  Time1=TimeEnvi(1)
  Time2=TimeEnvi(2)
  SELECT CASE(Environment)
! CALL InputWRFData(InputFileName)
! CALL ReadWRF(VecEnv1)
    CASE('WRF') 
      CALL InputWRFnc(InputFileName)
      Time1=StartTime
      Time2=Time1+TimeIncr
      CALL ReadWRFnc(VecEnv1,InputFileNameWRF,Time1,VecMet)
      CALL ExchangeCell(VecMet)
      CALL VelocityCellToFaceLR(VecMet,VelF1,VelF1,Time)
      CALL ReadWRFnc(VecEnv2,InputFileNameWRF,Time2)
    CASE('Init')
      IF ( uPosEnv>0) THEN
        CALL VectorInit(uPosEnv,VecEnv1,UStartE,Time1)
        PosE2Pos(uPosEnv)=uPosL
      END IF
      IF ( vPosEnv>0) THEN
        CALL VectorInit(vPosEnv,VecEnv1,vStartE,Time1)
        PosE2Pos(vPosEnv)=vPosL
      END IF
      IF ( wPosEnv>0) THEN
        CALL VectorInit(wPosEnv,VecEnv1,wStart,Time1)
        PosE2Pos(wPosEnv)=wPosL
      END IF
      IF (thPosEnv>0) THEN
        CALL VectorInit(thPosEnv,VecEnv1,ThStart,Time1)
        PosE2Pos(thPosEnv)=thPos
      END IF
      IF (tkePosEnv>0) THEN
        CALL VectorInit(tkePosEnv,VecEnv1,tkeStart,Time1)
        PosE2Pos(tkePosEnv)=tkePos
      END IF
      IF (disPosEnv>0) THEN
        CALL VectorInit(disPosEnv,VecEnv1,disStart,Time1)
        PosE2Pos(disPosEnv)=disPos
      END IF
      IF (RhoVPosEnv>0) THEN
        CALL VectorInit(RhoVPosEnv,VecEnv1,QvStart,Time1)
        PosE2Pos(RhoVPosEnv)=RhoVPos
      END IF
      IF ( uPosEnv>0) THEN
        CALL VectorInit(uPosEnv,VecEnv2,UStartE,Time2)
      END IF
      IF ( vPosEnv>0) THEN
        CALL VectorInit(vPosEnv,VecEnv2,vStartE,Time2)
      END IF
      IF ( wPosEnv>0) THEN
        CALL VectorInit(wPosEnv,VecEnv2,wStart,Time2)
      END IF
      IF (thPosEnv>0) THEN
        CALL VectorInit(thPosEnv,VecEnv2,ThStart,Time2)
      END IF
      IF (tkePosEnv>0) THEN
        CALL VectorInit(tkePosEnv,VecEnv2,tkeStart,Time2)
      END IF
      IF (disPosEnv>0) THEN
        CALL VectorInit(disPosEnv,VecEnv2,disStart,Time2)
      END IF
      IF (RhoVPosEnv>0) THEN
        CALL VectorInit(RhoVPosEnv,VecEnv2,QvStart,Time2)
      END IF
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        IF (uPosEnv>0) THEN
          VecEnv1(ibLoc)%Vec(uPosEnv)%c=VecEnv1(ibLoc)%Vec(uPosEnv)%c*VecMet(ibLoc)%Vec(RhoPos)%c
          VecEnv2(ibLoc)%Vec(uPosEnv)%c=VecEnv2(ibLoc)%Vec(uPosEnv)%c*VecMet(ibLoc)%Vec(RhoPos)%c
        END IF  
        IF (vPosEnv>0) THEN
          VecEnv1(ibLoc)%Vec(vPosEnv)%c=VecEnv1(ibLoc)%Vec(vPosEnv)%c*VecMet(ibLoc)%Vec(RhoPos)%c
          VecEnv2(ibLoc)%Vec(vPosEnv)%c=VecEnv2(ibLoc)%Vec(vPosEnv)%c*VecMet(ibLoc)%Vec(RhoPos)%c
        END IF  
        IF (wPosEnv>0) THEN
          VecEnv1(ibLoc)%Vec(wPosEnv)%c=VecEnv1(ibLoc)%Vec(wPosEnv)%c*VecMet(ibLoc)%Vec(RhoPos)%c
          VecEnv2(ibLoc)%Vec(wPosEnv)%c=VecEnv2(ibLoc)%Vec(wPosEnv)%c*VecMet(ibLoc)%Vec(RhoPos)%c
        END IF  
      END DO  
  END SELECT    
  CALL ExchangeCell(VecEnv1)
  CALL ExchangeCell(VecEnv2)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL DampingInit
  
  IF (Diffusion) THEN
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      CALL LenScaleGeomCompute
    END DO
    IF (TkeDis.OR.TkeDisRich) THEN
      CALL Allocate     (DiffKoeff)
      CALL ScalarInit(DiffKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffKoeff)
      CALL ExchangeCell (DiffKoeff)
    ELSE IF (TkeSGS.OR.TkeLen.OR.NoTke) THEN
      CALL Allocate     (DiffPotKoeff)
      CALL ScalarInit(DiffPotKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffPotKoeff)
      CALL ExchangeCell (DiffPotKoeff)
      CALL Allocate     (DiffMomKoeff)
      CALL ScalarInit(DiffMomKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffMomKoeff)
      CALL ExchangeCell (DiffMomKoeff)
      CALL Allocate     (LenKoeff)
      CALL ScalarInit(LenKoeff,LenStart,Time)
      CALL Mult(RhoCell,LenKoeff)
      CALL ExchangeCell (LenKoeff)
    ELSE IF (TkeSmag) THEN
      CALL Allocate     (DiffPotHKoeff)
      CALL Allocate     (DiffPotVKoeff)
      CALL ScalarInit(DiffPotHKoeff,DStart,Time)
      CALL ScalarInit(DiffPotVKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffPotHKoeff)
      CALL Mult(RhoCell,DiffPotVKoeff)
      CALL ExchangeCell (DiffPotHKoeff)
      CALL ExchangeCell (DiffPotVKoeff)
      CALL Allocate     (DiffMomHKoeff)
      CALL Allocate     (DiffMomVKoeff)
      CALL ScalarInit(DiffMomHKoeff,DStart,Time)
      CALL ScalarInit(DiffMomVKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffMomHKoeff)
      CALL Mult(RhoCell,DiffMomVKoeff)
      CALL ExchangeCell (DiffMomHKoeff)
      CALL ExchangeCell (DiffMomVKoeff)
      CALL Allocate     (LenKoeff)
      CALL ScalarInit(LenKoeff,LenStart,Time)
      CALL Mult(RhoCell,LenKoeff)
      CALL ExchangeCell (LenKoeff)
    ELSE IF (DynSmag) THEN
      CALL Allocate     (DiffPotHKoeff)
      CALL Allocate     (DiffPotVKoeff)
      CALL ScalarInit(DiffPotHKoeff,DStart,Time)
      CALL ScalarInit(DiffPotVKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffPotHKoeff)
      CALL Mult(RhoCell,DiffPotVKoeff)
      CALL ExchangeCell (DiffPotHKoeff)
      CALL ExchangeCell (DiffPotVKoeff)
      CALL Allocate     (DiffMomHKoeff)
      CALL Allocate     (DiffMomVKoeff)
      CALL ScalarInit(DiffMomHKoeff,DStart,Time)
      CALL ScalarInit(DiffMomVKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffMomHKoeff)
      CALL Mult(RhoCell,DiffMomVKoeff)
      CALL ExchangeCell (DiffMomHKoeff)
      CALL ExchangeCell (DiffMomVKoeff)
      CALL Allocate     (LenKoeff)
      CALL ScalarInit(LenKoeff,LenStart,Time)
      CALL Mult(RhoCell,LenKoeff)
      CALL ExchangeCell (LenKoeff)
    ELSE IF (TkeHVLen) THEN
      CALL Allocate     (DiffHKoeff)
      CALL ScalarInit(DiffHKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffHKoeff)
      CALL ExchangeCell (DiffHKoeff)
      CALL Allocate     (DiffVKoeff)
      CALL ScalarInit(DiffVKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffVKoeff)
      CALL ExchangeCell (DiffVKoeff)
    ELSE
      CALL Allocate     (DiffKoeff)
      CALL ScalarInit(DiffKoeff,DStart,Time)
      CALL Mult(RhoCell,DiffKoeff)
      CALL ExchangeCell (DiffKoeff)
    END IF
  END IF
   
  IF (ForcingExtern) THEN
    DO iForce=1,72+ForcingShift
      FinishTendency=.FALSE.
      INQUIRE(FILE=FileForceTendency,EXIST=FileExist)
      IF (FileExist) THEN
        CALL InputTendencyProfile(Time,FileForceTendency,FinishTendency)
      ELSE  
        ForcingExternTendency=.FALSE.
      END IF
      ForcingExternNudging=.TRUE.
      FinishNudging=.FALSE.
      CALL InputNudgingProfile(Time,FileForceNudging,FinishNudging)
      CALL DampProfileCompute
!     INQUIRE(FILE=FileForceSurface,EXIST=FileExist)
!     IF (FileExist) THEN
!       CALL InputSurface(Time,FileForceSurface)
!     ELSE  
!       ForcingExternSurface=.FALSE.
!     END IF
      IF (FinishNudging) THEN
        IF (MyId==0) WRITE(*,*) 'FinishNudging',FinishNudging,Time,NudgeTime1,NudgeTime2
        IF (MyId==0) WRITE(*,*) 'FinishTendency',FinishTendency,Time,TendTime1,TendTime2
        EXIT
      END IF
    END DO
    Time=Time+TendTime1
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        VecMet(ibLoc)%Vec(uPosL)%c(:,:,iz,1)=NudgingProfile1(ibLoc)%Vec(1)%c(iz)*RhoProfile1(ibLoc)%Vec(1)%c(iz)
        VecMet(ibLoc)%Vec(uPosR)%c(:,:,iz,1)=NudgingProfile1(ibLoc)%Vec(1)%c(iz)*RhoProfile1(ibLoc)%Vec(1)%c(iz)
        VecMet(ibLoc)%Vec(vPosL)%c(:,:,iz,1)=NudgingProfile1(ibLoc)%Vec(2)%c(iz)*RhoProfile1(ibLoc)%Vec(1)%c(iz)
        VecMet(ibLoc)%Vec(vPosR)%c(:,:,iz,1)=NudgingProfile1(ibLoc)%Vec(2)%c(iz)*RhoProfile1(ibLoc)%Vec(1)%c(iz)
        VecMet(ibLoc)%Vec(ThPos)%c(:,:,iz,1)=NudgingProfile1(ibLoc)%Vec(3)%c(iz)*RhoProfile1(ibLoc)%Vec(1)%c(iz)
        VecMet(ibLoc)%Vec(RhoVPos)%c(:,:,iz,1)=NudgingProfile1(ibLoc)%Vec(4)%c(iz)*RhoProfile1(ibLoc)%Vec(1)%c(iz)
        VecMet(ibLoc)%Vec(RhoPos)%c(:,:,iz,1)=RhoProfile1(ibLoc)%Vec(1)%c(iz)
      END DO  
    END DO
    CALL VelocityCellToFaceLR(VecMet,VelF1,VelF1,Time)
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO iz=iz0+1,iz1
        IF (zP(iz)<400.0d0) THEN
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1
              CALL Random_Number(rRand)
              VecMet(ibLoc)%Vec(ThPos)%c(ix,iy,iz,1)=VecMet(ibLoc)%Vec(ThPos)%c(ix,iy,iz,1)+1.0d-2*rRand
            END DO  
          END DO
        END IF  
      END DO  
    END DO
    CALL ExchangeCell(VecMet)
  END IF  
! Initialize BoundCells
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib)) 
    DO i=1,NumBoundCell
      CALL SetBoundCells(BoundCell(i))
    END DO
  END DO

  ! -- Project VelF1
  CALL PrepareF(VecMet,VelF1,Time)
  CALL IniJacAccGrav
  IF (Chemie) THEN
    IF (IniFile/='') THEN 
      CALL InitGas(VecChem,RhoCell,IniFile)
      CALL Allocate(VecAmb,VecChem)
      VecAmb=Zero
      CALL InitAmbientGas(VecAmb,IniFile)
      CALL InitBoundaryGas(IniFile)
    END IF
  END IF

  IF (Aerosol) THEN
    CALL InitGrid
    CALL InitAerosol(VecChem,IniFile)
    CALL InitAmbientAero(VecAmb,IniFile)
    CALL Limiter(VecChem,VecChem)
    CALL InitKoagulation
    CALL InitActivity
  END IF

  IF (Chemie.OR.Aerosol) THEN
    IF (Emiss) THEN
      CALL SetEmission(IniFile)
      CALL SetPointEmission(IniFile)
    END IF
    IF (EmissStreet) THEN 
      CALL SetStreetEmission(InputFileName)
    END IF
  END IF
  IF (ThetaKind=='PreEn') THEN
    CALL PrepareEn(VecMet,VelF1,Time)
  END IF  
  IF (dtP>0.0d0.AND.JacSound) THEN
    CALL AllocateVec4Chemie(VecMetP,VectorComponentsT)
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DEALLOCATE(VecMetP(ibLoc)%Vec(thPos)%cB)
      ALLOCATE(VecMetP(ibLoc)%Vec(thPos)%cB(1:NumBoundCell,1:1))
    END DO
    CALL Copy(VecMet,VecMetP)
    GravComp=0.0d0
    CALL JacAccGrav(VecMetP)
    VecMetP=Zero
    IF (Anelastic.OR.PseudoIn) THEN
      CALL ProjectVelface(dtP,VelF1,VecMetP,VecMet,VecG=VecG)
    ELSE  
      CALL ProjectVelface(dtP,VelF1,VecMetP,VecMet)
    END IF  
    CALL Deallocate(VecMetP)
    IF (ThetaKind=='PreEn') THEN
      CALL PrepareEn(VecMet,VelF1,Time)
    END IF  
  END IF  
  GravComp=Grav
  Time=StartTime
  CALL PrepareF(VecMet,VelF1,Time)

  !Time=0.0e0 ??

  CALL InputErrorControl(InputFileName)
  IF (Position('aNUMBER')>0) THEN
    ATol(Position('aNUMBER'))=1.0d10*aTolG
  END IF
  CALL ExchangeCell(VecMet)

  CALL MPI_Barrier(MPI_Comm_World,MPIerr)

  CALL InputModelOutput(InputFileName)
  CALL SpecialOutput(VecMet,VecG,VelF1,StartTime)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL Output(StartTime,VelF1,VecMet,VecMet)

  dtAct=dtStart
  dt=dtStart
  
! Output Profile
  IF (ProfOut) THEN
    IF (uPosl>0) THEN
      CALL VelocityFaceToCellLR(VelF1,VecMet)
    END IF  
    CALL OutputProfile(VecMet,VecMet,ProfileSTART)
  END IF

! Mean Profiles
  IF (BCVel%West  =='MeanFlow'.OR.BCVel%East =='MeanFlow'.OR. &
      BCVel%South =='MeanFlow'.OR.BCVel%North=='MeanFlow'.OR. &
      BCVel%Bottom=='MeanFlow'.OR.BCVel%Top  =='MeanFlow') THEN
    CALL VelocityFaceToCellLR(VelF1,VecMet)
    CALL MeanProfileCompute(VecMet)
  END IF

  IF (Anelastic) THEN
    VecG=Zero
  END IF  
  IF (Method(3:6)=='Peer') THEN
    ALLOCATE(VecVecMet(1:3))
    VecVecMet(3)%Vec=>VecMet
    DO i=1,2
      CALL Allocate(VecVecMet(i)%Vec,0,VectorComponentsM)
      CALL Copy(VecMet,VecVecMet(i)%Vec)
    END DO
    ALLOCATE(VecVeLF(1:3))
    VecVelF(3)%VecF=>VelF1
    DO i=1,2
      CALL Allocate(VecVelF(i)%VecF)
      CALL Copy(VeLF1,VecVelF(i)%VecF)
    END DO
  END IF

! Write SpeciesName
  IF (MyId==0) THEN
    WRITE(*,*) ' .chem-file ::  i  ,  SpeciesName  ,  SIZE(VecMet(1)%Vec(i)%c,4)  ,  Position(SpeciesName(i))'
    DO i=1,SIZE(SpeciesName)
      WRITE(*,*) i,SpeciesName(i),SIZE(VecMet(1)%Vec(i)%c,4),Position(SpeciesName(i))
    END DO
  END IF  

  Time=StartTime
  !OutputTime=StartTime !!ML.

  CheckMean=.FALSE.
  CALL MeanProfileNMLOutput(InputFileName,CheckMean) !!ML.
  IF (CheckMean) CALL MeanProfileDataTable(InputFileName) !!ML.

  CheckColumnX=.FALSE.
  CALL ColumnXProfileNMLOutput(InputFileName,CheckColumnX) !!ML.
  IF (CheckColumnX) THEN
    CALL ColumnXProfileDataTable(InputFileName) !!ML.
  END IF  
  CheckColumnZ=.FALSE.
  CALL ColumnZProfileNMLOutput(InputFileName,CheckColumnZ) !!ML.
  IF (CheckColumnZ) THEN
    CALL ColumnZProfileDataTable(InputFileName) !!ML.
  END IF  
  CheckSurfMean=.FALSE.
  CALL MeanSurfaceNMLOutput(InputFileName,CheckSurfMean)
  IF (CheckSurfMean) CALL MeanSurfaceDataTable(InputFileName)
  CheckSurfPoint=.FALSE.
  CALL PointSurfaceNMLOutput(InputFileName,CheckSurfPoint)
  IF (CheckSurfPoint) CALL PointSurfaceDataTable(InputFileName)

  IF (CheckMean) CALL OutputMeanProfile(VecMet,Time) !!ML.
  IF (CheckColumnX) THEN
    CALL OutputColumnXProfile(VecMet,Time) !!ML.
  END IF  
  IF (CheckColumnZ) THEN
    CALL OutputColumnZProfile(VecMet,Time) !!ML.
  END IF    
  IF (CheckSurfMean) CALL OutputMeanSurface(VecMet,Time)
  IF (CheckSurfPoint) CALL OutputPointSurface(VecMet,Time)

  CALL ReadRestart(VecMet,VelF1,Time,InputFileName)
  CALL AllocateVec4Chemie(VecMetP,VectorComponentsT)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DEALLOCATE(VecMetP(ibLoc)%Vec(thPos)%cB)
    ALLOCATE(VecMetP(ibLoc)%Vec(thPos)%cB(1:NumBoundCell,1:SIZE(VecMet(ibLoc)%Vec(thPos)%cB,2)))
  END DO
  CALL Copy(VecMet,VecMetP)

  IF (MyID==0) WRITE(*,*) 'Starting time integration loop, ',Method
  DO iInt=1,EndIter
    IF (MyID==0) THEN
      CALL CPU_TIME(Start)
    END IF
    IF (ForcingExtern.AND.Time>=NudgeTime2) THEN
      INQUIRE(FILE=FileForceTendency,EXIST=FileExist)
      FinishTendency=.FALSE.
      IF (FileExist) THEN
        CALL InputTendencyProfile(Time,FileForceTendency,FinishTendency)
      END IF
      FinishNudging=.FALSE.
      CALL InputNudgingProfile(Time,FileForceNudging,FinishNudging)
!     INQUIRE(FILE=FileForceSurface,EXIST=FileExist)
!     IF (FileExist) THEN
!       CALL InputSurface(Time,FileForceSurface)
!     END IF
    END IF
    TimeAct=Time
    IF (MyId==0) THEN
      WRITE(*,*);  WRITE(*,*) ' Number of Steps  ::  ',iInt, ' Time  ::  ', Time
    END IF
    IF (Wind) THEN
      CALL UpdateWind(VelF1,VecMet,dtAct)
    END IF
    CALL UpdateBoundary(Time) 
    IF ((ThetaKind=='PreEn'.OR.ThetaKind=='Exner').AND.EnPos>0.AND.PressureUpdate) THEN
      CALL PressureEnergyCompute(VecMet,VelF1) 
    END IF
    IF (dtAct<Time*1.d-12) THEN 
      IF (MyId==0) WRITE(*,*) 'dtAct/Time < 1.d-12' 
    END IF 
    dtAct=MIN(dtAct,EndTime-Time) 
    IF (Method(1:3)=='Ros') THEN
      IF (MyId==0.AND.PrintNameLists) THEN
        WRITE(*,*) 'dtAct',dtAct
      END IF
      CALL Ros3MetC(VecMet,VecChem,VelF1,VecG,dtAct,Time,ATol,RTol)
      Temp=TotalScalar(VecMet,'RHO')
      IF (MyId==0) THEN
        WRITE(*,*) 'TotalRho',Temp
      END IF  
    ELSE IF (Method(3:4)=='RK'.OR.Method(3:4)=='EB'.OR.Method(3:5)=='MIS') THEN
      CALL InitExpIntRk(VecMet)
      CALL TimeStepCFL(VelF1,dtAct,ns)
      IF (MyId==0.AND.PrintNameLists) THEN
        WRITE(*,*) 'dtAct',dtAct,Time
      END IF
!     CALL ExpIntRk(VelF1,VecMet,dtAct,Time,ATol,RTol)
    ELSE IF (Method(3:6)=='Peer') THEN
      IF (MyId==0.AND.PrintNameLists) THEN
        WRITE(*,*) 'Peer dtAct',dtAct,Time
      END IF
!     CALL InitExpIntPeer(VecVecMet(1)%Vec)
      CALL TimeStepCFL(VelF1,dtAct,ns)
!     CALL ExpIntPeer(VecVelF,VecVecMet,dtAct,Time,ATol,RTol)
    ELSE IF (Method(4:8)=='LPeer') THEN
!     CALL IntLinPeer(VelF1,VecMet,dtAct,Time,ATol,RTol)
    ELSE
      IF (MyId==0) WRITE(*,*) 'Falsche Methode'
      EXIT
    END IF
    IF ( Time>EndTime*(1.d0-1.d-12))THEN 
      Time=EndTime+1.d-12
      CALL PrepareF(VecMet,VelF1,Time)
      CALL SpecialOutput(VecMet,VecG,VelF1,Time)
      CALL Output(Time,VelF1,VecMet,VecMet)
      IF (CheckMean) CALL OutputMeanProfile(VecMet,Time) !ML.
      IF (CheckColumnX) THEN
        CALL OutputColumnXProfile(VecMet,Time) !ML.
      END IF  
      IF (CheckColumnZ) THEN
        CALL OutputColumnZProfile(VecMet,Time) !ML.
      END IF  
      IF (CheckSurfMean) CALL OutputMeanSurface(VecMet,Time)
      IF (CheckSurfPoint) CALL OutputPointSurface(VecMet,Time)
      EXIT
    END IF
    CALL PrepareF(VecMet,VelF1,Time)
    CALL SpecialOutput(VecMet,VecG,VelF1,Time)
    CALL Output(Time,VelF1,VecMet,VecMet)
    IF (CheckMean) CALL OutputMeanProfile(VecMet,Time)!ML.
    IF (CheckColumnX) CALL OutputColumnXProfile(VecMet,Time)!ML.
    IF (CheckColumnZ) THEN
      CALL OutputColumnZProfile(VecMet,Time)!ML.
    END IF
    IF (CheckSurfMean) CALL OutputMeanSurface(VecMet,Time)
    IF (CheckSurfPoint) CALL OutputPointSurface(VecMet,Time)
!   Mean Profiles
    IF (BCVel%West  =='MeanFlow'.OR.BCVel%East =='MeanFlow'.OR. &
        BCVel%South =='MeanFlow'.OR.BCVel%North=='MeanFlow'.OR. &
        BCVel%Bottom=='MeanFlow'.OR.BCVel%Top  =='MeanFlow') THEN
      CALL VelocityFaceToCellLR(VelF1,VecMet)
      CALL MeanProfileCompute(VecMet)
    END IF
    CALL WriteRestart(VecMet,VelF1,Time,InputFileName)
!   Check for NaN
    IF (ISNAN(DOT(VecMet,VecMet))) THEN
      CALL MPI_Finalize(MPIErr)
      IF (MyID==0) WRITE(*,*) 'NaN. Stopped'
      STOP
    END IF
    IF (MyID==0) THEN
      CALL CPU_TIME(Finish)
      WRITE(*,*) 'Passed time for one time step: ',finish-start,' s' 
      WRITE(*,*) 'Real time factor: ',(finish-start)/dtAct
    END IF
  END DO

! Output Profile
  IF (ProfOut) THEN
    CALL VelocityFaceToCellLR(VelF1,VecMet)
    CALL OutputProfile(VecMet,VecMet,ProfileEND)
  END IF
  CALL MPI_Finalize(MPIErr)

CONTAINS

SUBROUTINE MaxVelF(VelF)
  TYPE(VelocityFace_T), POINTER :: VelF(:)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1-1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          IF (ABS(VelF(ibLoc)%uF(ix,iy,iz))>15.0d0) THEN
            WRITE(*,*) ix,iy,iz,'VelF',VelF(ibLoc)%uF(ix,iy,iz)
            WRITE(*,*) 'Zelle links'
            WRITE(*,*) VolC(ix,iy,iz)
            WRITE(*,*) FU(ix-1,iy,iz),FU(ix,iy,iz)
            WRITE(*,*) FV(ix,iy-1,iz),FV(ix,iy,iz)
            WRITE(*,*) FW(ix,iy,iz-1),FW(ix,iy,iz)
            WRITE(*,*) 'Zelle rechts'
            WRITE(*,*) VolC(ix+1,iy,iz)
            WRITE(*,*) FU(ix+1-1,iy,iz),FU(ix+1,iy,iz)
            WRITE(*,*) FV(ix+1,iy-1,iz),FV(ix+1,iy,iz)
            WRITE(*,*) FW(ix+1,iy,iz-1),FW(ix+1,iy,iz)
          END IF  
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE

FUNCTION ThetaLoc(z)

  REAL(RealKind) :: ThetaLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: S,N=1.d-2,th0=300.0d0

  S=N*N/Grav
  ThetaLoc=th0*EXP(S*z)

END FUNCTION ThetaLoc

FUNCTION PressLoc(z)

  REAL(RealKind) :: PressLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: S,N=1.0d-2,th0=300.0d0

  S=N*N/Grav
  IF (N>Zero) THEN
    PressLoc=p0*(One-Grav/(Cpd*th0*S)*(One-EXP(-S*z)))**(Cpd/Rd)
  ELSE
    PressLoc=p0*(One-kappa*Grav*z/(Rd*th0))**(Cpd/Rd)
  END IF


END FUNCTION PressLoc

SUBROUTINE Hydrostatic

  INTEGER, PARAMETER :: nzz=10
  REAL(RealKind), PARAMETER :: H=10.0d3
  INTEGER :: i
  REAL(RealKind) :: pLoc 
  REAL(RealKind) :: zz(nzz),ThVert(nzz),RhoVert(nzz)
  REAL(RealKind) :: dzz

  dzz=H/nzz
  zz(1)=0.5d0*dzz
  DO i=2,nzz
    zz(i)=zz(i-1)+dzz
  END DO
  DO i=1,nzz
    ThVert(i)=ThetaLoc(zz(i))
  END DO
  pLoc=PressLoc(zz(1))
  RhoVert(1)=pLoc/((pLoc/p0)**kappa*Rd*ThVert(1))
  DO i=2,nzz
    Call RhoHydro(Rhovert(i),ThVert(i),zz(i) &
                 ,Rhovert(i-1),ThVert(i-1),zz(i-1))
  END DO
  OPEN(UNIT=10,FILE='Profile',STATUS='UNKNOWN')
  WRITE(10,*) 'RhoProf'
  WRITE(10,*) nzz,2 
  WRITE(10,*) 'Equal' 
  DO i=1,nzz
    WRITE(10,*) zz(i),Rhovert(i)
  END DO
  WRITE(10,*) 'thProf'
  WRITE(10,*) nzz,2 
  WRITE(10,*) 'Equal' 
  DO i=1,nzz
    WRITE(10,*) zz(i),Thvert(i)
  END DO
  CLOSE(10)
  
END SUBROUTINE Hydrostatic

SUBROUTINE RhoHydro(Rho,Th,z,Rho1,Th1,z1)

  REAL(RealKind) :: Rho,Th,z,Rho1,Th1,z1

  REAL(RealKind) :: p1,p
  REAL(RealKind) :: RhoNew,RhoL,RhoR
  REAL(RealKind) :: Rhs
  REAL(RealKind) :: F,FL,FR,DF
  REAL(RealKind) :: TolErr

  IF (RealKind==8) THEN
    TolErr=1.e-9_RealKind 
  ELSE  
    TolErr=1.e-4_RealKind 
  END IF  
  p1=PresLoc(Rho1,Th1)
  Rhs=-p1+0.5d0*(z-z1)*Grav*Rho1
  RhoL=0.0d0
  p=PresLoc(RhoL,Th)
  FL=p+0.5d0*(z-z1)*Grav*RhoL+Rhs
  RhoR=Rho1
  p=PresLoc(RhoR,Th)
  FR=p+0.5d0*(z-z1)*Grav*RhoR+Rhs
  F=FR
  Rho=RhoR
  DO 
    DF=DPreDRho(Rho,Th)+0.5d0*(z-z1)*Grav
    RhoNew=Rho-F/DF
    IF (RhoNew<RhoL.OR.RhoNew>RhoR) THEN
      RhoNew=0.5d0*(RhoL+RhoR)
    END IF
    p=PresLoc(RhoNew,Th)
    F=p+0.5d0*(z-z1)*Grav*RhoNew+Rhs
    IF (F<0.0d0) THEN
      FL=F
      RhoL=RhoNew
    ELSE
      FR=F
      RhoR=RhoNew
    END IF
    IF (ABS(F)<=TolErr.OR.ABS(Rho-RhoNew)<=TolErr*Rho) THEN
      EXIT
    ELSE
      Rho=RhoNew
    END IF
  END DO

END SUBROUTINE RhoHydro

FUNCTION PresLoc(RhoLoc,ThLoc)
  REAL(RealKind) :: PresLoc
  REAL(RealKind) :: RhoLoc,ThLoc
  PresLoc=p0*(Rd*RhoLoc*ThLoc/p0)**(One/(One-kappa))
END FUNCTION PresLoc

FUNCTION DPreDRho(RhoLoc,ThLoc)
  REAL(RealKind) :: DPreDrho
  REAL(RealKind) :: RhoLoc,ThLoc
  DPreDrho=p0*(One/(One-kappa))*(Rd*RhoLoc*ThLoc/p0)**(kappa/(One-kappa)) &
           *Rd*ThLoc/p0
END FUNCTION DPreDrho

SUBROUTINE UpdateBoundary(Time)  

  REAL(RealKind) :: Time
  INTEGER, SAVE :: PosTimeEnvi=3
  IF (Time>Time2) THEN
    Time1=Time2
    VecEnv2=Zero
    IF ( uPosEnv>0) THEN
      CALL VectorInit(uPosEnv,VecEnv1,UStartE,Time1)
    END IF
    IF ( vPosEnv>0) THEN
      CALL VectorInit(vPosEnv,VecEnv1,vStartE,Time1)
    END IF
    IF ( wPosEnv>0) THEN
      CALL VectorInit(wPosEnv,VecEnv1,wStart,Time1)
    END IF
    IF (thPosEnv>0) THEN
      CALL VectorInit(thPosEnv,VecEnv1,ThStart,Time1)
    END IF
    IF (tkePosEnv>0) THEN
      CALL VectorInit(tkePosEnv,VecEnv1,tkeStart,Time1)
    END IF
    IF (disPosEnv>0) THEN
      CALL VectorInit(disPosEnv,VecEnv1,disStart,Time1)
    END IF
    Time2=TimeEnvi(PosTimeEnvi)
    PosTimeEnvi = PosTimeEnvi + 1
    VecEnv2=Zero
    IF ( uPosEnv>0) THEN
      CALL VectorInit(uPosEnv,VecEnv2,UStartE,Time2)
    END IF
    IF ( vPosEnv>0) THEN
      CALL VectorInit(vPosEnv,VecEnv2,vStartE,Time2)
    END IF
    IF ( wPosEnv>0) THEN
      CALL VectorInit(wPosEnv,VecEnv2,wStart,Time2)
    END IF
    IF (thPosEnv>0) THEN
      CALL VectorInit(thPosEnv,VecEnv2,ThStart,Time2)
    END IF
    IF (tkePosEnv>0) THEN
      CALL VectorInit(tkePosEnv,VecEnv2,tkeStart,Time2)
    END IF
    IF (disPosEnv>0) THEN
      CALL VectorInit(disPosEnv,VecEnv2,disStart,Time2)
    END IF
  END IF
END SUBROUTINE UpdateBoundary

END PROGRAM MainProg
