PROGRAM MainProg

  USE Init_Mod
  USE JacAccGrav_Mod
  USE Int_Mod
  USE IntPeer_Mod
  USE IntLinPeer_Mod
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
  USE Restart_Mod
  USE Grid_Mod

  IMPLICIT NONE

  TYPE(VelocityFace_T), POINTER :: VelF1(:)
  TYPE(VecVelocityFace_T), POINTER :: VecVelF(:)
  TYPE(Vector4Cell_T), POINTER :: VecT(:)
  TYPE(Vector4Cell_T), POINTER :: VecG(:)
  TYPE(VecVector4Cell_T), POINTER :: VecVecT(:)
  TYPE(Vector4Cell_T), POINTER :: VecTP(:)
  INTEGER :: ib,ibLoc
  REAL(RealKind) :: Temp,RhoLoc
  INTEGER :: Iter
  INTEGER :: iInt,ix,iy,iz,ic
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
  REAL(RealKind) :: dtAct,Time
  REAL(RealKind) :: dtCFL=0.0d0
  INTEGER :: nsCFL

  REAL(RealKind) :: TotalStart(1:5),TotalCurrent(1:5)
  INTEGER :: SizeSeed=12
  INTEGER,Dimension(12) :: seed
  CHARACTER :: ArgRestart

  REAL(RealKind) :: ErrorL1,SumVol


  CALL start_MPI

  CALL Random_seed(size=SizeSeed)
! seed = 1 !myID
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
  WRITE(*,*)
  WRITE(*,*)
  IF (MyId==0) THEN
    IF (InputFileName == ' ') THEN
!      STOP 'Dateiname fehlt'
    END IF
  END IF
  OutputNameSoil=TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))//'.Soil'

  CALL MPI_Bcast(InputFileName,80,MPI_CHARACTER,0,MPI_COMM_WORLD,MPIErr)
  IF (MyID==0) WRITE(*,*) InputFileName
  WRITE(*,*) 'MyID',MyID

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
    CALL InputSystem(ChemieFile)
  END IF

  CALL MPI_Barrier(MPI_Comm_World,MPIerr)
  CALL SetIndices
  CALL Allocate(JacMet)
  JacMet=Zero
  CALL AllocateVec4Chemie(VecT,VectorComponentsT)

  CALL InitPrandtlNumber(InputFileName)

! FallVelocity
  ALLOCATE(FallVelocity(0:VectorComponentsT))
  FallVelocity=Zero
  IF (RhorPos>0) THEN
    FallVelocity(RhorPos)=0.0d0
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

  Vect=Zero
  IF (Chemie.AND.DataFile/='') THEN 
!   -- Eingabe Chemie --
    CALL InputChemicalData(DataFile)
    CALL OutputChemie(ChemOutFile)
  END IF

! -- WindFarm--
  IF (Wind.AND.WindFarmFileName/='') THEN 
! -- Input WindFarm --
    CALL InitWind(WindFarmFileName)
  END IF  
   

! -- Meteorologie --
  CALL Allocate(VelF1)
  VelF1=Zero
  VecT=Zero

  Time=StartTime
  IF (Radiation.AND..NOT.RadiationProfile) THEN
    CALL Allocate(RaddirCell)
    CALL Allocate(RaddifCell)
    CALL Allocate(RadinfCell)
    CALL Allocate(ShadowCell)
  END IF
  IF (Shallow) THEN
    CALL Allocate(HeightG)
    CALL ScalarInit(HeightG,HeightFun,Time)
  END IF
  CALL Allocate(RhoCell)
  ALLOCATE(DummyCell(nbLoc))
  RhoCell=Zero
  CALL ScalarInit(RhoCell,RhoFun,Time)
  CALL ExchangeCell(RhoCell)
  IF (RhoPos>0) THEN
    CALL Allocate(RhoProfG)
    CALL ScalarInit(RhoProfG,RhoProf,Time)
    CALL ExchangeCell(RhoProfG)
    CALL VectorInit(RhoPos,VecT,RhoStart,Time)
    CALL Mult(RhoCell,VecT,RhoPos)
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
    CALL Allocate(VecG,1)
    CALL Allocate(RhoProfG)
    CALL ScalarInit(RhoProfG,RhoFun,Time)
    CALL ExchangeCell(RhoProfG)
    CALL Allocate(RhoCell)
    CALL ScalarInit(RhoCell,RhoFun,Time)
    CALL ExchangeCell(RhoCell)
    CALL Allocate(PreCell)
  END IF
  IF ( uPosL>0)  THEN
    CALL VectorInit(uPosL,VecT,UStart,Time)
    CALL Mult(RhoCell,VecT,uPosL)
    CALL Allocate(uCell)
  END IF
  IF ( uPosR>0)  THEN
    CALL VectorInit(uPosR,VecT,UStart,Time)
    CALL Mult(RhoCell,VecT,uPosR)
  END IF
  IF ( vPosL>0) THEN
    CALL VectorInit(vPosL,VecT,vStart,Time)
    CALL Mult(RhoCell,VecT,vPosL)
    CALL Allocate(vCell)
  END IF
  IF ( vPosR>0) THEN
    CALL VectorInit(vPosR,VecT,vStart,Time)
    CALL Mult(RhoCell,VecT,vPosR)
  END IF
  IF ( wPosL>0) THEN
    CALL VectorInit(wPosL,VecT,wStart,Time)
    CALL Mult(RhoCell,VecT,wPosL)
    CALL Allocate(wCell)
  END IF
  IF ( wPosR>0) THEN
    CALL VectorInit(wPosR,VecT,wStart,Time)
    CALL Mult(RhoCell,VecT,wPosR)
  END IF
  IF ( thPos>0) THEN
    CALL VectorInit(thPos,VecT,ThStart,Time)
!   CALL BoundaryInit('Skin',TStart,Time)
    CALL Mult(RhoCell,VecT,thPos)
    CALL Allocate(thProfG)
    CALL ScalarInit(thProfG,thProfFun,Time)
    CALL Mult(RhoProfG,thProfG)
    CALL ExchangeCell(thProfG)
    CALL Allocate(TAbsCell,1,1)
    TAbsCell=300.0e0_RealKind
    CALL Allocate(KinEnCell)
    CALL Allocate(ECell)
    KinEnCell=Zero
    IF (DynamicSoil) THEN
      CALL InitSoil
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DEALLOCATE(VecT(ibLoc)%Vec(thPos)%cB)
        ALLOCATE(VecT(ibLoc)%Vec(thPos)%cB(1:NumBoundCell,1:Domain%nrsoillayers)) 
        DEALLOCATE(VecT(ibLoc)%Vec(RhoVPos)%cB)
        ALLOCATE(VecT(ibLoc)%Vec(RhoVPos)%cB(1:NumBoundCell,1:Domain%nrsoillayers+1))
      END DO
      CALL VectorInitSoilFunction(thPos,VecT,ThStartSoil,Time)
    ELSE IF (Canopy) THEN
      CALL InitCanopy
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DEALLOCATE(VecT(ibLoc)%Vec(thPos)%cB)
        ALLOCATE(VecT(ibLoc)%Vec(thPos)%cB(1:NumBoundCell,1:Domain%nrsoillayers)) 
        DEALLOCATE(VecT(ibLoc)%Vec(RhoVPos)%cB)
        ALLOCATE(VecT(ibLoc)%Vec(RhoVPos)%cB(1:NumBoundCell,1:Domain%nrsoillayers+1))
      END DO
      CALL VectorInitCanopyFunction(thPos,VecT,ThStartSoil,Time)
    END IF
  END IF
  IF (EnPos>0) THEN
    CALL VectorInit(EnPos,VecT,EnStart,Time)
    CALL Mult(RhoCell,VecT,EnPos)
  END IF  
  IF ( RhoVPos>0) THEN
    CALL Allocate(QvProfG)
    CALL ScalarInit(QvProfG,qvStart,Time)
    CALL ExchangeCell(QvProfG)
    CALL Mult(RhoProfG,QvProfG)
    CALL VectorInit(RhoVPos,VecT,qvStart,Time)
    CALL Mult(RhoCell,VecT,RhoVPos)
    IF (DynamicSoil) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
      CALL VectorInitSoilFunction(RhoVPos,VecT,QvStartSoil,Time)
    ELSE IF (Canopy) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
      CALL VectorInitCanopyFunction(RhoVPos,VecT,QvStartSoil,Time)
    END IF
  END IF
  IF ( RhoCPos>0) THEN
    CALL Allocate(QcProfG)
    QcProfG=Zero
    CALL ExchangeCell(QcProfG)
    CALL Mult(RhoProfG,QcProfG)
    CALL VectorInit(RhoCPos,VecT,qcStart,Time)
    CALL Mult(RhoCell,VecT,RhoCPos)
  ELSE 
    CALL Allocate(RhoLCell)
  END IF
  IF ( RhorPos>0) THEN
    CALL VectorInit(RhorPos,VecT,qrStart,Time)
    CALL Mult(RhoCell,VecT,RhorPos)
    IF (PrecipIce.AND.IceSurf) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
    END IF
  END IF
  IF ( RhoIPos>0) THEN
    CALL VectorInit(RhoIPos,VecT,qiStart,Time)
    CALL Mult(RhoCell,VecT,RhoIPos)
    IF (PrecipIce.AND.IceSurf) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
      END DO
    END IF
  END IF
  IF ( RhoSPos>0) THEN
    CALL VectorInit(RhoSPos,VecT,qsStart,Time)
    CALL Mult(RhoCell,VecT,RhoSPos)
  END IF
  IF ( nvPos>0) THEN
    CALL VectorInit(nvPos,VecT,nvStart,Time)
  END IF
  IF ( ncPos>0) THEN
    CALL VectorInit(ncPos,VecT,ncStart,Time)
  END IF
  IF ( nrPos>0) THEN
    CALL VectorInit(nrPos,VecT,nrStart,Time)
  END IF
  IF ( niPos>0) THEN
    CALL VectorInit(niPos,VecT,niStart,Time)
  END IF
  IF ( nsPos>0) THEN
    CALL VectorInit(nsPos,VecT,nsStart,Time)
  END IF
  IF (tkePos>0) THEN
    CALL VectorInit(tkePos,VecT,TkeStart,Time)
    CALL Mult(RhoCell,VecT,tkePos)
  END IF
  IF (tkeHPos>0) THEN
    CALL VectorInit(tkeHPos,VecT,TkeHStart,Time)
    CALL Mult(RhoCell,VecT,tkeHPos)
  END IF
  IF (tkeVPos>0) THEN
    CALL VectorInit(tkeVPos,VecT,TkeVStart,Time)
    CALL Mult(RhoCell,VecT,tkeVPos)
  END IF
  IF (disPos>0) THEN
    CALL VectorInit(disPos,VecT,DisStart,Time)
    CALL Mult(RhoCell,VecT,disPos)
  END IF
  IF (omePos>0) THEN
    CALL VectorInit(omePos,VecT,OmeStart,Time)
    CALL Mult(RhoCell,VecT,omePos)
  END IF
  IF (LenPos>0) THEN
    CALL VectorInit(LenPos,VecT,LenStart,Time)
    CALL Mult(RhoCell,VecT,LenPos)
  END IF
  IF (Tracer1Pos>0) THEN
    CALL VectorInit(Tracer1Pos,VecT,Tracer1Start,Time)
  END IF
  IF (Tracer2Pos>0) THEN
    CALL VectorInit(Tracer2Pos,VecT,Tracer2Start,Time)
  END IF
  CALL PerturbProfile(VecT)
  CALL ExchangeCell(VecT)
  IF (uPosL*uPosR>0) THEN
    CALL VelocityCellToFaceLR(VecT,VelF1,VelF1,Time)
    CALL RhoVelocityInit(VelF1,UStart,VStart,WStart,RhoFun,Time)
    CALL BoundaryVelocity(VelF1,Time)
  ELSE
    CALL VelocityInit(VelF1,UStart,VStart,WStart,Time)
  END IF
  Time= StartTime
  IF (.NOT.ALLOCATED(TimeEnvi)) THEN
    ALLOCATE(TimeEnvi(2))
    TimeEnvi(1)=StartTime
    TimeEnvi(2)=EndTime
  END IF
  Time1=TimeEnvi(1)
  Time2=TimeEnvi(2)
  IF (Parcel) THEN
    CALL InputAdiabatic(InputFileName,VecT)
    CALL Allocate(VecAmb,VecT)
    VecAmb=Zero
  END IF
  CALL PrepareF(VecT,VelF1,Time)
  IF (Chemie) THEN
    IF (IniFile/='') THEN 
      CALL InitGas(VecT,RhoCell,IniFile)
      CALL InitAmbientGas(VecAmb,IniFile)
      PosDummy=Position('Dummy1')
      IF (PosDummy>0) THEN
        CALL VectorInit(PosDummy,VecT,DummyStart1,Time)
        CALL Mult(RhoCell,VecT,PosDummy)
      END IF
      PosDummy=Position('Dummy2')
      IF (PosDummy>0) THEN
        CALL VectorInit(PosDummy,VecT,DummyStart2,Time)
        CALL Mult(RhoCell,VecT,PosDummy)
      END IF
      PosDummy=Position('Dummy3')
      IF (PosDummy>0) THEN
        CALL VectorInit(PosDummy,VecT,DummyStart3,Time)
        CALL Mult(RhoCell,VecT,PosDummy)
      END IF
      PosDummy=Position('Dummy4')
      IF (PosDummy>0) THEN
        CALL VectorInit(PosDummy,VecT,DummyStart4,Time)
        CALL Mult(RhoCell,VecT,PosDummy)
      END IF
    ELSE
      DO i=NumMet+1,VectorComponentsT 
        CALL VectorInit(i,VecT,DummyStart,Time)
        CALL Mult(RhoCell,VecT,i)
      END DO
    END IF
  END IF
  IF (Aerosol) THEN
    CALL InitGrid
    CALL InitAerosol(VecT,IniFile)
    CALL InitAmbientAero(VecAmb,IniFile)
    CALL Limiter(VecT,VecT)
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
  IF (Baum) THEN 
    CALL ReadBaum(InputFileName) 
  END IF
! Environment
  VectorComponentsME=MIN(VectorComponentsME,VectorComponentsM)
  CALL Allocate(VecEnv1,VectorComponentsME)
  CALL Allocate(VecEnv2,VectorComponentsME)
  VecEnv1=Zero
  VecEnv2=Zero
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
  CALL ExchangeCell(VecEnv1)
  CALL ExchangeCell(VecEnv2)
  CALL MPI_Barrier(MPI_Comm_World,MPIerr)

  CALL DampingInit

  IF (Diffusion) THEN
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      CALL LenScaleGeomCompute(ib)
    END DO
    IF (TkeDis) THEN
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


! -- Project VelF1
  CALL IniJacAccGrav
  IF (ThetaKind=='PreEn') THEN
    CALL PrepareEn(VecT,VelF1,Time)
  END IF  
  IF (dtP>0.0d0) THEN
    CALL AllocateVec4Chemie(VecTP,VectorComponentsT)
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DEALLOCATE(VecTP(ibLoc)%Vec(thPos)%cB)
      ALLOCATE(VecTP(ibLoc)%Vec(thPos)%cB(1:NumBoundCell,1:1))
    END DO
    CALL Copy(VecT,VecTP)
    GravComp=0.0d0
    CALL JacAccGrav(VecTP)
    VecTP=Zero
    IF (Anelastic.OR.PseudoIn) THEN
      CALL ProjectVelface(dtP,VelF1,VecTP,VecT,VecG=VecG)
    ELSE  
      CALL ProjectVelface(dtP,VelF1,VecTP,VecT)
    END IF  
    CALL Deallocate(VecTP)
    IF (ThetaKind=='PreEn') THEN
      CALL PrepareEn(VecT,VelF1,Time)
    END IF  
  END IF  
 
  GravComp=Grav

  Time=StartTime
  CALL PrepareF(VecT,VelF1,Time)

  Time=0.0e0

  CALL InputErrorControl(InputFileName)
  IF (Position('aNUMBER')>0) THEN
    ATol(Position('aNUMBER'))=1.0d10*aTolG
  END IF

  CALL MPI_Barrier(MPI_Comm_World,MPIerr)

  CALL InputModelOutput(InputFileName)
  CALL InitOutputSpecial(InputFileName)
  CALL SpecialOutput(VecT,VecG,VelF1,StartTime)
  CALL Output(StartTime,VelF1,VecT,VecT)

  dtAct=dtStart
  dt=dtStart

! Output Profile
  IF (ProfOut) THEN
    IF (uPosl>0) THEN
      CALL VelocityFaceToCellLR(VelF1,VecT)
    END IF  
    CALL OutputProfile(VecT,VecT,ProfileSTART)
  END IF

! Mean Profiles
  IF (BCVel%West  =='MeanFlow'.OR.BCVel%East =='MeanFlow'.OR. &
      BCVel%South =='MeanFlow'.OR.BCVel%North=='MeanFlow'.OR. &
      BCVel%Bottom=='MeanFlow'.OR.BCVel%Top  =='MeanFlow') THEN
    CALL VelocityFaceToCellLR(VelF1,VecT)
    CALL MeanProfileCompute(VecT)
  END IF

  IF (Anelastic) THEN
    VecG=Zero
  END IF  
  IF (Method(3:6)=='Peer') THEN
    ALLOCATE(VecVecT(1:3))
    VecVecT(3)%Vec=>VecT
    DO i=1,2
      CALL Allocate(VecVecT(i)%Vec,0,VectorComponentsM)
      CALL Copy(VecT,VecVecT(i)%Vec)
    END DO
    ALLOCATE(VecVeLF(1:3))
    VecVelF(3)%VecF=>VelF1
    DO i=1,2
      CALL Allocate(VecVelF(i)%VecF)
      CALL Copy(VeLF1,VecVelF(i)%VecF)
    END DO
  END IF

! Initialize BoundCells

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib)) 
    DO i=1,NumBoundCell
      CALL SetBoundCells(BoundCell(i))
    END DO
  END DO
! Write SpeciesName
  IF (MyId==0) THEN
    DO i=1,SIZE(SpeciesName)
      WRITE(*,*) i,SpeciesName(i)
    END DO
  END IF  

  Time=StartTime

  CALL ReadRestart(VecT,VelF1,Time,InputFileName)
  CALL AllocateVec4Chemie(VecTP,VectorComponentsT)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DEALLOCATE(VecTP(ibLoc)%Vec(thPos)%cB)
    ALLOCATE(VecTP(ibLoc)%Vec(thPos)%cB(1:NumBoundCell,1:1))
  END DO
  CALL Copy(VecT,VecTP)
  DO iInt=1,EndIter
    TimeAct=Time
    IF (MyId==0) THEN
      WRITE(*,*)
      WRITE(*,*)            '                               ',iInt,Time
    END IF
    IF (Wind) THEN
      CALL UpdateWind(VelF1,VecT,dtAct)
    END IF
    CALL UpdateBoundary(Time) 
    IF ((ThetaKind=='PreEn'.OR.ThetaKind=='Exner').AND.EnPos>0.AND.PressureUpdate) THEN
      CALL PressureEnergyCompute(VecT,VelF1) 
    END IF
    IF (dtAct<Time*1.d-12) THEN 
      WRITE(*,*) 'dtAct/Time < 1.d-12' 
    END IF 
    dtAct=MIN(dtAct,EndTime-Time) 
    IF (Method(1:3)=='Ros') THEN
      IF (MyId==0) THEN
        WRITE(*,*) 'dtAct',dtAct
      END IF
      CALL Ros3MetC(VecT,VelF1,VecG,dtAct,Time,ATol,RTol)
    ELSE IF (Method(3:4)=='RK'.OR.Method(3:4)=='EB'.OR.Method(3:5)=='MIS') THEN
      IF (MyId==0) THEN
        WRITE(*,*) 'dtAct',dtAct,Time
      END IF
      CALL TimeStepCFL(VelF1,dtCFL,nsCFL)
      IF (MyId==0) THEN
        WRITE(*,*) 'dtCFL ',dtCFL
        WRITE(*,*) 'ns    ',ns
        WRITE(*,*) 'nsCFL ',nsCFL
      END IF  
      dtAct=MIN(dtAct,dtMax)
      CALL ExpIntRk(VelF1,VecT,dtAct,Time,ATol,RTol)
    ELSE IF (Method(3:6)=='Peer') THEN
      CALL ExpIntPeer(VecVelF,VecVecT,dtAct,Time,ATol,RTol)
    ELSE IF (Method(4:8)=='LPeer') THEN
      CALL IntLinPeer(VelF1,VecT,dtAct,Time,ATol,RTol)
    ELSE
      WRITE(*,*) 'Falsche Methode'
      EXIT
    END IF
    IF ( Time>EndTime*(1.d0-1.d-12))THEN 
      Time=EndTime+1.d-12
      CALL PrepareF(VecT,VelF1,Time)
      CALL SpecialOutput(VecT,VecG,VelF1,Time)
      CALL Output(Time,VelF1,VecT,VecT)
      EXIT
    END IF
    CALL PrepareF(VecT,VelF1,Time)
    CALL SpecialOutput(VecT,VecG,VelF1,Time)
    CALL Output(Time,VelF1,VecT,VecT)
!   Mean Profiles
    IF (BCVel%West  =='MeanFlow'.OR.BCVel%East =='MeanFlow'.OR. &
        BCVel%South =='MeanFlow'.OR.BCVel%North=='MeanFlow'.OR. &
        BCVel%Bottom=='MeanFlow'.OR.BCVel%Top  =='MeanFlow') THEN
      CALL VelocityFaceToCellLR(VelF1,VecT)
      CALL MeanProfileCompute(VecT)
    END IF
    CALL WriteRestart(VecT,VelF1,Time,InputFileName)
!   Check for NaN
    IF (ISNAN(DOT(VecT,VecT))) THEN
      CALL MPI_Finalize(MPIErr)
      WRITE(*,*) 'NaN. Stopped'
      STOP
    END IF
  END DO
  ErrorL1=0.0d0
  SumVol=0.0d0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          ErrorL1=ErrorL1+ABS(VecTP(ib)%Vec(1)%c(ix,iy,iz,1)-VecTP(ib)%Vec(1)%c(ix,iy,iz,1))*VolC(ix,iy,iz)
        END DO
      END DO
    END DO
  END DO
  ErrorL1=ErrorL1/SumVol
  WRITE(*,*) 'ErrorL1',ErrorL1

! Output Profile
  IF (ProfOut) THEN
    CALL VelocityFaceToCellLR(VelF1,VecT)
    CALL OutputProfile(VecT,VecT,ProfileEND)
  END IF
  CALL MPI_Finalize(MPIErr)

CONTAINS

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
