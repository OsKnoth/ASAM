MODULE Output_Mod

  USE Kind_Mod
  USE DataType_Mod
  USE Physics_Mod
  USE Control_Mod
  USE Parallel_Mod
  USE Chemie_Mod
  USE Aerosol_Mod
  USE Soil_Mod
  USE ScalarVectorCellPar_Mod

  IMPLICIT NONE


!  INTEGER :: ix0Out,ix1Out
!  INTEGER :: iy0Out,iy1Out
!  INTEGER :: iz0Out,iz1Out
  INTEGER :: NumData
  INTEGER :: NumDataB
  INTEGER :: NumDataC ! Cut
  INTEGER(KIND=MPI_OFFSET_KIND) :: Offset0=0
  INTEGER(KIND=MPI_OFFSET_KIND) :: OffsetSoil0=0
  INTEGER(KIND=MPI_OFFSET_KIND) :: OffsetCut0=0
  INTEGER(KIND=MPI_OFFSET_KIND), PRIVATE :: Offset
  TYPE(MPI_Status) :: Status
  TYPE(MPI_File) :: fh,fh1,fh2,fh3
  INTEGER :: OffsetShift
  INTEGER :: OffsetShiftGMV
  INTEGER :: OffsetShiftSoil
  INTEGER :: OffsetShiftSoilGMV
  INTEGER :: OffsetShiftCut
  INTEGER :: OffsetShiftCutGMV

  LOGICAL :: VelOut
  LOGICAL :: thOut
  LOGICAL :: EnOut
  LOGICAL :: tkeOut
  LOGICAL :: disOut
  LOGICAL :: omeOut
  LOGICAL :: tkeHOut
  LOGICAL :: tkeVOut
  LOGICAL :: LenOut
  LOGICAL :: RhoVOut
  LOGICAL :: RhoCOut
  LOGICAL :: RhoROut
  LOGICAL :: RhoIOut
  LOGICAL :: RhoSOut
  LOGICAL :: nvOut
  LOGICAL :: ncOut
  LOGICAL :: nrOut
  LOGICAL :: niOut
  LOGICAL :: nsOut
  LOGICAL :: tracer1Out
  LOGICAL :: tracer2Out
  LOGICAL :: RhoOut
  LOGICAL :: PreOut
  LOGICAL :: DampOut
  LOGICAL :: ShadOut
  LOGICAL :: RadiationOut
  LOGICAL :: FluxOut
  LOGICAL :: LandClassOut
  LOGICAL :: DiffOut
  LOGICAL :: DiffHOut
  LOGICAL :: DiffVOut
  LOGICAL :: DiffPotOut
  LOGICAL :: DiffPotHOut
  LOGICAL :: DiffPotVOut
  LOGICAL :: DiffMomOut
  LOGICAL :: DiffMomHOut
  LOGICAL :: DiffMomVOut
  LOGICAL :: ChemieOut

  TYPE(Vector4Cell_T), ALLOCATABLE :: VectorOut(:)
  TYPE OutputControl_T
    CHARACTER(32) :: NameScalar='                                   '
    CHARACTER*8 :: DiffRho
    CHARACTER*12 :: Type='Gas'
    REAL(RealKind) :: UnitKonvFactor=1.d0
    INTEGER :: iPosC=0
  END TYPE OutputControl_T
  TYPE(OutputControl_T), ALLOCATABLE :: OutputC(:)
! CHARACTER(32), ALLOCATABLE :: NameScalars(:)
  CHARACTER(32), ALLOCATABLE :: NameScalarsB(:)
  INTEGER :: VectorComponentsOut
  INTEGER :: VectorComponentsOutB ! Soil
  INTEGER :: VectorComponentsOutC ! Cut
  INTEGER, ALLOCATABLE :: GasSpeciesC(:)
! CHARACTER*4, ALLOCATABLE :: DiffRho(:) 

  REAL(RealKind) :: OutputTimeStart
  REAL(RealKind) :: OutputTimeStep
  REAL(RealKind) :: OutputTimeEnd
  CHARACTER(4) :: OutputType
  CHARACTER(5) :: OutputGrid
  CHARACTER(5) :: OutputVelocity
  INTEGER :: OutputStart
  INTEGER :: OutputFrequ
  INTEGER :: OutputEnd
  !REAL(RealKind), PRIVATE :: OutputTime
  !INTEGER, PRIVATE :: OutputStep
  !INTEGER, PRIVATE :: Step
  REAL(RealKind) :: OutputTime
  INTEGER :: OutputStep
  INTEGER :: Step
  INTEGER, PRIVATE :: gmvStep=1000
  CHARACTER(160) :: OutputFileName
  CHARACTER(160) :: OutputVarFileName
  CHARACTER(160) :: OutputDirName=''
  CHARACTER(80) :: OutputFileSoil
  CHARACTER(80) :: OutputFileCut
  INTEGER, PRIVATE :: NumTime=200
  REAL(RealKind) :: NumberMin=1.0d0

! Profiles
  REAL(RealKind) :: xProf,yProf
  REAL(RealKind) :: xStart,yStart
  REAL(RealKind) :: xEnd,yEnd
  REAL(RealKind) :: ProfTimeStart,ProfTimeEnd
  INTEGER:: ProfgmvStep=1000
  INTEGER:: ProfgmvEnd=10000
  LOGICAL :: ProfOut
  LOGICAL :: ProfAverage=.TRUE.
  !LOGICAL :: ProfAverage=.FALSE.
  LOGICAL :: Deviation=.FALSE.
  LOGICAL :: ProfIslandOcean=.FALSE.
  REAL(RealKind) :: CloudCover=9999.d0
  REAL(RealKind) :: CloudBaseHeight=0.0d0
  REAL(RealKind) :: CloudTopHeight=0.0d0
  REAL(RealKind) :: CloudCover2=9999.d0
  REAL(RealKind) :: CloudBaseHeight2=0.0d0
  REAL(RealKind) :: CloudTopHeight2=0.0d0
  REAL(RealKind) :: CloudCover3=9999.d0
  REAL(RealKind) :: CloudBaseHeight3=0.0d0
  REAL(RealKind) :: CloudTopHeight3=0.0d0
  REAL(RealKind) :: MissValue=1.0d36

  LOGICAL :: ReadOut=.FALSE.

  CHARACTER(1) :: TimeOutputFormat='S' ! S=seconds, M=minutes, H=hours

! Time Profiles
  LOGICAL :: TimeProfile1=.FALSE.
  LOGICAL :: TimeProfile2=.FALSE.
  INTEGER :: TimeProfileIx1=1
  INTEGER :: TimeProfileIy1=1
  INTEGER :: TimeProfileIx2=2
  INTEGER :: TimeProfileIy2=1
  LOGICAL :: TimeProfileRH=.FALSE.
  LOGICAL :: TimeProfileLWC=.FALSE.
  LOGICAL :: TimeProfileRhoV=.FALSE.
  LOGICAL :: TimeProfileRhoC=.FALSE.
  LOGICAL :: TimeProfileRhoR=.FALSE.
  LOGICAL :: TimeProfileNc=.FALSE.
  LOGICAL :: TimeProfileNr=.FALSE.
  LOGICAL :: TimeProfileNCCN=.FALSE.
  LOGICAL :: TimeProfileNTotal=.FALSE.
  LOGICAL :: TimeProfileTracer1=.FALSE.
  LOGICAL :: TimeProfileTracer2=.FALSE.
  LOGICAL :: TimeProfileW=.FALSE.
  LOGICAL :: TimeProfileAbsVel=.FALSE.
  LOGICAL :: TimeProfileT=.FALSE.
  LOGICAL :: TimeProfileTC=.FALSE.
  LOGICAL :: TimeProfilePotTemp=.FALSE.
  LOGICAL :: TimeProfileThDens=.FALSE.
  LOGICAL :: TimeProfileThEquiv=.FALSE.


  INTEGER :: ibb
  REAL(RealKind), PRIVATE, POINTER :: c(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uC(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vC(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wC(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: f(:,:,:)
  INTEGER, PRIVATE :: NumberOfNodes
  INTEGER, PRIVATE :: NumberOfCells
  INTEGER, PARAMETER :: NumberOfMaterials=2
  CHARACTER(32) :: Materials(1:NumberOfMaterials) &
                  =(/'                            Wand' &
                    ,'                            Luft'/)
  CHARACTER(360) :: gmvFile
  CHARACTER(8),  PRIVATE, SAVE :: gmvOutputType
  INTEGER, PARAMETER :: NameByte=8
  CHARACTER(8), PARAMETER :: gmvinput ='gmvinput'
  CHARACTER(8), PARAMETER :: ascii    ='   ascii'
  CHARACTER(8), PARAMETER :: ieeei4r4 ='ieeei4r4'
  CHARACTER(8), PARAMETER :: ieeei4r8 ='ieeei4r8'
  CHARACTER(8), PARAMETER :: iecxi4r4 ='iecxi4r4'
  CHARACTER(8), PARAMETER :: iecxi4r8 ='iecxi4r8'
  CHARACTER(8), PARAMETER :: nodes    ='nodes   '
  CHARACTER(8), PARAMETER :: polygons ='polygons'
  CHARACTER(9), PARAMETER :: spolygons ='polygons '
  CHARACTER(9), PARAMETER :: probtime ='probtime '
  CHARACTER(8), PARAMETER :: nodev    ='nodev   '
  CHARACTER(8), PARAMETER :: cells    ='cells   '
  CHARACTER(8), PARAMETER :: material ='material'
  CHARACTER(8), PARAMETER :: velocity ='velocity'
  CHARACTER(8), PARAMETER :: hex      ='hex     '
  CHARACTER(8), PARAMETER :: endflag  ='endflag '
  CHARACTER(8), PARAMETER :: fromfile ='fromfile'
  CHARACTER(4), PARAMETER :: suffix_out    ='.out'
  CHARACTER(8), PARAMETER :: suffix_outgmvG='.out.gmvG'
  REAL(RealKind), PARAMETER :: RadOutput=1.0d4

  CHARACTER(8), ALLOCATABLE :: AeroOutUnit(:)
  REAL(RealKind), ALLOCATABLE :: AeroUnitKonvFactor(:)
  CHARACTER(8), ALLOCATABLE :: GasOutUnit(:)
  REAL(RealKind), ALLOCATABLE :: GasUnitKonvFactor(:)

  INTERFACE Output
    MODULE PROCEDURE OutputMet,OutputMetC,OutputVector
  END INTERFACE

  NAMELIST /ModelOutput/ VelOut, &
                         thOut,  &
                         EnOut,  &
                         tkeOut, &
                         disOut, &
                         omeOut, &
                         tkeHOut,&
                         tkeVOut,&
                         LenOut, &
                         RhoVOut,  &
                         RhoCOut,  &
                         RhoROut,  &
                         RhoIOut,  &
                         RhoSOut,  &
                         nvOut,  &
                         ncOut,  &
                         nrOut,  &
                         niOut,  &
                         nsOut,  &
                         tracer1Out,  &
                         tracer2Out,  &
                         RhoOut,&
                         PreOut,&
                         DampOut,&
                         ShadOut,&
                         RadiationOut,&
                         FluxOut,&
                         LandClassOut,&
                         DiffOut,&
                         DiffHOut,&
                         DiffVOut,&                         
                         DiffPotOut,&
                         DiffPotHOut,&
                         DiffPotVOut,&
                         DiffMomOut,&
                         DiffMomHOut,&
                         DiffMomVOut,&
                         ChemieOut,&
                         OutputTimeStart, &
                         OutputTimeStep, &
                         OutputTimeEnd,  &
                         OutputStart, &
                         OutputFrequ, &
                         OutputEnd,  &
                         OutputFileName, &
                         OutputVarFileName, &
                         OutputDirName, &
                         OutputType,    &
                         OutputGrid,    &
                         OutputVelocity,    &
                         ReadOut,    &
                         ProfOut,    &
                         ProfAverage,&
                         Deviation,  &
                         ProfIslandOcean,&
                         xStart,     &
                         xEnd,       &
                         yStart,     &
                         yEnd,       &
                         xProf,      &
                         yProf,      &
                         ProfTimeStart, &
                         ProfTimeEnd,   &
                         ProfgmvStep,   &
                         ProfgmvEnd,   &
                         MissValue,     &
                         NumberMin,     &
                         TimeOutputFormat, &
                         TimeProfile1, &
                         TimeProfile2, &  
                         TimeProfileIx1, &
                         TimeProfileIy1, &  
                         TimeProfileIx2, &
                         TimeProfileIy2, &  
                         TimeProfileRH,  &
                         TimeProfileRhoV,  &
                         TimeProfileRhoC,  &
                         TimeProfileRhoR,  &
                         TimeProfileNCCN,  &
                         TimeProfileNc,  &
                         TimeProfileNr,  &
                         TimeProfileTracer1,  &
                         TimeProfileTracer2,  &
                         TimeProfileNTotal,  &
                         TimeProfileLWC,     &
                         TimeProfileW, &
                         TimeProfileAbsVel,  &
                         TimeProfileT, &
                         TimeProfileTC, &
                         TimeProfilePotTemp, &
                         TimeProfileThDens, &
                         TimeProfileThEquiv

  TYPE(Vector4Cell_T), POINTER :: OutSpecial(:)
  INTEGER :: LenOutSpecial=1
  CHARACTER(32), ALLOCATABLE :: NameOutSpecial(:)
  CHARACTER(8), ALLOCATABLE :: UnitOutSpecial(:)
  CHARACTER(12), ALLOCATABLE :: SpecialType(:)
  REAL(RealKind), ALLOCATABLE :: UnitKonvSpecial(:)
  
CONTAINS

SUBROUTINE SetgmvStep(gmvStepAct)
  INTEGER :: gmvStepAct
  gmvStep=gmvStepAct
END SUBROUTINE SetgmvStep

FUNCTION GetgmvStep()
  INTEGER :: GetgmvStep
  GetgmvStep=gmvStep
END FUNCTION GetgmvStep

SUBROUTINE SetOutputTime(OutputTimeAct)
  REAL(RealKind) :: OutputTimeAct
  OutputTime=OutputTimeAct
END SUBROUTINE SetOutputTime

FUNCTION GetOutputTime()
  REAL(RealKind) :: GetOutputTime
  GetOutputTime=OutputTime
END FUNCTION GetOutputTime

SUBROUTINE InputOutputChemie(FileName)

  CHARACTER(*) :: FileName

  CHARACTER(20) :: S1,S2,End,SpeciesName
  CHARACTER(4) :: iFrac
  CHARACTER(8) :: OutUnit
  LOGICAL :: Back

  INTEGER :: i,j,Pos

  S1='BEGIN_GAS'
  S2='BEGIN_OUTPUT'
  End='END_OUTPUT'
  CALL OpenFile(FileName)


  NumGasOut=0
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,Name2=OutUnit)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      NumGasOut=NumGasOut+1
    END IF
  END DO
  CALL CloseFile

  ALLOCATE(GasOut(NumGasOut))
  ALLOCATE(GasSpecies(NumGasOut))
  ALLOCATE(GasOutUnit(NumGasOut))

  NumGasOut=0
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,Name2=OutUnit)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      NumGasOut=NumGasOut+1
      GasOut(NumGasOut)=Pos
      GasOutUnit(NumGasOut)=OutUnit
      GasSpecies(NumGasOut)=TRIM(SpeciesName)
    END IF
  END DO
  CALL CloseFile
! Check for right spelling of units:
  DO j=1,NumgasOut
    DO i=1,LEN_TRIM(GasOutUnit(j))
      IF (GasOutUnit(j)(i:i).EQ.'n') GasOutUnit(j)(i:i)='N'
      IF (GasOutUnit(j)(i:i).EQ.'k') GasOutUnit(j)(i:i)='K'
      IF (GasOutUnit(j)(i:i).EQ.'u') GasOutUnit(j)(i:i)='U'
      IF (GasOutUnit(j)(i:i).EQ.'g') GasOutUnit(j)(i:i)='G'
      IF (GasOutUnit(j)(i:i).EQ.'m') GasOutUnit(j)(i:i)='M'
      IF (GasOutUnit(j)(i:i).EQ.'o') GasOutUnit(j)(i:i)='O'
      IF (GasOutUnit(j)(i:i).EQ.'l') GasOutUnit(j)(i:i)='L'
      IF (GasOutUnit(j)(i:i).EQ.'p') GasOutUnit(j)(i:i)='P'
      IF (GasOutUnit(j)(i:i).EQ.'a') GasOutUnit(j)(i:i)='A'
      IF (GasOutUnit(j)(i:i).EQ.'r') GasOutUnit(j)(i:i)='R'
      IF (GasOutUnit(j)(i:i).EQ.'t') GasOutUnit(j)(i:i)='T'
    END DO
    GasSpecies(j)=TRIM(GasSpecies(j))//"_"//TRIM(GasOutUnit(j))
  END DO
  ALLOCATE(GasUnitKonvFactor(NumGasOut))
  DO j=1,NumGasOut
    IF (GasOutUnit(j)(1:2).EQ.'KG') GasUnitKonvFactor(j)=1.d0
    IF (GasOutUnit(j)(1:1).EQ.'G')  GasUnitKonvFactor(j)=1.d3
    IF (GasOutUnit(j)(1:2).EQ.'MG') GasUnitKonvFactor(j)=1.d6
    IF (GasOutUnit(j)(1:3).EQ.'MUG') GasUnitKonvFactor(j)=1.d9
    IF (GasOutUnit(j)(1:2).EQ.'NG') GasUnitKonvFactor(j)=1.d12
  END DO


  S1='BEGIN_AERO'
  S2='BEGIN_OUTPUT'
  End='END_OUTPUT'
  CALL OpenFile(FileName)
  NumAeroOut=0
  IF (Position('aNUMBER')>0)  NumAeroOut=1 ! set 'aNUMBER' as first value, when present
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,Name2=OutUnit)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      NumAeroOut=NumAeroOut+1
    END IF
  END DO
  CALL CloseFile
  ALLOCATE(AeroOut(NumAeroOut))
  ALLOCATE(AeroOutUnit(NumAeroOut))
  ALLOCATE(AeroSpecies(NumAeroOut,0:nFrac))
  IF (Position('aNUMBER')>0) THEN ! set 'aNUMBER' as first value, when present
    NumAeroOut=1
    AeroOut(NumAeroOut)=Position('aNUMBER')
    AeroSpecies(NumAeroOut,0)='aNUMBER'
  ELSE
    NumAeroOut=0
  END IF
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,Name2=OutUnit)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      NumAeroOut=NumAeroOut+1
      AeroOut(NumAeroOut)=Pos
      AeroOutUnit(NumAeroOut)=OutUnit
      AeroSpecies(NumAeroOut,0)=TRIM(SpeciesName)
    END IF
  END DO
  CALL CloseFile
! Check for right spelling of units:
  DO j=2,NumAeroOut
    DO i=1,LEN_TRIM(AeroOutUnit(j))
      IF (AeroOutUnit(j)(i:i).EQ.'n') AeroOutUnit(j)(i:i)='N'
      IF (AeroOutUnit(j)(i:i).EQ.'k') AeroOutUnit(j)(i:i)='K'
      IF (AeroOutUnit(j)(i:i).EQ.'u') AeroOutUnit(j)(i:i)='U'
      IF (AeroOutUnit(j)(i:i).EQ.'g') AeroOutUnit(j)(i:i)='G'
      IF (AeroOutUnit(j)(i:i).EQ.'m') AeroOutUnit(j)(i:i)='M'
      IF (AeroOutUnit(j)(i:i).EQ.'o') AeroOutUnit(j)(i:i)='O'
      IF (AeroOutUnit(j)(i:i).EQ.'l') AeroOutUnit(j)(i:i)='L'
      IF (AeroOutUnit(j)(i:i).EQ.'p') AeroOutUnit(j)(i:i)='P'
      IF (AeroOutUnit(j)(i:i).EQ.'a') AeroOutUnit(j)(i:i)='A'
      IF (AeroOutUnit(j)(i:i).EQ.'r') AeroOutUnit(j)(i:i)='R'
      IF (AeroOutUnit(j)(i:i).EQ.'t') AeroOutUnit(j)(i:i)='T'
    END DO
    AeroSpecies(j,0)=TRIM(AeroSpecies(j,0))//"_"//TRIM(AeroOutUnit(j))
  END DO

! Add bin number
  DO j=1,NumAeroOut
    IF (nFrac>1) THEN
      DO i=1,nFrac
        WRITE(iFrac,'(I4)') i
        iFrac = ADJUSTL(iFrac)
        AeroSpecies(j,i)=TRIM(AeroSpecies(j,0))//"_"//TRIM(iFrac)
      END DO
    ELSE
      AeroSpecies(j,nFrac)=TRIM(AeroSpecies(j,0))
    END IF
  END DO

  ALLOCATE(AeroUnitKonvFactor(NumAeroOut))
  DO j=1,NumAeroOut
    IF (AeroOutUnit(j)(1:2).EQ.'KG') AeroUnitKonvFactor(j)=1.d0
    IF (AeroOutUnit(j)(1:1).EQ.'G')  AeroUnitKonvFactor(j)=1.d3
    IF (AeroOutUnit(j)(1:2).EQ.'MG') AeroUnitKonvFactor(j)=1.d6
    IF (AeroOutUnit(j)(1:3).EQ.'MUG') AeroUnitKonvFactor(j)=1.d9
    IF (AeroOutUnit(j)(1:2).EQ.'NG') AeroUnitKonvFactor(j)=1.d12
  END DO

END SUBROUTINE InputOutputChemie

SUBROUTINE InputOutputSpecial(FileName)

  CHARACTER(*) :: FileName

  CHARACTER(20) :: S1,S2,End,SpecialName
  LOGICAL :: Back

  INTEGER :: i

  S1='BEGIN_SPECIAL'
  S2='BEGIN_OUTPUT'
  End='END_OUTPUT'
  CALL OpenFile(FileName)

  LenOutSpecial=0
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpecialName)
    IF (Back) THEN
      EXIT
    END IF
    IF (SpecialName.EQ.'PM1' .OR. SpecialName.EQ.'PM2.5' .OR.SpecialName.EQ.'PM10') THEN
      LenOutSpecial=LenOutSpecial+NumAeroOut
    ELSE
      LenOutSpecial=LenOutSpecial+1
    END IF
  END DO
  CALL CloseFile

  ALLOCATE(NameOutSpecial(LenOutSpecial))
  ALLOCATE(UnitOutSpecial(LenOutSpecial))
  ALLOCATE(UnitKonvSpecial(LenOutSpecial))
  ALLOCATE(SpecialType(LenOutSpecial))
  UnitKonvSpecial=1.d0
  CALL Allocate(OutSpecial,LenOutSpecial)

  LenOutSpecial=0
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpecialName)
    IF (Back) THEN
      EXIT
    END IF
    IF (SpecialName.EQ.'PM1' .OR. SpecialName.EQ.'PM2.5' .OR.SpecialName.EQ.'PM10') THEN
! Particel number
      LenOutSpecial=LenOutSpecial+1
      NameOutSpecial(LenOutSpecial)=TRIM(ADJUSTL(SpecialName))//'_'//TRIM(ADJUSTL(AeroSpecies(1,0)))
      UnitOutSpecial(LenOutSpecial)='Rho'
      UnitKonvSpecial(LenOutSpecial)=1.d0
      SpecialType(LenOutSpecial)='Gas'
! Species in Particel
      DO i=2,NumAeroOut
        LenOutSpecial=LenOutSpecial+1
        NameOutSpecial(LenOutSpecial)=TRIM(ADJUSTL(SpecialName))//'_'//TRIM(ADJUSTL(AeroSpecies(i,0)))
        UnitOutSpecial(LenOutSpecial)=AeroOutUnit(i)
        UnitKonvSpecial(LenOutSpecial)=AeroUnitKonvFactor(i)
        SpecialType(LenOutSpecial)='Aerosol_Spec'
      END DO
    ELSE
      LenOutSpecial=LenOutSpecial+1
      NameOutSpecial(LenOutSpecial)=SpecialName
      UnitOutSpecial(LenOutSpecial)='Unsc'
      SpecialType(LenOutSpecial)='Gas'
    END IF
  END DO
  CALL CloseFile

END SUBROUTINE InputOutputSpecial

SUBROUTINE gmvProbtimeWrite(ActTime)

  REAL(RealKind) :: ActTime

  REAL(RealKind) :: ActTimeOutput
  REAL(RealKind) :: Hour,Minutes,Seconds
  REAL(RealKind) :: Time,Temp
  REAL(Real4Kind) :: TimeOut
  

  IF (TimeOutputFormat=='M') THEN
    ActTimeOutput=ActTime/60.0d0
    Minutes=AINT(ActTimeOutput)
    Seconds=(ActTimeOutput-Minutes)*.6d0
    ActTimeOutput=Minutes+Seconds
  ELSE IF (TimeOutputFormat=='H') THEN
    ActTimeOutput=ActTime/3600.0d0
    Hour=AINT(ActTimeOutput)
    Minutes=(ActTimeOutput-Hour)*.6d0
    ActTimeOutput=Hour+Minutes
  ELSE IF (TimeOutputFormat=='S') THEN
    ActTimeOutput=ActTime
  ELSE
    ActTimeOutput=ActTime
  END IF

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a9,f16.2)') probtime,ActTimeOutput
  ELSE
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_WRITE_AT(fh,Offset,probtime(1:8),1, &
                        MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8

      TimeOut=ActTimeOutput

      CALL native_4byte_real(TimeOut,TimeOut)
      CALL MPI_FILE_WRITE_AT(fh,Offset,TimeOut,1, &
                        MPI_REAL,Status,MPIErr)
    END IF
    Offset0=Offset0+8+Real4Kind
  END IF

END SUBROUTINE gmvProbtimeWrite

SUBROUTINE gmvProbtimeWriteSoil(ActTime)

  REAL(RealKind) :: ActTime

  INTEGER :: MPIErr
  REAL(RealKind) :: ActTimeOutput
  REAL(RealKind) :: Hour,Minutes,Seconds
  REAL(RealKind) :: Time,Temp
  REAL(Real4Kind) :: TimeOut
  

  IF (TimeOutputFormat=='M') THEN
    ActTimeOutput=ActTime/60.0d0
    Minutes=AINT(ActTimeOutput)
    Seconds=(ActTimeOutput-Minutes)*.6d0
    ActTimeOutput=Minutes+Seconds
  ELSE IF (TimeOutputFormat=='H') THEN
    ActTimeOutput=ActTime/3600.0d0
    Hour=AINT(ActTimeOutput)
    Minutes=(ActTimeOutput-Hour)*.6d0
    ActTimeOutput=Hour+Minutes
  ELSE IF (TimeOutputFormat=='S') THEN
    ActTimeOutput=ActTime
  ELSE
    ActTimeOutput=ActTime
  END IF

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a9,f16.2)') probtime,ActTimeOutput
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetSoil0
      CALL MPI_FILE_WRITE_AT(fh2,Offset,probtime(1:8),1, &
                        MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8

      TimeOut=ActTimeOutput

      CALL native_4byte_real(TimeOut,TimeOut)
      CALL MPI_FILE_WRITE_AT(fh2,Offset,TimeOut,1, &
                        MPI_REAL,Status,MPIErr)
    END IF
    OffsetSoil0=OffsetSoil0+8+Real4Kind
  END IF

END SUBROUTINE gmvProbtimeWriteSoil

SUBROUTINE gmvProbtimeWriteCut(ActTime)

  REAL(RealKind) :: ActTime

  INTEGER :: MPIErr
  REAL(RealKind) :: ActTimeOutput
  REAL(RealKind) :: Hour,Minutes,Seconds
  REAL(RealKind) :: Time,Temp
  REAL(Real4Kind) :: TimeOut
  

  IF (TimeOutputFormat=='M') THEN
    ActTimeOutput=ActTime/60.0d0
    Minutes=AINT(ActTimeOutput)
    Seconds=(ActTimeOutput-Minutes)*.6d0
    ActTimeOutput=Minutes+Seconds
  ELSE IF (TimeOutputFormat=='H') THEN
    ActTimeOutput=ActTime/3600.0d0
    Hour=AINT(ActTimeOutput)
    Minutes=(ActTimeOutput-Hour)*.6d0
    ActTimeOutput=Hour+Minutes
  ELSE IF (TimeOutputFormat=='S') THEN
    ActTimeOutput=ActTime
  ELSE
    ActTimeOutput=ActTime
  END IF

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a9,f16.2)') probtime,ActTimeOutput
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetCut0
      CALL MPI_FILE_WRITE_AT(fh3,Offset,probtime(1:8),1, &
                        MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8

      TimeOut=ActTimeOutput

      CALL native_4byte_real(TimeOut,TimeOut)
      CALL MPI_FILE_WRITE_AT(fh3,Offset,TimeOut,1, &
                        MPI_REAL,Status,MPIErr)
    END IF
    OffsetCut0=OffsetCut0+8+Real4Kind
  END IF

END SUBROUTINE gmvProbtimeWriteCut

SUBROUTINE gmvProbtimeRead(ActTime)

  REAL(RealKind) :: ActTime

  REAL(RealKind) :: Hours,Minutes,Seconds
  REAL(RealKind) :: Time,Temp
  REAL(Real4Kind) :: GMVTime
  CHARACTER*8 :: Dummy
  

  IF (gmvOutputType==ascii) THEN
    READ(OutputUnit,'(a9,f16.2)') Dummy,ActTime
  ELSE
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_READ_AT(fh,Offset,Dummy(1:8),1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      WRITE(*,*) 'Dummy ',Dummy(1:8)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh,Offset,GMVTime,1, &
                      MPI_REAL,Status,MPIErr)
      CALL MPI_Bcast(GMVTime,1,MPI_REAL,0,MPI_COMM_WORLD,MPIErr)
    ELSE  
      CALL MPI_Bcast(GMVTime,1,MPI_REAL,0,MPI_COMM_WORLD,MPIErr)
    END IF
    Offset0=Offset0+8+Real4Kind
  END IF
  CALL native_4byte_real(GMVTime,GMVTime)

  WRITE(*,*) 'TimeOutputFormat ',TimeOutputFormat,GMVTime
  IF (TimeOutputFormat=='M') THEN 
    Minutes=AINT(GMVTime)
    Seconds=GMVTime-Minutes
    ActTime=Minutes*60.0d0+Seconds  
  ELSE IF (TimeOutputFormat=='H') THEN 
    Hours=AINT(GMVTime)
    Minutes=GMVTime-Hours
    ActTime=Hours*3600.0d0+Minutes*60.0d0
  ELSE 
    ActTime=GMVTime
  END IF
END SUBROUTINE gmvProbtimeRead

SUBROUTINE gmvOpenWrite(FileName)

  CHARACTER(*) :: FileName

  gmvFile=FileName
  IF (gmvOutputType==ascii) THEN
    OPEN(UNIT=OutputUnit,FILE=FileName,STATUS='UNKNOWN')
    WRITE(OutputUnit,'(a8,a8)') gmvinput,gmvOutputType
  ELSE
    Offset0=0
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,FileName, &
                       MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                       MPI_INFO_NULL,fh,MPIErr)
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_WRITE_AT(fh,Offset,gmvinput,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh,Offset,gmvOutputType,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
    END IF
    Offset0=Offset0+LEN(gmvinput)+LEN(gmvOutputType)
  END IF

END SUBROUTINE gmvOpenWrite


SUBROUTINE gmvOpenSoilWrite(FileName)

  CHARACTER(*) :: FileName

  gmvFile=FileName
  IF (gmvOutputType==ascii) THEN
    OPEN(UNIT=OutputUnitS,FILE=FileName,STATUS='UNKNOWN')
    WRITE(OutputUnitS,'(a8,a8)') gmvinput,gmvOutputType
  ELSE
    OffsetSoil0=0
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,FileName, &
                       MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                       MPI_INFO_NULL,fh2,MPIErr)
    IF (MyId==0) THEN
      Offset=OffsetSoil0
      CALL MPI_FILE_WRITE_AT(fh2,Offset,gmvinput,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh2,Offset,gmvOutputType,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
    END IF
    OffsetSoil0=OffsetSoil0+LEN(gmvinput)+LEN(gmvOutputType)
  END IF

END SUBROUTINE gmvOpenSoilWrite


SUBROUTINE gmvOpenCutWrite(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: MPIErr

  gmvFile=FileName
  IF (gmvOutputType==ascii) THEN
    OPEN(UNIT=OutputUnitS,FILE=FileName,STATUS='UNKNOWN')
    WRITE(OutputUnitS,'(a8,a8)') gmvinput,gmvOutputType
  ELSE
    OffsetCut0=0
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,FileName, &
                       MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                       MPI_INFO_NULL,fh3,MPIErr)
    IF (MyId==0) THEN
      Offset=OffsetCut0
      CALL MPI_FILE_WRITE_AT(fh3,Offset,gmvinput,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh3,Offset,gmvOutputType,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
    END IF
    OffsetCut0=OffsetCut0+LEN(gmvinput)+LEN(gmvOutputType)
  END IF

END SUBROUTINE gmvOpenCutWrite


SUBROUTINE gmvOpenRead(FileName)

  CHARACTER(*) :: FileName

  CHARACTER*8 :: Dummy

  INTEGER ERRORCLASS,LenString
  CHARACTER(MPI_MAX_ERROR_STRING) :: String

  gmvFile=FileName
  IF (gmvOutputType==ascii) THEN
    OPEN(UNIT=OutputUnit,FILE=FileName,STATUS='UNKNOWN')
    READ(OutputUnit,'(a8,a8)') Dummy,Dummy
  ELSE
    Offset0=0
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,FileName, &
!                      MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL,fh,MPIErr)
  CALL MPI_ERROR_CLASS(MPIErr, ERRORCLASS, MPIErr)
  WRITE(*,*) 'MPIErr nach MPI_FILE_CLOSE',MPIErr,ERRORCLASS
  CALL MPI_ERROR_STRING(ERRORCLASS,String,LenString,MPIErr)
  WRITE(*,*) 'StringOpenError   ',String
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_READ_AT(fh,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
    END IF
    Offset0=Offset0+LEN(gmvinput)+LEN(gmvOutputType)
  END IF

END SUBROUTINE gmvOpenRead

SUBROUTINE gmvOpenSoilRead(FileName)

  CHARACTER(*) :: FileName

  CHARACTER*8 :: Dummy

  gmvFile=FileName
  IF (gmvOutputType==ascii) THEN
    OPEN(UNIT=OutputUnitS,FILE=FileName,STATUS='UNKNOWN')
    READ(OutputUnitS,'(a8,a8)') Dummy,Dummy
  ELSE
    OffsetSoil0=0
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,FileName, &
                        MPI_MODE_RDONLY, &
                        MPI_INFO_NULL,fh2,MPIErr)
    IF (MyId==0) THEN
      Offset=OffsetSoil0
      CALL MPI_FILE_READ_AT(fh2,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh2,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
    END IF
    OffsetSoil0=OffsetSoil0+LEN(gmvinput)+LEN(gmvOutputType)
  END IF

END SUBROUTINE gmvOpenSoilRead


SUBROUTINE gmvOpenCutRead(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: MPIErr
  CHARACTER*8 :: Dummy

  gmvFile=FileName
  IF (gmvOutputType==ascii) THEN
    OPEN(UNIT=OutputUnitS,FILE=FileName,STATUS='UNKNOWN')
    READ(OutputUnitS,'(a8,a8)') Dummy,Dummy
  ELSE
    OffsetCut0=0
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,FileName, &
                        MPI_MODE_RDONLY, &
                        MPI_INFO_NULL,fh3,MPIErr)
    IF (MyId==0) THEN
      Offset=OffsetCut0
      CALL MPI_FILE_READ_AT(fh3,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh3,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
    END IF
    OffsetCut0=OffsetCut0+LEN(gmvinput)+LEN(gmvOutputType)
  END IF

END SUBROUTINE gmvOpenCutRead


SUBROUTINE gmvClose

  CHARACTER(8), PARAMETER :: endgmv  ='endgmv  '

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a8)') endgmv
    CLOSE(UNIT=OutputUnit)
  ELSE
    IF (MyId==0) THEN
      CALL MPI_FILE_WRITE_AT(fh,Offset0,endgmv,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    CALL MPI_FILE_CLOSE(fh,MPIErr)
  END IF

END SUBROUTINE gmvClose


SUBROUTINE gmvCloseSoil

  CHARACTER(8), PARAMETER :: endgmv  ='endgmv  '

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnitS,'(a8)') endgmv
    CLOSE(UNIT=OutputUnit)
  ELSE
    IF (MyId==0) THEN
      CALL MPI_FILE_WRITE_AT(fh2,OffsetSoil0,endgmv,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    CALL MPI_FILE_CLOSE(fh2,MPIErr)
  END IF

END SUBROUTINE gmvCloseSoil


SUBROUTINE gmvCloseCut

  CHARACTER(8), PARAMETER :: endgmv  ='endgmv  '

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnitS,'(a8)') endgmv
    CLOSE(UNIT=OutputUnit)
  ELSE
    IF (MyId==0) THEN
      CALL MPI_FILE_WRITE_AT(fh3,OffsetCut0,endgmv,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    CALL MPI_FILE_CLOSE(fh3,MPIErr)
  END IF

END SUBROUTINE gmvCloseCut


SUBROUTINE gmvGeometry


  INTEGER :: i,ix,iy,iz,j
  INTEGER :: iw
  INTEGER :: v000,v001,v010,v100,v011,v101,v110,v111
  INTEGER :: NumberOfNodesB
  INTEGER :: NumberOfCellsB
  INTEGER :: Node
  INTEGER :: Cell
  INTEGER, ALLOCATABLE :: iWork(:)
  REAL(4), ALLOCATABLE :: Work(:)
  REAL(4) :: xN,yN,zN

  IF (gmvOutputType==ascii) THEN
!   Number of Nodes and Cells
    NumberOfNodes=0
    NumberOfCells=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      NumberOfNodes=NumberOfNodes+(nx+1)*(ny+1)*(nz+1)
      NumberOfCells=NumberOfCells+nx*ny*nz
    END DO
    WRITE(OutputUnit,'(a8,i8)') nodev,NumberOfNodes
    IF (OutputGrid=='Plane') THEN
      Node=1
      DO ibb=1,nb
        CALL Set(Floor(ibb))
        DO iz=iz0,iz1
          DO iy=iy0,iy1
            DO ix=ix0,ix1
              xN=xP(ix)
              yN=yP(iy)
              zN=zP(iz)
              WRITE(OutputUnit,*) xN,yN,zN
            END DO
          END DO
        END DO
      END DO
    ELSE IF (OutputGrid=='Globe') THEN
      Node=1
      DO ibb=1,nb
        CALL Set(Floor(ibb))
        DO iz=iz0,iz1
          DO iy=iy0,iy1
            DO ix=ix0,ix1
              xN=(zP(iz)+RadOutput)*SIN(xP(ix))*COS(yp(iy))/RadOutput
              yN=(zP(iz)+RadOutput)*COS(xP(ix))*COS(yp(iy))/RadOutput
              zN=(zP(iz)+RadOutput)*SIN(yp(iy))/RadOutput
              WRITE(OutputUnit,*) xN,yN,zN
            END DO
          END DO
        END DO
      END DO
    END IF
    WRITE(OutputUnit,'(a8,i8)') cells,NumberOfCells
    NumberOfNodes=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      v000=NumberOfNodes
      NumberOfNodes=NumberOfNodes+(nx+1)*(ny+1)*(nz+1)
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            v000=v000+1
            v100=v000+1
            v010=v000+(nx+1)
            v110=v010+1 

            v001=v000+(nx+1)*(ny+1)
            v101=v001+1
            v011=v001+(nx+1)
            v111=v011+1 
            WRITE(OutputUnit,'(a8,I8)') hex,8 
            WRITE(OutputUnit,'(8I8)') v001,v101,v111,v011,v000,v100,v110,v010
          END DO
          v000=v000+1
        END DO
        v000=v000+(nx+1)
      END DO
    END DO
    WRITE(OutputUnit,'(a8,i8,i8)') material,NumberOfMaterials,0
    DO i=1,NumberOfMaterials
      WRITE(OutputUnit,'(a32)') Materials(i)
    END DO
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)==Zero) THEN
              WRITE(OutputUnit,*) 0
            ELSE
              WRITE(OutputUnit,*) 2
            END IF
          END DO
        END DO
      END DO
    END DO
  ELSE
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_WRITE_AT(fh,Offset,nodev,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh,Offset,NumberOfNodes,1, &
                      MPI_INTEGER,Status,MPIErr)
    END IF
    Offset0=Offset0+LEN(nodev)+IntKind
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      Offset=Offset0+3*Real4Kind*WriteOffsetN
      NumberOfNodesB=(ix1-ix0+1)*(iy1-iy0+1)*(iz1-iz0+1)
      IF (NumberOfNodesB>=8) THEN
        ALLOCATE(Work(3*NumberOfNodesB))
        IF (OutputGrid=='Plane') THEN
          iw=1
          DO iz=iz0,iz1
            DO iy=iy0,iy1
              DO ix=ix0,ix1
                xN=xP(ix)
                yN=yP(iy)
                zN=zP(iz)
                Work(iw)=xN
                iw=iw+1
                Work(iw)=yN
                iw=iw+1
                Work(iw)=zN
                iw=iw+1
              END DO
            END DO
          END DO
        ELSE IF (OutputGrid=='Globe') THEN
          iw=1
          DO iz=iz0,iz1
            DO iy=iy0,iy1
              DO ix=ix0,ix1
              xN=(zP(iz)+RadOutput)*SIN(xP(ix))*COS(yp(iy))/RadOutput
              yN=(zP(iz)+RadOutput)*COS(xP(ix))*COS(yp(iy))/RadOutput
              zN=(zP(iz)+RadOutput)*SIN(yp(iy))/RadOutput
                Work(iw)=xN
                iw=iw+1
                Work(iw)=yN
                iw=iw+1
                Work(iw)=zN
                iw=iw+1
              END DO
            END DO
          END DO
        END IF
        CALL MPI_FILE_WRITE_AT(fh,Offset,Work,3*NumberOfNodesB, &
                               MPI_REAL,Status,MPIErr)
        DEALLOCATE(Work)
      END IF
    END DO

    Offset0=Offset0+3*Real4Kind*NumberOfNodes
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_WRITE_AT(fh,Offset,cells,1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh,Offset,NumberOfCells,1, &
                      MPI_INTEGER,Status,MPIErr)
    END IF
    Offset0=Offset0+LEN(Cells)+IntKind
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      nx=ix1-ix0
      ny=iy1-iy0
      nz=iz1-iz0
      NumberOfCellsB=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
      Offset=Offset0+WriteOffsetC*(LEN(hex)+9*IntKind)
      v000=WriteOffsetN
      ALLOCATE(iWork(11*NumberOfCellsB))
      iw=1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            v000=v000+1
            v100=v000+1
            v010=v000+(nx+1)
            v110=v010+1 

            v001=v000+(nx+1)*(ny+1)
            v101=v001+1
            v011=v001+(nx+1)
            v111=v011+1 
            iWork(iw)=TRANSFER(hex(1:4),iWork(iw))
            iw=iw+1
            iWork(iw)=TRANSFER(hex(5:8),iWork(iw))
            iw=iw+1
            iWork(iw)=8
            iw=iw+1
            iWork(iw)=v001
            iw=iw+1
            iWork(iw)=v101
            iw=iw+1
            iWork(iw)=v111
            iw=iw+1
            iWork(iw)=v011
            iw=iw+1
            iWork(iw)=v000
            iw=iw+1
            iWork(iw)=v100
            iw=iw+1
            iWork(iw)=v110
            iw=iw+1
            iWork(iw)=v010
            iw=iw+1
          END DO
          v000=v000+1
        END DO
        v000=v000+(nx+1)
      END DO
      CALL MPI_FILE_WRITE_AT(fh,Offset,iWork, &
                             11*NumberOfCellsB, &
                             MPI_INTEGER,Status,MPIErr)
      DEALLOCATE(iWork)
    END DO
    Offset0=Offset0+NumberOfCells*(LEN(hex)+9*IntKind)
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_WRITE_AT(fh,Offset,material,1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh,Offset,NumberOfMaterials,1, &
                      MPI_INTEGER,Status,MPIErr)
      Offset=Offset+4
      CALL MPI_FILE_WRITE_AT(fh,Offset,0,1, &
                      MPI_INTEGER,Status,MPIErr)
      Offset=Offset+4
      ALLOCATE(iWork(NumberOfMaterials*8))
      iw=1
      DO i=1,NumberOfMaterials
        DO j=1,32,4
          iWork(iw)=TRANSFER(Materials(i)(j:j+3),iWork(iw))
          iw=iw+1
        END DO
      END DO
      CALL MPI_FILE_WRITE_AT(fh,Offset,iWork,NumberOfMaterials*8, &
                      MPI_INTEGER,Status,MPIErr)
      DEALLOCATE(iWork)
    END IF
    Offset0=Offset0+LEN(material)+2*IntKind+32*NumberOfMaterials

    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      nx=ix1-ix0
      ny=iy1-iy0
      nz=iz1-iz0
      NumberOfCellsB=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
      Offset=Offset0+WriteOffsetC*IntKind
      ALLOCATE(iWork(NumberOfCellsB))
      iw=1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)==Zero) THEN
              iWork(iw)=0
            ELSE
              iWork(iw)=1
            END IF
            iw=iw+1
          END DO
        END DO
      END DO
      CALL MPI_FILE_WRITE_AT(fh,Offset,iWork, &
                             NumberOfCellsB, &
                             MPI_INTEGER,Status,MPIErr)
      DEALLOCATE(iWork)
    END DO
    Offset0=Offset0+NumberOfCells*IntKind
  END IF

END SUBROUTINE gmvGeometry

SUBROUTINE gmvVelocityC(VecC,RhoCell)

  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(ScalarCell_T) :: RhoCell(:)

  INTEGER :: ix,iy,iz
  INTEGER :: Cell
  REAL(4), ALLOCATABLE :: xN(:),yN(:),zN(:)


  WRITE(*,*) 'in gmvVelocityC'
  WRITE(OutputUnit,'(a8,a8,a30)') nodes,fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
  WRITE(OutputUnit,'(a8,a8,a30)') cells,fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
  WRITE(OutputUnit,'(a8,a9,a30)') material,fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
  ALLOCATE(xN(NumberOfCells))
  ALLOCATE(yN(NumberOfCells))
  ALLOCATE(zN(NumberOfCells))
  WRITE(OutputUnit,'(a8,i8)') velocity,0
  Cell=1
  DO ibb=1,nb
    CALL Set(Floor(ibb))
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xN(Cell)=VecC(ibb)%Vec(uPosL)%c(ix,iy,iz,1) &
                  /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps)         
          yN(Cell)=VecC(ibb)%Vec(vPosL)%c(ix,iy,iz,1) &
                  /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps)         
          zN(Cell)=VecC(ibb)%Vec(wPosL)%c(ix,iy,iz,1) &
                  /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps)         
          Cell=Cell+1
        END DO
      END DO
    END DO
  END DO
  WRITE(OutputUnit,*) xN(1:NumberOfCells)
  WRITE(OutputUnit,*) yN(1:NumberOfCells)
  WRITE(OutputUnit,*) zN(1:NumberOfCells)
  DEALLOCATE(xN,yN,zN)

END SUBROUTINE gmvVelocityC

SUBROUTINE gmvVelocityFWrite(Vel,RhoCell)

  TYPE(VelocityFace_T) :: Vel(:)
  TYPE(ScalarCell_T) :: RhoCell(:)

  INTEGER :: ix,iy,iz
  INTEGER :: i,iw
  INTEGER :: NumberOfCellsB
  INTEGER :: Cell
  REAL(Real4Kind), ALLOCATABLE :: xN(:),yN(:),zN(:)
  CHARACTER*500 :: cWork
  INTEGER :: cWorkLen
  CHARACTER*4 :: cZero
  INTEGER :: iZero=0
  REAL(Real4Kind) :: lam,phi,Vx,Vy,Vz
  CHARACTER*100 :: OutputFileNameGMV
  INTEGER :: PosSlash,NumberSlash


  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a6,a1,a8,a80)') nodes,' ',fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
    WRITE(OutputUnit,'(a6,a1,a8,a80)') cells,' ',fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
    WRITE(OutputUnit,'(a8,a1,a8,a80)') polygons,' ',fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
    ALLOCATE(xN(NumberOfCells))
    ALLOCATE(yN(NumberOfCells))
    ALLOCATE(zN(NumberOfCells))
    WRITE(OutputUnit,'(a8,i8)') velocity,0
    Cell=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              Cell=Cell+1
              xN(Cell)=(Vel(ibb)%uF(ix-1,iy,iz)*FU(ix-1,iy,iz) &
                            +Vel(ibb)%uF(ix,iy,iz)*FU(ix,iy,iz))    &
                            /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps)         &
!                           /(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
                            /(Two*Max(FU(ix-1,iy,iz),FU(ix,iy,iz))+Eps)
              yN(Cell)=(Vel(ibb)%vF(ix,iy-1,iz)*FV(ix,iy-1,iz) &
                            +Vel(ibb)%vF(ix,iy,iz)*FV(ix,iy,iz))    &
                            /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps)         &
!                           /(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
                            /(Two*Max(FV(ix,iy-1,iz),FV(ix,iy,iz))+Eps)
              zN(Cell)=(Vel(ibb)%wF(ix,iy,iz-1)*FW(ix,iy,iz-1) &
                            +Vel(ibb)%wF(ix,iy,iz)*FW(ix,iy,iz))    &
                            /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps)         &
!                           /(FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)
                            /(Two*Max(FW(ix,iy,iz-1),FW(ix,iy,iz))+Eps)
              IF (OutputVelocity=='GCart') THEN
                lam=xP(ix-1)+Half*dx(ix)
                phi=yP(iy-1)+Half*dy(iy)
                Vx=SIN(phi)*COS(lam)*zN(Cell) &
                  +COS(phi)*SIN(lam)*yN(Cell) &
                  -         SIN(lam)*xN(Cell) 
                Vy=SIN(phi)*SIN(lam)*zN(Cell) &
                  +COS(phi)*SIN(lam)*yN(Cell) &
                  +         COS(lam)*xN(Cell) 
                Vz=COS(phi)         *zN(Cell) &
                  -SIN(phi)         *yN(Cell) 
                xN(Cell)=Vx
                yN(Cell)=Vy
                zN(Cell)=Vz
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO
    WRITE(OutputUnit,*) xN(1:Cell)
    WRITE(OutputUnit,*) yN(1:Cell)
    WRITE(OutputUnit,*) zN(1:Cell)
    DEALLOCATE(xN,yN,zN)
  ELSE
    OutputFileNameGMV=OutputFileName
    NumberSlash=0
    DO
      PosSlash=INDEX(OutputFileNameGMV,'/')
      IF (PosSlash>0) THEN
        OutputFileNameGMV=OutputFileNameGMV(PosSlash+1:)
        NumberSlash=NumberSlash+1
      ELSE
        EXIT
      END IF
    END DO
    DO i=1,NumberSlash
      OutputFileNameGMV='../'//OutputFileNameGMV
    END DO
    cZero=TRANSFER(iZero,cZero)
    cWorkLen=LEN(nodes)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(cells)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(polygons)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(velocity)+LEN(cZero)
    IF (MyId==0) THEN
      cWork(1:cWorkLen) = &
              nodes//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //cells//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //polygons//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //velocity//cZero
      CALL MPI_FILE_WRITE_AT(fh,Offset0,cWork(1:cWorkLen),cWorkLen, &
                      MPI_BYTE,Status,MPIErr)
    END IF
    Offset0=Offset0+cWorkLen
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      NumberOfCellsB=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
      ALLOCATE(xN(NumberOfCellsB))
      ALLOCATE(yN(NumberOfCellsB))
      ALLOCATE(zN(NumberOfCellsB))
      iw=0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              iw=iw+1
              xN(iw)=(Vel(ibLoc)%uF(ix-1,iy,iz)*FU(ix-1,iy,iz) &
                            +Vel(ibLoc)%uF(ix,iy,iz)*FU(ix,iy,iz))    &
                            /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)         &
                            /(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
!                           /(Two*Max(FU(ix-1,iy,iz),FU(ix,iy,iz))+Eps)
              yN(iw)=(Vel(ibLoc)%vF(ix,iy-1,iz)*FV(ix,iy-1,iz) &
                            +Vel(ibLoc)%vF(ix,iy,iz)*FV(ix,iy,iz))    &
                            /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)         &
                            /(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
!                           /(Two*Max(FV(ix,iy-1,iz),FV(ix,iy,iz))+Eps)
              zN(iw)=(Vel(ibLoc)%wF(ix,iy,iz-1)*FW(ix,iy,iz-1) &
                            +Vel(ibLoc)%wF(ix,iy,iz)*FW(ix,iy,iz))    &
                            /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)         &
                            /(FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)
!                           /(Two*Max(FW(ix,iy,iz-1),FW(ix,iy,iz))+Eps)
              IF (OutputVelocity=='GCart') THEN
                lam=xP(ix-1)+Half*dx(ix)
                phi=yP(iy-1)+Half*dy(iy)
                Vx=SIN(phi)*COS(lam)*zN(iw) &
                  +COS(phi)*SIN(lam)*yN(iw) &
                  -         SIN(lam)*xN(iw)
                Vy=SIN(phi)*SIN(lam)*zN(iw) &
                  +COS(phi)*SIN(lam)*yN(iw) &
                  +         COS(lam)*xN(iw)
                Vz=COS(phi)         *zN(iw) &
                  -SIN(phi)         *yN(iw)
                xN(iw)=Vx
                yN(iw)=Vy
                zN(iw)=Vz
              END IF
              CALL native_4byte_real(xN(iw),xN(iw))
              CALL native_4byte_real(yN(iw),yN(iw))
              CALL native_4byte_real(zN(iw),zN(iw))
            ELSE IF (FillValue) THEN
              iW=iW+1
              xN(iW)=Zero
              yN(iW)=Zero
              zN(iW)=Zero
              CALL native_4byte_real(xN(iw),xN(iw))
              CALL native_4byte_real(yN(iw),yN(iw))
              CALL native_4byte_real(zN(iw),zN(iw))
            END IF
          END DO
        END DO
      END DO
      NumberOfCellsB=iw
      Offset=Offset0+Real4Kind*WriteOffsetCgmv
      CALL MPI_FILE_WRITE_AT(fh,Offset,xN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      Offset=Offset0+Real4Kind*WriteOffsetCgmv+Real4Kind*OffSetShiftGMV
      CALL MPI_FILE_WRITE_AT(fh,Offset,yN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      Offset=Offset0+Real4Kind*WriteOffsetCgmv+2*Real4Kind*OffSetShiftGMV
      CALL MPI_FILE_WRITE_AT(fh,Offset,zN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      DEALLOCATE(xN)
      DEALLOCATE(yN)
      DEALLOCATE(zN)
    END DO
    Offset0=Offset0+3*Real4Kind*OffSetShiftGMV
  END IF
  
END SUBROUTINE gmvVelocityFWrite

SUBROUTINE gmvVelocityCWrite(VecC)

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: ix,iy,iz
  INTEGER :: i,iw
  INTEGER :: NumberOfCellsB
  INTEGER :: Cell
  REAL(Real4Kind), ALLOCATABLE :: xN(:),yN(:),zN(:)
  CHARACTER*500 :: cWork
  INTEGER :: cWorkLen
  CHARACTER*4 :: cZero
  INTEGER :: iZero=0
  REAL(Real4Kind) :: lam,phi,Vx,Vy,Vz
  CHARACTER*100 :: OutputFileNameGMV
  INTEGER :: PosSlash,NumberSlash


  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a6,a1,a8,a80)') nodes,' ',fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
    WRITE(OutputUnit,'(a6,a1,a8,a80)') cells,' ',fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
    WRITE(OutputUnit,'(a8,a1,a8,a80)') polygons,' ',fromfile,TRIM('"'//TRIM(OutputFileName)//'.gmvG'//'"')
    ALLOCATE(xN(NumberOfCells))
    ALLOCATE(yN(NumberOfCells))
    ALLOCATE(zN(NumberOfCells))
    WRITE(OutputUnit,'(a8,i8)') velocity,0
    Cell=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              Cell=Cell+1
              xN(Cell)=VecC(ibb)%Vec(uPosL)%c(ix,iy,iz,1)
              yN(Cell)=VecC(ibb)%Vec(vPosL)%c(ix,iy,iz,1)
              zN(Cell)=VecC(ibb)%Vec(wPosL)%c(ix,iy,iz,1)
              IF (OutputVelocity=='GCart') THEN
                lam=xP(ix-1)+Half*dx(ix)
                phi=yP(iy-1)+Half*dy(iy)
                Vx=SIN(phi)*COS(lam)*zN(Cell) &
                  +COS(phi)*SIN(lam)*yN(Cell) &
                  -         SIN(lam)*xN(Cell) 
                Vy=SIN(phi)*SIN(lam)*zN(Cell) &
                  +COS(phi)*SIN(lam)*yN(Cell) &
                  +         COS(lam)*xN(Cell) 
                Vz=COS(phi)         *zN(Cell) &
                  -SIN(phi)         *yN(Cell) 
                xN(Cell)=Vx
                yN(Cell)=Vy
                zN(Cell)=Vz
              END IF
            END IF
          END DO
        END DO
      END DO
    END DO
    WRITE(OutputUnit,*) xN(1:Cell)
    WRITE(OutputUnit,*) yN(1:Cell)
    WRITE(OutputUnit,*) zN(1:Cell)
    DEALLOCATE(xN,yN,zN)
  ELSE
    OutputFileNameGMV=OutputFileName
    NumberSlash=0
    DO
      PosSlash=INDEX(OutputFileNameGMV,'/')
      IF (PosSlash>0) THEN
        OutputFileNameGMV=OutputFileNameGMV(PosSlash+1:)
        NumberSlash=NumberSlash+1
      ELSE
        EXIT
      END IF
    END DO
    DO i=1,NumberSlash
      OutputFileNameGMV='../'//OutputFileNameGMV
    END DO
    cZero=TRANSFER(iZero,cZero)
    cWorkLen=LEN(nodes)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(cells)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(polygons)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(velocity)+LEN(cZero)
    IF (MyId==0) THEN
      cWork(1:cWorkLen) = &
              nodes//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //cells//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //polygons//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //velocity//cZero
      CALL MPI_FILE_WRITE_AT(fh,Offset0,cWork(1:cWorkLen),cWorkLen, &
                      MPI_BYTE,Status,MPIErr)
    END IF
    Offset0=Offset0+cWorkLen
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      NumberOfCellsB=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
      ALLOCATE(xN(NumberOfCellsB))
      ALLOCATE(yN(NumberOfCellsB))
      ALLOCATE(zN(NumberOfCellsB))
      iw=0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              iw=iw+1
              xN(iw)=VecC(ibLoc)%Vec(uPosL)%c(ix,iy,iz,1)
              yN(iw)=VecC(ibLoc)%Vec(vPosL)%c(ix,iy,iz,1)
              zN(iw)=VecC(ibLoc)%Vec(wPosL)%c(ix,iy,iz,1)
              IF (OutputVelocity=='GCart') THEN
                lam=xP(ix-1)+Half*dx(ix)
                phi=yP(iy-1)+Half*dy(iy)
                Vx=SIN(phi)*COS(lam)*zN(iw) &
                  +COS(phi)*SIN(lam)*yN(iw) &
                  -         SIN(lam)*xN(iw)
                Vy=SIN(phi)*SIN(lam)*zN(iw) &
                  +COS(phi)*SIN(lam)*yN(iw) &
                  +         COS(lam)*xN(iw)
                Vz=COS(phi)         *zN(iw) &
                  -SIN(phi)         *yN(iw)
                xN(iw)=Vx
                yN(iw)=Vy
                zN(iw)=Vz
              END IF
              CALL native_4byte_real(xN(iw),xN(iw))
              CALL native_4byte_real(yN(iw),yN(iw))
              CALL native_4byte_real(zN(iw),zN(iw))
            END IF
          END DO
        END DO
      END DO
      NumberOfCellsB=iw
      Offset=Offset0+Real4Kind*WriteOffsetCgmv
      CALL MPI_FILE_WRITE_AT(fh,Offset,xN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      Offset=Offset0+Real4Kind*WriteOffsetCgmv+Real4Kind*OffSetShiftGMV
      CALL MPI_FILE_WRITE_AT(fh,Offset,yN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      Offset=Offset0+Real4Kind*WriteOffsetCgmv+2*Real4Kind*OffSetShiftGMV
      CALL MPI_FILE_WRITE_AT(fh,Offset,zN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      DEALLOCATE(xN)
      DEALLOCATE(yN)
      DEALLOCATE(zN)
    END DO
    Offset0=Offset0+3*Real4Kind*OffSetShiftGMV
  END IF
  
END SUBROUTINE gmvVelocityCWrite

SUBROUTINE gmvVelocityFRead(VecC)

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: ix,iy,iz
  INTEGER :: i,iw
  INTEGER :: NumberOfCellsB
  INTEGER :: Cell
  REAL(Real4Kind), ALLOCATABLE :: xN(:),yN(:),zN(:)
  CHARACTER*500 :: cWork
  INTEGER :: cWorkLen
  CHARACTER*4 :: cZero
  INTEGER :: iZero=0
  CHARACTER*30 :: FileName
  CHARACTER*8 :: Dummy
  INTEGER :: iTemp
  CHARACTER*100 :: OutputFileNameGMV
  INTEGER :: PosSlash,NumberSlash


  IF (gmvOutputType==ascii) THEN
    READ(OutputUnit,'(a8,a8,a30)') Dummy,Dummy,FileName
    READ(OutputUnit,'(a8,a8,a30)') Dummy,Dummy,FileName
    READ(OutputUnit,'(a8,a8,a30)') Dummy,Dummy,FileName
    ALLOCATE(xN(NumberOfCells))
    ALLOCATE(yN(NumberOfCells))
    ALLOCATE(zN(NumberOfCells))
    READ(OutputUnit,'(a8,i8)') Dummy,iTemp
    Cell=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              Cell=Cell+1
            END IF
          END DO
        END DO
      END DO
    END DO
    WRITE(OutputUnit,*) xN(1:Cell)
    WRITE(OutputUnit,*) yN(1:Cell)
    WRITE(OutputUnit,*) zN(1:Cell)
    Cell=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              Cell=Cell+1
              VecC(ibLoc)%Vec(uPosL)%c(ix,iy,iz,1)=xN(Cell)
              VecC(ibLoc)%Vec(vPosL)%c(ix,iy,iz,1)=yN(Cell)
              VecC(ibLoc)%Vec(wPosL)%c(ix,iy,iz,1)=zN(Cell)
            ELSE
              VecC(ibLoc)%Vec(uPosL)%c(ix,iy,iz,1)=MissValue
              VecC(ibLoc)%Vec(vPosL)%c(ix,iy,iz,1)=MissValue
              VecC(ibLoc)%Vec(wPosL)%c(ix,iy,iz,1)=MissValue
            END IF
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(xN,yN,zN)
  ELSE
    OutputFileNameGMV=OutputFileName
    NumberSlash=0
    DO
      PosSlash=INDEX(OutputFileNameGMV,'/')
      IF (PosSlash>0) THEN
        OutputFileNameGMV=OutputFileNameGMV(PosSlash+1:)
        NumberSlash=NumberSlash+1
      ELSE
        EXIT
      END IF
    END DO
    DO i=1,NumberSlash
      OutputFileNameGMV='../'//OutputFileNameGMV
    END DO
    cZero=TRANSFER(iZero,cZero)
    cWorkLen=LEN(nodes)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(cells)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(polygons)+LEN(fromfile)+LEN('"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"') &
            +LEN(velocity)+LEN(cZero)
    IF (MyId==0) THEN
      cWork(1:cWorkLen) = &
              nodes//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //cells//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //polygons//fromfile//'"'//TRIM(OutputFileNameGMV)//'.gmvG'//'"' &
            //velocity//cZero
      CALL MPI_FILE_READ_AT(fh,Offset0,cWork(1:cWorkLen),cWorkLen, &
                      MPI_BYTE,Status,MPIErr)
    END IF
    Offset0=Offset0+cWorkLen
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      NumberOfCellsB=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
      ALLOCATE(xN(NumberOfCellsB))
      ALLOCATE(yN(NumberOfCellsB))
      ALLOCATE(zN(NumberOfCellsB))
      iw=0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              iw=iw+1
            END IF
          END DO
        END DO
      END DO
      NumberOfCellsB=iw
      Offset=Offset0+Real4Kind*WriteOffsetCgmv
      CALL MPI_FILE_READ_AT(fh,Offset,xN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      Offset=Offset0+Real4Kind*WriteOffsetCgmv+Real4Kind*OffSetShiftGMV
      CALL MPI_FILE_READ_AT(fh,Offset,yN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      Offset=Offset0+Real4Kind*WriteOffsetCgmv+2*Real4Kind*OffSetShiftGMV
      CALL MPI_FILE_READ_AT(fh,Offset,zN(1:NumberOfCellsB), &
                             NumberOfCellsB, &
                             MPI_REAL,Status,MPIErr)
      IF (uPosL>Zero) THEN
        iw=0
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1
              IF (VolC(ix,iy,iz)>Zero) THEN
                iw=iw+1
                Cell=Cell+1
                CALL native_4byte_real(xN(iw),xN(iw))
                VecC(ibLoc)%Vec(uPosL)%c(ix,iy,iz,1)=xN(iw)
                CALL native_4byte_real(yN(iw),yN(iw))
                VecC(ibLoc)%Vec(vPosL)%c(ix,iy,iz,1)=yN(iw)
                CALL native_4byte_real(zN(iw),zN(iw))
                VecC(ibLoc)%Vec(wPosL)%c(ix,iy,iz,1)=zN(iw)
              ELSE
                VecC(ibLoc)%Vec(uPosL)%c(ix,iy,iz,1)=MissValue
                VecC(ibLoc)%Vec(vPosL)%c(ix,iy,iz,1)=MissValue
                VecC(ibLoc)%Vec(wPosL)%c(ix,iy,iz,1)=MissValue
              END IF
            END DO
          END DO
        END DO
      END IF
      DEALLOCATE(xN)
      DEALLOCATE(yN)
      DEALLOCATE(zN)
    END DO
    Offset0=Offset0+3*Real4Kind*OffSetShiftGMV
  END IF
  
END SUBROUTINE gmvVelocityFRead


SUBROUTINE gmvOpenScalarWrite

  CHARACTER(8), PARAMETER :: variable='variable'

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a8)') variable
  ELSE
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_WRITE_AT(fh,Offset,variable,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    Offset0=Offset0+LEN(variable)
  END IF


END SUBROUTINE gmvOpenScalarWrite


SUBROUTINE gmvScalarVarSoilWrite

  CHARACTER(8), PARAMETER :: variable='variable'

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnitS,'(a8)') variable
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetSoil0
      CALL MPI_FILE_WRITE_AT(fh2,Offset,variable,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    OffsetSoil0=OffsetSoil0+LEN(variable)
  END IF

END SUBROUTINE gmvScalarVarSoilWrite


SUBROUTINE gmvScalarVarCutWrite

  CHARACTER(8), PARAMETER :: variable='variable'

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnitS,'(a8)') variable
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetCut0
      CALL MPI_FILE_WRITE_AT(fh3,Offset,variable,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    OffsetCut0=OffsetCut0+LEN(variable)
  END IF

END SUBROUTINE gmvScalarVarCutWrite


SUBROUTINE gmvScalarVarSoilRead

  CHARACTER(8), PARAMETER :: variable='variable'

  CHARACTER*8 :: Dummy

  IF (gmvOutputType==ascii) THEN
    READ(OutputUnit,'(a8)') Dummy
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetSoil0
      CALL MPI_FILE_READ_AT(fh2,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    OffsetSoil0=OffsetSoil0+LEN(variable)
  END IF

END SUBROUTINE gmvScalarVarSoilRead


SUBROUTINE gmvScalarVarCutRead

  CHARACTER(8), PARAMETER :: variable='variable'

  CHARACTER*8 :: Dummy

  IF (gmvOutputType==ascii) THEN
    READ(OutputUnit,'(a8)') Dummy
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetCut0
      CALL MPI_FILE_READ_AT(fh3,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    OffsetCut0=OffsetCut0+LEN(variable)
  END IF

END SUBROUTINE gmvScalarVarCutRead


SUBROUTINE gmvOpenScalarRead

  CHARACTER(8), PARAMETER :: variable='variable'
  CHARACTER*8 :: Dummy

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a8)') Dummy
  ELSE
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_READ_AT(fh,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF            
    Offset0=Offset0+LEN(variable)
  END IF
  
END SUBROUTINE gmvOpenScalarRead

SUBROUTINE gmvScalarWrite(Vec,iPos,name,it,RhoCell,Prof,VecC,Type,UnitKonv)

  TYPE(Vector4Cell_T) :: Vec(:)
  INTEGER :: iPos
  CHARACTER(32) :: name
  INTEGER :: it
  TYPE(ScalarCell_T), OPTIONAL :: RhoCell(:)
  TYPE(ScalarCell_T), OPTIONAL :: Prof(:)
  TYPE(Vector4Cell_T), OPTIONAL :: VecC(:)
  CHARACTER(*), OPTIONAL :: Type
  Real(RealKind), OPTIONAL :: UnitKonv

  INTEGER :: ix,iy,iz
  INTEGER :: NumberOfCellsB
  INTEGER :: Cell
  REAL(4), ALLOCATABLE :: xN(:)

  IF (gmvOutputType==ascii) THEN
    ALLOCATE(xN(NumberOfCells))
    WRITE(OutputUnit,'(a32,i8)') name,0
    Cell=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              Cell=Cell+1
              IF (PRESENT(VecC)) THEN   
                IF (PRESENT(Type)) THEN
                  IF (Type=='aH2O') THEN
                    xN(Cell)=Vec(ibb)%Vec(iPos)%c(ix,iy,iz,it) &
                                  /(VecC(ibb)%Vec(iWater)%c(ix,iy,iz,it)+Eps) &
                                  /MolMass(OutputC(iPos)%iPosC) 
                   ELSE IF (Type=='aNUMBER') THEN           
                     xN(Cell)=Vec(ibb)%Vec(iPos)%c(ix,iy,iz,it) &
                                  /(VecC(ibb)%Vec(iNC)%c(ix,iy,iz,it)+Eps) 
                   ELSE IF (Type=='KG_KG') THEN           
                     WRITE(*,*) 'KG_KG'
                     xN(Cell)=Vec(ibb)%Vec(iPos)%c(ix,iy,iz,it) &
                                  /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps)
                   END IF               
                ELSE   
                   xN(Cell)=Vec(ibb)%Vec(iPos)%c(ix,iy,iz,it) &
                                /(VecC(ibb)%Vec(iNC)%c(ix,iy,iz,it)+Eps) 
                END IF               
              ELSE IF (PRESENT(RhoCell)) THEN   
                IF (PRESENT(Prof)) THEN   
                  xN(Cell)=Vec(ibb)%Vec(iPos)%c(ix,iy,iz,it) &
                                /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps) &
                          -Prof(ibb)%c(ix,iy,iz,1)/(RhoProfG(ibb)%c(ix,iy,iz,1)+Eps) 
                ELSE
                  xN(Cell)=Vec(ibb)%Vec(iPos)%c(ix,iy,iz,it) &
                                /(RhoCell(ibb)%c(ix,iy,iz,1)+Eps)
                END IF
              ELSE
                IF (PRESENT(Prof)) THEN   
                  xN(Cell)=Vec(ibb)%Vec(iPos)%c(ix,iy,iz,it) &
                          -Prof(ibb)%c(ix,iy,iz,1)
                ELSE
                  xN(Cell)=Vec(ibb)%Vec(iPos)%c(ix,iy,iz,it)
                END IF
                IF (GasSpeciesC(iPos)==1) THEN 
                  xN(Cell)=xN(Cell)*MolMass(OutputC(iPos)%iPosC) &  !mol/l -> g/l bzw g/m3
                           *1.d9  ! Umrechnung in mug
                END IF
              END IF 
            END IF 
          END DO
        END DO
      END DO
    END DO
    WRITE(OutputUnit,*) xN(1:Cell)
    DEALLOCATE(xN)
  ELSE
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_WRITE_AT(fh,Offset,name(1:8),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh,Offset,name(9:16),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh,Offset,name(17:32),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh,Offset,0,1, &
                      MPI_INTEGER,Status,MPIErr)
    END IF
    Offset0=Offset0+LEN(name)+IntKind
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      nx=ix1-ix0
      ny=iy1-iy0
      nz=iz1-iz0
      NumberOfCellsB=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
      ALLOCATE(xN(NumberOfCellsB))
      Cell=0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              Cell=Cell+1
              IF (PRESENT(VecC)) THEN
                IF (PRESENT(Type)) THEN
                  IF (Type(1:3).EQ.'MOL') THEN
                    xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it) &
                             /MolMass(OutputC(iPos)%iPosC) 
                  ELSE
                    xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it)*UnitKonv
                  END IF
                  IF (Type(LEN_TRIM(Type):LEN_TRIM(Type)).EQ.'L') THEN
                    xN(Cell)=xN(Cell)/(VecC(ibLoc)%Vec(iWater)%c(ix,iy,iz,it)+Eps)
                  ELSE IF (Type(LEN_TRIM(Type)-3:LEN_TRIM(Type)).EQ.'PART') THEN
                    IF (VecC(ibLoc)%Vec(iNC)%c(ix,iy,iz,it)>NumberMin) THEN
                      xN(Cell)=xN(Cell)/(VecC(ibLoc)%Vec(iNC)%c(ix,iy,iz,it)+Eps)
                    ELSE
                      xN(Cell)=Zero
                    END IF  
                  ELSE IF (Type(LEN_TRIM(Type)-1:LEN_TRIM(Type)).EQ.'M3') THEN
                    xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it)
                  ELSE IF (Type(LEN_TRIM(Type)-1:LEN_TRIM(Type)).EQ.'KG') THEN
                    xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it)/(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)
                  ELSE
                    xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it)
                  END IF
                ELSE
                  IF (VecC(ibLoc)%Vec(iNC)%c(ix,iy,iz,it)>NumberMin) THEN
                    xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it) &
                              /(VecC(ibLoc)%Vec(iNC)%c(ix,iy,iz,it)+Eps)  
                  ELSE
                    xN(Cell)=Zero
                  END IF  
                END IF               
              ELSE IF (PRESENT(RhoCell)) THEN
                IF (PRESENT(Prof)) THEN
                  xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it) &
                                /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps) &
                          -Prof(ibLoc)%c(ix,iy,iz,1)/(RhoProfG(ibLoc)%c(ix,iy,iz,1)+Eps)
                ELSE
                  xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it) &
                                /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)
                END IF
              ELSE
                IF (PRESENT(Prof)) THEN
                  xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it) &
                          -Prof(ibLoc)%c(ix,iy,iz,1)
                ELSE
                  xN(Cell)=Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it)
                END IF
                IF (GasSpeciesC(iPos)==1) THEN 
                  xN(Cell)=xN(Cell)*MolMass(OutputC(iPos)%iPosC) &
                           *1.d9  ! Umrechnung in mug
                END IF
              END IF
              CALL native_4byte_real(xN(Cell),xN(Cell))
            ELSE IF (FillValue) THEN
              Cell=Cell+1
              xN(Cell)=Zero
              CALL native_4byte_real(xN(Cell),xN(Cell))
            END IF
          END DO
        END DO
      END DO
      NumberOfCellsB=Cell
      Offset=Offset0+WriteOffsetCgmv*IntKind
      CALL MPI_FILE_WRITE_AT(fh,Offset,xN, &
                             NumberOfCellsB, &
                             MPI_INTEGER,Status,MPIErr)
      DEALLOCATE(xN)
    END DO
    Offset0=Offset0+Real4Kind*OffSetShiftGMV
  END IF

END SUBROUTINE gmvScalarWrite

SUBROUTINE gmvScalarRead(Vec,iPos,name,it)

  TYPE(Vector4Cell_T) :: Vec(:)
  INTEGER :: iPos
  CHARACTER(32) :: name

  INTEGER :: ix,iy,iz
  INTEGER :: NumberOfCellsB
  INTEGER :: Cell
  REAL(4), ALLOCATABLE :: xN(:)
  INTEGER :: iTemp ,it,k,iPosMax,l

  IF (gmvOutputType==ascii) THEN
    ALLOCATE(xN(NumberOfCells))
    READ(OutputUnit,'(a32,i8)') name,iTemp
    Cell=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN
              Cell=Cell+1
            END IF 
          END DO
        END DO
      END DO
    END DO
    READ(OutputUnit,*) xN(1:Cell)
    Cell=0
    DO ibb=1,nb
      CALL Set(Floor(ibb))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            IF (VolC(ix,iy,iz)>Zero) THEN 
              Cell=Cell+1
              Vec(ibb)%Vec(iPos)%c(ix,iy,iz,LBOUND(Vec(ibb)%Vec(iPos)%c,4))=xN(Cell)
            END IF
          END DO
        END DO
      END DO
    END DO
    DEALLOCATE(xN)
  ELSE
        IF (MyId==0) THEN
          Offset=Offset0
          CALL MPI_FILE_READ_AT(fh,Offset,name(1:8),1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
          Offset=Offset+8
          CALL MPI_FILE_READ_AT(fh,Offset,name(9:16),1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
          Offset=Offset+8
          CALL MPI_FILE_READ_AT(fh,Offset,name(17:32),1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
          Offset=Offset+8
          CALL MPI_FILE_READ_AT(fh,Offset,iTemp,1, &
                      MPI_INTEGER,Status,MPIErr)
        END IF
        Offset0=Offset0+LEN(name)+IntKind
        DO ibLoc=1,nbLoc
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
          ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
          iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
          iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
          iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
          iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
          nx=ix1-ix0
          ny=iy1-iy0
          nz=iz1-iz0
          NumberOfCellsB=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
          ALLOCATE(xN(NumberOfCellsB))
          Cell=0
          DO iz=iz0+1,iz1
            DO iy=iy0+1,iy1
              DO ix=ix0+1,ix1
                IF (VolC(ix,iy,iz)>Zero) THEN
                  Cell=Cell+1
                END IF
              END DO
            END DO
          END DO
          NumberOfCellsB=Cell
          Offset=Offset0+WriteOffsetCgmv*IntKind
          CALL MPI_FILE_READ_AT(fh,Offset,xN, &
                                 NumberOfCellsB, &
                                 MPI_INTEGER,Status,MPIErr)
          Cell=0
          DO iz=iz0+1,iz1
            DO iy=iy0+1,iy1
              DO ix=ix0+1,ix1
                IF (VolC(ix,iy,iz)>Zero) THEN
                  Cell=Cell+1
                  CALL native_4byte_real(xN(Cell),xN(Cell))
                  Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it)=xN(Cell)
                ELSE
                  Vec(ibLoc)%Vec(iPos)%c(ix,iy,iz,it)=Zero
                END IF
              END DO
            END DO
          END DO
          DEALLOCATE(xN)
        END DO
        Offset0=Offset0+Real4Kind*OffSetShiftGMV
  END IF

END SUBROUTINE gmvScalarRead

SUBROUTINE gmvCloseScalar

  CHARACTER(8), PARAMETER :: endvars='endvars '


  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnit,'(a8)') endvars
  ELSE
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_WRITE_AT(fh,Offset,endvars,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    Offset0=Offset0+LEN(endvars)
  END IF

END SUBROUTINE gmvCloseScalar

SUBROUTINE gmvCloseReadScalar

  CHARACTER(8) :: Dummy
  CHARACTER(8), PARAMETER :: endvars='endvars '

  IF (gmvOutputType==ascii) THEN
  ELSE
    IF (MyId==0) THEN
      Offset=Offset0
      CALL MPI_FILE_READ_AT(fh,Offset,Dummy,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    Offset0=Offset0+LEN(endvars)
  END IF

END SUBROUTINE gmvCloseReadScalar

SUBROUTINE gmvCloseScalarVarSoil

  CHARACTER(8), PARAMETER :: endvars='endvars '


  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnitS,'(a8)') endvars
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetSoil0
      CALL MPI_FILE_WRITE_AT(fh2,Offset,endvars,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    OffsetSoil0=OffsetSoil0+LEN(endvars)
  END IF

END SUBROUTINE gmvCloseScalarVarSoil

 
SUBROUTINE gmvCloseScalarVarCut

  CHARACTER(8), PARAMETER :: endvars='endvars '


  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnitS,'(a8)') endvars
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetCut0
      CALL MPI_FILE_WRITE_AT(fh3,Offset,endvars,1, &
                      MPI_DOUBLE_PRECISION,Status,MPIErr)
    END IF
    OffsetCut0=OffsetCut0+LEN(endvars)
  END IF

END SUBROUTINE gmvCloseScalarVarCut


SUBROUTINE gmvIncludeSoilWrite(FileName)
  CHARACTER(*) :: FileName
  CHARACTER*500 :: cWork
  INTEGER :: cWorkLen
  CHARACTER*100 :: OutputNameSoilGMV
  INTEGER ::  PosSlash,NumberSlash,i

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnitS,'(a6,a1,a8,a80)') nodes,   ' ',fromfile,TRIM('"'//TRIM(OutputNameSoil)//'.out'//'.gmvG'//'"')
    WRITE(OutputUnitS,'(a6,a1,a8,a80)') cells,   ' ',fromfile,TRIM('"'//TRIM(OutputNameSoil)//'.out'//'.gmvG'//'"')
    WRITE(OutputUnitS,'(a8,a1,a8,a80)') polygons,' ',fromfile,TRIM('"'//TRIM(OutputNameSoil)//'.out'//'.gmvG'//'"')
  ELSE
    OutputNameSoilGMV=OutputNameSoil
    NumberSlash=0
    DO
      PosSlash=INDEX(OutputNameSoilGMV,'/')
      IF (PosSlash>0) THEN
        OutputNameSoilGMV=OutputNameSoilGMV(PosSlash+1:)
        NumberSlash=NumberSlash+1
      ELSE
        EXIT
      END IF
    END DO
    DO i=1,NumberSlash
      OutputNameSoilGMV='../'//OutputNameSoilGMV
    END DO
    cWorkLen=LEN(nodes)+LEN(fromfile)+LEN('"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"') &
            +LEN(cells)+LEN(fromfile)+LEN('"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"') &
            +LEN(polygons)+LEN(fromfile)+LEN('"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"')
    IF (MyId==0) THEN
      cWork(1:cWorkLen) = &
              nodes//fromfile//'"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"' &
            //cells//fromfile//'"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"' &
            //polygons//fromfile//'"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"'
      CALL MPI_FILE_WRITE_AT(fh2,OffsetSoil0,cWork(1:cWorkLen),cWorkLen, &
                      MPI_BYTE,Status,MPIErr)
    END IF
    OffsetSoil0=OffsetSoil0+cWorkLen
  END IF
END SUBROUTINE


SUBROUTINE gmvIncludeCutWrite(FileName)
  CHARACTER(*) :: FileName
  CHARACTER*500 :: cWork
  INTEGER :: cWorkLen
  CHARACTER*100 :: OutputNameCutGMV
  INTEGER ::  PosSlash,NumberSlash,i

  IF (gmvOutputType==ascii) THEN
    WRITE(OutputUnitS,'(a6,a1,a8,a80)') nodes,   ' ',fromfile,TRIM('"'//TRIM(OutputNameCut)//'.out'//'.gmvG'//'"')
    WRITE(OutputUnitS,'(a6,a1,a8,a80)') cells,   ' ',fromfile,TRIM('"'//TRIM(OutputNameCut)//'.out'//'.gmvG'//'"')
    WRITE(OutputUnitS,'(a8,a1,a8,a80)') polygons,' ',fromfile,TRIM('"'//TRIM(OutputNameCut)//'.out'//'.gmvG'//'"')
  ELSE
    OutputNameCutGMV=OutputNameCut
    NumberSlash=0
    DO
      PosSlash=INDEX(OutputNameCutGMV,'/')
      IF (PosSlash>0) THEN
        OutputNameCutGMV=OutputNameCutGMV(PosSlash+1:)
        NumberSlash=NumberSlash+1
      ELSE
        EXIT
      END IF
    END DO
    DO i=1,NumberSlash
      OutputNameCutGMV='../'//OutputNameCutGMV
    END DO
    cWorkLen=LEN(nodes)+LEN(fromfile)+LEN('"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"') &
            +LEN(cells)+LEN(fromfile)+LEN('"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"') &
            +LEN(polygons)+LEN(fromfile)+LEN('"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"')
    IF (MyId==0) THEN
      cWork(1:cWorkLen) = &
              nodes//fromfile//'"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"' &
            //cells//fromfile//'"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"' &
            //polygons//fromfile//'"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"'
      CALL MPI_FILE_WRITE_AT(fh3,OffsetCut0,cWork(1:cWorkLen),cWorkLen, &
                      MPI_BYTE,Status,MPIErr)
    END IF
    OffsetCut0=OffsetCut0+cWorkLen
  END IF
END SUBROUTINE


SUBROUTINE gmvIncludeSoilRead(Filename)
  CHARACTER(*) :: FileName
  CHARACTER*8 :: Dummy
  CHARACTER*500 :: cWork
  INTEGER :: cWorkLen
  CHARACTER*100 :: OutputNameSoilGMV
  INTEGER ::  PosSlash,NumberSlash,i


  IF (gmvOutputType==ascii) THEN
    READ(OutputUnitS,'(a6,a1,a8,a80)') Dummy,Dummy,Dummy,Filename
    READ(OutputUnitS,'(a6,a1,a8,a80)') Dummy,Dummy,Dummy,Filename
    READ(OutputUnitS,'(a8,a1,a8,a80)') Dummy,Dummy,Dummy,Filename
  ELSE
    OutputNameSoilGMV=OutputNameSoil
    NumberSlash=0
    DO
      PosSlash=INDEX(OutputNameSoilGMV,'/')
      IF (PosSlash>0) THEN
        OutputNameSoilGMV=OutputNameSoilGMV(PosSlash+1:)
        NumberSlash=NumberSlash+1
      ELSE
        EXIT
      END IF
    END DO
    DO i=1,NumberSlash
      OutputNameSoilGMV='../'//OutputNameSoilGMV
    END DO
    cWorkLen=LEN(nodes)+LEN(fromfile)+LEN('"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"') &
            +LEN(cells)+LEN(fromfile)+LEN('"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"') &
            +LEN(polygons)+LEN(fromfile)+LEN('"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"')
    IF (MyId==0) THEN
      cWork(1:cWorkLen) = &
              nodes//fromfile//'"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"' &
            //cells//fromfile//'"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"' &
            //polygons//fromfile//'"'//TRIM(OutputNameSoilGMV)//'.out'//'.gmvG'//'"'
      CALL MPI_FILE_READ_AT(fh2,OffsetSoil0,cWork(1:cWorkLen),cWorkLen, &
                      MPI_BYTE,Status,MPIErr)
    END IF
    OffsetSoil0=OffsetSoil0+cWorkLen
  END IF
END SUBROUTINE gmvIncludeSoilRead


SUBROUTINE gmvIncludeCutRead(Filename)
  CHARACTER(*) :: FileName
  CHARACTER*8 :: Dummy
  CHARACTER*500 :: cWork
  INTEGER :: cWorkLen
  CHARACTER*100 :: OutputNameCutGMV
  INTEGER ::  PosSlash,NumberSlash,i


  IF (gmvOutputType==ascii) THEN
    READ(OutputUnitS,'(a6,a1,a8,a80)') Dummy,Dummy,Dummy,Filename
    READ(OutputUnitS,'(a6,a1,a8,a80)') Dummy,Dummy,Dummy,Filename
    READ(OutputUnitS,'(a8,a1,a8,a80)') Dummy,Dummy,Dummy,Filename
  ELSE
    OutputNameCutGMV=OutputNameCut
    NumberSlash=0
    DO
      PosSlash=INDEX(OutputNameCutGMV,'/')
      IF (PosSlash>0) THEN
        OutputNameCutGMV=OutputNameCutGMV(PosSlash+1:)
        NumberSlash=NumberSlash+1
      ELSE
        EXIT
      END IF
    END DO
    DO i=1,NumberSlash
      OutputNameCutGMV='../'//OutputNameCutGMV
    END DO
    cWorkLen=LEN(nodes)+LEN(fromfile)+LEN('"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"') &
            +LEN(cells)+LEN(fromfile)+LEN('"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"') &
            +LEN(polygons)+LEN(fromfile)+LEN('"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"')
    IF (MyId==0) THEN
      cWork(1:cWorkLen) = &
              nodes//fromfile//'"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"' &
            //cells//fromfile//'"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"' &
            //polygons//fromfile//'"'//TRIM(OutputNameCutGMV)//'.out'//'.gmvG'//'"'
      CALL MPI_FILE_READ_AT(fh3,OffsetCut0,cWork(1:cWorkLen),cWorkLen, &
                      MPI_BYTE,Status,MPIErr)
    END IF
    OffsetCut0=OffsetCut0+cWorkLen
  END IF
END SUBROUTINE gmvIncludeCutRead


SUBROUTINE gmvSoilVecWrite(Vec,iPos,name,RhoCell,Prof)
  TYPE(Vector4Cell_T) :: Vec(:)
  INTEGER :: iPos
  CHARACTER(32) :: name
  TYPE(ScalarCell_T), OPTIONAL :: RhoCell(:)
  TYPE(ScalarCell_T), OPTIONAL :: Prof(:)

  INTEGER :: ix,iy,iz,i
  INTEGER :: NumberOfCellsB,NumBoundsAllSoil
  INTEGER :: iCell,iS
  REAL(4), ALLOCATABLE :: xN(:)
  REAL(RealKind), POINTER :: cB(:,:)
  INTEGER :: iTemp

! WRITE (*,*) 'gmvSoilVecWrite: iPos,name',iPos,name
  iTemp=0
  NumBoundsAllSoil=0
  DO ib=1,nb
      NumBoundsAllSoil=NumBoundsAllSoil+Floor(ib)%NumBoundCell*Domain%nrsoillayers !nzS
  END DO    
!  Write(*,*) "NumBoundsAllSoil=",NumBoundsAllSoil, "   aus gmvSoilVecWrite"

  IF (gmvOutputType==ascii) THEN
    ALLOCATE(xN(NumBoundsAllSoil))
    WRITE(OutputUnitS,'(a32,i8)') name,0
    iCell=0
    DO ib=1,nb
      CALL Set(Floor(ib))
      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers !nzS
          iCell=iCell+1
          xN(iCell)=cB(i,iS)
        END DO
      END DO
    END DO
    WRITE(OutputUnitS,*) xN(1:iCell)
    DEALLOCATE(xN)
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetSoil0
      CALL MPI_FILE_WRITE_AT(fh2,Offset,name(1:8),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh2,Offset,name(9:16),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh2,Offset,name(17:32),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh2,Offset,0,1, &
                      MPI_INTEGER,Status,MPIErr)
    END IF
    OffsetSoil0=OffsetSoil0+LEN(name)+IntKind
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      NumberOfCellsB=NumBoundCell*Domain%nrsoillayers !nzS
      iTemp=iTemp+NumberOfCellsB
      ALLOCATE(xN(NumberOfCellsB))
      iCell=0
      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers !nzS
          iCell=iCell+1
          xN(iCell)=cB(i,iS)
          CALL native_4byte_real(xN(iCell),xN(iCell))
        END DO
      END DO
      NumberOfCellsB=iCell
      Offset=OffsetSoil0+WriteOffsetCSoilgmv*IntKind
      CALL MPI_FILE_WRITE_AT(fh2,Offset,xN,NumberOfCellsB, &
                             MPI_INTEGER,Status,MPIErr)
      DEALLOCATE(xN)
    END DO
    OffsetSoil0=OffsetSoil0+Real4Kind*OffsetShiftSoilGMV
    !OffsetSoil0=OffsetSoil0+Real4Kind*NumBoundsAllSoil
  END IF
! Write(*,*) "NumBoundsAllSoil=",size(cB(:,:)),name

END SUBROUTINE gmvSoilVecWrite


SUBROUTINE gmvCutVecWrite(Vec,iPos,name,RhoCell,Prof)
  TYPE(Vector4Cell_T) :: Vec(:)
  INTEGER :: iPos
  CHARACTER(32) :: name
  TYPE(ScalarCell_T), OPTIONAL :: RhoCell(:)
  TYPE(ScalarCell_T), OPTIONAL :: Prof(:)

  INTEGER :: ix,iy,iz,i
  INTEGER :: NumberOfCellsC,NumBoundsAllCut
  INTEGER :: iCell,iS
  REAL(4), ALLOCATABLE :: xN(:)
  REAL(RealKind), POINTER :: cB(:,:)
  INTEGER :: iTemp

! WRITE (*,*) 'gmvCutVecWrite: iPos,name',iPos,name
  iTemp=0
  NumBoundsAllCut=0
  DO ib=1,nb
      NumBoundsAllCut=NumBoundsAllCut+Floor(ib)%NumBoundCell*Domain%nrsoillayers !nzS
  END DO    
! Write(*,*) "NumBoundsAllCut=",NumBoundsAllCut, "   aus gmvCutVecWrite"

  IF (gmvOutputType==ascii) THEN
    ALLOCATE(xN(NumBoundsAllCut))
    WRITE(OutputUnitS,'(a32,i8)') name,0
    iCell=0
    DO ib=1,nb
      CALL Set(Floor(ib))
      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        iCell=iCell+1
        xN(iCell)=cB(i,1)
      END DO
    END DO
    WRITE(OutputUnitS,*) xN(1:iCell)
    DEALLOCATE(xN)
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetCut0
      CALL MPI_FILE_WRITE_AT(fh3,Offset,name(1:8),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh3,Offset,name(9:16),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh3,Offset,name(17:32),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_WRITE_AT(fh3,Offset,0,1, &
                      MPI_INTEGER,Status,MPIErr)
    END IF
    OffsetCut0=OffsetCut0+LEN(name)+IntKind
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      NumberOfCellsC=NumBoundCell
      iTemp=iTemp+NumberOfCellsC
      ALLOCATE(xN(NumberOfCellsC))
      iCell=0
      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        iCell=iCell+1
        xN(iCell)=cB(i,1)
        CALL native_4byte_real(xN(iCell),xN(iCell))
      END DO
      NumberOfCellsC=iCell
      Offset=OffsetCut0+WriteOffsetCCutgmv*IntKind
      CALL MPI_FILE_WRITE_AT(fh3,Offset,xN,NumberOfCellsC, &
                             MPI_INTEGER,Status,MPIErr)
      DEALLOCATE(xN)
    END DO
    OffsetCut0=OffsetCut0+Real4Kind*OffsetShiftCutGMV
  END IF

END SUBROUTINE gmvCutVecWrite

SUBROUTINE gmvSoilVecRead(Vec,iPos,name,RhoCell,Prof)
  TYPE(Vector4Cell_T) :: Vec(:)
  INTEGER :: iPos
  CHARACTER(32) :: name
  CHARACTER*8 :: Dummy
  TYPE(ScalarCell_T), OPTIONAL :: RhoCell(:)
  TYPE(ScalarCell_T), OPTIONAL :: Prof(:)

  INTEGER :: ix,iy,iz,i
  INTEGER :: NumberOfCellsB,NumBoundsAllSoil
  INTEGER :: iCell,iS
  REAL(4), ALLOCATABLE :: xN(:)
  REAL(RealKind), POINTER :: cB(:,:)
  INTEGER :: iTemp

  NumBoundsAllSoil=0
  DO ib=1,nb
      NumBoundsAllSoil=NumBoundsAllSoil+Floor(ib)%NumBoundCell*Domain%nrsoillayers !nzS
  END DO   

  IF (gmvOutputType==ascii) THEN
    ALLOCATE(xN(NumBoundsAllSoil))
    READ(OutputUnitS,'(a32,i8)') name,iTemp !Dummy
    iCell=0
    DO ib=1,nb
      CALL Set(Floor(ib))
!      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers !nzS
          iCell=iCell+1
        END DO
      END DO
    END DO
    READ(OutputUnitS,*) xN(1:iCell)
    iCell=0
    DO ib=1,nb
      CALL Set(Floor(ib))
!      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers
          iCell=iCell+1
          Vec(ibLoc)%Vec(iPos)%cB(i,iS)=xN(iCell)
        END DO
      END DO
    END DO
    DEALLOCATE(xN)
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetSoil0
      CALL MPI_FILE_READ_AT(fh2,Offset,name(1:8),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh2,Offset,name(9:16),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh2,Offset,name(17:32),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh2,Offset,iTemp,1, &
                      MPI_INTEGER,Status,MPIErr)
    END IF
    OffsetSoil0=OffsetSoil0+LEN(name)+IntKind
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      NumberOfCellsB=NumBoundCell*Domain%nrsoillayers !nzS
      ALLOCATE(xN(NumberOfCellsB))
      iCell=0
      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers !nzS
          iCell=iCell+1
        END DO
      END DO
      NumberOfCellsB=iCell
      Offset=OffsetSoil0+WriteOffsetCSoilgmv*IntKind
      CALL MPI_FILE_READ_AT(fh2,Offset,xN,NumberOfCellsB, &
                             MPI_INTEGER,Status,MPIErr)
      iCell=0
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers !nzS
          iCell=iCell+1
          CALL native_4byte_real(xN(iCell),xN(iCell))
          Vec(ibLoc)%Vec(iPos)%cB(i,iS)=xN(iCell)
        END DO
      END DO

      DEALLOCATE(xN)
    END DO
    OffsetSoil0=OffsetSoil0+Real4Kind*OffsetShiftSoilGMV
    !OffsetSoil0=OffsetSoil0+Real4Kind*NumBoundsAllSoil
  END IF

END SUBROUTINE gmvSoilVecRead

SUBROUTINE gmvCutVecRead(Vec,iPos,name,RhoCell,Prof)
  TYPE(Vector4Cell_T) :: Vec(:)
  INTEGER :: iPos
  CHARACTER(32) :: name
  CHARACTER*8 :: Dummy
  TYPE(ScalarCell_T), OPTIONAL :: RhoCell(:)
  TYPE(ScalarCell_T), OPTIONAL :: Prof(:)

  INTEGER :: ix,iy,iz,i
  INTEGER :: NumberOfCellsB,NumBoundsAllCut
  INTEGER :: iCell,iS
  REAL(4), ALLOCATABLE :: xN(:)
  REAL(RealKind), POINTER :: cB(:,:)
  INTEGER :: iTemp

  NumBoundsAllCut=0
  DO ib=1,nb
      NumBoundsAllCut=NumBoundsAllCut+Floor(ib)%NumBoundCell*Domain%nrsoillayers !nzS
  END DO   

  IF (gmvOutputType==ascii) THEN
    ALLOCATE(xN(NumBoundsAllCut))
    READ(OutputUnitS,'(a32,i8)') name,iTemp !Dummy
    iCell=0
    DO ib=1,nb
      CALL Set(Floor(ib))
!      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers !nzS
          iCell=iCell+1
        END DO
      END DO
    END DO
    READ(OutputUnitS,*) xN(1:iCell)
    iCell=0
    DO ib=1,nb
      CALL Set(Floor(ib))
!      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers
          iCell=iCell+1
          Vec(ibLoc)%Vec(iPos)%cB(i,iS)=xN(iCell)
        END DO
      END DO
    END DO
    DEALLOCATE(xN)
  ELSE
    IF (MyId==0) THEN
      Offset=OffsetCut0
      CALL MPI_FILE_READ_AT(fh3,Offset,name(1:8),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh3,Offset,name(9:16),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh3,Offset,name(17:32),1, &
                  MPI_DOUBLE_PRECISION,Status,MPIErr)
      Offset=Offset+8
      CALL MPI_FILE_READ_AT(fh3,Offset,iTemp,1, &
                      MPI_INTEGER,Status,MPIErr)
    END IF
    OffsetCut0=OffsetCut0+LEN(name)+IntKind
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      NumberOfCellsB=NumBoundCell*Domain%nrsoillayers !nzS
      ALLOCATE(xN(NumberOfCellsB))
      iCell=0
      cB=>Vec(ibLoc)%Vec(iPos)%cb
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers !nzS
          iCell=iCell+1
        END DO
      END DO
      NumberOfCellsB=iCell
      Offset=OffsetCut0+WriteOffsetCCutgmv*IntKind
      CALL MPI_FILE_READ_AT(fh3,Offset,xN,NumberOfCellsB, &
                             MPI_INTEGER,Status,MPIErr)
      iCell=0
      DO i=1,NumBoundCell
        DO iS=1,Domain%nrsoillayers !nzS
          iCell=iCell+1
          CALL native_4byte_real(xN(iCell),xN(iCell))
          Vec(ibLoc)%Vec(iPos)%cB(i,iS)=xN(iCell)
        END DO
      END DO

      DEALLOCATE(xN)
    END DO
    OffsetCut0=OffsetCut0+Real4Kind*OffsetShiftCutGMV
  END IF

END SUBROUTINE gmvCutVecRead

SUBROUTINE InputModelOutput(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: Pos
  CHARACTER(300) :: Line

  VelOut=.FALSE.
  thOut=.FALSE.
  EnOut=.FALSE.
  tkeOut=.FALSE.
  disOut=.FALSE.
  omeOut=.FALSE.
  tkeHOut=.FALSE.
  tkeVOut=.FALSE.
  LenOut=.FALSE.
  RhoVOut=.FALSE.
  RhoCOut=.FALSE.
  RhoROut=.FALSE.
  RhoIOut=.FALSE.
  RhoSOut=.FALSE.
  nvOut=.FALSE.
  ncOut=.FALSE.
  nrOut=.FALSE.
  niOut=.FALSE.
  nsOut=.FALSE.
  tracer1Out=.FALSE.
  tracer2Out=.FALSE.
  RhoOut=.FALSE.
  PreOut=.FALSE.
  DampOut=.FALSE.
  ShadOut=.FALSE.
  RadiationOut=.FALSE.
  FluxOut=.FALSE.
  LandClassOut=.FALSE.
  DiffOut=.FALSE.
  DiffHOut=.FALSE.
  DiffVOut=.FALSE.
  DiffPotOut=.FALSE.
  DiffPotHOut=.FALSE.
  DiffPotVOut=.FALSE.
  DiffMomOut=.FALSE.
  DiffMomHOut=.FALSE. !marcelk
  DiffMomVOut=.FALSE. !marcelk
  ChemieOut=.FALSE.
  ProfOut=.FALSE.
  OutputTimeStart=1.0e30_RealKind
  ProfIslandOcean=.FALSE.
  OutputTimeStep=Zero
  OutputTimeEnd=Zero
  OutputStart=0
  OutputFrequ=0
  OutputEnd=0
  OutputFileName='VelocGlob'
  OutputType(1:3)='AVS'
  OutputGrid(1:5)='Plane'
  OutputVelocity(1:5)=''
  OutputVarFileName='Output'

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')

   DO
     READ(InputUnit,*,END=1) Line
     IF (INDEX(Line,'#OutputDomain')>0) THEN
       READ(InputUnit,*) Line
       IF (Line=='XYZ_Number') THEN 
         READ(InputUnit,*) ix0Out,ix1Out,iy0Out,iy1Out,iz0Out,iz1Out
       ELSE IF (Line=='XYZ_Coord') THEN 
       END IF  
       EXIT
     END IF
   END DO

! OSWALD
! DO ibLoc=1,nbLoc
!   ib=LocGlob(ibLoc)
!   CALL Set(Floor(ib))
!   IF (ib==LocGlob(1)) xStart=xP(ix0Out)
!   IF (ib==LocGlob(nbLoc)) xEnd=xP(ix1Out)
!   IF (ib==LocGlob(1)) yStart=yp(iy0Out)
!   IF (ib==LocGlob(nbLoc)) yEnd=yP(iy1Out)
!   xProf=(xStart+xEnd)/2
!   yProf=(yStart+yEnd)/2
! END DO

  xStart=Domain%xP(ix0Out)
  xEnd=Domain%xP(ix1Out)
  yStart=Domain%yP(iy0Out)
  yEnd=Domain%yP(iy1Out)

  REWIND(InputUnit)
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ModelOutput')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelOutput)
      EXIT
    END IF
  END DO

  OutputType=ADJUSTL(OutputType)
  IF (uPosL*vPosL*wPosL==0) THEN
    VelOut=.FALSE.
  END IF
  IF (thPos==0) THEN
    thOut=.FALSE.
  END IF
  IF (EnPos==0) THEN
    EnOut=.FALSE.
  END IF
  IF (tkePos==0) THEN
    tkeOut=.FALSE.
  END IF
  IF (disPos==0) THEN
    disOut=.FALSE.
  END IF
  IF (omePos==0) THEN
    omeOut=.FALSE.
  END IF
  IF (tkeHPos==0) THEN
    tkeHOut=.FALSE.
  END IF
  IF (tkeVPos==0) THEN
    tkeVOut=.FALSE.
  END IF
  IF (LenPos==0) THEN
    LenOut=.FALSE.
  END IF
  IF (RhoVPos==0) THEN
    RhoVOut=.FALSE.
  END IF
  IF (RhoCPos==0) THEN
    RhoCOut=.FALSE.
  END IF
  IF (RhoRPos==0) THEN
    RhoROut=.FALSE.
  END IF
  IF (RhoIPos==0) THEN
    RhoIOut=.FALSE.
  END IF
  IF (RhoSPos==0) THEN
    RhoSOut=.FALSE.
  END IF
  IF (nvPos==0) THEN
    nvOut=.FALSE.
  END IF
  IF (ncPos==0) THEN
    ncOut=.FALSE.
  END IF
  IF (nrPos==0) THEN
    nrOut=.FALSE.
  END IF
  IF (niPos==0) THEN
    niOut=.FALSE.
  END IF
  IF (nsPos==0) THEN
    nsOut=.FALSE.
  END IF
  IF (tracer1Pos==0) THEN
    tracer1Out=.FALSE.
  END IF
  IF (tracer2Pos==0) THEN
    tracer2Out=.FALSE.
  END IF
  IF (RhoPos==0) THEN
    RhoOut=.FALSE.
  END IF
  IF (.NOT.Diffusion) THEN ! Hinneburg
    DiffOut=.FALSE.
    DiffHOut=.FALSE.
    DiffVOut=.FALSE.
    DiffPotOut=.FALSE.
    DiffMomOut=.FALSE.
    DiffPotHOut=.FALSE.
    DiffPotVOut=.FALSE.
    DiffMomHOut=.FALSE.
    DiffMomVOut=.FALSE.
  END IF

  CLOSE(UNIT=InputUnit)
1 CONTINUE

  CALL InputOutputChemie(OutputVarFileName)
  CALL InputOutputSpecial(OutputVarFileName)

  !IF (MyID==0) THEN
  !  WRITE(*,*) 'Output region:           ',ix0Out,ix1Out,' ',iy0Out,iy1Out &
  !                                    ,' ',iz0Out,iz1Out
  !  WRITE(*,*)
  !END IF

  OutputStep=OutputStart
  Step=OutputStart
  OutputTime=OutputTimeStart
END SUBROUTINE InputModelOutput

SUBROUTINE OutputInit(VecC,VecT)

  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)
    
!---  local variables
    
  INTEGER :: i,NumOut,nbOut
  INTEGER :: ix,iy,iz
  INTEGER :: ic,icT
  INTEGER :: icB
  CHARACTER(2) :: nFracChar

  NumData=0
  NumDataB=0
  NumDataC=0
  VectorComponentsOut=0
  VectorComponentsOutB=0
  VectorComponentsOutC=0
  IF (DynamicSoil) THEN
    IF (thOut.AND.thPos>0) THEN
      VectorComponentsOutB=VectorComponentsOutB+1
    END IF
    IF (RhoVOut.AND.RhoVPos>0) THEN
      VectorComponentsOutB=VectorComponentsOutB+1
    END IF
  END IF
  IF (ShadOut.AND.Radiation) THEN
    VectorComponentsOutC=VectorComponentsOutC+1
  END IF
  IF (RadiationOut.AND.Radiation) THEN
    VectorComponentsOutC=VectorComponentsOutC+3
  END IF
  IF (FluxOut.AND.Radiation) THEN
    VectorComponentsOutC=VectorComponentsOutC+2
  END IF
  IF (LandClassOut.AND.Radiation) THEN
    VectorComponentsOutC=VectorComponentsOutC+6
  END IF
  IF (VelOut) THEN
    NumData=NumData+3
  END IF
  IF (thOut.AND.thPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (EnOut.AND.EnPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (tkeOut.AND.tkePos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (disOut.AND.disPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (omeOut.AND.omePos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (tkeHOut.AND.tkeHPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (tkeVOut.AND.tkeVPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (LenOut.AND.LenPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoVOut.AND.RhoVPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoCOut.AND.RhoCPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoROut.AND.RhoRPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoIOut.AND.RhoIPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoSOut.AND.RhoSPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (nvOut.AND.nvPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (ncOut.AND.ncPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (nrOut.AND.nrPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (niOut.AND.niPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (nsOut.AND.nsPos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (tracer1Out.AND.tracer1Pos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (tracer2Out.AND.tracer2Pos>0) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (PreOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DampOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffHOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffVOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffPotOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffPotHOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffPotVOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffMomOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffMomHOut) THEN !marcelk
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffMomVOut) THEN !marcelk
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (ChemieOut) THEN
    DO icT=1,NumAeroOut
      NumData=NumData+nFrac
      VectorComponentsOut=VectorComponentsOut+1
    END DO
    DO icT=1,NumGasOut
      NumData=NumData+1
      VectorComponentsOut=VectorComponentsOut+1
    END DO
  END IF
  DO icT=1,LenOutSpecial
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END DO
  
  ALLOCATE(VectorOut(nbLoc))
  ALLOCATE(OutputC(VectorComponentsOut))
! ALLOCATE(DiffRho(VectorComponentsOut))
! ALLOCATE(NameScalars(VectorComponentsOut))
  ALLOCATE(NameScalarsB(VectorComponentsOutB+VectorComponentsOutC))
  ALLOCATE(GasSpeciesC(VectorComponentsOut))
  GasSpeciesC(1:VectorComponentsOut)=0
  ic=0
  icB=0
  IF (DynamicSoil) THEN
    IF (thOut.AND.thPos>0) THEN
      icB=icB+1
      NameScalarsB(icB)='SoilTemp '
    END IF
    IF (RhoVOut.AND.RhoVPos>0) THEN
      icB=icB+1
      NameScalarsB(icB)='SoilWetness'
    END IF
  END IF
  IF (ShadOut.AND.Radiation) THEN
    icB=icB+1
    NameScalarsB(icB)='Shadow'
  END IF
  IF (RadiationOut.AND.Radiation) THEN
    icB=icB+1
    NameScalarsB(icB)='RadDirekt'
    icB=icB+1
    NameScalarsB(icB)='RadDiffus'
    icB=icB+1
    NameScalarsB(icB)='RadInfred'
  END IF
  IF (FluxOut.AND.Radiation) THEN
    icB=icB+1
    NameScalarsB(icB)='SensHeatFlux'
    icB=icB+1
    NameScalarsB(icB)='LatHeatFlux'
  END IF
  IF (LandClassOut.AND.Radiation) THEN
    icB=icB+1
    NameScalarsB(icB)='LandClass'
    icB=icB+1
    NameScalarsB(icB)='z0'
    icB=icB+1
    NameScalarsB(icB)='Albedo'
    icB=icB+1
    NameScalarsB(icB)='C_D'
    icB=icB+1
    NameScalarsB(icB)='C_H'
    icB=icB+1
    NameScalarsB(icB)='C_E'
  END IF
  IF (thOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='PotTemp '
  END IF
  IF (EnOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='TotEn '
  END IF
  IF (tkeOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='Tke     '
  END IF
  IF (disOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='Diss    '
  END IF
  IF (omeOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='Ome     '
  END IF
  IF (tkeHOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='TkeH    '
  END IF
  IF (tkeVOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='TkeV    '
  END IF
  IF (LenOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='Len     '
  END IF    
  IF (RhoVOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='RhoV      '
  END IF
  IF (RhoCOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='RhoC      '
  END IF
  IF (RhoROut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='RhoR      '
  END IF
  IF (RhoIOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='RhoI      '
  END IF
  IF (RhoSOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='RhoS      '
  END IF
  IF (nvOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='nv      '
  END IF
  IF (ncOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='nc      '
  END IF
  IF (nrOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='nr      '
  END IF
  IF (niOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='ni      '
  END IF
  IF (nsOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='ns      '
  END IF
  IF (tracer1Out) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='Tracer1 '
  END IF
  IF (tracer2Out) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='Tracer2 '
  END IF
  IF (RhoOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='Rho     '
  END IF
  IF (PreOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Unsc'
    OutputC(ic)%NameScalar='Pre     '
  END IF
  IF (DampOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DampKoef'
  END IF
  IF (DiffOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffKoef'
  END IF
  IF (DiffHOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffHKoef'
  END IF
  IF (DiffVOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffVKoef'
  END IF
  IF (DiffPotOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffPotKoef'
  END IF
  IF (DiffPotHOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffPotHKoef'
  END IF
  IF (DiffPotVOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffPotVKoef'
  END IF
  IF (DiffMomOut) THEN
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffMomKoef'
  END IF
  IF (DiffMomHOut) THEN !marcelk
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffMomHKoef'
  END IF
  IF (DiffMomVOut) THEN !marcelk
    ic=ic+1
    OutputC(ic)%DiffRho='Rho'
    OutputC(ic)%NameScalar='DiffMomVKoef'
  END IF
  IF (ChemieOut) THEN
    DO icT=1,NumAeroOut
      ic=ic+1
      OutputC(ic)%Type='Aerosol'
      IF (SpeciesName(AeroOut(icT))=='aNUMBER') THEN
        OutputC(ic)%DiffRho='Unsc'
      ELSE IF (SpeciesName(AeroOut(icT))=='aRELAX') THEN
        OutputC(ic)%DiffRho='Unsc'
      ELSE   
        OutputC(ic)%DiffRho=AeroOutUnit(icT)
        OutputC(ic)%UnitKonvFactor=AeroUnitKonvFactor(icT)
      END IF  
      OutputC(ic)%NameScalar=SpeciesName(AeroOut(icT))
    END DO
    DO icT=1,NumGasOut
      ic=ic+1
      GasSpeciesC(ic)=2
      OutputC(ic)%DiffRho=GasOutUnit(icT)
      OutputC(ic)%UnitKonvFactor=GasUnitKonvFactor(icT)
      OutputC(ic)%NameScalar=SpeciesName(GasOut(icT))
    END DO
  END IF
  DO icT=1,LenOutSpecial
    ic=ic+1
    OutputC(ic)%Type=SpecialType(icT)
    OutputC(ic)%DiffRho=UnitOutSpecial(icT)
    OutputC(ic)%NameScalar=NameOutSpecial(icT)
    OutputC(ic)%UnitKonvFactor=UnitKonvSpecial(icT)
  END DO
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    ALLOCATE(VectorOut(ibLoc)%Vec(VectorComponentsOut+VectorComponentsOutB+VectorComponentsOutC))
    ic=0
    icB=0
    IF (DynamicSoil) THEN
      IF (thOut.AND.thPos>0) THEN
        icB=icB+1
        VectorOut(ibLoc)%Vec(icB)%cb=>VecC(ibLoc)%Vec(thPos)%cb
      END IF  
      IF (RhoVOut.AND.RhoVPos>0) THEN
        icB=icB+1
        VectorOut(ibLoc)%Vec(icB)%cb=>VecC(ibLoc)%Vec(RhoVPos)%cb
      END IF
    END IF
    IF (ShadOut.AND.Radiation) THEN
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>ShadowCell(ibLoc)%cb 
    END IF
    IF (RadiationOut.AND.Radiation) THEN
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>RaddirCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>RaddifCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>RadinfCell(ibLoc)%cb 
    END IF
    IF (FluxOut.AND.Radiation) THEN
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>SensFluxCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>LatFluxCell(ibLoc)%cb 
    END IF
    IF (LandClassOut.AND.Radiation) THEN
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>LandClassCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>RoughnessLengthCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>AlbedoCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>BulkCoeffDragCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>BulkCoeffHeatCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>BulkCoeffMoistCell(ibLoc)%cb 
    END IF
    IF (thOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(thPos)%c
    END IF
    IF (EnOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(EnPos)%c
    END IF
    IF (tkeOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tkePos)%c
    END IF
    IF (disOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(disPos)%c
    END IF
    IF (omeOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(omePos)%c
    END IF
    IF (tkeHOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tkeHPos)%c
    END IF
    IF (tkeVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tkeVPos)%c
    END IF
    IF (LenOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(LenPos)%c
    END IF
    IF (RhoVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoVPos)%c
    END IF
    IF (RhoCOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoCPos)%c
    END IF
    IF (RhoROut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoRPos)%c
    END IF
    IF (RhoIOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoIPos)%c
    END IF
    IF (RhoSOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoSPos)%c
    END IF
    IF (nvOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nvPos)%c
    END IF
    IF (ncOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(ncPos)%c
    END IF
    IF (nrOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nrPos)%c
    END IF
    IF (niOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(niPos)%c
    END IF
    IF (nsOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nsPos)%c
    END IF
    IF (tracer1Out) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tracer1Pos)%c
    END IF
    IF (tracer2Out) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tracer2Pos)%c
    END IF
    IF (RhoOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoPos)%c
    END IF
    IF (PreOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>PreCell(ibLoc)%c
    END IF
    IF (DampOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DampKoeffCell(ibLoc)%c
    END IF
    IF (DiffOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffKoeff(ibLoc)%c
    END IF
    IF (DiffHOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffHKoeff(ibLoc)%c
    END IF
    IF (DiffVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffVKoeff(ibLoc)%c
    END IF
    IF (DiffPotOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffPotKoeff(ibLoc)%c
    END IF
    IF (DiffPotHOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffPotHKoeff(ibLoc)%c
    END IF
    IF (DiffPotVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffPotVKoeff(ibLoc)%c
    END IF
    IF (DiffMomOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffMomKoeff(ibLoc)%c
    END IF
    IF (DiffMomHOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffMomHKoeff(ibLoc)%c
    END IF
    IF (DiffMomVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffMomVKoeff(ibLoc)%c
    END IF
    IF (ChemieOut) THEN
      DO icT=1,NumAeroOut
        ic=ic+1
        VectorOut(ibLoc)%Vec(ic)%c=>VecT(ibLoc)%Vec(AeroOut(icT))%c 
        OutputC(ic)%iPosC=AeroOut(icT)
      END DO
      DO icT=1,NumGasOut
        ic=ic+1
        VectorOut(ibLoc)%Vec(ic)%c=>VecT(ibLoc)%Vec(GasOut(icT))%c 
        OutputC(ic)%iPosC=GasOut(icT)
      END DO
    END IF
    DO icT=1,LenOutSpecial
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>OutSpecial(ibLoc)%Vec(icT)%c 
    END DO
  END DO

!-- Determination of WriteOffsetC and WriteOffsetN
  Floor(1)%WriteOffsetC=0
  Floor(1)%WriteOffsetN=0
  Floor(1)%WriteOffsetCgmv=0
  Floor(1)%WriteOffsetCSoilgmv=0
  Floor(1)%WriteOffsetCCutgmv=0
  DO ib=1,nb-1
    CALL Set(Floor(ib))
    ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
    ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
    iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
    iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
    iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
    iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
    NumOut=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
    Floor(ib+1)%WriteOffsetC=Floor(ib)%WriteOffsetC+NumOut

    NumOut=(ix1-ix0+1)*(iy1-iy0+1)*(iz1-iz0+1)
    Floor(ib+1)%WriteOffsetN=Floor(ib)%WriteOffsetN+NumOut

    NumOut=Floor(ib)%FreeCells
    Floor(ib+1)%WriteOffsetCgmv=Floor(ib)%WriteOffsetCgmv+NumOut

    NumOut=Floor(ib)%NumBoundCell*Domain%nrsoillayers !nzS
    Floor(ib+1)%WriteOffsetCSoilgmv=Floor(ib)%WriteOffsetCSoilgmv+NumOut

    NumOut=Floor(ib)%NumBoundCell
    Floor(ib+1)%WriteOffsetCCutgmv=Floor(ib)%WriteOffsetCCutgmv+NumOut
  END DO
  CALL Set(Floor(nb))
  ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
  ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
  iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
  iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
  iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
  iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
  NumOut=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
  OffsetShift=Floor(nb)%WriteOffsetC+NumOut
  NumOut=Floor(ib)%FreeCells
  OffsetShiftGMV=Floor(nb)%WriteOffsetCgmv+NumOut
  NumOut=Floor(ib)%NumBoundCell*Domain%nrsoillayers !nzS
  OffsetShiftSoilGMV=Floor(nb)%WriteOffsetCSoilgmv+NumOut
  NumOut=Floor(ib)%NumBoundCell
  OffsetShiftCutGMV=Floor(nb)%WriteOffsetCCutgmv+NumOut
  NumberOfCells=OffsetShift
  NumOut=(ix1-ix0+1)*(iy1-iy0+1)*(iz1-iz0+1)
  NumberOfNodes=Floor(nb)%WriteOffsetN+NumOut

  IF (OutputType(1:3)=='AVS') THEN    
    IF (MyId == 0) THEN
!     OPEN(UNIT=90,FILE=TRIM(OutputFilename)//'Kopf',STATUS='UNKNOWN')
      OPEN(UNIT=90,FILE=TRIM(OutputFilename),STATUS='UNKNOWN')
      WRITE(90,*) 'Z-AXIS'
      WRITE(90,*) iz1Out-iz0Out+1
      WRITE(90,'(8(1X,1PE12.5))') (domain%zP(i),i=iz0Out,iz1Out)
      WRITE(90,*)
      WRITE(90,*) 'Y-AXIS'
      WRITE(90,*) iy1Out-iy0Out+1
      WRITE(90,'(8(1X,1PE12.5))') (domain%yP(i),i=iy0Out,iy1Out)
      WRITE(90,*)
      WRITE(90,*) 'X-AXIS'
      WRITE(90,*) ix1Out-ix0Out+1
      WRITE(90,'(8(1X,1PE12.5))') (domain%xP(i),i=ix0Out,ix1Out)
      WRITE(90,*)
      WRITE(90,*) 'R-AXIS'
      WRITE(90,*) nFrac
      WRITE(90,'(8(1X,1PE12.5))') (1.d6*r(i),i=1,nFrac)
      WRITE(90,*)
      WRITE(90,*) 'TIME'
      WRITE(90,*) NumTime
      WRITE(90,*) 0.0e0
      WRITE(90,*)
      WRITE(90,*) 'GRID'
      WRITE(90,'(8I10)') 0,ix1Out-ix0Out  
      WRITE(90,'(8I10)') 0,iy1Out-iy0Out  
      WRITE(90,'(8I10)') 0,iz1Out-iz0Out  
      nbOut=0
      DO ib=1,nb
        CALL Set(Floor(ib))
        IF (MIN(MAX(igx0,ix0Out),ix1Out)<     &
            MAX(MIN(igx1,ix1Out),ix0Out).AND. &
            MIN(MAX(igy0,iy0Out),iy1Out)<     &
            MAX(MIN(igy1,iy1Out),iy0Out).AND. &
            MIN(MAX(igz0,iz0Out),iz1Out)<     &
            MAX(MIN(igz1,iz1Out),iz0Out)) THEN
          nbOut=nbOut+1
        END IF
      END DO
      WRITE(90,*) nbOut
      DO ib=1,nb
        CALL Set(Floor(ib))
        IF (MIN(MAX(igx0,ix0Out),ix1Out)<     &
            MAX(MIN(igx1,ix1Out),ix0Out).AND. &
            MIN(MAX(igy0,iy0Out),iy1Out)<     &
            MAX(MIN(igy1,iy1Out),iy0Out).AND. &
            MIN(MAX(igz0,iz0Out),iz1Out)<     &
            MAX(MIN(igz1,iz1Out),iz0Out)) THEN
          WRITE(90,'(10I10)') MIN(MAX(igx0,ix0Out),ix1Out)-ix0Out, &
                              MAX(MIN(igx1,ix1Out),ix0Out)-ix0Out, &
                              MIN(MAX(igy0,iy0Out),iy1Out)-iy0Out, &
                              MAX(MIN(igy1,iy1Out),iy0Out)-iy0Out, &
                              MIN(MAX(igz0,iz0Out),iz1Out)-iz0Out, &
                              MAX(MIN(igz1,iz1Out),iz0Out)-iz0Out, &
                              RefineX,RefineY,RefineZ,ib
        END IF
      END DO
      WRITE(90,*)
      WRITE(90,*) 'MET-DATA'
      WRITE(90,*) VectorComponentsOut+1
      IF (VelOut)   WRITE(90,*) 'WIND(3)'
      IF (thOut)    WRITE(90,*) 'th'
      IF (EnOut)    WRITE(90,*) 'En'
      IF (tkeOut)   WRITE(90,*) 'tke'
      IF (disOut)   WRITE(90,*) 'dis'
      IF (omeOut)   WRITE(90,*) 'ome'
      IF (tkeHOut)  WRITE(90,*) 'tkeH'
      IF (tkeVOut)  WRITE(90,*) 'tkeV'
      IF (LenOut)   WRITE(90,*) 'Len'
      IF (RhoVOut)    WRITE(90,*) 'RhoV'
      IF (RhoCOut)    WRITE(90,*) 'RhoC'
      IF (RhoROut)    WRITE(90,*) 'RhoR'
      IF (RhoIOut)    WRITE(90,*) 'RhoI'
      IF (RhoSOut)    WRITE(90,*) 'RhoS'
      IF (nvOut)    WRITE(90,*) 'nv'
      IF (ncOut)    WRITE(90,*) 'nc'
      IF (nrOut)    WRITE(90,*) 'nr'
      IF (niOut)    WRITE(90,*) 'ni'
      IF (nsOut)    WRITE(90,*) 'ns'
      IF (tracer1Out)    WRITE(90,*) 'tracer1'
      IF (tracer2Out)    WRITE(90,*) 'tracer2'
      IF (RhoOut)  WRITE(90,*) 'Rho'
      IF (PreOut)  WRITE(90,*) 'Pre'
      IF (DampOut)  WRITE(90,*) 'Damp'
      IF (DiffOut)  WRITE(90,*) 'D'
      IF (DiffHOut)  WRITE(90,*) 'DH'
      IF (DiffVOut)  WRITE(90,*) 'DV'      
      IF (DiffPotOut)  WRITE(90,*) 'DPot'      
      IF (DiffPotHOut)  WRITE(90,*) 'DPot'      
      IF (DiffPotVOut)  WRITE(90,*) 'DPot'      
      IF (DiffMomOut)  WRITE(90,*) 'DMom'      
      IF (DiffMomHOut)  WRITE(90,*) 'DMomH'
      IF (DiffMomVOut)  WRITE(90,*) 'DMomV'
      IF (ChemieOut) THEN
        DO icT=1,NumAeroOut
          WRITE(nFracChar,'(I2)') nFrac
          WRITE(90,*) TRIM(SpeciesName(AeroOut(icT))) &
                      //'('//nFracChar//')'
        END DO
        DO icT=1,NumGasOut
          WRITE(90,*) SpeciesName(GasOut(icT))
        END DO
      END IF
      WRITE(90,*)
      WRITE(90,*) 'BINARY'
      CLOSE(90)

    END IF
    CALL MPI_Barrier(MPI_Comm_World,MPIerr)
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,OutputFileName, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL,fh,MPIErr)
    CALL MPI_FILE_GET_SIZE(fh,Offset0,MPIErr)
    CALL MPI_FILE_CLOSE(fh,MPIErr)
  ELSE
    IF (OutputType(4:4)=='A'.OR.OutputType(4:4)=='a') THEN
      gmvOutputType=ascii 
    ELSE
      gmvOutputType=iecxi4r4
    END IF
!   CALL gmvOpen(TRIM(OutputFilename)//'.gmvG')
!   CALL gmvGeometry
!   CALL gmvClose
  END IF

!-----------------------------------------------------------
END SUBROUTINE OutputInit

SUBROUTINE OutputSet(VecC,VecT)

  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)
    
!---  local variables
    
  INTEGER :: i,NumOut,nbOut
  INTEGER :: ic,icT
  INTEGER :: icB

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    ic=0
    icB=0
    IF (DynamicSoil) THEN
      IF (thOut.AND.thPos>0) THEN
        icB=icB+1
        VectorOut(ibLoc)%Vec(icB)%cb=>VecC(ibLoc)%Vec(thPos)%cb
      END IF
      IF (RhoVOut.AND.RhoVPos>0) THEN
        icB=icB+1
        VectorOut(ibLoc)%Vec(icB)%cb=>VecC(ibLoc)%Vec(RhoVPos)%cb
      END IF
    END IF
    IF (ShadOut.AND.Radiation) THEN
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>ShadowCell(ibLoc)%cb 
    END IF
    IF (RadiationOut.AND.Radiation) THEN
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>RaddirCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>RaddifCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>RadinfCell(ibLoc)%cb 
    END IF
    IF (FluxOut.AND.Radiation) THEN
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>SensFluxCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>LatFluxCell(ibLoc)%cb 
    END IF
    IF (LandClassOut.AND.Radiation) THEN
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>LandClassCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>RoughnessLengthCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>AlbedoCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>BulkCoeffDragCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>BulkCoeffHeatCell(ibLoc)%cb 
      icB=icB+1
      VectorOut(ibLoc)%Vec(icB)%cb=>BulkCoeffMoistCell(ibLoc)%cb 
    END IF
    IF (thOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(thPos)%c
    END IF
    IF (EnOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(EnPos)%c
    END IF
    IF (tkeOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tkePos)%c
    END IF
    IF (disOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(disPos)%c
    END IF
    IF (omeOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(omePos)%c
    END IF
    IF (tkeHOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tkeHPos)%c
    END IF
    IF (tkeVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tkeVPos)%c
    END IF
    IF (LenOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(LenPos)%c
    END IF    
    IF (RhoVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoVPos)%c
    END IF
    IF (RhoCOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoCPos)%c
    END IF
    IF (RhoROut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoRPos)%c
    END IF
    IF (RhoIOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoIPos)%c
    END IF
    IF (RhoSOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoSPos)%c
    END IF
    IF (nvOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nvPos)%c
    END IF
    IF (ncOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(ncPos)%c
    END IF
    IF (nrOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nrPos)%c
    END IF
    IF (niOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(niPos)%c
    END IF
    IF (nsOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nsPos)%c
    END IF
    IF (tracer1Out) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tracer1Pos)%c
    END IF
    IF (tracer2Out) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tracer2Pos)%c
    END IF
    IF (RhoOut) THEN
      ic=ic+1
!     VectorOut(ibLoc)%Vec(ic)%c=>RhoCell(ibLoc)%c
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoPos)%c
    END IF
    IF (PreOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>PreCell(ibLoc)%c
    END IF
    IF (DampOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DampKoeffCell(ibLoc)%c
    END IF
    IF (DiffOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffKoeff(ibLoc)%c
    END IF
    IF (DiffHOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffHKoeff(ibLoc)%c
    END IF
    IF (DiffVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffVKoeff(ibLoc)%c
    END IF    
    IF (DiffPotOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffPotKoeff(ibLoc)%c
    END IF
    IF (DiffPotHOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffPotHKoeff(ibLoc)%c
    END IF
    IF (DiffPotVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffPotVKoeff(ibLoc)%c
    END IF
    IF (DiffMomOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffMomKoeff(ibLoc)%c
    END IF
    IF (DiffMomHOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffMomHKoeff(ibLoc)%c
    END IF
    IF (DiffMomVOut) THEN
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>DiffMomVKoeff(ibLoc)%c
    END IF
    IF (ChemieOut) THEN
      DO icT=1,NumAeroOut
        ic=ic+1
        VectorOut(ibLoc)%Vec(ic)%c=>VecT(ibLoc)%Vec(AeroOut(icT))%c
      END DO
      DO icT=1,NumGasOut
        ic=ic+1
        VectorOut(ibLoc)%Vec(ic)%c=>VecT(ibLoc)%Vec(GasOut(icT))%c
      END DO
    END IF
    DO icT=1,LenOutSpecial
      ic=ic+1
      VectorOut(ibLoc)%Vec(ic)%c=>OutSpecial(ibLoc)%Vec(icT)%c 
    END DO
  END DO

END SUBROUTINE OutputSet

FUNCTION OutputControl(Time)

  LOGICAL :: OutputControl
  REAL(RealKind) :: Time
  REAL(RealKind) :: u1,u2,u3

  OutputControl=.FALSE.
  IF (Step==OutputStep.OR.Step==OutputEnd) THEN
    IF (Time>=OutputTimeStart) THEN
      OutputControl=.TRUE.
      OutputStep=OutputStep+OutputFrequ
      OutputTime=MAX(Time,OutputTime+OutputTimeStep)   
    END IF
  END IF
  Step=Step+1
  IF (Time>=OutputTime) THEN
    OutputControl=.TRUE.
    OutputTime=MAX(Time,OutputTime+OutputTimeStep)   
  END IF

END FUNCTION OutputControl

SUBROUTINE OutputMet(ActTime,Vel,VecC,VecT)
!
!==================================================
!----  Output of Velocity Values
!==================================================
!

  REAL(RealKind) :: ActTime
  TYPE(VelocityFace_T) :: Vel(:)
  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)


!---  local variables
  INTEGER :: ix,iy,iz,it,i
  INTEGER :: itL,itU
  INTEGER :: iScalar,NumOut
  
  REAL(4), ALLOCATABLE :: Work4(:)
  REAL(4) :: Time,Temp
  INTEGER :: iWork4
  INTEGER :: ncMax
  INTEGER :: nReals
  CHARACTER(10) :: iName, iFrac
  CHARACTER(40) :: DataRep
  CHARACTER(50) :: tempName
 

  IF (Step==OutputStart) THEN
    CALL OutputInit(VecC,VecT)
  ENDIF 
  IF (OutputControl(ActTime)) THEN
    !IF (MyID==0) WRITE(*,*) 'OutputTime',ActTime 
    CALL OutputSet(VecC,VecT)
    WRITE(iName,'(I8)') gmvStep
    gmvStep=gmvStep+1
    CALL gmvOpenWrite(TRIM(OutputDirName)//TRIM(OutputFilename)//TRIM(ADJUSTL(iName)))
    CALL gmvVelocityFWrite(Vel,RhoCell)
    CALL gmvOpenScalarWrite
    DO iScalar=1,VectorComponentsOut
      itL=1
      itU=1
      IF (OutputC(iScalar)%Type=='Aerosol') THEN
        itU=nFrac
      END IF  
      DO it=itL,itU
        IF (itU>1) THEN
          WRITE(iFrac,'(I8)') it
          iFrac='_'//ADJUSTL(iFrac)
        ELSE
          iFrac=''
        END IF
        tempName=TRIM(ADJUSTL(OutputC(iScalar)%NameScalar))//TRIM(iFrac)
        IF (OutputC(iScalar)%DiffRho=='Rho') THEN
          CALL gmvScalarWrite(VectorOut,iScalar,tempName,it,RhoCell)
        ELSE IF (OutputC(iScalar)%Type=='Aerosol' .OR. OutputC(iScalar)%Type=='Aerosol_Spec') THEN
          CALL gmvScalarWrite(VectorOut,iScalar,tempName,it,VecC=VecC,Type=OutputC(iScalar)%DiffRho,&
                              UnitKonv=OutputC(iScalar)%UnitKonvFactor)
        ELSE IF (OutputC(iScalar)%Type=='Gas') THEN
          CALL gmvScalarWrite(VectorOut,iScalar,tempName,it,RhoCell,VecC=VecC,Type=OutputC(iScalar)%DiffRho,&
                              UnitKonv=OutputC(iScalar)%UnitKonvFactor)
        ELSE
          CALL gmvScalarWrite(VectorOut,iScalar,tempName,it)
        END IF
      END DO
    END DO
    CALL gmvCloseScalar
    CALL gmvProbtimeWrite(ActTime)
    CALL gmvClose

    IF(DynamicSoil) THEN
       CALL gmvOpenSoilWrite(TRIM(OutputNameSoil)//suffix_out//TRIM(ADJUSTL(iName)))
       CALL gmvIncludeSoilWrite(TRIM(OutputNameSoil)//TRIM(suffix_outgmvG))
       CALL gmvScalarVarSoilWrite
       DO iScalar=1,VectorComponentsOutB
         CALL gmvSoilVecWrite(VectorOut ,iScalar,NameScalarsB(iScalar),RhoCell)
       END DO
       CALL gmvCloseScalarVarSoil
       CALL gmvProbtimeWriteSoil(ActTime)
       CALL gmvCloseSoil
    END IF ! (DynamicSoil)
    IF (ShadOut.OR.RadiationOut.OR.FluxOut) THEN
       CALL gmvOpenCutWrite(TRIM(OutputNameCut)//suffix_out//TRIM(ADJUSTL(iName)))
       CALL gmvIncludeCutWrite(TRIM(OutputNameCut)//TRIM(suffix_outgmvG))
       CALL gmvScalarVarCutWrite
       DO iScalar=1,VectorComponentsOutB+VectorComponentsOutC
         CALL gmvCutVecWrite(VectorOut ,iScalar,NameScalarsB(iScalar),RhoCell)
       END DO
       CALL gmvCloseScalarVarCut
       CALL gmvProbtimeWriteCut(ActTime)
       CALL gmvCloseCut
    END IF
  END IF   ! IF(ActTime)

!-----------------------------------------------------------
END SUBROUTINE OutputMet


SUBROUTINE InputMet(ActTime,VecC,VecT,InputFileName)
!
!==================================================
!----  Output of Velocity Values
!==================================================
!

  REAL(RealKind) :: ActTime
  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)
  CHARACTER(80) :: InputFileName


!---  local variables
  INTEGER :: ix, iy, iz, it
  INTEGER :: iScalar,NumOut
  
  REAL(4), ALLOCATABLE :: Work4(:)
  REAL(4) :: Time,Temp
  INTEGER :: iWork4,i
  INTEGER :: ncMax
  INTEGER :: nReals
  CHARACTER(10) :: iName
  CHARACTER(40) :: DataRep
  CHARACTER(50) :: SoilFileName
  LOGICAL, SAVE :: Init=.TRUE.

  INTEGER ERRORCLASS,LenString
  CHARACTER(MPI_MAX_ERROR_STRING) :: String



  
  IF (Init) THEN
    Init=.FALSE.
    CALL OutputInit(VecC,VecT)
    Step=Step+1
  END IF 
  CALL OutputSet(VecC,VecT)
  WRITE(iName,'(I8)') gmvStep
  gmvStep=gmvStep+1
  CALL gmvOpenRead(TRIM(InputFilename))
  CALL gmvVelocityFRead(VecC)
  CALL gmvOpenScalarRead
  DO iScalar=1,VectorComponentsOut
    DO it=LBOUND(VectorOut(1)%Vec(iScalar)%c,4),UBOUND(VectorOut(1)%Vec(iScalar)%c,4)
      CALL gmvScalarRead(VectorOut,iScalar,OutputC(iScalar)%NameScalar,it)
    END DO
  END DO

  IF (DynamicSoil) THEN
    SoilFileName=TRIM(InputFileName(1:INDEX(InputFileName,'.')-1))//'.Soil'&
                 //TRIM(InputFileName(INDEX(InputFileName,'.'):LEN(InputFileName)))
    CALL gmvOpenSoilRead(TRIM(SoilFileName))
    CALL gmvIncludeSoilRead(TRIM(SoilFileName))
    CALL gmvScalarVarSoilRead
    DO iScalar=1,VectorComponentsOutB
      IF (TRIM(NameScalarsB(iScalar))=='SoilTemp') THEN
        CALL gmvSoilVecRead(VectorOut,iScalar,NameScalarsB(iScalar),RhoCell)
      ELSE IF (TRIM(NameScalarsB(iScalar))=='SoilWetness') THEN
        CALL gmvSoilVecRead(VectorOut,iScalar,NameScalarsB(iScalar),RhoCell)
      END IF
    END DO
  END IF
  CALL gmvCloseReadScalar
  CALL gmvProbtimeRead(ActTime)
  WRITE(*,*) 'ActTime ',ActTime
  CALL MPI_FILE_CLOSE(fh,MPIErr)
  CALL MPI_ERROR_CLASS(MPIErr, ERRORCLASS, MPIErr)
  CALL MPI_ERROR_STRING(ERRORCLASS,String,LenString,MPIErr)

!-----------------------------------------------------------
END SUBROUTINE InputMet


SUBROUTINE OutputMetC(ActTime,VecC,VecT,FileName,FileStep)
!
!==================================================
!----  Output of Velocity Values
!==================================================
!

  REAL(RealKind) :: ActTime
  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)
  CHARACTER(*), OPTIONAL :: FileName
  INTEGER, OPTIONAL :: FileStep

!---  local variables
  INTEGER :: ix, iy, iz, it
  INTEGER :: iScalar,NumOut
  
  REAL(4), ALLOCATABLE :: Work4(:)
  REAL(4) :: Time,Temp
  INTEGER :: iWork4
  INTEGER :: ncMax
  INTEGER :: nReals
  CHARACTER(10) :: iName
  CHARACTER(40) :: DataRep

  
  IF (Step==OutputStart) THEN
    CALL OutputInit(VecC,VecT)
  ENDIF 
  IF (OutputControl(ActTime)) THEN
    CALL OutputSet(VecC,VecT)
    IF (PRESENT(FileStep)) THEN
      WRITE(iName,'(I8)') FileStep
    ELSE
      WRITE(iName,'(I8)') gmvStep
      gmvStep=gmvStep+1
    END IF
    IF (PRESENT(FileName)) THEN
      CALL gmvOpenWrite(TRIM(OutputDirName)//TRIM(Filename)//TRIM(ADJUSTL(iName)))
    ELSE
      CALL gmvOpenWrite(TRIM(OutputDirName)//TRIM(OutputFilename)//TRIM(ADJUSTL(iName)))
    END IF
    CALL gmvVelocityCWrite(VecC)
    CALL gmvOpenScalarWrite
    DO iScalar=1,SIZE(VectorOut(1)%Vec)
      DO it=LBOUND(VectorOut(1)%Vec(iScalar)%c,4),UBOUND(VectorOut(1)%Vec(iScalar)%c,4)
        CALL gmvScalarWrite(VectorOut,iScalar,OutputC(iScalar)%NameScalar,it)
      END DO
    END DO
    CALL gmvCloseScalar
    CALL gmvProbtimeWrite(ActTime)
    CALL gmvClose
  END IF

!-----------------------------------------------------------
END SUBROUTINE OutputMetC

SUBROUTINE OutputVector(ActTime,Vector,RhoCell)
!
!==================================================
!----  Output of Velocity Values
!==================================================
!

    IMPLICIT NONE

    REAL(RealKind) :: ActTime
    TYPE(Vector4Cell_T) :: Vector(:)
    TYPE(ScalarCell_T) :: RhoCell(:)


!---  local variables
    INTEGER :: ix, iy, iz
    INTEGER :: iScalar,NumOut
    
    REAL(4), ALLOCATABLE :: Work4(:)
    REAL(4) :: Time
    INTEGER :: iWork4
    INTEGER :: ncMax
    CHARACTER(80) :: FileName
    INTEGER :: nReals


    FileName='VelocGlob'
!-- Preparing Output
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,FileName, &
                       MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                       MPI_INFO_NULL,fh,MPIErr)
!-- Ausgabe Zeitpunkt
    Offset=Offset0
    Offset0=Offset0+1
    IF (MyId==0) THEN
      nReals=1
      Time=ActTime
      CALL MPI_FILE_WRITE_AT(fh,Offset,Time,nReals, &
                          MPI_REAL,Status,MPIErr)
    END IF
    ncMax=0
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      ncMax=MAX(ncMax,Floor(ib)%nc)
    END DO
!-- Ausgabe Geschwindigkeitsfeld
    ALLOCATE(Work4(NumData*ncMax))
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
      iWork4=1
      DO iScalar=1,SIZE(Vector(ibLoc)%Vec)
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            DO iz=iz0+1,iz1
              Work4(iWork4)=Vector(ibLoc)%Vec(iScalar)%c(ix,iy,iz &
                            ,LBOUND(Vector(ibLoc)%Vec(iScalar)%c,4)) &
                            /(RhoCell(ibLoc)%c(ix,iy,iz,1)+Eps)
              iWork4=iWork4+1
            END DO
          END DO
        END DO
      END DO
      Offset=NumData*Floor(ib)%WriteOffsetC*4+Offset0
      NumOut=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
      nReals=NumData*NumOut
      CALL MPI_FILE_WRITE_AT(fh,Offset,work4,nReals, &
                            MPI_REAL,Status,MPIErr)
    END DO
    CALL MPI_FILE_CLOSE(fh,MPIErr)
    DEALLOCATE(Work4)
    Offset0=Offset0+4*NumData*OffsetShift

!-----------------------------------------------------------
END SUBROUTINE OutputVector


SUBROUTINE OutputProfile(VecC,VecT,ProfName)
!
!==================================================
!----  Output of Velocity Values
!==================================================
!

  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)

  REAL(RealKind) :: GnuplotOut(0:15)
  INTEGER :: ix,iy,iz
  INTEGER :: ixProf,iyProf
  INTEGER :: i,ibProf,ibProfLoc
  LOGICAL :: OutProf=.FALSE.
  LOGICAL :: OutProfGnu=.FALSE.
  CHARACTER(80) :: ProfName

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        IF (xP(ix-1)<=xProf.AND. &
            xP(ix  )>=xProf.AND. &
            yP(iy-1)<=yProf.AND. &
            yP(iy  )>=yProf) THEN
          ixProf=ix
          iyProf=iy
          ibProf=ib
          ibProfLoc=ibLoc
          OutProf=.TRUE.
          OutProfGnu=.TRUE.
        END IF
      END DO
    END DO
  END DO
  IF (OutProf) THEN
    CALL Set(Floor(ibProf))
    OPEN(UNIT=OutputUnit,FILE=TRIM(ProfName),STATUS='UNKNOWN')
    WRITE(OutputUnit,*) 'uProf'
    WRITE(OutputUnit,*) nz,2
    DO iz=iz0+1,iz1
      WRITE(OutputUnit,*)  &
        0.5d0*(zP(iz-1)+zP(iz)) &
       ,VecC(ibProfLoc)%Vec(uPosL)%c(ixProf,iyProf,iz,1) &
                                /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
    END DO
    WRITE(OutputUnit,*) 'vProf'
    WRITE(OutputUnit,*) nz,2
    DO iz=iz0+1,iz1
      WRITE(OutputUnit,*)  &
        0.5d0*(zP(iz-1)+zP(iz)) &
       ,VecC(ibProfLoc)%Vec(vPosL)%c(ixProf,iyProf,iz,1) &
                                /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
    END DO
    WRITE(OutputUnit,*) 'wProf'
    WRITE(OutputUnit,*) nz,2
    DO iz=iz0+1,iz1
      WRITE(OutputUnit,*)  &
        0.5d0*(zP(iz-1)+zP(iz)) &
       ,0.5d0*(VecC(ibProfLoc)%Vec(wPosL)%c(ixProf,iyProf,iz,1) &
              +VecC(ibProfLoc)%Vec(wPosR)%c(ixProf,iyProf,iz,1)) &
                                /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
    END DO
    IF (thPos>0) THEN
      WRITE(OutputUnit,*) 'ThProf'
      WRITE(OutputUnit,*) nz,2
      DO iz=iz0+1,iz1
        WRITE(OutputUnit,*)  &
          0.5d0*(zP(iz-1)+zP(iz)) &
         ,VecC(ibProfLoc)%Vec(thPos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      END DO
    END IF
    IF (tkePos>0) THEN
      WRITE(OutputUnit,*) 'TkeProf'
      WRITE(OutputUnit,*) nz,2
      DO iz=iz0+1,iz1
        WRITE(OutputUnit,*)  &
          0.5d0*(zP(iz-1)+zP(iz)) &
         ,VecC(ibProfLoc)%Vec(tkePos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      END DO
    END IF
    IF (tkeHPos>0) THEN
      WRITE(OutputUnit,*) 'TkeHProf'
      WRITE(OutputUnit,*) nz,2
      DO iz=iz0+1,iz1
        WRITE(OutputUnit,*)  &
          0.5d0*(zP(iz-1)+zP(iz)) &
         ,VecC(ibProfLoc)%Vec(tkeHPos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      END DO
    END IF
    IF (tkeVPos>0) THEN
      WRITE(OutputUnit,*) 'TkeVProf'
      WRITE(OutputUnit,*) nz,2
      DO iz=iz0+1,iz1
        WRITE(OutputUnit,*)  &
          0.5d0*(zP(iz-1)+zP(iz)) &
         ,VecC(ibProfLoc)%Vec(tkeVPos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      END DO
    END IF
    IF (disPos>0) THEN
      WRITE(OutputUnit,*) 'DisProf'
      WRITE(OutputUnit,*) nz,2
      DO iz=iz0+1,iz1
        WRITE(OutputUnit,*)  &
          0.5d0*(zP(iz-1)+zP(iz)) &
         ,VecC(ibProfLoc)%Vec(disPos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      END DO
    END IF
    IF (omePos>0) THEN
      WRITE(OutputUnit,*) 'OmeProf'
      WRITE(OutputUnit,*) nz,2
      DO iz=iz0+1,iz1
        WRITE(OutputUnit,*)  &
          0.5d0*(zP(iz-1)+zP(iz)) &
         ,VecC(ibProfLoc)%Vec(omePos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      END DO
    END IF

    IF (Diffusion) THEN
      IF(TkeHVLen) THEN                      !!!FILAUS
        WRITE(OutputUnit,*) 'DHProf'
        WRITE(OutputUnit,*) nz,2
        DO iz=iz0+1,iz1
          WRITE(OutputUnit,*)  &
            0.5d0*(zP(iz-1)+zP(iz)) &
           ,DiffHKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        END DO
        WRITE(OutputUnit,*) 'DVProf'
        WRITE(OutputUnit,*) nz,2
        DO iz=iz0+1,iz1
          WRITE(OutputUnit,*)  &
            0.5d0*(zP(iz-1)+zP(iz)) &
           ,DiffVKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        END DO
      ELSE IF(TKESGS.OR.TKELEN) THEN
        WRITE(OutputUnit,*) 'DMomProf'
        WRITE(OutputUnit,*) nz,2
        DO iz=iz0+1,iz1
          WRITE(OutputUnit,*)  &
            0.5d0*(zP(iz-1)+zP(iz)) &
           ,DiffMomKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        END DO
        WRITE(OutputUnit,*) 'DPotProf'
        WRITE(OutputUnit,*) nz,2
        DO iz=iz0+1,iz1
          WRITE(OutputUnit,*)  &
            0.5d0*(zP(iz-1)+zP(iz)) &
           ,DiffPotKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        END DO
      ELSE IF(TKESMAG) THEN
        WRITE(OutputUnit,*) 'DMomHProf'
        WRITE(OutputUnit,*) nz,2
        DO iz=iz0+1,iz1
          WRITE(OutputUnit,*)  &
            0.5d0*(zP(iz-1)+zP(iz)) &
           ,DiffMomHKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        END DO
        WRITE(OutputUnit,*) 'DPotHProf'
        WRITE(OutputUnit,*) nz,2
        DO iz=iz0+1,iz1
          WRITE(OutputUnit,*)  &
            0.5d0*(zP(iz-1)+zP(iz)) &
           ,DiffPotHKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        END DO
      ELSE ! For example, TKEDIS and TKEOME
        WRITE(OutputUnit,*) 'DProf'
        WRITE(OutputUnit,*) nz,2
        DO iz=iz0+1,iz1
          WRITE(OutputUnit,*)  &
            0.5d0*(zP(iz-1)+zP(iz)) &
           ,DiffKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        END DO
      END IF
    END IF ! IF (Diffusion)
    WRITE(OutputUnit,*) 'RhoProf'
    WRITE(OutputUnit,*) nz,2
    DO iz=iz0+1,iz1
      WRITE(OutputUnit,*)  &
        0.5d0*(zP(iz-1)+zP(iz)) &
       ,RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1) 
    END DO
    CLOSE(OutputUnit)
  END IF
  IF (OutProfGnu) THEN                !!!FILAUS
    CALL Set(Floor(ibProf))
    OPEN(UNIT=OutputUnit,FILE='Gnu_'//TRIM(ProfName),STATUS='UNKNOWN')
    WRITE(OutputUnit,*) 'zHeight,uProf,vProf,ThProf,' &
                       ,'TkeProf,TkeHProf,TkeVProf,DisProf,OmeProf,DHProf,DVProf,' &
                       ,'DMomProf,DPotProf,DProf,RhoProf,LenProf'
    DO iz=iz0+1,iz1
      GnuplotOut(0)=0.5d0*(zP(iz-1)+zP(iz))
      GnuplotOut(1)=VecC(ibProfLoc)%Vec(uPosL)%c(ixProf,iyProf,iz,1) &
                                /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      GnuplotOut(2)=VecC(ibProfLoc)%Vec(vPosL)%c(ixProf,iyProf,iz,1) &
                                /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      IF (thPos>0) THEN
        GnuplotOut(3)=VecC(ibProfLoc)%Vec(thPos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      ELSE
        GnuplotOut(3)=0.0d0
      END IF
     
      IF (tkePos>0) THEN 
        GnuplotOut(4)=VecC(ibProfLoc)%Vec(tkePos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps) 
        GnuplotOut(5)=0.0d0
        GnuplotOut(6)=0.0d0
      ELSE IF (tkeHPos>0 .AND. tkeVPos>0) THEN
        GnuplotOut(4)=0.0d0
        GnuplotOut(5)=VecC(ibProfLoc)%Vec(tkeHPos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        GnuplotOut(6)=VecC(ibProfLoc)%Vec(tkeVPos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      ELSE
        GnuplotOut(4)=0.0d0
        GnuplotOut(5)=0.0d0
        GnuplotOut(6)=0.0d0
      END IF 

      IF (disPos>0) THEN
        GnuplotOut(7)=VecC(ibProfLoc)%Vec(disPos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      ELSE
        GnuplotOut(7)=0.0d0
      END IF

      IF (omePos>0) THEN
        GnuplotOut(15)=VecC(ibProfLoc)%Vec(omePos)%c(ixProf,iyProf,iz,1) &
                                  /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
      ELSE
        GnuplotOut(15)=0.0d0
      END IF
      
      IF (Diffusion) THEN
        IF(TkeHVLen) THEN
           GnuplotOut(8)=DiffHKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
           GnuplotOut(9)=DiffVKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
           GnuplotOut(10)=0.0d0
           GnuplotOut(11)=0.0d0
           GnuplotOut(12)=0.0d0
        ELSE IF (TKESGS.OR.TKELEN) THEN
           GnuplotOut(8)=0.0d0
           GnuplotOut(9)=0.0d0
           GnuplotOut(10)=DiffMomKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
           GnuplotOut(11)=DiffPotKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
           GnuplotOut(12)=0.0d0
        ELSE IF (TKESMAG) THEN
           GnuplotOut(8)=0.0d0
           GnuplotOut(9)=0.0d0
           GnuplotOut(10)=DiffMomHKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
           GnuplotOut(11)=DiffPotHKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
           GnuplotOut(12)=0.0d0
        ELSE IF (TKEDIS.OR.TKEDISRICH.OR.TKEOME) THEN
           GnuplotOut(8)=0.0d0
           GnuplotOut(9)=0.0d0
           GnuplotOut(10)=0.0d0
           GnuplotOut(11)=0.0d0
           GnuplotOut(12)=DiffKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1) &
                                    /(RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)+Eps)
        ELSE
           GnuplotOut(8)=0.0d0
           GnuplotOut(9)=0.0d0
           GnuplotOut(10)=0.0d0
           GnuplotOut(11)=0.0d0
           GnuplotOut(12)=0.0d0 
        END IF
      END IF
      GnuplotOut(13)=RhoCell(ibProfLoc)%c(ixProf,iyProf,iz,1)
      IF(TkeHVLen) THEN
         GnuplotOut(14)=VecC(ibProfLoc)%Vec(LenPos)%c(ixProf,iyProf,iz,1)
      ELSE IF(TKESGS) THEN
         GnuplotOut(14)=LenKoeff(ibProfLoc)%c(ixProf,iyProf,iz,1)
      ELSE IF(TkeLen) THEN
         GnuplotOut(14)=VecC(ibProfLoc)%Vec(LenPos)%c(ixProf,iyProf,iz,1)
      ELSE
         GnuplotOut(14)=0.0d0
      END IF
      WRITE(OutputUnit,'(16(D25.16, 1x))') (GnuplotOut(i),i=0,15)
    END DO
    CLOSE(OutputUnit)
  END IF

END SUBROUTINE OutputProfile

SUBROUTINE ReadOutputProfile(VecC,VecT,ProfName)
!
!==================================================
!----  Output of Velocity Values
!==================================================
!

  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)
  TYPE (ScalarCell_T), POINTER :: NoScalarCell(:) ! replacement character

  REAL(RealKind) :: GnuplotOut(0:40)
  INTEGER :: ix,iy,iz
  INTEGER :: ic,icT,icB
  INTEGER :: ixProf,iyProf
  INTEGER :: i,iS,ibProf,k
  INTEGER, PARAMETER :: noPos=9999
  CHARACTER(80)  :: ProfName
  CHARACTER(20)  :: ProfCase
  CHARACTER(400) :: ProfComp
  CHARACTER(400) :: SoilComp
  CHARACTER(80)  :: post
  REAL(RealKind), ALLOCATABLE :: GnuProfOut(:,:),GnuProfOut2(:,:),GnuProfOut3(:,:)
  REAL(RealKind), ALLOCATABLE :: SoilDataOut(:,:)

  INTEGER :: iScalar
  CHARACTER(10) :: iName

  CALL CheckValue
  
  IF (ProfAverage) THEN
    WRITE(*,*) 'xStart',xStart,'yStart',yStart
    WRITE(*,*) 'xEnd',xEnd,'yEnd',yEnd
  ELSE
    WRITE(*,*) 'xProf',xProf,'yProf',yProf
  END IF

  CALL getarg(2,post)
  IF (post=='Met') THEN
    ic=0
    icB=0
    ProfCase='Prof'
    ProfComp='zHeight'
    SoilComp='zHeight'
    ALLOCATE(GnuProfOut(2*VectorComponentsOut,Domain%iz1))
    IF (ProfIslandOcean) THEN
      ALLOCATE(GnuProfOut2(2*VectorComponentsOut,Domain%iz1))
      ALLOCATE(GnuProfOut3(2*VectorComponentsOut,Domain%iz1))
      OPEN(UNIT=OutputUnit+1,FILE=TRIM(ProfName)//'Island',STATUS='UNKNOWN')
      OPEN(UNIT=OutputUnit+2,FILE=TRIM(ProfName)//'Ocean',STATUS='UNKNOWN')
    END IF
    ALLOCATE(SoilDataOut(VectorComponentsOutB,Domain%nrsoillayers))
    OPEN(UNIT=OutputUnit,FILE=TRIM(ProfName),STATUS='UNKNOWN')
    IF (uPosL>0) THEN
      WRITE(OutputUnit,*) 'uProf'
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', uProf'
      CALL MakeProfile(NoScalarCell,VecC,uPosL,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,uPosL,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,uPosL,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (vPosL>0) THEN
      WRITE(OutputUnit,*) 'vProf'
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', vProf'
      CALL MakeProfile(NoScalarCell,VecC,vPosL,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,vPosL,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,vPosL,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (wPosL>0) THEN
      WRITE(OutputUnit,*) 'wProf'
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', wProf'
      CALL MakeProfile(NoScalarCell,VecC,wPosL,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,wPosL,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,wPosL,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (thPos>0) THEN
      ic=ic+1
      icB=icB+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', ThProf'
      SoilComp=SoilComp(1:LEN(TRIM(SoilComp)))//', SoilTemp'
      WRITE(OutputUnit,*) 'ThProf'
      CALL MakeProfile(NoScalarCell,VecC,thPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,thPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,thPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
      IF (DynamicSoil) THEN
        CALL MakeSoil(VecC,thPos,icB,SoilDataOut)
      END IF
    END IF
    IF (tkePos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', TkeProf'
      WRITE(OutputUnit,*) 'TkeProf'
      CALL MakeProfile(NoScalarCell,VecC,tkePos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,tkePos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,tkePos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (tkeHPos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', TkeHProf'
      WRITE(OutputUnit,*) 'TkeHProf'
      CALL MakeProfile(NoScalarCell,VecC,tkeHPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,tkeHPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,tkeHPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (tkeVPos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', TkeVProf'
      WRITE(OutputUnit,*) 'TkeVProf'
      CALL MakeProfile(NoScalarCell,VecC,tkeVPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,tkeVPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,tkeVPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (disPos>0) THEN
      ic=ic+1
      WRITE(OutputUnit,*) 'DisProf'
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DisProf'
      CALL MakeProfile(NoScalarCell,VecC,disPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,disPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,disPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (omePos>0) THEN
      ic=ic+1
      WRITE(OutputUnit,*) 'OmeProf'
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', OmeProf'
      CALL MakeProfile(NoScalarCell,VecC,omePos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,omePos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,omePos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (Diffusion) THEN
      IF(TkeHVLen) THEN                      !!!FILAUS
        IF (DiffHOut) THEN
          WRITE(OutputUnit,*) 'DHProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DHProf'
          CALL MakeProfile(DiffHKoeff,VecC,noPos,1,ic,GnuProfOut)
          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
          IF (ProfIslandOcean) THEN
            CALL MakeProfile(DiffHKoeff,VecC,noPos,1,ic,GnuProfOut2,Island=.TRUE.)
            CALL MakeProfile(DiffHKoeff,VecC,noPos,1,ic,GnuProfOut3,Island=.FALSE.)
            CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
            CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
          END IF
        END IF
        IF (DiffVOut) THEN
          WRITE(OutputUnit,*) 'DVProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DVProf'
          CALL MakeProfile(DiffVKoeff,VecC,noPos,1,ic,GnuProfOut)
          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
          IF (ProfIslandOcean) THEN
            CALL MakeProfile(DiffVKoeff,VecC,noPos,1,ic,GnuProfOut2,Island=.TRUE.)
            CALL MakeProfile(DiffVKoeff,VecC,noPos,1,ic,GnuProfOut3,Island=.FALSE.)
            CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
            CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
          END IF
        END IF
      ELSE IF(TKESGS.OR.TKELEN) THEN
        IF (DiffMomOut) THEN
          WRITE(OutputUnit,*) 'DMomProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DMomProf'
          CALL MakeProfile(DiffMomKoeff,VecC,noPos,1,ic,GnuProfOut)
          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
          IF (ProfIslandOcean) THEN
            CALL MakeProfile(DiffMomKoeff,VecC,noPos,1,ic,GnuProfOut2,Island=.TRUE.)
            CALL MakeProfile(DiffMomKoeff,VecC,noPos,1,ic,GnuProfOut3,Island=.FALSE.)
            CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
            CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
          END IF
        END IF
        IF (DiffPotOut) THEN
          WRITE(OutputUnit,*) 'DPotProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DPotProf'
          CALL MakeProfile(DiffPotKoeff,VecC,noPos,1,ic,GnuProfOut)
          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
          IF (ProfIslandOcean) THEN
            CALL MakeProfile(DiffPotKoeff,VecC,noPos,1,ic,GnuProfOut2,Island=.TRUE.)
            CALL MakeProfile(DiffPotKoeff,VecC,noPos,1,ic,GnuProfOut3,Island=.FALSE.)
            CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
            CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
          END IF
        END IF
      ELSE IF(TKESMAG.or.DynSmag) THEN
        IF (DiffMomHOut) THEN
          WRITE(OutputUnit,*) 'DMomHProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DMomHProf'
          CALL MakeProfile(DiffMomHKoeff,VecC,noPos,1,ic,GnuProfOut)
!          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
        END IF
        IF (DiffMomVOut) THEN
          WRITE(OutputUnit,*) 'DMomVProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DMomVProf'
          CALL MakeProfile(DiffMomVKoeff,VecC,noPos,1,ic,GnuProfOut)
!          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
        END IF
        IF (DiffPotHOut) THEN
          WRITE(OutputUnit,*) 'DPotHProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DPotHProf'
          CALL MakeProfile(DiffPotHKoeff,VecC,noPos,1,ic,GnuProfOut)
!          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
        END IF
        IF (DiffPotVOut) THEN
          WRITE(OutputUnit,*) 'DPotVProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DPotVProf'
          CALL MakeProfile(DiffPotVKoeff,VecC,noPos,1,ic,GnuProfOut)
!          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
        END IF
      ELSE 
        IF (DiffOut) THEN
          WRITE(OutputUnit,*) 'DProf'
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', DProf'
          CALL MakeProfile(DiffKoeff,VecC,noPos,1,ic,GnuProfOut)
          CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
          IF (ProfIslandOcean) THEN
            CALL MakeProfile(DiffKoeff,VecC,noPos,1,ic,GnuProfOut2,Island=.TRUE.)
            CALL MakeProfile(DiffKoeff,VecC,noPos,1,ic,GnuProfOut3,Island=.FALSE.)
            CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
            CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
          END IF
        END IF
      END IF
    END IF
    IF (RhoPos>0) THEN 
      WRITE(OutputUnit,*) 'RhoProf'
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', RhoProf'
      CALL MakeProfile(NoScalarCell,VecC,RhoPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,RhoPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,RhoPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (PreOut) THEN 
      WRITE(OutputUnit,*) 'PreProf'
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', PreProf'
      CALL MakeProfile(PreCell,VecC,noPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(PreCell,VecC,noPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(PreCell,VecC,noPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (RhoVPos>0) THEN
      ic=ic+1
      icB=icB+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', RhoVProf'
      SoilComp=SoilComp(1:LEN(TRIM(SoilComp)))//', SoilWetness'
      WRITE(OutputUnit,*) 'RhoVProf'
      CALL MakeProfile(NoScalarCell,VecC,RhoVPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,RhoVPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,RhoVPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
      IF (DynamicSoil) THEN
        CALL MakeSoil(VecC,RhoVPos,icB,SoilDataOut)
      END IF
    END IF
    IF (RhoCPos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', RhoCProf'
      WRITE(OutputUnit,*) 'RhoCProf'
      CALL MakeProfile(NoScalarCell,VecC,RhoCPos,1,ic,GnuProfOut,IsMet=.TRUE.) ! compute cloud cover
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp,IsMet=.TRUE.)        ! output cloud cover
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,RhoCPos,1,ic,GnuProfOut2,IsMet=.TRUE.,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,RhoCPos,1,ic,GnuProfOut3,IsMet=.TRUE.,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,IsMet=.TRUE.,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,IsMet=.TRUE.,Island=.FALSE.)
      END IF
    END IF
    IF (RhoRPos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', RhoRProf'
      WRITE(OutputUnit,*) 'RhoRProf'
      CALL MakeProfile(NoScalarCell,VecC,RhoRPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,RhoRPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,RhoRPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (nvPos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', NCCNProf'
      WRITE(OutputUnit,*) 'NCCNProf'
      CALL MakeProfile(NoScalarCell,VecC,nvPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,nvPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,nvPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (ncPos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', NCProf'
      WRITE(OutputUnit,*) 'NCProf'
      CALL MakeProfile(NoScalarCell,VecC,ncPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,ncPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,ncPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (nrPos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', NRProf'
      WRITE(OutputUnit,*) 'NRProf'
      CALL MakeProfile(NoScalarCell,VecC,nrPos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,nrPos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,nrPos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (tracer1Pos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', Tracer1Prof'
      WRITE(OutputUnit,*) 'Tracer1Prof'
      CALL MakeProfile(NoScalarCell,VecC,tracer1Pos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,tracer1Pos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,tracer1Pos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF
    IF (tracer2Pos>0) THEN
      ic=ic+1
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', Tracer2Prof'
      WRITE(OutputUnit,*) 'Tracer2Prof'
      CALL MakeProfile(NoScalarCell,VecC,tracer2Pos,1,ic,GnuProfOut)
      CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
      IF (ProfIslandOcean) THEN
        CALL MakeProfile(NoScalarCell,VecC,tracer2Pos,1,ic,GnuProfOut2,Island=.TRUE.)
        CALL MakeProfile(NoScalarCell,VecC,tracer2Pos,1,ic,GnuProfOut3,Island=.FALSE.)
        CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
        CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
      END IF
    END IF

    CLOSE(OutputUnit)
    IF (ProfIslandOcean) THEN
      CLOSE(OutputUnit+1)
      CLOSE(OutputUnit+2)
    END IF

    OPEN(UNIT=OutputUnit,FILE='Gnu_'//TRIM(ProfName),STATUS='UNKNOWN')
    IF (ProfIslandOcean) THEN
      OPEN(UNIT=OutputUnit+1,FILE='Gnu_'//TRIM(ProfName)//'Island',STATUS='UNKNOWN')
      OPEN(UNIT=OutputUnit+2,FILE='Gnu_'//TRIM(ProfName)//'Ocean',STATUS='UNKNOWN')
    END IF
    ProfCase='GnuProf'
    IF (Deviation) THEN
      ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', deviation: '//ProfComp(9:LEN(TRIM(ProfComp)))
      ic=ic+ic
    END IF
    CALL WriteProfile(GnuProfOut,ic,ProfCase,ProfComp)
    IF (ProfIslandOcean) THEN
      CALL WriteProfile(GnuProfOut2,ic,ProfCase,ProfComp,Island=.TRUE.)
      CALL WriteProfile(GnuProfOut3,ic,ProfCase,ProfComp,Island=.FALSE.)
    END IF
    CLOSE(OutputUnit)

    OPEN(UNIT=OutputUnit,FILE='Soil'//TRIM(ProfName),STATUS='UNKNOWN')
    CALL WriteSoil(SoilDataOut,icB,SoilComp)
    CLOSE(OutputUnit)
    IF (ProfIslandOcean) THEN
      CLOSE(OutputUnit+1)
      CLOSE(OutputUnit+2)
    END IF

    DEALLOCATE(GnuProfOut)
    DEALLOCATE(SoilDataOut)

  ELSE

    IF (post=='Special') THEN
      ALLOCATE(GnuProfOut(LenOutSpecial,iz1))
      ProfComp='zHeight   '
      DO i=1,LenOutSpecial
        ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//', '//NameOutSpecial(i)
      END DO
      ProfCase='Prof'
      OPEN(UNIT=OutputUnit,FILE=TRIM(ProfName),STATUS='UNKNOWN')
      DO i=1,LenOutSpecial
        CALL MakeProfile(NoScalarCell,OutSpecial,i,1,i,GnuProfOut)
        CALL WriteProfile(GnuProfOut,i,ProfCase,ProfComp)
      END DO
      CLOSE(OutputUnit)
      OPEN(UNIT=OutputUnit,FILE='Gnu_'//TRIM(ProfName),STATUS='UNKNOWN')
      ProfCase='GnuProf'
      CALL WriteProfile(GnuProfOut,LenOutSpecial,ProfCase,ProfComp)
      CLOSE(OutputUnit)
      DEALLOCATE(GnuProfOut)

    ELSE IF (post=='Besonderes') THEN
      WRITE(*,*) 'Noch nichts zu tun'

    ELSE IF (post=='aNUMBER') THEN
      ALLOCATE(GnuProfOut(nFrac,iz1))
      ProfCase='GnuProf'
      ProfComp='zHeight   '//'aNUMBER'
      OPEN(UNIT=OutputUnit,FILE='Gnu_'//TRIM(ProfName),STATUS='UNKNOWN')
      DO i=1,nFrac
        CALL MakeProfile(NoScalarCell,VecC,Position('aNUMBER'),i,i,GnuProfOut)
      END DO
      CALL WriteProfile(GnuProfOut,nFrac,ProfCase,ProfComp)
      DEALLOCATE(GnuProfOut)
      CLOSE(OutputUnit)
    ELSE
      OPEN(UNIT=OutputUnit,FILE='Gnu_'//TRIM(ProfName),STATUS='UNKNOWN')
      ALLOCATE(GnuProfOut(nFrac,iz1))
      ProfCase='GnuProf'
      ProfComp='zHeight   '
      ic=0
      DO k=1,80
        IF (post(k:k).EQ.'p'.OR.post(k:k)=='m'.OR.post(k:k)=='a'.OR.post(k:k)=='s') THEN ! Aerosolinhaltstoff
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//TRIM(post) 
          DO i=1,nFrac
            CALL MakeProfile(NoScalarCell,VecC,Position(trim(post)),i,i,GnuProfOut)
            DO iz=iz0+1,iz1
              GnuProfOut(i,iz)=GnuProfOut(i,iz)/Molmass(Position(trim(post)))  ! mol/m^3(Luft)
            END DO
          END DO
          CALL WriteProfile(GnuProfOut,nFrac,ProfCase,ProfComp)
          EXIT
        END IF
        IF (post(k:k).EQ.' ') THEN ! Gas
          ic=ic+1
          ProfComp=ProfComp(1:LEN(TRIM(ProfComp)))//TRIM(post)
          CALL MakeProfile(NoScalarCell,VecC,Position(trim(post)),1,ic,GnuProfOut) !kg/m^3
          CALL WriteProfile(GnuProfOut,nFrac,ProfCase,ProfComp)
          EXIT
        END IF
      END DO
      DEALLOCATE(GnuProfOut)
      CLOSE(OutputUnit)
    END IF

  END IF

END SUBROUTINE ReadOutputProfile

SUBROUTINE OutputPointProfile(VecC)

  TYPE(Vector4Cell_T) :: VecC(:)

  TYPE(ScalarCell_T) :: ScalarC(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  CHARACTER*20  :: ProfileName
  INTEGER :: Pos
  REAL(RealKind) :: Profile(Domain%nz,SIZE(VectorOut)+3)

WRITE(*,*) 'Im in OutputPointprof', myID

  IF (RhoPos>0) THEN
    CALL Assign(Rho,VecC,RhoPos)
  END IF  
  IF (uPosL>0) THEN

  WRITE(*,*) 'Im in OutputPointprof after U', myID

    CALL Assign(ScalarC,VecC,uPosl)
    IF (RhoPos>0) THEN
      CALL PointProfile(ScalarC,Profile(:,SIZE(VectorOut)+1),Rho)
    ELSE  
      CALL PointProfile(ScalarC,Profile(:,SIZE(VectorOut)+1))
    END IF  

    WRITE(*,*) 'Im in OutputPointprof after UU', myID

    ProfileName='UPoint'
    CALL WriteProfile1(Profile(:,SIZE(VectorOut)+1),Domain%zP,ProfileName)
  END IF  
  IF (vPosL>0) THEN
    CALL Assign(ScalarC,VecC,vPosl)
    IF (RhoPos>0) THEN
      CALL PointProfile(ScalarC,Profile(:,SIZE(VectorOut)+2),Rho)
    ELSE  
      CALL PointProfile(ScalarC,Profile(:,SIZE(VectorOut)+2))
    END IF  
    ProfileName='VPoint'
    CALL WriteProfile1(Profile(:,SIZE(VectorOut)+2),Domain%zP,ProfileName)

    WRITE(*,*) 'Im in OutputPointprof V', myID
  END IF  
  IF (wPosL>0) THEN
    CALL Assign(ScalarC,VecC,wPosl)
    IF (RhoPos>0) THEN
      CALL PointProfile(ScalarC,Profile(:,SIZE(VectorOut)+3),Rho)
    ELSE  
      CALL PointProfile(ScalarC,Profile(:,SIZE(VectorOut)+3))
    END IF  
    ProfileName='WPoint'
    CALL WriteProfile1(Profile(:,SIZE(VectorOut)+3),Domain%zP,ProfileName)
  END IF  
  DO Pos=1,SIZE(OutputC)
    CALL Assign(ScalarC,VectorOut,Pos)
    IF (RhoPos>0.AND.OutputC(Pos)%DiffRho=='Rho') THEN
      CALL PointProfile(ScalarC,Profile(:,Pos),Rho)
    ELSE  
      CALL PointProfile(ScalarC,Profile(:,Pos))
    END IF  
    ProfileName=TRIM(OutputC(Pos)%NameScalar)//'Point'
    CALL WriteProfile1(Profile(:,Pos),Domain%zP,ProfileName)
  END DO
  WRITE(*,*) 'Im in OutputPointprof End', myID

END SUBROUTINE OutputPointProfile

SUBROUTINE AreaProfiles(VecC,Profile)

  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind) :: Profile(:,:)

  TYPE(ScalarCell_T) :: ScalarC(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  CHARACTER*20  :: ProfileName
  INTEGER :: Pos,iProf
  
  iProf=0 
  IF (RhoPos>0) THEN
    CALL Assign(Rho,VecC,RhoPos)
  END IF  
  IF (uPosL>0.AND.RhoPos>0) THEN
    CALL Assign(ScalarC,VecC,uPosl)
    iProf=iProf+1
    CALL AreaProfile(ScalarC,Profile(:,iProf),Rho)
  END IF  
  IF (vPosL>0.AND.RhoPos>0) THEN
    CALL Assign(ScalarC,VecC,vPosl)
    iProf=iProf+1
    CALL AreaProfile(ScalarC,Profile(:,iProf),Rho)
  END IF  
  IF (wPosL>0.AND.RhoPos>0) THEN
    CALL Assign(ScalarC,VecC,wPosl)
    iProf=iProf+1
    CALL AreaProfile(ScalarC,Profile(:,iProf),Rho)
  END IF  
  DO Pos=1,SIZE(OutputC)
    CALL Assign(ScalarC,VectorOut,Pos)
    iProf=iProf+1
    IF (RhoPos>0.AND.OutputC(Pos)%DiffRho=='Rho') THEN
      CALL AreaProfile(ScalarC,Profile(:,iProf),Rho)
    ELSE
      CALL AreaProfile(ScalarC,Profile(:,iProf))
    END IF  
  END DO
END SUBROUTINE AreaProfiles

SUBROUTINE OutputAreaProfiles(Profile)

  REAL(RealKind) :: Profile(:,:)

  TYPE(ScalarCell_T) :: Rho(nbLoc)
  CHARACTER*20  :: ProfileName
  INTEGER :: Pos,iProf

  iProf=0
  IF (uPosL>0.AND.RhoPos>0) THEN
    ProfileName='UArea'
    iProf=iProf+1
    CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
  END IF  
  IF (vPosL>0.AND.RhoPos>0) THEN
    ProfileName='VArea'
    iProf=iProf+1
    CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
  END IF  
  IF (wPosL>0.AND.RhoPos>0) THEN
    ProfileName='WArea'
    iProf=iProf+1
    CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
  END IF  
  DO Pos=1,SIZE(OutputC)
    ProfileName=TRIM(OutputC(Pos)%NameScalar)//'Area'
    iProf=iProf+1
    CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
  END DO
END SUBROUTINE OutputAreaProfiles

SUBROUTINE VarianceProfiles(VecC,Profile)

  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind) :: Profile(:,:)

  TYPE(ScalarCell_T) :: ScalarC1(nbLoc)
  TYPE(ScalarCell_T) :: ScalarC2(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  CHARACTER*20  :: ProfileName
  INTEGER :: Pos,iProf
WRITE(*,*) 'Im in line 1'
  
  iProf=0 
  IF (RhoPos>0) THEN
    CALL Assign(Rho,VecC,RhoPos)

    WRITE(*,*) 'Im in line 2', uPosL
    IF (uPosL>0) THEN
    WRITE(*,*) 'Im in line 2aa', iProf

      CALL Assign(ScalarC1,VecC,uPosl)

      WRITE(*,*) 'Im in line 2aaa', iProf 
      DO Pos=1,SIZE(OutputC)
        CALL Assign(ScalarC2,VectorOut,Pos)
        iProf=iProf+1
        IF (RhoPos>0.AND.OutputC(Pos)%DiffRho=='Rho') THEN
          CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Rho,Profile(:,iProf))
        ELSE
          CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Profile=Profile(:,iProf))
        END IF  
WRITE(*,*) 'Im in line 2a', iProf

      END DO
      iProf=iProf+1
      CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC1,Rho,Profile(:,iProf))
      WRITE(*,*) 'Im in line 2b'
      IF (vPosL>0) THEN
        CALL Assign(ScalarC2,VecC,vPosl)
        iProf=iProf+1
        CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Rho,Profile(:,iProf))
      END IF  
      IF (wPosL>0) THEN
        CALL Assign(ScalarC2,VecC,wPosl)
        iProf=iProf+1
        CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Rho,Profile(:,iProf))
      END IF  
    END IF  

    WRITE(*,*) 'Im in line 3'

    IF (vPosL>0) THEN
      CALL Assign(ScalarC1,VecC,vPosl)


      WRITE(*,*) 'Im in line 4' 

      DO Pos=1,SIZE(OutputC)
        CALL Assign(ScalarC2,VectorOut,Pos)
        iProf=iProf+1
        IF (RhoPos>0.AND.OutputC(Pos)%DiffRho=='Rho') THEN
          CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Rho,Profile(:,iProf))
        ELSE
          CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Profile=Profile(:,iProf))
        END IF  
      END DO
      iProf=iProf+1
      CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC1,Rho,Profile(:,iProf))
      IF (wPosL>0) THEN
        CALL Assign(ScalarC2,VecC,wPosl)
        iProf=iProf+1
        CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Rho,Profile(:,iProf))
      END IF  
    END IF  
    IF (wPosL>0) THEN
      CALL Assign(ScalarC1,VecC,wPosl)
      DO Pos=1,SIZE(OutputC)
        CALL Assign(ScalarC2,VectorOut,Pos)
        iProf=iProf+1
        IF (RhoPos>0.AND.OutputC(Pos)%DiffRho=='Rho') THEN
          CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Rho,Profile(:,iProf))
        ELSE
          CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC2,Profile=Profile(:,iProf))
        END IF  
      END DO
      iProf=iProf+1
      CALL AreaVarianceProfile(ScalarC1,Rho,ScalarC1,Rho,Profile(:,iProf))
    END IF  
  END IF
END SUBROUTINE VarianceProfiles

SUBROUTINE OutputVarianceProfiles(Profile)

  REAL(RealKind) :: Profile(:,:)

  TYPE(ScalarCell_T) :: Rho(nbLoc)
  CHARACTER*20  :: ProfileName
  INTEGER :: Pos,iProf

  iProf=0

  WRITE(*,*) 'Im in line varB', myID


  IF (RhoPos>0) THEN
    IF (uPosL>0) THEN
      DO Pos=1,SIZE(OutputC)
        ProfileName='U'//TRIM(OutputC(Pos)%NameScalar)//'Var'
        iProf=iProf+1
        CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
      END DO

      WRITE(*,*) 'Im in line varB after U', myID


      ProfileName='UU'//'Var'
      iProf=iProf+1
      CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
      IF (vPosL>0) THEN
        ProfileName='UV'//'Var'
        iProf=iProf+1
        CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
      END IF  
      IF (wPosL>0) THEN
        ProfileName='UW'//'Var'
        iProf=iProf+1
        CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
      END IF  
    END IF  
    IF (vPosL>0) THEN
      DO Pos=1,SIZE(OutputC)
        ProfileName='V'//TRIM(OutputC(Pos)%NameScalar)//'Var'
        iProf=iProf+1
        CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
      END DO
      ProfileName='VV'//'Var'
      iProf=iProf+1
      CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
      IF (wPosL>0) THEN
        ProfileName='VW'//'Var'
        iProf=iProf+1
        CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)

         WRITE(*,*) 'Im in line varB after W', myID
      END IF  
    END IF  
    IF (wPosL>0) THEN
      DO Pos=1,SIZE(OutputC)
        ProfileName='W'//TRIM(OutputC(Pos)%NameScalar)//'Var'
        iProf=iProf+1
        CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
      END DO
      ProfileName='WW'//'Var'
      iProf=iProf+1
      CALL WriteProfile1(Profile(:,iProf),Domain%zP,ProfileName)
    END IF  
  END IF

   WRITE(*,*) 'Im in line varB after everyhing', myID
END SUBROUTINE OutputVarianceProfiles

SUBROUTINE PointProfile(ScalarC,Profile,Rho)
  TYPE(ScalarCell_T) :: ScalarC(:)
  REAL(RealKind) :: Profile(:)
  TYPE(ScalarCell_T), OPTIONAL :: Rho(:)

  INTEGER :: ix,iy,iz
  INTEGER :: igz
  REAL(RealKind) :: ProfileLoc(SIZE(Profile))

  ProfileLoc=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        IF (xP(ix-1)<=xProf.AND. &
            xP(ix  )> xProf.AND. &
            yP(iy-1)<=yProf.AND. &
            yP(iy  )> yProf) THEN
          DO iz=iz0+1,iz1
!           DO igz=iz**(-RefineZ+1),(iz+1)**(-RefineZ+1)-1 
            DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ)
              IF (PRESENT(Rho)) THEN
                ProfileLoc(igz)=ScalarC(ibLoc)%c(ix,iy,iz,1)/Rho(ibLoc)%c(ix,iy,iz,1)
              ELSE  
                ProfileLoc(igz)=ScalarC(ibLoc)%c(ix,iy,iz,1)
              END IF  
            END DO  
          END DO  
        END IF
      END DO
    END DO
  END DO  
  CALL MPI_Allreduce(ProfileLoc,Profile,SIZE(Profile),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
END SUBROUTINE PointProfile

SUBROUTINE AreaSurface(ScalarC,Values,Rho)
  TYPE(ScalarCell_T) :: ScalarC(:)
  REAL(RealKind) :: Values
  TYPE(ScalarCell_T), OPTIONAL :: Rho(:)

  INTEGER :: ix,iy !,iz
  INTEGER :: igz
  REAL(RealKind) :: xM,yM,Fac
  REAL(RealKind) :: ValuesLoc !
  REAL(RealKind) :: VolLoc !
  REAL(RealKind) :: Vol !

  ValuesLoc=Zero
  VolLoc=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      xM=Half*(xP(ix-1)+xP(ix))
      IF (xM<=xEnd.AND.xM>=xStart) THEN
        DO iy=iy0+1,iy1
          yM=Half*(yP(iy-1)+yP(iy))
          IF (yM<=yEnd.AND.yM>=yStart) THEN
                Fac=VolC(ix,iy,iz0+1) !*Domain%dz(igz)/dz(iz0) !
                VolLoc=VolLoc+Fac !
                IF (PRESENT(Rho)) THEN
                  ValuesLoc=ValuesLoc+Fac*ScalarC(ibLoc)%c(ix,iy,iz0+1,1)/Rho(ibLoc)%c(ix,iy,iz0+1,1)
                ELSE  
                  ValuesLoc=ValuesLoc+Fac*ScalarC(ibLoc)%c(ix,iy,iz0+1,1)
                END IF  
          END IF  
        END DO
      END IF
    END DO
  END DO  
  CALL MPI_Allreduce(ValuesLoc,Values,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  CALL MPI_Allreduce(VolLoc,Vol,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  Values=Values/(Vol+Eps)
END SUBROUTINE AreaSurface

SUBROUTINE AreaSurfaceBC(BoundCell,var,Val)
  TYPE(BoundCell_T) :: BoundCell(:)
  REAL(RealKind) :: Val
  INTEGER :: var

  INTEGER :: ix,iy,i,n,ixLoc,iyLoc
  INTEGER :: igz
  REAL(RealKind) :: xM,yM
  REAL(RealKind) :: ValLoc

  n=0
  Val=Zero
  ValLoc=Zero

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      xM=Half*(xP(ix-1)+xP(ix))
      IF (xM<=xEnd.AND.xM>=xStart) THEN 
        DO iy=iy0+1,iy1
          yM=Half*(yP(iy-1)+yP(iy))
          IF (yM<=yEnd.AND.yM>=yStart) THEN 
            DO i=1,NumBoundCell
              ixLoc=Floor(ib)%BoundCell(i)%ix
              iyLoc=Floor(ib)%BoundCell(i)%iy
              IF (ixLoc==ix.AND.iyLoc==iy) THEN
                IF (var==1) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%RadDirekt  
                  n=n+1
                ELSE IF (var==2) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%RadDiffus  
                  n=n+1
                ELSE IF (var==3) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%RadInfred  
                  n=n+1
                ELSE IF (var==4) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%FluxSens  
                  n=n+1
                ELSE IF (var==5) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%FluxLat  
                  n=n+1
                ELSE IF (var==6) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%DragM  
                  n=n+1
                ELSE IF (var==7) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%DragH  
                  n=n+1
                ELSE IF (var==8) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%DragQ  
                  n=n+1
                ELSE IF (var==9) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%TSoil  
                  n=n+1
                ELSE IF (var==10) THEN
                  ValLoc=ValLoc+Floor(ib)%BoundCell(i)%QVSoil  
                  n=n+1
                ELSE
                  WRITE(*,*) 'var not properly defined!'
                  EXIT
                END IF
              END IF
            END DO  
          END IF
        END DO
      END IF
    END DO
  END DO  
  CALL MPI_Allreduce(ValLoc,Val,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  Val=Val/(Real(n)+Eps)
END SUBROUTINE AreaSurfaceBC

SUBROUTINE AreaProfile(ScalarC,Profile,Rho,Subdomain,oxs,oxe,oys,oye)
  TYPE(ScalarCell_T) :: ScalarC(:)
  REAL(RealKind) :: Profile(:)
  TYPE(ScalarCell_T), OPTIONAL :: Rho(:)
  LOGICAL, OPTIONAL :: Subdomain
  REAL(RealKind), OPTIONAL :: oxs,oxe,oys,oye

  INTEGER :: ix,iy,iz
  INTEGER :: igz
  REAL(RealKind) :: xM,yM,Fac
  REAL(RealKind) :: ProfileLoc(SIZE(Profile))
  REAL(RealKind) :: VolLoc(SIZE(Profile))
  REAL(RealKind) :: Vol(SIZE(Profile))

  ProfileLoc=Zero
  VolLoc=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (PRESENT(Subdomain).AND.Subdomain) THEN
      xStart=oxs
      xEnd=oxe
      yStart=oys
      yEnd=oye
    ELSE
      xStart=xP(ix0)
      xEnd=xP(ix1)
      yStart=yP(iy0)
      yEnd=yP(iy1)
    END IF
    DO ix=ix0+1,ix1
      xM=Half*(xP(ix-1)+xP(ix))
      IF (xM<=xEnd.AND.xM>=xStart) THEN
        DO iy=iy0+1,iy1
          yM=Half*(yP(iy-1)+yP(iy))
          IF (yM<=yEnd.AND.yM>=yStart) THEN
            DO iz=iz0+1,iz1
!             DO igz=iz**(-RefineZ+1),(iz+1)**(-RefineZ+1)-1 
              DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ)
                Fac=VolC(ix,iy,iz)*Domain%dz(igz)/dz(iz)
                VolLoc(igz)=VolLoc(igz)+Fac
                IF (PRESENT(Rho)) THEN
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*ScalarC(ibLoc)%c(ix,iy,iz,1)/Rho(ibLoc)%c(ix,iy,iz,1)
                ELSE  
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*ScalarC(ibLoc)%c(ix,iy,iz,1)
                END IF  
              END DO  
            END DO  
          END IF  
        END DO
      END IF
    END DO
  END DO  
  CALL MPI_Allreduce(ProfileLoc,Profile,SIZE(Profile),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  CALL MPI_Allreduce(VolLoc,Vol,SIZE(Profile),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  Profile=Profile/(Vol+Eps)
END SUBROUTINE AreaProfile

SUBROUTINE AreaTotalHor(ScalarC,Total,Rho)
  TYPE(ScalarCell_T) :: ScalarC(:)
  REAL(RealKind) :: Total(Domain%igx0+1:Domain%igx1,Domain%igy0+1:Domain%igy1)
  TYPE(ScalarCell_T), OPTIONAL :: Rho(:)

  INTEGER :: ix,iy,iz
  INTEGER :: igx,igy,igz
  REAL(RealKind) :: xM,yM,Fac
  REAL(RealKind) :: TotalLoc(Domain%igx0+1:Domain%igx1,Domain%igy0+1:Domain%igy1)

  TotalLoc=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO igx=ix**(-RefineX+1),(ix+1)**(-RefineX+1)-1 
        xM=Half*(Domain%xP(igx-1)+Domain%xP(igx))
        IF (xM<=xEnd.AND.xM>=xStart) THEN
          DO iy=iy0+1,iy1
            DO igy=iy**(-RefineY+1),(iy+1)**(-RefineY+1)-1 
              yM=Half*(Domain%yP(igy-1)+Domain%yP(igy))
              IF (yM<=yEnd.AND.yM>=yStart) THEN
                DO iz=iz0+1,iz1
!                 DO igz=iz**(-RefineZ+1),(iz+1)**(-RefineZ+1)-1 
                  DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ)
                    Fac=VolC(ix,iy,iz)*(Domain%dx(igx)*Domain%dy(igy)*Domain%dz(igz))/(dx(ix)*dy(iy)*dz(iz))
                    IF (PRESENT(Rho)) THEN
                      TotalLoc(igx,igy)=TotalLoc(igx,igy)+Fac*ScalarC(ibLoc)%c(ix,iy,iz,1)/Rho(ibLoc)%c(ix,iy,iz,1)
                    ELSE  
                     TotalLoc(igx,igy)=TotalLoc(igx,igy)+Fac*ScalarC(ibLoc)%c(ix,iy,iz,1)
                    END IF  
                  END DO  
                END DO  
              END IF  
            END DO  
          END DO  
        END IF
      END DO
    END DO
  END DO  
  CALL MPI_Allreduce(TotalLoc,Total,SIZE(Total),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
END SUBROUTINE AreaTotalHor

SUBROUTINE AreaVarianceProfile(ScalarC1,Rho1,ScalarC2,Rho2,Profile,Subdomain,oxs,oxe,oys,oye)
  TYPE(ScalarCell_T) :: ScalarC1(:)
  TYPE(ScalarCell_T), OPTIONAL :: Rho1(:)
  TYPE(ScalarCell_T) :: ScalarC2(:)
  TYPE(ScalarCell_T), OPTIONAL :: Rho2(:)
  REAL(RealKind) :: Profile(:)
  LOGICAL, OPTIONAL :: Subdomain
  REAL(RealKind), OPTIONAL :: oxs,oxe,oys,oye

  INTEGER :: ix,iy,iz
  INTEGER :: igz
  REAL(RealKind) :: xM,yM,Fac
  REAL(RealKind) :: Profile1(SIZE(Profile)),Profile2(SIZE(Profile))
  REAL(RealKind) :: ProfileLoc(SIZE(Profile))
  REAL(RealKind) :: VolLoc(SIZE(Profile))
  REAL(RealKind) :: Vol(SIZE(Profile))

  IF (PRESENT(Rho1)) THEN
    CALL AreaProfile(ScalarC1,Profile1,Rho1)
  ELSE  
    CALL AreaProfile(ScalarC1,Profile1)
  END IF  
  IF (PRESENT(Rho2)) THEN
    CALL AreaProfile(ScalarC2,Profile2,Rho2)
  ELSE  
    CALL AreaProfile(ScalarC2,Profile2)
  END IF  
  ProfileLoc=Zero
  VolLoc=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (PRESENT(Subdomain).AND.Subdomain) THEN
      xStart=oxs
      xEnd=oxe
      yStart=oys
      yEnd=oye
    ELSE
      xStart=xP(ix0)
      xEnd=xP(ix1)
      yStart=yP(iy0)
      yEnd=yP(iy1)
    END IF
    DO ix=ix0+1,ix1
      xM=Half*(xP(ix-1)+xP(ix))
      IF (xM<=xEnd.AND.xM>=xStart) THEN
        DO iy=iy0+1,iy1
          yM=Half*(yP(iy-1)+yP(iy))
          IF (yM<=yEnd.AND.yM>=yStart) THEN
            DO iz=iz0+1,iz1
!             DO igz=iz*2**(-RefineZ+1),(iz+1)*2**(-RefineZ+1)-1 
              DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ)
                Fac=VolC(ix,iy,iz)*Domain%dz(igz)/dz(iz)
                VolLoc(igz)=VolLoc(igz)+Fac
                IF (PRESENT(Rho1).AND.PRESENT(Rho2)) THEN
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*(ScalarC1(ibLoc)%c(ix,iy,iz,1)/Rho1(ibLoc)%c(ix,iy,iz,1)-Profile1(igz)) &
                                                     *(ScalarC2(ibLoc)%c(ix,iy,iz,1)/Rho2(ibLoc)%c(ix,iy,iz,1)-Profile2(igz))
                ELSE IF (PRESENT(Rho1)) THEN
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*(ScalarC1(ibLoc)%c(ix,iy,iz,1)/Rho1(ibLoc)%c(ix,iy,iz,1)-Profile1(igz)) &
                                                     *(ScalarC2(ibLoc)%c(ix,iy,iz,1)-Profile2(igz))
                ELSE IF (PRESENT(Rho2)) THEN
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*(ScalarC1(ibLoc)%c(ix,iy,iz,1)-Profile1(igz)) &
                                                     *(ScalarC2(ibLoc)%c(ix,iy,iz,1)/Rho2(ibLoc)%c(ix,iy,iz,1)-Profile2(igz))
                ELSE  
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*(ScalarC1(ibLoc)%c(ix,iy,iz,1)-Profile1(igz)) &
                                                     *(ScalarC2(ibLoc)%c(ix,iy,iz,1)-Profile2(igz))
                END IF  
              END DO  
            END DO  
          END IF  
        END DO
      END IF
    END DO
  END DO  
  CALL MPI_Allreduce(ProfileLoc,Profile,SIZE(Profile),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  CALL MPI_Allreduce(VolLoc,Vol,SIZE(Profile),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  Profile=Profile/(Vol+Eps)
END SUBROUTINE AreaVarianceProfile

SUBROUTINE AreaVarianceTotalProfile(ScalarC1,Rho1,ScalarC2,Rho2,Profile,Subdomain,oxs,oxe,oys,oye)
  TYPE(ScalarCell_T) :: ScalarC1(:)
  TYPE(ScalarCell_T), OPTIONAL :: Rho1(:)
  TYPE(ScalarCell_T) :: ScalarC2(:)
  TYPE(ScalarCell_T), OPTIONAL :: Rho2(:)
  REAL(RealKind) :: Profile(:)
  LOGICAL, OPTIONAL :: Subdomain
  REAL(RealKind), OPTIONAL :: oxs,oxe,oys,oye

  INTEGER :: ix,iy,iz
  INTEGER :: igz
  REAL(RealKind) :: xM,yM,Fac
  REAL(RealKind) :: ProfileLoc(SIZE(Profile))
  REAL(RealKind) :: VolLoc(SIZE(Profile))
  REAL(RealKind) :: Vol(SIZE(Profile))

  ProfileLoc=Zero
  VolLoc=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (PRESENT(Subdomain).AND.Subdomain) THEN
      xStart=oxs
      xEnd=oxe
      yStart=oys
      yEnd=oye
    ELSE
      xStart=xP(ix0)
      xEnd=xP(ix1)
      yStart=yP(iy0)
      yEnd=yP(iy1)
    END IF
    DO ix=ix0+1,ix1
      xM=Half*(xP(ix-1)+xP(ix))
      IF (xM<=xEnd.AND.xM>=xStart) THEN 
        DO iy=iy0+1,iy1
          yM=Half*(yP(iy-1)+yP(iy))
          IF (yM<=yEnd.AND.yM>=yStart) THEN 
            DO iz=iz0+1,iz1
!             DO igz=iz*2**(-RefineZ+1),(iz+1)*2**(-RefineZ+1)-1 
              DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ)
                Fac=VolC(ix,iy,iz)*Domain%dz(igz)/dz(iz)
                VolLoc(igz)=VolLoc(igz)+Fac
                IF (PRESENT(Rho1).AND.PRESENT(Rho2)) THEN 
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*ScalarC1(ibLoc)%c(ix,iy,iz,1)/Rho1(ibLoc)%c(ix,iy,iz,1) &
                                                     *ScalarC2(ibLoc)%c(ix,iy,iz,1)/Rho2(ibLoc)%c(ix,iy,iz,1)
                ELSE IF (PRESENT(Rho1)) THEN
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*ScalarC1(ibLoc)%c(ix,iy,iz,1)/Rho1(ibLoc)%c(ix,iy,iz,1) &
                                                     *ScalarC2(ibLoc)%c(ix,iy,iz,1)
                ELSE IF (PRESENT(Rho2)) THEN
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*ScalarC1(ibLoc)%c(ix,iy,iz,1) &
                                                     *ScalarC2(ibLoc)%c(ix,iy,iz,1)/Rho2(ibLoc)%c(ix,iy,iz,1)
                ELSE
                  ProfileLoc(igz)=ProfileLoc(igz)+Fac*ScalarC1(ibLoc)%c(ix,iy,iz,1) &
                                                     *ScalarC2(ibLoc)%c(ix,iy,iz,1)
                END IF
              END DO
            END DO
          END IF
        END DO
      END IF
    END DO
  END DO
  CALL MPI_Allreduce(ProfileLoc,Profile,SIZE(Profile),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  CALL MPI_Allreduce(VolLoc,Vol,SIZE(Profile),MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)
  Profile=Profile/(Vol+Eps)
END SUBROUTINE AreaVarianceTotalProfile

SUBROUTINE WriteProfile1(Profile,zP,ProfileName)

  REAL(RealKind) :: Profile(:)
  REAL(RealKind) :: zP(:)
  CHARACTER(*) :: ProfileName

  INTEGER :: i

  IF (MyId==0) THEN
    OPEN(UNIT=OutputUnit,FILE=TRIM(ProfileName)//'.prof',STATUS='UNKNOWN')
    DO i=1,SIZE(Profile)
      WRITE(OutputUnit,*) Half*(zP(i)+zP(i+1)),Profile(i)
    END DO  
    CLOSE(OutputUnit)
  END IF
END SUBROUTINE WriteProfile1

SUBROUTINE MakeProfile(ScalarC,VecC,VecPos,it,ic,GnuplotOut,IsMet,Island) ! Vogelsberg ! MJ

  TYPE(ScalarCell_T) :: ScalarC(:)
  TYPE(Vector4Cell_T) :: VecC(:)
!  TYPE(Vector4Cell_T) :: VecT(:)
  REAL(RealKind), ALLOCATABLE :: GnuplotOut(:,:)
  LOGICAL, OPTIONAL :: IsMet
  LOGICAL, OPTIONAL :: Island

  REAL(RealKind)              :: SUMM(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM2(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMMVol(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM_ABWEICH(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM_ABWEICH2(Domain%iz0+1:Domain%iz1)

  REAL(RealKind)              :: SUMM_ISLAND(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM2_ISLAND(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMMVol_ISLAND(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM_ABWEICH_ISLAND(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM_ABWEICH2_ISLAND(Domain%iz0+1:Domain%iz1)

  REAL(RealKind)              :: SUMM_OCEAN(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM2_OCEAN(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMMVol_OCEAN(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM_ABWEICH_OCEAN(Domain%iz0+1:Domain%iz1)
  REAL(RealKind)              :: SUMM_ABWEICH2_OCEAN(Domain%iz0+1:Domain%iz1)

  INTEGER                     :: ix,iy,iz
  INTEGER                     :: ixProf,iyProf
  INTEGER                     :: ixStart,iyStart
  INTEGER                     :: ixEnd,iyEnd
  INTEGER                     :: i,ibProf,ibProfLoc,k
  INTEGER                     :: ic,VecPos,it
  INTEGER                     :: SumArea,SumCBH,SumCTH
  INTEGER                     :: NCells2D,CloudyCols
  LOGICAL                     :: CloudyCol
  Real(RealKind)              :: VolLoc

  CALL StartEndPoints(ixStart,ixEnd,iyStart,iyEnd,ixProf,iyProf,ibProf)

  IF (ProfAverage) THEN
    IF (ic==1) WRITE(*,*) 'ixStart',ixStart,'iyStart',iyStart
    IF (ic==1) WRITE(*,*) 'ixEnd',ixEnd,'iyEnd',iyEnd
  ELSE
    IF (ic==1) WRITE(*,*) 'ixProf',ixProf,'iyProf',iyProf
  END IF

    IF (ProfAverage) THEN
    IF (.NOT.PRESENT(Island)) THEN ! Domain Average
      SumArea=0
      SUMM = zero
      SUMM2 = zero
      SUMMVol = zero
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            DO iz=iz0+1,iz1
              IF ((ix>=ixStart.AND.ix<=ixEnd) &
              .AND.(iy>=iyStart.AND.iy<=iyEnd)) THEN
                IF (VecPos==9999) THEN
                  SUMM(iz) = SUMM(iz) +ScalarC(ibLoc)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                ELSE
                  IF(ic==3) THEN
                  SUMM(iz) = SUMM(iz) +VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                  SUMM2(iz) = SUMM2(iz) +VecC(ibLoc)%Vec(thPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                  ELSE
                  SUMM(iz) = SUMM(iz) +VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                  END IF
                END IF
                SumArea=SumArea+1
                SUMMVol(iz)=SUMMVol(iz) + VolC(ix,iy,iz)
              ELSE
              END IF
            END DO
          END DO
        END DO
      END DO
        SUMM = SUMM/SUMMVol 
        SUMM2 = SUMM2/SUMMVol 
        GnuplotOut(ic,:)=SUMM
      IF (Deviation) THEN
        SUMM_ABWEICH = zero
        SUMM_ABWEICH2 = zero
        SUMMVol = zero
        DO ibLoc=1,nbLoc
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          DO ix=ix0+1,ix1
            DO iy=iy0+1,iy1
              DO iz=iz0+1,iz1
                IF ((ix>=ixStart.AND.ix<=ixEnd) &
                .AND.(iy>=iyStart.AND.iy<=iyEnd)) THEN
                  IF (VecPos==9999) THEN
                    SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+(ABS(SUMM(iz)-ScalarC(ibLoc)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                  ELSE
                    IF(ic==3) THEN
!                    SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+((SUMM(iz)-VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                    SUMM_ABWEICH2(iz)=SUMM_ABWEICH2(iz)+((((VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)-SUMM(iz)))*&
                    & ((VecC(ibLoc)%Vec(thPos)%c(ix,iy,iz,it))-SUMM2(iz)))*VolC(ix,iy,iz))
                    ELSE
                    SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+(ABS(SUMM(iz)-VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                    END IF
                  END IF
                SUMMVol(iz)=SUMMVol(iz) + VolC(ix,iy,iz)
                ELSE
                END IF
              END DO
            END DO
          END DO
        END DO
        IF(ic==3) THEN
        SUMM_ABWEICH2 = SUMM_ABWEICH2/SUMMVol 
        GnuplotOut(ic+VectorComponentsOut,:)=SUMM_ABWEICH2
        ELSE
        SUMM_ABWEICH = SUMM_ABWEICH/SUMMVol 
        GnuplotOut(ic+VectorComponentsOut,:)=SUMM_ABWEICH
        END IF
      END IF
    ELSE IF (Island) THEN  ! Island Average 
      SumArea=0
      SUMM = zero
      SUMM2 = zero
      SUMMVol = zero
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
            IF (VolC(ix,iy,iz0+1)<VolLoc) THEN  ! Only apply loop for island colum
              DO iz=iz0+1,iz1
                IF ((ix>=ixStart.AND.ix<=ixEnd) &
                .AND.(iy>=iyStart.AND.iy<=iyEnd)) THEN
                  IF (VecPos==9999) THEN
                    SUMM(iz) = SUMM(iz) +ScalarC(ibLoc)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                  ELSE
                    IF(ic==3) THEN
                    SUMM(iz) = SUMM(iz) +VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                    SUMM2(iz) = SUMM2(iz) +VecC(ibLoc)%Vec(thPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                    ELSE
                    SUMM(iz) = SUMM(iz) +VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                    END IF
                  END IF
                  SumArea=SumArea+1
                  SUMMVol(iz)=SUMMVol(iz) + VolC(ix,iy,iz)
                ELSE
                END IF
              END DO
            END IF
          END DO
        END DO
      END DO
        SUMM = SUMM/SUMMVol 
        SUMM2 = SUMM2/SUMMVol 
        GnuplotOut(ic,:)=SUMM
      IF (Deviation) THEN
        SUMM_ABWEICH = zero
        SUMM_ABWEICH2 = zero
        SUMMVol = zero
        DO ibLoc=1,nbLoc
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          DO ix=ix0+1,ix1
            DO iy=iy0+1,iy1
              VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
              IF (VolC(ix,iy,iz0+1)<VolLoc) THEN  ! Only apply loop for island colum
                DO iz=iz0+1,iz1
                  IF ((ix>=ixStart.AND.ix<=ixEnd) &
                  .AND.(iy>=iyStart.AND.iy<=iyEnd)) THEN
                    IF (VecPos==9999) THEN
                      SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+(ABS(SUMM(iz)-ScalarC(ibLoc)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                    ELSE
                      IF(ic==3) THEN
!                      SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+((SUMM(iz)-VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                      SUMM_ABWEICH2(iz)=SUMM_ABWEICH2(iz)+((((VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)-&
                      & SUMM(iz)))*((VecC(ibLoc)%Vec(thPos)%c(ix,iy,iz,it))-SUMM2(iz)))*VolC(ix,iy,iz))
                      ELSE
                      SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+(ABS(SUMM(iz)-VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                      END IF
                    END IF
                  SUMMVol(iz)=SUMMVol(iz) + VolC(ix,iy,iz)
                  ELSE
                  END IF
                END DO
              END IF
            END DO
          END DO
        END DO
        IF(ic==3) THEN
        SUMM_ABWEICH2 = SUMM_ABWEICH2/SUMMVol 
        GnuplotOut(ic+VectorComponentsOut,:)=SUMM_ABWEICH2
        ELSE
        SUMM_ABWEICH = SUMM_ABWEICH/SUMMVol 
        GnuplotOut(ic+VectorComponentsOut,:)=SUMM_ABWEICH
        END IF
      END IF

    ELSE ! Ocean Average
      SumArea=0
      SUMM = zero
      SUMM2 = zero
      SUMMVol = zero
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
            IF (.NOT.(VolC(ix,iy,iz0+1)<VolLoc)) THEN  ! Only apply loop for ocean colum
              DO iz=iz0+1,iz1
                IF ((ix>=ixStart.AND.ix<=ixEnd) &
                .AND.(iy>=iyStart.AND.iy<=iyEnd)) THEN
                  IF (VecPos==9999) THEN
                    SUMM(iz) = SUMM(iz) +ScalarC(ibLoc)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                  ELSE
                    IF(ic==3) THEN
                    SUMM(iz) = SUMM(iz) +VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                    SUMM2(iz) = SUMM2(iz) +VecC(ibLoc)%Vec(thPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                    ELSE
                    SUMM(iz) = SUMM(iz) +VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)*VolC(ix,iy,iz)
                    END IF
                  END IF
                  SumArea=SumArea+1
                  SUMMVol(iz)=SUMMVol(iz) + VolC(ix,iy,iz)
                ELSE
                END IF
              END DO
            END IF
          END DO
        END DO
      END DO
        SUMM = SUMM/SUMMVol 
        SUMM2 = SUMM2/SUMMVol 
        GnuplotOut(ic,:)=SUMM
      IF (Deviation) THEN
        SUMM_ABWEICH = zero
        SUMM_ABWEICH2 = zero
        SUMMVol = zero
        DO ibLoc=1,nbLoc
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          DO ix=ix0+1,ix1
            DO iy=iy0+1,iy1
              VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
              IF (.NOT.(VolC(ix,iy,iz0+1)<VolLoc)) THEN  ! Only apply loop for ocean colum
                DO iz=iz0+1,iz1
                  IF ((ix>=ixStart.AND.ix<=ixEnd) &
                  .AND.(iy>=iyStart.AND.iy<=iyEnd)) THEN
                    IF (VecPos==9999) THEN
                      SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+(ABS(SUMM(iz)-ScalarC(ibLoc)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                    ELSE
                      IF(ic==3) THEN
!                      SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+((SUMM(iz)-VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                      SUMM_ABWEICH2(iz)=SUMM_ABWEICH2(iz)+((((VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)-&
                      & SUMM(iz)))*((VecC(ibLoc)%Vec(thPos)%c(ix,iy,iz,it))-SUMM2(iz)))*VolC(ix,iy,iz))
                      ELSE
                      SUMM_ABWEICH(iz)=SUMM_ABWEICH(iz)+(ABS(SUMM(iz)-VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)))*VolC(ix,iy,iz)
                      END IF
                    END IF
                  SUMMVol(iz)=SUMMVol(iz) + VolC(ix,iy,iz)
                  ELSE
                  END IF
                END DO
              END IF
            END DO
          END DO
        END DO
        IF(ic==3) THEN
        SUMM_ABWEICH2 = SUMM_ABWEICH2/SUMMVol 
        GnuplotOut(ic+VectorComponentsOut,:)=SUMM_ABWEICH2
        ELSE
        SUMM_ABWEICH = SUMM_ABWEICH/SUMMVol 
        GnuplotOut(ic+VectorComponentsOut,:)=SUMM_ABWEICH
        END IF
      END IF
    END IF ! .NOT. ProfIslandOcean
    ELSE   ! ProfAverage

    DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (ProfOut) THEN
      IF (.NOT.(ProfAverage)) THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            IF (xP(ix-1)<=xProf.AND. &
                xP(ix  )> xProf.AND. &
                yP(iy-1)<=yProf.AND. &
                yP(iy  )> yProf) THEN
            ixProf=ix
            iyProf=iy
            ibProf=ib
            ibProfLoc=ibLoc
              IF (VecPos==9999) THEN
                GnuplotOut(ic,iz) = ScalarC(ibProfLoc)%c(ixProf,iyProf,iz,it)
              ELSE
                GnuplotOut(ic,iz) = VecC(ibProfLoc)%Vec(VecPos)%c(ixProf,iyProf,iz,it)
              END IF
            END IF
          END DO
        END DO
      END DO
      ENDIF
    ENDIF
    ENDDO
    END IF
    
    IF (ProfOut.AND.ProfAverage.AND.VecPos==RhoCPos.AND.PRESENT(IsMet)) THEN
    IF (.NOT.PRESENT(Island)) THEN ! Domain Average
      NCells2D=0
      CloudyCols=0
      CloudCover=Zero
      SumCBH=Zero
      SumCTH=Zero
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            NCells2D=NCells2D+1
            CloudyCol=.FALSE.
            DO iz=iz0+1,iz1
              IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5.AND..NOT.CloudyCol) THEN
                CloudyCol=.TRUE.
                CloudyCols=CloudyCols+1
              END IF
            END DO
          END DO
        END DO
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            DO iz=iz0+1,iz1
              IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5) THEN
                SumCBH=SumCBH+0.5d0*(zP(iz-1)+zP(iz))             
                EXIT
              END IF
            END DO
          END DO
        END DO
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            DO iz=iz1,iz0+1,-1
              IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5) THEN
                SumCTH=SumCTH+0.5d0*(zP(iz-1)+zP(iz))             
                EXIT
              END IF
            END DO
          END DO
        END DO
      END DO
      CloudCover=REAL(CloudyCols)/REAL(NCells2D)
      CloudBaseHeight=SumCBH/REAL(CloudyCols+Eps)
      CloudTopHeight=SumCTH/REAL(CloudyCols+Eps)
    ELSE IF (Island) THEN ! Island average
      NCells2D=0
      CloudyCols=0
      CloudCover2=Zero
      SumCBH=Zero
      SumCTH=Zero
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
            IF (VolC(ix,iy,iz0+1)<VolLoc) THEN  ! Only apply loop for island colum
              CloudyCol=.FALSE.
              NCells2D=NCells2D+1
              DO iz=iz0+1,iz1
                IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5.AND..NOT.CloudyCol) THEN
                  CloudyCol=.TRUE.
                  CloudyCols=CloudyCols+1
                END IF
              END DO
            END IF
          END DO
        END DO
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
            IF (VolC(ix,iy,iz0+1)<VolLoc) THEN  ! Only apply loop for island colum
              DO iz=iz0+1,iz1
                IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5) THEN
                  SumCBH=SumCBH+0.5d0*(zP(iz-1)+zP(iz))             
                  EXIT
                END IF
              END DO
            END IF
          END DO
        END DO
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
            IF (VolC(ix,iy,iz0+1)<VolLoc) THEN  ! Only apply loop for island colum
              DO iz=iz1,iz0+1,-1
                IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5) THEN
                  SumCTH=SumCTH+0.5d0*(zP(iz-1)+zP(iz))             
                  EXIT
                END IF
              END DO
            END IF
          END DO
        END DO
      END DO
      CloudCover2=REAL(CloudyCols)/REAL(NCells2D)
      CloudBaseHeight2=SumCBH/REAL(CloudyCols+Eps)
      CloudTopHeight2=SumCTH/REAL(CloudyCols+Eps)
    ELSE ! Ocean average
      NCells2D=0
      CloudyCols=0
      CloudCover3=Zero
      SumCBH=Zero
      SumCTH=Zero
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
            IF (.NOT.(VolC(ix,iy,iz0+1)<VolLoc)) THEN  ! Only apply loop for ocean colum
              CloudyCol=.FALSE.
              NCells2D=NCells2D+1
              DO iz=iz0+1,iz1
                IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5.AND..NOT.CloudyCol) THEN
                  CloudyCol=.TRUE.
                  CloudyCols=CloudyCols+1
                END IF
              END DO
            END IF
          END DO
        END DO
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
            IF (.NOT.(VolC(ix,iy,iz0+1)<VolLoc)) THEN  ! Only apply loop for ocean colum
              DO iz=iz0+1,iz1
                IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5) THEN
                  SumCBH=SumCBH+0.5d0*(zP(iz-1)+zP(iz))             
                  EXIT
                END IF
              END DO
            END IF
          END DO
        END DO
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            VolLoc=dx(ix)*dy(iy)*dz(iz0+1)
            IF (.NOT.(VolC(ix,iy,iz0+1)<VolLoc)) THEN  ! Only apply loop for ocean colum
              DO iz=iz1,iz0+1,-1
                IF (VecC(ibLoc)%Vec(VecPos)%c(ix,iy,iz,it)>1.0d-5) THEN
                  SumCTH=SumCTH+0.5d0*(zP(iz-1)+zP(iz))             
                  EXIT
                END IF
              END DO
            END IF
          END DO
        END DO
      END DO
      CloudCover3=REAL(CloudyCols)/REAL(NCells2D)
      CloudBaseHeight3=SumCBH/REAL(CloudyCols+Eps)
      CloudTopHeight3=SumCTH/REAL(CloudyCols+Eps)
    END IF
    END IF 

END SUBROUTINE MakeProfile

SUBROUTINE MakeSoil(VecC,VecPos,icB,SoilOut) 

  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind), ALLOCATABLE :: SoilOut(:,:)
  REAL(RealKind)              :: SUMM(1:Domain%nrsoillayers)
  INTEGER                     :: ix,iy,iz
  INTEGER                     :: ixProf,iyProf
  INTEGER                     :: ixStart,iyStart
  INTEGER                     :: ixEnd,iyEnd
  INTEGER                     :: i,ibProf,ibProfLoc,k
  INTEGER                     :: icB,VecPos,it
  INTEGER                     :: SumArea

  CALL StartEndPoints(ixStart,ixEnd,iyStart,iyEnd,ixProf,iyProf,ibProf)

    DO iz=1,Domain%nrsoillayers
      IF (ProfAverage) THEN
        SumArea=0
        SUMM(iz) = zero
        DO ibLoc=1,nbLoc
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          DO i=1,Floor(ib)%NumBoundCell
            IF ((Floor(ib)%BoundCell(i)%ix>=ixStart.AND.Floor(ib)%BoundCell(i)%ix<=ixEnd) &
            .AND.(Floor(ib)%BoundCell(i)%iy>=iyStart.AND.Floor(ib)%BoundCell(i)%iy<=iyEnd)) THEN
              SUMM(iz) = SUMM(iz) + VecC(ibLoc)%Vec(VecPos)%cB(i,iz)
              SumArea=SumArea+1
            ELSE
              SUMM(iz)=SUMM(iz)
              SumArea=SumArea
            END IF
          END DO
        END DO
        SUMM(iz) = SUMM(iz)/SumArea
        SoilOut(icB,iz)=SUMM(iz)
      ELSE
        CALL Set(Floor(ibProf))
        DO i=1,Floor(ibProf)%NumBoundCell
          IF ((Floor(ibProf)%BoundCell(i)%ix==ixProf) &
            .AND.(Floor(ibProf)%BoundCell(i)%iy==iyProf)) THEN 
            SoilOut(icB,iz) = VecC(ibProfLoc)%Vec(VecPos)%cB(i,iz)
          END IF
        END DO
      END IF
    END DO

END SUBROUTINE MakeSoil

SUBROUTINE WriteProfile(GnuplotOut,ic,ProfCase,ProfComp,IsMet,Island) ! MJ

  REAL(RealKind), ALLOCATABLE :: GnuplotOut(:,:)
  CHARACTER(200)              :: ProfComp
  CHARACTER(20)               :: ProfCase
  LOGICAL, OPTIONAL           :: IsMet
  LOGICAL, OPTIONAL           :: Island
  INTEGER                     :: ic,i,iz

  IF (.NOT.PRESENT(Island)) THEN
    SELECT CASE(ProfCase)
    CASE('Prof')
      WRITE(OutputUnit,*) nz,2
      WRITE(OutputUnit,*) 'Equal'
    CASE('GnuProf')
      WRITE(OutputUnit,*) ProfComp
      WRITE(*,*) 'ProfComp: ', ProfComp
    END SELECT

    SELECT CASE(ProfCase)
    CASE('Prof')
      DO iz=iz0+1,iz1
        WRITE(OutputUnit,*)  &
              0.5d0*(zP(iz-1)+zP(iz)) &
              ,GnuplotOut(ic,iz)
      END DO
    CASE('GnuProf')
      DO iz=iz0+1,iz1
        WRITE(OutputUnit,*)  &
              0.5d0*(zP(iz-1)+zP(iz)) &
              ,GnuplotOut(1:ic,iz)
      END DO
      IF(RhoCPos>0.AND.CloudCover<=One) WRITE(OutputUnit,*) CloudCover,CloudBaseHeight,CloudTopHeight
    END SELECT
  ELSE IF (Island) THEN
    SELECT CASE(ProfCase)
    CASE('Prof')
      WRITE(OutputUnit+1,*) nz,2
      WRITE(OutputUnit+1,*) 'Equal'
    CASE('GnuProf')
      WRITE(OutputUnit+1,*) ProfComp
      WRITE(*,*) 'ProfComp: ', ProfComp
    END SELECT

    SELECT CASE(ProfCase)
    CASE('Prof')
      DO iz=iz0+1,iz1
        WRITE(OutputUnit+1,*)  &
              0.5d0*(zP(iz-1)+zP(iz)) &
              ,GnuplotOut(ic,iz)
      END DO
    CASE('GnuProf')
      DO iz=iz0+1,iz1
        WRITE(OutputUnit+1,*)  &
              0.5d0*(zP(iz-1)+zP(iz)) &
              ,GnuplotOut(1:ic,iz)
      END DO
      IF(RhoCPos>0.AND.CloudCover2<=One) WRITE(OutputUnit+1,*) CloudCover2,CloudBaseHeight2,CloudTopHeight2
    END SELECT
  ELSE ! Ocean
    SELECT CASE(ProfCase)
    CASE('Prof')
      WRITE(OutputUnit+2,*) nz,2
      WRITE(OutputUnit+2,*) 'Equal'
    CASE('GnuProf')
      WRITE(OutputUnit+2,*) ProfComp
      WRITE(*,*) 'ProfComp: ', ProfComp
    END SELECT

    SELECT CASE(ProfCase)
    CASE('Prof')
      DO iz=iz0+1,iz1
        WRITE(OutputUnit+2,*)  &
              0.5d0*(zP(iz-1)+zP(iz)) &
              ,GnuplotOut(ic,iz)
      END DO
    CASE('GnuProf')
      DO iz=iz0+1,iz1
        WRITE(OutputUnit+2,*)  &
              0.5d0*(zP(iz-1)+zP(iz)) &
              ,GnuplotOut(1:ic,iz)
      END DO
      IF(RhoCPos>0.AND.CloudCover3<=One) WRITE(OutputUnit+2,*) CloudCover3,CloudBaseHeight3,CloudTopHeight3
    END SELECT
  END IF

END SUBROUTINE WriteProfile

SUBROUTINE WriteSoil(SoilOut,ic,SoilComp)

  REAL(RealKind), ALLOCATABLE :: SoilOut(:,:)
  CHARACTER(200)              :: SoilComp
  INTEGER                     :: ic,i,iz

  WRITE(OutputUnit,*) SoilComp
  DO iz=1,Domain%nrsoillayers
    WRITE(OutputUnit,*)  &
          0.5d0*(Domain%zSDepth(iz-1)+Domain%zSDepth(iz)) &
          ,SoilOut(1:ic,iz)
  END DO

END SUBROUTINE WriteSoil

SUBROUTINE StartEndPoints(ixStart,ixEnd,iyStart,iyEnd,ixProf,iyProf,ibProf)
  INTEGER :: ix,iy,iz
  INTEGER :: ixProf,iyProf
  INTEGER :: ixStart,iyStart
  INTEGER :: ixEnd,iyEnd
  INTEGER :: i,ibProf,ibProfLoc,k

!  CALL CheckValue

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (ProfOut) THEN
      IF (.NOT.(ProfAverage)) THEN
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            IF (xP(ix-1)<=xProf.AND. &
                xP(ix  )> xProf.AND. &
                yP(iy-1)<=yProf.AND. &
                yP(iy  )> yProf) THEN
              ixProf=ix
              iyProf=iy
              ibProf=ib
              ibProfLoc=ibLoc
            END IF
          END DO
        END DO
      ELSE
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            IF (xP(ix-1)<=xStart.AND. &
                xP(ix  )>=xStart) THEN
              ixStart=ix-1
            END IF
            IF (yP(iy-1)<=yStart.AND. &
                yP(iy  )>=yStart) THEN
              iyStart=iy-1
            END IF
            IF (xP(ix-1)<=xEnd.AND. &
                xP(ix  )>=xEnd) THEN
              ixEnd=ix
            END IF
            IF (yP(iy-1)<=yEnd.AND. &
                yP(iy  )>=yEnd) THEN
              iyEnd=iy
            END IF
          END DO
        END DO
      END IF
    END IF
  END DO

END SUBROUTINE StartEndPoints

SUBROUTINE CheckValue

  IF (.NOT.(ProfAverage)) THEN
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      IF (ib==LocGlob(1)) xStart=xP(ix0Out)
      IF (ib==LocGlob(nbLoc)) xEnd=xP(ix1Out)
      IF (ib==LocGlob(1)) yStart=yp(iy0Out)
      IF (ib==LocGlob(nbLoc)) yEnd=yP(iy1Out)
    END DO
  END IF

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (ProfOut) THEN
      IF (.NOT.(ProfAverage)) THEN
        IF (((xProf<xP(ix0Out)).AND.(ibLoc==1)).OR. &
            ((xProf>xP(ix1Out)).AND.(ibLoc==nbLoc))) THEN
          xProf=(xStart+xEnd)/2
          WRITE(*,*) 'xStart',xStart,'xEnd',xEnd
          WRITE(*,*) 'xProf not within output domain' &
                     ,' and set to center point: ',xProf
        END IF
        IF (((yProf<yP(iy0Out)).AND.(ibLoc==1)).OR. &
            ((yProf>yP(iy1Out)).AND.(ibLoc==nbLoc))) THEN
          yProf=(yStart+yEnd)/2
          WRITE(*,*) 'yStart',yStart,'yEnd',yEnd
          WRITE(*,*) 'yProf not within output domain' &
                     ,' and set to center point: ',yProf
        END IF
      ELSE
        IF (((xStart<xP(ix0Out)).AND.(ibLoc==1)).OR. &
            ((xStart>xP(ix1Out)).AND.(ibLoc==nbLoc))) THEN
          CALL Set(Floor(LocGlob(1)))
          xStart=xP(ix0Out)
          WRITE(*,*) 'xStart not within output domain' &
                     ,' and set to initial value: ',xStart
        END IF
        IF (((yStart<yP(iy0Out)).AND.(ibLoc==1)).OR. &
            ((yStart>yP(iy1Out)).AND.(ibLoc==nbLoc))) THEN
          CALL Set(Floor(LocGlob(1)))
          yStart=yP(iy0Out)
          WRITE(*,*) 'yStart not within output domain' &
                     ,' and set to initial value: ',yStart
        END IF
        IF (((xEnd<xP(ix0Out)).AND.(ibLoc==1)).OR. &
            ((xEnd>xP(ix1Out)).AND.(ibLoc==nbLoc))) THEN
          CALL Set(Floor(LocGlob(nbLoc)))
          xEnd=xP(ix1Out)
          WRITE(*,*) 'xEnd not within output domain' &
                     ,' and set to limiting value: ',xEnd
        END IF
        IF (((yEnd<yP(iy0Out)).AND.(ibLoc==1)).OR. &
            ((yEnd>yP(iy1Out)).AND.(ibLoc==nbLoc))) THEN
          CALL Set(Floor(LocGlob(nbLoc)))
          yEnd=yP(iy1Out)
          WRITE(*,*) 'yEnd not within output domain' &
                     ,' and set to limiting value: ',yEnd
        END IF
      END IF
    END IF
  END DO

END SUBROUTINE CheckValue


END MODULE Output_Mod
