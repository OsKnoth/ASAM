PROGRAM MainReadProf

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
  USE Tools_Mod
  USE SpecialOutput_Mod
  USE Soil_Mod
  USE ReadWeights_Mod
  USE Output_Mod
  USE Grid_Mod

  IMPLICIT NONE

  TYPE(VelocityFace_T)  , POINTER :: VelF1(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecT(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecTP(:)
  TYPE(VecVector4Cell_T), POINTER :: VecVelR(:)
  TYPE(PressureVelocity), POINTER :: bbb(:)

  CHARACTER(80)                   :: ProfileSTART='ProfileSTART'
  CHARACTER(80)                   :: ProfileEND='ProfileEND'
  CHARACTER(80)                   :: ProfileMid='ProfileMid'
  CHARACTER(80)                   :: InputFileName
  CHARACTER(80)                   :: ReadFileName
  CHARACTER(8)                    :: ScalarName
  CHARACTER(4)                    :: record_string

  REAL(RealKind) :: xM,yM,zM
  REAL(RealKind) :: xLow,xUpp
  REAL(RealKind) :: yLow,yUpp
  REAL(RealKind) :: zLow,zUpp
  REAL(RealKind) :: Temp,RhoLoc
  REAL(RealKind) :: LocTime,Time,MaxTime

  REAL(RealKind) :: phiM,phi1,z,phi2,z1,lam,zHeight
  REAL(RealKind) :: p,th,rhoM,rho1,rho2
  REAL(RealKind) :: u1,u2,uM,uMax,Cor1,Cor2,Cor,Cur1,Cur2,Cur
  REAL(RealKind) :: p1,p2,dpdphi
  REAL(RealKind) :: p1Min,p2Min
  REAL(RealKind) :: grad1,grad2
  REAL(RealKind) :: dphi1,dphi2
  REAL(RealKind) :: TotalStart(1:5),TotalCurrent(1:5)

  INTEGER        :: iInt,ix,iy,iz,ic
  INTEGER        :: i,j,k,iShift
  INTEGER        :: iFrac
  INTEGER        :: iFort
  INTEGER        :: Iter
  INTEGER        :: unr

  external unlimit_stack

! NETCDF Output
  INTEGER        :: Pos, Ind, Pos0, Pos1, NTMAX
  INTEGER        :: FileIndex,NumberInputFiles, DtInputFiles
  INTEGER        :: NETCDF_TimeStep, record
  INTEGER        :: uIndex,vIndex,wIndex
  INTEGER        :: rhoIndex,thIndex,tkeIndex,disIndex
  INTEGER        :: tkeHIndex,tkeVIndex,LenIndex
  INTEGER        :: qvIndex,qrIndex,qiIndex,qcIndex
  INTEGER        :: PreIndex,DiffIndex,DiffHIndex,DiffVIndex,DiffPotIndex,DiffMomIndex
  REAL(RealKind), ALLOCATABLE :: ProfileVar(:,:)
  REAL(RealKind), ALLOCATABLE :: ProfileVarLoc(:,:)
  REAL(RealKind), ALLOCATABLE :: ProfileArea(:,:)
  REAL(RealKind), ALLOCATABLE :: ProfileAreaLoc(:,:)
  LOGICAL :: InitProfile=.TRUE.
  LOGICAL :: ExistFile=.TRUE.
  INTEGER :: NumberProfiles
  INTEGER :: gmvStep
  CHARACTER(10) :: iName


  Time=0.0d0

  CALL start_MPI

! -- Lesen der Gitterdatei --
  CALL getarg(1,InputFileName)
  IF (MyId==0) THEN
    IF (InputFileName == ' ') THEN
!      STOP 'Dateiname fehlt'
    END IF
  END IF
  OutputNameSoil=TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))//'.Soil'
  CALL MPI_Bcast(InputFileName,80,MPI_CHARACTER,0,MPI_COMM_WORLD,MPIErr)
  
  WRITE(*,*)
  WRITE(*,*) 'MyId',MyId,InputFileName

! -- Compute Parameter
  CALL ComputeParameter 
! -- Einlesen des Gitters (Block-Struktur) --
  CALL Allocate(Floor)
  CALL inp_part(InputFileName)

  CALL Allocate(Floor)
  CALL ReadWeights(InputFileName)

  CALL InputModelTransport(InputFileName)
  CALL InputSystem(ChemieFile)
  CALL SetIndices
  CALL SetPosition
  CALL AllocateVec4Chemie(VecT,VectorComponentsT)

! Prandtl number
  ALLOCATE(PrandtlNumber(VectorComponentsT))
  PrandtlNumber=One
  IF (disPos>0) THEN
    PrandtlNumber(disPos)=Karm*Karm/((Cmy2-Cmy1)*SQRT(Cmy0))
  END IF
! FallVelocity
  ALLOCATE(FallVelocity(VectorComponentsT))
  FallVelocity=Zero
  IF (RhoCpos>0) THEN
    FallVelocity(RhoCpos)=0.0d0
  END IF
  ALLOCATE(ScaleMat2(VectorComponentsT))
  ScaleMat2=FallVelocity

  VectorComponentsM=VectorComponentsT
  CALL InputModelBC(InputFileName)
  Vect=Zero
! IF (Chemie) THEN
!   -- Eingabe Chemie --

    !CALL InputChemicalData(DataFile)
    WRITE(*,*) ' nach InputChemicalData in mainreadprof'
    CALL OutputChemie('Output')
! END IF

! -- Meteorologie --
  CALL ALLOCATE(VelF1)
  VelF1=Zero
  VecT=Zero

! -- Allokieren der Jacobi-Matrix f�r Transport
  CALL InputExample(InputFileName)
  IF (RhoPos>0) THEN
    CALL ALLOCATE(RhoProfG)
    CALL ALLOCATE(RhoCell)
    CALL ALLOCATE(PreCell)
  END IF
  IF ( thPos>0) THEN
    CALL ALLOCATE(thProfG)
    CALL ALLOCATE(DiffKoeff)
  END IF
  IF (Diffusion) THEN
    IF (TkeDis) THEN
      CALL Allocate     (DiffKoeff)
    ELSE IF (TkeSGS.OR.TkeLen.OR.NoTke) THEN
      CALL ALLOCATE     (DiffPotKoeff)
      CALL ALLOCATE     (DiffMomKoeff)
      CALL ALLOCATE     (LenKoeff)
    ELSE IF (TkeSmag) THEN
      CALL ALLOCATE     (DiffPotHKoeff)
      CALL ALLOCATE     (DiffPotVKoeff)
      CALL ALLOCATE     (DiffMomHKoeff)
      CALL ALLOCATE     (DiffMomVKoeff)
      CALL ALLOCATE     (LenKoeff)
    ELSE IF (DynSmag) THEN
      CALL ALLOCATE     (DiffPotHKoeff)
      CALL ALLOCATE     (DiffPotVKoeff)
      CALL ALLOCATE     (DiffMomHKoeff)
      CALL ALLOCATE     (DiffMomVKoeff)
      CALL ALLOCATE     (LenKoeff)
    ELSE IF (TkeHVLen) THEN
      CALL ALLOCATE     (DiffHKoeff)
      CALL ALLOCATE     (DiffVKoeff)
    ELSE
      CALL ALLOCATE     (DiffKoeff)
    END IF
  END IF

  IF (Damping) THEN
    CALL ALLOCATE     (DampKoeffCell)
  END IF
  IF (Radiation) THEN
    CALL ALLOCATE     (ShadowCell)
    CALL ALLOCATE     (RaddirCell)
    CALL ALLOCATE     (RaddifCell)
    CALL ALLOCATE     (RadinfCell)
  END IF

  CALL InputModelOutput(InputFileName)
  CALL InitOutputSpecial(InputFileName)

  ReadFileName='Input.out'
  WRITE(*,*) 'OutputFileName ',OutputFileName
  gmvStep=ProfgmvStep
  DO 
    WRITE(iName,'(I8)') gmvStep
    ReadFileName=TRIM(OutputFileName)//TRIM(ADJUSTL(iName))
    INQUIRE(FILE=TRIM(ReadFileName),EXIST=ExistFile)

    IF (.NOT.ExistFile) EXIT
    CALL InputMet(Time,VecT,VecT,ReadFileName)
    IF (Time>ProfTimeEnd) EXIT
    IF (Time>=ProfTimeStart.AND.Time<=ProfTimeEnd) THEN
      IF (InitProfile) THEN
        ALLOCATE(ProfileVarLoc(Domain%nz,3*SIZE(VectorOut(1)%Vec)+6))
        ALLOCATE(ProfileVar(Domain%nz,3*SIZE(VectorOut(1)%Vec)+6))
        ALLOCATE(ProfileAreaLoc(Domain%nz,SIZE(VectorOut(1)%Vec)+3))
        ALLOCATE(ProfileArea(Domain%nz,SIZE(VectorOut(1)%Vec)+3))
        ProfileVar=0.0d0
        ProfileArea=0.0d0
        InitProfile=.FALSE.
        NumberProfiles=0
      END IF  
! CALL OutputPointProfile(VecT)
      WRITE(*,*) 'VarianceProfile'
      CALL VarianceProfiles(VecT,ProfileVarLoc)
      CALL AreaProfiles(VecT,ProfileVarLoc)
      ProfileVar=ProfileVar+ProfileVarLoc
      ProfileArea=ProfileArea+ProfileAreaLoc
      NumberProfiles=NumberProfiles+1
    END IF  
    gmvStep=gmvStep+1
  END DO  
  WRITE(*,*) 'NumberProfiles',NumberProfiles
  ProfileVar=ProfileVar/NumberProfiles
  ProfileArea=ProfileArea/NumberProfiles
  CALL OutputAreaProfiles(ProfileVar)
  CALL OutputVarianceProfiles(ProfileVar)
! IF (ProfOut) THEN
!   CALL ReadOutputProfile(VecT,VecT,ProfileMid)
! END IF
  CALL MPI_Finalize(MPIErr)

END PROGRAM MainReadProf
