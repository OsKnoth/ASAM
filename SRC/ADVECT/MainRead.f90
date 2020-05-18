PROGRAM MainRead

  USE Init_Mod
  USE JacAccGrav_Mod
  USE Int_Mod
  USE ReadOutput_Mod
  USE Emission_Mod
  USE Koagulation_Mod
  USE ScalarVectorCellPar_Mod
  USE Microphysics_Mod
  USE Activity_Mod
  USE RhoProf_Mod
  USE SpecialOutput_Mod
  USE ReadWeights_Mod
  USE Grid_Mod

  USE Grid_Mod

  USE Output_Mod
  USE typesizes
  USE netcdf

  IMPLICIT NONE

  TYPE(VelocityFace_T)  , POINTER :: VelF1(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecT(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecT1(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecTP(:)
  TYPE(VecVector4Cell_T), POINTER :: VecVelR(:)
  TYPE(PressureVelocity), POINTER :: bbb(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecG(:)

  CHARACTER(80)                   :: ProfileSTART='ProfileSTART'
  CHARACTER(80)                   :: ProfileEND='ProfileEND'
  CHARACTER(80)                   :: InputFileName
  CHARACTER(80)                   :: ReadFileName
  CHARACTER(80)                   :: ReadFileName1
  CHARACTER(8)                    :: ScalarName
  CHARACTER(4)                    :: record_string

  CHARACTER (20) :: TIME_FREQ = "" 
  CHARACTER (20) :: GRID_TYPE = "" 
  CHARACTER (20) :: NAT_GRIDTYPE = "" 
  CHARACTER (20) :: DESCRIPT = "" 
  CHARACTER (20) :: NAME_EQU = "nonhydro" 
  CHARACTER (20) :: HORIZ_RESOL = "" 
  CHARACTER (20) :: TEST_NAME = "" 
  CHARACTER (20) :: LEVEL_NAME = "" 
  CHARACTER (20) :: MODEL_NAME = "asam" 
  CHARACTER (20) :: X_NAME = "x" 
  CHARACTER (20) :: Y_NAME = "y" 
  CHARACTER (20) :: X_LONGNAME = "x" 
  CHARACTER (20) :: Y_LONGNAME = "y" 
  CHARACTER (20) :: X_UNITS = "" 
  CHARACTER (20) :: Y_UNITS = "" 
  CHARACTER (20) :: X_AXIS = "" 
  CHARACTER (20) :: Y_AXIS = "" 

  REAL(RealKind) :: xM,yM,zM
  REAL(RealKind) :: xLow,xUpp
  REAL(RealKind) :: yLow,yUpp
  REAL(RealKind) :: zLow,zUpp
  REAL(RealKind) :: Temp,RhoLoc
  REAL(RealKind) :: LocTime,Time,MaxTime

  REAL(RealKind) :: phiM,phi1,z,phi2,lam,zHeight
  REAL(RealKind) :: u1,u2,uM,uMax,Cor1,Cor2,Cor,Cur1,Cur2,Cur
  REAL(RealKind) :: p1,p2,dpdphi
  REAL(RealKind) :: p1Min,p2Min
  REAL(RealKind) :: grad1,grad2
  REAL(RealKind) :: dphi1,dphi2
  REAL(RealKind) :: TotalStart(1:5),TotalCurrent(1:5)

  INTEGER        :: iInt,ix,iy,iz
  INTEGER        :: it
  INTEGER        :: i,j,k,iShift
  INTEGER        :: iFrac
  INTEGER        :: iFort
  INTEGER        :: Iter
  INTEGER        :: unr
  LOGICAL        :: FileExist

  external unlimit_stack

! NETCDF Output
  INTEGER        :: Pos, Ind, Pos0, Pos1, NTMAX
  INTEGER        :: FileIndex,NumberInputFiles, DtInputFiles
  INTEGER        :: record
  REAL(RealKind) :: NETCDF_TimeStep
  INTEGER        :: gmvStartNum
  INTEGER        :: uIndex,vIndex,wIndex
  INTEGER        :: rhoIndex,thIndex,enIndex,tkeIndex,disIndex
  INTEGER        :: tkeHIndex,tkeVIndex,LenIndex
  INTEGER        :: RhoVIndex,RhoRIndex,RhoIIndex,RhoCIndex
  INTEGER        :: Tracer1Index,Dummy1Index,Dummy2Index,Dummy3Index,Dummy4Index
  INTEGER        :: PreIndex,DiffIndex,DiffHIndex,DiffVIndex,DiffPotIndex,DiffMomIndex,DiffMomHIndex
  INTEGER        :: DiffMomVIndex,AbsTIndex
  INTEGER, ALLOCATABLE :: GasIndex(:),AeroIndex(:,:)
  INTEGER, ALLOCATABLE :: SpecialIndex(:)
  REAL(RealKind), DIMENSION (:,:,:,:,:), ALLOCATABLE :: OutNetcdf

  Time=0.0d0

  CALL start_MPI

! -- Lesen der Gitterdatei --
  CALL getarg(1,InputFileName)
  WRITE(*,*)
  IF (MyId==0) THEN
    IF (InputFileName == ' ') THEN
!      STOP 'Dateiname fehlt'
    END IF
  END IF
  OutputNameSoil=TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))//'.Soil'
  CALL MPI_Bcast(InputFileName,80,MPI_CHARACTER,0,MPI_COMM_WORLD,MPIErr)
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
! CALL Allocate(JacMet)
! JacMet=Zero
  CALL AllocateVec4Chemie(VecT,VectorComponentsT)
  CALL AllocateVec4Chemie(VecT1,VectorComponentsT)

! Prandtl number
  ALLOCATE(PrandtlNumber(VectorComponentsT))
  PrandtlNumber=One
  IF (disPos>0) THEN
    PrandtlNumber(disPos)=Karm*Karm/((Cmy2-Cmy1)*SQRT(Cmy0))
  END IF
! FallVelocity
  ALLOCATE(FallVelocity(VectorComponentsT))
  FallVelocity=Zero
  IF (RhoCPos>0) THEN
    FallVelocity(RhoCPos)=0.0d0
  END IF
  ALLOCATE(ScaleMat2(VectorComponentsT))
  ScaleMat2=FallVelocity

  VectorComponentsM=VectorComponentsT
  CALL InputModelBC(InputFileName)
  Vect=Zero
! IF (Chemie) THEN
!   -- Eingabe Chemie --
    WRITE(*,*) ' vor inputchemicaldata ' 
    CALL InputChemicalData(DataFile)
!    CALL OutputChemie('Output')
!   CALL OutputChemie(ChemOutFile)
! END IF
 ALLOCATE(GasIndex(NumGasOut))
 ALLOCATE(AeroIndex(NumAeroOut,nFrac))

! -- Meteorologie --
  CALL ALLOCATE(VelF1)
  VelF1=Zero
  VecT=Zero
  VecT1=Zero

! -- Allokieren der Jacobi-Matrix f�r Transport
! CALL Allocate(JacMet)

! CALL InputExample(InputFileName)
  CALL SetMetVariables
  IF (RhoPos>0) THEN
    CALL ALLOCATE(RhoProfG)
    CALL ALLOCATE(RhoCell)
    CALL ALLOCATE(PreCell)
  END IF
  IF ( thPos>0) THEN
    CALL ALLOCATE(thProfG)
    CALL ALLOCATE(DiffKoeff)
    CALL ALLOCATE(TAbsCell,1,1)
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
! CALL InitOutputSpecial(InputFileName)
! CALL SpecialOutput(VecT,VecT,VelF1,Time)
  ALLOCATE(SpecialIndex(LenOutSpecial))

!  ReadFileName=TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))//'.out'
  ReadFileName=OutputFileName
  WRITE(*,*) 'Output Filename :   ',ReadFileName
!  READ(*,*) ReadFileName
! WRITE(*,*) 'Input ReadFileName1'
! READ(*,*) ReadFileName1
! CALL InputMet(Time,VecT,VecT,ReadFileName)

! ------------------------------------------   
! - Number of Files -
! ------------------------------------------   
  NumberInputFiles = INT((EndTime-StartTime)/OutputTimeStep)

! ------------------------------------------
! - User specific Number of Inputfiles - 
!
! - You can comment that part if You don't
!   like to change the OutputTimeStep!
! ------------------------------------------
  WRITE(*,*)'########################################################'
  WRITE(*,*)'User specific Number of Timesteps used for NETCDF-Output'
  WRITE(*,*)' '
  WRITE(*,*)'  StartTime      = ',StartTime
  WRITE(*,*)'  EndTime        = ',EndTime
  WRITE(*,*)'  OutputTimeStep = ',OutputTimeStep
  WRITE(*,*)' '
  WRITE(*,*)'NETCDF_TimeStep - A multiple integer of OutputTimeStep'
  WRITE(*,*)'NETCDF_TimeStep  >=  OutputTimeStep'
  WRITE(*,*)'########################################################'
  WRITE(*,*)' '
  WRITE(*,*)'GMV Output Start Number (e.g., 1000) =  ... [integer]'
  READ(*,*) gmvStartNum
  WRITE(*,*)' '
  WRITE(*,*)'NETCDF_TimeStep =  ... [s]'
  READ(*,*) NETCDF_TimeStep
  WRITE(*,*)'gmvStartNum = '    ,gmvStartNum 
  WRITE(*,*)'NETCDF_TimeStep = ',NETCDF_TimeStep,'[s]' 
 
  IF (NETCDF_TimeStep < OutputTimeStep) THEN
    WRITE(*,*)'########################################################'
    WRITE(*,*)' '
    WRITE(*,*)'Failure: NETCDF_TimeStep < OutputTimeStep'
    WRITE(*,*)'Program -STOP-'
    WRITE(*,*)' '
    WRITE(*,*)'########################################################'
    STOP
  END IF 
  WRITE(*,*) 'Filename gelesen 1'
! WRITE(*,*) 'Input ReadFileName1'
  DtInputFiles = INT(NETCDF_TimeStep/OutputTimeStep)
  NTMAX        = INT((EndTime-StartTime)/NETCDF_TimeStep)

! ------------------------------------------
! - Get Number of Variables for OutNetcdf -
! ------------------------------------------
  Pos0 = 0
  Pos1 = 0
  IF (uPosL   > 0) THEN
    Pos1      = Pos1 + 1
    uIndex    = Pos1
  END IF
  IF (vPosL   > 0) THEN
    Pos1      = Pos1 + 1
    vIndex    = Pos1
  END IF
  IF (wPosL   > 0) THEN
    Pos1      = Pos1 + 1
    wIndex    = Pos1
  END IF
  IF (thPos   > 0) THEN
    Pos1      = Pos1 + 1
    thIndex   = Pos1
  END IF
  IF (EnPos   > 0) THEN
    Pos1      = Pos1 + 1
    EnIndex   = Pos1
  END IF
  IF (rhoPos  > 0) THEN
    Pos1      = Pos1 + 1
    rhoIndex  = Pos1
  END IF
  IF (tkePos  > 0) THEN
    Pos1      = Pos1 + 1
    tkeIndex  = Pos1
  END IF
  IF (disPos  > 0) THEN
    Pos1      = Pos1 + 1
    disIndex  = Pos1
  END IF
  IF (tkeHPos > 0) THEN
    Pos1      = Pos1 + 1
    tkeHIndex = Pos1
  END IF
  IF (tkeVPos > 0) THEN
    Pos1      = Pos1 + 1
    tkeVIndex = Pos1
  END IF
  IF (LenPos  > 0) THEN
    Pos1      = Pos1 + 1
    LenIndex  = Pos1
  END IF
  IF (RhoVPos   > 0) THEN
    Pos1      = Pos1 + 1
    RhoVIndex   = Pos1
  END IF
  IF (RhoCPos   > 0) THEN
    Pos1      = Pos1 + 1
    RhoCIndex   = Pos1
  END IF
  IF (RhoIPos   > 0) THEN
    Pos1      = Pos1 + 1
    RhoIIndex   = Pos1
  END IF
  IF (RhoRPos   > 0) THEN
    Pos1      = Pos1 + 1
    RhoRIndex   = Pos1
  END IF
  IF (Position('Dummy1')   > 0) THEN
    Pos1      = Pos1 + 1
    Dummy1Index   = Pos1
  END IF
  IF (Position('Tracer1')   > 0) THEN
    Pos1      = Pos1 + 1
    Tracer1Index   = Pos1
  END IF
  IF (Position('Dummy2')   > 0) THEN
    Pos1      = Pos1 + 1
    Dummy2Index   = Pos1
  END IF
  IF (Position('Dummy3')   > 0) THEN
    Pos1      = Pos1 + 1
    Dummy3Index   = Pos1
  END IF
  IF (Position('Dummy4')   > 0) THEN
    Pos1      = Pos1 + 1
    Dummy4Index   = Pos1
  END IF
! IF (LenOutSpecial>0) THEN
!   DO i=1,LenOutSpecial
!     Pos1 = Pos1 + 1
!     SpecialIndex(i) = Pos1
!   END DO
! END IF
  IF (ChemieOut) THEN
    DO i=1,NumGasOut
      Pos1 = Pos1 + 1
      GasIndex(i) = Pos1
    END DO
    DO i=1,NumAeroOut
      DO it=1,nFrac
        Pos1 = Pos1 + 1
        AeroIndex(i,it) = Pos1
      END DO
    END DO
  END IF

  IF (PreOut     ) THEN
    Pos1          = Pos1 + 1
    PreIndex      = Pos1
  END IF
  IF (DiffOut    ) THEN
    Pos1          = Pos1 + 1
    DiffIndex     = Pos1
  END IF
  IF (DiffMomHOut    ) THEN !marcelk
    Pos1          = Pos1 + 1
    DiffMomHIndex     = Pos1
  END IF
  IF (DiffMomVOut    ) THEN !marcelk
    Pos1          = Pos1 + 1
    DiffMomVIndex     = Pos1
  END IF
  IF (DiffHOut   ) THEN
    Pos1          = Pos1 + 1
    DiffHIndex    = Pos1
  END IF
  IF (DiffVOut   ) THEN
    Pos1          = Pos1 + 1
    DiffVIndex    = Pos1
  END IF
  IF (DiffPotOut ) THEN
    Pos1          = Pos1 + 1
    DiffPotIndex  = Pos1
  END IF
  IF (DiffMomOut ) THEN
    Pos1          = Pos1 + 1
    DiffMomIndex  = Pos1
  END IF

! ------------------------------------------
! - Allocate OutNetcdf -
! ------------------------------------------

  CALL SetDomain(Domain)

  ALLOCATE(OutNetcdf(Pos0+1:Pos1, 1:NTMAX+1, ix0+1:ix1, &
                                             iy0+1:iy1, &
                                             iz0+1:iz1))

! ------------------------------------------
! - Set OutNetcdf by Loop over Inputfiles -
! ------------------------------------------
  FileIndex = 0
  DO i = 1, NumberInputFiles+1, DtInputFiles
    record = INT(gmvStartNum - 1 + i)
    FileIndex = FileIndex + 1
    write(*,*) 'Processing: ', record, FileIndex
! ------------------------------------------
! - Read Files -
! ------------------------------------------
    WRITE(record_string,'(I4)') record
    ReadFileName = TRIM(ReadFileName(1:INDEX(ReadFileName,'.out')-1)) &
                                    //'.out'//ADJUSTL(record_string)
    WRITE (*,*) 'Filename: ',ReadFileName
    INQUIRE(FILE=ReadFileName,EXIST=FileExist)
    IF (.NOT.FileExist) THEN
      NTMAX=i-2
      WRITE(*,*) 'EXIT',i
      EXIT
    END IF  
    WRITE (*,*) 'CALL InputMet'
    CALL InputMet(Time,VecT,VecT,ReadFileName)
    IF (ProfOut) THEN
      WRITE (*,*) 'CALL OutputProfile'
      CALL OutputProfile(VecT,VecT,ProfileEND)
    END IF

! ------------------------------------------
! - Set OutNetcdf -
! ------------------------------------------
    WRITE (*,*) 'Set OutNetcdf'
    Ind = 0
    IF (uPosL   > 0) THEN
      Pos = uPosL
      Ind = uIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (vPosL   > 0) THEN
      Pos = vPosL
      Ind = vIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (wPosL   > 0) THEN
      Pos = wPosL
      Ind = wIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (thPos   > 0) THEN
      Pos = thPos
      Ind = thIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (EnPos   > 0) THEN
      Pos = EnPos
      Ind = EnIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (rhoPos  > 0) THEN
      Pos = rhoPos
      Ind = rhoIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
!    TKE-Closure k-eps
    IF (tkePos  > 0) THEN
      Pos = tkePos
      Ind = tkeIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (disPos  > 0) THEN
      Pos = disPos
      Ind = disIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
!    TKE-Closure following Herzog (ATHAM)
    IF (tkeHPos > 0) THEN
      Pos = tkeHPos
      Ind = tkeHIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (tkeVPos > 0) THEN
      Pos = tkeVPos
      Ind = tkeVIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (LenPos  > 0) THEN
      Pos = LenPos
      Ind = LenIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
!    Feuchte
    IF (RhoVPos   > 0) THEN
      Pos = RhoVPos
      Ind = RhoVIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (RhoCPos   > 0) THEN
      Pos = RhoCPos
      Ind = RhoCIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (RhoIPos   > 0) THEN
      Pos = RhoIPos
      Ind = RhoIIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (RhoRPos   > 0) THEN
      Pos = RhoRPos
      Ind = RhoRIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (Position('Dummy1') > 0) THEN
      Pos = Position('Dummy1')
      Ind = Dummy1Index
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
IF (Position('Tracer1') > 0) THEN
      Pos = Position('Tracer1')
      Ind = Tracer1Index
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (Position('Dummy2') > 0) THEN
      Pos = Position('Dummy2')
      Ind = Dummy2Index
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (Position('Dummy3') > 0) THEN
      Pos = Position('Dummy3')
      Ind = Dummy3Index
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (Position('Dummy4') > 0) THEN
      Pos = Position('Dummy4')
      Ind = Dummy4Index
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF

!   Additional output
!!! not working due to SpecialOutput !!!
!!! DISABLED !!!
!   WRITE (*,*) 'Further diagnostic variables'
!   CALL SetOutNetcdf(OutNetcdf,FileIndex)

  END DO

! ------------------------------------------
! - Output Netcdf -
! ------------------------------------------
  WRITE (*,*) 'CALL SetDomain'
  CALL SetDomain(Domain)

  WRITE (*,*) 'CALL Netcdf_Output'
  CALL Netcdf_Output(OutNetcdf,InputFileName)

  DEALLOCATE(OutNetcdf)

CONTAINS

SUBROUTINE SetOutNetcdfVec(VecC,OutC,PosC,IndC,FileIndex)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: i,ix,iy,iz
  INTEGER :: PosC, IndC, time
  INTEGER :: FileIndex
  REAL(RealKind) :: OutC(:,:,:,:,:)
  
  WRITE (*,*) 'Routine SetOutNetcdfVec: PosC,IndC,FileIndex ',PosC,IndC,FileIndex

  DO ibLoc=1,nbLoc
     ib=LocGlob(ibLoc)
     CALL DomainSet(ib)
     DO iz=igz0+1,igz1
       DO iy=igy0+1,igy1
         DO ix=igx0+1,igx1
            OutC(IndC,FileIndex,ix,iy,iz)=VecC(ibLoc)%Vec(PosC)%c(ix,iy,iz,1)
        !  IF (PosC==5.AND.OutC(IndC,FileIndex,ix,iy,iz).gt.1e-08) write(*,*) OutC(IndC,FileIndex,ix,iy,iz)
         END DO
       END DO
     END DO
  END DO

END SUBROUTINE SetOutNetcdfVec

SUBROUTINE SetOutNetcdf(OutC,FileIndex)

  IMPLICIT NONE

  INTEGER :: i,ix,iy,iz
  INTEGER :: it
  INTEGER :: IndC
  INTEGER :: FileIndex
  REAL(RealKind) :: OutC(:,:,:,:,:)
    
  DO ibLoc=1,nbLoc
     ib=LocGlob(ibLoc)
     CALL DomainSet(ib)
     DO iz=igz0+1,igz1
       DO iy=igy0+1,igy1
         DO ix=igx0+1,igx1
         ! IF (PreOut     ) THEN
         !   IndC = PreIndex
         !   OutC(IndC,FileIndex,ix,iy,iz)=PreCell(ibLoc)%c(ix,iy,iz,1)
         ! END IF
           IF (DiffOut    ) THEN
             IndC = DiffIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffKoeff(ibLoc)%c(ix,iy,iz,1)
           END IF
           IF (DiffMomHOut    ) THEN !marcelk
             IndC = DiffMomHIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffMomHKoeff(ibLoc)%c(ix,iy,iz,1)
           END IF
           IF (DiffMomVOut    ) THEN !marcelk
             IndC = DiffMomVIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffMomVKoeff(ibLoc)%c(ix,iy,iz,1)
           END IF
           IF (DiffHOut   ) THEN
             IndC = DiffHIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffHKoeff(ibLoc)%c(ix,iy,iz,1)
           END IF
           IF (DiffVOut   ) THEN
             IndC = DiffVIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffVKoeff(ibLoc)%c(ix,iy,iz,1)
           END IF
           IF (DiffPotOut ) THEN
             IndC = DiffPotIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffPotKoeff(ibLoc)%c(ix,iy,iz,1)
           END IF
           IF (DiffMomOut ) THEN
             IndC = DiffMomIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffMomKoeff(ibLoc)%c(ix,iy,iz,1)
           END IF
       !   IF (LenOutSpecial>0) THEN
       !     DO j=1,LenOutSpecial
       !       Pos = j
       !       IndC = SpecialIndex(j)
       !       OutC(IndC,FileIndex,ix,iy,iz)=OutSpecial(ibLoc)%Vec(Pos)%c(ix,iy,iz,1)
       !     END DO
       !   END IF
           IF (ChemieOut) THEN
             DO j=1,NumGasOut
               Pos = Position(GasSpecies(j))
               IndC = GasIndex(j)
               OutC(IndC,FileIndex,ix,iy,iz)=VecT(ibLoc)%Vec(Pos)%c(ix,iy,iz,1)
             END DO
             DO j=1,NumAeroOut
               DO it=1,nFrac
                 Pos = Position(AeroSpecies(j,it))
                 IndC = AeroIndex(j,it)
                 OutC(IndC,FileIndex,ix,iy,iz)=VecT(ibLoc)%Vec(Pos)%c(ix,iy,iz,it)
               END DO
             END DO
           END IF

         END DO
       END DO
     END DO
  END DO

END SUBROUTINE SetOutNetcdf

SUBROUTINE Netcdf_Output(OutNetcdf,InputFileName)

  ! This is the name of the data file we will create.
  CHARACTER (80) :: FILE_NAME 
  CHARACTER (80) :: InputFileName
  INTEGER :: ncid

  INTEGER,             PARAMETER :: NDIMS  = 4
!  CHARACTER (LEN = *), PARAMETER :: X_NAME = "x"
!  CHARACTER (LEN = *), PARAMETER :: Y_NAME = "y"
!  CHARACTER (LEN = *), PARAMETER :: Z_NAME = "z"
  CHARACTER (LEN = *), PARAMETER :: Z_NAME = "lev" 
  CHARACTER (LEN = *), PARAMETER :: ZINT_NAME         = "ilev"
  CHARACTER (LEN = *), PARAMETER :: Z_LONGNAME        = "Height at midpoints" 
  CHARACTER (LEN = *), PARAMETER :: ZINT_LONGNAME     = "Height at interfaces" 
  CHARACTER (LEN = *), PARAMETER :: REC_NAME          = "time"
  CHARACTER (LEN = *), PARAMETER :: TIMESTEP_NAME     = "mdt"
  CHARACTER (LEN = *), PARAMETER :: REC_LONGNAME      = "time"
  CHARACTER (LEN = *), PARAMETER :: TIMESTEP_LONGNAME = "timestep"
  CHARACTER (LEN = *), PARAMETER :: PRES_LONGNAME     = "reference pressure"
  INTEGER :: x_dimid, y_dimid, z_dimid, zint_dimid, rec_dimid
  INTEGER :: rec

  ! The start and count arrays will tell the netCDF library where to
  ! write our data.
  INTEGER :: start(NDIMS), count(NDIMS)

  ! These program variables hold the cartesian grid
  INTEGER :: x_varid, y_varid, z_varid, zint_varid
  INTEGER :: rec_varid,dt_varid,time_varid
  INTEGER :: glob_varid

  ! We will create two netCDF variables, one each for temperature and
  ! pressure fields.
  CHARACTER (LEN = *), PARAMETER :: U_NAME        ="U"
  CHARACTER (LEN = *), PARAMETER :: V_NAME        ="V"
  CHARACTER (LEN = *), PARAMETER :: W_NAME        ="W"
  CHARACTER (LEN = *), PARAMETER :: TH_NAME       ="THETA"
  CHARACTER (LEN = *), PARAMETER :: EN_NAME       ="TOTE"
  CHARACTER (LEN = *), PARAMETER :: RHO_NAME      ="RHO"
  CHARACTER (LEN = *), PARAMETER :: TKE_NAME      ="TKE"
  CHARACTER (LEN = *), PARAMETER :: DIS_NAME      ="DIS"
  CHARACTER (LEN = *), PARAMETER :: TKEH_NAME     ="TKEH"
  CHARACTER (LEN = *), PARAMETER :: TKEV_NAME     ="TKEV"
  CHARACTER (LEN = *), PARAMETER :: LEN_NAME      ="LEN"
  CHARACTER (LEN = *), PARAMETER :: RhoV_NAME       ="RhoV"
  CHARACTER (LEN = *), PARAMETER :: RhoC_NAME       ="RhoC"
  CHARACTER (LEN = *), PARAMETER :: RhoI_NAME       ="RhoI"
  CHARACTER (LEN = *), PARAMETER :: RhoR_NAME       ="RhoR"
  CHARACTER (LEN = *), PARAMETER :: PRE_NAME      ="PRE"
  CHARACTER (LEN = *), PARAMETER :: PRES_NAME      ="P0"
  CHARACTER (LEN = *), PARAMETER :: DIFF_NAME     ="DIFF"
  CHARACTER (LEN = *), PARAMETER :: DIFFH_NAME    ="DIFFH"
  CHARACTER (LEN = *), PARAMETER :: DIFFV_NAME    ="DIFFV"
  CHARACTER (LEN = *), PARAMETER :: DIFFPOT_NAME  ="DIFFPOT"
  CHARACTER (LEN = *), PARAMETER :: DIFFMOM_NAME  ="DIFFMOM"
  CHARACTER (LEN = *), PARAMETER :: DUMMY1_NAME   ="DUMMY1"
  CHARACTER (LEN = *), PARAMETER :: TRACER1_NAME   ="TRACER1"
  CHARACTER (LEN = *), PARAMETER :: DUMMY2_NAME   ="DUMMY2"
  CHARACTER (LEN = *), PARAMETER :: DUMMY3_NAME   ="DUMMY3"
  CHARACTER (LEN = *), PARAMETER :: DUMMY4_NAME   ="DUMMY4"
  INTEGER :: u_varid, v_varid, w_varid
  INTEGER :: th_varid, en_varid, rho_varid
  INTEGER :: tke_varid, dis_varid
  INTEGER :: tkeh_varid, tkev_varid, len_varid
  INTEGER :: RhoV_varid, RhoC_varid, RhoI_varid, RhoR_varid
  INTEGER :: pre_varid, pres_varid, diff_varid, diffh_varid, diffv_varid
  INTEGER :: diffpot_varid, diffmom_varid
  INTEGER :: tracer1_varid,dummy1_varid,dummy2_varid,dummy3_varid,dummy4_varid
  INTEGER,ALLOCATABLE :: gas_varid(:),aero_varid(:,:)
  INTEGER,ALLOCATABLE :: special_varid(:)
  INTEGER :: dimids(NDIMS)
  INTEGER :: pres_dimid

  ! We recommend that each variable carry a "units" attribute.
  CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
!  CHARACTER (LEN = *), PARAMETER :: X_UNITS    = "m"
!  CHARACTER (LEN = *), PARAMETER :: Y_UNITS    = "m"
  CHARACTER (LEN = *), PARAMETER :: Z_UNITS    = "m"
  CHARACTER (LEN = *), PARAMETER :: ZINT_UNITS = "m"
  CHARACTER (LEN = *), PARAMETER :: REC_UNITS  = "seconds since 2000-01-01 00:00:00"
  CHARACTER (LEN = *), PARAMETER :: TIMESTEP_UNITS  = "s"
  CHARACTER (LEN = *), PARAMETER :: U_UNITS    = "m/s"
  CHARACTER (LEN = *), PARAMETER :: V_UNITS    = "m/s"
  CHARACTER (LEN = *), PARAMETER :: W_UNITS    = "m/s"
  CHARACTER (LEN = *), PARAMETER :: TH_UNITS   = "K"
  CHARACTER (LEN = *), PARAMETER :: EN_UNITS   = "K"
  CHARACTER (LEN = *), PARAMETER :: RHO_UNITS  = "kg/m**3"
  CHARACTER (LEN = *), PARAMETER :: TKE_UNITS  = "m**2/s**2"
  CHARACTER (LEN = *), PARAMETER :: DIS_UNITS  = "m**2/s**2"
  CHARACTER (LEN = *), PARAMETER :: TKEH_UNITS = "m**2/s**2"
  CHARACTER (LEN = *), PARAMETER :: TKEV_UNITS = "m**2/s**2"
  CHARACTER (LEN = *), PARAMETER :: LEN_UNITS  = "m"
  CHARACTER (LEN = *), PARAMETER :: RhoV_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: RhoC_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: RhoI_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: RhoR_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: PRE_UNITS       = " "
  CHARACTER (LEN = *), PARAMETER :: PRES_UNITS       = "Pa"
  CHARACTER (LEN = *), PARAMETER :: DIFF_UNITS      = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DIFFH_UNITS     = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DIFFV_UNITS     = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DIFFPOT_UNITS   = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DIFFMOM_UNITS   = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DUMMY1_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: TRACER1_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: DUMMY2_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: DUMMY3_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: DUMMY4_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: GAS_UNITS    = "mug/m3" !"g/mol "
  CHARACTER (LEN = *), PARAMETER :: AERO_UNITS   = "mol/l" !"g/mol "
  CHARACTER (20)                 :: LEVELnr   
  CHARACTER (20)                 :: TFREQ   

  ! -----------------------------------
  ! - ASAM -
  ! -----------------------------------
  REAL(RealKind) :: OutNetcdf(:,:,:,:,:)

  !---  local variables
  INTEGER :: ix, iy, iz, it, i,itime
  INTEGER :: iScalar,NumOut
  INTEGER :: NXX, NYY, NZZ, NTT
!  INTEGER, ALLOCATABLE :: xLoc(:),yLoc(:),zLoc(:)
  REAL(RealKind), ALLOCATABLE :: xLoc(:),yLoc(:),zLoc(:)
  REAL(RealKind), ALLOCATABLE :: z_int(:),timeLoc(:)

  ALLOCATE(aero_varid(NumAeroOut,nFrac))
  ALLOCATE(gas_varid(NumGasOut))
  ALLOCATE(special_varid(LenOutSpecial))

  ! Create the filename of the Netcdf-Outputfile
  FILE_NAME = TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))//'.nc' !'.cdf'

  ! Program variables to hold the data we will write out. 
  NXX = ix1Out - ix0Out
  NYY = iy1Out - iy0Out
  NZZ = iz1Out - iz0Out
! NTT = INT((EndTime-StartTime)/OutputTimeStep)
  NTT=NTMAX 
  WRITE(*,*) 'Number of time steps ',NTT+1
  WRITE(*,*) 'StartTime ',StartTime
  WRITE(*,*) 'EndTime ',EndTime

  WRITE(LEVELnr,'(I20)') NZZ
  WRITE(TFREQ,'(F10.2)') OutputTimeStep

  IF (Sphere) THEN
    X_NAME = "lon"
    Y_NAME = "lat"
    X_LONGNAME = "longitude"
    Y_LONGNAME = "latitude"
    X_UNITS = "degrees_east"
    Y_UNITS = "degrees_north"
    GRID_TYPE = "latlon"
    NAT_GRIDTYPE = "latlon"
    X_AXIS = "X"
    Y_AXIS = "Y"
  ELSE
    X_NAME = "xc"
    Y_NAME = "yc"
    X_LONGNAME = "x-coordinate in Cartesian system"
    Y_LONGNAME = "y-coordinate in Cartesian system"
    X_UNITS = "m"
    Y_UNITS = "m"
    GRID_TYPE = "cart"
    NAT_GRIDTYPE = "cart"
  END IF
  
  ALLOCATE(xLoc(1:NXX))
  ALLOCATE(yLoc(1:NYY)) 
  ALLOCATE(zLoc(1:NZZ))
  ALLOCATE(z_int(1:NZZ+1))
  ALLOCATE(timeLoc(0:NTT))

  DO ibLoc=1,nbLoc
   !WRITE(*,*) 'ibLoc',ibLoc,xP(0),xP(1),dx(1),Floor(ib)%dx(1)
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))

!   DO ix = 1, NXX
!     IF (Sphere) THEN
!       xLoc(ix)=Domain%x0+(ix-1)*(Domain%x1-Domain%x0)/Domain%nx
!       xLoc(ix)=xLoc(ix)*180.0d0/Pi
!     ELSE
!       xLoc(ix)=xP(ix-1)+0.5d0*dx(ix)
!     END IF
!   END DO
!   DO iy = 1, NYY
!     IF (Sphere) THEN
!       yLoc(iy)=Domain%y0+(iy-1)*(Domain%y1-Domain%y0)/Domain%ny
!       yLoc(iy)=yLoc(iy)*180.0d0/Pi
!     ELSE
!      yLoc(iy)=yP(iy-1)+0.5d0*dy(iy)
!    END IF
!   END DO

    DO itime = 0, NTT
      timeLoc(itime) = StartTime+OutputTimeStep*itime ! /3600.0d0
    END DO

!   CALL Set(Floor(1))
    DO ix=ix0+1,ix1
      xLoc(ix) = xP(ix-1)+0.5d0*dx(ix)
    !  WRITE(*,*) 'x',ix,xLoc(ix)
    END DO
    DO iy=iy0+1,iy1
      yLoc(iy) = yP(iy-1)+0.5d0*dy(iy)
    !  WRITE(*,*) 'y',iy,yLoc(iy)
    END DO
    DO iz=iz0+1,iz1
      zLoc(iz) = zP(iz-1)+0.5d0*dz(iz)
      z_int(iz) = zP(iz-1) 
    !  WRITE(*,*) 'z',iz,zLoc(iz)
    END DO
    z_int(iz1+1) = zP(iz1) 

  END DO

  TEST_NAME=TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))
  LEVEL_NAME=TRIM("L")//ADJUSTL(LEVELnr)
  TIME_FREQ=ADJUSTL(TRIM(TFREQ))//TRIM("s")

  IF (INDEX(InputFileName,'DCMIP')>0) THEN
    CALL InputDCMIP(InputFileName)
    IF (LEN(TRIM(DESCRIPT))>0) THEN
      FILE_NAME = TRIM(MODEL_NAME)//'.'//TRIM(TEST_NAME)//'.'//TRIM(HORIZ_RESOL)//'.'&
                //TRIM(LEVEL_NAME)//'.'//TRIM(GRID_TYPE)//'.'//TRIM(NAME_EQU)//'.' &
                //TRIM(DESCRIPT)//'.nc' 
    ELSE
      FILE_NAME = TRIM(MODEL_NAME)//'.'//TRIM(TEST_NAME)//'.'//TRIM(HORIZ_RESOL)//'.'&
                //TRIM(LEVEL_NAME)//'.'//TRIM(GRID_TYPE)//'.'//TRIM(NAME_EQU)//'.nc' 
    END IF
    DO iy = 1, NYY
      yLoc(iy)=Domain%y0+(iy-1)*(Domain%y1-Domain%y0)/Domain%ny
      yLoc(iy)=yLoc(iy)*180.0d0/Pi
    END DO
    DO ix = 1, NXX
      xLoc(ix)=Domain%x0+(ix-1)*(Domain%x1-Domain%x0)/Domain%nx
      xLoc(ix)=xLoc(ix)*180.0d0/Pi
    END DO
  END IF


  ! Cr
  CALL check( nf90_create(FILE_NAME, nf90_clobber, ncid) )

  ! Define the dimensions. The record dimension is defined to have
  ! unlimited length - it can grow as needed. In this example it is
  ! the time dimension.
  CALL check( nf90_def_dim(ncid, TRIM(X_NAME), NXX, x_dimid) )
  CALL check( nf90_def_dim(ncid, TRIM(Y_NAME), NYY, y_dimid) )
  CALL check( nf90_def_dim(ncid, Z_NAME, NZZ, z_dimid) )
  CALL check( nf90_def_dim(ncid, ZINT_NAME, NZZ+1, zint_dimid) )
  CALL check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

  ! Define the coordinate variables. We will only define coordinate
  ! variables for lat and lon.  Ordinarily we would need to provide
  ! an array of dimension IDs for each variable's dimensions, but
  ! since coordinate variables only have one dimension, we can
  ! simply provide the address of that dimension ID (lat_dimid) and
  ! similarly for (lon_dimid).
  CALL check( nf90_def_var(ncid, PRES_NAME, NF90_DOUBLE, varid=pres_varid) )
  CALL check( nf90_def_var(ncid, TRIM(X_NAME), NF90_DOUBLE, x_dimid, x_varid) )
  CALL check( nf90_def_var(ncid, TRIM(Y_NAME), NF90_DOUBLE, y_dimid, y_varid) )
  CALL check( nf90_def_var(ncid, Z_NAME, NF90_DOUBLE, z_dimid, z_varid) )
  CALL check( nf90_def_var(ncid, ZINT_NAME, NF90_DOUBLE, zint_dimid, zint_varid) )
  CALL check( nf90_def_var(ncid, REC_NAME, NF90_DOUBLE, rec_dimid, rec_varid) )
  CALL check( nf90_def_var(ncid, TIMESTEP_NAME, NF90_DOUBLE, varid=dt_varid) )

  ! Assign units attributes to coordinate variables.
  CALL check( nf90_put_att(ncid, pres_varid,"long_name",TRIM(PRES_LONGNAME)) )
  CALL check( nf90_put_att(ncid, pres_varid, UNITS, TRIM(PRES_UNITS)) )
  CALL check( nf90_put_att(ncid, x_varid,"long_name",TRIM(X_LONGNAME)) )
  CALL check( nf90_put_att(ncid, x_varid, UNITS, TRIM(X_UNITS)) )
  IF (.NOT.(Sphere)) THEN
  CALL check( nf90_put_att(ncid, x_varid, "axis",X_AXIS))
  END IF
  CALL check( nf90_put_att(ncid, y_varid,"long_name",TRIM(Y_LONGNAME)) )
  CALL check( nf90_put_att(ncid, y_varid, UNITS, TRIM(Y_UNITS)) )
  IF (.NOT.(Sphere)) THEN
  CALL check( nf90_put_att(ncid, y_varid, "axis",Y_AXIS))
  END IF
  CALL check( nf90_put_att(ncid, z_varid,"long_name",Z_LONGNAME) )
  CALL check( nf90_put_att(ncid, z_varid, UNITS, Z_UNITS) )
  CALL check( nf90_put_att(ncid, z_varid,"positive","up") )
  CALL check( nf90_put_att(ncid, zint_varid,"long_name",ZINT_LONGNAME) )
  CALL check( nf90_put_att(ncid, zint_varid, UNITS, ZINT_UNITS) )
  CALL check( nf90_put_att(ncid, zint_varid,"positive","up") )
  CALL check( nf90_put_att(ncid, rec_varid,"long_name",REC_NAME) )
  CALL check( nf90_put_att(ncid, rec_varid, UNITS, REC_UNITS) )
  CALL check( nf90_put_att(ncid, rec_varid,"calender","noleap") )
  CALL check( nf90_put_att(ncid, dt_varid,"long_name",TRIM(TIMESTEP_LONGNAME)) )
  CALL check( nf90_put_att(ncid, dt_varid, UNITS, TRIM(TIMESTEP_UNITS)) )

  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"Conventions","CF-1.0") )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"model","ASAM 2.8") )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"test_case",TEST_NAME) )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"horizontal_resolution",TRIM(HORIZ_RESOL)) )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"levels",TRIM(LEVEL_NAME)) )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"grid",TRIM(GRID_TYPE)) )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"native_grid",TRIM(NAT_GRIDTYPE)) )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"equation",TRIM(NAME_EQU)) )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"time_frequency",TRIM(TIME_FREQ)) )
  CALL check( nf90_put_att(ncid, NF90_GLOBAL,"description",TRIM(DESCRIPT)) )

  ! The dimids array is used to pass the dimids of the dimensions of
  ! the netCDF variables. Both of the netCDF variables we are creating
  ! share the same four dimensions. In Fortran, the unlimited
  ! dimension must come last on the list of dimids.
  dimids = (/ x_dimid, y_dimid, z_dimid, rec_dimid /)

  ! Define the netCDF variables for the pressure and temperature data.
  IF (uPosL > 0) THEN
    CALL check( nf90_def_var(ncid, U_NAME,    NF90_REAL, dimids, u_varid   ) )
  END IF
  IF (vPosL > 0) THEN
    CALL check( nf90_def_var(ncid, V_NAME,    NF90_REAL, dimids, v_varid   ) )
  END IF
  IF (wPosL > 0) THEN
    CALL check( nf90_def_var(ncid, W_NAME,    NF90_REAL, dimids, w_varid   ) )
  END IF
  IF (thPos > 0) THEN
    CALL check( nf90_def_var(ncid, TH_NAME,   NF90_REAL, dimids, th_varid  ) )
  END IF
  IF (EnPos > 0) THEN
    CALL check( nf90_def_var(ncid, EN_NAME,   NF90_REAL, dimids, En_varid  ) )
  END IF
  IF (rhoPos > 0) THEN
    CALL check( nf90_def_var(ncid, RHO_NAME,  NF90_REAL, dimids, rho_varid ) )
  END IF
  IF (tkePos > 0) THEN
    CALL check( nf90_def_var(ncid, TKE_NAME,  NF90_REAL, dimids, tke_varid ) )
  END IF
  IF (disPos > 0) THEN
    CALL check( nf90_def_var(ncid, DIS_NAME,  NF90_REAL, dimids, dis_varid ) )
  END IF
  IF (tkehPos > 0) THEN
    CALL check( nf90_def_var(ncid, TKEH_NAME, NF90_REAL, dimids, tkeh_varid) )
  END IF
  IF (tkevPos > 0) THEN
    CALL check( nf90_def_var(ncid, TKEV_NAME, NF90_REAL, dimids, tkev_varid) )
  END IF
  IF (lenPos > 0) THEN
    CALL check( nf90_def_var(ncid, LEN_NAME,  NF90_REAL, dimids, len_varid ) )
  END IF
  IF (RhoVPos > 0) THEN
    CALL check( nf90_def_var(ncid, RhoV_NAME,   NF90_REAL, dimids, RhoV_varid  ) )
  END IF
  IF (RhoCPos > 0) THEN
    CALL check( nf90_def_var(ncid, RhoC_NAME,   NF90_REAL, dimids, RhoC_varid  ) )
  END IF
  IF (RhoIPos > 0) THEN
    CALL check( nf90_def_var(ncid, RhoI_NAME,   NF90_REAL, dimids, RhoI_varid  ) )
  END IF
  IF (RhoRPos > 0) THEN
    CALL check( nf90_def_var(ncid, RhoR_NAME,    NF90_REAL, dimids, RhoR_varid     ) )
  END IF
  IF (PreOut   ) THEN
    CALL check( nf90_def_var(ncid, PRE_NAME,   NF90_REAL, dimids, pre_varid    ) )
  END IF
  IF (DiffOut  ) THEN
    CALL check( nf90_def_var(ncid, DIFF_NAME,  NF90_REAL, dimids, diff_varid   ) )
  END IF
  IF (DiffMomHOut  ) THEN !marcelk
    CALL check( nf90_def_var(ncid, DIFF_NAME,  NF90_REAL, dimids, diff_varid   ) )
  END IF
  IF (DiffMomVOut  ) THEN !marcelk
    CALL check( nf90_def_var(ncid, DIFF_NAME,  NF90_REAL, dimids, diff_varid   ) )
  END IF
  IF (DiffHOut ) THEN
    CALL check( nf90_def_var(ncid, DIFFH_NAME, NF90_REAL, dimids, diffh_varid  ) )
  END IF
  IF (DiffVOut ) THEN
    CALL check( nf90_def_var(ncid, DIFFV_NAME, NF90_REAL, dimids, diffv_varid  ) )
  END IF
  IF (DiffPotOut) THEN
    CALL check( nf90_def_var(ncid, DIFFPOT_NAME, NF90_REAL, dimids, diffpot_varid) )
  END IF
  IF (DiffMomOut) THEN
    CALL check( nf90_def_var(ncid, DIFFMOM_NAME, NF90_REAL, dimids, diffmom_varid) )
  END IF
  IF (Position('Dummy1')>0) THEN
    CALL check( nf90_def_var(ncid, DUMMY1_NAME, NF90_REAL, dimids, dummy1_varid) )
  END IF
  IF (Position('Tracer1')>0) THEN
    CALL check( nf90_def_var(ncid, TRACER1_NAME, NF90_REAL, dimids, tracer1_varid) )
  END IF
  IF (Position('Dummy2')>0) THEN
    CALL check( nf90_def_var(ncid, DUMMY2_NAME, NF90_REAL, dimids, dummy2_varid) )
  END IF
  IF (Position('Dummy3')>0) THEN
    CALL check( nf90_def_var(ncid, DUMMY3_NAME, NF90_REAL, dimids, dummy3_varid) )
  END IF
  IF (Position('Dummy4')>0) THEN
    CALL check( nf90_def_var(ncid, DUMMY4_NAME, NF90_REAL, dimids, dummy4_varid) )
  END IF

! IF (LenOutSpecial>0) THEN
!   DO j=1,LenOutSpecial
!     CALL check( nf90_def_var(ncid, NameOutSpecial(j), NF90_REAL, dimids, special_varid(j)) )
!   END DO
! END IF

  IF (ChemieOut) THEN
    DO j=1,NumAeroOut
      DO it=1,nFrac
        CALL check( nf90_def_var(ncid, AeroSpecies(j,it), NF90_REAL, dimids, aero_varid(j,it)) )
      END DO
    END DO
    DO j=1,NumGasOut
      CALL check( nf90_def_var(ncid, GasSpecies(j), NF90_REAL, dimids, gas_varid(j)) )
    END DO
  END IF

  ! Assign units attributes to the netCDF variables.
  IF (uPosL > 0) THEN
    CALL check( nf90_put_att(ncid, u_varid,    UNITS, U_UNITS   ) )
    CALL check( nf90_put_att(ncid, u_varid,"description", "zonal wind"            ) )
  END IF
  IF (vPosL > 0) THEN
    CALL check( nf90_put_att(ncid, v_varid,    UNITS, V_UNITS   ) )
    CALL check( nf90_put_att(ncid, v_varid,"description", "meridional wind"       ) )
  END IF
  IF (wPosL > 0) THEN
    CALL check( nf90_put_att(ncid, w_varid,    UNITS, W_UNITS   ) )
    CALL check( nf90_put_att(ncid, w_varid,"description", "vertical wind"       ) )
  END IF
  IF (thPos > 0) THEN
    CALL check( nf90_put_att(ncid, th_varid,   UNITS, TH_UNITS  ) )
    CALL check( nf90_put_att(ncid, th_varid,"description", "potential temperature") )
  END IF
  IF (EnPos > 0) THEN
    CALL check( nf90_put_att(ncid, th_varid,   UNITS, EN_UNITS  ) )
    CALL check( nf90_put_att(ncid, th_varid,"description", "total energy") )
  END IF
  IF (rhoPos > 0) THEN
    CALL check( nf90_put_att(ncid, rho_varid,  UNITS, RHO_UNITS ) )
    CALL check( nf90_put_att(ncid, rho_varid,"description", "density"             ) )
  END IF
  IF (tkePos > 0) THEN
    CALL check( nf90_put_att(ncid, tke_varid,  UNITS, TKE_UNITS ) )
    CALL check( nf90_put_att(ncid, tke_varid,"description", "turbulent kinetic energy"    ) )
  END IF
  IF (disPos > 0) THEN
    CALL check( nf90_put_att(ncid, dis_varid,  UNITS, DIS_UNITS ) )
    CALL check( nf90_put_att(ncid, dis_varid,"description", "dissipation of turbulent kinetic energy"   ) )
  END IF
  IF (tkehPos > 0) THEN
    CALL check( nf90_put_att(ncid, tkeh_varid, UNITS, TKEH_UNITS) )
    CALL check( nf90_put_att(ncid, tkeh_varid,"description", "horizontal component of tke" ) )
  END IF
  IF (tkevPos > 0) THEN
    CALL check( nf90_put_att(ncid, tkev_varid, UNITS, TKEV_UNITS) )
    CALL check( nf90_put_att(ncid, tkev_varid,"description", "vertical component of tke"   ) )
  END IF
  IF (lenPos > 0) THEN
    CALL check( nf90_put_att(ncid, len_varid,  UNITS, LEN_UNITS ) )
    CALL check( nf90_put_att(ncid, len_varid,"description", "turbulent lengthscale"       ) )
  END IF
  IF (RhoVPos > 0) THEN
    CALL check( nf90_put_att(ncid, RhoV_varid,   UNITS, RhoV_UNITS  ) )
    CALL check( nf90_put_att(ncid, RhoV_varid,"description", "water vapour content"        ) )
  END IF
  IF (RhoCPos > 0) THEN
    CALL check( nf90_put_att(ncid, RhoC_varid,   UNITS, RhoC_UNITS  ) )
    CALL check( nf90_put_att(ncid, RhoC_varid,"description", "cloud water content"         ) )
  END IF
  IF (RhoIPos > 0) THEN
    CALL check( nf90_put_att(ncid, RhoI_varid,   UNITS, RhoI_UNITS  ) )
    CALL check( nf90_put_att(ncid, RhoI_varid,"description", "cloud ice content"   ) )
  END IF
  IF (RhoRPos > 0) THEN
    CALL check( nf90_put_att(ncid, RhoR_varid,   UNITS, RhoR_UNITS  ) )
    CALL check( nf90_put_att(ncid, RhoR_varid,"description", ""   ) )
  END IF
  IF (PreOut   ) THEN
    CALL check( nf90_put_att(ncid, pre_varid,     UNITS, PRE_UNITS    ) )
    CALL check( nf90_put_att(ncid, pre_varid,"description", ""   ) )
  END IF
  IF (DiffOut  ) THEN
    CALL check( nf90_put_att(ncid, diff_varid,    UNITS, DIFF_UNITS   ) )
    CALL check( nf90_put_att(ncid, diff_varid,"description", "diffusion coefficient"            ) )
  END IF
 IF (DiffMomHOut  ) THEN !marcelk
    CALL check( nf90_put_att(ncid, diff_varid,    UNITS, DIFF_UNITS   ) )
    CALL check( nf90_put_att(ncid, diff_varid,"description", "diffusion coefficient"            ) )
  END IF
  IF (DiffMomVOut  ) THEN !marcelk
    CALL check( nf90_put_att(ncid, diff_varid,    UNITS, DIFF_UNITS   ) )
    CALL check( nf90_put_att(ncid, diff_varid,"description", "diffusion coefficient"            ) )
  END IF
  IF (DiffHOut ) THEN
    CALL check( nf90_put_att(ncid, diffh_varid,   UNITS, DIFFH_UNITS  ) )
    CALL check( nf90_put_att(ncid, diffh_varid,"description", "horizontal diffusion coefficient" ) )
  END IF
  IF (DiffVOut ) THEN
    CALL check( nf90_put_att(ncid, diffv_varid,   UNITS, DIFFV_UNITS  ) )
    CALL check( nf90_put_att(ncid, diffv_varid,"description", "vertical diffusion coefficient"   ) )
  END IF
  IF (DiffPotOut) THEN
    CALL check( nf90_put_att(ncid, diffpot_varid, UNITS, DIFFPOT_UNITS) )
    CALL check( nf90_put_att(ncid, diffpot_varid,"description", "diffusion coefficient of momentum") )
  END IF
  IF (DiffMomOut) THEN
    CALL check( nf90_put_att(ncid, diffmom_varid, UNITS, DIFFMOM_UNITS) )
    CALL check( nf90_put_att(ncid, diffmom_varid,"description", "diffusion coefficient of heat"  ) )
  END IF
  IF (Position('Dummy1')>0) THEN
    CALL check( nf90_put_att(ncid, dummy1_varid, UNITS, DUMMY1_UNITS) )
    CALL check( nf90_put_att(ncid, dummy1_varid,"description", "tracer 1"  ) )
  END IF
  IF (Position('Tracer1')>0) THEN
    CALL check( nf90_put_att(ncid, tracer1_varid, UNITS, TRACER1_UNITS) )
    CALL check( nf90_put_att(ncid, tracer1_varid,"description", "tracer 1"  ) )
  END IF
  IF (Position('Dummy2')>0) THEN
    CALL check( nf90_put_att(ncid, dummy2_varid, UNITS, DUMMY1_UNITS) )
    CALL check( nf90_put_att(ncid, dummy2_varid,"description", "tracer 2"  ) )
  END IF
  IF (Position('Dummy3')>0) THEN
    CALL check( nf90_put_att(ncid, dummy3_varid, UNITS, DUMMY1_UNITS) )
    CALL check( nf90_put_att(ncid, dummy3_varid,"description", "tracer 3"  ) )
  END IF
  IF (Position('Dummy4')>0) THEN
    CALL check( nf90_put_att(ncid, dummy4_varid, UNITS, DUMMY1_UNITS) )
    CALL check( nf90_put_att(ncid, dummy4_varid,"description", "tracer 4"  ) )
  END IF

! IF (LenOutSpecial>0) THEN
!   DO i=1,LenOutSpecial
!     CALL check( nf90_put_att(ncid, special_varid(i), UNITS, " ") )
!     CALL check( nf90_put_att(ncid, special_varid(i),"description", NameOutSpecial(i)  ) )
!   END DO
! END IF
  IF (ChemieOut) THEN
    DO j=1,NumGasOut
      CALL check( nf90_put_att(ncid, gas_varid(j), UNITS, GAS_UNITS) )
      CALL check( nf90_put_att(ncid, gas_varid(j),"description", "gas concentration"  ) )
    END DO
    DO j=1,NumAeroOut
      DO it=1,nFrac
        CALL check( nf90_put_att(ncid, aero_varid(j,it), UNITS, AERO_UNITS) )
        CALL check( nf90_put_att(ncid, aero_varid(j,it),"description", "aerosol concentration"  ) )
      END DO
    END DO
  END IF


  ! End define mode.
  CALL check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  CALL check( nf90_put_var(ncid, pres_varid, 1.d5) )
  DO i=1,NXX
    CALL check( nf90_put_var(ncid, x_varid, xLoc) )
  END DO
  DO i=1,NYY
    CALL check( nf90_put_var(ncid, y_varid, yLoc) )
  END DO
  DO i=1,NZZ
    CALL check( nf90_put_var(ncid, z_varid, zLoc) )
  END DO
  DO i=1,NZZ+1
    CALL check( nf90_put_var(ncid, zint_varid, z_int) )
  END DO
  DO i=0,NTT
    CALL check( nf90_put_var(ncid, rec_varid, timeLoc) )
  END DO
  CALL check( nf90_put_var(ncid, dt_varid, dtMax) )

  ! These settings tell netcdf to write one timestep of data. (The
  ! setting of start(4) inside the loop below tells netCDF which
  ! timestep to write.)
  count = (/ NXX, NYY, NZZ, 1 /)

  start = (/ 1, 1, 1, 1 /)

  ! Write the pretend data. This will write our surface pressure and
  ! surface temperature data. The arrays only hold one timestep worth
  ! of data. We will just rewrite the same data for each timestep. In
  ! a REAL :: application, the data would change between timesteps.

  DO rec = 1, NTMAX+1
    start(4) = rec
   IF (uPosL   > 0) THEN
    CALL check( nf90_put_var(ncid, u_varid,   &
                             OutNetcdf(uIndex,rec,ix0Out+1:ix1Out,  &
                                                  iy0Out+1:iy1Out,  &
                                                  iz0Out+1:iz1Out), & 
                             start = start,   &
                             count = count) )
   END IF
   IF (vPosL   > 0) THEN
    CALL check( nf90_put_var(ncid, v_varid,   &
                             OutNetcdf(vIndex,rec,ix0Out+1:ix1Out,  &
                                                  iy0Out+1:iy1Out,  &
                                                  iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (wPosL   > 0) THEN
    CALL check( nf90_put_var(ncid, w_varid,   &
                             OutNetcdf(wIndex,rec,ix0Out+1:ix1Out,  &
                                                  iy0Out+1:iy1Out,  &
                                                  iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (thPos   > 0) THEN
     CALL check( nf90_put_var(ncid, th_varid, &
                             OutNetcdf(thIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (EnPos   > 0) THEN
     CALL check( nf90_put_var(ncid, En_varid, &
                             OutNetcdf(EnIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (rhoPos  > 0) THEN
     CALL check( nf90_put_var(ncid, rho_varid, &
                             OutNetcdf(rhoIndex,rec,ix0Out+1:ix1Out,  &
                                                    iy0Out+1:iy1Out,  & 
                                                    iz0Out+1:iz1Out), &
                             start = start,    &
                             count = count) )
   END IF
!  TKE-Closure k-eps
   IF (tkePos  > 0) THEN
     CALL check( nf90_put_var(ncid, tke_varid, &
                             OutNetcdf(tkeIndex,rec,ix0Out+1:ix1Out,  &
                                                    iy0Out+1:iy1Out,  & 
                                                    iz0Out+1:iz1Out), &
                             start = start,    &
                             count = count) )
   END IF
   IF (disPos  > 0) THEN
     CALL check( nf90_put_var(ncid, dis_varid, &
                             OutNetcdf(disIndex,rec,ix0Out+1:ix1Out,  &
                                                    iy0Out+1:iy1Out,  & 
                                                    iz0Out+1:iz1Out), &
                             start = start,    &
                             count = count) )
   END IF
!  TKE-Closure following Herzog (ATHAM)
   IF (tkeHPos > 0) THEN
     CALL check( nf90_put_var(ncid, tkeh_varid, &
                             OutNetcdf(tkeHIndex,rec,ix0Out+1:ix1Out,  &
                                                     iy0Out+1:iy1Out,  & 
                                                     iz0Out+1:iz1Out), &
                             start = start,     &
                             count = count) )
   END IF
   IF (tkeVPos > 0) THEN
     CALL check( nf90_put_var(ncid, tkev_varid, &
                             OutNetcdf(tkeVIndex,rec,ix0Out+1:ix1Out,  &
                                                     iy0Out+1:iy1Out,  & 
                                                     iz0Out+1:iz1Out), &
                             start = start,     &
                             count = count) )
   END IF
   IF (LenPos  > 0) THEN
     CALL check( nf90_put_var(ncid, len_varid, &
                             OutNetcdf(LenIndex,rec,ix0Out+1:ix1Out,  &
                                                    iy0Out+1:iy1Out,  & 
                                                    iz0Out+1:iz1Out), &
                             start = start,    &
                             count = count) )
   END IF
!  Moisture
   IF (RhoVPos   > 0) THEN
     CALL check( nf90_put_var(ncid, RhoV_varid, &
                             OutNetcdf(RhoVIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (RhoCPos   > 0) THEN
     CALL check( nf90_put_var(ncid, RhoC_varid, &
                             OutNetcdf(RhoCIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (RhoIPos   > 0) THEN
     CALL check( nf90_put_var(ncid, RhoI_varid, &
                             OutNetcdf(RhoIIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (RhoRPos   > 0) THEN
     CALL check( nf90_put_var(ncid, RhoR_varid, &
                             OutNetcdf(RhoRIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
!  Additional Output
   IF (PreOut   ) THEN
     CALL check( nf90_put_var(ncid, pre_varid, &
                             OutNetcdf(PreIndex,rec,ix0Out+1:ix1Out,  &
                                                    iy0Out+1:iy1Out,  & 
                                                    iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (DiffOut  ) THEN
     CALL check( nf90_put_var(ncid, diff_varid, &
                             OutNetcdf(DiffIndex,rec,ix0Out+1:ix1Out,  &
                                                     iy0Out+1:iy1Out,  & 
                                                     iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (DiffMomHOut  ) THEN !marcelk
     CALL check( nf90_put_var(ncid, diff_varid, &
                             OutNetcdf(DiffMomHIndex,rec,ix0Out+1:ix1Out,  &
                                                     iy0Out+1:iy1Out,  &
                                                     iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (DiffMomVOut  ) THEN !marcelk
     CALL check( nf90_put_var(ncid, diff_varid, &
                             OutNetcdf(DiffMomVIndex,rec,ix0Out+1:ix1Out,  &
                                                     iy0Out+1:iy1Out,  &
                                                     iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (DiffHOut ) THEN
     CALL check( nf90_put_var(ncid, diffh_varid, &
                             OutNetcdf(DiffHIndex,rec,ix0Out+1:ix1Out,  &
                                                      iy0Out+1:iy1Out,  & 
                                                      iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (DiffVOut ) THEN
     CALL check( nf90_put_var(ncid, diffv_varid, &
                             OutNetcdf(DiffVIndex,rec,ix0Out+1:ix1Out,  &
                                                      iy0Out+1:iy1Out,  & 
                                                      iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (DiffPotOut) THEN
     CALL check( nf90_put_var(ncid, diffpot_varid, &
                             OutNetcdf(DiffPotIndex,rec,ix0Out+1:ix1Out,  &
                                                        iy0Out+1:iy1Out,  & 
                                                        iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (DiffMomOut) THEN
     CALL check( nf90_put_var(ncid, diffmom_varid, &
                             OutNetcdf(DiffMomIndex,rec,ix0Out+1:ix1Out,  &
                                                        iy0Out+1:iy1Out,  & 
                                                        iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (Position('Dummy1')>0) THEN
     CALL check( nf90_put_var(ncid, dummy1_varid, &
                             OutNetcdf(Dummy1Index,rec,ix0Out+1:ix1Out,  &
                                                        iy0Out+1:iy1Out,  & 
                                                        iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (Position('Tracer1')>0) THEN
     CALL check( nf90_put_var(ncid, tracer1_varid, &
                             OutNetcdf(Tracer1Index,rec,ix0Out+1:ix1Out,  &
                                                        iy0Out+1:iy1Out,  &
                                                        iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (Position('Dummy2')>0) THEN
     CALL check( nf90_put_var(ncid, dummy2_varid, &
                             OutNetcdf(Dummy2Index,rec,ix0Out+1:ix1Out,  &
                                                        iy0Out+1:iy1Out,  & 
                                                        iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (Position('Dummy3')>0) THEN
     CALL check( nf90_put_var(ncid, dummy3_varid, &
                             OutNetcdf(Dummy3Index,rec,ix0Out+1:ix1Out,  &
                                                        iy0Out+1:iy1Out,  & 
                                                        iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (Position('Dummy4')>0) THEN
     CALL check( nf90_put_var(ncid, dummy4_varid, &
                             OutNetcdf(Dummy4Index,rec,ix0Out+1:ix1Out,  &
                                                        iy0Out+1:iy1Out,  & 
                                                        iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
 ! IF (LenOutSpecial>0) THEN
 !   DO i=1,LenOutSpecial
 !     CALL check( nf90_put_var(ncid, special_varid(i), &
 !                             OutNetcdf(SpecialIndex(i),rec,ix0Out+1:ix1Out,  &
 !                                                        iy0Out+1:iy1Out,  &
 !                                                        iz0Out+1:iz1Out), &
 !                             start = start,   &
 !                             count = count) )
 !   END DO
 ! END IF
   IF (ChemieOut) THEN
     DO j=1,NumAeroOut
       DO it=1,nFrac
       CALL check( nf90_put_var(ncid, aero_varid(j,it), &
                               OutNetcdf(AeroIndex(j,it),rec,ix0Out+1:ix1Out,  &
                                                          iy0Out+1:iy1Out,  & 
                                                          iz0Out+1:iz1Out), &
                               start = start,   &
                               count = count) )
       END DO
     END DO
     DO j=1,NumGasOut
       CALL check( nf90_put_var(ncid, gas_varid(j), &
                               OutNetcdf(GasIndex(j),rec,ix0Out+1:ix1Out,  &
                                                         iy0Out+1:iy1Out,  & 
                                                         iz0Out+1:iz1Out), &
                               start = start,   &
                               count = count) )
     END DO
   END IF

  END DO

  ! Close the file. This causes netCDF to flush all buffers and make
  ! sure your data are REALly written to disk.
  CALL check( nf90_close(ncid) )

  DEALLOCATE(xLoc)
  DEALLOCATE(yLoc)
  DEALLOCATE(zLoc)

  print *,"*** SUCCESS writing netcdf-file ", FILE_NAME, "!"

END SUBROUTINE Netcdf_Output

  SUBROUTINE check(STATUS)
    INTEGER, INTENT ( in) :: STATUS

    IF(STATUS /= nf90_noerr) THEN
      PRINT *, trim(nf90_strerror(STATUS))
!      STOP "Stopped"
    END IF
  END SUBROUTINE check

SUBROUTINE InputDCMIP(FileName)
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line
  CHARACTER (20) :: Problem  
  INTEGER  :: ScalingFactor  
  REAL(RealKind) :: VelMax,N,phi0,lam0,nv,rpVo,DeltapVo,ProfFac,OmegaFac
  LOGICAL :: Perturbation
  NAMELIST /Example/    &

                    VelMax , &
                    ProfFac , &
                    Problem , &
                    Perturbation, &
                    OmegaFac, &
                    N, &
                    phi0, &
                    lam0, &
                    nv, &
                    rpVo, &
                    DeltapVo, &
                    ScalingFactor

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&Example')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=Example)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)

  DESCRIPT     = "" 
  NAME_EQU     = "nonhydro" 
  GRID_TYPE    = "latlon"
  NAT_GRIDTYPE = "latlon"
  X_NAME = "lon"
  Y_NAME = "lat"
  X_LONGNAME = "longitude"
  Y_LONGNAME = "latitude"
  X_UNITS = "degrees_east"
  Y_UNITS = "degrees_north"

!  CALL SetDomain(Domain)
  SELECT CASE(Problem)
    CASE ('BaroIn','DCMIP_4_2')
      TIME_FREQ    = "6hr" 
      IF (ix1==360) THEN
        HORIZ_RESOL  = "medium" 
      ELSE IF (ix1==720) THEN
        HORIZ_RESOL  = "high" 
      END IF
      LEVEL_NAME   = "L30"
      TEST_NAME    = "42"
    CASE ('DCMIP_1_2')
      TEST_NAME    = "12"
      TIME_FREQ    = "1hr" 
      IF (Domain%iz1==30) THEN
        HORIZ_RESOL  = "low" 
        LEVEL_NAME   = "L30"
      ELSE IF (Domain%iz1==60) THEN
        HORIZ_RESOL  = "medium" 
        LEVEL_NAME   = "L60"
      ELSE IF (Domain%iz1==120) THEN
        HORIZ_RESOL  = "high" 
        LEVEL_NAME   = "L120"
      END IF
    CASE ('DCMIP_1_3')
      TEST_NAME    = "13"
      TIME_FREQ    = "day" 
      HORIZ_RESOL  = "medium" 
      IF (Domain%iz1==120) THEN
        LEVEL_NAME   = "L120"
      ELSE IF (Domain%iz1==60) THEN
        LEVEL_NAME   = "L60"
      ELSE IF (Domain%iz1==30) THEN
        LEVEL_NAME   = "L30"
      END IF
    CASE ('DCMIP_2_1')
      TEST_NAME    = "21"
      TIME_FREQ    = "100 sec" 
      HORIZ_RESOL  = "medium" 
      LEVEL_NAME   = "L60"
    CASE ('DCMIP_2_2')
      TEST_NAME    = "22"
      TIME_FREQ    = "100 sec" 
      HORIZ_RESOL  = "medium" 
      LEVEL_NAME   = "L60"
    CASE ('DCMIP_1_1')
      TIME_FREQ    = "day" 
      HORIZ_RESOL  = "medium" 
      LEVEL_NAME   = "L60"
      TEST_NAME    = "11"
    CASE ('DCMIP_3')
      TIME_FREQ    = "day" 
      HORIZ_RESOL  = "medium" 
      LEVEL_NAME   = "L30"
      TEST_NAME    = "31"
    CASE ('DCMIP_4_1')
      TIME_FREQ    = "day" 
      HORIZ_RESOL  = "medium" 
      LEVEL_NAME   = "L30"
      IF (ScalingFactor==1) THEN
        TEST_NAME    = "410"
      ELSE IF (ScalingFactor==10) THEN
        TEST_NAME    = "411"
      ELSE IF (ScalingFactor==100) THEN
        TEST_NAME    = "412"
      ELSE IF (ScalingFactor==1000) THEN
        TEST_NAME    = "413"
      END IF
  END SELECT
END SUBROUTINE InputDCMIP


END PROGRAM MainRead

