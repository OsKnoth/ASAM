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

  USE Output_Mod
  USE netcdf

  IMPLICIT NONE

  TYPE(VelocityFace_T)  , POINTER :: VelF1(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecT(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecT1(:)
  TYPE(Vector4Cell_T)   , POINTER :: VecTP(:)
  TYPE(VecVector4Cell_T), POINTER :: VecVelR(:)
  TYPE(PressureVelocity), POINTER :: bbb(:)

  CHARACTER(80)                   :: ProfileSTART='ProfileSTART'
  CHARACTER(80)                   :: ProfileEND='ProfileEND'
  CHARACTER(80)                   :: InputFileName
  CHARACTER(80)                   :: ReadFileName
  CHARACTER(80)                   :: ReadFileName1
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

  INTEGER        :: ib,ibLoc
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
  CALL Allocate(JacMet)
  JacMet=Zero
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
    CALL InputChemicalData(DataFile)
    CALL OutputChemie('Output')
! END IF

! -- Meteorologie --
  CALL ALLOCATE(VelF1)
  VelF1=Zero
  VecT=Zero
  VecT1=Zero

! -- Allokieren der Jacobi-Matrix f�r Transport
! CALL Allocate(JacMet)

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
    ELSE IF (TkeHVLen) THEN
      CALL ALLOCATE     (DiffHKoeff)
      CALL ALLOCATE     (DiffVKoeff)
    ELSE
      CALL ALLOCATE     (DiffKoeff)
    END IF
  END IF


  CALL InputModelOutput(InputFileName)
  CALL InitOutputSpecial(InputFileName)
  WRITE(*,*) 'Input ReadFileName'
  READ(*,*) ReadFileName
  WRITE(*,*) 'Input ReadFileName1'
  READ(*,*) ReadFileName1
  CALL InputMet(Time,VecT,VecT,ReadFileName)
  CALL InputMet(Time,VecT1,VecT1,ReadFileName1)
  CALL AXPY(-One,VecT1,VecT)

! ------------------------------------------   
! - Number of Files -
! ------------------------------------------   
  NumberInputFiles = INT(EndTime/OutputTimeStep)

! ------------------------------------------
! - User specific Number of Inputfiles - 
!
! - You can comment that part if You don't
!   like to change the OutputTimeStep!
! ------------------------------------------
  WRITE(*,*)'########################################################'
  WRITE(*,*)'User specific Number of Timesteps used for NETCDF-Output'
  WRITE(*,*)' '
  WRITE(*,*)'  EndTime        = ',EndTime
  WRITE(*,*)'  OutputTimeStep = ',OutputTimeStep
  WRITE(*,*)' '
  WRITE(*,*)'NETCDF_TimeStep - A multiple integer of OutputTimeStep'
  WRITE(*,*)'NETCDF_TimeStep  >=  OutputTimeStep'
  WRITE(*,*)'########################################################'
  WRITE(*,*)' '
  WRITE(*,*)'NETCDF_TimeStep =  ... [s]'
  READ(*,'(I4)') NETCDF_TimeStep
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
  DtInputFiles = INT(NETCDF_TimeStep/OutputTimeStep)
  NTMAX        = INT(EndTime/NETCDF_TimeStep)

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
    qvIndex   = Pos1
  END IF
  IF (RhoCpos   > 0) THEN
    Pos1      = Pos1 + 1
    qcIndex   = Pos1
  END IF
  IF (RhoIPos   > 0) THEN
    Pos1      = Pos1 + 1
    qiIndex   = Pos1
  END IF
  IF (RhorPos   > 0) THEN
    Pos1      = Pos1 + 1
    qrIndex   = Pos1
  END IF

!  IF (PreOut     ) THEN
!    Pos1          = Pos1 + 1
!    PreIndex      = Pos1
!  END IF
  IF (DiffOut    ) THEN
    Pos1          = Pos1 + 1
    DiffIndex     = Pos1
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
    record = INT(999 + i)
    FileIndex = FileIndex + 1
! ------------------------------------------
! - Read Files -
! ------------------------------------------
    WRITE(record_string,'(I4)') record
    ReadFileName = TRIM(ReadFileName(1:INDEX(ReadFileName,'.out')-1)) &
                                    //'.out'//ADJUSTL(record_string)
    WRITE(*,*)'Reading ... ',ReadFileName
    
    CALL InputMet(Time,VecT,VecT,ReadFileName)
    IF (ProfOut) THEN
      CALL OutputProfile(VecT,VecT,ProfileEND)
    END IF

! ------------------------------------------
! - Set OutNetcdf -
! ------------------------------------------
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
    IF (rhoPos  > 0) THEN
      Pos = rhoPos
      Ind = rhoIndex
      WRITE(*,*) 'rhoPos',rhoPos
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
      Ind = qvIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (RhoCpos   > 0) THEN
      Pos = RhoCpos
      Ind = qcIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (RhoIPos   > 0) THEN
      Pos = RhoIPos
      Ind = qiIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF
    IF (RhorPos   > 0) THEN
      Pos = RhorPos
      Ind = qrIndex
      CALL SetOutNetcdfVec(VecT,OutNetcdf,Pos,Ind,FileIndex)
    END IF

!   Additional output
    CALL SetOutNetcdf(OutNetcdf,FileIndex)

  END DO

! ------------------------------------------
! - Output Netcdf -
! ------------------------------------------
  CALL SetDomain(Domain)

  CALL Netcdf_Output(OutNetcdf,InputFileName)

  DEALLOCATE(OutNetcdf)

CONTAINS

SUBROUTINE SetOutNetcdfVec(VecC,OutC,PosC,IndC,FileIndex)

  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  INTEGER :: PosC, IndC, time
  INTEGER :: FileIndex
  REAL(RealKind) :: OutC(:,:,:,:,:)
  
  DO ibLoc=1,nbLoc
     ib=LocGlob(ibLoc)
     CALL DomainSet(ib)
     DO iz=igz0+1,igz1
       DO iy=igy0+1,igy1
         DO ix=igx0+1,igx1
            OutC(IndC,FileIndex,ix,iy,iz)=VecC(ib)%Vec(PosC)%c(ix,iy,iz,1)
         END DO
       END DO
     END DO
  END DO

END SUBROUTINE SetOutNetcdfVec

SUBROUTINE SetOutNetcdf(OutC,FileIndex)

  IMPLICIT NONE

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  INTEGER :: IndC
  INTEGER :: FileIndex
  REAL(RealKind) :: OutC(:,:,:,:,:)
    
  DO ibLoc=1,nbLoc
     ib=LocGlob(ibLoc)
     CALL DomainSet(ib)
     DO iz=igz0+1,igz1
       DO iy=igy0+1,igy1
         DO ix=igx0+1,igx1
!           IF (PreOut     ) THEN
!             IndC = PreIndex
!             OutC(IndC,FileIndex,ix,iy,iz)=PreCell(ib)%c(ix,iy,iz,1)
!           END IF
           IF (DiffOut    ) THEN
             IndC = DiffIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffKoeff(ib)%c(ix,iy,iz,1)
           END IF
           IF (DiffHOut   ) THEN
             IndC = DiffHIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffHKoeff(ib)%c(ix,iy,iz,1)
           END IF
           IF (DiffVOut   ) THEN
             IndC = DiffVIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffVKoeff(ib)%c(ix,iy,iz,1)
           END IF
           IF (DiffPotOut ) THEN
             IndC = DiffPotIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffPotKoeff(ib)%c(ix,iy,iz,1)
           END IF
           IF (DiffMomOut ) THEN
             IndC = DiffMomIndex
             OutC(IndC,FileIndex,ix,iy,iz)=DiffMomKoeff(ib)%c(ix,iy,iz,1)
           END IF
         END DO
       END DO
     END DO
  END DO

END SUBROUTINE SetOutNetcdf

SUBROUTINE Netcdf_Output(OutNetcdf,InputFileName)
  USE Floor_mod
  USE netcdf
  IMPLICIT NONE

  ! This is the name of the data file we will create.
  CHARACTER (80) :: FILE_NAME 
  CHARACTER (80) :: InputFileName
  INTEGER :: ncid

  INTEGER,             PARAMETER :: NDIMS  = 4
  CHARACTER (LEN = *), PARAMETER :: X_NAME = "x"
  CHARACTER (LEN = *), PARAMETER :: Y_NAME = "y"
  CHARACTER (LEN = *), PARAMETER :: Z_NAME = "z"
  CHARACTER (LEN = *), PARAMETER :: REC_NAME = "time"
  INTEGER :: x_dimid, y_dimid, z_dimid, rec_dimid
  INTEGER :: rec

  ! The start and count arrays will tell the netCDF library where to
  ! write our data.
  INTEGER :: start(NDIMS), count(NDIMS)

  ! These program variables hold the cartesian grid
  INTEGER :: x_varid, y_varid, z_varid

  ! We will create two netCDF variables, one each for temperature and
  ! pressure fields.
  CHARACTER (LEN = *), PARAMETER :: U_NAME         ="U"
  CHARACTER (LEN = *), PARAMETER :: V_NAME         ="V"
  CHARACTER (LEN = *), PARAMETER :: W_NAME         ="W"
  CHARACTER (LEN = *), PARAMETER :: TH_NAME        ="THETA"
  CHARACTER (LEN = *), PARAMETER :: RHO_NAME       ="RHO"
  CHARACTER (LEN = *), PARAMETER :: TKE_NAME       ="TKE"
  CHARACTER (LEN = *), PARAMETER :: DIS_NAME       ="DIS"
  CHARACTER (LEN = *), PARAMETER :: TKEH_NAME      ="TKEH"
  CHARACTER (LEN = *), PARAMETER :: TKEV_NAME      ="TKEV"
  CHARACTER (LEN = *), PARAMETER :: LEN_NAME       ="LEN"
  CHARACTER (LEN = *), PARAMETER :: QV_NAME        ="QV"
  CHARACTER (LEN = *), PARAMETER :: QC_NAME        ="QC"
  CHARACTER (LEN = *), PARAMETER :: QI_NAME        ="QI"
  CHARACTER (LEN = *), PARAMETER :: QR_NAME        ="QR"
  CHARACTER (LEN = *), PARAMETER :: PRE_NAME       ="PRE"
  CHARACTER (LEN = *), PARAMETER :: DIFF_NAME      ="DIFF"
  CHARACTER (LEN = *), PARAMETER :: DIFFH_NAME     ="DIFFH"
  CHARACTER (LEN = *), PARAMETER :: DIFFV_NAME     ="DIFFV"
  CHARACTER (LEN = *), PARAMETER :: DIFFPOT_NAME   ="DIFFPOT"
  CHARACTER (LEN = *), PARAMETER :: DIFFMOM_NAME   ="DIFFMOM"
  INTEGER :: u_varid, v_varid, w_varid
  INTEGER :: th_varid, rho_varid
  INTEGER :: tke_varid, dis_varid
  INTEGER :: tkeh_varid, tkev_varid, len_varid
  INTEGER :: qv_varid, qc_varid, qi_varid, qr_varid
  INTEGER :: pre_varid, diff_varid, diffh_varid, diffv_varid
  INTEGER :: diffpot_varid, diffmom_varid
  INTEGER :: dimids(NDIMS)

  ! We recommend that each variable carry a "units" attribute.
  CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
  CHARACTER (LEN = *), PARAMETER :: X_UNITS    = "m"
  CHARACTER (LEN = *), PARAMETER :: Y_UNITS    = "m"
  CHARACTER (LEN = *), PARAMETER :: Z_UNITS    = "m"
  CHARACTER (LEN = *), PARAMETER :: U_UNITS    = "m/s"
  CHARACTER (LEN = *), PARAMETER :: V_UNITS    = "m/s"
  CHARACTER (LEN = *), PARAMETER :: W_UNITS    = "m/s"
  CHARACTER (LEN = *), PARAMETER :: TH_UNITS   = "K"
  CHARACTER (LEN = *), PARAMETER :: RHO_UNITS  = "kg/m**3"
  CHARACTER (LEN = *), PARAMETER :: TKE_UNITS  = "m**2/s**2"
  CHARACTER (LEN = *), PARAMETER :: DIS_UNITS  = "m**2/s**2"
  CHARACTER (LEN = *), PARAMETER :: TKEH_UNITS = "m**2/s**2"
  CHARACTER (LEN = *), PARAMETER :: TKEV_UNITS = "m**2/s**2"
  CHARACTER (LEN = *), PARAMETER :: LEN_UNITS  = "m"
  CHARACTER (LEN = *), PARAMETER :: QV_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: QC_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: QI_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: QR_UNITS   = " "
  CHARACTER (LEN = *), PARAMETER :: PRE_UNITS       = " "
  CHARACTER (LEN = *), PARAMETER :: DIFF_UNITS      = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DIFFH_UNITS     = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DIFFV_UNITS     = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DIFFPOT_UNITS   = "m**2/s "
  CHARACTER (LEN = *), PARAMETER :: DIFFMOM_UNITS   = "m**2/s "

  ! -----------------------------------
  ! - ASAM -
  ! -----------------------------------
  REAL(RealKind) :: OutNetcdf(:,:,:,:,:)

  !---  local variables
  INTEGER :: ix, iy, iz, it, i
  INTEGER :: ib, ibLoc,iScalar,NumOut
  INTEGER :: NXX, NYY, NZZ
  INTEGER, ALLOCATABLE :: xLoc(:),yLoc(:),zLoc(:)

  ! Create the filename of the Netcdf-Outputfile
  FILE_NAME = TRIM(InputFileName(1:INDEX(InputFileName,'.grid')-1))//'.cdf'

  ! Program variables to hold the data we will write out. 
  NXX = ix1Out - ix0Out
  NYY = iy1Out - iy0Out
  NZZ = iz1Out - iz0Out
  
  ALLOCATE(xLoc(1:NXX))
  ALLOCATE(yLoc(1:NYY))
  ALLOCATE(zLoc(1:NZZ))

  DO ix = 1, NXX
     xLoc(ix) = ix
  END DO
  DO iy = 1, NYY
     yLoc(iy) = iy
  END DO
  DO iz = 1, NZZ
     zLoc(iz) = iz
  END DO

  ! Cr
  CALL check( nf90_create(FILE_NAME, nf90_clobber, ncid) )

  ! Define the dimensions. The record dimension is defined to have
  ! unlimited length - it can grow as needed. In this example it is
  ! the time dimension.
  CALL check( nf90_def_dim(ncid, X_NAME, NXX, x_dimid) )
  CALL check( nf90_def_dim(ncid, Y_NAME, NYY, y_dimid) )
  CALL check( nf90_def_dim(ncid, Z_NAME, NZZ, z_dimid) )
  CALL check( nf90_def_dim(ncid, REC_NAME, NF90_UNLIMITED, rec_dimid) )

  ! Define the coordinate variables. We will only define coordinate
  ! variables for lat and lon.  Ordinarily we would need to provide
  ! an array of dimension IDs for each variable's dimensions, but
  ! since coordinate variables only have one dimension, we can
  ! simply provide the address of that dimension ID (lat_dimid) and
  ! similarly for (lon_dimid).
  CALL check( nf90_def_var(ncid, X_NAME, NF90_REAL, x_dimid, x_varid) )
  CALL check( nf90_def_var(ncid, Y_NAME, NF90_REAL, y_dimid, y_varid) )
  CALL check( nf90_def_var(ncid, Z_NAME, NF90_REAL, z_dimid, z_varid) )

  ! Assign units attributes to coordinate variables.
  CALL check( nf90_put_att(ncid, x_varid, UNITS, X_UNITS) )
  CALL check( nf90_put_att(ncid, y_varid, UNITS, Y_UNITS) )
  CALL check( nf90_put_att(ncid, z_varid, UNITS, Z_UNITS) )

  ! The dimids array is used to pass the dimids of the dimensions of
  ! the netCDF variables. Both of the netCDF variables we are creating
  ! share the same four dimensions. In Fortran, the unlimited
  ! dimension must come last on the list of dimids.
  dimids = (/ x_dimid, y_dimid, z_dimid, rec_dimid /)

  ! Define the netCDF variables for the pressure and temperature data.
  IF (uPosL > 0) THEN
  CALL check( nf90_def_var(ncid, U_NAME,    NF90_DOUBLE, dimids, u_varid   ) )
  END IF
  IF (vPosL > 0) THEN
  CALL check( nf90_def_var(ncid, V_NAME,    NF90_DOUBLE, dimids, v_varid   ) )
  END IF
  IF (wPosL > 0) THEN
  CALL check( nf90_def_var(ncid, W_NAME,    NF90_DOUBLE, dimids, w_varid   ) )
  END IF
  IF (thPos > 0) THEN
  CALL check( nf90_def_var(ncid, TH_NAME,   NF90_DOUBLE, dimids, th_varid  ) )
  END IF
  IF (rhoPos > 0) THEN
  CALL check( nf90_def_var(ncid, RHO_NAME,  NF90_DOUBLE, dimids, rho_varid ) )
  END IF
  IF (tkePos > 0) THEN
  CALL check( nf90_def_var(ncid, TKE_NAME,  NF90_DOUBLE, dimids, tke_varid ) )
  END IF
  IF (disPos > 0) THEN
  CALL check( nf90_def_var(ncid, DIS_NAME,  NF90_DOUBLE, dimids, dis_varid ) )
  END IF
  IF (tkehPos > 0) THEN
  CALL check( nf90_def_var(ncid, TKEH_NAME, NF90_DOUBLE, dimids, tkeh_varid) )
  END IF
  IF (tkevPos > 0) THEN
  CALL check( nf90_def_var(ncid, TKEV_NAME, NF90_DOUBLE, dimids, tkev_varid) )
  END IF
  IF (lenPos > 0) THEN
  CALL check( nf90_def_var(ncid, LEN_NAME,  NF90_DOUBLE, dimids, len_varid ) )
  END IF
  IF (RhoVPos > 0) THEN
  CALL check( nf90_def_var(ncid, QV_NAME,   NF90_DOUBLE, dimids, qv_varid  ) )
  END IF
  IF (RhoCpos > 0) THEN
  CALL check( nf90_def_var(ncid, QC_NAME,   NF90_DOUBLE, dimids, qc_varid  ) )
  END IF
  IF (RhoIPos > 0) THEN
  CALL check( nf90_def_var(ncid, QI_NAME,   NF90_DOUBLE, dimids, qi_varid  ) )
  END IF
  IF (RhorPos > 0) THEN
  CALL check( nf90_def_var(ncid, QR_NAME,      NF90_DOUBLE, dimids, qr_varid     ) )
  END IF
!  IF (PreOut   ) THEN
!  CALL check( nf90_def_var(ncid, PRE_NAME,     NF90_DOUBLE, dimids, pre_varid    ) )
!  END IF
  IF (DiffOut  ) THEN
  CALL check( nf90_def_var(ncid, DIFF_NAME,    NF90_DOUBLE, dimids, diff_varid   ) )
  END IF
  IF (DiffHOut ) THEN
  CALL check( nf90_def_var(ncid, DIFFH_NAME,   NF90_DOUBLE, dimids, diffh_varid  ) )
  END IF
  IF (DiffVOut ) THEN
  CALL check( nf90_def_var(ncid, DIFFV_NAME,   NF90_DOUBLE, dimids, diffv_varid  ) )
  END IF
  IF (DiffPotOut) THEN
  CALL check( nf90_def_var(ncid, DIFFPOT_NAME, NF90_DOUBLE, dimids, diffpot_varid) )
  END IF
  IF (DiffMomOut) THEN
  CALL check( nf90_def_var(ncid, DIFFMOM_NAME, NF90_DOUBLE, dimids, diffmom_varid) )
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
  CALL check( nf90_put_att(ncid, qv_varid,   UNITS, QV_UNITS  ) )
  CALL check( nf90_put_att(ncid, qv_varid,"description", "water vapour content"        ) )
  END IF
  IF (RhoCpos > 0) THEN
  CALL check( nf90_put_att(ncid, qc_varid,   UNITS, QC_UNITS  ) )
  CALL check( nf90_put_att(ncid, qc_varid,"description", "cloud water content"         ) )
  END IF
  IF (RhoIPos > 0) THEN
  CALL check( nf90_put_att(ncid, qi_varid,   UNITS, QI_UNITS  ) )
  CALL check( nf90_put_att(ncid, qi_varid,"description", "cloud ice content"   ) )
  END IF
  IF (RhorPos > 0) THEN
  CALL check( nf90_put_att(ncid, qr_varid,   UNITS, QR_UNITS  ) )
  CALL check( nf90_put_att(ncid, qr_varid,"description", ""   ) )
  END IF
!  IF (PreOut   ) THEN
!  CALL check( nf90_put_att(ncid, pre_varid,     UNITS, PRE_UNITS    ) )
!  CALL check( nf90_put_att(ncid, pre_varid,"description", ""   ) )
!  END IF
  IF (DiffOut  ) THEN
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
  CALL check( nf90_put_att(ncid, diffmom_varid,"description", "diffusion coefficient of heat"    ) )
  END IF


  ! End define mode.
  CALL check( nf90_enddef(ncid) )

  ! Write the coordinate variable data. This will put the latitudes
  ! and longitudes of our data grid into the netCDF file.
  DO i=1,NXX
     CALL check( nf90_put_var(ncid, x_varid, xLoc) )
  END DO
  DO i=1,NYY
     CALL check( nf90_put_var(ncid, y_varid, yLoc) )
  END DO
  DO i=1,NZZ
     CALL check( nf90_put_var(ncid, z_varid, zLoc) )
  END DO

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
     CALL check( nf90_put_var(ncid, qv_varid, &
                             OutNetcdf(qvIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (RhoCpos   > 0) THEN
     CALL check( nf90_put_var(ncid, qc_varid, &
                             OutNetcdf(qcIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (RhoIPos   > 0) THEN
     CALL check( nf90_put_var(ncid, qi_varid, &
                             OutNetcdf(qiIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
   IF (RhorPos   > 0) THEN
     CALL check( nf90_put_var(ncid, qr_varid, &
                             OutNetcdf(qrIndex,rec,ix0Out+1:ix1Out,  &
                                                   iy0Out+1:iy1Out,  & 
                                                   iz0Out+1:iz1Out), &
                             start = start,   &
                             count = count) )
   END IF
!  Additional Output
!   IF (PreOut   ) THEN
!     CALL check( nf90_put_var(ncid, pre_varid, &
!                             OutNetcdf(PreIndex,rec,ix0Out+1:ix1Out,  &
!                                                    iy0Out+1:iy1Out,  & 
!                                                    iz0Out+1:iz1Out), &
!                             start = start,   &
!                             count = count) )
!   END IF
   IF (DiffOut  ) THEN
     CALL check( nf90_put_var(ncid, diff_varid, &
                             OutNetcdf(DiffIndex,rec,ix0Out+1:ix1Out,  &
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

  END DO

  ! Close the file. This causes netCDF to flush all buffers and make
  ! sure your data are REALly written to disk.
  CALL check( nf90_close(ncid) )

  DEALLOCATE(xLoc)
  DEALLOCATE(yLoc)
  DEALLOCATE(zLoc)

  print *,"*** SUCCESS writing netcdf-file ", FILE_NAME, "!"

END SUBROUTINE Netcdf_Output

!CONTAINS
  SUBROUTINE check(STATUS)
    INTEGER, INTENT ( in) :: STATUS

    IF(STATUS /= nf90_noerr) THEN
      PRINT *, trim(nf90_strerror(STATUS))
      STOP "Stopped"
    END IF
  END SUBROUTINE check


END PROGRAM MainRead
