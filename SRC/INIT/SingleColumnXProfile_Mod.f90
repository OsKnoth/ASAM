MODULE SingleColumnXProfile_Mod

  USE Kind_Mod
  USE typesizes
  USE netcdf
  USE Domain_Mod
  USE Floor_Mod
  USE DataType_Mod
  USE Physics_Mod
  USE ScalarVectorCellPar_Mod
  USE Output_Mod
  USE Parallel_Mod
  IMPLICIT NONE

  INTEGER, PARAMETER :: NDIMS  = 2
  INTEGER, PRIVATE :: dimids(NDIMS)
  INTEGER, PRIVATE :: TimeStep
  INTEGER :: ScalarPos
  INTEGER :: MaxLenTable=0
  INTEGER, PRIVATE :: i
  CHARACTER(30)     :: NamePoint = " "
  CHARACTER(30)     :: NamePointFile = " "
  CHARACTER (LEN = *), PARAMETER :: TH_NAME           ="THETA"
  CHARACTER (LEN = *), PARAMETER :: Z_NAME            = "lev"
  CHARACTER (LEN = *), PARAMETER :: Z_LONGNAME        = "Height at midpoints"
  CHARACTER (LEN = *), PARAMETER :: REC_NAME          = "time"
  CHARACTER (LEN = *), PARAMETER :: REC_LONGNAME      = "time"
  CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
  CHARACTER (LEN = *), PARAMETER :: Z_UNITS    = "m"
  CHARACTER (LEN = *), PARAMETER :: REC_UNITS  = "seconds since 2000-01-01 00:00:00"
  REAL(RealKind), PRIVATE :: OutputTimeColumnX
  REAL(RealKind) :: dtimeColumnX=0.0d0
  REAL(RealKind) :: OutputTimeColumnXStart=1.0d40
  REAL(RealKind) :: OutputTimeColumnXEnd=-1.0d040

  TYPE ColumnXProfileTable_T
    INTEGER :: ibLoc
    CHARACTER(32) :: ScalarC=' '
    LOGICAL       :: ColumnXProfile=.FALSE.
    CHARACTER(20) :: Units=''
    CHARACTER(30) :: Description=''
    INTEGER       :: iPosC=0
  END TYPE ColumnXProfileTable_T
  TYPE ColumnXPoint_T
    REAL(RealKind) :: z
    REAL(RealKind) :: y
    REAL(RealKind) :: x
  END TYPE ColumnXPoint_T
  TYPE Profile_T
    REAL(RealKind), ALLOCATABLE :: Profile(:,:)
  END TYPE Profile_T  
  TYPE ColumnXProfile_T
    TYPE(ColumnXProfileTable_T), POINTER :: TableC(:)
    TYPE(ColumnXPoint_T) :: P
    TYPE(Profile_T), ALLOCATABLE :: ProfileBlock(:)
    INTEGER :: nbLoc=0
    INTEGER, ALLOCATABLE :: izProf(:)
    INTEGER, ALLOCATABLE :: iyProf(:)
    INTEGER, ALLOCATABLE :: ibLoc(:)
    LOGICAL, ALLOCATABLE :: OutBlock(:)
    LOGICAL :: OutHeader=.FALSE.
    CHARACTER*50 :: FileName
  END TYPE  
  TYPE ColumnXProfileWrite_T
    INTEGER :: ibLoc
    TYPE(ColumnXProfileTable_T), POINTER :: TableC(:)
    TYPE(ColumnXPoint_T) :: P
    REAL(RealKind), ALLOCATABLE :: Profile(:,:)
    INTEGER :: nb=0
    INTEGER :: nProc
    INTEGER, ALLOCATABLE :: izProf(:)
    INTEGER, ALLOCATABLE :: iyProf(:)
    INTEGER, ALLOCATABLE :: ib(:)
    INTEGER, ALLOCATABLE :: iProc(:)
    LOGICAL, ALLOCATABLE :: OutBlock(:)
    LOGICAL :: OutHeader=.FALSE.
    CHARACTER*50 :: FileName
    INTEGER :: ncid
    INTEGER :: x_varid, rec_varid, pz_varid, py_dimid
    INTEGER :: x_dimid, rec_dimid, pz_dimid, py_varid
    INTEGER, POINTER :: Scalar_VarId(:,:)
  END TYPE  
  TYPE(ColumnXProfile_T), ALLOCATABLE :: ColumnXProfile(:)
  TYPE(ColumnXProfileWrite_T), ALLOCATABLE :: ColumnXProfileWrite(:)

  NAMELIST /ColumnXProfileOutput/ OutputTimeColumnXStart, &
                                 OutputTimeColumnXEnd, &
                                 dtimeColumnX, &
                                 NamePointFile
CONTAINS

SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS

  IF(STATUS /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(STATUS))
  END IF
END SUBROUTINE check

SUBROUTINE OpenColumnXProfile(FileName)
  CHARACTER(*) :: FileName

  REAL(RealKind) :: xM(Domain%nx)
  INTEGER :: iz,iy,ix,izProf,iyProf,ixProf
  INTEGER :: iProf,ibProfLoc,ibProf
  TYPE(ColumnXProfileTable_T), POINTER :: TableC(:)
  INTEGER, POINTER :: Scalar_VarId(:,:)
  INTEGER :: iPz,iPy
  CHARACTER*12 :: NamePz,NamePy
  INTEGER :: ncid
  INTEGER :: x_varid, rec_varid, pz_varid, py_varid
  INTEGER :: x_dimid, rec_dimid, pz_dimid, py_dimid
  REAL(RealKind) :: yCProf,zCProf

! Determine number of local blocks
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iProf=1,SIZE(ColumnXProfile)
      zCProf=ColumnXProfile(iProf)%P%z
      yCProf=ColumnXProfile(iProf)%P%y
      S1:DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          IF (zP(iz-1)<=zCProf.AND. &
              zP(iz  )>=zCProf.AND. &
              yP(iy-1)<=yCProf.AND. &
              yP(iy  )>=yCProf) THEN
            ColumnXProfile(iProf)%nbLoc=ColumnXProfile(iProf)%nbLoc+1
            EXIT S1  
          END IF
        END DO 
      END DO S1
    END DO
  END DO
  DO iProf=1,SIZE(ColumnXProfile)
    IF (ColumnXProfile(iProf)%nbLoc>0) THEN
      ALLOCATE(ColumnXProfile(iProf)%izProf(ColumnXProfile(iProf)%nbLoc))
      ALLOCATE(ColumnXProfile(iProf)%iyProf(ColumnXProfile(iProf)%nbLoc))
      ALLOCATE(ColumnXProfile(iProf)%ibLoc(ColumnXProfile(iProf)%nbLoc))
      ALLOCATE(ColumnXProfile(iProf)%ProfileBlock(ColumnXProfile(iProf)%nbLoc))
      ColumnXProfile(iProf)%nbLoc=0
    END IF  
  END DO  

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iProf=1,SIZE(ColumnXProfile)
      zCProf=ColumnXProfile(iProf)%P%z
      yCProf=ColumnXProfile(iProf)%P%y
      S2:DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          IF (zP(iz-1)<=zCProf.AND. &
              zP(iz  )>=zCProf.AND. &
              yP(iy-1)<=yCProf.AND. &
              yP(iy  )>=yCProf) THEN
            ColumnXProfile(iProf)%nbLoc=ColumnXProfile(iProf)%nbLoc+1
            ColumnXProfile(iProf)%izProf(ColumnXProfile(iProf)%nbLoc)=iz
            ColumnXProfile(iProf)%iyProf(ColumnXProfile(iProf)%nbLoc)=iy
            ColumnXProfile(iProf)%ibLoc(ColumnXProfile(iProf)%nbLoc)=ibLoc
            ALLOCATE(ColumnXProfile(iProf)%ProfileBlock(ColumnXProfile(iProf)%nbLoc)%Profile(igx0+1:igx1,&
            &SIZE(ColumnXProfile(iProf)%TableC)))
            EXIT S2  
          END IF   
        END DO
      END DO S2
    END DO
  END DO

  IF (MyId==0) THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iProf=1,SIZE(ColumnXProfile)
        zCProf=ColumnXProfile(iProf)%P%z
        yCProf=ColumnXProfile(iProf)%P%y
        IF (z0<=zCProf.AND. &
          z1>=zCProf.AND. &
          y0<=yCProf.AND. &
          y1>=yCProf) THEN
          ColumnXProfileWrite(iProf)%nb=ColumnXProfileWrite(iProf)%nb+1
        END IF    
      END DO
    END DO
    DO iProf=1,SIZE(ColumnXProfile)
      ALLOCATE(ColumnXProfileWrite(iProf)%ib(ColumnXProfileWrite(iProf)%nb))
      ALLOCATE(ColumnXProfileWrite(iProf)%iProc(ColumnXProfileWrite(iProf)%nb))
      ColumnXProfileWrite(iProf)%nb=0
    END DO
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iProf=1,SIZE(ColumnXProfile)
        zCProf=ColumnXProfile(iProf)%P%z
        yCProf=ColumnXProfile(iProf)%P%y
        IF (z0<=zCProf.AND. &
          z1>=zCProf.AND. &
          y0<=yCProf.AND. &
          y1>=yCProf) THEN
          ColumnXProfileWrite(iProf)%nb=ColumnXProfileWrite(iProf)%nb+1
          ColumnXProfileWrite(iProf)%ib(ColumnXProfileWrite(iProf)%nb)=ib
          ColumnXProfileWrite(iProf)%iProc(ColumnXProfileWrite(iProf)%nb)=blMPI(ib)%Proc
          IF (.NOT.ALLOCATED(ColumnXProfileWrite(iProf)%Profile)) THEN
            ALLOCATE(ColumnXProfileWrite(iProf)%Profile(domain%ix0+1:domain%ix1,SIZE(ColumnXProfile(iProf)%TableC)))
          END IF
        END IF    
      END DO
    END DO

    DO iProf=1,SIZE(ColumnXProfileWrite)
      iPz=ColumnXProfileWrite(iProf)%P%z
      iPy=ColumnXProfileWrite(iProf)%P%y
      WRITE(NamePz,'(I10)') iPz
      WRITE(NamePy,'(I10)') iPy
      NamePoint='Py'//TRIM(ADJUSTL(NamePy))//'Pz'//TRIM(ADJUSTL(NamePz))
      ColumnXProfileWrite(iProf)%FileName=TRIM(FileName)//TRIM(ADJUSTL(NamePoint))//'.nc'
      CALL check(nf90_create(TRIM(FileName)//TRIM(ADJUSTL(NamePoint))//'.nc',nf90_clobber,ncid))           ! create netCDF dataset: enter define mode
 
      TableC=>ColumnXProfileWrite(iProf)%TableC
 
      CALL check(nf90_def_dim(ncid,TRIM('xM'),Domain%nx,x_dimid))                ! define dimensions: from name and length=30
      CALL check(nf90_def_dim(ncid,TRIM('time'),NF90_UNLIMITED,rec_dimid))     ! time
      CALL check(nf90_def_dim(ncid,TRIM('Py'),1,py_dimid))
      CALL check(nf90_def_dim(ncid,TRIM('Pz'),1,pz_dimid))

      CALL check(nf90_def_var(ncid,TRIM('xM'),NF90_DOUBLE,x_dimid,x_varid))      ! define variables: from name, type, dims
      CALL check(nf90_def_var(ncid,TRIM('time'),NF90_DOUBLE,rec_dimid,rec_varid))
      CALL check(nf90_def_var(ncid,TRIM('Py'),NF90_INT,py_dimid,py_varid))
      CALL check(nf90_def_var(ncid,TRIM('Pz'),NF90_INT,pz_dimid,pz_varid))
  
      CALL check(nf90_put_att(ncid,x_varid,UNITS,Z_UNITS))                 ! assign attribute values, a numeric rank 1 array
      CALL check(nf90_put_att(ncid,x_varid,'long_name',Z_LONGNAME))
      CALL check(nf90_put_att(ncid,x_varid,'positive','up'))
      CALL check(nf90_put_att(ncid,rec_varid,'long_name',REC_NAME))
      CALL check(nf90_put_att(ncid,rec_varid, UNITS, REC_UNITS))
      CALL check(nf90_put_att(ncid,rec_varid,'calender','noleap'))
      CALL check(nf90_put_att(ncid,pz_varid,UNITS,'m'))
      CALL check(nf90_put_att(ncid,py_varid,UNITS,'m'))
      !CALL check(nf90_put_att(ncid,pz_varid,'z',TRIM(ADJUSTL(NamePz))))
      !CALL check(nf90_put_att(ncid,py_varid,'y',TRIM(ADJUSTL(NamePy))))

      dimids=(/x_dimid,rec_dimid/)
      
      ALLOCATE(ColumnXProfileWrite(iProf)%Scalar_varid(SIZE(TableC),1))
      Scalar_VarId=>ColumnXProfileWrite(iProf)%Scalar_VarId
 
      DO i=1,SIZE(TableC)
        IF (TableC(i)%ColumnXProfile) THEN 
          CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'ColumnXProfile',NF90_DOUBLE,dimids,Scalar_varid(i,1)))
          CALL check(nf90_put_att(ncid,Scalar_varid(i,1),UNITS,TRIM(TableC(i)%Units)))
          CALL check(nf90_put_att(ncid,Scalar_varid(i,1),'Description',TRIM(TableC(i)%Description)))
        END IF  
      END DO  

      CALL check(nf90_enddef(ncid))                                        ! end definitions: leave define mode

      DO ix=Domain%ix0+1,Domain%ix1
        xM(ix)=0.5d0*(Domain%xP(ix-1)+Domain%xP(ix))
      END DO  
      CALL check(nf90_put_var(ncid,x_varid,xM))
      CALL check(nf90_put_var(ncid,pz_varid,iPz))
      CALL check(nf90_put_var(ncid,py_varid,iPy))
      CALL check(nf90_close(ncid))                                        ! close: save new netCDF dataset


      ColumnXProfileWrite(iProf)%ncid=ncid
      ColumnXProfileWrite(iProf)%x_dimid=x_dimid
      ColumnXProfileWrite(iProf)%rec_dimid=rec_dimid
      ColumnXProfileWrite(iProf)%x_varid=x_varid
      ColumnXProfileWrite(iProf)%rec_varid=rec_varid
    END DO     
  END IF
END SUBROUTINE OpenColumnXProfile

SUBROUTINE OutputColumnXProfile(VecC,Time)   
  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind) :: Time

  TYPE(ScalarCell_T) :: ScalarC(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  INTEGER :: ScalarPos
  INTEGER :: iProf
  INTEGER :: ncid
  TYPE(ColumnXProfileTable_T), POINTER :: TableC(:)
  INTEGER, POINTER :: Scalar_VarId(:,:)
  INTEGER :: rec_varid
  REAL(RealKind) :: SaeuleProfile(Domain%nx,MaxLenTable)
  REAL(RealKind) :: SaeuleProfile0(Domain%nx,MaxLenTable)
  INTEGER :: iz,iy,ix,igx
  INTEGER :: ibLoc1,ib1,ibW
  INTEGER :: iProc,iProc1
  TYPE(MPI_Status) :: Status

  LOGICAL, SAVE :: Init=.TRUE.
 
  IF (Init) THEN
    CALL OpenColumnXProfile(NamePointFile)
    Init=.FALSE.
    TimeStep=0
    IF (StartTime>OutputTimeColumnXStart) OutputTimeColumnXStart=StartTime
    OutputTimeColumnX=OutputTimeColumnXStart
  END IF  
  
  IF (OutputTimeColumnX<=Time.AND.OutputTimeColumnX<=OutputTimeColumnXEnd) THEN 
    OutputTimeColumnX=OutputTimeColumnX+dtimeColumnX 
    TimeStep=TimeStep+1
    DO iProf=1,SIZE(ColumnXProfile)
      TableC=>ColumnXProfile(iProf)%TableC
      IF (MyId==0) THEN
        ncid=ColumnXProfileWrite(iProf)%ncid
        rec_varid=ColumnXProfileWrite(iProf)%rec_varid
        Scalar_varid=>ColumnXProfileWrite(iProf)%Scalar_varid
        CALL check(nf90_open(ColumnXProfileWrite(iProf)%FileName,nf90_write,ncid))
        DO i=1,SIZE(TableC)
          IF (TableC(i)%ColumnXProfile) THEN 
            CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'ColumnXProfile',Scalar_varid(i,1)))
          END IF  
        END DO
      END IF  
   
      DO i=1,SIZE(TableC)
        ScalarPos=Position(TRIM(TableC(i)%ScalarC))
        IF (ScalarPos>0) THEN
          CALL Assign(ScalarC,VecC,ScalarPos)
        ELSE IF (TRIM(TableC(i)%ScalarC)=='U') THEN  
          CALL Assign(ScalarC,uCell)
        ! CALL Assign(ScalarC,VecC,uPosl)
        ELSE IF (TRIM(TableC(i)%ScalarC)=='V') THEN
          CALL Assign(ScalarC,vCell)
        ! CALL Assign(ScalarC,VecC,vPosl)
        ELSE IF (TRIM(TableC(i)%ScalarC)=='W') THEN
          CALL Assign(ScalarC,wCell)
        ! CALL Assign(ScalarC,VecC,wPosl)
        END IF  
        IF (TableC(i)%ColumnXProfile) THEN 
          IF (RhoPos>0) THEN
            CALL Assign(Rho,VecC,RhoPos)
            DO ibLoc1=1,ColumnXProfile(iProf)%nbLoc
              ibLoc=ColumnXProfile(iProf)%ibLoc(ibLoc1)
              ib=LocGlob(ibLoc)
              CALL Set(Floor(ib))
              iz=ColumnXProfile(iProf)%izProf(ibLoc1)
              iy=ColumnXProfile(iProf)%iyProf(ibLoc1)
              DO ix=ix0+1,ix1
                DO igx=(ix-1)*2**(-RefineZ)+1,ix*2**(-RefineZ) 
                  ColumnXProfile(iProf)%ProfileBlock(ibLoc1)%Profile(igx,i)=ScalarC(ibLoc)%c(ix,iy,iz,1)&
                  &/(Rho(ibLoc)%c(ix,iy,iz,1)+Eps)
                END DO  
              END DO  
            END DO  
          ELSE   
            DO ibLoc1=1,ColumnXProfile(iProf)%nbLoc
              ibLoc=ColumnXProfile(iProf)%ibLoc(ibLoc1)
              ib=LocGlob(ibLoc)
              CALL Set(Floor(ib))
              iz=ColumnXProfile(iProf)%izProf(ibLoc1)
              iy=ColumnXProfile(iProf)%iyProf(ibLoc1)
              DO ix=ix0+1,ix1
                DO igx=(ix-1)*2**(-RefineZ)+1,ix*2**(-RefineZ) 
                  ColumnXProfile(iProf)%ProfileBlock(ibLoc1)%Profile(igx,i)=ScalarC(ibLoc)%c(ix,iy,iz,1)
                END DO  
              END DO  
            END DO  
          END IF
        END IF
      END DO  
      IF (MyId==0) THEN
        DO ib1=1,ColumnXProfileWrite(iProf)%nb
          ib=ColumnXProfileWrite(iProf)%ib(ib1)
          iProc=ColumnXProfileWrite(iProf)%iProc(ib1)
          CALL Set(Floor(ib))
          IF (iProc>0) THEN
            CALL MPI_RECV(ColumnXProfileWrite(iProf)%Profile(igx0+1:igx1,:),nx*SIZE(ColumnXProfile(iProf)%TableC) &
                         ,MPI_RealKind,iProc,iProc+ib*nb,MPI_Comm_World,Status,MPIErr)
          END IF  
        END DO
      END IF  
      IF (ColumnXProfile(iProf)%nbLoc>0) THEN
        DO ibLoc1=1,ColumnXProfile(iProf)%nbLoc
          ibLoc=ColumnXProfile(iProf)%ibLoc(ibLoc1)
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          IF (MyId>0) THEN
            CALL MPI_SEND(ColumnXProfile(iProf)%ProfileBlock(ibLoc1)%Profile,&
            &nx*SIZE(ColumnXProfile(iProf)%TableC),MPI_RealKind,0,MyId+ib*nb,MPI_Comm_World,MPIErr)
          ELSE
            DO ib1=1,ColumnXProfileWrite(iProf)%nb
              ibW=ColumnXProfileWrite(iProf)%ib(ib1)
              IF (ibW==ib) THEN
                ColumnXProfileWrite(iProf)%Profile(igx0+1:igx1,:)=ColumnXProfile(iProf)%ProfileBlock(ibLoc1)%Profile
              END IF  
            END DO  
          END IF  
        END DO  
      END IF  
      IF (MyId==0) THEN
        DO i=1,SIZE(TableC)
          IF (TableC(i)%ColumnXProfile) THEN 
            CALL check(nf90_put_var(ncid,Scalar_varid(i,1),ColumnXProfileWrite(iProf)%Profile(:,i),start=(/1,TimeStep/)))
          END IF  
        END DO
        CALL check(nf90_inq_varid(ncid,TRIM('time'),rec_varid))
        CALL check(nf90_put_var(ncid,rec_varid,Time,start=(/TimeStep/)))
        
        !CALL start_MPI()!!
        !start=(/1,MyID+1/)!!
        !count=(/p,1/)!!
        !CALL check(nf90_put_var(ncid,Scalar_varid(i,1),SaeuleProfile(:,i),start=start,&
         !    count=count))!!
        ! Close file
        CALL check(nf90_sync(ncid))
        CALL check(nf90_close(ncid))
        ! Free local memory.
 !       DEALLOCATE(SaeuleProfile)!!
        ! MPI library must be shut down.
        !CALL MPI_Finalixe(MPIerr)!!
      END IF
    END DO
  END IF
END SUBROUTINE OutputColumnXProfile

SUBROUTINE CloseColumnXProfile
  INTEGER :: iProf
  INTEGER :: ncid
  TYPE(ColumnXProfileTable_T), POINTER :: TableC(:)
  
  IF (MyId==0) THEN
    DO iProf=1,SIZE(ColumnXProfileWrite)
      ncid=ColumnXProfileWrite(iProf)%ncid
      CALL check(nf90_close(ncid)) ! close: save new netCDF dataset
    END DO 
  END IF
END SUBROUTINE CloseColumnXProfile

SUBROUTINE ColumnXProfileNMLOutput(FileName,Check)
  CHARACTER(*) :: FileName
  LOGICAL :: Check
  CHARACTER(600) :: Line
 
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ColumnXProfileOutput')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ColumnXProfileOutput)
      Check=.TRUE.
      EXIT
    END IF
  END DO
1 CONTINUE 
  CLOSE(InputUnit)
END SUBROUTINE ColumnXProfileNMLOutput

SUBROUTINE ColumnXProfileDataTable(FileName)
  CHARACTER(*) :: FileName

  CHARACTER*200 :: Line
  INTEGER :: i,j,LenTable,NPoints 

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*) Line
    IF (INDEX(Line,'#ColumnXProfileTable')>0) THEN
      READ(InputUnit,*) NPoints 
      ALLOCATE(ColumnXProfile(NPoints))
      DO i=1,NPoints
        READ(InputUnit,*) Line
        READ(InputUnit,*) ColumnXProfile(i)%P%y,ColumnXProfile(i)%P%z
        READ(InputUnit,*) LenTable 
        MaxLenTable=MAX(MaxLenTable,LenTable)
        ALLOCATE(ColumnXProfile(i)%TableC(LenTable))
        READ(InputUnit,*) Line 
        DO j=1,LenTable
          READ(InputUnit,*) ColumnXProfile(i)%TableC(j)%ScalarC,ColumnXProfile(i)%TableC(j)%Units,&
          &ColumnXProfile(i)%TableC(j)%Description &
                           ,ColumnXProfile(i)%TableC(j)%ColumnXProfile
        END DO  
        READ(InputUnit,*) Line
      END DO  
      EXIT
    END IF
  END DO
  CLOSE(InputUnit) 
  IF (MyId==0) THEN
    ALLOCATE(ColumnXProfileWrite(NPoints))
    DO i=1,NPoints
      ColumnXProfileWrite(i)%P=ColumnXProfile(i)%P
      ALLOCATE(ColumnXProfileWrite(i)%TableC(LenTable))
      DO j=1,LenTable
        ColumnXProfileWrite(i)%TableC(j)=ColumnXProfile(i)%TableC(j)
      END DO  
    END DO  
  END IF  
END SUBROUTINE ColumnXProfileDataTable
END MODULE SingleColumnXProfile_Mod
