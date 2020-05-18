MODULE SingleColumnZProfile_Mod

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
  REAL(RealKind), PRIVATE :: OutputTimeColumnZ
  REAL(RealKind) :: dtimeColumnZ=0.0d0
  REAL(RealKind) :: OutputTimeColumnZStart=1.0d40
  REAL(RealKind) :: OutputTimeColumnZEnd=-1.0d40

  TYPE ColumnZProfileTable_T
    INTEGER :: ibLoc
    CHARACTER(32) :: ScalarC=' '
    LOGICAL       :: ColumnZProfile=.FALSE.
    CHARACTER(20) :: Units=''
    CHARACTER(30) :: Description=''
    INTEGER       :: iPosC=0
    LOGICAL       :: ScaledRho
  END TYPE ColumnZProfileTable_T
  TYPE ColumnZPoint_T
    REAL(RealKind) :: x
    REAL(RealKind) :: y
    REAL(RealKind) :: z
  END TYPE ColumnZPoint_T
  TYPE Profile_T
    REAL(RealKind), ALLOCATABLE :: Profile(:,:)
  END TYPE Profile_T  
  TYPE ColumnZProfile_T
    TYPE(ColumnZProfileTable_T), POINTER :: TableC(:)
    TYPE(ColumnZPoint_T) :: P
    TYPE(Profile_T), ALLOCATABLE :: ProfileBlock(:)
    INTEGER :: nbLoc=0
    INTEGER, ALLOCATABLE :: ixProf(:)
    INTEGER, ALLOCATABLE :: iyProf(:)
    INTEGER, ALLOCATABLE :: ibLoc(:)
    LOGICAL, ALLOCATABLE :: OutBlock(:)
    LOGICAL :: OutHeader=.FALSE.
    CHARACTER*50 :: FileName
!   INTEGER :: ncid
!   INTEGER :: z_varid, rec_varid, px_varid, py_dimid
!   INTEGER :: z_dimid, rec_dimid, px_dimid, py_varid
!   INTEGER, POINTER :: Scalar_VarId(:,:)
  END TYPE  
  TYPE ColumnZProfileWrite_T
    INTEGER :: ibLoc
    TYPE(ColumnZProfileTable_T), POINTER :: TableC(:)
    TYPE(ColumnZPoint_T) :: P
    REAL(RealKind), ALLOCATABLE :: Profile(:,:)
    INTEGER :: nb=0
    INTEGER :: nProc
    INTEGER, ALLOCATABLE :: ixProf(:)
    INTEGER, ALLOCATABLE :: iyProf(:)
    INTEGER, ALLOCATABLE :: ib(:)
    INTEGER, ALLOCATABLE :: iProc(:)
    LOGICAL, ALLOCATABLE :: OutBlock(:)
    LOGICAL :: OutHeader=.FALSE.
    CHARACTER*50 :: FileName
    INTEGER :: ncid
    INTEGER :: z_varid, rec_varid, px_varid, py_dimid
    INTEGER :: z_dimid, rec_dimid, px_dimid, py_varid
    INTEGER, POINTER :: Scalar_VarId(:,:)
  END TYPE  
  TYPE(ColumnZProfile_T), ALLOCATABLE :: ColumnZProfile(:)
  TYPE(ColumnZProfileWrite_T), ALLOCATABLE :: ColumnZProfileWrite(:)

  NAMELIST /ColumnZProfileOutput/ OutputTimeColumnZStart, &
                                 OutputTimeColumnZEnd, &
                                 dtimeColumnZ, &
                                 NamePointFile
CONTAINS

SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS

  IF(STATUS /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(STATUS))
  END IF
END SUBROUTINE check

SUBROUTINE OpenColumnZProfile(FileName)
  CHARACTER(*) :: FileName

  REAL(RealKind) :: zM(Domain%nz)
  INTEGER :: ix,iy,iz,ixProf,iyProf,izProf
  INTEGER :: iProf,ibProfLoc,ibProf
  TYPE(ColumnZProfileTable_T), POINTER :: TableC(:)
  INTEGER, POINTER :: Scalar_VarId(:,:)
  INTEGER :: iPx,iPy
  CHARACTER*12 :: NamePx,NamePy
  INTEGER :: ncid
  INTEGER :: z_varid, rec_varid, px_varid, py_varid
  INTEGER :: z_dimid, rec_dimid, px_dimid, py_dimid

! Determine number of local blocks
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iProf=1,UBOUND(ColumnZProfile,1)
      xProf=ColumnZProfile(iProf)%P%x
      yProf=ColumnZProfile(iProf)%P%y
      S1:DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          IF (xP(ix-1)<=xProf.AND. &
              xP(ix  )>=xProf.AND. &
              yP(iy-1)<=yProf.AND. &
              yP(iy  )>=yProf) THEN
            ColumnZProfile(iProf)%nbLoc=ColumnZProfile(iProf)%nbLoc+1
            EXIT S1  
          END IF
        END DO 
      END DO S1
    END DO
  END DO
  DO iProf=1,UBOUND(ColumnZProfile,1)
    IF (ColumnZProfile(iProf)%nbLoc>0) THEN
      ALLOCATE(ColumnZProfile(iProf)%ixProf(ColumnZProfile(iProf)%nbLoc))
      ALLOCATE(ColumnZProfile(iProf)%iyProf(ColumnZProfile(iProf)%nbLoc))
      ALLOCATE(ColumnZProfile(iProf)%ibLoc(ColumnZProfile(iProf)%nbLoc))
      ALLOCATE(ColumnZProfile(iProf)%ProfileBlock(ColumnZProfile(iProf)%nbLoc))
      ColumnZProfile(iProf)%nbLoc=0
    END IF  
  END DO  

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iProf=1,UBOUND(ColumnZProfile,1)
      xProf=ColumnZProfile(iProf)%P%x
      yProf=ColumnZProfile(iProf)%P%y
      S2:DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          IF (xP(ix-1)<=xProf.AND. &
              xP(ix  )>=xProf.AND. &
              yP(iy-1)<=yProf.AND. &
              yP(iy  )>=yProf) THEN
            ColumnZProfile(iProf)%nbLoc=ColumnZProfile(iProf)%nbLoc+1
            ColumnZProfile(iProf)%ixProf(ColumnZProfile(iProf)%nbLoc)=ix
            ColumnZProfile(iProf)%iyProf(ColumnZProfile(iProf)%nbLoc)=iy
            ColumnZProfile(iProf)%ibLoc(ColumnZProfile(iProf)%nbLoc)=ibLoc
            ALLOCATE(ColumnZProfile(iProf)%ProfileBlock(ColumnZProfile(iProf)%nbLoc)%Profile(igz0+1:igz1,&
            &SIZE(ColumnZProfile(iProf)%TableC)))
            EXIT S2  
          END IF   
        END DO
      END DO S2
    END DO
  END DO

  IF (MyId==0) THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iProf=1,UBOUND(ColumnZProfile,1)
        xProf=ColumnZProfile(iProf)%P%x
        yProf=ColumnZProfile(iProf)%P%y
        IF (x0<=xProf.AND. &
          x1>=xProf.AND. &
          y0<=yProf.AND. &
          y1>=yProf) THEN
          ColumnZProfileWrite(iProf)%nb=ColumnZProfileWrite(iProf)%nb+1
        END IF    
      END DO
    END DO
    DO iProf=1,SIZE(ColumnZProfile)
      ALLOCATE(ColumnZProfileWrite(iProf)%ib(ColumnZProfileWrite(iProf)%nb))
      ALLOCATE(ColumnZProfileWrite(iProf)%iProc(ColumnZProfileWrite(iProf)%nb))
      ColumnZProfileWrite(iProf)%nb=0
    END DO
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iProf=1,UBOUND(ColumnZProfile,1)
        xProf=ColumnZProfile(iProf)%P%x
        yProf=ColumnZProfile(iProf)%P%y
        IF (x0<=xProf.AND. &
          x1>=xProf.AND. &
          y0<=yProf.AND. &
          y1>=yProf) THEN
          ColumnZProfileWrite(iProf)%nb=ColumnZProfileWrite(iProf)%nb+1
          ColumnZProfileWrite(iProf)%ib(ColumnZProfileWrite(iProf)%nb)=ib
          ColumnZProfileWrite(iProf)%iProc(ColumnZProfileWrite(iProf)%nb)=blMPI(ib)%Proc
          IF (.NOT.ALLOCATED(ColumnZProfileWrite(iProf)%Profile)) THEN
            ALLOCATE(ColumnZProfileWrite(iProf)%Profile(domain%iz0+1:domain%iz1,SIZE(ColumnZProfile(iProf)%TableC)))
          END IF
        END IF    
      END DO
    END DO

    DO iProf=1,UBOUND(ColumnZProfile,1)
      iPx=ColumnZProfileWrite(iProf)%P%x
      iPy=ColumnZProfileWrite(iProf)%P%y
      WRITE(NamePx,'(I10)') iPx
      WRITE(NamePy,'(I10)') iPy
      NamePoint='Px'//TRIM(ADJUSTL(NamePx))//'Py'//TRIM(ADJUSTL(NamePy))
      ColumnZProfileWrite(iProf)%FileName=TRIM(FileName)//TRIM(ADJUSTL(NamePoint))//'.nc'
      CALL check(nf90_create(TRIM(FileName)//TRIM(ADJUSTL(NamePoint))//'.nc',nf90_clobber,ncid))           ! create netCDF dataset: enter define mode
 
      TableC=>ColumnZProfileWrite(iProf)%TableC
 
      CALL check(nf90_def_dim(ncid,TRIM('zM'),Domain%nz,z_dimid))                ! define dimensions: from name and length=30
      CALL check(nf90_def_dim(ncid,TRIM('time'),NF90_UNLIMITED,rec_dimid))     ! time
      CALL check(nf90_def_dim(ncid,TRIM('Px'),1,px_dimid))
      CALL check(nf90_def_dim(ncid,TRIM('Py'),1,py_dimid))

      CALL check(nf90_def_var(ncid,TRIM('zM'),NF90_DOUBLE,z_dimid,z_varid))      ! define variables: from name, type, dims
      CALL check(nf90_def_var(ncid,TRIM('time'),NF90_DOUBLE,rec_dimid,rec_varid))
      CALL check(nf90_def_var(ncid,TRIM('Px'),NF90_INT,px_dimid,px_varid))
      CALL check(nf90_def_var(ncid,TRIM('Py'),NF90_INT,py_dimid,py_varid))
  
      CALL check(nf90_put_att(ncid,z_varid,UNITS,Z_UNITS))                 ! assign attribute values, a numeric rank 1 array
      CALL check(nf90_put_att(ncid,z_varid,'long_name',Z_LONGNAME))
      CALL check(nf90_put_att(ncid,z_varid,'positive','up'))
      CALL check(nf90_put_att(ncid,rec_varid,'long_name',REC_NAME))
      CALL check(nf90_put_att(ncid,rec_varid, UNITS, REC_UNITS))
      CALL check(nf90_put_att(ncid,rec_varid,'calender','noleap'))
      CALL check(nf90_put_att(ncid,px_varid,UNITS,'m'))
      CALL check(nf90_put_att(ncid,py_varid,UNITS,'m'))
      !CALL check(nf90_put_att(ncid,px_varid,'x',TRIM(ADJUSTL(NamePx))))
      !CALL check(nf90_put_att(ncid,py_varid,'y',TRIM(ADJUSTL(NamePy))))

      dimids=(/z_dimid,rec_dimid/)
      
      ALLOCATE(ColumnZProfileWrite(iProf)%Scalar_varid(SIZE(TableC),1))
      Scalar_VarId=>ColumnZProfileWrite(iProf)%Scalar_VarId
 
      DO i=1,SIZE(TableC)
        IF (TableC(i)%ColumnZProfile) THEN 
          CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'ColumnZProfile',NF90_DOUBLE,dimids,Scalar_varid(i,1)))
          CALL check(nf90_put_att(ncid,Scalar_varid(i,1),UNITS,TRIM(TableC(i)%Units)))
          CALL check(nf90_put_att(ncid,Scalar_varid(i,1),'Description',TRIM(TableC(i)%Description)))
        END IF  
      END DO  

      CALL check(nf90_enddef(ncid))                                        ! end definitions: leave define mode

      DO iz=Domain%iz0+1,Domain%iz1
        zM(iz)=0.5d0*(Domain%zP(iz-1)+Domain%zP(iz))
      END DO  
      CALL check(nf90_put_var(ncid,z_varid,zM))
      CALL check(nf90_put_var(ncid,px_varid,iPx))
      CALL check(nf90_put_var(ncid,py_varid,iPy))
      CALL check(nf90_close(ncid))                                        ! close: save new netCDF dataset


      ColumnZProfileWrite(iProf)%ncid=ncid
      ColumnZProfileWrite(iProf)%z_dimid=z_dimid
      ColumnZProfileWrite(iProf)%rec_dimid=rec_dimid
      ColumnZProfileWrite(iProf)%z_varid=z_varid
      ColumnZProfileWrite(iProf)%rec_varid=rec_varid
    END DO     
  END IF
END SUBROUTINE OpenColumnZProfile

SUBROUTINE OutputColumnZProfile(VecC,Time)   
  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind) :: Time

  TYPE(ScalarCell_T) :: ScalarC(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  INTEGER :: ScalarPos
  INTEGER :: iProf
  INTEGER :: ncid
  TYPE(ColumnZProfileTable_T), POINTER :: TableC(:)
  INTEGER, POINTER :: Scalar_VarId(:,:)
  INTEGER :: rec_varid
  REAL(RealKind) :: SaeuleProfile(Domain%nz,MaxLenTable)
  REAL(RealKind) :: SaeuleProfile0(Domain%nz,MaxLenTable)
  INTEGER :: ix,iy,iz,igz
  INTEGER :: ibLoc1,ib1,ibW
  INTEGER :: iProc,iProc1
  TYPE(MPI_Status) :: Status

  LOGICAL, SAVE :: Init=.TRUE.
 
  IF (Init) THEN
    CALL OpenColumnZProfile(NamePointFile)
    Init=.FALSE.
    TimeStep=0
    IF (StartTime>OutputTimeColumnZStart) OutputTimeColumnZStart=StartTime
    OutputTimeColumnZ=OutputTimeColumnZStart
  END IF  
  
  IF (OutputTimeColumnZ<=Time.AND.OutputTimeColumnZ<=OutputTimeColumnZEnd) THEN 
    OutputTimeColumnZ=OutputTimeColumnZ+dtimeColumnZ 
    TimeStep=TimeStep+1
    DO iProf=1,SIZE(ColumnZProfile)
      TableC=>ColumnZProfile(iProf)%TableC
      IF (MyId==0) THEN
        ncid=ColumnZProfileWrite(iProf)%ncid
        rec_varid=ColumnZProfileWrite(iProf)%rec_varid
        Scalar_varid=>ColumnZProfileWrite(iProf)%Scalar_varid
        CALL check(nf90_open(ColumnZProfileWrite(iProf)%FileName,nf90_write,ncid))
        DO i=1,SIZE(TableC)
          IF (TableC(i)%ColumnZProfile) THEN 
            CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'ColumnZProfile',Scalar_varid(i,1)))
          END IF  
        END DO
      END IF  
   
      DO i=1,SIZE(TableC)
        ScalarPos=Position(TRIM(TableC(i)%ScalarC))
        IF (ScalarPos>0) THEN
          CALL Assign(ScalarC,VecC,ScalarPos)
        ELSE IF (TRIM(TableC(i)%ScalarC)=='U') THEN  
          CALL Assign(ScalarC,uCell)
!         CALL Assign(ScalarC,VecC,uPosl)
        ELSE IF (TRIM(TableC(i)%ScalarC)=='V') THEN
          CALL Assign(ScalarC,vCell)
!         CALL Assign(ScalarC,VecC,vPosl)
        ELSE IF (TRIM(TableC(i)%ScalarC)=='W') THEN
          CALL Assign(ScalarC,wCell)
!         CALL Assign(ScalarC,VecC,wPosl)
        END IF  
        IF (TableC(i)%ColumnZProfile) THEN 
          IF (RhoPos>0.AND.TableC(i)%ScaledRho) THEN
            CALL Assign(Rho,VecC,RhoPos)
            DO ibLoc1=1,ColumnZProfile(iProf)%nbLoc
              ibLoc=ColumnZProfile(iProf)%ibLoc(ibLoc1)
              ib=LocGlob(ibLoc)
              CALL Set(Floor(ib))
              ix=ColumnZProfile(iProf)%ixProf(ibLoc1)
              iy=ColumnZProfile(iProf)%iyProf(ibLoc1)
              DO iz=iz0+1,iz1
                DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ) 
                  ColumnZProfile(iProf)%ProfileBlock(ibLoc1)%Profile(igz,i)=ScalarC(ibLoc)%c(ix,iy,iz,1)&
                  &/(Rho(ibLoc)%c(ix,iy,iz,1)+Eps)
                END DO  
              END DO  
            END DO  
          ELSE   
            DO ibLoc1=1,ColumnZProfile(iProf)%nbLoc
              ibLoc=ColumnZProfile(iProf)%ibLoc(ibLoc1)
              ib=LocGlob(ibLoc)
              CALL Set(Floor(ib))
              ix=ColumnZProfile(iProf)%ixProf(ibLoc1)
              iy=ColumnZProfile(iProf)%iyProf(ibLoc1)
              DO iz=iz0+1,iz1
                DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ) 
                  ColumnZProfile(iProf)%ProfileBlock(ibLoc1)%Profile(igz,i)=ScalarC(ibLoc)%c(ix,iy,iz,1)
                END DO  
              END DO  
            END DO  
          END IF
        END IF
      END DO  
      IF (MyId==0) THEN
        DO ib1=1,ColumnZProfileWrite(iProf)%nb
          ib=ColumnZProfileWrite(iProf)%ib(ib1)
          iProc=ColumnZProfileWrite(iProf)%iProc(ib1)
          CALL Set(Floor(ib))
          IF (iProc>0) THEN
            CALL MPI_RECV(ColumnZProfileWrite(iProf)%Profile(igz0+1:igz1,:),nz*SIZE(ColumnZProfile(iProf)%TableC) &
                         ,MPI_RealKind,iProc,iProc+ib*nb,MPI_Comm_World,Status,MPIErr)
          END IF  
        END DO
      END IF  
      IF (ColumnZProfile(iProf)%nbLoc>0) THEN
        DO ibLoc1=1,ColumnZProfile(iProf)%nbLoc
          ibLoc=ColumnZProfile(iProf)%ibLoc(ibLoc1)
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          IF (MyId>0) THEN
            CALL MPI_SEND(ColumnZProfile(iProf)%ProfileBlock(ibLoc1)%Profile,nz*SIZE(ColumnZProfile(iProf)%TableC),&
            &MPI_RealKind,0,MyId+ib*nb,MPI_Comm_World,MPIErr)
          ELSE
            DO ib1=1,ColumnZProfileWrite(iProf)%nb
              ibW=ColumnZProfileWrite(iProf)%ib(ib1)
              IF (ibW==ib) THEN
                ColumnZProfileWrite(iProf)%Profile(igz0+1:igz1,:)=ColumnZProfile(iProf)%ProfileBlock(ibLoc1)%Profile
              END IF  
            END DO  
          END IF  
        END DO  
      END IF  
      IF (MyId==0) THEN
        DO i=1,SIZE(TableC)
          IF (TableC(i)%ColumnZProfile) THEN 
            CALL check(nf90_put_var(ncid,Scalar_varid(i,1),ColumnZProfileWrite(iProf)%Profile(:,i),start=(/1,TimeStep/)))
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
        !CALL MPI_Finalize(MPIerr)!!
      END IF
    END DO
  END IF
END SUBROUTINE OutputColumnZProfile

SUBROUTINE CloseColumnZProfile
  INTEGER :: iProf
  INTEGER :: ncid
  TYPE(ColumnZProfileTable_T), POINTER :: TableC(:)
  
  IF (MyId==0) THEN
    DO iProf=1,SIZE(ColumnZProfileWrite)
      ncid=ColumnZProfileWrite(iProf)%ncid
      CALL check(nf90_close(ncid)) ! close: save new netCDF dataset
    END DO 
  END IF
END SUBROUTINE CloseColumnZProfile

SUBROUTINE ColumnZProfileNMLOutput(FileName,Check)
  CHARACTER(*) :: FileName
  LOGICAL :: Check
  CHARACTER(600) :: Line
 
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ColumnZProfileOutput')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ColumnZProfileOutput)
      Check=.TRUE.
      EXIT
    END IF
  END DO
1 CONTINUE 
  CLOSE(InputUnit)
END SUBROUTINE ColumnZProfileNMLOutput

SUBROUTINE ColumnZProfileDataTable(FileName)
  CHARACTER(*) :: FileName

  CHARACTER*200 :: Line
  INTEGER :: i,j,LenTable,NPoints 

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*) Line
    IF (INDEX(Line,'#ColumnZProfileTable')>0) THEN
      READ(InputUnit,*) NPoints 
      ALLOCATE(ColumnZProfile(NPoints))
      DO i=1,NPoints
        READ(InputUnit,*) Line
        READ(InputUnit,*) ColumnZProfile(i)%P%x,ColumnZProfile(i)%P%y
        READ(InputUnit,*) LenTable 
        MaxLenTable=MAX(MaxLenTable,LenTable)
        ALLOCATE(ColumnZProfile(i)%TableC(LenTable))
        READ(InputUnit,*) Line 
        DO j=1,LenTable
          READ(InputUnit,*) ColumnZProfile(i)%TableC(j)%ScalarC,ColumnZProfile(i)%TableC(j)%Units,&
          &ColumnZProfile(i)%TableC(j)%Description &
                           ,ColumnZProfile(i)%TableC(j)%ColumnZProfile,ColumnZProfile(i)%TableC(j)%ScaledRho
        END DO  
        READ(InputUnit,*) Line
      END DO  
      EXIT
    END IF
  END DO
  CLOSE(InputUnit) 
  IF (MyId==0) THEN
    ALLOCATE(ColumnZProfileWrite(NPoints))
    DO i=1,NPoints
      ColumnZProfileWrite%P=ColumnZProfile%P
      ALLOCATE(ColumnZProfileWrite(i)%TableC(LenTable))
      DO j=1,LenTable
        ColumnZProfileWrite(i)%TableC(j)=ColumnZProfile(i)%TableC(j)
      END DO  
    END DO  
  END IF  
END SUBROUTINE ColumnZProfileDataTable
END MODULE SingleColumnZProfile_Mod
