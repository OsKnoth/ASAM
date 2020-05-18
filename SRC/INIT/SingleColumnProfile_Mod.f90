MODULE SingleColumnProfile_Mod

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
  REAL(RealKind), PRIVATE :: OutputTimeColumn
  REAL(RealKind) :: dtimeColumn
  REAL(RealKind) :: OutputTimeColumnStart
  REAL(RealKind) :: OutputTimeColumnEnd

  TYPE ColumnProfileTable_T
    INTEGER :: ibLoc
    CHARACTER(32) :: ScalarC=' '
    LOGICAL       :: ColumnProfile=.FALSE.
    CHARACTER(20) :: Units=''
    CHARACTER(30) :: Description=''
    INTEGER       :: iPosC=0
  END TYPE ColumnProfileTable_T
  TYPE ColumnPoint_T
    REAL(RealKind) :: x
    REAL(RealKind) :: y
    REAL(RealKind) :: z
  END TYPE ColumnPoint_T
  TYPE Profile_T
    REAL(RealKind), ALLOCATABLE :: Profile(:,:)
  END TYPE Profile_T  
  TYPE ColumnProfile_T
    TYPE(ColumnProfileTable_T), POINTER :: TableC(:)
    TYPE(ColumnPoint_T) :: P
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
  TYPE ColumnProfileWrite_T
    INTEGER :: ibLoc
    TYPE(ColumnProfileTable_T), POINTER :: TableC(:)
    TYPE(ColumnPoint_T) :: P
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
  TYPE(ColumnProfile_T), ALLOCATABLE :: ColumnProfile(:)
  TYPE(ColumnProfileWrite_T), ALLOCATABLE :: ColumnProfileWrite(:)

  NAMELIST /ColumnProfileOutput/ OutputTimeColumnStart, &
                                 OutputTimeColumnEnd, &
                                 dtimeColumn, &
                                 NamePointFile
CONTAINS

SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS

  IF(STATUS /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(STATUS))
  END IF
END SUBROUTINE check

SUBROUTINE OpenColumnProfile(FileName)
  CHARACTER(*) :: FileName

  REAL(RealKind) :: zM(Domain%nz)
  INTEGER :: ix,iy,iz,ixProf,iyProf,izProf
  INTEGER :: iProf,ibProfLoc,ibProf
  TYPE(ColumnProfileTable_T), POINTER :: TableC(:)
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
    DO iProf=1,SIZE(ColumnProfile)
      xProf=ColumnProfile(iProf)%P%x
      yProf=ColumnProfile(iProf)%P%y
      S1:DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          IF (xP(ix-1)<=xProf.AND. &
              xP(ix  )>=xProf.AND. &
              yP(iy-1)<=yProf.AND. &
              yP(iy  )>=yProf) THEN
            ColumnProfile(iProf)%nbLoc=ColumnProfile(iProf)%nbLoc+1
            EXIT S1  
          END IF
        END DO 
      END DO S1
    END DO
  END DO
  DO iProf=1,SIZE(ColumnProfile)
    IF (ColumnProfile(iProf)%nbLoc>0) THEN
      ALLOCATE(ColumnProfile(iProf)%ixProf(ColumnProfile(iProf)%nbLoc))
      ALLOCATE(ColumnProfile(iProf)%iyProf(ColumnProfile(iProf)%nbLoc))
      ALLOCATE(ColumnProfile(iProf)%ibLoc(ColumnProfile(iProf)%nbLoc))
      ALLOCATE(ColumnProfile(iProf)%ProfileBlock(ColumnProfile(iProf)%nbLoc))
      ColumnProfile(iProf)%nbLoc=0
    END IF  
  END DO  

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iProf=1,SIZE(ColumnProfile)
      xProf=ColumnProfile(iProf)%P%x
      yProf=ColumnProfile(iProf)%P%y
      S2:DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          IF (xP(ix-1)<=xProf.AND. &
              xP(ix  )>=xProf.AND. &
              yP(iy-1)<=yProf.AND. &
              yP(iy  )>=yProf) THEN
            ColumnProfile(iProf)%nbLoc=ColumnProfile(iProf)%nbLoc+1
            ColumnProfile(iProf)%ixProf(ColumnProfile(iProf)%nbLoc)=ix
            ColumnProfile(iProf)%iyProf(ColumnProfile(iProf)%nbLoc)=iy
            ColumnProfile(iProf)%ibLoc(ColumnProfile(iProf)%nbLoc)=ibLoc
            ALLOCATE(ColumnProfile(iProf)%ProfileBlock(ColumnProfile(iProf)%nbLoc)%&
            & Profile(igz0+1:igz1,SIZE(ColumnProfile(iProf)%TableC)))
            EXIT S2  
          END IF   
        END DO
      END DO S2
    END DO
  END DO

  IF (MyId==0) THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iProf=1,SIZE(ColumnProfile)
        xProf=ColumnProfile(iProf)%P%x
        yProf=ColumnProfile(iProf)%P%y
        IF (x0<=xProf.AND. &
          x1>=xProf.AND. &
          y0<=yProf.AND. &
          y1>=yProf) THEN
          ColumnProfileWrite(iProf)%nb=ColumnProfileWrite(iProf)%nb+1
        END IF    
      END DO
    END DO
    DO iProf=1,SIZE(ColumnProfile)
      ALLOCATE(ColumnProfileWrite(iProf)%ib(ColumnProfileWrite(iProf)%nb))
      ALLOCATE(ColumnProfileWrite(iProf)%iProc(ColumnProfileWrite(iProf)%nb))
      ColumnProfileWrite(iProf)%nb=0
    END DO
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO iProf=1,SIZE(ColumnProfile)
        xProf=ColumnProfile(iProf)%P%x
        yProf=ColumnProfile(iProf)%P%y
        IF (x0<=xProf.AND. &
          x1>=xProf.AND. &
          y0<=yProf.AND. &
          y1>=yProf) THEN
          ColumnProfileWrite(iProf)%nb=ColumnProfileWrite(iProf)%nb+1
          ColumnProfileWrite(iProf)%ib(ColumnProfileWrite(iProf)%nb)=ib
          ColumnProfileWrite(iProf)%iProc(ColumnProfileWrite(iProf)%nb)=blMPI(ib)%Proc
          IF (.NOT.ALLOCATED(ColumnProfileWrite(iProf)%Profile)) THEN
            ALLOCATE(ColumnProfileWrite(iProf)%Profile(domain%iz0+1:domain%iz1,SIZE(ColumnProfile(iProf)%TableC)))
          END IF
        END IF    
      END DO
    END DO

    DO iProf=1,SIZE(ColumnProfileWrite)
      iPx=ColumnProfileWrite(iProf)%P%x
      iPy=ColumnProfileWrite(iProf)%P%y
      WRITE(NamePx,'(I10)') iPx
      WRITE(NamePy,'(I10)') iPy
      NamePoint='Px'//TRIM(ADJUSTL(NamePx))//'Py'//TRIM(ADJUSTL(NamePy))
      ColumnProfileWrite(iProf)%FileName=TRIM(FileName)//TRIM(ADJUSTL(NamePoint))//'.nc'
      CALL check(nf90_create(TRIM(FileName)//TRIM(ADJUSTL(NamePoint))//'.nc',nf90_clobber,ncid))           ! create netCDF dataset: enter define mode
 
      TableC=>ColumnProfileWrite(iProf)%TableC
 
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
      
      ALLOCATE(ColumnProfileWrite(iProf)%Scalar_varid(SIZE(TableC),1))
      Scalar_VarId=>ColumnProfileWrite(iProf)%Scalar_VarId
 
      DO i=1,SIZE(TableC)
        IF (TableC(i)%ColumnProfile) THEN 
          CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'ColumnProfile',NF90_DOUBLE,dimids,Scalar_varid(i,1)))
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


      ColumnProfileWrite(iProf)%ncid=ncid
      ColumnProfileWrite(iProf)%z_dimid=z_dimid
      ColumnProfileWrite(iProf)%rec_dimid=rec_dimid
      ColumnProfileWrite(iProf)%z_varid=z_varid
      ColumnProfileWrite(iProf)%rec_varid=rec_varid
    END DO     
  END IF
END SUBROUTINE OpenColumnProfile

SUBROUTINE OutputColumnProfile(VecC,Time)   
  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind) :: Time

  TYPE(ScalarCell_T) :: ScalarC(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  INTEGER :: ScalarPos
  INTEGER :: iProf
  INTEGER :: ncid
  TYPE(ColumnProfileTable_T), POINTER :: TableC(:)
  INTEGER, POINTER :: Scalar_VarId(:,:)
  INTEGER :: rec_varid
  REAL(RealKind) :: SaeuleProfile(Domain%nz,MaxLenTable)
  REAL(RealKind) :: SaeuleProfile0(Domain%nz,MaxLenTable)
  INTEGER :: ix,iy,iz,igz
  INTEGER :: ibLoc1,ib1,ibW
  INTEGER :: iProc,iProc1
  INTEGER :: Status(MPI_STATUS_SIZE)

  LOGICAL, SAVE :: Init=.TRUE.
 
  IF (Init) THEN
    CALL OpenColumnProfile(NamePointFile)
    Init=.FALSE.
    TimeStep=0
    IF (StartTime>OutputTimeColumnStart) OutputTimeColumnStart=StartTime
    OutputTimeColumn=OutputTimeColumnStart
  END IF  
  
  IF (OutputTimeColumn<=Time.AND.OutputTimeColumn<=OutputTimeColumnEnd) THEN 
    OutputTimeColumn=OutputTimeColumn+dtimeColumn 
    TimeStep=TimeStep+1
    DO iProf=1,SIZE(ColumnProfileWrite)
      TableC=>ColumnProfile(iProf)%TableC
      IF (MyId==0) THEN
        ncid=ColumnProfileWrite(iProf)%ncid
        rec_varid=ColumnProfileWrite(iProf)%rec_varid
        Scalar_varid=>ColumnProfileWrite(iProf)%Scalar_varid
        CALL check(nf90_open(ColumnProfileWrite(iProf)%FileName,nf90_write,ncid))
        DO i=1,SIZE(TableC)
          IF (TableC(i)%ColumnProfile) THEN 
            CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'ColumnProfile',Scalar_varid(i,1)))
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
        IF (TableC(i)%ColumnProfile) THEN 
          IF (RhoPos>0) THEN
            CALL Assign(Rho,VecC,RhoPos)
            DO ibLoc1=1,ColumnProfile(iProf)%nbLoc
              ibLoc=ColumnProfile(iProf)%ibLoc(ibLoc1)
              ib=LocGlob(ibLoc)
              CALL Set(Floor(ib))
              ix=ColumnProfile(iProf)%ixProf(ibLoc1)
              iy=ColumnProfile(iProf)%iyProf(ibLoc1)
              DO iz=iz0+1,iz1
                DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ) 
                  ColumnProfile(iProf)%ProfileBlock(ibLoc1)%Profile(igz,i)=ScalarC(ibLoc)%c(ix,iy,iz,1)& 
                  & /(Rho(ibLoc)%c(ix,iy,iz,1)+Eps)
                END DO  
              END DO  
            END DO  
          ELSE   
            DO ibLoc1=1,ColumnProfile(iProf)%nbLoc
              ibLoc=ColumnProfile(iProf)%ibLoc(ibLoc1)
              ib=LocGlob(ibLoc)
              CALL Set(Floor(ib))
              ix=ColumnProfile(iProf)%ixProf(ibLoc1)
              iy=ColumnProfile(iProf)%iyProf(ibLoc1)
              DO iz=iz0+1,iz1
                DO igz=(iz-1)*2**(-RefineZ)+1,iz*2**(-RefineZ) 
                  ColumnProfile(iProf)%ProfileBlock(ibLoc1)%Profile(igz,i)=ScalarC(ibLoc)%c(ix,iy,iz,1)
                END DO  
              END DO  
            END DO  
          END IF
        END IF
      END DO  
      IF (MyId==0) THEN
        DO ib1=1,ColumnProfileWrite(iProf)%nb
          ib=ColumnProfileWrite(iProf)%ib(ib1)
          iProc=ColumnProfileWrite(iProf)%iProc(ib1)
          CALL Set(Floor(ib))
          IF (iProc>0) THEN
            CALL MPI_RECV(ColumnProfileWrite(iProf)%Profile(igz0+1:igz1,:),nz*SIZE(ColumnProfile(iProf)%TableC) &
                         ,MPI_RealKind,iProc,iProc+ib*nb,MPI_Comm_World,Status,MPIErr)
          END IF  
        END DO
      END IF  
      IF (ColumnProfile(iProf)%nbLoc>0) THEN
        DO ibLoc1=1,ColumnProfile(iProf)%nbLoc
          ibLoc=ColumnProfile(iProf)%ibLoc(ibLoc1)
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          IF (MyId>0) THEN
            CALL MPI_SEND(ColumnProfile(iProf)%ProfileBlock(ibLoc1)%Profile,nz*SIZE(ColumnProfile(iProf)%TableC),&
            &MPI_RealKind,0,MyId+ib*nb,MPI_Comm_World,MPIErr)
          ELSE
            DO ib1=1,ColumnProfileWrite(iProf)%nb
              ibW=ColumnProfileWrite(iProf)%ib(ib1)
              IF (ibW==ib) THEN
                ColumnProfileWrite(iProf)%Profile(igz0+1:igz1,:)=ColumnProfile(iProf)%ProfileBlock(ibLoc1)%Profile
              END IF  
            END DO  
          END IF  
        END DO  
      END IF  
      IF (MyId==0) THEN
        DO i=1,SIZE(TableC)
          IF (TableC(i)%ColumnProfile) THEN 
            CALL check(nf90_put_var(ncid,Scalar_varid(i,1),ColumnProfileWrite(iProf)%Profile(:,i),start=(/1,TimeStep/)))
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
END SUBROUTINE OutputColumnProfile

SUBROUTINE CloseColumnProfile
  INTEGER :: iProf
  INTEGER :: ncid
  TYPE(ColumnProfileTable_T), POINTER :: TableC(:)
  
  IF (MyId==0) THEN
    DO iProf=1,SIZE(ColumnProfileWrite)
      ncid=ColumnProfileWrite(iProf)%ncid
      CALL check(nf90_close(ncid)) ! close: save new netCDF dataset
    END DO 
  END IF
END SUBROUTINE CloseColumnProfile

SUBROUTINE ColumnProfileNMLOutput(FileName,Check)
  CHARACTER(*) :: FileName
  LOGICAL :: Check
  CHARACTER(600) :: Line
 
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ColumnProfileOutput')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ColumnProfileOutput)
      Check=.TRUE.
      EXIT
    END IF
  END DO
1 CONTINUE 
  CLOSE(InputUnit)
END SUBROUTINE ColumnProfileNMLOutput

SUBROUTINE ColumnProfileDataTable(FileName)
  CHARACTER(*) :: FileName

  CHARACTER*200 :: Line
  INTEGER :: i,j,LenTable,NPoints 

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*) Line
    IF (INDEX(Line,'#ColumnProfileTable')>0) THEN
      READ(InputUnit,*) NPoints 
      ALLOCATE(ColumnProfile(NPoints))
      DO i=1,NPoints
        READ(InputUnit,*) Line
        READ(InputUnit,*) ColumnProfile(i)%P%x,ColumnProfile(i)%P%y
        READ(InputUnit,*) LenTable 
        MaxLenTable=MAX(MaxLenTable,LenTable)
        ALLOCATE(ColumnProfile(i)%TableC(LenTable))
        READ(InputUnit,*) Line 
        DO j=1,LenTable
          READ(InputUnit,*) ColumnProfile(i)%TableC(j)%ScalarC,ColumnProfile(i)%TableC(j)%Units,&
          & ColumnProfile(i)%TableC(j)%Description &
                           ,ColumnProfile(i)%TableC(j)%ColumnProfile
        END DO  
        READ(InputUnit,*) Line
      END DO  
      EXIT
    END IF
  END DO
  CLOSE(InputUnit) 
  IF (MyId==0) THEN
    ALLOCATE(ColumnProfileWrite(NPoints))
    DO i=1,NPoints
      ALLOCATE(ColumnProfileWrite(i)%TableC(LenTable))
      DO j=1,LenTable
        ColumnProfileWrite(i)%TableC(j)=ColumnProfile(i)%TableC(j)
      END DO  
    END DO  
  END IF  
END SUBROUTINE ColumnProfileDataTable
END MODULE SingleColumnProfile_Mod
