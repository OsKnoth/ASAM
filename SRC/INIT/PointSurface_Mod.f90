MODULE PointSurface_Mod

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

  INTEGER, PRIVATE :: TimeStep
  INTEGER :: ScalarPos
  INTEGER :: MaxLenTable=0
  INTEGER, PRIVATE :: i
  CHARACTER(30)     :: NamePoint = " "
  CHARACTER(30)     :: NamePointFile = " "
  CHARACTER (LEN = *), PARAMETER :: TH_NAME           ="THETA"
  CHARACTER (LEN = *), PARAMETER :: REC_NAME          = "time"
  CHARACTER (LEN = *), PARAMETER :: REC_LONGNAME      = "time"
  CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
  CHARACTER (LEN = *), PARAMETER :: Z_UNITS    = "m"
  CHARACTER (LEN = *), PARAMETER :: REC_UNITS  = "seconds since 2000-01-01 00:00:00"
  REAL(RealKind), PRIVATE :: OutputTimePoint
  REAL(RealKind) :: dtimePoint
  REAL(RealKind) :: OutputTimePointStart
  REAL(RealKind) :: OutputTimePointEnd

  TYPE PointSurfaceTable_T
    CHARACTER(32) :: ScalarC=' '
    LOGICAL       :: PointSurface=.FALSE.
    CHARACTER(20) :: Units=''
    CHARACTER(30) :: Description=''
    INTEGER       :: iPosC=0
  END TYPE PointSurfaceTable_T
  TYPE ColumnPoint_T
    REAL(RealKind) :: x
    REAL(RealKind) :: y
  END TYPE ColumnPoint_T
  TYPE PointSurface_T
    INTEGER :: ibLoc
    TYPE(PointSurfaceTable_T), POINTER :: TableC(:)
    TYPE(ColumnPoint_T) :: P
    INTEGER :: ixProf
    INTEGER :: iyProf
    LOGICAL :: OutBlock=.FALSE.
    LOGICAL :: OutHeader=.FALSE.
    CHARACTER*50 :: FileName
    INTEGER :: ncid
    INTEGER :: rec_varid, px_varid, py_dimid
    INTEGER :: rec_dimid, px_dimid, py_varid
    INTEGER, POINTER :: Scalar_VarId(:,:)
  END TYPE  
  TYPE(PointSurface_T), ALLOCATABLE :: PointSurface(:)

  NAMELIST /PointSurfaceOutput/ OutputTimePointStart, &
                                 OutputTimePointEnd, &
                                 dtimePoint, &
                                 NamePointFile
CONTAINS

SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS

  IF(STATUS /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(STATUS))
  END IF
END SUBROUTINE check

SUBROUTINE OpenPointSurface(FileName)
  CHARACTER(*) :: FileName

  INTEGER :: ix,iy,ixProf,iyProf
  INTEGER :: iProf,ibProfLoc,ib,ibProf
  TYPE(PointSurfaceTable_T), POINTER :: TableC(:)
  INTEGER, POINTER :: Scalar_VarId(:,:)
  INTEGER :: iPx,iPy
  CHARACTER*12 :: NamePx,NamePy
  INTEGER :: ncid
  INTEGER :: rec_varid, px_varid, py_varid
  INTEGER :: rec_dimid, px_dimid, py_dimid

  S1:DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        DO iProf=1,SIZE(PointSurface)
          IF (.NOT.PointSurface(iProf)%OutBlock) THEN
            xProf=PointSurface(iProf)%P%x
            yProf=PointSurface(iProf)%P%y
            S2:DO ix=ix0+1,ix1
               S3:DO iy=iy0+1,iy1
                IF (xP(ix-1)<=xProf.AND. &
                  xP(ix  )>=xProf.AND. &
                  yP(iy-1)<=yProf.AND. &
                  yP(iy  )>=yProf) THEN
                  PointSurface(iProf)%ixProf=ix
                  PointSurface(iProf)%iyProf=iy
                  ibProf=ib
                  ibProfLoc=ibLoc
                  PointSurface(iProf)%OutBlock=.TRUE.
                  PointSurface(iProf)%ibLoc=ibLoc
                  IF (igz0==Domain%igz0) THEN
                    PointSurface(iProf)%OutHeader=.TRUE.
                  END IF   
                  EXIT S2  
                END IF
              END DO S3
            END DO S2
          END IF
        END DO
  END DO S1

  DO iProf=1,SIZE(PointSurface)
    IF (PointSurface(iProf)%OutHeader) THEN
      iPx=PointSurface(iProf)%P%x
      iPy=PointSurface(iProf)%P%y
      WRITE(NamePx,'(I10)') iPx
      WRITE(NamePy,'(I10)') iPy
      NamePoint='Px'//TRIM(ADJUSTL(NamePx))//'Py'//TRIM(ADJUSTL(NamePy))
      PointSurface(iProf)%FileName=TRIM(FileName)//TRIM(ADJUSTL(NamePoint))//'.nc'
    ! WRITE(*,*) 'Name ',PointSurface(iProf)%FileName
      CALL check(nf90_create(TRIM(FileName)//TRIM(ADJUSTL(NamePoint))//'.nc',nf90_clobber,ncid))           ! create netCDF dataset: enter define mode
 
      TableC=>PointSurface(iProf)%TableC
   
      CALL check(nf90_def_dim(ncid,TRIM('time'),NF90_UNLIMITED,rec_dimid))     ! time
      CALL check(nf90_def_dim(ncid,TRIM('Px'),1,px_dimid))
      CALL check(nf90_def_dim(ncid,TRIM('Py'),1,py_dimid))

      CALL check(nf90_def_var(ncid,TRIM('time'),NF90_DOUBLE,rec_dimid,rec_varid))
      CALL check(nf90_def_var(ncid,TRIM('Px'),NF90_INT,px_dimid,px_varid))
      CALL check(nf90_def_var(ncid,TRIM('Py'),NF90_INT,py_dimid,py_varid))
  
      CALL check(nf90_put_att(ncid,rec_varid,'long_name',REC_NAME))
      CALL check(nf90_put_att(ncid,rec_varid, UNITS, REC_UNITS))
      CALL check(nf90_put_att(ncid,rec_varid,'calender','noleap'))
      CALL check(nf90_put_att(ncid,px_varid,UNITS,'m'))
      CALL check(nf90_put_att(ncid,py_varid,UNITS,'m'))
      
      ALLOCATE(PointSurface(iProf)%Scalar_varid(SIZE(TableC),1))
      Scalar_VarId=>PointSurface(iProf)%Scalar_VarId
 
      DO i=1,SIZE(TableC)
        IF (TableC(i)%PointSurface) THEN 
          CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'PointSurface',NF90_DOUBLE &
            ,rec_dimid,Scalar_varid(i,1)))
          CALL check(nf90_put_att(ncid,Scalar_varid(i,1),UNITS,TRIM(TableC(i)%Units)))
          CALL check(nf90_put_att(ncid,Scalar_varid(i,1),'Description',TRIM(TableC(i)%Description)))
        END IF  
      END DO  

      CALL check(nf90_enddef(ncid))                                        ! end definitions: leave define mode

      CALL check(nf90_put_var(ncid,px_varid,iPx))
      CALL check(nf90_put_var(ncid,py_varid,iPy))
      CALL check(nf90_close(ncid))                                        ! close: save new netCDF dataset

      IF (StartTime>OutputTimePointStart) OutputTimePointStart=StartTime
      OutputTimePoint=OutputTimePointStart

      PointSurface(iProf)%ncid=ncid
      PointSurface(iProf)%rec_dimid=rec_dimid
      PointSurface(iProf)%rec_varid=rec_varid
    END IF
  END DO     
END SUBROUTINE OpenPointSurface

SUBROUTINE OutputPointSurface(VecC,Time)   
  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind) :: Time

  TYPE(ScalarCell_T) :: ScalarC(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  INTEGER :: ScalarPos
  INTEGER :: iProf
  INTEGER :: ncid
  TYPE(PointSurfaceTable_T), POINTER :: TableC(:)
  INTEGER, POINTER :: Scalar_VarId(:,:)
  INTEGER :: rec_varid
  REAL(RealKind) :: PointValues(MaxLenTable)
  INTEGER :: ix,iy
  INTEGER :: ixLoc,iyLoc
  LOGICAL :: SpecialOut

  LOGICAL, SAVE :: Init=.TRUE.
  INTEGER, SAVE :: ibound=0
 
! WRITE(*,*) MyId,'Anfang WritePointSurface'
  IF (Init) THEN
    CALL OpenPointSurface(NamePointFile)
    Init=.FALSE.
    TimeStep=0
!   Get BoundCellValue iBound
    ix=PointSurface(1)%ixProf
    iy=PointSurface(1)%iyProf
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
 !    WRITE(*,*) 'ix,iy',ix,iy,NumBoundCell
      DO i=1,NumBoundCell
        ixLoc=Floor(ib)%BoundCell(i)%ix
        iyLoc=Floor(ib)%BoundCell(i)%iy
 !      WRITE(*,*) 'boundcells ',i,NumBoundCell,ixLoc,iyLoc
        IF (ixLoc==ix.AND.iyLoc==iy) THEN
          iBound=i
          EXIT
        END IF
      END DO 
    END DO
  END IF  
  
  SpecialOut=.FALSE.

! WRITE(*,*) MyId,'Mitte WritePointSurface',OutputTimePoint,OutputTimePoint+dtimePoint
  IF (OutputTimePoint<=Time.AND.OutputTimePoint<=OutputTimePointEnd) THEN 
!   WRITE(*,*) MyId,'Mitte 1 WritePointSurface' 
    OutputTimePoint=OutputTimePoint+dtimePoint 
    TimeStep=TimeStep+1
!   WRITE(*,*) 'OutputTimePoint,TimeStep',OutputTimePoint,TimeStep
    DO iProf=1,SIZE(PointSurface)
      IF (PointSurface(iProf)%OutHeader) THEN
        ncid=PointSurface(iProf)%ncid
        rec_varid=PointSurface(iProf)%rec_varid
        TableC=>PointSurface(iProf)%TableC
        Scalar_varid=>PointSurface(iProf)%Scalar_varid
        CALL check(nf90_open(PointSurface(iProf)%FileName,nf90_write,ncid))
        DO i=1,SIZE(TableC)
          IF (TableC(i)%PointSurface) THEN 
            CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'PointSurface',Scalar_varid(i,1)))
          END IF  
        END DO
   
        DO i=1,SIZE(TableC)
          ScalarPos=Position(TRIM(TableC(i)%ScalarC))
!         WRITE(*,*) 'ScalarPos,ScalarC',ScalarPos,TRIM(TableC(i)%ScalarC)
!         WRITE(*,*) 'i,Size,iBound',i,SIZE(TableC),iBound
          IF (ScalarPos>0) THEN
            CALL Assign(ScalarC,VecC,ScalarPos)
          ELSE IF (TRIM(TableC(i)%ScalarC)=='U') THEN  
            CALL Assign(ScalarC,VecC,uPosl)
          ELSE IF (TRIM(TableC(i)%ScalarC)=='V') THEN
            CALL Assign(ScalarC,VecC,vPosl)
          ELSE IF (TRIM(TableC(i)%ScalarC)=='W') THEN
            CALL Assign(ScalarC,VecC,wPosl)
          ELSE IF (TRIM(TableC(i)%ScalarC)=='QDIR') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%RadDirekt
          ELSE IF (TRIM(TableC(i)%ScalarC)=='QDIF') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%RadDiffus
          ELSE IF (TRIM(TableC(i)%ScalarC)=='QINF') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%RadInfred
          ELSE IF (TRIM(TableC(i)%ScalarC)=='SHF') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%FluxSens
          ELSE IF (TRIM(TableC(i)%ScalarC)=='LHF') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%FluxLat
          ELSE IF (TRIM(TableC(i)%ScalarC)=='CD') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%DragM
          ELSE IF (TRIM(TableC(i)%ScalarC)=='CH') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%DragH
          ELSE IF (TRIM(TableC(i)%ScalarC)=='CE') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%DragQ
          ELSE IF (TRIM(TableC(i)%ScalarC)=='TS') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%TSoil
          ELSE IF (TRIM(TableC(i)%ScalarC)=='QVS') THEN
            SpecialOut=.TRUE.
            PointValues(i)=BoundCell(iBound)%QvSoil
          END IF  
          IF (TableC(i)%PointSurface.AND..NOT.SpecialOut) THEN 
        !   WRITE(*,*) 'not specialout ',TRIM(TableC(i)%ScalarC)
            IF (RhoPos>0) THEN
              CALL Assign(Rho,VecC,RhoPos)
              ibLoc=PointSurface(iProf)%ibLoc
              ib=LocGlob(ibLoc)
              CALL Set(Floor(ib))
              ix=PointSurface(iProf)%ixProf
              iy=PointSurface(iProf)%iyProf
              PointValues(i)=ScalarC(ibLoc)%c(ix,iy,iz0+1,1)/(Rho(ibLoc)%c(ix,iy,iz0+1,1)+Eps)
            ELSE   
              ib=LocGlob(PointSurface(iProf)%ibLoc)
              CALL Set(Floor(ib))
              ix=PointSurface(iProf)%ixProf
              iy=PointSurface(iProf)%iyProf
              PointValues(i)=ScalarC(ibLoc)%c(ix,iy,iz0+1,1)
           !  WRITE (*,*) 'PointProfile',ix,iy
            END IF
          END IF
        END DO  
        
        DO i=1,SIZE(TableC)
          IF (TableC(i)%PointSurface) THEN 
            CALL check(nf90_put_var(ncid,Scalar_varid(i,1),PointValues(i),start=(/TimeStep/)))
          END IF  
        END DO
        CALL check(nf90_inq_varid(ncid,TRIM('time'),rec_varid))
        CALL check(nf90_put_var(ncid,rec_varid,Time,start=(/TimeStep/)))
        CALL check(nf90_sync(ncid))
        CALL check(nf90_close(ncid))
      END IF
    END DO
  END IF
! WRITE(*,*) MyId,'Ende WritePointSurface' 
  ! SensFluxCell(ibLoc)%cB(i,1)=BoundCell(i)%FluxSens
  ! LatFluxCell(ibLoc)%cB(i,1)=BoundCell(i)%FluxLat
  ! LandClassCell(ibLoc)%cB(i,1)=BoundCell(i)%LandClass
  ! RoughnessLengthCell(ibLoc)%cB(i,1)=BoundCell(i)%zRauh
  ! AlbedoCell(ibLoc)%cB(i,1)=BoundCell(i)%alb
  ! BulkCoeffDragCell(ibLoc)%cB(i,1)=BoundCell(i)%DragM
  ! BulkCoeffHeatCell(ibLoc)%cB(i,1)=BoundCell(i)%DragH
  ! BulkCoeffMoistCell(ibLoc)%cB(i,1)=BoundCell(i)%DragQ
END SUBROUTINE OutputPointSurface

SUBROUTINE ClosePointSurface
  INTEGER :: iProf
  INTEGER :: ncid
  TYPE(PointSurfaceTable_T), POINTER :: TableC(:)
  
  DO iProf=1,SIZE(PointSurface)
    IF (PointSurface(iProf)%OutHeader) THEN
      ncid=PointSurface(iProf)%ncid
      CALL check(nf90_close(ncid)) ! close: save new netCDF dataset
    END IF
  END DO 
END SUBROUTINE ClosePointSurface

SUBROUTINE PointSurfaceNMLOutput(FileName,Check)
  CHARACTER(*) :: FileName
  LOGICAL :: Check
  CHARACTER(600) :: Line
 
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&PointSurfaceOutput')>0) THEN
      WRITE(*,*) Line
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=PointSurfaceOutput)
      Check=.TRUE.
      EXIT
    END IF
  END DO
1 CONTINUE 
  CLOSE(InputUnit)
END SUBROUTINE PointSurfaceNMLOutput

SUBROUTINE PointSurfaceDataTable(FileName)
  CHARACTER(*) :: FileName

  CHARACTER(1000) :: Line
  INTEGER :: i,j,LenTable,NPoints 

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*) Line
    IF (INDEX(Line,'#PointSurfaceTable')>0) THEN
      READ(InputUnit,*) NPoints 
      ALLOCATE(PointSurface(NPoints))
      DO i=1,NPoints
        READ(InputUnit,*) Line
        READ(InputUnit,*) PointSurface(i)%P%x,PointSurface(i)%P%y
        READ(InputUnit,*) LenTable 
        MaxLenTable=MAX(MaxLenTable,LenTable)
        ALLOCATE(PointSurface(i)%TableC(LenTable))
        READ(InputUnit,*) Line 
        DO j=1,LenTable
          READ(InputUnit,*) PointSurface(i)%TableC(j)%ScalarC,PointSurface(i)%TableC(j)%Units,&
          &PointSurface(i)%TableC(j)%Description &
                           ,PointSurface(i)%TableC(j)%PointSurface
        END DO  
        READ(InputUnit,*) Line
      END DO  
      EXIT
    END IF
  END DO
  CLOSE(InputUnit) 
END SUBROUTINE PointSurfaceDataTable
END MODULE PointSurface_Mod
