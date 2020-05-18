MODULE MeanSurface_Mod

  USE Kind_Mod
  USE typesizes
  USE netcdf
  USE Domain_Mod
  USE Floor_Mod
  USE DataType_Mod
  USE Physics_Mod
  USE ScalarVectorCellPar_Mod
  USE Output_Mod
  
  IMPLICIT NONE

  INTEGER, PRIVATE :: ncid
  INTEGER, PRIVATE :: rec_dimid
  INTEGER, PRIVATE :: TimeStep
  INTEGER :: ScalarPos
  INTEGER, PRIVATE :: i
  INTEGER, PRIVATE :: rec_varid, th_varid, w_varid 
  INTEGER, PRIVATE, ALLOCATABLE :: Scalar_varid(:,:)
  CHARACTER*30     :: FileNameMeanSurface = ""
  CHARACTER (LEN = *), PARAMETER :: TH_NAME       ="THETA"
  CHARACTER (LEN = *), PARAMETER :: REC_NAME          = "time"
  CHARACTER (LEN = *), PARAMETER :: REC_LONGNAME      = "time"
  CHARACTER (LEN = *), PARAMETER :: TH_LONGNAME      = "potential temperature"
  CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
  CHARACTER (LEN = *), PARAMETER :: REC_UNITS  = "seconds since 2000-01-01 00:00:00"
  CHARACTER (LEN = *), PARAMETER :: TH_UNITS   = "K"  
  CHARACTER (LEN = *), PARAMETER :: W_NAME     = "W"
  CHARACTER (LEN = *), PARAMETER :: W_UNITS    = "m/s"
  REAL(RealKind), PRIVATE :: OutputTimeMeanSurface
  REAL(RealKind) :: dtimeMeanSurface
  REAL(RealKind) :: OutputTimeMeanSurfaceStart
  REAL(RealKind) :: OutputTimeMeanSurfaceEnd

  TYPE MeanSurfaceTable_T
    CHARACTER(32) :: ScalarC=' '
    LOGICAL       :: MeanSurface=.FALSE.
    !LOGICAL       :: VertFlux=.FALSE.
    !LOGICAL       :: TotalVertFlux=.FALSE.
    !LOGICAL       :: Variance=.FALSE.
    CHARACTER(20) :: Units=''
    CHARACTER(30) :: Description=''
    INTEGER       :: iPosC=0
  END TYPE MeanSurfaceTable_T
  TYPE(MeanSurfaceTable_T), ALLOCATABLE :: TableC(:)

  NAMELIST /SurfaceOutput/ FileNameMeanSurface, &
                           OutputTimeMeanSurfaceStart, &
                           OutputTimeMeanSurfaceEnd, &
                           dtimeMeanSurface
CONTAINS

SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS

  IF(STATUS /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(STATUS))
  END IF
END SUBROUTINE check

SUBROUTINE OpenMeanSurface(FileName)
  CHARACTER(*) :: FileName

  CALL check(nf90_create(TRIM(FileName)//'.nc',nf90_clobber,ncid))           ! create netCDF dataset: enter define mode
 
  CALL check(nf90_def_dim(ncid,TRIM('time'),NF90_UNLIMITED,rec_dimid))     ! time

  CALL check(nf90_def_var(ncid,TRIM('time'),NF90_DOUBLE,rec_dimid,rec_varid)) 
  
  CALL check(nf90_put_att(ncid, rec_varid,"long_name",REC_NAME))
  CALL check(nf90_put_att(ncid, rec_varid, UNITS, REC_UNITS))
  CALL check(nf90_put_att(ncid, rec_varid,"calender","noleap"))
 
  ALLOCATE(Scalar_varid(SIZE(TableC),4))
  DO i=1,SIZE(TableC)
    IF (TableC(i)%MeanSurface) THEN 
      CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'MeanSurface',NF90_DOUBLE,rec_dimid,Scalar_varid(i,1)))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,1),UNITS,TRIM(TableC(i)%Units)))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,1),'Description',TRIM(TableC(i)%Description)))
    END IF  
    !IF (TableC(i)%VertFlux) THEN 
     ! CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'VertFlux',NF90_REAL,rec_dimid,Scalar_varid(i,2)))
      !CALL check(nf90_put_att(ncid,Scalar_varid(i,2),UNITS,TRIM(TableC(i)%Units)//' '//W_UNITS))
      !CALL check(nf90_put_att(ncid,Scalar_varid(i,2),'Description',TRIM('Vertical Flux of '//TableC(i)%Description)))
      !CALL check(nf90_put_att(ncid,Scalar_varid(i,2),'Description',TRIM('deviation from horizontal mean (double prime): resolved correlation')))
    !END IF  
    !IF (TableC(i)%TotalVertFlux) THEN
      !CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'TotalVertFlux',NF90_REAL,rec_dimid,Scalar_varid(i,3)))
      !CALL check(nf90_put_att(ncid,Scalar_varid(i,3),UNITS,TRIM(TableC(i)%Units)//' '//W_UNITS))
      !CALL check(nf90_put_att(ncid,Scalar_varid(i,3),'Description',TRIM('Total Vertical Flux of '//TableC(i)%Description)))
    !END IF
    !IF (TableC(i)%Variance) THEN
      !CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'Variance',NF90_REAL,rec_dimid,Scalar_varid(i,4)))
      !CALL check(nf90_put_att(ncid,Scalar_varid(i,4),UNITS,'('//TRIM(TableC(i)%Units)//')^2'))
      !CALL check(nf90_put_att(ncid,Scalar_varid(i,4),'Description',TRIM('Variance of '//TableC(i)%Description)))
    !END IF
  END DO  

  CALL check(nf90_enddef(ncid))                                        ! end definitions: leave define mode

  CALL check(nf90_close(ncid))                                        ! close: save new netCDF dataset


END SUBROUTINE OpenMeanSurface

SUBROUTINE OutputMeanSurface(VecC,Time)   
  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind) :: Time

  REAL(RealKind) :: MeanSurface(SIZE(TableC)) !
  !REAL(RealKind) :: MeanVarUProfile(SIZE(TableC))
  !REAL(RealKind) :: MeanVarVProfile(SIZE(TableC))
  !REAL(RealKind) :: MeanVarWProfile(SIZE(TableC))
  !REAL(RealKind) :: TotMeanVarUProfile(SIZE(TableC))
  !REAL(RealKind) :: TotMeanVarVProfile(SIZE(TableC))
  !REAL(RealKind) :: TotMeanVarWProfile(SIZE(TableC))
  !REAL(RealKind) :: MeanVarianceProfile(SIZE(TableC))

  TYPE(ScalarCell_T) :: ScalarC(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  INTEGER :: ScalarPos
  LOGICAL :: SpecialOut

  LOGICAL, SAVE :: Init=.TRUE.
 
  IF (Init) THEN
    IF (MyId==0) THEN
      CALL OpenMeanSurface(FileNameMeanSurface)
    END IF
    Init=.FALSE.
    TimeStep=0
    IF (StartTime>OutputTimeMeanSurfaceStart) OutputTimeMeanSurfaceStart=StartTime
    OutputTimeMeanSurface=OutputTimeMeanSurfaceStart
  END IF
 
  SpecialOut=.FALSE.

  IF (OutputTimeMeanSurface<=Time.AND.OutputTimeMeanSurface<=OutputTimeMeanSurfaceEnd) THEN 
    OutputTimeMeanSurface=OutputTimeMeanSurface+dtimeMeanSurface
    TimeStep=TimeStep+1
    IF (MyId==0) THEN 
      CALL check(nf90_open(TRIM(FileNameMeanSurface)//'.nc',nf90_write,ncid))
 
      !DO i=1,SIZE(TableC)
       ! IF (TableC(i)%MeanSurface) THEN 
       !   CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'MeanSurface',Scalar_varid(i,1)))
       ! END IF  
       ! IF (TableC(i)%VertFlux) THEN 
       !   CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'VertFlux',Scalar_varid(i,2)))
       ! END IF 
       ! IF (TableC(i)%TotalVertFlux) THEN
       !   CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'TotalVertFlux',Scalar_varid(i,3)))
       ! END IF
       ! IF (TableC(i)%Variance) THEN
       !   CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'Variance',Scalar_varid(i,4)))
       ! END IF
    END IF 
   
    IF (RhoPos>0) THEN
      CALL Assign(Rho,VecC,RhoPos)
    END IF  
    DO i=1,SIZE(TableC)
      ScalarPos=Position(TRIM(TableC(i)%ScalarC))
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
        CALL AreaSurfaceBC(BoundCell,1,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='QDIF') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,2,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='QINF') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,3,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='SHF') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,4,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='LHF') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,5,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='CD') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,6,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='CH') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,7,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='CE') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,8,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='TS') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,9,MeanSurface(i))
      ELSE IF (TRIM(TableC(i)%ScalarC)=='QVS') THEN
        SpecialOut=.TRUE.
        CALL AreaSurfaceBC(BoundCell,10,MeanSurface(i))
      END IF  
      IF (TableC(i)%MeanSurface.AND..NOT.SpecialOut) THEN 
        IF (RhoPos>0) THEN
          CALL AreaSurface(ScalarC,MeanSurface(i),Rho) !
        ELSE  
          CALL AreaSurface(ScalarC,MeanSurface(i)) !
        END IF
      END IF
      !IF (TableC(i)%VertFlux) THEN 
      !  IF (RhoPos>0) THEN
      !    CALL AreaVarianceProfile(uCell,Rho,ScalarC,Rho,MeanVarUProfile(i)) !
      !    CALL AreaVarianceProfile(vCell,Rho,ScalarC,Rho,MeanVarVProfile(i)) !
      !    CALL AreaVarianceProfile(wCell,Rho,ScalarC,Rho,MeanVarWProfile(i)) !
      !  END IF    
      !END IF
      !IF (TableC(i)%TotalVertFlux) THEN
      !  IF (RhoPos>0) THEN
      !    CALL AreaVarianceTotalProfile(uCell,Rho,ScalarC,Rho,TotMeanVarUProfile(i)) !
      !    CALL AreaVarianceTotalProfile(vCell,Rho,ScalarC,Rho,TotMeanVarVProfile(i)) !
      !    CALL AreaVarianceTotalProfile(wCell,Rho,ScalarC,Rho,TotMeanVarWProfile(i)) !
      !  END IF
      !END IF
      !IF (TableC(i)%Variance) THEN
      !  IF (RhoPos>0) THEN
      !    CALL AreaVarianceProfile(ScalarC,Rho,ScalarC,Rho,MeanVarianceProfile(i)) !
      !  END IF
      !END IF
    END DO  

    IF (MyID==0) THEN 
      WRITE(*,*) 'WriteMeanSurface' 
      DO i=1,SIZE(TableC)
        IF (TableC(i)%MeanSurface) THEN 
          CALL check(nf90_put_var(ncid,Scalar_varid(i,1),MeanSurface(i),start=(/TimeStep/))) !
        END IF  
      END DO
      CALL check(nf90_inq_varid(ncid,TRIM('time'),rec_varid))
      CALL check(nf90_put_var(ncid,rec_varid,Time,start=(/TimeStep/)))
      CALL check(nf90_sync(ncid))
      CALL check(nf90_close(ncid))
    END IF 
  END IF  
END SUBROUTINE OutputMeanSurface

SUBROUTINE CloseMeanSurface
  CALL check(nf90_close(ncid) ) ! close: save new netCDF dataset
END SUBROUTINE CloseMeanSurface

SUBROUTINE MeanSurfaceNMLOutput(FileName,Check)
  CHARACTER(*) :: FileName
  LOGICAL :: Check
  CHARACTER(600) :: Line
 
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&SurfaceOutput')>0) THEN
      WRITE(*,*) Line
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=SurfaceOutput)
      Check=.TRUE.
      EXIT
    END IF
  END DO
1 CONTINUE 
  CLOSE(InputUnit)
END SUBROUTINE MeanSurfaceNMLOutput

SUBROUTINE MeanSurfaceDataTable(FileName)
  CHARACTER(*) :: FileName

  CHARACTER(1000) :: Line
  INTEGER :: i,LenTable 

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*) Line
    IF (INDEX(Line,'#MeanSurfaceTable')>0) THEN
      READ(InputUnit,*) Lentable 
      ALLOCATE(TableC(LenTable))
      READ(InputUnit,*) Line 
      DO i=1,LenTable
        READ(InputUnit,*) TableC(i)%ScalarC,TableC(i)%Units,TableC(i)%Description,TableC(i)%MeanSurface
      END DO  
      EXIT
    END IF
  END DO
  CLOSE(InputUnit) 
END SUBROUTINE MeanSurfaceDataTable
END MODULE MeanSurface_Mod
