MODULE OutputMeanProfile_Mod

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

  INTEGER, PARAMETER :: NDIMS  = 2
  INTEGER, PRIVATE :: ncid, dimids(NDIMS)
  INTEGER, PRIVATE :: z_dimid, rec_dimid
  INTEGER, PRIVATE :: TimeStep
  INTEGER :: ScalarPos
  INTEGER, PRIVATE :: i
  INTEGER, PRIVATE :: z_varid, rec_varid, th_varid, w_varid 
  INTEGER, PRIVATE, ALLOCATABLE :: Scalar_varid(:,:)
  CHARACTER*128     :: FileNameMeanProfile = ""
  CHARACTER (LEN = *), PARAMETER :: TH_NAME       ="THETA"
  CHARACTER (LEN = *), PARAMETER :: Z_NAME            = "lev"
  CHARACTER (LEN = *), PARAMETER :: Z_LONGNAME        = "Height at midpoints"
  CHARACTER (LEN = *), PARAMETER :: REC_NAME          = "time"
  CHARACTER (LEN = *), PARAMETER :: REC_LONGNAME      = "time"
  CHARACTER (LEN = *), PARAMETER :: TH_LONGNAME      = "potential temperature"
  CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
  CHARACTER (LEN = *), PARAMETER :: Z_UNITS    = "m"
  CHARACTER (LEN = *), PARAMETER :: REC_UNITS  = "seconds since 2000-01-01 00:00:00"
  CHARACTER (LEN = *), PARAMETER :: TH_UNITS   = "K"  
  CHARACTER (LEN = *), PARAMETER :: W_NAME     = "W"
  CHARACTER (LEN = *), PARAMETER :: W_UNITS    = "m/s"
  REAL(RealKind), PRIVATE :: OutputTimeMean
  REAL(RealKind) :: dtimeMean
  REAL(RealKind) :: OutputTimeMeanStart
  REAL(RealKind) :: OutputTimeMeanEnd
  LOGICAL        :: OutputSubdomain=.FALSE.
  REAL(RealKind) :: OutputXStart
  REAL(RealKind) :: OutputXEnd
  REAL(RealKind) :: OutputYStart
  REAL(RealKind) :: OutputYEnd

  TYPE MeanProfileTable_T
    CHARACTER(32) :: ScalarC=' '
    LOGICAL       :: MeanProfile=.FALSE.
    LOGICAL       :: VertFlux=.FALSE.
    LOGICAL       :: TotalVertFlux=.FALSE.
    LOGICAL       :: Variance=.FALSE.
    CHARACTER(20) :: Units=''
    CHARACTER(30) :: Description=''
    INTEGER       :: iPosC=0
  END TYPE MeanProfileTable_T
  TYPE(MeanProfileTable_T), ALLOCATABLE :: TableC(:)

  NAMELIST /ProfileOutput/ FileNameMeanProfile, &
                           OutputTimeMeanStart, &
                           OutputTimeMeanEnd, &
                           OutputSubdomain, &
                           OutputXStart, &
                           OutputXEnd, &
                           OutputYStart, &
                           OutputYEnd, &
                           dtimeMean
CONTAINS

SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS

  IF(STATUS /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(STATUS))
      !      STOP "Stopped"
  END IF
END SUBROUTINE check

SUBROUTINE OpenMeanProfile(FileName)
  CHARACTER(*) :: FileName

  REAL(RealKind) :: zM(Domain%nz)
  INTEGER :: iz

  CALL check(nf90_create(TRIM(FileName)//'.nc',nf90_clobber,ncid))           ! create netCDF dataset: enter define mode
! WRITE (*,*) 'OpenMeanProfile: ncid',ncid,nf90_clobber,TRIM(FileName),Domain%nz
 
  CALL check(nf90_def_dim(ncid,TRIM('zM'),Domain%nz,z_dimid))                ! define dimensions: from name and length=30
  CALL check(nf90_def_dim(ncid,TRIM('time'), NF90_UNLIMITED, rec_dimid))     ! time

  CALL check(nf90_def_var(ncid,TRIM('zM'),NF90_DOUBLE,z_dimid,z_varid))      ! define variables: from name, type, dims
  CALL check(nf90_def_var(ncid,TRIM('time'),NF90_DOUBLE,rec_dimid,rec_varid)) 
  
  CALL check(nf90_put_att(ncid, z_varid,UNITS,Z_UNITS))                 ! assign attribute values, a numeric rank 1 array
  CALL check(nf90_put_att(ncid, z_varid,"long_name",Z_LONGNAME))
  CALL check(nf90_put_att(ncid, z_varid,"positive","up"))
  CALL check(nf90_put_att(ncid, rec_varid,"long_name",REC_NAME))
  CALL check(nf90_put_att(ncid, rec_varid, UNITS, REC_UNITS))
  CALL check(nf90_put_att(ncid, rec_varid,"calender","none"))

 
  dimids=(/z_dimid, rec_dimid/)
 
  ALLOCATE(Scalar_varid(SIZE(TableC),4))
  DO i=1,SIZE(TableC)
    IF (TableC(i)%MeanProfile) THEN 
      CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'MeanProfile',NF90_DOUBLE,dimids,Scalar_varid(i,1)))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,1),UNITS,TRIM(TableC(i)%Units)))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,1),'Description',TRIM(TableC(i)%Description)))
    END IF  
    IF (TableC(i)%VertFlux) THEN 
      CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'VertFlux',NF90_DOUBLE,dimids,Scalar_varid(i,2)))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,2),UNITS,TRIM(TableC(i)%Units)//' '//W_UNITS))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,2),'Description',TRIM('Vertical Flux of '//TableC(i)%Description)))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,2),'Description',&
      & TRIM('deviation from horizontal mean (double prime): resolved correlation')))
    END IF  
    IF (TableC(i)%TotalVertFlux) THEN
      CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'TotalVertFlux',NF90_DOUBLE,dimids,Scalar_varid(i,3)))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,3),UNITS,TRIM(TableC(i)%Units)//' '//W_UNITS))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,3),'Description',TRIM('Total Vertical Flux of '//TableC(i)%Description)))
    END IF
    IF (TableC(i)%Variance) THEN
      CALL check(nf90_def_var(ncid,TRIM(TableC(i)%ScalarC)//'Variance',NF90_DOUBLE,dimids,Scalar_varid(i,4)))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,4),UNITS,'('//TRIM(TableC(i)%Units)//')^2'))
      CALL check(nf90_put_att(ncid,Scalar_varid(i,4),'Description',TRIM('Variance of '//TableC(i)%Description)))
    END IF
  END DO  

  CALL check(nf90_enddef(ncid))                                        ! end definitions: leave define mode

  DO iz=Domain%iz0+1,Domain%iz1
    zM(iz)=0.5d0*(Domain%zP(iz-1)+Domain%zP(iz))
  END DO  
  CALL check(nf90_put_var(ncid,z_varid,zM))
  CALL check(nf90_close(ncid))                                        ! close: save new netCDF dataset

END SUBROUTINE OpenMeanProfile

SUBROUTINE OutputMeanProfile(VecC,Time)   
  TYPE(Vector4Cell_T) :: VecC(:)
  REAL(RealKind) :: Time

  REAL(RealKind) :: MeanProfile(Domain%nz,SIZE(TableC))
  REAL(RealKind) :: MeanVarUProfile(Domain%nz,SIZE(TableC))
  REAL(RealKind) :: MeanVarVProfile(Domain%nz,SIZE(TableC))
  REAL(RealKind) :: MeanVarWProfile(Domain%nz,SIZE(TableC))
  REAL(RealKind) :: TotMeanVarUProfile(Domain%nz,SIZE(TableC))
  REAL(RealKind) :: TotMeanVarVProfile(Domain%nz,SIZE(TableC))
  REAL(RealKind) :: TotMeanVarWProfile(Domain%nz,SIZE(TableC))
  REAL(RealKind) :: MeanVarianceProfile(Domain%nz,SIZE(TableC))

  TYPE(ScalarCell_T) :: ScalarC(nbLoc)
  TYPE(ScalarCell_T) :: Rho(nbLoc)
  INTEGER :: ScalarPos,j

  LOGICAL, SAVE :: Init=.TRUE.
 
  IF (Init) THEN
    IF (MyID==0) THEN
      CALL OpenMeanProfile(FileNameMeanProfile)
    END IF
    Init=.FALSE.
    TimeStep=0
    IF (StartTime>OutputTimeMeanStart) OutputTimeMeanStart=StartTime
    OutputTimeMean=OutputTimeMeanStart
  END IF  
 
  IF (OutputTimeMean<=Time.AND.OutputTimeMean<=OutputTimeMeanEnd) THEN 
  ! WRITE (*,*) MyID,Init,'OutputTime',OutputTimeMean,OutputTimeMean+dtimeMean
    OutputTimeMean=OutputTimeMean+dtimeMean 
    TimeStep=TimeStep+1
    IF (MyId==0) THEN 
      CALL check(nf90_open(TRIM(FileNameMeanProfile)//'.nc',nf90_write,ncid))
 
      DO i=1,SIZE(TableC)
        IF (TableC(i)%MeanProfile) THEN 
          CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'MeanProfile',Scalar_varid(i,1)))
        END IF  
        IF (TableC(i)%VertFlux) THEN 
          CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'VertFlux',Scalar_varid(i,2)))
        END IF 
        IF (TableC(i)%TotalVertFlux) THEN
          CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'TotalVertFlux',Scalar_varid(i,3)))
        END IF
        IF (TableC(i)%Variance) THEN
          CALL check(nf90_inq_varid(ncid,TRIM(TableC(i)%ScalarC)//'Variance',Scalar_varid(i,4)))
        END IF
      END DO
    END IF 
   
    IF (RhoPos>0) THEN
      CALL Assign(Rho,VecC,RhoPos)
    END IF  
    DO i=1,SIZE(TableC)
      ScalarPos=Position(TRIM(TableC(i)%ScalarC))
      IF (ScalarPos>0) THEN
        CALL Assign(ScalarC,VecC,ScalarPos)
      ELSE IF (TRIM(TableC(i)%ScalarC)=='U') THEN  
        CALL Assign(ScalarC,uCell)
      ELSE IF (TRIM(TableC(i)%ScalarC)=='V') THEN
        CALL Assign(ScalarC,vCell)
      ELSE IF (TRIM(TableC(i)%ScalarC)=='W') THEN
        CALL Assign(ScalarC,wCell)
      ELSE IF (TRIM(TableC(i)%ScalarC)=='T') THEN
        CALL Assign(ScalarC,TCell)
      END IF  
      IF (TableC(i)%MeanProfile) THEN 
        IF (RhoPos>0) THEN
          IF (TRIM(TableC(i)%ScalarC).NE.'RHO' &
              .AND.TRIM(TableC(i)%ScalarC).NE.'T'  &
              .AND.TRIM(TableC(i)%ScalarC).NE.'NV' &
              .AND.TRIM(TableC(i)%ScalarC).NE.'NC' &
              .AND.TRIM(TableC(i)%ScalarC).NE.'NR' &
              .AND.TRIM(TableC(i)%ScalarC).NE.'NI' &
              .AND.TRIM(TableC(i)%ScalarC).NE.'NS') THEN
            CALL AreaProfile(ScalarC,MeanProfile(:,i),Rho=Rho,Subdomain=OutputSubdomain,&
            & oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
          ELSE  
            CALL AreaProfile(ScalarC,MeanProfile(:,i),Subdomain=OutputSubdomain,& 
            & oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
          END IF
        END IF
      END IF
      IF (TableC(i)%VertFlux) THEN 
        IF (RhoPos>0) THEN
      !   WRITE(*,*) 'TableVertflux x3: ',TableC(i)%ScalarC
          CALL AreaVarianceProfile(uCell,Rho,ScalarC,Rho,MeanVarUProfile(:,i), &
          & Subdomain=OutputSubdomain,oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
          CALL AreaVarianceProfile(vCell,Rho,ScalarC,Rho,MeanVarVProfile(:,i), &
          & Subdomain=OutputSubdomain,oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
          CALL AreaVarianceProfile(wCell,Rho,ScalarC,Rho,MeanVarWProfile(:,i), & 
          & Subdomain=OutputSubdomain,oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
        END IF    
      END IF
      IF (TableC(i)%TotalVertFlux) THEN
        IF (RhoPos>0) THEN
      !   WRITE(*,*) 'TotalVertflux x3: ',TableC(i)%ScalarC
          CALL AreaVarianceTotalProfile(uCell,Rho,ScalarC,Rho,TotMeanVarUProfile(:,i),&
          & Subdomain=OutputSubdomain,oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
          CALL AreaVarianceTotalProfile(vCell,Rho,ScalarC,Rho,TotMeanVarVProfile(:,i),&
           & Subdomain=OutputSubdomain,oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
          CALL AreaVarianceTotalProfile(wCell,Rho,ScalarC,Rho,TotMeanVarWProfile(:,i),&
          & Subdomain=OutputSubdomain,oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
        END IF
      END IF
      IF (TableC(i)%Variance) THEN
        IF (RhoPos>0) THEN
      !   WRITE(*,*) 'TableVarianceProfile: ',TableC(i)%ScalarC
          CALL AreaVarianceProfile(ScalarC,Rho,ScalarC,Rho,MeanVarianceProfile(:,i),& 
          & Subdomain=OutputSubdomain,oxs=OutputXStart,oxe=OutputXEnd,oys=OutputYStart,oye=OutputYEnd)
        END IF
      END IF
    END DO  

    IF (MyID==0) THEN 
      WRITE(*,*) 'WriteMeanProfile',SIZE(TableC)
      DO i=1,SIZE(TableC)
        IF (TableC(i)%MeanProfile) THEN 
          CALL check(nf90_put_var(ncid,Scalar_varid(i,1),MeanProfile(:,i),start=(/1,TimeStep/)))
        END IF  
        IF (TableC(i)%VertFlux) THEN 
          CALL check(nf90_put_var(ncid,Scalar_varid(i,2),MeanVarUProfile(:,i),start=(/1,TimeStep/)))
          CALL check(nf90_put_var(ncid,Scalar_varid(i,2),MeanVarVProfile(:,i),start=(/1,TimeStep/)))
          CALL check(nf90_put_var(ncid,Scalar_varid(i,2),MeanVarWProfile(:,i),start=(/1,TimeStep/)))
        END IF  
        IF (TableC(i)%TotalVertFlux) THEN
          CALL check(nf90_put_var(ncid,Scalar_varid(i,3),TotMeanVarUProfile(:,i),start=(/1,TimeStep/)))
          CALL check(nf90_put_var(ncid,Scalar_varid(i,3),TotMeanVarVProfile(:,i),start=(/1,TimeStep/)))
          CALL check(nf90_put_var(ncid,Scalar_varid(i,3),TotMeanVarWProfile(:,i),start=(/1,TimeStep/)))
        END IF
        IF (TableC(i)%Variance) THEN
          CALL check(nf90_put_var(ncid,Scalar_varid(i,4),MeanVarianceProfile(:,i),start=(/1,TimeStep/)))
        END IF
      END DO
      CALL check(nf90_inq_varid(ncid,TRIM('time'),rec_varid))
      CALL check(nf90_put_var(ncid,rec_varid,Time,start=(/TimeStep/)))
      CALL check(nf90_sync(ncid))
      CALL check(nf90_close(ncid))
    END IF 
  END IF  
END SUBROUTINE OutputMeanProfile

SUBROUTINE CloseMeanProfile
  CALL check(nf90_close(ncid) ) ! close: save new netCDF dataset
END SUBROUTINE CloseMeanProfile

SUBROUTINE MeanProfileNMLOutput(FileName,Check)
  CHARACTER(*) :: FileName
  LOGICAL :: Check
  CHARACTER(600) :: Line
 
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ProfileOutput')>0) THEN
      WRITE(*,*) Line
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ProfileOutput)
      Check=.TRUE.
      EXIT
    END IF
  END DO
1 CONTINUE 
  CLOSE(InputUnit)
END SUBROUTINE MeanProfileNMLOutput

SUBROUTINE MeanProfileDataTable(FileName)
  CHARACTER(*) :: FileName

  CHARACTER(1000) :: Line
  INTEGER :: i,LenTable 

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*) Line
    IF (INDEX(Line,'#MeanProfileTable')>0) THEN
      READ(InputUnit,*) Lentable 
      ALLOCATE(TableC(LenTable))
      READ(InputUnit,*) Line 
      DO i=1,LenTable
        READ(InputUnit,*) TableC(i)%ScalarC,TableC(i)%Units,TableC(i)%Description,& 
        TableC(i)%MeanProfile,TableC(i)%VertFlux,TableC(i)%TotalVertFlux,TableC(i)%Variance
      ! WRITE(*,*) TableC(i)%ScalarC,TableC(i)%Units,TableC(i)%Description,TableC(i)%MeanProfile,TableC(i)%VertFlux,TableC(i)%TotalVertFlux,TableC(i)%Variance
      END DO  
      EXIT
    END IF
  END DO
  CLOSE(InputUnit) 
END SUBROUTINE MeanProfileDataTable

END MODULE OutputMeanProfile_Mod
