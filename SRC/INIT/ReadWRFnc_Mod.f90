MODULE ReadWRFnc_Mod
  USE Kind_Mod, ONLY : RealKind
  USE NETCDF
  USE Domain_Mod
  USE Parallel_Mod
  USE Names_Mod
  USE Parameter_Mod
  USE Floor_Mod, ONLY : Domain,ib,ibLoc
  USE Profile_Mod
  USE DataType_Mod, ONLY : Vector4Cell_T,LocGlob
  USE Control_Mod, ONLY : &
           thPosEnv,uPosEnv,vPosEnv,wPosEnv,RhoVPosEnv,RhoPosEnv &
          ,Year,Month,Day 

  IMPLICIT NONE

  CHARACTER(80)  :: InputFileNameWRF
  REAL(RealKind) :: TimeIncr

  NAMELIST /ReadWRFncControl/ InputFileNameWRF &
                             ,TimeIncr
                            
CONTAINS 
SUBROUTINE InputWRFnc(FileName)
  CHARACTER(*) :: FileName
  INTEGER      :: InputUnit=20

  OPEN(UNIT=InputUnit,FILE=FileName,STATUS='UNKNOWN')
  READ(InputUnit,NML=ReadWRFncControl)
  CLOSE(InputUnit)
END SUBROUTINE InputWRFnc

SUBROUTINE ReadWRFnc(VecEnv,InputFileNameWRF,Time,VecT)
  TYPE(Vector4Cell_T) :: VecEnv(:)
  REAL(RealKind) :: Time
  TYPE(Vector4Cell_T), OPTIONAL :: VecT(:)
  CHARACTER(*)                 :: InputFileNameWRF
  INTEGER                      :: ncid
  INTEGER                      :: mean_PreC_varid,pert_PreC_varid,pert_PotTemp_varid
  INTEGER                      :: qvC_varid,xlat_varid,xlong_varid,times_varid
  INTEGER                      :: uWindC_varid,vWindC_varid,wWindC_varid,pert_GeoPot_varid,mean_GeoPot_varid
  INTEGER                      :: t_dimid,x_stag_dimid,x_dimid,y_stag_dimid,y_dimid,z_stag_dimid,z_dimid 
  INTEGER                      :: datestr_dimid
  INTEGER                      :: ntime,nzC,nzC_stag,nLatC_stag,nLonC_stag
  INTEGER                      :: izC,iLatC,iLonC
  INTEGER                      :: DateStrLen
  REAL, ALLOCATABLE            :: mean_PreC(:,:,:,:)
  REAL, ALLOCATABLE            :: pert_PreC(:,:,:,:)
  REAL, ALLOCATABLE            :: PreC(:,:,:,:)
  REAL(RealKind), POINTER      :: PreFine(:,:,:,:)
  REAL, ALLOCATABLE            :: pert_PotTempC(:,:,:,:)
  REAL, ALLOCATABLE            :: PotTempC(:,:,:,:)
  REAL(RealKind), POINTER      :: PotTempFine(:,:,:,:)
  REAL(RealKind), POINTER      :: HeightFine(:,:,:,:)
  REAL, ALLOCATABLE            :: QvC(:,:,:,:)
  REAL(RealKind), POINTER      :: QvFine(:,:,:,:)
  REAL, ALLOCATABLE            :: uWindC(:,:,:,:)       !west_east_stag
  REAL(RealKind), POINTER      :: uWindFine(:,:,:,:)
  REAL, ALLOCATABLE            :: vWindC(:,:,:,:)       !south_north_stag 
  REAL(RealKind), POINTER      :: vWindFine(:,:,:,:)
  REAL, ALLOCATABLE            :: wWindC(:,:,:,:)       !bottom_top_stag 
  REAL(RealKind), POINTER      :: wWindFine(:,:,:,:)
  REAL, ALLOCATABLE            :: pert_GeoPot(:,:,:,:)  !bottom_top_stag 
  REAL, ALLOCATABLE            :: mean_GeoPot(:,:,:,:)  !bottom_top_stag 
  REAL, ALLOCATABLE            :: xlat(:,:)
  REAL, ALLOCATABLE            :: xlong(:,:)
  REAL, ALLOCATABLE            :: HeightWRF(:,:,:,:)
  REAL, ALLOCATABLE            :: dLonF(:)
  REAL, ALLOCATABLE            :: dLatF(:)
  REAL, ALLOCATABLE            :: LonFP(:)
  REAL, ALLOCATABLE            :: LatFP(:)
  REAL, ALLOCATABLE            :: dzF(:)
  REAL(RealKind), ALLOCATABLE  :: c(:,:)
  REAL(RealKind), ALLOCATABLE  :: HeightProf(:)
  REAL(RealKind), ALLOCATABLE  :: PreProf(:)
  REAL(RealKind), ALLOCATABLE  :: TempProf(:)
  REAL(RealKind), ALLOCATABLE  :: QvProf(:)
  REAL(RealKind), ALLOCATABLE  :: uWindProf(:)
  REAL(RealKind), ALLOCATABLE  :: vWindProf(:)
  REAL(RealKind), ALLOCATABLE  :: wWindProf(:)
  REAL(RealKind)               :: dLonC,dLatC
  REAL(RealKind)               :: Ttar
  REAL(RealKind)               :: lonC0,latC0
  REAL(RealKind)               :: fLon,fLat,dLonR,dLonL,dLatR,dLatL
  INTEGER                      :: nLonC,nLatC,nLonF,nLatF,iLonF,iLatF
  INTEGER                      :: nzF,izF
  INTEGER                      :: iLonR,iLonL,iLatR,iLatL
  INTEGER                      :: iLonLR,iLatLR
  INTEGER                      :: Hours, Minutes, Seconds
  INTEGER                      :: xtime,itime,y
  CHARACTER, ALLOCATABLE       :: ytime(:,:)

  CHARACTER(19) :: Date
  CHARACTER(19) :: testDate

  CALL check( nf90_open(InputFileNameWRF, nf90_nowrite, ncid)) !open
  CALL check( nf90_inq_dimid(ncid, "Time", t_dimid))
  CALL check( nf90_inquire_dimension(ncid, t_dimid, len=ntime))
  CALL check( nf90_inq_varid(ncid, "Times", times_varid))     
  CALL check( nf90_inq_dimid(ncid, "DateStrLen", datestr_dimid))
  CALL check( nf90_inquire_dimension(ncid, datestr_dimid, len=DateStrLen))
  ALLOCATE(ytime(DateStrLen,ntime))
  CALL check( nf90_get_var(ncid, times_varid, ytime(:,:) , &
                              start = (/ 1, 1 /),  &
                              count = (/ DateStrLen, ntime /)))


 ! 2018-11-29_00:00:00

  WRITE(Date(1:4),'(I4)') Year
  Date(5:5)='-'
  WRITE(Date(6:7),'(I2)') Month
  Date(8:8)='-'
  WRITE(Date(9:10),'(I2)') Day
  Date(11:11)='_'
  Hours=Time/3600.0d0
  IF (Hours<10) THEN
    Date(12:12)='0'
    WRITE(Date(13:13),'(I1)') Hours
  ELSE  
    WRITE(Date(12:13),'(I2)') Hours
  END IF
  Date(14:14)=':'
  Minutes=(Time-Hours*3600)/60
  IF (Minutes<10) THEN
    Date(15:15)='0'
    WRITE(Date(16:16),'(I1)') Minutes
  ELSE  
    WRITE(Date(15:16),'(I2)') Minutes
  END IF
  Date(17:17)=':'
  Seconds=Time-Hours*3600-Minutes*60
  IF (Seconds<10) THEN
    Date(18:18)='0'
    WRITE(Date(19:19),'(I1)') Seconds
  ELSE  
    WRITE(Date(18:19),'(I2)') Seconds
  END IF
   

  DO itime=1,ntime
    DO y=1,DateStrLen
      WRITE(testDate(y:y),'(A1)')ytime(y,itime)
    END DO 
    IF (testDate==Date) THEN
      xtime=itime
    END IF
  END DO
  
  DEALLOCATE(ytime)

  CALL check( nf90_inq_varid(ncid, "P", mean_PreC_varid))     
  CALL check( nf90_inq_varid(ncid, "PB", pert_PreC_varid))     !P+PB=real pressure (Pa)
  CALL check( nf90_inq_varid(ncid, "T", pert_PotTemp_varid))   !T=theta-300 (K)   !absolute Temp=(T+300)*((P+PB)/P00)**(Rdry/cpdry)
  CALL check( nf90_inq_varid(ncid, "QVAPOR", qvC_varid))       !water vapor mixing ratio (kg/kg)
  CALL check( nf90_inq_varid(ncid, "U", uWindC_varid))         !x-wind component (m/s) 
  CALL check( nf90_inq_varid(ncid, "V", vWindC_varid))         !y-wind component (m/s)
  CALL check( nf90_inq_varid(ncid, "W", wWindC_varid))         !z-wind component (m/s) 
  CALL check( nf90_inq_varid(ncid, "PH", pert_GeoPot_varid))   !perturbation geopotential (m2/s2)
  CALL check( nf90_inq_varid(ncid, "PHB", mean_GeoPot_varid))  !base state geopotential (m2/s2)
  CALL check( nf90_inq_dimid(ncid, "bottom_top", z_dimid))
  CALL check( nf90_inq_dimid(ncid, "bottom_top_stag", z_stag_dimid))
  CALL check( nf90_inq_dimid(ncid, "south_north", y_dimid))
  CALL check( nf90_inq_dimid(ncid, "south_north_stag", y_stag_dimid))
  CALL check( nf90_inq_dimid(ncid, "west_east", x_dimid))
  CALL check( nf90_inq_dimid(ncid, "west_east_stag", x_stag_dimid))
  CALL check( nf90_inq_varid(ncid, "XLONG", xlong_varid))  
  CALL check( nf90_inq_varid(ncid, "XLAT", xlat_varid)) 
  CALL check( nf90_inquire_dimension(ncid, z_dimid, len=nzC))
  CALL check( nf90_inquire_dimension(ncid, z_stag_dimid, len=nzC_stag))
  CALL check( nf90_inquire_dimension(ncid, y_dimid, len=nLatC))
  CALL check( nf90_inquire_dimension(ncid, y_stag_dimid, len=nLatC_stag))
  CALL check( nf90_inquire_dimension(ncid, x_dimid, len=nLonC))
  CALL check( nf90_inquire_dimension(ncid, x_stag_dimid, len=nLonC_stag))
  

  ALLOCATE(HeightWRF(nLonC,nLatC,nzC_stag,xtime:xtime))
  ALLOCATE(mean_PreC(nLonC,nLatC,nzC,xtime:xtime))
  ALLOCATE(pert_PreC(nLonC,nLatC,nzC,xtime:xtime))
  ALLOCATE(PreC(nLonC,nLatC,nzC,xtime:xtime))
  ALLOCATE(pert_PotTempC(nLonC,nLatC,nzC,xtime:xtime))
  ALLOCATE(PotTempC(nLonC,nLatC,nzC,xtime:xtime))
  ALLOCATE(QvC(nLonC,nLatC,nzC,xtime:xtime))
  ALLOCATE(uWindC(nLonC_stag,nLatC,nzC,xtime:xtime))
  ALLOCATE(vWindC(nLonC,nLatC_stag,nzC,xtime:xtime))
  ALLOCATE(wWindC(nLonC,nLatC,nzC_stag,xtime:xtime))
  ALLOCATE(pert_GeoPot(nLonC,nLatC,nzC_stag,xtime:xtime))
  ALLOCATE(mean_GeoPot(nLonC,nLatC,nzC_stag,xtime:xtime))
  ALLOCATE(xlat(nLonC,nLatC))
  ALLOCATE(xlong(nLonC,nLatC))


  CALL check( nf90_get_var(ncid, mean_PreC_varid, mean_PreC(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &              !in fortran you need to reverse the order of dimensions
                              count = (/ nLonC, nLatC, nzC, 1 /)))
  CALL check( nf90_get_var(ncid, pert_PreC_varid, pert_PreC(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &             
                              count = (/ nLonC, nLatC, nzC, 1 /)))
  CALL check( nf90_get_var(ncid, pert_PotTemp_varid, pert_PotTempC(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &             
                              count = (/ nLonC, nLatC, nzC, 1 /)))
  CALL check( nf90_get_var(ncid, qvC_varid, QvC(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &             
                              count = (/ nLonC, nLatC, nzC, 1 /)))
  CALL check( nf90_get_var(ncid, uWindC_varid, uWindC(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &             
                              count = (/ nLonC_stag, nLatC, nzC, 1 /)))
  CALL check( nf90_get_var(ncid, vWindC_varid, vWindC(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &             
                              count = (/ nLonC, nLatC_stag, nzC, 1 /)))
  CALL check( nf90_get_var(ncid, wWindC_varid, wWindC(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &             
                              count = (/ nLonC, nLatC, nzC_stag, 1 /)))
  CALL check( nf90_get_var(ncid, pert_GeoPot_varid, pert_GeoPot(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &             
                              count = (/ nLonC, nLatC, nzC_stag, 1 /)))
  CALL check( nf90_get_var(ncid, mean_GeoPot_varid, mean_GeoPot(:,:,:,:) , &
                              start = (/ 1, 1, 1, xtime /),  &             
                              count = (/ nLonC, nLatC, nzC_stag, 1 /)))
  CALL check( nf90_get_var(ncid, xlong_varid, xlong(:,:) , &
                              start = (/ 1, 1 /),  &
                              count = (/ nLonC, nLatC /)))
  CALL check( nf90_get_var(ncid, xlat_varid, xlat(:,:) , &
                              start = (/ 1, 1 /),  &
                              count = (/ nLonC, nLatC /)))




  PotTempC=pert_PotTempC+300.0d0
  PreC=mean_PreC+pert_PreC

  !Wind components at cell center
  DO iLonC=1,nLonC
    DO iLatC=1,nLatC
      DO izC=1,nzC
        uWindC(iLonC,iLatC,izC,xtime)=0.5*(uWindC(iLonC,iLatC,izC,xtime)+uWindC(iLonC+1,iLatC,izC,xtime))
        vWindC(iLonC,iLatC,izC,xtime)=0.5*(vWindC(iLonC,iLatC,izC,xtime)+vWindC(iLonC,iLatC+1,izC,xtime))
        wWindC(iLonC,iLatC,izC,xtime)=0.5*(wWindC(iLonC,iLatC,izC,xtime)+wWindC(iLonC,iLatC,izC+1,xtime))
      END DO
    END DO
  END DO

  !Height WRF at w-points (m)
  HeightWRF=(pert_GeoPot+mean_GeoPot)/Grav

  !Area
  lonC0=xlong(1,1)
  latC0=xlat(1,1)

  DO iLatC=2,nLatC
    DO iLonC=2,nLonC
      dLonC=xlong(iLonC,iLatC)-xlong(iLonC-1,iLatC)
      dLatC=xlat(iLonC,iLatC)-xlat(iLonC,iLatC-1)
    END DO  
  END DO

  dLonC=0.100000d0
  dLatC=0.100000d0

  ALLOCATE(HeightProf(nzC+1))
  ALLOCATE(TempProf(nzC))
  ALLOCATE(PreProf(nzC))
  ALLOCATE(QVProf(nzC))
  ALLOCATE(uWindProf(nzC))
  ALLOCATE(vWindProf(nzC))
  ALLOCATE(wWindProf(nzC))

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)

    nLonF=nx
    nLatF=ny
    nzF=nz

    ALLOCATE(dLonF(nLonF))
    dLonF=dx*180.0d0/Pi
    ALLOCATE(LonFP(0:nLonF+2))
    LonFP(1)=x0*180.0d0/Pi
    DO iLonF=1,nLonF
      LonFP(iLonF+1)=LonFP(iLonF)+dLonF(iLonF)
    END DO  
    LonFP(0)=LonFP(1)-dLonF(1)
    LonFP(nLonF+2)=LonFP(nLonF+1)+dLonF(nLonF)
    
    ALLOCATE(dLatF(nLatF))
    dLatF=dy*180.0d0/Pi
    ALLOCATE(LatFP(0:nLatF+2))
    LatFP(1)=y0*180.0d0/Pi
    DO iLatF=1,nLatF
      LatFP(iLatF+1)=LatFP(iLatF)+dLatF(iLatF)
    END DO  
    LatFP(0)=LatFP(1)-dLatF(1)
    LatFP(nLatF+2)=LatFP(nLatF+1)+dLatF(nLatF)
    
    ALLOCATE(PotTempFine(nLonF,nLatF,nzC,xtime:xtime))
    ALLOCATE(HeightFine(nLonF,nLatF,nzC+1,xtime:xtime))
    ALLOCATE(PreFine(nLonF,nLatF,nzC,xtime:xtime))
    ALLOCATE(QvFine(nLonF,nLatF,nzC,xtime:xtime))
    ALLOCATE(uWindFine(nLonF,nLatF,nzC,xtime:xtime))
    ALLOCATE(vWindFine(nLonF,nLatF,nzC,xtime:xtime))
    ALLOCATE(wWindFine(nLonF,nLatF,nzC,xtime:xtime))
    ALLOCATE(dzF(nzF))
    
    DO iLatC=1,nLatC
      DO iLonC=1,nLonC
        HeightProf(:)=HeightWRF(iLonC,iLatC,:,xtime)
        TempProf(:)=PotTempC(iLonC,iLatC,:,xtime)
        PreProf(:)=PreC(iLonC,iLatC,:,xtime)
        QVProf(:)=QvC(iLonC,iLatC,:,xtime)
      END DO
    END DO

    DO iLatF=1,nLatF
      DO iLonF=1,nLonF
        fLon=(LonFP(iLonF-1)+dLonF(iLonF)/2.0d0-lonC0)/dLonC
        fLat=(LatFP(iLatF-1)+dLatF(iLatF)/2.0d0-latC0)/dLatC
        iLonR=CEILING(fLon) 
        iLonL=iLonR-1 !FLOOR(fLon)
        iLonLR=NINT(fLon)
        IF (iLonLR==iLonR) THEN
          iLonL=iLonL+1
          iLonR=iLonR+1
        END IF
        dLonR=iLonR-0.5d0-fLon
        dLonL=fLon-iLonL+0.5d0

        iLatR=CEILING(fLat)
        iLatL=iLatR-1 !FLOOR(fLat)
        iLatLR=NINT(fLat)
        IF (iLatLR==iLatR) THEN
          iLatR=iLatR+1
          iLatL=iLatL+1
        END IF
        dLatR=iLatR-0.5d0-fLat
        dLatL=fLat-iLatL+0.5d0
        
        Ttar=(1.0d0-dLonL)*(1.0d0-dLatL)*HeightWRF(iLonL,iLatL,nzC+1,xtime) &
            +(1.0d0-dLonR)*(1.0d0-dLatL)*HeightWRF(iLonR,iLatL,nzC+1,xtime) &
            +(1.0d0-dLonL)*(1.0d0-dLatR)*HeightWRF(iLonL,iLatR,nzC+1,xtime) &
            +(1.0d0-dLonR)*(1.0d0-dLatR)*HeightWRF(iLonR,iLatR,nzC+1,xtime)
        HeightFine(iLonF,iLatF,nzC+1,xtime)=Ttar    
        DO izC=1,nzC
          Ttar=(1.0d0-dLonL)*(1.0d0-dLatL)*HeightWRF(iLonL,iLatL,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatL)*HeightWRF(iLonR,iLatL,izC,xtime) &
              +(1.0d0-dLonL)*(1.0d0-dLatR)*HeightWRF(iLonL,iLatR,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatR)*HeightWRF(iLonR,iLatR,izC,xtime)
          HeightFine(iLonF,iLatF,izC,xtime)=Ttar    
          Ttar=(1.0d0-dLonL)*(1.0d0-dLatL)*PotTempC(iLonL,iLatL,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatL)*PotTempC(iLonR,iLatL,izC,xtime) &
              +(1.0d0-dLonL)*(1.0d0-dLatR)*PotTempC(iLonL,iLatR,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatR)*PotTempC(iLonR,iLatR,izC,xtime)
          PotTempFine(iLonF,iLatF,izC,xtime)=Ttar    
          Ttar=(1.0d0-dLonL)*(1.0d0-dLatL)*QvC(iLonL,iLatL,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatL)*QvC(iLonR,iLatL,izC,xtime) &
              +(1.0d0-dLonL)*(1.0d0-dLatR)*QvC(iLonL,iLatR,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatR)*QvC(iLonR,iLatR,izC,xtime)
          QvFine(iLonF,iLatF,izC,xtime)=Ttar    
          Ttar=(1.0d0-dLonL)*(1.0d0-dLatL)*PreC(iLonL,iLatL,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatL)*PreC(iLonR,iLatL,izC,xtime) &
              +(1.0d0-dLonL)*(1.0d0-dLatR)*PreC(iLonL,iLatR,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatR)*PreC(iLonR,iLatR,izC,xtime)
          PreFine(iLonF,iLatF,izC,xtime)=Ttar    
          Ttar=(1.0d0-dLonL)*(1.0d0-dLatL)*uWindC(iLonL,iLatL,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatL)*uWindC(iLonR,iLatL,izC,xtime) &
              +(1.0d0-dLonL)*(1.0d0-dLatR)*uWindC(iLonL,iLatR,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatR)*uWindC(iLonR,iLatR,izC,xtime)
          uWindFine(iLonF,iLatF,izC,xtime)=Ttar    
          Ttar=(1.0d0-dLonL)*(1.0d0-dLatL)*vWindC(iLonL,iLatL,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatL)*vWindC(iLonR,iLatL,izC,xtime) &
              +(1.0d0-dLonL)*(1.0d0-dLatR)*vWindC(iLonL,iLatR,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatR)*vWindC(iLonR,iLatR,izC,xtime)
          vWindFine(iLonF,iLatF,izC,xtime)=Ttar    
          Ttar=(1.0d0-dLonL)*(1.0d0-dLatL)*wWindC(iLonL,iLatL,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatL)*wWindC(iLonR,iLatL,izC,xtime) &
              +(1.0d0-dLonL)*(1.0d0-dLatR)*wWindC(iLonL,iLatR,izC,xtime) &
              +(1.0d0-dLonR)*(1.0d0-dLatR)*wWindC(iLonR,iLatR,izC,xtime)
          wWindFine(iLonF,iLatF,izC,xtime)=Ttar    
        END DO  
    
      END DO
    END DO

    ALLOCATE(zMProf(nzF))
    nzProf=nz
    DO izF=1,nz
      zMProf(izF)=0.5d0*(zP(izF-1)+zP(izF))
    END DO  
    ProfileType='MoistQVPot'
    DO iLatF=1,nLatF
      DO iLonF=1,nLonF
        HeightProf(:)=HeightFine(iLonF,iLatF,:,xtime)
        TempProf(:)=PotTempFine(iLonF,iLatF,:,xtime)
        PreProf(:)=PreFine(iLonF,iLatF,:,xtime)
        QVProf(:)=QvFine(iLonF,iLatF,:,xtime)
        uWindProf(:)=uWindFine(iLonF,iLatF,:,xtime)
        vWindProf(:)=vWindFine(iLonF,iLatF,:,xtime)
        wWindProf(:)=wWindFine(iLonF,iLatF,:,xtime)
        CALL ComputeProfile(c,HeightProf,Pre=PreProf,ThD=TempProf,QV=QVProf,uWind=uWindProf,vWind=vWindProf,wWind=wWindProf)
        IF (ThPosEnv>0) THEN
          VecEnv(ibLoc)%Vec(ThPosEnv)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,3)*c(:,4)
        END IF  
        IF (RhoPosEnv>0) THEN
          VecEnv(ibLoc)%Vec(RhoPosEnv)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,4)
        END IF  
        IF (RhoVPosEnv>0) THEN
          VecEnv(ibLoc)%Vec(RhoVPosEnv)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,5)
        END IF  
        IF (uPosEnv>0) THEN
          VecEnv(ibLoc)%Vec(uPosEnv)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,8)*c(:,4)
        END IF  
        IF (vPosEnv>0) THEN
          VecEnv(ibLoc)%Vec(vPosEnv)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,9)*c(:,4)
        END IF  
        IF (wPosEnv>0) THEN
          VecEnv(ibLoc)%Vec(wPosEnv)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,10)*c(:,4)
        END IF  
        IF (PRESENT(VecT)) THEN
            IF (ThPos>0) THEN
            VecT(ibLoc)%Vec(ThPos)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,3)*c(:,4)
          END IF  
          IF (RhoPos>0) THEN
            VecT(ibLoc)%Vec(RhoPos)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,4)
          END IF  
          IF (RhoVPos>0) THEN
            VecT(ibLoc)%Vec(RhoVPos)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,5)
          END IF  
          IF (uPosL>0) THEN
            VecT(ibLoc)%Vec(uPosL)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,8)*c(:,4)
            VecT(ibLoc)%Vec(uPosR)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,8)*c(:,4)
          END IF  
          IF (vPosL>0) THEN
            VecT(ibLoc)%Vec(vPosL)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,9)*c(:,4)
            VecT(ibLoc)%Vec(vPosR)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,9)*c(:,4)
          END IF  
          IF (wPosL>0) THEN
            VecT(ibLoc)%Vec(wPosL)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,10)*c(:,4)
            VecT(ibLoc)%Vec(wPosR)%c(iLonF+ix0,iLatF+iy0,iz0+1:iz1,1)=c(:,10)*c(:,4)
          END IF  
        END IF
      END DO
    END DO
    DEALLOCATE(dLonF)
    DEALLOCATE(LonFP)
    DEALLOCATE(dLatF)
    DEALLOCATE(LatFP)
    DEALLOCATE(PotTempFine)
    DEALLOCATE(HeightFine)
    DEALLOCATE(PreFine)
    DEALLOCATE(QvFine)
    DEALLOCATE(uWindFine)
    DEALLOCATE(vWindFine)
    DEALLOCATE(wWindFine)
    DEALLOCATE(dzF)
    DEALLOCATE(zMProf)
    DEALLOCATE(c)
  END DO

  DEALLOCATE(HeightProf)
  DEALLOCATE(TempProf)
  DEALLOCATE(PreProf)
  DEALLOCATE(QVProf)
  DEALLOCATE(uWindProf)
  DEALLOCATE(vWindProf)
  DEALLOCATE(wWindProf)

END SUBROUTINE ReadWRFnc


SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS
     
  IF(STATUS /= nf90_noerr) THEN
    PRINT *, TRIM(nf90_strerror(STATUS))
    STOP "Stopped"
  END IF 
END SUBROUTINE check 
END MODULE ReadWRFnc_Mod
