PROGRAM ReadWRFnc 
   USE NETCDF
   USE Profile_Mod

IMPLICIT NONE

INTEGER                     :: ncid
INTEGER                     :: base_PreC_varid,mean_PreC_varid,pert_PreC_varid,pert_PotTemp_varid!,base_Temp_varid
INTEGER                     :: qvC_varid,xlat_varid,xlong_varid
INTEGER                     :: uWindC_varid,vWindC_varid,wWindC_varid,pert_GeoPot_varid,mean_GeoPot_varid
INTEGER                     :: ntime,nzC,nzC_stag,nLatC_stag,nLonC_stag
INTEGER                     :: t_dimid,x_stag_dimid,x_dimid,y_stag_dimid,y_dimid,z_stag_dimid,z_dimid, t,izC,iLatC,iLonC
REAL                        :: base_PreC
REAL, ALLOCATABLE           :: mean_PreC(:,:,:,:) 
REAL, ALLOCATABLE           :: pert_PreC(:,:,:,:) 
REAL, ALLOCATABLE           :: PreC(:,:,:,:) 
REAL, ALLOCATABLE           :: PreFine(:,:,:,:) 
REAL, ALLOCATABLE           :: pert_PotTempC(:,:,:,:) 
REAL, ALLOCATABLE           :: PotTempC(:,:,:,:) 
REAL, ALLOCATABLE           :: PotTempFine(:,:,:,:) 
REAL, ALLOCATABLE           :: QvC(:,:,:,:) 
REAL, ALLOCATABLE           :: QvFine(:,:,:,:) 
REAL, ALLOCATABLE           :: uWindC(:,:,:,:)     !west_east_stag
REAL, ALLOCATABLE           :: uWindFine(:,:,:,:) 
REAL, ALLOCATABLE           :: vWindC(:,:,:,:)     !south_north_stag 
REAL, ALLOCATABLE           :: vWindFine(:,:,:,:) 
REAL, ALLOCATABLE           :: wWindC(:,:,:,:)     !bottom_top_stag 
REAL, ALLOCATABLE           :: pert_GeoPot(:,:,:,:)     !bottom_top_stag 
REAL, ALLOCATABLE           :: mean_GeoPot(:,:,:,:)     !bottom_top_stag 
REAL, ALLOCATABLE           :: xlat(:,:)
REAL, ALLOCATABLE           :: xlong(:,:)
REAL, ALLOCATABLE           :: HeightWRF(:,:,:,:)
REAL, ALLOCATABLE           :: HeightWRF_theta(:,:,:,:)
REAL, ALLOCATABLE           :: HeightFine(:)
REAL, ALLOCATABLE           :: dLonF(:)
REAL, ALLOCATABLE           :: dLatF(:)
REAL, ALLOCATABLE           :: LonFP(:)
REAL, ALLOCATABLE           :: LatFP(:)
REAL, ALLOCATABLE           :: dzF(:)



!REAL(4) :: dLonL=19393.0         !Abweitung WGS84 80°
!REAL(4) :: dLatL=111660.0       
REAL(4) :: g=9.81     
REAL(4) :: dLonC,dLatC
REAL(4) :: dzR,dzL
REAL(4) :: Qvtar,Ttar,Ptar,Utar,Vtar
REAL(4) :: lonC0,latC0,lonF0,latF0,zF0,dLonFEqual,dLatFEqual
REAL(4) :: fLon,fLat,dLonR,dLonL,dLatR,dLatL
REAL(4) :: dLonR_stag,dLonL_stag
REAL(4) :: fz
INTEGER :: nLonC,nLatC,nLonF,nLatF,iLonF,iLatF
INTEGER :: nzF,izF
INTEGER :: iLonR,iLonL,iLatR,iLatL,izR,izL
INTEGER :: iLonLR,iLatLR


REAL(4), ALLOCATABLE :: xLoc(:),yLoc(:),zLoc(:)
REAL(8), ALLOCATABLE :: c(:,:)
REAL(8), ALLOCATABLE :: HeightProf(:)
REAL(8), ALLOCATABLE :: PreProf(:)
REAL(8), ALLOCATABLE :: TempProf(:)
REAL(8), ALLOCATABLE :: QvProf(:)
INTEGER :: i
CHARACTER (LEN = *), PARAMETER :: X_NAME = "xc"
CHARACTER (LEN = *), PARAMETER :: Y_NAME = "yc"
CHARACTER (LEN = *), PARAMETER :: Z_NAME = "lev"
CHARACTER (LEN = *), PARAMETER :: X_LONGNAME = "x-coordinate in Cartesian system"
CHARACTER (LEN = *), PARAMETER :: Y_LONGNAME = "y-coordinate in Cartesian system"
CHARACTER (LEN = *), PARAMETER :: Z_LONGNAME = "Height at midpoints"
CHARACTER (LEN = *), PARAMETER :: X_UNITS = "degrees_east"
CHARACTER (LEN = *), PARAMETER :: Y_UNITS = "degrees_north"
CHARACTER (LEN = *), PARAMETER :: Z_UNITS = "m"
CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
INTEGER :: x_varid, y_varid, z_varid,pottempfine_varid,pottemp_varid,prefine_varid,pre_varid
INTEGER :: qvCfine_varid
INTEGER :: dimids(3),count(3),start(3)
CHARACTER (LEN = *), PARAMETER :: POTTEMPFINE_NAME = "PotTempFine"
CHARACTER (LEN = *), PARAMETER :: POTTEMPFINE_UNITS = "K"
CHARACTER (LEN = *), PARAMETER :: PREFINE_NAME = "PreFine"
CHARACTER (LEN = *), PARAMETER :: PREFINE_UNITS = "Pa"
CHARACTER (LEN = *), PARAMETER :: POTTEMP_NAME = "PotTempCoarse"
CHARACTER (LEN = *), PARAMETER :: POTTEMP_UNITS = "K"
CHARACTER (LEN = *), PARAMETER :: PRE_NAME = "PreCoarse"
CHARACTER (LEN = *), PARAMETER :: PRE_UNITS = "Pa"
CHARACTER (LEN = *), PARAMETER :: QVFINE_NAME = "WVMRFine"
CHARACTER (LEN = *), PARAMETER :: QVFINE_UNITS = "kg/kg"
CHARACTER (LEN = *), PARAMETER :: QV_NAME = "WVMRCoarse"
CHARACTER (LEN = *), PARAMETER :: QV_UNITS = "kg/kg"

CALL check( nf90_open("/Users/urbanneck/WRFurbanneck/wrfout_d01_2018-11-29_00:00:00", nf90_nowrite, ncid))     !open
CALL check( nf90_inq_varid(ncid, "P00", base_PreC_varid))    !1000 hPa 
CALL check( nf90_inq_varid(ncid, "P", mean_PreC_varid))     
CALL check( nf90_inq_varid(ncid, "PB", pert_PreC_varid))     ! P+PB=real pressure (Pa)
!CALL check( nf90_inq_varid(ncid, "T00", base_Temp_varid))   !290 K   
!The base state value of 290 should not be changed. It's a parameter in WRF's hypothetical mean temperature profile. To get total potential temperature, it is equal to T + 300 (T is the perturbation potential temperature). Then use Poisson's equation to calculate temperature from potential temperature and pressure (pressure is also split into base and perturbation values). 
CALL check( nf90_inq_varid(ncid, "T", pert_PotTemp_varid))  !T=theta-300 (K)   !absolute Temp=(T+300)*((P+PB)/P00)**(Rdry/cpdry)
CALL check( nf90_inq_varid(ncid, "QVAPOR", qvC_varid))       !water vapor mixing ratio (kg/kg)
CALL check( nf90_inq_varid(ncid, "U", uWindC_varid))         !x-wind component (m/s) 
CALL check( nf90_inq_varid(ncid, "V", vWindC_varid))         !y-wind component (m/s)
CALL check( nf90_inq_varid(ncid, "W", wWindC_varid))         !z-wind component (m/s) 
CALL check( nf90_inq_varid(ncid, "PH", pert_GeoPot_varid))         !perturbation geopotential (m2/s2)
CALL check( nf90_inq_varid(ncid, "PHB", mean_GeoPot_varid))         !base state geopotential (m2/s2)
CALL check( nf90_inq_dimid(ncid, "Time", t_dimid))
CALL check( nf90_inq_dimid(ncid, "bottom_top", z_dimid))
CALL check( nf90_inq_dimid(ncid, "bottom_top_stag", z_stag_dimid))
CALL check( nf90_inq_dimid(ncid, "south_north", y_dimid))
CALL check( nf90_inq_dimid(ncid, "south_north_stag", y_stag_dimid))
CALL check( nf90_inq_dimid(ncid, "west_east", x_dimid))
CALL check( nf90_inq_dimid(ncid, "west_east_stag", x_stag_dimid))
CALL check( nf90_inq_varid(ncid, "XLONG", xlong_varid))  
CALL check( nf90_inq_varid(ncid, "XLAT", xlat_varid)) 
CALL check( nf90_inquire_dimension(ncid, t_dimid, len=ntime))
CALL check( nf90_inquire_dimension(ncid, z_dimid, len=nzC))
CALL check( nf90_inquire_dimension(ncid, z_stag_dimid, len=nzC_stag))
CALL check( nf90_inquire_dimension(ncid, y_dimid, len=nLatC))
CALL check( nf90_inquire_dimension(ncid, y_stag_dimid, len=nLatC_stag))
CALL check( nf90_inquire_dimension(ncid, x_dimid, len=nLonC))
CALL check( nf90_inquire_dimension(ncid, x_stag_dimid, len=nLonC_stag))


ALLOCATE(HeightWRF(nLonC,nLatC,nzC_stag,ntime))
ALLOCATE(HeightWRF_theta(nLonC,nLatC,nzC_stag,ntime))
ALLOCATE(mean_PreC(nLonC,nLatC,nzC,ntime))
ALLOCATE(pert_PreC(nLonC,nLatC,nzC,ntime))
ALLOCATE(PreC(nLonC,nLatC,nzC,ntime))
ALLOCATE(pert_PotTempC(nLonC,nLatC,nzC,ntime))
ALLOCATE(PotTempC(nLonC,nLatC,nzC,ntime))
ALLOCATE(QvC(nLonC,nLatC,nzC,ntime))
ALLOCATE(uWindC(nLonC_stag,nLatC,nzC,ntime))
ALLOCATE(vWindC(nLonC,nLatC_stag,nzC,ntime))
ALLOCATE(wWindC(nLonC,nLatC,nzC_stag,ntime))
ALLOCATE(pert_GeoPot(nLonC,nLatC,nzC_stag,ntime))
ALLOCATE(mean_GeoPot(nLonC,nLatC,nzC_stag,ntime))
ALLOCATE(xlat(nLonC,nLatC))
ALLOCATE(xlong(nLonC,nLatC))


CALL check( nf90_get_var(ncid, base_PreC_varid, base_PreC))
CALL check( nf90_get_var(ncid, mean_PreC_varid, mean_PreC(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &              !in fortran you need to reverse the order of dimensions
                              count = (/ nLonC, nLatC, nzC, ntime /)))
CALL check( nf90_get_var(ncid, pert_PreC_varid, pert_PreC(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &             
                              count = (/ nLonC, nLatC, nzC, ntime /)))
CALL check( nf90_get_var(ncid, pert_PotTemp_varid, pert_PotTempC(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &             
                              count = (/ nLonC, nLatC, nzC, ntime /)))
CALL check( nf90_get_var(ncid, qvC_varid, QvC(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &             
                              count = (/ nLonC, nLatC, nzC, ntime /)))
CALL check( nf90_get_var(ncid, uWindC_varid, uWindC(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &             
                              count = (/ nLonC_stag, nLatC, nzC, ntime /)))
CALL check( nf90_get_var(ncid, vWindC_varid, vWindC(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &             
                              count = (/ nLonC, nLatC_stag, nzC, ntime /)))
CALL check( nf90_get_var(ncid, wWindC_varid, wWindC(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &             
                              count = (/ nLonC, nLatC, nzC_stag, ntime /)))
CALL check( nf90_get_var(ncid, pert_GeoPot_varid, pert_GeoPot(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &             
                              count = (/ nLonC, nLatC, nzC_stag, ntime /)))
CALL check( nf90_get_var(ncid, mean_GeoPot_varid, mean_GeoPot(:,:,:,:) , &
                              start = (/ 1, 1, 1, 1 /),  &             
                              count = (/ nLonC, nLatC, nzC_stag, ntime /)))
CALL check( nf90_get_var(ncid, xlong_varid, xlong(:,:) , &
                              start = (/ 1, 1 /),  &
                              count = (/ nLonC, nLatC /)))
CALL check( nf90_get_var(ncid, xlat_varid, xlat(:,:) , &
                              start = (/ 1, 1 /),  &
                              count = (/ nLonC, nLatC /)))



WRITE(*,*) 'nzC_stag',nzC_stag
WRITE(*,*) 'base pressure',base_PreC



PotTempC=pert_PotTempC+300.0d0
PreC=mean_PreC+pert_PreC

!Height WRF at midpoints (m)
DO iLonC=1,nLonC
 DO iLatC=1,nLatC
  DO izC=1,nzC_stag-1
   DO t=1,ntime

    HeightWRF_theta(iLonC,iLatC,izC,t) =0.5*(mean_GeoPot(iLonC,iLatC,izC,t)+pert_GeoPot(iLonC,iLatC,izC,t) &
                                   +mean_GeoPot(iLonC,iLatC,izC+1,t)+pert_GeoPot(iLonC,iLatC,izC+1,t))/g

   END DO
  END DO
 END DO
END DO

!Height WRF at w-points (m)
HeightWRF=(pert_GeoPot+mean_GeoPot)/g


!Area
lonC0=xlong(1,1)
latC0=xlat(1,1)
lonF0=5.7d0
latF0=55.0d0


dLonC=0.100000d0
dLatC=0.100000d0
dLonFEqual=0.05d0                      
dLatFEqual=0.05d0               


nLonF=10
nLatF=9
nzF=10


write(*,*)'lonC0',lonC0,'latC0',latC0
write(*,*)'lonF0',lonF0,'latF0',latF0
write(*,*)'dLonC',dLonC,'dLatC',dLatC
write(*,*)'dLonFEqual',dLonFEqual,'dLatFEqual',dLatFEqual
write(*,*)'nLonC',nLonC,'nLatC',nLatC
write(*,*)'nLonF',nLonF,'nLatF',nLatF


ALLOCATE(dLonF(nLonF))
dLonF=dLonFEqual
ALLOCATE(LonFP(0:nLonF+2))
LonFP(1)=lonF0
DO iLonF=1,nLonF
  LonFP(iLonF+1)=LonFP(iLonF)+dLonF(iLonF)
END DO  
LonFP(0)=LonFP(1)-dLonF(1)
LonFP(nLonF+2)=LonFP(nLonF+1)+dLonF(nLonF)

ALLOCATE(dLatF(nLatF))
dLatF=dLatFEqual
ALLOCATE(LatFP(0:nLatF+2))
LatFP(1)=latF0
DO iLatF=1,nLatF
  LatFP(iLatF+1)=LatFP(iLatF)+dLatF(iLatF)
END DO  
LatFP(0)=LatFP(1)-dLatF(1)
LatFP(nLatF+2)=LatFP(nLatF+1)+dLatF(nLatF)

ALLOCATE(PotTempFine(nLonF,nLatF,nzF,ntime))
ALLOCATE(PreFine(nLonF,nLatF,nzF,ntime))
ALLOCATE(QvFine(nLonF,nLatF,nzF,ntime))
ALLOCATE(uWindFine(nLonF,nLatF,nzF,ntime))
ALLOCATE(vWindFine(nLonF,nLatF,nzF,ntime))
ALLOCATE(HeightFine(0:nzF))
ALLOCATE(dzF(nzF))

  WRITE(*,*) 'Höhe 1'
HeightFine(0)=600.0d0
dzF=50.0d0
DO izF=1,nzF
  HeightFine(izF)=HeightFine(izF-1)+dzF(izF)
END DO
  nzProf=nzF
  ALLOCATE(zPProf(0:nzF))
  ALLOCATE(dzProf(nzF))
  ALLOCATE(zMProf(nzF))
  dzProf=dzF
  zPProf=HeightFine
  WRITE(*,*) 'Höhe 2'
  DO izF=1,nzProf
    zMProf(izF)=0.5d0*(zPProf(izF-1)+zPProf(izF))
  END DO  
  ProfileType='MoistQVPot'
  ALLOCATE(HeightProf(nzC))
  ALLOCATE(TempProf(nzC))
  ALLOCATE(PreProf(nzC))
  ALLOCATE(QVProf(nzC))

DO iLatC=1,nLatC
  DO iLonC=1,nLonC
    HeightProf(:)=HeightWRF(iLonC,iLatC,:,1)
    TempProf(:)=PotTempC(iLonC,iLatC,:,1)
    PreProf(:)=PreC(iLonC,iLatC,:,1)
    QVProf(:)=QvC(iLonC,iLatC,:,1)
!   DO izC=1,nzC
!     WRITE(*,*) HeightWRF(iLonC,iLatC,izC,1),PreC(iLonC,iLatC,izC,1),PotTempC(iLonC,iLatC,izC,1),QVC(iLonC,iLatC,izC,1)
!   END DO  
    WRITE(*,*) 'iLatC',iLatC,'iLonC',iLonC
!   CALL ComputeProfile(c,HeightProf,Pre=PreProf,ThD=TempProf,QV=QVProf)
!   DO izF=1,nzF
!     WRITE(*,*) 'c',zMProf(izF),c(izF,:)
!   END DO  
  END DO
END DO
!STOP 'Ende '
DO iLatF=1,nLatF
  DO iLonF=1,nLonF
    fLon=(LonFP(iLonF-1)+dLonF(iLonF)/2.0d0-lonC0)/dLonC
    fLat=(LatFP(iLatF-1)+dLatF(iLatF)/2.0d0-latC0)/dLatC
    iLonR=CEILING(fLon)
    iLonL=FLOOR(fLon)
    iLonLR=NINT(fLon)
    IF (iLonLR==iLonR) THEN
      iLonL=iLonL+1
      iLonR=iLonR+1
    END IF
    dLonR=iLonR-0.5d0-fLon
    dLonL=fLon-iLonL+0.5d0
    !lon staggered, use for u-wind
    
    iLatR=CEILING(fLat)
    iLatL=FLOOR(fLat)
    iLatLR=NINT(fLat)
    IF (iLatLR==iLatR) THEN
      iLatR=iLatR+1
      iLatL=iLatL+1
    END IF
    dLatR=iLatR-0.5d0-fLat
    dLatL=fLat-iLatL+0.5d0
    !lat staggered, use for v-wind

    WRITE(*,*) 'iLonL',iLonL,'iLatL',iLatL
    DO izF=1,nzF
      Ttar=0.0d0
      Ptar=0.0d0
      Qvtar=0.0d0
      Utar=0.0d0
      Vtar=0.0d0
      fz=HeightFine(izF-1)+dzF(izF)/2.0d0
      DO izC=2,nzC+1
        IF (fz<0.5d0*(HeightWRF(iLonL,iLatL,izC-1,1)+HeightWRF(iLonL,iLatL,izC,1))) THEN
          izL=izC-2
          izR=izC-1
          !Interpolation zwischen izC-1 und izC-2
          dzR=0.5d0*(HeightWRF(iLonL,iLatL,izR,1)+HeightWRF(iLonL,iLatL,izR+1,1))-fz
          dzL=fz-0.5d0*(HeightWRF(iLonL,iLatL,izL,1)+HeightWRF(iLonL,iLatL,izR,1))
          dzL=dzL/(0.5d0*(HeightWRF(iLonL,iLatL,izR+1,1)-HeightWRF(iLonL,iLatL,izL,1)))
          dzR=dzR/(0.5d0*(HeightWRF(iLonL,iLatL,izR+1,1)-HeightWRF(iLonL,iLatL,izL,1)))
          Ttar=Ttar+(1.0d0-dLonL)*(1.0d0-dLatL)*(1.0d0-dzR)*PotTempC(iLonL,iLatL,izR,1)
          Ttar=Ttar+(1.0d0-dLonL)*(1.0d0-dLatL)*(1.0d0-dzL)*PotTempC(iLonL,iLatL,izL,1)
          Ptar=Ptar+(1.0d0-dLonL)*(1.0d0-dLatL)*(1.0d0-dzR)*PreC(iLonL,iLatL,izR,1)
          Ptar=Ptar+(1.0d0-dLonL)*(1.0d0-dLatL)*(1.0d0-dzL)*PreC(iLonL,iLatL,izL,1)
          Qvtar=Qvtar+(1.0d0-dLonL)*(1.0d0-dLatL)*(1.0d0-dzR)*QvC(iLonL,iLatL,izR,1)
          Qvtar=Qvtar+(1.0d0-dLonL)*(1.0d0-dLatL)*(1.0d0-dzL)*QvC(iLonL,iLatL,izL,1)
          Utar=Utar+(1.0d0-dLonL)*(1.0d0-dLatL)*(1.0d0-dzR)*uWindC(iLonL,iLatL,izR,1)
          Utar=Utar+(1.0d0-dLonL)*(1.0d0-dLatL)*(1.0d0-dzL)*uWindC(iLonL,iLatL,izL,1)
         ! WRITE(*,*) 'Coarse LL',PreC(iLonL,iLatL,izL,1),PreC(iLonL,iLatL,izR,1)
          EXIT
        ELSE IF (fz<0.5d0*HeightWRF(iLonL,iLatL,izC-1,1)  ) THEN
          Ttar=Ttar+(1.0d0-dLonL)*(1.0d0-dLatL)*PotTempC(iLonL,iLatL,izC-1,1)
          Ptar=Ptar+(1.0d0-dLonL)*(1.0d0-dLatL)*PreC(iLonL,iLatL,izC-1,1)
          Qvtar=Qvtar+(1.0d0-dLonL)*(1.0d0-dLatL)*QvC(iLonL,iLatL,izC-1,1)
          EXIT
        END IF  
      END DO  
      DO izC=2,nzC+1
        IF (fz<0.5d0*(HeightWRF(iLonR,iLatL,izC-1,1)+HeightWRF(iLonR,iLatL,izC,1))) THEN
          izL=izC-2
          izR=izC-1
          !Interpolation zwischen izR und izL
          dzR=0.5d0*(HeightWRF(iLonR,iLatL,izR,1)+HeightWRF(iLonR,iLatL,izR+1,1))-fz
          dzL=fz-0.5d0*(HeightWRF(iLonR,iLatL,izL,1)+HeightWRF(iLonR,iLatL,izR,1))
          dzL=dzL/(0.5d0*(HeightWRF(iLonR,iLatL,izR+1,1)-HeightWRF(iLonR,iLatL,izL,1)))
          dzR=dzR/(0.5d0*(HeightWRF(iLonR,iLatL,izR+1,1)-HeightWRF(iLonR,iLatL,izL,1)))
          Ttar=Ttar+(1.0d0-dLonR)*(1.0d0-dLatL)*(1.0d0-dzR)*PotTempC(iLonR,iLatL,izR,1)
          Ttar=Ttar+(1.0d0-dLonR)*(1.0d0-dLatL)*(1.0d0-dzL)*PotTempC(iLonR,iLatL,izL,1)
          Ptar=Ptar+(1.0d0-dLonR)*(1.0d0-dLatL)*(1.0d0-dzR)*PreC(iLonR,iLatL,izR,1)
          Ptar=Ptar+(1.0d0-dLonR)*(1.0d0-dLatL)*(1.0d0-dzL)*PreC(iLonR,iLatL,izL,1)
          Qvtar=Qvtar+(1.0d0-dLonR)*(1.0d0-dLatL)*(1.0d0-dzR)*QvC(iLonR,iLatL,izR,1)
          Qvtar=Qvtar+(1.0d0-dLonR)*(1.0d0-dLatL)*(1.0d0-dzL)*QvC(iLonR,iLatL,izL,1)
         ! WRITE(*,*) 'Coarse RL',PreC(iLonR,iLatL,izL,1),PreC(iLonR,iLatL,izR,1)
          EXIT
        ELSE IF (fz<0.5d0*HeightWRF(iLonR,iLatL,izC-1,1)  ) THEN
          Ttar=Ttar+(1.0d0-dLonR)*(1.0d0-dLatL)*PotTempC(iLonR,iLatL,izC-1,1)
          Ptar=Ptar+(1.0d0-dLonR)*(1.0d0-dLatL)*PreC(iLonR,iLatL,izC-1,1)
          Qvtar=Qvtar+(1.0d0-dLonR)*(1.0d0-dLatL)*QvC(iLonR,iLatL,izC-1,1)
          EXIT
        END IF  
      END DO  
      DO izC=2,nzC+1
        IF (fz<0.5d0*(HeightWRF(iLonL,iLatR,izC-1,1)+HeightWRF(iLonL,iLatR,izC,1))) THEN
          izL=izC-2
          izR=izC-1
          !Interpolation zwischen izR und izL
          dzR=0.5d0*(HeightWRF(iLonL,iLatR,izR,1)+HeightWRF(iLonL,iLatR,izR+1,1))-fz
          dzL=fz-0.5d0*(HeightWRF(iLonL,iLatR,izL,1)+HeightWRF(iLonL,iLatR,izR,1))
          dzL=dzL/(0.5d0*(HeightWRF(iLonL,iLatR,izR+1,1)-HeightWRF(iLonL,iLatR,izL,1)))
          dzR=dzR/(0.5d0*(HeightWRF(iLonL,iLatR,izR+1,1)-HeightWRF(iLonL,iLatR,izL,1)))
          Ttar=Ttar+(1.0d0-dLonL)*(1.0d0-dLatR)*(1.0d0-dzR)*PotTempC(iLonL,iLatR,izR,1)
          Ttar=Ttar+(1.0d0-dLonL)*(1.0d0-dLatR)*(1.0d0-dzL)*PotTempC(iLonL,iLatR,izL,1)
          Ptar=Ptar+(1.0d0-dLonL)*(1.0d0-dLatR)*(1.0d0-dzR)*PreC(iLonL,iLatR,izR,1)
          Ptar=Ptar+(1.0d0-dLonL)*(1.0d0-dLatR)*(1.0d0-dzL)*PreC(iLonL,iLatR,izL,1)
          Qvtar=Qvtar+(1.0d0-dLonL)*(1.0d0-dLatR)*(1.0d0-dzR)*QvC(iLonL,iLatR,izR,1)
          Qvtar=Qvtar+(1.0d0-dLonL)*(1.0d0-dLatR)*(1.0d0-dzL)*QvC(iLonL,iLatR,izL,1)
         ! WRITE(*,*) 'Coarse LR',PreC(iLonL,iLatR,izL,1),PreC(iLonL,iLatR,izR,1)
          EXIT
        ELSE IF (fz<0.5d0*HeightWRF(iLonL,iLatR,izC-1,1)  ) THEN
          Ttar=Ttar+(1.0d0-dLonL)*(1.0d0-dLatR)*PotTempC(iLonL,iLatR,izC-1,1)
          Ptar=Ptar+(1.0d0-dLonL)*(1.0d0-dLatR)*PreC(iLonL,iLatR,izC-1,1)
          Qvtar=Qvtar+(1.0d0-dLonL)*(1.0d0-dLatR)*QvC(iLonL,iLatR,izC-1,1)
          EXIT
        END IF  
      END DO  
      DO izC=2,nzC+1
        IF (fz<0.5d0*(HeightWRF(iLonR,iLatR,izC-1,1)+HeightWRF(iLonR,iLatR,izC,1))) THEN
          izL=izC-2
          izR=izC-1
          !Interpolation zwischen izR und izL
          dzR=0.5d0*(HeightWRF(iLonR,iLatR,izR,1)+HeightWRF(iLonR,iLatR,izR+1,1))-fz
          dzL=fz-0.5d0*(HeightWRF(iLonR,iLatR,izL,1)+HeightWRF(iLonR,iLatR,izR,1))
          dzL=dzL/(0.5d0*(HeightWRF(iLonR,iLatR,izR+1,1)-HeightWRF(iLonR,iLatR,izL,1)))
          dzR=dzR/(0.5d0*(HeightWRF(iLonR,iLatR,izR+1,1)-HeightWRF(iLonR,iLatR,izL,1)))
          Ttar=Ttar+(1.0d0-dLonR)*(1.0d0-dLatR)*(1.0d0-dzR)*PotTempC(iLonR,iLatR,izR,1)
          Ttar=Ttar+(1.0d0-dLonR)*(1.0d0-dLatR)*(1.0d0-dzL)*PotTempC(iLonR,iLatR,izL,1)
          Ptar=Ptar+(1.0d0-dLonR)*(1.0d0-dLatR)*(1.0d0-dzR)*PreC(iLonR,iLatR,izR,1)
          Ptar=Ptar+(1.0d0-dLonR)*(1.0d0-dLatR)*(1.0d0-dzL)*PreC(iLonR,iLatR,izL,1)
          Qvtar=Qvtar+(1.0d0-dLonR)*(1.0d0-dLatR)*(1.0d0-dzR)*QvC(iLonR,iLatR,izR,1)
          Qvtar=Qvtar+(1.0d0-dLonR)*(1.0d0-dLatR)*(1.0d0-dzL)*QvC(iLonR,iLatR,izL,1)
         ! WRITE(*,*) 'Coarse RR',PreC(iLonR,iLatR,izL,1),PreC(iLonR,iLatR,izR,1)
          EXIT
        ELSE IF (fz<0.5d0*HeightWRF(iLonR,iLatR,izC-1,1)  ) THEN
          Ttar=Ttar+(1.0d0-dLonR)*(1.0d0-dLatR)*PotTempC(iLonR,iLatR,izC-1,1)
          Ptar=Ptar+(1.0d0-dLonR)*(1.0d0-dLatR)*PreC(iLonR,iLatR,izC-1,1)
          Qvtar=Qvtar+(1.0d0-dLonR)*(1.0d0-dLatR)*QvC(iLonR,iLatR,izC-1,1)
          EXIT
        END IF  
      END DO  
      PotTempFine(iLonF,iLatF,izF,1)=Ttar
     ! WRITE(*,*) 'Fine',iLonF,iLatF,izF,Ptar
      PreFine(iLonF,iLatF,izF,1)=Ptar
      QvFine(iLonF,iLatF,izF,1)=Qvtar
      uWindFine(iLonF,iLatF,izF,1)=Utar
    END DO  
  END DO
END DO
  DEALLOCATE(HeightProf)
  DEALLOCATE(TempProf)
  DEALLOCATE(PreProf)
  DEALLOCATE(QVProf)
  ALLOCATE(HeightProf(nzF))
  ALLOCATE(TempProf(nzF))
  ALLOCATE(PreProf(nzF))
  ALLOCATE(QVProf(nzF))
 DO iLatF=1,nLatF
   DO iLonF=1,nLonF
!    HeightProf(:)=HeightFine(:)
     HeightProf(:)=zMProf(:)
     TempProf(:)=PotTempFine(iLonF,iLatF,:,1)
     PreProf(:)=PreFine(iLonF,iLatF,:,1)
     QVProf(:)=QvFine(iLonF,iLatF,:,1)
!    CALL ComputeProfile(c,HeightProf,Pre=PreProf,Temp=TempProf,QV=QVProf)
     WRITE(*,*) 'ComputeProfile',iLatF,iLonF
     CALL ComputeProfile(c,HeightProf,Pre=PreProf,ThD=TempProf,QV=QVProf)
     PotTempFine(iLonF,iLatF,:,1)=c(:,3)
!    VecEnv(1)%Vec(ThPos)=c(iLonF,iLatF,:,1)=c(:,3)     
!    DO izF=1,nzF
!      WRITE(*,*) 'c',c(izF,:)
!    END DO  
!    STOP
  END DO
END DO




!NetCDF FINE Output
ALLOCATE(xLoc(1:nLonF))
ALLOCATE(yLoc(1:nLatF))
ALLOCATE(zLoc(1:nzF))

xLoc(1)=LonF0
yLoc(1)=LatF0

DO i=2,nLonF
  xLoc(i)=xLoc(i-1)+dLonFEqual
END DO
 
DO i=2,nLatF
  yLoc(i)=yLoc(i-1)+dLatFEqual
END DO
 
DO izF=1,nzF
  zLoc(izF)=HeightFine(izF-1)+dzF(izF)
END DO

!Create .nc
CALL check( nf90_create('PotTempPreQvFine.nc', nf90_clobber, ncid) )

!Define the dimensions.
CALL check( nf90_def_dim(ncid, TRIM(X_NAME), nLonF, x_dimid) )
CALL check( nf90_def_dim(ncid, TRIM(Y_NAME), nLatF, y_dimid) )
CALL check( nf90_def_dim(ncid, Z_NAME, nzF, z_dimid) )

!Define the Coordinate variables.
CALL check( nf90_def_var(ncid, TRIM(X_NAME), NF90_DOUBLE, x_dimid, x_varid) )
CALL check( nf90_def_var(ncid, TRIM(Y_NAME), NF90_DOUBLE, y_dimid, y_varid) )
CALL check( nf90_def_var(ncid, Z_NAME, NF90_DOUBLE, z_dimid, z_varid) )

!Assign units attributes to coordinate variables.
CALL check( nf90_put_att(ncid, x_varid,"long_name",TRIM(X_LONGNAME)) )
CALL check( nf90_put_att(ncid, x_varid, UNITS, TRIM(X_UNITS)) )
CALL check( nf90_put_att(ncid, y_varid,"long_name",TRIM(Y_LONGNAME)) )
CALL check( nf90_put_att(ncid, y_varid, UNITS, TRIM(Y_UNITS)) )
CALL check( nf90_put_att(ncid, z_varid,"long_name",Z_LONGNAME) )
CALL check( nf90_put_att(ncid, z_varid, UNITS, Z_UNITS) )

!The dimids array is used to pass the dimids of the dimensions of
!the netCDF variables. Both of the netCDF variables we are creating
!share the same four dimensions. In Fortran, the unlimited
!dimension must come last on the list of dimids.
dimids = (/ x_dimid, y_dimid, z_dimid/)

!Define the netCDF variables for the pressure and temperature data.
  CALL check( nf90_def_var(ncid, POTTEMPFINE_NAME,    NF90_REAL, dimids, pottempfine_varid   ) )
  CALL check( nf90_def_var(ncid, PREFINE_NAME,    NF90_REAL, dimids, prefine_varid   ) )
  CALL check( nf90_def_var(ncid, QVFINE_NAME,    NF90_REAL, dimids, qvCfine_varid   ) )

!Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, pottempfine_varid,    UNITS, POTTEMPFINE_UNITS   ) )
    CALL check( nf90_put_att(ncid, pottempfine_varid,"description", "Pot. Temperature interpolated"            ) )
    CALL check( nf90_put_att(ncid, prefine_varid,    UNITS, PREFINE_UNITS   ) )
    CALL check( nf90_put_att(ncid, prefine_varid,"description", "PreCssure interpolated"            ) )
    CALL check( nf90_put_att(ncid, qvCfine_varid,    UNITS, QVFINE_UNITS   ) )
    CALL check( nf90_put_att(ncid, qvCfine_varid,"description", "WVMR interpolated"            ) )

!End define mode.
CALL check( nf90_enddef(ncid) )

!Write the coordinate variable data. This will put the latitudes
!and longitudes of our data grid into the netCDF file.
 DO i=1,nLonF
   CALL check( nf90_put_var(ncid, x_varid, xLoc) )
 END DO
 DO i=1,nLatF
   CALL check( nf90_put_var(ncid, y_varid, yLoc) )
 END DO
 DO i=1,nzF
   CALL check( nf90_put_var(ncid, z_varid, zLoc) )
END DO

! These settings tell netcdf to write one timestep of data.
! count = (/ nLonF, nLatF, nzF/)
! start = (/ 1, 1, 1/)
!Write data.
CALL check( nf90_put_var(ncid,pottempfine_varid,PotTempFine(1:nlonF,1:nLatF,1:nzF,1)))!,start=start,count=count) )
CALL check( nf90_put_var(ncid,prefine_varid,PreFine(1:nlonF,1:nLatF,1:nzF,1)))!,start=start,count=count) )
CALL check( nf90_put_var(ncid,qvCfine_varid,QvFine(1:nlonF,1:nLatF,1:nzF,1)))!,start=start,count=count) )



! Close the file. This causes netCDF to flush all buffers and make
! sure your data are REALly written to disk.
CALL check( nf90_close(ncid) )

DEALLOCATE(xLoc)
DEALLOCATE(yLoc)
DEALLOCATE(zLoc)



!NetCDF small COARSE Output
ALLOCATE(xLoc(1:7))
ALLOCATE(yLoc(94:99))
ALLOCATE(zLoc(1:nzC))

xLoc(1)=LonC0
yLoc(94)=LatC0+93*dLonC

DO i=2,7
  xLoc(i)=xLoc(i-1)+dLonC
END DO

DO i=95,99
  yLoc(i)=yLoc(i-1)+dLatC
END DO

DO i=1,nzC
  zLoc(i)=HeightWRF_theta(1,85,i,1) !zLoc(i-1)+1
END DO

!Create .nc
CALL check( nf90_create('PotTempPreQvCoarse_small.nc', nf90_clobber, ncid) )

!Define the dimensions.
CALL check( nf90_def_dim(ncid, TRIM(X_NAME), 7, x_dimid) )
CALL check( nf90_def_dim(ncid, TRIM(Y_NAME), 6, y_dimid) )
CALL check( nf90_def_dim(ncid, Z_NAME, nzC, z_dimid) )

!Define the Coordinate variables.
CALL check( nf90_def_var(ncid, TRIM(X_NAME), NF90_DOUBLE, x_dimid, x_varid) )
CALL check( nf90_def_var(ncid, TRIM(Y_NAME), NF90_DOUBLE, y_dimid, y_varid) )
CALL check( nf90_def_var(ncid, Z_NAME, NF90_DOUBLE, z_dimid, z_varid) )

!Assign units attributes to coordinate variables.
CALL check( nf90_put_att(ncid, x_varid,"long_name",TRIM(X_LONGNAME)) )
CALL check( nf90_put_att(ncid, x_varid, UNITS, TRIM(X_UNITS)) )
CALL check( nf90_put_att(ncid, y_varid,"long_name",TRIM(Y_LONGNAME)) )
CALL check( nf90_put_att(ncid, y_varid, UNITS, TRIM(Y_UNITS)) )
CALL check( nf90_put_att(ncid, z_varid,"long_name",Z_LONGNAME) )
CALL check( nf90_put_att(ncid, z_varid, UNITS, Z_UNITS) )

!The dimids array is used to pass the dimids of the dimensions of
!the netCDF variables. Both of the netCDF variables we are creating
!share the same four dimensions. In Fortran, the unlimited
!dimension must come last on the list of dimids.
dimids = (/ x_dimid, y_dimid, z_dimid/)

!Define the netCDF variables for the pressure and temperature data.
  CALL check( nf90_def_var(ncid, POTTEMP_NAME,    NF90_REAL, dimids, pottemp_varid   ) )
  CALL check( nf90_def_var(ncid, PRE_NAME,    NF90_REAL, dimids, pre_varid   ) )
  CALL check( nf90_def_var(ncid, QV_NAME,    NF90_REAL, dimids, qvC_varid   ) )

!Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, pottemp_varid,    UNITS, POTTEMP_UNITS   ) )
    CALL check( nf90_put_att(ncid, pottemp_varid,"description", "Pot. Temperature"            ) )
    CALL check( nf90_put_att(ncid, pre_varid,    UNITS, PRE_UNITS   ) )
    CALL check( nf90_put_att(ncid, pre_varid,"description", "PreCssure"            ) )
    CALL check( nf90_put_att(ncid, qvC_varid,    UNITS, QV_UNITS   ) )
    CALL check( nf90_put_att(ncid, qvC_varid,"description", "Water vapor mixing ratio"            ) )

!End define mode.
CALL check( nf90_enddef(ncid) )

!Write the coordinate variable data. This will put the latitudes
!and longitudes of our data grid into the netCDF file.
 DO i=1,5
   CALL check( nf90_put_var(ncid, x_varid, xLoc) )
 END DO
 DO i=94,99
   CALL check( nf90_put_var(ncid, y_varid, yLoc) )
 END DO
 DO i=1,nzC
   CALL check( nf90_put_var(ncid, z_varid, zLoc) )
END DO

! These settings tell netcdf to write one timestep of data.
! count = (/ nLonC, nLatC, nzC/)
! start = (/ 1, 1, 1/)
!Write  data.
CALL check( nf90_put_var(ncid,pottemp_varid,PotTempC(1:7,94:99,1:nzC,1)))!,start=start,count=count) )
CALL check( nf90_put_var(ncid,pre_varid,PreC(1:7,94:99,1:nzC,1)))!,start=start,count=count) )
CALL check( nf90_put_var(ncid,qvC_varid,QvC(1:7,94:99,1:nzC,1)))!,start=start,count=count) )



! Close the file. This causes netCDF to flush all buffers and make
! sure your data are REALly written to disk.
CALL check( nf90_close(ncid) )

DEALLOCATE(xLoc)
DEALLOCATE(yLoc)
DEALLOCATE(zLoc)



!NetCDF COARSE Output
ALLOCATE(xLoc(1:nLonC))
ALLOCATE(yLoc(1:nLatC))
ALLOCATE(zLoc(1:nzC))

xLoc(1)=LonC0
yLoc(1)=LatC0

DO i=2,nLonC
  xLoc(i)=xLoc(i-1)+dLonC
END DO
 
DO i=2,nLatC
  yLoc(i)=yLoc(i-1)+dLatC
END DO
 
DO i=1,nzC
  zLoc(i)=HeightWRF_theta(34,34,i,1) !zLoc(i-1)+1
END DO
 
!Create .nc
CALL check( nf90_create('PotTempPreQvCoarse.nc', nf90_clobber, ncid) )

!Define the dimensions.
CALL check( nf90_def_dim(ncid, TRIM(X_NAME), nLonC, x_dimid) )
CALL check( nf90_def_dim(ncid, TRIM(Y_NAME), nLatC, y_dimid) )
CALL check( nf90_def_dim(ncid, Z_NAME, nzC, z_dimid) )

!Define the Coordinate variables.
CALL check( nf90_def_var(ncid, TRIM(X_NAME), NF90_DOUBLE, x_dimid, x_varid) )
CALL check( nf90_def_var(ncid, TRIM(Y_NAME), NF90_DOUBLE, y_dimid, y_varid) )
CALL check( nf90_def_var(ncid, Z_NAME, NF90_DOUBLE, z_dimid, z_varid) )

!Assign units attributes to coordinate variables.
CALL check( nf90_put_att(ncid, x_varid,"long_name",TRIM(X_LONGNAME)) )
CALL check( nf90_put_att(ncid, x_varid, UNITS, TRIM(X_UNITS)) )
CALL check( nf90_put_att(ncid, y_varid,"long_name",TRIM(Y_LONGNAME)) )
CALL check( nf90_put_att(ncid, y_varid, UNITS, TRIM(Y_UNITS)) )
CALL check( nf90_put_att(ncid, z_varid,"long_name",Z_LONGNAME) )
CALL check( nf90_put_att(ncid, z_varid, UNITS, Z_UNITS) )

!The dimids array is used to pass the dimids of the dimensions of
!the netCDF variables. Both of the netCDF variables we are creating
!share the same four dimensions. In Fortran, the unlimited
!dimension must come last on the list of dimids.
dimids = (/ x_dimid, y_dimid, z_dimid/)

!Define the netCDF variables for the pressure and temperature data.
  CALL check( nf90_def_var(ncid, POTTEMP_NAME,    NF90_REAL, dimids, pottemp_varid   ) )
  CALL check( nf90_def_var(ncid, PRE_NAME,    NF90_REAL, dimids, pre_varid   ) )
  CALL check( nf90_def_var(ncid, QV_NAME,    NF90_REAL, dimids, qvC_varid   ) )

!Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, pottemp_varid,    UNITS, POTTEMP_UNITS   ) )
    CALL check( nf90_put_att(ncid, pottemp_varid,"description", "Pot. Temperature"            ) )
    CALL check( nf90_put_att(ncid, pre_varid,    UNITS, PRE_UNITS   ) )
    CALL check( nf90_put_att(ncid, pre_varid,"description", "PreCssure"            ) )
    CALL check( nf90_put_att(ncid, qvC_varid,    UNITS, QV_UNITS   ) )
    CALL check( nf90_put_att(ncid, qvC_varid,"description", "Water vapor mixing ratio"            ) )

!End define mode.
CALL check( nf90_enddef(ncid) )

!Write the coordinate variable data. This will put the latitudes
!and longitudes of our data grid into the netCDF file.
 DO i=1,nLonC
   CALL check( nf90_put_var(ncid, x_varid, xLoc) )
 END DO
 DO i=1,nLatC
   CALL check( nf90_put_var(ncid, y_varid, yLoc) )
 END DO
 DO i=1,nzC
   CALL check( nf90_put_var(ncid, z_varid, zLoc) )
END DO

! These settings tell netcdf to write one timestep of data.
! count = (/ nLonC, nLatC, nzC/)
! start = (/ 1, 1, 1/)
!Write  data.
CALL check( nf90_put_var(ncid,pottemp_varid,PotTempC(1:nLonC,1:nLatC,1:nzC,1)))!,start=start,count=count) )
CALL check( nf90_put_var(ncid,pre_varid,PreC(1:nLonC,1:nLatC,1:nzC,1)))!,start=start,count=count) )
CALL check( nf90_put_var(ncid,qvC_varid,QvC(1:nLonC,1:nLatC,1:nzC,1)))!,start=start,count=count) )



! Close the file. This causes netCDF to flush all buffers and make
! sure your data are REALly written to disk.
CALL check( nf90_close(ncid) )

DEALLOCATE(xLoc)
DEALLOCATE(yLoc)
DEALLOCATE(zLoc)




CONTAINS
SUBROUTINE check(STATUS)
         INTEGER, INTENT ( in) :: STATUS
     
         IF(STATUS /= nf90_noerr) THEN
           PRINT *, TRIM(nf90_strerror(STATUS))
           STOP "Stopped"
         END IF 
END SUBROUTINE check 
END PROGRAM ReadWRFnc
