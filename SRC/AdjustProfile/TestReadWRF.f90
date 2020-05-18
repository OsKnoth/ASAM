PROGRAM ReadWRF
  USE Profile_Mod
  USE typesizes
  USE netcdf
  IMPLICIT NONE
  INTEGER :: nxC,nyC,nzC 
  INTEGER :: nxF,nyF,nzF,wz,wzoffset
  INTEGER :: posx,posy,posz
  INTEGER :: wxr,wxl,wyr,wyl,wzr,wzl

  REAL(4) :: dlat,dlon,wlon0,wlat0,dxF,dyF,sdz,dxC,dyC,sx0,sy0,domaincenterx,domaincentery,wdsx,wdsy,fx,fy  !w for wrf, s for system
  REAL(4) :: dxr,dxl,dyr,dyl,tar,dzr,dzl
  REAL(4) :: HeightWRF(26)

  REAL(4), ALLOCATABLE :: T2(:,:)
  REAL(4), ALLOCATABLE :: PSFC(:,:)
  REAL(4), ALLOCATABLE :: U10(:,:)
  REAL(4), ALLOCATABLE :: V10(:,:)
  REAL(4), ALLOCATABLE :: Q2(:,:)
  REAL(4), ALLOCATABLE :: TSK(:,:)
  REAL(4), ALLOCATABLE :: HSFC(:,:)
  REAL(4), ALLOCATABLE :: RRTOT(:,:)
  REAL(4), ALLOCATABLE :: HFX(:,:)
  REAL(4), ALLOCATABLE :: LH(:,:)

  REAL(4), ALLOCATABLE :: TK(:,:,:)
  REAL(4), ALLOCATABLE :: PRESSURE(:,:,:)
  REAL(4), ALLOCATABLE :: UMET(:,:,:)
  REAL(4), ALLOCATABLE :: VMET(:,:,:)
  REAL(4), ALLOCATABLE :: W(:,:,:)
  REAL(4), ALLOCATABLE :: QVWRF(:,:,:)
  REAL(4), ALLOCATABLE :: QC(:,:,:)

  REAL(4), ALLOCATABLE :: T2tar(:,:)
  REAL(4), ALLOCATABLE :: PSFCtar(:,:)

  REAL(4), ALLOCATABLE :: TKtar(:,:,:)
  REAL(4), ALLOCATABLE :: PRESSUREtar(:,:,:)
  REAL(4), ALLOCATABLE :: QVWRFtar(:,:,:)


  INTEGER :: SizeOfReal=4
  INTEGER :: i,j,k
  INTEGER :: ix,iy,iz

  INTEGER :: nRec

  CHARACTER(40) :: FileNameWRF='prof/wrfout_wrf4km_asam_04.dat'  

  REAL(4) :: dlonL=78847.0
  REAL(4) :: dlatL=111132.0

! NetCDF
  INTEGER :: ncid
  CHARACTER (LEN = *), PARAMETER :: TK_NAME = "TK" 
  CHARACTER (LEN = *), PARAMETER :: PRES_NAME = "PRESSURE" 
  CHARACTER (LEN = *), PARAMETER :: QVWRF_NAME = "QVWRF" 
  CHARACTER (LEN = *), PARAMETER :: TK_UNITS = "K"
  CHARACTER (LEN = *), PARAMETER :: PRES_UNITS = "hPa"
  CHARACTER (LEN = *), PARAMETER :: QVWRF_UNITS = "kg/kg"
  CHARACTER (LEN = *), PARAMETER :: X_NAME = "xc"
  CHARACTER (LEN = *), PARAMETER :: Y_NAME = "yc"
  CHARACTER (LEN = *), PARAMETER :: Z_NAME = "lev"
  CHARACTER (LEN = *), PARAMETER :: X_LONGNAME = "x-coordinate in Cartesian system"
  CHARACTER (LEN = *), PARAMETER :: Y_LONGNAME = "y-coordinate in Cartesian system"
  CHARACTER (LEN = *), PARAMETER :: Z_LONGNAME = "Height at midpoints"
  CHARACTER (LEN = *), PARAMETER :: X_UNITS = "m"
  CHARACTER (LEN = *), PARAMETER :: Y_UNITS = "m"
  CHARACTER (LEN = *), PARAMETER :: Z_UNITS = "m"
  CHARACTER (LEN = *), PARAMETER :: GRID_TYPE = "cart"
  INTEGER :: x_dimid, y_dimid, z_dimid, tk_varid, pres_varid, qv_varid
  INTEGER :: x_varid, y_varid, z_varid 
  REAL(4), ALLOCATABLE :: xLoc(:),yLoc(:),zLoc(:)
  INTEGER :: dimids(3),start(3),count(3)
  CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
! read 
  REAL(RealKind), ALLOCATABLE :: c(:,:)
  CHARACTER(120) :: FileName


  nxC=136
  nyC=91
  nzC=26

  ALLOCATE(T2(nxC,nyC))
  ALLOCATE(PSFC(nxC,nyC))
  ALLOCATE(U10(nxC,nyC))
  ALLOCATE(V10(nxC,nyC))
  ALLOCATE(Q2(nxC,nyC))
  ALLOCATE(TSK(nxC,nyC))
  ALLOCATE(HSFC(nxC,nyC))
  ALLOCATE(RRTOT(nxC,nyC))
  ALLOCATE(HFX(nxC,nyC))
  ALLOCATE(LH(nxC,nyC))

  ALLOCATE(TK(nxC,nyC,nzC))
  ALLOCATE(PRESSURE(nxC,nyC,nzC))
  ALLOCATE(UMET(nxC,nyC,nzC))
  ALLOCATE(VMET(nxC,nyC,nzC))
  ALLOCATE(W(nxC,nyC,nzC))
  ALLOCATE(QVWRF(nxC,nyC,nzC))
  ALLOCATE(QC(nxC,nyC,nzC))


  OPEN(FILE=FileNameWRF,UNIT=10,STATUS='UNKNOWN'&
      ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)
  nRec=1
  DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) T2(i,j)
      nRec=nRec+1
    END DO  
  END DO  

  DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) PSFC(i,j)
      nRec=nRec+1
    END DO  
  END DO  

  DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) U10(i,j)
      nRec=nRec+1
    END DO  
  END DO  

  DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) V10(i,j)
      nRec=nRec+1
    END DO  
  END DO  

  DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) Q2(i,j)
      nRec=nRec+1
    END DO  
  END DO  

DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) TSK(i,j)
      nRec=nRec+1
    END DO  
END DO    

DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) HSFC(i,j)
      nRec=nRec+1
    END DO  
END DO    

DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) RRTOT(i,j)
      nRec=nRec+1
    END DO  
END DO    

DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) HFX(i,j)
      nRec=nRec+1
    END DO  
END DO    

DO j=1,nyC
    DO i=1,nxC
      READ(10,REC=nRec) LH(i,j)
      nRec=nRec+1
    END DO  
END DO    

DO k=1,nzC
    DO j=1,nyC
      DO i=1,nxC
        READ(10,REC=nRec) TK(i,j,k)
        nRec=nRec+1
      END DO  
    END DO  
END DO

DO k=1,nzC
    DO j=1,nyC
      DO i=1,nxC
        READ(10,REC=nRec) PRESSURE(i,j,k)
        nRec=nRec+1
      END DO  
    END DO  
END DO

DO k=1,nzC
    DO j=1,nyC
      DO i=1,nxC
        READ(10,REC=nRec) UMET(i,j,k)
        nRec=nRec+1
      END DO  
    END DO  
END DO

DO k=1,nzC
    DO j=1,nyC
      DO i=1,nxC
        READ(10,REC=nRec) VMET(i,j,k)
        nRec=nRec+1
      END DO  
    END DO  
END DO

DO k=1,nzC
    DO j=1,nyC
      DO i=1,nxC
        READ(10,REC=nRec) W(i,j,k)
        nRec=nRec+1
      END DO  
    END DO  
END DO

DO k=1,nzC
    DO j=1,nyC
      DO i=1,nxC
        READ(10,REC=nRec) QVWRF(i,j,k)
        nRec=nRec+1
      END DO  
    END DO  
END DO

DO k=1,nzC
    DO j=1,nyC
      DO i=1,nxC
        READ(10,REC=nRec) QC(i,j,k)
        nRec=nRec+1
      END DO  
    END DO  
END DO



!constructing area
dlat=0.02227027
dlon=0.02227027

wlon0=11.4946
wlat0=49.9988

  domaincenterx=12.927716
  domaincentery=51.525345

! wrf cell size in m
  dxC=dlon*dLonL
  dyC=dlat*dLatL

  HeightWRF(1)=125 !m
  HeightWRF(2)=250
  HeightWRF(3)=375
  HeightWRF(4)=500
  HeightWRF(5)=750
  HeightWRF(6)=1000
  HeightWRF(7)=1250
  HeightWRF(8)=1500
  HeightWRF(9)=1750
  HeightWRF(10)=2000
  HeightWRF(11)=2500
  HeightWRF(12)=3000
  HeightWRF(13)=3500
  HeightWRF(14)=4000
  HeightWRF(15)=4500
  HeightWRF(16)=5000
  HeightWRF(17)=6000
  HeightWRF(18)=7000
  HeightWRF(19)=8000
  HeightWRF(20)=9000
  HeightWRF(21)=10000
  HeightWRF(22)=11000
  HeightWRF(23)=12000
  HeightWRF(24)=14000
  HeightWRF(25)=16000
  HeightWRF(26)=18000

  nxF=100
  nyF=100
  nzF=50

  dxF=500
  dyF=500
  sdz=200
  !koordinatenursprung des zielsystems in wrf cells
  sx0 = CEILING(((domaincenterx-wlon0)*dLonL-nxF*dxF/2.0d0)/dxC)
  sy0 = CEILING(((domaincentery-wlat0)*dLatL-nyF*dyF/2.0d0)/dyC)

  WRITE(*,*) 'sx0 ',sx0
  WRITE(*,*) 'sy0 ',sy0

!target domain size (cells)
  wdsx=nxF*dxF/dxC
  wdsy=nyF*dyF/dyC

  ALLOCATE(T2tar(nxF,nyF))
  ALLOCATE(PSFCtar(nxF,nyF))

  ALLOCATE(TKtar(nxF,nyF,nzC))
  ALLOCATE(PRESSUREtar(nxF,nyF,nzC))
  ALLOCATE(QVWRFtar(nxF,nyF,nzC))

!2D-Interpolation

  CALL Interpolation2D(T2,T2Tar)

  CALL Interpolation3Dxy(TK,TKTar)
  CALL Interpolation3Dxy(Pressure,PressureTar)
  CALL Interpolation3Dxy(QVWRF,QVWRFTar)

  DO posz=1,nzF
    WRITE(*,*) posz,'TKTar',TKtar(50,50,posz)
  END DO
  DO posz=1,nzC
    WRITE(*,*) posz,'TK',TK(50,50,posz)
  END DO


! NetCDF
! Program variables to hold the data we will write out. 

  WRITE(*,*) 'nxC',nxC,'nxF',nxF
  WRITE(*,*) 'WDX',dxC,'dxF',dxF

  ALLOCATE(xLoc(1:nxC))
  ALLOCATE(yLoc(1:nyC)) 
  ALLOCATE(zLoc(1:nzC))

  xLoc(1)=wlon0+dxC/2.0
  yLoc(1)=wlat0+dyC/2.0
  zLoc(1)=HeightWRF(1)/2.0

  DO i=2,nxC
    xLoc(i)=xLoc(i-1)+dxC
  END DO

  DO i=2,nyC
    yLoc(i)=yLoc(i-1)+dyC
  END DO

  DO i=2,nzC
    zLoc(i)=HeightWRF(i-1)+(HeightWRF(i)-HeightWRF(i-1))/2.0
  END DO

! Create .nc
  CALL check( nf90_create('WRF.nc', nf90_clobber, ncid) )

! Define the dimensions.
  CALL check( nf90_def_dim(ncid, TRIM(X_NAME), nxC, x_dimid) )
  CALL check( nf90_def_dim(ncid, TRIM(Y_NAME), nyC, y_dimid) )
  CALL check( nf90_def_dim(ncid, TRIM(Z_NAME), nzC, z_dimid) )

! Define the Coordinate variables.
  CALL check( nf90_def_var(ncid, TRIM(X_NAME), NF90_DOUBLE, x_dimid, x_varid) )
  CALL check( nf90_def_var(ncid, TRIM(Y_NAME), NF90_DOUBLE, y_dimid, y_varid) )
  CALL check( nf90_def_var(ncid, TRIM(Z_NAME), NF90_DOUBLE, z_dimid, z_varid) )

! Assign units attributes to coordinate variables.
  CALL check( nf90_put_att(ncid, x_varid,"long_name",TRIM(X_LONGNAME)) )
  CALL check( nf90_put_att(ncid, x_varid, UNITS, TRIM(X_UNITS)) )
  CALL check( nf90_put_att(ncid, y_varid,"long_name",TRIM(Y_LONGNAME)) )
  CALL check( nf90_put_att(ncid, y_varid, UNITS, TRIM(Y_UNITS)) )
  CALL check( nf90_put_att(ncid, z_varid,"long_name",TRIM(Z_LONGNAME)) )
  CALL check( nf90_put_att(ncid, z_varid, UNITS, TRIM(Z_UNITS)) )

! The dimids array is used to pass the dimids of the dimensions of
! the netCDF variables. Both of the netCDF variables we are creating
! share the same four dimensions. In Fortran, the unlimited
! dimension must come last on the list of dimids.
  dimids = (/ x_dimid, y_dimid, z_dimid/)

! Define the netCDF variables for the pressure and temperature data.
  CALL check( nf90_def_var(ncid, TK_NAME,    NF90_REAL, dimids, tk_varid   ) )
  CALL check( nf90_def_var(ncid, PRES_NAME,    NF90_REAL, dimids, pres_varid   ) )
  CALL check( nf90_def_var(ncid, QVWRF_NAME,    NF90_REAL, dimids, qv_varid   ) )

! Assign units attributes to the netCDF variables.
  CALL check( nf90_put_att(ncid, tk_varid,    UNITS, TK_UNITS   ) )
  CALL check( nf90_put_att(ncid, tk_varid,"description", "Temperature"            ) )

  CALL check( nf90_put_att(ncid, pres_varid,    UNITS, PRES_UNITS   ) )
  CALL check( nf90_put_att(ncid, pres_varid,"description", "Pressure"            ) )
    
  CALL check( nf90_put_att(ncid, qv_varid,    UNITS, QVWRF_UNITS   ) )
  CALL check( nf90_put_att(ncid, qv_varid,"description", "Qvapor"            ) )
    
! End define mode.
  CALL check( nf90_enddef(ncid) )


! Write the coordinate variable data. This will put the latitudes
! and longitudes of our data grid into the netCDF file.
  DO i=1,nxC
    CALL check( nf90_put_var(ncid, x_varid, xLoc) )
  END DO
  DO i=1,nyC
    CALL check( nf90_put_var(ncid, y_varid, yLoc) )
  END DO
  DO i=1,nzC
    CALL check( nf90_put_var(ncid, z_varid, zLoc) )
  END DO

! Write TK data.
  CALL check( nf90_put_var(ncid,tk_varid,TK(1:nxC,1:nyC,1:nzC)) )
! Write PRESSURE data.
  CALL check( nf90_put_var(ncid,pres_varid,PRESSURE(1:nxC,1:nyC,1:nzC)) )

! Write QVWRF data.
  CALL check( nf90_put_var(ncid,qv_varid,QVWRFtar(1:nxC,1:nyC,1:nzC)) )
! Close the file. This causes netCDF to flush all buffers and make
! sure your data are REALly written to disk.
  CALL check( nf90_close(ncid) )

  DEALLOCATE(xLoc)
  DEALLOCATE(yLoc)
  DEALLOCATE(zLoc)

  WRITE(*,*) 'Adjust '
! Profile
  WRITE(*,*) "use: ./run outputfilename"
  ProfileType='Dry'
  CALL get_command_argument(1,FileName)
  CALL ComputeParameter
  NumSounding=26-wzOffset
  CALL ReadInput(FileName)
  WRITE(*,*) 'ProfileType ',ProfileType

  WRITE(*,*) 'wzOffset ',wzOffset,SIZE(Temp),SIZE(Pre),SIZE(QV)
  DO ix=1,nxC
    DO iy=1,nyC
      DO i=1,NumSounding
        Height(i)=HeightWRF(i+wzOffset)-(HeightWRF(1+wzOffset)+HeightWRF(i+wzOffset))/2.0d0
        Temp(i)=TK(ix,iy,i+wzOffset)
        Pre(i)=PRESSURE(ix,iy,i+wzOffset)*1.d2
        QV(i)=QVWRFTar(ix,iy,i+wzOffset)
      END DO  
      CALL ComputeProfile(c,Height,Pre=Pre,Temp=Temp,QV=QV)
    END DO
  END DO
  WRITE(*,*) 'Ende Compute Profile'
CONTAINS

SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS

  IF(STATUS /= nf90_noerr) THEN
    PRINT *, trim(nf90_strerror(STATUS))
      STOP "Stopped"
   END IF
END SUBROUTINE check

SUBROUTINE Interpolation2D(cCoarse,cFine)
  REAL(4) :: cCoarse(:,:)
  REAL(4) :: cFine(:,:)

  INTEGER :: ixF,iyF
  INTEGER :: ixR,iyR
  INTEGER :: ixL,iyL
  REAL(RealKind) :: fx,fy
  REAL(RealKind) :: dxR,dyR
  REAL(RealKind) :: dxL,dyL
  REAL(RealKind) :: tar

! 2D-Interpolation
  DO iyF=1,nyF
    DO ixF=1,nxF
      fx=(sx0*dxC+ixF*dxF+dxF/2.0d0)/dxC
      fy=(sy0*dyC+iyF*dyF+dyF/2.0d0)/dyC
      ixR=CEILING(fx)
      ixL=FLOOR(fx)
      IF (ixR==ixL) THEN
        ixL=ixL-1
      END IF
      dxR=ixR-fx
      dxL=fx-ixL
      iyR=CEILING(fy)
      iyL=FLOOR(fy)
      IF (iyR==iyL) THEN
        iyL=iyL-1
      END IF
      dyR=iyR-fy
      dyL=fy-iyL
      tar=0.0d0
      tar=tar+(1.0d0-dxR)*(1.0d0-dyR)*cCoarse(ixR,iyR)
      tar=tar+(1.0d0-dxR)*(1.0d0-dyL)*cCoarse(ixR,iyL)
      tar=tar+(1.0d0-dxL)*(1.0d0-dyR)*cCoarse(ixL,iyR)
      tar=tar+(1.0d0-dxL)*(1.0d0-dyL)*cCoarse(ixL,iyL)
      cFine(ixF,iyF)=tar
    END DO
  END DO
END SUBROUTINE Interpolation2D

SUBROUTINE Interpolation3Dxy(cCoarse,cFine)
  REAL(4) :: cCoarse(:,:,:)
  REAL(4) :: cFine(:,:,:)

  INTEGER :: ixF,iyF
  INTEGER :: izC
  INTEGER :: ixR,iyR
  INTEGER :: ixL,iyL
  REAL(RealKind) :: fx,fy,fz
  REAL(RealKind) :: dxR,dyR
  REAL(RealKind) :: dxL,dyL
  REAL(RealKind) :: tar

! 3D-Interpolation
  wzoffset=6
  DO izC=wzoffset,nzC
    DO iyF=1,nyF
      DO ixF=1,nxF
        fx=(sx0*dxC+ixF*dxF+dxF/2.0d0)/dxC
        fy=(sy0*dyC+iyF*dyF+dyF/2.0d0)/dyC
        ixR=CEILING(fx)
        ixL=FLOOR(fx)
        IF (ixR==ixL) THEN
          ixL=ixL-1
        END IF
        dxR=ixR-fx
        dxL=fx-ixL
        iyR=CEILING(fy)
        iyL=FLOOR(fy)
        IF (iyR==iyL) THEN
          iyL=iyL-1
        END IF
        dyR=iyR-fy
        dyL=fy-iyL
        tar=0.0d0
        tar=tar+(1.0d0-dxR)*(1.0d0-dyR)*cCoarse(ixR,iyR,izC)
        tar=tar+(1.0d0-dxR)*(1.0d0-dyL)*cCoarse(ixR,iyL,izC)
        tar=tar+(1.0d0-dxL)*(1.0d0-dyR)*cCoarse(ixL,iyR,izC)
        tar=tar+(1.0d0-dxL)*(1.0d0-dyL)*cCoarse(ixL,iyL,izC)
        cFine(ixF,iyF,izC)=tar
      END DO
    END DO
  END DO
END SUBROUTINE Interpolation3Dxy

SUBROUTINE Interpolation3D(cCoarse,cFine)
  REAL(4) :: cCoarse(:,:,:)
  REAL(4) :: cFine(:,:,:)

  INTEGER :: ixF,iyF,izF
  INTEGER :: izC
  INTEGER :: ixR,iyR,izR
  INTEGER :: ixL,iyL,izL
  REAL(RealKind) :: fx,fy,fz
  REAL(RealKind) :: dxR,dyR,dzR
  REAL(RealKind) :: dxL,dyL,dzL
  REAL(RealKind) :: tar

! 3D-Interpolation
  DO izF=1,nzF
    wzoffset=6
    izC=wzoffset
    DO WHILE (HeightWRF(izC)<=izF*sdz+HeightWRF(wzoffset)+sdz/2.0d0)
      izC=izC+1
    END DO
    izC=izC-1
    DO iyF=1,nyF
      DO ixF=1,nxF
        fx=(sx0*dxC+ixF*dxF+dxF/2.0d0)/dxC
        fy=(sy0*dyC+iyF*dyF+dyF/2.0d0)/dyC
        ixR=CEILING(fx)
        ixL=FLOOR(fx)
        IF (ixR==ixL) THEN
          ixL=ixL-1
        END IF
        dxR=ixR-fx
        dxL=fx-ixL
        iyR=CEILING(fy)
        iyL=FLOOR(fy)
        IF (iyR==iyL) THEN
          iyL=iyL-1
        END IF
        dyR=iyR-fy
        dyL=fy-iyL
        dzR=(HeightWRF(izC+1)-(izF*sdz+sdz/2.0d0+HeightWRF(wzoffset)))/(HeightWRF(izC+1)-HeightWRF(izC))
        dzL=(izF*sdz+sdz/2.0d0+HeightWRF(wzoffset)-HeightWRF(izC))/(HeightWRF(izC+1)-HeightWRF(izC))
        izR=izC+1
        izL=izC
        tar=0.0d0
        tar=tar+(1.0d0-dxR)*(1.0d0-dyR)*(1.0d0-dzR)*cCoarse(ixR,iyR,izR)
        tar=tar+(1.0d0-dxR)*(1.0d0-dyR)*(1.0d0-dzL)*cCoarse(ixR,iyR,izL)
        tar=tar+(1.0d0-dxR)*(1.0d0-dyL)*(1.0d0-dzR)*cCoarse(ixR,iyL,izR)
        tar=tar+(1.0d0-dxR)*(1.0d0-dyL)*(1.0d0-dzL)*cCoarse(ixR,iyL,izL)
        tar=tar+(1.0d0-dxL)*(1.0d0-dyR)*(1.0d0-dzR)*cCoarse(ixL,iyR,izR)
        tar=tar+(1.0d0-dxL)*(1.0d0-dyR)*(1.0d0-dzL)*cCoarse(ixL,iyR,izL)
        tar=tar+(1.0d0-dxL)*(1.0d0-dyL)*(1.0d0-dzR)*cCoarse(ixL,iyL,izR)
        tar=tar+(1.0d0-dxL)*(1.0d0-dyL)*(1.0d0-dzL)*cCoarse(ixL,iyL,izL)
        cFine(ixF,iyF,izF)=tar
      END DO
    END DO
  END DO
END SUBROUTINE Interpolation3D

END PROGRAM ReadWRF
