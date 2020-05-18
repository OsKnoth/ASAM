PROGRAM ReadWRFTKtar 
  USE typesizes
  USE netcdf
IMPLICIT NONE
INTEGER :: xdef,ydef,zdef 
INTEGER :: ssx,ssy,ssz,wz,wzoffset
INTEGER :: posx,posy,posz
INTEGER :: wxr,wxl,wyr,wyl,wzr,wzl

REAL(4) :: dlat,dlon,wlon0,wlat0,sdx,sdy,sdz,wdx,wdy,sx0,sy0,domaincenterx,domaincentery,wdsx,wdsy,fx,fy  !w for wrf, s for system
REAL(4) :: dxr,dxl,dyr,dyl,tar,dzr,dzl
REAL(4) :: h(26)

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
REAL(4), ALLOCATABLE :: Q(:,:,:)
REAL(4), ALLOCATABLE :: QC(:,:,:)

REAL(4), ALLOCATABLE :: T2tar(:,:)
REAL(4), ALLOCATABLE :: PSFCtar(:,:)

REAL(4), ALLOCATABLE :: TKtar(:,:,:)


INTEGER :: SizeOfReal=4
INTEGER :: i,j,k

INTEGER :: nRec

CHARACTER(40) :: FileName='wrfout_wrf4km_asam_04_TK.dat'  

REAL(4) :: dlonL=78847.0
REAL(4) :: dlatL=111132.0
! NetCDF
INTEGER :: ncid
CHARACTER (LEN = *), PARAMETER :: TK_NAME = "TK" 
CHARACTER (LEN = *), PARAMETER :: TK_UNITS = "K"
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
INTEGER :: NXX, NYY, NZZ 
INTEGER :: x_dimid, y_dimid, z_dimid, tk_varid
INTEGER :: x_varid, y_varid, z_varid 
REAL(4), ALLOCATABLE :: xLoc(:),yLoc(:),zLoc(:)
INTEGER :: dimids(3),count(3),start(3)
CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"

!read 

xdef=11 !# of gridpoints in file to be read
ydef=11
zdef=26

ALLOCATE(TK(xdef,ydef,zdef))


OPEN(FILE=FileName,UNIT=10,STATUS='UNKNOWN'&
    ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)
nRec=1
DO k=1,zdef
  DO j=1,ydef
    DO i=1,xdef
      READ(10,REC=nRec) TK(i,j,k)
      nRec=nRec+1
    END DO  
  END DO  
END DO



!constructing area
dlat=0.02227027
dlon=0.02227027

wlon0=11.4946
wlat0=49.9988

domaincenterx=wlon0+ceiling(xdef/2.0)*dlon
domaincentery=wlat0+ceiling(ydef/2.0)*dlat

!wrf cell size in m
wdx=dlon*dLonL
wdy=dlat*dLatL





h(1)=125 !m
h(2)=250
h(3)=375
h(4)=500
h(5)=750
h(6)=1000
h(7)=1250
h(8)=1500
h(9)=1750
h(10)=2000
h(11)=2500
h(12)=3000
h(13)=3500
h(14)=4000
h(15)=4500
h(16)=5000
h(17)=6000
h(18)=7000
h(19)=8000
h(20)=9000
h(21)=10000
h(22)=11000
h(23)=12000
h(24)=14000
h(25)=16000
h(26)=18000

!Anz. Gitterpunkte feineres Gebiet
ssx=100
ssy=100
ssz=90

!Gitterweite im feineren Gebiet in m
sdx=8.0*dlon*dLonL/ssx     !Anz. Gitterpunkte des großen, in denen das feine liegt * Gitterweite des großen in ' * Gitterweite des großen in m / Anz. Gitterpunkte des feinen 
sdy=8.0*dlat*dLatL/ssy
sdz=200


!koordinatenursprung des zielsystems in wrf cells
sx0 = ceiling(((domaincenterx-wlon0)*dLonL-ssx*sdx/2.0d0)/wdx)
sy0 = ceiling(((domaincentery-wlat0)*dLatL-ssy*sdy/2.0d0)/wdy)

WRITE(*,*) 'sx0 ',sx0
WRITE(*,*) 'sy0 ',sy0

!target domain size (cells)
wdsx=ssx*sdx/wdx
wdsy=ssy*sdy/wdy

ALLOCATE(TKtar(ssx,ssy,ssz))



!3D-Interpolation
DO posz=1,ssz

  wzoffset=6
  wz=wzoffset
  DO WHILE (h(wz)<=posz*sdz+h(wzoffset)+sdz/2.0)
    wz=wz+1
  ENDDO
  wz=wz-1


  DO posy=1,ssy
   DO posx=1,ssx

    fx=(sx0*wdx+posx*sdx+sdx/2.0)/wdx
    fy=(sy0*wdy+posy*sdy+sdy/2.0)/wdy

    wxr=ceiling(fx)
    wxl=floor(fx)
    IF (wxr==wxl) THEN
      wxl=wxl-1
    END IF
    dxr=wxr-fx
    dxl=fx-wxl

    wyr=ceiling(fy)
    wyl=floor(fy)
    IF (wyr==wyl) THEN
      wyl=wyl-1
    END IF
    dyr=wyr-fy
    dyl=fy-wyl



    dzr=(h(wz+1)-(posz*sdz+sdz/2.0+h(wzoffset)))/(h(wz+1)-h(wz))
    dzl=(posz*sdz+sdz/2.0+h(wzoffset)-h(wz))/(h(wz+1)-h(wz))

    wzr=wz+1
    wzl=wz


    tar=0
    tar=tar+(1.0-dxr)*(1.0-dyr)*(1.0-dzr)*TK(wxr,wyr,wzr)
    tar=tar+(1.0-dxr)*(1.0-dyr)*(1.0-dzl)*TK(wxr,wyr,wzl)
    tar=tar+(1.0-dxr)*(1.0-dyl)*(1.0-dzr)*TK(wxr,wyl,wzr)
    tar=tar+(1.0-dxr)*(1.0-dyl)*(1.0-dzl)*TK(wxr,wyl,wzl)
    tar=tar+(1.0-dxl)*(1.0-dyr)*(1.0-dzr)*TK(wxl,wyr,wzr)
    tar=tar+(1.0-dxl)*(1.0-dyr)*(1.0-dzl)*TK(wxl,wyr,wzl)
    tar=tar+(1.0-dxl)*(1.0-dyl)*(1.0-dzr)*TK(wxl,wyl,wzr)
    tar=tar+(1.0-dxl)*(1.0-dyl)*(1.0-dzl)*TK(wxl,wyl,wzl)

    TKtar(posx,posy,posz)=tar
   END DO
  END DO
END DO



! NetCDF
! Program variables to hold the data we will write out. 
NXX = ssx 
NYY = ssy
NZZ = ssz

ALLOCATE(xLoc(1:NXX))
ALLOCATE(yLoc(1:NYY)) 
ALLOCATE(zLoc(1:NZZ))

xLoc(1)=sx0*wdx+sdx/2.0
yLoc(1)=sy0*wdy+sdy/2.0
zLoc(1)=sdz/2.0


DO i=2,ssx
  xLoc(i)=xLoc(i-1)+sdx
END DO

DO i=2,ssy
  yLoc(i)=yLoc(i-1)+sdy
END DO

DO i=2,ssz
  zLoc(i)=zLoc(i-1)+sdz
END DO

!Create .nc
CALL check( nf90_create('TKtar.nc', nf90_clobber, ncid) )

!Define the dimensions.
CALL check( nf90_def_dim(ncid, TRIM(X_NAME), NXX, x_dimid) )
CALL check( nf90_def_dim(ncid, TRIM(Y_NAME), NYY, y_dimid) )
CALL check( nf90_def_dim(ncid, Z_NAME, NZZ, z_dimid) )

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
  CALL check( nf90_def_var(ncid, TK_NAME,    NF90_REAL, dimids, tk_varid   ) )

!Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, tk_varid,    UNITS, TK_UNITS   ) )
    CALL check( nf90_put_att(ncid, tk_varid,"description", "Temperature interpolated"            ) )

!End define mode.
CALL check( nf90_enddef(ncid) )

!Write the coordinate variable data. This will put the latitudes
!and longitudes of our data grid into the netCDF file.
 DO i=1,NXX
   CALL check( nf90_put_var(ncid, x_varid, xLoc) )
 END DO
 DO i=1,NYY
   CALL check( nf90_put_var(ncid, y_varid, yLoc) )
 END DO
 DO i=1,NZZ
   CALL check( nf90_put_var(ncid, z_varid, zLoc) )
 END DO

! These settings tell netcdf to write one timestep of data. 
 count = (/ NXX, NYY, NZZ/)
 start = (/ 1, 1, 1/)
!Write TKtar data.
CALL check( nf90_put_var(ncid,tk_varid,TKtar(1:posx,1:posy,1:posz),start=start,count=count) )



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
    PRINT *, trim(nf90_strerror(STATUS))
      STOP "Stopped"
   END IF
END SUBROUTINE check

END PROGRAM ReadWRFTKtar
