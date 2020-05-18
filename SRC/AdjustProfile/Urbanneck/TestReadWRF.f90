PROGRAM ReadWRF
  USE Profile_Mod
  USE typesizes
  USE netcdf
IMPLICIT NONE
INTEGER :: xdef,ydef,zdef 
INTEGER :: ssx,ssy,ssz,wz,wzoffset
INTEGER :: posx,posy,posz
INTEGER :: wxr,wxl,wyr,wyl,wzr,wzl

REAL(4) :: dlat,dlon,wlon0,wlat0,sdx,sdy,sdz,wdx,wdy,sx0,sy0,domaincenterx,domaincentery,wdsx,wdsy,fx,fy  !w for wrf, s for system
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
REAL(4), ALLOCATABLE :: QV(:,:,:)
REAL(4), ALLOCATABLE :: QC(:,:,:)

REAL(4), ALLOCATABLE :: T2tar(:,:)
REAL(4), ALLOCATABLE :: PSFCtar(:,:)

REAL(4), ALLOCATABLE :: TKtar(:,:,:)
REAL(4), ALLOCATABLE :: PRESSUREtar(:,:,:)
REAL(4), ALLOCATABLE :: QVtar(:,:,:)


INTEGER :: SizeOfReal=4
INTEGER :: i,j,k

INTEGER :: nRec

CHARACTER(40) :: FileNameWRF='prof/wrfout_wrf4km_asam_04.dat'  

REAL(4) :: dlonL=78847.0
REAL(4) :: dlatL=111132.0

! NetCDF
INTEGER :: ncid
CHARACTER (LEN = *), PARAMETER :: TK_NAME = "TK" 
CHARACTER (LEN = *), PARAMETER :: PRES_NAME = "PRESSURE" 
CHARACTER (LEN = *), PARAMETER :: QV_NAME = "QV" 
CHARACTER (LEN = *), PARAMETER :: TK_UNITS = "K"
CHARACTER (LEN = *), PARAMETER :: PRES_UNITS = "hPa"
CHARACTER (LEN = *), PARAMETER :: QV_UNITS = "kg/kg"
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
INTEGER :: x_dimid, y_dimid, z_dimid, tk_varid, pres_varid, qv_varid
INTEGER :: x_varid, y_varid, z_varid 
REAL(4), ALLOCATABLE :: xLoc(:),yLoc(:),zLoc(:)
INTEGER :: dimids(3),start(3),count(3)
CHARACTER (LEN = *), PARAMETER :: UNITS      = "units"
!read 
INTEGER  INFO(15), IDID, LRW, IWORK(100000), LIW, IPAR(1000)
REAL(RealKind) :: &
      T, Y(90), YPRIME(90), TOUT, RWORK(100000), &
      RPAR(1000), RTOL, ATOL
REAL(RealKind), ALLOCATABLE :: c(:,:)
REAL(RealKind) :: zu,zo,Tu,To
REAL(RealKind) :: z,zOut
INTEGER :: iHeight,iz
REAL(RealKind) :: TeStart     = 288.0d0
REAL(RealKind) :: PreStart    = 1013.25d2
CHARACTER(40) :: FileName=''


xdef=136
ydef=91
zdef=26

ALLOCATE(T2(xdef,ydef))
ALLOCATE(PSFC(xdef,ydef))
ALLOCATE(U10(xdef,ydef))
ALLOCATE(V10(xdef,ydef))
ALLOCATE(Q2(xdef,ydef))
ALLOCATE(TSK(xdef,ydef))
ALLOCATE(HSFC(xdef,ydef))
ALLOCATE(RRTOT(xdef,ydef))
ALLOCATE(HFX(xdef,ydef))
ALLOCATE(LH(xdef,ydef))

ALLOCATE(TK(xdef,ydef,zdef))
ALLOCATE(PRESSURE(xdef,ydef,zdef))
ALLOCATE(UMET(xdef,ydef,zdef))
ALLOCATE(VMET(xdef,ydef,zdef))
ALLOCATE(W(xdef,ydef,zdef))
ALLOCATE(QV(xdef,ydef,zdef))
ALLOCATE(QC(xdef,ydef,zdef))


OPEN(FILE=FileNameWRF,UNIT=10,STATUS='UNKNOWN'&
    ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)
nRec=1
DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) T2(i,j)
!    WRITE(*,*) 'T2  ',i,j,T2(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) PSFC(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) U10(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) V10(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) Q2(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) TSK(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) HSFC(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) RRTOT(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) HFX(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO j=1,ydef
  DO i=1,xdef
    READ(10,REC=nRec) LH(i,j)
    nRec=nRec+1
  END DO  
END DO  

DO k=1,zdef
  DO j=1,ydef
    DO i=1,xdef
      READ(10,REC=nRec) TK(i,j,k)
      nRec=nRec+1
    END DO  
  END DO  
END DO

DO k=1,zdef
  DO j=1,ydef
    DO i=1,xdef
      READ(10,REC=nRec) PRESSURE(i,j,k)
      PRESSURE(i,j,k)=PRESSURE(i,j,k)*100.d0 
      nRec=nRec+1
    END DO  
  END DO  
END DO

DO k=1,zdef
  DO j=1,ydef
    DO i=1,xdef
      READ(10,REC=nRec) UMET(i,j,k)
      nRec=nRec+1
    END DO  
  END DO  
END DO

DO k=1,zdef
  DO j=1,ydef
    DO i=1,xdef
      READ(10,REC=nRec) VMET(i,j,k)
      nRec=nRec+1
    END DO  
  END DO  
END DO

DO k=1,zdef
  DO j=1,ydef
    DO i=1,xdef
      READ(10,REC=nRec) W(i,j,k)
      nRec=nRec+1
    END DO  
  END DO  
END DO

DO k=1,zdef
  DO j=1,ydef
    DO i=1,xdef
      READ(10,REC=nRec) QV(i,j,k)
      nRec=nRec+1
    END DO  
  END DO  
END DO

DO k=1,zdef
  DO j=1,ydef
    DO i=1,xdef
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

!wrf cell size in m
wdx=dlon*dLonL
wdy=dlat*dLatL





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

ssx=100
ssy=100
ssz=50

sdx=500
sdy=500
sdz=200
!koordinatenursprung des zielsystems in wrf cells
sx0 = CEILING(((domaincenterx-wlon0)*dLonL-ssx*sdx/2.0d0)/wdx)
sy0 = CEILING(((domaincentery-wlat0)*dLatL-ssy*sdy/2.0d0)/wdy)

WRITE(*,*) 'sx0 ',sx0
WRITE(*,*) 'sy0 ',sy0

!target domain size (cells)
wdsx=ssx*sdx/wdx
wdsy=ssy*sdy/wdy

ALLOCATE(T2tar(ssx,ssy))
ALLOCATE(PSFCtar(ssx,ssy))

ALLOCATE(TKtar(ssx,ssy,ssz))
ALLOCATE(PRESSUREtar(ssx,ssy,ssz))
ALLOCATE(QVtar(ssx,ssy,ssz))






!2D-Interpolation

DO posy=1,ssy
  DO posx=1,ssx


    fx=(sx0*wdx+posx*sdx+sdx/2.0)/wdx
    fy=(sy0*wdy+posy*sdy+sdy/2.0)/wdy

    dxr=ceiling(fx)-fx
    dxl=fx-floor(fx)

    dyr=ceiling(fy)-fy
    dyl=fy-floor(fy)

    wxr=ceiling(fx)
    wxl=floor(fx)

    wyr=ceiling(fy) 
    wyl=floor(fy)



    tar=0
    tar=tar+(1.0-dxr)*(1.0-dyr)*T2(wxr,wyr)
    tar=tar+(1.0-dxr)*(1.0-dyl)*T2(wxr,wyl)
    tar=tar+(1.0-dxl)*(1.0-dyl)*T2(wxl,wyl)
    tar=tar+(1.0-dxl)*(1.0-dyr)*T2(wxl,wyr)

    T2tar(posx,posy)=tar

  ENDDO
ENDDO


!3D-Interpolation
DO posz=1,ssz

  wzoffset=6
  wz=wzoffset
  DO WHILE (HeightWRF(wz)<=posz*sdz+HeightWRF(wzoffset)+sdz/2.0)
    wz=wz+1
  ENDDO
  wz=wz-1

!    dzr=(h(wz+1)-(posz*sdz+sdz/2.0+h(wzoffset)))/(h(wz+1)-h(wz))
!    dzl=(posz*sdz+sdz/2.0+h(wzoffset)-h(wz))/(h(wz+1)-h(wz))

    wzr=wz+1
    wzl=wz

  DO posy=1,ssy
   DO posx=1,ssx

    fx=(sx0*wdx+posx*sdx+sdx/2.0)/wdx
    fy=(sy0*wdy+posy*sdy+sdy/2.0)/wdy

    dxr=ceiling(fx)-fx
    dxl=fx-floor(fx)

    dyr=ceiling(fy)-fy
    dyl=fy-floor(fy)

    wxr=ceiling(fx)
    wxl=floor(fx)

    wyr=ceiling(fy)
    wyl=floor(fy)




    tar=0
    tar=tar+(1.0-dxr)*(1.0-dyr)*TK(wxr,wyr,wzr)
    tar=tar+(1.0-dxr)*(1.0-dyr)*TK(wxr,wyr,wzl)
    tar=tar+(1.0-dxr)*(1.0-dyl)*TK(wxr,wyl,wzr)
    tar=tar+(1.0-dxr)*(1.0-dyl)*TK(wxr,wyl,wzl)
    tar=tar+(1.0-dxl)*(1.0-dyr)*TK(wxl,wyr,wzr)
    tar=tar+(1.0-dxl)*(1.0-dyr)*TK(wxl,wyr,wzl)
    tar=tar+(1.0-dxl)*(1.0-dyl)*TK(wxl,wyl,wzr)
    tar=tar+(1.0-dxl)*(1.0-dyl)*TK(wxl,wyl,wzl)

    TKtar(posx,posy,posz)=tar

    tar=0
    tar=tar+(1.0-dxr)*(1.0-dyr)*PRESSURE(wxr,wyr,wzr)
    tar=tar+(1.0-dxr)*(1.0-dyr)*PRESSURE(wxr,wyr,wzl)
    tar=tar+(1.0-dxr)*(1.0-dyl)*PRESSURE(wxr,wyl,wzr)
    tar=tar+(1.0-dxr)*(1.0-dyl)*PRESSURE(wxr,wyl,wzl)
    tar=tar+(1.0-dxl)*(1.0-dyr)*PRESSURE(wxl,wyr,wzr)
    tar=tar+(1.0-dxl)*(1.0-dyr)*PRESSURE(wxl,wyr,wzl)
    tar=tar+(1.0-dxl)*(1.0-dyl)*PRESSURE(wxl,wyl,wzr)
    tar=tar+(1.0-dxl)*(1.0-dyl)*PRESSURE(wxl,wyl,wzl)
    

    PRESSUREtar(posx,posy,posz)=tar

    tar=0
    tar=tar+(1.0-dxr)*(1.0-dyr)*QV(wxr,wyr,wzr)
    tar=tar+(1.0-dxr)*(1.0-dyr)*QV(wxr,wyr,wzl)
    tar=tar+(1.0-dxr)*(1.0-dyl)*QV(wxr,wyl,wzr)
    tar=tar+(1.0-dxr)*(1.0-dyl)*QV(wxr,wyl,wzl)
    tar=tar+(1.0-dxl)*(1.0-dyr)*QV(wxl,wyr,wzr)
    tar=tar+(1.0-dxl)*(1.0-dyr)*QV(wxl,wyr,wzl)
    tar=tar+(1.0-dxl)*(1.0-dyl)*QV(wxl,wyl,wzr)
    tar=tar+(1.0-dxl)*(1.0-dyl)*QV(wxl,wyl,wzl)

    QVtar(posx,posy,posz)=tar

      ENDDO
    ENDDO
ENDDO


! NetCDF
! Program variables to hold the data we will write out. 
NXX = ssx 
NYY = ssy
NZZ = ssz

ALLOCATE(xLoc(1:NXX))
ALLOCATE(yLoc(1:NYY)) 
ALLOCATE(zLoc(1:NZZ))

xLoc(1)=wlon0+wdx/2.0
yLoc(1)=wlat0+wdy/2.0
zLoc(1)=HeightWRF(1)/2.0

DO i=2,xdef
  xLoc(i)=xLoc(i-1)+wdx
END DO

DO i=2,ydef
  yLoc(i)=yLoc(i-1)+wdy
END DO

DO i=2,zdef
  zLoc(i)=HeightWRF(i-1)+(HeightWRF(i)-HeightWRF(i-1))/2.0
END DO

!Create .nc
CALL check( nf90_create('WRF.nc', nf90_clobber, ncid) )

!Define the dimensions.
CALL check( nf90_def_dim(ncid, TRIM(X_NAME), NXX, x_dimid) )
CALL check( nf90_def_dim(ncid, TRIM(Y_NAME), NYY, y_dimid) )
CALL check( nf90_def_dim(ncid, TRIM(Z_NAME), NZZ, z_dimid) )

!Define the Coordinate variables.
CALL check( nf90_def_var(ncid, TRIM(X_NAME), NF90_DOUBLE, x_dimid, x_varid) )
CALL check( nf90_def_var(ncid, TRIM(Y_NAME), NF90_DOUBLE, y_dimid, y_varid) )
CALL check( nf90_def_var(ncid, TRIM(Z_NAME), NF90_DOUBLE, z_dimid, z_varid) )

!Assign units attributes to coordinate variables.
CALL check( nf90_put_att(ncid, x_varid,"long_name",TRIM(X_LONGNAME)) )
CALL check( nf90_put_att(ncid, x_varid, UNITS, TRIM(X_UNITS)) )
CALL check( nf90_put_att(ncid, y_varid,"long_name",TRIM(Y_LONGNAME)) )
CALL check( nf90_put_att(ncid, y_varid, UNITS, TRIM(Y_UNITS)) )
CALL check( nf90_put_att(ncid, z_varid,"long_name",TRIM(Z_LONGNAME)) )
CALL check( nf90_put_att(ncid, z_varid, UNITS, TRIM(Z_UNITS)) )

!The dimids array is used to pass the dimids of the dimensions of
!the netCDF variables. Both of the netCDF variables we are creating
!share the same four dimensions. In Fortran, the unlimited
!dimension must come last on the list of dimids.
dimids = (/ x_dimid, y_dimid, z_dimid/)

!Define the netCDF variables for the pressure and temperature data.
  CALL check( nf90_def_var(ncid, TK_NAME,    NF90_REAL, dimids, tk_varid   ) )
  CALL check( nf90_def_var(ncid, PRES_NAME,    NF90_REAL, dimids, pres_varid   ) )
  CALL check( nf90_def_var(ncid, QV_NAME,    NF90_REAL, dimids, qv_varid   ) )

!Assign units attributes to the netCDF variables.
    CALL check( nf90_put_att(ncid, tk_varid,    UNITS, TK_UNITS   ) )
    CALL check( nf90_put_att(ncid, tk_varid,"description", "Temperature"            ) )

    CALL check( nf90_put_att(ncid, pres_varid,    UNITS, PRES_UNITS   ) )
    CALL check( nf90_put_att(ncid, pres_varid,"description", "Pressure"            ) )
    
    CALL check( nf90_put_att(ncid, qv_varid,    UNITS, QV_UNITS   ) )
    CALL check( nf90_put_att(ncid, qv_varid,"description", "Qvapor"            ) )
    


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


!Write TK data.
CALL check( nf90_put_var(ncid,tk_varid,TKtar(1:ssx,1:ssy,1:ssz)) )

!Write PRESSURE data.
CALL check( nf90_put_var(ncid,pres_varid,PRESSUREtar(1:ssx,1:ssy,1:ssz)) )

!Write QV data.
CALL check( nf90_put_var(ncid,qv_varid,QVtar(1:ssx,1:ssy,1:ssz)) )




! Close the file. This causes netCDF to flush all buffers and make
! sure your data are REALly written to disk.
CALL check( nf90_close(ncid) )

DEALLOCATE(xLoc)
DEALLOCATE(yLoc)
DEALLOCATE(zLoc)



!Profile
  WRITE(*,*) "use: ./run outputfilename"
  CALL get_command_argument(1,FileName)
  CALL ComputeParameter
  CALL ReadInput(FileName)
  iPar(1)=CaseProfile

  t=0.0
  info(1:15)=0
  lrw=100000
  liw=100000
  y=0.0d0

  NumSounding=26-wzOffset
  ALLOCATE(Height(NumSounding))
  DO i=1,NumSounding
    Height(i)=HeightWRF(i+wzOffset)-(HeightWRF(1+wzOffset)+HeightWRF(i+wzOffset))/2.0d0
  END DO  
  IF (ipar(1)==2.OR.ipar(1)==3) THEN
    rpar(7)=1.d-4
    rpar(8)=rpar(7)
    CALL InputRes(y,yprime,rpar,ipar,PreStart=REAL(Pressure(1,1,wzOffset+1),8),TeStart=REAL(TK(1,1,wzOffset+1),8) &
                 ,qvstart=rpar(7))
    DO  
      rwork=0.0d0
      info=0.0d0
      INFO(11)=1 ! DDASSL computes consisten initial values
      t=0.0d0
      tout=1.0d0
      WRITE(*,*) 'Vor DDASSL',RPAR(5)
      CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
                  IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
      IF (rpar(7)==REAL(QV(1,1,wzOffset+1),8)) THEN
        EXIT
      END IF  
      rpar(7)=MIN(rpar(7)+0.0001d0,REAL(QV(1,1,wzOffset+1),8))
      rpar(8)=rpar(7)
      CALL InputRes(y,yprime,rpar,ipar,PreStart=REAL(Pressure(1,1,wzOffset+1),8),TeStart=REAL(TK(1,1,wzOffset+1),8) &
                   ,qvstart=rpar(7))
    END DO  
  ELSE  
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart)
    rwork=0.0d0
    info=0.0d0
    info(11)=1 ! DDASSL computes consisten initial values
    t=0.0d0
    tout=1.0d0
    CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
                IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
  END IF  
  STOP
  iHeight=1
  iz=1
  ! Forward
  WRITE(*,*) 'Forward integration'
  S1:DO 
    zu=Height(iHeight)
    rpar(1)=zu
    zo=Height(iHeight+1)
    rpar(2)=zo
    z=zu
    rpar(3)=Temp(iHeight)
    rpar(4)=Temp(iHeight+1)
    rpar(5)=RH(iHeight)
    rpar(6)=RH(iHeight+1)
    c(1,:)=y(1:neq)
    S2:DO 
      IF (iz>nz) THEN
        EXIT S1
      END IF  
      IF (zu<zM(iz).AND.zM(iz)<=zO) THEN
        zOut=zM(iz)
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        c(iz,:)=y(1:neq)
        IF (zM(iz)==zO) THEN
          iHeight=iHeight+1
          iz=iz+1
          EXIT S2
        END IF 
        iz=iz+1
      ELSE IF (zM(iz)<=zu) THEN
        iz=iz+1
      ELSE 
        zOut=zo
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        iHeight=iHeight+1
        EXIT S2
      END IF  
    END DO S2  
  END DO S1 

  ! Rückwärts
!  =======
  ! Backward
  WRITE(*,*) 'Backward integration'
  yprime=0.0d0
  y=0.0d0
  IF (ipar(1)==2) THEN
    rpar(5)=0.1
    rpar(6)=rpar(5)
    y(6)=rpar(5)
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=rpar(5))
    DO  
      rwork=0.0d0
      info=0.0d0
      INFO(11)=1 ! DDASSL computes consisten initial values
      t=0.0d0
      tout=1.0d0
      CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
                  IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
      IF (rpar(5)==RH(1)) THEN
        EXIT
      END IF  
      rpar(5)=MIN(rpar(5)+0.1d0,RH(1))
      rpar(6)=rpar(5)
      CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart,rhstart=rpar(5))
    END DO  
  ELSE  
    CALL InputRes(y,yprime,rpar,ipar,PreStart,TeStart=TeStart)
    rwork=0.0d0
    info=0.0d0
    info(11)=1 ! DDASSL computes consisten initial values
    t=0.0d0
    tout=1.0d0
    CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL, &
                IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
  END IF  
  rwork=0.0d0
  info=0
  S3:DO 
    zu=Height(1)
    rpar(1)=zu
    zo=Height(2)
    rpar(2)=zo
    z=zu
    rpar(3)=Temp(1)
    rpar(4)=Temp(2)
    rpar(5)=RH(1)
    rpar(6)=RH(2)
    iz=nz
    S4:DO 
      IF (iz<1) THEN
        EXIT S3
      END IF  
      IF (zu>zM(iz)) THEN
        zOut=zM(iz)
        CALL DDASSL (RES, NEQ, z, Y, YPRIME, zOUT, INFO, RTOL, ATOL, &
                     IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
        c(iz,:)=y(1:neq)
        iz=iz-1
      ELSE IF (zM(iz)>=zu) THEN
        iz=iz-1
      END IF  
    END DO S4  
  END DO S3 


CONTAINS
SUBROUTINE check(STATUS)
  INTEGER, INTENT ( in) :: STATUS

  IF(STATUS /= nf90_noerr) THEN
    PRINT *, trim(nf90_strerror(STATUS))
      STOP "Stopped"
   END IF
END SUBROUTINE check


END PROGRAM ReadWRF
