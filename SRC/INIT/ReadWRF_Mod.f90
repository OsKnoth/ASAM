MODULE ReadWRF_Mod
  USE Kind_Mod, ONLY : RealKind
  USE Control_Mod, ONLY : InputUnit,ComputeParameter,ThPos
  USE Domain_Mod, ONLY : dx,dy,nx,ny,nz,x0,y0,x1,y1,Set
  USE Floor_Mod, ONLY : Domain,ib,ibLoc
  USE Parallel_Mod, ONLY : nbLoc
  USE DataType_Mod, ONLY : Vector4Cell_T,FloorDomain=>Floor,LocGlob 
  USE Profile_Mod, ONLY :  &
!           Variables  
            NumSounding,ProfileType,ProfTemp,Pre,QV,Height, &
!           SUBROUTINES
            ComputeProfile,ReadInput
  IMPLICIT NONE

  CHARACTER(120) :: FileNameWRF='prof/wrfout_wrf4km_asam_04.dat'  
  REAL(8) :: x0F,y0F
  REAL(8) :: wlon0,wlat0
  NAMELIST /ReadWRFControl/ FileNameWRF &
                           ,wlon0 &
                           ,wlat0 

CONTAINS  
SUBROUTINE InputWRFData(FileName)
  CHARACTER(*) :: FileName

  CHARACTER(300) :: Line

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ReadWRFControl')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ReadWRFControl)
      EXIT
    END IF
  END DO  
1 CONTINUE
  CLOSE(UNIT=InputUnit)
END SUBROUTINE InputWRFData

SUBROUTINE ReadWRF(VecEnv)
  TYPE(Vector4Cell_T), POINTER :: VecEnv(:)

  INTEGER :: nxC,nyC,nzC 
  INTEGER :: nxF,nyF,nzF,wz,wzoffset

  REAL(8) :: lenXF,lenYF
  REAL(4), ALLOCATABLE :: dxF(:)
  REAL(4), ALLOCATABLE :: dyF(:)
  REAL(4), ALLOCATABLE :: xPF(:)
  REAL(4), ALLOCATABLE :: yPF(:)

  REAL(4) :: dlat,dlon,sdz,dxC,dyC,domaincenterx,domaincentery,wdsx,wdsy,fx,fy  !w for wrf, s for system
  REAL(4) :: HeightWRF(0:26)

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
  INTEGER :: i,j,k,kFail
  INTEGER :: ix,iy,iz

  INTEGER :: nRec


  REAL(4) :: dlonL=78847.0
  REAL(4) :: dlatL=111132.0

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


  WRITE(*,*) 'FileNameWRF   ',FileNameWRF
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

  DO i=1,nxC
    DO j=1,nyC
      kFail=0
      DO k=1,nzC
        IF (QVWRF(i,j,k)<=-900.e0) THEN
          kFail=kFail+1
        ELSE  
          EXIT
        END IF  
      END DO  
      DO k=1,kFail
        QVWRF(i,j,k)=QVWRF(i,j,kFail+1)
        TK(i,j,k)=TK(i,j,kFail+1)
        Pressure(i,j,k)=Pressure(i,j,kFail+1)
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

  HeightWRF(0)=0 !m
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
  lenXF=domain%x1-domain%x0
  lenYF=domain%y1-domain%y0
  x0F=((domaincenterx-wlon0)*dLonL-lenXF/2.0d0)/dxC
  y0F=((domaincentery-wlat0)*dLatL-lenYF/2.0d0)/dyC
  x0F=(domaincenterx-wlon0)*dLonL
  y0F=(domaincentery-wlat0)*dLatL

  CALL get_command_argument(1,FileName)
  CALL ComputeParameter
  NumSounding=26
  CALL ReadInput(FileName)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(FloorDomain(ib))
    WRITE(*,*) 'domain%x0 ',domain%x0
    WRITE(*,*) 'domain%y0 ',domain%y0
    WRITE(*,*) 'FloorDomain(ib)%x0 ',FloorDomain(ib)%x0
    WRITE(*,*) 'FloorDomain(ib)%y0 ',FloorDomain(ib)%y0
    nxF=nx
    nyF=ny
    nzF=nz
    ALLOCATE(dxF(nxF))
    ALLOCATE(dyF(nyF))
    dxF=dx
    dyF=dy
    ALLOCATE(xPF(0:nxF))
    ALLOCATE(yPF(0:nyF))
    wdsx=nxF*(x1-x0)/dxC
    wdsy=nyF*(y1-y0)/dyC
    xPF(0)=x0F+FloorDomain(ib)%x0-domain%x0
    DO ix=1,nxF
      xPF(ix)=xPF(ix-1)+dxF(ix)
    END DO  
    yPF(0)=y0F+FloorDomain(ib)%y0-domain%y0
    DO iy=1,nyF
      yPF(iy)=yPF(iy-1)+dyF(iy)
    END DO  
    WRITE(*,*) 'domaincenterx,wlon0',domaincenterx,wlon0
    WRITE(*,*) 'domaincentery,wlat0',domaincentery,wlat0
    WRITE(*,*) 'dLonL',dLonL,'dLatL',dLatL
    WRITE(*,*) 'dxC',dxC,'dyC',dyC
    WRITE(*,*) 'x0F',x0F
    WRITE(*,*) 'y0F',y0F
    ALLOCATE(T2tar(nxF,nyF))
    ALLOCATE(PSFCtar(nxF,nyF))
    ALLOCATE(TKtar(nxF,nyF,nzC))
    ALLOCATE(PRESSUREtar(nxF,nyF,nzC))
    ALLOCATE(QVWRFtar(nxF,nyF,nzC))
    CALL Interpolation2D(T2,T2Tar)
    CALL Interpolation3Dxy(TK,TKTar)
    CALL Interpolation3Dxy(Pressure,PressureTar)
    CALL Interpolation3Dxy(QVWRF,QVWRFTar)
    DO ix=1,nxF
      DO iy=1,nyF
        DO i=1,NumSounding
          Height(i)=HeightWRF(i)
          ProfTemp(i)=TKTar(ix,iy,i)
          Pre(i)=PRESSURETar(ix,iy,i)*1.d2
          QV(i)=QVWRFTar(ix,iy,i)
        END DO  
        WRITE(*,*) 'Block ib,ixF,iyF',ib,ix,iy
        CALL ComputeProfile(c,Height,Pre=Pre,Temp=ProfTemp,QV=QV)
      END DO  
    END DO  
    DEALLOCATE(dxF)
    DEALLOCATE(dyF)
    DEALLOCATE(xPF)
    DEALLOCATE(yPF)
    DEALLOCATE(T2tar)
    DEALLOCATE(PSFCtar)
    DEALLOCATE(TKtar)
    DEALLOCATE(PRESSUREtar)
    DEALLOCATE(QVWRFtar)
  END DO  
  STOP 'Ende Weiter Programmieren'

  nxF=100
  nyF=100
  nzF=50

  dxF=500
  dyF=500
  sdz=200
  !koordinatenursprung des zielsystems in wrf cells
! x0F = CEILING(((domaincenterx-wlon0)*dLonL-nxF*dxF/2.0d0)/dxC)
! y0F = CEILING(((domaincentery-wlat0)*dLatL-nyF*dyF/2.0d0)/dyC)

  WRITE(*,*) 'x0F ',x0F
  WRITE(*,*) 'y0F ',y0F

!target domain size (cells)
  wdsx=nxF*dxF(1)/dxC
  wdsy=nyF*dyF(1)/dyC

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


  WRITE(*,*) 'nxC',nxC,'nxF',nxF
  WRITE(*,*) 'WDX',dxC,'dxF',dxF(1)

  WRITE(*,*) 'Adjust '
! Profile
  WRITE(*,*) "use: ./run outputfilename"
  ProfileType='Dry'
  CALL get_command_argument(1,FileName)
  CALL ComputeParameter
  NumSounding=26
  CALL ReadInput(FileName)
  WRITE(*,*) 'ProfileType ',ProfileType

  WRITE(*,*) 'wzOffset ',wzOffset,SIZE(ProfTemp),SIZE(Pre),SIZE(QV)
  DO ix=1,nxC
    DO iy=1,nyC
      DO i=1,NumSounding
        Height(i)=HeightWRF(i)
        ProfTemp(i)=TK(ix,iy,i)
        Pre(i)=PRESSURE(ix,iy,i)*1.d2
        QV(i)=QVWRFTar(ix,iy,i)
      END DO  
      CALL ComputeProfile(c,Height,Pre=Pre,Temp=ProfTemp,QV=QV)
    END DO
  END DO
  WRITE(*,*) 'Ende Compute Profile'
  VecEnv(1)%Vec(ThPos)%c(1,1,1,1)=0.0d0
CONTAINS

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
      fx=(xPF(ixF-1)+dxF(ixF)/2.0d0)/dxC
      fy=(yPF(iyF-1)+dyF(iyF)/2.0d0)/dyC
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
  DO izC=1,nzC
    DO iyF=1,nyF
      DO ixF=1,nxF
        fx=(xPF(ixF-1)+dxF(ixF)/2.0d0)/dxC
        fy=(yPF(iyF-1)+dyF(iyF)/2.0d0)/dyC
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
    DO iyF=1,nyF
      DO ixF=1,nxF
        fx=(xPF(ixF-1)+dxF(ixF)/2.0d0)/dxC
        fy=(yPF(iyF-1)+dyF(iyF)/2.0d0)/dyC
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
END SUBROUTINE ReadWRF
END MODULE ReadWRF_Mod
