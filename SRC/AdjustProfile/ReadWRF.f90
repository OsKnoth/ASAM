SUBROUTINE ReadWRF
  USE Profile_Mod

  IMPLICIT NONE
  INTEGER :: nxC,nyC,nzC 
  INTEGER :: nxF,nyF,nzF,wz,wzoffset

  REAL(4) :: dlat,dlon,wlon0,wlat0,dxF,dyF,sdz,dxC,dyC,sx0,sy0,domaincenterx,domaincentery,wdsx,wdsy,fx,fy  !w for wrf, s for system
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


  WRITE(*,*) 'nxC',nxC,'nxF',nxF
  WRITE(*,*) 'WDX',dxC,'dxF',dxF

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

END SUBROUTINE ReadWRF
