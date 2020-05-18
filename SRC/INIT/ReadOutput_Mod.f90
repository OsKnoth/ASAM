MODULE ReadOutput_mod

  USE Kind_Mod
  USE DataType_Mod
  USE Physics_Mod
  USE Parallel_Mod
  USE Chemie_Mod
  USE Output_Mod
  USE Aerosol_Mod

  IMPLICIT NONE

CONTAINS

SUBROUTINE ReadOutput(VecC,VecT)

  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)
  REAL(RealKind) :: TimeT=0.0d0 

  REAL(RealKind) :: xW,yW,zW
  INTEGER :: ixW,iyW,izW,ibW
  INTEGER :: ix,iy,iz
  INTEGER :: iTime
  INTEGER :: i,j
  CHARACTER*40 :: FileName

  CALL ReadOutputInit(VecC,VecT)
  CALL OutputSet(VecC,VecT) 
!
! Ausgabe Verteilung
!
  WRITE(*,*) 'Eingabe Punkt f�r Verteilung'
  READ(*,*) xW,yw,zW
  WRITE(*,*) 'Eingabe ZeitPunkt'
  READ(*,*) iTime
  WRITE(*,*) 'Eingabe FileName'
  READ(*,*) FileName
  ibW=-1
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    ixW=-1
    iyW=-1
    izW=-1
    Sx:DO ix=ix0+1,ix1
      IF (xP(ix-1)<=xW.AND.xW<=xP(ix)) THEN
        ixW=ix
        EXIT Sx
      END IF
    END DO Sx
    Sy:DO iy=iy0+1,iy1
      IF (yP(iy-1)<=yW.AND.yW<=yP(iy)) THEN
        iyW=iy
        EXIT Sy
      END IF
    END DO Sy
    Sz:DO iz=iz0+1,iz1
      IF (yP(iz-1)<=zW.AND.zW<=zP(iz)) THEN
        izW=iz
        EXIT Sz
      END IF
    END DO Sz
    IF (ixW>0.AND.iyW>0.AND.izW>0) THEN
      ibW=ib
      EXIT
    END IF
  END DO
  IF (ixW>0.AND.iyW>0.AND.izW>0.AND.ibW>0) THEN
    DO i=1,iTime
      CALL ReadOutputMetC(TimeT,VecC,VecT)
    END DO
    OPEN(UNIT=OutputUnit,FILE=TRIM(FileName),STATUS='UNKNOWN')
    DO i=1,nFrac
      WRITE(OutputUnit,*) r(i)*1.d6 &
              ,(VecT(ibW)%Vec(j)%c(ixW,iyW,izW,i)/LOG10(pFac)/1.d6,j=1,nAqua)
    END DO
  END IF

END SUBROUTINE ReadOutput

SUBROUTINE ReadOutputInit(VecC,VecT)

  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)
    
!---  local variables
    
  INTEGER :: i,NumOut,nbOut
  INTEGER :: ic,icT
  TYPE(MPI_File) :: fh
  CHARACTER(2) :: nFracChar

  NumData=0
  VectorComponentsOut=0
  IF (VelOut) THEN
    NumData=NumData+3
  END IF
  IF (thOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (tkeOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (disOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoVOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoCOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoIOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoROut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (RhoOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (DiffOut) THEN
    NumData=NumData+1
    VectorComponentsOut=VectorComponentsOut+1
  END IF
  IF (ChemieOut) THEN
    DO icT=1,NumAeroOut
      NumData=NumData+nFrac
      VectorComponentsOut=VectorComponentsOut+1
    END DO
    DO icT=1,NumGasOut
      NumData=NumData+1
      VectorComponentsOut=VectorComponentsOut+1
    END DO
  END IF

  ALLOCATE(VectorOut(nb))
!  ALLOCATE(DiffRho(VectorComponentsOut))
!  ALLOCATE(NameScalars(VectorComponentsOut))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    ALLOCATE(VectorOut(ibLoc)%Vec(VectorComponentsOut))
    ic=0
    IF (thOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(thPos)%c
      OutputC(ic)%NameScalar='PotTemp '
    END IF
    IF (tkeOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tkePos)%c
      OutputC(ic)%NameScalar='Tke     '
    END IF
    IF (disOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(disPos)%c
      OutputC(ic)%NameScalar='Diss    '
    END IF
    IF (RhoVOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoVPos)%c
      OutputC(ic)%NameScalar='RhoV      '
    END IF
    IF (RhoCOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoCPos)%c
      OutputC(ic)%NameScalar='RhoC      '
    END IF
    IF (RhoROut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoRPos)%c
      OutputC(ic)%NameScalar='RhoR      '
    END IF
    IF (RhoIOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoIPos)%c
      OutputC(ic)%NameScalar='RhoI      '
    END IF
    IF (RhoSOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(RhoSPos)%c
      OutputC(ic)%NameScalar='RhoS      '
    END IF
    IF (nvOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nvPos)%c
      OutputC(ic)%NameScalar='nv      '
    END IF
    IF (ncOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(ncPos)%c
      OutputC(ic)%NameScalar='nc      '
    END IF
    IF (nrOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nrPos)%c
      OutputC(ic)%NameScalar='nr      '
    END IF
    IF (niOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(niPos)%c
      OutputC(ic)%NameScalar='ni      '
    END IF
    IF (nsOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(nsPos)%c
      OutputC(ic)%NameScalar='ns      '
    END IF
    IF (tracer1Out) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tracer1Pos)%c
      OutputC(ic)%NameScalar='Tracer1 '
    END IF
    IF (tracer2Out) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      VectorOut(ibLoc)%Vec(ic)%c=>VecC(ibLoc)%Vec(tracer2Pos)%c
      OutputC(ic)%NameScalar='Tracer2 '
    END IF
    IF (RhoOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      VectorOut(ibLoc)%Vec(ic)%c=>RhoCell(ibLoc)%c
      OutputC(ic)%NameScalar='Rho     '
    END IF
    IF (DiffOut) THEN
      ic=ic+1
      OutputC(ic)%DiffRho='Rho'
      VectorOut(ibLoc)%Vec(ic)%c=>DiffKoeff(ibLoc)%c
      OutputC(ic)%NameScalar='DiffKoef'
    END IF
    IF (ChemieOut) THEN
      DO icT=1,NumAeroOut
        ic=ic+1
        IF (SpeciesName(AeroOut(icT))=='aNUMBER') THEN
          OutputC(ic)%DiffRho='Unsc'
        ELSE IF (SpeciesName(AeroOut(icT))=='aRELAX') THEN
          OutputC(ic)%DiffRho='Unsc'
        ELSE IF (SpeciesName(AeroOut(icT))=='aH2O') THEN
          OutputC(ic)%DiffRho='Unsc'
        ELSE
          OutputC(ic)%DiffRho='aH2O'
        END IF
        VectorOut(ibLoc)%Vec(ic)%c=>VecT(ibLoc)%Vec(AeroOut(icT))%c 
        OutputC(ic)%NameScalar=SpeciesName(AeroOut(icT))
      END DO
      DO icT=1,NumGasOut
        ic=ic+1
        OutputC(ic)%DiffRho='Unsc'
        VectorOut(ibLoc)%Vec(ic)%c=>VecT(ibLoc)%Vec(GasOut(icT))%c 
      END DO
    END IF
    DO icT=1,LenOutSpecial
      ic=ic+1
      OutputC(ic)%DiffRho='Unsc'
      OutputC(ic)%NameScalar=NameOutSpecial(icT)
    END DO

  END DO

  IF (OutputType=='AVS') THEN    
    IF (MyId == 0) THEN
!     OPEN(UNIT=90,FILE=TRIM(OutputFilename)//'Kopf',STATUS='UNKNOWN')
      OPEN(UNIT=90,FILE='Kopf',STATUS='UNKNOWN')
      WRITE(90,*) 'Z-AXIS'
      WRITE(90,*) iz1Out-iz0Out+1
      WRITE(90,'(8(1X,1PE12.5))') (domain%zP(i),i=iz0Out,iz1Out)
      WRITE(90,*)
      WRITE(90,*) 'Y-AXIS'
      WRITE(90,*) iy1Out-iy0Out+1
      WRITE(90,'(8(1X,1PE12.5))') (domain%yP(i),i=iy0Out,iy1Out)
      WRITE(90,*)
      WRITE(90,*) 'X-AXIS'
      WRITE(90,*) ix1Out-ix0Out+1
      WRITE(90,'(8(1X,1PE12.5))') (domain%xP(i),i=ix0Out,ix1Out)
      WRITE(90,*)
      WRITE(90,*) 'TIME'
      WRITE(90,*) 100
      WRITE(90,*) 0.0e0
      WRITE(90,*)
      WRITE(90,*) 'GRID'
      WRITE(90,'(8I10)') 0,ix1Out-ix0Out  
      WRITE(90,'(8I10)') 0,iy1Out-iy0Out  
      WRITE(90,'(8I10)') 0,iz1Out-iz0Out  
      nbOut=0
      DO ib=1,nb
        CALL Set(Floor(ib))
        IF (MIN(MAX(igx0,ix0Out),ix1Out)<     &
            MAX(MIN(igx1,ix1Out),ix0Out).AND. &
            MIN(MAX(igy0,iy0Out),iy1Out)<     &
            MAX(MIN(igy1,iy1Out),iy0Out).AND. &
            MIN(MAX(igz0,iz0Out),iz1Out)<     &
            MAX(MIN(igz1,iz1Out),iz0Out)) THEN
          nbOut=nbOut+1
        END IF
      END DO
      WRITE(90,*) nbOut
      DO ib=1,nb
        CALL Set(Floor(ib))
        IF (MIN(MAX(igx0,ix0Out),ix1Out)<     &
            MAX(MIN(igx1,ix1Out),ix0Out).AND. &
            MIN(MAX(igy0,iy0Out),iy1Out)<     &
            MAX(MIN(igy1,iy1Out),iy0Out).AND. &
            MIN(MAX(igz0,iz0Out),iz1Out)<     &
            MAX(MIN(igz1,iz1Out),iz0Out)) THEN
          WRITE(90,'(10I10)') MIN(MAX(igx0,ix0Out),ix1Out)-ix0Out, &
                              MAX(MIN(igx1,ix1Out),ix0Out)-ix0Out, &
                              MIN(MAX(igy0,iy0Out),iy1Out)-iy0Out, &
                              MAX(MIN(igy1,iy1Out),iy0Out)-iy0Out, &
                              MIN(MAX(igz0,iz0Out),iz1Out)-iz0Out, &
                              MAX(MIN(igz1,iz1Out),iz0Out)-iz0Out, &
                              RefineX,RefineY,RefineZ,ib
        END IF
      END DO
      WRITE(90,*)
      WRITE(90,*) 'MET-DATA'
      WRITE(90,*) VectorComponentsOut+1
      IF (VelOut)   WRITE(90,*) 'WIND(3)'
      IF (thOut)    WRITE(90,*) 'th'
      IF (tkeOut)   WRITE(90,*) 'tke'
      IF (disOut)   WRITE(90,*) 'dis'
      IF (RhoVOut)    WRITE(90,*) 'RhoV'
      IF (RhoCOut)    WRITE(90,*) 'RhoC'
      IF (RhoIOut)    WRITE(90,*) 'RhoI'
      IF (RhoROut)    WRITE(90,*) 'RhoR'
      IF (RhoOut)  WRITE(90,*) 'Rho'
      IF (DiffOut)  WRITE(90,*) 'D'
      IF (ChemieOut) THEN
        DO icT=1,NumAeroOut
          WRITE(nFracChar,'(I2)') nFrac
          WRITE(90,*) TRIM(SpeciesName(AeroOut(icT)))//'('//nFracChar//')'
        END DO
        DO icT=1,NumGasOut
          WRITE(90,*) SpeciesName(GasOut(icT))
        END DO
      END IF
      WRITE(90,*)
      WRITE(90,*) 'BINARY'
      CLOSE(90)

    END IF
!
!-- Determination of WriteOffsetC
    Floor(1)%WriteOffsetC=1
    DO ib=1,nb-1
      CALL Set(Floor(ib))
      ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
      ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
      iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
      iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
      iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
      iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
        NumOut=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
        Floor(ib+1)%WriteOffsetC=Floor(ib)%WriteOffsetC+ &
                                NumData*NumOut
    END DO 
    CALL Set(Floor(nb))
    ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
    ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
    iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
    iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
    iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
    iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
    NumOut=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
    OffsetShift=Floor(nb)%WriteOffsetC+NumData*NumOut
    CALL MPI_Barrier(MPI_Comm_World,MPIerr)
    CALL MPI_FILE_OPEN(MPI_COMM_WORLD,'Kopf', &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL,fh,MPIErr)
    CALL MPI_FILE_GET_SIZE(fh,Offset0,MPIErr)
    CALL MPI_FILE_CLOSE(fh,MPIErr)
  END IF

!-----------------------------------------------------------
END SUBROUTINE ReadOutputInit
SUBROUTINE ReadOutputMetC(ActTime,VecC,VecT)
!
!==================================================
!----  Output of Velocity Values
!==================================================
!

  REAL(RealKind) :: ActTime
  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecT(:)


!---  local variables
  INTEGER :: ix, iy, iz, it
  INTEGER :: iScalar,NumOut
  
  REAL(4), ALLOCATABLE :: Work4(:)
  REAL(4) :: Time
  INTEGER :: iWork4
  INTEGER :: ncMax
  TYPE(MPI_File) :: fh
  TYPE(MPI_Status) :: Status
  INTEGER(KIND=MPI_OFFSET_KIND) :: Offset
  INTEGER :: nReals
  CHARACTER(10) :: iName

!--   Preparing Output
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,OutputFileName, &
                       MPI_MODE_RDONLY, &
                       MPI_INFO_NULL,fh,MPIErr)
!--   Eingabe Zeitpunkt
  IF (MyId==0) THEN
    Offset=Offset0
    nReals=1
    Time=ActTime
    CALL MPI_FILE_READ_AT(fh,Offset,Time,nReals, &
                    MPI_REAL,Status,MPIErr)
  END IF
  ncMax=0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    ncMax=MAX(ncMax,Floor(ib)%nc)
  END DO
!-    Eingabe Geschwindigkeitsfeld
  ALLOCATE(Work4(NumData*ncMax))
  Work4=0.0e0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    ix0=MIN(MAX(igx0,ix0Out),ix1Out)/2**(-RefineX)
    ix1=MAX(MIN(igx1,ix1Out),ix0Out)/2**(-RefineX)
    iy0=MIN(MAX(igy0,iy0Out),iy1Out)/2**(-RefineY)
    iy1=MAX(MIN(igy1,iy1Out),iy0Out)/2**(-RefineY)
    iz0=MIN(MAX(igz0,iz0Out),iz1Out)/2**(-RefineZ)
    iz1=MAX(MIN(igz1,iz1Out),iz0Out)/2**(-RefineZ)
    Offset=Floor(ib)%WriteOffsetC*4+Offset0
    NumOut=(ix1-ix0)*(iy1-iy0)*(iz1-iz0)
    nReals=NumData*NumOut
    CALL MPI_FILE_READ_AT(fh,Offset,work4,nReals, &
                             MPI_REAL,Status,MPIErr)
    iWork4=1
    IF (VelOut) THEN
      DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          DO iz=iz0+1,iz1
            VecC(ibLoc)%Vec(uPosL)%c(ix,iy,iz,1)=Work4(iWork4)
            iWork4=iWork4+1 
          END DO
        END DO
      END DO
      DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          DO iz=iz0+1,iz1
            VecC(ibLoc)%Vec(vPosL)%c(ix,iy,iz,1)=Work4(iWork4)
            iWork4=iWork4+1 
          END DO
        END DO
      END DO
      DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          DO iz=iz0+1,iz1
            VecC(ibLoc)%Vec(wPosL)%c(ix,iy,iz,1)=Work4(iWork4)
            iWork4=iWork4+1
          END DO
        END DO
      END DO
    END IF
    DO iScalar=1,SIZE(VectorOut(ibLoc)%Vec)
      DO it=LBOUND(VectorOut(ibLoc)%Vec(iScalar)%c,4),&
            UBOUND(VectorOut(ibLoc)%Vec(iScalar)%c,4)
        DO ix=ix0+1,ix1
          DO iy=iy0+1,iy1
            DO iz=iz0+1,iz1
              VectorOut(ibLoc)%Vec(iScalar)%c(ix,iy,iz,it)=Work4(iWork4)
              iWork4=iWork4+1
            END DO
          END DO
        END DO
      END DO
    END DO
  END DO
  CALL MPI_FILE_CLOSE(fh,MPIErr)
  DEALLOCATE(Work4)
  Offset0=Offset0+4*OffsetShift

!-----------------------------------------------------------
END SUBROUTINE ReadOutputMetC

END MODULE ReadOutput_mod
