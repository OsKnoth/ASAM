MODULE Emission_Mod

  USE Chemie_Mod
  USE Distribution_Mod
  USE InitAerosol_Mod
  USE Physics_Mod
! New
  USE EmissDeposParameter_Mod
  USE MarineEmission_Mod
  USE FireEmission_Mod

  IMPLICIT NONE

CONTAINS

SUBROUTINE SetPointEmission(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: i,j
  INTEGER :: ix,iy,iz
  INTEGER :: NumPoints
  CHARACTER(30) :: SpeciesName,Type
  CHARACTER*300 :: Line
  REAL(RealKind), ALLOCATABLE :: Par(:),c(:)
  REAL(RealKind) :: x,y,z,Dummy
  INTEGER :: ibLocP,ixP,iyP,izP
  

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'TYPES')>0) THEN
      IF (MyID==0) WRITE (*,*) '***** Point emission(s) active *****'
      READ(InputUnit,*) NumberOfPointSnap
      ALLOCATE(PointSnap(NumberOfPointSnap))
      DO i=1,NumberOfPointSnap
        READ(InputUnit,*) PointSnap(i)%Type
      END DO
      EXIT
    END IF
  END DO
1 CONTINUE
  IF (NumberOfPointSnap==0) THEN ! Hinneburg
    IF (MyID==0) WRITE(*,*) '***** Attention: No point emissions *****'
    GOTO 5
  END IF
  REWIND(InputUnit)
  DO
    READ(InputUnit,*,END=2) Line
    IF (INDEX(Line,'BEGIN_SPAR')>0) THEN
      READ(InputUnit,*) Type
      IF (Type=='FOR_ALL') THEN
        READ(InputUnit,*) PointSnap(1)%NumPar
        ALLOCATE(PointSnap(1)%TypePar(PointSnap(1)%NumPar))
        DO i=1,PointSnap(1)%NumPar
          READ(InputUnit,*) PointSnap(1)%TypePar(i)
        END DO
        DO i=2,NumberOfPointSnap
          PointSnap(i)%NumPar=PointSnap(1)%NumPar
          ALLOCATE(PointSnap(i)%TypePar(PointSnap(i)%NumPar))
          PointSnap(i)%TypePar=PointSnap(1)%TypePar
        END DO
      ELSE
        DO i=1,NumberOfPointSnap
          READ(InputUnit,*) PointSnap(i)%NumPar 
          ALLOCATE(PointSnap(i)%TypePar(PointSnap(i)%NumPar))
          DO j=1,PointSnap(i)%NumPar
            READ(InputUnit,*) PointSnap(i)%TypePar(j)
          END DO
        END DO
      END IF
      EXIT
    END IF
  END DO
2 CONTINUE
  REWIND(InputUnit)
  DO
    READ(InputUnit,*,END=3) Line
    IF (INDEX(Line,'BEGIN_TRACER')>0) THEN
      READ(InputUnit,*) Type
      IF (Type=='FOR_ALL') THEN
        READ(InputUnit,*) PointSnap(1)%NumGas
        ALLOCATE(PointSnap(1)%Gas(PointSnap(1)%NumGas))
        DO i=1,PointSnap(1)%NumGas
          READ(InputUnit,*) SpeciesName
          PointSnap(1)%Gas(i)=Position(SpeciesName)
        END DO
        DO i=2,NumberOfPointSnap
          PointSnap(i)%NumGas=PointSnap(1)%NumGas
          ALLOCATE(PointSnap(i)%Gas(PointSnap(i)%NumGas))
          PointSnap(i)%Gas=PointSnap(1)%Gas
        END DO
      ELSE
        DO i=1,NumberOfPointSnap
          READ(InputUnit,*) PointSnap(i)%NumGas 
          ALLOCATE(PointSnap(i)%Gas(PointSnap(i)%NumGas))
          DO j=1,PointSnap(i)%NumGas
            READ(InputUnit,*) SpeciesName
            PointSnap(i)%Gas(j)=Position(SpeciesName)
          END DO
        END DO
      END IF
      EXIT
    END IF
  END DO
3 CONTINUE
  REWIND(InputUnit)
  DO
    READ(InputUnit,*,END=4) Line
    IF (INDEX(Line,'BEGIN_DATA')>0) THEN
      DO i=1,NumberOfPointSnap
        ALLOCATE(PointSnap(i)%EmiPointBlock(nbLoc))
        READ(InputUnit,*) Type
        READ(InputUnit,*) NumPoints
        DO j=1,NumPoints
          READ(InputUnit,*) x,y,z
          READ(InputUnit,*) Dummy
          CALL LocatePoint(x,y,z,ibLocP,ixP,iyP,izP)
          IF (ibLocP>0) THEN
            PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint= &
            PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint+1
          END IF
        END DO
      END DO
      EXIT
    END IF
  END DO
4 CONTINUE
  REWIND(InputUnit)
  DO
    READ(InputUnit,*,END=5) Line
    IF (INDEX(Line,'BEGIN_DATA')>0) THEN
      DO i=1,NumberOfPointSnap
        ALLOCATE(Par(1:PointSnap(i)%NumPar))
        ALLOCATE(c(1:PointSnap(i)%NumGas))
        DO ibLocP=1,nbLoc
          ALLOCATE(PointSnap(i)%EmiPointBlock(ibLocP) &
           %EmiPoint(PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint))
          PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint=0
        END DO
        READ(InputUnit,*) Type
        READ(InputUnit,*) NumPoints
        DO j=1,NumPoints
          READ(InputUnit,*) x,y,z,Par
          READ(InputUnit,*) c
          CALL LocatePoint(x,y,z,ibLocP,ixP,iyP,izP)
          IF (ibLocP>0) THEN
            PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint= &
            PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint+1
            ALLOCATE(PointSnap(i)%EmiPointBlock(ibLocP)% &
             EmiPoint(PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint)%Par(PointSnap(i)%NumPar))
            PointSnap(i)%EmiPointBlock(ibLocP)% &
             EmiPoint(PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint)%Par=Par
            ALLOCATE(PointSnap(i)%EmiPointBlock(ibLocP)% &
             EmiPoint(PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint)%c(PointSnap(i)%NumGas))
            PointSnap(i)%EmiPointBlock(ibLocP)% &
             EmiPoint(PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint)%c=c
            PointSnap(i)%EmiPointBlock(ibLocP)% &
             EmiPoint(PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint)%ix=ixP
            PointSnap(i)%EmiPointBlock(ibLocP)% &
             EmiPoint(PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint)%iy=iyP
            PointSnap(i)%EmiPointBlock(ibLocP)% &
             EmiPoint(PointSnap(i)%EmiPointBlock(ibLocP)%NumEmiPoint)%iz=izP
      END IF 
    END DO
    DEALLOCATE(Par)
    DEALLOCATE(c)
  END DO
  EXIT
END IF
  END DO
5 CONTINUE
  CLOSE(UNIT=InputUnit)

END SUBROUTINE SetPointEmission

SUBROUTINE SetEmission(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: i,j,iCell,Pos
  INTEGER :: ix,iy,iz
  INTEGER :: izStart
  INTEGER :: iSpecies,NumberOfSpecies
  INTEGER :: iMode,NumMode
  INTEGER :: iGas,NumGas,NumGasLoc,iAero,NumAero
  INTEGER :: DummyInt
  REAL(RealKind) :: Dummy
  CHARACTER*100 :: DummyChar
  REAL(RealKind) :: c1,c2,c3,c4,c5,c6,c7,c8
  REAL(RealKind) :: xQ0,xQ1,yQ0,yQ1,zQ0,zQ1
  REAL(RealKind) :: n1,n2
  REAL(RealKind) :: Temp=288.15d0-273.15d0
  REAL(RealKind) :: Num,D,Sigma 
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc,FacLoc 
  CHARACTER(20) :: S1,S2,End,SpeciesName
  CHARACTER*300 :: Line
  LOGICAL :: Back
  LOGICAL :: Intersect
  TYPE(AeroDistribution_T) :: DummyAeroDistribution


  NumGas=0
  NumAero=0
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'EmiDomain')>0) THEN
      IF (MyID==0) WRITE (*,*) '***** Area emission(s) active *****'
      READ(InputUnit,*) NumberOfEmiDomain
      NumberOfEmiDomainLoc=0
      DO i=1,NumberOfEmiDomain
        READ(InputUnit,*) xQ0,xQ1,yQ0,yQ1,zQ0,zQ1  
        DO ibLoc=1,nbLoc
          ib=LocGlob(ibLoc)
          CALL Set(Floor(ib))
          IF (.NOT.(x1<xQ0.OR.xQ1<x0.OR. &
              y1<yQ0.OR.yQ1<y0.OR. &
              z1<zQ0.OR.zQ1<z0)) THEN
            NumberOfEmiDomainLoc=NumberOfEmiDomainLoc+1  
            EXIT
          END IF  
        END DO  
        READ(InputUnit,*) n1,n2
        READ(InputUnit,*) NumGas                 
        IF (NumGas>0) THEN
          DO iGas=1,NumGas
            READ(InputUnit,*) DummyChar,Dummy
            IF (Position(TRIM(DummyChar))>0) THEN
              NumGas=NumGas+1
            END IF  
          END DO  
        END IF
        READ(InputUnit,*) NumAero                 
        IF (NumAero>0) THEN
          DO iAero=1,NumAero
            READ(InputUnit,*) DummyChar
          END DO  
          READ(InputUnit,*) NumMode                 
          DO iMode=1,NumMode
            READ(InputUnit,*) DummyAeroDistribution
          END DO  
          READ(InputUnit,*) DummyInt
          READ(InputUnit,*) Dummy
          READ(InputUnit,*) NumberOfSpecies
          DO iSpecies=1,NumberOfSpecies
            READ(InputUnit,*) DummyChar
          END DO  
        END IF
      END DO  
      EXIT
    END IF
  END DO
  1 CONTINUE
  IF (NumberOfEmiDomainLoc>0) THEN
    ALLOCATE(EmiDomain(NumberOfEmiDomainLoc))
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=2) Line
      IF (INDEX(Line,'EmiDomain')>0) THEN
        READ(InputUnit,*) NumberOfEmiDomain
        NumberOfEmiDomainLoc=0
        DO i=1,NumberOfEmiDomain
          READ(InputUnit,*) xQ0,xQ1,yQ0,yQ1,zQ0,zQ1  
          DO ibLoc=1,nbLoc
            ib=LocGlob(ibLoc)
            CALL Set(Floor(ib))
            Intersect=.FALSE.
            IF (.NOT.(x1<xQ0.OR.xQ1<x0.OR. &
                y1<yQ0.OR.yQ1<y0.OR. &
                z1<zQ0.OR.zQ1<z0)) THEN
              NumberOfEmiDomainLoc=NumberOfEmiDomainLoc+1  
              Intersect=.TRUE.
              EXIT
            END IF
          END DO  
          IF (Intersect) THEN
            EmiDomain(NumberOfEmiDomainLoc)%xQ0=xQ0
            EmiDomain(NumberOfEmiDomainLoc)%xQ1=xQ1
            EmiDomain(NumberOfEmiDomainLoc)%yQ0=yQ0
            EmiDomain(NumberOfEmiDomainLoc)%yQ1=yQ1
            EmiDomain(NumberOfEmiDomainLoc)%zQ0=zQ0
            EmiDomain(NumberOfEmiDomainLoc)%zQ1=zQ1
            READ(InputUnit,*) n1,n2
            EmiDomain(NumberOfEmiDomainLoc)%n1=n1
            EmiDomain(NumberOfEmiDomainLoc)%n2=n2
            READ(InputUnit,*) NumGasLoc                 
            IF (NumGasLoc>0) THEN
              IF (ChemieGas) THEN
                EmiDomain(NumberOfEmiDomainLoc)%NumGas=NumGasLoc
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%NameG(EmiDomain(NumberOfEmiDomainLoc)%NumGas))
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%EmiG(EmiDomain(NumberOfEmiDomainLoc)%NumGas))
                NumGas=0
                DO iGas=1,NumGasLoc
                  READ(InputUnit,*) DummyChar,Dummy
                  IF (Position(TRIM(DummyChar))>0) THEN
                    NumGas=NumGas+1
                    EmiDomain(NumberOfEmiDomainLoc)%NameG(NumGas)=DummyChar
                    EmiDomain(NumberOfEmiDomainLoc)%EmiG(NumGas)=Dummy
                  END IF  
                END DO  
                EmiDomain(NumberOfEmiDomainLoc)%NumGas=NumGas
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%PosG(EmiDomain(NumberOfEmiDomainLoc)%NumGas))
                DO iGas=1,EmiDomain(NumberOfEmiDomainLoc)%NumGas
                  EmiDomain(NumberOfEmiDomainLoc)%PosG(iGas)=Position(EmiDomain(NumberOfEmiDomainLoc)%NameG(iGas))
                END DO  
              ELSE  
                EmiDomain(NumberOfEmiDomainLoc)%NumGas=0
                DO iGas=1,NumGasLoc
                  READ(InputUnit,*) DummyChar,Dummy
                END DO  
              END IF
            ELSE  
              EmiDomain(NumberOfEmiDomainLoc)%NumGas=0
            END IF
            READ(InputUnit,*) NumAero                 
            IF (NumAero>0) THEN
              IF (Aerosol) THEN
                EmiDomain(NumberOfEmiDomainLoc)%NumAero=NumAero
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%NameA(EmiDomain(NumberOfEmiDomainLoc)%NumAero))
                DO iAero=1,EmiDomain(NumberOfEmiDomainLoc)%NumAero
                  READ(InputUnit,*) EmiDomain(NumberOfEmiDomainLoc)%NameA(iAero)
                END DO  
                READ(InputUnit,*) EmiDomain(NumberOfEmiDomainLoc)%NumberOfModes
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%AeroDistribution(EmiDomain(NumberOfEmiDomainLoc)%NumberOfModes))
                DO iMode=1,EmiDomain(NumberOfEmiDomainLoc)%NumberOfModes
                  READ(InputUnit,*) EmiDomain(NumberOfEmiDomainLoc)%AeroDistribution(iMode)
                END DO  
                READ(InputUnit,*) EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages%NumberOfStages
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages%Radius(EmiDomain(NumberOfEmiDomainLoc)%&
                & ImpactorStages%NumberOfStages))
                READ(InputUnit,*) EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages%Radius
                READ(InputUnit,*) EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages%NumberOfSpecies
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages% &
                         Fraction(EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages%NumberOfStages, &
                                  EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages%NumberOfSpecies))
                DO iSpecies=1,EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages%NumberOfSpecies
                  READ(InputUnit,*) DummyChar,EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages%Fraction(:,iSpecies)
                END DO
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%PosA(EmiDomain(NumberOfEmiDomainLoc)%NumAero))
                DO iAero=1,EmiDomain(NumberOfEmiDomainLoc)%NumAero
                  EmiDomain(NumberOfEmiDomainLoc)%PosA(iAero)=Position(EmiDomain(NumberOfEmiDomainLoc)%NameA(iAero))
                END DO  
                ALLOCATE(EmiDomain(NumberOfEmiDomainLoc)%EmiA(EmiDomain(NumberOfEmiDomainLoc)%NumAero+1,nFrac))
                CALL EmiAero(EmiDomain(NumberOfEmiDomainLoc)%AeroDistribution,EmiDomain(NumberOfEmiDomainLoc)%ImpactorStages &
                             ,EmiDomain(NumberOfEmiDomainLoc)%EmiA)
              ELSE    
                EmiDomain(NumberOfEmiDomainLoc)%NumAero=0
                DO iAero=1,NumAero
                  READ(InputUnit,*) DummyChar
                END DO  
                READ(InputUnit,*) NumMode                 
                DO iMode=1,NumMode
                  READ(InputUnit,*) DummyAeroDistribution
                END DO  
                READ(InputUnit,*) DummyInt
                READ(InputUnit,*) Dummy
                READ(InputUnit,*) NumberOfSpecies
                DO iSpecies=1,NumberOfSpecies
                  READ(InputUnit,*) DummyChar
                END DO  
              END IF
            ELSE  
              EmiDomain(NumberOfEmiDomainLoc)%NumAero=0
            END IF
          ELSE  
            READ(InputUnit,*) n1,n2
            READ(InputUnit,*) NumGas                 
            IF (NumGas>0) THEN
              DO iGas=1,NumGas
                READ(InputUnit,*) DummyChar,Dummy
              END DO  
            END IF
            READ(InputUnit,*) NumAero                 
            IF (NumAero>0) THEN
              DO iAero=1,NumAero
                READ(InputUnit,*) DummyChar
              END DO  
              READ(InputUnit,*) NumMode                 
              DO iMode=1,NumMode
                READ(InputUnit,*) DummyAeroDistribution
              END DO  
              READ(InputUnit,*) DummyInt
              READ(InputUnit,*) Dummy
              READ(InputUnit,*) NumberOfSpecies
              DO iSpecies=1,NumberOfSpecies
                READ(InputUnit,*) DummyChar
              END DO  
            END IF
          END IF  
        END DO  
      END IF
    END DO
    2 CONTINUE

    DO i=1,NumberOfEmiDomainLoc
      ALLOCATE(EmiDomain(i)%BlockCell(nbLoc))
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells=0
        DO iCell=1,SIZE(BoundCell)
          ix=BoundCell(iCell)%ix
          iy=BoundCell(iCell)%iy
          izStart=BoundCell(iCell)%iz
          DO iz=izStart,iz1
            dxLoc=MAX(MIN(xP(ix),EmiDomain(i)%xQ1)-MAX(xP(ix-1),EmiDomain(i)%xQ0),Zero)
            dyLoc=MAX(MIN(yP(iy),EmiDomain(i)%yQ1)-MAX(yP(iy-1),EmiDomain(i)%yQ0),Zero)
            dzLoc=MAX(MIN(zP(iz),EmiDomain(i)%zQ1)-MAX(zP(iz-1),EmiDomain(i)%zQ0),Zero)
            FacLoc=dxLoc*dyLoc*dzLoc/(dx(ix)*dy(iy)*dz(iz))
            IF (FacLoc>Zero.AND.VolC(ix,iy,iz)>0.0d0) THEN
              EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells=EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells+1
            END IF  
          END DO
        END DO
        ALLOCATE(EmiDomain(i)%BlockCell(ibLoc)%Cell(EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells))
      END DO
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells=0
        DO iCell=1,SIZE(BoundCell)
          ix=BoundCell(iCell)%ix
          iy=BoundCell(iCell)%iy
          izStart=BoundCell(iCell)%iz
          DO iz=izStart,iz1
            dxLoc=MAX(MIN(xP(ix),EmiDomain(i)%xQ1)-MAX(xP(ix-1),EmiDomain(i)%xQ0),Zero)
            dyLoc=MAX(MIN(yP(iy),EmiDomain(i)%yQ1)-MAX(yP(iy-1),EmiDomain(i)%yQ0),Zero)
            dzLoc=MAX(MIN(zP(iz),EmiDomain(i)%zQ1)-MAX(zP(iz-1),EmiDomain(i)%zQ0),Zero)
            FacLoc=dxLoc*dyLoc*dzLoc/(dx(ix)*dy(iy)*dz(iz))
            IF (FacLoc>Zero.AND.VolC(ix,iy,iz)>0.0d0) THEN
              EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells=EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells+1
              EmiDomain(i)%BlockCell(ibLoc)%Cell(EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells)%ix=ix
              EmiDomain(i)%BlockCell(ibLoc)%Cell(EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells)%iy=iy
              EmiDomain(i)%BlockCell(ibLoc)%Cell(EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells)%iz=iz
              EmiDomain(i)%BlockCell(ibLoc)%Cell(EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells)%Frac=FacLoc
              EmiDomain(i)%BlockCell(ibLoc)%Cell(EmiDomain(i)%BlockCell(ibLoc)%NumOfEmiCells)%FL=BoundCell(iCell)%FL
              BoundCell(iCell)%FacEmission=FacLoc
            END IF  
          END DO
        END DO
      END DO
    END DO
  END IF
  CLOSE(UNIT=InputUnit)

  S1='BEGIN_GAS'
  S2='BEGIN_EMISS'
  End='END_EMISS'
  CALL OpenFile(FileName)
  NumGasEmi=0
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      NumGasEmi=NumGasEmi+1
    END IF
  END DO
  CALL CloseFile
  ALLOCATE(GasEmi(NumGasEmi))
  NumGasEmi=0

  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,R1=c1)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      NumGasEmi=NumGasEmi+1
      GasEmi(NumGasEmi)%Pos=Position(SpeciesName)
      ALLOCATE(GasEmi(NumGasEmi)%Par(1))
      GasEmi(NumGasEmi)%Par(1)=c1
    END IF
  END DO
  CALL CloseFile

  S1='BEGIN_AERO'
  S2='BEGIN_EMISS'
  End='END_EMISS'
  CALL OpenFile(FileName)
  NumAerosolEmi=0
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      NumAerosolEmi=NumAerosolEmi+1
    END IF
  END DO
  CALL CloseFile
  ALLOCATE(AerosolEmi(NumAerosolEmi))
  NumAerosolEmi=0
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,R1=Num,R2=D,R3=Sigma)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      NumAerosolEmi=NumAerosolEmi+1
      AerosolEmi(NumAerosolEmi)%Pos=Position(SpeciesName)
      ALLOCATE(AerosolEmi(NumAerosolEmi)%AeroDistribution(1))
      AerosolEmi(NumAerosolEmi)%AeroDistribution(1)%Num=Num
      AerosolEmi(NumAerosolEmi)%AeroDistribution(1)%D=D
      AerosolEmi(NumAerosolEmi)%AeroDistribution(1)%Sigma=Sigma
      ALLOCATE(AerosolEmi(NumAerosolEmi)%Emi(2,nFrac))
      CALL EmiAero(AerosolEmi(NumAerosolEmi)%AeroDistribution,AerosolEmi(NumAerosolEmi)%ImpactorStages &
                  ,AerosolEmi(NumAerosolEmi)%Emi)
    END IF
  END DO
  CALL CloseFile

  IF (SeaEmiss) THEN
    CALL SetMarineEmission(FileName)
  END IF

END SUBROUTINE SetEmission

SUBROUTINE SetStreetEmission(FileName) ! Hinneburg

! Analogous to SetPointEmission

  CHARACTER(*) :: FileName
  CHARACTER(80) :: mi_wgts2
  INTEGER :: i,j
  INTEGER :: ix,iy,iz
  INTEGER :: NumPoints,SpeciesNum
  CHARACTER(5), ALLOCATABLE :: SpeciesName(:)
  REAL(RealKind), ALLOCATABLE :: c(:)


  mi_wgts2 = TRIM(FileName(1:INDEX(FileName,'.grid')-1))//'.Emission'
  OPEN(UNIT=InputUnit,FILE=mi_wgts2,STATUS='old')

  NumberOfStreetSnap=1
  ALLOCATE(StreetSnap(NumberOfStreetSnap))
  ALLOCATE(StreetSnap(1)%EmiStreetBlock(nbLoc))

  READ(InputUnit,*) SpeciesNum
  IF (SpeciesNum<=0) THEN
    WRITE(*,*) '***** Attention: SpeciesNum<=0 as read in from file',TRIM(mi_wgts2),' *****'
    NumberOfStreetSnap=0
    CLOSE(UNIT=InputUnit)
    RETURN
  END IF
  StreetSnap(1)%NumGas=SpeciesNum
  ALLOCATE(c(SpeciesNum))
  ALLOCATE(StreetSnap(1)%Gas(SpeciesNum))
  ALLOCATE(SpeciesName(SpeciesNum))

  WRITE(*,*) 'Check of emission file:  ',TRIM(mi_wgts2)
  j=0
  READ(InputUnit,*) SpeciesName
  DO i=1,SpeciesNum
    StreetSnap(1)%Gas(i)=Position(SpeciesName(i))
    IF      (StreetSnap(1)%Gas(i)/=0) THEN
      WRITE(*,*) '  ',SpeciesName(i),'  emission requested'
      j=j+1
    ELSE IF (StreetSnap(1)%Gas(i)==0) THEN
      WRITE(*,*) '  ',SpeciesName(i),'  emission not requested'
    END IF
  END DO
  IF      (j<VectorComponentsT-NumMet) THEN
    WRITE(*,*) '  -->    Some species remain without emissions'
  ELSE IF (j>VectorComponentsT-NumMet) THEN
    WRITE(*,*) '  -->    More emission entries than species'
  ELSE
    WRITE(*,*) '  -->    All species provided with emissions'
  END IF

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    READ(InputUnit,*) i,NumPoints
    IF (i/=ibLoc) WRITE(*,*) '*** Attention: Inconsistent block at SetStreetEmission input ***'
    StreetSnap(1)%EmiStreetBlock(i)%NumEmiStreet=NumPoints
    ALLOCATE(StreetSnap(1)%EmiStreetBlock(i)%EmiStreet(NumPoints))
    DO j=1,NumPoints
      ALLOCATE(StreetSnap(1)%EmiStreetBlock(i)%EmiStreet(j)%c(SpeciesNum))
      READ(InputUnit,*) ix,iy,iz
      READ(InputUnit,*) c ! µg/s
      StreetSnap(1)%EmiStreetBlock(i)%EmiStreet(j)%c=c
      StreetSnap(1)%EmiStreetBlock(i)%EmiStreet(j)%ix=ix
      StreetSnap(1)%EmiStreetBlock(i)%EmiStreet(j)%iy=iy
      StreetSnap(1)%EmiStreetBlock(i)%EmiStreet(j)%iz=iz
      IF (VolC(ix,iy,iz)<=Zero) THEN
        WRITE(*,1) '*** Inconsistency to GRID model: Emission in',i,ix,iy,iz,' meets closed cell (cancelled) ***'
1       FORMAT()
        StreetSnap(1)%EmiStreetBlock(i)%EmiStreet(j)%c=Zero
      END IF
    END DO
  END DO
  CLOSE(UNIT=InputUnit)
  DEALLOCATE(c,SpeciesName)

END SUBROUTINE SetStreetEmission

SUBROUTINE EmissionPointCompute(Time)

  REAL(RealKind) :: Time

  INTEGER :: i,j,k,LocNumGas0,LocNumGas1,t
  INTEGER :: Pos,NumGas,NumEmiPoint
  INTEGER :: ix,iy,iz
  INTEGER, POINTER :: Gas(:)
  TYPE(EmiPoint_T), POINTER :: EmiPoint(:)

  IF (Time<EmissionTimeStart.OR.NumberOfPointSnap==0) RETURN ! Hinneburg
  DO i=1,NumberOfPointSnap
    Gas=>PointSnap(i)%Gas   
    NumGas=PointSnap(i)%NumGas   
    NumEmiPoint=PointSnap(i)%EmiPointBlock(ibLoc)%NumEmiPoint
    EmiPoint=>PointSnap(i)%EmiPointBlock(ibLoc)%EmiPoint
    DO j=1,NumEmiPoint
      ix=EmiPoint(j)%ix
      iy=EmiPoint(j)%iy
      iz=EmiPoint(j)%iz
      DO k=1,NumGas
        Pos=Gas(k)  
        fVec(Pos)%c(ix,iy,iz,1)=fVec(Pos)%c(ix,iy,iz,1) &
               +EmiPoint(j)%c(k)/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE EmissionPointCompute

SUBROUTINE EmissionCompute(Time,VelF)

  REAL(RealKind) :: Time
  TYPE(VelocityFace_T) :: VelF(:)

  INTEGER :: i,iCell,ix,iy,iz
  INTEGER :: iGas,iAero
  INTEGER :: Pos,CountPos
  REAL(RealKind) :: n1,n2,n3,FL
  REAL(RealKind) :: FacEmission
  REAL(RealKind) :: EmissionWind
  REAL(RealKind) :: U10,TSurface
  REAL(RealKind) :: uLoc,vLoc,wLoc,VN,V,VT
  REAL(RealKind) :: uLoc2,vLoc2,wLoc2,VN2,V2,VT2
  REAL(RealKind) :: zPL
  REAL(RealKind) :: ustern, zRauh, DragM
  !REAL(RealKind), POINTER :: RhoLoc(:,:,:,:)
  REAL(RealKind) :: Frac

  !RhoLoc=>RhoC_T%c

  IF (Time<EmissionTimeStart.OR.NumberOfEmiDomain==0 &
       .OR.(NumberOfEmiDomain==0.AND.NumGasEmi==0.AND.NumAerosolEmi==0.AND..NOT.SeaEmiss)) RETURN ! Hinneburg

  DO iCell=1,NumBoundCell
    ix=BoundCell(iCell)%ix
    iy=BoundCell(iCell)%iy
    iz=BoundCell(iCell)%iz
    n1=BoundCell(iCell)%n1
    n2=BoundCell(iCell)%n2
    n3=BoundCell(iCell)%n3
    FL=BoundCell(iCell)%FL
    FacEmission=BoundCell(iCell)%FacEmission
    TSurface = BoundCell(iCell)%TeS
    U10 = BoundCell(iCell)%U10
    zRauh  = BoundCell(iCell)%zRauh
    DragM  = BoundCell(iCell)%DragM

    zPL = zP(iz-1)+0.5e0*dz(iz)

    IF (BoundCell(iCell)%LandClass==9.AND.SeaEmiss) THEN
      CALL SeaGasEmission(U10,FL,ix,iy,iz)
    ELSE  IF (BoundCell(iCell)%LandClass==4.AND.FireEmiss) THEN
      U10 = 0.5d0*(BoundCell(iCell)%VT1+BoundCell(iCell)%VT2)
      CALL FireEmission(U10,FL,ix,iy,iz)
    END IF

    IF (BoundCell(iCell)%LandClass==9.AND.SeaEmiss) THEN ! SeaSaltEmission
      ustern=SQRT(DragM)*U10*LOG(zPL/zRauh)/LOG(1.d1/zRauh)
      CALL SeaSaltEmission(U10,TSurface,ustern,zPL,zRauh,FL,ix,iy,iz)
    ELSE
      DO i=1,NumAerosolEmi
        Pos=AerosolEmi(i)%Pos
        fVec(iNc)%c(ix,iy,iz,:)=fVec(iNc)%c(ix,iy,iz,:) &
                               +AerosolEmi(i)%Emi(1,:)/3.6*FL/(VolC(ix,iy,iz)+Eps)*FacEmission
        fVec(Pos)%c(ix,iy,iz,:)=fVec(Pos)%c(ix,iy,iz,:) &
                               +AerosolEmi(i)%Emi(2,:)/3.6*FL/(VolC(ix,iy,iz)+Eps)*FacEmission
      END DO
    END IF
  END DO
  DO i=1,NumberOfEmiDomainLoc
    DO iCell=1,SIZE(EmiDomain(i)%BlockCell(ibLoc)%Cell)
      ix=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%ix
      iy=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%iy
      iz=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%iz
      Frac=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%Frac
      FL=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%FL
      DO iGas=1,EmiDomain(i)%NumGas
        Pos=EmiDomain(i)%PosG(iGas)
        fVec(Pos)%c(ix,iy,iz,:)=fVec(Pos)%c(ix,iy,iz,:) &
                               +EmiDomain(i)%EmiG(iGas)*FL/(VolC(ix,iy,iz)+Eps)*Frac
      END DO                         
      IF (EmiDomain(i)%NumAero>0) THEN
        fVec(iNC)%c(ix,iy,iz,:)=fVec(iNC)%c(ix,iy,iz,:) &
                               +EmiDomain(i)%EmiA(1,:)*FL/(VolC(ix,iy,iz)+Eps)*Frac
      END IF                         
      DO iAero=1,EmiDomain(i)%NumAero
        Pos=EmiDomain(i)%PosA(iAero)
        fVec(Pos)%c(ix,iy,iz,:)=fVec(Pos)%c(ix,iy,iz,:) &
                               +EmiDomain(i)%EmiA(iAero+1,:)*FL/(VolC(ix,iy,iz)+Eps)*Frac
      END DO                         
    END DO   
  END DO   
END SUBROUTINE EmissionCompute

SUBROUTINE EmissionStreetCompute(Time) ! Hinneburg

! Punktquellen (µg/s)

  REAL(RealKind) :: Time

  REAL(RealKind) :: FacEmission
  INTEGER :: i,j,k
  INTEGER :: Pos,NumGas,NumEmiStreets
  INTEGER :: ix,iy,iz
  INTEGER, POINTER :: Gas(:)
  TYPE(EmiStreet_T), POINTER :: EmiStreet(:)


  IF (Time<EmissionTimeStart.OR.NumberOfStreetSnap==0) RETURN

  DO i=1,NumberOfStreetSnap
    Gas=>StreetSnap(i)%Gas   
    NumGas=StreetSnap(i)%NumGas
    NumEmiStreets=StreetSnap(i)%EmiStreetBlock(ibLoc)%NumEmiStreet
    EmiStreet=>StreetSnap(i)%EmiStreetBlock(ibLoc)%EmiStreet
    DO j=1,NumEmiStreets
      ix=EmiStreet(j)%ix
      iy=EmiStreet(j)%iy
      iz=EmiStreet(j)%iz
      DO k=1,NumGas
        Pos=Gas(k)
        IF (Pos>0) THEN
          fVec(Pos)%c(ix,iy,iz,1)=fVec(Pos)%c(ix,iy,iz,1) &
                 +EmiStreet(j)%c(k)/(VolC(ix,iy,iz)+Eps)
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE EmissionStreetCompute

END MODULE Emission_Mod
