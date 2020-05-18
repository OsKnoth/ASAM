PROGRAM Main

  USE OSD_Mod
  USE DBase_Mod
  USE Shape_Mod
  USE Control_Mod

  IMPLICIT NONE

  CHARACTER(80) :: FileName=''
  CHARACTER(80) :: FileNameGrid=''
  INTEGER :: OutputUnitHaus=42
  CHARACTER(10) :: OutType,in_type
  INTEGER :: i,l,p, nr_wahlStruktur

  TYPE(PointOSD_T) :: Last_StoreP
  TYPE(LineOSD_T), POINTER :: Current
  TYPE(PolyLine_T) :: Line
  TYPE(PolyLineZ_T) :: PolyLineZ
  REAL(RealDouble) :: h_xmin,h_xmax,h_ymin,h_ymax,h_zmin,h_zmax
  TYPE(Point_T) :: Pzw 
  INTEGER :: nr_oj
  INTEGER :: LengthAllObjekte
  TYPE(BuildOSD_T) :: Build
  TYPE(BuildOSD_T), POINTER :: CurrBuild

  TYPE(Polygon_T), POINTER :: NewPoly
  TYPE(Polygon_T), POINTER :: CurrPoly
  TYPE(Polygon_T), POINTER :: StartPoly
  TYPE(Polygon_T), POINTER :: Temp
  !TYPE(RecordDbf_T) :: Record !dh
  INTEGER :: NumHaus,NumHaus1,NumberT,SpeciesNum
  LOGICAL :: RecordEnd,RecordErr
  CHARACTER(1)  :: q

! ================================== Erlaeuterungen ========================================
! Objekt-Arten:
!   Haus:   Mehreckige (anti-clockwise) Grundflaeche mit Parametern (Base,Height)
!           (im output 1 Ecke mehr: Startecke --> Endecke)
!   Street: Mehreckiger Linienzug mit Parametern (Base,Height,Width,Emissions),
!           Linienzug = Mittellinie der Strasse, Base = Strassenniveau, Width = Gesamtbreite,
!           Height = Hoehe der Fahrzeug- bzw. Emissionsschicht
!   Baum:   Linie mit Parametern (Base,Height,Width,Density[0...1])
!           Linie = Mittellinie der Grundflaeche eines Quaders mit Blattwerk, Base = Hoehe
!           der Grundflaeche, Width = Breite des Quaders
!
!
! Datei-Formate:
!   Input:   Shape (*.shp) oder Table (*.tbl), Verwendung mit Namens-Erweiterung
!            Table bietet Moeglichkeit des Auslassens von Objekten (ReadTable in Shape_Mod)
!   Output:  *.Haus  oder  *.Street  oder  *.Baum (Objekte jeweils nummeriert)
!
!
! ==========================================================================================

  CALL getarg(1,FileName)
  CALL getarg(2,FileNameGrid)
  WRITE(*,*) 'FileNameGrid  ',FileNameGrid
  CALL InputControl(FileNameGrid)
! Standard-Festlegung der Gebietseingrenzung P0 (links unten) und P1 (rechts oben), in Meter
! (--> keine Einschraenkung bei positiven Koordinaten < 1.d10):
  P0%x=0.0d0
  P1%x=1.d10
  P0%y=0.0d0
  P1%y=1.d10
! P0%x=x0DomainHaus
! P0%y=y0DomainHaus
! P1%x=x1DomainHaus
! P1%y=y1DomainHaus
! Standard-Festlegung der Skalierung Input-Koordinateneinheit --> Meter:
  Scale%x=1.d0
  Scale%y=1.d0
  write(*,*)
  WRITE(*,*) TRIM(FileName)
  write(*,*)

! Option: Shape or Table (taken from file extension)

  in_type  = TRIM(FileName(INDEX(FileName,'.')+1:INDEX(FileName,'.')+1))
  IF      (in_type=='s'.OR.in_type=='S') THEN
    in_type='Shape'
  ELSE IF (in_type=='t'.OR.in_type=='T') THEN
    in_type='Table'
  ELSE
    STOP '***** Stop: nor  *.shp  or  *.tbl *****'
  END IF
  FileName = TRIM(FileName(1:INDEX(FileName,'.')-1))
  write(*,*) 'in_type  = ',in_type
  write(*,*)

! Option Haus or Street, Input

  IF      (in_type=='Shape') THEN ! by inquiring
    write(*,*) 'Haus or Street or Baum ?'
    read(*,*) OutType
    OutType=TRIM(OutType)
    write(*,*) 'OutType = ',OutType
    write(*,*)
    IF (OutType/='Haus'.AND.OutType/='Street'.AND.OutType/='Baum') THEN
      STOP '***** Stop: nor Haus or Street or Baum *****'
    END IF
    CALL OpenFileShp(FileName)
    CALL ReadHeaderShp
    CALL OpenFileDbf(FileName)
    CALL ReadHeaderDbf(OutType)
  ELSE IF (in_type=='Table') THEN ! read from file
    CALL OpenFileTbl(FileName)
    CALL ReadHeaderTbl(OutType)
    write(*,*) 'OutType = ',OutType
    write(*,*)
    IF (OutType/='Haus'.AND.OutType/='Street'.AND.OutType/='Baum') THEN
      STOP '***** Stop: nor Haus or Street or Baum *****'
    END IF
  END IF
  write(*,*) '==========================================================================='
  write(*,*) 'THE FOLLOWING INPUT DATA CAN BE CHANGED IN THE MAIN PROGRAM'

!======================================================================================== INPUT for HAUS
  IF      (OutType=='Haus') THEN

    BaseIn=0.d0 ! standard base
    HeightIn=20.d0 ! standard height (above base)
    ValuesNum=1 ! do not change
    ALLOCATE(ValuesName(ValuesNum),ValuesCol(ValuesNum),ValuesIn(ValuesNum))
    ValuesCol=(/12/) ! Lutzener Straße
!   ValuesCol=(/7/)  ! Dresden Neustadt
    ! >0 (for Table): Base and Height taken from input file
    ! >0 (for Shape): = Position of Height (number of column) in the input file, Base=0
    ! =0:             Base and Height taken from standard values (BaseIn,HeightIn)

    CALL WriteHeaderHaus(in_type)
!======================================================================================== INPUT for STREET
  ELSE IF (OutType=='Street') THEN

    BaseIn=0.d0 ! standard base of emission layer
    HeightIn=2.d0 ! standard vertical extent of emission layer
    WidthIn=7.d0 ! standard width of emission layer (lanes)
    ValuesScale=1.d6/(365*24*3600) ! Factor for emission data from input file
                                   ! to convert g/ma --> µg/ms
    ValuesNum=2 ! Number of emission entries taken from input file
                ! (= number of species if each emission entry refers to another species)
    ALLOCATE(ValuesName(ValuesNum),ValuesCol(ValuesNum),ValuesIn(ValuesNum))
    ! These vectors have to correspond component-by-component
    ValuesName=(/'PM10','PM10'/)
    ! = species names (target names used in ASAM) for individual columns in input file
    !   (if several columns contribute to the same species --> multiple entry of name,
    !                                                          aggregation in output)
    ValuesIn=(/1.d0,0.d0/)
    ! = standard emission rates (µg/ms)
    ValuesCol=(/18,19/)
    ! >0 (for Table): = Position of species (number of column) in the input file,
    !                 Base,Height,Width also read from input file
    ! =0 (for Table): standard emission (ValuesIn) is used for that species,
    !                 Base,Height,Width read from input file
    ! >0 (for Shape): = Position of species (number of column) in the input file,
    !                 Base,Height,Width taken from standard values (BaseIn,HeightIn,WidthIn)
    ! =0 (for Shape): standard emission (ValuesIn) is used for that species,
    !                 Base,Height,Width taken from standard values (BaseIn,HeightIn,WidthIn)

    CALL WriteHeaderStreet(in_type)
!======================================================================================== INPUT for BAUM
  ELSE IF (OutType=='Baum') THEN

    BaseIn=0.d0 ! standard base of tree layer
    HeightIn=10.d0 ! standard vertical extent of tree layer
    WidthIn=5.d0 ! standard width of tree layer
    ValuesNum=1 ! do not change
    ALLOCATE(ValuesName(ValuesNum),ValuesCol(ValuesNum),ValuesIn(ValuesNum))
    ValuesIn=(/1.d0/) ! density of tree layer (0...1)
    ValuesCol=(/1/)
    ! >0 (for Table): Base,Height,Width and Density read from input file
    ! =0 (for Table): Base,Height,Width read from input file,
    !                 Density taken from standard values (ValuesIn)
    ! =0 (for Shape): All data taken from standard values (BaseIn,HeightIn,WidthIn,ValuesIn)
    ! >0 (for Shape): not realised

    CALL WriteHeaderBaum(in_type)
  END IF
!======================================================================================== INPUT for COORDINATES
!                                                                                               and BOUNDARIES
! Scaling factors for the input coordinates to transform them into meter
! (standard is set above to: 1.d0):
!   Scale%x=0.192870906d0 !(for Google-Earth coordinates)
!   Scale%y=0.308530319d0 !(for Google-Earth coordinates)

! Boundaries (P0,P1) to the coordinates with new origin P0 (in meter)
! (standard is set above to: P0%x=P0%y=0, P1%x=P1%y=1.d10)
! 
    !P0%x=4620490.0 !5410220.0 !+800.d0 !4621800.0
    !P1%x=4623330.0 !5413060.0 !-800.d0 !4622410.0
    !P0%y=5659290.0 !5658590.0 !+800.d0 !5660080.0
    !P1%y=5661340.0 !5660640.0 !-800.d0 !5660630.0

! Do also attend to different options in confining the Write loop at program end:
! either applying SR PartiallyInsideBox or TotallyInsideBox or not any
!
!======================================================================================== INPUT end
!
  write(*,*)
  write(*,*) 'Boundaries and scaling for coordinates:'
  write(*,4) 'P0 P1 Scal x:',P0%x,P1%x,Scale%x
  write(*,4) 'P0 P1 Scal y:',P0%y,P1%y,Scale%y
4 FORMAT(3x,a13,2f20.2,f10.5)
  WRITE(*,*) '==========================================================================='
  WRITE(*,*) 'Input correct ? y  or  n'
  READ(*,*) q
  IF (q=='n'.OR.q=='N') THEN
    WRITE(*,*) 'Please do change respective assignments in the Main program!'
    WRITE(*,*)
    STOP
  END IF

! READ

  IF (in_type=='Shape'.AND.MAXVAL(ValuesCol)>0) THEN
    write(*,*)
    write(*,*) '  Successive READ --> WRITE --> READ ...'
  END IF
  OPEN(UNIT=OutputUnitHaus,FILE=TRIM(FileName)//'.'//TRIM(OutType),STATUS='UNKNOWN')
  ALLOCATE(StartPoly)
  NumHaus=0
  CurrPoly=>StartPoly
  !WRITE(1,3) '  ','   NumHaus','    Number','  NumParts',' NumPoints'
3 FORMAT(a2,4a10)
  DO
    NumHaus=NumHaus+1
    IF      (in_type=='Shape') THEN
      CALL ReadPolygon(CurrPoly)
      IF (MAXVAL(ValuesCol)>0) THEN
        CALL ReadRecordDbf(RecordEnd,RecordErr)
        IF (RecordEnd.AND.RecordErr) THEN
          EXIT ! --> A
        ELSE IF (RecordErr) THEN
          RecordErr=.FALSE.
          CurrPoly%Number=-ABS(CurrPoly%Number)
          GOTO 10 ! --> D
        END IF
      END IF
      CALL CompletePolygon(CurrPoly,OutType)
10    CONTINUE
      IF (2*FileLengthShp<=nRec+4) THEN
        EXIT
      END IF
    ELSE IF (in_type=='Table') THEN
      CALL ReadTable(CurrPoly,OutType,RecordEnd,RecordErr)
      IF (RecordEnd) THEN
        EXIT ! --> C
      ELSE IF (RecordErr) THEN
        EXIT ! --> A
      END IF
      ! --> B
    END IF
    !WRITE(2,2) 'I',NumHaus,CurrPoly%Number,CurrPoly%NumParts,CurrPoly%NumPoints
2   FORMAT(a2,4i10)
!    IF (numhaus==3397) THEN ! Abbruch als Beispiel
!      exit
!    END IF
    ALLOCATE(NewPoly)
    CurrPoly%Next=>NewPoly
    CurrPoly=>NewPoly
  END DO
  NumHaus1=ABS(CurrPoly%Number)
  IF (RecordEnd) THEN    ! Objekt = Ende
    CurrPoly%NumPoints=0 ! Ende-Anzeiger
    NumHaus =NumHaus -1
    NumHaus1=NumHaus1-1
  END IF
  IF (in_type=='Shape') THEN
    CALL CloseFileDbf
  END IF
  write(*,*)
  IF (RecordErr) THEN
    write(*,*) 'NumHaus in    INPUT:     ',NumHaus1,OutType,'<-- ReadStop'
  ELSE
    write(*,*) 'NumHaus in    INPUT:     ',NumHaus1,OutType
  END IF
  write(*,*) 'NumHaus after READ:      ',NumHaus
  NumHaus1=NumHaus

! PARTITION

  NumHaus=0
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      NumHaus=NumHaus+CurrPoly%NumParts
      CALL PartitionPolygon(CurrPoly)
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO

  CALL WriteHausVTK(StartPoly,Filename)
  CALL WriteHausFaceVTK(StartPoly,Filename)
  write(*,*) 'NumHaus after PARTITION: ',NumHaus
  write(*,*)

! CHECK (coordinate shift, reorientation, revision)

  P1%x=P1%x*Scale%x
  P1%y=P1%y*Scale%y
  P0%x=P0%x*Scale%x
  P0%y=P0%y*Scale%y
  NumHaus=0
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      CALL Check(CurrPoly,P0,Scale,OutType,in_type)
      IF (CurrPoly%Number>0) THEN
        NumHaus=NumHaus+1
        !WRITE(1,2) 'C',NumHaus,CurrPoly%Number,CurrPoly%NumParts,CurrPoly%NumPoints
      END IF
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  write(*,*) 'NumHaus after CHECK:     ',NumHaus
  write(*,*)
  P1%x=P1%x -P0%x
  P1%y=P1%y -P0%y
  P0%x=0.d0
  P0%y=0.d0

! STREET and BAUM ignore Decompose and Simplify

  IF (OutType=='Street'.OR.OutType=='Baum') THEN
    GOTO 1
  END IF

! DECOMPOSE1 (coarse)

  NumHaus=0
  NumHaus1=0
  CurrPoly=>StartPoly
  DO
    Temp=>CurrPoly
    IF (CurrPoly%Number>0) THEN
      NumberT=CurrPoly%Number
      CALL DecomposePolygon(CurrPoly,NumberT,1.d-13)
      CurrPoly%Number=NumberT
      NumHaus=NumHaus+1
      IF (NumberT>0) THEN
        NumHaus1=NumHaus1+1
      END IF
      S1:DO
        Temp=>Temp%Next
        IF (.NOT.ASSOCIATED(Temp).OR.Temp%Number/=0) THEN ! Neu Hinzugekommene
          EXIT S1
        END IF
        NumHaus=NumHaus+1
        IF (NumberT>0) THEN
          NumHaus1=NumHaus1+1
        END IF
        Temp%Number=NumberT
      END DO S1
    ELSE
      Temp=>Temp%Next
    END IF
    IF (ASSOCIATED(Temp).AND.Temp%NumPoints>0) THEN
      CurrPoly=>Temp
    ELSE
      EXIT
    END IF
  END DO
  write(*,*) 'NumHaus by    DECOMPOSE1:',NumHaus
  write(*,*) 'NumHaus after DECOMPOSE1:',NumHaus1
  GO TO 1000

! SIMPLIFY1 (remove nearly lineal points)

  NumHaus=0
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      CALL SimplifyPolygon(CurrPoly,1.d-1)
      IF (CurrPoly%Number>0) THEN
        NumHaus=NumHaus+1
      END IF
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  write(*,*) 'NumHaus after SIMPLIFY1: ',NumHaus
  write(*,*)

! DECOMPOSE2 (fine)

  NumHaus=0
  NumHaus1=0
  CurrPoly=>StartPoly
  DO
    Temp=>CurrPoly
    IF (CurrPoly%Number>0) THEN
      NumberT=CurrPoly%Number
      CALL DecomposePolygon(CurrPoly,NumberT,1.d-5)
      CurrPoly%Number=NumberT
      NumHaus=NumHaus+1
      IF (NumberT>0) THEN
        NumHaus1=NumHaus1+1
      END IF
      S2:DO
        Temp=>Temp%Next
        IF (.NOT.ASSOCIATED(Temp).OR.Temp%Number/=0) THEN ! Neu Hinzugekommene
          EXIT S2
        END IF
        NumHaus=NumHaus+1
        IF (NumberT>0) THEN
          NumHaus1=NumHaus1+1
        END IF
        Temp%Number=NumberT
      END DO S2
    ELSE IF (CurrPoly%Number<0) THEN
      Temp=>Temp%Next
    END IF
    IF (ASSOCIATED(Temp).AND.Temp%NumPoints>0) THEN
      CurrPoly=>Temp
    ELSE
      EXIT
    END IF
  END DO
  write(*,*) 'NumHaus by    DECOMPOSE2:',NumHaus
  write(*,*) 'NumHaus after DECOMPOSE2:',NumHaus1

! SIMPLIFY2 (remove nearly lineal points)

  NumHaus=0
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      CALL SimplifyPolygon(CurrPoly,1.d-1)
      IF (CurrPoly%Number>0) THEN
        NumHaus=NumHaus+1
      END IF
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  write(*,*) 'NumHaus after SIMPLIFY2: ',NumHaus
  write(*,*)

1 CONTINUE
1000 CONTINUE

! STREET: RENUMBERING and WRITE of species in case of multiple entries

  IF (OutType=='Street') THEN
    SpeciesNum=1 ! numbers different species
    ValuesCol(1)=SpeciesNum
    DO i=2,ValuesNum ! comparing with the foregoing names
      DO l=1,i-1
        IF (ValuesName(i)==ValuesName(l)) THEN ! multiple entry found
          ValuesCol(i)=ValuesCol(l)
          EXIT
        END IF
        IF (l==i-1) THEN ! different species found
          SpeciesNum=SpeciesNum+1
          ValuesCol(i)=SpeciesNum ! contains number of (different) species
          ValuesName(SpeciesNum)=ValuesName(i) ! renumbered vector
        END IF
      END DO
    END DO
    WRITE(OutputUnitHaus,*) SpeciesNum,'!Species'
    WRITE(OutputUnitHaus,*) (ValuesName(i),' ',i=1,SpeciesNum)
  END IF

! COUNT and WRITE within BOUNDARIES (P0,P1)

  h_xmin= 1.d20
  h_xmax=-1.d20
  h_ymin= 1.d20
  h_ymax=-1.d20
  h_zmin= 1.d20
  h_zmax=-1.d20
  HeightIn=0.d0
  BaseIn  =0.d0
  WidthIn =0.d0
  ValuesIn=0.d0
  DO i=1,2 ! i=1 --> only count, i=2 --> count and write
    IF (i==2) THEN
      WRITE(OutputUnitHaus,*) NumHaus,'!',OutType
    END IF
    NumHaus=0
    CurrPoly=>StartPoly
    DO
      IF (CurrPoly%Number>0) THEN
        ! Einschränkung der Objekte nach ihrer Lage zwischen P0 und P1:
        !    TotallyInsideBox   --> grenzüberschreitende Objekte fallen heraus
        !    PartiallyInsideBox --> grenzüberschreitende Objekte berücksichtigt
        !    No condition       --> alle Objekte berücksichtigt
        ! In allen Fällen gilt: Vollkommen außerhalb liegende Objekte fallen heraus,
        !                       vollkommen innerhalb liegende Objekte berücksichtigt
        !
        WRITE(*,*) CurrPoly%Number,'Point 1 vor',CurrPoly%Points(1)
        IF (TotallyInsideBox(P0,P1,CurrPoly)) THEN ! condition is either TotallyInside ...
        !IF (PartiallyInsideBox(P0,P1,CurrPoly)) THEN ! ... or PartiallyInside ...
        !IF (0==0) THEN ! ... or not any
          NumHaus=NumHaus+1
          !WRITE(1,2) 'W',NumHaus,CurrPoly%Number,CurrPoly%NumParts,CurrPoly%NumPoints
          IF (i==2) THEN
            IF (OutType=='Haus') THEN
              WRITE(*,*) NumHaus,'Point 1',CurrPoly%Points(1)
!             CALL WriteHaus(CurrPoly,NumHaus)
              CALL WriteHaus1(CurrPoly,NumHaus)
            ELSE IF (OutType=='Haus1') THEN
              CALL WriteHaus1(CurrPoly,NumHaus)
            ELSE IF (OutType=='Street') THEN
              CALL WriteStreet(CurrPoly,NumHaus)
            ELSE IF (OutType=='Baum') THEN
              CALL WriteBaum(CurrPoly,NumHaus)
            END IF
            HeightIn   =MAX(CurrPoly%Height   ,HeightIn)
            BaseIn     =MAX(CurrPoly%Base     ,BaseIn)
            WidthIn    =MAX(CurrPoly%Width    ,WidthIn)
            ValuesIn(:)=MAX(CurrPoly%Values(:),ValuesIn(:))
          END IF
        END IF
      END IF
      IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
        CurrPoly=>CurrPoly%Next
      ELSE
        EXIT
      END IF
    END DO
  END DO
  WRITE(*,*) 'NumHaus in    BOUNDARIES:',NumHaus,'Written'
  WRITE(*,*)
  WRITE(*,*) 'Expansion x: ',h_xmin,h_xmax
  WRITE(*,*) '          y: ',h_ymin,h_ymax
  WRITE(*,*) '          z: ',h_zmin,h_zmax
  WRITE(*,7) 'Expansion x: ',h_xmin,h_xmax
  WRITE(*,7) '          y: ',h_ymin,h_ymax
  WRITE(*,7) '          z: ',h_zmin,h_zmax
  WRITE(*,7) 'Maxima       '
  IF      (OutType=='Haus') THEN
    WRITE(*,7) 'HausBase:    ',BaseIn
    WRITE(*,7) 'HausHeight:  ',HeightIn
  ELSE IF (OutType=='Street') THEN
    WRITE(*,7) 'StreetBase:  ',BaseIn
    WRITE(*,7) 'StreetHeight:',HeightIn
    WRITE(*,7) 'StreetWidth: ',WidthIn
    WRITE(*,5) 'SpeciesNum:  ',SpeciesNum
    DO i=1,SpeciesNum
      WRITE(*,6) ValuesName(i),' µg/ms ',ValuesIn(i)
    END DO
  ELSE IF (OutType=='Baum') THEN
    WRITE(*,7) 'BaumBase:    ',BaseIn
    WRITE(*,7) 'BaumHeight:  ',HeightIn
    WRITE(*,7) 'BaumWidth:   ',WidthIn
    WRITE(*,7) 'BaumDensity: ',ValuesIn(1)
  END IF
  WRITE(*,*)
  CLOSE(UNIT=22)
  CLOSE(OutputUnitHaus)
5 FORMAT(3x,a13,i15)
6 FORMAT(3x,a5,1x,a7,f16.7)
7 FORMAT(3x,a13,2f15.2,f10.5)

CONTAINS

FUNCTION TotallyInsideBox(P0Loc,P1Loc,Polygon)

  LOGICAL :: TotallyInsideBox
  TYPE(Polygon_T) :: Polygon
  TYPE(Point2_T) :: P0Loc,P1Loc
  REAL(8) :: eps=1.d-4
  INTEGER :: i

  TotallyInsideBox=.TRUE.
  DO i=0,Polygon%NumPoints-1
    IF (Polygon%Points(i)%x<P0Loc%x-eps.OR. &
        Polygon%Points(i)%y<P0Loc%y-eps.OR. &
        Polygon%Points(i)%x>P1Loc%x+eps.OR. &
        Polygon%Points(i)%y>P1Loc%y+eps) THEN ! partially outside --> exclude
      TotallyInsideBox=.FALSE.
      EXIT
    END IF
  END DO

END FUNCTION TotallyInsideBox

FUNCTION PartiallyInsideBox(P0Loc,P1Loc,Polygon)

  LOGICAL :: PartiallyInsideBox
  TYPE(Polygon_T) :: Polygon
  TYPE(Point2_T) :: P0Loc,P1Loc
  REAL(8) :: eps=1.d-4
  INTEGER :: i

  PartiallyInsideBox=.FALSE.
  DO i=0,Polygon%NumPoints-1
    IF (Polygon%Points(i)%x>P0Loc%x+eps.AND. &
        Polygon%Points(i)%y>P0Loc%y+eps.AND. &
        Polygon%Points(i)%x<P1Loc%x-eps.AND. &
        Polygon%Points(i)%y<P1Loc%y-eps) THEN ! partially inside --> include
      PartiallyInsideBox=.TRUE.
      EXIT
    END IF
  END DO

END FUNCTION PartiallyInsideBox

SUBROUTINE WriteHausVTK(StartPoly,FileName)
  TYPE(Polygon_T), POINTER :: StartPoly
  CHARACTER(*) :: FileName

  INTEGER :: OutputUnit=10

  INTEGER :: NumberOfPoints
  INTEGER :: NumHaus
  INTEGER :: i,j
  INTEGER :: iFaceOffsets
  TYPE(Polygon_T), POINTER :: CurrPoly

  CurrPoly=>StartPoly
  NumHaus=0
  NumberOfPoints=0
  DO
    IF (CurrPoly%Number>0) THEN
      NumHaus=NumHaus+1
      NumberOfPoints=NumberOfPoints+2*(CurrPoly%NumPoints-1)
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO

  OPEN(UNIT=OutputUnit,FILE=TRIM(FileName)//'Shape.vtu',STATUS='UNKNOWN')
  WRITE(OutputUnit,'(a21)') '<?xml version="1.0"?>'
  WRITE(OutputUnit,'(a33)') '<VTKFile type="UnstructuredGrid">'
  WRITE(OutputUnit,'(a18)') '<UnstructuredGrid>'
  WRITE(OutputUnit,'(a23,i8,a17,i8,a3)') '<Piece NumberOfPoints="',NumberOfPoints,'" NumberOfCells="',NumHaus,'">'
  WRITE(OutputUnit,'(a46)') '<CellData Scalars="scalars" Vector="Velocity">'
  WRITE(OutputUnit,*) '<DataArray type="Float64" Name="Number" Format="ascii">'
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      WRITE(OutputUnit,*) CurrPoly%Number
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a11)') '</CellData>'
  WRITE(OutputUnit,'(a8)') '<Points>'
  WRITE(OutputUnit,'(a64)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      DO j=1,CurrPoly%NumPoints-1
        WRITE(OutputUnit,*) CurrPoly%Points(j)%x-MoveXYGrid*xOffset,CurrPoly%Points(j)%y-MoveXYGrid*yOffset,CurrPoly%Base
      END DO  
      DO j=1,CurrPoly%NumPoints-1
        WRITE(OutputUnit,*) CurrPoly%Points(j)%x-MoveXYGrid*xOffset,CurrPoly%Points(j)%y-MoveXYGrid*yOffset,CurrPoly%Height
      END DO  
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a9)') '</Points>'
  WRITE(OutputUnit,'(a7)') '<Cells>'
  WRITE(OutputUnit,'(a59)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
  NumberOfPoints=0
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      WRITE(OutputUnit,*) (j+NumberOfPoints-1,j=1,CurrPoly%NumPoints-1) &
                         ,(j+NumberOfPoints-1+CurrPoly%NumPoints-1,j=1,CurrPoly%NumPoints-1)
      NumberOfPoints=NumberOfPoints+2*(CurrPoly%NumPoints-1)
                         
    END IF
    IF (ASSOCIATED(CurrPoly%Next)) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a54)') '<DataArray type="Int32" Name="offsets" format="ascii">'
  NumberOfPoints=0
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      NumberOfPoints=NumberOfPoints+2*(CurrPoly%NumPoints-1)
      WRITE(OutputUnit,*) NumberOfPoints
    END IF
    IF (ASSOCIATED(CurrPoly%Next)) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a52)') '<DataArray type="UInt8" Name="types" format="ascii">'
  DO i=1,NumHaus
    WRITE(OutputUnit,*) 42
  END DO  
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a52)') '<DataArray type="Int32" Name="faces" format="ascii">'
  NumberOfPoints=-1
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      WRITE(OutputUnit,*) CurrPoly%NumPoints+1
      DO i=1,CurrPoly%NumPoints-2
        WRITE(OutputUnit,*) 4 &
                           ,i+NumberOfPoints &
                           ,i+1+NumberOfPoints &
                           ,i+1+CurrPoly%NumPoints-1+NumberOfPoints &
                           ,i+CurrPoly%NumPoints-1+NumberOfPoints
      END DO
      i=CurrPoly%NumPoints-1
      WRITE(OutputUnit,*) 4 &
                           ,i+NumberOfPoints &
                           ,1+NumberOfPoints &
                           ,1+CurrPoly%NumPoints-1+NumberOfPoints &
                           ,i+CurrPoly%NumPoints-1+NumberOfPoints
      WRITE(OutputUnit,*) CurrPoly%NumPoints-1,(j+NumberOfPoints,j=1,CurrPoly%NumPoints-1) 
      WRITE(OutputUnit,*) CurrPoly%NumPoints-1,(j+NumberOfPoints+CurrPoly%NumPoints-1,j=1,CurrPoly%NumPoints-1)
      NumberOfPoints=NumberOfPoints+2*(CurrPoly%NumPoints-1)
    END IF
    IF (ASSOCIATED(CurrPoly%Next)) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a58)') '<DataArray type="Int32" Name="faceoffsets" format="ascii">'

  CurrPoly=>StartPoly
  iFaceOffsets=0
  DO
    IF (CurrPoly%Number>0) THEN
      iFaceOffsets=iFaceOffsets+1
      DO i=1,CurrPoly%NumPoints-2
        iFaceOffsets=iFaceOffsets+1+4
      END DO
      iFaceOffsets=iFaceOffsets+1+4
      iFaceOffsets=iFaceOffsets+1+CurrPoly%NumPoints-1
      iFaceOffsets=iFaceOffsets+1+CurrPoly%NumPoints-1
      WRITE(OutputUnit,*) iFaceOffsets
    END IF
    IF (ASSOCIATED(CurrPoly%Next)) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a8)') '</Cells>'
  WRITE(OutputUnit,'(a8)') '</Piece>'
  WRITE(OutputUnit,'(a19)') '</UnstructuredGrid>'
  WRITE(OutputUnit,'(a10)') '</VTKFile>'
  CLOSE(OutputUnit)
END SUBROUTINE WriteHausVTK

SUBROUTINE WriteHausFaceVTK(StartPoly,FileName)
  TYPE(Polygon_T), POINTER :: StartPoly
  CHARACTER(*) :: FileName

  INTEGER :: OutputUnit=10

  INTEGER :: NumberOfPoints
  INTEGER :: NumHaus
  INTEGER :: i,j
  INTEGER :: iFaceOffsets
  TYPE(Polygon_T), POINTER :: CurrPoly
  TYPE(Point2_T) :: P0,P1

  P0%x=x0DomainHaus
  P0%y=y0DomainHaus
  P1%x=x1DomainHaus
  P1%y=y1DomainHaus
  NumHaus=0
  NumberOfPoints=0
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0.AND.TotallyInsideBox(P0,P1,CurrPoly)) THEN
      NumHaus=NumHaus+1
      NumberOfPoints=NumberOfPoints+2*(CurrPoly%NumPoints-1)
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO

  OPEN(UNIT=OutputUnit,FILE=TRIM(FileName)//'Face.vtu',STATUS='UNKNOWN')
  WRITE(OutputUnit,'(a21)') '<?xml version="1.0"?>'
  WRITE(OutputUnit,'(a33)') '<VTKFile type="UnstructuredGrid">'
  WRITE(OutputUnit,'(a18)') '<UnstructuredGrid>'
  WRITE(OutputUnit,'(a23,i8,a17,i8,a3)') '<Piece NumberOfPoints="' &
                  ,NumberOfPoints,'" NumberOfCells="',NumHaus+NumberOfPoints/2,'">'
  WRITE(10,'(a28)') '<CellData Scalars="scalars">'
  WRITE(10,*) '<DataArray type="Float64" Name="'//'ColorFace'//'" Format="ascii">'
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0.AND.TotallyInsideBox(P0,P1,CurrPoly)) THEN
      DO i=1,CurrPoly%NumPoints-1
        WRITE(OutputUnit,*) 4 
      END DO
      WRITE(OutputUnit,*) 100
    END IF
    IF (ASSOCIATED(CurrPoly%Next)) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO  
  WRITE(10,'(a12)') '</DataArray>'
  WRITE(10,'(a11)') '</CellData>'

  WRITE(OutputUnit,'(a8)') '<Points>'
  WRITE(OutputUnit,'(a64)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0.AND.TotallyInsideBox(P0,P1,CurrPoly)) THEN
      DO j=1,CurrPoly%NumPoints-1
        WRITE(OutputUnit,*) CurrPoly%Points(j)%x-MoveXYGrid*xOffset,CurrPoly%Points(j)%y-MoveXYGrid*yOffset,CurrPoly%Base
      END DO  
      DO j=1,CurrPoly%NumPoints-1
        WRITE(OutputUnit,*) CurrPoly%Points(j)%x-MoveXYGrid*xOffset,CurrPoly%Points(j)%y-MoveXYGrid*yOffset,CurrPoly%Height
      END DO  
    END IF
    IF (ASSOCIATED(CurrPoly%Next).AND.CurrPoly%Next%NumPoints>0) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a9)') '</Points>'
  WRITE(OutputUnit,'(a7)') '<Cells>'
  WRITE(OutputUnit,'(a59)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
  NumberOfPoints=-1
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0.AND.TotallyInsideBox(P0,P1,CurrPoly)) THEN
      DO i=1,CurrPoly%NumPoints-2
        WRITE(OutputUnit,*) i+NumberOfPoints &
                           ,i+1+NumberOfPoints &
                           ,i+1+CurrPoly%NumPoints-1+NumberOfPoints &
                           ,i+CurrPoly%NumPoints-1+NumberOfPoints
      END DO
      i=CurrPoly%NumPoints-1
      WRITE(OutputUnit,*) i+NumberOfPoints &
                         ,1+NumberOfPoints &
                         ,1+CurrPoly%NumPoints-1+NumberOfPoints &
                         ,i+CurrPoly%NumPoints-1+NumberOfPoints
      WRITE(OutputUnit,*) (j+NumberOfPoints+CurrPoly%NumPoints-1,j=1,CurrPoly%NumPoints-1)
      NumberOfPoints=NumberOfPoints+2*(CurrPoly%NumPoints-1)
    END IF
    IF (ASSOCIATED(CurrPoly%Next)) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a54)') '<DataArray type="Int32" Name="offsets" format="ascii">'
  NumberOfPoints=0
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0.AND.TotallyInsideBox(P0,P1,CurrPoly)) THEN
      DO i=1,CurrPoly%NumPoints-1
        NumberOfPoints=NumberOfPoints+4
        WRITE(OutputUnit,*) NumberOfPoints
      END DO
      NumberOfPoints=NumberOfPoints+(CurrPoly%NumPoints-1)
      WRITE(OutputUnit,*) NumberOfPoints
    END IF
    IF (ASSOCIATED(CurrPoly%Next)) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a52)') '<DataArray type="UInt8" Name="types" format="ascii">'
  CurrPoly=>StartPoly
  DO
    IF (CurrPoly%Number>0) THEN
      DO i=1,CurrPoly%NumPoints
        WRITE(OutputUnit,*) 7
      END DO
    END IF
    IF (ASSOCIATED(CurrPoly%Next)) THEN
      CurrPoly=>CurrPoly%Next
    ELSE
      EXIT
    END IF
  END DO  
  WRITE(OutputUnit,'(a12)') '</DataArray>'
  WRITE(OutputUnit,'(a8)') '</Cells>'
  WRITE(OutputUnit,'(a8)') '</Piece>'
  WRITE(OutputUnit,'(a19)') '</UnstructuredGrid>'
  WRITE(OutputUnit,'(a10)') '</VTKFile>'
  CLOSE(OutputUnit)
END SUBROUTINE WriteHausFaceVTK

SUBROUTINE WriteHaus(PolygonSh,Number)

  TYPE(Polygon_T), POINTER :: PolygonSh
  INTEGER :: i,Number
  INTEGER :: NumEdges(0:PolygonSh%NumPoints)
  CHARACTER*2 :: TypeFaces(0:PolygonSh%NumPoints)

! Vom Punkt erst alle x, dann y, dann z

  DO i=0,PolygonSh%NumPoints-2
    NumEdges(i)=4
    TypeFaces(i)='w '
  END DO
  NumEdges(PolygonSh%NumPoints-1)=PolygonSh%NumPoints-1
  TypeFaces(PolygonSh%NumPoints-1)='r '
  IF (PolygonSh%Base<=0.d0) THEN
    WRITE(OutputUnitHaus,*) Number,PolygonSh%NumPoints,'!Number,Faces'
    WRITE(OutputUnitHaus,*) (NumEdges(i),i=0,PolygonSh%NumPoints-1),'!Points'
    WRITE(OutputUnitHaus,*) (TypeFaces(i),i=0,PolygonSh%NumPoints-1)
  ELSE
    NumEdges(PolygonSh%NumPoints)=PolygonSh%NumPoints-1
    TypeFaces(PolygonSh%NumPoints)='b '
    WRITE(OutputUnitHaus,*) Number,PolygonSh%NumPoints+1,'!Number,Faces'
    WRITE(OutputUnitHaus,*) (NumEdges(i),i=0,PolygonSh%NumPoints),'!Points'
    WRITE(OutputUnitHaus,*) (TypeFaces(i),i=0,PolygonSh%NumPoints)
  END IF
  DO i=0,PolygonSh%NumPoints-2
    WRITE(OutputUnitHaus,*) PolygonSh%Points(i)%x   &
                           ,PolygonSh%Points(i+1)%x &
                           ,PolygonSh%Points(i+1)%x &
                           ,PolygonSh%Points(i)%x   &
                           ,PolygonSh%Points(i)%y   &
                           ,PolygonSh%Points(i+1)%y &
                           ,PolygonSh%Points(i+1)%y &
                           ,PolygonSh%Points(i)%y   &
                           ,PolygonSh%Base        &
                           ,PolygonSh%Base        &
                           ,PolygonSh%Height+PolygonSh%Base        &
                           ,PolygonSh%Height+PolygonSh%Base,'!',TypeFaces(i)
    h_xmin=MIN(PolygonSh%Points(i)%x,PolygonSh%Points(i+1)%x,h_xmin)
    h_ymin=MIN(PolygonSh%Points(i)%y,PolygonSh%Points(i+1)%y,h_ymin)
    h_xmax=MAX(PolygonSh%Points(i)%x,PolygonSh%Points(i+1)%x,h_xmax)
    h_ymax=MAX(PolygonSh%Points(i)%y,PolygonSh%Points(i+1)%y,h_ymax)
    h_zmin=MIN(PolygonSh%Base,h_zmin)
    h_zmax=MAX(PolygonSh%Height+PolygonSh%Base,h_zmax)
  END DO
  IF (PolygonSh%Base>0.d0) THEN
    WRITE(OutputUnitHaus,*) (PolygonSh%Points(i)%x,i=PolygonSh%NumPoints-2,0,-1)   &
                           ,(PolygonSh%Points(i)%y,i=PolygonSh%NumPoints-2,0,-1)   &
                           ,(PolygonSh%Base,i=PolygonSh%NumPoints-2,0,-1),'!',TypeFaces(PolygonSh%NumPoints)
  END IF
  WRITE(OutputUnitHaus,*) (PolygonSh%Points(i)%x,i=0,PolygonSh%NumPoints-2)   &
                         ,(PolygonSh%Points(i)%y,i=0,PolygonSh%NumPoints-2)   &
                         ,(PolygonSh%Height+PolygonSh%Base,i=0,PolygonSh%NumPoints-2),'!',TypeFaces(PolygonSh%NumPoints-1)

END SUBROUTINE WriteHaus

SUBROUTINE WriteHaus1(PolygonSh,Number)

  TYPE(Polygon_T), POINTER :: PolygonSh
  INTEGER :: i,Number
  INTEGER :: NumEdges(0:PolygonSh%NumPoints)
  CHARACTER*2 :: TypeFaces(0:PolygonSh%NumPoints)

! Vom Punkt erst alle x, dann y, dann z

  DO i=0,PolygonSh%NumPoints-2
    NumEdges(i)=4
    TypeFaces(i)='w '
  END DO
  NumEdges(PolygonSh%NumPoints-1)=PolygonSh%NumPoints-1
  TypeFaces(PolygonSh%NumPoints-1)='r '

  WRITE(OutputUnitHaus,*) Number ! Number
  WRITE(OutputUnitHaus,*) 2*(PolygonSh%NumPoints-1)
  DO i=0,PolygonSh%NumPoints-2
    WRITE(OutputUnitHaus,*) PolygonSh%Points(i),PolygonSh%Base
  END DO  
  DO i=0,PolygonSh%NumPoints-2
    WRITE(OutputUnitHaus,*) PolygonSh%Points(i),PolygonSh%Base+PolygonSh%Height
    h_xmin=MIN(PolygonSh%Points(i)%x,h_xmin)
    h_ymin=MIN(PolygonSh%Points(i)%y,h_ymin)
    h_xmax=MAX(PolygonSh%Points(i)%x,h_xmax)
    h_ymax=MAX(PolygonSh%Points(i)%y,h_ymax)
    h_zmin=MIN(PolygonSh%Base,h_zmin)
    h_zmax=MAX(PolygonSh%Height+PolygonSh%Base,h_zmax)
  END DO  
  IF (PolygonSh%Base<=0.d0) THEN
    WRITE(OutputUnitHaus,*) PolygonSh%NumPoints,'!Number,Faces'
    WRITE(OutputUnitHaus,*) (NumEdges(i),i=0,PolygonSh%NumPoints-1),'!Points'
    WRITE(OutputUnitHaus,*) (TypeFaces(i),i=0,PolygonSh%NumPoints-1)
  ELSE
    NumEdges(PolygonSh%NumPoints)=PolygonSh%NumPoints-1
    TypeFaces(PolygonSh%NumPoints)='b '
    WRITE(OutputUnitHaus,*) Number,PolygonSh%NumPoints+1,'!Number,Faces'
    WRITE(OutputUnitHaus,*) (NumEdges(i),i=0,PolygonSh%NumPoints)
    WRITE(OutputUnitHaus,*) (TypeFaces(i),i=0,PolygonSh%NumPoints)
  END IF
  DO i=1,PolygonSh%NumPoints-2
    WRITE(OutputUnitHaus,*) i,i+1,i+1+PolygonSh%NumPoints-1,i+PolygonSh%NumPoints-1
  END DO
  i=PolygonSh%NumPoints-1
  WRITE(OutputUnitHaus,*) i,1,1+PolygonSh%NumPoints-1,i+PolygonSh%NumPoints-1
  IF (PolygonSh%Base>0.d0) THEN
    WRITE(OutputUnitHaus,*) (i,i=PolygonSh%NumPoints-1,1,-1) 
  END IF
  WRITE(OutputUnitHaus,*) (i+PolygonSh%NumPoints-1,i=1,PolygonSh%NumPoints-1) 

END SUBROUTINE WriteHaus1

SUBROUTINE WriteStreet(PolygonSh,Number)

  TYPE(Polygon_T), POINTER :: PolygonSh
  INTEGER :: i,j,Number
  INTEGER :: NumEdges(0:PolygonSh%NumPoints-2)
  CHARACTER*2 :: TypeFaces(0:PolygonSh%NumPoints-2)
  REAL(8) :: sum

! Vom Punkt erst alle x, dann y, dann z
  DO i=0,PolygonSh%NumPoints-2
    NumEdges(i)=2
    TypeFaces(i)='s '
  END DO
 !
 !
  WRITE(OutputUnitHaus,*) Number,PolygonSh%NumPoints-1,'!Number,Lines'
  WRITE(OutputUnitHaus,*) (NumEdges(i),i=0,PolygonSh%NumPoints-2),'!Points'
  WRITE(OutputUnitHaus,*) (TypeFaces(i),i=0,PolygonSh%NumPoints-2)
  DO i=0,PolygonSh%NumPoints-2
    WRITE(OutputUnitHaus,*) PolygonSh%Points(i)%x   &
                           ,PolygonSh%Points(i+1)%x &
                     !     ,PolygonSh%Points(i+1)%x &
                     !     ,PolygonSh%Points(i)%x   &
                           ,PolygonSh%Points(i)%y   &
                           ,PolygonSh%Points(i+1)%y &
                     !     ,PolygonSh%Points(i+1)%y &
                     !     ,PolygonSh%Points(i)%y   &
                           ,PolygonSh%Base        &
                           ,PolygonSh%Base,'!',TypeFaces(i)
                     !     ,Height                  &
                     !     ,Height
    h_xmin=MIN(PolygonSh%Points(i)%x,PolygonSh%Points(i+1)%x,h_xmin)
    h_ymin=MIN(PolygonSh%Points(i)%y,PolygonSh%Points(i+1)%y,h_ymin)
    h_xmax=MAX(PolygonSh%Points(i)%x,PolygonSh%Points(i+1)%x,h_xmax)
    h_ymax=MAX(PolygonSh%Points(i)%y,PolygonSh%Points(i+1)%y,h_ymax)
    h_zmin=MIN(PolygonSh%Base,h_zmin)
    h_zmax=MAX(PolygonSh%Base,h_zmax)
  END DO
!  WRITE(OutputUnitHaus,*) (PolygonSh%Points(i)%x,i=0,PolygonSh%NumPoints-2)   &
!                         ,(PolygonSh%Points(i)%y,i=0,PolygonSh%NumPoints-2)   &
!                         ,(Height,i=0,PolygonSh%NumPoints-2)
  WRITE(OutputUnitHaus,*) PolygonSh%Height,PolygonSh%Width,'!Height,Width'
  DO i=1,SpeciesNum
    sum=0.d0
    DO j=1,ValuesNum ! subsume multiple species entries
      IF (ValuesCol(j)==i) THEN
        sum=sum+PolygonSh%Values(j)
      END IF
    END DO
    PolygonSh%Values(i)=sum ! renumbering of the vector
  END DO
  WRITE(OutputUnitHaus,*) (PolygonSh%Values(i),i=1,SpeciesNum),'!Emiss µg/ms'

END SUBROUTINE WriteStreet

SUBROUTINE WriteBaum(PolygonSh,Number)

  TYPE(Polygon_T), POINTER :: PolygonSh
  INTEGER :: i,j,Number
  INTEGER :: NumEdges(0:PolygonSh%NumPoints-2)
  CHARACTER*2 :: TypeFaces(0:PolygonSh%NumPoints-2)

! Vom Punkt erst alle x, dann y, dann z
  DO i=0,PolygonSh%NumPoints-2
    NumEdges(i)=2
    TypeFaces(i)='t '
  END DO
 !
 !
  WRITE(OutputUnitHaus,*) Number,PolygonSh%NumPoints-1,'!Number,Lines'
  WRITE(OutputUnitHaus,*) (NumEdges(i),i=0,PolygonSh%NumPoints-2),'!Points'
  WRITE(OutputUnitHaus,*) (TypeFaces(i),i=0,PolygonSh%NumPoints-2)
  DO i=0,PolygonSh%NumPoints-2
    WRITE(OutputUnitHaus,*) PolygonSh%Points(i)%x   &
                           ,PolygonSh%Points(i+1)%x &
                     !     ,PolygonSh%Points(i+1)%x &
                     !     ,PolygonSh%Points(i)%x   &
                           ,PolygonSh%Points(i)%y   &
                           ,PolygonSh%Points(i+1)%y &
                     !     ,PolygonSh%Points(i+1)%y &
                     !     ,PolygonSh%Points(i)%y   &
                           ,PolygonSh%Base        &
                           ,PolygonSh%Base,'!',TypeFaces(i)
                     !     ,Height                  &
                     !     ,Height
    h_xmin=MIN(PolygonSh%Points(i)%x,PolygonSh%Points(i+1)%x,h_xmin)
    h_ymin=MIN(PolygonSh%Points(i)%y,PolygonSh%Points(i+1)%y,h_ymin)
    h_xmax=MAX(PolygonSh%Points(i)%x,PolygonSh%Points(i+1)%x,h_xmax)
    h_ymax=MAX(PolygonSh%Points(i)%y,PolygonSh%Points(i+1)%y,h_ymax)
    h_zmin=MIN(PolygonSh%Base,h_zmin)
    h_zmax=MAX(PolygonSh%Base,h_zmax)
  END DO
!  WRITE(OutputUnitHaus,*) (PolygonSh%Points(i)%x,i=0,PolygonSh%NumPoints-2)   &
!                         ,(PolygonSh%Points(i)%y,i=0,PolygonSh%NumPoints-2)   &
!                         ,(Height,i=0,PolygonSh%NumPoints-2)
  WRITE(OutputUnitHaus,*) PolygonSh%Height,PolygonSh%Width,'!Height,Width'
  WRITE(OutputUnitHaus,*) PolygonSh%Values(1),'!Density'

END SUBROUTINE WriteBaum

END PROGRAM Main
