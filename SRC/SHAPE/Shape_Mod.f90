MODULE Shape_Mod

  USE Geometry_Mod
  IMPLICIT NONE

  INTEGER, PARAMETER :: RealDouble=8

  TYPE Point2_T
   REAL(RealDouble) :: x                      !x Koordinate
   REAL(RealDouble) :: y                      !y Koordinate
  END TYPE Point2_T

  TYPE PolyLine_T
    INTEGER :: ContentLength=0
    INTEGER :: ShapeType=3
    REAL(RealDouble) :: xMin,yMin,xMax,yMax   !Bounding Box
    INTEGER :: NumParts                       !Number of Parts
    INTEGER :: NumPoints                      !Total Number of Points
    INTEGER , ALLOCATABLE :: Parts(:)         !Index to First Point in Part
    TYPE(Point2_T), ALLOCATABLE  :: Points(:) !Points for All Parts
    TYPE(PolyLine_T), POINTER :: Next=>NULL()
  END TYPE PolyLine_T

  TYPE(PolyLine_T), POINTER :: ListOfPoLine=>NULL()

  TYPE Polygon_T
    INTEGER :: Number=0
    INTEGER :: ContentLength=0
    INTEGER :: ShapeType=5
    REAL(RealDouble) :: xMin,yMin,xMax,yMax   !Bounding Box
    REAL(RealDouble) :: Base  =0.d0           !dh
    REAL(RealDouble) :: Height=0.d0           !dh
    REAL(RealDouble) :: Width =0.d0           !dh
    REAL(RealDouble), ALLOCATABLE :: Values(:)!dh
    INTEGER :: NumParts                       !Number of Parts
    INTEGER :: NumPoints                      !Total Number of Points
    INTEGER , ALLOCATABLE :: Parts(:)         !Index to First Point in Part
    TYPE(Point2_T), ALLOCATABLE  :: Points(:) !Points for All Parts
    TYPE(Polygon_T), POINTER :: Next=>NULL()
  END TYPE Polygon_T


  TYPE PolygonZ_T
    INTEGER :: ContentLength=0
    INTEGER :: ShapeType=15
    REAL(RealDouble) :: xMin,yMin,xMax,yMax   !Bounding Box
    INTEGER :: NumParts                       !Number of Parts
    INTEGER :: NumPoints                      !Total Number of Points
    INTEGER , ALLOCATABLE :: Parts(:)         !Index to First Point in Part
    TYPE(Point2_T), ALLOCATABLE  :: Points(:)                !Points for All Parts
    REAL(RealDouble) :: zMin,zMax             !Bounding Z Range
    REAL(RealDouble), ALLOCATABLE  :: zPoints(:)            !Z Values for All Points
    REAL(RealDouble) :: MMin,MMax             !Bounding Measure Range
    REAL(RealDouble), ALLOCATABLE  :: MPoints(:)            !Measures
  END TYPE PolygonZ_T

  TYPE PolyLineZ_T
    INTEGER :: ContentLength=0
    INTEGER :: ShapeType=13
    REAL(RealDouble) :: xMin,yMin,xMax,yMax   !Bounding Box
    INTEGER :: NumParts                       !Number of Parts
    INTEGER :: NumPoints                      !Total Number of Points
    INTEGER , ALLOCATABLE :: Parts(:)         !Index to First Point in Part
    TYPE(Point2_T), POINTER  :: Points(:)=>NULL()                !Points for All Parts
    REAL(RealDouble) :: zMin,zMax             !Bounding Z Range
    REAL(RealDouble), ALLOCATABLE  :: zPoints(:)            !Z Values for All Points
    REAL(RealDouble) :: MMin,MMax             !Bounding Measure Range
    REAL(RealDouble), ALLOCATABLE  :: MPoints(:)            !Measures
    TYPE(PolyLineZ_T), POINTER :: Next=>NULL()
  END TYPE PolyLineZ_T

  TYPE(PolyLineZ_T), POINTER :: ListOfPoLiZ=>NULL()


  INTEGER, PARAMETER :: MaxCol=100 !dh
  INTEGER :: OutputUnit=30
  INTEGER :: nRec=-3
  INTEGER :: iRec=4
  INTEGER :: RecordNumber=0
  INTEGER :: FileCode=9994
  INTEGER :: FileLengthShp=50
  INTEGER :: FileLengthShx=50
  INTEGER :: Version=1000
  INTEGER :: ShapeType=0
  INTEGER :: Offset=50
  REAL(RealDouble) :: xMin,yMin,xMax,yMax   !Global Bounding Box
  REAL(RealDouble) :: zMin=.0d0,zMax=0.0d0  !Global zBounding Box
  REAL(RealDouble) :: MMin=.0d0,MMax=0.0d0  !Global MBounding Box
  INTEGER :: NumberPolygon=0

  TYPE(Point2_T) :: P0,P1,Scale
  REAL(RealDouble) :: HeightIn,BaseIn,WidthIn
  INTEGER ::          ValuesNum
  CHARACTER(5)    , ALLOCATABLE :: ValuesName(:)
  REAL(RealDouble), ALLOCATABLE :: ValuesIn(:)
  INTEGER         , ALLOCATABLE :: ValuesCol(:)
  CHARACTER(10) :: ColName(MaxCol) ! from DBase_Mod
  REAL(RealDouble) :: ColValue(MaxCol)
  INTEGER ::       ncols=0 ! from DBase_Mod
  REAL(RealDouble) :: ValuesScale=1.d0

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE CopyPolyLineZ &
                    ,CopyPointPoint2 &
                    ,CopyPoint2Point
  END INTERFACE
  
CONTAINS
SUBROUTINE CopyPointPoint2(P,P2)
  TYPE(Point_T), INTENT(OUT) :: P
  TYPE(Point2_T), INTENT(IN) :: P2

  P%x=P2%x
  P%y=P2%y
  P%z=0.0d0
END SUBROUTINE CopyPointPoint2

SUBROUTINE CopyPoint2Point(P2,P)
  TYPE(Point2_T), INTENT(OUT) :: P2
  TYPE(Point_T), INTENT(IN) :: P

  P2%x=P%x
  P2%y=P%y
END SUBROUTINE CopyPoint2Point

SUBROUTINE OpenFileShp(FileName)
  CHARACTER(*) :: FileName
! OPEN(UNIT=OutputUnit,FILE=TRIM(FileName)//'.shp',STATUS='UNKNOWN'&
!   &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)
  OPEN(UNIT=OutputUnit,FILE=TRIM(FileName)//'.shp',STATUS='UNKNOWN' &
      ,ACCESS='STREAM')
  nRec=-3
END SUBROUTINE OpenFileShp

SUBROUTINE CloseFileShp
  CLOSE(OutputUnit)
END SUBROUTINE CloseFileShp

SUBROUTINE OpenFileShx(FileName)
  CHARACTER(*) :: FileName
! OPEN(UNIT=OutputUnit,FILE=TRIM(FileName)//'.shx',STATUS='UNKNOWN'&
!   &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=1)
  OPEN(UNIT=OutputUnit,FILE=TRIM(FileName)//'.shx',STATUS='UNKNOWN' &
      ,ACCESS='STREAM')

  nRec=-3
END SUBROUTINE OpenFileShx

SUBROUTINE CloseFileShx
  CLOSE(OutputUnit)
END SUBROUTINE CloseFileShx

SUBROUTINE WriteInteger(I)
  INTEGER :: I

  nRec=nRec+iRec
  WRITE(OutputUnit,POS=nRec) I
END SUBROUTINE WriteInteger

SUBROUTINE WriteIntegerB(I)
  INTEGER :: I

  INTEGER :: IB
  CALL MVBITS(I,24,8,IB, 0)
  CALL MVBITS(I,16,8,IB, 8)
  CALL MVBITS(I, 8,8,IB,16)
  CALL MVBITS(I, 0,8,IB,24)

  nRec=nRec+iRec
  WRITE(OutputUnit,POS=nRec) IB
END SUBROUTINE WriteIntegerB

SUBROUTINE WriteRealDouble(R)
  REAL(RealDouble) :: R

  INTEGER :: iReal(1:2)
  iReal=TRANSFER(R,(/0,0/))
  nRec=nRec+iRec
  WRITE(OutputUnit,POS=nRec) iReal(1)
  nRec=nRec+iRec
  WRITE(OutputUnit,POS=nRec) iReal(2)
END SUBROUTINE WriteRealDouble

SUBROUTINE WritePoint(P)

  TYPE(Point2_T) :: P
  CALL WriteRealDouble(P%x)
  CALL WriteRealDouble(P%y)

END SUBROUTINE WritePoint

SUBROUTINE WriteHeaderShp
  CALL WriteIntegerB(FileCode)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(FileLengthShp)
  CALL WriteInteger(Version)
  CALL WriteInteger(ShapeType)
  CALL WriteRealDouble(xMin)
  CALL WriteRealDouble(yMin)
  CALL WriteRealDouble(xMax)
  CALL WriteRealDouble(yMax)
  CALL WriteRealDouble(zMin)
  CALL WriteRealDouble(zMax)
  CALL WriteRealDouble(MMin)
  CALL WriteRealDouble(MMax)
END SUBROUTINE WriteHeaderShp

SUBROUTINE WriteHeaderShx
  CALL WriteIntegerB(FileCode)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(0)
  CALL WriteIntegerB(FileLengthShx)
  CALL WriteInteger(Version)
  CALL WriteInteger(ShapeType)
  CALL WriteRealDouble(xMin)
  CALL WriteRealDouble(yMin)
  CALL WriteRealDouble(xMax)
  CALL WriteRealDouble(yMax)
  CALL WriteRealDouble(zMin)
  CALL WriteRealDouble(zMax)
  CALL WriteRealDouble(MMin)
  CALL WriteRealDouble(MMax)
END SUBROUTINE WriteHeaderShx

SUBROUTINE WriteRecordShx(Offset,ContentLength)
  INTEGER :: Offset,ContentLength
  CALL WriteIntegerB(Offset)
  CALL WriteIntegerB(ContentLength)
END SUBROUTINE WriteRecordShx


SUBROUTINE WritePolygonZ(PolygonZ)

  TYPE(Polygonz_T) :: PolygonZ
  INTEGER :: i,ContentLength
  
! Header
  RecordNumber=RecordNumber+1
  CALL WriteIntegerB(RecordNumber)
  CALL WriteIntegerB(PolygonZ%ContentLength)
! Shape type
  CALL WriteInteger(PolygonZ%ShapeType)
! Bounding Box
  CALL WriteRealDouble(PolygonZ%xMin)
  CALL WriteRealDouble(PolygonZ%yMin)
  CALL WriteRealDouble(PolygonZ%xMax)
  CALL WriteRealDouble(PolygonZ%yMax)
! Number of Parts
  CALL WriteInteger(Polygonz%NumParts)
! Total Number of Points
  CALL WriteInteger(Polygonz%NumPoints)
! Parts
  DO i=0,Polygonz%NumParts-1
    CALL WriteInteger(Polygonz%Parts(i))
  END DO
! Points
  DO i=0,Polygonz%NumPoints-1
    CALL WritePoint(Polygonz%Points(i))
  END DO
! z Bounding Box
  CALL WriteRealDouble(PolygonZ%zMin)
  CALL WriteRealDouble(PolygonZ%zMax)
! z Points
  DO i=0,Polygonz%NumPoints-1
    CALL WriteRealDouble(Polygonz%zPoints(i))
  END DO
  IF (ALLOCATED(Polygonz%MPoints)) THEN
!   M Bounding Box
    CALL WriteRealDouble(PolygonZ%MMin)
    CALL WriteRealDouble(PolygonZ%MMax)
!   M Points
    DO i=0,Polygonz%NumPoints-1
      CALL WriteRealDouble(PolygonZ%MPoints(i))
    END DO
  END IF
END SUBROUTINE WritePolygonZ

SUBROUTINE CopyPolyLineZ(PolylineZ2,PolylineZ1)
  TYPE(PolyLineZ_T), INTENT(OUT) :: PolylineZ2
  TYPE(PolyLineZ_T), INTENT(IN) :: PolylineZ1

  PolylineZ2%ContentLength=PolylineZ1%ContentLength
  PolylineZ2%xMin=PolylineZ1%xMin
  PolylineZ2%yMin=PolylineZ1%yMin
  PolylineZ2%xMax=PolylineZ1%xMax
  PolylineZ2%yMax=PolylineZ1%yMax
  PolylineZ2%NumParts=PolylineZ1%NumParts
  PolylineZ2%NumPoints=PolylineZ1%NumPoints
  PolylineZ2%zMin=PolylineZ1%zMin
  PolylineZ2%zMax=PolylineZ1%zMax
  PolylineZ2%MMin=PolylineZ1%MMin
  PolylineZ2%MMax=PolylineZ1%MMax
  !Parts(:) ....................................
  IF (ALLOCATED(PolylineZ2%Parts)) THEN
    DEALLOCATE(PolylineZ2%Parts)
  END IF
  ALLOCATE(PolylineZ2%Parts(LBOUND(PolylineZ1%Parts,1):UBOUND(PolylineZ1%Parts,1)))
  PolylineZ2%Parts=PolylineZ1%Parts
  !Points(:).....................................
  IF (ASSOCIATED(PolylineZ2%Points)) THEN
    DEALLOCATE(PolylineZ2%Points)
  END IF
  ALLOCATE(PolylineZ2%Points(LBOUND(PolylineZ1%Points,1):UBOUND(PolylineZ1%Points,1)))
  PolylineZ2%Points=PolylineZ1%Points
  !zPoints(:)....................................
  IF (ALLOCATED(PolylineZ2%zPoints)) THEN
    DEALLOCATE(PolylineZ2%zPoints)
  END IF
  ALLOCATE(PolylineZ2%zPoints(LBOUND(PolylineZ1%zPoints,1):UBOUND(PolylineZ1%zPoints,1)))
  PolylineZ2%zPoints=PolylineZ1%zPoints 
  !MPoints(:).....................................
  IF (ALLOCATED(PolylineZ2%MPoints)) THEN
    DEALLOCATE(PolylineZ2%MPoints)
  END IF
  IF (ALLOCATED(PolylineZ1%MPoints)) THEN
    ALLOCATE(PolylineZ2%MPoints(LBOUND(PolylineZ1%MPoints,1):UBOUND(PolylineZ1%MPoints,1)))
    PolylineZ2%MPoints=PolylineZ1%MPoints
  END IF
  !Next...........................................
! PolylineZ2%Next=>PolylineZ1%Next

END SUBROUTINE CopyPolyLineZ


SUBROUTINE WritePolyLineZ(PolylineZ)

  TYPE(PolyLineZ_T) :: PolylineZ
  INTEGER :: i,ContentLength
  
! Header
  RecordNumber=RecordNumber+1
  CALL WriteIntegerB(RecordNumber)
  CALL WriteIntegerB(PolylineZ%ContentLength)
! Shape type
  CALL WriteInteger(PolylineZ%ShapeType)
! Bounding Box
  CALL WriteRealDouble(PolylineZ%xMin)
  CALL WriteRealDouble(PolylineZ%yMin)
  CALL WriteRealDouble(PolylineZ%xMax)
  CALL WriteRealDouble(PolylineZ%yMax)
! Number of Parts
  CALL WriteInteger(PolylineZ%NumParts)
! Total Number of Points
  CALL WriteInteger(PolylineZ%NumPoints)
! Parts
  DO i=0,PolylineZ%NumParts-1
    CALL WriteInteger(PolylineZ%Parts(i))
  END DO
! Points
  DO i=0,PolylineZ%NumPoints-1
    CALL WritePoint(PolylineZ%Points(i))
  END DO
! z Bounding Box
  CALL WriteRealDouble(PolylineZ%zMin)
  CALL WriteRealDouble(PolylineZ%zMax)
! z Points
  DO i=0,PolylineZ%NumPoints-1
    CALL WriteRealDouble(PolylineZ%zPoints(i))
  END DO
  IF (ALLOCATED(PolylineZ%MPoints)) THEN
!   M Bounding Box
    CALL WriteRealDouble(PolylineZ%MMin)
    CALL WriteRealDouble(PolylineZ%MMax)
!   M Points
    DO i=0,PolylineZ%NumPoints-1
      CALL WriteRealDouble(PolylineZ%MPoints(i))
    END DO
  END IF
END SUBROUTINE WritePolyLineZ

SUBROUTINE WritePolyLine(Polyline)

  TYPE(PolyLine_T) :: Polyline
  INTEGER :: i,ContentLength
  
! Header
  RecordNumber=RecordNumber+1
  CALL WriteIntegerB(RecordNumber)
  CALL WriteIntegerB(Polyline%ContentLength)
! Shape type
  CALL WriteInteger(Polyline%ShapeType)
! Bounding Box
  CALL WriteRealDouble(Polyline%xMin)
  CALL WriteRealDouble(Polyline%yMin)
  CALL WriteRealDouble(Polyline%xMax)
  CALL WriteRealDouble(Polyline%yMax)
! Number of Parts
  CALL WriteInteger(Polyline%NumParts)
! Total Number of Points
  CALL WriteInteger(Polyline%NumPoints)
! Parts
  DO i=0,Polyline%NumParts-1
    CALL WriteInteger(Polyline%Parts(i))
  END DO
! Points
  DO i=0,Polyline%NumPoints-1
    CALL WritePoint(Polyline%Points(i))
  END DO
END SUBROUTINE WritePolyLine


SUBROUTINE ReadInteger(I)
  INTEGER :: I

  nRec=nRec+iRec
  READ(OutputUnit,POS=nRec) I
END SUBROUTINE ReadInteger

SUBROUTINE ReadIntegerB(I)
  INTEGER :: I

  INTEGER :: IB

  nRec=nRec+iRec
  READ(OutputUnit,POS=nRec) IB
! CALL MVBITS(IB,24,8,I, 0)
! CALL MVBITS(IB,16,8,I, 8)
! CALL MVBITS(IB, 8,8,I,16)
! CALL MVBITS(IB, 0,8,I,24)

  CALL MVBITS(IB, 0,8,I,24)
  CALL MVBITS(IB, 8,8,I,16)
  CALL MVBITS(IB,16,8,I, 8)
  CALL MVBITS(IB,24,8,I, 0)
END SUBROUTINE ReadIntegerB

SUBROUTINE ReadRealDouble(R)
  REAL(RealDouble) :: R

  INTEGER :: iReal(1:2)
  nRec=nRec+iRec
  READ(OutputUnit,POS=nRec) iReal(1)
  nRec=nRec+iRec
  READ(OutputUnit,POS=nRec) iReal(2)
  R=TRANSFER(iReal,0.0d0)
END SUBROUTINE ReadRealDouble

SUBROUTINE ReadPoint(P)

  TYPE(Point2_T) :: P
  CALL ReadRealDouble(P%x)
  CALL ReadRealDouble(P%y)

END SUBROUTINE ReadPoint

SUBROUTINE ReadHeaderShp
  INTEGER :: Dummy
  CALL ReadIntegerB(FileCode)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(FileLengthShp)
  CALL ReadInteger(Version)
  CALL ReadInteger(ShapeType)
  CALL ReadRealDouble(xMin)
  CALL ReadRealDouble(yMin)
  CALL ReadRealDouble(xMax)
  CALL ReadRealDouble(yMax)
  CALL ReadRealDouble(zMin)
  CALL ReadRealDouble(zMax)
  CALL ReadRealDouble(MMin)
  CALL ReadRealDouble(MMax)
END SUBROUTINE ReadHeaderShp

SUBROUTINE ReadHeaderShx
  INTEGER :: Dummy
  CALL ReadIntegerB(FileCode)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(Dummy)
  CALL ReadIntegerB(FileLengthShx)
  CALL ReadInteger(Version)
  CALL ReadInteger(ShapeType)
  CALL ReadRealDouble(xMin)
  CALL ReadRealDouble(yMin)
  CALL ReadRealDouble(xMax)
  CALL ReadRealDouble(yMax)
  CALL ReadRealDouble(zMin)
  CALL ReadRealDouble(zMax)
  CALL ReadRealDouble(MMin)
  CALL ReadRealDouble(MMax)
END SUBROUTINE ReadHeaderShx

SUBROUTINE ReadRecordShx(Offset,ContentLength)
  INTEGER :: Offset,ContentLength
  CALL ReadIntegerB(Offset)
  CALL ReadIntegerB(ContentLength)
END SUBROUTINE ReadRecordShx

SUBROUTINE ReadPolyLineZ(PolylineZ)

  TYPE(PolyLineZ_T) :: PolylineZ
  INTEGER :: i,ContentLength
  
! Header
  CALL ReadIntegerB(RecordNumber)
  CALL ReadIntegerB(PolylineZ%ContentLength)
! Shape type
  CALL ReadInteger(PolylineZ%ShapeType)
! Bounding Box
  CALL ReadRealDouble(PolylineZ%xMin)
  CALL ReadRealDouble(PolylineZ%yMin)
  CALL ReadRealDouble(PolylineZ%xMax)
  CALL ReadRealDouble(PolylineZ%yMax)
! Number of Parts
  CALL ReadInteger(PolylineZ%NumParts)
! Total Number of Points
  CALL ReadInteger(PolylineZ%NumPoints)
! Parts
  ALLOCATE(PolylineZ%Parts(0:PolylineZ%NumParts-1))
  DO i=0,PolylineZ%NumParts-1
    CALL ReadInteger(PolylineZ%Parts(i))
  END DO
! Points
  ALLOCATE(PolylineZ%Points(0:PolylineZ%NumPoints-1))
  DO i=0,PolylineZ%NumPoints-1
    CALL ReadPoint(PolylineZ%Points(i))
  END DO
! z Bounding Box
  CALL ReadRealDouble(PolylineZ%zMin)
  CALL ReadRealDouble(PolylineZ%zMax)
! z Points
  ALLOCATE(PolylineZ%zPoints(0:PolylineZ%NumPoints-1))
  DO i=0,PolylineZ%NumPoints-1
    CALL ReadRealDouble(PolylineZ%zPoints(i))
  END DO
  STOP
  IF (ALLOCATED(PolylineZ%MPoints)) THEN
!   M Bounding Box
    CALL ReadRealDouble(PolylineZ%MMin)
    CALL ReadRealDouble(PolylineZ%MMax)
!   M Points
    DO i=0,PolylineZ%NumPoints-1
      CALL ReadRealDouble(PolylineZ%MPoints(i))
    END DO
  END IF
END SUBROUTINE ReadPolyLineZ

SUBROUTINE ReadPolyLine(Polyline)

  TYPE(PolyLine_T) :: Polyline
  INTEGER :: i,ContentLength
  
! Header
  CALL ReadIntegerB(RecordNumber)
  CALL ReadIntegerB(Polyline%ContentLength)
! Shape type
  CALL ReadInteger(Polyline%ShapeType)
! Bounding Box
  CALL ReadRealDouble(Polyline%xMin)
  CALL ReadRealDouble(Polyline%yMin)
  CALL ReadRealDouble(Polyline%xMax)
  CALL ReadRealDouble(Polyline%yMax)
! Number of Parts
  CALL ReadInteger(Polyline%NumParts)
! Total Number of Points
  CALL ReadInteger(Polyline%NumPoints)
! Parts
  DO i=0,Polyline%NumParts-1
    CALL ReadInteger(Polyline%Parts(i))
  END DO
! Points
  DO i=0,Polyline%NumPoints-1
    CALL ReadPoint(Polyline%Points(i))
  END DO
END SUBROUTINE ReadPolyLine

SUBROUTINE ReadPolygonZ(PolygonZ)

  TYPE(Polygonz_T) :: PolygonZ
  INTEGER :: i,ContentLength
  
! Header
  CALL ReadIntegerB(RecordNumber)
  CALL ReadIntegerB(PolygonZ%ContentLength)
! Shape type
  CALL ReadInteger(PolygonZ%ShapeType)
! Bounding Box
  CALL ReadRealDouble(PolygonZ%xMin)
  CALL ReadRealDouble(PolygonZ%yMin)
  CALL ReadRealDouble(PolygonZ%xMax)
  CALL ReadRealDouble(PolygonZ%yMax)
! Number of Parts
  CALL ReadInteger(Polygonz%NumParts)
! Total Number of Points
  CALL ReadInteger(Polygonz%NumPoints)
! Parts
  ALLOCATE(Polygonz%Parts(0:Polygonz%NumParts-1))
  DO i=0,Polygonz%NumParts-1
    CALL ReadInteger(Polygonz%Parts(i))
  END DO
! Points
  DO i=0,Polygonz%NumPoints-1
    CALL ReadPoint(Polygonz%Points(i))
  END DO
! z Bounding Box
  CALL ReadRealDouble(PolygonZ%zMin)
  CALL ReadRealDouble(PolygonZ%zMax)
! z Points
  DO i=0,Polygonz%NumPoints-1
    CALL ReadRealDouble(Polygonz%zPoints(i))
  END DO
  IF (ALLOCATED(Polygonz%MPoints)) THEN
!   M Bounding Box
    CALL ReadRealDouble(PolygonZ%MMin)
    CALL ReadRealDouble(PolygonZ%MMax)
!   M Points
    DO i=0,Polygonz%NumPoints-1
      CALL ReadRealDouble(PolygonZ%MPoints(i))
    END DO
  END IF
END SUBROUTINE ReadPolygonZ

SUBROUTINE ReadPolygon(Polygon)

  TYPE(Polygon_T) :: Polygon
  INTEGER :: i,ContentLength
  
  NumberPolygon=NumberPolygon+1
  Polygon%Number=NumberPolygon
! Header
  CALL ReadIntegerB(RecordNumber)
  CALL ReadIntegerB(Polygon%ContentLength)
! Shape type
  CALL ReadInteger(Polygon%ShapeType)
! Bounding Box
  CALL ReadRealDouble(Polygon%xMin)
  CALL ReadRealDouble(Polygon%yMin)
  CALL ReadRealDouble(Polygon%xMax)
  CALL ReadRealDouble(Polygon%yMax)
! Number of Parts
  CALL ReadInteger(Polygon%NumParts)
! Total Number of Points
  CALL ReadInteger(Polygon%NumPoints)
! Parts
  ALLOCATE(Polygon%Parts(0:Polygon%NumParts-1))
  DO i=0,Polygon%NumParts-1
    CALL ReadInteger(Polygon%Parts(i))
  END DO
! Points
  ALLOCATE(Polygon%Points(0:Polygon%NumPoints-1))
  DO i=0,Polygon%NumPoints-1
    CALL ReadPoint(Polygon%Points(i))
  END DO
! Species
  ALLOCATE(Polygon%Values(ValuesNum)) !dh
END SUBROUTINE ReadPolygon

SUBROUTINE CompletePolygon(Polygon,out_type)
  TYPE(Polygon_T) :: Polygon
  CHARACTER(10) :: out_type

  IF      (out_type=='Haus') THEN
    IF (MAXVAL(ValuesCol)>0) THEN
      BaseIn=0.d0 ! Attention: not generally
      HeightIn=ValuesIn(1)
      WRITE(*,*) 'ValuesIn(1)',ValuesIn(1)
    END IF
    Polygon%Base  =BaseIn
    Polygon%Height=HeightIn
  ELSE IF (out_type=='Street') THEN
    Polygon%Base  =BaseIn
    Polygon%Height=HeightIn
    Polygon%Width =WidthIn
    Polygon%Values=ValuesIn
  ELSE IF (out_type=='Baum') THEN ! Shape not realized
    Polygon%Base  =BaseIn
    Polygon%Height=HeightIn
    Polygon%Width =WidthIn
    Polygon%Values=ValuesIn
  END IF
END SUBROUTINE CompletePolygon

SUBROUTINE OpenFileTbl(FileName) !dh
  CHARACTER(*) :: FileName

  OPEN(UNIT=OutputUnit,FILE=TRIM(FileName)//'.tbl',STATUS='OLD')
END SUBROUTINE OpenFileTbl

SUBROUTINE ReadHeaderTbl(out_type) !dh
  CHARACTER(10) :: out_type
  INTEGER :: i

  READ(OutputUnit,*) out_type
  out_type=TRIM(out_type)
  IF      (out_type=='Haus') THEN
    ncols=1
    ColName(1)='height'
  ELSE IF (out_type=='Street') THEN
    READ(OutputUnit,*) ncols
    READ(OutputUnit,*) (ColName(i),i=1,ncols)
  ELSE IF (out_type=='Baum') THEN
    ncols=1
    ColName(1)='density'
  END IF
  READ(OutputUnit,*) ! Comment
END SUBROUTINE ReadHeaderTbl

SUBROUTINE WriteHeaderHaus(in_type) !dh
  CHARACTER(10) :: in_type

  WRITE(*,*)
  WRITE(*,*) 'Haus+',in_type
  WRITE(*,*)
  WRITE(*,*) 'ValuesCol',ValuesCol(1),'!position in input file'
  WRITE(*,*) '       >0  -->  Base,Height: read from input file (for Shape: Base=0)'
  WRITE(*,*) '       =0  -->  Base,Height: taken from standard values in Main program'
  WRITE(*,*) '                             (can be changed there)'
  WRITE(*,*)
  WRITE(*,*) 'Values    OUTPUT     <--    INPUT'
  IF (MINVAL(ValuesCol)<0.OR.MAXVAL(ValuesCol)>ncols) THEN
    WRITE(*,*) '***** Stop: ValuesCol wrong *****'
    STOP
  ELSE IF (MAXVAL(ValuesCol)==0) THEN
    WRITE(*,*) 1,'Base  ',BaseIn
    WRITE(*,*) 2,'Height',HeightIn
  ELSE IF (in_type=='Shape') THEN
    WRITE(*,*) 1,'Base  ',0.d0
    WRITE(*,*) 2,'Height',ValuesCol(1),ColName(ValuesCol(1))
  ELSE IF (in_type=='Table') THEN
    WRITE(*,*) 1,'Base  ',1,'. column in table'
    WRITE(*,*) 2,'Height',2,'. column in table'
  END IF
END SUBROUTINE WriteHeaderHaus

SUBROUTINE WriteHeaderStreet(in_type) !dh
  CHARACTER(10) :: in_type
  INTEGER :: i

  WRITE(*,*)
  IF (in_type=='Table') THEN
    WRITE(*,*) 'Street+Table'
    WRITE(*,*) '           -->  Base,Height,Width: read from input file'
  ELSE IF (in_type=='Shape') THEN
    WRITE(*,*) 'Street+Shape'
    WRITE(*,*) '           -->  Base,Height,Width: taken from standard values in Main program'
    WRITE(*,*) '                                   (can be changed there)'
  END IF
  WRITE(*,*)
  WRITE(*,*) 'ValuesNum ',ValuesNum,'!number of emissions'
  WRITE(*,*) 'ValuesName','           ',(ValuesName(i),' ',i=1,ValuesNum),'!target species'
  WRITE(*,*) 'ValuesCol ',ValuesCol,'!position in input file'
  WRITE(*,*) '       >0  -->  Emission: read from corresponding column in input file'
  WRITE(*,*) '       =0  -->  Emission: taken from corresponding ValuesIn in Main program'
  WRITE(*,*) '                          (can be changed there)'
  WRITE(*,*)
  IF (MAXVAL(ValuesCol)>0) THEN
    WRITE(*,*) 'Conversion factor for emissions (g/ma --> µg/ms):',ValuesScale
    WRITE(*,*) '(not applied to ValuesIn)'
  END IF
  WRITE(*,*)
  WRITE(*,*) 'Values    OUTPUT     <--    INPUT'
  IF (MINVAL(ValuesCol)<0.OR.MAXVAL(ValuesCol)>ncols) THEN
    WRITE(*,*) '***** Stop: ValuesCol wrong *****'
    STOP
  ELSE IF (in_type=='Shape') THEN
    WRITE(*,*) 1,'Base  ',BaseIn
    WRITE(*,*) 2,'Height',HeightIn
    WRITE(*,*) 3,'Width ',WidthIn
  ELSE IF (in_type=='Table') THEN
    WRITE(*,*) 1,'Base  ',1,'. column in table'
    WRITE(*,*) 2,'Height',2,'. column in table'
    WRITE(*,*) 3,'Width ',3,'. column in table'
  END IF
  WRITE(*,*)
  DO i=1,ValuesNum
    IF (ValuesCol(i)==0) THEN
      WRITE(*,*) i,ValuesName(i),ValuesIn(i)
    ELSE
      WRITE(*,*) i,ValuesName(i),' ',ValuesCol(i),ColName(ValuesCol(i))
    END IF
  END DO
END SUBROUTINE WriteHeaderStreet

SUBROUTINE WriteHeaderBaum(in_type) !dh
  CHARACTER(10) :: in_type
  INTEGER :: i

  WRITE(*,*)
  WRITE(*,*) 'Baum+',in_type
  WRITE(*,*)
  WRITE(*,*) 'ValuesCol',ValuesCol(1),'!position in input file'
  IF (in_type=='Table') THEN
    WRITE(*,*) '       >0  -->  Base,Height,Width: read from input file'
    WRITE(*,*) '                Density:           read from input file'
    WRITE(*,*) '       =0  -->  Base,Height,Width: read from input file'
    WRITE(*,*) '                Density:           taken from standard values in Main program'
    WRITE(*,*) '                                   (can be changed there)'
  ELSE IF (in_type=='Shape') THEN
    WRITE(*,*) '       =0  -->  Base,Height,Width: taken from standard values in Main program'
    WRITE(*,*) '                                   (can be changed there)'
    WRITE(*,*) '                Density:           taken from standard values in Main program'
    WRITE(*,*) '                                   (can be changed there)'
    WRITE(*,*) '       >0  -->  not realised'
  ELSE
    WRITE(*,*) '***** Stop: ValuesCol wrong *****'
    STOP
  END IF
  WRITE(*,*)
  WRITE(*,*) 'Values    OUTPUT      <--    INPUT'
  IF (MINVAL(ValuesCol)<0.OR.MAXVAL(ValuesCol)>ncols) THEN
    WRITE(*,*) '***** Stop: ValuesCol wrong *****'
    STOP
  ELSE IF (in_type=='Shape') THEN
    WRITE(*,*) 1,'Base   ',BaseIn
    WRITE(*,*) 2,'Height ',HeightIn
    WRITE(*,*) 3,'Width  ',WidthIn
  ELSE IF (in_type=='Table') THEN
    WRITE(*,*) 1,'Base   ',1,'. column in table'
    WRITE(*,*) 2,'Height ',2,'. column in table'
    WRITE(*,*) 3,'Width  ',3,'. column in table'
  END IF
  WRITE(*,*)
  DO i=1,ValuesNum
    IF (ValuesCol(i)==0) THEN
      WRITE(*,*) i,'Density',ValuesIn(i)
    ELSE
      WRITE(*,*) i,'Density',1,'. column in table (next row)'
    END IF
  END DO
END SUBROUTINE WriteHeaderBaum

SUBROUTINE ReadTable(Polygon,out_type,RecordEnd,RecordErr)

! Haus-Points anti-clockwise eingeben (sonst in CheckPolygon reorientation vornehmen!)
! Input-Ende: bei Datei-Ende oder Polygon%NumPoints=0
! Auslassung eines Objekts: bei Minuszeichen vor Polygon%NumPoints
!   (--> Polygon%Number erhaelt Minuszeichen)
  TYPE(Polygon_T) :: Polygon
  CHARACTER(10) :: out_type
  LOGICAL :: RecordEnd,RecordErr
  INTEGER :: i

  RecordEnd=.FALSE.
  RecordErr=.FALSE.
  NumberPolygon=NumberPolygon+1
  Polygon%Number=NumberPolygon
  Polygon%NumParts=1
  ALLOCATE(Polygon%Parts(0:0))
  Polygon%Parts(0)=0
  ALLOCATE(Polygon%Values(ValuesNum))
  READ(OutputUnit,*,END=1,ERR=1) Polygon%NumPoints ! Einlesen Input-Punktzahl
  IF      (Polygon%NumPoints==0) THEN ! Objekt markiert Ende
    GOTO 1
  ELSE IF (Polygon%NumPoints <0) THEN ! Objekt wird ueberlesen
    Polygon%NumPoints=-Polygon%NumPoints
    Polygon%Number=-NumberPolygon
  END IF
  IF     (out_type=='Haus') THEN
!   Einlesen von Haus: Points (unclosed ring), Base, Height
    Polygon%NumPoints=Polygon%NumPoints+1 ! Last point (no input) = first point
    ALLOCATE(Polygon%Points(0:Polygon%NumPoints-1))
    DO i=0,Polygon%NumPoints-2
      READ(OutputUnit,*,END=2,ERR=2) Polygon%Points(i)%x,Polygon%Points(i)%y
    END DO
    Polygon%Points(Polygon%NumPoints-1)%x=Polygon%Points(0)%x ! Closing of ring
    Polygon%Points(Polygon%NumPoints-1)%y=Polygon%Points(0)%y
    READ(OutputUnit,*,END=2,ERR=2) Polygon%Base,Polygon%Height
    IF (ValuesCol(1)==0) THEN
      Polygon%Base  =BaseIn
      Polygon%Height=HeightIn
    END IF
  ELSEIF (out_type=='Street') THEN
!   Einlesen von Street: Point series, Base, Height, Width, Emiss
    ALLOCATE(Polygon%Points(0:Polygon%NumPoints-1))
    DO i=0,Polygon%NumPoints-1
      READ(OutputUnit,*,END=2,ERR=2) Polygon%Points(i)%x,Polygon%Points(i)%y
    END DO
    READ(OutputUnit,*,END=2,ERR=2) Polygon%Base,Polygon%Height,Polygon%Width
    READ(OutputUnit,*,END=2,ERR=2) (ColValue(i),i=1,ncols)
    DO i=1,ValuesNum
      IF (ValuesCol(i)>0) THEN
        Polygon%Values(i)=ColValue(ValuesCol(i))*ValuesScale
      ELSE
        Polygon%Values(i)=ValuesIn(i)
      END IF
    END DO
  ELSEIF (out_type=='Baum') THEN
!   Einlesen von Baum: 2 Points, Base, Height, Width, Density
    ALLOCATE(Polygon%Points(0:Polygon%NumPoints-1))
    DO i=0,Polygon%NumPoints-1
      READ(OutputUnit,*,END=2,ERR=2) Polygon%Points(i)%x,Polygon%Points(i)%y
    END DO
    READ(OutputUnit,*,END=2,ERR=2) Polygon%Base,Polygon%Height,Polygon%Width
    IF (ValuesCol(1)==0) THEN
      Polygon%Values(1)=ValuesIn(1)
    ELSE
      READ(OutputUnit,*,END=2,ERR=2) Polygon%Values(1)
    END IF
  END IF
  RETURN
1 RecordEnd=.TRUE.
  RETURN
2 RecordEnd=.TRUE.
  RecordErr=.TRUE.
  RETURN

END SUBROUTINE ReadTable

SUBROUTINE PartitionPolygon(P)

  TYPE(Polygon_T), POINTER :: P

  TYPE(Polygon_T), POINTER :: NextPolygon,NewPolygon
  TYPE(Polygon_T), POINTER :: OldPolygon

  TYPE(Point2_T) :: Points(0:P%NumPoints-1)
  INTEGER :: Parts(0:P%NumParts)
  INTEGER :: NumParts,i,j,Number

  NumParts=P%NumParts
  HeightIn=P%Height
  BaseIn  =P%Base
  WidthIn =P%Width
  ValuesIn=P%Values

  IF (NumParts>1) THEN
    WRITE(*,*) 'Number Poly',P%Number,'NumParts',NumParts
    Number=P%Number
    DO i=0,P%NumPoints-1
      Points(i)=P%Points(i)
    END DO
    DO i=0,NumParts-1
      Parts(i)=P%Parts(i)
    END DO
    Parts(NumParts)=P%NumPoints
    DEALLOCATE(P%Points)
    DEALLOCATE(P%Parts)
    P%NumParts=1
    P%NumPoints=Parts(1)-Parts(0)
    ALLOCATE(P%Parts(0:0))
    P%Parts(0)=0
    ALLOCATE(P%Points(0:P%NumPoints-1))
    DO i=Parts(0),Parts(1)-1
      P%Points(i-Parts(0))=Points(i)
    END DO
    OldPolygon=>P
    DO j=1,NumParts-1
      ALLOCATE(NewPolygon)
      NextPolygon=>OldPolygon%Next
      OldPolygon%Next=>NewPolygon
      NewPolygon%Next=>NextPolygon
      NewPolygon%Number=Number
      NewPolygon%Height=HeightIn
      NewPolygon%Base=BaseIn
      NewPolygon%Width=WidthIn
      ALLOCATE(NewPolygon%Values(ValuesNum))
      NewPolygon%Values=ValuesIn
      NewPolygon%NumParts=j+1 ! Serial number (instead of total number before)
      NewPolygon%NumPoints=Parts(j+1)-Parts(j)
      ALLOCATE(NewPolygon%Parts(0:0))
      NewPolygon%Parts(0)=0
      ALLOCATE(NewPolygon%Points(0:NewPolygon%NumPoints-1))
      DO i=Parts(j),Parts(j+1)-1
        NewPolygon%Points(i-Parts(j))=Points(i)
      END DO
      OldPolygon=>NewPolygon
    END DO
    P=>NewPolygon
  END IF

END SUBROUTINE PartitionPolygon

SUBROUTINE Check(P,Pshift,Scale,out_type,in_type)

  TYPE(Polygon_T), TARGET :: P

  INTEGER :: i
  INTEGER :: NumPoints,NumPointsM
  CHARACTER(10) :: out_type,in_type
  TYPE(Point2_T) :: Pshift,Scale
  TYPE(Point2_T), ALLOCATABLE :: Points1(:),Points2(:)

! Minimal number of NumPoints for HAUS,STREET,BAUM
  NumPoints =P%NumPoints
  IF      (out_type=='Haus') THEN
    NumPointsM=4
  ELSE IF (out_type=='Street'.OR.out_type=='Baum') THEN
    NumPointsM=2
  END IF
! Objects with too few points are ignored
  IF (NumPoints<NumPointsM.OR.(NumPoints>NumPointsM.AND.out_type=='Baum')) THEN
    P%Number=-ABS(P%Number)
    RETURN
  END IF
! Objects with wrong parameters are ignored or corrected
  IF (P%Height<=P%Base) THEN
    P%Number=-ABS(P%Number)
  ELSE IF (       P%Width  <=0.d0.AND.(out_type=='Street'.OR.out_type=='Baum')) THEN
    P%Number=-ABS(P%Number)
  ELSE IF (MAXVAL(P%Values)<=0.d0.AND.(out_type=='Street'.OR.out_type=='Baum')) THEN
    P%Number=-ABS(P%Number)
  ELSE IF (MINVAL(P%Values)< 0.d0.AND.(out_type=='Street'.OR.out_type=='Baum')) THEN
    P%Number=-ABS(P%Number)
  END IF
  ALLOCATE(Points1(0:NumPoints-1),Points2(0:NumPoints-1))
! Shift + Scaling
  DO i=0,NumPoints-1
    Points2(i)%x=P%Points(i)%x*Scale%x
    Points2(i)%y=P%Points(i)%y*Scale%y
    Points2(i)%x=Points2(i)%x-Pshift%x
    Points2(i)%y=Points2(i)%y-Pshift%y
  END DO
! Reorientation clockwise (in Shape) --> anti-clockwise
  IF (in_type=='Shape') THEN
    DO i=0,NumPoints-1
      Points1(i)%x=Points2(NumPoints-1-i)%x
      Points1(i)%y=Points2(NumPoints-1-i)%y
    END DO
  ELSE
    DO i=0,NumPoints-1
      Points1(i)%x=Points2(i)%x
      Points1(i)%y=Points2(i)%y
    END DO
  END IF
! Correction for start=end point
  IF (out_type=='Haus') THEN
    IF (Points1(0)%x /= Points1(NumPoints-1)%x .OR. &
        Points1(0)%y /= Points1(NumPoints-1)%y ) THEN
      Points1(NumPoints-1)%x=Points1(0)%x
      Points1(NumPoints-1)%y=Points1(0)%y
    END IF
  END IF
! Restore
  DO i=0,P%NumPoints-1
    P%Points(i)%x=Points1(i)%x
    P%Points(i)%y=Points1(i)%y
  END DO
  DEALLOCATE(Points1,Points2)

END SUBROUTINE Check

RECURSIVE SUBROUTINE DecomposePolygon(P,Number,Tol)

  TYPE(Polygon_T), TARGET :: P

  INTEGER :: i,j,k,i1,j1,k1
  TYPE(Polygon_T), POINTER :: NextPolygon,NewPolygon
  TYPE(Point_T) :: P1,P2,P3,PCut
  TYPE(Point_T) :: P1P2,P2P3,P2PCut
  TYPE(Point_T) :: N
  INTEGER :: NumPoints,NumPoints1,NumPoints2,Number,NumParts
  TYPE(Point2_T), ALLOCATABLE :: Points1(:),Points2(:)
  REAL(8) :: ss,Tol

! Recursive decomposition (2 parts) in the case of concave corners
! --> final objects = convex only
! Also check for total convexity

  NumParts=P%NumParts
  NumPoints=P%NumPoints
  HeightIn=P%Height
  BaseIn  =P%Base
  WidthIn =P%Width
  ValuesIn=P%Values
  NextPolygon=>P%Next
  IF (Number<0) THEN
    RETURN
  END IF
  DO i=1,NumPoints-1
    P1%x=P%Points(i-1)%x
    P1%y=P%Points(i-1)%y
    P2%x=P%Points(i)%x
    P2%y=P%Points(i)%y
    i1=i+1
    IF (i1>NumPoints-1) THEN
      i1=i1-(NumPoints-1)
    END IF
    P3%x=P%Points(i1)%x
    P3%y=P%Points(i1)%y
    P1P2=P2-P1
    P2P3=P3-P2
    N=(P1P2.CROSS.P2P3) /MAX(Norm(P1P2),Norm(P2P3))
    ss=N%z
    IF (ss<-Tol) THEN
      ! Concave corner found --> search for cut point
      DO j=i+2,i+2+NumPoints-5
        j1=j
        IF (j1>NumPoints-1) THEN
          j1=j1-(NumPoints-1)
        END IF
        PCut%x=P%Points(j1)%x 
        PCut%y=P%Points(j1)%y 
        P2PCut=PCut-P2
        N=(P1P2.CROSS.P2PCut) /MAX(Norm(P1P2),Norm(P2PCut))
        ss=N%z
        IF (ss>=-Tol) THEN
          ! Cut point found --> division into 2 parts
          ! Points1 (remaining part, concave corner --> made convex)
          NumPoints1=i+1+(NumPoints-j)
          ALLOCATE(Points1(0:NumPoints1-1))
          DO k=j,j+i+NumPoints-j-1
            k1=k
            IF (k1>NumPoints-1) THEN
              k1=k1-(NumPoints-1)
            END IF
            Points1(k-j)=P%Points(k1)
          END DO
          Points1(NumPoints1-1)=P%Points(j1)
          ! Points2 (cut)
          NumPoints2=j-i+2
          ALLOCATE(Points2(0:NumPoints2-1))
          DO k=i,j
            k1=k
            IF (k1>NumPoints-1) THEN
              k1=k1-(NumPoints-1)
            END IF
            Points2(k-i)=P%Points(k1)
          END DO
          Points2(NumPoints2-1)=P%Points(i)
          ! Points1 --> P%Points
          DEALLOCATE(P%Points)
          ALLOCATE(P%Points(0:NumPoints1-1))
          P%NumPoints=NumPoints1
          P%Points=Points1
          !P%Number=Number
          P%NumParts=NumParts
          DEALLOCATE(Points1)
          ! Points2 --> NewPolygon%Points
          ALLOCATE(NewPolygon)
          ALLOCATE(NewPolygon%Points(0:NumPoints2-1))
          NewPolygon%NumPoints=NumPoints2
          NewPolygon%Points=Points2
          !NewPolygon%Number=Number
          NewPolygon%NumParts=NumParts
          ALLOCATE(NewPolygon%Parts(0:0))
          NewPolygon%Parts(0)=0
          NewPolygon%Height=HeightIn
          NewPolygon%Base=BaseIn
          NewPolygon%Width=WidthIn
          ALLOCATE(NewPolygon%Values(ValuesNum))
          NewPolygon%Values=ValuesIn
          DEALLOCATE(Points2)
          ! Rescan of both parts
          P%Next=>NewPolygon
          NewPolygon%Next=>NextPolygon
          CALL DecomposePolygon(P,Number,Tol)
          CALL DecomposePolygon(NewPolygon,Number,Tol)
          RETURN
        END IF
      END DO
      !WRITE(*,*) 'Attention! (Partial) object = concave, entire object ignored:' &
      !           ,Number,NumPoints
      Number=-ABS(Number)
      RETURN
    END IF
  END DO

END SUBROUTINE DecomposePolygon

SUBROUTINE SimplifyPolygon(P,Tol)

  TYPE(Polygon_T), TARGET :: P

  INTEGER :: i,j,i1,kv,kx
  TYPE(Point_T) :: P1,P2,P3
  TYPE(Point_T) :: P1P2,P2P3
  TYPE(Point_T) :: N
  INTEGER :: NumPoints,NumPoints1
  TYPE(Point2_T), ALLOCATABLE :: Points1(:),Points2(:)
  REAL(8) :: ss, Tol

  NumPoints =P%NumPoints
  NumPoints1=NumPoints
  ALLOCATE(Points1(0:NumPoints-1),Points2(0:NumPoints-1))
! P%Points --> Points1
  DO i=0,NumPoints-1
    Points1(i)%x=P%Points(i)%x
    Points1(i)%y=P%Points(i)%y
    Points2(i)%x=1.d0 ! default
  END DO
! Signing dispensable (lineal) points by Points2%x, numbering of convex corners
  P1%x=Points1(0)%x
  P1%y=Points1(0)%y
  kv=0 ! concave corners
  kx=0 ! convex  corners
  DO i=1,NumPoints-1
    P2%x=Points1(i)%x
    P2%y=Points1(i)%y
    i1=i+1
    IF (i1>NumPoints-1) THEN
      DO j=0,NumPoints-4
        i1=1+j
        IF (Points2(i1)%x>0.d0) THEN
          EXIT
        ELSE IF (j==NumPoints-4) THEN
          !WRITE(*,*) 'Attention! Too many dispensable points, object ignored:' &
          !           ,P%Number,NumPoints
          P%Number=-ABS(P%Number)
          GOTO 1
        END IF
      END DO
    END IF
    P3%x=Points1(i1)%x
    P3%y=Points1(i1)%y
    P1P2=P2-P1
    P2P3=P3-P2
    N=(P1P2.CROSS.P2P3) /MAX(Norm(P1P2),Norm(P2P3))
    ss=N%z ! = Lmin * Sinus = length deviation from lineal, /= vector product
    IF (ABS(ss)<Tol) THEN
      NumPoints1=NumPoints1-1
      Points2(i)%x=0.d0 ! dispensable
      IF (NumPoints1<4) THEN
        !WRITE(*,*) 'Attention! NumPoints < 4, simplified object ignored:' &
        !           ,P%Number,NumPoints
        P%Number=-ABS(P%Number)
        GOTO 1
      END IF
    ELSE
      P1%x=P2%x
      P1%y=P2%y
      Points2(i)%x=1.d0 ! indispensable
      IF (ss<0) THEN
        kv=kv+1
      ELSE
        kx=kx+1
      END IF
    END IF
  END DO
  IF (kx<3) THEN
!   Rejection of erroneous objects
    !WRITE(*,*) 'Attention! Number of convex corners < 3, object ignored:' &
    !           ,P%Number,NumPoints
    P%Number=-ABS(P%Number)
    GOTO 1
  ELSE IF (NumPoints1<NumPoints) THEN
! New dimensioning (shortening)
    i1=0
    DO i=1,NumPoints-1
      IF (Points2(i)%x>0.d0) THEN
        i1=i1+1
        Points1(i1)%x=Points1(i)%x
        Points1(i1)%y=Points1(i)%y
      ELSE IF (i==NumPoints-1) THEN
        Points1(0)%x=Points1(NumPoints1-1)%x
        Points1(0)%y=Points1(NumPoints1-1)%y
      END IF
    END DO
    DEALLOCATE(P%Points)
    ALLOCATE(P%Points(0:NumPoints1-1))
    P%NumPoints=NumPoints1
  END IF
! Restore
  DO i=0,P%NumPoints-1
    P%Points(i)%x=Points1(i)%x
    P%Points(i)%y=Points1(i)%y
  END DO

1 DEALLOCATE(Points1,Points2)

END SUBROUTINE SimplifyPolygon

END MODULE Shape_Mod
