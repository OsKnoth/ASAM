SUBROUTINE ReadInteger(I)
  INTEGER :: I

  nRec=nRec+iRec
  READ(OutputUnit,REC=nRec) I
END SUBROUTINE ReadInteger

SUBROUTINE ReadIntegerB(I)
  INTEGER :: I

  INTEGER :: IB

  nRec=nRec+iRec
  READ(OutputUnit,REC=nRec) IB
  CALL MVBITS(IB,24,8,I, 0)
  CALL MVBITS(IB,16,8,I, 8)
  CALL MVBITS(IB, 8,8,I,16)
  CALL MVBITS(IB, 0,8,I,24)
END SUBROUTINE ReadIntegerB

SUBROUTINE ReadRealDouble(R)
  REAL(RealDouble) :: R

  INTEGER :: iReal(1:2)
  nRec=nRec+iRec
  READ(OutputUnit,REC=nRec) iReal(1)
  nRec=nRec+iRec
  READ(OutputUnit,REC=nRec) iReal(2)
  R=TRANSFER(iReal,0.0d0)
END SUBROUTINE ReadRealDouble

SUBROUTINE ReadPoint(P)

  TYPE(Point_T) :: P
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
  DO i=0,PolylineZ%NumParts-1
    CALL ReadInteger(PolylineZ%Parts(i))
  END DO
! Points
  DO i=0,PolylineZ%NumPoints-1
    CALL ReadPoint(PolylineZ%Points(i))
  END DO
! z Bounding Box
  CALL ReadRealDouble(PolylineZ%zMin)
  CALL ReadRealDouble(PolylineZ%zMax)
! z Points
  DO i=0,PolylineZ%NumPoints-1
    CALL ReadRealDouble(PolylineZ%zPoints(i))
  END DO
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



