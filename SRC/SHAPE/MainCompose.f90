PROGRAM MainCompose

  USE OSD_Mod
  USE DBase_Mod
  USE Shape_Mod

  IMPLICIT NONE

  CHARACTER(80) :: FileName=''
  CHARACTER(80) :: FileHaus='Neustadt.Haus'
  INTEGER :: OutputUnitHaus=42
  CHARACTER :: out_type
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

  TYPE(Polygon_T), POINTER :: NewPoLi
  TYPE(Polygon_T), POINTER :: CurrPoLi
  TYPE(Polygon_T), POINTER :: StartPoLi
  TYPE(Polygon_T), POINTER :: Polygon
! TYPE(RecordDbf_T) :: Record
  REAL(RealDouble) :: Height

 
  ALLOCATE(Polygon)
  Polygon%NumPoints=7
  ALLOCATE(Polygon%Points(0:Polygon%NumPoints-1))
  Polygon%Points(0)%x=0.0d0
  Polygon%Points(0)%y=2.0d0

  Polygon%Points(1)%x=0.0d0
  Polygon%Points(1)%y=0.0d0

  Polygon%Points(2)%x=1.0d0
  Polygon%Points(2)%y=0.0d0

  Polygon%Points(3)%x=1.0d0
  Polygon%Points(3)%y=-1.0d0

  Polygon%Points(4)%x=2.0d0
  Polygon%Points(4)%y=-1.0d0

  Polygon%Points(5)%x=2.0d0
  Polygon%Points(5)%y=2.0d0

  Polygon%Points(6)%x=0.0d0
  Polygon%Points(6)%y=2.0d0

  CALL ComposePolygon(Polygon)

  OPEN(UNIT=OutputUnitHaus,FILE=FileHaus,STATUS='UNKNOWN')
  Height=20.0d0
  DO 
    CALL WriteHaus(Polygon,Height)
    IF (ASSOCIATED(Polygon%Next)) THEN
      Polygon=>Polygon%Next
    ELSE
      EXIT
    END IF 
  END DO

CONTAINS 

SUBROUTINE WriteHaus(PolygonSh,Height)
  TYPE(Polygon_T), POINTER :: PolygonSh
  REAL(8) :: Height
  INTEGER :: i
  INTEGER :: NumEdges(0:PolygonSh%NumPoints)
  CHARACTER*2 :: TypeFacecs(0:PolygonSh%NumPoints)

! Vom Punkt erst alle x, dann y, dann z

  DO i=0,PolygonSh%NumPoints-2
    NumEdges(i)=4
    TypeFacecs(i)='w '
  END DO
  NumEdges(PolygonSh%NumPoints-1)=PolygonSh%NumPoints
  TypeFacecs(PolygonSh%NumPoints-1)='r '
  WRITE(OutputUnitHaus,*) PolygonSh%NumPoints
  WRITE(OutputUnitHaus,*) (NumEdges(i),i=0,PolygonSh%NumPoints-1)
  WRITE(OutputUnitHaus,*) (TypeFacecs(i),i=0,PolygonSh%NumPoints-1)
  DO i=0,PolygonSh%NumPoints-2
    WRITE(OutputUnitHaus,*) PolygonSh%Points(i)%x   &
                     ,PolygonSh%Points(i+1)%x &
                     ,PolygonSh%Points(i+1)%x &
                     ,PolygonSh%Points(i)%x   &
                     ,PolygonSh%Points(i)%y   &
                     ,PolygonSh%Points(i+1)%y &
                     ,PolygonSh%Points(i+1)%y &
                     ,PolygonSh%Points(i)%y   &
                     ,0.0d0                   &
                     ,0.0d0                   &
                     ,Height                  &
                     ,Height
  END DO
  WRITE(OutputUnitHaus,*) (PolygonSh%Points(i)%x,i=0,PolygonSh%NumPoints-1)   &
                         ,(PolygonSh%Points(i)%y,i=0,PolygonSh%NumPoints-1)   &
                         ,(Height,i=0,PolygonSh%NumPoints-1)   
END SUBROUTINE WriteHaus


SUBROUTINE WriteHaus1(PolygonSh,Height)
  TYPE(Polygon_T), POINTER :: PolygonSh
  REAL(8) :: Height
  INTEGER :: i
  INTEGER :: NumEdges(0:PolygonSh%NumPoints)
  CHARACTER*2 :: TypeFacecs(0:PolygonSh%NumPoints)

! Vom Punkt erst alle x, dann y, dann z

  DO i=0,PolygonSh%NumPoints-2
    NumEdges(i)=4
    TypeFacecs(i)='w '
  END DO
  NumEdges(PolygonSh%NumPoints-1)=PolygonSh%NumPoints
  TypeFacecs(PolygonSh%NumPoints-1)='r '
  WRITE(OutputUnitHaus,*) PolygonSh%NumPoints
  WRITE(OutputUnitHaus,*) (NumEdges(i),i=0,PolygonSh%NumPoints-1)
  WRITE(OutputUnitHaus,*) (TypeFacecs(i),i=0,PolygonSh%NumPoints-1)
  DO i=0,PolygonSh%NumPoints-2
    WRITE(OutputUnitHaus,*) PolygonSh%Points(i)%x   &
                     ,PolygonSh%Points(i+1)%x &
                     ,PolygonSh%Points(i+1)%x &
                     ,PolygonSh%Points(i)%x   &
                     ,PolygonSh%Points(i)%y   &
                     ,PolygonSh%Points(i+1)%y &
                     ,PolygonSh%Points(i+1)%y &
                     ,PolygonSh%Points(i)%y   &
                     ,0.0d0                   &
                     ,0.0d0                   &
                     ,Height                  &
                     ,Height
  END DO
  WRITE(OutputUnitHaus,*) (PolygonSh%Points(i)%x,i=0,PolygonSh%NumPoints-1)   &
                         ,(PolygonSh%Points(i)%y,i=0,PolygonSh%NumPoints-1)   &
                         ,(Height,i=0,PolygonSh%NumPoints-1)   
END SUBROUTINE WriteHaus1

END PROGRAM MainCompose
