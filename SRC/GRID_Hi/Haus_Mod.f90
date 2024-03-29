MODULE Haus_Mod

  USE Domain_Mod 
  IMPLICIT NONE

  TYPE GFace_T
    INTEGER :: NumberOfPoints
    TYPE(Point_T), POINTER :: Points(:)
    INTEGER, POINTER :: ListOfPoints(:)
    TYPE(Point_T) :: Normal 
    CHARACTER :: Type         ! r=roof,w=wall
  END TYPE GFace_T

  TYPE Haus_T
    INTEGER :: NumberOfPoints
    INTEGER :: NumberOfFaces
    TYPE(Point_T), POINTER :: Points(:)
    TYPE(GFace_T), POINTER :: Faces(:)
    Type(Point_T) :: P0,P1
  END TYPE Haus_T
 

CONTAINS 

FUNCTION Dist(P,Haus)
  REAL(8) :: Dist
  TYPE(Point_T) :: P
  TYPE(Haus_T) :: Haus

  INTEGER :: i
  REAL(8) :: DistFace
  TYPE(Point_T) :: P1
  TYPE(GFace_T), POINTER :: Face

  IF (Haus%P0<=P.AND.P<=Haus%P1) THEN
    Dist=-1.0d99
    DO i=1,Haus%NumberOfFaces
      Face=>Haus%Faces(i)
      P1=Face%Points(1)
      DistFace=(P-P1)*Face%Normal
      Dist=MAX(Dist,DistFace)
    END DO
  ELSE
    Dist=1.0d0
  END IF 
END FUNCTION Dist

SUBROUTINE BoundingBox(Haus)
  TYPE(Haus_T) :: Haus

  INTEGER :: i,j
  TYPE(GFace_T), POINTER :: Face
  Haus%P0%x=1.d99
  Haus%P0%y=1.d99
  Haus%P0%z=1.d99
  Haus%P1%x=-1.d99
  Haus%P1%y=-1.d99
  Haus%P1%z=-1.d99
  DO i=1,Haus%NumberOfFaces
    Face=>Haus%Faces(i)
    DO j=1,Face%NumberOfPoints
      Haus%P0%x=MIN(Haus%P0%x,Face%Points(j)%x)
      Haus%P0%y=MIN(Haus%P0%y,Face%Points(j)%y)
      Haus%P0%z=MIN(Haus%P0%z,Face%Points(j)%z)
      Haus%P1%x=MAX(Haus%P1%x,Face%Points(j)%x)
      Haus%P1%y=MAX(Haus%P1%y,Face%Points(j)%y)
      Haus%P1%z=MAX(Haus%P1%z,Face%Points(j)%z)
    END DO
  END DO

  IF (Haus%P1%x<=Domain%x0.OR.Haus%P0%x>=Domain%x1.OR. &
      Haus%P1%y<=Domain%y0.OR.Haus%P0%y>=Domain%y1.OR. &
      Haus%P1%z<=Domain%z0.OR.Haus%P0%z>=Domain%z1) THEN
    Haus%Faces(Haus%NumberOfFaces)%Normal%z=-2.d0 ! sign for outside
    RETURN
  END IF

END SUBROUTINE BoundingBox

SUBROUTINE NormalForm(Haus)
  TYPE(Haus_T) :: Haus

  INTEGER :: i,i1,j,NumFaces,NumPoints
  TYPE(Point_T) :: P1,P2,P3
  TYPE(Point_T) :: P3mP1,P3mP2
  TYPE(GFace_T), POINTER :: Face
  REAL(8) :: Temp,Temp1

! All face normals
  NumFaces=Haus%NumberOfFaces
  DO i=1,Haus%NumberOfFaces
    Face=>Haus%Faces(i)
    P1=Face%Points(1)
    P2=Face%Points(2)
    P3=Face%Points(3)
    P3mP1=P3-P1
    P3mP2=P3-P2
    Face%Normal=P3mP1.CROSS.P3mP2
    Temp=Norm(Face%Normal)
    IF (Temp>1.d-4) THEN
      Face%Normal=Face%Normal/Temp
    ELSE
      Haus%Faces(NumFaces)%Normal%z=-1.d0 ! sign for faulty
      RETURN
    END IF
  END DO
  IF (Haus%Faces(NumFaces)%Normal%z<=0.d0) THEN
    Haus%Faces(NumFaces)%Normal%z=-1.d0 ! sign for faulty
    RETURN
  END IF

! Check of convexity of side walls (for roof points)
! Cancel of erroneous objects by Haus%Faces(NumFaces=Roof)%Normal%z=-1.d0
! (see SUBROUTINEs Read_Val_Haus...)
! Precondition for this check:  Haus%Faces(i)%Points(3) = roof points
  NumPoints=Haus%Faces(NumFaces)%NumberOfPoints
  DO j=1,NumPoints ! All roof points j
    P1=Haus%Faces(NumFaces)%Points(j) ! Roof point j
    i1=0 ! number of points on actual side wall i (i1=2 after i-loop!)
    Temp=-1.0d99
    DO i=1,NumPoints ! All side walls i
      P2=Haus%Faces(i)%Points(3) ! Upper point of side wall i (roof = 3 and 4)
      Temp1=(P1-P2)*Haus%Faces(i)%Normal
      IF (Temp1>=-1.d-2) THEN
        i1=i1+1
      END IF
      IF (Temp1<-1.d-2.OR.i1>2) THEN
        Temp=MAX(Temp,Temp1)
      END IF
    END DO
    IF (i1/=2.OR.Temp>=-1.d-2) THEN
      Haus%Faces(NumFaces)%Normal%z=-1.d0 ! sign for faulty
      RETURN
    END IF
  END DO

END SUBROUTINE NormalForm

END MODULE Haus_Mod

