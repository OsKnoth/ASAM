MODULE Haus_Mod

  USE Geometry_Mod 
  USE Control_Mod
  IMPLICIT NONE

  TYPE GFace_T
    INTEGER :: NumberOfPoints
    TYPE(Point_T), POINTER :: Points(:)
    INTEGER, POINTER :: ListOfPoints(:)
    TYPE(Point_T) :: Normal 
    CHARACTER :: Type         ! r=roof,w=wall
  END TYPE GFace_T

  TYPE Haus_T
    INTEGER :: Number
    INTEGER :: NumberOfPoints
    INTEGER :: NumberOfFaces
    TYPE(Point_T), POINTER :: Points(:)
    TYPE(GFace_T), POINTER :: Faces(:)
    Type(Point_T) :: P0,P1
    Type(Point_T) :: PM
  END TYPE Haus_T

  INTEGER :: NumberHaus
  TYPE(Haus_T),POINTER :: Haus(:)
 

CONTAINS 

FUNCTION Dist(P,Haus,DistSurf)
  REAL(8) :: Dist
  TYPE(Point_T) :: P
  TYPE(Haus_T) :: Haus
  REAL(8) :: DistSurf

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
  Dist=Dist-DistSurf
END FUNCTION Dist

SUBROUTINE BoundingBox(Haus)
  TYPE(Haus_T) :: Haus

  INTEGER :: i
  Haus%P0%x=1.d99
  Haus%P0%y=1.d99
  Haus%P0%z=1.d99
  Haus%P1%x=-1.d99
  Haus%P1%y=-1.d99
  Haus%P1%z=-1.d99
  DO i=1,Haus%NumberOfPoints
    Haus%P0%x=MIN(Haus%P0%x,Haus%Points(i)%x)
    Haus%P0%y=MIN(Haus%P0%y,Haus%Points(i)%y)
    Haus%P0%z=MIN(Haus%P0%z,Haus%Points(i)%z)
    Haus%P1%x=MAX(Haus%P1%x,Haus%Points(i)%x)
    Haus%P1%y=MAX(Haus%P1%y,Haus%Points(i)%y)
    Haus%P1%z=MAX(Haus%P1%z,Haus%Points(i)%z)
  END DO
  Haus%P0=Haus%P0-Point_T(1.0d-1,1.0d-1,1.0d-1)
  Haus%P1=Haus%P1+Point_T(1.0d-1,1.0d-1,1.0d-1)
  Haus%PM%x=0.0d0
  Haus%PM%y=0.0d0
  Haus%PM%z=0.0d0
  DO i=1,Haus%NumberOfPoints
    Haus%PM%x=Haus%PM%x+Haus%Points(i)%x
    Haus%PM%y=Haus%PM%y+Haus%Points(i)%y
    Haus%PM%z=Haus%PM%z+Haus%Points(i)%z
  END DO  
  Haus%PM%x=Haus%PM%x/Haus%NumberOfPoints
  Haus%PM%y=Haus%PM%y/Haus%NumberOfPoints
  Haus%PM%z=Haus%PM%z/Haus%NumberOfPoints
END SUBROUTINE BoundingBox

SUBROUTINE NormalForm(Haus)
  TYPE(Haus_T) :: Haus

  INTEGER :: i,j,jMax
  INTEGER :: NumberPoints
  TYPE(Point_T) :: P1,P2,P3
  TYPE(Point_T) :: P3mP1,P3mP2
  TYPE(GFace_T), POINTER :: Face
  REAL(8) :: AngleMax,Temp

  DO i=1,Haus%NumberOfFaces
    Face=>Haus%Faces(i)
    jMax=0
    AngleMax=-1.0d0
    NumberPoints=SIZE(Face%Points)
    DO j=1,NumberPoints
      IF (j==1) THEN
        P1=Face%Points(NumberPoints)
        P2=Face%Points(j)
        P3=Face%Points(j+1)
      ELSE IF (j==NumberPoints) THEN
        P1=Face%Points(j-1)
        P2=Face%Points(j)
        P3=Face%Points(1)
      ELSE  
        P1=Face%Points(j-1)
        P2=Face%Points(j)
        P3=Face%Points(j+1)
      END IF  
      P3mP1=P3-P1
      P3mP2=P3-P2
      Temp=Norm(P3mP1.CROSS.P3mP2)/Norm(P3mP1)/Norm(P3mP2)
      IF (ABS(TEMP)>AngleMax) THEN
        jMax=j
        AngleMax=Temp
      END IF  
    END DO  
    IF (jMax==1) THEN
      P1=Face%Points(NumberPoints)
      P2=Face%Points(jMax)
      P3=Face%Points(jMax+1)
    ELSE IF (jMax==NumberPoints) THEN
      P1=Face%Points(jMax-1)
      P2=Face%Points(jMax)
      P3=Face%Points(1)
    ELSE  
      P1=Face%Points(jMax-1)
      P2=Face%Points(jMax)
      P3=Face%Points(jMax+1)
    END IF  
    P3mP1=P3-P1
    P3mP2=P3-P2
    Face%Normal=P3mP1.CROSS.P3mP2
    Temp=Norm(Face%Normal)
    Face%Normal=Face%Normal/Temp
  END DO
END SUBROUTINE NormalForm

SUBROUTINE WriteHausVTK(FileName)
  CHARACTER(*) :: FileName
  INTEGER ::i,h
  INTEGER :: NumberOfPoints
  INTEGER :: iFaceOffsets=0
  INTEGER :: PointList(100)

  DO i=1,100
    PointList(i)=i-1
  END DO  

  NumberOfPoints=0
  DO h=1,SIZE(Haus)
    NumberOfPoints=NumberOfPoints+Haus(h)%NumberOfPoints
  END DO  


  OPEN(UNIT=10,FILE=TRIM(FileName)//'_'//'Haus'//'.vtu',STATUS='UNKNOWN')
  WRITE(10,'(a21)') '<?xml version="1.0"?>'
  WRITE(10,'(a33)') '<VTKFile type="UnstructuredGrid">'
  WRITE(10,'(a18)') '<UnstructuredGrid>'
  WRITE(10,'(a23,i8,a17,i8,a3)') '<Piece NumberOfPoints="',NumberOfPoints,'" NumberOfCells="',SIZE(Haus),'">'
  WRITE(10,'(a46)') '<CellData Scalars="scalars" Vector="Velocity">'
  WRITE(10,*) '<DataArray type="Float64" Name="Number" Format="ascii">'
  DO h=1,SIZE(Haus)
    WRITE(10,*) Haus(h)%Number
  END DO  
  WRITE(10,'(a12)') '</DataArray>'
  WRITE(10,'(a11)') '</CellData>'
  WRITE(10,'(a8)') '<Points>'
  WRITE(10,'(a64)') '<DataArray type="Float32" NumberOfComponents="3" format="ascii">'
  DO h=1,SIZE(Haus)
    DO i=1,Haus(h)%NumberOfPoints
      WRITE(10,*) Haus(h)%Points(i)%x-xOffset,Haus(h)%Points(i)%y-yOffset,Haus(h)%Points(i)%z
    END DO  
  END DO  
  WRITE(10,'(a12)') '</DataArray>'
  WRITE(10,'(a9)') '</Points>'
  WRITE(10,'(a7)') '<Cells>'
  WRITE(10,'(a59)') '<DataArray type="Int32" Name="connectivity" format="ascii">'
  NumberOfPoints=0
  DO h=1,SIZE(Haus)
    WRITE(10,*) PointList(1:Haus(h)%NumberOfPoints)+NumberOfPoints
    NumberOfPoints=NumberOfPoints+Haus(h)%NumberOfPoints
  END DO  
  WRITE(10,'(a12)') '</DataArray>'
  WRITE(10,'(a54)') '<DataArray type="Int32" Name="offsets" format="ascii">'
  NumberOfPoints=0
  DO h=1,SIZE(Haus)
    NumberOfPoints=NumberOfPoints+Haus(h)%NumberOfPoints
    WRITE(10,*) NumberOfPoints
  END DO  
  WRITE(10,'(a12)') '</DataArray>'
  WRITE(10,'(a52)') '<DataArray type="UInt8" Name="types" format="ascii">'
  DO h=1,SIZE(Haus)
    WRITE(10,*) 42
  END DO  
  WRITE(10,'(a12)') '</DataArray>'
  WRITE(10,'(a52)') '<DataArray type="Int32" Name="faces" format="ascii">'
  NumberOfPoints=-1
  DO h=1,SIZE(Haus)
    WRITE(10,*) Haus(h)%NumberOfFaces
    DO i=1,Haus(h)%NumberOfFaces
      WRITE(10,*) Haus(h)%Faces(i)%NumberOfPoints,Haus(h)%Faces(i)%ListOfPoints+NumberOfPoints
    END DO  
    NumberOfPoints=NumberOfPoints+Haus(h)%NumberOfPoints
  END DO  
  WRITE(10,'(a12)') '</DataArray>'
  WRITE(10,'(a58)') '<DataArray type="Int32" Name="faceoffsets" format="ascii">'
  iFaceOffsets=0
  DO h=1,SIZE(Haus)
    iFaceOffsets=iFaceOffsets+1
    DO i=1,Haus(h)%NumberOfFaces
      iFaceOffsets=iFaceOffsets+1+Haus(h)%Faces(i)%NumberOfPoints
    END DO  
    WRITE(10,*) iFaceOffsets
  END DO
  WRITE(10,'(a12)') '</DataArray>'
  WRITE(10,'(a8)') '</Cells>'
  WRITE(10,'(a8)') '</Piece>'
  WRITE(10,'(a19)') '</UnstructuredGrid>'
  WRITE(10,'(a10)') '</VTKFile>'
  CLOSE(10)
END SUBROUTINE WriteHausVTK

END MODULE Haus_Mod

