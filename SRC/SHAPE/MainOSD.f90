PROGRAM MainConvETYP

  USE OSD_Mod
  USE Shape_Mod

  IMPLICIT NONE

  CHARACTER(80) :: InputFileName
  CHARACTER(80) :: OutputFile
  CHARACTER(80) :: FileName='kiga'
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

  TYPE(PolyLineZ_T), POINTER :: NewPoLiZ
  TYPE(PolyLineZ_T), POINTER :: CurrPoLiZ
 
  CALL getarg(1,InputFileName)
  CALL getarg(2,FileName)
  OutputFile=FileName
  i=SCAN(OutputFile,'.')
  IF(i>0) THEN
    l=LEN(OutputFile)
    OutputFile(i:l)=" "
  END IF
  
  WRITE(*,*) InputFileName
  WRITE(*,*) "         ...Read_Line_Points..."
  CALL Read_Lpoints(InputFileName)
  WRITE(*,*) "     ... Ende  Read_Line_Points..."

  LengthAllObjekte=0 
  CurrBuild=>ListOfBuild
  ListOfPoLiZ=>NewPoliZ

  xMin=1.d99
  xMax=-1.d99
  yMin=1.d99
  yMax=-1.d99
  DO nr_oj=0,NumberOfBuild-1
    PolyLineZ%NumPoints=CurrBuild%NumberOfLines*2
    PolyLineZ%NumParts=1
    ALLOCATE(PolyLineZ%Parts(0:PolyLineZ%NumParts-1))
    ALLOCATE(PolyLineZ%Points(0:PolyLineZ%NumPoints-1))
    ALLOCATE(PolyLineZ%zPoints(0:PolyLineZ%NumPoints-1))
!   ALLOCATE(PolyLineZ%MPoints(0:PolyLineZ%NumPoints-1))
!   PolyLineZ%MPoints=0.0d0
    PolyLineZ%Parts(0)=0
    Current=>CurrBuild%ListOfLines 
    Last_StoreP=Current%P2
    p=0
    PolyLineZ%xMin=MIN(Current%P1%xP,Current%P2%xP)
    PolyLineZ%yMin=MIN(Current%P1%yP,Current%P2%yP)
    PolyLineZ%xMax=MAX(Current%P1%xP,Current%P2%xP)
    PolyLineZ%yMax=MAX(Current%P1%yP,Current%P2%yP)
    PolyLineZ%zMin=0.0d0
    PolyLineZ%zMax=1.0d0
!   xMin=PolyLineZ%xMin
!   yMin=PolyLineZ%yMin
!   xMax=PolyLineZ%xMax
!   yMax=PolyLineZ%yMax
    DO i=0,CurrBuild%NumberOfLines-1
      !Belege Struktur PolyLineZ
      h_xmin=MIN(Current%P1%xP,Current%P2%xP)
      h_ymin=MIN(Current%P1%yP,Current%P2%yP)
      !h_zmin=MIN()
      h_xmax=MAX(Current%P1%xP,Current%P2%xP)
      h_ymax=MAX(Current%P1%yP,Current%P2%yP)
      !h_zmax=MAX()
      PolyLineZ%xMin=MIN(PolyLineZ%xMin,h_xmin)
      PolyLineZ%yMin=MIN(PolyLineZ%yMin,h_ymin)
      PolyLineZ%zMin=MIN(PolyLineZ%zMin,h_zmin)
      PolyLineZ%xMax=MAX(PolyLineZ%xMax,h_xmax)
      PolyLineZ%yMax=MAX(PolyLineZ%yMax,h_ymax)
      PolyLineZ%zMax=MAX(PolyLineZ%zMax,h_zmax)
    
      PolyLineZ%Points(p)%x=Current%P1%xP
      PolyLineZ%Points(p)%y=Current%P1%yP
      PolyLineZ%zPoints(p)=0
     !LineZ%MPoints(p)=
      p=p+1
      PolyLineZ%Points(p)%x=Current%P2%xP
      PolyLineZ%Points(p)%y=Current%P2%yP
      PolyLineZ%zPoints(p)=0
     !LineZ%MPoints(p)=
      IF(i>0) THEN
        IF(Last_StoreP%xP==Current%P1%xP.AND.Last_StoreP%yP==Current%P1%yP) THEN
          Last_StoreP=Current%P2
        ELSE
          Pzw=PolyLineZ%Points(p-1)
          PolyLineZ%Points(p-1)=PolyLineZ%Points(p)
          PolyLineZ%Points(p)=Pzw 
          Last_StoreP=Current%P1
        END IF
      END IF
      p=p+1   
      IF (ASSOCIATED(Current%Next)) THEN
        Current=>Current%Next
      ELSE
        EXIT
      END IF
    END DO
    PolyLineZ%MMin=0.0d0
    PolyLineZ%MMax=1.0d0
    PolyLineZ%ContentLength=(44+4*PolyLineZ%NumParts+16*PolyLineZ%NumPoints+ &
                         16+8*PolyLineZ%NumPoints)/2
   
    xMin=MIN(PolyLineZ%xMin,xMin)
    yMin=MIN(PolyLineZ%yMin,yMin)
    xMax=MAX(PolyLineZ%xMax,xMax)
    yMax=MAX(PolyLineZ%yMax,yMax)
    zMin=PolyLineZ%zMin   !0.0d0   vordef.
    zMax=PolyLineZ%zMax   !1.0d0   vordef.
    MMin=PolyLineZ%MMin   !0.0d0   vordef.
    MMax=PolyLineZ%MMax   !1.0d0   vordef.
    LengthAllObjekte=LengthAllObjekte+(4+PolyLineZ%ContentLength)
   
    ALLOCATE(NewPoLiZ)
    IF (ASSOCIATED(ListOfPoLiZ)) THEN
      CurrPoLiZ%Next=>NewPoLiZ
      CurrPoLiZ=>NewPoLiZ
    ELSE
      ListOfPoLiZ=>NewPoLiZ
      CurrPoLiZ=>NewPoLiZ
    END IF
    CurrPoLiZ=PolyLineZ
    DEALLOCATE(PolyLineZ%Parts)
    DEALLOCATE(PolyLineZ%Points)
    DEALLOCATE(PolyLineZ%zPoints)
!   DEALLOCATE(PolyLineZ%MPoints)
     
    IF (ASSOCIATED(CurrBuild%Next)) THEN
      CurrBuild=>CurrBuild%Next
    ELSE
      EXIT
    END IF
  END DO  ! nr_oj

  ShapeType=13
  FileLengthShp=FileLengthShp+LengthAllObjekte
  CALL OpenFileShp(FileName) 
  CALL WriteHeaderShp

  CurrPoLiZ=>ListOfPoLiZ
  DO
    CALL WritePolyLineZ(CurrPoLiZ) 
    IF (ASSOCIATED(CurrPoLiZ%Next)) THEN
      CurrPoLiZ=>CurrPoLiZ%Next
    ELSE
      EXIT
    END IF
  END DO  !  ListOfPoLiZ,ListOfCurrBuild
  CALL CloseFileShp
  
  !Datei Shx
  FileLengthShx=FileLengthShx+NumberOfBuild*4 
  CALL OpenFileShx(FileName)
  CALL WriteHeaderShx
  CurrPoLiZ=>ListOfPoLiZ
  DO 
    CALL WriteRecordShx(Offset,CurrPoLiZ%ContentLength)
    Offset=OffSet+4+CurrPoLiZ%ContentLength
    IF (ASSOCIATED(CurrPoLiZ%Next)) THEN
      CurrPoLiZ=>CurrPoLiZ%Next
    ELSE
      EXIT
    END IF
  END DO
  CALL CloseFileShx

END PROGRAM MainConvETYP
