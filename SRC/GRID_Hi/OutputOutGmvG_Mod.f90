MODULE OutputOutGmvG_Mod

  USE F_Mod
  USE Floor_Mod
  USE IOControl_Mod
  USE Parametric_Mod

  USE Grid_Mod, ONLY: &
               VolFace_XY,VolFace_YZ,VolFace_ZX

  IMPLICIT NONE

CONTAINS

SUBROUTINE WriteCheckCellMaxVol(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k
  
  REAL(8) :: zw_rech
  IF (ASSOCIATED(Cell)) THEN
      Write(*,*)"!..... aus WriteCheckCellMaxVol(Cell,i,j,k) ....."
      Write(*,*)"Analyze Celle(",i,",",j,",",k,")"," :"
      Write(*,*)"        if (Cell%Face!%Vol/VolFace_!!(Cell%Face!)<1.d0-", &
                & dist_scMaxCell
      !..Face1..
      IF (zw_rech<1.d0-dist_scMaxCell) THEN  
         Write(*,*) &
           & "Cell%Face1%Vol/ VolFace_XY(Cell%Face1)=", zw_rech,"<",1.d0-dist_scMaxCell
      ELSE IF (zw_rech==1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face1%Vol/ VolFace_XY(Cell%Face1)=", zw_rech,"=",1.d0-dist_scMaxCell
      ELSE
         Write(*,*) &
           & "Cell%Face1%Vol/ VolFace_XY(Cell%Face1)=", zw_rech,">",1.d0-dist_scMaxCell
      END IF
      !..Face2..
      zw_rech=Cell%Face2%Vol/ VolFace_XY(Cell%Face2)
      IF (zw_rech<1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face2%Vol/ VolFace_XY(Cell%Face2)=", zw_rech,"<",1.d0-dist_scMaxCell
      ELSE IF (zw_rech==1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face2%Vol/ VolFace_XY(Cell%Face2)=", zw_rech,"=",1.d0-dist_scMaxCell
      ELSE
         Write(*,*) &
           & "Cell%Face2%Vol/ VolFace_XY(Cell%Face2)=", zw_rech,">",1.d0-dist_scMaxCell
      END IF
      !..Face3..
      zw_rech=Cell%Face3%Vol/ VolFace_ZX(Cell%Face3)
      IF (zw_rech<1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face3%Vol/ VolFace_ZX(Cell%Face3)=", zw_rech,"<",1.d0-dist_scMaxCell
      ELSE IF (zw_rech==1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face3%Vol/ VolFace_ZX(Cell%Face3)=", zw_rech,"=",1.d0-dist_scMaxCell
      ELSE
         Write(*,*) &
           & "Cell%Face3%Vol/ VolFace_ZX(Cell%Face3)=", zw_rech,">",1.d0-dist_scMaxCell
      END IF
      !..Face4..
      zw_rech=Cell%Face4%Vol/ VolFace_ZX(Cell%Face4)
      IF (zw_rech<1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face4%Vol/ VolFace_ZX(Cell%Face4)=", zw_rech,"<",1.d0-dist_scMaxCell
      ELSE IF (zw_rech==1.d0-dist_scMaxCell) THEN
         Write(*,*) & 
           & "Cell%Face4%Vol/ VolFace_ZX(Cell%Face4)=", zw_rech,"=",1.d0-dist_scMaxCell
      ELSE
         Write(*,*) &
           & "Cell%Face4%Vol/ VolFace_ZX(Cell%Face4)=", zw_rech,">",1.d0-dist_scMaxCell
      END IF
      !..Face5..
      zw_rech=Cell%Face5%Vol/ VolFace_YZ(Cell%Face5)
      IF (zw_rech<1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face5%Vol/ VolFace_YZ(Cell%Face5)=", zw_rech,"<",1.d0-dist_scMaxCell
      ELSE IF (zw_rech==1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face5%Vol/ VolFace_YZ(Cell%Face5)=", zw_rech,"=",1.d0-dist_scMaxCell
      ELSE
         Write(*,*) &
           & "Cell%Face5%Vol/ VolFace_YZ(Cell%Face5)=", zw_rech,">",1.d0-dist_scMaxCell
      END IF
      !..Face6..
      zw_rech=Cell%Face6%Vol/ VolFace_YZ(Cell%Face6)
      IF (zw_rech<1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face6%Vol/ VolFace_YZ(Cell%Face6)=", zw_rech,"<",1.d0-dist_scMaxCell
      ELSE IF (zw_rech==1.d0-dist_scMaxCell) THEN
         Write(*,*) &
           & "Cell%Face6%Vol/ VolFace_YZ(Cell%Face6)=", zw_rech,"=",1.d0-dist_scMaxCell
      ELSE
         Write(*,*) &
           & "Cell%Face6%Vol/ VolFace_YZ(Cell%Face6)=", zw_rech,">",1.d0-dist_scMaxCell
      END IF
      Write(*,*)"!...... Ende WriteCheckCellMaxVol ..................."
      !.....
  END IF
END SUBROUTINE WriteCheckCellMaxVol


SUBROUTINE WriteCellGMVAscii(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out

  IF (ASSOCIATED(Cell)) THEN

    IF ((Cell%in_out>-8.AND.Cell%Vol>0.0d0) &
        .OR.Cell%in_out==8) THEN
      IF (Cell%vc==0) THEN
       ! hex
        !WRITE(10,*) "Cell(i,j,k)", i, j, k," ASSOCIATED(Cell),Cell%vc==0"
        WRITE(10,'(a8,I8)') hex,8
        WRITE(10,'(8I8)') Cell%Face2%VertexList(1:4) &
                         ,Cell%Face1%VertexList(1:4)
      ELSE
        ! general
        nFaces=1
        IF (Cell%Face1%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face1%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face2%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face2%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face3%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face3%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face4%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face4%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face5%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face5%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face6%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face6%NumberVert
          nFaces=nFaces+1
        END IF
        nVerts(nFaces)=Cell%vc
        !WRITE(10,*) "Cell(i,j,k)", i, j, k," ASSOCIATED(Cell),vc>0"
        WRITE(10,'(a8,I8)') general,nfaces
        WRITE(10,*)(nVerts(iVert),iVert=1,nFaces)
        !IF(nfaces==7) THEN
        !    WRITE(10,'(7I4)') nVerts
        !ELSE 
        !    WRITE(10,*)(nVerts(iVert),iVert=1,nFaces)
        !END IF
        IF (Cell%Face1%NumberVert>2) THEN
          WRITE(10,*) Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
        END IF
        IF (Cell%Face2%NumberVert>2) THEN
          WRITE(10,*) Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
        END IF
        IF (Cell%Face3%NumberVert>2) THEN
          WRITE(10,*) Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
        END IF
        IF (Cell%Face4%NumberVert>2) THEN
          WRITE(10,*) Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
        END IF
        IF (Cell%Face5%NumberVert>2) THEN
          WRITE(10,*) Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
        END IF
        IF (Cell%Face6%NumberVert>2) THEN
          WRITE(10,*) Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
        END IF
        WRITE(10,*) Cell%VertCut(1:Cell%vc)
      END IF
    END IF
  ELSE !!! .NOT. ASSOCIATED(Cell)
    in_out=Vertices(i-1,j-1,k-1)%in_out &
          +Vertices(i,j-1,k-1)%in_out &
          +Vertices(i-1,j,k-1)%in_out &
          +Vertices(i-1,j-1,k)%in_out &
          +Vertices(i-1,j,k)%in_out &
          +Vertices(i,j-1,k)%in_out &
          +Vertices(i,j,k-1)%in_out &
          +Vertices(i,j,k)%in_out
    IF (in_out>=0) THEN
     ! hex
     !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def 
     !         /|       /|
     !        1--------2 |
     !        | 8------|-7
     !        |/       |/
     !        5--------6
      !WRITE(10,*) "Cell(i,j,k)", i, j, k," .NOT. ASSOCIATED(Cell)"
      WRITE(10,'(a8,I8)') hex,8
      WRITE(10,'(8I8)')  &
           Vertices(i-1,j-1,k)%nrP &
          ,Vertices(i,j-1,k)%nrP   &
          ,Vertices(i,j,k)%nrP   &
          ,Vertices(i-1,j,k)%nrP   &
          ,Vertices(i-1,j-1,k-1)%nrP     &
          ,Vertices(i,j-1,k-1)%nrP     &
          ,Vertices(i,j,k-1)%nrP     &
          ,Vertices(i-1,j,k-1)%nrP
    END IF
  END IF
END SUBROUTINE WriteCellGMVAscii


SUBROUTINE WriteCutPlanePolyGMVAscii(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out
  INTEGER :: mat1,mat2,mat3,nvert
  INTEGER :: li,ix,iy
  REAL(8) , POINTER :: xc(:)
  REAL(8) , POINTER :: yc(:)
  REAL(8) , POINTER :: zc(:)

  !! mat1=25 (rotbraun) ;mat2=26 (grün) ;mat3=27 (blau)
  !! (14=grau,16=lila-braun,25=rot-braun,27=blau,33=gruen)
  !! IF (ASSOCIATED(Cell)) THEN
  !!    IF (Cell%vc>0) THEN
  !!        mat1 -> alle Cell%vc
  !!        IF (k==1) und Cell%Face1%Edge(1 od.2 od. 3 od. 4)%ec
  !!             mat2  -> Schnitt Celle und am Boden
  !!    ELSE
  !!       If(k==1) und Cell%Face1%(Vert%in_out==0)
  !!              mat3 -> kein vc, maxVol, Celle am Boden 
  !!    END IF
  !! END IF
  mat1=25 !(rotbraun)
  mat2=26 !(grün)
  mat3=27 !(blau) Farbe mat1-3 -> für Var1-3
   
  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
    !~~~~~~~~~~~~~~~~~~
        ! polygons 
        mat1=25  !(rotbraun)
        nvert=Cell%vc
        !........................
        ALLOCATE(xc(1:nvert))
        ALLOCATE(yc(1:nvert))
        ALLOCATE(zc(1:nvert))
        DO li=1,nvert
          xc(li)=xParametricOut(VertOut(Cell%VertCut(li)))
          yc(li)=yParametricOut(VertOut(Cell%VertCut(li)))
          zc(li)=zParametricOut(VertOut(Cell%VertCut(li)))
        END DO
        !WRITE(10,*) "i,j,k=",i,j,k
        !WRITE(10,*) "Var1 Celle i,j,k=",i,j,k
        Write(10,*) mat1, nvert &
           ,(xc(li),li=1,nvert) &
           ,(yc(li),li=1,nvert) &
           ,(zc(li),li=1,nvert) 
        DEALLOCATE(xc)
        DEALLOCATE(yc)
        DEALLOCATE(zc)
        IF(k==1) THEN
          IF((Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
              Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) .OR. &
             (Cell%Face1%ec==-1.and.Cell%Face1%NumberVert>2)) THEN
             ! Nachtrag: diese Face muss/soll mit zur Celle-Darstellung eingebunden  
             ! werden,als Poly einzeln! 
             ! polygon: (Face1 mit Edge%ec=1 und Grundflaeche-Gelaende)
             !---1----- Specialfall 2Hügel, 
             !? ist vc>0, oder vc=0 aktuell gesetzt
             mat2=26    !(grün) 
             nvert=Cell%Face1%NumberVert
             ALLOCATE(xc(1:nvert))
             ALLOCATE(yc(1:nvert))
             ALLOCATE(zc(1:nvert))
             DO li=1,nvert
               xc(li)=xParametricOut(VertOut(Cell%Face1%VertexList(li)))
               yc(li)=yParametricOut(VertOut(Cell%Face1%VertexList(li)))
               zc(li)=zParametricOut(VertOut(Cell%Face1%VertexList(li)))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             !WRITE(10,*) "Var2 Celle i,j,k=",i,j,k
             Write(10,*) mat2, nvert &
                ,(xc(li),li=1,nvert) &
                ,(yc(li),li=1,nvert) &
                ,(zc(li),li=1,nvert)
             DEALLOCATE(xc)
             DEALLOCATE(yc)
             DEALLOCATE(zc)
          END IF
        END IF  ! cell%vc und k==1
    ELSE   ! Cell%vc==0        
    !~~~~~~~~~~~~~~~~~~
      !IF (Cell%in_out==4) THEN !Cell%Vol==maxVol oder >0.0d0 ausreichend
      !--> 23.03.2011, Betreff V8932 ,Bsp ContainerT_Verschiebe.grid geändert
      IF (k==1) THEN
         IF(Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
            Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
             ! polygon: (nur Face1)
             !---------
             mat3=27   !(blau)
             nvert=4
             !................................
             ALLOCATE(xc(1:nvert))
             ALLOCATE(yc(1:nvert))
             ALLOCATE(zc(1:nvert))
             DO li=1,nvert
               SELECT CASE (li)
                CASE (1)        !.. Vertices(i-1,j-1,k)
                    ix=i-1;iy=j-1
                CASE (2)        !.. Vertices(i,j-1,k
                    ix=i
                CASE (3)        !..Vertices(i,j,k)
                    iy=j
                CASE (4)        !..Vertices(i-1,j,k)
                    ix=i-1
               END SELECT
               xc(li)=xParametricOut(Vertices(ix,iy,k-1))
               yc(li)=yParametricOut(Vertices(ix,iy,k-1))
               zc(li)=zParametricOut(Vertices(ix,iy,k-1))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             !WRITE(10,*) "Var3 Celle i,j,k=",i,j,k
             Write(10,*) mat3, nvert &
                ,(xc(li),li=1,nvert) &
                ,(yc(li),li=1,nvert) &
                ,(zc(li),li=1,nvert)
             DEALLOCATE(xc)
             DEALLOCATE(yc)
             DEALLOCATE(zc)
         END IF
      END IF   ! (k==1)
    !ELSE IF (Cell%in_out==-4.AND.Cell%Vol==0.0d0) THEN
    !   !Cell%Vol=0.0d0 wurde in AnalyzeCell als 'special' gesetzt, Betreff Weight 
    !   IF (Cell%Face2%Edge1%Vert1%in_out==0.AND.Cell%Face2%Edge1%Vert2%in_out==0.AND. &
    !      Cell%Face2%Edge3%Vert1%in_out==0.AND.Cell%Face2%Edge3%Vert2%in_out==0.AND. &
    !      Cell%Face1%in_out==-4) THEN
    !     nFaces=1
    !     !fuer Check-Protokoll
    !     !WRITE(*,*) "i,j,k=",i,j,k,  " hex als general "
    !     WRITE(10,'(a8,I8)') general,nfaces
    !     WRITE(10,*) 4
    !     WRITE(10,*) Cell%VertCut(1:Cell%vc)
    !   END IF
    END IF  ! if-else cell%vc
    !~~~~~~~~~~~~~~~~~~~~~~~~~
  END IF   ! if associated celle

END SUBROUTINE WriteCutPlanePolyGMVAscii

SUBROUTINE WriteLandKlPlanePolyGMVAscii(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out
  INTEGER :: mat,nvert
  INTEGER :: li,ix,iy
  REAL(8) , POINTER :: xc(:)
  REAL(8) , POINTER :: yc(:)
  REAL(8) , POINTER :: zc(:)
  INTEGER :: Lk(1:9) 
  !LandKlassen: 
  !Land use:  1       2        3          4          5       6    
  !         urban savannah  deciduous  coniferous  mixed  shrubland 
  !         area            forest     forest      forest          
  !           7       8         9
  !         annual grassland  water
  !         crops
  !Zuordnung Lk:(1) (2) (3) (4) (5) (6) (7)     (8) (9)
  !             79  85  92  90  77  82  (81)89  96  76 
  Lk(:)=(/79,85,92,90,77,82,89,96,76/)
  !Lk(:)=(/74,77,82,87,90,92,96,97/)  ! Testfarben: Grün-Nuancen 
                                      ! für Wald-,Wiesen-,Sträucherdef.
  ! Im Nachtrag (7) 81 dark gray auf ( 89) braun 
  IF (ASSOCIATED(Cell)) THEN
   IF (Cell%LandClass /=0 ) THEN
    IF (Cell%vc>0) THEN
        ! polygons 
        mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
        IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
        nvert=Cell%vc
        ALLOCATE(xc(1:nvert))
        ALLOCATE(yc(1:nvert))
        ALLOCATE(zc(1:nvert))
        DO li=1,nvert
          xc(li)=xParametricOut(VertOut(Cell%VertCut(li)))
          yc(li)=yParametricOut(VertOut(Cell%VertCut(li)))
          zc(li)=zParametricOut(VertOut(Cell%VertCut(li)))
        END DO
        !WRITE(10,*) "i,j,k=",i,j,k
        Write(10,*) mat, nvert &
           ,(xc(li),li=1,nvert) &
           ,(yc(li),li=1,nvert) &
           ,(zc(li),li=1,nvert) 
        DEALLOCATE(xc)
        DEALLOCATE(yc)
        DEALLOCATE(zc)
        IF(k==1) THEN
          IF(Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
             Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) THEN
             ! polygon: (Face1 mit Edge%ec=1 und Grundflaeche-Gelaende)
             !--------- Specialfall 2Hügel, 
             !? ist vc>0, oder vc=0 aktuell gesetzt
             mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
             IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
             nvert=Cell%Face1%NumberVert
             ALLOCATE(xc(1:nvert))
             ALLOCATE(yc(1:nvert))
             ALLOCATE(zc(1:nvert))
             DO li=1,nvert
               xc(li)=xParametricOut(VertOut(Cell%Face1%VertexList(li)))
               yc(li)=yParametricOut(VertOut(Cell%Face1%VertexList(li)))
               zc(li)=zParametricOut(VertOut(Cell%Face1%VertexList(li)))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             Write(10,*) mat, nvert &
                ,(xc(li),li=1,nvert) &
                ,(yc(li),li=1,nvert) &
                ,(zc(li),li=1,nvert)
             DEALLOCATE(xc)
             DEALLOCATE(yc)
             DEALLOCATE(zc)
          END IF
        END IF
    ELSE IF & !(Cell%in_out==4) THEN !Cell%Vol==maxVol oder >0.0d0 ausreichend
       (k==1) THEN
       IF (Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
           Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
             ! polygon: (Face1)
             !---------
             mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
             IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
             nvert=4
             ALLOCATE(xc(1:nvert))
             ALLOCATE(yc(1:nvert))
             ALLOCATE(zc(1:nvert))
             DO li=1,nvert
               SELECT CASE (li)
                CASE (1)        !.. Vertices(i-1,j-1,k)
                    ix=i-1;iy=j-1
                CASE (2)        !.. Vertices(i,j-1,k
                    ix=i
                CASE (3)        !..Vertices(i,j,k)
                    iy=j
                CASE (4)        !..Vertices(i-1,j,k)
                    ix=i-1
               END SELECT
               xc(li)=xParametricOut(Vertices(ix,iy,k-1))
               yc(li)=yParametricOut(Vertices(ix,iy,k-1))
               zc(li)=zParametricOut(Vertices(ix,iy,k-1))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             Write(10,*) mat, nvert &
                ,(xc(li),li=1,nvert) &
                ,(yc(li),li=1,nvert) &
                ,(zc(li),li=1,nvert)
             DEALLOCATE(xc)
             DEALLOCATE(yc)
             DEALLOCATE(zc)
       END IF  ! IF Face1 Grenze ((V0-V4)%in_out==0)

    !1ELSE IF (Cell%in_out==-4.AND.Cell%Vol==0.0d0) THEN
    !   !Cell%Vol=0.0d0 wurde in AnalyzeCell als 'special' gesetzt, Betreff Weight 
    !   IF (Cell%Face2%Edge1%Vert1%in_out==0.AND.Cell%Face2%Edge1%Vert2%in_out==0.AND. &
    !      Cell%Face2%Edge3%Vert1%in_out==0.AND.Cell%Face2%Edge3%Vert2%in_out==0.AND. &
    !      Cell%Face1%in_out==-4) THEN
    !     nFaces=1
    !     !fuer Check-Protokoll
    !     !WRITE(*,*) "i,j,k=",i,j,k,  " hex als general "
    !     WRITE(10,'(a8,I8)') general,nfaces
    !     WRITE(10,*) 4
    !     WRITE(10,*) Cell%VertCut(1:Cell%vc)
    !   END IF
    END IF   ! IF (Cell%vc>0)  ELSE IF (k==1)
   END IF   ! IF (Cell%LandClass/=0)
  END IF   ! IF (ASSOCIATED(Cell))
END SUBROUTINE WriteLandKlPlanePolyGMVAscii

SUBROUTINE WriteCellCutPlaneGMVAscii(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out
  INTEGER :: tliste(1:4)
 
  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
        ! general
        nFaces=1
        nVerts(nFaces)=Cell%vc
        !WRITE(10,*) "i,j,k=",i,j,k
        !WRITE(10,*) "Var1 Celle i,j,k=",i,j,k
        WRITE(10,'(a8,I8)') general,nfaces
        WRITE(10,*) (nVerts(iVert),iVert=1,nFaces)
        WRITE(10,*) Cell%VertCut(1:Cell%vc)
        IF(k==1) THEN
          IF(Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
             Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) THEN
             ! polygon: (Face1,Grundflaeche-Gelaende)
             nFaces=1
             nVerts(nFaces)=Cell%Face1%NumberVert
             !WRITE(10,*) "i,j,k=",i,j,k
             !WRITE(10,*) "Var2 Celle i,j,k=",i,j,k
             WRITE(10,'(a8,I8)') general,nfaces
             WRITE(10,*) (nVerts(iVert),iVert=1,nFaces)
             WRITE(10,*) Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
          END IF
        END IF
    ELSE 
      IF (k==1) THEN 
        IF(Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
           Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0) THEN
          !Face1, Grenzflaeche-Gelaende
          nFaces=1
          !WRITE(10,*) "Var3 Celle i,j,k=",i,j,k
          WRITE(10,'(a8,I8)') general,nfaces
          WRITE(10,*) 4
          !tliste(1)=Cell%Face1%VertexList(4)
          !tliste(2)=Cell%Face1%VertexList(3)
          !tliste(3)=Cell%Face1%VertexList(2)
          !tliste(4)=Cell%Face1%VertexList(1)
          !WRITE(10,*) tliste(1:4)
          WRITE(10,*) Cell%Face1%VertexList(1:4)
        END IF
      END IF  ! IF(k==1) 
    END IF   ! IF (Cell%vc>0) ELSE
  END IF   ! IF (ASSOCIATED(Cell))
END SUBROUTINE WriteCellCutPlaneGMVAscii


SUBROUTINE WriteCellCutPlaneGMVAsciiTest2(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVertCutCell(8),iList
  INTEGER :: iVert,in_out,nrP,ic,nv
 
  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
        ! general
        nFaces=1
        nVerts(nFaces)=Cell%vc

        !WRITE(10,*) "i,j,k=",i,j,k
        WRITE(10,'(a8,I8)') general,nfaces
        WRITE(10,*) (nVerts(iVert),iVert=1,nFaces)

        DO ic=1,Cell%vc
          nrP=Cell%VertCut(ic)
          !Cell%VertCutOnlyP(ic)=VertNrCutPOut(nrP)
          Cell%VertCutOnlyP(ic)=VertOut(nrP)%nrCutP
        END DO
        WRITE(10,*) Cell%VertCutOnlyP(1:Cell%vc)
        IF(k==1) THEN
          IF(Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
             Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) THEN
             ! polygon: (Face1,Grundflaeche-Gelaende)
             nFaces=1
             nVerts(nFaces)=Cell%Face1%NumberVert
             WRITE(10,'(a8,I8)') general,nfaces
             WRITE(10,*) (nVerts(iVert),iVert=1,nFaces)
             DO nv=1,Cell%Face1%NumberVert
                nrP=Cell%Face1%VertexList(nv)
                ListVertCutCell(nv)=VertOut(nrP)%nrCutP
             END DO            
             WRITE(10,*) ListVertCutCell(1:Cell%Face1%NumberVert)
          END IF
        END IF
    ELSE 
      !IF (Cell%in_out==4.AND. &  ! Cell%in_out==4 schließt k==1 mit ein 
      ! V8.9.3.8.1.2_Entw) Kinzig_gk3_TSoilDat.grid für Soil-Ausgabe okey
      ! aber Cut/Cut2-Output muss k==1 bleiben, da i==73.and.j==75.and.k==12
      ! leere Celle entsteht obwohl VertexListe okey!!!! 
      ! auch Tauschvariante tliste und VertexList Ergebnis identisch

      IF (k==1  .and. &  ! Cell%in_out==4 schließt k==1 mit ein 
                           ! V8.9.3.8.1.2_Entw) Kinzig_gk3_TSoilDat.grid
         (Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
           Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0) )THEN
        !Face1, Grenzflaeche-Gelaende
        !WRITE(10,*) "Cell->(i,j,k)", i, j, k
        nFaces=1
        WRITE(10,'(a8,I8)') general,nfaces
        WRITE(10,*) 4   !Cell%Face1%NumberVert
        DO nv=1,Cell%Face1%NumberVert
           nrP=Cell%Face1%VertexList(nv)
           ListVertCutCell(nv)=VertOut(nrP)%nrCutP
        END DO            
        WRITE(10,*) ListVertCutCell(1:Cell%Face1%NumberVert)
      END IF
        
   !1ELSE IF (Cell%in_out==-4.AND.Cell%Vol==0.0d0) THEN 
   !   !Rand-Celle unter Berg: nur Flaeche Face2 ausgeben 
   !  IF (Cell%Face2%Edge1%Vert1%in_out==0.AND.Cell%Face2%Edge1%Vert2%in_out==0.AND. &
   !      Cell%Face2%Edge3%Vert1%in_out==0.AND.Cell%Face2%Edge3%Vert2%in_out==0.AND. &
   !      Cell%Face1%in_out==-4) THEN
   !      nFaces=1
   !      !fuer Check-Protokoll
   !      !WRITE(*,*) "i,j,k=",i,j,k,  " hex als general "
   !      WRITE(10,'(a8,I8)') general,nfaces
   !      WRITE(10,*) 4 
   !      WRITE(10,'(4I8)')  &
   !          Vertices(i-1,j-1,k)%nrCutP &
   !         ,Vertices(i-1,j,k)%nrCutP   &
   !         ,Vertices(i,j,k)%nrCutP     &
   !         ,Vertices(i,j-1,k)%nrCutP
   !  END IF
    END IF
  END IF
END SUBROUTINE WriteCellCutPlaneGMVAsciiTest2

SUBROUTINE WriteCellSoilSevenLayerGMVAsciiBin(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k
  !..............................
  INTEGER :: in_out,nv,ivc,is,isol,nrP,iNrP,anz,nrPCut,nrPF1,nrPFace
  INTEGER :: iVert,iFace,nFaces
  INTEGER :: fakt,aktfx,aktfy,aktfz,n,l,li,lix,liy,liz
  INTEGER :: nVertsF(8),ListVertCutCell(8)
  INTEGER :: ifvc
  INTEGER :: ListCut(1:2,1:10)
  INTEGER :: iv,jv,mj,mi,ii,maxcut,nrv
  INTEGER, POINTER :: fnvert(:),L_NrPCell(:,:),NrPFaces(:,:)
  INTEGER, POINTER :: LocCVertCut(:)=>NULL()
  INTEGER, POINTER :: LocCVertFace(:)=>NULL()
  INTEGER, POINTER :: LocCVertCutOnlyP(:)=>NULL()
  REAL(8) , POINTER :: cfxc(:,:)
  REAL(8) , POINTER :: cfyc(:,:)
  REAL(8) , POINTER :: cfzc(:,:)
  INTEGER :: F1V1,F1V2,F1V3,F1V4
  INTEGER :: F2V1,F2V2,F2V3,F2V4
  !..............................
  ifvc=0
  nrPFace=0
  
  IF (ASSOCIATED(Cell)) THEN   
    ! according to the selected NrB_Cells for Weight2
    IF (((Cell%in_out<5.AND.Cell%Vol>0.0d0) .OR. &
         (Cell%in_out==6.AND.Cell%vc>0) .OR. &
         (Cell%in_out==5.AND.Cell%vc>0)) &
         .and. &
        (Cell%Face1%Vol/ VolFace_XY(Cell%Face1)<1.d0-dist_scMaxCell .OR. &
         Cell%Face2%Vol/ VolFace_XY(Cell%Face2)<1.d0-dist_scMaxCell .OR. &
         Cell%Face3%Vol/ VolFace_ZX(Cell%Face3)<1.d0-dist_scMaxCell .OR. &
         Cell%Face4%Vol/ VolFace_ZX(Cell%Face4)<1.d0-dist_scMaxCell .OR. &
         Cell%Face5%Vol/ VolFace_YZ(Cell%Face5)<1.d0-dist_scMaxCell .OR. &
         Cell%Face6%Vol/ VolFace_YZ(Cell%Face6)<1.d0-dist_scMaxCell) ) THEN
         !WRITE(10,*) "i,j,k=",i,j,k, "Standard_GrenzeSoil"

        !-Soile-Boden unterhalb Cut-line
        IF(Cell%Face1%in_out==-1.and.Cell%Face1%ec==-1.and.Cell%Face1%NumberVert==3) THEN
           ! Cut_Face (Face1 und FCut)neu zusammenstellen, sollte aber schon in AnalyzeCell erfolgen
           ! Face1 und FaceCut einzeln im einem 'general' listen 
           nrPCut=Cell%vc
           nrPFace=Cell%Face1%NumberVert
           nFaces=nrPCut+2+nrPFace+2
           ALLOCATE(fnvert(1:nFaces))
           fnvert(1:2)=NrPCut
           fnvert(nrPCut+3:nrPCut+4)=nrPFace
           ifvc=1
           !WRITE(10,*) "i,j,k=",i,j,k
           !WRITE(10,*) "Var Cut unf Face1  Celle i,j,k=",i,j,k
         !ELSE IF(Cell%in_out==4.AND.Cell%Face2%in_out==4.AND.Cell%vc==0) THEN
         ELSE IF(Cell%Face1%numberVert==4.AND.Cell%Face2%numberVert==4.AND. &
                 Cell%Face1%in_out==0 .AND. Cell%vc==0) THEN
           !Face1, Grenzflaeche-Gelaende
           nrPCut=4      ! Nummer Point Face1, resultiert anz. Flächen (F3-nFace) 
           nFaces=nrPCut+2   ! nFaces = anz. Punkte(aus Celle-Face1=4) + (face1 und face2)
           ALLOCATE(fnvert(1:nFaces))
           IF(.NOT.ASSOCIATED(Cell%VertCut)) THEN
             ALLOCATE(Cell%VertCut(nrPCut))
           END IF
           IF(.NOT.ASSOCIATED(Cell%VertCutOnlyP)) THEN
             ALLOCATE(Cell%VertCutOnlyP(nrPCut))
           END IF
           fnvert(1:2)=nrPCut
           ! Face1 als Celle%VertCut nachtragen
           Cell%VertCut(1:nrPCut)=Cell%Face1%VertexList(1:nrPCut)
           !Cell%vc   !!!!nicht belegen da Orography-Abfrage noch folgt!!!!
           ifvc=1
           !WRITE(10,*) "i,j,k=",i,j,k
           !WRITE(10,*) "Celle%in_out==4.and.Cell%Face2%in_out==4, i,j,k=",i,j,k
        ELSE !Standard
          nrPCut=Cell%vc     ! Nummer Point Cell%vc, resultiert anz. Flächen (F3-nFace) 
          nFaces=nrPCut+2    ! nFaces = anz. Cell%vc + (face1 und face2)
          ALLOCATE(fnvert(1:nFaces))
          fnvert(1:2)=nrPCut
          ifvc=1
          !WRITE(10,*) "i,j,k=",i,j,k
          !WRITE(10,*) "Var nur Cut Std. Celle i,j,k=",i,j,k
        END IF
    ELSE  ! (vc==0 .OR. (vc>0.AND.Vol=0.0d0))
        !Face1, Grenzflaeche-Gelaende
        !IF (Cell%in_out==4.AND. &
        IF((k==1) .AND. & 
           (Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
            Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0) )THEN
            !Face1, Grenzflaeche-Gelaende
           nrPCut=4      ! Nummer Point Face1, resultiert anz. Flächen (F3-nFace) 
           nFaces=nrPCut+2   ! nFaces = anz. Punkte(aus Celle-Face1=4) + (face1 und face2)
           ALLOCATE(fnvert(1:nFaces))
           IF(.NOT.ASSOCIATED(Cell%VertCut)) THEN
             ALLOCATE(Cell%VertCut(nrPCut))
           END IF
           IF(.NOT.ASSOCIATED(Cell%VertCutOnlyP)) THEN
             ALLOCATE(Cell%VertCutOnlyP(nrPCut))
           END IF
           fnvert(1:2)=nrPCut
           ! Face1 als Celle%VertCut nachtragen
           Cell%VertCut(1:nrPCut)=Cell%Face1%VertexList(1:nrPCut)
           !Cell%vc   !!!!nicht belegen da Orography-Abfrage noch folgt!!!!
           ifvc=1
           !WRITE(10,*) "i,j,k=",i,j,k,"Face1_GrenzeSoil"
        ELSE
           !Wenn Check notwendig
           !Write(OutUnitProt,*) "Auswertung: \(Cell%vc==0 und nicht Cell%Face1==GrenzeCut\)"
           !Write(OutUnitProt,*) "             oder \(Cell%vc>0 und Cell%Vol==0.0\)"
           !Write(OutUnitProt,*) "i=",i," j=",j," k=",k," fakt=",fakt," nfaces=", nfaces
        END IF
        ! Speziell die nicht eingebunden werden:
        !       1) Cellen mit vc==0 .AND. (Cell%in_out==-6.OR.Cell%in_out==6)
        !          sind Cellen mit Edge als Kante
    END IF  ! IF (Cell%vc>0.AND.Cell%Vol>0.0d0.AND...<MAXVol... ELSE
    !------------------------------------------------------------------
    !------------------------------------------------------------------
    IF(ifvc==1) THEN
        ALLOCATE(LocCVertCut(nrPCut+1))      ! Local VertCut-Liste
        IF(nrPCut>nrPFace) THEN
             ALLOCATE(LocCVertCutOnlyP(nrPCut+1))  ! Local VertCutOnlyP
        ELSE 
             ALLOCATE(LocCVertCutOnlyP(nrPFace+1)) ! Local VertCutOnlyP
        END IF
        IF(nrPFace>0) THEN
             ALLOCATE(LocCVertFace(nrPFace+1))
        END IF
        ALLOCATE(L_NrPCell(1:nFaces,1:10)) ! max. 8 Points möglich
        ALLOCATE(NrPFaces(1:nFaces,1:10))  ! nrCutP je Face
        ALLOCATE(cfxc(1:nFaces,1:10))      ! max. 8 Points möglich
        ALLOCATE(cfyc(1:nFaces,1:10))
        ALLOCATE(cfzc(1:nFaces,1:10))
        L_NrPCell(:,:)=0
        NrPFaces(:,:)=0
        cfxc(:,:)=11111.11d0 !bessere Sicht Debugger
        cfyc(:,:)=11111.11d0
        cfzc(:,:)=11111.11d0

!1        IF(Cell%Face1%in_out==-1.and.Cell%Face1%ec==-1.and.Cell%Face1%NumberVert==3) THEN
!1          ! NummerVerts aus aus Face1 und CutFace in Liste
!1          DO  nrv=1, Cell%Face1%numberVert-1
!1            ListCut(1,nrv)=Cell%Face1%VertexList(nrv)
!1            ListCut(2,nrv)=Cell%Face1%VertexList(nrv+1)
!1          END DO
!1          ListCut(1,nrv)=Cell%Face1%VertexList(nrv)
!1          ListCut(2,nrv)=Cell%Face1%VertexList(1)
!1          DO ii=1,Cell%vc-1 
!1            nrv=nrv+1
!1            ListCut(1,nrv)=Cell%VertCut(ii)
!1            ListCut(2,nrv)=Cell%VertCut(ii+1)
!1          END DO
!1          nrv=nrv+1
!1          ListCut(1,nrv)=Cell%VertCut(Cell%vc)
!1          ListCut(2,nrv)=Cell%VertCut(1)
!1          !nrcut neu, suche angrenzendes Edge der beiden Faces
!1          maxcut=Cell%vc+Cell%Face1%NumberVert
!1          mj=0
!1          S1:DO iv=1,Cell%Face1%numberVert
!1           S2:DO jv=Cell%Face1%numberVert+1,maxcut
!1             IF ( (ListCut(1,iv)==ListCut(1,jv).AND.ListCut(2,iv)==ListCut(2,jv)).OR. &
!1                  (ListCut(1,iv)==ListCut(2,jv).AND.ListCut(2,iv)==ListCut(1,jv)) ) THEN
!1               ListCut(:,iv)=0   !ausblenden angrenzendes Edge
!1               ListCut(:,jv)=0   !ausblenden angrenzendes Edge
!1               mi=iv
!1               mj=jv
!1               EXIT S2
!1             END IF
!1           END DO S2
!1           if (mj/= 0) THEN
!1             EXIT S1
!1           end if
!1          END DO S1
!1          !Liste Cut aktuallisieren, NrPCuts um eine Position vorschieben
!1          !erste Null-Postion überschreiben
!1          DO iv=mi,maxcut-1
!1            ListCut(1,iv)=ListCut(1,iv+1)
!1            ListCut(2,iv)=ListCut(2,iv+1)
!1          END DO 
!1          maxcut=maxcut-1
!1          !zweite Null-Postion überschreiben
!1          DO jv=mj-1,maxcut-1
!1            ListCut(1,jv)=ListCut(1,jv+1)
!1            ListCut(2,jv)=ListCut(2,jv+1)
!1          END DO 
!1          maxcut=maxcut-1
!1          ! Celle-VertCut belegen, hier nur LocCVertCut (Local)
!1          LocCVertCut(1)= ListCut(1,1)
!1          LocCVertCut(2)= ListCut(2,1)
!1          iv=2
!1          S3:DO
!1           S4:DO is=1,maxcut
!1             IF (LocCVertCut(iv)==ListCut(1,is)) THEN
!1               iv=iv+1
!1               LocCVertCut(iv)=ListCut(2,is)
!1               ListCut(:,is)=0
!1               EXIT S4
!1             ELSE IF (LocCVertCut(iv)==ListCut(2,is)) THEN
!1               iv=iv+1
!1               LocCVertCut(iv)=ListCut(1,is)
!1               ListCut(:,is)=0
!1               EXIT S4
!1             END IF
!1           END DO S4
!1           IF (iv>=maxcut) THEN
!1             EXIT S3
!1           END IF
!1         END DO S3

        IF(nrPFace>0) THEN
          ! NummerVerts aus aus Face1 und CutFace in Liste
          LocCVertCut(1:nrPCut)=Cell%VertCut(1:nrPCut)
          LocCVertFace(1:nrPFace)=Cell%Face1%VertexList(1:nrPFace)
        ELSE
          LocCVertCut(1:nrPCut)=Cell%VertCut(1:nrPCut)
        END IF
        !Write(*,*) "i,j,k=",i,j,k
        !WRITE(10,*) "i,j,k=",i,j,k
        !IF(i==8.AND.j==1.AND.k==1) THEN
        !   !Check: Valley3D_1B_a_4LandClassDef.grid Version: V 8.9.2.3
        !IF (i==119.And.j==63.and.k==24) THEN
          !Check: Kinzig_gk3_TSoilDat.grid Version: V 8.9.3.8.1.2_Entw
          !Write(*,*) "i,j,k=",i,j,k
        !END IF
        DO isol=1,nr_lsoil
          DO ivc=1,nrPCut
            !wird zuvor in WriteAllCellsCutPlaneAsciiGMVTest2 beschrieben-Test-Variante
            !Achtung Cell%VertCutOnlyP  im oberen Teil bis Soil-Layer7 überschrieben-TestTeil
            nrP=LocCVertCut(ivc)      !nrP=Cell%VertCut(ivc)
            !!Cell%VertCutOnlyP(ivc)=VertNrCutPOut(nrP)
            !Cell%VertCutOnlyP(ivc)=VertOut(nrP)%nrCutP
            LocCVertCutOnlyP(ivc)=VertOut(nrP)%nrCutP
            !L_NrPCell(2,ivc)=Cell%VertCutOnlyP(ivc)+(isol-1)*nr_cutplane
            L_NrPCell(2,ivc)=LocCVertCutOnlyP(ivc)+(isol-1)*nr_cutplane
            !L_NrPCell(1,ivc)=Cell%VertCutOnlyP(ivc)+isol*nr_cutplane
            L_NrPCell(1,ivc)=LocCVertCutOnlyP(ivc)+isol*nr_cutplane
            NrPFaces(2,ivc)=L_NrPCell(2,ivc)
            NrPFaces(1,ivc)=L_NrPCell(1,ivc)
             !IF(i==8.AND.j==1.AND.k==1) THEN
             !   !Check: Valley3D_1B_a_4LandClassDef.grid Version: V 8.9.2.3
             !   !Problem mit Grid_intO hier negative Werte ->Abbruch
             !   Write(*,*) "       - ivc=", ivc 
             !   Write(*,*) "Cell%VertCutOnlyP(ivc)=", Cell%VertCutOnlyP(ivc) 
             !END IF
          END DO
          NrPFaces(2,ivc)=L_NrPCell(2,1)  !ivc = nrPCut+1
          NrPFaces(1,ivc)=L_NrPCell(1,1)  !ivc = nrPCut+1

          IF(nrPFace>0) THEN
            DO ivc=1,nrPFace
              nrP=LocCVertFace(ivc)
              LocCVertCutOnlyP(ivc)=VertOut(nrP)%nrCutP
              L_NrPCell(2,ivc)=LocCVertCutOnlyP(ivc)+(isol-1)*nr_cutplane
              L_NrPCell(1,ivc)=LocCVertCutOnlyP(ivc)+isol*nr_cutplane
              NrPFaces(nrPCut+2+2,ivc)=L_NrPCell(2,ivc)
              NrPFaces(nrPCut+2+1,ivc)=L_NrPCell(1,ivc)
            END DO
            NrPFaces(nrPCut+2+2,ivc)=L_NrPCell(2,1)
            NrPFaces(nrPCut+2+1,ivc)=L_NrPCell(1,1)
          END IF

          !Face3-Face(nFaces)
          li=1
          DO fakt=3,nrPCut+2
            cfxc(fakt,1)=cfxc(1,li)
            cfyc(fakt,1)=cfyc(1,li)
            cfzc(fakt,1)=cfzc(1,li)
            cfxc(fakt,2)=cfxc(1,li+1)
            cfyc(fakt,2)=cfyc(1,li+1)
            cfzc(fakt,2)=cfzc(1,li+1)
            cfxc(fakt,3)=cfxc(2,li+1)
            cfyc(fakt,3)=cfyc(2,li+1)
            cfzc(fakt,3)=cfzc(2,li+1)
            cfxc(fakt,4)=cfxc(2,li)
            cfyc(fakt,4)=cfyc(2,li)
            cfzc(fakt,4)=cfzc(2,li)
            NrPFaces(fakt,1)=NrPFaces(1,li)
            NrPFaces(fakt,2)=NrPFaces(1,li+1)
            NrPFaces(fakt,3)=NrPFaces(2,li+1)
            NrPFaces(fakt,4)=NrPFaces(2,li)
            fnvert(fakt)=4 !+1 !nur symbolisch; Anfangs-Endwert gleich
            li=li+1
          END DO
        
          IF(nrPFace>0) THEN
            li=1
            DO fakt=nrPCut+2+3,nFaces
              NrPFaces(fakt,1)=NrPFaces(nrPCut+2+1,li)
              NrPFaces(fakt,2)=NrPFaces(nrPCut+2+1,li+1)
              NrPFaces(fakt,3)=NrPFaces(nrPCut+2+2,li+1)
              NrPFaces(fakt,4)=NrPFaces(nrPCut+2+2,li)
              fnvert(fakt)=4 !+1 !nur symbolisch; Anfangs-Endwert gleich
              li=li+1
            END DO
          END IF 

          !Output CellSoil-Layer to file *.Soil.out.gmvG for ascii
          if(out_type=='a') then
             !WRITE(10,*) i,j,k, "Celle"
             WRITE(10,'(a8,I8)') general,nFaces
             WRITE(10,*) (fnvert(iFace),iFace=1,nFaces)
             DO iFace=1,nFaces
                !Write(10,*) "i=",i," j=",j," k=",k," fakt=",fakt, 'nfaces=', nfaces
                !WRITE(10,*) (L_NrPCell(iFace,iNrP),iNrP=1,fnvert(iFace))
                WRITE(10,*) (NrPFaces(iFace,iNrP),iNrP=1,fnvert(iFace)) !identisch, Ergebnis aus Analyze
             END DO
          else   !Output CellSoil-Layer to file *.Soil.out.gmvG for binary
             !.. WRITE(10) general,nfaces
             nRec=nRec+1
             WRITE(10,REC=nRec) general(1:4)
             nRec=nRec+1
             WRITE(10,REC=nRec) general(5:8)
             nRec=nRec+1
             WRITE(10,REC=nRec) nFaces
             !.. WRITE(10) (fnvert(iFace),iFace=1,nFaces)
             DO iFace=1,nFaces
               nRec=nRec+1
               WRITE(10,REC=nRec) fnvert(iFace)
             END DO
             !.. WRITE(10,*) (L_NrPCell(iFace,iNrP),iNrP=1,fnvert(iFace)) 
             DO iFace=1,nFaces
               DO iNrP=1,fnvert(iFace)
                  nRec=nRec+1
                  !WRITE(10,REC=nRec) L_NrPCell(iFace,iNrP)
                  WRITE(10,REC=nRec) NrPFaces(iFace,iNrP)
               END DO
             END DO
          end if

        END DO  !isol
        IF (ASSOCIATED(LocCVertCutOnlyP)) THEN
          DEALLOCATE(LocCVertCutOnlyP)
        END IF
        IF (ASSOCIATED(LocCVertCut)) THEN
          DEALLOCATE(LocCVertCut)
        END IF
        IF (ASSOCIATED(LocCVertFace)) THEN
          DEALLOCATE(LocCVertFace)
        END IF
        DEALLOCATE(fnvert)
        DEALLOCATE(L_NrPCell)
        DEALLOCATE(NrPFaces)
        DEALLOCATE(cfxc)
        DEALLOCATE(cfyc)
        DEALLOCATE(cfzc)

    END IF !IF (ifvc)
        !Speciell 2Hügel noch testen 30.09.2010
        !IF(k==1) THEN
        !  IF(Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
        !     Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) THEN
        !     ! polygon: (Face1,Grundflaeche-Gelaende)
        !     nFaces=1
        !     nVertsF(nFaces)=Cell%Face1%NumberVert
        !     WRITE(10,'(a8,I8)') general,nfaces
        !     WRITE(10,*) (nVertsF(iFace),iFace=1,nFaces)
        !     DO nv=1,Cell%Face1%NumberVert
        !        nrP=Cell%Face1%VertexList(nv)
        !        ListVertCutCell(nv)=VertOut(nrP)%nrCutP
        !     END DO            
        !     WRITE(10,*) ListVertCutCell(1:Cell%Face1%NumberVert)
        !  END IF
        !END IF
  END IF  ! (ASSOCIATED(Cell)) 
!TEST END IF
END SUBROUTINE WriteCellSoilSevenLayerGMVAsciiBin


SUBROUTINE WriteCellOroGMVAscii(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7),yFace(6)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out
  INTEGER :: tVert(10),m
 
  yFace(1:6)=0 
  IF (ASSOCIATED(Cell)) THEN
     IF (Cell%vc>0) THEN
        ! general
        !--------
     !WRITE(10,*)
     !WRITE(10,*) "Celle(i,j,k)= ",i,j,k, "                   Cell%vc=", Cell%vc
        nr_gen=nr_gen+1
        nr_cell_vc_oro=nr_cell_vc_oro+1
        nFaces=1
        IF ((MIN(Cell%Face1%Edge1%Vert1%in_out,Cell%Face1%Edge1%Vert2%in_out, &
                Cell%Face1%Edge3%Vert1%in_out,Cell%Face1%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face1%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face1%NumberVert
          nFaces=nFaces+1
          yFace(1)=1
        END IF
        IF ((MIN(Cell%Face2%Edge1%Vert1%in_out,Cell%Face2%Edge1%Vert2%in_out, &
                 Cell%Face2%Edge3%Vert1%in_out,Cell%Face2%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face2%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face2%NumberVert
          nFaces=nFaces+1
          yFace(2)=1
        END IF
        IF ((MIN(Cell%Face3%Edge1%Vert1%in_out,Cell%Face3%Edge1%Vert2%in_out, &
                 Cell%Face3%Edge3%Vert1%in_out,Cell%Face3%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face3%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face3%NumberVert
          nFaces=nFaces+1
          yFace(3)=1
        END IF
        IF ((MIN(Cell%Face4%Edge1%Vert1%in_out,Cell%Face4%Edge1%Vert2%in_out, &
                 Cell%Face4%Edge3%Vert1%in_out,Cell%Face4%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face4%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face4%NumberVert
          nFaces=nFaces+1
          yFace(4)=1
        END IF
        IF ((MIN(Cell%Face5%Edge1%Vert1%in_out,Cell%Face5%Edge1%Vert2%in_out, &
                 Cell%Face5%Edge3%Vert1%in_out,Cell%Face5%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face5%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face5%NumberVert
          nFaces=nFaces+1
          yFace(5)=1
        END IF
        IF ((MIN(Cell%Face6%Edge1%Vert1%in_out,Cell%Face6%Edge1%Vert2%in_out, &
                 Cell%Face6%Edge3%Vert1%in_out,Cell%Face6%Edge3%Vert2%in_out)<=0) &
            .AND.Cell%Face6%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face6%NumberVert
          nFaces=nFaces+1
          yFace(6)=1
        END IF
        nVerts(nFaces)=Cell%vc
        IF (Cell%Face1%Edge1%Vert1%in_out==0.AND. &
            Cell%Face1%Edge1%Vert2%in_out==0.AND. &
            Cell%Face1%Edge3%Vert1%in_out==0.AND. &
            Cell%Face1%Edge3%Vert2%in_out==0 ) THEN  
            ! Speziell Celle-Face1 Grenze Berg ist
              nfaces=nfaces-1
        END IF

        WRITE(10,'(a8,I8)') general,nfaces
        WRITE(10,*) (nVerts(iVert),iVert=1,nFaces)
        IF (yFace(1)==1) THEN
          WRITE(10,*) Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
        END IF
        IF (yFace(2)==1) THEN
          WRITE(10,*) Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
        END IF
        IF (yFace(3)==1) THEN
          WRITE(10,*) Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
        END IF
        IF (yFace(4)==1) THEN
          WRITE(10,*) Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
        END IF
        IF (yFace(5)==1) THEN
          WRITE(10,*) Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
        END IF
        IF (yFace(6)==1) THEN
          WRITE(10,*) Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
        END IF
        IF (nfaces/=1) THEN
          WRITE(10,*) Cell%VertCut(1:Cell%vc)
        END IF

     ELSE  ! Allokiert und vc==0 
        IF (Cell%Face1%in_out==0.AND.(Cell%Face2%in_out>=2.AND.Cell%Face2%in_out<=4)) THEN
          ! Schnittflaechen: Face1 komplett
          ! Face2%in_out=[2..4] 
          ! wird in SortVertCutCell vc=0 gesetzt, nur Face1 ausgeben! 
          ! general
          !--------
          IF(k==Domain%iz0+1) THEN
       !WRITE(10,*)
       !WRITE(10,*) "Celle(i,j,k)= ",i,j,k, "      no vc,  Face1%in_out==0.AND.Face2(in_out=2..4"
            nr_grenzeF1=nr_grenzeF1+1
            nr_cell_novc_oro=nr_cell_novc_oro+1
            WRITE(10,'(a8,I8)') general,1
            WRITE(10,*) Cell%Face1%NumberVert
            WRITE(10,*) Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
          END IF
        ELSE IF(MAX(Cell%Face2%Edge1%Vert1%in_out,Cell%Face2%Edge1%Vert2%in_out, &
                    Cell%Face2%Edge3%Vert1%in_out,Cell%Face2%Edge3%Vert2%in_out)==0) THEN
          ! Schnittflaechen: Face2 .OR. 1Edge .OR. 1Vertex
          ! wird in SortVertCutCell vc=0 gesetzt, Celle als hex ausgeben
          ! Face1%in_out=[-1 bis -4], Face2%in_out=[-3 bis 0]  
          ! hex
          !-----
          !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def
          !         /|       /|
          !        1--------2 |
          !        | 8------|-7
          !        |/       |/
          !        5--------6
          nr_grenzeF2=nr_grenzeF2+1
          nr_cell_novc_oro=nr_cell_novc_oro+1
       !WRITE(10,*)
       !WRITE(10,*) "Celle(i,j,k)= ",i,j,k,"     no vc,  MAX(Cell%Face2%Edge...%Vert..%in_out)==0"  
          WRITE(10,'(a8,I8)') hex,8
          WRITE(10,'(8I8)')  &
           Vertices(i-1,j-1,k)%nrInP &
           ,Vertices(i,j-1,k)%nrInP   &
           ,Vertices(i,j,k)%nrInP   &
           ,Vertices(i-1,j,k)%nrInP   &
           ,Vertices(i-1,j-1,k-1)%nrInP     &
           ,Vertices(i,j-1,k-1)%nrInP     &
           ,Vertices(i,j,k-1)%nrInP     &
           ,Vertices(i-1,j,k-1)%nrInP
        ELSE IF (Cell%Face2%in_out==-4 .AND. &
                 Cell%Face1%in_out==-1 .OR. Cell%Face1%in_out==-2) THEN
          nr_grenzeSp3=nr_grenzeSp3+1
          nr_cell_novc_oro=nr_cell_novc_oro+1
       !WRITE(10,*)
       !WRITE(10,*) "Celle(i,j,k)= ",i,j,k,"     no vc,  Special3-(Cell%Face2%in_out)==-4 und"
       !WRITE(10,*) "                                              Cell%Face1%in_out==-1 or. -2"  
          WRITE(10,'(a8,I8)') hex,8
          WRITE(10,'(8I8)')  &
           Vertices(i-1,j-1,k)%nrInP &
           ,Vertices(i,j-1,k)%nrInP   &
           ,Vertices(i,j,k)%nrInP   &
           ,Vertices(i-1,j,k)%nrInP   &
           ,Vertices(i-1,j-1,k-1)%nrInP     &
           ,Vertices(i,j-1,k-1)%nrInP     &
           ,Vertices(i,j,k-1)%nrInP     &
           ,Vertices(i-1,j,k-1)%nrInP
 
        END IF ! Spezialfälle, Celle allokiert, vc=0
     END IF  ! if( vc>0)  else

  ELSE !!! .NOT. ASSOCIATED(Cell) 
     in_out=Vertices(i-1,j-1,k-1)%in_out &
           +Vertices(i,j-1,k-1)%in_out &
           +Vertices(i-1,j,k-1)%in_out &
           +Vertices(i-1,j-1,k)%in_out &
           +Vertices(i-1,j,k)%in_out &
           +Vertices(i,j-1,k)%in_out &
           +Vertices(i,j,k-1)%in_out &
           +Vertices(i,j,k)%in_out
     IF (in_out==-8) THEN
       ! hex
       !----
       !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def
       !         /|       /|
       !        1--------2 |
       !        | 8------|-7
       !        |/       |/
       !        5--------6
     !WRITE(10,*)
     !WRITE(10,*) "i,j,k= ",i,j,k,"               .NOT. ASSOCIATED(Cell), Cell%in_out==-8"
       nr_insidehex=nr_insidehex+1
       WRITE(10,'(a8,I8)') hex,8
       WRITE(10,'(8I8)')  &
          Vertices(i-1,j-1,k)%nrInP &
           ,Vertices(i,j-1,k)%nrInP   &
           ,Vertices(i,j,k)%nrInP   &
           ,Vertices(i-1,j,k)%nrInP   &
           ,Vertices(i-1,j-1,k-1)%nrInP     &
           ,Vertices(i,j-1,k-1)%nrInP     &
           ,Vertices(i,j,k-1)%nrInP     &
           ,Vertices(i-1,j,k-1)%nrInP
     END IF
  END IF
END SUBROUTINE WriteCellOroGMVAscii


SUBROUTINE WriteOroCutPlanePolyGMVAscii(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  ! local:
  INTEGER :: mat,nvert
  INTEGER :: li,ix,iy
  REAL(8) , POINTER :: xc(:)
  REAL(8) , POINTER :: yc(:)
  REAL(8) , POINTER :: zc(:)
 
  mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
  IF (ASSOCIATED(Cell)) THEN
      IF (Cell%vc>0) THEN
         ! polygon: 
         !---------
         mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
         nvert=Cell%vc
         ALLOCATE(xc(1:nvert))
         ALLOCATE(yc(1:nvert))
         ALLOCATE(zc(1:nvert))
         DO li=1,nvert
           !xc(li)=VertIn(Cell%VertCut(li))%Point%x
           !yc(li)=VertIn(Cell%VertCut(li))%Point%y
           !zc(li)=VertIn(Cell%VertCut(li))%Point%z
           xc(li)=xParametricOut(VertIn(Cell%VertCut(li)))
           yc(li)=yParametricOut(VertIn(Cell%VertCut(li)))
           zc(li)=zParametricOut(VertIn(Cell%VertCut(li)))
         END DO
         !WRITE(10,*) "i,j,k=",i,j,k
         Write(10,*) mat, nvert &
            ,(xc(li),li=1,nvert) &
            ,(yc(li),li=1,nvert) &
            ,(zc(li),li=1,nvert) 
         DEALLOCATE(xc)
         DEALLOCATE(yc)
         DEALLOCATE(zc)
      ELSE IF (Cell%vc==0) THEN

          ALLOCATE(xc(1:4))
          ALLOCATE(yc(1:4))
          ALLOCATE(zc(1:4))
          IF(Cell%Face1%Edge1%Vert1%in_out==0.AND. &
             Cell%Face1%Edge1%Vert2%in_out==0.AND. &
             Cell%Face1%Edge3%Vert1%in_out==0.AND. &
             Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
             !mat=26   !mat=26 !zum check
             nvert=Cell%Face1%NumberVert !4
             !Write(*,*) "Cell%Face1%NumberVert=",Cell%Face1%NumberVert
             DO li=1,nvert
               !VertexList beinhaltet VertIn-Werte
               xc(li)=xParametricOut(VertIn(Cell%Face1%VertexList(li)))
               yc(li)=yParametricOut(VertIn(Cell%Face1%VertexList(li)))
               zc(li)=zParametricOut(VertIn(Cell%Face1%VertexList(li)))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             Write(10,*) mat, nvert &
                ,(xc(li),li=1,nvert) &
                ,(yc(li),li=1,nvert) &
                ,(zc(li),li=1,nvert) 
          END IF
          !ELSE IF(Cell%in_out==-4.AND. &
          !   (Cell%Face2%in_out==0.OR.Cell%Face3%in_out==0.OR. &
          !    Cell%Face4%in_out==0.OR.Cell%Face5%in_out==0.OR. &
          !    Cell%Face6%in_out==0) )THEN
          !   mat=27 !mat=27 ! zum check
          !   IF (Cell%Face2%in_out==0) THEN
          IF(Cell%Face2%Edge1%Vert1%in_out==0.AND.Cell%Face2%Edge1%Vert2%in_out==0.AND. &
             Cell%Face2%Edge3%Vert1%in_out==0.AND.Cell%Face2%Edge3%Vert2%in_out==0 ) THEN
             !mat=27 !mat=27 ! zum check
             nvert=Cell%Face2%NumberVert !4
             !Write(*,*) "Cell%Face2%NumberVert=",Cell%Face2%NumberVert
             DO li=1,nvert
               !VertexList beinhaltet VertIn-Werte
               xc(li)=xParametricOut(VertIn(Cell%Face2%VertexList(li)))
               yc(li)=yParametricOut(VertIn(Cell%Face2%VertexList(li)))
               zc(li)=zParametricOut(VertIn(Cell%Face2%VertexList(li)))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             Write(10,*) mat, nvert &
                ,(xc(li),li=1,nvert) &
                ,(yc(li),li=1,nvert) &
                ,(zc(li),li=1,nvert) 
          END IF
          !ELSE IF (Cell%Face3%in_out==0) THEN
          IF(Cell%Face3%Edge1%Vert1%in_out==0.AND. &
                  Cell%Face3%Edge1%Vert2%in_out==0.AND. &
                  Cell%Face3%Edge3%Vert1%in_out==0.AND. &
                  Cell%Face3%Edge3%Vert2%in_out==0 ) THEN
             !mat=28 !mat=28 ! zum check
             nvert=Cell%Face3%NumberVert !4
             !Write(*,*) "Cell%Face3%NumberVert=",Cell%Face3%NumberVert
             DO li=1,nvert
               !VertexList beinhaltet VertIn-Werte
               xc(li)=xParametricOut(VertIn(Cell%Face3%VertexList(li)))
               yc(li)=yParametricOut(VertIn(Cell%Face3%VertexList(li)))
               zc(li)=zParametricOut(VertIn(Cell%Face3%VertexList(li)))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             Write(10,*) mat, nvert &
               ,(xc(li),li=1,nvert) &
               ,(yc(li),li=1,nvert) &
               ,(zc(li),li=1,nvert) 
          END IF
          !ELSE IF (Cell%Face4%in_out==0) THEN
          IF(Cell%Face4%Edge1%Vert1%in_out==0.AND. &
                  Cell%Face4%Edge1%Vert2%in_out==0.AND. &
                  Cell%Face4%Edge3%Vert1%in_out==0.AND. &
                  Cell%Face4%Edge3%Vert2%in_out==0 ) THEN
             !mat=29 !mat=29 ! zum check
             nvert=Cell%Face4%NumberVert !4
             !Write(*,*) "Cell%Face4%NumberVert=",Cell%Face4%NumberVert
             DO li=1,nvert
               !VertexList beinhaltet VertIn-Werte
               xc(li)=xParametricOut(VertIn(Cell%Face4%VertexList(li)))
               yc(li)=yParametricOut(VertIn(Cell%Face4%VertexList(li)))
               zc(li)=zParametricOut(VertIn(Cell%Face4%VertexList(li)))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             Write(10,*) mat, nvert &
                ,(xc(li),li=1,nvert) &
                ,(yc(li),li=1,nvert) &
                ,(zc(li),li=1,nvert) 
          END IF
          !ELSE IF (Cell%Face5%in_out==0) THEN
          IF(Cell%Face5%Edge1%Vert1%in_out==0.AND. &
                  Cell%Face5%Edge1%Vert2%in_out==0.AND. &
                  Cell%Face5%Edge3%Vert1%in_out==0.AND. &
                  Cell%Face5%Edge3%Vert2%in_out==0 ) THEN
             !mat=30 !mat=30 ! zum check
             nvert=Cell%Face5%NumberVert !4
             !Write(*,*) "Cell%Face5%NumberVert=",Cell%Face5%NumberVert
             DO li=1,nvert
               !VertexList beinhaltet VertIn-Werte
               xc(li)=xParametricOut(VertIn(Cell%Face5%VertexList(li)))
               yc(li)=yParametricOut(VertIn(Cell%Face5%VertexList(li)))
               zc(li)=zParametricOut(VertIn(Cell%Face5%VertexList(li)))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             Write(10,*) mat, nvert &
                 ,(xc(li),li=1,nvert) &
                 ,(yc(li),li=1,nvert) &
                 ,(zc(li),li=1,nvert) 
          END IF
          !ELSE IF (Cell%Face6%in_out==0) THEN
          IF(Cell%Face6%Edge1%Vert1%in_out==0.AND. &
                  Cell%Face6%Edge1%Vert2%in_out==0.AND. &
                  Cell%Face6%Edge3%Vert1%in_out==0.AND. &
                  Cell%Face6%Edge3%Vert2%in_out==0 ) THEN
             !mat=31 !mat=31 ! zum check
             nvert=Cell%Face6%NumberVert !4
             !Write(*,*) "Cell%Face6%NumberVert=",Cell%Face6%NumberVert
             DO li=1,nvert
               !VertexList beinhaltet VertIn-Werte
               xc(li)=xParametricOut(VertIn(Cell%Face6%VertexList(li)))
               yc(li)=yParametricOut(VertIn(Cell%Face6%VertexList(li)))
               zc(li)=zParametricOut(VertIn(Cell%Face6%VertexList(li)))
             END DO
             !WRITE(10,*) "i,j,k=",i,j,k
             Write(10,*) mat, nvert &
                ,(xc(li),li=1,nvert) &
                ,(yc(li),li=1,nvert) &
                ,(zc(li),li=1,nvert) 
          END IF
          DEALLOCATE(xc)
          DEALLOCATE(yc)
          DEALLOCATE(zc)
                
      END IF
  END IF
END SUBROUTINE WriteOroCutPlanePolyGMVAscii

SUBROUTINE W_FaceSoilGMVAscii(Face,fmat,soiltype,i,j,k)
  TYPE (Face_T), POINTER :: Face
  INTEGER :: fmat,soiltype(1:8)
  INTEGER :: i,j,k

  INTEGER :: li,m,nvert,mat,i_st
  INTEGER :: posvert(1:8)
  REAL(8) , POINTER :: xc(:)
  REAL(8) , POINTER :: yc(:)
  REAL(8) , POINTER :: zc(:)
  REAL(8) :: dzi
  TYPE (Vertex_T) :: Fv1,Fvert(1:8)

           nvert=0
           posvert(:)=0
           IF (Face%Edge1%Vert1%in_out<=0) THEN
                   posvert(1)=1;nvert=nvert+1
           END IF
           IF (Face%Edge1%yes_sp==1) THEN
                   posvert(2)=1;nvert=nvert+1
           END IF
           IF (Face%Edge2%Vert1%in_out<=0) THEN
                   posvert(3)=1;nvert=nvert+1
           END IF
           IF (Face%Edge2%yes_sp==1) THEN
                   posvert(4)=1;nvert=nvert+1
           END IF
           IF (Face%Edge3%Vert2%in_out<=0) THEN
                   posvert(5)=1;nvert=nvert+1
           END IF
           IF (Face%Edge3%yes_sp==1) THEN
                   posvert(6)=1;nvert=nvert+1
           END IF
           IF (Face%Edge4%Vert2%in_out<=0) THEN
                   posvert(7)=1;nvert=nvert+1
           END IF
           IF (Face%Edge4%yes_sp==1) THEN
                   posvert(8)=1;nvert=nvert+1
           END IF
           nvert=nvert+1 !für Polygon EndVert (==AnfangsVert)
           IF (ASSOCIATED(Face%Edge1%Vert1)) Fvert(1)=Face%Edge1%Vert1
           IF (ASSOCIATED(Face%Edge1%VertS)) Fvert(2)=Face%Edge1%VertS
           IF (ASSOCIATED(Face%Edge2%Vert1)) Fvert(3)=Face%Edge2%Vert1
           IF (ASSOCIATED(Face%Edge2%VertS)) Fvert(4)=Face%Edge2%VertS
           IF (ASSOCIATED(Face%Edge3%Vert2)) Fvert(5)=Face%Edge3%Vert2
           IF (ASSOCIATED(Face%Edge3%VertS)) Fvert(6)=Face%Edge3%VertS
           IF (ASSOCIATED(Face%Edge4%Vert2)) Fvert(7)=Face%Edge4%Vert2
           IF (ASSOCIATED(Face%Edge4%VertS)) Fvert(8)=Face%Edge4%VertS
           DO m=1,8
             IF (posvert(m)==1) THEN
               Fv1=Fvert(m)
               EXIT
             END IF
           END DO
           DO li=1,nvert
             ALLOCATE(xc(1:nvert))
             ALLOCATE(yc(1:nvert))
             ALLOCATE(zc(1:nvert))
           END DO
           li=0
           IF (Face%Edge1%Vert1%in_out<=0) THEN
              li=li+1
              xc(li)=xParametricOut(Face%Edge1%Vert1)
              yc(li)=yParametricOut(Face%Edge1%Vert1)
              zc(li)=zParametricOut(Face%Edge1%Vert1)
           END IF
           IF (Face%Edge1%yes_sp==1)  THEN
              li=li+1
              xc(li)=xParametricOut(Face%Edge1%VertS)
              yc(li)=yParametricOut(Face%Edge1%VertS)
              zc(li)=zParametricOut(Face%Edge1%VertS)
           END IF
           IF (Face%Edge2%Vert1%in_out<=0) THEN
              li=li+1
              xc(li)=xParametricOut(Face%Edge2%Vert1)
              yc(li)=yParametricOut(Face%Edge2%Vert1)
              zc(li)=zParametricOut(Face%Edge2%Vert1)
           END IF
           IF (Face%Edge2%yes_sp==1) THEN
              li=li+1
              xc(li)=xParametricOut(Face%Edge2%VertS)
              yc(li)=yParametricOut(Face%Edge2%VertS)
              zc(li)=zParametricOut(Face%Edge2%VertS)
           END IF
           IF (Face%Edge3%Vert2%in_out<=0)  THEN
              li=li+1
              xc(li)=xParametricOut(Face%Edge3%Vert2)
              yc(li)=yParametricOut(Face%Edge3%Vert2)
              zc(li)=zParametricOut(Face%Edge3%Vert2)
           END IF
           IF (Face%Edge3%yes_sp==1) THEN
              li=li+1
              xc(li)=xParametricOut(Face%Edge3%VertS)
              yc(li)=yParametricOut(Face%Edge3%VertS)
              zc(li)=zParametricOut(Face%Edge3%VertS)
           END IF
           IF (Face%Edge4%Vert2%in_out<=0) THEN
              li=li+1
              xc(li)=xParametricOut(Face%Edge4%Vert2)
              yc(li)=yParametricOut(Face%Edge4%Vert2)
              zc(li)=zParametricOut(Face%Edge4%Vert2)
           END IF
           IF (Face%Edge4%yes_sp==1) THEN
              li=li+1
              xc(li)=xParametricOut(Face%Edge4%VertS)
              yc(li)=yParametricOut(Face%Edge4%VertS)
              zc(li)=zParametricOut(Face%Edge4%VertS)
           END IF
           li=li+1
           xc(li)=xParametricOut(Fv1)
           yc(li)=yParametricOut(Fv1)
           zc(li)=zParametricOut(Fv1)

           !Test Face3 mat=10   !mat=9
           !Write(10,*) fmat, nvert &
           !       ,(xc(li),li=1,nvert) &
           !       ,(yc(li),li=1,nvert) &
           !       ,(zc(li),li=1,nvert)
           !zc=0
           DO i_st=1,8         
              Write(10,*) soiltype(i_st), nvert &
                  ,(xc(li),li=1,nvert) &
                  ,(yc(li),li=1,nvert) &
                  ,(zc(li),li=1,nvert)
              zc=zc-20
           END DO
           DEALLOCATE(xc)
           DEALLOCATE(yc)
           DEALLOCATE(zc)
           
           ! Test für Face3 als Fz0-Darstellung
            !Fz0-Face3 Fvert(1)=Vertices(i-1,j-1,iz0)
            !Fz0-Face3 Fvert(2)=Vertices(i  ,j-1,iz0)
            !Fz0-Face3 Fvert(3)=Vertices(i  ,j-1,iz0)
            !Fz0-Face3 Fvert(4)=Vertices(i-1,j-1,iz0)
            !Fz0-Face3 dzi=(Vertices(i-1,j-1,iz0+1)%Point%z-Vertices(i-1,j-1,iz0)%Point%z)/7
            !Fz0-Face3 Fvert(1:2)%Point%z=Fvert(1:2)%Point%z-dzi
            !Fz0-Face3 nvert=4
            !Fz0-Face3 DO li=1,nvert
            !Fz0-Face3   ALLOCATE(xc(1:nvert))
            !Fz0-Face3   ALLOCATE(yc(1:nvert))
            !Fz0-Face3   ALLOCATE(zc(1:nvert))
            !Fz0-Face3 END DO
            !Fz0-Face3 !def. soil
            !Fz0-Face3 DO i_st=1,7         
            !Fz0-Face3    DO li=1,nvert
            !Fz0-Face3      xc(li)=xParametricOut(Fvert(li))
            !Fz0-Face3      yc(li)=yParametricOut(Fvert(li))
            !Fz0-Face3      zc(li)=zParametricOut(Fvert(li))
            !Fz0-Face3    END DO
            !Fz0-Face3    Write(10,*) soiltype(i_st), nvert &
            !Fz0-Face3        ,(xc(li),li=1,nvert) &
            !Fz0-Face3        ,(yc(li),li=1,nvert) &
            !Fz0-Face3        ,(zc(li),li=1,nvert)
            !Fz0-Face3    Fvert(:)%Point%z=Fvert(:)%Point%z-dzi
            !Fz0-Face3 END DO
            !Fz0-Face3 DEALLOCATE(xc)
            !Fz0-Face3 DEALLOCATE(yc)
            !Fz0-Face3 DEALLOCATE(zc)

END SUBROUTINE W_FaceSoilGMVAscii

SUBROUTINE W_Fz0SoilGMVAscii(NrF,i,j,k,soiltype)
  INTEGER :: NrF
  INTEGER :: i,j,k
  INTEGER :: soiltype(1:8)

  INTEGER :: li,nvert,i_st
  REAL(8) , POINTER :: xc(:)
  REAL(8) , POINTER :: yc(:)
  REAL(8) , POINTER :: zc(:)
  REAL(8) :: dzi
  TYPE (Vertex_T) :: Fvert(1:8),V(1:2)
        
  !dzi auf Grundgitter bezogen angepasst
  !7 Schichten im Boden
  !V(1:2)%Point%x=Domain%xP(0)
  !V(1:2)%Point%y=Domain%yP(0)
  !V(1)%Point%z=Domain%zP(0)
  !V(2)%Point%z=Domain%zP(1)
  !dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
  dzi=dzi_soil
  SELECT CASE (NrF)
       CASE (1)
          Fvert(1)=Vertices(i-1,j-1,iz0)
          Fvert(2)=Vertices(i  ,j-1,iz0)
          Fvert(3)=Vertices(i  ,j  ,iz0)
          Fvert(4)=Vertices(i-1,j  ,iz0)
          !iz0 Face1/Fvert(1:4) z0 verschieben
          Fvert(1:4)%Point%z=Fvert(1:4)%Point%z-dzi
       CASE (2)
          Fvert(1)=Vertices(i-1,j-1,iz0)
          Fvert(2)=Vertices(i  ,j-1,iz0)
          Fvert(3)=Vertices(i  ,j  ,iz0)
          Fvert(4)=Vertices(i-1,j  ,iz0)
          !iz0 Face2 z0 bleibt
       CASE (3)
          Fvert(1)=Vertices(i-1,j-1,iz0)
          Fvert(2)=Vertices(i  ,j-1,iz0)
          Fvert(3)=Vertices(i  ,j-1,iz0)
          Fvert(4)=Vertices(i-1,j-1,iz0)
          !iz0 Face3/Fvert(1:2) z0 verschieben
          Fvert(1:2)%Point%z=Fvert(1:2)%Point%z-dzi
       CASE (4)
          Fvert(1)=Vertices(i-1,j,iz0)
          Fvert(2)=Vertices(i  ,j,iz0)
          Fvert(3)=Vertices(i  ,j,iz0)
          Fvert(4)=Vertices(i-1,j,iz0)
          !iz0 Face4/Fvert(1:2) z0 verschieben
          Fvert(1:2)%Point%z=Fvert(1:2)%Point%z-dzi
       CASE (5)
          Fvert(1)=Vertices(i-1,j-1,iz0)
          Fvert(2)=Vertices(i-1,j  ,iz0)
          Fvert(3)=Vertices(i-1,j  ,iz0)
          Fvert(4)=Vertices(i-1,j-1,iz0)
          !iz0 Face4/Fvert(1:2) z0 verschieben
          Fvert(1:2)%Point%z=Fvert(1:2)%Point%z-dzi
       CASE (6)
          Fvert(1)=Vertices(i,j-1,iz0)
          Fvert(2)=Vertices(i,j  ,iz0)
          Fvert(3)=Vertices(i,j  ,iz0)
          Fvert(4)=Vertices(i,j-1,iz0)
          !iz0 Face4/Fvert(1:2) z0 verschieben
          Fvert(1:2)%Point%z=Fvert(1:2)%Point%z-dzi
   END SELECT
 
   nvert=4
   DO li=1,nvert
     ALLOCATE(xc(1:nvert))
     ALLOCATE(yc(1:nvert))
     ALLOCATE(zc(1:nvert))
   END DO
   !def. soil
   DO i_st=1,7         
      DO li=1,nvert
        xc(li)=xParametricOut(Fvert(li))
        yc(li)=yParametricOut(Fvert(li))
        zc(li)=zParametricOut(Fvert(li))
      END DO
      Write(10,*) soiltype(i_st), nvert &
          ,(xc(li),li=1,nvert) &
          ,(yc(li),li=1,nvert) &
          ,(zc(li),li=1,nvert)
      Fvert(:)%Point%z=Fvert(:)%Point%z-dzi
   END DO

   DEALLOCATE(xc)
   DEALLOCATE(yc)
   DEALLOCATE(zc)

END SUBROUTINE W_Fz0SoilGMVAscii


SUBROUTINE WritePolySoilSevenLayerGMVAsciiBin(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k
  !...............................
  
  INTEGER :: li,la,lix,liy,liz,ix,iy,iz
  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7),yFace(6)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out,nvert,posvert(1:8)
  INTEGER :: soiltype(1:8),i_st,mat,fmat(1:7)
  INTEGER :: st(1:9)
  INTEGER :: basis,allface,nrs,is
  INTEGER :: fv,nrf,aktfx,aktfy,aktfz,fakt,f1,f2,n,l
  INTEGER :: ifvc,ifgf,akt_vc,pos
  INTEGER, POINTER :: fnvert(:) 
  REAL(8), POINTER :: xc(:),cfxc(:,:)
  REAL(8), POINTER :: yc(:),cfyc(:,:)
  REAL(8), POINTER :: zc(:),cfzc(:,:)
  REAL(8) :: dzi
  TYPE (Vertex_T) :: V1,VCell(1:8),Fvert(1:4),V(1:2)

  !Specials: Celle k==11 ; Face2  CutFläche 3xin_out=0
         !   Celle k==12 ; Face1  Grenzfläche 3xin_out==0  und CutFace mit vc=3
         !     -> unter CellSoilSevenLayer k==12 dargestellt, muss eine Celle ergeben
         !     -> unter PolySoil  CutFläche-Face2 aus k=11 und CutFace aus k=12 dargestellt
         !        ist ok, da nicht eine Celle sein muss, Bsp Kinzig.grid V_8.9.3.8/TEST_NrB

  !soiltype: 1=ice  2=rock  3=sand  4=sandy loam  5=loam  6=clay loam
  !          7=clay  8=peat  9=water
  !          clay=Ton; peat=Torf
  !z.Zt: Zuordnung BodenType-->Farbe-View-GMV
  !      1=ice=hellblau(41) 2=rock=grau(49) 3=sand=gelb(36) 4=sandy loam=khaki(55)
  !      5=loam=rot(47) 6=clay loam(57) 7=clay=dunkelrot(43)  8=peat=dunkel_lila(13)
  !      9=water=blue(44)
  mat=16
  yFace(1:6)=0 
  fmat(:)   =(/8,9,10,11,12,13,14/)
  st(:)     =(/41,49,36,55,53,57,39,43,44/)
  !st(:)=(/74,77,82,87,90,92,96,97,97/)  ! Testfarben: Grün-Nuancen für Wald-,Wiesen-,Sträucherdef.

  !Testvarianten
  !soiltype(:)=(/41,49,132,55,47,89,43,13/)
  !soiltype(:)=(/41,17,53,132,89,47,46,13/)
  !soiltype(1:8)=25   ! 14=grau,16=lila-braun,25=rot-braun,33=gruen

  soiltype(:)=0
  ifvc=0
  IF (ASSOCIATED(Cell)) THEN
    ! according to the selected NrB_Cells for Weight2
    IF (Cell%vc>0) THEN  
        !IF (Cell%vc>0.AND.Cell%Vol>0.0d0.AND. & 
        !   (Cell%Face1%Vol/ VolFace_XY(Cell%Face1)<1.d0-1.d-2 &
        !.OR.Cell%Face2%Vol/ VolFace_XY(Cell%Face2)<1.d0-1.d-2 &
        !.OR.Cell%Face3%Vol/ VolFace_ZX(Cell%Face3)<1.d0-1.d-2 &
        !.OR.Cell%Face4%Vol/ VolFace_ZX(Cell%Face4)<1.d0-1.d-2 &
        !.OR.Cell%Face5%Vol/ VolFace_YZ(Cell%Face5)<1.d0-1.d-2 &
        !.OR.Cell%Face6%Vol/ VolFace_YZ(Cell%Face6)<1.d0-1.d-2)) THEN
           nvert=Cell%vc    ! Nummer Point Cell%vc, resultiert anz. Flächen (F3-nFace) 
           nrf=nvert+2      ! nFace = anz.Cell%vc + (face1 und face2)
           ALLOCATE(fnvert(1:nrf))
           ifvc=1
    ELSE  ! (vc==0 .OR. (vc>0.AND.Vol=0.0d0))
        IF (k==1 .and. &  !(Cell%in_out==4.AND. &
           (Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
            Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0) )THEN
           !Face1, Grenzflaeche-Gelaende
           nvert=4          ! Nummer Point Face1, resultiert anz. Flächen (F3-nFace)
           nrf=nvert+2      ! nFace = anz. Punkte(aus Celle-Face1=4) + (face1 und face2)
           ALLOCATE(fnvert(1:nrf))
           IF(.NOT.ASSOCIATED(Cell%VertCut)) THEN
                ALLOCATE(Cell%VertCut(nvert))
           END IF
           ! Face1 als Celle%VertCut nachtragen
           Cell%VertCut(1:nvert)=Cell%Face1%VertexList(1:nvert)
           !Cell%vc   !!!!nicht belegen da Orography-Abfrage noch folgt!!!!
           ifvc=1
        ELSE
           !Wenn Check notwendig
           !Write(OutUnitProt,*) "Auswertung: \(Cell%vc==0 und nicht Cell%Face1==GrenzeCut\)"
           !Write(OutUnitProt,*) "             oder \(Cell%vc>0 und Cell%Vol==0.0\)"
           !Write(OutUnitProt,*) "i=",i," j=",j," k=",k," fakt=",fakt," nfaces=", nfaces
        END IF
        ! Speziell die nicht eingebunden werden:
        !       1) Cellen mit vc==0 .AND. (Cell%in_out==-6.OR.Cell%in_out==6)
        !          sind Cellen mit Edge als Kante
    END IF
    !...........................................................
    ! 2.Variante:
    ! Beschreibt Bodenschichten entlang Orography mit Polygone
    !...........................................................
    !-Soile-Boden unterhalb Cut-line
    !-Face2 entspricht Cut-Fläche, Face2 kopiert auf Face1,
    !-Face1-zebene mit minus dzi_soil verschoben
    !-Face3-Face(nvert) entgegen Uhrzeigersinn aufgestellt
    !-Ermitteln der jeweiligen Faces, Checken gl. Points x-,bzw y-Richtung sind
    !-alle Polygone je Face erhalten Points mit Anzahl Face(nvert+1) 
    !-je Soil-Boden werden alle Faces je Ebene Soil(1-7)
    ! in z-Ebene komplett mit minus dzi_soil verschoben
    !...........................................................

    IF (ifvc==1) THEN            
      !Belege Soil-Liste mit akt.Soiltypen der akt Grund-Celle 
      !-----------------------------
      DO nrs=1,nr_soildef
        IF(i>(s_ixa(nrs)*2.e0**RefineX).AND.i<=(s_ixe(nrs)*2.e0**RefineX) .AND. &
           j>(s_iya(nrs)*2.e0**RefineY).AND.j<=(s_iye(nrs)*2.e0**RefineY)) THEN
           DO is=1,nr_lsoil !akt. 7 mit 8 deklariert 
             soiltype(is)=st(soil_type(nrs,is))
           END DO
           !WRITE(*,*) "def. soiltype : ", soiltype(1:nr_lsoil)
           EXIT
        END IF
      END DO
      dzi=dzi_soil

      !Points Face1-Face(nvert+2) zusammenstellen
      !-------------------------------------------
      ALLOCATE(cfxc(1:nrf,1:8)) ! max. 8 Points möglich
      ALLOCATE(cfyc(1:nrf,1:8))
      ALLOCATE(cfzc(1:nrf,1:8))
      cfxc(:,:)=11111.11d0 !bessere Sicht Debugger
      cfyc(:,:)=11111.11d0
      cfzc(:,:)=11111.11d0
      !-----
      !Face2
      !-----
      DO li=1,nvert
        cfxc(2,li)=xParametricOut(VertOut(Cell%VertCut(li)))
        cfyc(2,li)=yParametricOut(VertOut(Cell%VertCut(li)))
        cfzc(2,li)=zParametricOut(VertOut(Cell%VertCut(li)))
      END DO
      !für Polygone Anfangswert==Endwert; fnvert um 1 akumulieren
      cfxc(2,li)=cfxc(2,1)
      cfyc(2,li)=cfyc(2,1)
      cfzc(2,li)=cfzc(2,1)
      !-----
      !Face1
      !-----
      cfxc(1,:)=cfxc(2,:)
      cfyc(1,:)=cfyc(2,:)
      cfzc(1,:)=cfzc(2,:)-dzi_soil
      fnvert(1:2)=nvert+1  ! +1 für anhängen Anfangswert
! If(i==74.and.j==1.and.k==13) then
!    write(*,*) "cfxc(2:li) =",cfxc(2,:)
!    write(*,*) "cfxc(1:li) =",cfxc(1,:)
!    write(*,*) "cfyc(2:li) =",cfyc(2,:)
!    write(*,*) "cfyc(1:li) =",cfyc(1,:)
!    write(*,*) "cfzc(2:li) =",cfzc(2,:)
!    write(*,*) "cfzc(1:li) =",cfzc(1,:)
!  end if

      fakt=3
      f1=1;f2=2;pos=1
      DO While (fakt<=nrf)
        !.................................
        !immer 2Point(Vert) Face2/1 eine Face(fakt)
        cfxc(fakt,1) = cfxc(f1,pos) !AnfangPoint
        cfxc(fakt,2) = cfxc(f1,pos+1)
        cfxc(fakt,3) = cfxc(f2,pos+1)
        cfxc(fakt,4) = cfxc(f2,pos)
        cfxc(fakt,5) = cfxc(f1,pos) !EndPoint=AnfangPoint
        !.................................
        cfyc(fakt,1) = cfyc(f1,pos)
        cfyc(fakt,2) = cfyc(f1,pos+1)
        cfyc(fakt,3) = cfyc(f2,pos+1)
        cfyc(fakt,4) = cfyc(f2,pos)
        cfyc(fakt,5) = cfyc(f1,pos)
        !.................................
        cfzc(fakt,1) = cfzc(f1,pos)
        cfzc(fakt,2) = cfzc(f1,pos+1)
        cfzc(fakt,3) = cfzc(f2,pos+1)
        cfzc(fakt,4) = cfzc(f2,pos)
        cfzc(fakt,5) = cfzc(f1,pos)
! If(i==74.and.j==1.and.k==13) then
!    write(*,*) "cfxc(",fakt,":li) =",cfxc(fakt,:)
!    write(*,*) "cfyc(",fakt,":li) =",cfyc(fakt,:)
!    write(*,*) "cfzc(",fakt,":li) =",cfzc(fakt,:)
!  end if
        !.................................
        fnvert(fakt)=4+1 ! +1 für anhängen Anfangswert
        pos=pos+1;fakt=fakt+1
        !.................................
      END DO !DO While(fakt<=nrf)


      !Output Soil-Layer to file *out.gmvG for ascii
      if(out_type=='a') then
         DO i_st=1,nr_lsoil
           DO fv=1,nrf
              Write(10,*) soiltype(i_st), fnvert(fv) &
                 ,(cfxc(fv,li),li=1,fnvert(fv)) &
                ,(cfyc(fv,li),li=1,fnvert(fv)) &
                 ,(cfzc(fv,li),li=1,fnvert(fv))
           END DO
           DO fv=1,nrf
              cfzc(fv,:)=cfzc(fv,:)-dzi_soil
           END DO
         END DO
      else
         !Output Soil-Layer to file *out.gmvG for binary
         DO i_st=1,nr_lsoil
           DO fv=1,nrf
              nRec=nRec+1
              WRITE(10,REC=nRec) soiltype(i_st)
              nRec=nRec+1
              WRITE(10,REC=nRec) fnvert(fv)
              DO li=1,fnvert(fv)
                nRec=nRec+1
                tempw=cfxc(fv,li)
                WRITE(10,REC=nRec) tempw
              END DO
              DO li=1,fnvert(fv)
                nRec=nRec+1
                tempw=cfyc(fv,li)
                WRITE(10,REC=nRec) tempw
              END DO
              DO li=1,fnvert(fv)
                nRec=nRec+1
                tempw=cfzc(fv,li)
                WRITE(10,REC=nRec) tempw
              END DO
           END DO
           DO fv=1,nrf
              cfzc(fv,:)=cfzc(fv,:)-dzi_soil
           END DO
         END DO
      end if     

      DEALLOCATE(fnvert)
      DEALLOCATE(cfxc)
      DEALLOCATE(cfyc)
      DEALLOCATE(cfzc)

    END IF  ! IF(ifvc)
  END IF  ! If(ASSOCIATED(Cell))

END SUBROUTINE WritePolySoilSevenLayerGMVAsciiBin



SUBROUTINE WriteCellGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out

  IF (ASSOCIATED(Cell)) THEN
    IF ((Cell%in_out>-8.AND.Cell%Vol>0.0d0) &
         .OR.Cell%in_out==8) THEN
      IF (Cell%vc==0) THEN
       ! hex
        nRec=nRec+1
        WRITE(10,REC=nRec) hex(1:4)
        nRec=nRec+1
        WRITE(10,REC=nRec) hex(5:8)
        nRec=nRec+1
        WRITE(10,REC=nRec) 8
        DO iVert=1,4
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%Face2%VertexList(iVert)
        END DO
        DO iVert=1,4
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%Face1%VertexList(iVert)
        END DO
      ELSE
        ! general
        nFaces=1
        IF (Cell%Face1%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face1%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face2%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face2%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face3%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face3%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face4%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face4%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face5%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face5%NumberVert
          nFaces=nFaces+1
        END IF
        IF (Cell%Face6%NumberVert>2) THEN
          nVerts(nFaces)=Cell%Face6%NumberVert
          nFaces=nFaces+1
        END IF
        nVerts(nFaces)=Cell%vc

        !.. WRITE(10) general,nfaces
        nRec=nRec+1
        WRITE(10,REC=nRec) general(1:4)
        nRec=nRec+1
        WRITE(10,REC=nRec) general(5:8)
        nRec=nRec+1
        WRITE(10,REC=nRec) nfaces
        !.. WRITE(10) (nVerts(i),i=1,nFaces)
        DO iVert=1,nFaces
          nRec=nRec+1
          WRITE(10,REC=nRec) nVerts(iVert)
        END DO
        IF (Cell%Face1%NumberVert>2) THEN
          DO iVert=1,Cell%Face1%NumberVert
            nRec=nRec+1
            WRITE(10,REC=nRec)  Cell%Face1%VertexList(iVert)
          END DO
        END IF
        IF (Cell%Face2%NumberVert>2) THEN
          DO iVert=1,Cell%Face2%NumberVert
            nRec=nRec+1
            WRITE(10,REC=nRec) Cell%Face2%VertexList(iVert)
          END DO
        END IF
        IF (Cell%Face3%NumberVert>2) THEN
          DO iVert=1,Cell%Face3%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face3%VertexList(iVert)
          END DO
        END IF
        IF (Cell%Face4%NumberVert>2) THEN
          DO iVert=1,Cell%Face4%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face4%VertexList(iVert)
          END DO
        END IF
        IF (Cell%Face5%NumberVert>2) THEN
          DO iVert=1,Cell%Face5%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face5%VertexList(iVert)
          END DO
        END IF
        IF (Cell%Face6%NumberVert>2) THEN
          DO iVert=1,Cell%Face6%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face6%VertexList(iVert)
          END DO
        END IF
        DO iVert=1,Cell%vc
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%VertCut(iVert)
        END DO
      END IF
    END IF
  ELSE
    in_out=Vertices(i-1,j-1,k-1)%in_out &
          +Vertices(i,j-1,k-1)%in_out &
          +Vertices(i-1,j,k-1)%in_out &
          +Vertices(i-1,j-1,k)%in_out &
          +Vertices(i-1,j,k)%in_out &
          +Vertices(i,j-1,k)%in_out &
          +Vertices(i,j,k-1)%in_out &
          +Vertices(i,j,k)%in_out
    IF (in_out>=0) THEN
     ! hex
      !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def
      !         /|       /|
      !        1--------2 |
      !        | 8------|-7
      !        |/       |/
      !        5--------6
      nRec=nRec+1
      WRITE(10,REC=nRec) hex(1:4)
      nRec=nRec+1
      WRITE(10,REC=nRec) hex(5:8)
      nRec=nRec+1
      WRITE(10,REC=nRec) 8
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j-1,k)%nrP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j-1,k)%nrP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j,k)%nrP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j,k)%nrP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j-1,k-1)%nrP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j-1,k-1)%nrP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j,k-1)%nrP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j,k-1)%nrP
    END IF
  END IF
END SUBROUTINE WriteCellGMVBinary

SUBROUTINE WriteCutPlanePolyGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  ! loacal
  INTEGER :: mat1,mat2,mat3,nvert
  INTEGER :: li,ix,iy

  !! mat1=25 (rotbraun) ;mat2=26 (grün) ;mat3=27 (blau)
  !! (14=grau,16=lila-braun,25=rot-braun,27=blau,33=gruen)
  !! IF (ASSOCIATED(Cell)) THEN
  !!    IF (Cell%vc>0) THEN
  !!       mat1 -> alle Cell%vc
  !!       IF(k==1) und Cell%Face1%Edge(1 od.2 od. 3 od. 4)%ec
  !!           mat2  -> Schnitt Celle und am Boden
  !!    ELSE
  !!       If(k==1) und Cell%Face1%(Vert%in_out==0)
  !!           mat3 -> kein vc, maxVol, Celle am Boden 
  !!    END IF
  !! END IF

  mat1=25;mat2=26;mat3=27 ! Farbe mat1-3 -> für Var1-3
  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
    !~~~~~~~~~~~~~~~~~~
       ! polygons 
       mat1=25  !(rotbraun)
       nvert=Cell%vc
       !............
       nRec=nRec+1
       WRITE(10,REC=nRec) mat1
       nRec=nRec+1
       WRITE(10,REC=nRec) nvert
       DO li=1,nvert
         nRec=nRec+1
         tempw=xParametricOut(VertOut(Cell%VertCut(li)))
         WRITE(10,REC=nRec) tempw
       END DO
       DO li=1,nvert
         nRec=nRec+1
         tempw=yParametricOut(VertOut(Cell%VertCut(li)))
         WRITE(10,REC=nRec) tempw
       END DO
       DO li=1,nvert
         nRec=nRec+1
         tempw=zParametricOut(VertOut(Cell%VertCut(li)))
         WRITE(10,REC=nRec) tempw
       END DO
       IF(k==1) THEN
         IF((Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
             Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) .OR. &
            (Cell%Face1%ec==-1.and.Cell%Face1%NumberVert>2)) THEN
            ! polygon: (Face1 mit Edge%ec=1 und Grundflaeche-Gelaende)
            !--------- Specialfall 2Hügel, ContainerT.grid v8931 
            !? ist vc>0, oder vc=0 aktuell gesetzt
            mat2=26  !(grün)
            nvert=Cell%Face1%NumberVert
            !------------------------------------------------
            nRec=nRec+1
            WRITE(10,REC=nRec) mat2
            nRec=nRec+1
            WRITE(10,REC=nRec) nvert
            DO li=1,nvert
              nRec=nRec+1
              tempw=xParametricOut(VertOut(Cell%Face1%VertexList(li)))
              WRITE(10,REC=nRec) tempw
            END DO
            DO li=1,nvert
              nRec=nRec+1
              tempw=yParametricOut(VertOut(Cell%Face1%VertexList(li)))
              WRITE(10,REC=nRec) tempw
            END DO
            DO li=1,nvert
              nRec=nRec+1
              tempw=zParametricOut(VertOut(Cell%Face1%VertexList(li)))
              WRITE(10,REC=nRec) tempw
            END DO
         END IF
       END IF  ! cell%vc und k==1
    ELSE  !Cell%vc==0
    !~~~~~~~~~~~~~~~~~~~~
      !IF (Cell%in_out==4) THEN !Cell%Vol==maxVoloder >0.0d0 ausreichend
      IF (k==1) THEN
        IF(Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
           Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
          ! polygon: (nur Face1)
          !---------
          mat3=27  !(blau)
          nvert=4
          !-------------------------------------------------
          nRec=nRec+1
          WRITE(10,REC=nRec) mat3
          nRec=nRec+1
          WRITE(10,REC=nRec) nvert
          DO li=1,nvert
            SELECT CASE (li)
                 CASE (1)        !.. Vertices(i-1,j-1,k)
                     ix=i-1;iy=j-1
                 CASE (2)        !.. Vertices(i,j-1,k
                     ix=i
                 CASE (3)        !..Vertices(i,j,k)
                     iy=j
                 CASE (4)        !..Vertices(i-1,j,k)
                     ix=i-1
                END SELECT
            nRec=nRec+1
            tempw=xParametricOut(Vertices(ix,iy,k-1))
            WRITE(10,REC=nRec) tempw
          END DO
          DO li=1,nvert
            SELECT CASE (li)
                 CASE (1)        !.. Vertices(i-1,j-1,k)
                     ix=i-1;iy=j-1
                 CASE (2)        !.. Vertices(i,j-1,k
                     ix=i
                 CASE (3)        !..Vertices(i,j,k)
                     iy=j
                 CASE (4)        !..Vertices(i-1,j,k)
                     ix=i-1
                END SELECT
            nRec=nRec+1
            tempw=yParametricOut(Vertices(ix,iy,k-1))
            WRITE(10,REC=nRec) tempw
          END DO
          DO li=1,nvert
            SELECT CASE (li)
                 CASE (1)        !.. Vertices(i-1,j-1,k)
                     ix=i-1;iy=j-1
                 CASE (2)        !.. Vertices(i,j-1,k
                     ix=i
                 CASE (3)        !..Vertices(i,j,k)
                     iy=j
                 CASE (4)        !..Vertices(i-1,j,k)
                     ix=i-1
                END SELECT
            nRec=nRec+1
            tempw=zParametricOut(Vertices(ix,iy,k-1))
            WRITE(10,REC=nRec) tempw
          END DO
       END IF    ! IF Face1 Grenze ((V0-V4)%in_out==0)
      END IF ! (k==1)
    END IF  ! if-else cell%vc
    !~~~~~~~~~~~~~~~~~~~~~~~~~~~
  END IF  ! if associated celle

END SUBROUTINE WriteCutPlanePolyGMVBinary


SUBROUTINE WriteLandKlPlanePolyGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  ! loacal
  INTEGER :: mat,nvert
  INTEGER :: li,ix,iy
  INTEGER :: Lk(1:9)
  !LandKlassen: 
  !Land use:  1       2        3          4          5       6    
  !         urban savannah  deciduous  coniferous  mixed  shrubland 
  !         area            forest     forest      forest          
  !           7       8         9
  !         annual grassland  water
  !         crops
  !Zuordnung Lk:(1) (2) (3) (4) (5) (6)  (7)   (8) (9)
  !             79  85  92  90  77  82  81(89)  96  76 
  Lk(:)=(/79,85,92,90,77,82,89,96,76/)
  !Lk(:)=(/74,77,82,87,90,92,96,97/)  ! Testfarben: Grün-Nuancen 
                                      ! für Wald-,Wiesen-,Sträucherdef.
  ! Im Nachtrag (7) 81 dark gray auf ( 89) braun 
  IF (ASSOCIATED(Cell)) THEN
   IF (Cell%LandClass /=0 ) THEN
    IF (Cell%vc>0) THEN
        ! polygons 
        mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
        IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
        nvert=Cell%vc
        nRec=nRec+1
        WRITE(10,REC=nRec) mat
        nRec=nRec+1
        WRITE(10,REC=nRec) nvert
        DO li=1,nvert
          nRec=nRec+1
          !tempw=VertOut(Cell%VertCut(li))%Point%x
          tempw=xParametricOut(VertOut(Cell%VertCut(li)))
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          !tempw=VertOut(Cell%VertCut(li))%Point%y
          tempw=yParametricOut(VertOut(Cell%VertCut(li)))
          WRITE(10,REC=nRec) tempw
        END DO
        DO li=1,nvert
          nRec=nRec+1
          !tempw=VertOut(Cell%VertCut(li))%Point%z
          tempw=zParametricOut(VertOut(Cell%VertCut(li)))
          WRITE(10,REC=nRec) tempw
        END DO
        IF(k==1) THEN
          IF(Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
             Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) THEN
             ! polygon: (Face1 mit Edge%ec=1 und Grundflaeche-Gelaende)
             !--------- Specialfall 2Hügel, 
             !? ist vc>0, oder vc=0 aktuell gesetzt
             mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
             IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
             nvert=Cell%Face1%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) mat
             nRec=nRec+1
             WRITE(10,REC=nRec) nvert
             DO li=1,nvert
               nRec=nRec+1
               !tempw=VertOut(Cell%VertCut(li))%Point%x
               tempw=xParametricOut(VertOut(Cell%Face1%VertexList(li)))
               WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
               nRec=nRec+1
               !tempw=VertOut(Cell%VertCut(li))%Point%y
               tempw=yParametricOut(VertOut(Cell%Face1%VertexList(li)))
               WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
               nRec=nRec+1
               !tempw=VertOut(Cell%VertCut(li))%Point%z
               tempw=zParametricOut(VertOut(Cell%Face1%VertexList(li)))
               WRITE(10,REC=nRec) tempw
             END DO
          END IF
        END IF

    ELSE IF  & !(Cell%in_out==4) THEN !Cell%Vol==maxVoloder >0.0d0 ausreichend
      (k==1) THEN
      IF (Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
          Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
             ! polygon: (Face1)
             !---------
         mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
         IF(Cell%LandClass/=0) mat=Lk(Cell%LandClass)
         nvert=4
         nRec=nRec+1
         WRITE(10,REC=nRec) mat
         nRec=nRec+1
         WRITE(10,REC=nRec) nvert
         DO li=1,nvert
           SELECT CASE (li)
                CASE (1)        !.. Vertices(i-1,j-1,k)
                    ix=i-1;iy=j-1
                CASE (2)        !.. Vertices(i,j-1,k
                    ix=i
                CASE (3)        !..Vertices(i,j,k)
                    iy=j
                CASE (4)        !..Vertices(i-1,j,k)
                    ix=i-1
               END SELECT
           nRec=nRec+1
           tempw=xParametricOut(Vertices(ix,iy,k-1))
           WRITE(10,REC=nRec) tempw
         END DO
         DO li=1,nvert
           SELECT CASE (li)
                CASE (1)        !.. Vertices(i-1,j-1,k)
                    ix=i-1;iy=j-1
                CASE (2)        !.. Vertices(i,j-1,k
                    ix=i
                CASE (3)        !..Vertices(i,j,k)
                    iy=j
                CASE (4)        !..Vertices(i-1,j,k)
                    ix=i-1
               END SELECT
           nRec=nRec+1
           tempw=yParametricOut(Vertices(ix,iy,k-1))
           WRITE(10,REC=nRec) tempw
         END DO
         DO li=1,nvert
           SELECT CASE (li)
                CASE (1)        !.. Vertices(i-1,j-1,k)
                    ix=i-1;iy=j-1
                CASE (2)        !.. Vertices(i,j-1,k
                    ix=i
                CASE (3)        !..Vertices(i,j,k)
                    iy=j
                CASE (4)        !..Vertices(i-1,j,k)
                    ix=i-1
               END SELECT
           nRec=nRec+1
           tempw=zParametricOut(Vertices(ix,iy,k-1))
           WRITE(10,REC=nRec) tempw
         END DO
      END IF   ! IF Face1 Grenze ((V0-V4)%in_out==0)
    END IF   ! IF(Cell%vc>0) ELSE IF (k==1)
   END IF  !  IF (Cell%LandClass/=0)
  END IF  ! IF (ASSOCIATED(Cell)) 
END SUBROUTINE WriteLandKlPlanePolyGMVBinary


SUBROUTINE WriteCellCutPlaneGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out

  IF (ASSOCIATED(Cell)) THEN
    IF (Cell%vc>0) THEN
        ! general
        nFaces=1
        nVerts(nFaces)=Cell%vc

        !.. WRITE(10) general,nfaces
        nRec=nRec+1
        WRITE(10,REC=nRec) general(1:4)
        nRec=nRec+1
        WRITE(10,REC=nRec) general(5:8)
        nRec=nRec+1
        WRITE(10,REC=nRec) nfaces
        !.. WRITE(10) (nVerts(i),i=1,nFaces)
        nRec=nRec+1
        WRITE(10,REC=nRec) nVerts(1)
        DO iVert=1,Cell%vc
          nRec=nRec+1
          WRITE(10,REC=nRec) Cell%VertCut(iVert)
        END DO

        IF(k==1) THEN
          IF(Cell%Face1%Edge1%yes_sp==1.OR.Cell%Face1%Edge2%yes_sp==1.OR. &
             Cell%Face1%Edge3%yes_sp==1.OR.Cell%Face1%Edge4%yes_sp==1) THEN
             ! polygon: (Face1,Grundflaeche-Gelaende)
             nFaces=1
             nVerts(nFaces)=Cell%Face1%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) general(1:4)
             nRec=nRec+1
             WRITE(10,REC=nRec) general(5:8)
             nRec=nRec+1
             WRITE(10,REC=nRec) nfaces
             !.. WRITE(10) (nVerts(i),i=1,nFaces)
             nRec=nRec+1
             WRITE(10,REC=nRec) nVerts(1)
             DO iVert=1,nVerts(1)
               nRec=nRec+1
               WRITE(10,REC=nRec) Cell%Face1%VertexList(iVert)  
             END DO
          END IF
        END IF
    ELSE
      IF (k==1) THEN
        IF(Cell%Face1%Edge1%Vert1%in_out==0.AND.Cell%Face1%Edge1%Vert2%in_out==0.AND. &
           Cell%Face1%Edge3%Vert1%in_out==0.AND.Cell%Face1%Edge3%Vert2%in_out==0)THEN
         !Face1, Grenzflaeche-Gelaende
         nFaces=1
         !.. WRITE(10) general,nfaces
         nRec=nRec+1
         WRITE(10,REC=nRec) general(1:4)
         nRec=nRec+1
         WRITE(10,REC=nRec) general(5:8)
         nRec=nRec+1
         WRITE(10,REC=nRec) nfaces
         nRec=nRec+1
         WRITE(10,REC=nRec) 4   ! Punkte je Face
         DO iVert=1,4
           nRec=nRec+1
           WRITE(10,REC=nRec) Cell%Face1%VertexList(ivert)
         END DO
        END IF
      END IF ! k==1 
    END IF  ! IF (Cell%vc>0) ELSE
  END IF   ! IF (ASSOCIATED(Cell))
END SUBROUTINE WriteCellCutPlaneGMVBinary


SUBROUTINE WriteCellOroGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: nFaces,nVertices
  INTEGER :: nVerts(7),yFace(6)
  INTEGER :: ListVert(8),iList
  INTEGER :: iVert,in_out
 
  yFace(1:6)=0 
  IF (ASSOCIATED(Cell)) THEN
      IF (Cell%vc>0) THEN
         ! general
         !--------
     !WRITE(10,*)
     !WRITE(10,*) "Celle(i,j,k)= ",i,j,k, "                   Cell%vc=", Cell%vc
         nr_gen=nr_gen+1
         nr_cell_vc_oro=nr_cell_vc_oro+1
         nFaces=1
         IF ((MIN(Cell%Face1%Edge1%Vert1%in_out,Cell%Face1%Edge1%Vert2%in_out, &
                 Cell%Face1%Edge3%Vert1%in_out,Cell%Face1%Edge3%Vert2%in_out)<=0) &
             .AND.Cell%Face1%NumberVert>2) THEN
           nVerts(nFaces)=Cell%Face1%NumberVert
           nFaces=nFaces+1
           yFace(1)=1
         END IF
         IF ((MIN(Cell%Face2%Edge1%Vert1%in_out,Cell%Face2%Edge1%Vert2%in_out, &
                  Cell%Face2%Edge3%Vert1%in_out,Cell%Face2%Edge3%Vert2%in_out)<=0) &
             .AND.Cell%Face2%NumberVert>2) THEN
           nVerts(nFaces)=Cell%Face2%NumberVert
           nFaces=nFaces+1
           yFace(2)=1
         END IF
         IF ((MIN(Cell%Face3%Edge1%Vert1%in_out,Cell%Face3%Edge1%Vert2%in_out, &
                  Cell%Face3%Edge3%Vert1%in_out,Cell%Face3%Edge3%Vert2%in_out)<=0) &
             .AND.Cell%Face3%NumberVert>2) THEN
           nVerts(nFaces)=Cell%Face3%NumberVert
           nFaces=nFaces+1
           yFace(3)=1
         END IF
         IF ((MIN(Cell%Face4%Edge1%Vert1%in_out,Cell%Face4%Edge1%Vert2%in_out, &
                  Cell%Face4%Edge3%Vert1%in_out,Cell%Face4%Edge3%Vert2%in_out)<=0) &
             .AND.Cell%Face4%NumberVert>2) THEN
           nVerts(nFaces)=Cell%Face4%NumberVert
           nFaces=nFaces+1
           yFace(4)=1
         END IF
         IF ((MIN(Cell%Face5%Edge1%Vert1%in_out,Cell%Face5%Edge1%Vert2%in_out, &
                  Cell%Face5%Edge3%Vert1%in_out,Cell%Face5%Edge3%Vert2%in_out)<=0) &
             .AND.Cell%Face5%NumberVert>2) THEN
           nVerts(nFaces)=Cell%Face5%NumberVert
           nFaces=nFaces+1
           yFace(5)=1
         END IF
         IF ((MIN(Cell%Face6%Edge1%Vert1%in_out,Cell%Face6%Edge1%Vert2%in_out, &
                  Cell%Face6%Edge3%Vert1%in_out,Cell%Face6%Edge3%Vert2%in_out)<=0) &
             .AND.Cell%Face6%NumberVert>2) THEN
           nVerts(nFaces)=Cell%Face6%NumberVert
           nFaces=nFaces+1
           yFace(6)=1
         END IF
         nVerts(nFaces)=Cell%vc
         IF (Cell%Face1%Edge1%Vert1%in_out==0.AND. &
             Cell%Face1%Edge1%Vert2%in_out==0.AND. &
             Cell%Face1%Edge3%Vert1%in_out==0.AND. &
             Cell%Face1%Edge3%Vert2%in_out==0 ) THEN  
             ! Speziell Celle-Face1 Grenze Berg ist
               nfaces=nfaces-1
         END IF

         !WRITE(10,'(a8,I8)') general,nfaces
         nRec=nRec+1
         WRITE(10,REC=nRec) general(1:4)
         nRec=nRec+1
         WRITE(10,REC=nRec) general(5:8)
         nRec=nRec+1
         WRITE(10,REC=nRec) nfaces
         !WRITE(10,*) (nVerts(iVert),iVert=1,nFaces)
         DO iVert=1,nFaces
           nRec=nRec+1
           WRITE(10,REC=nRec) nVerts(iVert)
         END DO
         IF (yFace(1)==1) THEN
           !WRITE(10,*) Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
           DO iVert=1,Cell%Face1%NumberVert
              nRec=nRec+1
              WRITE(10,REC=nRec)  Cell%Face1%VertexList(iVert)
           END DO
         END IF
         IF (yFace(2)==1) THEN
           !WRITE(10,*) Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
           DO iVert=1,Cell%Face2%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face2%VertexList(iVert)
           END DO
         END IF
         IF (yFace(3)==1) THEN
           !WRITE(10,*) Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
           DO iVert=1,Cell%Face3%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face3%VertexList(iVert)
           END DO
         END IF
         IF (yFace(4)==1) THEN
           !WRITE(10,*) Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
           DO iVert=1,Cell%Face4%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face4%VertexList(iVert)
           END DO
         END IF
         IF (yFace(5)==1) THEN
           !WRITE(10,*) Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
           DO iVert=1,Cell%Face5%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face5%VertexList(iVert)
           END DO
         END IF
         IF (yFace(6)==1) THEN
           !WRITE(10,*) Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
           DO iVert=1,Cell%Face6%NumberVert
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%Face6%VertexList(iVert)
           END DO
         END IF
         !WRITE(10,*) Cell%VertCut(1:Cell%vc)
         IF (nfaces/=1) THEN
           DO iVert=1,Cell%vc
             nRec=nRec+1
             WRITE(10,REC=nRec) Cell%VertCut(iVert)
           END DO
         END IF
      
      ELSE  ! Allokiert und vc==0

        IF (Cell%Face1%in_out==0.AND.(Cell%Face2%in_out>=2.AND.Cell%Face2%in_out<=4)) THEN
          ! Schnittflaechen: Face1 komplett
          ! Face2%in_out=[2..4] 
          ! wird in SortVertCutCell vc=0 gesetzt, nur Face1 ausgeben! 
       !WRITE(10,*)
       !WRITE(10,*) "Celle(i,j,k)= ",i,j,k, "          Face1%in_out==0.AND.Face2(in_out=2..4"
          ! general
          !--------
          IF(k==Domain%iz0+1) THEN
            nr_grenzeF1=nr_grenzeF1+1
            nr_cell_novc_oro=nr_cell_novc_oro+1
            nRec=nRec+1
            WRITE(10,REC=nRec) general(1:4)
            nRec=nRec+1
            WRITE(10,REC=nRec) general(5:8)
            nRec=nRec+1
            WRITE(10,REC=nRec) 1 !Number Faces
            nRec=nRec+1
            WRITE(10,REC=nRec) Cell%Face1%NumberVert  
            DO iVert=1,Cell%Face1%NumberVert
              nRec=nRec+1
              WRITE(10,REC=nRec) Cell%Face1%VertexList(iVert)
            END DO
          END IF

        ELSE IF(MAX(Cell%Face2%Edge1%Vert1%in_out,Cell%Face2%Edge1%Vert2%in_out, &
                Cell%Face2%Edge3%Vert1%in_out,Cell%Face2%Edge3%Vert2%in_out)==0) THEN
          ! Schnittflaechen: Face2 .OR. 1Edge .OR. 1Vertex
          ! wird in SortVertCutCell vc=0 gesetzt, Celle als hex ausgeben
          ! Face1%in_out=[-1 bis -4], Face2%in_out=[-3 bis 0]  
          ! hex
          !-----
          !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def
          !         /|       /|
          !        1--------2 |
          !        | 8------|-7
          !        |/       |/
          !        5--------6
           nr_grenzeF2=nr_grenzeF2+1
           nr_cell_novc_oro=nr_cell_novc_oro+1
       !WRITE(10,*)
       !WRITE(10,*) "Celle(i,j,k)= ",i,j,k,"         MAX(Cell%Face2%Edge...%Vert..%in_out)==0" 

           nRec=nRec+1
           WRITE(10,REC=nRec) hex(1:4)
           nRec=nRec+1
           WRITE(10,REC=nRec) hex(5:8)
           nRec=nRec+1
           WRITE(10,REC=nRec) 8
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i-1,j-1,k)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i,j-1,k)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i,j,k)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i-1,j,k)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i-1,j-1,k-1)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i,j-1,k-1)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i,j,k-1)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i-1,j,k-1)%nrInP
        ELSE IF (Cell%Face2%in_out==-4 .AND. &
                 Cell%Face1%in_out==-1 .OR. Cell%Face1%in_out==-2) THEN
          nr_grenzeSp3=nr_grenzeSp3+1
          nr_cell_novc_oro=nr_cell_novc_oro+1
       !WRITE(10,*)
       !WRITE(10,*) "Celle(i,j,k)= ",i,j,k,"     no vc,  Special3-(Cell%Face2%in_out)==-4 und"
       !WRITE(10,*) "                                              Cell%Face1%in_out==-1 or. -2"  

           nRec=nRec+1
           WRITE(10,REC=nRec) hex(1:4)
           nRec=nRec+1
           WRITE(10,REC=nRec) hex(5:8)
           nRec=nRec+1
           WRITE(10,REC=nRec) 8
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i-1,j-1,k)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i,j-1,k)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i,j,k)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i-1,j,k)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i-1,j-1,k-1)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i,j-1,k-1)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i,j,k-1)%nrInP
           nRec=nRec+1
           WRITE(10,REC=nRec) Vertices(i-1,j,k-1)%nrInP

        END IF  ! Spezialfälle, Celle allokiert, vc=0
      END IF  ! if( vc>0)  else

  ELSE !!! .NOT. ASSOCIATED(Cell) 
    in_out=Vertices(i-1,j-1,k-1)%in_out &
          +Vertices(i,j-1,k-1)%in_out &
          +Vertices(i-1,j,k-1)%in_out &
          +Vertices(i-1,j-1,k)%in_out &
          +Vertices(i-1,j,k)%in_out &
          +Vertices(i,j-1,k)%in_out &
          +Vertices(i,j,k-1)%in_out &
          +Vertices(i,j,k)%in_out
    IF (in_out==-8) THEN
      ! hex
      !------
      !          4--------3    !Point-Folge Output-Hex-GMV-Vorlage-Def
      !         /|       /|
      !        1--------2 |
      !        | 8------|-7
      !        |/       |/
      !        5--------6
     !WRITE(10,*)
     !WRITE(10,*) "i,j,k= ",i,j,k,"               .NOT. ASSOCIATED(Cell), Cell%in_out==-8"
      nr_insidehex=nr_insidehex+1
      nRec=nRec+1
      WRITE(10,REC=nRec) hex(1:4)
      nRec=nRec+1
      WRITE(10,REC=nRec) hex(5:8)
      nRec=nRec+1
      WRITE(10,REC=nRec) 8
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j-1,k)%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j-1,k)%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j,k)%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j,k)%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j-1,k-1)%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j-1,k-1)%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i,j,k-1)%nrInP
      nRec=nRec+1
      WRITE(10,REC=nRec) Vertices(i-1,j,k-1)%nrInP
    END IF
  END IF
END SUBROUTINE WriteCellOroGMVBinary

SUBROUTINE WriteOroCutPlanePolyGMVBinary(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  ! local
  INTEGER :: mat,nvert
  INTEGER :: li

  mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
  IF (ASSOCIATED(Cell)) THEN
      IF (Cell%vc>0) THEN
         ! polygon
         !---------
         mat=25    ! 14=grau,16=lila-braun,25=rot-braun,33=gruen
         nvert=Cell%vc
         nRec=nRec+1
         WRITE(10,REC=nRec) mat
         nRec=nRec+1
         WRITE(10,REC=nRec) nvert
         DO li=1,nvert
           nRec=nRec+1
           !tempw=VertIn(Cell%VertCut(li))%Point%x
           tempw=xParametricOut(VertIn(Cell%VertCut(li)))
           WRITE(10,REC=nRec) tempw
         END DO
         DO li=1,nvert
           nRec=nRec+1
           !tempw=VertIn(Cell%VertCut(li))%Point%y
           tempw=yParametricOut(VertIn(Cell%VertCut(li)))
           WRITE(10,REC=nRec) tempw
         END DO
         DO li=1,nvert
           nRec=nRec+1
           !tempw=VertIn(Cell%VertCut(li))%Point%z
           tempw=zParametricOut(VertIn(Cell%VertCut(li)))
           WRITE(10,REC=nRec) tempw
         END DO
      ELSE IF (Cell%vc==0) THEN
          IF(Cell%Face1%Edge1%Vert1%in_out==0.AND. &
             Cell%Face1%Edge1%Vert2%in_out==0.AND. &
             Cell%Face1%Edge3%Vert1%in_out==0.AND. &
             Cell%Face1%Edge3%Vert2%in_out==0 ) THEN
             mat=25   !mat=26 !zum check
             nvert=Cell%Face1%NumberVert !4
             nRec=nRec+1
             WRITE(10,REC=nRec) mat
             nRec=nRec+1
             WRITE(10,REC=nRec) nvert
             DO li=1,nvert
               nRec=nRec+1
               tempw=xParametricOut(VertIn(Cell%Face1%VertexList(li)))
               WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
               nRec=nRec+1
               tempw=yParametricOut(VertIn(Cell%Face1%VertexList(li)))
               WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
               nRec=nRec+1
               tempw=zParametricOut(VertIn(Cell%Face1%VertexList(li)))
               WRITE(10,REC=nRec) tempw
             END DO
          END IF
          !IF (Cell%Face2%in_out==0) THEN
          IF(Cell%Face2%Edge1%Vert1%in_out==0.AND. &
             Cell%Face2%Edge1%Vert2%in_out==0.AND. &
             Cell%Face2%Edge3%Vert1%in_out==0.AND. &
             Cell%Face2%Edge3%Vert2%in_out==0 ) THEN
             mat=25 !mat=27 ! zum check
             nvert=Cell%Face2%NumberVert !4
             nRec=nRec+1
             WRITE(10,REC=nRec) mat
             nRec=nRec+1
             WRITE(10,REC=nRec) nvert
             DO li=1,nvert
                nRec=nRec+1
                tempw=xParametricOut(VertIn(Cell%Face2%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=yParametricOut(VertIn(Cell%Face2%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=zParametricOut(VertIn(Cell%Face2%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
          END IF
          !IF (Cell%Face3%in_out==0) THEN
          IF(Cell%Face3%Edge1%Vert1%in_out==0.AND. &
             Cell%Face3%Edge1%Vert2%in_out==0.AND. &
             Cell%Face3%Edge3%Vert1%in_out==0.AND. &
             Cell%Face3%Edge3%Vert2%in_out==0 ) THEN
             !mat=28 !mat=28 ! zum check
             nvert=Cell%Face3%NumberVert !4
             nRec=nRec+1
             WRITE(10,REC=nRec) mat
             nRec=nRec+1
             WRITE(10,REC=nRec) nvert
             DO li=1,nvert
               nRec=nRec+1
               tempw=xParametricOut(VertIn(Cell%Face3%VertexList(li)))
               WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=yParametricOut(VertIn(Cell%Face3%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=zParametricOut(VertIn(Cell%Face3%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
          END IF
          !IF (Cell%Face4%in_out==0) THEN
          IF(Cell%Face4%Edge1%Vert1%in_out==0.AND. &
             Cell%Face4%Edge1%Vert2%in_out==0.AND. &
             Cell%Face4%Edge3%Vert1%in_out==0.AND. &
             Cell%Face4%Edge3%Vert2%in_out==0 ) THEN
             !mat=29 !mat=29 ! zum check
             nvert=Cell%Face4%NumberVert !4
             nRec=nRec+1
             WRITE(10,REC=nRec) mat
             nRec=nRec+1
             WRITE(10,REC=nRec) nvert
             DO li=1,nvert
                nRec=nRec+1
                tempw=xParametricOut(VertIn(Cell%Face4%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
               nRec=nRec+1
               tempw=yParametricOut(VertIn(Cell%Face4%VertexList(li)))
               WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=zParametricOut(VertIn(Cell%Face4%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
          END IF
          !IF (Cell%Face5%in_out==0) THEN
          IF(Cell%Face5%Edge1%Vert1%in_out==0.AND. &
             Cell%Face5%Edge1%Vert2%in_out==0.AND. &
             Cell%Face5%Edge3%Vert1%in_out==0.AND. &
             Cell%Face5%Edge3%Vert2%in_out==0 ) THEN
             !mat=30 !mat=30 ! zum check
             nvert=Cell%Face5%NumberVert !4
             nRec=nRec+1
             WRITE(10,REC=nRec) mat
             nRec=nRec+1
             WRITE(10,REC=nRec) nvert
             DO li=1,nvert
                nRec=nRec+1
                tempw=xParametricOut(VertIn(Cell%Face5%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=yParametricOut(VertIn(Cell%Face5%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=zParametricOut(VertIn(Cell%Face5%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
          END IF
          !IF (Cell%Face6%in_out==0) THEN
          IF(Cell%Face6%Edge1%Vert1%in_out==0.AND. &
             Cell%Face6%Edge1%Vert2%in_out==0.AND. &
             Cell%Face6%Edge3%Vert1%in_out==0.AND. &
             Cell%Face6%Edge3%Vert2%in_out==0 ) THEN
             !mat=31 !mat=31 ! zum check
             nvert=Cell%Face6%NumberVert !4
             nRec=nRec+1
             WRITE(10,REC=nRec) mat
             nRec=nRec+1
             WRITE(10,REC=nRec) nvert
             DO li=1,nvert
                nRec=nRec+1
                tempw=xParametricOut(VertIn(Cell%Face6%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=yParametricOut(VertIn(Cell%Face6%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
             DO li=1,nvert
                nRec=nRec+1
                tempw=zParametricOut(VertIn(Cell%Face6%VertexList(li)))
                WRITE(10,REC=nRec) tempw
             END DO
          END IF

      END IF  ! if(vc) else vc==0
  END IF   ! allocated
END SUBROUTINE WriteOroCutPlanePolyGMVBinary


SUBROUTINE WriteBlockGMVAscii
  INTEGER :: i,j,k,schritt
  INTEGER :: mat,nvert

  ! polygons
  !WRITE(10,'(a8)') polygons
  mat=6
  nvert=4
  schritt=iy1-iy0
  DO j=iy0,iy1,schritt
    Write(10,*) mat, nvert  &
      ,xParametricOut(Vertices(ix0,j,iz0))   &
      ,xParametricOut(Vertices(ix1,j,iz0))   &
      ,xParametricOut(Vertices(ix1,j,iz1))   &
      ,xParametricOut(Vertices(ix0,j,iz1))   &
      ,yParametricOut(Vertices(ix0,j,iz0))   &
      ,yParametricOut(Vertices(ix1,j,iz0))   &
      ,yParametricOut(Vertices(ix1,j,iz1))   &
      ,yParametricOut(Vertices(ix0,j,iz1))   &
      ,zParametricOut(Vertices(ix0,j,iz0))   &
      ,zParametricOut(Vertices(ix1,j,iz0))   &
      ,zParametricOut(Vertices(ix1,j,iz1))   &
      ,zParametricOut(Vertices(ix0,j,iz1))
  END DO
  schritt=iz1-iz0
  DO k=iz0,iz1,schritt
    Write(10,*)  mat, nvert &
      ,xParametricOut(Vertices(ix0,iy0,k))   &
      ,xParametricOut(Vertices(ix1,iy0,k))   &
      ,xParametricOut(Vertices(ix1,iy1,k))   &
      ,xParametricOut(Vertices(ix0,iy1,k))   &
      ,yParametricOut(Vertices(ix0,iy0,k))   &
      ,yParametricOut(Vertices(ix1,iy0,k))   &
      ,yParametricOut(Vertices(ix1,iy1,k))   &
      ,yParametricOut(Vertices(ix0,iy1,k))   &
      ,zParametricOut(Vertices(ix0,iy0,k))   &
      ,zParametricOut(Vertices(ix1,iy0,k))   &
      ,zParametricOut(Vertices(ix1,iy1,k))   &
      ,zParametricOut(Vertices(ix0,iy1,k))
  END DO
  schritt=ix1-ix0
  DO i=ix0,ix1,schritt
    Write(10,*)  mat, nvert  &
      ,xParametricOut(Vertices(i,iy0,iz0))   &
      ,xParametricOut(Vertices(i,iy0,iz1))   &
      ,xParametricOut(Vertices(i,iy1,iz1))   &
      ,xParametricOut(Vertices(i,iy1,iz0))   &
      ,yParametricOut(Vertices(i,iy0,iz0))   &
      ,yParametricOut(Vertices(i,iy0,iz1))   &
      ,yParametricOut(Vertices(i,iy1,iz1))   &
      ,yParametricOut(Vertices(i,iy1,iz0))   &
      ,zParametricOut(Vertices(i,iy0,iz0))   &
      ,zParametricOut(Vertices(i,iy0,iz1))   &
      ,zParametricOut(Vertices(i,iy1,iz1))   &
      ,zParametricOut(Vertices(i,iy1,iz0))
  END DO

  !WRITE(10,'(a8)') endpoly

END SUBROUTINE WriteBlockGMVAscii


SUBROUTINE WriteBlockGMVBinary
  INTEGER :: i,j,k,schritt
  INTEGER :: mat,nvert

  ! polygons
  mat=6
  nvert=4
  schritt=iy1-iy0
  DO j=iy0,iy1,schritt
    nRec=nRec+1
    WRITE(10,REC=nRec) mat
    nRec=nRec+1
    WRITE(10,REC=nRec) nvert
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix0,j,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix1,j,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix1,j,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix0,j,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix0,j,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix1,j,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix1,j,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix0,j,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix0,j,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix1,j,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix1,j,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix0,j,iz1))
    WRITE(10,REC=nRec) tempw
  END DO

  schritt=iz1-iz0
  DO k=iz0,iz1,schritt
    nRec=nRec+1
    WRITE(10,REC=nRec) mat
    nRec=nRec+1
    WRITE(10,REC=nRec) nvert
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix0,iy0,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix1,iy0,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix1,iy1,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(ix0,iy1,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix0,iy0,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix1,iy0,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix1,iy1,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(ix0,iy1,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix0,iy0,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix1,iy0,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix1,iy1,k))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(ix0,iy1,k))
    WRITE(10,REC=nRec) tempw
  END DO

  schritt=ix1-ix0
  DO i=ix0,ix1,schritt
    nRec=nRec+1
    WRITE(10,REC=nRec) mat
    nRec=nRec+1
    WRITE(10,REC=nRec) nvert
    nRec=nRec+1
    tempw=xParametricOut(Vertices(i,iy0,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(i,iy0,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(i,iy1,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=xParametricOut(Vertices(i,iy1,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(i,iy0,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(i,iy0,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(i,iy1,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=yParametricOut(Vertices(i,iy1,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(i,iy0,iz0))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(i,iy0,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(i,iy1,iz1))
    WRITE(10,REC=nRec) tempw
    nRec=nRec+1
    tempw=zParametricOut(Vertices(i,iy1,iz0))
    WRITE(10,REC=nRec) tempw
  END DO

END SUBROUTINE WriteBlockGMVBinary


SUBROUTINE WriteBlockPartGMVAscii
  INTEGER :: i,j,k,schritt
  INTEGER :: mat,nvert

  ! polygons
  !WRITE(10,'(a8)') polygons
  mat=6
  nvert=4
  !block_plane+part
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  !...n,s
  IF(v_x0==ix0.AND.v_x1==ix1.AND.v_z0==iz0.AND.v_z1==iz1) THEN
      IF(v_y0==iy0) THEN   !'s'
        CALL WritePolyFaceXZ(mat,v_y0)
      ELSE IF(v_y1==iy1) THEN   ! 'n'
        CALL WritePolyFaceXZ(mat,v_y1)
      END IF
      CALL WritePolyEdgeY(mat,v_x0,v_z0)
      CALL WritePolyEdgeY(mat,v_x1,v_z0)
      CALL WritePolyEdgeY(mat,v_x1,v_z1)
      CALL WritePolyEdgeY(mat,v_x0,v_z1)
      ! jeweils LineY '---' verlängern
      CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
      CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
      CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
      CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
  !...e,w
  ELSE IF(v_y0==iy0.AND.v_y1==iy1.AND.v_z0==iz0.AND.v_z1==iz1) THEN   
      IF(v_x0==ix0) THEN   ! 'e'
        CALL WritePolyFaceYZ(mat,v_x0)
      ELSE IF(v_x1==ix1) THEN  !'w'
        CALL WritePolyFaceYZ(mat,v_x1)
      END IF
      CALL WritePolyEdgeX(mat,v_y0,v_z0)
      CALL WritePolyEdgeX(mat,v_y0,v_z1)
      CALL WritePolyEdgeX(mat,v_y1,v_z1)
      CALL WritePolyEdgeX(mat,v_y1,v_z0)
      ! jeweils LineX '---' verlängern
      CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
      CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
      CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
      CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
  !...t,b
  ELSE IF(v_x0==ix0.AND.v_x1==ix1.AND.v_y0==iy0.AND.v_y1==iy1) THEN
      IF(v_z0==iz0) THEN  ! 'b'
        CALL WritePolyFaceXY(mat,v_z0)
      ELSE IF(v_z1==iz1) THEN  !'t'
        CALL WritePolyFaceXY(mat,v_z1)
      END IF
      CALL WritePolyEdgeZ(mat,v_x0,v_y0)
      CALL WritePolyEdgeZ(mat,v_x1,v_y0)
      CALL WritePolyEdgeZ(mat,v_x1,v_y1)
      CALL WritePolyEdgeZ(mat,v_x0,v_y1)
      ! jeweils LineZ '---' verlängern
      CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
      CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
      CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
      CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
  ELSE  !block_edge+part
  !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     !......................................................
     IF(v_x0==ix0.AND.v_x1==ix1) THEN
        IF(v_z0==iz0) THEN
            IF(v_y0==iy0) THEN ! 'b,s'
              CALL WritePolyEdgeX(mat,v_y0,v_z0)
              CALL WritePolyEdgeY(mat,v_x0,v_z0)
              CALL WritePolyEdgeY(mat,v_x1,v_z0)
              CALL WritePolyEdgeZ(mat,v_x0,v_y0)
              CALL WritePolyEdgeZ(mat,v_x1,v_y0)
              !... LineY/-Z '---' verlängern
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
            ELSE IF(v_y1==iy1) THEN ! 'b,n'
              CALL WritePolyEdgeX(mat,v_y1,v_z0)
              CALL WritePolyEdgeY(mat,v_x0,v_z0)
              CALL WritePolyEdgeY(mat,v_x1,v_z0)
              CALL WritePolyEdgeZ(mat,v_x0,v_y1)
              CALL WritePolyEdgeZ(mat,v_x1,v_y1)
              !... LineY/-Z '---' verlängern
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
            ELSE  ! (v_x0==ix0.AND.v_x1==ix1)
                  ! (v_z0==iz0)  .and. (v_y0/=iy0 .and. v_y1/=iy1)
              !botton-direction
              CALL WritePolyEdgeY(mat,v_x0,v_z0)
              CALL WritePolyEdgeY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              !top-direction, inside block
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
            END IF
        ELSE IF(v_z1==iz1) THEN  
            IF(v_y0==iy0) THEN   ! 't,s'
              CALL WritePolyEdgeX(mat,v_y0,v_z1)
              CALL WritePolyEdgeY(mat,v_x0,v_z1)
              CALL WritePolyEdgeY(mat,v_x1,v_z1)
              CALL WritePolyEdgeZ(mat,v_x0,v_y0)
              CALL WritePolyEdgeZ(mat,v_x1,v_y0)
              !... LineY/-Z '---' verlängern
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
            ELSE IF(v_y1==iy1) THEN  ! 't,n'
              CALL WritePolyEdgeX(mat,v_y1,v_z1)
              CALL WritePolyEdgeY(mat,v_x0,v_z1)
              CALL WritePolyEdgeY(mat,v_x1,v_z1)
              CALL WritePolyEdgeZ(mat,v_x0,v_y1)
              CALL WritePolyEdgeZ(mat,v_x1,v_y1)
              !... LineY/-Z   '---' verlängern
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
            ELSE  ! (v_x0==ix0.AND.v_x1==ix1)
                  ! (v_z1==iz1)  .and. (v_y0/=iy0 .and. v_y1/=iy1)
              !top-direction
              CALL WritePolyEdgeY(mat,v_x0,v_z1)
              CALL WritePolyEdgeY(mat,v_x1,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
              !botton direction, inside block
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
            END IF
        ELSE  ! (v_x0==ix0 .AND. v_x1==ix1)
              ! (v_z0/=iz0 .and. v_z1/=iz1)
            IF(v_y0==iy0) THEN  
              !south-direction
              CALL WritePolyEdgeZ(mat,v_x0,v_y0)
              CALL WritePolyEdgeZ(mat,v_x1,v_y0)
              !... Line-Z '---' verlängern
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              !north-direction, inside block
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
            ELSE IF(v_y1==iy1) THEN 
              !north-direction
              CALL WritePolyEdgeZ(mat,v_x0,v_y1)
              CALL WritePolyEdgeZ(mat,v_x1,v_y1)
              !... Line-Z '---' verlängern
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              !south-direction-inside-block
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
            ELSE  ! (v_y0/=iy0 .and. v_y1/=iy1)
              !north-south-top-botton, inside block, 
              !x0-x1-flanking
              !... Line-Z/-Y  '---'  verlängern 
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
            END IF
        END IF
     !......................................................
     ELSE IF(v_y0==iy0.AND.v_y1==iy1) THEN
        IF(v_x0==ix0) THEN
           IF(v_z0==iz0) THEN      !'w,b'
              CALL WritePolyEdgeY(mat,v_x0,v_z0)
              CALL WritePolyEdgeZ(mat,v_x0,v_y0)
              CALL WritePolyEdgeZ(mat,v_x0,v_y1)
              CALL WritePolyEdgeX(mat,v_y0,v_z0)
              CALL WritePolyEdgeX(mat,v_y1,v_z0)
              !... LineZ/-X  '---'  verlängern
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
           ELSE IF(v_z1==iz1) THEN  !'w,t'
              CALL WritePolyEdgeY(mat,v_x0,v_z1)
              CALL WritePolyEdgeZ(mat,v_x0,v_y0)
              CALL WritePolyEdgeZ(mat,v_x0,v_y1)
              CALL WritePolyEdgeX(mat,v_y0,v_z1)
              CALL WritePolyEdgeX(mat,v_y1,v_z1)
              !... LineZ/-X  '---'  verlängern
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
           ELSE ! (v_y0==iy0.AND.v_y1==iy1) 
                ! (v_x0==ix0)  .and. (v_z0/=iz0 .and. v_z1/=iz1)
              !west-direction
              CALL WritePolyEdgeZ(mat,v_x0,v_y0)
              CALL WritePolyEdgeZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              !east-direction, inside block
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
           END IF
        ELSE IF(v_x1==ix1) THEN
           IF(v_z0==iz0) THEN      !'e,b'
              CALL WritePolyEdgeY(mat,v_x1,v_z0)
              CALL WritePolyEdgeZ(mat,v_x1,v_y0)
              CALL WritePolyEdgeZ(mat,v_x1,v_y1)
              CALL WritePolyEdgeX(mat,v_y0,v_z0)
              CALL WritePolyEdgeX(mat,v_y1,v_z0)
              !... LineZ/-X  '---'  verlängern
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
           ELSE IF(v_z1==iz1) THEN  !'e,t'
              CALL WritePolyEdgeY(mat,v_x1,v_z1)
              CALL WritePolyEdgeZ(mat,v_x1,v_y0)
              CALL WritePolyEdgeZ(mat,v_x1,v_y1)
              CALL WritePolyEdgeX(mat,v_y0,v_z1)
              CALL WritePolyEdgeX(mat,v_y1,v_z1)
              !... LineZ/-X  '---'  verlängern
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
           ELSE ! (v_y0==iy0.AND.v_y1==iy1) 
                ! (v_x1==ix1)  .and. (v_z0/=iz0 .and. v_z1/=iz1)
              !east-direction
              CALL WritePolyEdgeZ(mat,v_x1,v_y0)
              CALL WritePolyEdgeZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              !west-direction, inside block
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
           END IF
        ELSE  ! (v_y0==iy0 .AND. v_y1==iy1)
              ! (v_x0/=ix0 .and. v_x1/=ix1)
           IF(v_z0==iz0) THEN 
              !botton-direction 
              CALL WritePolyEdgeX(mat,v_y0,v_z0)
              CALL WritePolyEdgeX(mat,v_y1,v_z0)
              !... Line-X  '---'  verlängern
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              !top-direction, inside block
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
           ELSE IF(v_z1==iz1) THEN 
              !top-direction
              CALL WritePolyEdgeX(mat,v_y0,v_z1)
              CALL WritePolyEdgeX(mat,v_y1,v_z1)
              !... Line-X  '---'  verlängern
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
              !botton-direction, inside block
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
           ELSE  ! (v_z0/=iz0 .and. v_z1/=iz1)
              !top-botton-west-east, inside block, 
              !y0-y1-flanking
              !... Line-Z/-X  '---'  verlängern
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
              CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
           END IF
        END IF   
     !......................................................
     ELSE IF(v_z0==iz0.AND.v_z1==iz1) THEN 
        IF(v_y0==iy0) THEN
            IF(v_x0==ix0) THEN  !'s,w'
              CALL WritePolyEdgeZ(mat,v_x0,v_y0)
              CALL WritePolyEdgeX(mat,v_y0,v_z0)
              CALL WritePolyEdgeX(mat,v_y0,v_z1)
              CALL WritePolyEdgeY(mat,v_x0,v_z0)
              CALL WritePolyEdgeY(mat,v_x0,v_z1)
              !... Line-X/-Y  '---'  verlängern
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
            ELSE IF(v_x1==ix1) THEN  !'s,e'
              CALL WritePolyEdgeZ(mat,v_x1,v_y0)
              CALL WritePolyEdgeX(mat,v_y0,v_z0)
              CALL WritePolyEdgeX(mat,v_y0,v_z1)
              CALL WritePolyEdgeY(mat,v_x1,v_z0)
              CALL WritePolyEdgeY(mat,v_x1,v_z1)
              !... Line-X/-Y  '---'  verlängern
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
            ELSE ! (v_z0==iy0 .AND. v_z1==iy1) 
                 ! (v_y0==iy0)  .and. (v_x0/=ix0 .and. v_x1/=ix1)
              !south-direction
              CALL WritePolyEdgeX(mat,v_y0,v_z0)
              CALL WritePolyEdgeX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              !north-direction, inside block
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
            END IF
        ELSE IF(v_y1==iy1) THEN  
            IF(v_x0==ix0) THEN   !'n,w'
              CALL WritePolyEdgeZ(mat,v_x0,v_y1)
              CALL WritePolyEdgeX(mat,v_y1,v_z0)
              CALL WritePolyEdgeX(mat,v_y1,v_z1)
              CALL WritePolyEdgeY(mat,v_x0,v_z0)
              CALL WritePolyEdgeY(mat,v_x0,v_z1)
              !... Line-X/-Y  '---'  verlängern
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
            ELSE IF(v_x1==ix1) THEN !'n,e'
              CALL WritePolyEdgeZ(mat,v_x1,v_y1)
              CALL WritePolyEdgeX(mat,v_y1,v_z0)
              CALL WritePolyEdgeX(mat,v_y1,v_z1)
              CALL WritePolyEdgeY(mat,v_x1,v_z0)
              CALL WritePolyEdgeY(mat,v_x1,v_z1)
              !... Line-X/-Y  '---'  verlängern
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
            ELSE ! (v_z0==iy0 .AND. v_z1==iy1) 
                 ! (v_y1==iy1)  .and. (v_x0/=ix0 .and. v_x1/=ix1)
              !south-direction
              CALL WritePolyEdgeX(mat,v_y1,v_z0)
              CALL WritePolyEdgeX(mat,v_y1,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
              !north-direction, inside block
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
            END IF
        ELSE  ! (v_z0==iz0 .AND. v_z1==iz1)
              ! (v_y0/=iy0 .and. v_y1/=iy1)
            IF(v_x0==ix0) THEN   !'n,w'
              !west-direction
              CALL WritePolyEdgeY(mat,v_x0,v_z0)
              CALL WritePolyEdgeY(mat,v_x0,v_z1)
              !... Line-Y  '---'  verlängern
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              !east-deriction, inside block 
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
            ELSE IF(v_x1==ix1) THEN !'n,e'
              !east-direction
              CALL WritePolyEdgeY(mat,v_x1,v_z0)
              CALL WritePolyEdgeY(mat,v_x1,v_z1)
              !... Line-Y  '---'  verlängern
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
              !west-direction, inside block
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
            ELSE ! (v_z0==iy0 .AND. v_z1==iy1) 
              !east-west-north-south, inside block, 
              !z0-z1-flanking
              !... Line-X/-Y  '---'  verlängern
              CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
              CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
              CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
            END IF
        END IF
     ELSE  !block_node+part
     !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
        !......................................................
        IF(v_x0==ix0) THEN
            IF(v_y0==iy0) THEN
               IF(v_z0==iz0) THEN   !'w,s,b'
                 CALL WritePolyEdgeX(mat,v_y0,v_z0)
                 CALL WritePolyEdgeY(mat,v_x0,v_z0)
                 !... Line-X/-Y  '---'  verlängern
                 CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
                 CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
               ELSE IF(v_z1==iz1) THEN !'w,s,t'
                 CALL WritePolyEdgeX(mat,v_y0,v_z1)
                 CALL WritePolyEdgeY(mat,v_x0,v_z1)
                 !... Line-X/-Y  '---'  verlängern
                 CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
                 CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
               END IF
               CALL WritePolyEdgeZ(mat,v_x0,v_y0)
               !... Line-Z  '---'  verlängern
               CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
            ELSE IF(v_y1==iy1) THEN
               IF(v_z0==iz0) THEN    !'w,n,b'
                 CALL WritePolyEdgeX(mat,v_y1,v_z0)
                 CALL WritePolyEdgeY(mat,v_x0,v_z0)
                 !... Line-X/-Y  '---'  verlängern
                 CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
                 CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
               ELSE IF(v_z1==iz1) THEN  !'w,n,t'
                 CALL WritePolyEdgeX(mat,v_y1,v_z1)
                 CALL WritePolyEdgeY(mat,v_x0,v_z1)
                 !... Line-X/-Y  '---'  verlängern
                 CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
                 CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
               END IF
               CALL WritePolyEdgeZ(mat,v_x0,v_y1)
               !... Line-Z  '---'  verlängern
               CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
            ELSE  ! (v_x0==ix0) .and. (v_y0/=iy0) .and. (v_y1/=iy1)
               IF(v_z0==iz0) THEN
                 CALL WritePolyEdgeY(mat,v_x0,v_z0)
                 !CALL WriteThreePolyLittleY(mat,v_x0,v_z0) unten-> allg.
               ELSE IF(v_z1==iz1) THEN
                 CALL WritePolyEdgeY(mat,v_x0,v_z1)
               ! ELSE !--> überflüssig -->unten-> allg.
               END IF
            END IF
        ELSE IF(v_x1==ix1) THEN
            IF(v_y0==iy0) THEN
               IF(v_z0==iz0) THEN    !'e,s,b'
                 CALL WritePolyEdgeX(mat,v_y0,v_z0)
                 CALL WritePolyEdgeY(mat,v_x1,v_z0)
                 !... Line-X/-Y  '---'  verlängern
                 CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
                 CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
               ELSE IF(v_z1==iz1) THEN  !'e,s,t'
                 CALL WritePolyEdgeX(mat,v_y0,v_z1)
                 CALL WritePolyEdgeY(mat,v_x1,v_z1)
                 !... Line-X/-Y  '---'  verlängern
                 CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
                 CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
               END IF
               CALL WritePolyEdgeZ(mat,v_x1,v_y0)
               !... Line-Z  '---'  verlängern
               CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
            ELSE IF(v_y1==iy1) THEN
               IF(v_z0==iz0) THEN    !'e,n,b'
                 CALL WritePolyEdgeX(mat,v_y1,v_z0)
                 CALL WritePolyEdgeY(mat,v_x1,v_z0)
                 !... Line-X/-Y  '---'  verlängern
                 CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
                 CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
               ELSE IF(v_z1==iz1) THEN !'e,n,t'
                 CALL WritePolyEdgeX(mat,v_y1,v_z1)
                 CALL WritePolyEdgeY(mat,v_x1,v_z1)
                 !... Line-X/-Y  '---'  verlängern
                 CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
                 CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
               END IF
               CALL WritePolyEdgeZ(mat,v_x1,v_y1)
               !... Line-Z  '---'  verlängern
               CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
            ELSE  ! (v_x1==ix1) .and. (v_y0/=iy0) .and. (v_y1/=iy1)
               IF(v_z0==iz0) THEN
                 CALL WritePolyEdgeY(mat,v_x1,v_z0)
                 !CALL WriteThreePolyLittleY(mat,v_x1,v_z0) unten-> allg.
               ELSE IF(v_z1==iz1) THEN
                 CALL WritePolyEdgeY(mat,v_x1,v_z1)
               ! ELSE !--> überflüssig -->unten-> allg.
               END IF
            END IF
        ELSE IF(v_z0==iz0) THEN
            ! Speial Grenze-Block, Bedingung Reihenfolge x-y-z-IF's
            IF(v_y0==iy0) THEN
              CALL WritePolyEdgeX(mat,v_y0,v_z0)
            ELSE IF(v_y1==iy1) THEN
              CALL WritePolyEdgeX(mat,v_y1,v_z0)
            END IF
        ELSE IF(v_z1==iz1) THEN
            ! Speial Grenze-Block, Bedingung Reihenfolge x-y-z-IF's
            IF(v_y0==iy0) THEN
              CALL WritePolyEdgeX(mat,v_y0,v_z1)
            ELSE IF(v_y1==iy1) THEN
              CALL WritePolyEdgeX(mat,v_y1,v_z1)
            END IF
        ELSE ! "View -Part Block innerhalb Global-Block" 
        END IF
          CALL WriteThreePolyLittleZ(mat,v_x0,v_y0)
          CALL WriteThreePolyLittleZ(mat,v_x0,v_y1)
          CALL WriteThreePolyLittleZ(mat,v_x1,v_y0)
          CALL WriteThreePolyLittleZ(mat,v_x1,v_y1)
          !----------------------------------------
          CALL WriteThreePolyLittleY(mat,v_x0,v_z0)
          CALL WriteThreePolyLittleY(mat,v_x0,v_z1)
          CALL WriteThreePolyLittleY(mat,v_x1,v_z0)
          CALL WriteThreePolyLittleY(mat,v_x1,v_z1)
          !----------------------------------------
          CALL WriteThreePolyLittleX(mat,v_y0,v_z0)
          CALL WriteThreePolyLittleX(mat,v_y0,v_z1)
          CALL WriteThreePolyLittleX(mat,v_y1,v_z0)
          CALL WriteThreePolyLittleX(mat,v_y1,v_z1)
        !END IF
     END IF
  END IF   
END SUBROUTINE WriteBlockPartGMVAscii

SUBROUTINE WritePolyFaceYZ(mat,fx)
  INTEGER   :: mat,fx

  INTEGER   :: nvert=4
    Write(10,*) mat, nvert  &
      ,xParametricOut(Vertices(fx,v_y0,v_z0))   &
      ,xParametricOut(Vertices(fx,v_y0,v_z1))   &
      ,xParametricOut(Vertices(fx,v_y1,v_z1))   &
      ,xParametricOut(Vertices(fx,v_y1,v_z0))   &
      ,yParametricOut(Vertices(fx,v_y0,v_z0))   &
      ,yParametricOut(Vertices(fx,v_y0,v_z1))   &
      ,yParametricOut(Vertices(fx,v_y1,v_z1))   &
      ,yParametricOut(Vertices(fx,v_y1,v_z0))   &
      ,zParametricOut(Vertices(fx,v_y0,v_z0))   &
      ,zParametricOut(Vertices(fx,v_y0,v_z1))   &
      ,zParametricOut(Vertices(fx,v_y1,v_z1))   &
      ,zParametricOut(Vertices(fx,v_y1,v_z0))
END SUBROUTINE WritePolyFaceYZ

SUBROUTINE WritePolyFaceXZ(mat,fy)
  INTEGER   :: mat,fy

  INTEGER   :: nvert=4
    Write(10,*) mat, nvert  &
      ,xParametricOut(Vertices(v_x0,fy,v_z0))   &
      ,xParametricOut(Vertices(v_x1,fy,v_z0))   &
      ,xParametricOut(Vertices(v_x1,fy,v_z1))   &
      ,xParametricOut(Vertices(v_x0,fy,v_z1))   &
      ,yParametricOut(Vertices(v_x0,fy,v_z0))   &
      ,yParametricOut(Vertices(v_x1,fy,v_z0))   &
      ,yParametricOut(Vertices(v_x1,fy,v_z1))   &
      ,yParametricOut(Vertices(v_x0,fy,v_z1))   &
      ,zParametricOut(Vertices(v_x0,fy,v_z0))   &
      ,zParametricOut(Vertices(v_x1,fy,v_z0))   &
      ,zParametricOut(Vertices(v_x1,fy,v_z1))   &
      ,zParametricOut(Vertices(v_x0,fy,v_z1))
END SUBROUTINE WritePolyFaceXZ

SUBROUTINE WritePolyFaceXY(mat,fz)
  INTEGER   :: mat,fz

  INTEGER   :: nvert=4
    Write(10,*) mat, nvert  &
      ,xParametricOut(Vertices(v_x0,v_y0,fz))   &
      ,xParametricOut(Vertices(v_x1,v_y0,fz))   &
      ,xParametricOut(Vertices(v_x1,v_y1,fz))   &
      ,xParametricOut(Vertices(v_x0,v_y1,fz))   &
      ,yParametricOut(Vertices(v_x0,v_y0,fz))   &
      ,yParametricOut(Vertices(v_x1,v_y0,fz))   &
      ,yParametricOut(Vertices(v_x1,v_y1,fz))   &
      ,yParametricOut(Vertices(v_x0,v_y1,fz))   &
      ,zParametricOut(Vertices(v_x0,v_y0,fz))   &
      ,zParametricOut(Vertices(v_x1,v_y0,fz))   &
      ,zParametricOut(Vertices(v_x1,v_y1,fz))   &
      ,zParametricOut(Vertices(v_x0,v_y1,fz))
END SUBROUTINE WritePolyFaceXY

SUBROUTINE WritePolyEdgeX(mat,ey,ez)
  INTEGER   :: mat,ey,ez

  INTEGER   :: nvert=2
    Write(10,*) mat, nvert  &
      ,xParametricOut(Vertices(v_x0,ey,ez))   &
      ,xParametricOut(Vertices(v_x1,ey,ez))   &
      ,yParametricOut(Vertices(v_x0,ey,ez))   &
      ,yParametricOut(Vertices(v_x1,ey,ez))   &
      ,zParametricOut(Vertices(v_x0,ey,ez))   &
      ,zParametricOut(Vertices(v_x1,ey,ez))   
END SUBROUTINE WritePolyEdgeX

SUBROUTINE WritePolyEdgeY(mat,ex,ez)
  INTEGER   :: mat,ex,ez

  INTEGER   :: nvert=2
    Write(10,*) mat, nvert  &
      ,xParametricOut(Vertices(ex,v_y0,ez))   &
      ,xParametricOut(Vertices(ex,v_y1,ez))   &
      ,yParametricOut(Vertices(ex,v_y0,ez))   &
      ,yParametricOut(Vertices(ex,v_y1,ez))   &
      ,zParametricOut(Vertices(ex,v_y0,ez))   &
      ,zParametricOut(Vertices(ex,v_y1,ez))   
END SUBROUTINE WritePolyEdgeY

SUBROUTINE WritePolyEdgeZ(mat,ex,ey)
  INTEGER   :: mat,ex,ey

  INTEGER   :: nvert=2
    Write(10,*) mat, nvert  &
      ,xParametricOut(Vertices(ex,ey,v_z0))   &
      ,xParametricOut(Vertices(ex,ey,v_z1))   &
      ,yParametricOut(Vertices(ex,ey,v_z0))   &
      ,yParametricOut(Vertices(ex,ey,v_z1))   &
      ,zParametricOut(Vertices(ex,ey,v_z0))   &
      ,zParametricOut(Vertices(ex,ey,v_z1))   
END SUBROUTINE WritePolyEdgeZ

SUBROUTINE WriteThreePolyLittleX(mat,ey,ez)
  INTEGER   :: mat,ey,ez

  INTEGER   :: i, nvert=2
  REAL(8)   :: dxt
  TYPE (Vertex_T) :: VertA,VertE
  IF(v_x0==ix0) THEN   ! 'e'
    dxt=(Vertices(v_x1,ey,ez)%Point%x-Vertices(v_x1-1,ey,ez)%Point%x)/6.0d0
    VertA=Vertices(v_x1,ey,ez)
    VertE=VertA; VertE%Point%x=VertA%Point%x+dxt
    DO i=1,3
      Write(10,*) mat, nvert  &
        ,xParametricOut(VertA)   &
        ,xParametricOut(VertE)   &
        ,yParametricOut(VertA)   &
        ,yParametricOut(VertE)   &
        ,zParametricOut(VertA)   &
        ,zParametricOut(VertE)   
      VertA%Point%x=VertA%Point%x+2*dxt
      VertE%Point%x=VertE%Point%x+2*dxt 
    END DO 
  ELSE IF(v_x1==ix1) THEN  !'w'
    dxt=(Vertices(v_x0+1,ey,ez)%Point%x-Vertices(v_x0,ey,ez)%Point%x)/6.0d0
    VertA=Vertices(v_x0,ey,ez)
    VertE=VertA; VertE%Point%x=VertA%Point%x-dxt
    DO i=1,3
      Write(10,*) mat, nvert  &
        ,xParametricOut(VertA)   &
        ,xParametricOut(VertE)   &
        ,yParametricOut(VertA)   &
        ,yParametricOut(VertE)   &
        ,zParametricOut(VertA)   &
        ,zParametricOut(VertE)   
      VertA%Point%x=VertA%Point%x-2*dxt
      VertE%Point%x=VertE%Point%x-2*dxt 
    END DO 
  ELSE  ! 'e' and 'w'
    ! 'e'
    dxt=(Vertices(v_x1,ey,ez)%Point%x-Vertices(v_x1-1,ey,ez)%Point%x)/6.0d0
    VertA=Vertices(v_x1,ey,ez)
    VertE=VertA; VertE%Point%x=VertA%Point%x+dxt
    DO i=1,3
      Write(10,*) mat, nvert  &
        ,xParametricOut(VertA)   &
        ,xParametricOut(VertE)   &
        ,yParametricOut(VertA)   &
        ,yParametricOut(VertE)   &
        ,zParametricOut(VertA)   &
        ,zParametricOut(VertE)   
      VertA%Point%x=VertA%Point%x+2*dxt
      VertE%Point%x=VertE%Point%x+2*dxt 
    END DO 
    !'w'
    dxt=(Vertices(v_x0+1,ey,ez)%Point%x-Vertices(v_x0,ey,ez)%Point%x)/6.0d0
    VertA=Vertices(v_x0,ey,ez)
    VertE=VertA; VertE%Point%x=VertA%Point%x-dxt
    DO i=1,3
      Write(10,*) mat, nvert  &
        ,xParametricOut(VertA)   &
        ,xParametricOut(VertE)   &
        ,yParametricOut(VertA)   &
        ,yParametricOut(VertE)   &
        ,zParametricOut(VertA)   &
        ,zParametricOut(VertE)   
      VertA%Point%x=VertA%Point%x-2*dxt
      VertE%Point%x=VertE%Point%x-2*dxt 
    END DO 
  END IF
END SUBROUTINE WriteThreePolyLittleX

SUBROUTINE WriteThreePolyLittleY(mat,ex,ez)
  INTEGER   :: mat,ex,ez

  INTEGER   :: i, nvert=2
  REAL(8)   :: dyt
  TYPE (Vertex_T) :: VertA,VertE
  IF(v_y0==iy0) THEN   !'s'
     dyt=(Vertices(ex,v_y1,ez)%Point%y-Vertices(ex,v_y1-1,ez)%Point%y)/6.0d0
     VertA=Vertices(ex,v_y1,ez)
     VertE=VertA; VertE%Point%y=VertA%Point%y+dyt
     DO i=1,3
       Write(10,*) mat, nvert  &
         ,xParametricOut(VertA)   &
         ,xParametricOut(VertE)   &
         ,yParametricOut(VertA)   &
         ,yParametricOut(VertE)   &
         ,zParametricOut(VertA)   &
         ,zParametricOut(VertE)   
       VertA%Point%y=VertA%Point%y+2*dyt
       VertE%Point%y=VertE%Point%y+2*dyt 
     END DO
  ELSE IF(v_y1==iy1) THEN   ! 'n'
     dyt=(Vertices(ex,v_y0+1,ez)%Point%y-Vertices(ex,v_y0,ez)%Point%y)/6.0d0
     VertA=Vertices(ex,v_y0,ez)
     VertE=VertA; VertE%Point%y=VertA%Point%y-dyt
     DO i=1,3
       Write(10,*) mat, nvert  &
         ,xParametricOut(VertA)   &
         ,xParametricOut(VertE)   &
         ,yParametricOut(VertA)   &
         ,yParametricOut(VertE)   &
         ,zParametricOut(VertA)   &
         ,zParametricOut(VertE)   
       VertA%Point%y=VertA%Point%y-2*dyt
       VertE%Point%y=VertE%Point%y-2*dyt 
     END DO
  ELSE ! 's' and 'n'
     !'s'
     dyt=(Vertices(ex,v_y1,ez)%Point%y-Vertices(ex,v_y1-1,ez)%Point%y)/6.0d0
     VertA=Vertices(ex,v_y1,ez)
     VertE=VertA; VertE%Point%y=VertA%Point%y+dyt
     DO i=1,3
       Write(10,*) mat, nvert  &
         ,xParametricOut(VertA)   &
         ,xParametricOut(VertE)   &
         ,yParametricOut(VertA)   &
         ,yParametricOut(VertE)   &
         ,zParametricOut(VertA)   &
         ,zParametricOut(VertE)   
       VertA%Point%y=VertA%Point%y+2*dyt
       VertE%Point%y=VertE%Point%y+2*dyt 
     END DO
     ! 'n'
     dyt=(Vertices(ex,v_y0+1,ez)%Point%y-Vertices(ex,v_y0,ez)%Point%y)/6.0d0
     VertA=Vertices(ex,v_y0,ez)
     VertE=VertA; VertE%Point%y=VertA%Point%y-dyt
     DO i=1,3
       Write(10,*) mat, nvert  &
         ,xParametricOut(VertA)   &
         ,xParametricOut(VertE)   &
         ,yParametricOut(VertA)   &
         ,yParametricOut(VertE)   &
         ,zParametricOut(VertA)   &
         ,zParametricOut(VertE)   
       VertA%Point%y=VertA%Point%y-2*dyt
       VertE%Point%y=VertE%Point%y-2*dyt 
     END DO
  END IF
END SUBROUTINE WriteThreePolyLittleY


SUBROUTINE WriteThreePolyLittleZ(mat,ex,ey)
  INTEGER   :: mat,ex,ey

  INTEGER   :: i,nvert=2
  REAL(8)   :: dzt
  TYPE (Vertex_T) :: VertA,VertE
  IF(v_z0==iz0) THEN  ! 'b'
     dzt=(Vertices(ex,ey,v_z1)%Point%z-Vertices(ex,ey,v_z1-1)%Point%z)/6.0d0
     VertA=Vertices(ex,ey,v_z1)
     VertE=VertA; VertE%Point%z=VertA%Point%z+dzt
     DO i=1,3
       Write(10,*) mat, nvert  &
        ,xParametricOut(VertA)   &
        ,xParametricOut(VertE)   &
        ,yParametricOut(VertA)   &
        ,yParametricOut(VertE)   &
        ,zParametricOut(VertA)   &
        ,zParametricOut(VertE)
       VertA%Point%z=VertA%Point%z+2*dzt
       VertE%Point%z=VertE%Point%z+2*dzt 
     END DO
  ELSE IF(v_z1==iz1) THEN  !'t'
     dzt=(Vertices(ex,ey,v_z0+1)%Point%z-Vertices(ex,ey,v_z0)%Point%z)/6.0d0
     VertA=Vertices(ex,ey,v_z0)
     VertE=VertA; VertE%Point%z=VertA%Point%z-dzt
     DO i=1,3
       Write(10,*) mat, nvert  &
        ,xParametricOut(VertA)   &
        ,xParametricOut(VertE)   &
        ,yParametricOut(VertA)   &
        ,yParametricOut(VertE)   &
        ,zParametricOut(VertA)   &
        ,zParametricOut(VertE)
       VertA%Point%z=VertA%Point%z-2*dzt
       VertE%Point%z=VertE%Point%z-2*dzt 
     END DO
  ELSE  ! 'b' and 't'
     !'b'
     dzt=(Vertices(ex,ey,v_z1)%Point%z-Vertices(ex,ey,v_z1-1)%Point%z)/6.0d0
     VertA=Vertices(ex,ey,v_z1)
     VertE=VertA; VertE%Point%z=VertA%Point%z+dzt
     DO i=1,3
       Write(10,*) mat, nvert  &
        ,xParametricOut(VertA)   &
        ,xParametricOut(VertE)   &
        ,yParametricOut(VertA)   &
        ,yParametricOut(VertE)   &
        ,zParametricOut(VertA)   &
        ,zParametricOut(VertE)
       VertA%Point%z=VertA%Point%z+2*dzt
       VertE%Point%z=VertE%Point%z+2*dzt 
     END DO
     !'t'
     dzt=(Vertices(ex,ey,v_z0+1)%Point%z-Vertices(ex,ey,v_z0)%Point%z)/6.0d0
     VertA=Vertices(ex,ey,v_z0)
     VertE=VertA; VertE%Point%z=VertA%Point%z-dzt
     DO i=1,3
       Write(10,*) mat, nvert  &
        ,xParametricOut(VertA)   &
        ,xParametricOut(VertE)   &
        ,yParametricOut(VertA)   &
        ,yParametricOut(VertE)   &
        ,zParametricOut(VertA)   &
        ,zParametricOut(VertE)
       VertA%Point%z=VertA%Point%z-2*dzt
       VertE%Point%z=VertE%Point%z-2*dzt 
     END DO
  END IF
END SUBROUTINE WriteThreePolyLittleZ

SUBROUTINE WriteAllHausMultiColorAsciiGMV
  INTEGER   :: h,f,i,l,k,nr_p_faces
  INTEGER   :: nr_color_w,nr_farbe
  !WRITE(10,'(a8)') polygons
  !jedes Haus Waende verschiedene Farben
  !nr_farbe zum setzen von mat->gmv; 1 muss entfallen entspricht Farbe Gitterlinien!

    nr_color_w=1    !! nicht 0 initialisieren
    DO h=1,NumberHaus
      DO f=1,Haus(h)%NumberOfFaces
        nr_p_faces=Haus(h)%Faces(f)%NumberOfPoints
        !Farbe setzen, Standard 1-5 routieren, 6 für Blöcke
        IF (Haus(h)%Faces(f)%type=='r') THEN
           nr_farbe=5
        ELSE !Haus(h)%Faces(f)%type=='w'
          IF (nr_color_w<4) THEN
            nr_color_w=nr_color_w+1
            nr_farbe=nr_color_w
          ELSE
            nr_color_w=2
            nr_farbe=nr_color_w
          END IF
        END IF ! r,w
        !Output
        Write(10,*)  nr_farbe, nr_p_faces &
             ,(xPointParametricOut(Haus(h)%Faces(f)%Points(i)),i=1,nr_p_faces) &
             ,(yPointParametricOut(Haus(h)%Faces(f)%Points(l)),l=1,nr_p_faces) &
             ,(zPointParametricOut(Haus(h)%Faces(f)%Points(k)),k=1,nr_p_faces)
      END DO
    END DO  !h
  !WRITE(10,'(a8)') endpoly

END SUBROUTINE WriteAllHausMultiColorAsciiGMV


SUBROUTINE WriteAllHausOneColorAsciiGMV
  INTEGER   :: h,f,i,l,k,nr_p_faces
  INTEGER   :: nr_farbe

  !nr_farbe zum setzen von mat->gmv; 1 muss entfallen entspricht Farbe Gitterlinien!
    DO h=1,NumberHaus
      DO f=1,Haus(h)%NumberOfFaces
        nr_p_faces=Haus(h)%Faces(f)%NumberOfPoints

        !Farbe setzen, Dach->5 routieren, Wand->2, 6 fuer Bloecke vorbehalten
        IF (Haus(h)%Faces(f)%type=='r') THEN
           nr_farbe=5
        ELSE !Haus(h)%Faces(f)%type=='w'
           nr_farbe=2
        END IF ! r,w

        !Output
        Write(10,*)  nr_farbe, nr_p_faces &
             ,(xPointParametricOut(Haus(h)%Faces(f)%Points(i)),i=1,nr_p_faces) &
             ,(yPointParametricOut(Haus(h)%Faces(f)%Points(l)),l=1,nr_p_faces) &
             ,(zPointParametricOut(Haus(h)%Faces(f)%Points(k)),k=1,nr_p_faces)
      END DO
    END DO  !h
  !WRITE(10,'(a8)') endpoly

END SUBROUTINE WriteAllHausOneColorAsciiGMV

SUBROUTINE WriteAllHausColorAsciiGMV
  INTEGER   :: h,f,i,l,k,nr_p_faces
  INTEGER   :: nr_color_w,nr_farbe
  CHARACTER :: vor_type
  !WRITE(10,'(a8)') polygons
  !jedes Haus eine Farbe


    nr_color_w=1    !! nicht 0 initialisieren
    DO h=1,NumberHaus
      !Farbe setzen, Standard 2-5 routieren, 6 für Blöcke
      IF (nr_color_w<4) THEN
         nr_color_w=nr_color_w+1
         nr_farbe=nr_color_w
      ELSE
         nr_color_w=2
         nr_farbe=nr_color_w
      END IF

      DO f=1,Haus(h)%NumberOfFaces
        nr_p_faces=Haus(h)%Faces(f)%NumberOfPoints
        !Farbe setzen, Standard 2-5 routieren, 6 für Blöcke
        IF (Haus(h)%Faces(f)%type=='r') THEN
           nr_farbe=5
        ELSE
           nr_farbe=nr_color_w
        END IF
        !Output
        Write(10,*)  nr_farbe, nr_p_faces &
             ,(xPointParametricOut(Haus(h)%Faces(f)%Points(i)),i=1,nr_p_faces) &
             ,(yPointParametricOut(Haus(h)%Faces(f)%Points(l)),l=1,nr_p_faces) &
             ,(zPointParametricOut(Haus(h)%Faces(f)%Points(k)),k=1,nr_p_faces)
      END DO
    END DO  !h
  !WRITE(10,'(a8)') endpoly
END SUBROUTINE WriteAllHausColorAsciiGMV


SUBROUTINE WriteAllHausMultiColorBinGMV
  INTEGER   :: h,f,i,nr_p_faces
  INTEGER   :: nr_color_w,nr_farbe,nr_fvor
  CHARACTER :: vor_type
  !WRITE(10,'(a8)') polygons

    nr_color_w=1    !! nicht 0 initialisieren
    DO h=1,NumberHaus
      DO f=1,Haus(h)%NumberOfFaces
        nr_p_faces=Haus(h)%Faces(f)%NumberOfPoints
        !Farbe setzen, Standard 1-5 routieren, 6 für Blöcke
        IF (Haus(h)%Faces(f)%type=='r') THEN
           nr_farbe=5
        ELSE !Haus(h)%Faces(f)%type=='w'
          IF (nr_color_w<4) THEN
            nr_color_w=nr_color_w+1
            nr_farbe=nr_color_w
          ELSE
            nr_color_w=2
            nr_farbe=nr_color_w
          END IF
        END IF ! r,w

        !Output
        !Write(10,*)  nr_farbe, nr_p_faces &
        !     ,(Haus(h)%Faces(f)%Points(i)%x,i=1,nr_p_faces) &
        !     ,(Haus(h)%Faces(f)%Points(l)%y,l=1,nr_p_faces) &
        !     ,(Haus(h)%Faces(f)%Points(k)%z,k=1,nr_p_faces)
        nRec=nRec+1
        WRITE(10,REC=nRec) nr_farbe
        nRec=nRec+1
        WRITE(10,REC=nRec) nr_p_faces
        DO i=1,nr_p_faces
          nRec=nRec+1
          tempw=xPointParametricOut(Haus(h)%Faces(f)%Points(i))
          WRITE(10,REC=nRec) tempw
        END DO
        DO i=1,nr_p_faces
          nRec=nRec+1
          tempw=yPointParametricOut(Haus(h)%Faces(f)%Points(i))
          WRITE(10,REC=nRec) tempw
        END DO
        DO i=1,nr_p_faces
          nRec=nRec+1
          tempw=zPointParametricOut(Haus(h)%Faces(f)%Points(i))
          WRITE(10,REC=nRec) tempw
        END DO
        vor_type=Haus(h)%Faces(f)%type
      END DO
    END DO  !h
  !WRITE(10,'(a8)') endpoly

END SUBROUTINE WriteAllHausMultiColorBinGMV


SUBROUTINE WriteAllHausColorBinGMV
  INTEGER   :: h,f,i,l,k,nr_p_faces
  INTEGER   :: nr_color_w,nr_farbe
  CHARACTER :: vor_type
  !WRITE(10,'(a8)') polygon

    nr_color_w=1    !! nicht 0 initialisieren
    DO h=1,NumberHaus
      !Farbe setzen, Standard 2-5 routieren, 6 fuer Bloecke
      IF (nr_color_w<4) THEN
         nr_color_w=nr_color_w+1
         nr_farbe=nr_color_w
      ELSE
         nr_color_w=2
         nr_farbe=nr_color_w
      END IF

      DO f=1,Haus(h)%NumberOfFaces
        nr_p_faces=Haus(h)%Faces(f)%NumberOfPoints
        !Farbe setzen, Standard 2-5 routieren, 6 für Blöcke
        IF (Haus(h)%Faces(f)%type=='r') THEN
           nr_farbe=5
        ELSE
           nr_farbe=nr_color_w
        END IF

        !Output
        !Write(10,*)  nr_farbe, nr_p_faces &
        !     ,(Haus(h)%Faces(f)%Points(i)%x,i=1,nr_p_faces) &
        !     ,(Haus(h)%Faces(f)%Points(l)%y,l=1,nr_p_faces) &
        !     ,(Haus(h)%Faces(f)%Points(k)%z,k=1,nr_p_faces)
        nRec=nRec+1
        WRITE(10,REC=nRec) nr_farbe
        nRec=nRec+1
        WRITE(10,REC=nRec) nr_p_faces
        DO i=1,nr_p_faces
          nRec=nRec+1
          tempw=xPointParametricOut(Haus(h)%Faces(f)%Points(i))
          WRITE(10,REC=nRec) tempw
        END DO
        DO i=1,nr_p_faces
          nRec=nRec+1
          tempw=yPointParametricOut(Haus(h)%Faces(f)%Points(i))
          WRITE(10,REC=nRec) tempw
        END DO
        DO i=1,nr_p_faces
          nRec=nRec+1
          tempw=zPointParametricOut(Haus(h)%Faces(f)%Points(i))
          WRITE(10,REC=nRec) tempw
        END DO
        vor_type=Haus(h)%Faces(f)%type
      END DO
    END DO  !h
  !WRITE(10,'(a8)') endpoly
END SUBROUTINE WriteAllHausColorBinGMV




SUBROUTINE Set_ViewOutArea()
  IF(out_area=='n') THEN
     Domain%view_ixa=Domain%ix0
     Domain%view_ixe=Domain%ix1
     Domain%view_iya=Domain%iy0
     Domain%view_iye=Domain%iy1
     Domain%view_iza=Domain%iz0
     Domain%view_ize=Domain%iz1
  END IF
END SUBROUTINE Set_ViewOutArea


SUBROUTINE Search_OutArea()

    v_x0=Domain%view_ixe
    v_x1=Domain%view_ixa
    v_y0=Domain%view_iye
    v_y1=Domain%view_iya
    v_z0=Domain%view_ize
    v_z1=Domain%view_iza
    !...........................................................
    IF (igx0>=Domain%view_ixa.AND.igx0<Domain%view_ixe) THEN
        v_x0=igx0
    ELSE IF (igx0<Domain%view_ixa.AND.igx1>Domain%view_ixa) THEN
        v_x0=Domain%view_ixa
    END IF
    IF (igx1<=Domain%view_ixe.AND.igx1>Domain%view_ixa) THEN
         v_x1=igx1
    ELSE IF (igx1>Domain%view_ixe.AND.igx0<Domain%view_ixe) THEN
         v_x1=Domain%view_ixe
    END IF
    !............................................................
    IF (igy0>=Domain%view_iya.AND.igy0<Domain%view_iye) THEN
        v_y0=igy0
    ELSE IF (igy0<Domain%view_iya.AND.igy1>Domain%view_iya) THEN
        v_y0=Domain%view_iya
    END IF
    IF (igy1<=Domain%view_iye.AND.igy1>Domain%view_iya) THEN
         v_y1=igy1
    ELSE IF (igy1>Domain%view_iye.AND.igy0<Domain%view_iye) THEN
         v_y1=Domain%view_iye
    END IF
    !.............................................................
    IF (igz0>=Domain%view_iza.AND.igz0<Domain%view_ize) THEN
        v_z0=igz0
    ELSE IF (igz0<Domain%view_iza.AND.igz1>Domain%view_iza) THEN
        v_z0=Domain%view_iza
    END IF
    IF (igz1<=Domain%view_ize.AND.igz1>Domain%view_iza) THEN
         v_z1=igz1
    ELSE IF (igz1>Domain%view_ize.AND.igz0<Domain%view_ize) THEN
         v_z1=Domain%view_ize
    END IF
    !.............................................................
    IF(v_x0.NE.Domain%view_ixe .AND. v_x1.NE.Domain%view_ixa .AND. &
       v_y0.NE.Domain%view_iye .AND. v_y1.NE.Domain%view_iya .AND. &
       v_z0.NE.Domain%view_ize .AND. v_z1.NE.Domain%view_iza ) THEN
      IF(RefineX<0) THEN
        v_zw=v_x0
        v_x0=v_zw/2**(-RefineX)
        v_zw=v_x1
        v_x1=v_zw/2**(-RefineX)
      END IF
      IF(RefineY<0) THEN
        v_zw=v_y0
        v_y0=v_zw/2**(-RefineY)
        v_zw=v_y1
        v_y1=v_zw/2**(-RefineY)
      END IF
      IF(RefineZ<0) THEN
        v_zw=v_z0
        v_z0=v_zw/2**(-RefineZ)
        v_zw=v_z1
        v_z1=v_zw/2**(-RefineZ)
      END IF
      view_cell='y'
    ELSE
      view_cell='n'
    END IF

END SUBROUTINE Search_OutArea

SUBROUTINE CheckBlockBorder
     Write(*,*) "CheckBlockBorder\n"
    IF(v_x0.EQ.ix0 .AND. v_x1.EQ.ix1 .AND. &
       v_y0.EQ.iy0 .AND. v_y1.EQ.iy1 .AND. &
       v_z0.EQ.iz0 .AND. v_z1.EQ.iz1 ) THEN
       view_block='k'
    ELSE
       view_block='p'
    END IF
END SUBROUTINE CheckBlockBorder

SUBROUTINE WriteBlockVarBorderGMVAscii
    !Check_Block_Border
    IF(v_x0.EQ.ix0 .AND. v_x1.EQ.ix1 .AND. &
       v_y0.EQ.iy0 .AND. v_y1.EQ.iy1 .AND. &
       v_z0.EQ.iz0 .AND. v_z1.EQ.iz1 ) THEN
       view_block='k'
       CALL WriteBlockGMVAscii
    ELSE
       view_block='p'
       CALL WriteBlockPartGMVAscii
    END IF
END SUBROUTINE WriteBlockVarBorderGMVAscii 



SUBROUTINE Vertex_Parametric
  INTEGER :: i
  REAL(8) :: xTemp,yTemp,zTemp

  DO i=1,nr_out
    xTemp=xParametricOut(VertOut(i))
    yTemp=yParametricOut(VertOut(i))
    zTemp=zParametricOut(VertOut(i))
    VertOut(i)%Point%x=xTemp
    VertOut(i)%Point%y=yTemp
    VertOut(i)%Point%z=zTemp
  END DO

END SUBROUTINE Vertex_Parametric


SUBROUTINE WriteAllCellsAsciiGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k
  INTEGER :: lvertout
  TYPE (Vertex_T) :: V(1:2)
  REAL(8) :: sfac
  !WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.out.gmvG',STATUS='UNKNOWN')
  !WRITE(*,*)"> ",TRIM(FileName),'.out.gmvG'
  !WRITE(*,*)"              * ",TRIM(FileName),'.out.gmvG'
  CALL Display_OutGMVBlk(FileName)

  WRITE(10,'(a8,a8)') gmvinput,ascii
  WRITE(10,'(a8,i8)') nodes,nr_out
  lvertout=SIZE(VertOut)
  If(lvertout/=nr_out) THEN
    WRITE(*,*) "-------- Warnung! ----------"
    WRITE(*,*) "Grösse VertOut = ", lvertout
    WRITE(*,*) "        nr_out = ", nr_out
    WRITE(*,*) "----------------------------"
  END IF
  !VertOut(i)%Point%x
  !WRITE(10,*) (xParametricOut(VertOut(i)),i=1,SIZE(VertOut))
  WRITE(10,*) (xParametricOut(VertOut(i)),i=1,nr_out)
  !VertOut(i)%Point%y
  !WRITE(10,*) (yParametricOut(VertOut(i)),i=1,SIZE(VertOut))
  WRITE(10,*) (yParametricOut(VertOut(i)),i=1,nr_out)
  !VertOut(i)%Point%z
  !WRITE(10,*) (zParametricOut(VertOut(i)),i=1,SIZE(VertOut))
  WRITE(10,*) (zParametricOut(VertOut(i)),i=1,nr_out)

  !WRITE(10,'(a8,i8)') cells,nr_cells
  !WRITE(10,'(a8,i8)') cells,nr_viewcells
  IF(out_area=='y') THEN
    !WRITE(*,*) 'out1',cells,nr_viewcells
    WRITE(10,'(a8,i8)') cells,nr_viewcells
    DO ib=1,nb
      Write(*,*) "                                     Block : ",ib,"\/",nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              !WRITE(10,*) "Cell(i,j,k)", i, j, k,"   WriteAllCellsAsciiGMV"
              !WRITE(*,*) "Cell(i,j,k)", i, j, k,"   WriteAllCellsAsciiGMV"
              CALL WriteCellGMVAscii(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    !WRITE(*,*) 'out2',cells,nr_viewcells
    WRITE(10,'(a8,i8)') cells,nr_cells
    DO ib=1,nb
      Write(*,*) "                                     Block : ",ib,"\/",nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCellGMVAscii(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
  ! --------------------------
  ! ----- Polygone -----------
  ! CutPlane,Bloecke,Haeuser
  !---------------------------
  WRITE(10,'(a8)') polygons
  !Write(*,*) "WriteCutPlanePolyGMVAscii"
  ! CutPlane:
  !-----------
  IF(out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              !WRITE(*,*) "Cell(i,j,k)", i, j, k, "vc=", Cell(i,j,k)%Cell%vc
              !WRITE(10,*) "Cell(i,j,k)", i, j, k, "vc=", Cell(i,j,k)%Cell%vc
              CALL WriteCutPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCutPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
   END DO  !ib
  END IF
  !LandClass-Plane
  !---------------
  !Write(*,*) "WriteLandKlPlanePolyGMVAscii"
  IF(nr_landdef>0) THEN
     IF(out_area=='y') THEN
       DO ib=1,nb
         CALL Set(Floor(ib))
         CALL Search_OutArea
         IF(view_cell=='y') THEN
           DO k=v_z0+1,v_z1
             DO j=v_y0+1,v_y1
               DO i=v_x0+1,v_x1
                 CALL WriteLandKlPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
               END DO
             END DO
           END DO
         END IF
       END DO  !ib
     ELSE
       DO ib=1,nb
         CALL Set(Floor(ib))
         DO k=iz0+1,iz1
           DO j=iy0+1,iy1
             DO i=ix0+1,ix1
               CALL WriteLandKlPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
             END DO
           END DO
         END DO
      END DO  !ib
     END IF
  END IF ! if(nr_landdef)
  ! Soil-layering:
  !-----------
  !Write(*,*) "WritePolySoilSevenLayerGMVAsciiBin"
  !dzi_soil auf Grundgitter bezogen eingestellt
  IF(nr_sb>0) THEN
     V(1:2)%Point%x=Domain%xP(0)
     V(1:2)%Point%y=Domain%yP(0)
     V(1)%Point%z=Domain%zP(0)
     V(2)%Point%z=Domain%zP(1)
     dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
     sfac=(Domain%zP(nz)-Domain%zP(0))/(dzi_soil*100)
     dzi_soil=dzi_soil*sfac !skaliere dzi_soil auf 1% von Domain%zP(nz)
     dzi_soil=dzi_soil*ScaleSoil   !extern hinzu, 1% Standard
     IF(out_wahlgrid=='C') THEN
       IF(out_area=='y') THEN
         DO ib=1,nb
           CALL Set(Floor(ib))
           CALL Search_OutArea
           IF(view_cell=='y') THEN
             DO k=v_z0+1,v_z1
               DO j=v_y0+1,v_y1
                 DO i=v_x0+1,v_x1
                     !WRITE(*,*) "Cell(",i,",",j,",",k,")", "  Cell%vc=", &
                     !    &      Cell(i,j,k)%Cell%vc,"   PolySoilSevenLayer"
                     CALL WritePolySoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
                 END DO
               END DO
             END DO
           END IF
         END DO  !ib
       ELSE
         DO ib=1,nb
           CALL Set(Floor(ib))
           DO k=iz0+1,iz1
             DO j=iy0+1,iy1
               DO i=ix0+1,ix1
                   !WRITE(10,*) "Cell(i,j,k)", i, j, k
                   !WRITE(*,*) "Cell(",i,",",j,",",k,")", "  Cell%vc=", &
                   !    &      Cell(i,j,k)%Cell%vc,"   PolySoilSevenLayer"
                 CALL WritePolySoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
               END DO
             END DO
           END DO
         END DO  !ib
       END IF
     END IF
  END IF ! if (nr_sb)
  ! Bloecke:
  !---------
  DO ib=1,nb
    CALL Set(Floor(ib))
    IF(out_area=='y') THEN
       ! Output-View
       CALL Search_OutArea
       IF(view_cell=='y') THEN
          CALL WriteBlockVarBorderGMVAscii
       END IF
    ELSE  
      ! Output-Global
      CALL WriteBlockGMVAscii
    END IF
  END DO  !ib
  ! Haeuser:
  !----------
  CALL WriteAllHausMultiColorAsciiGMV
  !CALL WriteAllHausColorAsciiGMV
  !CALL WriteAllHausOneColorAsciiGMV
  WRITE(10,'(a8)') endpoly

  !.. WRITE endgmv
  WRITE(10,'(a8)') endgmv

  CLOSE(10)

END SUBROUTINE WriteAllCellsAsciiGMV


SUBROUTINE  WriteAllCellsCutPlaneAsciiGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k

  !WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Cut.out.gmvG',STATUS='UNKNOWN')
  !WRITE(*,*)"      ",TRIM(FileName),'.Cut.out.gmvG'
  !WRITE(*,*)"              * ",TRIM(FileName),'.Cut.out.gmvG'
  CALL Display_OutCutGMVBlk(FileName)

  WRITE(10,'(a8,a8)') gmvinput,ascii
  WRITE(10,'(a8,i8)') nodes,nr_out
! WRITE(10,*) VertOut%Point%x
  WRITE(10,*) (xParametricOut(VertOut(i)),i=1,SIZE(VertOut))
! WRITE(10,*) VertOut%Point%y
  WRITE(10,*) (yParametricOut(VertOut(i)),i=1,SIZE(VertOut))
! WRITE(10,*) VertOut%Point%z
  WRITE(10,*) (zParametricOut(VertOut(i)),i=1,SIZE(VertOut))
 
  IF(out_area=='y') THEN
    WRITE(10,'(a8,i8)') cells,nr_viewcutplanecells
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              !WRITE(10,*) "Cell(i,j,k)", i, j, k
              !WRITE(*,*) "Cell(i,j,k)", i, j, k
              CALL WriteCellCutPlaneGMVAscii(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    WRITE(10,'(a8,i8)') cells,nr_cutplanecells
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCellCutPlaneGMVAscii(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF

  ! Bloecke,Haeuser
  !.. WRITE polygons,...
  WRITE(10,'(a8)') polygons
  ! Bloecke
  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL WriteBlockGMVAscii
  END DO  !ib
  ! Haeuser
  CALL WriteAllHausMultiColorAsciiGMV
  !CALL WriteAllHausColorAsciiGMV
  !CALL WriteAllHausOneColorAsciiGMV
  WRITE(10,'(a8)') endpoly

  !.. WRITE endgmv
  WRITE(10,'(a8)') endgmv

  CLOSE(10)

END SUBROUTINE  WriteAllCellsCutPlaneAsciiGMV


SUBROUTINE  WriteAllCellsCutPlaneAsciiGMVTest2(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k

  !WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Cut2.out.gmvG',STATUS='UNKNOWN')
  !WRITE(*,*)"      ",TRIM(FileName),'.Cut2.out.gmvG'
  !WRITE(*,*)"              * ",TRIM(FileName),'.Cut2.out.gmvG'
  CALL Display_OutCut2GMVBlk(FileName)

  WRITE(10,'(a8,a8)') gmvinput,ascii
  WRITE(10,'(a8,i8)') nodes,nr_cutplane
  WRITE(10,*) (xParametricOut(VertCutPOut(i)),i=1,SIZE(VertCutPOut))
  WRITE(10,*) (yParametricOut(VertCutPOut(i)),i=1,SIZE(VertCutPOut))
  WRITE(10,*) (zParametricOut(VertCutPOut(i)),i=1,SIZE(VertCutPOut))
  
  IF(out_area=='y') THEN
    WRITE(10,'(a8,i8)') cells,nr_viewcutplanecells  !1+nr_viewcutplanecells2
    DO ib=1,nb
       CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              !WRITE(10,*) "Cell(i,j,k)", i, j, k,"   ----Cut2.out"
              !WRITE(*,*) "Cell(i,j,k)", i, j, k,"   ----Cut2.out"
              CALL WriteCellCutPlaneGMVAsciiTest2(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    WRITE(10,'(a8,i8)') cells,nr_cutplanecells   !1+nr_cutplanecells2
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCellCutPlaneGMVAsciiTest2(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF

  ! Bloecke,Haeuser
  !.. WRITE polygons,...
    !Bezug Koordinaten fuer GMV-Output,  VertCutPOut()
    !VertCutPOut nur Vert-Cut Koordinaten, 
    !Darstellung aller Bloecke, je nach ViewOutput-Definition
    !eventuell Groesser
    !
  WRITE(10,'(a8)') polygons
  ! Bloecke
  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL WriteBlockGMVAscii
  END DO  !ib
!  ! Haeuser
!  CALL WriteAllHausMultiColorAsciiGMV
!  !CALL WriteAllHausColorAsciiGMV
!  !CALL WriteAllHausOneColorAsciiGMV
  WRITE(10,'(a8)') endpoly
!  
  !.. WRITE endgmv
  WRITE(10,'(a8)') endgmv

  CLOSE(10)

END SUBROUTINE  WriteAllCellsCutPlaneAsciiGMVTest2



SUBROUTINE  WriteAllCellsSoilAsciiGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k
  TYPE (Vertex_T) :: V(1:2)
  REAL(8) :: sfac
 
  IF(nr_sb>0)THEN
     !WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
     OPEN(UNIT=10,FILE=TRIM(FileName)//'.Soil.out.gmvG',STATUS='UNKNOWN')
     !WRITE(*,*)"              * ",TRIM(FileName),'.Soil.out.gmvG'
     CALL Display_OutSoilGMVBlk(FileName)
     CALL WriteAuswSoilLayerToProt
   
     WRITE(10,'(a8,a8)') gmvinput,ascii
     WRITE(10,'(a8,i8)') nodes,nr_out_soil
     WRITE(10,*) (xParametricOut(VertSoilOut(i)),i=1,SIZE(VertSoilOut))
     WRITE(10,*) (yParametricOut(VertSoilOut(i)),i=1,SIZE(VertSoilOut))
     WRITE(10,*) (zParametricOut(VertSoilOut(i)),i=1,SIZE(VertSoilOut))
     !Write(*,*) "nr_soilplanecells=",nr_soilplanecells
     !Write(*,*) "nr_viewsoilplanecells=",nr_viewsoilplanecells     
     IF(out_area=='y') THEN
       WRITE(10,'(a8,i8)') cells,nr_viewsoilplanecells*nr_lsoil
       DO ib=1,nb
         Write(*,*) leerzei18,"Block : ", ib
         CALL WriteBlkNrToProt(ib)
         CALL Set(Floor(ib))
         CALL Search_OutArea
         IF(view_cell=='y') THEN
           DO k=v_z0+1,v_z1
             DO j=v_y0+1,v_y1
               DO i=v_x0+1,v_x1
                 CALL WriteCellSoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
               END DO 
             END DO
           END DO
         END IF
       END DO  !ib
     ELSE
       WRITE(10,'(a8,i8)') cells,nr_soilplanecells*nr_lsoil
       DO ib=1,nb
         Write(*,*) "                  inside Block : ", ib
         CALL WriteBlkNrToProt(ib)
         CALL Set(Floor(ib))
         DO k=iz0+1,iz1
           DO j=iy0+1,iy1
             DO i=ix0+1,ix1
               CALL WriteCellSoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
             END DO
           END DO
         END DO
       END DO  !ib
     END IF
  
     ! Polygons
     !----------- 
     WRITE(10,'(a8)') polygons

     !LandClass-Plane
     !---------------
     !Write(*,*) "WriteLandKlPlanePolyGMVAscii"
     IF(nr_landdef>0) THEN
        IF(out_area=='y') THEN
          DO ib=1,nb
            CALL Set(Floor(ib))
            CALL Search_OutArea
            IF(view_cell=='y') THEN
              DO k=v_z0+1,v_z1
                DO j=v_y0+1,v_y1
                  DO i=v_x0+1,v_x1
                   !WRITE(*,*) "Cell(i,j,k)", i, j, k, &
                   !   &       " WriteLandKlPlanePolyGMVAscii"
                    CALL WriteLandKlPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
                  END DO
                END DO
              END DO
            END IF
          END DO  !ib
        ELSE
          DO ib=1,nb
            CALL Set(Floor(ib))
            DO k=iz0+1,iz1
              DO j=iy0+1,iy1
                DO i=ix0+1,ix1
                  CALL WriteLandKlPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
                END DO
              END DO
            END DO
         END DO  !ib
        END IF
     END IF ! if(nr_landdef)

     ! Soil-layering:
     !-----------
     !Write(*,*) "WritePolySoilSevenLayerGMVAsciiBin"
     !dzi_soil auf Grundgitter bezogen eingestellt
     IF(nr_sb>0) THEN
        V(1:2)%Point%x=Domain%xP(0)
        V(1:2)%Point%y=Domain%yP(0)
        V(1)%Point%z=Domain%zP(0)
        V(2)%Point%z=Domain%zP(1)
        dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
        IF(out_wahlgrid=='C') THEN
          IF(out_area=='y') THEN
            DO ib=1,nb
              CALL Set(Floor(ib))
              CALL Search_OutArea
              IF(view_cell=='y') THEN
                DO k=v_z0+1,v_z1
                  DO j=v_y0+1,v_y1
                    DO i=v_x0+1,v_x1
                      !WRITE(*,*) "Cell(i,j,k)", i, j, k, "  vc =",Cell(i,j,k)%Cell%vc
                      CALL WritePolySoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
                    END DO
                  END DO
                END DO
              END IF
            END DO  !ib
          ELSE
            DO ib=1,nb
              CALL Set(Floor(ib))
              DO k=iz0+1,iz1
                DO j=iy0+1,iy1
                  DO i=ix0+1,ix1
                      !WRITE(10,*) "Cell(i,j,k)", i, j, k
                    CALL WritePolySoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
                  END DO
                END DO
              END DO
            END DO  !ib
          END IF
        END IF
     END IF ! if (nr_sb)

     ! Bloecke,Haeuser
       !Bezug Koordinaten fuer GMV-Output,  VertCutPOut()
       !VertCutPOut nur Vert-Cut Koordinaten, 
       !Darstellung aller Bloecke, je nach ViewOutput-Definition
       !eventuell Groesser
       !
     ! Bloecke
     !DO ib=1,nb
     !  CALL Set(Floor(ib))
     !  CALL WriteBlockGMVAscii
     !END DO  !ib
     !! Haeuser
     !CALL WriteAllHausMultiColorAsciiGMV
     !!CALL WriteAllHausColorAsciiGMV
     !!CALL WriteAllHausOneColorAsciiGMV
     WRITE(10,'(a8)') endpoly
   !  
!     WRITE(10,'(a8,a1,a8,a80)') polygons,' ',fromfile,TRIM('"'//TRIM(FileName)//'.out'//'.gmvG'//'"')
     !.. WRITE endgmv
     WRITE(10,'(a8)') endgmv
   
     CLOSE(10)
     CALL WriteEndSoilLayerToProt
   END IF
END SUBROUTINE  WriteAllCellsSoilAsciiGMV

SUBROUTINE Prot_CellsOroGMV
   Write(OutUnitProt,*)
   Write(OutUnitProt,*) "Print-Orography :"
   Write(OutUnitProt,*) "~~~~~~~~~~~~~~~~~"
   Write(OutUnitProt,*) "nr_gen         =  ",nr_gen
   Write(OutUnitProt,*) "nr_grenzeF2    = +",nr_grenzeF2
   Write(OutUnitProt,*) "nr_grenzeF1    = +",nr_grenzeF1
   Write(OutUnitProt,*) "nr_grenzeSp3   = +",nr_grenzeSp3
   Write(OutUnitProt,*) "nr_insidehex   = +",nr_insidehex
   Write(OutUnitProt,*) "                 -----------"
   Write(OutUnitProt,*) "nr_ges_oro_out =  ",nr_ges_oro_out
   Write(OutUnitProt,*) "                 -----------"
   Write(OutUnitProt,*) "........................................"
   Write(OutUnitProt,*) "nr_cell_vc_oro   (nr_gen)                  =",nr_cell_vc_oro
   Write(OutUnitProt,*) "nr_cell_novc_oro (nr_grenzeF1+nr_grenzeF2+  "
   Write(OutUnitProt,*) "                  nr_grenzeSp3)            =",nr_cell_novc_oro
   Write(OutUnitProt,*) "........................................"
   Write(OutUnitProt,*) "aus AnalyzeAllCell"
   Write(OutUnitProt,*) "nr_viewcutcells =   ",nr_viewcutcells
   Write(OutUnitProt,*) "nr_viewcutf1    = + ",nr_viewcutf1
   Write(OutUnitProt,*) "nr_viewincells  = + ",nr_viewincells
   Write(OutUnitProt,*) "                  -----------"
   Write(OutUnitProt,*) "nr_vieworocells =   ",nr_vieworocells
   Write(OutUnitProt,*) "                  -----------"
   Write(OutUnitProt,*) "........................................"
END SUBROUTINE Prot_CellsOroGMV

SUBROUTINE Screen_Prot_CellsOroGMV
   Write(*,*) "........................................"
   Write(*,*) "Print-Orography :"
   Write(*,*) "~~~~~~~~~~~~~~~~~"
   Write(*,*) "nr_gen         =  ",nr_gen
   Write(*,*) "nr_grenzeF2    = +",nr_grenzeF2
   Write(*,*) "nr_grenzeF1    = +",nr_grenzeF1
   Write(*,*) "nr_grenzeSp3   = +",nr_grenzeSp3
   Write(*,*) "nr_insidehex   = +",nr_insidehex
   Write(*,*) "                 -----------"
   Write(*,*) "nr_ges_oro_out =  ",nr_ges_oro_out
   Write(*,*) "                 -----------"
   Write(*,*) "........................................"
   Write(*,*) "nr_cell_vc_oro   (nr_gen)                  =",nr_cell_vc_oro
   Write(*,*) "nr_cell_novc_oro (nr_grenzeF1+nr_grenzeF2+  "
   Write(*,*) "                  nr_grenzeSp3)            =",nr_cell_novc_oro
   Write(*,*) "........................................"
   Write(*,*) "aus AnalyzeAllCell"
   Write(*,*) "nr_viewcutcells =   ",nr_viewcutcells
   Write(*,*) "nr_viewcutf1    = + ",nr_viewcutf1
   Write(*,*) "nr_viewincells  = + ",nr_viewincells
   Write(*,*) "                  -----------"
   Write(*,*) "nr_vieworocells =   ",nr_vieworocells
   Write(*,*) "                  -----------"
   Write(*,*) "........................................"
END SUBROUTINE Screen_Prot_CellsOroGMV


SUBROUTINE WriteAllCellsOroAsciiGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k

  !WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Oro.out.gmvG',STATUS='UNKNOWN')
  !WRITE(*,*) "> ",TRIM(FileName),'.Oro.out.gmvG'
  !WRITE(*,*) "              * ",TRIM(FileName),'.Oro.out.gmvG'
  CALL Display_OutOroGMVBlk(FileName)

  WRITE(10,'(a8,a8)') gmvinput,ascii
  WRITE(10,'(a8,i8)') nodes,nr_inside
  WRITE(10,*) (xParametricOut(VertIn(i)),i=1,SIZE(VertIn))
  WRITE(10,*) (yParametricOut(VertIn(i)),i=1,SIZE(VertIn))
  WRITE(10,*) (zParametricOut(VertIn(i)),i=1,SIZE(VertIn))
  IF(out_area=='y') THEN
    nr_vieworocells=nr_viewincells+nr_viewcutcells+nr_viewcutf1
    WRITE(10,'(a8,i8)') cells,nr_vieworocells
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
               !WRITE(10,*) "---------------------------"
               !WRITE(10,*) "Cell(i,j,k)", i, j, k
               CALL WriteCellOroGMVAscii(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    nr_orocells=nr_incells+nr_cutcells+nr_cutf1
    WRITE(10,'(a8,i8)') cells,nr_orocells
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCellOroGMVAscii(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
  nr_ges_oro_out=nr_gen+nr_grenzeF1+nr_grenzeF2+nr_grenzeSp3+nr_insidehex

  CALL WriteAuswCellsOroGMVToProt
  CALL Prot_CellsOroGMV
  CALL WriteEndCellsOroGMVToProt
  !CALL Screen_Prot_CellsOroGMV

  !  nr_vieworocells=nr_viewincells+nr_viewcutcells

  ! --------------------------
  ! ----- Polygone -----------
  ! CutPlane,Bloecke,Haeuser
  !---------------------------
  WRITE(10,'(a8)') polygons
  ! CutPlane
  !----------
  IF(out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
             !WRITE(10,*) "Cell(i,j,k)", i, j, k
               CALL WriteOroCutPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteOroCutPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
!  ! Bloecke:
   !----------
  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL WriteBlockGMVAscii
  END DO  !ib
  ! Haeuser
  CALL WriteAllHausMultiColorAsciiGMV
  !CALL WriteAllHausColorAsciiGMV
  !CALL WriteAllHausOneColorAsciiGMV

  WRITE(10,'(a8)') endpoly

  !.. WRITE endgmv
  WRITE(10,'(a8)') endgmv

  CLOSE(10)

END SUBROUTINE WriteAllCellsOroAsciiGMV

!!SUBROUTINE WriteAllSoilAsciiGMV(FileName)
!!  CHARACTER*50 :: FileName
!!  INTEGER :: ib,i,j,k
!!  nr_cell_vc_oro=0
!!  nr_cell_vc_oro_hex=0
!!  nr_insidehexprint=0
!!
!!  WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
!!  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Soil.out.gmvG',STATUS='UNKNOWN')
!!  WRITE(*,*) "> ",TRIM(FileName),'.Soil.out.gmvG'
!!
!!  WRITE(10,'(a8,a8)') gmvinput,ascii
!!  WRITE(10,'(a8,i8)') nodes,nr_inside
!!  WRITE(10,*) (xParametricOut(VertIn(i)),i=1,SIZE(VertIn))
!!  WRITE(10,*) (yParametricOut(VertIn(i)),i=1,SIZE(VertIn))
!!  WRITE(10,*) (zParametricOut(VertIn(i)),i=1,SIZE(VertIn))
!!  IF(out_area=='y') THEN
!!    !nr_vieworocells=nr_viewincells+nr_viewcutcells
!!    WRITE(10,'(a8,i8)') cells,nr_vieworocells
!!    DO ib=1,nb
!!      CALL Set(Floor(ib))
!!      CALL Search_OutArea
!!      IF(view_cell=='y') THEN
!!        DO k=v_z0+1,v_z1
!!          DO j=v_y0+1,v_y1
!!            DO i=v_x0+1,v_x1
!!               !WRITE(10,*) "Cell(i,j,k)", i, j, k
!!               CALL WriteCellSoilGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!            END DO
!!          END DO
!!        END DO
!!      END IF
!!    END DO  !ib
!!  ELSE
!!    nr_orocells=nr_incells+nr_cutcells
!!    WRITE(10,'(a8,i8)') cells,nr_orocells
!!    DO ib=1,nb
!!      CALL Set(Floor(ib))
!!      DO k=iz0+1,iz1
!!        DO j=iy0+1,iy1
!!          DO i=ix0+1,ix1
!!            CALL WriteCellSoilGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!          END DO
!!        END DO
!!      END DO
!!    END DO  !ib
!!  END IF
!!
!!   !Write(*,*) "nr_cell_vc_oro=",nr_cell_vc_oro
!!   !Write(*,*) "nr_cell_vc_oro_hex=",nr_cell_vc_oro_hex
!!   !Write(*,*) "nr_insidehexprint=",nr_insidehexprint
!!   !Write(*,*) "Oro gesamt Cellen=",nr_insidehexprint+nr_cell_vc_oro_hex+nr_cell_vc_oro
!!   !Write(*,*) "nr_viewcutplanecells=",nr_viewcutplanecells
!!
!!  ! --------------------------
!!  ! ----- Polygone -----------
!!  ! CutPlane,Bloecke,Haeuser
!!  !---------------------------
!!  WRITE(10,'(a8)') polygons
!!  ! CutPlane
!!  !----------
!!  IF(out_area=='y') THEN
!!    DO ib=1,nb
!!      CALL Set(Floor(ib))
!!      CALL Search_OutArea
!!      IF(view_cell=='y') THEN
!!        DO k=v_z0+1,v_z1
!!          DO j=v_y0+1,v_y1
!!            DO i=v_x0+1,v_x1
!!               !WRITE(10,*) "Cell(i,j,k)", i, j, k
!!              ! CALL WriteOroCutPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!            END DO
!!          END DO
!!        END DO
!!      END IF
!!    END DO  !ib
!!  ELSE
!!    DO ib=1,nb
!!      CALL Set(Floor(ib))
!!      DO k=iz0+1,iz1
!!        DO j=iy0+1,iy1
!!          DO i=ix0+1,ix1
!!           ! CALL WriteOroCutPlanePolyGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!          END DO
!!        END DO
!!      END DO
!!    END DO  !ib
!!  END IF
!!!  ! Bloecke:
!!   !----------
!!!  DO ib=1,nb
!!!    CALL Set(Floor(ib))
!!!    CALL WriteBlockGMVAscii
!!!  END DO  !ib
!!!  ! Haeuser
!!!  CALL WriteAllHausMultiColorAsciiGMV
!!!  !CALL WriteAllHausColorAsciiGMV
!!!  !CALL WriteAllHausOneColorAsciiGMV
!!!
!!  WRITE(10,'(a8)') endpoly
!!
!!!  !.. WRITE endgmv
!!  WRITE(10,'(a8)') endgmv
!!!
!!  CLOSE(10)
!!
!!END SUBROUTINE WriteAllSoilAsciiGMV



SUBROUTINE WriteAllCellsBinGMV(FileName)
  CHARACTER*50 :: FileName

  INTEGER :: ib,i,j,k,lvertout
  TYPE (Vertex_T) :: V(1:2)
  REAL(8) :: sfac
  !WRITE(*,*) 'Output-File im binary Format wird erzeugt'
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.out.gmvG',STATUS='UNKNOWN'&
    &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)
  !WRITE(*,*)"> ",TRIM(FileName),'.out.gmvG'
  !WRITE(*,*)"              * ",TRIM(FileName),'.out.gmvG'
  CALL Display_OutGMVBlk(FileName) 

  nRec=0
  !.. WRITE(10) gmvinput,ascii/iecxi4r4
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(5:8)
  !.. WRITE(10) nodes,nr_out
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) nr_out
  lvertout=SIZE(VertOut)
  If(lvertout/=nr_out) THEN
    WRITE(*,*) "-------- Warnung! ----------"
    WRITE(*,*) "Grösse VertOut = ", lvertout
    WRITE(*,*) "        nr_out = ", nr_out
    WRITE(*,*) "----------------------------"
  END IF
  !.. WRITE(10,*) VertOut%Point%x
  !WRITE(10,*) (xParametricOut(VertOut(i)),i=1,SIZE(VertOut))
  DO i=1,nr_out
    nRec=nRec+1
    tempw=xParametricOut(VertOut(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10,*) VertOut%Point%y
  DO i=1,nr_out
    nRec=nRec+1
    tempw=yParametricOut(VertOut(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10,*) VertOut%Point%z
  DO i=1,nr_out
    nRec=nRec+1
    tempw=zParametricOut(VertOut(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10) cells,nr_cells
  nRec=nRec+1
  WRITE(10,REC=nRec) cells(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) cells(5:8)
  nRec=nRec+1
  !WRITE(10,REC=nRec) nr_cells

  IF(out_area=='y') THEN
    WRITE(10,REC=nRec) nr_viewcells
    DO ib=1,nb
      Write(*,*) "                                     Block : ",ib,"\/",nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              !Write(*,*) "Celle-binary(1,j,k)", i, ,j, k
              CALL WriteCellGMVBinary(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    WRITE(10,REC=nRec) nr_cells
    DO ib=1,nb
      Write(*,*) "                                     Block : ",ib,"\/",nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCellGMVBinary(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
  ! --------------------------
  ! ----- Polygone -----------
  ! CutPlane,Bloecke,Haeuser
  !---------------------------
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(5:8)
  !Write(*,*) "WriteCutPlanePolyGMVBinary"
  ! CutPlane:
  !----------
  IF(out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              !WRITE(*,*) "Cell(i,j,k)", i,j,k, "vc=", Cell(i,j,k)%Cell%vc
              CALL WriteCutPlanePolyGMVBinary(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCutPlanePolyGMVBinary(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
  !LandClass-Plane
  !---------------
  !Write(*,*) "WriteLandKlPlanePolyGMVBinary"
  IF(nr_landdef>0) THEN
    IF(out_area=='y') THEN
      DO ib=1,nb
        CALL Set(Floor(ib))
        CALL Search_OutArea
        IF(view_cell=='y') THEN
          DO k=v_z0+1,v_z1
            DO j=v_y0+1,v_y1
              DO i=v_x0+1,v_x1
                CALL WriteLandKlPlanePolyGMVBinary(Cell(i,j,k)%Cell,i,j,k)
              END DO
            END DO
          END DO
       END IF
      END DO  !ib
    ELSE
      DO ib=1,nb
        CALL Set(Floor(ib))
        DO k=iz0+1,iz1
          DO j=iy0+1,iy1
            DO i=ix0+1,ix1
              CALL WriteLandKlPlanePolyGMVBinary(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END DO  !ib
    END IF
  END IF  !if(nr_landdef)
  ! Soil-layering:
  !-----------
  !Write(*,*) "WritePolySoilSevenLayerGMVAsciiBin"
  !dzi_soil auf Grundgitter bezogen eingestellt
  IF(nr_sb>0) THEN
    V(1:2)%Point%x=Domain%xP(0)
    V(1:2)%Point%y=Domain%yP(0)
    V(1)%Point%z=Domain%zP(0)
    V(2)%Point%z=Domain%zP(1)
    dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
    sfac=(Domain%zP(nz)-Domain%zP(0))/(dzi_soil*100)
    dzi_soil=dzi_soil*sfac !skaliere dzi_soil auf 1% von Domain%zP(nz)
    dzi_soil=dzi_soil*ScaleSoil   !extern hinzu, 1% Standard
    IF(out_wahlgrid=='C') THEN
      IF(out_area=='y') THEN
        DO ib=1,nb
          CALL Set(Floor(ib))
          CALL Search_OutArea
          IF(view_cell=='y') THEN
            DO k=v_z0+1,v_z1
              DO j=v_y0+1,v_y1
                DO i=v_x0+1,v_x1
                    !WRITE(*,*) "Cell(i,j,k)",i,j,k, "Cell%vc=", &
                    !    &      Cell(i,j,k)%Cell%vc,"   PolySoilSevenLayer"
                   CALL WritePolySoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
                END DO
              END DO
            END DO
          END IF
        END DO  !ib
      ELSE
        DO ib=1,nb
          CALL Set(Floor(ib))
          DO k=iz0+1,iz1
            DO j=iy0+1,iy1
              DO i=ix0+1,ix1
                CALL WritePolySoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
              END DO
            END DO
          END DO
        END DO  !ib
      END IF
    END IF
  END IF
  ! Bloecke:
  !----------
  DO ib=1,nb
     CALL Set(Floor(ib))
     IF(out_area=='y') THEN
       ! Output-View
       CALL Search_OutArea
       IF(view_cell=='y') THEN
          !CALL WriteBlockVarBorderGMVBinary  !noch erstellen
          CALL WriteBlockGMVBinary
       END IF
     ELSE  
       ! Output-Global
       CALL WriteBlockGMVBinary
     END IF
   END DO  !ib

  ! Haeuser:
  !----------
  CALL WriteAllHausMultiColorBinGMV
  !CALL WriteAllHausColorBinGMV
  !CALL WriteAllHausOneColorAsciiGMV
  !.. WRITE(10) endpoly
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(5:8)
 
  !.. WRITE(10) endgmv
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(5:8)

  CLOSE(10)

END SUBROUTINE WriteAllCellsBinGMV


SUBROUTINE WriteAllCellsCutPlaneBinGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k

  !WRITE(*,*) 'Output-File im binary Format wird erzeugt'
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Cut.out.gmvG',STATUS='UNKNOWN'&
    &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)
  !WRITE(*,*) "> ",TRIM(FileName),'.Cut.out.gmvG'
  !WRITE(*,*) "              * ",TRIM(FileName),'.Cut.out.gmvG'
  CALL Display_OutCutGMVBlk(FileName)

  nRec=0
  !.. WRITE(10) gmvinput,ascii/iecxi4r4
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(5:8)
  !.. WRITE(10) nodes,nr_out
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) nr_out
  !.. WRITE(10,*) VertOut%Point%x
  !WRITE(10,*) (xParametricOut(VertOut(i)),i=1,SIZE(VertOut))
  DO i=1,nr_out
    nRec=nRec+1
    tempw=xParametricOut(VertOut(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10,*) VertOut%Point%y
  DO i=1,nr_out
    nRec=nRec+1
    tempw=yParametricOut(VertOut(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10,*) VertOut%Point%z
  DO i=1,nr_out
    nRec=nRec+1
    tempw=zParametricOut(VertOut(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10) cells,nr_viewcutplanecells
  !.. WRITE(10) cells,nr_cutplanecells
  nRec=nRec+1
  WRITE(10,REC=nRec) cells(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) cells(5:8)
  nRec=nRec+1

  IF(out_area=='y') THEN
    WRITE(10,REC=nRec) nr_viewcutplanecells
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              CALL WriteCellCutPlaneGMVBinary(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    WRITE(10,REC=nRec) nr_cutplanecells
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCellCutPlaneGMVBinary(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF

  ! Bloecke,Haeuser
  !.. WRITE polygons,...
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(5:8)
  ! Bloecke
  DO ib=1,nb
    CALL Set(Floor(ib))
    CALL WriteBlockGMVBinary
  END DO  !ib
  ! Haeuser
  CALL WriteAllHausMultiColorBinGMV
  !CALL WriteAllHausColorBinGMV
  !CALL WriteAllHausOneColorAsciiGMV
  !.. WRITE(10) endpoly 
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(5:8)

  !.. WRITE(10) endgmv
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(5:8)

  CLOSE(10)

END SUBROUTINE WriteAllCellsCutPlaneBinGMV


SUBROUTINE  WriteAllCellsSoilBinGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k

  CHARACTER*120 :: workstr
  CHARACTER(4)  :: tempfrom
  CHARACTER, POINTER ::p
  INTEGER       :: lNameOut,lerg
  TYPE (Vertex_T) :: V(1:2)

  IF (nr_sb>0) THEN
     !WRITE(*,*) 'Output-File im binary Format wird erzeugt'
     OPEN(UNIT=10,FILE=TRIM(FileName)//'.Soil.out.gmvG',STATUS='UNKNOWN'&
       &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)
     !WRITE(*,*)"              * ",TRIM(FileName),'.Soil.out.gmvG'
     CALL Display_OutSoilGMVBlk(FileName)
     CALL WriteAuswSoilLayerToProt

     !1 !DateiString für polygons fromfile "Datei"   
     !1 lNameOut=LEN(TRIM(FileName)//'.out.gmvG'//'endgmv')
     !1 !2*Anführungsstriche,FilenameLänge,plus anz.Leerzeichen für füllen def. RECL 
     !1 IF(MOD((2+lNameOut),4)==0) THEN
     !1   lerg=2+lNameOut
     !1 ELSE
     !1   lerg=2+lNameOut+4-MOD((2+lNameOut),4)
     !1 END IF
     !1 workstr(1:lerg)='"'//TRIM(FileName)//'.out.gmvG'//'"'//'endgmv'
   
     nRec=0
     !.. WRITE(10) gmvinput,ascii/iecxi4r4
     nRec=nRec+1
     WRITE(10,REC=nRec) gmvinput(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) gmvinput(5:8)
     nRec=nRec+1
     WRITE(10,REC=nRec) iecxi4r4(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) iecxi4r4(5:8)
     !.. WRITE(10) nodes,nr_out_soil
     nRec=nRec+1
     WRITE(10,REC=nRec) nodes(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) nodes(5:8)
     nRec=nRec+1
     WRITE(10,REC=nRec) nr_out_soil
     !.. WRITE(10,*) (xParametricOut(VertSoilOut(i)),i=1,SIZE(VertSoilOut))
     DO i=1,nr_out_soil
       nRec=nRec+1
       tempw=xParametricOut(VertSoilOut(i))
       WRITE(10,REC=nRec) tempw
     END DO
     !.. WRITE(10,*) VertOut%Point%y
     DO i=1,nr_out_soil
       nRec=nRec+1
       tempw=yParametricOut(VertSoilOut(i))
       WRITE(10,REC=nRec) tempw
     END DO
     !.. WRITE(10,*) VertOut%Point%z
     DO i=1,nr_out_soil
       nRec=nRec+1
       tempw=zParametricOut(VertSoilOut(i))
       WRITE(10,REC=nRec) tempw
     END DO
     !.. WRITE(10) cells
     nRec=nRec+1
     WRITE(10,REC=nRec) cells(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) cells(5:8)
     nRec=nRec+1
    
     IF(out_area=='y') THEN
       !.. WRITE(10) nr_cells for soil
       WRITE(10,REC=nRec) nr_viewsoilplanecells*nr_lsoil
       DO ib=1,nb
         Write(*,*) leerzei18,"Block : ", ib
         CALL WriteBlkNrToProt(ib)
         CALL Set(Floor(ib))
         CALL Search_OutArea
         IF(view_cell=='y') THEN
           DO k=v_z0+1,v_z1
             DO j=v_y0+1,v_y1
               DO i=v_x0+1,v_x1
                 CALL WriteCellSoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
               END DO 
             END DO
           END DO
         END IF
       END DO  !ib
     ELSE
       !.. WRITE(10) nr_cells for soil
       WRITE(10,REC=nRec) nr_soilplanecells*nr_lsoil
       DO ib=1,nb
         Write(*,*) "                  inside Block : ", ib
         CALL WriteBlkNrToProt(ib)
         CALL Set(Floor(ib))
         DO k=iz0+1,iz1
           DO j=iy0+1,iy1
             DO i=ix0+1,ix1
               CALL WriteCellSoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
             END DO
           END DO
         END DO
       END DO  !ib
     END IF

     !.. WRITE(10) polygons
     nRec=nRec+1
     WRITE(10,REC=nRec) polygons(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) polygons(5:8)
   
     !LandClass-Plane
     !---------------
     IF(nr_landdef>0) THEN
       IF(out_area=='y') THEN
         DO ib=1,nb
           CALL Set(Floor(ib))
           CALL Search_OutArea
           IF(view_cell=='y') THEN
             DO k=v_z0+1,v_z1
               DO j=v_y0+1,v_y1
                 DO i=v_x0+1,v_x1
                   !WRITE(*,*) "Cell(i,j,k)", i, j, k, &
                   !   &       " WriteLandKlPlanePolyGMVBinary"
                   CALL WriteLandKlPlanePolyGMVBinary(Cell(i,j,k)%Cell,i,j,k)
                 END DO
               END DO
             END DO
          END IF
         END DO  !ib
       ELSE
         DO ib=1,nb
           CALL Set(Floor(ib))
           DO k=iz0+1,iz1
             DO j=iy0+1,iy1
               DO i=ix0+1,ix1
                 CALL WriteLandKlPlanePolyGMVBinary(Cell(i,j,k)%Cell,i,j,k)
               END DO
             END DO
           END DO
         END DO  !ib
       END IF
     END IF  !if(nr_landdef)

     ! Soil-layering:
     !-----------
     !Write(*,*) "WritePolySoilSevenLayerGMVAsciiBin"
     !dzi_soil auf Grundgitter bezogen eingestellt
     IF(nr_sb>0) THEN
       V(1:2)%Point%x=Domain%xP(0)
       V(1:2)%Point%y=Domain%yP(0)
       V(1)%Point%z=Domain%zP(0)
       V(2)%Point%z=Domain%zP(1)
       dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
       IF(out_wahlgrid=='C') THEN
         IF(out_area=='y') THEN
           DO ib=1,nb
             CALL Set(Floor(ib))
             CALL Search_OutArea
             IF(view_cell=='y') THEN
               DO k=v_z0+1,v_z1
                 DO j=v_y0+1,v_y1
                   DO i=v_x0+1,v_x1
                     CALL WritePolySoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
                   END DO
                 END DO
               END DO
             END IF
           END DO  !ib
         ELSE
           DO ib=1,nb
             CALL Set(Floor(ib))
             DO k=iz0+1,iz1
               DO j=iy0+1,iy1
                 DO i=ix0+1,ix1
                   CALL WritePolySoilSevenLayerGMVAsciiBin(Cell(i,j,k)%Cell,i,j,k)
                 END DO
               END DO
             END DO
           END DO  !ib
         END IF
       END IF
     END IF
     ! Bloecke,Haeuser
       !Bezug Koordinaten fuer GMV-Output,  VertCutPOut()
       !VertCutPOut nur Vert-Cut Koordinaten, 
       !Darstellung aller Bloecke, je nach ViewOutput-Definition
       !eventuell Groesser
       !
     !! Bloecke
     !DO ib=1,nb
     !  CALL Set(Floor(ib))
     !  CALL WriteBlockGMVBinary
     !END DO  !ib
     !! Haeuser
     !CALL WriteAllHausMultiColorBinGMV
     !!CALL WriteAllHausColorBinGMV
     !!CALL WriteAllHausOneColorAsciiGMV
   
     !.. WRITE(10) endpoly
     nRec=nRec+1
     WRITE(10,REC=nRec) endpoly(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) endpoly(5:8)
     !.. WRITE(10) endgmv
     nRec=nRec+1
     WRITE(10,REC=nRec) endgmv(1:4)
     nRec=nRec+1
     WRITE(10,REC=nRec) endgmv(5:8)
   
   !  
     !1 !WRITE(10,'(a8,a1,a8,a80)') polygons,' ',fromfile,TRIM('"'//TRIM(FileName)//'.out.gmvG'//'"')
     !1 nRec=nRec+1
     !1 WRITE(10,REC=nRec) polygons(1:4)
     !1 nRec=nRec+1
     !1 WRITE(10,REC=nRec) polygons(5:8)
     !1 nRec=nRec+1
     !1 WRITE(10,REC=nRec) fromfile(1:4)
     !1 nRec=nRec+1
     !1 WRITE(10,REC=nRec) fromfile(5:8)
     !1 DO i=1,lerg-3,4
     !1   nRec=nRec+1
     !1  !WRITE(*,*) "workstr(i:i+3)=",workstr(i:i+3), "   i=",i,"  i+3=",i+3
     !1   tempfrom=workstr(i:i+3) 
     !1   !WRITE(*,*) "tempfrom=",tempfrom,"'"
     !1   WRITE(10,REC=nRec) tempfrom 
     !1   !WRITE(10,REC=nRec) workstr(1:lerg)
     !1 END DO
   
     !1 !endgmv-String unter binary-Output,
     !1 !  direkt ohne Leerzeichen am FileName anhängen notwendig,
     !1 !  wenn polygons fromfile....Angabe gewünscht
     !1 !!.. WRITE(10) endgmv
     !1 !nRec=nRec+1
     !1 !WRITE(10,REC=nRec) endgmv(1:4)
     !1 !nRec=nRec+1
     !1 !WRITE(10,REC=nRec) endgmv(5:8)
   
     CLOSE(10)
     CALL WriteEndSoilLayerToProt
   END IF

END SUBROUTINE  WriteAllCellsSoilBinGMV


SUBROUTINE WriteAllCellsOroBinGMV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k
  nr_cell_vc_oro=0
  nr_cell_novc_oro=0
  nr_insidehex=0

  !WRITE(*,*) 'Output-File im binary Format wird erzeugt'
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.Oro.out.gmvG',STATUS='UNKNOWN'&
    &         ,FORM='UNFORMATTED',ACCESS='DIRECT',RECL=SizeOfReal)
  !WRITE(*,*) "> ",TRIM(FileName),'.Oro.out.gmvG'
  !WRITE(*,*) "              * ",TRIM(FileName),'.Oro.out.gmvG'
  CALL Display_OutOroGMVBlk(FileName)

  nRec=0
  !.. WRITE(10) gmvinput,ascii/iecxi4r4
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) gmvinput(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) iecxi4r4(5:8)
  !.. WRITE(10) nodes,nr_inside
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) nodes(5:8)
  nRec=nRec+1
  WRITE(10,REC=nRec) nr_inside
  !WRITE(10,*) (xParametricOut(VertIn(i)),i=1,SIZE(VertIn))
  DO i=1,nr_inside
    nRec=nRec+1
    tempw=xParametricOut(VertIn(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10,*) VertIn%Point%y
  DO i=1,nr_inside
    nRec=nRec+1
    tempw=yParametricOut(VertIn(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10,*) VertIn%Point%z
  DO i=1,nr_inside
    nRec=nRec+1
    tempw=zParametricOut(VertIn(i))
    WRITE(10,REC=nRec) tempw
  END DO
  !.. WRITE(10) cells,nr_vieworocells
  nRec=nRec+1
  WRITE(10,REC=nRec) cells(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) cells(5:8)
  nRec=nRec+1

  IF(out_area=='y') THEN
    nr_vieworocells=nr_viewincells+nr_viewcutcells+nr_viewcutf1
    WRITE(10,REC=nRec) nr_vieworocells
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              CALL WriteCellOroGMVBinary(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    nr_orocells=nr_incells+nr_cutcells+nr_cutf1
    WRITE(11,REC=nRec) nr_orocells
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteCellOroGMVBinary(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF
  nr_ges_oro_out=nr_gen+nr_grenzeF1+nr_grenzeF2+nr_grenzeSp3+nr_insidehex

  CALL WriteAuswCellsOroGMVToProt
  CALL Prot_CellsOroGMV
  CALL WriteEndCellsOroGMVToProt
  !CALL Screen_Prot_CellsOroGMV

  
  ! --------------------------
  ! ----- Polygone -----------
  ! CutPlane,Bloecke,Haeuser
  !---------------------------
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) polygons(5:8)
  ! CutPlane:
  !----------
  IF(out_area=='y') THEN
    DO ib=1,nb
      CALL Set(Floor(ib))
      CALL Search_OutArea
      IF(view_cell=='y') THEN
        DO k=v_z0+1,v_z1
          DO j=v_y0+1,v_y1
            DO i=v_x0+1,v_x1
              CALL WriteOroCutPlanePolyGMVBinary(Cell(i,j,k)%Cell,i,j,k)
            END DO
          END DO
        END DO
      END IF
    END DO  !ib
  ELSE
    DO ib=1,nb
      CALL Set(Floor(ib))
      DO k=iz0+1,iz1
        DO j=iy0+1,iy1
          DO i=ix0+1,ix1
            CALL WriteOroCutPlanePolyGMVBinary(Cell(i,j,k)%Cell,i,j,k)
          END DO
        END DO
      END DO
    END DO  !ib
  END IF

  ! Bloecke:
  !----------
  DO ib=1,nb
   CALL Set(Floor(ib))
   CALL WriteBlockGMVBinary
  END DO  !ib
  ! Haeuser:
  !----------
  CALL WriteAllHausMultiColorBinGMV
  !CALL WriteAllHausColorBinGMV
  !CALL WriteAllHausOneColorAsciiGMV
  !.. WRITE(10) endpoly
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endpoly(5:8)

  !.. WRITE(10) endgmv
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(1:4)
  nRec=nRec+1
  WRITE(10,REC=nRec) endgmv(5:8)

  CLOSE(10)

END SUBROUTINE WriteAllCellsOroBinGMV

END MODULE OutputOutGmvG_Mod
