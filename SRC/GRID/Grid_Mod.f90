MODULE Grid_Mod

  USE F_Mod
  USE Floor_Mod
  USE Parametric_Mod
  USE IOControl_Mod
  USE BoundaryCond_Mod
  USE GridInput_Mod
  USE GridNeu_Mod
  IMPLICIT NONE


CONTAINS



SUBROUTINE CheckVertexInOut(Vertex)
  TYPE (Vertex_T) :: Vertex
  REAL(8)  :: dist

  dist=sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y,Vertex%Point%z)
  IF (dist>=0.0d0) THEN
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x+dxLoc,Vertex%Point%y,Vertex%Point%z))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x-dxLoc,Vertex%Point%y,Vertex%Point%z))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y+dyLoc,Vertex%Point%z))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y-dyLoc,Vertex%Point%z))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y,Vertex%Point%z+dzLoc))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y,Vertex%Point%z-dzLoc))
    IF (dist>0) THEN
      Vertex%in_out=1
    ELSE
      Vertex%in_out=0
    END IF
  ELSE
    Vertex%in_out=-1
  END IF

END SUBROUTINE CheckVertexInOut


SUBROUTINE CheckVertexRand(Vertex)
  TYPE (Vertex_T) :: Vertex
  REAL(8)  :: dist
  dist=sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y,Vertex%Point%z)
  IF (dist>0.0d0) THEN
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x+dxLoc,Vertex%Point%y,Vertex%Point%z))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x-dxLoc,Vertex%Point%y,Vertex%Point%z))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y+dyLoc,Vertex%Point%z))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y-dyLoc,Vertex%Point%z))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y,Vertex%Point%z+dzLoc))
    dist=MIN(dist,sf(nr_wahlfkt,Vertex%Point%x,Vertex%Point%y,Vertex%Point%z-dzLoc))
    IF (dist>0) THEN
      Vertex%in_out=1
    ELSE
      Vertex%in_out=0
    END IF
  ELSE IF (dist>-dist_fscv) THEN  !Standard -1.0d-12
    Vertex%in_out=0
  ELSE
    Vertex%in_out=-1
  END IF
END SUBROUTINE CheckVertexRand

SUBROUTINE SearchMaxCalcZ
  INTEGER :: ix,iy,iz
    DO ix=ix0,ix1
      DO iy=iy0,iy1
        DO iz=iz0,iz1
          IF (Vertices(ix,iy,iz)%in_out==1) THEN
            ncz=MAX(ncz,iz)
            calc_z=MAX(calc_z,Vertices(ix,iy,iz)%Point%z)
            !Write(*,*)'ix=',ix,'iy=',iy,'iz=',iz,'cz=',cz,' calc_z=',calc_z
            EXIT
          END IF
        END DO
      END DO
    END DO
   !Write(*,*) "Ergebnis SearchMaxCalcZ"
   !Write(*,*) "calc_z=", calc_z
   !Write(*,*) "ncz=    ", ncz
END SUBROUTINE SearchMaxCalcZ

SUBROUTINE InitAnalyzeArea 
  domain%x0View=domain%xP(domain%ix0)
  domain%x1View=domain%xP(domain%ix1)
  domain%y0View=domain%yP(domain%iy0)
  domain%y1View=domain%yP(domain%iy1)
  domain%z0View=domain%zP(domain%iz0)
  domain%z1View=domain%zP(domain%iz1)
END SUBROUTINE InitAnalyzeArea 

SUBROUTINE Set_Domain_Fac_LonLat

  IF (TRIM(indata_type)==geographical.OR.TRIM(indata_type)==wgs84) THEN
         ! noch Absprache ob im grid-File minus für west bzw south angegeben 
         ! akt definition ja
         Domain%x0=Domain%x0   !*f_wlon
         Domain%x1=Domain%x1   !*f_elon
         Domain%y0=Domain%y0   !*f_slat
         Domain%y1=Domain%y1   !*f_nlat
         Domain%view_xa=Domain%view_xa   !*f_wlon
         Domain%view_xe=Domain%view_xe   !*f_elon
         Domain%view_ya=Domain%view_ya   !*f_slat
         Domain%view_ye=Domain%view_ye   !*f_nlat
  END IF ! TRIM(indata_type)....

END SUBROUTINE Set_Domain_Fac_LonLat

SUBROUTINE CheckVertRandBlock

  INTEGER :: ix,iy,iz
  INTEGER :: widthx,widthy,widthz

    ! for block borders 
    widthx=((ix1+1)-(ix0-1))
    widthy=((iy1+1)-(iy0-1))
    widthz=((iz1+1)-(iz0-1))
    !east,west
    DO ix=ix0-1,ix1+1,widthx
      DO iy=iy0,iy1
        DO iz=iz0,iz1
           CALL CheckVertexRand(Vertices(ix,iy,iz))
        END DO
     END DO
    END DO
    !south,north
    DO iy=iy0-1,iy1+1,widthy
      DO ix=ix0,ix1
        DO iz=iz0,iz1
           CALL CheckVertexRand(Vertices(ix,iy,iz))
        END DO
     END DO
    END DO
    !top,bottom
    DO iz=iz0-1,iz1+1,widthz
      DO ix=ix0,ix1
        DO iy=iy0,iy1
           CALL CheckVertexRand(Vertices(ix,iy,iz))
        END DO
      END DO
    END DO

END SUBROUTINE CheckVertRandBlock

SUBROUTINE AllocateEdge(Edge,Vert1,Vert2)

  TYPE(Edge_T), POINTER :: Edge
  TYPE(Vertex_T), TARGET :: Vert1,Vert2

  IF (.NOT.ASSOCIATED(Edge)) THEN
    ALLOCATE(Edge)
    Edge%Vert1=>Vert1
    Edge%Vert2=>Vert2
    Edge%yes_sp=-1
  END IF
END SUBROUTINE AllocateEdge

SUBROUTINE AllocateFace(Face,Edge1,Edge2,Edge3,Edge4,Vol)

  TYPE(Face_T), POINTER :: Face
  TYPE(Edge_T), TARGET :: Edge1,Edge2,Edge3,Edge4
  INTERFACE
    FUNCTION Vol(FaceV)
      USE Domain_Mod
      REAL(8) :: Vol
      TYPE(Face_T) :: FaceV
    END FUNCTION Vol
  END INTERFACE

  IF (.NOT.ASSOCIATED(Face)) THEN
    ALLOCATE(Face)
    Face%Edge1=>Edge1
    Face%Edge2=>Edge2
    Face%Edge3=>Edge3
    Face%Edge4=>Edge4
  END IF
  Face%Vol=Vol(Face)
END SUBROUTINE AllocateFace

SUBROUTINE AllocateCell(Cell,Face1,Face2,Face3,Face4,Face5,Face6,Vol)

  TYPE(Cell_T), POINTER :: Cell
  TYPE(Face_T), TARGET :: Face1,Face2,Face3,Face4,Face5,Face6
  INTERFACE
    FUNCTION Vol(CellV)
      USE Domain_Mod
      REAL(8) :: Vol
      TYPE(Cell_T) :: CellV
    END FUNCTION Vol
  END INTERFACE

  IF (.NOT.ASSOCIATED(Cell)) THEN
    ALLOCATE(Cell)
    Cell%Face1=>Face1
    Cell%Face2=>Face2
    Cell%Face3=>Face3
    Cell%Face4=>Face4
    Cell%Face5=>Face5
    Cell%Face6=>Face6
  END IF
  Cell%Vol=Vol(Cell)
END SUBROUTINE AllocateCell


SUBROUTINE InitAllEdges1
  INTEGER :: i,j,k
  INTEGER :: ii,jj,kk
  INTEGER :: is,js,ks
  INTEGER :: ib,Iter

  DO ib=1,nb
    CALL Set(Floor(ib))

    DO Iter=1,3
      IF (Iter==1) THEN
        is=1
        js=0
        ks=0
      ELSE IF (Iter==2) THEN
        is=0
        js=1
        ks=0
      ELSE
        is=0
        js=0
        ks=1
      END IF
      DO i=ix0+is,ix1
        DO j=iy0+js,iy1
          DO k=iz0+ks,iz1
            IF ((Vertices(i-is,j-js,k-ks)%in_out*Vertices(i,j,k)%in_out<=0 .AND.&
!                 MIN(Vertices(i-is,j-js,k-ks)%in_out,Vertices(i,j,k)%in_out)<0 .AND. &
                 MIN(Vertices(i-is,j-js,k-ks)%nrP,Vertices(i,j,k)%nrP)>-2)) THEN
!                MIN(Vertices(i-is,j-js,k-ks)%nrP,Vertices(i,j,k)%nrP)>-2)  .OR. &
!               (Vertices(i-is,j-js,k-ks)%in_out==0 .AND. Vertices(i,j-js,k-ks)%in_out==0 .AND. &
!                Vertices(i-is,j,k-ks)%in_out==0 .AND. Vertices(i,j,k-ks)%in_out==0) ) THEN
              DO ii=MAX(i,ix0+1),MIN(i+1-is,ix1)
                DO jj=MAX(j-1,iy0),MIN(j+is+ks,iy1)
                  DO kk=MAX(k-1,iz0),MIN(k+is+js,iz1)
                    CALL AllocateEdge(Edges_X(ii,jj,kk)%Edge,Vertices(ii-1,jj,kk),Vertices(ii,jj,kk))
                  END DO
                END DO
              END DO
              DO ii=MAX(i-1,ix0),MIN(i+js+ks,ix1)
                DO jj=MAX(j,iy0+1),MIN(j+1-js,iy1)
                  DO kk=MAX(k-1,iz0),MIN(k+js+is,iz1)
                    CALL AllocateEdge(Edges_Y(ii,jj,kk)%Edge,Vertices(ii,jj-1,kk),Vertices(ii,jj,kk))
                  END DO
                END DO
              END DO
              DO ii=MAX(i-1,ix0),MIN(i+ks+js,ix1)
                DO jj=MAX(j-1,iy0),MIN(j+ks+is,iy1)
                  DO kk=MAX(k,iz0+1),MIN(k+1-ks,iz1)
                    CALL AllocateEdge(Edges_Z(ii,jj,kk)%Edge,Vertices(ii,jj,kk-1),Vertices(ii,jj,kk))
                  END DO
                END DO
              END DO
              DO ii=MAX(i-1,ix0),MIN(i+js+ks,ix1)
                DO jj=MAX(j,iy0+1),MIN(j+ks+is,iy1)
                  DO kk=MAX(k,iz0+1),MIN(k+js+is,iz1)
                    CALL AllocateFace(Faces_YZ(ii,jj,kk)%Face &
                                     ,Edges_Y(ii,jj,kk-1)%Edge   &
                                     ,Edges_Z(ii,jj,kk)%Edge &
                                     ,Edges_Y(ii,jj,kk)%Edge   &
                                     ,Edges_Z(ii,jj-1,kk)%Edge &
                                     ,VolFace_YZ)
                  END DO
                END DO
              END DO
              DO ii=MAX(i,ix0+1),MIN(i+js+ks,ix1)
                DO jj=MAX(j,iy0+1),MIN(j+is+ks,iy1)
                  DO kk=MAX(k-1,iz0),MIN(k+is+js,iz1)
                    CALL AllocateFace(Faces_XY(ii,jj,kk)%Face &
                                     ,Edges_X(ii,jj-1,kk)%Edge   &
                                     ,Edges_Y(ii,jj,kk)%Edge &
                                     ,Edges_X(ii,jj,kk)%Edge   &
                                     ,Edges_Y(ii-1,jj,kk)%Edge &
                                     ,VolFace_XY)
                  END DO
                END DO
              END DO
              DO ii=MAX(i,ix0+1),MIN(i+ks+js,ix1)
                DO jj=MAX(j-1,iy0),MIN(j+is+ks,iy1)
                  DO kk=MAX(k,iz0+1),MIN(k+is+js,iz1)
                    CALL AllocateFace(Faces_ZX(ii,jj,kk)%Face &
                                     ,Edges_Z(ii-1,jj,kk)%Edge   &
                                     ,Edges_X(ii,jj,kk)%Edge &
                                     ,Edges_Z(ii,jj,kk)%Edge   &
                                     ,Edges_X(ii,jj,kk-1)%Edge &
                                     ,VolFace_ZX)
                  END DO
                END DO
              END DO
              DO ii=MAX(i,ix0+1),MIN(i+js+ks,ix1)
                DO jj=MAX(j,iy0+1),MIN(j+is+ks,iy1)
                  DO kk=MAX(k,iz0+1),MIN(k+is+js,iz1)
                    CALL AllocateCell(Cell(ii,jj,kk)%Cell &
                                     ,Faces_XY(ii,jj,kk-1)%Face &
                                     ,Faces_XY(ii,jj,kk)%Face &
                                     ,Faces_ZX(ii,jj-1,kk)%Face &
                                     ,Faces_ZX(ii,jj,kk)%Face &
                                     ,Faces_YZ(ii-1,jj,kk)%Face &
                                     ,Faces_YZ(ii,jj,kk)%Face &
                                     ,VolCell)
                  END DO
                END DO
              END DO
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO  ! ib

END SUBROUTINE InitAllEdges1

SUBROUTINE ProtInitAllEdges
  INTEGER :: i,j,k
  INTEGER :: ii,jj,kk,iia,iie,jja,jje,kka,kke
  INTEGER :: is,js,ks
  INTEGER :: ib,Iter

  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) trenn2 
  WRITE(OutUnitProt,*) leerzei2,'           >>>>   Prot AnalyzeInitAllEdges     <<<<'
  WRITE(OutUnitProt,*) trenn2 

  DO ib=1,nb
    CALL Set(Floor(ib))

    DO Iter=1,3
      IF (Iter==1) THEN
        is=1
        js=0
        ks=0
      ELSE IF (Iter==2) THEN
        is=0
        js=1
        ks=0
      ELSE
        is=0
        js=0
        ks=1
      END IF
      Write(OutUnitProt,*) "--| -------------------------------------------|--"
      Write(OutUnitProt,*) "--|               is "," js ","  ks            |--"
      Write(OutUnitProt,*) "--|                  ",is , js ,  ks 
      Write(OutUnitProt,*) "--| -------------------------------------------|--"
      Write(OutUnitProt,*) "          i-is,j-js,k-ks  "," i+is,j+js,k+ks" 
      DO i=ix0+is,ix1
        DO j=iy0+js,iy1
          DO k=iz0+ks,iz1
            Write(OutUnitProt,*)  i-is,j-js,k-ks,"  ", i,j,k, "  IF-Abfrage" 
            IF ((Vertices(i-is,j-js,k-ks)%in_out*Vertices(i,j,k)%in_out<=0 .AND.&
!                 MIN(Vertices(i-is,j-js,k-ks)%in_out,Vertices(i,j,k)%in_out)<0 .AND. &
                 MIN(Vertices(i-is,j-js,k-ks)%nrP,Vertices(i,j,k)%nrP)>-2)) THEN
!                MIN(Vertices(i-is,j-js,k-ks)%nrP,Vertices(i,j,k)%nrP)>-2)  .OR. &
!               (Vertices(i-is,j-js,k-ks)%in_out==0 .AND. Vertices(i,j-js,k-ks)%in_out==0 .AND. &
!                Vertices(i-is,j,k-ks)%in_out==0 .AND. Vertices(i,j,k-ks)%in_out==0) ) THEN
              ! Edge_X
              Write(OutUnitProt,*) "***************  Edge_X ***********************"
              Write(OutUnitProt,*) "IF-----  (i-is,j-js,k-ks)*(  i ,j ,k) <=0 und MIN(...nrP)>-2"
              Write(OutUnitProt,*) "-------","(", i-is,j-js ,k-ks,")"," (",i,j,k,")"
              Write(OutUnitProt,*) "Edges_X(ii,jj,kk)"
              iia=MAX(i,ix0+1);iie=MIN(i+1-is,ix1)
              jja=MAX(j-1,iy0);jje=MIN(j+is+ks,iy1)
              kka=MAX(k-1,iz0);kke=MIN(k+is+js,iz1)
              Write(OutUnitProt,*) "EdgeX-Schleife:", iia,iie, " iia,iie" 
              Write(OutUnitProt,*) "               ", jja,jje, " jja,jje" 
              Write(OutUnitProt,*) "               ", kka,kke, " kka,kke" 
              DO ii=iia,iie
                DO jj=jja,jje
                  DO kk=kka,kke
                    !CALL AllocateEdge(Edges_X(ii,jj,kk)%Edge,Vertices(ii-1,jj,kk),Vertices(ii,jj,kk))
                    Write(OutUnitProt,*) "Edges_X(",ii,jj,kk,")", "  (ii,jj,kk)"
                  END DO
                END DO
              END DO
              ! Edge_Y
              Write(OutUnitProt,*) "***************  Edge_Y ***********************"
              Write(OutUnitProt,*) "IF-----  (i-is,j-js,k-ks)*(  i ,j ,k) <=0 und MIN(...nrP)>-2"
              Write(OutUnitProt,*) "-------","(", i-is,j-js ,k-ks,")"," (",i,j,k,")"
              Write(OutUnitProt,*) "Edges_Y(ii,jj,kk)"
              iia=MAX(i-1,ix0);iie=MIN(i+js+ks,ix1)
              jja=MAX(j,iy0+1);jje=MIN(j+1-js,iy1)
              kka=MAX(k-1,iz0);kke=MIN(k+js+is,iz1)
              Write(OutUnitProt,*) "EdgeY-Schleife:", iia,iie, " iia,iie" 
              Write(OutUnitProt,*) "               ", jja,jje, " jja,jje" 
              Write(OutUnitProt,*) "               ", kka,kke, " kka,kke" 
              DO ii=MAX(i-1,ix0),MIN(i+js+ks,ix1)
                DO jj=MAX(j,iy0+1),MIN(j+1-js,iy1)
                  DO kk=MAX(k-1,iz0),MIN(k+js+is,iz1)
                    !CALL AllocateEdge(Edges_Y(ii,jj,kk)%Edge,Vertices(ii,jj-1,kk),Vertices(ii,jj,kk))
                    Write(OutUnitProt,*) "Edges_Y(",ii,jj,kk,")"
                 END DO
                END DO
              END DO
              Write(OutUnitProt,*) "***************  Edge_Z ***********************"
              Write(OutUnitProt,*) "IF-----  (i-is,j-js,k-ks)*(  i ,j ,k) <=0 und MIN(...nrP)>-2"
              Write(OutUnitProt,*) "-------","(", i-is,j-js ,k-ks,")"," (",i,j,k,")"
              Write(OutUnitProt,*) "Edges_Z(ii,jj,kk)"
              iia=MAX(i-1,ix0);iie=MIN(i+ks+js,ix1)
              jja=MAX(j-1,iy0);jje=MIN(j+ks+is,iy1)
              kka=MAX(k,iz0+1);kke=MIN(k+1-ks,iz1)
              Write(OutUnitProt,*) "EdgeZ-Schleife:", iia,iie, " iia,iie" 
              Write(OutUnitProt,*) "               ", jja,jje, " jja,jje" 
              Write(OutUnitProt,*) "               ", kka,kke, " kka,kke" 
              DO ii=MAX(i-1,ix0),MIN(i+ks+js,ix1)
                DO jj=MAX(j-1,iy0),MIN(j+ks+is,iy1)
                  DO kk=MAX(k,iz0+1),MIN(k+1-ks,iz1)
                    !CALL AllocateEdge(Edges_Z(ii,jj,kk)%Edge,Vertices(ii,jj,kk-1),Vertices(ii,jj,kk))
                    Write(OutUnitProt,*) "Edges_Z(",ii,jj,kk,")"
                  END DO
                END DO
              END DO
              Write(OutUnitProt,*) "***************  Face_YZ  ***********************"
              Write(OutUnitProt,*) "IF-----  (i-is,j-js,k-ks)*(  i ,j ,k) <=0 und MIN(...nrP)>-2"
              Write(OutUnitProt,*) "-------","(", i-is,j-js ,k-ks,")"," (",i,j,k,")"
              Write(OutUnitProt,*) "Face_YZ(ii,jj,kk)"
              iia=MAX(i-1,ix0);iie=MIN(i+js+ks,ix1)
              jja=MAX(j,iy0+1);jje=MIN(j+ks+is,iy1)
              kka=MAX(k,iz0+1);kke=MIN(k+js+is,iz1)
              Write(OutUnitProt,*) "Face_YZ-Schleife:", iia,iie, " iia,iie" 
              Write(OutUnitProt,*) "                 ", jja,jje, " jja,jje" 
              Write(OutUnitProt,*) "                 ", kka,kke, " kka,kke" 
              DO ii=MAX(i-1,ix0),MIN(i+js+ks,ix1)
                DO jj=MAX(j,iy0+1),MIN(j+ks+is,iy1)
                  DO kk=MAX(k,iz0+1),MIN(k+js+is,iz1)
                    !CALL AllocateFace(Faces_YZ(ii,jj,kk)%Face &
                    !                 ,Edges_Y(ii,jj,kk-1)%Edge   &
                    !                 ,Edges_Z(ii,jj,kk)%Edge &
                    !                 ,Edges_Y(ii,jj,kk)%Edge   &
                    !                 ,Edges_Z(ii,jj-1,kk)%Edge &
                    !                 ,VolFace_YZ)
                    Write(OutUnitProt,*) "Face_YZ(",ii,jj,kk,")"
                  END DO
                END DO
              END DO
              Write(OutUnitProt,*) "***************  Face_XY  ***********************"
              Write(OutUnitProt,*) "IF-----  (i-is,j-js,k-ks)*(  i ,j ,k) <=0 und MIN(...nrP)>-2"
              Write(OutUnitProt,*) "-------","(", i-is,j-js ,k-ks,")"," (",i,j,k,")"
              Write(OutUnitProt,*) "Face_XY(ii,jj,kk)"
              iia=MAX(i,ix0+1);iie=MIN(i+js+ks,ix1)
              jja=MAX(j,iy0+1);jje=MIN(j+is+ks,iy1)
              kka=MAX(k-1,iz0);kke=MIN(k+is+js,iz1)
              Write(OutUnitProt,*) "Face_XY-Schleife:", iia,iie, " iia,iie" 
              Write(OutUnitProt,*) "                 ", jja,jje, " jja,jje" 
              Write(OutUnitProt,*) "                 ", kka,kke, " kka,kke" 
              DO ii=MAX(i,ix0+1),MIN(i+js+ks,ix1)
                DO jj=MAX(j,iy0+1),MIN(j+is+ks,iy1)
                  DO kk=MAX(k-1,iz0),MIN(k+is+js,iz1)
                   !CALL AllocateFace(Faces_XY(ii,jj,kk)%Face &
                    !                 ,Edges_X(ii,jj-1,kk)%Edge   &
                    !                 ,Edges_Y(ii,jj,kk)%Edge &
                   !                 ,Edges_X(ii,jj,kk)%Edge   &
                   !                 ,Edges_Y(ii-1,jj,kk)%Edge &
                    !                 ,VolFace_XY)
                    Write(OutUnitProt,*) "Face_XY(",ii,jj,kk,")"
                  END DO
                END DO
              END DO
              Write(OutUnitProt,*) "***************  Face_ZX  ***********************"
              Write(OutUnitProt,*) "IF-----  (i-is,j-js,k-ks)*(  i ,j ,k) <=0 und MIN(...nrP)>-2"
              Write(OutUnitProt,*) "-------","(", i-is,j-js ,k-ks,")"," (",i,j,k,")"
              Write(OutUnitProt,*) "Face_ZX(ii,jj,kk)"
              iia=MAX(i,ix0+1);iie=MIN(i+js+ks,ix1)
              jja=MAX(j-1,iy0);jje=MIN(j+is+ks,iy1)
              kka=MAX(k,iz0+1);kke=MIN(k+is+js,iz1)
              Write(OutUnitProt,*) "Face_ZX-Schleife:", iia,iie, " iia,iie" 
              Write(OutUnitProt,*) "                 ", jja,jje, " jja,jje" 
              Write(OutUnitProt,*) "                 ", kka,kke, " kka,kke" 
              DO ii=MAX(i,ix0+1),MIN(i+ks+js,ix1)
                DO jj=MAX(j-1,iy0),MIN(j+is+ks,iy1)
                  DO kk=MAX(k,iz0+1),MIN(k+is+js,iz1)
                    !CALL AllocateFace(Faces_ZX(ii,jj,kk)%Face &
                    !                 ,Edges_Z(ii-1,jj,kk)%Edge   &
                    !                 ,Edges_X(ii,jj,kk)%Edge &
                    !                 ,Edges_Z(ii,jj,kk)%Edge   &
                    !                 ,Edges_X(ii,jj,kk-1)%Edge &
                    !                 ,VolFace_ZX)
                    Write(OutUnitProt,*) "Face_ZX(",ii,jj,kk,")"
                  END DO
               END DO
             END DO
             Write(OutUnitProt,*) "***************  Celle  ***********************"
             Write(OutUnitProt,*) "IF-----  (i-is,j-js,k-ks)*(  i ,j ,k) <=0 und MIN(...nrP)>-2"
             Write(OutUnitProt,*) "-------","(", i-is,j-js ,k-ks,")"," (",i,j,k,")"
             Write(OutUnitProt,*) "Celle(ii,jj,kk)"
             iia=MAX(i,ix0+1);iie=MIN(i+js+ks,ix1)
             jja=MAX(j,iy0+1);jje=MIN(j+is+ks,iy1)
             kka=MAX(k,iz0+1);kke=MIN(k+is+js,iz1)
             Write(OutUnitProt,*) "Celle-Schleife:", iia,iie, " iia,iie" 
             Write(OutUnitProt,*) "               ", jja,jje, " jja,jje" 
             Write(OutUnitProt,*) "               ", kka,kke, " kka,kke" 
             DO ii=MAX(i,ix0+1),MIN(i+js+ks,ix1)
               DO jj=MAX(j,iy0+1),MIN(j+is+ks,iy1)
                 DO kk=MAX(k,iz0+1),MIN(k+is+js,iz1)
                    !CALL AllocateCell(Cell(ii,jj,kk)%Cell &
                    !                 ,Faces_XY(ii,jj,kk-1)%Face &
                    !                 ,Faces_XY(ii,jj,kk)%Face &
                    !                 ,Faces_ZX(ii,jj-1,kk)%Face &
                    !                 ,Faces_ZX(ii,jj,kk)%Face &
                    !                 ,Faces_YZ(ii-1,jj,kk)%Face &
                    !                 ,Faces_YZ(ii,jj,kk)%Face &
                    !                 ,VolCell)
                    Write(OutUnitProt,*) "Cell(",ii,jj,kk,")"
                  END DO
                END DO
              END DO
            END IF
          END DO
        END DO
      END DO
    END DO
  END DO  ! ib

  WRITE(OutUnitProt,*)
  WRITE(OutUnitProt,*) leerzei2,'          >>>>  ProtInitAllEdges     <<<<'
  WRITE(OutUnitProt,*) leerzei2,'          >>>>>  End  Protokoll     <<<<<'
  WRITE(OutUnitProt,*) trenn2 
  WRITE(OutUnitProt,*)
END SUBROUTINE ProtInitAllEdges

SUBROUTINE SearchPoint(InterVert,Vertex1,Vertex2)
  TYPE (Vertex_T) :: InterVert
  TYPE (Vertex_T) :: Vertex1,Vertex2

  REAL(8) :: tL,xL,yL,zL,DistL,lHang
  REAL(8) :: tR,xR,yR,zR,DistR,rHang
  REAL(8) :: t,x,y,z,Dist

  REAL(8), PARAMETER :: EpsDist=1.d-12

  xL=Vertex1%Point%x
  yL=Vertex1%Point%y
  zL=Vertex1%Point%z
  xR=Vertex2%Point%x
  yR=Vertex2%Point%y
  zR=Vertex2%Point%z
  DistL=sf(nr_wahlfkt,xL,yL,zL)
  DistR=sf(nr_wahlfkt,xR,yR,zR)
  !IF (DistL==0.0d0 .OR. DistR==0.0d0) THEN  ! noch nicht gut wegen Doppelbelegungen
  ! t=1.0d0
  !ELSE

  IF (DistL>=0.0d0) THEN
    xL=Vertex2%Point%x
    yL=Vertex2%Point%y
    zL=Vertex2%Point%z
    xR=Vertex1%Point%x
    yR=Vertex1%Point%y
    zR=Vertex1%Point%z
    Dist=DistL
    DistL=DistR
    DistR=Dist
  END IF
  tL=0.0d0
  tR=1.0d0
  IF (ABS(DistL)<ABS(DistR)) THEN
    t=tL
    Dist=DistL
  ELSE
    t=tR
    Dist=DistR
  END IF
  lHang=1.0d0
  rHang=1.0d0
  IF (DistR>=0.0d0) THEN
    DO
      t=t-(tR-tL)/(DistR-DistL)*Dist
      t=MIN(t,tR-1.0d0/3.0d0*(tR-tL))
      t=MAX(t,tL+1.0d0/3.0d0*(tR-tL))
      x=t*xR+(1.0d0-t)*xL    !Sicht Gewichtsberechnung
      y=t*yR+(1.0d0-t)*yL
      z=t*zR+(1.0d0-t)*zL
      !x=xL+(xR-xL)*t        !Sicht Pointerberechnung
      !y=yL+(yR-yL)*t 
      !z=zL+(zR-zL)*t 
      Dist=sf(nr_wahlfkt,x,y,z)
      IF (Dist>=0.0d0) THEN
        tR=t
        DistR=Dist
        DistL=DistL/lHang
        lHang=lHang+1.0d0
        rHang=1.0d0
      ELSE IF (Dist<0.0d0) THEN
        tL=t
        DistL=Dist
        DistR=DistR/rHang
        rHang=rHang+1.0d0
        lHang=1.0d0
      ELSE
        IF (ABS(t-tL)>ABS(t-tR)) THEN
          IF (tL<t) THEN
            t=t-0.5d0*EpsDist
          ELSE
            t=t+0.5d0*EpsDist
          END IF
          tL=t
          x=t*xR+(1.0d0-t)*xL    !Sicht Gewichtsberechnung
          y=t*yR+(1.0d0-t)*yL
          z=t*zR+(1.0d0-t)*zL
          !x=xL+(xR-xL)*t        !Sicht Pointerberechnung
          !y=yL+(yR-yL)*t 
          !z=zL+(zR-zL)*t 
          DistL=sf(nr_wahlfkt,x,y,z)
        ELSE
          IF (tR<t) THEN
            t=t-0.5d0*EpsDist
          ELSE
            t=t+0.5d0*EpsDist
          END IF
          tR=t
          x=t*xR+(1.0d0-t)*xL  !Sicht Gewichtsberechnung
          y=t*yR+(1.0d0-t)*yL
          z=t*zR+(1.0d0-t)*zL
          !x=xL+(xR-xL)*t       !Sicht Pointerberechnung
          !y=yL+(yR-yL)*t 
          !z=zL+(zR-zL)*t 
          DistR=sf(nr_wahlfkt,x,y,z)
        END IF
      END IF
!      IF (tR-tL<=EpsDist.OR.MAX(ABS(DistL),ABS(DistR))<=EpsDist) THEN
      IF (ABS(tR-tL)<=EpsDist) THEN
        EXIT
      END IF
    END DO
  END IF
!END IF   ! DistL==0.0.d0.OR.DistR==0.0d0
  !InterVert%Point%x=t*xR+(1.0d0-t)*xL  !-Sicht Gewichtsberechnung
  !InterVert%Point%y=t*yR+(1.0d0-t)*yL
  !InterVert%Point%z=t*zR+(1.0d0-t)*zL
  !Write(*,*) " xL=", xL, " xR=",xR, " t=",t
  !Write(*,*) " yL=", yL, " yR=",yR, " t=",t
  !Write(*,*) " zL=", zL, " zR=",zR, " t=",t
  InterVert%Point%x=xL+(xR-xL)*t       !-Sicht Pointer-Berechnung
  InterVert%Point%y=yL+(yR-yL)*t  
  InterVert%Point%z=zL+(zR-zL)*t 

END SUBROUTINE SearchPoint


SUBROUTINE AnalyzeEdge(Edge)
  TYPE(Edge_T) :: Edge

  !local variable:
  TYPE (Vertex_T) :: VertS

  Edge%in_out=Edge%Vert1%in_out+Edge%Vert2%in_out
  IF (Edge%yes_sp==-1) THEN
    IF(Edge%Vert1%in_out*Edge%Vert2%in_out==-1 .AND. &
       Edge%Vert1%Shift*Edge%Vert2%Shift==1 .AND. &
       MIN(Edge%Vert1%nrP,Edge%Vert2%nrP)>-2) THEN
       !---> bleibt reines '-1' und '1' zur Analyze

    !1 IF (Edge%Vert1%in_out*Edge%Vert2%in_out<=0 .AND. & 
    !1    MIN(Edge%Vert1%in_out,Edge%Vert2%in_out)<0 .AND. &
    !1    Edge%Vert1%Shift*Edge%Vert2%Shift==1 .AND. &
    !1   MIN(Edge%Vert1%nrP,Edge%Vert2%nrP)>-2) THEN
        !---> bleibt Vert%in_out==-1 und Vert%in_out==0 dabei
        !---> Doppelbelegung  ein Vert,mit VertS 
        !---> fliegt dann in Soil-Koordinaten-Analyze raus,(V893712,Kinzig.grid)
        !---> IF zum Ausblenden Doppelbelegung
        !V893712 ALLOCATE(Edge%VertS)
        CALL SearchPoint(VertS,Edge%Vert1,Edge%Vert2)
       
    !1a   IF(.NOT.(PointEQ(Edge%Vert1%Point,VertS%Point)).AND. &
    !1a      .NOT.(PointEQ(Edge%Vert2%Point,VertS%Point)) ) THEN
    ! 1 und 1a ergibt gleiches Ergebnis im Weight2 
            ALLOCATE(Edge%VertS)
            Edge%VertS=VertS
            Edge%yes_sp=1
            !Vector Korrektur Bsp. statt 8000.0 Rückgabe 7999.999... (x||y||z)
            IF(Edge%Vert1%Point%x==Edge%Vert2%Point%x) THEN
                  Edge%VertS%Point%x=Edge%Vert1%Point%x
            END IF
            IF (Edge%Vert1%Point%y==Edge%Vert2%Point%y) THEN
                  Edge%VertS%Point%y=Edge%Vert1%Point%y
            END IF
            IF (Edge%Vert1%Point%z==Edge%Vert2%Point%z) THEN
                  Edge%VertS%Point%z=Edge%Vert1%Point%z
            END IF
            IF (Edge%VertS%Point%x>=domain%x0View-dxViewLoc .AND. &
                Edge%VertS%Point%x<=domain%x1View+dxViewLoc .AND. &
                Edge%VertS%Point%y>=domain%y0View-dyViewLoc .AND. &
                Edge%VertS%Point%y<=domain%y1View+dyViewLoc .AND. &
                Edge%VertS%Point%z>=domain%z0View-dzViewLoc .AND. &
                Edge%VertS%Point%z<=domain%z1View+dzViewLoc) THEN
                nr_out=nr_out+1
                Edge%VertS%nrP=nr_out
                nr_cutplane=nr_cutplane+1
                Edge%VertS%nrCutP=nr_cutplane
                nr_inside=nr_inside+1
                Edge%VertS%nrInP=nr_inside
            ELSE
              Edge%VertS%nrP=-2
              Edge%VertS%nrInP=-2
              Edge%VertS%nrCutP=-2
            END IF
    !  ELSE
    !    Edge%yes_sp=0   ! ein Vert%Point gleichzeitig PointS
    !  END IF
    ELSE
      Edge%yes_sp=0
    END IF
  END IF

END SUBROUTINE AnalyzeEdge


SUBROUTINE AnalyzeAllEdges

  INTEGER :: ix,iy,iz,jx,jy,jz,i,j,k
  INTEGER :: ib,in

! Suche Schnittpunkte
  DO ib=1,nb

    CALL Set(Floor(ib))

    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
            CALL AnalyzeEdge(Edges_X(i,j,k)%Edge)
          END IF
        END DO
      END DO
    END DO

    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
            CALL AnalyzeEdge(Edges_Y(i,j,k)%Edge)
          END IF
        END DO
      END DO
    END DO

    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
            CALL AnalyzeEdge(Edges_Z(i,j,k)%Edge)
          END IF
        END DO
      END DO
    END DO

!   Kopiere Schnittpunkte zu Nachbar-Edges
    DO in=1,AnzahlNachbar

      CALL Set(Nachbars(in))

      !IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw'.OR. &
      !    Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
      !  IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
      IF ((Nachbars(in)%nType=='iw'.and.TypeW=='iw') .OR. &
          (Nachbars(in)%nType=='ie'.and.TypeE=='ie')) THEN
        IF (Nachbars(in)%nType=='iw') THEN
          ix=Floor(ibn)%ix1
          jx=jx1
        ELSE
          ix=Floor(ibn)%ix0
          jx=jx0
        END IF
!       -------------------
!                |
!          N     |    D
!                |
!       -------------------
        IF (RefineY>RefineNachbarY) THEN
!         -------------
!         |     |  |  |  'iw,ie'
!         |     |-----|
!         |     |  |  |
!         -------------
          IF (RefineZ>RefineNachbarZ) THEN
            DO jy=jy0+1,jy1,IncrY
              iy=(jy+1)/IncrY
              DO jz=jz0,jz1,IncrZ
                iz=jz/Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_Y(jx,jy+IncrY-1,jz)%Edge).AND.Edges_Y(jx,jy+IncrY-1,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy+IncrY-1,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jy=jy0+1,jy1,IncrY
              iy=(jy+1)/IncrY
              DO jz=jz0,jz1
                iz=jz*Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_Y(jx,jy+IncrY-1,jz)%Edge).AND.Edges_Y(jx,jy+IncrY-1,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy+IncrY-1,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineY==RefineNachbarY) THEN
!         -------------
!         |     |     |  'iw,ie'
!         |---------- |
!         |     |     |
!         -------------
          IF (RefineZ>RefineNachbarZ) THEN
            DO jy=jy0+1,jy1
              iy=jy
              DO jz=jz0,jz1,IncrZ
                iz=jz/Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jy=jy0+1,jy1
              iy=jy
              DO jz=jz0,jz1
                iz=jz*Incrz
                !IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                ! wieder rückgesetzt V8.9.3.3.1.5-->6Entw
                IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge)) THEN
                IF (Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineY<RefineNachbarY) THEN
!         -------------
!         |  |  |     |   'iw,ie'
!         |---- |     |
!         |  |  |     |
!         -------------
          IF (RefineZ>RefineNachbarZ) THEN
            DO jy=jy0+1,jy1
              iy=IncrY*jy
              DO jz=jz0,jz1,IncrZ
                iz=jz/IncrZ
                IF(ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Y Vert1-/Vert2%in_out==-1(1) sind, und RefineY groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Y Vert1-/Vert2%in_out==-1(1) sind, und RefineYRefineY  groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          ELSE    ! RefineZ<=RefineNachbarZ
            DO jy=jy0+1,jy1
              iy=IncrY*jy
              DO jz=jz0,jz1
                iz=jz*IncrZ
                IF(ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Y Vert1-/Vert2%in_out==-1(1) sind, RefineY groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Y Vert1-/Vert2%in_out==-1(1) sind, RefineY groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          END IF
        END IF  ! RefineY

        IF (RefineZ>RefineNachbarZ) THEN
!         -------------
!         |     |  |  |  'iw,ie'
!         |     |-----|
!         |     |  |  |
!         -------------
          IF (RefineY>RefineNachbarY) THEN
            DO jy=jy0,jy1,IncrY
              iy=jy/IncrY
              DO jz=jz0+1,jz1,IncrZ
                iz=(jz+1)/Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz+IncrZ-1)%Edge).AND.Edges_Z(jx,jy,jz+IncrZ-1)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz+IncrZ-1)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          ELSE  ! RefineY<=RefineNachbarY
            DO jy=jy0,jy1
              iy=jy*IncrY
              DO jz=jz0+1,jz1,IncrZ
                iz=(jz+1)/Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz+IncrZ-1)%Edge).AND.Edges_Z(jx,jy,jz+IncrZ-1)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz+IncrZ-1)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineZ==RefineNachbarZ) THEN
!         -------------
!         |     |     |   'iw,ie'
!         |---------- |
!         |     |     |
!         -------------
          IF (RefineY>RefineNachbarY) THEN
            DO jy=jy0,jy1,IncrY
              iy=jy/IncrY
              DO jz=jz0+1,jz1
                iz=jz
                IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jy=jy0,jy1
              iy=jy*IncrY
              DO jz=jz0+1,jz1
                iz=jz
                !IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge)) THEN
                IF (Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineZ<RefineNachbarZ) THEN
!         -------------
!         |  |  |     |  'iw,ie'
!         |---- |     |
!         |  |  |     |
!         -------------
          IF (RefineY>RefineNachbarY) THEN
            DO jy=jy0,jy1,IncrY
              iy=jy/IncrY
              DO jz=jz0+1,jz1
                iz=IncrZ*jz
                IF(ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Z Vert1-/Vert2%in_out==-1(1) sind, und RefineZ groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=0
                     ! bedeutet das Edges_Z Vert1-/Vert2%in_out==-1(1) sind, und RefineZ groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jy=jy0,jy1
              iy=jy*IncrY
              DO jz=jz0+1,jz1
                iz=IncrZ*jz
                IF(ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Z Vert1-/Vert2%in_out==-1(1) sind, und RefineZ groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=0
                     ! bedeutet das Edges_Z Vert1-/Vert2%in_out==-1(1) sind, und RefineZ groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          END IF
        END IF ! RefineZ
      END IF  ! 'iw, ie'

      !IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps'.OR. &
      !    Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
      !  IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
      IF ((Nachbars(in)%nType=='is'.and.TypeS=='is').OR. &
          (Nachbars(in)%nType=='in'.and.TypeN=='in')) THEN
        IF (Nachbars(in)%nType=='is') THEN
          iy=Floor(ibn)%iy1
          jy=jy1
        ELSE
          iy=Floor(ibn)%iy0
          jy=jy0
        END IF
!       -------------------
!                |
!          N     |    D
!                |
!       -------------------
        IF (RefineX>RefineNachbarX) THEN
!         -------------
!         |     |  |  |  'is,in'
!         |     |-----|
!         |     |  |  |
!         -------------
          IF (RefineZ>RefineNachbarZ) THEN
            DO jx=jx0+1,jx1,IncrX
              ix=(jx+1)/IncrX
              DO jz=jz0,jz1,IncrZ
                iz=jz/Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_X(jx+IncrX-1,jy,jz)%Edge).AND.Edges_X(jx+IncrX-1,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx+IncrX-1,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          ELSE  ! RefineZ<=RefineNachbarZ
            DO jx=jx0+1,jx1,IncrX
              ix=(jx+1)/IncrX
              DO jz=jz0,jz1
                iz=jz*Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_X(jx+IncrX-1,jy,jz)%Edge).AND.Edges_X(jx+IncrX-1,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx+IncrX-1,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineX==RefineNachbarX) THEN
!         -------------
!         |     |     |  'is,in'
!         |---------- |
!         |     |     |
!         -------------
          IF (RefineZ>RefineNachbarZ) THEN
            DO jx=jx0+1,jx1
              ix=jx
              DO jz=jz0,jz1,IncrZ
                iz=jz/Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jx=jx0+1,jx1
              ix=jx
              DO jz=jz0,jz1
                iz=jz*Incrz
                !IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge)) THEN
                IF (Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineX<RefineNachbarX) THEN
!         -------------
!         |  |  |     |   'is,in'
!         |---- |     |
!         |  |  |     |
!         -------------
          IF (RefineZ>RefineNachbarZ) THEN
            DO jx=jx0+1,jx1
              ix=IncrX*jx
              DO jz=jz0,jz1,IncrZ
                iz=jz/IncrZ
                IF(ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_X Vert1-/Vert2%in_out==-1(1) sind, und RefineX groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_X Vert1-/Vert2%in_out==-1(1) sind, und RefineX groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          ELSE ! RefineZ<=RefineNachbarZ
            DO jx=jx0+1,jx1
              ix=IncrX*jx
              DO jz=jz0,jz1
                iz=jz*IncrZ
                IF(ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_X Vert1-/Vert2%in_out==-1(1) sind, und RefineX groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_X Vert1-/Vert2%in_out==-1(1) sind, und RefineX groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          END IF
        END IF  ! RefineX

        IF (RefineZ>RefineNachbarZ) THEN
!         -------------
!         |     |  |  |   'is,in'
!         |     |-----|
!         |     |  |  |
!         -------------
          IF (RefineX>RefineNachbarX) THEN
            DO jx=jx0,jx1,IncrX
              ix=jx/IncrX
              DO jz=jz0+1,jz1,IncrZ
                iz=(jz+1)/Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz-1+IncrZ)%Edge).AND.Edges_Z(jx,jy,jz-1+IncrZ)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz+IncrZ-1)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          ELSE  ! RefineX<=RefineNachbarX
            DO jx=jx0,jx1
              ix=jx*IncrX
              DO jz=jz0+1,jz1,IncrZ
                iz=(jz+1)/Incrz
                IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz-1+IncrZ)%Edge).AND.Edges_Z(jx,jy,jz-1+IncrZ)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz+IncrZ-1)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineZ==RefineNachbarZ) THEN
!         -------------
!         |     |     |   'is,in'
!         |---------- |
!         |     |     |
!         -------------
          IF (RefineX>RefineNachbarX) THEN
            DO jx=jx0,jx1,IncrX
              ix=jx/IncrX
              DO jz=jz0+1,jz1
                iz=jz
                IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
              END DO
            END DO
          ELSE ! RefineX<=RefineNachbarX
            DO jx=jx0,jx1
              ix=jx*IncrX
              DO jz=jz0+1,jz1
                iz=jz
                !IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                IF (ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge)) THEN
                IF (Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineZ<RefineNachbarZ) THEN
!         -------------
!         |  |  |     |   'is,in'
!         |---- |     |
!         |  |  |     |
!         -------------
          IF (RefineX>RefineNachbarX) THEN
            DO jx=jx0,jx1,IncrX
              ix=jx/IncrX
              DO jz=jz0+1,jz1
                iz=IncrZ*jz
                IF(ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Z Vert1-/Vert2%in_out==-1(1) sind, und RefineZ groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=0
                     ! bedeutet das Edges_Z Vert1-/Vert2%in_out==-1(1) sind, und RefineZ groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          ELSE ! RefineX<=RefineNachbarX
            DO jx=jx0,jx1
              ix=jx*IncrX
              DO jz=jz0+1,jz1
                iz=IncrZ*jz
                IF(ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Z(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Z Vert1-/Vert2%in_out==-1(1) sind, und RefineZ groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge).AND.Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%VertS)
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%VertS=Edges_Z(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Z(jx,jy,jz)%Edge).AND.Edges_Z(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Z(ix,iy,iz-1)%Edge%yes_sp=0
                     ! bedeutet das Edges_Z Vert1-/Vert2%in_out==-1(1) sind, und RefineZ groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          END IF
        END IF ! RefineZ
      END IF  ! 'is, in'

      !IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb'.OR. &
      !    Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
      !  IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
      IF ((Nachbars(in)%nType=='ib'.and.TypeB=='ib') .OR. &
          (Nachbars(in)%nType=='it'.and.TypeT=='it')) THEN
        IF (Nachbars(in)%nType=='ib') THEN
          iz=Floor(ibn)%iz1
          jz=jz1
        ELSE
          iz=Floor(ibn)%iz0
          jz=jz0
        END IF
!       -------------------
!                |
!          N     |    D
!                |
!       -------------------
        IF (RefineX>RefineNachbarX) THEN
!         -------------
!         |     |  |  |   'ib,it'
!         |     |-----|
!         |     |  |  |
!         -------------
          IF (RefineY>RefineNachbarY) THEN
            DO jx=jx0+1,jx1,IncrX
              ix=(jx+1)/IncrX
              DO jy=jy0,jy1,IncrY
                iy=jy/Incry
                IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                     IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp>0) THEN
                       ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                     ELSE IF (ASSOCIATED(Edges_X(jx-1+IncrX,jy,jz)%Edge).AND.Edges_X(jx-1+IncrX,jy,jz)%Edge%yes_sp>0) THEN
                       ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx-1+IncrX,jy,jz)%Edge%VertS
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                     ELSE
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                     END IF
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jx=jx0+1,jx1,IncrX
              ix=(jx+1)/IncrX
              DO jy=jy0,jy1
                iy=jy*Incry
                IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_X(jx-1+IncrX,jy,jz)%Edge).AND.Edges_X(jx-1+IncrX,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx-1+IncrX,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineX==RefineNachbarX) THEN
!         -------------
!         |     |     |  'ib,it'
!         |---------- |
!         |     |     |
!         -------------
          IF (RefineY>RefineNachbarY) THEN
            DO jx=jx0+1,jx1
              ix=jx
              DO jy=jy0,jy1,IncrY
                iy=jy/Incry
                IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jx=jx0+1,jx1
              ix=jx
              DO jy=jy0,jy1
                iy=jy*Incry
                !IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                IF (ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge)) THEN
                IF (Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  !IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1)  THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge)) THEN
                  IF(Edges_X(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                  END IF
                END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineX<RefineNachbarX) THEN
!         -------------
!         |  |  |     |  'ib,it'
!         |---- |     |
!         |  |  |     |
!         -------------
          IF (RefineY>RefineNachbarY) THEN
            DO jx=jx0+1,jx1
              ix=IncrX*jx
              DO jy=jy0,jy1,IncrY
                iy=jy/IncrY
                IF(ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_X Vert1-/Vert2%in_out==-1(1) sind, und RefineX groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_X Vert1-/Vert2%in_out==-1(1) sind, und RefineX groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          ELSE   ! RefineY<=RefineNachbarY
            DO jx=jx0+1,jx1
              ix=IncrX*jx
              DO jy=jy0,jy1
                iy=jy*IncrY
                IF(ASSOCIATED(Floor(ibn)%Edges_X(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_X(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_X(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_X Vert1-/Vert2%in_out==-1(1) sind, und RefineX groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge).AND.Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%VertS=Edges_X(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_X(jx,jy,jz)%Edge).AND.Edges_X(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_X(ix-1,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_X Vert1-/Vert2%in_out==-1(1) sind, und RefineX groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          END IF
        END IF  ! RefineX

        IF (RefineY>RefineNachbarY) THEN
!         -------------
!         |     |  |  |  'ib,it'
!         |     |-----|
!         |     |  |  |
!         -------------
          IF (RefineX>RefineNachbarX) THEN
            DO jx=jx0,jx1,IncrX
              ix=jx/IncrX
              DO jy=jy0+1,jy1,IncrY
                iy=(jy+1)/Incry
                IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_Y(jx,jy-1+IncrY,jz)%Edge).AND.Edges_Y(jx,jy-1+IncrY,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy-1+IncrY,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jx=jx0,jx1
              ix=jx*IncrX
              DO jy=jy0+1,jy1,IncrY
                iy=(jy+1)/Incry
                IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge) &
                    .AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert1%in_out &
                     *Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert2%in_out == 1) THEN  !Nachbar:Vert1-,Vert2%in_out==-1(1)
                       Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                  ELSE
                    IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE IF (ASSOCIATED(Edges_Y(jx,jy-1+IncrY,jz)%Edge).AND.Edges_Y(jx,jy-1+IncrY,jz)%Edge%yes_sp>0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy-1+IncrY,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineY==RefineNachbarY) THEN
!         -------------
!         |     |     |  'ib,it'
!         |---------- |
!         |     |     |
!         -------------
          IF (RefineX>RefineNachbarX) THEN
            DO jx=jx0,jx1,IncrX
              ix=jx/IncrX
              DO jy=jy0+1,jy1
                iy=jy
                IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jx=jx0,jx1
              ix=jx*IncrX
              DO jy=jy0+1,jy1
                iy=jy
                !IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                IF (ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge)) THEN
                IF (Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  !IF((ASSOCIATED(Edges_Y(jx,jy,jz)%Edge)).AND.(Edges_Y(jx,jy,jz)%Edge%yes_sp==1))  THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge)) THEN
                  IF(Edges_Y(jx,jy,jz)%Edge%yes_sp==1)  THEN
                    ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                    Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                  END IF
                  END IF
                END IF
                END IF
              END DO
            END DO
          END IF
        ELSE IF (RefineY<RefineNachbarY) THEN
!         -------------
!         |  |  |     |  'ib,it'
!         |---- |     |
!         |  |  |     |
!         -------------
          IF (RefineX>RefineNachbarX) THEN
            DO jx=jx0,jx1,IncrX
              ix=jx/IncrX
              DO jy=jy0+1,jy1
                iy=IncrY*jy
                IF(ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Y Vert1-/Vert2%in_out==-1(1) sind, und RefineY groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Y Vert1-/Vert2%in_out==-1(1) sind, und RefineY groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          ELSE
            DO jx=jx0,jx1
              ix=jx*IncrX
              DO jy=jy0+1,jy1
                iy=IncrY*jy
                IF(ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Y(ix,iy,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Y Vert1-/Vert2%in_out==-1(1) sind, und RefineY groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
                IF(ASSOCIATED(Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge).AND.Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp==-1) THEN
                  IF(ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==1) THEN
                    IF (Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%Vert1%in_out &
                       *Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%Vert2%in_out<0) THEN
                      ALLOCATE(Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%VertS)
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%VertS=Edges_Y(jx,jy,jz)%Edge%VertS
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=1
                    ELSE
                      Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=0
                    END IF
                  ELSE IF (ASSOCIATED(Edges_Y(jx,jy,jz)%Edge).AND.Edges_Y(jx,jy,jz)%Edge%yes_sp==0) THEN
                     !Floor(ibn)%Edges_Y(ix,iy-1,iz)%Edge%yes_sp=0
                     ! bedeutet das Edges_Y Vert1-/Vert2%in_out==-1(1) sind, und RefineY groeber als Nachbar
                     !                                                 ,Nachbar kann trotzdem Schnitt haben!!
                     ! Special-Fall getestet
                  END IF
                END IF
              END DO
            END DO
          END IF
        END IF ! RefineY
      END IF  ! 'ib, it'

!........................................................................................

    END DO  ! AnzahlNachbar
  END DO   ! ib

END SUBROUTINE AnalyzeAllEdges

SUBROUTINE SortVertex
  INTEGER :: ib,i,j,k
  INTEGER :: nrP

  ALLOCATE(VertOut(1:nr_out))
  nr_EgVerts=0
  DO ib=1,nb
    CALL Set(Floor(ib))

    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          nrP=Vertices(i,j,k)%nrP
          IF (nrP>0) THEN
            VertOut(nrP)=Vertices(i,j,k)
          END IF
        END DO
      END DO
    END DO
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
            IF (Edges_X(i,j,k)%Edge%yes_sp==1) THEN
               nrP=Edges_X(i,j,k)%Edge%VertS%nrP
               IF (nrP>0) THEN
                 VertOut(nrP)=Edges_X(i,j,k)%Edge%VertS
                 nr_EgVerts=nr_EgVerts+1
               END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
            IF (Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
              nrP=Edges_Y(i,j,k)%Edge%VertS%nrP
              IF (nrP>0) THEN
                VertOut(nrP)=Edges_Y(i,j,k)%Edge%VertS
                nr_EgVerts=nr_EgVerts+1
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
            IF (Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
              nrP=Edges_Z(i,j,k)%Edge%VertS%nrP
              IF (nrP>0) THEN
                VertOut(nrP)=Edges_Z(i,j,k)%Edge%VertS
                nr_EgVerts=nr_EgVerts+1
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO

  END DO   ! ib

END SUBROUTINE SortVertex

SUBROUTINE SortVertexView
  INTEGER :: ib,i,j,k
  INTEGER :: nrP,iout,nr_set

  ALLOCATE(VertViewScale(1:nr_out))
  VertViewScale(1:nr_out)=0
  nr_view_out=0
  nr_viewEgVerts=0

  DO ib=1,nb
    CALL Set(Floor(ib))

    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF ((i>=Domain%view_ixa.AND.i<=Domain%view_ixe) .AND. &
              (j>=Domain%view_iya.AND.j<=Domain%view_iye) .AND. &
              (k>=Domain%view_iza.AND.k<=Domain%view_ize)) THEN
            nrP=Vertices(i,j,k)%nrP
            IF (nrP>0) THEN
              IF (VertViewScale(nrP)==0) THEN
                nr_view_out=nr_view_out+1        ! noch Doppelbelegung Zylinder
                VertViewScale(nrP)=nr_view_out
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF ((i>=Domain%view_ixa+1.AND.i<=Domain%view_ixe) .AND. &
              (j>=Domain%view_iya.AND.j<=Domain%view_iye) .AND. &
              (k>=Domain%view_iza.AND.k<=Domain%view_ize)) THEN
            IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
              IF (Edges_X(i,j,k)%Edge%yes_sp==1) THEN
                 nrP=Edges_X(i,j,k)%Edge%VertS%nrP
                 IF (nrP>0) THEN
                   IF (VertViewScale(nrP)==0) THEN
                     nr_view_out=nr_view_out+1
                     VertViewScale(nrP)=nr_view_out
                     nr_viewEgVerts=nr_viewEgVerts+1
                   END IF
                 END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF ((i>=Domain%view_ixa.AND.i<=Domain%view_ixe) .AND. &
              (j>=Domain%view_iya+1.AND.j<=Domain%view_iye) .AND. &
              (k>=Domain%view_iza.AND.k<=Domain%view_ize)) THEN
            IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
              IF (Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
                nrP=Edges_Y(i,j,k)%Edge%VertS%nrP
                IF (nrP>0)  THEN
                   IF (VertViewScale(nrP)==0) THEN
                      nr_view_out=nr_view_out+1
                      VertViewScale(nrP)=nr_view_out
                      nr_viewEgVerts=nr_viewEgVerts+1
                   END IF
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF ((i>=Domain%view_ixa.AND.i<=Domain%view_ixe) .AND. &
              (j>=Domain%view_iya.AND.j<=Domain%view_iye) .AND. &
              (k>=Domain%view_iza+1.AND.k<=Domain%view_ize)) THEN
            IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
              IF (Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
                nrP=Edges_Z(i,j,k)%Edge%VertS%nrP
                IF (nrP>0) THEN
                  IF (VertViewScale(nrP)==0) THEN
                    nr_view_out=nr_view_out+1
                    VertViewScale(nrP)=nr_view_out
                    nr_viewEgVerts=nr_viewEgVerts+1
                  END IF
                END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO

  END DO   ! ib

  ALLOCATE(VertOutView(1:nr_view_out))
  nr_set=0
  DO iout=1,nr_out
   IF(VertViewScale(iout)>0) THEN
     nr_set=nr_set+1
     VertOutView(VertViewScale(iout))=VertOut(iout)
   END IF
  END DO
  IF(nr_set/=nr_view_out) THEN
     Write(OutUnitProt,*) "---SortVertexView---,nr_set/=nr_view_out"
  END IF
END SUBROUTINE SortVertexView


SUBROUTINE SortVertexOut2
  INTEGER :: ib,i,j,k
  INTEGER :: ioro

  ioro=nr_out
  !Vert-Zählung unterhalb Berg
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF (Vertices(i,j,k)%in_out==-1) THEN
            ioro=ioro+1 
            Vertices(i,j,k)%nrI=ioro
          END IF
        END DO
      END DO
    END DO
  END DO   ! ib
  nr_out2=ioro
  ALLOCATE(VertOut2(nr_out+1:nr_out2))

  DO ib=1,nb
    CALL Set(Floor(ib))
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          ioro=Vertices(i,j,k)%nrI
          IF (ioro>0) THEN
            VertOut2(ioro)=Vertices(i,j,k)
          END IF
        END DO
      END DO
    END DO
  END DO   ! ib

END SUBROUTINE SortVertexOut2


SUBROUTINE SortVertexOro
  INTEGER :: ib,i,j,k
  INTEGER :: nrInP

  ALLOCATE(VertIn(1:nr_inside))

  DO ib=1,nb
    CALL Set(Floor(ib))

    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          nrInP=Vertices(i,j,k)%nrInP
          IF (nrInP>0) THEN
            VertIn(nrInP)=Vertices(i,j,k)
          END IF
        END DO
      END DO
    END DO
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
            IF (Edges_X(i,j,k)%Edge%yes_sp==1) THEN
               nrInP=Edges_X(i,j,k)%Edge%VertS%nrInP
               IF (nrInP>0) THEN
                 VertIn(nrInP)=Edges_X(i,j,k)%Edge%VertS 
               END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
            IF (Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
              nrInP=Edges_Y(i,j,k)%Edge%VertS%nrInP
              IF (nrInP>0) THEN
                VertIn(nrInP)=Edges_Y(i,j,k)%Edge%VertS
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
            IF (Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
              nrInP=Edges_Z(i,j,k)%Edge%VertS%nrInP
              IF (nrInP>0) THEN
                VertIn(nrInP)=Edges_Z(i,j,k)%Edge%VertS
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO

  END DO   ! ib

END SUBROUTINE SortVertexOro


SUBROUTINE SortVertexCutPlane
  INTEGER :: ib,i,j,k
  INTEGER :: nrP
  INTEGER :: nrCutP

  !ALLOCATE(VertOut(1:nr_out))
  ALLOCATE(VertCutPOut(1:nr_cutplane))

  DO ib=1,nb
    CALL Set(Floor(ib))

    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          nrP=Vertices(i,j,k)%nrP
          nrCutP=Vertices(i,j,k)%nrCutP
          IF (nrCutP>0) THEN
            VertCutPOut(nrCutP)=Vertices(i,j,k)  !fuer Paramertic
          END IF
        END DO
      END DO
    END DO
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
            IF (Edges_X(i,j,k)%Edge%yes_sp==1) THEN
               nrP=Edges_X(i,j,k)%Edge%VertS%nrP
               nrCutP=Edges_X(i,j,k)%Edge%VertS%nrCutP
               IF (nrCutP>0) THEN
                  VertCutPOut(nrCutP)=Edges_X(i,j,k)%Edge%VertS
               END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
            IF (Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
              nrP=Edges_Y(i,j,k)%Edge%VertS%nrP
              nrCutP=Edges_Y(i,j,k)%Edge%VertS%nrCutP
              IF (nrCutP>0) THEN
                 VertCutPOut(nrCutP)=Edges_Y(i,j,k)%Edge%VertS
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
            IF (Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
              nrP=Edges_Z(i,j,k)%Edge%VertS%nrP
              nrCutP=Edges_Z(i,j,k)%Edge%VertS%nrCutP
              IF (nrCutP>0) THEN
                VertCutPOut(nrCutP)=Edges_Z(i,j,k)%Edge%VertS
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO

  END DO   ! ib

END SUBROUTINE SortVertexCutPlane


SUBROUTINE SortVertexSoil
  INTEGER :: ib,i,j,k
  INTEGER :: nrP
  INTEGER :: nrCutP
  TYPE (Vertex_T) :: V(1:2)
 
  IF(nr_sb>0) THEN
     !ALLOCATE(VertCutPOut(1:nr_cutplane))
     nr_out_soil=(1+nr_lsoil)*nr_cutplane
     ALLOCATE(VertSoilOut(1:nr_out_soil)) 
   
     !VertSoilOut(1:nr_cutplane)=VertCutPOut(1:nr_cutplane)
     
     DO ib=1,nb
       CALL Set(Floor(ib))
   
       DO i=ix0,ix1
         DO j=iy0,iy1
           DO k=iz0,iz1
             nrP=Vertices(i,j,k)%nrP
             nrCutP=Vertices(i,j,k)%nrCutP
             IF (nrCutP>0) THEN
               !VertCutPOut(nrCutP)=Vertices(i,j,k)  !fuer Paramertic
               VertSoilOut(nrCutP)=Vertices(i,j,k)
             END IF
           END DO
         END DO
       END DO
       DO i=ix0+1,ix1
         DO j=iy0,iy1
           DO k=iz0,iz1
             IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
               IF (Edges_X(i,j,k)%Edge%yes_sp==1) THEN
                  nrP=Edges_X(i,j,k)%Edge%VertS%nrP
                  nrCutP=Edges_X(i,j,k)%Edge%VertS%nrCutP
                  IF (nrCutP>0) THEN
                     VertSoilOut(nrCutP)=Edges_X(i,j,k)%Edge%VertS
                  END IF
               END IF
             END IF
           END DO
         END DO
       END DO
       DO i=ix0,ix1
         DO j=iy0+1,iy1
           DO k=iz0,iz1
             IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
               IF (Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
                 nrP=Edges_Y(i,j,k)%Edge%VertS%nrP
                 nrCutP=Edges_Y(i,j,k)%Edge%VertS%nrCutP
                 IF (nrCutP>0) THEN
                    VertSoilOut(nrCutP)=Edges_Y(i,j,k)%Edge%VertS
                 END IF
               END IF
             END IF
           END DO
         END DO
       END DO
       DO i=ix0,ix1
         DO j=iy0,iy1
           DO k=iz0+1,iz1
             IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
               IF (Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
                 nrP=Edges_Z(i,j,k)%Edge%VertS%nrP
                 nrCutP=Edges_Z(i,j,k)%Edge%VertS%nrCutP
                 IF (nrCutP>0) THEN
                   VertSoilOut(nrCutP)=Edges_Z(i,j,k)%Edge%VertS
                 END IF
               END IF
             END IF
           END DO
         END DO
       END DO
   
     END DO   ! ib
     !Write(*,*) "Domain%dz(1)=",Domain%dz(1),"   aus SortVertexSoil"
     !dzi_soil=zParametricOut(Domain%Vertices(i-1,j-1,iz0+1))-zParametricOut(Vertices(i-1,j-1,iz0))
     V(1:2)%Point%x=Domain%xP(0)
     V(1:2)%Point%y=Domain%yP(0)
     V(1)%Point%z=Domain%zP(0)
     V(2)%Point%z=Domain%zP(1)
     dzi_soil=zParametricOut(V(2))-zParametricOut(V(1))
     !Write(*,*) "dzi_soil=",dzi_soil,"   aus SortVertexSoil"
     DO i=2,nr_lsoil+1
       VertSoilOut(((i-1)*nr_cutplane)+1:(i*nr_cutplane))=VertCutPOut(1:nr_cutplane)
       VertSoilOut(((i-1)*nr_cutplane)+1:(i*nr_cutplane))%Point%z= &
          VertSoilOut(((i-1)*nr_cutplane)+1:(i*nr_cutplane))%Point%z-(i-1)*dzi_soil !dz(1) !dz(1)/2 !200.0d0 !dz(1)/7
       VertSoilOut(((i-1)*nr_cutplane)+1:(i*nr_cutplane))%nrCutP= &
          VertSoilOut(((i-1)*nr_cutplane)+1:(i*nr_cutplane))%nrCutP+((i-1)*nr_cutplane)
     END DO 
   END IF ! if(nr_sb>0)

END SUBROUTINE SortVertexSoil


FUNCTION VolFace_XY(Face)
   REAL(8) :: VolFace_XY
   TYPE (Face_T) :: Face
   REAL(8) :: r,phi1,phi2,lam1,lam2,r1,r2
   !phi-> [-Pi/2,Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. &
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
       indata_type=='geo' .OR. indata_type=='rad' .OR. &
       indata_type=='wgs84') THEN
     phi1=Face%Edge2%Vert1%Point%y
     phi2=Face%Edge2%Vert2%Point%y
     lam1=Face%Edge1%Vert1%Point%x
     lam2=Face%Edge1%Vert2%Point%x
     r=RadEarth+Face%Edge1%Vert1%Point%z
     VolFace_XY=r**2*(SIN(phi2)-SIN(phi1))*(lam2-lam1)
   ELSE If (indata_type=='cyl') THEN
     lam1=Face%Edge1%Vert1%Point%x
     lam2=Face%Edge1%Vert2%Point%x
     r1=Face%Edge2%Vert1%Point%y
     r2=Face%Edge2%Vert2%Point%y
     VolFace_XY=ABS((r2+r1)*(r2-r1)*(lam2-lam1)/2.0d0)
   ELSE  ! ('Cart',DEFAULT)
     VolFace_XY=ABS((Face%Edge1%Vert2%Point%x-Face%Edge1%Vert1%Point%x) &
                       *(Face%Edge2%Vert2%Point%y-Face%Edge2%Vert1%Point%y))
   END IF
END FUNCTION VolFace_XY

FUNCTION VolFace_YZ(Face)
   REAL(8) :: VolFace_YZ
   TYPE (Face_T) :: Face
   REAL(8) :: phi1,phi2,r1,r2
   !phi-> [-Pi/2;Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. & 
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
       indata_type=='geo' .OR. indata_type=='rad' .OR. &
       indata_type=='wgs84') THEN
     phi1=Face%Edge1%Vert1%Point%y
     phi2=Face%Edge1%Vert2%Point%y
     r1  =RadEarth+Face%Edge2%Vert1%Point%z
     r2  =RadEarth+Face%Edge2%Vert2%Point%z 
     VolFace_YZ=0.5*(phi2-phi1)*(r2**2-r1**2)
   ELSE If (indata_type=='cyl') THEN
     VolFace_YZ=ABS((Face%Edge1%Vert2%Point%y-Face%Edge1%Vert1%Point%y) &
                   *(Face%Edge2%Vert2%Point%z-Face%Edge2%Vert1%Point%z))
   ELSE !('Cart',DEFAULT)
     VolFace_YZ=ABS((Face%Edge1%Vert2%Point%y-Face%Edge1%Vert1%Point%y) &
                   *(Face%Edge2%Vert2%Point%z-Face%Edge2%Vert1%Point%z))
   END IF
END FUNCTION VolFace_YZ

FUNCTION VolFace_ZX(Face)
   REAL(8) :: VolFace_ZX
   TYPE (Face_T) :: Face
   REAL(8) :: phi,lam1,lam2,r1,r2
   !phi-> [-Pi/2,Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. &
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
      indata_type=='geo' .OR. indata_type=='rad' .OR. &
      indata_type=='wgs84') THEN
     phi =Face%Edge1%Vert1%Point%y
     lam1=Face%Edge2%Vert1%Point%x
     lam2=Face%Edge2%Vert2%Point%x
     r1  =RadEarth+Face%Edge1%Vert1%Point%z
     r2  =RadEarth+Face%Edge1%Vert2%Point%z
     VolFace_ZX=0.5*COS(phi)*(lam2-lam1)*(r2**2-r1**2)
   ELSE If (indata_type=='cyl') THEN
     r2 =Face%Edge1%Vert2%Point%y
     VolFace_ZX=ABS((Face%Edge1%Vert2%Point%z-Face%Edge1%Vert1%Point%z) &
                     *(Face%Edge2%Vert2%Point%x-Face%Edge2%Vert1%Point%x)&
                     *r2)
   ELSE !('Cart',DEFAULT)
     VolFace_ZX=ABS((Face%Edge1%Vert2%Point%z-Face%Edge1%Vert1%Point%z) &
                     *(Face%Edge2%Vert2%Point%x-Face%Edge2%Vert1%Point%x))
   END IF

END FUNCTION VolFace_ZX


SUBROUTINE InitAllFaces
  INTEGER :: ib,i,j,k
  REAl(8) :: xv1,yv1,zv1,xv2,yv2,zv2,Vol

  DO ib=1,nb
    CALL Set(Floor(ib))

  END DO  ! ib
END SUBROUTINE InitAllFaces


SUBROUTINE AnalyzeFace(Face,NrW,NrMP)
  TYPE (Face_T) :: Face
  INTEGER :: NrW,NrMP

  INTEGER :: i,j,v,w,ecut,iTemp,tdiff
  INTEGER :: ic,jc,kc,nxc,nyc,nzc
  REAL(8) :: xc0,xc1,yc0,yc1,zc0,zc1,dx,dy,dz
  TYPE (Vertex_T) :: InterVert
  TYPE (Vertex_T) :: Vertex1,Vertex2
  TYPE(Point_T) :: EdgP0,EdgP1,EdgP2,EdgP3
  TYPE(Point_T) :: P0,P1,P2,P3
  TYPE(Point_T) :: PMin,PMax,SumMidP,MidPoint
  REAL(8) :: xvc0,xvc1,xvc0u,xvc0o,xvc1u,xvc1o
  REAL(8) :: yvc0,yvc1,yvc0u,yvc0o,yvc1u,yvc1o
  REAL(8) :: zvc0,zvc1
  REAL(8) :: Vol,Len
  REAL(8) :: xsi,ys
  INTEGER :: n_egl
  INTEGER :: EdgeList(1:2,1:6) 

  Face%ec=-1
  Face%NumberVert=0
  Face%VertexList=0
  Face%EdgeCut(:)=0
  i=0
  ecut=0
  Face%in_out=Face%Edge1%in_out+Face%Edge3%in_out
  !...............................................
  IF (Face%Edge1%Vert1%in_out>=0) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge1%Vert1%nrP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge1%yes_sp==1) THEN
!D    IF(((Face%Edge1%Vert1%Point.PEQ.Face%Edge1%VertS%Point)==1) .OR. &
!D       ((Face%Edge1%Vert2%Point.PEQ.Face%Edge1%VertS%Point)==1) ) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge1%VertS%nrP
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge1%VertS%nrP
!D    END IF
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%Edge1%VertS%nrP
  ELSE IF (Face%Edge1%Vert1%in_out==0) THEN
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge1%Vert1%nrP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge1%Vert1%nrP
!D    END IF
  END IF
  !............................................
  IF (Face%Edge2%Vert1%in_out>=0) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge2%Vert1%nrP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge2%yes_sp==1) THEN
!D    IF(((Face%Edge2%Vert1%Point.PEQ.Face%Edge2%VertS%Point)==1) .OR. &
!D       ((Face%Edge2%Vert2%Point.PEQ.Face%Edge2%VertS%Point)==1) ) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge2%VertS%nrP
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge2%VertS%nrP
!D    END IF
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%Edge2%VertS%nrP
  ELSE IF (Face%Edge2%Vert1%in_out==0) THEN
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge2%Vert1%nrP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge2%Vert1%nrP
  END IF
  !..............................................
  IF (Face%Edge3%Vert2%in_out>=0) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge3%Vert2%nrP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge3%yes_sp==1) THEN
!D    IF(((Face%Edge3%Vert2%Point.PEQ.Face%Edge3%VertS%Point)==1) .OR. &
!D       ((Face%Edge3%Vert1%Point.PEQ.Face%Edge3%VertS%Point)==1) ) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge3%VertS%nrP
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge3%VertS%nrP
!D    END IF
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%Edge3%VertS%nrP
  ELSE IF (Face%Edge3%Vert2%in_out==0) THEN
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge3%Vert2%nrP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge3%Vert2%nrP
  END IF
  !................................................
  IF (Face%Edge4%Vert2%in_out>=0) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge4%Vert2%nrP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge4%yes_sp==1) THEN
!D    IF(((Face%Edge4%Vert1%Point.PEQ.Face%Edge4%VertS%Point)==1) .OR. &
!D       ((Face%Edge4%Vert2%Point.PEQ.Face%Edge4%VertS%Point)==1) )THEN
    i=i+1
    Face%VertexList(i)=Face%Edge4%VertS%nrP
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge4%VertS%nrP
      !IF(ecut<3) Face%EdgeCut(ecut)=Face%Edge4%VertS%nrP
!D    END IF
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%Edge4%VertS%nrP
  ELSE IF (Face%Edge4%Vert2%in_out==0) THEN
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge4%Vert2%nrP
    !IF(ecut<3) Face%EdgeCut(ecut)=Face%Edge4%Vert2%nrP
    Face%ec=1
    EdgeList(1,ecut)=0
    EdgeList(2,ecut)=Face%Edge4%Vert2%nrP
  END IF
  !Cut neu sortieren wenn ecut>2 ermittelt
  !Verts mit in_out==0 am Vert1||Vert2 filtern
  !richtiges ecut setzen, 
  !(aber dennoch Frage/offen, '1 x yes_sp und mehrfach Vert!%in_out==0'
  !je nachdem letzten 2 erfassten Cut in Liste? )
  IF(ecut>2) THEN
     j=0
     DO i=1,ecut
       IF(EdgeList(1,i)==1) THEN
         j=j+1
       END IF
     END DO 
     IF(j==2) THEN
       v=0
       Face%ec=1
       DO i=1,ecut
         IF(EdgeList(1,i)==1) THEN
           v=v+1
           Face%EdgeCut(v)=EdgeList(2,i)
         END IF
       END DO
     END IF
  END IF
  ! -------
  ! Special
  ! -------
  !IF (Face%Edge1%Vert1%in_out>=0 .AND. Face%Edge3%Vert2%in_out>=0 .AND. &
  !   Face%Edge1%Vert2%in_out<=0 .AND. Face%Edge3%Vert1%in_out<=0) THEN
  !   Face%ec=-1
  !   Face%EdgeCut(:)=0
  !ELSE IF (Face%Edge1%Vert1%in_out>=0 .AND. Face%Edge3%Vert2%in_out>=0 .AND. &
  !   Face%Edge1%Vert2%in_out<=0 .AND. Face%Edge3%Vert1%in_out<=0) THEN
  !   Face%ec=-1
  !   Face%EdgeCut(:)=0
  !  ! filtert auch  (1) ---*-- (-1)
  !  ! Schnitte:         |    *
  !                      |    |
  !                  (0) ------  (1)
  !END IF

  IF (Face%Edge1%Vert1%in_out<0 .AND. Face%Edge3%Vert2%in_out<0 .AND. &
     Face%Edge1%Vert2%in_out>0 .AND. Face%Edge3%Vert1%in_out>0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  ELSE IF (Face%Edge1%Vert1%in_out>0 .AND. Face%Edge3%Vert2%in_out>0 .AND. &
     Face%Edge1%Vert2%in_out<0 .AND. Face%Edge3%Vert1%in_out<0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  END IF
  !
  ! Face mit 1 Schnittpunkt -->EdgeCut ausblenden
  IF (ecut==1) THEN     ! zu Allgenein !?
    Face%ec=-1
  END IF

  !................................................
  ! Face-Volumen-Berechnung
  !IF (Face%in_out<=-4.OR.Face%ec==-1.OR. &
  IF (Face%Edge1%Vert1%in_out==0.AND.Face%Edge1%Vert2%in_out==0.AND. &
      Face%Edge3%Vert1%in_out==0.AND.Face%Edge3%Vert2%in_out==0) THEN
        Face%Vol=0.0d0
        EdgP0=Face%Edge1%Vert1%Point !Ossi
        EdgP1=Face%Edge1%Vert2%Point
        EdgP2=Face%Edge3%Vert2%Point
        EdgP3=Face%Edge3%Vert1%Point
        P0=PointParametricEarth(Face%Edge1%Vert1%Point)
        P1=PointParametricEarth(Face%Edge1%Vert2%Point)
        P2=PointParametricEarth(Face%Edge3%Vert2%Point)
        P3=PointParametricEarth(Face%Edge3%Vert1%Point)
        Face%MidPoint=FaceMidP4(P0,P1,P2,P3)
        Face%MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) ! Ossi
        !SumMidP=SumMidP+Vol*MidPoint
        !Face%MidPoint=SumMidP/Face%Vol
        NrMP=NrMP+1
        Face%mp=NrMP
  !Betreff Bsp.Agi_2Blk_H5000_test.grid in VGrid_GMV_8.2.1.1
  !Dürfte wenn in AnalyzeEdge ausgklinkt nicht nötig sein
  !In AnalyzeAllEdge sollche Flächen nicht allokieren
  !ELSE IF (Face%NumberVert==2.AND.ecut==2) THEN  ! Edge=Grenzfläche
  !     Face%Vol=0.0d0
  ELSE IF (Face%in_out==-4.OR. &
      MAX(Face%Edge1%Vert1%in_out,Face%Edge1%Vert2%in_out, &
          Face%Edge3%Vert1%in_out,Face%Edge3%Vert2%in_out)<=0 &
     ) THEN
    Face%Vol=0.0d0
  ELSE IF (Face%in_out==2.AND. &
           MIN(Face%Edge1%Vert1%in_out,Face%Edge1%Vert2%in_out, &
               Face%Edge3%Vert1%in_out,Face%Edge3%Vert2%in_out)==0 &
          ) THEN
    Face%ec=0   !! 1Edge Grenze dann Face=maxVol
  ELSE IF (Face%in_out==1.AND. &
           MIN(Face%Edge1%Vert1%in_out,Face%Edge1%Vert2%in_out, &
               Face%Edge3%Vert1%in_out,Face%Edge3%Vert2%in_out)==0 &
          ) THEN
    Face%ec=0   !! 2Edge Grenze dann Face=maxVol
  ELSE IF (Face%ec>0) THEN
    PMin=Face%Edge2%Vert2%Point
    PMax=Face%Edge1%Vert1%Point
    CALL SearchMinMaxFace(Face,PMin,PMax)
    Face%Vol=0.0d0
    nxc=IncrVol
    nyc=IncrVol
    nzc=IncrVol
    xc0=Face%Edge1%Vert1%Point%x
    yc0=Face%Edge1%Vert1%Point%y
    zc0=Face%Edge1%Vert1%Point%z
    xc1=Face%Edge2%Vert2%Point%x
    yc1=Face%Edge2%Vert2%Point%y
    zc1=Face%Edge2%Vert2%Point%z
    !        3-----2   ! generell Point-Folge Face-Vol-Berechnung
    !        | \   |
    !        |   \ |
    !        0-----1
    !---Face_YZ------------------------------------------------------------
    IF (xc0==xc1) THEN
      IF ((Face%Edge4%in_out==2).AND.(Face%Edge2%in_out==-2) .OR. &
          (Face%Edge4%in_out==-2).AND.(Face%Edge2%in_out==2))THEN
         !.........y-direction.............................................
         !     0 ---------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !     |          |   
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc0
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%z-PMin%z
         zvc0=PMin%z+Shrink*Len
         zvc1=PMax%z-Shrink*Len
         dz=(zvc1-zvc0)/nzc
         DO kc=1,nzc
           Vertex1%Point%z=zvc0
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeYaxis(EdgP0,EdgP3,Vertex1,Vertex2,yc0,yc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3)  !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3)  !nur je Teilflaeche
           SumMidP=SumMidP+Vol*MidPoint
           zvc0=MIN(zvc0+dz,zvc1)
         END DO
      ELSE
         !...........z-direction................................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc0
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%y-PMin%y
         yvc0=PMin%y+Shrink*Len
         yvc1=PMax%y-Shrink*Len
         dy=(yvc1-yvc0)/nyc
         DO jc=1,nyc
           Vertex1%Point%y=yvc0
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP0,EdgP3,Vertex1,Vertex2,zc0,zc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3)  !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3)  !nur je Teilflaeche
           SumMidP=SumMidP+Vol*MidPoint
           yvc0=MIN(yvc0+dy,yvc1)
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    !---Face_XZ------------------------------------------------------------
    ELSE IF (yc0==yc1) THEN
      IF ((Face%Edge1%in_out==2).AND.(Face%Edge3%in_out==-2) .OR. &
          (Face%Edge1%in_out==-2).AND.(Face%Edge3%in_out==2))THEN
         !........x-direction.......................................... 
         !     0----------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !     |          |
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         zvc0=PMin%z+Shrink*Len
         zvc1=PMax%z-Shrink*Len
         dz=(zvc1-zvc0)/nzc
         DO kc=1,nzc
           Vertex1%Point%z=zvc0
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeXaxis(EdgP0,EdgP3,Vertex1,Vertex2,xc0,xc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           zvc0=zvc0+dz
         END DO
      ELSE
         !..........z-direction..............................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc0
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         xvc0=PMin%x+Shrink*Len
         xvc1=PMax%x-Shrink*Len
         dx=(xvc1-xvc0)/nxc
         DO ic=1,nxc
           Vertex1%Point%x=xvc0
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeZaxis(EdgP0,EdgP3,Vertex1,Vertex2,zc0,zc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           xvc0=xvc0+dx
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    !---Face_XY------------------------------------------------------------
    ELSE IF (zc0==zc1) THEN
      IF ((Face%Edge4%in_out==2).AND.(Face%Edge2%in_out==-2) .OR. &
          (Face%Edge4%in_out==-2).AND.(Face%Edge2%in_out==2))THEN
         !....x-direction..................................................
         !     0----------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         ! dy  |          |
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc1
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%y-PMin%y
         yvc0=PMin%y+Shrink*Len
         yvc1=PMax%y-Shrink*Len
         dy=(yvc1-yvc0)/nyc
         DO jc=1,nyc
           Vertex1%Point%y=yvc0
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP0,EdgP3,Vertex1,Vertex2,xc0,xc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           yvc0=yvc0+dy
         END DO
      ELSE
         !....y-direction...................................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc1
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         xvc0=PMin%x+Shrink*Len
         xvc1=PMax%x-Shrink*Len
         dx=(xvc1-xvc0)/nxc
         DO ic=1,nxc
           Vertex1%Point%x=xvc0
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP0,EdgP3,Vertex1,Vertex2,yc0,yc1)
           P0=PointParametricEarth(EdgP0)
           P3=PointParametricEarth(EdgP3)
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           MidPoint=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           xvc0=xvc0+dx
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    END IF !!(xc0==xc1);(yc0==yc1);(zc0==zc1)
  END IF  !! ec

  !IF ((Face%in_out<4.AND.Face%ec>0).OR.(Face%in_out<=0.AND.Face%ec==-1)) THEN
  ! unter 7.7.1.3 eingearbeitet !! checken
  IF ((Face%in_out<4.AND.Face%ec>0).OR.(Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
    NrW=NrW+1
  END IF

END SUBROUTINE AnalyzeFace


SUBROUTINE AnalyzeAllFaces

  INTEGER :: ib,i,j,k
  INTEGER :: in_out

  INTEGER :: NrF_YZ_ASS,NrF_ZX_ASS,NrF_XY_ASS
  INTEGER :: NrF_YZ_UG, NrF_ZX_UG, NrF_XY_UG

  NrW_All_FXY=0; NrMP_All_FXY=0 
  NrW_All_FYZ=0; NrMP_All_FYZ=0
  NrW_All_FZX=0; NrMP_All_FZX=0
  
  DO ib=1,nb

    CALL Set(Floor(ib))
    NrF_YZ_ASS=0; NrF_ZX_ASS=0; NrF_XY_ASS=0
    NrF_YZ_UG=0;  NrF_ZX_UG=0;  NrF_XY_UG=0

    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
            !CALL AnalyzeFace(Faces_XY(i,j,k)%Face,NrW_FacesXY,NrMP_FacesXY)
            !Write(OutUnitProt,*)"Faces_XY...."
            !CALL WriteFaceEVProt(Faces_XY(i,j,k)%Face,i,j,k)
            CALL AnalyzeFace(Faces_XY(i,j,k)%Face,NrF_XY_ASS,NrMP_FacesXY)
          ELSE
            in_out=Vertices(i-1,j-1,k)%in_out &
                  +Vertices(i-1,j,k)%in_out &
                  +Vertices(i,j-1,k)%in_out &
                  +Vertices(i,j,k)%in_out
            IF (in_out<=0) THEN
              !NrW_FacesXY=NrW_FacesXY+1
              NrF_XY_UG=NrF_XY_UG+1
            END IF
          END IF
        END DO
      END DO
    END DO
    NrW_FacesXY=NrF_XY_ASS+NrF_XY_UG
    Floor(ib)%NrW_FacesXY=NrW_FacesXY
    Floor(ib)%NrMP_FacesXY=NrMP_FacesXY
    NrW_All_FXY=NrW_All_FXY+NrW_FacesXY
    NrMP_All_FXY=NrMP_All_FXY+NrMP_FacesXY

    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
            !CALL AnalyzeFace(Faces_YZ(i,j,k)%Face,NrW_FacesYZ,NrMP_FacesYZ)
            !Write(OutUnitProt,*)"Faces_YZ...."
            !CALL WriteFaceEVProt(Faces_YZ(i,j,k)%Face,i,j,k)
            CALL AnalyzeFace(Faces_YZ(i,j,k)%Face,NrF_YZ_ASS,NrMP_FacesYZ)
          ELSE
            in_out=Vertices(i,j-1,k-1)%in_out &
                  +Vertices(i,j,k-1)%in_out &
                  +Vertices(i,j-1,k)%in_out &
                  +Vertices(i,j,k)%in_out
            IF (in_out<=0) THEN
              !NrW_FacesYZ=NrW_FacesYZ+1
              NrF_YZ_UG=NrF_YZ_UG+1
            END IF
          END IF
        END DO
      END DO
    END DO
    NrW_FacesYZ=NrF_YZ_ASS+NrF_YZ_UG
    Floor(ib)%NrW_FacesYZ=NrW_FacesYZ
    Floor(ib)%NrMP_FacesYZ=NrMP_FacesYZ
    NrW_All_FYZ=NrW_All_FYZ+NrW_FacesYZ
    NrMP_All_FYZ=NrMP_All_FYZ+NrMP_FacesYZ
 
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
            !CALL AnalyzeFace(Faces_ZX(i,j,k)%Face,NrW_FacesZX,NrMP_FacesZX)
            !Write(OutUnitProt,*)"Faces_ZX...."
            !CALL WriteFaceEVProt(Faces_ZX(i,j,k)%Face,i,j,k)
            CALL AnalyzeFace(Faces_ZX(i,j,k)%Face,NrF_ZX_ASS,NrMP_FacesZX)
          ELSE
            in_out=Vertices(i-1,j,k-1)%in_out &
                  +Vertices(i,j,k-1)%in_out &
                  +Vertices(i-1,j,k)%in_out &
                  +Vertices(i,j,k)%in_out
            IF (in_out<=0) THEN
              !NrW_FacesZX=NrW_FacesZX+1
              NrF_ZX_UG=NrF_ZX_UG+1
            END IF
          END IF
        END DO
      END DO
    END DO
    NrW_FacesZX=NrF_ZX_ASS+NrF_ZX_UG
    Floor(ib)%NrW_FacesZX=NrW_FacesZX
    Floor(ib)%NrMP_FacesZX=NrMP_FacesZX
    NrW_All_FZX=NrW_All_FZX+NrW_FacesZX
    NrMP_All_FZX=NrMP_All_FZX+NrMP_FacesZX
  
  END DO   ! ib
END SUBROUTINE AnalyzeAllFaces

!.................................................................
! Rand Face
!...........

SUBROUTINE SearchMinMaxP4Face(Face,PMin,PMax)
  TYPE(FaceP4_T) :: Face 
  TYPE(Point_T) :: PMin,PMax
  TYPE(Point_T) :: P0,P1,P2,P3
   
  IF(Face%V0%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Face%V0%Point%x)
    PMin%y=MIN(PMin%y,Face%V0%Point%y)
    PMin%z=MIN(PMin%z,Face%V0%Point%z)
    PMax%x=MAX(PMax%x,Face%V0%Point%x)
    PMax%y=MAX(PMax%y,Face%V0%Point%y)
    PMax%z=MAX(PMax%z,Face%V0%Point%z)
  END IF
  IF(Face%V1%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Face%V1%Point%x)
    PMin%y=MIN(PMin%y,Face%V1%Point%y)
    PMin%z=MIN(PMin%z,Face%V1%Point%z)
    PMax%x=MAX(PMax%x,Face%V1%Point%x)
    PMax%y=MAX(PMax%y,Face%V1%Point%y)
    PMax%z=MAX(PMax%z,Face%V1%Point%z)
  END IF
  IF(Face%V2%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Face%V2%Point%x)
    PMin%y=MIN(PMin%y,Face%V2%Point%y)
    PMin%z=MIN(PMin%z,Face%V2%Point%z)
    PMax%x=MAX(PMax%x,Face%V2%Point%x)
    PMax%y=MAX(PMax%y,Face%V2%Point%y)
    PMax%z=MAX(PMax%z,Face%V2%Point%z)
  END IF
  IF(Face%V3%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Face%V3%Point%x)
    PMin%y=MIN(PMin%y,Face%V3%Point%y)
    PMin%z=MIN(PMin%z,Face%V3%Point%z)
    PMax%x=MAX(PMax%x,Face%V3%Point%x)
    PMax%y=MAX(PMax%y,Face%V3%Point%y)
    PMax%z=MAX(PMax%z,Face%V3%Point%z)
  END IF

  !................................
  ! Analyze Schnittpunkte
  ! EdgeX  -XY-Face und ZX -Face

  IF(Face%V0%Point%x/=Face%V1%Point%x) THEN
    CALL SearchPointsEdgeXaxis(P0,P1,Face%V0,Face%V1, &
                               Face%V0%Point%x,Face%V1%Point%x)
    CALL SearchPointsEdgeXaxis(P2,P3,Face%V2,Face%V3, &
                               Face%V2%Point%x,Face%V3%Point%x)
    IF(P0%x/=P1%x.AND.P0%x/=Face%V0%Point%x) THEN
      PMin%x=MIN(PMin%x,P0%x)
      PMax%x=MAX(PMax%x,P0%x)
    END IF
    IF(P0%x/=P1%x.AND.P1%x/=Face%V1%Point%x) THEN
      PMin%x=MIN(PMin%x,P1%x)
      PMax%x=MAX(PMax%x,P1%x)
    END IF
  
    IF(P2%x/=P3%x.AND.P2%x/=Face%V2%Point%x) THEN
      PMin%x=MIN(PMin%x,P2%x)
      PMax%x=MAX(PMax%x,P2%x)
    END IF
    IF(P2%x/=P3%x.AND.P3%x/=Face%V3%Point%x) THEN
      PMin%x=MIN(PMin%x,P3%x)
      PMax%x=MAX(PMax%x,P3%x)
    END IF
  END IF

  ! EdgeY -XY-Face und ZY-Face
  IF(Face%V0%Point%y/=Face%V2%Point%y.OR. &
     Face%V0%Point%y/=Face%V1%Point%y) THEN
    CALL SearchPointsEdgeYaxis(P0,P2,Face%V0,Face%V2, &
                               Face%V0%Point%y,Face%V2%Point%y)
    CALL SearchPointsEdgeYaxis(P1,P3,Face%V1,Face%V3, &
                               Face%V1%Point%y,Face%V3%Point%y)

    IF(P0%y/=P2%y.AND.P0%y/=Face%V0%Point%y) THEN
      PMin%y=MIN(PMin%y,P0%y)
      PMax%y=MAX(PMax%y,P0%y)
    END IF
    IF(P0%y/=P2%y.AND.P2%y/=Face%V2%Point%y) THEN
      PMin%y=MIN(PMin%y,P2%y)
      PMax%y=MAX(PMax%y,P2%y)
    END IF
  
    IF(P1%y/=P3%y.AND.P1%y/=Face%V1%Point%y) THEN
      PMin%y=MIN(PMin%y,P1%y)
      PMax%y=MAX(PMax%y,P1%y)
    END IF
    IF(P1%y/=P3%y.AND.P3%y/=Face%V3%Point%y) THEN
      PMin%y=MIN(PMin%y,P3%y)
      PMax%y=MAX(PMax%y,P3%y)
    END IF
  END IF


  ! EdgeZ  -ZX - Face  ZY - Face
  IF(Face%V0%Point%z/=Face%V2%Point%z) THEN
    CALL SearchPointsEdgeZaxis(P0,P2,Face%V0,Face%V2, &
                               Face%V0%Point%z,Face%V2%Point%z)
    CALL SearchPointsEdgeZaxis(P1,P3,Face%V1,Face%V3, &
                               Face%V1%Point%z,Face%V3%Point%z)
    IF(P0%z/=P2%z.AND.P0%z/=Face%V0%Point%z) THEN
      PMin%z=MIN(PMin%z,P0%z)
      PMax%z=MAX(PMax%z,P0%z)
    END IF
    IF(P0%z/=P2%z.AND.P2%z/=Face%V2%Point%z) THEN
      PMin%z=MIN(PMin%z,P2%z)
      PMax%z=MAX(PMax%z,P2%z)
    END IF
  
    IF(P1%z/=P3%z.AND.P1%z/=Face%V1%Point%z) THEN
      PMin%z=MIN(PMin%z,P1%z)
      PMax%z=MAX(PMax%z,P1%z)
    END IF
    IF(P1%z/=P3%z.AND.P3%z/=Face%V3%Point%z) THEN
      PMin%z=MIN(PMin%z,P3%z)
      PMax%z=MAX(PMax%z,P3%z)
    END IF
  END IF

END SUBROUTINE SearchMinMaxP4Face


SUBROUTINE AnalyzeFaceRand(Face,NrW,NrMP)
  TYPE (FaceP4_T) :: Face
  INTEGER :: NrW,NrMP

  INTEGER :: i,ecut,iTemp
  INTEGER :: ic,jc,kc,nxc,nyc,nzc
  REAL(8) :: xc0,xc1,yc0,yc1,zc0,zc1,dx,dy,dz
  TYPE (Vertex_T) :: InterVert
  TYPE (Vertex_T) :: Vertex1,Vertex2
  TYPE(Point_T) :: EdgP1,EdgP2,P0,P1,P2,P3
  TYPE(Point_T) :: PMin,PMax,SumMidP,MidPoint
  REAL(8) :: xvc0,xvc1,xvc0u,xvc0o,xvc1u,xvc1o
  REAL(8) :: yvc0,yvc1,yvc0u,yvc0o,yvc1u,yvc1o
  REAL(8) :: zvc0,zvc1
  REAL(8) :: Vol,Len
  REAL(8) :: xsi,ysi 

  !--------            v2 ----- v3
  ! Zählung FaceRand:     |   |
  ! -------            v0 ----- v1

  !Face%Edge[1-4]%in_out Berechnung extern da Reihenfolge nicht gleich (F3 u. F4)
  !Face%in_out=Face%Edge1%in_out+Face%Edge3%in_out
  !Face%in_out=Face%Edge2%in_out+Face%Edge4%in_out
  Face%in_out=Face%V0%in_out+Face%V1%in_out+Face%V2%in_out+Face%V3%in_out
  Face%ec=0
  ! -------
  ! Special
  ! -------
  ! 4 Schnittpunkte
  IF (Face%V0%in_out<=0 .AND. Face%V3%in_out<=0 .AND. &
      Face%V2%in_out>=0 .AND. Face%V1%in_out>=0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  ELSE IF (Face%V0%in_out>=0 .AND. Face%V3%in_out>=0 .AND. &
           Face%V2%in_out<=0 .AND. Face%V1%in_out<=0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  END IF
  ! Face der Cellen mit Cell%vc>0
  ! 1) Face-> Grenzfläche : Vol=0.0, Zählung NrW
  ! 2) Face-> in_out=4 : Vol=MaxVol, Keine Zählung, da maxVol std. gesetzt wird
  ! 3) Face-> in_out=-4, Vol=0.0
  ! 4) Face-> in_out=-2, MAX%[V0-V4]%in_out==0  :   Vol=0.0, Zählung NrW
  ! 5) Face-> Schnitt :  Vol=[>0.0;<maxVol], Zählung NrW
  !................................................
  ! Face-Volumen-Berechnung
  IF (Face%V0%in_out==0.AND.Face%V1%in_out==0.AND. &
      Face%V2%in_out==0.AND.Face%V3%in_out==0) THEN
        Face%Vol=0.0d0
        Face%ec=-1
        !P0=PointParametricEarth(Face%Edge1%Vert1%Point)
        !P1=PointParametricEarth(Face%Edge1%Vert2%Point)
        !P2=PointParametricEarth(Face%Edge3%Vert2%Point)
        !P3=PointParametricEarth(Face%Edge3%Vert1%Point)
        !Face%MidPoint=FaceMidP4(P0,P1,P2,P3)
        !!SumMidP=SumMidP+Vol*MidPoint
        !!Face%MidPoint=SumMidP/Face%Vol
        !NrMP=NrMP+1
        !Face%mp=NrMP
  ELSE IF(Face%in_out==4) THEN
        P0=PointParametricEarth(Face%V0%Point)
        P1=PointParametricEarth(Face%V1%Point)
        P2=PointParametricEarth(Face%V2%Point)
        P3=PointParametricEarth(Face%V3%Point)
        Face%Vol=ABS(FaceP4(P0,P1,P2,P3)) !symbolisch damit Vol!=0
        ! nicht in NrW einfliesen soll
        Face%MidPoint=FaceMidP4(P0,P1,P2,P3)
        NrMP=NrMP+1
        Face%mp=NrMP
  ELSE IF(Face%in_out>=2.AND. &
          MAX(Face%V0%in_out,Face%V1%in_out,Face%V2%in_out,Face%V3%in_out)==1 .AND. &
          MIN(Face%V0%in_out,Face%V1%in_out,Face%V2%in_out,Face%V3%in_out)==0) THEN
          ! 1 --- 1    1 --- 1
          !  |   |      |   |
          ! 0 --- 1    0 --- 0
        P0=PointParametricEarth(Face%V0%Point)
        P1=PointParametricEarth(Face%V1%Point)
        P2=PointParametricEarth(Face%V2%Point)
        P3=PointParametricEarth(Face%V3%Point)
        Face%Vol=ABS(FaceP4(P0,P1,P2,P3)) !symbolisch damit Vol!=0
        ! nicht in NrW einfliesen soll
        Face%MidPoint=FaceMidP4(P0,P1,P2,P3)
        NrMP=NrMP+1
        Face%mp=NrMP
  ELSE IF (Face%in_out==-4) THEN
        Face%Vol=0.0d0
  ELSE IF (Face%in_out==-2 .AND. &
           (MAX(Face%V0%in_out,Face%V1%in_out,Face%V2%in_out,Face%V3%in_out) &
            ==0) ) THEN  ! Schnitt-Grenze Edge
        Face%Vol=0.0d0
        NrW=NrW+1  ! da Face allokiert
  !ELSE IF (Face%ec>0) THEN  ! unter AnalyzeFace
  ELSE IF (MIN(Face%V0%in_out,Face%V1%in_out,Face%V2%in_out,Face%V3%in_out) &
            ==-1 .AND. &
           MAX(Face%V0%in_out,Face%V1%in_out,Face%V2%in_out,Face%V3%in_out) &
            ==1 ) THEN
    Face%ec=1
    PMin=Face%V3%Point  !Face%Edge2%Vert2%Point
    PMax=Face%V0%Point  !Face%Edge1%Vert1%Point
    CALL SearchMinMaxP4Face(Face,PMin,PMax)
    Face%Vol=0.0d0
    nxc=IncrVol
    nyc=IncrVol
    nzc=IncrVol
    xc0=Face%V0%Point%x     !Face%Edge1%Vert1%Point%x
    yc0=Face%V0%Point%y     !Face%Edge1%Vert1%Point%y
    zc0=Face%V0%Point%z     !Face%Edge1%Vert1%Point%z
    xc1=Face%V3%Point%x     !Face%Edge2%Vert2%Point%x
    yc1=Face%V3%Point%y     !Face%Edge2%Vert2%Point%y
    zc1=Face%V3%Point%z     !Face%Edge2%Vert2%Point%z
    !        3-----2   ! generell Point-Folge Face-Vol-Berechnung
    !        | \   |
    !        |   \ |
    !        0-----1
    !---Face_YZ------------------------------------------------------------
    IF (xc0==xc1) THEN
      IF ((Face%Edge4%in_out==2).AND.(Face%Edge2%in_out==-2) .OR. &
          (Face%Edge4%in_out==-2).AND.(Face%Edge2%in_out==2))THEN
         !.........y-direction.............................................
         !     0 ---------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !     |          |   
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc0
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%z-PMin%z
         zvc0=PMin%z+Shrink*Len
         zvc1=PMax%z-Shrink*Len
         dz=(zvc1-zvc0)/nzc
         DO kc=1,nzc
           Vertex1%Point%z=zvc0
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P0=PointParametricEarth(EdgP1)
           P3=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3)  !nur je Teilflaeche
           SumMidP=SumMidP+Vol*MidPoint
           zvc0=MIN(zvc0+dz,zvc1)
         END DO
      ELSE
         !...........z-direction................................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc0
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%y-PMin%y
         yvc0=PMin%y+Shrink*Len
         yvc1=PMax%y-Shrink*Len
         dy=(yvc1-yvc0)/nyc
         DO jc=1,nyc
           Vertex1%Point%y=yvc0
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
           P0=PointParametricEarth(EdgP1)
           P3=PointParametricEarth(EdgP2)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3)  !nur je Teilflaeche
           SumMidP=SumMidP+Vol*MidPoint
           yvc0=MIN(yvc0+dy,yvc1)
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    !---Face_XZ------------------------------------------------------------
    ELSE IF (yc0==yc1) THEN
      IF ((Face%Edge1%in_out==2).AND.(Face%Edge3%in_out==-2) .OR. &
          (Face%Edge1%in_out==-2).AND.(Face%Edge3%in_out==2))THEN
         !........x-direction.......................................... 
         !     0----------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !     |          |
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         zvc0=PMin%z+Shrink*Len
         zvc1=PMax%z-Shrink*Len
         dz=(zvc1-zvc0)/nzc
         DO kc=1,nzc
           Vertex1%Point%z=zvc0
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex2%Point%z=Vertex1%Point%z
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
           P0=PointParametricEarth(EdgP1)
           P3=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           zvc0=zvc0+dz
         END DO
      ELSE
         !..........z-direction..............................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc0
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc1
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         xvc0=PMin%x+Shrink*Len
         xvc1=PMax%x-Shrink*Len
         dx=(xvc1-xvc0)/nxc
         DO ic=1,nxc
           Vertex1%Point%x=xvc0
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
           P0=PointParametricEarth(EdgP1)
           P3=PointParametricEarth(EdgP2)
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           xvc0=xvc0+dx
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    !---Face_XY------------------------------------------------------------
    ELSE IF (zc0==zc1) THEN
      IF ((Face%Edge4%in_out==2).AND.(Face%Edge2%in_out==-2) .OR. &
          (Face%Edge4%in_out==-2).AND.(Face%Edge2%in_out==2))THEN
         !....x-direction..................................................
         !     0----------3   !Point-Folge Cut-Plane Face-Vol-Berechnung
         ! dy  |          |
         !     |          |
         !     1----------2
         Vertex1%Point%x=xc0
         Vertex2%Point%x=xc1
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%y-PMin%y
         yvc0=PMin%y+Shrink*Len
         yvc1=PMax%y-Shrink*Len
         dy=(yvc1-yvc0)/nyc
         DO jc=1,nyc
           Vertex1%Point%y=yvc0
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
           P0=PointParametricEarth(EdgP1)
           P3=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           yvc0=yvc0+dy
         END DO
      ELSE
         !....y-direction...................................................
         !        3-----2   !Point-Folge Cut-Plane Face-Vol-Berechnung
         !        | \   |
         !        |   \ |
         !        0-----1
         Vertex1%Point%y=yc0
         Vertex2%Point%y=yc1
         Vertex1%Point%z=zc0
         Vertex2%Point%z=zc0
         Vertex1%in_out=-1
         Vertex2%in_out=1
         SumMidP%x=0.0d0
         SumMidP%y=0.0d0
         SumMidP%z=0.0d0
         Len=PMax%x-PMin%x
         xvc0=PMin%x+Shrink*Len
         xvc1=PMax%x-Shrink*Len
         dx=(xvc1-xvc0)/nxc
         DO ic=1,nxc
           Vertex1%Point%x=xvc0
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P0=PointParametricEarth(EdgP1)
           P3=PointParametricEarth(EdgP2)
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P1=PointParametricEarth(EdgP1)
           P2=PointParametricEarth(EdgP2)
           Vol=ABS(FaceP4(P0,P1,P2,P3))
           Face%Vol=Face%Vol+Vol
           MidPoint=FaceMidP4(P0,P1,P2,P3) !nur je Teilflaeche
           SumMidP=Vol*MidPoint+SumMidP
           xvc0=xvc0+dx
         END DO
      END IF
      Face%MidPoint=SumMidP/Face%Vol
      NrMP=NrMP+1
      Face%mp=NrMP
    END IF !!(xc0==xc1);(yc0==yc1);(zc0==zc1)
    !NrW=NrW+1
  END IF  !! ec

  !IF ((Face%in_out<4.AND.Face%ec>0).OR.(Face%in_out<=0.AND.Face%ec==-1)) THEN
  ! unter 7.7.1.3 eingearbeitet !! checken
  !IF ((Face%in_out<4.AND.Face%ec>0).OR.(Face%in_out<0.AND.Face%ec==-1).OR.Face%Vol==0.0d0) THEN
  IF ((Face%in_out<4.AND.Face%ec>0)) THEN
    NrW=NrW+1
  END IF

END SUBROUTINE AnalyzeFaceRand


SUBROUTINE EvalBorderFaceXY_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
  INTEGER :: ib,ix,vix,iy,iz,ibn,vjx,fac

  !Übernehme Koordinate der akt RandCelle -> Berechne Face_XY (iz/iz-1)
  TYPE(Point_T) :: P0,P1,P2,P3
  TYPE(FaceP4_T) :: FaceRa
  INTEGER :: NrMP_FXY   ! zur Zt. nicht gezählt
  INTEGER :: NrR_FXY
  NrR_FXY=0;NrMP_FXY=0
  !.........................................
  ! Edge orientation Face_XY: (F1,F2)
  !    P2----------P3
  !     |  > e3 >   |
  !     |^         ^|   !--------                   P2 ----- P3
  !     |e4       e2|   ! Zählung FaceRand Vertex:     |   |
  !     |^          |   ! -------                   P0 ----- P1
  !     |  > e1 >  ^|
  !    P0----------P1
  !.........................................
  !Zuweisung Face-Rand-XY (WE)  für Übergabe (iz)-Ebene
  FaceRa%V0=Floor(ib)%Vertices(vix    ,iy-1,iz)
  FaceRa%V1=Floor(ib)%Vertices(vix+fac,iy-1,iz)
  FaceRa%V2=Floor(ib)%Vertices(vix    ,iy  ,iz)
  FaceRa%V3=Floor(ib)%Vertices(vix+fac,iy  ,iz)
  !.........................................
  !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
  FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
  FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
  FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
  FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
  !.........................................
  CALL AnalyzeFaceRand(FaceRa,NrR_FXY,NrMP_FXY)
  IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
    ALLOCATE(Floor(ib)%Faces_XY(ix,iy,iz)%Face)
    Floor(ib)%Faces_XY(ix,iy,iz)%Face%mp=FaceRa%mp
    Floor(ib)%Faces_XY(ix,iy,iz)%Face%MidPoint=FaceRa%MidPoint
    Floor(ib)%Faces_XY(ix,iy,iz)%Face%Vol=FaceRa%Vol
    Floor(ib)%Faces_XY(ix,iy,iz)%Face%in_out=FaceRa%in_out
    Floor(ib)%Faces_XY(ix,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
  END IF
  !.........................................
  !Zuweisung Face-Rand-XY (WE) für Übergabe (iz-1)-Ebene
  FaceRa%V0=Floor(ib)%Vertices(vix    ,iy-1,iz-1)
  FaceRa%V1=Floor(ib)%Vertices(vix+fac,iy-1,iz-1)
  FaceRa%V2=Floor(ib)%Vertices(vix    ,iy  ,iz-1)
  FaceRa%V3=Floor(ib)%Vertices(vix+fac,iy  ,iz-1)
  !.........................................
  !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
  FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
  FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
  FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
  FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
  !.........................................
  CALL AnalyzeFaceRand(FaceRa,NrR_FXY,NrMP_FXY)
  IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
    ALLOCATE(Floor(ib)%Faces_XY(ix,iy,iz-1)%Face)
    Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%mp=FaceRa%mp
    Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%MidPoint=FaceRa%MidPoint
    Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%Vol=FaceRa%Vol
    Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%in_out=FaceRa%in_out
    Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
  END IF
  Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+NrR_FXY
  Floor(ib)%NrW_FacesXY=Floor(ib)%NrW_FacesXY+NrR_FXY
END SUBROUTINE EvalBorderFaceXY_WE

SUBROUTINE EvalBorderFaceYZ_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
  INTEGER :: ib,ix,vix,iy,iz,ibn,vjx,fac

  !Übernehme Koordinate der akt RandCelle -> Berechne Face_YZ (vix+fac) 
  TYPE(Point_T) :: P0,P1,P2,P3
  TYPE(FaceP4_T) :: FaceRa
  INTEGER :: NrMP_FYZ  ! zur Zt. nicht gezählt
  INTEGER :: NrR_FYZ
  NrR_FYZ=0;NrMP_FYZ=0
  !.........................................
  ! Edge orientation Face_YZ: (F5,F6)
  !    P2----------P3
  !     |  > e3 >   |
  !     |^         ^|   !--------                   P2 ----- P3
  !     |e4       e2|   ! Zählung FaceRand Vertex:     |   |
  !     |^          |   ! -------                   P0 ----- P1
  !     |  > e1 >  ^|
  !    P0----------P1
  !.........................................
  !Zuweisung Face-Rand-YZ (iWE) für Übergabe (vix+fac) Ebene
  FaceRa%V0=Floor(ibn)%Vertices(jx0+fac,iy-1,iz-1)
  FaceRa%V1=Floor(ibn)%Vertices(jx0+fac,iy  ,iz-1)
  FaceRa%V2=Floor(ibn)%Vertices(jx0+fac,iy-1,iz)
  FaceRa%V3=Floor(ibn)%Vertices(jx0+fac,iy  ,iz)
  !.........................................
  !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
  FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
  FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
  FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
  FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
  !.........................................
  CALL AnalyzeFaceRand(FaceRa,NrR_FYZ,NrMP_FYZ)
  IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
    ALLOCATE(Floor(ib)%Faces_YZ(vix+fac,iy,iz)%Face)
    Floor(ib)%Faces_YZ(vix+fac,iy,iz)%Face%mp=FaceRa%mp
    Floor(ib)%Faces_YZ(vix+fac,iy,iz)%Face%MidPoint=FaceRa%MidPoint
    Floor(ib)%Faces_YZ(vix+fac,iy,iz)%Face%MidPoint%x= &
          & Floor(ib)%Vertices(vix+fac,iy,iz)%Point%x
    Floor(ib)%Faces_YZ(vix+fac,iy,iz)%Face%Vol=FaceRa%Vol
    Floor(ib)%Faces_YZ(vix+fac,iy,iz)%Face%in_out=FaceRa%in_out
    Floor(ib)%Faces_YZ(vix+fac,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
  END IF
  Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+NrR_FYZ
  Floor(ib)%NrW_FacesYZ=Floor(ib)%NrW_FacesYZ+NrR_FYZ
END SUBROUTINE EvalBorderFaceYZ_WE


SUBROUTINE EvalBorderFaceZX_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
  INTEGER :: ib,ix,vix,iy,iz,ibn,vjx,fac

  !Übernehme Koordinate der akt RandCelle -> Berechne Face_ZX iy/iy-1
  TYPE(Point_T) :: P0,P1,P2,P3
  TYPE(FaceP4_T) :: FaceRa
  INTEGER :: NrMP_FZX  ! zur Zeit nicht gezählt 
  INTEGER :: NrR_FZX
  NrR_FZX=0;NrMP_FZX=0
  !.........................................
  ! Edge orientation Face_ZX: (F3,F4)
  !    P2----------P3
  !     |  > e2 >   |
  !     |^         ^|   !--------                   P2 ----- P3
  !     |e1       e3|   ! Zählung FaceRand Vertex:     |   |
  !     |^          |   ! -------                   P0 ----- P1
  !     |  > e4 >  ^|
  !    P0----------P1
  !.........................................
  !Zuweisung Face-Rand-ZX (WE) für Übergabe (iy-1)Ebene
  FaceRa%V0=Floor(ib)%Vertices(vix    ,iy-1,iz-1)
  FaceRa%V1=Floor(ib)%Vertices(vix+fac,iy-1,iz-1)
  FaceRa%V2=Floor(ib)%Vertices(vix    ,iy-1,iz)
  FaceRa%V3=Floor(ib)%Vertices(vix+fac,iy-1,iz)
  !.........................................
  !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
  FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
  FaceRa%Edge2%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
  FaceRa%Edge3%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
  FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
  !.........................................
  CALL AnalyzeFaceRand(FaceRa,NrR_FZX,NrMP_FZX)
  IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
    ALLOCATE(Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face)
    Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%mp=FaceRa%mp
    Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%MidPoint=FaceRa%MidPoint
    Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%Vol=FaceRa%Vol
    Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%in_out=FaceRa%in_out
    Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
  END IF
  !.........................................
  !Zuweisung Face-Rand-ZX (WE) für Übergabe (iy)Ebene
  FaceRa%V0=Floor(ib)%Vertices(vix    ,iy,iz-1)
  FaceRa%V1=Floor(ib)%Vertices(vix+fac,iy,iz-1)
  FaceRa%V2=Floor(ib)%Vertices(vix    ,iy,iz)
  FaceRa%V3=Floor(ib)%Vertices(vix+fac,iy,iz)
  !.........................................
  !Orientierung Edge1-Edge4 in F3,F4(ZX) verschieden als F1,F2,F5,F6
  FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
  FaceRa%Edge2%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
  FaceRa%Edge3%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
  FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
  !.........................................
  CALL AnalyzeFaceRand(FaceRa,NrR_FZX,NrMP_FZX)
  IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
    ALLOCATE(Floor(ib)%Faces_ZX(ix,iy,iz)%Face)
    Floor(ib)%Faces_ZX(ix,iy,iz)%Face%mp=FaceRa%mp
    Floor(ib)%Faces_ZX(ix,iy,iz)%Face%MidPoint=FaceRa%MidPoint
    Floor(ib)%Faces_ZX(ix,iy,iz)%Face%Vol=FaceRa%Vol
    Floor(ib)%Faces_ZX(ix,iy,iz)%Face%in_out=FaceRa%in_out
    Floor(ib)%Faces_ZX(ix,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
  END IF
  Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+NrR_FZX
  Floor(ib)%NrW_FacesZX=Floor(ib)%NrW_FacesZX+NrR_FZX
END SUBROUTINE EvalBorderFaceZX_WE



SUBROUTINE EvalBorderFaceXY_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
  INTEGER :: ib,ix,iy,viy,iz,ibn,vjy,fac

  !Übernehme Koordinate der akt RandCelle -> Berechne Face_XY iz/iz-1
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(FaceP4_T) :: FaceRa
  INTEGER :: NrMP_FXY  ! zur Zt. nicht gezählt
  INTEGER :: NrR_FXY
  NrMP_FXY=0;NrR_FXY=0
     !.........................................
     ! Edge orientation Face_XY: (F1,F2)
     !    P2----------P3
     !     |  > e3 >   |
     !     |^         ^|   !--------                   P2 ----- P3
     !     |e4       e2|   ! Zählung FaceRand Vertex:     |   |
     !     |^          |   ! -------                   P0 ----- P1
     !     |  > e1 >  ^|
     !    P0----------P1
     !.........................................
     !Zuweisung Face-Rand-XY (NS) für Übergabe (iz-1)Ebene
     FaceRa%V0=Floor(ib)%Vertices(ix-1,viy    ,iz-1)
     FaceRa%V1=Floor(ib)%Vertices(ix  ,viy    ,iz-1)
     FaceRa%V2=Floor(ib)%Vertices(ix-1,viy+fac,iz-1)
     FaceRa%V3=Floor(ib)%Vertices(ix  ,viy+fac,iz-1)
     !.........................................
     !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
     FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
     FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
     FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
     FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
     !.........................................
     CALL AnalyzeFaceRand(FaceRa,NrR_FXY,NrMP_FXY)
     IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
       ALLOCATE(Floor(ib)%Faces_XY(ix,iy,iz-1)%Face)
       Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%mp=FaceRa%mp
       Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%MidPoint=FaceRa%MidPoint
       Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%Vol=FaceRa%Vol
       Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%in_out=FaceRa%in_out
       Floor(ib)%Faces_XY(ix,iy,iz-1)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
     END IF
     !............................................
     !Zuweisung Face-Rand-XY (NS) für Übergabe (iz) Ebene
     FaceRa%V0=Floor(ib)%Vertices(ix-1,viy    ,iz)
     FaceRa%V1=Floor(ib)%Vertices(ix  ,viy    ,iz)
     FaceRa%V2=Floor(ib)%Vertices(ix-1,viy+fac,iz)
     FaceRa%V3=Floor(ib)%Vertices(ix  ,viy+fac,iz)
     !.........................................
     !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
     FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
     FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
     FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
     FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
     !.........................................
     CALL AnalyzeFaceRand(FaceRa,NrR_FXY,NrMP_FXY)
     IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
       ALLOCATE(Floor(ib)%Faces_XY(ix,iy,iz)%Face)
       Floor(ib)%Faces_XY(ix,iy,iz)%Face%mp=FaceRa%mp
       Floor(ib)%Faces_XY(ix,iy,iz)%Face%MidPoint=FaceRa%MidPoint
       Floor(ib)%Faces_XY(ix,iy,iz)%Face%Vol=FaceRa%Vol
       Floor(ib)%Faces_XY(ix,iy,iz)%Face%in_out=FaceRa%in_out
       Floor(ib)%Faces_XY(ix,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
     END IF
     Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+NrR_FXY
     Floor(ib)%NrW_FacesXY=Floor(ib)%NrW_FacesXY+NrR_FXY
     !............................................
END SUBROUTINE EvalBorderFaceXY_NS

SUBROUTINE EvalBorderFaceYZ_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
  INTEGER :: ib,ix,iy,viy,iz,ibn,vjy,fac
  !Übernehme Koordinate der akt RandCelle -> Berechne Face_YZ ix/ix-1
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(FaceP4_T) :: FaceRa
  INTEGER :: NrMP_FYZ  ! zur Zt. nicht gezählt
  INTEGER :: NrR_FYZ
  NrMP_FYZ=0;NrR_FYZ=0
     !.........................................
     ! Edge orientation Face_YZ: (F5,F6)
     !    P2----------P3
     !     |  > e3 >   |
     !     |^         ^|   !--------                   P2 ----- P3
     !     |e4       e2|   ! Zählung FaceRand Vertex:     |   |
     !     |^          |   ! -------                   P0 ----- P1
     !     |  > e1 >  ^|
     !    P0----------P1
     !.........................................
     !Zuweisung Face-Rand-YZ (NS) für Übergabe (ix-1) Ebene
     FaceRa%V0=Floor(ib)%Vertices(ix-1,viy    ,iz-1)
     FaceRa%V1=Floor(ib)%Vertices(ix-1,viy+fac,iz-1)
     FaceRa%V2=Floor(ib)%Vertices(ix-1,viy    ,iz)
     FaceRa%V3=Floor(ib)%Vertices(ix-1,viy+fac,iz)
     !.........................................
     !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
     FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
     FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
     FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
     FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
     !.........................................
     CALL AnalyzeFaceRand(FaceRa,NrR_FYZ,NrMP_FYZ)
     IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
       ALLOCATE(Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face)
       Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%mp=FaceRa%mp
       Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%MidPoint=FaceRa%MidPoint
       Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%Vol=FaceRa%Vol
       Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%in_out=FaceRa%in_out
       Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
     END IF
     !............................................
     !Zuweisung Face-Rand-YZ (NS) für Übergabe (ix) Ebene
     FaceRa%V0=Floor(ib)%Vertices(ix,viy    ,iz-1)
     FaceRa%V1=Floor(ib)%Vertices(ix,viy+fac,iz-1)
     FaceRa%V2=Floor(ib)%Vertices(ix,viy    ,iz)
     FaceRa%V3=Floor(ib)%Vertices(ix,viy+fac,iz)
     !.........................................
     !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
     FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
     FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
     FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
     FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
     !.........................................
     CALL AnalyzeFaceRand(FaceRa,NrR_FYZ,NrMP_FYZ)
     IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
       ALLOCATE(Floor(ib)%Faces_YZ(ix,iy,iz)%Face)
       Floor(ib)%Faces_YZ(ix,iy,iz)%Face%mp=FaceRa%mp
       Floor(ib)%Faces_YZ(ix,iy,iz)%Face%MidPoint=FaceRa%MidPoint
       Floor(ib)%Faces_YZ(ix,iy,iz)%Face%Vol=FaceRa%Vol
       Floor(ib)%Faces_YZ(ix,iy,iz)%Face%in_out=FaceRa%in_out
       Floor(ib)%Faces_YZ(ix,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
     END IF
     Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+NrR_FYZ
     Floor(ib)%NrW_FacesYZ=Floor(ib)%NrW_FacesYZ+NrR_FYZ
     !............................................
END SUBROUTINE EvalBorderFaceYZ_NS

SUBROUTINE EvalBorderFaceZX_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
  INTEGER :: ib,ix,iy,viy,iz,ibn,vjy,fac
  !Übernehme Koordinate der akt RandCelle -> Berechne Face_ZX (viy+fac)
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(FaceP4_T) :: FaceRa
  INTEGER :: fiy
  INTEGER :: NrMP_FZX  ! zur Zt. nicht gezählt
  INTEGER :: NrR_FZX
  NrMP_FZX=0;NrR_FZX=0
     !.........................................
     ! Edge orientation Face_ZX: (F3.OR.F4)
     !    P2----------P3
     !     |  > e2 >   |
     !     |^         ^|   !--------                   V2 ----- V3
     !     |e1       e3|   ! Zählung FaceRand Vertex:     |   |
     !     |^          |   ! -------                   V0 ----- V1
     !     |  > e4 >  ^|
     !    P0----------P1
     !.........................................
     !Zuweisung Face-Rand-ZX (NS) für Übergabe (viy+fac)Ebene
     fiy=viy+fac
     FaceRa%V0=Floor(ib)%Vertices(ix-1,fiy,iz-1)
     FaceRa%V1=Floor(ib)%Vertices(ix  ,fiy,iz-1)
     FaceRa%V2=Floor(ib)%Vertices(ix-1,fiy,iz)
     FaceRa%V3=Floor(ib)%Vertices(ix  ,fiy,iz)
     !.........................................
     !Orientierung Edge1-Edge4 in F3,F4(ZX) verschieden als F1,F2,F5,F6
     FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
     FaceRa%Edge2%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
     FaceRa%Edge3%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
     FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
     !.........................................
     CALL AnalyzeFaceRand(FaceRa,NrR_FZX,NrMP_FZX)
     IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
       ALLOCATE(Floor(ib)%Faces_ZX(ix,fiy,iz)%Face)
       Floor(ib)%Faces_ZX(ix,fiy,iz)%Face%mp=FaceRa%mp
       Floor(ib)%Faces_ZX(ix,fiy,iz)%Face%MidPoint=FaceRa%MidPoint
       Floor(ib)%Faces_ZX(ix,fiy,iz)%Face%Vol=FaceRa%Vol
       Floor(ib)%Faces_ZX(ix,fiy,iz)%Face%in_out=FaceRa%in_out
       Floor(ib)%Faces_ZX(ix,fiy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
     END IF
     Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+NrR_FZX
     Floor(ib)%NrW_FacesZX=Floor(ib)%NrW_FacesZX+NrR_FZX
     !............................................
END SUBROUTINE EvalBorderFaceZX_NS

SUBROUTINE EvalBorderFaceXY_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
  INTEGER :: ib,ix,iy,iz,viz,ibn,vjz,fac
  !Übernehme Koordinate der akt RandCelle -> Berechne Face_XY (iz).OR.(iz-1) 
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(FaceP4_T) :: FaceRa
  INTEGER  :: fiz
  INTEGER :: NrMP_FXY  ! zur Zt. nicht gezählt
  INTEGER :: NrR_FXY
  NrMP_FXY=0;NrR_FXY=0
    !.........................................
    ! Edge orientation Face_XY: (F1,F2)
    !    P2----------P3
    !     |  > e3 >   |
    !     |^         ^|   !--------                   P2 ----- P3
    !     |e4       e2|   ! Zählung FaceRand Vertex:     |   |
    !     |^          |   ! -------                   P0 ----- P1
    !     |  > e1 >  ^|
    !    P0----------P1
    !.........................................
    !Zuweisung Face-Rand-XY (TB) für Übergabe (viz+fac)Ebene
    fiz=viz+fac
    FaceRa%V0=Floor(ib)%Vertices(ix-1,iy-1,fiz)
    FaceRa%V1=Floor(ib)%Vertices(ix  ,iy-1,fiz)
    FaceRa%V2=Floor(ib)%Vertices(ix-1,iy  ,fiz)
    FaceRa%V3=Floor(ib)%Vertices(ix  ,iy  ,fiz)
    !.........................................
    !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
    FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
    FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
    FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
    FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
    !.........................................
    CALL AnalyzeFaceRand(FaceRa,NrR_FXY,NrMP_FXY)
    IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
      ALLOCATE(Floor(ib)%Faces_XY(ix,iy,fiz)%Face)
      Floor(ib)%Faces_XY(ix,iy,fiz)%Face%mp=FaceRa%mp
      Floor(ib)%Faces_XY(ix,iy,fiz)%Face%MidPoint=FaceRa%MidPoint
      Floor(ib)%Faces_XY(ix,iy,fiz)%Face%Vol=FaceRa%Vol
      Floor(ib)%Faces_XY(ix,iy,fiz)%Face%in_out=FaceRa%in_out
      Floor(ib)%Faces_XY(ix,iy,fiz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
    END IF
    Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+NrR_FXY
    Floor(ib)%NrW_FacesXY=Floor(ib)%NrW_FacesXY+NrR_FXY
    !.........................................
END SUBROUTINE EvalBorderFaceXY_TB


SUBROUTINE EvalBorderFaceYZ_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
  INTEGER :: ib,ix,iy,iz,viz,ibn,vjz,fac
  !Übernehme Koordinate der akt RandCelle -> Berechne Face_YZ (ix/ix-1)
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(FaceP4_T) :: FaceRa
  INTEGER :: NrMP_FYZ     ! zur Zt. nicht gezählt 
  INTEGER :: NrR_FYZ
  NrMP_FYZ=0;NrR_FYZ=0
    !.........................................
    ! Edge orientation Face_YZ: (F5,F6)
    !    P2----------P3
    !     |  > e3 >   |
    !     |^         ^|   !--------                   P2 ----- P3
    !     |e4       e2|   ! Zählung FaceRand Vertex:     |   |
    !     |^          |   ! -------                   P0 ----- P1
    !     |  > e1 >  ^|
    !    P0----------P1
    !.........................................
    !Zuweisung Face-Rand-YZ (TB) für Übergabe (ix-1)Ebene
    FaceRa%V0=Floor(ib)%Vertices(ix-1,iy-1,viz)
    FaceRa%V1=Floor(ib)%Vertices(ix-1,iy  ,viz)
    FaceRa%V2=Floor(ib)%Vertices(ix-1,iy-1,viz+fac)
    FaceRa%V3=Floor(ib)%Vertices(ix-1,iy  ,viz+fac)
    !.........................................
    !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
    FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
    FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
    FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
    FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
    !.........................................
    CALL AnalyzeFaceRand(FaceRa,NrR_FYZ,NrMP_FYZ)
    IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
      ALLOCATE(Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face)
      Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%mp=FaceRa%mp
      Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%MidPoint=FaceRa%MidPoint
      Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%Vol=FaceRa%Vol
      Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%in_out=FaceRa%in_out
      Floor(ib)%Faces_YZ(ix-1,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
    END IF
    !.........................................
    !Zuweisung Face-Rand-YZ (TB) für Übergabe (ix)Ebene
    FaceRa%V0=Floor(ib)%Vertices(ix,iy-1,viz)
    FaceRa%V1=Floor(ib)%Vertices(ix,iy  ,viz)
    FaceRa%V2=Floor(ib)%Vertices(ix,iy-1,viz+fac)
    FaceRa%V3=Floor(ib)%Vertices(ix,iy  ,viz+fac)
    !.........................................
    !Orientierung Edge1-Edge4 in F3 u. F4 verschieden als F1,F2,F5,F6
    FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
    FaceRa%Edge2%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
    FaceRa%Edge3%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
    FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
    !.........................................
    CALL AnalyzeFaceRand(FaceRa,NrR_FYZ,NrMP_FYZ)
    IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
      ALLOCATE(Floor(ib)%Faces_YZ(ix,iy,iz)%Face)
      Floor(ib)%Faces_YZ(ix,iy,iz)%Face%mp=FaceRa%mp
      Floor(ib)%Faces_YZ(ix,iy,iz)%Face%MidPoint=FaceRa%MidPoint
      Floor(ib)%Faces_YZ(ix,iy,iz)%Face%Vol=FaceRa%Vol
      Floor(ib)%Faces_YZ(ix,iy,iz)%Face%in_out=FaceRa%in_out
      Floor(ib)%Faces_YZ(ix,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
    END IF
    Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+NrR_FYZ
    Floor(ib)%NrW_FacesYZ=Floor(ib)%NrW_FacesYZ+NrR_FYZ
    !...................................................
END SUBROUTINE EvalBorderFaceYZ_TB

SUBROUTINE EvalBorderFaceZX_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
  INTEGER :: ib,ix,iy,iz,viz,ibn,vjz,fac
  !Übernehme Koordinate der akt RandCelle -> Berechne Face_ZX (iy/iy-1)
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(FaceP4_T) :: FaceRa
  INTEGER :: NrMP_FZX   ! zur Zt. nicht gezählt
  INTEGER :: NrR_FZX
  NrMP_FZX=0;NrR_FZX=0
    !.........................................
    ! Edge orientation Face_ZX: (F3,F4)
    !    P2----------P3
    !     |  > e2 >   |
    !     |^         ^|   !--------                   V2 ----- V3
    !     |e1       e3|   ! Zählung FaceRand Vertex:     |   |
    !     |^          |   ! -------                   V0 ----- V1
    !     |  > e4 >  ^|
    !    P0----------P1
    !.........................................
    !Zuweisung Face-Rand-ZX (TB) für Übergabe (iy-1)Ebene
    FaceRa%V0=Floor(ib)%Vertices(ix-1,iy-1,viz)
    FaceRa%V1=Floor(ib)%Vertices(ix  ,iy-1,viz)
    FaceRa%V2=Floor(ib)%Vertices(ix-1,iy-1,viz+fac)
    FaceRa%V3=Floor(ib)%Vertices(ix  ,iy-1,viz+fac)
    !.........................................
    !Orientierung Edge1-Edge4 in F3,F4(ZX) verschieden als F1,F2,F5,F6
    FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
    FaceRa%Edge2%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
    FaceRa%Edge3%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
    FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
    !.........................................
    CALL AnalyzeFaceRand(FaceRa,NrR_FZX,NrMP_FZX)
    IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
      ALLOCATE(Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face)
      Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%mp=FaceRa%mp
      Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%MidPoint=FaceRa%MidPoint
      Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%Vol=FaceRa%Vol
      Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%in_out=FaceRa%in_out
      Floor(ib)%Faces_ZX(ix,iy-1,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
    END IF
    !.........................................
    !Zuweisung Face-Rand-ZX (TB) für Übergabe (iy)Ebene
    FaceRa%V0=Floor(ib)%Vertices(ix-1,iy,viz)
    FaceRa%V1=Floor(ib)%Vertices(ix  ,iy,viz)
    FaceRa%V2=Floor(ib)%Vertices(ix-1,iy,viz+fac)
    FaceRa%V3=Floor(ib)%Vertices(ix  ,iy,viz+fac)
    !.........................................
    !Orientierung Edge1-Edge4 in F3,F4(ZX) verschieden als F1,F2,F5,F6
    FaceRa%Edge1%in_out=FaceRa%V0%in_out+FaceRa%V2%in_out
    FaceRa%Edge2%in_out=FaceRa%V2%in_out+FaceRa%V3%in_out
    FaceRa%Edge3%in_out=FaceRa%V1%in_out+FaceRa%V3%in_out
    FaceRa%Edge4%in_out=FaceRa%V0%in_out+FaceRa%V1%in_out
    !.........................................
    CALL AnalyzeFaceRand(FaceRa,NrR_FZX,NrMP_FZX)
    IF(FaceRa%in_out>-4.AND.FaceRa%ec>0) THEN
      ALLOCATE(Floor(ib)%Faces_ZX(ix,iy,iz)%Face)
      Floor(ib)%Faces_ZX(ix,iy,iz)%Face%mp=FaceRa%mp
      Floor(ib)%Faces_ZX(ix,iy,iz)%Face%MidPoint=FaceRa%MidPoint
      Floor(ib)%Faces_ZX(ix,iy,iz)%Face%Vol=FaceRa%Vol
      Floor(ib)%Faces_ZX(ix,iy,iz)%Face%in_out=FaceRa%in_out
      Floor(ib)%Faces_ZX(ix,iy,iz)%Face%ec=FaceRa%ec !theoretisch/symbolisch für Weight-Ausgabe
    END IF
    Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+NrR_FZX
    Floor(ib)%NrW_FacesZX=Floor(ib)%NrW_FacesZX+NrR_FZX
    !...................................................
END SUBROUTINE EvalBorderFaceZX_TB

SUBROUTINE EvalBorderAllFace_OutW(ib,ix0,iy0,iy1,iz0,iz1)
  INTEGER :: ib,ix0,iy0,iy1,iz0,iz1
  
  INTEGER :: ix,iy,iz
  INTEGER :: NrR_FXY,NrR_FYZ,NrR_FZX,NrMP_FXY,NrMP_FYZ,NrMP_FZX
  !!1 !......................................
  !!1 !Für BorderWest -Face FU (YZ) aktiv
  !!1 !......................................
  !CALL Set(Floor(ib))
  CALL Set_EnvBlkFace(Floor(ib),ib)
  NrR_FXY=0;NrR_FYZ=0;NrR_FZX=0
  NrMP_FXY=0;NrMP_FYZ=0;NrMP_FZX=0  ! noch im Program umsetzen
  !!1 !----Faces_XY---------
  !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2
  !!1 DO iy=iy0+1,iy1
  !!1    iz=iz0
  !!1    DO WHILE(.NOT.ASSOCIATED(Faces_XY(ix0+1,iy,iz)%Face).AND.iz/=iz1)
  !!1        !           "  RandFace_XY unter Berg west"
  !!1        !Write(*,*) "ib=",ib, " Type-ow:  ix0=",ix0," iy=",iy," iz=",i, &
  !!1        iz=iz+1
  !!1    END DO 
  !!1    DO WHILE(ASSOCIATED(Faces_XY(ix0+1,iy,iz)%Face).AND.iz/=iz1+1)
  !!1        ALLOCATE(Floor(ib)%Faces_XY(ix0,iy,iz)%Face)
  !!1        Floor(ib)%Faces_XY(ix0,iy,iz)%Face%Vol=Floor(ib)%Faces_XY(ix0+1,iy,iz)%Face%Vol
  !!1        Floor(ib)%Faces_XY(ix0,iy,iz)%Face%MidPoint= &
  !!1              Floor(ib)%Faces_XY(ix0+1,iy,iz)%Face%MidPoint
  !!1        Floor(ib)%Faces_XY(ix0,iy,iz)%Face%MidPoint%x= &
  !!1              Floor(ib)%Faces_XY(ix0+1,iy,iz)%Face%MidPoint%x-dx(ix0+1)
  !!1        Floor(ib)%Faces_XY(ix0,iy,iz)%Face%in_out=Floor(ib)%Faces_XY(ix0+1,iy,iz)%Face%in_out
  !!1        Floor(ib)%Faces_XY(ix0,iy,iz)%Face%ec=Floor(ib)%Faces_XY(ix0+1,iy,iz)%Face%ec
  !!1        !IF(Floor(ib)%Faces_XY(ix0,iy,iz)%Face%Vol>0.0d0) THEN
  !!1        !Aenderung: Betreff Allocate EdgeX,EdgeY Face_XY,wenn iz==1 ?Laufvariable in InitAllEdge
  !!1        !?Fall in DurranHatS800.grid, wenn ix=1...,iy=1 iz=1 Vertex%in_out==1  sind 
  !!1        IF((Floor(ib)%Faces_XY(ix0,iy,iz)%Face%Vol>0.0d0 .AND. &
  !!1            Floor(ib)%Faces_XY(ix0,iy,iz)%Face%in_out<4)       .OR. &
  !!1           (Floor(ib)%Faces_XY(ix0,iy,iz)%Face%Vol==0.0d0 .AND. &
  !!1            Floor(ib)%Faces_XY(ix0,iy,iz)%Face%in_out==0) ) THEN
  !!1            !Floor(ib)%Faces_XY(ix0,iy,iz)%Face%numberVert==4) ) THEN
  !!1            !Beachte Grenzfall: Face_XY Grenzfläche
  !!1          NrR_FXY=NrR_FXY+1
  !!1          IF(Floor(ib)%Faces_XY(ix0,iy,iz)%Face%ec>0) THEN
  !!1            ! Faces_XY mit maxVol nicht unter NrMP_FacesXY
  !!1            IF(Floor(ib)%Faces_XY(ix0+1,iy,iz)%Face%mp>0) THEN
  !!1              !Zusatzabfrage notwendig,da allokierte Face unterhalb Berg sein können
  !!1              NrMP_FacesXY=NrMP_FacesXY+1
  !!1              Floor(ib)%Faces_XY(ix0,iy,iz)%Face%mp=NrMP_FacesXY
  !!1            ELSE
  !!1              Floor(ib)%Faces_XY(ix0,iy,iz)%Face%mp=0
  !!1            END IF
  !!1          END IF
  !!1        END IF
  !!1        !Write(*,*) "ib=",ib, " Type-ow:  ix0=",ix0," iy=",iy," iz=",iz, " Schnitt RandFace_XY west"
  !!1        iz=iz+1
  !!1    END DO 
  !!1 END DO
  !----Faces_YZ-------
  DO iy=iy0+1,iy1
    DO iz=iz0+1,iz1
       !DO WHILE(.NOT.ASSOCIATED(Faces_YZ(ix0,iy,iz)%Face).AND.iz/=iz1)
         !           "  RandFace_YZ unter Berg west"
         !Write(*,*) "ib=",ib, " Type-ow:  ix0=",ix0," iy=",iy," iz=",i, &
        ! iz=iz+1
       !END DO 
       !DO WHILE(ASSOCIATED(Faces_YZ(ix0,iy,iz)%Face).AND.iz/=iz1+1)
       IF(ASSOCIATED(Faces_YZ(ix0,iy,iz)%Face).AND.iz/=iz1+1) THEN
       !IF(Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%ec>0) THEN
       IF( (Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%in_out<4.AND.  &
            Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%ec>0) &
           .OR. &
           (Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%in_out<0.AND.  &
            Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%ec==-1)   &
           .OR. &
           Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%Vol==0.0d0 ) THEN
         ALLOCATE(Floor(ib)%Faces_YZ(ix0-1,iy,iz)%Face)
         Floor(ib)%Faces_YZ(ix0-1,iy,iz)%Face%Vol=Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%Vol
         Floor(ib)%Faces_YZ(ix0-1,iy,iz)%Face%MidPoint= &
               Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%MidPoint
         Floor(ib)%Faces_YZ(ix0-1,iy,iz)%Face%MidPoint%x= &
               Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%MidPoint%x-dx(ix0+1)
         Floor(ib)%Faces_YZ(ix0-1,iy,iz)%Face%in_out=Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%in_out
         Floor(ib)%Faces_YZ(ix0-1,iy,iz)%Face%ec=Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%ec
         NrR_FYZ=NrR_FYZ+1
         IF(Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%Vol>0.0d0) THEN
           ! Faces_YZ mit maxVol nicht unter NrMP_FacesYZ
           IF(Floor(ib)%Faces_YZ(ix0,iy,iz)%Face%mp>0) THEN
              !Zusatzabfrage notwendig,da allokierte Face unterhalb Berg sein können
              NrMP_FacesYZ=NrMP_FacesYZ+1
              Floor(ib)%Faces_YZ(ix0-1,iy,iz)%Face%mp=NrMP_FacesYZ
            ELSE
              Floor(ib)%Faces_YZ(ix0-1,iy,iz)%Face%mp=0
           END IF
         END IF
       END IF
       !Write(*,*) "ib=",ib, " Type-ow:  ix0=",ix0," iy=",iy," iz=",iz, " Schnitt RandFace_YZ west"
       !iz=iz+1
      END IF
     END DO  !iz
  END DO  !iy
  !!1 !----Faces_ZX-------
  !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2
  !!1 DO iy=iy0,iy1
  !!1    iz=iz0+1
  !!1    DO WHILE(.NOT.ASSOCIATED(Faces_ZX(ix0+1,iy,iz)%Face).AND.iz/=iz1)
  !!1        !           "  RandFace_ZX unter Berg west"
  !!1        !Write(*,*) "ib=",ib, " Type-ow:  ix0=",ix0," iy=",iy," iz=",i, &
  !!1        iz=iz+1
  !!1    END DO 
  !!1    DO WHILE(ASSOCIATED(Faces_ZX(ix0+1,iy,iz)%Face).AND.iz/=iz1+1)
  !!1        ALLOCATE(Floor(ib)%Faces_ZX(ix0,iy,iz)%Face)
  !!1        Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%Vol=Floor(ib)%Faces_ZX(ix0+1,iy,iz)%Face%Vol
  !!1        Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%MidPoint= &
  !!1              Floor(ib)%Faces_ZX(ix0+1,iy,iz)%Face%MidPoint
  !!1        Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%MidPoint%x= &
  !!1              Floor(ib)%Faces_ZX(ix0+1,iy,iz)%Face%MidPoint%x-dx(ix0+1)
  !!1        Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%in_out=Floor(ib)%Faces_ZX(ix0+1,iy,iz)%Face%in_out
  !!1        Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%ec=Floor(ib)%Faces_ZX(ix0+1,iy,iz)%Face%ec
  !!1        IF(Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%Vol>0.0d0) THEN
  !!1          NrR_FZX=NrR_FZX+1
  !!1          IF(Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%ec>0) THEN
  !!1            ! Faces_ZX mit maxVol nicht unter NrMP_FacesZX
  !!1            IF(Floor(ib)%Faces_ZX(ix0+1,iy,iz)%Face%mp>0) THEN
  !!1              !Zusatzabfrage notwendig,da allokierte Face unterhalb Berg sein können
  !!1              NrMP_FacesZX=NrMP_FacesZX+1
  !!1              Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%mp=NrMP_FacesZX
  !!1            ELSE
  !!1              Floor(ib)%Faces_ZX(ix0,iy,iz)%Face%mp=0
  !!1            END IF
  !!1          END IF
  !!1        END IF
  !!1        !Write(*,*) "ib=",ib, " Type-ow:  ix0=",ix0," iy=",iy," iz=",iz, " Schnitt RandFace_ZX west"
  !!1        iz=iz+1
  !!1    END DO 
  !!1 END DO
  !!1 Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+NrR_FXY
  Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+NrR_FYZ
  !!1 Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+NrR_FZX
  !!1 !Floor(ib)%NrMP_FacesXY=Floor(ib)%NrMP_FacesXY+NrMP_FXY
  !!1 !Floor(ib)%NrMP_FacesYZ=Floor(ib)%NrMP_FacesYZ+NrMP_FYZ
  !!1 !Floor(ib)%NrMP_FacesZX=Floor(ib)%NrMP_FacesZX+NrMP_FZX
  !!1 Floor(ib)%NrW_FacesXY=NrW_FacesXY+NrR_FXY
  Floor(ib)%NrW_FacesYZ=NrW_FacesYZ+NrR_FYZ
  !!1 Floor(ib)%NrW_FacesZX=NrW_FacesZX+NrR_FZX

END SUBROUTINE EvalBorderAllFace_OutW

SUBROUTINE EvalBorderAllFace_OutE(ib,ix1,iy0,iy1,iz0,iz1)
   INTEGER :: ib,ix1,iy0,iy1,iz0,iz1

   INTEGER :: ix,iy,iz
   INTEGER :: NrR_FXY,NrR_FYZ,NrR_FZX,NrMP_FXY,NrMP_FYZ,NrMP_FZX
   !!1 !......................................
   !!1 !Für BorderEast -Face FU (YZ) aktiv
   !!1 !......................................
   !CALL Set(Floor(ib))
   CALL Set_EnvBlkFace(Floor(ib),ib)
   NrR_FXY=0;NrR_FYZ=0;NrR_FZX=0
   NrMP_FXY=0;NrMP_FYZ=0;NrMP_FZX=0  ! noch im Program umsetzen
   !!1 !----Faces_XY---------
   !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2
   !!1 DO iy=iy0+1,iy1
   !!1   iz=iz0
   !!1   DO WHILE(.NOT.ASSOCIATED(Faces_XY(ix1,iy,iz)%Face).AND.iz/=iz1)
   !!1      !           "  RandFace_XY unter Berg east"
   !!1      !Write(*,*) "ib=",ib, " Type-oe:  ,ix1+1=",ix1+1, "iy=",iy,"iz=",i, &
   !!1      iz=iz+1
   !!1   END DO 
   !!1   DO While(ASSOCIATED(Faces_XY(ix1,iy,iz)%Face).AND.iz/=iz1+1)
   !!1      !Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face=>Floor(ib)%Faces_XY(ix1,iy,iz)%Face
   !!1      ALLOCATE(Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face)
   !!1      Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%Vol=Floor(ib)%Faces_XY(ix1,iy,iz)%Face%Vol
   !!1      Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%MidPoint= &
   !!1           Floor(ib)%Faces_XY(ix1,iy,iz)%Face%MidPoint
   !!1      Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%MidPoint%x= &
   !!1            Floor(ib)%Faces_XY(ix1,iy,iz)%Face%MidPoint%x+dx(ix1)
   !!1      Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%in_out=Floor(ib)%Faces_XY(ix1,iy,iz)%Face%in_out
   !!1      Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%ec=Floor(ib)%Faces_XY(ix1,iy,iz)%Face%ec
   !!1      !IF(Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%Vol>0.0d0) THEN
   !!1      !Aenderung: Betreff Allocate EdgeX,EdgeY Face_XY,wenn iz==1 ?Laufvariable in InitAllEdge
   !!1      !?Fall in DurranHatS800.grid, wenn ix=1...,iy=1 iz=1 Vertex%in_out==1  sind 
   !!1      IF((Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%Vol>0.0d0 .AND. &
   !!1        Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%in_out<4)       .OR. &
   !!1        (Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%Vol==0.0d0 .AND. &
   !!1        Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%in_out==0) ) THEN
   !!1        !Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%numberVert==4) ) THEN
   !!1        !Beachte Grenzfall: Face_XY Grenzfläche
   !!1        NrR_FXY=NrR_FXY+1
   !!1        IF(Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%ec>0) THEN
   !!1          IF(Floor(ib)%Faces_XY(ix1,iy,iz)%Face%mp>0) THEN
   !!1            !Zusatzabfrage notwendig,da allokierte Faces_XY unterhalb Berg sein können
   !!1            NrMP_FacesXY=NrMP_FacesXY+1
   !!1            Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%mp=NrMP_FacesXY
   !!1          ELSE
   !!1            Floor(ib)%Faces_XY(ix1+1,iy,iz)%Face%mp=0
   !!1          END IF
   !!1        END IF
   !!1      END IF
   !!1      !Write(*,*) "ib=",ib, " Type-oe:  ix1+1=",ix1+1,"iy=",iy, "iz=",iz, &
   !!1      !           " Schnitt RandFace_XY east"
   !!1      iz=iz+1
   !!1   END DO
   !!1 END DO
   !----Faces_YZ---------
   DO iy=iy0+1,iy1
     DO iz=iz0+1,iz1
        !DO WHILE(.NOT.ASSOCIATED(Faces_YZ(ix1,iy,iz)%Face).AND.iz/=iz1)
         !           "  RandFace_YZ unter Berg east"
         !Write(*,*) "ib=",ib, " Type-oe:  ,ix1+1=",ix1+1, "iy=",iy,"iz=",i, &
        ! iz=iz+1
        !END DO 
        !DO While(ASSOCIATED(Faces_YZ(ix1,iy,iz)%Face).AND.iz/=iz1+1)
        IF(ASSOCIATED(Faces_YZ(ix1,iy,iz)%Face).AND.iz/=iz1+1) THEN
        !IF(Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%ec>0) THEN
         IF((Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%in_out<4.AND. &
             Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%ec>0) &
            .OR. &
            (Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%in_out<0.AND. &
             Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%ec==-1) & 
            .OR. &
            Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%Vol==0.0d0 ) THEN
          ALLOCATE(Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face)
          Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face%Vol=Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%Vol
          Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face%MidPoint= &
                Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%MidPoint
          Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face%MidPoint%x= &
                Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%MidPoint%x+dx(ix1)
          Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face%in_out=Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%in_out
          Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face%ec=Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%ec
          NrR_FYZ=NrR_FYZ+1
          IF(Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face%Vol>0.0d0) THEN
            IF(Floor(ib)%Faces_YZ(ix1,iy,iz)%Face%mp>0) THEN
              !Zusatzabfrage notwendig,da allokierte Faces_YZ unterhalb Berg sein können
              NrMP_FacesYZ=NrMP_FacesYZ+1
              Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face%mp=NrMP_FacesYZ
            ELSE
              Floor(ib)%Faces_YZ(ix1+1,iy,iz)%Face%mp=0
            END IF
          END IF
        END IF
        !Write(*,*) "ib=",ib, " Type-oe:  ix1+1=",ix1+1,"iy=",iy, "iz=",iz, &
        !           " Schnitt RandFace_YZ east"
        !iz=iz+1
        END IF
      END DO !iz
   END DO  !iy
   !!1 !----Faces_ZX---------
   !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2
   !!1 DO iy=iy0,iy1
   !!1  iz=iz0+1
   !!1  DO WHILE(.NOT.ASSOCIATED(Faces_ZX(ix1,iy,iz)%Face).AND.iz/=iz1)
   !!1     !           "  RandFace_ZX unter Berg east"
   !!1     !Write(*,*) "ib=",ib, " Type-oe:  ,ix1+1=",ix1+1, "iy=",iy,"iz=",i, &
   !!1     iz=iz+1
   !!1  END DO 
   !!1  DO While(ASSOCIATED(Faces_ZX(ix1,iy,iz)%Face).AND.iz/=iz1+1)
   !!1      !Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face=>Floor(ib)%Faces_ZX(ix1,iy,iz)%Face
   !!1      ALLOCATE(Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face)
   !!1      Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%Vol=Floor(ib)%Faces_ZX(ix1,iy,iz)%Face%Vol
   !!1      Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%MidPoint= &
   !!1            Floor(ib)%Faces_ZX(ix1,iy,iz)%Face%MidPoint
   !!1      Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%MidPoint%x= &
   !!1            Floor(ib)%Faces_ZX(ix1,iy,iz)%Face%MidPoint%x+dx(ix1)
   !!1      Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%in_out=Floor(ib)%Faces_ZX(ix1,iy,iz)%Face%in_out
   !!1      Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%ec=Floor(ib)%Faces_ZX(ix1,iy,iz)%Face%ec
   !!1      !IF(Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%ec>0) THEN
   !!1      IF(Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%Vol>0.0d0) THEN
   !!1        NrR_FZX=NrR_FZX+1
   !!1        IF(Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%ec>0) THEN
   !!1          IF(Floor(ib)%Faces_ZX(ix1,iy,iz)%Face%mp>0) THEN
   !!1            !Zusatzabfrage notwendig,da allokierte Faces_ZX unterhalb Berg sein können
   !!1            NrMP_FacesZX=NrMP_FacesZX+1
   !!1            Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%mp=NrMP_FacesZX
   !!1          ELSE
   !!1            Floor(ib)%Faces_ZX(ix1+1,iy,iz)%Face%mp=0
   !!1          END IF
   !!1        END IF
   !!1      END IF
   !!1      !Write(*,*) "ib=",ib, " Type-oe:  ix1+1=",ix1+1,"iy=",iy, "iz=",iz, &
   !!1      !           " Schnitt RandFaces_ZX east"
   !!1      iz=iz+1
   !!1  END DO
   !!1 END DO
   !!1 Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+NrR_FXY
   Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+NrR_FYZ
   !!1 Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+NrR_FZX
   !!1 !Floor(ib)%NrMP_FacesXY=Floor(ib)%NrMP_FacesXY+NrMP_FXY
   !!1 !Floor(ib)%NrMP_FacesYZ=Floor(ib)%NrMP_FacesYZ+NrMP_FYZ
   !!1 !Floor(ib)%NrMP_FacesZX=Floor(ib)%NrMP_FacesZX+NrMP_FZX
   !!1 Floor(ib)%NrW_FacesXY=NrW_FacesXY+NrR_FXY
   Floor(ib)%NrW_FacesYZ=NrW_FacesYZ+NrR_FYZ
   !!1 Floor(ib)%NrW_FacesZX=NrW_FacesZX+NrR_FZX

END SUBROUTINE EvalBorderAllFace_OutE


SUBROUTINE EvalBorderAllFace_OutN(ib,ix0,ix1,iy1,iz0,iz1)
   INTEGER :: ib,ix0,ix1,iy1,iz0,iz1

   INTEGER :: ix,iy,iz
   INTEGER :: NrR_FXY,NrR_FYZ,NrR_FZX,NrMP_FXY,NrMP_FYZ,NrMP_FZX
   !!1 !......................................
   !!1 !Für BorderNorth -Face FV (ZX) aktiv
   !!1 !......................................
   !CALL Set(Floor(ib))
   CALL Set_EnvBlkFace(Floor(ib),ib)
   NrR_FXY=0;NrR_FYZ=0;NrR_FZX=0
   NrMP_FXY=0;NrMP_FYZ=0;NrMP_FZX=0  ! noch im Program umsetzen
   !!1 !----Faces_XY---------
   !!1 DO ix=ix0+1,ix1
   !!1  !!DO iz=iz0+1,iz1
   !!1  iz=iz0
   !!1  DO WHILE(.NOT.ASSOCIATED(Faces_XY(ix,iy1,iz)%Face).AND.iz/=iz1)
   !!1      !Write(*,*) "ib=",ib, " Type-on:  ix=",ix," iy1=",iy1," iz=",iz, &
   !!1      !           "  RandFaces_XY unter Berg north"
   !!1     !!ALLOCATE(Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face)
   !!1      !!Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%in_out=-8
   !!1      !!Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%Vol=0.0d0
   !!1      !NrRN_FXY=NrRN_FXY+1
   !!1      iz=iz+1
   !!1  END DO 
   !!1  DO WHILE(ASSOCIATED(Faces_XY(ix,iy1,iz)%Face).AND.iz/=iz1+1)
   !!1      !Floor(ib)%Faces_XY(ix,iy1,iz)%Face=>Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face
   !!1      ALLOCATE(Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face)
   !!1      Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%Vol=Floor(ib)%Faces_XY(ix,iy1,iz)%Face%Vol
   !!1      Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%MidPoint= &
   !!1                              Floor(ib)%Faces_XY(ix,iy1,iz)%Face%MidPoint
   !!1      Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%MidPoint%y= &
   !!1                              Floor(ib)%Faces_XY(ix,iy1,iz)%Face%MidPoint%y+dy(iy1)
   !!1      Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%in_out=Floor(ib)%Faces_XY(ix,iy1,iz)%Face%in_out
   !!1      Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%ec=Floor(ib)%Faces_XY(ix,iy1,iz)%Face%ec
   !!1      !IF(Floor(ib)%Faces_XY(ix,iy1,iz)%Face%Vol>0.0d0) THEN
   !!1      IF((Floor(ib)%Faces_XY(ix,iy1,iz)%Face%Vol>0.0d0 .AND. &
   !!1          Floor(ib)%Faces_XY(ix,iy1,iz)%Face%in_out<4)  .OR. &
   !!1         (Floor(ib)%Faces_XY(ix,iy1,iz)%Face%Vol==0.0d0 .AND. &
   !!1          Floor(ib)%Faces_XY(ix,iy1,iz)%Face%in_out==0) ) THEN
   !!1         NrR_FXY=NrR_FXY+1
   !!1         IF(Floor(ib)%Faces_XY(ix,iy1,iz)%Face%ec>0) THEN
   !!1           IF(Floor(ib)%Faces_XY(ix,iy1,iz)%Face%mp>0) THEN
   !!1             !Zusatzabfrage notwendig,da allokierte Faces_XY unterhalb Berg sein können
   !!1             NrMP_FXY=NrMP_FXY+1
   !!1             Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%mp=NrMP_FXY
   !!1           ELSE
   !!1             Floor(ib)%Faces_XY(ix,iy1+1,iz)%Face%mp=0
   !!1           END IF
   !!1         END IF
   !!1      END IF
   !!1      !Write(*,*) "ib=",ib, " Type-on:  ix=",ix," iy1=",iy1," iz=",iz, &
   !!1      !           " Schnitt RandFaces_XY north"
   !!1      iz=iz+1
   !!1  END DO 
   !!1 END DO
   !!1 !----Faces_YZ---------
   !!1 DO ix=ix0,ix1
   !!1  !!DO iz=iz0+1,iz1
   !!1  iz=iz0+1
   !!1  DO WHILE(.NOT.ASSOCIATED(Faces_YZ(ix,iy1,iz)%Face).AND.iz/=iz1)
   !!1      !Write(*,*) "ib=",ib, " Type-on:  ix=",ix," iy1=",iy1," iz=",iz, &
   !!1      !           "  RandFaces_YZ unter Berg north"
   !!1     !!ALLOCATE(Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face)
   !!1      !!Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%in_out=-8
   !!1      !!Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%Vol=0.0d0
   !!1      !NrRN_FYZ=NrRN_FYZ+1
   !!1      iz=iz+1
   !!1  END DO 
   !!1  DO WHILE(ASSOCIATED(Faces_YZ(ix,iy1,iz)%Face).AND.iz/=iz1+1)
   !!1      !Floor(ib)%Faces_YZ(ix,iy1,iz)%Face=>Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face
   !!1      ALLOCATE(Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face)
   !!1      Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%Vol=Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%Vol
   !!1      Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%MidPoint= &
   !!1                              Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%MidPoint
   !!1      Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%MidPoint%y= &
   !!1                              Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%MidPoint%y+dy(iy1)
   !!1      Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%in_out=Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%in_out
   !!1      Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%ec=Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%ec
   !!1      !IF(Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%Vol>0.0d0) THEN
   !!1      IF((Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%Vol>0.0d0 .AND. &
   !!1          Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%in_out<4)     .OR. &
   !!1         (Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%Vol==0.0d0 .AND. &
   !!1          Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%in_out==0) ) THEN 
   !!1         NrR_FYZ=NrR_FYZ+1
   !!1         IF(Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%ec>0) THEN
   !!1           IF(Floor(ib)%Faces_YZ(ix,iy1,iz)%Face%mp>0) THEN
   !!1             !Zusatzabfrage notwendig,da allokierte Faces_YZ unterhalb Berg sein können
   !!1             NrMP_FYZ=NrMP_FYZ+1
   !!1             Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%mp=NrMP_FYZ
   !!1           ELSE
   !!1             Floor(ib)%Faces_YZ(ix,iy1+1,iz)%Face%mp=0
   !!1           END IF
   !!1         END IF
   !!1      END IF
   !!1      !Write(*,*) "ib=",ib, " Type-on:  ix=",ix," iy1=",iy1," iz=",iz, &
   !!1      !           " Schnitt RandFaces_YZ north"
   !!1      iz=iz+1
   !!1  END DO 
   !!1 END DO
   !----Faces_ZX---------
   DO ix=ix0+1,ix1
      DO iz=iz0+1,iz1
      !iz=iz0+1
      !DO WHILE(.NOT.ASSOCIATED(Faces_ZX(ix,iy1,iz)%Face).AND.iz/=iz1)
      !    !Write(*,*) "ib=",ib, " Type-on:  ix=",ix," iy1=",iy1," iz=",iz, &
      !    !           "  RandFaces_ZX unter Berg north"
      !   !!ALLOCATE(Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face)
      !    !!Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%in_out=-8
      !    !!Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%Vol=0.0d0
      !    !NrRN_FZX=NrRN_FZX+1
      !    iz=iz+1
      !END DO 
      !DO WHILE(ASSOCIATED(Faces_ZX(ix,iy1,iz)%Face).AND.iz/=iz1+1)
      IF (ASSOCIATED(Faces_ZX(ix,iy1,iz)%Face).AND.iz/=iz1+1) THEN
          !IF(Faces_ZX(ix,iy1,iz)%Face%ec>0) THEN
          IF((Faces_ZX(ix,iy1,iz)%Face%in_out<4.AND.Faces_ZX(ix,iy1,iz)%Face%ec>0) &
             .OR. &
             (Faces_ZX(ix,iy1,iz)%Face%in_out<0.AND.Faces_ZX(ix,iy1,iz)%Face%ec==-1) &
             .OR. &
             Faces_ZX(ix,iy1,iz)%Face%Vol==0.0d0) THEN
          ALLOCATE(Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face)
          Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%Vol=Floor(ib)%Faces_ZX(ix,iy1,iz)%Face%Vol
          Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%MidPoint= &
                                  Floor(ib)%Faces_ZX(ix,iy1,iz)%Face%MidPoint
          Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%MidPoint%y= &
                                  Floor(ib)%Faces_ZX(ix,iy1,iz)%Face%MidPoint%y+dy(iy1)
          Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%in_out=Floor(ib)%Faces_ZX(ix,iy1,iz)%Face%in_out
          Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%ec=Floor(ib)%Faces_ZX(ix,iy1,iz)%Face%ec
          NrR_FZX=NrR_FZX+1
          IF(Floor(ib)%Faces_ZX(ix,iy1,iz)%Face%Vol>0.0d0) THEN
            IF(Floor(ib)%Faces_ZX(ix,iy1,iz)%Face%mp>0) THEN
              !Zusatzabfrage notwendig,da allokierte Faces_ZX unterhalb Berg sein können
              NrMP_FZX=NrMP_FZX+1
              Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%mp=NrMP_FZX
            ELSE
              Floor(ib)%Faces_ZX(ix,iy1+1,iz)%Face%mp=0
            END IF
          END IF
          !Write(*,*) "ib=",ib, " Type-on:  ix=",ix," iy1=",iy1," iz=",iz, &
          !           " Schnitt RandFaces_ZX north"
        END IF
        !iz=iz+1
      END IF
      END DO  !iz 
   END DO   !ix
   !!1 Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+NrR_FXY
   !!1 Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+NrR_FYZ
   Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+NrR_FZX
   !!1 !Floor(ib)%NrMP_FacesXY=Floor(ib)%NrMP_FacesXY+NrMP_FXY
   !!1 !Floor(ib)%NrMP_FacesYZ=Floor(ib)%NrMP_FacesYZ+NrMP_FYZ
   !!1 !Floor(ib)%NrMP_FacesZX=Floor(ib)%NrMP_FacesZX+NrMP_FZX
   !!1 Floor(ib)%NrW_FacesXY=NrW_FacesXY+NrR_FXY
   !!1 Floor(ib)%NrW_FacesYZ=NrW_FacesYZ+NrR_FYZ
   Floor(ib)%NrW_FacesZX=NrW_FacesZX+NrR_FZX
END SUBROUTINE EvalBorderAllFace_OutN

SUBROUTINE EvalBorderAllFace_OutS(ib,ix0,ix1,iy0,iz0,iz1)
   INTEGER :: ib,ix0,ix1,iy0,iz0,iz1

   INTEGER :: ix,iy,iz
   INTEGER :: NrR_FXY,NrR_FYZ,NrR_FZX,NrMP_FXY,NrMP_FYZ,NrMP_FZX
   !!1 !......................................
   !!1 !Für BorderSouth -Face FV (ZX) aktiv
   !!1 !......................................
   !CALL Set(Floor(ib))
   CALL Set_EnvBlkFace(Floor(ib),ib)
   NrR_FXY=0;NrR_FYZ=0;NrR_FZX=0
   NrMP_FXY=0;NrMP_FYZ=0;NrMP_FZX=0  ! noch im Program umsetzen
   !!1 !----Faces_XY---------
   !!1 DO ix=ix0+1,ix1
   !!1  !DO iz=iz0+1,iz1
   !!1  iz=iz0
   !!1  DO WHILE(.NOT.ASSOCIATED(Faces_XY(ix,iy0+1,iz)%Face).AND.iz/=iz1)
   !!1     !Write(*,*) "ib=",ib, " Type-os:  ,ix=",ix, "iy0+1=",iy0+1,"iz=",iz, &
   !!1     !           "  RandFaces_XY unter Berg south"
   !!1     !!ALLOCATE(Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face)
   !!1     !!Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%in_out=-8
   !!1     !!Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%Vol=0.0d0
   !!1     !NrRN_FXY=NrRN_FXY+1
   !!1     iz=iz+1
   !!1  END DO 
   !!1  DO While(ASSOCIATED(Faces_XY(ix,iy0+1,iz)%Face).AND.iz/=iz1+1)
   !!1      !Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face=>Floor(ib)%Faces_XY(ix,iy0,iz)%Face
   !!1      ALLOCATE(Floor(ib)%Faces_XY(ix,iy0,iz)%Face)
   !!1      Floor(ib)%Faces_XY(ix,iy0,iz)%Face%Vol=Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%Vol
   !!1      Floor(ib)%Faces_XY(ix,iy0,iz)%Face%MidPoint= &
   !!1                          Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%MidPoint
   !!1      Floor(ib)%Faces_XY(ix,iy0,iz)%Face%MidPoint%y= &
   !!1                          Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%MidPoint%y-dy(iy0+1)
   !!1      Floor(ib)%Faces_XY(ix,iy0,iz)%Face%in_out=Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%in_out
   !!1      Floor(ib)%Faces_XY(ix,iy0,iz)%Face%ec=Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%ec
   !!1      !IF(Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%Vol>0.0d0) THEN
   !!1      IF((Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%Vol>0.0d0 .AND. &
   !!1          Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%in_out<4)  .OR. &
   !!1         (Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%Vol==0.0d0 .AND. &
   !!1          Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%in_out==0) ) THEN ! .OR. &
   !!1         !(Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%numberVert==4 .AND. &
   !!1         ! Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%in_out==-2) ) THEN
   !!1        NrR_FXY=NrR_FXY+1
   !!1        IF(Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%ec>0) THEN
   !!1          IF(Floor(ib)%Faces_XY(ix,iy0+1,iz)%Face%mp>0) THEN
   !!1            !Zusatzabfrage notwendig,da allokierte Faces_XY unterhalb Berg sein können
   !!1            NrMP_FXY=NrMP_FXY+1
   !!1            Floor(ib)%Faces_XY(ix,iy0,iz)%Face%mp=NrMP_FXY
   !!1          ELSE
   !!1            Floor(ib)%Faces_XY(ix,iy0,iz)%Face%mp=0
   !!1          END IF
   !!1        END IF
   !!1      END IF
   !!1      !Write(*,*) "ib=",ib, " Type-os:  ix=",ix,"iy0+1=",iy0+1, "iz=",iz, &
   !!1      !           " Schnitt RandFaces_XY south"
   !!1      iz=iz+1
   !!1  END DO
   !!1 END DO
   !!1 !----Faces_YZ---------
   !!1 DO ix=ix0,ix1
   !!1  !DO iz=iz0+1,iz1
   !!1  iz=iz0+1
   !!1  DO WHILE(.NOT.ASSOCIATED(Faces_YZ(ix,iy0+1,iz)%Face).AND.iz/=iz1)
   !!1     !Write(*,*) "ib=",ib, " Type-os:  ,ix=",ix, "iy0+1=",iy0+1,"iz=",iz, &
   !!1     !           "  RandFaces_YZ unter Berg south"
   !!1     !!ALLOCATE(Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face)
   !!1     !!Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%in_out=-8
   !!1     !!Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%Vol=0.0d0
   !!1     !NrRN_FYZ=NrRN_FYZ+1
   !!1     iz=iz+1
   !!1  END DO 
   !!1  DO While(ASSOCIATED(Faces_YZ(ix,iy0+1,iz)%Face).AND.iz/=iz1+1)
   !!1      !Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face=>Floor(ib)%Faces_YZ(ix,iy0,iz)%Face
   !!1      ALLOCATE(Floor(ib)%Faces_YZ(ix,iy0,iz)%Face)
   !!1      Floor(ib)%Faces_YZ(ix,iy0,iz)%Face%Vol=Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%Vol
   !!1      Floor(ib)%Faces_YZ(ix,iy0,iz)%Face%MidPoint= &
   !!1                          Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%MidPoint
   !!1      Floor(ib)%Faces_YZ(ix,iy0,iz)%Face%MidPoint%y= &
   !!1                          Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%MidPoint%y-dy(iy0+1)
   !!1      Floor(ib)%Faces_YZ(ix,iy0,iz)%Face%in_out=Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%in_out
   !!1      Floor(ib)%Faces_YZ(ix,iy0,iz)%Face%ec=Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%ec
   !!1      !IF(Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%Vol>0.0d0) THEN
   !!1      IF((Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%Vol>0.0d0 .AND. &
   !!1          Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%in_out<4)     .OR. &
   !!1         (Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%Vol==0.0d0 .AND. &
   !!1          Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%in_out==0) ) THEN 
   !!1        NrR_FYZ=NrR_FYZ+1
   !!1        IF(Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%ec>0) THEN
   !!1          IF(Floor(ib)%Faces_YZ(ix,iy0+1,iz)%Face%mp>0) THEN
   !!1            !Zusatzabfrage notwendig,da allokierte Faces_YZ unterhalb Berg sein können
   !!1            NrMP_FYZ=NrMP_FYZ+1
   !!1            Floor(ib)%Faces_YZ(ix,iy0,iz)%Face%mp=NrMP_FYZ
   !!1          ELSE
   !!1            Floor(ib)%Faces_YZ(ix,iy0,iz)%Face%mp=0
   !!1          END IF
   !!1        END IF
   !!1      END IF
   !!1      !Write(*,*) "ib=",ib, " Type-os:  ix=",ix,"iy0+1=",iy0+1, "iz=",iz, &
   !!1      !           " Schnitt RandFaces_YZ south"
   !!1      iz=iz+1
   !!1  END DO
   !!1 END DO
   !----Faces_ZX---------
   DO ix=ix0+1,ix1
      DO iz=iz0+1,iz1
      !iz=iz0+1
      !DO WHILE(.NOT.ASSOCIATED(Faces_ZX(ix,iy0,iz)%Face).AND.iz/=iz1)
      !   !Write(*,*) "ib=",ib, " Type-os:  ,ix=",ix, "iy0+1=",iy0+1,"iz=",iz, &
      !   !           "  RandFaces_ZX unter Berg south"
      !   !!ALLOCATE(Floor(ib)%Faces_ZX(ix,iy0+1,iz)%Face)
      !   !!Floor(ib)%Faces_ZX(ix,iy0+1,iz)%Face%in_out=-8
      !   !!Floor(ib)%Faces_ZX(ix,iy0+1,iz)%Face%Vol=0.0d0
      !   !NrRN_FZX=NrRN_FZX+1
      !   iz=iz+1
      !END DO 
      !DO While(ASSOCIATED(Faces_ZX(ix,iy0,iz)%Face).AND.iz/=iz1+1)
      IF (ASSOCIATED(Faces_ZX(ix,iy0,iz)%Face).AND.iz/=iz1+1) THEN
        !IF(Faces_ZX(ix,iy0,iz)%Face%ec>0) THEN
        IF((Faces_ZX(ix,iy0,iz)%Face%in_out<4.AND.Faces_ZX(ix,iy0,iz)%Face%ec>0) &
           .OR. &
           (Faces_ZX(ix,iy0,iz)%Face%in_out<0.AND.Faces_ZX(ix,iy0,iz)%Face%ec==-1) & 
           .OR. &
           Faces_ZX(ix,iy0,iz)%Face%Vol==0.0d0) THEN
          ALLOCATE(Floor(ib)%Faces_ZX(ix,iy0-1,iz)%Face)
          Floor(ib)%Faces_ZX(ix,iy0-1,iz)%Face%Vol=Floor(ib)%Faces_ZX(ix,iy0,iz)%Face%Vol
          Floor(ib)%Faces_ZX(ix,iy0-1,iz)%Face%MidPoint= &
                              Floor(ib)%Faces_ZX(ix,iy0,iz)%Face%MidPoint
          Floor(ib)%Faces_ZX(ix,iy0-1,iz)%Face%MidPoint%y= &
                              Floor(ib)%Faces_ZX(ix,iy0,iz)%Face%MidPoint%y-dy(iy0+1)
          Floor(ib)%Faces_ZX(ix,iy0-1,iz)%Face%in_out=Floor(ib)%Faces_ZX(ix,iy0,iz)%Face%in_out
          Floor(ib)%Faces_ZX(ix,iy0-1,iz)%Face%ec=Floor(ib)%Faces_ZX(ix,iy0,iz)%Face%ec
          NrR_FZX=NrR_FZX+1 !alle allokierten auch wenn Vol=0 (Grenze Point/Edges)
          IF(Floor(ib)%Faces_ZX(ix,iy0,iz)%Face%Vol>0.0d0) THEN
            !NrR_FZX=NrR_FZX+1
            IF(Floor(ib)%Faces_ZX(ix,iy0,iz)%Face%mp>0) THEN
              !Zusatzabfrage notwendig,da allokierte Faces_ZX unterhalb Berg sein können
              NrMP_FZX=NrMP_FZX+1
              Floor(ib)%Faces_ZX(ix,iy0-1,iz)%Face%mp=NrMP_FZX
            ELSE
              Floor(ib)%Faces_ZX(ix,iy0-1,iz)%Face%mp=0
            END IF
          END IF
          !Write(*,*) "ib=",ib, " Type-os:  ix=",ix,"iy0+1=",iy0+1, "iz=",iz, &
          !           " Schnitt RandFaces_ZX south"
        END IF
        !iz=iz+1
      END IF
      END DO !iz
   END DO  !ix
   !!1 Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+NrR_FXY
   !!1 Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+NrR_FYZ
   Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+NrR_FZX
   !!1 !Floor(ib)%NrMP_FacesXY=Floor(ib)%NrMP_FacesXY+NrMP_FXY
   !!1 !Floor(ib)%NrMP_FacesYZ=Floor(ib)%NrMP_FacesYZ+NrMP_FYZ
   !!1 !Floor(ib)%NrMP_FacesZX=Floor(ib)%NrMP_FacesZX+NrMP_FZX
   !!1 Floor(ib)%NrW_FacesXY=NrW_FacesXY+NrR_FXY
   !!1 Floor(ib)%NrW_FacesYZ=NrW_FacesYZ+NrR_FYZ
   Floor(ib)%NrW_FacesZX=NrW_FacesZX+NrR_FZX

END SUBROUTINE EvalBorderAllFace_OutS


SUBROUTINE EvalBorderAllFace_OutT(ib,ix0,ix1,iy0,iy1,iz1)
   INTEGER :: ib,ix0,ix1,iy0,iy1,iz1

   INTEGER :: ix,iy,iz
   INTEGER :: NrR_FXY,NrR_FYZ,NrR_FZX,NrMP_FXY,NrMP_FYZ,NrMP_FZX
   !!1 !......................................
   !!1 !Für BorderTop -Face FW (XY) aktiv
   !!1 !......................................

   !CALL Set(Floor(ib))
   CALL Set_EnvBlkFace(Floor(ib),ib)
   NrR_FXY=0;NrR_FYZ=0;NrR_FZX=0
   NrMP_FXY=0;NrMP_FYZ=0;NrMP_FZX=0  ! noch im Program umsetzen
   !----Faces_XY---------
   DO ix=ix0+1,ix1
     DO iy=iy0+1,iy1
       IF (ASSOCIATED(Faces_XY(ix,iy,iz1)%Face)) THEN
          !NrR_FXY=NrR_FXY+1
          !IF(Faces_XY(ix,iy,iz1)%Face%ec>0) THEN
          !FolgeAbfrage dem Weight_Check angeglichen
           IF((Faces_XY(ix,iy,iz1)%Face%in_out<4.AND.Faces_XY(ix,iy,iz1)%Face%ec>0) &
              .OR. &
              (Faces_XY(ix,iy,iz1)%Face%in_out<0.AND.Faces_XY(ix,iy,iz1)%Face%ec==-1) &
              .OR. &
              Faces_XY(ix,iy,iz1)%Face%Vol==0.0)  THEN
            ALLOCATE(Floor(ib)%Faces_XY(ix,iy,iz1+1)%Face)
            Floor(ib)%Faces_XY(ix,iy,iz1+1)%Face%Vol=Floor(ib)%Faces_XY(ix,iy,iz1)%Face%Vol
            Floor(ib)%Faces_XY(ix,iy,iz1+1)%Face%MidPoint= &
                                    Floor(ib)%Faces_XY(ix,iy,iz1)%Face%MidPoint
            Floor(ib)%Faces_XY(ix,iy,iz1+1)%Face%MidPoint%z= &
                                    Floor(ib)%Faces_XY(ix,iy,iz1)%Face%MidPoint%z+dz(iz1)
            Floor(ib)%Faces_XY(ix,iy,iz1+1)%Face%in_out=Floor(ib)%Faces_XY(ix,iy,iz1)%Face%in_out
            Floor(ib)%Faces_XY(ix,iy,iz1+1)%Face%ec=Floor(ib)%Faces_XY(ix,iy,iz1)%Face%ec
            NrR_FXY=NrR_FXY+1
            NrMP_FXY=NrMP_FXY+1
            Floor(ib)%Faces_XY(ix,iy,iz1+1)%Face%mp=NrMP_FXY
            !Write(*,*) "ib=",ib, " Type-ot:  ix=",ix," iy=",iy," iz1=",iz1, &
            !           " Schnitt RandFacesXY top"
          END IF
       END IF
     END DO 
   END DO
   !!1 !----Faces_YZ---------
   !!1 DO ix=ix0,ix1
   !!1  DO iy=iy0+1,iy1
   !!1    IF (ASSOCIATED(Faces_YZ(ix,iy,iz1)%Face)) THEN
   !!1       !Floor(ib)%Faces_YZ(ix,iy,iz1+1)%Face=>Floor(ib)%Faces_YZ(ix,iy,iz1)%Face
   !!1       IF(Faces_YZ(ix,iy,iz1)%Face%ec>0) THEN
   !!1         ALLOCATE(Floor(ib)%Faces_YZ(ix,iy,iz1+1)%Face)
   !!1         Floor(ib)%Faces_YZ(ix,iy,iz1+1)%Face%Vol=Floor(ib)%Faces_YZ(ix,iy,iz1)%Face%Vol
   !!1         Floor(ib)%Faces_YZ(ix,iy,iz1+1)%Face%MidPoint= &
   !!1                                 Floor(ib)%Faces_YZ(ix,iy,iz1)%Face%MidPoint
   !!1         Floor(ib)%Faces_YZ(ix,iy,iz1+1)%Face%MidPoint%z= &
   !!1                                 Floor(ib)%Faces_YZ(ix,iy,iz1)%Face%MidPoint%z+dz(iz1)
   !!1         Floor(ib)%Faces_YZ(ix,iy,iz1+1)%Face%in_out=Floor(ib)%Faces_YZ(ix,iy,iz1)%Face%in_out
   !!1         Floor(ib)%Faces_YZ(ix,iy,iz1+1)%Face%ec=Floor(ib)%Faces_YZ(ix,iy,iz1)%Face%ec
   !!1         NrR_FYZ=NrR_FYZ+1
   !!1         NrMP_FYZ=NrMP_FYZ+1
   !!1         Floor(ib)%Faces_YZ(ix,iy,iz1+1)%Face%mp=NrMP_FYZ
   !!1         !Write(*,*) "ib=",ib, " Type-ot:  ix=",ix," iy=",iy," iz1=",iz1, &
   !!1         !           " Schnitt RandFacesYZ top"
   !!1       END IF
   !!1    END IF
   !!1  END DO 
   !!1 END DO
   !!1 !----Faces_ZX---------
   !!1 DO ix=ix0+1,ix1
   !!1  DO iy=iy0,iy1
   !!1    IF (ASSOCIATED(Faces_ZX(ix,iy,iz1)%Face)) THEN
   !!1       !Floor(ib)%Faces_ZX(ix,iy,iz1+1)%Face=>Floor(ib)%Faces_ZX(ix,iy,iz1)%Face
   !!1       IF(Faces_ZX(ix,iy,iz1)%Face%ec>0) THEN
   !!1         ALLOCATE(Floor(ib)%Faces_ZX(ix,iy,iz1+1)%Face)
   !!1         Floor(ib)%Faces_ZX(ix,iy,iz1+1)%Face%Vol=Floor(ib)%Faces_ZX(ix,iy,iz1)%Face%Vol
   !!1         Floor(ib)%Faces_ZX(ix,iy,iz1+1)%Face%MidPoint= &
   !!1                                 Floor(ib)%Faces_ZX(ix,iy,iz1)%Face%MidPoint
   !!1         Floor(ib)%Faces_ZX(ix,iy,iz1+1)%Face%MidPoint%z= &
   !!1                                 Floor(ib)%Faces_ZX(ix,iy,iz1)%Face%MidPoint%z+dz(iz1)
   !!1         Floor(ib)%Faces_ZX(ix,iy,iz1+1)%Face%in_out=Floor(ib)%Faces_ZX(ix,iy,iz1)%Face%in_out
   !!1         Floor(ib)%Faces_ZX(ix,iy,iz1+1)%Face%ec=Floor(ib)%Faces_ZX(ix,iy,iz1)%Face%ec
   !!1         NrR_FZX=NrR_FZX+1
   !!1         NrMP_FZX=NrMP_FZX+1
   !!1         Floor(ib)%Faces_ZX(ix,iy,iz1+1)%Face%mp=NrMP_FZX
   !!1         !Write(*,*) "ib=",ib, " Type-ot:  ix=",ix," iy=",iy," iz1=",iz1, &
   !!1         !           " Schnitt RandFacesZX top"
   !!1       END IF
   !!1    END IF
   !!1  END DO 
   !!1 END DO
   Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+NrR_FXY
   !!1 Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+NrR_FYZ
   !!1 Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+NrR_FZX
   !!1 !Floor(ib)%NrMP_FacesXY=Floor(ib)%NrMP_FacesXY+NrMP_FXY
   !!1 !Floor(ib)%NrMP_FacesYZ=Floor(ib)%NrMP_FacesYZ+NrMP_FYZ
   !!1 !Floor(ib)%NrMP_FacesZX=Floor(ib)%NrMP_FacesZX+NrMP_FZX
   Floor(ib)%NrW_FacesXY=NrW_FacesXY+NrR_FXY
   !!1 Floor(ib)%NrW_FacesYZ=NrW_FacesYZ+NrR_FYZ
   !!1 Floor(ib)%NrW_FacesZX=NrW_FacesZX+NrR_FZX
END SUBROUTINE EvalBorderAllFace_OutT


SUBROUTINE EvalBorderAllFace_OutB(ib,ix0,ix1,iy0,iy1,iz0)
   INTEGER :: ib,ix0,ix1,iy0,iy1,iz0

   INTEGER :: ix,iy,iz,in_out
   INTEGER :: NrRN_FXY,NrRN_FYZ,NrRN_FZX
   !!1 !......................................
   !!1 !Für BorderBottom -Face FW (XY) aktiv
   !!1 !......................................

   !CALL Set(Floor(ib))
   CALL Set_EnvBlkFace(Floor(ib),ib)
   NrRN_FXY=0;NrRN_FYZ=0;NrRN_FZX=0
   !----Faces_XY---------
!     iz=iz0
!     DO ix=ix0+1,ix1
!       DO iy=iy0+1,iy1
!           in_out=Vertices(ix-1,iy-1,iz)%in_out &
!                 +Vertices(ix  ,iy-1,iz)%in_out &
!                 +Vertices(ix-1,iy  ,iz)%in_out &
!                 +Vertices(ix  ,iy  ,iz)%in_out
!           !Berg-Grenze  unterhalb z0
!           IF (in_out==4) THEN
!             NrRN_FXY=NrRN_FXY+1   ! Erweiterung iz0
!                                   ! da im Block iz0 nur bis Grenze gezählt!
!           END IF
!         !Write(*,*) "ib=",ib, " Type-ob:  ,ix=",ix, "iy=",iy,"iz=",iz0, &
!         !           "  RandCell unter Berg bottom"
!         !!ALLOCATE(Floor(ib)%Cell(ix,iy,iz0)%Cell)
!         !!Floor(ib)%Cell(ix,iy,iz0)%Cell%in_out=-8
!         !!Floor(ib)%Cell(ix,iy,iz0)%Cell%Vol=0.0d0
!         NrRN_FXY=NrRN_FXY+1   ! iz0-1 Zählung
!       END DO 
!     END DO
   !!1 !----Faces_YZ---------
   !!1 DO ix=ix0,ix1
   !!1  DO iy=iy0+1,iy1
   !!1      !Write(*,*) "ib=",ib, " Type-ob:  ,ix=",ix, "iy=",iy,"iz=",iz0, &
   !!1      !           "  RandCell unter Berg bottom"
   !!1      !!ALLOCATE(Floor(ib)%Cell(ix,iy,iz0)%Cell)
   !!1      !!Floor(ib)%Cell(ix,iy,iz0)%Cell%in_out=-8
   !!1      !!Floor(ib)%Cell(ix,iy,iz0)%Cell%Vol=0.0d0
   !!1      NrRN_FYZ=NrRN_FYZ+1
   !!1  END DO 
   !!1 END DO
   !!1 !----Faces_ZX---------
   !!1 DO ix=ix0+1,ix1
   !!1  DO iy=iy0,iy1
   !!1      !Write(*,*) "ib=",ib, " Type-ob:  ,ix=",ix, "iy=",iy,"iz=",iz0, &
   !!1      !           "  RandCell unter Berg bottom"
   !!1      !!ALLOCATE(Floor(ib)%Cell(ix,iy,iz0)%Cell)
   !!1      !!Floor(ib)%Cell(ix,iy,iz0)%Cell%in_out=-8
   !!1      !!Floor(ib)%Cell(ix,iy,iz0)%Cell%Vol=0.0d0
   !!1      NrRN_FZX=NrRN_FZX+1
   !!1  END DO 
   !!1 END DO
   Floor(ib)%NrRN_FacesXY=Floor(ib)%NrRN_FacesXY+NrRN_FXY
   !!1 Floor(ib)%NrRN_FacesYZ=Floor(ib)%NrRN_FacesYZ+NrRN_FYZ
   !!1 Floor(ib)%NrRN_FacesZX=Floor(ib)%NrRN_FacesZX+NrRN_FZX
   Floor(ib)%NrW_FacesXY=Floor(ib)%NrW_FacesXY+NrRN_FXY
   !!1 Floor(ib)%NrW_FacesYZ=NrW_FacesYZ+NrRN_FYZ
   !!1 Floor(ib)%NrW_FacesZX=NrW_FacesZX+NrRN_FZX
END SUBROUTINE EvalBorderAllFace_OutB


SUBROUTINE CountFaceSubMontane_WE(ib)
   !Zählung Face Rand unterhalb Berg 'w,e'
   INTEGER :: ib

   INTEGER :: i,j,k,in_out,bli 
   INTEGER :: NrRN_FXY,NrRN_FYZ,NrRN_FZX
   !!1 !......................................
   !!1 !Für BorderWestEast -Face FU (YZ) aktiv
   !!1 !......................................

   !CALL Set(Floor(ib))
   CALL Set_EnvBlkFace(Floor(ib),ib)
   NrRN_FXY=0;NrRN_FYZ=0;NrRN_FZX=0
   !!1 IF (nType(2:2)=='w') THEN
   !!1   i=ix0
   !!1 ELSE
   !!1   i=ix1+1
   !!1 END IF
   !!1 !----Faces_XY---------
   !!1 DO j=iy0+1,iy1
   !!1   k=iz0;in_out=0
   !!1   DO WHILE(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face).AND.in_out/=4.AND.k/=iz1+1)
   !!1     in_out=Vertices(i-1,j-1, k  )%in_out &
   !!1          +Vertices(i  ,j-1, k  )%in_out &
   !!1          +Vertices(i-1,j  , k  )%in_out   &
   !!1          +Vertices(i  ,j  , k  )%in_out
   !!1     IF (.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face).AND.in_out<-2) THEN
   !!1         NrRN_FXY=NrRN_FXY+1
   !!1     END IF
   !!1     k=k+1
   !!1   END DO
   !!1 END DO
   !----Faces_YZ---------
   IF (nType(2:2)=='w') THEN
          i=ix0-1
         bli=ix0    !für Grenze Out
    ELSE
          i=ix1+1
          bli=ix1
   END IF
   DO j=iy0+1,iy1
     k=iz0+1
     in_out=Vertices(i,j-1, k-1)%in_out+Vertices(i,j-1, k  )%in_out &
           +Vertices(i,j  , k-1)%in_out+Vertices(i,j  , k  )%in_out
     !DO WHILE(.NOT.ASSOCIATED(Faces_YZ(i,j,k)%Face).AND.in_out/=4.AND.k/=iz1+1)
     DO k=iz0+1,iz1
      IF (.NOT.ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
       in_out=Vertices(i,j-1, k-1)%in_out+Vertices(i,j-1, k  )%in_out &
             +Vertices(i,j  , k-1)%in_out+Vertices(i,j  , k  )%in_out
       !IF ((.NOT.ASSOCIATED(Faces_YZ(i,j,k)%Face).AND.in_out<-2) THEN
       !+special (+1,-1,-1,-1) ?warum kein Schnitt def. in Analyze 
       !Valley3D_Blk7dy1400y10z12.grid unter 
       !/s1agi/reutgen/ASAM/ParExample/Test_Valley3D_V1_Blk, i=1,j=9,k=1
       !IF ((.NOT.ASSOCIATED(Faces_YZ(i,j,k)%Face).AND.in_out<-2) .OR. & 
       !    (in_out==-2 .AND. &
       !    (MAX(Vertices(i,j-1, k-1)%in_out,Vertices(i,j-1, k  )%in_out, &
       !         Vertices(i,j  , k-1)%in_out,Vertices(i,j  , k  )%in_out)==1).AND. &
       !    (MIN(Vertices(i,j-1, k-1)%in_out,Vertices(i,j-1, k  )%in_out, &
       !         Vertices(i,j  , k-1)%in_out,Vertices(i,j  , k  )%in_out)==-1)) ) THEN
       ! zu wenig gezählt Bsp: Valley3D_Blk7u8x5.grid ,Test_Valley3D_V1_Blk 
       !.......
       !IF(Vertices(bli,j-1, k-1)%in_out+Vertices(bli,j-1, k  )%in_out &
       !   +Vertices(bli,j  , k-1)%in_out+Vertices(bli,j  , k  )%in_out >=2) THEN
       !       Write (*,*) " EXIT aus  CountFaceSubMontane_WE"
       !       Write(*,*) " FaceYZ-Block(bli)=maxVol und FaceYZ-(i)Rand ist Vol==0.d0"
       !       Write(*,*) " Randzählung läuft über letztes FaceYZ-Block(blj) mit ec=1"
       !       Write(*,*) "i,j,k=", bli,j,k
       !       EXIT
       !ELSE  ! Valley3D_Blk7u8x4.grid viele zu wenig gezählt
         !IF (in_out<-2) THEN
       !.......
         IF (in_out<=0) THEN ! io Bsp:Valley3D_Blk7u8x5.grid,,Test_Valley3D_V1_Blk
                             ! io  Bsp:Valley3D_Blk7u8x4.grid,,Test_Valley3D_V1_Blk
                            ! zu viel gezählt auch Bsp.:
         !IF (in_out<0) THEN !zu wenig gezählt, Bsp:Valley3D_Blk7u8x5.grid,,Test_Valley3D_V1_Blk 
            NrRN_FYZ=NrRN_FYZ+1
         END IF
         !k=k+1
       END IF !.....
     END DO  !k
   END DO  !j
   !!1 !----Faces_ZX---------
   !!1 IF (nType(2:2)=='w') THEN
   !!1   i=ix0
   !!1 ELSE
   !!1   i=ix1+1
   !!1 END IF
   !!1 DO j=iy0,iy1
   !!1   k=iz0+1;in_out=0
   !!1   DO WHILE(.NOT.ASSOCIATED(Faces_ZX(i,j,k)%Face).AND.in_out/=4.AND.k/=iz1+1)
   !!1     in_out=Vertices(i-1,j, k-1)%in_out &
   !!1           +Vertices(i  ,j, k-1)%in_out &
   !!1           +Vertices(i-1,j, k  )%in_out   &
   !!1           +Vertices(i  ,j, k  )%in_out
   !!1     IF (.NOT.ASSOCIATED(Faces_ZX(i,j,k)%Face).AND.in_out<-2) THEN
   !!1         NrRN_FZX=NrRN_FZX+1
   !!1     END IF
   !!1     k=k+1
   !!1   END DO
   !!1 END DO
   !!1 Floor(ib)%NrRN_FacesXY=Floor(ib)%NrRN_FacesXY+NrRN_FXY
   Floor(ib)%NrRN_FacesYZ=Floor(ib)%NrRN_FacesYZ+NrRN_FYZ
   !!1 Floor(ib)%NrRN_FacesZX=Floor(ib)%NrRN_FacesZX+NrRN_FZX
   !!1 Floor(ib)%NrW_FacesXY=Floor(ib)%NrW_FacesXY+NrRN_FXY
   Floor(ib)%NrW_FacesYZ=Floor(ib)%NrW_FacesYZ+NrRN_FYZ
   !!1 Floor(ib)%NrW_FacesZX=Floor(ib)%NrW_FacesZX+NrRN_FZX
END SUBROUTINE CountFaceSubMontane_WE

SUBROUTINE CountFaceSubMontane_NS(ib)
   !Count Face Border Submontane 'n,s''i and o'
   INTEGER :: ib

   INTEGER :: i,j,k,in_out,blj 
   INTEGER :: NrRN_FXY,NrRN_FYZ,NrRN_FZX
   !!1 !......................................
   !!1 !Für BorderNorthSouth -Face FV (ZX) aktiv
   !!1 !......................................
   !CALL Set(Floor(ib))
   CALL Set_EnvBlkFace(Floor(ib),ib)
   NrRN_FXY=0;NrRN_FYZ=0;NrRN_FZX=0

   !!1 !----Faces_XY---------
   !!1 IF (nType(2:2)=='s') THEN
   !!1     j=iy0
   !!1 ELSE
   !!1     j=iy1+1
   !!1 END IF
   !!1 DO i=ix0+1,ix1
   !!1    k=iz0;in_out=0
   !!1    !DO WHILE(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face).AND.in_out/=4.AND.k/=iz1+1)
   !!1    !unter DurranHatS800.grid werden: 
   !!1    !        1) wenn z=1 und Face%in_out=4 Face_XY dennoch allokiet (GrenzZelle bedingt)
   !!1    !        2) sowohl Face_XY%in_out=(+4,0,-4) jeweils bis z=11 allokiert, an Grenze /Schnitt
   !!1    DO WHILE(k/=iz1+1)
   !!1      in_out=Vertices(i-1,j-1, k  )%in_out &
   !!1            +Vertices(i  ,j-1, k  )%in_out &
   !!1            +Vertices(i-1,j  , k  )%in_out   &
   !!1            +Vertices(i  ,j  , k  )%in_out
   !!1      !IF (.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face).AND.in_out<-2) THEN
   !!1      IF (in_out<-1) THEN
   !!1          NrRN_FXY=NrRN_FXY+1
   !!1      END IF
   !!1    k=k+1
   !!1    END DO
   !!1 END DO
   !!1 !----Faces_YZ---------
   !!1 IF (nType(2:2)=='s') THEN
   !!1     j=iy0
   !!1 ELSE
   !!1     j=iy1+1
   !!1 END IF
   !!1 DO i=ix0,ix1
   !!1    k=iz0+1;in_out=0
   !!1    !DO WHILE((.NOT.ASSOCIATED(Faces_YZ(i,j,k)%Face).AND.in_out/=4 .OR. &
   !!1    !         ASSOCIATED(Faces_YZ(i,j,k)%Face).AND.in_out<=-2) .AND. k/=iz1+1)
   !!1    DO WHILE(k/=iz1+1)
   !!1      in_out=Vertices(i,j-1, k-1)%in_out &
   !!1            +Vertices(i,j-1, k  )%in_out &
   !!1            +Vertices(i,j  , k-1)%in_out   &
   !!1            +Vertices(i,j  , k  )%in_out
   !!1      !IF ((.NOT.ASSOCIATED(Faces_YZ(i,j,k)%Face).AND.in_out<-2) .OR. &
   !!1      !     (ASSOCIATED(Faces_YZ(i,j,k)%Face).AND.in_out<=-2) ) THEN
   !!1      IF (in_out<-1) THEN
   !!1          NrRN_FYZ=NrRN_FYZ+1
   !!1      END IF
   !!1    k=k+1
   !!1    END DO
   !!1 END DO
   !----Faces_ZX---------
   IF (nType(2:2)=='s') THEN
      j=iy0-1
      blj=iy0
   ELSE
      j=iy1+1
      blj=iy1
   END IF
   DO i=ix0+1,ix1
     k=iz0+1
     in_out=Vertices(i-1,j, k-1)%in_out+Vertices(i  ,j, k-1)%in_out &
           +Vertices(i-1,j, k  )%in_out+Vertices(i  ,j, k  )%in_out
     DO k=iz0+1,iz1
     !DO WHILE(.NOT.ASSOCIATED(Faces_ZX(i,j,k)%Face).AND.in_out/=4.AND.k/=iz1+1)
     IF (.NOT.ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
       in_out=Vertices(i-1,j, k-1)%in_out+Vertices(i  ,j, k-1)%in_out &
             +Vertices(i-1,j, k  )%in_out+Vertices(i  ,j, k  )%in_out
      ! IF ((.NOT.ASSOCIATED(Faces_ZX(i,j,k)%Face).AND.in_out<-2) .OR. &
      !      (in_out==-2 .AND. &
      !      (MAX(Vertices(i-1,j, k-1)%in_out,Vertices(i  ,j, k-1)%in_out, &
      !          Vertices(i-1,j, k  )%in_out,Vertices(i  ,j, k  )%in_out)==1).AND. &
      !      (MIN(Vertices(i-1,j, k-1)%in_out,Vertices(i  ,j, k-1)%in_out, &
      !          Vertices(i-1,j, k  )%in_out,Vertices(i  ,j, k  )%in_out)==-1)) ) THEN
      !      !"Valley3D_Blk7u8.grid 'io'  ,Test_Valley3D_V1_Blk
       !IF(Vertices(i-1,blj, k-1)%in_out+Vertices(i  ,blj, k-1)%in_out &
       !   +Vertices(i-1,blj, k  )%in_out+Vertices(i  ,blj, k  )%in_out>=2) THEN
       !       !Write (*,*) " EXIT aus  CountFaceSubMontane_NS"
       !       !Write(*,*) " FaceZX-Block(blj)=maxVol und FaceZX-(j)Rand ist Vol==0.d0"
       !       !Write(*,*) " Randzählung läuft über letztes FaceZX-Block(blj) mit ec=1"
       !       EXIT
       !ELSE
       !  IF (in_out<-2) THEN
         !IF (in_out<=-2) THEN !"Valley3D_Blk7u8.grid 1 zuviel gezählt,Test_Valley3D_V1_Blk
         !IF (in_out<-1) THEN !"Valley3D_Blk7u8.grid 1 zuviel gezählt,Test_Valley3D_V1_Blk
         !IF (in_out<=-1) THEN !"Valley3D_Blk7u8.grid 1 zuviel gezählt,Test_Valley3D_V1_Blk
         !IF (in_out<0) THEN ! wenn Angleichung Weight2_Abfrage wird zu viel gezählt
         !                    !"Valley3D_Blk7u8.grid 1 zuviel gezählt,Test_Valley3D_V1_Blk
         IF (in_out<=0) THEN ! wenn Angleichung Weight2_Abfrage wird zu viel gezählt
                             !"Valley3D_Blk7u8.grid 1 zuviel gezählt,Test_Valley3D_V1_Blk
           NrRN_FZX=NrRN_FZX+1
         END IF
         !k=k+1
       !END IF
     END IF
     END DO
   END DO
   !!1 Floor(ib)%NrRN_FacesXY=Floor(ib)%NrRN_FacesXY+NrRN_FXY
   !!1 Floor(ib)%NrRN_FacesYZ=Floor(ib)%NrRN_FacesYZ+NrRN_FYZ
   Floor(ib)%NrRN_FacesZX=Floor(ib)%NrRN_FacesZX+NrRN_FZX
   !!1 Floor(ib)%NrW_FacesXY=Floor(ib)%NrW_FacesXY+NrRN_FXY
   !!1 Floor(ib)%NrW_FacesYZ=Floor(ib)%NrW_FacesYZ+NrRN_FYZ
   Floor(ib)%NrW_FacesZX=Floor(ib)%NrW_FacesZX+NrRN_FZX
END SUBROUTINE CountFaceSubMontane_NS

SUBROUTINE CountFaceSubMontane_TB(ib)
   !Count Face Border Submopntane 'it,ib,ob'
   INTEGER :: ib

   INTEGER :: i,j,k,in_out 
   INTEGER :: NrRN_FXY,NrRN_FYZ,NrRN_FZX
   !!1 !......................................
   !!1 !Für BorderTopBottom -Face FW (XY) aktiv
   !!1 !......................................

   !CALL Set(Floor(ib))
   CALL Set_EnvBlkFace(Floor(ib),ib)
   NrRN_FXY=0; NrRN_FYZ=0; NrRN_FZX=0
   k=-9999

   !----Faces_XY---------
   IF (nType(2:2)=='t') THEN
      k=iz1+1
   ELSE IF (nType=='ib') THEN
      k=iz0-1  ! 'ob' extra, da grenze z0 verschoben sein kann unter iz0
   END IF
   IF (k/=-9999) THEN
      DO i=ix0+1,ix1
         DO j=iy0+1,iy1
          IF(.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
             in_out=Vertices(i-1,j-1,k)%in_out+Vertices(i  ,j-1,k)%in_out &
                    +Vertices(i-1,j  ,k)%in_out+Vertices(i  ,j  ,k)%in_out
            !IF (in_out<-2 .OR. &
            !    (.NOT.ASSOCIATED(Faces_XY(i,j,k)%Face).AND.in_out<-2).OR. &
            !    (in_out==-2 .AND. &
            !     (MAX(Vertices(i-1,j-1,k)%in_out+Vertices(i  ,j-1,k)%in_out, &
            !          Vertices(i-1,j  ,k)%in_out+Vertices(i  ,j  ,k)%in_out)==1).AND. &
            !     (MIN(Vertices(i-1,j-1,k)%in_out+Vertices(i  ,j-1,k)%in_out, &
            !          Vertices(i-1,j  ,k)%in_out+Vertices(i  ,j  ,k)%in_out)==-1)) ) THEN
             IF (in_out<=0) THEN
                 NrRN_FXY=NrRN_FXY+1
             END IF
          END IF
        END DO
      END DO
      !!1 !----Faces_YZ---------
      !!1 IF (nType=='ib') THEN
      !!1   k=iz0
      !!1 ELSE IF (nType=='it') THEN
      !!1   k=iz1+1
      !!1 END IF
      !!1 DO i=ix0,ix1
      !!1   DO j=iy0+1,iy1
      !!1     in_out=Vertices(i,j-1, k-1)%in_out &
      !!1           +Vertices(i,j-1, k  )%in_out &
      !!1           +Vertices(i,j  , k-1)%in_out &
      !!1           +Vertices(i,j  , k  )%in_out
      !!1     IF (.NOT.ASSOCIATED(Faces_YZ(i,j,k)%Face).AND.in_out<-2) THEN
      !!1         NrRN_FYZ=NrRN_FYZ+1
      !!1     END IF
      !!1   END DO
      !!1 END DO
      !!1 !----Faces_ZX---------
      !!1 IF (nType=='ib') THEN
      !!1   k=iz0
      !!1 ELSE IF (nType=='it') THEN
      !!1   k=iz1+1
      !!1 END IF
      !!1 DO i=ix0+1,ix1
      !!1   DO j=iy0,iy1
      !!1     in_out=Vertices(i-1,j, k-1)%in_out &
      !!1           +Vertices(i  ,j, k-1)%in_out &
      !!1           +Vertices(i-1,j, k  )%in_out &
      !!1           +Vertices(i  ,j, k  )%in_out
      !!1     IF (.NOT.ASSOCIATED(Faces_ZX(i,j,k)%Face).AND.in_out<-2) THEN
      !!1         NrRN_FZX=NrRN_FZX+1
      !!1     END IF
      !!1   END DO
      !!1 END DO
      Floor(ib)%NrRN_FacesXY=Floor(ib)%NrRN_FacesXY+NrRN_FXY
      !!1 Floor(ib)%NrRN_FacesYZ=Floor(ib)%NrRN_FacesYZ+NrRN_FYZ
      !!1 Floor(ib)%NrRN_FacesZX=Floor(ib)%NrRN_FacesZX+NrRN_FZX
      Floor(ib)%NrW_FacesXY=Floor(ib)%NrW_FacesXY+NrRN_FXY
      !!1 Floor(ib)%NrW_FacesYZ=Floor(ib)%NrW_FacesYZ+NrRN_FYZ
      !!1 Floor(ib)%NrW_FacesZX=Floor(ib)%NrW_FacesZX+NrRN_FZX
   END IF  ! "it/ot/ib"
END SUBROUTINE CountFaceSubMontane_TB


FUNCTION VolCell(Cell)
   REAL(8) :: VolCell
   TYPE (Cell_T) :: Cell
   REAL(8) :: lam1,lam2,r1,r2,phi1,phi2,z2,z1
   !phi-> [-Pi/2,Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. & 
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
       indata_type=='geo' .OR. indata_type=='rad' .OR. &
       indata_type=='wgs84') THEN
     phi1=Cell%Face1%Edge2%Vert1%Point%y
     phi2=Cell%Face1%Edge2%Vert2%Point%y
     lam1=Cell%Face1%Edge1%Vert1%Point%x
     lam2=Cell%Face1%Edge1%Vert2%Point%x
     r1  =RadEarth+Cell%Face1%Edge1%Vert1%Point%z
     r2  =RadEarth+Cell%Face2%Edge1%Vert1%Point%z
     VolCell=1.0d0/3.0d0*((lam2-lam1)*(r2**3-r1**3)*(SIN(phi2)-SIN(phi1)))
   ELSE If (indata_type=='cyl') THEN
     lam1=Cell%Face1%Edge1%Vert1%Point%x
     lam2=Cell%Face1%Edge1%Vert2%Point%x
     r1=Cell%Face1%Edge2%Vert1%Point%y
     r2=Cell%Face1%Edge2%Vert2%Point%y
     z1=Cell%Face1%Edge1%Vert1%Point%z
     z2=Cell%Face2%Edge1%Vert1%Point%z
     VolCell=0.5d0*ABS((lam2-lam1) &
                      *(r2**2-r1**2) &
                      *(z2-z1))
   ELSE !('Cart',DEFAULT)
       VolCell=ABS((Cell%Face1%Edge1%Vert2%Point%x-Cell%Face1%Edge1%Vert1%Point%x) &
                  *(Cell%Face1%Edge2%Vert2%Point%y-Cell%Face1%Edge2%Vert1%Point%y) &
                  *(Cell%Face2%Edge1%Vert1%Point%z-Cell%Face1%Edge1%Vert1%Point%z))
   END IF
END FUNCTION VolCell

FUNCTION VolCellP8(Cell)
   REAL(8) :: VolCellP8
   TYPE(CellP8_T) :: Cell
   REAL(8) :: lam1,lam2,r1,r2,phi1,phi2,z2,z1
   !phi-> [-Pi/2,Pi/2]; lam ->[0,2Pi]

   !('Globe')
   IF((indata_type=='gk'.AND.conv_ingrid=='spherical') .OR. & 
      (indata_type=='gk'.AND.out_wahlgrid=='G') .OR. &
       indata_type=='geo' .OR. indata_type=='rad' .OR. &
       indata_type=='wgs84') THEN
     phi1=Cell%V1%Point%y  !Cell%Face1%Edge2%Vert1%Point%y
     phi2=Cell%V3%Point%y  !Cell%Face1%Edge2%Vert2%Point%y
     lam1=Cell%V0%Point%x  !Cell%Face1%Edge1%Vert1%Point%x
     lam2=Cell%V1%Point%x  !Cell%Face1%Edge1%Vert2%Point%x
     r1  =RadEarth+Cell%V0%Point%z  !Cell%Face1%Edge1%Vert1%Point%z
     r2  =RadEarth+Cell%V4%Point%z  !Cell%Face2%Edge1%Vert1%Point%z
     VolCellP8=1.0d0/3.0d0*((lam2-lam1)*(r2**3-r1**3)*(SIN(phi2)-SIN(phi1)))
   ELSE If (indata_type=='cyl') THEN
     lam1=Cell%V0%Point%x  !Cell%Face1%Edge1%Vert1%Point%x
     lam2=Cell%V1%Point%x  !Cell%Face1%Edge1%Vert2%Point%x
     r1=Cell%V1%Point%y    !Cell%Face1%Edge2%Vert1%Point%y
     r2=Cell%V3%Point%y    !Cell%Face1%Edge2%Vert2%Point%y
     z1=Cell%V0%Point%z    !Cell%Face1%Edge1%Vert1%Point%z
     z2=Cell%V4%Point%z    !Cell%Face2%Edge1%Vert1%Point%z
     VolCellP8=0.5d0*ABS((lam2-lam1) &
                      *(r2**2-r1**2) &
                      *(z2-z1))
   ELSE !('Cart',DEFAULT)
       VolCellP8=ABS((Cell%V1%Point%x-Cell%V0%Point%x) &
                    *(Cell%V3%Point%y-Cell%V1%Point%y) &
                    *(Cell%V4%Point%z-Cell%V0%Point%z))
   END IF
END FUNCTION VolCellP8

SUBROUTINE InitAllCells
  INTEGER :: ib,i,j,k

  DO ib=1,nb
    CALL Set(Floor(ib))

    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_XY(i,j,k-1)%Face) .AND. ASSOCIATED(Cell(i,j,k)%Cell) ) THEN
            Cell(i,j,k)%Cell%Face1=>Faces_XY(i,j,k-1)%Face
          END IF
          IF (ASSOCIATED(Faces_XY(i,j,k)%Face) .AND. ASSOCIATED(Cell(i,j,k)%Cell) ) THEN
            Cell(i,j,k)%Cell%Face2=>Faces_XY(i,j,k)%Face
          END IF
          IF (ASSOCIATED(Faces_ZX(i,j-1,k)%Face) .AND. ASSOCIATED(Cell(i,j,k)%Cell) ) THEN
            Cell(i,j,k)%Cell%Face3=>Faces_ZX(i,j-1,k)%Face
          END IF
          IF (ASSOCIATED(Faces_ZX(i,j,k)%Face) .AND. ASSOCIATED(Cell(i,j,k)%Cell) ) THEN
            Cell(i,j,k)%Cell%Face4=>Faces_ZX(i,j,k)%Face
          END IF
          IF (ASSOCIATED(Faces_YZ(i-1,j,k)%Face) .AND. ASSOCIATED(Cell(i,j,k)%Cell) ) THEN
            Cell(i,j,k)%Cell%Face5=>Faces_YZ(i-1,j,k)%Face
          END IF
          IF (ASSOCIATED(Faces_YZ(i,j,k)%Face) .AND. ASSOCIATED(Cell(i,j,k)%Cell) ) THEN
            Cell(i,j,k)%Cell%Face6=>Faces_YZ(i,j,k)%Face
          END IF
!          SELECT CASE (OutGrid)
!            CASE ('Globe')
!              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!                Cell(i,j,k)%Cell%Vol=VolCell(Cell(i,j,k)%Cell)
!              END IF
!            CASE ('Cylinder')
!              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!                Cell(i,j,k)%Cell%Vol=VolCell(Cell(i,j,k)%Cell)
!              END IF
!            CASE ('Cart')
!              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!                Cell(i,j,k)%Cell%Vol=dx(i)*dy(j)*dz(k)
!              END IF
!            CASE DEFAULT
!              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!                Cell(i,j,k)%Cell%Vol=dx(i)*dy(j)*dz(k)
!              END IF
!          END SELECT
        END DO
      END DO
    END DO
  END DO  ! ib
END SUBROUTINE InitAllCells


SUBROUTINE SearchPointsEdgeZaxis(P0,P1,Vertex1,Vertex2,zc0,zc1)
  TYPE (Vertex_T) :: Vertex1,Vertex2
  REAL(8) :: zc0,zc1
  TYPE(Point_T)   :: P0,P1
  REAL(8) :: zp,dist1,dist2
  TYPE (Vertex_T) :: InterVert

        dist1=sf(nr_wahlfkt,Vertex1%Point%x,Vertex1%Point%y,Vertex1%Point%z)
        dist2=sf(nr_wahlfkt,Vertex2%Point%x,Vertex2%Point%y,Vertex2%Point%z)
        IF (dist2>=0.0d0) THEN
          zp=zc0
          IF (dist1<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            zp=Intervert%Point%z
          END IF
          P0%x=Vertex1%Point%x
          P0%y=Vertex1%Point%y
          P0%z=zp
          P1%x=Vertex2%Point%x
          P1%y=Vertex2%Point%y
          P1%z=zc1
        ELSE IF (dist1>=0.0d0) THEN
          zp=zc1
          IF (dist2<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            zp=Intervert%Point%z
          END IF
          P0%x=Vertex1%Point%x
          P0%y=Vertex1%Point%y
          P0%z=zc0
          P1%x=Vertex2%Point%x
          P1%y=Vertex2%Point%y
          P1%z=zp
        ELSE
          !Fall,wenn letzter Streifen Vol-Berechnung
          !beide unter Null-Z-Ebene
          IF (dist1<dist2) THEN
            P0%x=Vertex2%Point%x
            P0%y=Vertex2%Point%y
            P0%z=Vertex2%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          ELSE IF (dist2<dist1) THEN
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex1%Point%x
            P1%y=Vertex1%Point%y
            P1%z=Vertex1%Point%z
          ELSE ! dist1==dist2
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          END IF
        END IF
END SUBROUTINE SearchPointsEdgeZaxis


SUBROUTINE SearchPointsEdgeYaxis(P0,P1,Vertex1,Vertex2,yc0,yc1)
  TYPE (Vertex_T) :: Vertex1,Vertex2
  REAL(8) :: yc0,yc1
  TYPE(Point_T)   :: P0,P1
  REAL(8) :: yp,dist1,dist2
  TYPE (Vertex_T) :: InterVert

        dist1=sf(nr_wahlfkt,Vertex1%Point%x,Vertex1%Point%y,Vertex1%Point%z)
        dist2=sf(nr_wahlfkt,Vertex2%Point%x,Vertex2%Point%y,Vertex2%Point%z)
        IF (dist2>=0.0d0) THEN
          yp=yc0
          IF (dist1<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            yp=Intervert%Point%y
          END IF
          P0%x=Vertex1%Point%x
          P0%y=yp
          P0%z=Vertex1%Point%z
          P1%x=Vertex2%Point%x
          P1%y=yc1
          P1%z=Vertex2%Point%z
        ELSE IF (dist1>=0.0d0) THEN
          yp=yc1
          IF (dist2<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            yp=Intervert%Point%y
          END IF
          P0%x=Vertex1%Point%x
          P0%y=yc0
          P0%z=Vertex1%Point%z
          P1%x=Vertex2%Point%x
          P1%y=yp
          P1%z=Vertex2%Point%z
        ELSE
          !Fall,wenn letzter Streifen Vol-Berechnung
          !beide unter Null-Z-Ebene
          IF (dist1<dist2) THEN
            P0%x=Vertex2%Point%x
            P0%y=Vertex2%Point%y
            P0%z=Vertex2%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          ELSE IF (dist2<dist1) THEN
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex1%Point%x
            P1%y=Vertex1%Point%y
            P1%z=Vertex1%Point%z
          ELSE !dist1==dist2
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          END IF
        END IF
END SUBROUTINE SearchPointsEdgeYaxis


SUBROUTINE SearchPointsEdgeXaxis(P0,P1,Vertex1,Vertex2,xc0,xc1)
  TYPE (Vertex_T) :: Vertex1,Vertex2
  REAL(8) :: xc0,xc1
  TYPE(Point_T)   :: P0,P1
  REAL(8) :: xp,dist1,dist2
  TYPE (Vertex_T) :: InterVert

        dist1=sf(nr_wahlfkt,Vertex1%Point%x,Vertex1%Point%y,Vertex1%Point%z)
        dist2=sf(nr_wahlfkt,Vertex2%Point%x,Vertex2%Point%y,Vertex2%Point%z)
        IF (dist2>=0.0d0) THEN
          xp=xc0
          IF (dist1<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            xp=Intervert%Point%x
          END IF
          P0%x=xp
          P0%y=Vertex1%Point%y
          P0%z=Vertex1%Point%z
          P1%x=xc1
          P1%y=Vertex2%Point%y
          P1%z=Vertex2%Point%z
        ELSE IF (dist1>=0.0d0) THEN
          xp=xc1
          IF (dist2<0.0d0) THEN
            CALL SearchPoint(InterVert,Vertex1,Vertex2)
            xp=Intervert%Point%x
          END IF
          P0%x=xc0
          P0%y=Vertex1%Point%y
          P0%z=Vertex1%Point%z
          P1%x=xp
          P1%y=Vertex2%Point%y
          P1%z=Vertex2%Point%z
        ELSE
          !Fall,wenn letzter Streifen Vol-Berechnung
          !beide unter Null-Z-Ebene
          IF (dist1<dist2) THEN
            P0%x=Vertex2%Point%x
            P0%y=Vertex2%Point%y
            P0%z=Vertex2%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          ELSE IF (dist2<dist1) THEN
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex1%Point%x
            P1%y=Vertex1%Point%y
            P1%z=Vertex1%Point%z
          ELSE
            P0%x=Vertex1%Point%x
            P0%y=Vertex1%Point%y
            P0%z=Vertex1%Point%z
            P1%x=Vertex2%Point%x
            P1%y=Vertex2%Point%y
            P1%z=Vertex2%Point%z
          END IF
        END IF
END SUBROUTINE SearchPointsEdgeXaxis


SUBROUTINE SearchMinMaxEdge(Edge,PMin,PMax)
  TYPE(Edge_T) :: Edge
  TYPE(Point_T) :: PMin,PMax

  IF (Edge%Vert1%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Edge%Vert1%Point%x)
    PMin%y=MIN(PMin%y,Edge%Vert1%Point%y)
    PMin%z=MIN(PMin%z,Edge%Vert1%Point%z)
    PMax%x=MAX(PMax%x,Edge%Vert1%Point%x)
    PMax%y=MAX(PMax%y,Edge%Vert1%Point%y)
    PMax%z=MAX(PMax%z,Edge%Vert1%Point%z)
  END IF
  IF (Edge%Vert2%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Edge%Vert2%Point%x)
    PMin%y=MIN(PMin%y,Edge%Vert2%Point%y)
    PMin%z=MIN(PMin%z,Edge%Vert2%Point%z)
    PMax%x=MAX(PMax%x,Edge%Vert2%Point%x)
    PMax%y=MAX(PMax%y,Edge%Vert2%Point%y)
    PMax%z=MAX(PMax%z,Edge%Vert2%Point%z)
  END IF
  IF (ASSOCIATED(Edge%VertS)) THEN
    PMin%x=MIN(PMin%x,Edge%VertS%Point%x)
    PMin%y=MIN(PMin%y,Edge%VertS%Point%y)
    PMin%z=MIN(PMin%z,Edge%VertS%Point%z)
    PMax%x=MAX(PMax%x,Edge%VertS%Point%x)
    PMax%y=MAX(PMax%y,Edge%VertS%Point%y)
    PMax%z=MAX(PMax%z,Edge%VertS%Point%z)
  END IF

END SUBROUTINE SearchMinMaxEdge

SUBROUTINE SearchMinMaxFace(Face,PMin,PMax)
  TYPE(Face_T) :: Face
  TYPE(Point_T) :: PMin,PMax
  CALL SearchMinMaxEdge(Face%Edge1,PMin,PMax)
  CALL SearchMinMaxEdge(Face%Edge2,PMin,PMax)
  CALL SearchMinMaxEdge(Face%Edge3,PMin,PMax)
  CALL SearchMinMaxEdge(Face%Edge4,PMin,PMax)

END SUBROUTINE SearchMinMaxFace

SUBROUTINE SearchMinMaxCell(Cell,PMin,PMax)
  TYPE(Cell_T) :: Cell
  TYPE(Point_T) :: PMin,PMax
  CALL SearchMinMaxFace(Cell%Face1,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face2,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face3,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face4,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face5,PMin,PMax)
  CALL SearchMinMaxFace(Cell%Face6,PMin,PMax)

END SUBROUTINE SearchMinMaxCell


SUBROUTINE SearchMinMaxXcFace(Face,xvc0,xvc1)
  TYPE(Face_T) ::Face
  REAL(8) :: xvc0,xvc1
  REAL(8) :: xc0,xc1,xvc0u,xvc0o,xvc1u,xvc1o

  xc0=Face%Edge1%Vert1%Point%x
  xc1=Face%Edge2%Vert2%Point%x
  IF (Face%Edge1%yes_sp==1 ) THEN
    IF(Face%Edge1%Vert1%in_out==1) THEN
       xvc0u=xc0
       xvc1u=Face%Edge1%VertS%Point%x
    ELSE
       xvc0u=Face%Edge1%VertS%Point%x
       xvc1u=xc1
    END IF
  ELSE
    xvc0u=xc0
    xvc1u=xc1
  END IF
  IF (Face%Edge3%yes_sp==1 ) THEN
    IF(Face%Edge3%Vert1%in_out==1) THEN
      xvc0o=xc0
      xvc1o=Face%Edge3%VertS%Point%x
    ELSE
      xvc0o=Face%Edge3%VertS%Point%x
      xvc1o=xc1
    END IF
  ELSE
    xvc0o=xc0
    xvc1o=xc1
  END IF
  xvc0=MIN(xvc0u,xvc0o)
  xvc1=MAX(xvc1u,xvc1o)

END SUBROUTINE SearchMinMaxXcFace


SUBROUTINE SearchMinMaxYcFace(Face,yvc0,yvc1)
  TYPE(Face_T) ::Face
  REAL(8) :: yvc0,yvc1
  REAL(8) :: yc0,yc1,yvc0u,yvc0o,yvc1u,yvc1o

  yc0=Face%Edge1%Vert1%Point%y
  yc1=Face%Edge2%Vert2%Point%y
  IF (Face%Edge2%yes_sp==1 ) THEN
    IF(Face%Edge2%Vert1%in_out==1) THEN
       yvc0u=yc0
       yvc1u=Face%Edge2%VertS%Point%y
    ELSE
       yvc0u=Face%Edge2%VertS%Point%y
       yvc1u=yc1
    END IF
  ELSE
    yvc0u=yc0
    yvc1u=yc1
  END IF
  IF (Face%Edge4%yes_sp==1 ) THEN
    IF(Face%Edge4%Vert1%in_out==1) THEN
      yvc0o=yc0
      yvc1o=Face%Edge4%VertS%Point%y
    ELSE
      yvc0o=Face%Edge4%VertS%Point%y
      yvc1o=yc1
    END IF
  ELSE
    yvc0o=yc0
    yvc1o=yc1
  END IF
  yvc0=MIN(yvc0u,yvc0o)
  yvc1=MAX(yvc1u,yvc1o)

END SUBROUTINE SearchMinMaxYcFace

SUBROUTINE SetLandClass(Cell,ci,cj,ck,ib)
  TYPE(Cell_T) :: Cell
  INTEGER :: ci,cj,ck,ib

  INTEGER :: ix,iy
  INTEGER :: nl
  !LandDef bezogen auf Grundgitter,
  !a) SetLandClassCell_CLC(Cell,ci,cj,ck,ib)
  !b) SetLandClass_allg
  !b1.abtasten welche Teile LandDef-Bereich mit Floor(ib)%[igx0..igx1]
  !  oder Floor(ib)%[igy0..igy1] überlappen
  !b2.abtasten spez. Cell-Index in LandDef-Bereich passt 
  !b3.Zuweisung LandClass-Def.
  IF(ifcorine) THEN
    ix=ci*2.e0**RefineX
    iy=cj*2.e0**RefineY
    !Cell%LandClass=LandCover(ix,iy)
    !Cell%LandClass=v_clc(LandCover(ix,iy))%clc_c
    Cell%LandClass=v_clc(LandCover(ix,iy))%mat_c
  ELSE
    IF(ASSOCIATED(LandDef)) THEN
      DO nl=1,nr_landdef
        IF((LandDef(1,1,nl)>igx0.AND.LandDef(1,1,nl)<=igx1).OR. &
           (LandDef(2,1,nl)>igx0.AND.LandDef(2,1,nl)<=igx1).OR. &
           (LandDef(1,2,nl)>igy0.AND.LandDef(1,2,nl)<=igy1).OR. &
           (LandDef(2,2,nl)>igy0.AND.LandDef(2,2,nl)<=igy1)) THEN
           IF(ci>LandDef(1,1,nl)*2.e0**RefineX .AND. &
              ci<=LandDef(2,1,nl)*2.e0**RefineX  .AND. &
              cj>LandDef(1,2,nl)*2.e0**RefineY .AND. &
              cj<=LandDef(2,2,nl)*2.e0**RefineY ) THEN
                 Cell%LandClass=LandDef(1,3,nl)
                 EXIT
           END IF
         END IF
      END DO
    END IF
  END IF
END SUBROUTINE SetLandClass

SUBROUTINE SetLandClassCell_CLC(Cell,ci,cj,ck,ib)
  TYPE(Cell_T) :: Cell
  INTEGER :: ci,cj,ck,ib
  INTEGER :: ix,iy
  IF(ifcorine) THEN
    ix=ci*2.e0**RefineX
    iy=cj*2.e0**RefineY
    !Cell%LandClass=LandCover(ix,iy)
    !Cell%LandClass=v_clc(LandCover(ix,iy))%clc_c
    Cell%LandClass=v_clc(LandCover(ix,iy))%mat_c
  END IF

END SUBROUTINE SetLandClassCell_CLC


SUBROUTINE SetLandCoverCellsCLC
  INTEGER :: i,j,k
  INTEGER :: ix,iy

  IF(ifcorine) THEN
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
             !..................................
             ix=i*2.e0**RefineX
             iy=j*2.e0**RefineY
             !Cell(i,j,k)%Cell%LandClass=LandCover(ix,iy)
             !Datei CLC.csv
             !Cell(i,j,k)%Cell%LandClass=v_clc(LandCover(ix,iy))%clc_c
             !Datei CLC_IfT.csv
             Cell(i,j,k)%Cell%LandClass=v_clc(LandCover(ix,iy))%mat_c
          END IF                              
        END DO  ! k
      END DO  ! j
    END DO  ! i
    nr_landdef=Nr_Land_Cl
  END IF

END SUBROUTINE SetLandCoverCellsCLC


SUBROUTINE AnalyzeCell(Cell)
  TYPE(Cell_T) :: Cell
  TYPE(Cell_T), POINTER :: pCell=>NULL()
  INTEGER :: ListCut(1:2,1:8)
  INTEGER :: nCut,nCutNeu,i,j,iv
  INTEGER :: ic,jc,kc,nxc,nyc,nzc
  INTEGER :: iTemp
  REAL(8) :: xc0,xc1,yc0,yc1,zc0,zc1,dx,dy,dz
  TYPE (Vertex_T) :: Vertex1,Vertex2
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(Point_T) :: EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7
  TYPE(Point_T) :: PMin,PMax,SumMidP,MidPoint
  TYPE(Point_T) :: PCut(1:8)
  REAL(8) :: xvc0,xvc1,xvu0,xvu1,xvo0,xvo1
  REAL(8) :: yvc0,yvc1,yvu0,yvu1,yvo0,yvo1
  REAL(8) :: zvc0,zvc1
  REAL(8) :: Vol,Len
  INTEGER :: VertexFace(1:6,1:9)
  LOGICAL :: Cut
  INTEGER :: iP1,iP2,jP1,jP2,iVert,jVert

  !Ossi
  EdgP0=Cell%Face1%Edge1%Vert1%Point
  EdgP1=Cell%Face1%Edge1%Vert2%Point
  EdgP2=Cell%Face1%Edge3%Vert1%Point
  EdgP3=Cell%Face1%Edge3%Vert2%Point
  EdgP4=Cell%Face2%Edge1%Vert1%Point
  EdgP5=Cell%Face2%Edge1%Vert2%Point
  EdgP6=Cell%Face2%Edge3%Vert1%Point
  EdgP7=Cell%Face2%Edge3%Vert2%Point
  P0=PointParametricEarth(EdgP0)
  P1=PointParametricEarth(EdgP1)
  P2=PointParametricEarth(EdgP2)
  P3=PointParametricEarth(EdgP3)
  P4=PointParametricEarth(EdgP4)
  P5=PointParametricEarth(EdgP5)
  P6=PointParametricEarth(EdgP6)
  P7=PointParametricEarth(EdgP7)
  Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
  Cell%MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
  !----------------------------------------------------------------------
  VertexFace=-1
  nCut=0
  ListCut=0
  Cell%in_out=Cell%Face1%in_out+Cell%Face2%in_out
  
  IF (Cell%Face1%NumberVert>2) THEN
    VertexFace(1,1:Cell%Face1%NumberVert)=Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
    VertexFace(1,Cell%Face1%NumberVert+1)=VertexFace(1,1)
  END IF
  IF (Cell%Face2%NumberVert>2) THEN
    VertexFace(2,1:Cell%Face2%NumberVert)=Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
    VertexFace(2,Cell%Face2%NumberVert+1)=VertexFace(2,1)
  END IF
  IF (Cell%Face3%NumberVert>2) THEN
    VertexFace(3,1:Cell%Face3%NumberVert)=Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
    VertexFace(3,Cell%Face3%NumberVert+1)=VertexFace(3,1)
  END IF
  IF (Cell%Face4%NumberVert>2) THEN
    VertexFace(4,1:Cell%Face4%NumberVert)=Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
    VertexFace(4,Cell%Face4%NumberVert+1)=VertexFace(4,1)
  END IF
  IF (Cell%Face5%NumberVert>2) THEN
    VertexFace(5,1:Cell%Face5%NumberVert)=Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
    VertexFace(5,Cell%Face5%NumberVert+1)=VertexFace(5,1)
  END IF
  IF (Cell%Face6%NumberVert>2) THEN
    VertexFace(6,1:Cell%Face6%NumberVert)=Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
    VertexFace(6,Cell%Face6%NumberVert+1)=VertexFace(6,1)
  END IF

  C1:DO i=1,6
    C2:DO iVert=1,8
      iP1=VertexFace(i,iVert)
      iP2=VertexFace(i,iVert+1)
      IF (ip2==-1) THEN
        EXIT C2
      END IF
      Cut=.TRUE.
      C3:DO j=1,6
        IF (j/=i) THEN
          C4:DO jVert=1,8
            jP1=VertexFace(j,jVert)
            jP2=VertexFace(j,jVert+1)
            IF (jp2==-1) THEN
              EXIT C4
            END IF
            IF ((iP1==jP1.AND.iP2==jP2).OR.  &
                (iP1==jP2.AND.iP2==jP1)) THEN
              Cut=.FALSE.
              EXIT C3
            END IF
          END DO C4
        END IF
      END DO C3
      IF (Cut) THEN
        ncut=ncut+1
        ListCut(1,nCut)=iP1
        ListCut(2,nCut)=iP2
      END IF
    END DO C2
  END DO C1
!  IF (Cell%Face1%ec==1) THEN
!    nCut=nCut+1
!    ListCut(1:2,nCut)=Cell%Face1%EdgeCut(1:2)
!  END IF
!  IF (Cell%Face2%ec==1) THEN
!    nCut=nCut+1
!    ListCut(2,nCut)=Cell%Face2%EdgeCut(1)
!    ListCut(1,nCut)=Cell%Face2%EdgeCut(2)
!  END IF
!  IF (Cell%Face3%ec==1) THEN
!    nCut=nCut+1
!    ListCut(1:2,nCut)=Cell%Face3%EdgeCut(1:2)
!  END IF
!  IF (Cell%Face4%ec==1) THEN
!    nCut=nCut+1
!    ListCut(1,nCut)=Cell%Face4%EdgeCut(2)
!    ListCut(2,nCut)=Cell%Face4%EdgeCut(1)
!  END IF
!  IF (Cell%Face5%ec==1) THEN
!    nCut=nCut+1
!    ListCut(1:2,nCut)=Cell%Face5%EdgeCut(1:2)
!  END IF
!  IF (Cell%Face6%ec==1) THEN
!    nCut=nCut+1
!    ListCut(1,nCut)=Cell%Face6%EdgeCut(2)
!    ListCut(2,nCut)=Cell%Face6%EdgeCut(1)
!  END IF

  IF (nCut>0) THEN
    nCutNeu=nCut
    DO i=1,ncut-1
      DO j=i+1,ncut
        IF (( (ListCut(1,i)==ListCut(1,j).AND.ListCut(2,i)==ListCut(2,j)).OR. &
              (ListCut(1,i)==ListCut(2,j).AND.ListCut(2,i)==ListCut(1,j)) ) &
           .AND.ListCut(1,i)>0) THEN
          ListCut(:,j)=0
          nCutNeu=nCutNeu-1
!         iTemP=ListCut(1,j)
!         ListCut(1,j)=ListCut(2,j)
!         ListCut(2,j)=iTemp
        END IF
      END DO
    END DO
    Cell%vc=nCutNeu
    ALLOCATE(Cell%VertCut(nCut))
    ALLOCATE(Cell%VertCutOnlyP(nCut))
    Cell%VertCut(1)= ListCut(1,1)
    Cell%VertCut(2)= ListCut(2,1)
    ListCut(1,1)=0
    ListCut(2,1)=0
    iv=2
    S2:DO
      S1:DO i=1,nCut
        IF (Cell%VertCut(iv)==ListCut(1,i)) THEN
          iv=iv+1
          Cell%VertCut(iv)=ListCut(2,i)
          ListCut(:,i)=0
          EXIT S1
        ELSE IF (Cell%VertCut(iv)==ListCut(2,i)) THEN
          iv=iv+1
          Cell%VertCut(iv)=ListCut(1,i)
          ListCut(:,i)=0
          EXIT S1
        END IF
      END DO S1
      IF (iv>=nCutNeu) THEN
        EXIT S2
      END IF
    END DO S2
  ELSE
    Cell%vc=0
    IF (Cell%in_out==-8) THEN
      Cell%Vol=0.0d0
    END IF
  END IF

  !Abfrage auf Vol-Berechnen umstellen
  ! - vc>0 kann auch (Grenzfläche-komplett) oder (3*in_out=0) sein, 
  !   dann max.Vol überschrieben, durch Vol-Analyze-Vektor-Berechnung
  ! - für Output Weight2 und *out.gmvG besser
  !   wenn Celle%Face1 Grenzfläche ist auf vc>0 setzen
  !   wenn Celle%Face2 Grenzfläche ist auf vc=0 setzen
  !   damit weniger Special-abfragen,
  !   z.B. wenn k==1 und Face1 Grenzfläche ist Ausgabe gewünscht
  ! z.Zt. Celle%Face1 Grenzfläche ist vc=0
  !       Celle%Face2 Grenzfläche ist vc>0 

  ! Volumen-Berechnung
  IF (Cell%vc>0) THEN
    Cell%Vol=0.0d0
    xc0=Cell%Face1%Edge1%Vert1%Point%x
    yc0=Cell%Face1%Edge1%Vert1%Point%y
    zc0=Cell%Face1%Edge1%Vert1%Point%z
    xc1=Cell%Face2%Edge2%Vert2%Point%x
    yc1=Cell%Face2%Edge2%Vert2%Point%y
    zc1=Cell%Face2%Edge2%Vert2%Point%z
    PMin=Cell%Face2%Edge3%Vert2%Point
    PMax=Cell%Face1%Edge1%Vert1%Point
    CALL SearchMinMaxCell(Cell,PMin,PMax)

    IF((Cell%face5%in_out==4.AND.Cell%face6%in_out==-4).OR. &
       (Cell%face5%in_out==-4.AND.Cell%face6%in_out==4)) THEN
       !...Plane x-direction.................................
       zvc0=PMin%z+Shrink*ABS(PMin%z)
       zvc1=PMax%z-Shrink*ABS(PMax%z)
       yvc0=PMin%y+Shrink*ABS(PMin%y)
       yvc1=PMax%y-Shrink*ABS(PMax%y)
       yc0=yvc0
       !...Vergleich zu InitCellVol mit nz=1
       !...IncrVol +10 (1Stelle) nach Komma gleich, +50 (2St.), +150 (3St.),+200 (4St.)
       nxc=IncrVol ! +50 , +200
       nyc=IncrVol ! +50 , +200
       nzc=IncrVol ! +50 , +200
       dy=(yvc1-yvc0)/nyc
       dz=(zvc1-zvc0)/nzc
       Vertex1%Point%x=xc0
       Vertex1%in_out=-1
       Vertex2%Point%x=xc1
       Vertex2%in_out=1
       SumMidP%x=0.0d0
       SumMidP%y=0.0d0
       SumMidP%z=0.0d0
       DO kc=1,nzc
         yvc0=yc0
         DO jc=1,nyc
           !      2--------------6    !Point-Folge Vol-Berechnung
           !     /|             /|
           !    0--------------4 |
           !    | 3------------|-7
           !    |/             |/
           !    1--------------5
           Vertex1%Point%z=zvc0
           Vertex1%Point%y=yvc0
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP1,EdgP5,Vertex1,Vertex2,xc0,xc1)
           P1=PointParametricEarth(EdgP1)
           P5=PointParametricEarth(EdgP5)
   
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex1%Point%y=yvc0
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP0,EdgP4,Vertex1,Vertex2,xc0,xc1)
           P0=PointParametricEarth(EdgP0)
           P4=PointParametricEarth(EdgP4)
   
           Vertex1%Point%z=zvc0
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP3,EdgP7,Vertex1,Vertex2,xc0,xc1)
           P3=PointParametricEarth(EdgP3)
           P7=PointParametricEarth(EdgP7)
   
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeXaxis(EdgP2,EdgP6,Vertex1,Vertex2,xc0,xc1)
           P2=PointParametricEarth(EdgP2)
           P6=PointParametricEarth(EdgP6)
   
           Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
           !auch Vol=ABS(VolP8(P1,P5,P3,P7,P0,P4,P2,P6))
           Cell%Vol=Cell%Vol+Vol
           MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
           MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
           SumMidP=Vol*MidPoint+SumMidP
         yvc0=yvc0+dy
         END DO
       zvc0=zvc0+dz
       END DO
    ELSE IF((Cell%face3%in_out==4.AND.Cell%face4%in_out==-4).OR. &
            (Cell%face3%in_out==-4.AND.Cell%face4%in_out==4)) THEN
       !...Plane y-direction..................................
       zvc0=PMin%z+Shrink*ABS(PMin%z)
       zvc1=PMax%z-Shrink*ABS(PMax%z)
       xvc0=PMin%x+Shrink*ABS(PMin%x)
       xvc1=PMax%x-Shrink*ABS(PMax%x)
       xc0=xvc0
       !...Vergleich zu InitCellVol mit nz=1
       !...IncrVol 1Stelle nach Komma gleich, +50 (2St.), +150 (3St.),+200 (4St.)
       nxc=IncrVol ! +50 , +200
       nyc=IncrVol ! +50 , +200
       nzc=IncrVol ! +50 , +200
       dx=(xvc1-xvc0)/nxc
       dz=(zvc1-zvc0)/nzc
       Vertex1%Point%y=yc0
       Vertex1%in_out=-1
       Vertex2%Point%y=yc1
       Vertex2%in_out=1
       SumMidP%x=0.0d0
       SumMidP%y=0.0d0
       SumMidP%z=0.0d0
       DO kc=1,nzc
         xvc0=xc0
         DO ic=1,nxc
           !          6----7    !Point-Folge Vol-Berechnung
           !         /|   /|
           !        / |  / |
           !       /  2 /  3 
           !      4----5  /
           !      | /  | / 
           !      |/   |/ 
           !      0----1
           Vertex1%Point%z=zvc0
           Vertex1%Point%x=xvc0
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP0,EdgP2,Vertex1,Vertex2,yc0,yc1)
           P0=PointParametricEarth(EdgP0)
           P2=PointParametricEarth(EdgP2)
   
           Vertex1%Point%z=zvc0
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP1,EdgP3,Vertex1,Vertex2,yc0,yc1)
           P1=PointParametricEarth(EdgP1)
           P3=PointParametricEarth(EdgP3)
   
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex1%Point%x=xvc0
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP4,EdgP6,Vertex1,Vertex2,yc0,yc1)
           P4=PointParametricEarth(EdgP4)
           P6=PointParametricEarth(EdgP6)
   
           Vertex1%Point%z=MIN(zvc0+dz,zvc1)
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex2%Point%z=Vertex1%Point%z
           Vertex2%Point%x=Vertex1%Point%x
           CALL SearchPointsEdgeYaxis(EdgP5,EdgP7,Vertex1,Vertex2,yc0,yc1)
           P5=PointParametricEarth(EdgP5)
           P7=PointParametricEarth(EdgP7)
   
           Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
           Cell%Vol=Cell%Vol+Vol
           MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
           MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
           SumMidP=Vol*MidPoint+SumMidP
         xvc0=xvc0+dx
         END DO
       zvc0=zvc0+dz
       END DO
    ELSE 
       !...Plane z-direction..................................
       xvc0=PMin%x+Shrink*ABS(PMin%x)
       xvc1=PMax%x-Shrink*ABS(PMax%x)
       yvc0=PMin%y+Shrink*ABS(PMin%y)
       yvc1=PMax%y-Shrink*ABS(PMax%y)
       yc0=yvc0
       !...Vergleich zu InitCellVol mit nz=1
       !...IncrVol 1Stelle nach Komma gleich, +50 (2St.), +150 (3St.),+200 (4St.)
       nxc=IncrVol ! +50 , +200
       nyc=IncrVol ! +50 , +200
       nzc=IncrVol ! +50 , +200
       dx=(xvc1-xvc0)/nxc
       dy=(yvc1-yvc0)/nyc
       Vertex1%Point%z=zc0
       Vertex1%in_out=-1
       Vertex2%Point%z=zc1
       Vertex2%in_out=1
       SumMidP%x=0.0d0
       SumMidP%y=0.0d0
       SumMidP%z=0.0d0
       DO ic=1,nxc
         yvc0=yc0
         DO jc=1,nyc
           !      6--------7    !Point-Folge Vol-Berechnung
           !     /|       /|
           !    4--------5 |
           !    | |      | |
           !    | 2------|-3
           !    |/       |/
           !    0--------1
           Vertex1%Point%x=xvc0
           Vertex1%Point%y=yvc0
           Vertex2%Point%x=Vertex1%Point%x
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP0,EdgP4,Vertex1,Vertex2,zc0,zc1)
           P0=PointParametricEarth(EdgP0)
           P4=PointParametricEarth(EdgP4)
   
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex1%Point%y=yvc0
           Vertex2%Point%x=Vertex1%Point%x
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP1,EdgP5,Vertex1,Vertex2,zc0,zc1)
           P1=PointParametricEarth(EdgP1)
           P5=PointParametricEarth(EdgP5)
   
           Vertex1%Point%x=xvc0
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%x=Vertex1%Point%x
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP2,EdgP6,Vertex1,Vertex2,zc0,zc1)
           P2=PointParametricEarth(EdgP2)
           P6=PointParametricEarth(EdgP6)
   
           Vertex1%Point%x=MIN(xvc0+dx,xvc1)
           Vertex1%Point%y=MIN(yvc0+dy,yvc1)
           Vertex2%Point%x=Vertex1%Point%x
           Vertex2%Point%y=Vertex1%Point%y
           CALL SearchPointsEdgeZaxis(EdgP3,EdgP7,Vertex1,Vertex2,zc0,zc1)
           P3=PointParametricEarth(EdgP3)
           P7=PointParametricEarth(EdgP7)
   
           Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
           Cell%Vol=Cell%Vol+Vol
           MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
           MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
           SumMidP=Vol*MidPoint+SumMidP
         yvc0=yvc0+dy
         END DO
       xvc0=xvc0+dx
       END DO
    END IF
    Cell%MidPoint=SumMidP/Cell%Vol 
    NrMP_Cells=NrMP_Cells+1
    Cell%mp=NrMP_Cells
    IF (ncut<=2) THEN
      IF (Cell%in_out<=0) THEN
        Cell%Vol=0.0d0
      ELSE
        Cell%in_out=8
      END IF
    END IF
    IF (Cell%vc>2) THEN
      i=1
      !write(*,*) "-------"
      IF(Cell%vc>6) THEN
        Write(OutUnitProt,*) "Cell%vc>6"
        pCell=Cell
        CALL WriteCellProt(pCell,99,99,99)
      END IF     
      DO i=1,Cell%vc
        PCut(i)=VertOut(Cell%VertCut(i))%Point
       ! write(*,*) PCut(i)
      END DO
      !write(*,*) (PCut(i),i=1,Cell%vc)
      SELECT CASE (Cell%vc)
        CASE(3)
          Cell%CutF_MidP=FaceMidP3(PCut(1),PCut(2),PCut(3))
          !relativer MidPoint der Celle;  xMidP-x0/dx; yMidP-y0/dy; zMidP-z0/dz
          !Cell%CutF_MidP%x=(Cell%CutF_MidP%x-Cell%Face1%Edge1%Vert1%Point%x)/ &
          !         (Cell%Face1%Edge1%Vert2%Point%x-Cell%Face1%Edge1%Vert1%Point%x)
          !Cell%CutF_MidP%y=(Cell%CutF_MidP%y-Cell%Face1%Edge1%Vert1%Point%y)/ &
          !         (Cell%Face1%Edge4%Vert2%Point%y-Cell%Face1%Edge4%Vert1%Point%y)
          !Cell%CutF_MidP%z=(Cell%CutF_MidP%z-Cell%Face1%Edge1%Vert1%Point%z)/ &
          !         (Cell%Face3%Edge1%Vert2%Point%z-Cell%Face3%Edge1%Vert1%Point%z)
        CASE(4)
          Cell%CutF_MidP=FaceMidP4(PCut(1),PCut(2),PCut(3),PCut(4))
         ! Write(*,*) "Cell%CutF_MidP=",Cell%CutF_MidP
          !Cell%CutF_MidP%x=(Cell%CutF_MidP%x-Cell%Face1%Edge1%Vert1%Point%x)/ &
          !         (Cell%Face1%Edge1%Vert2%Point%x-Cell%Face1%Edge1%Vert1%Point%x)
          !Cell%CutF_MidP%y=(Cell%CutF_MidP%y-Cell%Face1%Edge1%Vert1%Point%y)/ &
          !         (Cell%Face1%Edge4%Vert2%Point%y-Cell%Face1%Edge4%Vert1%Point%y)
          !Cell%CutF_MidP%z=(Cell%CutF_MidP%z-Cell%Face1%Edge1%Vert1%Point%z)/ &
          !         (Cell%Face3%Edge1%Vert2%Point%z-Cell%Face3%Edge1%Vert1%Point%z)
        CASE(5)
          Cell%CutF_MidP=FaceMidP5(PCut(1),PCut(2),PCut(3),PCut(4),PCut(5))
          !Cell%CutF_MidP%x=(Cell%CutF_MidP%x-Cell%Face1%Edge1%Vert1%Point%x)/ &
          !         (Cell%Face1%Edge1%Vert2%Point%x-Cell%Face1%Edge1%Vert1%Point%x)
          !Cell%CutF_MidP%y=(Cell%CutF_MidP%y-Cell%Face1%Edge1%Vert1%Point%y)/ &
          !         (Cell%Face1%Edge4%Vert2%Point%y-Cell%Face1%Edge4%Vert1%Point%y)
          !Cell%CutF_MidP%z=(Cell%CutF_MidP%z-Cell%Face1%Edge1%Vert1%Point%z)/ &
          !        (Cell%Face3%Edge1%Vert2%Point%z-Cell%Face3%Edge1%Vert1%Point%z)
        CASE(6)
          Cell%CutF_MidP=FaceMidP6(PCut(1),PCut(2),PCut(3),PCut(4),PCut(5),PCut(6))
          !Cell%CutF_MidP%x=(Cell%CutF_MidP%x-Cell%Face1%Edge1%Vert1%Point%x)/ &
          !         (Cell%Face1%Edge1%Vert2%Point%x-Cell%Face1%Edge1%Vert1%Point%x)
          !Cell%CutF_MidP%y=(Cell%CutF_MidP%y-Cell%Face1%Edge1%Vert1%Point%y)/ &
          !         (Cell%Face1%Edge4%Vert2%Point%y-Cell%Face1%Edge4%Vert1%Point%y)
          !Cell%CutF_MidP%z=(Cell%CutF_MidP%z-Cell%Face1%Edge1%Vert1%Point%z)/ &
          !         (Cell%Face3%Edge1%Vert2%Point%z-Cell%Face3%Edge1%Vert1%Point%z)
      END SELECT
    END IF
  ELSE  ! Cell%vc==0
    IF ((Cell%Face1%Edge1%Vert1%in_out==0 .AND. &
         Cell%Face1%Edge1%Vert2%in_out==0 .AND. &
         Cell%Face1%Edge3%Vert1%in_out==0 .AND. &
         Cell%Face1%Edge3%Vert2%in_out==0) &
        .AND. Cell%Face2%in_out==4) THEN 
       !Grenzflaeche Face1
       EdgP0=Cell%Face1%Edge1%Vert1%Point
       EdgP1=Cell%Face1%Edge1%Vert2%Point
       EdgP2=Cell%Face1%Edge3%Vert1%Point
       EdgP3=Cell%Face1%Edge3%Vert2%Point
       EdgP4=Cell%Face2%Edge1%Vert1%Point
       EdgP5=Cell%Face2%Edge1%Vert2%Point
       EdgP6=Cell%Face2%Edge3%Vert1%Point
       EdgP7=Cell%Face2%Edge3%Vert2%Point
       P0=PointParametricEarth(EdgP0)
       P1=PointParametricEarth(EdgP1)
       P2=PointParametricEarth(EdgP2)
       P3=PointParametricEarth(EdgP3)
       P4=PointParametricEarth(EdgP4)
       P5=PointParametricEarth(EdgP5)
       P6=PointParametricEarth(EdgP6)
       P7=PointParametricEarth(EdgP7)
       Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
       Cell%MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
       Cell%CutF_MidP=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3)
       !!relativer MidPoint der Celle;  xMidP-x0/dx; yMidP-y0/dy; zMidP-z0/dz
       !Cell%CutF_MidP%x=(Cell%CutF_MidP%x-Cell%Face1%Edge1%Vert1%Point%x)/ &
       !       (Cell%Face1%Edge1%Vert2%Point%x-Cell%Face1%Edge1%Vert1%Point%x)
       !Cell%CutF_MidP%y=(Cell%CutF_MidP%y-Cell%Face1%Edge1%Vert1%Point%y)/ &
       !       (Cell%Face1%Edge4%Vert2%Point%y-Cell%Face1%Edge4%Vert1%Point%y)
       !Cell%CutF_MidP%z=(Cell%CutF_MidP%z-Cell%Face1%Edge1%Vert1%Point%z)/ &
       !       (Cell%Face3%Edge1%Vert2%Point%z-Cell%Face3%Edge1%Vert1%Point%z)
    ELSE
       !fuer Cellen mit Grenzfläche Edge.OR.Vertex
       !diese MP-Berechnung wird benötigt wenn RandCellen berechnet
       !und Addition höherer Verfeinerung gebraucht wird
       EdgP0=Cell%Face1%Edge1%Vert1%Point
       EdgP1=Cell%Face1%Edge1%Vert2%Point
       EdgP2=Cell%Face1%Edge3%Vert1%Point
       EdgP3=Cell%Face1%Edge3%Vert2%Point
       EdgP4=Cell%Face2%Edge1%Vert1%Point
       EdgP5=Cell%Face2%Edge1%Vert2%Point
       EdgP6=Cell%Face2%Edge3%Vert1%Point
       EdgP7=Cell%Face2%Edge3%Vert2%Point
       P0=PointParametricEarth(EdgP0)
       P1=PointParametricEarth(EdgP1)
       P2=PointParametricEarth(EdgP2)
       P3=PointParametricEarth(EdgP3)
       P4=PointParametricEarth(EdgP4)
       P5=PointParametricEarth(EdgP5)
       P6=PointParametricEarth(EdgP6)
       P7=PointParametricEarth(EdgP7)
       Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
       Cell%MidPoint=VolMidP8(EdgP0,EdgP1,EdgP2,EdgP3,EdgP4,EdgP5,EdgP6,EdgP7)
       Cell%CutF_MidP=FaceMidP4(EdgP0,EdgP1,EdgP2,EdgP3)
       !!relativer MidPoint der Celle;  xMidP-x0/dx; yMidP-y0/dy; zMidP-z0/dz
       !Cell%CutF_MidP%x=(Cell%CutF_MidP%x-Cell%Face1%Edge1%Vert1%Point%x)/ &
       !       (Cell%Face1%Edge1%Vert2%Point%x-Cell%Face1%Edge1%Vert1%Point%x)
       !Cell%CutF_MidP%y=(Cell%CutF_MidP%y-Cell%Face1%Edge1%Vert1%Point%y)/ &
       !       (Cell%Face1%Edge4%Vert2%Point%y-Cell%Face1%Edge4%Vert1%Point%y)
       !Cell%CutF_MidP%z=(Cell%CutF_MidP%z-Cell%Face1%Edge1%Vert1%Point%z)/ &
       !       (Cell%Face3%Edge1%Vert2%Point%z-Cell%Face3%Edge1%Vert1%Point%z)
    END IF
  END IF
  !
  ! in 7.7.1.3 enthalten, mit gridfile.grid testen, wichtig fuer Cutplane-Check
  !Celle mit Grenzflaeche Face2, Cell%Vol =0 setzen
  !wird zuvor berechnet, da Face%ec -->Cell%vc=4 ->Vol berechnet
  !IF (MAX(Cell%Face1%Vol,Cell%Face2%Vol,Cell%Face3%Vol,Cell%Face4%Vol,Cell%Face5%Vol,Cell%Face6%Vol)<=0.0d0) THEN
  !schliesst folgende Anweisung mit ein
  !IF (Cell%Face2%Edge1%Vert1%in_out==0 .AND. Cell%Face2%Edge1%Vert2%in_out==0 .AND. &
  !    Cell%Face2%Edge3%Vert1%in_out==0 .AND. Cell%Face2%Edge3%Vert2%in_out==0 .AND. &
  !    Cell%Face1%in_out==-4) THEN
  !    Cell%Vol=0.0d0
  !END IF
 
  IF (Cell%Vol==0.0d0) THEN
    Cell%Face1%Vol=0.0d0
    Cell%Face2%Vol=0.0d0
    Cell%Face3%Vol=0.0d0
    Cell%Face4%Vol=0.0d0
    Cell%Face5%Vol=0.0d0
    Cell%Face6%Vol=0.0d0
    IF(Cell%mp>0)THEN
      NrMP_Cells=NrMP_Cells-1
      Cell%mp=0
    END IF
  END IF
  IF (MAX(Cell%Face1%Vol,Cell%Face2%Vol,Cell%Face3%Vol,Cell%Face4%Vol,Cell%Face5%Vol,Cell%Face6%Vol)<=0.0d0) THEN
    Cell%Vol=0.0d0  !Ossi
    Cell%mp=0  ! Celle unterhalb Berg angrenzt 
  END IF
 
END SUBROUTINE AnalyzeCell

SUBROUTINE Prot_AnalyzeAllCells
   Write(OutUnitProt,*) "........................................"
   Write(OutUnitProt,*) "aus AnalyzeAllCell:"
   Write(OutUnitProt,*) "~~~~~~~~~~~~~~~~~~~"
   Write(OutUnitProt,*) "NrW_Cells             =", NrW_Cells
   Write(OutUnitProt,*) "NrW_Cells1            =", NrW_Cells1
   !nicht global Write(OutUnitProt,*) "NrW_Cells1U           =", NrW_Cells1U
   Write(OutUnitProt,*) "NrW_Cells2            =", NrW_Cells2
   Write(OutUnitProt,*) "nr_cells              =", nr_cells
   Write(OutUnitProt,*) "nr_viewcells          =", nr_viewcells
   Write(OutUnitProt,*) "NrB_Cells             =", NrB_Cells
   Write(OutUnitProt,*) "nr_cutplanecells      =", nr_cutplanecells
   Write(OutUnitProt,*) "nr_cutplanecells1     =", nr_cutplanecells1
   Write(OutUnitProt,*) "nr_viewcutplanecells  =", nr_viewcutplanecells
   Write(OutUnitProt,*) "nr_viewcutplanecells1 =", nr_viewcutplanecells1
   Write(OutUnitProt,*) "nr_cutplanecells2     =", nr_cutplanecells2
   Write(OutUnitProt,*) "nr_viewcutplanecells2 =", nr_viewcutplanecells2
   Write(OutUnitProt,*) "nr_cutcells           =", nr_cutcells
   Write(OutUnitProt,*) "nr_viewcutcells       =", nr_viewcutcells
   Write(OutUnitProt,*) "NrAll_Cells           =", NrAll_Cells
   Write(OutUnitProt,*) "nr_incells            =", nr_incells
   Write(OutUnitProt,*) "nr_viewincells        =", nr_viewincells
   Write(OutUnitProt,*) "nr_orocells           =", nr_orocells
   Write(OutUnitProt,*) "nr_vieworocells       =", nr_vieworocells
   Write(OutUnitProt,*) "........................................"
END SUBROUTINE Prot_AnalyzeAllCells

SUBROUTINE Screen_Prot_AnalyzeAllCells
   Write(*,*) "........................................"
   Write(*,*) "aus AnalyzeAllCell:"
   Write(*,*) "~~~~~~~~~~~~~~~~~~~"
   Write(*,*) "NrW_Cells             =", NrW_Cells
   Write(*,*) "NrW_Cells1            =", NrW_Cells1
   !nicht global Write(*,*) "NrW_Cells1U           =", NrW_Cells1U
   Write(*,*) "NrW_Cells2            =", NrW_Cells2
   Write(*,*) "nr_cells              =", nr_cells
   Write(*,*) "nr_viewcells          =", nr_viewcells
   Write(*,*) "NrB_Cells             =", NrB_Cells
   Write(*,*) "nr_cutplanecells      =", nr_cutplanecells
   Write(*,*) "nr_cutplanecells1     =", nr_cutplanecells1
   Write(*,*) "nr_viewcutplanecells  =", nr_viewcutplanecells
   Write(*,*) "nr_viewcutplanecells1 =", nr_viewcutplanecells1
   Write(*,*) "nr_cutplanecells2     =", nr_cutplanecells2
   Write(*,*) "nr_viewcutplanecells2 =", nr_viewcutplanecells2
   Write(*,*) "nr_cutcells           =", nr_cutcells
   Write(*,*) "nr_viewcutcells       =", nr_viewcutcells
   Write(*,*) "NrAll_Cells           =", NrAll_Cells
   Write(*,*) "nr_incells            =", nr_incells
   Write(*,*) "nr_viewincells        =", nr_viewincells
   Write(*,*) "nr_orocells           =", nr_orocells
   Write(*,*) "nr_vieworocells       =", nr_vieworocells
   Write(*,*) "........................................"
END SUBROUTINE Screen_Prot_AnalyzeAllCells

SUBROUTINE AnalyzeAllCells
  INTEGER :: ib,i,j,k,liv1,liv2,gl
  INTEGER :: v_x,v_y,v_z,pix
  INTEGER :: in_out,in_out_view
  INTEGER :: NrW_Cells1U

  nr_cells=0
  nr_cutplanecells=0
  nr_cutplanecells1=0
  nr_cutplanecells2=0
  nr_cutcells=0
  nr_cutf1=0
  nr_soilplanecells=0
  nr_soilplanecells1=0
  nr_soilplanecells2=0
  nr_incells=0
  nr_viewcells=0
  nr_viewcutplanecells=0
  nr_viewcutplanecells1=0
  nr_viewcutplanecells2=0
  nr_viewcutcells=0
  nr_viewcutf1=0
  nr_viewsoilplanecells=0
  nr_viewsoilplanecells1=0
  nr_viewsoilplanecells2=0
  nr_viewincells=0
  NrAll_Cells=0
  NrW_All_Cells=0
  NrB_All_Cells=0
  NrMP_All_Cells=0

  !.....................................................................
  DO ib=1,nb
    CALL Set(Floor(ib))
    NrW_Cells1U=0
    !Write(*,*) "         ...AnalyzeAllCells inside Block : ", ib
    Write(*,*) "                                     Block : ", ib,"\/",nb
    Write(*,'(A7,$)') "       "
    if(mod(ix1,10)==0)then
       pix=ix1/10  
    else  
       pix=ABS(ix1/10)+1
    end if
    DO i=ix0+1,ix1
      ! Prozentausgabe
       !if (mod((100*i)/ix1,10)==0.and.(mod((100*(i+1))/ix1,10)/=0.or.i==ix1)) &
       !  !&  write(*,*) "                      ",(100*i)/ix1,"%"  !Detlef
       if(mod(i,pix)==0.and.i/=ix1) then 
        ! write(*,'(I5,1A,$)')(i*100)/ix1,"%"       ! mit Formatelement $-> Zeilenvorschub unterdrücken
        ! write(*,'(I5,1A)',advance='NO') (i*100)/ix1,"%"  ! mit expliziten Formatspezifizierer
        ! Kontrollformatelement tln-Zurücksetzen um n Stellen von der momentanen Position
        ! erfordert ein Formatelement vorraus Bsp. I3, tl3 
        ! wenn Formatspezifizierer advance='NO' mit angegeben ist funktioniert es nicht 
         write(*,'(I5,1A)',advance='NO')(i*100)/ix1,"%"
       else if (i==ix1)then 
         write(*,'(I5,1A)')(i*100)/ix1,"%"
       end if 
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
             !Test Compiler Print-Steuerung unter Linux noch problem
             !#ifdef _PANALYZE
               !WRITE(*,*) "Block : ",ib,"  Cell :  i=",i," j=",j," k=",k
             !#endif
             !Write(*,*) "AnalyzeCell(Cell(i,j,k)%Cell=",i,j,k
             !.................................
             CALL AnalyzeCell(Cell(i,j,k)%Cell)
             !IF(i==2.and.j==3.and.k==1)  THEN
                ! CALL WriteCellProt(Cell(i,j,k)%Cell,i,j,k)
             !END IF
             !.............................
             !Bereich View-Cells GMV zählen
             IF(RefineX<0 .AND. i>1) THEN
                v_x=i*RefineX*(-2)
             ELSE
                v_x=i
             END IF
             !...........................
             IF(RefineY<0 .AND. j>1) THEN
                v_y=j*RefineY*(-2)
             ELSE
                v_y=j
             END IF
             !...........................
             IF(RefineZ<0 .AND. k>1) THEN
                v_z=k*RefineZ*(-2)
             ELSE
                v_z=k
             END IF
             !...........................
             !geschnittene- sowie komplett angrenzende Cellen (Flaechen:F1,F2)
             !oberhalb, unterhalb am Berg
             !IF (Cell(i,j,k)%Cell%in_out<8) THEN
             IF (Cell(i,j,k)%Cell%in_out>6.OR. &
                 (Cell(i,j,k)%Cell%in_out==6.AND.Cell(i,j,k)%Cell%vc==0)) THEN
                 !(Cell(i,j,k)%Cell%in_out==6.AND.Cell(i,j,k)%Cell%vc==0).OR. &
                 !(Cell(i,j,k)%Cell%in_out==4.AND.Cell(i,j,k)%Cell%vc==0)) THEN
                 !Cellen mit Face od. Edge od. Point an Berg
                !Cellen mit Edge od. Point an Berg  
                !Cellen mit MaxVol, nicht unter NrW_ ausgeben
             ELSE IF (Cell(i,j,k)%Cell%in_out<6.AND.Cell(i,j,k)%Cell%Vol>0.0d0 .OR. &
                 (Cell(i,j,k)%Cell%in_out==6.AND.Cell(i,j,k)%Cell%vc>0) ) THEN
               NrW_Cells=NrW_Cells+1
               NrW_Cells1=NrW_Cells1+1
               !Write(OutUnitProt,*) i,j,k, "  NrW_Cells=",NrW_Cells," ----NrWCells1....."
             ELSE
                !Bsp./s1agi/reutgen/WRK_GITTER/TEST_GRID_GMV_8.1.1.2/Agi_2Blk_H5000_test.grid
                !    (2.Block Celle 6,1,1 unterhalb am Berg wird allokiert)
                ! Cellen mit Kanten unterhalb am Berg noch filtern in AnalyzeCell
               NrW_Cells=NrW_Cells+1
               NrW_Cells1U=NrW_Cells1U+1
               !Write(OutUnitProt,*) i,j,k, "  NrW_Cells=",NrW_Cells," ----NrWCells1U....."
             END IF

             !geschnittene- und angrenzende Cellen [Vol>0]
             !oberhalb Berg 
             IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
                 .OR.Cell(i,j,k)%Cell%in_out==8) THEN
               nr_cells=nr_cells+1
               !nr_cutplanecells=nr_cutplanecells+1
               IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                 IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                   IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                     nr_viewcells=nr_viewcells+1
                     !nr_viewcutplanecells=nr_viewcutplanecells+1                     
                   END IF
                 END IF
               END IF
             END IF
             !geschnittene Cellen [Vol>0,Vol<maxVol]
             !oberhalb Berg
             IF(Cell(i,j,k)%Cell%in_out==5.AND.Cell(i,j,k)%Cell%vc==0) THEN
               ! Special wird in Soil nicht ausgegeben
               ! bessere Lösung zur Anpassung ? 
             ELSE
               IF (Cell(i,j,k)%Cell%in_out<8.AND.Cell(i,j,k)%Cell%Vol>0.0d0.AND. &
                   (Cell(i,j,k)%Cell%Face1%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face1)<1.d0-1.d-2 &
                .OR.Cell(i,j,k)%Cell%Face2%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face2)<1.d0-1.d-2 &
                .OR.Cell(i,j,k)%Cell%Face3%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face3)<1.d0-1.d-2 &
                .OR.Cell(i,j,k)%Cell%Face4%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face4)<1.d0-1.d-2 &
                .OR.Cell(i,j,k)%Cell%Face5%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face5)<1.d0-1.d-2 &
                .OR.Cell(i,j,k)%Cell%Face6%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face6)<1.d0-1.d-2)) THEN
                 !NrB_Cells=NrB_Cells+1
               END IF
             END IF
             !angrenzende Cellen(Face2), [Vol==0]
             !unterhalb Berg
             IF (Cell(i,j,k)%Cell%in_out==-4.AND.Cell(i,j,k)%Cell%Vol==0.0d0 &
                 .AND.Cell(i,j,k)%Cell%vc==4) THEN
               !nr_cutplanecells=nr_cutplanecells+1
               IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                 IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                   IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                     !nr_viewcutplanecells=nr_viewcutplanecells+1                     
                   END IF
                 END IF
               END IF
             END IF
             IF (Cell(i,j,k)%Cell%in_out==-5.AND.Cell(i,j,k)%Cell%Vol==0.0d0 &
                 .AND.Cell(i,j,k)%Cell%vc==3) THEN
               !nr_cutplanecells=nr_cutplanecells+1
               IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                 IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                   IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                     !nr_viewcutplanecells=nr_viewcutplanecells+1                     
                   END IF
                 END IF
               END IF
             END IF
             !.....................................................................
             !alle Schnittflaechen Cellen
             IF (Cell(i,j,k)%Cell%vc>0) THEN
               !Standard-Schnittflaeche
               nr_cutplanecells=nr_cutplanecells+1
               nr_cutplanecells1=nr_cutplanecells1+1 
               IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                 IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                   IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                     nr_viewcutplanecells=nr_viewcutplanecells+1
                     nr_viewcutplanecells1=nr_viewcutplanecells1+1 
                   END IF
                 END IF
               END IF

               IF(k==1) THEN
                  gl=0
                  IF(Cell(i,j,k)%Cell%vc== & 
                     Cell(i,j,k)%Cell%Face1%NumberVert) THEN
                    DO liv1=1,Cell(i,j,k)%Cell%vc
                     DO liv2=1,Cell(i,j,k)%Cell%vc
                       IF(Cell(i,j,k)%Cell%VertCut(liv1)== &
                          Cell(i,j,k)%Cell%Face1%VertexList(liv2)) THEN
                         gl=gl+1
                       END IF
                     END DO
                    END DO
                  END IF
               IF(gl/=Cell(i,j,k)%Cell%vc) THEN

                 IF(Cell(i,j,k)%Cell%Face1%Edge1%yes_sp==1.OR. &
                    Cell(i,j,k)%Cell%Face1%Edge2%yes_sp==1.OR. &
                    Cell(i,j,k)%Cell%Face1%Edge3%yes_sp==1.OR. &
                    Cell(i,j,k)%Cell%Face1%Edge4%yes_sp==1) THEN
                    !wenn Face1 Schnitt und Grenzflaeche-Gelaende
                    !Bsp: Durran [139,162],1,1; Valley 41,[46...100],1
                    nr_cutplanecells=nr_cutplanecells+1
                    nr_cutplanecells2=nr_cutplanecells2+1
                    IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                      IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                        IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                          nr_viewcutplanecells=nr_viewcutplanecells+1
                          nr_viewcutplanecells2=nr_viewcutplanecells2+1
                        END IF
                      END IF
                    END IF
                 ELSE IF(Cell(i,j,k)%Cell%Face1%in_out==-1) THEN
                    nr_cutplanecells=nr_cutplanecells+1
                    nr_cutplanecells2=nr_cutplanecells2+1
                    IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                      IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                        IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                          nr_viewcutplanecells=nr_viewcutplanecells+1
                          nr_viewcutplanecells2=nr_viewcutplanecells2+1
                        END IF
                      END IF
                    END IF
                 END IF
               END IF ! Vert-Ungleich-Check
               END IF ! k==1
             ELSE
                !IF(Cell(i,j,k)%Cell%in_out==4.AND. &
                IF(k==1.AND. & ! in_out==4 schließt k==1 mit ein
                                ! fehlt unter Cut planecells wenn nicht k==1 ist
                   (Cell(i,j,k)%Cell%Face1%Edge1%Vert1%in_out==0.AND.&
                    Cell(i,j,k)%Cell%Face1%Edge1%Vert2%in_out==0.AND.&
                    Cell(i,j,k)%Cell%Face1%Edge3%Vert1%in_out==0.AND.&
                    Cell(i,j,k)%Cell%Face1%Edge3%Vert2%in_out==0) ) THEN
                   !Face1, Grenzflaeche-Gelaende
                   nr_cutplanecells=nr_cutplanecells+1
                   IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                     IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                       IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                          nr_viewcutplanecells=nr_viewcutplanecells+1
                          !Write(OutUnitProt,*) "i,j,k=",i,j,k
                       END IF
                     END IF
                   END IF
                END IF
             END IF
             !...........................................................................
             ! Soil-count: Übernommen aus Soil-Output !! mit alllen specials
             IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
                IF (((Cell(i,j,k)%Cell%in_out<5.AND.Cell(i,j,k)%Cell%Vol>0.0d0) .OR. &
                     (Cell(i,j,k)%Cell%in_out==6.AND.Cell(i,j,k)%Cell%vc>0) .OR. &
                     (Cell(i,j,k)%Cell%in_out==5.AND.Cell(i,j,k)%Cell%vc>0)) &
                    .AND. &
                    (Cell(i,j,k)%Cell%Face1%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face1) &
                       & <1.d0-dist_scMaxCell     .OR. &
                     Cell(i,j,k)%Cell%Face2%Vol/ VolFace_XY(Cell(i,j,k)%Cell%Face2) &
                       & <1.d0-dist_scMaxCell     .OR. &
                     Cell(i,j,k)%Cell%Face3%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face3) &
                       & <1.d0-dist_scMaxCell     .OR. &
                     Cell(i,j,k)%Cell%Face4%Vol/ VolFace_ZX(Cell(i,j,k)%Cell%Face4) &
                       & <1.d0-dist_scMaxCell     .OR. &
                     Cell(i,j,k)%Cell%Face5%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face5) &
                       & <1.d0-dist_scMaxCell     .OR. &
                     Cell(i,j,k)%Cell%Face6%Vol/ VolFace_YZ(Cell(i,j,k)%Cell%Face6) &
                       & <1.d0-dist_scMaxCell) &
                   ) THEN
                     NrB_Cells=NrB_Cells+1
                     nr_soilplanecells=nr_soilplanecells+1
                     nr_soilplanecells1=nr_soilplanecells1+1
                     IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                       IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                         IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                            nr_viewsoilplanecells=nr_viewsoilplanecells+1
                            nr_viewsoilplanecells1=nr_viewsoilplanecells1+1 
                             !Write(OutUnitProt,*) "i,j,k=",i,j,k, "Standard_GrenzeSoil"
                         END IF
                       END IF
                     END IF
                ELSE
                   !IF(Cell(i,j,k)%Cell%in_out==4.AND. &
                   IF(k==1.AND. & ! in_out==4 schließt k==1 mit ein
                                  ! fehlt unter Cut planecells wenn nicht k==1 ist
                      (Cell(i,j,k)%Cell%Face1%Edge1%Vert1%in_out==0.AND.&
                       Cell(i,j,k)%Cell%Face1%Edge1%Vert2%in_out==0.AND.&
                       Cell(i,j,k)%Cell%Face1%Edge3%Vert1%in_out==0.AND.&
                       Cell(i,j,k)%Cell%Face1%Edge3%Vert2%in_out==0) ) THEN
                      !Face1, Grenzflaeche-Gelaende
                      NrB_Cells=NrB_Cells+1
                      nr_soilplanecells=nr_soilplanecells+1
                      IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                        IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                          IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                             nr_viewsoilplanecells=nr_viewsoilplanecells+1
                             !Write(OutUnitProt,*) "i,j,k=",i,j,k, "Face1_GrenzeSoil"
                          END IF
                        END IF
                      END IF
                   END IF
                END IF  ! Soil-Std und Soil-Special
             END IF   ! IF (ASSOCIATED....Celle...
             !...........................................................................
             !geschnittene- und unterhalb angrenzende Cellen
             !am Berg fuer Oro-Darstellung
             ! (Ausblenden Cellen oberhalb Berg mit  Cell%in_out==[6,7])

             IF(Cell(i,j,k)%Cell%vc>0 ) THEN
               nr_cutcells=nr_cutcells+1
               IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                  IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                     nr_viewcutcells=nr_viewcutcells+1
                     !Write(OutUnitProt,*) "cut=  i j k : ", i, j,k
                  END IF
                END IF
               END IF
             ELSE ! Special 1-3, vc==0
               IF(k==1.AND. (Cell(i,j,k)%Cell%Face1%Edge1%Vert1%in_out==0.AND.&
                             Cell(i,j,k)%Cell%Face1%Edge1%Vert2%in_out==0.AND.&
                             Cell(i,j,k)%Cell%Face1%Edge3%Vert1%in_out==0.AND.&
                             Cell(i,j,k)%Cell%Face1%Edge3%Vert2%in_out==0) ) THEN
                 ! Special 1
                 nr_cutf1=nr_cutf1+1
                 IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                  IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                    IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                       nr_viewcutf1=nr_viewcutf1+1
                       !Write(OutUnitProt,*) "gf1=   i j k : ", i, j,k
                    END IF
                  END IF
                 END IF
               ELSE IF ((MAX(Cell(i,j,k)%Cell%Face2%Edge1%Vert1%in_out, &
                             Cell(i,j,k)%Cell%Face2%Edge1%Vert2%in_out, &
                             Cell(i,j,k)%Cell%Face2%Edge3%Vert1%in_out, &
                             Cell(i,j,k)%Cell%Face2%Edge3%Vert2%in_out)==0)&
                        .OR. &
                        (Cell(i,j,k)%Cell%Face2%in_out==-4 .AND. &
                         Cell(i,j,k)%Cell%Face1%in_out==-1 .OR. Cell(i,j,k)%Cell%Face1%in_out==-2)) THEN
                  ! Special 2 und 3
                  nr_cutcells=nr_cutcells+1
                  IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                   IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                     IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                        nr_viewcutcells=nr_viewcutcells+1
                        !Write(OutUnitProt,*) "cut=  i j k : ", i, j,k
                     END IF
                   END IF
                  END IF
               END IF  ! IF (Special 1-3)

             END IF ! allokiert und  (cell%vc>0 oder cell%vc==0)
             NrAll_Cells=NrAll_Cells+1
          !.....................................................................
          !LandDefinition setzen
            CALL SetLandClass(Cell(i,j,k)%Cell,i,j,k,ib)

          !.....................................................................
          ELSE  !!   NichtAllokierte Cellen
             in_out=Vertices(i-1,j-1,k-1)%in_out &
                   +Vertices(i,j-1,k-1)%in_out &
                   +Vertices(i-1,j-1,k)%in_out &
                   +Vertices(i-1,j,k-1)%in_out &
                   +Vertices(i-1,j,k)%in_out   &
                   +Vertices(i,j-1,k)%in_out   &
                   +Vertices(i,j,k-1)%in_out   &
                   +Vertices(i,j,k)%in_out
             !...........................
             IF(RefineX<0 .AND. i>1) THEN
                v_x=i*RefineX*(-2)
             ELSE
                v_x=i
             END IF
             !...........................
             IF(RefineY<0 .AND. j>1) THEN
                v_y=j*RefineY*(-2)
             ELSE
                v_y=j
             END IF
             !...........................
             IF(RefineZ<0 .AND. k>1) THEN
                v_z=k*RefineZ*(-2)
             ELSE
                v_z=k
             END IF
             !...........................
             IF (in_out<0) THEN         !inside Berg
               NrW_Cells=NrW_Cells+1
               NrW_Cells2=NrW_Cells2+1
               ! Write(OutUnitProt,*) i,j,k, "  NrW_Cells=",NrW_Cells," ---- inside Berg"
             END IF
             IF (in_out>0) THEN         !outside Berg
             nr_cells=nr_cells+1
               IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                 IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                   IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                     nr_viewcells=nr_viewcells+1
                   END IF
                 END IF
               END IF
             END IF
             IF (in_out==-8) THEN        !only 'hex' inside Berg
               nr_incells=nr_incells+1
               IF(v_x>Domain%view_ixa .AND. v_x<=Domain%view_ixe) THEN
                 IF(v_y>Domain%view_iya .AND. v_y<=Domain%view_iye) THEN
                   IF(v_z>Domain%view_iza .AND. v_z<=Domain%view_ize) THEN
                     nr_viewincells=nr_viewincells+1
                   END IF
                 END IF
               END IF
             END IF
             !..............................................................
          END IF  ! IF/ELSE (ASSOCIATED(Cell)) 
        END DO  ! k
      END DO   ! j
    END DO   ! i
    Floor(ib)%NrW_Cells=NrW_Cells
    Floor(ib)%NrMP_Cells=NrB_Cells   !Checken andere Variante, beteff RandCellen Block !?
    Floor(ib)%NrB_Cells=NrB_Cells
    NrW_All_Cells=NrW_All_Cells+NrW_Cells
    NrMP_All_Cells=NrMP_All_Cells+NrB_Cells
    NrB_All_Cells=NrB_All_Cells+NrB_Cells

    CALL SetLandCoverCellsCLC ! Test

  END DO  !ib
  nr_vieworocells=nr_viewincells+nr_viewcutcells+nr_viewcutf1
  nr_orocells=nr_incells+nr_cutcells+nr_cutf1


   CALL WriteAuswAnalyzeAllCellsToProt
   CALL Prot_AnalyzeAllCells
   CALL WriteEndAnalyzeAllCellsToProt
   !CALL Screen_Prot_AnalyzeAllCells

END SUBROUTINE AnalyzeAllCells


SUBROUTINE SearchMinMaxP8Cell(Cell,PMin,PMax)
  TYPE(CellP8_T) :: Cell
  TYPE(Point_T) :: PMin,PMax
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
   
  IF(Cell%V0%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Cell%V0%Point%x)
    PMin%y=MIN(PMin%y,Cell%V0%Point%y)
    PMin%z=MIN(PMin%z,Cell%V0%Point%z)
    PMax%x=MAX(PMax%x,Cell%V0%Point%x)
    PMax%y=MAX(PMax%y,Cell%V0%Point%y)
    PMax%z=MAX(PMax%z,Cell%V0%Point%z)
  END IF
  IF(Cell%V1%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Cell%V1%Point%x)
    PMin%y=MIN(PMin%y,Cell%V1%Point%y)
    PMin%z=MIN(PMin%z,Cell%V1%Point%z)
    PMax%x=MAX(PMax%x,Cell%V1%Point%x)
    PMax%y=MAX(PMax%y,Cell%V1%Point%y)
    PMax%z=MAX(PMax%z,Cell%V1%Point%z)
  END IF
  IF(Cell%V2%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Cell%V2%Point%x)
    PMin%y=MIN(PMin%y,Cell%V2%Point%y)
    PMin%z=MIN(PMin%z,Cell%V2%Point%z)
    PMax%x=MAX(PMax%x,Cell%V2%Point%x)
    PMax%y=MAX(PMax%y,Cell%V2%Point%y)
    PMax%z=MAX(PMax%z,Cell%V2%Point%z)
  END IF
  IF(Cell%V3%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Cell%V3%Point%x)
    PMin%y=MIN(PMin%y,Cell%V3%Point%y)
    PMin%z=MIN(PMin%z,Cell%V3%Point%z)
    PMax%x=MAX(PMax%x,Cell%V3%Point%x)
    PMax%y=MAX(PMax%y,Cell%V3%Point%y)
    PMax%z=MAX(PMax%z,Cell%V3%Point%z)
  END IF
  IF(Cell%V4%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Cell%V4%Point%x)
    PMin%y=MIN(PMin%y,Cell%V4%Point%y)
    PMin%z=MIN(PMin%z,Cell%V4%Point%z)
    PMax%x=MAX(PMax%x,Cell%V4%Point%x)
    PMax%y=MAX(PMax%y,Cell%V4%Point%y)
    PMax%z=MAX(PMax%z,Cell%V4%Point%z)
  END IF
  IF(Cell%V5%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Cell%V5%Point%x)
    PMin%y=MIN(PMin%y,Cell%V5%Point%y)
    PMin%z=MIN(PMin%z,Cell%V5%Point%z)
    PMax%x=MAX(PMax%x,Cell%V5%Point%x)
    PMax%y=MAX(PMax%y,Cell%V5%Point%y)
    PMax%z=MAX(PMax%z,Cell%V5%Point%z)
  END IF
  IF(Cell%V6%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Cell%V6%Point%x)
    PMin%y=MIN(PMin%y,Cell%V6%Point%y)
    PMin%z=MIN(PMin%z,Cell%V6%Point%z)
    PMax%x=MAX(PMax%x,Cell%V6%Point%x)
    PMax%y=MAX(PMax%y,Cell%V6%Point%y)
    PMax%z=MAX(PMax%z,Cell%V6%Point%z)
  END IF
  IF(Cell%V7%in_out>=0) THEN
    PMin%x=MIN(PMin%x,Cell%V7%Point%x)
    PMin%y=MIN(PMin%y,Cell%V7%Point%y)
    PMin%z=MIN(PMin%z,Cell%V7%Point%z)
    PMax%x=MAX(PMax%x,Cell%V7%Point%x)
    PMax%y=MAX(PMax%y,Cell%V7%Point%y)
    PMax%z=MAX(PMax%z,Cell%V7%Point%z)
  END IF
    !................................
    ! Analyze Schnittpunkte
    ! EdgeX
    CALL SearchPointsEdgeXaxis(P0,P1,Cell%V0,Cell%V1, &
                               Cell%V0%Point%x,Cell%V1%Point%x)
    CALL SearchPointsEdgeXaxis(P2,P3,Cell%V2,Cell%V3, &
                               Cell%V2%Point%x,Cell%V3%Point%x)
    CALL SearchPointsEdgeXaxis(P4,P5,Cell%V4,Cell%V5, &
                               Cell%V4%Point%x,Cell%V5%Point%x)
    CALL SearchPointsEdgeXaxis(P6,P7,Cell%V6,Cell%V7, &
                               Cell%V6%Point%x,Cell%V7%Point%x)
  IF(P0%x/=P1%x.AND.P0%x/=Cell%V0%Point%x) THEN
    PMin%x=MIN(PMin%x,P0%x)
    PMax%x=MAX(PMax%x,P0%x)
  END IF
  IF(P0%x/=P1%x.AND.P1%x/=Cell%V1%Point%x) THEN
    PMin%x=MIN(PMin%x,P1%x)
    PMax%x=MAX(PMax%x,P1%x)
  END IF

  IF(P2%x/=P3%x.AND.P2%x/=Cell%V2%Point%x) THEN
    PMin%x=MIN(PMin%x,P2%x)
    PMax%x=MAX(PMax%x,P2%x)
  END IF
  IF(P2%x/=P3%x.AND.P3%x/=Cell%V3%Point%x) THEN
    PMin%x=MIN(PMin%x,P3%x)
    PMax%x=MAX(PMax%x,P3%x)
  END IF

  IF(P4%x/=P5%x.AND.P4%x/=Cell%V4%Point%x) THEN
    PMin%x=MIN(PMin%x,P4%x)
    PMax%x=MAX(PMax%x,P4%x)
  END IF
  IF(P4%x/=P5%x.AND.P5%x/=Cell%V5%Point%x) THEN
    PMin%x=MIN(PMin%x,P5%x)
    PMax%x=MAX(PMax%x,P5%x)
  END IF

  IF(P6%x/=P7%x.AND.P6%x/=Cell%V6%Point%x) THEN
    PMin%x=MIN(PMin%x,P6%x)
    PMax%x=MAX(PMax%x,P6%x)
  END IF
  IF(P6%x/=P7%x.AND.P7%x/=Cell%V7%Point%x) THEN
    PMin%x=MIN(PMin%x,P7%x)
    PMax%x=MAX(PMax%x,P7%x)
  END IF

    ! EdgeY
    CALL SearchPointsEdgeYaxis(P0,P2,Cell%V0,Cell%V2, &
                               Cell%V0%Point%y,Cell%V2%Point%y)
    CALL SearchPointsEdgeYaxis(P1,P3,Cell%V1,Cell%V3, &
                               Cell%V1%Point%y,Cell%V3%Point%y)
    CALL SearchPointsEdgeYaxis(P4,P6,Cell%V4,Cell%V6, &
                               Cell%V4%Point%y,Cell%V6%Point%y)
    CALL SearchPointsEdgeYaxis(P5,P7,Cell%V5,Cell%V7, &
                               Cell%V5%Point%y,Cell%V7%Point%y)

  IF(P0%y/=P2%y.AND.P0%y/=Cell%V0%Point%y) THEN
    PMin%y=MIN(PMin%y,P0%y)
    PMax%y=MAX(PMax%y,P0%y)
  END IF
  IF(P0%y/=P2%y.AND.P2%y/=Cell%V2%Point%y) THEN
    PMin%y=MIN(PMin%y,P2%y)
    PMax%y=MAX(PMax%y,P2%y)
  END IF

  IF(P1%y/=P3%y.AND.P1%y/=Cell%V1%Point%y) THEN
    PMin%y=MIN(PMin%y,P1%y)
    PMax%y=MAX(PMax%y,P1%y)
  END IF
  IF(P1%y/=P3%y.AND.P3%y/=Cell%V3%Point%y) THEN
    PMin%y=MIN(PMin%y,P3%y)
    PMax%y=MAX(PMax%y,P3%y)
  END IF

  IF(P4%y/=P6%y.AND.P4%y/=Cell%V4%Point%y) THEN
    PMin%y=MIN(PMin%y,P4%y)
    PMax%y=MAX(PMax%y,P4%y)
  END IF
  IF(P4%y/=P6%y.AND.P6%y/=Cell%V6%Point%y) THEN
    PMin%y=MIN(PMin%y,P6%y)
    PMax%y=MAX(PMax%y,P6%y)
  END IF

  IF(P5%y/=P7%y.AND.P5%y/=Cell%V5%Point%y) THEN
    PMin%y=MIN(PMin%y,P5%y)
    PMax%y=MAX(PMax%y,P5%y)
  END IF
  IF(P5%y/=P7%y.AND.P7%y/=Cell%V7%Point%y) THEN
    PMin%y=MIN(PMin%y,P7%y)
    PMax%y=MAX(PMax%y,P7%y)
  END IF

    ! EdgeZ
    CALL SearchPointsEdgeZaxis(P0,P4,Cell%V0,Cell%V4, &
                               Cell%V0%Point%z,Cell%V4%Point%z)
    CALL SearchPointsEdgeZaxis(P1,P5,Cell%V1,Cell%V5, &
                               Cell%V1%Point%z,Cell%V5%Point%z)
    CALL SearchPointsEdgeZaxis(P2,P6,Cell%V2,Cell%V6, &
                               Cell%V2%Point%z,Cell%V6%Point%z)
    CALL SearchPointsEdgeZaxis(P3,P7,Cell%V3,Cell%V7, &
                               Cell%V3%Point%z,Cell%V7%Point%z)
  IF(P0%z/=P4%z.AND.P0%z/=Cell%V0%Point%z) THEN
    PMin%z=MIN(PMin%z,P0%z)
    PMax%z=MAX(PMax%z,P0%z)
  END IF
  IF(P0%z/=P4%z.AND.P4%z/=Cell%V4%Point%z) THEN
    PMin%z=MIN(PMin%z,P4%z)
    PMax%z=MAX(PMax%z,P4%z)
  END IF

  IF(P1%z/=P5%z.AND.P1%z/=Cell%V1%Point%z) THEN
    PMin%z=MIN(PMin%z,P1%z)
    PMax%z=MAX(PMax%z,P1%z)
  END IF
  IF(P1%z/=P5%z.AND.P5%z/=Cell%V5%Point%z) THEN
    PMin%z=MIN(PMin%z,P5%z)
    PMax%z=MAX(PMax%z,P5%z)
  END IF

  IF(P2%z/=P6%z.AND.P2%z/=Cell%V2%Point%z) THEN
    PMin%z=MIN(PMin%z,P2%z)
    PMax%z=MAX(PMax%z,P2%z)
  END IF
  IF(P2%z/=P6%z.AND.P6%z/=Cell%V6%Point%z) THEN
    PMin%z=MIN(PMin%z,P6%z)
    PMax%z=MAX(PMax%z,P6%z)
  END IF

  IF(P3%z/=P7%z.AND.P3%z/=Cell%V3%Point%z) THEN
    PMin%z=MIN(PMin%z,P3%z)
    PMax%z=MAX(PMax%z,P3%z)
  END IF
  IF(P3%z/=P7%z.AND.P7%z/=Cell%V7%Point%z) THEN
    PMin%z=MIN(PMin%z,P7%z)
    PMax%z=MAX(PMax%z,P7%z)
  END IF

END SUBROUTINE SearchMinMaxP8Cell


SUBROUTINE AnalyzeCellRand(Cell)
  TYPE(CellP8_T) :: Cell
  INTEGER :: ListCut(1:2,1:6)
  INTEGER :: nCut,nCutNeu,i,j,iv
  INTEGER :: ic,jc,kc,nxc,nyc,nzc
  INTEGER :: iTemp
  REAL(8) :: xc0,xc1,yc0,yc1,zc0,zc1,dx,dy,dz
  TYPE (Vertex_T) :: Vertex1,Vertex2
  TYPE(Point_T) :: EdgP1,EdgP2,P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(Point_T) :: PMin,PMax,SumMidP,MidPoint
  REAL(8) :: xvc0,xvc1,xvu0,xvu1,xvo0,xvo1
  REAL(8) :: yvc0,yvc1,yvu0,yvu1,yvo0,yvo1
  REAL(8) :: zvc0,zvc1
  REAL(8) :: Vol,Len
  INTEGER :: VertexFace(1:6,1:9)
  LOGICAL :: Cut
  INTEGER :: iP1,iP2,jP1,jP2,iVert,jVert

  !Face 1-6 statisch
  Cell%Face1%in_out=Cell%V0%in_out+Cell%V1%in_out+Cell%V2%in_out+Cell%V3%in_out
  Cell%Face2%in_out=Cell%V4%in_out+Cell%V5%in_out+Cell%V6%in_out+Cell%V7%in_out
  Cell%Face3%in_out=Cell%V0%in_out+Cell%V1%in_out+Cell%V4%in_out+Cell%V5%in_out
  Cell%Face4%in_out=Cell%V2%in_out+Cell%V3%in_out+Cell%V6%in_out+Cell%V7%in_out
  Cell%Face5%in_out=Cell%V0%in_out+Cell%V2%in_out+Cell%V4%in_out+Cell%V6%in_out
  Cell%Face6%in_out=Cell%V1%in_out+Cell%V3%in_out+Cell%V5%in_out+Cell%V7%in_out
  Cell%in_out=Cell%Face1%in_out+Cell%Face2%in_out
  Cell%in_out=Cell%V0%in_out+Cell%V1%in_out+Cell%V2%in_out+Cell%V3%in_out &
              +Cell%V4%in_out+Cell%V5%in_out+Cell%V6%in_out+Cell%V7%in_out

  ! Volumen-Berechnung
  !IF (Cell%in_out<7.AND.Cell%in_out>-8) THEN
  !IF ((Cell%in_out<6.AND.Cell%in_out>-8.AND.Cell%Vol>0.0d0) .OR. 
  IF ((Cell%in_out<6.AND.Cell%in_out>-8) .OR. &
       (Cell%in_out==6.AND.Cell%vc>0) ) THEN
    IF (Cell%V4%in_out==0.AND.Cell%V5%in_out==0.AND.Cell%V6%in_out==0.AND.Cell%V7%in_out==0) THEN
       !Special Face2 Grenzflaeche-Orography
       Cell%Vol=0.0d0 
       Cell%mp=0
       Cell%vc=0
       !NrRN_Cells=NrRN_Cells+1  ! jeweils unter Zählung unterhalb Berg
    ELSE IF (Cell%V0%in_out==0.AND.Cell%V1%in_out==0.AND.Cell%V2%in_out==0.AND.Cell%V3%in_out==0.AND. &
             Cell%V4%in_out==1.AND.Cell%V5%in_out==1.AND.Cell%V6%in_out==1.AND.Cell%V7%in_out==1) THEN
       !Special Face1 Grenzflaeche-Orography 
       !IF (Cell%V0%Point%z== Gitter z0 )
       !IF (Cell%V0%Point%z==0.0d0 ) THEN
          P0=PointParametricEarth(Cell%V0%Point)
          P1=PointParametricEarth(Cell%V1%Point)
          P2=PointParametricEarth(Cell%V2%Point)
          P3=PointParametricEarth(Cell%V3%Point)
          P4=PointParametricEarth(Cell%V4%Point)
          P5=PointParametricEarth(Cell%V5%Point)
          P6=PointParametricEarth(Cell%V6%Point)
          P7=PointParametricEarth(Cell%V7%Point)
          Vol=VolCellP8(Cell)
          Cell%Vol=Vol
          Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
          NrMP_Cells=NrMP_Cells+1
          Cell%mp=NrMP_Cells
          NrR_Cells=NrR_Cells+1
          Cell%vc=4
       !ELSE
       !   Cell%vc=0
       !END IF
    ELSE IF ((Cell%in_out==6.OR.Cell%in_out==-6).AND. &
             (Cell%V0%in_out==0.OR.Cell%V1%in_out==0.OR.Cell%V2%in_out==0.OR.Cell%V3%in_out==0.OR. &
              Cell%V4%in_out==0.OR.Cell%V5%in_out==0.OR.Cell%V6%in_out==0.OR.Cell%V7%in_out==0)) THEN
       !Special Edge, Grenzflaeche-Orography
       !Vol-Berechnung in beiden Fällen (6,-6)kann entfalllen
       P0=PointParametricEarth(Cell%V0%Point)
       P1=PointParametricEarth(Cell%V1%Point)
       P2=PointParametricEarth(Cell%V2%Point)
       P3=PointParametricEarth(Cell%V3%Point)
       P4=PointParametricEarth(Cell%V4%Point)
       P5=PointParametricEarth(Cell%V5%Point)
       P6=PointParametricEarth(Cell%V6%Point)
       P7=PointParametricEarth(Cell%V7%Point)
       Vol=VolCellP8(Cell)
       Cell%Vol=Vol
       Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
       Cell%vc=0
       Cell%mp=0
       !IF(Cell%in_out==-6) THEN
       !  NrRN_Cells=NrRN_Cells+1
       !END IF
    ELSE IF ( Cell%in_out==2 .and. &
              (Cell%V0%in_out==0.AND.Cell%V1%in_out==0.AND.Cell%V2%in_out==0.AND.Cell%V3%in_out==0).and. &
              (Cell%V4%in_out==0.OR.Cell%V5%in_out==0.OR.Cell%V6%in_out==0.OR.Cell%V7%in_out==0))then
          ! Special (Face1 und (F3.or.F4.or.F5.or.F6))-Grenze
          P0=PointParametricEarth(Cell%V0%Point)
          P1=PointParametricEarth(Cell%V1%Point)
          P2=PointParametricEarth(Cell%V2%Point)
          P3=PointParametricEarth(Cell%V3%Point)
          P4=PointParametricEarth(Cell%V4%Point)
          P5=PointParametricEarth(Cell%V5%Point)
          P6=PointParametricEarth(Cell%V6%Point)
          P7=PointParametricEarth(Cell%V7%Point)
          Vol=VolCellP8(Cell)
          Cell%Vol=Vol
          Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
          NrMP_Cells=NrMP_Cells+1
          Cell%mp=NrMP_Cells
          NrR_Cells=NrR_Cells+1
          Cell%vc=6
    ELSE IF (Cell%in_out==4 .and. &
             ((Cell%V0%in_out==0.AND.Cell%V1%in_out==0.AND.Cell%V4%in_out==0.AND.Cell%V5%in_out==0).OR. &
              (Cell%V1%in_out==0.AND.Cell%V3%in_out==0.AND.Cell%V5%in_out==0.AND.Cell%V7%in_out==0).OR. &
              (Cell%V3%in_out==0.AND.Cell%V2%in_out==0.AND.Cell%V7%in_out==0.AND.Cell%V6%in_out==0).OR. &
              (Cell%V2%in_out==0.AND.Cell%V0%in_out==0.AND.Cell%V6%in_out==0.AND.Cell%V4%in_out==0))) THEN
          ! Special (F3.or.F4.or.F5.or.F6) - Grenze
          P0=PointParametricEarth(Cell%V0%Point)
          P1=PointParametricEarth(Cell%V1%Point)
          P2=PointParametricEarth(Cell%V2%Point)
          P3=PointParametricEarth(Cell%V3%Point)
          P4=PointParametricEarth(Cell%V4%Point)
          P5=PointParametricEarth(Cell%V5%Point)
          P6=PointParametricEarth(Cell%V6%Point)
          P7=PointParametricEarth(Cell%V7%Point)
          Vol=VolCellP8(Cell)
          Cell%Vol=Vol
          Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
          NrMP_Cells=NrMP_Cells+1
          Cell%mp=NrMP_Cells
          NrR_Cells=NrR_Cells+1
          Cell%vc=4
    ELSE
       !Celle mit Schnitt 
       Cell%Vol=0.0d0
       xc0=Cell%V0%Point%x  !xc0=Cell%Face1%Edge1%Vert1%Point%x
       yc0=Cell%V0%Point%y  !yc0=Cell%Face1%Edge1%Vert1%Point%y
       zc0=Cell%V0%Point%z  !zc0=Cell%Face1%Edge1%Vert1%Point%z
       xc1=Cell%V7%Point%x  !xc1=Cell%Face2%Edge2%Vert2%Point%x
       yc1=Cell%V7%Point%y  !yc1=Cell%Face2%Edge2%Vert2%Point%y
       zc1=Cell%V7%Point%z  !zc1=Cell%Face2%Edge2%Vert2%Point%z
       IF(xc0<xc1) THEN
         PMin=Cell%V7%Point   !PMin=Cell%Face2%Edge3%Vert2%Point
         PMax=Cell%V0%Point   !PMax=Cell%Face1%Edge1%Vert1%Point
       ELSE
         PMin=Cell%V6%Point
         PMax=Cell%V1%Point
       END IF
       !Routine mit Point-Check verwenden
       CALL SearchMinMaxP8Cell(Cell,PMin,PMax)
   
       IF((Cell%face5%in_out==4.AND.Cell%face6%in_out==-4).OR. &
          (Cell%face5%in_out==-4.AND.Cell%face6%in_out==4)) THEN
          !...Plane x-direction.................................
          zvc0=PMin%z+Shrink*ABS(PMin%z)
          zvc1=PMax%z-Shrink*ABS(PMax%z)
          yvc0=PMin%y+Shrink*ABS(PMin%y)
          yvc1=PMax%y-Shrink*ABS(PMax%y)
          yc0=yvc0
          !...Vergleich zu InitCellVol mit nz=1
          !...IncrVol 1Stelle nach Komma gleich, +50 (2St.), +150 (3St.),+200 (4St.)
          nxc=IncrVol ! +50 , +200
          nyc=IncrVol ! +50 , +200
          nzc=IncrVol ! +50 , +200
          dy=(yvc1-yvc0)/nyc
          dz=(zvc1-zvc0)/nzc
          Vertex1%Point%x=xc0
          Vertex1%in_out=-1
          Vertex2%Point%x=xc1
          Vertex2%in_out=1
          SumMidP%x=0.0d0
          SumMidP%y=0.0d0
          SumMidP%z=0.0d0
          DO kc=1,nzc
            yvc0=yc0
            DO jc=1,nyc
              !      2--------------6    !Point-Folge Vol-Berechnung
              !     /|             /|
              !    0--------------4 |
              !    | 3------------|-7
              !    |/             |/
              !    1--------------5
              Vertex1%Point%z=zvc0
              Vertex1%Point%y=yvc0
              Vertex2%Point%z=Vertex1%Point%z
              Vertex2%Point%y=Vertex1%Point%y
              CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
              P1=PointParametricEarth(EdgP1)
              P5=PointParametricEarth(EdgP2)
      
              Vertex1%Point%z=MIN(zvc0+dz,zvc1)
              Vertex1%Point%y=yvc0
              Vertex2%Point%z=Vertex1%Point%z
              Vertex2%Point%y=Vertex1%Point%y
              CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
              P0=PointParametricEarth(EdgP1)
              P4=PointParametricEarth(EdgP2)
      
              Vertex1%Point%z=zvc0
              Vertex1%Point%y=MIN(yvc0+dy,yvc1)
              Vertex2%Point%z=Vertex1%Point%z
              Vertex2%Point%y=Vertex1%Point%y
              CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
              P3=PointParametricEarth(EdgP1)
              P7=PointParametricEarth(EdgP2)
      
              Vertex1%Point%z=MIN(zvc0+dz,zvc1)
              Vertex1%Point%y=MIN(yvc0+dy,yvc1)
              Vertex2%Point%z=Vertex1%Point%z
              Vertex2%Point%y=Vertex1%Point%y
              CALL SearchPointsEdgeXaxis(EdgP1,EdgP2,Vertex1,Vertex2,xc0,xc1)
              P2=PointParametricEarth(EdgP1)
              P6=PointParametricEarth(EdgP2)
      
              Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
              !auch Vol=ABS(VolP8(P1,P5,P3,P7,P0,P4,P2,P6))
              Cell%Vol=Cell%Vol+Vol
              MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
              SumMidP=Vol*MidPoint+SumMidP

              yvc0=yvc0+dy
            END DO
            zvc0=zvc0+dz
          END DO
       ELSE IF((Cell%face3%in_out==4.AND.Cell%face4%in_out==-4).OR. &
               (Cell%face3%in_out==-4.AND.Cell%face4%in_out==4)) THEN
          !...Plane y-direction..................................
          zvc0=PMin%z+Shrink*ABS(PMin%z)
          zvc1=PMax%z-Shrink*ABS(PMax%z)
          xvc0=PMin%x+Shrink*ABS(PMin%x)
          xvc1=PMax%x-Shrink*ABS(PMax%x)
          xc0=xvc0
          !...Vergleich zu InitCellVol mit nz=1
          !...IncrVol 1Stelle nach Komma gleich, +50 (2St.), +150 (3St.),+200 (4St.)
          nxc=IncrVol ! +50 , +200
          nyc=IncrVol ! +50 , +200
          nzc=IncrVol ! +50 , +200
          dx=(xvc1-xvc0)/nxc
          dz=(zvc1-zvc0)/nzc
          Vertex1%Point%y=yc0
          Vertex1%in_out=-1
          Vertex2%Point%y=yc1
          Vertex2%in_out=1
          SumMidP%x=0.0d0
          SumMidP%y=0.0d0
          SumMidP%z=0.0d0
          DO kc=1,nzc
            xvc0=xc0
            DO ic=1,nxc
              !          6----7    !Point-Folge Vol-Berechnung
              !         /|   /|
              !        / |  / |
              !       /  2 /  3 
              !      4----5  /
              !      | /  | / 
              !      |/   |/ 
              !      0----1
              Vertex1%Point%z=zvc0
              Vertex1%Point%x=xvc0
              Vertex2%Point%z=Vertex1%Point%z
              Vertex2%Point%x=Vertex1%Point%x
              CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
              P0=PointParametricEarth(EdgP1)
              P2=PointParametricEarth(EdgP2)
      
              Vertex1%Point%z=zvc0
              Vertex1%Point%x=MIN(xvc0+dx,xvc1)
              Vertex2%Point%z=Vertex1%Point%z
              Vertex2%Point%x=Vertex1%Point%x
              CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
              P1=PointParametricEarth(EdgP1)
              P3=PointParametricEarth(EdgP2)
      
              Vertex1%Point%z=MIN(zvc0+dz,zvc1)
              Vertex1%Point%x=xvc0
              Vertex2%Point%z=Vertex1%Point%z
              Vertex2%Point%x=Vertex1%Point%x
              CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
              P4=PointParametricEarth(EdgP1)
              P6=PointParametricEarth(EdgP2)
      
              Vertex1%Point%z=MIN(zvc0+dz,zvc1)
              Vertex1%Point%x=MIN(xvc0+dx,xvc1)
              Vertex2%Point%z=Vertex1%Point%z
              Vertex2%Point%x=Vertex1%Point%x
              CALL SearchPointsEdgeYaxis(EdgP1,EdgP2,Vertex1,Vertex2,yc0,yc1)
              P5=PointParametricEarth(EdgP1)
              P7=PointParametricEarth(EdgP2)
      
              Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
              Cell%Vol=Cell%Vol+Vol
              MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
              SumMidP=Vol*MidPoint+SumMidP
  
              xvc0=xvc0+dx
            END DO
            zvc0=zvc0+dz
          END DO
       ELSE 
          !...Plane z-direction..................................
          xvc0=PMin%x+Shrink*ABS(PMin%x)
          xvc1=PMax%x-Shrink*ABS(PMax%x)
          yvc0=PMin%y+Shrink*ABS(PMin%y)
          yvc1=PMax%y-Shrink*ABS(PMax%y)
          yc0=yvc0
          !...Vergleich zu InitCellVol mit nz=1
          !...IncrVol 1Stelle nach Komma gleich, +50 (2St.), +150 (3St.),+200 (4St.)
          nxc=IncrVol ! +50 , +200
          nyc=IncrVol ! +50 , +200
          nzc=IncrVol ! +50 , +200
          dx=(xvc1-xvc0)/nxc
          dy=(yvc1-yvc0)/nyc
          Vertex1%Point%z=zc0
          Vertex1%in_out=-1
          Vertex2%Point%z=zc1
          Vertex2%in_out=1
          SumMidP%x=0.0d0
          SumMidP%y=0.0d0
          SumMidP%z=0.0d0
          DO ic=1,nxc
            yvc0=yc0
            DO jc=1,nyc
              !      6--------7    !Point-Folge Vol-Berechnung
              !     /|       /|
              !    4--------5 |
              !    | |      | |
              !    | 2------|-3
              !    |/       |/
              !    0--------1
              Vertex1%Point%x=xvc0
              Vertex1%Point%y=yvc0
              Vertex2%Point%x=Vertex1%Point%x
              Vertex2%Point%y=Vertex1%Point%y
              CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
              P0=PointParametricEarth(EdgP1)
              P4=PointParametricEarth(EdgP2)

              Vertex1%Point%x=MIN(xvc0+dx,xvc1)
              Vertex1%Point%y=yvc0
              Vertex2%Point%x=Vertex1%Point%x
              Vertex2%Point%y=Vertex1%Point%y
              CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
              P1=PointParametricEarth(EdgP1)
              P5=PointParametricEarth(EdgP2)
      
              Vertex1%Point%x=xvc0
              Vertex1%Point%y=MIN(yvc0+dy,yvc1)
              Vertex2%Point%x=Vertex1%Point%x
              Vertex2%Point%y=Vertex1%Point%y
              CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
              P2=PointParametricEarth(EdgP1)
              P6=PointParametricEarth(EdgP2)
      
              Vertex1%Point%x=MIN(xvc0+dx,xvc1)
              Vertex1%Point%y=MIN(yvc0+dy,yvc1)
              Vertex2%Point%x=Vertex1%Point%x
              Vertex2%Point%y=Vertex1%Point%y
              CALL SearchPointsEdgeZaxis(EdgP1,EdgP2,Vertex1,Vertex2,zc0,zc1)
              P3=PointParametricEarth(EdgP1)
              P7=PointParametricEarth(EdgP2)
              
              Vol=ABS(VolP8(P0,P1,P2,P3,P4,P5,P6,P7))
               
              Cell%Vol=Cell%Vol+Vol
              MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
              SumMidP=Vol*MidPoint+SumMidP

              yvc0=yvc0+dy
            END DO
            xvc0=xvc0+dx
          END DO
       END IF
       Cell%MidPoint=SumMidP/Cell%Vol 
       NrMP_Cells=NrMP_Cells+1
       Cell%mp=NrMP_Cells
       NrR_Cells=NrR_Cells+1
       Cell%vc=5     !vc noch ermitteln wenn Face erarbeitet
    END IF
  ELSE IF (Cell%in_out>=7) THEN ! Cell%vc==0  ?!
       !nur Point angrenzend oder keine Grenze
       If(Cell%in_out>=7) THEN
         Write(*,*) "Celle mit in_out=",Cell%in_out
         Write(*,*) "festgestellt -> in AnalyzeRandCelle -> IF anschauen!"
       END IF
       P0=PointParametricEarth(Cell%V0%Point)
       P1=PointParametricEarth(Cell%V1%Point)
       P2=PointParametricEarth(Cell%V2%Point)
       P3=PointParametricEarth(Cell%V3%Point)
       P4=PointParametricEarth(Cell%V4%Point)
       P5=PointParametricEarth(Cell%V5%Point)
       P6=PointParametricEarth(Cell%V6%Point)
       P7=PointParametricEarth(Cell%V7%Point)
       Vol=VolCellP8(Cell)
       Cell%Vol=Vol
       Cell%MidPoint=VolMidP8(P0,P1,P2,P3,P4,P5,P6,P7)
       !NrMP_Cells=NrMP_Cells+1
       !Cell%mp=NrMP_Cells
       !NrR_Cells=NrR_Cells+1
       Cell%vc=0
       Cell%mp=0
  END IF
 
END SUBROUTINE AnalyzeCellRand

SUBROUTINE EvalBordCell_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
  INTEGER :: ib,ix,vix,iy,iz,ibn,vjx,fac

  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(CellP8_T) :: CellRa

     P0=Floor(ib)%Vertices(vix,iy-1,iz-1)%Point
     P1=P0
     P1%x=P0%x+fac*Floor(ibn)%dx(vjx)
     P2=Floor(ib)%Vertices(vix,iy,iz-1)%Point
     P3=P2
     P3%x=P2%x+fac*Floor(ibn)%dx(vjx)
     !.........................................
     P4=Floor(ib)%Vertices(vix,iy-1,iz)%Point
     P5=P4
     P5%x=P4%x+fac*Floor(ibn)%dx(vjx)
     P6=Floor(ib)%Vertices(vix,iy,iz)%Point
     P7=P6
     P7%x=P6%x+fac*Floor(ibn)%dx(vjx)
     !.........................................
     !Vertex-Zuweisung äußerer Rand + Analyze
     !eventuell als Teil AnalyzeVertices 
     Floor(ib)%Vertices(vix+fac,iy-1,iz-1)%Point=P1
     Floor(ib)%Vertices(vix+fac,iy  ,iz-1)%Point=P3
     Floor(ib)%Vertices(vix+fac,iy-1,iz  )%Point=P5
     Floor(ib)%Vertices(vix+fac,iy  ,iz  )%Point=P7
     CALL CheckVertexRand(Floor(ib)%Vertices(vix+fac,iy-1,iz-1)) !V-P1
     CALL CheckVertexRand(Floor(ib)%Vertices(vix+fac,iy  ,iz-1)) !V-P3
     CALL CheckVertexRand(Floor(ib)%Vertices(vix+fac,iy-1,iz  )) !V-P5
     CALL CheckVertexRand(Floor(ib)%Vertices(vix+fac,iy  ,iz  )) !V-P7
     !.........................................
     !Zuweisung Celle-Rand für Übergabe
     !Vertex-Folge entsprechend Point-Folge für VolP8-Berechnung
     !Orientation east-west-direction or coordinate system is identical
     if(fac==-1) THEN !west-direction
       !'iw / Periode V1 '-------------------------
       CellRa%V0=Floor(ibn)%Vertices(vjx-1,iy-1,iz-1) 
       CellRa%V1=Floor(ibn)%Vertices(vjx  ,iy-1,iz-1)
       CellRa%V2=Floor(ibn)%Vertices(vjx-1,iy  ,iz-1)
       CellRa%V3=Floor(ibn)%Vertices(vjx  ,iy  ,iz-1)
       CellRa%V4=Floor(ibn)%Vertices(vjx-1,iy-1,iz) 
       CellRa%V5=Floor(ibn)%Vertices(vjx  ,iy-1,iz)
       CellRa%V6=Floor(ibn)%Vertices(vjx-1,iy  ,iz)
       CellRa%V7=Floor(ibn)%Vertices(vjx  ,iy  ,iz)
       !'iw / Periode V2-------------------------
       CellRa%V0=Floor(ibn)%Vertices(jx0-1,iy-1,iz-1)
       CellRa%V1=Floor(ibn)%Vertices(jx0  ,iy-1,iz-1)
       CellRa%V2=Floor(ibn)%Vertices(jx0-1,iy  ,iz-1)
       CellRa%V3=Floor(ibn)%Vertices(jx0  ,iy  ,iz-1)
       CellRa%V4=Floor(ibn)%Vertices(jx0-1,iy-1,iz)
       CellRa%V5=Floor(ibn)%Vertices(jx0  ,iy-1,iz)
       CellRa%V6=Floor(ibn)%Vertices(jx0-1,iy  ,iz)
       CellRa%V7=Floor(ibn)%Vertices(jx0  ,iy  ,iz)
       
     else if(fac==1) THEN !east-direction
       ! 'ie / Periode V1 '-------------------------
       CellRa%V0=Floor(ibn)%Vertices(vjx-1,iy-1,iz-1)
       CellRa%V1=Floor(ibn)%Vertices(vjx  ,iy-1,iz-1)
       CellRa%V2=Floor(ibn)%Vertices(vjx-1,iy  ,iz-1)
       CellRa%V3=Floor(ibn)%Vertices(vjx  ,iy  ,iz-1)
       CellRa%V4=Floor(ibn)%Vertices(vjx-1,iy-1,iz)
       CellRa%V5=Floor(ibn)%Vertices(vjx  ,iy-1,iz)
       CellRa%V6=Floor(ibn)%Vertices(vjx-1,iy  ,iz)
       CellRa%V7=Floor(ibn)%Vertices(vjx  ,iy  ,iz)
      !'ie / Periode V2 '-------------------------
       CellRa%V0=Floor(ibn)%Vertices(jx0  ,iy-1,iz-1)
       CellRa%V1=Floor(ibn)%Vertices(jx0+1,iy-1,iz-1)
       CellRa%V2=Floor(ibn)%Vertices(jx0  ,iy  ,iz-1)
       CellRa%V3=Floor(ibn)%Vertices(jx0+1,iy  ,iz-1)
       CellRa%V4=Floor(ibn)%Vertices(jx0  ,iy-1,iz) 
       CellRa%V5=Floor(ibn)%Vertices(jx0+1,iy-1,iz)
       CellRa%V6=Floor(ibn)%Vertices(jx0  ,iy  ,iz)
       CellRa%V7=Floor(ibn)%Vertices(jx0+1,iy  ,iz)
     end if
     !.........................................
     CALL AnalyzeCellRand(CellRa)
     IF(CellRa%in_out>-8.AND.CellRa%vc>0) THEN
       ALLOCATE(Floor(ib)%Cell(ix,iy,iz)%Cell)
       Floor(ib)%Cell(ix,iy,iz)%Cell%mp=CellRa%mp
       Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint=CellRa%MidPoint
       Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint%x= &
            Floor(ib)%Vertices(vix,iy,iz)%Point%x + &
             & fac*ABS(CellRa%MidPoint%x-Floor(ibn)%Vertices(jx0,iy,iz)%Point%x) 
       Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=CellRa%Vol
       Floor(ib)%Cell(ix,iy,iz)%Cell%in_out=CellRa%in_out
       Floor(ib)%Cell(ix,iy,iz)%Cell%vc=CellRa%vc !theoretisch/symbolisch für Weight-Ausgabe
       IF(MAX(CellRa%V0%in_out,CellRa%V1%in_out,CellRa%V2%in_out,CellRa%V3%in_out, &
              CellRa%V4%in_out,CellRa%V5%in_out,CellRa%V6%in_out,CellRa%V7%in_out)==0)THEN
          Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=0
          Floor(ib)%Cell(ix,iy,iz)%Cell%mp=0
          ! Cellen mit Grenze, 
          ! aber Special, da Schnitt sehr sehr klein und vernachlässigt werden soll,
          ! durch dist(x,y,z)_coeff ausgeblendet und Face komplett Grenze wird
       END IF 
     END IF
     !............................................
     !Uebernehme Koordinate der akt RandCelle -> Berechne FacesRand
     CALL EvalBorderFaceYZ_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
     !............................................
END SUBROUTINE EvalBordCell_WE

SUBROUTINE EvalBordCell_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
  INTEGER :: ib,ix,iy,viy,iz,ibn,vjy,fac

  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(CellP8_T) :: CellRa

     P0=Floor(ib)%Vertices(ix,viy,iz-1)%Point
     P1=P0
     P1%y=P0%y+fac*Floor(ibn)%dy(vjy)
     P2=Floor(ib)%Vertices(ix-1,viy,iz-1)%Point
     P3=P2
     P3%y=P2%y+fac*Floor(ibn)%dy(vjy)
     !.........................................
     P4=Floor(ib)%Vertices(ix,viy,iz)%Point
     P5=P4
     P5%y=P4%y+fac*Floor(ibn)%dy(vjy)
     P6=Floor(ib)%Vertices(ix-1,viy,iz)%Point
     P7=P6
     P7%y=P6%y+fac*Floor(ibn)%dy(vjy)
     !.........................................
     !Vertex-Zuweisung äußerer Rand + Analyze
     !eventuell als Teil AnalyzeVertices 
     Floor(ib)%Vertices(ix  ,viy+fac,iz-1)%Point=P1
     Floor(ib)%Vertices(ix-1,viy+fac,iz-1)%Point=P3
     Floor(ib)%Vertices(ix  ,viy+fac,iz  )%Point=P5
     Floor(ib)%Vertices(ix-1,viy+fac,iz  )%Point=P7
     CALL CheckVertexRand(Floor(ib)%Vertices(ix  ,viy+fac,iz-1)) !V-P1
     CALL CheckVertexRand(Floor(ib)%Vertices(ix-1,viy+fac,iz-1)) !V-P3
     CALL CheckVertexRand(Floor(ib)%Vertices(ix  ,viy+fac,iz))   !V-P5
     CALL CheckVertexRand(Floor(ib)%Vertices(ix-1,viy+fac,iz))   !V-P7
     !.........................................
     !Zuweisung Celle-Rand für Übergabe
     !Vertex-Folge entsprechend Point-Folge für VolP8-Berechnung
     if(fac==-1) THEN !south-direction
       !'is / Periode V1' ------ mit vjy- Aktualparameter ----
       CellRa%V0=Floor(ibn)%Vertices(ix-1,vjy+fac,iz-1)  !V-P3
       CellRa%V1=Floor(ibn)%Vertices(ix  ,vjy+fac,iz-1)  !V-P1
       CellRa%V2=Floor(ibn)%Vertices(ix-1,vjy    ,iz-1)  !V-P2
       CellRa%V3=Floor(ibn)%Vertices(ix  ,vjy    ,iz-1)  !V-P0
       CellRa%V4=Floor(ibn)%Vertices(ix-1,vjy+fac,iz)   !V-P7
       CellRa%V5=Floor(ibn)%Vertices(ix  ,vjy+fac,iz)   !V-P5
       CellRa%V6=Floor(ibn)%Vertices(ix-1,vjy    ,iz)   !V-P6
       CellRa%V7=Floor(ibn)%Vertices(ix  ,vjy    ,iz)   !V-P4
       !'is / Periode V2' -----mit jy0-Nacbar-def -----------
       CellRa%V0=Floor(ibn)%Vertices(ix-1,jy0+fac,iz-1)  !V-P3
       CellRa%V1=Floor(ibn)%Vertices(ix  ,jy0+fac,iz-1)  !V-P1
       CellRa%V2=Floor(ibn)%Vertices(ix-1,jy0    ,iz-1)  !V-P2
       CellRa%V3=Floor(ibn)%Vertices(ix  ,jy0    ,iz-1)  !V-P0
       CellRa%V4=Floor(ibn)%Vertices(ix-1,jy0+fac,iz)   !V-P7
       CellRa%V5=Floor(ibn)%Vertices(ix  ,jy0+fac,iz)   !V-P5
       CellRa%V6=Floor(ibn)%Vertices(ix-1,jy0    ,iz)   !V-P6
       CellRa%V7=Floor(ibn)%Vertices(ix  ,jy0    ,iz)   !V-P4
     else if(fac==+1) THEN !north direction
       !'in / Period V1'  mit vjy-Aktualparameter -----------
       CellRa%V0=Floor(ibn)%Vertices(ix-1,vjy-1,iz-1) !V-P2
       CellRa%V1=Floor(ibn)%Vertices(ix  ,vjy-1,iz-1) !V-P0
       CellRa%V2=Floor(ibn)%Vertices(ix-1,vjy  ,iz-1) !V-P3
       CellRa%V3=Floor(ibn)%Vertices(ix  ,vjy  ,iz-1) !V-P1
       CellRa%V4=Floor(ibn)%Vertices(ix-1,vjy-1,iz)    !V-P6
       CellRa%V5=Floor(ibn)%Vertices(ix  ,vjy-1,iz)    !V-P4
       CellRa%V6=Floor(ibn)%Vertices(ix-1,vjy    ,iz)    !V-P7
       CellRa%V7=Floor(ibn)%Vertices(ix  ,vjy    ,iz)    !V-P5
       !'in / Period V2'  mit jy0 - Nachbar-def --------------
       CellRa%V0=Floor(ibn)%Vertices(ix-1,jy0  ,iz-1) !V-P2
       CellRa%V1=Floor(ibn)%Vertices(ix  ,jy0  ,iz-1) !V-P0
       CellRa%V2=Floor(ibn)%Vertices(ix-1,jy0+1,iz-1) !V-P3
       CellRa%V3=Floor(ibn)%Vertices(ix  ,jy0+1,iz-1) !V-P1
       CellRa%V4=Floor(ibn)%Vertices(ix-1,jy0  ,iz)    !V-P6
       CellRa%V5=Floor(ibn)%Vertices(ix  ,jy0  ,iz)    !V-P4
       CellRa%V6=Floor(ibn)%Vertices(ix-1,jy0+1,iz)    !V-P7
       CellRa%V7=Floor(ibn)%Vertices(ix  ,jy0+1,iz)    !V-P5
     end if
     !.........................................
     CALL AnalyzeCellRand(CellRa)
     !IF(CellRa%in_out>-8.AND.CellRa%in_out<6) THEN  ?NrMP nein NrW ja?
     IF(CellRa%in_out>-8.AND.CellRa%vc>0) THEN
       ALLOCATE(Floor(ib)%Cell(ix,iy,iz)%Cell)
       Floor(ib)%Cell(ix,iy,iz)%Cell%mp=CellRa%mp
       Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint=CellRa%MidPoint
       Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint%y= &
            Floor(ib)%Vertices(ix,viy,iz)%Point%y + &
             & fac*ABS(CellRa%MidPoint%y-Floor(ibn)%Vertices(ix,jy0,iz)%Point%y) 
       Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=CellRa%Vol
       Floor(ib)%Cell(ix,iy,iz)%Cell%in_out=CellRa%in_out
       Floor(ib)%Cell(ix,iy,iz)%Cell%vc=CellRa%vc !theoretisch/symbolisch für Weight-Ausgabe
       IF(MAX(CellRa%V0%in_out,CellRa%V1%in_out,CellRa%V2%in_out,CellRa%V3%in_out, &
              CellRa%V4%in_out,CellRa%V5%in_out,CellRa%V6%in_out,CellRa%V7%in_out)==0)THEN
          Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=0
          Floor(ib)%Cell(ix,iy,iz)%Cell%mp=0
       END IF 
     END IF
     !............................................
     !!Uebernehme Koordinate der akt RandCelle Berechne FacesRand
     CALL EvalBorderFaceZX_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
     !............................................
END SUBROUTINE EvalBordCell_NS

SUBROUTINE EvalBordCell_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
  INTEGER :: ib,ix,iy,iz,viz,ibn,vjz,fac

  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(CellP8_T) :: CellRa

    P0=Floor(ib)%Vertices(ix,iy,viz)%Point
    P1=P0
    P1%z=P0%z+fac*Floor(ibn)%dz(vjz)
    P2=Floor(ib)%Vertices(ix-1,iy,viz)%Point
    P3=P2
    P3%z=P2%z+fac*Floor(ibn)%dz(vjz)
    !.........................................
    P4=Floor(ib)%Vertices(ix,iy-1,viz)%Point
    P5=P4
    P5%z=P4%z+fac*Floor(ibn)%dz(vjz)
    P6=Floor(ib)%Vertices(ix-1,iy-1,viz)%Point
    P7=P6
    P7%z=P6%z+fac*Floor(ibn)%dz(vjz)
    !.........................................
    !Vertex-Zuweisung äußerer Rand + Analyze
    !eventuell als Teil AnalyzeVertices 
    Floor(ib)%Vertices(ix  ,iy  ,viz+fac)%Point=P1
    Floor(ib)%Vertices(ix-1,iy  ,viz+fac)%Point=P3
    Floor(ib)%Vertices(ix  ,iy-1,viz+fac)%Point=P5
    Floor(ib)%Vertices(ix-1,iy-1,viz+fac)%Point=P7
    CALL CheckVertexRand(Floor(ib)%Vertices(ix  ,iy  ,viz+fac)) !V-P1
    CALL CheckVertexRand(Floor(ib)%Vertices(ix-1,iy  ,viz+fac)) !V-P3
    CALL CheckVertexRand(Floor(ib)%Vertices(ix  ,iy-1,viz+fac)) !V-P5
    CALL CheckVertexRand(Floor(ib)%Vertices(ix-1,iy-1,viz+fac)) !V-P7
    !.........................................
    !Zuweisung Celle-Rand für Übergabe
    !Vertex-Folge entsprechend Point-Folge für VolP8-Berechnung
    !1) Orientation top-botton-direction, these assignment is as well right
    if(fac==-1) THEN !botton-direction
     !'ib / Period V1'  mit vjz-Aktualparameter -----------
     !'ib / Period V2'  mit jz0 - Nachbar-def --------------
     CellRa%V0=Floor(ibn)%Vertices(ix-1,iy-1,jz0+fac) !V-P7
     CellRa%V1=Floor(ibn)%Vertices(ix  ,iy-1,jz0+fac) !V-P5
     CellRa%V2=Floor(ibn)%Vertices(ix-1,iy  ,jz0+fac) !V-P3
     CellRa%V3=Floor(ibn)%Vertices(ix  ,iy  ,jz0+fac) !V-P1
     CellRa%V4=Floor(ibn)%Vertices(ix-1,iy-1,jz0)     !V-P6
     CellRa%V5=Floor(ibn)%Vertices(ix  ,iy-1,jz0)     !V-P4
     CellRa%V6=Floor(ibn)%Vertices(ix-1,iy  ,jz0)     !V-P2
     CellRa%V7=Floor(ibn)%Vertices(ix  ,iy  ,jz0)     !V-P0
    else if(fac==+1) THEN !top-direction
     !'it / Period V1'  mit vjz-Aktualparameter -----------
     !'it / Period V2'  mit jz0 - Nachbar-def --------------
     CellRa%V0=Floor(ibn)%Vertices(ix-1,iy-1,jz0)     !V-P6
     CellRa%V1=Floor(ibn)%Vertices(ix  ,iy-1,jz0)     !V-P4
     CellRa%V2=Floor(ibn)%Vertices(ix-1,iy  ,jz0)     !V-P2
     CellRa%V3=Floor(ibn)%Vertices(ix  ,iy  ,jz0)     !V-P0
     CellRa%V4=Floor(ibn)%Vertices(ix-1,iy-1,jz0+fac) !V-P7
     CellRa%V5=Floor(ibn)%Vertices(ix  ,iy-1,jz0+fac) !V-P5
     CellRa%V6=Floor(ibn)%Vertices(ix-1,iy  ,jz0+fac) !V-P3
     CellRa%V7=Floor(ibn)%Vertices(ix  ,iy  ,jz0+fac) !V-P1
    end if
    !.........................................
    CALL AnalyzeCellRand(CellRa)
    !IF(CellRa%in_out>-8.AND.CellRa%in_out<6) THEN  ?NrMP nein NrW ja?
    IF(CellRa%in_out>-8.AND.CellRa%vc>0) THEN
      ALLOCATE(Floor(ib)%Cell(ix,iy,iz)%Cell)
      Floor(ib)%Cell(ix,iy,iz)%Cell%mp=CellRa%mp
      Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint=CellRa%MidPoint
      Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint%z= &
            Floor(ib)%Vertices(ix,iy,viz)%Point%z + &
             & fac*ABS(CellRa%MidPoint%z-Floor(ibn)%Vertices(ix,iy,jz0)%Point%z) 
      Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=CellRa%Vol
      Floor(ib)%Cell(ix,iy,iz)%Cell%in_out=CellRa%in_out
      Floor(ib)%Cell(ix,iy,iz)%Cell%vc=CellRa%vc !theoretisch/symbolisch für Weight-Ausgabe
      IF(MAX(CellRa%V0%in_out,CellRa%V1%in_out,CellRa%V2%in_out,CellRa%V3%in_out, &
             CellRa%V4%in_out,CellRa%V5%in_out,CellRa%V6%in_out,CellRa%V7%in_out)==0)THEN
         Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=0
         Floor(ib)%Cell(ix,iy,iz)%Cell%mp=0
      END IF 
    END IF
     !............................................
     !Uebernehme Koordinate der akt RandCelle -> Berechne FacesRand
     !!1 !getestet 26.11.2009,VGrid_GMV_8.4.1.2
     CALL EvalBorderFaceXY_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
     !!1 CALL EvalBorderFaceYZ_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
     !!1 CALL EvalBorderFaceZX_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
    !.........................................
END SUBROUTINE EvalBordCell_TB


SUBROUTINE AnalyzeAllCellsRand
  !RandCellen der NachbarBloecke ermitteln, Volumen bzw. MidPoint berechnen,
  !                                         oder zuweisen bei gleichen Vergröberungsgrad
  !------------------------------------------------------------------------------------
  INTEGER :: i,j,k,bli,blj,ix,iy,iz,jx,jy,jz,vix,vjx,viy,vjy,viz,vjz,fac,facf
  INTEGER :: px,py,pz
  INTEGER :: in,ib
  INTEGER :: widthz,widthy,widthx,in_out
  TYPE(Point_T) :: P0,P1,P2,P3,P4,P5,P6,P7
  TYPE(CellP8_T) :: CellRa
  TYPE(Point_T) :: SumMidP
  REAL(8) :: Vol

  
  DO ib=1,nb
    Write(*,*) "                                     Block : ",ib,"\/",nb 
    CALL Set(Floor(ib))
    NrR_Cells=0
    NrRN_Cells=0
    NrR_FacesXY=0;NrR_FacesYZ=0;NrR_FacesZX=0
    NrRN_FacesXY=0;NrRN_FacesYZ=0;NrRN_FacesZX=0

    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))

      SELECT CASE (Nachbars(in)%nType(2:2))

      !......................................................................................
      !...........................TypeW/TypeE................................................
      !......................................................................................
      CASE ("w","e")      !? pe,pw !? läüft  unter else wie 'ie,iw'
        IF (nType=='ow') THEN
          DO iy=iy0+1,iy1
             DO iz=iz0+1,iz1
                !DO WHILE(.NOT.ASSOCIATED(Cell(ix0+1,iy,iz)%Cell).AND.iz/=iz1)
                  !Write(*,*) "ib=",ib, " Type-ow:  ix0=",ix0," iy=",iy," iz=",i, &
                  !           "  RandCell unter Berg west"
                  !!ALLOCATE(Floor(ib)%Cell(ix0,iy,iz)%Cell)
                  !!Floor(ib)%Cell(ix0,iy,iz)%Cell%in_out=-8
                  !!Floor(ib)%Cell(ix0,iy,iz)%Cell%Vol=0.0d0
                  !NrRN_Cells=NrRN_Cells+1
                  !IF(iz==Domain%iz0) THEN  
                  !   NrRN_Cells=NrRN_Cells-1
                  !END IF
                 !  iz=iz+1
                !END DO 
                !DO WHILE(ASSOCIATED(Cell(ix0+1,iy,iz)%Cell).AND.iz/=iz1+1)
                IF (ASSOCIATED(Cell(ix0+1,iy,iz)%Cell).AND.iz/=iz1+1) THEN
                  !IF(Cell(ix0+1,iy,iz)%Cell%vc>0) THEN
                  IF(Cell(ix0+1,iy,iz)%Cell%in_out<6.AND.Cell(ix0+1,iy,iz)%Cell%Vol>0.0d0.OR. &
                     (Cell(ix0+1,iy,iz)%Cell%in_out==6.AND.Cell(ix0+1,iy,iz)%Cell%vc>0) ) THEN
                    ALLOCATE(Floor(ib)%Cell(ix0,iy,iz)%Cell)
                    Floor(ib)%Cell(ix0,iy,iz)%Cell%Vol=Floor(ib)%Cell(ix0+1,iy,iz)%Cell%Vol
                    Floor(ib)%Cell(ix0,iy,iz)%Cell%MidPoint= &
                          Floor(ib)%Cell(ix0+1,iy,iz)%Cell%MidPoint
                    Floor(ib)%Cell(ix0,iy,iz)%Cell%MidPoint%x= &
                          Floor(ib)%Cell(ix0,iy,iz)%Cell%MidPoint%x-dx(ix0+1)
                    Floor(ib)%Cell(ix0,iy,iz)%Cell%in_out=Floor(ib)%Cell(ix0+1,iy,iz)%Cell%in_out
                    Floor(ib)%Cell(ix0,iy,iz)%Cell%vc=Floor(ib)%Cell(ix0+1,iy,iz)%Cell%vc
                    NrR_Cells=NrR_Cells+1
                    IF(Floor(ib)%Cell(ix0,iy,iz)%Cell%Vol>0.0d0) THEN
                      !Cell mit maxVol nicht unter NrMP_Cells
                      IF(Floor(ib)%Cell(ix0+1,iy,iz)%Cell%mp>0) THEN  
                        !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Grenze-Berg sein können
                         !.OR. &
                         !(MIN(Cell(ix0+1,iy,iz)%Cell%Face1%Edge1%Vert1%in_out, &
                         !     Cell(ix0+1,iy,iz)%Cell%Face1%Edge1%Vert2%in_out, &
                         !     Cell(ix0+1,iy,iz)%Cell%Face1%Edge3%Vert1%in_out, &
                         !     Cell(ix0+1,iy,iz)%Cell%Face1%Edge3%Vert2%in_out)==0) ) THEN  !Face1-Grenzfläche
                        !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Grenze-Berg sein können
                        NrMP_Cells=NrMP_Cells+1
                        Floor(ib)%Cell(ix0,iy,iz)%Cell%mp=NrMP_Cells
                      ELSE
                        Floor(ib)%Cell(ix0,iy,iz)%Cell%mp=0
                      END IF
                    END IF
                  END IF
                  !Write(*,*) "ib=",ib, " Type-ow:  ix0=",ix0," iy=",iy," iz=",iz, " Schnitt RandCelle west"
                  !iz=iz+1
                !END DO 
                END IF !ASSOCIATED(Cell....
              END DO  !iz
          END DO !iy
          CALL EvalBorderAllFace_OutW(ib,ix0,iy0,iy1,iz0,iz1)
        ELSE IF(nType=="oe") THEN
          DO iy=iy0+1,iy1
            DO iz=iz0+1,iz1
               !DO WHILE(.NOT.ASSOCIATED(Cell(ix1,iy,iz)%Cell).AND.iz/=iz1)
                 !Write(*,*) "ib=",ib, " Type-oe:  ,ix1+1=",ix1+1, "iy=",iy,"iz=",i, &
                 !           "  RandCell unter Berg east"
                  !!ALLOCATE(Floor(ib)%Cell(ix1+1,iy,iz)%Cell)
                  !!Floor(ib)%Cell(ix1+1,iy,iz)%Cell%in_out=-8
                  !!Floor(ib)%Cell(ix1+1,iy,iz)%Cell%Vol=0.0d0
                  !NrRN_Cells=NrRN_Cells+1
                  !IF(iz==Domain%iz0) THEN
                  !   NrRN_Cells=NrRN_Cells-1
                  !END IF
                  !iz=iz+1
               !END DO 
             !DO While(ASSOCIATED(Cell(ix1,iy,iz)%Cell).AND.iz/=iz1+1)
             IF(ASSOCIATED(Cell(ix1,iy,iz)%Cell).AND.iz/=iz1+1) THEN
               !IF(Cell(ix1,iy,iz)%Cell%vc>0) THEN reicht nicht + Cellen Face1-Grenzfläche
               IF(Cell(ix1,iy,iz)%Cell%in_out<6.AND.Cell(ix1,iy,iz)%Cell%Vol>0.0d0.OR. &
                  (Cell(ix1,iy,iz)%Cell%in_out==6.AND.Cell(ix1,iy,iz)%Cell%vc>0) ) THEN
                 ALLOCATE(Floor(ib)%Cell(ix1+1,iy,iz)%Cell)
                 Floor(ib)%Cell(ix1+1,iy,iz)%Cell%Vol=Floor(ib)%Cell(ix1,iy,iz)%Cell%Vol
                 Floor(ib)%Cell(ix1+1,iy,iz)%Cell%MidPoint= &
                       Floor(ib)%Cell(ix1,iy,iz)%Cell%MidPoint
                 Floor(ib)%Cell(ix1+1,iy,iz)%Cell%MidPoint%x= &
                       Floor(ib)%Cell(ix1,iy,iz)%Cell%MidPoint%x+dx(ix1)
                 Floor(ib)%Cell(ix1+1,iy,iz)%Cell%in_out=Floor(ib)%Cell(ix1,iy,iz)%Cell%in_out
                 Floor(ib)%Cell(ix1+1,iy,iz)%Cell%vc=Floor(ib)%Cell(ix1,iy,iz)%Cell%vc
                 NrR_Cells=NrR_Cells+1
                 IF(Floor(ib)%Cell(ix1+1,iy,iz)%Cell%Vol>0.0d0) THEN
                   IF(Floor(ib)%Cell(ix1,iy,iz)%Cell%mp>0 ) THEN 
                     !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Berg sein können
                     NrMP_Cells=NrMP_Cells+1
                     Floor(ib)%Cell(ix1+1,iy,iz)%Cell%mp=NrMP_Cells
                   ELSE
                     Floor(ib)%Cell(ix1+1,iy,iz)%Cell%mp=0
                   END IF
                 END IF
               END IF
               !Write(*,*) "ib=",ib, " Type-oe:  ix1+1=",ix1+1,"iy=",iy, "iz=",iz, &
               !           " Schnitt RandCelle east"
               !iz=iz+1
             !END DO
             END IF
            END DO !iz
          END DO !iy
          CALL EvalBorderAllFace_OutE(ib,ix1,iy0,iy1,iz0,iz1)
        !...oder 
        !SELECT CASE (Nachbars(in)%nType(1:1))
        !CASE ("o")
        !  IF (nType(2:2)=='w') THEN
        !    ix=ix0+1
        !    pos=-1
        !  ELSE !!"e"
        !    ix=ix1
        !    pos=+1
        !  END IF
        !  DO iz=iz0+1,iz1
        !    DO iy=iy0+1,iy1
        !      IF (ASSOCIATED(Cell(ix,iy,iz)%Cell)) THEN
        !        Floor(ib)%Cell(ix+pos,iy,iz)%Cell=>Floor(ib)%Cell(ix,iy,iz)%Cell
        !      END IF
        !    END DO
        !  END DO
        !CASE ("i")
        !   .......
        !END SELECT  !!'o','i' fuer 'w','e'

        ELSE  !!'iw'.OR.'ie' 
          IF (nType(2:2)=='w') THEN
            ix=ix0               !Plane-akt.RandCells
            jx=Floor(ibn)%ix1    !Plane-Nachbar-Cells
            dx(ix0)=Floor(ibn)%dx(jx)
            !........
            vix=ix0              !Plane-Point-P0 Celle
            vjx=Floor(ibn)%ix1   !Pos. dx-Nachbar
            fac=-1               !Faktor Position Point (P1,P3,P5,P7)
            facf=1
          ELSE !"e"
            ix=ix1+1             !Plane-akt.RandCells
            jx=Floor(ibn)%ix0+1  !Plane-Nachbar-Cells
            dx(ix1+1)=Floor(ibn)%dx(jx)  
            !........ 
            vix=ix1              !Plane-Point-P0 Celle
            vjx=Floor(ibn)%ix0+1 !Pos. dx-Nachbar
            fac=1                !Faktor Position Point (P1,P3,P5,P7)
            facf=0
          END IF

          !!!----
          IF (RefineZ>RefineNachbarZ) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!     -------------           -------------
          !!     |     |  |  |  'iw,ie'  |  |  |     |
          !!     |     |-----|           |-----|     |
          !!     |     |  |  |           |  |  |     |
          !!     -------------           -------------

            IF (RefineY>=RefineNachbarY) THEN
            !-------------------------------
            !!    ...7------6---....-----6------7...     !Point-Folge P0-P7 'w'/'e' allg.
            !!      /|     /|           /|     /|  
            !!  ...5------4-----....---4------5...
            !!    .|.3....|.2..........|.2....|.3...
            !!     |/ 'iw'|/   (ib)    |/ 'ie'|/
            !!  ...1------0----....----0------1...
            !!   (ibn)   /   (ib)     /  (ibn) 

              !Floor(ibn)%jy0,Floor(ibn)%jy1 
              !Floor(ibn)%jz0,Floor(ibn)%jz1
              !jx0,jx1,jy0,jy1,jz0,jz1: wird über SetNachbar gesetzt 
              !                         und def. GrenzFläche Akt-Block für Nachbar-Block
              DO jy=jy0+1,jy1
                jz=jz0+1
                DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND.jz/=jz1)
                   jz=jz+1
                END DO
                DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND.jz/=jz1+1)
                  iy=IncrY*jy
                  iz=IncrZ*jz
                  !............................
                  piz: DO pz=1,IncrZ
                    piy: DO py=1,IncrY
                      CALL EvalBordCell_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
                      iy=iy-1
                    END DO piy
                    iz=iz-1
                    iy=IncrY*jy
                  END DO piz
                  !............................
                  jz=jz+1
                END DO  !DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell...
              END DO   !jy, NachbarBlock-Cell             
            ELSE IF (RefineY<RefineNachbarY) THEN
            !-------------------------------
              DO jy=jy0+1,jy1,2
                jz=jz0+1
                DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND. &
                         .NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy+1,jz)%Cell).AND. &
                         jz/=jz1)
                   jz=jz+1
                END DO
                DO WHILE((ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                          ASSOCIATED(Floor(ibn)%Cell(jx,jy+1,jz)%Cell)) .AND. &
                          jz/=jz1+1)
                   iy=(jy+1)/IncrY
                   iz=IncrZ*jz
                   !............................
                   DO pz=1,IncrZ
                      CALL EvalBordCell_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
                      iz=iz-1
                   END DO
                   !............................
                   jz=jz+1
                END DO  !DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell... .OR.
                        !         ASSOCIATED(Floor(ibn)%Cell(jx,jy+1,jz)%Cell...) ...
              END DO  !jy, NachbarBlock-Cell             
            END IF  !(RefineY[>;==;<]RefineNachbarY)

          !!!----
          ELSE IF (RefineZ==RefineNachbarZ) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!      -------------
          !!      |     |     |   'iw,ie'
          !!      |---------- |
          !!      |     |     |
          !!      -------------
            !IF (RefineY>=RefineNachbarY) THEN
            IF (RefineY>RefineNachbarY) THEN
            !-------------------------------
              DO jy=jy0+1,jy1
                jz=jz0+1
                DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND.jz/=jz1)
                   jz=jz+1
                END DO
                DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND.jz/=jz1+1)
                  iy=IncrY*jy
                  iz=IncrZ*jz
                  !............................
                  DO py=1,IncrY
                    CALL EvalBordCell_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
                    iy=iy-1
                  END DO
                  !............................
                  jz=jz+1
                END DO  !DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell...
              END DO   !jy, NachbarBlock-Cell             
            ELSE IF (RefineY==RefineNachbarY) THEN
            !-------------------------------------
              DO iy=iy0+1,iy1
                DO iz=iz0+1,iz1
                 ! IF (ASSOCIATED(Floor(ib)%Cell(jx,iy,iz)%Cell)) THEN
                 !   !................................................. 
                 !   CALL EvalBordCell_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
                 !   !................................................. 
                 ! END IF
                  !Vol=0.0d0;SumMidP%x=0.0d0;SumMidP%y=0.0d0;SumMidP%z=0.0d0
                  IF (ASSOCIATED(Floor(ibn)%Cell(jx,iy,iz)%Cell)) THEN
                    ALLOCATE(Floor(ib)%Cell(ix,iy,iz)%Cell)
                    Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=Floor(ibn)%Cell(jx,iy,iz)%Cell%Vol
                    Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint= &
                                                 Floor(ibn)%Cell(jx,iy,iz)%Cell%MidPoint
                    Floor(ib)%Cell(ix,iy,iz)%Cell%in_out= &
                                                 Floor(ibn)%Cell(jx,iy,iz)%Cell%in_out
                    Floor(ib)%Cell(ix,iy,iz)%Cell%vc= &
                                                 Floor(ibn)%Cell(jx,iy,iz)%Cell%vc
                    !Floor(ib)%Cell(ix,iy,iz)%Cell=Floor(ibn)%Cell(jx,iy,iz)%Cell
    Write(OutUnitProt,*) "CellRandVergl iw/ew ibn,ib,jx,iy,iz", ibn,ib,jx,iy,iz," Celle gl. Refine"
                    IF(Floor(ibn)%Cell(jx,iy,iz)%Cell%Vol>0.0d0) THEN
                      NrR_Cells=NrR_Cells+1
                      IF(Floor(ibn)%Cell(jx,iy,iz)%Cell%vc>0) THEN
                        IF(Floor(ibn)%Cell(jx,iy,iz)%Cell%mp>0) THEN
                          !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Berg sein können
                          NrMP_Cells=NrMP_Cells+1
                          Floor(ib)%Cell(ix,iy,iz)%Cell%mp=NrMP_Cells
                        ELSE
                          Floor(ib)%Cell(ix,iy,iz)%Cell%mp=0
                        END IF  
                      END IF  
                    END IF
                  END IF
               END DO
              END DO
              DO iy=iy0+1,iy1
                DO iz=iz0+1,iz1
                   IF (ASSOCIATED(Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face)) THEN
                     IF(Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face%in_out>-4 .and. &
                        Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face%ec>0) THEN
                       ALLOCATE(Floor(ib)%Faces_YZ(ix-facf,iy,iz)%Face)
                       Floor(ib)%Faces_YZ(ix-facf,iy,iz)%Face%mp= &
                             Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face%mp
                       Floor(ib)%Faces_YZ(ix-facf,iy,iz)%Face%MidPoint= &
                             Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face%MidPoint
                       Floor(ib)%Faces_YZ(ix-facf,iy,iz)%Face%Vol= &
                             Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face%Vol
                       Floor(ib)%Faces_YZ(ix-facf,iy,iz)%Face%in_out= &
                             Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face%in_out
                       Floor(ib)%Faces_YZ(ix-facf,iy,iz)%Face%ec= &
                             Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face%ec
    Write(OutUnitProt,*) "FaceYZRandVergl iw/ew ibn,ib,jx-facf,iy,iz", ibn,ib,jx-facf,iy,iz," FaceYZ gl. Refine"
                       IF(Floor(ibn)%Faces_YZ(jx-facf,iy,iz)%Face%Vol>0) THEN
                         Floor(ib)%NrR_FacesYZ=Floor(ib)%NrR_FacesYZ+1 
                         Floor(ib)%NrW_FacesYZ=Floor(ib)%NrW_FacesYZ+1
                       END IF 
                     END IF
                   END IF
               END DO
              END DO
            ELSE  !!(RefineY<RefineNachbarY)
            !-------------------------------
              DO iy=iy0+1,iy1
                jy=iy*IncrY
                DO iz=iz0+1,iz1
                  jz=iz
                  IF(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                     ASSOCIATED(Floor(ibn)%Cell(jx,jy-1,jz)%Cell)) THEN
                    !................................................. 
                    CALL EvalBordCell_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
                    !................................................. 
                  END IF
                 !Vol=0.0d0;SumMidP%x=0.0d0;SumMidP%y=0.0d0;SumMidP%z=0.0d0
                 !IF (ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell)) THEN
                 !  Vol=Floor(ibn)%Cell(jx,jy,jz)%Cell%Vol
                 !  SumMidP=Floor(ibn)%Cell(jx,jy,jz)%Cell%Vol* &
                 !          Floor(ibn)%Cell(jx,jy,jz)%Cell%MidPoint
                 !END IF
                 !IF (ASSOCIATED(Floor(ibn)%Cell(jx,jy-1,jz)%Cell)) THEN
                 !  Vol=Vol+Floor(ibn)%Cell(jx,jy-1,jz)%Cell%Vol
                 !  SumMidP=SumMidP+Floor(ibn)%Cell(jx,jy-1,jz)%Cell%Vol* &
                 !          Floor(ibn)%Cell(jx,jy-1,jz)%Cell%MidPoint
                 !END IF
                 !IF (Vol>0.0d0) THEN  !!!! IF(Vol>0.0d0 .AND. Vol<maxVol(noch berechnen) THEN
                 !  ALLOCATE(Floor(ib)%Cell(ix,iy,iz)%Cell)
                 !  Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=Vol
                 !  Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint= SumMidP/Vol 
                 !  !für Weight_ausgabe gesetzt, für exakte Berechnung InOutCheck der GesamtCelle
                 !  ! nachtragen,geht aber nicht im Mechanismus ein 
                 !  Floor(ib)%Cell(ix,iy,iz)%Cell%in_out=5 !nur für Weight_ausgabe gesetzt
                 !  Floor(ib)%Cell(ix,iy,iz)%Cell%vc=1 !theoretisch/symbolisch für Weight-Ausgabe
                 !  NrMP_Cells=NrMP_Cells+1
                 !  Floor(ib)%Cell(ix,iy,iz)%Cell%mp=NrMP_Cells
                 !  NrR_Cells=NrR_Cells+1
                 !END IF
                END DO
              END DO
            END IF  !(RefineY[>;==;<]RefineNachbarY)

          !!!----
          ELSE IF (RefineZ<RefineNachbarZ) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!      -------------           ------------- 
          !!      |  |  |     |  'iw,ie'  |     |  |  |
          !!      |---- |     |           |     |-----|
          !!      |  |  |     |           |     |  |  |
          !!      -------------           -------------
            IF (RefineY>=RefineNachbarY) THEN
            !-------------------------------
              DO jy=jy0+1,jy1
                jz=jz0+1
                DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND. &
                         .NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz+1)%Cell).AND. &
                         jz<jz1)
                   jz=jz+2
                END DO
                DO WHILE((ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                          ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz+1)%Cell)) .AND. &
                          jz<jz1)
                  iy=IncrY*jy
                  iz=(jz+1)/IncrZ
                  !................................................. 
                  DO py=1,IncrY
                    CALL EvalBordCell_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
                    iy=iy-1
                  END DO
                  !................................................. 
                  jz=jz+2
                END DO  ! DO WHILE((ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR.
                        !           ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz+1)%Cell)) .AND. ...          
              END DO   !jy, NachbarBlock-Cell             

            ELSE  !!(RefineY<RefineNachbarY)
            !-------------------------------
              DO iy=iy0+1,iy1
                jy=iy*IncrY
                DO iz=iz0+1,iz1
                  jz=iz*IncrZ
                  IF(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                     ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz-1)%Cell).OR. &
                     ASSOCIATED(Floor(ibn)%Cell(jx,jy-1,jz)%Cell).OR. &
                     ASSOCIATED(Floor(ibn)%Cell(jx,jy-1,jz-1)%Cell)) THEN
                    !................................................. 
                    CALL EvalBordCell_WE(ib,ix,vix,iy,iz,ibn,vjx,fac)
                    !................................................. 
                  END IF
                END DO
              END DO
            END IF  ! (RefineY [>;==;<] RefineNachbarY)
            !-----
          END IF !! (RefineZ [>;==;<] RefineNachbarZ)

        END IF  !! west,east

        !Zählung Cellen Rand unterhalb Berg 'w,e' 'i and o'
        IF (nType(2:2)=='w') THEN
          i=ix0
          bli=ix0+1
        ELSE
          i=ix1+1
          bli=ix1
        END IF
        DO j=iy0+1,iy1
          !k=iz0+1
          !in_out=Vertices(i-1,j-1, k-1)%in_out+Vertices(i  ,j-1, k-1)%in_out &
          !      +Vertices(i-1,j-1, k  )%in_out+Vertices(i  ,j-1, k  )%in_out &
          !      +Vertices(i-1,j  , k-1)%in_out+Vertices(i  ,j  , k-1)%in_out &
          !      +Vertices(i-1,j  , k  )%in_out+Vertices(i  ,j  , k  )%in_out
          !DO WHILE(.NOT.ASSOCIATED(Cell(i,j,k)%Cell).AND.in_out/=8.AND.k/=iz1+1)
         DO k=iz0+1,iz1
          IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
            in_out=Vertices(i-1,j-1, k-1)%in_out+Vertices(i  ,j-1, k-1)%in_out &
                  +Vertices(i-1,j-1, k  )%in_out+Vertices(i  ,j-1, k  )%in_out &
                  +Vertices(i-1,j  , k-1)%in_out+Vertices(i  ,j  , k-1)%in_out &
                  +Vertices(i-1,j  , k  )%in_out+Vertices(i  ,j  , k  )%in_out
            !w IF (ASSOCIATED(Cell(bli,j,k)%Cell).AND.(Cell(bli,j,k)%Cell%in_out>=6)) THEN
            !w   !Write(*,*) " EXIT aus Zählung NrRN_Cells w/e AnalyzeAllCellsRand"
            !w   !Write(*,*) " Randzählung : i,j,k=",i,j,k
            !w   Write(OutUnitProt,*) " Info : EXIT aus Zählung NrRN_Cells w/e AnalyzeAllCellsRand"
            !w   Write(OutUnitProt,*) "        Randzählung : i,j,k=",i,j,k,"  Block : ",ib
            !w   EXIT
            !w ELSE 
              !IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell).AND.in_out<-6) THEN
            !w  IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell).AND.in_out<=0) THEN
            !w      NrRN_Cells=NrRN_Cells+1
            !w   END IF
              !k=k+1
              IF (in_out<=0) THEN
                  NrRN_Cells=NrRN_Cells+1
              END IF
            !w END IF
          END IF
          END DO
        END DO
        CALL CountFaceSubMontane_WE(ib)
      !......................................................................................
      !...........................TypeN/TypeS................................................
      !......................................................................................
      CASE ("n","s")      !? pn,ps !? läüft  unter else wie 'in,is'
        IF (nType=='on') THEN
          DO ix=ix0+1,ix1
            DO iz=iz0+1,iz1
              !DO WHILE(.NOT.ASSOCIATED(Cell(ix,iy1,iz)%Cell).AND.iz/=iz1)
              !   !Write(*,*) "ib=",ib, " Type-on:  ix=",ix," iy1=",iy1," iz=",iz, &
              !   !           "  RandCell unter Berg north"
              !  !!ALLOCATE(Floor(ib)%Cell(ix,iy1+1,iz)%Cell)
              !   !!Floor(ib)%Cell(ix,iy1+1,iz)%Cell%in_out=-8
              !   !!Floor(ib)%Cell(ix,iy1+1,iz)%Cell%Vol=0.0d0
              !   !NrRN_Cells=NrRN_Cells+1
              !   iz=iz+1
              !END DO 
             IF (ASSOCIATED(Cell(ix,iy1,iz)%Cell).AND.iz/=iz1+1 ) THEN
               ! Problem Celle-Edge-Kante-inout=6, dennoch allokiert, maxVol wurde berechnet 
               !IF(Cell(ix,iy1,iz)%Cell%vc>0) THEN
               IF(Cell(ix,iy1,iz)%Cell%in_out<6.AND.Cell(ix,iy1,iz)%Cell%Vol>0.0d0.OR. &
                  (Cell(ix,iy1,iz)%Cell%in_out==6.AND.Cell(ix,iy1,iz)%Cell%vc>0) ) THEN
                 ALLOCATE(Floor(ib)%Cell(ix,iy1+1,iz)%Cell)
                 Floor(ib)%Cell(ix,iy1+1,iz)%Cell%Vol=Floor(ib)%Cell(ix,iy1,iz)%Cell%Vol
                 Floor(ib)%Cell(ix,iy1+1,iz)%Cell%MidPoint= &
                                         Floor(ib)%Cell(ix,iy1,iz)%Cell%MidPoint
                 Floor(ib)%Cell(ix,iy1+1,iz)%Cell%MidPoint%y= &
                                         Floor(ib)%Cell(ix,iy1,iz)%Cell%MidPoint%y+dy(iy1)
                 Floor(ib)%Cell(ix,iy1+1,iz)%Cell%in_out=Floor(ib)%Cell(ix,iy1,iz)%Cell%in_out
                 Floor(ib)%Cell(ix,iy1+1,iz)%Cell%vc=Floor(ib)%Cell(ix,iy1,iz)%Cell%vc
                 NrR_Cells=NrR_Cells+1
                 IF(Floor(ib)%Cell(ix,iy1,iz)%Cell%Vol>0.0d0) THEN
                    IF(Floor(ib)%Cell(ix,iy1,iz)%Cell%mp>0) THEN
                      !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Berg sein können
                      NrMP_Cells=NrMP_Cells+1
                      Floor(ib)%Cell(ix,iy1+1,iz)%Cell%mp=NrMP_Cells
                    ELSE
                      Floor(ib)%Cell(ix,iy1+1,iz)%Cell%mp=0
                    END IF
                 END IF
               END IF
               !Write(*,*) "ib=",ib, " Type-on:  ix=",ix," iy1=",iy1," iz=",iz, &
               !           " Schnitt RandCelle north"
               !iz=iz+1
             END IF
            END DO !iz
          END DO  !ix
          CALL EvalBorderAllFace_OutN(ib,ix0,ix1,iy1,iz0,iz1)
        ELSE IF(nType=="os") THEN
          DO ix=ix0+1,ix1
             DO iz=iz0+1,iz1
              !DO WHILE(.NOT.ASSOCIATED(Cell(ix,iy0+1,iz)%Cell).AND.iz/=iz1)
              !  !Write(*,*) "ib=",ib, " Type-os:  ,ix=",ix, "iy0+1=",iy0+1,"iz=",iz, &
              !  !           "  RandCell unter Berg south"
              !  !!ALLOCATE(Floor(ib)%Cell(ix,iy0+1,iz)%Cell)
              !  !!Floor(ib)%Cell(ix,iy0+1,iz)%Cell%in_out=-8
              !  !!Floor(ib)%Cell(ix,iy0+1,iz)%Cell%Vol=0.0d0
              !  !NrRN_Cells=NrRN_Cells+1
              !  iz=iz+1
              !END DO 
             IF (ASSOCIATED(Cell(ix,iy0+1,iz)%Cell).AND.iz/=iz1+1) THEN
               ! Problem Celle-Edge-Kante-inout=6, dennoch allokiert, maxVol wurde berechnet 
               !IF(Cell(ix,iy0+1,iz)%Cell%vc>0) THEN
               IF(Cell(ix,iy0+1,iz)%Cell%in_out<6.AND.Cell(ix,iy0+1,iz)%Cell%Vol>0.0d0.OR. &
                  (Cell(ix,iy0+1,iz)%Cell%in_out==6.AND.Cell(ix,iy0+1,iz)%Cell%vc>0) ) THEN
                 ALLOCATE(Floor(ib)%Cell(ix,iy0,iz)%Cell)
                 Floor(ib)%Cell(ix,iy0,iz)%Cell%Vol=Floor(ib)%Cell(ix,iy0+1,iz)%Cell%Vol
                 Floor(ib)%Cell(ix,iy0,iz)%Cell%MidPoint= &
                                     Floor(ib)%Cell(ix,iy0+1,iz)%Cell%MidPoint
                 Floor(ib)%Cell(ix,iy0,iz)%Cell%MidPoint%y= &
                                     Floor(ib)%Cell(ix,iy0+1,iz)%Cell%MidPoint%y-dy(iy0+1)
                 Floor(ib)%Cell(ix,iy0,iz)%Cell%in_out=Floor(ib)%Cell(ix,iy0+1,iz)%Cell%in_out
                 Floor(ib)%Cell(ix,iy0,iz)%Cell%vc=Floor(ib)%Cell(ix,iy0+1,iz)%Cell%vc
                 NrR_Cells=NrR_Cells+1
                 IF(Floor(ib)%Cell(ix,iy0+1,iz)%Cell%Vol>0.0d0) THEN
                   IF(Floor(ib)%Cell(ix,iy0+1,iz)%Cell%mp>0) THEN
                     !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Berg sein können
                     NrMP_Cells=NrMP_Cells+1
                     Floor(ib)%Cell(ix,iy0,iz)%Cell%mp=NrMP_Cells
                   ELSE
                     Floor(ib)%Cell(ix,iy0,iz)%Cell%mp=0
                   END IF
                 END IF
               END IF
               !Write(*,*) "ib=",ib, " Type-os:  ix=",ix,"iy0+1=",iy0+1, "iz=",iz, &
               !           " Schnitt RandCelle south"
               !iz=iz+1
              END IF
            END DO  !iz
          END DO   !ix
             CALL EvalBorderAllFace_OutS(ib,ix0,ix1,iy0,iz0,iz1)

        ELSE  !!'in'.OR.'is' 
          IF (nType(2:2)=='s') THEN
            iy=iy0               !Plane-akt.RandCells
            jy=Floor(ibn)%iy1    !Plane-NachbarCells (jy1)
            dy(iy0)=Floor(ibn)%dy(jy)
            !........
            viy=iy0              !Plane-Point-P0 Celle
            vjy=Floor(ibn)%iy1   !Pos. dy-Nachbar
            fac=-1               !Faktor Position Point (P1,P3,P5,P7)
            facf=1
          ELSE !"n"
            iy=iy1+1             !Plane-akt.RandCells
            jy=Floor(ibn)%iy0+1  !Plane-NachbarCells (jy0+1)
            dy(iy1+1)=Floor(ibn)%dy(jy)  
            !........ 
            viy=iy1              !Plane-Point-P0 Celle
            vjy=Floor(ibn)%iy0+1 !Pos. dy-Nachbar
            fac=1                !Faktor Position Point (P1,P3,P5,P7)
            facf=0
          END IF
          !!!----
          IF (RefineZ>RefineNachbarZ) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!     -------------           -------------
          !!     |     |  |  |  'is,in'  |  |  |     |
          !!     |     |-----|           |-----|     |
          !!     |     |  |  |           |  |  |     |
          !!     -------------           -------------
   
            IF (RefineX>=RefineNachbarX) THEN
            !-------------------------------
            !!    ...7------6---....-----6------7...     !Point-Folge P0-P7 's'/'n' allg.
            !!      /|     /|           /|     /|  
            !!  ...5------4-----....---4------5...
            !!    .|.3....|.2..........|.2....|.3...
            !!     |/ 's' |/   (ib)    |/ 'n' |/
            !!  ...1------0----....----0------1...
            !!   (ibn)   /   (ib)     /  (ibn) 
   
              DO jx=jx0+1,jx1
                jz=jz0+1
                DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND.jz/=jz1)
                   jz=jz+1
                END DO
                DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND.jz/=jz1+1)
                   ix=IncrX*jx
                   iz=IncrZ*jz
                   !............................
                   DO pz=1,IncrZ
                     DO px=1,IncrX
                       CALL EvalBordCell_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
                       ix=ix-1
                     END DO
                     iz=iz-1
                     ix=IncrX*jx
                   END DO 
                   !............................
                   jz=jz+1
                END DO  !DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell...
              END DO   !jx, NachbarBlock-Cell             
            ELSE IF (RefineX<RefineNachbarX) THEN
            !-------------------------------
              DO jx=jx0+1,jx1,2
                jz=jz0+1
                DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND. &
                         .NOT.ASSOCIATED(Floor(ibn)%Cell(jx+1,jy,jz)%Cell).AND. &
                         jz/=jz1+1)
                   jz=jz+1
                END DO
                DO WHILE((ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                          ASSOCIATED(Floor(ibn)%Cell(jx+1,jy,jz)%Cell)) .AND. &
                          jz/=jz1+1)
                   ix=(jx+1)/IncrX
                   iz=IncrZ*jz
                   !............................
                   DO pz=1,IncrZ
                      CALL EvalBordCell_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
                      iz=iz-1
                   END DO
                   !............................
                   jz=jz+1
                END DO  !DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell... .OR.
                        !         ASSOCIATED(Floor(ibn)%Cell(jx,jy+1,jz)%Cell...) ...
              END DO  !jx, NachbarBlock-Cell             
            END IF  !(RefineX[>;==;<]RefineNachbarX)
   
          !!!----
          ELSE IF (RefineZ==RefineNachbarZ) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!      -------------
          !!      |     |     |   'is,in'
          !!      |---------- |
          !!      |     |     |
          !!      -------------
            IF (RefineX>RefineNachbarX) THEN
            !-------------------------------
              DO jx=jx0+1,jx1
                jz=jz0+1
                DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND.jz/=jz1)
                   jz=jz+1
                END DO
                DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND.jz/=jz1+1)
                  ix=IncrX*jx
                  iz=IncrZ*jz
                  !............................
                  DO px=1,IncrX
                    CALL EvalBordCell_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
                    ix=ix-1
                  END DO
                  !............................
                  jz=jz+1
                END DO  !DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell...
              END DO   !jx, NachbarBlock-Cell             
   
            ELSE IF (RefineX==RefineNachbarX) THEN
            !-------------------------------------
              DO ix=ix0+1,ix1
                DO iz=iz0+1,iz1
                 ! Vol=0.0d0;SumMidP%x=0.0d0;SumMidP%y=0.0d0;SumMidP%z=0.0d0
                 IF (ASSOCIATED(Floor(ibn)%Cell(ix,jy,iz)%Cell)) THEN
                   ALLOCATE(Floor(ib)%Cell(ix,iy,iz)%Cell)
                   Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=Floor(ibn)%Cell(ix,jy,iz)%Cell%Vol
                   Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint= &
                                                Floor(ibn)%Cell(ix,jy,iz)%Cell%MidPoint
                   Floor(ib)%Cell(ix,iy,iz)%Cell%in_out= &
                                                Floor(ibn)%Cell(ix,jy,iz)%Cell%in_out
                   Floor(ib)%Cell(ix,iy,iz)%Cell%vc= &
                                                Floor(ibn)%Cell(ix,jy,iz)%Cell%vc
    Write(OutUnitProt,*) "CellRandVergl in/is ibn,ib,ix,jy,iz", ibn,ib,ix,jy,iz," Celle gl. Refine"
                   !Floor(ib)%Cell(ix,iy,iz)%Cell=Floor(ibn)%Cell(ix,jy,iz)%Cell
                   IF(Floor(ibn)%Cell(ix,jy,iz)%Cell%Vol>0.0d0) THEN
                     NrR_Cells=NrR_Cells+1
                     IF(Floor(ibn)%Cell(ix,jy,iz)%Cell%vc>0) THEN
                       IF(Floor(ibn)%Cell(ix,jy,iz)%Cell%mp>0) THEN
                          !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Berg sein können
                          NrMP_Cells=NrMP_Cells+1
                          Floor(ib)%Cell(ix,iy,iz)%Cell%mp=NrMP_Cells
                       ELSE
                          Floor(ib)%Cell(ix,iy,iz)%Cell%mp=0
                       END IF
                     END IF
                   END IF
                 END IF
                END DO
              END DO
              DO ix=ix0+1,ix1
                DO iz=iz0+1,iz1
                   IF (ASSOCIATED(Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face)) THEN
                     IF(Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face%in_out>-4 .and. &
                        Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face%ec>0) THEN
                       ALLOCATE(Floor(ib)%Faces_ZX(ix,iy-facf,iz)%Face)
                       Floor(ib)%Faces_ZX(ix,iy-facf,iz)%Face%mp= &
                             Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face%mp
                       Floor(ib)%Faces_ZX(ix,iy-facf,iz)%Face%MidPoint= &
                             Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face%MidPoint
                       Floor(ib)%Faces_ZX(ix,iy-facf,iz)%Face%Vol= &
                             Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face%Vol
                       Floor(ib)%Faces_ZX(ix,iy-facf,iz)%Face%in_out= &
                             Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face%in_out
                       Floor(ib)%Faces_ZX(ix,iy-facf,iz)%Face%ec= &
                             Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face%ec
    Write(OutUnitProt,*) "FaceZXRandVergl in/is ibn,ib,ix,jy-facf,iz", ibn,ib,ix,jy-facf,iz," FaceZX gl. Refine"
                       IF(Floor(ibn)%Faces_ZX(ix,jy-facf,iz)%Face%Vol>0) THEN
                         Floor(ib)%NrR_FacesZX=Floor(ib)%NrR_FacesZX+1 
                         Floor(ib)%NrW_FacesZX=Floor(ib)%NrW_FacesZX+1
                       END IF 
                     END IF
                   END IF
                END DO
              END DO
            ELSE  !!(RefineX<RefineNachbarX)
            !-------------------------------
               DO ix=ix0+1,ix1
                 jx=ix*IncrX
                 jz=jz0+1
                 DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND. &
                          .NOT.ASSOCIATED(Floor(ibn)%Cell(jx-1,jy,jz)%Cell).AND.jz/=jz1)
                    jz=jz+1
                 END DO
                 DO WHILE((ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                           ASSOCIATED(Floor(ibn)%Cell(jx-1,jy,jz)%Cell)).AND.jz/=jz1+1)
                   iz=IncrZ*jz    !IncrZ==1
                   !............................
                   CALL EvalBordCell_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
                   !............................
                   jz=jz+1
                 END DO  !DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell...
               END DO   !ix-akt.Cell-ib             
            END IF  !(RefineX[>;==;<]RefineNachbarX)
   
          ELSE IF (RefineZ<RefineNachbarZ) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!      -------------           ------------- 
          !!      |  |  |     |  'is,in'  |     |  |  |
          !!      |---- |     |           |     |-----|
          !!      |  |  |     |           |     |  |  |
          !!      -------------           -------------
            IF (RefineX>=RefineNachbarX) THEN
            !-------------------------------
               DO jx=jx0+1,jx1
                 jz=jz0+1
                 DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND. &
                          .NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz+1)%Cell).AND. &
                          jz<jz1-1)
                    jz=jz+2
                 END DO
                 DO WHILE((ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                           ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz+1)%Cell)) .AND. &
                           jz<jz1)
                   ix=IncrX*jx
                   iz=(jz+1)/IncrZ
                   !............................
                   DO px=1,IncrX
                     CALL EvalBordCell_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
                     ix=ix-1
                   END DO
                   !............................
                   jz=jz+2
                 END DO  ! DO WHILE((ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR.
                         !           ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz+1)%Cell)) .AND. ...          
               END DO   !jx, NachbarBlock-Cell             
            ELSE  !!(RefineX<RefineNachbarX)
            !-------------------------------
               DO ix=ix0+1,ix1
                 jx=ix*IncrX
                 jz=jz0+1
                 DO WHILE(.NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).AND. &
                          .NOT.ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz+1)%Cell).AND. &
                          .NOT.ASSOCIATED(Floor(ibn)%Cell(jx-1,jy,jz)%Cell).AND. &
                          .NOT.ASSOCIATED(Floor(ibn)%Cell(jx-1,jy,jz+1)%Cell).AND. &
                          jz<jz1-1)
                    jz=jz+2
                 END DO
                 DO WHILE(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                          ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz+1)%Cell).OR. &
                          ASSOCIATED(Floor(ibn)%Cell(jx-1,jy,jz)%Cell).OR. &
                          ASSOCIATED(Floor(ibn)%Cell(jx-1,jy,jz+1)%Cell).OR. &
                          jz<jz1)
                   iz=(jz+1)/IncrZ
                   !............................
                   CALL EvalBordCell_NS(ib,ix,iy,viy,iz,ibn,vjy,fac)
                   !............................
                   jz=jz+2
                 END DO
               END DO
            END IF  ! (RefineX [>;==;<] RefineNachbarX)
            !-----
          END IF !! (RefineZ [>;==;<] RefineNachbarZ)

        END IF  !! south,north, (os,on,in-is) 

        !Zählung Cellen Rand unterhalb Berg 'n,s' 'i and o'
        IF (nType(2:2)=='s') THEN
          j=iy0
          blj=iy0+1
        ELSE
          j=iy1+1
          blj=jy1
        END IF
        DO i=ix0+1,ix1
          k=iz0+1
          in_out=Vertices(i-1,j-1, k-1)%in_out+Vertices(i  ,j-1, k-1)%in_out &
                +Vertices(i-1,j-1, k  )%in_out+Vertices(i  ,j-1, k  )%in_out &
                +Vertices(i-1,j  , k-1)%in_out+Vertices(i  ,j  , k-1)%in_out & 
                +Vertices(i-1,j  , k  )%in_out+Vertices(i  ,j  , k  )%in_out
          DO k=iz0+1,iz1
          !DO WHILE(.NOT.ASSOCIATED(Cell(i,j,k)%Cell).AND.in_out/=8.AND.k/=iz1+1)
          IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
            in_out=Vertices(i-1,j-1, k-1)%in_out+Vertices(i  ,j-1, k-1)%in_out &
                  +Vertices(i-1,j-1, k  )%in_out+Vertices(i  ,j-1, k  )%in_out &
                  +Vertices(i-1,j  , k-1)%in_out+Vertices(i  ,j  , k-1)%in_out & 
                  +Vertices(i-1,j  , k  )%in_out+Vertices(i  ,j  , k  )%in_out
           !w IF (ASSOCIATED(Cell(i,blj,k)%Cell).AND.Cell(i,blj,k)%Cell%in_out>=6) THEN
           !w    ! Problem Erste bzw. letzte Celle allokiert da Schnitt und 
           !w    ! dazugehörige RandCelle keinen Schnitt Zählung Überlauf

           !w  ! Write(*,*) " EXIT aus Zählung NrRN_Cells n/s AnalyzeAllCellsRand"
           !w  ! Write(*,*) " iandzählung : i,j,k=",i,j,k
           !w   Write(OutUnitProt,*) " Info : EXIT aus Zählung NrRN_Cells n/s AnalyzeAllCellsRand"
           !w   Write(OutUnitProt,*) "        Randzählung : i,j,k=",i,j,k,"  Block : ",ib
 
           !w    EXIT
           !w ELSE 
           !w   IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell).AND.in_out<=0) THEN
           !w      NrRN_Cells=NrRN_Cells+1
           !w   END IF
              !k=k+1
           !w END IF
            IF (in_out<=0) THEN
                 NrRN_Cells=NrRN_Cells+1
            END IF
          END IF
          END DO !k
        END DO !i
        CALL CountFaceSubMontane_NS(ib)

      !......................................................................................
      !...........................TypeT/TypeB................................................
      !......................................................................................
      CASE ("t","b")      !? pt,pb !? läüft  unter else wie 'it,ib'
        IF (nType=='ot') THEN
          DO ix=ix0+1,ix1
            DO iy=iy0+1,iy1
              IF (ASSOCIATED(Cell(ix,iy,iz1)%Cell)) THEN
                 !Floor(ib)%Cell(ix,iy,iz1+1)%Cell=>Floor(ib)%Cell(ix,iy,iz1)%Cell
                 !IF(Cell(ix,iy,iz1)%Cell%vc>0) THEN
                 IF((Cell(ix,iy,iz1)%Cell%in_out<6.AND.Cell(ix,iy,iz1)%Cell%Vol>0.0d0) .OR. & 
                    (Cell(ix,iy,iz1)%Cell%in_out==6.AND.Cell(ix,iy,iz1)%Cell%vc>0) ) THEN
                   ALLOCATE(Floor(ib)%Cell(ix,iy,iz1+1)%Cell)
                   Floor(ib)%Cell(ix,iy,iz1+1)%Cell%Vol=Floor(ib)%Cell(ix,iy,iz1)%Cell%Vol
                   Floor(ib)%Cell(ix,iy,iz1+1)%Cell%MidPoint= &
                                           Floor(ib)%Cell(ix,iy,iz1)%Cell%MidPoint
                   Floor(ib)%Cell(ix,iy,iz1+1)%Cell%MidPoint%z= &
                                           Floor(ib)%Cell(ix,iy,iz1)%Cell%MidPoint%z+dz(iz1)
                   Floor(ib)%Cell(ix,iy,iz1+1)%Cell%in_out=Floor(ib)%Cell(ix,iy,iz1)%Cell%in_out
                   Floor(ib)%Cell(ix,iy,iz1+1)%Cell%vc=Floor(ib)%Cell(ix,iy,iz1)%Cell%vc
                   NrR_Cells=NrR_Cells+1
                   IF(Floor(ib)%Cell(ix,iy,iz1)%Cell%Vol>0.0d0) THEN
                     IF(Floor(ib)%Cell(ix,iy,iz1)%Cell%mp>0) THEN
                       !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Berg sein können
                        NrMP_Cells=NrMP_Cells+1
                        Floor(ib)%Cell(ix,iy,iz1+1)%Cell%mp=NrMP_Cells
                     ELSE
                        Floor(ib)%Cell(ix,iy,iz1+1)%Cell%mp=0
                     END IF
                   END IF
                   !Write(*,*) "ib=",ib, " Type-ot:  ix=",ix," iy=",iy," iz1=",iz1, &
                   !           " Schnitt RandCelle top"
                 END IF
              END IF
            END DO 
          END DO
          CALL EvalBorderAllFace_OutT(ib,ix0,ix1,iy0,iy1,iz1) 
        ELSE IF(nType=="ob") THEN
          DO ix=ix0+1,ix1
            DO iy=iy0+1,iy1
                !Write(*,*) "ib=",ib, " Type-ob:  ,ix=",ix, "iy=",iy,"iz=",iz0, &
                !           "  RandCell unter Berg bottom"
                !!ALLOCATE(Floor(ib)%Cell(ix,iy,iz0)%Cell)
                !!Floor(ib)%Cell(ix,iy,iz0)%Cell%in_out=-8
                !!Floor(ib)%Cell(ix,iy,iz0)%Cell%Vol=0.0d0
!27.10.11                NrRN_Cells=NrRN_Cells+1
              IF (ASSOCIATED(Cell(ix,iy,iz0+1)%Cell)) THEN
                 !Floor(ib)%Cell(ix,iy,iz1+1)%Cell=>Floor(ib)%Cell(ix,iy,iz1)%Cell
                 !IF(Cell(ix,iy,iz1)%Cell%vc>0) THEN
                 IF((Cell(ix,iy,iz0+1)%Cell%in_out<6.AND.Cell(ix,iy,iz0+1)%Cell%Vol>0.0d0 &
                     .AND.Cell(ix,iy,iz0+1)%Cell%vc>0) .OR. &
                    (Cell(ix,iy,iz0+1)%Cell%in_out==6.AND.Cell(ix,iy,iz0+1)%Cell%vc>0) ) THEN
                   ALLOCATE(Floor(ib)%Cell(ix,iy,iz0)%Cell)
                   Floor(ib)%Cell(ix,iy,iz0)%Cell%Vol=Floor(ib)%Cell(ix,iy,iz0+1)%Cell%Vol
                   Floor(ib)%Cell(ix,iy,iz0)%Cell%MidPoint= &
                                           Floor(ib)%Cell(ix,iy,iz0+1)%Cell%MidPoint
                   Floor(ib)%Cell(ix,iy,iz0)%Cell%MidPoint%z= &
                                           Floor(ib)%Cell(ix,iy,iz0+1)%Cell%MidPoint%z+dz(iz1)
                   Floor(ib)%Cell(ix,iy,iz0)%Cell%in_out=Floor(ib)%Cell(ix,iy,iz0+1)%Cell%in_out
                   Floor(ib)%Cell(ix,iy,iz0)%Cell%vc=Floor(ib)%Cell(ix,iy,iz0+1)%Cell%vc
                   NrR_Cells=NrR_Cells+1
                   IF(Floor(ib)%Cell(ix,iy,iz0+1)%Cell%Vol>0.0d0) THEN
                     IF(Floor(ib)%Cell(ix,iy,iz0+1)%Cell%mp>0) THEN
                       !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Berg sein können
                        NrMP_Cells=NrMP_Cells+1
                        Floor(ib)%Cell(ix,iy,iz0)%Cell%mp=NrMP_Cells
                     ELSE
                        Floor(ib)%Cell(ix,iy,iz0)%Cell%mp=0
                     END IF
                   END IF
                   !Write(*,*) "ib=",ib, " Type-ot:  ix=",ix," iy=",iy," iz1=",iz1, &
                   !           " Schnitt RandCelle botton"
                 END IF
              END IF
            END DO 
          END DO
          CALL EvalBorderAllFace_OutB(ib,ix0,ix1,iy0,iy1,iz0) 

        ELSE  !!'ib'.OR.'it' 
           IF (nType(2:2)=='b') THEN
             iz=iz0               !Plane-akt.RandCells
             jz=Floor(ibn)%iz1    !Plane-NachbarCells
             dz(iz0)=Floor(ibn)%dz(jz)
             !........
             viz=iz0              !Plane-Point-P0 Celle
             vjz=Floor(ibn)%iz1   !Pos. dz-Nachbar
             fac=-1               !Faktor Position Point (P1,P3,P5,P7) 
             facf=1               !Faktor Position für Face
           ELSE !"t"
             iz=iz1+1             !Plane-akt.RandCells
             jz=Floor(ibn)%iz0+1  !Plane-NachbarCells
             dz(iz1+1)=Floor(ibn)%dz(jz)  
             !........ 
             viz=iz1              !Plane-Point-P0 Celle
             vjz=Floor(ibn)%iz0+1 !Pos. dz-Nachbar
             fac=1                !Faktor Position Point (P1,P3,P5,P7)
             facf=0               !Faktor Position für Face
           END IF

          !!!----
          IF (RefineY>RefineNachbarY) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!     ---------------------------
          !!     |      |    |     |      |
          !!     | 'ib' |--akt.ib--| 'it' |
          !!     |      |    |     |      |
          !!     ---------------------------
             IF (RefineX>=RefineNachbarX) THEN
             !-------------------------------
             !!    ...7------6---....-----6------7...     !Point-Folge P0-P7 'b'/'t' allg.
             !!      /|     /|           /|     /|  
             !!  ...5------4-----....---4------5...
             !!    .|.3....|.2..........|.2....|.3...
             !!     |/ 'b' |/   (ib)    |/ 't' |/
             !!  ...1------0----....----0------1...
             !!   (ibn)   /   (ib)     /  (ibn) 
 
               DO jx=jx0+1,jx1
                 DO jy=jy0+1,jy1
                   IF (ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell)) THEN
                     ix=IncrX*jx
                     iy=IncrY*jy
                     !............................
                     bt_piy: DO py=1,IncrY
                       bt_pix: DO px=1,IncrX
                         ! P0=Floor(ib)%Vertices(ix,iy,viz)%Point
                         ! P1=P0
                         ! P1%z=P0%z+fac*Floor(ibn)%dz(vjz)
                         ! P2=Floor(ib)%Vertices(ix-1,iy,viz)%Point
                         ! P3=P2
                         ! P3%z=P2%z+fac*Floor(ibn)%dz(vjz)
                         ! !.........................................
                         ! P4=Floor(ib)%Vertices(ix,iy-1,viz)%Point
                         ! P5=P4
                         ! P5%z=P4%z+fac*Floor(ibn)%dz(vjz)
                         ! P6=Floor(ib)%Vertices(ix-1,iy-1,viz)%Point
                         ! P7=P6
                         ! P7%z=P6%z+fac*Floor(ibn)%dz(vjz)
                         ! !.........................................
                         ! !Vertex-Zuweisung äußerer Rand + Analyze
                         ! !eventuell als Teil AnalyzeVertices 
                         ! Vertices(ix  ,iy  ,viz+fac)%Point=P1
                         ! Vertices(ix-1,iy  ,viz+fac)%Point=P3
                         ! Vertices(ix  ,iy-1,viz+fac)%Point=P5
                         ! Vertices(ix-1,iy-1,viz+fac)%Point=P7
                         ! CALL CheckVertexRand(Vertices(ix  ,iy  ,viz+fac)) !P1
                         ! CALL CheckVertexRand(Vertices(ix-1,iy  ,viz+fac)) !P3
                         ! CALL CheckVertexRand(Vertices(ix  ,iy-1,viz+fac)) !P5
                         ! CALL CheckVertexRand(Vertices(ix-1,iy-1,viz+fac)) !P7
                         ! !.........................................
                         ! !Zuweisung Celle-Rand für Übergabe
                         ! CellRa%V0=Floor(ib)%Vertices(ix  ,iy  ,viz)
                         ! CellRa%V1=Floor(ib)%Vertices(ix  ,iy  ,viz+fac)
                         ! CellRa%V2=Floor(ib)%Vertices(ix-1,iy  ,viz)
                         ! CellRa%V3=Floor(ib)%Vertices(ix-1,iy  ,viz+fac)
                         ! CellRa%V4=Floor(ib)%Vertices(ix  ,iy-1,viz)
                         ! CellRa%V5=Floor(ib)%Vertices(ix  ,iy-1,viz+fac)
                         ! CellRa%V6=Floor(ib)%Vertices(ix-1,iy-1,viz)
                         ! CellRa%V7=Floor(ib)%Vertices(ix-1,iy-1,viz+fac)
                         ! !.........................................
                         ! CALL AnalyzeCellRand(CellRa)
                         ! !IF(CellRa%in_out>-8.AND.CellRa%in_out<6) THEN  ?NrMP nein NrW ja?
                         ! IF(CellRa%in_out>-8.AND.CellRa%vc>0) THEN
                         !   ALLOCATE(Floor(ib)%Cell(ix,iy,iz)%Cell)
                         !   Cell(ix,iy,iz)%Cell%mp=CellRa%mp
                         !   Cell(ix,iy,iz)%Cell%MidPoint=CellRa%MidPoint
                         !   Cell(ix,iy,iz)%Cell%Vol=CellRa%Vol
                         !   Cell(ix,iy,iz)%Cell%in_out=CellRa%in_out
                         !   Cell(ix,iy,iz)%Cell%vc=CellRa%vc !theoretisch/symbolisch für Weight-Ausgabe
                         ! END IF
                         CALL EvalBordCell_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
                         ix=ix-1
                       END DO bt_pix
                       iy=iy-1
                       ix=IncrX*jx
                     END DO bt_piy
                     !............................
                   END IF    ! IF(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell)
                 END DO  ! jy, NachbarBlock-Cell
               END DO   !jx, NachbarBlock-Cell             
             ELSE IF (RefineX<RefineNachbarX) THEN
             !-------------------------------
               DO jx=jx0+1,jx1,2
                 DO jy=jy0+1,jy1
                   IF (ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                       ASSOCIATED(Floor(ibn)%Cell(jx+1,jy,jz)%Cell)) THEN
                     ix=(jx+1)/IncrX
                     iy=IncrY*jy
                     !............................
                     DO py=1,IncrY
                        CALL EvalBordCell_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
                        iy=iy-1
                     END DO
                     !............................
                   END IF ! IF(ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR.
                          !    ASSOCIATED(Floor(ibn)%Cell(jx+1,jy,jz)%Celli))
                 END DO  !jy, NachbarBlock-Cell             
               END DO  !jx, NachbarBlock-Cell             
             END IF  !(RefineX[>;==;<]RefineNachbarX)

          !!!----
          ELSE IF (RefineY==RefineNachbarY) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!     ---------------------------
          !!     |      |           |      |
          !!     | 'ib' | akt. ib   | 'it' |
          !!     |      |           |      |
          !!     ---------------------------
             IF (RefineX>RefineNachbarX) THEN
             !-------------------------------
               DO jx=jx0+1,jx1
                 DO jy=jy0+1,jy1
                   IF (ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell)) THEN
                      ix=IncrX*jx
                      iy=jy
                      !............................
                      DO px=1,IncrX
                        CALL EvalBordCell_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
                        ix=ix-1
                      END DO
                      !............................
                   END IF
                 END DO  !jy, NachbarBlock-Cell
               END DO   !jx, NachbarBlock-Cell             
             ELSE IF (RefineX==RefineNachbarX) THEN
             !-------------------------------------
               DO ix=ix0+1,ix1
                 DO iy=iy0+1,iy1
                   IF (ASSOCIATED(Floor(ibn)%Cell(ix,iy,jz)%Cell)) THEN
                      ! CALL EvalBordCell_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
                     Floor(ib)%Cell(ix,iy,iz)%Cell%Vol=Floor(ibn)%Cell(ix,iy,jz)%Cell%Vol
                     Floor(ib)%Cell(ix,iy,iz)%Cell%MidPoint= &
                                                  Floor(ibn)%Cell(ix,iy,jz)%Cell%MidPoint
                     Floor(ib)%Cell(ix,iy,iz)%Cell%in_out= &
                                                  Floor(ibn)%Cell(ix,iy,jz)%Cell%in_out
                     Floor(ib)%Cell(ix,iy,iz)%Cell%vc= &
                                                  Floor(ibn)%Cell(ix,iy,jz)%Cell%vc
    Write(OutUnitProt,*) "CellRandVergl it/ib ibn,ib,ix,iy,jz", ibn,ib,ix,iy,jz," Celle gl. Refine"
                     !Floor(ib)%Cell(ix,iy,iz)%Cell=Floor(ibn)%Cell(ix,jy,iz)%Cell
                     IF(Floor(ibn)%Cell(ix,iy,jz)%Cell%Vol>0.0d0) THEN
                       NrR_Cells=NrR_Cells+1
                       IF(Floor(ibn)%Cell(ix,iy,jz)%Cell%vc>0) THEN
                         IF(Floor(ibn)%Cell(ix,iy,jz)%Cell%mp>0) THEN
                            !Zusatzabfrage notwendig,da allokierte Cellen unterhalb Berg sein können
                            NrMP_Cells=NrMP_Cells+1
                            Floor(ib)%Cell(ix,iy,iz)%Cell%mp=NrMP_Cells
                         ELSE
                            Floor(ib)%Cell(ix,iy,iz)%Cell%mp=0
                         END IF
                       END IF
                     END IF
                   END IF
                 END DO
               END DO
               DO ix=ix0+1,ix1
                 DO iy=iy0+1,iy1
                   IF (ASSOCIATED(Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face)) THEN
                     IF(Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face%in_out>-4 .and. &
                        Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face%ec>0) THEN
                       ALLOCATE(Floor(ib)%Faces_XY(ix,iy,iz-facf)%Face)
                       Floor(ib)%Faces_XY(ix,iy,iz-facf)%Face%mp= &
                             Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face%mp
                       Floor(ib)%Faces_XY(ix,iy,iz-facf)%Face%MidPoint= &
                             Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face%MidPoint
                       Floor(ib)%Faces_XY(ix,iy,iz-facf)%Face%Vol= &
                             Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face%Vol
                       Floor(ib)%Faces_XY(ix,iy,iz-facf)%Face%in_out= &
                             Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face%in_out
                       Floor(ib)%Faces_XY(ix,iy,iz-facf)%Face%ec= &
                             Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face%ec
    Write(OutUnitProt,*) "FaceXYRandVergl in/is ibn,ib,ix,iy,jz-facf", ibn,ib,ix,iy,jz-facf," FaceXY gl. Refine"
                       IF(Floor(ibn)%Faces_XY(ix,iy,jz-facf)%Face%Vol>0) THEN
                         Floor(ib)%NrR_FacesXY=Floor(ib)%NrR_FacesXY+1 
                         Floor(ib)%NrW_FacesXY=Floor(ib)%NrW_FacesXY+1
                       END IF 
                     END IF
                   END IF
                 END DO
               END DO
             ELSE  !!(RefineX<RefineNachbarX)
             !-------------------------------
                DO ix=ix0+1,ix1
                  jx=ix*IncrX
                  DO iy=iy0+1,iy1
                     jy=iy
                     IF (ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                         ASSOCIATED(Floor(ibn)%Cell(jx-1,jy,jz)%Cell)) THEN
                        !............................
                        CALL EvalBordCell_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
                        !............................
                     END IF
                  END DO  !iy-akt.Celle-ib
                END DO   !ix-akt.Cell-ib             
             END IF  !(RefineX[>;==;<]RefineNachbarX)

          ELSE IF (RefineY<RefineNachbarY) THEN
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
          !!      ----------- -------------------- 
          !!      |        |           |         |
          !!      |--'ib'--| akt. ib   |--'it'---|
          !!      |        |           |         |
          !!      --------------------------------
          !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
             IF (RefineX>=RefineNachbarX) THEN
             !-------------------------------
               DO jx=jx0+1,jx1
                 DO jy=jy0+1,jy1,2
                   IF (ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                       ASSOCIATED(Floor(ibn)%Cell(jx,jy+1,jz)%Cell)) THEN 
                      ix=IncrX*jx
                      iy=(jy+1)/IncrY
                      !............................
                      DO px=1,IncrX
                        CALL EvalBordCell_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
                        ix=ix-1
                      END DO
                      !............................
                   END IF
                 END DO  !jy, NachbarBlock-Cell  
               END DO   !jx, NachbarBlock-Cell             
             ELSE  !!(RefineX<RefineNachbarX)
             !-------------------------------
               DO ix=ix0+1,ix1
                 jx=ix*IncrX
                 DO iy=iy0+1,iy1
                   jy=iy*IncrY
                   IF (ASSOCIATED(Floor(ibn)%Cell(jx,jy,jz)%Cell).OR. &
                       ASSOCIATED(Floor(ibn)%Cell(jx-1,jy,jz)%Cell).OR. &
                       ASSOCIATED(Floor(ibn)%Cell(jx,jy-1,jz)%Cell).OR. &
                       ASSOCIATED(Floor(ibn)%Cell(jx-1,jy-1,jz)%Cell)) THEN
                      !................................................
                      CALL EvalBordCell_TB(ib,ix,iy,iz,viz,ibn,vjz,fac)
                      !................................................
                   END IF
                 END DO !iy
               END DO !ix
             END IF  ! (RefineX [>;==;<] RefineNachbarX)
             !-----
          END IF !! (RefineY [>;==;<] RefineNachbarY)
        END IF  !! bottom,top, ob,ot,ib-it 

          !Zählung Cellen Rand unterhalb Berg 'it,ot,ib,ob'
          !IF(TypeT=="it") THEN
          IF (nType(2:2)=='t') THEN
             k=iz1+1
             DO i=ix0+1,ix1
               DO j=iy0+1,iy1
                  in_out=Vertices(i-1,j-1,k-1)%in_out+Vertices(i  ,j-1,k-1)%in_out &
                        +Vertices(i-1,j  ,k-1)%in_out+Vertices(i  ,j  ,k-1)%in_out &
                        +Vertices(i-1,j-1,k  )%in_out+Vertices(i  ,j-1,k  )%in_out &
                        +Vertices(i-1,j  ,k  )%in_out+Vertices(i  ,j  ,k  )%in_out
                  IF(.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
                   IF (in_out<=0) THEN
                  !IF (in_out<-6.OR.(in_out==-6 .AND.(.NOT.ASSOCIATED(Cell(i,j,k)%Cell)))) THEN
                      NrRN_Cells=NrRN_Cells+1
                   END IF
                  END IF
               END DO
             END DO
          END IF
          IF(TypeB=="ib") THEN
             k=iz0
             DO i=ix0+1,ix1
               DO j=iy0+1,iy1
                  in_out=Vertices(i-1,j-1,k-1)%in_out +Vertices(i  ,j-1,k-1)%in_out &
                        +Vertices(i-1,j  ,k-1)%in_out +Vertices(i  ,j  ,k-1)%in_out &
                        +Vertices(i-1,j-1,k  )%in_out +Vertices(i  ,j-1,k  )%in_out &
                        +Vertices(i-1,j  ,k  )%in_out +Vertices(i  ,j  ,k  )%in_out
                  IF(.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
                    IF (in_out<=0) THEN
                    !IF (in_out<-6.OR.(in_out==-6.AND.(.NOT.ASSOCIATED(Cell(i,j,k)%Cell)))) THEN
                      NrRN_Cells=NrRN_Cells+1
                    END IF
                  END IF
               END DO
             END DO
          END IF
        !END IF  !! bottom,top, ob,ot,ib-it 

        CALL CountFaceSubMontane_TB(ib)

      END SELECT  ! -w,-e,-s,-n,-b,-t Nachbar
    END DO  !...nachbars

    Floor(ib)%NrR_Cells=NrR_Cells
    Floor(ib)%NrRN_Cells=NrRN_Cells
    !Floor(ib)%NrMP_Cells=NrMP_Cells  !zur Zeit nicht erwünscht nur zu Cellen im Block
    Floor(ib)%NrW_Cells=Floor(ib)%NrW_Cells+NrR_Cells+NrRN_Cells
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrR_Cells = ", Floor(ib)%NrR_Cells
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrRN_Cells = ", Floor(ib)%NrRN_Cells
   !Write(*,*) "ib = ", ib, " NrMP_Cells = ", NrMP_Cells
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrMP_Cells = ", Floor(ib)%NrMP_Cells
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrW_Cells = ", Floor(ib)%NrW_Cells
   !Write(*,*) "....................................................."
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrR_FacesXY = ",Floor(ib)%NrR_FacesXY
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrR_FacesYZ = ",Floor(ib)%NrR_FacesYZ
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrR_FacesZX = ",Floor(ib)%NrR_FacesZX 
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrRN_FacesXY = ",Floor(ib)%NrRN_FacesXY
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrRN_FacesYZ = ",Floor(ib)%NrRN_FacesYZ
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrRN_FacesZX = ",Floor(ib)%NrRN_FacesZX 
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrW_FacesXY = ",Floor(ib)%NrW_FacesXY
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrW_FacesYZ = ",Floor(ib)%NrW_FacesYZ
   !Write(*,*) "ib = ", ib, " Floor(ib)%NrW_FacesZX = ",Floor(ib)%NrW_FacesZX 
   !Write(*,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
  END DO  !...nb

END SUBROUTINE AnalyzeAllCellsRand



SUBROUTINE SortVertFaceIn(Face)
  TYPE (Face_T) :: Face

  INTEGER :: i,j,v,ecut
  INTEGER :: EdgeList(1:2,1:6)

  Face%ec=-1
  Face%NumberVert=0
  Face%VertexList=0
  Face%EdgeCut(:)=0
  Face%in_out=Face%Edge1%in_out+Face%Edge3%in_out
  i=0
  ecut=0
  !...............................................
  IF (Face%Edge1%Vert1%in_out==-1) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge1%Vert1%nrInP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge1%yes_sp==1) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge1%VertS%nrInP
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge1%VertS%nrInP
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%Edge1%VertS%nrInP
  ELSE IF (Face%Edge1%Vert1%in_out==0) THEN
    IF(Face%Edge4%Vert2%in_out==1.and.Face%Edge1%Vert2%in_out==1) THEN
      !Point ohne Verbindung im Face-Inside -> Orography
      !      (1)---x--(-1)  
      !        |      |
      !        |      x
      !      (0)------( 1)
    ELSE
      i=i+1
      Face%VertexList(i)=Face%Edge1%Vert1%nrInP
      Face%NumberVert=Face%NumberVert+1
      ecut=ecut+1
      Face%EdgeCut(ecut)=Face%Edge1%Vert1%nrInP
      Face%ec=1
      EdgeList(1,ecut)=0
      EdgeList(2,ecut)=Face%Edge1%Vert1%nrInP
    END IF
  END IF
  !............................................
  IF (Face%Edge2%Vert1%in_out==-1) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge2%Vert1%nrInP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge2%yes_sp==1) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge2%VertS%nrInP
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge2%VertS%nrInP
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%Edge2%VertS%nrInP
  ELSE IF (Face%Edge2%Vert1%in_out==0) THEN
    IF(Face%Edge1%Vert1%in_out==1.and.Face%Edge2%Vert2%in_out==1) THEN
      !Point ohne Verbindung im Face
    ELSE
      i=i+1
      Face%VertexList(i)=Face%Edge2%Vert1%nrInP
      Face%NumberVert=Face%NumberVert+1
      ecut=ecut+1
      Face%EdgeCut(ecut)=Face%Edge2%Vert1%nrInP
      Face%ec=1
      EdgeList(1,ecut)=0
      EdgeList(2,ecut)=Face%Edge2%Vert1%nrInP
    END IF
  END IF
  !..............................................
  IF (Face%Edge3%Vert2%in_out==-1) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge3%Vert2%nrInP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge3%yes_sp==1) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge3%VertS%nrInP
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge3%VertS%nrInP
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%Edge3%VertS%nrInP
  ELSE IF (Face%Edge3%Vert2%in_out==0) THEN
    IF(Face%Edge3%Vert1%in_out==1.and.Face%Edge2%Vert1%in_out==1) THEN
      !Point ohne Verbindung im Face
    ELSE
      i=i+1
      Face%VertexList(i)=Face%Edge3%Vert2%nrInP
      Face%NumberVert=Face%NumberVert+1
      ecut=ecut+1
      Face%EdgeCut(ecut)=Face%Edge3%Vert2%nrInP
      Face%ec=1
      EdgeList(1,ecut)=0
      EdgeList(2,ecut)=Face%Edge3%Vert2%nrInP
    END IF
  END IF
  !................................................
  IF (Face%Edge4%Vert2%in_out==-1) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge4%Vert2%nrInP
    Face%NumberVert=Face%NumberVert+1
  END IF
  IF (Face%Edge4%yes_sp==1) THEN
    i=i+1
    Face%VertexList(i)=Face%Edge4%VertS%nrInP
    Face%NumberVert=Face%NumberVert+1
    ecut=ecut+1
    Face%EdgeCut(ecut)=Face%Edge4%VertS%nrInP
    Face%ec=1
    EdgeList(1,ecut)=1
    EdgeList(2,ecut)=Face%Edge4%VertS%nrInP
  ELSE IF (Face%Edge4%Vert2%in_out==0) THEN
    IF(Face%Edge4%Vert1%in_out==1.and.Face%Edge3%Vert2%in_out==1) THEN
      !Point ohne Verbindung im Face
    ELSE
      i=i+1
      Face%VertexList(i)=Face%Edge4%Vert2%nrInP
      Face%NumberVert=Face%NumberVert+1
      ecut=ecut+1
      Face%EdgeCut(ecut)=Face%Edge4%Vert2%nrInP
      Face%ec=1
      EdgeList(1,ecut)=0
      EdgeList(2,ecut)=Face%Edge4%Vert2%nrInP
    END IF
  END IF
  !Cut neu sortieren wenn ecut>2 ermittelt
  !Verts mit in_out==0 am Vert1||Vert2 filtern
  !richtiges ecut setzen, 
  !(aber dennoch Frage/offen, '1 x yes_sp und mehrfach Vert!%in_out==0'
  !je nachdem letzten 2 erfassten Cut in Liste? )
  IF(ecut>2) THEN
     j=0
     DO i=1,ecut
       IF(EdgeList(1,i)==1) THEN
         j=j+1
       END IF
     END DO 
     IF(j==2) THEN
       v=0
       Face%ec=1
       DO i=1,ecut
         IF(EdgeList(1,i)==1) THEN
           v=v+1
           Face%EdgeCut(v)=EdgeList(2,i)
         END IF
       END DO
     END IF
  END IF
 
  ! -------
  ! Special
  ! -------
  IF (Face%Edge1%Vert1%in_out<0 .AND. Face%Edge3%Vert2%in_out<0 .AND. &
     Face%Edge1%Vert2%in_out>0 .AND. Face%Edge3%Vert1%in_out>0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  ELSE IF (Face%Edge1%Vert1%in_out>0 .AND. Face%Edge3%Vert2%in_out>0 .AND. &
     Face%Edge1%Vert2%in_out<0 .AND. Face%Edge3%Vert1%in_out<0) THEN
     Face%ec=-1
     Face%EdgeCut(:)=0
  END IF
  !
  ! Face mit 1 Schnittpunkt -->EdgeCut ausblenden
  IF (ecut==1) THEN     ! zu Allgenein !?
    Face%ec=-1
  END IF

END SUBROUTINE SortVertFaceIn 


SUBROUTINE SortVertCutCell(Cell,ix,iy,iz)
  TYPE(Cell_T) :: Cell
  INTEGER :: ix,iy,iz
  INTEGER :: ListCut(1:2,1:8) !ListCut(1:2,1:6)
  INTEGER :: nCut,nCutNeu,i,j,iv
  INTEGER :: VertexFace(1:6,1:9)
  LOGICAL :: Cut
  INTEGER :: iP1,iP2,jP1,jP2,iVert,jVert
  INTEGER :: vc_out1
  TYPE(Face_T) :: Face

  IF(Cell%vc>0) THEN
     VertexFace=-1
     nCut=0
     ListCut=0
     Cell%in_out=Cell%Face1%in_out+Cell%Face2%in_out

     IF (Cell%Face1%NumberVert>2) THEN
       VertexFace(1,1:Cell%Face1%NumberVert)=Cell%Face1%VertexList(1:Cell%Face1%NumberVert)
       VertexFace(1,Cell%Face1%NumberVert+1)=VertexFace(1,1)
     END IF
     IF (Cell%Face2%NumberVert>2) THEN
       VertexFace(2,1:Cell%Face2%NumberVert)=Cell%Face2%VertexList(1:Cell%Face2%NumberVert)
       VertexFace(2,Cell%Face2%NumberVert+1)=VertexFace(2,1)
     END IF
     IF (Cell%Face3%NumberVert>2) THEN
       VertexFace(3,1:Cell%Face3%NumberVert)=Cell%Face3%VertexList(1:Cell%Face3%NumberVert)
       VertexFace(3,Cell%Face3%NumberVert+1)=VertexFace(3,1)
     END IF
     IF (Cell%Face4%NumberVert>2) THEN
       VertexFace(4,1:Cell%Face4%NumberVert)=Cell%Face4%VertexList(1:Cell%Face4%NumberVert)
       VertexFace(4,Cell%Face4%NumberVert+1)=VertexFace(4,1)
     END IF
     IF (Cell%Face5%NumberVert>2) THEN
       VertexFace(5,1:Cell%Face5%NumberVert)=Cell%Face5%VertexList(1:Cell%Face5%NumberVert)
       VertexFace(5,Cell%Face5%NumberVert+1)=VertexFace(5,1)
     END IF
     IF (Cell%Face6%NumberVert>2) THEN
       VertexFace(6,1:Cell%Face6%NumberVert)=Cell%Face6%VertexList(1:Cell%Face6%NumberVert)
       VertexFace(6,Cell%Face6%NumberVert+1)=VertexFace(6,1)
     END IF

     C1:DO i=1,6
       C2:DO iVert=1,8
         iP1=VertexFace(i,iVert)
         iP2=VertexFace(i,iVert+1)
         IF (ip2==-1) THEN
           EXIT C2
         END IF
         Cut=.TRUE.
         C3:DO j=1,6
           IF (j/=i) THEN
             C4:DO jVert=1,8
               jP1=VertexFace(j,jVert)
               jP2=VertexFace(j,jVert+1)
               IF (jp2==-1) THEN
                 EXIT C4
               END IF
               IF ((iP1==jP1.AND.iP2==jP2).OR.  &
                   (iP1==jP2.AND.iP2==jP1)) THEN
                 Cut=.FALSE.
                 EXIT C3
               END IF
             END DO C4
           END IF
         END DO C3
         IF (Cut) THEN
           ncut=ncut+1
           ListCut(1,nCut)=iP1
           ListCut(2,nCut)=iP2
         END IF
       END DO C2
     END DO C1

     IF (nCut>0) THEN
       nCutNeu=nCut
       DO i=1,ncut-1
         DO j=i+1,ncut
           IF (( (ListCut(1,i)==ListCut(1,j).AND.ListCut(2,i)==ListCut(2,j)).OR. &
                 (ListCut(1,i)==ListCut(2,j).AND.ListCut(2,i)==ListCut(1,j)) ) &
              .AND.ListCut(1,i)>0) THEN
             ListCut(:,j)=0
             nCutNeu=nCutNeu-1
   !         iTemP=ListCut(1,j)
   !         ListCut(1,j)=ListCut(2,j)
   !         ListCut(2,j)=iTemp
           END IF
         END DO
       END DO
       Cell%vc=nCutNeu
       IF (.NOT.ASSOCIATED(Cell%VertCut)) THEN
          ALLOCATE(Cell%VertCut(nCut))
          !Write(*,*) "Allocate new VertCut"
       ELSE
          DEALLOCATE(Cell%VertCut)
          ALLOCATE(Cell%VertCut(nCut))
       END IF
       Cell%VertCut(1)= ListCut(1,1)
       Cell%VertCut(2)= ListCut(2,1)
       ListCut(1,1)=0
       ListCut(2,1)=0
       iv=2
       S2:DO
         S1:DO i=1,nCut
           IF (Cell%VertCut(iv)==ListCut(1,i)) THEN
             iv=iv+1
             Cell%VertCut(iv)=ListCut(2,i)
             ListCut(:,i)=0
             EXIT S1
           ELSE IF (Cell%VertCut(iv)==ListCut(2,i)) THEN
             iv=iv+1
             Cell%VertCut(iv)=ListCut(1,i)
             ListCut(:,i)=0
             EXIT S1
           END IF
         END DO S1
         IF (iv>=nCutNeu) THEN
           EXIT S2
         END IF
       END DO S2
     ELSE  !ncut==0
       Cell%vc=0
       !IF (Cell%in_out==-8) THEN
       !  Cell%Vol=0.0d0
       !END IF
     END IF
  END IF

!Oro  ELSE  !vc==0
!Oro    IF (Cell%Face1%in_out==0.AND.(Cell%Face2%in_out>=2.AND.Cell%Face2%in_out<=4)) THEN 
!Oro       !Grenzflaeche Face1,und(Face2%in_out==4 od.
!Oro       !Grenzflaeche Face1,und Edge-Grenze , Grenzflaeche Face1,und Point-Grenze
!Oro       IF (.NOT.ASSOCIATED(Cell%VertCut)) THEN
!Oro          ALLOCATE(Cell%VertCut(Cell%Face1%NumberVert))
!Oro       END IF
!Oro       DO iv=1,Cell%Face1%NumberVert
!Oro         Cell%VertCut(iv)=Cell%Face1%VertexList(iv)
!Oro       END DO
!Oro    ELSE  !wenn Face 3-6  in_out==4
!Oro      IF (Cell%Face3%in_out==0) THEN
!Oro        IF (.NOT.ASSOCIATED(Cell%VertCut)) THEN
!Oro           ALLOCATE(Cell%VertCut(Cell%Face3%NumberVert))
!Oro        END IF
!Oro        DO iv=1,Cell%Face3%NumberVert
!Oro          Cell%VertCut(iv)=Cell%Face3%VertexList(iv)
!Oro        END DO
!Oro      ELSE IF (Cell%Face4%in_out==0) THEN
!Oro        IF (.NOT.ASSOCIATED(Cell%VertCut)) THEN
!Oro           ALLOCATE(Cell%VertCut(Cell%Face4%NumberVert))
!Oro        END IF
!Oro        DO iv=1,Cell%Face4%NumberVert
!Oro          Cell%VertCut(iv)=Cell%Face4%VertexList(iv)
!Oro        END DO
!Oro      ELSE IF (Cell%Face5%in_out==0) THEN
!Oro        IF (.NOT.ASSOCIATED(Cell%VertCut)) THEN
!Oro           ALLOCATE(Cell%VertCut(Cell%Face5%NumberVert))
!Oro        END IF
!Oro        DO iv=1,Cell%Face5%NumberVert
!Oro          Cell%VertCut(iv)=Cell%Face5%VertexList(iv)
!Oro        END DO
!Oro      ELSE IF (Cell%Face6%in_out==0) THEN
!Oro        IF (.NOT.ASSOCIATED(Cell%VertCut)) THEN
!Oro           ALLOCATE(Cell%VertCut(Cell%Face6%NumberVert))
!Oro        END IF
!Oro        DO iv=1,Cell%Face6%NumberVert
!Oro          Cell%VertCut(iv)=Cell%Face6%VertexList(iv)
!Oro        END DO
!Oro      END IF
!Oro    END IF
!Oro  END IF  ! If(vc>0) else 
!Oro  !Zaehlung Cut-Cellen (alle geschnittenen und GrenzCellen am Berg)
!Oro  IF (Cell%in_out>-8.AND.Cell%in_out<8 )THEN
!Oro    nr_cutinside=nr_cutinside+1
!Oro  END IF
!Oro
  !Specialfaelle {Cell%vc} ruecksetzen, wenn Schnitt nur Cell-Kante Streift
  !-----------------------------------  Betreff Abfrage Write general .or. hex
  !Cellen mit a) 1*(Edge%Vert1%in_out==0 und Edge%Vert2%in_out==0
  !           b) 1*Vertex an Grenzflaeche
  IF (Cell%vc==1) THEN
    Cell%vc=0
  END IF

!  !Celle mit Grenzflaeche Face2, Cell%Vol =0 setzen
!  !wird zuvor, da Face%ec -->Cell%vc ->Vol berechnet 
!  IF (Cell%Face2%Edge1%Vert1%in_out==0 .AND. Cell%Face2%Edge1%Vert2%in_out==0 .AND. &
!      Cell%Face2%Edge3%Vert1%in_out==0 .AND. Cell%Face2%Edge3%Vert2%in_out==0 .AND. &
!      Cell%Face1%in_out==-4) THEN
!      Cell%vc=0
!  END IF
END SUBROUTINE SortVertCutCell


SUBROUTINE SortVertAllFacesIn

  INTEGER :: ib,i,j,k
  INTEGER :: in_out
  !CALL Display_SortVertAllFacesIn
  DO ib=1,nb
    CALL Set(Floor(ib))

    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
            CALL SortVertFaceIn(Faces_XY(i,j,k)%Face)
          END IF
        END DO
      END DO
    END DO

    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
            CALL SortVertFaceIn(Faces_YZ(i,j,k)%Face)
          END IF
        END DO
      END DO
    END DO

    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
            CALL SortVertFaceIn(Faces_ZX(i,j,k)%Face)
          END IF
        END DO
      END DO
    END DO

  END DO   ! ib
END SUBROUTINE SortVertAllFacesIn


SUBROUTINE SortVertCutAllCells
  INTEGER :: ib,i,j,k
  INTEGER :: v_x,v_y,v_z
  INTEGER :: in_out,in_out_view
  !CALL Display_SortVertCutAllCells

  nr_cutinside=0
  DO ib=1,nb
    CALL Set(Floor(ib))
    !Write(*,*) "         ...SortVertCutAllCells inside Block : ", ib
    !Write(*,*) "               inside Block : ", ib 
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
           IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
             !#ifdef _PANALYZE
                 !WRITE(*,*) "Block : ",ib,"  SortVertCutCell :  i=",i," j=",j," k=",k
             !#endif
             CALL SortVertCutCell(Cell(i,j,k)%Cell,i,j,k)
           END IF
        END DO  ! k
      END DO   ! j
    END DO   ! i
  END DO  !ib
 !WRITE(*,*) "nr_cutinside =", nr_cutinside
END SUBROUTINE SortVertCutAllCells


SUBROUTINE Set_Domain_OutView()
INTEGER :: ix,iy,iz

TYPE OutMv_XYZ
  CHARACTER(3) :: who           ! ixa,ixe,iya,iye,iza,ize
  CHARACTER(5) :: scaling       ! upper,lower
  INTEGER      :: yes           ! 1 .or. 0
END TYPE
TYPE(OutMv_XYZ) :: zg(1:6)

  IF(out_area=='y') THEN 
  IF (TRIM(Domain%def_out_domain)==def_out_coord) THEN

  ! search between Domain%ix0,Domain%ix1
  !.....................................
  ix=Domain%ix0    ! lower ix
  DO WHILE (Domain%view_xa>Domain%xP(ix).and.ix<Domain%ix1)
     ix=ix+1
  END DO
  IF(Domain%view_xa==Domain%xP(ix)) THEN
     Domain%view_ixa=ix
  ELSE
     !Differenz-Check xa/xP(ix) und skalieren
     IF((Domain%xP(ix)-Domain%view_xa)<=(Domain%view_xa-Domain%xP(ix-1))) THEN
        Domain%view_ixa=ix
        !zg(1)=(/'ixa','upper',1/)
        !Compiler : All elements in an array constructor must have the same type and type parameters
        zg(1)%who="ixa";zg(1)%scaling="upper";zg(1)%yes=1
     ELSE
        Domain%view_ixa=ix-1
        !zg(1)=(/"ixa","lower",1/)
        zg(1)%who="ixa";zg(1)%scaling="lower";zg(1)%yes=1
     END IF
     Write(*,*) "                                     ",&
                &"View-OutGmV a must resized:"
     Write(*,*) "                                     ",&
                &"Domain%view_ixa next",zg(1)%scaling,"cell border!"
  END IF
  ix=Domain%ix1  !upper ix
  DO WHILE (Domain%view_xe<Domain%xP(ix).and.ix>Domain%ix0)
     ix=ix-1
  END DO
  IF(Domain%view_xe==Domain%xP(ix)) THEN
     Domain%view_ixe=ix
  ELSE
     !Differenz-Check xe/xP(ix) und skalieren
     IF((Domain%view_xe-Domain%xP(ix)) <= (Domain%xP(ix+1)-Domain%view_xe)) THEN
        Domain%view_ixe=ix
        !zg(2)=(/"ixe","lower",1/)
        zg(2)%who="ixe";zg(2)%scaling="lower";zg(2)%yes=1
     ELSE
        Domain%view_ixe=ix+1
        !zg(2)=(/"ixe","upper",1/)
        zg(2)%who="ixe";zg(2)%scaling="upper";zg(2)%yes=1
     END IF
     if(zg(1)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_ixe next",zg(2)%scaling,"cell border!"
  END IF

  ! search between Domain%iy0,Domain%iy1
  !.....................................
  iy=Domain%iy0    ! lower iy
  DO WHILE (Domain%view_ya>Domain%yP(iy).and.iy<Domain%iy1)
     iy=iy+1
  END DO
  IF(Domain%view_ya==Domain%yP(iy)) THEN
     Domain%view_iya=iy
  ELSE
     !Differenz-Check ya/yP(iy) und skalieren
     IF((Domain%yP(iy)-Domain%view_ya) <=(Domain%view_ya-Domain%yP(iy-1))) THEN
        Domain%view_iya=iy
        !zg(3)=(/"iya","upper",1/)
        zg(3)%who="iya";zg(3)%scaling="upper";zg(3)%yes=1
     ELSE
        Domain%view_iya=iy-1
        !zg(3)=(/"iya","lower",1/)
        zg(3)%who="iya";zg(3)%scaling="lower";zg(3)%yes=1
     END IF
     IF(zg(1)%yes==0.eqv.zg(2)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_iya next",zg(3)%scaling,"cell border!"
  END IF
  iy=Domain%iy1   ! upper iy
  DO WHILE (Domain%view_ye<Domain%yP(iy).and.iy>Domain%iy0)
     iy=iy-1
  END DO
  IF(Domain%view_ye==Domain%yP(iy)) THEN
     Domain%view_iye=iy
  ELSE
     !Differenz-Check ye/yP(iy) und skalieren
     IF((Domain%view_ye-Domain%yP(iy)) <= (Domain%yP(iy+1)-Domain%view_ye)) THEN
        Domain%view_iye=iy
        !zg(4)=(/"iye","lower",1/)
        zg(4)%who="iye";zg(4)%scaling="lower";zg(4)%yes=1
     ELSE
        Domain%view_iye=iy+1
        !zg(4)=(/"iye","upper,1"/)
        zg(4)%who="iye";zg(4)%scaling="upper";zg(4)%yes=1
     END IF
     if(zg(1)%yes==0.eqv.zg(2)%yes==0.eqv.zg(3)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_iye next",zg(4)%scaling,"cell border!"
  END IF

  ! search between Domain%iz0,Domain%iz1
  !.....................................
  iz=Domain%iz0   ! lower iz
  DO WHILE (Domain%view_za>Domain%zP(iz).and.iz<Domain%iz1)
     iz=iz+1
  END DO
  IF(Domain%view_za==Domain%zP(iz)) THEN
     Domain%view_iza=iz
  ELSE
     !Differenz-Check za/zP(zy) und skalieren
     IF((Domain%zP(iz)-Domain%view_za) <=(Domain%view_za-Domain%zP(iz-1))) THEN
        Domain%view_iza=iz
        !zg(5)=(/"iza","upper",1/)
        zg(5)%who="iza";zg(5)%scaling="upper";zg(5)%yes=1
     ELSE
        Domain%view_iza=iz-1
        !zg(5)=(/"iza","lower",1/)
        zg(5)%who="iza";zg(5)%scaling="lower";zg(5)%yes=1
     END IF
     if(zg(1)%yes==0.eqv.zg(2)%yes==0.eqv.zg(3)%yes==0.eqv.zg(4)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_iza next",zg(5)%scaling,"cell border!"
  END IF
  iz=Domain%iz1    ! upper iz
  DO WHILE (Domain%view_ze<Domain%zP(iz).and.iz>Domain%iz0)
     iz=iz-1
  END DO
  IF(Domain%view_ze==Domain%zP(iz)) THEN
     Domain%view_ize=iz
  ELSE
     !Differenz-Check ze/zP(iz) und skalieren
     IF((Domain%view_ze-Domain%zP(iz)) <= (Domain%zP(iz+1)-Domain%view_ze)) THEN
        Domain%view_ize=iz
        !zg(6)=(/"ize","lower",1/)
        zg(6)%who="ize";zg(6)%scaling="lower";zg(6)%yes=1
     ELSE
        Domain%view_ize=iz+1
        !zg(6)=(/"ize","upper",1/)
        zg(6)%who="ize";zg(6)%scaling="upper";zg(6)%yes=1
     END IF
     if(zg(1)%yes==0.eqv.zg(2)%yes==0.eqv.zg(3)%yes==0.eqv.zg(4)%yes==0.eqv.zg(5)%yes==0) THEN
       Write(*,*) "                                     ",&
                  &"View-OutGmV a must resized:"
     END IF
     Write(*,*) "                                     ",&
                &"Domain%view_iza next",zg(6)%scaling,"cell border!"
  END IF

  END IF  !IF (...def_out_coord)
  END IF  !IF(out_area=='y')
END SUBROUTINE Set_Domain_OutView


SUBROUTINE OutputCheckCellen
 Write(*,*)
 Write(*,*) "------------------------------------------------"
 Write(*,*) "| OutputCheckCellen()"
 Write(*,*) "| !Check aus WriteCellOroGMVAscii()"
 Write(*,*) "| nr_gen       = ",nr_gen
 Write(*,*) "| nr_grenzeF2  = ",nr_grenzeF2
 Write(*,*) "| nr_grenzeF1  = ",nr_grenzeF1
 Write(*,*) "| nr_insidehex = ",nr_insidehex
 Write(*,*) "| ----------------------------------------------"
 Write(*,*) "| !Check aus AnalyzeAllCells()"
 Write(*,*) "|                    ", "---Global---", "   ---View--- "
 Write(*,'(a22,i10,i12)') "| nr_cells         = ",nr_cells,nr_viewcells
 Write(*,'(a22,i10,i12)') "| NrAll_Cells      = ",NrAll_Cells
 Write(*,'(a22,i10,i12)') "| NrW_Cells        = ",NrW_Cells
 Write(*,'(a22,i10,i12)') "| NrW_Cells1       = ",NrW_Cells1
 Write(*,'(a22,i10,i12)') "| NrW_Cells2       = ",NrW_Cells2
 Write(*,'(a22,i10,i12)') "| NrB_Cells        = ",NrB_Cells
 Write(*,'(a22,i10,i12)') "| nr_cutcells      = ",nr_cutcells,nr_viewcutcells
 Write(*,'(a22,i10,i12)') "| nr_incells       = ",nr_incells,nr_viewincells
 Write(*,'(a22,i10,i12)') "| nr_orocells      = ",nr_orocells,nr_vieworocells
 Write(*,'(a22,i10,i12)') "| nr_cutplanecells = ",nr_cutplanecells,nr_viewcutplanecells
 Write(*,*) "| ----------------------------------------------"
 Write(*,*) "| !Check aus SortVertCutCell()"
 Write(*,*) "| nr_cutinside    = ",nr_cutinside
 Write(*,*) "------------------------------------------------"
 Write(*,*)
END SUBROUTINE OutputCheckCellen

SUBROUTINE WriteFaceEVProt(Face,i,j,k)
  TYPE(Face_T), POINTER :: Face
  INTEGER :: i,j,k
  WRITE(OutUnitProt,*) ' Face%--->',' \/',i,'\/',j,'\/',k,'\/'
  WRITE(OutUnitProt,*) "    !(TYPE Vertex_T:", "  Point(x,y,z)    in_out   nrP   nrInP   nrCutP   Shift)"
  WRITE(OutUnitProt,*) "    Edge1%Vert1:",Face%Edge1%Vert1
  WRITE(OutUnitProt,*) "    Edge1%Vert2:",Face%Edge1%Vert2
  WRITE(OutUnitProt,*) "    Edge3%Vert1:",Face%Edge3%Vert1
  WRITE(OutUnitProt,*) "    Edge3%Vert2:",Face%Edge3%Vert2
END SUBROUTINE WriteFaceEVProt
 
SUBROUTINE WriteFaceProt(Face,i,j,k)
  TYPE(Face_T), POINTER :: Face
  INTEGER :: i,j,k
  WRITE(OutUnitProt,*) ' Face%--->',' \/',i,'\/',j,'\/',k,'\/'
  WRITE(OutUnitProt,*) "    !(TYPE Vertex_T:", "  Point(x,y,z)    in_out   nrP   nrInP   nrCutP   Shift)"
  WRITE(OutUnitProt,*) "    Edge1%Vert1:",Face%Edge1%Vert1
  WRITE(OutUnitProt,*) "    Edge1%Vert2:",Face%Edge1%Vert2
  WRITE(OutUnitProt,*) "    Edge3%Vert1:",Face%Edge3%Vert1
  WRITE(OutUnitProt,*) "    Edge3%Vert2:",Face%Edge3%Vert2
  WRITE(OutUnitProt,*) "    ..........."
  WRITE(OutUnitProt,*) '    Vol        =',Face%Vol
  WRITE(OutUnitProt,*) '    NumberVert =',Face%NumberVert
  IF(Face%NumberVert>1) THEN
  WRITE(OutUnitProt,*) '    VertexList =',Face%VertexList(1:Face%NumberVert)
  ELSE
  WRITE(OutUnitProt,*) '    VertexList ='," no"
  END IF
  WRITE(OutUnitProt,*) '    in_out     =',Face%in_out
  WRITE(OutUnitProt,*) '    ec         =',Face%ec
  WRITE(OutUnitProt,*) '    EdgeCut    =',Face%EdgeCut
END SUBROUTINE WriteFaceProt

SUBROUTINE WriteCellProt(Cell,i,j,k)

  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: F1V1,F1V2,F1V3,F1V4
  INTEGER :: F2V1,F2V2,F2V3,F2V4

  IF (ASSOCIATED(Cell)) THEN
    !.................................
    F1V1=Cell%Face1%Edge1%Vert1%in_out
    F1V2=Cell%Face1%Edge1%Vert2%in_out
    F1V3=Cell%Face1%Edge3%Vert2%in_out
    F1V4=Cell%Face1%Edge3%Vert1%in_out 
    !.................................
    F2V1=Cell%Face2%Edge1%Vert1%in_out
    F2V2=Cell%Face2%Edge1%Vert2%in_out
    F2V3=Cell%Face2%Edge3%Vert2%in_out
    F2V4=Cell%Face2%Edge3%Vert1%in_out 
    !.................................
    WRITE(OutUnitProt,*) '-------------- WriteCellProt() ---------------------------------'
    WRITE(OutUnitProt,*) "        ",F2V4,"-----------",F2V3
    WRITE(OutUnitProt,*) "        /|  F2      /|" 
    WRITE(OutUnitProt,*) "      ",F2V1,"-----------",F2V2," |"  
    WRITE(OutUnitProt,*) "       | ",F1V4,"---------|-",F1V3
    WRITE(OutUnitProt,*) "       |/   F1     |/" 
    WRITE(OutUnitProt,*) "      ",F1V1,"-----------",F1V2   
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) '* Cell--> ',' \/',i,'\/',j,'\/',k,'\/'
    WRITE(OutUnitProt,*) '    Vol       =',Cell%Vol
    WRITE(OutUnitProt,*) '    MidPoint  =',Cell%MidPoint
    WRITE(OutUnitProt,*) '    vc        =',Cell%vc
    WRITE(OutUnitProt,*) '    in_out    =',Cell%in_out
    IF (ASSOCIATED(Cell%VertCut)) THEN
    WRITE(OutUnitProt,*) '    EdgeCut   =',Cell%VertCut
    ELSE 
    WRITE(OutUnitProt,*) '    EdgeCut   ='," no"
    END IF
    WRITE(OutUnitProt,*) '    CutF_MidP =',Cell%CutF_MidP
    WRITE(OutUnitProt,*) '    Face1%Vol =',Cell%Face1%Vol
    WRITE(OutUnitProt,*) '    Face2%Vol =',Cell%Face2%Vol
    WRITE(OutUnitProt,*) '    Face3%Vol =',Cell%Face3%Vol
    WRITE(OutUnitProt,*) '    Face4%Vol =',Cell%Face4%Vol
    WRITE(OutUnitProt,*) '    Face5%Vol =',Cell%Face5%Vol
    WRITE(OutUnitProt,*) '    Face6%Vol =',Cell%Face6%Vol
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceXY','  (F1)'
    CALL WriteFaceProt(Cell%Face1,i,j,k-1)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceXY','  (F2)'
    CALL WriteFaceProt(Cell%Face2,i,j,k)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceXZ','  (F3)'
    CALL WriteFaceProt(Cell%Face3,i,j-1,k)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceXZ','  (F4)'
    CALL WriteFaceProt(Cell%Face4,i,j,k)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceYZ','  (F5)'
    CALL WriteFaceProt(Cell%Face5,i-1,j,k)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceYZ','  (F6)'
    CALL WriteFaceProt(Cell%Face6,i,j,k)
    WRITE(OutUnitProt,*) 'Cell','(',i,',',j,',',k,')'
    WRITE(OutUnitProt,*) '-------------- WriteCellProt() Ende ----------------------------'
  END IF

END SUBROUTINE WriteCellProt

END MODULE Grid_Mod
