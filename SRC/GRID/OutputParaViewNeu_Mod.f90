MODULE OutputParaView_Mod

! USE F_Mod
  USE Floor_Mod
  USE IOControl_Mod
  USE Parametric_Mod
  !USE Grid_Mod !, ONLY: &
               ! VolFace_XY,VolFace_YZ,VolFace_ZX
  USE OutputOutGmvGNeu_Mod

  IMPLICIT NONE

CONTAINS

SUBROUTINE WriteTropoCellsParaV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k,n,ic,iec,is_eg,in_out
  INTEGER :: ecnr,exnr,eynr,eznr,enr,ec_set
  INTEGER :: fcnr,fxynr,fzxnr,fyznr,fnr
  INTEGER :: cnr,cnr_cut,ifn,nfn
  INTEGER :: lvertout
  INTEGER :: eg_listfxy(1:8)=0
  INTEGER :: eg_listfzx(1:8)=0
  INTEGER :: eg_listfyz(1:8)=0
  CHARACTER(2),PARAMETER :: zw_leerzei='  '
  CHARACTER(8),PARAMETER :: leerzei_acht='        '
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: liste_ec
  INTEGER :: nr_cells_write,nr_faces_write,nr_edges_write

  !CALL SetTropoCellsParaV(FileName)
  !CALL SetTropoOroCellsParaV(FileName)
  CALL WriteAuswOutParaVToProt()

  !WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.pva.tropo',STATUS='UNKNOWN')
  CALL Display_OutTroposParaV(FileName)

  WRITE(10,*) "Tropo  ",TRIM(FileName)
  Write(10,*) nb, "   ! Number of blocks"

  !WRITE(10,'(a8,i8)') nodes,nr_out
  WRITE(10,*) "Nr_Nodes",nr_out
  lvertout=SIZE(VertOut)
  If(lvertout/=nr_out) THEN
    WRITE(*,*) "-------- Warnung! ----------"
    WRITE(*,*) "Grösse VertOut = ", lvertout
    WRITE(*,*) "        nr_out = ", nr_out
    WRITE(*,*) "----------------------------"
  END IF
  !!VertOut(i)%Point%x,VertOut(i)%Point%y,VertOut(i)%Point%z
  DO i=1,nr_out
     WRITE(10,*)"n", i, zw_leerzei, xParametricOut(VertOut(i)), &
           &    yParametricOut(VertOut(i)),zParametricOut(VertOut(i))
  END DO

  !WRITE(*,*) ">....... Topo Cells for ParaView"
  !WRITE(10,'(a8,i8)') cells,nr_viewcells
  DO ib=1,nb
    Write(*,*) "                                     Block : ",ib,"\/",nb
    CALL Set(Floor(ib))

    Write(10,*)"............................................"
    Write(10,*)ib, "   ! block: ",ib,"/",nb
    Write(10,*)"............................................"

    !-- Output---EdgeX/EdgeY/EdgeZ/EdgeCut
    Write(10,*)"Nr_Edges",stopo_enr
    !-- Output--Edge_X
    nr_edges_write=0
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
             IF(MAX(Edges_X(i,j,k)%Edge%Vert1%in_out,Edges_X(i,j,k)%Edge%Vert2%in_out)>0 &
               .OR. &
               (Edges_X(i,j,k)%Edge%Vert1%in_out==Edges_X(i,j,k)%Edge%Vert2%in_out.AND. &
                Edges_X(i,j,k)%Edge%Vert1%in_out==0)) THEN
               IF(Edges_X(i,j,k)%Edge%yes_sp==1) THEN
                  IF(Edges_X(i,j,k)%Edge%Vert1%nrP==0) THEN
                    Write(10,*) "ex", Topo(i,j,k)%Topo%exnr, zw_leerzei, &
                     & Edges_X(i,j,k)%Edge%VertS%nrP,Edges_X(i,j,k)%Edge%Vert2%nrP
                    nr_edges_write=nr_edges_write+1
                  ELSE IF(Edges_X(i,j,k)%Edge%Vert2%nrP==0) THEN
                    Write(10,*) "ex", Topo(i,j,k)%Topo%exnr, zw_leerzei, &
                      & Edges_X(i,j,k)%Edge%Vert1%nrP,Edges_X(i,j,k)%Edge%VertS%nrP 
                    nr_edges_write=nr_edges_write+1
                  END IF
               ELSE
                 Write(10,*) "ex", Topo(i,j,k)%Topo%exnr, zw_leerzei, &
                   & Edges_X(i,j,k)%Edge%Vert1%nrP,Edges_X(i,j,k)%Edge%Vert2%nrP 
                 nr_edges_write=nr_edges_write+1
               END IF
             END IF
          ELSE
            IF(MAX(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)>0 &
               .OR. &
               (Vertices(i-1,j,k)%in_out==Vertices(i,j,k)%in_out .AND. &
                Vertices(i,j,k)%in_out==0)) THEN
              Write(10,*) "ex", Topo(i,j,k)%Topo%exnr, zw_leerzei, &
                & Vertices(i-1,j,k)%nrP,Vertices(i,j,k)%nrP
              nr_edges_write=nr_edges_write+1
            END IF 
          END IF
        END DO
      END DO
    END DO
    !-- Output--Edge_Y
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
            IF(MAX(Edges_Y(i,j,k)%Edge%Vert1%in_out,Edges_Y(i,j,k)%Edge%Vert2%in_out)>0 &
               .OR. &
               (Edges_Y(i,j,k)%Edge%Vert1%in_out==Edges_Y(i,j,k)%Edge%Vert2%in_out.AND. &
                Edges_Y(i,j,k)%Edge%Vert1%in_out==0)) THEN
               IF(Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
                 IF(Edges_Y(i,j,k)%Edge%Vert1%nrP==0) THEN
                   Write(10,*) "ey", Topo(i,j,k)%Topo%eynr, zw_leerzei, &
                    & Edges_Y(i,j,k)%Edge%VertS%nrP,Edges_Y(i,j,k)%Edge%Vert2%nrP
                   nr_edges_write=nr_edges_write+1
                 ELSE IF(Edges_Y(i,j,k)%Edge%Vert2%nrP==0) THEN
                   Write(10,*) "ey", Topo(i,j,k)%Topo%eynr, zw_leerzei, &
                    & Edges_Y(i,j,k)%Edge%Vert1%nrP,Edges_Y(i,j,k)%Edge%VertS%nrP
                   nr_edges_write=nr_edges_write+1
                 END IF
               ELSE
                  Write(10,*) "ey", Topo(i,j,k)%Topo%eynr, zw_leerzei, &
                   & Edges_Y(i,j,k)%Edge%Vert1%nrP,Edges_Y(i,j,k)%Edge%Vert2%nrP
                  nr_edges_write=nr_edges_write+1
               END IF
            END IF
          ELSE
            IF(MAX(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)>0 &
               .OR. &
               (Vertices(i,j-1,k)%in_out==Vertices(i,j,k)%in_out .AND. &
                Vertices(i,j,k)%in_out==0)) THEN
               Write(10,*) "ey", Topo(i,j,k)%Topo%eynr, zw_leerzei, &
                 & Vertices(i,j-1,k)%nrP,Vertices(i,j,k)%nrP
               nr_edges_write=nr_edges_write+1
            END IF
          END IF
        END DO
      END DO
    END DO
    !-- Output--Edge_Z
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
            IF(MAX(Edges_Z(i,j,k)%Edge%Vert1%in_out,Edges_Z(i,j,k)%Edge%Vert2%in_out)>0 &
               .OR. &
               (Edges_Z(i,j,k)%Edge%Vert1%in_out==Edges_Z(i,j,k)%Edge%Vert2%in_out.AND. &
                Edges_Z(i,j,k)%Edge%Vert1%in_out==0)) THEN
               IF(Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
                 IF(Edges_Z(i,j,k)%Edge%Vert1%nrP==0) THEN
                   Write(10,*) "ez", Topo(i,j,k)%Topo%eznr, zw_leerzei, &
                     & Edges_Z(i,j,k)%Edge%VertS%nrP,Edges_Z(i,j,k)%Edge%Vert2%nrP
                   nr_edges_write=nr_edges_write+1
                 ELSE IF(Edges_Z(i,j,k)%Edge%Vert2%nrP==0) THEN
                   Write(10,*) "ez", Topo(i,j,k)%Topo%eznr, zw_leerzei, &
                     & Edges_Z(i,j,k)%Edge%Vert1%nrP,Edges_Z(i,j,k)%Edge%VertS%nrP
                   nr_edges_write=nr_edges_write+1
                 END IF
               ELSE
                 Write(10,*) "ez", Topo(i,j,k)%Topo%eznr, zw_leerzei, &
                   & Edges_Z(i,j,k)%Edge%Vert1%nrP,Edges_Z(i,j,k)%Edge%Vert2%nrP
                 nr_edges_write=nr_edges_write+1
               END IF
            END IF
          ELSE
            IF(MAX(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)>0 &
               .OR. &
               (Vertices(i,j,k-1)%in_out==Vertices(i,j,k)%in_out.AND. &
                Vertices(i,j,k)%in_out==0)) THEN
               Write(10,*) "ez", Topo(i,j,k)%Topo%eznr, zw_leerzei, &
                 & Vertices(i,j,k-1)%nrP,Vertices(i,j,k)%nrP
               nr_edges_write=nr_edges_write+1
            END IF
          END IF
        END DO
      END DO
    END DO
    IF(stopo_ecnr>0) THEN
    ALLOCATE(liste_ec(topo_pos_ec1:topo_pos_ecn,1:6))
    DO i=ix0+1,ix1     !EdgeCut-Sortieren-für Ausgabe
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
            IF(Faces_XY(i,j,k)%Face%egc_nr>0) THEN
               !Write(10,*) "               aus Faces_XY--------",i,j,k
               !Write(10,*) "ec",Faces_XY(i,j,k)%Face%egc_nr, zw_leerzei, &
               !      & Faces_XY(i,j,k)%Face%EdgeCut(1),Faces_XY(i,j,k)%Face%EdgeCut(2)
               liste_ec(Topo(i,j,k)%Topo%ecnrxy,1)=Topo(i,j,k)%Topo%ecnrxy
               liste_ec(Topo(i,j,k)%Topo%ecnrxy,2)= &
                         &  Faces_XY(i,j,k)%Face%EdgeCut(1) 
               liste_ec(Topo(i,j,k)%Topo%ecnrxy,3)= &
                         &  Faces_XY(i,j,k)%Face%EdgeCut(2)
               liste_ec(Topo(i,j,k)%Topo%ecnrxy,4)=i
               liste_ec(Topo(i,j,k)%Topo%ecnrxy,5)=j
               liste_ec(Topo(i,j,k)%Topo%ecnrxy,6)=k
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
            IF(Faces_YZ(i,j,k)%Face%egc_nr>0) THEN
               !Write(10,*) "               aus Faces_YZ--------",i,j,k
               !Write(10,*) "ec",Faces_YZ(i,j,k)%Face%egc_nr, zw_leerzei, &
               !      & Faces_YZ(i,j,k)%Face%EdgeCut(1),Faces_YZ(i,j,k)%Face%EdgeCut(2)
               liste_ec(Topo(i,j,k)%Topo%ecnryz,1)=Topo(i,j,k)%Topo%ecnryz
               liste_ec(Topo(i,j,k)%Topo%ecnryz,2)=Faces_YZ(i,j,k)%Face%EdgeCut(1)
               liste_ec(Topo(i,j,k)%Topo%ecnryz,3)=Faces_YZ(i,j,k)%Face%EdgeCut(2)
               liste_ec(Topo(i,j,k)%Topo%ecnryz,4)=i
               liste_ec(Topo(i,j,k)%Topo%ecnryz,5)=j
               liste_ec(Topo(i,j,k)%Topo%ecnryz,6)=k
            END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
            IF(Faces_ZX(i,j,k)%Face%egc_nr>0) THEN
               !Write(10,*) "               aus Faces_ZX--------",i,j,k
               !Write(10,*) "ec",Faces_ZX(i,j,k)%Face%egc_nr, zw_leerzei, &
               !      & Faces_ZX(i,j,k)%Face%EdgeCut(1),Faces_ZX(i,j,k)%Face%EdgeCut(2)
               liste_ec(Topo(i,j,k)%Topo%ecnrzx,1)=Topo(i,j,k)%Topo%ecnrzx
               liste_ec(Topo(i,j,k)%Topo%ecnrzx,2)=Faces_ZX(i,j,k)%Face%EdgeCut(1)
               liste_ec(Topo(i,j,k)%Topo%ecnrzx,3)=Faces_ZX(i,j,k)%Face%EdgeCut(2)
               liste_ec(Topo(i,j,k)%Topo%ecnrzx,4)=i
               liste_ec(Topo(i,j,k)%Topo%ecnrzx,5)=j
               liste_ec(Topo(i,j,k)%Topo%ecnrzx,6)=k
            END IF
          END IF
        END DO
      END DO
    END DO
    !-- Output--EdgeCut
    DO i=topo_pos_ec1,topo_pos_ecn
       Write(10,*) "ec",liste_ec(i,1),zw_leerzei,liste_ec(i,2),liste_ec(i,3)
         !!&         liste_ec(i,4),liste_ec(i,5),liste_ec(i,6)
       nr_edges_write=nr_edges_write+1
    END DO
    DEALLOCATE(liste_ec)
    END IF
    !-- Output--FaceCut,Faces_XY,Faces_XZ,Faces_YZ ---------------------------------------
    !-- Output--FaceCut
    Write(10,*) "Nr_Faces", stopo_fnr 
    nr_faces_write=0
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
            IF(Cell(i,j,k)%Cell%vc>0) THEN
               !Write(10,*) "                  i =", i, " j = ",j, " k = ",k
               Write(10,*) "fc", Topo(i,j,k)%Topo%FCutNr, zw_leerzei,  &
                       &   Cell(i,j,k)%Cell%vc, zw_leerzei, tropos_border
               Write(10,*) leerzei_acht, &
                          & (Topo(i,j,k)%Topo%fceg_nrcut(iec),iec=1,Cell(i,j,k)%Cell%vc)
               nr_faces_write=nr_faces_write+1 
            END IF
          END IF
        END DO
      END DO
    END DO
    !-- Output--Faces_XY
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
          IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
            IF(Vertices(i-1,j-1,k)%in_out==0 .and. &
               Vertices(i  ,j-1,k)%in_out==0 .and. &
               Vertices(i-1,j  ,k)%in_out==0 .and. &
               Vertices(i  ,j  ,k)%in_out==0) THEN
               IF(k==iz0) THEN  ! Überarbeiten mehrere Blöcke 'igz0'
                 ! wenn Boden-Modell
                   Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, &
                            &   zw_leerzei, tropos_border 
                   Write(10,*) leerzei_acht, &
                           & Topo(i,j-1,k)%Topo%exnr, Topo(i,j,k)%Topo%eynr,  &
                           & Topo(i,j,k)%Topo%exnr, Topo(i-1,j,k)%Topo%eynr
                   nr_faces_write=nr_faces_write+1 
               ELSE
                 ! als Grenzfläche 'fc' ausgegeben
                 !Write(10,*) "Faces_XY()%Face als Grenzfläche fc ausgegeben (4xVert%in_out==0)",i,j,k
                 !Write(10,*) "                       FCutNr:", Topo(i,j,k)%Topo%FCutNr
               END IF
            ELSE IF(Faces_XY(i,j,k)%Face%in_out>-4.and.Faces_XY(i,j,k)%Face%numberVert>2) THEN
               !Reihenfolge:e1-e4, speciell je Face-Def.
               is_eg=0;ec_set=0
               IF(Faces_XY(i,j,k)%Face%ec==1) THEN
                  ! EdgeListe aufstellen
                 !.e1.....................................
                 IF(Topo(i,j-1,k)%Topo%exnr>0) THEN
                   is_eg=is_eg+1
                   eg_listfxy(is_eg)=Topo(i,j-1,k)%Topo%exnr
                   IF(Faces_XY(i,j,k)%Face%Edge1%yes_sp==1) THEN
                      IF(Faces_XY(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                         is_eg=is_eg+1
                         eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                         ec_set=1
                         !Write(10,*) i,j,k,", eg1 - Annulus Test F_XY"
                      END IF
                   END IF
                 ELSE 
                   IF(Faces_XY(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                      Faces_XY(i,j,k)%Face%Edge1%in_out==-1) THEN
                       is_eg=is_eg+1
                       eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                       ec_set=1
                   END IF
                 END IF
                 !.e2.....................................
                 IF(Topo(i,j,k)%Topo%eynr>0) THEN
                    is_eg=is_eg+1
                    eg_listfxy(is_eg)=Topo(i,j,k)%Topo%eynr
                    IF(Faces_XY(i,j,k)%Face%Edge2%yes_sp==1) THEN
                      IF(Faces_XY(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                          ec_set=1
                          !Write(10,*) i,j,k,",  eg2 - Annulus Test F_XY"
                        END IF
                      END IF
                    END IF
                 ELSE 
                    IF(Faces_XY(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                       Faces_XY(i,j,k)%Face%Edge2%in_out==-1) THEN
                      IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                        ec_set=1
                      END IF
                    END IF
                 END IF
                 !.e3.....................................
                 IF(Topo(i,j,k)%Topo%exnr>0) THEN
                   is_eg=is_eg+1
                   eg_listfxy(is_eg)=Topo(i,j,k)%Topo%exnr
                   IF(Faces_XY(i,j,k)%Face%Edge3%yes_sp==1) THEN
                     IF(Faces_XY(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                       IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                        ec_set=1
                        !Write(10,*) i,j,k,",  eg3 - Annulus Test F_XY"
                       END IF
                     END IF
                   END IF
                 ELSE
                   IF(Faces_XY(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                      Faces_XY(i,j,k)%Face%Edge3%in_out==-1) THEN
                     IF(ec_set==0) THEN
                       is_eg=is_eg+1
                       eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                       ec_set=1
                     END IF
                   END IF
                 END IF
                 !.e4.....................................
                 IF(Topo(i-1,j,k)%Topo%eynr>0) THEN
                   is_eg=is_eg+1
                   eg_listfxy(is_eg)=Topo(i-1,j,k)%Topo%eynr
                   IF(Faces_XY(i,j,k)%Face%Edge4%yes_sp==1) THEN
                     IF(Faces_XY(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                       IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                        ! Write(10,*) i,j,k,",  eg4 - Annulus Test F_XY"
                       END IF
                     END IF
                   END IF 
                 ELSE 
                    IF(Faces_XY(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                       Faces_XY(i,j,k)%Face%Edge4%in_out==-1) THEN
                      IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                      END IF
                    END IF
                 END IF
                 IF(is_eg==2 .and. Faces_XY(i,j,k)%Face%NumberVert>2) THEN
                    is_eg=is_eg+1
                    eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                 END IF
                 !........................................
                 !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==1"
                 IF(k==iz0.OR.k==iz1) THEN
                    Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, is_eg, &
                            &   zw_leerzei, tropos_border 
                 ELSE 
                    Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, is_eg, &
                            &   zw_leerzei, tropos_inside
                 END IF
                 Write(10,*) leerzei_acht,(eg_listfxy(n),n=1,is_eg)
                 nr_faces_write=nr_faces_write+1 
               ELSE 
                 ! Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==0"
                 IF(k==iz0.OR.k==iz1) THEN
                   Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, &
                            &   zw_leerzei, tropos_border 
                 ELSE 
                   Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, &
                            &   zw_leerzei, tropos_inside
                 END IF
                 Write(10,*) leerzei_acht, &
                           & Topo(i,j-1,k)%Topo%exnr, Topo(i,j,k)%Topo%eynr,  &
                           & Topo(i,j,k)%Topo%exnr, Topo(i-1,j,k)%Topo%eynr
                 nr_faces_write=nr_faces_write+1 
               END IF
            END IF
          ELSE
            in_out=Vertices(i-1,j-1,k)%in_out+Vertices(i,j-1,k)%in_out+ &
                   Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
            IF(in_out>=0) THEN
               !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, "(in_out>=0)"
               IF(k==iz0.OR.k==iz1) THEN
                 Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, & 
                         &   zw_leerzei, tropos_border 
               ELSE 
                 Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, & 
                         &   zw_leerzei, tropos_inside 
               END IF
               Write(10,*) leerzei_acht, &
                         & Topo(i,j-1,k)%Topo%exnr, Topo(i,j,k)%Topo%eynr, &
                         & Topo(i,j,k)%Topo%exnr, Topo(i-1,j,k)%Topo%eynr
               nr_faces_write=nr_faces_write+1 
            END IF
          END IF
        END DO
      END DO
    END DO
    !-- Output--Faces_ZX
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
            IF(Vertices(i-1,j,k-1)%in_out==0 .and. &
               Vertices(i  ,j,k-1)%in_out==0 .and. & 
               Vertices(i-1,j,k  )%in_out==0 .and. &
               Vertices(i  ,j,k  )%in_out==0 ) THEN
               !!als  Grenzfläche  'fc' ausgegeben 
               !Write(10,*) "Faces_ZX(i,j,k)%Face als Grenzfläche fc ausgegeben, 4xVert%in_out==0",i,j,k
               !Write(10,*) "                       FCutNr:", Topo(i,j,k)%Topo%FCutNr
            ELSE IF(Faces_ZX(i,j,k)%Face%in_out>-4.and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
               !Reihenfolge:e1-e4, speciell je Face-Def.
               is_eg=0;ec_set=0
               IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
                 ! EdgeListe aufstellen
                 !.e1.....................................
                 IF(Topo(i-1,j,k)%Topo%eznr>0) THEN
                   is_eg=is_eg+1
                   eg_listfzx(is_eg)=Topo(i-1,j,k)%Topo%eznr
                   IF(Faces_ZX(i,j,k)%Face%Edge1%yes_sp==1) THEN
                     IF(Faces_ZX(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                        is_eg=is_eg+1
                        eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                        ec_set=1
                        !Write(10,*) i,j,k,", eg1 - Annulus Test F_ZX"
                     END IF
                   END IF
                 ELSE 
                   IF(Faces_ZX(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                      Faces_ZX(i,j,k)%Face%Edge1%in_out==-1) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                     ec_set=1
                   END IF
                 END IF
                 !.e2.....................................
                 IF(Topo(i,j,k)%Topo%exnr>0) THEN
                    is_eg=is_eg+1
                    eg_listfzx(is_eg)=Topo(i,j,k)%Topo%exnr
                    IF(Faces_ZX(i,j,k)%Face%Edge2%yes_sp==1) THEN
                       IF(Faces_ZX(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                         IF(ec_set==0) THEN
                           is_eg=is_eg+1
                           eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                           ec_set=1
                           !Write(10,*) i,j,k,", eg2 - Annulus Test F_ZX"
                         END IF
                     END IF
                   END IF
                 ELSE 
                    IF(Faces_ZX(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                       Faces_ZX(i,j,k)%Face%Edge2%in_out==-1) THEN
                      IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                        ec_set=1
                      END IF
                    END IF
                 END IF
                 !.e3.....................................
                 IF(Topo(i,j,k)%Topo%eznr>0) THEN
                   is_eg=is_eg+1
                   eg_listfzx(is_eg)=Topo(i,j,k)%Topo%eznr
                   IF(Faces_ZX(i,j,k)%Face%Edge3%yes_sp==1) THEN
                      IF(Faces_ZX(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                          ec_set=1
                          !Write(10,*) i,j,k,", eg3 - Annulus Test F_ZX"
                        END IF
                      END IF
                   END IF
                 ELSE
                   IF(Faces_ZX(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                      Faces_ZX(i,j,k)%Face%Edge3%in_out==-1) THEN
                     IF(ec_set==0) THEN
                       is_eg=is_eg+1
                       eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                       ec_set=1
                     END IF
                   END IF
                 END IF
                 !.e4.....................................
                 IF(Topo(i,j,k-1)%Topo%exnr>0) THEN
                   is_eg=is_eg+1
                   eg_listfzx(is_eg)=Topo(i,j,k-1)%Topo%exnr
                   IF(Faces_ZX(i,j,k)%Face%Edge4%yes_sp==1) THEN
                      IF(Faces_ZX(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                          !Write(10,*) i,j,k,", eg4 - Annulus Test F_ZX"
                        END IF
                      END IF
                   END IF
                 ELSE 
                    IF(Faces_ZX(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                       Faces_ZX(i,j,k)%Face%Edge4%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                       END IF
                    END IF
                 END IF
                 IF(is_eg==2 .and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                 END IF
                 !........................................
                 !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==1"
                 IF(j==iy0.OR.j==iy1) THEN
                   Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, is_eg, &
                           &   zw_leerzei, tropos_border 
                 ELSE 
                   Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, is_eg, &
                           &   zw_leerzei, tropos_inside 
                 END IF
                 Write(10,*) leerzei_acht,(eg_listfzx(n),n=1,is_eg)
                 nr_faces_write=nr_faces_write+1 
               ELSE ! not. Faces_ZX(i,j,k)%Face%ec 
                 !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==0"
                 IF(j==iy0.OR.j==iy1) THEN
                   Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, 4, & 
                           &   zw_leerzei, tropos_border 
                 ELSE 
                   Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, 4, & 
                           &   zw_leerzei, tropos_inside 
                 END IF
                 Write(10,*) leerzei_acht, &
                           & Topo(i-1,j,k)%Topo%eznr,Topo(i,j,k)%Topo%exnr, &
                           & Topo(i,j,k)%Topo%eznr,Topo(i,j,k-1)%Topo%exnr
                 nr_faces_write=nr_faces_write+1 
               END IF !IF/Else Faces_ZX(i,j,k)%Face%ec
            END IF
          ELSE  ! .not. ASSOCIATED(Faces_ZX(i,j,k)%Face
            in_out=Vertices(i-1,j,k-1)%in_out+Vertices(i,j,k-1)%in_out+ &
                   Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
            IF(in_out>=0) THEN
               !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, "(in_out>=0)"
               IF(j==iy0.OR.j==iy1) THEN
                  Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, 4, &
                           &   zw_leerzei, tropos_border 
               ELSE 
                  Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, 4, &
                           &   zw_leerzei, tropos_inside 
               END IF
               Write(10,*) leerzei_acht, &
                         & Topo(i-1,j,k)%Topo%eznr,Topo(i,j,k)%Topo%exnr, &
                         & Topo(i,j,k)%Topo%eznr,Topo(i,j,k-1)%Topo%exnr
               nr_faces_write=nr_faces_write+1 
            END IF
          END IF  ! IF/Else ASSOCIATED(Faces_ZX(i,j,k)%Face
        END DO  ! k
      END DO   ! j
    END DO    ! i
    !-- Output--Faces_YZ
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
            IF(Vertices(i,j-1,k-1)%in_out==0 .and. &
               Vertices(i,j-1,k  )%in_out==0 .and. &
               Vertices(i,j  ,k  )%in_out==0 .and. &
               Vertices(i,j  ,k-1)%in_out==0 ) THEN
               ! als Grenzfläche 'fc' ausgegeben 
               !Write(10,*) "Faces_YZ(i,j,k)%Face als Grenzfläche fc ausgegeben, 4xVert%in_out==0",i,j,k
               !Write(10,*) "                       FCutNr:", Topo(i,j,k)%Topo%FCutNr
            ELSE IF(Faces_YZ(i,j,k)%Face%in_out>-4.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
               !Reihenfolge:e1-e4, speciell je Face-Def.
               is_eg=0;ec_set=0
               IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
                 ! EdgeListe aufstellen
                 !.e1.....................................
                 IF(Topo(i,j,k-1)%Topo%eynr>0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=Topo(i,j,k-1)%Topo%eynr
                   IF(Faces_YZ(i,j,k)%Face%Edge1%yes_sp==1) THEN
                     IF(Faces_YZ(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                       is_eg=is_eg+1
                       eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                       ec_set=1
                       !Write(10,*) i,j,k,",  eg1 - Annulus Test F_YZ"
                     END IF
                   END IF
                 ELSE 
                   IF(Faces_YZ(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                      Faces_YZ(i,j,k)%Face%Edge1%in_out==-1) THEN
                     is_eg=is_eg+1
                     eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                     ec_set=1
                   END IF
                 END IF
                 !.e2.....................................
                 IF(Topo(i,j,k)%Topo%eznr>0) THEN
                    is_eg=is_eg+1
                    eg_listfyz(is_eg)=Topo(i,j,k)%Topo%eznr
                    IF(Faces_YZ(i,j,k)%Face%Edge2%yes_sp==1) THEN
                      IF(Faces_YZ(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                          ec_set=1
                          !Write(10,*) i,j,k,",  eg2 - Annulus Test F_YZ"
                        END IF
                      END IF
                    END IF
                 ELSE 
                    IF(Faces_YZ(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                       Faces_YZ(i,j,k)%Face%Edge2%in_out==-1) THEN
                      IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                        ec_set=1
                      END IF
                    END IF
                 END IF
                 !.e3.....................................
                 IF(Topo(i,j,k)%Topo%eynr>0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=Topo(i,j,k)%Topo%eynr
                   IF(Faces_YZ(i,j,k)%Face%Edge3%yes_sp==1) THEN
                     IF(Faces_YZ(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                         ec_set=1
                         !Write(10,*) i,j,k,",  eg3 - Annulus Test F_YZ"
                       END IF
                     END IF
                   END IF
                 ELSE
                   IF(Faces_YZ(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                      Faces_YZ(i,j,k)%Face%Edge3%in_out==-1) THEN
                     IF(ec_set==0) THEN
                       is_eg=is_eg+1
                       eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                       ec_set=1
                     END IF
                   END IF
                 END IF
                 !.e4.....................................
                 IF(Topo(i,j-1,k)%Topo%eznr>0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=Topo(i,j-1,k)%Topo%eznr
                   IF(Faces_YZ(i,j,k)%Face%Edge4%yes_sp==1) THEN
                     IF(Faces_YZ(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                         !Write(10,*) i,j,k,",  eg4 - Annulus Test F_YZ"
                       END IF
                     END IF
                   END IF
                 ELSE 
                    IF(Faces_YZ(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                       Faces_YZ(i,j,k)%Face%Edge4%in_out==-1) THEN
                      IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                      END IF
                    END IF
                 END IF
                 IF(is_eg==2.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
                      is_eg=is_eg+1
                      eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                 END IF
                 !........................................
                 !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==1"
                 IF(i==ix0.OR.i==ix1) THEN
                   Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, is_eg, &
                           &   zw_leerzei, tropos_border 
                 ELSE 
                   Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, is_eg, &
                           &   zw_leerzei, tropos_inside 
                 END IF
                 Write(10,*) leerzei_acht,(eg_listfyz(n),n=1,is_eg)
                 nr_faces_write=nr_faces_write+1 
               ELSE 
                 !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==0"
                 IF(i==ix0.OR.i==ix1) THEN
                   Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, 4, & 
                            &   zw_leerzei, tropos_border 
                 ELSE 
                   Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, 4, & 
                            &   zw_leerzei, tropos_inside
                 END IF
                 Write(10,*) leerzei_acht, &
                           & Topo(i,j,k-1)%Topo%eynr,Topo(i,j,k)%Topo%eznr, &
                           & Topo(i,j,k)%Topo%eynr, Topo(i,j-1,k)%Topo%eznr
                 nr_faces_write=nr_faces_write+1 
               END IF
            END IF
          ELSE
            in_out=Vertices(i,j-1,k-1)%in_out+Vertices(i,j-1,k)%in_out+ &
                   Vertices(i,j,k)%in_out+Vertices(i,j,k-1)%in_out
            IF(in_out>=0) THEN
               !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, "(in_out>=0)"
               IF(i==ix0.OR.i==ix1) THEN
                 Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, 4, & 
                          &   zw_leerzei, tropos_border 
               ELSE 
                 Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, 4, & 
                          &   zw_leerzei, tropos_inside 
               END IF
               Write(10,*) leerzei_acht, &
                         & Topo(i,j,k-1)%Topo%eynr,Topo(i,j,k)%Topo%eznr, &
                         & Topo(i,j,k)%Topo%eynr, Topo(i,j-1,k)%Topo%eznr
               nr_faces_write=nr_faces_write+1 
            END IF
          END IF
        END DO
      END DO
    END DO
    !-- Output--Cells------------------------------------------------------------
    Write(10,*) "Nr_Cells", stopo_cnr
    nr_cells_write=0
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
             IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
                 & .OR. Cell(i,j,k)%Cell%in_out==8) THEN
                 IF (Cell(i,j,k)%Cell%vc==0) THEN
                   !Reihenfolge F1,F6,F4,F5,F3,F2; (F1,dann gegen Uhrzeiger,F2)
                   ! Write(10,*) "                    ",i,j,k,"-allok-vc=0"
                   ifn=1     !1 'F2'
                   IF(Topo(i,j,k)%Topo%fxynr>0) THEN
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fxynr 
                    ELSE
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                   END IF
                   ifn=ifn+1  !2 'F6'
                   IF(Topo(i,j,k)%Topo%fyznr>0) THEN
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fyznr 
                    ELSE
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                   END IF 
                   ifn=ifn+1  !3 'F4'
                   IF(Topo(i,j,k)%Topo%fzxnr>0) THEN
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fzxnr 
                    ELSE
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                   END IF
                   ifn=ifn+1  !4  'F1'
                   IF(Topo(i,j,k-1)%Topo%fxynr>0) THEN
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k-1)%Topo%fxynr 
                    ELSE
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                   END IF
                   ifn=ifn+1  !5  'F3'
                   IF(Topo(i,j-1,k)%Topo%fzxnr>0) THEN
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j-1,k)%Topo%fzxnr 
                    ELSE
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                   END IF
                   ifn=ifn+1  !6  'F5'
                   IF(Topo(i-1,j,k)%Topo%fyznr>0) THEN
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i-1,j,k)%Topo%fyznr 
                    ELSE
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                   END IF 
                   Write(10,*) "cv",Topo(i,j,k)%Topo%cnr, zw_leerzei, &
                             &     6, zw_leerzei, tropos_border 
                   Write(10,*) leerzei_acht,  &
                              & Topo(i,j,k)%Topo%cfn(4),Topo(i,j,k)%Topo%cfn(2), &
                              & Topo(i,j,k)%Topo%cfn(3),Topo(i,j,k)%Topo%cfn(6), &
                              & Topo(i,j,k)%Topo%cfn(5),Topo(i,j,k)%Topo%cfn(1)
                   nr_cells_write=nr_cells_write+1
                 ELSE  !  IF(Cell(i,j,k)%Cell%vc>0) THEN
                   Topo(i,j,k)%Topo%cfn(1)=Topo(i,j,k)%Topo%FCutNr 
                   ifn=1
                   IF(Topo(i,j,k)%Topo%fxynr>0) THEN
                     ifn=ifn+1
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fxynr 
                   END IF
                   IF(Topo(i,j,k)%Topo%fyznr>0) THEN
                     ifn=ifn+1
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fyznr 
                   END IF 
                   IF(Topo(i,j,k)%Topo%fzxnr>0) THEN
                     ifn=ifn+1
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fzxnr 
                   END IF
                   IF(Topo(i,j,k-1)%Topo%fxynr>0) THEN
                     ifn=ifn+1
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k-1)%Topo%fxynr 
                   END IF
                   IF(Topo(i,j-1,k)%Topo%fzxnr>0) THEN
                     ifn=ifn+1
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j-1,k)%Topo%fzxnr 
                   END IF
                   IF(Topo(i-1,j,k)%Topo%fyznr>0) THEN
                     ifn=ifn+1
                     Topo(i,j,k)%Topo%cfn(ifn)=Topo(i-1,j,k)%Topo%fyznr 
                   END IF 
                   Topo(i,j,k)%Topo%nfn=ifn
                   ! Write(10,*) "                       ",i,j,k," vc > 0 = ", &
                   !       &    Cell(i,j,k)%Cell%vc
                   Write(10,*) "cc",Topo(i,j,k)%Topo%cnr, zw_leerzei, &
                            &  Topo(i,j,k)%Topo%nfn, zw_leerzei, tropos_border
                   Write(10,*) leerzei_acht, &
                              & (Topo(i,j,k)%Topo%cfn(ifn),ifn=1,Topo(i,j,k)%Topo%nfn)
                   nr_cells_write=nr_cells_write+1
                   !CALL WriteCellProt2(Cell(i,j,k)%Cell,i,j,k)
                 END IF  !if(vc==0) then ... else
             END IF
          END IF
        END DO
      END DO
    END DO
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
          IF(.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN !!! .NOT. ASSOCIATED(Cell)
            in_out=Vertices(i-1,j-1,k-1)%in_out &
                  +Vertices(i,j-1,k-1)%in_out &
                  +Vertices(i-1,j,k-1)%in_out &
                  +Vertices(i-1,j-1,k)%in_out &
                  +Vertices(i-1,j,k)%in_out &
                  +Vertices(i,j-1,k)%in_out &
                  +Vertices(i,j,k-1)%in_out &
                  +Vertices(i,j,k)%in_out
            IF (in_out>=0) THEN
             !!IF (Cell(i,j,k)%Cell%in_out>5) THEN
              !Reihenfolge F1,F6,F4,F5,F3,F2; (F1,dann gegen Uhrzeiger,F2)
              !Write(10,*) "                       ",i,j,k," in_out >= 0  if"
              Write(10,*) "ch",Topo(i,j,k)%Topo%cnr, zw_leerzei, 6, & 
                       &  zw_leerzei, tropos_inside
              Write(10,*) leerzei_acht, &
                        & Topo(i,j,k-1)%Topo%fxynr,Topo(i,j,k)%Topo%fyznr, &
                        & Topo(i,j,k)%Topo%fzxnr,Topo(i-1,j,k)%Topo%fyznr, &
                        & Topo(i,j-1,k)%Topo%fzxnr,Topo(i,j,k)%Topo%fxynr
              nr_cells_write=nr_cells_write+1
            ELSE
              !Write(10,*) "                       ",i,j,k," in_out >= 0  else"
              !Topo(i,j,k)%Topo%ctp=-1  ! unterhalb Berg
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO  !ib
  Write(10,*) "Ende Tropo   ",TRIM(FileName), "  ",nb, "Blöcke"
  !------------------ Protokoll zur Ausgabe --------------------------------! 
  Write(OutUnitProt,*) leerzei3,"> Protokoll aus SUBROUTINE WriteTropoCellsParaV()"
  Write(OutUnitProt,*) leerzei3,"  Schreiben in  Datei  *.pva.tropo"
  Write(OutUnitProt,*) leerzei3,"-----------------------------------------------"
  IF (nr_edges_write /= topo_enr ) THEN
    Write(OutUnitProt,*) leerzei5,"analysierte Edges /= gelistete Faces", &
                &         topo_enr,"/=",nr_edges_write
  ELSE 
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Edges : analysierte == gelistete "
  END IF
  IF (nr_faces_write /= topo_fnr ) THEN
    Write(OutUnitProt,*) leerzei5,"analysierte Faces /= gelistete Faces", &
                &         topo_fnr,"/=",nr_faces_write
  ELSE
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Faces : analysierte == gelistete "
  END IF  
  IF (nr_cells_write /= topo_cnr ) THEN 
    Write(OutUnitProt,*) leerzei5,"analysierte Cellen /= gelistete Cellen ", &
                &        topo_cnr," /=",nr_cells_write        
  ELSE
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Cellen: analysierte == gelistete "
  END IF
  Write(OutUnitProt,*) ""
  !------------------ Ende Protokoll zur Ausgabe ----------------------------! 
!!    WRITE(10,'(a8,i8)') cells,nr_cells
!!    DO ib=1,nb
!!      Write(*,*) "                                     Block : ",ib,"\/",nb
!!      CALL Set(Floor(ib))
!!      DO k=iz0+1,iz1
!!        DO j=iy0+1,iy1
!!          DO i=ix0+1,ix1
!!            !CALL WriteCellGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!            IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!              Write(10,*) "...................................."
!!              Write(10,*) "Cell(",i,",",j,",",k,")%Cell% ->"
!!              Write(10,*) "  Face1%mp=",Cell(i,j,k)%Cell%Face1%mp
!!              Write(10,*) "  Face2%mp=",Cell(i,j,k)%Cell%Face2%mp
!!              Write(10,*) "  Face3%mp=",Cell(i,j,k)%Cell%Face3%mp
!!              Write(10,*) "  Face4%mp=",Cell(i,j,k)%Cell%Face4%mp
!!              Write(10,*) "  Face5%mp=",Cell(i,j,k)%Cell%Face5%mp
!!              Write(10,*) "  Face6%mp=",Cell(i,j,k)%Cell%Face6%mp
!!              Write(10,*) "  vc=",Cell(i,j,k)%Cell%vc
!!              WRITE(10,*) (Cell(i,j,k)%Cell%VertCut(ic),ic=1,Cell(i,j,k)%Cell%vc)
!!              Write(10,*) "...................................."
!!            END IF
!!          END DO
!!        END DO
!!      END DO
!!    END DO  !ib

  CLOSE(10)

END SUBROUTINE WriteTropoCellsParaV


SUBROUTINE WriteTropoOroCellsParaV(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib,i,j,k,n,ic,iec,is_eg,in_out
  INTEGER :: ecnr,exnr,eynr,eznr,enr,ec_set
  INTEGER :: fcnr,fxynr,fzxnr,fyznr,fnr
  INTEGER :: cnr,cnr_cut,ifn,nfn
  INTEGER :: eg_listfxy(1:8)=0
  INTEGER :: eg_listfzx(1:8)=0
  INTEGER :: eg_listfyz(1:8)=0
  REAL(8) :: sfac
  CHARACTER(2),PARAMETER :: zw_leerzei='  '
  CHARACTER(8),PARAMETER :: leerzei_acht='        '
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: liste_ec
  !Variablen Analyze/Auswertung
  INTEGER :: nr_cells_write,nr_faces_write,nr_edges_write

  !CALL WriteAuswOutParaVToProt()
  !CALL SetTropoOroCellsParaV(FileName)  ! in WriteTropoCellsParaV

  !WRITE(*,*) 'Output-File im ascii Format wird erzeugt: '
  OPEN(UNIT=10,FILE=TRIM(FileName)//'.pva.gall',STATUS='UNKNOWN')
  CALL Display_OutTroposOroParaV(FileName)
  
  WRITE(10,*) FileName
  !WRITE(10,'(a8,i8)') nodes,nr_out
  WRITE(10,*) "Nr_Nodes",nr_out2
  DO i=1,nr_out
     WRITE(10,*)"n", i, zw_leerzei, xParametricOut(VertOut(i)), &
           &    yParametricOut(VertOut(i)),zParametricOut(VertOut(i))
  END DO
  DO i=nr_out+1,nr_out2
     WRITE(10,*)"n", i, zw_leerzei, xParametricOut(VertOut2(i)), &
           &    yParametricOut(VertOut2(i)),zParametricOut(VertOut2(i))
  END DO

  !WRITE(10,'(a8,i8)') cells,nr_viewcells
  DO ib=1,1 !nb
    !Write(*,*) "                                     Block : ",ib,"\/",nb
    Write(*,*) "                                     z.Zt.: only with  1 Block"  
    CALL Set(Floor(ib))

      !-- Output---EdgeX/EdgeY/EdgeZ/EdgeCut
      Write(10,*)"Nr_Edges",para_grid_enr
      nr_edges_write=0
      DO i=ix0+1,ix1
        DO j=iy0,iy1
          DO k=iz0,iz1
            IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
               IF(MAX(Edges_X(i,j,k)%Edge%Vert1%in_out,Edges_X(i,j,k)%Edge%Vert2%in_out)>0 &
                 .OR. &
                 (Edges_X(i,j,k)%Edge%Vert1%in_out==Edges_X(i,j,k)%Edge%Vert2%in_out.AND. &
                  Edges_X(i,j,k)%Edge%Vert1%in_out==0)) THEN
                 IF(Edges_X(i,j,k)%Edge%yes_sp==1) THEN
                    IF(Edges_X(i,j,k)%Edge%Vert1%nrP==0) THEN
                      Write(10,*) "ex", Topo(i,j,k)%Topo%exnr, zw_leerzei, &
                       & Edges_X(i,j,k)%Edge%VertS%nrP,Edges_X(i,j,k)%Edge%Vert2%nrP
                      nr_edges_write=nr_edges_write+1
                    ELSE IF(Edges_X(i,j,k)%Edge%Vert2%nrP==0) THEN
                      Write(10,*) "ex", Topo(i,j,k)%Topo%exnr, zw_leerzei, &
                        & Edges_X(i,j,k)%Edge%Vert1%nrP,Edges_X(i,j,k)%Edge%VertS%nrP 
                      nr_edges_write=nr_edges_write+1
                    END IF
                 ELSE
                   Write(10,*) "ex", Topo(i,j,k)%Topo%exnr, zw_leerzei, &
                     & Edges_X(i,j,k)%Edge%Vert1%nrP,Edges_X(i,j,k)%Edge%Vert2%nrP 
                   nr_edges_write=nr_edges_write+1
                 END IF
               END IF
            ELSE
              IF(MAX(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)>0 &
                 .OR. &
                 (Vertices(i-1,j,k)%in_out==Vertices(i,j,k)%in_out .AND. &
                  Vertices(i,j,k)%in_out==0)) THEN
                Write(10,*) "ex", Topo(i,j,k)%Topo%exnr, zw_leerzei, &
                  & Vertices(i-1,j,k)%nrP,Vertices(i,j,k)%nrP
                nr_edges_write=nr_edges_write+1
              END IF 
            END IF
          END DO
        END DO
      END DO
      DO i=ix0,ix1
        DO j=iy0+1,iy1
          DO k=iz0,iz1
            IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
              IF(MAX(Edges_Y(i,j,k)%Edge%Vert1%in_out,Edges_Y(i,j,k)%Edge%Vert2%in_out)>0 &
                 .OR. &
                 (Edges_Y(i,j,k)%Edge%Vert1%in_out==Edges_Y(i,j,k)%Edge%Vert2%in_out.AND. &
                  Edges_Y(i,j,k)%Edge%Vert1%in_out==0)) THEN
                 IF(Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
                   IF(Edges_Y(i,j,k)%Edge%Vert1%nrP==0) THEN
                     Write(10,*) "ey", Topo(i,j,k)%Topo%eynr, zw_leerzei, &
                      & Edges_Y(i,j,k)%Edge%VertS%nrP,Edges_Y(i,j,k)%Edge%Vert2%nrP
                     nr_edges_write=nr_edges_write+1
                   ELSE IF(Edges_Y(i,j,k)%Edge%Vert2%nrP==0) THEN
                     Write(10,*) "ey", Topo(i,j,k)%Topo%eynr, zw_leerzei, &
                      & Edges_Y(i,j,k)%Edge%Vert1%nrP,Edges_Y(i,j,k)%Edge%VertS%nrP
                     nr_edges_write=nr_edges_write+1
                   END IF
                 ELSE
                    Write(10,*) "ey", Topo(i,j,k)%Topo%eynr, zw_leerzei, &
                     & Edges_Y(i,j,k)%Edge%Vert1%nrP,Edges_Y(i,j,k)%Edge%Vert2%nrP
                    nr_edges_write=nr_edges_write+1
                 END IF
              END IF
            ELSE
              IF(MAX(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)>0 &
                 .OR. &
                 (Vertices(i,j-1,k)%in_out==Vertices(i,j,k)%in_out .AND. &
                  Vertices(i,j,k)%in_out==0)) THEN
                 Write(10,*) "ey", Topo(i,j,k)%Topo%eynr, zw_leerzei, &
                   & Vertices(i,j-1,k)%nrP,Vertices(i,j,k)%nrP
                 nr_edges_write=nr_edges_write+1
              END IF
            END IF
          END DO
        END DO
      END DO
      DO i=ix0,ix1
        DO j=iy0,iy1
          DO k=iz0+1,iz1
            IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
              IF(MAX(Edges_Z(i,j,k)%Edge%Vert1%in_out,Edges_Z(i,j,k)%Edge%Vert2%in_out)>0 &
                 .OR. &
                 (Edges_Z(i,j,k)%Edge%Vert1%in_out==Edges_Z(i,j,k)%Edge%Vert2%in_out.AND. &
                  Edges_Z(i,j,k)%Edge%Vert1%in_out==0)) THEN
                 IF(Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
                   IF(Edges_Z(i,j,k)%Edge%Vert1%nrP==0) THEN
                     Write(10,*) "ez", Topo(i,j,k)%Topo%eznr, zw_leerzei, &
                        & Edges_Z(i,j,k)%Edge%VertS%nrP,Edges_Z(i,j,k)%Edge%Vert2%nrP
                     nr_edges_write=nr_edges_write+1
                   ELSE IF(Edges_Z(i,j,k)%Edge%Vert2%nrP==0) THEN
                     Write(10,*) "ez", Topo(i,j,k)%Topo%eznr, zw_leerzei, &
                       & Edges_Z(i,j,k)%Edge%Vert1%nrP,Edges_Z(i,j,k)%Edge%VertS%nrP
                     nr_edges_write=nr_edges_write+1
                   END IF
                 ELSE
                   Write(10,*) "ez", Topo(i,j,k)%Topo%eznr, zw_leerzei, &
                     & Edges_Z(i,j,k)%Edge%Vert1%nrP,Edges_Z(i,j,k)%Edge%Vert2%nrP
                   nr_edges_write=nr_edges_write+1
                 END IF
              END IF
            ELSE
              IF(MAX(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)>0 &
                 .OR. &
                 (Vertices(i,j,k-1)%in_out==Vertices(i,j,k)%in_out.AND. &
                  Vertices(i,j,k)%in_out==0)) THEN
                 Write(10,*) "ez", Topo(i,j,k)%Topo%eznr, zw_leerzei, &
                   & Vertices(i,j,k-1)%nrP,Vertices(i,j,k)%nrP
                 nr_edges_write=nr_edges_write+1
              END IF
            END IF
          END DO
        END DO
      END DO
      !EdgeCut-Liste Sortieren- für Ausgabe
      IF(stopo_ecnr>0) THEN
      ALLOCATE(liste_ec(stopo_enr-stopo_ecnr+1:stopo_enr,1:3))
      DO i=ix0+1,ix1
        DO j=iy0+1,iy1
          DO k=iz0,iz1
            IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
              IF(Faces_XY(i,j,k)%Face%egc_nr>0) THEN
                 !Write(10,*) "               aus Faces_XY--------",i,j,k
                 !Write(10,*) "ec",Faces_XY(i,j,k)%Face%egc_nr, zw_leerzei, &
                 !      & Faces_XY(i,j,k)%Face%EdgeCut(1),Faces_XY(i,j,k)%Face%EdgeCut(2)
                 liste_ec(Faces_XY(i,j,k)%Face%egc_nr,1)=Faces_XY(i,j,k)%Face%egc_nr
                 liste_ec(Faces_XY(i,j,k)%Face%egc_nr,2)=Faces_XY(i,j,k)%Face%EdgeCut(1)
                 liste_ec(Faces_XY(i,j,k)%Face%egc_nr,3)=Faces_XY(i,j,k)%Face%EdgeCut(2)    
              END IF
            END IF
          END DO
        END DO
      END DO
      DO i=ix0,ix1
        DO j=iy0+1,iy1
          DO k=iz0+1,iz1
            IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
              IF(Faces_YZ(i,j,k)%Face%egc_nr>0) THEN
                 !Write(10,*) "               aus Faces_YZ--------",i,j,k
                 !Write(10,*) "ec",Faces_YZ(i,j,k)%Face%egc_nr, zw_leerzei, &
                 !      & Faces_YZ(i,j,k)%Face%EdgeCut(1),Faces_YZ(i,j,k)%Face%EdgeCut(2)
                 liste_ec(Faces_YZ(i,j,k)%Face%egc_nr,1)=Faces_YZ(i,j,k)%Face%egc_nr
                 liste_ec(Faces_YZ(i,j,k)%Face%egc_nr,2)=Faces_YZ(i,j,k)%Face%EdgeCut(1)
                 liste_ec(Faces_YZ(i,j,k)%Face%egc_nr,3)=Faces_YZ(i,j,k)%Face%EdgeCut(2)
              END IF
            END IF
          END DO
        END DO
      END DO
      DO i=ix0+1,ix1
        DO j=iy0,iy1
          DO k=iz0+1,iz1
            IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
              IF(Faces_ZX(i,j,k)%Face%egc_nr>0) THEN
                 !Write(10,*) "               aus Faces_ZX--------",i,j,k
                 !Write(10,*) "ec",Faces_ZX(i,j,k)%Face%egc_nr, zw_leerzei, &
                 !      & Faces_ZX(i,j,k)%Face%EdgeCut(1),Faces_ZX(i,j,k)%Face%EdgeCut(2)
                 liste_ec(Faces_ZX(i,j,k)%Face%egc_nr,1)=Faces_ZX(i,j,k)%Face%egc_nr
                 liste_ec(Faces_ZX(i,j,k)%Face%egc_nr,2)=Faces_ZX(i,j,k)%Face%EdgeCut(1)
                 liste_ec(Faces_ZX(i,j,k)%Face%egc_nr,3)=Faces_ZX(i,j,k)%Face%EdgeCut(2)
              END IF
            END IF
          END DO
        END DO
      END DO
      !EdgeCut-Ausgabe
      DO i=stopo_enr-stopo_ecnr+1,stopo_enr 
         Write(10,*) "ec",liste_ec(i,1),zw_leerzei,liste_ec(i,2),liste_ec(i,3)
         nr_edges_write=nr_edges_write+1
      END DO
      DEALLOCATE(liste_ec)
      END IF
 
      !-- Output inside mountain (Orography) Edges
      !---Edge_X----inside mountain---
      DO i=ix0+1,ix1
        DO j=iy0,iy1
          DO k=iz0,iz1
             IF(Topo(i,j,k)%Topo%iexnr>topo_enr) THEN
                !(Edge%in_out==-2 -> Vertices([i-1:i],j,k)%in_out==-1)
                !(Edge%in_out==-1 und [yes_sp od. Vert!%in_out==0])
                IF(Vertices(i-1,j,k)%nrI>0.and.Vertices(i,j,k)%nrI>0) THEN
                  Write(10,*) "ex", Topo(i,j,k)%Topo%iexnr, zw_leerzei, &
                     & Vertices(i-1,j,k)%nrI,Vertices(i,j,k)%nrI
                  nr_edges_write=nr_edges_write+1
                ELSE
                  IF(Edges_X(i,j,k)%Edge%yes_sp==1) THEN
                    IF(Vertices(i-1,j,k)%nrI==-1) THEN
                      Write(10,*) "ex", Topo(i,j,k)%Topo%iexnr, zw_leerzei, &
                         & Edges_X(i,j,k)%Edge%VertS%nrP,Vertices(i,j,k)%nrI
                      nr_edges_write=nr_edges_write+1
                    ELSE
                      Write(10,*) "ex", Topo(i,j,k)%Topo%iexnr, zw_leerzei, &
                        & Vertices(i-1,j,k)%nrI,Edges_X(i,j,k)%Edge%VertS%nrP
                      nr_edges_write=nr_edges_write+1
                    END IF
                  ELSE IF(Edges_X(i,j,k)%Edge%yes_sp==0) THEN
                    IF(Vertices(i-1,j,k)%nrI==-1) THEN
                      Write(10,*) "ex", Topo(i,j,k)%Topo%iexnr, zw_leerzei, &
                         & Vertices(i-1,j,k)%nrP,Vertices(i,j,k)%nrI
                      nr_edges_write=nr_edges_write+1
                    ELSE
                      Write(10,*) "ex", Topo(i,j,k)%Topo%iexnr, zw_leerzei, &
                        & Vertices(i-1,j,k)%nrI,Vertices(i,j,k)%nrP
                      nr_edges_write=nr_edges_write+1
                    END IF
                  END IF
                END IF
             END IF
          END DO
        END DO
      END DO
      !---Edge_Y----inside mountain---
      DO i=ix0,ix1
        DO j=iy0+1,iy1
          DO k=iz0,iz1
             IF(Topo(i,j,k)%Topo%ieynr>topo_enr) THEN
               !Write(10,*) "ey","-----", i,j,k
                !(Edge%in_out==-2  -> Vertices(i,[j-1:j],k)%in_out==-1)
                !(Edge%in_out==-1 und [yes_sp od. Vert!%in_out==0])
                IF(Vertices(i,j-1,k)%nrI>0.and.Vertices(i,j,k)%nrI>0) THEN
                    Write(10,*) "ey", Topo(i,j,k)%Topo%ieynr, zw_leerzei, &
                            & Vertices(i,j-1,k)%nrI,Vertices(i,j,k)%nrI
                    nr_edges_write=nr_edges_write+1
                ELSE
                  IF(Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
                    IF(Vertices(i,j-1,k)%nrI==-1) THEN
                      Write(10,*) "ey", Topo(i,j,k)%Topo%ieynr, zw_leerzei, &
                            & Edges_Y(i,j,k)%Edge%VertS%nrP,Vertices(i,j,k)%nrI
                      nr_edges_write=nr_edges_write+1
                    ELSE
                      Write(10,*) "ey", Topo(i,j,k)%Topo%ieynr, zw_leerzei, &
                            & Vertices(i,j-1,k)%nrI,Edges_Y(i,j,k)%Edge%VertS%nrP
                      nr_edges_write=nr_edges_write+1
                    END IF
                  ELSE IF(Edges_Y(i,j,k)%Edge%yes_sp==0) THEN
                    IF(Vertices(i,j-1,k)%nrI==-1) THEN
                       Write(10,*) "ey", Topo(i,j,k)%Topo%ieynr, zw_leerzei, &
                            & Vertices(i,j-1,k)%nrP,Vertices(i,j,k)%nrI
                       nr_edges_write=nr_edges_write+1
                    ELSE
                       Write(10,*) "ey", Topo(i,j,k)%Topo%ieynr, zw_leerzei, &
                            & Vertices(i,j-1,k)%nrI,Vertices(i,j,k)%nrP
                       nr_edges_write=nr_edges_write+1
                    END IF
                  END IF
                END IF
             END IF
          END DO
        END DO
      END DO
      !---Edge_Z----inside mountain---
      DO i=ix0,ix1
        DO j=iy0,iy1
          DO k=iz0+1,iz1
             IF(Topo(i,j,k)%Topo%ieznr>topo_enr) THEN
               !Write(10,*) "ez","-----", i,j,k
                !(Edge%in_out==-2  -> Vertices(i,j,[k-1:k])%in_out==-1)
                !(Edge%in_out==-1 und [yes_sp od. Vert!%in_out==0])
                IF(Vertices(i,j,k-1)%nrI>0.and.Vertices(i,j,k)%nrI>0) THEN
                   Write(10,*) "ez", Topo(i,j,k)%Topo%ieznr, zw_leerzei, &
                            & Vertices(i,j,k-1)%nrI,Vertices(i,j,k)%nrI
                   nr_edges_write=nr_edges_write+1
                ELSE
                  IF(Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
                    IF(Vertices(i,j,k-1)%nrI==-1)THEN
                       Write(10,*) "ez", Topo(i,j,k)%Topo%ieznr, zw_leerzei, &
                          & Edges_Z(i,j,k)%Edge%VertS%nrP,Vertices(i,j,k)%nrI
                       nr_edges_write=nr_edges_write+1
                    ELSE
                       Write(10,*) "ez", Topo(i,j,k)%Topo%ieznr, zw_leerzei, &
                          & Vertices(i,j,k-1)%nrI,Edges_Z(i,j,k)%Edge%VertS%nrP
                       nr_edges_write=nr_edges_write+1
                    END IF
                  ELSE IF(Edges_Z(i,j,k)%Edge%yes_sp==0) THEN
                    IF(Vertices(i,j,k-1)%nrI==-1)THEN
                       Write(10,*) "ez", Topo(i,j,k)%Topo%ieznr, zw_leerzei, &
                               & Vertices(i,j,k-1)%nrP,Vertices(i,j,k)%nrI
                       nr_edges_write=nr_edges_write+1
                    ELSE
                       Write(10,*) "ez", Topo(i,j,k)%Topo%ieznr, zw_leerzei, &
                               & Vertices(i,j,k-1)%nrI,Vertices(i,j,k)%nrP
                       nr_edges_write=nr_edges_write+1
                    END IF
                  END IF
                END IF
             END IF
          END DO
        END DO
      END DO
      ! end inside  mountain Edges
      !-- Output--FaceCut,Faces_XY,Faces_XZ,Faces_YZ ---------------------------------------
      !Write(10,*) "Nr_Faces", topo_fnr 
      Write(10,*) "Nr_Faces", para_grid_fnr
      nr_faces_write=0
      DO i=ix0+1,ix1
        DO j=iy0+1,iy1
          DO k=iz0+1,iz1
            IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
              IF(Cell(i,j,k)%Cell%vc>0) THEN
                 Write(10,*) "fc", Topo(i,j,k)%Topo%FCutNr, zw_leerzei,  &
                         &   Cell(i,j,k)%Cell%vc, zw_leerzei, tropos_border
                 Write(10,*) leerzei_acht, &
                            & (Topo(i,j,k)%Topo%fceg_nrcut(iec),iec=1,Cell(i,j,k)%Cell%vc)
                 nr_faces_write=nr_faces_write+1
              END IF
            END IF
          END DO
        END DO
      END DO
      DO i=ix0+1,ix1
        DO j=iy0+1,iy1
          DO k=iz0,iz1
            IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
              IF(Vertices(i-1,j-1,k)%in_out==0 .and. &
                 Vertices(i  ,j-1,k)%in_out==0 .and. &
                 Vertices(i-1,j  ,k)%in_out==0 .and. &
                 Vertices(i  ,j  ,k)%in_out==0) THEN
                 IF(k==iz0) THEN  ! Überarbeiten mehrere Blöcke 'igz0'
                   ! wenn Boden-Modell
                     Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, &
                              &   zw_leerzei, tropos_border 
                     Write(10,*) leerzei_acht, &
                             & Topo(i,j-1,k)%Topo%exnr, Topo(i,j,k)%Topo%eynr,  &
                             & Topo(i,j,k)%Topo%exnr, Topo(i-1,j,k)%Topo%eynr
                     nr_faces_write=nr_faces_write+1
                 ELSE
                   ! als Grenzfläche 'fc' ausgegeben 
                 END IF
              ELSE IF(Faces_XY(i,j,k)%Face%in_out>-4.and.Faces_XY(i,j,k)%Face%numberVert>2) THEN
                 !Reihenfolge:e1-e4, speciell je Face-Def.
                 is_eg=0;ec_set=0
                 IF(Faces_XY(i,j,k)%Face%ec==1) THEN
                    ! EdgeListe aufstellen
                   !.e1.....................................
                   IF(Topo(i,j-1,k)%Topo%exnr>0) THEN
                     is_eg=is_eg+1
                     eg_listfxy(is_eg)=Topo(i,j-1,k)%Topo%exnr
                     IF(Faces_XY(i,j,k)%Face%Edge1%yes_sp==1) THEN
                        IF(Faces_XY(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                           is_eg=is_eg+1
                           eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                           ec_set=1
                           !Write(10,*) i,j,k,", eg1 - Annulus Test F_XY"
                        END IF
                     END IF
                   ELSE 
                     IF(Faces_XY(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                        Faces_XY(i,j,k)%Face%Edge1%in_out==-1) THEN
                         is_eg=is_eg+1
                         eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                         ec_set=1
                     END IF
                   END IF
                   !.e2.....................................
                   IF(Topo(i,j,k)%Topo%eynr>0) THEN
                      is_eg=is_eg+1
                      eg_listfxy(is_eg)=Topo(i,j,k)%Topo%eynr
                      IF(Faces_XY(i,j,k)%Face%Edge2%yes_sp==1) THEN
                        IF(Faces_XY(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                          IF(ec_set==0) THEN
                            is_eg=is_eg+1
                            eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                            ec_set=1
                            !Write(10,*) i,j,k,",  eg2 - Annulus Test F_XY"
                          END IF
                        END IF
                      END IF
                   ELSE 
                      IF(Faces_XY(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                         Faces_XY(i,j,k)%Face%Edge2%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                          ec_set=1
                        END IF
                      END IF
                   END IF
                   !.e3.....................................
                   IF(Topo(i,j,k)%Topo%exnr>0) THEN
                     is_eg=is_eg+1
                     eg_listfxy(is_eg)=Topo(i,j,k)%Topo%exnr
                     IF(Faces_XY(i,j,k)%Face%Edge3%yes_sp==1) THEN
                       IF(Faces_XY(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                         IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                          ec_set=1
                          !Write(10,*) i,j,k,",  eg3 - Annulus Test F_XY"
                         END IF
                       END IF
                     END IF
                   ELSE
                     IF(Faces_XY(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                        Faces_XY(i,j,k)%Face%Edge3%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                         ec_set=1
                       END IF
                     END IF
                   END IF
                   !.e4.....................................
                   IF(Topo(i-1,j,k)%Topo%eynr>0) THEN
                     is_eg=is_eg+1
                     eg_listfxy(is_eg)=Topo(i-1,j,k)%Topo%eynr
                     IF(Faces_XY(i,j,k)%Face%Edge4%yes_sp==1) THEN
                       IF(Faces_XY(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                         IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                          ! Write(10,*) i,j,k,",  eg4 - Annulus Test F_XY"
                         END IF
                       END IF
                     END IF 
                   ELSE 
                      IF(Faces_XY(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                         Faces_XY(i,j,k)%Face%Edge4%in_out==-1) THEN
                        IF(ec_set==0) THEN
                           is_eg=is_eg+1
                           eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                        END IF
                      END IF
                   END IF
                   IF(is_eg==2 .and. Faces_XY(i,j,k)%Face%NumberVert>2) THEN
                      is_eg=is_eg+1
                      eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                   END IF
                   !........................................
                   !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==1"
                   IF(k==iz0.OR.k==iz1) THEN
                      Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, is_eg, &
                              &   zw_leerzei, tropos_border 
                   ELSE 
                      Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, is_eg, &
                              &   zw_leerzei, tropos_inside
                   END IF
                   Write(10,*) leerzei_acht,(eg_listfxy(n),n=1,is_eg)
                   nr_faces_write=nr_faces_write+1
                 ELSE 
                   !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==0"
                   IF(k==iz0.OR.k==iz1) THEN
                     Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, &
                              &   zw_leerzei, tropos_border 
                   ELSE 
                     Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, &
                              &   zw_leerzei, tropos_inside
                   END IF
                   Write(10,*) leerzei_acht, &
                             & Topo(i,j-1,k)%Topo%exnr, Topo(i,j,k)%Topo%eynr,  &
                             & Topo(i,j,k)%Topo%exnr, Topo(i-1,j,k)%Topo%eynr
                   nr_faces_write=nr_faces_write+1
                 END IF
              END IF
            ELSE
              in_out=Vertices(i-1,j-1,k)%in_out+Vertices(i,j-1,k)%in_out+ &
                     Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
              IF(in_out>=0) THEN
                 !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, "(in_out>=0)"
                 IF(k==iz0.OR.k==iz1) THEN
                   Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, & 
                           &   zw_leerzei, tropos_border 
                 ELSE 
                   Write(10,*) "fxy", Topo(i,j,k)%Topo%fxynr, zw_leerzei, 4, & 
                           &   zw_leerzei, tropos_inside 
                 END IF
                 Write(10,*) leerzei_acht, &
                           & Topo(i,j-1,k)%Topo%exnr, Topo(i,j,k)%Topo%eynr, &
                           & Topo(i,j,k)%Topo%exnr, Topo(i-1,j,k)%Topo%eynr
                 nr_faces_write=nr_faces_write+1
              END IF
            END IF
          END DO
        END DO
      END DO
      DO i=ix0+1,ix1
        DO j=iy0,iy1
          DO k=iz0+1,iz1
            IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
              IF(Vertices(i-1,j,k-1)%in_out==0 .and. &
                 Vertices(i  ,j,k-1)%in_out==0 .and. & 
                 Vertices(i-1,j,k  )%in_out==0 .and. &
                 Vertices(i  ,j,k  )%in_out==0 ) THEN
                 ! als Grenzfläche 'fc' ausgegeben 
              ELSE IF(Faces_ZX(i,j,k)%Face%in_out>-4.and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
                 !Reihenfolge:e1-e4, speciell je Face-Def.
                 is_eg=0;ec_set=0
                 IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
                   ! EdgeListe aufstellen
                   !.e1.....................................
                   IF(Topo(i-1,j,k)%Topo%eznr>0) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=Topo(i-1,j,k)%Topo%eznr
                     IF(Faces_ZX(i,j,k)%Face%Edge1%yes_sp==1) THEN
                       IF(Faces_ZX(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                          is_eg=is_eg+1
                          eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                          ec_set=1
                          !Write(10,*) i,j,k,", eg1 - Annulus Test F_ZX"
                       END IF
                     END IF
                   ELSE 
                     IF(Faces_ZX(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                        Faces_ZX(i,j,k)%Face%Edge1%in_out==-1) THEN
                       is_eg=is_eg+1
                       eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                       ec_set=1
                     END IF
                   END IF
                   !.e2.....................................
                   IF(Topo(i,j,k)%Topo%exnr>0) THEN
                      is_eg=is_eg+1
                      eg_listfzx(is_eg)=Topo(i,j,k)%Topo%exnr
                      IF(Faces_ZX(i,j,k)%Face%Edge2%yes_sp==1) THEN
                         IF(Faces_ZX(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                           IF(ec_set==0) THEN
                             is_eg=is_eg+1
                             eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                             ec_set=1
                             !Write(10,*) i,j,k,", eg2 - Annulus Test F_ZX"
                           END IF
                       END IF
                     END IF
                   ELSE 
                      IF(Faces_ZX(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                         Faces_ZX(i,j,k)%Face%Edge2%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                          ec_set=1
                        END IF
                      END IF
                   END IF
                   !.e3.....................................
                   IF(Topo(i,j,k)%Topo%eznr>0) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=Topo(i,j,k)%Topo%eznr
                     IF(Faces_ZX(i,j,k)%Face%Edge3%yes_sp==1) THEN
                        IF(Faces_ZX(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                          IF(ec_set==0) THEN
                            is_eg=is_eg+1
                            eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                            ec_set=1
                            !Write(10,*) i,j,k,", eg3 - Annulus Test F_ZX"
                          END IF
                        END IF
                     END IF
                   ELSE
                     IF(Faces_ZX(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                        Faces_ZX(i,j,k)%Face%Edge3%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                         ec_set=1
                       END IF
                     END IF
                   END IF
                   !.e4.....................................
                   IF(Topo(i,j,k-1)%Topo%exnr>0) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=Topo(i,j,k-1)%Topo%exnr
                     IF(Faces_ZX(i,j,k)%Face%Edge4%yes_sp==1) THEN
                        IF(Faces_ZX(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                          IF(ec_set==0) THEN
                            is_eg=is_eg+1
                            eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                            !Write(10,*) i,j,k,", eg4 - Annulus Test F_ZX"
                          END IF
                        END IF
                     END IF
                   ELSE 
                      IF(Faces_ZX(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                         Faces_ZX(i,j,k)%Face%Edge4%in_out==-1) THEN
                         IF(ec_set==0) THEN
                           is_eg=is_eg+1
                           eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                         END IF
                      END IF
                   END IF
                   IF(is_eg==2 .and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
                       is_eg=is_eg+1
                       eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                   END IF
                   !........................................
                   !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==1"
                   IF(j==iy0.OR.j==iy1) THEN
                     Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, is_eg, &
                             &   zw_leerzei, tropos_border 
                   ELSE 
                     Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, is_eg, &
                             &   zw_leerzei, tropos_inside 
                   END IF
                   Write(10,*) leerzei_acht,(eg_listfzx(n),n=1,is_eg)
                   nr_faces_write=nr_faces_write+1
                 ELSE ! not. Faces_ZX(i,j,k)%Face%ec 
                   !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==0"
                   IF(j==iy0.OR.j==iy1) THEN
                     Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, 4, & 
                             &   zw_leerzei, tropos_border 
                   ELSE 
                     Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, 4, & 
                             &   zw_leerzei, tropos_inside 
                   END IF
                   Write(10,*) leerzei_acht, &
                             & Topo(i-1,j,k)%Topo%eznr,Topo(i,j,k)%Topo%exnr, &
                             & Topo(i,j,k)%Topo%eznr,Topo(i,j,k-1)%Topo%exnr
                   nr_faces_write=nr_faces_write+1
                 END IF !IF/Else Faces_ZX(i,j,k)%Face%ec
              END IF
            ELSE  ! .not. ASSOCIATED(Faces_ZX(i,j,k)%Face
              in_out=Vertices(i-1,j,k-1)%in_out+Vertices(i,j,k-1)%in_out+ &
                     Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
              IF(in_out>=0) THEN
                 ! Write(10,*) "                  i =", i, " j = ",j, " k = ",k, "(in_out>=0)"
                 IF(j==iy0.OR.j==iy1) THEN
                    Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, 4, &
                             &   zw_leerzei, tropos_border 
                 ELSE 
                    Write(10,*) "fzx", Topo(i,j,k)%Topo%fzxnr, zw_leerzei, 4, &
                             &   zw_leerzei, tropos_inside 
                 END IF
                 Write(10,*) leerzei_acht, &
                           & Topo(i-1,j,k)%Topo%eznr,Topo(i,j,k)%Topo%exnr, &
                           & Topo(i,j,k)%Topo%eznr,Topo(i,j,k-1)%Topo%exnr
                 nr_faces_write=nr_faces_write+1
              END IF
            END IF  ! IF/Else ASSOCIATED(Faces_ZX(i,j,k)%Face
          END DO  ! k
        END DO   ! j
      END DO    ! i
      DO i=ix0,ix1
        DO j=iy0+1,iy1
          DO k=iz0+1,iz1
            IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
              IF(Vertices(i,j-1,k-1)%in_out==0 .and. &
                 Vertices(i,j-1,k  )%in_out==0 .and. &
                 Vertices(i,j  ,k  )%in_out==0 .and. &
                 Vertices(i,j  ,k-1)%in_out==0 ) THEN
                 ! als Grenzfläche 'fc' ausgegeben 
              ELSE IF(Faces_YZ(i,j,k)%Face%in_out>-4.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
                 !Reihenfolge:e1-e4, speciell je Face-Def.
                 is_eg=0;ec_set=0
                 IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
                   ! EdgeListe aufstellen
                   !.e1.....................................
                   IF(Topo(i,j,k-1)%Topo%eynr>0) THEN
                     is_eg=is_eg+1
                     eg_listfyz(is_eg)=Topo(i,j,k-1)%Topo%eynr
                     IF(Faces_YZ(i,j,k)%Face%Edge1%yes_sp==1) THEN
                       IF(Faces_YZ(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                         is_eg=is_eg+1
                         eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                         ec_set=1
                         !Write(10,*) i,j,k,",  eg1 - Annulus Test F_YZ"
                       END IF
                     END IF
                   ELSE 
                     IF(Faces_YZ(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                        Faces_YZ(i,j,k)%Face%Edge1%in_out==-1) THEN
                       is_eg=is_eg+1
                       eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                       ec_set=1
                     END IF
                   END IF
                   !.e2.....................................
                   IF(Topo(i,j,k)%Topo%eznr>0) THEN
                      is_eg=is_eg+1
                      eg_listfyz(is_eg)=Topo(i,j,k)%Topo%eznr
                      IF(Faces_YZ(i,j,k)%Face%Edge2%yes_sp==1) THEN
                        IF(Faces_YZ(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                          IF(ec_set==0) THEN
                            is_eg=is_eg+1
                            eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                            ec_set=1
                            !Write(10,*) i,j,k,",  eg2 - Annulus Test F_YZ"
                          END IF
                        END IF
                      END IF
                   ELSE 
                      IF(Faces_YZ(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                         Faces_YZ(i,j,k)%Face%Edge2%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                          ec_set=1
                        END IF
                      END IF
                   END IF
                   !.e3.....................................
                   IF(Topo(i,j,k)%Topo%eynr>0) THEN
                     is_eg=is_eg+1
                     eg_listfyz(is_eg)=Topo(i,j,k)%Topo%eynr
                     IF(Faces_YZ(i,j,k)%Face%Edge3%yes_sp==1) THEN
                       IF(Faces_YZ(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                         IF(ec_set==0) THEN
                           is_eg=is_eg+1
                           eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                           ec_set=1
                           !Write(10,*) i,j,k,",  eg3 - Annulus Test F_YZ"
                         END IF
                       END IF
                     END IF
                   ELSE
                     IF(Faces_YZ(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                        Faces_YZ(i,j,k)%Face%Edge3%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                         ec_set=1
                       END IF
                     END IF
                   END IF
                   !.e4.....................................
                   IF(Topo(i,j-1,k)%Topo%eznr>0) THEN
                     is_eg=is_eg+1
                     eg_listfyz(is_eg)=Topo(i,j-1,k)%Topo%eznr
                     IF(Faces_YZ(i,j,k)%Face%Edge4%yes_sp==1) THEN
                       IF(Faces_YZ(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                         IF(ec_set==0) THEN
                           is_eg=is_eg+1
                           eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                           !Write(10,*) i,j,k,",  eg4 - Annulus Test F_YZ"
                         END IF
                       END IF
                     END IF
                   ELSE 
                      IF(Faces_YZ(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                         Faces_YZ(i,j,k)%Face%Edge4%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                        END IF
                      END IF
                   END IF
                   IF(is_eg==2.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
                        is_eg=is_eg+1
                        eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                   END IF
                   !........................................
                   !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==1"
                   IF(i==ix0.OR.i==ix1) THEN
                     Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, is_eg, &
                             &   zw_leerzei, tropos_border 
                   ELSE 
                     Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, is_eg, &
                             &   zw_leerzei, tropos_inside 
                   END IF
                   Write(10,*) leerzei_acht,(eg_listfyz(n),n=1,is_eg)
                   nr_faces_write=nr_faces_write+1
                 ELSE 
                   !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " ec==0"
                   IF(i==ix0.OR.i==ix1) THEN
                     Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, 4, & 
                              &   zw_leerzei, tropos_border 
                   ELSE 
                     Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, 4, & 
                              &   zw_leerzei, tropos_inside
                   END IF
                   Write(10,*) leerzei_acht, &
                             & Topo(i,j,k-1)%Topo%eynr,Topo(i,j,k)%Topo%eznr, &
                             & Topo(i,j,k)%Topo%eynr, Topo(i,j-1,k)%Topo%eznr
                   nr_faces_write=nr_faces_write+1
                 END IF
              END IF
            ELSE
              in_out=Vertices(i,j-1,k-1)%in_out+Vertices(i,j-1,k)%in_out+ &
                     Vertices(i,j,k)%in_out+Vertices(i,j,k-1)%in_out
              IF(in_out>=0) THEN
                 !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, "(in_out>=0)"
                 IF(i==ix0.OR.i==ix1) THEN
                   Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, 4, & 
                            &   zw_leerzei, tropos_border 
                 ELSE 
                   Write(10,*) "fyz", Topo(i,j,k)%Topo%fyznr, zw_leerzei, 4, & 
                            &   zw_leerzei, tropos_inside 
                 END IF
                 Write(10,*) leerzei_acht, &
                           & Topo(i,j,k-1)%Topo%eynr,Topo(i,j,k)%Topo%eznr, &
                           & Topo(i,j,k)%Topo%eynr, Topo(i,j-1,k)%Topo%eznr
                 nr_faces_write=nr_faces_write+1
              END IF
            END IF
          END DO
        END DO
      END DO
      !-- Output all inside mountain Faces
      !---Faces_XY----inside mountain-----
      !write(*,*)"---Faces_XY----inside mountain-----"
      DO i=ix0+1,ix1
        DO j=iy0+1,iy1
          DO k=iz0,iz1
            is_eg=0;ec_set=0
            IF(Topo(i,j,k)%Topo%ifxynr>topo_fnr) THEN
               !Write(10,*) "F_XY                  i =", i, " j = ",j, " k = ",k
               !.e1.................................................
               is_eg=is_eg+1
               eg_listfxy(is_eg)=Topo(i,j-1,k)%Topo%iexnr
               IF(eg_listfxy(is_eg)<0) THEN 
                  eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                  ec_set=1
               END IF
               !Position entsprechend 'Edges_X(i,j-1,k)%Edge%Vert2'
               IF(Vertices(i,j-1,k)%in_out==1.AND.ec_set==0) THEN
                  is_eg=is_eg+1
                  eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                  ec_set=1
               END IF
               !.e2................................................
               IF(Topo(i,j,k)%Topo%ieynr>0) THEN
                  is_eg=is_eg+1
                  eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ieynr
               ELSE
                  If(ec_set==0) THEN
                    is_eg=is_eg+1
                    eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                    ec_set=1
                  END IF
               END IF
               !Position entsprechend 'Edges_Y(i,j,k)%Edge%Vert2'
               IF(Vertices(i,j,k)%in_out==1.and.ec_set==0) THEN
                  is_eg=is_eg+1
                  eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                  ec_set=1
               END IF
               !.e3...............................................
               IF(Topo(i,j,k)%Topo%iexnr>0) THEN
                  is_eg=is_eg+1
                  eg_listfxy(is_eg)=Topo(i,j,k)%Topo%iexnr
               ELSE
                  IF(ec_set==0) THEN 
                    is_eg=is_eg+1
                    eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                    ec_set=1
                  END IF
               END IF
               !Position entsprechend 'Edges_X(i,j,k)%Edge%Vert1'
               IF(Vertices(i-1,j,k)%in_out==1.AND.ec_set==0) THEN
                  is_eg=is_eg+1
                  eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                  ec_set=1
               END IF
               !.e4...............................................
               IF(Topo(i-1,j,k)%Topo%ieynr>0) THEN
                 is_eg=is_eg+1
                 eg_listfxy(is_eg)=Topo(i-1,j,k)%Topo%ieynr 
               ELSE 
                 IF(ec_set==0) THEN 
                    is_eg=is_eg+1
                    eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
                    ec_set=1
                 END IF
               END IF
               !Position ensprechend 'Edges_Y(i-1,j,k)%Edge%Vert1'
               IF(Vertices(i-1,j-1,k)%in_out==1.and.ec_set==0) THEN
                  is_eg=is_eg+1
                  eg_listfxy(is_eg)=Topo(i,j,k)%Topo%ecnrxy
               END IF

             !IF(k==iz0.OR.k==iz1) THEN
             !ELSE 
               !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " inside"
               Write(10,*) "fxy", Topo(i,j,k)%Topo%ifxynr, zw_leerzei, is_eg, & 
                        &   zw_leerzei, soil_inside
             !END IF
               !Write(10,*) leerzei_acht, &
               !          & Topo(i,j-1,k)%Topo%iexnr, Topo(i,j,k)%Topo%ieynr, &
               !          & Topo(i,j,k)%Topo%iexnr, Topo(i-1,j,k)%Topo%ieynr
               Write(10,*) leerzei_acht,(eg_listfxy(n),n=1,is_eg)
               nr_faces_write=nr_faces_write+1
            END IF
          END DO
        END DO
      END DO
      !---Faces_ZX----inside mountain-----
      !write(*,*)"---Faces_ZX----inside mountain-----"
      DO i=ix0+1,ix1
        DO j=iy0,iy1
          DO k=iz0+1,iz1
             is_eg=0;ec_set=0
             IF (Topo(i,j,k)%Topo%ifzxnr>topo_fnr) THEN
               ! in Bearbeitung
               !IF(i == 25 .and. j ==  50 .and. k ==  1 ) THEN
               !   Write(10,*) "Analyze:i = 25  j =  50  k =  1 "
               !   Write(*,*) "Analyze:i = 25  j =  50  k =  1 "
               !END IF
                 !Write(10,*) "F_ZX                  i =", i, " j = ",j, " k = ",k
                 !Belegen eg_list + Auswertung inside or. border
                 !.e1...............................................
                 is_eg=is_eg+1
                 eg_listfzx(is_eg)=Topo(i-1,j,k)%Topo%ieznr
                 IF(eg_listfzx(is_eg)<0) THEN 
                    eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                    ec_set=1
                 END IF
                 !Position entsprechend 'Edges_Z(i-1,j,k)%Edge%Vert2'
                 IF(Vertices(i-1,j,k)%in_out==1.and.ec_set==0) THEN
                    is_eg=is_eg+1
                    eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                    ec_set=1  
                 END IF
                 !.e2................................................
                 IF(Topo(i,j,k)%Topo%iexnr>0) THEN
                   is_eg=is_eg+1
                   eg_listfzx(is_eg)=Topo(i,j,k)%Topo%iexnr
                 ELSE
                   IF(ec_set==0) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                     ec_set=1
                   END IF
                 END IF
                 !Position entsprechend 'Edges_X(i,j,k)%Edge%Vert2'
                 IF(Vertices(i,j,k)%in_out==1.and.ec_set==0) THEN
                   is_eg=is_eg+1
                   eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                   ec_set=1 
                 END IF
                 !.e3..............................................
                 IF(Topo(i,j,k)%Topo%ieznr>0) THEN
                    is_eg=is_eg+1
                    eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ieznr
                 ELSE
                   IF(ec_set==0) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                     ec_set=1
                   END IF
                 END IF      
                 !Position entsprechend 'Edges_Z(i,j,k)%Edge%Vert1'
                 IF(Vertices(i,j,k-1)%in_out==1.and.ec_set==0) THEN 
                    is_eg=is_eg+1
                    eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                    ec_set=1
                 END IF
                 !.e4...............................................
                 IF(Topo(i,j,k-1)%Topo%iexnr>0) THEN
                    is_eg=is_eg+1
                    eg_listfzx(is_eg)=Topo(i,j,k-1)%Topo%iexnr
                 ELSE
                    IF(ec_set==0) THEN
                      is_eg=is_eg+1
                      eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                      ec_set=1
                    END IF
                 END IF 
                 !Position entsprechend 'Edges_X(i,j,k-1)%Edge%Vert1'
                 IF(Vertices(i-1,j,k-1)%in_out==1.and.ec_set==0) THEN
                    is_eg=is_eg+1
                    eg_listfzx(is_eg)=Topo(i,j,k)%Topo%ecnrzx
                 END IF
                !IF(j==iy0.OR.j==iy1) THEN
                !ELSE 
                !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " inside"
                 Write(10,*) "fzx", Topo(i,j,k)%Topo%ifzxnr, zw_leerzei, is_eg, &
                            &   zw_leerzei, soil_inside 
                !END IF
                   !Write(10,*) leerzei_acht, &
                   !      & Topo(i-1,j,k)%Topo%ieznr,Topo(i,j,k)%Topo%iexnr, &
                   !      & Topo(i,j,k)%Topo%ieznr,Topo(i,j,k-1)%Topo%iexnr
                 Write(10,*) leerzei_acht,(eg_listfzx(n),n=1,is_eg)
                 nr_faces_write=nr_faces_write+1
             END IF
          END DO
        END DO
      END DO
      !---Faces_YZ----inside mountain-----
      !write(*,*)"---Faces_YZ----inside mountain-----"
      DO i=ix0,ix1
        DO j=iy0+1,iy1
          DO k=iz0+1,iz1
            is_eg=0;ec_set=0
            IF(Topo(i,j,k)%Topo%ifyznr>topo_fnr) THEN
               !Write(10,*) "F_YZ                  i =", i, " j = ",j, " k = ",k
               !.e1................................................
               is_eg=is_eg+1
               eg_listfyz(is_eg)=Topo(i,j,k-1)%Topo%ieynr
               IF(eg_listfyz(is_eg)<0) THEN
                   eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                   ec_set=1
               END IF
               !Position entsprechend 'Edges_Y(i,j,k-1)%Edge%Vert2'
               IF(Vertices(i,j,k-1)%in_out==1.AND.ec_set==0) THEN
                  is_eg=is_eg+1
                  eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                  ec_set=1
               END IF
               !.e2................................................
               IF(Topo(i,j,k)%Topo%ieznr>0) THEN
                  is_eg=is_eg+1
                  eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ieznr
               ELSE
                  IF(ec_set==0) THEN
                    is_eg=is_eg+1
                    eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                    ec_set=1
                  END IF
               END IF
               !Position entsprechend 'Edges_Z(i,j,k)%Edge%Vert2'
               IF(Vertices(i,j,k)%in_out==1.AND.ec_set==0) THEN 
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                   ec_set=1
               END IF
               !.e3...............................................
               IF(Topo(i,j,k)%Topo%ieynr>0) THEN
                 is_eg=is_eg+1
                 eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ieynr
               ELSE
                 IF(ec_set==0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                   ec_set=1
                 END IF
               END IF
               !Position entsprechend 'Edges_Y(i,j,k)%Edge%Vert1'
               IF(Vertices(i,j-1,k)%in_out==1.AND.ec_set==0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                   ec_set=1
               END IF
               !.e4...............................................
               IF(Topo(i,j-1,k)%Topo%ieznr>0) THEN 
                  is_eg=is_eg+1
                  eg_listfyz(is_eg)=Topo(i,j-1,k)%Topo%ieznr
               ELSE
                 IF(ec_set==0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                   ec_set=1
                 END IF
               END IF
               !Position entsprechend 'Edges_Z(i,j-1,k)%Edge%Vert1'
               IF(Vertices(i,j-1,k-1)%in_out==1.AND.ec_set==0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
               END IF
               !IF(i==ix0.OR.i==ix1) THEN
               !ELSE 
               !Write(10,*) "                  i =", i, " j = ",j, " k = ",k, " inside"
               Write(10,*) "fyz", Topo(i,j,k)%Topo%ifyznr, zw_leerzei, is_eg, & 
                         &   zw_leerzei, soil_inside 
               !END IF
               !Write(10,*) leerzei_acht, &
               !        & Topo(i,j,k-1)%Topo%ieynr,Topo(i,j,k)%Topo%ieznr, &
               !        & Topo(i,j,k)%Topo%ieynr, Topo(i,j-1,k)%Topo%ieznr
               Write(10,*) leerzei_acht,(eg_listfyz(n),n=1,is_eg)
               nr_faces_write=nr_faces_write+1
            END IF
          END DO
        END DO
      END DO
      !para_grid_fnr
      !para_oro_fnr 
      ! end inside  mountain faces 
      !-- Output--Cells
      Write(10,*) "Nr_Cells", para_grid_cnr !topo_cnr 
      nr_cells_write=0
      DO i=ix0+1,ix1
        DO j=iy0+1,iy1
          DO k=iz0+1,iz1
            IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
               IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
                   & .OR. Cell(i,j,k)%Cell%in_out==8) THEN
                   IF (Cell(i,j,k)%Cell%vc==0) THEN
                     !Special mit Face%in_out=0:  mit 4x in_out==0 je Vert1-4
                     !-------------------------   
                     !Reihenfolge F1,F6,F4,F5,F3,F2; (F1,dann gegen Uhrzeiger,F2)
                     ifn=1     !1 'F2'
                     IF(Topo(i,j,k)%Topo%fxynr>0) THEN
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fxynr 
                      ELSE
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     END IF
                     ifn=ifn+1  !2 'F6'
                     IF(Topo(i,j,k)%Topo%fyznr>0) THEN
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fyznr 
                      ELSE
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     END IF 
                     ifn=ifn+1  !3 'F4'
                     IF(Topo(i,j,k)%Topo%fzxnr>0) THEN
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fzxnr 
                      ELSE
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     END IF
                     ifn=ifn+1  !4  'F1'
                     IF(Topo(i,j,k-1)%Topo%fxynr>0) THEN
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k-1)%Topo%fxynr 
                      ELSE
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     END IF
                     ifn=ifn+1  !5  'F3'
                     IF(Topo(i,j-1,k)%Topo%fzxnr>0) THEN
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j-1,k)%Topo%fzxnr 
                      ELSE
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     END IF
                     ifn=ifn+1  !6  'F5'
                     IF(Topo(i-1,j,k)%Topo%fyznr>0) THEN
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i-1,j,k)%Topo%fyznr 
                      ELSE
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     END IF 
                !Write(10,*) "ijk                    ",i,j,k,"-allok-vc=0", &
                !       &  "tropos", Topo(i,j,k)%Topo%ctp
                     Write(10,*) "cv",Topo(i,j,k)%Topo%cnr, zw_leerzei, 6, &
                              &  zw_leerzei, tropos_border 
                     Write(10,*) leerzei_acht, &
                                & Topo(i,j,k)%Topo%cfn(4),Topo(i,j,k)%Topo%cfn(2), &
                                & Topo(i,j,k)%Topo%cfn(3),Topo(i,j,k)%Topo%cfn(6), &
                                & Topo(i,j,k)%Topo%cfn(5),Topo(i,j,k)%Topo%cfn(1)
                     nr_cells_write=nr_cells_write+1
                   ELSE  !  IF(Cell(i,j,k)%Cell%vc>0) THEN
                     Topo(i,j,k)%Topo%cfn(1)=Topo(i,j,k)%Topo%FCutNr 
                     ifn=1
                     IF(Topo(i,j,k)%Topo%fxynr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fxynr 
                     END IF
                     IF(Topo(i,j,k)%Topo%fyznr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fyznr 
                     END IF 
                     IF(Topo(i,j,k)%Topo%fzxnr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fzxnr 
                     END IF
                     IF(Topo(i,j,k-1)%Topo%fxynr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k-1)%Topo%fxynr 
                     END IF
                     IF(Topo(i,j-1,k)%Topo%fzxnr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j-1,k)%Topo%fzxnr 
                     END IF
                     IF(Topo(i-1,j,k)%Topo%fyznr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%cfn(ifn)=Topo(i-1,j,k)%Topo%fyznr 
                     END IF 
                     Topo(i,j,k)%Topo%nfn=ifn
                 ! Write(10,*) "ijk                       ",i,j,k," vc > 0 = ", &
                 !           &    Cell(i,j,k)%Cell%vc, "tropos", Topo(i,j,k)%Topo%ctp
                     Write(10,*) "cc",Topo(i,j,k)%Topo%cnr, zw_leerzei, &
                             &  Topo(i,j,k)%Topo%nfn, zw_leerzei, tropos_border
                     Write(10,*) leerzei_acht, &
                                & (Topo(i,j,k)%Topo%cfn(ifn),ifn=1,Topo(i,j,k)%Topo%nfn)
                     nr_cells_write=nr_cells_write+1
                     !CALL WriteCellProt2(Cell(i,j,k)%Cell,i,j,k)
                   END IF  !if(vc==0) then ... else
               END IF
            END IF
          END DO
        END DO
      END DO
      DO i=ix0+1,ix1
        DO j=iy0+1,iy1
          DO k=iz0+1,iz1
            IF(.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN !!! .NOT. ASSOCIATED(Cell)
              in_out=Vertices(i-1,j-1,k-1)%in_out &
                    +Vertices(i,j-1,k-1)%in_out &
                    +Vertices(i-1,j,k-1)%in_out &
                    +Vertices(i-1,j-1,k)%in_out &
                    +Vertices(i-1,j,k)%in_out &
                    +Vertices(i,j-1,k)%in_out &
                    +Vertices(i,j,k-1)%in_out &
                    +Vertices(i,j,k)%in_out
              IF (in_out>=0) THEN
               !!IF (Cell(i,j,k)%Cell%in_out>5) THEN
                !Reihenfolge F1,F6,F4,F5,F3,F2; (F1,dann gegen Uhrzeiger,F2)
     !  Write(10,*) "ijk                       ",i,j,k," not associated" , &
     !           &  "    tropos", Topo(i,j,k)%Topo%ctp 
                Write(10,*) "ch",Topo(i,j,k)%Topo%cnr, zw_leerzei, 6, &
                         &  zw_leerzei, tropos_inside 
                Write(10,*) leerzei_acht, &
                          & Topo(i,j,k-1)%Topo%fxynr,Topo(i,j,k)%Topo%fyznr, &
                          & Topo(i,j,k)%Topo%fzxnr,Topo(i-1,j,k)%Topo%fyznr, &
                          & Topo(i,j-1,k)%Topo%fzxnr,Topo(i,j,k)%Topo%fxynr
                nr_cells_write=nr_cells_write+1
              ELSE
                !Topo(i,j,k)%Topo%ctp=-1  ! unterhalb Berg
              END IF
            END IF
          END DO
        END DO
      END DO
      !-- Output all inside mountain cells
      DO i=ix0+1,ix1
        DO j=iy0+1,iy1
          DO k=iz0+1,iz1
            IF(.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN !!! .NOT. ASSOCIATED(Cell)
               IF(Topo(i,j,k)%Topo%icnr>0) THEN
          !Write(10,*) "ijk                    ",i,j,k," .not. associated", &
          !         &  "  inside mountine" , Topo(i,j,k)%Topo%ctp
                 Write(10,*) "ch",Topo(i,j,k)%Topo%icnr, zw_leerzei, 6, &
                          &  zw_leerzei, soil_inside 
                 Write(10,*) leerzei_acht, &
                          & Topo(i,j,k-1)%Topo%ifxynr,Topo(i,j,k)%Topo%ifyznr, &
                          & Topo(i,j,k)%Topo%ifzxnr,Topo(i-1,j,k)%Topo%ifyznr, &
                          & Topo(i,j-1,k)%Topo%ifzxnr,Topo(i,j,k)%Topo%ifxynr
                 nr_cells_write=nr_cells_write+1
               END IF
            ELSE
               IF(Topo(i,j,k)%Topo%icnr>0) THEN
                  IF(Cell(i,j,k)%Cell%vc>0) THEN
                     !Reihenfolge FCut,F2,F6,F4,F1,F3,F5; 
                     Topo(i,j,k)%Topo%icfn(1)=Topo(i,j,k)%Topo%FCutNr 
                     ifn=1
                     IF(Topo(i,j,k)%Topo%ifxynr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%icfn(ifn)=Topo(i,j,k)%Topo%ifxynr 
                     END IF
                     IF(Topo(i,j,k)%Topo%ifyznr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%icfn(ifn)=Topo(i,j,k)%Topo%ifyznr 
                     END IF 
                     IF(Topo(i,j,k)%Topo%ifzxnr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%icfn(ifn)=Topo(i,j,k)%Topo%ifzxnr 
                     END IF
                     IF(Topo(i,j,k-1)%Topo%ifxynr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%icfn(ifn)=Topo(i,j,k-1)%Topo%ifxynr 
                     END IF
                     IF(Topo(i,j-1,k)%Topo%ifzxnr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%icfn(ifn)=Topo(i,j-1,k)%Topo%ifzxnr 
                     END IF
                     IF(Topo(i-1,j,k)%Topo%ifyznr>0) THEN
                       ifn=ifn+1
                       Topo(i,j,k)%Topo%icfn(ifn)=Topo(i-1,j,k)%Topo%ifyznr 
                     END IF 
                     Topo(i,j,k)%Topo%infn=ifn
          !Write(10,*) "ijk                    ",i,j,k," vc > 0 = ", &
          !            &    Cell(i,j,k)%Cell%vc, "  inside mountine", Topo(i,j,k)%Topo%ctp
                     Write(10,*) "cc",Topo(i,j,k)%Topo%icnr, zw_leerzei, &
                              &  Topo(i,j,k)%Topo%infn, zw_leerzei, soil_inside
                     Write(10,*) leerzei_acht, &
                                & (Topo(i,j,k)%Topo%icfn(ifn),ifn=1,Topo(i,j,k)%Topo%infn)
                     nr_cells_write=nr_cells_write+1
                     !CALL WriteCellProt2(Cell(i,j,k)%Cell,i,j,k)

                  ELSE IF (Cell(i,j,k)%Cell%vc==0.AND.Cell(i,j,k)%Cell%in_out<=-4) THEN
                     !Write(10,*) "Celle vc==0 ch-special                " , i, j,k
                     !Write(10,*) "              Cell(i,j,k)%Cell%in_out=", &
                     !    & Cell(i,j,k)%Cell%in_out    
                     IF(Cell(i,j,k)%Cell%in_out==-7) THEN
          ! Write(10,*) "ijk                    ",i,j,k," vc = 0" , &
          !            &    Cell(i,j,k)%Cell%vc, "in_out=-7  inside mountine", Topo(i,j,k)%Topo%ctp
                       Write(10,*) "ch",Topo(i,j,k)%Topo%icnr, zw_leerzei, 6, &
                                &  zw_leerzei, soil_inside
                       Write(10,*) leerzei_acht, &
                          & Topo(i,j,k-1)%Topo%ifxynr,Topo(i,j,k)%Topo%ifyznr, &
                          & Topo(i,j,k)%Topo%ifzxnr,Topo(i-1,j,k)%Topo%ifyznr, &
                          & Topo(i,j-1,k)%Topo%ifzxnr,Topo(i,j,k)%Topo%ifxynr
                     nr_cells_write=nr_cells_write+1
                     ELSE IF(Cell(i,j,k)%Cell%in_out==-6) THEN
          ! Write(10,*) "ijk                    ",i,j,k," vc = 0" , &
          !            &    Cell(i,j,k)%Cell%vc, "in_out=-6  inside mountine", Topo(i,j,k)%Topo%ctp
                       Write(10,*) "ch",Topo(i,j,k)%Topo%icnr, zw_leerzei, 6, &
                                &  zw_leerzei, soil_inside
                       Write(10,*) leerzei_acht, &
                          & Topo(i,j,k-1)%Topo%ifxynr,Topo(i,j,k)%Topo%ifyznr, &
                          & Topo(i,j,k)%Topo%ifzxnr,Topo(i-1,j,k)%Topo%ifyznr, &
                          & Topo(i,j-1,k)%Topo%ifzxnr,Topo(i,j,k)%Topo%ifxynr
                     nr_cells_write=nr_cells_write+1
                     ELSE IF (Cell(i,j,k)%Cell%in_out==-4) THEN
                       ! Einfügen Face-grenzen  XY-F2,YZ-ZX-F3-F4 splitten  
          !Write(10,*) "ijk                    ",i,j,k," vc = 0" , &
          !            &    Cell(i,j,k)%Cell%vc, "in_out=-4  inside mountine", Topo(i,j,k)%Topo%ctp
                       Write(10,*) "ch",Topo(i,j,k)%Topo%icnr, zw_leerzei, 6, &
                                &  zw_leerzei, soil_inside
                       Write(10,*) leerzei_acht, &
                          & Topo(i,j,k-1)%Topo%ifxynr,Topo(i,j,k)%Topo%ifyznr, &
                          & Topo(i,j,k)%Topo%ifzxnr,Topo(i-1,j,k)%Topo%ifyznr, &
                          & Topo(i,j-1,k)%Topo%ifzxnr,Topo(i,j,k)%Topo%ifxynr
                     nr_cells_write=nr_cells_write+1
                     END IF
                  END IF ! vc>0 .or. (vc==0,Cell%in_out<=-4)
               END IF ! icnr>0
            END IF  ! if/else not associated cell
          END DO
        END DO
      END DO
      ! ende inside Cells
  END DO  !ib
  Write(10,*) "Ende"


  !------------------ Protokoll zur Ausgabe --------------------------------! 
  Write(OutUnitProt,*) leerzei3,"> Protokoll aus SUBROUTINE WriteTropoOroCellsParaV()"
  Write(OutUnitProt,*) leerzei3,"  Schreiben in  Datei  *.pva.gall"
  Write(OutUnitProt,*) leerzei3,"--------------------------------------------------"
  IF (nr_edges_write /= para_grid_enr ) THEN
    Write(OutUnitProt,*) leerzei5,"analysierte Edges /= gelistete Faces", &
                &         para_grid_enr,"/=",nr_edges_write
  ELSE 
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Edges : analysierte == gelistete "
  END IF
  IF (nr_faces_write /= para_grid_fnr ) THEN
    Write(OutUnitProt,*) leerzei5,"analysierte Faces /= gelistete Faces", &
                &         para_grid_fnr,"/=",nr_faces_write
  ELSE
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Faces : analysierte == gelistete "
  END IF  
  IF (nr_cells_write /= para_grid_cnr ) THEN 
    Write(OutUnitProt,*) leerzei5,"analysierte Cellen /= gelistete Cellen ", &
                &        para_grid_cnr," /=",nr_cells_write        
  ELSE
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Cellen: analysierte == gelistete "
  END IF
  Write(OutUnitProt,*) ""
    !------------------ Ende Protokoll zur Ausgabe --------------------------! 
!!    WRITE(10,'(a8,i8)') cells,nr_cells
!!    DO ib=1,nb
!!      Write(*,*) "                                     Block : ",ib,"\/",nb
!!      CALL Set(Floor(ib))
!!      DO k=iz0+1,iz1
!!        DO j=iy0+1,iy1
!!          DO i=ix0+1,ix1
!!            !CALL WriteCellGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!            IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!              Write(10,*) "...................................."
!!              Write(10,*) "Cell(",i,",",j,",",k,")%Cell% ->"
!!              Write(10,*) "  Face1%mp=",Cell(i,j,k)%Cell%Face1%mp
!!              Write(10,*) "  Face2%mp=",Cell(i,j,k)%Cell%Face2%mp
!!              Write(10,*) "  Face3%mp=",Cell(i,j,k)%Cell%Face3%mp
!!              Write(10,*) "  Face4%mp=",Cell(i,j,k)%Cell%Face4%mp
!!              Write(10,*) "  Face5%mp=",Cell(i,j,k)%Cell%Face5%mp
!!              Write(10,*) "  Face6%mp=",Cell(i,j,k)%Cell%Face6%mp
!!              Write(10,*) "  vc=",Cell(i,j,k)%Cell%vc
!!              WRITE(10,*) (Cell(i,j,k)%Cell%VertCut(ic),ic=1,Cell(i,j,k)%Cell%vc)
!!              Write(10,*) "...................................."
!!            END IF
!!          END DO
!!        END DO
!!      END DO
!!    END DO  !ib
  CALL WriteEndOutParaVToProt()

  CLOSE(10)
END SUBROUTINE WriteTropoOroCellsParaV

END MODULE OutputParaView_Mod
