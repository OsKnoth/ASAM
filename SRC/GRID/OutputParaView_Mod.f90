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

SUBROUTINE WriteFaceProt2(Face,i,j,k)
  TYPE(Face_T), POINTER :: Face
  INTEGER :: i,j,k
  WRITE(OutUnitProt,*) ' Face%--->'
  WRITE(OutUnitProt,*) "    !(TYPE Vertex_T:", "  Point(x,y,z)    in_out   nrP   nrInP   nrCutP   Shift)"
  WRITE(OutUnitProt,*) "    Edge1%Vert1:",Face%Edge1%Vert1
  WRITE(OutUnitProt,*) "    Edge1%Vert2:",Face%Edge1%Vert2
  WRITE(OutUnitProt,*) "    Edge3%Vert1:",Face%Edge3%Vert1
  WRITE(OutUnitProt,*) "    Edge3%Vert2:",Face%Edge3%Vert2
  WRITE(OutUnitProt,*) "    ..........."
  WRITE(OutUnitProt,*) '    Vol        =',Face%Vol
  WRITE(OutUnitProt,*) '    NumberVert =',Face%NumberVert
  WRITE(OutUnitProt,*) '    VertexList =',Face%VertexList(1:Face%NumberVert)
  WRITE(OutUnitProt,*) '    in_out     =',Face%in_out
  WRITE(OutUnitProt,*) '    ec         =',Face%ec
  WRITE(OutUnitProt,*) '    EdgeCut    =',Face%EdgeCut
END SUBROUTINE WriteFaceProt2

SUBROUTINE WriteCellProt2(Cell,i,j,k)

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
    WRITE(OutUnitProt,*) '-------------- WriteCellProt2() ---------------------------------'
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
    WRITE(OutUnitProt,*) '    EdgeCut   =',Cell%VertCut(1),Cell%VertCut(2)
    WRITE(OutUnitProt,*) '    CutF_MidP =',Cell%CutF_MidP
    WRITE(OutUnitProt,*) '    Face1%Vol =',Cell%Face1%Vol
    WRITE(OutUnitProt,*) '    Face2%Vol =',Cell%Face2%Vol
    WRITE(OutUnitProt,*) '    Face3%Vol =',Cell%Face3%Vol
    WRITE(OutUnitProt,*) '    Face4%Vol =',Cell%Face4%Vol
    WRITE(OutUnitProt,*) '    Face5%Vol =',Cell%Face5%Vol
    WRITE(OutUnitProt,*) '    Face6%Vol =',Cell%Face6%Vol
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceXY','  (F1)'
    CALL WriteFaceProt2(Cell%Face1,i,j,k-1)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceXY','  (F2)'
    CALL WriteFaceProt2(Cell%Face2,i,j,k)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceXZ','  (F3)'
    CALL WriteFaceProt2(Cell%Face3,i,j-1,k)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceXZ','  (F4)'
    CALL WriteFaceProt2(Cell%Face4,i,j,k)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceYZ','  (F5)'
    CALL WriteFaceProt2(Cell%Face5,i-1,j,k)
    WRITE(OutUnitProt,*) '............'
    WRITE(OutUnitProt,*) ' --> FaceYZ','  (F6)'
    CALL WriteFaceProt2(Cell%Face6,i,j,k)
    WRITE(OutUnitProt,*) 'Cell','(',i,',',j,',',k,')'
    WRITE(OutUnitProt,*) '-------------- WriteCellProt2() Ende ----------------------------'
  END IF

END SUBROUTINE WriteCellProt2



!!AusblendenSUBROUTINE SetTropoCellsParaV(FileName)
!!Ausblenden  CHARACTER*50 :: FileName
!!Ausblenden  INTEGER :: ib,i,j,k,in_out
!!Ausblenden  INTEGER :: iec,ecnr,exnr,eynr,eznr
!!Ausblenden  INTEGER :: fcnr,fxynr,fzxnr,fyznr,fnr
!!Ausblenden  INTEGER :: ic,cnr,cnr_cut
!!Ausblenden  TYPE (Edge_T) :: Edg(1:4)
!!Ausblenden  LOGICAL :: s_ec
!!Ausblenden  INTEGER :: li
!!Ausblenden
!!Ausblenden  WRITE(*,*) ">....... Set Troposphere for ParaView"
!!Ausblenden  IF(out_area=='y') THEN
!!Ausblenden    WRITE(OutUnitProt,'(a8,i8)') cells,nr_viewcells
!!Ausblenden    DO ib=1,1 !nb
!!Ausblenden!      Write(*,*) "                                     Block : ",ib,"\/",nb
!!Ausblenden      Write(*,*) "                                      z.Zt.: nur mit 1 Block"  
!!Ausblenden      CALL Set(Floor(ib))
!!Ausblenden      CALL Search_OutArea
!!Ausblenden      IF(view_cell=='y') THEN
!!Ausblenden        ! Für Analyze aktiv schalten, Ausgabe in *.prot-File
!!Ausblenden        !DO i=ix0+1,ix1
!!Ausblenden        !  DO j=iy0+1,iy1
!!Ausblenden        !    DO k=iz0+1,iz1
!!Ausblenden        !      !WRITE(OutUnitProt,*) "Cell(i,j,k)", i, j, k,"   WriteAllCellsAsciiGMV"
!!Ausblenden        !      !WRITE(*,*) "Cell(i,j,k)", i, j, k,"   WriteAllCellsAsciiGMV"
!!Ausblenden        !      !CALL WriteCellTopoAOut(Cell(i,j,k)%Cell,i,j,k)
!!Ausblenden        !      IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden        !       IF (i==2 .and. j==3 .and. k==1) THEN
!!Ausblenden        !        CALL WriteCellProt2(Cell(i,j,k)%Cell,i,j,k)
!!Ausblenden        !       END IF
!!Ausblenden        !      END IF
!!Ausblenden        !    END DO
!!Ausblenden        !  END DO
!!Ausblenden        !END DO
!!Ausblenden        cnr_cut=0
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
!!Ausblenden                    .OR.Cell(i,j,k)%Cell%in_out==8) THEN
!!Ausblenden                   !IF(Cell(i,j,k)%Cell%vc==0) THEN 
!!Ausblenden                   !ELSE ! vc>0 
!!Ausblenden                   !END IF   !beide eingeschossen
!!Ausblenden                  Topo(i,j,k)%Topo%ctp=0 
!!Ausblenden                  cnr_cut=cnr_cut+1
!!Ausblenden                  Topo(i,j,k)%Topo%cnr=cnr_cut
!!Ausblenden                END IF
!!Ausblenden              !ELSE !!! .NOT. ASSOCIATED(Cell)
!!Ausblenden              !  Betreff Ausgabe-Reihenfolge 
!!Ausblenden              !          extra Schleife für Cellen oberhalb
!!Ausblenden              !  END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        cnr=cnr_cut
!!Ausblenden        topo_cvcnr=cnr_cut
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                in_out=Vertices(i-1,j-1,k-1)%in_out &
!!Ausblenden                      +Vertices(i,j-1,k-1)%in_out &
!!Ausblenden                      +Vertices(i-1,j,k-1)%in_out &
!!Ausblenden                      +Vertices(i-1,j-1,k)%in_out &
!!Ausblenden                      +Vertices(i-1,j,k)%in_out &
!!Ausblenden                      +Vertices(i,j-1,k)%in_out &
!!Ausblenden                      +Vertices(i,j,k-1)%in_out &
!!Ausblenden                      +Vertices(i,j,k)%in_out
!!Ausblenden                IF (in_out>=0) THEN
!!Ausblenden                 !!IF (Cell(i,j,k)%Cell%in_out>5) THEN
!!Ausblenden                  Topo(i,j,k)%Topo%ctp=1
!!Ausblenden                  cnr=cnr+1
!!Ausblenden                  Topo(i,j,k)%Topo%cnr=cnr
!!Ausblenden                ELSE
!!Ausblenden                  Topo(i,j,k)%Topo%ctp=-1
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        topo_cnr=cnr 
!!Ausblenden        !--/EdgeX/EdgeY/EdgeZ/EdgeCut
!!Ausblenden        exnr=0            ! 1.Edge_X
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
!!Ausblenden                 IF(MAX(Edges_X(i,j,k)%Edge%Vert1%in_out,Edges_X(i,j,k)%Edge%Vert2%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Edges_X(i,j,k)%Edge%Vert1%in_out==Edges_X(i,j,k)%Edge%Vert2%in_out.AND. &
!!Ausblenden                    Edges_X(i,j,k)%Edge%Vert1%in_out==0)) THEN
!!Ausblenden                   exnr=exnr+1
!!Ausblenden                   Topo(i,j,k)%Topo%exnr=exnr
!!Ausblenden                   Edges_X(i,j,k)%Edge%eg_nr=exnr
!!Ausblenden                 END IF
!!Ausblenden              ELSE
!!Ausblenden                IF(MAX(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Vertices(i-1,j,k)%in_out==Vertices(i,j,k)%in_out .AND. &
!!Ausblenden                    Vertices(i,j,k)%in_out==0)) THEN
!!Ausblenden                  exnr=exnr+1
!!Ausblenden                  Topo(i,j,k)%Topo%exnr=exnr
!!Ausblenden                END IF 
!!Ausblenden              END IF
!!Ausblenden!              Write(*,*) "exnr= ",exnr, Topo(i,j,k)%Topo%exnr,"  "," i, j, k :", i,j,k 
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        eynr=exnr   ! 2.Edge_Y
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
!!Ausblenden                IF(MAX(Edges_Y(i,j,k)%Edge%Vert1%in_out,Edges_Y(i,j,k)%Edge%Vert2%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Edges_Y(i,j,k)%Edge%Vert1%in_out==Edges_Y(i,j,k)%Edge%Vert2%in_out.AND. &
!!Ausblenden                    Edges_Y(i,j,k)%Edge%Vert1%in_out==0)) THEN
!!Ausblenden                      eynr=eynr+1
!!Ausblenden                      Topo(i,j,k)%Topo%eynr=eynr
!!Ausblenden                      Edges_Y(i,j,k)%Edge%eg_nr=eynr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                IF(MAX(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Vertices(i,j-1,k)%in_out==Vertices(i,j,k)%in_out .AND. &
!!Ausblenden                    Vertices(i,j,k)%in_out==0)) THEN
!!Ausblenden                      eynr=eynr+1
!!Ausblenden                      Topo(i,j,k)%Topo%eynr=eynr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        eznr=eynr   ! 3.Edge_Z
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
!!Ausblenden                IF(MAX(Edges_Z(i,j,k)%Edge%Vert1%in_out,Edges_Z(i,j,k)%Edge%Vert2%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Edges_Z(i,j,k)%Edge%Vert1%in_out==Edges_Z(i,j,k)%Edge%Vert2%in_out.AND. &
!!Ausblenden                    Edges_Z(i,j,k)%Edge%Vert1%in_out==0)) THEN
!!Ausblenden                       eznr=eznr+1
!!Ausblenden                       Topo(i,j,k)%Topo%eznr=eznr
!!Ausblenden                       Edges_Z(i,j,k)%Edge%eg_nr=eznr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                IF(MAX(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Vertices(i,j,k-1)%in_out==Vertices(i,j,k)%in_out.AND. &
!!Ausblenden                    Vertices(i,j,k)%in_out==0)) THEN
!!Ausblenden                       eznr=eznr+1
!!Ausblenden                       Topo(i,j,k)%Topo%eznr=eznr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        ecnr=eznr   ! 4.EdgeCut
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                IF(Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                  !Edge-Cut zuordnen!!!
!!Ausblenden                  Topo(i,j,k)%Topo%fceg_nr=Cell(i,j,k)%Cell%vc
!!Ausblenden                  !----------------------------------------------------------------------------
!!Ausblenden                  !Analyze geschnittene Face_XY-/Face_YZ-/Face_ZX-; ecnr zuordnen                 
!!Ausblenden                  !Liste Topo(i,j,k)%Topo%fceg_nrcut(1:Cell(i,j,k)%Cell%vc) zuordnen
!!Ausblenden                  !----------------------------------------------------------------------------
!!Ausblenden                  DO iec=1,Cell(i,j,k)%Cell%vc-1
!!Ausblenden                    !ueber die VertCutListe die Vert-Nummern ermitteln
!!Ausblenden                    !im Output
!!Ausblenden                    !IF(Faces_XY(i,j,k-1)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1))  .OR. &
!!Ausblenden                       (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_XY(i,j,k-1)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_XY(i,j,k-1)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_XY(i,j,k-1)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_XY(i,j,k-1)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_XY(i,j,k-1)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k-1)%Topo%ecnrxy=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k-1)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_XY(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_XY(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_XY(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_XY(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_XY(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                           ecnr=ecnr+1
!!Ausblenden                           Faces_XY(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                           Topo(i,j,k)%Topo%ecnrxy=ecnr
!!Ausblenden                           Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                           Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_YZ(i-1,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) )THEN
!!Ausblenden                        Edg(1)=Faces_YZ(i-1,j,k)%Face%Edge1
!!Ausblenden                        Edg(2)=Faces_YZ(i-1,j,k)%Face%Edge2
!!Ausblenden                        Edg(3)=Faces_YZ(i-1,j,k)%Face%Edge3
!!Ausblenden                        Edg(4)=Faces_YZ(i-1,j,k)%Face%Edge4
!!Ausblenden                        s_ec=.TRUE.
!!Ausblenden                        DO li=1,4
!!Ausblenden                          IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                              Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                              Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                          END IF
!!Ausblenden                        END DO
!!Ausblenden                        IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                             ecnr=ecnr+1
!!Ausblenden                             Faces_YZ(i-1,j,k)%Face%egc_nr=ecnr
!!Ausblenden                             Topo(i-1,j,k)%Topo%ecnryz=ecnr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                        ELSE IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i-1,j,k)%Face%egc_nr
!!Ausblenden                        END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_YZ(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_YZ(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_YZ(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_YZ(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                              Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_YZ(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%ecnryz=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_ZX(i,j-1,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_ZX(i,j-1,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_ZX(i,j-1,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_ZX(i,j-1,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_ZX(i,j-1,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_ZX(i,j-1,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j-1,k)%Topo%ecnrzx=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j-1,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_ZX(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_ZX(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_ZX(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_ZX(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_ZX(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%ecnrzx=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                  END DO
!!Ausblenden                    !IF(Faces_XY(i,j,k-1)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
!!Ausblenden                       Edg(1)=Faces_XY(i,j,k-1)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_XY(i,j,k-1)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_XY(i,j,k-1)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_XY(i,j,k-1)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_XY(i,j,k-1)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k-1)%Topo%ecnrxy=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k-1)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_XY(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )    .OR. & 
!!Ausblenden                       (Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
!!Ausblenden                       Edg(1)=Faces_XY(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_XY(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_XY(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_XY(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                           ecnr=ecnr+1
!!Ausblenden                           Faces_XY(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                           Topo(i,j,k)%Topo%ecnrxy=ecnr
!!Ausblenden                           Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                           Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_YZ(i-1,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
!!Ausblenden                        Edg(1)=Faces_YZ(i-1,j,k)%Face%Edge1
!!Ausblenden                        Edg(2)=Faces_YZ(i-1,j,k)%Face%Edge2
!!Ausblenden                        Edg(3)=Faces_YZ(i-1,j,k)%Face%Edge3
!!Ausblenden                        Edg(4)=Faces_YZ(i-1,j,k)%Face%Edge4
!!Ausblenden                        s_ec=.TRUE.
!!Ausblenden                        DO li=1,4
!!Ausblenden                          IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                              Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                              Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                          END IF
!!Ausblenden                        END DO
!!Ausblenden                        IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                             ecnr=ecnr+1
!!Ausblenden                             Faces_YZ(i-1,j,k)%Face%egc_nr=ecnr
!!Ausblenden                             Topo(i-1,j,k)%Topo%ecnryz=ecnr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                        ELSE IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i-1,j,k)%Face%egc_nr
!!Ausblenden                        END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
!!Ausblenden                       Edg(1)=Faces_YZ(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_YZ(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_YZ(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_YZ(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_YZ(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%ecnryz=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_ZX(i,j-1,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
!!Ausblenden                       Edg(1)=Faces_ZX(i,j-1,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_ZX(i,j-1,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_ZX(i,j-1,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_ZX(i,j-1,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_ZX(i,j-1,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j-1,k)%Topo%ecnrzx=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j-1,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
!!Ausblenden                       Edg(1)=Faces_ZX(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_ZX(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_ZX(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_ZX(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_ZX(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%ecnrzx=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        topo_enr=ecnr
!!Ausblenden        topo_ecnr=topo_enr-eznr
!!Ausblenden
!!Ausblenden!        Write(*,*) "topo_enr:         ", "enr-summe"
!!Ausblenden!        Write(*,*) "exnr =",exnr,      "  ",exnr
!!Ausblenden!        Write(*,*) "eynr =",eynr-exnr, "  ",eynr
!!Ausblenden!        Write(*,*) "eznr =",eznr-eynr, "  ",eznr
!!Ausblenden!        Write(*,*) "ecnr =",topo_enr-eznr, "  ", topo_enr
!!Ausblenden        !FaceCut,Faces_XY,Faces_XZ,Faces_YZ
!!Ausblenden        fcnr=0
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                IF(Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                    fcnr=fcnr+1
!!Ausblenden                    Topo(i,j,k)%Topo%FCutNr=fcnr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden!!        Write(*,*) "fcnr=", fcnr, " --> anz. fcuts"
!!Ausblenden        fxynr=fcnr
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
!!Ausblenden                IF(Vertices(i-1,j-1,k)%in_out==0 .and. &
!!Ausblenden                   Vertices(i  ,j-1,k)%in_out==0 .and. &
!!Ausblenden                   Vertices(i-1,j  ,k)%in_out==0 .and. &
!!Ausblenden                   Vertices(i  ,j  ,k)%in_out==0) THEN
!!Ausblenden                   !Special-Grenzfläche als 'fc' ausgewertet,
!!Ausblenden                   !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
!!Ausblenden                   !Write(10,*) "Faces_XY(i,j,k)%Face Special Face Grenze als fc",  i,j,k
!!Ausblenden                   IF(Cell(i,j,k)%Cell%vc==0.and.Cell(i,j,k+1)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fxynr=Topo(i,j,k+1)%Topo%FCutNr
!!Ausblenden                      Topo(i,j,k)%Topo%FCutNr=Topo(i,j,k+1)%Topo%FCutNr
!!Ausblenden                   ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fxynr=Topo(i,j,k)%Topo%FCutNr
!!Ausblenden                   ELSE 
!!Ausblenden                      Write(OutUnitProt,*)"Warning aus SetTropoParaView(): "
!!Ausblenden                      Write(OutUnitProt,*)"  Face_XY(i,j,k)-Special-Neu gefunden!"  ,i,j,k
!!Ausblenden                      CALL WriteFaceProt2(Faces_XY(i,j,k)%Face)
!!Ausblenden                   END IF
!!Ausblenden                ELSE IF(Faces_XY(i,j,k)%Face%in_out>-4.and.Faces_XY(i,j,k)%Face%numberVert>2) THEN
!!Ausblenden                   fxynr=fxynr+1 
!!Ausblenden                   Topo(i,j,k)%Topo%fxynr=fxynr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                in_out=Vertices(i-1,j-1,k)%in_out+Vertices(i,j-1,k)%in_out+ &
!!Ausblenden                       Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
!!Ausblenden                IF(in_out>=0) THEN
!!Ausblenden                    fxynr=fxynr+1 
!!Ausblenden                    Topo(i,j,k)%Topo%fxynr=fxynr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        fzxnr=fxynr
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
!!Ausblenden                IF(Vertices(i-1,j,k-1)%in_out==0 .and. &
!!Ausblenden                   Vertices(i  ,j,k-1)%in_out==0 .and. & 
!!Ausblenden                   Vertices(i-1,j,k  )%in_out==0 .and. &
!!Ausblenden                   Vertices(i  ,j,k  )%in_out==0 ) THEN
!!Ausblenden                   !Special Grenzfläche als 'fc' ausgewertet,
!!Ausblenden                   !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
!!Ausblenden                   IF(Cell(i,j,k)%Cell%vc==0.and.Cell(i,j+1,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fzxnr=Topo(i,j+1,k)%Topo%FCutNr
!!Ausblenden                      Topo(i,j,k)%Topo%FCutNr=Topo(i,j+1,k)%Topo%FCutNr
!!Ausblenden                   ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fzxnr=Topo(i,j,k)%Topo%FCutNr
!!Ausblenden                   ELSE 
!!Ausblenden                      Write(OutUnitProt,*)"Warning aus SetTropoParaView(): "
!!Ausblenden                      Write(OutUnitProt,*)"  Face_ZX(i,j,k)-Special-Neu gefunden! "  ,i,j,k
!!Ausblenden                      CALL WriteFaceProt2(Faces_ZX(i,j,k)%Face)
!!Ausblenden                   END IF
!!Ausblenden                ELSE IF(Faces_ZX(i,j,k)%Face%in_out>-4.and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
!!Ausblenden                   fzxnr=fzxnr+1 
!!Ausblenden                   Topo(i,j,k)%Topo%fzxnr=fzxnr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                in_out=Vertices(i-1,j,k-1)%in_out+Vertices(i,j,k-1)%in_out+ &
!!Ausblenden                       Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
!!Ausblenden                IF(in_out>=0) THEN
!!Ausblenden                  fzxnr=fzxnr+1
!!Ausblenden                  Topo(i,j,k)%Topo%fzxnr=fzxnr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        fyznr=fzxnr
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
!!Ausblenden                IF(Vertices(i,j-1,k-1)%in_out==0 .and. &
!!Ausblenden                   Vertices(i,j-1,k  )%in_out==0 .and. &
!!Ausblenden                   Vertices(i,j  ,k  )%in_out==0 .and. &
!!Ausblenden                   Vertices(i,j  ,k-1)%in_out==0 ) THEN
!!Ausblenden                   !Special Grenzfläche als 'fc' ausgewertet,
!!Ausblenden                   !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
!!Ausblenden                   IF(Cell(i,j,k)%Cell%vc==0.and.Cell(i+1,j,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fyznr=Topo(i+1,j,k)%Topo%FCutNr
!!Ausblenden                      Topo(i,j,k)%Topo%FCutNr=Topo(i+1,j,k)%Topo%FCutNr
!!Ausblenden                   ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fyznr=Topo(i,j,k)%Topo%FCutNr
!!Ausblenden                   ELSE 
!!Ausblenden                      Write(OutUnitProt,*)"Warning aus SetTropoParaView(): "
!!Ausblenden                      Write(OutUnitProt,*)"  Face_YZ(i,j,k)-Special-Neu gefunden! "  ,i,j,k
!!Ausblenden                      CALL WriteFaceProt2(Faces_YZ(i,j,k)%Face)
!!Ausblenden                   END IF
!!Ausblenden                ELSE IF(Faces_YZ(i,j,k)%Face%in_out>-4.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
!!Ausblenden                  fyznr=fyznr+1 
!!Ausblenden                  Topo(i,j,k)%Topo%fyznr=fyznr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                in_out=Vertices(i,j-1,k-1)%in_out+Vertices(i,j-1,k)%in_out+ &
!!Ausblenden                       Vertices(i,j,k)%in_out+Vertices(i,j,k-1)%in_out
!!Ausblenden                IF(in_out>=0) THEN
!!Ausblenden                  fyznr=fyznr+1 
!!Ausblenden                  Topo(i,j,k)%Topo%fyznr=fyznr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        topo_fnr=fyznr
!!Ausblenden
!!Ausblenden!!        Write(*,*) "topo_fnr:        ", "fnr-summe"
!!Ausblenden!!        Write(*,*) "fcnr  =",fcnr,        "  ",fcnr
!!Ausblenden!!        Write(*,*) "fxynr =",fxynr-fcnr,  "  ",fxynr
!!Ausblenden!!        Write(*,*) "fzxnr =",fzxnr-fxynr, "  ",fzxnr
!!Ausblenden!!        Write(*,*) "fyznr =",fyznr-fzxnr, "  ",fyznr
!!Ausblenden      END IF   !IF(view_cell=='y')
!!Ausblenden    END DO  !ib
!!Ausblenden  ELSE
!!Ausblenden     Write(*,*) "z.Zt.  SetTropoCellsParaV ","  nur wenn  '#OutputDomain'  gesetzt ist!"  
!!Ausblenden!    WRITE(OutUnitProt,'(a8,i8)') cells,nr_cells
!!Ausblenden!    DO ib=1,nb
!!Ausblenden!      Write(*,*) "                                     Block : ",ib,"\/",nb
!!Ausblenden!      CALL Set(Floor(ib))
!!Ausblenden!      DO k=iz0+1,iz1
!!Ausblenden!        DO j=iy0+1,iy1
!!Ausblenden!          DO i=ix0+1,ix1
!!Ausblenden!            !CALL WriteCellGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!Ausblenden!            IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden!              Write(OutUnitProt,*) "...................................."
!!Ausblenden!              Write(OutUnitProt,*) "Cell(",i,",",j,",",k,")%Cell% ->"
!!Ausblenden!              Write(OutUnitProt,*) "  Face1%mp=",Cell(i,j,k)%Cell%Face1%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face2%mp=",Cell(i,j,k)%Cell%Face2%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face3%mp=",Cell(i,j,k)%Cell%Face3%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face4%mp=",Cell(i,j,k)%Cell%Face4%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face5%mp=",Cell(i,j,k)%Cell%Face5%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face6%mp=",Cell(i,j,k)%Cell%Face6%mp
!!Ausblenden!              Write(OutUnitProt,*) "  vc=",Cell(i,j,k)%Cell%vc
!!Ausblenden!              WRITE(OutUnitProt,*) (Cell(i,j,k)%Cell%VertCut(ic),ic=1,Cell(i,j,k)%Cell%vc)
!!Ausblenden!              Write(OutUnitProt,*) "...................................."
!!Ausblenden!            END IF
!!Ausblenden!          END DO
!!Ausblenden!        END DO
!!Ausblenden!      END DO
!!Ausblenden!    END DO  !ib
!!Ausblenden  END IF
!!Ausblenden
!!AusblendenEND SUBROUTINE SetTropoCellsParaV
!!Ausblenden
!!Ausblenden
!!AusblendenSUBROUTINE SetTropoOroCellsParaV(FileName)
!!Ausblenden  CHARACTER*50 :: FileName
!!Ausblenden  INTEGER :: ib,i,j,k,in_out
!!Ausblenden  INTEGER :: iec,ecnr,exnr,eynr,eznr
!!Ausblenden  INTEGER :: fcnr,fxynr,fzxnr,fyznr,fnr
!!Ausblenden  INTEGER :: ic,cnr,cnr_cut
!!Ausblenden  TYPE (Edge_T) :: Edg(1:4)
!!Ausblenden  LOGICAL :: s_ec
!!Ausblenden  INTEGER :: li
!!Ausblenden
!!Ausblenden  CALL Display_SetTroposOroParaV(FileName)
!!Ausblenden  IF(out_area=='y') THEN
!!Ausblenden    WRITE(OutUnitProt,'(a8,i8)') cells,nr_viewcells
!!Ausblenden    DO ib=1,1 !nb
!!Ausblenden      !Write(*,*) "                                     Block : ",ib,"\/",nb
!!Ausblenden      Write(*,*) "                                      z.Zt.: nur mit 1 Block"  
!!Ausblenden      CALL Set(Floor(ib))
!!Ausblenden      CALL Search_OutArea
!!Ausblenden      IF(view_cell=='y') THEN
!!Ausblenden        ! Für Analyze aktiv schalten, Ausgabe in *.prot-File
!!Ausblenden        !DO i=ix0+1,ix1
!!Ausblenden        !  DO j=iy0+1,iy1
!!Ausblenden        !    DO k=iz0+1,iz1
!!Ausblenden        !      !WRITE(OutUnitProt,*) "Cell(i,j,k)", i, j, k,"   WriteAllCellsAsciiGMV"
!!Ausblenden        !      !WRITE(*,*) "Cell(i,j,k)", i, j, k,"   WriteAllCellsAsciiGMV"
!!Ausblenden        !      !CALL WriteCellTopoAOut(Cell(i,j,k)%Cell,i,j,k)
!!Ausblenden        !      IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden        !       IF (i==2 .and. j==3 .and. k==1) THEN
!!Ausblenden        !        CALL WriteCellProt2(Cell(i,j,k)%Cell,i,j,k)
!!Ausblenden        !       END IF
!!Ausblenden        !      END IF
!!Ausblenden        !    END DO
!!Ausblenden        !  END DO
!!Ausblenden        !END DO
!!Ausblenden        cnr_cut=0
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                 IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
!!Ausblenden                     .OR.Cell(i,j,k)%Cell%in_out==8) THEN
!!Ausblenden                   !IF(Cell(i,j,k)%Cell%vc==0) THEN 
!!Ausblenden                   !ELSE ! vc>0 
!!Ausblenden                   !END IF   !beide eingeschossen
!!Ausblenden                   Topo(i,j,k)%Topo%ctp=0 
!!Ausblenden                   cnr_cut=cnr_cut+1
!!Ausblenden                   Topo(i,j,k)%Topo%cnr=cnr_cut
!!Ausblenden                 END IF
!!Ausblenden              !ELSE !!! .NOT. ASSOCIATED(Cell)
!!Ausblenden              !  Betreff Ausgabe-Reihenfolge 
!!Ausblenden              !          extra Schleife für Cellen oberhalb
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        cnr=cnr_cut
!!Ausblenden        topo_cvcnr=cnr_cut
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                in_out=Vertices(i-1,j-1,k-1)%in_out &
!!Ausblenden                      +Vertices(i,j-1,k-1)%in_out &
!!Ausblenden                      +Vertices(i-1,j,k-1)%in_out &
!!Ausblenden                      +Vertices(i-1,j-1,k)%in_out &
!!Ausblenden                      +Vertices(i-1,j,k)%in_out &
!!Ausblenden                      +Vertices(i,j-1,k)%in_out &
!!Ausblenden                      +Vertices(i,j,k-1)%in_out &
!!Ausblenden                      +Vertices(i,j,k)%in_out
!!Ausblenden                IF (in_out>=0) THEN
!!Ausblenden                 !!IF (Cell(i,j,k)%Cell%in_out>5) THEN
!!Ausblenden                  Topo(i,j,k)%Topo%ctp=1
!!Ausblenden                  cnr=cnr+1
!!Ausblenden                  Topo(i,j,k)%Topo%cnr=cnr
!!Ausblenden                ELSE
!!Ausblenden                  Topo(i,j,k)%Topo%ctp=-1
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        topo_cnr=cnr 
!!Ausblenden        !...................... 
!!Ausblenden        ! to Orography Cellen
!!Ausblenden        !......................
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden               !in_out checken damit nicht durch troposphäre abfrage läuft
!!Ausblenden               IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                 IF (Topo(i,j,k)%Topo%ctp==-1) THEN
!!Ausblenden                   cnr=cnr+1
!!Ausblenden                   Topo(i,j,k)%Topo%icnr=cnr
!!Ausblenden                 END IF
!!Ausblenden               ELSE   ! ASSOCIATED(..Celle)
!!Ausblenden                 IF(Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                   !Write(10,*) "Associated noch Auswertung Celle unterhalb Berg", i,j,k
!!Ausblenden                   !Write(10,*) "                 vc= ",Cell(i,j,k)%Cell%vc, &
!!Ausblenden                   !         &   "             in_out= ",Cell(i,j,k)%Cell%in_out
!!Ausblenden                   !CALL WriteCellProt2(Cell(i,j,k)%Cell,i,j,k)
!!Ausblenden                   cnr=cnr+1
!!Ausblenden                   Topo(i,j,k)%Topo%icnr=cnr
!!Ausblenden                 ELSE IF (Cell(i,j,k)%Cell%vc==0.AND.Cell(i,j,k)%Cell%in_out<=-4) THEN
!!Ausblenden                   !Beachtung Cellen mit Edge-GrenzKante und Face_Grenze
!!Ausblenden                   !Write(10,*) "Associated noch Auswertung Celle unterhalb Berg", i,j,k
!!Ausblenden                   !Write(10,*) "                 vc= ",Cell(i,j,k)%Cell%vc, &
!!Ausblenden                   !        &   "             in_out= ",Cell(i,j,k)%Cell%in_out
!!Ausblenden                   cnr=cnr+1
!!Ausblenden                   Topo(i,j,k)%Topo%icnr=cnr
!!Ausblenden                   SELECT CASE (Cell(i,j,k)%Cell%in_out)
!!Ausblenden                     CASE (-4)
!!Ausblenden                       Topo(i,j,k)%Topo%ctp=-4
!!Ausblenden                     CASE (-6)
!!Ausblenden                       Topo(i,j,k)%Topo%ctp=-6
!!Ausblenden                     CASE (-7)
!!Ausblenden                       Topo(i,j,k)%Topo%ctp=-7
!!Ausblenden                     CASE DEFAULT
!!Ausblenden                       Topo(i,j,k)%Topo%ctp=-9
!!Ausblenden                   END SELECT
!!Ausblenden                 END IF
!!Ausblenden               END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        para_grid_cnr=cnr
!!Ausblenden        para_oro_cnr=para_grid_cnr-topo_cnr
!!Ausblenden        ! End Orography Cells
!!Ausblenden        !----------------------------------------------------------------------------
!!Ausblenden        !--/EdgeX/EdgeY/EdgeZ/EdgeCut
!!Ausblenden        exnr=0            ! 1.Edge_X
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
!!Ausblenden                 IF(MAX(Edges_X(i,j,k)%Edge%Vert1%in_out,Edges_X(i,j,k)%Edge%Vert2%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Edges_X(i,j,k)%Edge%Vert1%in_out==Edges_X(i,j,k)%Edge%Vert2%in_out.AND. &
!!Ausblenden                    Edges_X(i,j,k)%Edge%Vert1%in_out==0)) THEN
!!Ausblenden                   exnr=exnr+1
!!Ausblenden                   Topo(i,j,k)%Topo%exnr=exnr
!!Ausblenden                   Edges_X(i,j,k)%Edge%eg_nr=exnr
!!Ausblenden                 END IF
!!Ausblenden              ELSE
!!Ausblenden                IF(MAX(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Vertices(i-1,j,k)%in_out==Vertices(i,j,k)%in_out .AND. &
!!Ausblenden                    Vertices(i,j,k)%in_out==0)) THEN
!!Ausblenden                  exnr=exnr+1
!!Ausblenden                  Topo(i,j,k)%Topo%exnr=exnr
!!Ausblenden                END IF 
!!Ausblenden              END IF
!!Ausblenden!              Write(*,*) "exnr= ",exnr, Topo(i,j,k)%Topo%exnr,"  "," i, j, k :", i,j,k 
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        eynr=exnr   ! 2.Edge_Y
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
!!Ausblenden                IF(MAX(Edges_Y(i,j,k)%Edge%Vert1%in_out,Edges_Y(i,j,k)%Edge%Vert2%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Edges_Y(i,j,k)%Edge%Vert1%in_out==Edges_Y(i,j,k)%Edge%Vert2%in_out.AND. &
!!Ausblenden                    Edges_Y(i,j,k)%Edge%Vert1%in_out==0)) THEN
!!Ausblenden                      eynr=eynr+1
!!Ausblenden                      Topo(i,j,k)%Topo%eynr=eynr
!!Ausblenden                      Edges_Y(i,j,k)%Edge%eg_nr=eynr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                IF(MAX(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Vertices(i,j-1,k)%in_out==Vertices(i,j,k)%in_out .AND. &
!!Ausblenden                    Vertices(i,j,k)%in_out==0)) THEN
!!Ausblenden                      eynr=eynr+1
!!Ausblenden                      Topo(i,j,k)%Topo%eynr=eynr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        eznr=eynr   ! 3.Edge_Z
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
!!Ausblenden                IF(MAX(Edges_Z(i,j,k)%Edge%Vert1%in_out,Edges_Z(i,j,k)%Edge%Vert2%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Edges_Z(i,j,k)%Edge%Vert1%in_out==Edges_Z(i,j,k)%Edge%Vert2%in_out.AND. &
!!Ausblenden                    Edges_Z(i,j,k)%Edge%Vert1%in_out==0)) THEN
!!Ausblenden                       eznr=eznr+1
!!Ausblenden                       Topo(i,j,k)%Topo%eznr=eznr
!!Ausblenden                       Edges_Z(i,j,k)%Edge%eg_nr=eznr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                IF(MAX(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)>0 &
!!Ausblenden                   .OR. &
!!Ausblenden                   (Vertices(i,j,k-1)%in_out==Vertices(i,j,k)%in_out.AND. &
!!Ausblenden                    Vertices(i,j,k)%in_out==0)) THEN
!!Ausblenden                       eznr=eznr+1
!!Ausblenden                       Topo(i,j,k)%Topo%eznr=eznr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        ecnr=eznr   ! 4.EdgeCut
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                IF(Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                  !Edge-Cut zuordnen!!!
!!Ausblenden                  Topo(i,j,k)%Topo%fceg_nr=Cell(i,j,k)%Cell%vc
!!Ausblenden                  !----------------------------------------------------------------------------
!!Ausblenden                  !Analyze geschnittene Face_XY-/Face_YZ-/Face_ZX-; ecnr zuordnen                 
!!Ausblenden                  !Liste Topo(i,j,k)%Topo%fceg_nrcut(1:Cell(i,j,k)%Cell%vc) zuordnen
!!Ausblenden                  !----------------------------------------------------------------------------
!!Ausblenden             ! in Bearbeitung
!!Ausblenden             !---------------
!!Ausblenden             !IF(i==25.and.j==51.and.k==1) THEN
!!Ausblenden             !  Write(*,*) "Analayze EdgeCut zuordnen Analyze geschnittene Face_XY-/Face_YZ-/Face_ZX"
!!Ausblenden             !END IF
!!Ausblenden                  DO iec=1,Cell(i,j,k)%Cell%vc-1
!!Ausblenden                    !ueber die VertCutListe die Vert-Nummern ermitteln
!!Ausblenden                    !im Output
!!Ausblenden                    !IF(Faces_XY(i,j,k-1)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1))  .OR. &
!!Ausblenden                       (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_XY(i,j,k-1)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_XY(i,j,k-1)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_XY(i,j,k-1)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_XY(i,j,k-1)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_XY(i,j,k-1)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k-1)%Topo%ecnrxy=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k-1)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_XY(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_XY(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_XY(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_XY(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_XY(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                           ecnr=ecnr+1
!!Ausblenden                           Faces_XY(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                           Topo(i,j,k)%Topo%ecnrxy=ecnr
!!Ausblenden                           Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                           Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_YZ(i-1,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) )THEN
!!Ausblenden                        Edg(1)=Faces_YZ(i-1,j,k)%Face%Edge1
!!Ausblenden                        Edg(2)=Faces_YZ(i-1,j,k)%Face%Edge2
!!Ausblenden                        Edg(3)=Faces_YZ(i-1,j,k)%Face%Edge3
!!Ausblenden                        Edg(4)=Faces_YZ(i-1,j,k)%Face%Edge4
!!Ausblenden                        s_ec=.TRUE.
!!Ausblenden                        DO li=1,4
!!Ausblenden                          IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                              Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                              Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                          END IF
!!Ausblenden                        END DO
!!Ausblenden                        IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                             ecnr=ecnr+1
!!Ausblenden                             Faces_YZ(i-1,j,k)%Face%egc_nr=ecnr
!!Ausblenden                             Topo(i-1,j,k)%Topo%ecnryz=ecnr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                        ELSE IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i-1,j,k)%Face%egc_nr
!!Ausblenden                        END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_YZ(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_YZ(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_YZ(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_YZ(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                              Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_YZ(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%ecnryz=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_ZX(i,j-1,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_ZX(i,j-1,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_ZX(i,j-1,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_ZX(i,j-1,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_ZX(i,j-1,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_ZX(i,j-1,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j-1,k)%Topo%ecnrzx=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j-1,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
!!Ausblenden                       (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
!!Ausblenden                       Edg(1)=Faces_ZX(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_ZX(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_ZX(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_ZX(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_ZX(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%ecnrzx=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                  END DO
!!Ausblenden                    !IF(Faces_XY(i,j,k-1)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
!!Ausblenden                       Edg(1)=Faces_XY(i,j,k-1)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_XY(i,j,k-1)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_XY(i,j,k-1)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_XY(i,j,k-1)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_XY(i,j,k-1)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k-1)%Topo%ecnrxy=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k-1)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_XY(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )    .OR. & 
!!Ausblenden                       (Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
!!Ausblenden                       Edg(1)=Faces_XY(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_XY(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_XY(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_XY(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE. 
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                           ecnr=ecnr+1
!!Ausblenden                           Faces_XY(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                           Topo(i,j,k)%Topo%ecnrxy=ecnr
!!Ausblenden                           Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                           Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_YZ(i-1,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
!!Ausblenden                        Edg(1)=Faces_YZ(i-1,j,k)%Face%Edge1
!!Ausblenden                        Edg(2)=Faces_YZ(i-1,j,k)%Face%Edge2
!!Ausblenden                        Edg(3)=Faces_YZ(i-1,j,k)%Face%Edge3
!!Ausblenden                        Edg(4)=Faces_YZ(i-1,j,k)%Face%Edge4
!!Ausblenden                        s_ec=.TRUE.
!!Ausblenden                        DO li=1,4
!!Ausblenden                          IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                              Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                              Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
!!Ausblenden                              Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                              s_ec=.FALSE. 
!!Ausblenden                          END IF
!!Ausblenden                        END DO
!!Ausblenden                        IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                             ecnr=ecnr+1
!!Ausblenden                             Faces_YZ(i-1,j,k)%Face%egc_nr=ecnr
!!Ausblenden                             Topo(i-1,j,k)%Topo%ecnryz=ecnr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                        ELSE IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i-1,j,k)%Face%egc_nr
!!Ausblenden                        END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
!!Ausblenden                       Edg(1)=Faces_YZ(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_YZ(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_YZ(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_YZ(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_YZ(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%ecnryz=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_ZX(i,j-1,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
!!Ausblenden                       Edg(1)=Faces_ZX(i,j-1,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_ZX(i,j-1,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_ZX(i,j-1,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_ZX(i,j-1,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_ZX(i,j-1,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j-1,k)%Topo%ecnrzx=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j-1,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                    !IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
!!Ausblenden                    IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
!!Ausblenden                       (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
!!Ausblenden                        Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
!!Ausblenden                       Edg(1)=Faces_ZX(i,j,k)%Face%Edge1
!!Ausblenden                       Edg(2)=Faces_ZX(i,j,k)%Face%Edge2
!!Ausblenden                       Edg(3)=Faces_ZX(i,j,k)%Face%Edge3
!!Ausblenden                       Edg(4)=Faces_ZX(i,j,k)%Face%Edge4
!!Ausblenden                       s_ec=.TRUE.
!!Ausblenden                       DO li=1,4
!!Ausblenden                         IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
!!Ausblenden                            (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                             Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                          ELSE IF(Edg(li)%yes_sp==1 .AND. &
!!Ausblenden                            ((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
!!Ausblenden                             (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
!!Ausblenden                              Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
!!Ausblenden                             Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
!!Ausblenden                             Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
!!Ausblenden                             s_ec=.FALSE.
!!Ausblenden                         END IF
!!Ausblenden                       END DO
!!Ausblenden                       IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr==0) THEN
!!Ausblenden                            ecnr=ecnr+1
!!Ausblenden                            Faces_ZX(i,j,k)%Face%egc_nr=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%ecnrzx=ecnr
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
!!Ausblenden                       ELSE IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr>0) THEN
!!Ausblenden                            Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j,k)%Face%egc_nr
!!Ausblenden                       END IF
!!Ausblenden                    END IF
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        topo_enr=ecnr
!!Ausblenden        topo_ecnr=topo_enr-eznr
!!Ausblenden!        Write(*,*) "Troposphäre:
!!Ausblenden!        Write(*,*) "topo_enr:         ", "enr-summe"
!!Ausblenden!        Write(*,*) "exnr =",exnr,      "  ",exnr
!!Ausblenden!        Write(*,*) "eynr =",eynr-exnr, "  ",eynr
!!Ausblenden!        Write(*,*) "eznr =",eznr-eynr, "  ",eznr
!!Ausblenden!        Write(*,*) "ecnr =",topo_enr-eznr, "  ", topo_enr
!!Ausblenden
!!Ausblenden        !.................................
!!Ausblenden        !--zu Orography /EdgeX/EdgeY/EdgeZ
!!Ausblenden        !.................................
!!Ausblenden        !                  ! 1.Edge_X
!!Ausblenden        exnr=topo_enr
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF(MIN(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
!!Ausblenden                  exnr=exnr+1
!!Ausblenden                  Topo(i,j,k)%Topo%iexnr=exnr
!!Ausblenden!W                !Write(*,*) "Orography: exnr= ",exnr, Topo(i,j,k)%Topo%exnr,"  "," i, j, k :", i,j,k 
!!Ausblenden              ELSE IF(MIN(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==0 .AND. &
!!Ausblenden                      (Vertices(i-1,j,k)%in_out==0.and.Vertices(i,j,k)%in_out==0) ) THEN
!!Ausblenden                     !Special Grenze Edge komplett)
!!Ausblenden                  Topo(i,j,k)%Topo%iexnr=Topo(i,j,k)%Topo%exnr
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        eynr=exnr   ! 2.Edge_Y
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF(MIN(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
!!Ausblenden                  eynr=eynr+1
!!Ausblenden                  Topo(i,j,k)%Topo%ieynr=eynr
!!Ausblenden!W                !Write(*,*) "Orography: eynr= ",eynr, Topo(i,j,k)%Topo%eynr,"  "," i, j, k :", i,j,k 
!!Ausblenden              ELSE IF(MIN(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)==0 .AND. &
!!Ausblenden                      (Vertices(i,j-1,k)%in_out==0.and.Vertices(i,j,k)%in_out==0) ) THEN
!!Ausblenden                     !Special Grenze Edge komplett)
!!Ausblenden                  Topo(i,j,k)%Topo%ieynr=Topo(i,j,k)%Topo%eynr
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        eznr=eynr   ! 3.Edge_Z
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF(MIN(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
!!Ausblenden                  eznr=eznr+1
!!Ausblenden                  Topo(i,j,k)%Topo%ieznr=eznr 
!!Ausblenden!W                !Write(*,*) "Orography: eznr= ",eznr, Topo(i,j,k)%Topo%eznr,"  "," i, j, k :", i,j,k 
!!Ausblenden              ELSE IF(MIN(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)==0 .AND. &
!!Ausblenden                      (Vertices(i,j,k-1)%in_out==0.and.Vertices(i,j,k)%in_out==0)) THEN
!!Ausblenden                     !Special Grenze Edge komplett)
!!Ausblenden                  Topo(i,j,k)%Topo%ieznr=Topo(i,j,k)%Topo%eznr 
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        para_grid_enr=eznr
!!Ausblenden        para_oro_enr=para_grid_enr-topo_enr
!!Ausblenden!        Write(*,*) "Orography_enr :",para_oro_enr
!!Ausblenden!        Write(*,*) "exnr =",topo_enr-exnr,      "  ",exnr
!!Ausblenden!        Write(*,*) "eynr =",eynr-exnr,          "  ",eynr
!!Ausblenden!        Write(*,*) "eznr =",eznr-eynr,          "  ",eznr
!!Ausblenden!        Write(*,*) "para_grid_enr =",para_grid_enr,"  "," grid_all_enr" 
!!Ausblenden!        Ende Orography-Edge
!!Ausblenden        !------------------------------------------------------------------------
!!Ausblenden        !FaceCut,Faces_XY,Faces_XZ,Faces_YZ
!!Ausblenden        fcnr=0
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden                IF(Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                    fcnr=fcnr+1
!!Ausblenden                    Topo(i,j,k)%Topo%FCutNr=fcnr
!!Ausblenden!W                  !Write(10,*) "Topo(i,j,k)%Topo%FCutNr =", fcnr, i,j,k
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden!W      ! Write(*,*) "fcnr=", fcnr, " --> anz. fcuts"
!!Ausblenden        fxynr=fcnr
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
!!Ausblenden                IF(Vertices(i-1,j-1,k)%in_out==0 .and. &
!!Ausblenden                   Vertices(i  ,j-1,k)%in_out==0 .and. &
!!Ausblenden                   Vertices(i-1,j  ,k)%in_out==0 .and. &
!!Ausblenden                   Vertices(i  ,j  ,k)%in_out==0) THEN
!!Ausblenden                   !Special-Grenzfläche als 'fc' ausgewertet,
!!Ausblenden                   !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
!!Ausblenden                   !Write(10,*) "Faces_XY(i,j,k)%Face Special Face Grenze als fc",  i,j,k
!!Ausblenden                   IF(Cell(i,j,k)%Cell%vc==0.and.Cell(i,j,k+1)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fxynr=Topo(i,j,k+1)%Topo%FCutNr
!!Ausblenden                      Topo(i,j,k)%Topo%FCutNr=Topo(i,j,k+1)%Topo%FCutNr
!!Ausblenden                   ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fxynr=Topo(i,j,k)%Topo%FCutNr
!!Ausblenden                   ELSE 
!!Ausblenden                      Write(OutUnitProt,*)"Warning aus SetTropoOroParaView(): "
!!Ausblenden                      Write(OutUnitProt,*)"  Face_XY(i,j,k)-Special-Neu gefunden!"  ,i,j,k
!!Ausblenden                      CALL WriteFaceProt2(Faces_XY(i,j,k)%Face)
!!Ausblenden                   END IF
!!Ausblenden                ELSE IF(Faces_XY(i,j,k)%Face%in_out>-4.and.Faces_XY(i,j,k)%Face%numberVert>2) THEN
!!Ausblenden                   fxynr=fxynr+1 
!!Ausblenden                   Topo(i,j,k)%Topo%fxynr=fxynr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                in_out=Vertices(i-1,j-1,k)%in_out+Vertices(i,j-1,k)%in_out+ &
!!Ausblenden                       Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
!!Ausblenden                IF(in_out>=0) THEN
!!Ausblenden                    fxynr=fxynr+1 
!!Ausblenden                    Topo(i,j,k)%Topo%fxynr=fxynr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        fzxnr=fxynr
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
!!Ausblenden                IF(Vertices(i-1,j,k-1)%in_out==0 .and. &
!!Ausblenden                   Vertices(i  ,j,k-1)%in_out==0 .and. & 
!!Ausblenden                   Vertices(i-1,j,k  )%in_out==0 .and. &
!!Ausblenden                   Vertices(i  ,j,k  )%in_out==0 ) THEN
!!Ausblenden                   !Special Grenzfläche als 'fc' ausgewertet,
!!Ausblenden                   !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
!!Ausblenden                   IF(Cell(i,j,k)%Cell%vc==0.and.Cell(i,j+1,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fzxnr=Topo(i,j+1,k)%Topo%FCutNr
!!Ausblenden                      Topo(i,j,k)%Topo%FCutNr=Topo(i,j+1,k)%Topo%FCutNr
!!Ausblenden                   ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fzxnr=Topo(i,j,k)%Topo%FCutNr
!!Ausblenden                   ELSE 
!!Ausblenden                      Write(OutUnitProt,*)"Warning aus SetTropoOroParaView(): "
!!Ausblenden                      Write(OutUnitProt,*)"  Face_ZX(i,j,k)-Special-Neu gefunden! "  ,i,j,k
!!Ausblenden                      CALL WriteFaceProt2(Faces_ZX(i,j,k)%Face)
!!Ausblenden                   END IF
!!Ausblenden                ELSE IF(Faces_ZX(i,j,k)%Face%in_out>-4.and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
!!Ausblenden                   fzxnr=fzxnr+1 
!!Ausblenden                   Topo(i,j,k)%Topo%fzxnr=fzxnr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                in_out=Vertices(i-1,j,k-1)%in_out+Vertices(i,j,k-1)%in_out+ &
!!Ausblenden                       Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
!!Ausblenden                IF(in_out>=0) THEN
!!Ausblenden                  fzxnr=fzxnr+1
!!Ausblenden                  Topo(i,j,k)%Topo%fzxnr=fzxnr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        fyznr=fzxnr
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
!!Ausblenden                IF(Vertices(i,j-1,k-1)%in_out==0 .and. &
!!Ausblenden                   Vertices(i,j-1,k  )%in_out==0 .and. &
!!Ausblenden                   Vertices(i,j  ,k  )%in_out==0 .and. &
!!Ausblenden                   Vertices(i,j  ,k-1)%in_out==0 ) THEN
!!Ausblenden                   !Special Grenzfläche als 'fc' ausgewertet,
!!Ausblenden                   !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
!!Ausblenden                   IF(Cell(i,j,k)%Cell%vc==0.and.Cell(i+1,j,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fyznr=Topo(i+1,j,k)%Topo%FCutNr
!!Ausblenden                      Topo(i,j,k)%Topo%FCutNr=Topo(i+1,j,k)%Topo%FCutNr
!!Ausblenden                   ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
!!Ausblenden                      Topo(i,j,k)%Topo%fyznr=Topo(i,j,k)%Topo%FCutNr
!!Ausblenden                   ELSE 
!!Ausblenden                      Write(OutUnitProt,*)"Warning aus SetTropoOroParaView(): "
!!Ausblenden                      Write(OutUnitProt,*)"  Face_YZ(i,j,k)-Special-Neu gefunden! "  ,i,j,k
!!Ausblenden                      CALL WriteFaceProt2(Faces_YZ(i,j,k)%Face)
!!Ausblenden                   END IF
!!Ausblenden                ELSE IF(Faces_YZ(i,j,k)%Face%in_out>-4.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
!!Ausblenden                  fyznr=fyznr+1 
!!Ausblenden                  Topo(i,j,k)%Topo%fyznr=fyznr
!!Ausblenden                END IF
!!Ausblenden              ELSE
!!Ausblenden                in_out=Vertices(i,j-1,k-1)%in_out+Vertices(i,j-1,k)%in_out+ &
!!Ausblenden                       Vertices(i,j,k)%in_out+Vertices(i,j,k-1)%in_out
!!Ausblenden                IF(in_out>=0) THEN
!!Ausblenden                  fyznr=fyznr+1 
!!Ausblenden                  Topo(i,j,k)%Topo%fyznr=fyznr
!!Ausblenden                END IF
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        topo_fnr=fyznr
!!Ausblenden
!!Ausblenden!!        Write(*,*) "topo_fnr:        ", "fnr-summe"
!!Ausblenden!!        Write(*,*) "fcnr  =",fcnr,        "  ",fcnr
!!Ausblenden!!        Write(*,*) "fxynr =",fxynr-fcnr,  "  ",fxynr
!!Ausblenden!!        Write(*,*) "fzxnr =",fzxnr-fxynr, "  ",fzxnr
!!Ausblenden!!        Write(*,*) "fyznr =",fyznr-fzxnr, "  ",fyznr
!!Ausblenden        !.............
!!Ausblenden        !Orography:
!!Ausblenden        !.............
!!Ausblenden        fxynr=topo_fnr
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0,iz1
!!Ausblenden              IF( MIN(Vertices(i-1,j-1,k)%in_out,Vertices(i,j-1,k)%in_out, &
!!Ausblenden                &  Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
!!Ausblenden                fxynr=fxynr+1 
!!Ausblenden                Topo(i,j,k)%Topo%ifxynr=fxynr
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        fzxnr=fxynr
!!Ausblenden        DO i=ix0+1,ix1
!!Ausblenden          DO j=iy0,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF( MIN(Vertices(i-1,j,k-1)%in_out,Vertices(i,j,k-1)%in_out, &
!!Ausblenden                &  Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
!!Ausblenden                 fzxnr=fzxnr+1 
!!Ausblenden                 Topo(i,j,k)%Topo%ifzxnr=fzxnr
!!Ausblenden               END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        fyznr=fzxnr
!!Ausblenden        DO i=ix0,ix1
!!Ausblenden          DO j=iy0+1,iy1
!!Ausblenden            DO k=iz0+1,iz1
!!Ausblenden              IF( MIN(Vertices(i,j-1,k-1)%in_out,Vertices(i,j,k-1)%in_out, &
!!Ausblenden                &  Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
!!Ausblenden                  fyznr=fyznr+1 
!!Ausblenden                  Topo(i,j,k)%Topo%ifyznr=fyznr
!!Ausblenden              END IF
!!Ausblenden            END DO
!!Ausblenden          END DO
!!Ausblenden        END DO
!!Ausblenden        para_grid_fnr=fyznr
!!Ausblenden        para_oro_fnr=para_grid_fnr-topo_fnr 
!!Ausblenden        ! End Orography Faces
!!Ausblenden!!        Write(*,*) "Oro_fnr:        ", "fnr-summe"
!!Ausblenden!!        Write(*,*) "fxynr =",fxynr-topo_fnr,  "  ",fxynr
!!Ausblenden!!        Write(*,*) "fzxnr =",fzxnr-fxynr, "  ",fzxnr
!!Ausblenden!!        Write(*,*) "fyznr =",fyznr-fzxnr, "  ",fyznr
!!Ausblenden!!        Write(*,*) "foronr=",para_oro_fnr, "   ",para_grid_fnr
!!Ausblenden      END IF   !IF(view_cell=='y')
!!Ausblenden    END DO  !ib
!!Ausblenden  ELSE
!!Ausblenden     Write(*,*) "z.Zt.  SetTropOroCellsParaV ","  nur wenn  '#OutputDomain'  gesetzt ist!"  
!!Ausblenden!    WRITE(OutUnitProt,'(a8,i8)') cells,nr_cells
!!Ausblenden!    DO ib=1,nb
!!Ausblenden!      Write(*,*) "                                     Block : ",ib,"\/",nb
!!Ausblenden!      CALL Set(Floor(ib))
!!Ausblenden!      DO k=iz0+1,iz1
!!Ausblenden!        DO j=iy0+1,iy1
!!Ausblenden!          DO i=ix0+1,ix1
!!Ausblenden!            !CALL WriteCellGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!Ausblenden!            IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!Ausblenden!              Write(OutUnitProt,*) "...................................."
!!Ausblenden!              Write(OutUnitProt,*) "Cell(",i,",",j,",",k,")%Cell% ->"
!!Ausblenden!              Write(OutUnitProt,*) "  Face1%mp=",Cell(i,j,k)%Cell%Face1%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face2%mp=",Cell(i,j,k)%Cell%Face2%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face3%mp=",Cell(i,j,k)%Cell%Face3%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face4%mp=",Cell(i,j,k)%Cell%Face4%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face5%mp=",Cell(i,j,k)%Cell%Face5%mp
!!Ausblenden!              Write(OutUnitProt,*) "  Face6%mp=",Cell(i,j,k)%Cell%Face6%mp
!!Ausblenden!              Write(OutUnitProt,*) "  vc=",Cell(i,j,k)%Cell%vc
!!Ausblenden!              WRITE(OutUnitProt,*) (Cell(i,j,k)%Cell%VertCut(ic),ic=1,Cell(i,j,k)%Cell%vc)
!!Ausblenden!              Write(OutUnitProt,*) "...................................."
!!Ausblenden!            END IF
!!Ausblenden!          END DO
!!Ausblenden!        END DO
!!Ausblenden!      END DO
!!Ausblenden!    END DO  !ib
!!Ausblenden  END IF
!!Ausblenden
!!AusblendenEND SUBROUTINE SetTropoOroCellsParaV


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
