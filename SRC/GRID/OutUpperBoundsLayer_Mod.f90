MODULE OutUpperBoundsLayer_Mod

  USE F_Mod
  USE Floor_Mod
  USE IOControl_Mod
  USE Parametric_Mod
  USE OutputOutGmvG_Mod, ONLY: Search_OutArea

  IMPLICIT NONE

CONTAINS

SUBROUTINE OpenFileForBound(FName)
  CHARACTER*50 :: FName
  OPEN(UNIT=OutUnitBpv,FILE=TRIM(FName)//'.pva.bound',STATUS='UNKNOWN')
  OPEN(UNIT=OutUnitB,FILE=TRIM(FName)//'.bound',STATUS='UNKNOWN')
  WRITE(OutUnitBpv,*) "Bound  ",TRIM(FName)
  WRITE(OutUnitB,*) "Bound  ",TRIM(FName)
END SUBROUTINE OpenFileForBound

SUBROUTINE EndeOutBound(FName,blkn)
  CHARACTER*50 :: FName
  INTEGER      :: blkn
  Write(OutUnitB,*)   "Ende Bound ----- ",TRIM(FName), " -----  ",blkn,"Blöcke"
  Write(OutUnitBpv,*) "Ende Bound ----- ",TRIM(FName), " -----  ",blkn,"Blöcke"
END SUBROUTINE EndeOutBound

SUBROUTINE CloseFileForBound
  CLOSE(OutUnitB)
  CLOSE(OutUnitBpv)
END SUBROUTINE CloseFileForBound

SUBROUTINE WOutNrBlkGrid(gb)
  INTEGER :: gb
  Write(OutUnitB,*) gb, "   ! Number of blocks"
  Write(OutUnitBpv,*) gb, "   ! Number of blocks"
END SUBROUTINE WOutNrBlkGrid

SUBROUTINE WOutVertBounds()
  INTEGER :: i
  INTEGER :: nrb
  nrb=b_upbounds-znr_upbounds
  IF(b_upbounds>0) THEN
    WRITE(OutUnitB,*) "NrB_Nodes",b_upbounds,leerzei2,"NrB_n",nrp_upb_gi, &
                 &    leerzei2,"NrB_nc",nrp_upb_nc,leerzei2,"NrB_nb",nrb
  ELSE   
    WRITE(OutUnitB,*) "NrB_Nodes",znr_upbounds, leerzei2,"NrB_n",nrp_upb_gi, &
                 &    leerzei2,"NrB_nc",nrp_upb_nc,"Nr_nb", 0
  END IF
  DO i=1,nrp_upb_gi
     WRITE(OutUnitB,*)"n ", i, leerzei2, & 
            &   xParametricOut(VertOutBounds(i)), &
            &   yParametricOut(VertOutBounds(i)), &
            &   zParametricOut(VertOutBounds(i)), &
            &   IndexOutBounds(i,2:4)
  END DO
  DO i=nrp_upb_gi+1,znr_upbounds
     WRITE(OutUnitB,*)"nc", i, leerzei2, & 
            &   xParametricOut(VertOutBounds(i)), &
            &   yParametricOut(VertOutBounds(i)), &
            &   zParametricOut(VertOutBounds(i)), &
            &   IndexOutBounds(i,2:4)
  END DO
  IF(b_upbounds>0) THEN
  DO i=znr_upbounds+1,b_upbounds
     WRITE(OutUnitB,*)"nb", i,leerzei2, & 
            &   xParametricOut(VertOutBoundsBord(i)), &
            &   yParametricOut(VertOutBoundsBord(i)), &
            &   zParametricOut(VertOutBoundsBord(i)), &
            &   IndexBordBounds(i,2:4)
  END DO
  END IF
  !..................  
  IF(b_upbounds>0) THEN
    WRITE(OutUnitBpv,*) "NrB_Nodes",b_upbounds,leerzei2,"NrB_n",nrp_upb_gi, &
                   &    leerzei2,"NrB_nc",nrp_upb_nc,leerzei2,"NrB_nb",nrb
  ELSE   
    WRITE(OutUnitBpv,*) "NrB_Nodes",znr_upbounds, leerzei2,"NrB_n",nrp_upb_gi, &
                   &    leerzei2,"NrB_nc",nrp_upb_nc,"Nr_nb", 0
  END IF
  DO i=1,nrp_upb_gi
     WRITE(OutUnitBpv,*)"n ", i, leerzei2, & 
            &   xParametricOut(VertOutBounds(i)), &
            &   yParametricOut(VertOutBounds(i)), &
            &   zParametricOut(VertOutBounds(i))
  END DO
  DO i=nrp_upb_gi+1,znr_upbounds
     WRITE(OutUnitBpv,*)"nc", i, leerzei2, & 
            &   xParametricOut(VertOutBounds(i)), &
            &   yParametricOut(VertOutBounds(i)), &
            &   zParametricOut(VertOutBounds(i))
  END DO
  IF(b_upbounds>0) THEN
  DO i=znr_upbounds+1,b_upbounds
     WRITE(OutUnitBpv,*)"nb", i,leerzei2, & 
            &   xParametricOut(VertOutBoundsBord(i)), &
            &   yParametricOut(VertOutBoundsBord(i)), &
            &   zParametricOut(VertOutBoundsBord(i))
  END DO
  END IF
END SUBROUTINE WOutVertBounds

SUBROUTINE WOutNrBEdge(ges_ibenr,ib_bexnr,ib_beynr,ib_beznr,ib_becnr)
  INTEGER :: ges_ibenr,ib_bexnr,ib_beynr,ib_beznr,ib_becnr

  !Reihenfolge Output:bexnr,beynr,beznr,becnr
  Write(OutUnitB,*)"NrB_Edges",ges_ibenr, &
           &  leerzei2,"NrB_ex",ib_bexnr,leerzei2,"NrB_ey",ib_beynr, &
           &  leerzei2,"NrB_ez",ib_beznr,leerzei2,"NrB_ec",ib_becnr
  !..................
  Write(OutUnitBpv,*)"Nr_Edges",ges_ibenr
END SUBROUTINE WOutNrBEdge

SUBROUTINE WOutBEdge(egi,is,js,ks,te,benr,bVert1,bVert2)
  INTEGER      :: egi,is,js,ks 
  CHARACTER(2) :: te
  INTEGER      :: benr,bVert1,bVert2

  !Write(OutUnitB,*) te,benr,leerzei2,bVert1,bVert2,"   ",is,js,ks, " ",egi 
  Write(OutUnitB,*) leerzei6,te,benr,leerzei2,bVert1,bVert2,leerzei6,egi 
  Write(OutUnitB,*) is,js,ks 
  !..................
  Write(OutUnitBpv,*)        te,benr,leerzei2,bVert1,bVert2
END SUBROUTINE WOutBEdge

SUBROUTINE WOutBEdgeCut(pos_beg_ecnr,max_upb_enr,bl_ec)
  INTEGER :: pos_beg_ecnr,max_upb_enr
  INTEGER :: bl_ec(pos_beg_ecnr:max_upb_enr,1:6)
  ! local...
  INTEGER :: i,nl
  !-- Output--EdgeCut
  nl=0
  DO i=pos_beg_ecnr,max_upb_enr
     nl=nl+1
     Write(OutUnitB,*) leerzei6,"ec",bl_ec(i,1),leerzei2,bl_ec(i,2),bl_ec(i,3),leerzei6,nl
     Write(OutUnitB,*) bl_ec(i,4),bl_ec(i,5),bl_ec(i,6)
  END DO
  !..................
  DO i=pos_beg_ecnr,max_upb_enr
     Write(OutUnitBpv,*) "ec",bl_ec(i,1),leerzei2,bl_ec(i,2),bl_ec(i,3)
  END DO
END SUBROUTINE WOutBEdgeCut

SUBROUTINE WOutNrBFace(ges_ibfnr,ib_bfcnr,ib_bfxynr,ib_bfzxnr,ib_bfyznr)
  INTEGER :: ges_ibfnr,ib_bfcnr,ib_bfxynr,ib_bfzxnr,ib_bfyznr
  !Reihenfolge Zählung/Output: fc,fxy,fzx,fyz
  Write(OutUnitB,*) "NrB_Faces", ges_ibfnr, &
           &  leerzei2,"NrB_fc",ib_bfcnr,leerzei2,"NrB_fxy",ib_bfxynr, &
           &  leerzei2,"NrB_fzx",ib_bfzxnr,leerzei2,"NrB_fyz",ib_bfyznr
  !..................
  Write(OutUnitBpv,*) "NrB_Faces", ges_ibfnr
END SUBROUTINE WOutNrBFace

SUBROUTINE WOutBFace(nfa,is,js,ks,tf,bfnr,be_nr,tid,feg_l)
  INTEGER      :: nfa,is,js,ks
  CHARACTER(3) :: tf
  INTEGER      :: bfnr,be_nr,tid
  INTEGER      :: feg_l(1:8)
  ! Local
  INTEGER :: ie
  Write(OutUnitB,*)   leerzei6,tf,bfnr,leerzei2,be_nr,leerzei2,tid,leerzei6,nfa
  Write(OutUnitB,*)   leerzei14,(feg_l(ie),ie=1,be_nr)
  Write(OutUnitB,*)   is,js,ks 
  !..................
  Write(OutUnitBpv,*)          tf,bfnr,leerzei2,be_nr,leerzei2,tid
  Write(OutUnitBpv,*) leerzei8,(feg_l(ie),ie=1,be_nr)
END SUBROUTINE WOutBFace

SUBROUTINE WOutNrBCell(ges_ibcnr,ib_bcvcnr,ib_bcnr,ib_bcnr_bdr)
  INTEGER :: ges_ibcnr,ib_bcvcnr,ib_bcnr,ib_bcnr_bdr(1:6)
  !Reihenfolge Zählung/Output: C-vc,C-hex
  Write(OutUnitB,*) "NrB_Cells", ges_ibcnr, &
           &  leerzei2,"NrB_cvc",ib_bcvcnr,leerzei2,"NrB_ch",ib_bcnr, leerzei2, &
           &  "NrB_w-e-s-n-b-t",leerzei2, ib_bcnr_bdr
  !..................
  Write(OutUnitBpv,*) "NrB_Cells", ges_ibcnr
END SUBROUTINE WOutNrBCell

SUBROUTINE WOutBCell(ncg,is,js,ks,tc,bcnr,bf_nr,tid,cf_l)
  INTEGER      :: ncg,is,js,ks
  CHARACTER(2) :: tc
  INTEGER      :: bcnr,bf_nr,tid
  INTEGER      :: cf_l(1:7)
  ! loacal
  INTEGER :: n
    Write(OutUnitB,*) leerzei6, tc,bcnr,leerzei2,bf_nr,leerzei2,tid,leerzei6,ncg 
    Write(OutUnitB,*) leerzei14,(cf_l(n),n=1,bf_nr) 
    Write(OutUnitB,*) is,js,ks
    !...................
    Write(OutUnitBpv,*)         tc,bcnr,leerzei2,bf_nr,leerzei2,tid 
    Write(OutUnitBpv,*) leerzei8,(cf_l(n),n=1,bf_nr) 
END SUBROUTINE WOutBCell

SUBROUTINE ProtCountNrPBoundLayerIdentEdgeBorders(ssx,ssy,ssz,itest)
  INTEGER:: ssx,ssy,ssz,itest
  Write(*,*) "----------------------------------------------"
  Write(*,*) "Ergebnis: CountEdgeBoundLayerIdentEdgeBorders()"
  Write(*,*) "sum nb sbexnr            =  ", ssx
  Write(*,*) "sum nb sbeynr            =  ", ssy
  Write(*,*) "sum nb sbeznr            =  ", ssz
  Write(*,*) "sum nb sbexnr+nr_EgVerts =  ", ssx+nr_EgVerts
  Write(*,*) "----------------------------------------------"
  Write(*,*) "all sbe*nr x,y,z         =  ", itest
  Write(*,*) "gleiche_benr i-border    = -", gl_benr
  Write(*,*) "Ergebnis:                =  ", itest-gl_benr
  Write(*,*) "----------------------------------------------"
  Write(*,*) "nr_upbounds (nrP-Gitter) =  ", nr_upbounds
  Write(*,*) "nr_EgVerts               = +", nr_EgVerts   !(aus Anlayze Verts)
  Write(*,*) "nr_upbounds+nr_EgVerts   = ",  nr_upbounds+nr_EgVerts
  Write(*,*) "----------------------------------------------"
END SUBROUTINE ProtCountNrPBoundLayerIdentEdgeBorders
  

SUBROUTINE CountNrPBoundLayerIdentEdgeBorders
  INTEGER :: ib,in,i,j,k,ix,jx,iy,biz,mx
  INTEGER :: fac_ra,fac_ne,fac_n
  INTEGER :: nrP,ix_next,itest,ssx,ssy,ssz
  TYPE (Domain_T) :: flibn
  !...........................................................................
  ! -count nr_upbounds (nrP-Grid) of UpperBoundLayer 
  ! -count sbenr (number bound edge) identical benr of inside border blocks
  ! -count sum sbexnr, sbeynr,sbeznr of all blocks wioth border
  !...........................................................................
  nr_upbounds=0; out_bounds=0; per_bounds=0
  ssx=0;ssy=0;ssz=0
  nr_upbounds=0; itest=0
  !...
  DO ib=1,nb
    !Write(*,*) "                                     Block : ",ib,"\/",nb
    CALL Set(Floor(ib))
  
    DO ix=ix0,ix1
      ix_next=MIN(ix+1,ix1)
      DO iy=iy0,iy1
        DO biz=DimUpBoundsArray(ix,iy,1),iz1
          IF((biz>=DimUpBoundsArray(ix,iy,1).AND.biz<=DimUpBoundsArray(ix,iy,2)) .OR. &
             (UpBoundsLayer(ix     ,iy,biz)%BCell%bexnr>0 .OR. &
              UpBoundsLayer(ix_next,iy,biz)%BCell%bexnr>0) )THEN
              ! --> schauen Debuger Inhalte, ob nrP Abfrage Grenze ausreicht
              ! zählt i-direction point doppelt
            nrP=Vertices(ix,iy,biz)%nrP
            IF(nrP>0) THEN 
              nr_upbounds=nr_upbounds+1
            END IF
          END IF
        END DO
      END DO
    END DO

    DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !-------------------------------------------------------------------------
        CASE ("w","e")
        !-------------
          SELECT CASE (Nachbars(in)%nType(1:1))
            !-----------------------TypeW=="iw".OR.TypeE=="ie"---------------------
            ! count out_bounds the direction ow,oe of border upper bounds layer 
            CASE ('o')
               IF(nType=="ow") THEN
                  ix=ix0+1               !NachbarPoint_Eigener
                  fac_n=0
               ELSE !'oe'
                  ix=ix1-1               !NachbarPoint_Eigener
                  fac_n=1
               END IF
               !nrP-Grid
               DO iy=iy0,iy1
                 DO biz=DimUpBoundsArray(ix,iy,1),DimUpBoundsArray(ix,iy,2)
                   !IF( UpBoundsLayer(ix,iy,biz)%BCell%bexnr>0) THEN
                       ! --> schauen Debuger Inhalte, ob nrP Abfrage Grenze ausreicht
                       ! zählt i-direction point doppelt
                     IF(Vertices(ix,iy,biz)%nrP>0) THEN
                        out_bounds=out_bounds+1
                     END IF
                   !END IF
                 END DO
               END DO
               !VertS-EdgeX
               DO iy=iy0,iy1
                 DO biz=DimUpBoundsArray(ix+fac_n,iy,1),DimUpBoundsArray(ix+fac_n,iy,2)
                   !IF( UpBoundsLayer(ix,iy,biz)%BCell%bexnr>0) THEN
                      IF (ASSOCIATED(Edges_X(ix+fac_n,iy,biz)%Edge)) THEN
                        IF (Edges_X(ix+fac_n,iy,biz)%Edge%yes_sp==1) THEN
                           out_bounds=out_bounds+1
                        END IF
                      END IF
                   !END IF
                 END DO
               END DO
               !VertS-EdgeY
               DO iy=iy0+1,iy1
                 DO biz=DimUpBoundsArray(ix,iy,1),DimUpBoundsArray(ix,iy,2)
                   !IF( UpBoundsLayer(ix,iy,biz)%BCell%beynr>0) THEN
                     IF (ASSOCIATED(Edges_Y(ix,iy,biz)%Edge)) THEN
                         IF (Edges_Y(ix,iy,biz)%Edge%yes_sp==1) THEN
                           out_bounds=out_bounds+1
                        END IF
                      END IF
                   !END IF
                 END DO
               END DO
               !VertS-EdgeZ
               DO iy=iy0,iy1
                 DO biz=DimUpBoundsArray(ix,iy,1)+1,DimUpBoundsArray(ix,iy,2)
                   !IF( UpBoundsLayer(ix,iy,biz)%BCell%beznr>0) THEN
                      IF (ASSOCIATED(Edges_Z(ix,iy,biz)%Edge)) THEN
                         IF (Edges_Z(ix,iy,biz)%Edge%yes_sp==1) THEN
                           out_bounds=out_bounds+1
                        END IF
                      END IF
                   !END IF
                 END DO
               END DO
            !-----------------------TypeW=="pw".AND.TypeE=="pe"---------------------
            ! count out_bounds the direction ow,oe of border upper bounds layer 
             IF(nType=="pw") THEN
                  jx=Floor(ibn)%ix1-1    !NachbarPoint
                  fac_n =1
               ELSE !'pe'
                  jx=Floor(ibn)%ix0+1    !NachbarPoint
                  fac_n =0
               END IF
               flibn=Floor(ibn)
               !nrP-Grid
               DO iy=iy0,iy1
                 DO biz=flibn%DimUpBoundsArray(jx,iy,1),flibn%DimUpBoundsArray(jx,iy,2)
                   !IF( UpBoundsLayer(ix,iy,biz)%BCell%bexnr>0) THEN
                       ! --> schauen Debuger Inhalte, ob nrP Abfrage Grenze ausreicht
                       ! zählt i-direction point doppelt
                     IF(flibn%Vertices(jx,iy,biz)%nrP>0) THEN
                        per_bounds=per_bounds+1
                     END IF
                   !END IF
                 END DO
               END DO
               !VertS-EdgeX
               mx=jx+fac_n
               DO iy=iy0,iy1
                 DO biz=flibn%DimUpBoundsArray(mx,iy,1),flibn%DimUpBoundsArray(mx,iy,2)
                   !IF( UpBoundsLayer(mx,iy,biz)%BCell%bexnr>0) THEN
                      IF (ASSOCIATED(flibn%Edges_X(mx,iy,biz)%Edge)) THEN
                        IF (flibn%Edges_X(mx,iy,biz)%Edge%yes_sp==1) THEN
                           per_bounds=per_bounds+1
                        END IF
                      END IF
                   !END IF
                 END DO
               END DO
               !VertS-EdgeY
               DO iy=iy0+1,iy1
                 DO biz=flibn%DimUpBoundsArray(jx,iy,1),flibn%DimUpBoundsArray(jx,iy,2)
                   !IF( UpBoundsLayer(jx,iy,biz)%BCell%beynr>0) THEN
                     IF (ASSOCIATED(flibn%Edges_Y(jx,iy,biz)%Edge)) THEN
                         IF (flibn%Edges_Y(jx,iy,biz)%Edge%yes_sp==1) THEN
                           per_bounds=per_bounds+1
                        END IF
                      END IF
                   !END IF
                 END DO
               END DO
               !VertS-EdgeZ
               DO iy=iy0,iy1
                 DO biz=flibn%DimUpBoundsArray(jx,iy,1)+1,flibn%DimUpBoundsArray(jx,iy,2)
                   !IF( UpBoundsLayer(jx,iy,biz)%BCell%beznr>0) THEN
                      IF (ASSOCIATED(flibn%Edges_Z(jx,iy,biz)%Edge)) THEN
                         IF (flibn%Edges_Z(jx,iy,biz)%Edge%yes_sp==1) THEN
                           per_bounds=per_bounds+1
                        END IF
                      END IF
                   !END IF
                 END DO
               END DO

            !-----------------------TypeW=="iw".OR.TypeE=="ie"---------------------
            ! 'i' -> west,east  z.Zt: nur Refine =0, Incr->X,Y,Z noch offen  
            CASE ('i')
               IF(nType=="iw") THEN
                  ix=ix0               !aktuellRandCelle
                  jx=Floor(ibn)%ix1    !NachbarCelle
                  fac_ra=1
                  fac_n =1
               ELSE !'ie'
                  ix=ix1+1             !aktuellRandCelle
                  jx=Floor(ibn)%ix0+1  !NachbarCelle
                  fac_ra=0
                  fac_n =0
               END IF
  
               DO j=iy0,iy1            !EdgeX   'i-we'
                 DO k=iz0,iz1
                   IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bexnr>0) THEN
                     IF(UpBoundsLayer(ix,j,k)%BCell%bexnr== &
                       & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bexnr) THEN
                         gl_benr=gl_benr+1
                     END IF
                   END IF
                 END DO
               END DO
               DO j=iy0+1,iy1          !EdgeY  'i_we'
                 DO k=iz0,iz1
                   IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beynr>0) THEN
                     IF(UpBoundsLayer(ix-fac_ra,j,k)%BCell%beynr== &
                       & Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beynr ) THEN
                         gl_benr=gl_benr+1
                     END IF
                   END IF
                 END DO
               END DO
               DO j=iy0,iy1            !EdgeZ   'i-we'
                 DO k=iz0+1,iz1
                   IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beznr>0) THEN
                     IF(UpBoundsLayer(ix-fac_ra,j,k)%BCell%beznr== &
                      & Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beznr) THEN
                         gl_benr=gl_benr+1
                     END IF
                   END IF
                 END DO
               END DO
          END SELECT   ! [i]-> 'w/e'
  
        !--------------------------------------------------------------------
        CASE ("n","s")
        !-------------
          SELECT CASE (Nachbars(in)%nType(1:1))
            ! 'i' -> north,south
            CASE ('i')
              IF (nType=="in") THEN
              ELSE IF (nType=="is") THEN
              END IF 
          END SELECT   ! [i]-> 'n/s'
  
        !--------------------------------------------------------------------
        CASE ("t","b")
        !-------------
          SELECT CASE (Nachbars(in)%nType(1:1))
            ! 'i' -> top, bottom
            CASE ('i')
               IF (nType=="it") THEN
               ELSE IF (nType=="ib") THEN
               END IF 
          END SELECT   ! [i]-> 't,b'
  
      END SELECT   ! 'w/e','n/s','t/b' -> to boundslayer
    END DO  ! AnzahlNachbar
    itest=itest+Floor(ib)%sbexnr+Floor(ib)%sbeynr+Floor(ib)%sbeznr
    ssx=ssx+Floor(ib)%sbexnr
    ssy=ssy+Floor(ib)%sbeynr
    ssz=ssz+Floor(ib)%sbeznr
  END DO  !ib
  !CALL ProtCountNrPBoundLayerIdentEdgeBorders(ssx,ssy,ssz,itest)

END SUBROUTINE CountNrPBoundLayerIdentEdgeBorders
  
SUBROUTINE Prot_SortIndexOutUpBounds
   Write(OutUnitProt,*) " "
   Write(OutUnitProt,*) leerzei3,"> SUBROUTINE SortIndexOutUpBounds()"
   Write(OutUnitProt,*) leerzei3,"  Allokierte-Aktuelle Anzahl zu IndexOutBounds(:)"
   Write(OutUnitProt,*) leerzei3,"--------------------------------------------------"
   IF(nr_upbounds>znr_upbounds) THEN
       Write(OutUnitProt,*) leerzei5, "i.o.  ", &
                      &     "nr_upbounds/=znr_upbounds", nr_upbounds,"/=",znr_upbounds
       Write(OutUnitProt,*) " "
   ELSE IF (nr_upbounds==znr_upbounds) THEN
       Write(OutUnitProt,*) leerzei5, "i.o.  ", &
                      &     "nr_upbounds=znr_upbounds", nr_upbounds,"=",znr_upbounds
       Write(OutUnitProt,*) " "
   ELSE
       Write(OutUnitProt,*) leerzei5,"Fehler!" &
                      &    ,"nr_upbounds<znr_upbounds", nr_upbounds,"<",znr_upbounds
       Write(OutUnitProt,*) " "
   END IF
END SUBROUTINE Prot_SortIndexOutUpBounds

   SUBROUTINE SortIndexOutUpBounds
     INTEGER :: ib,ix,iy,gz,biz,i,j,k,jx
     INTEGER :: nrP,ix_next
     REAL(8) :: mvx
     INTEGER :: in
     INTEGER :: itest,fac_ra,fac_n
     !.................................................................
     ! SortIndexOutUpBounds: Sort all Points in/of Upper Bounds Layer 
     !  1) Grid-Vertex-Points with in_out=0//1 in/of Upper Bounds Layer
     !  2) Nr. Verts/Points -> of Edge_X,_Y,_Z in/of upper border layer
     !.................................................................

     CALL CountNrPBoundLayerIdentEdgeBorders

     nr_upbounds=nr_upbounds+nr_EgVerts
     ALLOCATE(IndexOutBounds(1:nr_upbounds,1:4))
     ALLOCATE(VertOutBounds(1:nr_upbounds))
     ALLOCATE(NrPOutBounds(1:nr_out))
     IndexOutBounds(:,:)=0
     NrPOutBounds(:)=0
     znr_upbounds=0
     !Write(*,*) "per_bounds = ",per_bounds, &
     !       &   "  1.Analyze aus CountNrPBoundLayerIdentEdgeBorders Nachbar-Varainte"
     per_bounds=0
   
     DO ib=1,nb
   
       CALL Set(Floor(ib))
   
       !1) all Grid-Vertex-Points with in_out==0//1 the in/of upper border layer
       DO ix=ix0,ix1
         ix_next=MIN(ix+1,ix1)
         DO iy=iy0,iy1
           DO biz=DimUpBoundsArray(ix,iy,1),iz1
           !gz=MIN(DimUpBoundsArray(ix,iy,2)+1,iz1) !zu Special in_out=6,vc=0 wenn allok.
           !DO biz=DimUpBoundsArray(ix,iy,1),gz
             IF((biz>=DimUpBoundsArray(ix,iy,1).AND.biz<=DimUpBoundsArray(ix,iy,2)) .OR. &
                (UpBoundsLayer(ix  ,iy,biz)%BCell%bexnr>0 .OR. &
                 UpBoundsLayer(ix_next,iy,biz)%BCell%bexnr>0) ) THEN
                 ! oder ? check Debugger nach Vertices()%nrP, special Schnitte beachten!!
               nrP=Vertices(ix,iy,biz)%nrP
               IF(nrP>0) THEN 
                 IF( UpBoundsLayer(ix  ,iy,biz)%BCell%nrP==0) THEN
                   znr_upbounds=znr_upbounds+1
                   UpBoundsLayer(ix  ,iy,biz)%BCell%nrP=znr_upbounds
                   VertOutBounds(znr_upbounds)=Vertices(ix,iy,biz)
                   NrPOutBounds(Vertices(ix,iy,biz)%nrP)=znr_upbounds
                   IndexOutBounds(znr_upbounds,1)=Vertices(ix,iy,biz)%nrP
                   IndexOutBounds(znr_upbounds,2)=ix
                   IndexOutBounds(znr_upbounds,3)=iy
                   IndexOutBounds(znr_upbounds,4)=biz
                   If(ix==ix0+1.OR.ix==ix1-1) THEN
                     per_bounds=per_bounds+1
                   END IF
                 END IF
               END IF
             END IF
           END DO
         END DO
       END DO
       nrp_upb_gi=znr_upbounds
       !2) Nr. Verts/Points -> of Edge_X,_Y,_Z in/of upper border layer   
       !           ""       -> analogically with order of VertOut() 
       !                       and in the end the list
       DO i=ix0+1,ix1
         DO j=iy0,iy1
           DO k=iz0,iz1
             IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
               IF (Edges_X(i,j,k)%Edge%yes_sp==1) THEN
                 nrP=Edges_X(i,j,k)%Edge%VertS%nrP
                 IF(nrP>0) THEN
                   IF(UpBoundsLayer(i,j,k)%BCell%nrEcx==0) THEN
                      znr_upbounds=znr_upbounds+1
                      UpBoundsLayer(i,j,k)%BCell%nrEcx=znr_upbounds
                      !....
                      VertOutBounds(znr_upbounds)=Edges_X(i,j,k)%Edge%VertS
                      NrPOutBounds(Edges_X(i,j,k)%Edge%VertS%nrP)=znr_upbounds
                      !....
                      IndexOutBounds(znr_upbounds,1)=Edges_X(i,j,k)%Edge%VertS%nrP
                      IndexOutBounds(znr_upbounds,2)=i
                      IndexOutBounds(znr_upbounds,3)=j
                      IndexOutBounds(znr_upbounds,4)=k
                      If(i==ix0+1.OR.i==ix1) THEN
                        per_bounds=per_bounds+1
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
             IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
               IF (Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
                 nrP=Edges_Y(i,j,k)%Edge%VertS%nrP
                 IF(nrP>0) THEN
                   IF(UpBoundsLayer(i,j,k)%BCell%nrEcy==0) THEN
                      znr_upbounds=znr_upbounds+1
                      UpBoundsLayer(i,j,k)%BCell%nrEcy=znr_upbounds
                      !...
                      VertOutBounds(znr_upbounds)=Edges_Y(i,j,k)%Edge%VertS
                      NrPOutBounds(Edges_Y(i,j,k)%Edge%VertS%nrP)=znr_upbounds
                      !...
                      IndexOutBounds(znr_upbounds,1)=Edges_Y(i,j,k)%Edge%VertS%nrP
                      IndexOutBounds(znr_upbounds,2)=i
                      IndexOutBounds(znr_upbounds,3)=j
                      IndexOutBounds(znr_upbounds,4)=k
                      If(i==ix0+1.OR.i==ix1-1) THEN
                        per_bounds=per_bounds+1
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
             IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
               IF (Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
                 nrP=Edges_Z(i,j,k)%Edge%VertS%NrP
                 IF(nrP>0) THEN
                   IF(UpBoundsLayer(i,j,k)%BCell%nrEcz==0) THEN
                      znr_upbounds=znr_upbounds+1
                      UpBoundsLayer(i,j,k)%BCell%nrEcz=znr_upbounds
                      !....
                      VertOutBounds(znr_upbounds)=Edges_Z(i,j,k)%Edge%VertS
                      NrPOutBounds(Edges_Z(i,j,k)%Edge%VertS%NrP)=znr_upbounds
                      !...
                      IndexOutBounds(znr_upbounds,1)=Edges_Z(i,j,k)%Edge%VertS%NrP
                      IndexOutBounds(znr_upbounds,2)=i
                      IndexOutBounds(znr_upbounds,3)=j
                      IndexOutBounds(znr_upbounds,4)=k
                      If(i==ix0+1.OR.i==ix1-1) THEN
                        per_bounds=per_bounds+1
                      END IF
                   END IF
                 END IF
               END IF
             END IF
           END DO
         END DO
       END DO
       nrp_upb_nc=znr_upbounds-nrp_upb_gi
     END DO  !ib

     CALL Prot_SortIndexOutUpBounds
   
    ! Write(*,*) "out_bounds border new = ",out_bounds
    ! Write(*,*) "per_bounds border  = ",per_bounds
    ! Write(*,*) "znr_upbounds       = ",znr_upbounds

     !Ränder.......................................
     IF(per_bounds>0) THEN
     !ALLOCATE(VertOutBoundsBord(znr_upbounds+1:znr_upbounds+per_bounds))
     !ALLOCATE(IndexBordBounds(znr_upbounds+1:znr_upbounds+per_bounds,1:4))
     ALLOCATE(VertOutBoundsBord(znr_upbounds+1:200))
     ALLOCATE(IndexBordBounds(znr_upbounds+1:200,1:4))
     b_upbounds=znr_upbounds 

     DO ib=1,nb
       CALL Set(Floor(ib))
       DO in=1,AnzahlNachbar
         CALL Set(Nachbars(in))
         SELECT CASE (Nachbars(in)%nType(2:2))
           !----------------------------------
         CASE ("w","e")
             SELECT CASE (Nachbars(in)%nType(1:1)) 
             CASE ('o')
          
                IF (nType=="ow") THEN
                   !!! noch Anzahl vorrauszählen organisieren für ausen-Points
                   !spiegle Celle-Points ix0, damit gleiche Volume Celle-ix0-1    
                   ix=ix0-1   !aktuellRandPoint
                   jx=ix0+1   !grid-'owest'-Point, spiegle Points der Celle ix0+1
                   mvx=dx(ix0)+dx(ix0+1)
                   DO j=iy0,iy1
                     DO k=iz0,iz1
                       IF(UpBoundsLayer(jx,j,k)%BCell%nrP>0) THEN        
                         nrP=UpBoundsLayer(jx,j,k)%BCell%nrP
                         b_upbounds=b_upbounds+1
                         UpBoundsLayer(ix,j,k)%BCell%nrP=b_upbounds
                         VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                         VertOutBoundsBord(b_upbounds)%Point%x= &
                              & VertOutBounds(nrP)%Point%x-mvx
                         !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                         ! --> alle Nummern entsprechend Vertices-Gitter 
                         IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrP
                         IndexBordBounds(b_upbounds,2)=ix
                         IndexBordBounds(b_upbounds,3)=j
                         IndexBordBounds(b_upbounds,4)=k
                       END IF
                     END DO
                   END DO
                   ix=ix0     !aktuellRandPointEcx
                   jx=ix0+1   !grid-'owest'-Point, spiegle Points der Celle ix0+1
                   DO j=iy0,iy1
                     DO k=iz0,iz1
                       IF(UpBoundsLayer(jx,j,k)%BCell%nrEcx>0) THEN        
                         nrP=UpBoundsLayer(jx,j,k)%BCell%nrEcx
                         b_upbounds=b_upbounds+1
                         UpBoundsLayer(ix,j,k)%BCell%nrEcx=b_upbounds
                         VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                         VertOutBoundsBord(b_upbounds)%Point%x= &
                             & VertOutBounds(nrP)%Point%x-mvx
                         !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                         ! --> alle Nummern entsprechend Verts-EdgeX-GrundGitter 
                         IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrEcx
                         IndexBordBounds(b_upbounds,2)=ix
                         IndexBordBounds(b_upbounds,3)=j
                         IndexBordBounds(b_upbounds,4)=k
                       END IF
                     END DO
                   END DO
                   ix=ix0-1   !aktuellRandPointEcy
                   jx=ix0+1   !grid-'owest'-Point, spiegle Points der Celle ix0+1
                   DO j=iy0+1,iy1
                     DO k=iz0,iz1
                       IF(UpBoundsLayer(jx,j,k)%BCell%nrEcy>0) THEN        
                         nrP=UpBoundsLayer(jx,j,k)%BCell%nrEcy
                         b_upbounds=b_upbounds+1
                         UpBoundsLayer(ix,j,k)%BCell%nrEcy=b_upbounds
                         VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                         VertOutBoundsBord(b_upbounds)%Point%x= &
                              & VertOutBounds(nrP)%Point%x-mvx
                         !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                         ! --> alle Nummern entsprechend Verts-EdgeY-GrundGitter 
                         IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrEcy
                         IndexBordBounds(b_upbounds,2)=ix
                         IndexBordBounds(b_upbounds,3)=j
                         IndexBordBounds(b_upbounds,4)=k
                       END IF
                     END DO
                   END DO
                   ix=ix0-1   !aktuellRandPointEcz
                   jx=ix0+1   !grid-'owest'-Point, spiegle Points der Celle ix0+1 
                   DO j=iy0,iy1
                     DO k=iz0+1,iz1
                       IF(UpBoundsLayer(jx,j,k)%BCell%nrEcz>0) THEN        
                         nrP=UpBoundsLayer(jx,j,k)%BCell%nrEcz
                         b_upbounds=b_upbounds+1
                         UpBoundsLayer(ix,j,k)%BCell%nrEcz=b_upbounds
                         VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                         VertOutBoundsBord(b_upbounds)%Point%x= &
                             & VertOutBounds(nrP)%Point%x-mvx
                         !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                         ! --> alle Nummern entsprechend Verts-EdgeZ-GrundGitter 
                         IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrEcz
                         IndexBordBounds(b_upbounds,2)=ix
                         IndexBordBounds(b_upbounds,3)=j
                         IndexBordBounds(b_upbounds,4)=k
                       END IF
                     END DO
                   END DO
 
                ELSE IF (nType=="oe") THEN
                   !spiegle Celle-Points ix1, damit gleiche Volume Celle-ix1-1    
                   ix=ix1+1  !aktuellRandPoint
                   jx=ix1-1  !grid-'oeast'-Point, spiegle Points der Celle ix1
                   mvx=dx(ix1)+dx(ix1+1)
                   DO j=iy0,iy1
                     DO k=iz0,iz1
                       IF(UpBoundsLayer(jx,j,k)%BCell%nrP>0) THEN        
                         nrP=UpBoundsLayer(jx,j,k)%BCell%nrP
                         b_upbounds=b_upbounds+1
                         UpBoundsLayer(ix,j,k)%BCell%nrP=b_upbounds
                         VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                         VertOutBoundsBord(b_upbounds)%Point%x= &
                              & VertOutBounds(nrP)%Point%x+mvx
                         IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrP
                         IndexBordBounds(b_upbounds,2)=ix
                         IndexBordBounds(b_upbounds,3)=j
                         IndexBordBounds(b_upbounds,4)=k
                       END IF
                     END DO
                   END DO
                   ix=ix1+1  !aktuellRandPointEcx
                   jx=ix1    !grid-'oeast'-Point, spiegle Points der Celle ix1
                   DO j=iy0,iy1
                     DO k=iz0,iz1
                       IF(UpBoundsLayer(jx,j,k)%BCell%nrEcx>0) THEN        
                         nrP=UpBoundsLayer(jx,j,k)%BCell%nrEcx
                         b_upbounds=b_upbounds+1
                         UpBoundsLayer(ix,j,k)%BCell%nrEcx=b_upbounds
                         VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                         VertOutBoundsBord(b_upbounds)%Point%x= &
                             & VertOutBounds(nrP)%Point%x+mvx
                         IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrEcx
                         IndexBordBounds(b_upbounds,2)=ix
                         IndexBordBounds(b_upbounds,3)=j
                         IndexBordBounds(b_upbounds,4)=k
                       END IF
                     END DO
                   END DO
                   ix=ix1+1  !aktuellRandPointEcy
                   jx=ix1-1  !grid-'oeast'-Point, spiegle Points der Celle ix1
                   DO j=iy0+1,iy1
                     DO k=iz0,iz1
                       IF(UpBoundsLayer(jx,j,k)%BCell%nrEcy>0) THEN        
                         nrP=UpBoundsLayer(jx,j,k)%BCell%nrEcy
                         b_upbounds=b_upbounds+1
                         UpBoundsLayer(ix,j,k)%BCell%nrEcy=b_upbounds
                         VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                         VertOutBoundsBord(b_upbounds)%Point%x= &
                             & VertOutBounds(nrP)%Point%x+mvx
                         IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrEcy
                         IndexBordBounds(b_upbounds,2)=ix
                         IndexBordBounds(b_upbounds,3)=j
                         IndexBordBounds(b_upbounds,4)=k
                       END IF
                     END DO
                   END DO
                   ix=ix1+1  !aktuellRandPointEcz
                   jx=ix1-1  !grid-'oeast'-Point, spiegle Points der Celle ix1
                   DO j=iy0,iy1
                     DO k=iz0+1,iz1
                       IF(UpBoundsLayer(jx,j,k)%BCell%nrEcz>0) THEN        
                         nrP=UpBoundsLayer(jx,j,k)%BCell%nrEcz
                         b_upbounds=b_upbounds+1
                         UpBoundsLayer(ix,j,k)%BCell%nrEcz=b_upbounds
                         VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                         VertOutBoundsBord(b_upbounds)%Point%x= &
                              VertOutBounds(nrP)%Point%x+mvx
                         IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrEcz
                         IndexBordBounds(b_upbounds,2)=ix
                         IndexBordBounds(b_upbounds,3)=j
                         IndexBordBounds(b_upbounds,4)=k
                       END IF
                     END DO
                   END DO
                END IF  !"ow" "oe" 

             CASE ('p')
                IF(TypeW=="pw".AND.TypeE=="pe") THEN
                  !.................................
                  !pw und pe gleicher Funktionsverlauf am Rand
                  mvx=Vertices(ix1,0,0)%Point%x-Vertices(ix0,0,0)%Point%x
                  IF(nType=="pw") THEN
                     ix=ix0-1               !aktuellRandPoint
                     jx=Floor(ibn)%ix1-1    !Nachbar-Point
                     DO j=iy0,iy1
                       DO k=iz0,iz1
                         IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrP>0) THEN        
                           nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrP
                           b_upbounds=b_upbounds+1
                           UpBoundsLayer(ix,j,k)%BCell%nrP=b_upbounds
                           VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                           VertOutBoundsBord(b_upbounds)%Point%x= &
                                VertOutBounds(nrP)%Point%x-mvx
                           !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                           ! --> alle Nummern entsprechend Vertices-Gitter 
                           IndexBordBounds(b_upbounds,1)=0
                           IndexBordBounds(b_upbounds,2)=ix
                           IndexBordBounds(b_upbounds,3)=j
                           IndexBordBounds(b_upbounds,4)=k
                         END IF
                       END DO
                     END DO
                     ix=ix0    !aktuellRandPointEcx
                     jx=Floor(ibn)%ix1    !Nachbar-Point
                     DO j=iy0,iy1
                       DO k=iz0,iz1
                         IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcx>0) THEN        
                           nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcx
                           b_upbounds=b_upbounds+1
                           UpBoundsLayer(ix,j,k)%BCell%nrEcx=b_upbounds
                           VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                           VertOutBoundsBord(b_upbounds)%Point%x= &
                                VertOutBounds(nrP)%Point%x-mvx
                           !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                           ! --> alle Nummern entsprechend Verts-EdgeX-GrundGitter 
                           IndexBordBounds(b_upbounds,1)=0
                           IndexBordBounds(b_upbounds,2)=ix
                           IndexBordBounds(b_upbounds,3)=j
                           IndexBordBounds(b_upbounds,4)=k
                         END IF
                       END DO
                     END DO
                     ix=ix0-1    !aktuellRandPointEcy
                     jx=Floor(ibn)%ix1-1    !Nachbar-Point
                     DO j=iy0+1,iy1
                       DO k=iz0,iz1
                         IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcy>0) THEN        
                           nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcy
                           b_upbounds=b_upbounds+1
                           UpBoundsLayer(ix,j,k)%BCell%nrEcy=b_upbounds
                           VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                           VertOutBoundsBord(b_upbounds)%Point%x= &
                                VertOutBounds(nrP)%Point%x-mvx
                           !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                           ! --> alle Nummern entsprechend Verts-EdgeY-GrundGitter 
                           IndexBordBounds(b_upbounds,1)=0
                           IndexBordBounds(b_upbounds,2)=ix
                           IndexBordBounds(b_upbounds,3)=j
                           IndexBordBounds(b_upbounds,4)=k
                         END IF
                       END DO
                     END DO
                     ix=ix0-1    !aktuellRandPointEcz
                     jx=Floor(ibn)%ix1-1    !Nachbar-Point
                     DO j=iy0,iy1
                       DO k=iz0+1,iz1
                         IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcz>0) THEN        
                           nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcz
                           b_upbounds=b_upbounds+1
                           UpBoundsLayer(ix,j,k)%BCell%nrEcz=b_upbounds
                           VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                           VertOutBoundsBord(b_upbounds)%Point%x= &
                                VertOutBounds(nrP)%Point%x-mvx
                           !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                           ! --> alle Nummern entsprechend Verts-EdgeZ-GrundGitter 
                           IndexBordBounds(b_upbounds,1)=0
                           IndexBordBounds(b_upbounds,2)=ix
                           IndexBordBounds(b_upbounds,3)=j
                           IndexBordBounds(b_upbounds,4)=k
                         END IF
                       END DO
                     END DO
                  ELSE ! "pe"
                     ix=ix1+1  !aktuellRandPoint
                     jx=Floor(ibn)%ix0+1  !nachbar-Point
                     DO j=iy0,iy1
                       DO k=iz0,iz1
                         IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrP>0) THEN        
                           nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrP
                           b_upbounds=b_upbounds+1
                           UpBoundsLayer(ix,j,k)%BCell%nrP=b_upbounds
                           VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                           VertOutBoundsBord(b_upbounds)%Point%x= &
                                VertOutBounds(nrP)%Point%x+mvx
                           IndexBordBounds(b_upbounds,1)=0
                           IndexBordBounds(b_upbounds,2)=ix
                           IndexBordBounds(b_upbounds,3)=j
                           IndexBordBounds(b_upbounds,4)=k
                         END IF
                       END DO
                     END DO
                     ix=ix1+1  !aktuellRandPointEcx
                     jx=Floor(ibn)%ix0+1  !nachbar-Point
                     DO j=iy0,iy1
                       DO k=iz0,iz1
                         IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcx>0) THEN        
                           nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcx
                           b_upbounds=b_upbounds+1
                           UpBoundsLayer(ix,j,k)%BCell%nrEcx=b_upbounds
                           VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                           VertOutBoundsBord(b_upbounds)%Point%x= &
                                VertOutBounds(nrP)%Point%x+mvx
                           IndexBordBounds(b_upbounds,1)=0
                           IndexBordBounds(b_upbounds,2)=ix
                           IndexBordBounds(b_upbounds,3)=j
                           IndexBordBounds(b_upbounds,4)=k
                         END IF
                       END DO
                     END DO
                     ix=ix1+1  !aktuellRandPointEcy
                     jx=Floor(ibn)%ix0+1  !nachbar-Point
                     DO j=iy0+1,iy1
                       DO k=iz0,iz1
                         IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcy>0) THEN        
                           nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcy
                           b_upbounds=b_upbounds+1
                           UpBoundsLayer(ix,j,k)%BCell%nrEcy=b_upbounds
                           VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                           VertOutBoundsBord(b_upbounds)%Point%x= &
                                VertOutBounds(nrP)%Point%x+mvx
                           IndexBordBounds(b_upbounds,1)=0
                           IndexBordBounds(b_upbounds,2)=ix
                           IndexBordBounds(b_upbounds,3)=j
                           IndexBordBounds(b_upbounds,4)=k
                         END IF
                       END DO
                     END DO
                     ix=ix1+1  !aktuellRandPointEcz
                     jx=Floor(ibn)%ix0+1  !nachbar-Point
                     DO j=iy0,iy1
                       DO k=iz0+1,iz1
                         IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcz>0) THEN        
                           nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcz
                           b_upbounds=b_upbounds+1
                           UpBoundsLayer(ix,j,k)%BCell%nrEcz=b_upbounds
                           VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                           VertOutBoundsBord(b_upbounds)%Point%x= &
                                VertOutBounds(nrP)%Point%x+mvx
                           IndexBordBounds(b_upbounds,1)=0
                           IndexBordBounds(b_upbounds,2)=ix
                           IndexBordBounds(b_upbounds,3)=j
                           IndexBordBounds(b_upbounds,4)=k
                         END IF
                       END DO
                     END DO
                  END IF  ! If(w oder e)
                END IF  ! (TypeW=='pw'.AND.TypeE='pe')

             CASE ('i')  ! "iw","ie"
                ! z.Zt. gleiches Refine, 
                IF(nType=="iw") THEN
                   ix=ix0-1               !aktuellRandPoint
                   jx=Floor(ibn)%ix1-1    !Nachbar-Point
                   fac_ra=1
                   fac_n =1
                ELSE !'ie'
                   ix=ix1+1              !aktuellRandCelle
                   jx=Floor(ibn)%ix0+1   !NachbarCelle
                   fac_ra=0
                   fac_n =0
                END IF

                !RandPoint-Grid
                DO j=iy0,iy1
                  DO k=iz0,iz1
                    IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrP>0) THEN        
                      nrP=Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrP
                      UpBoundsLayer(ix,j,k)%BCell%nrP=nrP
                      !wenn refine ungleich Abfrage ob Grid-Point gleich
                      !b_upbounds=b_upbounds+1
                      !UpBoundsLayer(ix,j,k)%BCell%nrP=b_upbounds
                      !...
                      !VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP) !Point-Grid-Info
                      !VertOutBoundsBord(b_upbounds)%Point%x=  ! wenn Refine hinzu
                      !...
                      !IndexBordBounds(b_upbounds,1)=UpBoundsLayer(ix,j,k)%BCell%nrP
                      !IndexBordBounds(b_upbounds,2)=ix
                      !IndexBordBounds(b_upbounds,3)=j
                      !IndexBordBounds(b_upbounds,4)=k
                    END IF
                  END DO
                END DO
                !RandPoint-Ecx
                DO j=iy0,iy1
                  DO k=iz0,iz1
                    IF(Floor(ibn)%UpBoundsLayer(jx+fac_n,j,k)%BCell%nrEcx>0) THEN        
                      UpBoundsLayer(ix+fac_ra,j,k)%BCell%nrEcx= &
                        & Floor(ibn)%UpBoundsLayer(jx+fac_n,j,k)%BCell%nrEcx
                      ! wenn refine ungleich....
                      !b_upbounds=b_upbounds+1
                      !UpBoundsLayer(ix+fac_ra,j,k)%BCell%nrEcx=b_upbounds
                      !....
                      !VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                      !VertOutBoundsBord(b_upbounds)%Point%x=   ! wenn Refine hinzu
                      !VertOutBoundsBord(b_upbounds)%nrP
                      ! --> alle Nummern entsprechend Verts-EdgeX-GrundGitter 
                      !IndexBordBounds(b_upbounds,1)=0
                      !IndexBordBounds(b_upbounds,2)=ix
                      !IndexBordBounds(b_upbounds,3)=j
                      !IndexBordBounds(b_upbounds,4)=k
                    END IF
                  END DO
                END DO
                !RandPoint-Ecy
                DO j=iy0+1,iy1
                  DO k=iz0,iz1
                    IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcy>0) THEN        
                        UpBoundsLayer(ix,j,k)%BCell%nrEcy= &
                         & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcy
                      ! wenn refine ungleich....
                      !b_upbounds=b_upbounds+1
                      !UpBoundsLayer(ix,j,k)%BCell%nrEcy=b_upbounds
                      !VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                      !VertOutBoundsBord(b_upbounds)%Point%x= ! wenn Refine hinzu
                      !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                      ! --> alle Nummern entsprechend Verts-EdgeY-GrundGitter 
                      !IndexBordBounds(b_upbounds,1)=0
                      !IndexBordBounds(b_upbounds,2)=ix
                      !IndexBordBounds(b_upbounds,3)=j
                      !IndexBordBounds(b_upbounds,4)=k
                    END IF
                  END DO
                END DO
                !RandPoint-Ecz
                DO j=iy0,iy1
                  DO k=iz0+1,iz1
                    IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcz>0) THEN        
                       UpBoundsLayer(ix,j,k)%BCell%nrEcz= &
                         & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%nrEcz
                      ! wenn refine ungleich....
                      !b_upbounds=b_upbounds+1
                      !UpBoundsLayer(ix,j,k)%BCell%nrEcz=b_upbounds
                      !VertOutBoundsBord(b_upbounds)=VertOutBounds(nrP)
                      !VertOutBoundsBord(b_upbounds)%Point%x= !  wenn Refine hinzu
                      !VertOutBoundsBord(b_upbounds)%nrP,%nrInP,%nrCutP,%nrI 
                      ! --> alle Nummern entsprechend Verts-EdgeZ-GrundGitter 
                      !IndexBordBounds(b_upbounds,1)=0
                      !IndexBordBounds(b_upbounds,2)=ix
                      !IndexBordBounds(b_upbounds,3)=j
                      !IndexBordBounds(b_upbounds,4)=k
                    END IF
                  END DO
                END DO
 
             END SELECT  ! nType --> [o,p,i] -w,-e
   

         CASE ("n","s")
         !----------------------------------
             SELECT CASE (Nachbars(in)%nType(1:1)) 

             CASE ('o')  ! "os","on"
                  IF(nType=="os") THEN
                  ELSE ! "on"
                  END IF
             CASE ('p')  ! "pn","ps"
                IF(TypeS=="ps".AND.TypeN=="pn") THEN
                  IF(nType=="ps") THEN
                  ELSE ! "pn"
                  END IF
                END IF
             CASE ('i')  ! "in","is"

             END SELECT  ! nType --> [o,p,i] -n,-s
         CASE ("t","b")
         !----------------------------------
             SELECT CASE (Nachbars(in)%nType(1:1)) 

             CASE ('o')  ! "ob","ot"
                  IF(nType=="ob") THEN
                  ELSE ! "ot"
                  END IF
             CASE ('p')  ! "pb","pt"
                IF(TypeB=="pb".AND.TypeT=="pt") THEN
                  IF(nType=="pb") THEN
                  ELSE ! "pt"
                  END IF
                END IF
             CASE ('i')  ! "in","is"

             END SELECT  ! nType --> [o,p,i] -b,-t

         END SELECT  !'w/e','n/s','t/b' -> to boundslayer jeweils [o,p,i]
       END DO  ! AnzahlNachbar
     END DO  !ib
     !Write(*,*) "SortIndexOutUpBounds : " 
     !Write(*,*) "b_upbounds   = ",b_upbounds, " allokiert znr_upbounds+per_bounds = ",znr_upbounds+per_bounds 
     !Write(*,*) "znr_upbounds = ",znr_upbounds
     END IF 
END SUBROUTINE SortIndexOutUpBounds

SUBROUTINE WriteBoundBorderWE_EdgesX
   INTEGER :: in,ix,j,k,is,liType,li_max
   INTEGER :: wbnex,wbney,wbnez,wbnec,nrP,nrPec(1:2)
   wbnex=0;wbney=0;wbnez=0;wbnec=0
   li_max=0
   !-- Border-Output--Edge_X
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
          !..........................................
          !z.Zt. noch gleiches Refine wenn iw || ie
          !SELECT CASE (Nachbars(in)%nType(1:1))
          !   CASE ("o","p").....CASE("i")
          !..........................................
          IF(nType(2:2)=='w') THEN
              ix=ix0
              wbnex=0
          ELSE !'e' 
              ix=ix1+1
              wbnex=0
          END IF
          DO j=iy0,iy1
            DO k=iz0,iz1
              IF(UpBoundsLayer(ix,j,k)%BCell%bexnr>0) THEN
                 wbnex=wbnex+1
                 IF(UpBoundsLayer(ix,j,k)%BCell%nrEcx>0) THEN
                   IF(UpBoundsLayer(ix-1,j,k)%BCell%nrP==0) THEN
                      CALL WOutBEdge(wbnex,ix,j,k,"ex", UpBoundsLayer(ix,j,k)%BCell%bexnr, &
                             & UpBoundsLayer(ix  ,j,k)%BCell%nrEcx,  &
                             & UpBoundsLayer(ix,j,k)%BCell%nrP) 
                   ELSE IF (UpBoundsLayer(ix,j,k)%BCell%nrP==0) THEN
                      CALL WOutBEdge(wbnex,ix,j,k,"ex", UpBoundsLayer(ix,j,k)%BCell%bexnr, &
                             & UpBoundsLayer(ix-1,j,k)%BCell%nrP, &
                             & UpBoundsLayer(ix  ,j,k)%BCell%nrEcx) 
                   END IF
                 ELSE
                   CALL WOutBEdge(wbnex,ix,j,k,"ex", UpBoundsLayer(ix,j,k)%BCell%bexnr, &
                             & UpBoundsLayer(ix-1,j,k)%BCell%nrP, &
                             & UpBoundsLayer(ix  ,j,k)%BCell%nrP)
                 END IF
              END IF
            END DO 
          END DO 
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   ! [o,p,i]- 'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar

END SUBROUTINE WriteBoundBorderWE_EdgesX
       
SUBROUTINE WriteBoundBorderWE_EdgesY
   INTEGER :: in,ix,j,k,is,liType,li_max
   INTEGER :: wbney,nrP,nrPec(1:2)
   wbney=0
   li_max=0
   !-- Border-Output--Edge_Y
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
          !z.Zt. noch gleiches Refine wenn iw || ie
          IF(nType(2:2)=='w') THEN
             ix=ix0-1
             wbney=0
          ELSE !'e' 
             ix=ix1+1
             wbney=0
          END IF
          DO j=iy0+1,iy1
            DO k=iz0,iz1
              IF(UpBoundsLayer(ix,j,k)%BCell%beynr>0) THEN
                 wbney=wbney+1
                 IF(UpBoundsLayer(ix,j,k)%BCell%nrEcy>0) THEN
                    IF(UpBoundsLayer(ix,j-1,k)%BCell%nrP==0) THEN
                       CALL WOutBEdge(wbney,ix,j,k,"ey", UpBoundsLayer(ix,j,k)%BCell%beynr, &
                             & UpBoundsLayer(ix,j,k)%BCell%nrEcy, &
                             & UpBoundsLayer(ix,j,k)%BCell%nrP)
                    ELSE IF (UpBoundsLayer(ix,j,k)%BCell%nrP==0) THEN
                       CALL WOutBEdge(wbney,ix,j,k,"ey", UpBoundsLayer(ix,j,k)%BCell%beynr, &
                             & UpBoundsLayer(ix,j-1,k)%BCell%nrP, &
                             & UpBoundsLayer(ix,j,k)%BCell%nrEcy)
                    END IF
                 ELSE  
                    CALL WOutBEdge(wbney,ix,j,k,"ey", UpBoundsLayer(ix,j,k)%BCell%beynr, &
                             & UpBoundsLayer(ix,j-1,k)%BCell%nrP, &
                             & UpBoundsLayer(ix,j  ,k)%BCell%nrP)
                 END IF
              END IF
            END DO 
          END DO 
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   ! [o,p,i]- 'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar

END SUBROUTINE WriteBoundBorderWE_EdgesY


SUBROUTINE WriteBoundBorderWE_EdgesZ
   INTEGER :: in,ix,j,k,is,liType,li_max
   INTEGER :: wbnez,nrP,nrPec(1:2)
   wbnez=0
   li_max=0
   !-- Border-Output--Edge_Z
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
         !z.Zt. noch gleiches Refine wenn iw || ie
         IF(nType(2:2)=='w') THEN
            ix=ix0-1
            wbnez=0
         ELSE !'e' 
            ix=ix1+1
            wbnez=0
         END IF
         DO j=iy0,iy1
            DO k=iz0+1,iz1
              IF(UpBoundsLayer(ix,j,k)%BCell%beznr>0) THEN
                wbnez=wbnez+1
                IF(UpBoundsLayer(ix,j,k)%BCell%nrEcz>0) THEN
                  IF(UpBoundsLayer(ix,j,k-1)%BCell%nrP==0) THEN
                      CALL WOutBEdge(wbnez,ix,j,k,"ez", UpBoundsLayer(ix,j,k)%BCell%beznr, &
                            & UpBoundsLayer(ix,j,k)%BCell%nrEcz, &
                              UpBoundsLayer(ix,j,k)%BCell%beznr)
                  ELSE IF(UpBoundsLayer(ix,j,k)%BCell%nrP==0) THEN
                      CALL WOutBEdge(wbnez,ix,j,k,"ez", UpBoundsLayer(ix,j,k)%BCell%beznr, &
                            & UpBoundsLayer(ix,j,k-1)%BCell%nrP, &
                            & UpBoundsLayer(ix,j,k)%BCell%nrEcz)
                  END IF
                ELSE
                  CALL WOutBEdge(wbnez,ix,j,k,"ez", UpBoundsLayer(ix,j,k)%BCell%beznr, &
                            & UpBoundsLayer(ix,j,k-1)%BCell%nrP, &
                            & UpBoundsLayer(ix,j,k  )%BCell%nrP)
                END IF
              END IF
           END DO 
         END DO 
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   ! [o,p,i]- 'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar
END SUBROUTINE WriteBoundBorderWE_EdgesZ

SUBROUTINE WriteBoundBorderWE_EdgesEc
   INTEGER :: in,ix,j,k,is,liType,li_max
   INTEGER :: wbnec,nrP,nrPec(1:2)
   wbnec=0
   li_max=0
   !-- Border-Output--Edge_ec
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
        SELECT CASE (Nachbars(in)%nType(1:1))
        ! 'o,p,i' -> west,east  z.Zt: nur Refine =0, Incr->X,Y,Z noch offen 
        !------------------------TypeW=="ow".OR.TypeE=="oe"-------------------
        CASE ('o','p')
           !BCell%becnrxy.....
           IF (nType(2:2)=='w') THEN
               ix=ix0
               wbnec=0
           ELSE !'e' 
               ix=ix1+1
               wbnec=0
           END IF
           DO j=iy0+1,iy1
             DO k=iz0,iz1
               IF (UpBoundsLayer(ix,j,k)%BCell%becnrxy>0) THEN
                  !nrPec(1:2) ermitteln
                  is=1
                  IF(UpBoundsLayer(ix,j-1,k)%BCell%nrEcx>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j-1,k)%BCell%nrEcx
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix,j,k)%BCell%nrEcy>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j,k)%BCell%nrEcy
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix,j,k)%BCell%nrEcx>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j,k)%BCell%nrEcx
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix-1,j,k)%BCell%nrEcy>0) THEN
                     nrPec(is)=UpBoundsLayer(ix-1,j,k)%BCell%nrEcy
                     is=is+1
                  END IF
                  wbnec=wbnec+1
                  Write(OutUnitB,*) leerzei6,"ec",UpBoundsLayer(ix,j,k)%BCell%becnrxy,leerzei2, &
                         &   nrPec(1), nrPec(2),leerzei6,wbnec
                  Write(OutUnitB,*) ix,j,k       
                  Write(OutUnitBpv,*) "ec",UpBoundsLayer(ix,j,k)%BCell%becnrxy,leerzei2, &
                         &   nrPec(1), nrPec(2)
                  !Write(OutUnitBpv,*) "ec",bl_ec(i,1),leerzei2,bl_ec(i,2),bl_ec(i,3)
               END IF
             END DO 
           END DO 
           !BCell%becnryz.....
           IF(nType(2:2)=='w') THEN
              ix=ix0-1
           ELSE !'e' 
              ix=ix1+1
           END IF
           DO j=iy0+1,iy1
              DO k=iz0+1,iz1
               IF (UpBoundsLayer(ix,j,k)%BCell%becnryz>0) THEN
                  !nrPec(1:2) ermitteln
                  is=1
                  IF(UpBoundsLayer(ix,j,k-1)%BCell%nrEcy>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j,k-1)%BCell%nrEcy
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix,j,k)%BCell%nrEcz>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j,k)%BCell%nrEcz
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix,j,k)%BCell%nrEcy>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j,k)%BCell%nrEcy
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix,j-1,k)%BCell%nrEcz>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j-1,k)%BCell%nrEcz
                     is=is+1
                  END IF
                  wbnec=wbnec+1
                  Write(OutUnitB,*) leerzei6,"ec",UpBoundsLayer(ix,j,k)%BCell%becnryz,leerzei2, &
                         &   nrPec(1),nrPec(2),leerzei6,wbnec
                  Write(OutUnitB,*) ix,j,k       
                  Write(OutUnitBpv,*) "ec",UpBoundsLayer(ix,j,k)%BCell%becnryz,leerzei2, &
                         &   nrPec(1), nrPec(2)
               END IF
             END DO 
           END DO 
           !BCell%becnrzx..............................
           IF(nType(2:2)=='w') THEN
               ix=ix0
           ELSE !'e' 
               ix=ix1+1
           END IF
           DO j=iy0,iy1
             DO k=iz0+1,iz1
               IF (UpBoundsLayer(ix,j,k)%BCell%becnrzx>0) THEN
                  !nrPec(1:2) ermitteln andere Variante is kann überlaufen
                  ! Face%EdgeCut checken über diesen nrP in VetOut-Liste suchen!???
                  is=1
                  IF(UpBoundsLayer(ix,j,k-1)%BCell%nrEcx>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j,k-1)%BCell%nrEcx
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix,j,k)%BCell%nrEcz>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j,k)%BCell%nrEcz
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix,j,k)%BCell%nrEcx>0) THEN
                     nrPec(is)=UpBoundsLayer(ix,j,k)%BCell%nrEcx
                     is=is+1
                  END IF
                  IF(UpBoundsLayer(ix-1,j,k)%BCell%nrEcz>0) THEN
                     nrPec(is)=UpBoundsLayer(ix-1,j,k)%BCell%nrEcz
                     is=is+1
                  END IF
                  wbnec=wbnec+1
                  Write(OutUnitB,*) leerzei6,"ec",UpBoundsLayer(ix,j,k)%BCell%becnrzx,leerzei2, &
                         &   nrPec(1),nrPec(2),leerzei6,wbnec
                  Write(OutUnitB,*) ix,j,k
                  Write(OutUnitBpv,*) "ec",UpBoundsLayer(ix,j,k)%BCell%becnrzx,leerzei2, &
                         &   nrPec(1), nrPec(2)
               END IF
             END DO
           END DO
        CASE ('i')
         !z.Zt. noch gleiches Refine wenn iw || ie
           !BCell%becnrxy.....
           IF (nType(2:2)=='w') THEN
               ix=ix0
               wbnec=0
           ELSE !'e' 
               ix=ix1+1
               wbnec=0
           END IF
           !BCell%becnrxy.....
           DO j=iy0+1,iy1
             DO k=iz0,iz1
               IF (Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnrxy>0) THEN
                  !nrPec(1:2) ermitteln
                  nrPec(1)=NrPOutBounds(Floor(ibn)%Faces_XY(ix,j,k)%Face%EdgeCut(1))
                  nrPec(2)=NrPOutBounds(Floor(ibn)%Faces_XY(ix,j,k)%Face%EdgeCut(2))
                  wbnec=wbnec+1
                  Write(OutUnitB,*) leerzei6,"ec",Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnrxy,&
                         & leerzei2, nrPec(1), nrPec(2),leerzei6,wbnec
                  Write(OutUnitB,*) ix,j,k       
                  Write(OutUnitBpv,*) "ec",Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnrxy,&
                         & leerzei2, nrPec(1), nrPec(2)
               END IF
             END DO 
           END DO 
           !BCell%becnryz.....
           IF(nType(2:2)=='w') THEN
              ix=ix0-1
           ELSE !'e' 
              ix=ix1+1
           END IF
           DO j=iy0+1,iy1
              DO k=iz0+1,iz1
               IF (Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnryz>0) THEN
                  nrPec(1)=NrPOutBounds(Floor(ibn)%Faces_YZ(ix,j,k)%Face%EdgeCut(1))
                  nrPec(2)=NrPOutBounds(Floor(ibn)%Faces_YZ(ix,j,k)%Face%EdgeCut(2))
                  wbnec=wbnec+1
                  Write(OutUnitB,*) leerzei6,"ec",Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnryz,leerzei2, &
                         &   nrPec(1),nrPec(2),leerzei6,wbnec
                  Write(OutUnitB,*) ix,j,k       
                  Write(OutUnitBpv,*) "ec",Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnryz,leerzei2, &
                         &   nrPec(1),nrPec(2)
               END IF
             END DO 
           END DO 
           !BCell%becnrzx..............................
           IF(nType(2:2)=='w') THEN
               ix=ix0
           ELSE !'e' 
               ix=ix1+1
           END IF
           DO j=iy0,iy1
             DO k=iz0+1,iz1
               IF (Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnrzx>0) THEN
                  nrPec(1)=NrPOutBounds(Floor(ibn)%Faces_ZX(ix,j,k)%Face%EdgeCut(1))
                  nrPec(2)=NrPOutBounds(Floor(ibn)%Faces_ZX(ix,j,k)%Face%EdgeCut(2))
                  wbnec=wbnec+1
                  Write(OutUnitB,*) leerzei6,"ec",Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnrzx,leerzei2, &
                         &   nrPec(1),nrPec(2),leerzei6,wbnec
                  Write(OutUnitB,*) ix,j,k
                  Write(OutUnitBpv,*) "ec",Floor(ibn)%UpBoundsLayer(ix,j,k)%BCell%becnrzx,leerzei2, &
                         &   nrPec(1),nrPec(2)
               END IF
             END DO
           END DO
        END SELECT   ! [o,p,i]- 'w/e' -> to boundsborderlayer
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   ! 'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar

END SUBROUTINE WriteBoundBorderWE_EdgesEc


!#SUBROUTINE WriteBoundBorderWE_FacesFc
!#   INTEGER :: in,ix,i,j,k,jx,liType,li_max
!#   INTEGER :: wbfcnr,wbfxynr,wbfzxnr,wbfyznr
!#   INTEGER :: fegnr,fliste_egnr(1:8),ec_set,bfxynr,bfzxnr,bfyznr
!#   INTEGER :: fcegl_nb(1:8)
!#   INTEGER :: ai,aj,ak,li
!#   wbfcnr=0;wbfxynr=0;wbfzxnr=0;wbfyznr=0
!#   li_max=0
!#   !...%BCell%FCutNr
!#   DO in=1,AnzahlNachbar
!#      CALL Set(Nachbars(in))
!#      SELECT CASE (Nachbars(in)%nType(2:2))
!#        !----------------------------------
!#      CASE ("w","e")
!#      !-------------
!#        SELECT CASE (Nachbars(in)%nType(1:1))
!#        ! 'o,p,i' -> west,east  z.Zt: nur Refine =0, Incr->X,Y,Z noch offen 
!#        !------------------------TypeW=="ow".OR.TypeE=="oe"-------------------
!#        CASE ('o','p')
!#           IF(nType=="ow".OR.nType=="pw") THEN
!#               ix=ix0
!#               jx=ix1
!#               wbfcnr=0
!#           ELSE !'oe'.'pe' 
!#               ix=ix1+1
!#               jx=ix0+1
!#               wbfcnr=0
!#           END IF
!#           !...%BCell%FCutNr
!#           DO j=iy0+1,iy1
!#             DO k=iz0+1,iz1
!#               IF(UpBoundsLayer(ix,j,k)%BCell%FCutNr>0) THEN
!#                 i=0
!#                 !-ZX----
!#                 IF(UpBoundsLayer(ix,j,k)%BCell%becnrzx>0) THEN
!#                    i=i+1
!#                    UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
!#                           & UpBoundsLayer(ix,j,k)%BCell%becnrzx
!#                 ELSE 
!#                   IF(UpBoundsLayer(ix,j,k)%BCell%bexnr>0) THEN
!#                     i=i+1
!#                     UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
!#                           & UpBoundsLayer(ix,j,k)%BCell%bexnr
!#                   ELSE
!#                     IF(UpBoundsLayer(ix,j,k-1)%BCell%bexnr>0) THEN
!#                       i=i+1
!#                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
!#                           & UpBoundsLayer(ix,j,k-1)%BCell%bexnr
!#                     END IF
!#                   END IF
!#                 END IF
!#                 !-YZ----
!#                 IF(UpBoundsLayer(ix-1,j,k)%BCell%becnryz>0) THEN
!#                    i=i+1
!#                    UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
!#                         & UpBoundsLayer(ix-1,j,k)%BCell%becnryz 
!#                 ELSE
!#                   IF(UpBoundsLayer(ix-1,j,k)%BCell%beynr>0) THEN
!#                      i=i+1
!#                      UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
!#                         & UpBoundsLayer(ix-1,j,k)%BCell%beynr 
!#                   ELSE 
!#                     IF(UpBoundsLayer(ix-1,j,k-1)%BCell%beynr>0) THEN
!#                       i=i+1
!#                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
!#                         & UpBoundsLayer(ix-1,j,k-1)%BCell%beynr 
!#                     END IF
!#                   END IF
!#                 END IF
!#                 !--ZX---
!#                 IF(UpBoundsLayer(ix,j-1,k)%BCell%becnrzx>0) THEN
!#                    i=i+1
!#                    UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
!#                         & UpBoundsLayer(ix,j-1,k)%BCell%becnrzx
!#                 ELSE
!#                   IF(UpBoundsLayer(ix,j-1,k)%BCell%bexnr>0) THEN
!#                     i=i+1
!#                     UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
!#                         & UpBoundsLayer(ix,j-1,k)%BCell%bexnr
!#        END SELECT   ! [o,p,i]- 'w/e' -> to boundsborderlayer
!#      CASE ("n","s")
!#      !----------------------------------
!#          IF (nType=="on") THEN
!#          ELSE IF (nType=="os") THEN
!#          ELSE
!#          END IF
!#      CASE ("t","b")
!#      !----------------------------------
!#          IF (nType=="ot") THEN
!#          ELSE IF (nType=="ob") THEN
!#          ELSE
!#          END IF
!#      END SELECT   ! 'w/e','n/s','t/b' -> to boundsborderlayer
!#   END DO  ! AnzahlNachbar
!#
!#END SUBROUTINE WriteBoundBorderWE_EdgesEc


SUBROUTINE WriteBoundBorderWE_FacesFc
   INTEGER :: in,ix,i,j,k,jx,liType,li_max
   INTEGER :: wbfcnr,wbfxynr,wbfzxnr,wbfyznr
   INTEGER :: fegnr,fliste_egnr(1:8),ec_set,bfxynr,bfzxnr,bfyznr
   INTEGER :: fcegl_nb(1:8)
   INTEGER :: ai,aj,ak,li
   wbfcnr=0;wbfxynr=0;wbfzxnr=0;wbfyznr=0
   li_max=0
   !...%BCell%FCutNr
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
      !-------------
        SELECT CASE (Nachbars(in)%nType(1:1))
        ! 'o,p,i' -> west,east  z.Zt: nur Refine =0, Incr->X,Y,Z noch offen 
        !------------------------TypeW=="ow".OR.TypeE=="oe"-------------------
        CASE ('o','p')
           IF(nType=="ow".OR.nType=="pw") THEN
               ix=ix0
               jx=ix1
               wbfcnr=0
           ELSE !'oe'.'pe' 
               ix=ix1+1
               jx=ix0+1
               wbfcnr=0
           END IF
           !...%BCell%FCutNr
           DO j=iy0+1,iy1
             DO k=iz0+1,iz1
               IF(UpBoundsLayer(ix,j,k)%BCell%FCutNr>0) THEN
                 i=0
                 !-ZX----
                 IF(UpBoundsLayer(ix,j,k)%BCell%becnrzx>0) THEN
                    i=i+1
                    UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                           & UpBoundsLayer(ix,j,k)%BCell%becnrzx
                 ELSE 
                   IF(UpBoundsLayer(ix,j,k)%BCell%bexnr>0) THEN
                     i=i+1
                     UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                           & UpBoundsLayer(ix,j,k)%BCell%bexnr
                   ELSE
                     IF(UpBoundsLayer(ix,j,k-1)%BCell%bexnr>0) THEN
                       i=i+1
                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                           & UpBoundsLayer(ix,j,k-1)%BCell%bexnr
                     END IF
                   END IF
                 END IF
                 !-YZ----
                 IF(UpBoundsLayer(ix-1,j,k)%BCell%becnryz>0) THEN
                    i=i+1
                    UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                         & UpBoundsLayer(ix-1,j,k)%BCell%becnryz 
                 ELSE
                   IF(UpBoundsLayer(ix-1,j,k)%BCell%beynr>0) THEN
                      i=i+1
                      UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                         & UpBoundsLayer(ix-1,j,k)%BCell%beynr 
                   ELSE 
                     IF(UpBoundsLayer(ix-1,j,k-1)%BCell%beynr>0) THEN
                       i=i+1
                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                         & UpBoundsLayer(ix-1,j,k-1)%BCell%beynr 
                     END IF
                   END IF
                 END IF
                 !--ZX---
                 IF(UpBoundsLayer(ix,j-1,k)%BCell%becnrzx>0) THEN
                    i=i+1
                    UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                         & UpBoundsLayer(ix,j-1,k)%BCell%becnrzx
                 ELSE
                   IF(UpBoundsLayer(ix,j-1,k)%BCell%bexnr>0) THEN
                     i=i+1
                     UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                         & UpBoundsLayer(ix,j-1,k)%BCell%bexnr
                   ELSE 
                     IF(UpBoundsLayer(ix,j-1,k-1)%BCell%bexnr>0) THEN
                       i=i+1
                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                          & UpBoundsLayer(ix,j-1,k-1)%BCell%bexnr
                     END IF
                   END IF
                 END IF
                 !--YZ---
                 IF(UpBoundsLayer(ix,j,k)%BCell%becnryz>0) THEN
                    i=i+1
                    UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                         & UpBoundsLayer(ix,j,k)%BCell%becnryz
                 ELSE
                   IF(UpBoundsLayer(ix,j,k)%BCell%beynr>0) THEN
                      i=i+1
                      UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                         & UpBoundsLayer(ix,j,k)%BCell%beynr
                   ELSE 
                     IF(UpBoundsLayer(ix,j,k-1)%BCell%beynr>0) THEN
                       i=i+1
                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
                          & UpBoundsLayer(ix,j,k-1)%BCell%beynr
                     END IF
                   END IF
                 END IF
                 UpBoundsLayer(ix,j,k)%BCell%fceg_nr=i
                
                 wbfcnr=wbfcnr+1
                 CALL WOutBFace(wbfcnr,ix,j,k,"fc ", &
                       &  UpBoundsLayer(ix,j,k)%BCell%FCutNr, &
                       &  UpBoundsLayer(ix,j,k)%BCell%fceg_nr,tropos_border, &
                       &  UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(1:8))
               END IF          
             END DO
           END DO
        !-----------------------TypeW=="iw".OR.TypeE=="ie"---------------------
        CASE ('i')
           IF(nType=="iw") THEN
               ix=ix0
           ELSE !'ie' 
               ix=ix1+1
           END IF
           !z.Zt. noch gleiches Refine wenn iw || ie
           DO j=iy0+1,iy1
             DO k=iz0+1,iz1
               IF(UpBoundsLayer(ix,j,k)%BCell%FCutNr>0) THEN
   
                 wbfcnr=wbfcnr+1
                 CALL WOutBFace(wbfcnr,ix,j,k,"fc ", &
                       &  UpBoundsLayer(ix,j,k)%BCell%FCutNr, &
                       &  UpBoundsLayer(ix,j,k)%BCell%fceg_nr,tropos_border, &
                       &  UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(1:8))
               END IF          
             END DO
           END DO
         END SELECT  ! [o,i]- w/e
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   !'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar

   !!   
   !!            !nachbar-facecut edges suchen , entsprechend edge-positionen
   !!            !akt-edge-nr auf fc-rand-Cut setzen
   !!            fcegl_nb(1:8)=UpBoundsLayer(jx,j,k)%BCell%fceg_nrcut(1:8)
   !!            DO i=1,UpBoundsLayer(jx,j,k)%BCell%fceg_nr
   !!               ai=jx !suche bexnr-becnrxy-neighbor
   !!               DO ak=k-1,k
   !!                 DO aj=j-1,j
   !!                     IF (fcegl_nb(i)==UpBoundsLayer(ai,aj,ak)%BCell%bexnr) THEN
   !!                         UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                          & UpBoundsLayer(ix,aj,ak)%BCell%bexnr
   !!                     END IF    
   !!                 END DO
   !!                 IF (fcegl_nb(i)==UpBoundsLayer(ai,j,ak)%BCell%becnrxy) THEN
   !!                     UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                      & UpBoundsLayer(ix,j,ak)%BCell%becnrxy
   !!                 END IF    
   !!               END DO
   !!               aj=j !suche beynr-becnryz-neighbor
   !!               DO ai=jx-1,jx 
   !!                 DO ak=k-1,k
   !!                    IF (fcegl_nb(i)==UpBoundsLayer(ai,aj,ak)%BCell%beynr) THEN
   !!                        IF(ai==jx-1) THEN
   !!                          UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                           & UpBoundsLayer(ix-1,j,ak)%BCell%beynr
   !!                        ELSE !ai==jx
   !!                          UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                           & UpBoundsLayer(ix,j,ak)%BCell%beynr
   !!                        END IF
   !!                    END IF    
   !!                 END DO
   !!                 IF (fcegl_nb(i)==UpBoundsLayer(ai,j,k)%BCell%becnryz) THEN
   !!                     IF(ai==jx-1) THEN
   !!                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                        & UpBoundsLayer(ix-1,j,k)%BCell%becnryz
   !!                     ELSE !ai==jx
   !!                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                        & UpBoundsLayer(ix,j,k)%BCell%becnryz
   !!                     END IF
   !!                 END IF    
   !!               END DO
   !!               ak=k  !suche beznr-becnrzx-neighbor
   !!               DO aj=j-1,j
   !!                 DO ai=jx-1,jx
   !!                    IF (fcegl_nb(i)==UpBoundsLayer(ai,aj,k)%BCell%beznr) THEN
   !!                        IF(ai==jx-1) THEN
   !!                          UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                           & UpBoundsLayer(ix-1,aj,k)%BCell%beznr
   !!                        ELSE !ai==jx
   !!                          UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                          & UpBoundsLayer(ix,aj,k)%BCell%beznr
   !!                        END IF
   !!                    END IF    
   !!                 END DO
   !!                 IF (fcegl_nb(i)==UpBoundsLayer(jx,aj,k)%BCell%becnrzx) THEN
   !!                     IF(aj==j-1) THEN
   !!                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                        & UpBoundsLayer(ix,j-1,k)%BCell%becnrzx
   !!                     ELSE !aj==j
   !!                       UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(i)= &
   !!                        & UpBoundsLayer(ix,j,k)%BCell%becnrzx
   !!                     END IF
   !!                 END IF    
   !!               END DO
   !!            END DO
   !!            wbfcnr=wbfcnr+1
   !!            UpBoundsLayer(ix,j,k)%BCell%fceg_nr=UpBoundsLayer(jx,j,k)%BCell%fceg_nr
   !!            CALL WOutBFace(wbfcnr,ix,j,k,"fc ", &
   !!                  &  UpBoundsLayer(ix,j,k)%BCell%FCutNr, &
   !!                  &  UpBoundsLayer(ix,j,k)%BCell%fceg_nr,tropos_border, &
   !!                  &  UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(1:8))

END SUBROUTINE WriteBoundBorderWE_FacesFc

SUBROUTINE WriteBoundBorderWE_FacesXY
   INTEGER :: in,ix,i,j,k,jx,liType,li_max
   INTEGER :: wbfcnr,wbfxynr,wbfzxnr,wbfyznr
   INTEGER :: fegnr,fliste_egnr(1:8),ec_set,bfxynr,bfzxnr,bfyznr
   INTEGER :: fcegl_nb(1:8)
   INTEGER :: ai,aj,ak,li
   wbfxynr=0
   li_max=0
   !...%BCell%bfxynr
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
        !z.Zt. noch gleiches Refine wenn iw || ie
        IF(nType(2:2)=='w') THEN
            ix=ix0
            jx=ix1
            wbfxynr=0
        ELSE !'e' 
            ix=ix1+1
            jx=ix0+1
            wbfxynr=0
        END IF
        !...%BCell%bfxynr
        DO j=iy0+1,iy1
          DO k=iz0,iz1
            IF(UpBoundsLayer(ix,j,k)%BCell%bfxynr>0) THEN
              fegnr=0;fliste_egnr(1:8)=0;ec_set=0
              If( UpBoundsLayer(ix,j-1,k)%BCell%bexnr>0) THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix,j-1,k)%BCell%bexnr
              ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnrxy>0) THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrxy
                  ec_set=1
              END IF
              If( UpBoundsLayer(ix,j,k)%BCell%beynr>0) THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%beynr
              ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnrxy>0.and.ec_set==0) THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrxy
                  ec_set=1
              END IF
              If( UpBoundsLayer(ix,j,k)%BCell%bexnr>0)THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%bexnr
              ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnrxy>0.and.ec_set==0) THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrxy
                  ec_set=1
              END IF
              If( UpBoundsLayer(ix-1,j,k)%BCell%beynr>0) THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix-1,j,k)%BCell%beynr
              ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnrxy>0.and.ec_set==0) THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrxy
                  ec_set=1
              END IF
              IF(ec_set==0.and.UpBoundsLayer(ix,j,k)%BCell%becnrxy>0) THEN
                  fegnr=fegnr+1
                  fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrxy
              END IF   
              UpBoundsLayer(ix,j,k)%BCell%feg_nr=fegnr
              UpBoundsLayer(ix,j,k)%BCell%flisteg_nr(1:8)=fliste_egnr(1:8)   
              bfxynr=UpBoundsLayer(ix,j,k)%BCell%bfxynr
              wbfxynr=wbfxynr+1
              !tropo_id/bound_id noch offen
              CALL WOutBFace(wbfxynr,ix,j,k,"fxy", bfxynr, fegnr, 0, fliste_egnr(1:8))
            END IF          
          END DO
        END DO
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   ! [o,p,i]- 'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar
END SUBROUTINE WriteBoundBorderWE_FacesXY


SUBROUTINE WriteBoundBorderWE_FacesZX
   INTEGER :: in,ix,i,j,k,jx,liType,li_max
   INTEGER :: wbfzxnr
   INTEGER :: fegnr,fliste_egnr(1:8),ec_set,bfxynr,bfzxnr,bfyznr
   INTEGER :: fcegl_nb(1:8)
   INTEGER :: ai,aj,ak,li
   wbfzxnr=0
   li_max=0
   !...%BCell%bfzxnr
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
         !z.Zt. noch gleiches Refine wenn iw || ie
         IF(nType(2:2)=='w') THEN
             ix=ix0
             jx=ix1
             wbfzxnr=0
         ELSE !'e' 
             ix=ix1+1
             jx=ix0+1
             wbfzxnr=0
         END IF
         !...%BCell%bfzxnr
         DO j=iy0,iy1
           DO k=iz0+1,iz1
             IF(UpBoundsLayer(ix,j,k)%BCell%bfzxnr>0) THEN
               fegnr=0;fliste_egnr(1:8)=0;ec_set=0
               If( UpBoundsLayer(ix,j,k-1)%BCell%bexnr>0) THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k-1)%BCell%bexnr
               ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnrzx>0) THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrzx
                   ec_set=1
               END IF
               If( UpBoundsLayer(ix,j,k)%BCell%beznr>0) THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%beznr
               ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnrzx>0.and.ec_set==0) THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrzx
                   ec_set=1
               END IF
               If( UpBoundsLayer(ix,j,k)%BCell%bexnr>0) THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%bexnr
               ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnrzx>0.and.ec_set==0) THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrzx
                   ec_set=1
               END IF
               If( UpBoundsLayer(ix-1,j,k)%BCell%beznr>0) THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix-1,j,k)%BCell%beznr
               ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnrzx>0.and.ec_set==0) THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrzx
                   ec_set=1
               END IF
               IF(ec_set==0.and.UpBoundsLayer(ix,j,k)%BCell%becnrzx>0)THEN
                   fegnr=fegnr+1
                   fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnrzx
               END IF
               UpBoundsLayer(ix,j,k)%BCell%feg_nr=fegnr
               UpBoundsLayer(ix,j,k)%BCell%flisteg_nr(1:8)=fliste_egnr(1:8)   
               bfzxnr=UpBoundsLayer(ix,j,k)%BCell%bfzxnr
               wbfzxnr=wbfzxnr+1
               !tropo_id/bound_id noch offen
               CALL WOutBFace(wbfzxnr,ix,j,k,"fzx", bfzxnr, fegnr,0, fliste_egnr(1:8))
             END IF          
           END DO
         END DO
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   ! [o,p,i]- 'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar
END SUBROUTINE WriteBoundBorderWE_FacesZX

SUBROUTINE WriteBoundBorderWE_FacesYZ
   INTEGER :: in,ix,i,j,k,jx,liType,li_max
   INTEGER :: wbfyznr
   INTEGER :: fegnr,fliste_egnr(1:8),ec_set,bfxynr,bfzxnr,bfyznr
   INTEGER :: fcegl_nb(1:8)
   INTEGER :: ai,aj,ak,li
   wbfyznr=0
   li_max=0
   !...%BCell%bfyznr
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
          !z.Zt. noch gleiches Refine wenn iw || ie
          IF(nType(2:2)=='w') THEN
              ix=ix0-1
              jx=ix1-1
              wbfyznr=0
          ELSE !'e' 
              ix=ix1+1
              jx=ix0+1
              wbfyznr=0
          END IF
          !...%BCell%bfyznr
          DO j=iy0+1,iy1
            DO k=iz0+1,iz1
              IF(UpBoundsLayer(ix,j,k)%BCell%bfyznr>0) THEN
                fegnr=0
                fliste_egnr(1:8)=0
                ec_set=0
                If( UpBoundsLayer(ix,j,k-1)%BCell%beynr>0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k-1)%BCell%beynr
                ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnryz>0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnryz
                    ec_set=1
                END IF
                If( UpBoundsLayer(ix,j,k)%BCell%beznr>0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%beznr
                ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnryz>0.and.ec_set==0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnryz
                    ec_set=1
                END IF
                If( UpBoundsLayer(ix,j,k)%BCell%beynr>0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%beynr
                ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnryz>0.and.ec_set==0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnryz
                    ec_set=1
                END IF
                If( UpBoundsLayer(ix,j-1,k)%BCell%beznr>0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j-1,k)%BCell%beznr
                ELSE IF  (UpBoundsLayer(ix,j,k)%BCell%becnryz>0.and.ec_set==0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnryz
                    ec_set=1
                END IF
                IF(ec_set==0.and.UpBoundsLayer(ix,j,k)%BCell%becnryz>0) THEN
                    fegnr=fegnr+1
                    fliste_egnr(fegnr)=UpBoundsLayer(ix,j,k)%BCell%becnryz
                END IF 
                UpBoundsLayer(ix,j,k)%BCell%feg_nr=fegnr
                UpBoundsLayer(ix,j,k)%BCell%flisteg_nr(1:8)=fliste_egnr(1:8)
                bfyznr=UpBoundsLayer(ix,j,k)%BCell%bfyznr
                wbfyznr=wbfyznr+1
                !tropo_id/bound_id noch offen
                CALL WOutBFace(wbfyznr,ix,j,k,"fyz", bfyznr, fegnr, 0, fliste_egnr(1:8))
              END IF          
            END DO
          END DO
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   ! [o,p,i]- 'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar
END SUBROUTINE WriteBoundBorderWE_FacesYZ


SUBROUTINE WriteBoundBorderWE_Cells
   INTEGER :: in,ix,jx,j,k,liType,li_max
   INTEGER :: wbncpnr,ifn
   wbncpnr=0
   li_max=0
   !...%BCell%ub_cnr
   DO in=1,AnzahlNachbar
      CALL Set(Nachbars(in))
      SELECT CASE (Nachbars(in)%nType(2:2))
        !----------------------------------
      CASE ("w","e")
          !z.Zt. noch gleiches Refine wenn iw || ie
          IF(nType(2:2)=='w') THEN
              ix=ix0
              jx=ix1
              wbncpnr=0
          ELSE !'e' 
              ix=ix1+1
              jx=ix0+1
              wbncpnr=0
          END IF
          !...%BCell%ub_cnr
          DO j=iy0+1,iy1
            DO k=iz0+1,iz1
              IF(UpBoundsLayer(ix,j,k)%BCell%ub_cnr>0) THEN
                !Face-Bezeichner: bfxynr,bfzxnr,bfyznr,FCutNr
                !Reihenfolge: f1-xy,f6-yz,f4-zx,f5-yz,f3-zx,f2-xy,fc
                ifn=0
                IF(UpBoundsLayer(ix,j,k-1)%BCell%bfxynr>0) THEN
                   ifn=ifn+1
                   UpBoundsLayer(ix,j,k)%BCell%cfn(ifn)= &
                       & UpBoundsLayer(ix,j,k-1)%BCell%bfxynr
                END IF
                IF(UpBoundsLayer(ix,j,k)%BCell%bfyznr>0) THEN
                   ifn=ifn+1
                   UpBoundsLayer(ix,j,k)%BCell%cfn(ifn)= &
                       & UpBoundsLayer(ix,j,k)%BCell%bfyznr
                END IF
                IF(UpBoundsLayer(ix,j,k)%BCell%bfzxnr>0) THEN
                   ifn=ifn+1
                   UpBoundsLayer(ix,j,k)%BCell%cfn(ifn)= &
                       & UpBoundsLayer(ix,j,k)%BCell%bfzxnr
                END IF
                IF(UpBoundsLayer(ix-1,j,k)%BCell%bfyznr>0) THEN
                   ifn=ifn+1
                   UpBoundsLayer(ix,j,k)%BCell%cfn(ifn)= &
                       & UpBoundsLayer(ix-1,j,k)%BCell%bfyznr
                END IF
                IF(UpBoundsLayer(ix,j-1,k)%BCell%bfzxnr>0) THEN
                   ifn=ifn+1
                   UpBoundsLayer(ix,j,k)%BCell%cfn(ifn)= &
                       & UpBoundsLayer(ix,j-1,k)%BCell%bfzxnr
                END IF
                IF(UpBoundsLayer(ix,j,k)%BCell%bfxynr>0) THEN
                   ifn=ifn+1
                   UpBoundsLayer(ix,j,k)%BCell%cfn(ifn)= &
                       & UpBoundsLayer(ix,j,k)%BCell%bfxynr
                END IF
                IF(UpBoundsLayer(ix,j,k)%BCell%FCutNr>0) THEN
                   ifn=ifn+1
                   UpBoundsLayer(ix,j,k)%BCell%cfn(ifn)= &
                       & UpBoundsLayer(ix,j,k)%BCell%FCutNr
                END IF
                UpBoundsLayer(ix,j,k)%BCell%nfn=ifn
                wbncpnr=wbncpnr+1
                CALL WOutBCell(wbncpnr,ix,j,k,"cp",UpBoundsLayer(ix,j,k)%BCell%ub_cnr, &
                      &  ifn, tropos_border, UpBoundsLayer(ix,j,k)%BCell%cfn)
              ELSE
                 UpBoundsLayer(ix,j,k)%BCell%nfn=0
              END IF
            END DO
          END DO
      CASE ("n","s")
      !----------------------------------
          IF (nType=="on") THEN
          ELSE IF (nType=="os") THEN
          ELSE
          END IF
      CASE ("t","b")
      !----------------------------------
          IF (nType=="ot") THEN
          ELSE IF (nType=="ob") THEN
          ELSE
          END IF
      END SELECT   ! [o,p,i]- 'w/e','n/s','t/b' -> to boundsborderlayer
   END DO  ! AnzahlNachbar

END SUBROUTINE WriteBoundBorderWE_Cells


SUBROUTINE WriteUpperBoundsLayer(FileName)
  CHARACTER*50 :: FileName
  ! local...
  INTEGER :: ib,i,j,k,n,ic,iec,is_eg,in_out
  INTEGER :: becnr,ec_set,pos_ec1,pos_ecn,ifn,lvertout
  INTEGER :: eg_listfxy(1:8)=0
  INTEGER :: eg_listfzx(1:8)=0
  INTEGER :: eg_listfyz(1:8)=0
  INTEGER, DIMENSION(:,:), ALLOCATABLE :: bliste_ec
  INTEGER :: wbnex,wbney,wbnez,wbnec
  INTEGER :: wbnfxy,wbnfzx,wbnfyz,wbnfc
  INTEGER :: wbncnr,wbncvcnr
 
  CALL WriteAuswUpperBoundsLayerProt()

  CALL SortIndexOutUpBounds
  CALL Display_OutUpperBoundsLayer(FileName)
  !WRITE(*,*) 'Write File mit Struktur Upper-Bounds-Layer'
  !WRITE(*,*) ">....... Upper Bounds Layer for ASAM-Input"
  CALL OpenFileForBound(FileName)
  CALL WOutNrBlkGrid(nb)
  CALL WOutVertBounds

  DO ib=1,nb
    wbnex=0;wbney=0;wbnez=0;wbnec=0
    wbnfxy=0;wbnfyz=0;wbnfzx=0;wbnfc=0
    wbncnr=0;wbncvcnr=0
    Write(*,*) "                                     Block: ",ib,"\/",nb
    CALL Set(Floor(ib))

    Write(OutUnitB,*)"............................................"
    Write(OutUnitB,*)ib, "   ! block: ",ib,"/",nb
    Write(OutUnitB,*)"............................................"
    Write(OutUnitBpv,*)"............................................"
    Write(OutUnitBpv,*)ib, "   ! block: ",ib,"/",nb
    Write(OutUnitBpv,*)"............................................"

    IF(stopo_ecnr>0) THEN
       ALLOCATE(bliste_ec(upb_pos_ec1:upb_pos_ecn,1:6))
     END IF

    !-- Output---/EdgeX/EdgeY/EdgeZ/EdgeCut
    CALL WOutNrBEdge(sbenr,sbexnr,sbeynr,sbeznr,sbecnr)
    !-- Output--Edge_X
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
        !DO k=DimUpBoundsArray(i,j,1),DimUpBoundsArray(i,j,2)
          IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
             IF(MAX(Edges_X(i,j,k)%Edge%Vert1%in_out,Edges_X(i,j,k)%Edge%Vert2%in_out)>0 &
               .OR. &
               (Edges_X(i,j,k)%Edge%Vert1%in_out==Edges_X(i,j,k)%Edge%Vert2%in_out.AND. &
                Edges_X(i,j,k)%Edge%Vert1%in_out==0)) THEN
               IF(Edges_X(i,j,k)%Edge%yes_sp==1) THEN
                  IF(Edges_X(i,j,k)%Edge%Vert1%nrP==0) THEN
                    wbnex=wbnex+1
                    CALL WOutBEdge(wbnex,i,j,k,"ex", UpBoundsLayer(i,j,k)%BCell%bexnr, &
                           & NrPOutBounds(Edges_X(i,j,k)%Edge%VertS%nrP), &
                           & NrPOutBounds(Edges_X(i,j,k)%Edge%Vert2%nrP)) 
                  ELSE IF(Edges_X(i,j,k)%Edge%Vert2%nrP==0) THEN
                    wbnex=wbnex+1
                    CALL WOutBEdge(wbnex,i,j,k, "ex", UpBoundsLayer(i,j,k)%BCell%bexnr, &
                           & NrPOutBounds(Edges_X(i,j,k)%Edge%Vert1%nrP), &
                           & NrPOutBounds(Edges_X(i,j,k)%Edge%VertS%nrP))
                  END IF
               ELSE
                 wbnex=wbnex+1
                 CALL WOutBEdge(wbnex,i,j,k,"ex", UpBoundsLayer(i,j,k)%BCell%bexnr, &
                        & NrPOutBounds(Edges_X(i,j,k)%Edge%Vert1%nrP), &
                        & NrPOutBounds(Edges_X(i,j,k)%Edge%Vert2%nrP))
               END IF
             END IF
          ELSE
            IF(MAX(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)>0 &
               .OR. &
               (Vertices(i-1,j,k)%in_out==Vertices(i,j,k)%in_out .AND. &
                Vertices(i,j,k)%in_out==0)) THEN
               IF(UpBoundsLayer(i,j,k)%BCell%bexnr>0) THEN
                wbnex=wbnex+1
                CALL WOutBEdge(wbnex,i,j,k,"ex", UpBoundsLayer(i,j,k)%BCell%bexnr, &
                       & NrPOutBounds(Vertices(i-1,j,k)%nrP), &
                       & NrPOutBounds(Vertices(i,j,k)%nrP))
               END IF
            END IF 
          END IF
        END DO
      END DO
    END DO
    CALL WriteBoundBorderWE_EdgesX()
    !-- Output--Edge_Y
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
        !DO k=DimUpBoundsArray(i,j,1),DimUpBoundsArray(i,j,2)
          !IF(UpBoundsLayer(i,j,k)%BCell%beynr>0) THEN
          !      UpBoundsLayer(i,j,k)%BCell%beynr= &
          !            & UpBoundsLayer(i,j,k)%BCell%beynr+sbexnr
          !END IF
          IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
            IF(MAX(Edges_Y(i,j,k)%Edge%Vert1%in_out,Edges_Y(i,j,k)%Edge%Vert2%in_out)>0 &
               .OR. &
               (Edges_Y(i,j,k)%Edge%Vert1%in_out==Edges_Y(i,j,k)%Edge%Vert2%in_out.AND. &
                Edges_Y(i,j,k)%Edge%Vert1%in_out==0)) THEN
               IF(Edges_Y(i,j,k)%Edge%yes_sp==1) THEN
                 IF(Edges_Y(i,j,k)%Edge%Vert1%nrP==0) THEN
                   wbney=wbney+1
                   CALL WOutBEdge(wbney,i,j,k,"ey", UpBoundsLayer(i,j,k)%BCell%beynr,&
                          &  NrPOutBounds(Edges_Y(i,j,k)%Edge%VertS%nrP), &
                          &  NrPOutBounds(Edges_Y(i,j,k)%Edge%Vert2%nrP))
                 ELSE IF(Edges_Y(i,j,k)%Edge%Vert2%nrP==0) THEN
                   wbney=wbney+1
                   CALL WOutBEdge(wbney,i,j,k,"ey", UpBoundsLayer(i,j,k)%BCell%beynr, &
                          & NrPOutBounds(Edges_Y(i,j,k)%Edge%Vert1%nrP), &
                          & NrPOutBounds(Edges_Y(i,j,k)%Edge%VertS%nrP) )
                 END IF
               ELSE
                  wbney=wbney+1
                  CALL WOutBEdge(wbney,i,j,k,"ey", UpBoundsLayer(i,j,k)%BCell%beynr, &
                         & NrPOutBounds(Edges_Y(i,j,k)%Edge%Vert1%nrP), &
                         & NrPOutBounds(Edges_Y(i,j,k)%Edge%Vert2%nrP))
               END IF
            END IF
          ELSE
            IF(MAX(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)>0 &
               .OR. &
               (Vertices(i,j-1,k)%in_out==Vertices(i,j,k)%in_out .AND. &
                Vertices(i,j,k)%in_out==0)) THEN
              IF(UpBoundsLayer(i,j,k)%BCell%beynr>0) THEN
                wbney=wbney+1
                CALL WOutBEdge(wbney,i,j,k,"ey", UpBoundsLayer(i,j,k)%BCell%beynr, &
                       & NrPOutBounds(Vertices(i,j-1,k)%nrP), &
                       & NrPOutBounds(Vertices(i,j,k)%nrP))
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    CALL WriteBoundBorderWE_EdgesY()
    !-- Output--Edge_Z
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
        !DO k=DimUpBoundsArray(i,j,1)+1,DimUpBoundsArray(i,j,2)
          !IF(UpBoundsLayer(i,j,k)%BCell%beznr>0) THEN
          !    UpBoundsLayer(i,j,k)%BCell%beznr= &
          !       UpBoundsLayer(i,j,k)%BCell%beznr+sbexnr+sbeynr
          !END IF
          IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
            IF(MAX(Edges_Z(i,j,k)%Edge%Vert1%in_out,Edges_Z(i,j,k)%Edge%Vert2%in_out)>0 &
               .OR. &
               (Edges_Z(i,j,k)%Edge%Vert1%in_out==Edges_Z(i,j,k)%Edge%Vert2%in_out.AND. &
                Edges_Z(i,j,k)%Edge%Vert1%in_out==0)) THEN
               IF(Edges_Z(i,j,k)%Edge%yes_sp==1) THEN
                 IF(Edges_Z(i,j,k)%Edge%Vert1%nrP==0) THEN
                   wbnez=wbnez+1
                   CALL WOutBEdge(wbnez,i,j,k,"ez", UpBoundsLayer(i,j,k)%BCell%beznr, &
                          &  NrPOutBounds(Edges_Z(i,j,k)%Edge%VertS%nrP), &
                          &  NrPOutBounds(Edges_Z(i,j,k)%Edge%Vert2%nrP))
                 ELSE IF(Edges_Z(i,j,k)%Edge%Vert2%nrP==0) THEN
                   wbnez=wbnez+1
                   CALL WOutBEdge(wbnez,i,j,k,"ez", UpBoundsLayer(i,j,k)%BCell%beznr, &
                          & NrPOutBounds(Edges_Z(i,j,k)%Edge%Vert1%nrP), &
                          & NrPOutBounds(Edges_Z(i,j,k)%Edge%VertS%nrP))
                 END IF
               ELSE
                 wbnez=wbnez+1
                 CALL WOutBEdge(wbnez,i,j,k,"ez", UpBoundsLayer(i,j,k)%BCell%beznr, &
                        & NrPOutBounds(Edges_Z(i,j,k)%Edge%Vert1%nrP), &
                        & NrPOutBounds(Edges_Z(i,j,k)%Edge%Vert2%nrP))
               END IF
            END IF
          ELSE
            IF(MAX(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)>0 &
               .OR. &
               (Vertices(i,j,k-1)%in_out==Vertices(i,j,k)%in_out.AND. &
                Vertices(i,j,k)%in_out==0)) THEN
              IF(UpBoundsLayer(i,j,k)%BCell%beznr>0) THEN
               wbnez=wbnez+1
               CALL WOutBEdge(wbnez,i,j,k,"ez", UpBoundsLayer(i,j,k)%BCell%beznr, &
                      & NrPOutBounds(Vertices(i,j,k-1)%nrP), &
                      & NrPOutBounds(Vertices(i,j,k)%nrP))
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    CALL WriteBoundBorderWE_EdgesZ()
    !-- EdgeCut-Sortieren-Ausgabe
    becnr=0
    IF(topo_ecnr>0) THEN
      DO i=ix0+1,ix1     !EdgeCut-Sortieren-für Ausgabe
        DO j=iy0+1,iy1
          DO k=iz0,iz1
            IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
              IF(Faces_XY(i,j,k)%Face%egc_nr>0) THEN
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrxy,1)= &
                                      &  UpBoundsLayer(i,j,k)%BCell%becnrxy
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrxy,2)= &
                               & NrPOutBounds(Faces_XY(i,j,k)%Face%EdgeCut(1))
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrxy,3)= &
                               & NrPOutBounds(Faces_XY(i,j,k)%Face%EdgeCut(2))
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrxy,4)= i
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrxy,5)= j
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrxy,6)= k
                 becnr=becnr+1
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
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnryz,1)= &
                                     & UpBoundsLayer(i,j,k)%BCell%becnryz
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnryz,2)= &
                          & NrPOutBounds(Faces_YZ(i,j,k)%Face%EdgeCut(1))
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnryz,3)= &
                          & NrPOutBounds(Faces_YZ(i,j,k)%Face%EdgeCut(2))
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnryz,4)= i
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnryz,5)= j
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnryz,6)= k
                 becnr=becnr+1
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
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrzx,1)= &
                                    & UpBoundsLayer(i,j,k)%BCell%becnrzx
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrzx,2)= &
                         & NrPOutBounds(Faces_ZX(i,j,k)%Face%EdgeCut(1))
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrzx,3)= &
                         & NrPOutBounds(Faces_ZX(i,j,k)%Face%EdgeCut(2))
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrzx,4)= i
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrzx,5)= j
                 bliste_ec(UpBoundsLayer(i,j,k)%BCell%becnrzx,6)= k
                 becnr=becnr+1
              END IF
            END IF
          END DO
        END DO
      END DO
      !-- Output--EdgeCut

      CALL WOutBEdgeCut(upb_pos_ec1,upb_pos_ecn,bliste_ec)
      DEALLOCATE(bliste_ec)
    END IF
    CALL WriteBoundBorderWE_EdgesEc()

    !-- Output--FaceCut,Faces_XY,Faces_XZ,Faces_YZ ---------------------------------------
    !-- Output--FaceCut
    CALL WOutNrBFace(sbfnr,sbfcnr,sbfxynr,sbfzxnr,sbfyznr)
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
         IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
            !IF(Cell(i,j,k)%Cell%vc>0) THEN
            IF(Topo(i,j,k)%Topo%FCutNr>0) THEN
               wbnfc=wbnfc+1
               CALL WOutBFace(wbnfc,i,j,k,"fc ",Topo(i,j,k)%Topo%FCutNr, Cell(i,j,k)%Cell%vc, &
                      &  tropos_border, UpBoundsLayer(i,j,k)%BCell%fceg_nrcut )
            END IF
          END IF
        END DO
      END DO
    END DO
    CALL WriteBoundBorderWE_FacesFc
    !-- Output--Faces_XY
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0,iz1
        !DO k=DimUpBoundsArray(i,j,1),DimUpBoundsArray(i,j,2)
          IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
            IF(Vertices(i-1,j-1,k)%in_out==0 .and. &
               Vertices(i  ,j-1,k)%in_out==0 .and. &
               Vertices(i-1,j  ,k)%in_out==0 .and. &
               Vertices(i  ,j  ,k)%in_out==0) THEN
               IF(k==iz0) THEN
                 ! wenn Boden-Modell
                 !Face1-Grenzfläche ( Hill.grid, Cellen 23/1/1-79/1/1 bei iz0 Grenze)
                 !(Cell(i,j,k)%Cell%vc==0.and.Cell(i,j,k)%Cell%in_out==4)
                 ! kein F-Cut!
                 eg_listfxy(1)=UpBoundsLayer(i,j-1,k)%BCell%bexnr
                 eg_listfxy(2)=UpBoundsLayer(i,j  ,k)%BCell%beynr
                 eg_listfxy(3)=UpBoundsLayer(i,j  ,k)%BCell%bexnr
                 eg_listfxy(4)=UpBoundsLayer(i-1,j,k)%BCell%beynr
                   wbnfxy=wbnfxy+1
                   CALL WOutBFace(wbnfxy,i,j,k,"fxy", UpBoundsLayer(i,j,k)%BCell%bfxynr, 4, &
                            &  tropos_border,eg_listfxy) 
               ELSE
                 ! als Grenzfläche 'fc' ausgegeben
                 !Write(OutUnitBpv,*) "Faces_XY()%Face als Grenzfläche fc ausgegeben &
                 !         &  (4xVert%in_out==0)",i,j,k
                 !Write(OutUnitBpv,*) "               FCutNr:", Topo(i,j,k)%Topo%FCutNr
                 ! --> Special suchen , betrifft Berg <----
               END IF
            ELSE IF(Faces_XY(i,j,k)%Face%in_out>-4.and.Faces_XY(i,j,k)%Face%numberVert>2) THEN
               !Reihenfolge:e1-e4, speciell je Face-Def.
               is_eg=0;ec_set=0
               IF(Faces_XY(i,j,k)%Face%ec==1) THEN
                  ! EdgeListe aufstellen
                 !.e1.....................................
                 IF(Topo(i,j-1,k)%Topo%exnr>0) THEN
                   is_eg=is_eg+1
                   eg_listfxy(is_eg)=UpBoundsLayer(i,j-1,k)%BCell%bexnr !Topo(i,j-1,k)%Topo%exnr
                   IF(Faces_XY(i,j,k)%Face%Edge1%yes_sp==1) THEN
                      IF(Faces_XY(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                         is_eg=is_eg+1
                         eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                         ec_set=1
                         !Write(OutUnitBpv,*) i,j,k,", eg1 - Annulus Test F_XY"
                      END IF
                   END IF
                 ELSE 
                   IF(Faces_XY(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                      Faces_XY(i,j,k)%Face%Edge1%in_out==-1) THEN
                       is_eg=is_eg+1
                       eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                       ec_set=1
                   END IF
                 END IF
                 !.e2.....................................
                 IF(Topo(i,j,k)%Topo%eynr>0) THEN
                    is_eg=is_eg+1
                    eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%beynr !Topo(i,j,k)%Topo%eynr
                    IF(Faces_XY(i,j,k)%Face%Edge2%yes_sp==1) THEN
                      IF(Faces_XY(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                          ec_set=1
                          !Write(OutUnitBpv,*) i,j,k,",  eg2 - Annulus Test F_XY"
                        END IF
                      END IF
                    END IF
                 ELSE 
                    IF(Faces_XY(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                       Faces_XY(i,j,k)%Face%Edge2%in_out==-1) THEN
                      IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                        ec_set=1
                      END IF
                    END IF
                 END IF
                 !.e3.....................................
                 IF(Topo(i,j,k)%Topo%exnr>0) THEN
                   is_eg=is_eg+1
                   eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%bexnr !Topo(i,j,k)%Topo%exnr
                   IF(Faces_XY(i,j,k)%Face%Edge3%yes_sp==1) THEN
                     IF(Faces_XY(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                       IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                        ec_set=1
                        !Write(OutUnitBpv,*) i,j,k,",  eg3 - Annulus Test F_XY"
                       END IF
                     END IF
                   END IF
                 ELSE
                   IF(Faces_XY(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                      Faces_XY(i,j,k)%Face%Edge3%in_out==-1) THEN
                     IF(ec_set==0) THEN
                       is_eg=is_eg+1
                       eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                       ec_set=1
                     END IF
                   END IF
                 END IF
                 !.e4.....................................
                 IF(Topo(i-1,j,k)%Topo%eynr>0) THEN
                   is_eg=is_eg+1
                   eg_listfxy(is_eg)=UpBoundsLayer(i-1,j,k)%BCell%beynr  !Topo(i-1,j,k)%Topo%eynr
                   IF(Faces_XY(i,j,k)%Face%Edge4%yes_sp==1) THEN
                     IF(Faces_XY(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                       IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                        ! Write(OutUnitBpv,*) i,j,k,",  eg4 - Annulus Test F_XY"
                       END IF
                     END IF
                   END IF 
                 ELSE 
                    IF(Faces_XY(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                       Faces_XY(i,j,k)%Face%Edge4%in_out==-1) THEN
                      IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                      END IF
                    END IF
                 END IF
                 IF(is_eg==2 .and. Faces_XY(i,j,k)%Face%NumberVert>2) THEN
                    is_eg=is_eg+1
                    eg_listfxy(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrxy
                 END IF
                 !........................................
                 IF(k==iz0.OR.k==iz1) THEN
                    wbnfxy=wbnfxy+1
                    CALL WOutBFace(wbnfxy,i,j,k,"fxy", UpBoundsLayer(i,j,k)%BCell%bfxynr, is_eg, &
                            &   tropos_border,eg_listfxy)
                 ELSE 
                    wbnfxy=wbnfxy+1
                    CALL WOutBFace(wbnfxy,i,j,k,"fxy", UpBoundsLayer(i,j,k)%BCell%bfxynr, is_eg, &
                            &   tropos_inside,eg_listfxy)
                 END IF
               ELSE 
                 eg_listfxy(1)=UpBoundsLayer(i,j-1,k)%BCell%bexnr
                 eg_listfxy(2)=UpBoundsLayer(i,j  ,k)%BCell%beynr
                 eg_listfxy(3)=UpBoundsLayer(i,j  ,k)%BCell%bexnr
                 eg_listfxy(4)=UpBoundsLayer(i-1,j,k)%BCell%beynr
                 IF(k==iz0.OR.k==iz1) THEN
                    wbnfxy=wbnfxy+1
                   CALL WOutBFace(wbnfxy,i,j,k,"fxy", UpBoundsLayer(i,j,k)%BCell%bfxynr, 4, &
                            &  tropos_border,eg_listfxy) 
                 ELSE 
                    wbnfxy=wbnfxy+1
                   CALL WOutBFace(wbnfxy,i,j,k,"fxy", UpBoundsLayer(i,j,k)%BCell%bfxynr, 4, &
                            &  tropos_inside,eg_listfxy)
                 END IF
               END IF
            END IF
          ELSE
            !IF(k<=DimUpBoundsArray(i-1,j,2)) THEN
            IF(UpBoundsLayer(i,j,k)%BCell%bfxynr>0) THEN
              in_out=Vertices(i-1,j-1,k)%in_out+Vertices(i,j-1,k)%in_out+ &
                     Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
              IF(in_out>=0) THEN
                 eg_listfxy(1)=UpBoundsLayer(i,j-1,k)%BCell%bexnr
                 eg_listfxy(2)=UpBoundsLayer(i,j  ,k)%BCell%beynr
                 eg_listfxy(3)=UpBoundsLayer(i,j  ,k)%BCell%bexnr
                 eg_listfxy(4)=UpBoundsLayer(i-1,j,k)%BCell%beynr
                 IF(k==iz0.OR.k==iz1) THEN
                    wbnfxy=wbnfxy+1
                    CALL WOutBFace(wbnfxy,i,j,k,"fxy", UpBoundsLayer(i,j,k)%BCell%bfxynr, 4, & 
                           &   tropos_border,eg_listfxy) 
                 ELSE 
                    wbnfxy=wbnfxy+1
                    CALL WOutBFace(wbnfxy,i,j,k,"fxy", UpBoundsLayer(i,j,k)%BCell%bfxynr, 4, & 
                           &   tropos_inside,eg_listfxy) 
                 END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    CALL WriteBoundBorderWE_FacesXY
    !-- Output--Faces_ZX
    DO i=ix0+1,ix1
      DO j=iy0,iy1
        DO k=iz0+1,iz1
        !DO k=DimUpBoundsArray(i,j,1)+1,DimUpBoundsArray(i,j,2)
          IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
            IF(Vertices(i-1,j,k-1)%in_out==0 .and. &
               Vertices(i  ,j,k-1)%in_out==0 .and. & 
               Vertices(i-1,j,k  )%in_out==0 .and. &
               Vertices(i  ,j,k  )%in_out==0 ) THEN
               !!als  Grenzfläche  'fc' ausgegeben 
               !Write(OutUnitBpv,*) "Faces_ZX(i,j,k)%Face als Grenzfläche fc ausgegeben, &
               !           &  4xVert%in_out==0",i,j,k
               !Write(OutUnitBpv,*) "                  FCutNr:", Topo(i,j,k)%Topo%FCutNr
            ELSE IF(Faces_ZX(i,j,k)%Face%in_out>-4.and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
               !Reihenfolge:e1-e4, speciell je Face-Def.
               is_eg=0;ec_set=0
               IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
                 ! EdgeListe aufstellen
                 !.e1.....................................
                 IF(Topo(i-1,j,k)%Topo%eznr>0) THEN
                   is_eg=is_eg+1
                   eg_listfzx(is_eg)=UpBoundsLayer(i-1,j,k)%BCell%beznr !Topo(i-1,j,k)%Topo%eznr
                   IF(Faces_ZX(i,j,k)%Face%Edge1%yes_sp==1) THEN
                     IF(Faces_ZX(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                        is_eg=is_eg+1
                        eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                        ec_set=1
                        !Write(OutUnitBpv,*) i,j,k,", eg1 - Annulus Test F_ZX"
                     END IF
                   END IF
                 ELSE 
                   IF(Faces_ZX(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                      Faces_ZX(i,j,k)%Face%Edge1%in_out==-1) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                     ec_set=1
                   END IF
                 END IF
                 !.e2.....................................
                 IF(Topo(i,j,k)%Topo%exnr>0) THEN
                    is_eg=is_eg+1
                    eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%bexnr !Topo(i,j,k)%Topo%exnr
                    IF(Faces_ZX(i,j,k)%Face%Edge2%yes_sp==1) THEN
                       IF(Faces_ZX(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                         IF(ec_set==0) THEN
                           is_eg=is_eg+1
                           eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                           ec_set=1
                           !Write(OutUnitBpv,*) i,j,k,", eg2 - Annulus Test F_ZX"
                         END IF
                     END IF
                   END IF
                 ELSE 
                    IF(Faces_ZX(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                       Faces_ZX(i,j,k)%Face%Edge2%in_out==-1) THEN
                      IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                        ec_set=1
                      END IF
                    END IF
                 END IF
                 !.e3.....................................
                 IF(Topo(i,j,k)%Topo%eznr>0) THEN
                   is_eg=is_eg+1
                   eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%beznr  !Topo(i,j,k)%Topo%eznr
                   IF(Faces_ZX(i,j,k)%Face%Edge3%yes_sp==1) THEN
                      IF(Faces_ZX(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                          ec_set=1
                          !Write(OutUnitBpv,*) i,j,k,", eg3 - Annulus Test F_ZX"
                        END IF
                      END IF
                   END IF
                 ELSE
                   IF(Faces_ZX(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                      Faces_ZX(i,j,k)%Face%Edge3%in_out==-1) THEN
                     IF(ec_set==0) THEN
                       is_eg=is_eg+1
                       eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                       ec_set=1
                     END IF
                   END IF
                 END IF
                 !.e4.....................................
                 IF(Topo(i,j,k-1)%Topo%exnr>0) THEN
                   is_eg=is_eg+1
                   eg_listfzx(is_eg)=UpBoundsLayer(i,j,k-1)%BCell%bexnr  !Topo(i,j,k-1)%Topo%exnr
                   IF(Faces_ZX(i,j,k)%Face%Edge4%yes_sp==1) THEN
                      IF(Faces_ZX(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                          !Write(OutUnitBpv,*) i,j,k,", eg4 - Annulus Test F_ZX"
                        END IF
                      END IF
                   END IF
                 ELSE 
                    IF(Faces_ZX(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                       Faces_ZX(i,j,k)%Face%Edge4%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                       END IF
                    END IF
                 END IF
                 IF(is_eg==2 .and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
                     is_eg=is_eg+1
                     eg_listfzx(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnrzx
                 END IF
                 !........................................
                 IF(j==iy0.OR.j==iy1) THEN
                    wbnfzx=wbnfzx+1
                    CALL WOutBFace(wbnfzx,i,j,k,"fzx", UpBoundsLayer(i,j,k)%BCell%bfzxnr, is_eg, &
                           &   tropos_border, eg_listfzx) 
                 ELSE 
                    CALL WOutBFace(wbnfzx,i,j,k,"fzx", UpBoundsLayer(i,j,k)%BCell%bfzxnr, is_eg, &
                           &   tropos_inside,eg_listfzx) 
                 END IF
               ELSE ! not. Faces_ZX(i,j,k)%Face%ec 
                 eg_listfzx(1)=UpBoundsLayer(i-1,j,k)%BCell%beznr
                 eg_listfzx(2)=UpBoundsLayer(i  ,j,k)%BCell%bexnr
                 eg_listfzx(3)=UpBoundsLayer(i  ,j,k)%BCell%beznr
                 eg_listfzx(4)=UpBoundsLayer(i,j,k-1)%BCell%bexnr
                 IF(j==iy0.OR.j==iy1) THEN
                    wbnfzx=wbnfzx+1
                   CALL WOutBFace(wbnfzx,i,j,k,"fzx", UpBoundsLayer(i,j,k)%BCell%bfzxnr, 4, & 
                           &   tropos_border,eg_listfzx) 
                 ELSE 
                    wbnfzx=wbnfzx+1
                   CALL WOutBFace(wbnfzx,i,j,k,"fzx", UpBoundsLayer(i,j,k)%BCell%bfzxnr, 4, & 
                           &   tropos_inside,eg_listfzx) 
                 END IF
               END IF !IF/Else Faces_ZX(i,j,k)%Face%ec
            END IF
          ELSE  ! .not. ASSOCIATED(Faces_ZX(i,j,k)%Face
            !IF(k<=DimUpBoundsArray(i-1,j,2)) THEN
            IF(UpBoundsLayer(i,j,k)%BCell%bfzxnr>0) THEN
              in_out=Vertices(i-1,j,k-1)%in_out+Vertices(i,j,k-1)%in_out+ &
                     Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
              IF(in_out>=0) THEN
                 eg_listfzx(1)=UpBoundsLayer(i-1,j,k)%BCell%beznr
                 eg_listfzx(2)=UpBoundsLayer(i  ,j,k)%BCell%bexnr
                 eg_listfzx(3)=UpBoundsLayer(i  ,j,k)%BCell%beznr
                 eg_listfzx(4)=UpBoundsLayer(i,j,k-1)%BCell%bexnr
                 IF(j==iy0.OR.j==iy1) THEN
                    wbnfzx=wbnfzx+1
                    CALL WOutBFace(wbnfzx,i,j,k,"fzx", UpBoundsLayer(i,j,k)%BCell%bfzxnr, 4, &
                             &   tropos_border,eg_listfzx) 
                 ELSE 
                    wbnfzx=wbnfzx+1
                    CALL WOutBFace(wbnfzx,i,j,k,"fzx", UpBoundsLayer(i,j,k)%BCell%bfzxnr, 4, &
                             &   tropos_inside,eg_listfzx) 
                 END IF
              END IF
            END IF
          END IF  ! IF/Else ASSOCIATED(Faces_ZX(i,j,k)%Face
        END DO  ! k
      END DO   ! j
    END DO    ! i
    CALL WriteBoundBorderWE_FacesZX
    !-- Output--Faces_YZ
    DO i=ix0,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
        !DO k=DimUpBoundsArray(i,j,1)+1,DimUpBoundsArray(i,j,2)
          IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
            IF(Vertices(i,j-1,k-1)%in_out==0 .and. &
               Vertices(i,j-1,k  )%in_out==0 .and. &
               Vertices(i,j  ,k  )%in_out==0 .and. &
               Vertices(i,j  ,k-1)%in_out==0 ) THEN
               ! als Grenzfläche 'fc' ausgegeben 
               !Write(OutUnitBpv,*) "Faces_YZ(i,j,k)%Face als Grenzfläche fc ausgegeben, &
               !          &  4xVert%in_out==0",i,j,k
               !Write(OutUnitBpv,*) "                       FCutNr:", Topo(i,j,k)%Topo%FCutNr
            ELSE IF(Faces_YZ(i,j,k)%Face%in_out>-4.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
               !Reihenfolge:e1-e4, speciell je Face-Def.
               is_eg=0;ec_set=0
               IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
                 ! EdgeListe aufstellen
                 !.e1.....................................
                 IF(Topo(i,j,k-1)%Topo%eynr>0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=UpBoundsLayer(i,j,k-1)%BCell%beynr !Topo(i,j,k-1)%Topo%eynr
                   IF(Faces_YZ(i,j,k)%Face%Edge1%yes_sp==1) THEN
                     IF(Faces_YZ(i,j,k)%Face%Edge1%Vert2%in_out==-1) THEN
                       is_eg=is_eg+1
                       eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnryz
                       ec_set=1
                       !Write(OutUnitBpv,*) i,j,k,",  eg1 - Annulus Test F_YZ"
                     END IF
                   END IF
                 ELSE 
                   IF(Faces_YZ(i,j,k)%Face%Edge1%yes_sp==1 .OR. &
                      Faces_YZ(i,j,k)%Face%Edge1%in_out==-1) THEN
                     is_eg=is_eg+1
                     eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnryz
                     ec_set=1
                   END IF
                 END IF
                 !.e2.....................................
                 IF(Topo(i,j,k)%Topo%eznr>0) THEN
                    is_eg=is_eg+1
                    eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%beznr  !Topo(i,j,k)%Topo%eznr
                    IF(Faces_YZ(i,j,k)%Face%Edge2%yes_sp==1) THEN
                      IF(Faces_YZ(i,j,k)%Face%Edge2%Vert2%in_out==-1) THEN
                        IF(ec_set==0) THEN
                          is_eg=is_eg+1
                          eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnryz
                          ec_set=1
                          !Write(OutUnitBpv,*) i,j,k,",  eg2 - Annulus Test F_YZ"
                        END IF
                      END IF
                    END IF
                 ELSE 
                    IF(Faces_YZ(i,j,k)%Face%Edge2%yes_sp==1 .OR. &
                       Faces_YZ(i,j,k)%Face%Edge2%in_out==-1) THEN
                      IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnryz
                        ec_set=1
                      END IF
                    END IF
                 END IF
                 !.e3.....................................
                 IF(Topo(i,j,k)%Topo%eynr>0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%beynr  !Topo(i,j,k)%Topo%eynr
                   IF(Faces_YZ(i,j,k)%Face%Edge3%yes_sp==1) THEN
                     IF(Faces_YZ(i,j,k)%Face%Edge3%Vert1%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnryz
                         ec_set=1
                         !Write(OutUnitBpv,*) i,j,k,",  eg3 - Annulus Test F_YZ"
                       END IF
                     END IF
                   END IF
                 ELSE
                   IF(Faces_YZ(i,j,k)%Face%Edge3%yes_sp==1 .OR. &
                      Faces_YZ(i,j,k)%Face%Edge3%in_out==-1) THEN
                     IF(ec_set==0) THEN
                       is_eg=is_eg+1
                       eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnryz
                       ec_set=1
                     END IF
                   END IF
                 END IF
                 !.e4.....................................
                 IF(Topo(i,j-1,k)%Topo%eznr>0) THEN
                   is_eg=is_eg+1
                   eg_listfyz(is_eg)=UpBoundsLayer(i,j-1,k)%BCell%beznr  !Topo(i,j-1,k)%Topo%eznr
                   IF(Faces_YZ(i,j,k)%Face%Edge4%yes_sp==1) THEN
                     IF(Faces_YZ(i,j,k)%Face%Edge4%Vert1%in_out==-1) THEN
                       IF(ec_set==0) THEN
                         is_eg=is_eg+1
                         eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnryz
                         !Write(OutUnitBpv,*) i,j,k,",  eg4 - Annulus Test F_YZ"
                       END IF
                     END IF
                   END IF
                 ELSE 
                    IF(Faces_YZ(i,j,k)%Face%Edge4%yes_sp==1 .OR. &
                       Faces_YZ(i,j,k)%Face%Edge4%in_out==-1) THEN
                      IF(ec_set==0) THEN
                        is_eg=is_eg+1
                        eg_listfyz(is_eg)=UpBoundsLayer(i,j,k)%BCell%becnryz
                      END IF
                    END IF
                 END IF
                 IF(is_eg==2.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
                      is_eg=is_eg+1
                      eg_listfyz(is_eg)=Topo(i,j,k)%Topo%ecnryz
                 END IF
                 !........................................
                 !Write(OutUnitBpv,*) "                  i =", i, " j = ",j, " k = ",k, " ec==1"
                 IF(i==ix0.OR.i==ix1) THEN
                   wbnfyz=wbnfyz+1
                   CALL WOutBFace(wbnfyz,i,j,k,"fyz", UpBoundsLayer(i,j,k)%BCell%bfyznr, is_eg, &
                           &   tropos_border,eg_listfyz) 
                 ELSE 
                   wbnfyz=wbnfyz+1
                   CALL WOutBFace(wbnfyz,i,j,k,"fyz", UpBoundsLayer(i,j,k)%BCell%bfyznr, is_eg, &
                           &   tropos_inside,eg_listfyz) 
                 END IF
               ELSE 
                 !Write(OutUnitBpv,*) "                  i =", i, " j = ",j, " k = ",k, " ec==0"
                 eg_listfyz(1)=UpBoundsLayer(i,j,k-1)%BCell%beynr
                 eg_listfyz(2)=UpBoundsLayer(i,j,k  )%BCell%beznr
                 eg_listfyz(3)=UpBoundsLayer(i,j,k  )%BCell%beynr
                 eg_listfyz(4)=UpBoundsLayer(i,j-1,k)%BCell%beznr
                 IF(i==ix0.OR.i==ix1) THEN
                   wbnfyz=wbnfyz+1
                   CALL WOutBFace(wbnfyz,i,j,k,"fyz", UpBoundsLayer(i,j,k)%BCell%bfyznr, 4, & 
                            &   tropos_border,eg_listfyz) 
                 ELSE 
                   wbnfyz=wbnfyz+1
                   CALL WOutBFace(wbnfyz,i,j,k,"fyz", UpBoundsLayer(i,j,k)%BCell%bfyznr, 4, & 
                            &   tropos_inside,eg_listfyz)
                 END IF
               END IF
            END IF
          ELSE
            !IF(k<=DimUpBoundsArray(i,j-1,2)) THEN
            IF(UpBoundsLayer(i,j,k)%BCell%bfyznr>0) THEN
              in_out=Vertices(i,j-1,k-1)%in_out+Vertices(i,j-1,k)%in_out+ &
                     Vertices(i,j,k)%in_out+Vertices(i,j,k-1)%in_out
              IF(in_out>=0) THEN
                 !Write(OutUnitBpv,*) "                  i =", i, " j = ",j, " k = ",k, "(in_out>=0)"
                 eg_listfyz(1)=UpBoundsLayer(i,j,k-1)%BCell%beynr
                 eg_listfyz(2)=UpBoundsLayer(i,j,k  )%BCell%beznr
                 eg_listfyz(3)=UpBoundsLayer(i,j,k  )%BCell%beynr
                 eg_listfyz(4)=UpBoundsLayer(i,j-1,k)%BCell%beznr
                 IF(i==ix0.OR.i==ix1) THEN
                   wbnfyz=wbnfyz+1
                   CALL WOutBFace(wbnfyz,i,j,k,"fyz", UpBoundsLayer(i,j,k)%BCell%bfyznr, 4, & 
                            &   tropos_border,eg_listfyz) 
                 ELSE 
                   wbnfyz=wbnfyz+1
                   CALL WOutBFace(wbnfyz,i,j,k,"fyz", UpBoundsLayer(i,j,k)%BCell%bfyznr, 4, & 
                            &   tropos_inside,eg_listfyz) 
                 END IF
              END IF
            END IF
          END IF
        END DO
      END DO
    END DO
    !CALL WriteBoundBorderWE_Faces
    CALL WriteBoundBorderWE_FacesYZ

    !-- Output--Cells------------------------------------------------------------
    CALL WOutNrBCell(sbcnr,sbcvcnr,sbchnr,sbcnr_bdr)
    DO i=ix0+1,ix1
      DO j=iy0+1,iy1
        DO k=iz0+1,iz1
        !DO k=DimUpBoundsArray(i,j,1)+1,DimUpBoundsArray(i,j,2)
          IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
             IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
                 & .OR. Cell(i,j,k)%Cell%in_out==8) THEN
                 IF (Cell(i,j,k)%Cell%vc==0) THEN
                   !Reihenfolge F1,F6,F4,F5,F3,F2; (F1,dann gegen Uhrzeiger,F2)
                   !Write(OutUnitBpv,*) "                    ",i,j,k,"-allok-vc=0"
                   ifn=1     !1 'F2'
                   IF(Topo(i,j,k)%Topo%fxynr>0) THEN
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fxynr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j,k)%BCell%bfxynr
                    ELSE
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                                           & Topo(i,j,k)%Topo%FCutNr
                   END IF
                   ifn=ifn+1  !2 'F6'
                   IF(Topo(i,j,k)%Topo%fyznr>0) THEN
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fyznr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j,k)%BCell%bfyznr 
                    ELSE
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                                      & Topo(i,j,k)%Topo%FCutNr
                   END IF 
                   ifn=ifn+1  !3 'F4'
                   IF(Topo(i,j,k)%Topo%fzxnr>0) THEN
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fzxnr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j,k)%BCell%bfzxnr
                    ELSE
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                                      & Topo(i,j,k)%Topo%FCutNr
                   END IF
                   ifn=ifn+1  !4  'F1'
                   IF(Topo(i,j,k-1)%Topo%fxynr>0) THEN
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k-1)%Topo%fxynr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j,k-1)%BCell%bfxynr
                    ELSE
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                                      & Topo(i,j,k)%Topo%FCutNr 
                   END IF
                   ifn=ifn+1  !5  'F3'
                   IF(Topo(i,j-1,k)%Topo%fzxnr>0) THEN
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j-1,k)%Topo%fzxnr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j-1,k)%BCell%bfzxnr
                    ELSE
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                                      & Topo(i,j,k)%Topo%FCutNr
                   END IF
                   ifn=ifn+1  !6  'F5'
                   IF(Topo(i-1,j,k)%Topo%fyznr>0) THEN
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i-1,j,k)%Topo%fyznr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i-1,j,k)%BCell%bfyznr
                    ELSE
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%FCutNr
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                                      & Topo(i,j,k)%Topo%FCutNr
                   END IF 
                   wbncvcnr=wbncvcnr+1
                   CALL WOutBCell(wbncvcnr,i,j,k,"cc",UpBoundsLayer(i,j,k)%BCell%ub_cnr, &
                          &   6, tropos_border, UpBoundsLayer(i,j,k)%BCell%cfn)
                 ELSE  !  IF(Cell(i,j,k)%Cell%vc>0) THEN
                   !Topo(i,j,k)%Topo%cfn(1)=Topo(i,j,k)%Topo%FCutNr 
                   UpBoundsLayer(i,j,k)%BCell%cfn(1)=Topo(i,j,k)%Topo%FCutNr
                   ifn=1
                   IF(Topo(i,j,k)%Topo%fxynr>0) THEN
                     ifn=ifn+1
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fxynr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j,k)%BCell%bfxynr
                   END IF
                   IF(Topo(i,j,k)%Topo%fyznr>0) THEN
                     ifn=ifn+1
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fyznr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j,k)%BCell%bfyznr 
                   END IF 
                   IF(Topo(i,j,k)%Topo%fzxnr>0) THEN
                     ifn=ifn+1
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k)%Topo%fzxnr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j,k)%BCell%bfzxnr
                   END IF
                   IF(Topo(i,j,k-1)%Topo%fxynr>0) THEN
                     ifn=ifn+1
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j,k-1)%Topo%fxynr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j,k-1)%BCell%bfxynr
                   END IF
                   IF(Topo(i,j-1,k)%Topo%fzxnr>0) THEN
                     ifn=ifn+1
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i,j-1,k)%Topo%fzxnr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                            & UpBoundsLayer(i,j-1,k)%BCell%bfzxnr
                   END IF
                   IF(Topo(i-1,j,k)%Topo%fyznr>0) THEN
                     ifn=ifn+1
                     !Topo(i,j,k)%Topo%cfn(ifn)=Topo(i-1,j,k)%Topo%fyznr 
                     UpBoundsLayer(i,j,k)%BCell%cfn(ifn)= &
                          & UpBoundsLayer(i-1,j,k)%BCell%bfyznr
                   END IF 
                   !Topo(i,j,k)%Topo%nfn=ifn
                   !Write(OutUnitBpv,*) "                       ",i,j,k," vc > 0 = ", &
                   !       &    Cell(i,j,k)%Cell%vc
                   UpBoundsLayer(i,j,k)%BCell%nfn=ifn
                   wbncvcnr=wbncvcnr+1
                   CALL WOutBCell(wbncvcnr,i,j,k,"cv",UpBoundsLayer(i,j,k)%BCell%ub_cnr, &
                          &  ifn, tropos_border,UpBoundsLayer(i,j,k)%BCell%cfn)
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
        !DO k=DimUpBoundsArray(i,j,1)+1,DimUpBoundsArray(i,j,2)
          !IF(k<=DimUpBoundsArray(i-1,j-1,2)) THEN
          IF(UpBoundsLayer(i,j,k)%BCell%ub_cnr>0) THEN 
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
                !Write(OutUnitBpv,*) "                       ",i,j,k," in_out >= 0  if"
                UpBoundsLayer(i,j,k)%BCell%cfn(1)=UpBoundsLayer(i,j,k-1)%BCell%bfxynr
                UpBoundsLayer(i,j,k)%BCell%cfn(2)=UpBoundsLayer(i,j,k  )%BCell%bfyznr
                UpBoundsLayer(i,j,k)%BCell%cfn(3)=UpBoundsLayer(i,j,k  )%BCell%bfzxnr
                UpBoundsLayer(i,j,k)%BCell%cfn(4)=UpBoundsLayer(i-1,j,k)%BCell%bfyznr 
                UpBoundsLayer(i,j,k)%BCell%cfn(5)=UpBoundsLayer(i,j-1,k)%BCell%bfzxnr
                UpBoundsLayer(i,j,k)%BCell%cfn(6)=UpBoundsLayer(i,j,k  )%BCell%bfxynr
                wbncnr=wbncnr+1
                CALL WOutBCell(wbncnr,i,j,k,"ch",UpBoundsLayer(i,j,k)%BCell%ub_cnr, 6, & 
                         &  tropos_inside,UpBoundsLayer(i,j,k)%BCell%cfn)
              ELSE
                !Write(OutUnitBpv,*) "                       ",i,j,k," in_out >= 0  else"
                !Topo(i,j,k)%Topo%ctp=-1  ! unterhalb Berg
              END IF
            END IF
          END IF 
        END DO
      END DO
    END DO
    CALL WriteBoundBorderWE_Cells
  END DO  !ib
  CALL EndeOutBound(FileName,nb)
  !------------------ Protokoll zur Ausgabe --------------------------------! 
  Write(OutUnitProt,*) leerzei3,"> Protokoll aus SUBROUTINE WriteUpperBoundsLayer()"
  Write(OutUnitProt,*) leerzei3,"  Schreiben in  Datei  *.pva.bound  *.bound"
  Write(OutUnitProt,*) leerzei3,"--------------------------------------------------"
  IF ((wbnex+wbney+wbnez+becnr) /= upbound_enr ) THEN
    Write(OutUnitProt,*) leerzei5,"analysierte Edges /= gelistete Edges", &
                &         upbound_enr,"/=",(wbnex+wbney+wbnez+becnr)
  ELSE
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Edges : analysierte == gelistete "
  END IF
  IF ((wbnfxy+wbnfyz+wbnfzx+wbnfc) /= upbound_fnr ) THEN
    Write(OutUnitProt,*) leerzei5,"analysierte Faces /= gelistete Faces", &
                &         upbound_fnr,"/=",(wbnfxy+wbnfyz+wbnfzx+wbnfc)
  ELSE
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Faces : analysierte == gelistete "
  END IF
  IF ((wbncnr+wbncvcnr) /= upbound_cnr ) THEN
    Write(OutUnitProt,*) leerzei5,"analysierte Cellen /= gelistete Cellen ", &
                &        upbound_cnr," /=",(wbncnr+wbncvcnr)
  ELSE
    Write(OutUnitProt,*) leerzei5,"i.o.  anz. Cellen: analysierte == gelistete "
  END IF
  Write(OutUnitProt,*) ""
  !------------------ Ende Protokoll zur Ausgabe --------------------------! 
!!    WRITE(OutUnitBpv,'(a8,i8)') cells,nr_cells
!!    DO ib=1,nb
!!      Write(*,*) "                                     Block : ",ib,"\/",nb
!!      CALL Set(Floor(ib))
!!      DO k=iz0+1,iz1
!!        DO j=iy0+1,iy1
!!          DO i=ix0+1,ix1
!!            !CALL WriteCellGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!!            IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!!              Write(OutUnitBpv,*) "...................................."
!!              Write(OutUnitBpv,*) "Cell(",i,",",j,",",k,")%Cell% ->"
!!              Write(OutUnitBpv,*) "  Face1%mp=",Cell(i,j,k)%Cell%Face1%mp
!!              Write(OutUnitBpv,*) "  Face2%mp=",Cell(i,j,k)%Cell%Face2%mp
!!              Write(OutUnitBpv,*) "  Face3%mp=",Cell(i,j,k)%Cell%Face3%mp
!!              Write(OutUnitBpv,*) "  Face4%mp=",Cell(i,j,k)%Cell%Face4%mp
!!              Write(OutUnitBpv,*) "  Face5%mp=",Cell(i,j,k)%Cell%Face5%mp
!!              Write(OutUnitBpv,*) "  Face6%mp=",Cell(i,j,k)%Cell%Face6%mp
!!              Write(OutUnitBpv,*) "  vc=",Cell(i,j,k)%Cell%vc
!!              WRITE(OutUnitBpv,*) (Cell(i,j,k)%Cell%VertCut(ic),ic=1,Cell(i,j,k)%Cell%vc)
!!              Write(OutUnitBpv,*) "...................................."
!!            END IF
!!          END DO
!!        END DO
!!      END DO
!!    END DO  !ib
  CALL WriteEndUpperBoundsLayerProt
  CALL CloseFileForBound 
END SUBROUTINE WriteUpperBoundsLayer


END MODULE OutUpperBoundsLayer_Mod
