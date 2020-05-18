MODULE InitTropoOroStruct_Mod
  
  USE F_Mod
  USE Floor_Mod
  USE IOControl_Mod
  USE Parametric_Mod
  USE OutputOutGmvG_Mod


  IMPLICIT NONE
CONTAINS

SUBROUTINE InitAllTopo
  INTEGER :: ib,i,j,k

  topo_cnr=0;topo_cvcnr=0
  topo_enr=0;topo_ecnr=0
  topo_fnr=0;topo_fcnr=0

  DO ib=1,nb
    CALL Set(Floor(ib))
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
            Topo(i,j,k)%Topo%ctp=-99
            Topo(i,j,k)%Topo%cnr=0
            Topo(i,j,k)%Topo%nfn=0
            !inside mountine
            Topo(i,j,k)%Topo%icnr=0
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
            Topo(i,j,k)%Topo%Vert=>Vertices(i,j,k)
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
            Topo(i,j,k)%Topo%ecnrxy=-99
            Topo(i,j,k)%Topo%ecnryz=-99
            Topo(i,j,k)%Topo%ecnrzx=-99
            Topo(i,j,k)%Topo%exnr=-99
            Topo(i,j,k)%Topo%eynr=-99
            Topo(i,j,k)%Topo%eznr=-99
            !inside mountine
            Topo(i,j,k)%Topo%iexnr=-99
            Topo(i,j,k)%Topo%ieynr=-99
            Topo(i,j,k)%Topo%ieznr=-99
            Topo(i,j,k)%Topo%iecnrxy=-99
            Topo(i,j,k)%Topo%iecnryz=-99
            Topo(i,j,k)%Topo%iecnrzx=-99
        END DO
      END DO
    END DO
    DO i=ix0,ix1
      DO j=iy0,iy1
        DO k=iz0,iz1
            Topo(i,j,k)%Topo%fxynr=-999
            Topo(i,j,k)%Topo%fzxnr=-999
            Topo(i,j,k)%Topo%fyznr=-999
            Topo(i,j,k)%Topo%FCutNr=-999
            Topo(i,j,k)%Topo%fceg_nr=0
            Topo(i,j,k)%Topo%fceg_nrcut(:)=0
            !inside mountine
            Topo(i,j,k)%Topo%ifxynr=-999
            Topo(i,j,k)%Topo%ifzxnr=-999
            Topo(i,j,k)%Topo%ifyznr=-999
        END DO
      END DO
    END DO
  END DO  ! ib
END SUBROUTINE InitAllTopo


SUBROUTINE InitUpBoundsLayer
  INTEGER :: ib,i,j,k

  !...................................................................
  upbound_ps=0       ! counter grid cut points of upper bounds
  upbound_p =0       ! counter grid points of upper bounds
                     ! akt.area upper bound with direction side east, -west 
  upbound_cnr=0      ! summe all area cell bounds
  upbound_cvc=0      ! summe all area cell-vc bounds
  upbound_enr=0      ! ""             edge
  upbound_ecnr=0     ! ""             ec
  upbound_fnr=0      ! ""             face
  upbound_fcnr=0     ! ""             fcnr
  gl_benr=0          ! counter, identical edge-number of borders-boundslayer 
                     ! in all direction w,e,s,n,b,t when border inside blocks

  ! allocation with 'w-e-n-s-border'
  !...................................................................
  DO ib=1,nb
    CALL Set(Floor(ib))
    DO i=ix0-1,ix1+1
      DO j=iy0-1,iy1+1
        DO k=iz0,iz1
            UpBoundsLayer(i,j,k)%BCell%i=-9
            UpBoundsLayer(i,j,k)%BCell%j=-9
            UpBoundsLayer(i,j,k)%BCell%k=-9
            UpBoundsLayer(i,j,k)%BCell%nrP=0
            UpBoundsLayer(i,j,k)%BCell%nrEcx=0
            UpBoundsLayer(i,j,k)%BCell%nrEcy=0
            UpBoundsLayer(i,j,k)%BCell%nrEcz=0
            !...................................
            UpBoundsLayer(i,j,k)%BCell%bexnr=-99
            UpBoundsLayer(i,j,k)%BCell%beynr=-99
            UpBoundsLayer(i,j,k)%BCell%beznr=-99
            UpBoundsLayer(i,j,k)%BCell%becnrxy=-99
            UpBoundsLayer(i,j,k)%BCell%becnryz=-99
            UpBoundsLayer(i,j,k)%BCell%becnrzx=-99
            !...................................
            UpBoundsLayer(i,j,k)%BCell%bfxynr=-999
            UpBoundsLayer(i,j,k)%BCell%bfzxnr=-999
            UpBoundsLayer(i,j,k)%BCell%bfyznr=-999
            UpBoundsLayer(i,j,k)%BCell%flisteg_nr(1:8)=0
            UpBoundsLayer(i,j,k)%BCell%feg_nr=0
            !...................................
            UpBoundsLayer(i,j,k)%BCell%FCutNr=-999
            UpBoundsLayer(i,j,k)%BCell%fceg_nr=0
            UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(:)=0
            !...................................
            UpBoundsLayer(i,j,k)%BCell%ub_cnr=0
            UpBoundsLayer(i,j,k)%BCell%ubc='n'
            UpBoundsLayer(i,j,k)%BCell%nfn=0
            UpBoundsLayer(i,j,k)%BCell%cfn(:)=0
        END DO
      END DO
    END DO
    DimUpBoundsArray(:,:,:)=0
    DimUpBoundsN(:,:,:)=0
    sbcnr_bdr(1:6)=0
    sbcvcnr=0
    sbchnr =0
    sbcnr  =0
    sbexnr=0
    sbeynr=0
    sbeznr=0
    sbecnr=0
    sbenr =0
    sbfcnr =0
    sbfxynr=0
    sbfzxnr=0
    sbfyznr=0
    sbfnr  =0
  END DO  ! ib
END SUBROUTINE InitUpBoundsLayer


SUBROUTINE ProveTropoCells(Flid)
  TYPE (Domain_T) :: Flid
  INTEGER :: i,j,k,in_out
  INTEGER :: cnr,cnr_cut
  INTEGER :: stbound_cnr,begtopo_cnr
  !......................................................
  !Tropoposphäre Cells
  !     -alle Cut-Cells: Troposphäre == UpperBoundsLayer
  !......................................................
  ! Für Analyze aktiv schalten, Ausgabe in *.prot-File
  !DO i=ix0+1,ix1
  !  DO j=iy0+1,iy1
  !    DO k=iz0+1,iz1
  !      !WRITE(OutUnitProt,*) "Cell(i,j,k)", i, j, k,"   WriteAllCellsAsciiGMV"
  !      !WRITE(*,*) "Cell(i,j,k)", i, j, k,"   WriteAllCellsAsciiGMV"
  !      !CALL WriteCellTopoAOut(Cell(i,j,k)%Cell,i,j,k)
  !      IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
  !       IF (i==2 .and. j==3 .and. k==1) THEN
  !        CALL WriteCellProt2(Cell(i,j,k)%Cell,i,j,k)
  !       END IF
  !      END IF
  !    END DO
  !  END DO
  !END DO
  !......................................................
  stbound_cnr=upbound_cnr
  begtopo_cnr=topo_cnr
  cnr_cut=begtopo_cnr
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      DO k=iz0+1,iz1
        IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
           IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
               .OR.Cell(i,j,k)%Cell%in_out==8) THEN
             !IF(Cell(i,j,k)%Cell%vc==0) THEN 
             !ELSE ! vc>0 
             !END IF   !beide eingeschossen
             Topo(i,j,k)%Topo%ctp=0 
             cnr_cut=cnr_cut+1
             Topo(i,j,k)%Topo%cnr=cnr_cut
             upbound_cnr=upbound_cnr+1
             UpBoundsLayer(i,j,k)%BCell%ub_cnr=upbound_cnr
           END IF
        !ELSE !!! .NOT. ASSOCIATED(Cell)
        !  Betreff Ausgabe-Reihenfolge 
        !          extra Schleife für Cellen oberhalb
        END IF
      END DO
    END DO
  END DO
  cnr=cnr_cut
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      DO k=iz0+1,iz1
        IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
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
            Topo(i,j,k)%Topo%ctp=1
            cnr=cnr+1
            Topo(i,j,k)%Topo%cnr=cnr
          ELSE
            Topo(i,j,k)%Topo%ctp=-1
          END IF
        END IF
      END DO
    END DO
  END DO
  !.......................................
  topo_cnr=cnr
  topo_cvcnr=topo_cvcnr+cnr_cut-begtopo_cnr
  Flid%stopo_cnr=topo_cnr-begtopo_cnr
  Flid%stopo_cv=cnr_cut-begtopo_cnr
  !.................................
  !bound only cuts zugewiesen
  Flid%sbcvcnr=upbound_cnr-stbound_cnr
  upbound_cvc=upbound_cvc+Flid%sbcvcnr
  !.................................
END SUBROUTINE ProveTropoCells

SUBROUTINE ProveTOroCells(Flid)
  TYPE (Domain_T) :: Flid

  INTEGER :: i,j,k,cnr
  !.............................................. 
  ! Orography Cellen (TOroCells)
  !   - Zählung Cellen fortlaufend der Liste Tropo
  !   - Counting Cells consecutively of the list Tropo
  !..............................................
  cnr=topo_cnr
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      DO k=iz0+1,iz1
         !in_out checken damit nicht durch troposphäre abfrage läuft
         IF (.NOT.ASSOCIATED(Cell(i,j,k)%Cell)) THEN
           IF (Topo(i,j,k)%Topo%ctp==-1) THEN
             cnr=cnr+1
             Topo(i,j,k)%Topo%icnr=cnr
           END IF
         ELSE   ! ASSOCIATED(..Celle)
           IF(Cell(i,j,k)%Cell%vc>0) THEN
             !Write(10,*) "Associated noch Auswertung Celle unterhalb Berg", i,j,k
             !Write(10,*) "                 vc= ",Cell(i,j,k)%Cell%vc, &
             !         &   "             in_out= ",Cell(i,j,k)%Cell%in_out
             !CALL WriteCellProt2(Cell(i,j,k)%Cell,i,j,k)
             cnr=cnr+1
             Topo(i,j,k)%Topo%icnr=cnr
           ELSE IF (Cell(i,j,k)%Cell%vc==0.AND.Cell(i,j,k)%Cell%in_out<=-4) THEN
             !Beachtung Cellen mit Edge-GrenzKante und Face_Grenze
             !Write(10,*) "Associated noch Auswertung Celle unterhalb Berg", i,j,k
             !Write(10,*) "                 vc= ",Cell(i,j,k)%Cell%vc, &
             !        &   "             in_out= ",Cell(i,j,k)%Cell%in_out
             cnr=cnr+1
             Topo(i,j,k)%Topo%icnr=cnr
             SELECT CASE (Cell(i,j,k)%Cell%in_out)
               CASE (-4)
                 Topo(i,j,k)%Topo%ctp=-4
               CASE (-6)
                 Topo(i,j,k)%Topo%ctp=-6
               CASE (-7)
                 Topo(i,j,k)%Topo%ctp=-7
               CASE DEFAULT
                 Topo(i,j,k)%Topo%ctp=-9
             END SELECT
           END IF
         END IF
      END DO
    END DO
  END DO
  para_grid_cnr=cnr
  para_oro_cnr=para_grid_cnr-topo_cnr
END SUBROUTINE ProveTOroCells


SUBROUTINE ProveDimUpBounds(Flid)
  TYPE (Domain_T) :: Flid

  INTEGER :: i,j,k,is,js,ks,ke,ie,je
  INTEGER :: ch_k1,ch_k2,ch_k_isj,ch_k_ijs,ch_k_isjs
  !.....................................................................
  ! DimUpBoundsArray(ix0:ix1,iy0:iy1,1:2)
  !    - beschreibt Index Z-plane der UpperBoundsLayer
  !    - 'iz' Grid-Vertex-Points und Schnitt-Vertex-Points suchen 
  !.........................................
  ! 'iz' Vertex-Points%in_out=0//1 der Bounds Layer suchen 
  ! an allen Knoten einer Schnitt-Celle angrenzende Cellen ohne Schnitt,
  ! als UpperBoundLayer-Celle kennzeichnen 
  !.....................................................................

    DimUpBoundsArray(:,:,1)=999999
    DimUpBoundsArray(:,:,2)=0
    DO k=iz0+1,iz1
      ke=MIN(k+1,iz1)
      DO i=ix0+1,ix1
        ie=MIN(i+1,ix1)
        DO j=iy0+1,iy1
          je=MIN(j+1,iy1)
          ch_k1=DimUpBoundsArray(i,j,1)
          ch_k2=DimUpBoundsArray(i,j,2)
          IF(Topo(i,j,k)%Topo%ctp==0) THEN
             UpBoundsLayer(i,j,k)%BCell%ubc='b'
             DimUpBoundsArray(i,j,1)=MIN(ch_k1,k-1)
             !"we","we,t"
             DO ks=k,ke         !k+1
              DO is=i-1,ie,2    !i+1,2
                IF(Topo(is,j,ks)%Topo%ctp==1) THEN
                 UpBoundsLayer(is,j,ks)%BCell%ubc='b'
                 ch_k_isj=DimUpBoundsArray(is,j,2)
                 DimUpBoundsArray(is,j,2)=MAX(ch_k_isj,ks)
                END IF
              END DO
             END DO
             !"sn", "sn,t"
             DO ks=k,ke         !k+1
              DO js=j-1,je,2    !j+1,2
                IF(Topo(i,js,ks)%Topo%ctp==1) THEN
                 UpBoundsLayer(i,js,ks)%BCell%ubc='b'
                 ch_k_ijs=DimUpBoundsArray(i,js,2)
                 DimUpBoundsArray(i,js,2)=MAX(ch_k_ijs,ks)
                END IF
              END DO
             END DO
             !"nw","ne","sw","se"
             !"nw,t","ne,t","sw,t","se,t"
             DO ks=k,ke         !k+1
              DO is=i-1,ie,2    !i+1,2
               DO js=j-1,je,2   !j+1,2
                IF(Topo(is,js,ks)%Topo%ctp==1) THEN
                 UpBoundsLayer(is,js,ks)%BCell%ubc='b'
                 ch_k_isjs=DimUpBoundsArray(is,js,2)
                 DimUpBoundsArray(is,js,2)=MAX(ch_k_isjs,ks)
                END IF
               END DO
              END DO
             END DO
             UpBoundsLayer(i,j,k+1)%BCell%ubc='b'
          ELSE IF(Topo(i,j,k)%Topo%ctp==1) THEN 
            IF(ch_k2==0) THEN
              UpBoundsLayer(i,j,k)%BCell%ubc='b'
              DimUpBoundsArray(i,j,2)=k
            END IF
          END IF
          !
        END DO
      END DO
    END DO
    ! ix0-belegen, Area ix0:ix1, erforderlich für all edges and faces der Celle-ix0+1
    i=ix0+1 
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        IF(UpBoundsLayer(i,j,k)%BCell%ubc=='b') THEN
            DimUpBoundsArray(ix0,j,:)=DimUpBoundsArray(i,j,:)
        END IF
      END DO
    END DO
    ! iy0-belegen, area iy0:iy1, erforderlich für edges and faces der celle-iy0+1
    j=iy0+1 
    DO k=iz0+1,iz1
      DO i=ix0+1,ix1
        IF(UpBoundsLayer(i,j,k)%BCell%ubc=='b') THEN
            DimUpBoundsArray(i,iy0,:)=DimUpBoundsArray(i,iy0+1,:)
        END IF
      END DO
    END DO
   ! ix1+1-belegen, Area ix0:ix1, erforderlich für all edges and faces der Celle-ix1+1
    i=ix1 
    DO k=iz0+1,iz1
      DO j=iy0+1,iy1
        IF(UpBoundsLayer(i,j,k)%BCell%ubc=='b') THEN
            DimUpBoundsArray(ix1+1,j,:)=DimUpBoundsArray(i,j,:)
        END IF
      END DO
    END DO
   ! iy1+1-belegen, area iy0:iy1, erforderlich für edges and faces der celle-iy1+1
    j=iy1 
    DO k=iz0+1,iz1
      DO i=ix0+1,ix1
        IF(UpBoundsLayer(i,j,k)%BCell%ubc=='b') THEN
            DimUpBoundsArray(i,iy1+1,:)=DimUpBoundsArray(i,iy0+1,:)
        END IF
      END DO
    END DO
  
    ! Point ix0,iy0 belegen
    DimUpBoundsArray(ix0,iy0,:)=DimUpBoundsArray(ix0,iy0+1,:)
  
    ! Point ix1+1,iy0
    DimUpBoundsArray(ix1+1,iy0,:)=DimUpBoundsArray(ix1,iy0,:)
  
    IF(DimUpBoundsArray(ix1+1,iy0,2)==0) THEN
       DimUpBoundsArray(ix1+1,iy0,2)=iz1  ! Schnitt bis maxZ-Gitter (Zylinder-Bsp.)
    END IF
  
    minZ_UBL=MINVAL(DimUpBoundsArray(ix0:ix1,iy0:iy1,1))
    maxZ_UBL=MAXVAL(DimUpBoundsArray(ix0:ix1,iy0:iy1,2))
    DO i=ix0,ix1
      DO j=iy0,iy1
        DimUpBoundsN(i,j,1:2)=DimUpBoundsArray(i,j,1:2)
      END DO
    END DO
  

END SUBROUTINE ProveDimUpBounds


SUBROUTINE ProveDimUpBounds_V1(Flid)
  TYPE (Domain_T) :: Flid

  INTEGER :: i,j,k,ix,iy,iz,ib
  INTEGER :: zminUpBL,zmaxUpBL
  INTEGER :: ctp_erf
  !........................................
  ! DimUpBoundsArray(ix0:ix1,iy0:iy1,1:2)
  !    - beschreibt Index Z-plane der UpperBoundsLayer
  !    - 'iz' Grid-Vertex-Points suchen 
  !.........................................
  ! 'iz' Grid-Vertex-Points%in_out=0//1 der Bounds Layer suchen 
  DimUpBoundsArray(:,:,:)=0
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      ctp_erf=0
      DO k=iz0+1,iz1
        IF(Topo(i,j,k)%Topo%ctp==0.AND.ctp_erf==0) THEN
          DimUpBoundsArray(i,j,1)=k-1
          ctp_erf=1
        ELSE IF (Topo(i,j,k)%Topo%ctp==1.AND.ctp_erf==1) THEN
          DimUpBoundsArray(i,j,2)=k
        Write(*,*) "DimUpBoundsArray\(",i,",",j,",2\)=",k
          EXIT
        END IF
      END DO
      IF(DimUpBoundsArray(i,j,2)==0) THEN
          DimUpBoundsArray(i,j,2)=iz1  ! Schnitt bis maxZ-Gitter (Zylinder-Bsp.)
      END IF
    END DO
  END DO
  ! ix0-belegen, Area ix0:ix1, erforderlich für all edges and faces der Celle-ix0+1 
  i=ix0+1
  DO j=iy0+1,iy1
    ctp_erf=0
    DO k=iz0+1,iz1
      IF(Topo(i,j,k)%Topo%ctp==0.AND.ctp_erf==0) THEN
        DimUpBoundsArray(ix0,j,1)=k-1
        ctp_erf=1
      ELSE IF (Topo(i,j,k)%Topo%ctp==1.AND.ctp_erf==1) THEN
        DimUpBoundsArray(ix0,j,2)=k
        EXIT
      END IF
    END DO
    IF(DimUpBoundsArray(ix0,j,2)==0) THEN
        DimUpBoundsArray(ix0,j,2)=iz1  ! Schnitt bis maxZ-Gitter (Zylinder-Bsp.)
    END IF
  END DO
  ! ix1+1-belegen, Area ix0:ix1, erforderlich für all edges and faces der Celle-ix0+1 
  i=ix1+1
  DO j=iy0+1,iy1
    ctp_erf=0
    DO k=iz0+1,iz1
      IF(Topo(ix1,j,k)%Topo%ctp==0.AND.ctp_erf==0) THEN
        DimUpBoundsArray(i,j,1)=k-1
        ctp_erf=1
      ELSE IF (Topo(ix1,j,k)%Topo%ctp==1.AND.ctp_erf==1) THEN
        DimUpBoundsArray(i,j,2)=k
        EXIT
      END IF
    END DO
    IF(DimUpBoundsArray(i,j,2)==0) THEN
        DimUpBoundsArray(i,j,2)=iz1  ! Schnitt bis maxZ-Gitter (Zylinder-Bsp.)
    END IF
  END DO
  
  ! iy0-belegen, area iy0:iy1, erforderlich für edges and faces der celle-iy0+1
  j=iy0+1
  DO i=ix0+1,ix1
    ctp_erf=0
    DO k=iz0+1,iz1
      IF(Topo(i,j,k)%Topo%ctp==0.AND.ctp_erf==0) THEN
        DimUpBoundsArray(i,iy0,1)=k-1
        ctp_erf=1
      ELSE IF (Topo(i,j,k)%Topo%ctp==1.AND.ctp_erf==1) THEN
        DimUpBoundsArray(i,iy0,2)=k
        Write(*,*) "DimUpBoundsArray\(",i,",",j-1,",2\)=",k
        EXIT
      END IF
    END DO
    IF(DimUpBoundsArray(i,iy0,2)==0) THEN
        DimUpBoundsArray(i,iy0,2)=iz1  ! Schnitt bis maxZ-Gitter (Zylinder-Bsp.)
    END IF
  END DO
  ! iy1+1-belegen, area iy0:iy1, erforderlich für edges and faces der celle-iy0+1
  j=iy1+1
  DO i=ix0+1,ix1
    ctp_erf=0
    DO k=iz0+1,iz1
      IF(Topo(i,iy1,k)%Topo%ctp==0.AND.ctp_erf==0) THEN
        DimUpBoundsArray(i,j,1)=k-1
        ctp_erf=1
      ELSE IF (Topo(i,iy1,k)%Topo%ctp==1.AND.ctp_erf==1) THEN
        DimUpBoundsArray(i,j,2)=k
        Write(*,*) "DimUpBoundsArray\(",i,",",j-1,",2\)=",k
        EXIT
      END IF
    END DO
    IF(DimUpBoundsArray(i,j,2)==0) THEN
        DimUpBoundsArray(i,j,2)=iz1  ! Schnitt bis maxZ-Gitter (Zylinder-Bsp.)
    END IF
  END DO
  !Point ix0,iy0
  i=ix0+1;j=iy0+1
  ctp_erf=0
  DO k=iz0+1,iz1
    IF(Topo(i,j,k)%Topo%ctp==0.AND.ctp_erf==0) THEN
       DimUpBoundsArray(ix0,iy0,1)=k-1
       ctp_erf=1
    ELSE IF (Topo(i,j,k)%Topo%ctp==1.AND.ctp_erf==1) THEN
       DimUpBoundsArray(ix0,iy0,2)=k
       EXIT
    END IF
  END DO
  IF(DimUpBoundsArray(ix0,iy0,2)==0) THEN
     DimUpBoundsArray(ix0,iy0,2)=iz1  ! Schnitt bis maxZ-Gitter (Zylinder-Bsp.)
  END IF
  !Point ix1+1,iy0
  i=ix1;j=iy0+1
  ctp_erf=0
  DO k=iz0+1,iz1
    IF(Topo(i,j,k)%Topo%ctp==0.AND.ctp_erf==0) THEN
       DimUpBoundsArray(ix1+1,iy0,1)=k-1
       ctp_erf=1
    ELSE IF (Topo(i,j,k)%Topo%ctp==1.AND.ctp_erf==1) THEN
       DimUpBoundsArray(ix1+1,iy0,2)=k
       EXIT
    END IF
  END DO
  IF(DimUpBoundsArray(ix1+1,iy0,2)==0) THEN
     DimUpBoundsArray(ix1+1,iy0,2)=iz1  ! Schnitt bis maxZ-Gitter (Zylinder-Bsp.)
  END IF
  minZ_UBL=MINVAL(DimUpBoundsArray(ix0:ix1,iy0:iy1,1))
  maxZ_UBL=MAXVAL(DimUpBoundsArray(ix0:ix1,iy0:iy1,2)) 
  DO i=ix0,ix1
    DO j=iy0,iy1
      DimUpBoundsN(i,j,1:2)=DimUpBoundsArray(i,j,1:2)
    END DO
  END DO

Write(*,*) (DimUpBoundsArray(i,0,2) ,i=ix0,ix1)
Write(*,*) (DimUpBoundsArray(i,0,1) ,i=ix0,ix1)
Write(*,*)
Write(*,*) (DimUpBoundsArray(i,0,2) ,i=ix0,ix1)
Write(*,*) (DimUpBoundsArray(i,0,1) ,i=ix0,ix1)

END SUBROUTINE ProveDimUpBounds_V1

SUBROUTINE ProveTropoEdges(Flid)
  TYPE (Domain_T) :: Flid

  INTEGER       :: i,j,k,li,iec
  INTEGER       :: ecnr,exnr,eynr,eznr,st_tenr
  INTEGER       :: bexnr,beynr,beznr,becnr,st_benr
  TYPE (Edge_T) :: Edg(1:4)
  INTEGER       :: Edg_bnr(1:4)
  LOGICAL       :: s_ec
  !.....................................................
  !Troposphäre Edges
  !--/EdgeX/EdgeY/EdgeZ/EdgeCut
  !  ReihenFolge: EdgeCut als lezten Block zwingend!
  !.....................................................
  exnr=topo_enr               ! 1.Edge_X
  st_tenr=topo_enr
  bexnr=upbound_enr
  st_benr=upbound_enr
  DO i=ix0+1,ix1
    DO j=iy0,iy1
      DO k=iz0,iz1
        IF (ASSOCIATED(Edges_X(i,j,k)%Edge)) THEN
           IF(MAX(Edges_X(i,j,k)%Edge%Vert1%in_out,Edges_X(i,j,k)%Edge%Vert2%in_out)>0 &
              .OR. &
              (Edges_X(i,j,k)%Edge%Vert1%in_out==Edges_X(i,j,k)%Edge%Vert2%in_out.AND. &
               Edges_X(i,j,k)%Edge%Vert1%in_out==0)) THEN
             exnr=exnr+1
             Topo(i,j,k)%Topo%exnr=exnr
             Edges_X(i,j,k)%Edge%eg_nr=exnr
             !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
             IF(UpBoundsLayer(i,j,k)%BCell%bexnr<=0)THEN
               bexnr=bexnr+1
               UpBoundsLayer(i,j,k)%BCell%bexnr=bexnr
             END IF
             !END IF
           END IF
        ELSE
          IF(MAX(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)>0 ) THEN
            exnr=exnr+1
            Topo(i,j,k)%Topo%exnr=exnr
            !IF(j/=iy0.and.k/=iz0) THEN
            !  CALL ProveBoundEdgX(i,j,k,bexnr)
            !  bexnr=sbexnr
            !END IF 
          ELSE IF (Vertices(i-1,j,k)%in_out==Vertices(i,j,k)%in_out .AND. &
                   Vertices(i,j,k)%in_out==0) THEN
            exnr=exnr+1
            Topo(i,j,k)%Topo%exnr=exnr
            !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
            IF(UpBoundsLayer(i,j,k)%BCell%bexnr<=0)THEN
                bexnr=bexnr+1
                UpBoundsLayer(i,j,k)%BCell%bexnr=bexnr
            END IF
            !END IF
          END IF 
        END IF
!              Write(*,*) "exnr= ",exnr, Topo(i,j,k)%Topo%exnr,"  "," i, j, k :", i,j,k 
      END DO
    END DO
  END DO
  eynr=exnr   ! 2.Edge_Y
  beynr=bexnr
  DO i=ix0,ix1
    DO j=iy0+1,iy1
      DO k=iz0,iz1
        IF (ASSOCIATED(Edges_Y(i,j,k)%Edge)) THEN
          IF(MAX(Edges_Y(i,j,k)%Edge%Vert1%in_out,Edges_Y(i,j,k)%Edge%Vert2%in_out)>0 &
             .OR. &
             (Edges_Y(i,j,k)%Edge%Vert1%in_out==Edges_Y(i,j,k)%Edge%Vert2%in_out.AND. &
              Edges_Y(i,j,k)%Edge%Vert1%in_out==0)) THEN
            eynr=eynr+1
            Topo(i,j,k)%Topo%eynr=eynr
            Edges_Y(i,j,k)%Edge%eg_nr=eynr
            !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
              !IF((Topo(i,j,k)%Topo%ctp==0).AND. &
              !    & UpBoundsLayer(i-1,j,k)%BCell%beynr==0) THEN
              !    ! Special celle%in_out==6 .and. vc==0 
              !    beynr=beynr+1
              !    UpBoundsLayer(i-1,j,k)%BCell%beynr=beynr=beynr+1
              !END IF
             IF(UpBoundsLayer(i,j,k)%BCell%beynr<=0) THEN
              beynr=beynr+1
              UpBoundsLayer(i,j,k)%BCell%beynr=beynr
             END IF
            !END IF
          END IF
        ELSE
          IF(MAX(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)>0) THEN 
             eynr=eynr+1
             Topo(i,j,k)%Topo%eynr=eynr
          ELSE IF (Vertices(i,j-1,k)%in_out==Vertices(i,j,k)%in_out .AND. &
                  Vertices(i,j,k)%in_out==0) THEN
             eynr=eynr+1
             Topo(i,j,k)%Topo%eynr=eynr
             !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
             IF(UpBoundsLayer(i,j,k)%BCell%beynr<=0) THEN
                beynr=beynr+1
                UpBoundsLayer(i,j,k)%BCell%beynr=beynr
             END IF
             !END IF
          END IF
        END IF
      END DO
    END DO
  END DO
  eznr=eynr   ! 3.Edge_Z
  beznr=beynr
  DO i=ix0,ix1
    DO j=iy0,iy1
      DO k=iz0+1,iz1
        IF (ASSOCIATED(Edges_Z(i,j,k)%Edge)) THEN
          IF(MAX(Edges_Z(i,j,k)%Edge%Vert1%in_out,Edges_Z(i,j,k)%Edge%Vert2%in_out)>0 &
             .OR. &
             (Edges_Z(i,j,k)%Edge%Vert1%in_out==Edges_Z(i,j,k)%Edge%Vert2%in_out.AND. &
              Edges_Z(i,j,k)%Edge%Vert1%in_out==0)) THEN
            eznr=eznr+1
            Topo(i,j,k)%Topo%eznr=eznr
            Edges_Z(i,j,k)%Edge%eg_nr=eznr
            !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
              !IF((Topo(i,j,k)%Topo%ctp==0) THEN
              !  IF( UpBoundsLayer(i-1,j-1,k)%BCell%beznr==0) THEN
              !   beznr=beznr+1
              !   UpBoundsLayer(i-1,j-1,k)%BCell%beznr=beznr
              !  END IF
              !  IF( UpBoundsLayer(i-1,j,k)%BCell%beznr==0) THEN
              !   beznr=beznr+1
              !   UpBoundsLayer(i-1,j,k)%BCell%beznr=beznr
              !  END IF
              !END IF    
             IF(UpBoundsLayer(i,j,k)%BCell%beznr<=0) THEN
              beznr=beznr+1
              UpBoundsLayer(i,j,k)%BCell%beznr=beznr
             END IF
            !END IF
          END IF
        ELSE
          IF(MAX(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)>0) THEN 
             eznr=eznr+1
             Topo(i,j,k)%Topo%eznr=eznr
          ELSE IF (Vertices(i,j,k-1)%in_out==Vertices(i,j,k)%in_out.AND. &
                   Vertices(i,j,k)%in_out==0) THEN
             eznr=eznr+1
             Topo(i,j,k)%Topo%eznr=eznr
             !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
             IF(UpBoundsLayer(i,j,k)%BCell%beznr<=0) THEN
                beznr=beznr+1
                UpBoundsLayer(i,j,k)%BCell%beznr=beznr
             END IF
             !END IF
          END IF
        END IF
      END DO
    END DO
  END DO
  ecnr=eznr   ! 4.EdgeCut
  becnr=beznr
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      DO k=iz0+1,iz1
        IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
          IF(Cell(i,j,k)%Cell%vc>0) THEN
            !Edge-Cut zuordnen!!!
            Topo(i,j,k)%Topo%fceg_nr=Cell(i,j,k)%Cell%vc
            UpBoundsLayer(i,j,k)%BCell%fceg_nr=Cell(i,j,k)%Cell%vc 
            !----------------------------------------------------------------------------
            !Analyze geschnittene Face_XY-/Face_YZ-/Face_ZX-; ecnr zuordnen                 
            !Liste Topo(i,j,k)%Topo%fceg_nrcut(1:Cell(i,j,k)%Cell%vc) zuordnen
            !----------------------------------------------------------------------------
            DO iec=1,Cell(i,j,k)%Cell%vc-1
              !ueber die VertCutListe die Vert-Nummern ermitteln
              !im Output
              !IF(Faces_XY(i,j,k-1)%Face%ec==1) THEN
              IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1))  .OR. &
                 (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
                 Edg(1)=Faces_XY(i,j,k-1)%Face%Edge1
                 Edg(2)=Faces_XY(i,j,k-1)%Face%Edge2
                 Edg(3)=Faces_XY(i,j,k-1)%Face%Edge3
                 Edg(4)=Faces_XY(i,j,k-1)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i,j-1,k-1)%BCell%bexnr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j  ,k-1)%BCell%beynr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j  ,k-1)%BCell%bexnr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i-1,j,k-1)%BCell%beynr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j,k-1)%Topo%ecnrxy         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE. 
                       !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                         ! UpBoundsLayer(i,j,k-1)%BCell%becnrxy      =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE. 
                       !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k-1)%BCell%becnrxy        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr==0) THEN
                      ecnr=ecnr+1
                      Faces_XY(i,j,k-1)%Face%egc_nr=ecnr
                      Topo(i,j,k-1)%Topo%ecnrxy=ecnr
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr>0) THEN
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k-1)%Face%egc_nr
                 END IF
                 !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                   IF(s_ec .and. UpBoundsLayer(i,j,k-1)%BCell%becnrxy<=0) THEN
                        becnr=becnr+1
                        UpBoundsLayer(i,j,k-1)%BCell%becnrxy=becnr
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                   ELSE IF(s_ec .and. UpBoundsLayer(i,j,k-1)%BCell%becnrxy>0) THEN
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                            & UpBoundsLayer(i,j,k-1)%BCell%becnrxy
                   END IF
                 !END IF
              END IF
              !IF(Faces_XY(i,j,k)%Face%ec==1) THEN
              IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
                 (Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
                 Edg(1)=Faces_XY(i,j,k)%Face%Edge1
                 Edg(2)=Faces_XY(i,j,k)%Face%Edge2
                 Edg(3)=Faces_XY(i,j,k)%Face%Edge3
                 Edg(4)=Faces_XY(i,j,k)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i,j-1,k)%BCell%bexnr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j  ,k)%BCell%beynr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j  ,k)%BCell%bexnr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i-1,j,k)%BCell%beynr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                         ! UpBoundsLayer(i,j,k)%BCell%becnrxy        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE. 
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                         ! UpBoundsLayer(i,j,k)%BCell%becnrxy        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr==0) THEN
                     ecnr=ecnr+1
                     Faces_XY(i,j,k)%Face%egc_nr=ecnr
                     Topo(i,j,k)%Topo%ecnrxy=ecnr
                     Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr>0) THEN
                     Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k)%Face%egc_nr
                 END IF
                 !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                   IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnrxy<=0) THEN
                        becnr=becnr+1
                        UpBoundsLayer(i,j,k)%BCell%becnrxy=becnr
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                   ELSE IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnrxy>0) THEN
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                            & UpBoundsLayer(i,j,k)%BCell%becnrxy
                   END IF
                 !END IF
              END IF
              !IF(Faces_YZ(i-1,j,k)%Face%ec==1) THEN
              IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
                 (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) )THEN
                  Edg(1)=Faces_YZ(i-1,j,k)%Face%Edge1
                  Edg(2)=Faces_YZ(i-1,j,k)%Face%Edge2
                  Edg(3)=Faces_YZ(i-1,j,k)%Face%Edge3
                  Edg(4)=Faces_YZ(i-1,j,k)%Face%Edge4
                  s_ec=.TRUE.
                  Edg_bnr(1)=UpBoundsLayer(i-1,j,k-1)%BCell%beynr !..Face%Edge1-benr
                  Edg_bnr(2)=UpBoundsLayer(i-1,j,k)%BCell%beznr !..Face%Edge2-benr
                  Edg_bnr(3)=UpBoundsLayer(i-1,j,k)%BCell%beynr !..Face%Edge3-benr
                  Edg_bnr(4)=UpBoundsLayer(i-1,j-1,k)%BCell%beznr !..Face%Edge4-benr
                  DO li=1,4
                    IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                        !Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
                        Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                        s_ec=.FALSE. 
                        !IF(k>=DimUpBoundsArray(i-1,j,1).AND.k<=DimUpBoundsArray(i-1,j,2)) THEN
                           !UpBoundsLayer(i-1,j,k)%BCell%becnryz      =Edg_bnr(li)
                           UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                        !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                        !Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
                        Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                        s_ec=.FALSE. 
                        !IF(k>=DimUpBoundsArray(i-1,j,1).AND.k<=DimUpBoundsArray(i-1,j,2)) THEN
                           !UpBoundsLayer(i-1,j,k)%BCell%becnryz      =Edg_bnr(li)
                           UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                        !END IF
                    END IF
                  END DO
                  IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr==0) THEN
                       ecnr=ecnr+1
                       Faces_YZ(i-1,j,k)%Face%egc_nr=ecnr
                       Topo(i-1,j,k)%Topo%ecnryz=ecnr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                  ELSE IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr>0) THEN
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i-1,j,k)%Face%egc_nr
                  END IF
                  !IF(k>=DimUpBoundsArray(i-1,j,1).AND.k<=DimUpBoundsArray(i-1,j,2)) THEN
                    IF(s_ec .and. UpBoundsLayer(i-1,j,k)%BCell%becnryz<=0) THEN
                         becnr=becnr+1
                         UpBoundsLayer(i-1,j,k)%BCell%becnryz=becnr 
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                    ELSE IF(s_ec .and. UpBoundsLayer(i-1,j,k)%BCell%becnryz>0) THEN
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                               & UpBoundsLayer(i-1,j,k)%BCell%becnryz
                    END IF
                  !END IF
              END IF
              !IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
              IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
                 (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
                 Edg(1)=Faces_YZ(i,j,k)%Face%Edge1
                 Edg(2)=Faces_YZ(i,j,k)%Face%Edge2
                 Edg(3)=Faces_YZ(i,j,k)%Face%Edge3
                 Edg(4)=Faces_YZ(i,j,k)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i,j,k-1)%BCell%beynr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j,k)%BCell%beznr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j,k)%BCell%beynr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i,j-1,k)%BCell%beznr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k)%BCell%becnryz        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                        !Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
                        Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                        s_ec=.FALSE. 
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k)%BCell%becnryz        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr==0) THEN
                      ecnr=ecnr+1
                      Faces_YZ(i,j,k)%Face%egc_nr=ecnr
                      Topo(i,j,k)%Topo%ecnryz=ecnr
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr>0) THEN
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i,j,k)%Face%egc_nr
                 END IF
                 !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                    IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnryz<=0) THEN
                         becnr=becnr+1
                         UpBoundsLayer(i,j,k)%BCell%becnryz=becnr 
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                    ELSE IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnryz>0) THEN
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                               & UpBoundsLayer(i,j,k)%BCell%becnryz
                    END IF
                 !END IF
              END IF
              !IF(Faces_ZX(i,j-1,k)%Face%ec==1) THEN
              IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
                 (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
                 Edg(1)=Faces_ZX(i,j-1,k)%Face%Edge1
                 Edg(2)=Faces_ZX(i,j-1,k)%Face%Edge2
                 Edg(3)=Faces_ZX(i,j-1,k)%Face%Edge3
                 Edg(4)=Faces_ZX(i,j-1,k)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i-1,j-1,k)%BCell%beznr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j-1,k)%BCell%bexnr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j-1,k)%BCell%beznr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i,j-1,k-1)%BCell%bexnr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j-1,1).AND.k<=DimUpBoundsArray(i,j-1,2)) THEN
                         !UpBoundsLayer(i,j-1,k)%BCell%becnrzx      =Edg_bnr(li)
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j-1,1).AND.k<=DimUpBoundsArray(i,j-1,2)) THEN
                         !UpBoundsLayer(i,j-1,k)%BCell%becnrzx      =Edg_bnr(li)
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    END IF
                 END DO
                 IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr==0) THEN
                      ecnr=ecnr+1
                      Faces_ZX(i,j-1,k)%Face%egc_nr=ecnr
                      Topo(i,j-1,k)%Topo%ecnrzx=ecnr
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr>0) THEN
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j-1,k)%Face%egc_nr
                 END IF
                 !IF(k>=DimUpBoundsArray(i,j-1,1).AND.k<=DimUpBoundsArray(i,j-1,2)) THEN
                   IF(s_ec .and. UpBoundsLayer(i,j-1,k)%BCell%becnrzx<=0) THEN
                        becnr=becnr+1
                        UpBoundsLayer(i,j-1,k)%BCell%becnrzx      =becnr 
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                   ELSE IF(s_ec .and. UpBoundsLayer(i,j-1,k)%BCell%becnrzx>0) THEN
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                             &  UpBoundsLayer(i,j-1,k)%BCell%becnrzx
                   END IF
                 !END IF
              END IF
              !IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
              IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec+1)) .OR. &
                 (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec+1)) ) THEN
                 Edg(1)=Faces_ZX(i,j,k)%Face%Edge1
                 Edg(2)=Faces_ZX(i,j,k)%Face%Edge2
                 Edg(3)=Faces_ZX(i,j,k)%Face%Edge3
                 Edg(4)=Faces_ZX(i,j,k)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i-1,j,k)%BCell%beznr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j,k)%BCell%bexnr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j,k)%BCell%beznr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i,j,k-1)%BCell%bexnr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                         !UpBoundsLayer(i,j,k)%BCell%becnrzx        =Edg_bnr(li)
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                         !UpBoundsLayer(i,j,k)%BCell%becnrzx        =Edg_bnr(li)
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr==0) THEN
                      ecnr=ecnr+1
                      Faces_ZX(i,j,k)%Face%egc_nr=ecnr
                      Topo(i,j,k)%Topo%ecnrzx=ecnr
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr>0) THEN
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j,k)%Face%egc_nr
                 END IF
                 !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                   IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnrzx<=0) THEN
                        becnr=becnr+1
                        UpBoundsLayer(i,j,k)%BCell%becnrzx        =becnr 
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                   ELSE IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnrzx>0) THEN
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                             &  UpBoundsLayer(i,j,k)%BCell%becnrzx
                   END IF
                 !END IF
              END IF
            END DO
              !IF(Faces_XY(i,j,k-1)%Face%ec==1) THEN
              IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
                 (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
                 Edg(1)=Faces_XY(i,j,k-1)%Face%Edge1
                 Edg(2)=Faces_XY(i,j,k-1)%Face%Edge2
                 Edg(3)=Faces_XY(i,j,k-1)%Face%Edge3
                 Edg(4)=Faces_XY(i,j,k-1)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i,j-1,k-1)%BCell%bexnr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j  ,k-1)%BCell%beynr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j  ,k-1)%BCell%bexnr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i-1,j,k-1)%BCell%beynr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE. 
                       !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k-1)%BCell%becnrxy      =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k-1)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_XY(i,j,k-1)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j,k-1)%Topo%ecnrxy       =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE. 
                       !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k-1)%BCell%becnrxy      =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr==0) THEN
                      ecnr=ecnr+1
                      Faces_XY(i,j,k-1)%Face%egc_nr=ecnr
                      Topo(i,j,k-1)%Topo%ecnrxy=ecnr
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_XY(i,j,k-1)%Face%egc_nr>0) THEN
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k-1)%Face%egc_nr
                 END IF
                 !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                   IF(s_ec .and. UpBoundsLayer(i,j,k-1)%BCell%becnrxy<=0) THEN
                        becnr=becnr+1
                        UpBoundsLayer(i,j,k-1)%BCell%becnrxy=becnr
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                   ELSE IF(s_ec .and. UpBoundsLayer(i,j,k-1)%BCell%becnrxy>0) THEN
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                            & UpBoundsLayer(i,j,k-1)%BCell%becnrxy
                   END IF
                 !END IF
              END IF
              !IF(Faces_XY(i,j,k)%Face%ec==1) THEN
              IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )    .OR. & 
                 (Faces_XY(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_XY(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
                 Edg(1)=Faces_XY(i,j,k)%Face%Edge1
                 Edg(2)=Faces_XY(i,j,k)%Face%Edge2
                 Edg(3)=Faces_XY(i,j,k)%Face%Edge3
                 Edg(4)=Faces_XY(i,j,k)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i,j-1,k)%BCell%bexnr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j  ,k)%BCell%beynr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j  ,k)%BCell%bexnr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i-1,j,k)%BCell%beynr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE. 
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k)%BCell%becnrxy        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_XY(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_XY(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j,k)%Topo%ecnrxy         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE. 
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k)%BCell%becnrxy        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr==0) THEN
                     ecnr=ecnr+1
                     Faces_XY(i,j,k)%Face%egc_nr=ecnr
                     Topo(i,j,k)%Topo%ecnrxy=ecnr
                     Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_XY(i,j,k)%Face%egc_nr>0) THEN
                     Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_XY(i,j,k)%Face%egc_nr
                 END IF
                 !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                   IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnrxy<=0) THEN
                        becnr=becnr+1
                        UpBoundsLayer(i,j,k)%BCell%becnrxy=becnr
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                   ELSE IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnrxy>0) THEN
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                            & UpBoundsLayer(i,j,k)%BCell%becnrxy
                   END IF
                 !END IF
              END IF
              !IF(Faces_YZ(i-1,j,k)%Face%ec==1) THEN
              IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
                 (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
                  Edg(1)=Faces_YZ(i-1,j,k)%Face%Edge1
                  Edg(2)=Faces_YZ(i-1,j,k)%Face%Edge2
                  Edg(3)=Faces_YZ(i-1,j,k)%Face%Edge3
                  Edg(4)=Faces_YZ(i-1,j,k)%Face%Edge4
                  s_ec=.TRUE.
                  Edg_bnr(1)=UpBoundsLayer(i-1,j,k-1)%BCell%beynr !..Face%Edge1-benr
                  Edg_bnr(2)=UpBoundsLayer(i-1,j,k)%BCell%beznr !..Face%Edge2-benr
                  Edg_bnr(3)=UpBoundsLayer(i-1,j,k)%BCell%beynr !..Face%Edge3-benr
                  Edg_bnr(4)=UpBoundsLayer(i-1,j-1,k)%BCell%beznr !..Face%Edge4-benr
                  DO li=1,4
                    IF((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                        !Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
                        Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                        s_ec=.FALSE. 
                        !IF(k>=DimUpBoundsArray(i-1,j,1).AND.k<=DimUpBoundsArray(i-1,j,2)) THEN
                           !UpBoundsLayer(i-1,j,k)%BCell%becnryz      =Edg_bnr(li)
                           UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                        !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i-1,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_YZ(i-1,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                        !Topo(i-1,j,k)%Topo%ecnryz       =Edg(li)%eg_nr
                        Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                        s_ec=.FALSE. 
                        !IF(k>=DimUpBoundsArray(i-1,j,1).AND.k<=DimUpBoundsArray(i-1,j,2)) THEN
                           !UpBoundsLayer(i-1,j,k)%BCell%becnryz      =Edg_bnr(li)
                           UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                        !END IF
                    END IF
                  END DO
                  IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr==0) THEN
                       ecnr=ecnr+1
                       Faces_YZ(i-1,j,k)%Face%egc_nr=ecnr
                       Topo(i-1,j,k)%Topo%ecnryz=ecnr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                  ELSE IF(s_ec .and. Faces_YZ(i-1,j,k)%Face%egc_nr>0) THEN
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i-1,j,k)%Face%egc_nr
                  END IF
                  !IF(k>=DimUpBoundsArray(i-1,j,1).AND.k<=DimUpBoundsArray(i-1,j,2)) THEN
                    IF(s_ec .and. UpBoundsLayer(i-1,j,k)%BCell%becnryz<=0) THEN
                         becnr=becnr+1
                         UpBoundsLayer(i-1,j,k)%BCell%becnryz=becnr 
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                    ELSE IF(s_ec .and. UpBoundsLayer(i-1,j,k)%BCell%becnryz>0) THEN
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                               & UpBoundsLayer(i-1,j,k)%BCell%becnryz
                    END IF
                  !END IF
              END IF
              !IF(Faces_YZ(i,j,k)%Face%ec==1) THEN
              IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
                 (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_YZ(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
                 Edg(1)=Faces_YZ(i,j,k)%Face%Edge1
                 Edg(2)=Faces_YZ(i,j,k)%Face%Edge2
                 Edg(3)=Faces_YZ(i,j,k)%Face%Edge3
                 Edg(4)=Faces_YZ(i,j,k)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i,j,k-1)%BCell%beynr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j,k)%BCell%beznr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j,k)%BCell%beynr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i,j-1,k)%BCell%beznr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k)%BCell%becnryz        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_YZ(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_YZ(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j,k)%Topo%ecnryz         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                          !UpBoundsLayer(i,j,k)%BCell%becnryz        =Edg_bnr(li)
                          UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr==0) THEN
                      ecnr=ecnr+1
                      Faces_YZ(i,j,k)%Face%egc_nr=ecnr
                      Topo(i,j,k)%Topo%ecnryz=ecnr
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_YZ(i,j,k)%Face%egc_nr>0) THEN
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_YZ(i,j,k)%Face%egc_nr
                 END IF
                 !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                    IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnryz<=0) THEN
                         becnr=becnr+1
                         UpBoundsLayer(i,j,k)%BCell%becnryz=becnr 
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                    ELSE IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnryz>0) THEN
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                               & UpBoundsLayer(i,j,k)%BCell%becnryz
                    END IF
                 !END IF
              END IF
              !IF(Faces_ZX(i,j-1,k)%Face%ec==1) THEN
              IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
                 (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) )THEN
                 Edg(1)=Faces_ZX(i,j-1,k)%Face%Edge1
                 Edg(2)=Faces_ZX(i,j-1,k)%Face%Edge2
                 Edg(3)=Faces_ZX(i,j-1,k)%Face%Edge3
                 Edg(4)=Faces_ZX(i,j-1,k)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i-1,j-1,k)%BCell%beznr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j-1,k)%BCell%bexnr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j-1,k)%BCell%beznr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i,j-1,k-1)%BCell%bexnr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j-1,1).AND.k<=DimUpBoundsArray(i,j-1,2)) THEN
                         !UpBoundsLayer(i,j-1,k)%BCell%becnrzx      =Edg_bnr(li)
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j-1,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_ZX(i,j-1,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j-1,k)%Topo%ecnrzx       =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j-1,1).AND.k<=DimUpBoundsArray(i,j-1,2)) THEN
                         !UpBoundsLayer(i,j-1,k)%BCell%becnrzx      =Edg_bnr(li)
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr==0) THEN
                      ecnr=ecnr+1
                      Faces_ZX(i,j-1,k)%Face%egc_nr=ecnr
                      Topo(i,j-1,k)%Topo%ecnrzx=ecnr
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_ZX(i,j-1,k)%Face%egc_nr>0) THEN
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j-1,k)%Face%egc_nr
                 END IF
                 !IF(k>=DimUpBoundsArray(i,j-1,1).AND.k<=DimUpBoundsArray(i,j-1,2)) THEN
                   IF(s_ec .and. UpBoundsLayer(i,j-1,k)%BCell%becnrzx<=0) THEN
                        becnr=becnr+1
                        UpBoundsLayer(i,j-1,k)%BCell%becnrzx      =becnr 
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                   ELSE IF(s_ec .and. UpBoundsLayer(i,j-1,k)%BCell%becnrzx>0) THEN
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                             &  UpBoundsLayer(i,j-1,k)%BCell%becnrzx
                   END IF
                 !END IF
              END IF
              !IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
              IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(1) )   .OR. & 
                 (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Cell(i,j,k)%Cell%VertCut(iec) .AND. &
                  Faces_ZX(i,j,k)%Face%EdgeCut(1)==Cell(i,j,k)%Cell%VertCut(1) ) ) THEN
                 Edg(1)=Faces_ZX(i,j,k)%Face%Edge1
                 Edg(2)=Faces_ZX(i,j,k)%Face%Edge2
                 Edg(3)=Faces_ZX(i,j,k)%Face%Edge3
                 Edg(4)=Faces_ZX(i,j,k)%Face%Edge4
                 s_ec=.TRUE.
                 Edg_bnr(1)=UpBoundsLayer(i-1,j,k)%BCell%beznr !..Face%Edge1-benr
                 Edg_bnr(2)=UpBoundsLayer(i,j,k)%BCell%bexnr !..Face%Edge2-benr
                 Edg_bnr(3)=UpBoundsLayer(i,j,k)%BCell%beznr !..Face%Edge3-benr
                 Edg_bnr(4)=UpBoundsLayer(i,j,k-1)%BCell%bexnr !..Face%Edge4-benr
                 DO li=1,4
                   IF((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                       Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP)  .OR. &
                      (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                       Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP) ) THEN
                       !Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                         !UpBoundsLayer(i,j,k)%BCell%becnrzx        =Edg_bnr(li)
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                    ELSE IF(Edg(li)%yes_sp==1 .AND. &
                      ((Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert1%nrP .AND. &
                        Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert1%nrP .AND. &
                        Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%Vert2%nrP .AND. &
                        Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%VertS%nrP)  .OR. &
                       (Faces_ZX(i,j,k)%Face%EdgeCut(2)==Edg(li)%Vert2%nrP .AND. &
                        Faces_ZX(i,j,k)%Face%EdgeCut(1)==Edg(li)%VertS%nrP)) ) THEN
                       !Topo(i,j,k)%Topo%ecnrzx         =Edg(li)%eg_nr
                       Topo(i,j,k)%Topo%fceg_nrcut(iec)=Edg(li)%eg_nr
                       s_ec=.FALSE.
                       !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                         !UpBoundsLayer(i,j,k)%BCell%becnrzx        =Edg_bnr(li)
                         UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=Edg_bnr(li)
                       !END IF
                   END IF
                 END DO
                 IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr==0) THEN
                      ecnr=ecnr+1
                      Faces_ZX(i,j,k)%Face%egc_nr=ecnr
                      Topo(i,j,k)%Topo%ecnrzx=ecnr
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=ecnr
                 ELSE IF(s_ec .and. Faces_ZX(i,j,k)%Face%egc_nr>0) THEN
                      Topo(i,j,k)%Topo%fceg_nrcut(iec)=Faces_ZX(i,j,k)%Face%egc_nr
                 END IF
                 !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                   IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnrzx<=0) THEN
                        becnr=becnr+1
                        UpBoundsLayer(i,j,k)%BCell%becnrzx        =becnr 
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)=becnr
                   ELSE IF(s_ec .and. UpBoundsLayer(i,j,k)%BCell%becnrzx>0) THEN
                        UpBoundsLayer(i,j,k)%BCell%fceg_nrcut(iec)= &
                             &  UpBoundsLayer(i,j,k)%BCell%becnrzx
                   END IF
                 !END IF
              END IF   ! IF(Faces_ZX(i,j,k)%Face%ec==1) THEN
          !ELSE IF (Cell(i,j,k)%Cell%vc==0) THEN

          END IF    ! IF(Cell(i,j,k)%Cell%vc>0 || vc==0 Special Faces-kompl.Grenze)
        END IF   ! IF (ASSOCIATED(Cell(i,j,k)%Cell))
      END DO    ! k
    END DO    ! j
  END DO    ! i
  !..............................
  topo_enr=ecnr
  topo_ecnr=topo_ecnr+ecnr-eznr
  Flid%topo_pos_ec1=eznr+1
  Flid%topo_pos_ecn=ecnr
  Flid%stopo_ecnr=ecnr-eznr
  Flid%stopo_enr=ecnr-st_tenr
  !..............................
  upbound_enr=becnr
  upbound_ecnr=upbound_ecnr+becnr-beznr
  Flid%upb_pos_ec1=beznr+1
  Flid%upb_pos_ecn=becnr
  Flid%sbexnr=bexnr-st_benr
  Flid%sbeynr=beynr-bexnr
  Flid%sbeznr=beznr-beynr
  Flid%sbecnr=becnr-beznr
  Flid%sbenr=becnr-st_benr

END SUBROUTINE ProveTropoEdges

SUBROUTINE ProveTOroEdges(Flid)
  TYPE (Domain_T) :: Flid

  INTEGER :: i,j,k
  INTEGER :: exnr,eynr,eznr
  !............................................................
  ! Orography /EdgeX/EdgeY/EdgeZ  (TOroEdges)
  !   - Zählung Kanten-Celle fortlaufend der Liste Tropo
  !   - Counting Edges of Cells consecutively of the list Tropo
  !............................................................
  !            1.Edge_X
  exnr=topo_enr
  DO i=ix0+1,ix1
    DO j=iy0,iy1
      DO k=iz0,iz1
        IF(MIN(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
            exnr=exnr+1
            Topo(i,j,k)%Topo%iexnr=exnr
!W                !Write(*,*) "Orography: exnr= ",exnr, Topo(i,j,k)%Topo%exnr,"  "," i, j, k :", i,j,k 
        ELSE IF(MIN(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==0 .AND. &
                (Vertices(i-1,j,k)%in_out==0.and.Vertices(i,j,k)%in_out==0) ) THEN
               !Special Grenze Edge komplett)
            Topo(i,j,k)%Topo%iexnr=Topo(i,j,k)%Topo%exnr
        END IF
      END DO
    END DO
  END DO
  eynr=exnr   ! 2.Edge_Y
  DO i=ix0,ix1
    DO j=iy0+1,iy1
      DO k=iz0,iz1
        IF(MIN(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
            eynr=eynr+1
            Topo(i,j,k)%Topo%ieynr=eynr
!W                !Write(*,*) "Orography: eynr= ",eynr, Topo(i,j,k)%Topo%eynr,"  "," i, j, k :", i,j,k 
        ELSE IF(MIN(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)==0 .AND. &
                (Vertices(i,j-1,k)%in_out==0.and.Vertices(i,j,k)%in_out==0) ) THEN
               !Special Grenze Edge komplett)
            Topo(i,j,k)%Topo%ieynr=Topo(i,j,k)%Topo%eynr
        END IF
      END DO
    END DO
  END DO
  eznr=eynr   ! 3.Edge_Z
  DO i=ix0,ix1
    DO j=iy0,iy1
      DO k=iz0+1,iz1
        IF(MIN(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
            eznr=eznr+1
            Topo(i,j,k)%Topo%ieznr=eznr 
!W                !Write(*,*) "Orography: eznr= ",eznr, Topo(i,j,k)%Topo%eznr,"  "," i, j, k :", i,j,k 
        ELSE IF(MIN(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)==0 .AND. &
                (Vertices(i,j,k-1)%in_out==0.and.Vertices(i,j,k)%in_out==0)) THEN
               !Special Grenze Edge komplett)
            Topo(i,j,k)%Topo%ieznr=Topo(i,j,k)%Topo%eznr 
        END IF
      END DO
    END DO
  END DO
  para_grid_enr=eznr
  para_oro_enr=para_grid_enr-topo_enr
!        Write(*,*) "Orography_enr :",para_oro_enr
!        Write(*,*) "exnr =",topo_enr-exnr,      "  ",exnr
!        Write(*,*) "eynr =",eynr-exnr,          "  ",eynr
!        Write(*,*) "eznr =",eznr-eynr,          "  ",eznr
!        Write(*,*) "para_grid_enr =",para_grid_enr,"  "," grid_all_enr" 
!        Ende Orography-Edge
        !------------------------------------------------------------------------
END SUBROUTINE ProveTOroEdges

SUBROUTINE ProveTropoFaces(Flid)
  TYPE (Domain_T) :: Flid

  INTEGER :: i,j,k,in_out
  INTEGER :: fcnr,fxynr,fzxnr,fyznr
  INTEGER :: bfxynr,bfzxnr,bfyznr
  INTEGER :: st_topofnr,st_bfnr
  !......................................
  ! Troposphäre Faces
  ! FaceCut,Faces_XY,Faces_XZ,Faces_YZ
  ! Bound-Face:
  !   FaceCut,Faces_XY,Faces_XZ,Faces_YZ
  !......................................
  fcnr=topo_fnr
  st_topofnr=topo_fnr 
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      DO k=iz0+1,iz1
        IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
          IF(Cell(i,j,k)%Cell%vc>0) THEN
              fcnr=fcnr+1
              Topo(i,j,k)%Topo%FCutNr=fcnr
!W                  !Write(10,*) "Topo(i,j,k)%Topo%FCutNr =", fcnr, i,j,k
          !ELSE IF(Cell(i,j,k)%Cell%vc==0.AND.Cell(i,j,k)%Cell%in_out==4) THEN
          !    !Face1 Grenzfläche
          !    fcnr=fcnr+1
          !    Topo(i,j,k)%Topo%FCutNr=fcnr
          END IF
        END IF
      END DO
    END DO
  END DO
!W      ! Write(*,*) "fcnr=", fcnr, " --> anz. fcuts"
  fxynr=fcnr
  bfxynr=fcnr
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      DO k=iz0,iz1
        IF (ASSOCIATED(Faces_XY(i,j,k)%Face)) THEN
          IF(Vertices(i-1,j-1,k)%in_out==0 .and. &
             Vertices(i  ,j-1,k)%in_out==0 .and. &
             Vertices(i-1,j  ,k)%in_out==0 .and. &
             Vertices(i  ,j  ,k)%in_out==0) THEN
             !Special-Grenzfläche als 'fc' ausgewertet,
             !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
             !Write(10,*) "Faces_XY(i,j,k)%Face Special Face Grenze als fc",  i,j,k
             IF (k==iz0) THEN
                !Face1-Grenzfläche ( Hill.grid, Cellen 23/1/1-79/1/1 bei iz0 Grenze)
                !(Cell(i,j,k)%Cell%vc==0.and.Cell(i,j,k)%Cell%in_out==4)
                ! kein F-Cut!
                fxynr=fxynr+1
                Topo(i,j,k)%Topo%fxynr=fxynr
                bfxynr=bfxynr+1
                UpBoundsLayer(i,j,k)%BCell%bfxynr=bfxynr
             ELSE IF (Cell(i,j,k)%Cell%vc==0.and.Cell(i,j,k+1)%Cell%vc>0) THEN
                Topo(i,j,k)%Topo%fxynr=Topo(i,j,k+1)%Topo%FCutNr
                Topo(i,j,k)%Topo%FCutNr=Topo(i,j,k+1)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%bfxynr=Topo(i,j,k+1)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%FCutNr=Topo(i,j,k+1)%Topo%FCutNr
             ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
                Topo(i,j,k)%Topo%fxynr=Topo(i,j,k)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%bfxynr=Topo(i,j,k)%Topo%FCutNr
             ELSE
                Write(OutUnitProt,*)"Warning aus ProveTropoFaces(): "
                Write(OutUnitProt,*)"  Face_XY(i,j,k)-Special-Neu gefunden!"  ,i,j,k
             !   !CALL WriteFaceProt2(Faces_XY(i,j,k)%Face)
             END IF
          ELSE IF(Faces_XY(i,j,k)%Face%in_out>-4.and.Faces_XY(i,j,k)%Face%numberVert>2) THEN
             fxynr=fxynr+1 
             Topo(i,j,k)%Topo%fxynr=fxynr
             !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
                bfxynr=bfxynr+1 
                UpBoundsLayer(i,j,k)%BCell%bfxynr=bfxynr
             !END IF
          END IF
        ELSE
          in_out=Vertices(i-1,j-1,k)%in_out+Vertices(i,j-1,k)%in_out+ &
                 Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
          IF(in_out>=0) THEN
              fxynr=fxynr+1 
              Topo(i,j,k)%Topo%fxynr=fxynr
              !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
              !  IF(k<=DimUpBoundsArray(i-1,j,2)) THEN
              !    bfxynr=bfxynr+1 
              !    UpBoundsLayer(i,j,k)%BCell%bfxynr=bfxynr
              !  END IF
              !END IF
          END IF
        END IF
      END DO
    END DO
  END DO
  fzxnr=fxynr
  bfzxnr=bfxynr
  DO i=ix0+1,ix1
    DO j=iy0,iy1
      DO k=iz0+1,iz1
        IF (ASSOCIATED(Faces_ZX(i,j,k)%Face)) THEN
          IF(Vertices(i-1,j,k-1)%in_out==0 .and. &
             Vertices(i  ,j,k-1)%in_out==0 .and. & 
             Vertices(i-1,j,k  )%in_out==0 .and. &
             Vertices(i  ,j,k  )%in_out==0 ) THEN
             !Special Grenzfläche als 'fc' ausgewertet,
             !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
             IF(Cell(i,j,k)%Cell%vc==0.and.Cell(i,j+1,k)%Cell%vc>0) THEN
                Topo(i,j,k)%Topo%fzxnr=Topo(i,j+1,k)%Topo%FCutNr
                Topo(i,j,k)%Topo%FCutNr=Topo(i,j+1,k)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%bfzxnr=Topo(i,j+1,k)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%FCutNr=Topo(i,j+1,k)%Topo%FCutNr
             ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
                Topo(i,j,k)%Topo%fzxnr=Topo(i,j,k)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%bfzxnr=Topo(i,j,k)%Topo%FCutNr
             ELSE 
                Write(OutUnitProt,*)"Warning aus SetTropoOroParaView(): "
                Write(OutUnitProt,*)"  Face_ZX(i,j,k)-Special-Neu gefunden! "  ,i,j,k
                !CALL WriteFaceProt2(Faces_ZX(i,j,k)%Face)
             END IF
          ELSE IF(Faces_ZX(i,j,k)%Face%in_out>-4.and.Faces_ZX(i,j,k)%Face%numberVert>2) THEN
             fzxnr=fzxnr+1 
             Topo(i,j,k)%Topo%fzxnr=fzxnr
             !IF(k>=DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
               bfzxnr=bfzxnr+1 
               UpBoundsLayer(i,j,k)%BCell%bfzxnr=bfzxnr
             !END IF
          END IF
        ELSE
          in_out=Vertices(i-1,j,k-1)%in_out+Vertices(i,j,k-1)%in_out+ &
                 Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
          IF(in_out>=0) THEN
            fzxnr=fzxnr+1
            Topo(i,j,k)%Topo%fzxnr=fzxnr
            !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
            !  IF(k<=DimUpBoundsArray(i-1,j,2)) THEN
            !    bfzxnr=bfzxnr+1 
            !    UpBoundsLayer(i,j,k)%BCell%bfzxnr=bfzxnr
            !  END IF
            !END IF
          END IF
        END IF
      END DO
    END DO
  END DO
  fyznr=fzxnr
  bfyznr=bfzxnr
  DO i=ix0,ix1
    DO j=iy0+1,iy1
      DO k=iz0+1,iz1
        IF (ASSOCIATED(Faces_YZ(i,j,k)%Face)) THEN
          IF(Vertices(i,j-1,k-1)%in_out==0 .and. &
             Vertices(i,j-1,k  )%in_out==0 .and. &
             Vertices(i,j  ,k  )%in_out==0 .and. &
             Vertices(i,j  ,k-1)%in_out==0 ) THEN
             !Special Grenzfläche als 'fc' ausgewertet,
             !Copieren fc-nr der entsprechenden Topo%FCutNr der Celle oberhalb Berg
             IF(Cell(i,j,k)%Cell%vc==0.and.Cell(i+1,j,k)%Cell%vc>0) THEN
                Topo(i,j,k)%Topo%fyznr=Topo(i+1,j,k)%Topo%FCutNr
                Topo(i,j,k)%Topo%FCutNr=Topo(i+1,j,k)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%bfyznr=Topo(i+1,j,k)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%FCutNr=Topo(i+1,j,k)%Topo%FCutNr
             ELSE IF (Cell(i,j,k)%Cell%vc>0) THEN
                Topo(i,j,k)%Topo%fyznr=Topo(i,j,k)%Topo%FCutNr
                UpBoundsLayer(i,j,k)%BCell%bfyznr=Topo(i,j,k)%Topo%FCutNr
             ELSE 
                Write(OutUnitProt,*)"Warning aus SetTropoOroParaView(): "
                Write(OutUnitProt,*)"  Face_YZ(i,j,k)-Special-Neu gefunden! "  ,i,j,k
                !CALL WriteFaceProt2(Faces_YZ(i,j,k)%Face)
             END IF
          ELSE IF(Faces_YZ(i,j,k)%Face%in_out>-4.and.Faces_YZ(i,j,k)%Face%numberVert>2) THEN
            fyznr=fyznr+1 
            Topo(i,j,k)%Topo%fyznr=fyznr
            !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
              bfyznr=bfyznr+1 
              UpBoundsLayer(i,j,k)%BCell%bfyznr=bfyznr
            !END IF
          END IF
        ELSE
          in_out=Vertices(i,j-1,k-1)%in_out+Vertices(i,j-1,k)%in_out+ &
                 Vertices(i,j,k)%in_out+Vertices(i,j,k-1)%in_out
          IF(in_out>=0) THEN
            fyznr=fyznr+1 
            Topo(i,j,k)%Topo%fyznr=fyznr
            !IF(k>DimUpBoundsArray(i,j,1).AND.k<=DimUpBoundsArray(i,j,2)) THEN
            !  IF(k<=DimUpBoundsArray(i,j-1,2)) THEN
            !    bfyznr=bfyznr+1 
            !    UpBoundsLayer(i,j,k)%BCell%bfyznr=bfyznr
            !  END IF
            !END IF
          END IF
        END IF
      END DO
    END DO
  END DO
  !.........................................
  topo_fnr=fyznr
  topo_fcnr=topo_fcnr+fcnr-st_topofnr
  Flid%stopo_fnr=fyznr-st_topofnr
  Flid%stopo_fcnr=fcnr-st_topofnr
  !.........................................
  upbound_fnr=bfyznr
  upbound_fcnr=upbound_fcnr+fcnr-st_topofnr
  Flid%sbfcnr=fcnr-st_topofnr
  Flid%sbfxynr=bfxynr-fcnr
  Flid%sbfzxnr=bfzxnr-bfxynr
  Flid%sbfyznr=bfyznr-bfzxnr
  Flid%sbfnr=bfyznr-st_topofnr
  !.........................................
!!        Write(*,*) "topo_fnr:        ", "fnr-summe"
!!        Write(*,*) "fcnr  =",fcnr,        "  ",fcnr
!!        Write(*,*) "fxynr =",fxynr-fcnr,  "  ",fxynr
!!        Write(*,*) "fzxnr =",fzxnr-fxynr, "  ",fzxnr
!!        Write(*,*) "fyznr =",fyznr-fzxnr, "  ",fyznr
!!        Write(*,*) "upbound_fnr:     ", "fnr-summe bound (allok.+fcut)
!!        Write(*,*) "fcnr    =",fcnr,   "   ",fcnr, " identisch Tropo  
!!        Write(*,*) "sbfxynr =",sbfxynr,"   ",bfxynr
!!        Write(*,*) "sbfzxnr =",sbfzxnr,"   ",bfzxnr
!!        Write(*,*) "sbfyznr =",sbfyznr,"   ",bfyznr
END SUBROUTINE ProveTropoFaces

SUBROUTINE ProveTOroFaces(Flid)
  TYPE (Domain_T) :: Flid

  INTEGER :: i,j,k
  INTEGER :: fxynr,fzxnr,fyznr
  !.............................................................
  ! Orography: Faces (TOroFaces)
  !   - Zählung Fläche-Celle fortlaufend der Liste Tropo
  !   - Counting Faces of Cells consecutively of the list Tropo
  !.............................................................
  fxynr=topo_fnr
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      DO k=iz0,iz1
        IF( MIN(Vertices(i-1,j-1,k)%in_out,Vertices(i,j-1,k)%in_out, &
          &  Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
          fxynr=fxynr+1
          Topo(i,j,k)%Topo%ifxynr=fxynr
        END IF
      END DO
    END DO
  END DO
  fzxnr=fxynr
  DO i=ix0+1,ix1
    DO j=iy0,iy1
      DO k=iz0+1,iz1
        IF( MIN(Vertices(i-1,j,k-1)%in_out,Vertices(i,j,k-1)%in_out, &
          &  Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
           fzxnr=fzxnr+1
           Topo(i,j,k)%Topo%ifzxnr=fzxnr
         END IF
      END DO
    END DO
  END DO
  fyznr=fzxnr
  DO i=ix0,ix1
    DO j=iy0+1,iy1
      DO k=iz0+1,iz1
        IF( MIN(Vertices(i,j-1,k-1)%in_out,Vertices(i,j,k-1)%in_out, &
          &  Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)==-1) THEN
            fyznr=fyznr+1
            Topo(i,j,k)%Topo%ifyznr=fyznr
        END IF
      END DO
    END DO
  END DO
  para_grid_fnr=fyznr
  para_oro_fnr=para_grid_fnr-topo_fnr
!!        Write(*,*) "Oro_fnr:        ", "fnr-summe"
!!        Write(*,*) "fxynr =",fxynr-topo_fnr,  "  ",fxynr
!!        Write(*,*) "fzxnr =",fzxnr-fxynr, "  ",fzxnr
!!        Write(*,*) "fyznr =",fyznr-fzxnr, "  ",fyznr
!!        Write(*,*) "foronr=",para_oro_fnr, "   ",para_grid_fnr
END SUBROUTINE ProveTOroFaces


! Begin Bounds-Layer
SUBROUTINE ProveBoundsEdgeX(is,js,ks,Flid)
 INTEGER is,js,ks
 TYPE (Domain_T) :: Flid
 !local
 !Point wenn Celle%vc
 !EdgesX: (k-1,k,k+1)  | i,j+1  |  i+1,j+1 |
 !                     | i,j    |  i+1,j   |
 !                     | i,j-1  |  i+1,j-1 |
 !     ------    --------
 !    /7/8/9/   /16/17/18/     Reihenfolge jeweils 'i,j,k'
 !    ------ * ---------       Point '*'
 !   /4/5/6/  /13/14/15/
 !   ------   ---------
 !  /1/2/3/  /10/11/12/ 
 ! 
 INTEGER i,j,k
 DO i=is,is+1
   DO j=js-1,js+1
     DO k=ks-1,ks+1
       IF(k/=iz0-1.and.k/=iz1+1.and. &
          j/=iy0-1.and.j/=iy1+1.and. &
          i/=ix0  .and.i/=ix1+1) THEN
         IF(UpBoundsLayer(i,j,k)%BCell%bexnr<=0.AND. &
            (MAX(Vertices(i-1,j,k)%in_out,Vertices(i,j,k)%in_out)==1) ) THEN
           upbound_enr=upbound_enr+1
           Flid%sbexnr=Flid%sbexnr+1
           UpBoundsLayer(i,j,k)%BCell%bexnr=upbound_enr
         END IF
       END IF
     END DO
   END DO
 END DO
END SUBROUTINE ProveBoundsEdgeX

SUBROUTINE ProveBoundsEdgeY(is,js,ks,Flid)
 INTEGER is,js,ks
 TYPE (Domain_T) :: Flid
 !local
 !Point wenn Celle%vc
 !EdgesY: (k-1,k,k+1)   | i-1,j+1 | i,j+1  | i+1,j+1  |
 !                      | i-1,j   | i,j    | i+1,j    |
 !      ----  ----
 !    2/   4/   6/    Reihenfolge jeweils 'i,j,k'
 !    ---- *----      Point '*'
 !  1/   3/   5/
 !  ----  ----
 INTEGER i,j,k
 DO i=is-1,is+1
   DO j=js,js+1
     DO k=ks-1,ks+1
       IF(k/=iz0-1.and.k/=iz1+1.and. &
          j/=iy0  .and.j/=iy1+1.and. &
          i/=ix0-1.and.i/=ix1+1) THEN
          IF(UpBoundsLayer(i,j,k)%BCell%beynr<=0 .and. &
            (MAX(Vertices(i,j-1,k)%in_out,Vertices(i,j,k)%in_out)==1) ) THEN
            upbound_enr=upbound_enr+1
            Flid%sbeynr=Flid%sbeynr+1
            UpBoundsLayer(i,j,k)%BCell%beynr=upbound_enr
          END IF
       END IF
     END DO
   END DO
 END DO
END SUBROUTINE ProveBoundsEdgeY

SUBROUTINE ProveBoundsEdgeZ(is,js,ks,Flid)
 INTEGER is,js,ks
 TYPE (Domain_T) :: Flid
 !local
 !Point wenn Celle%vc
 !EdgesZ: (k,k+1)   | i-1,j+1  | i,j+1  | i+1,j+1  |
 !                  | i-1,j    | i,j    | i+1,j    |
 !                  | i-1,j-1  | i,j-1  | i+1,j-1  |
 !           -------- --------
 !      5/6 |  11/12 |  17/18 |
 !         -------- *--------
 !    3/4 |   9/10 |  15/16  |    Reihenfolge jeweils 'i,j,k'
 !       -------- ----------      Point '*'
 !  1/2 |    7/8 |   13/14 |
 INTEGER i,j,k
 DO i=is-1,is+1
   DO j=js-1,js+1
     DO k=ks,ks+1
       IF(k/=iz0  .and.k/=iz1+1.and. &
          j/=iy0-1.and.j/=iy1+1.and. &
          i/=ix0-1.and.i/=ix1+1) THEN
          IF(UpBoundsLayer(i,j,k)%BCell%beznr<=0 .and. &
            (MAX(Vertices(i,j,k-1)%in_out,Vertices(i,j,k)%in_out)==1) ) THEN
            upbound_enr=upbound_enr+1
            Flid%sbeznr=Flid%sbeznr+1
            UpBoundsLayer(i,j,k)%BCell%beznr=upbound_enr
          END IF
       END IF
     END DO
   END DO
 END DO
END SUBROUTINE ProveBoundsEdgeZ

SUBROUTINE ProveBoundsFaceXY(is,js,ks,Flid)
 INTEGER is,js,ks
 TYPE (Domain_T) :: Flid
 !local
 !Point wenn Celle%vc
 !FaceXY: ( k-1,k,k+1) | i,j+1 | i+1,j+1 | 
 !                     | i,j   | i+1,j   |
 !  ---------------------
 !  | 4/5/6  | 10/11/12 |  Reihenfolge jeweils 'i,j,k'
 !  ---------*-----------  Point '*'
 !  | 1/2/3  | 7/8/9    |
 !  ---------------------
 INTEGER i,j,k,in_out
 DO i=is,is+1
   DO j=js,js+1
     DO k=ks-1,ks+1
       IF(k/=iz0-1.and.k/=iz1+1.and. &
          j/=iy0  .and.j/=iy1+1.and. &
          i/=ix0  .and.i/=ix1+1) THEN
         IF(UpBoundsLayer(i,j,k)%BCell%bfxynr<=0) THEN
           in_out=Vertices(i-1,j-1,k)%in_out+Vertices(i,j-1,k)%in_out+ &
                  Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
           IF(in_out>0) THEN
             upbound_fnr=upbound_fnr+1
             Flid%sbfxynr=Flid%sbfxynr+1 
             UpBoundsLayer(i,j,k)%BCell%bfxynr=upbound_fnr 
             !write(11,*) "face_xy",i,j,k
             !write(11,*) "                  ",Vertices(i-1,j,k)%in_out,"  ",Vertices(i,j,k)%in_out
             !write(11,*) "                  ",Vertices(i-1,j-1,k)%in_out,"  ",Vertices(i,j-1,k)%in_out
           END IF
         END IF
       END IF
     END DO
   END DO
 END DO
END SUBROUTINE ProveBoundsFaceXY


SUBROUTINE ProveBoundsFaceYZ(is,js,ks,Flid)
 INTEGER is,js,ks
 TYPE (Domain_T) :: Flid
 !local
 !Point wenn Celle%vc
 !FaceYZ: ( k,k+1) | i-1,j+1 | i,j+1 | i+1,j+1 |
 !                 | i-1,j   | i,j   | i+1,j   |
 !      -------------
 !    2/    4/    6/  Reihenfolge jeweils 'i,j,k'
 !     ------*-----   Point '*'
 !   1/    3/   5/
 !  --------------
 !
 INTEGER i,j,k,in_out
 DO i=is-1,is+1
   DO j=js,js+1
     DO k=ks,ks+1
       IF(k/=iz0  .and.k/=iz1+1.and. &
          j/=iy0  .and.j/=iy1+1.and. &
          i/=ix0-1.and.i/=ix1+1) THEN
         IF(UpBoundsLayer(i,j,k)%BCell%bfyznr<=0) THEN
           in_out=Vertices(i,j-1,k-1)%in_out+Vertices(i,j-1,k)%in_out+ &
                  Vertices(i,j,k)%in_out+Vertices(i,j,k-1)%in_out
           IF(in_out>0) THEN
             upbound_fnr=upbound_fnr+1
             Flid%sbfyznr=Flid%sbfyznr+1 
             UpBoundsLayer(i,j,k)%BCell%bfyznr=upbound_fnr
             !write(11,*) "face_yz",i,j,k
             !write(11,*) "                  ",Vertices(i,j-1,k)%in_out,"  ",Vertices(i,j,k)%in_out
             !write(11,*) "                  ",Vertices(i,j-1,k-1)%in_out,"  ",Vertices(i,j,k-1)%in_out
           END IF
         END IF
       END IF 
     END DO
   END DO
 END DO
END SUBROUTINE ProveBoundsFaceYZ

SUBROUTINE ProveBoundsFaceZX(is,js,ks,Flid)
 INTEGER is,js,ks
 TYPE (Domain_T) :: Flid
 !local
 !Point wenn Celle%vc
 !FaceZX: ( k,k+1) | i,j+1 | i+1,j+1 | 
 !                 | i,j   | i+1,j   | 
 !                 | i,j-1 | i+1,j-1 |
 !          ---- ----
 !         | 3  |  6 |  Reihenfolge jeweils 'i,j,k'
 !          ---- ----
 !        ----*----     Point '*'
 !       | 2  |  5 |
 !        ---- ----
 !      ---- ----
 !     | 1  |  4 |
 !      ---- ----
 INTEGER i,j,k,in_out
 DO i=is,is+1
   DO j=js-1,js+1
     DO k=ks,ks+1
       IF(k/=iz0  .and.k/=iz1+1.and. &
          j/=iy0-1.and.j/=iy1+1.and. &
          i/=ix0  .and.i/=ix1+1) THEN
         IF(UpBoundsLayer(i,j,k)%BCell%bfzxnr<=0) THEN
           in_out=Vertices(i-1,j,k-1)%in_out+Vertices(i,j,k-1)%in_out+ &
                  Vertices(i-1,j,k)%in_out+Vertices(i,j,k)%in_out
           IF(in_out>0) THEN
             upbound_fnr=upbound_fnr+1
             Flid%sbfzxnr=Flid%sbfzxnr+1
             UpBoundsLayer(i,j,k)%BCell%bfzxnr=upbound_fnr
             !write(11,*) "face_zx",i,j,k
             !write(11,*) "                  ",Vertices(i-1,j,k)%in_out,"  ",Vertices(i,j,k)%in_out
             !write(11,*) "                  ",Vertices(i-1,j,k-1)%in_out,"  ",Vertices(i,j,k-1)%in_out
           END IF
         END IF
       END IF
     END DO
   END DO
 END DO
END SUBROUTINE ProveBoundsFaceZX

SUBROUTINE ProveBoundEdgX(is,js,ks,bexnr)
  INTEGER :: is,js,ks,bexnr
  INTEGER :: i,j,k,p_vc

  p_vc=0
  DO i=is,is+1
    DO j=js,js+1
      DO k=ks,ks+1
      IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
         IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
            .OR.Cell(i,j,k)%Cell%in_out==8) THEN
          !IF(Cell(i,j,k)%Cell%vc==0) THEN 
          !ELSE ! vc>0 
          !END IF   !beide eingeschossen
          p_vc=p_vc+1 
         END IF 
      END IF
      END DO
    END DO
  END DO
  sbexnr=bexnr 
  IF (p_vc>0) THEN
    DO i=is,is+1
      DO j=js-1,js+1
        DO k=ks-1,ks+1
           IF(k/=iz0-1.and.k/=iz1+1.and. &
              j/=iy0-1.and.j/=iy1+1.and. &
              i/=ix1+1) THEN
             IF(UpBoundsLayer(i,j,k)%BCell%bexnr<=0) THEN
               sbexnr=sbexnr+1
               UpBoundsLayer(i,j,k)%BCell%bexnr=sbexnr
             END IF
           END IF
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE ProveBoundEdgX 

SUBROUTINE ProveBoundEdgY(is,js,ks,beynr)
  INTEGER :: is,js,ks,beynr
  INTEGER :: i,j,k,p_vc

  p_vc=0
  DO i=is,is+1
    DO j=js,js+1
      DO k=ks,ks+1
      IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
         IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
            .OR.Cell(i,j,k)%Cell%in_out==8) THEN
          !IF(Cell(i,j,k)%Cell%vc==0) THEN 
          !ELSE ! vc>0 
          !END IF   !beide eingeschossen
          p_vc=p_vc+1 
         END IF 
      END IF
      END DO
    END DO
  END DO
  sbeynr=beynr
  DO i=is-1,is+1
    DO j=js,js+1
      DO k=ks-1,ks+1
        IF(k/=iz0-1.and.k/=iz1+1.and. &
           j/=iy1+1.and. &
           i/=ix0-1.and.i/=ix1+1) THEN
           IF(UpBoundsLayer(i,j,k)%BCell%beynr<=0) THEN
             sbeynr=sbeynr+1
             UpBoundsLayer(i,j,k)%BCell%beynr=sbeynr 
           END IF
        END IF
      END DO
    END DO
  END DO
 
END SUBROUTINE ProveBoundEdgY

SUBROUTINE ProveBoundEdgZ(is,js,ks,beznr)
  INTEGER :: is,js,ks,beznr
  INTEGER :: i,j,k,p_vc

  p_vc=0
  DO i=is,is+1
    DO j=js,js+1
      DO k=ks,ks+1
      IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
         IF ((Cell(i,j,k)%Cell%in_out>-8.AND.Cell(i,j,k)%Cell%Vol>0.0d0) &
            .OR.Cell(i,j,k)%Cell%in_out==8) THEN
          !IF(Cell(i,j,k)%Cell%vc==0) THEN 
          !ELSE ! vc>0 
          !END IF   !beide eingeschossen
          p_vc=p_vc+1 
         END IF 
      END IF
      END DO
    END DO
  END DO
  sbeznr=beznr
  DO i=is-1,is+1
    DO j=js-1,js+1
      DO k=ks,ks+1
        IF(k/=iz1+1.and. &
           j/=iy0-1.and.j/=iy1+1.and. &
           i/=ix0-1.and.i/=ix1+1) THEN
           IF(UpBoundsLayer(i,j,k)%BCell%beznr<=0) THEN
             sbeznr=sbeznr+1
             UpBoundsLayer(i,j,k)%BCell%beznr=sbeznr
           END IF
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE ProveBoundEdgZ


SUBROUTINE ProveBoundsCells(Flid)
  TYPE (Domain_T) :: Flid

  INTEGER :: i,j,k,in_out
  INTEGER :: ib,is,js,ks
  INTEGER :: p_vc,bcnr,stbound_cnr,ai,aj,defdimz
  !1) Zählung Cut-Cells für Bounds in ProveBoundsCells()
  !2) Ermitteln Bounds, 1 Celle über Point Celle mit vc existiert  
  !  
  !Cell: (k) | i,j+1,k  | i+1,j+1,k |
  !          | i,j  ,k  | i+1,j  ,k |
  !      -------------
  !     /  2  /  4  /   Reihenfolge jeweils 'i,j,k'
  !     -----*------    Point '*'
  !   /  1  /  3  /
  !  --------------
  stbound_cnr=upbound_cnr
  bcnr=0
  DO i=ix0+1,ix1
    DO j=iy0+1,iy1
      DO k=iz0+1,iz1
         p_vc=0 
         DO is=i,i+1
           DO js=j,j+1
             IF (ASSOCIATED(Cell(is,js,k)%Cell)) THEN
                IF ((Cell(is,js,k)%Cell%in_out>-8.AND.Cell(is,js,k)%Cell%Vol>0.0d0) &
                   .OR.Cell(is,js,k)%Cell%in_out==8) THEN
                 !IF(Cell(i,j,k)%Cell%vc==0) THEN 
                 !ELSE ! vc>0 
                 !END IF   !beide eingeschossen
                 p_vc=p_vc+1 
                END IF 
             END IF
           END DO
         END DO
         IF(p_vc>0) THEN
           CALL ProveBoundsEdgeX(i,j,k,Flid)
           CALL ProveBoundsEdgeY(i,j,k,Flid)
           CALL ProveBoundsEdgeZ(i,j,k,Flid)
           CALL ProveBoundsFaceXY(i,j,k,Flid)
           CALL ProveBoundsFaceYZ(i,j,k,Flid)
           CALL ProveBoundsFaceZX(i,j,k,Flid)
           ai=MIN(i+1,ix1)
           aj=MIN(j+1,iy1)
           defdimz=MIN(k+1,iz1)
           DimUpBoundsN(i-1:ai,j-1:aj,2)=defdimz
           DO is=i,i+1
             DO js=j,j+1
               DO ks=k,k+1
                 IF(ks/=iz0.and.ks/=iz1+1.and. &
                    js/=iy0.and.js/=iy1+1.and. &
                    is/=ix0.and.is/=ix1+1) THEN
                   IF (.NOT.ASSOCIATED(Cell(is,js,ks)%Cell)) THEN
                     in_out=Vertices(is-1,js-1,ks-1)%in_out &
                           +Vertices(is-1,js  ,ks-1)%in_out &
                           +Vertices(is  ,js-1,ks-1)%in_out &
                           +Vertices(is  ,js  ,ks-1)%in_out &
                           +Vertices(is-1,js-1,ks  )%in_out &
                           +Vertices(is-1,js  ,ks  )%in_out &
                           +Vertices(is  ,js-1,ks  )%in_out &
                           +Vertices(is  ,js  ,ks  )%in_out
                     IF (in_out>=0) THEN
                       IF (UpBoundsLayer(is,js,ks)%BCell%ub_cnr<=0) THEN
                         bcnr=bcnr+1
                         upbound_cnr=upbound_cnr+1
                         UpBoundsLayer(is,js,ks)%BCell%ub_cnr=upbound_cnr
                       END IF
                     END IF
                   END IF
                 END IF
               END DO ! ks=k,k+1
             END DO ! js=j,j+1
           END DO ! is=i,i+1
         END IF  ! IF(p_vc>0) 
      END DO  ! k=iz0+1,iz1-1
    END DO   ! j=iy0+1,iy1-1
  END DO   ! i=ix0+1,ix1-1
  Flid%sbchnr=bcnr
  Flid%sbcnr=Flid%sbcvcnr+Flid%sbchnr
  Flid%sbfnr=Flid%sbfxynr+Flid%sbfzxnr+Flid%sbfyznr+Flid%sbfcnr
  Flid%sbenr=Flid%sbexnr+Flid%sbeynr+Flid%sbeznr+Flid%sbecnr

END SUBROUTINE ProveBoundsCells

SUBROUTINE ProveBordersBoundsLayer
  INTEGER :: ib,in,i,j,k,ix,jx
  INTEGER :: fac_ra,fac_ne,fac_n
  INTEGER :: st_upe,st_upf,st_upc,zw_cnr,zw_iedg,zw_iface
  !..............................
  upboundgi_enr=upbound_enr
  upboundgi_fnr=upbound_fnr
  upboundgi_cnr=upbound_cnr
  !Zählung:  upbound_* -> aufsteigend der Reihe nach über alle blöcke, mit border
  !          sbe*-edges,sbf*-faces -> je block, binhaltet auch border-zählung
  !                                   border-liste jeweils am ende special gelistet
  !          sbc*                  -> je block, border-zählung-liste zusätzlich
  !                                   border-liste am ende block gelistet
  DO ib=1,nb
    Write(*,*) "                                     Block : ",ib,"\/",nb
    CALL Set(Floor(ib))
  
  st_upe=upbound_enr
  st_upf=upbound_fnr
  st_upc=upbound_cnr
  DO in=1,AnzahlNachbar
    CALL Set(Nachbars(in))
    SELECT CASE (Nachbars(in)%nType(2:2))
      !-------------------------------------------------------------------------
      CASE ("w","e")
      !-------------
        SELECT CASE (Nachbars(in)%nType(1:1))
          ! 'o,p,i' -> west,east  z.Zt: nur Refine =0, Incr->X,Y,Z noch offen  
          !------------------------TypeW=="ow".OR.TypeE=="oe"------------------- 
          CASE ('o')
             ! Spieglung der cellen, damit gleiches Vol-Rand
             IF (nType=="ow") THEN
                ix=ix0        !aktuellRandCelle
                jx=ix0+1      !NachbarEigenCelle
                fac_ra=1      !Anw. minus
                fac_ne=0 
             ELSE IF (nType=="oe") THEN
                ix=ix1+1      !aktuellRandCelle
                jx=ix1        !NachbarEigenCelle
                fac_ra=0
                fac_ne=1      !Anw. minus
             END IF
             zw_cnr=upbound_cnr
             DO j=iy0+1,iy1          !Cellen 'o-we'
               DO k=iz0+1,iz1
                 IF(UpBoundsLayer(jx,j,k)%BCell%ub_cnr>0) THEN
                    sbcnr=sbcnr+1
                    upbound_cnr=upbound_cnr+1
                    UpBoundsLayer(ix,j,k)%BCell%ub_cnr=upbound_cnr
                 END IF
               END DO
             END DO
             IF(nType=="ow") THEN
                  sbcnr_bdr(1)=upbound_cnr-zw_cnr
             ELSE !"oe"
                  sbcnr_bdr(2)=upbound_cnr-zw_cnr
             END IF         
             DO j=iy0,iy1            !EdgeX  'o-we'
               DO k=iz0,iz1
                 IF(UpBoundsLayer(jx,j,k)%BCell%bexnr>0) THEN
                   sbexnr=sbexnr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix,j,k)%BCell%bexnr=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeY  'o-we'
               DO k=iz0,iz1
                 IF(UpBoundsLayer(jx-fac_ne,j,k)%BCell%beynr>0) THEN
                   sbeynr=sbeynr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%beynr=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !EdgeZ  'o-we'
               DO k=iz0+1,iz1
                 IF(UpBoundsLayer(jx-fac_ne,j,k)%BCell%beznr>0) THEN
                   sbeznr=sbeznr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%beznr=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeCut-becnrxy  'o-we'
               DO k=iz0,iz1
                 IF(UpBoundsLayer(jx,j,k)%BCell%becnrxy>0) THEN
                   sbecnr=sbecnr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix,j,k)%BCell%becnrxy=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeCut-becnryz  'o-we'
               DO k=iz0+1,iz1
                 IF(UpBoundsLayer(jx-fac_ne,j,k)%BCell%becnryz>0) THEN
                   sbecnr=sbecnr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%becnryz=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !EdgeCut-becnrzx 'o-we'
               DO k=iz0+1,iz1
                 IF(UpBoundsLayer(jx,j,k)%BCell%becnrzx>0) THEN
                   sbecnr=sbecnr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix,j,k)%BCell%becnrzx=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !FCutNr  'o-we'
               DO k=iz0+1,iz1
                 IF(Topo(jx,j,k)%Topo%FCutNr>0) THEN
                   sbfcnr=sbfcnr+1
                   upbound_fnr=upbound_fnr+1
                   UpBoundsLayer(ix,j,k)%BCell%FCutNr=upbound_fnr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !bfxynr  'o-we'
               DO k=iz0,iz1
                 IF(UpBoundsLayer(jx,j,k)%BCell%bfxynr>0) THEN
                   sbfxynr=sbfxynr+1
                   upbound_fnr=upbound_fnr+1
                   UpBoundsLayer(ix,j,k)%BCell%bfxynr=upbound_fnr
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !bfzxnr  'o-we'
               DO k=iz0+1,iz1
                 IF(UpBoundsLayer(jx,j,k)%BCell%bfzxnr>0) THEN
                   sbfzxnr=sbfzxnr+1
                   upbound_fnr=upbound_fnr+1
                   UpBoundsLayer(ix,j,k)%BCell%bfzxnr=upbound_fnr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !bfyznr  'o-we'
               DO k=iz0+1,iz1
                 IF(UpBoundsLayer(jx-fac_ne,j,k)%BCell%bfyznr>0) THEN
                   sbfyznr=sbfyznr+1
                   upbound_fnr=upbound_fnr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%bfyznr=upbound_fnr
                 END IF
               END DO
             END DO
          !-----------------------TypeW=="pw".OR.TypeE=="pe"---------------------
          CASE ('p')
             IF(nType=="pw") THEN
                ix=ix0               !aktuellRandCelle
                jx=Floor(ibn)%ix1    !NachbarCelle
                fac_ra=1
                fac_n =1
             ELSE !'pe'
                ix=ix1+1             !aktuellRandCelle
                jx=Floor(ibn)%ix0+1  !NachbarCelle
                fac_ra=0
                fac_n =0
             END IF

             zw_cnr=upbound_cnr
             DO j=iy0+1,iy1          !Celle   'p-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%ub_cnr>0) THEN
                   sbcnr=sbcnr+1  
                   upbound_cnr=upbound_cnr+1
                   UpBoundsLayer(ix,j,k)%BCell%ub_cnr=upbound_cnr
                 END IF
               END DO
             END DO
             IF(nType=="pw") THEN
                  sbcnr_bdr(1)=upbound_cnr-zw_cnr
             ELSE !"pe"
                  sbcnr_bdr(2)=upbound_cnr-zw_cnr
             END IF         
             DO j=iy0,iy1            !EdgeX   'p-we'
               DO k=iz0,iz1       
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bexnr>0) THEN
                   sbexnr=sbexnr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix,j,k)%BCell%bexnr=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeY   'p-we'
               DO k=iz0,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beynr>0) THEN
                   sbeynr=sbeynr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%beynr=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !EdgeZ   'p-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beznr>0) THEN
                   sbeznr=sbeznr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%beznr=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeCut-becnrxy  'p-we'
               DO k=iz0,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%becnrxy>0) THEN
                   sbecnr=sbecnr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix,j,k)%BCell%becnrxy=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeCut-becnryz  'p-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%becnryz>0) THEN
                   sbecnr=sbecnr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%becnryz=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !EdgeCut-becnrzx   'p-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%becnrzx>0) THEN
                   sbecnr=sbecnr+1
                   upbound_enr=upbound_enr+1
                   UpBoundsLayer(ix,j,k)%BCell%becnrzx=upbound_enr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !FCutNr   'p-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%Topo(jx,j,k)%Topo%FCutNr>0) THEN
                   sbfcnr=sbfcnr+1
                   upbound_fnr=upbound_fnr+1
                   UpBoundsLayer(ix,j,k)%BCell%FCutNr=upbound_fnr
                   !UpBoundsLayer(ix,j,k)%BCell%fceg_nr= &
                   !  & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%fceg_nr
                   !UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(:)= &
                   !  & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%fceg_nrcut(:)
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !bfxynr   'p-we'
               DO k=iz0,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bfxynr>0) THEN
                   sbfxynr=sbfxynr+1
                   upbound_fnr=upbound_fnr+1
                   UpBoundsLayer(ix,j,k)%BCell%bfxynr=upbound_fnr
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !bfzxnr   'p-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bfzxnr>0) THEN
                   sbfzxnr=sbfzxnr+1
                   upbound_fnr=upbound_fnr+1
                   UpBoundsLayer(ix,j,k)%BCell%bfzxnr=upbound_fnr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !bfyznr   'p-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%bfyznr>0) THEN
                   sbfyznr=sbfyznr+1
                   upbound_fnr=upbound_fnr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%bfyznr=upbound_fnr
                 END IF
               END DO
             END DO
          !-----------------------TypeW=="iw".OR.TypeE=="ie"---------------------
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
             zw_iedg=0;zw_iface=0;
             zw_cnr=sbcnr
             DO j=iy0+1,iy1          !Celle   'i_we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%ub_cnr>0) THEN
                   sbcnr=sbcnr+1
                   UpBoundsLayer(ix,j,k)%BCell%ub_cnr= &
                    & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%ub_cnr
                 END IF
               END DO
             END DO
             IF(nType=="iw") THEN
                  sbcnr_bdr(1)=sbcnr-zw_cnr
             ELSE !"ie"
                  sbcnr_bdr(2)=sbcnr-zw_cnr
             END IF         
             DO j=iy0,iy1            !EdgeX   'i-we'
               DO k=iz0,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bexnr>0) THEN
                   sbexnr=sbexnr+1  !zw_iedg=zw_iedg+1 wenn extra gezählt
                   UpBoundsLayer(ix,j,k)%BCell%bexnr= &
                     & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bexnr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeY  'i_we'
               DO k=iz0,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beynr>0) THEN
                   sbeynr=sbeynr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%beynr= &
                     & Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beynr
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !EdgeZ   'i-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beznr>0) THEN
                   sbeznr=sbeznr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%beznr= &
                    & Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%beznr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeCut-becnrxy   'i-we'
               DO k=iz0,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%becnrxy>0) THEN
                   sbecnr=sbecnr+1
                   UpBoundsLayer(ix,j,k)%BCell%becnrxy= &
                    & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%becnrxy
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !EdgeCut-becnryz    'i-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%becnryz>0) THEN
                   sbecnr=sbecnr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%becnryz= &
                     & Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%becnryz
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !EdgeCut-becnrzx   'i-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%becnrzx>0) THEN
                   sbecnr=sbecnr+1
                   UpBoundsLayer(ix,j,k)%BCell%becnrzx= &
                    & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%becnrzx
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !FCutNr    'i-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%Topo(jx,j,k)%Topo%FCutNr>0) THEN
                   sbfcnr=sbfcnr+1
                   UpBoundsLayer(ix,j,k)%BCell%FCutNr= &
                     & Floor(ibn)%Topo(jx,j,k)%Topo%FCutNr
                   UpBoundsLayer(ix,j,k)%BCell%fceg_nr= &
                     & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%fceg_nr
                   UpBoundsLayer(ix,j,k)%BCell%fceg_nrcut(:)= &
                     & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%fceg_nrcut(:)
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !bfxynr   'i-we'
               DO k=iz0,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bfxynr>0) THEN
                   sbfxynr=sbfxynr+1
                   UpBoundsLayer(ix,j,k)%BCell%bfxynr= &
                    & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bfxynr
                 END IF
               END DO
             END DO
             DO j=iy0,iy1            !bfzxnr  'i-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bfzxnr>0) THEN
                   sbfzxnr=sbfzxnr+1
                   UpBoundsLayer(ix,j,k)%BCell%bfzxnr= &
                    & Floor(ibn)%UpBoundsLayer(jx,j,k)%BCell%bfzxnr
                 END IF
               END DO
             END DO
             DO j=iy0+1,iy1          !bfyznr   'i-we'
               DO k=iz0+1,iz1
                 IF(Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%bfyznr>0) THEN
                   sbfyznr=sbfyznr+1
                   UpBoundsLayer(ix-fac_ra,j,k)%BCell%bfyznr= &
                    &  Floor(ibn)%UpBoundsLayer(jx-fac_n,j,k)%BCell%bfyznr
                 END IF
               END DO
             END DO
        END SELECT   ! [o,p,i]-> 'w/e'

      !--------------------------------------------------------------------
      CASE ("n","s")
      !-------------
        SELECT CASE (Nachbars(in)%nType(1:1))
          ! 'o,p,i' -> north,south
          CASE ('o')
            IF (nType=="on") THEN
            ELSE IF (nType=="os") THEN
            END IF 
          CASE ('p')
            IF (nType=="pn") THEN
            ELSE IF (nType=="ps") THEN
            END IF 
          CASE ('i')
            IF (nType=="in") THEN
            ELSE IF (nType=="is") THEN
            END IF 
        END SELECT   ! [o,p,i]-> 'n/s'

      !--------------------------------------------------------------------
      CASE ("t","b")
      !-------------
        SELECT CASE (Nachbars(in)%nType(1:1))
          ! 'o,p,i' -> top, bottom
          CASE ('o')
             IF (nType=="ot") THEN
             ELSE IF (nType=="ob") THEN
             END IF 
          CASE ('p')
             IF (nType=="pt") THEN
             ELSE IF (nType=="pb") THEN
             END IF 
          CASE ('i')
             IF (nType=="it") THEN
             ELSE IF (nType=="ib") THEN
             END IF 
        END SELECT   ! [o,p,i]-> 't,b'

    END SELECT   ! 'w/e','n/s','t/b' -> to boundslayer
  END DO  ! AnzahlNachbar
  Floor(ib)%sbexnr=sbexnr
  Floor(ib)%sbeynr=sbeynr
  Floor(ib)%sbeznr=sbeznr
  Floor(ib)%sbecnr=sbecnr
  Floor(ib)%sbenr=sbexnr+sbeynr+sbeznr+sbecnr
  !Floor(ib)%sbenr=sbenr+upbound_enr-st_upe+zw_iedg
  Floor(ib)%sbfxynr=sbfxynr
  Floor(ib)%sbfzxnr=sbfzxnr
  Floor(ib)%sbfyznr=sbfyznr
  Floor(ib)%sbfcnr=sbfcnr
  Floor(ib)%sbfnr=sbfxynr+sbfzxnr+sbfyznr+sbfcnr
  !Floor(ib)%sbfnr=sbfnr+upbound_fnr-st_upf+zw_iface
  Floor(ib)%sbcnr=sbcnr
  Floor(ib)%sbcnr_bdr(:)=sbcnr_bdr(:)
 
  END DO  !ib

END SUBROUTINE ProveBordersBoundsLayer

SUBROUTINE InitTropoBoundOroStruct(FileName)
  CHARACTER*50 :: FileName
  INTEGER :: ib

  !-------------------------------------------
  ! Tropo/BoundsLayer/Oro - Struct
  !-------------------------------------------
  CALL Display_InitTropoOroOutStruct(FileName)
  CALL InitAllTopo
  CALL InitUpBoundsLayer

  DO ib=1,nb
    Write(*,*) "                                     Block : ",ib,"\/",nb
    CALL Set(Floor(ib))
    !............................................
    CALL ProveTropoCells(Floor(ib))
    CALL ProveTOroCells(Floor(ib))
    !.........................................
    CALL ProveDimUpBounds(Floor(ib))
    !.........................................
    CALL ProveTropoEdges(Floor(ib))
    CALL ProveTOroEdges(Floor(ib))
    !.........................................
    CALL ProveTropoFaces(Floor(ib))
    CALL ProveTOroFaces(Floor(ib))
    CALL ProveBoundsCells(Floor(ib))
  END DO  !ib

  !CALL Display_InitBoundsBorder_AllDirection(FileName)

  CALL ProveBordersBoundsLayer

!     Write(*,*) "z.Zt.  InitTropoOroOutStruct ","  nur wenn  '#OutputDomain'  gesetzt ist!"  
!    WRITE(OutUnitProt,'(a8,i8)') cells,nr_cells
!    DO ib=1,nb
!      Write(*,*) "                                     Block : ",ib,"\/",nb
!      CALL Set(Floor(ib))
!      DO k=iz0+1,iz1
!        DO j=iy0+1,iy1
!          DO i=ix0+1,ix1
!            !CALL WriteCellGMVAscii(Cell(i,j,k)%Cell,i,j,k)
!            IF (ASSOCIATED(Cell(i,j,k)%Cell)) THEN
!              Write(OutUnitProt,*) "...................................."
!              Write(OutUnitProt,*) "Cell(",i,",",j,",",k,")%Cell% ->"
!              Write(OutUnitProt,*) "  Face1%mp=",Cell(i,j,k)%Cell%Face1%mp
!              Write(OutUnitProt,*) "  Face2%mp=",Cell(i,j,k)%Cell%Face2%mp
!              Write(OutUnitProt,*) "  Face3%mp=",Cell(i,j,k)%Cell%Face3%mp
!              Write(OutUnitProt,*) "  Face4%mp=",Cell(i,j,k)%Cell%Face4%mp
!              Write(OutUnitProt,*) "  Face5%mp=",Cell(i,j,k)%Cell%Face5%mp
!              Write(OutUnitProt,*) "  Face6%mp=",Cell(i,j,k)%Cell%Face6%mp
!              Write(OutUnitProt,*) "  vc=",Cell(i,j,k)%Cell%vc
!              WRITE(OutUnitProt,*) (Cell(i,j,k)%Cell%VertCut(ic),ic=1,Cell(i,j,k)%Cell%vc)
!              Write(OutUnitProt,*) "...................................."
!            END IF
!          END DO
!        END DO
!      END DO
!    END DO  !ib

END SUBROUTINE InitTropoBoundOroStruct

END MODULE InitTropoOroStruct_Mod
