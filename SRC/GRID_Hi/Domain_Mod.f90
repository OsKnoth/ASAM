MODULE Domain_Mod

  USE Kind_Mod
  USE Geometry_Mod

  IMPLICIT NONE

  TYPE Vertex_T
    TYPE (Point_T) :: Point
    INTEGER :: in_out
    INTEGER :: nrP=0
    INTEGER :: nrInP=0
    INTEGER :: nrCutP=0
    INTEGER :: Shift=1
  END TYPE

  TYPE Edge_T
    TYPE (Vertex_T), POINTER :: Vert1=>NULL()
    TYPE (Vertex_T), POINTER :: Vert2=>NULL()
    INTEGER :: in_out
    INTEGER :: yes_sp=-1
    TYPE (Vertex_T), POINTER :: VertS=>NULL()
  END TYPE

  TYPE EdgeP_T
    TYPE(Edge_T), POINTER :: Edge=>NULL()
  END TYPE

  TYPE Face_T
    TYPE (Edge_T), POINTER :: Edge1=>NULL()
    TYPE (Edge_T), POINTER :: Edge2=>NULL()
    TYPE (Edge_T), POINTER :: Edge3=>NULL()
    TYPE (Edge_T), POINTER :: Edge4=>NULL()
    INTEGER :: in_out
    INTEGER :: NumberVert
    INTEGER :: VertexList(1:8)
    INTEGER :: ec
    INTEGER :: mp=0
    INTEGER :: EdgeCut(1:2)
    REAL(8) :: Vol=-1111111.111
    TYPE (Point_T) :: MidPoint
  END TYPE

  TYPE FaceP_T
    TYPE(Face_T), POINTER :: Face=>NULL()
  END TYPE

  TYPE Cell_T
    TYPE (Face_T), POINTER :: Face1=>NULL()
    TYPE (Face_T), POINTER :: Face2=>NULL()
    TYPE (Face_T), POINTER :: Face3=>NULL()
    TYPE (Face_T), POINTER :: Face4=>NULL()
    TYPE (Face_T), POINTER :: Face5=>NULL()
    TYPE (Face_T), POINTER :: Face6=>NULL()
    INTEGER :: in_out
    INTEGER :: vc
    INTEGER :: mp=0
    INTEGER, POINTER :: VertCut(:)=>NULL()
    INTEGER, POINTER :: VertCutOnlyP(:)=>NULL()
    REAL(8) :: Vol=-1111111.111
    TYPE (Point_T) :: MidPoint
    TYPE (Point_T) :: CutF_MidP
    INTEGER, POINTER :: layer_soil(:)=>NULL()
    INTEGER :: LandClass=0
  END TYPE

  TYPE CellP_T
    TYPE(Cell_T), POINTER :: Cell=>NULL()
  END TYPE

  TYPE FaceP4_T
    TYPE (Vertex_T) :: V0,V1,V2,V3
    TYPE (Edge_T) :: Edge1
    TYPE (Edge_T) :: Edge2
    TYPE (Edge_T) :: Edge3
    TYPE (Edge_T) :: Edge4
    INTEGER :: in_out
    INTEGER :: NumberVert
    INTEGER :: VertexList(1:8)
    INTEGER :: ec
    INTEGER :: mp=0
    INTEGER :: EdgeCut(1:2)
    REAL(8) :: Vol=-1111111.111
    TYPE (Point_T) :: MidPoint
  END TYPE

  TYPE CellP8_T
    TYPE (Vertex_T) :: V0,V1,V2,V3,V4,V5,V6,V7
    TYPE (Face_T) :: Face1
    TYPE (Face_T) :: Face2
    TYPE (Face_T) :: Face3
    TYPE (Face_T) :: Face4
    TYPE (Face_T) :: Face5
    TYPE (Face_T) :: Face6
    INTEGER :: in_out    
    INTEGER :: vc
    INTEGER :: mp=0
    REAL(8) :: Vol=-1111111.111
    TYPE (Point_T) :: MidPoint
  END TYPE
!.....................................................  
  CHARACTER(8), PARAMETER :: gmvinput ='gmvinput'
  CHARACTER(8), PARAMETER :: ascii    ='   ascii'
  CHARACTER(8), PARAMETER :: ieeei4r4 ='ieeei4r4'
  CHARACTER(8), PARAMETER :: ieeei4r8 ='ieeei4r8'
  CHARACTER(8), PARAMETER :: iecxi4r4 ='iecxi4r4'
  CHARACTER(8), PARAMETER :: iecxi4r8 ='iecxi4r8'
  CHARACTER(8), PARAMETER :: nodes    ='nodes   '
  CHARACTER(8), PARAMETER :: nodev    ='nodev   '
  CHARACTER(8), PARAMETER :: cells    ='cells   '
  CHARACTER(8), PARAMETER :: material ='material'
  CHARACTER(8), PARAMETER :: velocity ='velocity'
  CHARACTER(8), PARAMETER :: hex      ='hex     '
  CHARACTER(8), PARAMETER :: general  ='general '
  CHARACTER(8), PARAMETER :: polygons ='polygons'
  CHARACTER(8), PARAMETER :: endpoly  ='endpoly '
  CHARACTER(8), PARAMETER :: endflag  ='endflag '
  CHARACTER(8), PARAMETER :: fromfile ='fromfile'
  CHARACTER(8), PARAMETER :: endgmv   ='endgmv  '
  !..................................................
  CHARACTER(10), PARAMETER :: def_out_index='XYZ_Number' !'Gr_Index'
  CHARACTER( 9), PARAMETER :: def_out_coord='XYZ_Coord'  !'Gr_Coord'
  !..................................................
  TYPE Boundary_T
    CHARACTER(10) :: West=''
    CHARACTER(10) :: East=''
    CHARACTER(10) :: South=''
    CHARACTER(10) :: North=''
    CHARACTER(10) :: Bottom=''
    CHARACTER(10) :: Top=''
  END TYPE Boundary_T

  TYPE(Boundary_T), SAVE :: BC

  TYPE Nachbar_T
    CHARACTER(2) :: nTYPE
    INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1    ! Coordinates of border in block coordinates
    INTEGER :: ib
    INTEGER :: Refine
    INTEGER :: RefineX
    INTEGER :: RefineY
    INTEGER :: RefineZ
    INTEGER :: IncrX
    INTEGER :: IncrY
    INTEGER :: IncrZ
    INTEGER :: CopyCase
  END TYPE Nachbar_T

  TYPE Domain_T
    CHARACTER(2) :: TypeE,TypeW,TypeS,TypeN,TypeB,TypeT
    REAL(RealKind) :: Boundary
    INTEGER :: nx,ny,nz            ! numbers def. x-,y-,z-direction the domain-grid
    INTEGER :: nc                  ! number of cells the domain-grid
    INTEGER :: nr_out              ! number points output for gmv 
    INTEGER :: WriteOffsetC        ! Cell Offset for parallel I/O
    INTEGER :: WriteOffsetN        ! Node Offset for parallel I/O
    INTEGER :: igx0,igx1,igy0,igy1,igz0,igz1
    INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1,ncz
    !......................................................................
    !Output-Area-ViewGMV
    CHARACTER(10)   :: def_out_domain
    INTEGER        :: view_ixa,view_ixe,view_iya,view_iye,view_iza,view_ize
    REAL(RealKind) :: view_xa,view_xe,view_ya,view_ye,view_za,view_ze
    !......................................................................
    INTEGER :: ib
    INTEGER :: RefLevel
    INTEGER :: Refine
    INTEGER :: RefineX
    INTEGER :: RefineY
    INTEGER :: RefineZ
    INTEGER :: Xshift
    INTEGER :: Yshift
    INTEGER :: Zshift
    REAL(RealKind) :: x0,x1,y0,y1,z0,z1,calc_z
    REAL(RealKind) :: x0View,x1View,y0View,y1View,z0View,z1View !Domain
    REAL(RealKind), POINTER :: dx(:),dy(:),dz(:)
    REAL(RealKind), POINTER :: xP(:),yP(:),zP(:)
    TYPE(Vertex_T), POINTER :: Vertices(:,:,:)
    TYPE (EdgeP_T), POINTER :: Edges_X(:,:,:)
    TYPE (EdgeP_T), POINTER :: Edges_Y(:,:,:)
    TYPE (EdgeP_T), POINTER :: Edges_Z(:,:,:)
    TYPE (FaceP_T), POINTER :: Faces_XY(:,:,:)
    TYPE (FaceP_T), POINTER :: Faces_YZ(:,:,:)
    TYPE (FaceP_T), POINTER :: Faces_ZX(:,:,:)
    TYPE (CellP_T), POINTER :: Cell(:,:,:)
    INTEGER          :: nr_soildef
    INTEGER, POINTER :: soil_type(:,:)
    INTEGER, POINTER :: s_ixa(:),s_ixe(:),s_iya(:),s_iye(:),s_iza(:),s_ize(:)
    REAL(RealKind), POINTER :: WeiFU(:,:,:),WeiFV(:,:,:),WeiFW(:,:,:)
    REAL(RealKind), POINTER :: VolC(:,:,:)
    REAL(RealKind), POINTER :: Emiss(:,:,:,:),Emiss0(:,:,:,:)
    REAL(RealKind), POINTER :: Mask(:,:,:),Mask0(:,:,:)
    INTEGER, POINTER :: Nemiss(:,:,:)

    INTEGER :: NrW_Cells,NrR_Cells,NrRN_Cells,NrRW_Cells
    !INTEGER :: NrW_Cells1,NrW_Cells2
    INTEGER :: NrMP_Cells,NrB_Cells
    INTEGER :: NrW_FacesXY,NrW_FacesYZ,NrW_FacesZX
    INTEGER :: NrR_FacesXY,NrR_FacesYZ,NrR_FacesZX
    INTEGER :: NrRN_FacesXY,NrRN_FacesYZ,NrRN_FacesZX
    INTEGER :: NrRW_FacesXY,NrRW_FacesYZ,NrRW_FacesZX
    INTEGER :: NrMP_FacesXY,NrMP_FacesYZ,NrMP_FacesZX,NrMP_F
    !INTEGER :: sub_face_gr(1:6)  !Index Grenze Face subMountain
    INTEGER :: AnzahlNachbar
    TYPE (Nachbar_T), POINTER :: Nachbars(:)
  END TYPE Domain_T
  
  TYPE (Domain_T) :: Domain
  INTEGER :: nx,ny,nz 
  INTEGER :: WriteOffsetC
  INTEGER :: WriteOffsetN
  INTEGER :: igx0,igx1,igy0,igy1,igz0,igz1
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1,ncz
  INTEGER :: jx0,jx1,jy0,jy1,jz0,jz1
  REAL(RealKind) :: calc_z
  REAL(RealKind), POINTER :: dx(:),dy(:),dz(:)
  REAL(RealKind), POINTER :: xP(:),yP(:),zP(:)
  TYPE(Vertex_T), POINTER :: Vertices(:,:,:)
  TYPE (EdgeP_T), POINTER :: Edges_X(:,:,:)
  TYPE (EdgeP_T), POINTER :: Edges_Y(:,:,:)
  TYPE (EdgeP_T), POINTER :: Edges_Z(:,:,:)

  TYPE (FaceP_T), POINTER :: Faces_XY(:,:,:)
  TYPE (FaceP_T), POINTER :: Faces_YZ(:,:,:)
  TYPE (FaceP_T), POINTER :: Faces_ZX(:,:,:)
 
  TYPE (CellP_T), POINTER :: Cell(:,:,:)

  REAL(RealKind), POINTER :: FU(:,:,:),FV(:,:,:),FW(:,:,:)
  REAL(RealKind), POINTER :: VolC(:,:,:)

  REAL(RealKind), POINTER :: Emiss(:,:,:,:),Emiss0(:,:,:,:)
  REAL(RealKind), POINTER :: Mask(:,:,:),Mask0(:,:,:)
  INTEGER, POINTER :: Nemiss(:,:,:)
  INTEGER :: SpeciesNum=1

  INTEGER :: nr_sb, struct_bound, std_bound=3
  INTEGER :: nr_lsoil=7   !akt. def. Bodenschichten Cell-Soil-Darstellung
  INTEGER :: nr_soildef
  REAL(RealKind) :: dzi_soil
  INTEGER, POINTER :: soil_type(:,:)
  INTEGER, POINTER :: s_ixa(:),s_ixe(:),s_iya(:),s_iye(:),s_iza(:),s_ize(:)
  
  INTEGER :: nr_landdef
  INTEGER, POINTER :: LandDef(:,:,:)

  INTEGER :: AnzahlNachbar
  TYPE (Nachbar_T), POINTER :: Nachbars(:)
  TYPE(Nachbar_T), POINTER :: Nachbar
  INTEGER :: RefLevel
  INTEGER :: Refine
  INTEGER :: RefineX
  INTEGER :: RefineY
  INTEGER :: RefineZ
  REAL(RealKind) :: Boundary
  CHARACTER(2) :: TypeE,TypeW,TypeS,TypeN,TypeB,TypeT
  
  INTEGER :: ibc,ibn
  INTEGER :: RefineNachbar
  INTEGER :: RefineNachbarX
  INTEGER :: RefineNachbarY
  INTEGER :: RefineNachbarZ
  INTEGER :: IncrX,IncrY,IncrZ
  INTEGER :: CopyCase
  CHARACTER(2) :: nType

  INTEGER :: nr_out
  TYPE (Vertex_T), POINTER :: VertOut(:)=>NULL()
  INTEGER :: nr_cutplane
  TYPE (Vertex_T), POINTER :: VertCutPOut(:)=>NULL()
  INTEGER , POINTER        :: VertNrCutPOut(:)=>NULL()
  INTEGER :: nr_out_soil
  TYPE (Vertex_T), POINTER :: VertSoilOut(:)=>NULL()

  INTEGER :: nr_inside
  TYPE (Vertex_T), POINTER :: VertIn(:)=>NULL()

  !Cells counter for GMV
  INTEGER :: nr_cutinside
  INTEGER :: nr_cells
  INTEGER :: nr_viewcells                              ! outside Topogr.
  INTEGER :: nr_cutplanecells,nr_viewcutplanecells     ! cut Topogr.
  INTEGER :: nr_cutplanecells1,nr_viewcutplanecells1   ! only 'cut' Topogr.
  INTEGER :: nr_cutplanecells2,nr_viewcutplanecells2   ! only Celle k==1 'cut' Topogr.
  INTEGER :: nr_soilplanecells,nr_viewsoilplanecells   ! cut=soilplane Topogr.
  INTEGER :: nr_soilplanecells1,nr_viewsoilplanecells1 !  
  INTEGER :: nr_soilplanecells2,nr_viewsoilplanecells2 ! only Celle k==1 'soilplane'
  INTEGER :: nr_cutcells,nr_cutf1,nr_incells,nr_orocells        ! inside Topogr.
  INTEGER :: nr_viewcutcells,nr_viewcutf1,nr_viewincells,nr_vieworocells
  ! into  Oro.out
  INTEGER :: nr_gen=0
  INTEGER :: nr_grenzeF2=0 
  INTEGER :: nr_grenzeF1=0
  INTEGER :: nr_grenzeSp3=0
  INTEGER :: nr_insidehex=0
  INTEGER :: nr_cell_vc_oro=0
  INTEGER :: nr_cell_novc_oro=0
  INTEGER :: nr_ges_oro_out=0

  !Cells counter for Weight
  INTEGER :: NrAll_Cells,NrAll_FacesXY,NrAll_FacesYZ,NrAll_FacesZX
  INTEGER :: NrW_Cells,NrW_FacesXY,NrW_FacesYZ,NrW_FacesZX
  INTEGER :: NrRW_Cells,NrRW_FacesXY,NrRW_FacesYZ,NrRW_FacesZX
  INTEGER :: NrW_All_Cells,NrW_All_FXY,NrW_All_FYZ,NrW_All_FZX
  INTEGER :: NrW_Cells1=0
  INTEGER :: NrW_Cells2=0
  INTEGER :: NrMP_Cells,NrMP_FacesXY,NrMP_FacesYZ,NrMP_FacesZX,NrMP_F
  INTEGER :: NrMP_All_Cells,NrMP_All_FXY,NrMP_All_FYZ,NrMP_All_FZX
  INTEGER :: NrB_Cells,NrB_All_Cells
  INTEGER :: NrR_Cells,NrR_FacesXY,NrR_FacesYZ,NrR_FacesZX
  INTEGER :: NrRN_Cells,NrRN_FacesXY,NrRN_FacesYZ,NrRN_FacesZX
  INTEGER :: sub_cell_gr,sub_face_gr
  !INTEGER :: sub_face_gr(1:6)  !Index Grenze Face subMountain
  !Face Output Weight
  CHARACTER(6), DIMENSION(6):: FName_MP=(/"MP_FB","MP_FT",&
                                          "MP_FS","MP_FN",&
                                          "MP_FW","MP_FE"/)
  CHARACTER(6) :: CName_MP="C_MP"
  TYPE CF_MP_T
    CHARACTER(6)   :: FN_MP
    TYPE (Point_T) :: F_MP
  END TYPE CF_MP_T 
  TYPE(CF_MP_T) :: CF_MP(1:6)
  INTEGER :: nfmp

  !
  CHARACTER(8) :: OutputType
  INTEGER :: nr_wahlfkt,nr_wahlemi,nr_wahlmask
  INTEGER :: v_x0,v_x1,v_y0,v_y1,v_z0,v_z1,v_zw
  CHARACTER :: out_area,view_cell,view_block

  REAL(4) :: tempw
  INTEGER :: nRec
  INTEGER, PARAMETER :: SizeOfReal=4
  
  REAL(8) :: zRauh=0.01d0
  REAL(8) :: Albedo=0.0d0
  REAL(8) :: Emissivity=0.0d0

  !Glob.  Distance Point(x,y,z) coefficient
  REAL(8) :: dxLoc           ! (Domain%x1-Domain%x0)/Domain%nx*distx_coeff 
  REAL(8) :: dyLoc           ! (Domain%y1-Domain%y0)/Domain%ny*disty_coeff
  REAL(8) :: dzLoc           ! (Domain%z1-Domain%z0)/Domain%nz*distz_coeff
  !Input GridDistsCtrl
  REAL(8) :: dist_fscv=1.0d-12     ! for fine scaling dist to in_out-def, CheckVertex 
  REAL(8) :: distx_coeff=0.01      ! Distance Point(x) coefficient
  REAL(8) :: disty_coeff=0.01      ! Distance Point(y) coefficient
  REAL(8) :: distz_coeff=0.01      ! Distance Point(z) coefficient
  REAL(8) :: dxViewLoc=1.0d-8      ! Distance Point(x) coefficient of model border
  REAL(8) :: dyViewLoc=1.0d-8      ! Distance Point(y) coefficient of model border
  REAL(8) :: dzViewLoc=1.0d-8      ! Distance Point(z) coefficient of model border
  REAL(8) :: Shrink=1.0d-6         ! shrinks points of volume-analysis (0 für Einheitswürfel)
  INTEGER :: IncrVol=10            ! counter for cuts of volume-analysis (1 für Einheitswürfel 
  REAL(8) :: dist_scMaxCell=1.0d-12 ! adjustment value to the filter(screen) of cells 
                                    ! with roughly max. Vol


  INTERFACE Set
    MODULE PROCEDURE SetDomain,SetNachbar
    !MODULE PROCEDURE SetDomain,SetNachbar,Set_EnvBlkFace
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE DomainAllocate
  END INTERFACE
  INTERFACE DeAllocate
    MODULE PROCEDURE DomainDeAllocate
  END INTERFACE

CONTAINS

SUBROUTINE WriteFace(Face)
  TYPE(Face_T), POINTER :: Face
  WRITE(*,*) ' Face%--->'
  WRITE(*,*) "    !(TYPE Vertex_T:", "  Point(x,y,z)    in_out   nrP   nrInP   nrCutP   Shift)"
  WRITE(*,*) "    Edge1%Vert1:",Face%Edge1%Vert1 
  WRITE(*,*) "    Edge1%Vert2:",Face%Edge1%Vert2
  WRITE(*,*) "    Edge3%Vert1:",Face%Edge3%Vert1
  WRITE(*,*) "    Edge3%Vert2:",Face%Edge3%Vert2
  WRITE(*,*) "    ..........."
  WRITE(*,*) '    Vol        =',Face%Vol
  WRITE(*,*) '    NumberVert =',Face%NumberVert
  WRITE(*,*) '    VertexList =',Face%VertexList(1:Face%NumberVert)
  WRITE(*,*) '    in_out     =',Face%in_out
  WRITE(*,*) '    ec         =',Face%ec
  WRITE(*,*) '    EdgeCut    =',Face%EdgeCut 
END SUBROUTINE WriteFace

SUBROUTINE WriteCell(Cell,i,j,k)

  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  INTEGER :: F1V1,F1V2,F1V3,F1V4
  INTEGER :: F2V1,F2V2,F2V3,F2V4
  !.............................. 
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

    WRITE(*,*) '-------------- WriteCell() ---------------------------------'
    WRITE(*,*) '............'
    WRITE(*,*) "        ",F2V4,"-----------",F2V3
    WRITE(*,*) "        /|  F2      /|" 
    WRITE(*,*) "      ",F2V1,"-----------",F2V2," |"  
    WRITE(*,*) "       | ",F1V4,"---------|-",F1V3
    WRITE(*,*) "       |/   F1     |/" 
    WRITE(*,*) "      ",F1V1,"-----------",F1V2   
    WRITE(*,*) '* Cell--> ',' \/',i,'\/',j,'\/',k,'\/'
    WRITE(*,*) '    Vol       =',Cell%Vol
    WRITE(*,*) '    MidPoint  =',Cell%MidPoint
    WRITE(*,*) '    vc        =',Cell%vc
    WRITE(*,*) '    in_out    =',Cell%in_out
    WRITE(*,*) '    EdgeCut   =',Cell%VertCut
    WRITE(*,*) '    CutF_MidP =',Cell%CutF_MidP
    WRITE(*,*) '    Face1%Vol =',Cell%Face1%Vol 
    WRITE(*,*) '    Face2%Vol =',Cell%Face2%Vol 
    WRITE(*,*) '    Face3%Vol =',Cell%Face3%Vol 
    WRITE(*,*) '    Face4%Vol =',Cell%Face4%Vol 
    WRITE(*,*) '    Face5%Vol =',Cell%Face5%Vol 
    WRITE(*,*) '    Face6%Vol =',Cell%Face6%Vol 
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceXY','  (F1)'
    CALL WriteFace(Cell%Face1)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceXY','  (F2)'
    CALL WriteFace(Cell%Face2)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceXZ','  (F3)'
    CALL WriteFace(Cell%Face3)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceXZ','  (F4)'
    CALL WriteFace(Cell%Face4)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceYZ','  (F5)'
    CALL WriteFace(Cell%Face5)
    WRITE(*,*) '............'
    WRITE(*,*) ' --> FaceYZ','  (F6)'
    CALL WriteFace(Cell%Face6)
    WRITE(*,*) 'Cell','(',i,',',j,',',k,')'
    WRITE(*,*) '-------------- WriteCell() Ende ----------------------------'
  END IF

END SUBROUTINE WriteCell


SUBROUTINE CheckWriteCellFace(Cell)
  TYPE(Cell_T), POINTER :: Cell
  nfmp=0
  CName_MP="MP_C"
  IF(Cell%Face1%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(1)
     CF_MP(nfmp)%F_MP=Cell%Face1%MidPoint
  END IF  
  IF(Cell%Face2%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(2)
     CF_MP(nfmp)%F_MP=Cell%Face2%MidPoint
  END IF  
  IF(Cell%Face3%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(3)
     CF_MP(nfmp)%F_MP=Cell%Face3%MidPoint
  END IF  
  IF(Cell%Face4%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(4)
     CF_MP(nfmp)%F_MP=Cell%Face4%MidPoint
  END IF  
  IF(Cell%Face5%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(5)
     CF_MP(nfmp)%F_MP=Cell%Face5%MidPoint
  END IF  
  IF(Cell%Face6%mp>0) THEN
     nfmp=nfmp+1
     CF_MP(nfmp)%FN_MP=FName_MP(6)
     CF_MP(nfmp)%F_MP=Cell%Face6%MidPoint
  END IF  
END SUBROUTINE CheckWriteCellFace

SUBROUTINE Display_CellVol(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  Write(*,*)
  Write(*,*) "Check :"
  Write(*,*) "Cell%Vol "," ix=",i," iy=",j," iz=",k  
  IF (ASSOCIATED(Cell)) THEN
     Write(*,*) "Cell%Vol=",Cell%Vol
  ELSE
    Write(*,*) "Celle not associated"
  END IF
END SUBROUTINE Display_CellVol


SUBROUTINE Display_Cell_vc(Cell,i,j,k)
  TYPE(Cell_T), POINTER :: Cell
  INTEGER :: i,j,k

  Write(*,*)
  Write(*,*) "Check :"
  IF (ASSOCIATED(Cell)) THEN
    Write(*,*) "Cell", "(", i, ",", j, ",", k, ")", &
             & " Cell%vc=",Cell%vc," Cell%in_out=",Cell%in_out  
  ELSE
    Write(*,*) "Cell", "(", i, ",", j, ",", k, ")", & 
             & " Celle not associated "
  END IF
END SUBROUTINE Display_Cell_vc


SUBROUTINE SetDomain(Domain)

  TYPE (Domain_T) :: Domain

  nx=Domain%nx
  ny=Domain%ny
  nz=Domain%nz
  WriteOffsetC=Domain%WriteOffsetC
  WriteOffsetN=Domain%WriteOffsetN
  ix0=Domain%ix0
  ix1=Domain%ix1
  iy0=Domain%iy0
  iy1=Domain%iy1
  iz0=Domain%iz0
  iz1=Domain%iz1
  igx0=Domain%igx0
  igx1=Domain%igx1
  igy0=Domain%igy0
  igy1=Domain%igy1
  igz0=Domain%igz0
  igz1=Domain%igz1
  dx=>Domain%dx
  dy=>Domain%dy
  dz=>Domain%dz
  xP=>Domain%xP
  yP=>Domain%yP
  zP=>Domain%zP
  !.......................
  Vertices=>Domain%Vertices
  Edges_X=>Domain%Edges_X
  Edges_Y=>Domain%Edges_Y
  Edges_Z=>Domain%Edges_Z
  Faces_XY=>Domain%Faces_XY
  Faces_YZ=>Domain%Faces_YZ
  Faces_ZX=>Domain%Faces_ZX
  Cell=>Domain%Cell
  !........................
  nr_soildef=Domain%nr_soildef
  soil_type=>Domain%soil_type
  s_ixa=>Domain%s_ixa
  s_ixe=>Domain%s_ixe
  s_iya=>Domain%s_iya
  s_iye=>Domain%s_iye
  s_iza=>Domain%s_iza
  s_ize=>Domain%s_ize
  !........................
  NrW_FacesXY=Domain%NrW_FacesXY
  NrW_FacesYZ=Domain%NrW_FacesYZ
  NrW_FacesZX=Domain%NrW_FacesZX
  NrMP_FacesXY=Domain%NrMP_FacesXY
  NrMP_FacesYZ=Domain%NrMP_FacesYZ
  NrMP_FacesZX=Domain%NrMP_FacesZX
  NrMP_F=Domain%NrMP_F ! Allg. für Rand->NrMP_F* Splittung später 
  NrR_FacesXY=Domain%NrR_FacesXY
  NrR_FacesYZ=Domain%NrR_FacesYZ
  NrR_FacesZX=Domain%NrR_FacesZX
  NrRN_FacesXY=Domain%NrRN_FacesXY
  NrRN_FacesYZ=Domain%NrRN_FacesYZ
  NrRN_FacesZX=Domain%NrRN_FacesZX
  NrRW_FacesXY=Domain%NrRW_FacesXY
  NrRW_FacesYZ=Domain%NrRW_FacesYZ
  NrRW_FacesZX=Domain%NrRW_FacesZX
  !sub_face_gr=Domain%sub_face_gr
  !...............................
  NrW_Cells=Domain%NrW_Cells
  NrB_Cells=Domain%NrB_Cells
  NrMP_Cells=Domain%NrMP_Cells
  NrR_Cells=Domain%NrR_Cells
  NrRN_Cells=Domain%NrRN_Cells
  NrRW_Cells=Domain%NrRW_Cells
  !............................
  FU=>Domain%WeiFU
  FV=>Domain%WeiFV
  FW=>Domain%WeiFW
  VolC=>Domain%VolC
  !................................
  Emiss=>Domain%Emiss
  Emiss0=>Domain%Emiss0
  Nemiss=>Domain%Nemiss
  Mask=>Domain%Mask
  Mask0=>Domain%Mask0
  !................................
  AnzahlNachbar=Domain%AnzahlNachbar
  Nachbars=>Domain%Nachbars
  RefLevel=Domain%RefLevel
  Refine=Domain%Refine
  RefineX=Domain%RefineX
  RefineY=Domain%RefineY
  RefineZ=Domain%RefineZ
  Boundary=Domain%Boundary
  TypeW=Domain%TypeW
  TypeE=Domain%TypeE
  TypeS=Domain%TypeS
  TypeN=Domain%TypeN
  TypeT=Domain%TypeT
  TypeB=Domain%TypeB
END SUBROUTINE SetDomain

SUBROUTINE SetNachbar(Nachbar)

  TYPE (Nachbar_T) :: Nachbar

  ibn=Nachbar%ib
  jx0=Nachbar%ix0
  jx1=Nachbar%ix1
  jy0=Nachbar%iy0
  jy1=Nachbar%iy1
  jz0=Nachbar%iz0
  jz1=Nachbar%iz1
  RefineNachbar=Nachbar%Refine
  RefineNachbarX=Nachbar%RefineX
  RefineNachbarY=Nachbar%RefineY
  RefineNachbarZ=Nachbar%RefineZ
  IncrX=Nachbar%IncrX
  IncrY=Nachbar%IncrY
  IncrZ=Nachbar%IncrZ
  CopyCase=Nachbar%CopyCase
  nType=Nachbar%nType

END SUBROUTINE SetNachbar

!SUBROUTINE Set_EnvBlkFace(Domain,Face)
!  TYPE (Domain_T) :: Domain 
!  TYPE FaceP_T :: Face

SUBROUTINE Set_EnvBlkFace(Domain,ib)

  TYPE (Domain_T) :: Domain 
  INTEGER :: ib

  ix0=Domain%ix0
  ix1=Domain%ix1
  iy0=Domain%iy0
  iy1=Domain%iy1
  iz0=Domain%iz0
  iz1=Domain%iz1
  dx=>Domain%dx
  dy=>Domain%dy
  dz=>Domain%dz
  xP=>Domain%xP
  yP=>Domain%yP
  zP=>Domain%zP
  !.......................
  Vertices=>Domain%Vertices
  Faces_XY=>Domain%Faces_XY
  Faces_YZ=>Domain%Faces_YZ
  Faces_ZX=>Domain%Faces_ZX
  !........................
  NrW_FacesXY=Domain%NrW_FacesXY
  NrW_FacesYZ=Domain%NrW_FacesYZ
  NrW_FacesZX=Domain%NrW_FacesZX
  NrMP_FacesXY=Domain%NrMP_FacesXY
  NrMP_FacesYZ=Domain%NrMP_FacesYZ
  NrMP_FacesZX=Domain%NrMP_FacesZX
  NrMP_F=Domain%NrMP_F ! Allg. für Rand->NrMP_F* Splittung später 
  NrR_FacesXY=Domain%NrR_FacesXY
  NrR_FacesYZ=Domain%NrR_FacesYZ
  NrR_FacesZX=Domain%NrR_FacesZX
  NrRN_FacesXY=Domain%NrRN_FacesXY
  NrRN_FacesYZ=Domain%NrRN_FacesYZ
  NrRN_FacesZX=Domain%NrRN_FacesZX
  !...............................
END SUBROUTINE Set_EnvBlkFace


SUBROUTINE DomainAllocate(Domain)
   TYPE (Domain_T)           ::  &
      Domain
   ! Local variables
   INTEGER                   ::  &
      ix,iy,iz
   INTEGER                   ::  & 
      stat_P, stat_d, stat_V,    & ! error status variables Points,d[x,y,z],Verts
      stat_E, stat_F, stat_C,    & ! error status variables Edges,Faces,Cells
      stat_W, stat_VolC            ! error status variables Wei*,VolC
   CHARACTER (LEN=65)        ::  &
      nerrmsg                      ! error message
   CHARACTER (LEN=15)        ::  &
      nroutine                     ! name of this subroutine

   CALL Set(Domain)

!  common definition area is [ix0:ix1,iy0:iy1,iz0:iz1]
!  mit Rand
   ALLOCATE(Domain%xP(ix0-1:ix1+1),STAT=stat_P)
   ALLOCATE(Domain%yP(iy0-1:iy1+1),STAT=stat_P)
   ALLOCATE(Domain%zP(iz0-1:iz1+1),STAT=stat_P)
!  mit Rand
   ALLOCATE(Domain%dx(ix0:ix1+1),STAT=stat_d)
   ALLOCATE(Domain%dy(iy0:iy1+1),STAT=stat_d)
   ALLOCATE(Domain%dz(iz0:iz1+1),STAT=stat_d)

!  mit Rand
   ALLOCATE(Domain%Vertices(ix0-1:ix1+1,iy0-1:iy1+1,iz0-1:iz1+1),STAT=stat_V)

!  ohne Rand
   ALLOCATE(Domain%Edges_X(ix0+1:ix1,iy0:iy1,iz0:iz1),STAT=stat_E)
   ALLOCATE(Domain%Edges_Y(ix0:ix1,iy0+1:iy1,iz0:iz1),STAT=stat_E)
   ALLOCATE(Domain%Edges_Z(ix0:ix1,iy0:iy1,iz0+1:iz1),STAT=stat_E)

!  mit Rand
   ALLOCATE(Domain%Faces_XY(ix0:ix1+1,iy0:iy1+1,iz0-1:iz1+1),STAT=stat_F)
   ALLOCATE(Domain%Faces_ZX(ix0:ix1+1,iy0-1:iy1+1,iz0:iz1+1),STAT=stat_F)
   ALLOCATE(Domain%Faces_YZ(ix0-1:ix1+1,iy0:iy1+1,iz0:iz1+1),STAT=stat_F)

!  mit Rand
   ALLOCATE(Domain%Cell(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1),STAT=stat_C)
   !!Cells-Block
   !ALLOCATE(Domain%Cell(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
   !!RandCellen-Block, when not allocate corner cells
   !ALLOCATE(Domain%Cell(ix0,iy0+1:iy1,iz0+1:iz1))   !w
   !ALLOCATE(Domain%Cell(ix1+1,iy0+1:iy1,iz0+1:iz1))  !e
   !ALLOCATE(Domain%Cell(ix0+1:ix1,iy0,iz0+1:iz1))   !s
   !ALLOCATE(Domain%Cell(ix0+1:ix1,iy1+1,iz0+1:iz1))  !n
   !ALLOCATE(Domain%Cell(ix0+1:ix1,iy0+1:iy1,iz0))   !b
   !ALLOCATE(Domain%Cell(ix0+1:ix1,iy0+1:iy1,iz1+1))  !t

!  mit Rand
   ALLOCATE(Domain%WeiFU(ix0-1:ix1+1,iy0+1:iy1,iz0+1:iz1),STAT=stat_W) ! [x]-+
   ALLOCATE(Domain%WeiFV(ix0+1:ix1,iy0-1:iy1+1,iz0+1:iz1),STAT=stat_W) ! [y]-+
   ALLOCATE(Domain%WeiFW(ix0+1:ix1,iy0+1:iy1,iz0-1:iz1+1),STAT=stat_W) ! [z]-+

!  mit Rand
   ALLOCATE(Domain%VolC(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1),STAT=stat_VolC)
   !!VolC-Block
   !ALLOCATE(Domain%VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
   !!RandCellsVolC, when not allocate corner cells
   !ALLOCATE(Domain%VolC(ix0,iy0+1:iy1,iz0+1:iz1))  !w
   !ALLOCATE(Domain%VolC(ix1+1,iy0+1:iy1,iz0+1:iz1))  !e
   !ALLOCATE(Domain%VolC(ix0+1:ix1,iy0,iz0+1:iz1))  !s
   !ALLOCATE(Domain%VolC(ix0+1:ix1,iy1+1,iz0+1:iz1))  !n
   !ALLOCATE(Domain%VolC(ix0+1:ix1,iy0+1:iy1,iz0))  !b
   !ALLOCATE(Domain%VolC(ix0+1:ix1,iy0+1:iy1,iz1+1))  !t
   IF(stat_P/=0 .OR. stat_d/=0 .OR. stat_V/=0) THEN
     nroutine= 'DomainAllocate'
     nerrmsg = 'Allocate-Error: "stat_P, stat_d, stat_V"  aus DomainAllocate'
     Write(*,*) nerrmsg
     !CALL tool_break (my_cart_id, 1001, nerrmsg, nroutine)
   END IF
   IF(stat_E/=0 .OR. stat_F/=0 .OR. stat_C/=0) THEN
     nroutine= 'DomainAllocate'
     nerrmsg = 'Allocate-Error:  "stat_E, stat_F, stat_C" aus DomainAllocate'
     Write(*,*) nerrmsg
     !CALL tool_break (my_cart_id, 1002, nerrmsg, nroutine)
   END IF
   IF(stat_W/=0 .OR. stat_VolC/=0 ) THEN                        
     nroutine= 'DomainAllocate'
     nerrmsg = 'Allocate-Error:  "stat_W, stat_VolC"  aus DomainAllocate'
     Write(*,*) nerrmsg
     !CALL tool_break (my_cart_id, 1003, nerrmsg, nroutine)
   END IF
!  ohne Rand
   ALLOCATE(Domain%Emiss (ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1:SpeciesNum))
   ALLOCATE(Domain%Emiss0(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1:SpeciesNum))
   ALLOCATE(Domain%Nemiss(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
   ALLOCATE(Domain%Mask (ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
   ALLOCATE(Domain%Mask0(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))

END SUBROUTINE DomainAllocate


SUBROUTINE DEALLOC_VertP(Edge)
  TYPE (Edge_T) Edge

  IF (ASSOCIATED(Edge%VertS)) THEN
    DEALLOCATE(Edge%VertS) 
  END IF
END SUBROUTINE DEALLOC_VertP

SUBROUTINE DEALLOC_EdgeP(Face)
  TYPE (Face_T) :: Face

   IF (ASSOCIATED(Face%Edge1)) THEN
     CALL DEALLOC_VertP(Face%Edge1) 
     DEALLOCATE(Face%Edge1)
   END IF
   IF (ASSOCIATED(Face%Edge2)) THEN
     CALL DEALLOC_VertP(Face%Edge2) 
     DEALLOCATE(Face%Edge2)
   END IF
   IF (ASSOCIATED(Face%Edge3)) THEN
     CALL DEALLOC_VertP(Face%Edge3) 
     DEALLOCATE(Face%Edge3)
   END IF
   IF (ASSOCIATED(Face%Edge4)) THEN
     CALL DEALLOC_VertP(Face%Edge4)
     DEALLOCATE(Face%Edge4)
   END IF
END SUBROUTINE DEALLOC_EdgeP

SUBROUTINE DEALLOC_Cell_FaceP(Cell)
   TYPE (Cell_T) ::Cell

   IF (ASSOCIATED(Cell%Face1)) THEN
     CALL DEALLOC_EdgeP(Cell%Face1)
     DEALLOCATE(Cell%Face1)
   END IF
   IF (ASSOCIATED(Cell%Face2)) THEN
     CALL DEALLOC_EdgeP(Cell%Face2)
     DEALLOCATE(Cell%Face2)
   END IF
   IF (ASSOCIATED(Cell%Face3)) THEN
     CALL DEALLOC_EdgeP(Cell%Face3)
     DEALLOCATE(Cell%Face3)
   END IF
   IF (ASSOCIATED(Cell%Face4)) THEN
     CALL DEALLOC_EdgeP(Cell%Face4)
     DEALLOCATE(Cell%Face4)
   END IF
   IF (ASSOCIATED(Cell%Face5)) THEN
     CALL DEALLOC_EdgeP(Cell%Face5)
     DEALLOCATE(Cell%Face5)
   END IF
   IF (ASSOCIATED(Cell%Face6)) THEN
     CALL DEALLOC_EdgeP(Cell%Face6)
     DEALLOCATE(Cell%Face6)
   END IF
END SUBROUTINE DEALLOC_Cell_FaceP

SUBROUTINE DEALLOC_Cell_CutsListP(Cell)
   TYPE (Cell_T) ::Cell

   IF (ASSOCIATED(Cell%VertCut)) THEN
      DEALLOCATE(Cell%VertCut)
   END IF
   IF (ASSOCIATED(Cell%VertCutOnlyP)) THEN
      DEALLOCATE(Cell%VertCutOnlyP)
   END IF
END SUBROUTINE DEALLOC_Cell_CutsListP


SUBROUTINE DEALLOC_Cell_Parts(Cell)
   TYPE (Cell_T) ::Cell
   
   CALL DEALLOC_Cell_CutsListP(Cell)
   CALL DEALLOC_Cell_FaceP(Cell)
END SUBROUTINE DEALLOC_Cell_Parts


SUBROUTINE DomainDeAllocate(Domain)

   TYPE (Domain_T) :: Domain

   INTEGER :: i,j,k
   CALL Set(Domain)

   DEALLOCATE(Domain%xP)
   DEALLOCATE(Domain%yP)
   DEALLOCATE(Domain%zP)

   DEALLOCATE(Domain%dx)
   DEALLOCATE(Domain%dy)
   DEALLOCATE(Domain%dz)

   DEALLOCATE(Domain%Vertices)

   DO i=ix0+1,ix1
     DO j=iy0,iy1
       DO k=iz0,iz1
         IF (ASSOCIATED(Domain%Edges_X(i,j,k)%Edge)) THEN
           DEALLOCATE(Domain%Edges_X(i,j,k)%Edge)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Edges_X)
   DO i=ix0,ix1
     DO j=iy0+1,iy1
       DO k=iz0,iz1
         IF (ASSOCIATED(Domain%Edges_Y(i,j,k)%Edge)) THEN
           DEALLOCATE(Domain%Edges_Y(i,j,k)%Edge)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Edges_Y)
   DO i=ix0,ix1
     DO j=iy0,iy1
       DO k=iz0+1,iz1
         IF (ASSOCIATED(Domain%Edges_Z(i,j,k)%Edge)) THEN
           DEALLOCATE(Domain%Edges_Z(i,j,k)%Edge)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Edges_Z)

   DO i=ix0+1,ix1
     DO j=iy0+1,iy1
       DO k=iz0,iz1
         IF (ASSOCIATED(Domain%Faces_XY(i,j,k)%Face)) THEN
           DEALLOCATE(Domain%Faces_XY(i,j,k)%Face)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Faces_XY)
   DO i=ix0+1,ix1
     DO j=iy0,iy1
       DO k=iz0+1,iz1
         IF (ASSOCIATED(Domain%Faces_ZX(i,j,k)%Face)) THEN
           DEALLOCATE(Domain%Faces_ZX(i,j,k)%Face)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Faces_ZX)
   DO i=ix0,ix1
     DO j=iy0+1,iy1
       DO k=iz0+1,iz1
         IF (ASSOCIATED(Domain%Faces_YZ(i,j,k)%Face)) THEN
           DEALLOCATE(Domain%Faces_YZ(i,j,k)%Face)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Faces_YZ)  

   DO i=ix0,ix1+1
     DO j=iy0,iy1+1
       DO k=iz0,iz1+1
         IF (ASSOCIATED(Domain%Cell(i,j,k)%Cell)) THEN
          ! CALL DEALLOC_Cell_Parts(Domain%Cell(i,j,k)%Cell)
           !If (i /= ix0 .AND. i /= ix1+1 ) THEN
           !CALL DEALLOC_Cell_CutsListP(Domain%Cell(i,j,k)%Cell)
           !CALL DEALLOC_Cell_FaceP(Domain%Cell(i,j,k)%Cell)
           !END IF
           DEALLOCATE(Domain%Cell(i,j,k)%Cell)
         END IF
       END DO
     END DO
   END DO
   DEALLOCATE(Domain%Cell)

   DEALLOCATE(Domain%WeiFU)
   DEALLOCATE(Domain%WeiFV)
   DEALLOCATE(Domain%WeiFW)
   DEALLOCATE(Domain%VolC)

   DEALLOCATE(Domain%Emiss)
   DEALLOCATE(Domain%Emiss0)
   DEALLOCATE(Domain%Nemiss)
   DEALLOCATE(Domain%Mask)
   DEALLOCATE(Domain%Mask0)
  
   DEALLOCATE(Domain%Nachbars)

END SUBROUTINE DomainDeAllocate

END MODULE Domain_Mod

