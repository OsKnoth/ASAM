MODULE Domain_Mod

  USE Kind_Mod

  IMPLICIT NONE

  TYPE Nachbar_T
    CHARACTER(2) :: nTYPE
    INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1       ! Coordinates of border in block coordinates
    INTEGER :: iNx0,iNx1,iNy0,iNy1,iNz0,iNz1 ! Coordinates of border in neighbour coordinates
    INTEGER :: ixO,ixI,iyO,iyI,izO,izI
    INTEGER :: iNxO,iNxI,iNyO,iNyI,iNzO,iNzI
    INTEGER :: ib
    INTEGER :: ibLoc
    INTEGER :: Refine
    INTEGER :: RefineX
    INTEGER :: RefineY
    INTEGER :: RefineZ
    INTEGER :: IncrX
    INTEGER :: IncrY
    INTEGER :: IncrZ
    INTEGER :: CopyCase
    REAL(RealKind) :: dLoc
  END TYPE Nachbar_T
  TYPE gliedT
    INTEGER :: OwnBlock         ! block receiving data
    INTEGER :: OwnBlockLoc         ! block receiving data
    INTEGER :: NeiBlock(2)      ! block sending data
    INTEGER :: ix0U,ix1U,iy0U,iy1U,iz0U,iz1U
    CHARACTER(2) ::   side
    INTEGER :: RefOwn
    INTEGER :: RefNei
    INTEGER :: RefSide
    INTEGER :: IncrX
    INTEGER :: IncrY
    INTEGER :: IncrZ
    TYPE(Nachbar_T), POINTER :: Nachbar
  END TYPE gliedT
  
  TYPE SoilCell_T
    REAL(RealKind) :: T1             ! First soil layer temperature
    REAL(RealKind) :: T2             ! Second soil layer temperature
    REAL(RealKind) :: D              ! Exchange coefficient
    REAL(RealKind) :: k,Cp,RhoBeton !Waermeleitfähigkeit,WärmeKap und Dichte der Wand
  END TYPE SoilCell_T

  TYPE CanopyCell_T
    INTEGER                 :: NrCanopyLayers
    REAL(RealKind), POINTER :: LAD(:) ! Leaf Area Density
    REAL(RealKind)          :: LAI    ! Leaf Area Index, the summation of local LAD
    REAL(RealKind)          :: RadPtrRate=1.0d0 ! Penetration rate of radiation after going through the canopy
    REAL(RealKind)          :: GapFunc=1.0d0 ! gap function, how much of the sky will be seen under the canopy
    REAL(RealKind)          :: AvgTem=285.0d0 ! average temperature of the leaves inside the canopy, used to calculate the downward long wave radiation of canopy, it varies with time
    REAL(RealKind)          :: emissivity=0.990d0 ! emissivity of foreast, it should depends on the type of the forest
    REAL(RealKind)          :: Alb=0.0d0 ! albedo of the canopy
    REAL(RealKind)          :: AlbSoil=0.1d0 ! albedo of the soil under the canopy
    REAL(RealKind), ALLOCATABLE :: SunlitLeafT(:), ShadedLeafT(:) ! Leaf temperature
    REAL(RealKind), ALLOCATABLE :: Hc(:), LEc(:) ! Sensible heat flux and latent heat flux from canopy to air
    REAL(RealKind)          :: Rbs=0.0d0 ! direct radiation under the canopy
    REAL(RealKind)          :: Rds=0.0d0 ! diffuse radiation under the canopy
    REAL(RealKind)          :: RLs=0.0d0 ! longwave radiation under the canopy
  END TYPE CanopyCell_T

  TYPE BoundCell_T
    INTEGER :: ix,iy,iz
    INTEGER :: shad=1
    REAL(RealKind) :: xS,yS,zS
    REAL(RealKind) :: n1=0.0d0,n2=0.0d0,n3=0.0d0
    REAL(RealKind) :: n1G=0.0d0,n2G=0.0d0,n3G=0.0d0
    REAL(RealKind) :: xFLS,yFLS,zFLS
    REAL(RealKind) :: skyviewfactor=1.0d0    
    REAL(RealKind) :: FL=0.0d0
    REAL(RealKind) :: dL=0.0d0
    REAL(RealKind) :: z=0.0d0
    REAL(RealKind) :: zRauh=0.0d0
    REAL(RealKind) :: zRauhT=0.0d0
    REAL(RealKind) :: alb=0.0d0
    REAL(RealKind) :: ee=0.0d0
    REAL(RealKind) :: DragM=0.0d0  ! Drag coefficient for momentum
    REAL(RealKind) :: DragH=0.0d0  ! Drag coefficient for heat 
    REAL(RealKind) :: DragQ=0.0d0  ! Drag coefficient for Moisture 
    REAL(RealKind) :: FacEmission=0.0d0
    REAL(RealKind) :: raddirekt=0.0d0
    REAL(RealKind) :: raddiffus=0.0d0
    REAL(RealKind) :: radinfred=0.0d0
    REAL(RealKind) :: FluxSens=0.0d0
    REAL(RealKind) :: FluxLat=0.0d0
    REAL(RealKind) :: TSoil=0.0d0
    REAL(RealKind) :: QvSoil=0.0d0
    INTEGER        :: LandClass=9
    INTEGER, POINTER :: SoilType(:)
!   TYPE(SoilCell_T), Pointer :: SoilCell 
    TYPE(CanopyCell_T) :: CanopyCell
    REAL(RealKind) :: TeS                ! Surface temperature
    REAL(RealKind) :: qv=0.0d0           ! Surface humidity
    REAL(RealKind) :: ThetaS             ! Surface potential temperature
    REAL(RealKind) :: U10                ! Wind in 10m height ! Barthel
    REAL(RealKind) :: VT1                ! Tangential Wind in first level over ground
    REAL(RealKind) :: VT2                ! Tangential Wind in second level over ground
    REAL(RealKind) :: T1             ! First soil layer temperature
    REAL(RealKind) :: D              ! Exchange coefficient
    LOGICAL :: Land=.FALSE.
    LOGICAL :: LandCheck=.FALSE.
  END TYPE BoundCell_T 

  TYPE BoundFace_T
    INTEGER :: ix,iy,iz
    REAL(RealKind) :: xS,yS,zS
  END TYPE BoundFace_T

! Masking(Baum): analogous to Point-Emission ! Hinneburg
  TYPE TreePoint_T 
    INTEGER :: ix,iy,iz
    REAL(RealKind) :: c
  END TYPE TreePoint_T 
  TYPE TreePointBlock_T
    INTEGER :: NumTreePoint=0
    TYPE(TreePoint_T), POINTER :: TreePoint(:)
  END TYPE TreePointBlock_T
  TYPE PointTree_T
    TYPE(TreePointBlock_T), POINTER :: TreePointBlock(:)
  END TYPE PointTree_T
  TYPE(PointTree_T) :: PointTree
  REAL(RealKind) :: zRauhBaum=0.1d0,distBaum=0.2d0 ! Hinneburg: gilt innerhalb Baumschicht (Testwerte)

  TYPE Domain_T
    CHARACTER(3) :: TypeE='',TypeW='',TypeS='',TypeN='',TypeB='',TypeT=''
    REAL(RealKind) :: Boundary
    INTEGER :: nx,ny,nz            ! number of cells in x-,y-,z-direction
    INTEGER :: nc                  ! number of cells of the block
    INTEGER :: WriteOffsetC        ! Cell Offset for parallel I/O
    INTEGER :: WriteOffsetCgmv     ! Cell Offset for parallel I/O for gmv
    INTEGER :: WriteOffsetN        ! Node Offset for parallel I/O
    INTEGER :: WriteOffsetCSoilgmv ! Cell Offset Soil for parallel I/O for gmv
    INTEGER :: WriteOffsetCCutgmv  ! Cell Offset Cut for parallel I/O for gmv
    INTEGER :: igx0,igy0,igx1,igy1,igz0,igz1
    INTEGER :: ix0,iy0,ix1,iy1,iz0,iz1
    INTEGER :: ib
    LOGICAL :: JacAdvection
    LOGICAL :: JacDiffusion
    INTEGER :: RefLevel
    INTEGER :: Refine
    INTEGER :: RefineX
    INTEGER :: RefineY
    INTEGER :: RefineZ
    INTEGER :: Xshift
    INTEGER :: Yshift
    INTEGER :: Zshift
    REAL(RealKind) :: x0,x1,y0,y1,z0,z1
    REAL(RealKind) :: LenScaleGeom
    REAL(RealKind), POINTER :: dx(:),dy(:),dz(:)
    REAL(RealKind), POINTER :: zSDepth(:)
    REAL(RealKind), POINTER :: MetrXY(:),MetrXZ(:)
    REAL(RealKind), POINTER :: MetrYX(:),MetrYZ(:)
    REAL(RealKind), POINTER :: MetrZX(:),MetrZY(:)
    REAL(RealKind), POINTER :: xP(:),yP(:),zP(:)
    REAL(RealKind), POINTER :: xPG(:,:,:),yPG(:,:,:),zPG(:,:,:)
    REAL(RealKind), POINTER :: zH(:,:)
    REAL(RealKind), POINTER :: fCor(:,:)
    REAL(RealKind), POINTER :: WeiFU(:,:,:),WeiFV(:,:,:),WeiFW(:,:,:)
    REAL(RealKind), POINTER :: WeiFUT(:,:,:)=>NULL(),WeiFVT(:,:,:)=>NULL(),WeiFWT(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: WeiFUG(:,:,:)=>NULL(),WeiFVG(:,:,:)=>NULL(),WeiFWG(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: VolC(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: VolCD(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: VolCP(:,:,:)=>NULL()
    LOGICAL, POINTER :: Mask(:,:,:)=>NULL()
    INTEGER :: NumBoundCell
    TYPE (BoundCell_T), POINTER :: BoundCell(:)
    INTEGER :: NumBoundFaceU
    TYPE (BoundCell_T), POINTER :: BoundFaceU(:)
    INTEGER :: NumBoundFaceV
    TYPE (BoundCell_T), POINTER :: BoundFaceV(:)
    INTEGER :: NumBoundFaceW
    TYPE (BoundCell_T), POINTER :: BoundFaceW(:)
    INTEGER :: FreeCells
    INTEGER :: AnzahlNachbar
    TYPE (Nachbar_T), POINTER :: Nachbars(:)
    INTEGER, POINTER :: BoundCell3d(:,:,:)
    INTEGER :: nrsoillayers
  END TYPE Domain_T
  
  INTEGER :: ix0Out,ix1Out,iy0Out,iy1Out,iz0Out,iz1Out
  INTEGER :: nx,ny,nz 
  INTEGER :: WriteOffsetC
  INTEGER :: WriteOffsetCgmv
  INTEGER :: WriteOffsetN
  INTEGER :: WriteOffsetCSoilgmv
  INTEGER :: WriteOffsetCCutgmv
  INTEGER :: igx0,igy0,igx1,igy1,igz0,igz1
  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  INTEGER :: jx0,jx1,jy0,jy1,jz0,jz1
  INTEGER :: jNx0,jNx1,jNy0,jNy1,jNz0,jNz1
  INTEGER :: jxO,jxI,jyO,jyI,jzO,jzI
  INTEGER :: jNxO,jNxI,jNyO,jNyI,jNzO,jNzI
  REAL(RealKind) :: x0,x1,y0,y1,z0,z1
  REAL(RealKind) :: LenScaleGeom
  REAL(RealKind), POINTER :: dx(:),dy(:),dz(:)
  REAL(RealKind), POINTER :: zSDepth(:)
  REAL(RealKind), POINTER :: MetrXY(:),MetrXZ(:)
  REAL(RealKind), POINTER :: MetrYX(:),MetrYZ(:)
  REAL(RealKind), POINTER :: MetrZX(:),MetrZY(:)
  REAL(RealKind), POINTER :: xP(:),yP(:),zP(:)
  REAL(RealKind), POINTER :: xPG(:,:,:),yPG(:,:,:),zPG(:,:,:)
  REAL(RealKind), POINTER :: zH(:,:)
  REAL(RealKind), POINTER :: fCor(:,:)
  REAL(RealKind), POINTER :: FU(:,:,:),FV(:,:,:),FW(:,:,:)
  REAL(RealKind), POINTER :: FUT(:,:,:),FVT(:,:,:),FWT(:,:,:)
  REAL(RealKind), POINTER :: FUG(:,:,:),FVG(:,:,:),FWG(:,:,:)
  REAL(RealKind), POINTER :: WeiVol(:,:,:)
  REAL(RealKind), POINTER :: VolC(:,:,:)
  REAL(RealKind), POINTER :: VolCD(:,:,:)  
  REAL(RealKind), POINTER :: VolCP(:,:,:)
  REAL(RealKind), POINTER :: VolB(:,:,:)
  LOGICAL, POINTER :: Mask(:,:,:)
  INTEGER :: NumBoundCell
  TYPE (BoundCell_T), POINTER :: BoundCell(:)
  INTEGER :: NumBoundFaceU
  TYPE (BoundCell_T), POINTER :: BoundFaceU(:)
  INTEGER :: NumBoundFaceV
  TYPE (BoundCell_T), POINTER :: BoundFaceV(:)
  INTEGER :: NumBoundFaceW
  TYPE (BoundCell_T), POINTER :: BoundFaceW(:)
  INTEGER :: FreeCells
  INTEGER :: AnzahlNachbar
  TYPE (Nachbar_T), POINTER :: Nachbars(:)
  TYPE(Nachbar_T), POINTER :: Nachbar
  LOGICAL :: JacAdvection
  LOGICAL :: JacDiffusion
  INTEGER :: RefLevel
  INTEGER :: Refine
  INTEGER :: RefineX
  INTEGER :: RefineY
  INTEGER :: RefineZ
  REAL(RealKind) :: Boundary
  INTEGER :: nrsoillayers
  CHARACTER(3) :: TypeE,TypeW,TypeS,TypeN,TypeB,TypeT
  REAL(RealKind) :: dphi
  REAL(RealKind) :: MaxHeight

  INTEGER, POINTER :: BoundCell3d(:,:,:)
  
  INTEGER :: ibC,ibCLoc,ibN,ibNLoc
  INTEGER :: RefineNachbar
  INTEGER :: RefineNachbarX
  INTEGER :: RefineNachbarY
  INTEGER :: RefineNachbarZ
  INTEGER :: IncrX,IncrY,IncrZ
  INTEGER :: CopyCase
  REAL(RealKind) :: dLoc
  CHARACTER(2) :: nType

  INTERFACE Set
    MODULE PROCEDURE SetDomain,SetNachbar,SetChain
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE DomainAllocate
  END INTERFACE

CONTAINS

SUBROUTINE SetDomain(Domain)

  TYPE (Domain_T) :: Domain

  nx=Domain%nx
  ny=Domain%ny
  nz=Domain%nz
  WriteOffsetC=Domain%WriteOffsetC
  WriteOffsetCgmv=Domain%WriteOffsetCgmv
  WriteOffsetN=Domain%WriteOffsetN
  WriteOffsetCSoilgmv=Domain%WriteOffsetCSoilgmv
  WriteOffsetCCutgmv=Domain%WriteOffsetCCutgmv
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
  x0=Domain%x0
  x1=Domain%x1
  y0=Domain%y0
  y1=Domain%y1
  z0=Domain%z0
  z1=Domain%z1
  dx=>Domain%dx
  dy=>Domain%dy
  dz=>Domain%dz
  MetrXY=>Domain%MetrXY
  MetrXZ=>Domain%MetrXZ
  MetrYX=>Domain%MetrYX
  MetrYZ=>Domain%MetrYZ
  MetrZX=>Domain%MetrZX
  MetrZY=>Domain%MetrZY
  zSDepth=>Domain%zSDepth
  xP=>Domain%xP
  xPG=>Domain%xPG
  yP=>Domain%yP
  yPG=>Domain%yPG
  zP=>Domain%zP
  zPG=>Domain%zPG
  zH=>Domain%zH
  LenScaleGeom=Domain%LenScaleGeom
  fCor=>Domain%fCor
  FU=>Domain%WeiFU
  FV=>Domain%WeiFV
  FW=>Domain%WeiFW
  FUT=>Domain%WeiFUT
  FVT=>Domain%WeiFVT
  FWT=>Domain%WeiFWT
  FUG=>Domain%WeiFUG
  FVG=>Domain%WeiFVG
  FWG=>Domain%WeiFWG
  VolC=>Domain%VolC
  Mask=>Domain%Mask
  VolCD=>Domain%VolCD
  VolCP=>Domain%VolCP
  IF (ASSOCIATED(Domain%VolC)) THEN
    VolB=>Domain%VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)
  END IF
  NumBoundCell=Domain%NumBoundCell
  BoundCell=>Domain%BoundCell
  NumBoundFaceU=Domain%NumBoundFaceU
  BoundFaceU=>Domain%BoundFaceU
  NumBoundFaceV=Domain%NumBoundFaceV
  BoundFaceV=>Domain%BoundFaceV
  NumBoundFaceW=Domain%NumBoundFaceW
  BoundFaceW=>Domain%BoundFaceW
  FreeCells=Domain%FreeCells
  AnzahlNachbar=Domain%AnzahlNachbar
  Nachbars=>Domain%Nachbars
  JacAdvection=Domain%JacAdvection
  JacDiffusion=Domain%JacDiffusion
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
  BoundCell3d=>Domain%BoundCell3d
  nrsoillayers=Domain%nrsoillayers
END SUBROUTINE SetDomain

SUBROUTINE SetNachbar(Nachbar)

  TYPE (Nachbar_T) :: Nachbar

  ibn=Nachbar%ib
  ibnLoc=Nachbar%ibLoc
  jx0=Nachbar%ix0
  jx1=Nachbar%ix1
  jy0=Nachbar%iy0
  jy1=Nachbar%iy1
  jz0=Nachbar%iz0
  jz1=Nachbar%iz1
  jNx0=Nachbar%iNx0
  jNx1=Nachbar%iNx1
  jNy0=Nachbar%iNy0
  jNy1=Nachbar%iNy1
  jNz0=Nachbar%iNz0
  jNz1=Nachbar%iNz1
  jxO=Nachbar%ixO
  jxI=Nachbar%ixI
  jyO=Nachbar%iyO
  jyI=Nachbar%iyI
  jzO=Nachbar%izO
  jzI=Nachbar%izI
  jNxO=Nachbar%iNxO
  jNxI=Nachbar%iNxI
  jNyO=Nachbar%iNyO
  jNyI=Nachbar%iNyI
  jNzO=Nachbar%iNzO
  jNzI=Nachbar%iNzI
  RefineNachbar=Nachbar%Refine
  RefineNachbarX=Nachbar%RefineX
  RefineNachbarY=Nachbar%RefineY
  RefineNachbarZ=Nachbar%RefineZ
  IncrX=Nachbar%IncrX
  IncrY=Nachbar%IncrY
  IncrZ=Nachbar%IncrZ
  CopyCase=Nachbar%CopyCase
  dLoc=Nachbar%dLoc
  nType=Nachbar%nType

END SUBROUTINE SetNachbar

SUBROUTINE SetChain(Chain)

  TYPE (gliedT) :: Chain

  ibC=Chain%OwnBlock
  ibCLoc=Chain%OwnBlockLoc
  jx0=Chain%ix0U
  jx1=Chain%ix1U
  jy0=Chain%iy0U
  jy1=Chain%iy1U
  jz0=Chain%iz0U
  jz1=Chain%iz1U
  IncrX=Chain%IncrX
  IncrY=Chain%IncrY
  IncrZ=Chain%IncrZ
  nType=Chain%side
  Refine=Chain%RefOwn
  RefineNachbar=Chain%RefNei
  Nachbar=>Chain%Nachbar

END SUBROUTINE SetChain

SUBROUTINE DomainAllocate(Domain)

   TYPE (Domain_T) :: Domain
   CALL Set(Domain)

   ALLOCATE(Domain%xP(ix0:ix1))
   ALLOCATE(Domain%yP(iy0:iy1))
   ALLOCATE(Domain%zP(iz0:iz1))

   ALLOCATE(Domain%xPG(ix0:ix1,iy0:iy1,iz0:iz1))
   ALLOCATE(Domain%yPG(ix0:ix1,iy0:iy1,iz0:iz1))
   ALLOCATE(Domain%zPG(ix0:ix1,iy0:iy1,iz0:iz1))

   ALLOCATE(Domain%dx(ix0+1:ix1))
   ALLOCATE(Domain%dy(iy0+1:iy1))
   ALLOCATE(Domain%dz(iz0+1:iz1))

   ALLOCATE(Domain%MetrXY(iy0+1:iy1))
   ALLOCATE(Domain%MetrXZ(iz0+1:iz1))
   ALLOCATE(Domain%MetrYX(ix0+1:ix1))
   ALLOCATE(Domain%MetrYZ(iz0+1:iz1))
   ALLOCATE(Domain%MetrZX(ix0+1:ix1))
   ALLOCATE(Domain%MetrZY(iy0+1:iy1))
   ALLOCATE(Domain%zH(ix0:ix1+1,iy0:iy1+1))
   ALLOCATE(Domain%fCor(ix0+1:ix1,iy0+1:iy1))

   ! mit Rand
   ALLOCATE(Domain%WeiFU(ix0-1:ix1+1,iy0+1:iy1,iz0+1:iz1))
   ALLOCATE(Domain%WeiFV(ix0+1:ix1,iy0-1:iy1+1,iz0+1:iz1))
   ALLOCATE(Domain%WeiFW(ix0+1:ix1,iy0+1:iy1,iz0-1:iz1+1))

   ALLOCATE(Domain%WeiFUT(ix0:ix1,iy0+1:iy1,iz0+1:iz1))
   ALLOCATE(Domain%WeiFVT(ix0+1:ix1,iy0:iy1,iz0+1:iz1))
   ALLOCATE(Domain%WeiFWT(ix0+1:ix1,iy0+1:iy1,iz0:iz1))

   ALLOCATE(Domain%WeiFUG(ix0:ix1,iy0+1:iy1,iz0+1:iz1))
   ALLOCATE(Domain%WeiFVG(ix0+1:ix1,iy0:iy1,iz0+1:iz1))
   ALLOCATE(Domain%WeiFWG(ix0+1:ix1,iy0+1:iy1,iz0:iz1))

   ! mit Rand
   ALLOCATE(Domain%VolC(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1))
   ALLOCATE(Domain%Mask(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1))
   Domain%Mask=.TRUE.
   ALLOCATE(Domain%VolCD(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1))   
   ALLOCATE(Domain%VolCP(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))

END SUBROUTINE DomainAllocate

END MODULE Domain_Mod



