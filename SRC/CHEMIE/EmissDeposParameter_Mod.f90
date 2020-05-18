MODULE EmissDeposParameter_Mod

  USE Chemie_Mod
  USE Distribution_Mod
  USE InitAerosol_Mod
  USE Physics_Mod

  IMPLICIT NONE

! New keyword for all emissions:
!   EmissionStartTime in Namelist ModelTransport

! Point emission:
! Keyword = Emiss in Namelist ModelTransport
! Data file = IniFile in Namelist ModelTransport
  TYPE EmiPoint_T 
    INTEGER :: ix,iy,iz
    REAL(RealKind), POINTER :: c(:)
    REAL(RealKind), POINTER :: Par(:)
  END TYPE EmiPoint_T 
  TYPE EmiPointBlock_T
    INTEGER :: NumEmiPoint=0
    TYPE(EmiPoint_T), POINTER :: EmiPoint(:)
  END TYPE EmiPointBlock_T
  INTEGER :: NumberOfPointSnap=0
  TYPE PointSnap_T
    CHARACTER*6 :: Type
    INTEGER :: NumGas=0
    INTEGER, POINTER :: Gas(:)
    INTEGER :: NumPar=0
    CHARACTER*10, POINTER :: TypePar(:)
    TYPE(EmiPointBlock_T), POINTER :: EmiPointBlock(:)
  END TYPE PointSnap_T
  TYPE(PointSnap_T), ALLOCATABLE :: PointSnap(:)

! Area emission:
! Keyword = Emiss in Namelist ModelTransport
! Data file = IniFile in Namelist ModelTransport
  TYPE EmiDomain_T    
    REAL(RealKind) :: xQ0,xQ1,yQ0,yQ1,zQ0,zQ1
    REAL(RealKind) :: n1,n2
    INTEGER :: NumGas
    CHARACTER(20), POINTER :: NameG(:)
    REAL(RealKind), POINTER :: EmiG(:)
    INTEGER, POINTER :: PosG(:)
    INTEGER :: NumAero
    CHARACTER(20), POINTER :: NameA(:)
    REAL(RealKind), POINTER :: EmiA(:,:)
    INTEGER :: NumberOfModes
    TYPE(AeroDistribution_T), POINTER :: AeroDistribution(:)
    TYPE(ImpactorStages_T) :: ImpactorStages
    INTEGER, POINTER :: PosA(:)
    TYPE(BlockCellSource_T), POINTER :: BlockCell(:)
  END TYPE EmiDomain_T
  TYPE BlockCellSource_T 
    INTEGER :: NumOfEmiCells
    TYPE(CellSource_T), POINTER :: Cell(:) 
  END TYPE BlockCellSource_T 
  
  TYPE CellSource_T 
    INTEGER :: ix,iy,iz
    REAL(RealKind) :: Frac
    REAL(RealKind) :: FL
  END TYPE CellSource_T 
  INTEGER :: NumberOfEmiDomainLoc
  TYPE (EmiDomain_T), POINTER :: EmiDomain(:)
  TYPE GasEmi_T    
    INTEGER :: Pos
    REAL(RealKind), POINTER :: Par(:) 
  END TYPE GasEmi_T
  TYPE AerosolEmi_T    
    INTEGER :: Pos
    INTEGER :: NumberOfModes
    TYPE(AeroDistribution_T), POINTER :: AeroDistribution(:)
    TYPE(ImpactorStages_T) :: ImpactorStages
    REAL(RealKind), POINTER :: Emi(:,:) 
  END TYPE AerosolEmi_T
  INTEGER :: NumberOfEmiDomain=0
  INTEGER :: NumGasEmi=0
  TYPE(GasEmi_T), ALLOCATABLE :: GasEmi(:)
  INTEGER :: NumAerosolEmi=0
  TYPE(AerosolEmi_T), ALLOCATABLE :: AerosolEmi(:)

! Street emission (analogous to point emission): ! Hinneburg
! Keyword = EmissStreet in Namelist ModelTransport
! Input file = InputFileName.Emission (supplied by GRID model
!              analogously to Weight files)
  TYPE EmiStreet_T 
    INTEGER :: ix,iy,iz
    REAL(RealKind), POINTER :: c(:)
    REAL(RealKind), POINTER :: Par(:)
  END TYPE EmiStreet_T 
  TYPE EmiStreetBlock_T
    INTEGER :: NumEmiStreet=0
    TYPE(EmiStreet_T), POINTER :: EmiStreet(:)
  END TYPE EmiStreetBlock_T
  INTEGER :: NumberOfStreetSnap=0
  TYPE StreetSnap_T
    CHARACTER*6 :: Type
    INTEGER :: NumGas=0
    INTEGER, POINTER :: Gas(:)
    INTEGER :: NumPar=0
    CHARACTER*10, POINTER :: TypePar(:)
    TYPE(EmiStreetBlock_T), POINTER :: EmiStreetBlock(:)
  END TYPE StreetSnap_T
  TYPE(StreetSnap_T), ALLOCATABLE :: StreetSnap(:)

! Marine values
  TYPE SeaGasEmi_T
    INTEGER :: Pos
    REAL(RealKind):: Par
    REAL(RealKind):: Konz
    REAL(RealKind):: Schmidt
    REAL(RealKind):: Henry
    REAL(RealKind):: Depos
  END TYPE SeaGasEmi_T
  TYPE(SeaGasEmi_T), ALLOCATABLE :: SeaGasEmi(:)
  INTEGER :: NumSeaTrans=0

  REAL(RealKind), PRIVATE :: Brownexp(0:9)
  REAL(RealKind), PRIVATE :: Impact_alpha(0:9)
  REAL(RealKind), PRIVATE :: Intercept_A(0:9)

!                        (0)     (1)     (2)       (3)        (4)        (5)     (6)        (7)     (8)     (9)
!  Land use:            bare    urban   savannah  deciduous  coniferous mixed              annual  grass   sea
!                       soil    area              forest     forest     forest  shrubland  crops   land
  DATA Brownexp     /   5.4d-1, 5.6d-1, 5.4d-1,   5.8d-1,    5.6d-1,    5.6d-1, 5.4d-1,    5.4d-1, 5.4d-1, 5.d-1 /
  DATA Impact_alpha /   5.0d1,  1.5d0,  1.2d0,    6.d-1,     1.d0,      8.d-1,  1.3d0,     1.2d0,  1.2d0,  1.d2 /
  DATA Intercept_A  /  -9.0d9,  1.d-2,  5.d-3,    5.d-3,     2.d-3,     5.d-3,  1.d-2,     5.d-3,  5.d-3, -9.d9 /

! Fire values  
  TYPE FireAeroEmiss_T
    INTEGER :: Pos
    REAL(RealKind):: Konz
  END TYPE FireAeroEmiss_T
  INTEGER :: NumFireAeroEmiss=0
  TYPE(FireAeroEmiss_T), ALLOCATABLE :: FireAeroEmiss(:)
  

CONTAINS

FUNCTION SchmidtF(Diff)

  REAL(RealKind) :: SchmidtF
  REAL(RealKind) :: Diff

  SchmidtF=KinViscAir/Diff ! Schmidtzahl

END FUNCTION SchmidtF

! Sedimentationsgeschwindigkeit

FUNCTION VFinalF(RhoPart,Rho,Rad)

  REAL(RealKind) :: VFinalF
  REAL(RealKind) :: RhoPart,Rho,Rad

  VFinalF=(2.d0*Rad**2*(RhoPart-Rho)*SLIP(Rad))/(9.d0*VISC_air)

END FUNCTION VFinalF

FUNCTION BrownF(Schmidt,Landuse)

  INTEGER :: Landuse
  REAL(RealKind) :: Schmidt
  REAL(RealKind) :: BrownF

  BrownF=Schmidt**(-Brownexp(Landuse))

END FUNCTION BrownF

FUNCTION StokeF(VFinal,ustern,Landuse)

  INTEGER :: Landuse
  REAL(RealKind) :: VFinal,ustern
  REAL(RealKind) :: StokeF

  IF (Intercept_A(Landuse).GT.0.d0) THEN
    StokeF=VFinal*ustern/(Grav*Intercept_A(Landuse))
  ELSE
    StokeF=VFinal*ustern**2/(Grav*KinViscAir)
  END IF

END FUNCTION StokeF

FUNCTION ImpactF(Stoke,Landuse)

  INTEGER :: Landuse
  REAL(RealKind) :: ImpactF
  REAL(RealKind) :: Stoke

  ImpactF=(Stoke/(Impact_alpha(Landuse)+Stoke))**2

END FUNCTION ImpactF

FUNCTION InterceptionF(Diff,Landuse)

  INTEGER :: Landuse
  REAL(RealKind) :: InterceptionF
  REAL(RealKind) :: Diff

  IF (Intercept_A(Landuse).GT.0.d0) THEN
    InterceptionF=2.d0*(Diff/Intercept_A(Landuse))**2 ! A LandUseFactor
  ELSE
    InterceptionF=0.d0
  END IF

END FUNCTION InterceptionF

FUNCTION SchmidtGas(SST,a0,a1,a2,a3)

  REAL(RealKind) :: SchmidtGas
  REAL(RealKind) :: SST,a0,a1,a2,a3

  SchmidtGas=a0-a1*SST+a2*SST**2-a3*SST**3 ! Wannkinhof 1992

END FUNCTION SchmidtGas

FUNCTION TransCoeff(U10,Schmidt)

  REAL(RealKind) :: TransCoeff
  REAL(RealKind) :: U10,Schmidt

    TransCoeff=0.31d0*U10**2/(Schmidt/660d0)**(0.5) ! cm/h ! Wannkinhof92
    TransCoeff=TransCoeff/3.6d5 ! m/s

END FUNCTION TransCoeff


END MODULE EmissDeposParameter_Mod
