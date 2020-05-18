MODULE Aerosol_Mod

  USE Kind_Mod
  USE Transport_Mod

! REAL(RealKind) :: Numin=1.d-10
! REAL(RealKind) :: bMax=10.0d0
! REAL(RealKind) :: MassEps=1.d-4
! REAL(RealKind) :: rStart=1.d-8
! REAL(RealKind) :: pFac=2.0d0
! INTEGER :: nFrac=40

  INTEGER :: iNC,iWater,iRelax,iClm,iNap
  REAL(RealKind), POINTER :: m(:)  ! Grid points of mass
  REAL(RealKind), POINTER :: r(:)  ! cell centered radius
!  REAL(RealKind) :: RhoAerosol=1.5d3 ! in kg/m^3
  REAL(RealKind) :: RhoAerosol=2.2d3 ! in kg/m^3

END MODULE Aerosol_Mod
