MODULE Kind_Mod

  IMPLICIT NONE

! INTEGER, PARAMETER :: RealKind=8
  INTEGER, PARAMETER :: RealKind=4
  INTEGER, PARAMETER :: Real4Kind=4
  INTEGER, PARAMETER :: IntKind=4


! ---Numerical Parameter-----------------------------------------------------------------------------------------
  REAL(RealKind), PARAMETER :: Eps        = 1.e-30_RealKind
  REAL(RealKind), PARAMETER :: Zero       = 0.0e0_RealKind
  REAL(RealKind), PARAMETER :: Half       = 0.5e0_RealKind
  REAL(RealKind), PARAMETER :: One        = 1.0e0_RealKind
  REAL(RealKind), PARAMETER :: NegOne     =-1.0e0_RealKind
  REAL(RealKind), PARAMETER :: Two        = 2.0e0_RealKind
  REAL(RealKind), PARAMETER :: Three      = 3.0e0_RealKind
  REAL(RealKind), PARAMETER :: Four       = 4.0e0_RealKind

CONTAINS


  SUBROUTINE native_4byte_real( realIn, realOut )
  IMPLICIT NONE
  REAL(4), INTENT(IN)                              :: realIn
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element
  REAL(4), INTENT(OUT)                             :: realOut
                                                   ! a single 32 bit, 4 byte
                                                   ! REAL data element, with
                                                   ! reverse byte order to
                                                   ! that of realIn
  realOut=realIn
  END SUBROUTINE native_4byte_real

END MODULE Kind_Mod
