MODULE Kind_Mod

  INTEGER, PARAMETER :: RealKind=8
  INTEGER, PARAMETER :: Real4Kind=4
  INTEGER, PARAMETER :: IntKind=4

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
