MODULE MatSpInd_Mod

  USE Kind_Mod

  IMPLICIT NONE

  TYPE SpInd
    INTEGER :: m,n
    INTEGER :: NumNonZero
    INTEGER, POINTER :: RowInd(:)
    INTEGER, POINTER :: ColInd(:)
    REAL(RealKind), POINTER :: Val(:)
  END TYPE SpInd

CONTAINS

SUBROUTINE SpAllocate_SpInd(A)

!
! Allocate A
!
  IMPLICIT NONE

  TYPE(SpInd) :: A

  IF (.NOT.ASSOCIATED(A%RowInd)) THEN
    ALLOCATE(A%RowInd(A%NumNonZero))
  END IF
  IF (.NOT.ASSOCIATED(A%ColInd)) THEN
    ALLOCATE(A%ColInd(A%NumNonZero))
  END IF
  IF (.NOT.ASSOCIATED(A%Val)) THEN
    ALLOCATE(A%Val(A%NumNonZero))
  END IF

END SUBROUTINE SpAllocate_SpInd

SUBROUTINE SpDeallocate_SpInd(A)

!
! Allocate A
!

  TYPE(SpInd) :: A

  A%m=0
  A%n=0
  A%NumNonZero=0
  IF (ASSOCIATED(A%RowInd)) THEN
    DEALLOCATE(A%RowInd)
  END IF
  IF (ASSOCIATED(A%ColInd)) THEN
    DEALLOCATE(A%ColInd)
  END IF
  IF (ASSOCIATED(A%Val)) THEN
    DEALLOCATE(A%Val)
  END IF

END SUBROUTINE SpDeallocate_SpInd

SUBROUTINE SpNullify_SpInd(A)

!
! Nullify Pointers of A 
!

  TYPE(SpInd) :: A

  NULLIFY(A%RowInd)
  NULLIFY(A%ColInd)
  NULLIFY(A%Val)

END SUBROUTINE SpNullify_SpInd

END MODULE MatSpInd_Mod
