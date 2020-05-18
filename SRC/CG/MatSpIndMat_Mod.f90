MODULE MatSpIndMat_Mod

  USE Kind_Mod
  USE Matrix_Mod

  IMPLICIT NONE

  TYPE SpIndMat
    INTEGER :: m,n
    INTEGER :: NumNonZero
    INTEGER, POINTER :: RowInd(:)
    INTEGER, POINTER :: ColInd(:)
    TYPE(A22_T), POINTER :: Val(:)
  END TYPE SpIndMat

CONTAINS

SUBROUTINE SpAllocate_SpIndMat(A)

!
! Allocate A
!
  IMPLICIT NONE

  TYPE(SpIndMat) :: A

  IF (.NOT.ASSOCIATED(A%RowInd)) THEN
    ALLOCATE(A%RowInd(A%NumNonZero))
  END IF
  IF (.NOT.ASSOCIATED(A%ColInd)) THEN
    ALLOCATE(A%ColInd(A%NumNonZero))
  END IF
  IF (.NOT.ASSOCIATED(A%Val)) THEN
    ALLOCATE(A%Val(A%NumNonZero))
  END IF

END SUBROUTINE SpAllocate_SpIndMat

SUBROUTINE SpDeallocate_SpIndMat(A)

!
! Allocate A
!

  TYPE(SpIndMat) :: A

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

END SUBROUTINE SpDeallocate_SpIndMat

SUBROUTINE SpNullify_SpIndMat(A)

!
! Nullify Pointers of A 
!

  TYPE(SpIndMat) :: A

  NULLIFY(A%RowInd)
  NULLIFY(A%ColInd)
  NULLIFY(A%Val)

END SUBROUTINE SpNullify_SpIndMat

END MODULE MatSpIndMat_Mod
