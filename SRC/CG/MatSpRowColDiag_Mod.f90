MODULE MatSpRowColDiag_Mod

  USE Kind_Mod

  IMPLICIT NONE

  TYPE SpRowColDiag
    INTEGER :: m,n
    INTEGER, POINTER :: RowPtr(:)=>NULL()
    INTEGER, POINTER :: DiagPtr(:)=>NULL()
    INTEGER, POINTER :: ColInd(:)=>NULL()
    INTEGER, POINTER :: Permu(:)=>NULL()
    INTEGER, POINTER :: InvPer(:)=>NULL()
    REAL(RealKind), POINTER :: Val(:)=>NULL()
  END TYPE SpRowColDiag

CONTAINS

SUBROUTINE SpNullify_SpRowColDiag(A)

!
! Nullify Pointers of A 
!

  TYPE(SpRowColDiag) :: A

  NULLIFY(A%RowPtr)
  NULLIFY(A%DiagPtr)
  NULLIFY(A%ColInd)
  NULLIFY(A%Permu)
  NULLIFY(A%InvPer)
  NULLIFY(A%Val)

END SUBROUTINE SpNullify_SpRowColDiag

SUBROUTINE SpDeallocate_SpRowColDiag(A)

!
! Deallocate A
!

  TYPE(SpRowColDiag) :: A

  A%m=0
  A%n=0
  IF (ASSOCIATED(A%RowPtr)) THEN
    DEALLOCATE(A%RowPtr)
  END IF
  IF (ASSOCIATED(A%DiagPtr)) THEN
    DEALLOCATE(A%DiagPtr)
  END IF
  IF (ASSOCIATED(A%ColInd)) THEN
    DEALLOCATE(A%ColInd)
  END IF
  IF (ASSOCIATED(A%Permu)) THEN
    DEALLOCATE(A%Permu)
  END IF
  IF (ASSOCIATED(A%InvPer)) THEN
    DEALLOCATE(A%InvPer)
  END IF
  IF (ASSOCIATED(A%Val)) THEN
    DEALLOCATE(A%Val)
  END IF

END SUBROUTINE SpDeallocate_SpRowColDiag


END MODULE MatSpRowColDiag_Mod
