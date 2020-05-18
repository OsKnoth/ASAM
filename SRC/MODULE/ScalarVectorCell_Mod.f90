MODULE ScalarVectorCellPar_Mod

  USE Floor_Mod
  USE Vector4Cell_Mod
  USE ScalarCell_Mod

  IMPLICIT NONE

  INTERFACE Assign
    MODULE PROCEDURE Vector2ScalarCell,Scalar2ScalarCell
  END INTERFACE
  INTERFACE Mult
    MODULE PROCEDURE ScalarVectorMult
  END INTERFACE
CONTAINS 

SUBROUTINE Scalar2ScalarCell(ScalarCell1,ScalarCell2)

  TYPE(ScalarCell_T) :: ScalarCell1(:)
  TYPE(ScalarCell_T), INTENT(IN)  :: ScalarCell2(:)

  DO ibLoc=1,nbLoc
    ScalarCell1(ibLoc)%c=>ScalarCell2(ibLoc)%c
  END DO

END SUBROUTINE Scalar2ScalarCell

SUBROUTINE Vector2ScalarCell(ScalarCell,VectorCell,Pos)

  TYPE(ScalarCell_T) :: ScalarCell(:)
  TYPE(Vector4Cell_T), INTENT(IN)  :: VectorCell(:)
  INTEGER :: Pos

  DO ibLoc=1,nbLoc
    ScalarCell(ibLoc)%c=>VectorCell(ibLoc)%Vec(Pos)%c
  END DO

END SUBROUTINE Vector2ScalarCell

SUBROUTINE ScalarVectorMult(ScalarCell,VectorCell,Pos)

  TYPE(ScalarCell_T), INTENT(IN) :: ScalarCell(:)
  TYPE(Vector4Cell_T), INTENT(INOUT)  :: VectorCell(:)
  INTEGER :: Pos

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    VectorCell(ibLoc)%Vec(Pos)%c=ScalarCell(ibLoc)%c*VectorCell(ibLoc)%Vec(Pos)%c
  END DO

END SUBROUTINE ScalarVectorMult


END MODULE ScalarVectorCellPar_Mod

