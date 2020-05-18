MODULE Vector34Cell_Mod

  USE Vector3Cell_Mod
  USE Vector4Cell_Mod

  IMPLICIT NONE

  INTERFACE Assign
    MODULE PROCEDURE Assign_Vector34Cell
  END INTERFACE
  TYPE (Vector3Cell_T), POINTER :: Vec34(:)

CONTAINS

SUBROUTINE Assign_Vector34Cell(Vector3Cell,Vector4Cell)

  TYPE (Vector3Cell_T) :: Vector3Cell
  TYPE (Vector4Cell_T) :: Vector4Cell
  INTEGER :: Vector3Components
  INTEGER :: Vector4Components
  INTEGER :: ic,ic3,it

  IF (.NOT.ASSOCIATED(Vector3Cell%Vec)) THEN
    Vector3Components=0
    DO ic=1,SIZE(Vector4Cell%Vec)
      Vector3Components=Vector3Components+SIZE(Vector4Cell%Vec(ic)%c,4) 
    END DO
    ALLOCATE(Vector3Cell%Vec(Vector3Components))
  END IF
  ic3=0
  DO ic=1,SIZE(Vector4Cell%Vec)
    DO it=LBOUND(Vector4Cell%Vec(ic)%c,4),UBOUND(Vector4Cell%Vec(ic)%c,4)
      ic3=ic3+1
      Vector3Cell%Vec(ic3)%c=>Vector4Cell%Vec(ic)%c(:,:,:,it)
    END DO
  END DO

END SUBROUTINE Assign_Vector34Cell
END MODULE Vector34Cell_Mod
