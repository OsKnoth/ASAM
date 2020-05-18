MODULE Vector1FacePar_Mod
  USE Parallel_Mod
  USE Exchange_Mod
  USE Floor_Mod

  USE Vector1Face_Mod

  IMPLICIT NONE

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_Vec1Face
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_Vec1Face
  END INTERFACE
  INTERFACE Axpby
    MODULE PROCEDURE Axpby_Vec1Face
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_Vec1Face
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE Deallocate_Vec1Face
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_Vec1Face1 &
                    ,Allocate_Vec1Face2 
  END INTERFACE
CONTAINS 

SUBROUTINE Scale_Vec1Face(alpha,y)


  REAL(RealKind), INTENT(IN) :: alpha
  TYPE(Vector1Face_T), INTENT(INOUT) :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL ScaleV(alpha,y(ibLoc))
  END DO

END SUBROUTINE SCALE_Vec1Face

SUBROUTINE  Axpby_Vec1Face(alpha,x,beta,y)

  IMPLICIT NONE

  REAL(RealKind) :: alpha,beta
  TYPE(Vector1Face_T), INTENT(IN)  :: x(:)
  TYPE(Vector1Face_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Axpby(alpha,x(ibLoc),beta,y(ibLoc))
  END DO

END SUBROUTINE  Axpby_Vec1Face

SUBROUTINE  Copy_Vec1Face(x,y)

  IMPLICIT NONE

  TYPE(Vector1Face_T), INTENT(IN)  :: x(:)
  TYPE(Vector1Face_T), INTENT(INOUT)  :: y(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Copy(x(ibLoc),y(ibLoc))
  END DO

END SUBROUTINE  Copy_Vec1Face

SUBROUTINE Value_Vec1Face(x,Value)

  IMPLICIT NONE

  TYPE(Vector1Face_T), INTENT(OUT) :: x(:)
  REAL(RealKind), INTENT(IN) :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    x(ibLoc)=Value
  END DO

END SUBROUTINE Value_Vec1Face

SUBROUTINE  Allocate_Vec1Face1(x,VectorComponents)

  TYPE(Vector1Face_T), POINTER :: x(:)
  INTEGER :: VectorComponents

  ALLOCATE(x(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(x(ibLoc),VectorComponents)
  END DO

END SUBROUTINE  Allocate_Vec1Face1

SUBROUTINE  Allocate_Vec1Face2(x,VectorComponents1,VectorComponents2)

  TYPE(Vector1Face_T), POINTER :: x(:)
  INTEGER :: VectorComponents1,VectorComponents2

  ALLOCATE(x(nbLoc))
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Allocate(x(ibLoc),VectorComponents1,VectorComponents2)
  END DO

END SUBROUTINE  Allocate_Vec1Face2

SUBROUTINE  Deallocate_Vec1Face(x)

  IMPLICIT NONE

  TYPE(Vector1Face_T), POINTER :: x(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Deallocate(x(ibLoc))
  END DO
  DEALLOCATE(x)

END SUBROUTINE  Deallocate_Vec1Face

END MODULE Vector1FacePar_Mod

