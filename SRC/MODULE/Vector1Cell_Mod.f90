MODULE Vector1Cell_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Sp_Mod
  USE Parameter_Mod
  USE Control_Mod
  USE Parallel_Mod

  IMPLICIT NONE


  TYPE Vec1_T
    REAL(RealKind), POINTER :: c(:)=>NULL()
  END TYPE Vec1_T

  TYPE Vector1Cell_T
    TYPE (Vec1_T), POINTER :: Vec(:)
  END TYPE Vector1Cell_T
  TYPE Vec1Vec1tor1Cell_T
    TYPE (Vector1Cell_T), POINTER :: Vec(:)
  END TYPE Vec1Vec1tor1Cell_T

  TYPE (Vector1Cell_T), POINTER :: Vector1CellAct(:)



  INTERFACE Add
    MODULE PROCEDURE Add_Vector,Add_Vector_Vector
  END INTERFACE Add
  INTERFACE AddScalar
    MODULE PROCEDURE AddScalar_MatVector
  END INTERFACE


  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_Vector1Cell &
                    ,Copy_Vec1Vec1 &
                    ,Copy_Vector1Scalar &
                    ,Copy_Vector_Vector1Scalar &
                    ,Copy_Vector_Vector1Vector1 

  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_Vector1Cell
  END INTERFACE
  INTERFACE Axpby
    MODULE PROCEDURE Axpby_Vector1Cell
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_Vector1Cell
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_Vector1Cell &
                    ,Scale_Vector
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE DeallocateVector1Cell
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_Vector,AllocateVector1Cell1,AllocateVector1Cell2
  END INTERFACE
  INTERFACE Assign
    MODULE PROCEDURE Assign_Vector1Vector1
  END INTERFACE
CONTAINS

SUBROUTINE Assign_Vector1Vector1(Vec1,Vec2)

  TYPE(Vec1_T), INTENT(OUT) :: Vec1
  TYPE(Vec1_T), INTENT(IN) :: Vec2


  Vec1%c=>Vec2%c

END SUBROUTINE Assign_Vector1Vector1

SUBROUTINE Add_Vector(alpha,Vec1,Vec2)

  REAL(RealKind) :: alpha
  TYPE(Vec1_T), INTENT(IN) :: Vec1
  TYPE(Vec1_T), INTENT(INOUT) :: Vec2

  Vec2%c=alpha*Vec1%c+Vec2%c

END SUBROUTINE Add_Vector

SUBROUTINE Add_Vector_Vector(alpha,Vec1,Vec2)

  REAL(RealKind) :: alpha
  TYPE(Vec1_T), INTENT(IN) :: Vec1(:)
  TYPE(Vec1_T), INTENT(INOUT) :: Vec2(:)

  INTEGER :: i

  DO i=1,SIZE(Vec1(:))
    Vec2(i)%c=alpha*Vec1(i)%c+Vec2(i)%c
  END DO

END SUBROUTINE Add_Vector_Vector

SUBROUTINE Neg_MatVector(Vec)

  TYPE(Vec1_T), INTENT(INOUT) :: Vec

  Vec%c=-Vec%c

END SUBROUTINE Neg_MatVector

SUBROUTINE AddScalar_MatVector(Vec1,Vec2,Value)

  TYPE(Vec1_T) :: Vec1,Vec2
  REAL(RealKind)  :: Value

! Vec1=Vec2+Value

  Vec1%c=Vec2%c+Value

END SUBROUTINE AddScalar_MatVector

SUBROUTINE Scale_Vector(Value,Vec)

  REAL(RealKind), INTENT(IN) :: Value
  TYPE(Vec1_T), INTENT(INOUT) :: Vec

  Vec%c=Value*Vec%c
END SUBROUTINE Scale_Vector

SUBROUTINE Deallocate_Vector(Vec)

  TYPE(Vec1_T), INTENT(INOUT) :: Vec

  IF (ASSOCIATED(Vec%c)) THEN
    DEALLOCATE(Vec%c)
  END IF

END SUBROUTINE Deallocate_Vector

SUBROUTINE Copy_Vector1Scalar(Vec,Scalar)

  TYPE(Vec1_T), INTENT(INOUT) :: Vec
  REAL(RealKind), INTENT(IN) :: Scalar

  Vec%c=Scalar

END SUBROUTINE Copy_Vector1Scalar

SUBROUTINE Copy_Vector_Vector1Scalar(Vec,Scalar)

  TYPE(Vec1_T), INTENT(INOUT) :: Vec(:)
  REAL(RealKind), INTENT(IN) :: Scalar

  INTEGER :: i
  
  DO i=1,SIZE(Vec)
    Vec(i)%c=Scalar
  END DO 

END SUBROUTINE Copy_Vector_Vector1Scalar

SUBROUTINE Copy_Vector_Vector1Vector1(Vec1,Vec2)

  TYPE(Vec1_T), INTENT(INOUT) :: Vec1(:)
  TYPE(Vec1_T), INTENT(IN) :: Vec2(:)

  INTEGER :: i
 
  DO i=1,SIZE(Vec1)
    Vec1(i)%c=Vec2(i)%c
  END DO

END SUBROUTINE Copy_Vector_Vector1Vector1

SUBROUTINE Copy_Vec1Vec1(Vec1,Vec2)

  TYPE(Vec1_T), INTENT(INOUT) :: Vec1
  TYPE(Vec1_T), INTENT(IN) :: Vec2

  Vec1%c=Vec2%c
END SUBROUTINE Copy_Vec1Vec1

SUBROUTINE AllocateVector1Cell1(Vector1Cell,VectorComponents)

  TYPE (Vector1Cell_T) :: Vector1Cell
  INTEGER :: VectorComponents

  INTEGER :: ic
  ALLOCATE(Vector1Cell%Vec(VectorComponents))
  DO ic=1,VectorComponents
    ALLOCATE(Vector1Cell%Vec(ic)%c(iz0:iz1+1))
  END DO

END SUBROUTINE AllocateVector1Cell1

SUBROUTINE AllocateVector1Cell2(Vector1Cell,VectorComponents1,VectorComponents2)

  TYPE (Vector1Cell_T) :: Vector1Cell
  INTEGER :: VectorComponents1,VectorComponents2

  INTEGER :: ic
  ALLOCATE(Vector1Cell%Vec(VectorComponents1:VectorComponents2))
  DO ic=LBOUND(Vector1Cell%Vec,1),UBOUND(Vector1Cell%Vec,1)
    ALLOCATE(Vector1Cell%Vec(ic)%c(iz0:iz1+1))
    Vector1Cell%Vec(ic)%c=0.0e0
  END DO

END SUBROUTINE AllocateVector1Cell2

SUBROUTINE AllocateVector1CellVector1Cell(Vector1Cell,Vector1CellO,LBound0,UBound0)

  TYPE (Vector1Cell_T) :: Vector1Cell
  TYPE (Vector1Cell_T) :: Vector1CellO
  INTEGER, OPTIONAL :: LBound0,UBound0

  INTEGER :: LBound0Loc,UBound0Loc
  INTEGER :: ic

  IF (PRESENT(LBound0)) THEN
    LBound0Loc=LBound0
  ELSE
    LBound0Loc=LBOUND(Vector1CellO%Vec,1)
  END IF
  IF (PRESENT(UBound0)) THEN
    UBound0Loc=UBound0
  ELSE
    UBound0Loc=UBOUND(Vector1CellO%Vec,1)
  END IF
  ALLOCATE(Vector1Cell%Vec(LBound0Loc:UBound0Loc))
  DO ic=LBound0Loc,UBound0Loc
    CALL Allocate(Vector1Cell%Vec(ic),Vector1CellO%Vec(ic))
    Vector1Cell%Vec(ic)%c=0.0e0
  END DO

END SUBROUTINE AllocateVector1CellVector1Cell

SUBROUTINE Allocate_Vector(Vec1,Vec2)

  TYPE(Vec1_T), INTENT(INOUT) :: Vec1
  TYPE(Vec1_T), INTENT(IN) :: Vec2

  ALLOCATE(Vec1%c(LBOUND(Vec2%c,1):UBOUND(Vec2%c,1)))
  Vec1%c=Zero

END SUBROUTINE Allocate_Vector

SUBROUTINE DeallocateVector1Cell(Vector1Cell)

  TYPE (Vector1Cell_T) :: Vector1Cell

  INTEGER :: ic

  DO ic=LBOUND(Vector1Cell%Vec,1),UBOUND(Vector1Cell%Vec,1)
    DEALLOCATE(Vector1Cell%Vec(ic)%c)
  END DO
  DEALLOCATE(Vector1Cell%Vec)

END SUBROUTINE DeallocateVector1Cell

SUBROUTINE Axpy_Vector1Cell(alpha,Vector1Cell1,Vector1Cell2)

  REAL(RealKind) :: alpha
  TYPE (Vector1Cell_T) :: Vector1Cell1,Vector1Cell2

  INTEGER :: ic
  INTEGER :: n

  IF (alpha/=Zero) THEN
    DO ic=MAX(LBOUND(Vector1Cell1%Vec,1),LBOUND(Vector1Cell2%Vec,1)), &
          MIN(UBOUND(Vector1Cell1%Vec,1),UBOUND(Vector1Cell2%Vec,1))
      Vector1Cell2%Vec(ic)%c=alpha*Vector1Cell1%Vec(ic)%c+Vector1Cell2%Vec(ic)%c
    END DO
  END IF

END SUBROUTINE Axpy_Vector1Cell

SUBROUTINE Axpby_Vector1Cell(alpha,Vector1Cell1,beta,Vector1Cell2)

  REAL(RealKind) :: alpha,beta
  TYPE (Vector1Cell_T) :: Vector1Cell1,Vector1Cell2

  INTEGER :: ic

  DO ic=MAX(LBOUND(Vector1Cell1%Vec,1),LBOUND(Vector1Cell2%Vec,1)), &
        MIN(UBOUND(Vector1Cell1%Vec,1),UBOUND(Vector1Cell2%Vec,1))
    Vector1Cell2%Vec(ic)%c=alpha*Vector1Cell1%Vec(ic)%c+beta*Vector1Cell2%Vec(ic)%c
  END DO

END SUBROUTINE Axpby_Vector1Cell

SUBROUTINE Copy_Vector1Cell(Vector1Cell1,Vector1Cell2)

  TYPE (Vector1Cell_T) :: Vector1Cell1,Vector1Cell2

  INTEGER :: ic

  DO ic=MAX(LBOUND(Vector1Cell1%Vec,1),LBOUND(Vector1Cell2%Vec,1)), &
        MIN(UBOUND(Vector1Cell1%Vec,1),UBOUND(Vector1Cell2%Vec,1))
    Vector1Cell2%Vec(ic)%c=Vector1Cell1%Vec(ic)%c
  END DO

END SUBROUTINE Copy_Vector1Cell

SUBROUTINE Scale_Vector1Cell(alpha,Vector1Cell1)

  REAL(RealKind) :: alpha
  TYPE (Vector1Cell_T) :: Vector1Cell1

  INTEGER :: ic

  DO ic=LBOUND(Vector1Cell1%Vec,1),UBOUND(Vector1Cell1%Vec,1)
    Vector1Cell1%Vec(ic)%c=alpha*Vector1Cell1%Vec(ic)%c
  END DO

END SUBROUTINE Scale_Vector1Cell

SUBROUTINE Value_Vector1Cell(Vector1Cell1,Value)

  TYPE(Vector1Cell_T), INTENT(OUT)  :: Vector1Cell1
  REAL(RealKind), INTENT(IN)  :: Value

  INTEGER :: ic

  DO ic=LBOUND(Vector1Cell1%Vec,1),UBOUND(Vector1Cell1%Vec,1)
    Vector1Cell1%Vec(ic)%c=Value
  END DO

END SUBROUTINE Value_Vector1Cell

END MODULE Vector1Cell_Mod

