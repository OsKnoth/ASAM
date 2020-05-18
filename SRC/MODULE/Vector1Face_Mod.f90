MODULE Vector1Face_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Sp_Mod
  USE Parameter_Mod
  USE Control_Mod
  USE Parallel_Mod

  IMPLICIT NONE


  TYPE VecF1_T
    REAL(RealKind), POINTER :: c(:)=>NULL()
  END TYPE VecF1_T

  TYPE Vector1Face_T
    TYPE (VecF1_T), POINTER :: Vec(:)
  END TYPE Vector1Face_T
  TYPE Vec1Vec1tor1Face_T
    TYPE (Vector1Face_T), POINTER :: Vec(:)
  END TYPE Vec1Vec1tor1Face_T

  TYPE (Vector1Face_T), POINTER :: Vector1FaceAct(:)



  INTERFACE Add
    MODULE PROCEDURE Add_Vector,Add_Vector_Vector
  END INTERFACE Add
  INTERFACE AddScalar
    MODULE PROCEDURE AddScalar_MatVector
  END INTERFACE


  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_Vector1Face &
                    ,Copy_Vec1Vec1 &
                    ,Copy_Vector1Scalar &
                    ,Copy_Vector_Vector1Scalar &
                    ,Copy_Vector_Vector1Vector1 

  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_Vector1Face
  END INTERFACE
  INTERFACE Axpby
    MODULE PROCEDURE Axpby_Vector1Face
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_Vector1Face
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_Vector1Face &
                    ,Scale_Vector
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE DeallocateVector1Face
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_Vector,AllocateVector1Face1,AllocateVector1Face2
  END INTERFACE
  INTERFACE Assign
    MODULE PROCEDURE Assign_Vector1Vector1
  END INTERFACE
CONTAINS

SUBROUTINE Assign_Vector1Vector1(Vec1,Vec2)

  TYPE(VecF1_T), INTENT(OUT) :: Vec1
  TYPE(VecF1_T), INTENT(IN) :: Vec2


  Vec1%c=>Vec2%c

END SUBROUTINE Assign_Vector1Vector1

SUBROUTINE Add_Vector(alpha,Vec1,Vec2)

  REAL(RealKind) :: alpha
  TYPE(VecF1_T), INTENT(IN) :: Vec1
  TYPE(VecF1_T), INTENT(INOUT) :: Vec2

  Vec2%c=alpha*Vec1%c+Vec2%c

END SUBROUTINE Add_Vector

SUBROUTINE Add_Vector_Vector(alpha,Vec1,Vec2)

  REAL(RealKind) :: alpha
  TYPE(VecF1_T), INTENT(IN) :: Vec1(:)
  TYPE(VecF1_T), INTENT(INOUT) :: Vec2(:)

  INTEGER :: i

  DO i=1,SIZE(Vec1(:))
    Vec2(i)%c=alpha*Vec1(i)%c+Vec2(i)%c
  END DO

END SUBROUTINE Add_Vector_Vector

SUBROUTINE Neg_MatVector(Vec)

  TYPE(VecF1_T), INTENT(INOUT) :: Vec

  Vec%c=-Vec%c

END SUBROUTINE Neg_MatVector

SUBROUTINE AddScalar_MatVector(Vec1,Vec2,Value)

  TYPE(VecF1_T) :: Vec1,Vec2
  REAL(RealKind)  :: Value

! Vec1=Vec2+Value

  Vec1%c=Vec2%c+Value

END SUBROUTINE AddScalar_MatVector

SUBROUTINE Scale_Vector(Value,Vec)

  REAL(RealKind), INTENT(IN) :: Value
  TYPE(VecF1_T), INTENT(INOUT) :: Vec

  Vec%c=Value*Vec%c
END SUBROUTINE Scale_Vector

SUBROUTINE Deallocate_Vector(Vec)

  TYPE(VecF1_T), INTENT(INOUT) :: Vec

  IF (ASSOCIATED(Vec%c)) THEN
    DEALLOCATE(Vec%c)
  END IF

END SUBROUTINE Deallocate_Vector

SUBROUTINE Copy_Vector1Scalar(Vec,Scalar)

  TYPE(VecF1_T), INTENT(INOUT) :: Vec
  REAL(RealKind), INTENT(IN) :: Scalar

  Vec%c=Scalar

END SUBROUTINE Copy_Vector1Scalar

SUBROUTINE Copy_Vector_Vector1Scalar(Vec,Scalar)

  TYPE(VecF1_T), INTENT(INOUT) :: Vec(:)
  REAL(RealKind), INTENT(IN) :: Scalar

  INTEGER :: i
  
  DO i=1,SIZE(Vec)
    Vec(i)%c=Scalar
  END DO 

END SUBROUTINE Copy_Vector_Vector1Scalar

SUBROUTINE Copy_Vector_Vector1Vector1(Vec1,Vec2)

  TYPE(VecF1_T), INTENT(INOUT) :: Vec1(:)
  TYPE(VecF1_T), INTENT(IN) :: Vec2(:)

  INTEGER :: i
 
  DO i=1,SIZE(Vec1)
    Vec1(i)%c=Vec2(i)%c
  END DO

END SUBROUTINE Copy_Vector_Vector1Vector1

SUBROUTINE Copy_Vec1Vec1(Vec1,Vec2)

  TYPE(VecF1_T), INTENT(INOUT) :: Vec1
  TYPE(VecF1_T), INTENT(IN) :: Vec2

  Vec1%c=Vec2%c
END SUBROUTINE Copy_Vec1Vec1

SUBROUTINE AllocateVector1Face1(Vector1Face,VectorComponents)

  TYPE (Vector1Face_T) :: Vector1Face
  INTEGER :: VectorComponents

  INTEGER :: ic
  ALLOCATE(Vector1Face%Vec(VectorComponents))
  DO ic=1,VectorComponents
    ALLOCATE(Vector1Face%Vec(ic)%c(iz0:iz1))
  END DO

END SUBROUTINE AllocateVector1Face1

SUBROUTINE AllocateVector1Face2(Vector1Face,VectorComponents1,VectorComponents2)

  TYPE (Vector1Face_T) :: Vector1Face
  INTEGER :: VectorComponents1,VectorComponents2

  INTEGER :: ic
  ALLOCATE(Vector1Face%Vec(VectorComponents1:VectorComponents2))
  DO ic=LBOUND(Vector1Face%Vec,1),UBOUND(Vector1Face%Vec,1)
    ALLOCATE(Vector1Face%Vec(ic)%c(iz0:iz1))
    Vector1Face%Vec(ic)%c=0.0e0
  END DO

END SUBROUTINE AllocateVector1Face2

SUBROUTINE AllocateVector1FaceVector1Face(Vector1Face,Vector1FaceO,LBound0,UBound0)

  TYPE (Vector1Face_T) :: Vector1Face
  TYPE (Vector1Face_T) :: Vector1FaceO
  INTEGER, OPTIONAL :: LBound0,UBound0

  INTEGER :: LBound0Loc,UBound0Loc
  INTEGER :: ic

  IF (PRESENT(LBound0)) THEN
    LBound0Loc=LBound0
  ELSE
    LBound0Loc=LBOUND(Vector1FaceO%Vec,1)
  END IF
  IF (PRESENT(UBound0)) THEN
    UBound0Loc=UBound0
  ELSE
    UBound0Loc=UBOUND(Vector1FaceO%Vec,1)
  END IF
  ALLOCATE(Vector1Face%Vec(LBound0Loc:UBound0Loc))
  DO ic=LBound0Loc,UBound0Loc
    CALL Allocate(Vector1Face%Vec(ic),Vector1FaceO%Vec(ic))
    Vector1Face%Vec(ic)%c=0.0e0
  END DO

END SUBROUTINE AllocateVector1FaceVector1Face

SUBROUTINE Allocate_Vector(Vec1,Vec2)

  TYPE(VecF1_T), INTENT(INOUT) :: Vec1
  TYPE(VecF1_T), INTENT(IN) :: Vec2

  ALLOCATE(Vec1%c(LBOUND(Vec2%c,1):UBOUND(Vec2%c,1)))
  Vec1%c=Zero

END SUBROUTINE Allocate_Vector

SUBROUTINE DeallocateVector1Face(Vector1Face)

  TYPE (Vector1Face_T) :: Vector1Face

  INTEGER :: ic

  DO ic=LBOUND(Vector1Face%Vec,1),UBOUND(Vector1Face%Vec,1)
    DEALLOCATE(Vector1Face%Vec(ic)%c)
  END DO
  DEALLOCATE(Vector1Face%Vec)

END SUBROUTINE DeallocateVector1Face

SUBROUTINE Axpy_Vector1Face(alpha,Vector1Face1,Vector1Face2)

  REAL(RealKind) :: alpha
  TYPE (Vector1Face_T) :: Vector1Face1,Vector1Face2

  INTEGER :: ic
  INTEGER :: n

  IF (alpha/=Zero) THEN
    DO ic=MAX(LBOUND(Vector1Face1%Vec,1),LBOUND(Vector1Face2%Vec,1)), &
          MIN(UBOUND(Vector1Face1%Vec,1),UBOUND(Vector1Face2%Vec,1))
      Vector1Face2%Vec(ic)%c=alpha*Vector1Face1%Vec(ic)%c+Vector1Face2%Vec(ic)%c
    END DO
  END IF

END SUBROUTINE Axpy_Vector1Face

SUBROUTINE Axpby_Vector1Face(alpha,Vector1Face1,beta,Vector1Face2)

  REAL(RealKind) :: alpha,beta
  TYPE (Vector1Face_T) :: Vector1Face1,Vector1Face2

  INTEGER :: ic

  DO ic=MAX(LBOUND(Vector1Face1%Vec,1),LBOUND(Vector1Face2%Vec,1)), &
        MIN(UBOUND(Vector1Face1%Vec,1),UBOUND(Vector1Face2%Vec,1))
    Vector1Face2%Vec(ic)%c=alpha*Vector1Face1%Vec(ic)%c+beta*Vector1Face2%Vec(ic)%c
  END DO

END SUBROUTINE Axpby_Vector1Face

SUBROUTINE Copy_Vector1Face(Vector1Face1,Vector1Face2)

  TYPE (Vector1Face_T) :: Vector1Face1,Vector1Face2

  INTEGER :: ic

  DO ic=MAX(LBOUND(Vector1Face1%Vec,1),LBOUND(Vector1Face2%Vec,1)), &
        MIN(UBOUND(Vector1Face1%Vec,1),UBOUND(Vector1Face2%Vec,1))
    Vector1Face2%Vec(ic)%c=Vector1Face1%Vec(ic)%c
  END DO

END SUBROUTINE Copy_Vector1Face

SUBROUTINE Scale_Vector1Face(alpha,Vector1Face1)

  REAL(RealKind) :: alpha
  TYPE (Vector1Face_T) :: Vector1Face1

  INTEGER :: ic

  DO ic=LBOUND(Vector1Face1%Vec,1),UBOUND(Vector1Face1%Vec,1)
    Vector1Face1%Vec(ic)%c=alpha*Vector1Face1%Vec(ic)%c
  END DO

END SUBROUTINE Scale_Vector1Face

SUBROUTINE Value_Vector1Face(Vector1Face1,Value)

  TYPE(Vector1Face_T), INTENT(OUT)  :: Vector1Face1
  REAL(RealKind), INTENT(IN)  :: Value

  INTEGER :: ic

  DO ic=LBOUND(Vector1Face1%Vec,1),UBOUND(Vector1Face1%Vec,1)
    Vector1Face1%Vec(ic)%c=Value
  END DO

END SUBROUTINE Value_Vector1Face

END MODULE Vector1Face_Mod

