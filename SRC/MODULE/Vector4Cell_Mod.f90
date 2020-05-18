MODULE Vector4Cell_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Floor_Mod
  USE Sp_Mod
  USE Parameter_Mod
  USE Control_Mod
  USE Parallel_Mod

  IMPLICIT NONE

  TYPE Vec34_T
    REAL(RealKind), POINTER :: c(:,:,:,:)=>NULL()
  END TYPE Vec34_T

  TYPE Vec4_T
    REAL(RealKind), POINTER :: c(:,:,:,:)=>NULL()
    REAL(RealKind), POINTER :: cInt(:,:,:,:)=>NULL()
    REAL(RealKind), POINTER :: cB(:,:)=>NULL()
    TYPE(Vec34_T) :: Vec3
    INTEGER :: Rank(2)
  END TYPE Vec4_T

  TYPE Vector4Cell_T
    TYPE (Vec4_T), POINTER :: Vec(:)
  END TYPE Vector4Cell_T
  TYPE VecVector4Cell_T
    TYPE (Vector4Cell_T), POINTER :: Vec(:)
  END TYPE VecVector4Cell_T

  TYPE Matrix4Cell_T
    TYPE (Vec4_T), POINTER :: Mat(:,:)
  END TYPE Matrix4Cell_T

  TYPE JacMatrix4_T
    TYPE (Matrix4Cell_T), POINTER :: JacSLU=>NULL()
    TYPE (SpDiag), POINTER :: JacT=>NULL()
  END TYPE JacMatrix4_T


  TYPE SpMatrix4Cell_T
    INTEGER :: n
    INTEGER :: nMaxVec
    TYPE(Vec4_T), POINTER :: Mat(:)=>NULL()
    TYPE(SpRowColDiag), POINTER :: Struct=>NULL()
    LOGICAL :: Factor=.FALSE.
  END TYPE SpMatrix4Cell_T

  TYPE JacSpMatrix4_T
    TYPE (SpMatrix4Cell_T), POINTER :: JacSLU=>NULL()
    TYPE (SpDiag), POINTER :: JacTMom=>NULL()
    TYPE (SpDiag), POINTER :: JacTPot=>NULL()
    TYPE (SpDiag), POINTER :: JacFall=>NULL()
    TYPE (SpDiag), POINTER :: JacFallRhoL=>NULL()
  END TYPE JacSpMatrix4_T

  TYPE (Vector4Cell_T), POINTER :: VectorCellAct(:)
  TYPE (Matrix4Cell_T), POINTER :: MatrixCellAct(:)
  TYPE (SpMatrix4Cell_T), POINTER :: SpMatrixCellAct(:)

  INTEGER, PRIVATE :: ic,jc,kc,it
  INTEGER, PUBLIC :: ixM,iyM,izM,icM,itM
  INTEGER, PRIVATE :: ixc0,ixc1,iyc0,iyc1,izc0,izc1,itc0,itc1
  INTEGER, PUBLIC :: VectorComponents
  INTEGER, PUBLIC :: ixS,iyS,izS
  REAL(RealKind), ALLOCATABLE :: ScaleMat2(:)

  INTERFACE Mult
    MODULE PROCEDURE Mult_VecVec,Mult_Vec4,Mult_Vec5,Mult_Vec6
  END INTERFACE Mult
  INTERFACE MultD
    MODULE PROCEDURE MultD_Vec6
  END INTERFACE MultD
  INTERFACE Div
    MODULE PROCEDURE Div_VecVec
  END INTERFACE Div
  INTERFACE SetZero
    MODULE PROCEDURE SetZero_VecVec
  END INTERFACE SetZero
  INTERFACE Add
    MODULE PROCEDURE Add_Vector,Add_Vector_Vector           &
                    ,Add3_Vector_Vector,Add_Vector_Partiell
  END INTERFACE Add
  INTERFACE AddScalar
    MODULE PROCEDURE AddScalar_MatVector,AddScalarPos_MatVector
  END INTERFACE


  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_VectorCell &
                    ,Value_MatrixCell&
                    ,Value_SpMatrixCell&
                    ,Copy_ScalarVector &
                    ,Copy_VecVec &
                    ,Copy_VectorScalar &
                    ,Copy_Vector_VectorScalar &
                    ,Copy_Vector_VectorVector &
                    ,Copy_LowHigh 

  END INTERFACE
  INTERFACE Dot
    MODULE PROCEDURE Dot_VectorCell &
                    ,Dot_Vec4
  END INTERFACE
  INTERFACE Sum
    MODULE PROCEDURE Sum_VectorCell
  END INTERFACE
  INTERFACE Error
    MODULE PROCEDURE Error_VectorCell
  END INTERFACE
  INTERFACE ErrorMax
    MODULE PROCEDURE ErrorMax_VectorCell
  END INTERFACE
  INTERFACE Xpy
    MODULE PROCEDURE Xpy_VectorCell
  END INTERFACE
  INTERFACE Axpy
    MODULE PROCEDURE Axpy_VectorCell
  END INTERFACE
  INTERFACE Ax1mx2py
    MODULE PROCEDURE Ax1mx2py_VectorCell
  END INTERFACE
  INTERFACE Ax1px2py
    MODULE PROCEDURE Ax1px2py_VectorCell
  END INTERFACE
  INTERFACE Axpby
    MODULE PROCEDURE Axpby_VectorCell
  END INTERFACE
  INTERFACE Xpay
    MODULE PROCEDURE Xpay_VectorCell
  END INTERFACE
  INTERFACE Copy
    MODULE PROCEDURE Copy_VectorCell
  END INTERFACE
  INTERFACE Random
    MODULE PROCEDURE Random_VectorCell
  END INTERFACE
  INTERFACE ScaleV
    MODULE PROCEDURE Scale_VectorCell &
                    ,Scale_Vector
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE DeallocateVectorCell, &
                     DeallocateMatrixCell , &
                     Deallocate_Vector
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE AllocateVectorCell1, &
                     AllocateVectorCell2, &
                     AllocateMatrixCell, &
                     AllocateVectorCellVectorCell, &
                     AllocateMatrixCellVectorCell, &
                     SpDiagAllocate, &
                     Allocate_Vector,&
                     Allocate_Rank,&
                     Allocate_Dim,& 
                     Allocate_Dim1 
  END INTERFACE
  INTERFACE Gefa
    MODULE PROCEDURE GefaMatrix,GefaSpMatrix
  END INTERFACE
  INTERFACE MatVec
    MODULE PROCEDURE MatVecSp,MatVecSp1,MatVecMat
  END INTERFACE
  INTERFACE PSolve
    MODULE PROCEDURE PSolveSpMatrix,PSolveMatrix,PSolveSp
  END INTERFACE
  INTERFACE Assign
    MODULE PROCEDURE Assign_VectorVector
  END INTERFACE
CONTAINS

SUBROUTINE Assign_VectorVector(Vec1,Vec2)

  TYPE(Vec4_T), INTENT(OUT) :: Vec1
  TYPE(Vec4_T), INTENT(IN) :: Vec2


  Vec1%c=>Vec2%c
  Vec1%cInt=>Vec2%cInt
  Vec1%cB=>Vec2%cB
  Vec1%Rank=Vec2%Rank

END SUBROUTINE Assign_VectorVector

SUBROUTINE Mult_VecVec(Vec1,Vec2)

  TYPE(Vec4_T), INTENT(IN) :: Vec1
  TYPE(Vec4_T), INTENT(INOUT) :: Vec2

  Vec2%c=Vec1%c*Vec2%c

END SUBROUTINE Mult_VecVec

SUBROUTINE Mult_Vec4(Vec1,Scalar,Vec2,Vec3)
 
  REAL(RealKind) :: Scalar
  TYPE(Vec4_T), INTENT(IN) :: Vec1,Vec2
  TYPE(Vec4_T), INTENT(INOUT) :: Vec3
 
  WHERE (Vec2%c>Zero) 
    Vec3%c=Vec3%c*Vec1%c/(Scalar*Vec2%c)
  ELSEWHERE
    Vec3%c=Zero
  END WHERE 
 
END SUBROUTINE Mult_Vec4

SUBROUTINE Mult_Vec5(Vec1,Vec2,Scalar,Vec3,Vec4)

  REAL(RealKind) :: Scalar
  TYPE(Vec4_T), INTENT(IN) :: Vec1,Vec2,Vec3
  TYPE(Vec4_T), INTENT(INOUT) :: Vec4

  WHERE (Vec3%c>Zero)
    Vec4%c=Vec4%c*Vec1%c*Vec2%c/(Scalar*Vec3%c)
  ELSEWHERE
    Vec4%c=Zero
  END WHERE

END SUBROUTINE Mult_Vec5

SUBROUTINE Mult_Vec6(Vec1,Vec2,Scalar1,Scalar2,Vec3,Vec4)

  REAL(RealKind) :: Scalar1,Scalar2
  TYPE(Vec4_T), INTENT(IN) :: Vec1,Vec2,Vec3
  TYPE(Vec4_T), INTENT(INOUT) :: Vec4

  WHERE (Vec3%c>Zero)
    Vec4%c=Vec4%c*MAX(Vec1%c*Vec2%c/(Scalar2*Vec3%c),Zero)**Scalar1
  ELSEWHERE
    Vec4%c=Zero
  END WHERE

END SUBROUTINE Mult_Vec6

SUBROUTINE MultD_Vec6(Vec1,Vec2,Scalar1,Scalar2,Vec3,Vec4)

  REAL(RealKind) :: Scalar1,Scalar2
  TYPE(Vec4_T), INTENT(IN) :: Vec1,Vec2,Vec3
  TYPE(Vec4_T), INTENT(INOUT) :: Vec4

  WHERE (Vec3%c>Zero)
    Vec4%c=Vec4%c*MAX(Vec1%c*Vec2%c/(Scalar2*Vec3%c),1.d-10)**Scalar1
  ELSEWHERE
    Vec4%c=Zero
  END WHERE

END SUBROUTINE MultD_Vec6

SUBROUTINE Div_VecVec(Vec1,Vec2)

  TYPE(Vec4_T), INTENT(IN) :: Vec1
  TYPE(Vec4_T), INTENT(INOUT) :: Vec2

  WHERE (Vec1%c>Zero) 
    Vec2%c=Vec2%c/Vec1%c
  ELSEWHERE
    Vec2%c=Zero
  END WHERE 

END SUBROUTINE Div_VecVec

SUBROUTINE SetZero_VecVec(Vec1,Vec2)

  TYPE(Vec4_T), INTENT(IN) :: Vec1
  TYPE(Vec4_T), INTENT(INOUT) :: Vec2

  WHERE (Vec1%c<=Zero)
    Vec2%c=Zero
  END WHERE

END SUBROUTINE SetZero_VecVec

SUBROUTINE Add_Vector(alpha,Vec1,Vec2)

  REAL(RealKind) :: alpha
  TYPE(Vec4_T), INTENT(IN) :: Vec1
  TYPE(Vec4_T), INTENT(INOUT) :: Vec2

  Vec2%cInt=alpha*Vec1%cInt+Vec2%cInt
  Vec2%cB=alpha*Vec1%cB+Vec2%cB

END SUBROUTINE Add_Vector

SUBROUTINE Add_Vector_Vector(alpha,Vec1,Vec2)

  REAL(RealKind) :: alpha
  TYPE(Vec4_T), INTENT(IN) :: Vec1(:)
  TYPE(Vec4_T), INTENT(INOUT) :: Vec2(:)

  INTEGER :: i

  DO i=1,SIZE(Vec1(:))
    Vec2(i)%c=alpha*Vec1(i)%c+Vec2(i)%c
    Vec2(i)%cB=alpha*Vec1(i)%cB+Vec2(i)%cB
  END DO

END SUBROUTINE Add_Vector_Vector

SUBROUTINE Add3_Vector_Vector(alpha,Vec1,Vec2,Vec3)

  REAL(RealKind) :: alpha
  TYPE(Vec4_T), INTENT(IN) :: Vec1(:)
  TYPE(Vec4_T), INTENT(IN) :: Vec2(:)
  TYPE(Vec4_T), INTENT(INOUT) :: Vec3(:)

  INTEGER :: i

  DO i=1,SIZE(Vec1(:))
    Vec3(i)%c=alpha*Vec1(i)%c+Vec2(i)%c
    Vec3(i)%cB=alpha*Vec1(i)%cB+Vec2(i)%cB
  END DO

END SUBROUTINE Add3_Vector_Vector

SUBROUTINE Add_Vector_Partiell(alpha,Vec1,Vec2,i1,i2)

  REAL(RealKind) :: alpha
  TYPE(Vec4_T), INTENT(IN) :: Vec1
  TYPE(Vec4_T), INTENT(INOUT) :: Vec2
  INTEGER :: i1,i2

  Vec2%cInt(:,:,:,i1:i2)=alpha*Vec1%cInt(:,:,:,i1:i2)+Vec2%cInt(:,:,:,i1:i2)

END SUBROUTINE Add_Vector_Partiell


SUBROUTINE DivLeft_MatVector(Div,Vec1,Vec2)

  TYPE(Vec4_T), INTENT(INOUT) :: Div
  TYPE(Vec4_T), INTENT(IN) :: Vec1,Vec2

! Div=Vec1**-1*Vec2

  Div%Rank(1)=Vec1%Rank(2)
  Div%Rank(2)=Vec2%Rank(2)
! Div%c%c=Vec2%c%c/Vec1%c%c
  Div%cInt(:,:,:,1:MAX(Div%Rank(1),Div%Rank(2))) &
  =Vec2%cInt(:,:,:,1:MAX(Div%Rank(1),Div%Rank(2))) &
  /Vec1%cInt(:,:,:,1:MAX(Div%Rank(1),Div%Rank(2)))

END SUBROUTINE DivLeft_MatVector

SUBROUTINE DivRight_MatVector(Div,Vec1,Vec2)
 
  TYPE(Vec4_T), INTENT(INOUT) :: Div
  TYPE(Vec4_T), INTENT(IN) :: Vec1,Vec2
 
! Div=Vec1*Vec2**-1
 
  Div%Rank(1)=Vec1%Rank(1)
  Div%Rank(2)=Vec2%Rank(1)
  Div%cInt(:,:,:,1:MAX(Div%Rank(1),Div%Rank(2))) &
  =Vec1%cInt(:,:,:,1:MAX(Div%Rank(1),Div%Rank(2))) &
  /Vec2%cInt(:,:,:,1:MAX(Div%Rank(1),Div%Rank(2)))
 
END SUBROUTINE DivRight_MatVector

SUBROUTINE AddMult_MatVector(AddMult,Vec1,Vec2)

  TYPE(Vec4_T), INTENT(INOUT) :: AddMult
  TYPE(Vec4_T), INTENT(IN) :: Vec1,Vec2

! Mult=Vec1*Vec2

  INTEGER :: i,i1
  INTEGER :: Shift1

  AddMult%Rank(1)=Vec1%Rank(1)
  AddMult%Rank(2)=Vec2%Rank(2)

  Shift1=MIN(1,Vec1%Rank(1)-1)
  i1=1
  DO i=1,Vec1%Rank(2)
    AddMult%cInt(:,:,:,i1)=AddMult%cInt(:,:,:,i1)+Vec1%cInt(:,:,:,i)*Vec2%cInt(:,:,:,i)
    i1=i1+Shift1
  END DO

END SUBROUTINE AddMult_MatVector

SUBROUTINE SubMult_MatVector(SubMult,Vec1,Vec2)
 
  TYPE(Vec4_T), INTENT(INOUT) :: SubMult
  TYPE(Vec4_T), INTENT(IN) :: Vec1,Vec2
 
! Mult=Vec1*Vec2
 
  INTEGER :: i,i1,i2,i3
  INTEGER :: Shift1,Shift2,Shift3
  INTEGER :: MaxRank1,MaxRank2,MaxRank3,MaxRank
 
  SubMult%Rank(1)=Vec1%Rank(1)
  SubMult%Rank(2)=Vec2%Rank(2)
  MaxRank1=MAX(Vec1%Rank(1),Vec1%Rank(2))
  Shift1=MIN(1,MaxRank1-1)
  MaxRank2=MAX(Vec2%Rank(1),Vec2%Rank(2))
  Shift2=MIN(1,MaxRank2-1)
  MaxRank3=MAX(SubMult%Rank(1),SubMult%Rank(2))
  Shift3=MIN(1,MaxRank3-1)
  i1=1
  i2=1
  i3=1
  MaxRank=MAX(MaxRank1,MaxRank2,MaxRank3)
  DO i=1,MaxRank
    SubMult%cInt(:,:,:,i3)=SubMult%cInt(:,:,:,i3)-Vec1%cInt(:,:,:,i1)*Vec2%cInt(:,:,:,i2)
    i1=i1+Shift1
    i2=i2+Shift2
    i3=i3+Shift3
  END DO
 
END SUBROUTINE SubMult_MatVector

SUBROUTINE Neg_MatVector(Vec)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec


  Vec%cInt=-Vec%cInt

END SUBROUTINE Neg_MatVector

SUBROUTINE AddScalar_MatVector(Vec1,Vec2,Value)

  TYPE(Vec4_T) :: Vec1,Vec2
  REAL(RealKind)  :: Value

! Vec1=Vec2+Value

  Vec1%cInt=Vec2%cInt+Value

END SUBROUTINE AddScalar_MatVector

SUBROUTINE AddScalarPos_MatVector(Vec1,Vec2,Value,i)

  TYPE(Vec4_T) :: Vec1,Vec2
  REAL(RealKind)  :: Value
  INTEGER :: i

  INTEGER :: iPos

! Vec1=Value

  iPos=MIN(i,SIZE(Vec1%cInt(ixS,iyS,izS,:)))
  Vec1%c(ixS,iyS,izS,iPos)=Vec2%c(ixS,iyS,izS,iPos)+Value

END SUBROUTINE AddScalarPos_MatVector
SUBROUTINE Scale_Vector(Value,Vec)

  REAL(RealKind), INTENT(IN) :: Value
  TYPE(Vec4_T), INTENT(INOUT) :: Vec

  Vec%cInt=Value*Vec%cInt
  Vec%cB=Value*Vec%cB
END SUBROUTINE Scale_Vector

SUBROUTINE Deallocate_Vector(Vec)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec

  IF (ASSOCIATED(Vec%c)) THEN
    DEALLOCATE(Vec%c)
  END IF
  IF (ASSOCIATED(Vec%cB)) THEN
    DEALLOCATE(Vec%cB)
  END IF

END SUBROUTINE Deallocate_Vector

SUBROUTINE Allocate_Rank(Vec,Rank)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec
  INTEGER :: Rank(2)

  Vec%Rank=Rank
  ALLOCATE(Vec%c(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1:MAX(Rank(1),Rank(2))))
  Vec%c=Zero 
  Vec%cInt=>Vec%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1:MAX(Rank(1),Rank(2)))
  Vec%Vec3%c=>Vec%c(:,:,:,LBOUND(Vec%c,4):LBOUND(Vec%c,4))
  ALLOCATE(Vec%cB(0:0,0:0))

END SUBROUTINE Allocate_Rank

SUBROUTINE Allocate_Dim(Vec,it0,it1)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec
  INTEGER :: it0,it1

  ALLOCATE(Vec%c(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,it0+1:it1))
  Vec%c=Zero 
  Vec%cInt=>Vec%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,it0+1:it1)
  Vec%Rank(1)=1
  Vec%Rank(2)=it1-it0
  Vec%Vec3%c=>Vec%c(:,:,:,LBOUND(Vec%c,4):LBOUND(Vec%c,4))
  ALLOCATE(Vec%cB(0:0,0:0))

END SUBROUTINE Allocate_Dim

SUBROUTINE Allocate_Dim1(Vec,it)
                                                                                
  TYPE(Vec4_T), INTENT(INOUT) :: Vec
  INTEGER :: it
                                                                                
  ALLOCATE(Vec%c(1,1,1,it))
  Vec%c=Zero
  Vec%cInt=>Vec%c
  Vec%Rank(1)=1
  Vec%Rank(2)=1
  Vec%Vec3%c=>Vec%c(:,:,:,1:1)
  ALLOCATE(Vec%cB(0:0,0:0))
                                                                                
END SUBROUTINE Allocate_Dim1


SUBROUTINE Copy_VectorScalar(Vec,Scalar)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec
  REAL(RealKind), INTENT(IN) :: Scalar

  Vec%c=Scalar
  Vec%cB=Scalar

END SUBROUTINE Copy_VectorScalar

SUBROUTINE Copy_Vector_VectorScalar(Vec,Scalar)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec(:)
  REAL(RealKind), INTENT(IN) :: Scalar

  INTEGER :: i
  
  DO i=LBOUND(Vec,1),UBOUND(Vec,1)
    Vec(i)%c=Scalar
    Vec(i)%cB=Scalar
  END DO 

END SUBROUTINE Copy_Vector_VectorScalar

SUBROUTINE Copy_Vector_VectorVector(Vec1,Vec2)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec1(:)
  TYPE(Vec4_T), INTENT(IN) :: Vec2(:)

  INTEGER :: i
 
  DO ic=MAX(LBOUND(Vec1,1),LBOUND(Vec2,1)), &
        MIN(UBOUND(Vec1,1),UBOUND(Vec2,1))
    Vec1(i)%c=Vec2(i)%c
    Vec1(i)%cB=Vec2(i)%cB
    Vec1(i)%Rank=Vec2(i)%Rank
  END DO

END SUBROUTINE Copy_Vector_VectorVector

SUBROUTINE Copy_ScalarVector(Scalar,Vec)

  REAL(RealKind), INTENT(OUT) :: Scalar
  TYPE(Vec4_T), INTENT(IN) :: Vec

  Scalar=Vec%c(1,1,1,1)
END SUBROUTINE Copy_ScalarVector

SUBROUTINE Copy_VecVec(Vec1,Vec2)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec1
  TYPE(Vec4_T), INTENT(IN) :: Vec2

  Vec1%c=Vec2%c(:,:,:,LBOUND(Vec1%c,4):UBOUND(Vec1%c,4))
  Vec1%cB=Vec2%cB
  Vec1%Rank=Vec2%Rank
END SUBROUTINE Copy_VecVec

SUBROUTINE Copy_LowHigh(High,Low)

  TYPE(Vec4_T) ,INTENT(INOUT) :: High
  TYPE(Vec34_T) ,INTENT(IN) :: Low

  INTEGER ::i

  DO i=LBOUND(High%c,4)+1,UBOUND(High%c,4)
    High%c(:,:,:,i)=Low%c(:,:,:,LBOUND(Low%c,4))
  END DO

END SUBROUTINE Copy_LowHigh

FUNCTION NumberVec(x)

  INTEGER :: NumberVec
  TYPE(Vector4Cell_T) :: x

  NumberVec=0
  DO ic=LBOUND(x%Vec,1),UBOUND(x%Vec,1)
    NumberVec=NumberVec+SIZE(x%Vec(ic)%c,4)
  END DO
END FUNCTION NumberVec

SUBROUTINE MatVecSp(Ax,A,x)

  TYPE(Vector4Cell_T) :: Ax
  TYPE(JacSpMatrix4_T) :: A
! TYPE(SpDiag) :: A1,A2,A3
  TYPE(Vector4Cell_T) :: x

  REAL(RealKind) :: alpha

  DO ic=LBOUND(Ax%Vec,1),UBOUND(Ax%Vec,1)
    CALL VectorCellBounds(Ax%Vec(ic))
    IF (ic==uPosL.OR.ic==uPosR.OR. &
        ic==vPosL.OR.ic==vPosR.OR. &
        ic==wPosL.OR.ic==wPosR) THEN
      DO it=LBOUND(Ax%Vec(ic)%c,4),UBOUND(Ax%Vec(ic)%c,4)
        CALL SpAVec(Ax%Vec(ic)%c(:,iyc0,izc0,it),A%JactMom,x%Vec(ic)%c(:,iyc0,izc0,it))
      END DO
    ELSE IF (ic/=RhoRPos.OR.ic/=nrPos) THEN
      alpha=ScaleMat2(ic) 
      IF (alpha<Zero) THEN
        A%JacFall%Val=A%JacFall%Val*alpha
        DO it=LBOUND(Ax%Vec(ic)%c,4),UBOUND(Ax%Vec(ic)%c,4)
          CALL SpAVec(Ax%Vec(ic)%c(:,iyc0,izc0,it),A%JacTPot,A%JacFall &
                     ,x%Vec(ic)%c(:,iyc0,izc0,it))
        END DO
        A%JacFall%Val=A%JacFall%Val/alpha
      ELSE
        DO it=LBOUND(Ax%Vec(ic)%c,4),UBOUND(Ax%Vec(ic)%c,4)
          CALL SpAVec(Ax%Vec(ic)%c(:,iyc0,izc0,it),A%JacTPot,x%Vec(ic)%c(:,iyc0,izc0,it))
        END DO
      END IF
    ELSE IF (ic==RhoRpos.OR.ic==nrPos) THEN
      DO it=LBOUND(Ax%Vec(ic)%c,4),UBOUND(Ax%Vec(ic)%c,4)
        CALL SpAVec(Ax%Vec(ic)%c(:,iyc0,izc0,it),A%JacTPot,A%JacFallRhoL,x%Vec(ic)%c(:,iyc0,izc0,it))
      END DO
    ELSE
      Ax%Vec(ic)%c=x%Vec(ic)%c
    END IF
    Ax%Vec(ic)%cB=x%Vec(ic)%cB
  END DO
END SUBROUTINE MatVecSp

SUBROUTINE MatVecSp1(Ax,A)

  TYPE(Vector4Cell_T) :: Ax
  TYPE(SpDiag) :: A

  DO ic=LBOUND(Ax%Vec,1),UBOUND(Ax%Vec,1)
    CALL VectorCellBounds(Ax%Vec(ic))
    DO it=LBOUND(Ax%Vec(ic)%c,4),UBOUND(Ax%Vec(ic)%c,4)
      CALL SpAVec(Ax%Vec(ic)%c(:,iyc0,izc0,it),A)
    END DO
  END DO
END SUBROUTINE MatVecSp1


SUBROUTINE MatVecMat(Ax,A,x)

  TYPE(Vector4Cell_T) :: Ax
  TYPE(Matrix4Cell_T) :: A
  TYPE(Vector4Cell_T) :: x

  CALL VectorCellBounds(Ax%Vec(LBOUND(Ax%Vec,1)))
  DO ic=LBOUND(Ax%Vec,1),UBOUND(Ax%Vec,1)
    DO jc=LBOUND(Ax%Vec,1),UBOUND(Ax%Vec,1)
      Ax%Vec(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,:)= &
      Ax%Vec(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,:)- &
      dt*beta0*A%Mat(ic,jc)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,:)* &
      x%Vec(jc)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,:)
    END DO
  END DO
END SUBROUTINE MatVecMat

SUBROUTINE PSolveSp(Px,P,x)

  TYPE(Vector4Cell_T) :: Px
  TYPE(JacSpMatrix4_T) :: P
  TYPE(Vector4Cell_T) :: x

  REAL(RealKind) :: alpha

  CALL VectorCellBounds(Px%Vec(1))

  DO ic=LBOUND(Px%Vec,1),UBOUND(Px%Vec,1)
    CALL VectorCellBounds(Px%Vec(ic))
    IF (ic==uPosL.OR.ic==uPosR.OR. &
        ic==vPosL.OR.ic==vPosR.OR. &
        ic==wPosL.OR.ic==wPosR) THEN
      DO it=LBOUND(Px%Vec(ic)%c,4),UBOUND(Px%Vec(ic)%c,4)
        CALL GaussSeidel(Px%Vec(ic)%c(:,iyc0,izc0,it),P%JacTMom,x%Vec(ic)%c(:,iyc0,izc0,it))
      END DO
    ELSE IF (ic/=RhoRPos.OR.ic/=nrPos) THEN
      alpha=ScaleMat2(ic)
      IF (alpha<Zero) THEN
        P%JacFall%Val=P%JacFall%Val*alpha
        DO it=LBOUND(Px%Vec(ic)%c,4),UBOUND(Px%Vec(ic)%c,4)
          CALL GaussSeidel(Px%Vec(ic)%c(:,iyc0,izc0,it),P%JacTPot,P%JacFall, &
                           x%Vec(ic)%c(:,iyc0,izc0,it))
        END DO
        P%JacFall%Val=P%JacFall%Val/alpha
      ELSE
        DO it=LBOUND(Px%Vec(ic)%c,4),UBOUND(Px%Vec(ic)%c,4)
          CALL GaussSeidel(Px%Vec(ic)%c(:,iyc0,izc0,it),P%JacTPot,x%Vec(ic)%c(:,iyc0,izc0,it))
        END DO
      END IF
    ELSE IF (ic==RhoRpos.OR.ic==nrPos) THEN
      DO it=LBOUND(Px%Vec(ic)%c,4),UBOUND(Px%Vec(ic)%c,4)
        CALL GaussSeidel(Px%Vec(ic)%c(:,iyc0,izc0,it),P%JacTPot,P%JacFallRhoL, &
                         x%Vec(ic)%c(:,iyc0,izc0,it))
      END DO   
    END IF
    Px%Vec(ic)%c(:,:,izc0,:)=Zero
    Px%Vec(ic)%c(:,:,izc1+1,:)=Zero
    Px%Vec(ic)%c(:,iyc0,:,:)=Zero
    Px%Vec(ic)%c(:,iyc1+1,:,:)=Zero
    Px%Vec(ic)%c(ixc0,:,:,:)=Zero
    Px%Vec(ic)%c(ixc1+1,:,:,:)=Zero
  END DO

END SUBROUTINE PSolveSp

SUBROUTINE PSolveMatrix(P,x)

  TYPE(Matrix4Cell_T) :: P
  TYPE(Vector4Cell_T) :: x

  REAL(RealKind) :: Fac
  TYPE (Vec4_T), POINTER :: A(:,:)
  TYPE (Vec4_T), POINTER :: b(:)

  A=>P%Mat
  b=>x%Vec
  CALL VectorCellBounds(b(LBOUND(b,1)))
  DO kc=LBOUND(b,1),UBOUND(b,1)
    DO ic=kc+1,UBOUND(b,1)
      DO it=LBOUND(b(ic)%c,4),UBOUND(b(ic)%c,4)
        b(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,it:it)= &
        b(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,it:it)- &
        A(ic,kc)%c*b(kc)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,it:it)
      END DO
    END DO
  END DO
  DO kc=UBOUND(b,1),LBOUND(b,1),-1
    DO it=LBOUND(b(kc)%c,4),UBOUND(b(kc)%c,4)
      b(kc)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,it:it)= &
      b(kc)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,it:it)/A(kc,kc)%c
    END DO
    DO ic=LBOUND(b,1),kc-1
      DO it=LBOUND(b(ic)%c,4),UBOUND(b(ic)%c,4)
        b(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,it:it)= &
        b(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,it:it)- &
        A(ic,kc)%c*b(kc)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,it:it)
      END DO
    END DO
  END DO
END SUBROUTINE PSolveMatrix

SUBROUTINE GefaSpMatrix(MatrixMatVectorLU)

  TYPE(SpMatrix4Cell_T) :: MatrixMatVectorLU

  INTEGER :: n
  TYPE (Vec4_T), POINTER :: LU(:)
  TYPE (Vec4_T) :: w(MatrixMatVectorLU%n)
  TYPE (Vec4_T) :: alpha
  INTEGER, POINTER :: RowPtr(:),ColInd(:),DiagPtr(:)
  INTEGER :: i,j,jj,jjLU,kk
  INTEGER :: Rank(2)

  MatrixMatVectorLU%Factor=.TRUE.
  n=MatrixMatVectorLU%n
  LU=>MatrixMatVectorLU%Mat
  RowPtr=>MatrixMatVectorLU%Struct%RowPtr(:)
  DiagPtr=>MatrixMatVectorLU%Struct%DiagPtr(:)
  ColInd=>MatrixMatVectorLU%Struct%ColInd(:)

  Rank(1:2)=(/MatrixMatVectorLU%nMaxVec,MatrixMatVectorLU%nMaxVec/)
  CALL Allocate(alpha,Rank)
  DO i=1,SIZE(ColInd)
    CALL ScaleV(-dt*beta0,LU(i))
  END DO
  DO i=1,n
    CALL AddScalar(LU(DiagPtr(i)),LU(DiagPtr(i)),One)
  END DO

  DO i=1,n
    DO jj=RowPtr(i),RowPtr(i+1)-1
      CALL Assign(w(ColInd(jj)),LU(jj))
    END DO
    DO jj=RowPtr(i),DiagPtr(i)-1
      j=ColInd(jj)
      CALL DivRight_MatVector(alpha,w(j),LU(DiagPtr(j)))
!     CALL Copy_MatVector(w(j),alpha)
      w(j)=alpha
      DO kk=DiagPtr(j)+1,RowPtr(j+1)-1
        CALL SubMult_MatVector(w(ColInd(kk)),alpha,LU(kk))
      END DO
    END DO
  END DO

  CALL Deallocate(alpha)

END SUBROUTINE GefaSpMatrix

SUBROUTINE PSolveSpMatrix(MatrixMatVectorLU,bVec)

  TYPE(SpMatrix4Cell_T) :: MatrixMatVectorLU
  TYPE(Vector4Cell_T) :: bVec

  TYPE (Vec4_T) :: w(MatrixMatVectorLU%n)

  INTEGER :: n
  TYPE (Vec4_T), POINTER :: LU(:)
  TYPE(Vec4_T), POINTER :: b(:)
  INTEGER, POINTER :: RowPtr(:),ColInd(:),DiagPtr(:)
  INTEGER, POINTER :: Permu(:),InvPer(:)
  INTEGER :: i,ii,j,jj,jjLU,kk

 
  IF (.NOT.MatrixMatVectorLU%Factor) THEN
    CALL GefaSpMatrix(MatrixMatVectorLU)
  END IF
    
  n=MatrixMatVectorLU%n
  LU=>MatrixMatVectorLU%Mat
  b=>bVec%Vec
  RowPtr=>MatrixMatVectorLU%Struct%RowPtr(:)
  DiagPtr=>MatrixMatVectorLU%Struct%DiagPtr(:)
  ColInd=>MatrixMatVectorLU%Struct%ColInd(:)
  Permu=>MatrixMatVectorLU%Struct%Permu(:)
  InvPer=>MatrixMatVectorLU%Struct%InvPer(:)

  DO i=1,n
    CALL Assign(w(Permu(i)),b(i))
  END DO
  DO i=2,n
    DO jj=RowPtr(i),DiagPtr(i)-1
      j=ColInd(jj)
      CALL SubMult_MatVector(w(i),LU(jj),w(j))
    END DO
  END DO
  DO i=n,1,-1
    DO jj=DiagPtr(i)+1,RowPtr(i+1)-1
      j=ColInd(jj)
      CALL SubMult_MatVector(w(i),LU(jj),w(j))
    END DO
    CALL DivLeft_MatVector(w(i),LU(DiagPtr(i)),w(i))
  END DO
END SUBROUTINE PSolveSpMatrix


SUBROUTINE Mat2(MatrixMatVectorLU,xVec,yVec)

  TYPE(SpMatrix4Cell_T) :: MatrixMatVectorLU
  TYPE(Vector4Cell_T) :: xVec,yVec

  TYPE (Vec4_T) :: w(MatrixMatVectorLU%n)
  TYPE (Vec4_T) :: w1(MatrixMatVectorLU%n)

  INTEGER :: n
  TYPE (Vec4_T), POINTER :: LU(:)
  TYPE(Vec4_T), POINTER :: x(:),y(:)
  INTEGER, POINTER :: RowPtr(:),ColInd(:),DiagPtr(:)
  INTEGER, POINTER :: Permu(:),InvPer(:)
  INTEGER :: i,ii,j,jj,jjLU,kk

 
  n=MatrixMatVectorLU%n
  LU=>MatrixMatVectorLU%Mat
  x=>xVec%Vec
  y=>yVec%Vec
  RowPtr=>MatrixMatVectorLU%Struct%RowPtr(:)
  DiagPtr=>MatrixMatVectorLU%Struct%DiagPtr(:)
  ColInd=>MatrixMatVectorLU%Struct%ColInd(:)
  Permu=>MatrixMatVectorLU%Struct%Permu(:)
  InvPer=>MatrixMatVectorLU%Struct%InvPer(:)

  yVec=Zero
  DO i=1,n
    CALL Assign(w(Permu(i)),x(i))
    CALL Assign(w1(Permu(i)),y(i))
  END DO
  DO i=1,n
    DO jj=RowPtr(i),RowPtr(i+1)-1
      j=ColInd(jj)
      CALL SubMult_MatVector(w1(j),LU(jj),w(i))
    END DO
  END DO
END SUBROUTINE Mat2



SUBROUTINE AllocateVectorCell1(VectorCell,VectorComponents)

  TYPE (Vector4Cell_T) :: VectorCell
  INTEGER :: VectorComponents
  ALLOCATE(VectorCell%Vec(VectorComponents))
  DO ic=1,VectorComponents
    ALLOCATE(VectorCell%Vec(ic)%c(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1))
    VectorCell%Vec(ic)%cInt=>VectorCell%Vec(ic)%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1:1)
    VectorCell%Vec(ic)%c=0.0e0
    ALLOCATE(VectorCell%Vec(ic)%cB(0,0))
  END DO

END SUBROUTINE AllocateVectorCell1

SUBROUTINE AllocateVectorCell2(VectorCell,VectorComponents1,VectorComponents2)

  TYPE (Vector4Cell_T) :: VectorCell
  INTEGER :: VectorComponents1,VectorComponents2
  ALLOCATE(VectorCell%Vec(VectorComponents1:VectorComponents2))
  DO ic=LBOUND(VectorCell%Vec,1),UBOUND(VectorCell%Vec,1)
    ALLOCATE(VectorCell%Vec(ic)%c(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1))
    VectorCell%Vec(ic)%cInt=>VectorCell%Vec(ic)%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1:1)
    ALLOCATE(VectorCell%Vec(ic)%cB(0:0,0:0))
    VectorCell%Vec(ic)%c=0.0e0
  END DO

END SUBROUTINE AllocateVectorCell2

SUBROUTINE AllocateVectorCellVectorCell(VectorCell,VectorCellO,LBound0,UBound0)

  TYPE (Vector4Cell_T) :: VectorCell
  TYPE (Vector4Cell_T) :: VectorCellO
  INTEGER, OPTIONAL :: LBound0,UBound0

  INTEGER :: LBound0Loc,UBound0Loc

  IF (PRESENT(LBound0)) THEN
    LBound0Loc=LBound0
  ELSE
    LBound0Loc=LBOUND(VectorCellO%Vec,1)
  END IF
  IF (PRESENT(UBound0)) THEN
    UBound0Loc=UBound0
  ELSE
    UBound0Loc=UBOUND(VectorCellO%Vec,1)
  END IF
  ALLOCATE(VectorCell%Vec(LBound0Loc:UBound0Loc))
  DO ic=LBound0Loc,UBound0Loc
    CALL Allocate(VectorCell%Vec(ic),VectorCellO%Vec(ic))
    VectorCell%Vec(ic)%c=0.0e0
  END DO

END SUBROUTINE AllocateVectorCellVectorCell

SUBROUTINE Allocate_Vector(Vec1,Vec2)

  TYPE(Vec4_T), INTENT(INOUT) :: Vec1
  TYPE(Vec4_T), INTENT(IN) :: Vec2

  ALLOCATE(Vec1%c(LBOUND(Vec2%c,1):UBOUND(Vec2%c,1) &
                 ,LBOUND(Vec2%c,2):UBOUND(Vec2%c,2) &
                 ,LBOUND(Vec2%c,3):UBOUND(Vec2%c,3) &
                 ,LBOUND(Vec2%c,4):UBOUND(Vec2%c,4)))
  Vec1%cInt=>Vec1%c(LBOUND(Vec2%c,1)+1:UBOUND(Vec2%c,1)-1 &
                 ,LBOUND(Vec2%c,2)+1:UBOUND(Vec2%c,2)-1 &
                 ,LBOUND(Vec2%c,3)+1:UBOUND(Vec2%c,3)-1 &
                 ,LBOUND(Vec2%c,4):UBOUND(Vec2%c,4))
  Vec1%c=Zero
  Vec1%Rank=Vec2%Rank
  Vec1%Vec3%c=>Vec1%c(:,:,:,LBOUND(Vec1%c,4):LBOUND(Vec1%c,4))
  IF (ASSOCIATED(Vec2%cB)) THEN
    ALLOCATE(Vec1%cB(LBOUND(Vec2%cB,1):UBOUND(Vec2%cB,1) &
                 ,LBOUND(Vec2%cB,2):UBOUND(Vec2%cB,2)))
    Vec1%cB=Zero
  END IF 

END SUBROUTINE Allocate_Vector

SUBROUTINE AllocateMatrixCell(MatrixCell,VectorComponents)

  TYPE (Matrix4Cell_T) :: MatrixCell
  INTEGER :: VectorComponents
  ALLOCATE(MatrixCell%Mat(VectorComponents, &
                           VectorComponents))
  DO ic=1,VectorComponents
    DO jc=1,VectorComponents
      ALLOCATE(MatrixCell%Mat(ic,jc)%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1))
      MatrixCell%Mat(ic,jc)%c=0.0e0
    END DO
  END DO

END SUBROUTINE AllocateMatrixCell

SUBROUTINE AllocateMatrixCellVectorCell(MatrixCell,VectorCellO)

  TYPE (Matrix4Cell_T) :: MatrixCell
  TYPE (Vector4Cell_T) :: VectorCellO
  ALLOCATE(MatrixCell%Mat(LBOUND(VectorCellO%Vec,1):UBOUND(VectorCellO%Vec,1)  &
                         ,LBOUND(VectorCellO%Vec,1):UBOUND(VectorCellO%Vec,1))) 
  DO ic=LBOUND(VectorCellO%Vec,1),UBOUND(VectorCellO%Vec,1)
    DO jc=LBOUND(VectorCellO%Vec,1),UBOUND(VectorCellO%Vec,1)
      CALL Allocate(MatrixCell%Mat(ic,jc),VectorCellO%Vec(1))
    END DO
  END DO

END SUBROUTINE AllocateMatrixCellVectorCell

SUBROUTINE SpDiagAllocate(A,NumDiag)

  TYPE(SpDiag), POINTER :: A
  INTEGER, OPTIONAL :: NumDiag

  IF (.NOT.PRESENT(NumDiag)) THEN
    NumDiag=7
  END IF
  IF (NumDiag==7) THEN
    A%m=(nx+2)*(ny+2)*(nz+2)
    A%n=(nx+2)*(ny+2)*(nz+2)
    A%NumDiag=7
    IF (.NOT.ASSOCIATED(A%DiagPtr)) THEN
      ALLOCATE(A%DiagPtr(A%NumDiag))
      A%DiagPtr(1)=-(nx+2)*(ny+2)
      A%DiagPtr(2)=-(nx+2)
      A%DiagPtr(3)=-1
      A%DiagPtr(4)=0
      A%DiagPtr(5)=1
      A%DiagPtr(6)=(nx+2)
      A%DiagPtr(7)=(nx+2)*(ny+2)
    END IF
    IF (.NOT.ASSOCIATED(A%Val)) THEN
      ALLOCATE(A%Val(A%m,A%NumDiag))
      A%Val=Zero
    END IF
    A%Diag=4
  ELSE IF (NumDiag==2) THEN
    A%m=(nx+2)*(ny+2)*(nz+2)
    A%n=(nx+2)*(ny+2)*(nz+2)
    A%NumDiag=2
    IF (.NOT.ASSOCIATED(A%DiagPtr)) THEN
      ALLOCATE(A%DiagPtr(A%NumDiag))
      A%DiagPtr(1)=0
      A%DiagPtr(2)=(nx+2)*(ny+2)
    END IF
    IF (.NOT.ASSOCIATED(A%Val)) THEN
      ALLOCATE(A%Val(A%m,A%NumDiag))
      A%Val=Zero
    END IF
    A%Diag=1
  END IF

END SUBROUTINE SpDiagAllocate

SUBROUTINE DeallocateVectorCell(VectorCell)

  TYPE (Vector4Cell_T) :: VectorCell

  DO ic=LBOUND(VectorCell%Vec,1),UBOUND(VectorCell%Vec,1)
    DEALLOCATE(VectorCell%Vec(ic)%c)
    IF (ASSOCIATED(VectorCell%Vec(ic)%cB)) THEN
      DEALLOCATE(VectorCell%Vec(ic)%cB)
    END IF
  END DO
  DEALLOCATE(VectorCell%Vec)

END SUBROUTINE DeallocateVectorCell

SUBROUTINE DeallocateMatrixCell(MatrixCell)

  TYPE (Matrix4Cell_T) :: MatrixCell

  DO ic=LBOUND(MatrixCell%Mat,1),UBOUND(MatrixCell%Mat,1)
    DO jc=LBOUND(MatrixCell%Mat,2),UBOUND(MatrixCell%Mat,2)
      DEALLOCATE(MatrixCell%Mat(ic,jc)%c)
    END DO
  END DO
  DEALLOCATE(MatrixCell%Mat)

END SUBROUTINE DeallocateMatrixCell
SUBROUTINE VectorCellBounds(Vec4)

  TYPE (Vec4_T) :: Vec4

  ixc0=LBOUND(Vec4%c,1)
  ixc1=UBOUND(Vec4%c,1)-1
  iyc0=LBOUND(Vec4%c,2)
  iyc1=UBOUND(Vec4%c,2)-1
  izc0=LBOUND(Vec4%c,3)
  izc1=UBOUND(Vec4%c,3)-1
  itc0=LBOUND(Vec4%c,4)
  itc1=UBOUND(Vec4%c,4)

END SUBROUTINE VectorCellBounds

FUNCTION Dot_Vec4(Vector41,Vector42)

  REAL(RealKind) :: Dot_Vec4
  TYPE (Vec4_T) :: Vector41(:),Vector42(:)

  Dot_Vec4=Zero
  DO ic=MAX(LBOUND(Vector41,1),1),UBOUND(Vector41,1)
    CALL VectorCellBounds(Vector41(ic))
    Dot_Vec4=Dot_Vec4+ &
       SUM(Vector41(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,:)* &
           Vector42(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,:)) 
  END DO
END FUNCTION Dot_Vec4

FUNCTION Dot_VectorCell(VectorCell1,VectorCell2)

  REAL(RealKind) :: Dot_VectorCell
  TYPE (Vector4Cell_T) :: VectorCell1,VectorCell2

  Dot_VectorCell=Zero
  DO ic=MAX(LBOUND(VectorCell1%Vec,1),1),UBOUND(VectorCell1%Vec,1)
    CALL VectorCellBounds(VectorCell1%Vec(ic))
    Dot_VectorCell=Dot_VectorCell+ &
       SUM(VectorCell1%Vec(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,:)* &
           VectorCell2%Vec(ic)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,:)) 
  END DO
END FUNCTION Dot_VectorCell

FUNCTION Sum_VectorCell(VectorCell1,Pos)

  REAL(RealKind) :: Sum_VectorCell
  TYPE (Vector4Cell_T) :: VectorCell1
  INTEGER :: Pos

  CALL VectorCellBounds(VectorCell1%Vec(Pos))
  Sum_VectorCell= SUM(VolC(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1)*VectorCell1%Vec(Pos)%c(ixc0+1:ixc1,iyc0+1:iyc1,izc0+1:izc1,1))
END FUNCTION Sum_VectorCell

FUNCTION Error_VectorCell(ErrorInd,VectorCell,VectorCell1,ATol,RTol,LenVec,V2)

  REAL(RealKind) :: Error_VectorCell
  TYPE (Vector4Cell_T) :: ErrorInd,VectorCell,VectorCell1,V2
  REAL(RealKind) :: ATol(:),RTol(:)
  REAL(RealKind) :: LenVec

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp
  REAL(RealKind) :: MaxTemp

  Error_VectorCell=Zero
  MaxTemp=0.0d0
  ixM=ixc0+1
  iyM=iyc0+1
  izM=izc0+1
  icM=1
  DO ic=MAX(LBOUND(VectorCell%Vec,1),1),UBOUND(VectorCell%Vec,1)
    CALL VectorCellBounds(VectorCell%Vec(ic))
    DO iz=izc0+1,izc1
      DO iy=iyc0+1,iyc1
        DO ix=ixc0+1,ixc1
          DO it=LBOUND(VectorCell%Vec(ic)%c,4),UBOUND(VectorCell%Vec(ic)%c,4)
            LenVec=LenVec+1
            Temp=VolC(ix,iy,iz)*ErrorInd%Vec(ic)%c(ix,iy,iz,it)/ &
             (ATol(ic)+RTol(ic)*ABS(VectorCell%Vec(ic)%c(ix,iy,iz,it)))
            IF (ABS(Temp)>MaxTemp) THEN
              ixM=ix
              iyM=iy
              izM=iz  
              icM=ic
              itM=it
              MaxTemp=ABS(Temp)
            END IF
            Error_VectorCell=Error_VectorCell+ABS(Temp) 
            LenVec=LenVec+VolC(ix,iy,iz) 
         END DO
       END DO
      END DO
    END DO
  END DO
   WRITE(*,*) 'Error',ixM,iyM,izM,icM
   WRITE(*,*) ErrorInd%Vec(icM)%c(ixM,iyM,izM,itM)
   WRITE(*,*) VectorCell%Vec(icM)%c(ixM,iyM,izM,itM),VectorCell1%Vec(icM)%c(ixM,iyM,izM,itM)
   WRITE(*,*) V2%Vec(icM)%c(ixM,iyM,izM,itM)
END FUNCTION Error_VectorCell

FUNCTION ErrorMax_VectorCell(ErrorInd,VectorCell,VectorCell1,ATol,RTol,V2)

  REAL(RealKind) :: ErrorMax_VectorCell
  TYPE (Vector4Cell_T) :: ErrorInd,VectorCell,VectorCell1,V2
  REAL(RealKind) :: ATol(:),RTol(:)

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp
  REAL(RealKind) :: MaxTemp,MaxVol

  ErrorMax_VectorCell=Zero
  MaxTemp=0.0d0
  MaxVol=0.0d0
  ixM=ixc0+1
  iyM=iyc0+1
  izM=izc0+1
  icM=1
  DO ic=LBOUND(VectorCell%Vec,1),UBOUND(VectorCell%Vec,1)
    CALL VectorCellBounds(VectorCell%Vec(ic))
    DO iz=izc0+1,izc1
      DO iy=iyc0+1,iyc1
        DO ix=ixc0+1,ixc1
          MaxVol=MAX(MaxVol,VolC(ix,iy,iz))
          DO it=LBOUND(VectorCell%Vec(ic)%c,4),UBOUND(VectorCell%Vec(ic)%c,4)
            Temp=VolC(ix,iy,iz)*ErrorInd%Vec(ic)%c(ix,iy,iz,it)/ &
             (ATol(ic)+RTol(ic)*ABS(VectorCell%Vec(ic)%c(ix,iy,iz,it)))
            IF (ABS(Temp)>MaxTemp) THEN
              ixM=ix
              iyM=iy
              izM=iz  
              icM=ic
              itM=it
              MaxTemp=ABS(Temp)
            END IF
         END DO
       END DO
      END DO
    END DO
  END DO
  ErrorMax_VectorCell=MaxTemp/MaxVol
   WRITE(*,*) 'Error',ixM,iyM,izM,icM
   WRITE(*,*) ErrorInd%Vec(icM)%c(ixM,iyM,izM,itM)
   WRITE(*,*) VectorCell%Vec(icM)%c(ixM,iyM,izM,itM),VectorCell1%Vec(icM)%c(ixM,iyM,izM,itM)
   WRITE(*,*) V2%Vec(icM)%c(ixM,iyM,izM,itM)
END FUNCTION ErrorMax_VectorCell

SUBROUTINE Axpy_VectorCell(alpha,VectorCell1,VectorCell2)

  REAL(RealKind) :: alpha
  TYPE (Vector4Cell_T) :: VectorCell1,VectorCell2
  INTEGER :: n

  IF (alpha/=Zero) THEN
    DO ic=MAX(LBOUND(VectorCell1%Vec,1),LBOUND(VectorCell2%Vec,1)), &
          MIN(UBOUND(VectorCell1%Vec,1),UBOUND(VectorCell2%Vec,1))
      VectorCell2%Vec(ic)%c=alpha*VectorCell1%Vec(ic)%c+VectorCell2%Vec(ic)%c
      VectorCell2%Vec(ic)%cB=alpha*VectorCell1%Vec(ic)%cB+VectorCell2%Vec(ic)%cB
    END DO
  END IF

END SUBROUTINE Axpy_VectorCell

SUBROUTINE Ax1mx2py_VectorCell(alpha,VectorCell1,VectorCell2,VectorCell3)

  REAL(RealKind) :: alpha
  TYPE (Vector4Cell_T) :: VectorCell1,VectorCell2,VectorCell3
  INTEGER :: n

  IF (alpha/=Zero) THEN
    DO ic=MAX(LBOUND(VectorCell1%Vec,1),LBOUND(VectorCell2%Vec,1),LBOUND(VectorCell3%Vec,1)), &
          MIN(UBOUND(VectorCell1%Vec,1),UBOUND(VectorCell2%Vec,1),UBOUND(VectorCell3%Vec,1))
      VectorCell3%Vec(ic)%c=alpha*(VectorCell1%Vec(ic)%c-VectorCell2%Vec(ic)%c) &
                           +VectorCell3%Vec(ic)%c
      VectorCell3%Vec(ic)%cB=alpha*(VectorCell1%Vec(ic)%cB-VectorCell2%Vec(ic)%cB) &
                            +VectorCell3%Vec(ic)%cB
    END DO
  END IF

END SUBROUTINE Ax1mx2py_VectorCell

SUBROUTINE Ax1px2py_VectorCell(alpha,VectorCell1,VectorCell2,VectorCell3)

  REAL(RealKind) :: alpha
  TYPE (Vector4Cell_T) :: VectorCell1,VectorCell2,VectorCell3
  INTEGER :: n

  IF (alpha/=Zero) THEN
    DO ic=MAX(LBOUND(VectorCell1%Vec,1),LBOUND(VectorCell2%Vec,1),LBOUND(VectorCell3%Vec,1)), &
          MIN(UBOUND(VectorCell1%Vec,1),UBOUND(VectorCell2%Vec,1),UBOUND(VectorCell3%Vec,1))
      VectorCell3%Vec(ic)%c=alpha*(VectorCell1%Vec(ic)%c+VectorCell2%Vec(ic)%c) &
                           +VectorCell3%Vec(ic)%c
      VectorCell3%Vec(ic)%cB=alpha*(VectorCell1%Vec(ic)%cB+VectorCell2%Vec(ic)%cB) &
                            +VectorCell3%Vec(ic)%cB
    END DO
  END IF

END SUBROUTINE Ax1px2py_VectorCell


SUBROUTINE Xpy_VectorCell(VectorCell1,VectorCell2)

  TYPE (Vector4Cell_T) :: VectorCell1,VectorCell2

  DO ic=MAX(LBOUND(VectorCell1%Vec,1),LBOUND(VectorCell2%Vec,1)), &
        MIN(UBOUND(VectorCell1%Vec,1),UBOUND(VectorCell2%Vec,1))
    VectorCell2%Vec(ic)%c=VectorCell1%Vec(ic)%c+VectorCell2%Vec(ic)%c
    VectorCell2%Vec(ic)%cB=VectorCell1%Vec(ic)%cB+VectorCell2%Vec(ic)%cB
  END DO

END SUBROUTINE Xpy_VectorCell


SUBROUTINE Axpby_VectorCell(alpha,VectorCell1,beta,VectorCell2)

  REAL(RealKind) :: alpha,beta
  TYPE (Vector4Cell_T) :: VectorCell1,VectorCell2

  DO ic=MAX(LBOUND(VectorCell1%Vec,1),LBOUND(VectorCell2%Vec,1)), &
        MIN(UBOUND(VectorCell1%Vec,1),UBOUND(VectorCell2%Vec,1))
    VectorCell2%Vec(ic)%c=alpha*VectorCell1%Vec(ic)%c+beta*VectorCell2%Vec(ic)%c
    VectorCell2%Vec(ic)%cB=alpha*VectorCell1%Vec(ic)%cB+beta*VectorCell2%Vec(ic)%cB
  END DO

END SUBROUTINE Axpby_VectorCell

SUBROUTINE Xpay_VectorCell(VectorCell1,alpha,VectorCell2)

  REAL(RealKind) :: alpha
  TYPE (Vector4Cell_T) :: VectorCell1,VectorCell2

  DO ic=MAX(LBOUND(VectorCell1%Vec,1),LBOUND(VectorCell2%Vec,1)), &
        MIN(UBOUND(VectorCell1%Vec,1),UBOUND(VectorCell2%Vec,1))
    VectorCell2%Vec(ic)%c=alpha*VectorCell2%Vec(ic)%c+VectorCell1%Vec(ic)%c
    VectorCell2%Vec(ic)%cB=alpha*VectorCell2%Vec(ic)%cB+VectorCell1%Vec(ic)%cB
  END DO

END SUBROUTINE Xpay_VectorCell
SUBROUTINE Copy_VectorCell(VectorCell1,VectorCell2)

  TYPE (Vector4Cell_T) :: VectorCell1,VectorCell2

  DO ic=MAX(LBOUND(VectorCell1%Vec,1),LBOUND(VectorCell2%Vec,1)), &
        MIN(UBOUND(VectorCell1%Vec,1),UBOUND(VectorCell2%Vec,1))
    VectorCell2%Vec(ic)%c=VectorCell1%Vec(ic)%c
    IF (SIZE(VectorCell2%Vec(ic)%cB)==SIZE(VectorCell1%Vec(ic)%cB)) THEN
      VectorCell2%Vec(ic)%cB=VectorCell1%Vec(ic)%cB
    END IF
  END DO

END SUBROUTINE Copy_VectorCell
SUBROUTINE Scale_VectorCell(alpha,VectorCell1)

  REAL(RealKind) :: alpha
  TYPE (Vector4Cell_T) :: VectorCell1

  DO ic=LBOUND(VectorCell1%Vec,1),UBOUND(VectorCell1%Vec,1)
    VectorCell1%Vec(ic)%c=alpha*VectorCell1%Vec(ic)%c
    VectorCell1%Vec(ic)%cB=alpha*VectorCell1%Vec(ic)%cB
  END DO

END SUBROUTINE Scale_VectorCell
SUBROUTINE Value_VectorCell(VectorCell1,Value)

  TYPE(Vector4Cell_T), INTENT(OUT)  :: VectorCell1
  REAL(RealKind), INTENT(IN)  :: Value

  DO ic=LBOUND(VectorCell1%Vec,1),UBOUND(VectorCell1%Vec,1)
    VectorCell1%Vec(ic)%c=Value
    VectorCell1%Vec(ic)%cB=Value
  END DO

END SUBROUTINE Value_VectorCell

SUBROUTINE Value_MatrixCell(MatrixCell,Value)

  TYPE(Matrix4Cell_T), INTENT(INOUT)  :: MatrixCell
  REAL(RealKind), INTENT(IN)  :: Value

  DO ic=LBOUND(MatrixCell%Mat,1),UBOUND(MatrixCell%Mat,1)
    DO jc=LBOUND(MatrixCell%Mat,2),UBOUND(MatrixCell%Mat,2)
      MatrixCell%Mat(ic,jc)%c=Value
    END DO
  END DO

END SUBROUTINE Value_MatrixCell


SUBROUTINE Value_SpMatrixCell(SpMatrixCell,Value)

  TYPE(SpMatrix4Cell_T), INTENT(INOUT) :: SpMatrixCell
  REAL(RealKind), INTENT(IN) :: Value

  INTEGER :: i

  DO i=1,SIZE(SpMatrixCell%Mat)
    SpMatrixCell%Mat(i)=Value
  END DO
  SpMatrixCell%Factor=.FALSE.

END SUBROUTINE Value_SpMatrixCell

SUBROUTINE Random_VectorCell(VectorCell1)

  IMPLICIT NONE

  TYPE(Vector4Cell_T), INTENT(OUT)  :: VectorCell1

  DO ic=LBOUND(VectorCell1%Vec,1),UBOUND(VectorCell1%Vec,1)
    CALL RANDOM_NUMBER(VectorCell1%Vec(ic)%c)
  END DO
END SUBROUTINE Random_VectorCell


SUBROUTINE GefaMatrix(MatrixCellLU)

  TYPE(Matrix4Cell_T) :: MatrixCellLU

  REAL(RealKind) :: Fac
  TYPE (Vec4_T), POINTER :: LU(:,:)

  LU=>MatrixCellLU%Mat
  Fac=dt*beta0
  DO ic=LBOUND(LU,1),UBOUND(LU,1)
    DO jc=1,UBOUND(LU,2)
      LU(ic,jc)%c=-Fac*LU(ic,jc)%c
    END DO
    LU(ic,ic)%c=One+LU(ic,ic)%c
  END DO
  DO kc=LBOUND(LU,2),UBOUND(LU,2)
    DO ic=kc+1,UBOUND(LU,1)
      LU(ic,kc)%c=LU(ic,kc)%c/LU(kc,kc)%c
      DO jc=kc+1,UBOUND(LU,2)
        LU(ic,jc)%c=LU(ic,jc)%c-LU(ic,kc)%c*LU(kc,jc)%c
      END DO
    END DO
  END DO

END SUBROUTINE GefaMatrix
END MODULE Vector4Cell_Mod

