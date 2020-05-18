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
