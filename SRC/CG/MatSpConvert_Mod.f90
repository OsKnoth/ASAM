MODULE MatSpConvert_Mod

  USE Kind_Mod
  USE MatSpDiag_Mod
  USE MatSpInd_Mod
  USE MatSpRowCol_Mod  
  USE MatSpIndMat_Mod
  USE MatSpRowColMat_Mod  
  USE MatSpRowColL_Mod  
  USE MatSpRowColD_Mod  
  USE MatSpRowColDiag_Mod  

  IMPLICIT NONE

CONTAINS

SUBROUTINE RowColD_To_RowColDiag(A,B)

  TYPE(SpRowColDiag), INTENT(INOUT) :: A
  TYPE(SpRowColD), INTENT(IN) :: B

  INTEGER :: i,j,jj,nzrA

  A%n=B%n
  A%m=B%m

  IF (.NOT.ASSOCIATED(A%RowPtr)) THEN
    ALLOCATE(A%RowPtr(A%n+1))
  END IF
  IF (.NOT.ASSOCIATED(A%DiagPtr)) THEN
    ALLOCATE(A%DiagPtr(A%n+1))
  END IF
  A%RowPtr(1)=1
  NzrA=0
  DO i=1,A%n
    A%RowPtr(i+1)=A%RowPtr(i)
    DO jj=B%RowPtr(1,i),B%RowPtr(2,i)
      j=B%ColInd(jj) 
      IF (i==j) THEN
        A%DiagPtr(i)=A%RowPtr(i+1)
      END IF
      A%RowPtr(i+1)=A%RowPtr(i+1)+1
      NzrA=NzrA+1
    END DO
  END DO
  IF (.NOT.ASSOCIATED(A%ColInd)) THEN
    ALLOCATE(A%ColInd(NzrA))
  END IF 
  NzrA=0
  DO i=1,A%n
    DO jj=B%RowPtr(1,i),B%RowPtr(2,i)
      j=B%ColInd(jj)
      NzrA=NzrA+1
      A%ColInd(NzrA)=j
    END DO
  END DO
  IF (ASSOCIATED(B%Permu)) THEN
    IF (.NOT.ASSOCIATED(A%Permu)) THEN
      ALLOCATE(A%Permu(A%n))
    END IF
    A%Permu=B%Permu
  END IF
  IF (ASSOCIATED(B%InvPer)) THEN
    IF (.NOT.ASSOCIATED(A%InvPer)) THEN
      ALLOCATE(A%InvPer(A%n))
    END IF
    A%InvPer=B%InvPer
  END IF
END SUBROUTINE RowColD_To_RowColDiag

SUBROUTINE RowCol_To_RowColL(A,B)

!
! Converts
!

  TYPE(SpRowColL), INTENT(INOUT) :: A
  TYPE(SpRowCol), INTENT(IN) :: B

  INTEGER :: i,jj,ind

  A%n=B%n
  A%m=B%m

  IF (.NOT.ASSOCIATED(A%RowPtr)) THEN
    ALLOCATE(A%RowPtr(A%m+1))
    A%RowPtr(1)=1
    DO i=1,B%m
      A%RowPtr(i+1)=A%RowPtr(i)
      DO jj=B%RowPtr(i),B%RowPtr(i+1)-1
        IF (B%ColInd(jj)<i) THEN
          A%RowPtr(i+1)=A%RowPtr(i+1)+1
        END IF
      END DO
    END DO
  END IF
  IF (.NOT.ASSOCIATED(A%ColInd)) THEN
    ALLOCATE(A%ColInd(A%RowPtr(A%n+1)-1))
  END IF
  IF (.NOT.ASSOCIATED(A%Val)) THEN
    ALLOCATE(A%Val(A%RowPtr(A%n+1)-1))
  END IF
  IF (.NOT.ASSOCIATED(A%Diag)) THEN
    ALLOCATE(A%Diag(A%n))
  END IF
  ind=1
  DO i=1,B%m
    DO jj=B%RowPtr(i),B%RowPtr(i+1)-1
      IF (B%ColInd(jj)<i) THEN
        A%ColInd(ind)=B%ColInd(jj)
        A%Val(ind)=B%Val(jj)
        ind=ind+1
      ELSE IF (B%ColInd(jj)==i) THEN
        A%Diag(i)=B%Val(jj)
      END IF
    END DO
  END DO
END SUBROUTINE RowCol_To_RowColL

SUBROUTINE Ind_To_RowCol(A,B)

!
! Converts 
!
  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(INOUT) :: A
  TYPE(SpInd), INTENT(IN) :: B

  INTEGER :: i,jj
  INTEGER :: RowPtr(B%m)

  A%n=B%n
  A%m=B%m

  IF (.NOT.ASSOCIATED(A%RowPtr)) THEN
    ALLOCATE(A%RowPtr(A%m+1))
    A%RowPtr=0
    DO i=1,B%NumNonZero
      A%RowPtr(B%RowInd(i)+1)=A%RowPtr(B%RowInd(i)+1)+1
    END DO
    A%RowPtr(1)=1
    DO i=1,A%m
      A%RowPtr(i+1)=A%RowPtr(i+1)+A%RowPtr(i)
    END DO
  END IF
  IF (.NOT.ASSOCIATED(A%ColInd)) THEN
    ALLOCATE(A%ColInd(B%NumNonZero))
  END IF
  IF (.NOT.ASSOCIATED(A%Val)) THEN
     ALLOCATE(A%Val(B%NumNonZero))
  END IF
  RowPtr(1:A%m)=A%RowPtr(1:A%m)
  DO i=1,B%NumNonZero
     jj=RowPtr(B%RowInd(i))
     RowPtr(B%RowInd(i))=RowPtr(B%RowInd(i))+1
     A%ColInd(jj)=B%ColInd(i)
     A%Val(jj)=B%Val(i)
  END DO
  IF (A%n==A%m) THEN
    IF (.NOT.ASSOCIATED(A%DiagPtr)) THEN
      ALLOCATE(A%DiagPtr(A%n))
    END IF
    DO i=1,A%n
      CALL Sort(A%ColInd(A%RowPtr(i):A%RowPtr(i+1)-1) &
               ,A%Val(A%RowPtr(i):A%RowPtr(i+1)-1))
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        IF (A%ColInd(jj)==i) THEN
          A%DiagPtr(i)=jj
        END IF  
      END DO   
    END DO   
  END IF  
CONTAINS
SUBROUTINE sort(vec,Val)
  INTEGER :: vec(:)
  REAL(RealKind) :: Val(:)

  INTEGER :: i,itemp,j,n
  REAL(RealKind) :: Temp
  n=SIZE(Vec)
  DO i=1,n
    DO j=1,n-i
      IF (vec(j).gt.vec(j+1)) THEN
        itemp=vec(j)
        vec(j)=vec(j+1)
        vec(j+1)=itemp
        temp=Val(j)
        Val(j)=Val(j+1)
        Val(j+1)=Temp
     END IF
    END DO
  END DO
END SUBROUTINE sort

END SUBROUTINE Ind_To_RowCol

SUBROUTINE Ind_To_RowColMat(A,B)

!
! Converts 
!
  IMPLICIT NONE

  TYPE(SpRowColMat), INTENT(INOUT) :: A
  TYPE(SpIndMat), INTENT(IN) :: B

  INTEGER :: i,jj
  INTEGER :: RowPtr(B%m)

  A%n=B%n
  A%m=B%m

  IF (.NOT.ASSOCIATED(A%RowPtr)) THEN
    ALLOCATE(A%RowPtr(A%m+1))
    A%RowPtr=0
    DO i=1,B%NumNonZero
      A%RowPtr(B%RowInd(i)+1)=A%RowPtr(B%RowInd(i)+1)+1
    END DO
    A%RowPtr(1)=1
    DO i=1,A%m
      A%RowPtr(i+1)=A%RowPtr(i+1)+A%RowPtr(i)
    END DO
  END IF
  IF (.NOT.ASSOCIATED(A%ColInd)) THEN
    ALLOCATE(A%ColInd(B%NumNonZero))
  END IF
  IF (.NOT.ASSOCIATED(A%Val)) THEN
     ALLOCATE(A%Val(B%NumNonZero))
  END IF
  RowPtr(1:A%m)=A%RowPtr(1:A%m)
  DO i=1,B%NumNonZero
     jj=RowPtr(B%RowInd(i))
     RowPtr(B%RowInd(i))=RowPtr(B%RowInd(i))+1
     A%ColInd(jj)=B%ColInd(i)
     A%Val(jj)=B%Val(i)
  END DO
  IF (A%n==A%m) THEN
    IF (.NOT.ASSOCIATED(A%DiagPtr)) THEN
      ALLOCATE(A%DiagPtr(A%n))
    END IF
    DO i=1,A%n
      CALL Sort(A%ColInd(A%RowPtr(i):A%RowPtr(i+1)-1) &
               ,A%Val(A%RowPtr(i):A%RowPtr(i+1)-1))
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        IF (A%ColInd(jj)==i) THEN
          A%DiagPtr(i)=jj
        END IF  
      END DO   
    END DO   
  END IF  
CONTAINS
SUBROUTINE sort(vec,Val)
  INTEGER :: vec(:)
  TYPE(A22_T):: Val(:)

  INTEGER :: i,itemp,j,n
  TYPE(A22_T):: Temp
  n=SIZE(Vec)
  DO i=1,n
    DO j=1,n-i
      IF (vec(j).gt.vec(j+1)) THEN
        itemp=vec(j)
        vec(j)=vec(j+1)
        vec(j+1)=itemp
        temp=Val(j)
        Val(j)=Val(j+1)
        Val(j+1)=Temp
     END IF
    END DO
  END DO
END SUBROUTINE sort

END SUBROUTINE Ind_To_RowColMat

SUBROUTINE RowCol_To_Diag(A,B)

!
! Converts 
!
  IMPLICIT NONE

  TYPE(SpDiag), INTENT(INOUT) :: A
  TYPE(SpRowCol), INTENT(IN) :: B

  INTEGER :: i,jj,iNumDiag
  INTEGER :: Ptr(-B%m+1:B%m-1)

  A%n=B%n
  A%m=B%m

  IF (.NOT.ASSOCIATED(A%DiagPtr)) THEN
    Ptr=0
    DO i=1,B%m
      DO jj=B%RowPtr(i),B%RowPtr(i+1)-1
        Ptr(B%ColInd(jj)-i)=Ptr(B%ColInd(jj)-i)+1
      END DO
    END DO
    A%NumDiag=0
    DO i=-A%m+1,A%m-1
      IF (Ptr(i)>0) THEN
        A%NumDiag=A%NumDiag+1
      END IF
    END DO
    ALLOCATE(A%DiagPtr(A%NumDiag))
    iNumDiag=1
    DO i=-A%m+1,A%m-1
      IF (Ptr(i)>0) THEN
        A%DiagPtr(iNumDiag)=i
        Ptr(i)=iNumDiag
        iNumDiag=iNumDiag+1
      END IF
    END DO
  ELSE
    DO i=1,A%NumDiag
      Ptr(A%DiagPtr(i))=i
    END DO
  END IF
  IF (.NOT.ASSOCIATED(A%Val)) THEN
    ALLOCATE(A%Val(A%m,A%NumDiag))
    A%Val=0.0e0
  END IF
  DO i=1,B%m
    DO jj=B%RowPtr(i),B%RowPtr(i+1)-1
      A%Val(B%ColInd(jj),Ptr(B%ColInd(jj)-i))=B%Val(jj)
    END DO
  END DO

END SUBROUTINE RowCol_To_Diag

END MODULE MatSpConvert_Mod

  

