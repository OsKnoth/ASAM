MODULE MatSpRowCol_Mod

  USE Kind_Mod
  USE Index_Mod

  IMPLICIT NONE

  TYPE :: SpRowCol
    INTEGER :: m,n
    INTEGER, POINTER :: RowPtr(:)=>NULL()
    INTEGER, POINTER :: ColInd(:)=>NULL()
    INTEGER, POINTER :: DiagPtr(:)=>NULL()
    REAL(RealKind), POINTER :: Val(:)=>NULL()
  END TYPE SpRowCol

  REAL(RealKind) :: Shift

CONTAINS

SUBROUTINE AddDiag3_SpRowCol(alpha,x3,A)

  REAL(RealKind) :: alpha,x3(:,:,:)
  TYPE(SpRowCol) :: A

  INTEGER :: i,ii,j,k

  ii=0
  DO k=1,SIZE(x3,3)
    DO j=1,SIZE(x3,2)
      DO i=1,SIZE(x3,1)
        ii=ii+1
        A%Val(A%DiagPtr(ii))=A%Val(A%DiagPtr(ii))+alpha*x3(i,j,k)+Eps
      END DO  
    END DO  
  END DO  
END SUBROUTINE AddDiag3_SpRowCol

SUBROUTINE AddDiag_SpRowCol(alpha,x3,A,Shift)

  REAL(RealKind) :: alpha,x3(:,:,:)
  TYPE(SpRowCol) :: A
  INTEGER :: Shift

  INTEGER :: i,ii,iShift,j,k

  ii=0
  DO k=1,SIZE(x3,3)
    DO j=1,SIZE(x3,2)
      DO i=1,SIZE(x3,1)
        DO iShift=1,Shift
          ii=ii+1
          A%Val(A%DiagPtr(ii))=A%Val(A%DiagPtr(ii))+alpha*x3(i,j,k)+Eps
        END DO    
      END DO  
    END DO  
  END DO  
END SUBROUTINE AddDiag_SpRowCol

SUBROUTINE SpOutput_SpRowCol(A,FileName)

  TYPE(SpRowCol) :: A
  CHARACTER(*) :: FileName

  INTEGER :: i,jj

  OPEN(UNIT=10,FILE=TRIM(FileName)//'Mat',STATUS='UNKNOWN',POSITION='APPEND')
  WRITE(10,*) 'm=',A%m
  WRITE(10,*) 'n=',A%n
  DO i=1,A%m
    IF (ASSOCIATED(A%DiagPtr)) THEN
      WRITE(10,*) 'Diag',A%DiagPtr(i)
    END IF  
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
       WRITE(10,*) i,A%ColInd(jj),A%Val(jj),jj
    END DO
  END DO
  CLOSE(10)
END SUBROUTINE SpOutput_SpRowCol

SUBROUTINE SpNullify_SpRowCol(A)

!
! Nullify Pointers of A 
!

  TYPE(SpRowCol) :: A

  NULLIFY(A%RowPtr)
  NULLIFY(A%ColInd)
  NULLIFY(A%DiagPtr)
  NULLIFY(A%Val)

END SUBROUTINE SpNullify_SpRowCol

SUBROUTINE SpDeallocate_SpRowCol(A)

!
! Deallocate A
!

  TYPE(SpRowCol) :: A

  A%m=0
  A%n=0
  IF (ASSOCIATED(A%RowPtr)) THEN
    DEALLOCATE(A%RowPtr)
  END IF
  IF (ASSOCIATED(A%ColInd)) THEN
    DEALLOCATE(A%ColInd)
  END IF
  IF (ASSOCIATED(A%DiagPtr)) THEN
    DEALLOCATE(A%DiagPtr)
  END IF
  IF (ASSOCIATED(A%Val)) THEN
    DEALLOCATE(A%Val)
  END IF

END SUBROUTINE SpDeallocate_SpRowCol

SUBROUTINE SpAVecI3_SpRowCol(y3,A,x3)

!
! Computes y=A*x
!

  TYPE(SpRowCol) :: A
  REAL(RealKind) :: y3(:,:,:)
  REAL(RealKind) :: x3(:,:,:)

  INTEGER :: i,j,k,jj,n1,n2,n3
  REAL(RealKind) :: y(SIZE(y3))
  REAL(RealKind) :: x(SIZE(x3))

  n1=SIZE(x3,1)
  n2=SIZE(x3,2)
  n3=SIZE(x3,3)
  DO k=1,n3
    DO j=1,n2
      DO i=1,n1
        x(IndexC(i,j,k,n1,n2,n3))=x3(i,j,k)
      END DO  
    END DO  
  END DO  
  DO i=1,A%m
    y(i)=0.0e0
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)+A%Val(jj)*x(A%ColInd(jj))
    END DO
  END DO
  IF (Shift>Zero) THEN
    DO i=1,A%m
      y(i)=y(i)+Shift*x(i)
    END DO  
  END IF
  DO k=1,n3
    DO j=1,n2
      DO i=1,n1
        y3(i,j,k)=y(IndexC(i,j,k,n1,n2,n3))
      END DO  
    END DO  
  END DO  
END SUBROUTINE SpAVecI3_SpRowCol

SUBROUTINE SpAVecI4_SpRowCol(y4,A,x4)

!
! Computes y=A*x
!

  TYPE(SpRowCol) :: A
  REAL(RealKind) :: y4(:,:,:,:)
  REAL(RealKind) :: x4(:,:,:,:)

  INTEGER :: i0,i,j,k,jj,n0,n1,n2,n3
  REAL(RealKind) :: y(SIZE(y4))
  REAL(RealKind) :: x(SIZE(x4))

  n0=SIZE(x4,1)
  n1=SIZE(x4,2)
  n2=SIZE(x4,3)
  n3=SIZE(x4,4)
  DO k=1,n3
    DO j=1,n2
      DO i=1,n1
        DO i0=1,n0
          x(IndexC(i0,i,j,k,n0,n1,n2,n3))=x4(i0,i,j,k)
        END DO  
      END DO  
    END DO  
  END DO  
  DO i=1,A%m
    y(i)=0.0e0
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)+A%Val(jj)*x(A%ColInd(jj))
    END DO
  END DO
  DO k=1,n3
    DO j=1,n2
      DO i=1,n1
        DO i0=1,n0
          y4(i0,i,j,k)=y(IndexC(i0,i,j,k,n0,n1,n2,n3))
        END DO  
      END DO  
    END DO  
  END DO  
END SUBROUTINE SpAVecI4_SpRowCol

SUBROUTINE GaussSeidel_SpRowCol(x,A,b,NumIter)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A = L + D + U
!


  TYPE(SpRowCol) :: A
  REAL(RealKind) :: x(:)
  REAL(RealKind) :: b(:)
  INTEGER :: NumIter

  INTEGER :: i,j,k,jj,Iter
  REAL(RealKind) :: y(A%n)

  DO Iter=1,NumIter
    DO i=1,A%n
      y(i)=b(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
      END DO
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(i),A%DiagPtr(i)-1
        y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
      END DO
      y(i)=y(i)/(A%Val(A%DiagPtr(i))+Shift)
    END DO
    DO i=1,A%n 
      y(i)=y(i)*(A%Val(A%DiagPtr(i))+Shift)
    END DO
    DO i=A%n,1,-1
      DO jj=A%DiagPtr(i)+1,A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
      END DO
      y(i)=y(i)/(A%Val(A%DiagPtr(i))+Shift)
    END DO
    x(1:A%n)=x(1:A%n)+y
  END DO
      
END SUBROUTINE GaussSeidel_SpRowCol


SUBROUTINE Overrelaxation_lexico(x, A, b, NumIter, omega)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(IN) :: A
  REAL(Realkind), INTENT(IN) :: b(:)
  REAL(Realkind), INTENT(INOUT) :: x(:)
  REAL, INTENT(IN) :: omega
  INTEGER, INTENT(IN) :: NumIter
  
  INTEGER :: i, j, iter
  REAL(Realkind) :: ax


  DO iter = 1, NumIter
    DO i = 1, A%n
      ax = 0.0
      DO j = A%Rowptr(i), A%Rowptr(i + 1) - 1
        ax = ax + A%Val(j) * x(A%Colind(j))
      END DO
      x(i) = x(i) + (b(i) - ax) / A%Val(A%Diagptr(i)) * omega
    END DO
  END DO

END SUBROUTINE

SUBROUTINE Overrelaxation_redblack(x, A, b, rd_ind, bl_ind, NumIter, omega)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(IN) :: A
  REAL(Realkind), INTENT(IN) :: b(:)
  REAL(Realkind), INTENT(INOUT) :: x(:)
  INTEGER, INTENT(IN) :: rd_ind(:), bl_ind(:)
  REAL, INTENT(IN) :: omega
  INTEGER, INTENT(IN) :: NumIter

  INTEGER :: i, j, ii, iter
  REAL(Realkind) :: ax
  
  DO iter = 1, NumIter
    DO ii = 1, SIZE(rd_ind)
      i = rd_ind(ii)
      ax = 0.0
      DO j = A%Rowptr(i), A%Rowptr(i + 1) - 1
        ax = ax + A%Val(j) * x(A%Colind(j))
      END DO
      x(i) = x(i) + (b(i) - ax) / A%Val(A%Diagptr(i)) * omega
    END DO
    DO ii = 1, SIZE(bl_ind)
      i = bl_ind(ii)
      ax = 0.0
      DO j = A%Rowptr(i), A%Rowptr(i + 1) - 1
        ax = ax + A%Val(j) * x(A%Colind(j))
      END DO
      x(i) = x(i) + (b(i) - ax) / A%Val(A%Diagptr(i)) * omega
    END DO
  END DO

END SUBROUTINE 

SUBROUTINE GaussSeidelI3_SpRowCol(x3,A,b3,NumIter)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A = L + D + U
!


  TYPE(SpRowCol) :: A
  REAL(RealKind) :: x3(:,:,:)
  REAL(RealKind) :: b3(:,:,:)
  INTEGER :: NumIter

  INTEGER :: i,j,k,jj,n1,n2,n3,Iter
  REAL(RealKind) :: y(SIZE(x3))
  REAL(RealKind) :: x(SIZE(x3))
  REAL(RealKind) :: b(SIZE(x3))

  n1=SIZE(x3,1)
  n2=SIZE(x3,2)
  n3=SIZE(x3,3)
  DO k=1,n3
    DO j=1,n2
      DO i=1,n1
        x(IndexC(i,j,k,n1,n2,n3))=x3(i,j,k)
        b(IndexC(i,j,k,n1,n2,n3))=b3(i,j,k)
      END DO  
    END DO  
  END DO  
  DO Iter=1,NumIter
    DO i=1,A%n
      y(i)=b(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
      END DO
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(i),A%DiagPtr(i)-1
        y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
      END DO
      y(i)=y(i)/(A%Val(A%DiagPtr(i))+Shift)
    END DO
    DO i=1,A%n 
      y(i)=y(i)*(A%Val(A%DiagPtr(i))+Shift)
    END DO
    DO i=A%n,1,-1
      DO jj=A%DiagPtr(i)+1,A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
      END DO
      y(i)=y(i)/(A%Val(A%DiagPtr(i))+Shift)
    END DO
    x=x+y
  END DO
      
  DO k=1,n3
    DO j=1,n2
      DO i=1,n1
        x3(i,j,k)=x(IndexC(i,j,k,n1,n2,n3))
      END DO  
    END DO  
  END DO  
END SUBROUTINE GaussSeidelI3_SpRowCol

SUBROUTINE GaussSeidelI4_SpRowCol(x4,A,b4,NumIter)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A = L + D + U
!


  TYPE(SpRowCol) :: A
  REAL(RealKind) :: x4(:,:,:,:)
  REAL(RealKind) :: b4(:,:,:,:)
  INTEGER :: NumIter

  INTEGER :: i0,i,j,k,jj,n0,n1,n2,n3,Iter
  REAL(RealKind) :: y(SIZE(x4))
  REAL(RealKind) :: x(SIZE(x4))
  REAL(RealKind) :: b(SIZE(x4))

  n0=SIZE(x4,1)
  n1=SIZE(x4,2)
  n2=SIZE(x4,3)
  n3=SIZE(x4,4)
  DO k=1,n3
    DO j=1,n2
      DO i=1,n1
        DO i0=1,n0
          x(IndexC(i0,i,j,k,n0,n1,n2,n3))=x4(i0,i,j,k)
          b(IndexC(i0,i,j,k,n0,n1,n2,n3))=b4(i0,i,j,k)
        END DO  
      END DO  
    END DO  
  END DO  
  DO Iter=1,NumIter
    DO i=1,A%n
      y(i)=b(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
      END DO
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(i),A%DiagPtr(i)-1
        y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
      END DO
      y(i)=y(i)/A%Val(A%DiagPtr(i))
    END DO
    DO i=1,A%n
      y(i)=y(i)*A%Val(A%DiagPtr(i))
    END DO
    DO i=A%n,1,-1
      DO jj=A%DiagPtr(i)+1,A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
      END DO
      y(i)=y(i)/A%Val(A%DiagPtr(i))
    END DO
    x=x+y
  END DO
      
  DO k=1,n3
    DO j=1,n2
      DO i=1,n1
        DO i0=1,n0
          x4(i0,i,j,k)=x(IndexC(i0,i,j,k,n0,n1,n2,n3))
        END DO  
      END DO  
    END DO  
  END DO  
END SUBROUTINE GaussSeidelI4_SpRowCol

SUBROUTINE SpAVec_SpRowCol(y,A,x)

!
! Computes y=A*x
!

  TYPE(SpRowCol) :: A
  REAL(RealKind) :: y(:)
  REAL(RealKind) :: x(:)

  INTEGER :: i,jj

  DO i=1,A%m
    y(i)=0.0e0
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)+A%Val(jj)*x(A%ColInd(jj))
    END DO
  END DO
END SUBROUTINE SpAVec_SpRowCol

SUBROUTINE Axpy_SpRowCol(A,x,y)

!
! Computes y=A*x+y
!

  TYPE(SpRowCol) :: A
  REAL(RealKind) :: x(:)
  REAL(RealKind) :: y(:)

  INTEGER :: i,jj

  DO i=1,A%m
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)+A%Val(jj)*x(A%ColInd(jj))
    END DO
  END DO
END SUBROUTINE Axpy_SpRowCol

SUBROUTINE ATxpy_SpRowCol(A,x,y)

!
! Computes y=(A**T)*x+y
!

  TYPE(SpRowCol) :: A
  REAL(RealKind) :: x(:)
  REAL(RealKind) :: y(:)

  INTEGER :: i,jj

  DO i=1,A%m
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(A%ColInd(jj))=y(A%ColInd(jj))+A%Val(jj)*x(i)
    END DO
  END DO
END SUBROUTINE ATxpy_SpRowCol


SUBROUTINE SpTrans_SpRowCol(AT,A,allocate_mat)

!
! AT = A**T
!
  TYPE(SpRowCol), INTENT(inout) :: AT 
  TYPE(SpRowCol), INTENT(in) :: A
  LOGICAL, OPTIONAL, INTENT(in) :: allocate_mat

  LOGICAL :: allocate_mat_def

  INTEGER :: iAT,jAT,jjAT
  INTEGER :: iA,jA,jjA
  INTEGER :: iTemp1,iTemp2,iTemp3


  allocate_mat_def = .TRUE.

  IF (PRESENT(allocate_mat)) allocate_mat_def = allocate_mat

  AT%m=A%n
  AT%n=A%m

  IF (allocate_mat_def) THEN
    ALLOCATE(AT%RowPtr(AT%m+1))
    ALLOCATE(AT%ColInd(A%RowPtr(A%m+1)-1))
    ALLOCATE(AT%Val(A%RowPtr(A%m+1)-1))
  END IF

! Determine AT%RowPtr structure

  DO iAT=1,AT%m
    AT%RowPtr(iAT)=0
  END DO
  DO jjA=1,A%RowPtr(A%m+1)-1
    jA=A%ColInd(jjA)
    AT%RowPtr(jA)=AT%RowPtr(jA)+1
  END DO
  iTemp1=AT%RowPtr(1)

  IF(SIZE(AT%RowPtr) .GT. 1) THEN
    iTemp2=AT%RowPtr(2)
    AT%RowPtr(2)=1
    DO iAT=2,AT%m
      iTemp3=AT%RowPtr(iAT+1) 
      AT%RowPtr(iAT+1)=AT%RowPtr(iAT)+iTemp1
      iTemp1=iTemp2
      iTemp2=iTemp3
    END DO
  END IF

  DO iA=1,A%m
    DO jjA=A%RowPtr(iA),A%RowPtr(iA+1)-1
      jA=A%ColInd(jjA)
      jjAT=AT%RowPtr(jA+1)
      AT%Colind(jjAT)=iA
      AT%Val(jjAT)=A%Val(jjA)
      AT%RowPtr(jA+1)=jjAT+1
    END DO
  END DO
  AT%RowPtr(1)=1

END SUBROUTINE SpTrans_SpRowCol
      
SUBROUTINE SpMm_SpRowCol(C,A,B, restruct_mat, withdiagptr)

!
! C=A*B
! 

  TYPE(SpRowCol), INTENT(in) :: A,B
  TYPE(SpRowCol), INTENT(inout) :: C

  INTEGER :: iA,jA,jjA
  INTEGER :: iB,jB,jjB
  INTEGER :: iC,jC,jjC,jjjC
  INTEGER :: NzrIndC(B%n)
  INTEGER :: Status
  INTEGER :: n_old
  REAL(RealKind) :: Val(B%n)
  LOGICAL, OPTIONAL, INTENT(in) :: withdiagptr, &
                                   restruct_mat
  LOGICAL :: withdiagptr_def, &
             restruct_mat_def 

  LOGICAL :: comp_nnz

  withdiagptr_def = .True.
  IF (PRESENT(withdiagptr)) withdiagptr_def = withdiagptr

  restruct_mat_def = .True.
  IF (PRESENT(restruct_mat)) restruct_mat_def = restruct_mat

  comp_nnz = .TRUE.

  n_old = C%n
  C%m=A%m
  C%n=B%n 
 
  IF (.NOT.ASSOCIATED(C%RowPtr)) THEN
    ALLOCATE(C%RowPtr(A%m+1))
    restruct_mat_def = .True.      
    NzrIndC(1:C%n)=0

!   Determine Number of Nonzeroes in C

    jjC=1
    DO iA=1,A%m
      C%RowPtr(iA)=jjC
      DO jjA=A%RowPtr(iA),A%RowPtr(iA+1)-1
        jA=A%ColInd(jjA)
        iB=jA
        DO jjB=B%RowPtr(iB),B%RowPtr(iB+1)-1
          jB=B%ColInd(jjB)
          IF (NzrIndC(jB)/=iA) THEN
            jjC=jjC+1
            NzrIndC(jB)=iA
          END IF
        END DO
      END DO
    END DO
    C%RowPtr(C%m+1)=jjC
    ALLOCATE(C%ColInd(jjC-1))
    ALLOCATE(C%Val(jjC-1))
    comp_nnz = .FALSE.
  END IF

  IF (restruct_mat_def) THEN
    IF (comp_nnz) THEN
    !   Determine Number of Nonzeroes in C
      NzrIndC(1:C%n)=0
      jjC=1
      DO iA=1,A%m
        C%RowPtr(iA)=jjC
        DO jjA=A%RowPtr(iA),A%RowPtr(iA+1)-1
          jA=A%ColInd(jjA)
          iB=jA
          DO jjB=B%RowPtr(iB),B%RowPtr(iB+1)-1
            jB=B%ColInd(jjB)
            IF (NzrIndC(jB)/=iA) THEN
              jjC=jjC+1
              NzrIndC(jB)=iA
            END IF
          END DO
        END DO
      END DO
      C%RowPtr(C%m+1)=jjC
    END IF

    NzrIndC(1:C%n)=0 
    jjC=1
    DO iA=1,A%m
      DO jjA=A%RowPtr(iA),A%RowPtr(iA+1)-1
        jA=A%ColInd(jjA)
        iB=jA
        DO jjB=B%RowPtr(iB),B%RowPtr(iB+1)-1
          jB=B%ColInd(jjB)
          IF (NzrIndC(jB)/=iA) THEN
            C%ColInd(jjC)=jB
            jjC=jjC+1
            NzrIndC(jB)=iA
            Val(jB)=A%Val(jjA)*B%Val(jjB)
          ELSE
            Val(jB)=Val(jB)+A%Val(jjA)*B%Val(jjB)
          END IF
        END DO
      END DO
      DO jjjC=C%RowPtr(iA),C%RowPtr(iA+1)-1
        jC=C%ColInd(jjjC)
        C%Val(jjjC)=Val(jC)
      END DO
    END DO
  ELSE
    jjC=1
    Val=0.0e0
    DO iA=1,A%m
      DO jjA=A%RowPtr(iA),A%RowPtr(iA+1)-1
        jA=A%ColInd(jjA)
        iB=jA
        DO jjB=B%RowPtr(iB),B%RowPtr(iB+1)-1
          jB=B%ColInd(jjB)
          Val(jB)=Val(jB)+A%Val(jjA)*B%Val(jjB)
        END DO
      END DO
      DO jjjC=C%RowPtr(iA),C%RowPtr(iA+1)-1
        jC=C%ColInd(jjjC)
        C%Val(jjjC)=Val(jC)
        Val(jC)=0.0e0
      END DO
    END DO
  END IF
  IF ((C%n==C%m) .AND. withdiagptr_def .EQV. .True.) THEN
    IF (.NOT.ASSOCIATED(C%DiagPtr)) THEN
      ALLOCATE(C%DiagPtr(C%n))
    END IF
    C%DiagPtr(:)=-1
    DO iC=1,C%n
      CALL Sort(C%ColInd(C%RowPtr(iC):C%RowPtr(iC+1)-1) &
               ,C%Val(C%RowPtr(iC):C%RowPtr(iC+1)-1))
      DO jjC=C%RowPtr(iC),C%RowPtr(iC+1)-1
        IF (C%ColInd(jjC)==iC) THEN
          C%DiagPtr(iC)=jjC
        END IF  
      END DO   
    END DO   
  END IF  

END SUBROUTINE SpMm_SpRowCol

SUBROUTINE SpA_SpRowCol(C,A,B, allocate_mat, restruct_mat, withdiagptr)
!
! C=A+B
! 

  TYPE(SpRowCol), INTENT(in) :: A,B
  TYPE(SpRowCol), INTENT(inout) :: C
  LOGICAL, OPTIONAL, INTENT(in) :: allocate_mat, restruct_mat, withdiagptr

  INTEGER :: iA,jA,jjA
  INTEGER :: iB,jB,jjB
  INTEGER :: iC,jC,jjC,jjjC
  INTEGER :: NzrIndC(A%n)
  REAL(Realkind) :: Vals_tmp(A%n)
  LOGICAL :: isin
  LOGICAL :: allocate_mat_def, restruct_mat_def, withdiagptr_def

  allocate_mat_def = .TRUE.
  restruct_mat_def = .TRUE.
  withdiagptr_def = .TRUE.

  IF (PRESENT(allocate_mat)) allocate_mat_def = allocate_mat
  IF (PRESENT(restruct_mat)) restruct_mat_def = restruct_mat
  IF (PRESENT(withdiagptr)) withdiagptr_def = withdiagptr


  C%m=A%m
  C%n=A%n

  IF (allocate_mat_def) THEN

    restruct_mat_def = .True.

    IF (.NOT.ASSOCIATED(C%RowPtr)) THEN
!   Determine Number of Nonzeroes in C

      NzrIndC(1:C%n)=0

      ALLOCATE(C%RowPtr(C%m+1))
      jjC=1
      DO iA=1,A%m
        C%RowPtr(iA)=jjC
        DO jjA=A%RowPtr(iA),A%RowPtr(iA+1)-1
          jA=A%ColInd(jjA)
          IF (NzrIndC(jA)/=iA) THEN
            jjC=jjC+1
            NzrIndC(jA)=iA
          END IF
        END DO
        DO jjB=B%RowPtr(iA),B%RowPtr(iA+1)-1
          jB=B%ColInd(jjB)
          IF (NzrIndC(jB)/=iA) THEN
            jjC=jjC+1
            NzrIndC(jB)=iA
          END IF
        END DO
      END DO

      C%RowPtr(C%m+1)=jjC
      ALLOCATE(C%ColInd(jjC-1))
      ALLOCATE(C%Val(jjC-1))

    ELSE
      WRITE(*,*) "Cannot allocate matrix. RowPtr is already associated"
      STOP 'error'
    END IF
  END IF

  IF (restruct_mat_def) THEN
    IF (.NOT. allocate_mat_def) THEN
      NzrIndC(1:C%n)=0
      jjC=1
      DO iA=1,A%m
        C%RowPtr(iA)=jjC
        DO jjA=A%RowPtr(iA),A%RowPtr(iA+1)-1
          jA=A%ColInd(jjA)
          IF (NzrIndC(jA)/=iA) THEN
            jjC=jjC+1
            NzrIndC(jA)=iA
          END IF
        END DO
        DO jjB=B%RowPtr(iA),B%RowPtr(iA+1)-1
          jB=B%ColInd(jjB)
          IF (NzrIndC(jB)/=iA) THEN
            jjC=jjC+1
            NzrIndC(jB)=iA
          END IF
        END DO
      END DO
      C%RowPtr(C%m+1)=jjC
    END IF
  END IF

  DO iA=1,A%m
    jjC = 0
    DO jjA=A%RowPtr(iA),A%RowPtr(iA+1)-1
      jjC=jjC+1
      jA=A%ColInd(jjA)
      Vals_tmp(jA)=A%Val(jjA)
      NzrIndC(jjC)=jA
    END DO

    DO jjB=B%RowPtr(iA),B%RowPtr(iA+1)-1
      jB=B%ColInd(jjB)

      isin = .FALSE.
      DO jjjC = 1, jjC 
        IF (NzrIndC(jjjC) .EQ. jB) THEN
          Vals_tmp(jB) = Vals_tmp(jB) + B%Val(jjB)
          isin = .TRUE.
        END IF
      END DO
      IF (isin .EQV. .FALSE.) THEN
        jjC=jjC+1
        Vals_tmp(jB)=B%Val(jjB)
        NzrIndC(jjC)=jB
      END IF
    END DO

    DO jjc=C%RowPtr(iA),C%RowPtr(iA+1)-1
      C%ColInd(jjc)=NzrIndC(jjc-C%RowPtr(iA)+1)
      C%Val(jjc)=Vals_tmp(C%ColInd(jjc))
    END DO
    CALL Sort(C%ColInd(C%RowPtr(iA):C%RowPtr(iA+1)-1) &
             ,C%Val(C%RowPtr(iA):C%RowPtr(iA+1)-1))
  END DO

  IF ((C%n==C%m) .AND. withdiagptr_def) THEN
    IF (.NOT.ASSOCIATED(C%DiagPtr)) THEN
      ALLOCATE(C%DiagPtr(C%n))
    END IF
    C%DiagPtr(:)=-1
    DO iC=1,C%n
      CALL Sort(C%ColInd(C%RowPtr(iC):C%RowPtr(iC+1)-1) &
               ,C%Val(C%RowPtr(iC):C%RowPtr(iC+1)-1))
      DO jjC=C%RowPtr(iC),C%RowPtr(iC+1)-1
        IF (C%ColInd(jjC)==iC) THEN
          C%DiagPtr(iC)=jjC
        END IF
      END DO
    END DO
  END IF

END SUBROUTINE SpA_SpRowCol

SUBROUTINE SpAScalRow_SpRowCol(A, scal)
  !
  ! computes A(i,:) = scal(i) * A(i,:) 
  !
  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: A
  REAL(Realkind), INTENT(in) :: scal(:)
  INTEGER :: irow, ival

  DO irow = 1, A%m
    DO ival = A%RowPtr(irow), A%RowPtr(irow + 1) - 1 
      A%Val(ival) = A%Val(ival) * scal(irow)
    END DO 
  END DO

END SUBROUTINE SpAScalRow_SpRowCol


SUBROUTINE SpAScalinvRow_SpRowCol(A, scal)
  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: A
  REAL(Realkind), INTENT(in) :: scal(:)
  INTEGER :: irow, ival
  REAl(Realkind), PARAMETER :: zero = 0.0
  REAL(Realkind), PARAMETER :: eps = 1e-20

  DO irow = 1, A%m
    DO ival = A%RowPtr(irow), A%RowPtr(irow + 1) - 1
!      IF (scal(irow) .NE. zero) THEN
         A%Val(ival) = A%Val(ival) / (scal(irow) + eps)
!      ELSE
!         A%Val(ival) = zero
!      END IF
    END DO
  END DO
END SUBROUTINE SpAScalinvRow_SpRowCol

SUBROUTINE SpAScalCol_SpRowCol(A, scal)
  !
  ! computes A(:,i) = A(:,i) * scal(i)
  !
  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: A
  REAL(Realkind), INTENT(in) :: scal(:)

  INTEGER :: irow, ival, icol

  DO irow = 1, A%m
    DO ival = A%RowPtr(irow), A%RowPtr(irow + 1) - 1
      icol = A%ColInd(ival)
      A%Val(ival) = A%Val(ival) * scal(icol)
    END DO
  END DO

END SUBROUTINE SpAScalCol_SpRowCol


SUBROUTINE SpMmBAC(RES, B, A, C, nnz)
  IMPLICIT NONE
  !
  ! Computes RES = B * A * C 
  !
  TYPE(SpRowCol), INTENT(inout) :: RES
  TYPE(SpRowCol), INTENT(in) :: A, B, C
  INTEGER, INTENT(in) :: nnz

  TYPE(SpRowCol) :: AC
  INTEGER, TARGET :: AC_RowPtr(A%m + 1)
  INTEGER, TARGET :: AC_ColInd(nnz)
  REAL(Realkind), TARGET :: AC_val(nnz)

  AC%RowPtr => AC_RowPtr
  AC%ColInd => AC_ColInd
  AC%Val => AC_val

  CALL SpMm_SpRowCol(AC, A, C, restruct_mat=.TRUE., withdiagptr=.FALSE.)

!  WRITE(*,*) "B", B%m, SIZE(B%RowPtr, 1) - 1, B%RowPtr(B%m + 1) - 1, SIZE(B%ColInd, 1), &
!             B%n, MINVAL(B%ColInd, 1), MAXVAL(B%ColInd, 1)
!  WRITE(*,*) "AC", AC%m, SIZE(AC%RowPtr, 1) - 1, AC%RowPtr(AC%m + 1) - 1, SIZE(AC%ColInd, 1), &
!             AC%n, MINVAL(AC%ColInd, 1), MAXVAL(AC%ColInd, 1)
  CALL SpMm_SpRowCol(RES, B, AC)

  CALL SpNullify_SpRowCol(AC)

END SUBROUTINE SpMmBAC


SUBROUTINE SpMmBAC_tmp(RES, B, A, C, nnz)
  IMPLICIT NONE
  !
  ! Computes RES = B * A * C
  !
  TYPE(SpRowCol), INTENT(inout) :: RES
  TYPE(SpRowCol), INTENT(in) :: A, B, C
  INTEGER, INTENT(in) :: nnz

  TYPE(SpRowCol) :: AC
  INTEGER, TARGET :: AC_RowPtr(A%m + 1)
  INTEGER, TARGET :: AC_ColInd(nnz)
  REAL(Realkind), TARGET :: AC_val(nnz)

  AC%RowPtr => AC_RowPtr
  AC%ColInd => AC_ColInd
  AC%Val => AC_val

  CALL SpMm_SpRowCol(AC, A, C, restruct_mat=.TRUE., withdiagptr=.FALSE.)

  WRITE(*,*) "B", B%m, SIZE(B%RowPtr, 1) - 1, B%RowPtr(B%m + 1) - 1, SIZE(B%ColInd, 1), &
             B%n, MINVAL(B%ColInd, 1), MAXVAL(B%ColInd, 1)
  WRITE(*,*) "AC", AC%m, SIZE(AC%RowPtr, 1) - 1, AC%RowPtr(AC%m + 1) - 1, SIZE(AC%ColInd, 1), &
             AC%n, MINVAL(AC%ColInd, 1), MAXVAL(AC%ColInd, 1)
  CALL SpMm_SpRowCol(RES, B, AC)

  WRITE(*,*) "RES", RES%m, SIZE(RES%RowPtr, 1) - 1, RES%RowPtr(RES%m + 1) - 1, SIZE(RES%ColInd, 1), &
             RES%n, MINVAL(RES%ColInd, 1), MAXVAL(RES%ColInd, 1)
  
  CALL SpNullify_SpRowCol(AC)

END SUBROUTINE SpMmBAC_tmp


SUBROUTINE SpMmBACT(RES, B, A, C)
  IMPLICIT NONE
  !
  ! Computes RES = B * A * C^T
  !

  TYPE(SpRowCol), INTENT(inout) :: RES
  TYPE(SpRowCol), INTENT(in) :: A, B, C

  TYPE(SpRowCol) :: CT
  INTEGER :: nnz, iA, jjA, jA, jC, jjC

  INTEGER, TARGET :: CT_RowPtr(C%n + 1), &
                     CT_ColInd(SIZE(C%ColInd))
  REAL(Realkind), TARGET :: CT_Val(SIZE(C%Val))

  CT%RowPtr => CT_RowPtr
  CT%ColInd => CT_ColInd
  CT%Val => CT_Val

  CALL SpTrans_SpRowCol(CT, C, allocate_mat=.FALSE.)

  ! Compute nnz of AC^T (To keep fast automatic memory allocation from stack in the next step)
  nnz = SpMm_nnz(A, CT)
  CALL SpMmBAC(RES, B, A, CT, nnz)
!  WRITE(*,*) "done"
  CALL SpNullify_SpRowCol(CT)

END SUBROUTINE 


SUBROUTINE make_eye(eye_op, cellisopen)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(INOUT) :: eye_op
  REAL(Realkind), INTENT(in) :: cellisopen(:)
  INTEGER :: i, n

  n = SIZE(cellisopen, 1)

  eye_op%n = n
  eye_op%m = n

  ALLOCATE(eye_op%RowPtr(n + 1))
  ALLOCATE(eye_op%ColInd(n))
  ALLOCATE(eye_op%Val(n))

  eye_op%Val(:) = cellisopen

  eye_op%RowPtr(1) = 1

  DO i = 1, n
    eye_op%ColInd(i) = i
    eye_op%RowPtr(i + 1) = i + 1
  END DO

END SUBROUTINE make_eye


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

SUBROUTINE resize_sprowcol_dataflds(spcr, nmax)
  TYPE(SpRowCol), INTENT(INOUT) :: spcr
  INTEGER, POINTER :: int_new(:)
  REAL(Realkind), POINTER :: real_new(:)
  INTEGER, INTENT(IN) :: nmax

  ALLOCATE(int_new(nmax))
  ALLOCATE(real_new(nmax))

  int_new(:) = spcr%ColInd(1:nmax)
  real_new(:) = spcr%val(1:nmax)
  DEALLOCATE(spcr%ColInd)
  DEALLOCATE(spcr%val)
  spcr%ColInd => int_new
  spcr%val => real_new

END SUBROUTINE resize_sprowcol_dataflds


FUNCTION SpMm_nnz(A, B) RESULT(nnz)
  IMPLICIT NONE
  TYPE(SpRowCol), INTENT(in) :: A, B

  INTEGER :: NzrIndC(B%n)
  INTEGER :: nnz, iA, jA, jjA, jjB, jB, iB
  
  nnz = 0
  NzrIndC(1:B%n) = 0

  nnz = 0
  DO iA = 1, A%m
    DO jjA = A%RowPtr(iA), A%RowPtr(iA + 1) - 1
      jA = A%ColInd(jjA)
      iB = jA
      DO jjB = B%RowPtr(iB), B%RowPtr(iB + 1) - 1
        jB = B%ColInd(jjB)
        IF (NzrIndC(jB) /= iA) THEN
          nnz = nnz + 1
          NzrIndC(jB) = iA
        END IF
      END DO
    END DO
  END DO

END FUNCTION SpMm_nnz




END MODULE MatSpRowCol_Mod
  

