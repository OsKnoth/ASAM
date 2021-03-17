MODULE MatSpRowCol_Mod

  USE Kind_Mod
  USE Index_Mod

  IMPLICIT NONE

  TYPE SpRowCol
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


SUBROUTINE SpTrans_SpRowCol(AT,A)

!
! AT = A**T
!

  TYPE(SpRowCol) :: AT,A


  INTEGER :: iAT,jAT,jjAT
  INTEGER :: iA,jA,jjA
  INTEGER :: iTemp1,iTemp2,iTemp3



  AT%m=A%n
  AT%n=A%m
  ALLOCATE(AT%RowPtr(AT%m+1))
  ALLOCATE(AT%ColInd(A%RowPtr(A%m+1)-1))
  ALLOCATE(AT%Val(A%RowPtr(A%m+1)-1))

! Determine AT%RowPtr structure

  DO iAT=1,AT%m
    AT%RowPtr(iAT)=0
  END DO
  DO jjA=1,A%RowPtr(A%m+1)-1
    jA=A%ColInd(jjA)
    AT%RowPtr(jA)=AT%RowPtr(jA)+1
  END DO
  iTemp1=AT%RowPtr(1)
  iTemp2=AT%RowPtr(2)
  AT%RowPtr(2)=1
  DO iAT=2,AT%m
    iTemp3=AT%RowPtr(iAT+1) 
    AT%RowPtr(iAT+1)=AT%RowPtr(iAT)+iTemp1
    iTemp1=iTemp2
    iTemp2=iTemp3
  END DO

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
      
SUBROUTINE SpMm_SpRowCol(C,A,B)

!
! C=A*B
! 

  TYPE(SpRowCol) :: A,B,C


  INTEGER :: iA,jA,jjA
  INTEGER :: iB,jB,jjB
  INTEGER :: iC,jC,jjC,jjjC
  INTEGER :: NzrIndC(B%n)
  REAL(RealKind) :: Val(B%n)

  WRITE(*,*) 'Inside SpMm_SpRowCol',ASSOCIATED(C%RowPtr)
  IF (.NOT.ASSOCIATED(C%RowPtr)) THEN
    C%m=A%m
    C%n=B%n
    NzrIndC(1:C%n)=0

!   Determine Number of Nonzeroes in C

    ALLOCATE(C%RowPtr(C%m+1))
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
  IF (C%n==C%m) THEN
    IF (.NOT.ASSOCIATED(C%DiagPtr)) THEN
      ALLOCATE(C%DiagPtr(C%n))
    END IF
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
END SUBROUTINE SpMm_SpRowCol 

END MODULE MatSpRowCol_Mod
  

