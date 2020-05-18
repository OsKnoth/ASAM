MODULE MatSpRowColL_Mod

  USE Kind_Mod

  USE Index_Mod

  IMPLICIT NONE

  TYPE SpRowColL
    INTEGER :: m,n
    INTEGER, POINTER :: RowPtr(:)=>NULL()
    INTEGER, POINTER :: ColInd(:)=>NULL()
    REAL(RealKind), POINTER :: Val(:)=>NULL()
    REAL(RealKind), POINTER :: Diag(:)=>NULL()
    REAL(RealKind), POINTER :: PDiag(:)=>NULL()
  END TYPE SpRowColL
  INTEGER, PRIVATE :: nx,ny,nz
  INTEGER, PRIVATE :: ix,iy,iz

CONTAINS

SUBROUTINE SpNullify_SpRowColL(A)

!
! Nullify Pointers of A 
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A

  NULLIFY(A%RowPtr)
  NULLIFY(A%ColInd)
  NULLIFY(A%Val)
  NULLIFY(A%Diag)
  NULLIFY(A%PDiag)

END SUBROUTINE SpNullify_SpRowColL

SUBROUTINE SpDeallocate_SpRowColL(A)

!
! Deallocate A
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A

  A%m=0
  A%n=0
  IF (ASSOCIATED(A%RowPtr)) THEN
    DEALLOCATE(A%RowPtr)
  END IF
  IF (ASSOCIATED(A%ColInd)) THEN
    DEALLOCATE(A%ColInd)
  END IF
  IF (ASSOCIATED(A%Val)) THEN
    DEALLOCATE(A%Val)
  END IF
  IF (ASSOCIATED(A%Diag)) THEN
    DEALLOCATE(A%Diag)
  END IF
  IF (ASSOCIATED(A%PDiag)) THEN
    DEALLOCATE(A%PDiag)
  END IF

END SUBROUTINE SpDeallocate_SpRowColL

SUBROUTINE IC_SpRowColL(A)

!
! Computes Incomplete Cholesky
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A

  INTEGER :: i,jj

  IF (.NOT.ASSOCIATED(A%PDiag)) THEN
    ALLOCATE(A%PDiag(A%n))
  END IF
  DO i=1,A%n
    A%PDiag(i)=A%Diag(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      A%PDiag(i)=A%PDiag(i)-A%Val(jj)**2/A%PDiag(A%ColInd(jj))
    END DO
    A%PDiag(i)=ABS(A%PDiag(i))+Eps
  END DO

END SUBROUTINE IC_SpRowColL

SUBROUTINE ICSolveMI3_SpRowColL(x3,A,b3,NumIter)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A = L + D + U
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: x3(:,:,:)
  REAL(RealKind) :: b3(:,:,:)
  INTEGER :: NumIter

  INTEGER :: i,jj,Iter
  REAL(RealKind) :: y(A%n)
  REAL(RealKind) :: x(SIZE(x3))
  REAL(RealKind) :: b(SIZE(b3))

  nx=SIZE(x3,1)
  ny=SIZE(x3,2)
  nz=SIZE(x3,3)
  DO iz=1,nz
    DO iy=1,ny
      DO ix=1,nx
        b(IndexC(ix,iy,iz,nx,ny,nz))=b3(ix,iy,iz)
        x(IndexC(ix,iy,iz,nx,ny,nz))=x3(ix,iy,iz)
      END DO
    END DO
  END DO

  DO Iter=1,NumIter
    DO i=1,A%n
      y(i)=b(i)-A%Diag(i)*x(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
        y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*x(i)
      END DO
    END DO

    DO i=1,A%n
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
      END DO
      y(i)=y(i)/A%PDiag(i)
    END DO
    y=A%PDiag*y
    DO i=A%n,1,-1
      y(i)=y(i)/A%PDiag(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*y(i)
      END DO
    END DO
    x=x+y
  END DO

  DO iz=1,nz
    DO iy=1,ny
      DO ix=1,nx
        x3(ix,iy,iz)=x(IndexC(ix,iy,iz,nx,ny,nz))
      END DO
    END DO
  END DO

END SUBROUTINE ICSolveMI3_SpRowColL

SUBROUTINE ICSolveM_SpRowColL(x,A,b,NumIter)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A = L + D + U
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: x(:)
  REAL(RealKind) :: b(:)
  INTEGER :: NumIter

  INTEGER :: i,jj,Iter
  REAL(RealKind) :: y(A%n)

  DO Iter=1,NumIter
    DO i=1,A%n
      y(i)=b(i)-A%Diag(i)*x(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
        y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*x(i)
      END DO
    END DO

    DO i=1,A%n
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
      END DO
      y(i)=y(i)/A%PDiag(i)
    END DO
    y=A%PDiag*y
    DO i=A%n,1,-1
      y(i)=y(i)/A%PDiag(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*y(i)
      END DO
    END DO
    x(A%n)=x(A%n)+y(A%n)
  END DO


END SUBROUTINE ICSolveM_SpRowColL

SUBROUTINE ICSolve_SpRowColL(x,A)

!
! Solves (L+D)*D**(-1)*(D+U)x=x
! where A = L + D + U
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: x(:)

  INTEGER :: i,jj
  REAL(RealKind) :: b(1:A%n)

  DO i=1,A%n
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      x(i)=x(i)-A%Val(jj)*x(A%ColInd(jj))
    END DO
    x(i)=x(i)/A%PDiag(i)
  END DO
  b(1:A%n)=A%PDiag(1:A%n)*x(1:A%n)
  DO i=A%n,1,-1
    x(i)=b(i)/A%PDiag(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      b(A%ColInd(jj))=b(A%ColInd(jj))-A%Val(jj)*x(i)
    END DO
  END DO


END SUBROUTINE ICSolve_SpRowColL

SUBROUTINE SpAVecI3_SpRowColL(y3,A,x3)

!
! Computes y=A*x
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: y3(:,:,:)
  REAL(RealKind) :: x3(:,:,:)

  INTEGER :: i,jj
  REAL(RealKind) :: y(SIZE(y3))
  REAL(RealKind) :: x(SIZE(x3))

  nx=SIZE(x3,1)
  ny=SIZE(x3,2)
  nz=SIZE(x3,3)
  DO iz=1,nz
    DO iy=1,ny
      DO ix=1,nx
        x(IndexC(ix,iy,iz,nx,ny,nz))=x3(ix,iy,iz)
      END DO
    END DO
  END DO

  DO i=1,A%n
    y(i)=A%Diag(i)*x(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)+A%Val(jj)*x(A%ColInd(jj))
      y(A%ColInd(jj))=y(A%ColInd(jj))+A%Val(jj)*x(i)
    END DO
  END DO
  DO iz=1,nz
    DO iy=1,ny
      DO ix=1,nx
        y3(ix,iy,iz)=y(IndexC(ix,iy,iz,nx,ny,nz))
      END DO
    END DO
  END DO

END SUBROUTINE SpAVecI3_SpRowColL

SUBROUTINE SpAVec_SpRowColL(y,A,x)

!
! Computes y=A*x
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: y(:)
  REAL(RealKind) :: x(:)

  INTEGER :: i,jj

  DO i=1,A%n
    y(i)=A%Diag(i)*x(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)+A%Val(jj)*x(A%ColInd(jj))
      y(A%ColInd(jj))=y(A%ColInd(jj))+A%Val(jj)*x(i)
    END DO
  END DO

END SUBROUTINE SpAVec_SpRowColL

SUBROUTINE Axpy_SpRowColL(A,x,y)

!
! Computes y=A*x+y
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: x(:)
  REAL(RealKind) :: y(:)

  INTEGER :: i,jj

  DO i=1,A%n
    y(i)=y(i)+A%Diag(i)*x(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)+A%Val(jj)*x(A%ColInd(jj))
      y(A%ColInd(jj))=y(A%ColInd(jj))+A%Val(jj)*x(i)
    END DO
  END DO

END SUBROUTINE Axpy_SpRowColL


SUBROUTINE GaussSeidel_SpRowColL(x,A,b)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A = L + D + U
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: x(A%n)
  REAL(RealKind) :: b(A%n)

  INTEGER :: i,jj
  REAL(RealKind) :: y(A%n)

  DO i=1,A%n
    y(i)=b(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)-A%Val(jj)*y(A%ColInd(jj))
    END DO
    y(i)=y(i)/A%Diag(i)
  END DO
  y=A%Diag*y
  DO i=A%n,1,-1
    x(i)=y(i)/A%Diag(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*x(i)
    END DO
  END DO
END SUBROUTINE GaussSeidel_SpRowColL

SUBROUTINE GaussSeidelF_SpRowColL(x3,A,b3,NumIter)

!
! Performs the iteration
!  x = (L+D)**-1*b
! where A = L + D + U
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: x3(:,:,:)
  REAL(RealKind) :: b3(:,:,:)
  INTEGER :: NumIter

  INTEGER :: i,jj
  INTEGER :: Iter
  REAL(RealKind) :: Temp
  REAL(RealKind) :: y(A%n)
  REAL(RealKind) :: x(SIZE(x3))
  REAL(RealKind) :: b(SIZE(b3))

  nx=SIZE(x3,1)
  ny=SIZE(x3,2)
  nz=SIZE(x3,3)
  DO iz=1,nz
    DO iy=1,ny
      DO ix=1,nx
        b(IndexC(ix,iy,iz,nx,ny,nz))=b3(ix,iy,iz)
      END DO
    END DO
  END DO
  DO i=1,A%n
    Temp=b(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      Temp=Temp-A%Val(jj)*x(A%ColInd(jj))
    END DO
    x(i)=Temp/A%Diag(i)
  END DO
  y(1:A%n)=b(1:A%n)
  DO i=A%n,1,-1
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
    END DO
    x(i)=y(i)/A%Diag(i)
    DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
      y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*x(i)
    END DO
  END DO

  DO Iter=1,NumIter-1
    DO i=1,A%n
      y(i)=b(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*x(i)
      END DO
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
      END DO
      x(i)=y(i)/A%Diag(i)
    END DO

    y(1:A%n)=b(1:A%n)
    DO i=A%n,1,-1
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
      END DO
      x(i)=y(i)/A%Diag(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*x(i)
      END DO
    END DO
  END DO
  DO iz=1,nz
    DO iy=1,ny
      DO ix=1,nx
        x3(ix,iy,iz)=x(IndexC(ix,iy,iz,nx,ny,nz))
      END DO
    END DO
  END DO
END SUBROUTINE GaussSeidelF_SpRowColL
SUBROUTINE GaussSeidelB_SpRowColL(x3,A,b3,NumIter)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A = L + D + U
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  REAL(RealKind) :: x3(:,:,:)
  REAL(RealKind) :: b3(:,:,:)
  INTEGER :: NumIter

  INTEGER :: i,jj
  INTEGER :: Iter
  REAL(RealKind) :: y(A%n)
  REAL(RealKind) :: x(SIZE(x3))
  REAL(RealKind) :: b(SIZE(b3))

  nx=SIZE(x3,1)
  ny=SIZE(x3,2)
  nz=SIZE(x3,3)
  DO iz=1,nz
    DO iy=1,ny
      DO ix=1,nx
        b(IndexC(ix,iy,iz,nx,ny,nz))=b3(ix,iy,iz)
      END DO
    END DO
  END DO

  DO Iter=1,NumIter
    DO i=1,A%n
      y(i)=b(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*x(i)
      END DO
    END DO
    DO i=1,A%n
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
      END DO
      x(i)=y(i)/A%Diag(i)
    END DO

    y(1:A%n)=b(1:A%n)
    DO i=A%n,1,-1
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(i)=y(i)-A%Val(jj)*x(A%ColInd(jj))
      END DO
      x(i)=y(i)/A%Diag(i)
      DO jj=A%RowPtr(i),A%RowPtr(i+1)-1
        y(A%ColInd(jj))=y(A%ColInd(jj))-A%Val(jj)*x(i)
      END DO
    END DO
  END DO
  DO iz=1,nz
    DO iy=1,ny
      DO ix=1,nx
        x3(ix,iy,iz)=x(IndexC(ix,iy,iz,nx,ny,nz))
      END DO
    END DO
  END DO
END SUBROUTINE GaussSeidelB_SpRowColL

SUBROUTINE SpOutput_SpRowColL(A,FileName)

!
! Output 
!
  IMPLICIT NONE

  TYPE(SpRowColL) :: A
  CHARACTER(80) :: FileName

  INTEGER :: i,jj

  OPEN(UNIT=10,FILE=TRIM(FileName)//'Mat',STATUS='UNKNOWN')
  DO i=1,A%m
    WRITE(10,*) i
    WRITE(10,'(8D12.4)') (A%Val(jj),jj=A%RowPtr(i),A%RowPtr(i+1)-1),A%Diag(i)
  END DO
  CLOSE(10)
END SUBROUTINE SpOutput_SpRowColL

END MODULE MatSpRowColL_Mod
  

