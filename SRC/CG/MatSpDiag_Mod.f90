MODULE MatSpDiag_Mod

  USE Kind_Mod

  IMPLICIT NONE

  TYPE SpDiag
    INTEGER :: m,n
    INTEGER :: NumDiag
    INTEGER :: Diag
    INTEGER, POINTER :: DiagPtr(:)=>NULL()
    REAL(RealKind), POINTER :: Val(:,:)=>NULL()
  END TYPE SpDiag
  INTEGER :: IterG=5

CONTAINS

SUBROUTINE Copy_MatSpDiag(A,B)

  TYPE(SpDiag), INTENT(INOUT) :: A
  TYPE(SpDiag), INTENT(IN) :: B

  A%n=B%n
  A%m=B%m
  A%NumDiag=B%NumDiag
  A%Diag=B%Diag
  IF (.NOT.ASSOCIATED(A%DiagPtr)) THEN
    ALLOCATE(A%DiagPtr(A%NumDiag))
  END IF
  A%DiagPtr=B%DiagPtr
  IF (.NOT.ASSOCIATED(A%Val)) THEN
    ALLOCATE(A%Val(A%m,A%NumDiag))
  END IF
  A%Val=B%Val
END SUBROUTINE Copy_MatSpDiag

SUBROUTINE SpAllocate_SpDiag(A)

!
! Allocate A
!

  TYPE(SpDiag) :: A

  IF (.NOT.ASSOCIATED(A%DiagPtr)) THEN
    ALLOCATE(A%DiagPtr(A%NumDiag))
  END IF
  IF (.NOT.ASSOCIATED(A%Val)) THEN
    ALLOCATE(A%Val(A%m,A%NumDiag))
  END IF

END SUBROUTINE SpAllocate_SpDiag


SUBROUTINE SpOutput_SpDiag(A,FileName)

!
! Nullify Pointers of A
!

  TYPE(SpDiag) :: A
  CHARACTER(80) :: FileName

  INTEGER :: i,j

  OPEN(UNIT=10,FILE=TRIM(FileName)//'Mat',STATUS='UNKNOWN')
  DO i=1,A%n
    WRITE(10,*) (A%Val(i,j),j=1,A%NumDiag)
  END DO
  WRITE(10,*) (A%DiagPtr(j), j=1,A%NumDiag)
END SUBROUTINE SpOutput_SpDiag


SUBROUTINE SpNullify_SpDiag(A)

!
! Nullify Pointers of A
!

  TYPE(SpDiag) :: A

  NULLIFY(A%DiagPtr)
  NULLIFY(A%Val)

END SUBROUTINE SpNullify_SpDiag

SUBROUTINE SpAVec_SpDiag(y,A,x)

!
! Computes y=A*x
!

  TYPE(SpDiag) :: A
  REAL(RealKind) :: y(:)
  REAL(RealKind) :: x(:)

  INTEGER :: i,iDiag,Shift
! REAL(RealKind) :: Temp(1:A%n)

! WRITE(*,*) 'SpAVec_SpDiag: x,y,n',x,y,A%n
  y(1:A%n)=0.0e0
  DO iDiag=1,A%NumDiag
    Shift=A%DiagPtr(iDiag)
    DO i=1+MAX(0,Shift),A%n+MIN(0,Shift)
      y(i-Shift)=y(i-Shift)+A%Val(i,iDiag)*x(i)
    END DO
  END DO
END SUBROUTINE SpAVec_SpDiag


SUBROUTINE SpABVec_SpDiag(y,A,B,x)

!
! Computes y=A*x+B*x
!

  TYPE(SpDiag) :: A,B
  REAL(RealKind) :: y(:)
  REAL(RealKind) :: x(:)

  INTEGER :: i,iDiag,Shift
! REAL(RealKind) :: Temp(1:A%n)

  y(1:A%n)=0.0e0
  DO iDiag=1,A%NumDiag
    Shift=A%DiagPtr(iDiag)
    DO i=1+MAX(0,Shift),A%n+MIN(0,Shift)
      y(i-Shift)=y(i-Shift)+A%Val(i,iDiag)*x(i)
    END DO
  END DO
  DO iDiag=1,B%NumDiag
    Shift=B%DiagPtr(iDiag)
    DO i=1+MAX(0,Shift),B%n+MIN(0,Shift)
      y(i-Shift)=y(i-Shift)+B%Val(i,iDiag)*x(i)
    END DO
  END DO
END SUBROUTINE SpABVec_SpDiag


SUBROUTINE SpAVec1_SpDiag(y,A)

!
! Computes y=A*y-y
!

  TYPE(SpDiag) :: A
  REAL(RealKind) :: y(:)

  INTEGER :: i,iDiag,Shift
  REAL(RealKind) :: Temp(A%n)

  Temp(1:A%n)=y(1:A%n)
  y(1:A%n)=-y(1:A%n)
  DO iDiag=1,A%NumDiag
    Shift=A%DiagPtr(iDiag)
    DO i=1+MAX(0,Shift),A%n+MIN(0,Shift)
      y(i-Shift)=y(i-Shift)+A%Val(i,iDiag)*Temp(i)
    END DO
  END DO
END SUBROUTINE SpAVec1_SpDiag

SUBROUTINE SpABVec1_SpDiag(y,A,B)

!
! Computes y=A*y+B*y-y
!

  TYPE(SpDiag) :: A,B
  REAL(RealKind) :: y(:)

  INTEGER :: i,iDiag,Shift
  REAL(RealKind) :: Temp(A%n)

  Temp(1:A%n)=y(1:A%n)
  y(1:A%n)=-y(1:A%n)
  DO iDiag=1,A%NumDiag
    Shift=A%DiagPtr(iDiag)
    DO i=1+MAX(0,Shift),A%n+MIN(0,Shift)
      y(i-Shift)=y(i-Shift)+A%Val(i,iDiag)*Temp(i)
    END DO
  END DO
  DO iDiag=1,B%NumDiag
    Shift=B%DiagPtr(iDiag)
    DO i=1+MAX(0,Shift),B%n+MIN(0,Shift)
      y(i-Shift)=y(i-Shift)+B%Val(i,iDiag)*Temp(i)
    END DO
  END DO
END SUBROUTINE SpABVec1_SpDiag


SUBROUTINE GaussSeidel_SpDiag(x,A,b)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A = L + D + U
!


  TYPE(SpDiag) :: A
  REAL(RealKind) :: x(:)
  REAL(RealKind) :: b(:)

  INTEGER :: i,iDiag,iShift,Iter,imax
  REAL(RealKind), ALLOCATABLE :: Temp(:)

  ALLOCATE(Temp(1+A%DiagPtr(1):A%n-A%DiagPtr(1)))
  Temp(1:A%n)=b(1:A%n)
  DO i=1,A%n
    Temp(i)=Temp(i)/A%Val(i,A%Diag)
    DO iDiag=1,A%Diag-1 
      Temp(i-A%DiagPtr(iDiag))= &
       Temp(i-A%DiagPtr(iDiag))-A%Val(i,iDiag)*Temp(i)
    END DO
  END DO

  DO Iter=1,IterG 
    DO i=1,A%n
      DO  iDiag=A%Diag+1,A%NumDiag
        Temp(i-A%DiagPtr(iDiag))= &
         Temp(i-A%DiagPtr(iDiag))-A%Val(i,iDiag)*Temp(i)
      END DO
      Temp(i)=b(i)
    END DO
    DO i=1,A%n
      Temp(i)=Temp(i)/A%Val(i,A%Diag)
      DO iDiag=1,A%Diag-1
        Temp(i-A%DiagPtr(iDiag))= &
         Temp(i-A%DiagPtr(iDiag))-A%Val(i,iDiag)*Temp(i)
      END DO
    END DO
  END DO
  x(1:A%n)=Temp(1:A%n)
  DEALLOCATE(Temp)

END SUBROUTINE GaussSeidel_SpDiag

SUBROUTINE GaussSeidel_SpDiagAB(x,A,C,b)

!
! Solves (L+D)*D**(-1)*(D+U)x=b
! where A + C = L + D + U
!


  TYPE(SpDiag) :: A,C
  REAL(RealKind) :: x(:)
  REAL(RealKind) :: b(:)

  INTEGER :: i,iDiag,iShift,Iter,imax
  REAL(RealKind), ALLOCATABLE :: Temp(:)


  ALLOCATE(Temp(1+A%DiagPtr(1):A%n-A%DiagPtr(1)))
  Temp(1:A%n)=b(1:A%n)
  DO i=1,A%n
    Temp(i)=Temp(i)/(A%Val(i,A%Diag)+C%Val(i,C%Diag))
    DO iDiag=1,A%Diag-1
      Temp(i-A%DiagPtr(iDiag))= &
       Temp(i-A%DiagPtr(iDiag))-A%Val(i,iDiag)*Temp(i)
    END DO
    DO iDiag=1,C%Diag-1
      Temp(i-C%DiagPtr(iDiag))= &
       Temp(i-C%DiagPtr(iDiag))-C%Val(i,iDiag)*Temp(i)
    END DO
  END DO

  DO Iter=1,IterG
    DO i=1,A%n
      DO  iDiag=A%Diag+1,A%NumDiag
        Temp(i-A%DiagPtr(iDiag))= &
         Temp(i-A%DiagPtr(iDiag))-A%Val(i,iDiag)*Temp(i)
      END DO
      DO  iDiag=C%Diag+1,C%NumDiag
        Temp(i-C%DiagPtr(iDiag))= &
         Temp(i-C%DiagPtr(iDiag))-C%Val(i,iDiag)*Temp(i)
      END DO
      Temp(i)=b(i)
    END DO
    DO i=1,A%n
      Temp(i)=Temp(i)/(A%Val(i,A%Diag)+C%Val(i,C%Diag))
      DO iDiag=1,A%Diag-1
        Temp(i-A%DiagPtr(iDiag))= &
         Temp(i-A%DiagPtr(iDiag))-A%Val(i,iDiag)*Temp(i)
      END DO
      DO iDiag=1,C%Diag-1
        Temp(i-C%DiagPtr(iDiag))= &
         Temp(i-C%DiagPtr(iDiag))-C%Val(i,iDiag)*Temp(i)
      END DO
    END DO
  END DO
  x(1:A%n)=Temp(1:A%n)
  DEALLOCATE(Temp)

END SUBROUTINE GaussSeidel_SpDiagAB


END MODULE MatSpDiag_Mod


