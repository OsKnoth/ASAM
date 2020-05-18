MODULE LinearSolver_Mod

  USE DataType_Mod
  IMPLICIT NONE

CONTAINS

SUBROUTINE MultVec(dt,LU,Pivot,VecVelC)

  REAL(8) :: dt
  REAL(8) :: LU(:,:)
  INTEGER :: Pivot(:)
  TYPE (VecVector4Cell_T), POINTER :: VecVelC(:)

  INTEGER :: i,j,n
  TYPE (VecVector4Cell_T), POINTER :: VecVelCC(:)

  n=UBOUND(A,1)
  DO i=1,n
    CALL ScaleV(dt*LU(iStage,iStage),VecVelC(i)%Vec)
    DO j=i+1,n
      CALL Axpy(dt*LU(i,j),VecVelC(jStage)%Vec,VecVelC(iStage)%Vec)
    END DO
  END DO
  DO i=n,1,-1
    DO j=1,i-1
      CALL Axpy(dt*LU(i,j),VecVelC(jStage)%Vec,VecVelC(iStage)%Vec)
    END DO
  END DO

END SUBROUTINE MultVec

SUBROUTINE MultVecF(dt,LU,Pivot,VecVelC)

  REAL(8) :: dt
  REAL(8) :: LU(:,:)
  INTEGER :: Pivot(:)
  TYPE (VecVelocityFace_T), POINTER :: VecVelF(:)

  INTEGER :: i,j,n
  TYPE (VecVelocityFace_T), POINTER :: VecVelFF(:)

  n=UBOUND(A,1)
  DO i=1,n
    CALL ScaleV(dt*LU(iStage,iStage),VecVelF(i)%Vec)
    DO j=i+1,n
      CALL Axpy(dt*LU(i,j),,VecVelF(jStage)%VecF,VecVelF(iStage)%VecF)
    END DO
  END DO
  DO i=n,1,-1
    DO j=1,i-1
      CALL Axpy(dt*LU(i,j),VecVelF(jStage)%VecF,VecVelF(iStage)%VecF)
    END DO
  END DO

END SUBROUTINE MultVecF

SUBROUTINE DGEFAP(A,IPVT)
  REAL(8) :: A(:,:) 
  INTEGER :: IPVT(:)
 
  REAL(8) :: temp 
  INTEGER :: i,j,l,k,n

  n=SIZE(A,1)
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
  WRITE(*,*) 'n',n
  DO k=1,n
    temp=A(k,k)  
    l=k
    DO i=k+1,n
      IF (ABS(temp)<ABS(A(i,k))) THEN
        l=i
        temp=A(i,k)
      END IF
    END DO
    IPVT(k)=l
    IF (A(l,k)/=0.0d0) THEN
      IF (l/=k) THEN
        DO j=k,n
          temp=A(l,j)
          A(l,j)=A(k,j)
          A(k,j)=temp
        END DO
      END IF
      temp=1.0d0/A(k,k)
      DO i=k+1,n
        A(i,k)=temp*A(i,k)
      END DO
      DO j=k+1,n
        temp=A(k,j)
        DO i=k+1,n
          A(i,j)=A(i,j)-temp*A(i,k)
        END DO
      END DO
    END IF
  END DO
END SUBROUTINE DGEFAP
END MODULE LinearSolver_Mod

PROGRAM tt
  USE LinearSolver_Mod
  IMPLICIT NONE
  REAL(8) :: A(4,4)
  REAL(8) :: LU(4,4)
  INTEGER :: ipvt(4)
  TYPE R_T
    REAL(8), POINTER :: x(:)=>NULL()
  END TYPE R_T
  TYPE(R_T), POINTER :: w(:)
  INTEGER :: i,n

  A=0.0d0
  A(2,1)=1.0d0
  A(3,1)=2.0d0
  A(4,1)=3.0d0
  ALLOCATE(w(4))
  n=4
  DO i=1,n
    ALLOCATE(w(i)%x(10))
    w(i)%x=1
  END DO
  LU=A
  CALL DGEFAP(LU,ipvt)
  CALL Mult(A,LU,ipvt,w)

END PROGRAM tt


