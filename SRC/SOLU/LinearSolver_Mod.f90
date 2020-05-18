MODULE LinearSolver_Mod

  USE DataType_Mod
  IMPLICIT NONE

CONTAINS

SUBROUTINE DGEFAP(A,IPVT)
  REAL(RealKind) :: A(:,:) 
  INTEGER :: IPVT(:)
 
  REAL(RealKind) :: temp 
  INTEGER :: i,j,l,k,n,itemp

  n=SIZE(A,1)
!
!     GAUSSIAN ELIMINATION WITH PARTIAL PIVOTING
!
  DO k=1,n
    IPVT(k)=k
  END DO
  DO k=1,n
    temp=A(k,k)  
    l=k
    DO i=k+1,n
      IF (ABS(temp)<ABS(A(i,k))) THEN
        l=i
        temp=A(i,k)
      END IF
    END DO
    IF (A(l,k)/=0.0d0) THEN
      IF (l/=k) THEN
        itemp=IPVT(l)
        IPVT(l)=IPVT(k)
        IPVT(k)=itemp
!       DO j=k,n
        DO j=1,n
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



