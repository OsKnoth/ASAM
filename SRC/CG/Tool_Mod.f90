MODULE Tool_Mod

  USE Kind_Mod
  USE Index_Mod

  INTERFACE Axpy
    MODULE PROCEDURE Axpy31
  END INTERFACE 

CONTAINS

SUBROUTINE Axpy31(alpha,c3,c1)

  REAL(RealKind) :: alpha
  REAL(RealKind), TARGET :: c3(:,:,:)
  REAL(RealKind) :: c1(:)

  INTEGER :: i,j,k,i1,l,m,n
  REAL(RealKind), POINTER :: c3P(:,:,:)

  c3P=>c3
  m=SIZE(c3P,1)
  n=SIZE(c3P,2)
  l=SIZE(c3P,3)
  DO k=1,l
    DO j=1,n
      DO i=1,m
        i1=IndexC(i,j,k,m,n,l)
        c1(i1)=alpha*c1(i1)+c3P(i,j,k)
      END DO
    END DO
  END DO
END SUBROUTINE Axpy31
END MODULE Tool_Mod

