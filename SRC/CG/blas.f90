FUNCTION DNRM2(n,x,incx)

  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DNRM2
  INTEGER :: n,incx
  REAL(RealKind) :: x(n)

  DNRM2=SQRT(SUM(x*x))

END FUNCTION DNRM2

FUNCTION DDOT(n,x,incx,y,incy)

  USE Kind_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DDOT
  INTEGER :: n,incx,incy
  REAL(RealKind) :: x(n),y(n) 

  DDOT=SUM(x*y)

END FUNCTION DDOT   
SUBROUTINE DAXPY(n,alpha,x,incx,y,incy)

  USE Kind_Mod
  IMPLICIT NONE
  INTEGER :: n,incx,incy
  REAL(RealKind) :: alpha,x(n),y(n) 

  y=alpha*x+y

END SUBROUTINE DAXPY
SUBROUTINE DCOPY(n,x,incx,y,incy)

  USE Kind_Mod
  IMPLICIT NONE
  INTEGER :: n,incx,incy
  REAL(RealKind) :: alpha,x(n),y(n) 

  y=x

END SUBROUTINE DCOPY
   
