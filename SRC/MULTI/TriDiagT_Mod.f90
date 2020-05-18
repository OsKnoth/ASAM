MODULE TriT_1D_Mod
  USE Matrix_Mod
  USE Parameter_Mod
  USE Control_Mod
  USE Reverse_Mod

  IMPLICIT NONE

  TYPE TriTDiag1D_T
    INTEGER :: n1=0
    REAL(RealKind), POINTER :: d(:)=>NULL()
    REAL(RealKind), POINTER :: u(:)=>NULL()
    REAL(RealKind), POINTER :: l(:)=>NULL()
    TYPE(TriTDiag1D_T), POINTER :: Coarse=>NULL()
    TYPE(TriTDiag1D_T), POINTER :: Fine=>NULL()
  END TYPE TriTDiag1D_T

  INTERFACE Allocate
    MODULE PROCEDURE Allocate_1D
  END INTERFACE Allocate
  INTERFACE Solve
    MODULE PROCEDURE Solve_1D
  END INTERFACE Solve
  INTERFACE MatVecT
    MODULE PROCEDURE MatVecT_1D
  END INTERFACE MatVecT
  INTERFACE Compute
    MODULE PROCEDURE ComputeDual_1D,ComputeFull_1D
  END INTERFACE Compute
CONTAINS

SUBROUTINE ComputeDual_1D(A,F,Vol,DTU,DUU,DUT)

  TYPE (TriTDiag1D_T) :: A
  REAL(RealKind) :: F(:),Vol(:)
  REAL(RealKind) :: DTU(:),DUU(:),DUT(:)

  INTEGER :: i,n1

! a_1 u_1
! l_1 a_2 u_2
!     l_2 a_3 u_3


  n1=A%n1
  A%d(1)=A%d(1)+(beta0*dtp)**2*Fac1_L*F(1)*DTU(1)*DUU(1)*F(1)*DUT(1)/(Vol(1)+Eps)
  DO i=1,n1-1
    A%u(i)=-(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*F(i+1)*DUT(i+1)/(Vol(i)+Vol(i+1)+Eps)
    A%l(i)=-(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*F(i+1)*DUT(i)/(Vol(i)+Vol(i+1)+Eps)
!   A%d(i)=A%d(i)-A%l(i)
    A%d(i)=A%d(i)+(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*F(i+1)*DUT(i)/(Vol(i)+Vol(i+1)+Eps)
!   A%d(i+1)=A%d(i+1)-A%u(i)
    A%d(i+1)=A%d(i+1)+(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*F(i+1)*DUT(i+1)/(Vol(i)+Vol(i+1)+Eps)
  END DO
  A%d(n1)=A%d(n1)+(beta0*dtp)**2*Fac1_R*F(n1+1)*DTU(n1+1)*DUU(n1+1)*F(n1+1)*DUT(n1)/(Vol(n1)+Eps)
  DO i=1,n1
    A%d(i)=A%d(i)+FacAnela*Vol(i)
  END DO

END SUBROUTINE ComputeDual_1D

SUBROUTINE ComputeFull_1D(A,F,dd,Metr1,Metr2,Vol,DTU,DUU,DUT)

  TYPE (TriTDiag1D_T) :: A
  REAL(RealKind) :: F(:),dd(:),Metr1,Metr2,Vol(:)
  REAL(RealKind) :: DTU(:),DUU(:),DUT(:)

  INTEGER :: i,n1

! a_1 u_1
! l_1 a_2 u_2
!     l_2 a_3 u_3


  n1=A%n1
  A%d(1)=A%d(1)+(beta0*dtp)**2*Fac1_L*F(1)*DTU(1)*DUU(1)*DUT(1)/(dd(1)*Metr1*Metr2)
  DO i=1,n1-1
    A%u(i)=-(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*DUT(i+1)/((dd(i)+dd(i+1))*Metr1*Metr2)
    A%l(i)=-(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*DUT(i)/((dd(i)+dd(i+1))*Metr1*Metr2)
!   A%d(i)=A%d(i)-A%l(i)
    A%d(i)=A%d(i)+(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*DUT(i)/((dd(i)+dd(i+1))*Metr1*Metr2)
!   A%d(i+1)=A%d(i+1)-A%u(i)
    A%d(i+1)=A%d(i+1)+(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*DUT(i+1)/((dd(i)+dd(i+1))*Metr1*Metr2)
  END DO
  A%d(n1)=A%d(n1)+(beta0*dtp)**2*Fac1_R*F(n1+1)*DTU(n1+1)*DUU(n1+1)*DUT(n1)/(dd(n1)*Metr1*Metr2)
  DO i=1,n1
    A%d(i)=A%d(i)+FacAnela*Vol(i)
  END DO

END SUBROUTINE ComputeFull_1D

SUBROUTINE Allocate_1D(A,n)

  TYPE (TriTDiag1D_T) :: A
  INTEGER :: i,n

  A%n1=n
  IF (.NOT.ASSOCIATED(A%d)) THEN
    ALLOCATE(A%d(A%n1))
  END IF
  A%d=1.d-40
  IF (.NOT.ASSOCIATED(A%u)) THEN
    ALLOCATE(A%u(A%n1-1))
  END IF
  A%u=Zero
  IF (.NOT.ASSOCIATED(A%l)) THEN
    ALLOCATE(A%l(A%n1-1))
  END IF
  A%l=Zero

END SUBROUTINE Allocate_1D

SUBROUTINE MatVecT_1D(y,A,x)

  REAL(RealKind) :: y(:)
  TYPE(TriTDiag1D_T) :: A
  REAL(RealKind) :: x(:)
 
  REAL(RealKind) :: ss

  INTEGER :: i
  
  IF (A%n1>1) THEN
    y(1)=A%d(1)*x(1)+A%u(1)*x(2)
    DO i=2,A%n1-1
      y(i)=A%l(i-1)*x(i-1)+A%d(i)*x(i)+A%u(i)*x(i+1)
    END DO
    y(A%n1)=A%l(A%n1-1)*x(A%n1-1)+A%d(A%n1)*x(A%n1)
  ELSE
    y(1)=A%d(1)*x(1)
  END IF

END SUBROUTINE MatVecT_1D

SUBROUTINE Solve_1D(A,b)

  TYPE(TriTDiag1D_T) :: A
  REAL(RealKind) :: b(:)

  REAL(RealKind) :: w(A%n1)
  REAL(RealKind) :: t
  INTEGER :: i

  t=A%d(1)
  b(1)=b(1)/t
  DO i=2,A%n1
    w(i-1)=A%u(i-1)/t
    t=A%d(i)-A%l(i-1)*w(i-1)
    b(i)=(b(i)-A%l(i-1)*b(i-1))/t
  END DO
  DO i=A%n1-1,1,-1
    b(i)=b(i)-w(i)*b(i+1)
  END DO

END SUBROUTINE Solve_1D

END MODULE TriT_1D_Mod

MODULE TriT_2D_Mod

  USE TriT_1D_Mod

  IMPLICIT NONE

  TYPE TriTDiag2D_T
    INTEGER :: n1=0
    INTEGER :: n2=0
    TYPE(TriTDiag1D_T), POINTER :: d(:)=>NULL()
    TYPE(TriTDiag2D_T), POINTER :: Coarse=>NULL()
    TYPE(TriTDiag2D_T), POINTER :: Fine=>NULL()
    REAL(RealKind), POINTER :: u(:,:)=>NULL()
    REAL(RealKind), POINTER :: l(:,:)=>NULL()
    REAL(RealKind), POINTER :: Res(:,:)=>NULL()
    REAL(RealKind), POINTER :: ResCoarse(:,:)=>NULL()
    REAL(RealKind), POINTER :: xCoarse(:,:)=>NULL()
  END TYPE TriTDiag2D_T
  TYPE TriTDiag2DP_T
    TYPE(TriTDiag2D_T), POINTER :: P=>NULL()
  END TYPE TriTDiag2DP_T


  INTERFACE Allocate
    MODULE PROCEDURE Allocate_2D
  END INTERFACE Allocate
  INTERFACE GaussSeidel
    MODULE PROCEDURE GaussSeidel_2D
  END INTERFACE GaussSeidel
  INTERFACE MatVecT
    MODULE PROCEDURE MatVecT_2D
  END INTERFACE MatVecT
  INTERFACE RestrictT
    MODULE PROCEDURE RestrictT_2D
  END INTERFACE RestrictT
  INTERFACE ProlongateT
    MODULE PROCEDURE ProlongateT_2D
  END INTERFACE ProlongateT
  INTERFACE Compute
    MODULE PROCEDURE ComputeDual_2D,ComputeFull_2D
  END INTERFACE Compute
  INTERFACE Multigrid
    MODULE PROCEDURE Multigrid_2D
  END INTERFACE Multigrid
  INTERFACE Solve
    MODULE PROCEDURE Solve_2D
  END INTERFACE Solve


CONTAINS

SUBROUTINE Allocate_2D(A,n2,n1)

  TYPE (TriTDiag2D_T), POINTER :: A
  INTEGER, OPTIONAL :: n1
  INTEGER, OPTIONAL :: n2

  INTEGER :: j

  IF (.NOT.ASSOCIATED(A)) THEN
    ALLOCATE(A)
  END IF
  IF (PRESENT(n2)) THEN
    A%n2=n2
    IF (.NOT.ASSOCIATED(A%d)) THEN
      ALLOCATE(A%d(n2))
    END IF
    IF (PRESENT(n1)) THEN
      A%n1=n1
      DO j=1,n2
        CALL ALLOCATE(A%d(j),n1)
      END DO
      IF (.NOT.ASSOCIATED(A%u)) THEN
        ALLOCATE(A%u(n1,n2-1))
      END IF
      A%u=Zero
      IF (.NOT.ASSOCIATED(A%l)) THEN
        ALLOCATE(A%l(n1,n2-1))
      END IF
      A%l=Zero
      IF (.NOT.ASSOCIATED(A%Res)) THEN
        ALLOCATE(A%Res(n1,n2))
      END IF
      IF (.NOT.ASSOCIATED(A%ResCoarse)) THEN
        ALLOCATE(A%ResCoarse((n1+Coarse1-1)/Coarse1,(n2+Coarse2-1)/Coarse2))
      END IF
      IF (.NOT.ASSOCIATED(A%xCoarse)) THEN
        ALLOCATE(A%xCoarse((n1+Coarse1-1)/Coarse1,(n2+Coarse2-1)/Coarse2))
      END IF
    END IF
  END IF

END SUBROUTINE Allocate_2D

SUBROUTINE ComputeDual_2D(A,F1,F2,Vol,DTU1,DTU2,DUU1,DUU2,DUT)

  TYPE (TriTDiag2D_T) :: A
  REAL(RealKind) :: F1(:,:)
  REAL(RealKind) :: F2(:,:)
  REAL(RealKind) :: Vol(:,:)
  REAL(RealKind) :: DTU1(:,:)
  REAL(RealKind) :: DTU2(:,:)
  REAL(RealKind) :: DUU1(:,:)
  REAL(RealKind) :: DUU2(:,:)
  REAL(RealKind) :: DUT(:,:)

  INTEGER :: i,j,n1,n2

  n1=A%n1
  n2=A%n2
  IF (ASSOCIATED(A%u).AND.ASSOCIATED(A%l)) THEN
    CALL Compute(A%d(1),F1(:,1),Vol(:,1),DTU1(:,1),DUU1(:,1),DUT(:,1))
    DO i=1,n1
      A%d(1)%d(i)=A%d(1)%d(i)+(beta0*dtp)**2*Fac2_L*F2(i,1)*DTU2(i,1)*DUU2(i,1)*F2(i,1)*DUT(i,1)/(Vol(i,1)+Eps)
    END DO
    DO j=1,n2-1
      DO i=1,n1
        A%u(i,j)=-(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
        A%l(i,j)=-(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
!       A%d(j)%d(i)=A%d(j)%d(i)-A%e(i,j)
        A%d(j)%d(i)=A%d(j)%d(i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
      END DO
      CALL Compute(A%d(j+1),F1(:,j+1),Vol(:,j+1),DTU1(:,j+1),DUU1(:,j+1),DUT(:,j+1))
      DO i=1,n1
!       A%d(j+1)%d(i)=A%d(j+1)%d(i)-A%e(i,j)
        A%d(j+1)%d(i)=A%d(j+1)%d(i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
      END DO
    END DO
    DO i=1,n1
      A%d(n2)%d(i)=A%d(n2)%d(i) &
                         +(beta0*dtp)**2*Fac2_R*F2(i,n2+1)*DTU2(i,n2+1)*DUU2(i,n2+1)*F2(i,n2+1)*DUT(i,n2)/(Vol(i,n2)+Eps)
    END DO
  ELSE
    CALL Compute(A%d(1),F1(:,1),Vol(:,1),DTU1(:,1),DUU1(:,1),DUT(:,1))
    DO i=1,n1
      A%d(1)%d(i)=A%d(1)%d(i)+(beta0*dtp)**2*Fac2_L*F2(i,1)*DTU2(i,1)*DUU2(i,1)*F2(i,1)*DUT(i,1)/(Vol(i,1)+Eps)
    END DO
    DO j=1,n2-1
      DO i=1,n1
!       A%d(j)%d(i)=A%d(j)%d(i)-e(i)
        A%d(j)%d(i)=A%d(j)%d(i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
      END DO
      CALL Compute(A%d(j+1),F1(:,j+1),Vol(:,j+1),DTU1(:,j+1),DUU1(:,j+1),DUT(:,j+1))
      DO i=1,n1
!       A%d(j+1)%d(i)=A%d(j+1)%d(i)-e(i)
        A%d(j+1)%d(i)=A%d(j+1)%d(i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
      END DO
    END DO
    DO i=1,n1
      A%d(n2)%d(i)=A%d(n2)%d(i) &
                         +(beta0*dtp)**2*Fac2_R*F2(i,n2+1)*DTU2(i,n2+1)*DUU2(i,n2+1)*F2(i,n2+1)*DUT(i,n2)/(Vol(i,n2)+Eps)
    END DO
  END IF

END SUBROUTINE ComputeDual_2D

SUBROUTINE ComputeFull_2D(A,F1,F2,dd1,Metr12,Metr13,dd2,Metr21,Metr23,Vol,DTU1,DTU2,DUU1,DUU2,DUT)

  TYPE (TriTDiag2D_T) :: A
  REAL(RealKind) :: F1(:,:)
  REAL(RealKind) :: F2(:,:)
  REAL(RealKind) :: dd1(:)
  REAL(RealKind) :: Metr12(:)
  REAL(RealKind) :: Metr13
  REAL(RealKind) :: dd2(:)
  REAL(RealKind) :: Metr21(:)
  REAL(RealKind) :: Metr23
  REAL(RealKind) :: Vol(:,:)
  REAL(RealKind) :: DTU1(:,:)
  REAL(RealKind) :: DTU2(:,:)
  REAL(RealKind) :: DUU1(:,:)
  REAL(RealKind) :: DUU2(:,:)
  REAL(RealKind) :: DUT(:,:)

  INTEGER :: i,j,n1,n2

  n1=A%n1
  n2=A%n2
  IF (ASSOCIATED(A%u).AND.ASSOCIATED(A%l)) THEN
    CALL Compute(A%d(1),F1(:,1),dd1,Metr12(1),Metr13,Vol(:,1),DTU1(:,1),DUU1(:,1),DUT(:,1))
    DO i=1,n1
      A%d(1)%d(i)=A%d(1)%d(i)+(beta0*dtp)**2*Fac2_L*F2(i,1)*DTU2(i,1)*DUU2(i,1)*DUT(i,1) &
                 /(dd2(1)*Metr21(i)*Metr23)
    END DO
    DO j=1,n2-1
      DO i=1,n1
        A%u(i,j)=-(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*DUT(i,j+1) &
                 /((dd2(j)+dd2(j+1))*Metr21(i)*Metr23)
        A%l(i,j)=-(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*DUT(i,j) &
                 /((dd2(j)+dd2(j+1))*Metr21(i)*Metr23)
!       A%d(j)%d(i)=A%d(j)%d(i)-A%e(i,j)
        A%d(j)%d(i)=A%d(j)%d(i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*DUT(i,j) &
            /((dd2(j)+dd2(j+1))*Metr21(i)*Metr23)
      END DO
      CALL Compute(A%d(j+1),F1(:,j+1),dd1(:),Metr12(j+1),Metr13,Vol(:,j+1),DTU1(:,j+1),DUU1(:,j+1),DUT(:,j+1))
      DO i=1,n1
!       A%d(j+1)%d(i)=A%d(j+1)%d(i)-A%e(i,j)
        A%d(j+1)%d(i)=A%d(j+1)%d(i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*DUT(i,j+1) &
            /((dd2(j)+dd2(j+1))*Metr21(i)*Metr23)
      END DO
    END DO
    DO i=1,n1
      A%d(n2)%d(i)=A%d(n2)%d(i) &
                         +(beta0*dtp)**2*Fac2_R*F2(i,n2+1)*DTU2(i,n2+1)*DUU2(i,n2+1)*DUT(i,n2) &
                         /(dd2(n2)*Metr21(i)*Metr23)
    END DO
  ELSE
    CALL Compute(A%d(1),F1(:,1),dd1(:),Metr12(1),Metr13,Vol(:,1),DTU1(:,1),DUU1(:,1),DUT(:,1))
    DO i=1,n1
      A%d(1)%d(i)=A%d(1)%d(i)+(beta0*dtp)**2*Fac2_L*F2(i,1)*DTU2(i,1)*DUU2(i,1)*DUT(i,1) &
                             /(dd2(1)*Metr21(i)*Metr23)
    END DO
    DO j=1,n2-1
      DO i=1,n1
!       A%d(j)%d(i)=A%d(j)%d(i)-e(i)
        A%d(j)%d(i)=A%d(j)%d(i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*DUT(i,j) &
            /((dd2(j)+dd2(j+1))*Metr21(i)*Metr23)
      END DO
      CALL Compute(A%d(j+1),F1(:,j+1),dd1(:),Metr12(j+1),Metr13,Vol(:,j+1),DTU1(:,j+1),DUU1(:,j+1),DUT(:,j+1))
      DO i=1,n1
!       A%d(j+1)%d(i)=A%d(j+1)%d(i)-e(i)
        A%d(j+1)%d(i)=A%d(j+1)%d(i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*DUT(i,j+1) &
            /((dd2(j)+dd2(j+1))*Metr21(i)*Metr23)
      END DO
    END DO
    DO i=1,n1
      A%d(n2)%d(i)=A%d(n2)%d(i) &
                         +(beta0*dtp)**2*Fac2_R*F2(i,n2+1)*DTU2(i,n2+1)*DUU2(i,n2+1)*DUT(i,n2) &
                         /(dd2(n2)*Metr21(i)*Metr23)
    END DO
  END IF

END SUBROUTINE ComputeFull_2D

SUBROUTINE MatVecT_2D(y,A,x)

  REAL(RealKind) :: y(:,:)
  TYPE(TriTDiag2D_T) :: A
  REAL(RealKind) :: x(:,:)

  INTEGER :: i,j,n1,n2

  n1=A%n1
  n2=A%n2
  CALL MatVecT(y(:,1),A%d(1),x(:,1))
  IF (n2>1) THEN
    DO i=1,n1
      y(i,1)=y(i,1)+A%u(i,1)*x(i,2)    
    END DO
    DO j=2,n2-1
      CALL MatVecT(y(:,j),A%d(j),x(:,j))
      DO i=1,n1
        y(i,j)=y(i,j)+A%l(i,j-1)*x(i,j-1)+A%u(i,j)*x(i,j+1)    
      END DO
    END DO
    CALL MatVecT(y(:,n2),A%d(n2),x(:,n2))
    DO i=1,n1
      y(i,n2)=y(i,n2)+A%l(i,n2-1)*x(i,n2-1)    
    END DO
  END IF

END SUBROUTINE MatVecT_2D

SUBROUTINE GaussSeidel_2D(A,b,x,GaussIter)

  TYPE(TriTDiag2D_T) :: A
  REAL(RealKind) :: b(:,:)
  REAL(RealKind) :: x(:,:)
  INTEGER :: GaussIter

  REAL(RealKind) :: t(A%n1)
  INTEGER :: i,j,iter

  IF (A%n2>1) THEN
    DO iter=1,GaussIter
      t(:)=b(:,1)
      DO j=1,A%n2-1
        DO i=1,A%n1
          x(i,j)=t(i)-A%u(i,j)*x(i,j+1)
        END DO 
        CALL Solve(A%d(j),x(:,j))
        DO i=1,A%n1
          t(i)=b(i,j+1)-A%l(i,j)*x(i,j)
        END DO 
      END DO
      x(:,A%n2)=t(:)
      CALL Solve(A%d(A%n2),x(:,A%n2))
      DO i=1,A%n1
        t(i)=b(i,A%n2-1)-A%u(i,A%n2-1)*x(i,A%n2)
      END DO  
      DO j=A%n2-1,2,-1
        DO i=1,A%n1
          x(i,j)=t(i)-A%l(i,j-1)*x(i,j-1)
        END DO 
        CALL Solve(A%d(j),x(:,j))
        DO i=1,A%n1
          t(i)=b(i,j-1)-A%u(i,j-1)*x(i,j)
        END DO 
      END DO
    END DO
    x(:,1)=t(:)
    CALL Solve(A%d(1),x(:,1))
  ELSE
    x(:,1)=b(:,1)
    CALL Solve(A%d(1),x(:,1))
  END IF
   
END SUBROUTINE GaussSeidel_2D
SUBROUTINE Solve_2D(A,b,x)

  TYPE(TriTDiag2D_T) :: A
  REAL(RealKind) :: b(:,:)
  REAL(RealKind) :: x(:,:)

  INTEGER :: Iter

  
  DO Iter=1,MultIter2
    CALL Multigrid(A,b,x)
  END DO

END SUBROUTINE Solve_2D

SUBROUTINE Multigrid_2D(A,b,x)

  TYPE(TriTDiag2D_T), TARGET :: A
  REAL(RealKind), TARGET :: b(:,:)
  REAL(RealKind), TARGET :: x(:,:)

  INTEGER :: i,j,n1,n2
  TYPE(TriTDiag2D_T), POINTER :: ALoc
  REAL(RealKind), POINTER :: bLoc(:,:)
  REAL(RealKind), POINTER :: xLoc(:,:)
  INTEGER :: GaussIter

  bLoc=>b
  xLoc=>x
  ALoc=>A
  GaussIter=1
  DO 
    n1=ALoc%n1
    n2=ALoc%n2
!   Relax
    CALL GaussSeidel(ALoc,bLoc,xLoc,GaussIter)
    GaussIter=GaussIter+GaussIter2

    IF (ASSOCIATED(ALoc%Coarse)) THEN
!     RestrictTR
      CALL MatVecT(ALoc%Res,ALoc,xLoc)
      ALoc%Res=bLoc-ALoc%Res
      CALL RestrictT(ALoc%Res,ALoc%ResCoarse)
!     Down
      bLoc=>ALoc%ResCoarse
      xLoc=>ALoc%xCoarse
      xLoc=Zero
      ALoc=>ALoc%Coarse
    ELSE
      GaussIter=GaussIter-GaussIter2
      EXIT
    END IF
  END DO
  DO
!   Smooth
    CALL GaussSeidel(ALoc,bLoc,xLoc,GaussIter)
    GaussIter=GaussIter-GaussIter2
    IF (ASSOCIATED(ALoc%Fine)) THEN
!     Up
      ALoc=>ALoc%Fine
      n1=ALoc%n1
      n2=ALoc%n2
      CALL ProlongateT(xLoc,ALoc%Res)
      IF (ASSOCIATED(ALoc%Fine)) THEN
        xLoc=>Aloc%Fine%xCoarse
        bLoc=>ALoc%Fine%ResCoarse
      ELSE
        xLoc=>x
        bLoc=>b
      END IF
      xLoc=xLoc+ALoc%Res
    ELSE
      EXIT
    END IF 
  END DO

END SUBROUTINE Multigrid_2D

SUBROUTINE RestrictT_2D(Res,ResCoarse)

  REAL(RealKind) :: Res(:,:)
  REAL(RealKind) :: ResCoarse(:,:)

  INTEGER :: n1,n2
  INTEGER :: i1,i2,i1C,i2C

  n1=SIZE(Res,1)
  n2=SIZE(Res,2)
  ResCoarse=0.0d0
  DO i1=1,n1
    i1C=(i1+Coarse1-1)/Coarse1
    DO i2=1,n2
      i2C=(i2+Coarse2-1)/Coarse2
      ResCoarse(i1C,i2C)=ResCoarse(i1C,i2C)+Res(i1,i2)
    END DO
  END DO
END SUBROUTINE RestrictT_2D

SUBROUTINE ProlongateT_2D(xCoarse,x)
  REAL(RealKind) :: x(:,:),xCoarse(:,:)

  INTEGER :: n1,n2
  INTEGER :: i1,i2,i1C,i2C

  n1=SIZE(x,1)
  n2=SIZE(x,2)
  DO i1=1,n1
    i1C=(i1+Coarse1-1)/Coarse1
    DO i2=1,n2
      i2C=(i2+Coarse2-1)/Coarse2
      x(i1,i2)=xCoarse(i1C,i2C)
    END DO
  END DO

END SUBROUTINE ProlongateT_2D

END MODULE TriT_2D_Mod

MODULE TriT_3D_Mod

  USE TriT_2D_Mod

  IMPLICIT NONE

  TYPE TriTDiag3D_T
    INTEGER :: n1=0
    INTEGER :: n2=0
    INTEGER :: n3=0
    TYPE(TriTDiag2DP_T), POINTER :: d(:)=>NULL()
    REAL(RealKind), POINTER :: u(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: l(:,:,:)=>NULL()
    TYPE(TriTDiag3D_T), POINTER :: Coarse=>NULL()
    TYPE(TriTDiag3D_T), POINTER :: Fine=>NULL()
    REAL(RealKind), POINTER :: Res(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: ResCoarse(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: xCoarse(:,:,:)=>NULL()
  END TYPE TriTDiag3D_T



  INTERFACE Allocate
    MODULE PROCEDURE Allocate_3D
  END INTERFACE Allocate
  INTERFACE GaussSeidel
    MODULE PROCEDURE GaussSeidel_3D
  END INTERFACE GaussSeidel
  INTERFACE MatVecT
    MODULE PROCEDURE MatVecT_3D
  END INTERFACE MatVecT
  INTERFACE RestrictT
    MODULE PROCEDURE RestrictT_3D
  END INTERFACE RestrictT
  INTERFACE ProlongateT
    MODULE PROCEDURE ProlongateT_3D
  END INTERFACE ProlongateT
  INTERFACE Compute
    MODULE PROCEDURE ComputeDual_3D,ComputeFull_3D
  END INTERFACE Compute
  INTERFACE Multigrid
    MODULE PROCEDURE Multigrid_3D
  END INTERFACE Multigrid
  INTERFACE Solve
    MODULE PROCEDURE Solve_3D
  END INTERFACE Solve


CONTAINS

SUBROUTINE Allocate_3D(A,n3,n2,n1)

  TYPE(TriTDiag3D_T), POINTER :: A
  INTEGER, OPTIONAL :: n1
  INTEGER, OPTIONAL :: n2
  INTEGER, OPTIONAL :: n3

  INTEGER :: i
  
  IF (.NOT.ASSOCIATED(A)) THEN
    ALLOCATE(A)
  END IF
  IF (PRESENT(n3)) THEN
    A%n3=n3
    IF (.NOT.ASSOCIATED(A%d)) THEN
      ALLOCATE(A%d(n3))
    END IF
    IF (PRESENT(n2)) THEN
      A%n2=n2
      DO i=1,n3
        IF (.NOT.ASSOCIATED(A%d(i)%P)) THEN
          ALLOCATE(A%d(i)%P)
        END IF
        CALL ALLOCATE(A%d(i)%P,n2,n1)
      END DO
      IF (PRESENT(n1)) THEN
        A%n1=n1
        IF (.NOT.ASSOCIATED(A%u)) THEN
          ALLOCATE(A%u(n1,n2,n3-1))
        END IF
        A%u=Zero
        IF (.NOT.ASSOCIATED(A%l)) THEN
          ALLOCATE(A%l(n1,n2,n3-1))
        END IF
        A%l=Zero
        IF (.NOT.ASSOCIATED(A%Res)) THEN
          ALLOCATE(A%Res(n1,n2,n3))
        END IF
        IF (.NOT.ASSOCIATED(A%ResCoarse)) THEN
          ALLOCATE(A%ResCoarse((n1+Coarse1-1)/Coarse1,(n2+Coarse2-1)/Coarse2,(n3+Coarse3-1)/Coarse3))
        END IF
        IF (.NOT.ASSOCIATED(A%xCoarse)) THEN
          ALLOCATE(A%xCoarse((n1+Coarse1-1)/Coarse1,(n2+Coarse2-1)/Coarse2,(n3+Coarse3-1)/Coarse3))
        END IF
      END IF
    END IF
  END IF

END SUBROUTINE Allocate_3D

SUBROUTINE ComputeDual_3D(A,F1,F2,F3,Vol,DTU1,DTU2,DTU3,DUU1,DUU2,DUU3,DUT)

  TYPE (TriTDiag3D_T) :: A
  REAL(RealKind) :: F1(:,:,:)
  REAL(RealKind) :: F2(:,:,:)
  REAL(RealKind) :: F3(:,:,:)
  REAL(RealKind) :: Vol(:,:,:)
  REAL(RealKind) :: DTU1(:,:,:)
  REAL(RealKind) :: DTU2(:,:,:)
  REAL(RealKind) :: DTU3(:,:,:)
  REAL(RealKind) :: DUU1(:,:,:)
  REAL(RealKind) :: DUU2(:,:,:)
  REAL(RealKind) :: DUU3(:,:,:)
  REAL(RealKind) :: DUT(:,:,:)

  INTEGER :: i,j,k,n1,n2,n3

  n1=A%n1
  n2=A%n2
  n3=A%n3
  IF (ASSOCIATED(A%u).AND.ASSOCIATED(A%l)) THEN
    CALL Compute(A%d(1)%P,F1(:,:,1),F2(:,:,1),Vol(:,:,1), &
                 DTU1(:,:,1),DTU2(:,:,1),DUU1(:,:,1),DUU2(:,:,1),DUT(:,:,1))
    DO j=1,n2
      DO i=1,n1
        A%d(1)%P%d(j)%d(i)=A%d(1)%P%d(j)%d(i) &
           +(beta0*dtp)**2*Fac3_L*F3(i,j,1)*DTU3(i,j,1)*DUU3(i,j,1)*F3(i,j,1)*DUT(i,j,1)/(Vol(i,j,1)+Eps)
      END DO
    END DO
    DO k=1,n3-1
      DO j=1,n2
        DO i=1,n1
          A%u(i,j,k)=-(beta0*dtp)**2*2.0d0*F3(i,j,k+1) &
              *DTU3(i,j,k+1)*F3(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
!         A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*F3(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
        END DO
      END DO
      CALL Compute(A%d(k+1)%P,F1(:,:,k+1),F2(:,:,k+1),Vol(:,:,k+1), &
                   DTU1(:,:,k+1),DTU2(:,:,k+1),DUU1(:,:,k+1),DUU2(:,:,k+1),DUT(:,:,k+1))
      DO j=1,n2
        DO i=1,n1
          A%l(i,j,k)=-(beta0*dtp)**2*2.0d0*F3(i,j,k+1) &
              *DTU3(i,j,k+1)*F3(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
!         A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*F3(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
        END DO
      END DO
    END DO
    DO j=1,n2
      DO i=1,n1
        A%d(n3)%P%d(j)%d(i)=A%d(n3)%P%d(j)%d(i) &
           +(beta0*dtp)**2*Fac3_R*F3(i,j,n3+1)*DTU3(i,j,n3+1)*F3(i,j,n3+1)*DUT(i,j,n3)/(Vol(i,j,n3)+Eps)
      END DO
    END DO
  ELSE
    n2=A%d(1)%P%n2
    n1=A%d(1)%P%n1
    CALL Compute(A%d(1)%P,F1(:,:,1),F2(:,:,1),Vol(:,:,1), &
                 DTU1(:,:,1),DTU2(:,:,1),DUU1(:,:,1),DUU2(:,:,1),DUT(:,:,1))
    DO j=1,n2
      DO i=1,n1
        A%d(1)%P%d(j)%d(i)=A%d(1)%P%d(j)%d(i) &
           +(beta0*dtp)**2*Fac3_L*F3(i,j,1)*DTU3(i,j,1)*DUU3(i,j,1)*F3(i,j,1)*DUT(i,j,1)/(Vol(i,j,1)+Eps)
      END DO
    END DO
    DO k=1,n3-1
      DO j=1,n2
        DO i=1,n1
!         A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*DUU3(i,j,k+1)*F3(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
        END DO
      END DO
      CALL Compute(A%d(k+1)%P,F1(:,:,k+1),F2(:,:,k+1),Vol(:,:,k+1), &
                   DTU1(:,:,k+1),DTU2(:,:,k+1),DUU1(:,:,k+1),DUU2(:,:,k+1),DUT(:,:,k+1))
      DO j=1,n2
        DO i=1,n1
!         A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*DUU3(i,j,k+1)*F3(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
        END DO
      END DO
    END DO
    DO j=1,n2
      DO i=1,n1
        A%d(n3)%P%d(j)%d(i)=A%d(n3)%P%d(j)%d(i) &
           +(beta0*dtp)**2*Fac3_R*F3(i,j,n3+1)*DTU3(i,j,n3+1)*DUU3(i,j,n3+1)*F3(i,j,n3+1)*DUT(i,j,n3)/(Vol(i,j,n3)+Eps)
      END DO
    END DO
  END IF
  
END SUBROUTINE ComputeDual_3D

SUBROUTINE ComputeFull_3D(A,F1,F2,F3,dd1,Metr12,Metr13,dd2,Metr21,Metr23,dd3,Metr31,Metr32 &
                         ,Vol,DTU1,DTU2,DTU3,DUU1,DUU2,DUU3,DUT)

  TYPE (TriTDiag3D_T) :: A
  REAL(RealKind) :: F1(:,:,:)
  REAL(RealKind) :: F2(:,:,:)
  REAL(RealKind) :: F3(:,:,:)
  REAL(RealKind) :: dd1(:)
  REAL(RealKind) :: Metr12(:)
  REAL(RealKind) :: Metr13(:)
  REAL(RealKind) :: dd2(:)
  REAL(RealKind) :: Metr21(:)
  REAL(RealKind) :: Metr23(:)
  REAL(RealKind) :: dd3(:)
  REAL(RealKind) :: Metr31(:)
  REAL(RealKind) :: Metr32(:)
  REAL(RealKind) :: Vol(:,:,:)
  REAL(RealKind) :: DTU1(:,:,:)
  REAL(RealKind) :: DTU2(:,:,:)
  REAL(RealKind) :: DTU3(:,:,:)
  REAL(RealKind) :: DUU1(:,:,:)
  REAL(RealKind) :: DUU2(:,:,:)
  REAL(RealKind) :: DUU3(:,:,:)
  REAL(RealKind) :: DUT(:,:,:)

  INTEGER :: i,j,k,n1,n2,n3

  n1=A%n1
  n2=A%n2
  n3=A%n3
  IF (ASSOCIATED(A%u).AND.ASSOCIATED(A%l)) THEN
    CALL Compute(A%d(1)%P,F1(:,:,1),F2(:,:,1),dd1,Metr12,Metr13(1),dd2,Metr21,Metr23(1),Vol(:,:,1), &
                 DTU1(:,:,1),DTU2(:,:,1),DUU1(:,:,1),DUU2(:,:,1),DUT(:,:,1))
    DO j=1,n2
      DO i=1,n1
        A%d(1)%P%d(j)%d(i)=A%d(1)%P%d(j)%d(i) &
           +(beta0*dtp)**2*Fac3_L*F3(i,j,1)*DTU3(i,j,1)*DUU3(i,j,1)*DUT(i,j,1)/(dd3(1)*Metr31(i)*Metr32(j))
      END DO
    END DO
    DO k=1,n3-1
      DO j=1,n2
        DO i=1,n1
          A%u(i,j,k)=-(beta0*dtp)**2*2.0d0*F3(i,j,k+1) &
              *DTU3(i,j,k+1)*DUT(i,j,k+1)/((dd3(k)+dd3(k+1))*Metr31(i)*Metr32(j))
!         A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*DUT(i,j,k) &
           /((dd3(k)+dd3(k+1))*Metr31(i)*Metr32(j))
        END DO
      END DO
      CALL Compute(A%d(k+1)%P,F1(:,:,k+1),F2(:,:,k+1),dd1,Metr12,Metr13(k+1),dd2,Metr21,Metr23(k+1),Vol(:,:,k+1), &
                   DTU1(:,:,k+1),DTU2(:,:,k+1),DUU1(:,:,k+1),DUU2(:,:,k+1),DUT(:,:,k+1))
      DO j=1,n2
        DO i=1,n1
          A%l(i,j,k)=-(beta0*dtp)**2*2.0d0*F3(i,j,k+1) &
              *DTU3(i,j,k+1)*DUT(i,j,k)/((dd3(k)+dd3(k+1))*Metr31(i)*Metr32(j))
!         A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*DUT(i,j,k+1) &
           /((dd3(k)+dd3(k+1))*Metr31(i)*Metr32(j))
        END DO
      END DO
    END DO
    DO j=1,n2
      DO i=1,n1
        A%d(n3)%P%d(j)%d(i)=A%d(n3)%P%d(j)%d(i) &
           +(beta0*dtp)**2*Fac3_R*F3(i,j,n3+1)*DTU3(i,j,n3+1)*DUT(i,j,n3) &
           /(dd3(n3)*Metr31(i)*Metr32(j))
      END DO
    END DO
  ELSE
    n2=A%d(1)%P%n2
    n1=A%d(1)%P%n1
    CALL Compute(A%d(1)%P,F1(:,:,1),F2(:,:,1),dd1,Metr12,Metr13(1),dd2,Metr21,Metr23(1),Vol(:,:,1), &
                 DTU1(:,:,1),DTU2(:,:,1),DUU1(:,:,1),DUU2(:,:,1),DUT(:,:,1))
    DO j=1,n2
      DO i=1,n1
        A%d(1)%P%d(j)%d(i)=A%d(1)%P%d(j)%d(i) &
           +(beta0*dtp)**2*Fac3_L*F3(i,j,1)*DTU3(i,j,1)*DUU3(i,j,1)*DUT(i,j,1) &
           /(dd3(1)*Metr31(i)*Metr32(j))
      END DO
    END DO
    DO k=1,n3-1
      DO j=1,n2
        DO i=1,n1
!         A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*DUU3(i,j,k+1)*DUT(i,j,k) &
           /((dd3(k)+dd3(k+1))*Metr31(i)*Metr32(j))
        END DO
      END DO
      CALL Compute(A%d(k+1)%P,F1(:,:,k+1),F2(:,:,k+1),dd1,Metr12,Metr13(k+1),dd2,Metr21,Metr23(k+1),Vol(:,:,k+1), &
                   DTU1(:,:,k+1),DTU2(:,:,k+1),DUU1(:,:,k+1),DUU2(:,:,k+1),DUT(:,:,k+1))
      DO j=1,n2
        DO i=1,n1
!         A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*DUU3(i,j,k+1)*DUT(i,j,k+1) &
           /((dd3(k)+dd3(k+1))*Metr31(i)*Metr32(j))
        END DO
      END DO
    END DO
    DO j=1,n2
      DO i=1,n1
        A%d(n3)%P%d(j)%d(i)=A%d(n3)%P%d(j)%d(i) &
           +(beta0*dtp)**2*Fac3_R*F3(i,j,n3+1)*DTU3(i,j,n3+1)*DUU3(i,j,n3+1)*DUT(i,j,n3) &
           /(dd3(n3)*Metr31(i)*Metr32(j))
      END DO
    END DO
  END IF
  
END SUBROUTINE ComputeFull_3D

SUBROUTINE MatVecT_3D(y,A,x)

  REAL(RealKind) :: y(:,:,:)
  TYPE(TriTDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:)

  INTEGER :: i,j,k,n1,n2,n3

  n1=A%n1
  n2=A%n2
  n3=A%n3
  CALL MatVecT(y(:,:,1),A%d(1)%P,x(:,:,1))
  IF (n3>1) THEN
    DO j=1,n2
      DO i=1,n1
        y(i,j,1)=y(i,j,1)+A%u(i,j,1)*x(i,j,2)
      END DO   
    END DO   
    DO k=2,n3-1
      CALL MatVecT(y(:,:,k),A%d(k)%P,x(:,:,k))
      DO j=1,n2
        DO i=1,n1
          y(i,j,k)=y(i,j,k)+A%l(i,j,k-1)*x(i,j,k-1)+A%u(i,j,k)*x(i,j,k+1)
        END DO   
      END DO   
    END DO
    CALL MatVecT(y(:,:,n3),A%d(n3)%P,x(:,:,n3))
    DO j=1,n2
      DO i=1,n1
        y(i,j,n3)=y(i,j,n3)+A%l(i,j,n3-1)*x(i,j,n3-1)
      END DO   
    END DO   
  END IF

END SUBROUTINE MatVecT_3D

FUNCTION Res_3D(A,x,b)
  REAL(8) :: Res_3D
  TYPE(TriTDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:)
  REAL(RealKind) :: b(:,:,:)
  REAL(RealKind) :: Res(SIZE(b,1),SIZE(b,2),SIZE(b,3))

  Res=0.0d0
  CALL MatVecT_3D(Res,A,x)
  Res=Res-b
  Res_3D=SUM(ABS(Res))
END FUNCTION Res_3D

SUBROUTINE GaussSeidel_3D(A,b,x,GaussIter)

  TYPE(TriTDiag3D_T) :: A
  REAL(RealKind) :: b(:,:,:)
  REAL(RealKind) :: x(:,:,:)
  INTEGER :: GaussIter

  REAL(RealKind) :: t(A%n1,A%n2)
  INTEGER :: i,j,k,iter

  IF (A%n3>1) THEN
    DO iter=1,GaussIter
      t(:,:)=b(:,:,1)
      DO k=1,A%n3-1
        DO j=1,A%n2
          DO i=1,A%n1
            t(i,j)=t(i,j)-A%u(i,j,k)*x(i,j,k+1)
          END DO   
        END DO   
        CALL Solve(A%d(k)%P,t,x(:,:,k))
        DO j=1,A%n2
          DO i=1,A%n1
            t(i,j)=b(i,j,k+1)-A%l(i,j,k)*x(i,j,k)
          END DO   
        END DO   
      END DO
      CALL Solve(A%d(A%n3)%P,t,x(:,:,A%n3))
      DO j=1,A%n2
        DO i=1,A%n1
          t(i,j)=b(i,j,A%n3-1)-A%u(i,j,A%n3-1)*x(i,j,A%n3)
        END DO   
      END DO   
      DO k=A%n3-1,2,-1
        DO j=1,A%n2
          DO i=1,A%n1
            t(i,j)=t(i,j)-A%l(i,j,k-1)*x(i,j,k-1)
          END DO   
        END DO   
        CALL Solve(A%d(k)%P,t,x(:,:,k))
        DO j=1,A%n2
          DO i=1,A%n1
            t(i,j)=b(i,j,k-1)-A%u(i,j,k-1)*x(i,j,k)
          END DO   
        END DO   
      END DO
    END DO
    CALL Solve(A%d(1)%P,t,x(:,:,1))
!   zusaetzlich 
    t(:,:)=b(:,:,A%n3)-A%l(:,:,A%n3-1)*x(:,:,A%n3-1)
!   t(2,:,:)=b(2,:,:,A%n3)-A%e(2,:,:,A%n3-1)*x(1,:,:,A%n3-1)
    CALL Solve(A%d(A%n3)%P,t,x(:,:,A%n3))
!   zusaetzlich 
  ELSE
    CALL Solve(A%d(1)%P,b(:,:,1),x(:,:,1))
  END IF
   
END SUBROUTINE GaussSeidel_3D

SUBROUTINE Multigrid_3D(A,b,x)

  TYPE(TriTDiag3D_T), TARGET :: A
  REAL(RealKind), TARGET :: b(:,:,:)
  REAL(RealKind), TARGET :: x(:,:,:)

  INTEGER :: i,j,k,n1,n2,n3
  TYPE(TriTDiag3D_T), POINTER :: ALoc
  REAL(RealKind), POINTER :: bLoc(:,:,:)
  REAL(RealKind), POINTER :: xLoc(:,:,:)
  INTEGER :: GaussIter

  bLoc=>b(:,:,:)
  xLoc=>x(:,:,:)
  ALoc=>A
  GaussIter=1
  DO 
    n1=ALoc%n1
    n2=ALoc%n2
    n3=ALoc%n3
!   Relax
    CALL GaussSeidel(ALoc,bLoc,xLoc,GaussIter)
    GaussIter=GaussIter+GaussIter3
    IF (ASSOCIATED(ALoc%Coarse)) THEN
!     RestrictTR
      CALL MatVecT(ALoc%Res,ALoc,xLoc)
      ALoc%Res=bLoc-ALoc%Res
      CALL RestrictT(ALoc%Res,ALoc%ResCoarse)
!     Down
      bLoc=>ALoc%ResCoarse
      xLoc=>ALoc%xCoarse
      xLoc=Zero
      ALoc=>ALoc%Coarse
    ELSE
      GaussIter=GaussIter-GaussIter3
      EXIT
    END IF
  END DO
  DO
!   Smooth
    CALL GaussSeidel(ALoc,bLoc,xLoc,GaussIter)
    GaussIter=GaussIter-GaussIter3
    IF (ASSOCIATED(ALoc%Fine)) THEN
!     Up
      ALoc=>ALoc%Fine
      n1=ALoc%n1
      n2=ALoc%n2
      n3=ALoc%n3
      CALL ProlongateT(xLoc,ALoc%Res)
      IF (ASSOCIATED(ALoc%Fine)) THEN
        xLoc=>Aloc%Fine%xCoarse
        bLoc=>ALoc%Fine%ResCoarse
      ELSE
        xLoc=>x(:,:,:)
        bLoc=>b(:,:,:)
      END IF
      xLoc=xLoc+ALoc%Res
    ELSE
      EXIT
    END IF 
  END DO

END SUBROUTINE Multigrid_3D

SUBROUTINE Solve_3D(A,b,x)

  TYPE(TriTDiag3D_T) :: A
  REAL(RealKind) :: b(:,:,:)
  REAL(RealKind) :: x(:,:,:)

  INTEGER :: Iter

  DO Iter=1,MultIter3
    CALL Multigrid(A,b,x)
  END DO

END SUBROUTINE Solve_3D

SUBROUTINE RestrictT_3D(Res,ResCoarse)

  REAL(RealKind) :: Res(:,:,:)
  REAL(RealKind) :: ResCoarse(:,:,:)

  INTEGER :: n1,n2,n3
  INTEGER :: i1,i2,i3,i1C,i2C,i3C

  n1=SIZE(Res,1)
  n2=SIZE(Res,2)
  n3=SIZE(Res,3)
  ResCoarse=0.0d0
  DO i1=1,n1
    i1C=(i1+Coarse1-1)/Coarse1
    DO i2=1,n2
      i2C=(i2+Coarse2-1)/Coarse2
      DO i3=1,n3
        i3C=(i3+Coarse3-1)/Coarse3
        ResCoarse(i1C,i2C,i3C)=ResCoarse(i1C,i2C,i3C)+Res(i1,i2,i3)
      END DO
    END DO
  END DO
END SUBROUTINE RestrictT_3D

SUBROUTINE ProlongateT_3D(xCoarse,x)
  REAL(RealKind) :: x(:,:,:),xCoarse(:,:,:)

  INTEGER :: n1,n2,n3
  INTEGER :: i1,i2,i3,i1C,i2C,i3C

  n1=SIZE(x,1)
  n2=SIZE(x,2)
  n3=SIZE(x,3)
  DO i1=1,n1
    i1C=(i1+Coarse1-1)/Coarse1
    DO i2=1,n2
      i2C=(i2+Coarse2-1)/Coarse2
      DO i3=1,n3
        i3C=(i3+Coarse3-1)/Coarse3
        x(i1,i2,i3)=xCoarse(i1C,i2C,i3C)
      END DO
    END DO
  END DO

END SUBROUTINE ProlongateT_3D

END MODULE TriT_3D_Mod

