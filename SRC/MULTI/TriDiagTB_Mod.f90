MODULE TriTB_1D_Mod
  USE Matrix_Mod
  USE Parameter_Mod
  USE Control_Mod
  USE Reverse_Mod

  IMPLICIT NONE

  TYPE TriTBDiag1D_T
    INTEGER :: n1=0
    REAL(RealKind), POINTER :: d(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: u(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: l(:,:,:)=>NULL()
    TYPE(TriTBDiag1D_T), POINTER :: Coarse=>NULL()
    TYPE(TriTBDiag1D_T), POINTER :: Fine=>NULL()
  END TYPE TriTBDiag1D_T

  INTERFACE Allocate
    MODULE PROCEDURE Allocate_1D
  END INTERFACE Allocate
  INTERFACE Solve
    MODULE PROCEDURE Solve_1D
  END INTERFACE Solve
  INTERFACE MatVecTB
    MODULE PROCEDURE MatVecTB_1D
  END INTERFACE MatVecTB
  INTERFACE Compute
    MODULE PROCEDURE Compute_1D
  END INTERFACE Compute
CONTAINS

SUBROUTINE Compute_1D(A,F,FG,Vol,DTU,DUU,DUT)

  TYPE (TriTBDiag1D_T) :: A
  REAL(RealKind) :: F(:),FG(:),Vol(:)
  REAL(RealKind) :: DTU(:),DUU(:),DUT(:)

  INTEGER :: i,n1

! a_1 u_1
! l_1 a_2 u_2
!     l_2 a_3 u_3


  n1=A%n1
  A%d(1,1,1)=A%d(1,1,1)+(beta0*dtp)**2*Fac1_L*F(1)*DTU(1)*DUU(1)*FG(1)*DUT(1)/(Vol(1)+Eps)
  A%d(2,1,1)=A%d(2,1,1)+(beta0*dtp)**2*Fac1_L*F(1)*DUU(1)*FG(1)*DUT(1)/(Vol(1)+Eps)
  A%d(1,2,1)=A%d(1,2,1)-Half*GravComp*(beta0*dtp)**2*Fac1_L*F(1)*DTU(1)*Vol(1)/(Vol(1)+Eps)
  A%d(2,2,1)=A%d(2,2,1)-Half*GravComp*(beta0*dtp)**2*Fac1_L*F(1)*Vol(1)/(Vol(1)+Eps)
  DO i=1,n1-1
    A%u(1,1,i)=-(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*FG(i+1)*DUT(i+1)/(Vol(i)+Vol(i+1)+Eps)
    A%u(2,1,i)=-(beta0*dtp)**2*2.0d0*F(i+1)*DUU(i+1)*FG(i+1)*DUT(i+1)/(Vol(i)+Vol(i+1)+Eps)
    A%l(1,1,i)=-(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*FG(i+1)*DUT(i)/(Vol(i)+Vol(i+1)+Eps)
    A%l(2,1,i)=-(beta0*dtp)**2*2.0d0*F(i+1)*DUU(i+1)*FG(i+1)*DUT(i)/(Vol(i)+Vol(i+1)+Eps)
    IF (.NOT.GradFull) THEN
      A%u(1,2,i)=A%u(1,2,i)-GravComp*(beta0*dtp)**2*F(i+1)*DTU(i+1)*Vol(i+1)/(Vol(i)+Vol(i+1)+Eps)
      A%u(2,2,i)=A%u(2,2,i)-GravComp*(beta0*dtp)**2*F(i+1)*Vol(i+1)/(Vol(i)+Vol(i+1)+Eps)
      A%l(1,2,i)=A%l(1,2,i)+GravComp*(beta0*dtp)**2*F(i+1)*DTU(i+1)*Vol(i)/(Vol(i)+Vol(i+1)+Eps)
      A%l(2,2,i)=A%l(2,2,i)+GravComp*(beta0*dtp)**2*F(i+1)*Vol(i)/(Vol(i)+Vol(i+1)+Eps)
    ELSE
      A%u(1,2,i)=A%u(1,2,i)-GravComp*(beta0*dtp)**2*F(i+1)*DTU(i+1)*Vol(i+1)/(Vol(i)+Vol(i+1)+Eps)
      A%u(2,2,i)=A%u(2,2,i)-GravComp*(beta0*dtp)**2*F(i+1)*Vol(i+1)/(Vol(i)+Vol(i+1)+Eps)
      A%l(1,2,i)=A%l(1,2,i)+GravComp*(beta0*dtp)**2*F(i+1)*DTU(i+1)*Vol(i)/(Vol(i)+Vol(i+1)+Eps)
      A%l(2,2,i)=A%l(2,2,i)+GravComp*(beta0*dtp)**2*F(i+1)*Vol(i)/(Vol(i)+Vol(i+1)+Eps)
    END IF
!   A%d(i)=A%d(i)-A%l(i)
    A%d(1,1,i)=A%d(1,1,i)+(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*FG(i+1)*DUT(i)/(Vol(i)+Vol(i+1)+Eps)
    A%d(2,1,i)=A%d(2,1,i)+(beta0*dtp)**2*2.0d0*F(i+1)*DUU(i+1)*FG(i+1)*DUT(i)/(Vol(i)+Vol(i+1)+Eps)
!   A%d(i+1)=A%d(i+1)-A%u(i)
    A%d(1,1,i+1)=A%d(1,1,i+1)+(beta0*dtp)**2*2.0d0*F(i+1)*DTU(i+1)*DUU(i+1)*FG(i+1)*DUT(i+1)/(Vol(i)+Vol(i+1)+Eps)
    A%d(2,1,i+1)=A%d(2,1,i+1)+(beta0*dtp)**2*2.0d0*F(i+1)*DUU(i+1)*FG(i+1)*DUT(i+1)/(Vol(i)+Vol(i+1)+Eps)
    IF (.NOT.GradFull) THEN
      A%d(1,2,i)=A%d(1,2,i)-GravComp*(beta0*dtp)**2*F(i+1)*DTU(i+1)*Vol(i)/(Vol(i)+Vol(i+1)+Eps)
      A%d(2,2,i)=A%d(2,2,i)-GravComp*(beta0*dtp)**2*F(i+1)*Vol(i)/(Vol(i)+Vol(i+1)+Eps)
      A%d(1,2,i+1)=A%d(1,2,i+1)+GravComp*(beta0*dtp)**2*F(i+1)*DTU(i+1)*Vol(i+1)/(Vol(i)+Vol(i+1)+Eps)
      A%d(2,2,i+1)=A%d(2,2,i+1)+GravComp*(beta0*dtp)**2*F(i+1)*Vol(i+1)/(Vol(i)+Vol(i+1)+Eps)
    ELSE
      A%d(1,2,i)=A%d(1,2,i)-GravComp*(beta0*dtp)**2*F(i+1)*DTU(i+1)*Vol(i)/(Vol(i)+Vol(i+1)+Eps)
      A%d(2,2,i)=A%d(2,2,i)-GravComp*(beta0*dtp)**2*F(i+1)*Vol(i)/(Vol(i)+Vol(i+1)+Eps)
      A%d(1,2,i+1)=A%d(1,2,i+1)+GravComp*(beta0*dtp)**2*F(i+1)*DTU(i+1)*Vol(i+1)/(Vol(i)+Vol(i+1)+Eps)
      A%d(2,2,i+1)=A%d(2,2,i+1)+GravComp*(beta0*dtp)**2*F(i+1)*Vol(i+1)/(Vol(i)+Vol(i+1)+Eps)
    END IF
  END DO
  A%d(1,1,n1)=A%d(1,1,n1)+(beta0*dtp)**2*Fac1_R*F(n1+1)*DTU(n1+1)*DUU(n1+1)*FG(n1+1)*DUT(n1)/(Vol(n1)+Eps)
  A%d(2,1,n1)=A%d(2,1,n1)+(beta0*dtp)**2*Fac1_R*F(n1+1)*DUU(n1+1)*FG(n1+1)*DUT(n1)/(Vol(n1)+Eps)
  A%d(1,2,n1)=A%d(1,2,n1)+Half*GravComp*(beta0*dtp)**2*Fac1_R*F(n1+1)*DTU(n1+1)*Vol(n1)/(Vol(n1)+Eps)
  A%d(2,2,n1)=A%d(2,2,n1)+Half*GravComp*(beta0*dtp)**2*Fac1_R*F(n1+1)*Vol(n1)/(Vol(n1)+Eps)
  DO i=1,n1
    A%d(1,1,i)=A%d(1,1,i)+FacAnela*Vol(i)
    A%d(2,2,i)=A%d(2,2,i)+Vol(i)
  END DO

END SUBROUTINE Compute_1D

SUBROUTINE Allocate_1D(A,n)

  TYPE (TriTBDiag1D_T) :: A
  INTEGER :: i,n

  A%n1=n
  IF (.NOT.ASSOCIATED(A%d)) THEN
    ALLOCATE(A%d(2,2,A%n1))
  END IF
  DO i=1,A%n1
    A%d(1,1,i)=1.d-40
    A%d(2,1,i)=0.0d0
    A%d(1,2,i)=0.0d0
    A%d(2,2,i)=1.d-40
  END DO
  IF (.NOT.ASSOCIATED(A%u)) THEN
    ALLOCATE(A%u(2,2,A%n1-1))
  END IF
  A%u=Zero
  IF (.NOT.ASSOCIATED(A%l)) THEN
    ALLOCATE(A%l(2,2,A%n1-1))
  END IF
  A%l=Zero

END SUBROUTINE Allocate_1D

SUBROUTINE MatVecTB_1D(y,A,x)

  REAL(RealKind) :: y(:,:)
  TYPE(TriTBDiag1D_T) :: A
  REAL(RealKind) :: x(:,:)
 
  REAL(RealKind) :: ss

  INTEGER :: i
  
  IF (A%n1>1) THEN
    y(1:2,1)=Mult_2(A%d(1:2,1:2,1),x(1:2,1))+Mult_2(A%u(1:2,1:2,1),x(1:2,2))
    DO i=2,A%n1-1
      y(1:2,i)=Mult_2(A%l(1:2,1:2,i-1),x(1:2,i-1))+Mult_2(A%d(1:2,1:2,i),x(1:2,i))+Mult_2(A%u(1:2,1:2,i),x(1:2,i+1))
    END DO
    y(1:2,A%n1)=Mult_2(A%l(1:2,1:2,A%n1-1),x(1:2,A%n1-1))+Mult_2(A%d(1:2,1:2,A%n1),x(1:2,A%n1))
  ELSE
    y(1:2,1)=Mult_2(A%d(1:2,1:2,1),x(1:2,1))
  END IF

END SUBROUTINE MatVecTB_1D

SUBROUTINE Solve_1D(A,b)

  TYPE(TriTBDiag1D_T) :: A
  REAL(RealKind) :: b(:,:)

  REAL(RealKind) :: w(1:2,1:2,A%n1)
  REAL(RealKind) :: t(1:2,1:2)
  INTEGER :: i

  t=A%d(1:2,1:2,1)
  b(1:2,1)=Div_2(b(1:2,1),t(1:2,1:2))
  DO i=2,A%n1
    w(1:2,1:2,i-1)=Div_22(A%u(1:2,1:2,i-1),t(1:2,1:2))
    t=A%d(1:2,1:2,i)-Mult_22(A%l(1:2,1:2,i-1),w(1:2,1:2,i-1))
    b(1:2,i)=Div_2(b(1:2,i)-Mult_2(A%l(1:2,1:2,i-1),b(1:2,i-1)),t(1:2,1:2))
  END DO
  DO i=A%n1-1,1,-1
    b(1:2,i)=b(1:2,i)-Mult_2(w(1:2,1:2,i),b(1:2,i+1))
  END DO

END SUBROUTINE Solve_1D

END MODULE TriTB_1D_Mod

MODULE TriTB_2D_Mod

  USE TriTB_1D_Mod

  IMPLICIT NONE

  TYPE TriTBDiag2D_T
    INTEGER :: n1=0
    INTEGER :: n2=0
    TYPE(TriTBDiag1D_T), POINTER :: d(:)=>NULL()
    TYPE(TriTBDiag2D_T), POINTER :: Coarse=>NULL()
    TYPE(TriTBDiag2D_T), POINTER :: Fine=>NULL()
    REAL(RealKind), POINTER :: u(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: l(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: Res(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: ResCoarse(:,:,:)=>NULL()
    REAL(RealKind), POINTER :: xCoarse(:,:,:)=>NULL()
  END TYPE TriTBDiag2D_T
  TYPE TriTBDiag2DP_T
    TYPE(TriTBDiag2D_T), POINTER :: P=>NULL()
  END TYPE TriTBDiag2DP_T


  INTERFACE Allocate
    MODULE PROCEDURE Allocate_2D
  END INTERFACE Allocate
  INTERFACE GaussSeidel
    MODULE PROCEDURE GaussSeidel_2D
  END INTERFACE GaussSeidel
  INTERFACE MatVecTB
    MODULE PROCEDURE MatVecTB_2D
  END INTERFACE MatVecTB
  INTERFACE RestrictTB
    MODULE PROCEDURE RestrictTB_2D
  END INTERFACE RestrictTB
  INTERFACE ProlongateTB
    MODULE PROCEDURE ProlongateTB_2D
  END INTERFACE ProlongateTB
  INTERFACE Compute
    MODULE PROCEDURE Compute_2D
  END INTERFACE Compute
  INTERFACE Multigrid
    MODULE PROCEDURE Multigrid_2D
  END INTERFACE Multigrid
  INTERFACE Solve
    MODULE PROCEDURE Solve_2D
  END INTERFACE Solve


CONTAINS

SUBROUTINE Allocate_2D(A,n2,n1)

  TYPE (TriTBDiag2D_T), POINTER :: A
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
        ALLOCATE(A%u(1:2,n1,n2-1))
      END IF
      A%u=Zero
      IF (.NOT.ASSOCIATED(A%l)) THEN
        ALLOCATE(A%l(1:2,n1,n2-1))
      END IF
      A%l=Zero
      IF (.NOT.ASSOCIATED(A%Res)) THEN
        ALLOCATE(A%Res(1:2,n1,n2))
      END IF
      IF (.NOT.ASSOCIATED(A%ResCoarse)) THEN
        ALLOCATE(A%ResCoarse(1:2,(n1+1)/2,(n2+1)/2))
      END IF
      IF (.NOT.ASSOCIATED(A%xCoarse)) THEN
        ALLOCATE(A%xCoarse(1:2,(n1+1)/2,(n2+1)/2))
      END IF
    END IF
  END IF

END SUBROUTINE Allocate_2D

SUBROUTINE Compute_2D(A,F1,F2,F1G,F2G,Vol,DTU1,DTU2,DUU1,DUU2,DUT)

  TYPE (TriTBDiag2D_T) :: A
  REAL(RealKind) :: F1(:,:)
  REAL(RealKind) :: F2(:,:)
  REAL(RealKind) :: F1G(:,:)
  REAL(RealKind) :: F2G(:,:)
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
    CALL Compute(A%d(1),F1(:,1),F1G(:,1),Vol(:,1),DTU1(:,1),DUU1(:,1),DUT(:,1))
    DO i=1,n1
      A%d(1)%d(1,1,i)=A%d(1)%d(1,1,i)+(beta0*dtp)**2*Fac2_L*F2(i,1)*DTU2(i,1)*DUU2(i,1)*F2G(i,1)*DUT(i,1)/(Vol(i,1)+Eps)
      A%d(1)%d(2,1,i)=A%d(1)%d(2,1,i)+(beta0*dtp)**2*Fac2_L*F2(i,1)*DUU2(i,1)*F2G(i,1)*DUT(i,1)/(Vol(i,1)+Eps)
    END DO
    DO j=1,n2-1
      DO i=1,n1
        A%u(1,i,j)=-(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
        A%u(2,i,j)=-(beta0*dtp)**2*2.0d0*F2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
        A%l(1,i,j)=-(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
        A%l(2,i,j)=-(beta0*dtp)**2*2.0d0*F2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
!       A%d(j)%d(i)=A%d(j)%d(i)-A%e(i,j)
        A%d(j)%d(1,1,i)=A%d(j)%d(1,1,i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
        A%d(j)%d(2,1,i)=A%d(j)%d(2,1,i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
      END DO
      CALL Compute(A%d(j+1),F1(:,j+1),F1G(:,j+1),Vol(:,j+1),DTU1(:,j+1),DUU1(:,j+1),DUT(:,j+1))
      DO i=1,n1
!       A%d(j+1)%d(i)=A%d(j+1)%d(i)-A%e(i,j)
        A%d(j+1)%d(1,1,i)=A%d(j+1)%d(1,1,i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
        A%d(j+1)%d(2,1,i)=A%d(j+1)%d(2,1,i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
      END DO
    END DO
    DO i=1,n1
      A%d(n2)%d(1,1,i)=A%d(n2)%d(1,1,i) &
                         +(beta0*dtp)**2*Fac2_R*F2(i,n2+1)*DTU2(i,n2+1)*DUU2(i,n2+1)*F2G(i,n2+1)*DUT(i,n2)/(Vol(i,n2)+Eps)
      A%d(n2)%d(2,1,i)=A%d(n2)%d(2,1,i) &
                         +(beta0*dtp)**2*Fac2_R*F2(i,n2+1)*DUU2(i,n2+1)*F2G(i,n2+1)*DUT(i,n2)/(Vol(i,n2)+Eps)
    END DO
  ELSE
    CALL Compute(A%d(1),F1(:,1),F1G(:,1),Vol(:,1),DTU1(:,1),DUU1(:,1),DUT(:,1))
    DO i=1,n1
      A%d(1)%d(1,1,i)=A%d(1)%d(1,1,i)+(beta0*dtp)**2*Fac2_L*F2(i,1)*DTU2(i,1)*DUU2(i,1)*F2G(i,1)*DUT(i,1)/(Vol(i,1)+Eps)
      A%d(1)%d(2,1,i)=A%d(1)%d(2,1,i)+(beta0*dtp)**2*Fac2_L*F2(i,1)*DUU2(i,1)*F2G(i,1)*DUT(i,1)/(Vol(i,1)+Eps)
    END DO
    DO j=1,n2-1
      DO i=1,n1
!       A%d(j)%d(i)=A%d(j)%d(i)-e(i)
        A%d(j)%d(1,1,i)=A%d(j)%d(1,1,i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
        A%d(j)%d(2,1,i)=A%d(j)%d(2,1,i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j)/(Vol(i,j)+Vol(i,j+1)+Eps)
      END DO
      CALL Compute(A%d(j+1),F1(:,j+1),F1G(:,j+1),Vol(:,j+1),DTU1(:,j+1),DUU1(:,j+1),DUT(:,j+1))
      DO i=1,n1
!       A%d(j+1)%d(i)=A%d(j+1)%d(i)-e(i)
        A%d(j+1)%d(1,1,i)=A%d(j+1)%d(1,1,i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DTU2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
        A%d(j+1)%d(2,1,i)=A%d(j+1)%d(2,1,i) &
            +(beta0*dtp)**2*2.0d0*F2(i,j+1)*DUU2(i,j+1)*F2G(i,j+1)*DUT(i,j+1)/(Vol(i,j)+Vol(i,j+1)+Eps)
      END DO
    END DO
    DO i=1,n1
      A%d(n2)%d(1,1,i)=A%d(n2)%d(1,1,i) &
                         +(beta0*dtp)**2*Fac2_R*F2(i,n2+1)*DTU2(i,n2+1)*DUU2(i,n2+1)*F2G(i,n2+1)*DUT(i,n2)/(Vol(i,n2)+Eps)
      A%d(n2)%d(2,1,i)=A%d(n2)%d(2,1,i) &
                         +(beta0*dtp)**2*Fac2_R*F2(i,n2+1)*DUU2(i,n2+1)*F2G(i,n2+1)*DUT(i,n2)/(Vol(i,n2)+Eps)
    END DO
  END IF

END SUBROUTINE Compute_2D

SUBROUTINE MatVecTB_2D(y,A,x)

  REAL(RealKind) :: y(:,:,:)
  TYPE(TriTBDiag2D_T) :: A
  REAL(RealKind) :: x(:,:,:)

  INTEGER :: i,j,n1,n2

  n1=A%n1
  n2=A%n2
  CALL MatVecTB(y(:,:,1),A%d(1),x(:,:,1))
  IF (n2>1) THEN
    DO i=1,n1
      y(1:2,i,1)=y(1:2,i,1)+A%u(1:2,i,1)*x(1,i,2)    
    END DO
    DO j=2,n2-1
      CALL MatVecTB(y(:,:,j),A%d(j),x(:,:,j))
      DO i=1,n1
        y(1:2,i,j)=y(1:2,i,j)+A%l(1:2,i,j-1)*x(1,i,j-1)+A%u(1:2,i,j)*x(1,i,j+1)   
      END DO
    END DO
    CALL MatVecTB(y(:,:,n2),A%d(n2),x(:,:,n2))
    DO i=1,n1
      y(1:2,i,n2)=y(1:2,i,n2)+A%l(1:2,i,n2-1)*x(1,i,n2-1)    
    END DO
  END IF

END SUBROUTINE MatVecTB_2D

SUBROUTINE GaussSeidel_2D(A,b,x,GaussIter)

  TYPE(TriTBDiag2D_T) :: A
  REAL(RealKind) :: b(:,:,:)
  REAL(RealKind) :: x(:,:,:)
  INTEGER :: GaussIter

  REAL(RealKind) :: t(1:2,A%n1)
  INTEGER :: i,j,iter

  IF (A%n2>1) THEN
    DO iter=1,GaussIter
      t(:,:)=b(:,:,1)
      DO j=1,A%n2-1
        DO i=1,A%n1
          x(1:2,i,j)=t(1:2,i)-A%u(1:2,i,j)*x(1,i,j+1)
        END DO 
        CALL Solve(A%d(j),x(:,:,j))
        DO i=1,A%n1
          t(1:2,i)=b(1:2,i,j+1)-A%l(1:2,i,j)*x(1,i,j)
        END DO 
      END DO
      x(:,:,A%n2)=t(:,:)
      CALL Solve(A%d(A%n2),x(:,:,A%n2))
      DO i=1,A%n1
        t(1:2,i)=b(1:2,i,A%n2-1)-A%u(1:2,i,A%n2-1)*x(1,i,A%n2)
      END DO  
      DO j=A%n2-1,2,-1
        DO i=1,A%n1
          x(1:2,i,j)=t(1:2,i)-A%l(1:2,i,j-1)*x(1,i,j-1)
        END DO 
        CALL Solve(A%d(j),x(:,:,j))
        DO i=1,A%n1
          t(1:2,i)=b(1:2,i,j-1)-A%u(1:2,i,j-1)*x(1,i,j)
        END DO 
      END DO
    END DO
    x(:,:,1)=t(:,:)
    CALL Solve(A%d(1),x(:,:,1))
  ELSE
    x(:,:,1)=b(:,:,1)
    CALL Solve(A%d(1),x(:,:,1))
  END IF
   
END SUBROUTINE GaussSeidel_2D
SUBROUTINE Solve_2D(A,b,x)

  TYPE(TriTBDiag2D_T) :: A
  REAL(RealKind) :: b(:,:,:)
  REAL(RealKind) :: x(:,:,:)

  INTEGER :: Iter
  
  DO Iter=1,MultIter2
    CALL Multigrid(A,b,x)
  END DO

END SUBROUTINE Solve_2D

SUBROUTINE Multigrid_2D(A,b,x)

  TYPE(TriTBDiag2D_T), TARGET :: A
  REAL(RealKind), TARGET :: b(:,:,:)
  REAL(RealKind), TARGET :: x(:,:,:)

  INTEGER :: i,j,n1,n2
  TYPE(TriTBDiag2D_T), POINTER :: ALoc
  REAL(RealKind), POINTER :: bLoc(:,:,:)
  REAL(RealKind), POINTER :: xLoc(:,:,:)
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
!     RestrictTB
      CALL MatVecTB(ALoc%Res,ALoc,xLoc)
      ALoc%Res=bLoc-ALoc%Res
      CALL RestrictTB(ALoc%Res,ALoc%ResCoarse)
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
      CALL ProlongateTB(xLoc,ALoc%Res)
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
SUBROUTINE RestrictTB_2D(Res,ResCoarse)

  REAL(RealKind) :: Res(:,:,:),ResCoarse(:,:,:)

  INTEGER :: iC,iF,jC,jF,n1,n2,n1C,n2C

  n1=SIZE(Res,2)
  n1C=(n1+1)/2
  n2=SIZE(Res,3)
  n2C=(n2+1)/2
  ResCoarse(:,n1C,n2C)=Res(:,n1,n2)
  iF=1
  DO iC=1,n1/2
    ResCoarse(:,iC,n2C)=Res(:,iF  ,n2) &
                     +Res(:,iF+1,n2) 
    iF=iF+2
  END DO
  jF=1
  DO jC=1,n2/2
    ResCoarse(:,n1C,jC)=Res(:,n1,jF)+Res(:,n1,jF+1)
    iF=1
    DO iC=1,n1/2
      ResCoarse(:,iC,jC)=Res(:,iF  ,jF  ) &
                      +Res(:,iF+1,jF  ) &
                      +Res(:,iF  ,jF+1) &
                      +Res(:,iF+1,jF+1) 
      iF=iF+2
    END DO
    jF=jF+2
  END DO
END SUBROUTINE RestrictTB_2D

SUBROUTINE ProlongateTB_2D(xCoarse,x)

  REAL(RealKind) :: x(:,:,:),xCoarse(:,:,:)

  INTEGER :: iC,iF,jC,jF,n1,n2,n1C,n2C

  n1=SIZE(x,2)
  n1C=(n1+1)/2
  n2=SIZE(x,3)
  n2C=(n2+1)/2

  jF=1
  DO jC=1,n2/2
    iF=1
    DO iC=1,n1/2
      x(:,iF  ,jF  )=xCoarse(:,iC,jC)
      x(:,iF+1,jF  )=xCoarse(:,iC,jC)
      x(:,iF  ,jF+1)=xCoarse(:,iC,jC)
      x(:,iF+1,jF+1)=xCoarse(:,iC,jC)
      iF=iF+2
    END DO
    x(:,n1,jF  )=xCoarse(:,n1C,jC)
    x(:,n1,jF+1)=xCoarse(:,n1C,jC)
    jF=jF+2
  END DO
! CALL ProlongateTB(XCoarse(:,n2C),x(:,n2))
  iF=1
  DO iC=1,n1/2
    x(:,iF  ,n2)=xCoarse(:,iC,n2C)
    x(:,iF+1,n2)=xCoarse(:,iC,n2C)
    iF=iF+2
  END DO
  x(:,n1,n2)=xCoarse(:,n1C,n2C)
END SUBROUTINE ProlongateTB_2D
END MODULE TriTB_2D_Mod

MODULE TriTB_3D_Mod

  USE TriTB_2D_Mod

  IMPLICIT NONE

  TYPE TriTBDiag3D_T
    INTEGER :: n1=0
    INTEGER :: n2=0
    INTEGER :: n3=0
    TYPE(TriTBDiag2DP_T), POINTER :: d(:)=>NULL()
    REAL(RealKind), POINTER :: u(:,:,:,:)=>NULL()
    REAL(RealKind), POINTER :: l(:,:,:,:)=>NULL()
    TYPE(TriTBDiag3D_T), POINTER :: Coarse=>NULL()
    TYPE(TriTBDiag3D_T), POINTER :: Fine=>NULL()
    REAL(RealKind), POINTER :: Res(:,:,:,:)=>NULL()
    REAL(RealKind), POINTER :: ResCoarse(:,:,:,:)=>NULL()
    REAL(RealKind), POINTER :: xCoarse(:,:,:,:)=>NULL()
  END TYPE TriTBDiag3D_T



  INTERFACE Allocate
    MODULE PROCEDURE Allocate_3D
  END INTERFACE Allocate
  INTERFACE GaussSeidel
    MODULE PROCEDURE GaussSeidel_3D
  END INTERFACE GaussSeidel
  INTERFACE MatVecTB
    MODULE PROCEDURE MatVecTB_3D
  END INTERFACE MatVecTB
  INTERFACE RestrictTB
    MODULE PROCEDURE RestrictTB_3D
  END INTERFACE RestrictTB
  INTERFACE ProlongateTB
    MODULE PROCEDURE ProlongateTB_3D
  END INTERFACE ProlongateTB
  INTERFACE Compute
    MODULE PROCEDURE Compute_3D
  END INTERFACE Compute
  INTERFACE Multigrid
    MODULE PROCEDURE Multigrid_3D
  END INTERFACE Multigrid
  INTERFACE Solve
    MODULE PROCEDURE Solve_3D
  END INTERFACE Solve


CONTAINS

SUBROUTINE Allocate_3D(A,n3,n2,n1)

  TYPE(TriTBDiag3D_T), POINTER :: A
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
          ALLOCATE(A%u(1:2,n1,n2,n3-1))
        END IF
        A%u=Zero
        IF (.NOT.ASSOCIATED(A%l)) THEN
          ALLOCATE(A%l(1:2,n1,n2,n3-1))
        END IF
        A%l=Zero
        IF (.NOT.ASSOCIATED(A%Res)) THEN
          ALLOCATE(A%Res(1:2,n1,n2,n3))
        END IF
        IF (.NOT.ASSOCIATED(A%ResCoarse)) THEN
          ALLOCATE(A%ResCoarse(1:2,(n1+1)/2,(n2+1)/2,(n3+1)/2))
        END IF
        IF (.NOT.ASSOCIATED(A%xCoarse)) THEN
          ALLOCATE(A%xCoarse(1:2,(n1+1)/2,(n2+1)/2,(n3+1)/2))
        END IF
      END IF
    END IF
  END IF

END SUBROUTINE Allocate_3D

SUBROUTINE Compute_3D(A,F1,F2,F3,F1G,F2G,F3G,Vol,DTU1,DTU2,DTU3,DUU1,DUU2,DUU3,DUT)

  TYPE (TriTBDiag3D_T) :: A
  REAL(RealKind) :: F1(:,:,:)
  REAL(RealKind) :: F2(:,:,:)
  REAL(RealKind) :: F3(:,:,:)
  REAL(RealKind) :: F1G(:,:,:)
  REAL(RealKind) :: F2G(:,:,:)
  REAL(RealKind) :: F3G(:,:,:)
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
    CALL Compute(A%d(1)%P,F1(:,:,1),F2(:,:,1),F1G(:,:,1),F2G(:,:,1),Vol(:,:,1), &
                 DTU1(:,:,1),DTU2(:,:,1),DUU1(:,:,1),DUU2(:,:,1),DUT(:,:,1))
    DO j=1,n2
      DO i=1,n1
        A%d(1)%P%d(j)%d(1,1,i)=A%d(1)%P%d(j)%d(1,1,i) &
           +(beta0*dtp)**2*Fac3_L*F3(i,j,1)*DTU3(i,j,1)*DUU3(i,j,1)*F3G(i,j,1)*DUT(i,j,1)/(Vol(i,j,1)+Eps)
        A%d(1)%P%d(j)%d(2,1,i)=A%d(1)%P%d(j)%d(2,1,i) &
           +(beta0*dtp)**2*Fac3_L*F3(i,j,1)*DUU3(i,j,1)*F3G(i,j,1)*DUT(i,j,1)/(Vol(i,j,1)+Eps)
      END DO
    END DO
    DO k=1,n3-1
      DO j=1,n2
        DO i=1,n1
          A%u(1,i,j,k)=-(beta0*dtp)**2*2.0d0*F3(i,j,k+1) &
              *DTU3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
          A%u(2,i,j,k)=-(beta0*dtp)**2*2.0d0*F3(i,j,k+1) &
              *F3G(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
!         A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k)%P%d(j)%d(1,1,i)=A%d(k)%P%d(j)%d(1,1,i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
          A%d(k)%P%d(j)%d(2,1,i)=A%d(k)%P%d(j)%d(2,1,i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
        END DO
      END DO
      CALL Compute(A%d(k+1)%P,F1(:,:,k+1),F2(:,:,k+1),F1G(:,:,k+1),F2G(:,:,k+1),Vol(:,:,k+1), &
                   DTU1(:,:,k+1),DTU2(:,:,k+1),DUU1(:,:,k+1),DUU2(:,:,k+1),DUT(:,:,k+1))
      DO j=1,n2
        DO i=1,n1
          A%l(1,i,j,k)=-(beta0*dtp)**2*2.0d0*F3(i,j,k+1) &
              *DTU3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
          A%l(2,i,j,k)=-(beta0*dtp)**2*2.0d0*F3(i,j,k+1) &
              *F3G(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
!         A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k+1)%P%d(j)%d(1,1,i)=A%d(k+1)%P%d(j)%d(1,1,i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
          A%d(k+1)%P%d(j)%d(2,1,i)=A%d(k+1)%P%d(j)%d(2,1,i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
        END DO
      END DO
    END DO
    DO j=1,n2
      DO i=1,n1
        A%d(n3)%P%d(j)%d(1,1,i)=A%d(n3)%P%d(j)%d(1,1,i) &
           +(beta0*dtp)**2*Fac3_R*F3(i,j,n3+1)*DTU3(i,j,n3+1)*F3G(i,j,n3+1)*DUT(i,j,n3)/(Vol(i,j,n3)+Eps)
        A%d(n3)%P%d(j)%d(2,1,i)=A%d(n3)%P%d(j)%d(2,1,i) &
           +(beta0*dtp)**2*Fac3_R*F3(i,j,n3+1)*F3G(i,j,n3+1)*DUT(i,j,n3)/(Vol(i,j,n3)+Eps)
      END DO
    END DO
  ELSE
    n2=A%d(1)%P%n2
    n1=A%d(1)%P%n1
    CALL Compute(A%d(1)%P,F1(:,:,1),F2(:,:,1),F1G(:,:,1),F2G(:,:,1),Vol(:,:,1), &
                 DTU1(:,:,1),DTU2(:,:,1),DUU1(:,:,1),DUU2(:,:,1),DUT(:,:,1))
    DO j=1,n2
      DO i=1,n1
        A%d(1)%P%d(j)%d(1,1,i)=A%d(1)%P%d(j)%d(1,1,i) &
           +(beta0*dtp)**2*Fac3_L*F3(i,j,1)*DTU3(i,j,1)*DUU3(i,j,1)*F3G(i,j,1)*DUT(i,j,1)/(Vol(i,j,1)+Eps)
        A%d(1)%P%d(j)%d(2,1,i)=A%d(1)%P%d(j)%d(2,1,i) &
           +(beta0*dtp)**2*Fac3_L*F3(i,j,1)*DUU3(i,j,1)*F3G(i,j,1)*DUT(i,j,1)/(Vol(i,j,1)+Eps)
      END DO
    END DO
    DO k=1,n3-1
      DO j=1,n2
        DO i=1,n1
!         A%d(k)%P%d(j)%d(i)=A%d(k)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k)%P%d(j)%d(1,1,i)=A%d(k)%P%d(j)%d(1,1,i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*DUU3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
          A%d(k)%P%d(j)%d(2,1,i)=A%d(k)%P%d(j)%d(2,1,i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DUU3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
        END DO
      END DO
      CALL Compute(A%d(k+1)%P,F1(:,:,k+1),F2(:,:,k+1),F1G(:,:,k+1),F2G(:,:,k+1),Vol(:,:,k+1), &
                   DTU1(:,:,k+1),DTU2(:,:,k+1),DUU1(:,:,k+1),DUU2(:,:,k+1),DUT(:,:,k+1))
      DO j=1,n2
        DO i=1,n1
!         A%d(k+1)%P%d(j)%d(i)=A%d(k+1)%P%d(j)%d(i)-A%e(:,i,j,k)
          A%d(k+1)%P%d(j)%d(1,1,i)=A%d(k+1)%P%d(j)%d(1,1,i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DTU3(i,j,k+1)*DUU3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
          A%d(k+1)%P%d(j)%d(2,1,i)=A%d(k+1)%P%d(j)%d(2,1,i) &
           +(beta0*dtp)**2*2.0d0*F3(i,j,k+1)*DUU3(i,j,k+1)*F3G(i,j,k+1)*DUT(i,j,k+1)/(Vol(i,j,k)+Vol(i,j,k+1)+Eps)
        END DO
      END DO
    END DO
    DO j=1,n2
      DO i=1,n1
        A%d(n3)%P%d(j)%d(1,1,i)=A%d(n3)%P%d(j)%d(1,1,i) &
           +(beta0*dtp)**2*Fac3_R*F3(i,j,n3+1)*DTU3(i,j,n3+1)*DUU3(i,j,n3+1)*F3G(i,j,n3+1)*DUT(i,j,n3)/(Vol(i,j,n3)+Eps)
        A%d(n3)%P%d(j)%d(2,1,i)=A%d(n3)%P%d(j)%d(2,1,i) &
           +(beta0*dtp)**2*Fac3_R*F3(i,j,n3+1)*DUU3(i,j,n3+1)*F3G(i,j,n3+1)*DUT(i,j,n3)/(Vol(i,j,n3)+Eps)
      END DO
    END DO
  END IF
  
END SUBROUTINE Compute_3D

SUBROUTINE MatVecTB_3D(y,A,x)

  REAL(RealKind) :: y(:,:,:,:)
  TYPE(TriTBDiag3D_T) :: A
  REAL(RealKind) :: x(:,:,:,:)

  INTEGER :: i,j,k,n1,n2,n3

  n1=A%n1
  n2=A%n2
  n3=A%n3
  CALL MatVecTB(y(:,:,:,1),A%d(1)%P,x(:,:,:,1))
  IF (n3>1) THEN
    DO j=1,n2
      DO i=1,n1
        y(1:2,i,j,1)=y(1:2,i,j,1)+A%u(1:2,i,j,1)*x(1,i,j,2)
      END DO   
    END DO   
    DO k=2,n3-1
      CALL MatVecTB(y(:,:,:,k),A%d(k)%P,x(:,:,:,k))
      DO j=1,n2
        DO i=1,n1
          y(1:2,i,j,k)=y(1:2,i,j,k)+A%l(1:2,i,j,k-1)*x(1,i,j,k-1)+A%u(1:2,i,j,k)*x(1,i,j,k+1)
        END DO   
      END DO   
    END DO
    CALL MatVecTB(y(:,:,:,n3),A%d(n3)%P,x(:,:,:,n3))
    DO j=1,n2
      DO i=1,n1
        y(1:2,i,j,n3)=y(1:2,i,j,n3)+A%l(1:2,i,j,n3-1)*x(1,i,j,n3-1)
      END DO   
    END DO   
  END IF

END SUBROUTINE MatVecTB_3D

SUBROUTINE GaussSeidel_3D(A,b,x,GaussIter)

  TYPE(TriTBDiag3D_T) :: A
  REAL(RealKind) :: b(:,:,:,:)
  REAL(RealKind) :: x(:,:,:,:)
  INTEGER :: GaussIter

  REAL(RealKind) :: t(2,A%n1,A%n2)
  INTEGER :: i,j,k,iter

  IF (A%n3>1) THEN
    DO iter=1,GaussIter
      t(:,:,:)=b(:,:,:,1)
      DO k=1,A%n3-1
        DO j=1,A%n2
          DO i=1,A%n1
            t(1:2,i,j)=t(1:2,i,j)-A%u(1:2,i,j,k)*x(1,i,j,k+1)
          END DO   
        END DO   
        CALL Solve(A%d(k)%P,t,x(:,:,:,k))
        DO j=1,A%n2
          DO i=1,A%n1
            t(1:2,i,j)=b(1:2,i,j,k+1)-A%l(1:2,i,j,k)*x(1,i,j,k)
          END DO   
        END DO   
      END DO
      CALL Solve(A%d(A%n3)%P,t,x(:,:,:,A%n3))
      DO j=1,A%n2
        DO i=1,A%n1
          t(1:2,i,j)=b(1:2,i,j,A%n3-1)-A%u(1:2,i,j,A%n3-1)*x(1,i,j,A%n3)
        END DO   
      END DO   
      DO k=A%n3-1,2,-1
        DO j=1,A%n2
          DO i=1,A%n1
            t(1:2,i,j)=t(1:2,i,j)-A%l(1:2,i,j,k-1)*x(1,i,j,k-1)
          END DO   
        END DO   
        CALL Solve(A%d(k)%P,t,x(:,:,:,k))
        DO j=1,A%n2
          DO i=1,A%n1
            t(1:2,i,j)=b(1:2,i,j,k-1)-A%u(1:2,i,j,k-1)*x(1,i,j,k)
          END DO   
        END DO   
      END DO
    END DO
    CALL Solve(A%d(1)%P,t,x(:,:,:,1))
!   zusaetzlich 
    t(1,:,:)=b(1,:,:,A%n3)-A%l(1,:,:,A%n3-1)*x(1,:,:,A%n3-1)
    t(2,:,:)=b(2,:,:,A%n3)-A%l(2,:,:,A%n3-1)*x(1,:,:,A%n3-1)
    CALL Solve(A%d(A%n3)%P,t,x(:,:,:,A%n3))
!   zusaetzlich 
  ELSE
    CALL Solve(A%d(1)%P,b(:,:,:,1),x(:,:,:,1))
  END IF
   
END SUBROUTINE GaussSeidel_3D

SUBROUTINE Multigrid_3D(A,b,x)

  TYPE(TriTBDiag3D_T), TARGET :: A
  REAL(RealKind), TARGET :: b(:,:,:,:)
  REAL(RealKind), TARGET :: x(:,:,:,:)

  INTEGER :: i,j,k,n1,n2,n3
  TYPE(TriTBDiag3D_T), POINTER :: ALoc
  REAL(RealKind), POINTER :: bLoc(:,:,:,:)
  REAL(RealKind), POINTER :: xLoc(:,:,:,:)
  INTEGER :: GaussIter

  bLoc=>b(:,:,:,:)
  xLoc=>x(:,:,:,:)
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
!     RestrictTB
      CALL MatVecTB(ALoc%Res,ALoc,xLoc)
      ALoc%Res=bLoc-ALoc%Res
      CALL RestrictTB(ALoc%Res,ALoc%ResCoarse)
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
      CALL ProlongateTB(xLoc,ALoc%Res)
      IF (ASSOCIATED(ALoc%Fine)) THEN
        xLoc=>Aloc%Fine%xCoarse
        bLoc=>ALoc%Fine%ResCoarse
      ELSE
        xLoc=>x(:,:,:,:)
        bLoc=>b(:,:,:,:)
      END IF
      xLoc=xLoc+ALoc%Res
    ELSE
      EXIT
    END IF 
  END DO

END SUBROUTINE Multigrid_3D

SUBROUTINE Solve_3D(A,b,x)

  TYPE(TriTBDiag3D_T) :: A
  REAL(RealKind) :: b(:,:,:,:)
  REAL(RealKind) :: x(:,:,:,:)

  INTEGER :: Iter

  DO Iter=1,MultIter3
    CALL Multigrid(A,b,x)
  END DO

END SUBROUTINE Solve_3D

SUBROUTINE RestrictTB_3D(Res,ResCoarse)

  REAL(RealKind) :: Res(:,:,:,:),ResCoarse(:,:,:,:)

  INTEGER :: iF,jF,kF,n1,n2,n3
  INTEGER :: iC,jC,kC,n1C,n2C,n3C

  n1=SIZE(Res,2)
  n1C=(n1+1)/2
  n2=SIZE(Res,3)
  n2C=(n2+1)/2
  n3=SIZE(Res,4)
  n3C=(n3+1)/2


  ResCoarse(:,n1C,n2C,n3C)=Res(:,n1,n2,n3)
  iF=1
  DO iC=1,n1/2
    ResCoarse(:,iC,n2C,n3C)=Res(:,iF  ,n2,n3) &
                         +Res(:,iF+1,n2,n3)  
    iF=iF+2
  END DO

  jF=1
  DO jC=1,n2/2
    ResCoarse(:,n1C,jC,n3C)=Res(:,n1,jF  ,n3) &
                         +Res(:,n1,jF+1,n3) 
    iF=1
    DO iC=1,n1/2
      ResCoarse(:,iC,jC,n3C)=Res(:,iF  ,jF  ,n3) &
                          +Res(:,iF+1,jF  ,n3) &
                          +Res(:,iF  ,jF+1,n3) &
                          +Res(:,iF+1,jF+1,n3)  
      iF=iF+2
    END DO
    jF=jF+2
  END DO

  kF=1
  DO kC=1,n3/2
    ResCoarse(:,n1C,n2C,kC)=Res(:,n1,n2,kF  ) &
                         +Res(:,n1,n2,kF+1)  
    iF=1
    DO iC=1,n1/2
      ResCoarse(:,iC,n2C,kC)=Res(:,iF  ,n2,kF  ) &
                          +Res(:,iF+1,n2,kF  ) &
                          +Res(:,iF  ,n2,kF+1) &
                          +Res(:,iF+1,n2,kF+1)  
      iF=iF+2
    END DO
    jF=1
    DO jC=1,n2/2
      ResCoarse(:,n1C,jC,kC)=Res(:,n1,jF  ,kF  ) &
                          +Res(:,n1,jF+1,kF  ) &
                          +Res(:,n1,jF  ,kF+1) &
                          +Res(:,n1,jF+1,kF+1)  
      iF=1
      DO iC=1,n1/2
        ResCoarse(:,iC,jC,kC)=Res(:,iF  ,jF  ,kF  ) &
                           +Res(:,iF+1,jF  ,kF  ) &
                           +Res(:,iF  ,jF+1,kF  ) &
                           +Res(:,iF  ,jF  ,kF+1) &
                           +Res(:,iF+1,jF+1,kF  ) &
                           +Res(:,iF+1,jF  ,kF+1) &
                           +Res(:,iF  ,jF+1,kF+1) &
                           +Res(:,iF+1,jF+1,kF+1)
        iF=iF+2
      END DO
      jF=jF+2
    END DO
    kF=kF+2
  END DO
END SUBROUTINE RestrictTB_3D

SUBROUTINE ProlongateTB_3D(xCoarse,x)

  REAL(RealKind) :: x(:,:,:,:),xCoarse(:,:,:,:)

  INTEGER :: iF,jF,kF,n1,n2,n3
  INTEGER :: iC,jC,kC,n1C,n2C,n3C

  n1=SIZE(x,2)
  n1C=(n1+1)/2 
  n2=SIZE(x,3)
  n2C=(n2+1)/2 
  n3=SIZE(x,4)
  n3C=(n3+1)/2 

  kF=1
  DO kC=1,n3/2
    jF=1
    DO jC=1,n2/2
      iF=1
      DO iC=1,n1/2
        x(:,iF  ,jF  ,kF  )=xCoarse(:,iC,jC,kC)
        x(:,iF+1,jF  ,kF  )=xCoarse(:,iC,jC,kC)
        x(:,iF  ,jF+1,kF  )=xCoarse(:,iC,jC,kC)
        x(:,iF  ,jF  ,kF+1)=xCoarse(:,iC,jC,kC)
        x(:,iF+1,jF+1,kF  )=xCoarse(:,iC,jC,kC)
        x(:,iF+1,jF  ,kF+1)=xCoarse(:,iC,jC,kC)
        x(:,iF  ,jF+1,kF+1)=xCoarse(:,iC,jC,kC)
        x(:,iF+1,jF+1,kF+1)=xCoarse(:,iC,jC,kC)
        iF=iF+2
      END DO
      x(:,n1,jF  ,kF  )=xCoarse(:,n1C,jC,kC)
      x(:,n1,jF+1,kF  )=xCoarse(:,n1C,jC,kC)
      x(:,n1,jF  ,kF+1)=xCoarse(:,n1C,jC,kC)
      x(:,n1,jF+1,kF+1)=xCoarse(:,n1C,jC,kC)
      jF=jF+2
    END DO
    iF=1
    DO iC=1,n1/2
      x(:,iF  ,n2,kF  )=xCoarse(:,iC,n2C,kC)
      x(:,iF+1,n2,kF  )=xCoarse(:,iC,n2C,kC)
      x(:,iF  ,n2,kF+1)=xCoarse(:,iC,n2C,kC)
      x(:,iF+1,n2,kF+1)=xCoarse(:,iC,n2C,kC)
      iF=iF+2
    END DO
    x(:,n1,n2,kF  )=xCoarse(:,n1C,n2C,kC)
    x(:,n1,n2,kF+1)=xCoarse(:,n1C,n2C,kC)
    kF=kF+2
  END DO

  jF=1
  DO jC=1,n2/2
    iF=1
    DO iC=1,n1/2
      x(:,iF  ,jF  ,n3)=xCoarse(:,iC,jC,n3C)
      x(:,iF+1,jF  ,n3)=xCoarse(:,iC,jC,n3C)
      x(:,iF  ,jF+1,n3)=xCoarse(:,iC,jC,n3C)
      x(:,iF+1,jF+1,n3)=xCoarse(:,iC,jC,n3C)
      iF=iF+2
    END DO
    x(:,n1,jF  ,n3)=xCoarse(:,n1C,jC,n3C)
    x(:,n1,jF+1,n3)=xCoarse(:,n1C,jC,n3C)
    x(:,n1,jF  ,n3)=xCoarse(:,n1C,jC,n3C)
    x(:,n1,jF+1,n3)=xCoarse(:,n1C,jC,n3C)
    jF=jF+2
  END DO
  iF=1
  DO iC=1,n1/2
    x(:,iF  ,n2,n3)=xCoarse(:,iC,n2C,n3C)
    x(:,iF+1,n2,n3)=xCoarse(:,iC,n2C,n3C)
    iF=iF+2
  END DO
  x(:,n1,n2,n3)=xCoarse(:,n1C,n2C,n3C)
END SUBROUTINE ProlongateTB_3D
END MODULE TriTB_3D_Mod

