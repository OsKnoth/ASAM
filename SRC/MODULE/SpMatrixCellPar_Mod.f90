MODULE SpMatrixCellPar_Mod
  USE Vector4CellPar_Mod

  IMPLICIT NONE

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_JacMatrixSp
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_MatCellSp
  END INTERFACE
  INTERFACE BICGStabSp
    MODULE PROCEDURE BICGStab_VecCellSp
  END INTERFACE
  INTERFACE SolveSourceSp
    MODULE PROCEDURE SolveSource_VecCellSp
  END INTERFACE
  INTERFACE MatVec
    MODULE PROCEDURE MatVec_VecCellSp
  END INTERFACE
  INTERFACE PSolve
    MODULE PROCEDURE PSolve_VecCellSp
  END INTERFACE
CONTAINS 

SUBROUTINE Value_JacMatrixSp(x,Value)

  IMPLICIT NONE

  TYPE(JacSpMatrix_T), INTENT(OUT)  :: x(:)
  REAL(RealKind), INTENT(IN)  :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    x(ibLoc)%JacT%Val=Value
    x(ibLoc)%JacSLU=Value
  END DO

END SUBROUTINE Value_JacMatrixSp

SUBROUTINE MatVec_VecCellSp(Ax,A,x)

  TYPE(Vector4Cell_T) :: Ax(:)
  TYPE(JacSpMatrix_T) :: A(:)
  TYPE(Vector4Cell_T) :: x(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL MatVec(Ax(ibLoc),A(ibLoc)%JacT,x(ibLoc))
    CALL Axpy(-One,x(ibLoc),Ax(ibLoc))
    CALL PSolve(A(ibLoc)%JacSLU,Ax(ibLoc))
    CALL Xpy(x(ibLoc),Ax(ibLoc))
  END DO
END SUBROUTINE MatVec_VecCellSp

SUBROUTINE PSolve_VecCellSp(Px,P,x)

  TYPE(Vector4Cell_T) :: Px(:)
  TYPE(JacSpMatrix_T) :: P(:)
  TYPE(Vector4Cell_T) :: x(:)

  INTEGER :: ii

  CALL Copy(x,Px)
  CALL ExchangeCell(Px) 
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL PSolve(Px(ibLoc),P(ibLoc)%JacT,Px(ibLoc))
  END DO
  CALL ExchangeCell(Px) 
END SUBROUTINE PSolve_VecCellSp

SUBROUTINE SolveSource_VecCellSp(Px,P)

  TYPE(Vector4Cell_T) :: Px(:)
  TYPE(JacSpMatrix_T) :: P(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL PSolve(P(ibLoc)%JacSLU,Px(ibLoc))
  END DO
END SUBROUTINE SolveSource_VecCellSp

SUBROUTINE Mat1(P,x,y)

  TYPE(Vector4Cell_T) :: x(:),y(:)
  TYPE(JacSpMatrix_T) :: P(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL Mat2(P(ibLoc)%JacSLU,x(ibLoc),y(ibLoc))
  END DO
END SUBROUTINE Mat1

SUBROUTINE  Allocate_MatCellSp(x)

  IMPLICIT NONE

  TYPE(JacSpMatrix_T), POINTER :: x(:)

  IF (.NOT.ASSOCIATED(x)) THEN
    ALLOCATE(x(nbLoc))
    DO ibLoc=1,nb
      x(ibLoc)%JacT=>NULL()
      x(ibLoc)%JacSLU=>NULL()
    END DO
  END IF
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    IF (.NOT.ASSOCIATED(x(ibLoc)%JacT)) THEN
      ALLOCATE(x(ibLoc)%JacT)
      x(ibLoc)%JacT%Val=>NULL()
      x(ibLoc)%JacT%DiagPtr=>NULL()
    END IF
    CALL Allocate(x(ibLoc)%JacT)
    IF (.NOT.ASSOCIATED(x(ibLoc)%JacSLU)) THEN
      ALLOCATE(x(ibLoc)%JacSLU)
      x(ibLoc)%JacSLU%A=>NULL()
      x(ibLoc)%JacSLU%Mat=>NULL()
    END IF
    CALL Allocate(x(ibLoc)%JacSLU)
  END DO

END SUBROUTINE  Allocate_MatCellSp

SUBROUTINE BiCGStab_VecCellSp(A,x,b,MaxIter,Tol)

  TYPE(JacSpMatrix_T) :: A(:)
  TYPE(Vector4Cell_T), POINTER :: b(:),x(:) ! TYPE(Vector) :: b,x
  INTEGER :: MaxIter
  REAL(RealKind) :: Tol


  TYPE(Vector4Cell_T), POINTER :: r(:),rTilde(:),rHat(:),p(:),pHat(:),t(:),v(:) 
  REAL(RealKind) :: alpha,beta,omega,omegaOld,rho,rhoOld
  REAL(RealKind) :: rNorm,rNorm0
  INTEGER :: Iter

  REAL(RealKind) :: Res(MaxIter)

  REAL(RealKind), PARAMETER :: NegOne=-1.0e0
  
  CALL Allocate(r,x)
  CALL Allocate(rTilde,x)
  CALL Allocate(rHat,x)
  CALL Allocate(p,x)
  CALL Allocate(pHat,x)
  CALL Allocate(t,x)
  CALL Allocate(v,x)

! -- Compute initial Residual r=b-A*x --
  CALL MatVec(r,A,x)

! -- r=b-r --
  CALL Xpay(b,NegOne,r)

! -- Choose rTilde, for example rTilde=r --
! CALL Copy(r,rTilde)
  CALL Random(rTilde)

  omegaOld=0.0e0
  rNorm0=SQRT(Dot(r,r))

  DO Iter=1,MaxIter
    rhoOld=rho
    rho=Dot(rTilde,r)
    IF (rho==0.0e0) THEN
      MaxIter=Iter
      Tol=0.0d0
      EXIT
    END IF
    IF (Iter==1) THEN
      CALL Copy(r,p)
    ELSE
!      beta=(rho/rhoOld)*(alpha/omega)
! --  p=r+beta*(p-omega*v) --
      CALL Axpy(-omega,v,p)
      CALL Xpay(r,beta,p)
    END IF
! -- Solve pHat from preconditioning
    CALL PSolve(pHat,A,p)
    CALL MatVec(v,A,pHat)
    alpha=rho/Dot(rTilde,v)
! -- s=r-alpha*v
! -- r=r-alpha*v
    CALL Axpy(-alpha,v,r)

    rNorm=SQRT(Dot(r,r))
    IF (rNorm<=Tol*rNorm0) THEN
      CALL Axpy(alpha,pHat,x)
      MaxIter=-Iter
      Tol=rNorm
      EXIT
    END IF
    CALL PSolve(rHat,A,r)
    CALL MatVec(t,A,rHat)
    omegaOld=omega
    omega=Dot(t,r)/Dot(t,t)
    CALL Axpy(alpha,pHat,x)
    CALL Axpy(omega,rHat,x)
    CALL Axpy(-omega,t,r)
    rNorm=SQRT(Dot(r,r))
    WRITE(*,*) 'BiCG',Iter,rNorm
    IF (rNorm<=Tol*rNorm0) THEN
      MaxIter=Iter
      Tol=rNorm
      EXIT
    END IF
    IF (Iter==MaxIter) THEN
      MaxIter=2*Iter
      Tol=rNorm
      EXIT
    END IF
  END DO

  CALL Deallocate(r)
  CALL Deallocate(rTilde)
  CALL Deallocate(rHat)
  CALL Deallocate(p)
  CALL Deallocate(pHat)
  CALL Deallocate(t)
  CALL Deallocate(v)
 
END SUBROUTINE BiCGStab_VecCellSp
    
END MODULE SpMatrixCellPar_Mod

