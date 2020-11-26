MODULE SpMatrix4CellPar_Mod
  USE Vector4CellPar_Mod
  USE Control_Mod

  IMPLICIT NONE

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Value_JacMatrix4Sp
  END INTERFACE
  INTERFACE BiCGStabChemSp
    MODULE PROCEDURE BiCGStabChem_Vec4CellSp
  END INTERFACE
  INTERFACE BiCGStabMetSp
    MODULE PROCEDURE BiCGStabMet_Vec4CellSp
  END INTERFACE
  INTERFACE SolveSourceChemSp
    MODULE PROCEDURE SolveSourceChem_Vec4CellSp
  END INTERFACE
  INTERFACE SolveSourceMetSp
    MODULE PROCEDURE SolveSourceMet_Vec4CellSp
  END INTERFACE
  INTERFACE MatVecChem
    MODULE PROCEDURE MatVecChem_Vec4CellSp
  END INTERFACE
  INTERFACE MatVecMet
    MODULE PROCEDURE MatVecMet_Vec4CellSp
  END INTERFACE
  INTERFACE PSolve
    MODULE PROCEDURE PSolve_Vec4CellSp
  END INTERFACE
CONTAINS 

SUBROUTINE Value_JacMatrix4Sp(x,Value)

  IMPLICIT NONE

  TYPE(JacSpMatrix4_T), INTENT(INOUT)  :: x(:)
  REAL(RealKind), INTENT(IN)  :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    x(ibLoc)%JacTMom%Val=Value
    x(ibLoc)%JacTPot%Val=Value
    x(ibLoc)%JacFall%Val=Value
    IF (ASSOCIATED(x(ibLoc)%JacFallRhoL%Val)) THEN
      x(ibLoc)%JacFallRhoL%Val=Value
    END IF
    x(ibLoc)%JacSChemLU=Value
    x(ibLoc)%JacSMetLU=Value
  END DO

END SUBROUTINE Value_JacMatrix4Sp

SUBROUTINE MatVecChem_Vec4CellSp(Ax,A,x)

  TYPE(Vector4Cell_T) :: Ax(:)
  TYPE(JacSpMatrix4_T) :: A(:)
  TYPE(Vector4Cell_T) :: x(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL MatVecChem(Ax(ibLoc),A(ibLoc),x(ibLoc))
    CALL Axpy(-One,x(ibLoc),Ax(ibLoc))
    CALL PSolve(A(ibLoc)%JacSChemLU,Ax(ibLoc))
    CALL Xpy(x(ibLoc),Ax(ibLoc)) 
  END DO
END SUBROUTINE MatVecChem_Vec4CellSp

SUBROUTINE MatVecMet_Vec4CellSp(Ax,A,x)

  TYPE(Vector4Cell_T) :: Ax(:)
  TYPE(JacSpMatrix4_T) :: A(:)
  TYPE(Vector4Cell_T) :: x(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL MatVecMet(Ax(ibLoc),A(ibLoc),x(ibLoc))
    CALL Axpy(-One,x(ibLoc),Ax(ibLoc))
    CALL PSolve(A(ibLoc)%JacSMetLU,Ax(ibLoc))
    CALL Xpy(x(ibLoc),Ax(ibLoc)) 
  END DO
END SUBROUTINE MatVecMet_Vec4CellSp

SUBROUTINE PSolve_Vec4CellSp(Px,P,x)

  TYPE(Vector4Cell_T) :: Px(:)
  TYPE(JacSpMatrix4_T) :: P(:)
  TYPE(Vector4Cell_T) :: x(:)

  INTEGER :: ii

  CALL Copy(x,Px)
  CALL ExchangeCell(Px) 
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL PSolve(Px(ibLoc),P(ibLoc),Px(ibLoc))
  END DO
  CALL ExchangeCell(Px) 
END SUBROUTINE PSolve_Vec4CellSp

SUBROUTINE SolveSourceChem_Vec4CellSp(Px,P)

  TYPE(Vector4Cell_T) :: Px(:)
  TYPE(JacSpMatrix4_T) :: P(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL PSolve(P(ibLoc)%JacSChemLU,Px(ibLoc))
  END DO
END SUBROUTINE SolveSourceChem_Vec4CellSp

SUBROUTINE SolveSourceMet_Vec4CellSp(Px,P)

  TYPE(Vector4Cell_T) :: Px(:)
  TYPE(JacSpMatrix4_T) :: P(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL PSolve(P(ibLoc)%JacSMetLU,Px(ibLoc))
  END DO
END SUBROUTINE SolveSourceMet_Vec4CellSp

!SUBROUTINE Mat1(P,x,y)
!
!  TYPE(Vector4Cell_T) :: x(:),y(:)
!  TYPE(JacSpMatrix4_T) :: P(:)
!
!  DO ibLoc=1,nbLoc
!    ib=LocGlob(ibLoc)
!    CALL Set(Floor(ib))
!    CALL Mat2(P(ibLoc)%JacSLU,x(ibLoc),y(ibLoc))
!  END DO
!END SUBROUTINE Mat1

SUBROUTINE BiCGStabChem_Vec4CellSp(A,x,b,MaxIter,Tol)

  TYPE(JacSpMatrix4_T) :: A(:)
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
  CALL MatVecChem(r,A,x)

! -- r=b-r --
  CALL Xpay(b,NegOne,r)

! -- Choose rTilde, for example rTilde=r --
  CALL Copy(r,rTilde)
! CALL Random(rTilde)

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
      beta=(rho/rhoOld)*(alpha/omega)
! --  p=r+beta*(p-omega*v) --
      CALL Axpy(-omega,v,p)
      CALL Xpay(r,beta,p)
    END IF
! -- Solve pHat from preconditioning
    CALL PSolve(pHat,A,p)
    CALL MatVecChem(v,A,pHat)
    alpha=rho/Dot(rTilde,v)
! -- s=r-alpha*v
! -- r=r-alpha*v
    CALL Axpy(-alpha,v,r)

    rNorm=SQRT(Dot(r,r))
    IF (MyId==0.AND.PrintNameLists) THEN
      WRITE(TermUnit,*) 'BiCGS',rNorm,Tol,rNorm0
    END IF
    IF (rNorm<=Tol*rNorm0) THEN
      CALL Axpy(alpha,pHat,x)
      MaxIter=-Iter
      Tol=rNorm
      EXIT
    END IF
    CALL PSolve(rHat,A,r)
    CALL MatVecChem(t,A,rHat)
    omegaOld=omega
    omega=Dot(t,r)/Dot(t,t)
    CALL Axpy(alpha,pHat,x)
    CALL Axpy(omega,rHat,x)
    CALL Axpy(-omega,t,r)
    rNorm=SQRT(Dot(r,r))
    IF (MyId==0.AND.PrintNameLists) THEN
      WRITE(TermUnit,*) 'BiCGS',rNorm,Tol,rNorm0
    END IF
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
 
END SUBROUTINE BiCGStabChem_Vec4CellSp
    
SUBROUTINE BiCGStabMet_Vec4CellSp(A,x,b,MaxIter,Tol)

  TYPE(JacSpMatrix4_T) :: A(:)
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
  CALL MatVecChem(r,A,x)

! -- r=b-r --
  CALL Xpay(b,NegOne,r)

! -- Choose rTilde, for example rTilde=r --
  CALL Copy(r,rTilde)
! CALL Random(rTilde)

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
      beta=(rho/rhoOld)*(alpha/omega)
! --  p=r+beta*(p-omega*v) --
      CALL Axpy(-omega,v,p)
      CALL Xpay(r,beta,p)
    END IF
! -- Solve pHat from preconditioning
    CALL PSolve(pHat,A,p)
    CALL MatVecMet(v,A,pHat)
    alpha=rho/Dot(rTilde,v)
! -- s=r-alpha*v
! -- r=r-alpha*v
    CALL Axpy(-alpha,v,r)

    rNorm=SQRT(Dot(r,r))
    IF (MyId==0.AND.PrintNameLists) THEN
      WRITE(TermUnit,*) 'BiCGS',rNorm,Tol,rNorm0
    END IF
    IF (rNorm<=Tol*rNorm0) THEN
      CALL Axpy(alpha,pHat,x)
      MaxIter=-Iter
      Tol=rNorm
      EXIT
    END IF
    CALL PSolve(rHat,A,r)
    CALL MatVecMet(t,A,rHat)
    omegaOld=omega
    omega=Dot(t,r)/Dot(t,t)
    CALL Axpy(alpha,pHat,x)
    CALL Axpy(omega,rHat,x)
    CALL Axpy(-omega,t,r)
    rNorm=SQRT(Dot(r,r))
    IF (MyId==0.AND.PrintNameLists) THEN
      WRITE(TermUnit,*) 'BiCGS',rNorm,Tol,rNorm0
    END IF
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
 
END SUBROUTINE BiCGStabMet_Vec4CellSp
    
END MODULE SpMatrix4CellPar_Mod

