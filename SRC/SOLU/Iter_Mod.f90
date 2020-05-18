MODULE Iter_Mod
  USE VecPrVel_Mod
  USE MatVec_Mod
  USE PreCond_Mod
  USE Control_Mod
  IMPLICIT NONE
CONTAINS

SUBROUTINE SolveSound(x,b,MaxIter,Tol)
  TYPE(PressureVelocity), DIMENSION(:) :: b,x
  INTEGER :: MaxIter
  REAL(RealKind) :: Tol
  
  SELECT CASE(MethSound)
    CASE('GMRES')
      CALL GMRES(x,b,MaxIter,Tol)
    CASE('FGMRES')
      CALL FGMRES(x,b,MaxIter,Tol)
    CASE('FGMRESV')
      CALL FGMRESV(x,b,MaxIter,Tol)
    CASE('QMR')
      CALL QMR(x,b,MaxIter,Tol)
    CASE('MINRES')
      CALL MinRes(x,b,MaxIter,Tol)
    CASE('CGR')
      CALL CGR(x,b,MaxIter,Tol)
    CASE DEFAULT
      CALL GMRES(x,b,MaxIter,Tol)
  END SELECT
END SUBROUTINE SolveSound

SUBROUTINE GMRES(x,b,MaxIter,Tol)
  TYPE(PressureVelocity), DIMENSION(:) :: b,x
  INTEGER :: MaxIter
  REAL(RealKind) :: Tol
  TYPE(PressureVelocity), POINTER, DIMENSION(:) ::  &
       w
  TYPE(PressureVelocity), POINTER, DIMENSION(:) ::  &
       w1
  TYPE(PressureVelocity), POINTER, DIMENSION(:,:) ::  &
       v
  REAL(RealKind) :: beta,delta,c(StepGMRES),s(StepGMRES),Rho,Rho0,Temp
  REAL(RealKind) :: H(StepGMRES+1,StepGMRES)
  REAL(RealKind) :: bHat(StepGMRES+1)
  REAL(RealKind) :: y(StepGMRES)
  INTEGER :: Iter,i,k,nR
  REAL(RealKind) :: Temp1
  CALL Allocate(v,StepGMRES+1)
  CALL Allocate(w)
  CALL Allocate(w1)
! Compute initial Residual w=b-A*x
  DO Iter=1,MaxIter
    CALL MatVec1(w,x)
    CALL DXPBY(b,-One,w)
    CALL Exchange(w) !OSSI
    CALL PSolve1(v(:,1),w)
!   beta=SQRT(v(:,1)*v(:,1))
    beta=SQRT(DOT2(v(:,1),v(:,1)))
    IF (Iter==1) THEN
      Rho0=beta
    END IF
    IF (beta<=Zero) THEN
      EXIT
    END IF
    CALL DSCALE(One/beta,v(:,1)) !
    bHat(1)=beta*One
    S2:DO i=1,StepGMRES
      CALL MatVec1(w,v(:,i))
      CALL Exchange(w) !OSSI
      CALL PSolve1(v(:,i+1),w)
      DO k=1,i
!       H(k,i)=v(:,k)*v(:,i+1)
        H(k,i)=DOT2(v(:,k),v(:,i+1))
        CALL DAXPY(-H(k,i),v(:,k),v(:,i+1))
      END DO
!     H(i+1,i)=SQRT(v(:,i+1)*v(:,i+1))
      H(i+1,i)=SQRT(DOT2(v(:,i+1),v(:,i+1)))
      CALL DSCALE(One/H(i+1,i),v(:,i+1))
      DO k=2,i
        Temp=c(k-1)*H(k-1,i)+s(k-1)*H(k,i)
        H(k,i)=-s(k-1)*H(k-1,i)+c(k-1)*H(k,i)
        H(k-1,i)=Temp
      END DO
      delta=SQRT(H(i,i)*H(i,i)+H(i+1,i)*H(i+1,i))
      c(i)=H(i,i)/delta
      s(i)=H(i+1,i)/delta
      H(i,i)=c(i)*H(i,i)+s(i)*H(i+1,i)
      bHat(i+1)=-s(i)*bHat(i)
      bHat(i)=c(i)*bHat(i)
      rho=ABS(bHat(i+1))
      IF (MyId==0.AND.PrintNameLists) THEN
        WRITE(*,*) 'GMRES Rho',i,rho,rho0,Tol,Tol*rho0
      END IF
      IF (rho<=Tol*Rho0) THEN
        nR=i
        EXIT S2
      END IF
      nr=StepGMRES
    END DO S2
    y(nR)=bHat(nR)/H(nR,nR)
    DO k=nR-1,1,-1
      y(k)=bHat(k)
      DO i=k+1,nR
        y(k)=y(k)-H(k,i)*y(i)
      END DO
      y(k)=y(k)/H(k,k)
    END DO
    DO i=1,nR 
      CALL DAXPY(y(i),v(:,i),x)
    END DO
    IF (Rho<=Tol*Rho0) THEN
      EXIT
    END IF
  END DO
  CALL Deallocate(w)
  CALL Deallocate(w1)
  CALL Deallocate(v)
    
END SUBROUTINE GMRES

SUBROUTINE FGMRES(x,b,MaxIter,Tol)
  TYPE(PressureVelocity), DIMENSION(:) :: b,x
  INTEGER :: MaxIter
  REAL(RealKind) :: Tol
  TYPE(PressureVelocity), POINTER, DIMENSION(:,:) ::  &
       v,z
  REAL(RealKind) :: beta,delta,c(StepGMRES),s(StepGMRES),rho,rho0,temp
  REAL(RealKind) :: H(StepGMRES+1,StepGMRES)
  REAL(RealKind) :: bHat(StepGMRES+1)
  REAL(RealKind) :: y(StepGMRES)
  INTEGER :: Iter,i,k,nR
  INTEGER :: ii,ib
  CALL Allocate(v,StepGMRES+1)
  CALL Allocate(z,StepGMRES)
! Compute initial Residual w=b-A*x
  DO Iter=1,MaxIter
    CALL MatVec1(v(:,1),x)
    CALL DXPBY(b,-One,v(:,1))
    CALL Exchange(v(:,1))
    beta=SQRT(DOT2(v(:,1),v(:,1)))
    IF (Iter==1) THEN
      rho0=beta
    END IF
    IF (beta<=Zero) THEN
      EXIT
    END IF
    CALL DSCALE(One/beta,v(:,1))
    bHat(1)=beta*One
    S2:DO i=1,StepGMRES
      CALL PSolve1(z(:,i),v(:,i))
      CALL MatVec1(v(:,i+1),z(:,i))
      CALL Exchange(v(:,i+1))
      DO k=1,i
        H(k,i)=DOT2(v(:,k),v(:,i+1))
        CALL DAXPY(-H(k,i),v(:,k),v(:,i+1))
      END DO
      H(i+1,i)=SQRT(DOT2(v(:,i+1),v(:,i+1)))
      CALL DSCALE(One/H(i+1,i),v(:,i+1))
      DO k=2,i
        Temp=c(k-1)*H(k-1,i)+s(k-1)*H(k,i)
        H(k,i)=-s(k-1)*H(k-1,i)+c(k-1)*H(k,i)
        H(k-1,i)=Temp 
      END DO
      delta=SQRT(H(i,i)*H(i,i)+H(i+1,i)*H(i+1,i))
      c(i)=H(i,i)/delta
      s(i)=H(i+1,i)/delta
      H(i,i)=c(i)*H(i,i)+s(i)*H(i+1,i)
      bHat(i+1)=-s(i)*bHat(i)
      bHat(i)=c(i)*bHat(i)
      rho=ABS(bHat(i+1))
      IF (MyId==0.AND.PrintNameLists) THEN
        WRITE(*,*) 'FGMRES Rho',i,rho,rho0,Tol,Tol*rho0
      END IF
      IF (rho<=Tol*rho0) THEN
!     IF (rho<=Tol) THEN
        nR=i
        EXIT S2
      END IF
      nr=StepGMRES
    END DO S2
    y(nR)=bHat(nR)/H(nR,nR)
    DO k=nR-1,1,-1
      y(k)=bHat(k)
      DO i=k+1,nR
        y(k)=y(k)-H(k,i)*y(i)
      END DO
      y(k)=y(k)/H(k,k)
    END DO
    DO i=1,nR 
      CALL DAXPY(y(i),z(:,i),x)
    END DO
    IF (rho<=Tol*rho0) THEN
      EXIT 
    END IF
  END DO
  CALL Deallocate(v)
  CALL Deallocate(z)
    
END SUBROUTINE FGMRES

SUBROUTINE FGMRESV(x,b,MaxIter,Tol)
  TYPE(PressureVelocity), DIMENSION(:) :: b,x
  INTEGER :: MaxIter
  REAL(RealKind) :: Tol
  TYPE(PressureVelocity), POINTER, DIMENSION(:) :: w
  TYPE(PressureVelocity), POINTER, DIMENSION(:,:) ::  &
       v,z
  REAL(RealKind) :: beta,delta,c(StepGMRES),s(StepGMRES),rho,rho0,temp
  REAL(RealKind) :: H(StepGMRES+1,StepGMRES)
  REAL(RealKind) :: bHat(StepGMRES+1)
  REAL(RealKind) :: y(StepGMRES)
  INTEGER :: Iter,i,k,nR
  INTEGER :: ii
  CALL Allocate(v,StepGMRES+1)
  CALL Allocate(z,StepGMRES)
  CALL Allocate(w)
! Compute initial Residual w=b-A*x
  DO Iter=1,MaxIter
    CALL MatVec1(v(:,1),x)
    CALL DXPBY(b,-One,v(:,1))
    w=v(:,1)
    CALL Exchange(w)
    beta=SQRT(w*v(:,1))
    IF (Iter==1) THEN
      rho0=beta
    END IF
    IF (beta<=Zero) THEN
      EXIT
    END IF
    CALL DSCALE(One/beta,v(:,1))
    bHat(1)=beta*One
    S2:DO i=1,StepGMRES
      CALL PSolve1(z(:,i),v(:,i))
      CALL MatVec1(v(:,i+1),z(:,i))
      w=v(:,i+1)
      CALL Exchange(w)
      DO k=1,i
        H(k,i)=v(:,k)*w
      END DO
      DO k=1,i
        CALL DAXPY(-H(k,i),v(:,k),v(:,i+1))
      END DO
      w=v(:,i+1)
      CALL Exchange(w)
      H(i+1,i)=SQRT(w*v(:,i+1))
      CALL DSCALE(One/H(i+1,i),v(:,i+1))
      DO k=2,i
        Temp=c(k-1)*H(k-1,i)+s(k-1)*H(k,i)
        H(k,i)=-s(k-1)*H(k-1,i)+c(k-1)*H(k,i)
        H(k-1,i)=Temp 
      END DO
      delta=SQRT(H(i,i)*H(i,i)+H(i+1,i)*H(i+1,i))
      c(i)=H(i,i)/delta
      s(i)=H(i+1,i)/delta
      H(i,i)=c(i)*H(i,i)+s(i)*H(i+1,i)
      bHat(i+1)=-s(i)*bHat(i)
      bHat(i)=c(i)*bHat(i)
      rho=ABS(bHat(i+1))
      IF (MyId==0) THEN
        WRITE(*,*) 'Rho',rho,rho0,Tol,Tol*rho0
      END IF
      IF (rho<=Tol*rho0) THEN
!     IF (rho<=Tol) THEN
        nR=i
        EXIT S2
      END IF
      nr=StepGMRES
    END DO S2
    y(nR)=bHat(nR)/H(nR,nR)
    DO k=nR-1,1,-1
      y(k)=bHat(k)
      DO i=k+1,nR
        y(k)=y(k)-H(k,i)*y(i)
      END DO
      y(k)=y(k)/H(k,k)
    END DO
    DO i=1,nR 
      CALL DAXPY(y(i),z(:,i),x)
    END DO
    IF (rho<=Tol*rho0) THEN
      EXIT 
    END IF
  END DO
  CALL Deallocate(v)
  CALL Deallocate(w)
  CALL Deallocate(z)
    
END SUBROUTINE FGMRESV

SUBROUTINE QMR(x,b,MaxIter,Tol)
  TYPE(PressureVelocity), DIMENSION(:) :: b,x
  INTEGER :: MaxIter
  REAL(RealKind) :: Tol
  TYPE(PressureVelocity), POINTER, DIMENSION(:) ::  &
       v,Pv,p,d,t
  REAL(RealKind) :: alpha,beta,eta,theta,thetaOld
  REAL(RealKind) :: tau,rho,rhoOld,psi,sigma
  REAL(RealKind) :: NormRes0,NormRes
  REAL(RealKind) :: NormRes0T,NormResT
  INTEGER :: Iter
  CALL Allocate(v)
  CALL Allocate(Pv)
  CALL Allocate(p)
  CALL Allocate(d)
  CALL Allocate(t)

! Compute initial Residual t=b-A*x
  x=Zero  !OSSI
  CALL MatVec1(t,x)
  CALL DXPBY(b,-One,t)
  v=t
! Compute norm of initial Residual
  CALL Exchange(t) 
! CALL Output_PressureVelocity(t)
! NormRes0=NormV(t)
  NormRes0=SQRT(Dot2(t,t))
  IF (NormRes0<=Tol) THEN
    IF (MyID==0) THEN
      WRITE(TermUnit,*) 'QMR Iter',0,NormRes0
    END IF
  ELSE
!   Preconditioning
    CALL PSolve1(Pv,v)  
    p=Pv
    d=Zero
    theta=Zero
    tau=Dot2(Pv,Pv)
    rho=p*v
  
    DO Iter=1,MaxIter
!     Compute t=A*p
      CALL MatVec1(t,p)
      sigma=p*t
      IF (sigma==Zero) THEN
        STOP 'sigma gleich Null'
      END IF
      alpha=rho/sigma
!     Compute v=-alpha*t+v
      CALL DAXPY(-alpha,t,v)
!     Preconditioning
      CALL PSolve1(Pv,v)  
      thetaOld=theta
      theta=Dot2(Pv,Pv)/tau
      psi=one/(One+theta)
      tau=tau*theta*psi
      NormResT=SQRT(tau)
!     Compute d=(psi*thetaOld)*d+(psi*alpha)*p
      CALL DAXPBY(psi*alpha,p,psi*thetaOld,d)
!     Compute x=x+d
      CALL DXPY(d,x)
!     Compute Residual t=b-A*x
      CALL MatVec1(t,x)
      CALL DXPBY(b,-One,t)
!     Compute norm of Residual
      CALL Exchange(t) 
!     CALL Output_PressureVelocity(t)
!     NormRes=NormV(t)
      NormRes=SQRT(Dot2(t,t))
!     Convergence check
!     x-Convergence OR rho=0
      IF (MyID==0) THEN
        WRITE(TermUnit,*) 'QMR Iter',Iter,NormRes,NormRes0,Tol
      END IF
      IF (NormRes<=Tol*NormRes0) THEN
        MaxIter=Iter
        EXIT
      END IF
      rhoOld=rho
      rho=v*Pv
      beta=rho/rhoOld
!     Compute p=Pv+beta*p
      CALL DXPBY(Pv,beta,p)
    END DO
  END IF
    
  CALL Deallocate(v)
  CALL Deallocate(Pv)
  CALL Deallocate(p)
  CALL Deallocate(d)
  CALL Deallocate(t)


END SUBROUTINE QMR

SUBROUTINE MinRes(x,b,MaxIter,Tol)


  TYPE(PressureVelocity), DIMENSION(:) :: b,x
  INTEGER :: MaxIter
  REAL(RealKind) :: Tol


  TYPE(PressureVelocity), POINTER, DIMENSION(:) ::  &
       v,v_hat,v_hat_old,Av,w,w_old,w_oold,y
  REAL(RealKind) :: alpha,beta,beta_old,eta
  REAL(RealKind) :: r1,r1_hat,r2,r3
  REAL(RealKind) :: c,c_old,c_oold
  REAL(RealKind) :: s,s_old,s_oold
  REAL(RealKind) :: NormRes0,NormRes
  REAL(RealKind) :: Fac
  INTEGER :: Iter

  REAL(RealKind) :: Res(MaxIter)

  REAL(RealKind), PARAMETER :: NegOne=-1.0e0

  CALL Allocate(v)
  CALL Allocate(v_hat)
  CALL Allocate(v_hat_old)
  CALL Allocate(Av)
  CALL Allocate(w)
  CALL Allocate(w_old)
  CALL Allocate(w_oold)
  CALL Allocate(y)

! Compute initial Residual v_hat=b-A*x
  CALL MatVec1(v_hat,x)
  CALL DXPBY(b,NegOne,v_hat)
! Compute norm of initial Residual
! Preconditioning
  CALL PSolve1(y,v_hat)  
 
  beta=SQRT(y*v_hat)
  NormRes0=beta
  NormRes=beta
  beta_old=1.0e0

  c=1.0e0
  c_old=1.0e0
  s=0.0e0
  s_old=0.0e0

  eta=beta
  
  DO Iter=1,MaxIter
!   Lanczos
!   v=(1.0e0/beta)*y
    Fac=1.0e0/beta
    CALL DSCALE(Fac,y,v)
    y=v_hat
    CALL MatVec1(Av,v)
    alpha=v*Av
!   v_hat=Av-(alpha/beta)*v_hat-(beta/beta_old)*v_hat_old
    CALL DXPBY(Av,-(alpha/beta),v_hat)
    CALL DAXPY(-(beta/beta_old),v_hat_old,v_hat) 
    v_hat_old=y
    CALL PSolve1(y,v_hat)
    beta_old=beta
    beta=SQRT(MAX(y*v_hat,Zero)) 
    IF (beta < 1.e-12) THEN
      WRITE(*,*) 'MinRes beta<0',y*v_hat
      MaxIter=Iter
      STOP   
      EXIT
    END IF

!   QR-Factorization
    c_oold=c_old
    c_old=c
    s_oold=s_old
    s_old=s
     
    r1_hat=c_old*alpha-c_oold*s_old*beta_old;
    r1    =SQRT(r1_hat**2+beta**2)
    r2    =s_old*alpha+c_oold*c_old*beta_old
    r3    =s_oold*beta_old

!   Givens rotation
    c=r1_hat/r1
    s=beta/r1

!   Update
    w_oold=w_old
    w_old=w
!   w=(v-r3*w_oold-r2*w_old)/r1    
    w=v
    CALL DAXPY(-r3,w_oold,w)
    CALL DAXPY(-r2,w_old,w)
!   w=(1.0e0/r1)*w
    Fac=1.0e0/r1
    CALL DSCALE(Fac,w)

!   x=x+c*eta*w
    CALL DAXPY(c*eta,w,x)
    NormRes=NormRes*s

    eta=-s*eta 

!   Convergence check
    WRITE(*,*) 'MinRes Iter',Iter,NormRes,NormRes0
    IF (NormRes<=Tol*NormRes0) THEN
!   IF (NormRes<=Tol) THEN
      MaxIter=Iter
      EXIT
    END IF
  END DO
  CALL Deallocate(v)
  CALL Deallocate(v_hat)
  CALL Deallocate(v_hat_old)
  CALL Deallocate(Av)
  CALL Deallocate(w)
  CALL Deallocate(w_old)
  CALL Deallocate(w_oold)
  CALL Deallocate(y)


END SUBROUTINE MinRes
 

SUBROUTINE CGR(x,b,MaxIter,Tol)

  TYPE(PressureVelocity), DIMENSION(:) :: b,x
  INTEGER :: MaxIter
  REAL(RealKind) :: Tol


  TYPE(PressureVelocity), POINTER, DIMENSION(:) :: &
          r,z,Az,p,Ap,PAp,t
  REAL(RealKind) :: alpha,beta,rho,rhoOld
  REAL(RealKind) :: NormRes,NormRes0
  INTEGER :: Iter

  REAL(RealKind) :: Res(MaxIter)

  REAL(RealKind), PARAMETER :: NegOne=-1.0e0

  CALL Allocate(r)
  CALL Allocate(t)
  CALL Allocate(z)
  CALL Allocate(Az)
  CALL Allocate(p)
  CALL Allocate(Ap)
  CALL Allocate(PAp)
  

! Compute initial Residual r=b-A*x
  CALL MatVec1(r,x)
! r=b-r
  CALL DXPBY(b,NegOne,r)
  t=r
! Compute norm of initial Residual
  CALL Exchange(t)
  NormRes0=NormV(t)
  IF (NormRes0<=Tol) THEN
    IF (MyID==0) THEN
      WRITE(*,*) 'QMR Iter',0,NormRes0
    END IF
  ELSE
    CALL PSolve1(p,r)
    z=p
    NormRes0=SQRT(r*z)
    CALL MatVec1(Ap,p)
    Az=Ap
    rho=Az*z

    DO Iter=1,MaxIter
      CALL PSolve1(PAp,Ap)
      alpha=rho/(PAp*Ap)
!     x=x+alpha*p
      CALL DAXPY(alpha,p,x)
!     r=r-alpha*Ap
      CALL DAXPY(-alpha,Ap,r)
!     z=z-alpha*PAp
      CALL DAXPY(-alpha,PAp,z)
      NormRes=SQRT(r*z)
      t=r
      CALL Exchange(t)
      NormRes=NormV(t)
!
      CALL MatVec1(Az,z)
      rhoOld=rho
      rho=z*Az
      beta=rho/rhoOld
!     p=z+beta*p
      CALL DXPBY(z,beta,p)    
!     Ap=Az+beta*Ap
      CALL DXPBY(Az,beta,Ap)    

!     Convergence check
      IF (MyID==0) THEN
        WRITE(*,*) 'Iter',Iter,'NormRes',NormRes,'NormRes0',NormRes0
      END IF
      IF (NormRes<=Tol) THEN
        MaxIter=Iter
        EXIT
      END IF
    END DO   
  END IF
  CALL Deallocate(r)
  CALL Deallocate(t)
  CALL Deallocate(z)
  CALL Deallocate(Az)
  CALL Deallocate(p)
  CALL Deallocate(Ap)
  CALL Deallocate(PAp)

END SUBROUTINE CGR

 
END MODULE Iter_Mod

