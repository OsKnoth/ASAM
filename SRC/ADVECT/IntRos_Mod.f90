MODULE IntRos_Mod

  USE Parallel_Mod
  USE Floor_Mod
  USE Control_Mod
  USE DataType_Mod
! USE Iter_Mod,CG=>CGR
  USE Iter_Mod,CG=>FGMRES
  USE JacAccGrav_Mod
  USE Function_Mod
  USE Rhs_Mod
  USE Output_Mod
  USE RungeKuttaMethods_Mod
  USE IntSub_Mod

  IMPLICIT NONE


  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsVec(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecIncrF(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecIncrC(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VelCOld(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: VelFOld(:)

  INTEGER, PRIVATE :: nStage
  INTEGER, PRIVATE, SAVE :: InitR=0

  REAL(RealKind) :: alpha(1:4,1:4)
  REAL(RealKind) :: Gamma1(0:4)
  REAL(RealKind) :: beta(1:4,1:4)

CONTAINS

SUBROUTINE InitRos(VecT,Method)

  TYPE (Vector4Cell_T), POINTER :: VecT(:)
  CHARACTER*20 :: Method
  INTEGER :: i

  IF (InitR==0) THEN

    alpha=Zero
    beta=Zero
    gamma1=Zero
    IF (Method=='Ros3AMF') THEN
!     Ros2AMF    
      nStage=2
      beta0=0.5e0+SQRT(3.0e0)/6.0e0
      alpha(1,1)=2.0d0/3.0d0
      alpha(1,2)=0.0d0
      alpha(2,1)=5.0d0/4.0d0
      alpha(2,2)=3.0d0/4.0d0
      beta=Zero
      beta(2,1)=-4.0d0/3.0d0
      gamma1(0)=Zero
      gamma1(1)=Zero
      gamma1(2)=1.0d0
    ELSE IF (Method=='RosEul') THEN
!     RosEul
      beta0=1.0d0
      nStage=1
      alpha(1,1)=1.0d0
    ELSE IF (Method=='Ros3P') THEN
!     Ros3P
      nStage=3
      beta0=7.886751345948129d-1
      alpha(1,1)=1.267949192431123d0
      alpha(2,1)=1.267949192431123d0
      alpha(3,1)=2.000000000000000d+0
      alpha(3,2)=5.773502691896258d-1
      alpha(3,3)=4.226497308103742d-1
    
      beta(2,1)=-1.607695154586736d0
      beta(3,1)=-3.464101615137755d0
      beta(3,2)=-1.732050807568877d0
      alpha=alpha*beta0
      beta=beta*beta0
    ELSE IF (Method=='Ros3Pw') THEN
      WRITE(*,*) 'Method',Method
      nStage=3
      beta0=7.88675134594812865529e-01
      alpha(1,1)=2.0d0
      alpha(2,1)=6.33974596215561403412e-01
      alpha(3,1)=1.63397459621556140341e+00;
      alpha(3,2)=2.94228634059947813384e-01;
      alpha(3,3)=1.07179676972449078320e+00;

      beta(2,1)=-2.53589838486224561365e+00;
      beta(3,1)=-1.62740473580835520728e+00;
      beta(3,2)=-2.74519052838329002952e-01;
      alpha=alpha*beta0
      beta=beta*beta0
    ELSE IF (Method=='RosRK3') THEN
      WRITE(*,*) 'Method',Method
      nStage=3
      beta0=1.0d0
      alpha(1,1)=1.0d0/(3.0d0*beta0)
      alpha(2,1)=(-1.0d0 + 12.0d0*beta0**2.0d0)/(18.0d0*beta0**2.0d0*(-1.0d0 + 4.0d0*beta0))
      alpha(2,2)=1.0d0/(2.0d0*beta0)
      alpha(3,1)=   (1.0d0 - 3.0d0*beta0*(7.0d0 + 16.0d0*beta0*(-2.0d0 + 3.0d0*beta0)))/(36.0d0*beta0**3.0d0*(-1.0d0 + 4.0d0*beta0))
      alpha(3,2)=   (-1.0d0 + 12.0d0*beta0)/(4.0d0*beta0**2.0d0)
      alpha(3,3)=   1.0d0/beta0
      beta(2,1)= (1.0d0 - 12.0d0*beta0**2.0d0)/(beta0**2.0d0*(-9.0d0 + 36.0d0*beta0))
      beta(3,1)= (-1.0d0 + 3.0d0*beta0*(7.0d0 + 16.0d0*beta0*(-2.0d0 + 3.0d0*beta0)))/(36.*beta0**3.0d0*(-1.0d0 + 4.0d0*beta0))
      beta(3,2)= (0.250d0 - 3.0d0*beta0)/beta0**2.0d0

      alpha=alpha*beta0
      beta=beta*beta0
    ELSE IF (Method=='RosTVD') THEN
      WRITE(*,*) 'Method',Method
      nStage=3
      beta0=0.5d0
      alpha(1,1)= 1.0d0/beta0
      alpha(2,1)= (1.0d0 + beta0*(-9.0d0 + 16.0d0*beta0))/(4.*beta0**2.0d0*(-1.0d0 + 4.0d0*beta0))
      alpha(2,2)= 1.0d0/(4.0d0*beta0)
      alpha(3,1)= (1.0d0 + beta0*(-3.0d0 + 2.0d0*beta0*(-17.0d0 + 44.0d0*beta0)))/(12.*beta0**3.0d0*(-1.0d0 +  4.0d0*beta0))
      alpha(3,2)= (1.0d0 + 6.0d0*beta0)/(12.0d0*beta0**2.0d0)
      alpha(3,3)= 2.0d0/(3.*beta0)
      beta(2,1)= beta0**(-2.0d0) - 4.0d0/beta0 + 4.0d0/(-1.0d0 + 4.0d0*beta0)
      beta(3,1)= (-1.0d0 + beta0*(5.0d0 + 8.0d0*(2.0d0 - 7.0d0*beta0)*beta0))/(8.0d0*beta0**3.0d0*(-1.0d0 +  4.0d0*beta0))
      beta(3,2)= -(1.0d0 + 4.0d0*beta0)/(8.*beta0**2.0d0)
      alpha=alpha*beta0
      beta=beta*beta0
    ELSE IF (Method=='RosRK4') THEN
      WRITE(*,*) 'Method',Method
      nStage=4
      beta0=1.0d0/4.0d0
      alpha(1,1)=2.0d0
      alpha(2,1)=2.0d0
      alpha(2,2)=2.0d0
      alpha(3,1)=6.0d0
      alpha(3,2)=10.0d0
      alpha(3,3)=4.0d0
      alpha(4,1)=14.0d0/3.0d0
      alpha(4,2)=20.0d0/3.0d0
      alpha(4,3)=4.0d0/3.0d0
      alpha(4,4)=2.0d0/3.0d0
    
      beta(2,1)=-4.0d0
      beta(3,1)=-6.0d0
      beta(3,2)=-10.0d0
      beta(4,1)=-4.0d0
      beta(4,2)=-12.0d0
      alpha=alpha*beta0
      beta=beta*beta0
    END IF
    CALL Allocate(RhsVec,VecT)
    ALLOCATE(VecIncrC(nStage))
    DO i=1,nStage
      CALL Allocate(VecIncrC(i)%Vec,VecT)
    END DO
    ALLOCATE(VecIncrF(nStage))
    DO i=1,nStage
      CALL Allocate(VecIncrF(i)%VecF)
      VecIncrF(i)%VecF=Zero
    END DO
    CALL Allocate(VelCOld,VecT)
    CALL Allocate(VelFOld)
    InitR=1
  END IF

END SUBROUTINE InitRos

SUBROUTINE FinalRos

  INTEGER :: i

  IF (InitR==1) THEN 
    CALL Deallocate(RhsVec)
    DO i=1,nStage
      CALL Deallocate(VecIncrC(i)%Vec)
    END DO
    DEALLOCATE(VecIncrC)
    DO i=1,nStage
      CALL Deallocate(VecIncrF(i)%VecF)
    END DO
    DEALLOCATE(VecIncrF)
    CALL Deallocate(VelCOld)
    CALL Deallocate(VelFOld)
    InitR=0
  END IF

END SUBROUTINE FinalRos

SUBROUTINE StageRos(VelC,VelF,dtAct,Time)
  TYPE (Vector4Cell_T), POINTER :: VelC(:)
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  REAL(RealKind) :: dtAct
  REAL(RealKind) :: Time

  REAL(RealKind) :: TimeF,dtLoc,Tol
  INTEGER :: iStage
  INTEGER :: Iter,jStage

  DO iStage=1,nStage
    dtLoc=(Gamma1(iStage)-Gamma1(iStage-1))*dtAct
    TimeF=Time+Gamma1(iStage)*dtAct
    IF (uPosL>0) THEN
      CALL BoundaryVelocity(VelF,TimeF)
      CALL VelocityFaceToCellLR(VelF,VelC)
      CALL BoundaryVelocity(VelC,TimeF)
    END IF
    CALL BoundaryCondition(VelC,TimeF)
    CALL ExchangeCell(VelC)
    CALL PrepareF(VelC,VelF,TimeF)
    IF (iStage==1) THEN
      JacMet=Zero
      dt=dtAct
      CALL Jac(VelC,VelF,JacMet,Time)
      IF (uPosL>0.AND.PGradient) THEN
        dtP=dtAct
        CALL JacAccGrav(VelC)
      END IF
    END IF
    RhsVec=Zero
    dtAdvec=dtAct
    CALL Fcn(VelC,VelF,RhsVec,TimeF,dtAct)
    CALL ScaleV(dtAct,RhsVec)
    DO jStage=1,iStage-1
      CALL Axpy(beta(iStage,jStage),VecIncrC(jStage)%Vec,RhsVec) 
    END DO
    Iter=BiCGStabMaxIter
    Tol=BiCGStabTol
    VecIncrC(iStage)%Vec=Zero
    IF (Source) THEN
      CALL SolveSourceSp(RhsVec,JacMet)
    END IF
    IF (Transport) THEN
      CALL BICGStabSp(JacMet,VecIncrC(iStage)%Vec,RhsVec,Iter,Tol)
    ELSE
      Iter=0
      CALL Copy(RhsVec,VecIncrC(iStage)%Vec)
    END IF
    IF (MyId==0) THEN
      WRITE(TermUnit,*) Iter,Tol
    END IF
    CALL ExchangeCell(VecIncrC(iStage)%Vec)
    IF (uPosL>0) THEN
      CALL VelocityCellToFaceLR(VecIncrC(iStage)%Vec,VecIncrF(iStage)%VecF,VelF,TimeF)
      CALL ProjectVelFace(dtAct,VecIncrF(iStage)%VecF,VecIncrC(iStage)%Vec,VelCOld) 
      CALL VelocityFaceToCellLR(VecIncrF(iStage)%VecF,VecIncrC(iStage)%Vec)
    END IF
    CALL ExchangeCell(VecIncrC(iStage)%Vec)
    CALL Copy(VelCOld,VelC)
    CALL Copy(VelFOld,VelF)
    DO jStage=1,iStage
      CALL Axpy(alpha(iStage,jStage),VecIncrF(jStage)%VecF,VelF)
      CALL Axpy(alpha(iStage,jStage),VecIncrC(jStage)%Vec,VelC)
    END DO
  END DO

END SUBROUTINE StageRos

SUBROUTINE IntRos(VelF,VelC,dtAct,Time,Method)

  TYPE (VelocityFace_T), POINTER :: VelF(:)
  TYPE (Vector4Cell_T), POINTER :: VelC(:)
  REAL(RealKind) :: dtAct,Time
  CHARACTER*20 :: Method

  CALL InitRos(VelC,Method)
  CALL Copy(VelF,VelFOld)
  CALL Copy(VelC,VelCOld)
  CALL StageRos(VelC,VelF,dtAct,Time)
  Time=Time+dtAct
END SUBROUTINE IntRos

END MODULE IntRos_Mod

