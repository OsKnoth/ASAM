MODULE IntLinPeer_Mod

  USE Parallel_Mod
  USE Floor_Mod
  USE Control_Mod
  USE DataType_Mod
  USE Iter_Mod,CG=>FGMRES
  USE JacAccGrav_Mod
  USE BoundaryCondition_Mod
  USE Function_Mod
  USE Rhs_Mod
  USE LinPeerMethods_Mod
  USE IntRos_Mod

  IMPLICIT NONE

  INTEGER, PRIVATE :: Iter
  INTEGER, PRIVATE :: iStage
  INTEGER, PRIVATE :: nStage
  INTEGER, PRIVATE, SAVE :: InitLinP=0
  REAL(RealKind), PRIVATE :: beta0Loc

  REAL(RealKind), PRIVATE, POINTER :: GammaLinPeer(:)
  REAL(RealKind), PRIVATE, POINTER :: cLinPeer(:)
  REAL(RealKind), PRIVATE, POINTER :: A0LinPeer(:,:)
  REAL(RealKind), PRIVATE, POINTER :: ALinPeer(:,:)
  REAL(RealKind), PRIVATE, POINTER :: LUULinPeer(:,:)
  REAL(RealKind), PRIVATE, POINTER :: LUU0LinPeer(:,:)
  INTEGER, POINTER :: PULinPeer(:)
  INTEGER, POINTER :: PU0LinPeer(:)
  REAL(RealKind), POINTER :: GLinPeer(:,:)
  REAL(RealKind), POINTER :: HLinPeer(:,:)
  REAL(RealKind), POINTER :: SLinPeer(:,:)
  REAL(RealKind), POINTER :: BLinPeer(:,:)
  REAL(RealKind), POINTER :: RLinPeer(:,:)
  REAL(RealKind), POINTER :: LUGLinPeer(:,:)
  REAL(RealKind), POINTER :: LUALinPeer(:,:)
  REAL(RealKind), POINTER :: LUBLinPeer(:,:)
  INTEGER, POINTER :: PGLinPeer(:)
  INTEGER, POINTER :: PALinPeer(:)

  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecVelC(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecRhs(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecVelF(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecWC(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecWF(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VecW(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsVec(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: IncrC(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: IncrF(:)
 


CONTAINS


SUBROUTINE AssignMethP2(Method)

  TYPE(LinP2Method_T) :: Method

  nStage=Method%nStage
  beta0Loc=Method%GammaLinPeer(1)
  GammaLinPeer=>Method%GammaLinPeer
  cLinPeer=>Method%cLinPeer
  GLinPeer=>Method%GLinPeer
  HLinPeer=>Method%HLinPeer
  SLinPeer=>Method%SLinPeer
  BLinPeer=>Method%BLinPeer
  ALinPeer=>Method%ALinPeer
  RLinPeer=>Method%RLinPeer
  LUGLinPeer=>Method%LUGLinPeer
  LUALinPeer=>Method%LUALinPeer 
  LUBLinPeer=>Method%LUALinPeer
  PGLinPeer=>Method%PGLinPeer
  PALinPeer=>Method%PALinPeer
END SUBROUTINE AssignMethP2

SUBROUTINE InitLinPeer(VelC,VelF,dtAct,Time)

  TYPE (Vector4Cell_T), POINTER :: VelC(:)
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  REAL(RealKind) :: dtAct,Time

  INTEGER :: iStage
  REAL(RealKind) :: dtLoc
  INTEGER :: iSmall,nSmall=3
  INTEGER :: ComponentsVel
  CHARACTER*20 :: MethodRos='Ros3AMF'

  IF (uPosL>0) THEN
    ComponentsVel=6
  ELSE
    ComponentsVel=0
  END IF
  IF (InitLinP==0) THEN
    CALL MethodsLinP
    IF (Method=='GL5LPeer') THEN
!     Gerisch Lang    
      nStage=GeriLang5%nStage
      beta0Loc=GeriLang5%GammaLinPeer(1)
      GammaLinPeer=>GeriLang5%GammaLinPeer
      cLinPeer=>GeriLang5%cLinPeer
      ALinPeer=>GeriLang5%ALinPeer
      A0LinPeer=>GeriLang5%A0LinPeer
      LUULinPeer=>GeriLang5%LUULinPeer
      LUU0LinPeer=>GeriLang5%LUU0LinPeer
      PULinPeer=>GeriLang5%PULinPeer
      PU0LinPeer=>GeriLang5%PU0LinPeer
    ELSE IF (Method=='GL4LPeer') THEN
      WRITE(*,*) Method
      nStage=GeriLang4%nStage
      beta0Loc=GeriLang4%GammaLinPeer(1)
      GammaLinPeer=>GeriLang4%GammaLinPeer
      cLinPeer=>GeriLang4%cLinPeer
      ALinPeer=>GeriLang4%ALinPeer
      A0LinPeer=>GeriLang4%A0LinPeer
      LUULinPeer=>GeriLang4%LUULinPeer
      LUU0LinPeer=>GeriLang4%LUU0LinPeer
      PULinPeer=>GeriLang4%PULinPeer
      PU0LinPeer=>GeriLang4%PU0LinPeer
    ELSE IF (Method=='JebLPeer') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(Jebens)
    ELSE IF (Method=='JebLPeer0') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebensNeu1)
    ELSE IF (Method=='JebLPeer1') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebensNeu2)
    ELSE IF (Method=='JebLPeer2') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(ImplPeer2)
    ELSE IF (Method=='JebLPeer3') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(ImplPeer3)
    ELSE IF (Method=='JeRLPeer3') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebRosRK3)
    ELSE IF (Method=='JeRLPeer2') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebRos)
      nStage=JebRos%nStage
    ELSE IF (Method=='JeELPeer') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebIEul)
    ELSE IF (Method=='JeWLPeer') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebWei)
    ELSE IF (Method=='JeCLPeer') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebCon)
    END IF
    IF (Method(1:2)=='GL') THEN
      ALLOCATE(VecVelC(nStage))
      ALLOCATE(VecVelF(nStage))
      ALLOCATE(VecWC(nStage))
      ALLOCATE(VecWF(nStage))
  
      DO iStage=1,nStage
        CALL Allocate(VecVelC(iStage)%Vec,1,VectorComponentsM-ComponentsVel)
        CALL Allocate(VecWC(iStage)%Vec,1,VectorComponentsM-ComponentsVel)
        CALL Copy(VelC,VecVelC(iStage)%Vec)
        CALL Copy(VelC,VecWC(iStage)%Vec)
        CALL Allocate(VecVelF(iStage)%VecF)
        CALL Allocate(VecWF(iStage)%VecF)
        CALL Copy(VelF,VecVelF(iStage)%VecF)
        CALL Copy(VelF,VecWF(iStage)%VecF)
      END DO
      CALL Allocate(VecW,1,VectorComponentsM)
      CALL Allocate(RhsVec,VelC)
      CALL Allocate(IncrC,VelC)
      CALL Allocate(IncrF)
    ELSE IF(Method(1:2)=='Je') THEN
      ALLOCATE(VecVelC(nStage))
      IF (uPosL>0) THEN
        ALLOCATE(VecVelF(nStage))
      END IF
      ALLOCATE(VecRhs(nStage))
      CALL Allocate(IncrC,VelC)
    ELSE IF (Method=='JeRLPeer') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebRos)
      nStage=JebRos%nStage
    ELSE IF (Method=='JeELPeer') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebIEul)
    ELSE IF (Method=='JeWLPeer') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebWei)
    ELSE IF (Method=='JeCLPeer') THEN
      WRITE(*,*) Method
      CALL AssignMethP2(JebCon)
    END IF
    IF (Method(1:2)=='GL') THEN
      ALLOCATE(VecVelC(nStage))
      ALLOCATE(VecVelF(nStage))
      ALLOCATE(VecWC(nStage))
      ALLOCATE(VecWF(nStage))
  
      DO iStage=1,nStage
        CALL Allocate(VecVelC(iStage)%Vec,1,VectorComponentsM-ComponentsVel)
        CALL Allocate(VecWC(iStage)%Vec,1,VectorComponentsM-ComponentsVel)
        CALL Copy(VelC,VecVelC(iStage)%Vec)
        CALL Copy(VelC,VecWC(iStage)%Vec)
        CALL Allocate(VecVelF(iStage)%VecF)
        CALL Allocate(VecWF(iStage)%VecF)
        CALL Copy(VelF,VecVelF(iStage)%VecF)
        CALL Copy(VelF,VecWF(iStage)%VecF)
      END DO
      CALL Allocate(VecW,1,VectorComponentsM)
      CALL Allocate(RhsVec,VelC)
      CALL Allocate(IncrC,VelC)
      CALL Allocate(IncrF)
    ELSE IF(Method(1:2)=='Je') THEN
      ALLOCATE(VecVelC(nStage))
      IF (uPosL>0) THEN
        ALLOCATE(VecVelF(nStage))
      END IF
      ALLOCATE(VecRhs(nStage))
      CALL Allocate(IncrC,VelC)
      DO iStage=1,nStage
        CALL Allocate(VecVelC(iStage)%Vec,0,VectorComponentsM-ComponentsVel)
        IF (uPosL>0) THEN
          CALL Allocate(VecVelF(iStage)%VecF)
          CALL Copy(VelF,VecVelF(iStage)%VecF)
        END IF
        CALL Copy(VelC,IncrC)
        IF (iStage>1) THEN
          dtLoc=(cLinPeer(iStage)-cLinPeer(1))*dtAct/nSmall
          IF (uPosL>0) THEN
            CALL Copy(VecVelF(iStage)%VecF,VelF)
          END IF
          IF (InitPeer) THEN
            DO iSmall=1,nSmall
              CALL IntRos(VelF,IncrC,dtLoc,Time,MethodRos)
            END DO
          END IF    
        END IF
        IF (uPosL>0) THEN
          CALL Copy(VelF,VecVelF(iStage)%VecF)
        END IF
        CALL Copy(IncrC,VecVelC(iStage)%Vec)
      END DO
      IF (InitPeer) THEN
        CALL FinalRos
      END IF  
      DO iStage=1,nStage
        CALL Copy(VecVelC(iStage)%Vec,VelC)
        CALL Allocate(VecRhs(iStage)%Vec,1,VectorComponentsM)
        IF (uPosL>0) THEN
          CALL BoundaryVelocity(VecVelF(iStage)%VecF,Time)
          CALL VelocityFaceToCellLR(VecVelF(iStage)%VecF,VelC)
          CALL BoundaryVelocity(VelC,Time)
          CALL Copy(VecVelF(iStage)%VecF,VelF)
        END IF
        CALL BoundaryCondition(VelC,Time)
        CALL ExchangeCell(VelC)
        CALL PrepareF(VelC,VelF,Time)
        VecRhs(iStage)%Vec=Zero
        dtAdvec=dtAct
        CALL Fcn(VelC,VelF,VecRhs(iStage)%Vec,Time,dtAct)
      END DO
      IF (uPosL>0) THEN
        CALL Copy(VecVelF(nStage)%VecF,VelF)
      END IF
      beta0=beta0Loc
      IF (JacTransport) THEN
        JacMet=Zero
        dt=dtAct
        CALL Jac(VelC,VelF,JacMet,Time)
      END IF  
      IF (uPosL>0.AND.PGradient.AND.JacSound) THEN
        dtP=dtAct
        CALL JacAccGrav(VelC)
      END IF

      CALL Allocate(RhsVec,VelC)
      IF (uPosL>0) THEN
        CALL Allocate(IncrF)
      END IF
      CALL Copy(VecVelC(nStage)%Vec,VelC)
    END IF
    InitLinP=1
    beta0=beta0Loc
  END IF

END SUBROUTINE InitLinPeer

SUBROUTINE StageLinPeerGL(VelC,VelF,dtAct,Time)
  TYPE (Vector4Cell_T), POINTER :: VelC(:)
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  REAL(RealKind) :: dtAct
  REAL(RealKind) :: Time

  REAL(RealKind) :: TimeF,dtLoc,Tol
  INTEGER :: iStage,jStage
  INTEGER :: ipStage,jpStage
  INTEGER :: Iter

  DO ipStage=1,nStage
    iStage=PULinPeer(ipStage)
    CALL ScaleV(LUULinPeer(iStage,iStage),VecWF(iStage)%VecF)
    CALL ScaleV(LUULinPeer(iStage,iStage),VecWC(iStage)%Vec)
    DO jpStage=ipStage+1,nStage
      jStage=PULinPeer(jpStage)
      CALL Axpy(LUULinPeer(iStage,jStage),VecWF(jStage)%VecF,VecWF(iStage)%VecF)
      CALL Axpy(LUULinPeer(iStage,jStage),VecWC(jStage)%Vec,VecWC(iStage)%Vec)
    END DO
  END DO
  DO ipStage=nStage,1,-1
    iStage=PULinPeer(ipStage)
    DO jpStage=1,ipStage-1
      jStage=PULinPeer(jpStage)
      CALL Axpy(LUULinPeer(iStage,jStage),VecWF(jStage)%VecF,VecWF(iStage)%VecF)
      CALL Axpy(LUULinPeer(iStage,jStage),VecWC(jStage)%Vec,VecWC(iStage)%Vec)
    END DO
  END DO

  DO ipStage=1,nStage
    iStage=PU0LinPeer(ipStage)
    CALL ScaleV(LUU0LinPeer(iStage,iStage),VecVelF(iStage)%VecF)
    CALL ScaleV(LUU0LinPeer(iStage,iStage),VecVelC(iStage)%Vec)
    DO jpStage=ipStage+1,nStage
      jStage=PU0LinPeer(jpStage)
      CALL Axpy(LUU0LinPeer(iStage,jStage),VecVelF(jStage)%VecF,VecVelF(iStage)%VecF)
      CALL Axpy(LUU0LinPeer(iStage,jStage),VecVelC(jStage)%Vec,VecVelC(iStage)%Vec)
    END DO
  END DO
  DO ipStage=nStage,1,-1
    iStage=PU0LinPeer(ipStage)
    DO jpStage=1,ipStage-1
      jStage=PU0LinPeer(jpStage)
      CALL Axpy(LUU0LinPeer(iStage,jStage),VecVelF(jStage)%VecF,VecVelF(iStage)%VecF)
      CALL Axpy(LUU0LinPeer(iStage,jStage),VecVelC(jStage)%Vec,VecVelC(iStage)%Vec)
    END DO
  END DO

  DO iStage=1,nStage
    DO jStage=1,iStage-1
      CALL Ax1mx2py(ALinPeer(iStage,jStage),VecVelC(jStage)%Vec,VecWC(jStage)%Vec,VecWC(iStage)%Vec)
      CALL Ax1mx2py(ALinPeer(iStage,jStage),VecVelF(jStage)%VecF,VecWF(jStage)%VecF,VecWF(iStage)%VecF)
      CALL Ax1mx2py(A0LinPeer(iStage,jStage),VecVelC(jStage)%Vec,VecWC(jStage)%Vec,VecVelC(iStage)%Vec)
      CALL Ax1mx2py(A0LinPeer(iStage,jStage),VecVelF(jStage)%VecF,VecWF(jStage)%VecF,VecVelF(iStage)%VecF)
    END DO
    CALL Copy(VecVelC(iStage)%Vec,VelC)
    IF (uPosl>0) THEN
      CALL VelocityFaceToCellLR(VecVelF(iStage)%VecF,VelC)
    END IF
    CALL Copy(VecWC(iStage)%Vec,VecW)
    IF (uPosl>0) THEN
      CALL VelocityFaceToCellLR(VecWF(iStage)%VecF,VecW)
    END IF
    dtLoc=(cLinPeer(iStage)-cLinPeer(iStage-1))*dtAct
    dtLoc=dtAct
    TimeF=Time+cLinPeer(iStage)*dtAct
    IF (uPosL>0) THEN
      CALL BoundaryVelocity(VecVelF(iStage)%VecF,TimeF)
      CALL VelocityFaceToCellLR(VecVelF(iStage)%VecF,VelC)
      CALL BoundaryVelocity(VelC,TimeF)
    END IF
    CALL BoundaryCondition(VelC,TimeF)
    CALL ExchangeCell(VelC)
    CALL ExchangeCell(VecW)
    CALL PrepareF(VelC,VecVelF(iStage)%VecF,TimeF)
    IF (iStage==1) THEN
      dt=dtAct
      IF (JacTransport) THEN
        JacMet=Zero
        CALL Jac(VelC,VecVelF(iStage)%VecF,JacMet,Time)
      END IF  
      IF (uPosL>0.AND.PGradient.AND.JacSound) THEN
        dtP=dtAct
        CALL JacAccGrav(VelC)
      END IF
    END IF

    RhsVec=Zero
    dtAdvec=dtAct
    CALL Fcn(VelC,VecVelF(iStage)%VecF,RhsVec,TimeF,dt)
    CALL ScaleV(beta0*dtAct,RhsVec)
    CALL Ax1mx2py(One,VecW,VelC,RhsVec)
    Iter=BiCGStabMaxIter
    Tol=BiCGStabTol
    IncrC=Zero
    IF (JacTransport) THEN
      IF (Source) THEN
        CALL SolveSourceSp(RhsVec,JacMet)
      END IF
      IF (Transport) THEN
        CALL BICGStabSp(JacMet,IncrC,RhsVec,Iter,Tol)
      ELSE
        Iter=0
        CALL Copy(RhsVec,IncrC)
      END IF
    ELSE
      CALL Copy(RhsVec,IncrC)
    END IF
    IF (MyId==0) THEN
      WRITE(TermUnit,*) Iter,Tol
    END IF
    CALL ExchangeCell(IncrC)
    IF (uPosL>0) THEN
      CALL VelocityCellToFaceLR(IncrC,IncrF,VelF,TimeF)
      CALL ProjectVelFace(dtAct,IncrF,IncrC,VelC) 
    END IF
    CALL ExchangeCell(IncrC)
    CALL Axpy(One,IncrF,VecVelF(iStage)%VecF)
    CALL Axpy(One,IncrC,VecVelC(iStage)%Vec)
  END DO
  CALL Copy(VecVelF(nStage)%VecF,VelF)
  CALL Copy(VecVelC(nStage)%Vec,VelC)
  DO iStage=1,nStage
    CALL Copy(VecVelC(iStage)%Vec,VecWC(iStage)%Vec)
    CALL Copy(VecVelF(iStage)%VecF,VecWF(iStage)%VecF)
  END DO

END SUBROUTINE StageLinPeerGL

SUBROUTINE StageLinPeerJeb(VelC,VelF,dtAct,Time)
  TYPE (Vector4Cell_T), POINTER :: VelC(:)
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  REAL(RealKind) :: dtAct
  REAL(RealKind) :: Time

  REAL(RealKind) :: TimeF,dtLoc,Tol
  INTEGER :: iStage,jStage
  INTEGER :: ipStage,jpStage
  INTEGER :: Iter

! dt*A*F_{m-1}
  CALL MultVec(dtAct,LUALinPeer,PALinPeer,VecRhs)
! B*Y_{m-1}+dt*A*F_{m-1}
  DO iStage=1,nStage
    VelC=Zero
    DO jStage=1,nStage
      CALL Axpy(BLinPeer(iStage,jStage),VecVelC(jStage)%Vec,VelC)
    END DO
    IF (uPosl>0) THEN
      VelF=Zero
      DO jStage=1,nStage
        CALL Axpy(BLinPeer(iStage,jStage),VecVelF(jStage)%VecF,VelF)
      END DO
      CALL VelocityFaceToCellLR(VelF,VelC)
    END IF
    CALL ExchangeCell(VelC)
    CALL Axpy(One,VelC,VecRhs(iStage)%Vec)
  END DO
! -1/gamma*G*Y_{m-1}
! B*Y_{m-1}+dt*A*F_{m-1}
  CALL MultVec(-One,LUGLinPeer,PGLinPeer,VecVelC)
  IF (uPosl>0) THEN
    CALL MultVecF(-One,LUGLinPeer,PGLinPeer,VecVelF)
  END IF

  DO iStage=1,nStage
    TimeF=Time+cLinPeer(iStage)*dtAct
!   -1/gamma*H*Y_m-1/gamma*G*Y_{m-1}
!   dt*R*F_m+B*Y_{m-1}+dt*A*F_{m-1}
    DO jStage=1,iStage-1
      CALL Axpy(-HLinPeer(iStage,jStage),VecVelC(jStage)%Vec,VecVelC(iStage)%Vec)
      CALL Axpy(dtAct*RLinPeer(iStage,jStage),VecRhs(jStage)%Vec,VecRhs(iStage)%Vec)
    END DO
    CALL Copy(VecVelC(iStage)%Vec,VelC)
    IF (uPosl>0) THEN
      DO jStage=1,iStage-1
        CALL Axpy(-HLinPeer(iStage,jStage),VecVelF(jStage)%VecF,VecVelF(iStage)%VecF)
      END DO
      CALL VelocityFaceToCellLR(VecVelF(iStage)%VecF,VelC)
    END IF
    CALL ExchangeCell(VelC)
    CALL Axpy(-One,VelC,VecRhs(iStage)%Vec)

    IncrC=Zero
    DO jStage=1,iStage-1
      CALL Axpy(SLinPeer(iStage,jStage),VecVelC(jStage)%Vec,IncrC) !OSSI
    END DO
    IF (uPosL>0) THEN
      VelF=Zero
      DO jStage=1,iStage-1
        CALL Axpy(SLinPeer(iStage,jStage),VecVelF(jStage)%VecF,VelF) !OSSI
      END DO
      CALL VelocityFaceToCellLR(VelF,IncrC)
    END IF  
    CALL ExchangeCell(IncrC)
    CALL Axpy(One,IncrC,VecRhs(iStage)%Vec)

    CALL Copy(VecRhs(iStage)%Vec,RhsVec)
    
    Iter=BiCGStabMaxIter
    Tol=BiCGStabTol
    IncrC=Zero
    IF (JacTransport) THEN
      IF (Source) THEN
        CALL SolveSourceSp(RhsVec,JacMet)
      END IF
      IF (Transport) THEN
        CALL BICGStabSp(JacMet,IncrC,RhsVec,Iter,Tol)
      ELSE
        Iter=0
        CALL Copy(RhsVec,IncrC)
      END IF
    ELSE  
      CALL Copy(RhsVec,IncrC)
    END IF  
    IF (MyId==0) THEN
      WRITE(TermUnit,*) Iter,Tol
    END IF
    CALL ExchangeCell(IncrC)
    IF (uPosL>0) THEN
      CALL VelocityCellToFaceLR(IncrC,IncrF,VelF,TimeF)
      CALL ProjectVelFace(dtAct,IncrF,IncrC,VelC) 
      CALL Axpy(One,IncrF,VecVelF(iStage)%VecF)
    END IF
    CALL ExchangeCell(IncrC)
    CALL Axpy(One,IncrC,VecVelC(iStage)%Vec)
    IF (uPosL>0) THEN
      CALL BoundaryVelocity(VecVelF(iStage)%VecF,TimeF)
      CALL VelocityFaceToCellLR(VecVelF(iStage)%VecF,VelC)
      CALL BoundaryVelocity(VelC,TimeF)
      CALL Copy(VecVelF(iStage)%VecF,VelF)
    END IF
    CALL BoundaryCondition(VecVelC(iStage)%Vec,TimeF)
    CALL Copy(VecVelC(iStage)%Vec,VelC)
    CALL ExchangeCell(VelC)
    CALL PrepareF(VelC,VelF,TimeF)
    IF (iStage==nStage) THEN
      IF (iStage<nStage) THEN
        beta0=GammaLinPeer(iStage+1)
      ELSE  
        beta0=GammaLinPeer(1)
      END IF  
      IF (JacTransport) THEN
        JacMet=Zero
        dt=dtAct
        CALL Jac(VelC,VelF,JacMet,Time)
      END IF  
      IF (uPosL>0.AND.PGradient.AND.JacSound) THEN
        dtP=dtAct
        CALL JacAccGrav(VelC)
      END IF  
    END IF
    VecRhs(iStage)%Vec=Zero
    CALL Fcn(VelC,VelF,VecRhs(iStage)%Vec,TimeF,dtAct)
  END DO
  IF (uPosL>0) THEN
    CALL Copy(VecVelF(nStage)%VecF,VelF)
  END IF 
  CALL Copy(VecVelC(nStage)%Vec,VelC)

END SUBROUTINE StageLinPeerJeb

SUBROUTINE IntLinPeer(VelF,VelC,dtAct,Time,ATol,RTol)

  TYPE (VelocityFace_T), POINTER :: VelF(:)
  TYPE (Vector4Cell_T), POINTER :: VelC(:)
  REAL(RealKind) :: dtAct,Time,ATol(:),RTol(:)

  CALL InitLinPeer(VelC,VelF,dtAct,Time)

  IF (Method(1:2)=='GL') THEN
    CALL StageLinPeerGL(VelC,VelF,dtAct,Time)
  ELSE IF (Method(1:2)=='Je') THEN
    CALL StageLinPeerJeb(VelC,VelF,dtAct,Time)
  END IF
  Time=Time+dtAct

END SUBROUTINE IntLinPeer

END MODULE IntLinPeer_Mod
