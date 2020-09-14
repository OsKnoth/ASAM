MODULE IntPeer_Mod

  USE Parallel_Mod
  USE Floor_Mod
  USE Control_Mod
  USE DataType_Mod
  USE Rhs_Mod
  USE Output_Mod
  USE PeerMethods_Mod
  USE Function_Mod
  USE TimeStep_Mod

  IMPLICIT NONE
  INTEGER, PRIVATE :: iStage
  INTEGER, PRIVATE :: InitP=0
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsVecS(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsVecF(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: RhsVecFF(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecRhsVelFS(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecRhsVecCS(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: VelFFast(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsCS(:)
  TYPE (VecVectorSFace_T), POINTER, PRIVATE, SAVE :: VecThetaF(:)
  TYPE (VectorSFace_T), POINTER, PRIVATE, SAVE :: ThetaF(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecPreFacF(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: PreFacF(:)
  TYPE (VecScalarCell_T), POINTER, PRIVATE, SAVE :: VecSoundFac(:)
  TYPE (ScalarCell_T), POINTER, PRIVATE, SAVE :: SoundFac(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VecU(:)

  INTEGER, PRIVATE :: nStage
  REAL(RealKind), PRIVATE, POINTER :: cPeer(:)
  REAL(RealKind), PRIVATE, POINTER :: alphaPeer(:)
  REAL(RealKind), PRIVATE, POINTER :: APeer(:,:)
  REAL(RealKind), PRIVATE, POINTER :: BPeer(:,:)
  REAL(RealKind), PRIVATE, POINTER :: RPeer(:,:)
  REAL(RealKind), PRIVATE, POINTER :: SPeer(:,:)
  REAL(RealKind), PRIVATE, POINTER :: ALUPeer(:,:)
  REAL(RealKind), PRIVATE, POINTER :: BLUPeer(:,:)


CONTAINS

SUBROUTINE InitExpIntPeer(VecT)

  TYPE (Vector4Cell_T), POINTER :: VecT(:)
  INTEGER :: i

  IF (InitP==0) THEN
    Call MethodsP
    IF (TRIM(Method(3:))=='PeerJ') THEN
      cPeer=>PJebens%cPeer
      alphaPeer=>PJebens%alphaPeer
      APeer=>PJebens%APeer
      BPeer=>PJebens%BPeer
      RPeer=>PJebens%RPeer
      SPeer=>PJebens%SPeer
      ALUPeer=>PJebens%ALUPeer
      BLUPeer=>PJebens%BLUPeer
      nStage=PJebens%nStage
      CFLNumber=PJebens%CFLNumber
    ELSE
      WRITE(*,*) 'Falsche Methode'
      WRITE(*,*) Method
      STOP
    END IF
    CALL Allocate(RhsVecS,VecT,0,VectorComponentsM)
    CALL Allocate(RhsVecF,VecT,1,VectorComponentsM-6)
    CALL Allocate(RhsVecFF)
    ALLOCATE(VecRhsVelFS(nStage))
    DO i=1,nStage
      CALL Allocate(VecRhsVelFS(i)%VecF)
    END DO
    ALLOCATE(VecRhsVecCS(nStage))
    DO i=1,nStage
      CALL Allocate(VecRhsVecCS(i)%Vec,VecT,1,VectorComponentsM-6)
    END DO
    ALLOCATE(VecThetaF(nStage))
    DO i=1,nStage
      CALL Allocate(VecThetaF(i)%VecF,VecT,1,VectorComponentsM-6)
    END DO
    CALL Allocate(ThetaF,VecT,1,VectorComponentsM-6)

    ALLOCATE(VecPreFacF(nStage))
    DO i=1,nStage
      CALL Allocate(VecPreFacF(i)%VecF)
    END DO
    CALL Allocate(PreFacF)

    ALLOCATE(VecSoundFac(nStage))
    DO i=1,nStage
      CALL Allocate(VecSoundFac(i)%Vec)
    END DO
    CALL Allocate(SoundFac)

    CALL Allocate(VelFFast)
    CALL Allocate(RhsCS,VecT,1,VectorComponentsM-6)
    RhsCS=Zero
    CALL Allocate(VecU,uPosl,wPosR)
  END IF

END SUBROUTINE InitExpIntPeer

SUBROUTINE ExpIntPeer(VecVelF,VecVelC,dtAct,Time,ATol,RTol)
  
  TYPE (VecVelocityFace_T), POINTER :: VecVelF(:)
  TYPE (VecVector4Cell_T), POINTER :: VecVelC(:)
  REAL(RealKind) :: dtAct,Time,ATol(:),RTol(:)
  
  INTEGER :: i,is,nsLoc,jStage
  REAL(RealKind) :: dTau,dtLoc

  TimeAct=Time
! CALL InitExpIntPeer(VecVelC(1)%Vec)
  IF (InitP==0) THEN
    DO iStage=1,nStage
      CALL BoundaryVelocity(VecVelF(iStage)%VecF,Time)
      CALL VelocityFaceToCellLR(VecVelF(iStage)%VecF,VecU)
      CALL BoundaryVelocity(VecU,Time)
      CALL ExchangeCell(VecU)
      CALL BoundaryCondition(VecVelC(iStage)%Vec,VecVelF(iStage)%VecF,Time)
!----------------------------------------
!   Part of PrepareF
    CALL PrepareFEx(VecVelC(iStage)%Vec,VecVelF(iStage)%VecF,VecU,Time)

!   Part of PrepareF
!----------------------------------------
      RhsVecS=Zero 
      CALL FcnMetSlow(VecU,VecVelC(iStage)%Vec,VecVelF(iStage)%VecF,RhsVecS &
                     ,VecThetaF(iStage)%VecF,VecPreFacF(iStage)%VecF,VecSoundFac(iStage)%Vec,Time,dt)
      CALL ExchangeCell(RhsVecS)
      CALL VelocityCellToFaceLR(RhsVecS,VecRhsVelFS(iStage)%VecF,VecVelF(iStage)%VecF,Time)
      CALL Copy(RhsVecS,VecRhsVecCS(iStage)%Vec)
    END DO
    InitP=1
  END IF
! A*F, a is decomposed by A=L*R, where L is a permuted lower triangular matrix
! and U an upper triangular matrix with 1 in the diagonal
  DO iStage=1,nStage 
    CALL ScaleV(ALUPeer(iStage,iStage),VecRhsVelFS(iStage)%VecF)
    CALL ScaleV(ALUPeer(iStage,iStage),VecRhsVecCS(iStage)%Vec)
    CALL ScaleV(ALUPeer(iStage,iStage),VecThetaF(iStage)%VecF)
    CALL ScaleV(ALUPeer(iStage,iStage),VecPreFacF(iStage)%VecF)
    CALL ScaleV(ALUPeer(iStage,iStage),VecSoundFac(iStage)%Vec)
    CALL ScaleV(BLUPeer(iStage,iStage),VecVelF(iStage)%VecF)
    CALL ScaleV(BLUPeer(iStage,iStage),VecVelC(iStage)%Vec)
    DO jStage=iStage+1,nStage
      CALL Axpy(ALUPeer(iStage,jStage),VecRhsVelFS(jStage)%VecF,VecRhsVelFS(iStage)%VecF)
      CALL Axpy(ALUPeer(iStage,jStage),VecRhsVecCS(jStage)%Vec,VecRhsVecCS(iStage)%Vec)
      CALL Axpy(ALUPeer(iStage,jStage),VecThetaF(jStage)%VecF,VecThetaF(iStage)%VecF)
      CALL Axpy(ALUPeer(iStage,jStage),VecPreFacF(jStage)%VecF,VecPreFacF(iStage)%VecF)
      CALL Axpy(ALUPeer(iStage,jStage),VecSoundFac(jStage)%Vec,VecSoundFac(iStage)%Vec)
      CALL Axpy(BLUPeer(iStage,jStage),VecVelF(jStage)%VecF,VecVelF(iStage)%VecF)
      CALL Axpy(BLUPeer(iStage,jStage),VecVelC(jStage)%Vec,VecVelC(iStage)%Vec)
    END DO
  END DO
  DO iStage=nStage,1,-1
    DO jStage=1,iStage-1
      CALL Axpy(ALUPeer(iStage,jStage),VecRhsVelFS(jStage)%VecF,VecRhsVelFS(iStage)%VecF)
      CALL Axpy(ALUPeer(iStage,jStage),VecRhsVecCS(jStage)%Vec,VecRhsVecCS(iStage)%Vec)
      CALL Axpy(ALUPeer(iStage,jStage),VecThetaF(jStage)%VecF,VecThetaF(iStage)%VecF)
      CALL Axpy(ALUPeer(iStage,jStage),VecPreFacF(jStage)%VecF,VecPreFacF(iStage)%VecF)
      CALL Axpy(ALUPeer(iStage,jStage),VecSoundFac(jStage)%Vec,VecSoundFac(iStage)%Vec)
      CALL Axpy(BLUPeer(iStage,jStage),VecVelF(jStage)%VecF,VecVelF(iStage)%VecF)
      CALL Axpy(BLUPeer(iStage,jStage),VecVelC(jStage)%Vec,VecVelC(iStage)%Vec)
    END DO
  END DO

  dTau=dtAct/ns
  DO iStage=1,nStage
    DO jStage=1,iStage-1
      CALL Axpy(RPeer(iStage,jStage),VecRhsVelFS(jStage)%VecF,VecRhsVelFS(iStage)%VecF)
      CALL Axpy(RPeer(iStage,jStage),VecRhsVecCS(jStage)%Vec,VecRhsVecCS(iStage)%Vec)
      CALL Axpy(RPeer(iStage,jStage),VecThetaF(jStage)%VecF,VecThetaF(iStage)%VecF)
      CALL Axpy(RPeer(iStage,jStage),VecPreFacF(jStage)%VecF,VecPreFacF(iStage)%VecF)
      CALL Axpy(RPeer(iStage,jStage),VecSoundFac(jStage)%Vec,VecSoundFac(iStage)%Vec)
      CALL Axpy(SPeer(iStage,jStage),VecVelF(jStage)%VecF,VecVelF(iStage)%VecF)
      CALL Axpy(SPeer(iStage,jStage),VecVelC(jStage)%Vec,VecVelC(iStage)%Vec)
    END DO

    nsLoc=CEILING(ns*ABS(alphaPeer(iStage)))
    dtLoc=dtAct*alphaPeer(iStage)
    dTau=dtLoc/nsLoc
    WRITE(*,*) 'nsLoc',nsLoc,'ns',ns
    CALL ForBack1Peer(dTau,nsLoc,VecVelF(iStage)%VecF,VecVelC(iStage)%Vec &
!   CALL SToermerVerlet(dTau,nsLoc,VecVelF(iStage)%VecF,VecVelC(iStage)%Vec &
                     ,VecThetaF(iStage)%VecF &
                     ,VecPreFacF(iStage)%VecF &
                     ,VecSoundFac(iStage)%Vec &
                     ,VecRhsVelFS(iStage)%VecF &
                     ,VecRhsVecCS(iStage)%Vec,Time,dtAct)
  ! Compute slow tendencies
    CALL ExchangeCell(VecVelC(iStage)%Vec) !OSSI
    CALL BoundaryVelocity(VecVelF(iStage)%VecF,Time)
    CALL VelocityFaceToCellLR(VecVelF(iStage)%VecF,VecU)
    CALL BoundaryVelocity(VecU,Time)
    CALL ExchangeCell(VecU)
    CALL BoundaryCondition(VecVelC(iStage)%Vec,VecVelF(iStage)%VecF,Time)
!----------------------------------------
!   Part of PrepareF
    CALL PrepareFEx(VecVelC(iStage)%Vec,VecVelF(iStage)%VecF,VecU,Time)

!   Part of PrepareF
!----------------------------------------
    RhsVecS=Zero 
    CALL FcnMetSlow(VecU,VecVelC(iStage)%Vec,VecVelF(iStage)%VecF,RhsVecS &
                   ,VecThetaF(iStage)%VecF,VecPreFacF(iStage)%VecF,VecSoundFac(iStage)%Vec,Time,dt)
    CALL ExchangeCell(RhsVecS)
    CALL VelocityCellToFaceLR(RhsVecS,VecRhsVelFS(iStage)%VecF,VecVelF(iStage)%VecF,Time)
    CALL Copy(RhsVecS,VecRhsVecCS(iStage)%Vec)
    TimeAct=Time+cPeer(iStage)*dtAct
  END DO
  Time=Time+dtAct

END SUBROUTINE ExpIntPeer

SUBROUTINE ForBack1Peer(dTau,ns,VelF,VelC, &
                    ThetaF,         &
                    PreFacF,         &
                    SoundFac,         &
                    RhsVelFS,         &
                    RhsCS,Time,dT)

  REAL(RealKind) :: dTau,Time,dT
  INTEGER :: ns
  TYPE (VelocityFace_T), POINTER :: VelF(:) ! Stufen für Geschwindikeit
  TYPE (Vector4Cell_T), POINTER :: VelC(:)  ! Stufen für Skalare
  TYPE (VectorSFace_T), POINTER :: ThetaF(:)
  TYPE (VelocityFace_T), POINTER :: PreFacF(:)
  TYPE (ScalarCell_T), POINTER :: SoundFac(:)
  TYPE (VelocityFace_T), POINTER :: RhsVelFS(:)
  TYPE (Vector4Cell_T), POINTER :: RhsCS(:)

  INTEGER :: iPhi,is,jStage
  REAL(RealKind) :: Tau,Fac

  DO is=1,ns
    RhsVecFF=Zero
    CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
    CALL Axpy(dTau,RhsVecFF,VelF)
    CALL Axpy(dTau,RhsVelFS,VelF)
    RhsVecF=Zero
    CALL FcnMetFastScalar(VelC,VelF,ThetaF,SoundFac,RhsVecF,Time)
    CALL Ax1px2py(dTau,RhsVecF,RhsCS,VelC)
  END DO  
END SUBROUTINE ForBack1Peer

SUBROUTINE StoermerVerlet(dTau,ns,VelF,VelC, &
                    ThetaF,         &
                    PreFacF,         &
                    SoundFac,         &
                    RhsVelFS,         &
                    RhsCS,Time,dT)

  REAL(RealKind) :: dTau,Time,dT
  INTEGER :: ns
  TYPE (VelocityFace_T), POINTER :: VelF(:) ! Stufen für Geschwindikeit
  TYPE (Vector4Cell_T), POINTER :: VelC(:)  ! Stufen für Skalare
  TYPE (VectorSFace_T), POINTER :: ThetaF(:)
  TYPE (VelocityFace_T), POINTER :: PreFacF(:)
  TYPE (ScalarCell_T), POINTER :: SoundFac(:)
  TYPE (VelocityFace_T), POINTER :: RhsVelFS(:)
  TYPE (Vector4Cell_T), POINTER :: RhsCS(:)

  INTEGER :: iPhi,is,jStage
  REAL(RealKind) :: Tau,Fac
  REAL(RealKind) :: Summe

  RhsVecFF=Zero
  CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
  CALL Axpy(Half*dTau,RhsVecFF,VelF)
  CALL Axpy(Half*dTau,RhsVelFS,VelF)
  DO is=1,ns-1
    RhsVecF=Zero
    CALL FcnMetFastScalar(VelC,VelF,ThetaF,SoundFac,RhsVecF,Time)
    CALL Ax1px2py(dTau,RhsVecF,RhsCS,VelC)
!   CALL Axpy(dTau,RhsVecF,VelC)
!   CALL Axpy(dTau,RhsCS,VelC)
    RhsVecFF=Zero
    CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
    CALL Axpy(dTau,RhsVecFF,VelF)
    CALL Axpy(dTau,RhsVelFS,VelF)
  END DO
  RhsVecF=Zero
  CALL FcnMetFastScalar(VelC,VelF,ThetaF,SoundFac,RhsVecF,Time)
  CALL Ax1px2py(dTau,RhsVecF,RhsCS,VelC)
! CALL Axpy(dTau,RhsVecF,VelC)
! CALL Axpy(dTau,RhsCS,VelC)
  RhsVecFF=Zero
  CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
  CALL Axpy(Half*dTau,RhsVecFF,VelF)
  CALL Axpy(Half*dTau,RhsVelFS,VelF)
END SUBROUTINE StoermerVerlet

END MODULE IntPeer_Mod

