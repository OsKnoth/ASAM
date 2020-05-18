MODULE PeerMethods_Mod
  USE Parameter_Mod

  IMPLICIT NONE

  TYPE PMethod_T
    INTEGER :: nStage
    REAL(RealKind), POINTER :: cPeer(:)
    REAL(RealKind), POINTER :: alphaPeer(:)
    REAL(RealKind), POINTER :: APeer(:,:)
    REAL(RealKind), POINTER :: BPeer(:,:)
    REAL(RealKind), POINTER :: RPeer(:,:)
    REAL(RealKind), POINTER :: SPeer(:,:)
    REAL(RealKind), POINTER :: ALUPeer(:,:)
    REAL(RealKind), POINTER :: BLUPeer(:,:)
    REAL(RealKind) :: CFLNumber
  END TYPE PMethod_T

  TYPE(PMethod_T) :: PJebens

CONTAINS

SUBROUTINE MethodsP

  INTEGER :: i,j,Info
  INTEGER :: nStage
  REAL(RealKind), POINTER :: cPeer(:)
  REAL(RealKind), POINTER :: alphaPeer(:)
  REAL(RealKind), POINTER :: APeer(:,:)
  REAL(RealKind), POINTER :: BPeer(:,:)
  REAL(RealKind), POINTER :: RPeer(:,:)
  REAL(RealKind), POINTER :: SPeer(:,:)
  REAL(RealKind), POINTER :: ALUPeer(:,:)
  REAL(RealKind), POINTER :: BLUPeer(:,:)

  nStage=3
  PJebens%nStage=nStage
  PJebens%CFLNumber=1.5d0
  ALLOCATE(PJebens%cPeer(nStage))
  ALLOCATE(PJebens%alphaPeer(nStage))
  ALLOCATE(PJebens%APeer(nStage,nStage))
  ALLOCATE(PJebens%BPeer(nStage,nStage))
  ALLOCATE(PJebens%RPeer(nStage,nStage))
  ALLOCATE(PJebens%SPeer(nStage,nStage))
  ALLOCATE(PJebens%ALUPeer(nStage,nStage))
  ALLOCATE(PJebens%BLUPeer(nStage,nStage))
  cPeer=>PJebens%cPeer
  alphaPeer=>PJebens%alphaPeer
  APeer=>PJebens%APeer
  BPeer=>PJebens%BPeer
  RPeer=>PJebens%RPeer
  SPeer=>PJebens%SPeer
  ALUPeer=>PJebens%ALUPeer
  BLUPeer=>PJebens%BLUPeer

  cPeer(1:3)=(/-8.9953162787855270e-002 &
              , 4.6764288306976509e-001 &
              , 1.0000000000000000e+000/)

  APeer(1,:)=(/ 7.2100732200857512e-002,-1.3228042883312888e-001, 1.2650691731921049e-001/)
  APeer(2,:)=(/ 4.7823871966525877e-002,-4.8313723987222795e-001,-1.1633321060462617e-001/)
  APeer(3,:)=(/ 3.2590697144031344e-002, 7.0244009589084278e-002, 1.2867615058926479e-001/)

  BPeer(1,:)=(/-9.6705998384565617e-002, 4.9155986452023448e-001, 6.0514613386433114e-001/)
  BPeer(2,:)=(/-4.7092982628159370e-002, 2.1699465817029365e-001, 5.7208159637221156e-001/)
  BPeer(3,:)=(/-8.9143731284548000e-002, 1.5738303158840131e-001, 1.9732333925866852e-001/)

  RPeer(1,:)=(/ 0.0000000000000000e+000, 0.0000000000000000e+000, 0.0000000000000000e+000/)
  RPeer(2,:)=(/ 1.1066883875756954e+000, 0.0000000000000000e+000, 0.0000000000000000e+000/)
  RPeer(3,:)=(/-5.0202716737489572e-001, 1.0959786066300778e+000, 0.0000000000000000e+000/)

  SPeer(1,:)=(/ 0.0000000000000000e+000, 0.0000000000000000e+000, 0.0000000000000000e+000/)
  SPeer(2,:)=(/ 2.5801672808565412e-001, 0.0000000000000000e+000, 0.0000000000000000e+000/)
  SPeer(3,:)=(/ 3.2693061133974349e-001, 4.0750674909773471e-001, 0.0000000000000000e+000/)

  alphaPeer(1)=SUM(APeer(1,:))+SUM(RPeer(1,:))
  APeer(1,:)=APeer(1,:)/alphaPeer(1)
  RPeer(1,:)=RPeer(1,:)/alphaPeer(1)
  alphaPeer(2)=SUM(APeer(2,:))+SUM(RPeer(2,:))
  APeer(2,:)=APeer(2,:)/alphaPeer(2)
  RPeer(2,:)=RPeer(2,:)/alphaPeer(2)
  alphaPeer(3)=SUM(APeer(3,:))+SUM(RPeer(3,:))
  APeer(3,:)=APeer(3,:)/alphaPeer(3)
  RPeer(3,:)=RPeer(3,:)/alphaPeer(3)
  ALUPeer=APeer
  CALL DGEFA (ALUPeer,3,3,Info)
  BLUPeer=BPeer
  CALL DGEFA (BLUPeer,3,3,Info)

END SUBROUTINE MethodsP

END MODULE PeerMethods_Mod

