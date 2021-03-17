MODULE Names_Mod
  USE DataType_Mod
  USE BoundaryCond_Mod

  IMPLICIT NONE

  INTEGER :: ic,icE
  INTEGER :: uPosJac,vPosJac,wPosJac
  REAL(RealKind), POINTER :: c(:,:,:,:)
  REAL(RealKind), POINTER :: cB(:,:)
  REAL(RealKind), POINTER :: Acc(:,:)
  REAL(RealKind), POINTER :: ce1(:,:,:,:)
  REAL(RealKind), POINTER :: ce2(:,:,:,:)
  REAL(RealKind), POINTER :: th(:,:,:,:)
  REAL(RealKind), POINTER :: thProf(:,:,:,:)
  REAL(RealKind), POINTER :: WaterLiq(:,:,:,:) 
  REAL(RealKind), POINTER :: the(:,:,:,:)
  REAL(RealKind), POINTER :: tke(:,:,:,:)
  REAL(RealKind), POINTER :: dis(:,:,:,:)
  REAL(RealKind), POINTER :: ome(:,:,:,:)
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoC(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: RhoI(:,:,:,:)
  REAL(RealKind), POINTER :: RhoS(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: qd(:,:,:,:)
  REAL(RealKind), POINTER :: ql(:,:,:,:)
  REAL(RealKind), POINTER :: Nvap(:,:,:,:)
  REAL(RealKind), POINTER :: Ncloud(:,:,:,:)
  REAL(RealKind), POINTER :: Nrain(:,:,:,:)
  REAL(RealKind), POINTER :: Nice(:,:,:,:)
  REAL(RealKind), POINTER :: Nsnow(:,:,:,:)
  REAL(RealKind), POINTER :: ThetaV(:)
  REAL(RealKind), POINTER :: p(:,:,:,:)
  REAL(RealKind), POINTER :: Sound(:,:,:,:)
  REAL(RealKind), POINTER :: Weight(:,:,:,:)
  REAL(RealKind), POINTER :: T(:,:,:,:)
  REAL(RealKind), POINTER :: KinEn(:,:,:,:)
  REAL(RealKind), POINTER :: E(:,:,:,:)
  REAL(RealKind), POINTER :: RhoEn(:,:,:,:)
  REAL(RealKind), POINTER :: EnergyStart(:,:,:,:)
  REAL(RealKind), POINTER :: PressureStart(:,:,:,:)
  REAL(RealKind), POINTER :: SoS(:,:,:,:)
  REAL(RealKind), POINTER :: HeatRate(:,:,:,:)
  REAL(RealKind), POINTER :: ThForcing(:,:,:,:)
  REAL(RealKind), POINTER :: RhoVForcing(:,:,:,:)
  REAL(RealKind), POINTER :: HeightC(:,:,:,:)
  REAL(RealKind), POINTER :: DivC(:,:,:,:)
  REAL(RealKind), POINTER :: dpdTheta(:,:,:,:)
  REAL(RealKind), POINTER :: Cpml(:,:,:,:)
  REAL(RealKind), POINTER :: DH(:,:,:,:)
  REAL(RealKind), POINTER :: DV(:,:,:,:)  
  REAL(RealKind), POINTER :: uF(:,:,:)
  REAL(RealKind), POINTER :: uFOld(:,:,:)
  REAL(RealKind), POINTER :: vFOld(:,:,:)
  REAL(RealKind), POINTER :: wFOld(:,:,:)
  REAL(RealKind), POINTER :: vF(:,:,:)
  REAL(RealKind), POINTER :: wF(:,:,:)
  REAL(RealKind), POINTER :: cFU(:,:,:,:)
  REAL(RealKind), POINTER :: cFV(:,:,:,:)
  REAL(RealKind), POINTER :: cFW(:,:,:,:)
  REAL(RealKind), POINTER :: pFU(:,:,:)
  REAL(RealKind), POINTER :: pFV(:,:,:)
  REAL(RealKind), POINTER :: pFW(:,:,:)
  REAL(RealKind), POINTER :: DUU(:,:,:)
  REAL(RealKind), POINTER :: DUV(:,:,:)
  REAL(RealKind), POINTER :: DUW(:,:,:)
  REAL(RealKind), POINTER :: uC(:,:,:,:)
  REAL(RealKind), POINTER :: vC(:,:,:,:)
  REAL(RealKind), POINTER :: wC(:,:,:,:)
  REAL(RealKind), POINTER :: uCL(:,:,:,:)
  REAL(RealKind), POINTER :: vCL(:,:,:,:)
  REAL(RealKind), POINTER :: wCL(:,:,:,:)
  REAL(RealKind), POINTER :: uCR(:,:,:,:)
  REAL(RealKind), POINTER :: vCR(:,:,:,:)
  REAL(RealKind), POINTER :: wCR(:,:,:,:)
  REAL(RealKind), POINTER :: uE(:,:,:,:)
  REAL(RealKind), POINTER :: vE(:,:,:,:)
  REAL(RealKind), POINTER :: wE(:,:,:,:)
  REAL(RealKind), POINTER :: uEProf(:)
  REAL(RealKind), POINTER :: vEProf(:)
  REAL(RealKind), POINTER :: f(:,:,:,:),fL(:,:,:,:),fR(:,:,:,:)
  REAL(RealKind), POINTER :: uFPGrad(:,:,:)
  REAL(RealKind), POINTER :: vFPGrad(:,:,:)
  REAL(RealKind), POINTER :: wFPGrad(:,:,:)
  REAL(RealKind), POINTER :: uRhs(:,:,:,:)
  REAL(RealKind), POINTER :: vRhs(:,:,:,:)
  REAL(RealKind), POINTER :: wRhs(:,:,:,:)
  REAL(RealKind), POINTER :: uFRhs(:,:,:)
  REAL(RealKind), POINTER :: vFRhs(:,:,:)
  REAL(RealKind), POINTER :: wFRhs(:,:,:)
  REAL(RealKind), POINTER :: uRhsL(:,:,:,:)
  REAL(RealKind), POINTER :: vRhsL(:,:,:,:)
  REAL(RealKind), POINTER :: wRhsL(:,:,:,:)
  REAL(RealKind), POINTER :: uRhsR(:,:,:,:)
  REAL(RealKind), POINTER :: vRhsR(:,:,:,:)
  REAL(RealKind), POINTER :: wRhsR(:,:,:,:)
  REAL(RealKind), POINTER :: thRhs(:,:,:,:)
  REAL(RealKind), POINTER :: tkeRhs(:,:,:,:)
  REAL(RealKind), POINTER :: disRhs(:,:,:,:)
  REAL(RealKind), POINTER :: omeRhs(:,:,:,:)
 ! REAL(RealKind), POINTER :: RhoVRhs(:,:,:,:)
  REAL(RealKind), POINTER :: RhoVRhs(:,:,:,:)
  REAL(RealKind), POINTER :: RhoCRhs(:,:,:,:)
  REAL(RealKind), POINTER :: RhoLRhs(:,:,:,:)
  REAL(RealKind), POINTER :: RhoIRhs(:,:,:,:)
  REAL(RealKind), POINTER :: RhoRRhs(:,:,:,:)
  REAL(RealKind), POINTER :: RhoRhs(:,:,:,:)
  REAL(RealKind), POINTER :: DampKoeff(:,:,:,:)
  REAL(RealKind), POINTER :: VelSedi(:,:,:,:)
  TYPE(SpDiag), POINTER :: AT
  TYPE(SpDiag), POINTER :: ATMom
  TYPE(SpDiag), POINTER :: ATPot
  TYPE(SpDiag), POINTER :: AFall
  TYPE(SpDiag), POINTER :: AFallRhoL
  TYPE(Vec4_T), POINTER :: AS(:)
  TYPE(BoundaryCon_T), POINTER :: BC 
  REAL(RealKind) :: vFall=0.0d0
  REAL(RealKind) :: dtAdvec
  INTEGER, POINTER :: DiagP(:)
  INTEGER, POINTER :: Permu(:)
  INTEGER :: Diag
  REAL(RealKind), POINTER :: cAve(:),cProf(:),RhoAve(:),wSubs(:),DampProf(:)
  Real(RealKind), TARGET, PUBLIC, ALLOCATABLE :: AccRain(:,:)

CONTAINS

SUBROUTINE DomainSet(ib)

  INTEGER, INTENT(IN) :: ib

  CALL Set(Floor(ib))

END SUBROUTINE DomainSet


SUBROUTINE SetVelocityFace(ibLoc)

  INTEGER, INTENT(IN) :: ibLoc

  uF=>VelocityFaceAct(ibLoc)%uF
  vF=>VelocityFaceAct(ibLoc)%vF
  wF=>VelocityFaceAct(ibLoc)%wF

END SUBROUTINE SetVelocityFace

END MODULE Names_Mod
