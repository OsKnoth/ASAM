MODULE Int_Mod

  USE Parallel_Mod
  USE Floor_Mod
  USE Control_Mod
  USE DataType_Mod
  USE Iter_Mod,CG=>FGMRES
  USE JacAccGrav_Mod
  USE BoundaryCondition_Mod
  USE Function_Mod
  USE Rhs_Mod
  USE Output_Mod
  USE RungeKuttaMethods_Mod
  USE IntSub_Mod
  USE Tools_Mod
  USE Init_Mod
  USE TimeStep_Mod

  IMPLICIT NONE


  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsVecMet(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsVecChem(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsVecS(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsVecF(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecIncrF(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecIncrMetC(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecIncrChemC(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecIncrCG(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: RhsVecFF(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: RhsVecFS(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecRhsVecFS(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecRhsVecCS(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VecMetOld(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VecChemOld(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VelGOld(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VelC1(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VelC2(:)
  TYPE (VecVector4Cell_T), POINTER, PRIVATE, SAVE :: VecVelC(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: VelFOld(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: VelF1(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: VelF2(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: VelFFast(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: RhsVelFS(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: RhsCS(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecVelF(:)
  TYPE (VecVectorSFace_T), POINTER, PRIVATE, SAVE :: VecMethetaF(:)
  TYPE (VectorSFace_T), POINTER, PRIVATE, SAVE :: ThetaF(:)
  TYPE (VecVelocityFace_T), POINTER, PRIVATE, SAVE :: VecPreFacF(:)
  TYPE (VelocityFace_T), POINTER, PRIVATE, SAVE :: PreFacF(:)
  TYPE (VecScalarCell_T), POINTER, PRIVATE, SAVE :: VecSoundFac(:)
  TYPE (ScalarCell_T), POINTER, PRIVATE, SAVE :: SoundFac(:)
  TYPE (Vector4Cell_T), POINTER, PRIVATE, SAVE :: VecU(:)


  INTEGER, PRIVATE :: Iter
  INTEGER, PRIVATE :: iStage
  INTEGER, PRIVATE :: nStage,nPhi
  INTEGER, PRIVATE :: ii,jj
  INTEGER, PRIVATE, SAVE :: InitR=0

  INTEGER, PRIVATE :: Stage=1
  TYPE(Koeff_T), PRIVATE, POINTER :: aRunge(:,:)
  REAL(RealKind), PRIVATE, POINTER :: cRunge(:)
  REAL(RealKind), PRIVATE, POINTER :: dRunge(:,:)
  REAL(RealKind), PRIVATE, POINTER :: gRunge(:,:)
  REAL(RealKind), PRIVATE, POINTER :: dtRunge(:)

  REAL(RealKind) :: alpha(1:4,1:4)
  REAL(RealKind) :: Gamma1(0:4)
  REAL(RealKind) :: beta(1:4,1:4)
  REAL(RealKind) :: delta=-Half
  LOGICAL :: ChemieLoc=.FALSE.
  INTEGER :: kx=90

CONTAINS

SUBROUTINE InitRos3MetC(VecMet,VecChem)

  TYPE (Vector4Cell_T), POINTER :: VecMet(:)
  TYPE (Vector4Cell_T), POINTER :: VecChem(:)
  INTEGER :: i

  REAL(RealKind) :: theta
  REAL(RealKind) :: tt(1:4,1:4)
  REAL(RealKind) :: gamma21,gamma31,gamma32

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
      gamma1(0)=0.0d0
      gamma1(1)=0.666666666666667d0
      gamma1(2)=1.0d0
      delta=-Half
    ELSE IF (Method=='RosRK32') THEN
!     RosRK32    
      nStage=3
      theta=1.0d0/3.0d0*ATAN(SQRT(2.0d0)/4.0d0)
      beta0=1.0d0-Half*SQRT(2.0d0)*COS(theta)+Half*SQRT(6.0d0)*SIN(theta)
      gamma32=Half-3.0d0*beta0
      gamma31=-(6.0d0*beta0*beta0*beta0-12.0d0*beta0*beta0 &
              +6.0d0*(1.0d0+gamma32)*beta0+2.0d0*gamma32*gamma32-Half) &
             /(1.0d0+2.0d0*gamma32)
      gamma21=-(3.0d0*beta0+gamma31+gamma32)
      
      tt=Zero
      tt(1,1)=One
      tt(2,1)=-gamma21/beta0
      tt(2,2)=One
      tt(3,1)=-gamma31/beta0+gamma32*gamma21/beta0**2
      tt(3,2)=-gamma32/beta0
      tt(3,3)=One
      alpha=Zero
      alpha(1,1)=Half
      alpha(2,1)=Half
      alpha(2,2)=Half
      alpha(3,1)=1.0d0/3.0d0
      alpha(3,2)=1.0d0/3.0d0
      alpha(3,3)=1.0d0/3.0d0
      alpha(1:3,1:3)=MATMUL(alpha(1:3,1:3),tt(1:3,1:3))
      beta=Zero
      beta(2,1)=gamma21/beta0
      beta(3,1)=gamma31/beta0
      beta(3,2)=gamma32/beta0
      beta(1:3,1:3)=MATMUL(beta(1:3,1:3),tt(1:3,1:3))
      gamma1(0)=Zero
      gamma1(1)=Zero
      gamma1(2)=1.0d0
      delta=-Half
    ELSE IF (Method=='RosEul') THEN
!     RosEul
      beta0=1.0d0
      nStage=1
      alpha(1,1)=1.0d0
      Gamma1(0)=0.0d0
      Gamma1(1)=1.0d0
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
      nStage=3
      beta0=1.0d0
!     beta0=1.0d0/2.0d0
!     beta0=(3.0d0+SQRT(3.0d0))/6.0d0
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
      Gamma1(0)=0.0d0
      Gamma1(1)=1.0d0/3.0d0
      Gamma1(2)=1.0d0/2.0d0
      Gamma1(3)=1.0d0
    ELSE IF (Method=='RosTVD') THEN
      nStage=3
      beta0=Half
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
    CALL Allocate(RhsVecMet,VecMet)
    CALL Allocate(RhsVecChem,VecChem)
    ALLOCATE(VecIncrMetC(nStage))
    ALLOCATE(VecIncrChemC(nStage))
    DO i=1,nStage
      CALL Allocate(VecIncrMetC(i)%Vec,VecMet)
      CALL Allocate(VecIncrChemC(i)%Vec,VecMet)
    END DO
    ALLOCATE(VecIncrF(nStage))
    DO i=1,nStage
      CALL Allocate(VecIncrF(i)%VecF)
      VecIncrF(i)%VecF=Zero
    END DO
    CALL Allocate(VecMetOld,VecMet)
    CALL Allocate(VecChemOld,VecChem)
    CALL Allocate(VelFOld)
    IF (Anelastic.OR.PseudoIn) THEN
      ALLOCATE(VecIncrCG(nStage))
      DO i=1,nStage
        CALL Allocate(VecIncrCG(i)%Vec,1,1)
      END DO
    CALL Allocate(VelGOld,1,1)
    END IF
    InitR=1
  END IF

END SUBROUTINE InitRos3MetC

SUBROUTINE StageRos(VecMet,VecChem,VelF,VelG,dt,Time)
  TYPE (Vector4Cell_T), POINTER :: VecMet(:)
  TYPE (Vector4Cell_T), POINTER :: VecChem(:)
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  TYPE (Vector4Cell_T), POINTER :: VelG(:)
  REAL(RealKind) :: dt
  REAL(RealKind) :: Time
  INTEGER, SAVE :: JacStart=1

  REAL(RealKind) :: TimeF,dtLoc,Tol
  INTEGER :: iStage
  INTEGER :: Iter,jStage
  INTEGER :: i,j,k
  INTEGER :: ix,iy
  REAL(RealKind) :: Temp

  DO iStage=1,nStage
    TimeF=Time+Gamma1(iStage-1)*dt
    CALL ExchangeCell(VecMet)
    CALL ExchangeCell(VecChem)
    IF (uPosL>0) THEN
      CALL VelocityFaceToCellLR(VelF,VecMet)
      CALL BoundaryVelocity(VecMet,TimeF)
    END IF
    CALL BoundaryCondition(VecMet,VecChem,VelF,TimeF)
    CALL PrepareF(VecMet,VelF,TimeF)
    IF (iStage==1) THEN
      IF (JacTransport) THEN
        JacTrans=Zero
        CALL Jac(VecMet,VecChem,VelF,JacTrans,Time)
      END IF
      IF ((JacSound.AND.uPosL>0.AND.PGradient).OR.Anelastic.OR.PseudoIn) THEN
        dtP=dt
        CALL JacAccGrav(VecMet)
      END IF
    END IF
    RhsVecMet=Zero
    RhsVecChem=Zero
    dtAdvec=dt
    CALL Fcn(VecMet,VecChem,VelF,RhsVecMet,RhsVecChem,TimeF,dt,VelG)
    CALL ScaleV(dt,RhsVecMet)
    CALL ScaleV(dt,RhsVecChem)
    DO jStage=1,iStage-1
      CALL Axpy(beta(iStage,jStage),VecIncrMetC(jStage)%Vec,RhsVecMet) 
      CALL Axpy(beta(iStage,jStage),VecIncrChemC(jStage)%Vec,RhsVecChem) 
    END DO
    IF (Anelastic.OR.PseudoIn) THEN
      VecIncrCG(iStage)%Vec=Zero
      CALL FcnG(VecMet,VelF,VecIncrCG(iStage)%Vec,TimeF,dt)
      CALL ScaleV(dt,VecIncrCG(iStage)%Vec)
    END IF  
    IF (JacTransport) THEN
      VecIncrMetC(iStage)%Vec=Zero
      VecIncrChemC(iStage)%Vec=Zero
      IF (Source) THEN
        CALL SolveSourceMetSp(RhsVecMet,JacTrans)
        CALL SolveSourceChemSp(RhsVecChem,JacTrans)
      END IF
      IF (Transport) THEN
        Tol=BiCGStabTol
        Iter=BiCGStabMaxIter
        CALL BICGStabMetSp(JacTrans,VecIncrMetC(iStage)%Vec,RhsVecMet,Iter,Tol)
        Tol=BiCGStabTol
        Iter=BiCGStabMaxIter
        CALL BICGStabChemSp(JacTrans,VecIncrChemC(iStage)%Vec,RhsVecChem,Iter,Tol)
      ELSE
        Iter=0
        CALL Copy(RhsVecMet,VecIncrMetC(iStage)%Vec)
        CALL Copy(RhsVecChem,VecIncrChemC(iStage)%Vec)
      END IF
    ELSE
      CALL Copy(RhsVecMet,VecIncrMetC(iStage)%Vec)
      CALL Copy(RhsVecChem,VecIncrChemC(iStage)%Vec)
    END IF
    IF (MyId==0) THEN
      WRITE(TermUnit,*) Iter,Tol
    END IF
    CALL ExchangeCell(VecIncrMetC(iStage)%Vec)
    CALL ExchangeCell(VecIncrChemC(iStage)%Vec)
    IF (uPosL>0) THEN
      CALL VelocityCellToFaceLR(VecIncrMetC(iStage)%Vec,VecIncrF(iStage)%VecF,VelF,TimeF)
      IF ((JacSound.AND.uPosL>0.AND.PGradient).OR.Anelastic.OR.PseudoIn) THEN
        IF (Anelastic.OR.PseudoIn) THEN
          CALL ProjectVelFace(dt,VecIncrF(iStage)%VecF,VecIncrMetC(iStage)%Vec,VecMetOld &
                                                      ,VecIncrChemC(iStage)%Vec,VecChemOld,VecG=VecIncrCG(iStage)%Vec) 
        ELSE  
          CALL ProjectVelFace(dt,VecIncrF(iStage)%VecF,VecIncrMetC(iStage)%Vec,VecMetOld, &
                                                       VecIncrChemC(iStage)%Vec,VecChemOld) 
        END IF  
      END IF  
      CALL VelocityFaceToCellLR(VecIncrF(iStage)%VecF,VecIncrMetC(iStage)%Vec)
      CALL Copy(VelFOld,VelF)
      DO jStage=1,iStage
        CALL Axpy(alpha(iStage,jStage),VecIncrF(jStage)%VecF,VelF)
      END DO
    ELSE
      CALL VelocityInit(VelF,UStart,VStart,WStart,Time+Gamma1(iStage)*dt)   
    END IF
    CALL ExchangeCell(VecIncrMetC(iStage)%Vec)
    CALL Copy(VecMetOld,VecMet)
    CALL Copy(VecChemOld,VecChem)
    DO jStage=1,iStage
      CALL Axpy(alpha(iStage,jStage),VecIncrMetC(jStage)%Vec,VecMet)
      CALL Axpy(alpha(iStage,jStage),VecIncrChemC(jStage)%Vec,VecChem)
    END DO
    IF (Anelastic.OR.PseudoIn) THEN
      CALL Copy(VelGOld,VelG)
      DO jStage=1,iStage
        CALL Axpy(alpha(iStage,jStage),VecIncrCG(jStage)%Vec,VelG)
      END DO
    END IF
  END DO
  CONTAINS

SUBROUTINE Div(ix0,ix1,iy0,iy1,iz0,iz1,uF,vF,wF,FU,FV,FW)

  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  REAL(RealKind) :: uF(ix0:ix1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: vF(ix0+1:ix1,iy0:iy1,iz0+1:iz1)
  REAL(RealKind) :: wF(ix0+1:ix1,iy0+1:iy1,iz0:iz1)
  REAL(RealKind) :: FU(ix0-1:ix1+1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: FV(ix0+1:ix1,iy0-1:iy1+1,iz0+1:iz1)
  REAL(RealKind) :: FW(ix0+1:ix1,iy0+1:iy1,iz0-1:iz1+1)

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: tt

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        tt= &
                (uF(ix,iy,iz)*FU(ix,iy,iz)-uF(ix-1,iy,iz)*FU(ix-1,iy,iz) &
                +vF(ix,iy,iz)*FV(ix,iy,iz)-vF(ix,iy-1,iz)*FV(ix,iy-1,iz) &
                +wF(ix,iy,iz)*FW(ix,iy,iz)-wF(ix,iy,iz-1)*FW(ix,iy,iz-1))
      END DO
    END DO
  END DO

END SUBROUTINE Div


END SUBROUTINE StageRos

SUBROUTINE StageRosParcel(VecMet,VecChem,VelF,dt,Time)
  TYPE (Vector4Cell_T), POINTER :: VecMet(:)
  TYPE (Vector4Cell_T), POINTER :: VecChem(:)
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  REAL(RealKind) :: dt
  REAL(RealKind) :: Time

  REAL(RealKind) :: TimeF,dtLoc,Tol
  INTEGER :: iStage
  INTEGER :: Iter,jStage

  DO iStage=1,nStage
    dtLoc=(Gamma1(Stage)-Gamma1(Stage-1))*dt
    TimeF=Time+Gamma1(Stage)*dt
    CALL PrepareF(VecMet,VelF,TimeF)
    IF (iStage==1) THEN
      JacTrans=Zero
      CALL Jac(VecMet,VecChem,VelF,JacTrans,Time)
      dtP=dt
    END IF
    RhsVecMet=Zero
    RhsVecChem=Zero
    dtAdvec=dt
    CALL Fcn(VecMet,VecChem,VelF,RhsVecMet,RhsVecChem,Time,dt)
    CALL ScaleV(dt,RhsVecMet)
    CALL ScaleV(dt,RhsVecChem)
    DO jStage=1,iStage-1
      CALL Axpy(beta(iStage,jStage),VecIncrMetC(jStage)%Vec,RhsVecMet) 
      CALL Axpy(beta(iStage,jStage),VecIncrChemC(jStage)%Vec,RhsVecChem) 
    END DO
    Iter=BiCGStabMaxIter
    Tol=BiCGStabTol
    VecIncrMetC(iStage)%Vec=Zero
    VecIncrChemC(iStage)%Vec=Zero
    IF (Source) THEN
      CALL SolveSourceMetSp(RhsVecMet,JacTrans)
      CALL SolveSourceChemSp(RhsVecMet,JacTrans)
    END IF
    CALL Copy(RhsVecMet,VecIncrMetC(iStage)%Vec)
    CALL Copy(RhsVecChem,VecIncrChemC(iStage)%Vec)
    CALL ExchangeCell(VecIncrMetC(iStage)%Vec)
    CALL ExchangeCell(VecIncrChemC(iStage)%Vec)
    CALL Copy(VecMetOld,VecMet)
    CALL Copy(VecChemOld,VecChem)
    DO jStage=1,iStage
      CALL Axpy(alpha(iStage,jStage),VecIncrMetC(jStage)%Vec,VecMet)
      CALL Axpy(alpha(iStage,jStage),VecIncrChemC(jStage)%Vec,VecChem)
    END DO
  END DO

END SUBROUTINE StageRosParcel

SUBROUTINE InitExpIntRK(VecMet)

  TYPE (Vector4Cell_T), POINTER :: VecMet(:)
  INTEGER :: i,j

  IF (InitR==0) THEN
    Call MethodsRK
    WRITE(*,*) 'Method   ',Method
    IF (TRIM(Method(3:))=='RK1') THEN
      nStage=RK1%nStage
      nPhi=RK1%nPhi
      CFLNumber=RK1%CFLNumber
      aRunge=>RK1%aRunge
      cRunge=>RK1%cRunge
      dRunge=>RK1%dRunge
      gRunge=>RK1%gRunge
      dtRunge=>RK1%dtRunge
    ELSE IF (TRIM(Method(3:))=='RK2a') THEN
      nStage=RK2a%nStage
      nPhi=RK2a%nPhi
      aRunge=>RK2a%aRunge
      cRunge=>RK2a%cRunge
      dRunge=>RK2a%dRunge
      gRunge=>RK2a%gRunge
      dtRunge=>RK2a%dtRunge
    ELSE IF (TRIM(Method(3:))=='RK2Ros') THEN
      nStage=RK2Ros%nStage
      nPhi=RK2Ros%nPhi
      aRunge=>RK2Ros%aRunge
      cRunge=>RK2Ros%cRunge
      dRunge=>RK2Ros%dRunge
      gRunge=>RK2Ros%dRunge
      dtRunge=>RK2Ros%dtRunge
    ELSE IF (TRIM(Method(3:))=='RK2b') THEN
      nStage=RK2b%nStage
      nPhi=RK2b%nPhi
      aRunge=>RK2b%aRunge
      cRunge=>RK2b%cRunge
      dRunge=>RK2b%dRunge
      gRunge=>RK2b%gRunge
      dtRunge=>RK2b%dtRunge
    ELSE IF (TRIM(Method(3:))=='RK2bT') THEN
      nStage=RK2bT%nStage
      nPhi=RK2bT%nPhi
      aRunge=>RK2bT%aRunge
      cRunge=>RK2bT%cRunge
      dRunge=>RK2bT%dRunge
      gRunge=>RK2bT%gRunge
      dtRunge=>RK2bT%dtRunge
    ELSE IF (TRIM(Method(3:))=='RK3') THEN
      WRITE(*,*) 'Method RK3'
      nStage=RK3%nStage
      nPhi=RK3%nPhi
      CFLNumber=RK3%CFLNumber
      aRunge=>RK3%aRunge
      cRunge=>RK3%cRunge
      dRunge=>RK3%dRunge
      gRunge=>RK3%gRunge
      dtRunge=>RK3%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKJ') THEN
      nStage=RKJeb%nStage
      nPhi=RKJeb%nPhi
      CFLNumber=RKJeb%CFLNumber
      aRunge=>RKJeb%aRunge
      cRunge=>RKJeb%cRunge
      dRunge=>RKJeb%dRunge
      gRunge=>RKJeb%gRunge
      dtRunge=>RKJeb%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKJ1') THEN
      nStage=RKJeb1%nStage
      nPhi=RKJeb1%nPhi
      CFLNumber=RKJeb1%CFLNumber
      aRunge=>RKJeb1%aRunge
      cRunge=>RKJeb1%cRunge
      dRunge=>RKJeb1%dRunge
      gRunge=>RKJeb1%gRunge
      dtRunge=>RKJeb1%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKL') THEN
      nStage=RKLan%nStage
      nPhi=RKLan%nPhi
      aRunge=>RKLan%aRunge
      cRunge=>RKLan%cRunge
      dRunge=>RKLan%dRunge
      gRunge=>RKLan%gRunge
      dtRunge=>RKLan%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKNl') THEN
      WRITE(*,*) 'Method RKNl'
      nStage=RKNlopt%nStage
      nPhi=RKNlopt%nPhi
      aRunge=>RKNlopt%aRunge
      cRunge=>RKNlopt%cRunge
      dRunge=>RKNlopt%dRunge
      gRunge=>RKNlopt%gRunge
      dtRunge=>RKNlopt%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS4') THEN
      WRITE(*,*) 'Method MIS4'
      nStage=MIS4%nStage
      nPhi=MIS4%nPhi
      CFLNumber=MIS4%CFLNumber
      aRunge=>MIS4%aRunge
      cRunge=>MIS4%cRunge
      dRunge=>MIS4%dRunge
      gRunge=>MIS4%gRunge
      dtRunge=>MIS4%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS_SV3_3_4_13_A_fixed') THEN
      nStage=MIS_SV3_3_4_13_A_fixed%nStage
      nPhi=MIS_SV3_3_4_13_A_fixed%nPhi
      CFLNumber=MIS_SV3_3_4_13_A_fixed%CFLNumber
      aRunge=>MIS_SV3_3_4_13_A_fixed%aRunge
      cRunge=>MIS_SV3_3_4_13_A_fixed%cRunge
      dRunge=>MIS_SV3_3_4_13_A_fixed%dRunge
      gRunge=>MIS_SV3_3_4_13_A_fixed%gRunge
      dtRunge=>MIS_SV3_3_4_13_A_fixed%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS_SV2_2_3_10_A_fixed') THEN
      nStage=MIS_SV2_2_3_10_A_fixed%nStage
      nPhi=MIS_SV2_2_3_10_A_fixed%nPhi
      CFLNumber=MIS_SV2_2_3_10_A_fixed%CFLNumber
      WRITE(*,*) 'CFLNumber Set',CFLNumber
      aRunge=>MIS_SV2_2_3_10_A_fixed%aRunge
      cRunge=>MIS_SV2_2_3_10_A_fixed%cRunge
      dRunge=>MIS_SV2_2_3_10_A_fixed%dRunge
      gRunge=>MIS_SV2_2_3_10_A_fixed%gRunge
      dtRunge=>MIS_SV2_2_3_10_A_fixed%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS_SV2_2_4_14_A_fixed') THEN
      nStage=MIS_SV2_2_4_14_A_fixed%nStage
      nPhi=MIS_SV2_2_4_14_A_fixed%nPhi
      CFLNumber=MIS_SV2_2_4_14_A_fixed%CFLNumber
      aRunge=>MIS_SV2_2_4_14_A_fixed%aRunge
      cRunge=>MIS_SV2_2_4_14_A_fixed%cRunge
      dRunge=>MIS_SV2_2_4_14_A_fixed%dRunge
      gRunge=>MIS_SV2_2_4_14_A_fixed%gRunge
      dtRunge=>MIS_SV2_2_4_14_A_fixed%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS_SV2_2_5_19_A_fixed') THEN
      nStage=MIS_SV2_2_5_19_A_fixed%nStage
      nPhi=MIS_SV2_2_5_19_A_fixed%nPhi
      CFLNumber=MIS_SV2_2_5_19_A_fixed%CFLNumber
      aRunge=>MIS_SV2_2_5_19_A_fixed%aRunge
      cRunge=>MIS_SV2_2_5_19_A_fixed%cRunge
      dRunge=>MIS_SV2_2_5_19_A_fixed%dRunge
      gRunge=>MIS_SV2_2_5_19_A_fixed%gRunge
      dtRunge=>MIS_SV2_2_5_19_A_fixed%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS_SV3_3_5_17_A_fixed') THEN
      nStage=MIS_SV3_3_5_17_A_fixed%nStage
      nPhi=MIS_SV3_3_5_17_A_fixed%nPhi
      CFLNumber=MIS_SV3_3_5_17_A_fixed%CFLNumber
      WRITE(*,*) 'MIS_SV3_3_5_17_A_fixed%CFLNumber',CFLNumber
      aRunge=>MIS_SV3_3_5_17_A_fixed%aRunge
      cRunge=>MIS_SV3_3_5_17_A_fixed%cRunge
      dRunge=>MIS_SV3_3_5_17_A_fixed%dRunge
      gRunge=>MIS_SV3_3_5_17_A_fixed%gRunge
      dtRunge=>MIS_SV3_3_5_17_A_fixed%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS4_4') THEN
      WRITE(*,*) 'Method MIS4_4'
      nStage=MIS4_4%nStage
      nPhi=MIS4_4%nPhi
      CFLNumber=MIS4_4%CFLNumber
      aRunge=>MIS4_4%aRunge
      cRunge=>MIS4_4%cRunge
      dRunge=>MIS4_4%dRunge
      gRunge=>MIS4_4%gRunge
      dtRunge=>MIS4_4%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS2') THEN
      WRITE(*,*) 'Method MIS2'
      nStage=MIS2%nStage
      nPhi=MIS2%nPhi
      CFLNumber=MIS2%CFLNumber
      aRunge=>MIS2%aRunge
      cRunge=>MIS2%cRunge
      dRunge=>MIS2%dRunge
      gRunge=>MIS2%gRunge
      dtRunge=>MIS2%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKN3E4') THEN
      WRITE(*,*) 'Method RKN3E4'
      nStage=RKN3E4%nStage
      nPhi=RKN3E4%nPhi
      aRunge=>RKN3E4%aRunge
      cRunge=>RKN3E4%cRunge
      dRunge=>RKN3E4%dRunge
      gRunge=>RKN3E4%gRunge
      dtRunge=>RKN3E4%dtRunge
    ELSE IF (TRIM(Method(3:))=='MIS3C') THEN
      WRITE(*,*) 'Method MIS3C'
      nStage=MIS3C%nStage
      nPhi=MIS3C%nPhi
      aRunge=>MIS3C%aRunge
      cRunge=>MIS3C%cRunge
      dRunge=>MIS3C%dRunge
      gRunge=>MIS3C%gRunge
      dtRunge=>MIS3C%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKtvdA') THEN
      WRITE(*,*) 'Method RKtvdA'
      nStage=RKtvdA%nStage
      nPhi=RKtvdA%nPhi
      aRunge=>RKtvdA%aRunge
      cRunge=>RKtvdA%cRunge
      dRunge=>RKtvdA%dRunge
      gRunge=>RKtvdA%gRunge
      dtRunge=>RKtvdA%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKtvdB') THEN
      WRITE(*,*) 'Method RKtvdB'
      nStage=RKtvdB%nStage
      nPhi=RKtvdB%nPhi
      aRunge=>RKtvdB%aRunge
      cRunge=>RKtvdB%cRunge
      dRunge=>RKtvdB%dRunge
      gRunge=>RKtvdB%gRunge
      dtRunge=>RKtvdB%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKO') THEN
      nStage=RK3%nStage
      nPhi=RK3%nPhi
      aRunge=>RKW23%aRunge
      cRunge=>RKW23%cRunge
      dRunge=>RKW23%dRunge
      gRunge=>RKW23%gRunge
      dtRunge=>RKW23%dtRunge
    ELSE IF (TRIM(Method(3:))=='RK3TVDF') THEN
      nStage=RK3TVDF%nStage
      nPhi=RK3TVDF%nPhi
      aRunge=>RK3TVDF%aRunge
      cRunge=>RK3TVDF%cRunge
      dRunge=>RK3TVDF%dRunge
      gRunge=>RK3TVDF%gRunge
      dtRunge=>RK3TVDF%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKC') THEN
      nStage=RKCel%nStage
      nPhi=RKCel%nPhi
      aRunge=>RKCel%aRunge
      cRunge=>RKCel%cRunge
      dRunge=>RKCel%dRunge
      gRunge=>RKCel%gRunge
      dtRunge=>RKCel%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKI') THEN
      nStage=RKIMEX%nStage
      nPhi=RKIMEX%nPhi
      aRunge=>RKIMEX%aRunge
      cRunge=>RKIMEX%cRunge
      dRunge=>RKIMEX%dRunge
      gRunge=>RKIMEX%gRunge
      dtRunge=>RKIMEX%dtRunge
    ELSE IF (TRIM(Method(3:))=='RKW') THEN
      nStage=RKWen%nStage
      nPhi=RKWen%nPhi
      aRunge=>RKWen%aRunge
      cRunge=>RKWen%cRunge
      dRunge=>RKWen%dRunge
      gRunge=>RKWen%gRunge
      dtRunge=>RKWen%dtRunge
    ELSE IF (TRIM(Method(3:))=='RK4') THEN
      nStage=RK4%nStage
      nPhi=RK4%nPhi
      aRunge=>RK4%aRunge
      cRunge=>RK4%cRunge
      dRunge=>RK4%dRunge
      gRunge=>RK43K%gRunge
      dtRunge=>RK43K%dtRunge
    ELSE IF (TRIM(Method(3:))=='RK43K') THEN
      nStage=RK43K%nStage
      nPhi=RK43K%nPhi
      aRunge=>RK43K%aRunge
      cRunge=>RK43K%cRunge
      dRunge=>RK43K%dRunge
      gRunge=>RK43K%gRunge
      dtRunge=>RK43K%dtRunge
    ELSE IF (TRIM(Method(3:))=='EB4') THEN
      nStage=RKEB4%nStage
      nPhi=RKEB4%nPhi
      aRunge=>RKEB4%aRunge
      cRunge=>RKEB4%cRunge
      dRunge=>RKEB4%dRunge
      gRunge=>RKEB4%gRunge
      dtRunge=>RKEB4%dtRunge
    ELSE
      WRITE(*,*) 'Falsche Methode'
      WRITE(*,*) Method
      STOP
    END IF
    IF (MyId==0) THEN
      WRITE(*,*) Method
      DO i=1,SIZE(aRunge,1)
        DO j=1,i-1
          WRITE(*,*) 'aRunge',i,j,aRunge(i,j)%Koeff(0)
        END DO  
      END DO  
    END IF  
    CALL Allocate(RhsVecS,VecMet,0,VectorComponentsM)
    CALL Allocate(RhsVecF,VecMet,1,VectorComponentsM-6)
    ALLOCATE(VecVelC(nStage+1))
    DO i=1,nStage+1
      CALL Allocate(VecVelC(i)%Vec,VecMet,0,VectorComponentsM-6)
    END DO
    CALL Allocate(RhsVecFF)
    ALLOCATE(VecRhsVecFS(nStage))
    DO i=1,nStage
      CALL Allocate(VecRhsVecFS(i)%VecF)
    END DO
    ALLOCATE(VecRhsVecCS(nStage))
    DO i=1,nStage
      CALL Allocate(VecRhsVecCS(i)%Vec,VecMet,1,VectorComponentsM-6)
    END DO
    ALLOCATE(VecVelF(nStage+1))
    DO i=1,nStage+1
      CALL Allocate(VecVelF(i)%VecF)
    END DO

    ALLOCATE(VecMethetaF(nStage))
    DO i=1,nStage
      CALL Allocate(VecMethetaF(i)%VecF,VecMet,1,VectorComponentsM-6)
    END DO
    CALL Allocate(ThetaF,VecMet,1,VectorComponentsM-6)

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

    CALL Allocate(VecMetOld,VecMet)
    CALL Allocate(VelFOld)

    CALL Allocate(VelFFast)
    CALL Allocate(RhsVelFS)
    CALL Allocate(RhsCS,VecMet,1,VectorComponentsM-6)
    CALL Allocate(VecU,uPosL,wPosR)
    InitR=1
  END IF

END SUBROUTINE InitExpIntRK

SUBROUTINE ExpIntRk(VelF,VelC,dtAct,Time,ATol,RTol)
  
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  TYPE (Vector4Cell_T), POINTER :: VelC(:)
  REAL(RealKind) :: dtAct,Time,ATol(:),RTol(:)
  
  INTEGER :: i,is,nsLoc,jStage
  REAL(RealKind) :: dTau,dtLoc,TimeFast
  REAL(RealKind) :: Temp

  dTau=dtAct/ns
  DO iStage=1,nStage+1
    CALL Copy(VelF,VecVelF(iStage)%VecF)
    CALL Copy(VelC,VecVelC(iStage)%Vec)
  END DO
  DO iStage=1,nStage
    TimeAct=Time+cRunge(iStage)*dtAct
!   Computation of the initial value 
    DO jStage=1,iStage
      CALL Ax1mx2py(dRunge(iStage+1,jStage),VecVelF(jStage)%VecF,VelF,VecVelF(iStage+1)%VecF)
      CALL Ax1mx2py(dRunge(iStage+1,jStage),VecVelC(jStage)%Vec,VelC,VecVelC(iStage+1)%Vec) 
    END DO
    CALL BoundaryVelocity(VecVelF(iStage)%VecF,Time+cRunge(iStage)*dtAct)
    CALL VelocityFaceToCellLR(VecVelF(iStage)%VecF,VecU)
    CALL BoundaryVelocity(VecU,Time+cRunge(iStage)*dtAct)
    CALL ExchangeCell(VecU)
    CALL BoundaryCondition(VecVelC(iStage)%Vec,VelF=VecVelF(iStage)%VecF,Time=Time+cRunge(iStage)*dtAct)
!----------------------------------------
!   Part of PrepareF
    CALL PrepareFEx(VecVelC(iStage)%Vec,VecVelF(iStage)%VecF,VecU,Time+cRunge(iStage)*dtAct)
!   Part of PrepareF
!----------------------------------------
    RhsVecS=Zero
    CALL FcnMetSlow(VecU,VecVelC(iStage)%Vec,VecVelF(iStage)%VecF,RhsVecS &
                   ,VecMethetaF(iStage)%VecF,VecPreFacF(iStage)%VecF,VecSoundFac(iStage)%Vec &
                   ,Time+cRunge(iStage)*dtAct,dt)
    CALL ExchangeCell(RhsVecS)
    CALL VelocityCellToFaceLR(RhsVecS,VecRhsVecFS(iStage)%VecF,VelF,Time+cRunge(iStage)*dtAct)
    CALL Copy(RhsVecS,VecRhsVecCS(iStage)%Vec)
    nsLoc=CEILING(ns*dtRunge(iStage+1))
    dtLoc=dtAct*dtRunge(iStage+1)
    dTau=dtLoc/nsLoc
    TimeFast=Time+cRunge(iStage+1)*dtAct-dtLoc
    IF (nsLoc>0) THEN
      SELECT CASE(IntFast)
        CASE ('FB')
          CALL ForBack1(dTau,nsLoc,VecVelF(iStage+1)%VecF,VecVelC(iStage+1)%Vec,VecVelF,VecVelC,VecMethetaF,VecPreFacF &
                       ,VecSoundFac,VecRhsVecFS,VecRhsVecCS,TimeFast,dtAct)
        CASE ('Verlet')
          CALL StoermerVerlet(dTau,nsLoc,VecVelF(iStage+1)%VecF,VecVelC(iStage+1)%Vec,VecVelF,VecVelC,VecMethetaF,VecPreFacF &
                       ,VecSoundFac,VecRhsVecFS,VecRhsVecCS,TimeFast,dtAct)
      END SELECT    
      CALL ExchangeCell(VecVelC(iStage+1)%Vec) !OSSI
    ELSE
!     CALL ConstInt(VecVelF(iStage+1)%VecF,VecVelC(iStage+1)%Vec,VecRhsVecFS,VecRhsVecS,dtAct)
    END IF
  END DO
  CALL Copy(VecVelC(nStage+1)%Vec,VelC)
  CALL Copy(VecVelF(nStage+1)%VecF,VelF)
  
  Time=Time+dtAct

END SUBROUTINE ExpIntRk

SUBROUTINE ConstInt(VelF,VelC,VecRhsVecFS,VecRhsVecS,dT)

  REAL(RealKind) :: dT
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  TYPE (Vector4Cell_T), POINTER :: VelC(:)
  TYPE (VecVelocityFace_T), POINTER :: VecRhsVecFS(:)
  TYPE (VecVector4Cell_T), POINTER :: VecRhsVecS(:)

  INTEGER :: iPhi,is,jStage
  REAL(RealKind) :: Tau,Fac
  REAL(RealKind) :: Summe

  DO jStage=1,iStage
    Fac=aRunge(iStage+1,jStage)%Koeff(0)*dT
    CALL Axpy(Fac,VecRhsVecFS(jStage)%VecF,VelF)
    CALL Axpy(Fac,VecRhsVecS(jStage)%Vec,VelC)
  END DO
  CALL ExchangeCell(VelC)

END SUBROUTINE ConstInt

SUBROUTINE ForBack1(dTau,ns,VelF,VelC, &
                    VecVelFS,VecVelCS, &
                    VecMethetaF,         &
                    VecPreFacF,         &
                    VecSoundFac,         &
                    VecRhsVecFS,        &
                    VecRhsVecCS,Time,dT)

  REAL(RealKind) :: dTau,Time,dT
  INTEGER :: ns
  TYPE (VelocityFace_T), POINTER :: VelF(:) ! Stufen für Geschwindikeit
  TYPE (Vector4Cell_T), POINTER :: VelC(:)  ! Stufen für Skalare
  TYPE (VecVelocityFace_T), POINTER :: VecVelFS(:)
  TYPE (VecVector4Cell_T), POINTER :: VecVelCS(:)
  TYPE (VecVectorSFace_T), POINTER :: VecMethetaF(:)
  TYPE (VecVelocityFace_T), POINTER :: VecPreFacF(:)
  TYPE (VecScalarCell_T), POINTER :: VecSoundFac(:)
  TYPE (VecVelocityFace_T), POINTER :: VecRhsVecFS(:)
  TYPE (VecVector4Cell_T), POINTER :: VecRhsVecCS(:)

  INTEGER :: iPhi,is,jStage
  REAL(RealKind) :: Tau,Fac

  RhsVelFS=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecRhsVecFS(jStage)%VecF,RhsVelFS)
    CALL Ax1mx2py(gRunge(iStage+1,jStage)/dT,VecVelFS(jStage)%VecF,VecVelFS(1)%VecF,RhsVelFS)
  END DO
  ThetaF=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecMethetaF(jStage)%VecF,ThetaF)
  END DO
  PreFacF=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecPreFacF(jStage)%VecF,PreFacF)
  END DO

  SoundFac=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecSoundFac(jStage)%Vec,SoundFac)
  END DO

  RhsCS=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecRhsVecCS(jStage)%Vec,RhsCS)
    CALL Ax1mx2py(gRunge(iStage+1,jStage)/dT,VecVelCS(jStage)%Vec,VecVelCS(1)%Vec,RhsCS)
  END DO
  IF (ns>0) THEN
    DO is=1,ns
!     CALL Copy(VelC,VecMetOld)
      CALL Copy(VelF,VelFOld)
      RhsVecFF=Zero
      CALL BoundaryCondition(VelC,VelF=VelF,Time=Time)
      CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
      CALL Axpy(dTau,RhsVecFF,VelF)
      CALL Axpy(dTau,RhsVelFS,VelF)
      RhsVecF=Zero
      CALL Axpby((1.0d0+GammaDiv),VelF,-GammaDiv,VelFOld)
      CALL BoundaryVelocity(VelFOld,Time+dTau)
      CALL FcnMetFastScalar(VelC,VelFOld,ThetaF,SoundFac,RhsVecF,Time)
!     CALL BoundaryVelocity(VelF,Time+dTau)
!     CALL FcnMetFastScalar(VelC,VelF,ThetaF,SoundFac,RhsVecF,Time)
      CALL Ax1px2py(dTau,RhsVecF,RhsCS,VelC)
      Time=Time+dTau
    END DO
  ELSE
    RhsVecFF=Zero
    CALL BoundaryCondition(VelC,VelF=VelF,Time=Time)
    CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
    CALL Axpy(dTau,RhsVecFF,VelF)
    CALL Axpy(dTau,RhsVelFS,VelF)
    RhsVecF=Zero
    CALL BoundaryVelocity(VelF,Time)
    CALL FcnMetFastScalar(VelC,VelF,ThetaF,SoundFac,RhsVecF,Time)
    CALL Ax1px2py(dTau,RhsVecF,RhsCS,VelC)
  END IF
END SUBROUTINE ForBack1

SUBROUTINE StoermerVerlet(dTau,ns,VelF,VelC, &
                          VecVelFS,VecVelCS, &
                          VecMethetaF,         &
                          VecPreFacF,         &
                          VecSoundFac,         &
                          VecRhsVecFS,        &
                          VecRhsVecCS,Time,dT)

  REAL(RealKind) :: dTau,Time,dT
  INTEGER :: ns
  TYPE (VelocityFace_T), POINTER :: VelF(:) ! Stufen für Geschwindikeit
  TYPE (Vector4Cell_T), POINTER :: VelC(:)  ! Stufen für Skalare
  TYPE (VecVelocityFace_T), POINTER :: VecVelFS(:)
  TYPE (VecVector4Cell_T), POINTER :: VecVelCS(:)
  TYPE (VecVectorSFace_T), POINTER :: VecMethetaF(:)
  TYPE (VecVelocityFace_T), POINTER :: VecPreFacF(:)
  TYPE (VecScalarCell_T), POINTER :: VecSoundFac(:)
  TYPE (VecVelocityFace_T), POINTER :: VecRhsVecFS(:)
  TYPE (VecVector4Cell_T), POINTER :: VecRhsVecCS(:)

  INTEGER :: iPhi,is,jStage
  REAL(RealKind) :: Tau,Fac
  REAL(RealKind) :: Summe

  RhsVelFS=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecRhsVecFS(jStage)%VecF,RhsVelFS)
    CALL Ax1mx2py(gRunge(iStage+1,jStage)/dT,VecVelFS(jStage)%VecF,VecVelFS(1)%VecF,RhsVelFS)
  END DO
  ThetaF=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecMethetaF(jStage)%VecF,ThetaF)
  END DO
  PreFacF=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecPreFacF(jStage)%VecF,PreFacF)
  END DO

  SoundFac=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecSoundFac(jStage)%Vec,SoundFac)
  END DO

  RhsCS=Zero
  DO jStage=1,iStage
    CALL Axpy(aRunge(iStage+1,jStage)%Koeff(0),VecRhsVecCS(jStage)%Vec,RhsCS)
    CALL Ax1mx2py(gRunge(iStage+1,jStage)/dT,VecVelCS(jStage)%Vec,VecVelCS(1)%Vec,RhsCS)
  END DO
  IF (ns>0) THEN
    RhsVecFF=Zero
    CALL BoundaryCondition(VelC,VelF=VelF,Time=Time)
    CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
    CALL Axpy(Half*dTau,RhsVecFF,VelF)
    CALL Axpy(Half*dTau,RhsVelFS,VelF)
    DO is=1,ns-1
      RhsVecF=Zero
      CALL BoundaryVelocity(VelF,Time+Half*dTau)
      CALL FcnMetFastScalar(VelC,VelF,ThetaF,SoundFac,RhsVecF,Time)
      CALL Ax1px2py(dTau,RhsVecF,RhsCS,VelC)
      RhsVecFF=Zero
      CALL BoundaryCondition(VelC,VelF=VelF,Time=Time)
      CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
      CALL Axpy(dTau,RhsVecFF,VelF)
      CALL Axpy(dTau,RhsVelFS,VelF)
      Time=Time+dTau
    END DO
    RhsVecF=Zero
    CALL BoundaryVelocity(VelF,Time-Half*dTau)
    CALL FcnMetFastScalar(VelC,VelF,ThetaF,SoundFac,RhsVecF,Time)
    CALL Ax1px2py(dTau,RhsVecF,RhsCS,VelC)
    RhsVecFF=Zero
    CALL BoundaryCondition(VelC,VelF=VelF,Time=Time)
    CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
    CALL Axpy(Half*dTau,RhsVecFF,VelF)
    CALL Axpy(Half*dTau,RhsVelFS,VelF)
  ELSE
    RhsVecFF=Zero
    CALL BoundaryCondition(VelC,VelF=VelF,Time=Time)
    CALL FcnMetFastU(VelC,VelF,PreFacF,RhsVecFF,Time)
    CALL Axpy(dTau,RhsVecFF,VelF)
    CALL Axpy(dTau,RhsVelFS,VelF)
    RhsVecF=Zero
    CALL BoundaryVelocity(VelF,Time+dTau)
    CALL FcnMetFastScalar(VelC,VelF,ThetaF,SoundFac,RhsVecF,Time)
    CALL Ax1px2py(dTau,RhsVecF,RhsCS,VelC)
  END IF
END SUBROUTINE StoermerVerlet

FUNCTION Horner(Koeff,tau)

  REAL(RealKind) :: Horner
  REAL(RealKind) :: Koeff(:)
  REAL(RealKind) :: tau

  INTEGER :: i,n

  n=SIZE(Koeff)
  Horner=Koeff(n)
  DO i=n-1,1,-1
    Horner=Horner*Tau+Koeff(i)
  END DO

END FUNCTION Horner

SUBROUTINE Ros3MetC(VecMet,VecChem,VelF,VecG,dtAct,Time,ATol,RTol)

  TYPE (Vector4Cell_T), POINTER :: VecMet(:)
  TYPE (Vector4Cell_T), POINTER :: VecChem(:)
  TYPE (VelocityFace_T), POINTER :: VelF(:)
  TYPE (Vector4Cell_T), POINTER :: VecG(:)
  REAL(RealKind) :: dtAct,Time,ATol(:),RTol(:)

  INTEGER :: iTer1
  REAL(RealKind) :: ErrLoc
  REAL(RealKind) :: dt1,dt2
  INTEGER :: SignCheck,SignCheckLoc
  INTEGER :: ix,iy,iz

  CALL InitRos3MetC(VecMet,VecChem)
  CALL Copy(VelF,VelFOld)
  CALL Copy(VecMet,VecMetOld)
  CALL Copy(VecChem,VecChemOld)
  IF (Anelastic.OR.PseudoIn) THEN
    CALL Copy(VecG,VelGOld)
  END IF

  iTer1=0
  DO
    dt=dtAct
    CALL StageRos(VecMet,VecChem,VelF,VecG,dt,Time)
!   Error estimation
    IF (ErrControl) THEN
      CALL Xpay(VecIncrMetC(2)%Vec,Delta,VecIncrMetC(1)%Vec)
      CALL Xpay(VecIncrF(2)%VecF,Delta,VecIncrF(1)%VecF)
      IF (uPosL>0) THEN
        CALL VelocityFaceToCellLR(VecIncrF(1)%VecF,VecIncrMetC(1)%Vec)
        CALL VelocityFaceToCellLR(VecIncrF(2)%VecF,VecIncrMetC(2)%Vec)
      END IF
      ErrLoc=Error(VecIncrMetC(1)%Vec,VecMet,VecMetOld,ATol,RTol,VecIncrMetC(2)%Vec)
!     ErrLoc=ErrorMax(VecIncrC(1)%Vec,VelC,VecMetOld,ATol,RTol,VecIncrC(2)%Vec)
      IF (MyId==0) THEN
        WRITE(TermUnit,*) ' ErrLoc',ErrLoc,dtAct
        WRITE(*,*) ' ErrLoc',ErrLoc,dtAct
      END IF
      IF (ErrLoc>One) THEN
        dtAct=MAX(0.5e0_RealKind,MIN(2.0e0_RealKind,0.8e0_RealKind/SQRT(ErrLoc)))*dtAct
        dtAct=MIN(dtAct,dtMax)
        CALL Copy(VelFOld,VelF)
        CALL Copy(VecMetOld,VecMet)
        iTer1=iTer1+1
      ELSE
        IF (MyId==0) THEN
          WRITE(110,*) Time,dtAct
        END IF
        IF (iTer1==0) THEN
          dtAct=MAX(0.5e0_RealKind,MIN(2.0e0_RealKind,0.8e0_RealKind/SQRT(ErrLoc)))*dtAct
          dtAct=MIN(dtAct,dtMax)
        END IF
        EXIT
      END IF
    ELSE
      EXIT
    END IF
  END DO
  !CALL BoundaryVelocity(VelF,Time+dt) OSSI
  Time=Time+dt
  dtAct=MIN(EndTime-Time,dtAct)*(1.d0+1.d-13)
END SUBROUTINE Ros3MetC

SUBROUTINE RosParcel(VelF,VecMet,VecChem,dtAct,Time,ATol,RTol)

  TYPE (VelocityFace_T), POINTER :: VelF(:)
  TYPE (Vector4Cell_T), POINTER :: VecMet(:)
  TYPE (Vector4Cell_T), POINTER :: VecChem(:)
  REAL(RealKind) :: dtAct,Time,ATol(:),RTol(:)

  INTEGER :: iTer1
  REAL(RealKind) :: ErrLoc
  REAL(RealKind) :: dt1,dt2
  INTEGER :: SignCheck,SignCheckLoc
  INTEGER :: ix,iy,iz

  CALL InitRos3MetC(VecMet,VecChem)
  CALL Copy(VecMet,VecMetOld)
  CALL Copy(VecChem,VecChemOld)

  iTer1=0
  DO
    dt=dtAct
    CALL StageRosParcel(VecMet,VecChem,VelF,dt,Time)

    IF (ErrControl) THEN
      CALL Xpay(VecIncrMetC(2)%Vec,Delta,VecIncrMetC(1)%Vec)
      ErrLoc=Error(VecIncrMetC(1)%Vec,VecMet,VecMetOld,ATol,RTol,VecIncrMetC(2)%Vec)
      IF (MyId==0) THEN
        WRITE(TermUnit,*) ' ErrLoc',ErrLoc,dtAct
        WRITE(*,*) ' ErrLoc',ErrLoc,dtAct
      END IF
      IF (ErrLoc>One) THEN
        dtAct=MAX(0.5e0_RealKind,MIN(2.0e0_RealKind,0.8e0_RealKind/SQRT(ErrLoc)))*dtAct
        dtAct=MIN(dtAct,dtMax)
        CALL Copy(VecMetOld,VecMet)
        iTer1=iTer1+1
      ELSE
        IF (MyId==0) THEN
          WRITE(110,*) Time,dtAct
        END IF
        IF (iTer1==0) THEN
          dtAct=MAX(0.5e0_RealKind,MIN(2.0e0_RealKind,0.8e0_RealKind/SQRT(ErrLoc)))*dtAct
          dtAct=MIN(dtAct,dtMax)
        END IF
        EXIT
      END IF
    ELSE
      EXIT
    END IF
  END DO
  Time=Time+dt
  dtAct=MIN(EndTime-Time,dtAct)*(1.d0+1.d-13)
END SUBROUTINE RosParcel

END MODULE Int_Mod
