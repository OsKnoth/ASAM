MODULE Turbulence_Mod
  USE DataType_Mod
  USE Domain_Mod
  USE Floor_Mod
  USE Physics_Mod
  USE Operator_Mod
  USE Parameter_mod
  USE Output_Mod
  USE Emission_Mod

  IMPLICIT NONE

  INTEGER, PRIVATE :: i,j,k
  INTEGER, PRIVATE :: ic,icE
  INTEGER, PRIVATE :: uPosJac,vPosJac,wPosJac
  REAL(RealKind), PRIVATE :: LenScale
  REAL(RealKind), PRIVATE, POINTER :: c(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: th(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: tke(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: dis(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: ome(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: tkeH(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: tkeV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: LenC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: qv(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: qc(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: qi(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: qr(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Rho(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: D(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DH(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DPot(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DPotH(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DPotV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DMom(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DMomH(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: DMomV(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wF(:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wC(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCL(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wCR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: f(:,:,:,:),fL(:,:,:,:),fR(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: tkeRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: disRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: omeRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: tkeHRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: tkeVRhs(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: LenRhs(:,:,:,:)
  TYPE(Vec4_T), PRIVATE, POINTER :: AS(:)
  TYPE(Vector4Cell_T), POINTER, PRIVATE :: VecAve(:)
  TYPE(Vector4Cell_T), POINTER, PRIVATE :: VecAve2(:)

  REAL(RealKind), PRIVATE, POINTER :: uHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uuHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vvHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wwHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uvHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uwHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vwHat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau11Hat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau22Hat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau33Hat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau12Hat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau13Hat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau23Hat(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uuHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vvHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: wwHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uvHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: uwHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: vwHat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau11Hat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau22Hat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau33Hat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau12Hat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau13Hat2(:,:,:,:)
  REAL(RealKind), PRIVATE, POINTER :: Tau23Hat2(:,:,:,:)

  INTEGER, PRIVATE :: InitDynSmag=0
  INTEGER, PARAMETER :: uHatP=1
  INTEGER, PARAMETER :: vHatP=2
  INTEGER, PARAMETER :: wHatP=3
  INTEGER, PARAMETER :: uuHatP=4
  INTEGER, PARAMETER :: vvHatP=5
  INTEGER, PARAMETER :: wwHatP=6
  INTEGER, PARAMETER :: uvHatP=7
  INTEGER, PARAMETER :: uwHatP=8
  INTEGER, PARAMETER :: vwHatP=9
  INTEGER, PARAMETER :: Tau11HatP=10
  INTEGER, PARAMETER :: Tau22HatP=11
  INTEGER, PARAMETER :: Tau33HatP=12
  INTEGER, PARAMETER :: Tau12HatP=13
  INTEGER, PARAMETER :: Tau13HatP=14
  INTEGER, PARAMETER :: Tau23HatP=15

CONTAINS

SUBROUTINE Turbulence(Vector,Rhs,UVec)

  TYPE(Vector4Cell_T) :: Vector,Rhs
  TYPE (Vector4Cell_T), OPTIONAL :: UVec

  IF (TkeDis) THEN
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tke=>Vector%Vec(tkePos)%c
    dis=>Vector%Vec(disPos)%c
    tkeRhs=>Rhs%Vec(tkePos)%c
    disRhs=>Rhs%Vec(disPos)%c
    Rho=>Vector%Vec(RhoPos)%c
    D=>DiffKoeff(ibLoc)%c
    th=>Vector%Vec(thPos)%c
    CALL TkeDisWallCompute
    CALL TkeDisCompute
    IF (Baum) THEN
      CALL TkeDisBaumCompute
    END IF  
  ELSE IF (TkeOme) THEN
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tke=>Vector%Vec(tkePos)%c
    ome=>Vector%Vec(omePos)%c
    tkeRhs=>Rhs%Vec(tkePos)%c
    omeRhs=>Rhs%Vec(omePos)%c
    Rho=>Vector%Vec(RhoPos)%c
    D=>DiffKoeff(ibLoc)%c
    CALL TkeOmeWallCompute
    CALL TkeOmeCompute
    ! IF (Baum) CALL TkeDisBaumCompute(ibLoc) 
  ELSE IF (TkeDisRich) THEN
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tke=>Vector%Vec(tkePos)%c
    dis=>Vector%Vec(disPos)%c
    tkeRhs=>Rhs%Vec(tkePos)%c
    disRhs=>Rhs%Vec(disPos)%c
    Rho=>Vector%Vec(RhoPos)%c
    D=>DiffKoeff(ibLoc)%c
    CALL TkeDisWallCompute
    th=>Vector%Vec(thPos)%c
    CALL TkeDisRichCompute
    IF (Baum) THEN
      CALL TkeDisBaumCompute
    END IF  
  ELSE IF (TkeSGS) THEN
    IF (PRESENT(UVec)) THEN
      uCL=>UVec%Vec(uPosL)%c
      vCL=>UVec%Vec(vPosL)%c
      wCL=>UVec%Vec(wPosL)%c
      uCR=>UVec%Vec(uPosR)%c
      vCR=>UVec%Vec(vPosR)%c
      wCR=>UVec%Vec(wPosR)%c
    ELSE  
      uCL=>Vector%Vec(uPosL)%c
      vCL=>Vector%Vec(vPosL)%c
      wCL=>Vector%Vec(wPosL)%c
      uCR=>Vector%Vec(uPosR)%c
      vCR=>Vector%Vec(vPosR)%c
      wCR=>Vector%Vec(wPosR)%c
    END IF  
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tke=>Vector%Vec(tkePos)%c
    tkeRhs=>Rhs%Vec(tkePos)%c
    th=>Vector%Vec(thPos)%c
    Rho=>Vector%Vec(RhoPos)%c
    DPot=>DiffPotKoeff(ibLoc)%c
    DMom=>DiffMomKoeff(ibLoc)%c
    LenC=>LenKoeff(ibLoc)%c
    CALL TkeSGSCompute
  ELSE IF (TkeLen) THEN
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    th=>Vector%Vec(thPos)%c
    Rho=>Vector%Vec(RhoPos)%c
    tke=>Vector%Vec(tkePos)%c
    LenC=>Vector%Vec(LenPos)%c
    tkeRhs=>Rhs%Vec(tkePos)%c
    LenRhs=>Rhs%Vec(LenPos)%c
    DPot=>DiffPotKoeff(ibLoc)%c
    DMom=>DiffMomKoeff(ibLoc)%c
    CALL TkeLenCompute
  ELSE IF (TkeHVLen) THEN
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tkeH=>Vector%Vec(tkeHPos)%c
    tkeV=>Vector%Vec(tkeVPos)%c
    LenC=>Vector%Vec(LenPos)%c
    tkeHRhs=>Rhs%Vec(tkeHPos)%c
    tkeVRhs=>Rhs%Vec(tkeVPos)%c
    LenRhs=>Rhs%Vec(LenPos)%c
    DH=>DiffHKoeff(ibLoc)%c
    DV=>DiffVKoeff(ibLoc)%c
    th=>Vector%Vec(thPos)%c
    Rho=>Vector%Vec(RhoPos)%c
    CALL TkeHVLenCompute
  ENDIF
END SUBROUTINE Turbulence

SUBROUTINE JacTurbulence(Vector,Jac)

  TYPE(Vector4Cell_T) :: Vector
  TYPE(Vec4_T), POINTER :: Jac(:)

  IF (TkeDis) THEN
    AS=>Jac
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tke=>Vector%Vec(tkePos)%c
    dis=>Vector%Vec(disPos)%c
    D=>DiffKoeff(ibLoc)%c
    Rho=>Vector%Vec(RhoPos)%c
    th=>Vector%Vec(thPos)%c
    CALL JacTkeDisWallCompute
    CALL JacTkeDisCompute
    IF (Baum) THEN
      CALL JacTkeDisBaumCompute
    END IF  
  ELSE IF (TkeOme) THEN
    AS=>Jac
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tke=>Vector%Vec(tkePos)%c
    ome=>Vector%Vec(omePos)%c
    D=>DiffKoeff(ibLoc)%c
    Rho=>Vector%Vec(RhoPos)%c
    CALL JacTkeOmeWallCompute
    CALL JacTkeOmeCompute
    ! IF (Baum) CALL JacTkeDisBaumCompute(ibLoc) 
  ELSE IF (TkeDisRich) THEN
    AS=>Jac
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tke=>Vector%Vec(tkePos)%c
    dis=>Vector%Vec(disPos)%c
    D=>DiffKoeff(ibLoc)%c
    Rho=>Vector%Vec(RhoPos)%c
    CALL JacTkeDisWallCompute
    th=>Vector%Vec(thPos)%c
    CALL JacTkeDisRichCompute
    IF (Baum) THEN
      CALL JacTkeDisBaumCompute
    END IF  
  ELSE IF (TkeSGS) THEN
    AS=>Jac
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    th=>Vector%Vec(thPos)%c
    tke=>Vector%Vec(tkePos)%c
    th=>Vector%Vec(thPos)%c
    Rho=>Vector%Vec(RhoPos)%c
    DPot=>DiffPotKoeff(ibLoc)%c
    DMom=>DiffMomKoeff(ibLoc)%c
    LenC=>LenKoeff(ibLoc)%c
    CALL JacTkeSGSCompute
  ELSE IF (TkeLen) THEN
    AS=>Jac
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tke=>Vector%Vec(tkePos)%c
    LenC=>Vector%Vec(LenPos)%c
    DPot=>DiffPotKoeff(ibLoc)%c
    DMom=>DiffMomKoeff(ibLoc)%c
    th=>Vector%Vec(thPos)%c
    Rho=>Vector%Vec(RhoPos)%c
    CALL JacTkeLenCompute
  ELSE IF (TkeHVLen) THEN
    AS=>Jac
    uCL=>Vector%Vec(uPosL)%c
    vCL=>Vector%Vec(vPosL)%c
    wCL=>Vector%Vec(wPosL)%c
    uCR=>Vector%Vec(uPosR)%c
    vCR=>Vector%Vec(vPosR)%c
    wCR=>Vector%Vec(wPosR)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    tkeH=>Vector%Vec(tkeHPos)%c
    tkeV=>Vector%Vec(tkeVPos)%c
    LenC=>Vector%Vec(LenPos)%c
    DH=>DiffHKoeff(ibLoc)%c
    DV=>DiffVKoeff(ibLoc)%c
    th=>Vector%Vec(thPos)%c
    Rho=>Vector%Vec(RhoPos)%c
    CALL JacTkeHVLenCompute
  END IF

END SUBROUTINE JacTurbulence

SUBROUTINE TkeDisCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: tkeC,disC
  REAL(RealKind) :: S,O,dudz,uStar,z0,z,dzLoc

  INTEGER :: i,j,nrcl
  REAL(RealKind) :: dthdz,thC,rhoC
  REAL(RealKind) :: TempU,TempV,TempW ! velocities
  REAL(RealKind) :: TempVM ! mean velocity
  REAL(RealKind) :: FB,FT,VB,VCe,VT
  REAL(RealKind) :: Sd ! the dissipation of TKE by the interaction of the flow with the canopy elements
  REAL(RealKind) :: N2,B ! Brunt-Vasala-Frequency and buoyance production of TKE induced by canopy
  REAL(RealKind), POINTER :: LAD(:)
  REAL(RealKind), PARAMETER :: Prt=0.85d0
  REAL(RealKind), PARAMETER :: Cd=0.15d0 ! drag coefficient of canopy
  REAL(RealKind) :: FL,Frac
  REAL(RealKind) :: vTraffic=30000.0d0/3600.0d0    ! [m/s]
  REAL(RealKind) :: DensTraffic=600.0d0/3600.0d0  ! [ /s]
  REAL(RealKind) :: cDTraffic=0.45d0    ! []

! Quellen und Senken der Tke-Dis-Gleichungen
! Zentrale Gradienten linear interpoliert zwischen den prim�en Gradienten

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),        &
                         FV(ix:ix,iy-1:iy,iz:iz),        &
                         FW(ix:ix,iy:iy,iz-1:iz),        &
                         S,O)

        tkeC = tke(ix,iy,iz,1)
        disC = dis(ix,iy,iz,1)
        tkeRhs(ix,iy,iz,1) = tkeRhs(ix,iy,iz,1)  &           !!!FILAUS
                             +D(ix,iy,iz,1)*S*O-Cmy0*ABS(tkeC)*tkeC/(D(ix,iy,iz,1)+Eps)
        disRhs(ix,iy,iz,1) = disRhs(ix,iy,iz,1)  &
                             +Cmy1*Cmy0*tkeC*S*O-Cmy0*Cmy2 & 
                             *SQRT(D(ix,iy,iz,1)*ABS(disC)/Cmy0)*disC/(D(ix,iy,iz,1)+Eps)
      END DO
    END DO
  END DO

  IF (Traffic) THEN
    DO i=1,NumberOfEmiDomainLoc
      DO iCell=1,SIZE(EmiDomain(i)%BlockCell(ibLoc)%Cell)
        ix=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%ix
        iy=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%iy
        iz=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%iz
        Frac=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%Frac
        FL=EmiDomain(i)%BlockCell(ibLoc)%Cell(iCell)%FL
        tkeRhs(ix,iy,iz,1)=tkeRhs(ix,iy,iz,1)+ &
           Frac*(0.5d0*cDTraffic*vTraffic*vTraffic)*DensTraffic
      END DO                         
    END DO   
  END IF

  IF (Canopy) THEN
    DO i=1,NumBoundCell
      ix = BoundCell(i)%ix
      iy = BoundCell(i)%iy
      iz = BoundCell(i)%iz
      nrcl = BoundCell(i)%CanopyCell%NrCanopyLayers
      LAD => BoundCell(i)%CanopyCell%LAD

      DO j=iz, iz+nrcl-1
        CALL DeformVortN(uC(ix-1:ix+1,iy-1:iy+1,j-1:j+1,1),  &
                         vC(ix-1:ix+1,iy-1:iy+1,j-1:j+1,1),  &
                         wC(ix-1:ix+1,iy-1:iy+1,j-1:j+1,1),  &
                         Rho(ix-1:ix+1,iy-1:iy+1,j-1:j+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,j-1:j+1),  &
                         FU(ix-1:ix,iy:iy,j:j),              &
                         FV(ix:ix,iy-1:iy,j:j),              &
                         FW(ix:ix,iy:iy,j-1:j),              &
                         S,O,dudz)
        VCe   = VolC(ix  ,iy  ,j  )
        VB    = VolC(ix  ,iy  ,j-1)
        VT    = VolC(ix  ,iy  ,j+1)
        FB    = FW(ix  ,iy  ,j-1)
        FT    = FW(ix  ,iy  ,j  )
        rhoC  = Rho(ix,iy,j,1)
        thC   = th(ix,iy,j,1)/(rhoC+Eps)
        tkeC  = MAX(tke(ix,iy,j,1),Zero)
        dthdz = GradCentr(th(ix,iy,j-1,1)/(Rho(ix,iy,j-1,1)+Eps) &
                          ,thC &
                          ,th(ix,iy,j+1,1)/(Rho(ix,iy,j+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)
        N2 = Grav/(thC+Eps)*dthdz/(S*O+1.d-3)*S*O
        B  = -D(ix,iy,j,1)/Prt*N2

        TempU=(uCR(ix,iy,j,1)*FU(ix,iy,j)+uCL(ix,iy,j,1)*FU(ix-1,iy,j))/(FU(ix-1,iy,j)+FU(ix,iy,j)+Eps)
        TempV=(vCR(ix,iy,j,1)*FV(ix,iy,j)+vCL(ix,iy,j,1)*FV(ix,iy-1,j))/(FV(ix,iy-1,j)+FV(ix,iy,j)+Eps)
        TempW=(wCR(ix,iy,j,1)*FW(ix,iy,j)+wCL(ix,iy,j,1)*FW(ix,iy,j-1))/(FW(ix,iy,j-1)+FW(ix,iy,j)+Eps)
        TempVM=SQRT(TempU*TempU+TempV*TempV+TempW*TempW)
        Sd = 12.d0*SQRT(Cmy0)*Cd*0.5d0*LAD(j-iz+1)*TempVM*tkeC/(rhoC+Eps)

        tkeRhs(ix,iy,j,1) = tkeRhs(ix,iy,j,1) + B
        disRhs(ix,iy,j,1) = disRhs(ix,iy,j,1) + Cmy0*tkeC/(D(ix,iy,j,1)+Eps) * (Cmy1-Cmy2)*(B-Sd) ! dis/tke = Cmy0*tke/D
      END DO ! j
    END DO ! i
  END IF ! Canopy
END SUBROUTINE TkeDisCompute

SUBROUTINE JacTkeDisCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: tkeC,disC
  REAL(RealKind) :: S,O,dudz,uStar,z0,z,dzLoc

  INTEGER :: i,j,nrcl
  REAL(RealKind) :: dthdz,thC,rhoC
  REAL(RealKind) :: TempU,TempV,TempW ! velocities
  REAL(RealKind) :: TempVM ! mean velocity
  REAL(RealKind) :: FB,FT,VB,VCe,VT
  REAL(RealKind) :: N2,B ! Brunt-Vasala-Frequency and buoyance production of TKE induced by canopy
  REAL(RealKind) :: FUL,FUR,FVL,FVR,FWL,FWR
  REAL(RealKind) :: KSd,Kc,Kv
  REAL(RealKind), POINTER :: LAD(:)
  REAL(RealKind), PARAMETER :: Prt=0.85d0
  REAL(RealKind), PARAMETER :: Cd=0.15d0 ! drag coefficient of canopy

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),        &
                         FV(ix:ix,iy-1:iy,iz:iz),        &
                         FW(ix:ix,iy:iy,iz-1:iz),        &
                         S,O)

        tkeC = MAX(tke(ix,iy,iz,1),Zero)
        disC = MAX(dis(ix,iy,iz,1),Zero)
        tkeC = tke(ix,iy,iz,1)
        disC = dis(ix,iy,iz,1)

        AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) &
                                                          -Cmy0*Two*tkeC/(D(ix,iy,iz,1)+Eps)
        AS(IndexMet(disPosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(disPosJac,tkePosJac))%c(ix,iy,iz,1) &
                                                          +Cmy1*Cmy0*S*S  
        AS(IndexMet(disPosJac,disPosJac))%c(ix,iy,iz,1) = AS(IndexMet(disPosJac,disPosJac))%c(ix,iy,iz,1) &
                                                          -Cmy0*Cmy2 &
                                                          *SQRT(D(ix,iy,iz,1)/Cmy0)*&
                                                          &1.5d0*SQRT(ABS(disC))*SIGN(disC,1.0d0)&
                                                          &/(D(ix,iy,iz,1)+Eps)


      END DO
    END DO
  END DO

  IF (Canopy) THEN
    DO i=1,NumBoundCell
      ix = BoundCell(i)%ix
      iy = BoundCell(i)%iy
      iz = BoundCell(i)%iz
      nrcl = BoundCell(i)%CanopyCell%NrCanopyLayers
      LAD => BoundCell(i)%CanopyCell%LAD

      DO j=iz, iz+nrcl-1
        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,j-1:j+1,1), &
                         vC(ix-1:ix+1,iy-1:iy+1,j-1:j+1,1), &
                         wC(ix-1:ix+1,iy-1:iy+1,j-1:j+1,1), &
                         Rho(ix-1:ix+1,iy-1:iy+1,j-1:j+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,j-1:j+1),  &
                         FU(ix-1:ix,iy:iy,j:j),        &
                         FV(ix:ix,iy-1:iy,j:j),        &
                         FW(ix:ix,iy:iy,j-1:j),        &
                         S,O)
        FUL = FU(ix-1,iy,j)
        FUR = FU(ix,iy,j)
        FVL = FV(ix,iy-1,j)
        FVR = FV(ix,iy,j)
        FWL = FW(ix,iy,j-1)
        FWR = FW(ix,iy,j)

        VCe   = VolC(ix  ,iy  ,j  )
        VB    = VolC(ix  ,iy  ,j-1)
        VT    = VolC(ix  ,iy  ,j+1)
        FB    = FW(ix  ,iy  ,j-1)
        FT    = FW(ix  ,iy  ,j  )
        rhoC  = Rho(ix,iy,j,1)
        thC   = th(ix,iy,j,1)/(rhoC+Eps)
        tkeC = MAX(tke(ix,iy,j,1),Zero)
        dthdz = GradCentr(th(ix,iy,j-1,1)/(Rho(ix,iy,j-1,1)+Eps) &
                          ,thC &
                          ,th(ix,iy,j+1,1)/(Rho(ix,iy,j+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)
        N2 = Grav/(thC+Eps)*dthdz/(S*O+1.d-3)*S*O
        B  = -D(ix,iy,j,1)/Prt*N2

        TempU=(uCR(ix,iy,j,1)*FUR+uCL(ix,iy,j,1)*FUL)/(FUL+FUR+Eps)
        TempV=(vCR(ix,iy,j,1)*FVR+vCL(ix,iy,j,1)*FVL)/(FVL+FVR+Eps)
        TempW=(wCR(ix,iy,j,1)*FWR+wCL(ix,iy,j,1)*FWL)/(FWL+FWR+Eps)
        TempVM=SQRT(TempU*TempU+TempV*TempV+TempW*TempW)
        KSd= 12.d0*SQRT(Cmy0)*Cd*0.5d0*LAD(j-iz+1)/(rhoC+Eps)
        Kc = Cmy0*(Cmy1-Cmy2)/(D(ix,iy,j,1)+Eps)
        Kv = -Kc*KSd*tkeC*tkeC/(TempVM+Eps)

        AS(IndexMet(disPosJac,tkePosJac))%c(ix,iy,j,1) = AS(IndexMet(disPosJac,tkePosJac))%c(ix,iy,j,1) &
                                                          +Kc*B -2*Kc*KSd*tkeC*TempVM
        AS(IndexMet(disPosJac,uPosLJac))%c(ix,iy,j,1) = AS(IndexMet(disPosJac,uPosLJac))%c(ix,iy,j,1) &
                                                         +Kv*TempU*FUL/(FUL+FUR+Eps)
        AS(IndexMet(disPosJac,uPosRJac))%c(ix,iy,j,1) = AS(IndexMet(disPosJac,uPosRJac))%c(ix,iy,j,1) &
                                                         +Kv*TempU*FUR/(FUL+FUR+Eps)
        AS(IndexMet(disPosJac,vPosLJac))%c(ix,iy,j,1) = AS(IndexMet(disPosJac,vPosLJac))%c(ix,iy,j,1) &
                                                         +Kv*TempV*FVL/(FVL+FVR+Eps)
        AS(IndexMet(disPosJac,vPosRJac))%c(ix,iy,j,1) = AS(IndexMet(disPosJac,vPosRJac))%c(ix,iy,j,1) &
                                                         +Kv*TempV*FVR/(FVL+FVR+Eps)
        AS(IndexMet(disPosJac,wPosLJac))%c(ix,iy,j,1) = AS(IndexMet(disPosJac,wPosLJac))%c(ix,iy,j,1) &
                                                         +Kv*TempW*FWL/(FWL+FWR+Eps)
        AS(IndexMet(disPosJac,wPosRJac))%c(ix,iy,j,1) = AS(IndexMet(disPosJac,wPosRJac))%c(ix,iy,j,1) &
                                                         +Kv*TempW*FWR/(FWL+FWR+Eps)
      END DO ! j
    END DO ! i
  END IF ! Canopy
END SUBROUTINE JacTkeDisCompute

SUBROUTINE TkeOmeCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: tkeC,omeC,rhoC
  REAL(RealKind) :: S,O,dudz,uStar,z0,z,dzLoc
! Quellen und Senken der Tke-Dis-Gleichungen
! Zentrale Gradienten linear interpoliert zwischen den prim�en Gradienten

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        CALL DeformVortN(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),              &
                         FV(ix:ix,iy-1:iy,iz:iz),              &
                         FW(ix:ix,iy:iy,iz-1:iz),              &
                         S,O,dudz)

        rhoC = Rho(ix,iy,iz,1)
        tkeC = MAX(tke(ix,iy,iz,1),Zero)
        omeC = MAX(ome(ix,iy,iz,1),Zero)
        tkeRhs(ix,iy,iz,1) = tkeRhs(ix,iy,iz,1)  &
                             +D(ix,iy,iz,1)*S*O-omeC*tkeC/(rhoC+Eps)
        omeRhs(ix,iy,iz,1) = omeRhs(ix,iy,iz,1)  &
                            + Come1*omeC/(tkeC+Eps)*D(ix,iy,iz,1)*S*O &
                            -Come2*omeC*omeC/(rhoC+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE TkeOmeCompute

SUBROUTINE JacTkeOmeCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: tkeC,omeC,rhoC
  REAL(RealKind) :: S,O,dudz,uStar,z0,z,dzLoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        
        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),        &
                         FV(ix:ix,iy-1:iy,iz:iz),        &
                         FW(ix:ix,iy:iy,iz-1:iz),        &
                         S,O)

        tkeC = MAX(tke(ix,iy,iz,1),Zero)
        omeC = MAX(ome(ix,iy,iz,1),Zero)
        rhoC = Rho(ix,iy,iz,1)

        AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) &
                                                          -omeC/(rhoC+Eps)
        AS(IndexMet(tkePosJac,omePosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,omePosJac))%c(ix,iy,iz,1) &
                                                          -tkeC/(rhoC+Eps)
        AS(IndexMet(omePosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(omePosJac,tkePosJac))%c(ix,iy,iz,1) &
                                                          -omeC*Come1*D(ix,iy,iz,1)*S*O/(tkeC**2+Eps)
        AS(IndexMet(omePosJac,omePosJac))%c(ix,iy,iz,1) = AS(IndexMet(omePosJac,omePosJac))%c(ix,iy,iz,1) &
                                                          +Come1*D(ix,iy,iz,1)*S*O/(tkeC+Eps) -2*Come2*omeC/(rhoC+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE JacTkeOmeCompute

SUBROUTINE TkeSGSCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: S,O
  REAL(RealKind) :: N2
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC,rhoC
  REAL(RealKind) :: tkeC,Ptke,Rich,DMomLoc,DPotLoc,LenLoc
  REAL(RealKind) :: CDis

  INTEGER :: i,j,nrcl
  REAL(RealKind) :: TempU,TempV,TempW ! velocities
  REAL(RealKind) :: TempVM ! mean velocity
  REAL(RealKind), POINTER :: LAD(:)
  REAL(RealKind), PARAMETER :: Cd=0.15d0 ! drag coefficient of canopy

! inserted by Christian M
  REAL(RealKind), PARAMETER :: Cd_Tree=0.2d0 ! drag coefficient of canopy
  TYPE(TreePoint_T), POINTER :: TreePoint(:)
  INTEGER :: NumPoints
  REAL(RealKind) :: c

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe     = VolC(ix  ,iy  ,iz  )
        VB      = VolC(ix  ,iy  ,iz-1)
        VT      = VolC(ix  ,iy  ,iz+1  )
        FB      = FW(ix  ,iy  ,iz-1)
        FT      = FW(ix  ,iy  ,iz  )
        tkeC    = MAX(tke(ix,iy,iz,1),Zero)
        rhoC    = Rho(ix,iy,iz,1)
        thC     = th(ix,iy,iz,1)/(rhoC+Eps)
        DMomLoc = DMom(ix,iy,iz,1)
        DPotLoc = DPot(ix,iy,iz,1)
        LenLoc  = LenScaleGeom

        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),        &
                         FV(ix:ix,iy-1:iy,iz:iz),        &
                         FW(ix:ix,iy:iy,iz-1:iz),        &
                         S,O)

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                          ,thC &
                          ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)
 
!       Brunt-Vasala-Frequency
        N2 = Grav/(thC+Eps)*dthdz
        N2 = Grav/(thC+Eps)*dthdz/(S*O+1.d-3)*S*O

!       Richardson-Zahl und Faktoren fuer Produktionsterme:
        Rich = N2/(S*O+1.d-3)
        Rich = MAX(Rich,-Four/Three)
        Ptke = MAX(1.d0-DPotLoc*Rich/( SigT*DMomLoc+Eps ),ftke)              !Faktor fuer tke-Produktion
!       Local Length Scale
        tkeC = tke(ix,iy,iz,1)
        IF (N2>Zero) THEN
          LenLoc=MIN(LenScaleGeom,0.76*SQRT(ABS(tkeC)/(RhoC+Eps)/N2))
        ELSE 
          LenLoc=LenScaleGeom
        END IF
        
        CDis=1.9d0*Cs+(0.93d0-1.9d0*Cs)*LenLoc/LenScaleGeom
        tkeRhs(ix,iy,iz,1) = tkeRhs(ix,iy,iz,1)  &
                             +DMomLoc*S*O-DPotLoc*N2 &
                             -CDis*tkeC*SQRT(ABS(tkeC))/(LenLoc*SQRT(rhoC)+Eps)


      END DO
    END DO
  END DO

  IF (Canopy) THEN
    DO i=1, NumBoundCell
      ix = BoundCell(i)%ix
      iy = BoundCell(i)%iy
      iz = BoundCell(i)%iz
      nrcl = BoundCell(i)%CanopyCell%NrCanopyLayers
      LAD => BoundCell(i)%CanopyCell%LAD

      DO j=iz, iz+nrcl-1
        tkeC    = MAX(tke(ix,iy,j,1),Zero)
        rhoC    = Rho(ix,iy,j,1)

        TempU=(uCR(ix,iy,j,1)*FU(ix,iy,j)+uCL(ix,iy,j,1)*FU(ix-1,iy,j))/(FU(ix-1,iy,j)+FU(ix,iy,j)+Eps)
        TempV=(vCR(ix,iy,j,1)*FV(ix,iy,j)+vCL(ix,iy,j,1)*FV(ix,iy-1,j))/(FV(ix,iy-1,j)+FV(ix,iy,j)+Eps)
        TempW=(wCR(ix,iy,j,1)*FW(ix,iy,j)+wCL(ix,iy,j,1)*FW(ix,iy,j-1))/(FW(ix,iy,j-1)+FW(ix,iy,j)+Eps)
        TempVM=SQRT(TempU*TempU+TempV*TempV+TempW*TempW)

        tkeRhs(ix,iy,j,1) = tkeRhs(ix,iy,j,1) - Cd*LAD(j-iz+1)*TempVM*tkeC/(rhoC+Eps)
      END DO ! j
    END DO ! i
  END IF ! Canopy

  IF (Baum) THEN ! Christian M.
  NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
  TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint
  DO i=1, NumPoints
    ix = TreePoint(i)%ix
    iy = TreePoint(i)%iy
    iz = TreePoint(i)%iz
    c = 2.0d0*TreePoint(i)%c

    tkeC    = MAX(tke(ix,iy,iz,1),Zero)
    rhoC    = Rho(ix,iy,iz,1)

    TempU=(uCR(ix,iy,iz,1)*FU(ix,iy,iz)+uCL(ix,iy,iz,1)*FU(ix-1,iy,iz))/(FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
    TempV=(vCR(ix,iy,iz,1)*FV(ix,iy,iz)+vCL(ix,iy,iz,1)*FV(ix,iy-1,iz))/(FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
    TempW=(wCR(ix,iy,iz,1)*FW(ix,iy,iz)+wCL(ix,iy,iz,1)*FW(ix,iy,iz-1))/(FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)
    TempVM=SQRT(TempU*TempU+TempV*TempV+TempW*TempW)

    tkeRhs(ix,iy,iz,1) = tkeRhs(ix,iy,iz,1) - c*TempVM*tkeC/(rhoC+Eps)
  END DO ! i
  END IF
END SUBROUTINE TkeSGSCompute

SUBROUTINE JacTkeSGSCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: S,O
  REAL(RealKind) :: N2
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC,rhoC
  REAL(RealKind) :: tkeC,Ptke,Rich,DMomLoc,DPotLoc,LenLoc
  REAL(RealKind) :: CDis

  INTEGER :: i,j,nrcl
  REAL(RealKind), POINTER :: LAD(:)
  REAL(RealKind), PARAMETER :: Cd=0.15d0 ! drag coefficient of canopy
  REAL(RealKind) :: TempU,TempV,TempW,TempVM
  REAL(RealKind) :: K1,K2
  REAL(RealKind) :: FUL,FUR,FVL,FVR,FWL,FWR

  ! inserted by Christian Markwitz
  REAL(RealKind), PARAMETER :: Cd_Tree=0.2d0 ! drag coefficient of canopy
  TYPE(TreePoint_T), POINTER :: TreePoint(:)
  INTEGER :: NumPoints
  REAL(RealKind) :: c

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe     = VolC(ix  ,iy  ,iz  )
        VB      = VolC(ix  ,iy  ,iz-1)
        VT      = VolC(ix  ,iy  ,iz+1  )
        FB      = FW(ix  ,iy  ,iz-1)
        FT      = FW(ix  ,iy  ,iz  )
        tkeC    = tke(ix,iy,iz,1)
        rhoC    = Rho(ix,iy,iz,1)
        thC     = th(ix,iy,iz,1)/(rhoC+Eps)
        DMomLoc = DMom(ix,iy,iz,1)
        DPotLoc = DPot(ix,iy,iz,1)
        LenLoc  = LenScaleGeom
         
        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),        &
                         FV(ix:ix,iy-1:iy,iz:iz),        &
                         FW(ix:ix,iy:iy,iz-1:iz),        &
                         S,O)

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                          ,thC &
                          ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)

!       Brunt-Vasala-Frequency
        N2 = Grav/(thC+Eps)*dthdz
        N2 = Grav/(thC+Eps)*dthdz/(S*O+1.d-3)*S*O
!       Local Length Scale
        IF (N2>Zero) THEN
          LenLoc=MIN(LenScaleGeom,0.76*SQRT(ABS(tkeC)/(RhoC+Eps)/N2))
        ELSE 
          LenLoc=LenScaleGeom
        END IF
        CDis=1.9d0*Cs+(0.93d0-1.9d0*Cs)*LenLoc/LenScaleGeom
        AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) &
                                                         -CDis*(Three/Two)*SQRT(ABS(tkeC)) &
                                                         /(LenLoc*SQRT(rhoC)+Eps) 


      END DO
    END DO
  END DO

  IF (Canopy) THEN
    DO i=1, NumBoundCell
      ix = BoundCell(i)%ix
      iy = BoundCell(i)%iy
      iz = BoundCell(i)%iz
      nrcl = BoundCell(i)%CanopyCell%NrCanopyLayers
      LAD => BoundCell(i)%CanopyCell%LAD

      DO j=iz, iz+nrcl-1
        tkeC    = MAX(tke(ix,iy,j,1),Zero)
        rhoC    = Rho(ix,iy,j,1)

        FUL = FU(ix-1,iy,j)
        FUR = FU(ix,iy,j)
        FVL = FV(ix,iy-1,j)
        FVR = FV(ix,iy,j)
        FWL = FW(ix,iy,j-1)
        FWR = FW(ix,iy,j)

        TempU=(uCR(ix,iy,j,1)*FUR+uCL(ix,iy,j,1)*FUL)/(FUL+FUR+Eps)
        TempV=(vCR(ix,iy,j,1)*FVR+vCL(ix,iy,j,1)*FVL)/(FVL+FVR+Eps)
        TempW=(wCR(ix,iy,j,1)*FWR+wCL(ix,iy,j,1)*FWL)/(FWL+FWR+Eps)
        TempVM=SQRT(TempU*TempU+TempV*TempV+TempW*TempW)
        K1 = Cd*LAD(j-iz+1)
        K2 = K1*tkeC/(TempVM*rhoC+Eps)

        AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,j,1) = AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,j,1) &
                                                          -K1*TempVM/(rhoC+Eps)
        AS(IndexMet(tkePosJac,uPosLJac))%c(ix,iy,j,1) = AS(IndexMet(tkePosJac,uPosLJac))%c(ix,iy,j,1) &
                                                         -K2*TempU*FUL/(FUL+FUR+Eps)
        AS(IndexMet(tkePosJac,uPosRJac))%c(ix,iy,j,1) = AS(IndexMet(tkePosJac,uPosRJac))%c(ix,iy,j,1) &
                                                         -K2*TempU*FUR/(FUL+FUR+Eps)
        AS(IndexMet(tkePosJac,vPosLJac))%c(ix,iy,j,1) = AS(IndexMet(tkePosJac,vPosLJac))%c(ix,iy,j,1) &
                                                         -K2*TempV*FVL/(FVL+FVR+Eps)
        AS(IndexMet(tkePosJac,vPosRJac))%c(ix,iy,j,1) = AS(IndexMet(tkePosJac,vPosRJac))%c(ix,iy,j,1) &
                                                         -K2*TempV*FVR/(FVL+FVR+Eps)
        AS(IndexMet(tkePosJac,wPosLJac))%c(ix,iy,j,1) = AS(IndexMet(tkePosJac,wPosLJac))%c(ix,iy,j,1) &
                                                         -K2*TempW*FWL/(FWL+FWR+Eps)
        AS(IndexMet(tkePosJac,wPosRJac))%c(ix,iy,j,1) = AS(IndexMet(tkePosJac,wPosRJac))%c(ix,iy,j,1) &
                                                         -K2*TempW*FWR/(FWL+FWR+Eps)
      END DO ! j
    END DO ! i
  END IF ! Canopy

  IF (Baum) THEN
  NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
  TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint
    DO i=1, NumPoints
      ix = TreePoint(i)%ix
      iy = TreePoint(i)%iy
      iz = TreePoint(i)%iz
      c = 2.0d0*TreePoint(i)%c

      tkeC    = MAX(tke(ix,iy,iz,1),Zero)
      rhoC    = Rho(ix,iy,iz,1)

      FUL = FU(ix-1,iy,iz)
      FUR = FU(ix,iy,iz)
      FVL = FV(ix,iy-1,iz)
      FVR = FV(ix,iy,iz)
      FWL = FW(ix,iy,iz-1)
      FWR = FW(ix,iy,iz)

      TempU=(uCR(ix,iy,iz,1)*FUR+uCL(ix,iy,iz,1)*FUL)/(FUL+FUR+Eps)
      TempV=(vCR(ix,iy,iz,1)*FVR+vCL(ix,iy,iz,1)*FVL)/(FVL+FVR+Eps)
      TempW=(wCR(ix,iy,iz,1)*FWR+wCL(ix,iy,iz,1)*FWL)/(FWL+FWR+Eps)
      TempVM=SQRT(TempU*TempU+TempV*TempV+TempW*TempW)
      K1 = c
      K2 = K1*tkeC/(TempVM*rhoC+Eps)

      AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) &
                                                          -K1*TempVM/(rhoC+Eps)
      AS(IndexMet(tkePosJac,uPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,uPosLJac))%c(ix,iy,iz,1) &
                                                         -K2*TempU*FUL/(FUL+FUR+Eps)
      AS(IndexMet(tkePosJac,uPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,uPosRJac))%c(ix,iy,iz,1) &
                                                         -K2*TempU*FUR/(FUL+FUR+Eps)
      AS(IndexMet(tkePosJac,vPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,vPosLJac))%c(ix,iy,iz,1) &
                                                         -K2*TempV*FVL/(FVL+FVR+Eps)
      AS(IndexMet(tkePosJac,vPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,vPosRJac))%c(ix,iy,iz,1) &
                                                         -K2*TempV*FVR/(FVL+FVR+Eps)
      AS(IndexMet(tkePosJac,wPosLJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,wPosLJac))%c(ix,iy,iz,1) &
                                                         -K2*TempW*FWL/(FWL+FWR+Eps)
      AS(IndexMet(tkePosJac,wPosRJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,wPosRJac))%c(ix,iy,iz,1) &
                                                         -K2*TempW*FWR/(FWL+FWR+Eps)
    END DO 
  END IF
END SUBROUTINE JacTkeSGSCompute

SUBROUTINE TkeLenCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: S,O
  REAL(RealKind) :: N2,Len0,Len0s
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC,rhoC
  REAL(RealKind) :: tkeC,Ptke,Rich,DMomLoc,DPotLoc
  REAL(RealKind) :: dzLoc,ds,zPL

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe     = VolC(ix  ,iy  ,iz  )
        VB      = VolC(ix  ,iy  ,iz-1)
        VT      = VolC(ix  ,iy  ,iz+1  )
        FB      = FW(ix  ,iy  ,iz-1)
        FT      = FW(ix  ,iy  ,iz  )
        tkeC    = MAX(tke(ix,iy,iz,1),Zero)
        rhoC    = Rho(ix,iy,iz,1)
        thC     = th(ix,iy,iz,1)/(rhoC+Eps)
        DMomLoc = DMom(ix,iy,iz,1)
        DPotLoc = DPot(ix,iy,iz,1)
        dzLoc   = MAXVAL(dz)

        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),              &
                         FV(ix:ix,iy-1:iy,iz:iz),              &
                         FW(ix:ix,iy:iy,iz-1:iz),              &
                         S,O)

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps)  &
                          ,thC                                      &
                          ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)
 
!	Brunt-Vaisala-Frequency N**2
        N2 = Grav/(thC+Eps)*dthdz

!	Calculate height z
        IF (FW(ix,iy,iz-1)>FW(ix,iy,iz)) THEN
          dzLoc = MIN(VolC(ix,iy,iz)/(FW(ix,iy,iz-1)+Eps),dz(iz))
          zPL   = zP(iz-1)+0.5e0*dzLoc
        ELSE
          dzLoc = MIN(VolC(ix,iy,iz)/(FW(ix,iy,iz)+Eps),dz(iz))
          zPL   = zP(iz)-0.5e0*dzLoc
        END IF
!       Mean grid resolution ds
        ds = (two*VolC(ix,iy,iz)*(One/(FU(ix+1,iy,iz)+FU(ix,iy,iz)+Eps) &
                  +One/(FV(ix,iy+1,iz)+FV(ix,iy,iz)+Eps)+One/(FW(ix,iy,iz+1)+FW(ix,iy,iz)+Eps)))/Three

!       Turbulent Length Scale
        Len0s = 0.54d0*SQRT(tkeC/(ABS(N2)+Eps))

!       Equilibrium turbulent Length Scale Len0
        Len0 = Rho(ix,iy,iz,1)    &
                  *MIN(0.67d0*zPL &
                  ,ds             &
                  ,Len0s)

        LenRhs(ix,iy,iz,1) = LenRhs(ix,iy,iz,1)-c2*SQRT(rhoC*tkeC)           &
                             *(LenC(ix,iy,iz,1)-Len0)/(LenC(ix,iy,iz,1)+Eps)

!       Richardson-Zahl und Faktoren fuer Produktionsterme:
        Rich = N2/(S*O+1.d-3)
        Rich = MAX(Rich,-Four/Three)
        Ptke = MAX(1.d0-DPotLoc*Rich/( SigT*DMomLoc+Eps ),ftke)              !Faktor fuer tke-Produktion
        
        tkeC = tke(ix,iy,iz,1)
        tkeRhs(ix,iy,iz,1) = tkeRhs(ix,iy,iz,1)      &
                             +DMomLoc*S*O-DPotLoc*N2 &
                             -Cmy0**(Three/Four)*tkeC*SQRT(ABS(tkeC))/(LenC(ix,iy,iz,1)*SQRT(rhoC)+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE TkeLenCompute

SUBROUTINE JacTkeLenCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: N2,Len0,Len0s
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC,rhoC
  REAL(RealKind) :: tkeC
  REAL(RealKind) :: dzLoc,zPL,ds

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe     = VolC(ix  ,iy  ,iz  )
        VB      = VolC(ix  ,iy  ,iz-1)
        VT      = VolC(ix  ,iy  ,iz+1  )
        FB      = FW(ix  ,iy  ,iz-1)
        FT      = FW(ix  ,iy  ,iz  )
        tkeC    = MAX(tke(ix,iy,iz,1),Zero)
        rhoC    = Rho(ix,iy,iz,1)
        thC     = th(ix,iy,iz,1)/(rhoC+Eps)
        dzLoc   = MAXVAL(dz)

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps)  &
                          ,thC                                      &
                          ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)
 
!	Brunt-Vaisala-Frequency N**2
        N2 = Grav/(thC+Eps)*dthdz

!	Calculate height z
        IF (FW(ix,iy,iz-1)>FW(ix,iy,iz)) THEN
          dzLoc = MIN(VolC(ix,iy,iz)/(FW(ix,iy,iz-1)+Eps),dz(iz))
          zPL   = zP(iz-1)+0.5e0*dzLoc
        ELSE
          dzLoc = MIN(VolC(ix,iy,iz)/(FW(ix,iy,iz)+Eps),dz(iz))
          zPL   = zP(iz)-0.5e0*dzLoc
        END IF
!       Mean grid resolution ds
        ds = (two*VolC(ix,iy,iz)*(One/(FU(ix+1,iy,iz)+FU(ix,iy,iz)+Eps) &
                  +One/(FV(ix,iy+1,iz)+FV(ix,iy,iz)+Eps)+One/(FW(ix,iy,iz+1)+FW(ix,iy,iz)+Eps)))/Three

!       Turbulent Length Scale
        Len0s = 0.54d0*SQRT(tkeC/(ABS(N2)+Eps))

!       Equilibrium turbulent Length Scale Len0
        Len0 = Rho(ix,iy,iz,1)    &
                  *MIN(0.67d0*zPL &
                  ,ds             &
                  ,Len0s)

!       ---Derivatives tkeRhs tke,Len---
!                                   tkeHRhs(ix,iy,iz,1) = tkeHRhs(ix,iy,iz,1) &
!                                                         +DMomLoc*S*O-DPotLoc*N2 &
!                                                         -Cmy0**(Three/Four)*tkeC*SQRT(ABS(tkeC))
!                                                         /(LenC*SQRT(rhoC)+Eps)
        AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1)               &
!                                                         -Cmy4*(Three/Two)*SIGN(One,tkeC)*SQRT(ABS(tkeC)) &
                                                         -Cmy0**(Three/Four)*(Three/Two)*SIGN(One,tkeC)*SQRT(ABS(tkeC)) &
                                                         /(LenC(ix,iy,iz,1)*SQRT(rhoC)+Eps)                  !!!FILAUS
        AS(IndexMet(tkePosJac,LenPosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,LenPosJac))%c(ix,iy,iz,1)  &
!                                                         -Cmy4*(Three/Two)*SIGN(One,tkeC)*SQRT(ABS(tkeC)) &
                                                         +Cmy0**(Three/Four)*tkeC*SQRT(ABS(tkeC))          &
                                                         /(LenC(ix,iy,iz,1)**Two*SQRT(rhoC)+Eps)                  !!!FILAUS
!       ---Derivatives LenRhs tke,Len---
!                                    LenRhs(ix,iy,iz,1) = LenRhs(ix,iy,iz,1)-c2*SQRT(rhoC*tkeC) &
!                                                         *(LenC(ix,iy,iz,1)-Len0)/(LenC(ix,iy,iz,1)+Eps)
        AS(IndexMet(LenPosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(LenPosJac,tkePosJac))%c(ix,iy,iz,1) &
                                                          -c2*SQRT(rhoC)*(LenC(ix,iy,iz,1)-Len0)          &
                                                          !/(SQRT(tkeC)*LenC(ix,iy,iz,1)+Eps) ! = falsch nach tke abgeleitet
                                                          /(SQRT(tkeC)*LenC(ix,iy,iz,1)*Two+Eps) ! Hinneburg
        AS(IndexMet(LenPosJac,LenPosJac))%c(ix,iy,iz,1) = AS(IndexMet(LenPosJac,LenPosJac))%c(ix,iy,iz,1) &
                                                          !-c2*SQRT(rhoC)*Len0/(LenC(ix,iy,iz,1)**Two+Eps) ! = falsch nach Len abgeleitet
                                                          -c2*SQRT(rhoC)*SQRT(tkeC)*Len0/(LenC(ix,iy,iz,1)**Two+Eps) ! Hinneburg
      END DO
    END DO
  END DO
END SUBROUTINE JacTkeLenCompute

SUBROUTINE TkeDisRichCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FB,FT
  REAL(RealKind) :: VB,VT,VCe
  REAL(RealKind) :: tkeC,disC,rhoC
  REAL(RealKind) :: S,O
  REAL(RealKind) :: Rich,dthdz,thC,Ptke,Pdis
  REAL(RealKind) :: uStar
  REAL(RealKind) :: GammaPos,TkeCPos,DisCPos

! Quellen und Senken der Tke-Dis-Gleichungen
! Zentrale Gradienten linear interpoliert zwischen den prim�en Gradienten

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       Zentrale Groessen der betrachteten Box
        VCe  = VolC(ix  ,iy  ,iz  )
        VB   = VolC(ix  ,iy  ,iz-1)
        VT   = VolC(ix  ,iy  ,iz+1  )
        FB   = FW(ix  ,iy  ,iz-1)
        FT   = FW(ix  ,iy  ,iz  )
        rhoC = Rho(ix,iy,iz,1)
        thC  = th(ix,iy,iz,1)/(rhoC+Eps)

        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),        &
                         FV(ix:ix,iy-1:iy,iz:iz),        &
                         FW(ix:ix,iy:iy,iz-1:iz),        &
                         S,O)

!       Richardson-Zahl und Faktoren fr Produktionsterme:
        dthdz =GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                       ,thC &
                       ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                       ,FB,FT,VB,VCe,VT)
        Rich  = (Grav/(thC+Eps)) * dthdz/(S*S+1.d-3)

        Rich     = MAX(Rich,-4.0e0_RealKind/3.0e0_RealKind)
        Ptke     = MAX(1.d0-Rich/SigT,ftke)              !Faktor fr tke-Produktion
        Pdis     = MAX(1.d0-Rich*Cmy3/(Cmy1*SigT),ftke)  !Faktor fr dis-Produktion
        tkeC     = tke(ix,iy,iz,1)
        tkeCPos  = MAX(tkeC,1.d-3)
        disC     = dis(ix,iy,iz,1)
        disCPos  = MAX(disC,1.d-6)
        GammaPos = MAX(DisCPos/tkeCPos,Zero)

        tkeRhs(ix,iy,iz,1) = tkeRhs(ix,iy,iz,1)      &
                             +D(ix,iy,iz,1)*S*O*Ptke &
                             -GammaPos*tkeC
        disRhs(ix,iy,iz,1) = disRhs(ix,iy,iz,1)  &
                             +GammaPos*(Cmy1*D(ix,iy,iz,1)*S*O*PTke-Cmy2*DisC)

      END DO
    END DO
  END DO

END SUBROUTINE TkeDisRichCompute

SUBROUTINE JacTkeDisRichCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: tkeC,disC,rhoC
  REAL(RealKind) :: GammaPos,TkeCPos,DisCPos

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        tkeC     = tke(ix,iy,iz,1)
        tkeCPos  = MAX(tkeC,1.d-3)
        disC     = dis(ix,iy,iz,1)
        disCPos  = MAX(disC,1.d-6)
        GammaPos = MAX(DisCPos/tkeCPos,Zero)

        AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkePosJac,tkePosJac))%c(ix,iy,iz,1) &
                                                          -GammaPos
        AS(IndexMet(disPosJac,disPosJac))%c(ix,iy,iz,1) = AS(IndexMet(disPosJac,disPosJac))%c(ix,iy,iz,1) &
                                                          -Cmy2*GammaPos

      END DO
    END DO
  END DO
END SUBROUTINE JacTkeDisRichCompute

SUBROUTINE TkeDisWallCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,zRauh,z
  REAL(RealKind) :: V,VT,VN
  REAL(RealKind) :: tkeC,rhoC

  DO i=1,NumBoundCell
    ix    = BoundCell(i)%ix
    iy    = BoundCell(i)%iy
    iz    = BoundCell(i)%iz
    n1    = BoundCell(i)%n1
    n2    = BoundCell(i)%n2
    n3    = BoundCell(i)%n3
    FL    = BoundCell(i)%FL+Eps
    zRauh = BoundCell(i)%zRauh
    z     = BoundCell(i)%dL
    rhoC  = Rho(ix,iy,iz,1)
    tkeC  = MAX(tke(ix,iy,iz,1),Zero)
    tkeC  = tke(ix,iy,iz,1)
    disRhs(ix,iy,iz,1) = disRhs(ix,iy,iz,1)                                     &
                        +(Cmy0)**0.75d0*FL*D(ix,iy,iz,1)/PrandtlNumber(disPos)  &       
                        *(ABS(tkeC)/(Rho(ix,iy,iz,1)+Eps))**0.5d0                    &
                        *(tkeC/(Rho(ix,iy,iz,1)+Eps))                           &
                        /(Karm*(z+zRauh)**2)/(VolC(ix,iy,iz)+Eps) 
  END DO

END SUBROUTINE TkeDisWallCompute

SUBROUTINE JacTkeDisWallCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,zRauh,z,logz
  REAL(RealKind) :: tkeC,rhoC,Factor
  REAL(RealKind) ,Parameter :: logmin=0.1d0,logmax=1.d10,zero=0.d0

  DO i=1,NumBoundCell
    ix    = BoundCell(i)%ix
    iy    = BoundCell(i)%iy
    iz    = BoundCell(i)%iz
    FL    = BoundCell(i)%FL+Eps
    zRauh = BoundCell(i)%zRauh
    z     = BoundCell(i)%dL
    tkeC  = MAX(tke(ix,iy,iz,1),Zero)
    tkeC  = tke(ix,iy,iz,1)
    rhoC  = Rho(ix,iy,iz,1)
    AS(IndexMet(disPosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(disPosJac,tkePosJac))%c(ix,iy,iz,1)        &
                                                      +(Cmy0)**0.75d0*FL*D(ix,iy,iz,1)/PrandtlNumber(disPos) &
                                                      *1.5d0*SQRT(ABS(tkeC))*SIGN(tkeC,1.0d0)/(Rho(ix,iy,iz,1)+Eps)**1.5d0
  END DO

END SUBROUTINE JacTkeDisWallCompute

SUBROUTINE TkeDisBaumCompute 

  INTEGER :: NumPoints
  INTEGER :: ix,iy,iz
  TYPE(TreePoint_T), POINTER :: TreePoint(:)

  REAL(RealKind) :: FL
  REAL(RealKind) :: tkeC

  NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
  TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint
  DO i=1,NumPoints
    ix=TreePoint(i)%ix
    iy=TreePoint(i)%iy
    iz=TreePoint(i)%iz
    FL=TreePoint(i)%c*FW(ix,iy,iz)+Eps
    tkeC  = MAX(tke(ix,iy,iz,1),Zero)
    disRhs(ix,iy,iz,1) = disRhs(ix,iy,iz,1)                                     &
                        +(Cmy0)**0.75d0*FL*D(ix,iy,iz,1)/PrandtlNumber(disPos)  &       
                        *(tkeC/(Rho(ix,iy,iz,1)+Eps))**1.5d0                    &
                        /(Karm*(distBaum+zRauhBaum)**2)/(VolC(ix,iy,iz)+Eps)
  END DO

END SUBROUTINE TkeDisBaumCompute

SUBROUTINE JacTkeDisBaumCompute

  INTEGER :: NumPoints
  INTEGER :: ix,iy,iz
  TYPE(TreePoint_T), POINTER :: TreePoint(:)

  REAL(RealKind) :: FL
  REAL(RealKind) :: tkeC

  NumPoints=PointTree%TreePointBlock(ibLoc)%NumTreePoint
  TreePoint=>PointTree%TreePointBlock(ibLoc)%TreePoint
  DO i=1,NumPoints
    ix=TreePoint(i)%ix
    iy=TreePoint(i)%iy
    iz=TreePoint(i)%iz
    FL=TreePoint(i)%c*FW(ix,iy,iz)+Eps
    tkeC  = MAX(tke(ix,iy,iz,1),Zero)

    AS(IndexMet(disPosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(disPosJac,tkePosJac))%c(ix,iy,iz,1)        &
                                                      +(Cmy0)**0.75d0*FL*D(ix,iy,iz,1)/PrandtlNumber(disPos) &
                                                      *1.5d0*tkeC**0.5d0/(Rho(ix,iy,iz,1)+Eps)**1.5d0        &
                                                      /(Karm*(distBaum+zRauhBaum)**2)/(VolC(ix,iy,iz)+Eps)
  END DO

END SUBROUTINE JacTkeDisBaumCompute

SUBROUTINE TkeOmeWallCompute
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,zRauh,z
  REAL(RealKind) :: V,VT,VN
  REAL(RealKind) :: tkeC,omeC,rhoC

  DO i=1,NumBoundCell
    ix    = BoundCell(i)%ix
    iy    = BoundCell(i)%iy
    iz    = BoundCell(i)%iz
    n1    = BoundCell(i)%n1
    n2    = BoundCell(i)%n2
    n3    = BoundCell(i)%n3
    FL    = BoundCell(i)%FL+Eps
    zRauh = BoundCell(i)%zRauh
    z     = BoundCell(i)%dL
    rhoC  = Rho(ix,iy,iz,1)
    tkeC  = MAX(tke(ix,iy,iz,1),Zero)
    omeC  = MAX(ome(ix,iy,iz,1),Zero)
    omeRhs(ix,iy,iz,1) = omeRhs(ix,iy,iz,1)                                     &
                        +(Cmy0)**0.75d0*FL*D(ix,iy,iz,1)/PrandtlNumber(omePos)  &
                        *(tkeC/(Rho(ix,iy,iz,1)+Eps))**0.5d0                    &
                        /(Karm*(z+zRauh)**2)/(VolC(ix,iy,iz)+Eps)
  END DO
END SUBROUTINE TkeOmeWallCompute

SUBROUTINE JacTkeOmeWallCompute
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: FL,zRauh,z,logz
  REAL(RealKind) :: tkeC,rhoC,Factor
  REAL(RealKind) ,Parameter :: logmin=0.1d0,logmax=1.d10,zero=0.d0

  DO i=1,NumBoundCell
    ix    = BoundCell(i)%ix
    iy    = BoundCell(i)%iy
    iz    = BoundCell(i)%iz
    FL    = BoundCell(i)%FL+Eps
    zRauh = BoundCell(i)%zRauh
    z     = BoundCell(i)%dL
    tkeC  = MAX(tke(ix,iy,iz,1),Zero)
    rhoC  = Rho(ix,iy,iz,1)
    AS(IndexMet(omePosJac,tkePosJac))%c(ix,iy,iz,1) = AS(IndexMet(omePosJac,tkePosJac))%c(ix,iy,iz,1)        &
                                                      +(Cmy0)**0.75d0*FL*D(ix,iy,iz,1)/PrandtlNumber(omePos) &
                                                      *0.5d0/(tkeC*Rho(ix,iy,iz,1)+Eps)**0.5d0        &
                                                      /(Karm*(z+zRauh)**2)/(VolC(ix,iy,iz)+Eps)
  END DO
END SUBROUTINE JacTkeOmeWallCompute

SUBROUTINE TkeHVLenCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: FWe,FE,FS,FN,FB,FT
  REAL(RealKind) :: VW,VE,VS,VN,VB,VT,VCe
  REAL(RealKind) :: u,v,w
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dudx,dudy,dudz,dvdx,dvdy,dvdz,dwdx,dwdy,dwdz,dthdz,ds
  REAL(RealKind) :: thC,rhoC,tkeHC,tkeVC,tkeG
  REAL(RealKind) :: N2,Len0,Len0s
  REAL(RealKind) :: SH,SV,RH,RV,DissiH,DissiV
  REAL(RealKind) :: c3,zRauh
  REAL(RealKind) :: Bou,zPL,dzLoc
! Quellen und Senken der Tke-Dis-Gleichungen
! Zentrale Gradienten linear interpoliert zwischen den primaeren Gradienten

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       Zentrale Groessen der betrachteten Box
        VCe   = VolC(ix  ,iy  ,iz  )
        VW    = VolC(ix-1,iy  ,iz  )
        VE    = VolC(ix+1,iy  ,iz  )
        VS    = VolC(ix  ,iy-1,iz  )
        VN    = VolC(ix  ,iy+1,iz  )
        VB    = VolC(ix  ,iy  ,iz-1)
        VT    = VolC(ix  ,iy  ,iz+1  )
        FWe   = FU(ix-1,iy  ,iz  )
        FE    = FU(ix  ,iy  ,iz  )
        FS    = FV(ix  ,iy-1,iz  )
        FN    = FV(ix  ,iy  ,iz  )
        FB    = FW(ix  ,iy  ,iz-1)
        FT    = FW(ix  ,iy  ,iz  )
        rhoC  = Rho(ix,iy,iz,1)
        u     = Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(rhoC+Eps)
        v     = Half*(vCL(ix,iy,iz,1)+vCR(ix,iy,iz,1))/(rhoC+Eps)
        w     = Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(rhoC+Eps)
        thC   = th(ix,iy,iz,1)/(rhoC+Eps)
        tkeHC = tkeH(ix,iy,iz,1)
        tkeVC = tkeV(ix,iy,iz,1)
        tkeG  = tkeH(ix,iy,iz,1)+tkeV(ix,iy,iz,1)
        dzLoc = MaxVal(dz)
        zRauh = 0.1d0

        dudx = GradCentrT(Half*(uCL(ix-1,iy,iz,1)+uCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                         ,u &
                         ,Half*(uCL(ix+1,iy,iz,1)+uCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps) &
                         ,DH(ix-1,iy,iz,1),DH(ix,iy,iz,1),DH(ix+1,iy,iz,1)                   &
                         ,FWe,FE,VW,VCe,VE)
        dudy = GradCentrT(Half*(uCL(ix,iy-1,iz,1)+uCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                        ,u &
                        ,Half*(uCL(ix,iy+1,iz,1)+uCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps)  &
                        ,DH(ix,iy-1,iz,1),DH(ix,iy,iz,1),DH(ix,iy+1,iz,1)                    &
                        ,FS,FN,VS,VCe,VN)
        dudz = GradCentrT(Half*(uCL(ix,iy,iz-1,1)+uCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                        ,u &
                        ,Half*(uCL(ix,iy,iz+1,1)+uCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps)  &
                        ,DV(ix,iy,iz-1,1),DV(ix,iy,iz,1),DV(ix,iy,iz+1,1)                    &
                        ,FB,FT,VB,VCe,VT)
        dvdx = GradCentrT(Half*(vCL(ix-1,iy,iz,1)+vCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                        ,v &
                        ,Half*(vCL(ix+1,iy,iz,1)+vCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps)  &
                        ,DH(ix-1,iy,iz,1),DH(ix,iy,iz,1),DH(ix+1,iy,iz,1)                    &
                        ,FWe,FE,VW,VCe,VE)
        dvdy = GradCentrT(Half*(vCL(ix,iy-1,iz,1)+vCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                        ,v &
                        ,Half*(vCL(ix,iy+1,iz,1)+vCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps)  &
                        ,DH(ix,iy-1,iz,1),DH(ix,iy,iz,1),DH(ix,iy+1,iz,1)                    &
                        ,FS,FN,VS,VCe,VN)
        dvdz = GradCentrT(Half*(vCL(ix,iy,iz-1,1)+vCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                        ,v &
                        ,Half*(vCL(ix,iy,iz+1,1)+vCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps)  &
                        ,DV(ix,iy,iz-1,1),DV(ix,iy,iz,1),DV(ix,iy,iz+1,1)                    &
                        ,FB,FT,VB,VCe,VT)
        dwdx = GradCentrT(Half*(wCL(ix-1,iy,iz,1)+wCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                        ,w &
                        ,Half*(wCL(ix+1,iy,iz,1)+wCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps)  &
                        ,DH(ix-1,iy,iz,1),DH(ix,iy,iz,1),DH(ix+1,iy,iz,1)                    &
                        ,FWe,FE,VW,VCe,VE)
        dwdy = GradCentrT(Half*(wCL(ix,iy-1,iz,1)+wCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                        ,w &
                        ,Half*(wCL(ix,iy+1,iz,1)+wCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps)  &
                        ,DH(ix,iy-1,iz,1),DH(ix,iy,iz,1),DH(ix,iy+1,iz,1)                    &
                        ,FS,FN,VS,VCe,VN)
        dwdz = GradCentrT(Half*(wCL(ix,iy,iz-1,1)+wCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                        ,w &
                        ,Half*(wCL(ix,iy,iz+1,1)+wCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps)  &
                        ,DV(ix,iy,iz-1,1),DV(ix,iy,iz,1),DV(ix,iy,iz+1,1)                    &
                        ,FB,FT,VB,VCe,VT)
        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                         ,thC &
                         ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                         ,FB,FT,VB,VCe,VT)
!	Constants
!	c0 = 0.32d0
!	c2 = 0.43d0
!        c3 = 0.067d0+0.18d0*LenC(ix,iy,iz,1)/((two*VolC(ix,iy,iz)*(One/(FU(ix+1,iy,iz)+FU(ix,iy,iz)+Eps) &
!             +One/(FV(ix,iy+1,iz)+FV(ix,iy,iz)+Eps)+One/(FW(ix,iy,iz+1)+FW(ix,iy,iz)+Eps)))/Three+Eps)
        c3 = 2*Cmy0**0.75d0

!	Brunt-Vaisala-Frequency N**2
        N2 = Grav/(thC+Eps)*dthdz

!	Calculate height z
        IF (FW(ix,iy,iz-1)>FW(ix,iy,iz)) THEN
          dzLoc = MIN(VolC(ix,iy,iz)/(FW(ix,iy,iz-1)+Eps),dz(iz))
          zPL   = zP(iz-1)+0.5e0*dzLoc
        ELSE
          dzLoc = MIN(VolC(ix,iy,iz)/(FW(ix,iy,iz)+Eps),dz(iz))
          zPL   = zP(iz)-0.5e0*dzLoc
        END IF
!       Mean grid resolution ds
        ds = (two*VolC(ix,iy,iz)*(One/(FU(ix+1,iy,iz)+FU(ix,iy,iz)+Eps) &
                                 +One/(FV(ix,iy+1,iz)+FV(ix,iy,iz)+Eps) &
                                 +One/(FW(ix,iy,iz+1)+FW(ix,iy,iz)+Eps)))/Three

!       Turbulent Length Scale
        IF (N2>0.0d0) THEN
          Len0s = 0.54d0*SQRT(tkeG/(rhoC*ABS(N2)+Eps))
        ELSE
          Len0s=ds
        END IF 
!        Len0s = 0.54d0*SQRT((tkeHC+tkeVC)/(rhoC*N2+Eps))

!       Equilibrium turbulent Length Scale Len0
         Len0 = Rho(ix,iy,iz,1)    &
                   *MIN(0.67d0*zPL &
                   ,ds             &
                   ,Len0s)
!        Len0 = Rho(ix,iy,iz,1)    &
!                  *MIN(Karm*(zPL+zRauh)*(One-Karm**Two*Cmy0**0.25d0/c2) &
!                  ,ds             &
!                  ,Len0s)
!       Len0 = Rho(ix,iy,iz,1)*Karm*(zPL+zRauh)*(One-Karm**Two*Cmy0**0.25d0/c2)

        LenRhs(ix,iy,iz,1) = LenRhs(ix,iy,iz,1)-c2*SQRT(rhoC*tkeG) &
                             *(LenC(ix,iy,iz,1)-Len0)/(LenC(ix,iy,iz,1)+Eps)

!	Shear Production
        SH  = two*(DH(ix,iy,iz,1)*(dudx*dudx+dvdx*dvdx+dudy*dudy+dvdy*dvdy)+DV(ix,iy,iz,1)* &  !!!c3=2*Cmy0**0.75
                                  (dudz*dudz+dvdz*dvdz))
        SV  = two*(DH(ix,iy,iz,1)*(dwdx*dwdx+dwdy*dwdy)+DV(ix,iy,iz,1)*(dwdz*dwdz))

!	Production of Buoyancy
!        Bou = two*(One/PrandtlDis)*c0*LenC(ix,iy,iz,1)*SQRT(3.0d0*ABS(tkeVC))*N2
!        Bou = two*(One/PrandtlDis)*c0*LenC(ix,iy,iz,1)*SQRT(3.0d0*tkeVC/(rhoC+Eps))*N2
        Bou = -two*(One/PrandtlDis)*DV(ix,iy,iz,1)*N2 
        IF (-Bou/(SH+SV)>0.25d0) THEN
          SH=0.0d0
          SV=0.0d0
          Bou=0.0d0
        END IF

!	Redistribution
        RH  = c2*SQRT(rhoC)*SQRT(tkeG)                        &
              *(tkeHC-Two*tkeG/Three)/(LenC(ix,iy,iz,1)+Eps)
        RV  = c2*SQRT(rhoC)*SQRT(tkeG)                        &
              *(tkeVC-    tkeG/Three)/(LenC(ix,iy,iz,1)+Eps)

!	Dissipation
        DissiH = c3*SQRT(rhoC)*SQRT(tkeG)*tkeHC &
                 /(LenC(ix,iy,iz,1)+Eps)
        DissiV = c3*SQRT(rhoC)*SQRT(tkeG)*tkeVC &
                 /(LenC(ix,iy,iz,1)+Eps)

!	RHS of TKE-Equation
        tkeHRhs(ix,iy,iz,1) = tkeHRhs(ix,iy,iz,1) &
                              +SH-RH-DissiH
        tkeVRhs(ix,iy,iz,1) = tkeVRhs(ix,iy,iz,1) &
                              +SV+Bou-RV-DissiV
      END DO
    END DO
  END DO

END SUBROUTINE TkeHVLenCompute

SUBROUTINE JacTkeHVLenCompute

INTEGER :: ix,iy,iz
  REAL(RealKind) :: FB,FT
  REAL(RealKind) :: VB,VT,VCe
  REAL(RealKind) :: Vol
  REAL(RealKind) :: dthdz,thC,rhoC,tkeHC,tkeVC,tkeG
  REAL(RealKind) :: N2,c3
  REAL(RealKind) :: zPL,dzLoc,Len0,Len0s,zRauh,ds
! Quellen und Senken der Tke-Dis-Gleichungen
! Zentrale Gradienten linear interpoliert zwischen den primaeren Gradienten

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       Zentrale Groessen der betrachteten Box
        VCe   = VolC(ix  ,iy  ,iz  )
        VB    = VolC(ix  ,iy  ,iz-1)
        VT    = VolC(ix  ,iy  ,iz+1  )   
        FB    = FW(ix  ,iy  ,iz-1)
        FT    = FW(ix  ,iy  ,iz  )
        rhoC  = Rho(ix,iy,iz,1)
        thC   = th(ix,iy,iz,1)/(rhoC+Eps)
        tkeHC = tkeH(ix,iy,iz,1)
        tkeVC = tkeV(ix,iy,iz,1)
        tkeG  = tkeH(ix,iy,iz,1)+tkeV(ix,iy,iz,1)
        dzLoc = MaxVal(dz)
        zRauh = 0.1d0

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                         ,thC                                      &
                         ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                         ,FB,FT,VB,VCe,VT)

!	Constants
!	c0 = 0.32d0
!	c2 = 0.43d0
!        c3 = 0.067d0+0.18d0*LenC(ix,iy,iz,1)/((two*VolC(ix,iy,iz)*(One/(FU(ix+1,iy,iz)+FU(ix,iy,iz)+Eps) &
!             +One/(FV(ix,iy+1,iz)+FV(ix,iy,iz)+Eps)+One/(FW(ix,iy,iz+1)+FW(ix,iy,iz)+Eps)))/Three+Eps)
        c3 = 2*Cmy0**0.75d0

!	Brunt-Vasala Frequency N**2
        N2 = Grav/(thC+Eps)*dthdz

!	Calculate height z
        IF (FW(ix,iy,iz-1)>FW(ix,iy,iz)) THEN
          dzLoc = MIN(VolC(ix,iy,iz)/(FW(ix,iy,iz-1)+Eps),dz(iz))
          zPL   = zP(iz-1)+0.5e0*dzLoc
        ELSE
          dzLoc = MIN(VolC(ix,iy,iz)/(FW(ix,iy,iz)+Eps),dz(iz))
          zPL   = zP(iz)-0.5e0*dzLoc
        END IF
!       Mean grid resolution ds
        ds = (two*VolC(ix,iy,iz)*(One/(FU(ix+1,iy,iz)+FU(ix,iy,iz)+Eps) &
                                 +One/(FV(ix,iy+1,iz)+FV(ix,iy,iz)+Eps) &
                                 +One/(FW(ix,iy,iz+1)+FW(ix,iy,iz)+Eps)))/Three

!       Turbulent Length Scale                          !!!FILAUS
        IF (N2>0.0d0) THEN
          Len0s = 0.54d0*SQRT(tkeG/(rhoC*ABS(N2)+Eps))
        ELSE
          Len0s=ds
        END IF 
        Len0s = 0.54d0*SQRT(tkeG/(rhoC*ABS(N2)+Eps))
!        Len0s = 0.54d0*SQRT((tkeHC+tkeVC)/(rhoC*N2+Eps))

!       Equilibrium turbulent Length Scale Len0
         Len0 = Rho(ix,iy,iz,1)    &
                   *MIN(0.67d0*zPL &
                   ,ds             &
                   ,Len0s)
!       Len0 = Rho(ix,iy,iz,1)*Karm*(zPL+zRauh)*(One-Karm**Two*Cmy0**0.25d0/c2)

!       ---Derivatives tkeHRhs tkeH,tkeV,Len---
!                                     tkeHRhs(ix,iy,iz,1) = tkeHRhs(ix,iy,iz,1) &
!                                                           +SH-RH-DissiH	

        AS(IndexMet(tkeHPosJac,tkeHPosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkeHPosJac,tkeHPosJac))%c(ix,iy,iz,1)     &
                                                            -c2*SQRT(rhoC)*SQRT(tkeG)                             &
                                                             *((0.5d0*tkeHC-tkeVC)/(tkeG+Eps) + One)              &
                                                             /(Three*LenC(ix,iy,iz,1)+Eps)                        &
                                                            -c3*SQRT(rhoC)*SQRT(tkeG)                             &
                                                             *(0.5d0*tkeHC/(tkeG+Eps) + One)                      &
                                                             /(LenC(ix,iy,iz,1)+Eps)
        AS(IndexMet(tkeHPosJac,tkeVPosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkeHPosJac,tkeVPosJac))%c(ix,iy,iz,1)     &
                                                            -c2*SQRT(rhoC)*SQRT(tkeG)                             &
                                                             *((0.5d0*tkeHC-tkeVC)/(tkeG+Eps) - Two)              &
                                                             /(Three*LenC(ix,iy,iz,1)+Eps)                        &
                                                            -c3*SQRT(rhoC)*tkeHC                                  &
                                                             /(Two*SQRT(tkeG)*LenC(ix,iy,iz,1)+Eps)
        AS(IndexMet(tkeHPosJac,LenPosJac))%c(ix,iy,iz,1)  = AS(IndexMet(tkeHPosJac,LenPosJac))%c(ix,iy,iz,1)      &
                                                            +c2*SQRT(rhoC)*SQRT(tkeG)                             &
                                                             *(tkeHC-Two*tkeG/Three)                              &
                                                             /((LenC(ix,iy,iz,1))**Two+Eps)                       &
                                                            +c3*SQRT(rhoC)*SQRT(tkeG)*tkeHC                       &
                                                             /((LenC(ix,iy,iz,1))**Two+Eps)

!       ---Derivatives tkeVRhs tkeH,tkeV,Len---
!                                     tkeVRhs(ix,iy,iz,1) = tkeVRhs(ix,iy,iz,1) &
!                                                           +Bou-RV-DissiV

        AS(IndexMet(tkeVPosJac,tkeHPosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkeVPosJac,tkeHPosJac))%c(ix,iy,iz,1)  &
                                                            -c2*SQRT(rhoC)*SQRT(tkeG)                          &
                                                             *((tkeVC-0.5d0*tkeHC)/(tkeG+Eps) - One)           &
                                                             /(Three*LenC(ix,iy,iz,1)+Eps)                     &
                                                            -c3*SQRT(rhoC)*tkeVC                               &
                                                             /(Two*SQRT(tkeG)*LenC(ix,iy,iz,1)+Eps)
        AS(IndexMet(tkeVPosJac,tkeVPosJac))%c(ix,iy,iz,1) = AS(IndexMet(tkeVPosJac,tkeVPosJac))%c(ix,iy,iz,1)  &
!                                                           +Three*c0*N2*LenC(ix,iy,iz,1)                      &
!                                                            /(PrandtlDis*SQRT(Three*tkeVC/(rhoC+Eps))+Eps)    &
                                                            -c2*SQRT(rhoC)*SQRT(tkeG)                          &
                                                             *((tkeVC-0.5d0*tkeHC)/(tkeG+Eps) + Two)           &
                                                             /(Three*LenC(ix,iy,iz,1)+Eps)                     &
                                                            -c3*SQRT(rhoC)*SQRT(tkeG)                          &
                                                             *(0.5d0*tkeVC/(tkeG+Eps) + One)                   &
                                                             /(LenC(ix,iy,iz,1)+Eps)
        AS(IndexMet(tkeVPosJac,LenPosJac))%c(ix,iy,iz,1)  = AS(IndexMet(tkeVPosJac,LenPosJac))%c(ix,iy,iz,1)   &
!                                                           +Two*c0*N2*SQRT(Three*tkeVC/(rhoC+Eps))            &
!                                                            /(PrandtlDis+Eps)                                 &
                                                            +c2*SQRT(rhoC)*SQRT(tkeG)                          &
                                                             *(tkeVC-tkeG/Three)/((LenC(ix,iy,iz,1))**Two+Eps) &
                                                            +c3*SQRT(rhoC)*SQRT(tkeG)*tkeVC                    &
                                                             /((LenC(ix,iy,iz,1))**Two+Eps)

!       ---Derivatives LenRhs tkeH,tkeV,Len---
!                                      LenRhs(ix,iy,iz,1) = LenRhs(ix,iy,iz,1)-c2*SQRT(rhoC*(tkeH(ix,iy,iz,1)  &
!	                                                    +tkeV(ix,iy,iz,1)))*(LenC(ix,iy,iz,1)-Len0)         &
!                                                            /(LenC(ix,iy,iz,1)+Eps)

        AS(IndexMet(LenPosJac,TkeHPosJac))%c(ix,iy,iz,1)  = AS(IndexMet(LenPosJac,TkeHPosJac))%c(ix,iy,iz,1) &
                                                            -0.5d0*c2*SQRT(rhoC)/(SQRT(tkeG)+Eps)            &
                                                            *(One-Len0/(LenC(ix,iy,iz,1)+Eps))
        AS(IndexMet(LenPosJac,TkeVPosJac))%c(ix,iy,iz,1)  = AS(IndexMet(LenPosJac,TkeVPosJac))%c(ix,iy,iz,1) &
                                                            -0.5d0*c2*SQRT(rhoC)/(SQRT(tkeG)+Eps)            &
                                                            *(One-Len0/(LenC(ix,iy,iz,1)+Eps))
        AS(IndexMet(LenPosJac,LenPosJac))%c(ix,iy,iz,1)   = AS(IndexMet(LenPosJac,LenPosJac))%c(ix,iy,iz,1)  &
                                                            -c2*SQRT(rhoC)*SQRT(tkeG)*Len0                   &
                                                             /(LenC(ix,iy,iz,1)**Two+Eps)

      END DO
    END DO
  END DO
END SUBROUTINE JacTkeHVLenCompute

SUBROUTINE ReadBaum(FileName)

! Analogous to ReadStreet

  CHARACTER(*) :: FileName
  CHARACTER(80) :: mi_wgts2
  INTEGER :: ix,iy,iz
  INTEGER :: NumPoints
  REAL(RealKind) :: c
  INTEGER :: ibInput,j
  CHARACTER(300) :: Line

  mi_wgts2 = TRIM(FileName(1:INDEX(FileName,'.grid')-1))//'.Tree'
  OPEN(UNIT=InputUnit,FILE=mi_wgts2,STATUS='old')
  ALLOCATE(PointTree%TreePointBlock(nbLoc))

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    S:DO
      READ(InputUnit,'(a300)') Line
      IF (INDEX(Line,'#.... Block')>0) THEN
        READ(Line(INDEX(Line,'#.... Block')+LEN('#.... Block'):),*) ibInput,NumPoints ! the indices of blocks
        IF (ib==ibInput) THEN 
          WRITE(*,*) 'Input Mask',ib,ibInput,NumPoints
          PointTree%TreePointBlock(ibLoc)%NumTreePoint=NumPoints
          ALLOCATE(PointTree%TreePointBlock(ibLoc)%TreePoint(NumPoints))
          DO j=1,NumPoints
            READ(InputUnit,*) ix,iy,iz
            READ(InputUnit,*) c
            PointTree%TreePointBlock(ibLoc)%TreePoint(j)%c=c
            PointTree%TreePointBlock(ibLoc)%TreePoint(j)%ix=ix
            PointTree%TreePointBlock(ibLoc)%TreePoint(j)%iy=iy
            PointTree%TreePointBlock(ibLoc)%TreePoint(j)%iz=iz
            IF (VolC(ix,iy,iz)<=Zero) THEN
              WRITE(*,*) '*** Baum in',i,ix,iy,iz,'  cancelled (closed cell) ***'
              PointTree%TreePointBlock(ibLoc)%TreePoint(j)%c=Zero
            END IF
          END DO
          EXIT S
        END IF   ! ibInput
      END IF   ! ibInput
    END DO S
  END DO
  CLOSE(UNIT=InputUnit)

END SUBROUTINE ReadBaum

END MODULE Turbulence_Mod
