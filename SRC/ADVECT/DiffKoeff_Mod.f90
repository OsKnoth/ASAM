MODULE DiffKoeff_Mod
  USE DataType_Mod
  USE Domain_Mod
  USE Floor_Mod
  USE Physics_Mod
  USE Operator_Mod
  USE Parameter_mod
  USE Output_Mod

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

SUBROUTINE DiffKoeffSelect(Vector,UVec)
!
!==================================================
!----  Diagnostic Computation (from Velocity Values)
!==================================================
!
  IMPLICIT NONE

  TYPE (Vector4Cell_T) :: Vector(:)
  TYPE (Vector4Cell_T), OPTIONAL :: UVec(:)

  IF (Diffusion) THEN
    IF (DynSmag) THEN
      CALL DFromDynSmag(Vector,UVec)
    ELSE
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        Rho=>Vector(ibLoc)%Vec(RhoPos)%c
        IF (TkeDis) THEN
          D=>DiffKoeff(ibLoc)%c
          tke=>Vector(ibLoc)%Vec(tkePos)%c
          dis=>Vector(ibLoc)%Vec(disPos)%c
          CALL DFromTkeDisCompute
        ELSE IF (TkeDisRich) THEN
          D=>DiffKoeff(ibLoc)%c
          tke=>Vector(ibLoc)%Vec(tkePos)%c
          dis=>Vector(ibLoc)%Vec(disPos)%c
          CALL DFromTkeDisCompute
        ELSE IF (TkeOme) THEN
          D=>DiffKoeff(ibLoc)%c
          tke=>Vector(ibLoc)%Vec(tkePos)%c
          ome=>Vector(ibLoc)%Vec(omePos)%c
          CALL DFromTkeOmeCompute
        ELSE IF (TkeSGS) THEN
          th=>Vector(ibLoc)%Vec(thPos)%c
          tke=>Vector(ibLoc)%Vec(tkePos)%c
          IF (PRESENT(UVec)) THEN
            uCL=>UVec(ibLoc)%Vec(uPosL)%c
            vCL=>UVec(ibLoc)%Vec(vPosL)%c
            wCL=>UVec(ibLoc)%Vec(wPosL)%c
            uCR=>UVec(ibLoc)%Vec(uPosR)%c
            vCR=>UVec(ibLoc)%Vec(vPosR)%c
            wCR=>UVec(ibLoc)%Vec(wPosR)%c
          ELSE
            uCL=>Vector(ibLoc)%Vec(uPosL)%c
            vCL=>Vector(ibLoc)%Vec(vPosL)%c
            wCL=>Vector(ibLoc)%Vec(wPosL)%c
            uCR=>Vector(ibLoc)%Vec(uPosR)%c
            vCR=>Vector(ibLoc)%Vec(vPosR)%c
            wCR=>Vector(ibLoc)%Vec(wPosR)%c
          END IF
          uC=>uCell(ibLoc)%c
          vC=>vCell(ibLoc)%c
          wC=>wCell(ibLoc)%c
          DPot=>DiffPotKoeff(ibLoc)%c
          DMom=>DiffMomKoeff(ibLoc)%c
          LenC=>LenKoeff(ibLoc)%c
          CALL DFromTkeSGSCompute
        ELSE IF (TkeSmag) THEN
          th=>Vector(ibLoc)%Vec(thPos)%c
          uC=>uCell(ibLoc)%c
          vC=>vCell(ibLoc)%c
          wC=>wCell(ibLoc)%c
          DPotH=>DiffPotHKoeff(ibLoc)%c
          DPotV=>DiffPotVKoeff(ibLoc)%c
          DMomH=>DiffMomHKoeff(ibLoc)%c
          DMomV=>DiffMomVKoeff(ibLoc)%c
          LenC=>LenKoeff(ibLoc)%c
          CALL DFromTkeSmagCompute
        ELSE IF (TkeLen) THEN
          th=>Vector(ibLoc)%Vec(thPos)%c
          tke=>Vector(ibLoc)%Vec(tkePos)%c
          uCL=>Vector(ibLoc)%Vec(uPosL)%c
          vCL=>Vector(ibLoc)%Vec(vPosL)%c
          wCL=>Vector(ibLoc)%Vec(wPosL)%c
          uCR=>Vector(ibLoc)%Vec(uPosR)%c
          vCR=>Vector(ibLoc)%Vec(vPosR)%c
          wCR=>Vector(ibLoc)%Vec(wPosR)%c
          uC=>uCell(ibLoc)%c
          vC=>vCell(ibLoc)%c
          wC=>wCell(ibLoc)%c
          DPot=>DiffPotKoeff(ibLoc)%c
          DMom=>DiffMomKoeff(ibLoc)%c
          LenC=>Vector(ibLoc)%Vec(LenPos)%c
          CALL DFromTkeLenCompute
        ELSE IF (TkeHVLen) THEN
          DH=>DiffHKoeff(ibLoc)%c
          DV=>DiffVKoeff(ibLoc)%c
          tkeH=>Vector(ibLoc)%Vec(tkeHPos)%c
          tkeV=>Vector(ibLoc)%Vec(tkeVPos)%c
          LenC=>Vector(ibLoc)%Vec(LenPos)%c
          CALL DFromTkeHVLenCompute
        ELSE IF (NoTke) THEN
          th=>Vector(ibLoc)%Vec(thPos)%c
          tke=>Vector(ibLoc)%Vec(tkePos)%c
          uCL=>Vector(ibLoc)%Vec(uPosL)%c
          vCL=>Vector(ibLoc)%Vec(vPosL)%c
          wCL=>Vector(ibLoc)%Vec(wPosL)%c
          uCR=>Vector(ibLoc)%Vec(uPosR)%c
          vCR=>Vector(ibLoc)%Vec(vPosR)%c
          wCR=>Vector(ibLoc)%Vec(wPosR)%c
          uC=>uCell(ibLoc)%c
          vC=>vCell(ibLoc)%c
          wC=>wCell(ibLoc)%c
          DPot=>DiffPotKoeff(ibLoc)%c
          DMom=>DiffMomKoeff(ibLoc)%c
          LenC=>LenKoeff(ibLoc)%c
          CALL DCompute
        ELSE
          DiffKoeff(ibLoc)%c=Rho*DiffMin
        END IF
      END DO
    END IF
    IF (TkeDis.OR.TkeDisRich) THEN
      CALL ExchangeCell(DiffKoeff)
    ELSE IF (TkeSGS.OR.TkeLen.OR.NoTke) THEN
      CALL ExchangeCell(DiffPotKoeff)
      CALL ExchangeCell(DiffMomKoeff)
    ELSE IF (TkeSmag) THEN
      CALL ExchangeCell(DiffPotHKoeff)
      CALL ExchangeCell(DiffPotVKoeff)
      CALL ExchangeCell(DiffMomHKoeff)
      CALL ExchangeCell(DiffMomVKoeff)
    ELSE IF (DynSmag) THEN
      CALL ExchangeCell(DiffPotHKoeff)
      CALL ExchangeCell(DiffPotVKoeff)
      CALL ExchangeCell(DiffMomHKoeff)
      CALL ExchangeCell(DiffMomVKoeff)
    ELSE IF (TkeHVLen) THEN
      CALL ExchangeCell(DiffHKoeff)
      CALL ExchangeCell(DiffVKoeff)
    ELSE 
      CALL ExchangeCell(DiffKoeff)
    END IF
  END IF

!-----------------------------------------------------------
END SUBROUTINE DiffKoeffSelect

SUBROUTINE DCompute 
!
!==================================================
!----  Diffusion Coefficient without TKE
!==================================================
! 
  IMPLICIT NONE

!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: S,O
  REAL(RealKind) :: N2
  REAL(RealKind) :: Cm,Ch1,Ch2,Ceps1,Ceps2,Ccs,delta
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC
  REAL(RealKind) :: RhoC,LenLoc

!--- local constants
  dxLoc = MaxVal(dx(ix0+1:ix1))
  dyLoc = MaxVal(dy(iy0+1:iy1))
  dzLoc = MaxVal(dz(iz0+1:iz1))
  
  IF (ix1-(ix0+1)==0) THEN
    dxLoc=dyLoc
  ELSE IF (iy1-(iy0+1)==0) THEN 
    dyLoc=dxLoc
  END IF

  Cm    = 0.0856d0
  Ch1   = Cm
  Ch2   = 0.1184d0
  Ceps1 = 0.1911d0
  Ceps2 = 0.6539d0
  
  Ccs    = (Cm**Three/(Ceps1+Ceps2))**(One/Four)
  delta  = (dxLoc*dyLoc*dzLoc)**(One/Three)
  
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe  = VolC(ix  ,iy  ,iz  )
        VB   = VolC(ix  ,iy  ,iz-1)
        VT   = VolC(ix  ,iy  ,iz+1)
        FB   = FW(ix  ,iy  ,iz-1)
        FT   = FW(ix  ,iy  ,iz  )
        RhoC = Rho(ix,iy,iz,1)
        thC  = th(ix,iy,iz,1)/(RhoC+Eps)

        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),              &
                         FV(ix:ix,iy-1:iy,iz:iz),              &
                         FW(ix:ix,iy:iy,iz-1:iz),              &
                         S,O)

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                         ,thC                                      &
                         ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                         ,FB,FT,VB,VCe,VT)

!       Richardson Number
        N2  = (Grav/(thC+Eps))*dthdz

!       Local Length Scale
          LenLoc=Ccs*MIN(delta,dzLoc)

!       Diffusion Coefficients Compute
        DMom(ix,iy,iz,1) = MIN(RhoC*LenLoc**Two*SQRT(MAX(S*S-(Ch1+Ch2)*N2/Cm,Zero)),RhoC*DiffMax)
        DPot(ix,iy,iz,1) = (Ch1+Ch2)/Cm*DMom(ix,iy,iz,1)
      END DO
    END DO
  END DO

!-----------------------------------------------------------
END SUBROUTINE DCompute

SUBROUTINE DFromTkeDisCompute
!
!==================================================
!----  Diffusion Coefficient from Tke*Tke/Dis
!==================================================
!
  IMPLICIT NONE

!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: tkeC,disC,rhoC

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        tkeC = HALF*(ABS(tke(ix,iy,iz,1))+tke(ix,iy,iz,1))
        disC = MAX(dis(ix,iy,iz,1),Zero)
        rhoC = Rho(ix,iy,iz,1)
        IF (dis(ix,iy,iz,1) > 1.d-8) THEN
          D(ix,iy,iz,1) = MIN(MAX(Cmy0*tkeC*tkeC/(disC+Eps),rhoC*DiffMin),rhoC*DiffMax) &
                          *VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        ELSE
          D(ix,iy,iz,1) = RhoC*DiffMin
        END IF 
      END DO
    END DO
  END DO

END SUBROUTINE DFromTkeDisCompute
!-----------------------------------------------------------

SUBROUTINE DFromTkeOmeCompute
  IMPLICIT NONE
!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: tkeC,omeC,rhoC

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        tkeC = HALF*(ABS(tke(ix,iy,iz,1))+tke(ix,iy,iz,1))
        omeC = MAX(ome(ix,iy,iz,1),Zero)
        rhoC = Rho(ix,iy,iz,1)
        D(ix,iy,iz,1) = MIN(MAX(Cmiu*rhoC**1.5d0*SQRT(tkeC)/(omeC+Eps),rhoC*DiffMin),rhoC*DiffMax) &
                          *VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE DFromTkeOmeCompute
!-----------------------------------------------------------

SUBROUTINE DFromTkeSGSCompute 
!
!==================================================
!----  Diffusion Coefficient depending on TKE
!==================================================
!
  IMPLICIT NONE

!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: S,O
  REAL(RealKind) :: N2,Rich,Pr
  REAL(RealKind) :: FM,FH,PhiM
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC
  REAL(RealKind) :: tkeC,rhoC
  REAL(RealKind) :: LenLoc,Brunt

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe  = VolC(ix  ,iy  ,iz  )
        VB   = VolC(ix  ,iy  ,iz-1)
        VT   = VolC(ix  ,iy  ,iz+1)
        FB   = FW(ix  ,iy  ,iz-1)
        FT   = FW(ix  ,iy  ,iz  )
        rhoC = Rho(ix,iy,iz,1)
        thC  = th(ix,iy,iz,1)/(RhoC+Eps)
        tkeC = ABS(tke(ix,iy,iz,1))
        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),              &
                         FV(ix:ix,iy-1:iy,iz:iz),              &
                         FW(ix:ix,iy:iy,iz-1:iz),              &
                         S,O)

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                         ,thC                                      &
                         ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                         ,FB,FT,VB,VCe,VT)
 
!       Richardson Number
        Brunt  = (Grav/(thC+Eps))*dthdz
        Rich  = (Grav/(thC+Eps))*dthdz/(S*O+1.d-3)
  
!       Parametric Functions
        IF (Rich < Zero) THEN
!         Unstable Conditions
          FM = (One-16.0d0*Rich)**0.5d0
          FH = (One-40.0d0*Rich)**0.5d0/(Prn+Eps)
        ELSE IF (Rich >= RichC) THEN
!         Very Stable Conditions
          FM = 0.0d0 
          FH = 0.0d0 
        ELSE  
!         Stable Conditions
          FM = (One-Rich/RichC)**Four
          FH = (One-Rich/RichC)**Four*(One-1.2d0*Rich)/(Prn+Eps)
        END IF

!       Local Length Scale
        IF (Brunt>Zero) THEN
          LenLoc=MIN(LenScaleGeom,0.76*SQRT(tkeC/(RhoC+Eps)/Brunt))
        ELSE 
          LenLoc=LenScaleGeom
        END IF

!       Turbulent Kinetic Energy and Diffusion Coefficients
        IF (Rich >= RichC) THEN
          PhiM             = 0.0d0
          DMom(ix,iy,iz,1) = RhoC*DiffMin
          DPot(ix,iy,iz,1) = RhoC*DiffMin
        ELSE
!         Prandtl Number
          Pr   = FM/( FH+Eps )
!         Stability Function for Momentum 
          PhiM = (Cs**(Four/Three)) * (Ceps**(One/Three)) * FM**(Two/Three) &
                 / (One-Rich/(Pr+Eps))**(One/Three)
!         Diffusion Coefficients Compute
          DMom(ix,iy,iz,1)=MAX(MIN(PhiM*LenLoc*SQRT(tkeC*rhoC),RhoC*DiffMax),RhoC*DiffMin)
          DPot(ix,iy,iz,1)=MAX(MIN(PhiM/(Pr+Eps)*LenLoc*SQRT(tkeC*rhoC),RhoC*DiffMax),RhoC*DiffMin)
        END IF

      END DO
    END DO
  END DO

!-----------------------------------------------------------
END SUBROUTINE DFromTkeSGSCompute

SUBROUTINE DFromTkeSmagCompute
!
!==================================================
!----  Diffusion Coefficient depending on Smagorinsky
!==================================================
!
  IMPLICIT NONE

!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: S,O
  REAL(RealKind) :: N2,Rich,Pr
  REAL(RealKind) :: FM,FH,PhiM
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC
  REAL(RealKind) :: tkeC,rhoC
  REAL(RealKind) :: LenLoc,Brunt,a1,a2
  REAL(RealKind) :: Csmag,alpha,phi,Lennull,LenScale
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc,dxdz,dydz,dzdx,dxdy,delta
  REAL(RealKind),DIMENSION(1:6) :: Sij

!  dxLoc = MAXVAL(dx(ix0+1:ix1))
!  dyLoc = MAXVAL(dy(iy0+1:iy1))
!  dzLoc = MAXVAL(dz(iz0+1:iz1))
!        a1   = (MIN(dxLoc,dyLoc,dzLoc))/ &
!               ((dxLoc+dyLoc+dzLoc)-MIN(dxLoc,dyLoc,dzLoc)-MAX(dxLoc,dyLoc,dzLoc))
!        a2   = ((dxLoc+dyLoc+dzLoc)-MIN(dxLoc,dyLoc,dzLoc)-MAX(dxLoc,dyLoc,dzLoc))/ &
!               (MAX(dxLoc,dyLoc,dzLoc))
!        LenLoc = (dxLoc*dyLoc*dzLoc)**(1.0d0/3.0d0)
!        LenLoc = LenLoc* &
!                 cosh(SQRT((4.0d0/27.0d0)*((log(a1))**Two-log(a1)*log(a2)+(log(a2))**Two)))

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe  = VolC(ix  ,iy  ,iz  )
        VB   = VolC(ix  ,iy  ,iz-1)
        VT   = VolC(ix  ,iy  ,iz+1)
        FB   = FW(ix  ,iy  ,iz-1)
        FT   = FW(ix  ,iy  ,iz  )
        rhoC = Rho(ix,iy,iz,1)
        thC  = th(ix,iy,iz,1)/(RhoC+Eps)
        LenLoc = (dx(ix)*dy(iy)*dz(iz))**(1.0d0/3.0d0)

        CALL StrainRateRho(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),              &
                         FV(ix:ix,iy-1:iy,iz:iz),              &
                         FW(ix:ix,iy:iy,iz-1:iz),              &
                         S,Sij)
!       dthdz = GradCentrFull(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
!                        ,thC                                      &
!                        ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
!                        ,dz(iz-1),dz(iz),dz(iz+1))
         dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                          ,thC                                      &
                          ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)
!       Richardson Number
        Brunt  = (Grav/(thC+Eps))*dthdz
        Rich  = (Grav/(thC+Eps))*dthdz/(S*S+1.d-3)

!    use this for depending LenScale on atmos. stability
!        IF (Rich<Zero) THEN
!!         Unstable Conditions
!          phi=(One-Four*alpha*Rich)**0.25d0
!          Lennull = LenLoc
!          LenScale=One/ (((((One-Rich)**0.25d0)*phi) / (0.4d0*(zP(iz)+0.1d0)))+((One) / (Lennull)))
!        ELSE IF (Rich >= RichC) THEN
!!         Very Stable Conditions
!          LenScale= LenLoc/LenLoc  ! eigentlich null
!        ELSE
!!         Neutral or slightly stable Conditions
!          Lennull = LenLoc * (One-alpha*Rich)
!          LenScale=One/ (((((One-Rich)**0.25d0)) / (0.4d0*(zP(iz)+0.1d0)))+((One) / (Lennull)))
!        END IF

!       Parametric Functions
        IF (Rich < Zero) THEN
!         Unstable Conditions
          FM = (One-16.0d0*Rich)**0.5d0
          FH = (One-40.0d0*Rich)**0.5d0/(Prn+Eps)
        ELSE IF (Rich >= RichC) THEN
!         Very Stable Conditions
          FM = 0.0d-5
          FH = 0.0d-5
        ELSE
!         Neutral or slightly Stable Conditions
          FM = (One-Rich/RichC)**Four
          FH = (One-Rich/RichC)**Four*(One-1.2d0*Rich)/(Prn+Eps)
        END IF

!       Prandtl Number
        Pr   = MAX((FM/(FH+Eps)),3.3d-1) 

!       Local Length Scale Wall function
!       Decreasing of LengthScale near the wall aus Kleissl Dis.
        LenScale=((1.0d0/(0.4d0*(zP(iz)+0.1d0))**One)+(1.0d0/(LenLoc**One)))**(-1.0d0/One)

!       Output as SpecialOutput   -Len_turb-
        LenC(ix,iy,iz,1)= LenScale

!       Diffusion Coefficients Compute
        DMomH(ix,iy,iz,1)=RhoC*((CsmagK*LenScale)**Two)*S*FM
        DMomV(ix,iy,iz,1)=RhoC*((CsmagK*LenScale)**Two)*S*FM
        DPotH(ix,iy,iz,1)=RhoC*((CsmagK*LenScale)**Two)*S*FH
        DPotV(ix,iy,iz,1)=RhoC*((CsmagK*LenScale)**Two)*S*FH 
!!         Very Stable Conditions                 ! eigentlich nur bei
!vertikalen D, denn der horizontale D ist nicht von der stabilitaet betroffen

      END DO
    END DO
  END DO

!-----------------------------------------------------------
END SUBROUTINE DFromTkeSmagCompute

SUBROUTINE DFromDynSmag(Vector,UVec) !!!marcelk

  IMPLICIT NONE
  TYPE (Vector4Cell_T) :: Vector(:)
  TYPE (Vector4Cell_T), OPTIONAL :: UVec(:)
!---  local variables
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: S
  REAL(RealKind) :: LenLoc
  REAL(RealKind) :: FM,FH,PhiM
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind),DIMENSION(1:6) :: Sij

  IF (InitDynSmag==0) THEN
    InitDynSmag=1
    CALL Allocate(VecAve,1,17)
  END IF
  CALL ExchangeCell(VecAve)
  CALL BoundaryZeroGrad(VecAve)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Rho=>Vector(ibLoc)%Vec(RhoPos)%c
    th=>Vector(ibLoc)%Vec(thPos)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c

    uHat => VecAve(ibLoc)%Vec(uHatP)%c
    vHat => VecAve(ibLoc)%Vec(vHatP)%c
    wHat => VecAve(ibLoc)%Vec(wHatP)%c
    uuHat=>VecAve(ibLoc)%Vec(uuHatP)%c
    vvHat=>VecAve(ibLoc)%Vec(vvHatP)%c
    wwHat=>VecAve(ibLoc)%Vec(wwHatP)%c
    uvHat=>VecAve(ibLoc)%Vec(uvHatP)%c
    uwHat=>VecAve(ibLoc)%Vec(uwHatP)%c
    vwHat=>VecAve(ibLoc)%Vec(vwHatP)%c
    Tau11Hat=>VecAve(ibLoc)%Vec(Tau11HatP)%c
    Tau22Hat=>VecAve(ibLoc)%Vec(Tau22HatP)%c
    Tau33Hat=>VecAve(ibLoc)%Vec(Tau33HatP)%c
    Tau12Hat=>VecAve(ibLoc)%Vec(Tau12HatP)%c
    Tau13Hat=>VecAve(ibLoc)%Vec(Tau13HatP)%c
    Tau23Hat=>VecAve(ibLoc)%Vec(Tau23HatP)%c

    CALL FilterFeinCompute
  END DO

  CALL ExchangeCell(VecAve)
  CALL BoundaryZeroGrad(VecAve)
! ---------------------Filtern 1
! -------------------x-direction
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    uHat => VecAve(ibLoc)%Vec(uHatP)%c
    vHat => VecAve(ibLoc)%Vec(vHatP)%c
    wHat => VecAve(ibLoc)%Vec(wHatP)%c
    uuHat=>VecAve(ibLoc)%Vec(uuHatP)%c
    vvHat=>VecAve(ibLoc)%Vec(vvHatP)%c
    wwHat=>VecAve(ibLoc)%Vec(wwHatP)%c
    uvHat=>VecAve(ibLoc)%Vec(uvHatP)%c
    uwHat=>VecAve(ibLoc)%Vec(uwHatP)%c
    vwHat=>VecAve(ibLoc)%Vec(vwHatP)%c
    Tau11Hat=>VecAve(ibLoc)%Vec(Tau11HatP)%c
    Tau22Hat=>VecAve(ibLoc)%Vec(Tau22HatP)%c
    Tau33Hat=>VecAve(ibLoc)%Vec(Tau33HatP)%c
    Tau12Hat=>VecAve(ibLoc)%Vec(Tau12HatP)%c
    Tau13Hat=>VecAve(ibLoc)%Vec(Tau13HatP)%c
    Tau23Hat=>VecAve(ibLoc)%Vec(Tau23HatP)%c
    DO ic=LBOUND(VecAve(ibLoc)%Vec,1),UBOUND(VecAve(ibLoc)%Vec,1)
      CALL FilterX(VecAve(ibLoc)%Vec(ic)%c(:,:,:,1),VolC(:,:,:))
    END DO
  END DO
  CALL ExchangeCell(VecAve)
  CALL BoundaryZeroGrad(VecAve)
! -------------------y-direction
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    uHat => VecAve(ibLoc)%Vec(uHatP)%c
    vHat => VecAve(ibLoc)%Vec(vHatP)%c
    wHat => VecAve(ibLoc)%Vec(wHatP)%c
    uuHat=>VecAve(ibLoc)%Vec(uuHatP)%c
    vvHat=>VecAve(ibLoc)%Vec(vvHatP)%c
    wwHat=>VecAve(ibLoc)%Vec(wwHatP)%c
    uvHat=>VecAve(ibLoc)%Vec(uvHatP)%c
    uwHat=>VecAve(ibLoc)%Vec(uwHatP)%c
    vwHat=>VecAve(ibLoc)%Vec(vwHatP)%c
    Tau11Hat=>VecAve(ibLoc)%Vec(Tau11HatP)%c
    Tau22Hat=>VecAve(ibLoc)%Vec(Tau22HatP)%c
    Tau33Hat=>VecAve(ibLoc)%Vec(Tau33HatP)%c
    Tau12Hat=>VecAve(ibLoc)%Vec(Tau12HatP)%c
    Tau13Hat=>VecAve(ibLoc)%Vec(Tau13HatP)%c
    Tau23Hat=>VecAve(ibLoc)%Vec(Tau23HatP)%c
    DO ic=LBOUND(VecAve(ibLoc)%Vec,1),UBOUND(VecAve(ibLoc)%Vec,1)
      CALL FilterY(VecAve(ibLoc)%Vec(ic)%c(:,:,:,1),VolC(:,:,:))
    END DO
  END DO
  CALL ExchangeCell(VecAve)
  CALL BoundaryZeroGrad(VecAve)
! -------------------z-direction
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    uHat => VecAve(ibLoc)%Vec(uHatP)%c
    vHat => VecAve(ibLoc)%Vec(vHatP)%c
    wHat => VecAve(ibLoc)%Vec(wHatP)%c
    uuHat=>VecAve(ibLoc)%Vec(uuHatP)%c
    vvHat=>VecAve(ibLoc)%Vec(vvHatP)%c
    wwHat=>VecAve(ibLoc)%Vec(wwHatP)%c
    uvHat=>VecAve(ibLoc)%Vec(uvHatP)%c
    uwHat=>VecAve(ibLoc)%Vec(uwHatP)%c
    vwHat=>VecAve(ibLoc)%Vec(vwHatP)%c
    Tau11Hat=>VecAve(ibLoc)%Vec(Tau11HatP)%c
    Tau22Hat=>VecAve(ibLoc)%Vec(Tau22HatP)%c
    Tau33Hat=>VecAve(ibLoc)%Vec(Tau33HatP)%c
    Tau12Hat=>VecAve(ibLoc)%Vec(Tau12HatP)%c
    Tau13Hat=>VecAve(ibLoc)%Vec(Tau13HatP)%c
    Tau23Hat=>VecAve(ibLoc)%Vec(Tau23HatP)%c
    DO ic=LBOUND(VecAve(ibLoc)%Vec,1),UBOUND(VecAve(ibLoc)%Vec,1)
      CALL FilterZ(VecAve(ibLoc)%Vec(ic)%c(:,:,:,1),VolC(:,:,:))
    END DO
  END DO
  CALL BoundaryZeroGrad(VecAve)
  CALL ExchangeCell(VecAve)
! ---------------------Filtern 1 ENDE

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    th=>Vector(ibLoc)%Vec(thPos)%c
    Rho=>Vector(ibLoc)%Vec(RhoPos)%c

    LenC=>LenKoeff(ibLoc)%c
    uHat => VecAve(ibLoc)%Vec(uHatP)%c
    vHat => VecAve(ibLoc)%Vec(vHatP)%c
    wHat => VecAve(ibLoc)%Vec(wHatP)%c
    uuHat=>VecAve(ibLoc)%Vec(uuHatP)%c
    vvHat=>VecAve(ibLoc)%Vec(vvHatP)%c
    wwHat=>VecAve(ibLoc)%Vec(wwHatP)%c
    uvHat=>VecAve(ibLoc)%Vec(uvHatP)%c
    uwHat=>VecAve(ibLoc)%Vec(uwHatP)%c
    vwHat=>VecAve(ibLoc)%Vec(vwHatP)%c
    Tau11Hat=>VecAve(ibLoc)%Vec(Tau11HatP)%c
    Tau22Hat=>VecAve(ibLoc)%Vec(Tau22HatP)%c
    Tau33Hat=>VecAve(ibLoc)%Vec(Tau33HatP)%c
    Tau12Hat=>VecAve(ibLoc)%Vec(Tau12HatP)%c
    Tau13Hat=>VecAve(ibLoc)%Vec(Tau13HatP)%c
    Tau23Hat=>VecAve(ibLoc)%Vec(Tau23HatP)%c

    CALL  DynSmagKoeffCompute
  END DO
  CALL ExchangeCell(VecAve)
  CALL BoundaryZeroGrad(VecAve)

!  Volumen-Average for LM und MM for Csmag computation
! -------------------x-direction
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL FilterX3(VecAve(ibLoc)%Vec(16)%c(:,:,:,1),VolC(:,:,:))                                                                                                                                                              
    CALL FilterX3(VecAve(ibLoc)%Vec(17)%c(:,:,:,1),VolC(:,:,:)) 
  END DO
  CALL ExchangeCell(VecAve)
  CALL BoundaryZeroGrad(VecAve)
! -------------------y-direction
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL FilterY3(VecAve(ibLoc)%Vec(16)%c(:,:,:,1),VolC(:,:,:))                                                                                                                                                              
    CALL FilterY3(VecAve(ibLoc)%Vec(17)%c(:,:,:,1),VolC(:,:,:)) 
  END DO
  CALL ExchangeCell(VecAve)
  CALL BoundaryZeroGrad(VecAve)
! -------------------z-direction
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL FilterZ3(VecAve(ibLoc)%Vec(16)%c(:,:,:,1),VolC(:,:,:))                                                                                                                                                              
    CALL FilterZ3(VecAve(ibLoc)%Vec(17)%c(:,:,:,1),VolC(:,:,:)) 
  END DO
  CALL BoundaryZeroGrad(VecAve)
  CALL ExchangeCell(VecAve)
! ENDE Volumen-Average

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    th=>Vector(ibLoc)%Vec(thPos)%c
    Rho=>Vector(ibLoc)%Vec(RhoPos)%c
    uC=>uCell(ibLoc)%c
    vC=>vCell(ibLoc)%c
    wC=>wCell(ibLoc)%c
    DPotH=>DiffPotHKoeff(ibLoc)%c
    DPotV=>DiffPotVKoeff(ibLoc)%c
    DMomH=>DiffMomHKoeff(ibLoc)%c
    DMomV=>DiffMomVKoeff(ibLoc)%c

    LenC=>LenKoeff(ibLoc)%c
    uHat => VecAve(ibLoc)%Vec(uHatP)%c
    vHat => VecAve(ibLoc)%Vec(vHatP)%c
    wHat => VecAve(ibLoc)%Vec(wHatP)%c
    uuHat=>VecAve(ibLoc)%Vec(uuHatP)%c
    vvHat=>VecAve(ibLoc)%Vec(vvHatP)%c
    wwHat=>VecAve(ibLoc)%Vec(wwHatP)%c
    uvHat=>VecAve(ibLoc)%Vec(uvHatP)%c
    uwHat=>VecAve(ibLoc)%Vec(uwHatP)%c
    vwHat=>VecAve(ibLoc)%Vec(vwHatP)%c
    Tau11Hat=>VecAve(ibLoc)%Vec(Tau11HatP)%c
    Tau22Hat=>VecAve(ibLoc)%Vec(Tau22HatP)%c
    Tau33Hat=>VecAve(ibLoc)%Vec(Tau33HatP)%c
    Tau12Hat=>VecAve(ibLoc)%Vec(Tau12HatP)%c
    Tau13Hat=>VecAve(ibLoc)%Vec(Tau13HatP)%c
    Tau23Hat=>VecAve(ibLoc)%Vec(Tau23HatP)%c

    CALL  DFromDynSmagCompute
  END DO
  CALL ExchangeCell(DiffPotHKoeff)
  CALL ExchangeCell(DiffPotVKoeff)
  CALL ExchangeCell(DiffMomHKoeff)
  CALL ExchangeCell(DiffMomVKoeff)

END SUBROUTINE DFromDynSmag !!!marcelk

SUBROUTINE FilterFeinCompute

  INTEGER :: ix, iy, iz
  REAL(RealKind) :: S
  REAL(RealKind),DIMENSION(1:6) :: Sij
  REAL(RealKind) :: LenLoc,Brunt,LenScale,Rich
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC
  REAL(RealKind) :: tkeC,rhoC
  REAL(RealKind) :: a1,a2
  REAL(RealKind) :: alpha,phi,Lennull
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc,dxdz,dydz,dzdx,dxdy,delta

  uuHat = uC/(Rho+Eps)*uC/(Rho+Eps)
  vvHat = vC/(Rho+Eps)*vC/(Rho+Eps)
  wwHat = wC/(Rho+Eps)*wC/(Rho+Eps)
  uvHat = uC/(Rho+Eps)*vC/(Rho+Eps)
  uwHat = uC/(Rho+Eps)*wC/(Rho+Eps)
  vwHat = vC/(Rho+Eps)*wC/(Rho+Eps)
  uHat  = uC/(Rho+Eps)
  vHat  = vC/(Rho+Eps)
  wHat  = wC/(Rho+Eps)
  DO iz=iz0+1,iz1 ! u,v,w not yet filtered
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe  = VolC(ix  ,iy  ,iz  )
        VB   = VolC(ix  ,iy  ,iz-1)
        VT   = VolC(ix  ,iy  ,iz+1)
        FB   = FW(ix  ,iy  ,iz-1)
        FT   = FW(ix  ,iy  ,iz  )
        rhoC = Rho(ix,iy,iz,1)
        thC  = th(ix,iy,iz,1)/(RhoC+Eps)
!        LenLoc = (dx(ix)*dy(iy)*dz(iz))**(1.0d0/3.0d0)

        CALL StrainRateRho(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                    VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                    FU(ix-1:ix,iy:iy,iz:iz),              &
                    FV(ix:ix,iy-1:iy,iz:iz),              &
                    FW(ix:ix,iy:iy,iz-1:iz),              &
                    S,Sij)

!!        dthdz = GradCentrFull(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
!                         ,thC                                      &
!                         ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
!                         ,dz(iz-1),dz(iz),dz(iz+1))
         dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                          ,thC                                      &
                          ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)

!         Local Length Scale Wall function
!         Decreasing of LengthScale near the wall aus Kleissl Dis.
!          LenScale=((1.0d0/(0.4d0*(zP(iz)+0.1d0))**One)+(1.0d0/(LenLoc**One)))**(-1.0d0/One)

!        Tau11Hat(ix,iy,iz,1) = (LenScale**1.0)*ABS(S)*Sij(1)
!        Tau22Hat(ix,iy,iz,1) = (LenScale**1.0)*ABS(S)*Sij(2)
!        Tau33Hat(ix,iy,iz,1) = (LenScale**1.0)*ABS(S)*Sij(3)
!        Tau12Hat(ix,iy,iz,1) = (LenScale**1.0)*ABS(S)*Sij(4)
!        Tau13Hat(ix,iy,iz,1) = (LenScale**1.0)*ABS(S)*Sij(5)
!        Tau23Hat(ix,iy,iz,1) = (LenScale**1.0)*ABS(S)*Sij(6)
        Tau11Hat(ix,iy,iz,1) = ABS(S)*Sij(1)
        Tau22Hat(ix,iy,iz,1) = ABS(S)*Sij(2)
        Tau33Hat(ix,iy,iz,1) = ABS(S)*Sij(3)
        Tau12Hat(ix,iy,iz,1) = ABS(S)*Sij(4)
        Tau13Hat(ix,iy,iz,1) = ABS(S)*Sij(5)
        Tau23Hat(ix,iy,iz,1) = ABS(S)*Sij(6)
      END DO
    END DO
  END DO

END SUBROUTINE FilterFeinCompute

SUBROUTINE DynSmagKoeffCompute !!!marcelk nach Tejada-Martinez or Froehlich Buch
!
!==================================================
!----  Diffusion Coefficient depending on dynamic Smagorinsky
!==================================================
!
!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: S
  REAL(RealKind) :: N2,Rich,Pr
  REAL(RealKind) :: FM,FH,PhiM
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC
  REAL(RealKind) :: LenLoc,Brunt,LenScale
  REAL(RealKind) :: rhoC
  REAL(RealKind) :: alpha,phi,Lennull,LM,MM,QN,NN
  REAL(RealKind) :: BETA,Beta1,Csmag,Csmag2,CsmagZae,CsmagNen
  REAL(RealKind) :: T11,T22,T33,T12,T13,T23,a1,a2,S11
  REAL(RealKind) :: TT11,TT22,TT33,TT12,TT13,TT23,a,b,c
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc,dxdz,dydz,dzdx,dxdy,delta
  REAL(RealKind) :: aa1,bb1,cc1,dd1,ee1,aa2,bb2,cc2, dd2, ee2
  REAL(RealKind) :: P1,P2,P3,P4,P5,P0,PP,dPPdx
  REAL(RealKind),DIMENSION(1:6) :: Sij
  REAL(RealKind),DIMENSION(1:6) :: Lij
  REAL(RealKind),DIMENSION(1:6) :: Mij
  REAL(RealKind),DIMENSION(1:6) :: Qij
  REAL(RealKind),DIMENSION(1:6) :: Nij
  REAL(RealKind) :: x0, x1, epsilon 
  epsilon = 1.0d-10 

!  dxLoc = MAXVAL(dx(ix0+1:ix1))
!  dyLoc = MAXVAL(dy(iy0+1:iy1))
!  dzLoc = MAXVAL(dz(iz0+1:iz1))
!        a1   = (MIN(dxLoc,dyLoc,dzLoc))/ &
!               ((dxLoc+dyLoc+dzLoc)-MIN(dxLoc,dyLoc,dzLoc)-MAX(dxLoc,dyLoc,dzLoc))
!        a2   = ((dxLoc+dyLoc+dzLoc)-MIN(dxLoc,dyLoc,dzLoc)-MAX(dxLoc,dyLoc,dzLoc))/ &
!               (MAX(dxLoc,dyLoc,dzLoc))
!        LenLoc = LenLoc* &
!                 cosh(SQRT((4.0d0/27.0d0)*((log(a1))**Two-log(a1)*log(a2)+(log(a2))**Two)))

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe  = VolC(ix  ,iy  ,iz  )
        VB   = VolC(ix  ,iy  ,iz-1)
        VT   = VolC(ix  ,iy  ,iz+1)
        FB   = FW(ix  ,iy  ,iz-1)
        FT   = FW(ix  ,iy  ,iz  )
        rhoC = Rho(ix,iy,iz,1)
        thC  = th(ix,iy,iz,1)/(RhoC+Eps)
        LenLoc = (dx(ix)*dy(iy)*dz(iz))**(1.0d0/3.0d0)

    CALL StrainRate(VecAve(ibLoc)%Vec(uHatP)%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    VecAve(ibLoc)%Vec(vHatP)%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    VecAve(ibLoc)%Vec(wHatP)%c(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                    FU(ix-1:ix,iy:iy,iz:iz),              &
                    FV(ix:ix,iy-1:iy,iz:iz),              &
                    FW(ix:ix,iy:iy,iz-1:iz),              &
                    S,Sij)

!!        dthdz = GradCentrFull(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
!                         ,thC                                      &
!                         ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
!                         ,dz(iz-1),dz(iz),dz(iz+1))
         dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                          ,thC                                      &
                          ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                          ,FB,FT,VB,VCe,VT)

!         Local Length Scale Wall function
!         Decreasing of LengthScale near the wall aus Kleissl Dis.
!          LenScale=((1.0d0/(0.4d0*(zP(iz)+0.1d0))**One)+(1.0d0/(LenLoc**One)))**(-1.0d0/One)

    T11 = ABS(S)*Sij(1) ! Einheit: (m/sm)*(m/sm)=(1/s)*(1/s) normal noch mit LenScale**2 multipliziert
    T22 = ABS(S)*Sij(2)
    T33 = ABS(S)*Sij(3)
    T12 = ABS(S)*Sij(4)
    T13 = ABS(S)*Sij(5)
    T23 = ABS(S)*Sij(6)

    Lij(1) = uuHat(ix,iy,iz,1) - (uHat(ix,iy,iz,1)*uHat(ix,iy,iz,1))
    Lij(2) = vvHat(ix,iy,iz,1) - (vHat(ix,iy,iz,1)*vHat(ix,iy,iz,1))
    Lij(3) = wwHat(ix,iy,iz,1) - (wHat(ix,iy,iz,1)*wHat(ix,iy,iz,1))
    Lij(4) = uvHat(ix,iy,iz,1) - (uHat(ix,iy,iz,1)*vHat(ix,iy,iz,1))
    Lij(5) = uwHat(ix,iy,iz,1) - (uHat(ix,iy,iz,1)*wHat(ix,iy,iz,1))
    Lij(6) = vwHat(ix,iy,iz,1) - (vHat(ix,iy,iz,1)*wHat(ix,iy,iz,1))

    Mij(1) = Tau11Hat(ix,iy,iz,1) - (((2.0)**Two)*T11)
    Mij(2) = Tau22Hat(ix,iy,iz,1) - (((2.0)**Two)*T22)
    Mij(3) = Tau33Hat(ix,iy,iz,1) - (((2.0)**Two)*T33)
    Mij(4) = Tau12Hat(ix,iy,iz,1) - (((2.0)**Two)*T12)
    Mij(5) = Tau13Hat(ix,iy,iz,1) - (((2.0)**Two)*T13)
    Mij(6) = Tau23Hat(ix,iy,iz,1) - (((2.0)**Two)*T23)

!______________________________________________

     LM = Lij(1)*Mij(1) + Lij(2)*Mij(2) + Lij(3)*Mij(3) + Two*Lij(4)*Mij(4) + Two*Lij(5)*Mij(5) + Two*Lij(6)*Mij(6)
     MM = Mij(1)*Mij(1) + Mij(2)*Mij(2) + Mij(3)*Mij(3) + Two*Mij(4)*Mij(4) + Two*Mij(5)*Mij(5) + Two*Mij(6)*Mij(6)

     VecAve(ibLoc)%Vec(16)%c(ix,iy,iz,1) = LM
     VecAve(ibLoc)%Vec(17)%c(ix,iy,iz,1) = MM

!      Output as SpecialOutput   -Len_turb-
!     LenC(ix,iy,iz,1)=  SQRT(ABS((LM)/(2.0*MM+Eps)))*((LM)/(2.0*MM+Eps)/ABS(((LM)/(2.0*MM+Eps))+Eps))/LenLoc

      END DO
    END DO
  END DO
END SUBROUTINE DynSmagKoeffCompute

SUBROUTINE DFromDynSmagCompute !!!marcelk nach Tejada-Martinez or Froehlich Buch
!
!==================================================
!----  Diffusion Coefficient depending on dynamic Smagorinsky
!==================================================
!
!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: S
  REAL(RealKind) :: N2,Rich,Pr
  REAL(RealKind) :: FM,FH,PhiM
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC
  REAL(RealKind) :: LenLoc,Brunt,LenScale
  REAL(RealKind) :: rhoC
  REAL(RealKind) :: alpha,phi,Lennull,LM,MM,QN,NN
  REAL(RealKind) :: BETA,Beta1,Csmag,Csmag2,CsmagZae,CsmagNen
  REAL(RealKind) :: dxLoc,dyLoc,dzLoc,dxdz,dydz,dzdx,dxdy,delta
  REAL(RealKind),DIMENSION(1:6) :: Sij
  REAL(RealKind) :: x0, x1, epsilon 
  epsilon = 1.0d-10 

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe  = VolC(ix  ,iy  ,iz  )
        VB   = VolC(ix  ,iy  ,iz-1)
        VT   = VolC(ix  ,iy  ,iz+1)
        FB   = FW(ix  ,iy  ,iz-1)
        FT   = FW(ix  ,iy  ,iz  )
        rhoC = Rho(ix,iy,iz,1)
        thC  = th(ix,iy,iz,1)/(RhoC+Eps)
        LenLoc = (dx(ix)*dy(iy)*dz(iz))**(1.0d0/3.0d0)
!         Local Length Scale Wall function
!         Decreasing of LengthScale near the wall aus Kleissl Dis.
!          LenScale=((1.0d0/(0.4d0*(zP(iz)+0.1d0))**One)+(1.0d0/(LenLoc**One)))**(-1.0d0/One)

        LM = VecAve(ibLoc)%Vec(16)%c(ix,iy,iz,1)
        MM = VecAve(ibLoc)%Vec(17)%c(ix,iy,iz,1)

      Csmag = ((LM)/(2.0*MM+Eps))  !Toda_Cabrit_Balarac_Bose_Lee_Choi_Nicoud_PAPER

         IF(Csmag.LT.(-(5.0d-2**Two)*(LenLoc**Two))) THEN ! -0.05 < Cs < 1.0
         Csmag = -(5.0d-2**Two)*(LenLoc**Two)
         END IF
         IF(Csmag.GT.(LenLoc**Two)) THEN !2.0d-1
         Csmag = (LenLoc**Two) 
         ENDIF

!        Output as SpecialOutput   -Len_turb-
        LenC(ix,iy,iz,1)=  SQRT(ABS(Csmag))*(Csmag/ABS(Csmag+Eps))/LenLoc

    CALL StrainRateRho(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                    Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                    VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                    FU(ix-1:ix,iy:iy,iz:iz),              &
                    FV(ix:ix,iy-1:iy,iz:iz),              &
                    FW(ix:ix,iy:iy,iz-1:iz),              &
                    S,Sij)

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                         ,thC                                      &
                         ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                         ,FB,FT,VB,VCe,VT)

!       Richardson Number
        Brunt  = (Grav/(thC+Eps))*dthdz
        Rich  = (Grav/(thC+Eps))*dthdz/(S*S+1.d-3)

!    use this for depending LenScale on atmos. stability
!        IF (Rich<Zero) THEN
!!         Unstable Conditions
!          phi=(One-Four*alpha*Rich)**0.25d0
!          Lennull = LenLoc
!          LenScale=One/ (((((One-Rich)**0.25d0)*phi) / (0.4d0*(zP(iz)+0.1d0)))+((One) / (Lennull)))
!        ELSE IF (Rich >= RichC) THEN
!!         Very Stable Conditions
!          LenScale= LenLoc/LenLoc  ! eigentlich null
!        ELSE
!!         Neutral or slightly stable Conditions
!          Lennull = LenLoc * (One-alpha*Rich)
!          LenScale=One/ (((((One-Rich)**0.25d0)) / (0.4d0*(zP(iz)+0.1d0)))+((One) / (Lennull)))
!        END IF

!       Parametric Functions
        IF (Rich < Zero) THEN
!         Unstable Conditions
          FM = (One-16.0d0*Rich)**0.5d0
          FH = (One-40.0d0*Rich)**0.5d0/(Prn+Eps)
        ELSE IF (Rich >= RichC) THEN
!         Very Stable Conditions
          FM = 0.0d-5
          FH = 0.0d-5 
        ELSE
!         Stable Conditions
          FM = (One-Rich/RichC)**Four
          FH = (One-Rich/RichC)**Four*(One-1.2d0*Rich)/(Prn+Eps)
        END IF


!       Turbulent Kinetic Energy and Diffusion Coefficients
!         Prandtl Number
         Pr = MAX((FM/(FH+Eps)),3.3d-1)

!         Diffusion Coefficients Compute
!          DMomH(ix,iy,iz,1)=(Csmag)*(LenScale**Two)*S*FM
!          DMomV(ix,iy,iz,1)=(Csmag)*(LenScale**Two)*S*FM
!          DPotH(ix,iy,iz,1)=(Csmag)*(LenScale**Two)*S*FH
!          DPotV(ix,iy,iz,1)=(Csmag)*(LenScale**Two)*S*FH
          DMomH(ix,iy,iz,1)=(Csmag)*S*FM
          DMomV(ix,iy,iz,1)=(Csmag)*S*FM
          DPotH(ix,iy,iz,1)=(Csmag)*S*FH
          DPotV(ix,iy,iz,1)=(Csmag)*S*FH
      END DO
    END DO
  END DO

!-----------------------------------------------------------
END SUBROUTINE DFromDynSmagCompute

FUNCTION Filter(cL,cC,cR,VL,VC,VR)

  REAL(RealKind) :: Filter
  REAL(RealKind) :: cl,cC,cR,VL,VC,VR

!  Filter=(cL*VL+cC*VC+cR*VR)/(VL+VC+VR+Eps)
  Filter=(0.5*cL*VL+cC*VC+0.5*cR*VR)/(0.5*VL+VC+0.5*VR+Eps)

END FUNCTION Filter

FUNCTION Filter1(cL,cC,cR,VL,VC,VR)

  REAL(RealKind) :: Filter1
  REAL(RealKind) :: cl,cC,cR,VL,VC,VR

  Filter1=(cL*VL+cC*VC+cR*VR)/(VL+VC+VR+Eps)

END FUNCTION Filter1

FUNCTION Filter2(cLL,cL,cC,cR,cRR,VLL,VL,VC,VR,VRR)

  REAL(RealKind) :: Filter2
  REAL(RealKind) :: cl,cC,cR,VL,VC,VR
  REAL(RealKind) :: cLL,cRR,VLL,VRR

  Filter2=(0.5*cLL*VLL+cL*VL+cC*VC+cR*VR+0.5*cRR*VRR)/(0.5*VLL+VL+VC+VR+0.5*VRR+Eps)

END FUNCTION Filter2

SUBROUTINE FilterX(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
        c(i,j,k)=Filter(cLoc(i-1,j,k),cLoc(i,j,k),cLoc(i+1,j,k),VCeLoc(i-1,j,k),VCeLoc(i,j,k),VCeLoc(i+1,j,k))
      END DO
    END DO
  END DO
END SUBROUTINE FilterX

SUBROUTINE FilterY(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
        c(i,j,k)=Filter(cLoc(i,j-1,k),cLoc(i,j,k),cLoc(i,j+1,k),VCeLoc(i,j-1,k),VCeLoc(i,j,k),VCeLoc(i,j+1,k))
      END DO
    END DO
  END DO
END SUBROUTINE FilterY

SUBROUTINE FilterZ(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
!        VCeLoc(i,j,1)= VCeLoc(i,j,2)
!        cLoc(i,j,LBOUND(c,3)+1-1) = 0.0d0 ! marcelk Wind im Boden = 0
        c(i,j,k)=Filter(cLoc(i,j,k-1),cLoc(i,j,k),cLoc(i,j,k+1),VCeLoc(i,j,k-1),VCeLoc(i,j,k),VCeLoc(i,j,k+1))
      END DO
    END DO
  END DO
END SUBROUTINE FilterZ

SUBROUTINE FilterX2(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
        IF(i==(LBOUND(c,1)+1).or.i==(UBOUND(c,1)-1)) THEN
        c(i,j,k)=Filter1(cLoc(i-1,j,k),cLoc(i,j,k),cLoc(i+1,j,k),VCeLoc(i-1,j,k),VCeLoc(i,j,k),VCeLoc(i+1,j,k))
        ELSE
        c(i,j,k)=Filter2(cLoc(i-2,j,k),cLoc(i-1,j,k),cLoc(i,j,k),cLoc(i+1,j,k),&
        &cLoc(i+2,j,k),VCeLoc(i-2,j,k),VCeLoc(i-1,j,k),VCeLoc(i,j,k),VCeLoc(i+1,j,k),VCeLoc(i+2,j,k))
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE FilterX2

SUBROUTINE FilterY2(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
        IF(j==(LBOUND(c,2)+1).or.j==(UBOUND(c,2)-1)) THEN
        c(i,j,k)=Filter1(cLoc(i,j-1,k),cLoc(i,j,k),cLoc(i,j+1,k),VCeLoc(i,j-1,k),VCeLoc(i,j,k),VCeLoc(i,j+1,k))
        ELSE
        c(i,j,k)=Filter2(cLoc(i,j-2,k),cLoc(i,j-1,k),cLoc(i,j,k),cLoc(i,j+1,k),&
        &cLoc(i,j+2,k),VCeLoc(i,j-2,k),VCeLoc(i,j-1,k),VCeLoc(i,j,k),VCeLoc(i,j+1,k),VCeLoc(i,j+2,k))
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE FilterY2

SUBROUTINE FilterZ2(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
!        VCeLoc(i,j,1)= VCeLoc(i,j,2)
!        cLoc(i,j,LBOUND(c,3)+1-1) = 0.0d0 ! marcelk Wind im Boden = 0
        IF(k==(LBOUND(c,3)+1).or.k==(UBOUND(c,3)-1)) THEN
        c(i,j,k)=Filter1(cLoc(i,j,k-1),cLoc(i,j,k),cLoc(i,j,k+1),VCeLoc(i,j,k-1),VCeLoc(i,j,k),VCeLoc(i,j,k+1))
        ELSE
        c(i,j,k)=Filter2(cLoc(i,j,k-2),cLoc(i,j,k-1),cLoc(i,j,k),cLoc(i,j,k+1),&
        &cLoc(i,j,k+2),VCeLoc(i,j,k-2),VCeLoc(i,j,k-1),VCeLoc(i,j,k),VCeLoc(i,j,k+1),VCeLoc(i,j,k+2))
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE FilterZ2

SUBROUTINE FilterX3(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
        c(i,j,k)=Filter1(cLoc(i-1,j,k),cLoc(i,j,k),cLoc(i+1,j,k),VCeLoc(i-1,j,k),VCeLoc(i,j,k),VCeLoc(i+1,j,k))
      END DO
    END DO
  END DO
END SUBROUTINE FilterX3

SUBROUTINE FilterY3(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
        c(i,j,k)=Filter1(cLoc(i,j-1,k),cLoc(i,j,k),cLoc(i,j+1,k),VCeLoc(i,j-1,k),VCeLoc(i,j,k),VCeLoc(i,j+1,k))
      END DO
    END DO
  END DO
END SUBROUTINE FilterY3

SUBROUTINE FilterZ3(c,VCe)

  INTEGER i,j,k
  REAL(RealKind) :: c(:,:,:)
  REAL(RealKind) :: VCe(:,:,:)
  REAL(RealKind) :: VCeLoc(LBOUND(VCe,1):UBOUND(VCe,1),LBOUND(VCe,2):UBOUND(VCe,2),LBOUND(VCe,3):UBOUND(VCe,3))
  REAL(RealKind) :: cLoc(LBOUND(c,1):UBOUND(c,1),LBOUND(c,2):UBOUND(c,2),LBOUND(c,3):UBOUND(c,3))

  cLoc=c
  VCeLoc=VCe
  DO i=LBOUND(c,1)+1,UBOUND(c,1)-1
    DO j=LBOUND(c,2)+1,UBOUND(c,2)-1
      DO k=LBOUND(c,3)+1,UBOUND(c,3)-1
        c(i,j,k)=Filter1(cLoc(i,j,k-1),cLoc(i,j,k),cLoc(i,j,k+1),VCeLoc(i,j,k-1),VCeLoc(i,j,k),VCeLoc(i,j,k+1))
      END DO
    END DO
  END DO
END SUBROUTINE FilterZ3

SUBROUTINE DFromTkeLenCompute !!!FILAUS
!
!============================================================
!----  Diffusion Coefficient depending on TKE and Lengthscale
!============================================================
!
  IMPLICIT NONE

!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: S,O
  REAL(RealKind) :: N2,Rich,Pr
  REAL(RealKind) :: FM,FH,PhiM
  REAL(RealKind) :: FB,FT,VB,VT,VCe
  REAL(RealKind) :: dthdz,thC
  REAL(RealKind) :: tkeC,rhoC

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VCe  = VolC(ix  ,iy  ,iz  )
        VB   = VolC(ix  ,iy  ,iz-1)
        VT   = VolC(ix  ,iy  ,iz+1)
        FB   = FW(ix  ,iy  ,iz-1)
        FT   = FW(ix  ,iy  ,iz  )
        rhoC = Rho(ix,iy,iz,1)
        thC  = th(ix,iy,iz,1)/(RhoC+Eps)
        tkeC = MAX(tke(ix,iy,iz,1),Zero)
        
        CALL DeformVortC(uC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         vC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         wC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1),  &
                         Rho(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1,1), &
                         VolC(ix-1:ix+1,iy-1:iy+1,iz-1:iz+1),  &
                         FU(ix-1:ix,iy:iy,iz:iz),              &
                         FV(ix:ix,iy-1:iy,iz:iz),              &
                         FW(ix:ix,iy:iy,iz-1:iz),              &
                         S,O)

        dthdz = GradCentr(th(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps) &
                         ,thC                                      &
                         ,th(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps) &
                         ,FB,FT,VB,VCe,VT)
 
!       Richardson Number
        Rich  = (Grav/(thC+Eps))*dthdz/(S*S+1.d-3)
  
!       Parametric Functions
        IF (Rich < Zero) THEN
!         Unstable Conditions
          FM = (One-16.0d0*Rich)**0.5d0
          FH = (One-40.0d0*Rich)**0.5d0/(Prn+Eps)
        ELSE IF (Rich >= RichC) THEN
!         Very Stable Conditions
          FM = 0.0d0 
          FH = 0.0d0 
        ELSE  
!         Stable Conditions
          FM = (One-Rich/RichC)**Four
          FH = (One-Rich/RichC)**Four*(One-1.2d0*Rich)/(Prn+Eps)
        END IF

!       Turbulent Kinetic Energy and Diffusion Coefficients
        IF (Rich >= RichC) THEN
          PhiM             = 0.0d0
          DMom(ix,iy,iz,1) = RhoC*DiffMin
          DPot(ix,iy,iz,1) = RhoC*DiffMin
        ELSE
!         Prandtl Number
          Pr   = FM/( FH+Eps )
!         Stability Function for Momentum 
          PhiM = (Cs**(Four/Three)) * (Ceps**(One/Three)) * FM**(Two/Three) &
                 / (One-Rich/(Pr+Eps))**(One/Three)
!         Diffusion Coefficients Compute
          DMom(ix,iy,iz,1)=MAX(MIN(PhiM*LenC(ix,iy,iz,1)*SQRT(tkeC*rhoC),RhoC*DiffMax),RhoC*DiffMin)
          DPot(ix,iy,iz,1)=MAX(MIN(PhiM/(Pr+Eps)*LenC(ix,iy,iz,1)*SQRT(tkeC*rhoC),RhoC*DiffMax),RhoC*DiffMin)
!          DMom(ix,iy,iz,1) = MAX(MIN(Cmy0**(One/Four)*LenC(ix,iy,iz,1)*SQRT(ABS(tkeC*rhoC)) &
!                           ,RhoC*DiffMax),RhoC*DiffMin)
!         DPot(ix,iy,iz,1) = DMom(ix,iy,iz,1)                                                !!!FILAUS
!          DMom(ix,iy,iz,1)=MAX(MIN(PhiM*LenC*SQRT(ABS(tkeC*rhoC)) &
!                          ,RhoC*DiffMax),RhoC*DiffMin)
!          DPot(ix,iy,iz,1)=MAX(MIN(DMom(ix,iy,iz,1)*FH/(FM+Eps) &     
!                          ,RhoC*DiffMax),RhoC*DiffMin)
!          DMom(ix,iy,iz,1)=MAX(MIN((Cs*LenC)**Two*S*FM,RhoC*DiffMax),RhoC*DiffMin)
!          DPot(ix,iy,iz,1)=MAX(MIN((Cs*LenC)**Two*S*FM,RhoC*DiffMax),RhoC*DiffMin)

        END IF

      END DO
    END DO
  END DO

END SUBROUTINE DFromTkeLenCompute

SUBROUTINE DFromTkeHVLenCompute
!
!==================================================
!----  Diffusion Coefficient depending on Herzog
!==================================================
!
  IMPLICIT NONE

!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: tkeHC,tkeVC,rhoC

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        rhoC           = Rho(ix,iy,iz,1)
        tkeHC          = MAX(tkeH(ix,iy,iz,1),Zero)
        tkeVC          = MAX(tkeV(ix,iy,iz,1),Zero)
!        DH(ix,iy,iz,1) = c0*LenC(ix,iy,iz,1)*SQRT((Three/Two)*tkeHC*rhoC)   !!!Original
!        DV(ix,iy,iz,1) = c0*LenC(ix,iy,iz,1)*SQRT(Three*tkeVC*rhoC)
!        DH(ix,iy,iz,1) = Cmy0**0.25d0*LenC(ix,iy,iz,1)*SQRT((Three/Two)*tkeHC*rhoC)
!        DV(ix,iy,iz,1) = Cmy0**0.25d0*LenC(ix,iy,iz,1)*SQRT(Three*tkeVC*rhoC)
        DH(ix,iy,iz,1) = (Cmy0**0.25d0*SQRT(1+2*Cmy0**0.75d0/c2))*LenC(ix,iy,iz,1)*SQRT((Three/Two)*tkeHC*rhoC)  !!!c3=2*Cmy0**0.75
        DV(ix,iy,iz,1) = (Cmy0**0.25d0*SQRT(1+2*Cmy0**0.75d0/c2))*LenC(ix,iy,iz,1)*SQRT(Three*tkeVC*rhoC)
      END DO
    END DO
  END DO

!-----------------------------------------------------------
END SUBROUTINE DFromTkeHVLenCompute

SUBROUTINE LenScaleGeomCompute

  REAL(RealKind) :: dxLoc,dyLoc,dzLoc,dxdz,dydz,dzdx,dxdy,delta
  REAL(RealKind) :: aspec1,aspec2
  REAL(RealKind) :: FCorr

  dxLoc = MAXVAL(dx(ix0+1:ix1))
  dyLoc = MAXVAL(dy(iy0+1:iy1))
  dzLoc = MAXVAL(dz(iz0+1:iz1))
  IF (ix1-(ix0+1)==0) THEN
    dxLoc=dyLoc
  ELSE IF (iy1-(iy0+1)==0) THEN
    dyLoc=dxLoc
  END IF
  IF (GridType=='Globe') THEN 
    dxLoc=dxLoc*RadEarth
    dyLoc=dyLoc*RadEarth
  END IF  
  delta = dxLoc*dyLoc*dzLoc

! Aspect Ratio of the given Grid
  IF (dxLoc < dyLoc .AND. dyLoc < dzLoc) THEN
    aspec1 = dxLoc/dzLoc
    aspec2 = dyLoc/dzLoc
  ELSE IF (dxLoc < dzLoc .AND. dzLoc < dyLoc) THEN
    aspec1 = dxLoc/dyLoc
    aspec2 = dzLoc/dyLoc
  ELSE IF (dyLoc < dxLoc .AND. dxLoc < dzLoc) THEN
    aspec1 = dyLoc/dzLoc
    aspec2 = dxLoc/dzLoc
  ELSE IF (dyLoc < dzLoc .AND. dzLoc < dxLoc) THEN
    aspec1 = dyLoc/dxLoc
    aspec2 = dzLoc/dxLoc
  ELSE IF (dzLoc < dxLoc .AND. dxLoc < dyLoc) THEN
    aspec1 = dzLoc/dyLoc
    aspec2 = dxLoc/dyLoc
  ELSE IF (dzLoc < dyLoc .AND. dyLoc < dxLoc) THEN
    aspec1 = dzLoc/dxLoc
    aspec2 = dyLoc/dxLoc
  ELSE IF (dxLoc == dyLoc .AND. dyLoc < dzLoc) THEN
    aspec1 = dxLoc/dzLoc
    aspec2 = dyLoc/dzLoc
  ELSE IF (dxLoc == dzLoc .AND. dzLoc < dyLoc) THEN
    aspec1 = dxLoc/dyLoc
    aspec2 = dzLoc/dyLoc
  ELSE IF (dyLoc == dzLoc .AND. dzLoc < dxLoc) THEN
    aspec1 = dyLoc/dxLoc
    aspec2 = dzLoc/dxLoc
  ELSE IF (dxLoc < dyLoc .AND. dyLoc == dzLoc) THEN
    aspec1 = dxLoc/dzLoc
    aspec2 = dyLoc/dzLoc
  ELSE IF (dxLoc < dzLoc .AND. dzLoc == dyLoc) THEN
    aspec1 = dxLoc/dyLoc
    aspec2 = dzLoc/dyLoc
  ELSE IF (dyLoc < dxLoc .AND. dxLoc == dzLoc) THEN
    aspec1 = dyLoc/dzLoc
    aspec2 = dxLoc/dzLoc
  ELSE IF (dyLoc < dzLoc .AND. dzLoc == dxLoc) THEN
    aspec1 = dyLoc/dxLoc
    aspec2 = dzLoc/dxLoc
  ELSE IF (dzLoc < dxLoc .AND. dxLoc == dyLoc) THEN
    aspec1 = dzLoc/dyLoc
    aspec2 = dxLoc/dyLoc
  ELSE IF (dzLoc < dyLoc .AND. dyLoc == dxLoc) THEN
    aspec1 = dzLoc/dxLoc
    aspec2 = dyLoc/dxLoc
  ELSE IF (dxLoc == dyLoc .AND. dyLoc == dzLoc) THEN
    aspec1 = dxLoc/dzLoc
    aspec2 = dyLoc/dzLoc
  END IF

! Corrective Factor
    FCorr = COSH(((Four/27.0d0)*(log(aspec1)**Two-log(aspec1)*log(aspec2)+log(aspec2)**Two))**0.5d0)
! Geometric Length Scale
    LenScaleGeom           = delta**(One/Three)*FCorr
    Floor(ib)%LenScaleGeom = LenScaleGeom
END SUBROUTINE LenScaleGeomCompute

SUBROUTINE BoundaryZeroGrad(VectorChange)

  TYPE(Vector4Cell_T), TARGET :: VectorChange(:)

  INTEGER :: ic

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ic=1,SIZE(VectorChange(ibLoc)%Vec)
      c=>VectorChange(ibLoc)%Vec(ic)%c
      CALL BoundaryZeroGradCompute
    END DO
  END DO
END SUBROUTINE BoundaryZeroGrad

SUBROUTINE BoundaryZeroGradCompute

  INTEGER :: ix,iy,iz

  IF (TypeW=='ow') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          c(ix0,iy,iz,:)=c(ix0+1,iy,iz,:)
        END DO
      END DO
  END IF

  IF (TypeE=='oe') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          c(ix1+1,iy,iz,:)=c(ix1,iy,iz,:)
        END DO
      END DO
  END IF

  IF (TypeS=='os') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          c(ix,iy0,iz,:)=c(ix,iy0+1,iz,:)
        END DO
      END DO
  END IF

  IF (TypeN=='on') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          c(ix,iy1+1,iz,:)=c(ix,iy1,iz,:)
        END DO
      END DO
  END IF

  IF (TypeB=='ob') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          c(ix,iy,iz0,:)=c(ix,iy,iz0+1,:)
        END DO
      END DO
  END IF
  IF (TypeT=='ot') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          c(ix,iy,iz1+1,:)=c(ix,iy,iz1,:)
        END DO
      END DO
  END IF

END SUBROUTINE BoundaryZeroGradCompute

END MODULE DiffKoeff_Mod
