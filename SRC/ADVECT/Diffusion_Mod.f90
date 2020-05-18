MODULE Diffusion_Mod
  USE DataType_Mod
  USE Names_Mod
  USE DiffKoeff_Mod
  USE Operator_Mod

  IMPLICIT NONE

CONTAINS

SUBROUTINE DiffusionCompute

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: cL,cR
  REAL(RealKind) :: VL,VR
  REAL(RealKind) :: DL,DR
  REAL(RealKind) :: FF,DRhoF
  REAL(RealKind) :: Flux

  DO it=LBOUND(c,4),UBOUND(c,4)
! Fluxes
! x-Direction
  IF (TypeW(1:1)/='o') THEN
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,1)+Eps)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        Flux=  &
             DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix+1,iy,iz,it)=f(ix+1,iy,iz,it)-Flux
        f(ix  ,iy,iz,it)=f(ix  ,iy,iz,it)+Flux
      END DO
    END DO
  ELSE
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,1)+Eps)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DH(ix+1,iy,iz,1)
        FF=FU(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix+1,iy,iz,it)=f(ix+1,iy,iz,it)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,1)+Eps)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL) 
        f(ix+1,iy,iz,it)=f(ix+1,iy,iz,it)-Flux
        f(ix  ,iy,iz,it)=f(ix  ,iy,iz,it)+Flux
      END DO
    END DO
  END DO
  IF (TypeE(1:1)/='o') THEN
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,1)+Eps)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        Flux=  &
             DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix+1,iy,iz,it)=f(ix+1,iy,iz,it)-Flux
        f(ix  ,iy,iz,it)=f(ix  ,iy,iz,it)+Flux
      END DO
    END DO
  ELSE
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,1)+Eps)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DH(ix,iy,iz,1)
        FF=FU(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix  ,iy,iz,it)=f(ix  ,iy,iz,it)+Flux
      END DO
    END DO
  END IF

! y-Direction
  IF (TypeS(1:1)/='o') THEN
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy+1,iz,it)/(Rho(ix,iy+1,iz,1)+Eps)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        Flux=  &
             DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy+1,iz,it)=f(ix,iy+1,iz,it)-Flux
        f(ix,iy  ,iz,it)=f(ix,iy  ,iz,it)+Flux
      END DO
    END DO
  ELSE
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy+1,iz,it)/(Rho(ix,iy+1,iz,1)+Eps)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DH(ix,iy+1,iz,1)
        FF=FV(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy+1,iz,it)=f(ix,iy+1,iz,it)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy+1,iz,it)/(Rho(ix,iy+1,iz,1)+Eps)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy+1,iz,it)=f(ix,iy+1,iz,it)-Flux
        f(ix,iy  ,iz,it)=f(ix,iy  ,iz,it)+Flux
      END DO
   END DO
  END DO
  IF (TypeN(1:1)/='o') THEN
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy+1,iz,it)/(Rho(ix,iy+1,iz,1)+Eps)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        Flux=  &
             DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy+1,iz,it)=f(ix,iy+1,iz,it)-Flux
        f(ix,iy  ,iz,it)=f(ix,iy  ,iz,it)+Flux
      END DO
    END DO
  ELSE
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy+1,iz,it)/(Rho(ix,iy+1,iz,1)+Eps)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DH(ix,iy,iz,1)
        FF=FV(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy  ,iz,it)=f(ix,iy  ,iz,it)+Flux
      END DO
    END DO
  END IF
! z-Direction
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy,iz+1,it)/(Rho(ix,iy,iz+1,1)+Eps)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        Flux=  &
             DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy,iz+1,it)=f(ix,iy,iz+1,it)-Flux
        f(ix,iy,iz  ,it)=f(ix,iy,iz  ,it)+Flux
      END DO
    END DO
  ELSE
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy,iz+1,it)/(Rho(ix,iy,iz+1,1)+Eps)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DV(ix,iy,iz+1,1)
        FF=FW(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy,iz+1,it)=f(ix,iy,iz+1,it)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy,iz+1,it)/(Rho(ix,iy,iz+1,1)+Eps)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy,iz+1,it)=f(ix,iy,iz+1,it)-Flux
        f(ix,iy,iz  ,it)=f(ix,iy,iz  ,it)+Flux
      END DO
    END DO
  END DO
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy,iz+1,it)/(Rho(ix,iy,iz+1,1)+Eps)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        Flux=  &
             DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy,iz+1,it)=f(ix,iy,iz+1,it)-Flux
        f(ix,iy,iz  ,it)=f(ix,iy,iz  ,it)+Flux
      END DO
    END DO
  ELSE
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cL=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VL=VolC(ix,iy,iz)
        cR=c(ix,iy,iz+1,it)/(Rho(ix,iy,iz+1,1)+Eps)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DV(ix,iy,iz,1)
        FF=FW(ix,iy,iz)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        f(ix,iy,iz  ,it)=f(ix,iy,iz  ,it)+Flux
      END DO
    END DO
  END IF

  DO iz=iz0,iz1+1
    DO iy=iy0,iy1+1
      DO ix=ix0,ix1+1
        f(ix,iy,iz,it)=f(ix,iy,iz,it)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  END DO
END SUBROUTINE DiffusionCompute

SUBROUTINE DiffusionComputeU

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: cL,cR
  REAL(RealKind) :: VL,VR
  REAL(RealKind) :: DL,DR
  REAL(RealKind) :: FF,DRhoF
  REAL(RealKind) :: Flux,FluxCr
  REAL(RealKind) :: GradL,GradR

! Fluxes
! x-Direction
  IF (TypeW(1:1)/='o') THEN
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        FF=FU(ix,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)-Flux
        uRhsL(ix  ,iy,iz,1)=uRhsL(ix  ,iy,iz,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        uRhsR(ix+1,iy,iz,1)=uRhsR(ix+1,iy,iz,1)-Flux
        uRhsR(ix  ,iy,iz,1)=uRhsR(ix  ,iy,iz,1)+Flux
      END DO
    END DO
  ELSE
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DH(ix+1,iy,iz,1)
        FF=FU(ix,iy,iz)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)-Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        uRhsR(ix+1,iy,iz,1)=uRhsR(ix+1,iy,iz,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL) 
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)-Flux
        uRhsL(ix  ,iy,iz,1)=uRhsL(ix  ,iy,iz,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL) 
        uRhsR(ix+1,iy,iz,1)=uRhsR(ix+1,iy,iz,1)-Flux
        uRhsR(ix  ,iy,iz,1)=uRhsR(ix  ,iy,iz,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeE(1:1)/='o') THEN
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)-Flux
        uRhsL(ix  ,iy,iz,1)=uRhsL(ix  ,iy,iz,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        uRhsR(ix+1,iy,iz,1)=uRhsR(ix+1,iy,iz,1)-Flux
        uRhsR(ix  ,iy,iz,1)=uRhsR(ix  ,iy,iz,1)+Flux
      END DO
    END DO
  ELSE
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DH(ix,iy,iz,1)
        FF=FU(ix,iy,iz)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        uRhsL(ix  ,iy,iz,1)=uRhsL(ix  ,iy,iz,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        uRhsR(ix  ,iy,iz,1)=uRhsR(ix  ,iy,iz,1)+Flux
      END DO
    END DO
  END IF

! y-Direction
  IF (TypeS(1:1)/='o') THEN
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((vCL(ix-1,iy+1,iz,1)+vCR(ix-1,iy+1,iz,1))/(Rho(ix-1,iy+1,iz,1)+Eps) &
                            ,(vCL(ix  ,iy+1,iz,1)+vCR(ix  ,iy+1,iz,1))/(Rho(ix  ,iy+1,iz,1)+Eps) &
                            ,(vCL(ix+1,iy+1,iz,1)+vCR(ix+1,iy+1,iz,1))/(Rho(ix+1,iy+1,iz,1)+Eps) &
                             ,FU(ix-1,iy+1,iz),FU(ix,iy+1,iz)          & 
                             ,VolC(ix-1,iy+1,iz),VolC(ix,iy+1,iz),VolC(ix+1,iy+1,iz)) 
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy+1,iz,1)=uRhsL(ix,iy+1,iz,1)-Flux
        uRhsL(ix,iy  ,iz,1)=uRhsL(ix,iy  ,iz,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy+1,iz,1)=uRhsR(ix,iy+1,iz,1)-Flux
        uRhsR(ix,iy  ,iz,1)=uRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  ELSE
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DH(ix,iy+1,iz,1)
        FF=FV(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((vCL(ix-1,iy+1,iz,1)+vCR(ix-1,iy+1,iz,1))/(Rho(ix-1,iy+1,iz,1)+Eps) &
                            ,(vCL(ix  ,iy+1,iz,1)+vCR(ix  ,iy+1,iz,1))/(Rho(ix  ,iy+1,iz,1)+Eps) &
                            ,(vCL(ix+1,iy+1,iz,1)+vCR(ix+1,iy+1,iz,1))/(Rho(ix+1,iy+1,iz,1)+Eps) &
                             ,FU(ix-1,iy+1,iz),FU(ix,iy+1,iz)          &
                             ,VolC(ix-1,iy+1,iz),VolC(ix,iy+1,iz),VolC(ix+1,iy+1,iz))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy+1,iz,1)=uRhsL(ix,iy+1,iz,1)-Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy+1,iz,1)=uRhsR(ix,iy+1,iz,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        GradL=Zero
        GradL=Half*GradCentr((vCL(ix-1,iy,iz,1)+vCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                            ,(vCL(ix  ,iy,iz,1)+vCR(ix  ,iy,iz,1))/(Rho(ix  ,iy,iz,1)+Eps) &
                            ,(vCL(ix+1,iy,iz,1)+vCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps) &
                             ,FU(ix-1,iy,iz),FU(ix,iy,iz)          &
                             ,VolC(ix-1,iy,iz),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
        GradR=Half*GradCentr((vCL(ix-1,iy+1,iz,1)+vCR(ix-1,iy+1,iz,1))/(Rho(ix-1,iy+1,iz,1)+Eps) &
                            ,(vCL(ix  ,iy+1,iz,1)+vCR(ix  ,iy+1,iz,1))/(Rho(ix  ,iy+1,iz,1)+Eps) &
                            ,(vCL(ix+1,iy+1,iz,1)+vCR(ix+1,iy+1,iz,1))/(Rho(ix+1,iy+1,iz,1)+Eps) &
                             ,FU(ix-1,iy+1,iz),FU(ix,iy+1,iz)          &
                             ,VolC(ix-1,iy+1,iz),VolC(ix,iy+1,iz),VolC(ix+1,iy+1,iz))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy+1,iz,1)=uRhsL(ix,iy+1,iz,1)-Flux
        uRhsL(ix,iy  ,iz,1)=uRhsL(ix,iy  ,iz,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy+1,iz,1)=uRhsR(ix,iy+1,iz,1)-Flux
        uRhsR(ix,iy  ,iz,1)=uRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeN(1:1)/='o') THEN
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        GradL=Zero
        GradL=Half*GradCentr((vCL(ix-1,iy,iz,1)+vCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                            ,(vCL(ix  ,iy,iz,1)+vCR(ix  ,iy,iz,1))/(Rho(ix  ,iy,iz,1)+Eps) &
                            ,(vCL(ix+1,iy,iz,1)+vCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps) &
                             ,FU(ix-1,iy,iz),FU(ix,iy,iz)          &
                             ,VolC(ix-1,iy,iz),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy+1,iz,1)=uRhsL(ix,iy+1,iz,1)-Flux
        uRhsL(ix,iy  ,iz,1)=uRhsL(ix,iy  ,iz,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy+1,iz,1)=uRhsR(ix,iy+1,iz,1)-Flux
        uRhsR(ix,iy  ,iz,1)=uRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  ELSE
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DH(ix,iy,iz,1)
        FF=FV(ix,iy,iz)
        GradL=Zero
        GradL=Half*GradCentr((vCL(ix-1,iy,iz,1)+vCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                            ,(vCL(ix  ,iy,iz,1)+vCR(ix  ,iy,iz,1))/(Rho(ix  ,iy,iz,1)+Eps) &
                            ,(vCL(ix+1,iy,iz,1)+vCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps) &
                             ,FU(ix-1,iy,iz),FU(ix,iy,iz)          &
                             ,VolC(ix-1,iy,iz),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy  ,iz,1)=uRhsL(ix,iy  ,iz,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy  ,iz,1)=uRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  END IF
! z-Direction
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR) 
        FF=FW(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((wCL(ix-1,iy,iz+1,1)+wCR(ix-1,iy,iz+1,1))/(Rho(ix-1,iy,iz+1,1)+Eps) &
                            ,(wCL(ix  ,iy,iz+1,1)+wCR(ix  ,iy,iz+1,1))/(Rho(ix  ,iy,iz+1,1)+Eps) &
                            ,(wCL(ix+1,iy,iz+1,1)+wCR(ix+1,iy,iz+1,1))/(Rho(ix+1,iy,iz+1,1)+Eps) &
                             ,FU(ix-1,iy,iz+1),FU(ix,iy,iz+1)          &
                             ,VolC(ix-1,iy,iz+1),VolC(ix,iy,iz+1),VolC(ix+1,iy,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy,iz+1,1)=uRhsL(ix,iy,iz+1,1)-Flux
        uRhsL(ix,iy,iz  ,1)=uRhsL(ix,iy,iz  ,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy,iz+1,1)=uRhsR(ix,iy,iz+1,1)-Flux
        uRhsR(ix,iy,iz  ,1)=uRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  ELSE
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DV(ix,iy,iz+1,1)  
        FF=FW(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((wCL(ix-1,iy,iz+1,1)+wCR(ix-1,iy,iz+1,1))/(Rho(ix-1,iy,iz+1,1)+Eps) &
                            ,(wCL(ix  ,iy,iz+1,1)+wCR(ix  ,iy,iz+1,1))/(Rho(ix  ,iy,iz+1,1)+Eps) &
                            ,(wCL(ix+1,iy,iz+1,1)+wCR(ix+1,iy,iz+1,1))/(Rho(ix+1,iy,iz+1,1)+Eps) &
                             ,FU(ix-1,iy,iz+1),FU(ix,iy,iz+1)          &
                             ,VolC(ix-1,iy,iz+1),VolC(ix,iy,iz+1),VolC(ix+1,iy,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy,iz+1,1)=uRhsL(ix,iy,iz+1,1)-Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy,iz+1,1)=uRhsR(ix,iy,iz+1,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR) 
        FF=FW(ix,iy,iz)
        GradL=Half*GradCentr((wCL(ix-1,iy,iz,1)+wCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                            ,(wCL(ix  ,iy,iz,1)+wCR(ix  ,iy,iz,1))/(Rho(ix  ,iy,iz,1)+Eps) &
                            ,(wCL(ix+1,iy,iz,1)+wCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps) &
                             ,FU(ix-1,iy,iz),FU(ix,iy,iz)          &
                             ,VolC(ix-1,iy,iz),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
        GradR=Half*GradCentr((wCL(ix-1,iy,iz+1,1)+wCR(ix-1,iy,iz+1,1))/(Rho(ix-1,iy,iz+1,1)+Eps) &
                            ,(wCL(ix  ,iy,iz+1,1)+wCR(ix  ,iy,iz+1,1))/(Rho(ix  ,iy,iz+1,1)+Eps) &
                            ,(wCL(ix+1,iy,iz+1,1)+wCR(ix+1,iy,iz+1,1))/(Rho(ix+1,iy,iz+1,1)+Eps) &
                             ,FU(ix-1,iy,iz+1),FU(ix,iy,iz+1)          &
                             ,VolC(ix-1,iy,iz+1),VolC(ix,iy,iz+1),VolC(ix+1,iy,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy,iz+1,1)=uRhsL(ix,iy,iz+1,1)-Flux
        uRhsL(ix,iy,iz  ,1)=uRhsL(ix,iy,iz  ,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy,iz+1,1)=uRhsR(ix,iy,iz+1,1)-Flux
        uRhsR(ix,iy,iz  ,1)=uRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)  
        FF=FW(ix,iy,iz)
        GradL=Half*GradCentr((wCL(ix-1,iy,iz,1)+wCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                            ,(wCL(ix  ,iy,iz,1)+wCR(ix  ,iy,iz,1))/(Rho(ix  ,iy,iz,1)+Eps) &
                            ,(wCL(ix+1,iy,iz,1)+wCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps) &
                             ,FU(ix-1,iy,iz),FU(ix,iy,iz)          &
                             ,VolC(ix-1,iy,iz),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy,iz+1,1)=uRhsL(ix,iy,iz+1,1)-Flux
        uRhsL(ix,iy,iz  ,1)=uRhsL(ix,iy,iz  ,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy,iz+1,1)=uRhsR(ix,iy,iz+1,1)-Flux
        uRhsR(ix,iy,iz  ,1)=uRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  ELSE
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DV(ix,iy,iz,1) 
        FF=FW(ix,iy,iz)
        GradL=Half*GradCentr((wCL(ix-1,iy,iz,1)+wCR(ix-1,iy,iz,1))/(Rho(ix-1,iy,iz,1)+Eps) &
                            ,(wCL(ix  ,iy,iz,1)+wCR(ix  ,iy,iz,1))/(Rho(ix  ,iy,iz,1)+Eps) &
                            ,(wCL(ix+1,iy,iz,1)+wCR(ix+1,iy,iz,1))/(Rho(ix+1,iy,iz,1)+Eps) &
                             ,FU(ix-1,iy,iz),FU(ix,iy,iz)          &
                             ,VolC(ix-1,iy,iz),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=ucL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsL(ix,iy,iz  ,1)=uRhsL(ix,iy,iz  ,1)+Flux
        cL=ucR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=ucR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        uRhsR(ix,iy,iz  ,1)=uRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  END IF

  DO iz=iz0,iz1+1
    DO iy=iy0,iy1+1
      DO ix=ix0,ix1+1
        uRhsL(ix,iy,iz,1)=uRhsL(ix,iy,iz,1)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE DiffusionComputeU

SUBROUTINE DiffusionComputeV

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: cL,cR
  REAL(RealKind) :: VL,VR
  REAL(RealKind) :: DL,DR
  REAL(RealKind) :: FF,DRhoF
  REAL(RealKind) :: Flux,FluxCr
  REAL(RealKind) :: GradL,GradR

! Fluxes
! x-Direction
  IF (TypeW(1:1)/='o') THEN
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((uCL(ix+1,iy-1,iz,1)+uCR(ix+1,iy-1,iz,1))/(Rho(ix+1,iy-1,iz,1)+Eps) &
                            ,(uCL(ix+1,iy  ,iz,1)+uCR(ix+1,iy  ,iz,1))/(Rho(ix+1,iy  ,iz,1)+Eps) &
                            ,(uCL(ix+1,iy+1,iz,1)+uCR(ix+1,iy+1,iz,1))/(Rho(ix+1,iy+1,iz,1)+Eps) &
                             ,FV(ix+1,iy-1,iz),FV(ix+1,iy,iz)          &
                             ,VolC(ix+1,iy-1,iz),VolC(ix+1,iy,iz),VolC(ix+1,iy+1,iz))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix+1,iy,iz,1)=vRhsL(ix+1,iy,iz,1)-Flux
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy  ,iz,1)+Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix+1,iy,iz,1)=vRhsR(ix+1,iy,iz,1)-Flux
        vRhsR(ix,iy  ,iz,1)=vRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  ELSE
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DH(ix+1,iy,iz,1)
        FF=FU(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((uCL(ix+1,iy-1,iz,1)+uCR(ix+1,iy-1,iz,1))/(Rho(ix+1,iy-1,iz,1)+Eps) &
                            ,(uCL(ix+1,iy  ,iz,1)+uCR(ix+1,iy  ,iz,1))/(Rho(ix+1,iy  ,iz,1)+Eps) &
                            ,(uCL(ix+1,iy+1,iz,1)+uCR(ix+1,iy+1,iz,1))/(Rho(ix+1,iy+1,iz,1)+Eps) &
                             ,FV(ix+1,iy-1,iz),FV(ix+1,iy,iz)          &
                             ,VolC(ix+1,iy-1,iz),VolC(ix+1,iy,iz),VolC(ix+1,iy+1,iz))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix+1,iy,iz,1)=vRhsL(ix+1,iy,iz,1)-Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix+1,iy,iz,1)=vRhsR(ix+1,iy,iz,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        GradL=Half*GradCentr((uCL(ix,iy-1,iz,1)+uCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                            ,(uCL(ix,iy  ,iz,1)+uCR(ix,iy  ,iz,1))/(Rho(ix,iy  ,iz,1)+Eps) &
                            ,(uCL(ix,iy+1,iz,1)+uCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps) &
                             ,FV(ix,iy-1,iz),FV(ix,iy,iz)          &
                             ,VolC(ix,iy-1,iz),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
        GradR=Half*GradCentr((uCL(ix+1,iy-1,iz,1)+uCR(ix+1,iy-1,iz,1))/(Rho(ix+1,iy-1,iz,1)+Eps) &
                            ,(uCL(ix+1,iy  ,iz,1)+uCR(ix+1,iy  ,iz,1))/(Rho(ix+1,iy  ,iz,1)+Eps) &
                            ,(uCL(ix+1,iy+1,iz,1)+uCR(ix+1,iy+1,iz,1))/(Rho(ix+1,iy+1,iz,1)+Eps) &
                             ,FV(ix+1,iy-1,iz),FV(ix+1,iy,iz)          &
                             ,VolC(ix+1,iy-1,iz),VolC(ix+1,iy,iz),VolC(ix+1,iy+1,iz))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix+1,iy,iz,1)=vRhsL(ix+1,iy,iz,1)-Flux
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy  ,iz,1)+Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix+1,iy,iz,1)=vRhsR(ix+1,iy,iz,1)-Flux
        vRhsR(ix,iy  ,iz,1)=vRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeE(1:1)/='o') THEN
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        GradL=Half*GradCentr((uCL(ix,iy-1,iz,1)+uCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                            ,(uCL(ix,iy  ,iz,1)+uCR(ix,iy  ,iz,1))/(Rho(ix,iy  ,iz,1)+Eps) &
                            ,(uCL(ix,iy+1,iz,1)+uCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps) &
                             ,FV(ix,iy-1,iz),FV(ix,iy,iz)          &
                             ,VolC(ix,iy-1,iz),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix+1,iy,iz,1)=vRhsL(ix+1,iy,iz,1)-Flux
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy  ,iz,1)+Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix+1,iy,iz,1)=vRhsR(ix+1,iy,iz,1)-Flux
        vRhsR(ix,iy  ,iz,1)=vRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  ELSE
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DH(ix,iy,iz,1)
        FF=FU(ix,iy,iz)
        GradL=Half*GradCentr((uCL(ix,iy-1,iz,1)+uCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                            ,(uCL(ix,iy  ,iz,1)+uCR(ix,iy  ,iz,1))/(Rho(ix,iy  ,iz,1)+Eps) &
                            ,(uCL(ix,iy+1,iz,1)+uCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps) &
                             ,FV(ix,iy-1,iz),FV(ix,iy,iz)          &
                             ,VolC(ix,iy-1,iz),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy  ,iz,1)+Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix,iy  ,iz,1)=vRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  END IF

! y-Direction
  IF (TypeS(1:1)/='o') THEN
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)-Flux
        vRhsL(ix,iy  ,iz,1)=vRhsL(ix,iy  ,iz,1)+Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsR(ix,iy+1,iz,1)=vRhsR(ix,iy+1,iz,1)-Flux
        vRhsR(ix,iy  ,iz,1)=vRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  ELSE
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DH(ix,iy+1,iz,1)
        FF=FV(ix,iy,iz)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)-Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsR(ix,iy+1,iz,1)=vRhsR(ix,iy+1,iz,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)-Flux
        vRhsL(ix,iy  ,iz,1)=vRhsL(ix,iy  ,iz,1)+Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsR(ix,iy+1,iz,1)=vRhsR(ix,iy+1,iz,1)-Flux
        vRhsR(ix,iy  ,iz,1)=vRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeN(1:1)/='o') THEN
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)-Flux
        vRhsL(ix,iy  ,iz,1)=vRhsL(ix,iy  ,iz,1)+Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsR(ix,iy+1,iz,1)=vRhsR(ix,iy+1,iz,1)-Flux
        vRhsR(ix,iy  ,iz,1)=vRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  ELSE
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DH(ix,iy,iz,1)
        FF=FV(ix,iy,iz)
        cL=vcL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsL(ix,iy  ,iz,1)=vRhsL(ix,iy  ,iz,1)+Flux
        cL=vcR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vcR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        vRhsR(ix,iy  ,iz,1)=vRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  END IF
! z-Direction
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        FF=FW(ix,iy,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR) 
        GradL=Zero
        GradR=Half*GradCentr((wCL(ix,iy-1,iz+1,1)+wCR(ix,iy-1,iz+1,1))/(Rho(ix,iy-1,iz+1,1)+Eps) &
                            ,(wCL(ix,iy  ,iz+1,1)+wCR(ix,iy  ,iz+1,1))/(Rho(ix,iy  ,iz+1,1)+Eps) &
                            ,(wCL(ix,iy+1,iz+1,1)+wCR(ix,iy+1,iz+1,1))/(Rho(ix,iy+1,iz+1,1)+Eps) &
                            ,FV(ix,iy-1,iz+1),FV(ix,iy,iz+1)          & 
                            ,VolC(ix,iy-1,iz+1),VolC(ix,iy,iz+1),VolC(ix,iy+1,iz+1)) 
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr) 
        vRhsL(ix,iy,iz+1,1)=vRhsL(ix,iy,iz+1,1)-Flux
        vRhsL(ix,iy,iz  ,1)=vRhsL(ix,iy,iz  ,1)+Flux
        cL=vCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr) 
        vRhsR(ix,iy,iz+1,1)=vRhsR(ix,iy,iz+1,1)-Flux
        vRhsR(ix,iy,iz  ,1)=vRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  ELSE
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        FF=FW(ix,iy,iz)
        DRhoF=DV(ix,iy,iz+1,1) 
        GradL=Zero
        GradR=Half*GradCentr((wCL(ix,iy-1,iz+1,1)+wCR(ix,iy-1,iz+1,1))/(Rho(ix,iy-1,iz+1,1)+Eps) &
                            ,(wCL(ix,iy  ,iz+1,1)+wCR(ix,iy  ,iz+1,1))/(Rho(ix,iy  ,iz+1,1)+Eps) &
                            ,(wCL(ix,iy+1,iz+1,1)+wCR(ix,iy+1,iz+1,1))/(Rho(ix,iy+1,iz+1,1)+Eps) &
                            ,FV(ix,iy-1,iz+1),FV(ix,iy,iz+1)          &
                            ,VolC(ix,iy-1,iz+1),VolC(ix,iy,iz+1),VolC(ix,iy+1,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix,iy,iz+1,1)=vRhsL(ix,iy,iz+1,1)-Flux
        cL=vCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix,iy,iz+1,1)=vRhsR(ix,iy,iz+1,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        FF=FW(ix,iy,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)  
        GradL=Half*GradCentr((wCL(ix,iy-1,iz,1)+wCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                            ,(wCL(ix,iy  ,iz,1)+wCR(ix,iy  ,iz,1))/(Rho(ix,iy  ,iz,1)+Eps) &
                            ,(wCL(ix,iy+1,iz,1)+wCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps) &
                            ,FV(ix,iy-1,iz),FV(ix,iy,iz)          &
                            ,VolC(ix,iy-1,iz),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
        GradR=Half*GradCentr((wCL(ix,iy-1,iz+1,1)+wCR(ix,iy-1,iz+1,1))/(Rho(ix,iy-1,iz+1,1)+Eps) &
                            ,(wCL(ix,iy  ,iz+1,1)+wCR(ix,iy  ,iz+1,1))/(Rho(ix,iy  ,iz+1,1)+Eps) &
                            ,(wCL(ix,iy+1,iz+1,1)+wCR(ix,iy+1,iz+1,1))/(Rho(ix,iy+1,iz+1,1)+Eps) &
                            ,FV(ix,iy-1,iz+1),FV(ix,iy,iz+1)          &
                            ,VolC(ix,iy-1,iz+1),VolC(ix,iy,iz+1),VolC(ix,iy+1,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix,iy,iz+1,1)=vRhsL(ix,iy,iz+1,1)-Flux
        vRhsL(ix,iy,iz  ,1)=vRhsL(ix,iy,iz  ,1)+Flux
        cL=vCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix,iy,iz+1,1)=vRhsR(ix,iy,iz+1,1)-Flux
        vRhsR(ix,iy,iz  ,1)=vRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        FF=FW(ix,iy,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR) 
        GradL=Half*GradCentr((wCL(ix,iy-1,iz,1)+wCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                            ,(wCL(ix,iy  ,iz,1)+wCR(ix,iy  ,iz,1))/(Rho(ix,iy  ,iz,1)+Eps) &
                            ,(wCL(ix,iy+1,iz,1)+wCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps) &
                            ,FV(ix,iy-1,iz),FV(ix,iy,iz)          &
                            ,VolC(ix,iy-1,iz),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix,iy,iz+1,1)=vRhsL(ix,iy,iz+1,1)-Flux
        vRhsL(ix,iy,iz  ,1)=vRhsL(ix,iy,iz  ,1)+Flux
        cL=vCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix,iy,iz+1,1)=vRhsR(ix,iy,iz+1,1)-Flux
        vRhsR(ix,iy,iz  ,1)=vRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  ELSE
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        FF=FW(ix,iy,iz)
        DRhoF=DV(ix,iy,iz,1) 
        GradL=Half*GradCentr((wCL(ix,iy-1,iz,1)+wCR(ix,iy-1,iz,1))/(Rho(ix,iy-1,iz,1)+Eps) &
                            ,(wCL(ix,iy  ,iz,1)+wCR(ix,iy  ,iz,1))/(Rho(ix,iy  ,iz,1)+Eps) &
                            ,(wCL(ix,iy+1,iz,1)+wCR(ix,iy+1,iz,1))/(Rho(ix,iy+1,iz,1)+Eps) &
                            ,FV(ix,iy-1,iz),FV(ix,iy,iz)          &
                            ,VolC(ix,iy-1,iz),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=vCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsL(ix,iy,iz  ,1)=vRhsL(ix,iy,iz  ,1)+Flux
        cL=vCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=vCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        vRhsR(ix,iy,iz  ,1)=vRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  END IF

  DO iz=iz0,iz1+1
    DO iy=iy0,iy1+1
      DO ix=ix0,ix1+1
        vRhsL(ix,iy,iz,1)=vRhsL(ix,iy,iz,1)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE DiffusionComputeV

SUBROUTINE DiffusionComputeW

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: cL,cR
  REAL(RealKind) :: VL,VR
  REAL(RealKind) :: DL,DR
  REAL(RealKind) :: FF,DRhoF
  REAL(RealKind) :: Flux,FluxCr
  REAL(RealKind) :: GradL,GradR

! Fluxes
! x-Direction
  IF (TypeW(1:1)/='o') THEN
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((uCL(ix+1,iy,iz-1,1)+uCR(ix+1,iy,iz-1,1))/(Rho(ix+1,iy,iz-1,1)+Eps) &
                            ,(uCL(ix+1,iy,iz  ,1)+uCR(ix+1,iy,iz  ,1))/(Rho(ix+1,iy,iz  ,1)+Eps) &
                            ,(uCL(ix+1,iy,iz+1,1)+uCR(ix+1,iy,iz+1,1))/(Rho(ix+1,iy,iz+1,1)+Eps) &
                            ,FW(ix+1,iy,iz-1),FW(ix+1,iy,iz)          & 
                            ,VolC(ix+1,iy,iz-1),VolC(ix+1,iy,iz),VolC(ix+1,iy,iz+1)) 
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr) 
        wRhsL(ix+1,iy,iz,1)=wRhsL(ix+1,iy,iz,1)-Flux
        wRhsL(ix  ,iy,iz,1)=wRhsL(ix  ,iy,iz,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr) 
        wRhsR(ix+1,iy,iz,1)=wRhsR(ix+1,iy,iz,1)-Flux
        wRhsR(ix  ,iy,iz,1)=wRhsR(ix  ,iy,iz,1)+Flux
      END DO
    END DO
  ELSE
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DV(ix+1,iy,iz,1)
        FF=FU(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((uCL(ix+1,iy,iz-1,1)+uCR(ix+1,iy,iz-1,1))/(Rho(ix+1,iy,iz-1,1)+Eps) &
                            ,(uCL(ix+1,iy,iz  ,1)+uCR(ix+1,iy,iz  ,1))/(Rho(ix+1,iy,iz  ,1)+Eps) &
                            ,(uCL(ix+1,iy,iz+1,1)+uCR(ix+1,iy,iz+1,1))/(Rho(ix+1,iy,iz+1,1)+Eps) &
                            ,FW(ix+1,iy,iz-1),FW(ix+1,iy,iz)          &
                            ,VolC(ix+1,iy,iz-1),VolC(ix+1,iy,iz),VolC(ix+1,iy,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix+1,iy,iz,1)=wRhsL(ix+1,iy,iz,1)-Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix+1,iy,iz,1)=wRhsR(ix+1,iy,iz,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        GradL=Half*GradCentr((uCL(ix,iy,iz-1,1)+uCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                            ,(uCL(ix,iy,iz  ,1)+uCR(ix,iy,iz  ,1))/(Rho(ix,iy,iz  ,1)+Eps) &
                            ,(uCL(ix,iy,iz+1,1)+uCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps) &
                            ,FW(ix,iy,iz-1),FW(ix,iy,iz)          &
                            ,VolC(ix,iy,iz-1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
        GradR=Half*GradCentr((uCL(ix+1,iy,iz-1,1)+uCR(ix+1,iy,iz-1,1))/(Rho(ix+1,iy,iz-1,1)+Eps) &
                            ,(uCL(ix+1,iy,iz  ,1)+uCR(ix+1,iy,iz  ,1))/(Rho(ix+1,iy,iz  ,1)+Eps) &
                            ,(uCL(ix+1,iy,iz+1,1)+uCR(ix+1,iy,iz+1,1))/(Rho(ix+1,iy,iz+1,1)+Eps) &
                            ,FW(ix+1,iy,iz-1),FW(ix+1,iy,iz)          &
                            ,VolC(ix+1,iy,iz-1),VolC(ix+1,iy,iz),VolC(ix+1,iy,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix+1,iy,iz,1)=wRhsL(ix+1,iy,iz,1)-Flux
        wRhsL(ix  ,iy,iz,1)=wRhsL(ix  ,iy,iz,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix+1,iy,iz,1)=wRhsR(ix+1,iy,iz,1)-Flux
        wRhsR(ix  ,iy,iz,1)=wRhsR(ix  ,iy,iz,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeE(1:1)/='o') THEN
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        GradL=Half*GradCentr((uCL(ix,iy,iz-1,1)+uCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                            ,(uCL(ix,iy,iz  ,1)+uCR(ix,iy,iz  ,1))/(Rho(ix,iy,iz  ,1)+Eps) &
                            ,(uCL(ix,iy,iz+1,1)+uCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps) &
                            ,FW(ix,iy,iz-1),FW(ix,iy,iz)          &
                            ,VolC(ix,iy,iz-1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix+1,iy,iz,1)=wRhsL(ix+1,iy,iz,1)-Flux
        wRhsL(ix  ,iy,iz,1)=wRhsL(ix  ,iy,iz,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix+1,iy,iz,1)=wRhsR(ix+1,iy,iz,1)-Flux
        wRhsR(ix  ,iy,iz,1)=wRhsR(ix  ,iy,iz,1)+Flux
      END DO
    END DO
  ELSE
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DV(ix,iy,iz,1)
        FF=FU(ix,iy,iz)
        GradL=Half*GradCentr((uCL(ix,iy,iz-1,1)+uCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                            ,(uCL(ix,iy,iz  ,1)+uCR(ix,iy,iz  ,1))/(Rho(ix,iy,iz  ,1)+Eps) &
                            ,(uCL(ix,iy,iz+1,1)+uCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps) &
                            ,FW(ix,iy,iz-1),FW(ix,iy,iz)          &
                            ,VolC(ix,iy,iz-1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix  ,iy,iz,1)=wRhsL(ix  ,iy,iz,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix+1,iy,iz,1)/(Rho(ix+1,iy,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix  ,iy,iz,1)=wRhsR(ix  ,iy,iz,1)+Flux
      END DO
    END DO
  END IF

! y-Direction
  IF (TypeS(1:1)/='o') THEN
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((vCL(ix,iy+1,iz-1,1)+vCR(ix,iy+1,iz-1,1))/(Rho(ix,iy+1,iz-1,1)+Eps) &
                            ,(vCL(ix,iy+1,iz  ,1)+vCR(ix,iy+1,iz  ,1))/(Rho(ix,iy+1,iz  ,1)+Eps) &
                            ,(vCL(ix,iy+1,iz+1,1)+vCR(ix,iy+1,iz+1,1))/(Rho(ix,iy+1,iz+1,1)+Eps) &
                            ,FW(ix,iy+1,iz-1),FW(ix,iy+1,iz)          & 
                            ,VolC(ix,iy+1,iz-1),VolC(ix,iy+1,iz),VolC(ix,iy+1,iz+1)) 
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix,iy+1,iz,1)=wRhsL(ix,iy+1,iz,1)-Flux
        wRhsL(ix,iy  ,iz,1)=wRhsL(ix,iy  ,iz,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix,iy+1,iz,1)=wRhsR(ix,iy+1,iz,1)-Flux
        wRhsR(ix,iy  ,iz,1)=wRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  ELSE
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DV(ix,iy+1,iz,1)
        FF=FV(ix,iy,iz)
        GradL=Zero
        GradR=Half*GradCentr((vCL(ix,iy+1,iz-1,1)+vCR(ix,iy+1,iz-1,1))/(Rho(ix,iy+1,iz-1,1)+Eps) &
                            ,(vCL(ix,iy+1,iz  ,1)+vCR(ix,iy+1,iz  ,1))/(Rho(ix,iy+1,iz  ,1)+Eps) &
                            ,(vCL(ix,iy+1,iz+1,1)+vCR(ix,iy+1,iz+1,1))/(Rho(ix,iy+1,iz+1,1)+Eps) &
                            ,FW(ix,iy+1,iz-1),FW(ix,iy+1,iz)          &
                            ,VolC(ix,iy+1,iz-1),VolC(ix,iy+1,iz),VolC(ix,iy+1,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix,iy+1,iz,1)=wRhsL(ix,iy+1,iz,1)-Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix,iy+1,iz,1)=wRhsR(ix,iy+1,iz,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        GradL=Half*GradCentr((vCL(ix,iy,iz-1,1)+vCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                            ,(vCL(ix,iy,iz  ,1)+vCR(ix,iy,iz  ,1))/(Rho(ix,iy,iz  ,1)+Eps) &
                            ,(vCL(ix,iy,iz+1,1)+vCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps) &
                            ,FW(ix,iy,iz-1),FW(ix,iy,iz)          &
                            ,VolC(ix,iy,iz-1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
        GradR=Half*GradCentr((vCL(ix,iy+1,iz-1,1)+vCR(ix,iy+1,iz-1,1))/(Rho(ix,iy+1,iz-1,1)+Eps) &
                            ,(vCL(ix,iy+1,iz  ,1)+vCR(ix,iy+1,iz  ,1))/(Rho(ix,iy+1,iz  ,1)+Eps) &
                            ,(vCL(ix,iy+1,iz+1,1)+vCR(ix,iy+1,iz+1,1))/(Rho(ix,iy+1,iz+1,1)+Eps) &
                            ,FW(ix,iy+1,iz-1),FW(ix,iy+1,iz)          &
                            ,VolC(ix,iy+1,iz-1),VolC(ix,iy+1,iz),VolC(ix,iy+1,iz+1))
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix,iy+1,iz,1)=wRhsL(ix,iy+1,iz,1)-Flux
        wRhsL(ix,iy  ,iz,1)=wRhsL(ix,iy  ,iz,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(Two*FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix,iy+1,iz,1)=wRhsR(ix,iy+1,iz,1)-Flux
        wRhsR(ix,iy  ,iz,1)=wRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeN(1:1)/='o') THEN
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        GradL=Half*GradCentr((vCL(ix,iy,iz-1,1)+vCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                            ,(vCL(ix,iy,iz  ,1)+vCR(ix,iy,iz  ,1))/(Rho(ix,iy,iz  ,1)+Eps) &
                            ,(vCL(ix,iy,iz+1,1)+vCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps) &
                            ,FW(ix,iy,iz-1),FW(ix,iy,iz)          &
                            ,VolC(ix,iy,iz-1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix,iy+1,iz,1)=wRhsL(ix,iy+1,iz,1)-Flux
        wRhsL(ix,iy  ,iz,1)=wRhsL(ix,iy  ,iz,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix,iy+1,iz,1)=wRhsR(ix,iy+1,iz,1)-Flux
        wRhsR(ix,iy  ,iz,1)=wRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  ELSE
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DV(ix,iy,iz,1)
        FF=FV(ix,iy,iz)
        GradL=Half*GradCentr((vCL(ix,iy,iz-1,1)+vCR(ix,iy,iz-1,1))/(Rho(ix,iy,iz-1,1)+Eps) &
                            ,(vCL(ix,iy,iz  ,1)+vCR(ix,iy,iz  ,1))/(Rho(ix,iy,iz  ,1)+Eps) &
                            ,(vCL(ix,iy,iz+1,1)+vCR(ix,iy,iz+1,1))/(Rho(ix,iy,iz+1,1)+Eps) &
                            ,FW(ix,iy,iz-1),FW(ix,iy,iz)          &
                            ,VolC(ix,iy,iz-1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
        GradR=Zero
        FluxCr=IntCellToFace(GradL,GradR,VL,VR)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsL(ix,iy  ,iz,1)=wRhsL(ix,iy  ,iz,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy+1,iz,1)/(Rho(ix,iy+1,iz,1)+Eps)
        Flux=Two*DRhoF*FF*(FF/(VL+VR+Eps)*(cR-cL)+FluxCr)
        wRhsR(ix,iy  ,iz,1)=wRhsR(ix,iy  ,iz,1)+Flux
      END DO
    END DO
  END IF
! z-Direction
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsL(ix,iy,iz+1,1)=wRhsL(ix,iy,iz+1,1)-Flux
        wRhsL(ix,iy,iz  ,1)=wRhsL(ix,iy,iz  ,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsR(ix,iy,iz+1,1)=wRhsR(ix,iy,iz+1,1)-Flux
        wRhsR(ix,iy,iz  ,1)=wRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  ELSE
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
       VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DV(ix,iy,iz+1,1)
        FF=FW(ix,iy,iz)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsL(ix,iy,iz+1,1)=wRhsL(ix,iy,iz+1,1)-Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsR(ix,iy,iz+1,1)=wRhsR(ix,iy,iz+1,1)-Flux
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
       VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsL(ix,iy,iz+1,1)=wRhsL(ix,iy,iz+1,1)-Flux
        wRhsL(ix,iy,iz  ,1)=wRhsL(ix,iy,iz  ,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsR(ix,iy,iz+1,1)=wRhsR(ix,iy,iz+1,1)-Flux
        wRhsR(ix,iy,iz  ,1)=wRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  END DO
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsL(ix,iy,iz+1,1)=wRhsL(ix,iy,iz+1,1)-Flux
        wRhsL(ix,iy,iz  ,1)=wRhsL(ix,iy,iz  ,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Two*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsR(ix,iy,iz+1,1)=wRhsR(ix,iy,iz+1,1)-Flux
        wRhsR(ix,iy,iz  ,1)=wRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  ELSE
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DV(ix,iy,iz,1)
        FF=FW(ix,iy,iz)
        cL=wCL(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCL(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsL(ix,iy,iz  ,1)=wRhsL(ix,iy,iz  ,1)+Flux
        cL=wCR(ix,iy,iz,1)/(Rho(ix,iy,iz,1)+Eps)
        cR=wCR(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)
        Flux=  &
             Four*DRhoF*FF*FF/(VL+VR+Eps)*(cR-cL)
        wRhsR(ix,iy,iz  ,1)=wRhsR(ix,iy,iz  ,1)+Flux
      END DO
    END DO
  END IF

  DO iz=iz0,iz1+1
    DO iy=iy0,iy1+1
      DO ix=ix0,ix1+1
        wRhsL(ix,iy,iz,1)=wRhsL(ix,iy,iz,1)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        wRhsR(ix,iy,iz,1)=wRhsR(ix,iy,iz,1)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE DiffusionComputeW

SUBROUTINE JacDiffusionCompute

  INTEGER :: i
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: VL,VR
  REAL(RealKind) :: DL,DR
  REAL(RealKind) :: FF,DRhoF
  REAL(RealKind) :: Temp
  INTEGER :: nxP2,nyP2,nzP2

  nxP2=ix1-ix0+2
  nyP2=iy1-iy0+2
  nzP2=iz1-iz0+2

! Dichte Einarbeiten -------------

! x-Direction
  IF (TypeW(1:1)/='o') THEN
    ix=ix0
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_i
        i=Index(ix,iy,iz)
        AT%Val(i,3)=AT%Val(i,3)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_i+1
        i=Index(ix+1,iy,iz)
        AT%Val(i,5)=AT%Val(i,5)+Temp/(VL+Eps)/(Rho(ix+1,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix+1,iy,iz,1)+Eps)
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_i
        i=Index(ix,iy,iz)
        AT%Val(i,3)=AT%Val(i,3)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_i+1
        i=Index(ix+1,iy,iz)
        AT%Val(i,5)=AT%Val(i,5)+Temp/(VL+Eps)/(Rho(ix+1,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix+1,iy,iz,1)+Eps)
      END DO
    END DO
  END DO
  IF (TypeE(1:1)/='o') THEN
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix+1,iy,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix+1,iy,iz,1),VR)
        FF=FU(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_i
        i=Index(ix,iy,iz)
        AT%Val(i,3)=AT%Val(i,3)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_i+1
        i=Index(ix+1,iy,iz)
        AT%Val(i,5)=AT%Val(i,5)+Temp/(VL+Eps)/(Rho(ix+1,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix+1,iy,iz,1)+Eps)
      END DO
    END DO
  END IF
! y-Direction
  IF (TypeS(1:1)/='o') THEN
    iy=iy0
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_j
        i=Index(ix,iy,iz)
        AT%Val(i,2)=AT%Val(i,2)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_j+1
        i=Index(ix,iy+1,iz)
        AT%Val(i,6)=AT%Val(i,6)+Temp/(VL+Eps)/(Rho(ix,iy+1,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix,iy+1,iz,1)+Eps)
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_j
        i=Index(ix,iy,iz)
        AT%Val(i,2)=AT%Val(i,2)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_j+1
        i=Index(ix,iy+1,iz)
        AT%Val(i,6)=AT%Val(i,6)+Temp/(VL+Eps)/(Rho(ix,iy+1,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix,iy+1,iz,1)+Eps)
      END DO
    END DO
  END DO
  IF (TypeN(1:1)/='o') THEN
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy+1,iz)
        DRhoF=DF(DH(ix,iy,iz,1),VL,DH(ix,iy+1,iz,1),VR)
        FF=FV(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_j
        i=Index(ix,iy,iz)
        AT%Val(i,2)=AT%Val(i,2)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_j+1
        i=Index(ix,iy+1,iz)
        AT%Val(i,6)=AT%Val(i,6)+Temp/(VL+Eps)/(Rho(ix,iy+1,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix,iy+1,iz,1)+Eps)
      END DO
    END DO
  END IF

! z-Direction
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_k
        i=Index(ix,iy,iz)
        AT%Val(i,1)=AT%Val(i,1)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_k+1
        i=Index(ix,iy,iz+1)
        AT%Val(i,7)=AT%Val(i,7)+Temp/(VL+Eps)/(Rho(ix,iy,iz+1,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix,iy,iz+1,1)+Eps)
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        i=Index(ix,iy,iz)
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_k
        i=Index(ix,iy,iz)
        AT%Val(i,1)=AT%Val(i,1)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_k+1
        i=Index(ix,iy,iz+1)
        AT%Val(i,7)=AT%Val(i,7)+Temp/(VL+Eps)/(Rho(ix,iy,iz+1,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix,iy,iz+1,1)+Eps)
      END DO
    END DO
  END DO
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VL=VolC(ix,iy,iz)
        VR=VolC(ix,iy,iz+1)
        DRhoF=DF(DV(ix,iy,iz,1),VL,DV(ix,iy,iz+1,1),VR)
        FF=FW(ix,iy,iz)
        Temp=Two*DRhoF*FF*FF/(VL+VR+Eps)
!       Ableitung nach c_k
        i=Index(ix,iy,iz)
        AT%Val(i,1)=AT%Val(i,1)+Temp/(VR+Eps)/(Rho(ix,iy,iz,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VL+Eps)/(Rho(ix,iy,iz,1)+Eps)
!       Ableitung nach c_k+1
        i=Index(ix,iy,iz+1)
        AT%Val(i,7)=AT%Val(i,7)+Temp/(VL+Eps)/(Rho(ix,iy,iz+1,1)+Eps)
        AT%Val(i,4)=AT%Val(i,4)-Temp/(VR+Eps)/(Rho(ix,iy,iz+1,1)+Eps)
      END DO
    END DO
  END IF

CONTAINS
  FUNCTION Index(ix,iy,iz)
    INTEGER :: Index,ix,iy,iz
    Index=ix-ix0+1+nxP2*(iy-iy0)+nxP2*nyP2*(iz-iz0)  
  END FUNCTION Index
END SUBROUTINE JacDiffusionCompute

END MODULE Diffusion_Mod
