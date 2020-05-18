MODULE Advection_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Operator_Mod
  USE Thermodynamic_Mod
  USE ParameterMicrophys_Mod
  USE Control_Mod

  IMPLICIT NONE

CONTAINS

SUBROUTINE AdvectionCompute(Phi,FallF,Species)

  CHARACTER(2),OPTIONAL :: Species
  INTERFACE
    FUNCTION Phi1(r,k1,k2)
      USE Kind_Mod
      REAL(RealKind) :: Phi1,r,k1,k2
    END FUNCTION Phi1
    FUNCTION Phi(r,k1,k2,cL,cC,cR,hC)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2,cL,cC,cR,hC
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  IF (Primitive) THEN
    CALL AdvectionComputeP(Phi,FallF,Species)
  ELSE 
    CALL AdvectionComputeC(Phi,FallF,Species)
  END IF

END SUBROUTINE AdvectionCompute

SUBROUTINE AdvectionPreCompute(Phi,FallF)
  INTERFACE
    FUNCTION Phi(r,k1,k2,cL,cC,cR,hC)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2,cL,cC,cR,hC
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  IF (AdvPr=='V1') THEN
    CALL AdvectionPreComputeC(FallF)
  ELSE IF (AdvPr=='V2') THEN
    CALL AdvectionPreComputeV2(Phi,FallF)
  ELSE IF (AdvPr=='V3') THEN
    CALL AdvectionPreComputeV3(Phi,FallF)
  END IF  
END SUBROUTINE AdvectionPreCompute

SUBROUTINE AdvectionPreComputeV2(Phi,FallF)

  INTERFACE
    FUNCTION Phi(r,k1,k2,cR,cC,cL,hC)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2,cR,cC,cL,hC
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cL,cC,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux
  REAL(RealKind) :: cDiff,cLeft,cRight
  REAL(RealKind) :: uLe,uRi,vLe,vRi,wLe,wRi
  REAL(RealKind) :: uLoc,vLoc,wLoc
  REAL(RealKind) :: RhoL,RhoC,RhoR
  REAL(RealKind) :: RhoLF

  DO it=LBOUND(c,4),UBOUND(c,4)
! x-Direction
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0,ix1+1
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        uLoc=Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(Rho(ix,iy,iz,1)+Eps)
!       Left Value        
!       uLoc>0 take for interpolation cLL,cL,cC
!       uLoc<0 take for interpolation cL,cC,cR
!       Right Value        
!       uLoc>0 take for interpolation cL,cC,cR
!       uLoc<0 take for interpolation cC,cR,cRR
        cLeft=Zero
        cRight=Zero
        IF (uLoc>Zero.AND.ix>ix0) THEN
!         uLoc>0 take for interpolation cLL,cL,cC
!         Left boundary
          IF (ix>ix0+1) THEN
            VL=VolC(ix-2,iy,iz)
            cL=c(ix-2,iy,iz,it)
            VC=VolC(ix-1,iy,iz)
            cC=c(ix-1,iy,iz,it)
            VR=VolC(ix,iy,iz)
            cR=c(ix,iy,iz,it)
            FL=FU(ix-2,iy,iz)
            FR=FU(ix-1,iy,iz)
            hL=hLF(VL,FL)
            hC=hCF(VC,FL,FR)
            hR=hRF(VR,FR)
            k1=k1F(VL,VC,VR)
            k2=k2F(VL,VC,VR)
            cDiff=cC-cL
            r=(cR-cC+Eps)/(cDiff+Eps)
            cLeft=(cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
          END IF  
!         Right boundary
          IF (ix<ix1+1) THEN
            VL=VolC(ix-1,iy,iz)
            cL=c(ix-1,iy,iz,it)
            VC=VolC(ix,iy,iz)
            cC=c(ix,iy,iz,it)
            VR=VolC(ix+1,iy,iz)
            cR=c(ix+1,iy,iz,it)
            FL=FU(ix-1,iy,iz)
            FR=FU(ix,iy,iz)
            hL=hLF(VL,FL)
            hC=hCF(VC,FL,FR)
            hR=hRF(VR,FR)
            k1=k1F(VL,VC,VR)
            k2=k2F(VL,VC,VR)
            cDiff=cC-cL
            r=(cR-cC+Eps)/(cDiff+Eps)
            cRight=cC+FR/(FR+Eps)*phi(r,k1,k2,cR,cC,cL,hC)
          END IF  
        END IF  
        IF (uLoc<Zero.AND.ix<ix1+1) THEN
!         uLoc<0 take for interpolation cL,cC,cR
!         Left boundary          
          IF (ix>ix0) THEN
            VL=VolC(ix-1,iy,iz)
            cL=c(ix-1,iy,iz,it)
            VC=VolC(ix,iy,iz)
            cC=c(ix,iy,iz,it)
            VR=VolC(ix+1,iy,iz)
            cR=c(ix+1,iy,iz,it)
            FL=FU(ix-1,iy,iz)
            FR=FU(ix,iy,iz)
            hL=hLF(VL,FL)
            hC=hCF(VC,FL,FR)
            hR=hRF(VR,FR)
            k1=k1F(VR,VC,VL)
            k2=k2F(VR,VC,VL)
            cDiff=cR-cC
            r=(cC-cL+Eps)/(cDiff+Eps)
            cLeft=cC+FL/(FL+Eps)*phi(r,k1,k2,cL,cC,cR,hC)       
          END IF  
!         Right boundary
          IF (ix<ix1) THEN
            VL=VolC(ix,iy,iz)
            cL=c(ix,iy,iz,it)
            VC=VolC(ix+1,iy,iz)
            cC=c(ix+1,iy,iz,it)
            VR=VolC(ix+2,iy,iz)
            cR=c(ix+2,iy,iz,it)
            FL=FU(ix,iy,iz)
            FR=FU(ix+1,iy,iz)
            hL=hLF(VL,FL)
            hC=hCF(VC,FL,FR)
            hR=hRF(VR,FR)
            k1=k1F(VR,VC,VL)
            k2=k2F(VR,VC,VL)
            cDiff=cR-cC
            r=(cC-cL+Eps)/(cDiff+Eps)
            cRight=cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC)       
          END IF  
        END IF  
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-uLoc*(cRight-cLeft)*Half*(FL+FR)
      END DO  
    END DO  
  END DO  
! z-Direction
  DO iz=iz0,iz1+1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        wLoc=Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(Rho(ix,iy,iz,1)+Eps)
!       Left Value        
!       wLoc>0 take for interpolation cLL,cL,cC
!       wLoc<0 take for interpolation cL,cC,cR
!       Right Value        
!       wLoc>0 take for interpolation cL,cC,cR
!       wLoc<0 take for interpolation cC,cR,cRR
        cLeft=Zero
        cRight=Zero
        IF (wLoc>Zero.AND.iz>iz0) THEN
!         wLoc>0 take for interpolation cLL,cL,cC
!         Left boundary
          IF (iz>iz0+1) THEN
            VL=VolC(ix,iy,iz-2)
            cL=c(ix,iy,iz-2,it)
            VC=VolC(ix,iy,iz-1)
            cC=c(ix,iy,iz-1,it)
            VR=VolC(ix,iy,iz)
            cR=c(ix,iy,iz,it)
            FL=FW(ix,iy,iz-2)
            FR=FW(ix,iy,iz-1)
            hL=hLF(VL,FL)
            hC=hCF(VC,FL,FR)
            hR=hRF(VR,FR)
            k1=k1F(VL,VC,VR)
            k2=k2F(VL,VC,VR)
            cDiff=cC-cL
            r=(cR-cC+Eps)/(cDiff+Eps)
            cLeft=(cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
          END IF  
!         Right boundary
          IF (iz<iz1+1) THEN
            VL=VolC(ix,iy,iz-1)
            cL=c(ix,iy,iz-1,it)
            VC=VolC(ix,iy,iz)
            cC=c(ix,iy,iz,it)
            VR=VolC(ix,iy,iz+1)
            cR=c(ix,iy,iz+1,it)
            FL=FW(ix,iy,iz-1)
            FR=FW(ix,iy,iz)
            hL=hLF(VL,FL)
            hC=hCF(VC,FL,FR)
            hR=hRF(VR,FR)
            k1=k1F(VL,VC,VR)
            k2=k2F(VL,VC,VR)
            cDiff=cC-cL
            r=(cR-cC+Eps)/(cDiff+Eps)
            cRight=cC+FR/(FR+Eps)*phi(r,k1,k2,cR,cC,cL,hC)
          END IF  
        END IF  
        IF (wLoc<Zero.AND.iz<iz1+1) THEN
!         wLoc<0 take for interpolation cL,cC,cR
!         Left boundary          
          IF (iz>iz0) THEN
            VL=VolC(ix,iy,iz-1)
            cL=c(ix,iy,iz-1,it)
            VC=VolC(ix,iy,iz)
            cC=c(ix,iy,iz,it)
            VR=VolC(ix,iy,iz+1)
            cR=c(ix,iy,iz+1,it)
            FL=FW(ix,iy,iz-1)
            FR=FW(ix,iy,iz)
            hL=hLF(VL,FL)
            hC=hCF(VC,FL,FR)
            hR=hRF(VR,FR)
            k1=k1F(VR,VC,VL)
            k2=k2F(VR,VC,VL)
            cDiff=cR-cC
            r=(cC-cL+Eps)/(cDiff+Eps)
            cLeft=cC+FL/(FL+Eps)*phi(r,k1,k2,cL,cC,cR,hC)       
          END IF  
!         Right boundary
          IF (iz<iz1) THEN
            VL=VolC(ix,iy,iz)
            cL=c(ix,iy,iz,it)
            VC=VolC(ix,iy,iz+1)
            cC=c(ix,iy,iz+1,it)
            VR=VolC(ix,iy,iz+2)
            cR=c(ix,iy,iz+2,it)
            FL=FW(ix,iy,iz)
            FR=FW(ix,iy,iz+1)
            hL=hLF(VL,FL)
            hC=hCF(VC,FL,FR)
            hR=hRF(VR,FR)
            k1=k1F(VR,VC,VL)
            k2=k2F(VR,VC,VL)
            cDiff=cR-cC
            r=(cC-cL+Eps)/(cDiff+Eps)
            cRight=cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC)       
          END IF  
        END IF  
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-wLoc*(cRight-cLeft)*Half*(FL+FR)
      END DO  
    END DO  
  END DO  
  IF (TypeW=='ow') THEN
    ix=ix0+1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        uLoc=Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(Rho(ix,iy,iz,1)+Eps)
        cLeft=c(ix-1,iy,iz,it)
        IF (uLoc>Zero) THEN
          f(ix,iy,iz,it)=f(ix,iy,iz,it)-uLoc*(-cLeft)*Half*(FL+FR)
        END IF  
      END DO
    END DO
  END IF
  IF (TypeW=='oe') THEN
    ix=ix1
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        uLoc=Half*(uCL(ix,iy,iz,1)+uCR(ix,iy,iz,1))/(Rho(ix,iy,iz,1)+Eps)
        cRight=c(ix+1,iy,iz,it)
        IF (uLoc<Zero) THEN
          f(ix,iy,iz,it)=f(ix,iy,iz,it)-uLoc*(cRight)*Half*(FL+FR)
        END IF  
      END DO
    END DO
  END IF
  IF (TypeB=='ob') THEN
    iz=iz0+1
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        wLoc=Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(Rho(ix,iy,iz,1)+Eps)
        cLeft=c(ix,iy,iz-1,it)
        IF (FL==Zero) THEN
          cLeft=c(ix,iy,iz,it)
        END IF  
        IF (wLoc>Zero) THEN
          f(ix,iy,iz,it)=f(ix,iy,iz,it)-wLoc*(-cLeft)*Half*(FL+FR)
        END IF  
      END DO
    END DO
  END IF
  IF (TypeT=='ot') THEN
    iz=iz1
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        wLoc=Half*(wCL(ix,iy,iz,1)+wCR(ix,iy,iz,1))/(Rho(ix,iy,iz,1)+Eps)
        cRight=c(ix,iy,iz+1,it)
        IF (wLoc<Zero) THEN
          f(ix,iy,iz,it)=f(ix,iy,iz,it)-wLoc*(cRight)*Half*(FL+FR)
        END IF  
      END DO
    END DO
  END IF
  END DO !it

END SUBROUTINE AdvectionPreComputeV2

SUBROUTINE AdvectionPreComputeV3(Phi,FallF)

  INTERFACE
    FUNCTION Phi(r,k1,k2,cR,cC,cL,hC)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2,cR,cC,cL,hC
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cL,cC,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux
  REAL(RealKind) :: cDiff,cLeft,cRight
  REAL(RealKind) :: uLe,uRi,vLe,vRi,wLe,wRi
  REAL(RealKind) :: uLoc,vLoc,wLoc
  REAL(RealKind) :: RhoLF,wfLoc


  DO it=LBOUND(c,4),UBOUND(c,4)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cC=c(ix,iy,iz,it)
        VC=VolC(ix,iy,iz)

!       x-Direction
        cL=c(ix-1,iy,iz,it)
        cR=c(ix+1,iy,iz,it)
        VL=VolC(ix-1,iy,iz)
        VR=VolC(ix+1,iy,iz)
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(uF(ix-1,iy,iz))-uF(ix-1,iy,iz))* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))       
        f(ix-1,iy,iz,it)=f(ix-1,iy,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary 
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(uF(ix,iy,iz))+uF(ix,iy,iz))* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix+1,iy,iz,it)=f(ix+1,iy,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux

!       y-Direction
        cL=c(ix,iy-1,iz,it)
        cR=c(ix,iy+1,iz,it)
        VL=VolC(ix,iy-1,iz)
        VR=VolC(ix,iy+1,iz)
        FR=FV(ix,iy,iz)
        FL=FV(ix,iy-1,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(vF(ix,iy-1,iz))-vF(ix,iy-1,iz))* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))
        f(ix,iy-1,iz,it)=f(ix,iy-1,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(vF(ix,iy,iz))+vF(ix,iy,iz))* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix,iy+1,iz,it)=f(ix,iy+1,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux

!       z-Direction
        cL=c(ix,iy,iz-1,it)
        cR=c(ix,iy,iz+1,it)
        VL=VolC(ix,iy,iz-1)
        VR=VolC(ix,iy,iz+1)
        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        IF (PRESENT(FallF)) THEN
           RhoLF=CellToFaceVol(RhoR(ix,iy,iz-1,1),RhoR(ix,iy,iz,1) &
                           ,VL,VC)
          vFall=FallF(RhoLF)
        END IF
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz-1) &
             +CellToFaceVol(Rho(ix,iy,iz-1,1),Rho(ix,iy,iz,1) &
                           ,VolC(ix,iy,iz-1),VolC(ix,iy,iz))*(vFall+vSub(iz-1))
        Flux=Half*(ABS(wFLoc)-wFLoc)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))
        f(ix,iy,iz-1,it)=f(ix,iy,iz-1,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary
        IF (PRESENT(FallF)) THEN
           RhoLF=CellToFaceVol(RhoR(ix,iy,iz+1,1),RhoR(ix,iy,iz,1) &
                               ,VR,VC)
          vFall=FallF(RhoLF)
        END IF
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz) &
             +CellToFaceVol(Rho(ix,iy,iz+1,1),Rho(ix,iy,iz,1) &
                           ,VolC(ix,iy,iz+1),VolC(ix,iy,iz))*(vFall+vSub(iz))
        Flux=Half*(ABS(wFLoc)+wFLoc)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix,iy,iz+1,it)=f(ix,iy,iz+1,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
      END DO
    END DO
  END DO
  IF (TypeW=='ow') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix0,iy,iz))+uF(ix0,iy,iz))* &
              c(ix0,iy,iz,it)*FU(ix0,iy,iz)
        f(ix0+1,iy,iz,it)=f(ix0+1,iy,iz,it)+Flux
        f(ix0,iy,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeE=='oe') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix1,iy,iz))-uF(ix1,iy,iz))* &
              c(ix1+1,iy,iz,it)*FU(ix1,iy,iz) 
        f(ix1,iy,iz,it)=f(ix1,iy,iz,it)+Flux
        f(ix1+1,iy,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeS=='os') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy0,iz))+vF(ix,iy0,iz))* &
              c(ix,iy0,iz,it)*FV(ix,iy0,iz)
        f(ix,iy0+1,iz,it)=f(ix,iy0+1,iz,it)+Flux
        f(ix,iy0,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeN=='on') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy1,iz))-vF(ix,iy1,iz))* &
              c(ix,iy1+1,iz,it)*FV(ix,iy1,iz)
        f(ix,iy1,iz,it)=f(ix,iy1,iz,it)+Flux
        f(ix,iy1+1,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeB=='ob') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        IF (PRESENT(FallF)) THEN
          vFall=FallF(RhoR(ix,iy,iz0+1,1))
        END IF
        wFLoc=wF(ix,iy,iz0) &
             +Rho(ix,iy,iz0+1,1)*(vFall+vSub(iz0))
        Flux=Half*(ABS(wFLoc)+wFLoc)* &
              c(ix,iy,iz0,it)*FW(ix,iy,iz0)
        f(ix,iy,iz0+1,it)=f(ix,iy,iz0+1,it)+Flux 
        f(ix,iy,iz0,it)=Zero 
      END DO
    END DO
  END IF
  IF (TypeT=='ot') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        IF (PRESENT(FallF)) THEN
          vFall=FallF(RhoR(ix,iy,iz1,1))
        END IF
        wFLoc=wF(ix,iy,iz1) &
             +Rho(ix,iy,iz1,1)*(vFall+vSub(iz1))
        Flux=Half*(ABS(wFLoc)-wFLoc)* &
              c(ix,iy,iz1+1,it)*FW(ix,iy,iz1)
        f(ix,iy,iz1,it)=f(ix,iy,iz1,it)+Flux
        f(ix,iy,iz1+1,it)=Zero
      END DO
    END DO
  END IF
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        f(ix,iy,iz,it)=f(ix,iy,iz,it)+c(ix,iy,iz,it) &
                      *(FU(ix,iy,iz)*uF(ix,iy,iz) &
                       -FU(ix-1,iy,iz)*uF(ix-1,iy,iz) &
                       +FV(ix,iy,iz)*vF(ix,iy,iz) &
                       -FV(ix,iy-1,iz)*vF(ix,iy-1,iz) &
                       +FW(ix,iy,iz)*wF(ix,iy,iz) &
                       -FW(ix,iy,iz-1)*wF(ix,iy,iz-1)) 
      END DO
    END DO
  END DO
  DO iz=iz0,iz1+1
    DO iy=iy0,iy1+1
      DO ix=ix0,ix1+1
        f(ix,iy,iz,it)=f(ix,iy,iz,it)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)/(Rho(ix,iy,iz,1)+Eps)
      END DO
    END DO
  END DO
  END DO

END SUBROUTINE AdvectionPreComputeV3

FUNCTION RhoFace(RhoL,RhoR,VL,VR)

  REAL(RealKind) :: RhoFace
  REAL(RealKind) :: RhoL,RhoR,VL,VR

  RhoFace=(RhoL*VL+RhoR*VR)/(VR+VL+Eps)
END FUNCTION RhoFace

SUBROUTINE AdvectionPreComputeC(FallF)

  INTERFACE
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cC,cL,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux
  REAL(RealKind) :: cDiff
  REAL(RealKind) :: uLe,uRi,vLe,vRi,wLe,wRi
  REAL(RealKind) :: RhoLF

  DO it=LBOUND(c,4),UBOUND(c,4)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cC=c(ix,iy,iz,it)
        VC=VolC(ix,iy,iz)

!       x-Direction
        cL=c(ix-1,iy,iz,it)
        cR=c(ix+1,iy,iz,it)
        VL=VolC(ix-1,iy,iz)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix+1,iy,iz)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)
        uLe=uF(ix-1,iy,iz)
        uRi=uF(ix,iy,iz)
        Flux=FaceToCell(uLe*Gradient(cL,cC,FL,VL,VC),uRi*Gradient(cC,cR,FR,VC,VR),FL,FR) &
             *VC/(Rho(ix,iy,iz,1)+Eps)

!       y-Direction
        cL=c(ix,iy-1,iz,it)
        cR=c(ix,iy+1,iz,it)
        VL=VolC(ix,iy-1,iz)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix,iy+1,iz)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FV(ix,iy,iz)
        FL=FV(ix,iy-1,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)
        vLe=vF(ix,iy-1,iz)
        vRi=vF(ix,iy,iz)
        Flux=FaceToCell(vLe*Gradient(cL,cC,FL,VL,VC),vRi*Gradient(cC,cR,FR,VC,VR),FL,FR) &
             *VC/(Rho(ix,iy,iz,1)+Eps)
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux

!       z-Direction
        cL=c(ix,iy,iz-1,it)
        cR=c(ix,iy,iz+1,it)
        VL=VolC(ix,iy,iz-1)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix,iy,iz+1)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)
        wLe=wF(ix,iy,iz-1)
        wRi=wF(ix,iy,iz)
        Flux=FaceToCell(wLe*Gradient(cL,cC,FL,VL,VC),wRi*Gradient(cC,cR,FR,VC,VR),FL,FR) &
             *VC/(Rho(ix,iy,iz,1)+Eps)
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux

      END DO 
    END DO 
  END DO 
  END DO 

END SUBROUTINE AdvectionPreComputeC


SUBROUTINE AdvectionCompute4
  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: cLL,cL,cC,cR,cRR,cF
  REAL(RealKind) :: hLL,hL,hC,hR,hRR
  DO it=LBOUND(c,4),UBOUND(c,4)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+2,ix1-1
        cLL=c(ix-2,iy,iz,it)/(Rho(ix-2,iy,iz,it)+Eps)
        cL=c(ix-1,iy,iz,it)/(Rho(ix-1,iy,iz,it)+Eps)
        cC=c(ix,iy,iz,it)/(Rho(ix,iy,iz,it)+Eps)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,it)+Eps)
        cRR=c(ix+2,iy,iz,it)/(Rho(ix+2,iy,iz,it)+Eps)
        cF=(cLL+5.0d0*cL)
        f(ix,iy,iz,it)=f(ix,iy,iz,it)
        f(ix+1,iy,iz,it)=f(ix+1,iy,iz,it)
      END DO
    END DO
  END DO
  END DO
END SUBROUTINE AdvectionCompute4

SUBROUTINE AdvectionComputeP(Phi,FallF,Species)

  INTERFACE
    FUNCTION Phi1(r,k1,k2)
      USE Kind_Mod
      REAL(RealKind) :: Phi1,r,k1,k2
    END FUNCTION Phi1
    FUNCTION Phi(r,k1,k2,cL,cC,cR,hC)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2,cL,cC,cR,hC
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF
  CHARACTER(2),OPTIONAL :: Species

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cC,cL,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux
  REAL(RealKind) :: cDiff
  REAL(RealKind) :: wFLoc
  REAL(RealKind) :: RhoLF
  REAL(RealKind) :: RhoRF,RhoF,NRF,RhoIF,NIF
  REAL(RealKind) :: RhoSF,NSF
  REAL(RealKind) :: wFall(ix0+1:ix1,iy0+1:iy1,iz0:iz1,MAX(nFrac,1))

  wFall=0.0d0
  IF (PRESENT(Species)) THEN
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          SELECT CASE(Species)
          CASE('RhoR')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoRF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                               ,VL,VR)
            NRF=CellToFaceVol(Nrain(ix,iy,iz,1),Nrain(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz,1)=fall_rain_q(RhoRF,NRF,RhoF)!+vFall+vSub(zP(iz)) 
          CASE('NR')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoRF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                               ,VL,VR)
            NRF=CellToFaceVol(Nrain(ix,iy,iz,1),Nrain(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz,1)=fall_rain_n(RhoRF,NRF,RhoF)!+vFall+vSub(zP(iz))
          CASE('RhoI')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoIF=CellToFaceVol(RhoI(ix,iy,iz,1),RhoI(ix,iy,iz+1,1) &
                               ,VL,VR)
            NIF=CellToFaceVol(Nice(ix,iy,iz,1),Nice(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz,1)=fall_ice_q(RhoIF,NIF,RhoF)!+vFall+vSub(zP(iz))
          CASE('NI')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoIF=CellToFaceVol(RhoI(ix,iy,iz,1),RhoI(ix,iy,iz+1,1) &
                                ,VL,VR)
            NIF=CellToFaceVol(Nice(ix,iy,iz,1),Nice(ix,iy,iz+1,1) &
                                ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                ,VL,VR)
            wFall(ix,iy,iz,1)=fall_ice_n(RhoIF,NIF,RhoF)!+vFall+vSub(zP(iz))
          CASE('RhoS')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoSF=CellToFaceVol(RhoS(ix,iy,iz,1),RhoS(ix,iy,iz+1,1) &
                                ,VL,VR)
            NSF=CellToFaceVol(Nsnow(ix,iy,iz,1),Nsnow(ix,iy,iz+1,1) &
                                ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                ,VL,VR)
            wFall(ix,iy,iz,1)=fall_snow_q(RhoSF,NSF,RhoF)!+vFall+vSub(zP(iz))
          CASE('NS')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoSF=CellToFaceVol(RhoS(ix,iy,iz,1),RhoS(ix,iy,iz+1,1) &
                                 ,VL,VR)
            NSF=CellToFaceVol(Nsnow(ix,iy,iz,1),Nsnow(ix,iy,iz+1,1) &
                                 ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                 ,VL,VR)
            wFall(ix,iy,iz,1)=fall_snow_n(RhoSF,NSF,RhoF)!+vFall+vSub(zP(iz))
          CASE('AE') ! Aerosol
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            DO it=LBOUND(c,4),UBOUND(c,4)
              wFall(ix,iy,iz,it)=CellToFaceVol(VelSedi(ix,iy,iz,it),VelSedi(ix,iy,iz+1,it) &
                           ,VL,VR)
            END DO
          END SELECT                                                                                                                                                   
        END DO
      END DO
    END DO                                                                                                                                                            
  ELSE
    wFall=0.0d0
  END IF
  
  DO it=LBOUND(c,4),UBOUND(c,4)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cC=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VC=VolC(ix,iy,iz)

!       x-Direction
        cL=c(ix-1,iy,iz,it)/(Rho(ix-1,iy,iz,1)+Eps)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,1)+Eps)
        VL=VolC(ix-1,iy,iz)
        VR=VolC(ix+1,iy,iz)
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        IF (VL<VC) cL=((VC-VL)*cC+VL*cL)/VC
        IF (VR<VC) cR=((VC-VR)*cC+VR*cR)/VC
        IF (FL<FR) cL=((FR-FL)*cC+FL*cL)/FR
        IF (FR<FL) cR=((FL-FR)*cC+FR*cR)/FL
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(uF(ix-1,iy,iz))-uF(ix-1,iy,iz))* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))       
        f(ix-1,iy,iz,it)=f(ix-1,iy,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary 
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(uF(ix,iy,iz))+uF(ix,iy,iz))* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix+1,iy,iz,it)=f(ix+1,iy,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       y-Direction
        cL=c(ix,iy-1,iz,it)/(Rho(ix,iy-1,iz,1)+Eps)
        cR=c(ix,iy+1,iz,it)/(Rho(ix,iy+1,iz,1)+Eps)
        VL=VolC(ix,iy-1,iz)
        VR=VolC(ix,iy+1,iz)
        FR=FV(ix,iy,iz)
        FL=FV(ix,iy-1,iz)
!       IF (VL<VC) cL=((VC-VL)*cC+VL*cL)/VC
!       IF (VR<VC) cR=((VC-VR)*cC+VR*cR)/VC
        IF (FL<FR) cL=((FR-FL)*cC+FL*cL)/FR
        IF (FR<FL) cR=((FL-FR)*cC+FR*cR)/FL
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(vF(ix,iy-1,iz))-vF(ix,iy-1,iz))* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))
        f(ix,iy-1,iz,it)=f(ix,iy-1,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(vF(ix,iy,iz))+vF(ix,iy,iz))* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix,iy+1,iz,it)=f(ix,iy+1,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux

!       z-Direction
        cL=c(ix,iy,iz-1,it)/(Rho(ix,iy,iz-1,1)+Eps)
        cR=c(ix,iy,iz+1,it)/(Rho(ix,iy,iz+1,1)+Eps)
        VL=VolC(ix,iy,iz-1)
        VR=VolC(ix,iy,iz+1)
        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        IF (FL<FR) cL=((FR-FL)*cC+FL*cL)/FR
        IF (FR<FL) cR=((FL-FR)*cC+FR*cR)/FL
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Precipitation part
!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz-1) &
             +CellToFaceVol(Rho(ix,iy,iz-1,1),Rho(ix,iy,iz,1) &
                           ,VolC(ix,iy,iz-1),VolC(ix,iy,iz))*wFall(ix,iy,iz-1,it)
        Flux=Half*(ABS(wFLoc)-wFLoc)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))
        f(ix,iy,iz-1,it)=f(ix,iy,iz-1,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz) &
             +CellToFaceVol(Rho(ix,iy,iz+1,1),Rho(ix,iy,iz,1) &
                           ,VolC(ix,iy,iz+1),VolC(ix,iy,iz))*wFall(ix,iy,iz,it)
        Flux=Half*(ABS(wFLoc)+wFLoc)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix,iy,iz+1,it)=f(ix,iy,iz+1,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
      END DO
    END DO
  END DO
  IF (TypeW=='ow') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix0,iy,iz))+uF(ix0,iy,iz))* &
              c(ix0,iy,iz,it)/(Rho(ix0,iy,iz,1)+Eps)*FU(ix0,iy,iz)
        f(ix0+1,iy,iz,it)=f(ix0+1,iy,iz,it)+Flux
        f(ix0,iy,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeE=='oe') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix1,iy,iz))-uF(ix1,iy,iz))* &
              c(ix1+1,iy,iz,it)/(Rho(ix1+1,iy,iz,1)+Eps)*FU(ix1,iy,iz) 
        f(ix1,iy,iz,it)=f(ix1,iy,iz,it)+Flux
        f(ix1+1,iy,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeS=='os') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy0,iz))+vF(ix,iy0,iz))* &
              c(ix,iy0,iz,it)/(Rho(ix,iy0,iz,1)+Eps)*FV(ix,iy0,iz)
        f(ix,iy0+1,iz,it)=f(ix,iy0+1,iz,it)+Flux
        f(ix,iy0,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeN=='on') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy1,iz))-vF(ix,iy1,iz))* &
              c(ix,iy1+1,iz,it)/(Rho(ix,iy1+1,iz,1)+Eps)*FV(ix,iy1,iz)
        f(ix,iy1,iz,it)=f(ix,iy1,iz,it)+Flux
        f(ix,iy1+1,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeB=='ob') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       IF (PRESENT(FallF)) THEN
!         vFall=FallF(RhoR(ix,iy,iz0+1,1))
!       END IF
        wFLoc=wF(ix,iy,iz0) &
             +Rho(ix,iy,iz0+1,1)*wFall(ix,iy,iz0+1,it)
        Flux=Half*(ABS(wFLoc)+wFLoc)* &
              c(ix,iy,iz0,it)/(Rho(ix,iy,iz0+1,1)+Eps)*FW(ix,iy,iz0)
        f(ix,iy,iz0+1,it)=f(ix,iy,iz0+1,it)+Flux 
        f(ix,iy,iz0,it)=Zero 
      END DO
    END DO
  END IF
  IF (TypeT=='ot') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       IF (PRESENT(FallF)) THEN
!         vFall=FallF(RhoR(ix,iy,iz1,1))
!       END IF
        wFLoc=wF(ix,iy,iz1) &
             +Rho(ix,iy,iz1,1)*wFall(ix,iy,iz1,it)
        Flux=Half*(ABS(wFLoc)-wFLoc)* &
              c(ix,iy,iz1+1,it)/(Rho(ix,iy,iz1+1,1)+Eps)*FW(ix,iy,iz1)
        f(ix,iy,iz1,it)=f(ix,iy,iz1,it)+Flux
        f(ix,iy,iz1+1,it)=Zero
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

END SUBROUTINE AdvectionComputeP

SUBROUTINE AdvectionComputeC(Phi,FallF,Species)

  INTERFACE
    FUNCTION Phi1(r,k1,k2)
      USE Kind_Mod
      REAL(RealKind) :: Phi1,r,k1,k2
    END FUNCTION Phi1
    FUNCTION Phi(r,k1,k2,cL,cC,cR,hC)
      USE KIND_Mod
      REAL(RealKind) :: Phi,r,k1,k2,cL,cC,cR,hC
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF
  CHARACTER(2),OPTIONAL :: Species

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cC,cL,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux
  REAL(RealKind) :: cDiff
  REAL(RealKind) :: uFLoc,vFloc,wFLoc
  REAL(RealKind) :: RhoLe,RhoC,RhoRi
  REAL(RealKind) :: RhoLF,RhoRF,NRF,RhoF,RhoIF,NIF,RhoSF,NSF
  REAL(RealKind) :: wFall(ix0+1:ix1,iy0+1:iy1,iz0:iz1)

  IF (PRESENT(Species)) THEN
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          SELECT CASE(Species)
          CASE('RhoR')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoRF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                               ,VL,VR)
            NRF=CellToFaceVol(Nrain(ix,iy,iz,1),Nrain(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_rain_q(RhoRF,NRF,RhoF)!+vFall+vSub(zP(iz)) 
          CASE('NR')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoRF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                               ,VL,VR)
            NRF=CellToFaceVol(Nrain(ix,iy,iz,1),Nrain(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_rain_n(RhoRF,NRF,RhoF)!+vFall+vSub(zP(iz))
          CASE('RhoI')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoIF=CellToFaceVol(RhoI(ix,iy,iz,1),RhoI(ix,iy,iz+1,1) &
                               ,VL,VR)
            NIF=CellToFaceVol(Nice(ix,iy,iz,1),Nice(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_ice_q(RhoIF,NIF,RhoF)!+vFall+vSub(zP(iz))
          CASE('NI')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoIF=CellToFaceVol(RhoI(ix,iy,iz,1),RhoI(ix,iy,iz+1,1) &
                                ,VL,VR)
            NIF=CellToFaceVol(Nice(ix,iy,iz,1),Nice(ix,iy,iz+1,1) &
                                ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                ,VL,VR)
            wFall(ix,iy,iz)=fall_ice_n(RhoIF,NIF,RhoF)!+vFall+vSub(zP(iz))
          CASE('RhoS')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoSF=CellToFaceVol(RhoS(ix,iy,iz,1),RhoS(ix,iy,iz+1,1) &
                                ,VL,VR)
            NSF=CellToFaceVol(Nsnow(ix,iy,iz,1),Nsnow(ix,iy,iz+1,1) &
                                ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                ,VL,VR)
            wFall(ix,iy,iz)=fall_snow_q(RhoSF,NSF,RhoF)!+vFall+vSub(zP(iz))
          CASE('NS')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoSF=CellToFaceVol(RhoS(ix,iy,iz,1),RhoS(ix,iy,iz+1,1) &
                                 ,VL,VR)
            NSF=CellToFaceVol(Nsnow(ix,iy,iz,1),Nsnow(ix,iy,iz+1,1) &
                                 ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                 ,VL,VR)
            wFall(ix,iy,iz)=fall_snow_n(RhoSF,NSF,RhoF)!+vFall+vSub(zP(iz))
          END SELECT                                                                                                                                                   
        END DO
      END DO
    END DO                                                                                                                                                               
  ELSE
    wFall=0.0d0
  END IF

  DO it=LBOUND(c,4),UBOUND(c,4)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoC=Rho(ix,iy,iz,1)
        cC=c(ix,iy,iz,it)
        VC=VolC(ix,iy,iz)

!       x-Direction
        RhoLe=Rho(ix-1,iy,iz,1)
        cL=c(ix-1,iy,iz,it)
        RhoRi=Rho(ix+1,iy,iz,1)
        cR=c(ix+1,iy,iz,it)
        VL=VolC(ix-1,iy,iz)
        VR=VolC(ix+1,iy,iz)
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        uFLoc=uF(ix-1,iy,iz)*(VC+VL)/(VC*RhoC+VL*RhoLe+Eps)
        Flux=Half*(ABS(uFLoc)-uFLoc)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))
        f(ix-1,iy,iz,it)=f(ix-1,iy,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary 
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        uFLoc=uF(ix,iy,iz)*(VC+VR)/(VC*RhoC+VR*RhoRi+Eps)
        Flux=Half*(ABS(uFLoc)+uFLoc)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix+1,iy,iz,it)=f(ix+1,iy,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux

!       y-Direction
        RhoLe=Rho(ix,iy-1,iz,1)
        cL=c(ix,iy-1,iz,it)
        RhoRi=Rho(ix,iy+1,iz,1)
        cR=c(ix,iy+1,iz,it)
        VL=VolC(ix,iy-1,iz)
        VR=VolC(ix,iy+1,iz)
        FR=FV(ix,iy,iz)
        FL=FV(ix,iy-1,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        vFLoc=vF(ix,iy-1,iz)*(VC+VL)/(VC*RhoC+VL*RhoLe+Eps)
        Flux=Half*(ABS(vFLoc)-vFLoc)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))
        f(ix,iy-1,iz,it)=f(ix,iy-1,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        vFLoc=vF(ix,iy,iz)*(VC+VR)/(VC*RhoC+VR*RhoRi+Eps)
        Flux=Half*(ABS(vFLoc)+vFLoc)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix,iy+1,iz,it)=f(ix,iy+1,iz,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux

!       z-Direction
        RhoLe=Rho(ix,iy,iz-1,1)
        cL=c(ix,iy,iz-1,it)
        RhoRi=Rho(ix,iy,iz+1,1)
        cR=c(ix,iy,iz+1,it)
        VL=VolC(ix,iy,iz-1)
        VR=VolC(ix,iy,iz+1)
        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cR-cC
        r=(cC-cL+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz-1)*(VC+VL)/(VC*RhoC+VL*RhoLe+Eps)+wFall(ix,iy,iz-1)
        Flux=Half*(ABS(wFLoc)-wFLoc)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))
        f(ix,iy,iz-1,it)=f(ix,iy,iz-1,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
!       Right boundary
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz)*(VC+VR)/(VC*RhoC+VR*RhoRi+Eps)+wFall(ix,iy,iz)
        Flux=Half*(ABS(wFLoc)+wFLoc)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        f(ix,iy,iz+1,it)=f(ix,iy,iz+1,it)+Flux
        f(ix,iy,iz,it)=f(ix,iy,iz,it)-Flux
      END DO
    END DO
  END DO
  IF (TypeW=='ow') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix0,iy,iz))+uF(ix0,iy,iz))* &
              c(ix0,iy,iz,it)/(Rho(ix0,iy,iz,1)+Eps)*FU(ix0,iy,iz)
        f(ix0+1,iy,iz,it)=f(ix0+1,iy,iz,it)+Flux
        f(ix0,iy,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeE=='oe'.AND.BCvel%East/='OutFlow') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix1,iy,iz))-uF(ix1,iy,iz))* &
              c(ix1+1,iy,iz,it)/(Rho(ix1+1,iy,iz,1)+Eps)*FU(ix1,iy,iz) 
        f(ix1,iy,iz,it)=f(ix1,iy,iz,it)+Flux
        f(ix1+1,iy,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeS=='os') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy0,iz))+vF(ix,iy0,iz))* &
              c(ix,iy0,iz,it)/(Rho(ix,iy0,iz,1)+Eps)*FV(ix,iy0,iz)
        f(ix,iy0+1,iz,it)=f(ix,iy0+1,iz,it)+Flux
        f(ix,iy0,iz,it)=Zero
      END DO
    END DO
  END IF
  IF (TypeN=='on') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy1,iz))-vF(ix,iy1,iz))* &
              c(ix,iy1+1,iz,it)/(Rho(ix,iy1+1,iz,1)+Eps)*FV(ix,iy1,iz)
        f(ix,iy1,iz,it)=f(ix,iy1,iz,it)+Flux
        f(ix,iy1+1,iz,it)=Zero
      END DO
    END DO
  END IF 

  IF (TypeB=='ob') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       IF (PRESENT(FallF)) THEN
!         vFall=FallF(RhoR(ix,iy,iz0+1,1))
!       END IF
        wFLoc=wF(ix,iy,iz0) &
             +Rho(ix,iy,iz0+1,1)*wFall(ix,iy,iz0+1)
        Flux=Half*(ABS(wFLoc)+wFLoc)* &
              c(ix,iy,iz0,it)/(Rho(ix,iy,iz0+1,1)+Eps)*FW(ix,iy,iz0)
        f(ix,iy,iz0+1,it)=f(ix,iy,iz0+1,it)+Flux 
        f(ix,iy,iz0,it)=Zero 
      END DO
    END DO
  END IF
  IF (TypeT=='ot') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       IF (PRESENT(FallF)) THEN
!         vFall=FallF(RhoR(ix,iy,iz1,1))
!       END IF
        wFLoc=wF(ix,iy,iz1) &
             +Rho(ix,iy,iz1,1)*wFall(ix,iy,iz1)
        Flux=Half*(ABS(wFLoc)-wFLoc)* &
              c(ix,iy,iz1+1,it)/(Rho(ix,iy,iz1+1,1)+Eps)*FW(ix,iy,iz1)
        f(ix,iy,iz1,it)=f(ix,iy,iz1,it)+Flux
        f(ix,iy,iz1+1,it)=Zero
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

END SUBROUTINE AdvectionComputeC

SUBROUTINE AdvectionQECompute(Phi,FallF)

  INTERFACE
    FUNCTION Phi(r,k1,k2)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cC,cL,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux,Flux1,Flux2
  REAL(RealKind) :: PreFac 
  REAL(RealKind) :: cDiff
  REAL(RealKind) :: uFLoc,vFloc,wFLoc
  REAL(RealKind) :: RhoL,RhoC,RhoR
  REAL(RealKind) :: RhoLL,RhoCL,RhoRL
  REAL(RealKind) :: RhoLF

END SUBROUTINE AdvectionQECompute


SUBROUTINE AdvectionQFallCompute(Phi,FallF)

  INTERFACE
    FUNCTION Phi(r,k1,k2)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cC,cL,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux,Flux1,Flux2
  REAL(RealKind) :: cDiff
  REAL(RealKind) :: uFLoc,vFloc,wFLoc
  REAL(RealKind) :: RhoLe,RhoC,RhoRi
  REAL(RealKind) :: RhoLF
  REAL(RealKind) :: RhoLoc,RhoDLoc,RhoVLoc,RhoLLoc,TLoc,pLoc,ExnerPre
  REAL(RealKind) :: Cpml,Rm,Cp_eff,ThDens,ThEquiv,TotalEnergy 
  REAL(RealKind) :: PreFac 

  DO it=LBOUND(c,4),UBOUND(c,4)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoC=Rho(ix,iy,iz,1)
        cC=RhoR(ix,iy,iz,it)
        VC=VolC(ix,iy,iz)

!       z-Direction
        RhoLe=Rho(ix,iy,iz-1,1)
        cL=RhoR(ix,iy,iz-1,it)
        VL=VolC(ix,iy,iz-1)
        IF (VL==Zero) THEN
          cL=cC
        END IF

        RhoRi=Rho(ix,iy,iz+1,1)
        cR=RhoR(ix,iy,iz+1,it)
        VR=VolC(ix,iy,iz+1)
        IF (VR==Zero) THEN
          cR=cC
        END IF

        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
!       Compute d/dz([W_f+w]*RhoR) = Flux1
        IF (PRESENT(FallF)) THEN
          RhoLF=CellToFaceVol(RhoR(ix,iy,iz-1,1),RhoR(ix,iy,iz,1) &
                           ,VL,VC)
          vFall=FallF(RhoLF)
        END IF
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cC-cR
        r=(cL-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz-1)*(VC+VL)/(VC*RhoC+VL*RhoLe+Eps)+(vFall+vSub(iz-1))
        Flux1=Half*(ABS(wFLoc)-wFLoc)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2)*(cC-cR))
!       Compute d/dz(w*RhoR) = Flux2
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cC-cR
        r=(cL-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz-1)*(VC+VL)/(VC*RhoC+VL*RhoLe+Eps)
        Flux2=Half*(ABS(wFLoc)-wFLoc)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2)*(cC-cR))
!       Compute Flux
        Flux=Flux1-Flux2
        IF (ic==RhoPos) THEN
          RhoRhs(ix,iy,iz-1,it)=RhoRhs(ix,iy,iz-1,it)+Flux
          RhoRhs(ix,iy,iz,it)=RhoRhs(ix,iy,iz,it)-Flux
        ELSE IF (ic==ThPos) THEN
          RhoLoc=Rho(ix,iy,iz,1)+Eps
          RhoVLoc=RhoV(ix,iy,iz,1)+Eps
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)+Eps
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Three*Eps
          TLoc=T(ix,iy,iz,1)+Eps
          pLoc=p(ix,iy,iz,1)+Eps
          Rm=Rd+RhoVLoc/RhoDLoc*Rv
          Cpml=Cpd+RhoVLoc/RhoDLoc*Cpv+RhoLLoc/RhoDLoc*Cpl
          ExnerPre=(pLoc/p0)**(Rm/Cpml)
          SELECT CASE(ThetaKind)
          CASE('Density')
            ThDens=TLoc/ExnerPre*(1+Rv/Rd*RhoVLoc/RhoDLoc)/(1+RhoVLoc/RhoDLoc+RhoLLoc/RhoDLoc)
            PreFac=RhoLoc*ThDens*(Cpl/(RhoDLoc*Cpml)-LOG(ExnerPre)/RhoDLoc*(Cpl/Cpml)+1.0d0/RhoLoc)
          CASE('Equiv')
            Cp_eff=Cpd+(RhoVLoc/RhoDLoc+RhoLLoc/RhoDLoc)*Cpl
            ThEquiv=TLoc*(p0/(RhoDLoc*Rd*TLoc))**(Rd/Cp_eff) &
                     *RelHumidity(TLoc,RhoVLoc)**(-RhoVLoc/RhoDLoc*Rv/Cp_eff) &
                     *EXP(RhoDLoc*LatHeat(TLoc)/(RhoDLoc*Cp_eff*TLoc)) 
            PreFac=-RhoLoc*ThEquiv*Cpl/(RhoDLoc*Cp_eff**2) &
                    *(Rd*LOG(p0/pLoc*(1+Rv/Rd*RhoVLoc/RhoDLoc)-RhoVLoc) &
                      -RhoVLoc/RhoDLoc*Rv*LOG(RelHumidity(TLoc,RhoVLoc)) &
                      +RhoVLoc/RhoDLoc*LatHeat(TLoc)/TLoc &
                     ) &
                   -ThEquiv
          END SELECT
          ThRhs(ix,iy,iz+1,it)=ThRhs(ix,iy,iz+1,it)+PreFac*Flux
          ThRhs(ix,iy,iz,it)=ThRhs(ix,iy,iz,it)-PreFac*Flux
        END IF

!       Right boundary
!       Compute d/dz([W_f+w]*RhoR) = Flux1
        IF (PRESENT(FallF)) THEN
           RhoLF=CellToFaceVol(RhoR(ix,iy,iz+1,1),RhoR(ix,iy,iz,1) &
                               ,VR,VC)
          vFall=FallF(RhoLF)
        END IF
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz)*(VC+VR)/(VC*RhoC+VR*RhoRi+Eps)+(vFall+vSub(iz))
        Flux1=Half*(ABS(wFLoc)+wFLoc)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2)*(cC-cL))
!       Compute d/dz(w*RhoR) = Flux2
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz)*(VC+VR)/(VC*RhoC+VR*RhoRi+Eps)
        Flux2=Half*(ABS(wFLoc)+wFLoc)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2)*(cC-cL))
!       Compute Flux
        Flux=Flux1-Flux2
        IF (ic==RhoPos) THEN
          RhoRhs(ix,iy,iz+1,it)=RhoRhs(ix,iy,iz+1,it)+Flux
          RhoRhs(ix,iy,iz,it)=RhoRhs(ix,iy,iz,it)-Flux
        ELSE IF (ic==ThPos) THEN
          RhoLoc=Rho(ix,iy,iz,1)+Eps
          RhoVLoc=RhoV(ix,iy,iz,1)+Eps
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)+Eps
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Three*Eps
          TLoc=T(ix,iy,iz,1)+Eps
          pLoc=p(ix,iy,iz,1)+Eps
          Rm=Rd+RhoVLoc/RhoDLoc*Rv
          Cpml=Cpd+RhoVLoc/RhoDLoc*Cpv+RhoLLoc/RhoDLoc*Cpl
          ExnerPre=(pLoc/p0)**(Rm/Cpml)
          SELECT CASE(ThetaKind)
          CASE('Density')
            ThDens=TLoc/ExnerPre*(1+Rv/Rd*RhoVLoc/RhoDLoc)/(1+RhoVLoc/RhoDLoc+RhoLLoc/RhoDLoc)
            PreFac=RhoLoc*ThDens*(Cpl/(RhoDLoc*Cpml)-LOG(ExnerPre)/RhoDLoc*(Cpl/Cpml)+1.0d0/RhoLoc)
          CASE('Equiv')
            Cp_eff=Cpd+(RhoVLoc/RhoDLoc+RhoLLoc/RhoDLoc)*Cpl
            ThEquiv=TLoc*(p0/(RhoDLoc*Rd*TLoc))**(Rd/Cp_eff) &
                     *RelHumidity(TLoc,RhoVLoc)**(-RhoVLoc/RhoDLoc*Rv/Cp_eff) &
                     *EXP(RhoDLoc*LatHeat(TLoc)/(RhoDLoc*Cp_eff*TLoc)) 
            PreFac=-ThEquiv*Cpl/(RhoDLoc*Cp_eff**2) &
                    *(Rd*LOG(p0/pLoc*(1+Rv/Rd*RhoVLoc/RhoDLoc)-RhoVLoc) &
                      -RhoVLoc/RhoDLoc*Rv*LOG(RelHumidity(TLoc,RhoVLoc)) &
                      +RhoVLoc/RhoDLoc*LatHeat(TLoc)/TLoc &
                     )
          END SELECT
          ThRhs(ix,iy,iz+1,it)=ThRhs(ix,iy,iz+1,it)+PreFac*Flux
          ThRhs(ix,iy,iz,it)=ThRhs(ix,iy,iz,it)-PreFac*Flux
        END IF
      END DO
    END DO
  END DO 
  END DO !it

END SUBROUTINE AdvectionQFallCompute

SUBROUTINE AdvectionFaceCompute1(Phi,FallF)

  INTERFACE
    FUNCTION Phi(r,k1,k2)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2
    END FUNCTION Phi
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cC,cL,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux
  REAL(RealKind) :: cDiff
  REAL(RealKind) :: wFLoc
  REAL(RealKind) :: RhoLF

  cFU=Zero
  cFV=Zero
  cFW=Zero
! Fluxes
  DO it=LBOUND(c,4),UBOUND(c,4)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cC=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VC=VolC(ix,iy,iz)

!       x-Direction
        cL=c(ix-1,iy,iz,it)/(Rho(ix-1,iy,iz,1)+Eps)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,1)+Eps)
        VL=VolC(ix-1,iy,iz)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix+1,iy,iz)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cC-cR
        r=(cL-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(uF(ix-1,iy,iz)+Eps)-uF(ix-1,iy,iz)-Eps)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2)*(cC-cR))
        cFU(ix-1,iy,iz,it)=cFU(ix-1,iy,iz,it)+Flux/(FL*ABS(uF(ix-1,iy,iz)+Eps)+Eps**2)
!       Right boundary 
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(uF(ix,iy,iz)+Eps)+uF(ix,iy,iz)+Eps)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2)*(cC-cL))
        cFU(ix,iy,iz,it)=cFU(ix,iy,iz,it)+Flux/(FR*ABS(uF(ix,iy,iz)+Eps)+Eps**2)

!       y-Direction
        cL=c(ix,iy-1,iz,it)/(Rho(ix,iy-1,iz,1)+Eps)
        cR=c(ix,iy+1,iz,it)/(Rho(ix,iy+1,iz,1)+Eps)
        VL=VolC(ix,iy-1,iz)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix,iy+1,iz)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FV(ix,iy,iz)
        FL=FV(ix,iy-1,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cC-cR
        r=(cL-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(vF(ix,iy-1,iz)+Eps)-vF(ix,iy-1,iz)-Eps)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2)*(cC-cR))
        cFV(ix,iy-1,iz,it)=cFV(ix,iy-1,iz,it)+Flux/(FL*ABS(vF(ix,iy-1,iz)+Eps)+Eps**2)
!       Right boundary
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(vF(ix,iy,iz)+Eps)+vF(ix,iy,iz)+Eps)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2)*(cC-cL))
        cFV(ix,iy,iz,it)=cFV(ix,iy,iz,it)+Flux/(FR*ABS(vF(ix,iy,iz)+Eps)+Eps**2)

!       z-Direction
        cL=c(ix,iy,iz-1,it)/(Rho(ix,iy,iz-1,1)+Eps)
        cR=c(ix,iy,iz+1,it)/(Rho(ix,iy,iz+1,1)+Eps)
        VL=VolC(ix,iy,iz-1)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix,iy,iz+1)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cC-cR
        r=(cL-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz-1) 
        Flux=Half*(ABS(wFLoc+Eps)-wFLoc-Eps)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2)*(cC-cR))
        cFW(ix,iy,iz-1,it)=cFW(ix,iy,iz-1,it)+Flux/(FL*ABS(wFLoc+Eps)+Eps**2)
!       Right boundary
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz) 
        Flux=Half*(ABS(wFLoc+Eps)+wFLoc+Eps)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2)*(cC-cL))
        cFW(ix,iy,iz,it)=cFW(ix,iy,iz,it)+Flux/FR*(ABS(wFLoc+Eps)+Eps**2)
      END DO
    END DO
  END DO
  IF (TypeW=='ow') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix0,iy,iz)+Eps)+uF(ix0,iy,iz)+Eps)* &
              c(ix0,iy,iz,it)/(Rho(ix0,iy,iz,1)+Eps)*FU(ix0,iy,iz)
        cFU(ix0,iy,iz,it)=cFU(ix0,iy,iz,it)+Flux/(FU(ix0,iy,iz)*ABS(uF(ix0,iy,iz)+Eps)+Eps**2)
      END DO
    END DO
  END IF
  IF (TypeE=='oe') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix1,iy,iz)+Eps)-uF(ix1,iy,iz)-Eps)* &
              c(ix1+1,iy,iz,it)/(Rho(ix1+1,iy,iz,1)+Eps)*FU(ix1,iy,iz) 
        cFU(ix1,iy,iz,it)=cFU(ix1,iy,iz,it)+Flux/(FU(ix1,iy,iz)*ABS(uF(ix1,iy,iz)+Eps)+Eps**2)
      END DO
    END DO
  END IF
  IF (TypeS=='os') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy0,iz)+Eps)+vF(ix,iy0,iz)+Eps)* &
              c(ix,iy0,iz,it)/(Rho(ix,iy0,iz,1)+Eps)*FV(ix,iy0,iz)
        cFV(ix,iy0,iz,it)=cFV(ix,iy0,iz,it)+Flux/(FV(ix,iy0,iz)*ABS(vF(ix,iy0,iz)+Eps)+Eps**2)
      END DO
    END DO
  END IF
  IF (TypeN=='on') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy1,iz)+Eps)-vF(ix,iy1,iz)-Eps)* &
              c(ix,iy1+1,iz,it)/(Rho(ix,iy1+1,iz,1)+Eps)*FV(ix,iy1,iz)
        cFV(ix,iy1,iz,it)=cFV(ix,iy1,iz,it)+Flux/(FV(ix,iy1,iz)*ABS(vF(ix,iy1,iz)+Eps)+Eps**2)
      END DO
    END DO
  END IF
  IF (TypeB=='ob') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        IF (PRESENT(FallF)) THEN
          vFall=FallF(RhoR(ix,iy,iz0+1,1))
        END IF
        wFLoc=wF(ix,iy,iz0) &
             +Rho(ix,iy,iz0+1,1)*(vFall+vSub(iz0))
        Flux=Half*(ABS(wFLoc+Eps)+wFLoc+Eps)* &
              c(ix,iy,iz0,it)/(Rho(ix,iy,iz0+1,1)+Eps)*FW(ix,iy,iz0)
        cFW(ix,iy,iz0,it)=cFW(ix,iy,iz0,it)+Flux/(FW(ix,iy,iz0)*ABS(wFLoc+Eps)+Eps**2)
      END DO
    END DO
  END IF
  IF (TypeT=='ot') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        IF (PRESENT(FallF)) THEN
          vFall=FallF(RhoR(ix,iy,iz1,1))
        END IF
        wFLoc=wF(ix,iy,iz1) &
             +Rho(ix,iy,iz1,1)*(vFall+vSub(iz1))
        Flux=Half*(ABS(wFLoc+Eps)-wFLoc-Eps)* &
              c(ix,iy,iz1+1,it)/(Rho(ix,iy,iz1+1,1)+Eps)*FW(ix,iy,iz1)
        cFW(ix,iy,iz1,it)=cFW(ix,iy,iz1,it)+Flux/(FW(ix,iy,iz1)*ABS(wFLoc+Eps)+Eps**2)
      END DO
    END DO
  END IF
  END DO

END SUBROUTINE AdvectionFaceCompute1

SUBROUTINE AdvectionFaceCompute(Phi,Species)

  INTERFACE
    FUNCTION Phi(r,k1,k2,cL,cC,cR,hC)
      USE Kind_Mod
      REAL(RealKind) :: Phi,r,k1,k2,cL,cC,cR,hC
    END FUNCTION Phi
  END INTERFACE

  CHARACTER(2),OPTIONAL :: Species

  INTEGER :: ix,iy,iz,it
  REAL(RealKind) :: r
  REAL(RealKind) :: hL,hC,hR
  REAL(RealKind) :: k1,k2
  REAL(RealKind) :: cx,cy,cz
  REAL(RealKind) :: cC,cL,cR
  REAL(RealKind) :: VC,VL,VR
  REAL(RealKind) :: FL,FR
  REAL(RealKind) :: Flux
  REAL(RealKind) :: cDiff
  REAL(RealKind) :: wFLoc
  REAL(RealKind) :: RhoLF
  REAL(RealKind) :: RhoRF,RhoF,NRF,RhoIF,NIF
  REAL(RealKind) :: RhoSF,NSF
  REAL(RealKind) :: wFall(ix0+1:ix1,iy0+1:iy1,iz0:iz1)

  IF (PRESENT(Species)) THEN
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          SELECT CASE(Species)
          CASE('RhoR')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoRF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                               ,VL,VR)
            NRF=CellToFaceVol(Nrain(ix,iy,iz,1),Nrain(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_rain_q(RhoRF,NRF,RhoF)!+vFall+vSub(zP(iz)) 
          CASE('NR')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoRF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                               ,VL,VR)
            NRF=CellToFaceVol(Nrain(ix,iy,iz,1),Nrain(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_rain_n(RhoRF,NRF,RhoF)!+vFall+vSub(zP(iz))
          CASE('RhoI')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoIF=CellToFaceVol(RhoI(ix,iy,iz,1),RhoI(ix,iy,iz+1,1) &
                               ,VL,VR)
            NIF=CellToFaceVol(Nice(ix,iy,iz,1),Nice(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_ice_q(RhoIF,NIF,RhoF)!+vFall+vSub(zP(iz))
          CASE('NI')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoIF=CellToFaceVol(RhoI(ix,iy,iz,1),RhoI(ix,iy,iz+1,1) &
                                ,VL,VR)
            NIF=CellToFaceVol(Nice(ix,iy,iz,1),Nice(ix,iy,iz+1,1) &
                                ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                ,VL,VR)
            wFall(ix,iy,iz)=fall_ice_n(RhoIF,NIF,RhoF)!+vFall+vSub(zP(iz))
          CASE('RhoS')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoSF=CellToFaceVol(RhoS(ix,iy,iz,1),RhoS(ix,iy,iz+1,1) &
                                ,VL,VR)
            NSF=CellToFaceVol(Nsnow(ix,iy,iz,1),Nsnow(ix,iy,iz+1,1) &
                                ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                ,VL,VR)
            wFall(ix,iy,iz)=fall_snow_q(RhoSF,NSF,RhoF)!+vFall+vSub(zP(iz))
          CASE('NS')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoSF=CellToFaceVol(RhoS(ix,iy,iz,1),RhoS(ix,iy,iz+1,1) &
                                 ,VL,VR)
            NSF=CellToFaceVol(Nsnow(ix,iy,iz,1),Nsnow(ix,iy,iz+1,1) &
                                 ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                 ,VL,VR)
            wFall(ix,iy,iz)=fall_snow_n(RhoSF,NSF,RhoF)!+vFall+vSub(zP(iz))
          END SELECT                                                                                                                                                   
        END DO
      END DO
    END DO                                                                                                                                                               
  ELSE
    wFall=0.0d0
  END IF

  cFU=Zero
  cFV=Zero
  cFW=Zero
! Fluxes
  DO it=LBOUND(c,4),UBOUND(c,4)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        cC=c(ix,iy,iz,it)/(Rho(ix,iy,iz,1)+Eps)
        VC=VolC(ix,iy,iz)

!       x-Direction
        cL=c(ix-1,iy,iz,it)/(Rho(ix-1,iy,iz,1)+Eps)
        cR=c(ix+1,iy,iz,it)/(Rho(ix+1,iy,iz,1)+Eps)
        VL=VolC(ix-1,iy,iz)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix+1,iy,iz)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FU(ix,iy,iz)
        FL=FU(ix-1,iy,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cC-cR
        r=(cL-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(uF(ix-1,iy,iz))-uF(ix-1,iy,iz))* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))       
        cFU(ix-1,iy,iz,it)=cFU(ix-1,iy,iz,it)+Flux/(FL*(ABS(uF(ix-1,iy,iz)))+Eps*Eps)
!       Right boundary 
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(uF(ix,iy,iz))+uF(ix,iy,iz))* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        cFU(ix,iy,iz,it)=cFU(ix,iy,iz,it)+Flux/(FR*(ABS(uF(ix,iy,iz)))+Eps*Eps)

!       y-Direction
        cL=c(ix,iy-1,iz,it)/(Rho(ix,iy-1,iz,1)+Eps)
        cR=c(ix,iy+1,iz,it)/(Rho(ix,iy+1,iz,1)+Eps)
        VL=VolC(ix,iy-1,iz)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix,iy+1,iz)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FV(ix,iy,iz)
        FL=FV(ix,iy-1,iz)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cC-cR
        r=(cL-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(vF(ix,iy-1,iz))-vF(ix,iy-1,iz))* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))       
        cFV(ix,iy-1,iz,it)=cFV(ix,iy-1,iz,it)+Flux/(FL*(ABS(vF(ix,iy-1,iz)))+Eps*Eps)
!       Right boundary
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        Flux=Half*(ABS(vF(ix,iy,iz))+vF(ix,iy,iz))* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        cFV(ix,iy,iz,it)=cFV(ix,iy,iz,it)+Flux/(FR*(ABS(vF(ix,iy,iz)))+Eps*Eps)

!       z-Direction
        cL=c(ix,iy,iz-1,it)/(Rho(ix,iy,iz-1,1)+Eps)
        cR=c(ix,iy,iz+1,it)/(Rho(ix,iy,iz+1,1)+Eps)
        VL=VolC(ix,iy,iz-1)
        IF (VL==Zero) THEN
          cL=cC
        END IF
        VR=VolC(ix,iy,iz+1)
        IF (VR==Zero) THEN
          cR=cC
        END IF
        FR=FW(ix,iy,iz)
        FL=FW(ix,iy,iz-1)
        hL=hLF(VL,FL)
        hC=hCF(VC,FL,FR)
        hR=hRF(VR,FR)

!       Left boundary
!       IF (PRESENT(FallF)) THEN
!         RhoLF=CellToFaceVol(RhoR(ix,iy,iz-1,1),RhoR(ix,iy,iz,1) &
!                          ,VL,VC)
!         vFall=FallF(RhoLF)
!       END IF
        k1=k1F(VR,VC,VL)
        k2=k2F(VR,VC,VL)
        cDiff=cC-cR
        r=(cL-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz-1)+CellToFaceVol(Rho(ix,iy,iz-1,1),Rho(ix,iy,iz,1) &
                                        ,VolC(ix,iy,iz-1),VolC(ix,iy,iz)) &
                            *wFall(ix,iy,iz-1)
        Flux=Half*(ABS(wFLoc)-wFLoc)* &
             FL* &
             (cC+FR/(FR+Eps)*phi(r,k1,k2,cL,cC,cR,hC))       
        cFW(ix,iy,iz-1,it)=cFW(ix,iy,iz-1,it)+Flux/(FL*(ABS(wFLoc))+Eps*Eps)
!       Right boundary
!       IF (PRESENT(FallF)) THEN
!         RhoLF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
!                          ,VL,VC)
!         vFall=FallF(RhoLF)
!       END IF
        k1=k1F(VL,VC,VR)
        k2=k2F(VL,VC,VR)
        cDiff=cC-cL
        r=(cR-cC+Eps)/(cDiff+Eps)
        wFLoc=wF(ix,iy,iz)+CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                        ,VolC(ix,iy,iz),VolC(ix,iy,iz+1)) &
                          *wFall(ix,iy,iz)
        Flux=Half*(ABS(wFLoc)+wFLoc)* &
             FR* &
             (cC+FL/(FL+Eps)*phi(r,k1,k2,cR,cC,cL,hC))
        cFW(ix,iy,iz,it)=cFW(ix,iy,iz,it)+Flux/(FR*(ABS(wFLoc))+Eps*Eps)
      END DO
    END DO
  END DO
  IF (TypeW=='ow') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix0,iy,iz))+uF(ix0,iy,iz))* &
              c(ix0,iy,iz,it)/(Rho(ix0,iy,iz,1)+Eps)*FU(ix0,iy,iz)
        cFU(ix0,iy,iz,it)=cFU(ix0,iy,iz,it)+Flux/(FU(ix0,iy,iz)*(ABS(uF(ix0,iy,iz))+Eps)+Eps)
      END DO
    END DO
  END IF
  IF (TypeE=='oe') THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Flux=Half*(ABS(uF(ix1,iy,iz))-uF(ix1,iy,iz))* &
              c(ix1+1,iy,iz,it)/(Rho(ix1+1,iy,iz,1)+Eps)*FU(ix1,iy,iz) 
        cFU(ix1,iy,iz,it)=cFU(ix1,iy,iz,it)+Flux/(FU(ix1,iy,iz)*(ABS(uF(ix1,iy,iz))+Eps)+Eps)
      END DO
    END DO
  END IF
  IF (TypeS=='os') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy0,iz))+vF(ix,iy0,iz))* &
              c(ix,iy0,iz,it)/(Rho(ix,iy0,iz,1)+Eps)*FV(ix,iy0,iz)
        cFV(ix,iy0,iz,it)=cFV(ix,iy0,iz,it)+Flux/(FV(ix,iy0,iz)*(ABS(vF(ix,iy0,iz))+Eps)+Eps)
      END DO
    END DO
  END IF
  IF (TypeN=='on') THEN
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Flux=Half*(ABS(vF(ix,iy1,iz))-vF(ix,iy1,iz))* &
              c(ix,iy1+1,iz,it)/(Rho(ix,iy1+1,iz,1)+Eps)*FV(ix,iy1,iz)
        cFV(ix,iy1,iz,it)=cFV(ix,iy1,iz,it)+Flux/(FV(ix,iy1,iz)*(ABS(vF(ix,iy1,iz))+Eps)+Eps)
      END DO
    END DO
  END IF
  IF (TypeB=='ob') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       IF (PRESENT(FallF)) THEN
!         vFall=FallF(RhoR(ix,iy,iz0+1,1))
!       END IF
        wFLoc=wF(ix,iy,iz0) &
             +Rho(ix,iy,iz0+1,1)*wFall(ix,iy,iz0+1)
        Flux=Half*(ABS(wFLoc)+wFLoc+Eps)* &
              c(ix,iy,iz0,it)/(Rho(ix,iy,iz0+1,1)+Eps)*FW(ix,iy,iz0)
        cFW(ix,iy,iz0,it)=cFW(ix,iy,iz0,it)+Flux/(FW(ix,iy,iz0)*(ABS(wFLoc)+Eps)+Eps)
      END DO
    END DO
  END IF
  IF (TypeT=='ot') THEN
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       IF (PRESENT(FallF)) THEN
!         vFall=FallF(RhoR(ix,iy,iz1,1))
!       END IF
        wFLoc=wF(ix,iy,iz1) &
             +Rho(ix,iy,iz1,1)*wFall(ix,iy,iz1)
        Flux=Half*(ABS(wFLoc)-wFLoc)* &
              c(ix,iy,iz1+1,it)/(Rho(ix,iy,iz1+1,1)+Eps)*FW(ix,iy,iz1)
        cFW(ix,iy,iz1,it)=cFW(ix,iy,iz1,it)+Flux/(FW(ix,iy,iz1)*(ABS(wFLoc)+Eps)+Eps)
      END DO
    END DO
  END IF
  END DO

END SUBROUTINE AdvectionFaceCompute

FUNCTION hLF(VL,FL)
  REAL(RealKind) :: hLF,VL,FL
  hLF=VL/(FL+Eps)
END FUNCTION hLF
FUNCTION hCF(VC,FL,FR)
  REAL(RealKind) :: hCF,VC,FL,FR
  hCF=2.0e0*VC/(FL+FR+Eps)
END FUNCTION hCF
FUNCTION hRF(VR,FR)
  REAL(RealKind) :: hRF,VR,FR
  hRF=VR/(FR+Eps)
END FUNCTION hRF
FUNCTION k1F(hL,hC,hR)
  REAL(RealKind) :: k1F,hL,hC,hR
  k1F=(hC/(hC+hR+Eps))*((hL+hC)/(hL+hC+hR+Eps))
END FUNCTION k1F
FUNCTION k2F(hL,hC,hR)
  REAL(RealKind) :: k2F,hL,hC,hR
  k2F=(hC/(hL+hC+Eps))*(hR/(hL+hC+hR+Eps))
END FUNCTION k2F

FUNCTION PhiLim(r,k1,k2,cL,cC,cR,hC)
  REAL(RealKind) :: PhiLim,r,k1,k2,cL,cC,cR,hC
  REAL(RealKind) :: rad,eps,eta
  REAL(RealKind) :: q,TOL,a,b,c,d,d1,d2
  SELECT CASE(MethAdv) 
    CASE('Koren')
      PhiLim=(MAX(Zero,MIN(r,MIN(k1*r+k2,One))))*(cC-cR)
    CASE('Third')
      PhiLim=(k1*r+k2)*(cC-cR)
    CASE('Upwind')  
      PhiLim=Zero
    CASE('Leer1')
      PhiLim=(Half*(r+r*r)/(One+r*r))*(cC-cR)
    CASE('Leer2')
      PhiLim=(Half*(r+ABS(r))/(One+ABS(r)))*(cC-cR)
    CASE('Central')
      PhiLim=(Half*r)*(cC-cR)
    CASE('Recovery')
      PhiLim=(.25d0*r+.25d0)*(cC-cR)
      !cL=0.25*cL+cC-.25cR=cC+.25*(cL-cR)=cR+0.25*(cL-cC)+(cC-cR))
    CASE('Limo3')
       eps=0.1d0
       eta=(((cC-cR)**2+(cL-cC)**2)/((0.01d0*hC)**2))
       IF (eta<=(1-eps)) THEN
          PhiLim=(k1*r+k2)*(cC-cR)
       ELSE IF (eta>=(1+eps)) THEN
          PhiLim=(MAX(Zero,MIN((k1*r+k2),MAX(-0.25d0*r,MIN(r,(k1*r+k2),0.8d0)))))*(cC-cR)
       ELSE
          PhiLim=0.5d0*((1-(eta-1)/eps)*((k1*r+k2))+(1+(eta-1)/eps)*MAX(Zero,MIN(((k1*r+k2))&
          &,MAX(-0.25d0*r,MIN(r,(k1*r+k2),0.8d0)))))*(cC-cR)
       END IF
    CASE('LDLR')
       IF (hC<=0.0d0) THEN
         PhiLim=0.0d0
       ELSE  
         q=1.4d0
         TOL=0.1d0*(1.0d0/160.0d0)**q
         d1=(cC-cR)/hC
         d2=(cL-cC)/hC
         a=1.0d0+TOL-(2.0d0*ABS(d1)**q*ABS(d2)**q+TOL)/(ABS(d1)**(2.0d0*q)+ABS(d2)**(2.0d0*q)+TOL)
         b=a/(a-1)
         c=((a-1.0d0)*(d2*(1.0d0-b)-d1))/(b-a)
         d=d1-c
         PhiLim=(c*hC*Etap(a)+d*hC*Etap(b))
       END IF  
    CASE DEFAULT 
      PhiLim=(k1*r+k2)*(cC-cR)
  END SELECT  
END FUNCTION PhiLim

FUNCTION Etap(k)
  REAL(RealKind) :: Etap,k
  Etap=-(LOG(ABS(1.0d0-k))+k)/(k**2.0d0)
END FUNCTION Etap

FUNCTION phiUnlim(r,k1,k2,cL,cC,cR,hC)
  REAL(RealKind) :: phiUnlim,r,k1,k2
  REAL(RealKind) :: cL,cC,cR,hC
  phiUnlim=(k1*r+k2)*(cR-cC)
END FUNCTION phiUnlim

FUNCTION phiLimSym(r,k1,k2,cL,cC,cR,hC)
  REAL(RealKind) :: phiLimSym,r,k1,k2
  REAL(RealKind) :: cL,cC,cR,hC
  phiLimSym=(Half*r)*(cR-cC)
END FUNCTION phiLimSym


FUNCTION phiLim1(r,k1,k2)
  REAL(RealKind) :: phiLim1,r,k1,k2
  SELECT CASE(MethAdv) 
    CASE('Koren')
      phiLim1=MAX(Zero,Min(r,MIN(k1*r+k2,One)))
    CASE('Third')
      phiLim1=k1*r+k2
    CASE('Upwind')  
      phiLim1=zero
    CASE('Leer1')
      philim1=Half*(r+r*r)/(One+r*r)
    CASE('Leer2')
      philim1=Half*(r+ABS(r))/(One+ABS(r))
    CASE('Central')
      philim1=Half*r
    CASE DEFAULT
      phiLim1=zero
  END SELECT  
END FUNCTION phiLim1

FUNCTION phiLim2(r,k1,k2)
  REAL(RealKind) :: phiLim2,r,k1,k2
  philim2=Half*r
END FUNCTION phiLim2

SUBROUTINE JacAdvectionCompute(FallF)

  INTERFACE
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF

  INTEGER :: i,ix,iy,iz
  INTEGER :: in,jx,jy,jz
  REAL(RealKind) :: VL,VR
  REAL(RealKind) :: F
  REAL(RealKind) :: RhoLF,vFallRhoL
  REAL(RealKind) :: RhoCoarse,RhoFine,VolCoarse,VolFine
  INTEGER :: nxP2,nyP2,nzP2


  nxP2=ix1-ix0+2
  nyP2=iy1-iy0+2
  nzP2=iz1-iz0+2
  vFall=One

  IF (JacTransportX) THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1-1
          VL=VolC(ix,iy,iz)+Eps
          VR=VolC(ix+1,iy,iz)+Eps
          F=FU(ix,iy,iz)
!         uF > 0
          i=Index(ix,iy,iz)
          AT%Val(i,3)=AT%Val(i,3)+Half*(ABS(uF(ix,iy,iz))+uF(ix,iy,iz)) &
                     *F/VR/(Rho(ix,iy,iz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix,iy,iz))+uF(ix,iy,iz)) &
                     *F/VL/(Rho(ix,iy,iz,1)+Eps)
!         uF < 0
          i=Index(ix+1,iy,iz)
          AT%Val(i,5)=AT%Val(i,5)+Half*(ABS(uF(ix,iy,iz))-uF(ix,iy,iz)) &
                     *F/VL/(Rho(ix+1,iy,iz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix,iy,iz))-uF(ix,iy,iz)) &
                     *F/VR/(Rho(ix+1,iy,iz,1)+Eps)
        END DO
      END DO
    END DO
  END IF
  IF (JacTransportY) THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1-1
        DO ix=ix0+1,ix1
          VL=VolC(ix,iy,iz)+Eps
          VR=VolC(ix,iy+1,iz)+Eps
          F=FV(ix,iy,iz)
!         vF > 0
          i=Index(ix,iy,iz)
          AT%Val(i,2)=AT%Val(i,2)+Half*(ABS(vF(ix,iy,iz))+vF(ix,iy,iz)) &
                     *F/VR/(Rho(ix,iy,iz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(ix,iy,iz))+vF(ix,iy,iz)) &
                     *F/VL/(Rho(ix,iy,iz,1)+Eps)
!         vF < 0
          i=Index(ix,iy+1,iz)
          AT%Val(i,6)=AT%Val(i,6)+Half*(ABS(vF(ix,iy,iz))-vF(ix,iy,iz)) &
                     *F/VL/(Rho(ix,iy+1,iz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(ix,iy,iz))-vF(ix,iy,iz)) &
                     *F/VR/(Rho(ix,iy+1,iz,1)+Eps)
        END DO
      END DO
    END DO
  END IF
  IF (JacTransportZ) THEN
    DO iz=iz0+1,iz1-1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          i=Index(ix,iy,iz)
          VL=VolC(ix,iy,iz)+Eps
          VR=VolC(ix,iy,iz+1)+Eps
          F=FW(ix,iy,iz)
!         wF > 0
          i=Index(ix,iy,iz)
          AT%Val(i,1)=AT%Val(i,1)+Half*(ABS(wF(ix,iy,iz))+wF(ix,iy,iz)) &
                     *F/VR/(Rho(ix,iy,iz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(ix,iy,iz))+wF(ix,iy,iz)) &
                     *F/VL/(Rho(ix,iy,iz,1)+Eps)
!         wF < 0
          i=Index(ix,iy,iz+1)
          AT%Val(i,7)=AT%Val(i,7)+Half*(ABS(wF(ix,iy,iz))-wF(ix,iy,iz)) &
                     *F/VL/(Rho(ix,iy,iz+1,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(ix,iy,iz))-wF(ix,iy,iz)) &
                     *F/VR/(Rho(ix,iy,iz+1,1)+Eps)
          IF (PRESENT(FallF)) THEN
            RhoLF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                            ,VL,VR)
            vFallRhoL=FallF(RhoLF)
            AFallRhoL%Val(i,2)=AFallRhoL%Val(i,2)-vFallRhoL*F/VL/(Rho(ix,iy,iz+1,1)+Eps)
            AFallRhoL%Val(i,1)=AFallRhoL%Val(i,1)+vFallRhoL*F/VR/(Rho(ix,iy,iz+1,1)+Eps)
          END IF
          AFall%Val(i,2)=AFall%Val(i,2)-F/VL/(Rho(ix,iy,iz+1,1)+Eps)
          AFall%Val(i,1)=AFall%Val(i,1)+F/VR/(Rho(ix,iy,iz+1,1)+Eps)
        END DO
      END DO
    END DO
  END IF

  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (JacTransportX) THEN
    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            RhoFine= &
               SUM(Rho(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* & 
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            RhoCoarse= &
               SUM(Rho(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            DO iz=jz,jz+IncrZ-1
              DO iy=jy,jy+IncrY-1
                VL=VolC(ix0,iy,iz)+Eps
                VR=VolC(ix0+1,iy,iz)+Eps
                F=FU(ix0,iy,iz)
!               uF > 0
                i=Index(ix0,iy,iz)
                AT%Val(i,3)=AT%Val(i,3)+Half*(ABS(uF(ix0,iy,iz))+uF(ix0,iy,iz)) &
                   *F/VR/(RhoCoarse+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix0,iy,iz))+uF(ix0,iy,iz)) &
                           *F/VL/(RhoCoarse+Eps)
!               uF < 0
                i=Index(ix0+1,iy,iz)
                AT%Val(i,5)=AT%Val(i,5)+Half*(ABS(uF(ix0,iy,iz))-uF(ix0,iy,iz)) &
                           *F/VL/(RhoFine+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix0,iy,iz))-uF(ix0,iy,iz)) &
                           *F/VR/(RhoFine+Eps)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            VL=VolC(ix0,jy,jz)+Eps
            VR=VolC(ix0+1,jy,jz)+Eps
            F=FU(ix0,jy,jz)
!           uF > 0
            i=Index(ix0,jy,jz)
            AT%Val(i,3)=AT%Val(i,3)+Half*(ABS(uF(ix0,jy,jz))+uF(ix0,jy,jz)) &
                       *F/VR/(Rho(ix0,jy,jz,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix0,jy,jz))+uF(ix0,jy,jz)) &
                       *F/VL/(Rho(ix0,jy,jz,1)+Eps)
!           uF < 0
            i=Index(ix0+1,jy,jz)
            AT%Val(i,5)=AT%Val(i,5)+Half*(ABS(uF(ix0,jy,jz))-uF(ix0,jy,jz)) &
                       *F/VL/(Rho(ix0+1,jy,jz,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix0,jy,jz))-uF(ix0,jy,jz)) &
                       *F/VR/(Rho(ix0+1,jy,jz,1)+Eps)
          END DO
        END DO
      END IF  
    ELSE IF (Nachbars(in)%nType=='ow') THEN
      DO jz=jz0+1,jz1
        DO jy=jy0+1,jy1
          VL=VolC(ix0,jy,jz)+Eps
          VR=VolC(ix0+1,jy,jz)+Eps
          F=FU(ix0,jy,jz)
!         uF > 0
!         i=Index(ix0,jy,jz)
!         AT%Val(i,3)=AT%Val(i,3)+Half*(ABS(uF(ix0,jy,jz))+uF(ix0,jy,jz)) &
!                    *F/VR/(Rho(ix0,jy,jz,1)+Eps)
!         AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix0,jy,jz))+uF(ix0,jy,jz)) &
!                    *F/VL/(Rho(ix0,jy,jz,1)+Eps)
!         uF < 0
          i=Index(ix0+1,jy,jz)
!         AT%Val(i,5)=AT%Val(i,5)+Half*(ABS(uF(ix0,jy,jz))-uF(ix0,jy,jz)) &
!                    *F/VL/(Rho(ix0+1,jy,jz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix0,jy,jz))-uF(ix0,jy,jz)) &
                     *F/VR/(Rho(ix0+1,jy,jz,1)+Eps)
        END DO
      END DO
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jy=jy0+1,jy1,IncrY
            F=SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            RhoFine= &
               SUM(Rho(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* & 
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            RhoCoarse= &
               SUM(Rho(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            DO iz=jz,jz+IncrZ-1
              DO iy=jy,jy+IncrY-1
                VL=VolC(ix1,iy,iz)+Eps
                VR=VolC(ix1+1,iy,iz)+Eps
                F=FU(ix1,iy,iz)
!               uF > 0
                i=Index(ix1,iy,iz)
                AT%Val(i,3)=AT%Val(i,3)+Half*(ABS(uF(ix1,iy,iz))+uF(ix1,iy,iz)) &
                   *F/VR/(RhoFine+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix1,iy,iz))+uF(ix1,iy,iz)) &
                           *F/VL/(RhoFine+Eps)
!               uF < 0
                i=Index(ix1+1,iy,iz)
                AT%Val(i,5)=AT%Val(i,5)+Half*(ABS(uF(ix1,iy,iz))-uF(ix1,iy,iz)) &
                           *F/VL/(RhoCoarse+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix1,iy,iz))-uF(ix1,iy,iz)) &
                           *F/VR/(RhoCoarse+Eps)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jy=jy0+1,jy1
            VL=VolC(ix1,jy,jz)+Eps
            VR=VolC(ix1+1,jy,jz)+Eps
            F=FU(ix1,jy,jz)
!           uF > 0
            i=Index(ix1,jy,jz)
            AT%Val(i,3)=AT%Val(i,3)+Half*(ABS(uF(ix1,jy,jz))+uF(ix1,jy,jz)) &
                       *F/VR/(Rho(ix1,jy,jz,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix1,jy,jz))+uF(ix1,jy,jz)) &
                       *F/VL/(Rho(ix1,jy,jz,1)+Eps)
!           uF < 0
            i=Index(ix1+1,jy,jz)
            AT%Val(i,5)=AT%Val(i,5)+Half*(ABS(uF(ix1,jy,jz))-uF(ix1,jy,jz)) &
                       *F/VL/(Rho(ix1+1,jy,jz,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix1,jy,jz))-uF(ix1,jy,jz)) &
                       *F/VR/(Rho(ix1+1,jy,jz,1)+Eps)
          END DO
        END DO
      END IF
    ELSE IF (Nachbars(in)%nType=='oe') THEN
      DO jz=jz0+1,jz1
        DO jy=jy0+1,jy1
          VL=VolC(ix1,jy,jz)+Eps
          VR=VolC(ix1+1,jy,jz)+Eps
          F=FU(ix1,jy,jz)
!         uF > 0
          i=Index(ix1,jy,jz)
!         AT%Val(i,3)=AT%Val(i,3)+Half*(ABS(uF(ix1,jy,jz))+uF(ix1,jy,jz)) &
!                    *F/VR/(Rho(ix1,jy,jz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix1,jy,jz))+uF(ix1,jy,jz)) &
                     *F/VL/(Rho(ix1,jy,jz,1)+Eps)
!         uF < 0
!         i=Index(ix1+1,jy,jz)
!         AT%Val(i,5)=AT%Val(i,5)+Half*(ABS(uF(ix1,jy,jz))-uF(ix1,jy,jz)) &
!                    *F/VL/(Rho(ix1+1,jy,jz,1)+Eps)
!         AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(uF(ix1,jy,jz))-uF(ix1,jy,jz)) &
!                    *F/VR/(Rho(ix1+1,jy,jz,1)+Eps)
        END DO
      END DO
    END IF
    END IF
    IF (JacTransportY) THEN
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            RhoFine= &
               SUM(Rho(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            RhoCoarse= &
               SUM(Rho(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            DO iz=jz,jz+IncrZ-1
              DO ix=jx,jx+IncrX-1
                VL=VolC(ix,iy0,iz)+Eps
                VR=VolC(ix,iy0+1,iz)+Eps
                F=FV(ix,iy0,iz)
!               vF > 0
                i=Index(ix,iy0,iz)
                AT%Val(i,2)=AT%Val(i,2)+Half*(ABS(vF(ix,iy0,iz))+vF(ix,iy0,iz)) &
                           *F/VR/(RhoCoarse+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(ix,iy0,iz))+vF(ix,iy0,iz)) &
                           *F/VL/(RhoCoarse+Eps)
!               vF < 0
                i=Index(ix,iy0+1,iz)
                AT%Val(i,6)=AT%Val(i,6)+Half*(ABS(vF(ix,iy0,iz))-vF(ix,iy0,iz)) &
                           *F/VL/(RhoFine+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(ix,iy0,iz))-vF(ix,iy0,iz)) &
                           *F/VR/(RhoFine+Eps)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            VL=VolC(jx,iy0,jz)+Eps
            VR=VolC(jx,iy0+1,jz)+Eps
            F=FV(jx,iy0,jz)
!           vF > 0
            i=Index(jx,iy0,jz)
            AT%Val(i,2)=AT%Val(i,2)+Half*(ABS(vF(jx,iy0,jz))+vF(jx,iy0,jz)) &
                       *F/VR/(Rho(jx,iy0,jz,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(jx,iy0,jz))+vF(jx,iy0,jz)) &
                       *F/VL/(Rho(jx,iy0,jz,1)+Eps)
!           vF < 0
            i=Index(jx,iy0+1,jz)
            AT%Val(i,6)=AT%Val(i,6)+Half*(ABS(vF(jx,iy0,jz))-vF(jx,iy0,jz)) &
                       *F/VL/(Rho(jx,iy0+1,jz,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(jx,iy0,jz))-vF(jx,iy0,jz)) &
                       *F/VR/(Rho(jx,iy0+1,jz,1)+Eps)
          END DO
        END DO
      END IF
    ELSE IF (Nachbars(in)%nType=='os') THEN
      DO jz=jz0+1,jz1
        DO jx=jx0+1,jx1
          VL=VolC(jx,iy0,jz)+Eps
          VR=VolC(jx,iy0+1,jz)+Eps
          F=FV(jx,iy0,jz)
!         vF > 0
!         i=Index(jx,iy0,jz)
!         AT%Val(i,2)=AT%Val(i,2)+Half*(ABS(vF(jx,iy0,jz))+vF(jx,iy0,jz)) &
!                    *F/VR/(Rho(jx,iy0,jz,1)+Eps)
!         AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(jx,iy0,jz))+vF(jx,iy0,jz)) &
!                    *F/VL/(Rho(jx,iy0,jz,1)+Eps)
!         vF < 0
          i=Index(jx,iy0+1,jz)
!         AT%Val(i,6)=AT%Val(i,6)+Half*(ABS(vF(jx,iy0,jz))-vF(jx,iy0,jz)) &
!                    *F/VL/(Rho(jx,iy0+1,jz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(jx,iy0,jz))-vF(jx,iy0,jz)) &
                     *F/VR/(Rho(jx,iy0+1,jz,1)+Eps)
        END DO
      END DO
    END IF

    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
      IF (Refine>RefineNachbar) THEN
        DO jz=jz0+1,jz1,IncrZ
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            RhoFine= &
               SUM(Rho(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            RhoCoarse= &
               SUM(Rho(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))/ &
               (VolCoarse+Eps)
            DO iz=jz,jz+IncrZ-1
              DO ix=jx,jx+IncrX-1
                VL=VolC(ix,iy1,iz)+Eps
                VR=VolC(ix,iy1+1,iz)+Eps
                F=FV(ix,iy1,iz)
!               vF > 0
                i=Index(ix,iy1,iz)
                AT%Val(i,2)=AT%Val(i,2)+Half*(ABS(vF(ix,iy1,iz))+vF(ix,iy1,iz)) &
                           *F/VR/(RhoFine+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(ix,iy1,iz))+vF(ix,iy1,iz)) &
                           *F/VL/(RhoFine+Eps)
!               vF < 0
                i=Index(ix,iy1+1,iz)
                AT%Val(i,6)=AT%Val(i,6)+Half*(ABS(vF(ix,iy1,iz))-vF(ix,iy1,iz)) &
                           *F/VL/(RhoCoarse+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(ix,iy1,iz))-vF(ix,iy1,iz)) &
                           *F/VR/(RhoCoarse+Eps)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            VL=VolC(jx,iy1,jz)+Eps
            VR=VolC(jx,iy1+1,jz)+Eps
            F=FV(jx,iy1,jz)
!           vF > 0
            i=Index(jx,iy1,jz)
            AT%Val(i,2)=AT%Val(i,2)+Half*(ABS(vF(jx,iy1,jz))+vF(jx,iy1,jz)) &
                       *F/VR/(Rho(jx,iy1,jz,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(jx,iy1,jz))+vF(jx,iy1,jz)) &
                       *F/VL/(Rho(jx,iy1,jz,1)+Eps)
!           vF < 0
            i=Index(jx,iy1+1,jz)
            AT%Val(i,6)=AT%Val(i,6)+Half*(ABS(vF(jx,iy1,jz))-vF(jx,iy1,jz)) &
                       *F/VL/(Rho(jx,iy1+1,jz,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(jx,iy1,jz))-vF(jx,iy1,jz)) &
                       *F/VR/(Rho(jx,iy1+1,jz,1)+Eps)
          END DO
        END DO
      END IF
    ELSE IF (Nachbars(in)%nType=='on') THEN
      DO jz=jz0+1,jz1
        DO jx=jx0+1,jx1
          VL=VolC(jx,iy1,jz)+Eps
          VR=VolC(jx,iy1+1,jz)+Eps
          F=FV(jx,iy1,jz)
!           vF > 0
          i=Index(jx,iy1,jz)
!         AT%Val(i,2)=AT%Val(i,2)+Half*(ABS(vF(jx,iy1,jz))+vF(jx,iy1,jz)) &
!                    *F/VR/(Rho(jx,iy1,jz,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(jx,iy1,jz))+vF(jx,iy1,jz)) &
                     *F/VL/(Rho(jx,iy1,jz,1)+Eps)
!           vF < 0
!         i=Index(jx,iy1+1,jz)
!         AT%Val(i,6)=AT%Val(i,6)+Half*(ABS(vF(jx,iy1,jz))-vF(jx,iy1,jz)) &
!                    *F/VL/(Rho(jx,iy1+1,jz,1)+Eps)
!         AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(vF(jx,iy1,jz))-vF(jx,iy1,jz)) &
!                    *F/VR/(Rho(jx,iy1+1,jz,1)+Eps)
        END DO
      END DO
    END IF
    END IF
    IF (JacTransportZ) THEN
    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            RhoFine= &
               SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps
            RhoCoarse= &
               SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))/ &
               (VolCoarse+Eps)
            DO iy=jy,jy+IncrY-1
              DO ix=jx,jx+IncrX-1
                VL=VolC(ix,iy,iz0)+Eps
                VR=VolC(ix,iy,iz0+1)+Eps
                F=FW(ix,iy,iz0)
!               wF > 0
                i=Index(ix,iy,iz0)
                AT%Val(i,1)=AT%Val(i,1)+Half*(ABS(wF(ix,iy,iz0))+wF(ix,iy,iz0)) &
                           *F/VR/(RhoCoarse+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(ix,iy,iz0))+wF(ix,iy,iz0)) &
                           *F/VL/(RhoCoarse+Eps)
!               wF < 0
                i=Index(ix,iy,iz0+1)
                AT%Val(i,7)=AT%Val(i,7)+Half*(ABS(wF(ix,iy,iz0))-wF(ix,iy,iz0)) &
                           *F/VL/(RhoFine+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(ix,iy,iz0))-wF(ix,iy,iz0)) &
                           *F/VR/(RhoFine+Eps)
                IF (PRESENT(FallF)) THEN
                  RhoLF=CellToFaceVol(RhoR(ix,iy,iz0,1),RhoR(ix,iy,iz0+1,1) &
                                  ,VL,VR)
                  vFallRhoL=FallF(RhoLF)
                  AFallRhoL%Val(i,2)=AFallRhoL%Val(i,2)-vFallRhoL*F/VL/(RhoCoarse+Eps)
                  AFallRhoL%Val(i,1)=AFallRhoL%Val(i,1)+vFallRhoL*F/VR/(RhoCoarse+Eps)
                END IF
                AFall%Val(i,2)=AFall%Val(i,2)-F/VL/(RhoCoarse+Eps)
                AFall%Val(i,1)=AFall%Val(i,1)+F/VR/(RhoCoarse+Eps)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            i=Index(jx,jy,iz0)
            VL=VolC(jx,jy,iz0)+Eps
            VR=VolC(jx,jy,iz0+1)+Eps
            F=FW(jx,jy,iz0)
!           wF > 0
            i=Index(jx,jy,iz0)
            AT%Val(i,1)=AT%Val(i,1)+Half*(ABS(wF(jx,jy,iz0))+wF(jx,jy,iz0)) &
                       *F/VR/(Rho(jx,jy,iz0,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(jx,jy,iz0))+wF(jx,jy,iz0)) &
                       *F/VL/(Rho(jx,jy,iz0,1)+Eps)
!           wF < 0
            i=Index(jx,jy,iz0+1)
            AT%Val(i,7)=AT%Val(i,7)+Half*(ABS(wF(jx,jy,iz0))-wF(jx,jy,iz0)) &
                       *F/VL/(Rho(jx,jy,iz0+1,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(jx,jy,iz0))-wF(jx,jy,iz0)) &
                       *F/VR/(Rho(jx,jy,iz0+1,1)+Eps)
            IF (PRESENT(FallF)) THEN 
              RhoLF=CellToFaceVol(RhoR(jx,jy,iz0,1),RhoR(jx,jy,iz0+1,1) &
                              ,VL,VR)
              vFallRhoL=FallF(RhoLF)
              AFallRhoL%Val(i,2)=AFallRhoL%Val(i,2)-vFallRhoL*F/VL/(Rho(jx,jy,iz0+1,1)+Eps)
              AFallRhoL%Val(i,1)=AFallRhoL%Val(i,1)+vFallRhoL*F/VR/(Rho(jx,jy,iz0+1,1)+Eps)
            END IF
            AFall%Val(i,2)=AFall%Val(i,2)-F/VL/(Rho(jx,jy,iz0+1,1)+Eps)
            AFall%Val(i,1)=AFall%Val(i,1)+F/VR/(Rho(jx,jy,iz0+1,1)+Eps)
          END DO
        END DO
      END IF
    ELSE IF (Nachbars(in)%nType=='ob') THEN
      DO jy=jy0+1,jy1
        DO jx=jx0+1,jx1
          i=Index(jx,jy,iz0)
          VL=VolC(jx,jy,iz0)+Eps
          VR=VolC(jx,jy,iz0+1)+Eps
          F=FW(jx,jy,iz0)
!         wF > 0
!         i=Index(jx,jy,iz0)
!         AT%Val(i,1)=AT%Val(i,1)+Half*(ABS(wF(jx,jy,iz0))+wF(jx,jy,iz0)) &
!                    *F/VR/(Rho(jx,jy,iz0,1)+Eps)
!         AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(jx,jy,iz0))+wF(jx,jy,iz0)) &
!                    *F/VL/(Rho(jx,jy,iz0,1)+Eps)
!         wF < 0
          i=Index(jx,jy,iz0+1)
!         AT%Val(i,7)=AT%Val(i,7)+Half*(ABS(wF(jx,jy,iz0))-wF(jx,jy,iz0)) &
!                    *F/VL/(Rho(jx,jy,iz0+1,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(jx,jy,iz0))-wF(jx,jy,iz0)) &
                     *F/VR/(Rho(jx,jy,iz0+1,1)+Eps)
          IF (PRESENT(FallF)) THEN 
            RhoLF=CellToFaceVol(RhoR(jx,jy,iz0,1),RhoR(jx,jy,iz0+1,1) &
                            ,VL,VR)
            vFallRhoL=FallF(RhoLF)
            AFallRhoL%Val(i,1)=AFallRhoL%Val(i,1)+vFallRhoL*F/VR/(Rho(jx,jy,iz0+1,1)+Eps)
          END IF
          AFall%Val(i,1)=AFall%Val(i,1)+F/VR/(Rho(jx,jy,iz0+1,1)+Eps)
        END DO
      END DO
    END IF

    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
      IF (Refine>RefineNachbar) THEN
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            F=SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            RhoFine= &
               SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            RhoCoarse= &
               SUM(Rho(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))/ &
               (VolCoarse+Eps)
            DO iy=jy,jy+IncrY-1
              DO ix=jx,jx+IncrX-1
                i=Index(ix,iy,iz1)
                VL=VolC(ix,iy,iz1)+Eps
                VR=VolC(ix,iy,iz1+1)+Eps
                F=FW(ix,iy,iz1)
!               wF > 0
                i=Index(ix,iy,iz1)
                AT%Val(i,1)=AT%Val(i,1)+Half*(ABS(wF(ix,iy,iz1))+wF(ix,iy,iz1)) &
                           *F/VR/(RhoFine+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(ix,iy,iz1))+wF(ix,iy,iz1)) &
                           *F/VL/(RhoFine+Eps)
!               wF < 0
                i=Index(ix,iy,iz1+1)
                AT%Val(i,7)=AT%Val(i,7)+Half*(ABS(wF(ix,iy,iz1))-wF(ix,iy,iz1)) &
                           *F/VL/(RhoCoarse+Eps)
                AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(ix,iy,iz1))-wF(ix,iy,iz1)) &
                           *F/VR/(RhoCoarse+Eps)
                IF (PRESENT(FallF)) THEN
                  RhoLF=CellToFaceVol(RhoR(ix,iy,iz1,1),RhoR(ix,iy,iz1+1,1) &
                                  ,VL,VR)
                  vFallRhoL=FallF(RhoLF)
                  AFallRhoL%Val(i,2)=AFallRhoL%Val(i,2)-vFallRhoL*F/VL/(RhoFine+Eps)
                  AFallRhoL%Val(i,1)=AFallRhoL%Val(i,1)+vFallRhoL*F/VR/(RhoFine+Eps)
                END IF
                AFall%Val(i,2)=AFall%Val(i,2)-F/VL/(RhoFine+Eps)
                AFall%Val(i,1)=AFall%Val(i,1)+F/VR/(RhoFine+Eps)
              END DO
            END DO
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            i=Index(jx,jy,iz1)
            VL=VolC(jx,jy,iz1)+Eps
            VR=VolC(jx,jy,iz1+1)+Eps
            F=FW(jx,jy,iz1)
!           wF > 0
            i=Index(jx,jy,iz1)
            AT%Val(i,1)=AT%Val(i,1)+Half*(ABS(wF(jx,jy,iz1))+wF(jx,jy,iz1)) &
                       *F/VR/(Rho(jx,jy,iz1,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(jx,jy,iz1))+wF(jx,jy,iz1)) &
                       *F/VL/(Rho(jx,jy,iz1,1)+Eps)
!           wF < 0
            i=Index(jx,jy,iz1+1)
            AT%Val(i,7)=AT%Val(i,7)+Half*(ABS(wF(jx,jy,iz1))-wF(jx,jy,iz1)) &
                       *F/VL/(Rho(jx,jy,iz1+1,1)+Eps)
            AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(jx,jy,iz1))-wF(jx,jy,iz1)) &
                       *F/VR/(Rho(jx,jy,iz1+1,1)+Eps)
            IF (PRESENT(FallF)) THEN 
              RhoLF=CellToFaceVol(RhoR(jx,jy,iz1,1),RhoR(jx,jy,iz1+1,1) &
                              ,VL,VR)
              vFallRhoL=FallF(RhoLF)
              AFallRhoL%Val(i,2)=AFallRhoL%Val(i,2)-vFallRhoL*F/VL/(Rho(jx,jy,iz1+1,1)+Eps)
              AFallRhoL%Val(i,1)=AFallRhoL%Val(i,1)+vFallRhoL*F/VR/(Rho(jx,jy,iz1+1,1)+Eps)
            END IF
            AFall%Val(i,2)=AFall%Val(i,2)-F/VL/(Rho(jx,jy,iz1+1,1)+Eps)
            AFall%Val(i,1)=AFall%Val(i,1)+F/VR/(Rho(jx,jy,iz1+1,1)+Eps)
          END DO
        END DO
      END IF
    ELSE IF (Nachbars(in)%nType=='ot') THEN
      DO jy=jy0+1,jy1
        DO jx=jx0+1,jx1
          i=Index(jx,jy,iz1)
          VL=VolC(jx,jy,iz1)+Eps
          VR=VolC(jx,jy,iz1+1)+Eps
          F=FW(jx,jy,iz1)
!           wF > 0
          i=Index(jx,jy,iz1)
!         AT%Val(i,1)=AT%Val(i,1)+Half*(ABS(wF(jx,jy,iz1))+wF(jx,jy,iz1)) &
!                    *F/VR/(Rho(jx,jy,iz1,1)+Eps)
          AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(jx,jy,iz1))+wF(jx,jy,iz1)) &
                     *F/VL/(Rho(jx,jy,iz1,1)+Eps)
!           wF < 0
!         i=Index(jx,jy,iz1+1)
!         AT%Val(i,7)=AT%Val(i,7)+Half*(ABS(wF(jx,jy,iz1))-wF(jx,jy,iz1)) &
!                    *F/VL/(Rho(jx,jy,iz1+1,1)+Eps)
!         AT%Val(i,4)=AT%Val(i,4)-Half*(ABS(wF(jx,jy,iz1))-wF(jx,jy,iz1)) &
!                    *F/VR/(Rho(jx,jy,iz1+1,1)+Eps)
        END DO
      END DO
    END IF
    END IF
  END DO

CONTAINS
FUNCTION Index(ix,iy,iz)
  INTEGER :: Index,ix,iy,iz
  Index=ix-ix0+1+nxP2*(iy-iy0)+nxP2*nyP2*(iz-iz0)  
END FUNCTION Index
END SUBROUTINE JacAdvectionCompute

SUBROUTINE AdvectionUpdate(alpha)

  REAL(RealKind) :: alpha
  INTEGER :: it

  DO it=LBOUND(c,4),UBOUND(c,4)
    c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,it)=c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,it)+ &
           alpha*f(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,it)/ &
           (VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1) &
            +Eps)  
  END DO

END SUBROUTINE AdvectionUpdate

SUBROUTINE AdvectionScale

  INTEGER :: ix,iy,iz,it

  DO it=LBOUND(f,4),UBOUNd(f,4)
    f(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,it)=                                 &
           f(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,it)/ &
           (VolC(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1) &
            +Eps)  
  END DO

END SUBROUTINE AdvectionScale

SUBROUTINE DivCompute

  INTEGER :: ix,iy,iz

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoRhs(ix,iy,iz,1)=RhoRhs(ix,iy,iz,1)- &
                (uF(ix,iy,iz)*FU(ix,iy,iz)-uF(ix-1,iy,iz)*FU(ix-1,iy,iz) &
                +vF(ix,iy,iz)*FV(ix,iy,iz)-vF(ix,iy-1,iz)*FV(ix,iy-1,iz) &
                +wF(ix,iy,iz)*FW(ix,iy,iz)-wF(ix,iy,iz-1)*FW(ix,iy,iz-1))
      END DO
    END DO
  END DO

END SUBROUTINE DivCompute

SUBROUTINE DivPreCompute

  INTEGER :: ix,iy,iz

  REAL(RealKind) :: DivLoc(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)
  REAL(RealKind) :: PreFac,FluxLoc

  DivLoc=Zero

  ix=ix0
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      FluxLoc=uf(ix,iy,iz)*FU(ix,iy,iz)*(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)) &
              /(VolC(ix+1,iy,iz)*Rho(ix+1,iy,iz,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
      DivLoc(ix+1,iy,iz)=DivLoc(ix+1,iy,iz)-FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        FluxLoc=uf(ix,iy,iz)*FU(ix,iy,iz)*(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)) &
                /(VolC(ix+1,iy,iz)*Rho(ix+1,iy,iz,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        DivLoc(ix+1,iy,iz)=DivLoc(ix+1,iy,iz)-FluxLoc
        DivLoc(ix,iy,iz)=DivLoc(ix,iy,iz)+FluxLoc
      END DO
    END DO
  END DO
  ix=ix1
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      FluxLoc=uf(ix,iy,iz)*FU(ix,iy,iz)*(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)) &
              /(VolC(ix+1,iy,iz)*Rho(ix+1,iy,iz,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
      DivLoc(ix,iy,iz)=DivLoc(ix,iy,iz)+FluxLoc
    END DO
  END DO
 
  iy=iy0
  DO iz=iz0+1,iz1
    DO ix=ix0+1,ix1
      FluxLoc=vf(ix,iy,iz)*FV(ix,iy,iz)*(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)) &
              /(VolC(ix,iy+1,iz)*Rho(ix,iy+1,iz,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
      DivLoc(ix,iy+1,iz)=DivLoc(ix,iy+1,iz)-FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        FluxLoc=vf(ix,iy,iz)*FV(ix,iy,iz)*(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)) &
                /(VolC(ix,iy+1,iz)*Rho(ix,iy+1,iz,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        DivLoc(ix,iy+1,iz)=DivLoc(ix,iy+1,iz)-FluxLoc
        DivLoc(ix,iy,iz)=DivLoc(ix,iy,iz)+FluxLoc
      END DO
    END DO
  END DO
  iy=iy1
  DO iz=iz0+1,iz1
    DO ix=ix0+1,ix1
      FluxLoc=vf(ix,iy,iz)*FV(ix,iy,iz)*(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)) &
              /(VolC(ix,iy+1,iz)*Rho(ix,iy+1,iz,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
      DivLoc(ix,iy,iz)=DivLoc(ix,iy,iz)+FluxLoc
    END DO
  END DO
 
  iz=iz0
  DO iy=iy0+1,iy1
    DO ix=ix0+1,ix1
      FluxLoc=wf(ix,iy,iz)*FW(ix,iy,iz)*(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)) & 
              /(VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
      DivLoc(ix,iy,iz+1)=DivLoc(ix,iy,iz+1)-FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        FluxLoc=wf(ix,iy,iz)*FW(ix,iy,iz)*(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)) & 
                /(VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        DivLoc(ix,iy,iz+1)=DivLoc(ix,iy,iz+1)-FluxLoc
        DivLoc(ix,iy,iz)=DivLoc(ix,iy,iz)+FluxLoc
      END DO
    END DO
  END DO
  iz=iz1
  DO iy=iy0+1,iy1
    DO ix=ix0+1,ix1
      FluxLoc=wf(ix,iy,iz)*FW(ix,iy,iz)*(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)) & 
              /(VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)+VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
      DivLoc(ix,iy,iz)=DivLoc(ix,iy,iz)+FluxLoc
    END DO
  END DO
 
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        PreFac=Sound(ix,iy,iz,1)
        thRhs(ix,iy,iz,1)=thRhs(ix,iy,iz,1)-DivLoc(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)*PreFac  
      END DO
    END DO
  END DO

END SUBROUTINE DivPreCompute

SUBROUTINE DivEnCompute

  INTEGER :: ix,iy,iz

  REAL(RealKind) :: FluxLoc,Fac

  ix=ix0
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Fac=PreFace(p(ix,iy,iz,1),p(ix+1,iy,iz,1),Rho(ix,iy,iz,1)     &
                 ,Rho(ix+1,iy,iz,1),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
      FluxLoc=uFace(uf(ix,iy,iz),Rho(ix,iy,iz,1),Rho(ix+1,iy,iz,1) &
                   ,VolC(ix,iy,iz),VolC(ix+1,iy,iz)) &
             *FU(ix,iy,iz)*Fac
      f(ix+1,iy,iz,1)=f(ix+1,iy,iz,1)+FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        Fac=PreFace(p(ix,iy,iz,1),p(ix+1,iy,iz,1),Rho(ix,iy,iz,1)     &
                   ,Rho(ix+1,iy,iz,1),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
        FluxLoc=uFace(uf(ix,iy,iz),Rho(ix,iy,iz,1),Rho(ix+1,iy,iz,1) &
                     ,VolC(ix,iy,iz),VolC(ix+1,iy,iz)) &
               *FU(ix,iy,iz)*Fac
        f(ix+1,iy,iz,1)=f(ix+1,iy,iz,1)+FluxLoc
        f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
      END DO
    END DO
  END DO
  ix=ix1
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Fac=PreFace(p(ix,iy,iz,1),p(ix+1,iy,iz,1),Rho(ix,iy,iz,1)     &
                 ,Rho(ix+1,iy,iz,1),VolC(ix,iy,iz),VolC(ix+1,iy,iz))
      FluxLoc=uFace(uf(ix,iy,iz),Rho(ix,iy,iz,1),Rho(ix+1,iy,iz,1) &
                   ,VolC(ix,iy,iz),VolC(ix+1,iy,iz)) &
             *FU(ix,iy,iz)*Fac
      f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
    END DO
  END DO
 
  iy=iy0
  DO iz=iz0+1,iz1
    DO ix=ix0+1,ix1
      Fac=PreFace(p(ix,iy,iz,1),p(ix,iy+1,iz,1),Rho(ix,iy,iz,1)     &
                 ,Rho(ix,iy+1,iz,1),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
      FluxLoc=vf(ix,iy,iz)*FV(ix,iy,iz)*Fac
      f(ix,iy+1,iz,1)=f(ix,iy+1,iz,1)+FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        Fac=PreFace(p(ix,iy,iz,1),p(ix,iy+1,iz,1),Rho(ix,iy,iz,1)     &
                   ,Rho(ix,iy+1,iz,1),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
        FluxLoc=vf(ix,iy,iz)*FV(ix,iy,iz)*Fac
        f(ix,iy+1,iz,1)=f(ix,iy+1,iz,1)+FluxLoc
        f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
      END DO
    END DO
  END DO
  iy=iy1
  DO iz=iz0+1,iz1
    DO ix=ix0+1,ix1
      Fac=PreFace(p(ix,iy,iz,1),p(ix,iy+1,iz,1),Rho(ix,iy,iz,1)     &
                 ,Rho(ix,iy+1,iz,1),VolC(ix,iy,iz),VolC(ix,iy+1,iz))
      FluxLoc=vf(ix,iy,iz)*FV(ix,iy,iz)*Fac
      f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
    END DO
  END DO
 
  iz=iz0
  DO iy=iy0+1,iy1
    DO ix=ix0+1,ix1
      Fac=PreFace(p(ix,iy,iz,1),p(ix,iy,iz+1,1),Rho(ix,iy,iz,1)     &
                 ,Rho(ix,iy,iz+1,1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
      FluxLoc=wf(ix,iy,iz)*FW(ix,iy,iz)*Fac
      f(ix,iy,iz+1,1)=f(ix,iy,iz+1,1)+FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Fac=PreFace(p(ix,iy,iz,1),p(ix,iy,iz+1,1),Rho(ix,iy,iz,1)     &
                   ,Rho(ix,iy,iz+1,1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
        FluxLoc=wf(ix,iy,iz)*FW(ix,iy,iz)*Fac
        f(ix,iy,iz+1,1)=f(ix,iy,iz+1,1)+FluxLoc
        f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
      END DO
    END DO
  END DO
  iz=iz1
  DO iy=iy0+1,iy1
    DO ix=ix0+1,ix1
      Fac=PreFace(p(ix,iy,iz,1),p(ix,iy,iz+1,1),Rho(ix,iy,iz,1)     &
                 ,Rho(ix,iy,iz+1,1),VolC(ix,iy,iz),VolC(ix,iy,iz+1))
      FluxLoc=wf(ix,iy,iz)*FW(ix,iy,iz)*Fac
      f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
    END DO
  END DO
 
END SUBROUTINE DivEnCompute

SUBROUTINE DivWeightCompute

  INTEGER :: ix,iy,iz

  REAL(RealKind) :: FluxLoc,Fac

  ix=ix0
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Fac=WFace(Weight(ix,iy,iz,1),Weight(ix+1,iy,iz,1)     &
               ,Rho(ix,iy,iz,1),Rho(ix+1,iy,iz,1)     &
               ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
      FluxLoc=uF(ix,iy,iz)*FU(Ix,iy,iz)*Fac
      f(ix+1,iy,iz,1)=f(ix+1,iy,iz,1)+FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        Fac=WFace(Weight(ix,iy,iz,1),Weight(ix+1,iy,iz,1)     &
                 ,Rho(ix,iy,iz,1),Rho(ix+1,iy,iz,1)     &
                 ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
        FluxLoc=uF(ix,iy,iz)*FU(Ix,iy,iz)*Fac
        f(ix+1,iy,iz,1)=f(ix+1,iy,iz,1)+FluxLoc
        f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
      END DO
    END DO
  END DO
  ix=ix1
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      Fac=WFace(Weight(ix,iy,iz,1),Weight(ix+1,iy,iz,1)     &
               ,Rho(ix,iy,iz,1),Rho(ix+1,iy,iz,1)     &
               ,VolC(ix,iy,iz),VolC(ix+1,iy,iz))
      FluxLoc=uF(ix,iy,iz)*FU(Ix,iy,iz)*Fac
      f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
    END DO
  END DO
 
  iy=iy0
  DO iz=iz0+1,iz1
    DO ix=ix0+1,ix1
      Fac=WFace(Weight(ix,iy,iz,1),Weight(ix,iy+1,iz,1)     &
               ,Rho(ix,iy,iz,1),Rho(ix,iy+1,iz,1)     &
               ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
      FluxLoc=vF(ix,iy,iz)*FV(ix,iy,iz)*Fac
      f(ix,iy+1,iz,1)=f(ix,iy+1,iz,1)+FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        Fac=WFace(Weight(ix,iy,iz,1),Weight(ix,iy+1,iz,1)     &
                 ,Rho(ix,iy,iz,1),Rho(ix,iy+1,iz,1)     &
                 ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
        FluxLoc=vF(ix,iy,iz)*FV(ix,iy,iz)*Fac
        f(ix,iy+1,iz,1)=f(ix,iy+1,iz,1)+FluxLoc
        f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
      END DO
    END DO
  END DO
  iy=iy1
  DO iz=iz0+1,iz1
    DO ix=ix0+1,ix1
      Fac=WFace(Weight(ix,iy,iz,1),Weight(ix,iy+1,iz,1)     &
               ,Rho(ix,iy,iz,1),Rho(ix,iy+1,iz,1)     &
               ,VolC(ix,iy,iz),VolC(ix,iy+1,iz))
      FluxLoc=vF(ix,iy,iz)*FV(ix,iy,iz)*Fac
      f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
    END DO
  END DO
 
  iz=iz0
  DO iy=iy0+1,iy1
    DO ix=ix0+1,ix1
      Fac=WFace(Weight(ix,iy,iz,1),Weight(ix,iy,iz+1,1)     &
               ,Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1)     &
               ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
      FluxLoc=wF(ix,iy,iz)*FW(ix,iy,iz)*Fac
      f(ix,iy,iz+1,1)=f(ix,iy,iz+1,1)+FluxLoc
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Fac=WFace(Weight(ix,iy,iz,1),Weight(ix,iy,iz+1,1)     &
                 ,Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1)     &
                 ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
        FluxLoc=wF(ix,iy,iz)*FW(ix,iy,iz)*Fac
        f(ix,iy,iz+1,1)=f(ix,iy,iz+1,1)+FluxLoc
        f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
      END DO
    END DO
  END DO
  iz=iz1
  DO iy=iy0+1,iy1
    DO ix=ix0+1,ix1
      Fac=WFace(Weight(ix,iy,iz,1),Weight(ix,iy,iz+1,1)     &
               ,Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1)     &
               ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
      FluxLoc=wF(ix,iy,iz)*FW(ix,iy,iz)*Fac
      f(ix,iy,iz,1)=f(ix,iy,iz,1)-FluxLoc
    END DO
  END DO
 
END SUBROUTINE DivWeightCompute

SUBROUTINE DivFCompute

  INTEGER :: ix,iy,iz

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        DivC(ix,iy,iz,1)= &
                (uF(ix,iy,iz)*FU(ix,iy,iz)-uF(ix-1,iy,iz)*FU(ix-1,iy,iz) &
                +vF(ix,iy,iz)*FV(ix,iy,iz)-vF(ix,iy-1,iz)*FV(ix,iy-1,iz) &
                +wF(ix,iy,iz)*FW(ix,iy,iz)-wF(ix,iy,iz-1)*FW(ix,iy,iz-1))/(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE DivFCompute

SUBROUTINE DivScalarCompute(FallF,Species)

  INTERFACE
    FUNCTION FallF(qc)
      USE Kind_Mod
      REAL(RealKind) :: FallF,qc
    END FUNCTION FallF
  END INTERFACE

  OPTIONAL :: FallF
  CHARACTER(2),OPTIONAL :: Species

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Flux(SIZE(cFU,4)),FacLoc
  REAL(RealKind) :: VL,VR,VC
  REAL(RealKind) :: wFLoc
  REAL(RealKind) :: RhoLF
  REAL(RealKind) :: RhoRF,RhoF,NRF,RhoIF,NIF
  REAL(RealKind) :: RhoSF,NSF
  REAL(RealKind) :: wFall(ix0+1:ix1,iy0+1:iy1,iz0:iz1)

  IF (PRESENT(Species)) THEN
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          SELECT CASE(Species)
          CASE('RhoR')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoRF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                               ,VL,VR)
            NRF=CellToFaceVol(Nrain(ix,iy,iz,1),Nrain(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_rain_q(RhoRF,NRF,RhoF)!+vFall+vSub(zP(iz)) 
          CASE('NR')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoRF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
                               ,VL,VR)
            NRF=CellToFaceVol(Nrain(ix,iy,iz,1),Nrain(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_rain_n(RhoRF,NRF,RhoF)!+vFall+vSub(zP(iz))
          CASE('RhoI')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoIF=CellToFaceVol(RhoI(ix,iy,iz,1),RhoI(ix,iy,iz+1,1) &
                               ,VL,VR)
            NIF=CellToFaceVol(Nice(ix,iy,iz,1),Nice(ix,iy,iz+1,1) &
                               ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                               ,VL,VR)
            wFall(ix,iy,iz)=fall_ice_q(RhoIF,NIF,RhoF)!+vFall+vSub(zP(iz))
          CASE('NI')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoIF=CellToFaceVol(RhoI(ix,iy,iz,1),RhoI(ix,iy,iz+1,1) &
                                ,VL,VR)
            NIF=CellToFaceVol(Nice(ix,iy,iz,1),Nice(ix,iy,iz+1,1) &
                                ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                ,VL,VR)
            wFall(ix,iy,iz)=fall_ice_n(RhoIF,NIF,RhoF)!+vFall+vSub(zP(iz))
          CASE('RhoS')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoSF=CellToFaceVol(RhoS(ix,iy,iz,1),RhoS(ix,iy,iz+1,1) &
                                ,VL,VR)
            NSF=CellToFaceVol(Nsnow(ix,iy,iz,1),Nsnow(ix,iy,iz+1,1) &
                                ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                ,VL,VR)
            wFall(ix,iy,iz)=fall_snow_q(RhoSF,NSF,RhoF)!+vFall+vSub(zP(iz))
          CASE('NS')
            VL=VolC(ix,iy,iz)
            VR=VolC(ix,iy,iz+1)
            RhoSF=CellToFaceVol(RhoS(ix,iy,iz,1),RhoS(ix,iy,iz+1,1) &
                                 ,VL,VR)
            NSF=CellToFaceVol(Nsnow(ix,iy,iz,1),Nsnow(ix,iy,iz+1,1) &
                                 ,VL,VR)
            RhoF=CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                 ,VL,VR)
            wFall(ix,iy,iz)=fall_snow_n(RhoSF,NSF,RhoF)!+vFall+vSub(zP(iz))
          END SELECT                                                                                                                                                   
        END DO
      END DO
    END DO                                                                                                                                                               
  ELSE
    wFall=0.0d0
  END IF

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0,ix1
        Flux=FU(ix,iy,iz)*cFU(ix,iy,iz,:)*uF(ix,iy,iz) 
        f(ix,iy,iz,:)=f(ix,iy,iz,:)-Flux 
        f(ix+1,iy,iz,:)=f(ix+1,iy,iz,:)+Flux 
      END DO
    END DO
  END DO
  
  DO iz=iz0+1,iz1
    DO iy=iy0,iy1
      DO ix=ix0+1,ix1
        Flux=FV(ix,iy,iz)*cFV(ix,iy,iz,:)*vF(ix,iy,iz) 
        f(ix,iy,iz,:)=f(ix,iy,iz,:)-Flux 
        f(ix,iy+1,iz,:)=f(ix,iy+1,iz,:)+Flux 
      END DO
    END DO
  END DO

! IF (PRESENT(FallF)) THEN
!   DO iz=iz0,iz1
!     DO iy=iy0+1,iy1
!       DO ix=ix0+1,ix1
!         RhoLF=CellToFaceVol(RhoR(ix,iy,iz,1),RhoR(ix,iy,iz+1,1) &
!                            ,VolC(ix,iy,iz),VolC(ix,iy,iz+1))
!         vFall=FallF(RhoLF)
!         wFLoc=wF(ix,iy,iz)+CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
!                                         ,VolC(ix,iy,iz),VolC(ix,iy,iz+1)) &
!                           *(vFall+vSub(iz))
!         Flux=FW(ix,iy,iz)*cFW(ix,iy,iz,:)*wFLoc
!         f(ix,iy,iz,:)=f(ix,iy,iz,:)-Flux 
!         f(ix,iy,iz+1,:)=f(ix,iy,iz+1,:)+Flux 
!       END DO
!     END DO
!   END DO
! ELSE  
    DO iz=iz0,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wFLoc=wF(ix,iy,iz)+CellToFaceVol(Rho(ix,iy,iz,1),Rho(ix,iy,iz+1,1) &
                                          ,VolC(ix,iy,iz),VolC(ix,iy,iz+1)) &
                            *wFall(ix,iy,iz)
          Flux=FW(ix,iy,iz)*cFW(ix,iy,iz,:)*wFLoc
          f(ix,iy,iz,:)=f(ix,iy,iz,:)-Flux 
          f(ix,iy,iz+1,:)=f(ix,iy,iz+1,:)+Flux 
        END DO
      END DO
    END DO
! END IF
END SUBROUTINE DivScalarCompute

FUNCTION vSub(iz)

  REAL(RealKind) :: vSub
  INTEGER :: iz

  REAL(RealKind) :: zHeight

  vSub=0.0d0

END FUNCTION vSub

FUNCTION vSub1(iz)

  REAL(RealKind) :: vSub1
  INTEGER :: iz

  REAL(RealKind) :: zHeight

  zHeight=zP(iz)
  SELECT CASE(vSubProfile)
    CASE('BOMEX')
      IF(zHeight<1500.0d0) THEN
        vSub1=(vSubConst/1500.0d0)*zHeight
      ELSE IF(zHeight>=1500.0d0.AND.zHeight<2100.0d0) THEN
        vSub1=vSubConst-vSubConst/(2100.0d0-1500.0d0)*(zHeight-1500.0d0)
      ELSE
        vSub1=Zero
      END IF
    CASE('RICO')
      IF(zHeight<2260.0d0) THEN
        vSub1=(vSubConst/2260.0d0)*zHeight
      ELSE IF(zHeight>=2260.0d0.AND.zHeight<4000.0d0) THEN
        vSub1=vSubConst
      ELSE
        vSub1=Zero
      END IF
    CASE('ISDAC')
      IF(zHeight<825.0d0) THEN
        vSub1=-5.0d-6*zHeight
      ELSE
        vSub1=-0.412510d-2
      END IF
    CASE('FORCING')
      vSub1=wSubs(iz)
    CASE DEFAULT
      vSub1=vSubConst
  END SELECT

END FUNCTION vSub1

SUBROUTINE SubsidenceScalar(Vec,Rhs)

  TYPE(Vec4_T), POINTER :: Vec(:)
  TYPE(Vec4_T), POINTER :: Rhs(:)

  LOGICAL :: DivRho=.TRUE.
  LOGICAL :: RhoTend=.TRUE.

  Rho=>Vec(RhoPos)%c
  RhoRhs=>Rhs(RhoPos)%c

! Applying subsidence source term to all scalars
  IF (thPos>0) THEN
    c=>Vec(ThPos)%c
    f=>Rhs(ThPos)%c
    CALL SubsidenceCompute(DivRho)
  END IF
  IF (RhoVPos>0) THEN
    c=>Vec(RhoVPos)%c
    f=>Rhs(RhoVPos)%c
    CALL SubsidenceCompute(DivRho,RhoTend)
  END IF
  IF (RhoCPos>0) THEN
    c=>Vec(RhoCPos)%c
    f=>Rhs(RhoCPos)%c
    CALL SubsidenceCompute(DivRho,RhoTend)
  END IF
  IF (RhoRPos>0) THEN
    c=>Vec(RhoRPos)%c
    f=>Rhs(RhoRPos)%c
    CALL SubsidenceCompute(DivRho,RhoTend)
  END IF
  IF (RhoIPos>0) THEN
    c=>Vec(RhoIPos)%c
    f=>Rhs(RhoIPos)%c
    CALL SubsidenceCompute(DivRho,RhoTend)
  END IF
  IF (RhoSPos>0) THEN
    c=>Vec(RhoSPos)%c
    f=>Rhs(RhoSPos)%c
    CALL SubsidenceCompute(DivRho,RhoTend)
  END IF
  IF (NvPos>0) THEN
    c=>Vec(NvPos)%c
    f=>Rhs(NvPos)%c
    CALL SubsidenceCompute
  END IF
  IF (NcPos>0) THEN
    c=>Vec(NcPos)%c
    f=>Rhs(NcPos)%c
    CALL SubsidenceCompute
  END IF
  IF (NrPos>0) THEN
    c=>Vec(NrPos)%c
    f=>Rhs(NrPos)%c
    CALL SubsidenceCompute
  END IF
  IF (NiPos>0) THEN
    c=>Vec(NiPos)%c
    f=>Rhs(NiPos)%c
    CALL SubsidenceCompute
  END IF
  IF (NsPos>0) THEN
    c=>Vec(NsPos)%c
    f=>Rhs(NsPos)%c
    CALL SubsidenceCompute
  END IF

END SUBROUTINE SubsidenceScalar

SUBROUTINE SubsidenceCompute(DivRho,RhoTend)

  LOGICAL, OPTIONAL :: DivRho
  LOGICAL, OPTIONAL :: RhoTend

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp,ws

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        ws=0.5d0*(vSub1(iz)+vSub1(iz-1))
        IF (PRESENT(DivRho)) THEN ! Th, RhoV, RhoC, RhoR, RhoI, RhoS
          IF (ws<0.0d0) THEN 
            Temp=-ws*(c(ix,iy,iz+1,1)/(Rho(ix,iy,iz+1,1)+Eps)-c(ix,iy,iz  ,1)/(Rho(ix,iy,iz  ,1)+Eps))/dz(iz)
          ELSE 
            Temp=-ws*(c(ix,iy,iz,1  )/(Rho(ix,iy,iz,1  )+Eps)-c(ix,iy,iz-1,1)/(Rho(ix,iy,iz-1,1)+Eps))/dz(iz)
          END IF
          f(ix,iy,iz,1)=f(ix,iy,iz,1)+Temp*Rho(ix,iy,iz,1)
          IF (PRESENT(RhoTend)) THEN ! RhoV, RhoC, RhoR, RhoI, RhoS
            RhoRhs(ix,iy,iz,1)=RhoRhs(ix,iy,iz,1)+Temp*Rho(ix,iy,iz,1)
          END IF
        ELSE ! number concentrations Nv, Nc, Nr, Ni, Ns 
          IF (ws<0.0d0) THEN 
            Temp=-ws*(c(ix,iy,iz+1,1)-c(ix,iy,iz  ,1))/dz(iz)
          ELSE 
            Temp=-ws*(c(ix,iy,iz  ,1)-c(ix,iy,iz-1,1))/dz(iz)
          END IF
          f(ix,iy,iz,1)=f(ix,iy,iz,1)+Temp
        END IF
      END DO
    END DO
  END DO
END SUBROUTINE SubsidenceCompute

END MODULE Advection_Mod
