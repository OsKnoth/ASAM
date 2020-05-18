MODULE Buoyancy_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Control_Mod

  IMPLICIT NONE
  !REAL(RealKind), POINTER:: RhoC(:,:,:,:)

CONTAINS

SUBROUTINE AuftriebDryCompute

  INTEGER :: ix,iy,iz

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        f(ix,iy,iz,1)=f(ix,iy,iz,1)- &
            Rho(ix,iy,iz,1)*Grav 
      END DO
    END DO
  END DO

END SUBROUTINE AuftriebDryCompute

SUBROUTINE AuftriebDryComputeLR
  IF (.NOT.Anelastic) THEN
    IF (GradFull) THEN
     CALL AuftriebDryComputeLRFull
    ELSE
      CALL AuftriebDryComputeLRVol
    END IF
  ELSE
    IF (GradFull) THEN
      CALL BousinesqComputeLRFull
    ELSE  
      CALL BousinesqComputeLRVol
    END IF
  END IF

END SUBROUTINE AuftriebDryComputeLR

SUBROUTINE AuftriebDryComputeLRFull

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: RhoLoc,F
  REAL(RealKind) :: RhoFine,VolFine
  REAL(RealKind) :: RhoCoarse,VolCoarse

  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)*(dz(iz)*Rho(ix,iy,iz,1) &
               +dz(iz+1)*Rho(ix,iy,iz+1,1)) &
               /(dz(iz)+dz(iz+1)+Eps) 
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1)- &
            RhoLoc*Grav
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1)- &
            RhoLoc*Grav
      END DO
    END DO
  END DO
  DO in=1,AnzahlNachbar
    CALL Set(Nachbars(in))
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
            RhoLoc=F/(F+Eps)*(RhoFine*dz(iz0+1)+RhoCoarse*dLoc)/(dz(iz0+1)+dLoc+Eps)

            fL(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)=fL(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1) &
               -Grav*RhoLoc
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            RhoLoc=FW(jx,jy,iz0)/(FW(jx,jy,iz0)+Eps)*(Rho(jx,jy,iz0+1,1)*dz(iz0+1)+Rho(jx,jy,iz0,1)*dLoc)/(dz(iz0+1)+dLoc+Eps)
            fL(jx,jy,iz0+1,1)=fL(jx,jy,iz0+1,1)-Grav*RhoLoc
          END DO
        END DO
      END IF
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
            RhoLoc=F/(F+Eps)*(RhoCoarse*dLoc+RhoFine*dz(iz1))/(dLoc+dz(iz1)+Eps)
            fR(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)=fR(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1) &
               -Grav*RhoLoc
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            RhoLoc=FW(jx,jy,iz1)/(FW(jx,jy,iz1)+Eps)*(Rho(jx,jy,iz1,1)*dz(iz1)+Rho(jx,jy,iz1+1,1)*dLoc)/(dz(iz1)+dLoc+Eps)
            fR(jx,jy,iz1,1)=fR(jx,jy,iz1,1)-Grav*RhoLoc
          END DO
        END DO
      END IF
    END IF
  END DO  

END SUBROUTINE AuftriebDryComputeLRFull

SUBROUTINE AuftriebDryComputeLRVol

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc

  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*Rho(ix,iy,iz,1) &
               +VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1)- &
            RhoLoc*Grav
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1)- &
            RhoLoc*Grav
      END DO
    END DO
  END DO
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*Rho(ix,iy,iz,1) &
               +VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1)- &
            RhoLoc*Grav
      END DO
    END DO
  END IF
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*Rho(ix,iy,iz,1) &
               +VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1)- &
            RhoLoc*Grav
      END DO
    END DO
  END IF

END SUBROUTINE AuftriebDryComputeLRVol

SUBROUTINE HeightComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: HeightLoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        HeightLoc=(Rho(ix,iy,iz,1)+Rho(ix+1,iy,iz,1)) &
               *FU(ix,iy,iz)/(FU(ix,iy,iz)+Eps)*(HeightC(ix+1,iy,iz,1)-HeightC(ix,iy,iz,1)) &
               /(MetrXY(iy)*(dx(ix)+dx(ix+1)+Eps))
        uRhsL(ix+1,iy,iz,1)=uRhsL(ix+1,iy,iz,1)- &
            HeightLoc*Grav
        uRhsR(ix,iy,iz,1)=uRhsR(ix,iy,iz,1)- &
            HeightLoc*Grav
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        HeightLoc=(Rho(ix,iy,iz,1)+Rho(ix,iy+1,iz,1)) &
               *FV(ix,iy,iz)/(FV(ix,iy,iz)+Eps)*(HeightC(ix,iy+1,iz,1)-HeightC(ix,iy,iz,1)) &
               /(MetrYX(ix)*(dy(iy+1)+dy(iy)+Eps))
        vRhsL(ix,iy+1,iz,1)=vRhsL(ix,iy+1,iz,1)- &
            HeightLoc*Grav
        vRhsR(ix,iy,iz,1)=vRhsR(ix,iy,iz,1)- &
            HeightLoc*Grav
      END DO
    END DO
  END DO

END SUBROUTINE HeightComputeLR

SUBROUTINE AuftriebComputeF
  SELECT CASE(ThetaKind)
    CASE ('Pseudo')
      CALL BousinesqComputeF
    CASE DEFAULT
      IF (GradFull) THEN
        CALL AuftriebDryComputeFFull
      ELSE  
        CALL AuftriebDryComputeFVol
      END IF  
  END SELECT    
END SUBROUTINE AuftriebComputeF

SUBROUTINE AuftriebDryComputeFFull

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: RhoLoc,F
  REAL(RealKind) :: RhoFine,VolFine
  REAL(RealKind) :: RhoCoarse,VolCoarse

  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)*(dz(iz)*Rho(ix,iy,iz,1) &
               +dz(iz+1)*Rho(ix,iy,iz+1,1)) &
               /(dz(iz)+dz(iz+1)+Eps) 
        wfRhs(ix,iy,iz)=wfRhs(ix,iy,iz)- &
            RhoLoc*Grav
      END DO
    END DO
  END DO
  DO in=1,AnzahlNachbar
    CALL Set(Nachbars(in))
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
            RhoLoc=F/(F+Eps)*(RhoFine*dz(iz0+1)+RhoCoarse*dLoc)/(dz(iz0+1)+dLoc+Eps)
            wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)=wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0) &
               -Grav*RhoLoc
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            RhoLoc=FW(jx,jy,iz0)/(FW(jx,jy,iz0)+Eps)*(Rho(jx,jy,iz0+1,1)*dz(iz0+1)+Rho(jx,jy,iz0,1)*dLoc)/(dz(iz0+1)+dLoc+Eps)
            wFRhs(jx,jy,iz0)=wFRhs(jx,jy,iz0)-Grav*RhoLoc
          END DO
        END DO
      END IF
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
            RhoLoc=F/(F+Eps)*(RhoCoarse*dLoc+RhoFine*dz(iz1))/(dLoc+dz(iz1)+Eps)
            wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)=wFRhs(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1) &
               -Grav*RhoLoc
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            RhoLoc=FW(jx,jy,iz1)/(FW(jx,jy,iz1)+Eps)*(Rho(jx,jy,iz1,1)*dz(iz1)+Rho(jx,jy,iz1+1,1)*dLoc)/(dz(iz1)+dLoc+Eps)
            wFRhs(jx,jy,iz1)=wFRhs(jx,jy,iz1)-Grav*RhoLoc
          END DO
        END DO
      END IF
    END IF
  END DO  

END SUBROUTINE AuftriebDryComputeFFull

SUBROUTINE AuftriebDryComputeFVol

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc

  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*Rho(ix,iy,iz,1) &
               +VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)- &
            RhoLoc*Grav
      END DO
    END DO
  END DO
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*Rho(ix,iy,iz,1) &
               +VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)- &
            RhoLoc*Grav
      END DO
    END DO
  END IF
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*Rho(ix,iy,iz,1) &
               +VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1)) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)- &
            RhoLoc*Grav
      END DO
    END DO
  END IF

END SUBROUTINE AuftriebDryComputeFVol

SUBROUTINE JacAuftriebDryComputeLR

  IF (Anelastic) THEN
    CALL JacBousinesqComputeLR
  ELSE
    CALL JacAuftriebDryComComputeLR
  END IF
END SUBROUTINE JacAuftriebDryComputeLR

SUBROUTINE JacAuftriebDryCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: thL,thC,thR
  REAL(RealKind) :: GradL,GradR
  REAL(RealKind) :: VC,VR,VL
  REAL(RealKind) :: FR,FL

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        AS(IndexMet(wPosLJac,rhoPosJac))%c(ix,iy,iz,1)=-Grav 
        AS(IndexMet(wPosRJac,rhoPosJac))%c(ix,iy,iz,1)=-Grav 
      END DO
    END DO
  END DO

END SUBROUTINE JacAuftriebDryCompute

SUBROUTINE JacAuftriebDryComComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: thL,thC,thR
  REAL(RealKind) :: GradL,GradR
  REAL(RealKind) :: VC,VR,VL
  REAL(RealKind) :: FR,FL
  REAL(RealKind) :: Temp

  IF (.NOT.MultiTriTB.AND..NOT.MultiTriTR.AND..NOT.MultiMuTR) THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          AS(IndexMet(wPosLJac,rhoPosJac))%c(ix,iy,iz,1)=-Grav
          AS(IndexMet(wPosRJac,rhoPosJac))%c(ix,iy,iz,1)=-Grav
        END DO
      END DO
    END DO
  END IF

END SUBROUTINE JacAuftriebDryComComputeLR

SUBROUTINE AuftriebCloudComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc

  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*(Rho(ix,iy,iz,1) &
               +RhoV(ix,iy,iz,1)   &
               +RhoC(ix,iy,iz,1)   &
                             )   &
               +VolC(ix,iy,iz+1)*(Rho(ix,iy,iz+1,1) &
               +RhoV(ix,iy,iz+1,1)   &
               +RhoC(ix,iy,iz+1,1)   &
                             ))    &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1) &
                        -RhoLoc*Grav
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1) &
                      -RhoLoc*Grav
      END DO
    END DO
  END DO
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*(Rho(ix,iy,iz,1) &
               +RhoV(ix,iy,iz,1)   &
               +RhoC(ix,iy,iz,1))   &
               +VolC(ix,iy,iz+1)*(Rho(ix,iy,iz+1,1) &
               +RhoV(ix,iy,iz+1,1)   &
               +RhoC(ix,iy,iz+1,1)))   &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1) &
                      -RhoLoc*Grav
      END DO
    END DO
  END IF
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*(Rho(ix,iy,iz,1) &
               +RhoV(ix,iy,iz,1)   &
               +RhoC(ix,iy,iz,1))   &
               +VolC(ix,iy,iz+1)*(Rho(ix,iy,iz+1,1) &
               +RhoV(ix,iy,iz+1,1)   &
               +RhoC(ix,iy,iz+1,1)))   &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1) &
                        -RhoLoc*Grav
      END DO
    END DO
  END IF

END SUBROUTINE AuftriebCloudComputeLR

SUBROUTINE JacAuftriebCloudComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: thL,thC,thR
  REAL(RealKind) :: GradL,GradR
  REAL(RealKind) :: VC,VR,VL
  REAL(RealKind) :: FR,FL
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        AS(IndexMet(wPosLJac,rhoPosJac))%c(ix,iy,iz,1)=-Grav
        AS(IndexMet(wPosLJac,RhoVPosJac))%c(ix,iy,iz,1)=-Grav
        AS(IndexMet(wPosLJac,RhoCPosJac))%c(ix,iy,iz,1)=-Grav
        AS(IndexMet(wPosRJac,RhoPosJac))%c(ix,iy,iz,1)=-Grav
        AS(IndexMet(wPosRJac,RhoVPosJac))%c(ix,iy,iz,1)=-Grav
        AS(IndexMet(wPosRJac,RhoCPosJac))%c(ix,iy,iz,1)=-Grav
      END DO
    END DO
  END DO

END SUBROUTINE JacAuftriebCloudComputeLR

SUBROUTINE BousinesqComputeLRFull

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        RhoLoc=FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps) &
              *(dz(iz)*BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +dz(iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(dz(iz)+dz(iz+1)+Eps)
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1)+ &
            RhoLoc*Grav
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1)+ &
            RhoLoc*Grav
      END DO
    END DO
  END DO
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps) &
              *(dz(iz)*BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +dz(iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(dz(iz)+dz(iz+1)+Eps)
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1)+ &
            RhoLoc*Grav
      END DO
    END DO
  END IF
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps) &
              *(dz(iz)**BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +dz(iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(dz(iz)+dz(iz+1)+Eps)
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1)+ &
            RhoLoc*Grav
      END DO
    END DO
  END IF

END SUBROUTINE BousinesqComputeLRFull

SUBROUTINE BousinesqComputeLRVol

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc

  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +VolC(ix,iy,iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1)+ &
            RhoLoc*Grav
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1)+ &
            RhoLoc*Grav
      END DO
    END DO
  END DO
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +VolC(ix,iy,iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fL(ix,iy,iz+1,1)=fL(ix,iy,iz+1,1)+ &
            RhoLoc*Grav
      END DO
    END DO
  END IF
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +VolC(ix,iy,iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        fR(ix,iy,iz,1)=fR(ix,iy,iz,1)+ &
            RhoLoc*Grav
      END DO
    END DO
  END IF

END SUBROUTINE BousinesqComputeLRVol

SUBROUTINE BousinesqComputeF

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc

  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +VolC(ix,iy,iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)+ &
            RhoLoc*Grav
      END DO
    END DO
  END DO
  IF (TypeB(1:1)/='o') THEN
    iz=iz0
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +VolC(ix,iy,iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)+ &
            RhoLoc*Grav
      END DO
    END DO
  END IF
  IF (TypeT(1:1)/='o') THEN
    iz=iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=(VolC(ix,iy,iz)*BousinesqF(Th(ix,iy,iz,1),ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1)) &
               +VolC(ix,iy,iz+1)*BousinesqF(Th(ix,iy,iz+1,1),ThProf(ix,iy,iz+1,1),Rho(ix,iy,iz+1,1))) &
               /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
        wFRhs(ix,iy,iz)=wFRhs(ix,iy,iz)+ &
            RhoLoc*Grav
      END DO
    END DO
  END IF

END SUBROUTINE BousinesqComputeF

FUNCTION BousinesqF(T,TProf,Rho)
  REAL(RealKind) BousinesqF,T,TProf,Rho
  
  BousinesqF=Rho
  SELECT CASE(BousinesqType) 
    CASE('Scaled')
      BousinesqF=Rho*(T-TProf)/(TProf+Eps)
    CASE('UnScaled')
      BousinesqF=alphaBous*(T-TProf)
  END SELECT    
END FUNCTION BousinesqF

FUNCTION BousinesqF1(T,TProf,Rho)
  REAL(RealKind) BousinesqF1,T,TProf,Rho
  
  BousinesqF1=Rho
  SELECT CASE(BousinesqType) 
    CASE('Scaled')
      BousinesqF1=Rho*(T-TProf)/(Tprof+Eps)+Rho
    CASE('UnScaled')
      BousinesqF1=-alphaBous*(T-TProf)+Rho
  END SELECT    
END FUNCTION BousinesqF1


FUNCTION BousinesqF2(T,TProf,Rho)
  REAL(RealKind) BousinesqF2,T,TProf,Rho
  
  BousinesqF2=Rho
  SELECT CASE(BousinesqType) 
    CASE('Scaled')
      BousinesqF2=Rho*(T-TProf)/(Tprof+Eps)+Rho
    CASE('UnScaled')
      BousinesqF2=-alphaBous*(T-TProf)
  END SELECT    
END FUNCTION BousinesqF2
SUBROUTINE JacBousinesqComputeLR

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: thL,thC,thR
  REAL(RealKind) :: GradL,GradR
  REAL(RealKind) :: VC,VR,VL
  REAL(RealKind) :: FR,FL
  REAL(RealKind) :: Temp

  IF (.NOT.MultiTriTB.AND..NOT.MultiTriTR.AND..NOT.MultiMuTr) THEN
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          AS(IndexMet(wPosLJac,thPosJac))%c(ix,iy,iz,1) &
            =Grav*BousinesqJac(ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1))
          AS(IndexMet(wPosRJac,thPosJac))%c(ix,iy,iz,1) &
            =Grav*BousinesqJac(ThProf(ix,iy,iz,1),Rho(ix,iy,iz,1))
        END DO
      END DO
    END DO
  END IF

END SUBROUTINE JacBousinesqComputeLR

FUNCTION BousinesqJac(TProf,Rho)
  REAL(RealKind) BousinesqJac,TProf,Rho
  
  BousinesqJac=One
  SELECT CASE(BousinesqType) 
    CASE('Scaled')
      BousinesqJac=Rho/(TProf+Eps)
    CASE('UnScaled')
      BousinesqJac=alphaBous
  END SELECT    
END FUNCTION BousinesqJac

FUNCTION BousinesqJac2(TProf,Rho)
  REAL(RealKind) BousinesqJac2,TProf,Rho
  
  BousinesqJac2=One
  SELECT CASE(BousinesqType) 
    CASE('Scaled')
      BousinesqJac2=Rho/(TProf+Eps)
    CASE('UnScaled')
      BousinesqJac2=-alphaBous
  END SELECT    
END FUNCTION BousinesqJac2

END MODULE Buoyancy_Mod