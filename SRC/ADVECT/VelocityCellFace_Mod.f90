MODULE VelocityCellFace_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Example_Mod

  IMPLICIT NONE


CONTAINS

SUBROUTINE VelocityCellToFace(VectorCell,VelocityFace)

  TYPE(Vector4Cell_T), TARGET :: VectorCell(:)
  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)

  VectorCellAct=>VectorCell
  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVectorCell(ibLoc)
    CALL SetVelocityFace(ibLoc)
    CALL VelCellToFaceCompute
  END DO
END SUBROUTINE VelocityCellToFace

SUBROUTINE VelocityCellToFaceLR(VectorCell,VelocityFace,VelFOld,Time)

  TYPE(Vector4Cell_T), TARGET :: VectorCell(:)
  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)
  TYPE(VelocityFace_T), TARGET :: VelFOld(:)
  REAL(RealKind) :: Time

  INTEGER :: iz

  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    uFOld=>VelFOld(ibLoc)%uF
    vFOld=>VelFOld(ibLoc)%vF
    wFOld=>VelFOld(ibLoc)%wF
    uCL=>VectorCell(ibLoc)%Vec(uPosL)%c
    vCL=>VectorCell(ibLoc)%Vec(vPosL)%c
    wCL=>VectorCell(ibLoc)%Vec(wPosL)%c
    uCR=>VectorCell(ibLoc)%Vec(uPosR)%c
    vCR=>VectorCell(ibLoc)%Vec(vPosR)%c
    wCR=>VectorCell(ibLoc)%Vec(wPosR)%c
    CALL VelCellToFaceComputeLR(Time)
  END DO
END SUBROUTINE VelocityCellToFaceLR

SUBROUTINE VelocityFaceToCell(VelocityFace,VectorCell)

  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)
  TYPE(Vector4Cell_T), TARGET :: VectorCell(:)

  INTEGER :: ic

  VectorCellAct=>VectorCell
  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVectorCell(ibLoc)
    CALL SetVelocityFace(ibLoc)
    CALL VelFaceToCellCompute
  END DO
END SUBROUTINE VelocityFaceToCell

SUBROUTINE VelocityFaceToCellLR(VelocityFace,VectorCell)

  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)
  TYPE(Vector4Cell_T), TARGET :: VectorCell(:)

  INTEGER :: ic

  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    uCL=>VectorCell(ibLoc)%Vec(uPosL)%c
    vCL=>VectorCell(ibLoc)%Vec(vPosL)%c
    wCL=>VectorCell(ibLoc)%Vec(wPosL)%c
    uCR=>VectorCell(ibLoc)%Vec(uPosR)%c
    vCR=>VectorCell(ibLoc)%Vec(vPosR)%c
    wCR=>VectorCell(ibLoc)%Vec(wPosR)%c
    Rho=>RhoCell(ibLoc)%c
    CALL VelFaceToCellComputeLR
!   CALL VelFaceToCellComputeLRC 
  END DO
  CALL ExchangeCell(VectorCell)
END SUBROUTINE VelocityFaceToCellLR

SUBROUTINE SetVectorCell(ibLoc)

  INTEGER, INTENT(IN) :: ibLoc

  uC=>VectorCellAct(ibLoc)%Vec(uPosL)%c
  vC=>VectorCellAct(ibLoc)%Vec(vPosL)%c
  wC=>VectorCellAct(ibLoc)%Vec(wPosL)%c

END SUBROUTINE SetVectorCell

SUBROUTINE VelCellToFaceCompute

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: VolCoarse,VolFine
  REAL(RealKind) :: uCoarse,uFine,uFace
  REAL(RealKind) :: vCoarse,vFine,vFace
  REAL(RealKind) :: wCoarse,wFine,wFace

  uF=Zero
  vF=Zero
  wF=Zero
  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1)) 
            uFine= &
               SUM(uC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps 
             uCoarse=uC(ix0,jy,jz,1)
             uFace=(uCoarse*VolCoarse+uFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps) 
             uF(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                uFace*FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
                      (FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            uF(ix0,jy,jz)= &
                          (uC(ix0,jy,jz,1)*VolC(ix0,jy,jz)+ &
                           uC(ix0+1,jy,jz,1)*VolC(ix0+1,jy,jz))/ &
                          (VolC(ix0+1,jy,jz)+VolC(ix0,jy,jz)+Eps)* &
                           FU(ix0,jy,jz)/(FU(ix0,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF

    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)) 
            uFine= &
               SUM(uC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps 
             uCoarse=uC(ix1+1,jy,jz,1)
             uFace=(uCoarse*VolCoarse+uFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps) 
             uF(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                uFace*FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
                      (FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            uF(ix1,jy,jz)= &
                          (uC(ix1,jy,jz,1)*VolC(ix1,jy,jz)+ &
                           uC(ix1+1,jy,jz,1)*VolC(ix1+1,jy,jz))/ &
                          (VolC(ix1+1,jy,jz)+VolC(ix1,jy,jz)+Eps)* &
                           FU(ix1,jy,jz)/(FU(ix1,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF

    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)) 
            vFine= &
               SUM(vC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps 
             vCoarse=vC(jx,iy0,jz,1)
             vFace=(vCoarse*VolCoarse+vFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps) 
             vF(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)= &
                vFace*FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)/ &
                      (FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            vF(jx,iy0,jz)= &
                          (vC(jx,iy0,jz,1)*VolC(jx,iy0,jz)+ &
                           vC(jx,iy0+1,jz,1)*VolC(jx,iy0+1,jz))/ &
                          (VolC(jx,iy0+1,jz)+VolC(jx,iy0,jz)+Eps)* &
                           FV(jx,iy0,jz)/(FV(jx,iy0,jz)+Eps)
          END DO
        END DO
      END IF
    END IF

    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)) 
            vFine= &
               SUM(vC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps 
             vCoarse=vC(jx,iy1+1,jz,1)
             vFace=(vCoarse*VolCoarse+vFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps) 
             vF(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
                vFace*FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)/ &
                      (FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jz=jz0+1,jz1
          DO jx=jx0+1,jx1
            vF(jx,iy1,jz)= &
                          (vC(jx,iy1,jz,1)*VolC(jx,iy1,jz)+ &
                           vC(jx,iy1+1,jz,1)*VolC(jx,iy1+1,jz))/ &
                          (VolC(jx,iy1+1,jz)+VolC(jx,iy1,jz)+Eps)* &
                           FV(jx,iy1,jz)/(FV(jx,iy1,jz)+Eps)
          END DO
        END DO
      END IF
    END IF

    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jy=jy0+1,jy1,IncrY
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)) 
            wFine= &
               SUM(wC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps 
             wCoarse=wC(jx,jy,iz0,1)
             wFace=(wCoarse*VolCoarse+wFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps) 
             wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)= &
                wFace*FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)/ &
                      (FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            wF(jx,jy,iz0)= &
                          (wC(jx,jy,iz0,1)*VolC(jx,jy,iz0)+ &
                           wC(jx,jy,iz0+1,1)*VolC(jx,jy,iz0+1))/ &
                          (VolC(jx,jy,iz0+1)+VolC(jx,jy,iz0)+Eps)* &
                           FW(jx,jy,iz0)/(FW(jx,jy,iz0)+Eps)
          END DO
        END DO
      END IF
    END IF

    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jy=jy0+1,jy1,IncrY
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)) 
            wFine= &
               SUM(wC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps 
             wCoarse=wC(jx,jy,iz1+1,1)
             wFace=(wCoarse*VolCoarse+wFine*VolFine)/ &
                   (VolFine+VolCoarse+Eps) 
             wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
                wFace*FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)/ &
                      (FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            wF(jx,jy,iz1)= &
                          (wC(jx,jy,iz1,1)*VolC(jx,jy,iz1)+ &
                           wC(jx,jy,iz1+1,1)*VolC(jx,jy,iz1+1))/ &
                          (VolC(jx,jy,iz1+1)+VolC(jx,jy,iz1)+Eps)* &
                           FW(jx,jy,iz1)/(FW(jx,jy,iz1)+Eps)
          END DO
        END DO
      END IF
    END IF
  END DO


  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        uF(ix,iy,iz)= &
                     (uC(ix,iy,iz,1)*VolC(ix,iy,iz)+ &
                      uC(ix+1,iy,iz,1)*VolC(ix+1,iy,iz))/ &
                     (VolC(ix+1,iy,iz)+VolC(ix,iy,iz)+Eps)* &
                      FU(ix,iy,iz)/(FU(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        vF(ix,iy,iz)= &
                     (vC(ix,iy,iz,1)*VolC(ix,iy,iz)+ &
                      vC(ix,iy+1,iz,1)*VolC(ix,iy+1,iz))/ &
                     (VolC(ix,iy+1,iz)+VolC(ix,iy,iz)+Eps)* &
                      FV(ix,iy,iz)/(FV(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        wF(ix,iy,iz)= &
                     (wC(ix,iy,iz,1)*VolC(ix,iy,iz)+ &
                      wC(ix,iy,iz+1,1)*VolC(ix,iy,iz))/ &
                     (VolC(ix,iy,iz+1)+VolC(ix,iy,iz)+Eps)* &
                      FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO

  IF (TypeW=='ow') THEN
    IF (BCMetVec(uPosL)%West=='ZeroGrad'.OR.BCMetVec(uPosL)%West=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix0,iy,iz)=uC(ix0+1,iy,iz,1) &
                                *FU(ix0,iy,iz)/(FU(ix0,iy,iz)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeE=='oe') THEN
    IF (BCMetVec(uPosL)%East=='ZeroGrad'.OR.BCMetVec(uPosL)%East=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix1,iy,iz)=uC(ix1,iy,iz,1) &
                                *FU(ix1,iy,iz)/(FU(ix1,iy,iz)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeS=='os') THEN
    IF (BCMetVec(vPosL)%South=='ZeroGrad'.OR.BCMetVec(vPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy0,iz)=vC(ix,iy0+1,iz,1) &
                                *FV(ix,iy0,iz)/(FV(ix,iy0,iz)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeN=='on') THEN
    IF (BCMetVec(vPosL)%North=='ZeroGrad'.OR.BCMetVec(vPosL)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy1,iz)=vC(ix,iy1,iz,1) &
                                *FV(ix,iy1,iz)/(FV(ix,iy1,iz)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeB=='ob') THEN
    IF (BCMetVec(wPosL)%Bottom=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz0)=wC(ix,iy,iz0+1,1) &
                                *FW(ix,iy,iz0)/(FW(ix,iy,iz0)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeT=='ot') THEN
    IF (BCMetVec(wPosL)%Top=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz1)=wC(ix,iy,iz1,1) &
                                *FW(ix,iy,iz1)/(FW(ix,iy,iz1)+eps)
        END DO
      END DO
    END IF
  END IF
  
END SUBROUTINE VelCellToFaceCompute

SUBROUTINE VelFaceToCellCompute

  REAL(RealKind) :: Time 

  INTEGER :: ix,iy,iz

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        uC(ix,iy,iz,1)=          &
                (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
                 FU(ix,iy,iz)*uF(ix,iy,iz)) / &
                (FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        vC(ix,iy,iz,1)=          &
                (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
                 FV(ix,iy,iz)*vF(ix,iy,iz)) / &
                (FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        wC(ix,iy,iz,1)=          &
                (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
                 FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
                (FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE VelFaceToCellCompute

SUBROUTINE VelFaceToCellComputeLRC

  REAL(RealKind) :: Time 

  INTEGER :: ix,iy,iz

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        uCL(ix,iy,iz,1)=          &
                (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
                 FU(ix,iy,iz)*uF(ix,iy,iz)) / &
                (FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)
        uCR(ix,iy,iz,1)=uCL(ix,iy,iz,1)
        vCL(ix,iy,iz,1)=          &
                (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
                 FV(ix,iy,iz)*vF(ix,iy,iz)) / &
                (FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)
        vCR(ix,iy,iz,1)=vCL(ix,iy,iz,1)
        wCL(ix,iy,iz,1)=          &
                (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
                 FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
                (FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)
        wCR(ix,iy,iz,1)=wCL(ix,iy,iz,1)
      END DO
    END DO
  END DO
END SUBROUTINE VelFaceToCellComputeLRC

SUBROUTINE VelFaceToCellComputeLR

  INTEGER :: i
  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: uFine,vFine,wFine,RhoFine
  REAL(RealKind) :: FUL,FUR
  REAL(RealKind) :: FVL,FVR
  REAL(RealKind) :: FWL,FWR
  REAL(RealKind) :: uLoc,vLoc,wLoc
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: nu,nn

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        uCL(ix,iy,iz,1)=Rho(ix,iy,iz,1) &
                       *(VolC(ix-1,iy,iz)+VolC(ix,iy,iz)) &
                       /(Rho(ix-1,iy,iz,1)*VolC(ix-1,iy,iz) &
                        +Rho(ix,iy,iz,1)*VolC(ix,iy,iz)+Eps) &
                        *uF(ix-1,iy,iz)          
        uCR(ix,iy,iz,1)=Rho(ix,iy,iz,1) &
                       *(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)) &
                       /(Rho(ix+1,iy,iz,1)*VolC(ix+1,iy,iz) &
                        +Rho(ix,iy,iz,1)*VolC(ix,iy,iz)+Eps) &
                        *uF(ix,iy,iz)          
        vCL(ix,iy,iz,1)=Rho(ix,iy,iz,1) &
                       *(VolC(ix,iy-1,iz)+VolC(ix,iy,iz)) &
                       /(Rho(ix,iy-1,iz,1)*VolC(ix,iy-1,iz) &
                        +Rho(ix,iy,iz,1)*VolC(ix,iy,iz)+Eps) &
                        *vF(ix,iy-1,iz)          
        vCR(ix,iy,iz,1)=Rho(ix,iy,iz,1) &
                       *(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)) &
                       /(Rho(ix,iy+1,iz,1)*VolC(ix,iy+1,iz) &
                        +Rho(ix,iy,iz,1)*VolC(ix,iy,iz)+Eps) &
                        *vF(ix,iy,iz)          
        wCL(ix,iy,iz,1)=Rho(ix,iy,iz,1) &
                       *(VolC(ix,iy,iz-1)+VolC(ix,iy,iz)) &
                       /(Rho(ix,iy,iz-1,1)*VolC(ix,iy,iz-1) &
                        +Rho(ix,iy,iz,1)*VolC(ix,iy,iz)+Eps) &
                        *wF(ix,iy,iz-1)          
        wCR(ix,iy,iz,1)=Rho(ix,iy,iz,1) &
                       *(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)) &
                       /(Rho(ix,iy,iz+1,1)*VolC(ix,iy,iz+1) &
                        +Rho(ix,iy,iz,1)*VolC(ix,iy,iz)+Eps) &
                        *wF(ix,iy,iz)          
        uCL(ix,iy,iz,1)=uF(ix-1,iy,iz)          
        uCR(ix,iy,iz,1)=uF(ix,iy,iz)          
        vCL(ix,iy,iz,1)=vF(ix,iy-1,iz)          
        vCR(ix,iy,iz,1)=vF(ix,iy,iz)          
        wCL(ix,iy,iz,1)=wF(ix,iy,iz-1)          
        wCR(ix,iy,iz,1)=wF(ix,iy,iz)          
      END DO
    END DO
  END DO
  IF (Sphere.AND.igy0==domain%igy0) THEN
    iy=iy0+1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        vCL(ix,iy,iz,1)=vF(ix,iy,iz)          
      END DO
    END DO
  END IF  
  IF (Sphere.AND.igy1==domain%igy1) THEN
    iy=iy1
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        vCR(ix,iy,iz,1)=vF(ix,iy-1,iz)          
      END DO
    END DO
  END IF  

  IF (FreeSlip) THEN
    DO i=1,NumBoundCell
      ix    =BoundCell(i)%ix
      iy    =BoundCell(i)%iy
      iz    =BoundCell(i)%iz
      CALL ShiftVel1(VolC(ix,iy,iz), &  
                    FU(ix-1,iy,iz),FU(ix,iy,iz), &
                    xP(ix-1),xP(ix), &
                    FV(ix,iy-1,iz),FV(ix,iy,iz), &
                    yP(iy-1),yP(iy), &
                    FW(ix,iy,iz-1),FW(ix,iy,iz), &
                    zP(iz-1),zP(iz), &
                    uF(ix-1,iy,iz),uF(ix,iy,iz), &
                    vF(ix,iy-1,iz),vF(ix,iy,iz), &
                    wF(ix,iy,iz-1),wF(ix,iy,iz), &
                    uCL(ix,iy,iz,1),uCR(ix,iy,iz,1), &
                    vCL(ix,iy,iz,1),vCR(ix,iy,iz,1), &
                    wCL(ix,iy,iz,1),wCR(ix,iy,iz,1))
    END DO
  ELSE
    DO i=1,NumBoundCell
      ix    =BoundCell(i)%ix
      iy    =BoundCell(i)%iy
      iz    =BoundCell(i)%iz
      CALL ShiftNoSlip(VolC(ix,iy,iz), &
                    FU(ix-1,iy,iz),FU(ix,iy,iz), &
                    xP(ix-1),xP(ix), &
                    FV(ix,iy-1,iz),FV(ix,iy,iz), &
                    yP(iy-1),yP(iy), &
                    FW(ix,iy,iz-1),FW(ix,iy,iz), &
                    zP(iz-1),zP(iz), &
                    uF(ix-1,iy,iz),uF(ix,iy,iz), &
                    vF(ix,iy-1,iz),vF(ix,iy,iz), &
                    wF(ix,iy,iz-1),wF(ix,iy,iz), &
                    uCL(ix,iy,iz,1),uCR(ix,iy,iz,1), &
                    vCL(ix,iy,iz,1),vCR(ix,iy,iz,1), &
                    wCL(ix,iy,iz,1),wCR(ix,iy,iz,1))
    END DO
  END IF

END SUBROUTINE VelFaceToCellComputeLR

SUBROUTINE ShiftVel(Vol, &
                    FUL,FUR, &
                    xL,xR, &
                    FVL,FVR, &
                    yL,yR, &
                    FWL,FWR, &
                    zL,zR, &
                    uL,uR,vL,vR,wL,wR, &
                    uCL,uCR,vCL,vCR,wCL,wCR)
                    

  REAL(RealKind) :: Vol
  REAL(RealKind) :: FUL,FUR
  REAL(RealKind) :: xL,xR
  REAL(RealKind) :: FVL,FVR
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: FWL,FWR
  REAL(RealKind) :: zL,zR
  REAL(RealKind) :: uL,uR,vL,vR,wL,wR
  REAL(RealKind) :: uCL,uCR,vCL,vCR,wCL,wCR

  REAL(RealKind) :: FL,xFL,yFL,zFL
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: b(7)
  REAL(RealKind) :: A(7,4)
  INTEGER :: m,n
  INTEGER :: info

  INTEGER :: lwork=200
  REAL(RealKind) :: Work(200)
  
  m=7
  n=4

  IF (ABS(FUL-FUR)>1.d-7*MAX(FUL,FUR)) THEN
    xFL=MAX(MIN((Vol-FUR*xR+FUL*xL)/(FUL-FUR),xR),xL)
  ELSE
    xFL=Half*(xL+xR)
  END IF 
  IF (ABS(FVL-FVR)>1.d-7*MAX(FVL,FVR)) THEN
    yFL=MAX(MIN((Vol-FVR*yR+FVL*zL)/(FVL-FVR),yR),yL)
  ELSE
    yFL=Half*(yL+yR)
  END IF 
  IF (ABS(FWL-FWR)>1.d-7*MAX(FWL,FWR)) THEN
    zFL=MAX(MIN((Vol-FWR*zR+FWL*zL)/(FWL-FWR),zR),zL)
  ELSE
    zFL=Half*(zL+zR)
  END IF 

  n1=FUR-FUL 
  n2=FVR-FVL 
  n3=FWR-FWL 
  FL=SQRT(n1*n1+n2*n2+n3*n3) 

! FUL
  b(1)=FUL*uL
  A(1,1)=FUL
  A(1,2)=0.0d0
  A(1,3)=0.0d0
  A(1,4)=FUL*(xL-xFL)
! FUR
  b(2)=FUR*uR
  A(2,1)=FUR
  A(2,2)=0.0d0
  A(2,3)=0.0d0
  A(2,4)=FUR*(xR-xFL)
! FVL
  b(3)=FVL*vL
  A(3,1)=0.0d0
  A(3,2)=FVL
  A(3,3)=0.0d0
  A(3,4)=FVL*(yL-yFL)
! FVR
  b(4)=FVR*vR
  A(4,1)=0.0d0
  A(4,2)=FVR
  A(4,3)=0.0d0
  A(4,4)=FVR*(yR-yFL)
! FWL
  b(5)=FWL*wL
  A(5,1)=0.0d0
  A(5,2)=0.0d0
  A(5,3)=FWL
  A(5,4)=FWL*(zL-zFL)
! FWR
  b(6)=FWR*wR
  A(6,1)=0.0d0
  A(6,2)=0.0d0
  A(6,3)=FWR
  A(6,4)=FWR*(zR-zFL)
! FL 
  b(7)=0.0d0
  A(7,1)=n1
  A(7,2)=n2
  A(7,3)=n3
  A(7,4)=0.0d0
  A(7,:)=1.0d7*A(7,:)

  CALL DGELS('N',m,n,1,A,m,b,m,work,lwork,info) 
  IF (info>0) THEN
    WRITE(*,*) 'info aus DGELS',info
    STOP
  END IF

  IF (FUL>FUR) THEN
    uCL=uL 
    uCR=(FUR*uR+(FUL-FUR)*b(1))/FUL
  ELSE
    uCL=(FUL*uL+(FUR-FUL)*b(1))/FUR
    uCR=uR 
  END IF
  IF (FVL>FVR) THEN
    vCL=vL 
    vCR=(FVR*vR+(FVL-FVR)*b(2))/FVL
  ELSE
    vCL=(FVL*vL+(FVR-FVL)*b(2))/FVR
    vCR=vR 
  END IF
  IF (FWL>FWR) THEN
    wCL=wL 
    wCR=(FWR*wR+(FWL-FWR)*b(3))/FWL
  ELSE
    wCL=(FWL*wL+(FWR-FWL)*b(3))/FWR
    wCR=wR 
  END IF
    
END SUBROUTINE ShiftVel

SUBROUTINE ShiftVel1(Vol, &
                    FUL,FUR, &
                    xL,xR, &
                    FVL,FVR, &
                    yL,yR, &
                    FWL,FWR, &
                    zL,zR, &
                    uL,uR,vL,vR,wL,wR, &
                    uCL,uCR,vCL,vCR,wCL,wCR)
                    

  REAL(RealKind) :: Vol
  REAL(RealKind) :: FUL,FUR
  REAL(RealKind) :: xL,xR
  REAL(RealKind) :: FVL,FVR
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: FWL,FWR
  REAL(RealKind) :: zL,zR
  REAL(RealKind) :: uL,uR,vL,vR,wL,wR
  REAL(RealKind) :: uCL,uCR,vCL,vCR,wCL,wCR

  REAL(RealKind) :: FL,xFL,yFL,zFL
  REAL(RealKind) :: n1,n2,n3

  IF (FUL>FUR) THEN
    uCL=uL 
    uCR=(FUR*uR+(FUL-FUR)*uL)/(FUL+Eps)
  ELSE
    uCL=(FUL*uL+(FUR-FUL)*uR)/(FUR+Eps)
    uCR=uR 
  END IF
  IF (FVL>FVR) THEN
    vCL=vL 
    vCR=(FVR*vR+(FVL-FVR)*vL)/(FVL+Eps)
  ELSE
    vCL=(FVL*vL+(FVR-FVL)*vR)/(FVR+Eps)
    vCR=vR 
  END IF
  IF (FWL>FWR) THEN
    wCL=wL 
    wCR=(FWR*wR+(FWL-FWR)*wL)/(FWL+Eps)
  ELSE
    wCL=(FWL*wL+(FWR-FWL)*wR)/(FWR+Eps)
    wCR=wR 
  END IF
END SUBROUTINE ShiftVel1

SUBROUTINE ShiftVel2(Vol, &
                    FUL,FUR, &
                    xL,xR, &
                    FVL,FVR, &
                    yL,yR, &
                    FWL,FWR, &
                    zL,zR, &
                    uL,uR,vL,vR,wL,wR, &
                    uCL,uCR,vCL,vCR,wCL,wCR)


  REAL(RealKind) :: Vol
  REAL(RealKind) :: FUL,FUR
  REAL(RealKind) :: xL,xR
  REAL(RealKind) :: FVL,FVR
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: FWL,FWR
  REAL(RealKind) :: zL,zR
  REAL(RealKind) :: uL,uR,vL,vR,wL,wR
  REAL(RealKind) :: uCL,uCR,vCL,vCR,wCL,wCR

  REAL(RealKind) :: FL,xFL,yFL,zFL
  REAL(RealKind) :: n1,n2,n3
  REAL(RealKind) :: uLoc,vLoc,wLoc
  REAL(RealKind) :: nu,nn



  uLoc=(FUL*uL+FUR*uR)/(FUL+FUR+Eps)
  vLoc=(FVL*vL+FVR*vR)/(FVL+FVR+Eps)
  wLoc=(FWL*wL+FWR*wR)/(FWL+FWR+Eps)
! uLoc=(MAX(FUR-FUL,Zero)*uR+MAX(FUL-FUR,Zero)*uL)/(ABS(FUR-FUL)+Eps)
! vLoc=(MAX(FVR-FVL,Zero)*vR+MAX(FVL-FVR,Zero)*vL)/(ABS(FVR-FVL)+Eps)
! wLoc=(MAX(FWR-FWL,Zero)*wR+MAX(FWL-FWR,Zero)*wL)/(ABS(FWR-FWL)+Eps)
  n1=FUR-FUL
  n2=FVR-FVL
  n3=FWR-FWL
  nu=n1*uLoc+n2*vLoc+n3*wLoc
  nn=n1*n1+n2*n2+n3*n3
  uCL=(FUL*uL+MAX(FUR-FUL,Zero)*(uLoc-nu/nn*n1))/(MAX(FUL,FUR)+Eps)
  uCR=(FUR*uR+MAX(FUL-FUR,Zero)*(uLoc-nu/nn*n1))/(MAX(FUL,FUR)+Eps)
  vCL=(FVL*vL+MAX(FVR-FVL,Zero)*(vLoc-nu/nn*n2))/(MAX(FVL,FVR)+Eps)
  vCR=(FVR*vR+MAX(FVL-FVR,Zero)*(vLoc-nu/nn*n2))/(MAX(FVL,FVR)+Eps)
  wCL=(FWL*wL+MAX(FWR-FWL,Zero)*(wLoc-nu/nn*n3))/(MAX(FWL,FWR)+Eps)
  wCR=(FWR*wR+MAX(FWL-FWR,Zero)*(wLoc-nu/nn*n3))/(MAX(FWL,FWR)+Eps)

END SUBROUTINE ShiftVel2

SUBROUTINE ShiftNoSlip(Vol, &
                    FUL,FUR, &
                    xL,xR, &
                    FVL,FVR, &
                    yL,yR, &
                    FWL,FWR, &
                    zL,zR, &
                    uL,uR,vL,vR,wL,wR, &
                    uCL,uCR,vCL,vCR,wCL,wCR)


  REAL(RealKind) :: Vol
  REAL(RealKind) :: FUL,FUR
  REAL(RealKind) :: xL,xR
  REAL(RealKind) :: FVL,FVR
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: FWL,FWR
  REAL(RealKind) :: zL,zR
  REAL(RealKind) :: uL,uR,vL,vR,wL,wR
  REAL(RealKind) :: uCL,uCR,vCL,vCR,wCL,wCR


  uCL=FUL*uL/MAX(FUR,FUL)
  uCR=FUR*uR/MAX(FUR,FUL)
  vCL=FVL*vL/MAX(FVR,FVL)
  vCR=FVR*vR/MAX(FVR,FVL)
  wCL=FWL*wL/MAX(FWR,FWL)
  wCR=FWR*wR/MAX(FWR,FWL)
    
END SUBROUTINE ShiftNoSlip

SUBROUTINE VelCellToFaceComputeLR2(Time)

  REAL(RealKind) :: Time

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: VolCoarse,VolFine
  REAL(RealKind) :: FL,FC,FR
  REAL(RealKind) :: uCoarse,uFine,uFace
  REAL(RealKind) :: vCoarse,vFine,vFace
  REAL(RealKind) :: wCoarse,wFine,wFace
  REAL(RealKind) :: VL,VR

  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            FR= &
               SUM(FU(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            uFine= &
               SUM(uCL(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            FC= &
               SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            FL= &
               SUM(FU(ix0-1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            uCoarse=uCR(ix0,jy,jz,1)
            VL=VolCoarse
            VR=VolFine
            uFace=(uCoarse*VL+uFine*VR)/ &
                  (VL+VR+Eps)
            uF(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
               uFace*FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
                     (FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            VR=VolC(ix0+1,jy,jz)
            VL=VolC(ix0,jy,jz)
            uF(ix0,jy,jz)=(VR*uCL(ix0+1,jy,jz,1) &
                         +VL*uCR(ix0,jy,jz,1)) &
                         /(VL+VR+Eps)*FU(ix0,jy,jz)/(FU(ix0,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            FL= &
               SUM(FU(ix1-1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            uFine= &
               SUM(uCR(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            FC= &
               SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            FR= &
               SUM(FU(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            VL=VolFine
            VR=VolCoarse
            uCoarse=uCL(ix1+1,jy,jz,1)
            uFace=(uCoarse*VR+uFine*VL)/ &
                  (VL+VR+Eps)
            uF(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
               uFace*FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
                     (FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            VR=VolC(ix1+1,jy,jz)
            VL=VolC(ix1,jy,jz)
            uF(ix1,jy,jz)=(VR*uCL(ix1+1,jy,jz,1) &
                         +VL*uCR(ix1,jy,jz,1)) &
                         /(VL+VR+Eps)*FU(ix1,jy,jz)/(FU(ix1,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            FR= &
               SUM(FV(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            vFine= &
               SUM(vCL(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            FC= &
               SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            FL= &
               SUM(FV(jx:jx+IncrX-1,iy0-1,jz:jz+IncrZ-1))+Eps
            vCoarse=vCR(jx,iy0,jz,1)
            VL=VolCoarse
            VR=VolFine
            vFace=(vCoarse*VL+vFine*VR)/ &
                  (VL+VR+Eps)
            vF(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)= &
                vFace*FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)/ &
                      (FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            VR=VolC(jx,iy0+1,jz)
            VL=VolC(jx,iy0,jz)
            vF(jx,iy0,jz)=(VR*vCL(jx,iy0+1,jz,1) &
                          +VL*vCR(jx,iy0,jz,1)) &
                         /(VL+VR+Eps)*FV(jx,iy0,jz)/(FV(jx,iy0,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            FL= &
               SUM(FV(jx:jx+IncrX-1,iy1-1,jz:jz+IncrZ-1))
            vFine= &
               SUM(vCR(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            FC= &
               SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))+Eps
            FR= &
               SUM(FV(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            VL=VolFine
            VR=VolCoarse
            vCoarse=vCL(jx,iy1+1,jz,1)
            vFace=(vCoarse*VR+vFine*VL)/ &
                  (VL+VR+Eps)
            vF(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
               vFace*FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)/ &
                     (FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            VR=VolC(jx,iy1+1,jz)
            VL=VolC(jx,iy1,jz)
            vF(jx,iy1,jz)=(VR*vCL(jx,iy1+1,jz,1) &
                          +VL*vCR(jx,iy1,jz,1)) &
                         /(VL+VR+Eps)*FV(jx,iy1,jz)/(FV(jx,iy1,jz)+Eps)
          END DO
        END DO
      END IF
    END IF


    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            FR= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            wFine= &
               SUM(wCL(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            FC= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            FL= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0-1))
            wCoarse=wCR(jx,jy,iz0,1)
            VL=VolCoarse
            VR=VolFine
            wFace=(wCoarse*VL+wFine*VR)/ &
                  (VL+VR+Eps)
            wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)= &
                wFace*FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)/ &
                      (FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            VR=VolC(jx,jy,iz0+1)
            VL=VolC(jx,jy,iz0)
            wF(jx,jy,iz0)=(VR*wCL(jx,jy,iz0+1,1) &
                          +VL*wCR(jx,jy,iz0,1)) &
                         /(VL+VR+Eps)*FW(jx,jy,iz0)/(FW(jx,jy,iz0)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            FL= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1-1))
            wFine= &
               SUM(wCR(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            FC= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))+Eps
            FR= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            VL=VolFine
            VR=VolCoarse
            wCoarse=wCL(jx,jy,iz1+1,1)
            wFace=(wCoarse*VR+wFine*VL)/ &
                  (VL+VR+Eps)
            wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
                wFace*FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)/ &
                      (FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            VR=VolC(jx,jy,iz1+1)
            VL=VolC(jx,jy,iz1)
            wF(jx,jy,iz1)=(VR*wCL(jx,jy,iz1+1,1) &
                          +VL*wCR(jx,jy,iz1,1)) &
                         /(VL+VR+Eps)*FW(jx,jy,iz1)/(FW(jx,jy,iz1)+Eps)
          END DO
        END DO
      END IF
    END IF
  END DO

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        VR=VolC(ix+1,iy,iz)
        VL=VolC(ix,iy,iz)
        uF(ix,iy,iz)=(VR*uCL(ix+1,iy,iz,1) &
                     +VL*uCR(ix,iy,iz,1)) &
                     /(VL+VR+Eps)*FU(ix,iy,iz)/(FU(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        VR=VolC(ix,iy+1,iz)
        VL=VolC(ix,iy,iz)
        vF(ix,iy,iz)=(VR*vCL(ix,iy+1,iz,1) &
                     +VL*vCR(ix,iy,iz,1)) &
                     /(VL+VR+Eps)*FV(ix,iy,iz)/(FV(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VR=VolC(ix,iy,iz+1)
        VL=VolC(ix,iy,iz)
        wF(ix,iy,iz)=(VR*wCL(ix,iy,iz+1,1) &
                     +VL*wCR(ix,iy,iz,1)) &
                     /(VL+VR+Eps)*FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  IF (TypeW(1:1)=='o') THEN
    !IF (BCVel%West/='OutFlow') THEN
    IF (BCMetVec(uPosL)%West=='ZeroGrad'.OR.BCMetVec(uPosL)%West=='MeanValue') THEN ! Hinneburg
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix,iy,iz)=Zero
        END DO
      END DO
    ELSE
      ix=ix0
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix,iy,iz)=uF(ix+1,iy,iz)
        END DO
      END DO
    END IF
  END IF
  IF (TypeE(1:1)=='o') THEN
    IF (BCMetVec(uPosL)%East=='ZeroGrad'.OR.BCMetVec(uPosL)%East=='MeanValue') THEN 
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix,iy,iz)=Zero
        END DO
      END DO
    ELSE
      ix=ix1
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix,iy,iz)=uF(ix-1,iy,iz)
        END DO
      END DO
    END IF
  END IF
  IF (TypeS(1:1)=='o') THEN
    iy=iy0
    IF (BCMetVec(vPosL)%South=='ZeroGrad'.OR.BCMetVec(vPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy,iz)=Zero
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy,iz)=vF(ix,iy+1,iz)
        END DO
      END DO
    END IF
  END IF
  IF (TypeN(1:1)=='o') THEN
    iy=iy1
    IF (BCMetVec(vPosL)%North=='ZeroGrad'.OR.BCMetVec(vPosL)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy,iz)=Zero
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy,iz)=vF(ix,iy-1,iz)
        END DO
      END DO
    END IF
  END IF
  IF (TypeB(1:1)=='o') THEN
    iz=iz0
    IF (BCMetVec(wPosL)%Bottom=='ZeroGrad'.OR.BCMetVec(wPosL)%Bottom=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz)=Zero
        END DO
      END DO
    ELSE
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz)=wF(ix,iy,iz+1)
        END DO
      END DO
    END IF
  END IF
  IF (TypeT(1:1)=='o') THEN
    IF (BCMetVec(wPosL)%Top=='ZeroGrad'.OR.BCMetVec(wPosL)%Top=='MeanValue') THEN 
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz)=Zero
        END DO
      END DO
    ELSE 
      iz=iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz)=wF(ix,iy,iz-1)
        END DO
      END DO
    END IF
  END IF

END SUBROUTINE VelCellToFaceComputeLR2

SUBROUTINE VelCellToFaceComputeLR(Time)

  REAL(RealKind) :: Time

  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  INTEGER :: in
  REAL(RealKind) :: VolCoarse,VolFine
  REAL(RealKind) :: FL,FC,FR
  REAL(RealKind) :: uCoarse,uFine,uFace
  REAL(RealKind) :: vCoarse,vFine,vFace
  REAL(RealKind) :: wCoarse,wFine,wFace
  REAL(RealKind) :: VL,VR
  REAL(RealKind) :: xPL,yPL,zPL

  DO in=1,AnzahlNachbar

    CALL Set(Nachbars(in))

    IF (Nachbars(in)%nType=='iw'.OR.Nachbars(in)%nType=='pw') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            FR= &
               SUM(FU(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            uFine= &
               SUM(uCL(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            FC= &
               SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            FL= &
               SUM(FU(ix0-1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            uCoarse=uCR(ix0,jy,jz,1)
            VL=VolCoarse*FC/(FC+FL+Eps)
            VR=VolFine*FC/(FC+FR+Eps)
            uFace=(uCoarse*VL+uFine*VR)/ &
                  (VL+VR+Eps)
            uF(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
               uFace*FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
                     (FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            VR=VolC(ix0+1,jy,jz)*FU(ix0,jy,jz)/(FU(ix0,jy,jz)+FU(ix0+1,jy,jz)+Eps)
            VL=VolC(ix0,jy,jz)*FU(ix0,jy,jz)/(FU(ix0,jy,jz)+FU(ix0-1,jy,jz)+Eps)
            uF(ix0,jy,jz)=(VR*uCL(ix0+1,jy,jz,1) &
                         +VL*uCR(ix0,jy,jz,1)) &
                         /(VL+VR+Eps)*FU(ix0,jy,jz)/(FU(ix0,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='ie'.OR.Nachbars(in)%nType=='pe') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            FL= &
               SUM(FU(ix1-1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
            uFine= &
               SUM(uCR(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1,1)* &
                   VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            FC= &
               SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            FR= &
               SUM(FU(ix1+1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps
            VL=VolFine*FC/(FC+FL+Eps)
            VR=VolCoarse*FC/(FC+FR+Eps)
            uCoarse=uCL(ix1+1,jy,jz,1)
            uFace=(uCoarse*VR+uFine*VL)/ &
                  (VL+VR+Eps)
            uF(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
               uFace*FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
                     (FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jz=jz0+1,jz1
            VR=VolC(ix1+1,jy,jz)*FU(ix1,jy,jz)/(FU(ix1,jy,jz)+FU(ix1+1,jy,jz)+Eps)
            VL=VolC(ix1,jy,jz)*FU(ix1,jy,jz)/(FU(ix1,jy,jz)+FU(ix1-1,jy,jz)+Eps)
            uF(ix1,jy,jz)=(VR*uCL(ix1+1,jy,jz,1) &
                         +VL*uCR(ix1,jy,jz,1)) &
                         /(VL+VR+Eps)*FU(ix1,jy,jz)/(FU(ix1,jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='is'.OR.Nachbars(in)%nType=='ps') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            FR= &
               SUM(FV(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))
            vFine= &
               SUM(vCL(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            FC= &
               SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps
            FL= &
               SUM(FV(jx:jx+IncrX-1,iy0-1,jz:jz+IncrZ-1))+Eps
            vCoarse=vCR(jx,iy0,jz,1)
            VL=VolCoarse*FC/(FC+FL+Eps)
            VR=VolFine*FC/(FC+FR+Eps)
            vFace=(vCoarse*VL+vFine*VR)/ &
                  (VL+VR+Eps)
            vF(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)= &
                vFace*FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)/ &
                      (FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            VR=VolC(jx,iy0+1,jz)*FV(jx,iy0,jz)/(FV(jx,iy0,jz)+FV(jx,iy0+1,jz)+Eps)
            VL=VolC(jx,iy0,jz)*FV(jx,iy0,jz)/(FV(jx,iy0,jz)+FV(jx,iy0-1,jz)+Eps)
            vF(jx,iy0,jz)=(VR*vCL(jx,iy0+1,jz,1) &
                          +VL*vCR(jx,iy0,jz,1)) &
                         /(VL+VR+Eps)*FV(jx,iy0,jz)/(FV(jx,iy0,jz)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='in'.OR.Nachbars(in)%nType=='pn') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jx=jx0+1,jx1,IncrX
          DO jz=jz0+1,jz1,IncrZ
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
            FL= &
               SUM(FV(jx:jx+IncrX-1,iy1-1,jz:jz+IncrZ-1))
            vFine= &
               SUM(vCR(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1,1)* &
                   VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            FC= &
               SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))+Eps
            FR= &
               SUM(FV(jx:jx+IncrX-1,iy1+1,jz:jz+IncrZ-1))+Eps
            VL=VolFine*FC/(FC+FL+Eps)
            VR=VolCoarse*FC/(FC+FR+Eps)
            vCoarse=vCL(jx,iy1+1,jz,1)
            vFace=(vCoarse*VR+vFine*VL)/ &
                  (VL+VR+Eps)
            vF(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
               vFace*FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)/ &
                     (FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)+Eps)
          END DO
        END DO
      ELSE
        DO jx=jx0+1,jx1
          DO jz=jz0+1,jz1
            VR=VolC(jx,iy1+1,jz)*FV(jx,iy1,jz)/(FV(jx,iy1,jz)+FV(jx,iy1+1,jz)+Eps)
            VL=VolC(jx,iy1,jz)*FV(jx,iy1,jz)/(FV(jx,iy1,jz)+FV(jx,iy1-1,jz)+Eps)
            vF(jx,iy1,jz)=(VR*vCL(jx,iy1+1,jz,1) &
                          +VL*vCR(jx,iy1,jz,1)) &
                         /(VL+VR+Eps)*FV(jx,iy1,jz)/(FV(jx,iy1,jz)+Eps)
          END DO
        END DO
      END IF
    END IF


    IF (Nachbars(in)%nType=='ib'.OR.Nachbars(in)%nType=='pb') THEN
!     -------------------
!              |
!        N     |    D
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |     |  |  |
!       |     |------
!       |     |  |  |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            FR= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))
            wFine= &
               SUM(wCL(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            FC= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
            FL= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0-1))
            wCoarse=wCR(jx,jy,iz0,1)
            VL=VolCoarse*FC/(FC+FL+Eps)
            VR=VolFine*FC/(FC+FR+Eps)
            wFace=(wCoarse*VL+wFine*VR)/ &
                  (VL+VR+Eps)
            wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)= &
                wFace*FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)/ &
                      (FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            VR=VolC(jx,jy,iz0+1)*FW(jx,jy,iz0)/(FW(jx,jy,iz0)+FW(jx,jy,iz0+1)+Eps)
            VL=VolC(jx,jy,iz0)*FW(jx,jy,iz0)/(FW(jx,jy,iz0)+FW(jx,jy,iz0-1)+Eps)
            wF(jx,jy,iz0)=(VR*wCL(jx,jy,iz0+1,1) &
                          +VL*wCR(jx,jy,iz0,1)) &
                         /(VL+VR+Eps)*FW(jx,jy,iz0)/(FW(jx,jy,iz0)+Eps)
          END DO
        END DO
      END IF
    END IF
    IF (Nachbars(in)%nType=='it'.OR.Nachbars(in)%nType=='pt') THEN
!     -------------------
!              |
!        D     |    N
!              |
!     -------------------
      IF (Refine>RefineNachbar) THEN
!       -------------
!       |  |  |     |
!       |------     |
!       |  |  |     |
!       -------------
        DO jy=jy0+1,jy1,IncrY
          DO jx=jx0+1,jx1,IncrX
            VolFine= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
            FL= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1-1))
            wFine= &
               SUM(wCR(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1,1)* &
                   VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))/ &
               (VolFine+Eps)
            VolCoarse= &
               SUM(VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            FC= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))+Eps
            FR= &
               SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1+1))+Eps
            VL=VolFine*FC/(FC+FL+Eps)
            VR=VolCoarse*FC/(FC+FR+Eps)
            wCoarse=wCL(jx,jy,iz1+1,1)
            wFace=(wCoarse*VR+wFine*VL)/ &
                  (VL+VR+Eps)
            wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
                wFace*FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)/ &
                      (FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)+Eps)
          END DO
        END DO
      ELSE
        DO jy=jy0+1,jy1
          DO jx=jx0+1,jx1
            VR=VolC(jx,jy,iz1+1)*FW(jx,jy,iz1)/(FW(jx,jy,iz1)+FW(jx,jy,iz1+1)+Eps)
            VL=VolC(jx,jy,iz1)*FW(jx,jy,iz1)/(FW(jx,jy,iz1)+FW(jx,jy,iz1-1)+Eps)
            wF(jx,jy,iz1)=(VR*wCL(jx,jy,iz1+1,1) &
                          +VL*wCR(jx,jy,iz1,1)) &
                         /(VL+VR+Eps)*FW(jx,jy,iz1)/(FW(jx,jy,iz1)+Eps)
          END DO
        END DO
      END IF
    END IF
  END DO

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1-1
        VR=VolC(ix+1,iy,iz)*FU(ix,iy,iz)/(FU(ix,iy,iz)+FU(ix+1,iy,iz)+Eps)
        VL=VolC(ix,iy,iz)*FU(ix,iy,iz)/(FU(ix,iy,iz)+FU(ix-1,iy,iz)+Eps)
        uF(ix,iy,iz)=(VR*uCL(ix+1,iy,iz,1) &
                     +VL*uCR(ix,iy,iz,1)) &
                     /(VL+VR+Eps)*FU(ix,iy,iz)/(FU(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1-1
      DO ix=ix0+1,ix1
        VR=VolC(ix,iy+1,iz)*FV(ix,iy,iz)/(FV(ix,iy,iz)+FV(ix,iy+1,iz)+Eps)
        VL=VolC(ix,iy,iz)*FV(ix,iy,iz)/(FV(ix,iy,iz)+FV(ix,iy-1,iz)+Eps)
        vF(ix,iy,iz)=(VR*vCL(ix,iy+1,iz,1) &
                     +VL*vCR(ix,iy,iz,1)) &
                     /(VL+VR+Eps)*FV(ix,iy,iz)/(FV(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1-1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        VR=VolC(ix,iy,iz+1)*FW(ix,iy,iz)/(FW(ix,iy,iz)+FW(ix,iy,iz+1)+Eps)
        VL=VolC(ix,iy,iz)*FW(ix,iy,iz)/(FW(ix,iy,iz)+FW(ix,iy,iz-1)+Eps)
        wF(ix,iy,iz)=(VR*wCL(ix,iy,iz+1,1) &
                     +VL*wCR(ix,iy,iz,1)) &
                     /(VL+VR+Eps)*FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
  IF (TypeW(1:1)=='o') THEN
    ix=ix0
    IF (BCMetVec(uPosL)%West=='Function') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix0)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uF(ix,iy,iz)=RhoFun(xPL,yPL,zPl,Half*(zH(ix,iy)+zh(ix+1,iy)),Time) &
                      *UStart(xPL,yPL,zPl,Half*(zH(ix,iy)+zh(ix+1,iy)),Time) &
                      *FU(ix,iy,iz)/(FU(ix,iy,iz)+eps)-uFOld(ix,iy,iz)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%West=='ZeroGrad' &
         .OR.BCMetVec(uPosL)%West=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix,iy,iz)=uCL(ix+1,iy,iz,1) 
        END DO
      END DO
    ELSE
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix,iy,iz)=Zero
        END DO
      END DO
    END IF
  END IF
  IF (TypeE(1:1)=='o') THEN
    ix=ix1
    IF (BCMetVec(uPosR)%East=='Function') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uF(ix,iy,iz)=RhoFun(xPL,yPL,zPl,Half*(zH(ix,iy)+zh(ix+1,iy)),Time) &
                      *UStart(xPL,yPL,zPl,Half*(zH(ix,iy)+zh(ix+1,iy)),Time) &
                      *FU(ix,iy,iz)/(FU(ix,iy,iz)+eps)-uFOld(ix,iy,iz)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosR)%East=='ZeroGrad' &
         .OR.BCMetVec(uPosR)%East=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix,iy,iz)=uCR(ix,iy,iz,1) 
        END DO
      END DO
    ELSE  
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix,iy,iz)=Zero
        END DO
      END DO
    END IF
  END IF
  IF (TypeS(1:1)=='o') THEN
    iy=iy0
    IF (BCMetVec(vPosL)%South=='Function') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vF(ix,iy,iz)=RhoFun(xPL,yPL,zPl,Half*(zH(ix,iy)+zh(ix,iy+1)),Time) &
                      *VStart(xPL,yPL,zPl,Half*(zH(ix,iy)+zh(ix,iy+1)),Time) &
                      *FV(ix,iy,iz)/(FV(ix,iy,iz)+eps)-vFOld(ix,iy,iz)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%South=='ZeroGrad' &
         .OR.BCMetVec(vPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy,iz)=vCL(ix,iy+1,iz,1) 
        END DO
      END DO
    ELSE  
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy,iz)=Zero
        END DO
      END DO
    END IF
  END IF
  IF (TypeN(1:1)=='o') THEN
    iy=iy1
    IF (BCMetVec(vPosR)%North=='Function') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vF(ix,iy,iz)=RhoFun(xPL,yPL,zPl,Half*(zH(ix,iy)+zh(ix,iy+1)),Time) &
                      *VStart(xPL,yPL,zPl,Half*(zH(ix,iy)+zh(ix,iy+1)),Time) &
                      *FV(ix,iy,iz)/(FV(ix,iy,iz)+eps)-vFOld(ix,iy,iz)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosR)%North=='ZeroGrad' &
         .OR.BCMetVec(vPosR)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy,iz)=vCR(ix,iy,iz,1) 
        END DO
      END DO
    ELSE  
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy,iz)=Zero
        END DO
      END DO
    END IF
  END IF
  IF (TypeB(1:1)=='o') THEN
    iz=iz0
    IF (BCMetVec(wPosL)%Bottom=='Function') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz)
          wF(ix,iy,iz)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time) &
                      *WStart(xPL,yPL,zPl,zH(ix,iy),Time) &
                      *FW(ix,iy,iz)/(FW(ix,iy,iz)+eps)-wFOld(ix,iy,iz)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Bottom=='ZeroGrad' &
         .OR.BCMetVec(wPosL)%Bottom=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz)=wCL(ix,iy,iz+1,1) 
        END DO
      END DO
    ELSE  
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz)=Zero
        END DO
      END DO
    END IF
  END IF
  IF (TypeT(1:1)=='o') THEN
    iz=iz1
    IF (BCMetVec(wPosR)%Top=='Function') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz1)
          wF(ix,iy,iz)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time) &
                      *WStart(xPL,yPL,zPl,zH(ix,iy),Time) &
                      *FW(ix,iy,iz)/(FW(ix,iy,iz)+eps)-wFOld(ix,iy,iz)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosR)%Top=='ZeroGrad' &
         .OR.BCMetVec(wPosR)%Top=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz)=wCR(ix,iy,iz,1) 
        END DO
      END DO
    ELSE  
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz)=Zero
        END DO
      END DO
    END IF
  END IF

END SUBROUTINE VelCellToFaceComputeLR

END MODULE VelocityCellFace_Mod
