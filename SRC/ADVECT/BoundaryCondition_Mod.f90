MODULE BoundaryCondition_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Example_Mod
  USE Control_Mod
  USE Chemie_Mod

  IMPLICIT NONE

  INTERFACE BoundaryVelocity
    MODULE PROCEDURE BoundaryVelocityFace,BoundaryVelocityCell
  END INTERFACE

CONTAINS

SUBROUTINE BoundaryComputeZeroGrad

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
          c(ix,iy,iz0,:)=Zero
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

END SUBROUTINE BoundaryComputeZeroGrad


SUBROUTINE BoundaryCompute(Val,Time,ic) 

  REAL(RealKind), EXTERNAL :: Val
  REAL(RealKind) :: Time

  INTEGER :: ix,iy,iz,ic
  REAL(RealKind) :: xPL,yPL,zPl

  IF (TypeW(1:2)=='ow'.OR.TypeW=='iwo') THEN
    ix=ix0
    IF (BC%West=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          c(ix0,iy,iz,:)=c(ix0+1,iy,iz,:)
        END DO
      END DO
    ELSE IF (BC%West=='InFlowOutFlow') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          IF (uF(ix0,iy,iz)>=0.0d0) THEN
            c(ix0,iy,iz,:)=Rho(ix0,iy,iz,1)*VecAmb(ibLoc)%Vec(ic)%c(ix0,iy,iz,:)
          ELSE
            c(ix0,iy,iz,:)=c(ix0+1,iy,iz,:)
          END IF
        END DO
      END DO
    ELSE IF (BC%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          c(ix0,iy,iz,:)=-c(ix0+1,iy,iz,:)
        END DO
      END DO
    ELSE IF (BC%West=='MeanValue') THEN 
       DO iz=iz0+1,iz1
         DO iy=iy0+1,iy1
           c(ix0,iy,iz,:) = 2.d0*MeanProfile(ic,iz)-c(ix0+1,iy,iz,:)
         END DO
       END DO
    ELSE IF (BC%West=='Function') THEN
       DO iz=iz0+1,iz1
         DO iy=iy0+1,iy1
           xPL=xP(ix0)-0.5d0*dx(ix0+1)
           yPL=yP(iy-1)+0.5e0*dy(iy)
           zPL=zP(iz-1)+0.5e0*dz(iz)
           IF (ic==nvPos.OR.ic==ncPos.OR.ic==nrPos.OR.ic==nsPos.OR.ic==niPos) THEN 
             c(ix0,iy,iz,:)=2.d0*Val(xPL,yPL,zPL,zH(ix0,iy),Time)-c(ix0+1,iy,iz,:)
           ELSE
             c(ix0,iy,iz,:)=2.d0*RhoFun(xPL,yPL,zPL,zH(ix0,iy),Time)*Val(xPL,yPL,zPL,zH(ix0,iy),Time)-c(ix0+1,iy,iz,:)
           END IF
         END DO
       END DO
    END IF
  END IF

  IF (TypeE(1:2)=='oe'.OR.TypeE=='ieo') THEN
    ix=ix1
    IF (BC%East=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          c(ix1+1,iy,iz,:)=c(ix1,iy,iz,:)
        END DO
      END DO
    ELSE IF (BC%East=='InFlowOutFlow') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          IF (uF(ix1,iy,iz)<=0.0d0) THEN
            c(ix1+1,iy,iz,:)=Rho(ix1+1,iy,iz,1)*VecAmb(ibLoc)%Vec(ic)%c(ix1+1,iy,iz,:)
          ELSE
            c(ix1+1,iy,iz,:)=c(ix1,iy,iz,:)
          END IF
        END DO
      END DO
    ELSE IF (BC%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          c(ix1+1,iy,iz,:)=-c(ix1,iy,iz,:)
        END DO
      END DO
    ELSE IF (BC%East=='MeanValue') THEN 
       DO iz=iz0+1,iz1
         DO iy=iy0+1,iy1
           c(ix1+1,iy,iz,:) = 2.d0*MeanProfile(ic,iz)-c(ix1,iy,iz,:)
         END DO
       END DO
    ELSE IF (BC%East=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix1)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          IF (ic==nvPos.OR.ic==ncPos.OR.ic==nrPos.OR.ic==nsPos.OR.ic==niPos) THEN 
            c(ix1+1,iy,iz,:)=2.d0*Val(xPL,yPL,zPL,zH(ix1+1,iy),Time)-c(ix1,iy,iz,:)
          ELSE
            c(ix1+1,iy,iz,:)=2.d0*RhoFun(xPL,yPL,zPL,zH(ix1+1,iy),Time)*Val(xPL,yPL,zPL,zH(ix1+1,iy),Time)-c(ix1,iy,iz,:)
          END IF
        END DO
      END DO
    END IF
  END IF

  IF (TypeS(1:2)=='os'.OR.TypeS=='iso') THEN
    iy=iy0
    IF (BC%South=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          c(ix,iy0,iz,:)=c(ix,iy0+1,iz,:)
        END DO
      END DO
    ELSE IF (BC%South=='InFlowOutFlow') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          IF (vF(ix,iy0,iz)>=0.0d0) THEN
            c(ix,iy0,iz,:)=Rho(ix,iy0,iz,1)*VecAmb(ibLoc)%Vec(ic)%c(ix,iy0,iz,:)
          ELSE
            c(ix,iy0,iz,:)=c(ix,iy0+1,iz,:)
          END IF
        END DO
      END DO
    ELSE IF (BC%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          c(ix,iy0,iz,:)=-c(ix,iy0+1,iz,:)
        END DO
      END DO
    ELSE IF (BC%South=='MeanValue') THEN 
       DO iz=iz0+1,iz1
         DO ix=ix0+1,ix1
           c(ix,iy0,iz,:) = 2.d0*MeanProfile(ic,iz)-c(ix,iy0+1,iz,:)
         END DO
       END DO
    ELSE IF (BC%South=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy0)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          IF (ic==nvPos.OR.ic==ncPos.OR.ic==nrPos.OR.ic==nsPos.OR.ic==niPos) THEN 
            c(ix,iy0,iz,:)=2.d0*Val(xPL,yPL,zPL,zH(ix,iy0),Time)-c(ix,iy0+1,iz,:)
          ELSE
            c(ix,iy0,iz,:)=2.d0*RhoFun(xPL,yPL,zPL,zH(ix,iy0),Time)*Val(xPL,yPL,zPL,zH(ix,iy0),Time)-c(ix,iy0+1,iz,:)
          END IF
        END DO
      END DO
    END IF
  END IF

  IF (TypeN(1:2)=='on'.OR.TypeN=='ino') THEN
    iy=iy1
    IF (BC%North=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          c(ix,iy1+1,iz,:)=c(ix,iy1,iz,:)
        END DO
      END DO
    ELSE IF (BC%North=='InFlowOutFlow') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          IF (vF(ix,iy1,iz)<=0.0d0) THEN
            c(ix,iy1+1,iz,:)=Rho(ix,iy1+1,iz,1)*VecAmb(ibLoc)%Vec(ic)%c(ix,iy1+1,iz,:)
          ELSE
            c(ix,iy1+1,iz,:)=c(ix,iy1,iz,:)
          END IF
        END DO
      END DO
    ELSE IF (BC%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          c(ix,iy1+1,iz,:)=-c(ix,iy1,iz,:)
        END DO
      END DO
    ELSE IF (BC%North=='MeanValue') THEN 
       DO iz=iz0+1,iz1
         DO ix=ix0+1,ix1
           c(ix,iy1+1,iz,:) = 2.d0*MeanProfile(ic,iz)-c(ix,iy1,iz,:)
         END DO
       END DO
    ELSE IF (BC%North=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy1)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          IF (ic==nvPos.OR.ic==ncPos.OR.ic==nrPos.OR.ic==nsPos.OR.ic==niPos) THEN 
            c(ix,iy1+1,iz,:)=2.d0*Val(xPL,yPL,zPL,zH(ix,iy1+1),Time)-c(ix,iy1,iz,:)
          ELSE
            c(ix,iy1+1,iz,:)=2.d0*RhoFun(xPL,yPL,zPL,zH(ix,iy1+1),Time)*Val(xPL,yPL,zPL,zH(ix,iy1+1),Time)-c(ix,iy1,iz,:)
          END IF
        END DO
      END DO
    END IF
  END IF

  IF (TypeB(1:2)=='ob'.OR.TypeB=='ibo') THEN
    iz=iz0
    IF (BC%Bottom=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          c(ix,iy,iz0,:)=c(ix,iy,iz0+1,:)
        END DO
      END DO
    ELSE IF (BC%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          c(ix,iy,iz0,:)=-c(ix,iy,iz0+1,:)
        END DO
      END DO
    ELSE IF (BC%Bottom=='Extrapolation') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          c(ix,iy,iz0,:)=Two*c(ix,iy,iz0+1,:)-c(ix,iy,iz0+2,:)
        END DO
      END DO
    ELSE IF (BC%Bottom=='MeanValue') THEN 
       DO iy=iy0+1,iy1
         DO ix=ix0+1,ix1
           c(ix,iy,iz0,:) = 2.d0*MeanProfile(ic,iz0+1)-c(ix,iy,iz0+1,:)
         END DO
       END DO
    ELSE IF (BC%Bottom=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz0)
          IF (ic==nvPos.OR.ic==ncPos.OR.ic==nrPos.OR.ic==nsPos.OR.ic==niPos) THEN 
            c(ix,iy,iz0,:)=2.d0*Val(xPL,yPL,zPL,zH(ix,iy),Time)-c(ix,iy,iz0+1,:)
          ELSE
            c(ix,iy,iz0,:)=2.d0*RhoFun(xPL,yPL,zPL,zH(ix,iy),Time)*Val(xPL,yPL,zPL,zH(ix,iy),Time)-c(ix,iy,iz0+1,:)
          END IF
        END DO
      END DO
    END IF
  END IF

  IF (TypeT(1:2)=='ot'.OR.TypeT=='ito') THEN
    iz=iz1
    IF (BC%Top=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          c(ix,iy,iz1+1,:)=c(ix,iy,iz1,:)
        END DO
      END DO
    ELSE IF (BC%Top=='Extrapolation') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          c(ix,iy,iz1+1,:)=Two*c(ix,iy,iz1,:)-c(ix,iy,iz1-1,:)
        END DO
      END DO
    ELSE IF (BC%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          c(ix,iy,iz1+1,:)=-c(ix,iy,iz1,:)
        END DO
      END DO
    ELSE IF (BC%Top=='MeanValue') THEN 
       DO iy=iy0+1,iy1
         DO ix=ix0+1,ix1
           c(ix,iy,iz1+1,:) = 2.d0*MeanProfile(ic,iz1)-c(ix,iy,iz1,:)
         END DO
       END DO
    ELSE IF (BC%Top=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz1)
          IF (ic==nvPos.OR.ic==ncPos.OR.ic==nrPos.OR.ic==nsPos.OR.ic==niPos) THEN 
            c(ix,iy,iz1+1,:)=2.d0*Val(xPL,yPL,zPL,zH(ix,iy),Time)-c(ix,iy,iz1,:)
          ELSE
            c(ix,iy,iz1+1,:)=2.d0*RhoFun(xPL,yPL,zPL,zH(ix,iy),Time)*Val(xPL,yPL,zPL,zH(ix,iy),Time)-c(ix,iy,iz1,:)
          END IF
        END DO
      END DO
    END IF
  END IF

END SUBROUTINE BoundaryCompute

SUBROUTINE BoundaryFluxCompute(ic) 
  INTEGER :: ic

  IF (TypeW(1:2)=='ow'.OR.TypeW=='iwo') THEN
    IF (BC%West=='InFlowOutFlow') THEN
      c(ix0,:,:,:)=0.0d0
    END IF
  END IF

  IF (TypeE(1:2)=='oe'.OR.TypeE=='ieo') THEN
    IF (BC%East=='InFlowOutFlow') THEN
      c(ix1+1,:,:,:)=0.0d0
    END IF
  END IF

  IF (TypeS(1:2)=='os'.OR.TypeS=='iso') THEN
    IF (BC%South=='InFlowOutFlow') THEN
      c(:,iy0,:,:)=0.0d0
    END IF
  END IF

  IF (TypeN(1:2)=='on'.OR.TypeN=='ino') THEN
    IF (BC%North=='InFlowOutFlow') THEN
      c(:,iy1+1,:,:)=0.0d0
    END IF
  END IF

  IF (TypeB(1:2)=='ob'.OR.TypeB=='ibo') THEN
    IF (BC%Bottom=='InFlowOutFlow') THEN
      c(:,:,iz0,:)=0.0d0
    END IF
  END IF

  IF (TypeT(1:2)=='ot'.OR.TypeT=='ito') THEN
    IF (BC%Top=='InFlowOutFlow') THEN
      c(:,:,iz1+1,:)=0.0d0
    END IF
  END IF

END SUBROUTINE BoundaryFluxCompute

SUBROUTINE BoundaryVelocityFace(VelocityFace,Time)
 
  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)
  REAL(RealKind) :: Time

  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    CALL BoundaryVelocityFaceCompute(Time)
  END DO

END SUBROUTINE BoundaryVelocityFace

SUBROUTINE BoundaryVelocityCell(VelocityCell,Time)
 
  TYPE(Vector4Cell_T), TARGET :: VelocityCell(:)
  REAL(RealKind) :: Time

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    uC=>VelocityCell(ibLoc)%Vec(uPosL)%c
    vC=>VelocityCell(ibLoc)%Vec(vPosL)%c
    wC=>VelocityCell(ibLoc)%Vec(wPosL)%c
    CALL BoundaryVelocityCellCompute(Time)
    uC=>VelocityCell(ibLoc)%Vec(uPosR)%c
    vC=>VelocityCell(ibLoc)%Vec(vPosR)%c
    wC=>VelocityCell(ibLoc)%Vec(wPosR)%c
    CALL BoundaryVelocityCellCompute(Time)
  END DO

END SUBROUTINE BoundaryVelocityCell

SUBROUTINE BoundaryVelocityCellCompute(Time) 

  REAL(RealKind) :: Time 

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: xPL,yPL,zPl

  IF (TypeW=='ow') THEN
    IF (BCMetVec(uPosL)%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix0,iy,iz,1)=-uC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%West=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix0,iy,iz,1)=uC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%West=='MeanValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix0,iy,iz,1)=0.5d0*(MeanProfile(uPosL,iz)+MeanProfile(uPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%West=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix0)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uC(ix0,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix0,iy),Time)*UStart(xPL,yPL,zPl,zH(ix0,iy),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(vPosL)%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix0,iy,iz,1)=-vC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%West=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix0,iy,iz,1)=vC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%West=='MeanValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix0,iy,iz,1)=0.5d0*(MeanProfile(vPosL,iz)+MeanProfile(vPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%West=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix0)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vC(ix0,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix0,iy),Time)*VStart(xPL,yPL,zPl,zH(ix0,iy),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(wPosL)%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix0,iy,iz,1)=-wC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%West=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix0,iy,iz,1)=wC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%West=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix0,iy,iz,1)=0.5d0*(MeanProfile(wPosL,iz)+MeanProfile(wPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%West=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix0)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          wC(ix0,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix0,iy),Time)*WStart(xPL,yPL,zPl,zH(ix0,iy),Time)
        END DO
      END DO
    END IF
  END IF

  IF (TypeE=='oe') THEN
    IF (BCMetVec(uPosL)%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix1+1,iy,iz,1)=-uC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%East=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix1+1,iy,iz,1)=uC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%East=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix1+1,iy,iz,1)=0.5d0*(MeanProfile(uPosL,iz)+MeanProfile(uPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%East=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix1)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uC(ix1+1,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix1+1,iy),Time)*UStart(xPL,yPL,zPl,zH(ix1+1,iy),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(vPosL)%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix1+1,iy,iz,1)=-vC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%East=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix1+1,iy,iz,1)=vC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%East=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix1+1,iy,iz,1)=0.5d0*(MeanProfile(vPosL,iz)+MeanProfile(vPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%East=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix1)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vC(ix1+1,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix1+1,iy),Time)*VStart(xPL,yPL,zPl,zH(ix1+1,iy),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(wPosL)%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix1+1,iy,iz,1)=-wC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%East=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix1+1,iy,iz,1)=wC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%East=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix1+1,iy,iz,1)=0.5d0*(MeanProfile(wPosL,iz)+MeanProfile(wPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%East=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix1)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          wC(ix1+1,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix1+1,iy),Time)*WStart(xPL,yPL,zPl,zH(ix1+1,iy),Time)
        END DO
      END DO
    END IF

  END IF

  IF (TypeS=='os') THEN
    IF (BCMetVec(uPosL)%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy0,iz,1)=-uC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%South=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy0,iz,1)=uC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy0,iz,1)=0.5d0*(MeanProfile(uPosL,iz)+MeanProfile(uPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%South=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy0)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uC(ix,iy0,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy0),Time)*UStart(xPL,yPL,zPl,zH(ix,iy0),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(vPosL)%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy0,iz,1)=-vC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%South=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy0,iz,1)=vC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy0,iz,1)=0.5d0*(MeanProfile(vPosL,iz)+MeanProfile(vPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%South=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy0)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vC(ix,iy0,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy0),Time)*VStart(xPL,yPL,zPl,zH(ix,iy0),Time)
        END DO
      END DO
    END IF

    IF      (BCMetVec(wPosL)%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy0,iz,1)=-wC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%South=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy0,iz,1)=wC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy0,iz,1)=0.5d0*(MeanProfile(wPosL,iz)+MeanProfile(wPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%South=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy0)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          wC(ix,iy0,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy0),Time)*WStart(xPL,yPL,zPl,zH(ix,iy0),Time)
        END DO
      END DO
    END IF

  END IF

  IF (TypeN=='on') THEN
    IF (BCMetVec(uPosL)%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy1+1,iz,1)=-uC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%North=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy1+1,iz,1)=uC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy1+1,iz,1)=0.5d0*(MeanProfile(uPosL,iz)+MeanProfile(uPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%North=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy1)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uC(ix,iy1+1,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy1+1),Time)*UStart(xPL,yPL,zPl,zH(ix,iy1+1),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(vPosL)%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy1+1,iz,1)=-vC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%North=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy1+1,iz,1)=vC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy1+1,iz,1)=0.5d0*(MeanProfile(vPosL,iz)+MeanProfile(vPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%North=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy1)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vC(ix,iy1+1,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy1+1),Time)*VStart(xPL,yPL,zPl,zH(ix,iy1+1),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(wPosL)%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy1+1,iz,1)=-wC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%North=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy1+1,iz,1)=wC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy1+1,iz,1)=0.5d0*(MeanProfile(wPosL,iz)+MeanProfile(wPosR,iz))
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%North=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy1)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          wC(ix,iy1+1,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy1+1),Time)*WStart(xPL,yPL,zPl,zH(ix,iy1+1),Time)
        END DO
      END DO
    END IF
  END IF

  IF (TypeB=='ob') THEN

    IF (BCMetVec(uPosL)%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz0,1)=-uC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%Bottom=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz0,1)=uC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%Bottom=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz0,1)=0.5d0*(MeanProfile(uPosL,iz0+1)+MeanProfile(uPosR,iz0+1))
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%Bottom=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz0)
          uC(ix,iy,iz0,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*UStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(vPosL)%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz0,1)=-vC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%Bottom=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz0,1)=vC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%Bottom=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz0,1)=0.5d0*(MeanProfile(vPosL,iz0+1)+MeanProfile(vPosR,iz0+1))
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%Bottom=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz0)
          vC(ix,iy,iz0,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*VStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(wPosL)%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz0,1)=-wC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Bottom=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz0,1)=wC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Bottom=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz0,1)=0.5d0*(MeanProfile(wPosL,iz0+1)+MeanProfile(wPosR,iz0+1))
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Bottom=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz0)
          wC(ix,iy,iz0,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*WStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

  END IF

  IF (TypeT=='ot') THEN
    IF (BCMetVec(uPosL)%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz1+1,1)=-uC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%Top=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz1+1,1)=uC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%Top=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz1+1,1)=0.5d0*(MeanProfile(uPosL,iz1)+MeanProfile(uPosR,iz1))
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%Top=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz1)
          uC(ix,iy,iz1+1,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*UStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(vPosL)%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz1+1,1)=-vC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%Top=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz1+1,1)=vC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%Top=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz1+1,1)=0.5d0*(MeanProfile(vPosL,iz1)+MeanProfile(vPosR,iz1))
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%Top=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz1)
          vC(ix,iy,iz1+1,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*VStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

    IF (BCMetVec(wPosL)%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz1+1,1)=-wC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Top=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz1+1,1)=wC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Top=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz1+1,1)=0.5d0*(MeanProfile(wPosL,iz1)+MeanProfile(wPosR,iz1))
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Top=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz1)
          wC(ix,iy,iz1+1,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*WStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

  END IF

END SUBROUTINE BoundaryVelocityCellCompute

SUBROUTINE BoundaryVelocityFaceCompute(Time)

  REAL(RealKind) :: Time

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: xPL,yPL,zPl

  IF (TypeW=='ow') THEN
    ix=ix0
    IF (BCMetVec(uPosL)%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix0,iy,iz)=0.d0
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%West=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix0)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uF(ix0,iy,iz)=RhoFun(xPL,yPL,zPl,Half*(zH(ix0,iy)+zh(ix0+1,iy)),Time) &
                        *UStart(xPL,yPL,zPl,Half*(zH(ix0,iy)+zh(ix0+1,iy)),Time) &
                        *FU(ix0,iy,iz)/(FU(ix0,iy,iz)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeE=='oe') THEN
    ix=ix1
    IF (BCMetVec(uPosL)%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix1,iy,iz)=0.d0
        END DO
      END DO
    ELSE IF (BCMetVec(uPosL)%East=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix1)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uF(ix1,iy,iz)=RhoFun(xPL,yPL,zPl,Half*(zH(ix1,iy)+zh(ix1+1,iy)),Time) &
                       *UStart(xPL,yPL,zPl,Half*(zH(ix1,iy)+zh(ix1+1,iy)),Time) &
                       *FU(ix1,iy,iz)/(FU(ix1,iy,iz)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeS=='os') THEN
    iy=iy0
    IF (BCMetVec(vPosL)%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy0,iz)=0.d0
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%South=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy0)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vF(ix,iy0,iz) = RhoFun(xPL,yPL,zPl,Half*(zH(ix,iy0)+zh(ix,iy0+1)),Time)*&
          &VStart(xPL,yPL,zPl,Half*(zH(ix,iy0)+zh(ix,iy0+1)),Time) &
                                *FV(ix,iy0,iz)/(FV(ix,iy0,iz)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeN=='on') THEN
    iy=iy1
    IF (BCMetVec(vPosL)%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy1,iz)=0.d0
        END DO
      END DO
    ELSE IF (BCMetVec(vPosL)%North=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy1)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vF(ix,iy1,iz) = RhoFun(xPL,yPL,zPl,Half*(zH(ix,iy1)+zh(ix,iy1+1)),Time)*&
          &VStart(xPL,yPL,zPl,Half*(zH(ix,iy1)+zh(ix,iy1+1)),Time) &
                                *FV(ix,iy1,iz)/(FV(ix,iy1,iz)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeB=='ob') THEN
    iz=iz0
    IF (BCMetVec(wPosL)%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz0)=0.d0
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Bottom=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz0)
          wF(ix,iy,iz0) = RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*WStart(xPL,yPL,zPl,zH(ix,iy),Time) &
                                *FW(ix,iy,iz0)/(FW(ix,iy,iz0)+eps)
        END DO
      END DO
    END IF
  END IF

  IF (TypeT=='ot') THEN
    iz=iz1
    IF (BCMetVec(wPosL)%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz1)=0.d0
        END DO
      END DO
    ELSE IF (BCMetVec(wPosL)%Top=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz1)
          wF(ix,iy,iz1) = RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*WStart(xPL,yPL,zPl,zH(ix,iy),Time) &
                                *FW(ix,iy,iz1)/(FW(ix,iy,iz1)+eps)
        END DO
      END DO
    END IF
  END IF

END SUBROUTINE BoundaryVelocityFaceCompute

SUBROUTINE BoundaryConditionT(VectorCell,Time)
 
  TYPE(Vector4Cell_T), TARGET :: VectorCell(:)
  REAL(RealKind) :: Time
 
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    DO ic=1,UBOUND(VectorCell(ibLoc)%Vec,1)
      c=>VectorCell(ibLoc)%Vec(ic)%c
      CALL BoundaryComputeZeroGrad
    END DO
  END DO
END SUBROUTINE BoundaryConditionT

SUBROUTINE BoundaryFluxCondition(VectorCellMet,VectorCellChem)
 
  TYPE(Vector4Cell_T), TARGET :: VectorCellMet(:)
  TYPE(Vector4Cell_T), TARGET :: VectorCellChem(:)
 
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    DO ic=LBOUND(VectorCellChem(ibLoc)%Vec,1)+1,UBOUND(VectorCellChem(ibLoc)%Vec,1)
      c=>VectorCellChem(ibLoc)%Vec(ic)%c
      BC=>BCChemVec(ic)
      CALL BoundaryFluxCompute(ic)
    END DO
  END DO

END SUBROUTINE BoundaryFluxCondition

SUBROUTINE BoundaryCondition(VectorCellMet,VectorCellChem,VelF,Time)
 
  TYPE(Vector4Cell_T), TARGET :: VectorCellMet(:)
  TYPE(Vector4Cell_T), OPTIONAL, TARGET :: VectorCellChem(:)
  TYPE (VelocityFace_T), TARGET :: VelF(:)
  REAL(RealKind) :: Time
 
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    Rho=>VectorCellMet(ibLoc)%Vec(RhoPos)%c
    uF=>VelF(ibLoc)%uF
    vF=>VelF(ibLoc)%vF
    wF=>VelF(ibLoc)%wF
    DO ic=LBOUND(VectorCellMet(ibLoc)%Vec,1)+1,UBOUND(VectorCellMet(ibLoc)%Vec,1)
      c=>VectorCellMet(ibLoc)%Vec(ic)%c
      BC=>BCMetVec(ic)
      IF      (ic==  uPosL) THEN
      ELSE IF (ic==  vPosL) THEN
      ELSE IF (ic==  wPosL) THEN
      ELSE IF (ic==  uPosR) THEN
      ELSE IF (ic==  vPosR) THEN
      ELSE IF (ic==  wPosR) THEN
      ELSE IF (ic== thPos) THEN
        CALL BoundaryCompute(   ThStart,Time,ic)
      ELSE IF (ic== EnPos) THEN
        CALL BoundaryCompute(   EnStart,Time,ic)
      ELSE IF (ic== RhoPos) THEN
        CALL BoundaryCompute(   RhoStart,Time,ic)
      ELSE IF (ic==tkePos) THEN
        CALL BoundaryCompute(  TkeStart,Time,ic)
      ELSE IF (ic==disPos) THEN
        CALL BoundaryCompute(  DisStart,Time,ic)
      ELSE IF (ic==omePos) THEN
        CALL BoundaryCompute(  OmeStart,Time,ic)
      ELSE IF (ic==tkeHPos) THEN
        CALL BoundaryCompute(  TkeHStart,Time,ic)
      ELSE IF (ic==tkeVPos) THEN
        CALL BoundaryCompute(  TkeVStart,Time,ic)
      ELSE IF (ic==LenPos) THEN
        CALL BoundaryCompute(  LenStart,Time,ic)
      ELSE IF (ic== RhoVPos) THEN
        CALL BoundaryCompute(   QvStart,Time,ic)
      ELSE IF (ic== RhoCPos) THEN
        CALL BoundaryCompute(   QcStart,Time,ic)
      ELSE IF (ic== RhoRPos) THEN
        CALL BoundaryCompute(   QrStart,Time,ic)
      ELSE IF (ic== RhoIPos) THEN
        CALL BoundaryCompute(   QiStart,Time,ic)
      ELSE IF (ic== RhoSPos) THEN
        CALL BoundaryCompute(   QsStart,Time,ic)
      ELSE IF (ic== nvPos) THEN
        CALL BoundaryCompute(   NvStart,Time,ic)
      ELSE IF (ic== ncPos) THEN
        CALL BoundaryCompute(   NcStart,Time,ic)
      ELSE IF (ic== nrPos) THEN
        CALL BoundaryCompute(   NrStart,Time,ic)
      ELSE IF (ic== niPos) THEN
        CALL BoundaryCompute(   NiStart,Time,ic)
      ELSE IF (ic== nsPos) THEN
        CALL BoundaryCompute(   NsStart,Time,ic)
      ELSE IF (ic== tracer1Pos) THEN
        CALL BoundaryCompute(Tracer1Start,Time,ic)
      ELSE IF (ic== tracer2Pos) THEN
        CALL BoundaryCompute(Tracer2Start,Time,ic)
      ELSE
        CALL BoundaryCompute(DummyStart,Time,ic)
      END IF
    END DO
    IF (PRESENT(VectorCellChem)) THEN
      DO ic=LBOUND(VectorCellChem(ibLoc)%Vec,1)+1,UBOUND(VectorCellChem(ibLoc)%Vec,1)
        c=>VectorCellChem(ibLoc)%Vec(ic)%c
        BC=>BCChemVec(ic)
        CALL BoundaryCompute(DummyStart,Time,ic)
      END DO
    END IF
  END DO

END SUBROUTINE BoundaryCondition

END MODULE BoundaryCondition_Mod
