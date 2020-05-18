MODULE BoundaryCondition_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Example_Mod
  USE Control_Mod

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

  IF (TypeW=='ow') THEN
    ix=ix0
    IF (BC%West=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          c(ix0,iy,iz,:)=c(ix0+1,iy,iz,:)
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

  IF (TypeE=='oe') THEN
    ix=ix1
    IF (BC%East=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          c(ix1+1,iy,iz,:)=c(ix1,iy,iz,:)
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

  IF (TypeS=='os') THEN
    iy=iy0
    IF (BC%South=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          c(ix,iy0,iz,:)=c(ix,iy0+1,iz,:)
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

  IF (TypeN=='on') THEN
    iy=iy1
    IF (BC%North=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          c(ix,iy1+1,iz,:)=c(ix,iy1,iz,:)
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

  IF (TypeB=='ob') THEN
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

  IF (TypeT=='ot') THEN
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
    IF (BCVec(uPosL)%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix0,iy,iz,1)=-uC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%West=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix0,iy,iz,1)=uC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%West=='MeanValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix0,iy,iz,1)=0.5d0*(MeanProfile(uPosL,iz)+MeanProfile(uPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%West=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix0)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uC(ix0,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix0,iy),Time)*UStart(xPL,yPL,zPl,zH(ix0,iy),Time)
        END DO
      END DO
    END IF

    IF (BCVec(vPosL)%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix0,iy,iz,1)=-vC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%West=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix0,iy,iz,1)=vC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%West=='MeanValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix0,iy,iz,1)=0.5d0*(MeanProfile(vPosL,iz)+MeanProfile(vPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%West=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix0)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vC(ix0,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix0,iy),Time)*VStart(xPL,yPL,zPl,zH(ix0,iy),Time)
        END DO
      END DO
    END IF

    IF (BCVec(wPosL)%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix0,iy,iz,1)=-wC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%West=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix0,iy,iz,1)=wC(ix0+1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%West=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix0,iy,iz,1)=0.5d0*(MeanProfile(wPosL,iz)+MeanProfile(wPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%West=='Function') THEN
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
    IF (BCVec(uPosL)%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix1+1,iy,iz,1)=-uC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%East=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix1+1,iy,iz,1)=uC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%East=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uC(ix1+1,iy,iz,1)=0.5d0*(MeanProfile(uPosL,iz)+MeanProfile(uPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%East=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix1)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uC(ix1+1,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix1+1,iy),Time)*UStart(xPL,yPL,zPl,zH(ix1+1,iy),Time)
        END DO
      END DO
    END IF

    IF (BCVec(vPosL)%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix1+1,iy,iz,1)=-vC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%East=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix1+1,iy,iz,1)=vC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%East=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          vC(ix1+1,iy,iz,1)=0.5d0*(MeanProfile(vPosL,iz)+MeanProfile(vPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%East=='Function') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          xPL=xP(ix1)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vC(ix1+1,iy,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix1+1,iy),Time)*VStart(xPL,yPL,zPl,zH(ix1+1,iy),Time)
        END DO
      END DO
    END IF

    IF (BCVec(wPosL)%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix1+1,iy,iz,1)=-wC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%East=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix1+1,iy,iz,1)=wC(ix1,iy,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%East=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          wC(ix1+1,iy,iz,1)=0.5d0*(MeanProfile(wPosL,iz)+MeanProfile(wPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%East=='Function') THEN
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
    IF (BCVec(uPosL)%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy0,iz,1)=-uC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%South=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy0,iz,1)=uC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy0,iz,1)=0.5d0*(MeanProfile(uPosL,iz)+MeanProfile(uPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%South=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy0)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uC(ix,iy0,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy0),Time)*UStart(xPL,yPL,zPl,zH(ix,iy0),Time)
        END DO
      END DO
    END IF

    IF (BCVec(vPosL)%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy0,iz,1)=-vC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%South=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy0,iz,1)=vC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy0,iz,1)=0.5d0*(MeanProfile(vPosL,iz)+MeanProfile(vPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%South=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy0)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vC(ix,iy0,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy0),Time)*VStart(xPL,yPL,zPl,zH(ix,iy0),Time)
        END DO
      END DO
    END IF

    IF      (BCVec(wPosL)%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy0,iz,1)=-wC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%South=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy0,iz,1)=wC(ix,iy0+1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%South=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy0,iz,1)=0.5d0*(MeanProfile(wPosL,iz)+MeanProfile(wPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%South=='Function') THEN
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
    IF (BCVec(uPosL)%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy1+1,iz,1)=-uC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%North=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy1+1,iz,1)=uC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          uC(ix,iy1+1,iz,1)=0.5d0*(MeanProfile(uPosL,iz)+MeanProfile(uPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%North=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy1)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          uC(ix,iy1+1,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy1+1),Time)*UStart(xPL,yPL,zPl,zH(ix,iy1+1),Time)
        END DO
      END DO
    END IF

    IF (BCVec(vPosL)%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy1+1,iz,1)=-vC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%North=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy1+1,iz,1)=vC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vC(ix,iy1+1,iz,1)=0.5d0*(MeanProfile(vPosL,iz)+MeanProfile(vPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%North=='Function') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy1)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          vC(ix,iy1+1,iz,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy1+1),Time)*VStart(xPL,yPL,zPl,zH(ix,iy1+1),Time)
        END DO
      END DO
    END IF

    IF (BCVec(wPosL)%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy1+1,iz,1)=-wC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%North=='ZeroGrad') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy1+1,iz,1)=wC(ix,iy1,iz,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%North=='MeanValue') THEN 
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          wC(ix,iy1+1,iz,1)=0.5d0*(MeanProfile(wPosL,iz)+MeanProfile(wPosR,iz))
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%North=='Function') THEN
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

    IF (BCVec(uPosL)%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz0,1)=-uC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%Bottom=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz0,1)=uC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%Bottom=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz0,1)=0.5d0*(MeanProfile(uPosL,iz0+1)+MeanProfile(uPosR,iz0+1))
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%Bottom=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz0)
          uC(ix,iy,iz0,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*UStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

    IF (BCVec(vPosL)%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz0,1)=-vC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%Bottom=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz0,1)=vC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%Bottom=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz0,1)=0.5d0*(MeanProfile(vPosL,iz0+1)+MeanProfile(vPosR,iz0+1))
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%Bottom=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz0)
          vC(ix,iy,iz0,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*VStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

    IF (BCVec(wPosL)%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz0,1)=-wC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%Bottom=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz0,1)=wC(ix,iy,iz0+1,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%Bottom=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz0,1)=0.5d0*(MeanProfile(wPosL,iz0+1)+MeanProfile(wPosR,iz0+1))
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%Bottom=='Function') THEN
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
    IF (BCVec(uPosL)%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz1+1,1)=-uC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%Top=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz1+1,1)=uC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%Top=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          uC(ix,iy,iz1+1,1)=0.5d0*(MeanProfile(uPosL,iz1)+MeanProfile(uPosR,iz1))
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%Top=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz1)
          uC(ix,iy,iz1+1,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*UStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

    IF (BCVec(vPosL)%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz1+1,1)=-vC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%Top=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz1+1,1)=vC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%Top=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          vC(ix,iy,iz1+1,1)=0.5d0*(MeanProfile(vPosL,iz1)+MeanProfile(vPosR,iz1))
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%Top=='Function') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz1)
          vC(ix,iy,iz1+1,1)=RhoFun(xPL,yPL,zPl,zH(ix,iy),Time)*VStart(xPL,yPL,zPl,zH(ix,iy),Time)
        END DO
      END DO
    END IF

    IF (BCVec(wPosL)%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz1+1,1)=-wC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%Top=='ZeroGrad') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz1+1,1)=wC(ix,iy,iz1,1)
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%Top=='MeanValue') THEN 
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wC(ix,iy,iz1+1,1)=0.5d0*(MeanProfile(wPosL,iz1)+MeanProfile(wPosR,iz1))
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%Top=='Function') THEN
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
    IF (BCVec(uPosL)%West=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix0,iy,iz)=0.d0
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%West=='Function') THEN
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
    IF (BCVec(uPosL)%East=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          uF(ix1,iy,iz)=0.d0
        END DO
      END DO
    ELSE IF (BCVec(uPosL)%East=='Function') THEN
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
    IF (BCVec(vPosL)%South=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy0,iz)=0.d0
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%South=='Function') THEN
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
    IF (BCVec(vPosL)%North=='ZeroValue') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          vF(ix,iy1,iz)=0.d0
        END DO
      END DO
    ELSE IF (BCVec(vPosL)%North=='Function') THEN
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
    IF (BCVec(wPosL)%Bottom=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz0)=0.d0
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%Bottom=='Function') THEN
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
    IF (BCVec(wPosL)%Top=='ZeroValue') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          wF(ix,iy,iz1)=0.d0
        END DO
      END DO
    ELSE IF (BCVec(wPosL)%Top=='Function') THEN
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
 
SUBROUTINE BoundaryCondition(VectorCell,Time)
 
  TYPE(Vector4Cell_T), TARGET :: VectorCell(:)
  REAL(RealKind) :: Time
 
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    DO ic=LBOUND(VectorCell(ibLoc)%Vec,1)+1,UBOUND(VectorCell(ibLoc)%Vec,1)
      c=>VectorCell(ibLoc)%Vec(ic)%c
      BC=>BCVec(ic)
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
  END DO

END SUBROUTINE BoundaryCondition

END MODULE BoundaryCondition_Mod
