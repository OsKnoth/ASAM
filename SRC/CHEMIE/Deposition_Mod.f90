MODULE Deposition_Mod

  USE EmissDeposParameter_Mod

  IMPLICIT NONE

  TYPE(ScalarCell_T), POINTER :: SediVelC


CONTAINS

SUBROUTINE Deposition(VelF)

  TYPE(VelocityFace_T) :: VelF(:)

  INTEGER :: i,j,iCell,is
  INTEGER :: ix,iy,iz
  INTEGER :: LandUse
  REAL(RealKind) :: n1,n2,n3,FL
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: Schmidt,Brown,Stoke,Impact,Interception,Diff
  REAL(RealKind) :: Vdepos,Rb,Ra,R1,Rc,Kw
  REAL(RealKind) :: qv,RhoPart
  REAL(RealKind) :: DragH,DragM,ustern
  REAL(RealKind) :: U10,TSurface
  REAL(RealKind) :: nC,Rad
  REAL(RealKind) :: MM(0:nFrac+1,nAqua)
  REAL(RealKind) :: mFrac(nAqua)
  REAL(RealKind), POINTER :: VelSedi(:,:,:,:), RhoLoc(:,:,:,:)
  REAL(RealKind), POINTER :: uF(:,:,:)
  REAL(RealKind), POINTER :: vF(:,:,:)
  REAL(RealKind), POINTER :: wF(:,:,:)
  REAL(RealKind) :: VC,VL
  REAL(RealKind) :: vFall
  REAL(RealKind) :: RH
     
  uF=>VelF(ibLoc)%uF
  vF=>VelF(ibLoc)%vF
  wF=>VelF(ibLoc)%wF
  VelSedi=>SediVelC%c
  RhoLoc=>cVec(RhoPos)%c

  DO iCell=1,NumBoundCell
    ix   = BoundCell(iCell)%ix
    iy   = BoundCell(iCell)%iy
    iz   = BoundCell(iCell)%iz
    n1   = BoundCell(iCell)%n1
    n2   = BoundCell(iCell)%n2
    n3   = BoundCell(iCell)%n3
    FL   = BoundCell(iCell)%FL
    DragH= BoundCell(iCell)%DragH
    DragM= BoundCell(iCell)%DragM
    TSurface = BoundCell(iCell)%TeS
    U10 = BoundCell(iCell)%U10
    LandUse = BoundCell(iCell)%LandClass
    qv=BoundCell(iCell)%qv

    RhoPart=1.2d3 ! kg/m³     ! Aus ASAM
    IF (LandUse.EQ.9) THEN
      RH=1.d0
    ELSE
      RH=qv*RV*TSurface/SaturVapor(TSurface) ! Calculation of rel hum at surface
    END IF

!   Evaluate velocity tangential and normal  
    uLoc = &
           (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
            FU(ix,iy,iz)*uF(ix,iy,iz)) / &
          ((FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*(RhoLoc(ix,iy,iz,1)+Eps))
    vLoc = &
          (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
           FV(ix,iy,iz)*vF(ix,iy,iz)) / &
         ((FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)*(RhoLoc(ix,iy,iz,1)+Eps))
    wLoc = &
          (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
           FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
         ((FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)*(RhoLoc(ix,iy,iz,1)+Eps))

    VN   = uLoc*n1+vLoc*n2+wLoc*n3
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero)) ! Wind

    VC=VolC(ix,iy,iz)
    VL=VolC(ix,iy,iz-1)

    ustern=SQRT(DragM)*VT

    Ra=1/(VT*DragH) ! Turbulent Resistenz

    IF (Aerosol) THEN
      DO i=1,nFrac
        nC=cVec(iNC)%c(ix,iy,iz,i)
        IF (nC.GT.0.0d0) THEN
          DO is=1,nAqua
            MM(i,is)=cVec(is)%c(ix,iy,iz,i)
          END DO
          mFrac=MM(i,:)/nC
          iF=iz
          Rad=Radius(mFrac)

          IF (Rad>0.0d0) THEN
            vFall=(VL*VelSedi(ix,iy,iz-1,i)+VC*VelSedi(ix,iy,iz,i))/(VL+VC+Eps) ! Calculation of sedimentation velocity at surface

            Diff=DIFF_coef(Rad)
            Schmidt=SchmidtF(Diff)

            Brown=BrownF(Schmidt,LandUse)
            Stoke=StokeF(vFall,ustern,LandUse)
            Impact=ImpactF(Stoke,LandUse)
            Interception=InterceptionF(Diff,LandUse)
            IF (RH.GE.9.5d-1) THEN
              R1=1.d0 ! wet surface 
            ELSE
              R1=EXP(-SQRT(Stoke))
            END IF
            Rb=1.d0/(3d0*ustern*(Brown+Impact+Interception)*R1)

            Vdepos=(Ra+Rb+Ra*Rb*vFall)**(-1) ! [m/s] Addition von vFall erfolgt in AdvectionCompute

            ! Anzahl 
            fVec(1)%c(ix,iy,iz,i)=fVec(1)%c(ix,iy,iz,i)-cVec(1)%c(ix,iy,iz,i)*Vdepos*FL/(VolC(ix,iy,iz)+Eps)
            ! Massen 
            DO j=3,nAqua
              fVec(j)%c(ix,iy,iz,i)=fVec(j)%c(ix,iy,iz,i)-cVec(j)%c(ix,iy,iz,i)*Vdepos*FL/(VolC(ix,iy,iz)+Eps)
            END DO
          END IF
        END IF
      END DO
    END IF
    IF (ChemieGas) THEN
      DO i=1,NumSeaTrans

        IF (SeaGasEmi(i)%Pos>0) THEN
          IF (SeaGasEmi(i)%Schmidt.NE.0) THEN
            Rb=5*SeaGasEmi(i)%Schmidt**(2.d0/3.d0)/ustern

            Kw=TransCoeff(U10,SeaGasEmi(i)%Schmidt)

            Rc=1/(Kw*SeaGasEmi(i)%Henry)
            Vdepos=(Ra+Rb+Rc)**(-1)
            fVec(SeaGasEmi(i)%Pos)%c(ix,iy,iz,1)=fVec(SeaGasEmi(i)%Pos)%c(ix,iy,iz,1) &
                          -1*cVec(SeaGasEmi(i)%Pos)%c(ix,iy,iz,1)*Vdepos*FL/(VolC(ix,iy,iz)+Eps)
          ELSE
            Vdepos=SeaGasEmi(i)%Depos
            fVec(SeaGasEmi(i)%Pos)%c(ix,iy,iz,1)=fVec(SeaGasEmi(i)%Pos)%c(ix,iy,iz,1) &
                          -1*cVec(SeaGasEmi(i)%Pos)%c(ix,iy,iz,1)*Vdepos*FL/(VolC(ix,iy,iz)+Eps)
          END IF
        END IF

      END DO
    END IF
  END DO

END SUBROUTINE Deposition

SUBROUTINE DepositionJac(MatrixMatVectorLU,VelF,cVec,Jac) 

  TYPE(SpMatrix4Cell_T) :: MatrixMatVectorLU
  TYPE(VelocityFace_T) :: VelF(:)
  TYPE(Vec4_T), POINTER :: cVec(:)
  TYPE(Vec4_T), POINTER :: Jac(:)

  INTEGER, POINTER :: DiagP(:),Permu(:)
  REAL(RealKind), POINTER :: VelSedi(:,:,:,:), RhoLoc(:,:,:,:)
  REAL(RealKind), POINTER :: uF(:,:,:)
  REAL(RealKind), POINTER :: vF(:,:,:)
  REAL(RealKind), POINTER :: wF(:,:,:)

  INTEGER :: i,j,iCell,diag
  INTEGER :: ix,iy,iz,is
  INTEGER :: LandUse
  REAL(RealKind) :: U10,TSurface
  REAL(RealKind) :: n1,n2,n3,FL
  REAL(RealKind) :: uLoc,vLoc,wLoc,V,VT,VN
  REAL(RealKind) :: Schmidt,Brown,Stoke,Impact,Interception,Diff
  REAL(RealKind) :: Vdepos,Rb,Ra,R1,Rc,Kw
  REAL(RealKind) :: qv,RhoPart
  REAL(RealKind) :: DragH,DragM,ustern
  REAL(RealKind) :: nC,Rad
  REAL(RealKind) :: MM(0:nFrac+1,nAqua)
  REAL(RealKind) :: mFrac(nAqua)
  REAL(RealKind) :: VC,VL
  REAL(RealKind) :: vFall
  REAL(RealKind) :: RH

  RhoPart=1.2d3 ! kg/m³     ! Aus ASAM

  uF=>VelF(ibLoc)%uF
  vF=>VelF(ibLoc)%vF
  wF=>VelF(ibLoc)%wF
  VelSedi=>SediVelC%c
  RhoLoc=>cVec(RhoPos)%c
  DiagP=>MatrixMatVectorLU%Struct%DiagPtr(:)
  Permu=>MatrixMatVectorLU%Struct%Permu(:)

  DO iCell=1,NumBoundCell
    ix   = BoundCell(iCell)%ix
    iy   = BoundCell(iCell)%iy
    iz   = BoundCell(iCell)%iz
    n1   = BoundCell(iCell)%n1
    n2   = BoundCell(iCell)%n2
    n3   = BoundCell(iCell)%n3
    FL   = BoundCell(iCell)%FL
    DragH= BoundCell(iCell)%DragH
    DragM= BoundCell(iCell)%DragM
    TSurface = BoundCell(iCell)%TeS
    U10 = BoundCell(iCell)%U10
    LandUse = BoundCell(iCell)%LandClass
    qv=BoundCell(iCell)%qv

    RhoPart=1.2d3 ! kg/m³     ! Aus ASAM
    IF (LandUse.EQ.9) THEN
      RH=1.d0
    ELSE
      RH=qv*RV*TSurface/SaturVapor(TSurface)
    END IF

!   Evaluate velocity tangential and normal  
    uLoc = &
           (FU(ix-1,iy,iz)*uF(ix-1,iy,iz)+  &
            FU(ix,iy,iz)*uF(ix,iy,iz)) / &
          ((FU(ix-1,iy,iz)+FU(ix,iy,iz)+Eps)*(RhoLoc(ix,iy,iz,1)+Eps))
    vLoc = &
          (FV(ix,iy-1,iz)*vF(ix,iy-1,iz)+  &
           FV(ix,iy,iz)*vF(ix,iy,iz)) / &
         ((FV(ix,iy-1,iz)+FV(ix,iy,iz)+Eps)*(RhoLoc(ix,iy,iz,1)+Eps))
    wLoc = &
          (FW(ix,iy,iz-1)*wF(ix,iy,iz-1)+  &
           FW(ix,iy,iz)*wF(ix,iy,iz)) /  &
         ((FW(ix,iy,iz-1)+FW(ix,iy,iz)+Eps)*(RhoLoc(ix,iy,iz,1)+Eps))

    VN   = uLoc*n1+vLoc*n2+wLoc*n3
    V    = uLoc*uLoc+ &
           vLoc*vLoc+ &
           wLoc*wLoc
    VT   = SQRT(MAX(V-VN*VN,Zero)) ! Wind

    VC=VolC(ix,iy,iz)
    VL=VolC(ix,iy,iz-1)

    ustern=SQRT(DragM)*VT

    Ra=1.d0/(VT*DragH+Eps) ! Turbulent Resistenz ! Verglichen passt


    IF (Aerosol.AND.BoundCell(iCell)%LandClass.EQ.9) THEN
      DO i=1,nFrac

        nC=cVec(iNC)%c(ix,iy,iz,i)
        IF (nC.NE.0.0d0) THEN
          DO is=1,nAqua
            MM(i,is)=cVec(is)%c(ix,iy,iz,i)
          END DO
          mFrac=MM(i,:)/nC
          iF=iz
          Rad=Radius(mFrac)

          IF (Rad>0.0d0) THEN
            vFall=(VL*VelSedi(ix,iy,iz-1,i)+VC*VelSedi(ix,iy,iz,i))/(VL+VC+Eps)

            Diff=DIFF_coef(Rad)
            Schmidt=SchmidtF(Diff)

            Brown=BrownF(Schmidt,LandUse)
            Stoke=StokeF(vFall,ustern,LandUse)
            Impact=ImpactF(Stoke,LandUse)
            Interception=InterceptionF(Diff,LandUse)
            IF (RH.GE.9.d-1) THEN
              R1=1.d0 ! wet surface 
            ELSE
              R1=EXP(-SQRT(Stoke))
            END IF
            Rb=1.d0/(3.d0*ustern*(Brown+Impact+Interception)*R1+Eps)

            Vdepos=(Ra+Rb+Ra*Rb*vFall)**(-1.d0) ! [m/s]
            diag=DiagP(Permu(1))
            Jac(diag)%c(ix,iy,iz,i)=Jac(diag)%c(ix,iy,iz,i)-Vdepos*FL/(VolC(ix,iy,iz)+Eps)
            DO j=3,nAqua
              diag=DiagP(Permu(j))
              Jac(diag)%c(ix,iy,iz,i)=Jac(diag)%c(ix,iy,iz,i)-Vdepos*FL/(VolC(ix,iy,iz)+Eps)
            END DO
          END IF
        END IF
      END DO
    END IF
    IF (ChemieGas) THEN
      DO i=1,NumSeaTrans

        IF (SeaGasEmi(i)%Pos>0) THEN
          diag=DiagP(Permu(SeaGasEmi(i)%Pos))
          IF (SeaGasEmi(i)%Schmidt.NE.0.d0) THEN
            Rb=5*SeaGasEmi(i)%Schmidt**(2.d0/3.d0)/ustern

            Kw=TransCoeff(U10,SeaGasEmi(i)%Schmidt)

            Rc=1/(Kw*SeaGasEmi(i)%Henry)
            Vdepos=(Ra+Rb+Rc)**(-1)
            Jac(diag)%c(ix,iy,iz,1)=Jac(diag)%c(ix,iy,iz,1)-Vdepos*FL/(VolC(ix,iy,iz)+Eps)
          ELSE
            Vdepos=SeaGasEmi(i)%Depos
            Jac(diag)%c(ix,iy,iz,1)=Jac(diag)%c(ix,iy,iz,1)-Vdepos*FL/(VolC(ix,iy,iz)+Eps)
          END IF
        END IF

      END DO
    END IF
  END DO

END SUBROUTINE DepositionJac

SUBROUTINE SedimentVelocity(RhoZelle,SediZelle)

  TYPE(ScalarCell_T) :: RhoZelle, SediZelle
  INTEGER :: ix, iy, iz, i, is
  REAL(RealKind) :: nC,Rad
  REAL(RealKind) :: MM(0:nFrac+1,nAqua)
  REAL(RealKind) :: mFrac(nAqua)
  REAL(RealKind) :: RhoPart
  REAL(RealKind), POINTER :: VelSedi(:,:,:,:), Rhoair(:,:,:,:)

  Rhoair => RhoZelle%c
  VelSedi => SediZelle%c
  VelSedi = 0.E0
  RhoPart = 1.2d3 ! kg/m³     ! Aus ASAM

  DO ix=ix0,ix1+1
    DO iy=iy0,iy1+1
      DO iz=iz0,iz1+1
        DO i=1,nFrac
          nC=cVec(iNC)%c(ix,iy,iz,i)
          IF (nC.GT.0.0d0) THEN
            DO is=1,nAqua
              MM(i,is)=cVec(is)%c(ix,iy,iz,i)
            END DO
            mFrac=MM(i,:)/nC
            iCell=i
            iF=iz
            Rad=Radius(mFrac) ! m

            IF (Rad.GT.0.0d0) THEN
              VelSedi(ix,iy,iz,i) = VFinalF(RhoPart,Rhoair(ix,iy,iz,1),Rad)
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE SedimentVelocity

END MODULE Deposition_Mod
