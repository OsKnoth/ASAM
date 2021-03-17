MODULE Diagnostic_Mod
  USE DataType_Mod
  USE Names_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE Operator_Mod
  USE Chemie_Mod

  IMPLICIT NONE

CONTAINS

FUNCTION TotalScalar(VectorCell,Name)
  REAL(RealKind) :: TotalScalar
  TYPE(Vector4Cell_T) :: VectorCell(:)
  CHARACTER(*) :: Name

  INTEGER :: Pos
  REAL(RealKind) :: TotalScalarLoc

  Pos=Position(Name)
  TotalScalarLoc=0.0d0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    TotalScalarLoc=TotalScalarLoc+SUM(VectorCell(ibLoc)%Vec(Pos)%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1)*VolB)
  END DO
  CALL MPI_Allreduce(TotalScalarLoc,TotalScalar,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,MPIErr)

END FUNCTION TotalScalar

FUNCTION TotalValues(VectorCell,Velocityface)
  REAL(RealKind) :: TotalValues(1:6)
  TYPE(Vector4Cell_T) :: VectorCell(:)
  TYPE(VelocityFace_T) :: VelocityFace(:)

  INTEGER :: ix,iy,iz 
  REAL(RealKind) :: PotMLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoILoc,RhoDloc,pLoc,TLoc
  REAL(RealKind) :: Cvml,Lv
  REAL(RealKind) :: TotalEnergy,KineticEnergy,RhoTotal,PotentialEnergy
  REAL(RealKind) :: TotalEntropy
  REAL(RealKind) :: InternalEnergy


  TotalValues=Zero
  TotalEnergy=Zero
  KineticEnergy=Zero
  RhoTotal=Zero
  PotentialEnergy=Zero
  InternalEnergy=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
    IF (RhoCPos>0) THEN 
      RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
    ELSE
      RhoL=>RhoLCell(ibLoc)%c
    END IF
    RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
    RhoI=>VectorCell(ibLoc)%Vec(RhoIPos)%c
    RhoS=>VectorCell(ibLoc)%Vec(RhoSPos)%c
    Th=>VectorCell(ibLoc)%Vec(ThPos)%c
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          PotMLoc=Th(ix,iy,iz,1)
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoILoc=RhoI(ix,iy,iz,1)+RhoS(ix,iy,iz,1)
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc-RhoILoc
          pLoc=PressureTheta(RhoDLoc,RhoVLoc,RhoLLoc,RhoILoc,PotMLoc)
          TLoc=AbsTemp(RhoDLoc,RhoVLoc,pLoc)
          Cvml=CvmlFun(RhoDLoc,RhoVLoc,RhoLLoc,RhoILoc)
          Lv=LatHeat(TLoc)
          RhoVLoc=RhoV(ix,iy,iz,1)
          TotalEnergy=TotalEnergy &
          +Cvml*TLoc-Lv*RhoLLoc &
          +Half*(FaceToCell(uF(ix-1,iy,iz),uF(ix,iy,iz),FU(ix-1,iy,iz),FU(ix,iy,iz))**2 &
           +FaceToCell(vF(ix,iy-1,iz),vF(ix,iy,iz),FV(ix,iy-1,iz),FV(ix,iy,iz))**2 &
           +FaceToCell(wF(ix,iy,iz-1),wF(ix,iy,iz),FW(ix,iy,iz-1),FW(ix,iy,iz))**2 &
           )/(RhoLoc+Eps) &
          +Half*RhoLoc*(zP(iz-1)+zP(iz))*Grav
          InternalEnergy=InternalEnergy&
           +Cvml*TLoc-Lv*RhoLLoc
          KineticEnergy=KineticEnergy &
            +Half*(FaceToCell(uF(ix-1,iy,iz),uF(ix,iy,iz),FU(ix-1,iy,iz),FU(ix,iy,iz))**2 &
           +FaceToCell(vF(ix,iy-1,iz),vF(ix,iy,iz),FV(ix,iy-1,iz),FV(ix,iy,iz))**2 &
           +FaceToCell(wF(ix,iy,iz-1),wF(ix,iy,iz),FW(ix,iy,iz-1),FW(ix,iy,iz))**2 &
           )/(RhoLoc+Eps)
          PotentialEnergy=PotentialEnergy &
           +Half*RhoLoc*(zP(iz-1)+zP(iz))*Grav
          RhoTotal=RhoTotal &
           +RhoLoc   
          TotalEntropy=TotalEntropy+PotMLoc
        END DO
      END DO
    END DO
  END DO
  TotalValues(1)=TotalEnergy
  TotalValues(2)=KineticEnergy
  TotalValues(3)=InternalEnergy
  TotalValues(4)=PotentialEnergy
  TotalValues(5)=RhoTotal
  TotalValues(6)=TotalEntropy
END FUNCTION TotalValues

FUNCTION KineticEnergy(VectorCell)
  REAL(RealKind) :: KineticEnergy
  TYPE(Vector4Cell_T) :: VectorCell(:)

  KineticEnergy=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    DO ic=1,3
      c=>VectorCell(ibLoc)%Vec(ic)%c
      KineticEnergy=KineticEnergy &
         +SUM(c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1) &          
         *c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1) &          
         /(Rho(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1)+Eps))           
    END DO
  END DO
END FUNCTION KineticEnergy

FUNCTION Momentum(VelocityFace)
  REAL(RealKind) :: Momentum
  TYPE(VelocityFace_T) :: VelocityFace(:)

  Momentum=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    Momentum=Momentum+Half*SUM(VelocityFace(ibLoc)%uF(ix0:ix1-1,iy0+1:iy1,iz0+1:iz1) &
            *VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
    Momentum=Momentum+Half*SUM(VelocityFace(ibLoc)%uF(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1) &
            *VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
    Momentum=Momentum+Half*SUM(VelocityFace(ibLoc)%vF(ix0:ix1,iy0:iy1-1,iz0+1:iz1) &
            *VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
    Momentum=Momentum+Half*SUM(VelocityFace(ibLoc)%vF(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1) &
            *VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
    Momentum=Momentum+Half*SUM(VelocityFace(ibLoc)%wF(ix0:ix1,iy0:iy1,iz0:iz1-1) &
            *VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
    Momentum=Momentum+Half*SUM(VelocityFace(ibLoc)%wF(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1) &
            *VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
  END DO
END FUNCTION Momentum

FUNCTION TotalDensity(VectorCell)
  REAL(RealKind) :: TotalDensity
  TYPE(Vector4Cell_T) :: VectorCell(:)

  TotalDensity=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    TotalDensity=TotalDensity+SUM(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1) &
                *VectorCell(ibLoc)%Vec(RhoPos)%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1))
  END DO
END FUNCTION TotalDensity


FUNCTION TotalEntropy(VectorCell)
  REAL(RealKind) :: TotalEntropy
  TYPE(Vector4Cell_T) :: VectorCell(:)

  TotalEntropy=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    TotalEntropy=TotalEntropy+SUM(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1) &
                *VectorCell(ibLoc)%Vec(ThPos)%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1))
  END DO
END FUNCTION TotalEntropy

SUBROUTINE PrepareLWP(LWPCell,VectorCell)
  TYPE(ScalarCell_T), POINTER :: LWPCell(:)
  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
END SUBROUTINE PrepareLWP

SUBROUTINE PrepareEn(VectorCell,VelocityFace,Time)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE(VelocityFace_T), TARGET :: VelocityFace(:)
  REAL(RealKind) :: Time

  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
    RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
    IF (RhoCPos>0) THEN 
      RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
    ELSE
      RhoL=>RhoLCell(ibLoc)%c
    END IF  
    CALL Set(Floor(ib))
    KinEn=>KinEnCell(ibLoc)%c
    CALL KinEnCompute
    T=>TAbsCell(ibLoc)%Vec(1)%c
    CALL AbsTCompute
    E=>VectorCell(ibLoc)%Vec(enPos)%c
    CALL ECompute
  END DO
  CALL ExchangeCell(VectorCell)

END SUBROUTINE PrepareEn

SUBROUTINE PressureEnergyCompute(VectorCell,VelocityFace)

  TYPE (Vector4Cell_T), POINTER :: VectorCell(:)
  TYPE (VelocityFace_T), TARGET :: VelocityFace(:)
   
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDLoc
  REAL(RealKind) :: Kinetic,Energy,zLoc,pLoc

  INTEGER :: ix,iy,iz

  VelocityFaceAct=>VelocityFace
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    p=>PreCell(ibLoc)%c
    KinEn=>KinEnCell(ibLoc)%c
    Rho=>VectorCell(ibLoc)%Vec(RhoPos)%c
    E=>VectorCell(ibLoc)%Vec(enPos)%c
    RhoV=>VectorCell(ibLoc)%Vec(RhoVPos)%c
    RhoR=>VectorCell(ibLoc)%Vec(RhoRPos)%c
    IF (RhoCPos>0) THEN 
       RhoL=>VectorCell(ibLoc)%Vec(RhoCPos)%c
    ELSE
       RhoL=>RhoLCell(ibLoc)%c
    END IF  
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc
          Energy=E(ix,iy,iz,1)
          Kinetic=KinEn(ix,iy,iz,1)
          zLoc=Half*(zp(iz-1)+zp(iz))
          pLoc=(RhoDLoc*Rd+RhoVLoc*Rv) &
                        *(Energy &
                         -RhoLoc*Kinetic &
                         -RhoLoc*zLoc*Grav &
                         -RhoVLoc*L00) &
                        /(RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps)
          p(ix,iy,iz,1)=pLoc
        END DO 
      END DO 
    END DO 
  END DO 
  CALL ExchangeCell(PreCell)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    CALL SetVelocityFace(ibLoc)
    VectorCell(ibLoc)%Vec(thPos)%c=PreCell(ibLoc)%c
  END DO  


END SUBROUTINE PressureEnergyCompute

SUBROUTINE MeanProfileCompute(VecT) 

! Profil wird aus allen Zellen berechnet
! Voraussetzung: alle Blocks weisen die originale z-Einteilung auf

  TYPE(Vector4Cell_T), POINTER :: VecT(:)
  INTEGER :: ic,ix,iy,iz
  REAL(RealKind), ALLOCATABLE :: VolS(:)

  ALLOCATE(VolS(Domain%iz0+1:Domain%iz1))
  VolS=0.d0
  MeanProfile=0.d0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL DomainSet(ib)
    IF (iz0<Domain%iz0.OR.iz1>Domain%iz1) THEN
      STOP 'SR MeanProfileCompute qualified only for blocks with original z-grid'
    END IF
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          VolS(iz)=VolS(iz)+VolC(ix,iy,iz)
          MeanProfile(0,iz)=MeanProfile(0,iz)+VecT(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1)*VolC(ix,iy,iz)
          DO ic=1,UBOUND(VecT(ibLoc)%Vec,1)
            MeanProfile(ic,iz)=MeanProfile(ic,iz)+VecT(ibLoc)%Vec(ic)%c(ix,iy,iz,1)*VolC(ix,iy,iz)
          END DO
        END DO
      END DO
    END DO
  END DO
  DO iz=Domain%iz0+1,Domain%iz1
    MeanProfile(0,iz)=MeanProfile(0,iz)/(ABS(VolS(iz))+Eps)
    DO ic=1,UBOUND(MeanProfile,1)
      MeanProfile(ic,iz)=MeanProfile(ic,iz)/(ABS(VolS(iz))+Eps)
    END DO
  END DO

  DEALLOCATE(VolS)

END SUBROUTINE MeanProfileCompute

SUBROUTINE SoundCompute

  INTEGER :: ix,iy,iz

  REAL(RealKind) :: RhoLoc,RhoDLoc,RhoVLoc,RhoLLoc
  REAL(RealKind) :: qdLoc,qvLoc,qlLoc,Cvml
  REAL(RealKind) :: TLoc,eLoc,pLoc,PotEn,KinEnergy,SoSLoc,SoSLoc1,PreFac
  REAL(RealKind) :: DpDRho,DpDe,DpDRhoV,DpDRhoL

  SELECT CASE (ThetaKind)
    CASE ('PreEn')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            RhoLoc=Rho(ix,iy,iz,1)
            pLoc=Th(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc
            TLoc=pLoc/(RhoDLoc*Rd+RhoVLoc*Rv+Eps)
            Cvml=RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps
            eLoc=(Cvml*TLoc+RhoVLoc*L00)/(RhoLoc+Eps)
            DpDRho=Rd*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                   -Cvd*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2 &
                   +eLoc*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml
            DpDe=RhoLoc*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml
            DpDRhoV=(Rv-Rd)*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                    -L00*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml &
                    -(Cvv-Cvd)*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2 
            DpDRhoL=-Rd*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                    -(Cpl-Cvd)*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2
            SoSLoc=DpDRho+pLoc/(RhoLoc+Eps)**2*DpDe+RhoVLoc/(RhoLoc+Eps)*DpDRhoV+RhoLLoc/(RhoLoc+Eps)*DpDRhoL 
            PreFac=RhoLoc*SoSLoc
            Sound(ix,iy,iz,1)=PreFac
          END DO
        END DO
      END DO
    CASE('Density')  
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            RhoLoc=Rho(ix,iy,iz,1)
            pLoc=p(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc
            TLoc=T(ix,iy,iz,1)
            Cvml=RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps
            eLoc=(Cvml*TLoc+RhoVLoc*L00)/(RhoLoc+Eps)
            DpDRho=Rd*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                   -Cvd*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2 &
                   +eLoc*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml
            DpDe=RhoLoc*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml
            DpDRhoV=(Rv-Rd)*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                    -L00*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml &
                    -(Cvv-Cvd)*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2 
            DpDRhoL=-Rd*(RhoLoc*eLoc-RhoVLoc*L00)/Cvml &
                    -(Cpl-Cvd)*(RhoLoc*eLoc-RhoVLoc*L00)*(RhoDLoc*Rd+RhoVLoc*Rv)/Cvml**2
            SoSLoc=DpDRho+pLoc/(RhoLoc+Eps)**2*DpDe+RhoVLoc/(RhoLoc+Eps)*DpDRhoV+RhoLLoc/(RhoLoc+Eps)*DpDRhoL 
            PreFac=RhoLoc*SoSLoc
            Sound(ix,iy,iz,1)=PreFac
          END DO
        END DO
      END DO
    CASE ('Exner')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            Sound(ix,iy,iz,1)=Th(ix,iy,iz,1)*Rd/Cvd  
          END DO
        END DO
      END DO
  END SELECT

END SUBROUTINE SoundCompute

SUBROUTINE KinEnCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: uFL,uFR
  REAL(RealKind) :: vFL,vFR
  REAL(RealKind) :: wFL,wFR
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        uFL=uF(ix-1,iy,iz)
        uFR=uF(ix,iy,iz)
        vFL=vF(ix,iy-1,iz)
        vFR=vF(ix,iy,iz)
        wFL=wF(ix,iy,iz-1)
        wFR=wF(ix,iy,iz)
        RhoLoc=Rho(ix,iy,iz,1)
        KinEn(ix,iy,iz,1)= &  ! ohne Rho (KinEn=v^2/2)
                           0.25d0*(uFL*uFL+uFR*uFR+vFL*vFL+vFR*vFR+wFL*wFL+wFR*wFR) &
                           /(RhoLoc+Eps)/(RhoLoc+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE KinEnCompute 

SUBROUTINE KinEnCompute2

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: uFL,uFR
  REAL(RealKind) :: vFL,vFR
  REAL(RealKind) :: wFL,wFR
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        uFL=uF(ix-1,iy,iz)*(VolC(ix-1,iy,iz)+VolC(ix,iy,iz)) &
                          /(VolC(ix-1,iy,iz)*Rho(ix-1,iy,iz,1) &
                           +VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        uFR=uF(ix,iy,iz)*(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)) &
                          /(VolC(ix+1,iy,iz)*Rho(ix+1,iy,iz,1) &
                           +VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        vFL=vF(ix,iy-1,iz)*(VolC(ix,iy-1,iz)+VolC(ix,iy,iz)) &
                          /(VolC(ix,iy-1,iz)*Rho(ix,iy-1,iz,1) &
                           +VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        vFR=vF(ix,iy,iz)*(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)) &
                          /(VolC(ix,iy+1,iz)*Rho(ix,iy+1,iz,1) &
                           +VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        wFL=wF(ix,iy,iz-1)*(VolC(ix,iy,iz-1)+VolC(ix,iy,iz)) &
                          /(VolC(ix,iy,iz-1)*Rho(ix,iy,iz-1,1) &
                           +VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        wFR=wF(ix,iy,iz)*(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)) &
                          /(VolC(ix,iy,iz+1)*Rho(ix,iy,iz+1,1) &
                           +VolC(ix,iy,iz)*Rho(ix,iy,iz,1)+Eps)
        RhoLoc=Rho(ix,iy,iz,1)
        KinEn(ix,iy,iz,1)= &  ! ohne Rho (KinEn=v^2/2)
                           0.25d0*(uFL*uFL+uFR*uFR+vFL*vFL+vFR*vFR+wFL*wFL+wFR*wFR)
        KinEn(ix,iy,iz,1)= &  ! ohne Rho (KinEn=v^2/2)
                           0.125d0*(uFL+uFR)**2.0d0+(vFL+vFR)**2.0d0+(wFL+wFR)*2.0d0
      END DO
    END DO
  END DO

END SUBROUTINE KinEnCompute2 

SUBROUTINE KinEnCompute3

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc
  REAL(RealKind) :: uFL,uFR
  REAL(RealKind) :: vFL,vFR
  REAL(RealKind) :: wFL,wFR
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        uFL=uF(ix-1,iy,iz)
        uFR=uF(ix,iy,iz)
        vFL=vF(ix,iy-1,iz)
        vFR=vF(ix,iy,iz)
        wFL=wF(ix,iy,iz-1)
        wFR=wF(ix,iy,iz)
        RhoLoc=Rho(ix,iy,iz,1)
        KinEn(ix,iy,iz,1)=0.125d0*((uFL+uFR)**2.0d0+(vFL+vFR)**2.0d0+(wFL+wFR)**2.0d0)/(RhoLoc+Eps)**2.0d0
      END DO
    END DO
  END DO

END SUBROUTINE KinEnCompute3 

SUBROUTINE KinEnCompute4

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: uFLS,vFLS,wFLS
  REAL(RealKind) :: uFRS,vFRS,wFRS
  REAL(RealKind) :: Temp
  REAL(RealKind) :: KinV(ix0:ix1,iy0:iy1,iz0:iz1)
  REAL(RealKind) :: VolNode(ix0:ix1,iy0:iy1,iz0:iz1)

  KinV=0.0d0
  VoLNode=0.0d0
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=Two*uF(ix-1,iy,iz)/(Rho(ix-1,iy,iz,1)+Rho(ix,iy,iz,1)+Eps)
        uFLS=Temp*Temp
        Temp=Two*uF(ix,iy,iz)/(Rho(ix+1,iy,iz,1)+Rho(ix,iy,iz,1)+Eps)
        uFRS=Temp*Temp
        Temp=Two*vF(ix,iy-1,iz)/(Rho(ix,iy-1,iz,1)+Rho(ix,iy,iz,1)+Eps)
        vFLS=Temp*Temp
        Temp=Two*vF(ix,iy,iz)/(Rho(ix,iy+1,iz,1)+Rho(ix,iy,iz,1)+Eps)
        vFRS=Temp*Temp
        Temp=Two*wF(ix,iy,iz-1)/(Rho(ix,iy,iz-1,1)+Rho(ix,iy,iz,1)+Eps)
        wFLS=Temp*Temp
        Temp=Two*wF(ix,iy,iz)/(Rho(ix,iy,iz+1,1)+Rho(ix,iy,iz,1)+Eps)
        wFRS=Temp*Temp
        KinV(ix-1,iy-1,iz-1)=KinV(ix-1,iy-1,iz-1)+VolC(ix,iy,iz)*(uFLS+vFLS+wFLS) 
        KinV(ix  ,iy-1,iz-1)=KinV(ix  ,iy-1,iz-1)+VolC(ix,iy,iz)*(uFRS+vFLS+wFLS) 
        KinV(ix-1,iy  ,iz-1)=KinV(ix-1,iy  ,iz-1)+VolC(ix,iy,iz)*(uFLS+vFRS+wFLS) 
        KinV(ix-1,iy-1,iz  )=KinV(ix-1,iy-1,iz  )+VolC(ix,iy,iz)*(uFLS+vFLS+wFRS) 
        KinV(ix-1,iy  ,iz  )=KinV(ix-1,iy  ,iz  )+VolC(ix,iy,iz)*(uFLS+vFRS+wFRS) 
        KinV(ix  ,iy-1,iz  )=KinV(ix  ,iy-1,iz  )+VolC(ix,iy,iz)*(uFRS+vFLS+wFRS) 
        KinV(ix  ,iy  ,iz-1)=KinV(ix  ,iy  ,iz-1)+VolC(ix,iy,iz)*(uFRS+vFRS+wFLS) 
        KinV(ix  ,iy  ,iz  )=KinV(ix  ,iy  ,iz  )+VolC(ix,iy,iz)*(uFRS+vFRS+wFRS) 
        VolNode(ix-1,iy-1,iz-1)=VolNode(ix-1,iy-1,iz-1)+VolC(ix,iy,iz)
        VolNode(ix  ,iy-1,iz-1)=VolNode(ix  ,iy-1,iz-1)+VolC(ix,iy,iz)
        VolNode(ix-1,iy  ,iz-1)=VolNode(ix-1,iy  ,iz-1)+VolC(ix,iy,iz)
        VolNode(ix-1,iy-1,iz  )=VolNode(ix-1,iy-1,iz  )+VolC(ix,iy,iz)
        VolNode(ix-1,iy  ,iz  )=VolNode(ix-1,iy  ,iz  )+VolC(ix,iy,iz)
        VolNode(ix  ,iy-1,iz  )=VolNode(ix  ,iy-1,iz  )+VolC(ix,iy,iz)
        VolNode(ix  ,iy  ,iz-1)=VolNode(ix  ,iy  ,iz-1)+VolC(ix,iy,iz)
        VolNode(ix  ,iy  ,iz  )=VolNode(ix  ,iy  ,iz  )+VolC(ix,iy,iz)
      END DO
    END DO
  END DO
  KinV=KinV/(VolNode+Eps)
  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
         KinEn(ix,iy,iz,1)=1.0d0/8.0d0*( &
                           KinV(ix-1,iy-1,iz-1) &
                          +KinV(ix  ,iy-1,iz-1) &
                          +KinV(ix-1,iy  ,iz-1) &
                          +KinV(ix-1,iy-1,iz  ) &
                          +KinV(ix  ,iy  ,iz-1) &
                          +KinV(ix  ,iy-1,iz  ) &
                          +KinV(ix-1,iy  ,iz  ) &
                          +KinV(ix  ,iy  ,iz  )) 
      END DO
    END DO
  END DO

END SUBROUTINE KinEnCompute4 

SUBROUTINE AbsTPreCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Rm,Cpml,KappaLoc
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoILoc,RhoDLoc
  REAL(RealKind) :: TLoc


  SELECT CASE(ThetaKind)
    CASE('Density')
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoILoc=RhoI(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc-RhoILoc+Eps
            Rm=Rd*RhoDLoc+Rv*RhoVLoc
            Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc+Cpi*RhoILoc
            KappaLoc=Rm/Cpml
            T(ix,iy,iz,1)=(Rd*Th(ix,iy,iz,1)/p0**KappaLoc)**(One/(One-KappaLoc))/Rm
            p(ix,iy,iz,1)=Rm*T(ix,iy,iz,1)
          END DO
        END DO
      END DO
    CASE('Energy','EnergyBryan')
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
            TLoc=(Th(ix,iy,iz,1) &
                   -RhoLoc*KinEn(ix,iy,iz,1) &
                   -RhoLoc*Half*(zP(iz-1)+zP(iz))*Grav &
                   -RhoVLoc*L00) &
                  /(RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps)
            T(ix,iy,iz,1)=TLoc+Eps
            Rm=Rd*RhoDLoc+Rv*RhoVLoc
            p(ix,iy,iz,1)=Rm*T(ix,iy,iz,1)
          END DO
        END DO
      END DO      
    CASE('PreEn')
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
            Rm=Rd*RhoDLoc+Rv*RhoVLoc
            TLoc=(RhoEn(ix,iy,iz,1) &
                   -RhoLoc*KinEn(ix,iy,iz,1) &
                   -RhoLoc*Half*(zP(iz-1)+zP(iz))*Grav &
                   -RhoVLoc*L00) &
                  /(RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps)
            T(ix,iy,iz,1)=TLoc+Eps
            p(ix,iy,iz,1)=Th(ix,iy,iz,1)
          END DO
        END DO
      END DO      
  END SELECT 
END SUBROUTINE AbsTPreCompute

SUBROUTINE CpmlCompute
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDLoc,RhoILoc

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=Rho(ix,iy,iz,1)
        RhoVLoc=RhoV(ix,iy,iz,1)
        RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
        RhoILoc=RhoI(ix,iy,iz,1)
        RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc-RhoILoc+Eps
        Cpml(ix,iy,iz,1)=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc+Cpi*RhoILoc
      END DO
    END DO
  END DO
END SUBROUTINE CpmlCompute


SUBROUTINE AbsTCompute(TotalEn,RhoLocEn)

  INTEGER :: ix,iy,iz
  REAL(RealKind), OPTIONAL :: TotalEn,RhoLocEn
  REAL(RealKind) :: ThLoc,RhoLoc,RhoVLoc,RhoLLoc,RhoDLoc,RhoILoc,TLoc
  REAL(RealKind) :: TLoc1
  REAL(RealKind) :: Rm,Cpml,KappaLoc,Cp_eff
  REAL(RealKind) :: Acc,NewtonError,NewtonFunc,AbsT
  REAL(RealKind) :: pLoc

  SELECT CASE(ThetaKind)
    CASE('Density')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoILoc=RhoI(ix,iy,iz,1)!+RhoS(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc-RhoILoc+Eps
            Rm=Rd*RhoDLoc+Rv*RhoVLoc+Eps
            Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc+Cpi*RhoILoc+Eps
            KappaLoc=Rm/Cpml
            pLoc=(Rd*Th(ix,iy,iz,1)/p0**KappaLoc)**(One/(One-KappaLoc))
            T(ix,iy,iz,1)=(Rd*Th(ix,iy,iz,1)/p0**KappaLoc)**(One/(One-KappaLoc))/Rm+Eps
          END DO
        END DO
      END DO
    CASE('Equiv')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoILoc=RhoI(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc-RhoILoc+Eps
            ThLoc=Th(ix,iy,iz,1)/(RhoLoc+Eps)
            Rm=Rd*RhoDLoc+Rv*RhoVLoc
            Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc+Cpi*RhoILoc
            TLoc=T(ix,iy,iz,1)
            CALL AbsTNewton(TLoc,ThLoc,RhoDLoc,RhoVLoc,RhoLLoc)
            T(ix,iy,iz,1)=TLoc+Eps
          END DO
        END DO
      END DO
    CASE('Energy')
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
            TLoc=(Th(ix,iy,iz,1) &
                   -RhoLoc*KinEn(ix,iy,iz,1) &
                   -RhoLoc*Half*(zP(iz-1)+zP(iz))*Grav &
                   -RhoVLoc*L00) &
                  /(RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps)
            T(ix,iy,iz,1)=TLoc+Eps
          END DO
        END DO
      END DO
    CASE('PreEn') !p
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
            TLoc=Th(ix,iy,iz,1)/(RhoDLoc*Rd+RhoVLoc*Rv+Eps)
            T(ix,iy,iz,1)=TLoc+Eps
          END DO
        END DO
      END DO
    CASE('Exner') 
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
            TLoc=p0*Th(ix,iy,iz,1)**(1.0d0/kappa)/(RhoDLoc*Rd)
            T(ix,iy,iz,1)=TLoc+Eps
          END DO
        END DO
      END DO
  END SELECT

END SUBROUTINE AbsTCompute

SUBROUTINE ECompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: RhoLoc,TLoc,RhoDLoc,RhoVLoc,RhoLLoc,KinEnLoc
  REAL(RealKind) :: Cvml 

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        RhoLoc=Rho(ix,iy,iz,1)
        RhoVLoc=RhoV(ix,iy,iz,1)
        RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
        RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc
        Cvml=Cvd*RhoDLoc+Cvv*RhoVLoc+Cpl*RhoLLoc
        TLoc=T(ix,iy,iz,1)
        KinEnLoc=KinEn(ix,iy,iz,1)
        E(ix,iy,iz,1)=Cvml*TLoc & 
                      +RhoLoc*KinEnLoc &
                      +RhoVLoc*L00 &
                      +Half*RhoLoc*(zP(iz-1)+zP(iz))*Grav
      END DO
    END DO
  END DO

END SUBROUTINE ECompute

SUBROUTINE PreCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Arg,ThLoc,RhoLoc,RhoVLoc,RhoLLoc,KappaLoc
  REAL(RealKind) :: Kinetic,Energy,zLoc
  REAL(RealKind) :: RhoDLoc,TLoc

  IF (Shallow) THEN
    DO iz=iz0,iz1+1
      DO iy=iy0,iy1+1
        DO ix=ix0,ix1+1
          p(ix,iy,iz,1)=Half*Grav*Th(ix,iy,iz,1)*Th(ix,iy,iz,1)
        END DO
      END DO
    END DO
  ELSE IF (Liquid) THEN
    DO iz=iz0,iz1+1
      DO iy=iy0,iy1+1
        DO ix=ix0,ix1+1
          RhoLoc=Rho(ix,iy,iz,1)
          ThLoc=Th(ix,iy,iz,1)/(RhoLoc+Eps)
          arg=(EXP(-alpha_liq*(ThLoc-t0_liq))- &
              alpha_liq/(rho0_liq*kappa_liq*cpl)) &
              /(RhoLoc/Rho0_liq-alpha_liq/(rho0_liq*kappa_liq*cpl))
          IF (arg>Zero) THEN
            p(ix,iy,iz,1)=-1.d0/kappa_liq*LOG(arg)+p0_liq
          END IF
        END DO
      END DO
    END DO
  ELSE IF (LiquidTam) THEN
    DO iz=iz0,iz1+1
      DO iy=iy0,iy1+1
        DO ix=ix0,ix1+1
          p(ix,iy,iz,1)=(Th(ix,iy,iz,1)*(1.0d0-1.0d0/gamma_tam)*Cpl_tam)**(gamma_tam) &
                        *(1.0d0/(p0+gamma_tam*p0_tam))**(gamma_tam-1.0d0)
        END DO
      END DO
    END DO
  ELSE  
    SELECT CASE(ThetaKind)
    CASE ('Pseudo')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            p(ix,iy,iz,1)=Half*cS0*Cs0/Rho0*Rho(ix,iy,iz,1)*Rho(ix,iy,iz,1)
          END DO
        END DO
      END DO
    CASE('PreEn') ! p (thPos) and E (enPos)
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            p(ix,iy,iz,1)=Th(ix,iy,iz,1)
          END DO
        END DO
      END DO
    CASE('EnergyBryan') ! En (thPos) 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
            TLoc=(Th(ix,iy,iz,1) &
                   -RhoLoc*Half*(zP(iz-1)+zP(iz))*Grav &
                   -RhoVLoc*L00) &
                  /(RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps)
            p(ix,iy,iz,1)=(Rd*(Rho(ix,iy,iz,1)-RhoL(ix,iy,iz,1)-RhoR(ix,iy,iz,1)) &
                          +(Rv-Rd)*RhoV(ix,iy,iz,1))*TLoc
          END DO
        END DO
      END DO
    CASE('EnergyBryanSlow')
      !Kinetic Energy times Rm/Cvm for the slow part
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
            TLoc=(-RhoLoc*KinEn(ix,iy,iz,1)) &
                  /(RhoDLoc*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl+Eps)
            p(ix,iy,iz,1)=(Rd*(Rho(ix,iy,iz,1)-RhoL(ix,iy,iz,1)-RhoR(ix,iy,iz,1)) &
                          +(Rv-Rd)*RhoV(ix,iy,iz,1))*TLoc
          END DO
        END DO
      END DO
    CASE('Exner') ! p (thPos) and E (enPos)
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            p(ix,iy,iz,1)=p0*Th(ix,iy,iz,1)**(1.0d0/kappa)
          END DO
        END DO
      END DO
    CASE DEFAULT
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            p(ix,iy,iz,1)=(Rd*(Rho(ix,iy,iz,1)-RhoL(ix,iy,iz,1)-RhoR(ix,iy,iz,1)) &
                          +(Rv-Rd)*RhoV(ix,iy,iz,1))*T(ix,iy,iz,1)
          END DO
        END DO
      END DO
    END SELECT
  END IF

IF (TypeW=='ow') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          p(ix0,iy,iz,:)=p(ix0+1,iy,iz,:)
        END DO
      END DO
  END IF

  IF (TypeE=='oe') THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          p(ix1+1,iy,iz,:)=p(ix1,iy,iz,:)
        END DO
      END DO
  END IF

  IF (TypeS=='os') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          p(ix,iy0,iz,:)=p(ix,iy0+1,iz,:)
        END DO
      END DO
  END IF
  IF (TypeN=='on') THEN
      DO iz=iz0+1,iz1
        DO ix=ix0+1,ix1
          p(ix,iy1+1,iz,:)=p(ix,iy1,iz,:)
        END DO
      END DO
  END IF

  IF (TypeB=='ob') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          p(ix,iy,iz0,:)=p(ix,iy,iz0+1,:)
        END DO
      END DO
  END IF
  IF (TypeT=='ot') THEN
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          p(ix,iy,iz1+1,:)=p(ix,iy,iz1,:)
        END DO
      END DO
  END IF

END SUBROUTINE PreCompute

SUBROUTINE RhoLCompute 

  INTEGER :: ix,iy,iz,i

  DO iz=iz0,iz1+1
    DO iy=iy0,iy1+1
      DO ix=ix0,ix1+1
        RhoL(ix,iy,iz,1)=Sum(WaterLiq(ix,iy,iz,:))
      END DO
    END DO
  END DO

END SUBROUTINE RhoLCompute

SUBROUTINE PreFacCompute

  INTEGER :: ix,iy,iz

  DO iz=iz0,iz1+1
    DO iy=iy0,iy1+1
      DO ix=ix0,ix1+1
        p(ix,iy,iz,1)=Rd*(Rd*Th(ix,iy,iz,1)/p0)**(Kappa/(One-Kappa))
      END DO
    END DO
  END DO

END SUBROUTINE PreFacCompute

SUBROUTINE PreFacFCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: dpdthetaL,dpdthetaR

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0,ix1
        dpdthetaL=Rd/(One-Kappa)*(Rd*Th(ix,iy,iz,1)/p0)**(Kappa/(One-Kappa))
        dpdthetaR=Rd/(One-Kappa)*(Rd*Th(ix+1,iy,iz,1)/p0)**(Kappa/(One-Kappa))
        pFU(ix,iy,iz)=(VolC(ix,iy,iz)*dpdthetaL+VolC(ix+1,iy,iz)*dpdthetaR) &
                       /(VolC(ix,iy,iz)+VolC(ix+1,iy,iz)+Eps)
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0,iy1
      DO ix=ix0+1,ix1
        dpdthetaL=Rd/(One-Kappa)*(Rd*Th(ix,iy,iz,1)/p0)**(Kappa/(One-Kappa))
        dpdthetaR=Rd/(One-Kappa)*(Rd*Th(ix,iy+1,iz,1)/p0)**(Kappa/(One-Kappa))
        pFV(ix,iy,iz)=(VolC(ix,iy,iz)*dpdthetaL+VolC(ix,iy+1,iz)*dpdthetaR) &
                       /(VolC(ix,iy,iz)+VolC(ix,iy+1,iz)+Eps)
      END DO
    END DO
  END DO
  DO iz=iz0+1,iz1
    DO iy=iy0,iy1
      DO ix=ix0+1,ix1
        dpdthetaL=Rd/(One-Kappa)*(Rd*Th(ix,iy,iz,1)/p0)**(Kappa/(One-Kappa))
        dpdthetaR=Rd/(One-Kappa)*(Rd*Th(ix,iy,iz+1,1)/p0)**(Kappa/(One-Kappa))
        pFW(ix,iy,iz)=(VolC(ix,iy,iz)*dpdthetaL+VolC(ix,iy,iz+1)*dpdthetaR) &
                       /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps)
      END DO
    END DO
  END DO

END SUBROUTINE PreFacFCompute

SUBROUTINE DpDThetaCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Rm,Cpml,Cp_eff,KappaLoc
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDLoc,RhoILoc

  SELECT CASE (ThetaKind)
    CASE('Density')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
            RhoILoc=RhoI(ix,iy,iz,1)
            RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc-RhoILoc+Eps
            Rm=Rd+Rv*RhoV(ix,iy,iz,1)
            Cpml=Cpd+Cpv*RhoVLoc+Cpl*RhoLLoc+Cpi*RhoILoc
            KappaLoc=Rm/Cpml
            dpdtheta(ix,iy,iz,1)=Rd/(One-KappaLoc)*(Rd*Th(ix,iy,iz,1)/p0)**(KappaLoc/(One-Kappa))
          END DO
        END DO
      END DO
    CASE('Equiv')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            ! Fehler OSSI
            Cp_eff=Cpd+(RhoV(ix,iy,iz,1)+RhoL(ix,iy,iz,1))/Rho(ix,iy,iz,1)*Cpl
            Rm=Rd+Rv*RhoV(ix,iy,iz,1)/Rho(ix,iy,iz,1)
            Cpml=Cpd+Cpv*RhoV(ix,iy,iz,1)/Rho(ix,iy,iz,1) &
                +Cpl*(RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1))/Rho(ix,iy,iz,1) &
                +Cpi*RhoI(ix,iy,iz,1)/Rho(ix,iy,iz,1)
            dpdtheta(ix,iy,iz,1)=T(ix,iy,iz,1)/Th(ix,iy,iz,1)*((Cpml-Rm)/Cp_eff) &
                                  *(Rho(ix,iy,iz,1)*Rd+RhoV(ix,iy,iz,1)*Rv)
          END DO
        END DO
      END DO
    CASE('Energy')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            dpdtheta(ix,iy,iz,1)=(Rho(ix,iy,iz,1)*Rd+RhoV(ix,iy,iz,1)*Rv) &
                                  /(Rho(ix,iy,iz,1)*Cvd+RhoV(ix,iy,iz,1)*Cvv+(RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1))*Cpl)
          END DO
        END DO
      END DO
    CASE('Pressure')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            dpdtheta(ix,iy,iz,1)=One
          END DO
        END DO
      END DO
    CASE('PreEn')
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            dpdtheta(ix,iy,iz,1)=One
          END DO
        END DO
      END DO
    CASE DEFAULT
      DO iz=iz0,iz1+1
        DO iy=iy0,iy1+1
          DO ix=ix0,ix1+1
            dpdtheta(ix,iy,iz,1)=Rd/(One-Kappa)*(Rd*Th(ix,iy,iz,1)/p0)**(Kappa/(One-Kappa))
          END DO
        END DO
      END DO
  END SELECT
END SUBROUTINE DpDThetaCompute


END MODULE Diagnostic_Mod
