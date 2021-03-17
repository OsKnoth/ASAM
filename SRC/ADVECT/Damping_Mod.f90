MODULE Damping_Mod
  USE DataType_Mod
  USE Domain_Mod
  USE Floor_Mod
  USE Physics_Mod
  USE Operator_Mod
  USE Parameter_mod
  USE Output_Mod
  USE Names_Mod
  USE Forcing_Mod

  IMPLICIT NONE


CONTAINS

SUBROUTINE SetPosE2Pos
  IF (uPosEnv>0) THEN
     PosE2Pos(uPosEnv)=uPosL
  END IF
  IF (vPosEnv>0) THEN
    PosE2Pos(vPosEnv)=vPosL
  END IF
  IF (wPosEnv>0) THEN
    PosE2Pos(wPosEnv)=wPosL
  END IF
  IF (thPosEnv>0) THEN
    PosE2Pos(thPosEnv)=thPos
  END IF
  IF (tkePosEnv>0) THEN
    PosE2Pos(tkePosEnv)=tkePos
  END IF
  IF (disPosEnv>0) THEN
    PosE2Pos(disPosEnv)=disPos
  END IF
  IF (RhoVPosEnv>0) THEN
    PosE2Pos(RhoVPosEnv)=RhoVPos
  END IF
  IF (RhoPosEnv>0) THEN
    PosE2Pos(RhoPosEnv)=RhoPos
  END IF
END SUBROUTINE SetPosE2Pos

SUBROUTINE ForceDamp(VectorCell,Rhs,Time,UVec)

  TYPE(Vector4Cell_T) :: VectorCell,Rhs
  REAL(RealKind) :: Time
  TYPE (Vector4Cell_T), OPTIONAL :: UVec

  INTEGER :: ic,ic0

  Rho=>VectorCell%Vec(RhoPos)%c
  SELECT CASE(DampType)
    CASE ('Interpol')
      DampKoeff=>DampKoeffCell(ibLoc)%c
      IF (PRESENT(UVec)) THEN
        DO ic=1,SIZE(VecEnv1(ibLoc)%Vec)
          ce1=>VecEnv1(ibLoc)%Vec(ic)%c
          ce2=>VecEnv2(ibLoc)%Vec(ic)%c
          IF (PosE2Pos(ic)==uPosL) THEN
            f=>Rhs%Vec(uPosL)%c
            c=>UVec%Vec(uPosL)%c
            CALL DampingInterpolCompute(Time)
            f=>Rhs%Vec(uPosR)%c
            c=>UVec%Vec(uPosR)%c
            CALL DampingInterpolCompute(Time)
          ELSE IF (PosE2Pos(ic)==vPosL) THEN
            f=>Rhs%Vec(vPosL)%c
            c=>UVec%Vec(vPosL)%c
            CALL DampingInterpolCompute(Time)
            f=>Rhs%Vec(vPosR)%c
            c=>UVec%Vec(vPosR)%c
            CALL DampingInterpolCompute(Time)
          ELSE IF (PosE2Pos(ic)==wPosL) THEN
            f=>Rhs%Vec(wPosL)%c
            c=>UVec%Vec(wPosL)%c
            CALL DampingInterpolCompute(Time)
            f=>Rhs%Vec(wPosR)%c
            c=>UVec%Vec(wPosR)%c
            CALL DampingInterpolCompute(Time)
          ELSE
            f=>Rhs%Vec(PosE2Pos(ic))%c
            c=>VectorCell%Vec(PosE2Pos(ic))%c
            CALL DampingInterpolCompute(Time)
          END IF
        END DO
      ELSE
        DO ic=1,SIZE(VecEnv1(ibLoc)%Vec)
          ce1=>VecEnv1(ibLoc)%Vec(ic)%c
          ce2=>VecEnv2(ibLoc)%Vec(ic)%c
          IF (PosE2Pos(ic)==uPosL) THEN
            f=>Rhs%Vec(uPosL)%c
            c=>VectorCell%Vec(uPosL)%c
            CALL DampingInterpolCompute(Time)
            f=>Rhs%Vec(uPosR)%c
            c=>VectorCell%Vec(uPosR)%c
            CALL DampingInterpolCompute(Time)
          ELSE IF (PosE2Pos(ic)==vPosL) THEN
            f=>Rhs%Vec(vPosL)%c
            c=>VectorCell%Vec(vPosL)%c
            CALL DampingInterpolCompute(Time)
            f=>Rhs%Vec(vPosR)%c
            c=>VectorCell%Vec(vPosR)%c
            CALL DampingInterpolCompute(Time)
          ELSE IF (PosE2Pos(ic)==wPosL) THEN
            f=>Rhs%Vec(wPosL)%c
            c=>VectorCell%Vec(wPosL)%c
            CALL DampingInterpolCompute(Time)
            f=>Rhs%Vec(wPosR)%c
            c=>VectorCell%Vec(wPosR)%c
            CALL DampingInterpolCompute(Time)
          ELSE
            f=>Rhs%Vec(PosE2Pos(ic))%c
            c=>VectorCell%Vec(PosE2Pos(ic))%c
            CALL DampingInterpolCompute(Time)
          END IF
        END DO
      END IF
    CASE ('Const')
      DampKoeff=>DampKoeffCell(ibLoc)%c
      IF (PRESENT(UVec)) THEN
        DO ic=1,SIZE(VecEnv1(ibLoc)%Vec)
          ce1=>VecEnv1(ibLoc)%Vec(ic)%c
          IF (PosE2Pos(ic)==uPosL) THEN
            f=>Rhs%Vec(uPosL)%c
            c=>UVec%Vec(uPosL)%c
            CALL DampingConstCompute(Time)
            f=>Rhs%Vec(uPosR)%c
            c=>UVec%Vec(uPosR)%c
            CALL DampingConstCompute(Time)
          ELSE IF (PosE2Pos(ic)==vPosL) THEN
            f=>Rhs%Vec(vPosL)%c
            c=>UVec%Vec(vPosL)%c
            CALL DampingConstCompute(Time)
            f=>Rhs%Vec(vPosR)%c
            c=>UVec%Vec(vPosR)%c
            CALL DampingConstCompute(Time)
          ELSE IF (PosE2Pos(ic)==wPosL) THEN
            f=>Rhs%Vec(wPosL)%c
            c=>UVec%Vec(wPosL)%c
            CALL DampingConstCompute(Time)
            f=>Rhs%Vec(wPosR)%c
            c=>UVec%Vec(wPosR)%c
            CALL DampingConstCompute(Time)
          ELSE
            f=>Rhs%Vec(PosE2Pos(ic))%c
            c=>VectorCell%Vec(PosE2Pos(ic))%c
            CALL DampingConstCompute(Time)
          END IF
        END DO
      ELSE
        DO ic=1,SIZE(VecEnv1(ibLoc)%Vec)
          ce1=>VecEnv1(ibLoc)%Vec(ic)%c
          IF (PosE2Pos(ic)==uPosL) THEN
            f=>Rhs%Vec(uPosL)%c
            c=>VectorCell%Vec(uPosL)%c
            CALL DampingConstCompute(Time)
            f=>Rhs%Vec(uPosR)%c
            c=>VectorCell%Vec(uPosR)%c
            CALL DampingConstCompute(Time)
          ELSE IF (PosE2Pos(ic)==vPosL) THEN
            f=>Rhs%Vec(vPosL)%c
            c=>VectorCell%Vec(vPosL)%c
            CALL DampingConstCompute(Time)
            f=>Rhs%Vec(vPosR)%c
            c=>VectorCell%Vec(vPosR)%c
            CALL DampingConstCompute(Time)
          ELSE IF (PosE2Pos(ic)==wPosL) THEN
            f=>Rhs%Vec(wPosL)%c
            c=>VectorCell%Vec(wPosL)%c
            CALL DampingConstCompute(Time)
            f=>Rhs%Vec(wPosR)%c
            c=>VectorCell%Vec(wPosR)%c
            CALL DampingConstCompute(Time)
          ELSE
            f=>Rhs%Vec(PosE2Pos(ic))%c
            c=>VectorCell%Vec(PosE2Pos(ic))%c
            CALL DampingConstCompute(Time)
          END IF
        END DO
      END IF
    CASE ('Average')
      DampKoeff=>DampKoeffCell(ibLoc)%c
      ic=1
      RhoAve=>RhoProfileMean(ibLoc)%Vec(1)%c
      DO ic0=1,VectorComponentsME
        IF (PosE2Pos(ic0)==uPosL) THEN
          f=>Rhs%Vec(uPosL)%c
          cAve=>ProfileMean(ibLoc)%Vec(ic)%c
          CALL DampingAverageCompute
          ic=ic+1
          f=>Rhs%Vec(uPosR)%c
          cAve=>ProfileMean(ibLoc)%Vec(ic)%c
          CALL DampingAverageCompute
          ic=ic+1
        ELSE IF (PosE2Pos(ic0)==vPosL) THEN
          f=>Rhs%Vec(vPosL)%c
          cAve=>ProfileMean(ibLoc)%Vec(ic)%c
          cProf=>NudgingProfile(ibLoc)%Vec(ic0)%c
          DampProf=>DampProfile(ibLoc)%Vec(ic0)%c
          CALL DampingAverageCompute
          ic=ic+1
          f=>Rhs%Vec(vPosR)%c
          cAve=>ProfileMean(ibLoc)%Vec(ic)%c
          CALL DampingAverageCompute
          ic=ic+1
        ELSE  
          f=>Rhs%Vec(PosE2Pos(ic0))%c
          cAve=>ProfileMean(ibLoc)%Vec(ic)%c
          cProf=>NudgingProfile(ibLoc)%Vec(ic0)%c
          DampProf=>DampProfile(ibLoc)%Vec(ic0)%c
          CALL DampingAverageCompute
          ic=ic+1
        END IF  
      END DO
  END SELECT
END SUBROUTINE ForceDamp

SUBROUTINE DampingInterpolCompute(Time)

  REAL(RealKind) :: Time
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: TempEnv

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
!       TempEnv=Rho(ix,iy,iz,1)*((Time2-Time)*ce1(ix,iy,iz,1) &
        TempEnv=((Time2-Time)*ce1(ix,iy,iz,1) &
               +(Time-Time1)*ce2(ix,iy,iz,1)) &
               /(Time2-Time1)
        f(ix,iy,iz,1)=f(ix,iy,iz,1)- &
               DampKoeff(ix,iy,iz,1)*(c(ix,iy,iz,1)-TempEnv)
      END DO
    END DO
  END DO

END SUBROUTINE DampingInterpolCompute

SUBROUTINE DampingConstCompute(Time)

  REAL(RealKind) :: Time
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: TempEnv

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        f(ix,iy,iz,1)=f(ix,iy,iz,1)- &
               DampKoeff(ix,iy,iz,1)*(c(ix,iy,iz,1)-ce1(ix,iy,iz,1))
      END DO
    END DO
  END DO

END SUBROUTINE DampingConstCompute

SUBROUTINE DampingAverageCompute

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: Temp

  DO iz=iz0+1,iz1
    Temp=DampProf(iz)*RelaxTime*(cAve(iz)-RhoAve(iz)*cProf(iz))
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        f(ix,iy,iz,1)=f(ix,iy,iz,1)-Temp
      END DO
    END DO
  END DO

END SUBROUTINE DampingAverageCompute

SUBROUTINE JacForceDamp(Jac)
  TYPE(SpMatrix4Cell_T) :: Jac

  INTEGER :: ic,icE

  SELECT CASE(DampType)
    CASE ('Interpol','Const')
      AS=>Jac%Mat
      DiagP=>Jac%Struct%DiagPtr(:)
      Permu=>Jac%Struct%Permu(:)
      DampKoeff=>DampKoeffCell(ibLoc)%c
      DO icE=1,SIZE(VecEnv1(ibLoc)%Vec)
        ic=PosE2Pos(icE)
        Diag=DiagP(Permu(ic))
        CALL JacDampingCompute 
        IF (ic==uPosL) THEN
          ic=uPosR
          Diag=DiagP(Permu(ic))
          CALL JacDampingCompute 
        ELSE IF (ic==vPosL) THEN
          ic=vPosR
          Diag=DiagP(Permu(ic))
          CALL JacDampingCompute 
        ELSE IF (ic==wPosL) THEN
          ic=wPosR
          Diag=DiagP(Permu(ic))
          CALL JacDampingCompute 
        END IF
      END DO
  END SELECT  
END SUBROUTINE JacForceDamp

SUBROUTINE JacDampingCompute

  INTEGER :: ix,iy,iz

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        AS(Diag)%c(ix,iy,iz,1)=AS(Diag)%c(ix,iy,iz,1)- &
               DampKoeff(ix,iy,iz,1)
      END DO
    END DO
  END DO

END SUBROUTINE JacDampingCompute

SUBROUTINE DampingInit

  INTEGER :: ix,iy,iz
  REAL(RealKind) :: xPLoc,yPLoc,zPLoc

  IF (Damping) THEN
    CALL Allocate(DampKoeffCell)
    DampKoeffCell=Zero
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO ix=ix0+1,ix1
        xPLoc=xP(ix-1)+Half*dx(ix)
!       West Boundary
        IF (xPLoc<=Domain%x0+StrideDamp%R_W) THEN
          DampKoeffCell(ibLoc)%c(ix,iy0+1:iy1,iz0+1:iz1,1) &
            =MAX(DampKoeffCell(ibLoc)%c(ix,iy0+1:iy1,iz0+1:iz1,1) &
                ,Damp%R_W*DampF(One-(xPLoc-Domain%x0)/StrideDamp%R_W))
        END IF
!       East Boundary
        IF (xPLoc>=Domain%x1-StrideDamp%R_E) THEN
          DampKoeffCell(ibLoc)%c(ix,iy0+1:iy1,iz0+1:iz1,1) &
            =MAX(DampKoeffCell(ibLoc)%c(ix,iy0+1:iy1,iz0+1:iz1,1) &
                ,Damp%R_E*DampF(One-(Domain%x1-xPLoc)/StrideDamp%R_E))
        END IF
      END DO
      DO iy=iy0+1,iy1
        yPLoc=yP(iy-1)+Half*dy(iy)
!       South Boundary
        IF (yPLoc<=Domain%y0+StrideDamp%R_S) THEN
          DampKoeffCell(ibLoc)%c(ix0+1:ix1,iy,iz0+1:iz1,1) &
            =MAX(DampKoeffCell(ibLoc)%c(ix0+1:ix1,iy,iz0+1:iz1,1) &
                ,Damp%R_S*DampF(One-(yPLoc-Domain%y0)/StrideDamp%R_S))
        END IF
!       North Boundary
        IF (yPLoc>=Domain%y1-StrideDamp%R_N) THEN
          DampKoeffCell(ibLoc)%c(ix0+1:ix1,iy,iz0+1:iz1,1) &
            =MAX(DampKoeffCell(ibLoc)%c(ix0+1:ix1,iy,iz0+1:iz1,1) &
                ,Damp%R_N*DampF(One-(Domain%y1-yPLoc)/StrideDamp%R_N))
        END IF
      END DO

      DO iz=iz0+1,iz1
        zPLoc=zP(iz-1)+Half*dz(iz)
!       Top Boundary
        IF (zPLoc>=Domain%z1-StrideDamp%R_T) THEN
          DampKoeffCell(ibLoc)%c(ix0+1:ix1,iy0+1:iy1,iz,1) &
            =MAX(DampKoeffCell(ibLoc)%c(ix0+1:ix1,iy0+1:iy1,iz,1) &
                ,Damp%R_T*DampF(One-(Domain%z1-zpLoc)/StrideDamp%R_T))
        END IF
      END DO
    END DO
  END IF
CONTAINS 
FUNCTION DampF(alpha)

  REAL(RealKind) :: DampF
  REAL(RealKind) :: alpha

  DampF=alpha
  DampF=SIN(Half*Pi*alpha)**2

END FUNCTION DampF

END SUBROUTINE DampingInit

END MODULE Damping_Mod
