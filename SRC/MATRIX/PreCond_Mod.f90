MODULE PreCond_Mod

  USE Control_Mod
  USE JacAccGrav_Mod

  IMPLICIT NONE

  REAL(RealKind), POINTER :: DTU(:,:,:)
  REAL(RealKind), POINTER :: DTV(:,:,:)
  REAL(RealKind), POINTER :: DTW(:,:,:)
  REAL(RealKind), POINTER :: DUT(:,:,:)
  REAL(RealKind), POINTER :: DUR(:,:,:)
  REAL(RealKind), POINTER :: DTT(:,:,:)
  REAL(RealKind), POINTER :: DUU(:,:,:)
  REAL(RealKind), POINTER :: DUV(:,:,:)
  REAL(RealKind), POINTER :: DUW(:,:,:)

CONTAINS

SUBROUTINE SchurPreFull(b)

  TYPE(PressureVelocity), TARGET :: b(:)

  TYPE(PressureVelocity), POINTER :: BVec
  TYPE(Nachbar_T), POINTER :: Nachbar

  INTEGER :: jx,jy,jz
  REAL(RealKind) :: TolRel,TolAbs
  INTEGER :: MaxIter,nCG
  INTEGER :: in
  REAL(RealKind) :: Temp


  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    BVec => b(ibLoc)
    CALL Set(Floor(ib))
    DTU=>DTUG(ibLoc)%uF
    DTV=>DTUG(ibLoc)%vF
    DTW=>DTUG(ibLoc)%wF
    DUT=>DUTG(ibLoc)%p
    DUR=>DURG(ibLoc)%p
    DTT=>DTTG(ibLoc)%p
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    CALL SchurPreBlock
  END DO  
CONTAINS   

SUBROUTINE SchurPreBlock  

  REAL(RealKind) :: p(1:UBOUND(BVec%P,1),ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)
  INTEGER :: ix,iy,iz

    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
!  -- Westlicher Rand --
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar<Refine) THEN
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*   & 
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
              (SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF 
      END IF 
!  -- Oestlicher Rand --
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar<Refine) THEN
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*   & 
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
              (SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Suedlicher Rand --
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
              FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*   &  
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)/ &
              (SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Noerdlicher Rand --
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
              FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*   &  
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)/ &
              (SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Unterer Rand --
      IF (Nachbar%nType(2:2) == 'b') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
              FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*   &  
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)/ &
              (SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Oberer Rand --
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
              FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*   &  
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)/ &
              (SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))+Eps) 
            END DO
          END DO
        END IF
      END IF
    END DO

    p=BVec%p
    BVec%p=0.0e0
    IF (MultiMu.OR.MultiTriT.OR.MultiEx) THEN
      IF (TypeW/='ow') THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%p(1,ix0+1,iy,iz)=BVec%p(1,ix0+1,iy,iz)   &
                            +beta0*dtP*FU(ix0,iy,iz)*DTU(ix0,iy,iz)*DTT(ix0+1,iy,iz)*BVec%u_w(iy,iz)  &
                            /(0.5d0*dx(ix0+1)*MetrXY(iy)*MetrXZ(iz)*FU(ix0,iy,iz)+Eps)
          END DO                  
        END DO                  
      ELSE IF (TypeW=='ow'.AND.BCP%West==1) THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%p(1,ix0+1,iy,iz)=BVec%p(1,ix0+1,iy,iz)   &
                            +Two*beta0*dtP*FU(ix0,iy,iz)*DTU(ix0,iy,iz)*DTT(ix0+1,iy,iz)*BVec%u_w(iy,iz)  &
                            /(dx(ix0+1)*MetrXY(iy)*MetrXZ(iz)*FU(ix0,iy,iz)+Eps)
          END DO                  
        END DO                  
      END IF
      IF (TypeE/='oe') THEN
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%p(1,ix1,iy,iz)=BVec%p(1,ix1,iy,iz)   &
                          -beta0*dtP*FU(ix1,iy,iz)*DTU(ix1,iy,iz)*DTT(ix1,iy,iz)*BVec%u_e(iy,iz)  &
                          /(0.5d0*dx(ix1)*MetrXY(iy)*MetrXZ(iz)*FU(ix1,iy,iz)+Eps)
          END DO                  
        END DO                  
      ELSE IF (TypeE=='oe'.AND.BCP%East==1) THEN
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%p(1,ix1,iy,iz)=BVec%p(1,ix1,iy,iz)   &
                          -Two*beta0*dtP*FU(ix1,iy,iz)*DTU(ix1,iy,iz)*DTT(ix1,iy,iz)*BVec%u_e(iy,iz)  &
                          /(dx(ix1)*MetrXY(iy)*MetrXZ(iz)*FU(ix1,iy,iz)+Eps)
          END DO                  
        END DO                  
      END IF
      IF (TypeS/='os') THEN
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%p(1,ix,iy0+1,iz)=BVec%p(1,ix,iy0+1,iz)   &
                            +beta0*dtP*FV(ix,iy0,iz)*DTV(ix,iy0,iz)*DTT(ix,iy0+1,iz)*BVec%v_s(ix,iz) & 
                            /(0.5d0*dy(iy0+1)*MetrYX(ix)*MetrYZ(iz)*FV(ix,iy0,iz)+Eps) 
          END DO                  
        END DO                  
      ELSE IF (TypeS=='os'.AND.BCP%South==1) THEN
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%p(1,ix,iy0+1,iz)=BVec%p(1,ix,iy0+1,iz)   &
                            +Two*beta0*dtP*FV(ix,iy0,iz)*DTV(ix,iy0,iz)*DTT(ix,iy0+1,iz)*BVec%v_s(ix,iz) & 
                            /(dy(iy0+1)*MetrYX(ix)*MetrYZ(iz)*FV(ix,iy0,iz)+Eps) 
          END DO                  
        END DO                  
      END IF
      IF (TypeN/='on') THEN     
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%p(1,ix,iy1,iz)=BVec%p(1,ix,iy1,iz)   &
                          -beta0*dtP*FV(ix,iy1,iz)*DTV(ix,iy1,iz)*DTT(ix,iy1,iz)*BVec%v_n(ix,iz) &
                          /(0.5d0*dy(iy1)*MetrYX(ix)*MetrYZ(iz)*FV(ix,iy1,iz)+Eps)
          END DO                  
        END DO                  
      ELSE IF (TypeN=='on'.AND.BCP%North==1) THEN     
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%p(1,ix,iy1,iz)=BVec%p(1,ix,iy1,iz)   &
                          -Two*beta0*dtP*FV(ix,iy1,iz)*DTV(ix,iy1,iz)*DTT(ix,iy1,iz)*BVec%v_n(ix,iz) &
                          /(dy(iy1)*MetrYX(ix)*MetrYZ(iz)*FV(ix,iy1,iz)+Eps)
          END DO                  
        END DO                  
      END IF
      IF (TypeB/='ob') THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%p(1,ix,iy,iz0+1)=BVec%p(1,ix,iy,iz0+1)   &
                            +beta0*dtP*FW(ix,iy,iz0)*DTW(ix,iy,iz0)*DTT(ix,iy,iz0+1)*BVec%w_b(ix,iy) &
                            /(0.5d0*dz(iz0+1)*FW(ix,iy,iz0)+Eps)
          END DO                  
        END DO                  
      ELSE IF (TypeB=='ob'.AND.BCP%Bottom==1) THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%p(1,ix,iy,iz0+1)=BVec%p(1,ix,iy,iz0+1)   &
                            +Two*beta0*dtP*FW(ix,iy,iz0)*DTW(ix,iy,iz0)*DTT(ix,iy,iz0+1)*BVec%w_b(ix,iy) &
                            /(dz(iz0+1)*FW(ix,iy,iz0)+Eps)
          END DO                  
        END DO                  
      END IF
      IF (TypeT/='ot') THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%p(1,ix,iy,iz1)=BVec%p(1,ix,iy,iz1)   &
                          -beta0*dtP*FW(ix,iy,iz1)*DTW(ix,iy,iz1)*DTT(ix,iy,iz1)*BVec%w_t(ix,iy) &
                          /(0.5d0*dz(iz1)*FW(ix,iy,iz1)+Eps)
          END DO                  
        END DO                  
      ELSE IF (TypeT=='ot'.AND.BCP%Top==1) THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%p(1,ix,iy,iz1)=BVec%p(1,ix,iy,iz1)   &
                              -Two*beta0*dtP*FW(ix,iy,iz1)*DTW(ix,iy,iz1)*DTT(ix,iy,iz1)*BVec%w_t(ix,iy) &
                              /(dz(iz1)*FW(ix,iy,iz1)+Eps)
          END DO                
        END DO                
      END IF
    ELSE IF (MultiTriTB.OR.MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:)   &
                            +Two*beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:) &
                            /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)
        BVec%p(2,ix0+1,:,:) = BVec%p(2,ix0+1,:,:)   &
                            +Two*beta0*dtP*FU(ix0,:,:)*BVec%u_w(:,:) &
                            /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:)   &
                          -Two*beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:) &
                          /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
        BVec%p(2,ix1,:,:) = BVec%p(2,ix1,:,:)   &
                          -Two*beta0*dtP*FU(ix1,:,:)*BVec%u_e(:,:) &
                          /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:)   &
                            +Two*beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)
        BVec%p(2,:,iy0+1,:)=BVec%p(2,:,iy0+1,:)   &
                            +Two*beta0*dtP*FV(:,iy0,:)*BVec%v_s(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:)   &
                          -Two*beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:) &
                          /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)
        BVec%p(2,:,iy1,:)=BVec%p(2,:,iy1,:)   &
                          -Two*beta0*dtP*FV(:,iy1,:)*BVec%v_n(:,:) &
                          /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1)   &
                            +Two*beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)
        BVec%p(2,:,:,iz0+1)=BVec%p(2,:,:,iz0+1)   &
                            +beta0*dtP*FW(:,:,iz0)*BVec%w_b(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1)   &
                          -Two*beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:) &
                          /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)
        BVec%p(2,:,:,iz1)=BVec%p(2,:,:,iz1)   &
                          -Two*beta0*dtP*FW(:,:,iz1)*BVec%w_t(:,:) &
                          /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)
      END IF
    END IF

    MaxIter=CGMaxIterSch
    TolRel=CGTolSch
    TolAbs=1.d-12
    nCG=SIZE(BVec%p)
    MaActual=>LaplDirichlet
    xx=>b
    bb=>b
    CALL VCycleSet(ibLoc)
    CALL VCycleCG(MaxIter=MaxIter,TolRel=TolRel,TolAbs=TolAbs)
    
    IF (MultiMu.OR.MultiTriT.OR.MultiEx) THEN
      IF (TypeW/='ow') THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_w(iy,iz)=PreFac2*(BVec%u_w(iy,iz) &
                      -beta0*dtP*FU(ix0,iy,iz)*DUT(ix0+1,iy,iz)*DUU(ix0,iy,iz)*BVec%p(1,ix0+1,iy,iz)) &
                      /(0.5d0*dx(ix0+1)*MetrXY(iy)*MetrXZ(iz)*FU(ix0,iy,iz)+Eps) 
          END DO            
        END DO            
      ELSE IF (TypeW=='ow'.AND.BCP%West==1) THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_w(iy,iz)=PreFac3*(BVec%u_w(iy,iz) &
                      -beta0*dtP*FU(ix0,iy,iz)*DUT(ix0+1,iy,iz)*DUU(ix0,iy,iz)*BVec%p(1,ix0+1,iy,iz)) &
                      /(dx(ix0+1)*MetrXY(iy)*MetrXZ(iz)*FU(ix0,iy,iz)+Eps) 
          END DO            
        END DO            
      END IF
      IF (TypeE/='oe') THEN
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_e(iy,iz)=PreFac2*(BVec%u_e(iy,iz) &
                      +beta0*dtP*FU(ix1,iy,iz)*DUT(ix1,iy,iz)*DUU(ix1,iy,iz)*BVec%p(1,ix1,iy,iz)) &
                      /(0.5d0*dx(ix1)*MetrXY(iy)*MetrXZ(iz)*FU(ix1,iy,iz)+Eps) 
          END DO            
        END DO            
      ELSE IF (TypeE=='oe'.AND.BCP%East==1) THEN
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_e(iy,iz)=PreFac3*(BVec%u_e(iy,iz) &
                      +beta0*dtP*FU(ix1,iy,iz)*DUT(ix1,iy,iz)*DUU(ix1,iy,iz)*BVec%p(1,ix1,iy,iz)) &
                      /(dx(ix1)*MetrXY(iy)*MetrXZ(iz)*FU(ix1,iy,iz)+Eps) 
          END DO            
        END DO            
      END IF
      IF (TypeS/='os') THEN
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_s(ix,iz)=PreFac2*(BVec%v_s(ix,iz) &
                      -beta0*dtP*FV(ix,iy0,iz)*DUT(ix,iy0+1,iz)*DUV(ix,iy0,iz)*BVec%p(1,ix,iy0+1,iz)) &
                      /(0.5d0*dy(iy0+1)*MetrYX(ix)*MetrYZ(iz)*FV(ix,iy0,iz)+Eps) 
          END DO            
        END DO            
      ELSE IF (TypeS=='os'.AND.BCP%South==1) THEN
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_s(ix,iz)=PreFac3*(BVec%v_s(ix,iz) &
                      -beta0*dtP*FV(ix,iy0,iz)*DUT(ix,iy0+1,iz)*DUV(ix,iy0,iz)*BVec%p(1,ix,iy0+1,iz)) &
                      /(dy(iy0+1)*MetrYX(ix)*MetrYZ(iz)*FV(ix,iy0,iz)+Eps) 
          END DO            
        END DO            
      END IF
      IF (TypeN/='on') THEN     
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_n(ix,iz)=PreFac2*(BVec%v_n(ix,iz) &
                      +beta0*dtP*FV(ix,iy1,iz)*DUT(ix,iy1,iz)*DUV(ix,iy1,iz)*BVec%p(1,ix,iy1,iz)) &
                      /(0.5d0*dy(iy1)*MetrYX(ix)*MetrYZ(iz)*FV(ix,iy1,iz)+Eps) 
          END DO            
        END DO            
      ELSE IF (TypeN=='on'.AND.BCP%North==1) THEN     
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_n(ix,iz)=PreFac3*(BVec%v_n(ix,iz) &
                      +beta0*dtP*FV(ix,iy1,iz)*DUT(ix,iy1,iz)*DUV(ix,iy1,iz)*BVec%p(1,ix,iy1,iz)) &
                      /(dy(iy1)*MetrYX(ix)*MetrYZ(iz)*FV(ix,iy1,iz)+Eps) 
          END DO            
        END DO            
      END IF
      IF (TypeB/='ob') THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%w_b(ix,iy)=PreFac2*(BVec%w_b(ix,iy) &
                      -beta0*dtP*FW(ix,iy,iz0)*DUT(ix,iy,iz0+1)*DUW(ix,iy,iz0)*BVec%p(1,ix,iy,iz0+1)) &
                      /(0.5d0*dz(iz0+1)*FW(ix,iy,iz0)+Eps)
          END DO            
        END DO            
      ELSE IF (TypeB=='ob'.AND.BCP%Bottom==1) THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%w_b(ix,iy)=PreFac3*(BVec%w_b(ix,iy) &
                      -beta0*dtP*FW(ix,iy,iz0)*DUT(ix,iy,iz0+1)*DUW(ix,iy,iz0)*BVec%p(1,ix,iy,iz0+1)) &
                      /(dz(iz0+1)*FW(ix,iy,iz0)+Eps)
          END DO            
        END DO            
      END IF
      IF (TypeT/='ot') THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%w_t(ix,iy)=PreFac2*(BVec%w_t(ix,iy) &
                      +beta0*dtP*FW(ix,iy,iz1)*DUT(ix,iy,iz1)*DUW(ix,iy,iz1)*BVec%p(1,ix,iy,iz1)) &
                      /(0.5d0*dz(iz1)*FW(ix,iy,iz1)+Eps)
          END DO            
        END DO            
      ELSE IF (TypeT=='ot'.AND.BCP%Top==1) THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%w_t(ix,iy)=PreFac3*(BVec%w_t(ix,iy) &
                      +beta0*dtP*FW(ix,iy,iz1)*DUT(ix,iy,iz1)*DUW(ix,iy,iz1)*BVec%p(1,ix,iy,iz1)) &
                      /(dz(iz1)*FW(ix,iy,iz1)+Eps)
          END DO            
        END DO            
      END IF
    ELSE IF (MultiTriTB) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_w(iy,iz)=PreFac2*(BVec%u_w(iy,iz) &
                      -beta0*dtP*DUT(ix0+1,iy,iz)*DUU(ix0,iy,iz)*BVec%p(1,ix0+1,iy,iz)) &
                      /(dx(ix0+1)*MetrXY(iy)*MetrXZ(iz))
          END DO
        END DO  
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_e(iy,iz)=PreFac2*(BVec%u_e(iy,iz) &
                      +beta0*dtP*DUT(ix1,iy,iz)*DUU(ix1,iy,iz)*BVec%p(1,ix1,iy,iz)) &
                      /(dx(ix1)*MetrXY(iy)*MetrXZ(iz))
          END DO
        END DO  
      END IF
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_s(ix,iz)=PreFac2*(BVec%v_s(ix,iz) &
                      -beta0*dtP*DUT(ix,iy0+1,iz)*DUV(ix,iy0,iz)*BVec%p(1,ix,iy0+1,iz)) &
                      /(dy(iy0+1)*MetrYX(ix)*MetrYZ(iz))
          END DO
        END DO  
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_n(ix,iz)=PreFac2*(BVec%v_n(ix,iz) &
                      +beta0*dtP*DUT(ix,iy1,iz)*DUV(ix,iy1,iz)*BVec%p(1,ix,iy1,iz)) &
                      /(dy(iy1)*MetrYX(ix)*MetrYZ(iz))
          END DO
        END DO  
      END IF
!!!!!!! Rho terme !!!!!!
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%w_b(:,:)=PreFac2*(BVec%w_b(:,:) &
                      -beta0*dtP*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)) &
                      /dz(iz0+1)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%w_t(:,:)=PreFac2*(BVec%w_t(:,:) &
                      +beta0*dtP*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)) &
                      /dz(iz1)
      END IF
    ELSE IF (MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_w(iy,iz)=PreFac2*(BVec%u_w(iy,iz) &
                      -beta0*dtP*DUU(ix0,iy,iz)*DUT(ix0+1,iy,iz)*BVec%p(1,ix0+1,iy,iz)  &
                      -beta0*dtP*DUU(ix0,iy,iz)*DUR(ix0+1,iy,iz)*BVec%p(2,ix0+1,iy,iz)) &
                      /(dx(ix0+1)*MetrXY(iy)*MetrXZ(iz))
          END DO
        END DO  
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_e(iy,iz)=PreFac2*(BVec%u_e(iy,iz) &
                      +beta0*dtP*DUU(ix1,iy,iz)*DUT(ix1,iy,iz)*BVec%p(1,ix1,iy,iz)  &
                      +beta0*dtP*DUU(ix1,iy,iz)*DUR(ix1,iy,iz)*BVec%p(2,ix1,iy,iz)) &
                      /(dx(ix1)*MetrXY(iy)*MetrXZ(iz))  
          END DO
        END DO  
      END IF
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_s(ix,iz)=PreFac2*(BVec%v_s(ix,iz) &
                      -beta0*dtP*DUV(ix,iy0,iz)*DUT(ix,iy0+1,iz)*BVec%p(1,ix,iy0+1,iz)  &
                      -beta0*dtP*DUV(ix,iy0,iz)*DUR(ix,iy0+1,iz)*BVec%p(2,ix,iy0+1,iz)) &
                      /(dy(iy0+1)*MetrYX(ix)*MetrYZ(iz))
          END DO
        END DO  
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_n(ix,iz)=PreFac2*(BVec%v_n(ix,iz) &
                      +beta0*dtP*DUV(ix,iy1,iz)*DUT(ix,iy1,iz)*BVec%p(1,ix,iy1,iz)  &
                      +beta0*dtP*DUV(ix,iy1,iz)*DUR(ix,iy1,iz)*BVec%p(2,ix,iy1,iz)) &
                      /(dy(iy1)*MetrYX(ix)*MetrYZ(iz))
          END DO
        END DO  
      END IF
!!!!!!! Rho terme !!!!!!
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%w_b(:,:)=PreFac2*(BVec%w_b(:,:) &
                      -beta0*dtP*DUW(:,:,iz0)*DUT(:,:,iz0+1)*BVec%p(1,:,:,iz0+1)  &
                      -beta0*dtP*DUW(:,:,iz0)*DUR(:,:,iz0+1)*BVec%p(2,:,:,iz0+1)) &
                      /dz(iz0+1)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%w_t(:,:)=PreFac2*(BVec%w_t(:,:) &
                      +beta0*dtP*DUW(:,:,iz1)*DUT(:,:,iz1)*BVec%p(1,:,:,iz1)  &
                      +beta0*dtP*DUW(:,:,iz1)*DUR(:,:,iz1)*BVec%p(2,:,:,iz1)) &
                      /dz(iz1)
      END IF
    END IF

    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
!  -- Westlicher Rand --
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar<Refine) THEN
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*   &  
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
              (SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF 
      END IF 
!  -- Oestlicher Rand --
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar<Refine) THEN
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*   &  
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
              (SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Suedlicher Rand --
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
              FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*   &  
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)/ &
              (SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Noerdlicher Rand --
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
              FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*   &  
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)/ &
              (SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Unterer Rand --
      IF (Nachbar%nType(2:2) == 'b') THEN
        jz0 = iz0
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
              FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*   &  
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)/ &
              (SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Oberer Rand --
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
              FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*   &  
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)/ &
              (SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))+Eps) 
            END DO
          END DO
        END IF
      END IF
    END DO
    BVec%p=p
END SUBROUTINE SchurPreBlock  
END SUBROUTINE SchurPreFull

SUBROUTINE PressurePre(Pr,r)


  TYPE(PressureVelocity), TARGET :: Pr(:),r(:) 
  INTEGER :: MaxIter 
  REAL(RealKind) :: TolRel,TolAbs,SumP

  xx=>Pr
  bb=>r
  MaActual=>LaplNeumann
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    MaxIter=CGMaxIterPre
    TolRel=CGTolPre
    TolAbs=1.d-12
    Pr(ibLoc)%p=r(ibLoc)%p
    CALL VCycleSet(ibLoc)
    Shift=ShiftLoc
    CALL VCycleCG(MaxIter=MaxIter,TolRel=TolRel,TolAbs=TolAbs)
    Shift=0.0d0
  END DO

END SUBROUTINE PressurePre

SUBROUTINE b1pBTx2Full(b,x)

  TYPE(PressureVelocity), TARGET :: b(:),x(:) 

  TYPE(PressureVelocity), POINTER :: BVec, XVec

  INTEGER :: ix,iz
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    BVec => b(ibLoc)
    XVec => x(ibLoc)
    DTU=>DTUG(ibLoc)%uF
    DTV=>DTUG(ibLoc)%vF
    DTW=>DTUG(ibLoc)%wF
    DUT=>DUTG(ibLoc)%p
    DUR=>DURG(ibLoc)%p
    DTT=>DTTG(ibLoc)%p
    CALL Set(Floor(ib))
    Bvec%p=-XVec%p
    IF (MultiMu.OR.MultiTriT.OR.MultiEx) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:) &
                            +beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:) &
                          -beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:)
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:) &
                            +beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:) &
                          -beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1) &
                            +beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1) &
                          -beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:)
      END IF
    ELSE IF (MultiTriTB.OR.MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:) &
                            +beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:)
        BVec%p(2,ix0+1,:,:)=BVec%p(2,ix0+1,:,:) &
                            +beta0*dtP*FU(ix0,:,:)*BVec%u_w(:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:) &
                          -beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:)
        BVec%p(2,ix1,:,:)=BVec%p(2,ix1,:,:) &
                          -beta0*dtP*FU(ix1,:,:)*BVec%u_e(:,:)
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:) &
                            +beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:)
        BVec%p(2,:,iy0+1,:)=BVec%p(2,:,iy0+1,:) &
                            +beta0*dtP*FV(:,iy0,:)*BVec%v_s(:,:)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:) &
                          -beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:)
        BVec%p(2,:,iy1,:)=BVec%p(2,:,iy1,:) &
                          -beta0*dtP*FV(:,iy1,:)*BVec%v_n(:,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1) &
                            +beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:)
        BVec%p(2,:,:,iz0+1)=BVec%p(2,:,:,iz0+1) &
                            +beta0*dtP*FW(:,:,iz0)*BVec%w_b(:,:)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1) &
                          -beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:)
        BVec%p(2,:,:,iz1)=BVec%p(2,:,:,iz1) &
                          -beta0*dtP*FW(:,:,iz1)*BVec%w_t(:,:)
      END IF
    END IF 
  END DO

END SUBROUTINE b1pBTx2Full

SUBROUTINE b2pBx1Full(b,x)

  TYPE(PressureVelocity), TARGET :: b(:),x(:) 

  INTEGER :: ix,iy,iz
  TYPE(PressureVelocity), POINTER :: BVec, XVec
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    
    BVec => b(ibLoc)
    XVec => x(ibLoc)
    DTU=>DTUG(ibLoc)%uF
    DTV=>DTUG(ibLoc)%vF
    DTW=>DTUG(ibLoc)%wF
    DUT=>DUTG(ibLoc)%p
    DUR=>DURG(ibLoc)%p
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    IF (MultiMu.OR.MultiTriT.OR.MultiEx) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_w(iy,iz)=XVec%u_w(iy,iz)           &
                          +beta0*dtP*FU(ix0,iy,iz)*DUT(ix0+1,iy,iz)*DUU(ix0,iy,iz)*BVec%p(1,ix0+1,iy,iz)
          END DO            
        END DO            
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            BVec%u_e(iy,iz)=XVec%u_e(iy,iz)           &
                          -beta0*dtP*FU(ix1,iy,iz)*DUT(ix1,iy,iz)*DUU(ix1,iy,iz)*BVec%p(1,ix1,iy,iz)
          END DO            
        END DO            
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_s(ix,iz)=XVec%v_s(ix,iz)           &
                          +beta0*dtP*FV(ix,iy0,iz)*DUT(ix,iy0+1,iz)*DUV(ix,iy0,iz)*BVec%p(1,ix,iy0+1,iz)
          END DO            
        END DO            
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        DO iz=iz0+1,iz1
          DO ix=ix0+1,ix1
            BVec%v_n(ix,iz)=XVec%v_n(ix,iz)           &
                          -beta0*dtP*FV(ix,iy1,iz)*DUT(ix,iy1,iz)*DUV(ix,iy1,iz)*BVec%p(1,ix,iy1,iz)
          END DO            
        END DO            
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%w_b(ix,iy)=XVec%w_b(ix,iy)           &
                          +beta0*dtP*FW(ix,iy,iz0)*DUT(ix,iy,iz0+1)*DUW(ix,iy,iz0)*BVec%p(1,ix,iy,iz0+1)
          END DO            
        END DO            
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            BVec%w_t(ix,iy)=XVec%w_t(ix,iy)           &
                          -beta0*dtP*FW(ix,iy,iz1)*DUT(ix,iy,iz1)*DUW(ix,iy,iz1)*BVec%p(1,ix,iy,iz1)
          END DO            
        END DO            
      END IF
    ELSE IF (MultiTriTB) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%u_w(:,:)=XVec%u_w(:,:)           &
                      +beta0*dtP*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%u_e(:,:)=XVec%u_e(:,:)           &
                      -beta0*dtP*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%v_s(:,:)=XVec%v_s(:,:)           &
                      +beta0*dtP*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%v_n(:,:)=XVec%v_n(:,:)           &
                      -beta0*dtP*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%w_b(:,:)=XVec%w_b(:,:)           &
                      +beta0*dtP*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1) &
                      -Half*beta0*dtP*GravComp*dz(iz0+1)*BVec%p(2,:,:,iz0+1)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%w_t(:,:)=XVec%w_t(:,:)           &
                      -beta0*dtP*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1) &
                      -Half*beta0*dtP*GravComp*dz(iz1)*BVec%p(2,:,:,iz1)
      END IF
    ELSE IF (MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%u_w(:,:)=XVec%u_w(:,:)           &
                      +beta0*dtP*DUU(ix0,:,:)*DUT(ix0+1,:,:)*BVec%p(1,ix0+1,:,:) &
                      +beta0*dtP*DUU(ix0,:,:)*DUR(ix0+1,:,:)*BVec%p(2,ix0+1,:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%u_e(:,:)=XVec%u_e(:,:)           &
                      -beta0*dtP*DUU(ix1,:,:)*DUT(ix1,:,:)*BVec%p(1,ix1,:,:) &
                      -beta0*dtP*DUU(ix1,:,:)*DUR(ix1,:,:)*BVec%p(2,ix1,:,:)                      
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%v_s(:,:)=XVec%v_s(:,:)           &
                      +beta0*dtP*DUV(:,iy0,:)*DUT(:,iy0+1,:)*BVec%p(1,:,iy0+1,:) &
                      +beta0*dtP*DUV(:,iy0,:)*DUR(:,iy0+1,:)*BVec%p(2,:,iy0+1,:)
     END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%v_n(:,:)=XVec%v_n(:,:)           &
                      -beta0*dtP*DUV(:,iy1,:)*DUT(:,iy1,:)*BVec%p(1,:,iy1,:) &
                      -beta0*dtP*DUV(:,iy1,:)*DUR(:,iy1,:)*BVec%p(2,:,iy1,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%w_b(:,:)=XVec%w_b(:,:)           &
                      +beta0*dtP*DUW(:,:,iz0)*DUT(:,:,iz0+1)*BVec%p(1,:,:,iz0+1) &
                      +beta0*dtP*DUW(:,:,iz0)*DUR(:,:,iz0+1)*BVec%p(2,:,:,iz0+1) &
                      -Half*beta0*dtP*GravComp*dz(iz0+1)*BVec%p(2,:,:,iz0+1)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%w_t(:,:)=XVec%w_t(:,:)           &
                      -beta0*dtP*DUW(:,:,iz1)*DUT(:,:,iz1)*BVec%p(1,:,:,iz1) &
                      -beta0*dtP*DUW(:,:,iz1)*DUR(:,:,iz1)*BVec%p(2,:,:,iz1) &
                      -Half*beta0*dtP*GravComp*dz(iz1)*BVec%p(2,:,:,iz1)
      END IF
    END IF
  END DO  

END SUBROUTINE b2pBx1Full

SUBROUTINE SchurPreDual(b)

  TYPE(PressureVelocity), TARGET :: b(:)

  TYPE(PressureVelocity), POINTER :: BVec
  TYPE(Nachbar_T), POINTER :: Nachbar

  INTEGER :: jx,jy,jz
  REAL(RealKind) :: TolRel,TolAbs
  INTEGER :: MaxIter,nCG
  INTEGER :: in
  REAL(RealKind) :: Temp


  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    BVec => b(ibLoc)
    CALL Set(Floor(ib))
    DTU=>DTUG(ibLoc)%uF
    DTV=>DTUG(ibLoc)%vF
    DTW=>DTUG(ibLoc)%wF
    DUT=>DUTG(ibLoc)%p
    DUR=>DURG(ibLoc)%p
    DTT=>DTTG(ibLoc)%p
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    CALL SchurPreBlock
  END DO  
CONTAINS   

SUBROUTINE SchurPreBlock  

  REAL(RealKind) :: p(1:UBOUND(BVec%P,1),ix0+1:ix1,iy0+1:iy1,iz0+1:iz1)

    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
!  -- Westlicher Rand --
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar<Refine) THEN
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*   & 
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
              (SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF 
      END IF 
!  -- Oestlicher Rand --
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar<Refine) THEN
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*   & 
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
              (SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Suedlicher Rand --
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
              FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*   &  
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)/ &
              (SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Noerdlicher Rand --
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
              FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*   &  
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)/ &
              (SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Unterer Rand --
      IF (Nachbar%nType(2:2) == 'b') THEN
        jz0 = iz0
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
              FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*   &  
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)/ &
              (SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Oberer Rand --
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
              FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*   &  
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)/ &
              (SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))+Eps) 
            END DO
          END DO
        END IF
      END IF
    END DO

    p=BVec%p
    BVec%p=0.0e0
    IF (MultiMu.OR.MultiTriT.OR.MultiEx) THEN
      IF (TypeW/='ow') THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:)   &
                            +beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:) &
                            /(0.5d0*VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)
      ELSE IF (TypeW=='ow'.AND.BCP%West==1) THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:)   &
                            +Two*beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:) &
                            /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeE/='oe') THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:)   &
                          -beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:) &
                          /(0.5d0*VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
      ELSE IF (TypeE=='oe'.AND.BCP%East==1) THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:)   &
                          -Two*beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:) &
                          /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeS/='os') THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:)   &
                            +beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:) &
                            /(0.5d0*VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)
      ELSE IF (TypeS=='os'.AND.BCP%South==1) THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:)   &
                            +Two*beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeN/='on') THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:)   &
                          -beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:) &
                          /(0.5d0*VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)
      ELSE IF (TypeN=='on'.AND.BCP%North==1) THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:)   &
                          -Two*beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:) &
                          /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeB/='ob') THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1)   &
                            +beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:) &
                            /(0.5d0*VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)
      ELSE IF (TypeB=='ob'.AND.BCP%Bottom==1) THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1)   &
                            +Two*beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)
      END IF
      IF (TypeT/='ot') THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1)   &
                          -beta0*beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:) &
                          /(0.5d0*VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)
      ELSE IF (TypeT=='ot'.AND.BCP%Top==1) THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1)   &
                          -Two*beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:) &
                          /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)
      END IF
    ELSE IF (MultiTriTB.OR.MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:)   &
                            +Two*beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:) &
                            /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)
        BVec%p(2,ix0+1,:,:) = BVec%p(2,ix0+1,:,:)   &
                            +Two*beta0*dtP*FU(ix0,:,:)*BVec%u_w(:,:) &
                            /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:)   &
                          -Two*beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:) &
                          /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
        BVec%p(2,ix1,:,:) = BVec%p(2,ix1,:,:)   &
                          -Two*beta0*dtP*FU(ix1,:,:)*BVec%u_e(:,:) &
                          /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:)   &
                            +Two*beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)
        BVec%p(2,:,iy0+1,:)=BVec%p(2,:,iy0+1,:)   &
                            +Two*beta0*dtP*FV(:,iy0,:)*BVec%v_s(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:)   &
                          -Two*beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:) &
                          /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)
        BVec%p(2,:,iy1,:)=BVec%p(2,:,iy1,:)   &
                          -Two*beta0*dtP*FV(:,iy1,:)*BVec%v_n(:,:) &
                          /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1)   &
                            +Two*beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)
        BVec%p(2,:,:,iz0+1)=BVec%p(2,:,:,iz0+1)   &
                            +beta0*dtP*FW(:,:,iz0)*BVec%w_b(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1)   &
                          -Two*beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:) &
                          /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)
        BVec%p(2,:,:,iz1)=BVec%p(2,:,:,iz1)   &
                          -Two*beta0*dtP*FW(:,:,iz1)*BVec%w_t(:,:) &
                          /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)
      END IF
    END IF

    MaxIter=CGMaxIterSch
    TolRel=CGTolSch
    TolAbs=1.d-12
    nCG=SIZE(BVec%p)
    MaActual=>LaplDirichlet
    xx=>b
    bb=>b
    CALL VCycleSet(ibLoc)
    CALL VCycleCG(MaxIter=MaxIter,TolRel=TolRel,TolAbs=TolAbs)
    
    IF (MultiMu.OR.MultiTriT.OR.MultiEx) THEN
      IF (TypeW/='ow') THEN 
        BVec%u_w(:,:)=PreFac2*(BVec%u_w(:,:) &
                      -beta0*dtP*FU(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)) &
                      /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+VolC(ix0,iy0+1:iy1,iz0+1:iz1)+Eps)/Two 
      ELSE IF (TypeW=='ow'.AND.BCP%West==1) THEN 
        BVec%u_w(:,:)=PreFac3*(BVec%u_w(:,:) &
                      -beta0*dtP*FU(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)) &
                      /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps) 
      END IF
      IF (TypeE/='oe') THEN
        BVec%u_e(:,:)=PreFac2*(BVec%u_e(:,:) &
                      +beta0*dtP*FU(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)) &
                      /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+VolC(ix1+1,iy0+1:iy1,iz0+1:iz1)+Eps)/Two 
      ELSE IF (TypeE=='oe'.AND.BCP%East==1) THEN
        BVec%u_e(:,:)=PreFac3*(BVec%u_e(:,:) &
                      +beta0*dtP*FU(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)) &
                      /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)   
      END IF
      IF (TypeS/='os') THEN
        BVec%v_s(:,:)=PreFac2*(BVec%v_s(:,:) &
                      -beta0*dtP*FV(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)) &
                      /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+VolC(ix0+1:ix1,iy0,iz0+1:iz1)+Eps)/Two  
      ELSE IF (TypeS=='os'.AND.BCP%South==1) THEN
        BVec%v_s(:,:)=PreFac3*(BVec%v_s(:,:) &
                      -beta0*dtP*FV(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)) &
                      /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeN/='on') THEN     
        BVec%v_n(:,:)=PreFac2*(BVec%v_n(:,:) &
                      +beta0*dtP*FV(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)) &
                      /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+VolC(ix0+1:ix1,iy1+1,iz0+1:iz1)+Eps)/Two  
      ELSE IF (TypeN=='on'.AND.BCP%North==1) THEN     
        BVec%v_n(:,:)=PreFac3*(BVec%v_n(:,:) &
                      +beta0*dtP*FV(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)) &
                      /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeB/='ob') THEN     
        BVec%w_b(:,:)=PreFac2*(BVec%w_b(:,:) &
                      -beta0*dtP*FW(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)  
      ELSE IF (TypeB=='ob'.AND.BCP%Bottom==1) THEN     
        BVec%w_b(:,:)=PreFac3*(BVec%w_b(:,:) &
                      -beta0*dtP*FW(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)  
      END IF
      IF (TypeT/='ot') THEN     
        BVec%w_t(:,:)=PreFac2*(BVec%w_t(:,:) &
                      +beta0*dtP*FW(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)  
      ELSE IF (TypeT=='ot'.AND.BCP%Top==1) THEN     
        BVec%w_t(:,:)=PreFac3*(BVec%w_t(:,:) &
                      +beta0*dtP*FW(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)  
      END IF
    ELSE IF (MultiTriTB) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
!       BVec%u_w(:,:)=(0.5e0*BVec%u_w(:,:) &
        BVec%u_w(:,:)=0.25d0*(BVec%u_w(:,:) &
                      -beta0*dtP*FU(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)) &
                      /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
!       BVec%u_e(:,:)=(0.5e0*BVec%u_e(:,:) &
        BVec%u_e(:,:)=0.25d0*(BVec%u_e(:,:) &
                      +beta0*dtP*FU(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)) &
                      /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
!       BVec%v_s(:,:)=(0.5e0*BVec%v_s(:,:) &
        BVec%v_s(:,:)=0.25d0*(BVec%v_s(:,:) &
                      -beta0*dtP*FV(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)) &
                      /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
!       BVec%v_n(:,:)=(0.5e0*BVec%v_n(:,:) &
        BVec%v_n(:,:)=0.25d0*(BVec%v_n(:,:) &
                      +beta0*dtP*FV(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)) &
                      /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)  
      END IF
!!!!!!! Rho terme !!!!!!
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
!       BVec%w_b(:,:)=(0.5e0*BVec%w_b(:,:) &
        BVec%w_b(:,:)=0.25d0*(BVec%w_b(:,:) &
                      -beta0*dtP*FW(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)  
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
!       BVec%w_t(:,:)=(0.5e0*BVec%w_t(:,:) &
        BVec%w_t(:,:)=0.25d0*(BVec%w_t(:,:) &
                      +beta0*dtP*FW(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)  
      END IF
    ELSE IF (MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
!       BVec%u_w(:,:)=(0.5e0*BVec%u_w(:,:) &
        BVec%u_w(:,:)=PreFac2*(BVec%u_w(:,:) &
                      -beta0*dtP*DUU(ix0,:,:)*FU(ix0,:,:)*DUT(ix0+1,:,:)*BVec%p(1,ix0+1,:,:)  &
                      -beta0*dtP*DUU(ix0,:,:)*FU(ix0,:,:)*DUR(ix0+1,:,:)*BVec%p(2,ix0+1,:,:)) &
                      /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps) 
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
!       BVec%u_e(:,:)=(0.5e0*BVec%u_e(:,:) &
        BVec%u_e(:,:)=PreFac2*(BVec%u_e(:,:) &
                      +beta0*dtP*DUU(ix1,:,:)*FU(ix1,:,:)*DUT(ix1,:,:)*BVec%p(1,ix1,:,:)  &
                      +beta0*dtP*DUU(ix1,:,:)*FU(ix1,:,:)*DUR(ix1,:,:)*BVec%p(2,ix1,:,:)) &
                      /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
!       BVec%v_s(:,:)=(0.5e0*BVec%v_s(:,:) &
        BVec%v_s(:,:)=PreFac2*(BVec%v_s(:,:) &
                      -beta0*dtP*DUV(:,iy0,:)*FV(:,iy0,:)*DUT(:,iy0+1,:)*BVec%p(1,:,iy0+1,:)  &
                      -beta0*dtP*DUV(:,iy0,:)*FV(:,iy0,:)*DUR(:,iy0+1,:)*BVec%p(2,:,iy0+1,:)) &
                      /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
!       BVec%v_n(:,:)=(0.5e0*BVec%v_n(:,:) &
        BVec%v_n(:,:)=PreFac2*(BVec%v_n(:,:) &
                      +beta0*dtP*DUV(:,iy1,:)*FV(:,iy1,:)*DUT(:,iy1,:)*BVec%p(1,:,iy1,:)  &
                      +beta0*dtP*DUV(:,iy1,:)*FV(:,iy1,:)*DUR(:,iy1,:)*BVec%p(2,:,iy1,:)) &
                      /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)  
      END IF
!!!!!!! Rho terme !!!!!!
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
!       BVec%w_b(:,:)=(0.5e0*BVec%w_b(:,:) &
        BVec%w_b(:,:)=PreFac2*(BVec%w_b(:,:) &
                      -beta0*dtP*DUW(:,:,iz0)*FW(:,:,iz0)*DUT(:,:,iz0+1)*BVec%p(1,:,:,iz0+1)  &
                      -beta0*dtP*DUW(:,:,iz0)*FW(:,:,iz0)*DUR(:,:,iz0+1)*BVec%p(2,:,:,iz0+1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)  
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
!       BVec%w_t(:,:)=(0.5e0*BVec%w_t(:,:) &
        BVec%w_t(:,:)=PreFac2*(BVec%w_t(:,:) &
                      +beta0*dtP*DUW(:,:,iz1)*FW(:,:,iz1)*DUT(:,:,iz1)*BVec%p(1,:,:,iz1)  &
                      +beta0*dtP*DUW(:,:,iz1)*FW(:,:,iz1)*DUR(:,:,iz1)*BVec%p(2,:,:,iz1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)  
      END IF
    END IF

    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
!  -- Westlicher Rand --
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar<Refine) THEN
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*   &  
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
              (SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF 
      END IF 
!  -- Oestlicher Rand --
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar<Refine) THEN
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*   &  
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)/ &
              (SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Suedlicher Rand --
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
              FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*   &  
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)/ &
              (SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Noerdlicher Rand --
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
              FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*   &  
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)/ &
              (SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Unterer Rand --
      IF (Nachbar%nType(2:2) == 'b') THEN
        jz0 = iz0
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
              FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*   &  
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)/ &
              (SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))+Eps) 
            END DO
          END DO
        END IF
      END IF
!  -- Oberer Rand --
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar<Refine) THEN
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
              FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*   &  
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)/ &
              (SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))+Eps) 
            END DO
          END DO
        END IF
      END IF
    END DO
    BVec%p=p
END SUBROUTINE SchurPreBlock  
END SUBROUTINE SchurPreDual

SUBROUTINE b1pBTx2Dual(b,x)

  TYPE(PressureVelocity), TARGET :: b(:),x(:) 

  TYPE(PressureVelocity), POINTER :: BVec, XVec

  INTEGER :: ix,iz
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    BVec => b(ibLoc)
    XVec => x(ibLoc)
    DTU=>DTUG(ibLoc)%uF
    DTV=>DTUG(ibLoc)%vF
    DTW=>DTUG(ibLoc)%wF
    DUT=>DUTG(ibLoc)%p
    DUR=>DURG(ibLoc)%p
    DTT=>DTTG(ibLoc)%p
    CALL Set(Floor(ib))
    Bvec%p=-XVec%p      !!!!!!!!!!!!!!!!!!!!
    IF (MultiMu.OR.MultiTriT.OR.MultiEx) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:) &
                            +beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:) &
                          -beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:)
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:) &
                            +beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:) &
                          -beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1) &
                            +beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1) &
                          -beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:)
      END IF
    ELSE IF (MultiTriTB.OR.MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:) &
                            +beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:)
        BVec%p(2,ix0+1,:,:)=BVec%p(2,ix0+1,:,:) &
                            +beta0*dtP*FU(ix0,:,:)*BVec%u_w(:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:) &
                          -beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:)
        BVec%p(2,ix1,:,:)=BVec%p(2,ix1,:,:) &
                          -beta0*dtP*FU(ix1,:,:)*BVec%u_e(:,:)
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:) &
                            +beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:)
        BVec%p(2,:,iy0+1,:)=BVec%p(2,:,iy0+1,:) &
                            +beta0*dtP*FV(:,iy0,:)*BVec%v_s(:,:)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:) &
                          -beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:)
        BVec%p(2,:,iy1,:)=BVec%p(2,:,iy1,:) &
                          -beta0*dtP*FV(:,iy1,:)*BVec%v_n(:,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1) &
                            +beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:)
        BVec%p(2,:,:,iz0+1)=BVec%p(2,:,:,iz0+1) &
                            +beta0*dtP*FW(:,:,iz0)*BVec%w_b(:,:)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1) &
                          -beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:)
        BVec%p(2,:,:,iz1)=BVec%p(2,:,:,iz1) &
                          -beta0*dtP*FW(:,:,iz1)*BVec%w_t(:,:)
      END IF
    END IF 
  END DO

END SUBROUTINE b1pBTx2Dual

SUBROUTINE b2pBx1Dual(b,x)

  TYPE(PressureVelocity), TARGET :: b(:),x(:) 

  TYPE(PressureVelocity), POINTER :: BVec, XVec
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    
    BVec => b(ibLoc)
    XVec => x(ibLoc)
    DTU=>DTUG(ibLoc)%uF
    DTV=>DTUG(ibLoc)%vF
    DTW=>DTUG(ibLoc)%wF
    DUT=>DUTG(ibLoc)%p
    DUR=>DURG(ibLoc)%p
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    IF (MultiMu.OR.MultiTriT.OR.MultiEx) THEN
      IF (TypeW/='ow') THEN 
        BVec%u_w(:,:)=PreFac5*XVec%u_w(:,:)           &
                      +beta0*dtP*FU(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)
      ELSE IF (TypeW=='ow'.AND.BCP%West==1) THEN 
        BVec%u_w(:,:)=XVec%u_w(:,:)           &
                      +beta0*dtP*FU(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)
      END IF
      IF (TypeE/='oe') THEN
        BVec%u_e(:,:)=PreFac5*XVec%u_e(:,:)           &
                      -beta0*dtP*FU(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)
      ELSE IF (TypeE=='oe'.AND.BCP%East==1) THEN
        BVec%u_e(:,:)=XVec%u_e(:,:)           &
                      -beta0*dtP*FU(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)
      END IF    
      IF (TypeS/='os') THEN
        BVec%v_s(:,:)=PreFac5*XVec%v_s(:,:)           &
                      +beta0*dtP*FV(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)
      ELSE IF (TypeS=='os'.AND.BCP%South==1) THEN
        BVec%v_s(:,:)=XVec%v_s(:,:)           &
                      +beta0*dtP*FV(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)
      END IF
      IF (TypeN/='on') THEN     
        BVec%v_n(:,:)=PreFac5*XVec%v_n(:,:)           &
                      -beta0*dtP*FV(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)
      ELSE IF (TypeN=='on'.AND.BCP%North==1) THEN     
        BVec%v_n(:,:)=XVec%v_n(:,:)           &
                      -beta0*dtP*FV(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)
      END IF
      IF (TypeB/='ob') THEN     
        BVec%w_b(:,:)=PreFac5*XVec%w_b(:,:)           &
                      +beta0*dtP*FW(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)
      ELSE IF (TypeB=='ob'.AND.BCP%Bottom==1) THEN     
        BVec%w_b(:,:)=XVec%w_b(:,:)           &
                      +beta0*dtP*FW(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)
      END IF
      IF (TypeT/='ot') THEN     
        BVec%w_t(:,:)=PreFac5*XVec%w_t(:,:)           &
                      -beta0*dtP*FW(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)
      ELSE IF (TypeT=='ot'.AND.BCP%Top==1) THEN     
        BVec%w_t(:,:)=XVec%w_t(:,:)           &
                      -beta0*dtP*FW(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)
      END IF
    ELSE IF (MultiTriTB) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%u_w(:,:)=XVec%u_w(:,:)           &
                      +beta0*dtP*FU(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%u_e(:,:)=XVec%u_e(:,:)           &
                      -beta0*dtP*FU(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%v_s(:,:)=XVec%v_s(:,:)           &
                      +beta0*dtP*FV(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%v_n(:,:)=XVec%v_n(:,:)           &
                      -beta0*dtP*FV(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%w_b(:,:)=XVec%w_b(:,:)           &
                      +beta0*dtP*FW(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1) &
                      -Half*beta0*dtP*GravComp*VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)*BVec%p(2,:,:,iz0+1)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%w_t(:,:)=XVec%w_t(:,:)           &
                      -beta0*dtP*FW(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1) &
                      -Half*beta0*dtP*GravComp*VolC(ix0+1:ix1,iy0+1:iy1,iz1)*BVec%p(2,:,:,iz1)
      END IF
    ELSE IF (MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%u_w(:,:)=XVec%u_w(:,:)           &
                      +beta0*dtP*DUU(ix0,:,:)*FU(ix0,:,:)*DUT(ix0+1,:,:)*BVec%p(1,ix0+1,:,:) &
                      +beta0*dtP*DUU(ix0,:,:)*FU(ix0,:,:)*DUR(ix0+1,:,:)*BVec%p(2,ix0+1,:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%u_e(:,:)=XVec%u_e(:,:)           &
                      -beta0*dtP*DUU(ix1,:,:)*FU(ix1,:,:)*DUT(ix1,:,:)*BVec%p(1,ix1,:,:) &
                      -beta0*dtP*DUU(ix1,:,:)*FU(ix1,:,:)*DUR(ix1,:,:)*BVec%p(2,ix1,:,:)                      
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%v_s(:,:)=XVec%v_s(:,:)           &
                      +beta0*dtP*DUV(:,iy0,:)*FV(:,iy0,:)*DUT(:,iy0+1,:)*BVec%p(1,:,iy0+1,:) &
                      +beta0*dtP*DUV(:,iy0,:)*FV(:,iy0,:)*DUR(:,iy0+1,:)*BVec%p(2,:,iy0+1,:)
     END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%v_n(:,:)=XVec%v_n(:,:)           &
                      -beta0*dtP*DUV(:,iy1,:)*FV(:,iy1,:)*DUT(:,iy1,:)*BVec%p(1,:,iy1,:) &
                      -beta0*dtP*DUV(:,iy1,:)*FV(:,iy1,:)*DUR(:,iy1,:)*BVec%p(2,:,iy1,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%w_b(:,:)=XVec%w_b(:,:)           &
                      +beta0*dtP*DUW(:,:,iz0)*FW(:,:,iz0)*DUT(:,:,iz0+1)*BVec%p(1,:,:,iz0+1) &
                      +beta0*dtP*DUW(:,:,iz0)*FW(:,:,iz0)*DUR(:,:,iz0+1)*BVec%p(2,:,:,iz0+1) &
                      -Half*beta0*dtP*GravComp*VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)*BVec%p(2,:,:,iz0+1)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%w_t(:,:)=XVec%w_t(:,:)           &
                      -beta0*dtP*DUW(:,:,iz1)*FW(:,:,iz1)*DUT(:,:,iz1)*BVec%p(1,:,:,iz1) &
                      -beta0*dtP*DUW(:,:,iz1)*FW(:,:,iz1)*DUR(:,:,iz1)*BVec%p(2,:,:,iz1) &
                      -Half*beta0*dtP*GravComp*VolC(ix0+1:ix1,iy0+1:iy1,iz1)*BVec%p(2,:,:,iz1)
      END IF
    END IF
  END DO   ! ibLoc

END SUBROUTINE b2pBx1Dual

SUBROUTINE SchurPre(b)

  TYPE(PressureVelocity), TARGET :: b(:)

  IF (GradFull) THEN
    CALL SchurPreFull(b)
  ELSE  
    CALL SchurPreDual(b)
  END IF
END SUBROUTINE SchurPre

SUBROUTINE b1pBTx2(b,x)

  TYPE(PressureVelocity), TARGET :: b(:),x(:) 
  IF (GradFull) THEN
    CALL b1pBTx2Full(b,x)
  ELSE  
    CALL b1pBTx2Dual(b,x)
  END IF
END SUBROUTINE b1pBTx2

SUBROUTINE b2pBx1(b,x)
  TYPE(PressureVelocity), TARGET :: b(:),x(:) 

  IF (GradFull) THEN
    CALL b2pBx1Full(b,x)
  ELSE  
    CALL b2pBx1Dual(b,x)
  END IF
END SUBROUTINE b2pBx1

SUBROUTINE PSolve1(b,x)

  TYPE(PressureVelocity) :: b(:),x(:) 
  REAL(RealKind) :: Temp

  b=x
  CALL PressurePre(b,x)
  CALL b2pBx1(b,x)
! -- Randwert-Austausch 
  CALL Exchange(b)
  CALL SchurPre(b)
!   -- Randwert-Austausch 
  CALL Exchange(b)
  CALL b1pBTx2(b,x) 
  CALL PressurePre(b,b)

END SUBROUTINE PSolve1

END MODULE PreCond_Mod
