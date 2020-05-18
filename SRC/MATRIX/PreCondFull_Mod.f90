MODULE PreCond_Mod

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

  REAL(RealKind) :: Fac2=0.25d0 !0.5d0
  REAL(RealKind) :: Fac3=0.5d0 !2.0d0
  REAL(RealKind) :: Fac4=2.0d0
  REAL(RealKind) :: Fac5=0.5d0
CONTAINS

SUBROUTINE SchurPre(b)

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
    IF (MultiMU.OR.MultiTriT) THEN
      IF (TypeW/='ow') THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:)   &
                            +Fac4*beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:) &
                            /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)
      ELSE IF (TypeW=='ow'.AND.BCP%West==1) THEN 
        BVec%p(1,ix0+1,:,:)=BVec%p(1,ix0+1,:,:)   &
                            +Two*beta0*dtP*FU(ix0,:,:)*DTU(ix0,:,:)*DTT(ix0+1,:,:)*BVec%u_w(:,:) &
                            /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeE/='oe') THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:)   &
                          -Fac4*beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:) &
                          /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
      ELSE IF (TypeE=='oe'.AND.BCP%East==1) THEN
        BVec%p(1,ix1,:,:)=BVec%p(1,ix1,:,:)   &
                          -Two*beta0*dtP*FU(ix1,:,:)*DTU(ix1,:,:)*DTT(ix1,:,:)*BVec%u_e(:,:) &
                          /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeS/='os') THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:)   &
                            +Fac4*beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)
      ELSE IF (TypeS=='os'.AND.BCP%South==1) THEN
        BVec%p(1,:,iy0+1,:)=BVec%p(1,:,iy0+1,:)   &
                            +Two*beta0*dtP*FV(:,iy0,:)*DTV(:,iy0,:)*DTT(:,iy0+1,:)*BVec%v_s(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeN/='on') THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:)   &
                          -Fac4*beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:) &
                          /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)
      ELSE IF (TypeN=='on'.AND.BCP%North==1) THEN     
        BVec%p(1,:,iy1,:)=BVec%p(1,:,iy1,:)   &
                          -Two*beta0*dtP*FV(:,iy1,:)*DTV(:,iy1,:)*DTT(:,iy1,:)*BVec%v_n(:,:) &
                          /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)
      END IF
      IF (TypeB/='ob') THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1)   &
                            +Two*beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)
      ELSE IF (TypeB=='ob'.AND.BCP%Bottom==1) THEN     
        BVec%p(1,:,:,iz0+1)=BVec%p(1,:,:,iz0+1)   &
                            +Two*beta0*dtP*FW(:,:,iz0)*DTW(:,:,iz0)*DTT(:,:,iz0+1)*BVec%w_b(:,:) &
                            /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)
      END IF
      IF (TypeT/='ot') THEN     
        BVec%p(1,:,:,iz1)=BVec%p(1,:,:,iz1)   &
                          -Two*beta0*dtP*FW(:,:,iz1)*DTW(:,:,iz1)*DTT(:,:,iz1)*BVec%w_t(:,:) &
                          /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)
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
    
    IF (MultiMU.OR.MultiTriT) THEN
      IF (TypeW/='ow') THEN 
        BVec%u_w(:,:)=Fac2*(BVec%u_w(:,:) &
                      -beta0*dtP*FUG(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)) &
                      /dx(ix0+1) 
!                     /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+VolC(ix0,iy0+1:iy1,iz0+1:iz1)+Eps)/Two 
      ELSE IF (TypeW=='ow'.AND.BCP%West==1) THEN 
        BVec%u_w(:,:)=Fac3*(BVec%u_w(:,:) &
                      -beta0*dtP*FUG(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)) &
                      /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps) 
      END IF
      IF (TypeE/='oe') THEN
        BVec%u_e(:,:)=Fac2*(BVec%u_e(:,:) &
                      +beta0*dtP*FUG(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)) &
                      /dx(ix1) 
!                     /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+VolC(ix1+1,iy0+1:iy1,iz0+1:iz1)+Eps)/Two 
      ELSE IF (TypeE=='oe'.AND.BCP%East==1) THEN
        BVec%u_e(:,:)=Fac3*(BVec%u_e(:,:) &
                      +beta0*dtP*FUG(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)) &
                      /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)   
      END IF
      IF (TypeS/='os') THEN
        BVec%v_s(:,:)=Fac2*(BVec%v_s(:,:) &
                      -beta0*dtP*FVG(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)) &
                      /dy(iy0+1)
!                     /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+VolC(ix0+1:ix1,iy0,iz0+1:iz1)+Eps)/Two  
      ELSE IF (TypeS=='os'.AND.BCP%South==1) THEN
        BVec%v_s(:,:)=Fac3*(BVec%v_s(:,:) &
                      -beta0*dtP*FVG(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)) &
                      /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeN/='on') THEN     
        BVec%v_n(:,:)=Fac2*(BVec%v_n(:,:) &
                      +beta0*dtP*FVG(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)) &
                      /dy(iy1+1)
!                     /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+VolC(ix0+1:ix1,iy1+1,iz0+1:iz1)+Eps)/Two  
      ELSE IF (TypeN=='on'.AND.BCP%North==1) THEN     
        BVec%v_n(:,:)=Fac3*(BVec%v_n(:,:) &
                      +beta0*dtP*FVG(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)) &
                      /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeB/='ob') THEN     
        BVec%w_b(:,:)=Fac2*(BVec%w_b(:,:) &
                      -beta0*dtP*FWG(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)  
      ELSE IF (TypeB=='ob'.AND.BCP%Bottom==1) THEN     
        BVec%w_b(:,:)=Fac3*(BVec%w_b(:,:) &
                      -beta0*dtP*FWG(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)  
      END IF
      IF (TypeT/='ot') THEN     
        BVec%w_t(:,:)=Fac2*(BVec%w_t(:,:) &
                      +beta0*dtP*FWG(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)  
      ELSE IF (TypeT=='ot'.AND.BCP%Top==1) THEN     
        BVec%w_t(:,:)=Fac3*(BVec%w_t(:,:) &
                      +beta0*dtP*FWG(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)  
      END IF
    ELSE IF (MultiTriTB) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
!       BVec%u_w(:,:)=(0.5e0*BVec%u_w(:,:) &
        BVec%u_w(:,:)=0.25d0*(BVec%u_w(:,:) &
                      -beta0*dtP*FUG(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)) &
                      /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
!       BVec%u_e(:,:)=(0.5e0*BVec%u_e(:,:) &
        BVec%u_e(:,:)=0.25d0*(BVec%u_e(:,:) &
                      +beta0*dtP*FUG(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)) &
                      /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
!       BVec%v_s(:,:)=(0.5e0*BVec%v_s(:,:) &
        BVec%v_s(:,:)=0.25d0*(BVec%v_s(:,:) &
                      -beta0*dtP*FVG(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)) &
                      /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
!       BVec%v_n(:,:)=(0.5e0*BVec%v_n(:,:) &
        BVec%v_n(:,:)=0.25d0*(BVec%v_n(:,:) &
                      +beta0*dtP*FVG(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)) &
                      /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)  
      END IF
!!!!!!! Rho terme !!!!!!
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
!       BVec%w_b(:,:)=(0.5e0*BVec%w_b(:,:) &
        BVec%w_b(:,:)=0.25d0*(BVec%w_b(:,:) &
                      -beta0*dtP*FWG(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)  
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
!       BVec%w_t(:,:)=(0.5e0*BVec%w_t(:,:) &
        BVec%w_t(:,:)=0.25d0*(BVec%w_t(:,:) &
                      +beta0*dtP*FWG(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz1)+Eps)  
      END IF
    ELSE IF (MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
!       BVec%u_w(:,:)=(0.5e0*BVec%u_w(:,:) &
        BVec%u_w(:,:)=Fac2*(BVec%u_w(:,:) &
                      -beta0*dtP*DUU(ix0,:,:)*FUG(ix0,:,:)*DUT(ix0+1,:,:)*BVec%p(1,ix0+1,:,:)  &
                      -beta0*dtP*DUU(ix0,:,:)*FUG(ix0,:,:)*DUR(ix0+1,:,:)*BVec%p(2,ix0+1,:,:)) &
                      /(VolC(ix0+1,iy0+1:iy1,iz0+1:iz1)+Eps) 
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
!       BVec%u_e(:,:)=(0.5e0*BVec%u_e(:,:) &
        BVec%u_e(:,:)=Fac2*(BVec%u_e(:,:) &
                      +beta0*dtP*DUU(ix1,:,:)*FUG(ix1,:,:)*DUT(ix1,:,:)*BVec%p(1,ix1,:,:)  &
                      +beta0*dtP*DUU(ix1,:,:)*FUG(ix1,:,:)*DUR(ix1,:,:)*BVec%p(2,ix1,:,:)) &
                      /(VolC(ix1,iy0+1:iy1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
!       BVec%v_s(:,:)=(0.5e0*BVec%v_s(:,:) &
        BVec%v_s(:,:)=Fac2*(BVec%v_s(:,:) &
                      -beta0*dtP*DUV(:,iy0,:)*FVG(:,iy0,:)*DUT(:,iy0+1,:)*BVec%p(1,:,iy0+1,:)  &
                      -beta0*dtP*DUV(:,iy0,:)*FVG(:,iy0,:)*DUR(:,iy0+1,:)*BVec%p(2,:,iy0+1,:)) &
                      /(VolC(ix0+1:ix1,iy0+1,iz0+1:iz1)+Eps)  
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
!       BVec%v_n(:,:)=(0.5e0*BVec%v_n(:,:) &
        BVec%v_n(:,:)=Fac2*(BVec%v_n(:,:) &
                      +beta0*dtP*DUV(:,iy1,:)*FVG(:,iy1,:)*DUT(:,iy1,:)*BVec%p(1,:,iy1,:)  &
                      +beta0*dtP*DUV(:,iy1,:)*FVG(:,iy1,:)*DUR(:,iy1,:)*BVec%p(2,:,iy1,:)) &
                      /(VolC(ix0+1:ix1,iy1,iz0+1:iz1)+Eps)  
      END IF
!!!!!!! Rho terme !!!!!!
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
!       BVec%w_b(:,:)=(0.5e0*BVec%w_b(:,:) &
        BVec%w_b(:,:)=Fac2*(BVec%w_b(:,:) &
                      -beta0*dtP*DUW(:,:,iz0)*FWG(:,:,iz0)*DUT(:,:,iz0+1)*BVec%p(1,:,:,iz0+1)  &
                      -beta0*dtP*DUW(:,:,iz0)*FWG(:,:,iz0)*DUR(:,:,iz0+1)*BVec%p(2,:,:,iz0+1)) &
                      /(VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)+Eps)  
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
!       BVec%w_t(:,:)=(0.5e0*BVec%w_t(:,:) &
        BVec%w_t(:,:)=Fac2*(BVec%w_t(:,:) &
                      +beta0*dtP*DUW(:,:,iz1)*FWG(:,:,iz1)*DUT(:,:,iz1)*BVec%p(1,:,:,iz1)  &
                      +beta0*dtP*DUW(:,:,iz1)*FWG(:,:,iz1)*DUR(:,:,iz1)*BVec%p(2,:,:,iz1)) &
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
END SUBROUTINE SchurPre

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

SUBROUTINE Project(b)

  TYPE(PressureVelocity) :: b(:)
  
  REAL(RealKind) :: y(nb)
  INTEGER :: i

! Projektion
  IF (CoarsePre) THEN
    CALL DivInterface(b,y)
    DO ib=1,nb
      y(ib)=y(ib)*Floor(ib)%Boundary 
    END DO
    CALL dposl(LaplNeumGlob,nb,nb,y)
    CALL ProjectInterface(y,b)
  END IF  
END SUBROUTINE Project

SUBROUTINE ProjectT(b)

  TYPE(PressureVelocity) :: b(:)
  
  REAL(RealKind) :: y(nb)
  INTEGER :: i

! Projektion
  IF (CoarsePre) THEN
    CALL DivInterfaceT(b,y)
    DO ib=1,nb
      y(ib)=-y(ib)*Floor(ib)%Boundary 
    END DO
    CALL dposl(LaplNeumGlob,nb,nb,y)
    CALL ProjectInterfaceT(y,b)
  END IF  
END SUBROUTINE ProjectT

SUBROUTINE ProjectInterface(y,x)
  
  TYPE(PressureVelocity), TARGET :: x(:)
  REAL(RealKind) :: y(nb)

  INTEGER :: in
  INTEGER :: jx,jy,jz
  REAL(RealKind) :: Temp
  TYPE(PressureVelocity), POINTER :: XVec
  TYPE(Nachbar_T), POINTER   :: Nachbar
  REAL(RealKind) :: GradY
  
! -- Aufdatieren von x = ugamma = ugamma - Dtrans * y auf dem Interface --

! -- alle lokalen Bloecke, alle Nachbarn --
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    XVec => x(ibLoc)

    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      GradY=y(ib)-y(ibn)
!  -- Westlicher Rand --     
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              XVec%u_w(jy,jz)=(XVec%u_w(jy,jz)+dtP*beta0*GradY*FUG(ix0,jy,jz)*DUU(ix0,jy,jz))/(VolFace(ibLoc)%u_w(jy,jz)+Eps) 
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%West==1) THEN 
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              XVec%u_w(jy,jz)=(XVec%u_w(jy,jz)+dtP*beta0*y(ib)*FUG(ix0,jy,jz)*DUU(ix0,jy,jz))/(VolFace(ibLoc)%u_w(jy,jz)+Eps) 
            END DO
          END DO
        ELSE 
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FUG(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
              XVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)=(XVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)+Temp)/ &
                                                    (VolFace(ibLoc)%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)+Eps)
            END DO
          END DO
        END IF
      END IF
           
!  -- Oestlicher Rand --     
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              XVec%u_e(jy,jz)=(XVec%u_e(jy,jz)-dtP*beta0*GradY*FUG(ix1,jy,jz)*DUU(ix1,jy,jz))/(VolFace(ibLoc)%u_e(jy,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%East==1) THEN 
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              XVec%u_e(jy,jz)=(XVec%u_e(jy,jz)-dtP*beta0*y(ib)*FUG(ix1,jy,jz)*DUU(ix1,jy,jz))/(VolFace(ibLoc)%u_e(jy,jz)+Eps) 
            END DO
          END DO
        ELSE
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FUG(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
              XVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)=XVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)-Temp
            END DO
          END DO
        END IF
      END IF
!  -- Suedlicher Rand --     
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              XVec%v_s(jx,jz)=(XVec%v_s(jx,jz)+dtP*beta0*GradY*FVG(jx,iy0,jz)*DUV(jx,iy0,jz))/(VolFace(ibLoc)%v_s(jx,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%South==1) THEN 
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              XVec%v_s(jx,jz)=(XVec%v_s(jx,jz)+dtP*beta0*y(ib)*FVG(jx,iy0,jz)*DUV(jx,iy0,jz))/(VolFace(ibLoc)%v_s(jx,jz)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FVG(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
              XVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)=XVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)+Temp
            END DO
          END DO
        END IF
      END IF
! --  Noerdlicher Rand --     
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              XVec%v_n(jx,jz)=(XVec%v_n(jx,jz)-dtP*beta0*GradY*FVG(jx,iy1,jz)*DUV(jx,iy1,jz))/(VolFace(ibLoc)%v_n(jx,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%North==1) THEN 
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              XVec%v_n(jx,jz)=(XVec%v_n(jx,jz)-dtP*beta0*y(ib)*FVG(jx,iy1,jz)*DUV(jx,iy1,jz))/(VolFace(ibLoc)%v_n(jx,jz)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FVG(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
              XVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)=XVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)-Temp
            END DO
          END DO
        END IF
      END IF
! --  Unterer Rand --     
      IF (Nachbar%nType(2:2) == 'b') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_b(jx,jy)=(XVec%w_b(jx,jy)+dtP*beta0*GradY*FWG(jx,jy,iz0)*DUW(jx,jy,iz0))/(VolFace(ibLoc)%w_b(jx,jy)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Bottom==1) THEN 
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_b(jx,jy)=(XVec%w_b(jx,jy)+dtP*beta0*y(ib)*FWG(jx,jy,iz0)*DUW(jx,jy,iz0))/(VolFace(ibLoc)%w_b(jx,jy)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              Temp=GradY*SUM(FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
              XVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)=XVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)+Temp
            END DO
          END DO
        END IF
      END IF
! --  Oberer Rand --     
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_t(jx,jy)=(XVec%w_t(jx,jy)-dtP*beta0*GradY*FWG(jx,jy,iz1)*DUW(jx,jy,iz1))/(VolFace(ibLoc)%w_t(jx,jy)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Top==1) THEN 
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_t(jx,jy)=(XVec%w_t(jx,jy)-dtP*beta0*y(ib)*FWG(jx,jy,iz1)*DUW(jx,jy,iz1))/(VolFace(ibLoc)%w_t(jx,jy)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              Temp=GradY*SUM(FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
              XVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)=XVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)-Temp
            END DO
          END DO
        END IF
      END IF
    END DO     ! in
  END DO       ! ibLoc
    
END SUBROUTINE ProjectInterface
  
SUBROUTINE ProjectInterfaceT(y,x)
  
  TYPE(PressureVelocity), TARGET :: x(:)
  REAL(RealKind) :: y(nb)

  INTEGER :: in
  INTEGER :: jx,jy,jz
  REAL(RealKind) :: Temp
  TYPE(PressureVelocity), POINTER :: XVec
  TYPE(Nachbar_T), POINTER   :: Nachbar
  REAL(RealKind) :: GradY
  
! -- Aufdatieren von x = ugamma = ugamma - Dtrans * y auf dem Interface --

! -- alle lokalen Bloecke, alle Nachbarn --
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    XVec => x(ibLoc)

    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      GradY=y(ib)-y(ibn)
!  -- Westlicher Rand --     
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              XVec%u_w(jy,jz) = XVec%u_w(jy,jz) +      &
                  dtP*beta0*GradY * FUG(ix0,jy,jz)*DUU(ix0,jy,jz)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%West==1) THEN 
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              XVec%u_w(jy,jz) = XVec%u_w(jy,jz) +      &
                  dtP*beta0*y(ib) * FUG(ix0,jy,jz)*DUU(ix0,jy,jz)
            END DO
          END DO
        ELSE 
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FUG(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
              XVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)=XVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)+Temp
            END DO
          END DO
        END IF
      END IF
           
!  -- Oestlicher Rand --     
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              XVec%u_e(jy,jz) = XVec%u_e(jy,jz) -      &
                  dtP*beta0*GradY * FUG(ix1,jy,jz)*DUU(ix1,jy,jz)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%East==1) THEN 
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              XVec%u_e(jy,jz) = XVec%u_e(jy,jz) -      &
                  dtP*beta0*y(ib) * FUG(ix1,jy,jz)*DUU(ix1,jy,jz)
            END DO
          END DO
        ELSE
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FUG(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
              XVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)=XVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)-Temp
            END DO
          END DO
        END IF
      END IF
!  -- Suedlicher Rand --     
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              XVec%v_s(jx,jz) = XVec%v_s(jx,jz) +      &
                  dtP*beta0*GradY * FVG(jx,iy0,jz)*DUV(jx,iy0,jz)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%South==1) THEN 
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              XVec%v_s(jx,jz) = XVec%v_s(jx,jz) +      &
                  dtP*beta0*y(ib) * FVG(jx,iy0,jz)*DUV(jx,iy0,jz)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FVG(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
              XVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)=XVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)+Temp
            END DO
          END DO
        END IF
      END IF
! --  Noerdlicher Rand --     
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              XVec%v_n(jx,jz) = XVec%v_n(jx,jz) -      &
                  dtP*beta0*GradY * FVG(jx,iy1,jz)*DUV(jx,iy1,jz)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%North==1) THEN 
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              XVec%v_n(jx,jz) = XVec%v_n(jx,jz) -      &
                  dtP*beta0*y(ib) * FVG(jx,iy1,jz)*DUV(jx,iy1,jz)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FVG(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
              XVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)=XVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)-Temp
            END DO
          END DO
        END IF
      END IF
! --  Unterer Rand --     
      IF (Nachbar%nType(2:2) == 'b') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_b(jx,jy) = XVec%w_b(jx,jy) +      &
                  dtP*beta0*GradY * FWG(jx,jy,iz0)*DUW(jx,jy,iz0)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Bottom==1) THEN 
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_b(jx,jy) = XVec%w_b(jx,jy) +      &
                  dtP*beta0*y(ib) * FWG(jx,jy,iz0)*DUW(jx,jy,iz0)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              Temp=GradY*SUM(FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
              XVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)=XVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)+Temp
            END DO
          END DO
        END IF
      END IF
! --  Oberer Rand --     
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_t(jx,jy) = XVec%w_t(jx,jy) -      &
                  dtP*beta0*GradY * FWG(jx,jy,iz1)*DUW(jx,jy,iz1)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Top==1) THEN 
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_t(jx,jy) = XVec%w_t(jx,jy) -      &
                  dtP*beta0*y(ib) * FWG(jx,jy,iz1)*DUW(jx,jy,iz1)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              Temp=GradY*SUM(FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
              XVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)=XVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)-Temp
            END DO
          END DO
        END IF
      END IF
    END DO     ! in
  END DO       ! ibLoc
    
END SUBROUTINE ProjectInterfaceT
  
SUBROUTINE b1pBTx2(b,x)

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
    IF (MultiMu.OR.MultiTriT) THEN
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

END SUBROUTINE b1pBTx2

SUBROUTINE b2pBx1(b,x)

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
    IF (MultiMu.OR.MultiTriT) THEN
      IF (TypeW/='ow') THEN 
        BVec%u_w(:,:)=Fac5*XVec%u_w(:,:)           &
                      +beta0*dtP*FUG(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)
      ELSE IF (TypeW=='ow'.AND.BCP%West==1) THEN 
        BVec%u_w(:,:)=XVec%u_w(:,:)           &
                      +beta0*dtP*FUG(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)
      END IF
      IF (TypeE/='oe') THEN
        BVec%u_e(:,:)=Fac5*XVec%u_e(:,:)           &
                      -beta0*dtP*FUG(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)
      ELSE IF (TypeE=='oe'.AND.BCP%East==1) THEN
        BVec%u_e(:,:)=XVec%u_e(:,:)           &
                      -beta0*dtP*FUG(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)
      END IF    
      IF (TypeS/='os') THEN
        BVec%v_s(:,:)=Fac5*XVec%v_s(:,:)           &
                      +beta0*dtP*FVG(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)
      ELSE IF (TypeS=='os'.AND.BCP%South==1) THEN
        BVec%v_s(:,:)=XVec%v_s(:,:)           &
                      +beta0*dtP*FVG(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)
      END IF
      IF (TypeN/='on') THEN     
        BVec%v_n(:,:)=Fac5*XVec%v_n(:,:)           &
                      -beta0*dtP*FVG(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)
      ELSE IF (TypeN=='on'.AND.BCP%North==1) THEN     
        BVec%v_n(:,:)=XVec%v_n(:,:)           &
                      -beta0*dtP*FVG(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)
      END IF
      IF (TypeB/='ob') THEN     
        BVec%w_b(:,:)=Fac5*XVec%w_b(:,:)           &
                      +beta0*dtP*FWG(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)
      ELSE IF (TypeB=='ob'.AND.BCP%Bottom==1) THEN     
        BVec%w_b(:,:)=XVec%w_b(:,:)           &
                      +beta0*dtP*FWG(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1)
      END IF
      IF (TypeT/='ot') THEN     
        BVec%w_t(:,:)=Fac5*XVec%w_t(:,:)           &
                      -beta0*dtP*FWG(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)
      ELSE IF (TypeT=='ot'.AND.BCP%Top==1) THEN     
        BVec%w_t(:,:)=XVec%w_t(:,:)           &
                      -beta0*dtP*FWG(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1)
      END IF
    ELSE IF (MultiTriTB) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%u_w(:,:)=XVec%u_w(:,:)           &
                      +beta0*dtP*FUG(ix0,:,:)*DUT(ix0+1,:,:)*DUU(ix0,:,:)*BVec%p(1,ix0+1,:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%u_e(:,:)=XVec%u_e(:,:)           &
                      -beta0*dtP*FUG(ix1,:,:)*DUT(ix1,:,:)*DUU(ix1,:,:)*BVec%p(1,ix1,:,:)
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%v_s(:,:)=XVec%v_s(:,:)           &
                      +beta0*dtP*FVG(:,iy0,:)*DUT(:,iy0+1,:)*DUV(:,iy0,:)*BVec%p(1,:,iy0+1,:)
      END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%v_n(:,:)=XVec%v_n(:,:)           &
                      -beta0*dtP*FVG(:,iy1,:)*DUT(:,iy1,:)*DUV(:,iy1,:)*BVec%p(1,:,iy1,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%w_b(:,:)=XVec%w_b(:,:)           &
                      +beta0*dtP*FWG(:,:,iz0)*DUT(:,:,iz0+1)*DUW(:,:,iz0)*BVec%p(1,:,:,iz0+1) &
                      -Half*beta0*dtP*GravComp*VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)*BVec%p(2,:,:,iz0+1)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%w_t(:,:)=XVec%w_t(:,:)           &
                      -beta0*dtP*FWG(:,:,iz1)*DUT(:,:,iz1)*DUW(:,:,iz1)*BVec%p(1,:,:,iz1) &
                      -Half*beta0*dtP*GravComp*VolC(ix0+1:ix1,iy0+1:iy1,iz1)*BVec%p(2,:,:,iz1)
      END IF
    ELSE IF (MultiTriTR.OR.MultiMuTR) THEN
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        BVec%u_w(:,:)=XVec%u_w(:,:)           &
                      +beta0*dtP*DUU(ix0,:,:)*FUG(ix0,:,:)*DUT(ix0+1,:,:)*BVec%p(1,ix0+1,:,:) &
                      +beta0*dtP*DUU(ix0,:,:)*FUG(ix0,:,:)*DUR(ix0+1,:,:)*BVec%p(2,ix0+1,:,:)
      END IF
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        BVec%u_e(:,:)=XVec%u_e(:,:)           &
                      -beta0*dtP*DUU(ix1,:,:)*FUG(ix1,:,:)*DUT(ix1,:,:)*BVec%p(1,ix1,:,:) &
                      -beta0*dtP*DUU(ix1,:,:)*FUG(ix1,:,:)*DUR(ix1,:,:)*BVec%p(2,ix1,:,:)                      
      END IF    
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        BVec%v_s(:,:)=XVec%v_s(:,:)           &
                      +beta0*dtP*DUV(:,iy0,:)*FVG(:,iy0,:)*DUT(:,iy0+1,:)*BVec%p(1,:,iy0+1,:) &
                      +beta0*dtP*DUV(:,iy0,:)*FVG(:,iy0,:)*DUR(:,iy0+1,:)*BVec%p(2,:,iy0+1,:)
     END IF
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        BVec%v_n(:,:)=XVec%v_n(:,:)           &
                      -beta0*dtP*DUV(:,iy1,:)*FVG(:,iy1,:)*DUT(:,iy1,:)*BVec%p(1,:,iy1,:) &
                      -beta0*dtP*DUV(:,iy1,:)*FVG(:,iy1,:)*DUR(:,iy1,:)*BVec%p(2,:,iy1,:)
      END IF
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        BVec%w_b(:,:)=XVec%w_b(:,:)           &
                      +beta0*dtP*DUW(:,:,iz0)*FWG(:,:,iz0)*DUT(:,:,iz0+1)*BVec%p(1,:,:,iz0+1) &
                      +beta0*dtP*DUW(:,:,iz0)*FWG(:,:,iz0)*DUR(:,:,iz0+1)*BVec%p(2,:,:,iz0+1) &
                      -Half*beta0*dtP*GravComp*VolC(ix0+1:ix1,iy0+1:iy1,iz0+1)*BVec%p(2,:,:,iz0+1)
      END IF
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        BVec%w_t(:,:)=XVec%w_t(:,:)           &
                      -beta0*dtP*DUW(:,:,iz1)*FWG(:,:,iz1)*DUT(:,:,iz1)*BVec%p(1,:,:,iz1) &
                      -beta0*dtP*DUW(:,:,iz1)*FWG(:,:,iz1)*DUR(:,:,iz1)*BVec%p(2,:,:,iz1) &
                      -Half*beta0*dtP*GravComp*VolC(ix0+1:ix1,iy0+1:iy1,iz1)*BVec%p(2,:,:,iz1)
      END IF
    END IF
  END DO   ! ibLoc

END SUBROUTINE b2pBx1

SUBROUTINE PSolve11(b,x)

  TYPE(PressureVelocity) :: b(:),x(:) 

  b=x
  CALL PressurePre(b,x)
  CALL b2pBx1(b,x)
! -- Randwert-Austausch 
  CALL Exchange(b)
! Projektion
! CALL ProjectT(b)
  CALL Project(b)
  CALL SchurPre(b)
!   -- Randwert-Austausch 
  CALL Exchange(b)
!   Projektion
  CALL Project(b)
  CALL b1pBTx2(b,x) 
  CALL PressurePre(b,b) 

END SUBROUTINE PSolve11

SUBROUTINE PSolve1(b,x)

  TYPE(PressureVelocity) :: b(:),x(:) 

  b=x
  CALL PressurePre(b,x)
  CALL b2pBx1(b,x)
! -- Randwert-Austausch 
  CALL Exchange(b)
  CALL SchurPre(b)
!   -- Randwert-Austausch 
  CALL Exchange(b)
  CALL b1pBTx2(b,x) 
  CALL PressurePre(b,b) !OSWALD
! CALL Coarse(b,x)

END SUBROUTINE PSolve1

SUBROUTINE Coarse(b,x)
  TYPE(PressureVelocity) :: b(:),x(:) 

  REAL(RealKind) :: y(nb)

  CALL CoarseInterface(x,y)
  CALL dposl(LaplNeumGlob,nb,nb,y)
  CALL CoarseProject(b,x,y)

END SUBROUTINE Coarse

SUBROUTINE CoarseProject(b,x,y)
  
  TYPE(PressureVelocity), TARGET :: b(:),x(:)
  REAL(RealKind) :: y(nb)

  INTEGER :: in
  INTEGER :: ix,iy,iz
  INTEGER :: jx,jy,jz
  REAL(RealKind) :: Temp
  TYPE(PressureVelocity), POINTER :: XVec
  TYPE(PressureVelocity), POINTER :: BVec
  TYPE(Nachbar_T), POINTER   :: Nachbar
  REAL(RealKind) :: GradY
  
! -- Aufdatieren von x = ugamma = ugamma - Dtrans * y auf dem Interface --

! -- alle lokalen Bloecke, alle Nachbarn --
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    XVec => x(ibLoc)
    BVec => b(ibLoc)

    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1-1
          BVec%p(1,ix,iy,iz)=BVec%p(1,ix,iy,iz)+y(ib)/(DUTG(ibLoc)%p(ix,iy,iz)+Eps)
        END DO
      END DO
    END DO
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      GradY=y(ib)-y(ibn)
!  -- Westlicher Rand --     
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jz=jz0+1,jz1
            DO jy=jy0+1,jy1
              BVec%u_w(jy,jz)=BVec%u_w(jy,jz)+(XVec%u_w(jy,jz)+dtP*beta0*GradY*FUG(ix0,jy,jz)*DUU(ix0,jy,jz))/(VolFace(ibLoc)%u_w(jy,jz)+Eps) 
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%West==1) THEN 
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              BVec%u_w(jy,jz)=BVec%u_w(jy,jz)+(XVec%u_w(jy,jz)+dtP*beta0*y(ib)*FUG(ix0,jy,jz)*DUU(ix0,jy,jz))/(VolFace(ibLoc)%u_w(jy,jz)+Eps) 
            END DO
          END DO
        ELSE 
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FUG(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)=BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)+ &
                                                    (XVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)+Temp)/ &
                                                    (VolFace(ibLoc)%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)+Eps)
            END DO
          END DO
        END IF
      END IF
           
!  -- Oestlicher Rand --     
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              BVec%u_e(jy,jz)=BVec%u_e(jy,jz)+(XVec%u_e(jy,jz)-dtP*beta0*GradY*FUG(ix1,jy,jz)*DUU(ix1,jy,jz))/(VolFace(ibLoc)%u_e(jy,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%East==1) THEN 
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              BVec%u_e(jy,jz)=BVec%u_e(jy,jz)+(XVec%u_e(jy,jz)-dtP*beta0*y(ib)*FUG(ix1,jy,jz)*DUU(ix1,jy,jz))/(VolFace(ibLoc)%u_e(jy,jz)+Eps) 
            END DO
          END DO
        ELSE
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FUG(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)=XVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)-Temp
            END DO
          END DO
        END IF
      END IF
!  -- Suedlicher Rand --     
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              BVec%v_s(jx,jz)=BVec%v_s(jx,jz)+(XVec%v_s(jx,jz)+dtP*beta0*GradY*FVG(jx,iy0,jz)*DUV(jx,iy0,jz))/(VolFace(ibLoc)%v_s(jx,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%South==1) THEN 
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              BVec%v_s(jx,jz)=BVec%v_s(jx,jz)+(XVec%v_s(jx,jz)+dtP*beta0*y(ib)*FVG(jx,iy0,jz)*DUV(jx,iy0,jz))/(VolFace(ibLoc)%v_s(jx,jz)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FVG(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)=XVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)+Temp
            END DO
          END DO
        END IF
      END IF
! --  Noerdlicher Rand --     
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              BVec%v_n(jx,jz)=BVec%v_n(jx,jz)+(XVec%v_n(jx,jz)-dtP*beta0*GradY*FVG(jx,iy1,jz)*DUV(jx,iy1,jz))/(VolFace(ibLoc)%v_n(jx,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%North==1) THEN 
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              BVec%v_n(jx,jz)=BVec%v_n(jx,jz)+(XVec%v_n(jx,jz)-dtP*beta0*y(ib)*FVG(jx,iy1,jz)*DUV(jx,iy1,jz))/(VolFace(ibLoc)%v_n(jx,jz)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              Temp=GradY*SUM(FVG(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)=XVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)-Temp
            END DO
          END DO
        END IF
      END IF
! --  Unterer Rand --     
      IF (Nachbar%nType(2:2) == 'b') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_b(jx,jy)=(XVec%w_b(jx,jy)+dtP*beta0*GradY*FWG(jx,jy,iz0)*DUW(jx,jy,iz0))/(VolFace(ibLoc)%w_b(jx,jy)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Bottom==1) THEN 
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_b(jx,jy)=(XVec%w_b(jx,jy)+dtP*beta0*y(ib)*FWG(jx,jy,iz0)*DUW(jx,jy,iz0))/(VolFace(ibLoc)%w_b(jx,jy)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              Temp=GradY*SUM(FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))
              XVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)=XVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)+Temp
            END DO
          END DO
        END IF
      END IF
! --  Oberer Rand --     
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_t(jx,jy)=(XVec%w_t(jx,jy)-dtP*beta0*GradY*FWG(jx,jy,iz1)*DUW(jx,jy,iz1))/(VolFace(ibLoc)%w_t(jx,jy)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Top==1) THEN 
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              XVec%w_t(jx,jy)=(XVec%w_t(jx,jy)-dtP*beta0*y(ib)*FWG(jx,jy,iz1)*DUW(jx,jy,iz1))/(VolFace(ibLoc)%w_t(jx,jy)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              Temp=GradY*SUM(FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))
              XVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)=XVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)-Temp
            END DO
          END DO
        END IF
      END IF
    END DO     ! in
  END DO       ! ibLoc
    
END SUBROUTINE CoarseProject

SUBROUTINE CoarseInterface(u,y)
  
  INTEGER :: in,jx,jy,jz

  TYPE(PressureVelocity) :: u(:) 
  REAL(RealKind) :: y(nb)

  DO ib=1,nb
    IF (MyId==blMPI(ib)%Proc) THEN
      y(ib)=0.0d0
      CALL Set(Floor(ib))
      ibLoc=blMPI(ib)%ibLoc
      DTU=>DTUG(ibLoc)%uF
      DTV=>DTUG(ibLoc)%vF
      DTW=>DTUG(ibLoc)%wF
      DTT=>DTTG(ibLoc)%p
! --  Westlicher Rand --     
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        DO jy=iy0+1,iy1
          DO jz=iz0+1,iz1
            y(ib)=y(ib)-u(ibLoc)%u_w(jy,jz)*FU(ix0,jy,jz)*DTU(ix0,jy,jz)*DTT(ix0+1,jy,jz)/ &
                        (VolFace(ibLoc)%u_w(jy,jz)+Eps)
          END DO
        END DO
      END IF
! --  Oestlicher Rand --     
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        DO jy=iy0+1,iy1
          DO jz=iz0+1,iz1
            y(ib)=y(ib)+u(ibLoc)%u_e(jy,jz)*FU(ix1,jy,jz)*DTU(ix1,jy,jz)*DTT(ix1,jy,jz)/ &
                        (VolFace(ibLoc)%u_w(jy,jz)+Eps)
          END DO
        END DO
      END IF
! --  Suedlicher Rand --     
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        DO jx=ix0+1,ix1
          DO jz=iz0+1,iz1
            y(ib)=y(ib)-u(ibLoc)%v_s(jx,jz)*FV(jx,iy0,jz)*DTV(jx,iy0,jz)*DTT(jx,iy0+1,jz)/ &
                        (VolFace(ibLoc)%v_s(jx,jz)+Eps)
           END DO
        END DO
      END IF
! --  Noerdlicher Rand --     
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        DO jx=ix0+1,ix1
          DO jz=iz0+1,iz1
            y(ib)=y(ib)+u(ibLoc)%v_n(jx,jz)*FV(jx,iy1,jz)*DTV(jx,iy1,jz)*DTT(jx,iy1,jz)/ &
                        (VolFace(ibLoc)%v_n(jx,jz)+Eps)
          END DO
        END DO
      END IF
! --  Unterer Rand --     
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        DO jx=ix0+1,ix1
          DO jy=iy0+1,iy1
            y(ib)=y(ib)-u(ibLoc)%w_b(jx,jy)*FW(jx,jy,iz0)*DTW(jx,jy,iz0)*DTT(jx,jy,iz0+1)/ &
                        (VolFace(ibLoc)%w_b(jy,jz)+Eps)
          END DO
        END DO
      END IF
           
! --  Oberer Rand --     
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        DO jx=ix0+1,ix1
          DO jy=iy0+1,iy1
            y(ib)=y(ib)+u(ibLoc)%w_t(jx,jy)*FW(jx,jy,iz1)*DTW(jx,jy,iz1)*DTT(jx,jy,iz1)/ &
                        (VolFace(ibLoc)%w_t(jy,jz)+Eps)
          END DO
        END DO
      END IF
      y(ib)=SUM(u(ibLoc)%p)-beta0*dtP*y(ib)
    END IF  ! MyId
! -- Senden der Komponente ib an die anderen Prozesse --
    CALL MPI_Bcast(y(ib),1,MPI_RealKind,blMPI(ib)%Proc, &
 &                MPI_COMM_WORLD,MPIErr)
  END DO     ! ib
    
END SUBROUTINE CoarseInterface

SUBROUTINE DivInterface(u,y)
  
  INTEGER :: in,jx,jy,jz

  TYPE(PressureVelocity) :: u(:) 
  REAL(RealKind) :: y(nb)

  
! -- Berechnung des Matrix-Vektor-Produkts y = D * x auf dem Interface Gamma --

  DO ib=1,nb
    IF (MyId == blMPI(ib)%Proc) THEN
      y(ib) = 0.0e0
      CALL Set(Floor(ib))
      ibLoc=blMPI(ib)%ibLoc
      DTU=>DTUG(ibLoc)%uF
      DTV=>DTUG(ibLoc)%vF
      DTW=>DTUG(ibLoc)%wF
      DTT=>DTTG(ibLoc)%p

! --  Westlicher Rand --     
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
        DO jy=iy0+1,iy1
          DO jz=iz0+1,iz1
            y(ib)=y(ib)-u(ibLoc)%u_w(jy,jz)*FU(ix0,jy,jz)*DTU(ix0,jy,jz)*DTT(ix0+1,jy,jz)/ &
                        (VolFace(ibLoc)%u_w(jy,jz)+Eps)
          END DO
        END DO
      END IF
           
! --  Oestlicher Rand --     
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
        DO jy=iy0+1,iy1
          DO jz=iz0+1,iz1
            y(ib)=y(ib)+u(ibLoc)%u_e(jy,jz)*FU(ix1,jy,jz)*DTU(ix1,jy,jz)*DTT(ix1,jy,jz)/ &
                        (VolFace(ibLoc)%u_w(jy,jz)+Eps)
          END DO
        END DO
      END IF
           
! --  Suedlicher Rand --     
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
        DO jx=ix0+1,ix1
          DO jz=iz0+1,iz1
            y(ib)=y(ib)-u(ibLoc)%v_s(jx,jz)*FV(jx,iy0,jz)*DTV(jx,iy0,jz)*DTT(jx,iy0+1,jz)/ &
                        (VolFace(ibLoc)%v_s(jx,jz)+Eps)
           END DO
        END DO
      END IF
           
! --  Noerdlicher Rand --     
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
        DO jx=ix0+1,ix1
          DO jz=iz0+1,iz1
            y(ib)=y(ib)+u(ibLoc)%v_n(jx,jz)*FV(jx,iy1,jz)*DTV(jx,iy1,jz)*DTT(jx,iy1,jz)/ &
                        (VolFace(ibLoc)%v_n(jx,jz)+Eps)
          END DO
        END DO
      END IF
            
! --  Unterer Rand --     
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
        DO jx=ix0+1,ix1
          DO jy=iy0+1,iy1
            y(ib)=y(ib)-u(ibLoc)%w_b(jx,jy)*FW(jx,jy,iz0)*DTW(jx,jy,iz0)*DTT(jx,jy,iz0+1)/ &
                        (VolFace(ibLoc)%w_b(jy,jz)+Eps)
          END DO
        END DO
      END IF
           
! --  Oberer Rand --     
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
        DO jx=ix0+1,ix1
          DO jy=iy0+1,iy1
            y(ib)=y(ib)+u(ibLoc)%w_t(jx,jy)*FW(jx,jy,iz1)*DTW(jx,jy,iz1)*DTT(jx,jy,iz1)/ &
                        (VolFace(ibLoc)%w_t(jy,jz)+Eps)
          END DO
        END DO
      END IF
    END IF  ! MyId

! -- Senden der Komponente ib an die anderen Prozesse --
    CALL MPI_Bcast(y(ib),1,MPI_RealKind,blMPI(ib)%Proc, &
 &                MPI_COMM_WORLD,MPIErr)
  END DO     ! ib
  y=y*dtP*beta0
    
END SUBROUTINE DivInterface

SUBROUTINE DivInterfaceT(u,y)
  
  INTEGER :: in,jx,jy,jz

  TYPE(PressureVelocity) :: u(:) 
  REAL(RealKind) :: y(nb)

  
! -- Berechnung des Matrix-Vektor-Produkts y = D * x auf dem Interface Gamma --

  DO ib=1,nb
    IF (MyId == blMPI(ib)%Proc) THEN
      y(ib) = 0.0e0
      CALL Set(Floor(ib))
      ibLoc=blMPI(ib)%ibLoc
      DTU=>DTUG(ibLoc)%uF
      DTV=>DTUG(ibLoc)%vF
      DTW=>DTUG(ibLoc)%wF
      DTT=>DTTG(ibLoc)%p

! -- Westlicher Rand --     
      IF (TypeW/='ow'.OR.(TypeW=='ow'.AND.BCP%West==1)) THEN 
       DO jy=iy0+1,iy1
         DO jz=iz0+1,iz1
           y(ib)=y(ib)-u(ibLoc)%u_w(jy,jz)*FU(ix0,jy,jz)*DTU(ix0,jy,jz)*DTT(ix0+1,jy,jz)/(VolFace(ibLoc)%u_w(jy,jz)+Eps)
         END DO
       END DO
     END IF
           
! -- Oestlicher Rand --     
      IF (TypeE/='oe'.OR.(TypeE=='oe'.AND.BCP%East==1)) THEN
       DO jy=iy0+1,iy1
         DO jz=iz0+1,iz1
           y(ib)=y(ib)+u(ibLoc)%u_e(jy,jz)*FU(ix1,jy,jz)*DTU(ix1,jy,jz)*DTT(ix1,jy,jz)/(VolFace(ibLoc)%u_e(jy,jz)+Eps)
         END DO
       END DO
     END IF
           
! -- Suedlicher Rand --     
      IF (TypeS/='os'.OR.(TypeS=='os'.AND.BCP%South==1)) THEN
       DO jx=ix0+1,ix1
         DO jz=iz0+1,iz1
           y(ib)=y(ib)-u(ibLoc)%v_s(jx,jz)*FV(jx,iy0,jz)*DTV(jx,iy0,jz)*DTT(jx,iy0+1,jz)/(VolFace(ibLoc)%v_s(jx,jz)+Eps)
          END DO
       END DO
     END IF
           
! -- Noerdlicher Rand --     
      IF (TypeN/='on'.OR.(TypeN=='on'.AND.BCP%North==1)) THEN     
       DO jx=ix0+1,ix1
         DO jz=iz0+1,iz1
           y(ib)=y(ib)+u(ibLoc)%v_n(jx,jz)*FV(jx,iy1,jz)*DTV(jx,iy1,jz)*DTT(jx,iy1,jz)/(VolFace(ibLoc)%v_n(jx,jz)+Eps)
         END DO
       END DO
     END IF
           
! -- Unterer Rand --     
      IF (TypeB/='ob'.OR.(TypeB=='ob'.AND.BCP%Bottom==1)) THEN     
           DO jx=ix0+1,ix1
              DO jy=iy0+1,iy1
                 y(ib) = y(ib) - u(ibLoc)%w_b(jx,jy) *  &
                                   FW(jx,jy,iz0)*DTW(jx,jy,iz0)*DTT(jx,jy,iz0+1)/(VolFace(ibLoc)%w_b(jx,jy)+Eps)
              END DO
           END DO
        END IF
           
! -- Oberer Rand --     
      IF (TypeT/='ot'.OR.(TypeT=='ot'.AND.BCP%Top==1)) THEN     
           DO jx=ix0+1,ix1
              DO jy=iy0+1,iy1
                 y(ib) = y(ib) + u(ibLoc)%w_t(jx,jy) *  &
                                   FW(jx,jy,iz1)*DTW(jx,jy,iz1)*DTT(jx,jy,iz1)/(VolFace(ibLoc)%w_t(jx,jy)+Eps)
              END DO
           END DO
        END IF
           
     END IF  ! MyId

! -- Senden der Komponente ib an die anderen Prozesse --
     CALL MPI_Bcast(y(ib),1,MPI_RealKind,blMPI(ib)%Proc, &
 &                  MPI_COMM_WORLD,MPIErr)

  END DO     ! ib
  y=y*dtP*beta0
    
END SUBROUTINE DivInterfaceT

END MODULE PreCond_Mod
