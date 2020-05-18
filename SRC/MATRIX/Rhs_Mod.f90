MODULE Rhs_Mod
  USE Domain_Mod
  USE DataType_Mod
  USE BoundaryCond_Mod
  USE JacAccGrav_Mod
  USE Tools_Mod
  USE Chemie_Mod
  
  IMPLICIT NONE

CONTAINS

SUBROUTINE Divergence(Div,n,m,l,WeiFU,WeiFV,WeiFW,u,v,w,DTT) !MJ

  IMPLICIT NONE

  INTEGER :: n,m,l
  REAL(RealKind) :: Div(1:n,1:m,1:l)
  REAL(RealKind) :: WeiFU(0:n,1:m,1:l)
  REAL(RealKind) :: WeiFV(1:n,0:m,1:l)
  REAL(RealKind) :: WeiFW(1:n,1:m,0:l)
  REAL(RealKind) :: u(0:n,1:m,1:l)
  REAL(RealKind) :: v(1:n,0:m,1:l)
  REAL(RealKind) :: w(1:n,1:m,0:l)
  REAL(RealKind) :: DTT(1:n,1:m,1:l)

  INTEGER :: i,j,k

  DO k=1,l
    DO j=1,m
      DO i=1,n
        Div(i,j,k)=(-WeiFU(i-1,j,k)*u(i-1,j,k) &
                   +WeiFU(i  ,j,k)*u(i  ,j,k) &
                   -WeiFV(i,j-1,k)*v(i,j-1,k) &
                   +WeiFV(i  ,j,k)*v(i  ,j,k) &
                   -WeiFW(i,j,k-1)*w(i,j,k-1) &
                   +WeiFW(i  ,j,k)*w(i  ,j,k))*DTT(i,j,k) 
                   
      END DO
    END DO
  END DO  
  
END SUBROUTINE Divergence


SUBROUTINE DivergenceS(Div,n,m,l, &
                       WeiFU,WeiFV,WeiFW, &
                       DTU,DTV,DTW, &
                       u,v,w,DTT)

  IMPLICIT NONE

  INTEGER :: n,m,l
  REAL(RealKind) :: Div(1:n,1:m,1:l)
  REAL(RealKind) :: WeiFU(0:n,1:m,1:l)
  REAL(RealKind) :: WeiFV(1:n,0:m,1:l)
  REAL(RealKind) :: WeiFW(1:n,1:m,0:l)
  REAL(RealKind) :: DTU(0:n,1:m,1:l)
  REAL(RealKind) :: DTV(1:n,0:m,1:l)
  REAL(RealKind) :: DTW(1:n,1:m,0:l)
  REAL(RealKind) :: u(0:n,1:m,1:l)
  REAL(RealKind) :: v(1:n,0:m,1:l)
  REAL(RealKind) :: w(1:n,1:m,0:l)
  REAL(RealKind) :: DTT(1:n,1:m,1:l)

  INTEGER :: i,j,k

  DO k=1,l
    DO j=1,m
      DO i=1,n
        Div(i,j,k)=(-WeiFU(i-1,j,k)*DTU(i-1,j,k)*u(i-1,j,k) &
                   +WeiFU(i  ,j,k)*DTU(i  ,j,k)*u(i  ,j,k) &
                   -WeiFV(i,j-1,k)*DTV(i,j-1,k)*v(i,j-1,k) &
                   +WeiFV(i  ,j,k)*DTV(i,  j,k)*v(i  ,j,k) &
                   -WeiFW(i,j,k-1)*DTW(i,j,k-1)*w(i,j,k-1) &
                   +WeiFW(i  ,j,k)*DTW(i,j,  k)*w(i  ,j,k) &
                   )*DTT(i,j,k)
      END DO
    END DO
  END DO  
  
END SUBROUTINE DivergenceS
  
SUBROUTINE DivInter(Div,n,m,l, &
                    WeiFU,WeiFV,WeiFW, &
                    DTU,DTV,DTW, &
                    u,v,w)

  IMPLICIT NONE

  INTEGER :: n,m,l
  REAL(RealKind) :: Div
  REAL(RealKind) :: WeiFU(0:n,1:m,1:l)
  REAL(RealKind) :: WeiFV(1:n,0:m,1:l)
  REAL(RealKind) :: WeiFW(1:n,1:m,0:l)
  REAL(RealKind) :: DTU(0:n,1:m,1:l)
  REAL(RealKind) :: DTV(1:n,0:m,1:l)
  REAL(RealKind) :: DTW(1:n,1:m,0:l)
  REAL(RealKind) :: u(0:n,1:m,1:l)
  REAL(RealKind) :: v(1:n,0:m,1:l)
  REAL(RealKind) :: w(1:n,1:m,0:l)

  
  Div = SUM(-WeiFU(0,:,:)*DTU(0,:,:)*u(0,:,:)) + SUM(WeiFU(n,:,:)*DTU(n,:,:)*u(n,:,:)) & 
 &    + SUM(-WeiFV(:,0,:)*DTV(:,0,:)*v(:,0,:)) + SUM(WeiFV(:,m,:)*DTV(:,m,:)*v(:,m,:)) &  
 &    + SUM(-WeiFW(:,:,0)*DTW(:,:,0)*w(:,:,0)) + SUM(WeiFW(:,:,l)*DTW(:,:,l)*w(:,:,l)) 
 
END SUBROUTINE DivInter

SUBROUTINE Rhs(b,Vel,VecC,VecG)
  
  TYPE(PressureVelocity) :: b(:)
  TYPE(VelocityFace_T) :: Vel(:)
  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VecG(:)

  REAL(RealKind) :: y(nb)
  INTEGER :: info
  INTEGER :: ix,iy,iz
  INTEGER :: i


! -- Berechnung der rechten Seite --
!
! -- Blockelemente 1 bis N=nb --
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    CALL DivergenceS(b(ibLoc)%p(1,:,:,:),nx,ny,nz, &
                    FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1), &
                    FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1), &
                    FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1), &
                    DTUG(ibLoc)%uF,DTUG(ibLoc)%vF,DTUG(ibLoc)%wF,      &
                    Vel(ibLoc)%uF,Vel(ibLoc)%vF,Vel(ibLoc)%wF,DTTG(ibLoc)%p)
    IF (FacAnela>=1.0d0) THEN !OSSI
      b(ibLoc)%p(1,:,:,:)=beta0*dtP*b(ibLoc)%p(1,:,:,:) &
                        -VolB*VecC(ibLoc)%Vec(thPos)%cInt(:,:,:,1)    
    ELSE
      b(ibLoc)%p(1,:,:,:)=beta0*dtP*b(ibLoc)%p(1,:,:,:)
!                       -VolB*VecG(ibLoc)%Vec(1)%cInt(:,:,:,1)
    END IF  
    IF (MultiTriTB.OR.MultiTriTR.OR.MultiMuTR) THEN  
      CALL Divergence(b(ibLoc)%p(2,:,:,:),nx,ny,nz, &
                     FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1), &
                     FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1), &
                     FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1), &
                     Vel(ibLoc)%uF,Vel(ibLoc)%vF,Vel(ibLoc)%wF,DTTG(ibLoc)%p)
      b(ibLoc)%p(2,:,:,:)=beta0*dtP*b(ibLoc)%p(2,:,:,:) &
                      -VolB*VecC(ibLoc)%Vec(RhoPos)%cInt(:,:,:,1)    
    END IF
  END DO  
      
  IF (Anelastic.OR.PseudoIn) THEN
! -- Blockelement N+1 (Interface) = -d --
    ibLoc=0
    DO ib=1,nb
      IF (MyId == blMPI(ib)%Proc) THEN
        ibLoc=ibLoc+1
        CALL Set(Floor(ib))
        CALL DivInter(y(ib),nx,ny,nz,   &
                      FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1), &
                      FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1), &
                      FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1), &
                      DTUG(ibLoc)%uF,DTUG(ibLoc)%vF,DTUG(ibLoc)%wF,      &
                      Vel(ibLoc)%uF, Vel(ibLoc)%vF, Vel(ibLoc)%wF)
        y(ib)=beta0*dtP*y(ib)
        y(ib)=y(ib)*Boundary
      END IF
! -- Senden der Komponente ib an die anderen Prozesse --
     CALL MPI_Bcast(y(ib),1,MPI_RealKind,blMPI(ib)%Proc, &
 &                  MPI_COMM_WORLD,MPIErr)
    END DO 
  

! -- Loesen des Systems D * vgamma_s = d  <==>  vgamma_s = D*(D(+)*D)inv * d --
    CALL dposl(LaplNeumGlob,nb,nb,y) ! y = (D(+)*D)inv * d
    DO ib=1,nb
      y(ib)=y(ib)*Floor(ib)%Boundary
    END DO
    CALL UpdateRhsInterface(b,y)  ! b(N+1) = {b(N+1)} - Dtrans*y = -D(+)*(-d)
  
!   -- Addieren von C_j * vgamma_s zur rechten Seite (Komponenten 1 bis N) --

    CALL UpdateRhsPressure(b,y)                   ! b(j) = r_j - C_j * vgamma_s
  END IF  

END SUBROUTINE rhs

SUBROUTINE UpdateRhsInterface(b,y)  
  
  TYPE(PressureVelocity), TARGET :: b(:)
  REAL(RealKind) :: y(nb)
  
  TYPE(PressureVelocity), POINTER :: BVec

  TYPE(Nachbar_T), POINTER   :: Nachbar
  
  INTEGER :: in
  INTEGER :: jx,jy,jz
  REAL(RealKind) :: Temp
  REAL(RealKind) :: GradY
  REAL(RealKind), POINTER :: DUU(:,:,:)
  REAL(RealKind), POINTER :: DUV(:,:,:)
  REAL(RealKind), POINTER :: DUW(:,:,:)

!-- Aufdatieren von b = b - Dtrans * y auf dem Interface --
!--  (verteilte Vektoren) 
!-- alle lokalen Bloecke --
  
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    BVec => b(ibLoc)
    CALL Set(Floor(ib))
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      GradY=y(ibLoc)-y(ibn)
!  -- Westlicher Rand --
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (Nachbar%nType(1:1)/='o'.AND.RefineNachbar>=Refine) THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              BVec%u_w(jy,jz)=-beta0*dtP*GradY*FUG(ix0,jy,jz)*DUU(ix0,jy,jz)/(VolFace(ibLoc)%u_w(jy,jz)+Eps)* &
                              0.5e0*VolC(ix0+1,jy,jz)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%West==1) THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              BVec%u_w(jy,jz)=-beta0*dtP*y(ibLoc)*FUG(ix0,jy,jz)*DUU(ix0,jy,jz)/(VolFace(ibLoc)%u_w(jy,jz)+Eps)* &
                              0.5e0*VolC(ix0+1,jy,jz)
            END DO
          END DO
        ELSE 
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_w(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
               -GradY*SUM(FUG(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))* &
                0.5e0*VolC(ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1)
            END DO
          END DO
        END IF
      END IF
! --  Oestlicher Rand --     
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (Nachbar%nType(1:1)/='o'.AND.RefineNachbar>=Refine) THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              BVec%u_e(jy,jz)=beta0*dtP*GradY*FUG(ix1,jy,jz)*DUU(ix1,jy,jz)/(VolFace(ibLoc)%u_e(jy,jz)+Eps)* &
                              0.5d0*VolC(ix1,jy,jz)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%East==1) THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              BVec%u_e(jy,jz)=+beta0*dtP*y(ibLoc)*FUG(ix1,jy,jz)*DUU(ix1,jy,jz)/(VolFace(ibLoc)%u_e(jy,jz)+Eps)* &
                              0.5d0*VolC(ix1,jy,jz)
            END DO
          END DO
        ELSE
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              BVec%u_e(jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
               GradY*SUM(FUG(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))*   &
               0.5e0*VolC(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)
            END DO
          END DO
        END IF
      END IF
!  -- Suedlicher Rand --     
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (Nachbar%nType(1:1)/='o'.AND.RefineNachbar>=Refine) THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              BVec%v_s(jx,jz)=-beta0*dtP*GradY*FVG(jx,iy0,jz)*DUV(jx,iy0,jz)/(VolFace(ibLoc)%v_s(jx,jz)+Eps)* &
                              0.5d0*VolC(jx,iy0+1,jz)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%South==1) THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              BVec%v_s(jx,jz)=-beta0*dtP*y(ibLoc)*FVG(jx,iy0,jz)*DUV(jx,iy0,jz)/(VolFace(ibLoc)%v_s(jx,jz)+Eps)* &
                              0.5d0*VolC(jx,iy0+1,jz)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_s(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
               -GradY*SUM(FVG(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))* &
                0.5e0*VolC(jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)
            END DO
          END DO
        END IF
      END IF
!  -- Noerdlicher Rand --     
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (Nachbar%nType(1:1)/='o'.AND.RefineNachbar>=Refine) THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              BVec%v_n(jx,jz)=beta0*dtP*GradY*FVG(jx,iy1,jz)*DUV(jx,iy1,jz)/(VolFace(ibLoc)%v_n(jx,jz)+Eps)* &
                              0.5d0*VolC(jx,iy1,jz)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%North==1) THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              BVec%v_n(jx,jz)=+beta0*dtP*y(ibLoc)*FVG(jx,iy1,jz)*DUV(jx,iy1,jz)/(VolFace(ibLoc)%v_n(jx,jz)+Eps)* &
                              0.5d0*VolC(jx,iy1,jz)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              BVec%v_n(jx:jx+IncrX-1,jz:jz+IncrZ-1)= &
               GradY*SUM(FVG(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))* &
               0.5e0*VolC(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)
            END DO
          END DO
        END IF
      END IF
! --  Unterer Rand --     
      IF (Nachbar%nType(2:2) == 'b') THEN
        IF (Nachbar%nType(1:1)/='o'.AND.RefineNachbar>=Refine) THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              BVec%w_b(jx,jy)=-beta0*dtP*GradY*FWG(jx,jy,iz0)*DUW(jx,jy,iz0)/(VolFace(ibLoc)%w_b(jx,jy)+Eps)* &
                              0.5d0*VolC(jx,jy,iz0+1)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Bottom==1) THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              BVec%w_b(jx,jy)=-beta0*dtP*y(ibLoc)*FWG(jx,jy,iz0)*DUW(jx,jy,iz0)/(VolFace(ibLoc)%w_b(jx,jy)+Eps)* &
                              0.5d0*VolC(jx,jy,iz0+1)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_b(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
               -GradY*SUM(FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))* &
                0.5e0*VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)
            END DO
          END DO
        END IF
      END IF
! --  Oberer Rand --     
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (Nachbar%nType(1:1)/='o'.AND.RefineNachbar>=Refine) THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              BVec%w_t(jx,jy)=beta0*dtP*GradY*FWG(jx,jy,iz1)*DUW(jx,jy,iz1)/(VolFace(ibLoc)%w_t(jx,jy)+Eps)* &
                              0.5d0*VolC(jx,jy,iz1)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Top==1) THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1
              BVec%w_t(jx,jy)=+beta0*dtP*y(ibLoc)*FWG(jx,jy,iz1)*DUW(jx,jy,iz1)/(VolFace(ibLoc)%w_t(jx,jy)+Eps)* &
                              0.5d0*VolC(jx,jy,iz1)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              BVec%w_t(jx:jx+IncrX-1,jy:jy+IncrY-1)= &
               GradY*SUM(FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))* &
               0.5e0*VolC(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)
            END DO
          END DO
        END IF
      END IF
    END DO    ! in
  END DO      ! ibLoc
    
END SUBROUTINE UpdateRhsInterface  
  
SUBROUTINE UpdateRhsPressure(b,y)
  
  TYPE(PressureVelocity) :: b(:)
  REAL(RealKind) :: y(nb)
  
  INTEGER :: in
  INTEGER :: jx,jy,jz
  Real(8) :: Temp
  Real(8) :: GradY

  REAL(RealKind), POINTER :: bp(:,:,:,:)
  TYPE(Nachbar_T), POINTER :: Nachbar
  REAL(RealKind), POINTER :: DUU(:,:,:)
  REAL(RealKind), POINTER :: DUV(:,:,:)
  REAL(RealKind), POINTER :: DUW(:,:,:)
  REAL(RealKind), POINTER :: DTU(:,:,:)
  REAL(RealKind), POINTER :: DTV(:,:,:)
  REAL(RealKind), POINTER :: DTW(:,:,:)
  
! -- Aufdatieren von r_j = b(j)%p = b(j)%p - C_j * vgamma_s --

! -- alle lokalen Bloecke, alle Nachbarn --
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Call Set(Floor(ib))
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    DTU=>DTUG(ibLoc)%uF
    DTV=>DTUG(ibLoc)%vF
    DTW=>DTUG(ibLoc)%wF
    bp => b(ibLoc)%p
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      GradY=y(ibLoc)-y(ibn)
! --  Westlicher Rand --     
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              bp(1,ix0+1,jy,jz)=bp(1,ix0+1,jy,jz)-beta0*dtP*beta0*dtP*(y(ibLoc)-y(ibn))* &
                            FU(ix0,jy,jz)*DTU(ix0,jy,jz)*FUG(ix0,jy,jz)*DUU(ix0,jy,jz)/(VolFace(ibLoc)%u_w(jy,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%West==1) THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              bp(1,ix0+1,jy,jz)=bp(1,ix0+1,jy,jz)-beta0*dtP*beta0*dtP*y(ibLoc)* &
                            FU(ix0,jy,jz)*DTU(ix0,jy,jz)*FUG(ix0,jy,jz)*DUU(ix0,jy,jz)/(VolFace(ibLoc)%u_w(jy,jz)+Eps)
            END DO
          END DO
        ELSE
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              bp(1,ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              bp(1,ix0+1,jy:jy+IncrY-1,jz:jz+IncrZ-1)-beta0*dtP* &
              Grady*SUM(FU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DTU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1))* &
              FUG(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)          
            END DO
          END DO
        END IF
      END IF
! -- Oestlicher Rand --     
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              bp(1,ix1,jy,jz)=bp(1,ix1,jy,jz)-beta0*dtP*beta0*dtP*(y(ibLoc)-y(ibn))* &
                            FU(ix1,jy,jz)*DTU(ix1,jy,jz)*FUG(ix1,jy,jz)*DUU(ix1,jy,jz)/(VolFace(ibLoc)%u_e(jy,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%East==1) THEN
          DO jy=jy0+1,jy1
            DO jz=jz0+1,jz1
              bp(1,ix1,jy,jz)=bp(1,ix1,jy,jz)-beta0*dtP*beta0*dtP*y(ibLoc)* &
                            FU(ix1,jy,jz)*DTU(ix1,jy,jz)*FUG(ix1,jy,jz)*DUU(ix1,jy,jz)/(VolFace(ibLoc)%u_e(jy,jz)+Eps)
            END DO
          END DO
        ELSE
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              bp(1,ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
              bp(1,ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)-beta0*dtP* &
              Grady*SUM(FU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DTU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1))* &
              FUG(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)*DUU(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)          
            END DO
          END DO
        END IF
      END IF
           
! -- Suedlicher Rand --     
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              bp(1,jx,iy0+1,jz)=bp(1,jx,iy0+1,jz)-beta0*dtP*beta0*dtP*(y(ibLoc)-y(ibn))* &
                              FV(jx,iy0,jz)*DTV(jx,iy0,jz)*FVG(jx,iy0,jz)*DUV(jx,iy0,jz)/(VolFace(ibLoc)%v_s(jx,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%South==1) THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              bp(1,jx,iy0+1,jz)=bp(1,jx,iy0+1,jz)-beta0*dtP*beta0*dtP*y(ibLoc)* &
                              FV(jx,iy0,jz)*DTV(jx,iy0,jz)*FVG(jx,iy0,jz)*DUV(jx,iy0,jz)/(VolFace(ibLoc)%v_s(jx,jz)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              bp(1,jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)= &
              bp(1,jx:jx+IncrX-1,iy0+1,jz:jz+IncrZ-1)-beta0*dtP* &
              GradY*SUM(FV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*DTV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1))* &
              FVG(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)          
            END DO
          END DO
        END IF
      END IF
           
! -- Noerdlicher Rand --     
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              bp(1,jx,iy1,jz)=bp(1,jx,iy1,jz)-beta0*dtP*beta0*dtP*(y(ibLoc)-y(ibn))* &
                            FV(jx,iy1,jz)*DTV(jx,iy1,jz)*FVG(jx,iy1,jz)*DUV(jx,iy1,jz)/(VolFace(ibLoc)%v_n(jx,jz)+Eps)
            END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%North==1) THEN
          DO jx=jx0+1,jx1
            DO jz=jz0+1,jz1
              bp(1,jx,iy1,jz)=bp(1,jx,iy1,jz)-beta0*dtP*beta0*dtP*y(ibLoc)* &
                            FV(jx,iy1,jz)*DTV(jx,iy1,jz)*FVG(jx,iy1,jz)*DUV(jx,iy1,jz)/(VolFace(ibLoc)%v_n(jx,jz)+Eps)
            END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              bp(1,jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
              bp(1,jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)-beta0*dtP* &
              GradY*SUM(FV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*DTV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1))* &
              FVG(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)*DUV(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)          
            END DO
          END DO
        END IF
      END IF
           
! -- Unterer Rand --     
      IF (Nachbar%nType(2:2) == 'b') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1 
               bp(1,jx,jy,iz0+1)=bp(1,jx,jy,iz0+1) -beta0*dtP*beta0*dtP* (y(ibLoc)-y(ibn)) * &
                                 FW(jx,jy,iz0)*DTW(jx,jy,iz0)*FWG(jx,jy,iz0)*DUW(jx,jy,iz0)/(VolFace(ibLoc)%w_b(jx,jy)+Eps)
             END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Bottom==1) THEN
          DO jx=jx0+1,jx1
            DO jy=jy0+1,jy1 
               bp(1,jx,jy,iz0+1)=bp(1,jx,jy,iz0+1) -beta0*dtP*beta0*dtP*y(ibLoc)* &
                                 FW(jx,jy,iz0)*DTW(jx,jy,iz0)*FWG(jx,jy,iz0)*DUW(jx,jy,iz0)/(VolFace(ibLoc)%w_b(jx,jy)+Eps)
             END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              bp(1,jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)= &
              bp(1,jx:jx+IncrX-1,jy:jy+IncrY-1,iz0+1)-beta0*dtP* &
              GradY*SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*DTW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0))* &
              FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)          
            END DO
          END DO
        END IF
      END IF
           
! -- Oberer Rand --     
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar>=Refine.AND.Nachbar%nType(1:1)/='o') THEN
          DO jx=jx0+1,jx1
             DO jy=jy0+1,jy1
                bp(1,jx,jy,iz1)=bp(1,jx,jy,iz1)-beta0*dtP*beta0*dtP*(y(ibLoc)-y(ibn))* &
                                FW(jx,jy,iz1)*DTW(jx,jy,iz1)*FWG(jx,jy,iz1)*DUW(jx,jy,iz1)/(VolFace(ibLoc)%w_t(jx,jy)+Eps)
             END DO
          END DO
        ELSE IF (Nachbar%nType(1:1)=='o'.AND.BCP%Top==1) THEN
          DO jx=jx0+1,jx1
             DO jy=jy0+1,jy1
                bp(1,jx,jy,iz1)=bp(1,jx,jy,iz1)-beta0*dtP*beta0*dtP*y(ibLoc)* &
                                FW(jx,jy,iz1)*DTW(jx,jy,iz1)*FWG(jx,jy,iz1)*DUW(jx,jy,iz1)/(VolFace(ibLoc)%w_t(jx,jy)+Eps)
             END DO
          END DO
        ELSE
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              bp(1,jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
              bp(1,jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)-beta0*dtP* &
              GradY*SUM(FW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*DTW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1))* &
              FWG(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)*DUW(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)          
            END DO
          END DO
        END IF
      END IF
    END DO    ! in
  END DO       ! ibLoc

END SUBROUTINE UpdateRhsPressure

SUBROUTINE UpdateVelocity(x,b,Vel,VecC,VelC,VecG)
  
  TYPE(PressureVelocity), TARGET :: x(:),b(:)
  TYPE(VelocityFace_T) :: Vel(:)
  TYPE(Vector4Cell_T) :: VecC(:)
  TYPE(Vector4Cell_T) :: VelC(:)
  TYPE(Vector4Cell_T), OPTIONAL :: VecG(:)

  INTEGER :: ix,iy,iz,ic,it
  REAL(RealKind) :: Temp 
  TYPE(PressureVelocity), POINTER :: XVec, BVec
  REAL(RealKind), POINTER :: DUU(:,:,:)
  REAL(RealKind), POINTER :: DUV(:,:,:)
  REAL(RealKind), POINTER :: DUW(:,:,:)
  REAL(RealKind), POINTER :: DTT(:,:,:)

! -- Aufdatieren von x%U = x_0%U + x%V (verteilter Vektor!) --


! -- Aufdatieren --
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)

    XVec => x(ibLoc)
    BVec => b(ibLoc)
    DUU=>DUUG(ibLoc)%uF
    DUV=>DUUG(ibLoc)%vF
    DUW=>DUUG(ibLoc)%wF
    DTT=>DTTG(ibLoc)%p
    CALL Set(Floor(ib))
    XVec%p(1,:,:,:)=XVec%p(1,:,:,:)*DUTG(ibLoc)%p
    IF (MultiTriTR.OR.MultiTriTB.OR.MultiMuTR) THEN 
      DO iz=iz0+1,iz1-1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1                       
            Temp=beta0*dtp*GravComp*(VolC(ix,iy,iz)*XVec%p(2,ix,iy,iz) &
                +VolC(ix,iy,iz+1)*XVec%p(2,ix,iy,iz+1)) &
                /(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)+Eps)
            Vel(ibLoc)%wF(ix,iy,iz)=Vel(ibLoc)%wF(ix,iy,iz)-Temp
          END DO
        END DO
      END DO
      XVec%p(2,:,:,:)=XVec%p(2,:,:,:)*DURG(ibLoc)%p
    END IF
     
    IF (GradFull) THEN 
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          Temp=XVec%u_w(iy,iz)-BVec%u_w(iy,iz)/(0.5e0*dx(ix0+1)*MetrXY(iy)*MetrXZ(iz))
          Vel(ibLoc)%uF(ix0,iy,iz)=Vel(ibLoc)%uF(ix0,iy,iz)+Temp
          Temp=XVec%u_e(iy,iz)-BVec%u_e(iy,iz)/(0.5e0*dx(ix1)*MetrXY(iy)*MetrXZ(iz))
          Vel(ibLoc)%uF(ix1,iy,iz)=Vel(ibLoc)%uF(ix1,iy,iz)+Temp
        END DO
      END DO
      DO ix=ix0+1,ix1
        DO iz=iz0+1,iz1
          Temp=XVec%v_s(ix,iz)-BVec%v_s(ix,iz)/(0.5e0*dy(iy0+1)*MetrYX(ix)*MetrYZ(iz))
          Vel(ibLoc)%vF(ix,iy0,iz)=Vel(ibLoc)%vF(ix,iy0,iz)+Temp
          Temp=XVec%v_n(ix,iz)-BVec%v_n(ix,iz)/(0.5e0*dy(iy1)*MetrYX(ix)*MetrYZ(iz))
          Vel(ibLoc)%vF(ix,iy1,iz)=Vel(ibLoc)%vF(ix,iy1,iz)+Temp
        END DO
      END DO
      DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          Temp=XVec%w_b(ix,iy)-BVec%w_b(ix,iy)/(0.5e0*dz(iz0+1)*MetrZX(ix)*MetrZY(iy))
          Vel(ibLoc)%wF(ix,iy,iz0)=Vel(ibLoc)%wF(ix,iy,iz0)+Temp
          Temp=XVec%w_t(ix,iy)-BVec%w_t(ix,iy)/(0.5e0*dz(iz1)*MetrZX(ix)*MetrZY(iy))
          Vel(ibLoc)%wF(ix,iy,iz1)=Vel(ibLoc)%wF(ix,iy,iz1)+FW(ix,iy,iz1)/(FW(ix,iy,iz1)+Eps)*Temp
        END DO
      END DO
! -- Geschwindigkeiten im Innern des Blocks: u = utilde - grad p --
      IF (MultiTriT.OR.MultiMu.OR.MultiEx) THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1-1                       
              Temp=Two*beta0*dtp*DUU(ix,iy,iz)*(XVec%p(1,ix+1,iy,iz)-XVec%p(1,ix,iy,iz)) &
                /((dx(ix+1)+dx(ix))*MetrXY(iy)*MetrXZ(iz))
              Vel(ibLoc)%uF(ix,iy,iz)=Vel(ibLoc)%uF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1-1
            DO ix=ix0+1,ix1                       
              Temp=Two*beta0*dtp*DUV(ix,iy,iz)*(XVec%p(1,ix,iy+1,iz)-XVec%p(1,ix,iy,iz)) &
                /((dy(iy+1)+dy(iy))*MetrYX(ix)*MetrYZ(iz))
              Vel(ibLoc)%vF(ix,iy,iz)=Vel(ibLoc)%vF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
        DO iz=iz0+1,iz1-1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1                       
              Temp=Two*beta0*dtp*DUW(ix,iy,iz)*(XVec%p(1,ix,iy,iz+1)-XVec%p(1,ix,iy,iz)) &
                /(dz(iz+1)+dz(iz))
              Vel(ibLoc)%wF(ix,iy,iz)=Vel(ibLoc)%wF(ix,iy,iz)-FW(ix,iy,iz)/(FW(ix,iy,iz)+Eps)*Temp
            END DO
          END DO
        END DO
      ELSE
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1-1
              Temp=Two*beta0*dtp*DUU(ix,iy,iz) &
                *((XVec%p(1,ix+1,iy,iz)-XVec%p(1,ix,iy,iz)) &
                 +(XVec%p(2,ix+1,iy,iz)-XVec%p(2,ix,iy,iz))) &
                /(dx(ix+1)+dx(ix))
              Vel(ibLoc)%uF(ix,iy,iz)=Vel(ibLoc)%uF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1-1
            DO ix=ix0+1,ix1
              Temp=Two*beta0*dtp*DUV(ix,iy,iz) &
                  *((XVec%p(1,ix,iy+1,iz)-XVec%p(1,ix,iy,iz)) &
                   +(XVec%p(2,ix,iy+1,iz)-XVec%p(2,ix,iy,iz))) &
                /(dy(iy+1)+dy(iy))
              Vel(ibLoc)%vF(ix,iy,iz)=Vel(ibLoc)%vF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
        DO iz=iz0+1,iz1-1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1
              Temp=Two*beta0*dtp*DUW(ix,iy,iz) &
                  *((XVec%p(1,ix,iy,iz+1)-XVec%p(1,ix,iy,iz)) &
                   +(XVec%p(2,ix,iy,iz+1)-XVec%p(2,ix,iy,iz))) &
                /(dz(iz+1)+dz(iz))
              Vel(ibLoc)%wF(ix,iy,iz)=Vel(ibLoc)%wF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
      END IF
    ELSE  
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          Temp=XVec%u_w(iy,iz)-BVec%u_w(iy,iz)/(0.5e0*VolC(ix0+1,iy,iz)+Eps)
          Vel(ibLoc)%uF(ix0,iy,iz)=Vel(ibLoc)%uF(ix0,iy,iz)+Temp
          Temp=XVec%u_e(iy,iz)-BVec%u_e(iy,iz)/(0.5e0*VolC(ix1,iy,iz)+Eps)
          Vel(ibLoc)%uF(ix1,iy,iz)=Vel(ibLoc)%uF(ix1,iy,iz)+Temp
        END DO
      END DO
      DO ix=ix0+1,ix1
        DO iz=iz0+1,iz1
          Temp=XVec%v_s(ix,iz)-BVec%v_s(ix,iz)/(0.5e0*VolC(ix,iy0+1,iz)+Eps)
          Vel(ibLoc)%vF(ix,iy0,iz)=Vel(ibLoc)%vF(ix,iy0,iz)+Temp
          Temp=XVec%v_n(ix,iz)-BVec%v_n(ix,iz)/(0.5e0*VolC(ix,iy1,iz)+Eps)
          Vel(ibLoc)%vF(ix,iy1,iz)=Vel(ibLoc)%vF(ix,iy1,iz)+Temp
        END DO
      END DO
      DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          Temp=XVec%w_b(ix,iy)-BVec%w_b(ix,iy)/(0.5e0*VolC(ix,iy,iz0+1)+Eps)
          Vel(ibLoc)%wF(ix,iy,iz0)=Vel(ibLoc)%wF(ix,iy,iz0)+Temp
          Temp=XVec%w_t(ix,iy)-BVec%w_t(ix,iy)/(0.5e0*VolC(ix,iy,iz1)+Eps)
          Vel(ibLoc)%wF(ix,iy,iz1)=Vel(ibLoc)%wF(ix,iy,iz1)+Temp
        END DO
      END DO
      IF (MultiTriT.OR.MultiMu.OR.MultiEx) THEN 
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1-1                       
              Temp=Two*beta0*dtp*FU(ix,iy,iz)*DUU(ix,iy,iz)*(XVec%p(1,ix+1,iy,iz)-XVec%p(1,ix,iy,iz)) &
                /(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)+Eps)
              Vel(ibLoc)%uF(ix,iy,iz)=Vel(ibLoc)%uF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1-1
            DO ix=ix0+1,ix1                       
              Temp=Two*beta0*dtp*FV(ix,iy,iz)*DUV(ix,iy,iz)*(XVec%p(1,ix,iy+1,iz)-XVec%p(1,ix,iy,iz)) &
                /(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)+Eps)
              Vel(ibLoc)%vF(ix,iy,iz)=Vel(ibLoc)%vF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
        DO iz=iz0+1,iz1-1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1                       
              Temp=Two*beta0*dtp*FW(ix,iy,iz)*DUW(ix,iy,iz)*(XVec%p(1,ix,iy,iz+1)-XVec%p(1,ix,iy,iz)) &
                /(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)+Eps)
              Vel(ibLoc)%wF(ix,iy,iz)=Vel(ibLoc)%wF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
      ELSE
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1-1
              Temp=Two*beta0*dtp*FU(ix,iy,iz)*DUU(ix,iy,iz) &
                *((XVec%p(1,ix+1,iy,iz)-XVec%p(1,ix,iy,iz)) &
                 +(XVec%p(2,ix+1,iy,iz)-XVec%p(2,ix,iy,iz))) &
                /(VolC(ix+1,iy,iz)+VolC(ix,iy,iz)+Eps)
              Vel(ibLoc)%uF(ix,iy,iz)=Vel(ibLoc)%uF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
        DO iz=iz0+1,iz1
          DO iy=iy0+1,iy1-1
            DO ix=ix0+1,ix1
              Temp=Two*beta0*dtp*FV(ix,iy,iz)*DUV(ix,iy,iz) &
                  *((XVec%p(1,ix,iy+1,iz)-XVec%p(1,ix,iy,iz)) &
                   +(XVec%p(2,ix,iy+1,iz)-XVec%p(2,ix,iy,iz))) &
                /(VolC(ix,iy+1,iz)+VolC(ix,iy,iz)+Eps)
              Vel(ibLoc)%vF(ix,iy,iz)=Vel(ibLoc)%vF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
        DO iz=iz0+1,iz1-1
          DO iy=iy0+1,iy1
            DO ix=ix0+1,ix1
              Temp=Two*beta0*dtp*FW(ix,iy,iz)*DUW(ix,iy,iz) &
                  *((XVec%p(1,ix,iy,iz+1)-XVec%p(1,ix,iy,iz)) &
                   +(XVec%p(2,ix,iy,iz+1)-XVec%p(2,ix,iy,iz))) &
                /(VolC(ix,iy,iz+1)+VolC(ix,iy,iz)+Eps)
              Vel(ibLoc)%wF(ix,iy,iz)=Vel(ibLoc)%wF(ix,iy,iz)-Temp
            END DO
          END DO
        END DO
      END IF
    END IF
    IF (ThPos>0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1                       
            VecC(ibLoc)%Vec(thPos)%c(ix,iy,iz,1)=VecC(ibLoc)%Vec(thPos)%c(ix,iy,iz,1) &
            -dtP*beta0*DTT(ix,iy,iz) &
            *(FU(ix,iy,iz)*DTUG(ibLoc)%uF(ix,iy,iz)*Vel(ibLoc)%uF(ix,iy,iz) &
             -FU(ix-1,iy,iz)*DTUG(ibLoc)%uF(ix-1,iy,iz)*Vel(ibLoc)%uF(ix-1,iy,iz) &
             +FV(ix,iy,iz)*DTUG(ibLoc)%vF(ix,iy,iz)*Vel(ibLoc)%vF(ix,iy,iz) &
             -FV(ix,iy-1,iz)*DTUG(ibLoc)%vF(ix,iy-1,iz)*Vel(ibLoc)%vF(ix,iy-1,iz) &
             +FW(ix,iy,iz)*DTUG(ibLoc)%wF(ix,iy,iz)*Vel(ibLoc)%wF(ix,iy,iz) &
             -FW(ix,iy,iz-1)*DTUG(ibLoc)%wF(ix,iy,iz-1)*Vel(ibLoc)%wF(ix,iy,iz-1)) & 
             /(VolC(ix,iy,iz)+Eps)
          END DO
        END DO
      END DO
    END IF
    IF (RhoPos>0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1                       
            VecC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1)=VecC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1) &
            -dtP*beta0 & 
            *(FU(ix,iy,iz)*Vel(ibLoc)%uF(ix,iy,iz) &
             -FU(ix-1,iy,iz)*Vel(ibLoc)%uF(ix-1,iy,iz) &
             +FV(ix,iy,iz)*Vel(ibLoc)%vF(ix,iy,iz) &
             -FV(ix,iy-1,iz)*Vel(ibLoc)%vF(ix,iy-1,iz) &
             +FW(ix,iy,iz)*Vel(ibLoc)%wF(ix,iy,iz) &
             -FW(ix,iy,iz-1)*Vel(ibLoc)%wF(ix,iy,iz-1)) &
             /(VolC(ix,iy,iz)+Eps)
          END DO
        END DO
      END DO
    END IF
    IF (EnPos>0) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1                       
            VecC(ibLoc)%Vec(enPos)%c(ix,iy,iz,1)=VecC(ibLoc)%Vec(enPos)%c(ix,iy,iz,1) &
            -dtP*beta0 &
            *(FU(ix,iy,iz)*Vel(ibLoc)%uF(ix,iy,iz) &
             *PreFace(VelC(ibLoc)%Vec(thPos)%c(ix,iy,iz,1),VelC(ibLoc)%Vec(thPos)%c(ix+1,iy,iz,1) &
                     ,VelC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1),VelC(ibLoc)%Vec(RhoPos)%c(ix+1,iy,iz,1) &
                     ,VolC(ix,iy,iz),VolC(ix+1,iy,iz)) &
             -FU(ix-1,iy,iz)*Vel(ibLoc)%uF(ix-1,iy,iz) &
             *PreFace(VelC(ibLoc)%Vec(thPos)%c(ix-1,iy,iz,1),VelC(ibLoc)%Vec(thPos)%c(ix,iy,iz,1) &
                     ,VelC(ibLoc)%Vec(RhoPos)%c(ix-1,iy,iz,1),VelC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1) &
                     ,VolC(ix-1,iy,iz),VolC(ix,iy,iz)) &
             +FV(ix,iy,iz)*Vel(ibLoc)%vF(ix,iy,iz) &
             *PreFace(VelC(ibLoc)%Vec(thPos)%c(ix,iy,iz,1),VelC(ibLoc)%Vec(thPos)%c(ix,iy+1,iz,1) &
                     ,VelC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1),VelC(ibLoc)%Vec(RhoPos)%c(ix,iy+1,iz,1) &
                     ,VolC(ix,iy,iz),VolC(ix,iy+1,iz)) &
             -FV(ix,iy-1,iz)*Vel(ibLoc)%vF(ix,iy-1,iz) &
             *PreFace(VelC(ibLoc)%Vec(thPos)%c(ix,iy-1,iz,1),VelC(ibLoc)%Vec(thPos)%c(ix,iy,iz,1) &
                     ,VelC(ibLoc)%Vec(RhoPos)%c(ix,iy-1,iz,1),VelC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1) &
                     ,VolC(ix,iy-1,iz),VolC(ix,iy,iz)) &
             +FW(ix,iy,iz)*Vel(ibLoc)%wF(ix,iy,iz) &
             *PreFace(VelC(ibLoc)%Vec(thPos)%c(ix,iy,iz,1),VelC(ibLoc)%Vec(thPos)%c(ix,iy,iz+1,1) &
                     ,VelC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1),VelC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz+1,1) &
                     ,VolC(ix,iy,iz),VolC(ix,iy,iz+1)) &
             -FW(ix,iy,iz-1)*Vel(ibLoc)%wF(ix,iy,iz-1) &
             *PreFace(VelC(ibLoc)%Vec(thPos)%c(ix,iy,iz-1,1),VelC(ibLoc)%Vec(thPos)%c(ix,iy,iz,1) &
                     ,VelC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz-1,1),VelC(ibLoc)%Vec(RhoPos)%c(ix,iy,iz,1) &
                     ,VolC(ix,iy,iz-1),VolC(ix,iy,iz))) &
             /(VolC(ix,iy,iz)+Eps)
          END DO
        END DO
      END DO
    END IF
    DO ic=1,uPosL-1
      IF (ic/=thPos.AND.ic/=RhoPos) THEN
        DO it=1,SIZE(VelC(ibLoc)%Vec(iC)%c,4)
          CALL DivScalar(VecC(ibLoc)%Vec(iC)%c(:,:,:,it),VelC(ibLoc)%Vec(iC)%c(:,:,:,it),VelC(ibLoc)%Vec(RhoPos)%c)
        END DO  
      END IF
    END DO  
    IF (PRESENT(VecG)) THEN
      VecG(ibLoc)%Vec(1)%c(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1,1)=Xvec%p(1,:,:,:)
    END IF  
  END DO
CONTAINS

FUNCTION PreFace(pL,pR,RhoL,RhoR,VolL,VolR)

  REAL(RealKind) :: PreFace
  REAL(RealKind) :: pL,pR,RhoL,RhoR,VolL,VolR

  PreFace=(pL/(RhoL+Eps)*RhoR+pR/(RhoR+Eps)*RhoL)/(RhoL+RhoR+Eps)
! PreFace=(pL*RhoR+pR*RhoL)/(RhoL+RhoR)

END FUNCTION PreFace  
SUBROUTINE DivScalar(cDiv,c,RhoC)
  REAL(RealKind) :: cDiv(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1)
  REAL(RealKind) :: c(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1)
  REAL(RealKind) :: RhoC(ix0:ix1+1,iy0:iy1+1,iz0:iz1+1,1)

  DO iz=iz0+1,iz1
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1                       
        cDiv(ix,iy,iz,1)=cDiv(ix,iy,iz,1) &
        -dtP*beta0 &
        *(FU(ix,iy,iz)*Vel(ibLoc)%uF(ix,iy,iz) &
         *(VolC(ix,iy,iz)*c(ix,iy,iz,1)/(RhoC(ix,iy,iz,1)+Eps) &
          +VolC(ix+1,iy,iz)*c(ix+1,iy,iz,1)/(RhoC(ix+1,iy,iz,1)+Eps)) &
          /(VolC(ix,iy,iz)+VolC(ix+1,iy,iz)+Eps) &
         -FU(ix-1,iy,iz)*Vel(ibLoc)%uF(ix-1,iy,iz) &
         *(VolC(ix,iy,iz)*c(ix,iy,iz,1)/(RhoC(ix,iy,iz,1)+Eps) &
          +VolC(ix-1,iy,iz)*c(ix-1,iy,iz,1)/(RhoC(ix-1,iy,iz,1)+Eps)) &
         /(VolC(ix,iy,iz)+VolC(ix-1,iy,iz)+Eps) &
         +FV(ix,iy,iz)*Vel(ibLoc)%vF(ix,iy,iz) &
         *(VolC(ix,iy,iz)*c(ix,iy,iz,1)/(RhoC(ix,iy,iz,1)+Eps) &
          +VolC(ix,iy+1,iz)*c(ix,iy+1,iz,1)/(RhoC(ix,iy+1,iz,1)+Eps)) &
         /(VolC(ix,iy,iz)+VolC(ix,iy+1,iz)+Eps) &
         -FV(ix,iy-1,iz)*Vel(ibLoc)%vF(ix,iy-1,iz) &
         *(VolC(ix,iy,iz)*c(ix,iy,iz,1)/(RhoC(ix,iy,iz,1)+Eps) &
          +VolC(ix,iy-1,iz)*c(ix,iy-1,iz,1)/(RhoC(ix,iy-1,iz,1)+Eps)) &
         /(VolC(ix,iy,iz)+VolC(ix,iy-1,iz)+Eps) &
         +FW(ix,iy,iz)*Vel(ibLoc)%wF(ix,iy,iz) &
         *(VolC(ix,iy,iz)*c(ix,iy,iz,1)/(RhoC(ix,iy,iz,1)+Eps) &
          +VolC(ix,iy,iz+1)*c(ix,iy,iz+1,1)/(RhoC(ix,iy,iz+1,1)+Eps)) &
         /(VolC(ix,iy,iz)+VolC(ix,iy,iz+1)+Eps) &
         -FW(ix,iy,iz-1)*Vel(ibLoc)%wF(ix,iy,iz-1) &
         *(VolC(ix,iy,iz)*c(ix,iy,iz,1)/(RhoC(ix,iy,iz,1)+Eps) &
          +VolC(ix,iy,iz-1)*c(ix,iy,iz-1,1)/(RhoC(ix,iy,iz-1,1)+Eps)) &
         /(VolC(ix,iy,iz)+VolC(ix,iy,iz-1)+Eps)) &
         /(VolC(ix,iy,iz)+Eps)
      END DO
    END DO
  END DO
END SUBROUTINE DivScalar
      
END SUBROUTINE UpdateVelocity
  
END MODULE Rhs_Mod
