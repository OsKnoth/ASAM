MODULE VecPrVel_Mod

  USE Kind_Mod
  USE Floor_Mod
  USE Parameter_Mod
  USE Parallel_Mod
  USE Exchange_Mod

  IMPLICIT NONE

  TYPE PressureVelocity
    INTEGER :: nx,ny,nz
    INTEGER :: refine
    REAL(RealKind), POINTER :: p(:,:,:,:)
    REAL(RealKind), POINTER :: u_e(:,:),u_w(:,:)
    REAL(RealKind), POINTER :: v_n(:,:),v_s(:,:)
    REAL(RealKind), POINTER :: w_t(:,:),w_b(:,:)
  END TYPE PressureVelocity

  TYPE Pressure_T
    REAL(RealKind), POINTER :: p(:,:,:)
  END TYPE Pressure_T

  INTEGER, PRIVATE :: np1=1
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE DDOT_PressureVelocity
  END INTERFACE
  INTERFACE DOT2
    MODULE PROCEDURE DDOT2_PressureVelocity
  END INTERFACE
  INTERFACE DOTR2
    MODULE PROCEDURE DDOTR2_PressureVelocity
  END INTERFACE
  INTERFACE DOTI2
    MODULE PROCEDURE DDOTI2_PressureVelocity
  END INTERFACE
  INTERFACE NormV
    MODULE PROCEDURE NormV_PressureVelocity
  END INTERFACE
  INTERFACE DDOTR
    MODULE PROCEDURE DDOTR_PressureVelocity
  END INTERFACE
  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE COPY_PressureVelocity, &
                     Value_PressureVelocity
  END INTERFACE
  INTERFACE DSCALE
    MODULE PROCEDURE SCALE_PressureVelocity, &
                     SCALE2_PressureVelocity
  END INTERFACE
  INTERFACE DXPY
    MODULE PROCEDURE DXPY_PressureVelocity
  END INTERFACE
  INTERFACE DAXPY
    MODULE PROCEDURE DAXPY_PressureVelocity
  END INTERFACE
  INTERFACE DXPBY
    MODULE PROCEDURE DXPBY_PressureVelocity
  END INTERFACE
  INTERFACE DAXPBY
    MODULE PROCEDURE DAXPBY_PressureVelocity
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_PressureVelocity, &
                     Allocate_Pressure, &
                     Allocate2_PressureVelocity
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE Deallocate_PressureVelocity, &
                     Deallocate2_PressureVelocity
  END INTERFACE
  INTERFACE Output
    MODULE PROCEDURE Output_PressureVelocity, &
                     Output2_PressureVelocity
  END INTERFACE

CONTAINS

SUBROUTINE SetDim1P(n)
  INTEGER :: n
  np1=n
END SUBROUTINE SetDim1P

FUNCTION DDOT_PressureVelocity(xV,yV)

  REAL(RealKind) :: DDOT_PressureVelocity
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:),yV(:)

  REAL(RealKind) :: SUM,Temp
  INTEGER :: Ierr
  INTEGER :: nx,ny,nz


  Temp=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Temp = Temp + SUM(xV(ibLoc)%p   * yV(ibLoc)%p)
    Temp = Temp + SUM(xV(ibLoc)%u_e * yV(ibLoc)%u_e)
    Temp = Temp + SUM(xV(ibLoc)%u_w * yV(ibLoc)%u_w)
    Temp = Temp + SUM(xV(ibLoc)%v_n * yV(ibLoc)%v_n)
    Temp = Temp + SUM(xV(ibLoc)%v_s * yV(ibLoc)%v_s)
    Temp = Temp + SUM(xV(ibLoc)%w_t * yV(ibLoc)%w_t)
    Temp = Temp + SUM(xV(ibLoc)%w_b * yV(ibLoc)%w_b)
  END DO

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,DDOT_PressureVelocity,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,Ierr)

END FUNCTION DDOT_PressureVelocity

FUNCTION DDOT2_PressureVelocity(xV,yV)

  REAL(RealKind) :: DDOT2_PressureVelocity
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:),yV(:)

  REAL(RealKind) :: Temp
  INTEGER :: in, Ierr
  INTEGER :: nx,ny,nz

  Temp=Zero
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc) 
    CALL Set(Floor(ib))
    Temp = Temp + SUM(xV(ibLoc)%p   * yV(ibLoc)%p)
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      IF (Nachbar%nType(2:2)=='w') THEN
        IF (Refine>RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%u_w(jy0+1:jy1,jz0+1:jz1)* &
                         yV(ibLoc)%u_w(jy0+1:jy1,jz0+1:jz1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='e') THEN
        IF (Refine>=RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%u_e(jy0+1:jy1,jz0+1:jz1)* &
                         yV(ibLoc)%u_e(jy0+1:jy1,jz0+1:jz1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='s') THEN
        IF (Refine>RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%v_s(jx0+1:jx1,jz0+1:jz1)* &
                         yV(ibLoc)%v_s(jx0+1:jx1,jz0+1:jz1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='n') THEN
        IF (Refine>=RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%v_n(jx0+1:jx1,jz0+1:jz1)* &
                         yV(ibLoc)%v_n(jx0+1:jx1,jz0+1:jz1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='b') THEN
        IF (Refine>RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%w_b(jx0+1:jx1,jy0+1:jy1)* &
                         yV(ibLoc)%w_b(jx0+1:jx1,jy0+1:jy1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='t') THEN
        IF (Refine>=RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%w_t(jx0+1:jx1,jy0+1:jy1)* &
                         yV(ibLoc)%w_t(jx0+1:jx1,jy0+1:jy1))
        END IF
      END IF
    END DO
  END DO

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,DDOT2_PressureVelocity,1,MPI_RealKind, &
                     MPI_SUM,MPI_Comm_World,Ierr)

END FUNCTION DDOT2_PressureVelocity

FUNCTION DDOTR2_PressureVelocity(xV,yV)

  REAL(RealKind) :: DDOTR2_PressureVelocity
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:),yV(:)

  REAL(RealKind) :: SUM,Temp
  INTEGER :: in, Ierr
  INTEGER :: nx,ny,nz

  Temp=0.0e0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)
      IF (Nachbar%nType(2:2)=='w') THEN
        IF (Refine>RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%u_w(jy0+1:jy1,jz0+1:jz1)* &
                         yV(ibLoc)%u_w(jy0+1:jy1,jz0+1:jz1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='e') THEN
        IF (Refine>=RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%u_e(jy0+1:jy1,jz0+1:jz1)* &
                         yV(ibLoc)%u_e(jy0+1:jy1,jz0+1:jz1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='s') THEN
        IF (Refine>RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%v_s(jx0+1:jx1,jz0+1:jz1)* &
                         yV(ibLoc)%v_s(jx0+1:jx1,jz0+1:jz1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='n') THEN
        IF (Refine>=RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%v_n(jx0+1:jx1,jz0+1:jz1)* &
                         yV(ibLoc)%v_n(jx0+1:jx1,jz0+1:jz1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='b') THEN
        IF (Refine>RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%w_b(jx0+1:jx1,jy0+1:jy1)* &
                         yV(ibLoc)%w_b(jx0+1:jx1,jy0+1:jy1))
        END IF
      END IF
      IF (Nachbar%nType(2:2)=='t') THEN
        IF (Refine>=RefineNachbar.OR.Nachbar%nType(1:1)=='o') THEN
           Temp=Temp+SUM(xV(ibLoc)%w_t(jx0+1:jx1,jy0+1:jy1)* &
                         yV(ibLoc)%w_t(jx0+1:jx1,jy0+1:jy1))
        END IF
      END IF
    END DO
  END DO

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,DDOTR2_PressureVelocity,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,Ierr)

END FUNCTION DDOTR2_PressureVelocity

FUNCTION DDOTI2_PressureVelocity(xV,yV)

  REAL(RealKind) :: DDOTI2_PressureVelocity
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:),yV(:)

  REAL(RealKind) :: SUM,Temp
  INTEGER :: in, Ierr
  INTEGER :: nx,ny,nz

  Temp=0.0e0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Temp = Temp + SUM(xV(ibLoc)%p   * yV(ibLoc)%p)
  END DO

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,DDOTI2_PressureVelocity,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,Ierr)

END FUNCTION DDOTI2_PressureVelocity

FUNCTION NormV_PressureVelocity(xV)

  REAL(RealKind) :: NormV_PressureVelocity
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:)

  REAL(RealKind) :: SUM,Temp
  INTEGER :: Ierr
  INTEGER :: ix,iy,iz
  INTEGER :: nx,ny,nz

  Temp=0.0e0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          Temp=Temp+SUM((xV(ibLoc)%p(:,ix,iy,iz)/(VolC(ix,iy,iz)+Eps))**2)
        END DO
      END DO
    END DO
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        Temp=Temp+(xV(ibLoc)%u_w(iy,iz)*FU(ix0,iy,iz) &
            /((VolC(ix0,iy,iz)+VolC(ix0+1,iy,iz)+Eps)/Two)**2)**2/Two
        Temp=Temp+(xV(ibLoc)%u_e(iy,iz)*FU(ix1,iy,iz) &
            /((VolC(ix1,iy,iz)+VolC(ix1+1,iy,iz)+Eps)/Two)**2)**2/Two
      END DO
    END DO
    DO iz=iz0+1,iz1
      DO ix=ix0+1,ix1
        Temp=Temp+(xV(ibLoc)%v_s(ix,iz)*FV(ix,iy0,iz) &
            /((VolC(ix,iy0,iz)+VolC(ix,iy0+1,iz)+Eps)/Two)**2)**2/Two
        Temp=Temp+(xV(ibLoc)%v_n(ix,iz)*FV(ix,iy1,iz) &
            /((VolC(ix,iy1,iz)+VolC(ix,iy1+1,iz)+Eps)/Two)**2)**2/Two
      END DO
    END DO
    DO iy=iy0+1,iy1
      DO ix=ix0+1,ix1
        Temp=Temp+(xV(ibLoc)%w_b(ix,iy)*FW(ix,iy,iz0) &
            /((VolC(ix,iy,iz0)+VolC(ix,iy,iz0+1)+Eps)/Two)**2)**2/Two
        Temp=Temp+(xV(ibLoc)%w_t(ix,iy)*FW(ix,iy,iz1) &
            /((VolC(ix,iy,iz1)+VolC(ix,iy,iz1+1)+Eps)/Two)**2)**2/Two
      END DO
    END DO
  END DO

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,NormV_PressureVelocity,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,Ierr)
  NormV_PressureVelocity=SQRT(NormV_PressureVelocity)

END FUNCTION NormV_PressureVelocity


FUNCTION DDOTR_PressureVelocity(xV,yV)

  REAL(RealKind) :: DDOTR_PressureVelocity
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:),yV(:)

  REAL(RealKind) :: SUM,Temp
  INTEGER :: Ierr
  INTEGER :: nx,ny,nz


  Temp=0.0e0
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    Temp = Temp + SUM(xV(ibLoc)%u_e * yV(ibLoc)%u_e)
    Temp = Temp + SUM(xV(ibLoc)%u_w * yV(ibLoc)%u_w)
    Temp = Temp + SUM(xV(ibLoc)%v_n * yV(ibLoc)%v_n)
    Temp = Temp + SUM(xV(ibLoc)%v_s * yV(ibLoc)%v_s)
    Temp = Temp + SUM(xV(ibLoc)%w_t * yV(ibLoc)%w_t)
    Temp = Temp + SUM(xV(ibLoc)%w_b * yV(ibLoc)%w_b)
  END DO

! -- Global scalar product over all processors --
  CALL MPI_Allreduce(Temp,DDOTR_PressureVelocity,1,MPI_RealKind, &
&                    MPI_SUM,MPI_Comm_World,Ierr)

END FUNCTION DDOTR_PressureVelocity

SUBROUTINE SCALE_PressureVelocity(alpha,yV)

  REAL(RealKind), INTENT(IN) :: alpha
  TYPE(PressureVelocity), INTENT(INOUT) :: yV(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    yV(ibLoc)%p   = alpha*yV(ibLoc)%p
    yV(ibLoc)%u_e = alpha*yV(ibLoc)%u_e
    yV(ibLoc)%u_w = alpha*yV(ibLoc)%u_w
    yV(ibLoc)%v_n = alpha*yV(ibLoc)%v_n
    yV(ibLoc)%v_s = alpha*yV(ibLoc)%v_s
    yV(ibLoc)%w_t = alpha*yV(ibLoc)%w_t
    yV(ibLoc)%w_b = alpha*yV(ibLoc)%w_b
  END DO

END SUBROUTINE SCALE_PressureVelocity

SUBROUTINE SCALE2_PressureVelocity(alpha,xV,yV)

  REAL(RealKind), INTENT(IN) :: alpha
  TYPE(PressureVelocity), INTENT(IN) :: xV(:)
  TYPE(PressureVelocity), INTENT(OUT) :: yV(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    yV(ibLoc)%p   = alpha*xV(ibLoc)%p
    yV(ibLoc)%u_e = alpha*xV(ibLoc)%u_e
    yV(ibLoc)%u_w = alpha*xV(ibLoc)%u_w
    yV(ibLoc)%v_n = alpha*xV(ibLoc)%v_n
    yV(ibLoc)%v_s = alpha*xV(ibLoc)%v_s
    yV(ibLoc)%w_t = alpha*xV(ibLoc)%w_t
    yV(ibLoc)%w_b = alpha*xV(ibLoc)%w_b
  END DO

END SUBROUTINE SCALE2_PressureVelocity

SUBROUTINE  DAXPY_PressureVelocity(alpha,xV,yV)

  REAL(RealKind) :: alpha
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:)
  TYPE(PressureVelocity), INTENT(INOUT)  :: yV(:)

  INTEGER :: nx,ny,nz

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    yV(ibLoc)%p   = alpha*xV(ibLoc)%p   + yV(ibLoc)%p
    yV(ibLoc)%u_e = alpha*xV(ibLoc)%u_e + yV(ibLoc)%u_e
    yV(ibLoc)%u_w = alpha*xV(ibLoc)%u_w + yV(ibLoc)%u_w
    yV(ibLoc)%v_n = alpha*xV(ibLoc)%v_n + yV(ibLoc)%v_n
    yV(ibLoc)%v_s = alpha*xV(ibLoc)%v_s + yV(ibLoc)%v_s
    yV(ibLoc)%w_t = alpha*xV(ibLoc)%w_t + yV(ibLoc)%w_t
    yV(ibLoc)%w_b = alpha*xV(ibLoc)%w_b + yV(ibLoc)%w_b
  END DO

END SUBROUTINE  DAXPY_PressureVelocity

SUBROUTINE  DXPBY_PressureVelocity(xV,alpha,yV)

  REAL(RealKind) :: alpha
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:)
  TYPE(PressureVelocity), INTENT(INOUT)  :: yV(:)

  INTEGER :: nx,ny,nz

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    yV(ibLoc)%p   = xV(ibLoc)%p   + alpha*yV(ibLoc)%p
    yV(ibLoc)%u_e = xV(ibLoc)%u_e + alpha*yV(ibLoc)%u_e
    yV(ibLoc)%u_w = xV(ibLoc)%u_w + alpha*yV(ibLoc)%u_w
    yV(ibLoc)%v_n = xV(ibLoc)%v_n + alpha*yV(ibLoc)%v_n
    yV(ibLoc)%v_s = xV(ibLoc)%v_s + alpha*yV(ibLoc)%v_s
    yV(ibLoc)%w_t = xV(ibLoc)%w_t + alpha*yV(ibLoc)%w_t
    yV(ibLoc)%w_b = xV(ibLoc)%w_b + alpha*yV(ibLoc)%w_b
  END DO

END SUBROUTINE  DXPBY_PressureVelocity

SUBROUTINE DXPY_PressureVelocity(xV,yV)

  TYPE(PressureVelocity), INTENT(IN)  :: xV(:)
  TYPE(PressureVelocity), INTENT(INOUT)  :: yV(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    yV(ibLoc)%p   = xV(ibLoc)%p   + yV(ibLoc)%p
    yV(ibLoc)%u_e = xV(ibLoc)%u_e + yV(ibLoc)%u_e
    yV(ibLoc)%u_w = xV(ibLoc)%u_w + yV(ibLoc)%u_w
    yV(ibLoc)%v_n = xV(ibLoc)%v_n + yV(ibLoc)%v_n
    yV(ibLoc)%v_s = xV(ibLoc)%v_s + yV(ibLoc)%v_s
    yV(ibLoc)%w_t = xV(ibLoc)%w_t + yV(ibLoc)%w_t
    yV(ibLoc)%w_b = xV(ibLoc)%w_b + yV(ibLoc)%w_b
  END DO

END SUBROUTINE DXPY_PressureVelocity

SUBROUTINE  DAXPBY_PressureVelocity(alpha,xV,beta,yV)

  REAL(RealKind) :: alpha,beta
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:)
  TYPE(PressureVelocity), INTENT(INOUT)  :: yV(:)

  INTEGER :: nx,ny,nz

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    yV(ibLoc)%p   = alpha*xV(ibLoc)%p   + beta*yV(ibLoc)%p
    yV(ibLoc)%u_e = alpha*xV(ibLoc)%u_e + beta*yV(ibLoc)%u_e
    yV(ibLoc)%u_w = alpha*xV(ibLoc)%u_w + beta*yV(ibLoc)%u_w
    yV(ibLoc)%v_n = alpha*xV(ibLoc)%v_n + beta*yV(ibLoc)%v_n
    yV(ibLoc)%v_s = alpha*xV(ibLoc)%v_s + beta*yV(ibLoc)%v_s
    yV(ibLoc)%w_t = alpha*xV(ibLoc)%w_t + beta*yV(ibLoc)%w_t
    yV(ibLoc)%w_b = alpha*xV(ibLoc)%w_b + beta*yV(ibLoc)%w_b
  END DO

END SUBROUTINE  DAXPBY_PressureVelocity

SUBROUTINE Deallocate_PressureVelocity(xV)

  TYPE(PressureVelocity), POINTER :: xV(:)

  IF (ASSOCIATED(xV)) THEN
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)

      IF (ASSOCIATED(xV(ibLoc)%p)) THEN
        DEALLOCATE(xV(ibLoc)%p)
      END IF

      IF (ASSOCIATED(xV(ibLoc)%u_e)) THEN
        DEALLOCATE(xV(ibLoc)%u_e)
      END IF
      IF (ASSOCIATED(xV(ibLoc)%u_w)) THEN
        DEALLOCATE(xV(ibLoc)%u_w)
      END IF

      IF (ASSOCIATED(xV(ibLoc)%v_s)) THEN
        DEALLOCATE(xV(ibLoc)%v_s)
      END IF
      IF (ASSOCIATED(xV(ibLoc)%v_n)) THEN
        DEALLOCATE(xV(ibLoc)%v_n)
      END IF

      IF (ASSOCIATED(xV(ibLoc)%w_t)) THEN
        DEALLOCATE(xV(ibLoc)%w_t)
      END IF
      IF (ASSOCIATED(xV(ibLoc)%w_b)) THEN
        DEALLOCATE(xV(ibLoc)%w_b)
      END IF
    END DO
    DEALLOCATE(xV)
    NULLIFY(xV)
  END IF
END SUBROUTINE Deallocate_PressureVelocity

SUBROUTINE Deallocate2_PressureVelocity(xV)

  TYPE(PressureVelocity), POINTER :: xV(:,:)

  INTEGER :: iDim,Dim

  IF (ASSOCIATED(xV)) THEN
    Dim=SIZE(xV,2)
    DO iDim=1,Dim
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        IF (ASSOCIATED(xV(ibLoc,iDim)%p)) THEN
          DEALLOCATE(xV(ibLoc,iDim)%p)
        END IF
        IF (ASSOCIATED(xV(ibLoc,iDim)%u_e)) THEN
          DEALLOCATE(xV(ibLoc,iDim)%u_e)
        END IF
        IF (ASSOCIATED(xV(ibLoc,iDim)%u_w)) THEN
          DEALLOCATE(xV(ibLoc,iDim)%u_w)
        END IF
        IF (ASSOCIATED(xV(ibLoc,iDim)%v_s)) THEN
          DEALLOCATE(xV(ibLoc,iDim)%v_s)
        END IF
        IF (ASSOCIATED(xV(ibLoc,iDim)%v_n)) THEN
          DEALLOCATE(xV(ibLoc,iDim)%v_n)
        END IF
        IF (ASSOCIATED(xV(ibLoc,iDim)%w_t)) THEN
          DEALLOCATE(xV(ibLoc,iDim)%w_t)
        END IF
        IF (ASSOCIATED(xV(ibLoc,iDim)%w_b)) THEN
          DEALLOCATE(xV(ibLoc,iDim)%w_b)
        END IF
      END DO
    END DO
    DEALLOCATE(xV)
    NULLIFY(xV)
  END IF
END SUBROUTINE Deallocate2_PressureVelocity

SUBROUTINE COPY_PressureVelocity(yV,xV)

  TYPE(PressureVelocity), INTENT(OUT)  :: yV(:)
  TYPE(PressureVelocity), INTENT(IN)  :: xV(:)

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    yV(ibLoc)%p   = xV(ibLoc)%p
    yV(ibLoc)%u_e = xV(ibLoc)%u_e
    yV(ibLoc)%u_w = xV(ibLoc)%u_w
    yV(ibLoc)%v_n = xV(ibLoc)%v_n
    yV(ibLoc)%v_s = xV(ibLoc)%v_s
    yV(ibLoc)%w_t = xV(ibLoc)%w_t
    yV(ibLoc)%w_b = xV(ibLoc)%w_b
  END DO

END SUBROUTINE COPY_PressureVelocity

SUBROUTINE Value_PressureVelocity(yV,Value)

  TYPE(PressureVelocity), INTENT(OUT)  :: yV(:)
  REAL(RealKind), INTENT(IN)  :: Value

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    yV(ibLoc)%p   = Value
    yV(ibLoc)%u_e = Value
    yV(ibLoc)%u_w = Value
    yV(ibLoc)%v_n = Value
    yV(ibLoc)%v_s = Value
    yV(ibLoc)%w_t = Value
    yV(ibLoc)%w_b = Value
  END DO

END SUBROUTINE Value_PressureVelocity

SUBROUTINE Output_PressureVelocity(yV,FileName)

  INTEGER :: iz

  TYPE(PressureVelocity), INTENT(IN)  :: yV(:)
  CHARACTER(80), INTENT(IN), OPTIONAL :: FileName

  IF (PRESENT(Filename)) THEN
    OPEN(UNIT=10,FILE=FileName,STATUS='UNKNOWN')
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      WRITE(10,*) 'BlockNumber',ib
      WRITE(10,*) 'c' 
      WRITE(10,*) SUM(ABS(yV(ibLoc)%p)) 
      WRITE(10,*) 'e' 
      WRITE(10,*) SUM(ABS(yV(ibLoc)%u_e))
      WRITE(10,*) 'w' 
      WRITE(10,*) SUM(ABS(yV(ibLoc)%u_w))
      WRITE(10,*) 'n' 
      WRITE(10,*) SUM(ABS(yV(ibLoc)%v_n))
      WRITE(10,*) 's' 
      WRITE(10,*) SUM(ABS(yV(ibLoc)%v_s))
      WRITE(10,*) 't' 
      WRITE(10,*) SUM(ABS(yV(ibLoc)%w_t))
      WRITE(10,*) 'b' 
      WRITE(10,*) SUM(ABS(yV(ibLoc)%w_b)) 
    END DO
    CLOSE(10)
  ELSE  
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      WRITE(*,*) 'BlockNumber',ib
      WRITE(*,*) 'c' 
      WRITE(*,*) (SUM(ABS(yV(ibLoc)%p(1,:,:,iz))),iz=1,SIZE(yV(ibLoc)%p,4))
      WRITE(*,*) 'e' 
      WRITE(*,*) (SUM(ABS(yV(ibLoc)%u_e(:,iz))),iz=1,SIZE(yV(ibLoc)%u_e,2))
      WRITE(*,*) 'w' 
      WRITE(*,*) SUM(ABS(yV(ibLoc)%u_w))
      WRITE(*,*) 'n' 
      WRITE(*,*) SUM(ABS(yV(ibLoc)%v_n))
      WRITE(*,*) 's' 
      WRITE(*,*) SUM(ABS(yV(ibLoc)%v_s))
      WRITE(*,*) 't' 
      WRITE(*,*) SUM(ABS(yV(ibLoc)%w_t))
      WRITE(*,*) 'b' 
      WRITE(*,*) SUM(ABS(yV(ibLoc)%w_b)) 
    END DO
  END IF  

END SUBROUTINE Output_PressureVelocity

SUBROUTINE Output_Pressure(yV,FileName)

  TYPE(PressureVelocity), INTENT(IN)  :: yV(:)
  CHARACTER(*), INTENT(IN)  :: FileName

  INTEGER :: ix,iy,iz

  OPEN(UNIT=10,FILE=FileName,STATUS='UNKNOWN')
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          IF (yV(ibLoc)%p(1,ix,iy,iz)/=Zero) THEN
            WRITE(10,*) ix,iy,iz,yV(ibLoc)%p(1,ix,iy,iz)
          END IF  
        END DO
      END DO
    END DO
  END DO
  CLOSE(10)

END SUBROUTINE Output_Pressure

SUBROUTINE Output2_PressureVelocity(xV,yV,FileName)

  TYPE(PressureVelocity), INTENT(IN)  :: xV(:)
  TYPE(PressureVelocity), INTENT(IN)  :: yV(:)
  CHARACTER(80), INTENT(IN)  :: FileName

  OPEN(UNIT=10,FILE=FileName,STATUS='UNKNOWN')
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    WRITE(10,*) 'BlockNumber',ib
    WRITE(10,*) 'c'
    WRITE(10,*) SUM(yV(ibLoc)%p*xV(ibLoc)%p)
    WRITE(10,*) 'e'
    WRITE(10,*) SUM(yV(ibLoc)%u_e*xV(ibLoc)%u_e)
    WRITE(10,*) 'w'
    WRITE(10,*) SUM(yV(ibLoc)%u_w*xV(ibLoc)%u_w)
    WRITE(10,*) 'n'
    WRITE(10,*) SUM(yV(ibLoc)%v_n*xV(ibLoc)%v_n)
    WRITE(10,*) 's'
    WRITE(10,*) SUM(yV(ibLoc)%v_s*xV(ibLoc)%v_s)
    WRITE(10,*) 't'
    WRITE(10,*) SUM(yV(ibLoc)%w_t*xV(ibLoc)%w_t)
    WRITE(10,*) 'b'
    WRITE(10,*) SUM(yV(ibLoc)%w_b*xV(ibLoc)%w_b)
  END DO
  CLOSE(10)

END SUBROUTINE Output2_PressureVelocity

SUBROUTINE Allocate_PressureVelocity(xV)

  TYPE(PressureVelocity), POINTER :: xV(:)

  ALLOCATE(xV(nbLoc))

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    
    xV(ibLoc)%refine=Refine
    xV(ibLoc)%nx=nx
    xV(ibLoc)%ny=ny
    xV(ibLoc)%nz=nz

    ALLOCATE(xV(ibLoc)%p(1:np1,ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
    xV(ibLoc)%p=0.0e0

    ALLOCATE(xV(ibLoc)%u_w(iy0+1:iy1,iz0+1:iz1))
    xV(ibLoc)%u_w=0.0e0

    ALLOCATE(xV(ibLoc)%u_e(iy0+1:iy1,iz0+1:iz1))
    xV(ibLoc)%u_e=0.0e0

    ALLOCATE(xV(ibLoc)%v_s(ix0+1:ix1,iz0+1:iz1))
    xV(ibLoc)%v_s=0.0e0

    ALLOCATE(xV(ibLoc)%v_n(ix0+1:ix1,iz0+1:iz1))
    xV(ibLoc)%v_n=0.0e0

    ALLOCATE(xV(ibLoc)%w_b(ix0+1:ix1,iy0+1:iy1))
    xV(ibLoc)%w_b=0.0e0

    ALLOCATE(xV(ibLoc)%w_t(ix0+1:ix1,iy0+1:iy1))
    xV(ibLoc)%w_t=0.0e0
  END DO

END SUBROUTINE Allocate_PressureVelocity

SUBROUTINE Allocate_Velocity(xV)

  TYPE(PressureVelocity), POINTER :: xV(:)

  ALLOCATE(xV(nbLoc))

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    
    xV(ibLoc)%refine=Refine
    xV(ibLoc)%nx=nx
    xV(ibLoc)%ny=ny
    xV(ibLoc)%nz=nz

    ALLOCATE(xV(ibLoc)%u_w(iy0+1:iy1,iz0+1:iz1))
    xV(ibLoc)%u_w=0.0e0

    ALLOCATE(xV(ibLoc)%u_e(iy0+1:iy1,iz0+1:iz1))
    xV(ibLoc)%u_e=0.0e0

    ALLOCATE(xV(ibLoc)%v_s(ix0+1:ix1,iz0+1:iz1))
    xV(ibLoc)%v_s=0.0e0

    ALLOCATE(xV(ibLoc)%v_n(ix0+1:ix1,iz0+1:iz1))
    xV(ibLoc)%v_n=0.0e0

    ALLOCATE(xV(ibLoc)%w_b(ix0+1:ix1,iy0+1:iy1))
    xV(ibLoc)%w_b=0.0e0

    ALLOCATE(xV(ibLoc)%w_t(ix0+1:ix1,iy0+1:iy1))
    xV(ibLoc)%w_t=0.0e0
  END DO

END SUBROUTINE Allocate_Velocity

SUBROUTINE Allocate_Pressure(xV)

  TYPE(Pressure_T), POINTER :: xV(:)

  ALLOCATE(xV(nbLoc))

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    
    ALLOCATE(xV(ibLoc)%p(ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
    xV(ibLoc)%p=Zero

  END DO

END SUBROUTINE Allocate_Pressure

SUBROUTINE Allocate2_PressureVelocity(xV,Dim)

  TYPE(PressureVelocity), POINTER :: xV(:,:)
  INTEGER :: Dim

  INTEGER :: iDim

  ALLOCATE(xV(nbLoc,Dim))

  DO iDim=1,Dim
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))

      xV(ibLoc,iDim)%refine=Refine
      xV(ibLoc,iDim)%nx=nx
      xV(ibLoc,iDim)%ny=ny
      xV(ibLoc,iDim)%nz=nz

      ALLOCATE(xV(ibLoc,iDim)%p(1:np1,ix0+1:ix1,iy0+1:iy1,iz0+1:iz1))
      xV(ibLoc,iDim)%p=0.0e0

      ALLOCATE(xV(ibLoc,iDim)%u_w(iy0+1:iy1,iz0+1:iz1))
      xV(ibLoc,iDim)%u_w=0.0e0

      ALLOCATE(xV(ibLoc,iDim)%u_e(iy0+1:iy1,iz0+1:iz1))
      xV(ibLoc,iDim)%u_e=0.0e0

      ALLOCATE(xV(ibLoc,iDim)%v_s(ix0+1:ix1,iz0+1:iz1))
      xV(ibLoc,iDim)%v_s=0.0e0

      ALLOCATE(xV(ibLoc,iDim)%v_n(ix0+1:ix1,iz0+1:iz1))
      xV(ibLoc,iDim)%v_n=0.0e0

      ALLOCATE(xV(ibLoc,iDim)%w_b(ix0+1:ix1,iy0+1:iy1))
      xV(ibLoc,iDim)%w_b=0.0e0

      ALLOCATE(xV(ibLoc,iDim)%w_t(ix0+1:ix1,iy0+1:iy1))
      xV(ibLoc,iDim)%w_t=0.0e0
    END DO
  END DO
END SUBROUTINE Allocate2_PressureVelocity

SUBROUTINE RandCopy(vec)

!---------------------------------------------------------------
!---  Copy Concentrations from Neighboring Blocks
!---------------------------------------------------------------

  TYPE (PressureVelocity), TARGET :: vec(:)

  INTEGER :: in
  INTEGER :: jx,jy,jz
  TYPE (PressureVelocity), POINTER :: subblock
  TYPE (Nachbar_T), POINTER :: Nachbar

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    subblock=>vec(ibLoc)
    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)

      IF (Nachbar%nType=='iw'.AND.ibLoc<=ibnLoc) THEN           ! west neighbor
        CALL CopyFluxVel(subblock%u_w(jy0+1:jy1,jz0+1:jz1), &
                         vec(ibnLoc)%u_e(jNy0+1:jNy1,jNz0+1:jNz1))
      ELSE IF (Nachbar%nType=='ie'.AND.ibLoc<ibnLoc) THEN       ! east neighbor
        CALL CopyFluxVel(subblock%u_e(jy0+1:jy1,jz0+1:jz1), &
                         vec(ibnLoc)%u_w(jNy0+1:jNy1,jNz0+1:jNz1))
      ELSE IF (Nachbar%nType=='is'.AND.ibLoc<=ibnLoc) THEN      ! south neighbor
        CALL CopyFluxVel(subblock%v_s(jx0+1:jx1,jz0+1:jz1),  &
                         vec(ibnLoc)%v_n(jNx0+1:jNx1,jNz0+1:jNz1))
      ELSE IF (Nachbar%nType=='in'.AND.ibLoc<ibnLoc) THEN       ! north neighbor
        CALL CopyFluxVel(subblock%v_n(jx0+1:jx1,jz0+1:jz1),  &
                         vec(ibnLoc)%v_s(jNx0+1:jNx1,jNz0+1:jNz1))
      ELSE IF (Nachbar%nType=='ib'.AND.ibLoc<=ibnLoc) THEN      ! bottom neighbor
        CALL CopyFluxVel(subblock%w_b(jx0+1:jx1,jy0+1:jy1), &
                         vec(ibnLoc)%w_t(jNx0+1:jNx1,jNy0+1:jNy1))
      ELSE IF (Nachbar%nType=='it'.AND.ibLoc<ibnLoc) THEN       ! top neighbor
        CALL CopyFluxVel(subblock%w_t(jx0+1:jx1,jy0+1:jy1), &
                         vec(ibnLoc)%w_b(jNx0+1:jNx1,jNy0+1:jNy1))
      END IF
    END DO ! in
  END DO   ! ibLoc

 END SUBROUTINE RandCopy


 SUBROUTINE Exchange(b)
   
   TYPE (PressureVelocity) :: b(:) 

   INTEGER :: i,ip,j   
   INTEGER :: tag
   TYPE(MPI_Request) :: req(2*NumberNeiProc)
   TYPE(MPI_Status) :: StatusArray(2*NumberNeiProc)

   REAL(RealKind) :: PutBuf(MaxPutBuf,NumberNeiProc)
   REAL(RealKind) :: GetBuf(MaxGetBuf,NumberNeiProc)
   
   CALL RandCopy(b)
   i=0
   DO ip=0,NumProcs-1
      IF (chain(ip)%nlink > 0) THEN
         i=i+1
         CALL PackBuff(chain(ip)%glied,PutBuf(:,i),b)
      END IF
   END DO
                         
   i=0
   DO ip=0,NumProcs-1
      IF (chain(ip)%nlink > 0) THEN
         i=i+1
         tag=MyId+NumProcs*ip
         CALL MPI_IRECV(GetBuf(1,i),chain(ip)%LenUnPack, &
                       MPI_RealKind,ip,tag,      &
                       MPI_COMM_WORLD,req(i),MPIErr)
      END IF
   END DO
   i=0
   DO ip=0,NumProcs-1
      IF (chain(ip)%nlink > 0) THEN
         i=i+1
         tag=ip+NumProcs*MyId
         CALL MPI_ISEND(PutBuf(1,i),chain(ip)%LenPack, &
                       MPI_RealKind,ip,tag,      &
                       MPI_COMM_WORLD,req(NumberNeiProc+i),MPIErr)
      END IF
   END DO

   CALL MPI_WAITALL(2*NumberNeiProc,req,StatusArray,MPIErr)

   i=0
      
   DO ip=0,NumProcs-1
      IF (chain(ip)%nlink > 0) THEN
         i=i+1
         CALL UnPackBuff(chain(ip)%glied,GetBuf(:,i),b)
      END IF
   END DO

CONTAINS 

SUBROUTINE PackBuff(chain,Buffer,b)

  TYPE (gliedT) :: chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (PressureVelocity) :: b(:) ! TYPE (Vector) b

  INTEGER :: iBuf
  INTEGER :: i, len

  iBuf=0
  Buffer(:)=Zero
  len=SIZE(chain)
  DO i=1,len
    CALL Set(Chain(i))
    CALL Set(Floor(ibC))
    CALL Set(Nachbar)
    IF (nType=='pw') THEN
      CALL PackFluxVel(Buffer,iBuf,b(ibCLoc)%u_w(jy0+1:jy1,jz0+1:jz1),CopyCase)
    ELSE IF (nType=='pe') THEN
      CALL PackFluxVel(Buffer,iBuf,b(ibCLoc)%u_e(jy0+1:jy1,jz0+1:jz1),CopyCase)
    ELSE IF (nType=='ps') THEN
      CALL PackFluxVel(Buffer,iBuf,b(ibCLoc)%v_s(jx0+1:jx1,jz0+1:jz1),CopyCase)
    ELSE IF (nType=='pn') THEN
      CALL PackFluxVel(Buffer,iBuf,b(ibCLoc)%v_n(jx0+1:jx1,jz0+1:jz1),CopyCase)
    ELSE IF (nType=='pb') THEN
      CALL PackFluxVel(Buffer,iBuf,b(ibCLoc)%w_b(jx0+1:jx1,jy0+1:jy1),CopyCase)
    ELSE IF (nType=='pt') THEN
      CALL PackFluxVel(Buffer,iBuf,b(ibCLoc)%w_t(jx0+1:jx1,jy0+1:jy1),CopyCase)
    END IF
  END DO   ! i=1,SIZE(chain)
END SUBROUTINE PackBuff

SUBROUTINE UnPackBuff(chain,Buffer,b)

  TYPE (gliedT) :: chain(:)
  REAL(RealKind) :: Buffer(:)
  TYPE (PressureVelocity) :: b(:) ! TYPE (Vector) b

  INTEGER :: iBuf
  INTEGER :: i,len

  iBuf=0
  len=SIZE(chain)
  DO i=1,len
    CALL Set(Chain(i))
    CALL Set(Floor(ibC))
    CALL Set(Nachbar)

    IF (nType=='pw') THEN
      CALL UnpackFluxVel(b(ibCLoc)%u_w(jy0+1:jy1,jz0+1:jz1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='pe') THEN
      CALL UnpackFluxVel(b(ibCLoc)%u_e(jy0+1:jy1,jz0+1:jz1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='ps') THEN
      CALL UnpackFluxVel(b(ibCLoc)%v_s(jx0+1:jx1,jz0+1:jz1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='pn') THEN
      CALL UnpackFluxVel(b(ibCLoc)%v_n(jx0+1:jx1,jz0+1:jz1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='pb') THEN
      CALL UnpackFluxVel(b(ibCLoc)%w_b(jx0+1:jx1,jy0+1:jy1),Buffer,iBuf,CopyCase)
    ELSE IF (nType=='pt') THEN
      CALL UnpackFluxVel(b(ibCLoc)%w_t(jx0+1:jx1,jy0+1:jy1),Buffer,iBuf,CopyCase)
    END IF
  END DO   ! i=1,SIZE(chain)
END SUBROUTINE UnPackBuff

END SUBROUTINE Exchange 
END MODULE VecPrVel_Mod

