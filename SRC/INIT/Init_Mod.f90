MODULE Init_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE Control_Mod
  USE Chemie_Mod
  USE DataType_Mod

  IMPLICIT NONE

  INTEGER, PRIVATE :: i,in
  REAL(RealKind), PRIVATE :: xPL,yPL,zPL,dzLoc
  REAL(RealKind), POINTER, PRIVATE :: c(:,:,:,:) 

  INTERFACE VectorInit
    MODULE PROCEDURE VectorInitFunction,VectorInitValue
  END INTERFACE VectorInit
  INTERFACE ScalarInit
    MODULE PROCEDURE ScalarInit
  END INTERFACE ScalarInit

CONTAINS

SUBROUTINE InputInital(FileName,Vector)

  CHARACTER(*) :: FileName
  TYPE (Vector4Cell_T), TARGET :: Vector(:)

  INTEGER :: Pos
  REAL(RealKind) :: Value
  CHARACTER(20) :: Species
  INTEGER, PARAMETER :: Unit=10

  OPEN(UNIT=Unit,FILE=TRIM(FileName),STATUS='OLD')

  DO 
    READ(Unit,*,END=100) Species
    READ(Unit,*,END=100) Value
    Pos=Position(Species)
    CALL VectorInit(Pos,Vector,Value)
    READ(Unit,*,END=100) 
  END DO

100 CONTINUE

END SUBROUTINE InputInital


SUBROUTINE ScalarInit(Scalar,Val,Time)
!
!==================================================
!----  Initialization of Scalar Values
!==================================================

  TYPE (ScalarCell_T), TARGET :: Scalar(:)
  REAL(RealKind) :: Val
  EXTERNAL :: Val
  REAL(RealKind) :: Time
  INTEGER :: ix, iy, iz

  DO ibLoc=1,nbLoc
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))
    c=>Scalar(ibLoc)%c
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          c(ix,iy,iz,1) = Val(xPL,yPL,zPL,zH(ix,iy),Time) &
            *VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        END DO
      END DO
    END DO
    IF (.NOT.GradFull) THEN
      DO i=1,NumBoundCell
        ix=BoundCell(i)%ix 
        iy=BoundCell(i)%iy 
        iz=BoundCell(i)%iz 
        xPL=BoundCell(i)%xS
        yPL=BoundCell(i)%yS
        zPL=BoundCell(i)%zS
        c(ix,iy,iz,1) = Val(xPL,yPL,zPL,zH(ix,iy),Time) &
          *VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END IF
    c(ix0,:,:,:)=c(ix0+1,:,:,:)
    c(ix1+1,:,:,:)=c(ix1,:,:,:)
    c(:,iy0,:,:)=c(:,iy0+1,:,:)
    c(:,iy1+1,:,:)=c(:,iy1,:,:)
    c(:,:,iz0,:)=c(:,:,iz0+1,:)
    c(:,:,iz1+1,:)=c(:,:,iz1,:)

  END DO   ! ibLoc

!-----------------------------------------------------------
END SUBROUTINE ScalarInit

SUBROUTINE BoundaryInit(Type,Val,Time)
!
!==================================================
!----  Initialization of Scalar Values
!==================================================

  CHARACTER(*) :: Type
  REAL(RealKind) :: Val
  EXTERNAL :: Val
  REAL(RealKind) :: Time
  INTEGER :: ix, iy, iz, i

  IF(Type=='Skin') THEN
    DO ibLoc=1,nbLoc
      ib = LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO i=1,NumBoundCell
        ix     = BoundCell(i)%ix
        iy     = BoundCell(i)%iy
        iz     = BoundCell(i)%iz
        xPL=xP(ix-1)+0.5e0*dx(ix)
        yPL=yP(iy-1)+0.5e0*dy(iy)
        zPL=zP(iz-1)
        BoundCell(i)%TeS=Val(xPL,yPL,zPL,zH(ix,iy),Time)
      END DO
    END DO   ! ibLoc
  END IF
!-----------------------------------------------------------
END SUBROUTINE BoundaryInit

SUBROUTINE VectorInitFunction(Pos,Vector,Val,Time)
!
!==================================================
!----  Initialization of Scalar Values
!==================================================
!

  INTEGER :: Pos
  TYPE (Vector4Cell_T), TARGET :: Vector(:)
  REAL(RealKind) :: Time
  REAL(RealKind) :: Val
  EXTERNAL Val
  INTEGER :: ix, iy, iz


    
  DO ibLoc=1,nbLoc
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))
    c=>Vector(ibLoc)%Vec(Pos)%c
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          c(ix,iy,iz,1) = Val(xPL,yPL,zPL,zH(ix,iy),Time) &
            *VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        END DO
      END DO
    END DO
    IF (.NOT.GradFull) THEN
      DO i=1,NumBoundCell
        ix=BoundCell(i)%ix
        iy=BoundCell(i)%iy
        iz=BoundCell(i)%iz
        xPL=BoundCell(i)%xS
        yPL=BoundCell(i)%yS
        zPL=BoundCell(i)%zS
        c(ix,iy,iz,1) = Val(xPL,yPL,zPL,zH(ix,iy),Time) &
          *VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
      END DO
    END IF
    c(ix0,:,:,:)=c(ix0+1,:,:,:)
    c(ix1+1,:,:,:)=c(ix1,:,:,:)
    c(:,iy0,:,:)=c(:,iy0+1,:,:)
    c(:,iy1+1,:,:)=c(:,iy1,:,:)
    c(:,:,iz0,:)=c(:,:,iz0+1,:)
    c(:,:,iz1+1,:)=c(:,:,iz1,:)

  END DO   ! ibLoc

!-----------------------------------------------------------
END SUBROUTINE VectorInitFunction

SUBROUTINE VectorInitValue(Pos,Vector,Vel)
!
!==================================================
!----  Initialization of Scalar Values
!==================================================
!

  INTEGER :: Pos
  TYPE (Vector4Cell_T) :: Vector(:)
  REAL(RealKind) :: Vel
  REAL(RealKind) :: Time

  INTEGER :: ix, iy, iz

  DO ibLoc=1,nbLoc
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))
    c=>Vector(ibLoc)%Vec(Pos)%c
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          c(ix,iy,iz,1) = Vel &
            *VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        END DO
      END DO
    END DO
    c(ix0,:,:,:)=c(ix0+1,:,:,:)
    c(ix1+1,:,:,:)=c(ix1,:,:,:)
    c(:,iy0,:,:)=c(:,iy0+1,:,:)
    c(:,iy1+1,:,:)=c(:,iy1,:,:)
    c(:,:,iz0,:)=c(:,:,iz0+1,:)
    c(:,:,iz1+1,:)=c(:,:,iz1,:)
  END DO   ! ibLoc

!-----------------------------------------------------------
END SUBROUTINE VectorInitValue

SUBROUTINE VelocityInit(Vel,U,V,W,Time)
!
!==================================================
!----  Initialization of Velocity Values
!==================================================
!

  TYPE (VelocityFace_T) :: Vel(:)
  REAL(RealKind) :: U,V,W
  EXTERNAL U,V,W
  REAL(RealKind) :: Time

!---  local variables
  INTEGER :: ix, iy, iz
  INTEGER :: jx, jy, jz
  TYPE(Nachbar_T), POINTER :: Nachbar
    

!------------------------------------------------------------
!+++++++++  Allocation of Velocity Variables ++++++++++++++++
!------------------------------------------------------------
!--- local block structure --


  DO ibLoc=1,nbLoc
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))

    DO ix=ix0,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          xPL=xP(ix) 
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          Vel(ibLoc)%uF(ix,iy,iz) = U(xPL,yPL,zPL,zH(ix,iy),Time)
        END DO
      END DO
    END DO
    DO i=1,NumBoundFaceU
      ix=BoundFaceU(i)%ix
      iy=BoundFaceU(i)%iy
      iz=BoundFaceU(i)%iz
      xPL=BoundFaceU(i)%xS
      yPL=BoundFaceU(i)%yS
      zPL=BoundFaceU(i)%zS
      Vel(ibLoc)%uF(ix,iy,iz) = U(xPL,yPL,zPL,zH(ix,iy),Time)
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0,iy1
        DO iz=iz0+1,iz1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy) 
          zPL=zP(iz-1)+0.5e0*dz(iz)
          Vel(ibLoc)%vF(ix,iy,iz) = V(xPL,yPL,zPL,zH(ix,iy),Time)
        END DO
      END DO
    END DO
    DO i=1,NumBoundFaceV
      ix=BoundFaceV(i)%ix
      iy=BoundFaceV(i)%iy
      iz=BoundFaceV(i)%iz
      xPL=BoundFaceV(i)%xS
      yPL=BoundFaceV(i)%yS
      zPL=BoundFaceV(i)%zS
      Vel(ibLoc)%vF(ix,iy,iz) = V(xPL,yPL,zPL,zH(ix,iy),Time)
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0,iz1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz) 
          Vel(ibLoc)%wF(ix,iy,iz) = W(xPL,yPL,zPL,zH(ix,iy),Time)
        END DO
      END DO
    END DO
    DO i=1,NumBoundFaceW
      ix=BoundFaceW(i)%ix
      iy=BoundFaceW(i)%iy
      iz=BoundFaceW(i)%iz
      xPL=BoundFaceW(i)%xS
      yPL=BoundFaceW(i)%yS
      zPL=BoundFaceW(i)%zS
      Vel(ibLoc)%wF(ix,iy,iz) = W(xPL,yPL,zPL,zH(ix,iy),Time)
    END DO

    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)

! --  Westlicher Rand --
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar<Refine) THEN
          xPL=xP(ix0)
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              yPL=yP(jy)
              zPL=zP(jz)
              Vel(ibLoc)%uF(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                U(xPL,yPL,zPL,zH(ix0,jy),Time)
            END DO
          END DO
        END IF
      END IF
! --  Oestlicher Rand --
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar<Refine) THEN
          xPL=xP(ix1)
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              yPL=yP(jy)
              zPL=zP(jz)
              Vel(ibLoc)%uF(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                U(xPL,yPL,zPL,zH(ix1,jy),Time)
            END DO
          END DO
        END IF
      END IF
! --  Suedlicher Rand --
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar<Refine) THEN
          yPL=yP(iy0)
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              xPL=xP(jx)
              zPL=zP(jz)
              Vel(ibLoc)%vF(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)= &
                V(xPL,yPL,zPL,zH(jx,iy0),Time)
            END DO
          END DO
        END IF
      END IF
! --  Noerdlicher Rand --
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar<Refine) THEN
          yPL=yP(iy1)
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              xPL=xP(jx)
              zPL=zP(jz)
              Vel(ibLoc)%vF(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
                V(xPL,yPL,zPL,zH(jx,iy1),Time)
            END DO
          END DO
        END IF
      END IF
! --  Unterer Rand --
      IF (Nachbar%nType(2:2) == 'b') THEN
        IF (RefineNachbar<Refine) THEN
          zPL=zP(iz0)
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              xPL=xP(jx)
              yPL=yP(jy)
              Vel(ibLoc)%wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)= &
                W(xPL,yPL,zPL,zH(jx,jy),Time)
            END DO
          END DO
        END IF
      END IF
! --  Oberer Rand --
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar<Refine) THEN
          zPL=zP(iz1)
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              xPL=xP(jx)
              yPL=yP(jy)
              Vel(ibLoc)%wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
                W(xPL,yPL,zPL,zH(jx,jy),Time)
            END DO
          END DO
        END IF
      END IF
    END DO
    Vel(ibLoc)%uF=Vel(ibLoc)%uF*FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1) &
              /(FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
    Vel(ibLoc)%vF=Vel(ibLoc)%vF*FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1) &
              /(FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)+Eps)
    Vel(ibLoc)%wF=Vel(ibLoc)%wF*FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1) &
              /(FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)+Eps)

  END DO   ! ibLoc

!-----------------------------------------------------------
END SUBROUTINE VelocityInit

SUBROUTINE RhoVelocityInit(Vel,U,V,W,Rho,Time)
!
!==================================================
!----  Initialization of Velocity Values
!==================================================
!

  TYPE (VelocityFace_T) :: Vel(:)
  REAL(RealKind) :: U,V,W,Rho
  EXTERNAL U,V,W,Rho
  REAL(RealKind) :: Time

!---  local variables
  INTEGER :: ix, iy, iz
  INTEGER :: jx, jy, jz
  TYPE(Nachbar_T), POINTER :: Nachbar
    

!------------------------------------------------------------
!+++++++++  Allocation of Velocity Variables ++++++++++++++++
!------------------------------------------------------------
!--- local block structure --


  DO ibLoc=1,nbLoc
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))

    DO ix=ix0,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          xPL=xP(ix) 
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz-1)+0.5e0*dz(iz)
          Vel(ibLoc)%uF(ix,iy,iz) = U(xPL,yPL,zPL,zH(ix,iy),Time)*Rho(xPL,yPL,zPL,zH(ix,iy),Time)
        END DO
      END DO
    END DO
    DO i=1,NumBoundFaceU
      ix=BoundFaceU(i)%ix
      iy=BoundFaceU(i)%iy
      iz=BoundFaceU(i)%iz
      xPL=BoundFaceU(i)%xS
      yPL=BoundFaceU(i)%yS
      zPL=BoundFaceU(i)%zS
!     WRITE(*,*) 'uF',LBOUND(Vel(ibLoc)%uF,1)
!     WRITE(*,*) 'uF',ix,iy,iz
      Vel(ibLoc)%uF(ix,iy,iz) = U(xPL,yPL,zPL,zH(ix,iy),Time)*Rho(xPL,yPL,zPL,zH(ix,iy),Time)
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0,iy1
        DO iz=iz0+1,iz1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy) 
          zPL=zP(iz-1)+0.5e0*dz(iz)
          Vel(ibLoc)%vF(ix,iy,iz) = V(xPL,yPL,zPL,zH(ix,iy),Time)*Rho(xPL,yPL,zPL,zH(ix,iy),Time)
        END DO
      END DO
    END DO
    DO i=1,NumBoundFaceV
      ix=BoundFaceV(i)%ix
      iy=BoundFaceV(i)%iy
      iz=BoundFaceV(i)%iz
      xPL=BoundFaceV(i)%xS
      yPL=BoundFaceV(i)%yS
      zPL=BoundFaceV(i)%zS
      Vel(ibLoc)%vF(ix,iy,iz) = V(xPL,yPL,zPL,zH(ix,iy),Time)*Rho(xPL,yPL,zPL,zH(ix,iy),Time)
    END DO
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0,iz1
          xPL=xP(ix-1)+0.5e0*dx(ix)
          yPL=yP(iy-1)+0.5e0*dy(iy)
          zPL=zP(iz) 
          Vel(ibLoc)%wF(ix,iy,iz) = W(xPL,yPL,zPL,zH(ix,iy),Time)*Rho(xPL,yPL,zPL,zH(ix,iy),Time)
        END DO
      END DO
    END DO
    DO i=1,NumBoundFaceW
      ix=BoundFaceW(i)%ix
      iy=BoundFaceW(i)%iy
      iz=BoundFaceW(i)%iz
      xPL=BoundFaceW(i)%xS
      yPL=BoundFaceW(i)%yS
      zPL=BoundFaceW(i)%zS
      Vel(ibLoc)%wF(ix,iy,iz) = W(xPL,yPL,zPL,zH(ix,iy),Time)*Rho(xPL,yPL,zPL,zH(ix,iy),Time)
    END DO

    DO in=1,AnzahlNachbar
      Nachbar=>Nachbars(in)
      CALL Set(Nachbar)

! --  Westlicher Rand --
      IF (Nachbar%nType(2:2) == 'w') THEN
        IF (RefineNachbar<Refine) THEN
          xPL=xP(ix0)
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              yPL=yP(jy)
              zPL=zP(jz)
              Vel(ibLoc)%uF(ix0,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                U(xPL,yPL,zPL,zH(ix0,jy),Time)*Rho(xPL,yPL,zPL,zH(ix0,jy),Time)
            END DO
          END DO
        END IF
      END IF
! --  Oestlicher Rand --
      IF (Nachbar%nType(2:2) == 'e') THEN
        IF (RefineNachbar<Refine) THEN
          xPL=xP(ix1)
          DO jy=jy0+1,jy1,IncrY
            DO jz=jz0+1,jz1,IncrZ
              yPL=yP(jy)
              zPL=zP(jz)
              Vel(ibLoc)%uF(ix1,jy:jy+IncrY-1,jz:jz+IncrZ-1)= &
                U(xPL,yPL,zPL,zH(ix1,jy),Time)*Rho(xPL,yPL,zPL,zH(ix1,jy),Time)
            END DO
          END DO
        END IF
      END IF
! --  Suedlicher Rand --
      IF (Nachbar%nType(2:2) == 's') THEN
        IF (RefineNachbar<Refine) THEN
          yPL=yP(iy0)
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              xPL=xP(jx)
              zPL=zP(jz)
              Vel(ibLoc)%vF(jx:jx+IncrX-1,iy0,jz:jz+IncrZ-1)= &
                V(xPL,yPL,zPL,zH(jx,iy0),Time)*Rho(xPL,yPL,zPL,zH(jx,iy0),Time)
            END DO
          END DO
        END IF
      END IF
! --  Noerdlicher Rand --
      IF (Nachbar%nType(2:2) == 'n') THEN
        IF (RefineNachbar<Refine) THEN
          yPL=yP(iy1)
          DO jx=jx0+1,jx1,IncrX
            DO jz=jz0+1,jz1,IncrZ
              xPL=xP(jx)
              zPL=zP(jz)
              Vel(ibLoc)%vF(jx:jx+IncrX-1,iy1,jz:jz+IncrZ-1)= &
                V(xPL,yPL,zPL,zH(jx,iy1),Time)*Rho(xPL,yPL,zPL,zH(jx,iy1),Time)
            END DO
          END DO
        END IF
      END IF
! --  Unterer Rand --
      IF (Nachbar%nType(2:2) == 'b') THEN
        IF (RefineNachbar<Refine) THEN
          zPL=zP(iz0)
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              xPL=xP(jx)
              yPL=yP(jy)
              Vel(ibLoc)%wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz0)= &
                W(xPL,yPL,zPL,zH(jx,jy),Time)*Rho(xPL,yPL,zPL,zH(jx,jy),Time)
            END DO
          END DO
        END IF
      END IF
! --  Oberer Rand --
      IF (Nachbar%nType(2:2) == 't') THEN
        IF (RefineNachbar<Refine) THEN
          zPL=zP(iz1)
          DO jx=jx0+1,jx1,IncrX
            DO jy=jy0+1,jy1,IncrY
              xPL=xP(jx)
              yPL=yP(jy)
              Vel(ibLoc)%wF(jx:jx+IncrX-1,jy:jy+IncrY-1,iz1)= &
                W(xPL,yPL,zPL,zH(jx,jy),Time)*Rho(xPL,yPL,zPL,zH(jx,jy),Time)
            END DO
          END DO
        END IF
      END IF
    END DO
    Vel(ibLoc)%uF=Vel(ibLoc)%uF*FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1) &
              /(FU(ix0:ix1,iy0+1:iy1,iz0+1:iz1)+Eps)
    Vel(ibLoc)%vF=Vel(ibLoc)%vF*FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1) &
              /(FV(ix0+1:ix1,iy0:iy1,iz0+1:iz1)+Eps)
    Vel(ibLoc)%wF=Vel(ibLoc)%wF*FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1) &
              /(FW(ix0+1:ix1,iy0+1:iy1,iz0:iz1)+Eps)

  END DO   ! ibLoc

!-----------------------------------------------------------
END SUBROUTINE RhoVelocityInit

SUBROUTINE QvsCompute(Scalar,Val,Vector,Time)
!
!==================================================
!----  Initialization of Qvs
!==================================================
!
  TYPE (ScalarCell_T), TARGET :: Scalar(:)
  TYPE (Vector4Cell_T) :: Vector(:)
  REAL(RealKind) :: Val
  EXTERNAL :: Val
  REAL(RealKind) :: Time

!---  local variables
  INTEGER :: ix, iy, iz
  REAL(RealKind) :: rh,thP,th
  REAL(RealKind), POINTER :: c(:,:,:,:) 
    
  DO ibLoc=1,nbLoc
    ib = LocGlob(ibLoc)
    CALL Set(Floor(ib))
    c=>Scalar(ibLoc)%c
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          rh =RhoCell(ibLoc)           %c(ix,iy,iz,1)
          thP=thProfG(ibLoc)           %c(ix,iy,iz,1)
          th =Vector (ibLoc)%Vec(thPos)%c(ix,iy,iz,1)
          c(ix,iy,iz,1) = Val(rh,thP,th) &
            *VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
        END DO
      END DO
    END DO
    c(ix0,:,:,:)=c(ix0+1,:,:,:)
    c(ix1+1,:,:,:)=c(ix1,:,:,:)
    c(:,iy0,:,:)=c(:,iy0+1,:,:)
    c(:,iy1+1,:,:)=c(:,iy1,:,:)
    c(:,:,iz0,:)=c(:,:,iz0+1,:)
    c(:,:,iz1+1,:)=c(:,:,iz1,:)

  END DO   ! ibLoc

!-----------------------------------------------------------
END SUBROUTINE QvsCompute

END MODULE Init_Mod
