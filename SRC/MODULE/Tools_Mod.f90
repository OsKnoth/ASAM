MODULE Tools_Mod
  
  USE Domain_Mod
  USE Floor_Mod
  USE DataType_Mod
  USE Control_Mod,  ONLY: PrintNameLists

  IMPLICIT NONE

CONTAINS 

SUBROUTINE MinVol

  INTEGER :: ixVol,iyVol,izVol
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: VolMin ! ,VolMax

  VolMin=1.0e20_RealKind
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          IF (VolC(ix,iy,iz)>0.0d0) THEN
            IF (VolC(ix,iy,iz)<VolMin) THEN
              ixVol=ix
              iyVol=iy
              izVol=iz
              VolMin=VolC(ix,iy,iz)
            END IF
          END IF
        END DO
      END DO
    END DO
  END DO
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    DO ix=ix0+1,ix1
      DO iy=iy0+1,iy1
        DO iz=iz0+1,iz1
          IF (ix==ixVol.AND.iy==iyVol.AND.iz==izVol .AND. PrintNameLists) THEN
            WRITE(*,*) 'Volume'
            WRITE(*,*) ixVol,iyVol,izVol
            WRITE(*,*) VolMin,VolMin/(dx(ixVol)*dy(iyVol)*dz(izVol)) ! ,VolMax,VolMin/VolMax
            WRITE(*,*) FU(ixVol-1,iyVol,izVol),FU(ixVol,iyVol,izVol)
            WRITE(*,*) FV(ixVol,iyVol-1,izVol),FV(ixVol,iyVol,izVol)
            WRITE(*,*) FW(ixVol,iyVol,izVol-1),FW(ixVol,iyVol,izVol)
            WRITE(*,*) 1.0d0/6.0d0*SQRT(8.0d0*MAX(FU(ixVol-1,iyVol,izVol),FU(ixVol,iyVol,izVol)) &
                                    *MAX(FV(ixVol,iyVol-1,izVol),FV(ixVol,iyVol,izVol)) &
                                    *MAX(FW(ixVol,iyVol,izVol-1),FW(ixVol,iyVol,izVol)))
          END IF
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE MinVol

SUBROUTINE WriteVec4(C,Name,ix0,ix1,iy0,iy1,iz0,iz1)

  INTEGER :: ix0,ix1,iy0,iy1,iz0,iz1
  REAL(RealKind):: C(ix0:ix1,iy0:iy1,iz0:iz1)
  CHARACTER(LEN=*) :: Name
  
  INTEGER :: ix,iy,iz

  WRITE(*,*) Name
  DO ix=ix0,ix1
    DO iy=iy0,iy1
      DO iz=iz0,iz1
        WRITE(*,*) ix,iy,iz,c(ix,iy,iz)
      END DO
    END DO
  END DO

END SUBROUTINE WriteVec4

SUBROUTINE WriteSomeC(VectorCell,Name,ixP,iyP,izP,ibP)

  TYPE(Vector4Cell_T), POINTER :: VectorCell(:)
  CHARACTER(LEN=*) :: Name
  INTEGER, OPTIONAL :: ixP,iyP,izP,ibP
  
  INTEGER :: ix,iy,iz

  IF (PRESENT(ixP)) THEN
    WRITE(*,*) Name,ixP,iyP,izP
    WRITE(*,*) VectorCell(ibP)%Vec(1)%c(ixP-1:ixP+1,iyP,izP,1)
    WRITE(*,*) VectorCell(ibP)%Vec(2)%c(ixP-1:ixP+1,iyP,izP,1)
    DO ix=ixp-5,ixp+5
      WRITE(*,*) ix,VectorCell(ibP)%Vec(5)%c(ix,iyP,izP:izp+3,1)
      WRITE(*,*) ix,VectorCell(ibP)%Vec(6)%c(ix,iyP,izP:izp+3,1)
    END DO
  ELSE
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          DO iz=iz0+1,iz1
            IF (VolC(ix,iy,iz)>Zero.AND.VolC(ix,iy,iz)<=5.d-2*dx(ix)*dy(iy)*dz(iz)) THEN
              WRITE(*,*) Name,ix,iy,iz 
              WRITE(*,*) VectorCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
              WRITE(*,*) VectorCell(ibLoc)%Vec(2)%c(ix,iy,iz,1)
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE WriteSomeC

SUBROUTINE WriteSomeS(ScalarCell,Name,ixP,iyP,izP,ibP)

  TYPE(ScalarCell_T), POINTER :: ScalarCell(:)
  CHARACTER(LEN=*) :: Name
  INTEGER, OPTIONAL :: ixP,iyP,izP,ibP

  INTEGER :: ix,iy,iz

  IF (PRESENT(ixP)) THEN
    WRITE(*,*) Name,ixP,iyP,izP
    WRITE(*,*) ScalarCell(ibP)%c(ixP-1:ixP+1,iyP,izP,1)
    DO ix=ixp-5,ixp+5
      WRITE(*,*) ix,ScalarCell(ibP)%c(ix,iyP,izP:izp+3,1)
      WRITE(*,*) VolC(ix,iyp,izP:izp+3)
      WRITE(*,*) FW(ix,iyp,izP-1:izp+3)
    END DO
  ELSE
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          DO iz=iz0+1,iz1
            IF (VolC(ix,iy,iz)>Zero.AND.VolC(ix,iy,iz)<=5.d-2*dx(ix)*dy(iy)*dz(iz)) THEN
              WRITE(*,*) Name,ix,iy,iz
              WRITE(*,*) ScalarCell(ibLoc)%c(ix,iy,iz,1)
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE WriteSomeS

SUBROUTINE WriteSomeF(VelocityFace,Name,ixP,iyP,izP,ibP)

  TYPE(VelocityFace_T), POINTER :: VelocityFace(:)
  CHARACTER(LEN=*) :: Name
  INTEGER, OPTIONAL :: ixP,iyP,izP,ibP

  INTEGER :: ix,iy,iz

  IF (PRESENT(ixP)) THEN
    WRITE(*,*) Name,ixP,iyP,izP
    WRITE(*,*) VelocityFace(ibP)%uF(ixP-1,iyP,izP),VelocityFace(ibLoc)%uF(ixP,iyP,izP)
    WRITE(*,*) VelocityFace(ibP)%wF(ixP,iyP,izP-1),VelocityFace(ibLoc)%wF(ixP,iyP,izP)
  ELSE
    VelocityFaceAct=>VelocityFace
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      DO ix=ix0+1,ix1
        DO iy=iy0+1,iy1
          DO iz=iz0+1,iz1
            IF (VolC(ix,iy,iz)>Zero.AND.VolC(ix,iy,iz)<=5.d-2*dx(ix)*dy(iy)*dz(iz)) THEN
              WRITE(*,*) Name,ix,iy,iz
              WRITE(*,*) VelocityFace(ibLoc)%uF(ix-1,iy,iz),VelocityFace(ibLoc)%uF(ix,iy,iz)
              WRITE(*,*) VelocityFace(ibLoc)%wF(ix,iy,iz-1),VelocityFace(ibLoc)%wF(ix,iy,iz)
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF
END SUBROUTINE WriteSomeF

END MODULE Tools_Mod
