MODULE Floor_Mod

  USE Kind_Mod
  USE Domain_Mod
  USE Parallel_Mod

  IMPLICIT NONE

  INTEGER :: nb
  INTEGER :: ib,ibLoc
  TYPE (Domain_T) :: Domain
  TYPE (Domain_T), POINTER :: Floor(:)

  INTERFACE Allocate
    MODULE PROCEDURE FloorAllocate
  END INTERFACE  
CONTAINS

SUBROUTINE FloorAllocate(Floor)

   TYPE (Domain_T), POINTER :: Floor(:)

   DO ibLoc=1,nbLoc
     ib=LocGlob(ibLoc)
     CALL Allocate(Floor(ib))
   END DO

END SUBROUTINE FloorAllocate

SUBROUTINE LocatePoint(x,y,z,ibLocP,ixP,iyP,izP)
  REAL(RealKind) :: x,y,z
  INTEGER :: ibLocP,ib1,ixP,iyP,izP

  INTEGER :: ix,iy,iz

  ibLocP=0
  SB:DO ib1=1,nbLoc
    ib=LocGlob(ib1)
    IF (Floor(ib)%x0<=x.AND.x<=Floor(ib)%x1.AND. &
        Floor(ib)%y0<=y.AND.y<=Floor(ib)%y1.AND. &
        Floor(ib)%z0<=z.AND.z<=Floor(ib)%z1) THEN
     !ibLocP=ibLoc
      ibLocP=ib1
      SX:DO ix=Floor(ib)%ix0+1,Floor(ib)%ix1
          IF (Floor(ib)%xP(ix-1)<=x.AND.x<=Floor(ib)%xP(ix)) THEN 
            ixP=ix
            EXIT SX
          END IF
      END DO SX
      SY:DO iy=Floor(ib)%iy0+1,Floor(ib)%iy1
          IF (Floor(ib)%yP(iy-1)<=y.AND.y<=Floor(ib)%yP(iy)) THEN 
            iyP=iy
            EXIT SY
          END IF
      END DO SY
      SZ:DO iz=Floor(ib)%iz0+1,Floor(ib)%iz1
          IF (Floor(ib)%zP(iz-1)<=z.AND.z<=Floor(ib)%zP(iz)) THEN 
            izP=iz
            EXIT SZ
          END IF
      END DO SZ
      EXIT SB
    END IF
  END DO SB

END SUBROUTINE LocatePoint

END MODULE Floor_Mod





