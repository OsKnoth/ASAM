MODULE ReadProfile3D_Mod

  USE Kind_Mod

  IMPLICIT NONE 

  REAL(RealKind), ALLOCATABLE, PRIVATE :: xP(:),yP(:),zP(:)
  LOGICAL, PRIVATE :: Grid=.FALSE.
  INTEGER, PRIVATE :: nx,ny,nz

CONTAINS 

FUNCTION Profile(c,x,y,z)

  REAL(RealKind) :: Profile
  REAL(RealKind) :: c(:,:,:),x,y,z

  INTEGER :: ix,iy,iz
  

  DO ix=1,nx
    IF (x<=xP(ix+1)) THEN
      EXIT
    END IF
  END DO 
  DO iy=1,ny
    IF (y<=yP(iy+1)) THEN
      EXIT
    END IF
  END DO 
  DO iz=1,nz
    IF (z<=zP(iz+1)) THEN
      EXIT
    END IF
  END DO 
  Profile=c(ix,iy,iz)

END FUNCTION Profile

SUBROUTINE ReadGrid

  INTEGER :: ix,iy,iz,InputUnit
  CHARACTER(128) :: Line

  IF (.NOT.Grid) THEN
!   Find line
    InputUnit=1
    OPEN(UNIT=InputUnit,FILE='Profile',STATUS='OLD')
    DO
      READ(InputUnit,*,END=1) Line
      IF (INDEX(Line,'Grid')>0) THEN
        EXIT
      END IF
    END DO
    READ(InputUnit,*) nx,ny,nz
    ALLOCATE(xP(nx+1))
    ALLOCATE(yP(ny+1))
    ALLOCATE(zP(nz+1))
    DO ix=1,nx+1
      READ(InputUnit,*) xP(ix)
    END DO
    DO iy=1,ny+1
      READ(InputUnit,*) yP(iy)
    END DO
    DO iz=1,nz+1
      READ(InputUnit,*) zP(iz)
    END DO
    CLOSE(InputUnit)
1   CONTINUE
    Grid=.TRUE.
  END IF

END SUBROUTINE ReadGrid


SUBROUTINE ReadProfile(c,Name)

  REAl(RealKind), POINTER :: c(:,:,:)
  CHARACTER(*) :: Name

  INTEGER :: ix,iy,iz,InputUnit
  CHARACTER(128) :: Line

  CALL ReadGrid
! Find line
  InputUnit=1
  OPEN(UNIT=InputUnit,FILE='Profile',STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,TRIM(Name))>0) THEN
      EXIT
    END IF
  END DO
  READ(InputUnit,*) nx,ny,nz
  Allocate(c(nx,ny,nz))
  DO ix=1,nx
    DO iy=1,ny
      DO iz=1,nz
        READ(InputUnit,*) c(ix,iy,iz)
      END DO
    END DO
  END DO
  CLOSE(InputUnit)
1 CONTINUE

END SUBROUTINE ReadProfile

FUNCTION TkeStart(x,y,z,Time)

  REAL(RealKind) :: TkeStart

  REAL(RealKind) :: x,y,z,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:,:)
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Load) THEN
    CALL ReadProfile(cInt,'TkeProf')
    Load=.FALSE.
  END IF
  TkeStart=Profile(cInt,x,y,z)

END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,Time)

  REAL(RealKind) :: DisStart

  REAL(RealKind) :: x,y,z,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:,:)
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Load) THEN
    CALL ReadProfile(cInt,'DisProf')
    Load=.FALSE.
  END IF
  DisStart=Profile(cInt,x,y,z)

END FUNCTION DisStart

END MODULE ReadProfile3D_Mod

