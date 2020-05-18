MODULE Read3D_Mod

  USE Kind_Mod

  IMPLICIT NONE 

  INTEGER, PRIVATE :: nx,ny,nz

CONTAINS 

FUNCTION Set3D(c,x,y,z)

  REAL(RealKind) :: Set3D
  REAL(RealKind) :: c(:,:,:),x,y,z

  INTEGER :: i,n
  REAL(RealKind) :: zLoc
  
  n=SIZE(c,1)

  zLoc=MIN(z,c(n,1))
  IF (zLoc<=c(1,1)) THEN
    3D=c(1,2)
  ELSE
    DO i=2,n
      IF (zLoc<=c(i,1)) THEN
        3D=Int(zLoc,c(i-1,1),c(i-1,2),c(i,1),c(i,2))
        EXIT
      END IF
    END DO
  END IF

CONTAINS
  FUNCTION Int(zM,zL,cL,zR,cR)
    REAL(RealKind) :: Int,zM,zL,zR,cL,cR
    Int=((zM-zL)*cR+(zR-zM)*cL)/(zR-zL)
   END FUNCTION Int
END FUNCTION 3D

SUBROUTINE Read3D(c,Name)

  REAl(8), POINTER :: c(:,:,:)
  CHARACTER(*) :: Name

  INTEGER :: i,n,InputUnit
  CHARACTER(128) :: Line

! Find line
  InputUnit=1
  OPEN(UNIT=InputUnit,FILE='3D',STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,TRIM(Name))>0) THEN
      EXIT
    END IF
  END DO
  READ(InputUnit,*) nx,ny,nz
  Allocate(c(nx,ny,nz))
  DO i=1,n
    READ(InputUnit,*) c(i,:)
  END DO
  CLOSE(InputUnit)
1 CONTINUE

END SUBROUTINE Read3D

FUNCTION UStart(x,y,z,Time)

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,Time
  REAL(RealKind), POINTER, SAVE :: uInt(:,:,:)) 
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Load) THEN
    CALL Read3D(uInt,'uProf')
    Load=.FALSE.
  END IF
  UStart=Set3D(uInt,x,y,z)

END FUNCTION UStart

FUNCTION VStart(x,y,z,Time)

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,Time
  REAL(RealKind), POINTER, SAVE :: vInt(:,:)
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Load) THEN
    CALL Read3D(vInt,'vProf')
    Load=.FALSE.
  END IF
  VStart=3D(vInt,z)

END FUNCTION VStart

FUNCTION ThStart(x,y,z,Time)

  REAL(RealKind) :: ThStart

  REAL(RealKind) :: x,y,z,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Load) THEN
    CALL Read3D(cInt,'ThProf')
    Load=.FALSE.
  END IF
  ThStart=3D(cInt,z)

END FUNCTION ThStart

FUNCTION TkeStart(x,y,z,Time)

  REAL(RealKind) :: TkeStart

  REAL(RealKind) :: x,y,z,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Load) THEN
    CALL Read3D(cInt,'TkeProf')
    Load=.FALSE.
  END IF
  TkeStart=3D(cInt,z)

END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,Time)

  REAL(RealKind) :: DisStart

  REAL(RealKind) :: x,y,z,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  LOGICAL, SAVE :: Load=.TRUE.

  IF (Load) THEN
    CALL Read3D(cInt,'DisProf')
    Load=.FALSE.
  END IF
  DisStart=3D(cInt,z)

END FUNCTION DisStart

END MODULE Read3D_Mod

