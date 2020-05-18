MODULE ReadProfileRad_Mod

  USE Physics_Mod

  IMPLICIT NONE 

CONTAINS 

FUNCTION ProfileRad(r,z,c,rP,zP)

  REAL(RealKind) :: ProfileRad
  REAL(RealKind) :: r,z,c(:,:),rP(:),zP(:)

  INTEGER :: ir,iz,nr,nz
  REAL(RealKind) :: zLoc,rLoc
  
  nr=SIZE(rP)
  nz=SIZE(zP)

  zLoc=MAX(MIN(z,zP(nz)),zP(1))
  rLoc=MAX(MIN(r,rP(nr)),rP(1))
  DO iz=2,nz
    IF (zLoc<=zP(iz)) THEN
      EXIT
    END IF
  END DO
  DO ir=2,nr
    IF (rLoc<=rP(ir)) THEN
      EXIT
    END IF
  END DO
  ProfileRad=Bilinear(rLoc,zLoc,c(ir-1,iz-1),c(ir,iz-1),c(ir-1,iz),c(ir,iz) &
                       ,rP(ir-1),rP(ir),zP(iz-1),zP(iz))


END FUNCTION ProfileRad


FUNCTION Bilinear(x,y,c11,c21,c12,c22 &
                 ,x1,x2,y1,y2)
  REAL(RealKind) :: Bilinear,x,y
  REAL(RealKind) :: c11,c21,c12,c22
  REAL(RealKind) :: x1,x2,y1,y2
  Bilinear=(c11*(x2-x)*(y2-y) &
           +c21*(x-x1)*(y2-y) &
           +c12*(x2-x)*(y-y1) &
           +c22*(x-x1)*(y-y1))&
           /((x2-x1)*(y2-y1)) 
END FUNCTION Bilinear


SUBROUTINE ReadProfileRad(zP,rP,c,Type,Name)

  REAl(RealKind), POINTER :: zP(:)
  REAl(RealKind), POINTER :: rP(:)
  REAl(RealKind), POINTER :: c(:,:)
  CHARACTER(*) :: Type
  CHARACTER(*) :: Name

  INTEGER :: i,k,nz,nr,InputUnit
  CHARACTER(128) :: Line

! Find line
  InputUnit=1
  OPEN(UNIT=InputUnit,FILE=Profile,STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,TRIM(Name))>0) THEN
      EXIT
    END IF
  END DO
  READ(InputUnit,*) nz,nr
  ALLOCATE(zP(nz))
  DO k=1,nz
    READ(InputUnit,*) zP(k)
  END DO
  ALLOCATE(rP(nr))
  DO i=1,nr
    READ(InputUnit,*) rP(i)
  END DO
  ALLOCATE(c(nr,nz))
  DO i=1,nr
    DO k=1,nz
      READ(InputUnit,*) c(i,k)
    END DO
  END DO
  CLOSE(InputUnit)
1 CONTINUE

END SUBROUTINE ReadProfileRad

END MODULE ReadProfileRad_Mod
