MODULE ReadProfile_Mod

  USE Physics_Mod

  IMPLICIT NONE 

  INTERFACE ProfileEqual
    MODULE PROCEDURE ProfileEqual,ProfileEqual2
  END INTERFACE

CONTAINS 

FUNCTION ProfileEqual(c,z)

  REAL(RealKind) :: ProfileEqual
  REAL(RealKind) :: c(:,:),z

  INTEGER :: i,n
  REAL(RealKind) :: zLoc
  
  n=SIZE(c,1)

  zLoc=MIN(z,c(n,1))
  IF (zLoc<=c(1,1)) THEN
    ProfileEqual=c(1,2)
  ELSE
    DO i=2,n
      IF (zLoc<=c(i,1)) THEN
        ProfileEqual=Interp(zLoc,c(i-1,1),c(i-1,2),c(i,1),c(i,2))
        EXIT
      END IF
    END DO
  END IF

END FUNCTION ProfileEqual

FUNCTION ProfileEqual2(c,zC,z)

  REAL(RealKind) :: ProfileEqual2
  REAL(RealKind) :: c(:),zC(:),z

  INTEGER :: i,n
  REAL(RealKind) :: zLoc
  
  n=SIZE(c)

  zLoc=MIN(z,zC(n))
  IF (zLoc<=zC(1)) THEN
    ProfileEqual2=c(1)
  ELSE
    DO i=2,n
      IF (zLoc<=zC(i)) THEN
        ProfileEqual2=Interp(zLoc,zC(i-1),c(i-1),zC(i),c(i))
        EXIT
      END IF
    END DO
  END IF

END FUNCTION ProfileEqual2

FUNCTION ProfileStretch(c,z,zSurf)

  REAL(RealKind) :: ProfileStretch
  REAL(RealKind) :: c(:,:),z,zSurf

  INTEGER :: i,n
  REAL(RealKind) :: zLoc,zEnd
  REAL(RealKind) :: zeta(SIZE(c,1))
  REAL(RealKind) :: alpha,beta  ! zeta=alpha*z+beta
 
  n=SIZE(c,1)
  zeta(n)=c(n,1)
  zeta(1)=zSurf
  beta=(zeta(1)-zeta(n))/(c(1,1)-c(n,1))
  alpha=zeta(1)-beta*c(1,1)
  DO i=2,n-1
    zeta(i)=alpha+beta*c(i,1)
  END DO

  zLoc=MIN(z,zeta(n))
  DO i=2,n
    IF (zLoc<=zeta(i)) THEN
      ProfileStretch=Interp(zLoc,zeta(i-1),c(i-1,2),zeta(i),c(i,2))
      EXIT
    END IF
  END DO

END FUNCTION ProfileStretch

FUNCTION Interp(zM,zL,cL,zR,cR)
  REAL(RealKind) :: Interp,zM,zL,zR,cL,cR
  Interp=((zM-zL)*cR+(zR-zM)*cL)/(zR-zL)
END FUNCTION Interp


SUBROUTINE ReadProfile(c,Type,Name)

  REAl(RealKind), POINTER :: c(:,:)
  CHARACTER(*) :: Type
  CHARACTER(*) :: Name

  INTEGER :: i,n,nc,InputUnit
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
  READ(InputUnit,*) n,nc
  READ(InputUnit,*) Type
  Allocate(c(n,nc))
  DO i=1,n
    READ(InputUnit,*) c(i,:)
  END DO
1 CONTINUE
  CLOSE(InputUnit)

END SUBROUTINE ReadProfile

SUBROUTINE ReadRProfile(c,Name)

  REAl(RealKind), POINTER :: c(:,:)
  CHARACTER(*) :: Name

  INTEGER :: i,n,nc,InputUnit
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
  READ(InputUnit,*) n,nc
  Allocate(c(n,nc))
  DO i=n,1,-1
    READ(InputUnit,*) c(i,:)
  END DO
  CLOSE(InputUnit)
1 CONTINUE

END SUBROUTINE ReadRProfile

END MODULE ReadProfile_Mod

