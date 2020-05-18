MODULE Vortex_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: x,y,z
  REAL(RealKind) :: xr,yr,zr
  REAL(RealKind) :: Xi0=2.34d-3      !in s^-1
  REAL(RealKind) :: b=53.5d+3        !in m
!  REAL(RealKind) :: v0               !=v0(r,z) vertical wind field 
  REAL(RealKind) :: alpha=2.5d0        !vertical decay rate of the wind in the baroclinic region
  REAL(RealKind) :: Lz=6.0d+3          !in m, indicates the depth of the barotropic part of the vortex
!  REAL(RealKind) :: Thetai           !potential temperature perturbation
  REAL(RealKind) :: rb=5.0d+3         !in m
  REAL(RealKind) :: zb=6.0d+3          !in m
  REAL(RealKind) :: Sigmar=15.0d+3     !in m
  REAL(RealKind) :: Sigmaz=3000d0      !in m
  REAL(RealKind) :: A=0.5d0            !in K
!  REAL(RealKind) :: Lambda=2*Pi 
  
  NAMELIST /Example/ x &
                    ,y &
                    ,z &
                    ,Xi0 &
                    ,b &
                    ,alpha &
                    ,Lz &
                    ,rb &
                    ,zb &
                    ,Sigmar &
                    ,Sigmaz &
                    ,A &
                    ,Lambda 


END MODULE Vortex_Mod

SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile


SUBROUTINE InputExample(FileName)
  USE Vortex_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&Example')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=Example)
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(UNIT=InputUnit)
END SUBROUTINE InputExample

FUNCTION UStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: Vr               !azimuthal wind field
  REAL(RealKind) :: r,Lambda
!  Pi=4.0d0*ATAN(1.0d0)
  r=SQRT((x*x)+(y*y))
  Vr=(Xi0*(b*b)*(1.0d0-EXP(-(r/b)*(r/b))))
  Lambda=ATAN(y/x)
  UStart=Vr*SIN(Lambda)

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: Vr               !azimuthal wind field
  REAL(RealKind) :: r,Lambda
  r=SQRT((x*x)+(y*y))
  Vr=(Xi0*(b*b)*(1.0d0-EXP(-(r/b)*(r/b))))
  Lambda=ATAN(y/x)
  VStart=Vr*COS(Lambda)

END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  ThProfFun=ThInit
END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad,pLoc

  pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
  Rad=SQRT(((x-x0)/xr)**2 &
          +((z-z0)/zr)**2 &
            ) 
  ThStart=ThInit
  IF (Rad<1.0d0) THEN
    ThStart=ThStart+DeltaTh*COS(Pi*Rad/2.0d0)**2 
  END IF

END FUNCTION ThStart

FUNCTION QvProfFun(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QvStart=Zero

END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad

  Rad=SQRT(((x-x0Qc)/xr)**2 &
          +((y-y0Qc)/yr)**2 &
          +((z-z0Qc)/zr)**2 &
            ) 
  QcStart=Zero
  IF (Rad<1.0d0) THEN
    QcStart=QcStart+DeltaQc*(COS(Pi*Rad)+1.0d0)/2.0d0 
  END IF


END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=Zero

END FUNCTION QiStart


FUNCTION DStart(x,y,z,zHeight,Time)

  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=DMax

END FUNCTION DStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,ThLoc,Rad

  pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
  Rad=SQRT(((x-x0)/xr)**2 &
          +((z-z0)/zr)**2 &
            )
  ThLoc=ThInit
  IF (Rad<1.0d0) THEN
    ThLoc=ThLoc+DeltaTh*COS(Pi*Rad/2.0d0)**2 
  END IF
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Vortex_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc

  pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThInit)
  RhoProf=Zero

END FUNCTION RhoProf


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Vortex_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Vortex_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart


FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(cpD/Rd)
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Vortex_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart





