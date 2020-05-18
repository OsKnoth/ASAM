MODULE Durran_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 
  REAL(RealKind) :: NL=1.0d-2
  REAL(RealKind) :: NU=1.0d-2
  REAL(RealKind) :: H=1.0d-2
  REAL(RealKind) :: RhoL=1.15d0
  REAL(RealKind) :: RhoU=0.0d0
  REAL(RealKind) :: th0=293.16d0
  REAL(RealKind) :: D0=10.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind) :: uMax=10.0d0,vMax=0.0d0
  REAL(RealKind) :: TkeMax=1.0d-2,DisMax=1.0d-4
  REAL(RealKind) :: TkeHMax=1.0d-2,TkeVMax=1.0d-2,LenMax=1.0d3
  REAL(RealKind) :: uStar=1.00d0, z0=0.1d0                          !!!FILAUS

  NAMELIST /Example/    &
                    uMax , &
                    TkeMax , &
                    NL,    &
                    NU,    &
                    th0,   &
                    D0,    &
                    uStar, &
                    z0, &
                    H

CONTAINS


FUNCTION ThetaLoc(z)

  REAL(RealKind) :: ThetaLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: SL,SU

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    ThetaLoc=th0*EXP(SL)
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
    ThetaLoc=th0*EXP(SL+SU)
  END IF

END FUNCTION ThetaLoc

FUNCTION PressLoc(z)

  REAL(RealKind) :: PressLoc
  REAL(RealKind) :: z

  REAL(RealKind) :: SL,SU
  REAL(RealKind) :: t1,t2,pH
  REAL(RealKind) :: s1,s2

  IF (z<=H) THEN
    SL=NL*NL*z/Grav
    t1=Grav*Grav*kappa/(Rd*th0*NL*NL)
    t2=Grav*Grav*kappa/(Rd*th0*NL*NL)*EXP(-SL)
    PressLoc=p0 &
            *(t2-t1+1.0d0)**(1.0d0/kappa)
           
  ELSE
    SL=NL*NL*H/Grav
    SU=NU*NU*(z-H)/Grav
    t1=Grav*Grav*kappa/(Rd*th0*NL*NL)
    t2=Grav*Grav*kappa/(Rd*th0*NL*NL)*EXP(-SL)
    pH=p0 &
       *(t2-t1+1.0d0)**(1.0d0/kappa)
    s1=Grav*Grav*kappa/(Rd*th0*NU*NU)*EXP(-SL)
    s2=Grav*Grav*kappa/(Rd*th0*NU*NU)*EXP(-SL-SU)
    PressLoc=p0 &
       *(t2-t1+s2-s1+1.0d0)**(1.0d0/kappa)
  END IF

END FUNCTION PressLoc

END MODULE Durran_Mod

SUBROUTINE PerturbProfile(VecC)

  USE DataType_Mod
  IMPLICIT NONE
  TYPE(Vector4Cell_T) :: VecC(:)

END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE Durran_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

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

FUNCTION RhoFun(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  pLoc=PressLoc(z)
  ThLoc=ThetaLoc(z)
  RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
!  RhoFun=1.0d0         !!!FILAUS

END FUNCTION RhoFun

FUNCTION RhoProf(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc

  pLoc=PressLoc(z)
  ThLoc=ThetaLoc(z)
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
!  RhoProf=1.0d0         !!!FILAUS

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc
  ThLoc=ThetaLoc(z)
  ThProfFun=ThLoc
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvProfFun=Zero
END FUNCTION QvProfFun

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: yL,yR
  REAL(RealKind) :: U0
!  UStart=uStar/karm*LOG(z/z0+1)
  UStart=uMax
END FUNCTION UStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=vMax
END FUNCTION VStart

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThStart=ThetaLoc(z)
END FUNCTION ThStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThLoc,pLoc,Fpot,RhoLoc
  ThLoc  =ThetaLoc(z)
  pLoc   =PressLoc(z)
  RhoLoc =pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  Fpot   = (Rd*RhoLoc*ThLoc/p0)**(-kappa/(One-kappa))
  TStart =ThLoc/Fpot
END FUNCTION TStart
 
!FUNCTION RhoStart(x,y,z,zHeight,Time)
!  USE Durran_Mod
!  USE Rho_Mod
!  USE RhoProf_Mod
!  IMPLICIT NONE
!  REAL(RealKind) :: RhoStart
!  REAL(RealKind) :: x,y,z,zHeight,Time
!  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
!END FUNCTION RhoStart


FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
!  TkeStart=uStar**2.0d0/SQRT(Cmy0)
  TkeStart=TkeMax
END FUNCTION TkeStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
!  TkeHStart=(2.0d0/3.0d0)*uStar**2.0d0/SQRT(Cmy0)
!  TkeHStart=TkeHMax
  TkeHStart=One/(0.43d0+Two*Cmy0**0.75d0)*((Two/Three)*0.43d0*uStar**Two/SQRT(Cmy0) &
                                           +Two*uStar**Two*Cmy0**0.25d0) 
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
!  TkeVStart=(1.0d0/3.0d0)*uStar**2.0d0/SQRT(Cmy0)
!  TkeVStart=TkeVMax
  TkeVStart=(One/Three)*0.43d0*uStar**Two/(SQRT(Cmy0)*(0.43d0+Two*Cmy0**0.75d0))
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=karm*(z+z0)
END FUNCTION LenStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
!  DisStart=uStar**3.0d0/(karm*(z+z0))
  DisStart=DisMax
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QcStart=0.0d0
END FUNCTION QcStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=D0
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: S,RadLoc,ThBack
  PreStart=PressLoc(z)
END FUNCTION PreStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Durran_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.d0
END FUNCTION QiStart


