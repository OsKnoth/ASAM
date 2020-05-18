MODULE Cosine_Mod

! USE ReadProfile_Mod, ONLY : 

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind), PARAMETER :: VelMax=1.0d0
  REAL(RealKind), PARAMETER :: N=1.0d-2
  REAL(RealKind), PARAMETER :: Rho0=1.0d0
  REAL(RealKind), PARAMETER :: th0=293.16d0
  REAL(RealKind), PARAMETER :: p0=101300.0d0 ! Pa
  REAL(RealKind), PARAMETER :: D0=1.0d0
  REAL(RealKind), PARAMETER :: qvr=0.95d0
  REAL(RealKind), PRIVATE :: uMax=0.0d0,vMax=0.0d0
  REAL(RealKind), PRIVATE :: TkeMax=0.0d0,DisMax=0.0d0

CONTAINS 

SUBROUTINE InputExample

  WRITE(*,*) 'Eingabe uMax'
  READ(*,*) uMax

END SUBROUTINE InputExample

FUNCTION RhoFun(x,y,z)

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z


  RhoFun=p0/(Rd*th0)*(1.0d0-kappa*Grav*z/(Rd*th0))**(1.0d0/kappa-1.0d0)
  RhoFun=1.0d0

END FUNCTION RhoFun

FUNCTION ThProfFun(x,y,z)

  REAL(RealKind) :: ThProfFun

  REAL(RealKind) :: x,y,z

  ThProfFun=th0*RhoFun(x,y,z)

END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z)

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z

  REAL(RealKind) :: RhoLoc,ThLoc

  RhoLoc=RhoFun(x,y,z)
  ThLoc=thProfFun(x,y,z)
  QvProfFun=qvr*qvs(RhoLoc,ThLoc,ThLoc)

END FUNCTION QvProfFun

FUNCTION UStart(x,y,z,Time)

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,Time
  REAL(RealKind) :: yL,yR

  REAL(RealKind) :: U0

  UStart=uMax*RhoFun(x,y,z)

END FUNCTION UStart

FUNCTION UStartE(x,y,z,Time)

  REAL(RealKind) :: UStartE

  REAL(RealKind) :: x,y,z,Time

  UStartE=UStart(x,y,z,Time)

END FUNCTION UStartE

FUNCTION VStart(x,y,z,Time)

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,Time

  VStart=vMax

END FUNCTION VStart

FUNCTION VStartE(x,y,z,Time)

  REAL(RealKind) :: VStartE

  REAL(RealKind) :: x,y,z,Time

  VStartE=VStart(x,y,z,Time)

END FUNCTION VStartE

FUNCTION WStart(x,y,z,Time)

  REAL(RealKind) :: WStart

  REAL(RealKind) :: x,y,z,Time

  WStart=0.0d0*RhoFun(x,y,z)

END FUNCTION WStart

FUNCTION ThStart(x,y,z,Time)

  REAL(RealKind) :: ThStart

  REAL(RealKind) :: x,y,z,Time

  ThStart=ThProfFun(x,y,z)

END FUNCTION ThStart

FUNCTION TkeStart(x,y,z,Time)

  REAL(RealKind) :: TkeStart

  REAL(RealKind) :: x,y,z,Time

  TkeStart=TkeMax
  
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,Time)

  REAL(RealKind) :: DisStart

  REAL(RealKind) :: x,y,z,Time

  DisStart=DisMax
  
END FUNCTION DisStart

FUNCTION QvStart(x,y,z,Time)

  REAL(RealKind) :: QvStart

  REAL(RealKind) :: x,y,z,Time

  QvStart=QvProfFun(x,y,z)
  
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,Time)

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,Time

  QcStart=0.0d0

END FUNCTION QcStart

FUNCTION QrStart(x,y,z,Time)

  REAL(RealKind) :: QrStart

  REAL(RealKind) :: x,y,z,Time

  QrStart=0.d0

END FUNCTION QrStart

FUNCTION DStart(x,y,z,Time)

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,Time

  DStart=D0*RhoFun(x,y,z)
  
END FUNCTION DStart

FUNCTION DummyStart(x,y,z,Time)

  REAL(RealKind) :: DummyStart

  REAL(RealKind) :: x,y,z,Time

  DummyStart=0.0d0

END FUNCTION DummyStart

END MODULE Cosine_Mod

