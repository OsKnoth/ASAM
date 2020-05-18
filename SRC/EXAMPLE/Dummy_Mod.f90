MODULE Dummy_Mod

  USE Kind_Mod

  IMPLICIT NONE 

CONTAINS 
FUNCTION RhoProf(x,y,z,zHeight,Time)
  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoProf=0.0d0
END FUNCTION RhoProf
FUNCTION ThProfFun(x,y,z,zHeight,Time)
  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThProfFun=0.0d0
END FUNCTION ThProfFun
FUNCTION RhoFun(x,y,z,zHeight,Time)
  REAL(RealKind) :: RhoFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoFun=0.0d0
END FUNCTION RhoFun
FUNCTION HeightFun(x,y,z,zHeight,Time)
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  HeightFun=0.0d0
END FUNCTION HeightFun
FUNCTION ThStart(x,y,z,zHeight,Time)
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThStart=0.0d0
END FUNCTION ThStart
FUNCTION QvStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  qvStart=0.0d0
END FUNCTION QvStart
FUNCTION QvProfFun(x,y,z,zHeight,Time)
  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  qvProfFun=0.0d0
END FUNCTION QvProfFun
FUNCTION QcStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: QcStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  qcStart=0.0d0
END FUNCTION QcStart
FUNCTION QiStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=0.0d0
END FUNCTION QiStart
FUNCTION ForceU(x,y,z,zHeight,Time)
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=0.0d0
END FUNCTION ForceU
FUNCTION ForceV(x,y,z,zHeight,Time)
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=0.0d0
END FUNCTION ForceV
FUNCTION ForceW(x,y,z,zHeight,Time)
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=0.0d0
END FUNCTION ForceW
FUNCTION ForceRho(x,y,z,zHeight,Time)
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceRho=0.0d0
END FUNCTION ForceRho
FUNCTION UStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStart=0.0d0
END FUNCTION UStart
FUNCTION VStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=0.0d0
END FUNCTION VStart
FUNCTION DStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=0.0d0
END FUNCTION DStart
FUNCTION UStartE(x,y,z,zHeight,Time)
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=0.0d0
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=0.0d0
END FUNCTION VStartE
FUNCTION WStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart
FUNCTION TkeStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=0.0d0
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=0.0d0
END FUNCTION DisStart
FUNCTION TkeHStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=0.0d0
END FUNCTION TkeHStart
FUNCTION TkeVStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=0.0d0
END FUNCTION TkeVStart
FUNCTION LenStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=0.0d0
END FUNCTION LenStart
FUNCTION QrStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart
FUNCTION DummyStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart
FUNCTION TStart(x,y,z,zHeight,Time)
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=0.0d0
END FUNCTION TStart
END MODULE Dummy_Mod
