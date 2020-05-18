MODULE Rho_Mod
  INTERFACE 
  FUNCTION RhoFun(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: RhoFun
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION RhoFun
  END INTERFACE
END MODULE Rho_Mod

MODULE RhoProf_Mod
  INTERFACE 
  FUNCTION RhoProf(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: RhoProf
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION RhoProf
  END INTERFACE
END MODULE RhoProf_Mod

MODULE RhoStart_Mod
  INTERFACE 
  FUNCTION RhoStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE RhoProf_Mod
    USE Rho_Mod
    REAL(RealKind) :: RhoStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION RhoStart
  END INTERFACE
END MODULE RhoStart_Mod

MODULE ThProf_Mod
  INTERFACE 
  FUNCTION ThProfFun(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    REAL(RealKind) :: ThProfFun
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION ThProfFun
  END INTERFACE
END MODULE ThProf_Mod

MODULE PreProf_Mod
  INTERFACE 
  FUNCTION PreProfFun(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    REAL(RealKind) :: PreProfFun
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION PreProfFun
  END INTERFACE
END MODULE PreProf_Mod

MODULE QvProf_Mod
  INTERFACE 
  FUNCTION QvProfFun(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    USE ThProf_Mod
    REAL(RealKind) :: QvProfFun
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION QvProfFun
  END INTERFACE
END MODULE QvProf_Mod

MODULE QcProf_Mod
  INTERFACE 
  FUNCTION QcProfFun(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    USE ThProf_Mod
    REAL(RealKind) :: QcProfFun
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION QcProfFun
  END INTERFACE
END MODULE QcProf_Mod

MODULE NvProf_Mod
  INTERFACE 
  FUNCTION NvProfFun(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    USE ThProf_Mod
    REAL(RealKind) :: NvProfFun
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION NvProfFun
  END INTERFACE
END MODULE NvProf_Mod

MODULE NcProf_Mod
  INTERFACE 
  FUNCTION NcProfFun(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    USE ThProf_Mod
    REAL(RealKind) :: NcProfFun
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION NcProfFun
  END INTERFACE
END MODULE NcProf_Mod

MODULE UVW_Mod
  INTERFACE 
  FUNCTION UStart(x,y,z,zH,Time)
    USE Kind_Mod
    USE Rho_Mod
    REAL(RealKind) :: UStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION UStart
  END INTERFACE

  INTERFACE 
  FUNCTION VStart(x,y,z,zH,Time)
    USE Kind_Mod
    USE Rho_Mod
    REAL(RealKind) :: VStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION VStart
  END INTERFACE

  INTERFACE 
  FUNCTION WStart(x,y,z,zH,Time)
    USE Kind_Mod
    USE Rho_Mod
    REAL(RealKind) :: WStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION WStart
  END INTERFACE
END MODULE UVW_Mod

MODULE Start_Mod
  INTERFACE 
  FUNCTION ForceU(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: ForceU
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION ForceU
  END INTERFACE

  INTERFACE 
  FUNCTION ForceV(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: ForceV
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION ForceV
  END INTERFACE

  INTERFACE 
  FUNCTION ForceW(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: ForceW
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION ForceW
  END INTERFACE

  INTERFACE
  FUNCTION ForceRho(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: ForceRho
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION ForceRho
  END INTERFACE

  INTERFACE
  FUNCTION ForceTh(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: ForceTh
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION ForceTh
  END INTERFACE

  INTERFACE 
  FUNCTION UStartE(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE UVW_Mod
    REAL(RealKind) :: UStartE
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION UStartE
  END INTERFACE

  INTERFACE 
  FUNCTION VStartE(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE UVW_Mod
    REAL(RealKind) :: VStartE
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION VStartE
  END INTERFACE

  INTERFACE 
  FUNCTION ThStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE ThProf_Mod
    REAL(RealKind) :: ThStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION ThStart
  END INTERFACE

  INTERFACE 
  FUNCTION EnStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE ThProf_Mod
    REAL(RealKind) :: EnStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION EnStart
  END INTERFACE

  INTERFACE 
  FUNCTION TStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE ThProf_Mod
    REAL(RealKind) :: TStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION TStart
  END INTERFACE

!  INTERFACE
!  FUNCTION RhoStart(x,y,z,zHeight,Time)
!    USE Kind_Mod
!    USE Rho_Mod
!    REAL(RealKind) :: RhoStart
!    REAL(RealKind) :: x,y,z,zHeight,Time
!  END FUNCTION RhoStart
!  END INTERFACE

  INTERFACE
  FUNCTION PreStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: PreStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION PreStart
  END INTERFACE

  INTERFACE
  FUNCTION PreProf(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: PreProf
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION PreProf
  END INTERFACE

  INTERFACE 
  FUNCTION TkeStart(x,y,z,zH,Time)
    USE Kind_Mod
    REAL(RealKind) :: TkeStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION TkeStart
  END INTERFACE

  INTERFACE 
  FUNCTION DisStart(x,y,z,zH,Time)
    USE Kind_Mod
    REAL(RealKind) :: DisStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION DisStart
  END INTERFACE

  INTERFACE 
  FUNCTION OmeStart(x,y,z,zH,Time)
    USE Kind_Mod
    REAL(RealKind) :: OmeStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION OmeStart
  END INTERFACE

  INTERFACE
  FUNCTION TkeHStart(x,y,z,zH,Time)
    USE Kind_Mod
    REAL(RealKind) :: TkeHStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION TkeHStart
  END INTERFACE

  INTERFACE
  FUNCTION TkeVStart(x,y,z,zH,Time)
    USE Kind_Mod
    REAL(RealKind) :: TkeVStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION TkeVStart
  END INTERFACE

  INTERFACE 
  FUNCTION LenStart(x,y,z,zH,Time)
    USE Kind_Mod
    REAL(RealKind) :: LenStart
    REAL(RealKind) :: x,y,z,zH,Time
  END FUNCTION LenStart
  END INTERFACE
  
  INTERFACE 
  FUNCTION QvStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE QvProf_Mod
    REAL(RealKind) :: QvStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION QvStart
  END INTERFACE

  INTERFACE 
  FUNCTION QcStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: QcStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION QcStart
  END INTERFACE

  INTERFACE 
  FUNCTION QrStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: QrStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION QrStart
  END INTERFACE

  INTERFACE
  FUNCTION QiStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: QiStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION QiStart
  END INTERFACE

  INTERFACE 
  FUNCTION QsStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: QsStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION QsStart
  END INTERFACE
  
  INTERFACE 
  FUNCTION NvStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: NvStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION NvStart
  END INTERFACE

  INTERFACE 
  FUNCTION NcStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: NcStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION NcStart
  END INTERFACE

  INTERFACE 
  FUNCTION NrStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: NrStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION NrStart
  END INTERFACE

  INTERFACE
  FUNCTION NiStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: NiStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION NiStart
  END INTERFACE

  INTERFACE 
  FUNCTION NsStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: NsStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION NsStart
  END INTERFACE

  INTERFACE 
  FUNCTION DStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    REAL(RealKind) :: DStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION DStart
  END INTERFACE

  INTERFACE 
  FUNCTION DHStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    REAL(RealKind) :: DHStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION DHStart
  END INTERFACE

  INTERFACE 
  FUNCTION DVStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    USE Rho_Mod
    REAL(RealKind) :: DVStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION DVStart
  END INTERFACE
  
  INTERFACE 
  FUNCTION HeightFun(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: HeightFun
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION HeightFun
  END INTERFACE

  INTERFACE 
  FUNCTION DummyStart(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: DummyStart
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION DummyStart
  END INTERFACE

  INTERFACE 
  FUNCTION Tracer1Start(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: Tracer1Start
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION Tracer1Start
  END INTERFACE

  INTERFACE 
  FUNCTION Tracer2Start(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: Tracer2Start
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION Tracer2Start
  END INTERFACE

  INTERFACE 
  FUNCTION DummyStart1(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: DummyStart1
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION DummyStart1
  END INTERFACE

  INTERFACE 
  FUNCTION DummyStart2(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: DummyStart2
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION DummyStart2
  END INTERFACE

  INTERFACE 
  FUNCTION DummyStart3(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: DummyStart3
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION DummyStart3
  END INTERFACE

  INTERFACE 
  FUNCTION DummyStart4(x,y,z,zHeight,Time)
    USE Kind_Mod
    REAL(RealKind) :: DummyStart4
    REAL(RealKind) :: x,y,z,zHeight,Time
  END FUNCTION DummyStart4
  END INTERFACE

  INTERFACE 
  SUBROUTINE PerturbProfile(VecC)
    USE DataType_Mod
    IMPLICIT NONE
    TYPE(Vector4Cell_T) :: VecC(:)
  END SUBROUTINE PerturbProfile
  END INTERFACE

  INTERFACE 
  SUBROUTINE SetBoundCells(BoundCellLoc)
    USE DataType_Mod
    IMPLICIT NONE
    TYPE(BoundCell_T) :: BoundCellLoc
  END SUBROUTINE SetBoundCells
  END INTERFACE

  INTERFACE
  SUBROUTINE SetParameterExample
    USE DataType_Mod
    IMPLICIT NONE
  END SUBROUTINE SetParameterExample
  END INTERFACE

  INTERFACE
  FUNCTION DampFun(z,Name)
    USE DataType_Mod
    IMPLICIT NONE
    REAL(RealKind) :: DampFun
    REAL(RealKind) :: z
    CHARACTER(*) :: Name
  END FUNCTION DampFun
  END INTERFACE

END MODULE Start_Mod

MODULE StartSoil_Mod
  INTERFACE 
  FUNCTION ThStartSoil(x,y,z,zHeight,zS,LandClass,SoilType,Time)
    USE Kind_Mod
    REAL(RealKind) :: ThStartSoil
    REAL(RealKind) :: x,y,z,zHeight,zS,Time
    INTEGER :: LandClass,SoilType
  END FUNCTION ThStartSoil
  END INTERFACE

  INTERFACE 
  FUNCTION QvStartSoil(x,y,z,zHeight,zS,LandClass,SoilType,Time)
    USE Kind_Mod
    REAL(RealKind) :: QvStartSoil
    REAL(RealKind) :: x,y,z,zHeight,zS,Time
    INTEGER :: LandClass,SoilType
  END FUNCTION QvStartSoil
  END INTERFACE

END MODULE StartSoil_Mod

MODULE Example_Mod
  USE Kind_Mod
  USE Parameter_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  USE RhoStart_Mod
  USE ThProf_Mod
  USE QvProf_Mod
  USE UVW_Mod
  USE Start_Mod
  USE StartSoil_Mod
  IMPLICIT NONE
CONTAINS
  
!FUNCTION RhoStart(x,y,z,zHeight,Time)
!  USE RhoProf_Mod
!  REAL(RealKind) :: RhoStart
!  REAL(RealKind) :: RhoLoc
!  REAL(RealKind) :: x,y,z,zHeight,Time
!  RhoLoc=RhoFun(x,y,z,zHeight,Time)
!  RhoStart=(RhoLoc-RhoProf(x,y,z,zHeight,Time))/(RhoLoc+Eps)
!  RhoStart=One
!END FUNCTION RhoStart

END MODULE Example_Mod

