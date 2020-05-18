MODULE Cloud_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: xC,yC,zC
  REAL(RealKind) :: xr,yr,zr
  REAL(RealKind) :: DeltaTh,uMax,vMax,ThInit
  REAL(RealKind) :: qc0,qv0,qt0
  REAL(RealKind) :: RelHum,N=0.0d0
  LOGICAL :: ProfIn=.FALSE.
  
  NAMELIST /Example/ xC &
                    ,yC &
                    ,zC &
                    ,xr &
                    ,yr &
                    ,zr &
                    ,qv0 &
                    ,qc0 &
                    ,qt0 &
                    ,RelHum &
                    ,N      &
                    ,uMax &
                    ,vMax &
                    ,ThInit &
                    ,DeltaTh &
                    ,ProFIn


END MODULE Cloud_Mod


SUBROUTINE PerturbProfile(VecC)

  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Parameter_Mod
  USE Floor_Mod
  USE Cloud_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)
END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE Cloud_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line

  REAL(RealKind) :: p1,p2,z1,z2,RhoProf

  qv0=1.d-2
  qc0=1.d-6
  qt0=1.d-2
  RelHum=0.0d0
! Find line
  xC=10000.d0
  yC=0.0d0
  zC=2000.0d0
  xr=2000.0d0
  yr=2000.0d0
  zr=2000.0d0
  uMax=0.0d0
  vMax=0.0d0
  ThInit=300.0d0
  DeltaTh=2.0d0
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

FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Cloud_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,S,ThLoc

  S=N*N/Grav
  ThLoc=thInit*exp(z*S)
  IF (N>Zero) THEN
    pLoc=p0*(One-Grav/(Cpd*ThInit*S)*(One-EXP(-S*z)))**(Cpd/Rd)
  ELSE
    pLoc=p0*(One-kappa*Grav*z/(Rd*thInit))**(Cpd/Rd)
  END IF
  RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE Cloud_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: S

  S=N*N/Grav
  ThProfFun=thInit*exp(z*S)
END FUNCTION ThProfFun

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Cloud_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: RhoType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S,pLoc,ThLoc,L


  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,RhoType,'RhoProf')
      Load=.FALSE.
    END IF
    RhoFun=ProfileEqual(cInt,z)
  ELSE
    S=N*N/Grav
    ThLoc=ThInit*exp(z*S)
    IF (N>Zero) THEN
      pLoc=p0*(One-Grav/(Cpd*ThInit*S)*(One-EXP(-S*z)))**(Cpd/Rd)
    ELSE
      pLoc=p0*(One-kappa*Grav*z/(Rd*ThInit))**(Cpd/Rd)
    END IF
    RhoFun=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END IF

END FUNCTION RhoFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE Cloud_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThProfFun
  REAL(RealKind), POINTER, SAVE :: cIntT(:,:)
  REAL(RealKind), POINTER, SAVE :: cIntqv(:,:)
  REAL(RealKind), POINTER, SAVE :: cIntRho(:,:)
  CHARACTER*10, SAVE :: tType,qvType,RhoType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: qvLoc,TLoc,pLoc,RhoLoc,RhoVLoc,kappaLoc,S

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cIntT,tType,'TProf')
      CALL ReadProfile(cIntqv,qvType,'QvProf')
      CALL ReadProfile(cIntRho,RhoType,'RhoProf')
      Load=.FALSE.
    END IF
    TLoc=ProfileEqual(cIntT,z)
    qvLoc=ProfileEqual(cIntqv,z)
    RhoLoc=ProfileEqual(cIntRho,z)
    RhoVLoc=RhoLoc*qvLoc
    pLoc=(Rd*RhoLoc+(Rv-Rd)*RhoVLoc)*TLoc
    KappaLoc=(Rd*RhoLoc+(Rv-Rd)*RhoVLoc) &
             /(Cpd*RhoLoc+(Cpv-Cpd)*RhoVLoc) 
    ThStart=Tloc*(pLoc/p0)**(-KappaLoc)*(RhoLoc+(Rv-Rd)/Rd*RhoVLoc)/RhoLoc
  ELSE
    S=N*N/Grav
    ThStart=thInit*exp(z*S)
  END IF
! L=SQRT(((x-xC)/xr)**2 &
!       +((z-zC)/zr)**2)
! IF (L<1.0d0) THEN
!   ThStart=ThStart*(1.0d0+DeltaTh*COS(Pi*L/2.0d0)**2/300.0d0)
! END IF

END FUNCTION ThStart

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE Cloud_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qvType,'QvProf')
      Load=.FALSE.
    END IF
    qvStart=ProfileEqual(cInt,z)
  ELSE
    IF (z>=2000.0d0.AND.((x<=20000.AND.x>=-20000).OR.(x>=60000.AND.x<=80000))) THEN
      qvStart=1.0d-4
    ELSE
      qvStart=1.0d-8
    END IF
  END IF

END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE Cloud_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

! IF (ProfIn) THEN
!   IF (Load) THEN
!     CALL ReadProfile(cInt,qcType,'QcProf')
!     Load=.FALSE.
!   END IF
!   qcStart=ProfileEqual(cInt,z)
! ELSE
!  qcStart=0.0d0
! END IF

   qcStart=1.0d-6
    IF (z>=2000.0d0.AND.((x<=20000.AND.x>=-20000).OR.(x>=60000.AND.x<=80000))) THEN
      qcStart=1.0d-4
    ELSE
      qcStart=1.0d-8
    END IF
END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=Zero

END FUNCTION QiStart

FUNCTION UStart(x,y,z,zHeight,Time)

  USE Cloud_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  UStart=uMax

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE Cloud_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  VStart=VMax

END FUNCTION VStart



FUNCTION DStart(x,y,z,zHeight,Time)

  USE Cloud_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=Zero

END FUNCTION DStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Cloud_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Cloud_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: Rho1,Rho2
  Rho1=RhoFun(x,y,z,zHeight,Time)
  Rho2=RhoProf(x,y,z,zHeight,Time)
! RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
  RhoStart=Rho1-Rho2
END FUNCTION RhoStart


FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(Cpd/Rd)
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE Cloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(Cpd/Rd)
END FUNCTION PreProf




