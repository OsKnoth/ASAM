MODULE StabCloud_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: xC,yC,zC
  REAL(RealKind) :: xLC,xRC
  REAL(RealKind) :: xr,yr,zr
  REAL(RealKind) :: DeltaTh,uMax,vMax,ThInit
  REAL(RealKind) :: qc0,qv0,qt0
  REAL(RealKind) :: RelHum,N=0.0d0
  LOGICAL :: ProfIn=.FALSE.
  LOGICAL :: ProfInP=.FALSE.
  CHARACTER*50 :: FileProf
  CHARACTER*50 :: FileProfP
  REAL(RealKind) :: PiLoc
  REAL(RealKind) :: TimePeriod=600.0d0
  
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
                    ,ProfIn &
                    ,ProfInP &
                    ,FileProf &
                    ,FileProfP &
                    ,xLC &
                    ,xRC  &
                    ,TimePeriod 

CONTAINS
FUNCTION TimeMix(Time)

  REAL(RealKind) :: TimeMix
  REAL(RealKind) :: Time
  
  TimeMix=MAX(SIN(Two*PiLoc*Time/TimePeriod),Zero)
END FUNCTION TimeMix
END MODULE StabCloud_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE StabCloud_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc

END SUBROUTINE SetBoundCells

SUBROUTINE PerturbProfile(VecC)

  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Parameter_Mod
  USE Floor_Mod
  USE StabCloud_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)
  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: ThVirt(:,:,:,:)
  REAL(RealKind) :: RhoLoc,RhoNew,RhoVLoc,RhoLLoc,RhoDry,RhoDryNew
  REAL(RealKind) :: KappaLoc,pLoc,pvs,rt,rvs
  REAL(RealKind) :: ThDensLoc,ThDensNew,ThLoc,ThNew,TLoc
  REAL(RealKind) :: xPLoc,zPloc,L,Delta
  REAL(RealKind) :: qv

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    ThVirt=>VecC(ib)%Vec(thPos)%c
    call srand(ib) 
    DO ix=ix0+1,ix1
     DO iy=iy0+1,iy1    
      DO iz=iz0+1,10    
       ThVirt(ix,iy,iz,1)=ThVirt(ix,iy,iz,1)+DeltaTh*(0.5-rand())
      END DO
     END DO
    END DO
  END DO 
END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE StabCloud_Mod
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
  uMax=10.0d0
  vMax=0.0d0
  ThInit=300.0d0
  DeltaTh=0.0d0
  PiLoc=Pi
  WRITE(*,*) 'PiLoc',PiLoc
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
  USE ReadProfile_Mod
  USE StabCloud_Mod
  USE Rho_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: RhoType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S,pLoc,ThLoc,L


  IF (ProfIn) THEN
    IF (Load) THEN
      Profile=FileProf
      CALL ReadProfile(cInt,RhoType,'RhoProf')
      Load=.FALSE.
    END IF
    RhoProf=ProfileEqual(cInt,z)
  ELSE
    S=N*N/Grav
    ThLoc=ThInit*exp(z*S)
    IF (N>Zero) THEN
      pLoc=p0*(One-Grav/(Cpd*ThInit*S)*(One-EXP(-S*z)))**(Cpd/Rd)
    ELSE
      pLoc=p0*(One-kappa*Grav*z/(Rd*ThInit))**(Cpd/Rd)
    END IF
    RhoProf=pLoc/((pLoc/p0)**kappa*Rd*ThLoc)
  END IF

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE ReadProfile_Mod
  USE StabCloud_Mod
  USE Start_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      Profile=FileProf
      CALL ReadProfile(cInt,thType,'ThProf')
      Load=.FALSE.
    END IF
    ThProfFun=ProfileEqual(cInt,z)
  ELSE
    S=N*N/Grav
    ThProfFun=thInit*exp(z*S)
  END IF

END FUNCTION ThProfFun

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE StabCloud_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  REAL(RealKind), POINTER, SAVE :: cIntP(:,:)
  CHARACTER*10, SAVE :: RhoType
  LOGICAL, SAVE :: Load=.TRUE.
  LOGICAL, SAVE :: LoadP=.TRUE.
  REAL(RealKind) :: S,pLoc,ThLoc,L


  IF (ProfIn) THEN
    IF (Load) THEN
      Profile=FileProf
      CALL ReadProfile(cInt,RhoType,'RhoProf')
      Load=.FALSE.
    END IF
    IF (LoadP) THEN
      Profile=FileProfP
      CALL ReadProfile(cIntP,RhoType,'RhoProf')
      LoadP=.FALSE.
    END IF
    IF (x<xLC.OR.x>xRC) THEN
      RhoFun=ProfileEqual(cInt,z)
    ELSE
      RhoFun=ProfileEqual(cIntP,z)
    END IF
    IF (Time>Zero) THEN
      RhoFun=(One-TimeMix(Time))*ProfileEqual(cInt,z) &
                  +TimeMix(Time)*ProfileEqual(cIntP,z)
    END IF
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

  USE StabCloud_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: ThProfFun
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  REAL(RealKind), POINTER, SAVE :: cIntP(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  LOGICAL, SAVE :: LoadP=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      Profile=FileProf 
      CALL ReadProfile(cInt,thType,'ThProf')
      Load=.FALSE.
    END IF
    IF (LoadP) THEN
      Profile=FileProfP 
      CALL ReadProfile(cIntP,thType,'ThProf')
      LoadP=.FALSE.
    END IF
    IF (x<xLC.OR.x>xRC) THEN
      ThStart=ProfileEqual(cInt,z)
    ELSE
      ThStart=ProfileEqual(cIntP,z)
    END IF
    IF (Time>Zero) THEN
      ThStart=(One-TimeMix(Time))*ProfileEqual(cInt,z) &
                  +TimeMix(Time)*ProfileEqual(cIntP,z)
    END IF
  ELSE
    S=N*N/Grav
    ThStart=thInit*exp(z*S)
  END IF
! L=SQRT(((x-xC)/xr)**2 &
!       +((z-zC)/zr)**2)
! IF (L<1.0d0) THEN
!   ThStart=ThStart*(1.0d0+DeltaTh*COS(PiLoc*L/2.0d0)**2/300.0d0)
! END IF

END FUNCTION ThStart

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE StabCloud_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  REAL(RealKind), POINTER, SAVE :: cIntP(:,:)
  CHARACTER*10, SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.
  LOGICAL, SAVE :: LoadP=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      Profile=FileProf
      CALL ReadProfile(cInt,qvType,'QvProf')
      Load=.FALSE.
    END IF
    IF (LoadP) THEN
      Profile=FileProfP
      CALL ReadProfile(cIntP,qvType,'QvProf')
      LoadP=.FALSE.
    END IF
    IF (x<xLC.OR.x>xRC) THEN
      qvStart=ProfileEqual(cInt,z)
    ELSE
      qvStart=ProfileEqual(cIntP,z)
    END IF
    IF (Time>Zero) THEN
      qvStart=(One-TimeMix(Time))*ProfileEqual(cInt,z) &
                  +TimeMix(Time)*ProfileEqual(cIntP,z)
    END IF
  ELSE
   qvStart=0.0d0
  END IF

END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE StabCloud_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  REAL(RealKind), POINTER, SAVE :: cIntP(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.
  LOGICAL, SAVE :: LoadP=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      Profile=FileProf
      CALL ReadProfile(cInt,qcType,'QcProf')
      Load=.FALSE.
    END IF
    IF (LoadP) THEN
      Profile=FileProfP
      CALL ReadProfile(cIntP,qcType,'QcProf')
      LoadP=.FALSE.
    END IF
    IF (x<xLC.OR.x>xRC) THEN
      qcStart=ProfileEqual(cInt,z)
    ELSE
      qcStart=ProfileEqual(cIntP,z)
    END IF
    IF (Time>Zero) THEN
      qcStart=(One-TimeMix(Time))*ProfileEqual(cInt,z) &
                  +TimeMix(Time)*ProfileEqual(cIntP,z)
    END IF
  ELSE
   qcStart=0.0d0
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

  USE StabCloud_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  UStart=uMax

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE StabCloud_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  VStart=VMax

END FUNCTION VStart



FUNCTION DStart(x,y,z,zHeight,Time)

  USE StabCloud_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=Zero

END FUNCTION DStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=p0*(One-kappa*Grav*z/(Rd*thInit))**(Cpd/Rd)
END FUNCTION PreStart

FUNCTION PreProf(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=p0*(One-kappa*Grav*z/(Rd*thInit))**(Cpd/Rd)
END FUNCTION PreProf

FUNCTION TStart(x,y,z,zHeight,Time)
  USE StabCloud_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart=0.0d0
END FUNCTION TStart




