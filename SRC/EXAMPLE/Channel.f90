MODULE Channel_Mod

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

END MODULE Channel_Mod

SUBROUTINE PerturbProfile(VecC)

  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Parameter_Mod
  USE Floor_Mod
  USE Channel_Mod

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
  REAL(RealKind) :: xPLoc,yPLoc,zPloc,L,Delta
  REAL(RealKind) :: qv

  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Rho=>RhoCell(ib)%c
    RhoV=>VecC(ib)%Vec(RhoVPos)%c
    RhoL=>VecC(ib)%Vec(RhoCPos)%c
    ThVirt=>VecC(ib)%Vec(thPos)%c
    DO iz=iz0+1,iz1
      zpLoc=zP(iz-1)+0.5d0*dz(iz)
      DO iy=iy0+1,iy1
        ypLoc=yP(iy-1)+0.5d0*dy(iy)
        DO ix=ix0+1,ix1
          xpLoc=xP(ix-1)+0.5d0*dx(ix)
          L=SQRT(((xPLoc-xC)/xr)**2 &
                +((zPLoc-zC)/zr)**2)
          RhoLLoc=RhoL(ix,iy,iz,1)
          IF (L<1.0d0.AND.DeltaTh>0.0d0) THEN
            Delta=DeltaTh*COS(Pi*L/2.0d0)**2/300.0d0
            RhoLoc=Rho(ix,iy,iz,1)
            RhoVLoc=RhoV(ix,iy,iz,1)
            RhoDry=RhoLoc-RhoVLoc-RhoLLoc
            KappaLoc=(Rd*(RhoLoc-RhoLLoc) &
                    +(Rv-Rd)*RhoVLoc) &
                    /(Cpd*RhoLoc+(Cpv-Cpd)*RhoVLoc &
                     +(Cpl-Cpd)*RhoLLoc)
            pLoc=p0*(Rd*ThVirt(ix,iy,iz,1)/p0)**(One/(One-KappaLoc))
            ThDensLoc=ThVirt(ix,iy,iz,1)/(RhoLoc+Eps)
            ThDensNew=ThDensLoc*(One+Delta)
            RhoLoc=(pLoc/p0)**(One-KappaLoc)*p0/(Rd*ThDensNew)
            VecC(ib)%Vec(RhoPos)%c(ix,iy,iz,1)=RhoLoc-RhoProfG(ib)%c(ix,iy,iz,1)
            write(*,*) ThVirt(ix,iy,iz,1),ThDensNew*RhoLoc
            ThVirt(ix,iy,iz,1)=ThDensNew*RhoLoc
          END IF
        END DO
      END DO
    END DO
  END DO
END SUBROUTINE PerturbProfile

SUBROUTINE InputExample(FileName)
  USE Channel_Mod
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
  xC=0.0d0
  yC=15000.0d0
  zC=1400.0d0
  xr=5000.0d0
  yr=5000.0d0
  zr=700.0d0
  uMax=0.0d0
  vMax=0.0d0
  ThInit=303.0d0
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
  USE ReadProfile_Mod
  USE Channel_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: RhoProfType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,RhoProfType,'RhoProf')
      Load=.FALSE.
    END IF
    RhoProf=ProfileEqual(cInt,z)
  ELSE
    RhoProf=1.14d0
  END IF

END FUNCTION RhoProf

FUNCTION ThProfFun(x,y,z,zHeight,Time)
  USE Parameter_Mod
  USE ReadProfile_Mod
  USE Channel_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,thType,'ThProf')
      Load=.FALSE.
    END IF
    ThProfFun=ProfileEqual(cInt,z)
  ELSE
    S=N*N/Grav
    ThProfFun=thInit
  END IF
END FUNCTION ThProfFun

FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE ReadProfile_Mod
  USE RhoProf_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qvType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qvType,'QvProf')
      Load=.FALSE.
    END IF
    qvProfFun=ProfileEqual(cInt,z)
  ELSE
   qvProfFun=0.0d0
  END IF
END FUNCTION QvProfFun

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE Channel_Mod
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
    RhoFun=1.14d0
  END IF

END FUNCTION RhoFun

FUNCTION ThStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE ThProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  ThStart=ThProfFun(x,y,z,zHeight,Time)
END FUNCTION ThStart

FUNCTION QvStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE QvProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QvStart=QvProfFun(x,y,z,zHeight,Time)
END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

  qcStart=1.0d-23
END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)
  USE Parameter_Mod

  IMPLICIT NONE
  REAL(RealKind) :: QiStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QiStart=Zero
END FUNCTION QiStart

FUNCTION UStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStart=umax
END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStart=0.0
END FUNCTION VStart

FUNCTION DStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DStart=1.0d+0
END FUNCTION DStart

FUNCTION UStartE(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=0.1
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=Zero
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=Zero
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=Zero
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=0.0d0
END FUNCTION DummyStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: Rho1,Rho2
  Rho1=RhoFun(x,y,z,zHeight,Time)
  Rho2=RhoProf(x,y,z,zHeight,Time)
  RhoStart=Rho1-Rho2
END FUNCTION RhoStart

FUNCTION PreStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE PreProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreStart=PreProfFun(x,y,z,zHeight,Time)
END FUNCTION PreStart

FUNCTION PreProfFun(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: PrePType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,PrePType,'PreProf')
      Load=.FALSE.
    END IF
    PreProfFun=ProfileEqual(cInt,z)
  ELSE
    PreProfFun=p0*(One-kappa*Grav*z/(Rd*thInit))**(Cpd/Rd)
  END IF
END FUNCTION PreProfFun

FUNCTION TStart(x,y,z,zHeight,Time)
  USE Channel_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: TStartType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,TStartType,'TProf')
      Load=.FALSE.
    END IF
    TStart=ProfileEqual(cInt,z)
  ELSE
    TStart=thInit
  END IF

END FUNCTION TStart
