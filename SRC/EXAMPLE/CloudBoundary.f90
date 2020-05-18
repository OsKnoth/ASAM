MODULE CloudBoundary_Mod

  USE Kind_Mod
  USE Parameter_Mod
  USE Domain_Mod
  USE Physics_Mod
  USE Thermodynamic_Mod
  USE SoilData_Mod

  IMPLICIT NONE 

  REAL(RealKind) :: uMax,vMax
  REAL(RealKind) :: uMin,vMin
  REAL(RealKind) :: zFree=250.0d0
  REAL(RealKind) :: z2=1500.0d0
  REAL(RealKind) :: z3=4000.0d0
  REAL(RealKind) :: TkeHMax=1.0d-2
  REAL(RealKind) :: TkeVMax=1.0d-2
  REAL(RealKind) :: LenMax=1.0d0
  REAL(RealKind) :: AeroStart=0.0d0
  REAL(RealKind) :: HeightTracer1Top=0.0d0
  REAL(RealKind) :: HeightTracer1Bottom=0.0d0
  REAL(RealKind) :: HeightTracer2Top=0.0d0
  REAL(RealKind) :: HeightTracer2Bottom=0.0d0
  REAL(RealKind) :: DeltaT=0.02d0
  REAL(RealKind) :: HeightDeltaT=1000.0d0
  REAL(RealKind) :: DeltaQv=0.05d-3
  REAL(RealKind) :: HeightDeltaQv=1000.0d0
  REAL(RealKind) :: offset_y=1.0d0
  REAL(RealKind) :: offset_y1=1.0d0
  REAL(RealKind) :: intenz=0.1d0
  REAL(RealKind) :: phase(1000)
  REAL(RealKind) :: phase1(1000)
  REAL(RealKind) :: TurbHeight=1000.0d0
  LOGICAL :: ProfIn=.FALSE.
  LOGICAL :: Perturb=.FALSE.
  LOGICAL :: TurbInFlow=.FALSE.
  LOGICAL :: LogWindProf=.FALSE.
  LOGICAL :: WindProf=.FALSE.
  LOGICAL :: BomexWind=.FALSE.
  
  NAMELIST /Example/ uMax &
                    ,vMax &
                    ,uMin &
                    ,vMin &
                    ,zFree &
                    ,z2 &
                    ,z3 &
                    ,AeroStart &
                    ,HeightTracer1Top &
                    ,HeightTracer1Bottom &
                    ,HeightTracer2Top &
                    ,HeightTracer2Bottom &
                    ,DeltaT &
                    ,HeightDeltaQv &
                    ,DeltaQv &
                    ,HeightDeltaT &
                    ,ProfIn &
                    ,TurbInFlow &
                    ,LogWindProf &
                    ,WindProf &
                    ,BomexWind &
                    ,TurbHeight &
                    ,Perturb 

END MODULE CloudBoundary_Mod

SUBROUTINE SetBoundCells(BoundCellLoc)

  USE CloudBoundary_Mod
  USE DataType_Mod
  IMPLICIT NONE
  TYPE(BoundCell_T) :: BoundCellLoc
  SELECT CASE (Problem)
    CASE ('RICO')
      BoundCellLoc%DragH=0.001094d0
      BoundCellLoc%DragQ=0.001133d0
      BoundCellLoc%DragM=0.001229d0
      BoundCellLoc%TeS=299.8d0
      BoundCellLoc%ThetaS=299.8d0
      BoundCellLoc%qv=SaturVapor(BoundCellLoc%TeS)/(Rv*BoundCellLoc%TeS)
    CASE DEFAULT  
  END SELECT    

END SUBROUTINE SetBoundCells


SUBROUTINE PerturbProfile(VecC)

  USE Physics_Mod
  USE Thermodynamic_Mod
  USE DataType_Mod
  USE Parameter_Mod
  USE Floor_Mod
  USE CloudBoundary_Mod
  IMPLICIT NONE

  TYPE(Vector4Cell_T) :: VecC(:)

  INTEGER :: i,ib,ibLoc,ix,iy,iz
  REAL(RealKind), POINTER :: Rho(:,:,:,:)
  REAL(RealKind), POINTER :: RhoV(:,:,:,:)
  REAL(RealKind), POINTER :: RhoL(:,:,:,:)
  REAL(RealKind), POINTER :: RhoR(:,:,:,:)
  REAL(RealKind), POINTER :: ThVirt(:,:,:,:)
  REAL(RealKind), POINTER :: En(:,:,:,:)
  REAL(RealKind) :: RhoLoc,RhoVLoc,RhoLLoc,RhoDry
  REAL(RealKind) :: KappaLoc,pLoc,rvLoc,rlLoc,rtLoc
  REAL(RealKind) :: ThDensLoc,ThLoc,TLoc,ThEquivLoc
  REAL(RealKind) :: zPloc,r

  IF (Perturb) THEN
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    Rho=>RhoCell(ib)%c
    RhoV=>VecC(ib)%Vec(RhoVPos)%c
    RhoL=>VecC(ib)%Vec(RhoCPos)%c
    RhoR=>VecC(ib)%Vec(RhoRPos)%c
    ThVirt=>VecC(ib)%Vec(thPos)%c
    IF (EnPos>0) THEN
      En=>VecC(ib)%Vec(EnPos)%c
    END IF 
    DO iz=iz0+1,iz1
      zpLoc=zP(iz-1)+0.5d0*dz(iz)
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          RhoLLoc=RhoL(ix,iy,iz,1)+RhoR(ix,iy,iz,1)
          RhoLoc=Rho(ix,iy,iz,1)
          RhoVLoc=RhoV(ix,iy,iz,1)
          IF (zpLoc<HeightDeltaQv) THEN
            CALL Random_number(r)
            RhoVLoc=RhoVLoc+(r-0.5d0)*DeltaQv
          END IF
          RhoDry=RhoLoc-RhoVLoc-RhoLLoc
          rvLoc=RhoVLoc/RhoDry+Eps
          rlLoc=RhoLLoc/RhoDry+Eps
          rtLoc=rvLoc+rlLoc
          KappaLoc=(Rd*(RhoLoc-RhoLLoc) &
                   +(Rv-Rd)*RhoVLoc) &
                   /(Cpd*RhoLoc+(Cpv-Cpd)*RhoVLoc &
                   +(Cpl-Cpd)*RhoLLoc)
          pLoc=p0*(Rd*ThVirt(ix,iy,iz,1)/p0)**(One/(One-KappaLoc))
          TLoc=pLoc/(Rd*RhoDry+Rv*RhoVLoc)
          IF (zpLoc<HeightDeltaT) THEN
            CALL Random_number(r)
            TLoc=TLoc+(r-0.5d0)*DeltaT
          END IF
          SELECT CASE(ThetaKind)
          CASE('Density')
            ThVirt(ix,iy,iz,1)=RhoLoc*TLoc*(p0/pLoc)**KappaLoc*(One+rvLoc*Rv/Rd)/(One+rtLoc)
          CASE('PotTemp')
            ThVirt(ix,iy,iz,1)=TLoc*(p0/pLoc)**(Rd/Cpd)*RhoLoc
          CASE('Equiv')
            ThEquivLoc=ThEquiv(TLoc,RhoDry,RhoVLoc,RhoLLoc)
            ThVirt(ix,iy,iz,1)=ThEquivLoc*RhoLoc
          CASE('Energy')
            ThVirt(ix,iy,iz,1)= &
              (RhoDry*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl)*TLoc+RhoVLoc*L00 &
              +Half*RhoLoc*(zP(iz-1)+zP(iz))*Grav
          CASE('PreEn')
            ThVirt(ix,iy,iz,1)=pLoc
            En(ix,iy,iz,1)= &
              (RhoDry*Cvd+RhoVLoc*Cvv+RhoLLoc*Cpl)*TLoc+RhoVLoc*L00 &
              +Half*RhoLoc*(zP(iz-1)+zP(iz))*Grav+Half*RhoLoc*uMax*uMax
          END SELECT
        END DO 
      END DO 
    END DO 
  END DO 
  END IF

END SUBROUTINE PerturbProfile


SUBROUTINE InputExample(FileName)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  CHARACTER(*) :: FileName
  INTEGER :: Pos
  CHARACTER(300) :: Line
  INTEGER :: i
  REAL(RealKind) :: r

  offset_y=domain%y0
  offset_y1=domain%y1

   DO i=1,1000
    Call Random_NUmber(r)
    phase(i)=r
!WRITE(*,*)'phase',i,phase(i)
    phase1(i)=2.0*r-5.0
!WRITE(*,*)i,(phase1(i))
    phase1(i)=3.0*10.0**(phase1(i))
  ENDDO 

! Find line
  uMax=0.0d0
  vMax=0.0d0
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

  USE CloudBoundary_Mod
  IMPLICIT NONE

  REAL(RealKind) :: UStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: z0

  REAL(REALKIND) :: SpeedStartU1
  REAL(RealKind) :: SpeedStartU,Speed
  REAL(RealKind) :: wy,kmean,AMP=0.8 !0.75
  INTEGER        :: s

  IF (TurbInFlow) THEN
    SpeedStartU=0.0d0
    Speed = (SQRT(uMax**Two+vMax**Two))
    DO s = 1,100 !100 !10 !k_end
       wy = phase1(s)   !kmean*(Speed)
       kmean = wy/speed !(k/(offset_y1-offset_y)) !Einheit jetzt [1/m]
!       AMP = 0.55*0.35*((79.0*((wy*4.0/5.0)/(1.0+4.7*(0.6))))/ &
!            (1.0+263.0*(((wy*4.0/5.0)/(1.0+4.7*(0.6)))**(5.0/3.0)))) &
!           * (((1.0+2.5*((0.6)**(0.6)))/(1.0+4.7*(0.6)))**(2.0/3.0))
    SpeedStartU1 =  AMP*SIN(kmean*((y-offset_y))*2.0*3.1415 + phase(s)*offset_y1)* &
                 (COS(Time*(wy)*2.0*3.1415))
    SpeedStartU=SpeedStartU + SpeedStartU1
    ENDDO
    SpeedStartU=intenz*SpeedStartU/1.0
    SpeedStartU=SpeedStartU + uMax ! Speed

    IF(Time.lt.1.0d-4) SpeedStartU= uMax ! Speed
    IF((y-offset_y).lt.1.0d4.or.(y-offset_y).gt.1.0d5.or.z>TurbHeight) SpeedStartU=uMax
  ELSE
    SpeedStartU=uMax
  END IF

  IF (LogWindProf.AND..NOT.WindProf) THEN
    IF (zHeight>Zero) THEN
      z0=0.5d0 ! land
    ELSE
      z0=0.0002d0 ! sea
    END IF

    IF (uMax<Zero) THEN
      UStart=MAX(SpeedStartU,SpeedStartU*LOG((z-zHeight)/z0)/LOG(zFree/z0))
    ELSE
      UStart=MIN(SpeedStartU,SpeedStartU*LOG((z-zHeight)/z0)/LOG(zFree/z0))
    END IF
  ELSE
    UStart=SpeedStartU
  END IF

  IF (WindProf) THEN
    IF (zHeight>Zero) THEN
      z0=0.5d0 ! land
    ELSE
      z0=0.0002d0 ! sea
    END IF

    IF ((z-zHeight)<=zFree) THEN
    IF (uMax<Zero) THEN
      UStart=MAX(SpeedStartU,SpeedStartU*LOG((z-zHeight)/z0)/LOG(zFree/z0))
    ELSE
      UStart=MIN(SpeedStartU,SpeedStartU*LOG((z-zHeight)/z0)/LOG(zFree/z0))
    END IF
    ELSE IF ((z-zHeight)<=z2) THEN
      UStart=SpeedStartU
    ELSE IF ((z-zHeight)<=z3) THEN
      UStart=SpeedStartU+(uMin-uMax)/(z3-z2)*((z-zHeight)-z2)
    ELSE
      UStart=uMin
    END IF 
  END IF

  IF (BomexWind) THEN
    IF (z-zHeight<700.0d0) THEN
      UStart=SpeedStartU
    ELSE
      UStart=SpeedStartU+1.8d-3*(z-zHeight-700.0d0) 
    END IF 
  END IF
>>>>>>> d1af75040bfa212804a901250156e2f33b403426

END FUNCTION UStart

FUNCTION VStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  IMPLICIT NONE

  REAL(RealKind) :: VStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: z0

  REAL(REALKIND) :: SpeedStartV1
  REAL(RealKind) :: SpeedStartV,Speed
  REAL(RealKind) :: wy,kmean,AMP=0.8 !0.75
  INTEGER        :: s

  IF (TurbInFlow) THEN
  SpeedStartV=0.0d0
  Speed = (SQRT(uMax**Two+vMax**Two))
  DO s = 1,100 !100 !k_end
     wy = phase1(s) !kmean*(Speed)
     kmean = wy/speed !(k/(offset_y1-offset_y)) !Einheit jetzt [1/m]
!       AMP = 0.55*0.35*((79.0*((wy*4.0/5.0)/(1.0+4.7*(0.6))))/ &
!            (1.0+263.0*(((wy*4.0/5.0)/(1.0+4.7*(0.6)))**(5.0/3.0)))) &
!           * (((1.0+2.5*((0.6)**(0.6)))/(1.0+4.7*(0.6)))**(2.0/3.0))
SpeedStartV1 =  AMP*COS(kmean*((y-offset_y))*2.0*3.1415 + phase(s)*offset_y1)* &
               (SIN(Time*(wy)*2.0*3.1415))
  SpeedStartV=SpeedStartV + SpeedStartV1
  ENDDO
  SpeedStartV=intenz*SpeedStartV/1.0
  SpeedStartV=SpeedStartV + vMax

  IF(Time.lt.1.0d-4) SpeedStartV=0.0
  IF((y-offset_y).lt.1.0d4.or.(y-offset_y).gt.1.0d5.or.z>TurbHeight) SpeedStartV=vMax
  ELSE
    SpeedStartV=vMax
  END IF

  IF (LogWindProf) THEN
    IF (zHeight>Zero) THEN
      z0=0.5d0 ! land
    ELSE
      z0=0.0002 ! sea
    END IF

    IF (vMax<Zero) THEN
      VStart=MAX(SpeedStartV,SpeedStartV*LOG((z-zHeight)/z0)/LOG(zFree/z0))
    ELSE
      VStart=MIN(SpeedStartV,SpeedStartV*LOG((z-zHeight)/z0)/LOG(zFree/z0))
    END IF
  ELSE
    VStart=SpeedStartV
  END IF

  IF (WindProf) THEN
    IF (zHeight>Zero) THEN
      z0=0.5d0 ! land
    ELSE
      z0=0.0002d0 ! sea
    END IF
    
    IF ((z-zHeight)<zFree) THEN
    IF (vMax<Zero) THEN
      VStart=MAX(SpeedStartV,SpeedStartV*LOG((z-zHeight)/z0)/LOG(zFree/z0))
    ELSE
      VStart=MIN(SpeedStartV,SpeedStartV*LOG((z-zHeight)/z0)/LOG(zFree/z0))
    END IF
    ELSE IF ((z-zHeight)<z2) THEN
      VStart=SpeedStartV
    ELSE IF ((z-zHeight)<z3) THEN
      VStart=SpeedStartV+(vMin-vMax)/(z3-z2)*((z-zHeight)-z2)
    ELSE
      VStart=vMin
    END IF
  END IF
>>>>>>> d1af75040bfa212804a901250156e2f33b403426

END FUNCTION VStart

FUNCTION ThProfFun(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  USE ReadProfile_Mod

  IMPLICIT NONE

  REAL(RealKind) :: ThProfFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.

  REAL(RealKind) :: pLoc,TLoc

  ThProfFun=Zero

END FUNCTION ThProfFun

FUNCTION ThStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  USE Rho_Mod
  USE QvProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: ThStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: Rad,pLoc,TLoc
  REAL(RealKind) :: RhoLoc,RhoDLoc,RhoVLoc,RhoLLoc
  REAL(RealKind) :: rvLoc,rlLoc 
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: thType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      SELECT CASE(ThetaKind)
      CASE('Density')
        CALL ReadProfile(cInt,thType,'ThDensProf')
      CASE('Equiv')
        CALL ReadProfile(cInt,thType,'ThetaEProf')
      CASE('Energy')
        CALL ReadProfile(cInt,thType,'EnergyProf')
      CASE('PreEn')
        CALL ReadProfile(cInt,thType,'PreProf')
      CASE DEFAULT
        WRITE(*,*) 'ThetaKind pruefen.' 
      END SELECT
      Load=.FALSE.
    END IF
    ThStart=ProfileEqual(cInt,z)
  ELSE
    ThStart=0.0d0
  END IF

END FUNCTION ThStart

FUNCTION RhoFun(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE CloudBoundary_Mod
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
    RhoFun=One
  END IF
END FUNCTION RhoFun


FUNCTION PreStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: PreStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: tType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,tType,'PreProf')
      Load=.FALSE.
    END IF
    PreStart=ProfileEqual(cInt,z)
  ELSE
    PreStart=p0*(One-kappa*Grav*z/(Rd*300.0d0))**(Cpd/Rd)
  END IF

END FUNCTION PreStart


FUNCTION PreProf(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: PreProf
  REAL(RealKind) :: x,y,z,zHeight,Time
  PreProf=Zero
END FUNCTION PreProf


FUNCTION TFun(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  USE Rho_Mod
  USE PreProf_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: TFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  REAL(RealKind) :: pLoc,RhoLoc,GammaDry
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: tType
  LOGICAL, SAVE :: Load=.TRUE.
  REAL(RealKind) :: S

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,tType,'TProf')
      Load=.FALSE.
    END IF
    TFun=ProfileEqual(cInt,z)
  ELSE
    TFun=300.0d0
  END IF

END FUNCTION TFun


FUNCTION QvProfFun(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  USE Rho_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QvProfFun

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: TLoc,RhoLoc

  QvProfFun=0.0d0

END FUNCTION QvProfFun

FUNCTION QvStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  USE QvProf_Mod
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
    QvStart=0.0d0
  END IF

END FUNCTION QvStart

FUNCTION QcStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QcStart

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qcType,'QcProf')
      Load=.FALSE.
    END IF
    qcStart=ProfileEqual(cInt,z)
  ELSE
    QcStart=0.0d0
  END IF

END FUNCTION QcStart

FUNCTION QiStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QiStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QiStart=Zero

END FUNCTION QiStart

FUNCTION QsStart(x,y,z,zHeight,Time)

  USE Parameter_Mod
  IMPLICIT NONE

  REAL(RealKind) :: QsStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  QsStart=Zero

END FUNCTION QsStart

FUNCTION NvStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NvStart
  REAL(RealKind) :: QcLoc

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

  QcLoc=Zero
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qcType,'QcProf')
      Load=.FALSE.
    END IF
    QcLoc=ProfileEqual(cInt,z)
  END IF
  NvStart=AeroStart

END FUNCTION NvStart

FUNCTION NcStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NcStart
  REAL(RealKind) :: QcLoc

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind), POINTER, SAVE :: cInt(:,:)
  CHARACTER*10, SAVE :: qcType
  LOGICAL, SAVE :: Load=.TRUE.

  QcLoc=Zero
  IF (ProfIn) THEN
    IF (Load) THEN
      CALL ReadProfile(cInt,qcType,'QcProf')
      Load=.FALSE.
    END IF
    QcLoc=ProfileEqual(cInt,z)
  END IF
  NcStart=Zero

END FUNCTION NcStart

FUNCTION NrStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NrStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NrStart=Zero

END FUNCTION NrStart

FUNCTION NiStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NiStart

  REAL(RealKind) :: x,y,z,zHeight,Time


END FUNCTION NiStart

FUNCTION NsStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  IMPLICIT NONE

  REAL(RealKind) :: NsStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  NsStart=Zero

END FUNCTION NsStart

FUNCTION OmeStart(lam,phi,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE 
  REAL(RealKind) :: OmeStart
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  OmeStart = 0.0
END FUNCTION OmeStart

FUNCTION DStart(x,y,z,zHeight,Time)

  USE CloudBoundary_Mod
  IMPLICIT NONE

  REAL(RealKind) :: DStart

  REAL(RealKind) :: x,y,z,zHeight,Time

  DStart=Zero

END FUNCTION DStart


FUNCTION RhoProf(x,y,z,zHeight,Time)

  USE Parameter_Mod
  USE CloudBoundary_Mod
  IMPLICIT NONE

  REAL(RealKind) :: RhoProf

  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: pLoc,TLoc,ThLoc

  RhoProf=Zero

END FUNCTION RhoProf


FUNCTION UStartE(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: UStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  UStartE=UStart(x,y,z,zHeight,Time)
END FUNCTION UStartE

FUNCTION VStartE(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  USE UVW_Mod
  IMPLICIT NONE
  REAL(RealKind) :: VStartE
  REAL(RealKind) :: x,y,z,zHeight,Time
  VStartE=VStart(x,y,z,zHeight,Time)
END FUNCTION VStartE

FUNCTION WStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: WStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  WStart=0.0d0
END FUNCTION WStart

FUNCTION TkeStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeStart=Zero
END FUNCTION TkeStart

FUNCTION DisStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DisStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DisStart=Zero
END FUNCTION DisStart

FUNCTION TkeHStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeHStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeHStart=TkeHMax
END FUNCTION TkeHStart

FUNCTION TkeVStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TkeVStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TkeVStart=TkeVMax
END FUNCTION TkeVStart

FUNCTION LenStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: LenStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  LenStart=LenMax
END FUNCTION LenStart

FUNCTION QrStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QrStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  QrStart=0.d0
END FUNCTION QrStart

FUNCTION RhoStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  USE Rho_Mod
  USE RhoProf_Mod
  IMPLICIT NONE
  REAL(RealKind) :: RhoStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  RhoStart=RhoFun(x,y,z,zHeight,Time)-RhoProf(x,y,z,zHeight,Time)
END FUNCTION RhoStart

FUNCTION TStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: TStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  TStart =0.0d0
END FUNCTION TStart

FUNCTION DummyStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: DummyStart
  REAL(RealKind) :: x,y,z,zHeight,Time
  DummyStart=Zero
END FUNCTION DummyStart

FUNCTION Tracer1Start(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer1Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  IF (z<=HeightTracer1Top.AND.z>HeightTracer1Bottom) THEN
    Tracer1Start=One
  ELSE
    Tracer1Start=Zero
  END IF
END FUNCTION Tracer1Start

FUNCTION Tracer2Start(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: Tracer2Start
  REAL(RealKind) :: x,y,z,zHeight,Time
  IF (z<=HeightTracer2Top.AND.z>HeightTracer2Bottom) THEN
    Tracer2Start=One
  ELSE
    Tracer2Start=Zero
  END IF
END FUNCTION Tracer2Start

FUNCTION HeightFun(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: HeightFun
  REAL(RealKind) :: x,y,z,zHeight,Time

  HeightFun=Zero

END FUNCTION HeightFun

FUNCTION ForceU(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceU
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceU=0.0d0
END FUNCTION ForceU

FUNCTION ForceV(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceV
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceV=0.0d0
END FUNCTION ForceV

FUNCTION ForceW(x,y,z,zHeight,Time)
  USE Kind_Mod
  REAL(RealKind) :: ForceW
  REAL(RealKind) :: x,y,z,zHeight,Time
  ForceW=0.0d0
END FUNCTION ForceW

FUNCTION ForceRho(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  USE Kind_Mod
  USE Rho_Mod
  REAL(RealKind) :: ForceRho
  REAL(RealKind) :: x,y,z,zHeight,Time
  REAL(RealKind) :: fTime,fSpace 
  ForceRho=0.0d0
END FUNCTION ForceRho

FUNCTION ThStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE CloudBoundary_Mod
  USE Start_Mod
  USE Rho_Mod
  IMPLICIT NONE
  REAL(RealKind) :: ThStartSoil,ThAir
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  REAL(RealKind) :: Rm,Cpml,KappaLoc
  REAL(RealKind) :: RhoDLoc,RhoLoc,RhoVLoc,RhoLLoc
  RhoLoc=RhoFun(x,y,z,zHeight,Time)
  RhoVLoc=QvStart(x,y,z,zHeight,Time)
  RhoLLoc=QcStart(x,y,z,zHeight,Time)
  RhoDLoc=RhoLoc-RhoVLoc-RhoLLoc+Eps
  Rm=Rd*RhoDLoc+Rv*RhoVLoc
  Cpml=Cpd*RhoDLoc+Cpv*RhoVLoc+Cpl*RhoLLoc
  KappaLoc=Rm/Cpml
  ThAir=ThStart(x,y,z,zHeight,Time)
  ThAir=(Rd*RhoLoc*ThAir/p0**KappaLoc)**(One/(One-KappaLoc))/(Rd*RhoDLoc+Rv*RhoVLoc)
  IF (zSoil<=0.05d0) THEN
!    ThStartSoil=280.50d0
    ThStartSoil=ThAir
  ELSE IF (zSoil<=0.5d0) THEN
    ThStartSoil=ThAir-0.5d0
!    ThStartSoil=280.00d0
  ELSE IF (zSoil<=1.0d0) THEN
    ThStartSoil=ThAir-2.0d0
!    ThStartSoil=279.50d0
  ELSE
    ThStartSoil=ThAir-3.5d0
!    ThStartSoil=278.00d0
  END IF
END FUNCTION ThStartSoil

FUNCTION QvStartSoil(x,y,z,zHeight,zSoil,LandClass,SoilType,Time)
  USE CloudBoundary_Mod
  IMPLICIT NONE
  REAL(RealKind) :: QvStartSoil
  REAL(RealKind) :: x,y,zSoil,z,zHeight
  INTEGER :: LandClass, SoilType
  REAL(RealKind) :: Time
  IF (zSoil==0.0d0) THEN
    QvStartSoil=0.0d0  ! interception reservoir
    QvStartSoil=MIN(QvStartSoil,5.d-4)
  ELSE 
    QvStartSoil=0.15d0
    QvStartSoil=MIN(QvStartSoil,cporv(SoilType))
  END IF
END FUNCTION QvStartSoil

FUNCTION DummyStart1(lam,phi,z,zHeight,Time)
  USE CloudBoundary_Mod
  REAL(RealKind) :: DummyStart1
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart1=0.0d0
END FUNCTION DummyStart1

FUNCTION DummyStart2(lam,phi,z,zHeight,Time)
  USE CloudBoundary_Mod
  REAL(RealKind) :: DummyStart2
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart2=0.0d0
END FUNCTION DummyStart2

FUNCTION DummyStart3(lam,phi,z,zHeight,Time)
  USE CloudBoundary_Mod
  REAL(RealKind) :: DummyStart3
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart3=0.0d0
END FUNCTION DummyStart3

FUNCTION DummyStart4(lam,phi,z,zHeight,Time)
  USE CloudBoundary_Mod
  REAL(RealKind) :: DummyStart4
  REAL(RealKind) :: lam,phi,z,zHeight,Time
  DummyStart4=0.0d0
END FUNCTION DummyStart4

FUNCTION EnStart(x,y,z,zHeight,Time)
  USE CloudBoundary_Mod
  USE ReadProfile_Mod
  IMPLICIT NONE
  REAL(RealKind) :: EnStart
  REAL(RealKind) :: x,y,z,zHeight,Time

  EnStart=Zero

END FUNCTION EnStart

