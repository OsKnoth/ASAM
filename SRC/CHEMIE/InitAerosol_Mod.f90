MODULE Distribution_Mod

  USE Kind_Mod
  USE Parameter_Mod

  IMPLICIT NONE

  INTERFACE Moment0
     MODULE PROCEDURE Moment0_Log,Moment0_Exp
  END INTERFACE 
  INTERFACE Moment1
     MODULE PROCEDURE Moment1_Log,Moment1_Exp
  END INTERFACE 
  INTERFACE Moment3
     MODULE PROCEDURE Moment3_Log
  END INTERFACE 

CONTAINS

FUNCTION Mass2Number(M,Dg,Sigma,RhoSalt)
  REAL(RealKind) :: Mass2Number
  REAL(RealKind) :: M,Dg,Sigma,RhoSalt
  
  IF (M>Zero) THEN
    Mass2Number=M/(Dg**3*EXP(Half*9.0d0*LOG(Sigma)**2)*4.0d0/3.0d0*Pi*RhoSalt)
  ELSE
    Mass2Number=Zero
  END IF
END FUNCTION Mass2Number
  

FUNCTION Num(M,DGN,SIGMAG)
  REAL(RealKind):: Num,M,DGN,SIGMAG
  IF(SIGMAG > 0.D0) THEN
    Num= (6.0d0*M/PI)/(DGN**3*EXP(4.5d0*LOG(SIGMAG)**2.0d0))
  ELSE
    Num = 0.d0
  ENDIF
END FUNCTION Num

FUNCTION MGN(DGN)
  REAL(RealKind)::MGN,DGN
  MGN=Pi*(DGN**3)/6
END FUNCTION MGN


FUNCTION SIGMAGM(SIGMAG)
  REAL(RealKind)::SIGMAGM,SIGMAG 
  SIGMAGM=SIGMAG**3
END FUNCTION SIGMAGM

FUNCTION DGN(DGV,SIGMAG)
  REAL(RealKind)::DGN,DGV,SIGMAG
  IF (SIGMAG.NE.0.0d0)THEN
     DGN=DGV/EXP(3.0d0*LOG(SIGMAG)**2.0d0)
   ELSE
     DGN=0.0d0
   ENDIF
END FUNCTION DGN

FUNCTION Nlog(x,Num,MGN,SIGMAGM)
  REAL(RealKind):: Nlog
  REAL(RealKind):: x, Num,MGN,SIGMAGM
  IF (MGN==0.0d0.OR.SIGMAGM==0.0d0.OR.x==0.0d0)THEN
     Nlog=0.0d0
  ELSE 
    Nlog=Num*EXP(-Half*(LOG(x/MGN)**2.0d0/LOG(SIGMAGM)**2.0d0))&
        /(SQRT(2.0d0*PI)*x*LOG(SIGMAGM))
   ENDIF
END FUNCTION Nlog

FUNCTION Moment0_Log(xk,xkp1,Num,MGN,SIGMAGM)
  REAL(RealKind)::Moment0_Log

  REAL(RealKind)::xk, xkp1
  REAL(RealKind)::Num, MGN, SIGMAGM

  Moment0_Log= &
        Nlog(Half*(xk+xkp1),Num,MGN,SIGMAGM)*(xkp1-xk)
 
END FUNCTION Moment0_Log

FUNCTION Moment1_Log(xk,xkp1,Num,MGN,SIGMAGM)
  REAL(RealKind)::Moment1_Log

  REAL(RealKind)::xk, xkp1
  REAL(RealKind)::Num, MGN, SIGMAGM

  Moment1_Log= &
           Nlog(Half*(xk+xkp1),Num,MGN,SIGMAGM)*Half*(xk+xkp1)*(xkp1-xk)
END FUNCTION Moment1_Log

FUNCTION Moment3_Log(xk,xkp1,Num,MGN,SIGMAGM)
  REAL(RealKind)::Moment3_Log
 
  REAL(RealKind)::xk, xkp1
  REAL(RealKind)::Num, MGN, SIGMAGM
 
  Moment3_Log= &
           Nlog(Half*(xk+xkp1),Num,MGN,SIGMAGM)*(Half*(xk+xkp1))**3*(xkp1-xk)
  Moment3_Log=Moment3_Log*4.0d0/3.0d0*Pi*1.0d3
END FUNCTION Moment3_Log

FUNCTION Moment0_Exp(xk,xkp1,N0,m0)
 
  REAL(RealKind) :: Moment0_Exp
 
  REAL(RealKind) :: xk,xkp1
  REAL(RealKind) :: N0,m0
  REAL(RealKind) :: a1,a2
 
  REAL(RealKind) :: xBar
  xBar=m0/N0
  IF (-xk/xBar<=-100.0d0) THEN
    a1=0.0d0
  ELSE
    a1=EXP(-xk/xBar)
  END IF
  IF (-xkp1/xBar<=-100.0d0) THEN
    a2=0.0d0
  ELSE
    a2=EXP(-xkp1/xBar)
  END IF
  Moment0_Exp=N0/xBar*(xkp1-xk)*Half*(a1+a2)
END FUNCTION Moment0_Exp
 
FUNCTION Moment1_Exp(xk,xkp1,N0,m0)
  REAL(RealKind) :: Moment1_Exp
  REAL(RealKind) :: xk,xkp1
  REAL(RealKind) :: N0,m0
  REAL(RealKind) :: a1,a2

  REAL(RealKind) :: Temp,xBar
  xBar=m0/N0
  IF (-xk/xBar<=-100.0d0) THEN
    a1=0.0d0
  ELSE
    a1=EXP(-xk/xBar)
  END IF
  IF (-xkp1/xBar<=-100.0d0) THEN
    a2=0.0d0
  ELSE
    a2=EXP(-xkp1/xBar)
  END IF
  Moment1_Exp=N0/xBar*(xkp1-xk)*Half*(xk*a1+xkp1*a2)
END FUNCTION Moment1_Exp

END MODULE Distribution_Mod

MODULE InitAerosol_Mod

  USE Control_Mod
  USE Parameter_Mod
  USE Physics_Mod
  USE Physics_Mod
  USE InputTool_Mod
  USE Distribution_Mod
  USE Aerosol_Mod
  USE Chemie_Mod
  USE Koehler_Mod

  IMPLICIT NONE

  TYPE AeroDistribution_T
    REAL(RealKind) :: Num,D,Sigma
  END TYPE AeroDistribution_T
  TYPE ImpactorStages_T
    INTEGER :: NumberOfStages
    INTEGER :: NumberOfSpecies
    REAL(RealKind), ALLOCATABLE :: Radius(:)
    REAL(RealKind), ALLOCATABLE :: Fraction(:,:)
  END TYPE ImpactorStages_T
  REAL(RealKind), PARAMETER :: EpsInit=1.0d-8

CONTAINS 

SUBROUTINE SetIndices

  iNC=PositionAero('aNUMBER')
  iWater=PositionAero('aH2O')
  iRelax=PositionAero('aRELAX')

END SUBROUTINE SetIndices

SUBROUTINE InitAerosol(VecT,FileName)

  TYPE(Vector4Cell_T), POINTER :: VecT(:)
  CHARACTER(*) :: FileName

  INTEGER :: ix,iy,iz
  REAL(RealKind), ALLOCATABLE :: dn_ap(:),dp_ap(:),sig_ap(:)
  INTEGER :: i
  INTEGER :: n,ns
  TYPE (Reaction_T), POINTER :: Current
  INTEGER :: nD,nSpek
  REAL(RealKind), ALLOCATABLE :: mD(:)
  REAL(RealKind), ALLOCATABLE :: NumD(:),MassD(:,:)
  REAL(RealKind), ALLOCATABLE :: Frac(:,:)
  REAL(RealKind), ALLOCATABLE :: c(:,:)
  REAL(RealKind), ALLOCATABLE :: numspek(:)
  REAL(RealKind), ALLOCATABLE :: radspek(:)
  REAL(RealKind), ALLOCATABLE :: densspek(:)
  REAL(RealKind) :: NumLoc
  REAL(RealKind) :: rMin,rMax
  REAL(RealKind) :: rL,rR,rC
  REAL(RealKind) :: epsi
  REAL(RealKind) :: eSatt, TAbs
  REAL(RealKind), POINTER :: RhoLoc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoVLoc(:,:,:,:)
  REAL(RealKind), POINTER :: RhoLLoc(:,:,:,:)
  INTEGER :: iMode,iSpek,is,iPos
  CHARACTER(30) :: SpeciesName
  REAL(RealKind) :: c1,c2,c3,c4,c5,c6,c7,c8
  LOGICAL :: Back
  REAL(RealKind) :: MassGes,ChargeGes
  INTEGER :: NumberImpactor=0
  INTEGER :: NumberModes=0
  INTEGER :: NumberSpek=0
  REAL(RealKind), ALLOCATABLE :: rMode(:)
  REAL(RealKind), ALLOCATABLE :: Work(:)
  REAL(RealKind) :: SattW
  REAL(RealKind) :: pH
  LOGICAL :: InputMass,InputNumber

  iVec1=1
  iVec2=nFrac
  ns=nAqua
  InputMass=.FALSE.
  InputNumber=.FALSE.
  ALLOCATE(LatentHeat(ns))
  LatentHeat=0.0d0
  IF (iWater>0) THEN
    LatentHeat(iWater)=lv/Cpd
  END IF

  CALL OpenFile(FileName)
  CALL LineFile(Back,'BEGIN_AERO','BEGIN_MASS','END_MASS', &
                R1=c1)
  IF (.NOT.Back) THEN
    InputMass=.TRUE.
    NumberImpactor=c1
    ALLOCATE(rMode(NumberImpactor+1))
    ALLOCATE(Frac(NumberImpactor,ns+1))
    ALLOCATE(Work(NumberImpactor))
    Frac=Zero
    CALL LineFile(Back,END='END_MASS', &
                  R=rMode)
    DO
      CALL LineFile(Back,END='END_MASS', &
                    Name1=SpeciesName, &
                    R=Work)
      IF (Back) THEN
        EXIT
      END IF
      iPos=PositionAero(SpeciesName)
      IF (iPos>0) THEN
        Frac(:,iPos)=Work 
      END IF 
    END DO
    DO iMode=1,NumberImpactor
      Frac(iMode,:)=Frac(iMode,:)/(SUM(Frac(iMode,:))+Eps)
    END DO
  END IF
  CALL CloseFile

  CALL OpenFile(FileName)
  CALL LineFile(Back,'BEGIN_AERO','BEGIN_MOL','END_MOL', &
                R1=c1)
  IF (.NOT.Back) THEN
    InputMass=.TRUE.
    NumberImpactor=c1
    ALLOCATE(rMode(NumberImpactor+1))
    ALLOCATE(Frac(NumberImpactor,ns+1))
    ALLOCATE(Work(NumberImpactor))
    Frac=Zero
    CALL LineFile(Back,END='END_MOL', &
                  R=rMode)
    DO
      CALL LineFile(Back,END='END_MOL', &
                    Name1=SpeciesName, &
                    R=Work)
      IF (Back) THEN
        EXIT
      END IF
      iPos=PositionAero(SpeciesName)
      IF (iPos>0) THEN
        Frac(:,iPos)=Work*MolMass(iPos)
      END IF
    END DO
    DO iMode=1,NumberImpactor
      Frac(iMode,:)=Frac(iMode,:)/(SUM(Frac(iMode,:))+Eps)
    END DO
  END IF
  CALL CloseFile

  CALL OpenFile(FileName)
  CALL LineFile(Back,'BEGIN_AERO','BEGIN_SPEK','END_SPEK', &
                  R1=c1)
  IF (.NOT.Back) THEN
    InputNumber=.TRUE.
    NumberSpek=c1
    ALLOCATE(radspek(NumberSpek))
    ALLOCATE(numspek(NumberSpek))
    ALLOCATE(densspek(NumberSpek))
    DO iSpek=1,NumberSpek
      CALL LineFile(Back,END='END_SPEK', &
                  R1=numspek(iSpek),        &
                  R2=radspek(iSpek),        &
                  R3=densspek(iSpek))
    END DO

    nD=NumberSpek
    ALLOCATE(mD(nD))
    ALLOCATE(NumD(nD))
    NumD=0.0d0
    ALLOCATE(MassD(nD,nAqua))
    MassD=0.0d0
    DO i=1,nD
      NumD(i)=numspek(i)
      MassD(i,1)=densspek(i)*4.0d0/3.0d0*Pi*radspek(i)**3*numspek(i)
      mD(i)=MassD(i,1)
    END DO
    iMode=1
    DO i=1,nD
      MassGes=MassD(i,1)
      ChargeGes=Zero
      NumD(i)=NumD(i)*SUM(Frac(iMode,:))
      DO is=1,nAqua
        MassD(i,is)=MassGes*Frac(iMode,is)
        ChargeGes=ChargeGes+Charge(is)*MassD(i,is)/MolMass(is)
      END DO
      IF (Neutral) THEN
        IF (ChargeGes>Zero) THEN
          MassD(i,Position('OHm'))=ChargeGes*MolMass(Position('OHm'))
        ELSE
          MassD(i,Position('Hp'))=-ChargeGes*MolMass(Position('Hp'))
        END IF
      ELSE  
        MassD(i,Position('OHm'))=EpsInit*MolMass(Position('OHm')) 
        MassD(i,Position('Hp'))=EpsInit*MolMass(Position('Hp')) 
      END IF
    END DO
    DEALLOCATE(radspek)
    DEALLOCATE(numspek)
    DEALLOCATE(densspek)
  END IF
  CALL CloseFile

  CALL OpenFile(FileName)
  CALL LineFile(Back,'BEGIN_AERO','BEGIN_MODE','END_MODE', &
                  R1=c1)
  IF (.NOT.Back) THEN
    InputNumber=.TRUE.
    NumberModes=c1
    ALLOCATE(dp_ap(NumberModes))
    ALLOCATE(sig_ap(NumberModes))
    ALLOCATE(dn_ap(NumberModes))
    DO iMode=1,NumberModes
      CALL LineFile(Back,END='END_MODE', &
                  R1=dn_ap(iMode),        & 
                  R2=dp_ap(iMode),        & 
                  R3=sig_ap(iMode)) 
      dn_ap(iMode)=dn_ap(iMode)
    END DO

    rMax=0.0d0
    rMin=1.0d20
    DO iMode=1,NumberModes
      rMax=MAX(rMax,dp_ap(iMode)*sig_ap(iMode)**6.0d0)
      rMin=MIN(rMin,dp_ap(iMode)/sig_ap(iMode)**6.0d0)
    END DO
    nD=LOG10(rMax/rMin)/LOG10(pFacInit)+1
    ALLOCATE(mD(nD))
    ALLOCATE(NumD(nD))
    NumD=0.0d0
    ALLOCATE(MassD(nD,nAqua))
    MassD=0.0d0
    rL=rMin
!  mD(1)=RhoAerosol*4.0d0/3.0d0*Pi*rL**3
    DO i=1,nD
      rR=rL*pFacInit
      rC=Half*(rL+rR)
      MassD(i,1)=Zero
      DO iMode=1,NumberModes
        NumLoc=Moment0(Two*rL,Two*rR,dn_ap(iMode),dp_ap(iMode),sig_ap(iMode)) 
        NumD(i)=NumD(i)+NumLoc
        MassD(i,1)=MassD(i,1)+RhoAerosol*4.0d0/3.0d0*Pi*rC**3*NumLoc 
      END DO
      rL=rR
      mD(i)=RhoAerosol*4.0d0/3.0d0*Pi*rC**3
    END DO
    iMode=1
    rL=rMin
    DO i=1,nD
      IF (rL>=rMode(iMode+1).AND.iMode<NumberImpactor) THEN
        iMode=iMode+1
      END IF
      rR=rL*pFacInit
      rC=Half*(rL+rR)
      MassGes=MassD(i,1) 
      ChargeGes=Zero
      NumD(i)=NumD(i)*SUM(Frac(iMode,:))  
      DO is=1,nAqua
        MassD(i,is)=MassGes*Frac(iMode,is) 
        ChargeGes=ChargeGes+Charge(is)*MassD(i,is)/MolMass(is) 
      END DO
      IF (Position('Hp')>0) THEN
        IF (Neutral.AND.Position('Hp')>0) THEN
          IF (ChargeGes>Zero) THEN
            MassD(i,Position('OHm'))=ChargeGes*MolMass(Position('OHm')) 
          ELSE
            MassD(i,Position('Hp'))=-ChargeGes*MolMass(Position('Hp')) 
          END IF
        ELSE  
          MassD(i,Position('OHm'))=EpsInit*MolMass(Position('OHm')) 
          MassD(i,Position('Hp'))=EpsInit*MolMass(Position('Hp')) 
        END IF
      END IF
      rL=rR
    END DO
    DEALLOCATE(dp_ap)
    DEALLOCATE(sig_ap)
    DEALLOCATE(dn_ap)
  END IF
  CALL CloseFile

  IF (InputMass.AND.InputNumber) THEN
    ALLOCATE(c(nFrac,ns))
    DO ibLoc=1,nbLoc
      ib=LocGlob(ibLoc)
      CALL Set(Floor(ib))
      RhoLoc=>RhoCell(ibLoc)%c
      RhoVLoc=>VecT(ibLoc)%Vec(RhoVPos)%c
      IF (RhoCpos>0) THEN
        RhoLLoc=>VecT(ibLoc)%Vec(RhoCpos)%c
      ELSE
        RhoLLoc=>RhoLCell(ibLoc)%c
      END IF
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            TAbs=TAbsCell(ibLoc)%Vec(1)%c(ix,iy,iz,1)
            eSatt=SaturVapor(TAbs)
            SattW=RhoVLoc(ix,iy,iz,1)*Rv*TAbs/eSatt-1.0d0
            IF (Position('aH2O')>0) THEN
              CALL InitDroplet(nFrac,m,c,nD,mD,NumD,MassD,SattW,TAbs)
              DO is=1,ns
                VecT(ibLoc)%Vec(is)%c(ix,iy,iz,:)=c(:,is)*VolC(ix,iy,iz)/(VolC(ix,iy,iz)+Eps)
              END DO
            END IF
          END DO
        END DO
      END DO
    END DO
  END IF  

  IF (ALLOCATED(c)) DEALLOCATE(c)
  IF (ALLOCATED(mD)) DEALLOCATE(mD)
  IF (ALLOCATED(NumD)) DEALLOCATE(NumD)
  IF (ALLOCATED(MassD)) DEALLOCATE(MassD)
  IF (ALLOCATED(Frac)) DEALLOCATE(Frac)
  IF (ALLOCATED(rMode)) DEALLOCATE(rMode)
  IF (ALLOCATED(Work)) DEALLOCATE(Work)

END SUBROUTINE InitAerosol

SUBROUTINE InitAmbientAero(VecT,FileName)

  TYPE(Vector4Cell_T), POINTER :: VecT(:)
  CHARACTER(*) :: FileName

  REAL(RealKind), ALLOCATABLE :: AmbientConc(:)
  REAL(RealKind) :: c1
  INTEGER :: i,k,iCell,Pos
  INTEGER :: ix,iy,iz
  CHARACTER(20) :: S1,S2,End,SpeciesName
  CHARACTER*300 :: Line
  LOGICAL :: Back


  S1='BEGIN_AEROSOL'
  S2='BEGIN_AMBIENT'
  End='END_AMBIENT'
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,R1=c1)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        VecAmb(ibLoc)%Vec(Pos)%c(:,:,:,1)=c1
      END DO
    END IF
  END DO
  CALL CloseFile

END SUBROUTINE InitAmbientAero

SUBROUTINE InitGrid

  INTEGER :: i
  REAL(RealKind) :: mC

  ALLOCATE(m(nFrac+1))
  ALLOCATE(r(nFrac))
  m(1)=4.0d0/3.0d0*Pi*RhoW*rStart**3
  DO i=1,nFrac
    m(i+1)=m(i)*pFac
    mC=Half*(m(i)+m(i+1))
    r(i)=(3.0d0*mC/4.0d0/Pi/RhoW)**(1.0d0/3.0d0) ! in m
  END DO 
! m(1)=0.0d0  OSSI
END SUBROUTINE InitGrid

SUBROUTINE EmiAero(AeroDistribution,ImpactorStages,Emi)

  TYPE(AeroDistribution_T) :: AeroDistribution(:)
  TYPE(ImpactorStages_T) :: ImpactorStages
  REAL(RealKInd) :: Emi(:,:)

  INTEGER :: i,iMode,iStage,iFrac,nD
  REAL(RealKind) :: D,Sigma,Num
  REAL(RealKind) :: NumLoc
  REAL(RealKind) :: rL,rR
  REAL(RealKind) :: rC
  REAL(RealKind) :: rMax,rMin
  REAL(RealKind), ALLOCATABLE :: mD(:)
  REAL(RealKind), ALLOCATABLE :: NumD(:)
  REAL(RealKind), ALLOCATABLE :: MassD(:,:)

  rMax=0.0d0
  rMin=1.0d20
  DO iMode=1,SIZE(AeroDistribution)
    D=AeroDistribution(iMode)%D
    Sigma=AeroDistribution(iMode)%Sigma
    rMax=MAX(rMax,D*Sigma**6.0d0)
    rMin=MIN(rMin,D/Sigma**6.0d0)
  END DO
  nD=LOG10(rMax/rMin)/LOG10(pFacInit)+1
  ALLOCATE(mD(nD))
  mD=0.0d0
  ALLOCATE(NumD(nD))
  NumD=0.0d0
  ALLOCATE(MassD(nD,SIZE(ImpactorStages%Fraction,2)))
  MassD=0.0d0
  rL=rMin
  DO i=1,nD
    rR=rL*pFacInit
    rC=Half*(rL+rR)
    DO iMode=1,SIZE(AeroDistribution)
      D=AeroDistribution(iMode)%D
      Sigma=AeroDistribution(iMode)%Sigma
      Num=AeroDistribution(iMode)%Num
      NumLoc=Moment0(Two*rL,Two*rR,Num,D,Sigma) 
      NumD(i)=NumD(i)+NumLoc
      mD(i)=mD(i)+RhoAerosol*4.0d0/3.0d0*Pi*rC**3*NumLoc 
    END DO
    rL=rR
  END DO
  iStage=1
  rL=rMin
  DO i=1,nD
    IF (rL>=ImpactorStages%Radius(iStage+1).AND.iStage<ImpactorStages%NumberOfStages) THEN
      iStage=iStage+1
    END IF
    rR=rL*pFacInit
    rC=Half*(rL+rR)
    MassD(i,:)=mD(i)*ImpactorStages%Fraction(iStage,:) 
  END DO  

  Emi=0.0d0
  iFrac=1
  DO i=1,nD
    IF (mD(i)<m(iFrac+1)) THEN
      Emi(1,iFrac)=Emi(1,iFrac)+NumD(i)
      Emi(2:,iFrac)=Emi(2:,iFrac)+MassD(i,:)
    ELSE IF (iFrac<nFrac) THEN
      iFrac=iFrac+1
      Emi(1,iFrac)=Emi(1,iFrac)+NumD(i)
      Emi(2:,iFrac)=Emi(2:,iFrac)+MassD(i,:)
    ELSE
      EXIT
    END IF
  END DO  
  DEALLOCATE(NumD)
  DEALLOCATE(mD)
  DEALLOCATE(MassD)

END SUBROUTINE EmiAero

END MODULE InitAerosol_Mod
