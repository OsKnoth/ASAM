MODULE Chemie_Mod

  USE Sp_Mod
  USE Physics_Mod
  USE Control_Mod
  USE Reaction_Mod
  USE Rates_Mod
  USE InputTool_Mod
  USE Transport_Mod
  USE Thermodynamic_Mod
  USE Variable_Mod

  IMPLICIT NONE

  REAL(RealKind) :: ChargeGes
  INTEGER :: nSpc,nGas,nSoluble,nPartic,nSolid,nKat,nreak
  INTEGER :: nAqua
  INTEGER, PRIVATE :: NumHenry,NumHenryPas
  INTEGER :: nReakGas    &
            ,nReakGPhoto &
            ,nReakGConst &
            ,nReakGTemp  &
            ,nReakGTroe  &
            ,nReakGSpec  &
            ,nReakHenry  &
            ,nReakDissoc &
            ,nReakAqua   &
            ,nReakAPhoto &
            ,nReakAConst &
            ,nReakATemp  &
            ,nReakASpec  &
            ,nReakSolid  & 
            ,nReakPartic & 
            ,nReakSTemp  &
            ,nReakSEqui  &
            ,nReakSSpec  &
            ,nReakMicro

  INTEGER :: nges

  TYPE (Reaction_T), POINTER :: ReactionFirst
  TYPE (Reaction_T), POINTER :: GasPhotoFirst &
                               ,GasConstFirst &
                               ,GasOtherFirst &
                               ,AquaPhotoFirst &
                               ,AquaConstFirst &
                               ,AquaOtherFirst &
                               ,DissocFirst &
                               ,HenryFirst &
                               ,HenryPasFirst &
                               ,SolidFirst &
                               ,WaterFirst &
                               ,RelaxFirst &
                               ,AdiabaticFirst &
                               ,MicroFirst

  CHARACTER*20, ALLOCATABLE :: SpeciesNameGas(:)
  CHARACTER*20, ALLOCATABLE :: SpeciesNameAqua(:)
  CHARACTER*20, ALLOCATABLE, PRIVATE :: SpeciesNameKat(:)
  REAL(RealKind), ALLOCATABLE :: MolMass(:)
  REAL(RealKind), ALLOCATABLE :: Accomod(:)
  REAL(RealKind), ALLOCATABLE :: Dg(:)
  REAL(RealKind), ALLOCATABLE :: RhoSpecies(:)
  INTEGER, ALLOCATABLE :: Charge(:)
  INTEGER, ALLOCATABLE :: ActIndex(:)

  TYPE (SpRowColDiag), TARGET, SAVE :: LUChemie
  TYPE (SpRowColDiag), TARGET, SAVE :: AChemie
! TYPE(Vec4_T) :: TVec
  INTEGER :: iVec1=15,iVec2=56 
  TYPE(Vec4_T), POINTER :: cVec(:)
  TYPE(Vec4_T), POINTER :: cVec0(:)
  TYPE(Vec4_T), POINTER :: fVec(:)
  TYPE(Vec4_T), POINTER :: AVec(:)
  TYPE(Vec4_T), POINTER :: dfVec(:)
! TYPE(Vec4_T), SAVE :: TAbs
  TYPE(Vec4_T), POINTER :: TAbs
  TYPE(ScalarCell_T), POINTER :: RhoLC ! Barthel
  INTEGER :: NumAeroOut
  INTEGER :: NumGasOut
  CHARACTER(20),ALLOCATABLE :: GasSpecies(:)
  CHARACTER(20),ALLOCATABLE :: AeroSpecies(:,:)
  INTEGER, ALLOCATABLE :: AeroOut(:)
  INTEGER, ALLOCATABLE :: GasOut(:)

  REAL(RealKind), PARAMETER :: EpsChemie=0.d-20
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_Jac
  END INTERFACE Allocate

  REAL(RealKind), PARAMETER :: Lambda=1.0d0
CONTAINS
!
!SUBROUTINE SetMetVariables
!  IF (Position('UCL')>0) THEN
!    uPosL=Position('UCL')
!  END IF
!  IF (Position('UCR')>0) THEN
!    uPosR=Position('UCR')
!  END IF
!  IF (Position('VCL')>0) THEN
!    vPosL=Position('VCL')
!  END IF
!  IF (Position('VCR')>0) THEN
!    vPosR=Position('VCR')
!  END IF
!  IF (Position('WCL')>0) THEN
!    wPosL=Position('WCL')
!  END IF
!  IF (Position('WCR')>0) THEN
!    wPosR=Position('WCR')
!  END IF
!  IF (Position('TE')>0) THEN
!    thPos=Position('TE')
!  END IF
!  IF (Position('EN')>0) THEN
!    EnPos=Position('EN')
!  END IF
!  IF (Position('TKE')>0) THEN
!    tkePos=Position('TKE')
!  END IF
!  IF (Position('DIS')>0) THEN
!    disPos=Position('DIS')
!  END IF
!  IF (Position('TKEH')>0) THEN
!    tkeHPos=Position('TKEH')
!  END IF
!  IF (Position('TKEV')>0) THEN
!    tkeVPos=Position('TKEV')
!  END IF
!  IF (Position('LEN')>0) THEN
!    LenPos=Position('LEN')
!  END IF
!  IF (Position('RHO')>0) THEN
!    rhoPos=Position('RHO')
!  END IF
!  IF (Position('PRE')>0) THEN
!    prePos=Position('PRE')
!  END IF
!  IF (Position('RHOV')>0) THEN
!    RhoVPos=Position('RHOV')
!  END IF
!  IF (Position('RHOC')>0) THEN
!    RhoCPos=Position('RHOC')
!  END IF
!  IF (Position('RHOR')>0) THEN
!    RhoRPos=Position('RHOR')
!  END IF
!  IF (Position('NV')>0) THEN
!    nvPos=Position('NV')
!  END IF
!  IF (Position('NC')>0) THEN
!    ncPos=Position('NC')
!  END IF
!  IF (Position('NR')>0) THEN
!    nrPos=Position('NR')
!  END IF
!  IF (Position('TRACER1')>0) THEN
!    tracer1Pos=Position('TRACER1')
!  END IF
!  IF (Position('TRACER2')>0) THEN
!    tracer2Pos=Position('TRACER2')
!  END IF
!END SUBROUTINE SetMetVariables
!SUBROUTINE OutputChemie(FileName)

!  CHARACTER(*) :: FileName

!  CHARACTER(20) :: S1,S2,End,SpeciesName
!  LOGICAL :: Back

!  INTEGER :: i,Pos

!  S1='BEGIN_GAS'
!  S2='BEGIN_OUTPUT'
!  End='END_OUTPUT'
!  CALL OpenFile(FileName)
!  NumGasOut=0
!  DO
!    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName)
!    IF (Back) THEN
!      EXIT
!    END IF
!    Pos=Position(SpeciesName)
!    IF (Pos>0) THEN
!      NumGasOut=NumGasOut+1
!    END IF
!  END DO
!  CALL CloseFile
!  ALLOCATE(GasOut(NumGasOut))
!  NumGasOut=0
!  CALL OpenFile(FileName)
!  DO
!    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName)
!    IF (Back) THEN
!      EXIT
!    END IF
!    Pos=Position(SpeciesName)
!    IF (Pos>0) THEN
!      NumGasOut=NumGasOut+1
!      GasOut(NumGasOut)=Pos
!    END IF
!  END DO
!  CALL CloseFile

!  S1='BEGIN_AERO'
!  S2='BEGIN_OUTPUT'
!  End='END_OUTPUT'
!  CALL OpenFile(FileName)
!  NumAeroOut=0
!  DO
!    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName)
!    IF (Back) THEN
!      EXIT
!    END IF
!    Pos=Position(SpeciesName)
!    IF (Pos>0) THEN
!      NumAeroOut=NumAeroOut+1
!    END IF
!  END DO
!  CALL CloseFile
!  ALLOCATE(AeroOut(NumAeroOut))
!  NumAeroOut=0
!  CALL OpenFile(FileName)
!  DO
!    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName)
!    IF (Back) THEN
!      EXIT
!    END IF
!    Pos=Position(SpeciesName)
!    IF (Pos>0) THEN
!      NumAeroOut=NumAeroOut+1
!      AeroOut(NumAeroOut)=Pos
!    END IF
!  END DO
!  CALL CloseFile
!END SUBROUTINE OutputChemie


SUBROUTINE AllocateVec4Chemie(Vec,VecComponents)

  TYPE(Vector4Cell_T), POINTER :: Vec(:)
  INTEGER :: VecComponents

  INTEGER :: i
  INTEGER :: Rank(2) 

  ALLOCATE(Vec(nbLoc))
  VecComponents=nGas+nAqua
  Rank(2)=1
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    ALLOCATE(Vec(ibLoc)%Vec(0:nAqua+nGas))
!   ALLOCATE(Vec(ibLoc)%Vec(nAqua+nGas))
    DO i=1,nAqua
      Rank(1)=nFrac
      CALL Allocate(Vec(ibLoc)%Vec(i),Rank)
    END DO  
    DO i=nAqua+1,nGas+nAqua
      Rank(1)=1
      CALL Allocate(Vec(ibLoc)%Vec(i),Rank)
    END DO  
    Rank(1)=1
    CALL Allocate(Vec(ibLoc)%Vec(0),Rank)
  END DO  

END SUBROUTINE AllocateVec4Chemie


SUBROUTINE Allocate_Jac(JacChemie)

  TYPE(JacSpMatrix4_T), POINTER :: JacChemie(:)

  INTEGER :: i,iAct,j,jj,jAct,n,nzr
  INTEGER :: Rank(2)
  TYPE(SpRowColDiag), POINTER :: StructLU
  TYPE(SpMatrix4Cell_T), POINTER :: JacSLU

  ALLOCATE(StructLU)
  ALLOCATE(JacChemie(nbLoc))
  CALL Jacstr_LU(StructLU)
  DO ibLoc=1,nbLoc
    ib=LocGlob(ibLoc)
    CALL Set(Floor(ib))
    ALLOCATE(JacChemie(ibLoc)%JacTMom)
    CALL SpNullify(JacChemie(ibLoc)%JacTMom)
    CALL Allocate(JacChemie(ibLoc)%JacTMom,7)
    IF (TkeSGS.OR.NoTke.OR.TkeSmag.OR.DynSmag) THEN
      ALLOCATE(JacChemie(ibLoc)%JacTPot)
      CALL SpNullify(JacChemie(ibLoc)%JacTPot)
      CALL Allocate(JacChemie(ibLoc)%JacTPot,7)
    ELSE
      JacChemie(ibLoc)%JacTPot=>JacChemie(ibLoc)%JacTMom
    END IF
    ALLOCATE(JacChemie(ibLoc)%JacFall)
    CALL SpNullify(JacChemie(ibLoc)%JacFall)
    CALL Allocate(JacChemie(ibLoc)%JacFall,2)
    ALLOCATE(JacChemie(ibLoc)%JacFallRhoL)
    CALL SpNullify(JacChemie(ibLoc)%JacFallRhoL)
    CALL Allocate(JacChemie(ibLoc)%JacFallRhoL,2)
    ALLOCATE(JacChemie(ibLoc)%JacSLU)
    JacSLU=>JacChemie(ibLoc)%JacSLU
    JacSLU%Struct=>StructLU
    n=StructLU%n
    nzr=StructLU%RowPtr(n+1)-1
    JacSLU%n=n
    JacSLU%nMaxVec=MAX(nFrac,1)
    ALLOCATE(JacSLU%Mat(nzr))
    DO i=1,n
      iAct=StructLU%InvPer(i)
      IF (TypeSpecies(iAct)=='GAS') THEN
        Rank(1)=1
      ELSE
        Rank(1)=nFrac
      END IF
      DO jj=StructLU%RowPtr(i),StructLU%RowPtr(i+1)-1
        j=StructLU%ColInd(jj)
        jAct=StructLU%InvPer(j)
        IF (TypeSpecies(jAct)=='GAS') THEN
          Rank(2)=1
        ELSE
          Rank(2)=nFrac
        END IF
        CALL ALLOCATE(JacSLU%Mat(jj),Rank)
      END DO
    END DO
  END DO

END SUBROUTINE Allocate_Jac 

FUNCTION TypeSpecies(i)

  CHARACTER*4 TypeSpecies
  INTEGER :: i
  IF (i<=nAqua) THEN
    TypeSpecies='AQUA'
  ELSE
    TypeSpecies='GAS'
  END IF

END FUNCTION TypeSpecies


FUNCTION PositionGas(Species)

  INTEGER :: PositionGas
  CHARACTER(*) :: Species

  INTEGER :: i

  PositionGas=0
  DO i=1,nGas
    IF (TRIM(Species)==TRIM(SpeciesNameGas(i))) THEN
      PositionGas=i
      EXIT
    END IF
  END DO

END FUNCTION PositionGas

FUNCTION PositionAero(Species)

  INTEGER :: PositionAero
  CHARACTER(*) :: Species

  INTEGER :: i

  PositionAero=0
  DO i=1,nAqua
    IF (TRIM(Species)==TRIM(SpeciesNameAqua(i))) THEN
      PositionAero=i
      EXIT
    END IF
  END DO

END FUNCTION PositionAero

SUBROUTINE Jacstr_LU(LU)

  TYPE (SpRowColDiag), POINTER :: LU

  TYPE (SpRowColD) :: A

  INTEGER :: i,k,l,istr,Shift
  INTEGER, ALLOCATABLE :: Struct(:)
  LOGICAL :: ins
  TYPE (Reaction_T), POINTER :: Current
  INTEGER :: PosMet(20)
  INTEGER :: iPos

! Initialisierung von  rwptrj und clindj
 

  A%n=nges
  A%m=nges
  A%len=20*nges
  CALL SpNullify(A)
  CALL Allocate(A)

! Insert diagonal
  DO i=1,nges
    CALL SpInsert(A,i,i,ins)   
  END DO
  A%Restr=0

! Reaktionen 

  Current=>ReactionFirst
  DO WHILE(ASSOCIATED(Current))
    DO i=1,Current%NumSpeciesLeftAktiv
      DO l=1,Current%NumSpecies
        CALL SpInsert(A,Current%Species(l)%Species &
                     ,Current%SpeciesLeft(i),ins)
      END DO
    END DO
    IF(TRIM(Current%ClassR)=='DISS'.OR. &
       TRIM(Current%ClassR)=='SOLID'.OR. &
       TRIM(Current%ClassR)=='HENRY'.OR. &
       TRIM(Current%ClassR)=='MICROPHYS') THEN
      DO i=1,Current%NumSpeciesRightAktiv
        DO l=1,Current%NumSpecies
          CALL SpInsert(A,Current%Species(l)%Species &
                       ,Current%SpeciesRight(i),ins)
        END DO
      END DO
    END IF
    IF (TRIM(Current%ClassR)=="HENRY".AND.Current%NumSpeciesLeftAktiv>0) THEN
      A%Restr(Current%SpeciesLeft(1))=-1
    END IF
    Current=>Current%Next
  END DO
! Restriction for Water
  IF (Position('RHOV')>0) THEN
    A%Restr(Position('RHOV'))=-1
  END IF
  IF (Position('TE')>0) THEN
    A%Restr(Position('TE'))=-1
  END IF
  IF (Position('EN')>0) THEN
    A%Restr(Position('EN'))=-1
  END IF
  IF (Position('PRE')>0) THEN
    A%Restr(Position('PRE'))=-1
  END IF
  IF (Position('RHO')>0) THEN
    A%Restr(Position('RHO'))=-1
  END IF

! Insert meteorological variables
  NumMet=0
  IF (Position('UCL')>0) THEN
    NumMet=NumMet+1
    uPosLJac=NumMet
    uPosL=Position('UCL')
    PosMet(NumMet)=uPosL
  END IF
  IF (Position('UCR')>0) THEN
    NumMet=NumMet+1
    uPosRJac=NumMet
    uPosR=Position('UCR')
    PosMet(NumMet)=uPosR
  END IF
  IF (Position('VCL')>0) THEN
    NumMet=NumMet+1
    vPosLJac=NumMet
    vPosL=Position('VCL')
    PosMet(NumMet)=vPosL
  END IF
  IF (Position('VCR')>0) THEN
    NumMet=NumMet+1
    vPosRJac=NumMet
    vPosR=Position('VCR')
    PosMet(NumMet)=vPosR
  END IF
  IF (Position('WCL')>0) THEN
    NumMet=NumMet+1
    wPosLJac=NumMet
    wPosL=Position('WCL')
    PosMet(NumMet)=wPosL
  END IF
  IF (Position('WCR')>0) THEN
    NumMet=NumMet+1
    wPosRJac=NumMet
    wPosR=Position('WCR')
    PosMet(NumMet)=wPosR
  END IF
  IF (Position('TE')>0) THEN
    NumMet=NumMet+1
    thPosJac=NumMet
    thPos=Position('TE')
    PosMet(NumMet)=thPos
  END IF
  IF (Position('EN')>0) THEN
    NumMet=NumMet+1
    EnPosJac=NumMet
    EnPos=Position('EN')
    PosMet(NumMet)=EnPos
  END IF
  IF (Position('TKE')>0) THEN
    NumMet=NumMet+1
    tkePosJac=NumMet
    tkePos=Position('TKE')
    PosMet(NumMet)=tkePos
  END IF
  IF (Position('DIS')>0) THEN
    NumMet=NumMet+1
    disPosJac=NumMet
    disPos=Position('DIS')
    PosMet(NumMet)=disPos
  END IF
  IF (Position('OME')>0) THEN
    NumMet=NumMet+1
    omePosJac=NumMet
    omePos=Position('OME')
    PosMet(NumMet)=omePos
  END IF
  IF (Position('TKEH')>0) THEN
    NumMet=NumMet+1
    tkeHPosJac=NumMet
    tkeHPos=Position('TKEH')
    PosMet(NumMet)=tkeHPos
  END IF
  IF (Position('TKEV')>0) THEN
    NumMet=NumMet+1
    tkeVPosJac=NumMet
    tkeVPos=Position('TKEV')
    PosMet(NumMet)=tkeVPos
  END IF
  IF (Position('LEN')>0) THEN
    NumMet=NumMet+1
    LenPosJac=NumMet
    LenPos=Position('LEN')
    PosMet(NumMet)=LenPos
  END IF
  IF (Position('RHO')>0) THEN
    NumMet=NumMet+1
    rhoPosJac=NumMet
    rhoPos=Position('RHO')
    PosMet(NumMet)=rhoPos
  END IF
  IF (Position('PRE')>0) THEN
    NumMet=NumMet+1
    prePosJac=NumMet
    prePos=Position('PRE')
    PosMet(NumMet)=prePos
  END IF
  IF (Position('RHOV')>0) THEN
    NumMet=NumMet+1
    RhoVPosJac=NumMet
    RhoVPos=Position('RHOV')
    PosMet(NumMet)=RhoVPos
  END IF
  IF (Position('RHOC')>0) THEN
    NumMet=NumMet+1
    RhoCPosJac=NumMet
    RhoCPos=Position('RHOC')
    PosMet(NumMet)=RhoCPos
  END IF
  IF (Position('RHOR')>0) THEN
    NumMet=NumMet+1
    RhoRPosJac=NumMet
    RhoRPos=Position('RHOR')
    PosMet(NumMet)=RhoRPos
  END IF
  IF (Position('RHOI')>0) THEN
    NumMet=NumMet+1
    RhoIPosJac=NumMet
    RhoIPos=Position('RHOI')
    PosMet(NumMet)=RhoIPos
  END IF
  IF (Position('RHOS')>0) THEN
    NumMet=NumMet+1
    RhoSPosJac=NumMet
    RhoSPos=Position('RHOS')
    PosMet(NumMet)=RhoSPos
  END IF
  IF (Position('NV')>0) THEN
    NumMet=NumMet+1
    nvPosJac=NumMet
    nvPos=Position('NV')
    PosMet(NumMet)=nvPos
  END IF
  IF (Position('NC')>0) THEN
    NumMet=NumMet+1
    ncPosJac=NumMet
    ncPos=Position('NC')
    PosMet(NumMet)=ncPos
  END IF
  IF (Position('NR')>0) THEN
    NumMet=NumMet+1
    nrPosJac=NumMet
    nrPos=Position('NR')
    PosMet(NumMet)=nrPos
  END IF
  IF (Position('NI')>0) THEN
    NumMet=NumMet+1
    niPosJac=NumMet
    niPos=Position('NI')
    PosMet(NumMet)=niPos
  END IF
  IF (Position('NS')>0) THEN
    NumMet=NumMet+1
    nsPosJac=NumMet
    nsPos=Position('NS')
    PosMet(NumMet)=nsPos
  END IF
  IF (Position('TRACER1')>0) THEN 
    NumMet=NumMet+1
    tracer1PosJac=NumMet
    tracer1Pos=Position('TRACER1')
    PosMet(NumMet)=tracer1Pos
  END IF
  IF (Position('TRACER2')>0) THEN 
    NumMet=NumMet+1
    tracer2PosJac=NumMet
    tracer2Pos=Position('TRACER2')
    PosMet(NumMet)=tracer2Pos
  END IF
  ALLOCATE(IndexMet(NumMet,NumMet))
  DO i=1,NumMet
    DO l=1,NumMet
      CALL SpInsert(A,PosMet(i)&
                     ,PosMet(l),ins)
    END DO
  END DO

! Symbolic Factorization

  CALL SymbLU(A)
  CALL SpNullify(LU)
  LU=A
  CALL Deallocate(A)

! Belegung der lokalen Anteile 

  Current=>ReactionFirst
  DO WHILE(ASSOCIATED(Current))
    istr=0
    ALLOCATE(Struct(Current%NumSpecies &
                   *(Current%NumSpeciesLeftAktiv+Current%NumSpeciesRightAktiv)))
    DO i=1,Current%NumSpeciesLeftAktiv
      DO l=1,Current%NumSpecies
        DO k=LU%RowPtr(LU%Permu(Current%Species(l)%Species)),&
             LU%RowPtr(LU%Permu(Current%Species(l)%Species)+1)-1
          IF (LU%Permu(Current%SpeciesLeft(i)).eq.LU%ColInd(k)) THEN
            istr=istr+1
            Struct(istr)=k
          END IF
        END DO
      END DO
    END DO
    IF(TRIM(Current%ClassR)=='DISS'.OR. &
       TRIM(Current%ClassR)=='SOLID'.OR. &
       TRIM(Current%ClassR)=='HENRY'.OR. &
       TRIM(Current%ClassR)=='MICROPHYS') THEN
      DO i=1,Current%NumSpeciesRightAktiv
        DO l=1,Current%NumSpecies
          DO k=LU%RowPtr(LU%Permu(Current%Species(l)%Species)),&
               LU%RowPtr(LU%Permu(Current%Species(l)%Species)+1)-1
            IF (LU%Permu(Current%SpeciesRight(i)).eq.LU%ColInd(k)) THEN
              istr=istr+1
              Struct(istr)=k
            END IF
          END DO
        END DO
      END DO
    END IF
    ALLOCATE(Current%struct(istr))
    DO i=1,istr
       Current%struct(i)=Struct(i)
    END DO
    DEALLOCATE(struct)
    Current=>Current%Next
  END DO
! Belegung der Meteorologie
  DO i=1,NumMet
    DO l=1,NumMet
      DO k=LU%RowPtr(LU%Permu(PosMet(i))),&
           LU%RowPtr(LU%Permu(PosMet(i))+1)-1
        IF (LU%Permu(PosMet(l))==LU%ColInd(k)) THEN
          IndexMet(i,l)=k
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE Jacstr_LU

SUBROUTINE SetPosition

  INTEGER :: NumMet
! Insert meteorological variables
  NumMet=0
  IF (Position('UCL')>0) THEN
    NumMet=NumMet+1
    uPosL=Position('UCL')
  END IF
  IF (Position('UCR')>0) THEN
    NumMet=NumMet+1
    uPosR=Position('UCR')
  END IF
  IF (Position('VCL')>0) THEN
    NumMet=NumMet+1
    vPosL=Position('VCL')
  END IF
  IF (Position('VCR')>0) THEN
    NumMet=NumMet+1
    vPosR=Position('VCR')
  END IF
  IF (Position('WCL')>0) THEN
    NumMet=NumMet+1
    wPosL=Position('WCL')
  END IF
  IF (Position('WCR')>0) THEN
    NumMet=NumMet+1
    wPosR=Position('WCR')
  END IF
  IF (Position('TE')>0) THEN
    NumMet=NumMet+1
    thPos=Position('TE')
  END IF
  IF (Position('EN')>0) THEN
    NumMet=NumMet+1
    EnPos=Position('EN')
  END IF
  IF (Position('TKE')>0) THEN
    NumMet=NumMet+1
    tkePos=Position('TKE')
  END IF
  IF (Position('DIS')>0) THEN
    NumMet=NumMet+1
    disPos=Position('DIS')
  END IF
  IF (Position('OME')>0) THEN
    NumMet=NumMet+1
    omePos=Position('OME')
  END IF
  IF (Position('TKEH')>0) THEN
    NumMet=NumMet+1
    tkeHPos=Position('TKEH')
  END IF
  IF (Position('TKEV')>0) THEN
    NumMet=NumMet+1
    tkeVPos=Position('TKEV')
  END IF
  IF (Position('LEN')>0) THEN
    NumMet=NumMet+1
    LenPos=Position('LEN')
  END IF
  IF (Position('RHO')>0) THEN
    NumMet=NumMet+1
    rhoPos=Position('RHO')
  END IF
  IF (Position('PRE')>0) THEN
    NumMet=NumMet+1
    prePos=Position('PRE')
  END IF
  IF (Position('RHOV')>0) THEN
    NumMet=NumMet+1
    RhoVPos=Position('RHOV')
  END IF
  IF (Position('RHOC')>0) THEN
    NumMet=NumMet+1
    RhoCPos=Position('RHOC')
  END IF
  IF (Position('RHOR')>0) THEN
    NumMet=NumMet+1
    RhoRPos=Position('RHOR')
  END IF
  IF (Position('RHOI')>0) THEN
    NumMet=NumMet+1
    RhoIPos=Position('RHOI')
  END IF
  IF (Position('RHOS')>0) THEN
    NumMet=NumMet+1
    RhoSPos=Position('RHOS')
  END IF
  IF (Position('NV')>0) THEN
    NumMet=NumMet+1
    nvPos=Position('NV')
  END IF
  IF (Position('NC')>0) THEN
    NumMet=NumMet+1
    ncPos=Position('NC')
  END IF
  IF (Position('NR')>0) THEN
    NumMet=NumMet+1
    nrPos=Position('NR')
  END IF
  IF (Position('TRACER1')>0) THEN
    NumMet=NumMet+1
    tracer1Pos=Position('TRACER1')
  END IF
  IF (Position('TRACER2')>0) THEN
    NumMet=NumMet+1
    tracer2Pos=Position('TRACER2')
  END IF
END SUBROUTINE SetPosition

SUBROUTINE InputChemicalData(FileName)

  CHARACTER(*) :: FileName

!---------------------------------------------------------------
!---  Read Chemical Data
!---------------------------------------------------------------

  INTEGER :: iPos
  REAL(RealKind) :: c1,c2,c3
  CHARACTER*300 :: Line
  CHARACTER*20 :: SpeciesName
  LOGICAL :: Back


  ALLOCATE(MolMass(nspc))
  ALLOCATE(RhoSpecies(nspc))
  MolMass=Zero
  ALLOCATE(Accomod(nGas))
  Accomod=Zero
  ALLOCATE(Dg(nGas))
  Dg=Zero
  ALLOCATE(Charge(nAqua))
  ALLOCATE(ActIndex(nAqua))
  Charge=0
  ActIndex=1


  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,'BEGIN_GAS','BEGIN_DATA','END_DATA', &
                  Name1=SpeciesName,R1=c1,R2=c2,R3=c3)
    IF (Back) THEN
      EXIT
    END IF
    iPos=PositionGas(SpeciesName)
    IF (iPos>0) THEN
      Accomod(iPos)=c2
      Dg(iPos)=c3
    END IF
    iPos=Position(SpeciesName)
    IF (iPos>0) THEN
      MolMass(iPos)=c1*1.d-3
    END IF
  END DO
  CALL CloseFile
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,'BEGIN_AERO','BEGIN_DATA','END_DATA', &
                  Name1=SpeciesName,R1=c1,R2=c2,R3=c3)
    IF (Back) THEN
      EXIT
    END IF
    iPos=Position(SpeciesName)
    IF (iPos>0) THEN
      MolMass(iPos)=c1*1.d-3
      RhoSpecies(iPos)=c3
    END IF
    iPos=PositionAero(SpeciesName)
    IF (iPos>0) THEN
      Charge(iPos)=c2
    END IF
  END DO
  CALL CloseFile

END SUBROUTINE InputChemicalData

SUBROUTINE InputSystem(FileName)

  CHARACTER(*) :: FileName

!---------------------------------------------------------------
!---  Read Chemical Reaction Mechanism
!---------------------------------------------------------------


  INTEGER :: i,i_special,j,l
  INTEGER :: Unit
  INTEGER :: iWork(10)
  TYPE (Reaction_T), POINTER :: New,Next
  TYPE (Reaction_T), POINTER :: GasPhoto &
                               ,GasConst &
                               ,GasOther &
                               ,AquaPhoto &
                               ,AquaConst &
                               ,AquaOther &
                               ,Dissoc &
                               ,Henry &
                               ,HenryPas &
                               ,Solid &
                               ,Water &
                               ,Relax &
                               ,Adiabatic &
                               ,Micro

  TYPE PointerReaction_T
    TYPE(Reaction_T), POINTER :: PoiReac
  END TYPE PointerReaction_T
  TYPE (PointerReaction_T), ALLOCATABLE :: First(:),Last(:)
  TYPE (Reaction_T), POINTER :: Current
  CHARACTER(30) :: Dummy

  NULLIFY(ReactionFirst)
  NULLIFY(GasPhotoFirst)
  NULLIFY(GasConstFirst)
  NULLIFY(GasOtherFirst)
  NULLIFY(AquaPhotoFirst)
  NULLIFY(AquaConstFirst)
  NULLIFY(AquaOtherFirst)
  NULLIFY(DissocFirst)
  NULLIFY(HenryFirst)
  NULLIFY(HenryPasFirst)
  NULLIFY(SolidFirst)
  NULLIFY(MicroFirst)
  NULLIFY(GasPhoto)
  NULLIFY(GasConst)
  NULLIFY(GasOther)
  NULLIFY(AquaPhoto)
  NULLIFY(AquaConst)
  NULLIFY(AquaOther)
  NULLIFY(Dissoc)
  NULLIFY(Henry)
  NULLIFY(HenryPas)
  NULLIFY(Solid)
  NULLIFY(Micro)

  i_special = 0
  NumHenry=0
  NumHenryPas=0
  Unit=10

  OPEN(UNIT=Unit,FILE=TRIM(FileName),STATUS='OLD')

  DO i=1,15
    READ(Unit,*)
  END DO

  READ(Unit,*) nspc
  READ(Unit,*) nGas
  READ(Unit,*) nSoluble
  READ(Unit,*) nPartic
  READ(Unit,*) nSolid
  READ(Unit,*) nKat
  nges=nspc-nkat
  nAqua=nSoluble+nSolid
  VectorComponentsT=nges

  DO i=1,3
    READ(Unit,*) 
  END DO
  NumberVariables=nspc
  ALLOCATE(Variables(NumberVariables))
  ALLOCATE(SpeciesName(nspc))
  ALLOCATE(SpeciesNameGas(nGas))
  DO i=1,nGas
    READ(Unit,*) SpeciesNameGas(i)
    SpeciesName(i+nAqua)=SpeciesNameGas(i)
    Variables(i+nAqua)%Name=SpeciesNameGas(i)
    Variables(i+nAqua)%Type='GAS'
    Variables(i+nAqua)%Unit='Mol/m^3'
  END DO
  ALLOCATE(SpeciesNameAqua(nAqua))
  DO i=1,nAqua
    READ(Unit,*) SpeciesNameAqua(i)
    SpeciesName(i)=SpeciesNameAqua(i)
    Variables(i)%Name=SpeciesNameAqua(i)
    Variables(i)%Type='AERO'
    Variables(i+nAqua)%Unit='Kg/l'
  END DO
  ALLOCATE(SpeciesNameKat(nKat))
  DO i=1,nKat
    READ(Unit,*) SpeciesNameKat(i)
    SpeciesName(i+nGas+nAqua)=SpeciesNameKat(i)
  END DO

  DO i=1,3
    READ(Unit,*) 
  END DO
  READ(Unit,*) nReak    
  READ(Unit,*) nReakGas    
  READ(Unit,*) nReakGPhoto 
  READ(Unit,*) nReakGConst 
  READ(Unit,*) nReakGTemp  
  READ(Unit,*) nReakGTroe  
  READ(Unit,*) nReakGSpec  
  READ(Unit,*) nReakHenry  
  READ(Unit,*) nReakDissoc 
  READ(Unit,*) nReakAqua   
  READ(Unit,*) nReakAPhoto 
  READ(Unit,*) nReakAConst 
  READ(Unit,*) nReakATemp  
  READ(Unit,*) nReakASpec  
  READ(Unit,*) nReakPartic  
  READ(Unit,*) nReakSolid   
  READ(Unit,*) nReakSTemp   
  READ(Unit,*) nReakSEqui   
  READ(Unit,*) nReakSSpec   
  READ(Unit,*) nReakMicro   
  READ(Unit,*)
  READ(Unit,*)
  READ(Unit,*)

  ALLOCATE(Current)
  NULLIFY(Current%next)

  DO j=1,nreak
    READ(Unit,*,END=100)
    CALL ReactionInput(Current,Unit)
    IF (Current%ClassR=='GAS') THEN
      Current%SpeciesLeft=Current%SpeciesLeft+nAqua
      Current%SpeciesRight=Current%SpeciesRight+nAqua
      Current%Species(:)%Species=Current%Species(:)%Species+nAqua
      IF (Current%TypeR(1:3)=='PHO') THEN
        IF ((.NOT.ASSOCIATED(GasPhotoFirst))) THEN 
          GasPhotoFirst=>Current
        ELSE
          GasPhoto%Next=>Current
        END IF
        GasPhoto=>Current
      ELSE IF (Current%TypeR(1:4)=='CONS') THEN
        IF ((.NOT.ASSOCIATED(GasConstFirst))) THEN 
          GasConstFirst=>Current
        ELSE
          GasConst%Next=>Current
        END IF
        GasConst=>Current
      ELSE
        IF ((.NOT.ASSOCIATED(GasOtherFirst))) THEN 
          GasOtherFirst=>Current
        ELSE
          GasOther%Next=>Current
        END IF
        GasOther=>Current
      END IF
    ELSE IF (Current%ClassR=='AQUA') THEN
      iWork(1:Current%NumSpeciesLeft)=Current%SpeciesLeft
      DEALLOCATE(Current%SpeciesLeft)
      Current%NumSpeciesLeftAktiv=Current%NumSpeciesLeftAktiv+1
      Current%NumSpeciesLeft=Current%NumSpeciesLeft+1
      ALLOCATE(Current%SpeciesLeft(Current%NumSpeciesLeft))
      Current%SpeciesLeft(1:Current%NumSpeciesLeftAktiv-1) &
         =iWork(1:Current%NumSpeciesLeftAktiv-1)
      Current%SpeciesLeft(Current%NumSpeciesLeftAktiv) &
         =Position('aH2O')
      IF (Current%NumSpeciesLeftAktiv< &
          Current%NumSpeciesLeft) THEN
        Current%SpeciesLeft(Current%NumSpeciesLeftAktiv+1:) &
           =iWork(Current%NumSpeciesLeftAktiv &
                 :Current%NumSpeciesLeft-1)
      END IF
      IF (Current%TypeR(1:3)=='PHO') THEN
        IF ((.NOT.ASSOCIATED(AquaPhotoFirst))) THEN
          AquaPhotoFirst=>Current
        ELSE
          AquaPhoto%Next=>Current
        END IF
        AquaPhoto=>Current
      ELSE IF (Current%TypeR(1:4)=='CONS') THEN
        IF ((.NOT.ASSOCIATED(AquaConstFirst))) THEN
          AquaConstFirst=>Current
        ELSE
          AquaConst%Next=>Current
        END IF
        AquaConst=>Current
      ELSE
        IF ((.NOT.ASSOCIATED(AquaOtherFirst))) THEN
          AquaOtherFirst=>Current
        ELSE
          AquaOther%Next=>Current
        END IF
        AquaOther=>Current
      END IF
    ELSE IF (Current%ClassR=='HENRY') THEN
      iWork(1:Current%NumSpeciesRight)=Current%SpeciesRight
      DEALLOCATE(Current%SpeciesRight)
      Current%NumSpeciesRightAktiv=Current%NumSpeciesRightAktiv+2
      Current%NumSpeciesRight=Current%NumSpeciesRight+2
      ALLOCATE(Current%SpeciesRight(Current%NumSpeciesRight))
      Current%SpeciesRight(1:Current%NumSpeciesRightAktiv-2) &
         =iWork(1:Current%NumSpeciesRightAktiv-2)
      Current%SpeciesRight(Current%NumSpeciesRightAktiv-1) &
         =Position('aNUMBER')
      Current%SpeciesRight(Current%NumSpeciesRightAktiv) &
         =Position('aH2O')
      IF (Current%NumSpeciesRightAktiv< &
          Current%NumSpeciesRight) THEN
        Current%SpeciesRight(Current%NumSpeciesRightAktiv+1:) &
           =iWork(Current%NumSpeciesRightAktiv-1 &
                 :Current%NumSpeciesRight-2)
      END IF
      Current%SpeciesLeft=Current%SpeciesLeft+nAqua
      IF (Current%NumSpeciesRightAktiv/=0) THEN
        IF (Current%NumSpeciesLeftAktiv/=0) THEN
          Current%Species(1)%Species=Current%Species(1)%Species+nAqua
        END IF
      END IF
      IF ((Current%NumSpeciesLeftAktiv==0).OR.(Current%NumSpeciesRightAktiv==0)) THEN
        NumHenryPas=NumHenryPas+1
        IF ((.NOT.ASSOCIATED(HenryPasFirst))) THEN
          HenryPasFirst=>Current
        ELSE
          HenryPas%Next=>Current
        END IF
        HenryPas=>Current
      ELSE
        NumHenry=NumHenry+1
        IF ((.NOT.ASSOCIATED(HenryFirst))) THEN
          HenryFirst=>Current
        ELSE
          Henry%Next=>Current
        END IF
        Henry=>Current
      END IF
    ELSE IF (Current%ClassR=='DISS') THEN
      iWork(1:Current%NumSpeciesRight)=Current%SpeciesRight
      DEALLOCATE(Current%SpeciesRight)
      Current%NumSpeciesRightAktiv=Current%NumSpeciesRightAktiv+1
      Current%NumSpeciesRight=Current%NumSpeciesRight+1
      ALLOCATE(Current%SpeciesRight(Current%NumSpeciesRight))
      Current%SpeciesRight(1:Current%NumSpeciesRightAktiv-1) &
         =iWork(1:Current%NumSpeciesRightAktiv-1)
      Current%SpeciesRight(Current%NumSpeciesRightAktiv) &
         =Position('aH2O')
      IF (Current%NumSpeciesRightAktiv< &
          Current%NumSpeciesRight) THEN
        Current%SpeciesRight(Current%NumSpeciesRightAktiv+1:) &
           =iWork(Current%NumSpeciesRightAktiv &
                 :Current%NumSpeciesRight-1)
      END IF
      IF ((.NOT.ASSOCIATED(DissocFirst))) THEN
        DissocFirst=>Current
      ELSE
        Dissoc%Next=>Current
      END IF
      Dissoc=>Current
    ELSE IF (Current%ClassR=='SOLID') THEN
      IF (Current%NumSpeciesLeftAktiv/=0) THEN
!       Current%Species(1)%Species=Current%Species(1)%Species+nSoluble
!       Current%SpeciesLeft=Current%SpeciesLeft+nSoluble
      END IF
      iWork(1:Current%NumSpeciesRight)=Current%SpeciesRight
      DEALLOCATE(Current%SpeciesRight)
      Current%NumSpeciesRightAktiv=Current%NumSpeciesRightAktiv+1
      Current%NumSpeciesRight=Current%NumSpeciesRight+1
      ALLOCATE(Current%SpeciesRight(Current%NumSpeciesRight))
      Current%SpeciesRight(1:Current%NumSpeciesRightAktiv-1) &
         =iWork(1:Current%NumSpeciesRightAktiv-1)
      Current%SpeciesRight(Current%NumSpeciesRightAktiv) &
         =Position('aH2O')
      IF (Current%NumSpeciesRightAktiv< &
          Current%NumSpeciesRight) THEN
        Current%SpeciesRight(Current%NumSpeciesRightAktiv+1:) &
           =iWork(Current%NumSpeciesRightAktiv &
                 :Current%NumSpeciesRight-1)
      END IF
      IF ((.NOT.ASSOCIATED(SolidFirst))) THEN
        SolidFirst=>Current
      ELSE
        Solid%Next=>Current
      END IF
      Solid=>Current
    ELSE IF (Current%ClassR=='MICROPHYS') THEN
      IF ((.NOT.ASSOCIATED(MicroFirst))) THEN
        MicroFirst=>Current
      ELSE
        Micro%Next=>Current
      END IF
      Micro=>Current
    END IF

    IF (j<nreak)  THEN
      ALLOCATE(New)
      NULLIFY(New%next)
      Current=>New
    END IF
  END DO 

! Append waterflux
  IF (Position('aH2O')>0) THEN
    IF (Position('RHOV')==0) THEN
      STOP 'No RHOV in .chem-file present, but necessary'
    ELSE IF (Position('TE')==0) THEN
      STOP 'No TE in .chem-file present, but necessary'
    END IF

    ALLOCATE(Current)
    NULLIFY(Current%Next)
    WaterFirst=>Current
    Water=>Current
    Current%ClassR='WATER'
    Current%NumSpeciesLeft=nAqua+2
    Current%NumSpeciesLeftAktiv=nAqua+2
    Current%NumSpeciesRight=0
    Current%NumSpeciesRightAktiv=0
    ALLOCATE(Current%SpeciesLeft(Current%NumSpeciesLeft))
    DO i=1,nAqua
      Current%SpeciesLeft(i)=i
    END DO
    Current%SpeciesLeft(nAqua+1)=Position('RHOV')
    Current%SpeciesLeft(nAqua+2)=Position('TE')
!   Current%SpeciesLeft(nAqua+3)=Position('PRE')
    Current%NumSpecies=1
    ALLOCATE(Current%Species(Current%NumSpecies))
    Current%Species(:)%Koeff=RelFac1  !OSSI
    Current%Species(1)%Species=Position('aRELAX')

    ALLOCATE(Current)
    NULLIFY(Current%Next)
    Water%Next=>Current
    Water=>Current
    Current%ClassR='WATER'
    Current%NumSpeciesLeft=2
    Current%NumSpeciesLeftAktiv=2
    Current%NumSpeciesRight=0
    Current%NumSpeciesRightAktiv=0
    ALLOCATE(Current%SpeciesLeft(Current%NumSpeciesLeft))
    Current%SpeciesLeft(1)=Position('aNUMBER')
    Current%SpeciesLeft(2)=Position('aRELAX')
    Current%NumSpecies=3
    ALLOCATE(Current%Species(Current%NumSpecies))
    Current%Species(1)%Koeff=1.0d0
    Current%Species(1)%Species=Position('aH2O')
    Current%Species(2)%Koeff=-1.0d0
    Current%Species(2)%Species=Position('RHOV')
    Current%Species(3)%Koeff=lv/Cpd
!   Current%Species(3)%Koeff=Zero !OSSI
    Current%Species(3)%Species=Position('TE')
  END IF

! Append Relax
  IF (Position('aRELAX')>0) THEN
    ALLOCATE(Current)
    NULLIFY(Current%Next)
    RelaxFirst=>Current
    Relax=>Current
    Current%ClassR='RELA'
    Current%NumSpeciesLeft=1
    Current%NumSpeciesLeftAktiv=1
    Current%NumSpeciesRight=0
    Current%NumSpeciesRightAktiv=0
    ALLOCATE(Current%SpeciesLeft(Current%NumSpeciesLeft))
    Current%SpeciesLeft(1)=Position('aRELAX')
    Current%NumSpecies=1
    ALLOCATE(Current%Species(Current%NumSpecies))
    Current%Species(:)%Koeff=-1.0d0
    Current%Species(1)%Species=Position('aRELAX')
  END IF

! Append Adiabatic
  IF (Parcel) THEN
    ALLOCATE(Current)
    NULLIFY(Current%Next)
    AdiabaticFirst=>Current
    Adiabatic=>Current
    Current%ClassR='HEAT'
    Current%NumSpeciesLeft=3
    Current%NumSpeciesLeftAktiv=3
    Current%NumSpeciesRight=0
    Current%NumSpeciesRightAktiv=0
    ALLOCATE(Current%SpeciesLeft(Current%NumSpeciesLeft))
    Current%SpeciesLeft(1)=Position('RHOV')
    Current%SpeciesLeft(2)=Position('TE')
    Current%SpeciesLeft(3)=Position('RHO')
    Current%NumSpecies=1
    ALLOCATE(Current%Species(Current%NumSpecies))
    Current%Species(:)%Koeff=1.0d0
    Current%Species(1)%Species=Position('TE')
  END IF

100 continue
! CALL AllocateStack_calculate

! Create whole list
  ALLOCATE(First(NumReactionTypes))
  ALLOCATE(Last(NumReactionTypes))
  First(1)%PoiReac=>GasPhotoFirst
  First(2)%PoiReac=>GasConstFirst
  First(3)%PoiReac=>GasOtherFirst
  First(4)%PoiReac=>AquaPhotoFirst
  First(5)%PoiReac=>AquaConstFirst
  First(6)%PoiReac=>AquaOtherFirst
  First(7)%PoiReac=>DissocFirst
  First(8)%PoiReac=>HenryFirst
  First(9)%PoiReac=>HenryPasFirst
  First(10)%PoiReac=>SolidFirst  
  First(11)%PoiReac=>WaterFirst  
  First(12)%PoiReac=>AdiabaticFirst  
  First(13)%PoiReac=>RelaxFirst  
  First(14)%PoiReac=>MicroFirst  

  Last(1)%PoiReac=>GasPhoto
  Last(2)%PoiReac=>GasConst
  Last(3)%PoiReac=>GasOther
  Last(4)%PoiReac=>AquaPhoto
  Last(5)%PoiReac=>AquaConst
  Last(6)%PoiReac=>AquaOther
  Last(7)%PoiReac=>Dissoc
  Last(8)%PoiReac=>Henry
  Last(9)%PoiReac=>HenryPas
  Last(10)%PoiReac=>Solid  
  Last(11)%PoiReac=>Water  
  Last(12)%PoiReac=>Adiabatic  
  Last(13)%PoiReac=>Relax  
  Last(14)%PoiReac=>Micro  

  DO i=1,NumReactionTypes
    IF (ASSOCIATED(First(i)%PoiReac)) THEN
      IF (.NOT.ASSOCIATED(ReactionFirst)) THEN
        ReactionFirst=>First(i)%PoiReac
      ELSE
        Next%Next=>First(i)%PoiReac
      END IF
      Next=>Last(i)%PoiReac
    END IF
  END DO
END SUBROUTINE InputSystem

SUBROUTINE GasChemie(cGas,f,Temp)

  TYPE (Vec4_T), TARGET :: cGas(0:)
  TYPE (Vec4_T) :: f(0:)
  TYPE (Vec4_T), TARGET :: Temp

  INTEGER :: l,iTAbs,iReak,nreakGLoc
  TYPE (Vec4_T) :: w 
  INTEGER :: iReac,i
  TYPE (Reaction_T), POINTER :: Reaction
  TYPE (Vec4_T), POINTER :: Pres
  TYPE (Vec4_T) :: mAir, mH2O
  TYPE (Vec4_T), POINTER :: rHum
  REAL(RealKind), POINTER :: Constants(:)

  cVecReac=>cGas
  Temperature=>Temp
  rHum=>cGas(RhoVPos)

  CALL Allocate(w,Temp)
  CALL Allocate(mAir,Temp)
  CALL Allocate(mH2O,Temp)
  CALL Init_StackVector(Temp)
  Chi=acos(cosSun)
  !mAir%c=2.46E19*298.15E0*Rd*RhoC%c/1.01325d5*ConvAir ! [mol m-3]
  !Hier für mol/m^3  und molec/cm^3 Faktoren vorab berechnen in InitChemie
  mAir%c=2.46E19/Temp%c*298.15E0*PreCell(ibLoc)%c*ConvAir
  mH2O%c=rHum%c/18.01528d-3*Avga/1.d6 ! [molec cm-3], Added by Tian
  IF (chi<=PiHalf.AND.ASSOCIATED(GasPhotoFirst)) THEN
    Reaction=>GasPhotoFirst
    nreakGLoc=nReakGas
  ELSE IF(ASSOCIATED(GasConstFirst)) THEN
    Reaction=>GasConstFirst
    nreakGLoc=nReakGas-nReakGPhoto
  ELSE
    Reaction=>GasOtherFirst
    nreakGLoc=nReakGas-nReakGPhoto-nReakGConst
  END IF
  DO iReak=1,nreakGLoc
    Constants=>Reaction%Constants(1,1,1,:)
    SELECT CASE (Reaction%TypeR)
      CASE ('PHOTABC')
        CALL PhoABCCompute(w,Constants)
      CASE ('PHOTAB')
        CALL PhoABCompute(w,Constants)
      CASE ('PHOTMCM')
        CALL PhoMCMCompute(w,Constants)
      CASE ('PHOTMCM1')
        CALL PhoMCM1Compute(w,Constants,Temp)
      CASE ('PHOTMCM_S')
        CALL PhoMCMSCompute(w,Constants)
      CASE ('SPECIAL')
        CALL Calculate(w,Reaction%Infix,Reaction%Constants)
      CASE ('CONST')
        CALL ConstCompute(w,Constants)
      CASE ('CONST_S')
        CALL ConstSCompute(w,Constants)
      CASE ('TEMP1')
        CALL Temp1Compute(w,Constants,Temp)
      CASE ('TEMP2')
        CALL Temp2Compute(w,Constants,Temp)
      CASE ('TEMP3')
        CALL Temp3Compute(w,Constants,Temp)
      CASE ('TROE') 
        CALL TroeCompute(w,Constants,Temp,mAir)
      CASE ('TROEQ') 
        CALL TroeEqCompute(w,Constants,Temp,mAir)
      CASE ('TROEF') ! Barthel
        CALL TroeFCompute(w,Constants,Temp,mAir)
      CASE ('TROEXP') ! Barthel
        CALL TroeXPCompute(w,Constants,Temp,mAir)
      CASE ('TROEQF') ! Barthel
        CALL TroeEqfCompute(w,Constants,Temp,mAir)
      CASE ('SPEC1') 
        CALL Spec1Compute(w,Constants,mAir)
      CASE ('SPEC2') 
        CALL Spec2Compute(w,Constants,Temp,mAir)
      CASE ('SPEC3') 
        CALL Spec3Compute(w,Constants,Temp,mAir)
      CASE ('SPEC4') 
        CALL Spec4Compute(w,Constants,Temp,mAir)
      CASE ('S4H2O') 
        CALL S4H2OCompute(w,Constants,Temp,mAir,rHum)
      CASE DEFAULT
        w=Zero
        STOP 'Falsche ReaktionsTyp'
    END SELECT
    IF (ASSOCIATED(Reaction%Factor)) THEN  
      SELECT CASE (Reaction%Factor)
        CASE ('$O2')               
          w%c=w%c*((mO2*mair%c)**Reaction%fac_exp)*Reaction%fac_A
        CASE ('$H2O')                 
          w%c=w%c*((mH2O%c)**Reaction%fac_exp)*Reaction%fac_A
        CASE ('$N2')                 
          w%c=w%c*((mN2*mair%c)**Reaction%fac_exp)*Reaction%fac_A
        CASE ('$H2')                 
          w%c=w%c*((mH2*mair%c)**Reaction%fac_exp)*Reaction%fac_A
        CASE DEFAULT
          WRITE(*,*) 'GasPhase: WRONG FACTOR FOR PASSIVE SPECIES ',Reaction%Factor
          STOP       'GasPhase: WRONG FACTOR FOR PASSIVE SPECIES'
      END SELECT
    END IF
    DO l=1,Reaction%NumSpeciesLeftAktiv
      CALL Mult(cGas(Reaction%SpeciesLeft(l)),w)
    END DO
    DO l=1,Reaction%NumSpecies
      CALL Add(Reaction%Species(l)%Koeff,w,f(Reaction%Species(l)%Species))
    END DO
    Reaction=>Reaction%next
  END DO
  CALL Deallocate(w)
  CALL Deallocate(mAir)
  CALL Deallocate(mH2O)
  CALL Close_StackVector
END SUBROUTINE GasChemie

SUBROUTINE GasChemieJac(cGas,Jac,Temp)

  TYPE(Vec4_T), TARGET :: cGas(0:)
  TYPE(Vec4_T), POINTER :: Jac(:)
  TYPE(Vec4_T), TARGET :: Temp

  INTEGER :: i,idf,istr,l
  INTEGER :: iReak,nreakGLoc
  TYPE(Vec4_T) :: w,Work
  INTEGER :: iReac
  TYPE (Reaction_T), POINTER :: Reaction
  TYPE (Vec4_T) Pres
  TYPE (Vec4_T) :: mAir, mH2O
  TYPE (Vec4_T), POINTER :: rHum
  REAL(RealKind), POINTER :: Constants(:)

  cVecReac=>cGas
  Temperature=>Temp
!  Pres=>cGas(POSITION('PRE'))
! Pres=>PreCell(ibLoc)
  rHum=>cGas(RhoVPos)

!  IF (Position('H2O')>0) THEN ! Barthel
!    cGas(Position('H2O'))%c=rHum%c/Molmass(Position('H2O'))
!  ENDIF

  CALL Allocate(w,Temp)
  CALL Allocate(Work,Temp)
  CALL Allocate(mAir,Temp)
  CALL Allocate(mH2O,Temp)
  CALL Init_StackVector(Temp)

  Chi=acos(cosSun)
  mAir%c=2.46d19/Temp%c*298.15d0*PreCell(ibLoc)%c*ConvAir
! mAir%c=2.46E19*298.15E0*Rd*RhoC%c/1.01325d5*ConvAir ! [mol m-3]
  mH2O%c=rHum%c/18.01528d-3 ! [mol m-3], Added by Tian
  IF (chi<=PiHalf.AND.ASSOCIATED(GasPhotoFirst)) THEN
    Reaction=>GasPhotoFirst
    nreakGLoc=nReakGas
  ELSE IF(ASSOCIATED(GasConstFirst)) THEN
    Reaction=>GasConstFirst
    nreakGLoc=nReakGas-nReakGPhoto
  ELSE
    Reaction=>GasOtherFirst
    nreakGLoc=nReakGas-nReakGPhoto-nReakGConst
  END IF
  DO iReak=1,nreakGLoc
    Constants=>Reaction%Constants(1,1,1,:)
    SELECT CASE (Reaction%TypeR)
      CASE ('PHOTABC')
        CALL PhoABCCompute(w,Constants)
      CASE ('PHOTAB')
        CALL PhoABCompute(w,Constants)
      CASE ('PHOTMCM')
        CALL PhoMCMCompute(w,Constants)
      CASE ('PHOTMCM1')
        CALL PhoMCM1Compute(w,Constants,Temp)
      CASE ('PHOTMCM_S')
        CALL PhoMCMSCompute(w,Constants)
      CASE ('SPECIAL')
        CALL Calculate(w,Reaction%Infix,Reaction%Constants)
      CASE ('CONST')
        CALL ConstCompute(w,Constants)
      CASE ('CONST_S')
        CALL ConstSCompute(w,Constants)
      CASE ('TEMP1')
        CALL Temp1Compute(w,Constants,Temp)
      CASE ('TEMP2')
        CALL Temp2Compute(w,Constants,Temp)
      CASE ('TEMP3')
        CALL Temp3Compute(w,Constants,Temp)
      CASE ('TROE') 
        CALL TroeCompute(w,Constants,Temp,mAir)
      CASE ('TROEQ') 
        CALL TroeEqCompute(w,Constants,Temp,mAir)
      CASE ('TROEF') ! Barthel
        CALL TroeFCompute(w,Constants,Temp,mAir)
      CASE ('TROEXP') ! Barthel
        CALL TroeXPCompute(w,Constants,Temp,mAir)
      CASE ('TROEQF') ! Barthel
        CALL TroeEqfCompute(w,Constants,Temp,mAir)
      CASE ('SPEC1') 
        CALL Spec1Compute(w,Constants,mAir)
      CASE ('SPEC2') 
        CALL Spec2Compute(w,Constants,Temp,mAir)
      CASE ('SPEC3') 
        CALL Spec3Compute(w,Constants,Temp,mAir)
      CASE ('SPEC4') 
        CALL Spec4Compute(w,Constants,Temp,mAir)
      CASE ('S4H2O') 
        CALL S4H2OCompute(w,Constants,Temp,mAir,rHum)
      CASE DEFAULT
        WRITE(*,*) 'STOP Falscher ReaktionsTyp',Reaction%TypeR
        w=Zero
        STOP 
    END SELECT
    IF (ASSOCIATED(Reaction%Factor)) THEN  
      SELECT CASE (Reaction%Factor)
        CASE ('$O2')               
          w%c=w%c*((mO2*mair%c)**Reaction%fac_exp)*Reaction%fac_A
        CASE ('$H2O')                 
          w%c=w%c*((mH2O%c)**Reaction%fac_exp)*Reaction%fac_A
          ! Rate(:,1) = Rate(:,1)*(mH2O(:)**Reaction%fac_exp)*Reaction%fac_A
        CASE ('$N2')                 
          w%c=w%c*((mN2*mair%c)**Reaction%fac_exp)*Reaction%fac_A
        CASE ('$H2')                 
          w%c=w%c*((mH2*mair%c)**Reaction%fac_exp)*Reaction%fac_A
        CASE DEFAULT
          WRITE(*,*) 'GasPhase: WRONG FACTOR FOR PASSIVE SPECIES ',Reaction%Factor
          STOP       'GasPhase: WRONG FACTOR FOR PASSIVE SPECIES'
      END SELECT
    END IF
    istr=0
    DO i=1,Reaction%NumSpeciesLeftAktiv
      Work=w
      DO l=1,Reaction%NumSpeciesLeftAktiv
        IF (l/=i) THEN
          CALL Mult(cGas(Reaction%SpeciesLeft(l)),Work)
        END IF
      END DO
      DO l=1,Reaction%NumSpecies
        istr=istr+1
        idf=Reaction%Struct(istr)
        CALL Add(Reaction%Species(l)%Koeff,Work,Jac(idf))
      END DO
    END DO
    Reaction=>Reaction%next
  END DO
  CALL Deallocate(w)
  CALL Deallocate(mAir)
  CALL Deallocate(mH2O)
  CALL Deallocate(Work)
  CALL Close_StackVector
END SUBROUTINE GasChemieJac

SUBROUTINE AquaChemie(cAqua,f,Temp,Act)

  TYPE (Vec4_T), TARGET :: cAqua(0:)
  TYPE (Vec4_T) :: f(0:)
  TYPE (Vec4_T), TARGET :: Temp
  TYPE (Vec4_T) :: Act(0:)

  INTEGER :: l,is,isW
  INTEGER :: iReak,nreakALoc
  TYPE (Vec4_T) :: w
  REAL(RealKind) :: Koeff
  INTEGER :: iReac
  TYPE (Reaction_T), POINTER :: Reaction
  REAL(RealKind), POINTER :: Constants(:)

  cVecReac=>cAqua
  Temperature=>Temp
  isW=Position('aH2O')
  CALL Allocate(w,cAqua(1))
  CALL Init_StackVector(cAqua(1))
  Chi=acos(cosSun)
  IF (chi<=PiHalf.AND.ASSOCIATED(AquaPhotoFirst)) THEN
    Reaction=>AquaPhotoFirst
    nreakALoc=nReakAqua
  ELSE IF(ASSOCIATED(AquaConstFirst)) THEN
    Reaction=>AquaConstFirst
    nreakALoc=nReakAqua-nReakAPhoto
  ELSE
    Reaction=>AquaOtherFirst
    nreakALoc=nReakAqua-nReakAPhoto-nReakAConst
  END IF
  DO  iReak=1,nreakALoc
    Constants=>Reaction%Constants(1,1,1,:)
    SELECT CASE (Reaction%TypeR)
      CASE ('PHOTABC')
        CALL PhoABCCompute(w,Constants)
      CASE ('PHOTAB')
        CALL PhoABCompute(w,Constants)
      CASE ('PHOTMCM')
        CALL PhoMCMCompute(w,Constants)
      CASE ('SPECIAL')
        CALL Calculate(w,Reaction%Infix,Reaction%Constants)
      CASE ('CONST')
        CALL ConstCompute(w,Constants)
      CASE ('TEMP1')
        CALL Temp1Compute(w,Constants,Temp)
      CASE ('TEMP2')
        CALL Temp2Compute(w,Constants,Temp)
      CASE ('TEMP3')
        CALL Temp3Compute(w,Constants,Temp)
      CASE DEFAULT
        w=Zero
        STOP 'Falscher ReaktionsTyp in AquaChemie'
    END SELECT
    DO l=1,Reaction%NumSpeciesLeft-1
      is=Reaction%SpeciesLeft(l)
      CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),w)
    END DO
    CALL Mult(cAqua(isW),w) ! w wird zu mol/m^3
    DO l=1,Reaction%NumSpecies
      is=Reaction%Species(l)%Species
      Koeff=Reaction%Species(l)%Koeff
      CALL Add(Koeff*MolMass(is),w,f(is),iVec1,iVec2) ! wenn Koeff einheitenlos dann f in kg/m^3
    END DO
    Reaction=>Reaction%next
  END DO
  CALL Deallocate(w)
  CALL Close_StackVector
END SUBROUTINE AquaChemie

SUBROUTINE AquaChemieJac(cAqua,Jac,Temp,Act)

  TYPE (Vec4_T), TARGET :: cAqua(0:)
  TYPE(Vec4_T), POINTER :: Jac(:)
  TYPE (Vec4_T), TARGET :: Temp
  TYPE (Vec4_T) :: Act(0:)

  INTEGER :: i,is,isD,isW,idf,istr,l,iTAbs
  INTEGER :: iReak,nreakALoc
  REAL(RealKind) :: Fac,Koeff
  TYPE (Vec4_T) :: w,Work
  INTEGER :: iReac
  TYPE (Reaction_T), POINTER :: Reaction
  REAL(RealKind), POINTER :: Constants(:)

  cVecReac=>cAqua
  Temperature=>Temp
  isW=Position('aH2O')
  iTAbs=PositionGas('TE')
  CALL Allocate(w,cAqua(1))
  CALL Allocate(Work,cAqua(1))
  CALL Init_StackVector(cAqua(1))
  Chi=acos(cosSun)
  IF (chi<=PiHalf.AND.ASSOCIATED(AquaPhotoFirst)) THEN
    Reaction=>AquaPhotoFirst
    nreakALoc=nReakAqua
  ELSE IF(ASSOCIATED(AquaConstFirst)) THEN
    Reaction=>AquaConstFirst
    nreakALoc=nReakAqua-nReakAPhoto
  ELSE
    Reaction=>AquaOtherFirst
    nreakALoc=nReakAqua-nReakAPhoto-nReakAConst
  END IF
  DO  iReak=1,nreakALoc
    Constants=>Reaction%Constants(1,1,1,:)
    SELECT CASE (Reaction%TypeR)
      CASE ('PHOTABC')
        CALL PhoABCCompute(w,Constants)
      CASE ('PHOTAB')
        CALL PhoABCompute(w,Constants)
      CASE ('PHOTMCM')
        CALL PhoMCMCompute(w,Constants)
      CASE ('SPECIAL')
        CALL Calculate(w,Reaction%Infix,Reaction%Constants)
      CASE ('CONST')
        CALL ConstCompute(w,Constants)
      CASE ('TEMP1')
        CALL Temp1Compute(w,Constants,Temp)
      CASE ('TEMP2')
        CALL Temp2Compute(w,Constants,Temp)
      CASE ('TEMP3')
        CALL Temp3Compute(w,Constants,Temp)
      CASE DEFAULT
        w=Zero
        STOP 'Falscher ReaktionsTyp in AquaChemieJac'
    END SELECT
    istr=0
    DO i=1,Reaction%NumSpeciesLeftAktiv-1
      isD=Reaction%SpeciesLeft(i)
      Work=w
      CALL Mult(Act(isD),Work)
      DO l=1,Reaction%NumSpeciesLeft-1
        IF (l/=i) THEN
          is=Reaction%SpeciesLeft(l)
          CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),Work)
        END IF
      END DO
      DO l=1,Reaction%NumSpecies
        is=Reaction%Species(l)%Species
        Koeff=Reaction%Species(l)%Koeff
        istr=istr+1
        idf=Reaction%Struct(istr)
        CALL Add(Koeff*MolMass(is)/MolMass(isD),Work,Jac(idf),iVec1,iVec2)
      END DO
    END DO
    is=Reaction%SpeciesLeft(Reaction%NumSpeciesLeftAktiv-1)
    CALL Mult(cAqua(is),MolMass(is),cAqua(isW),Work)
    Fac=-(Reaction%NumSpeciesLeftAktiv-2)
    CALL ScaleV(Fac,Work)
    DO l=1,Reaction%NumSpecies
      is=Reaction%Species(l)%Species
      Koeff=Reaction%Species(l)%Koeff
      istr=istr+1
      idf=Reaction%Struct(istr)
      CALL Add(Koeff*MolMass(is),Work,Jac(idf),iVec1,iVec2)
    END DO
    Reaction=>Reaction%next
  END DO
  CALL Deallocate(w)
  CALL Deallocate(Work)
  CALL Close_StackVector
END SUBROUTINE AquaChemieJac


SUBROUTINE  DissChemie(cAqua,f,Temp,Act)

  TYPE (Vec4_T), TARGET :: cAqua(0:)
  TYPE (Vec4_T) :: f(0:)
  TYPE (Vec4_T), TARGET :: Temp
  TYPE (Vec4_T) :: Act(0:)

  INTEGER :: l,il,is,isW,iTAbs
  INTEGER :: iReak
  REAL(RealKind) :: Koeff
  TYPE (Reaction_T), POINTER :: Reaction
  TYPE(Vec4_T) :: w,v
  REAL(RealKind), POINTER :: Constants(:)

  cVecReac=>cAqua
  Temperature=>Temp
  isW=Position('aH2O')
  CALL Allocate(w,cAqua(1))
  CALL Allocate(v,cAqua(1))
  Reaction=>DissocFirst
  DO iReak=1,nReakDissoc
    Constants=>Reaction%Constants(1,1,1,:)
    SELECT CASE (Reaction%TypeR)
      CASE ('DCONST')
        CALL DConstCompute(w,v,Constants)
      CASE ('DTEMP')
        CALL DTempCompute(w,v,Constants,Temp)
      CASE ('DTEMP2')
        CALL DTemp2Compute(w,v,Constants,Temp)
      CASE ('DTEMP3')
        CALL DTemp3Compute(w,v,Constants,Temp)
      CASE ('DTEMP4')
        CALL DTemp4Compute(w,v,Constants,Temp)
      CASE DEFAULT
        w=Zero
        STOP 'Falscher ReaktionsTyp in DissChemie'
    END SELECT
    CALL Mult(v,w)
    IF (Reaction%TypeE) THEN
      DO l=1,Reaction%NumSpeciesLeftAktiv
        is=Reaction%SpeciesLeft(l)
        il=Reaction%SpeciesLeftLoc(l)
        CALL Mult(Act(is),cAqua(is),-Reaction%Species(il)%Koeff,MolMass(is),cAqua(isW),w)
      END DO
      DO l=1,Reaction%NumSpeciesRightAktiv-1
        is=Reaction%SpeciesRight(l)
        il=Reaction%SpeciesRightLoc(l)
        CALL Mult(Act(is),cAqua(is),Reaction%Species(il)%Koeff,MolMass(is),cAqua(isW),v)
      END DO
    ELSE
      DO l=1,Reaction%NumSpeciesLeftAktiv
        is=Reaction%SpeciesLeft(l)
        CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),w)
      END DO
      DO l=1,Reaction%NumSpeciesRightAktiv-1
        is=Reaction%SpeciesRight(l)
        CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),v)
      END DO
    END IF
    CALL Add(-One,v,w)
    CALL Mult(cAqua(isW),w)
    DO l=1,Reaction%NumSpecies
      is=Reaction%Species(l)%Species
      Koeff=Reaction%Species(l)%Koeff
      CALL Add(Koeff*MolMass(is),w,f(is),iVec1,iVec2)
    END DO
    Reaction=>Reaction%Next
  END DO
  
  CALL Deallocate(w)
  CALL Deallocate(v)

END SUBROUTINE DissChemie

SUBROUTINE  DissChemieJac(cAqua,Jac,Temp,Act)

  TYPE(Vec4_T), TARGET :: cAqua(0:)
  TYPE(Vec4_T), POINTER :: Jac(:)
  TYPE(Vec4_T), TARGET :: Temp
  TYPE (Vec4_T) :: Act(0:)

  INTEGER :: i,il,ilLast,idf,istr,l,is,isW,isD,iTAbs
  INTEGER :: iReak
  REAL(RealKind) :: Fac,Koeff
  TYPE(Reaction_T), POINTER :: Reaction
  TYPE(Vec4_T) :: w,v,Work
  REAL(RealKind), POINTER :: Constants(:)

  cVecReac=>cAqua
  Temperature=>Temp
  isW=Position('aH2O')
  iTAbs=PositionGas('TE')
  CALL Allocate(w,cAqua(1))
  CALL Allocate(v,cAqua(1))
  CALL Allocate(Work,cAqua(1))
  Reaction=>DissocFirst
  DO iReak=1,nReakDissoc
    Constants=>Reaction%Constants(1,1,1,:)
    SELECT CASE (Reaction%TypeR)
      CASE ('DCONST')
        CALL DConstCompute(w,v,Constants)
      CASE ('DTEMP')
        CALL DTempCompute(w,v,Constants,Temp)
      CASE ('DTEMP2')
        CALL DTemp2Compute(w,v,Constants,Temp)
      CASE ('DTEMP3')
        CALL DTemp3Compute(w,v,Constants,Temp)
      CASE ('DTEMP4')
        CALL DTemp4Compute(w,v,Constants,Temp)
      CASE DEFAULT
        w=Zero
        STOP 'Falscher ReaktionsTyp in DissChemieJac'
    END SELECT
    CALL Mult(v,w)
    istr=0
    DO i=1,Reaction%NumSpeciesLeftAktiv
      Work=w
      isD=Reaction%SpeciesLeft(i)
      CALL Mult(Act(isD),Work)
      IF (Reaction%TypeE) THEN
        DO l=1,Reaction%NumSpeciesLeftAktiv
          is=Reaction%SpeciesLeft(l)
          il=Reaction%SpeciesLeftLoc(l)
          IF (l/=i) THEN
            CALL Mult(Act(is),cAqua(is),-Reaction%Species(il)%Koeff,MolMass(is),cAqua(isW),Work)
          ELSE
            ilLast=il
            CALL MultD(Act(is),cAqua(is),-Reaction%Species(il)%Koeff-One,  &
                      MolMass(is),cAqua(isW),Work)
            Work%c=(-Reaction%Species(il)%Koeff)*Work%c
          END IF
        END DO
      ELSE
        DO l=1,Reaction%NumSpeciesLeftAktiv
          IF (l/=i) THEN
            is=Reaction%SpeciesLeft(l)
            CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),Work)
          END IF
        END DO
      END IF
      DO l=1,Reaction%NumSpecies
        is=Reaction%Species(l)%Species
        Koeff=Reaction%Species(l)%Koeff
        istr=istr+1
        idf=Reaction%Struct(istr)
        CALL Add(Koeff*MolMass(is)/MolMass(isD),Work,Jac(idf),iVec1,iVec2)
      END DO
    END DO
!   Derivative with respect to water    
    is=Reaction%SpeciesLeft(Reaction%NumSpeciesLeftAktiv)
    IF (Reaction%TypeE) THEN
       Work%c=1.0d0/(-Reaction%Species(ilLast)%Koeff)*Work%c
    END IF   
    CALL Mult(cAqua(is),MolMass(is),cAqua(isW),Work)
    w=Work
!-----------------------------------------------------------    
    DO i=1,Reaction%NumSpeciesRightAktiv-1
      Work=v
      isD=Reaction%SpeciesRight(i)
      CALL Mult(Act(isD),Work)
      IF (Reaction%TypeE) THEN
        DO l=1,Reaction%NumSpeciesRightAktiv-1
          is=Reaction%SpeciesRight(l)
          il=Reaction%SpeciesRightLoc(l)
          IF (l/=i) THEN
            CALL Mult(Act(is),cAqua(is),Reaction%Species(il)%Koeff,MolMass(is),cAqua(isW),Work)
          ELSE
            ilLast=il
            CALL MultD(Act(is),cAqua(is),Reaction%Species(il)%Koeff-One,  &
                      MolMass(is),cAqua(isW),Work)
            Work%c=Reaction%Species(il)%Koeff*Work%c
          END IF
        END DO
      ELSE
        DO l=1,Reaction%NumSpeciesRightAktiv-1
          IF (l/=i) THEN
            is=Reaction%SpeciesRight(l)
            CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),Work)
          END IF
        END DO
      END IF
      DO l=1,Reaction%NumSpecies
        is=Reaction%Species(l)%Species
        Koeff=Reaction%Species(l)%Koeff
        istr=istr+1
        idf=Reaction%Struct(istr)
        CALL Add(-Koeff*MolMass(is)/MolMass(isD),Work,Jac(idf),iVec1,iVec2)
      END DO
    END DO
!   Derivative with respect to water    
    is=Reaction%SpeciesRight(Reaction%NumSpeciesRightAktiv-1)
    IF (Reaction%TypeE) THEN
      Work%c=1.0d0/Reaction%Species(ilLast)%Koeff*Work%c
    END IF  
    CALL Mult(cAqua(is),MolMass(is),cAqua(isW),Work)
    IF (Reaction%TypeE) THEN
      Fac=-1.0d0
      DO l=1,Reaction%NumSpecies
        Koeff=Reaction%Species(l)%Koeff
        Fac=Fac+MAX(Koeff,0.0d0)
      END DO  
    ELSE
      Fac=Reaction%NumSpeciesRightAktiv-2
    END IF  
    CALL ScaleV(Fac,Work)
    IF (Reaction%TypeE) THEN
      Fac=-1.0d0
      DO l=1,Reaction%NumSpecies
        Koeff=Reaction%Species(l)%Koeff
        Fac=Fac-MIN(Koeff,0.0d0)
      END DO  
      Fac=-Fac
    ELSE
      Fac=-(Reaction%NumSpeciesLeftAktiv-1)
    END IF  
    CALL Add(Fac,w,Work)
    DO l=1,Reaction%NumSpecies
      is=Reaction%Species(l)%Species
      Koeff=Reaction%Species(l)%Koeff
      istr=istr+1
      idf=Reaction%Struct(istr)
      CALL Add(Koeff*MolMass(is),Work,Jac(idf),iVec1,iVec2)
    END DO
    Reaction=>Reaction%Next
  END DO

  CALL Deallocate(w)
  CALL Deallocate(v)
  CALL Deallocate(Work)

END SUBROUTINE DissChemieJac


SUBROUTINE SolidChemie(cAqua,f,Temp,Act)

  TYPE (Vec4_T), TARGET :: cAqua(0:)
  TYPE (Vec4_T) :: f(0:)
  TYPE (Vec4_T), TARGET :: Temp
  TYPE (Vec4_T) :: Act(0:)

  INTEGER :: l,is,isW,iTAbs
  INTEGER :: iReak,il
  REAL(RealKind) :: Koeff
  TYPE (Reaction_T), POINTER :: Reaction
  TYPE(Vec4_T) :: w,v
  TYPE(Vec4_T) :: fL,fR
  REAL(RealKind), POINTER :: Constants(:)


  cVecReac=>cAqua
  Temperature=>Temp
  isW=Position('aH2O')
  CALL Init_StackVector(cAqua(1))
  CALL Allocate(w,cAqua(1))
  CALL Allocate(v,cAqua(1))
  CALL Allocate(fL,cAqua(1))
  CALL Allocate(fR,cAqua(1))
  Reaction=>SolidFirst
  DO iReak=1,nReakSolid
    Constants=>Reaction%Constants(1,1,1,:)
    SELECT CASE (Reaction%TypeR)
      CASE ('DCONST')
        CALL DConstCompute(w,v,Constants)
      CASE ('DTEMP')
        CALL DTempCompute(w,v,Constants,Temp)
      CASE ('DTEMP2')
        CALL DTemp2Compute(w,v,Constants,Temp)
      CASE ('DTEMP3')
        CALL DTemp3Compute(w,v,Constants,Temp)
        v=1.e-4_RealKind !OSSI
      CASE ('SPECIAL')
        CALL Calculate(w,Reaction%Infix,Reaction%Constants)
        v=One
        v=1.e-4_RealKind !OSSI
      CASE ('MESKHIDZE')  
        CALL MeskhidzeCompute(w,v,Constants,Temp)
      CASE DEFAULT
        w=Zero
        WRITE(*,*) 'Falscher ReaktionsTyp in SolidChemie',Reaction%TypeR
        STOP 'Falscher ReaktionsTyp in SolidChemie'
    END SELECT
    fR=One
    IF (Reaction%TypeE) THEN
      DO l=1,Reaction%NumSpeciesRightAktiv-1
        is=Reaction%SpeciesRight(l)
        il=Reaction%SpeciesRightLoc(l)
        CALL Mult(Act(is),cAqua(is),Reaction%Species(il)%Koeff,MolMass(is),cAqua(isW),fR)
      END DO
    ELSE
      DO l=1,Reaction%NumSpeciesRightAktiv-1
        is=Reaction%SpeciesRight(l)
        CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),fR)
      END DO
    END IF  
    fL=One
    is=Reaction%SpeciesLeft(1)
    CALL Mult(cAqua(is),MolMass(is),cAqua(isW),fL)
!   fR%cInt=w%cInt-fR%cInt
    fR%cInt=One-fR%cInt/w%cInt
    w%cInt=v%cInt* &
           (fR%cInt+fL%cInt-SQRT(fR%cInt*fR%cInt+fL%cInt*fL%cInt+1.d-80)) 
    CALL Mult(cAqua(isW),w)
    DO l=1,Reaction%NumSpecies
      is=Reaction%Species(l)%Species
      Koeff=Reaction%Species(l)%Koeff
      CALL Add(Koeff*MolMass(is),w,f(is),iVec1,iVec2)
    END DO
    Reaction=>Reaction%Next
  END DO
  
  CALL Deallocate(w)
  CALL Deallocate(v)
  CALL Deallocate(fL)
  CALL Deallocate(fR)
  CALL Close_StackVector

END SUBROUTINE SolidChemie

SUBROUTINE  SolidChemieJac(cAqua,Jac,Temp,Act)

  TYPE(Vec4_T), TARGET :: cAqua(0:)
  TYPE(Vec4_T), POINTER :: Jac(:)
  TYPE(Vec4_T), TARGET :: Temp
  TYPE (Vec4_T) :: Act(0:)

  INTEGER :: i,idf,istr,l,is,isW,isD,iTAbs
  INTEGER :: iReak,il
  REAL(RealKind) :: Fac,Koeff
  TYPE(Reaction_T), POINTER :: Reaction
  TYPE(Vec4_T) :: w,v,Work,Work1
  TYPE(Vec4_T) :: fL,fR
  REAL(RealKind), POINTER :: Constants(:)

  cVecReac=>cAqua
  Temperature=>Temp
  isW=Position('aH2O')
  iTAbs=PositionGas('TE')
  CALL Init_StackVector(cAqua(1))
  CALL Allocate(w,cAqua(1))
  CALL Allocate(v,cAqua(1))
  CALL Allocate(Work,cAqua(1))
  CALL Allocate(Work1,cAqua(1))
  CALL Allocate(fL,cAqua(1))
  CALL Allocate(fR,cAqua(1))
  Reaction=>SolidFirst
  DO iReak=1,nReakSolid
    Constants=>Reaction%Constants(1,1,1,:)
    SELECT CASE (Reaction%TypeR)
      CASE ('DCONST')
        CALL DConstCompute(w,v,Constants)
      CASE ('DTEMP')
        CALL DTempCompute(w,v,Constants,Temp)
      CASE ('DTEMP2')
        CALL DTemp2Compute(w,v,Constants,Temp)
      CASE ('DTEMP3')
        CALL DTemp3Compute(w,v,Constants,Temp)
        v=1.e-4_RealKind !OSSI
      CASE ('SPECIAL')
        CALL Calculate(w,Reaction%Infix,Reaction%Constants)
        v=One
        v=1.e-4_RealKind !OSSI
      CASE ('MESKHIDZE')  
        CALL MeskhidzeCompute(w,v,Constants,Temp)
      CASE DEFAULT
        w=Zero
        WRITE(*,*) 'Falscher ReaktionsTyp in SolidChemieJac',Reaction%TypeR
        STOP 'Falscher ReaktionsTyp in SolidChemieJac'
    END SELECT
    fR=One
    IF (Reaction%TypeE) THEN
      DO l=1,Reaction%NumSpeciesRightAktiv-1
        is=Reaction%SpeciesRight(l)
        il=Reaction%SpeciesRightLoc(l)
        cAqua(is)%c=cAqua(is)%c+EpsChemie !OSSI
        CALL Mult(Act(is),cAqua(is),Reaction%Species(il)%Koeff,MolMass(is),cAqua(isW),fR)
      END DO
    ELSE
      DO l=1,Reaction%NumSpeciesRightAktiv-1
        is=Reaction%SpeciesRight(l)
        CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),fR)
      END DO
    END IF
!   fR%cInt=w%cInt-fR%cInt
    fR%cInt=One-fR%cInt/w%cInt
    fL=One
    is=Reaction%SpeciesLeft(1)
    CALL Mult(cAqua(is),MolMass(is),cAqua(isW),fL)

    Work%cInt=v%cInt*(One-fL%cInt/SQRT(fL%cInt*fL%cInt+fR%cInt*fR%cInt+1.d-80)) 
    istr=0
    DO i=1,Reaction%NumSpeciesLeftAktiv
      isD=Reaction%SpeciesLeft(i)
      DO l=1,Reaction%NumSpecies
        is=Reaction%Species(l)%Species
        Koeff=Reaction%Species(l)%Koeff
        istr=istr+1
        idf=Reaction%Struct(istr)
        CALL Add(Koeff*MolMass(is)/MolMass(isD),Work,Jac(idf),iVec1,iVec2)
      END DO
    END DO

    Work1%cInt=v%cInt*(-(One-fR%cInt/SQRT(fL%cInt*fL%cInt+fR%cInt*fR%cInt+1.d-80)))/w%cInt
    DO i=1,Reaction%NumSpeciesRightAktiv-1
      Work=Work1
      isD=Reaction%SpeciesRight(i)
      CALL Mult(Act(isD),Work)
      IF (Reaction%TypeE) THEN
        DO l=1,Reaction%NumSpeciesRightAktiv-1
          is=Reaction%SpeciesRight(l)
          il=Reaction%SpeciesRightLoc(l)
          cAqua(is)%c=cAqua(is)%c+EpsChemie !OSSI
          IF (l/=i) THEN
            CALL Mult(Act(is),cAqua(is),Reaction%Species(il)%Koeff,MolMass(is),cAqua(isW),Work)
          ELSE
            CALL MultD(Act(is),cAqua(is),Reaction%Species(il)%Koeff-One,  &
                      MolMass(is),cAqua(isW),Work)
            Work%c=Reaction%Species(il)%Koeff*Work%c
          END IF
        END DO
      ELSE
        DO l=1,Reaction%NumSpeciesRightAktiv-1
          IF (l/=i) THEN
            is=Reaction%SpeciesRight(l)
            CALL Mult(Act(is),cAqua(is),MolMass(is),cAqua(isW),Work)
          END IF
        END DO
      END IF
      DO l=1,Reaction%NumSpecies
        is=Reaction%Species(l)%Species
        Koeff=Reaction%Species(l)%Koeff
        istr=istr+1
        idf=Reaction%Struct(istr)
        CALL Add(Koeff*MolMass(is)/MolMass(isD),Work,Jac(idf),iVec1,iVec2)
      END DO
    END DO

    Work%cInt=v%cInt &
              *((fL%cInt+fR%cInt-SQRT(fL%cInt*fL%cInt+fR%cInt*fR%cInt+1.d-80)) &
              -((1.0d0-fL%cInt/SQRT(fL%cInt*fL%cInt+fR%cInt*fR%cInt+1.d-80)))*fL%cInt &
              +2.0d0*((1.0d0-fR%cInt/SQRT(fL%cInt*fL%cInt+fR%cInt*fR%cInt+1.d-80)))*(fR%cInt-w%cInt))
    DO l=1,Reaction%NumSpecies
      is=Reaction%Species(l)%Species
      Koeff=Reaction%Species(l)%Koeff
      istr=istr+1
      idf=Reaction%Struct(istr)
    END DO
  END DO
  CALL Deallocate(fR)
  CALL Deallocate(Work)
  CALL Deallocate(Work1)
  CALL DeAllocate(w)
  CALL DeAllocate(v)
  CALL DeAllocate(fL)
  CALL Close_StackVector

END SUBROUTINE SolidChemieJac

SUBROUTINE InitGas(VecT,RhoCell,FileName)

  TYPE(Vector4Cell_T), POINTER :: VecT(:)
  TYPE(ScalarCell_T), POINTER :: RhoCell(:)
  CHARACTER(*) :: FileName

  INTEGER :: Pos,iz
  REAL(RealKind) :: c1,tAbs,p,cValue
  CHARACTER(20) :: S1,S2,End,SpeciesName,GasUnit
  CHARACTER*300 :: Line
  LOGICAL :: Back


  S1='BEGIN_GAS'
  S2='BEGIN_INIT'
  End='END_INIT'

  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,Name2=GasUnit,R1=c1)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        tAbs=VecT(ibLoc)%Vec(thPos)%c(ix0+1,iy0+1,iz0+1,1)
!       p=VecT(ibLoc)%Vec(prePos)%c(ix0+1,iy0+1,iz0+1,1)   OSSI
        IF (TRIM(GasUnit)=='pp') THEN
!         ppb --> molec/cm^3
          cValue=c1*2.4615d10*(p/101300.0d0)*(298.0d0/tAbs) ! ppb --> molec/cm^3
          cValue=cValue/Avga*1.0d06 ! molec/cm^3 --> Mol/m^3
        ELSE IF(TRIM(GasUnit)=='molec') THEN
          cValue=c1/Avga*1.0d06 ! molec/cm^3 --> Mol/m^3
          cValue=c1
        ELSE IF(TRIM(GasUnit)=='kg') THEN
          cValue=c1
        END IF
!        VecT(ibLoc)%Vec(Pos)%c=c1
!        VecT(ibLoc)%Vec(Pos)%c(:,:,:,1)=cValue*VolC/(VolC+Eps)
         VecT(ibLoc)%Vec(Pos)%c(:,:,:,1)=RhoCell(ibLoc)%c(:,:,:,1)*cValue*VolC/(VolC+Eps)
      END DO
    END IF
  END DO
  CALL CloseFile

END SUBROUTINE InitGas

SUBROUTINE InitAmbientGas(VecT,FileName)

  CHARACTER(*) :: FileName
  TYPE(Vector4Cell_T), POINTER :: VecT(:)

  INTEGER :: i,k,iCell,Pos
  INTEGER :: ix,iy,iz
  REAL(RealKind) :: c1,tAbs,p,cValue
  CHARACTER(20) :: S1,S2,End,SpeciesName,GasUnit
  CHARACTER*300 :: Line
  LOGICAL :: Back


  S1='BEGIN_GAS'
  S2='BEGIN_AMBIENT'
  End='END_AMBIENT'
  CALL OpenFile(FileName)
  DO
    CALL LineFile(Back,S1,S2,End,Name1=SpeciesName,Name2=GasUnit,R1=c1)
    IF (Back) THEN
      EXIT
    END IF
    Pos=Position(SpeciesName)
    IF (Pos>0) THEN
      DO ibLoc=1,nbLoc
        ib=LocGlob(ibLoc)
        CALL Set(Floor(ib))
        tAbs=VecT(ibLoc)%Vec(thPos)%c(ix0+1,iy0+1,iz0+1,1)
!       p=VecT(ibLoc)%Vec(prePos)%c(ix0+1,iy0+1,iz0+1,1)   OSSI
        IF (TRIM(GasUnit)=='pp') THEN
!         ppb --> molec/cm^3
          cValue=c1*2.4615d10*(p/101300.0d0)*(298.0d0/tAbs) ! ppb --> molec/cm^3
          cValue=cValue/Avga*1.0d06 ! molec/cm^3 --> Mol/m^3
        ELSE IF(TRIM(GasUnit)=='molec') THEN
          cValue=c1/Avga*1.0d06 ! molec/cm^3 --> Mol/m^3
        ELSE IF(TRIM(GasUnit)=='kg') THEN
          cValue=c1
        END IF
         VecAmb(ibLoc)%Vec(Pos)%c(:,:,:,1)=cValue*VolC/(VolC+Eps)
      END DO
    END IF
  END DO
  CALL CloseFile

END SUBROUTINE InitAmbientGas

subroutine To_upper(str)
  character(*), intent(in out) :: str
  integer :: i

  do i = 1, len(str)
    select case(str(i:i))
      case("a":"z")
        str(i:i) = achar(iachar(str(i:i))-32)
    end select
  end do 
end subroutine To_upper

subroutine To_lower(str)
  character(*), intent(in out) :: str
  integer :: i

  do i = 1, len(str)
    select case(str(i:i))
      case("A":"Z")
        str(i:i) = achar(iachar(str(i:i))+32)
    end select
  end do  
end subroutine To_Lower

END MODULE Chemie_Mod
