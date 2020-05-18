MODULE Chemsys1_Mod
  USE String_Mod
  USE hashtbl
  IMPLICIT NONE

  INTEGER, PARAMETER :: RealKind=8
  INTEGER, PARAMETER :: LenLine=400
  INTEGER, PARAMETER :: LenName=60

  TYPE Duct_T
    CHARACTER(LenName) :: Species
    CHARACTER(20) :: Type
    REAL(RealKind) :: Coeff
  END TYPE Duct_T
  TYPE Reaction_T
    CHARACTER*20 :: Type
    CHARACTER*20 :: TypeConstant
    CHARACTER(LenLine) :: Line1
    CHARACTER(LenLine) :: Line2
    TYPE(Duct_T), POINTER :: Educt(:)=>NULL()
    TYPE(Duct_T), POINTER :: Product(:)=>NULL()
    REAL(RealKind), POINTER :: Constants(:)=>NULL()
    TYPE(Reaction_T), POINTER :: Next=>NULL()
    LOGICAL :: Factor=.FALSE.
    CHARACTER*50 :: NameFactor=''
    REAL(8) :: FacCoeff1=1.0d0,FacCoeff2=1.0d0
  END TYPE Reaction_T
  TYPE ListReaction_T
    TYPE(Reaction_T), POINTER :: Start=>NULL()
    TYPE(Reaction_T), POINTER :: End=>NULL()
    INTEGER :: LenList=0
  END TYPE ListReaction_T
  TYPE Species_T
    CHARACTER(LenName) :: Species=''
    REAL(8) :: Hf=0.0d0,Gf=0.0d0,Cp=0.0d0
  END TYPE Species_T
  TYPE Element_T
    CHARACTER*5 :: Element=''
  END TYPE Element_T
  TYPE(Element_T) :: Elements(11)=(/Element_T('(') &
                                  ,Element_T(')') & 
                                  ,Element_T('exp') & 
                                  ,Element_T('+') & 
                                  ,Element_T('-') & 
                                  ,Element_T('*') & 
                                  ,Element_T('/') & 
                                  ,Element_T('**') & 
                                  ,Element_T('abs') & 
                                  ,Element_T('sqrt') & 
                                  ,Element_T('log') & 
                                  /)
    

  TYPE(Reaction_T), POINTER :: System
  LOGICAL, PRIVATE :: InitSystem=.TRUE.
  TYPE(ListReaction_T), SAVE :: ListRGas,ListRHenry,ListRAqua,ListRDiss,ListRSolid,ListRPartic,ListRMicro
  TYPE(hash_tbl_sll), PRIVATE :: ListAqua,ListGas,ListSolid,ListPartic,ListNonReac
  TYPE(Species_T), ALLOCATABLE, TARGET :: ListAqua2(:),ListGas2(:),ListSolid2(:) &
                                        ,ListPartic2(:),ListNonReac2(:)
  INTEGER :: InputUnit=10
  INTEGER, PARAMETER :: MaxEduct=10
  INTEGER, PARAMETER :: MaxProduct=10
  CHARACTER(33), PARAMETER :: SetSpecies='ABCDEFGHIJKLMNOPQRSTUVWXYZapsc[]()=+*'
  CHARACTER(14), PARAMETER :: SetConstants='ABCDEFGKINMOR/'
  CHARACTER(12), PARAMETER :: SetExponent='0123456789+-'
  INTEGER :: NumberSpeciesGas=0
  INTEGER :: NumberSpeciesAqua=0
  INTEGER :: NumberSpeciesSolid=0
  INTEGER :: NumberSpeciesSolub=0 ! counter not implemented (willi)
  INTEGER :: NumberSpeciesPartic=0
  INTEGER :: NumberSpeciesMicro=0
  INTEGER :: NumberSpeciesNonReac=0

  INTEGER :: NumberReactionsGas=0
  INTEGER :: NumberReactionsGasPhoto=0
  INTEGER :: NumberReactionsGasConst=0
  INTEGER :: NumberReactionsGasTemp=0
  INTEGER :: NumberReactionsGasTroe=0
  INTEGER :: NumberReactionsGasSpecial=0
  INTEGER :: NumberReactionsHenry=0
  INTEGER :: NumberReactionsAqua=0
  INTEGER :: NumberReactionsAquaPhoto=0
  INTEGER :: NumberReactionsAquaConst=0
  INTEGER :: NumberReactionsAquaTemp=0
  INTEGER :: NumberReactionsAquaSpecial=0
  INTEGER :: NumberReactionsDiss=0
  INTEGER :: NumberReactionsSolid=0
  INTEGER :: NumberReactionsSolidEqui=0
  INTEGER :: NumberReactionsSolidDtemp3=0
  INTEGER :: NumberReactionsSolidSpecial=0
  INTEGER :: NumberReactionsPartic=0
  INTEGER :: NumberReactionsMicro=0
  INTEGER :: UnitGas=0
  INTEGER :: UnitAqua=0
  CHARACTER*20 :: Filename='Salt'
  REAL(8), PARAMETER :: RGas=8.3145d0 
  REAL(8), PARAMETER :: TRef=298.15d0  

CONTAINS

SUBROUTINE ReadSpecies(Out)
  LOGICAL :: Out

  CHARACTER(100) :: Species
  CHARACTER(20) :: Type
  INTEGER :: Pos

  READ(InputUnit,'(a100)',END=1) Species
  DO 
    Pos=SCAN(Species,"'")
    IF (Pos>0) THEN
      Species(Pos:)=Species(Pos+1:)
    ELSE
      EXIT
    END IF
  END DO  
  IF (Species/='') THEN
    CALL InsertSpecies(Species,Type)
  END IF  
  Out=.FALSE.
  GO TO 999
1 CONTINUE
  Out=.TRUE.
999 CONTINUE
END SUBROUTINE ReadSpecies

SUBROUTINE ReadReaction(Out)

  LOGICAL :: Out

  INTEGER :: iLine,PosColon,Pos,is
  CHARACTER(LenLine) :: LocLine
  CHARACTER(LenLine) :: Line(1:5)
  CHARACTER*20 :: Type
  CHARACTER*40 :: TypeR

  iLine=0
  S1:DO 
    Line=''
    READ(InputUnit,'(a400)',IOSTAT=is) LocLine
    IF (ABS(is)>0) THEN
      EXIT
    END IF  
    IF (INDEX(LocLine,'CLASS')>0) THEN
      iLine=1
      Line(iLine)=LocLine
      DO
        READ(InputUnit,'(a400)',IOSTAT=is) LocLine
        WRITE(*,*) 'LocLine 1',LocLine(1:30)
        IF (ABS(is)>0) THEN
          EXIT S1
        END IF  
        IF (INDEX(LocLine,'CLASS')>0) THEN
          BACKSPACE(InputUnit)
          EXIT S1
        END IF  
        IF (ADJUSTL(LocLine(1:1))/='#'.AND.ADJUSTL(LocLine(1:7))/='COMMENT'.AND.LEN(TRIM(LocLine))>0) THEN
          iLine=iLine+1
          Line(iLine)=LocLine
        END IF
      END DO  
    END IF 
  END DO S1
  IF (iLine>=3) THEN
    Pos=SCAN(Line(1),'#')
    IF (Pos>0) THEN
      Line(1)=Line(1)(:Pos-1)
    END IF  
    PosColon=Index(Line(1),':')
    Type=ADJUSTL(Line(1)(PosColon+1:))
    SELECT CASE (Type)
      CASE ('GAS')
        NumberReactionsGas=NumberReactionsGas+1
        CALL InsertReaction(ListRGas,Line,TypeR) 
        SELECT CASE (TypeR)
          CASE ('PHOTO','PHOTMCM','PHOTMCM1')
            NumberReactionsGasPhoto=NumberReactionsGasPhoto+1
          CASE ('CONST')
            NumberReactionsGasConst=NumberReactionsGasConst+1
          CASE ('TEMP','TEMP1')
            NumberReactionsGasTemp=NumberReactionsGasTemp+1
          CASE ('TROE','TROEF','TROEQ')
            NumberReactionsGasTroe=NumberReactionsGasTroe+1
          CASE ('SPECIAL')
            NumberReactionsGasSpecial=NumberReactionsGasSpecial+1
        END SELECT
      CASE ('HENRY')
        NumberReactionsHenry=NumberReactionsHenry+1
        CALL InsertReaction(ListRHenry,Line,TypeR) 
      CASE ('AQUA')
        NumberReactionsAqua=NumberReactionsAqua+1
        CALL InsertReaction(ListRAqua,Line,TypeR) 
        SELECT CASE (TypeR)
          CASE ('PHOTO','PHOTMCM','PHOTMCM1')
            NumberReactionsAquaPhoto=NumberReactionsAquaPhoto+1
          CASE ('CONST')
            NumberReactionsAquaConst=NumberReactionsAquaConst+1
          CASE ('TEMP','TEMP1','TEMP3')
            NumberReactionsAquaTemp=NumberReactionsAquaTemp+1
          CASE ('SPECIAL')
            NumberReactionsAquaSpecial=NumberReactionsAquaSpecial+1
        END SELECT
      CASE ('DISS')
        NumberReactionsDiss=NumberReactionsDiss+1
        CALL InsertReaction(ListRDiss,Line,TypeR) 
      CASE ('SOLID')
        NumberReactionsSolid=NumberReactionsSolid+1
        CALL InsertReaction(ListRSolid,Line,TypeR) 
        SELECT CASE (TypeR)
          CASE ('EQUI')
            NumberReactionsSolidEqui=NumberReactionsSolidEqui+1
          CASE ('DTEMP3')
            NumberReactionsSolidDtemp3=NumberReactionsSolidDtemp3+1
          CASE ('SPECIAL')
            NumberReactionsSolidSpecial=NumberReactionsSolidSpecial+1
        END SELECT
      CASE ('PARTI')
        NumberReactionsPartic=NumberReactionsPartic+1
        CALL InsertReaction(ListRPartic,Line,TypeR) 
      CASE ('MICROPHYS')
        NumberReactionsMicro=NumberReactionsMicro+1
        CALL InsertReaction(ListRMicro,Line,TypeR) 
    END SELECT
    Out=.FALSE.
  ELSE  
    Out=.TRUE.
  END IF  

END SUBROUTINE ReadReaction

SUBROUTINE PrintSpecies(ListName,Unit)
  TYPE(Species_T) :: ListName(:)
  INTEGER :: Unit

  INTEGER :: i

  DO i=1,SIZE(ListName)
    WRITE(Unit,*) "'"//TRIM(ListName(i)%Species)//"'"
  END DO
END SUBROUTINE PrintSpecies

SUBROUTINE PrintHeadSpecies(Unit)

  INTEGER :: Unit

  CHARACTER*8 :: Date
  CHARACTER*10 :: Time
  INTEGER(8) :: Value(8)

  CALL DATE_AND_TIME(Date,Time,VALUES=Value) 

  WRITE(Unit,*) ' ==========================================================='
  WRITE(Unit,*) ' ========                 MODMEP                    ========'
  WRITE(Unit,*) ' ========     Output -  Chemical Reaction Data      ========'
  WRITE(Unit,*) ' ==========================================================='
  WRITE(Unit,*) ''
  WRITE(Unit,*) ' Created:             ',Date(7:8),'.',Date(5:6),'.',Date(1:4)
  WRITE(Unit,*) ' Chemical Mechanism:  ',TRIM(ADJUSTL(FileName)),'.chem'
  WRITE(Unit,*) ''
  WRITE(Unit,*) ' =================     Units         ======================='
  WRITE(Unit,*) ''
  IF (UnitGas==0) THEN
    WRITE(Unit,*) ' Gas Phase Units:     molec/cm3'
  ELSE  
    WRITE(Unit,*) ' Gas Phase Units:     mol/m3'
  END IF  
  IF (UnitAqua==0) THEN
    WRITE(Unit,*) ' Aqueous Phase Units: mol/l'
  END IF  
  WRITE(Unit,*) ''
  WRITE(Unit,*) ' =================    Numbers        ======================='
  WRITE(Unit,*) ''
  WRITE(Unit,*) NumberSpeciesGas &
               +NumberSpeciesAqua &
               +NumberSpeciesSolub &
               +NumberSpeciesPartic &
               +NumberSpeciesNonReac &
               +NumberSpeciesSolid,  '     Number of Species' 
  WRITE(Unit,*) NumberSpeciesGas,    '     No. of gaseous species'           
  WRITE(Unit,*) NumberSpeciesAqua,   '     No. of aqueous species'           
  WRITE(Unit,*) NumberSpeciesSolub,  '     No. of soluble species'           
  WRITE(Unit,*) NumberSpeciesPartic, '     No. of particular species'
  WRITE(Unit,*) NumberSpeciesSolid,  '     No. of solid   species'           
  WRITE(Unit,*) NumberSpeciesNonReac,'     Number of Non-reactive Species '
  WRITE(Unit,*) ''
  WRITE(Unit,*) ' =================   Species Names   ======================='
  WRITE(Unit,*) ''
END SUBROUTINE PrintHeadSpecies

SUBROUTINE PrintHeadReactions(Unit)

  INTEGER :: Unit

  INTEGER :: NumberReactions
  CHARACTER*8 :: Date
  CHARACTER*10 :: Time
  INTEGER(8) :: Value(8)

  NumberReactions=NumberReactionsGas &
                 +NumberReactionsHenry &
                 +NumberReactionsAqua &
                 +NumberReactionsDiss &
                 +NumberReactionsSolid & 
                 +NumberReactionsPartic &
                 +NumberReactionsMicro 

  WRITE(Unit,*) ''
  WRITE(Unit,*) ' ================   Description of Reactions   =============='
  WRITE(Unit,*) ''
  WRITE(Unit,*) NumberReactions,            '        NREAK   : Number of Reactions' 
  WRITE(Unit,*) NumberReactionsGas,         '        NGAS   : Gas phase reactions'  
  WRITE(Unit,*) NumberReactionsGasPhoto,    '           Gaseous PHOTO - type reactions'
  WRITE(Unit,*) NumberReactionsGasConst,    '           Gaseous CONST - type reactions'
  WRITE(Unit,*) NumberReactionsGasTemp,     '           Gaseous TEMP - type reactions'
  WRITE(Unit,*) NumberReactionsGasTroe,     '           Gaseous TROE - type reactions'
  WRITE(Unit,*) NumberReactionsGasSpecial,  '           Gaseous SPECIAL - type reactions'
  WRITE(Unit,*) NumberReactionsHenry,       '        NHENRY : Henry Equilib. reactions'
  WRITE(Unit,*) NumberReactionsDiss,        '        NDISS  : Dissociation reactions'
  WRITE(Unit,*) NumberReactionsAqua,        '        NAQUA  : Aquatic Equilib. reactions'
  WRITE(Unit,*) NumberReactionsAquaPhoto,   '           Aqueous PHOTO - type reactions'
  WRITE(Unit,*) NumberReactionsAquaConst,   '           Aqueous CONST - type reactions'
  WRITE(Unit,*) NumberReactionsAquaTemp,    '           Aqueous TEMP - type reactions'
  WRITE(Unit,*) NumberReactionsAquaSpecial, '           Aqueous SPECIAL - type reactions'
  WRITE(Unit,*) NumberReactionsPartic,      '        NPARTI  : Particulare reactions   '
  WRITE(Unit,*) NumberReactionsSolid,       '        NSOLID  : Solid Equilib. reactions'
  WRITE(Unit,*) NumberReactionsSolidDTemp3, '           Solid DTEMP3 - type reactions'
  WRITE(Unit,*) NumberReactionsSolidEqui,   '           Solid EQUI - type reactions'
  WRITE(Unit,*) NumberReactionsSolidSpecial,'           Solid SPECIAL - type reactions'
  WRITE(Unit,*) NumberReactionsMicro,       '        NMICRO  : Microphysical reactions'
  WRITE(Unit,*)
  WRITE(Unit,*) ' ======================  Reactions   ========================'
  WRITE(Unit,*) ''
END SUBROUTINE PrintHeadReactions

SUBROUTINE PrintReaction(List,Unit)

  TYPE(ListReaction_T) :: List
  INTEGER :: Unit

  INTEGER :: i
  TYPE(Reaction_T), POINTER :: Current
  INTEGER, SAVE :: NumberReaction=0
  INTEGER :: NumActiveEduct,NumActiveProduct
  TYPE(Duct_T) :: ActiveEduct(30)
  TYPE(Duct_T) :: ActiveProduct(30)

  Current=>List%Start
  DO
    IF (ASSOCIATED(Current)) THEN
      NumActiveEduct=0
      DO i=1,SIZE(Current%Educt)
        SELECT CASE(Current%Educt(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic')
            NumActiveEduct=NumActiveEduct+1
            ActiveEduct(NumActiveEduct)=Current%Educt(i)
        END SELECT
      END DO  
      NumActiveProduct=0
      DO i=1,SIZE(Current%Product)
        SELECT CASE(Current%Product(i)%Type)
          CASE ('Gas','Aqua','Solid','Partic')
            NumActiveProduct=NumActiveProduct+1
            ActiveProduct(NumActiveProduct)=Current%Product(i)
        END SELECT
      END DO  
      NumberReaction=NumberReaction+1
      WRITE(Unit,'(a12,i5,a23)')   '#-----------',NumberReaction,'. Reaction ----------- '
      WRITE(Unit,*) TRIM(Current%Type)//'   '//TRIM(Current%TypeConstant)
      WRITE(Unit,555) SIZE(Current%Educt) , SIZE(Current%Product) , NumActiveEduct , NumActiveProduct 
      WRITE(Unit,555) ( PositionSpecies(Current%Educt(i)%Species) , i=1,SIZE(Current%Educt)   ) &
                   ,( PositionSpecies(Current%Product(i)%Species) , i=1,SIZE(Current%Product) ) &
                   , NumActiveEduct + NumActiveProduct
      WRITE(Unit,666)  ( PositionSpecies(ActiveEduct(i)%Species)   , -ActiveEduct(i)%Coeff  , i=1,NumActiveEduct   ) &
                     , ( PositionSpecies(ActiveProduct(i)%Species) , ActiveProduct(i)%Coeff , i=1,NumActiveProduct ) 
      IF (Current%TypeConstant=='SPECIAL') THEN
        WRITE(Unit,*) Current%Line2(:LEN(TRIM(Current%Line2))-1)
      END IF
      WRITE(Unit,*) SIZE(Current%Constants),Current%Constants
      IF (Current%Factor) THEN
        WRITE(Unit,*) 'FACTOR ',TRIM(Current%NameFactor),Current%FacCoeff1,Current%FacCoeff2
      END IF
      Current=>Current%Next
    ELSE
      EXIT
    END IF
  END DO
  
  555 FORMAT( *(2X,I0) )
  666 FORMAT( *(2X,I0,1X,Es10.2) )
  
END SUBROUTINE PrintReaction

SUBROUTINE InsertReaction(List,Line,TypeR)

  TYPE(ListReaction_T) :: List
  CHARACTER(*) :: Line(1:4)
  CHARACTER(*)  :: TypeR

  INTEGER :: PosColon,PosEqual,PosPlus,PosSpecies,PosDol,Pos
  CHARACTER(LenLine) :: Left,Right
  TYPE(Reaction_T), POINTER :: Reaction

  IF (ASSOCIATED(List%Start)) THEN
    ALLOCATE(List%End%Next)
    List%End=>List%End%Next
  ELSE
    ALLOCATE(List%End)
    List%Start=>List%End
  END IF
  List%LenList=List%LenList+1
  Reaction=>List%End
  PosColon=INDEX(Line(1),':')
  Reaction%Type=ADJUSTL(Line(1)(PosColon+1:))
  Reaction%Line1=Line(2)
  PosEqual=INDEX(Reaction%Line1,' = ')
  IF (PosEqual==0) THEN
    WRITE(*,*) 'Trennzeichen Falsch'
    STOP
  END IF  
  Left=Reaction%Line1(1:PosEqual-1)
  CALL ExtractSpecies(Left,Reaction%Educt)
  Right=Reaction%Line1(PosEqual+3:)
  CALL ExtractSpecies(Right,Reaction%Product)
  CALL ExtractConstants(Line(3),Reaction%Constants,Reaction%TypeConstant)
  Reaction%Line2=Line(3)
  TypeR=Reaction%TypeConstant
  IF (INDEX(Line(4),'FACTOR')>0) THEN
    Reaction%Factor=.TRUE.
    PosDol=INDEX(Line(4),'$')
    IF (PosDol>0) THEN
      Pos=INDEX(Line(4)(PosDol:),' ')+PosDol-1 
      Reaction%NameFactor=Line(4)(PosDol:Pos)
      IF (LEN(TRIM(Line(4)(Pos+1:)))>0) THEN
        READ(Line(4)(Pos+1:),*) Reaction%FacCoeff1,Reaction%FacCoeff2
      END IF
    END IF  
  END IF  
END SUBROUTINE InsertReaction

SUBROUTINE ReadUnits
!
INTEGER :: Pos
CHARACTER(LenLine) :: LocLine
INTEGER :: ios
CHARACTER(400) :: iom
LOGICAL :: comments
!
REWIND(InputUnit)
DO 
  READ(InputUnit,'(A200)',IOSTAT=ios,IOMSG=iom) LocLine
  IF ( ios /= 0 ) THEN
    WRITE(*,*) ' Error Reading Units' 
    WRITE(*,*) ' Error Message  ::  ',TRIM(iom)
  ELSE
    comments = (ADJUSTL(LocLine(1:1)) == '#'        .OR. &
    &           ADJUSTL(LocLine(1:7)) == 'COMMENT'  .OR. &
    &           LEN_TRIM(LocLine)     == 0 )
    IF ( comments ) CYCLE

    Pos = INDEX(LocLine,'GAS')
    IF ( Pos > 0 ) READ(LocLine(Pos+3:),*) UnitGas 
    Pos = INDEX(LocLine,'AQUA')
    IF ( Pos > 0 ) THEN
      READ(LocLine(Pos+4:),*) UnitAQUA 
      EXIT
    END IF
  END IF
END DO

END SUBROUTINE ReadUnits

SUBROUTINE ExtractConstants(String,Constants,Type)

  CHARACTER(*) :: String
  REAL(RealKind), POINTER :: Constants(:)
  CHARACTER(*) :: Type

  INTEGER :: NumColon,PosColon,PosName,PosComment
  INTEGER :: i,PosNum1,PosNum2,PosNum3,NumNum,PosElem
  CHARACTER*4 :: NameNumNum
  CHARACTER*10 :: DummyString
  CHARACTER(LEN(String)) :: LocString
  CHARACTER(LEN(String)) :: LocString1
  CHARACTER(LEN(String)) :: NameConstant
  REAL(RealKind) :: Dummy
  INTEGER :: is

  LocString=String
  String=''
  PosColon=Index(LocString,':')
  Type=LocString(:PosColon-1)
  LocString=ADJUSTL(LocString(PosColon+1:))
  PosComment=INDEX(LocString,'#')
  IF (PosComment>0) THEN
    LocString=LocString(:PosComment-1)
  END IF  
  NumColon=0
  LocString1=LocString
  IF (Type/='SPECIAL') THEN
    DO  
      PosColon=Index(LocString1,':')
      IF (PosColon>0) THEN
        LocString1=ADJUSTL(LocString1(PosColon+1:))
        READ(LocString1,*,IOSTAT=is) Dummy,DummyString
        NumColon=NumColon+1
        PosName=INDEX(LocString1,TRIM(DummyString))
        IF (PosName>0) THEN
          LocString1=LocString1(PosName:)
        END IF  
      ELSE
        EXIT
      END IF
    END DO
    ALLOCATE(Constants(NumColon))
    NumColon=0
    DO  
      PosColon=Index(LocString,':')
      IF (PosColon>0) THEN
        LocString=ADJUSTL(LocString(PosColon+1:))
        NumColon=NumColon+1
        READ(LocString,*,IOSTAT=is) Constants(NumColon),DummyString
        PosName=INDEX(LocString,TRIM(DummyString))
        IF (PosName>0) THEN
          LocString=LocString(PosName:)
        END IF  
      ELSE
        EXIT
      END IF
    END DO
  ELSE
    NumNum=0
    DO
      PosNum1=SCAN(LocString1,'1234567890.')
      PosNum2=LEN(LocString1(MAX(PosNum1,1):))
      DO i=1,SIZE(Elements)
        PosElem=INDEX(LocString1(MAX(PosNum1,1):),TRIM(Elements(i)%Element))
        IF (PosElem>0) THEN
          PosNum2=MIN(PosNum2,PosElem-1)
        END IF
      END DO
      PosNum2=PosNum2+PosNum1-1
      IF (PosNum2>0) THEN
        IF ((LocString1(PosNum2:PosNum2)=='e'.OR.LocString1(PosNum2:PosNum2)=='E').AND. &
             (LocString1(PosNum2+1:PosNum2+1)=='-'.OR. &
              LocString1(PosNum2+1:PosNum2+1)=='+')) THEN
          PosNum3=PosNum2+2
          PosNum2=LEN(LocString(MAX(PosNum1,1):))
          DO i=1,SIZE(Elements)
            PosElem=INDEX(LocString1(MAX(PosNum3,1):),TRIM(Elements(i)%Element))
            IF (PosElem>0) THEN
              PosNum2=MIN(PosNum2,PosElem-1)
            END IF
          END DO
          PosNum2=PosNum2+PosNum3-1
        ELSE
          PosNum2=PosNum2-PosNum1+1
        END IF          
      ELSE
        PosNum2=PosNum2-PosNum1+1
      END IF          
      PosNum2=PosNum2+PosNum1-1
      IF (PosNum1>0) THEN
        LocString1=LocString1(PosNum2+1:)
        NumNum=NumNum+1
      ELSE
        EXIT
      END IF
    END DO
    ALLOCATE(Constants(NumNum))
    NumNum=0
    DO
      PosNum1=SCAN(LocString,'1234567890.')
      PosNum2=LEN(LocString(MAX(PosNum1,1):))
      DO i=1,SIZE(Elements)
        PosElem=INDEX(LocString(MAX(PosNum1,1):),TRIM(Elements(i)%Element))
        IF (PosElem>0) THEN
          PosNum2=MIN(PosNum2,PosElem-1)
        END IF
      END DO
      PosNum2=PosNum2+PosNum1-1
      IF (PosNum2>0) THEN
        IF ((LocString(PosNum2:PosNum2)=='e'.OR.LocString(PosNum2:PosNum2)=='E').AND. &
             (LocString(PosNum2+1:PosNum2+1)=='-'.OR. &
              LocString(PosNum2+1:PosNum2+1)=='+')) THEN
          PosNum3=PosNum2+2
          PosNum2=LEN(LocString(MAX(PosNum1,1):))
          DO i=1,SIZE(Elements)
            PosElem=INDEX(LocString(MAX(PosNum3,1):),TRIM(Elements(i)%Element))
            IF (PosElem>0) THEN
              PosNum2=MIN(PosNum2,PosElem-1)
            END IF
          END DO
          PosNum2=PosNum2+PosNum3-1
        ELSE
          PosNum2=PosNum2-PosNum1+1
        END IF          
      ELSE
        PosNum2=PosNum2-PosNum1+1
      END IF          
      PosNum2=PosNum2+PosNum1-1
      IF (PosNum1>0) THEN
        NameConstant=LocString(PosNum1:PosNum2)
        NumNum=NumNum+1
        READ(NameConstant,*) Constants(NumNum)
        !WRITE(NameNumNum,'(I2)') NumNum
        String=TRIM(String)//LocString(:PosNum1-1)//'$'//TRIM(ADJUSTL(NameNumNum))
        LocString=LocString(PosNum2+1:)
      ELSE
        String=TRIM(String)//TRIM(LocString)
        EXIT
      END IF
    END DO
  END IF
END SUBROUTINE ExtractConstants

SUBROUTINE ExtractSpecies(String,Duct)

  CHARACTER(*) :: String
  TYPE(Duct_T), POINTER :: Duct(:)

  INTEGER :: PosMinus,PosPlus,NumSpec,PosSpecies
  REAL(RealKind) :: NumberSpecies,PreFac
  CHARACTER(LenLine) :: Species
  CHARACTER(LEN(String)) :: LocString

  LocString=String
  NumSpec=1
  DO
    IF (INDEX(LocString,' + ')>0) THEN
      LocString=LocString(INDEX(LocString,' + ')+3:)
      NumSpec=NumSpec+1
    ELSE IF (INDEX(LocString,' - ')>0) THEN
      LocString=LocString(INDEX(LocString,' - ')+3:)
      NumSpec=NumSpec+1
    ELSE
      EXIT
    END IF
  END DO
  ALLOCATE(Duct(NumSpec))
  LocString=String
  NumSpec=0
  DO
    PosPlus=INDEX(LocString,' + ')
    PosMinus=INDEX(LocString,' - ')
    IF (PosMinus>0.AND.PosMinus<PosPlus) THEN
      PreFac=-1.0d0
    ELSE
      PreFac=1.0d0
    END IF  
    IF (PosPlus>0.AND.(PosPlus<PosMinus.OR.PosMinus==0)) THEN
      Species=ADJUSTL(LocString(:PosPlus-1))
      LocString=LocString(PosPlus+3:)
    ELSE IF (PosMinus>0.AND.(PosMinus<PosPlus.OR.PosPlus==0)) THEN
      Species=ADJUSTL(LocString(:PosMinus-1))
      LocString=LocString(PosMinus+3:)
    ELSE
      Species=ADJUSTL(LocString)
    END IF
    PosSpecies=SCAN(Species,SetSpecies)
    NumSpec=NumSpec+1
    IF (PosSpecies==1) THEN
      Duct(NumSpec)%Coeff=PreFac
      Duct(NumSpec)%Species=Species
    ELSE
      READ(Species(1:PosSpecies-1),*) Duct(NumSpec)%Coeff
      Duct(NumSpec)%Coeff=PreFac*Duct(NumSpec)%Coeff
      Duct(NumSpec)%Species=Species(PosSpecies:)
    END IF
    CALL InsertSpecies(Duct(NumSpec)%Species,Duct(NumSpec)%Type)
    IF (PosPlus==0.AND.PosMinus==0) THEN
      EXIT
    END IF
  END DO
END SUBROUTINE ExtractSpecies

SUBROUTINE InsertSpecies(Species,Type)

  CHARACTER(*) :: Species
  CHARACTER(*) :: Type

  IF (Species(1:1)=='p') THEN
    CALL InsertHash(ListPartic,TRIM(ADJUSTL(Species)),NumberSpeciesPartic)
    Type='Partic'
  ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
    CALL InsertHash(ListAqua,TRIM(ADJUSTL(Species)),NumberSpeciesAqua)
    Type='Aqua'
  ELSE IF (Species(1:1)=='s') THEN
    CALL InsertHash(ListSolid,TRIM(ADJUSTL(Species)),NumberSpeciesSolid)
    Type='Solid'
  ELSE IF (Species(1:1)=='['.AND.Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
    CALL InsertHash(ListNonReac,TRIM(ADJUSTL(Species)),NumberSpeciesNonReac)
    Type='Inert'
  ELSE IF (Species(1:1)=='(') THEN
  ELSE
    CALL InsertHash(ListGas,TRIM(ADJUSTL(Species)),NumberSpeciesGas)
    Type='Gas'
  END IF
END SUBROUTINE InsertSpecies

FUNCTION PositionListSpecies(Species)

  TYPE(Species_T), POINTER :: PositionListSpecies
  CHARACTER(*) :: Species

  INTEGER :: PositionSpecies

  PositionSpecies=0
  PositionListSpecies=>NULL()
  IF (Species(1:1)=='p') THEN
    PositionSpecies=GetHash(ListPartic,TRIM(ADJUSTL(Species)))
    IF (PositionSpecies>0) PositionListSpecies=>ListPartic2(PositionSpecies)
  ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
    PositionSpecies=GetHash(ListAqua,TRIM(ADJUSTL(Species)))
    IF (PositionSpecies>0) PositionListSpecies=>ListAqua2(PositionSpecies)
  ELSE IF (Species(1:1)=='s') THEN
    PositionSpecies=GetHash(ListSolid,TRIM(ADJUSTL(Species)))
    IF (PositionSpecies>0) PositionListSpecies=>ListSolid2(PositionSpecies)
  ELSE IF (Species(1:1)=='['.AND.Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
    PositionSpecies=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))
    IF (PositionSpecies>0) PositionListSpecies=>ListNonReac2(PositionSpecies)
  ELSE
    PositionSpecies=GetHash(ListGas,TRIM(ADJUSTL(Species)))
    IF (PositionSpecies>0) PositionListSpecies=>ListGas2(PositionSpecies)
  END IF
END FUNCTION PositionListSpecies

FUNCTION PositionSpecies(Species)

  INTEGER :: PositionSpecies

  CHARACTER(*) :: Species

  PositionSpecies=0
  IF (Species(1:1)=='p') THEN
    PositionSpecies=GetHash(ListPartic,TRIM(ADJUSTL(Species)))
  ELSE IF (Species(1:1)=='a'.OR.SCAN(Species,'pm')>0) THEN
    PositionSpecies=GetHash(ListAqua,TRIM(ADJUSTL(Species)))
  ELSE IF (Species(1:1)=='s') THEN
    PositionSpecies=GetHash(ListSolid,TRIM(ADJUSTL(Species)))
  ELSE IF (Species(1:1)=='['.AND.Species(LEN(TRIM(Species)):LEN(TRIM(Species)))==']') THEN
    PositionSpecies=GetHash(ListNonReac,TRIM(ADJUSTL(Species)))
  ELSE
    PositionSpecies=GetHash(ListGas,TRIM(ADJUSTL(Species)))
  END IF
END FUNCTION PositionSpecies

SUBROUTINE OpenFile(FileName,Type)

  CHARACTER(*) :: Filename
  CHARACTER(*) :: Type

  LOGICAL :: ExistFile

  INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
  IF (ExistFile) THEN
    OPEN(UNIT=InputUnit,FILE=TRIM(Filename)//'.'//TRIM(Type),STATUS='UNKNOWN')
  END IF  

END SUBROUTINE OpenFile

SUBROUTINE CloseFile(FileName,Type)
  CHARACTER(*) :: Filename
  CHARACTER(*) :: Type
  LOGICAL :: ExistFile

  INQUIRE(FILE=TRIM(Filename)//'.'//TRIM(Type),EXIST=ExistFile)
  IF (ExistFile) THEN
    CLOSE(UNIT=InputUnit)
  END IF  
END SUBROUTINE CloseFile

SUBROUTINE ReadThermoData(FileName)

  CHARACTER(*) :: Filename

  CHARACTER(LenLine) :: LocLine
  CHARACTER(LenName) :: Name

  INTEGER :: is,iLine,i
  REAL(8) :: Hf,Gf,Cp
  TYPE(Species_T), POINTER :: Species
  TYPE(Reaction_T), POINTER :: Current

  CALL OpenFile(FileName,'dat')
  iLine=0
  DO 
    READ(InputUnit,'(a400)',IOSTAT=is) LocLine
    IF (ABS(is)>0) THEN
      EXIT
    END IF  
    IF (ADJUSTL(LocLine(1:1))/='#'.AND.LEN(TRIM(LocLine))>0) THEN
      READ(LocLine,*) Name
      IF (PositionSpecies(Name)>0) THEN
        iLine=iLine+1
      END IF  
    END IF
  END DO  
  REWIND(InputUnit)
  DO 
    READ(InputUnit,'(a400)',IOSTAT=is) LocLine
    IF (ABS(is)>0) THEN
      EXIT
    END IF  
    IF (ADJUSTL(LocLine(1:1))/='#'.AND.LEN(TRIM(LocLine))>0) THEN
      READ(LocLine,*) Name,Hf,Gf,Cp
      Species=>PositionListSpecies(Name)
      IF (ASSOCIATED(Species)) THEN
        Species%Hf=Hf
        Species%Gf=Gf
        Species%Cp=Cp
      END IF  
    END IF
  END DO  
  CLOSE(InputUnit)

  Current=>ListRSolid%Start
  DO
    IF (ASSOCIATED(Current)) THEN
      IF (Current%TypeConstant=='DTEMP3') THEN
        Hf=0.0d0
        Gf=0.0d0
        Cp=0.0d0
        DO i=1,SIZE(Current%Educt)
          Species=>PositionListSpecies(Current%Educt(i)%Species)
          Hf=Hf-Current%Educt(i)%Coeff*Species%Hf
          Gf=Gf-Current%Educt(i)%Coeff*Species%Gf
          Cp=Cp-Current%Educt(i)%Coeff*Species%Cp
        END DO
        DO i=1,SIZE(Current%Product)
          Species=>PositionListSpecies(Current%Product(i)%Species)
          Hf=Hf+Current%Product(i)%Coeff*Species%Hf
          Gf=Gf+Current%Product(i)%Coeff*Species%Gf
          Cp=Cp+Current%Product(i)%Coeff*Species%Cp
        END DO
        Hf=Hf*1000.0d0
        Gf=Gf*1000.0d0
        IF (ASSOCIATED(Current%Constants)) THEN
          DEALLOCATE(Current%Constants)
        END IF  
        ALLOCATE(Current%Constants(3))
        Current%Constants(1)=EXP(-Gf/(RGas*TRef))
        Current%Constants(2)=Cp/RGas
        Current%Constants(3)=-Hf/RGas+Cp*TRef/RGas
      END IF  
      Current=>Current%Next
    ELSE
      EXIT
    END IF
  END DO  
END SUBROUTINE ReadThermoData

SUBROUTINE ReadSystem(FileName)

  CHARACTER(*) :: Filename

  LOGICAL :: Out

  CALL InitHashTable(ListAqua,100)
  CALL InitHashTable(ListGas,100)
  CALL InitHashTable(ListSolid,100)
  CALL InitHashTable(ListPartic,100)
  CALL InitHashTable(ListNonReac,100)
  CALL OpenFile(FileName,'sys')
  CALL ReadUnits
  DO
    CALL ReadReaction(Out)
    IF (Out) THEN
      EXIT
    END IF
  END DO
  CALL CloseFile(FileName,'sys')
  CALL OpenFile(FileName,'spc')
  DO
    CALL ReadSpecies(Out)
    IF (Out) THEN
      EXIT
    END IF
  END DO
  CALL CloseFile(FileName,'spc')
  ALLOCATE(ListGas2(NumberSpeciesGas))
  CALL HashTableToList(ListGas,ListGas2)
! CALL SortList(ListGas2)
  CALL ListToHashTable(ListGas2,ListGas)
  ALLOCATE(ListAqua2(NumberSpeciesAqua))
  CALL HashTableToList(ListAqua,ListAqua2)
  CALL SortList(ListAqua2)
  CALL ListToHashTable(ListAqua2,ListAqua)
  ALLOCATE(ListSolid2(NumberSpeciesSolid))
  CALL HashTableToList(ListSolid,ListSolid2)
  CALL SortList(ListSolid2)
  CALL ListToHashTable(ListSolid2,ListSolid)
  ALLOCATE(ListPartic2(NumberSpeciesPartic))
  CALL HashTableToList(ListPartic,ListPartic2)
  CALL SortList(ListPartic2)
  CALL ListToHashTable(ListPartic2,ListPartic)
  ALLOCATE(ListNonReac2(NumberSpeciesNonReac))
  CALL HashTableToList(ListNonReac,ListNonReac2)
  CALL SortList(ListNonReac2)
  CALL ListToHashTable(ListNonReac2,ListNonReac)
END SUBROUTINE ReadSystem

SUBROUTINE SortList(List)
  TYPE(Species_T) :: List(:)

  TYPE(Species_T) :: Temp
  INTEGER :: i,j

  DO i=1,SIZE(List)
    DO j=1,SIZE(List)-i
      IF (List(j+1)%Species<List(j)%Species) THEN
        Temp=List(j+1)
        List(j+1)=List(j)
        List(j)=Temp
      END IF  
    END DO
  END DO  
  DO i=1,SIZE(List)
    IF (List(i)%Species=='OHm') THEN
      Temp=List(i)
      IF (i==SIZE(List)) EXIT
      DO j=i,SIZE(List)-1
        List(j)=List(j+1)
      END DO  
      List(SIZE(List))=Temp
      EXIT
    END IF  
  END DO  
  DO i=1,SIZE(List)
    IF (List(i)%Species=='Hp') THEN
      Temp=List(i)
      IF (i==SIZE(List)) EXIT
      DO j=i,SIZE(List)-1
        List(j)=List(j+1)
      END DO  
      List(SIZE(List))=Temp
      EXIT
    END IF  
  END DO  

END SUBROUTINE SortList

SUBROUTINE ListToHashTable(List,Table)
  TYPE(Species_T) :: List(:)
  TYPE(hash_tbl_sll) :: Table

  INTEGER :: i
  DO i=1,SIZE(List)
    CALL Table%put(TRIM(ADJUSTL(List(i)%Species)),i)
  END DO  
END SUBROUTINE ListToHashTable

SUBROUTINE HashTableToList(Table,List)
  TYPE(hash_tbl_sll) :: Table
  TYPE(Species_T) :: List(:)

  INTEGER :: i,j
  TYPE(sllist), POINTER :: child => NULL()

  DO i=LBOUND(table%vec,dim=1), UBOUND(table%vec,dim=1)
    IF (ALLOCATED(table%vec(i)%key)) THEN 
      DO j=1,SIZE(table%vec(i)%key)
        List(table%vec(i)%Val)%Species(j:j)=table%vec(i)%key(j)
      END DO
    END IF
    Child=>table%vec(i)%Child
    DO 
      IF (ASSOCIATED(Child)) THEN
        DO j=1,SIZE(Child%key)
          List(Child%Val)%Species(j:j)=Child%key(j)
        END DO
        Child=>Child%Child
      ELSE  
        EXIT
      END IF  
    END DO
  END DO
END SUBROUTINE HashTableToList
   
END MODULE Chemsys1_Mod

PROGRAM TChemsys
  USE Chemsys1_Mod
  IMPLICIT NONE

  Filename='Salt'
  CALL getarg(1,FileName)
  CALL ReadSystem(FileName)
  CALL ReadThermoData('Thermo')
  OPEN(UNIT=89,FILE=ADJUSTL(TRIM(FileName))//'.chem',STATUS='UNKNOWN')
  CALL PrintHeadSpecies(89)
  CALL PrintSpecies(ListGas2,89)
  CALL PrintSpecies(ListAqua2,89)
  CALL PrintSpecies(ListSolid2,89)
  CALL PrintSpecies(ListPartic2,89)
  CALL PrintSpecies(ListNonReac2,89)
! WRITE(89,*) "'UCL'"
! WRITE(89,*) "'UCR'"
! WRITE(89,*) "'VCL'"
! WRITE(89,*) "'VCR'"
! WRITE(89,*) "'WCL'"
! WRITE(89,*) "'WCR'"
  CALL PrintHeadReactions(89)
  CALL PrintReaction(ListRGas,89)
  CALL PrintReaction(ListRHenry,89)
  CALL PrintReaction(ListRDiss,89)
  CALL PrintReaction(ListRAqua,89)
  CALL PrintReaction(ListRSolid,89)
  CALL PrintReaction(ListRPartic,89)
  CALL PrintReaction(ListRMicro,89)
  CLOSE(89)
END PROGRAM TChemsys
