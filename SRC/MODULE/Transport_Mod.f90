MODULE Transport_Mod

  USE Kind_Mod
  USE Control_Mod

  IMPLICIT NONE
  INTEGER :: it0,it1
  INTEGER :: nFrac
  REAL(RealKind) :: rStart
  REAL(RealKind) :: pFac
  REAL(RealKind) :: pFacInit
  REAL(RealKind) :: Numin
  REAL(RealKind) :: bMax
  REAL(RealKind) :: RelFac1,RelFac2
  REAL(RealKind) :: MassEps
  REAL(RealKind) :: EmissionTimeStart ! Hinneburg
  CHARACTER*40 :: ChemieFile
  CHARACTER*40 :: ChemOutFile
  CHARACTER*40 :: TrajFile
  CHARACTER*40 :: DataFile
  CHARACTER*40 :: IniFile
  CHARACTER*20 :: GasUnit
  LOGICAL :: Neutral
  LOGICAL :: Chemie
  LOGICAL :: ChemieGas
  LOGICAL :: ChemieAqua
  LOGICAL :: Aerosol
  LOGICAL :: Koagul
  LOGICAL :: Condens
  LOGICAL :: Depos
  LOGICAL :: Emiss
  LOGICAL :: SeaEmiss
  LOGICAL :: EmissStreet ! Hinneburg
  LOGICAL :: AdvMass
  LOGICAL :: HenryTrans
  LOGICAL :: Limit
  LOGICAL :: Dry
  CHARACTER*4 :: MethodAct='Aim'
  NAMELIST /ModelTransport/ & 
                   it0, &
                   it1, & 
                   nFrac,  &
                   rStart, &
                   pFac,   &
                   pFacInit,   &
                   Numin, &
                   bMax, &
                   MassEps, &
                   RelFac1, &
                   RelFac2, &
                   ChemieFile,&
                   ChemOutFile,&
                   TrajFile,&
                   DataFile,&
                   IniFile,    &
                   GasUnit, &
                   Neutral, &
                   Chemie, &
                   ChemieGas, &
                   ChemieAqua, &
                   Aerosol, &
                   Koagul, &
                   Condens, &
                   Depos, &
                   Emiss, &
                   SeaEmiss, &
                   EmissStreet, & ! Hinneburg
                   EmissionTimeStart, & ! Hinneburg
                   AdvMass, &
                   HenryTrans, &
                   Limit, &
                   MethodAct,    &
                   Dry
CONTAINS

SUBROUTINE InputModelTransport(FileName)

  CHARACTER(*) :: FileName
  
  INTEGER :: Pos
  CHARACTER(300) :: Line

  it0=0
  it1=1
  NuMin=1.d-10
  ChemieFile=''
  ChemOutFile='Output'
  TrajFile=''
  DataFile=''
  IniFile=''
  GasUnit=''
  Neutral =.FALSE.
  Chemie =.FALSE.
  ChemieGas =.FALSE.
  ChemieAqua =.FALSE.
  Aerosol =.FALSE.
  Koagul =.FALSE.
  Condens =.FALSE.
  Depos =.FALSE.
  Emiss =.FALSE.
  SeaEmiss =.FALSE.
  EmissStreet =.FALSE. ! Hinneburg
  EmissionTimeStart=StartTime ! Hinneburg
  AdvMass =.FALSE.
  HenryTrans =.FALSE.
  Limit =.FALSE.
  Dry =.FALSE.

! Find line
  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ModelTransport')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelTransport)
      EXIT
    END IF
  END DO
  CLOSE(UNIT=InputUnit)
1 CONTINUE

  CLOSE(UNIT=InputUnit)

END SUBROUTINE InputModelTransport

END MODULE Transport_Mod

