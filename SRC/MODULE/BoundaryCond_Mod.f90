MODULE BoundaryCond_Mod

  USE Kind_Mod
  USE Control_Mod
  USE Physics_Mod
  USE Control_Mod,  ONLY: PrintNameLists

  IMPLICIT NONE

  TYPE BoundaryP_T
    INTEGER :: West
    INTEGER :: East
    INTEGER :: South
    INTEGER :: North
    INTEGER :: Bottom
    INTEGER :: Top
  END TYPE BoundaryP_T
!  INTEGER :: WestI
!  INTEGER :: EastI
!  INTEGER :: SouthI
!  INTEGER :: NorthI
!  INTEGER :: BottomI
!  INTEGER :: TopI

  TYPE BoundaryCon_T
    CHARACTER(20) :: West
    CHARACTER(20) :: East
    CHARACTER(20) :: South
    CHARACTER(20) :: North
    CHARACTER(20) :: Bottom
    CHARACTER(20) :: Top
  END TYPE BoundaryCon_T
!  CHARACTER(10) :: WestC
!  CHARACTER(10) :: EastC
!  CHARACTER(10) :: SouthC
!  CHARACTER(10) :: NorthC
!  CHARACTER(10) :: BottomC
!  CHARACTER(10) :: TopC

  TYPE (BoundaryP_T)   :: BCP
  TYPE (BoundaryCon_T) :: BCVel
  TYPE (BoundaryCon_T) :: BCScal
  TYPE (BoundaryCon_T),POINTER :: BCMetVec(:)
  TYPE (BoundaryCon_T),POINTER :: BCChemVec(:)

  NAMELIST /ModelBCP/   BCP
  NAMELIST /ModelBCVel/ BCVel
  NAMELIST /ModelBCu/   BCScal
  NAMELIST /ModelBCv/   BCScal
  NAMELIST /ModelBCw/   BCScal
  NAMELIST /ModelBCth/  BCScal
  NAMELIST /ModelBCRho/  BCScal
  NAMELIST /ModelBCtke/ BCScal
  NAMELIST /ModelBCdis/ BCScal
  NAMELIST /ModelBCtkeH/ BCScal
  NAMELIST /ModelBCtkeV/ BCScal
  NAMELIST /ModelBCLen/ BCScal
  NAMELIST /ModelBCEn/ BCScal
  NAMELIST /ModelBCRhoV/  BCScal
  NAMELIST /ModelBCRhoC/  BCScal
  NAMELIST /ModelBCRhoR/  BCScal
  NAMELIST /ModelBCRhoI/  BCScal
  NAMELIST /ModelBCgas/  BCScal ! Hinneburg: alle species

! NAMELIST /ModelBCP/ WestI,   &
!                     EastI,   &
!                     SouthI,  &
!                     NorthI,  &
!                     BottomI, &
!                     TopI
!  NAMELIST /ModelBCVel/ WestC,   &
!                        EastC,   &
!                        SouthC,  &
!                        NorthC,  &
!                        BottomC, &
!                        TopC
!  NAMELIST /ModelBCu/ WestC,   &
!                      EastC,   &
!                      SouthC,  &
!                      NorthC,  &
!                      BottomC, &
!                      TopC
!  NAMELIST /ModelBCv/ WestC,   &
!                      EastC,   &
!                      SouthC,  &
!                      NorthC,  &
!                      BottomC, &
!                      TopC
!  NAMELIST /ModelBCw/ WestC,   &
!                      EastC,   &
!                      SouthC,  &
!                      NorthC,  &
!                      BottomC, &
!                      TopC
! NAMELIST /ModelBCth/ WestC,   &
!                      EastC,   &
!                      SouthC,  &
!                      NorthC,  &
!                      BottomC, &
!                      TopC
!  NAMELIST /ModelBCtke/ WestC,   &
!                        EastC,   &
!                        SouthC,  &
!                        NorthC,  &
!                        BottomC, &
!                        TopC
!  NAMELIST /ModelBCdis/ WestC,   &
!                        EastC,   &
!                        SouthC,  &
!                        NorthC,  &
!                        BottomC, &
!                        TopC
!  NAMELIST /ModelBCtkeH/ WestC,   &
!                        EastC,   &
!                        SouthC,  &
!                        NorthC,  &
!                        BottomC, &
!                        TopC
!  NAMELIST /ModelBCtkeV/ WestC,   &
!                        EastC,   &
!                        SouthC,  &
!                        NorthC,  &
!                        BottomC, &
!                        TopC
!  NAMELIST /ModelBCLen/ WestC,   &
!                        EastC,   &
!                        SouthC,  &
!                        NorthC,  &
!                        BottomC, &
!                        TopC
!  NAMELIST /ModelBCEn/  WestC,   &
!                        EastC,   &
!                        SouthC,  &
!                        NorthC,  &
!                        BottomC, &
!                        TopC
!  NAMELIST /ModelBCRhoV/ WestC,   &
!                       EastC,   &
!                       SouthC,  &
!                       NorthC,  &
!                       BottomC, &
!                       TopC
!  NAMELIST /ModelBCRhoC/ WestC,   &
!                       EastC,   &
!                       SouthC,  &
!                       NorthC,  &
!                       BottomC, &
!                       TopC
!  NAMELIST /ModelBCRhoI/ WestC,   &
!                       EastC,   &
!                       SouthC,  &
!                       NorthC,  &
!                       BottomC, &
!                       TopC
!  NAMELIST /ModelBCRhoR/ WestC,   &
!                       EastC,   &
!                       SouthC,  &
!                       NorthC,  &
!                       BottomC, &
!                       TopC
! NAMELIST /ModelBCRho/ WestC,   &
!                      EastC,   &
!                      SouthC,  &
!                      NorthC,  &
!                      BottomC, &
!                      TopC


CONTAINS

SUBROUTINE InputModelBCVel(FileName)

  CHARACTER(*) :: FileName

  INTEGER :: ic
  CHARACTER(300) :: Line

! Additional boundary condition ! Hinneburg
! (homogenized periodic boundary condition using actual mean profile)
! Keywords: MeanFlow and MeanValue in Namelists ModelBC...
!
! BCVel     = InFlow(=Function), OutFlow(=ZeroGrad), MeanFlow, FreeSlip, or NoSlip
! MeanFlow  = OutFlow, but instead of ZeroGrad:
!               MeanValue for meteorological variables,
!               Function  for other variables
! MeanValue = value of mean profile, uniform for x and y, only dependent on z
!             (actual avarage of all cells with same z)

  BCVel%West  ='OutFlow'
  BCVel%East  ='OutFlow'
  BCVel%South ='OutFlow'
  BCVel%North ='OutFlow'
  BCVel%Bottom='OutFlow'
  BCVel%Top   ='OutFlow'

! Find line 'BCVel' (first appearance) and modify

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ModelBCVel')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelBCVel)
      !BCVel%West  =WestC
      !BCVel%East  =EastC
      !BCVel%South =SouthC
      !BCVel%North =NorthC
      !BCVel%Bottom=BottomC
      !BCVel%Top   =TopC
      EXIT
    END IF
  END DO
1 CONTINUE
  CLOSE(InputUnit)

END SUBROUTINE InputModelBCVel


SUBROUTINE InputModelBC(FileName)


  CHARACTER(*) :: FileName

  INTEGER :: ic
  CHARACTER(300) :: Line


  ALLOCATE(BCMetVec(VectorComponentsMet))
  ALLOCATE(BCChemVec(VectorComponentsChem))

! Standard BCVel, BCMetVec and BCP (corresponding to OutFlow)

  BCVel%West  ='OutFlow'
  BCVel%East  ='OutFlow'
  BCVel%South ='OutFlow'
  BCVel%North ='OutFlow'
  BCVel%Bottom='OutFlow'
  BCVel%Top   ='OutFlow'

  DO ic=1,VectorComponentsMet
    BCMetVec(ic)%West=  'ZeroGrad'
    BCMetVec(ic)%East=  'ZeroGrad'
    BCMetVec(ic)%South= 'ZeroGrad'
    BCMetVec(ic)%North= 'ZeroGrad'
    BCMetVec(ic)%Bottom='ZeroGrad'
    BCMetVec(ic)%Top=   'ZeroGrad'
  END DO
  DO ic=1,VectorComponentsChem
    BCChemVec(ic)%West=  'ZeroGrad'
    BCChemVec(ic)%East=  'ZeroGrad'
    BCChemVec(ic)%South= 'ZeroGrad'
    BCChemVec(ic)%North= 'ZeroGrad'
    BCChemVec(ic)%Bottom='ZeroGrad'
    BCChemVec(ic)%Top=   'ZeroGrad'
  END DO

  BCP%West  =1
  BCP%East  =1
  BCP%South =1
  BCP%North =1
  BCP%Bottom=1
  BCP%Top   =1

! Find line 'BCVel' (first appearance) and modify

  OPEN(UNIT=InputUnit,FILE=TRIM(FileName),STATUS='OLD')
  REWIND(InputUnit)
  DO
    READ(InputUnit,*,END=1) Line
    IF (INDEX(Line,'&ModelBCVel')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelBCVel)
      EXIT
    END IF
  END DO
1 CONTINUE

! Derive BCMetVec from (modified) BCVel ('Period' after all modifications)

  IF (uPosL*vPosL*wPosL>0) THEN
    IF (BCVel%West=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCMetVec(ic)%West='Function'
      END DO
    ELSE IF (BCVel%West=='MeanFlow') THEN ! Hinneburg
      DO ic=1,NumMet
        BCMetVec(ic)%West='MeanValue'
      END DO
      DO ic=NumMet+1,VectorComponentsM
        BCMetVec(ic)%West='Function'
      END DO
    ELSE IF (BCVel%West=='FreeSlip') THEN
      BCMetVec(uPosL)%West='ZeroValue'
      BCMetVec(uPosR)%West='ZeroValue'
    ELSE IF (BCVel%West=='NoSlip') THEN
      BCMetVec(uPosL)%West='ZeroValue'
      BCMetVec(uPosR)%West='ZeroValue'
      BCMetVec(vPosL)%West='ZeroValue'
      BCMetVec(vPosR)%West='ZeroValue'
      BCMetVec(wPosL)%West='ZeroValue'
      BCMetVec(wPosR)%West='ZeroValue'
    END IF

    IF (BCVel%East=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCMetVec(ic)%East='Function'
      END DO
    ELSE IF (BCVel%East=='MeanFlow') THEN ! Hinneburg
      DO ic=1,NumMet
        BCMetVec(ic)%East='MeanValue'
      END DO
      DO ic=NumMet+1,VectorComponentsM
        BCMetVec(ic)%East='Function'
      END DO
    ELSE IF (BCVel%East=='FreeSlip') THEN
      BCMetVec(uPosL)%East='ZeroValue'
      BCMetVec(uPosR)%East='ZeroValue'
    ELSE IF (BCVel%East=='NoSlip') THEN
      BCMetVec(uPosL)%East='ZeroValue'
      BCMetVec(uPosR)%East='ZeroValue'
      BCMetVec(vPosL)%East='ZeroValue'
      BCMetVec(vPosR)%East='ZeroValue'
      BCMetVec(wPosL)%East='ZeroValue'
      BCMetVec(wPosR)%East='ZeroValue'
    END IF
   
    IF (BCVel%South=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCMetVec(ic)%South='Function'
      END DO
    ELSE IF (BCVel%South=='MeanFlow') THEN ! Hinneburg
      DO ic=1,NumMet
        BCMetVec(ic)%South='MeanValue'
      END DO
      DO ic=NumMet+1,VectorComponentsM
        BCMetVec(ic)%South='Function'
      END DO
    ELSE IF (BCVel%South=='FreeSlip') THEN
      BCMetVec(vPosL)%South='ZeroValue'
      BCMetVec(vPosR)%South='ZeroValue'
    ELSE IF (BCVel%South=='NoSlip') THEN
      BCMetVec(uPosL)%South='ZeroValue'
      BCMetVec(uPosR)%South='ZeroValue'
      BCMetVec(vPosL)%South='ZeroValue'
      BCMetVec(vPosR)%South='ZeroValue'
      BCMetVec(wPosL)%South='ZeroValue'
      BCMetVec(wPosR)%South='ZeroValue'
    END IF

    IF (BCVel%North=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCMetVec(ic)%North='Function'
      END DO
    ELSE IF (BCVel%North=='MeanFlow') THEN ! Hinneburg
      DO ic=1,NumMet
        BCMetVec(ic)%North='MeanValue'
      END DO
      DO ic=NumMet+1,VectorComponentsM
        BCMetVec(ic)%North='Function'
      END DO
    ELSE IF (BCVel%North=='FreeSlip') THEN
      BCMetVec(vPosL)%North='ZeroValue'
      BCMetVec(vPosR)%North='ZeroValue'
    ELSE IF (BCVel%North=='NoSlip') THEN
      BCMetVec(uPosL)%North='ZeroValue'
      BCMetVec(uPosR)%North='ZeroValue'
      BCMetVec(vPosL)%North='ZeroValue'
      BCMetVec(vPosR)%North='ZeroValue'
      BCMetVec(wPosL)%North='ZeroValue'
      BCMetVec(wPosR)%North='ZeroValue'
    END IF

    IF (BCVel%Bottom=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCMetVec(ic)%Bottom='Function'
      END DO
    ELSE IF (BCVel%Bottom=='MeanFlow') THEN ! Hinneburg 
      BCVel%Bottom='OutFlow'
    ELSE IF (BCVel%Bottom=='FreeSlip') THEN
      BCMetVec(wPosL)%Bottom='ZeroValue'
      BCMetVec(wPosR)%Bottom='ZeroValue'
    ELSE IF (BCVel%Bottom=='NoSlip') THEN
      BCMetVec(uPosL)%Bottom='ZeroValue'
      BCMetVec(uPosR)%Bottom='ZeroValue'
      BCMetVec(vPosL)%Bottom='ZeroValue'
      BCMetVec(vPosR)%Bottom='ZeroValue'
      BCMetVec(wPosL)%Bottom='ZeroValue'
      BCMetVec(wPosR)%Bottom='ZeroValue'
    END IF

    IF (BCVel%Top=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCMetVec(ic)%Top='Function'
      END DO
    ELSE IF (BCVel%Top=='MeanFlow') THEN ! Hinneburg
      BCVel%Top='OutFlow' 
    ELSE IF (BCVel%Top=='FreeSlip') THEN
      BCMetVec(wPosL)%Top='ZeroValue'
      BCMetVec(wPosR)%Top='ZeroValue'
    ELSE IF (BCVel%Top=='NoSlip') THEN
      BCMetVec(uPosL)%Top='ZeroValue'
      BCMetVec(uPosR)%Top='ZeroValue'
      BCMetVec(vPosL)%Top='ZeroValue'
      BCMetVec(vPosR)%Top='ZeroValue'
      BCMetVec(wPosL)%Top='ZeroValue'
      BCMetVec(wPosR)%Top='ZeroValue'
    END IF
  END IF

! Find lines 'BCMetVec' (first appearance) and modify

  IF (uPosL>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=11) Line
      IF (INDEX(Line,'&ModelBCu')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(uPosL)
        READ(InputUnit,NML=ModelBCu)
        BCMetVec(uPosL)=BCScal
        BCMetVec(uPosR)=BCScal
        EXIT
      END IF
    END DO
11  CONTINUE
  END IF

  IF (vPosL>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=12) Line
      IF (INDEX(Line,'&ModelBCv')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(vPosL)
        READ(InputUnit,NML=ModelBCv)
        BCMetVec(vPosL)=BCScal
        BCMetVec(vPosR)=BCScal
        EXIT
      END IF
    END DO
12  CONTINUE
  END IF

  IF (wPosL>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=13) Line
      IF (INDEX(Line,'&ModelBCw')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(wPosL)
        READ(InputUnit,NML=ModelBCw)
        BCMetVec(wPosL)=BCScal
        BCMetVec(wPosR)=BCScal
        EXIT
      END IF
    END DO
13  CONTINUE
  END IF

  IF (thPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=14) Line
      IF (INDEX(Line,'&ModelBCth')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(thPos)
        READ(InputUnit,NML=ModelBCth)
        BCMetVec(thPos)=BCScal
        EXIT
      END IF
    END DO
14  CONTINUE
  END IF

  IF (tkePos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=15) Line
      IF (INDEX(Line,'&ModelBCtke')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(tkePos)
        READ(InputUnit,NML=ModelBCtke)
        BCMetVec(tkePos)=BCScal
        EXIT
      END IF
    END DO
15  CONTINUE
  END IF

  IF (disPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=16) Line
      IF (INDEX(Line,'&ModelBCdis')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(disPos)
        READ(InputUnit,NML=ModelBCdis)
        BCMetVec(disPos)=BCScal
        EXIT
      END IF
    END DO
16  CONTINUE
  END IF

  IF (RhoVPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=17) Line
      IF (INDEX(Line,'&ModelBCRhoV')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(RhoVPos)
        READ(InputUnit,NML=ModelBCRhoV)
        BCMetVec(RhoVPos)=BCScal
        EXIT
      END IF
    END DO
17  CONTINUE
  END IF

  IF (RhoCPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=18) Line
      IF (INDEX(Line,'&ModelBCRhoC')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(RhoCPos)
        READ(InputUnit,NML=ModelBCRhoC)
        BCMetVec(RhoCPos)=BCScal
        EXIT
      END IF
    END DO
18  CONTINUE
  END IF

  IF (RhoIPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=19) Line
      IF (INDEX(Line,'&ModelBCRhoI')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(RhoIPos)
        READ(InputUnit,NML=ModelBCRhoI)
        BCMetVec(RhoIPos)=BCScal
        EXIT
      END IF
    END DO
19  CONTINUE
  END IF

  IF (RhoRPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=20) Line
      IF (INDEX(Line,'&ModelBCRhoR')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(RhoRPos)
        READ(InputUnit,NML=ModelBCRhoR)
        BCMetVec(RhoRPos)=BCScal
        EXIT
      END IF
    END DO
20  CONTINUE
  END IF

  IF (RhoPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=21) Line
      IF (INDEX(Line,'&ModelBCRho')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(RhoPos)
        READ(InputUnit,NML=ModelBCRho)
        BCMetVec(RhoPos)=BCScal
        EXIT
      END IF
    END DO
21 CONTINUE
  END IF

  IF (tkeHPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=22) Line
      IF (INDEX(Line,'&ModelBCtkeH')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(tkeHPos)
        READ(InputUnit,NML=ModelBCtkeH)
        BCMetVec(tkeHPos)=BCScal
        EXIT
      END IF
    END DO
22  CONTINUE
  END IF

  IF (tkeVPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=23) Line
      IF (INDEX(Line,'&ModelBCtkeV')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(tkeVPos)
        READ(InputUnit,NML=ModelBCtkeV)
        BCMetVec(tkeVPos)=BCScal
        EXIT
      END IF
    END DO
23  CONTINUE
  END IF

  IF (LenPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=24) Line
      IF (INDEX(Line,'&ModelBCLen')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(LenPos)
        READ(InputUnit,NML=ModelBCLen)
        BCMetVec(LenPos)=BCScal
        EXIT
      END IF
    END DO
24  CONTINUE
  END IF

  IF (EnPos>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=25) Line
      IF (INDEX(Line,'&ModelBCEn')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCMetVec(EnPos)
        READ(InputUnit,NML=ModelBCEn)
        BCMetVec(EnPos)=BCScal
        EXIT
      END IF
    END DO
25  CONTINUE
  END IF

! Signing BCVel='MeanFlow' from BCMetVec
! (one 'MeanValue' --> 'MeanFlow') ! Hinneburg (for CALL MeanProfileCompute)

    IF (BCVel%West/='MeanFlow') THEN
      DO ic=1,VectorComponentsMet
        IF (BCMetVec(ic)%West=='MeanValue') THEN
          BCVel%West='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%East/='MeanFlow') THEN
      DO ic=1,VectorComponentsMet
        IF (BCMetVec(ic)%East=='MeanValue') THEN
          BCVel%East='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%South/='MeanFlow') THEN
      DO ic=1,VectorComponentsMet
        IF (BCMetVec(ic)%South=='MeanValue') THEN
          BCVel%South='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%North/='MeanFlow') THEN
      DO ic=1,VectorComponentsMet
        IF (BCMetVec(ic)%North=='MeanValue') THEN
          BCVel%North='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%Bottom/='MeanFlow') THEN
      DO ic=1,VectorComponentsMet
        IF (BCMetVec(ic)%Bottom=='MeanValue') THEN
          BCVel%Bottom='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%Top/='MeanFlow') THEN
      DO ic=1,VectorComponentsMet
        IF (BCMetVec(ic)%Top=='MeanValue') THEN
          BCVel%Top='MeanFlow'
          EXIT
        END IF
      END DO
    END IF

! Fixing of all BCMetVec=Period in case of BCVel=Period

  IF (BCVel%West=='Period') THEN
    DO ic=1,VectorComponentsMet
      BCMetVec(ic)%West='Period'
    END DO
  END IF
  IF (BCVel%East=='Period') THEN
    DO ic=1,VectorComponentsMet
      BCMetVec(ic)%East='Period'
    END DO
  END IF
  IF (BCVel%South=='Period') THEN
    DO ic=1,VectorComponentsMet
      BCMetVec(ic)%South='Period'
    END DO
  END IF
  IF (BCVel%North=='Period') THEN
    DO ic=1,VectorComponentsMet
      BCMetVec(ic)%North='Period'
    END DO
  END IF
  IF (BCVel%Bottom=='Period') THEN
    DO ic=1,VectorComponentsMet
      BCMetVec(ic)%Bottom='Period'
    END DO
  END IF
  IF (BCVel%Top=='Period') THEN
    DO ic=1,VectorComponentsMet
      BCMetVec(ic)%Top='Period'
    END DO
  END IF

! Derive BCP from (modified) BCMetVec

  IF (uPosL*vPosL*wPosL>0) THEN
    IF (BCMetVec(uPosL)%West  .NE.'ZeroGrad' .AND. BCMetVec(uPosL)%West  .NE.'MeanValue') THEN ! Hinneburg
      BCP%West  =-1
    END IF
    IF (BCMetVec(uPosL)%East  .NE.'ZeroGrad' .AND. BCMetVec(uPosL)%East  .NE.'MeanValue') THEN
      BCP%East  =-1
    END IF
    IF (BCMetVec(vPosL)%South .NE.'ZeroGrad' .AND. BCMetVec(vPosL)%South .NE.'MeanValue') THEN
      BCP%South =-1
    END IF
    IF (BCMetVec(vPosL)%North .NE.'ZeroGrad' .AND. BCMetVec(vPosL)%North .NE.'MeanValue') THEN
      BCP%North =-1
    END IF
    IF (BCMetVec(wPosL)%Bottom.NE.'ZeroGrad' .AND. BCMetVec(wPosL)%Bottom.NE.'MeanValue') THEN
      BCP%Bottom=-1
    END IF
    IF (BCMetVec(wPosL)%Top   .NE.'ZeroGrad' .AND. BCMetVec(wPosL)%Top   .NE.'MeanValue') THEN
      BCP%Top   =-1
    END IF
  END IF

! Find line 'BCP' (first appearance) and modify

  REWIND(InputUnit)
  DO
    READ(InputUnit,*,END=31) Line
    IF (INDEX(Line,'&ModelBCP')>0) THEN
      BACKSPACE(InputUnit)
      READ(InputUnit,NML=ModelBCP)
      EXIT
    END IF
  END DO
31 CONTINUE
  IF (VectorComponentsChem>0) THEN ! Hinneburg
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=32) Line
      IF (INDEX(Line,'&ModelBCGas')>0) THEN ! Hinneburg
        BACKSPACE(InputUnit)
        READ(InputUnit,NML=ModelBCgas)
        DO ic=1,VectorComponentsChem
          BCChemVec(ic)=BCScal
        END DO
        EXIT
      END IF
    END DO
32  CONTINUE
  END IF
  CLOSE(InputUnit)

  IF (MyId==0.AND.PrintNameLists) THEN
    WRITE(TermUnit,*) '         ','WEST      EAST      ','  ' &
                          ,'SOUTH     NORTH     ','  ' &
                          ,'BOTTOM    TOP       '
    WRITE(TermUnit,*)
    IF (  uPosL.GT. 0) THEN
      WRITE(TermUnit,*) '  uPosL  ',BCMetVec(  uPosL)%West  ,BCMetVec(  uPosL)%East  ,'  ' &
                            ,BCMetVec(  uPosL)%South ,BCMetVec(  uPosL)%North ,'  ' &
                            ,BCMetVec(  uPosL)%Bottom,BCMetVec(  uPosL)%Top
    END IF
    IF (  vPosL .GT. 0) THEN
      WRITE(TermUnit,*) '  vPosL   ',BCMetVec(  vPosL)%West  ,BCMetVec(  vPosL)%East  ,'  ' &
                            ,BCMetVec(  vPosL)%South ,BCMetVec(  vPosL)%North ,'  ' &
                            ,BCMetVec(  vPosL)%Bottom,BCMetVec(  vPosL)%Top
    END IF
    IF (  wPosL .GT. 0) THEN
      WRITE(TermUnit,*) '  wPosL   ',BCMetVec(  wPosL)%West  ,BCMetVec(  wPosL)%East  ,'  ' &
                            ,BCMetVec(  wPosL)%South ,BCMetVec(  wPosL)%North ,'  ' &
                            ,BCMetVec(  wPosL)%Bottom,BCMetVec(  wPosL)%Top
    END IF
    WRITE(TermUnit,*)
    IF ( thPos .GT. 0) THEN
      WRITE(TermUnit,*) ' thPos   ',BCMetVec( thPos)%West  ,BCMetVec( thPos)%East  ,'  ' &
                            ,BCMetVec( thPos)%South ,BCMetVec( thPos)%North ,'  ' &
                            ,BCMetVec( thPos)%Bottom,BCMetVec( thPos)%Top
    END IF
    IF (tkePos .GT. 0) THEN
      WRITE(TermUnit,*) 'tkePos   ',BCMetVec(tkePos)%West  ,BCMetVec(tkePos)%East  ,'  ' &
                            ,BCMetVec(tkePos)%South ,BCMetVec(tkePos)%North ,'  ' &
                            ,BCMetVec(tkePos)%Bottom,BCMetVec(tkePos)%Top
    END IF
    IF (disPos .GT. 0) THEN
      WRITE(TermUnit,*) 'disPos   ',BCMetVec(disPos)%West  ,BCMetVec(disPos)%East  ,'  ' &
                            ,BCMetVec(disPos)%South ,BCMetVec(disPos)%North ,'  ' &
                            ,BCMetVec(disPos)%Bottom,BCMetVec(disPos)%Top
    END IF
    IF (tkeHPos .GT. 0) THEN
      WRITE(TermUnit,*) 'tkeHPos  ',BCMetVec(tkeHPos)%West  ,BCMetVec(tkeHPos)%East  ,'  ' &
                            ,BCMetVec(tkeHPos)%South ,BCMetVec(tkeHPos)%North ,'  ' &
                            ,BCMetVec(tkeHPos)%Bottom,BCMetVec(tkeHPos)%Top
    END IF
    IF (tkeVPos .GT. 0) THEN
      WRITE(TermUnit,*) 'tkeVPos  ',BCMetVec(tkeVPos)%West  ,BCMetVec(tkeVPos)%East  ,'  ' &
                            ,BCMetVec(tkeVPos)%South ,BCMetVec(tkeVPos)%North ,'  ' &
                            ,BCMetVec(tkeVPos)%Bottom,BCMetVec(tkeVPos)%Top
    END IF
    IF (LenPos .GT. 0) THEN
      WRITE(TermUnit,*) 'LenPos   ',BCMetVec(LenPos)%West  ,BCMetVec(LenPos)%East  ,'  ' &
                            ,BCMetVec(LenPos)%South ,BCMetVec(LenPos)%North ,'  ' &
                            ,BCMetVec(LenPos)%Bottom,BCMetVec(LenPos)%Top
    END IF
    IF (EnPos .GT. 0) THEN
      WRITE(TermUnit,*) 'EnPos   ',BCMetVec(EnPos)%West  ,BCMetVec(EnPos)%East  ,'  ' &
                            ,BCMetVec(EnPos)%South ,BCMetVec(EnPos)%North ,'  ' &
                            ,BCMetVec(EnPos)%Bottom,BCMetVec(EnPos)%Top
    END IF
    IF ( RhoVPos .GT. 0) THEN
      WRITE(TermUnit,*) ' RhoVPos   ',BCMetVec( RhoVPos)%West  ,BCMetVec( RhoVPos)%East  ,'  ' &
                            ,BCMetVec( RhoVPos)%South ,BCMetVec( RhoVPos)%North ,'  ' &
                            ,BCMetVec( RhoVPos)%Bottom,BCMetVec( RhoVPos)%Top
    END IF
    IF ( RhoCPos .GT. 0) THEN
      WRITE(TermUnit,*) ' RhoCPos   ',BCMetVec( RhoCPos)%West  ,BCMetVec( RhoCPos)%East  ,'  ' &
                            ,BCMetVec( RhoCPos)%South ,BCMetVec( RhoCPos)%North ,'  ' &
                            ,BCMetVec( RhoCPos)%Bottom,BCMetVec( RhoCPos)%Top
    END IF
    IF ( RhoIPos .GT. 0) THEN
      WRITE(TermUnit,*) ' RhoIPos   ',BCMetVec( RhoIPos)%West  ,BCMetVec( RhoIPos)%East  ,'  ' &
                            ,BCMetVec( RhoIPos)%South ,BCMetVec( RhoIPos)%North ,'  ' &
                            ,BCMetVec( RhoIPos)%Bottom,BCMetVec( RhoIPos)%Top
    END IF
    IF ( RhoRPos .GT. 0) THEN
      WRITE(TermUnit,*) ' RhoRPos   ',BCMetVec( RhoRPos)%West  ,BCMetVec( RhoRPos)%East  ,'  ' &
                            ,BCMetVec( RhoRPos)%South ,BCMetVec( RhoRPos)%North ,'  ' &
                            ,BCMetVec( RhoRPos)%Bottom,BCMetVec( RhoRPos)%Top
    END IF
    IF (RhoPos .GT. 0) THEN
      WRITE(TermUnit,*) 'RhoPos   ',BCMetVec(RhoPos)%West  ,BCMetVec(RhoPos)%East  ,'  ' &
                            ,BCMetVec(RhoPos)%South ,BCMetVec(RhoPos)%North ,'  ' &
                            ,BCMetVec(RhoPos)%Bottom,BCMetVec(RhoPos)%Top
    END IF
    DO ic=1,UBOUND(BCChemVec,1)
      WRITE(TermUnit,*)
      WRITE(TermUnit,*) 'gasPos   ',BCChemVec(ic)%West  ,BCChemVec(ic)%East  ,'  ' &
                                   ,BCChemVec(ic)%South ,BCChemVec(ic)%North ,'  ' &
                                   ,BCChemVec(ic)%Bottom,BCChemVec(ic)%Top
    END DO
    WRITE(TermUnit,*)
    WRITE(TermUnit,40) BCP%West,BCP%East,BCP%South,BCP%North,BCP%Bottom,BCP%Top
    WRITE(TermUnit,*)
  END IF
40 FORMAT(' P (boundaryCond_mod)        ',3(2(i2,8x),2x))

END SUBROUTINE InputModelBC

END MODULE BoundaryCond_Mod
