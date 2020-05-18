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
  TYPE (BoundaryCon_T),POINTER :: BCVec(:)

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


  ALLOCATE(BCVec(VectorComponentsM))

! Standard BCVel, BCVec and BCP (corresponding to OutFlow)

  BCVel%West  ='OutFlow'
  BCVel%East  ='OutFlow'
  BCVel%South ='OutFlow'
  BCVel%North ='OutFlow'
  BCVel%Bottom='OutFlow'
  BCVel%Top   ='OutFlow'

  DO ic=1,VectorComponentsM
    BCVec(ic)%West=  'ZeroGrad'
    BCVec(ic)%East=  'ZeroGrad'
    BCVec(ic)%South= 'ZeroGrad'
    BCVec(ic)%North= 'ZeroGrad'
    BCVec(ic)%Bottom='ZeroGrad'
    BCVec(ic)%Top=   'ZeroGrad'
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

! Derive BCVec from (modified) BCVel ('Period' after all modifications)

  IF (uPosL*vPosL*wPosL>0) THEN
    IF (BCVel%West=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCVec(ic)%West='Function'
      END DO
    ELSE IF (BCVel%West=='MeanFlow') THEN ! Hinneburg
      DO ic=1,NumMet
        BCVec(ic)%West='MeanValue'
      END DO
      DO ic=NumMet+1,VectorComponentsM
        BCVec(ic)%West='Function'
      END DO
    ELSE IF (BCVel%West=='FreeSlip') THEN
      BCVec(uPosL)%West='ZeroValue'
      BCVec(uPosR)%West='ZeroValue'
    ELSE IF (BCVel%West=='NoSlip') THEN
      BCVec(uPosL)%West='ZeroValue'
      BCVec(uPosR)%West='ZeroValue'
      BCVec(vPosL)%West='ZeroValue'
      BCVec(vPosR)%West='ZeroValue'
      BCVec(wPosL)%West='ZeroValue'
      BCVec(wPosR)%West='ZeroValue'
    END IF

    IF (BCVel%East=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCVec(ic)%East='Function'
      END DO
    ELSE IF (BCVel%East=='MeanFlow') THEN ! Hinneburg
      DO ic=1,NumMet
        BCVec(ic)%East='MeanValue'
      END DO
      DO ic=NumMet+1,VectorComponentsM
        BCVec(ic)%East='Function'
      END DO
    ELSE IF (BCVel%East=='FreeSlip') THEN
      BCVec(uPosL)%East='ZeroValue'
      BCVec(uPosR)%East='ZeroValue'
    ELSE IF (BCVel%East=='NoSlip') THEN
      BCVec(uPosL)%East='ZeroValue'
      BCVec(uPosR)%East='ZeroValue'
      BCVec(vPosL)%East='ZeroValue'
      BCVec(vPosR)%East='ZeroValue'
      BCVec(wPosL)%East='ZeroValue'
      BCVec(wPosR)%East='ZeroValue'
    END IF
   
    IF (BCVel%South=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCVec(ic)%South='Function'
      END DO
    ELSE IF (BCVel%South=='MeanFlow') THEN ! Hinneburg
      DO ic=1,NumMet
        BCVec(ic)%South='MeanValue'
      END DO
      DO ic=NumMet+1,VectorComponentsM
        BCVec(ic)%South='Function'
      END DO
    ELSE IF (BCVel%South=='FreeSlip') THEN
      BCVec(vPosL)%South='ZeroValue'
      BCVec(vPosR)%South='ZeroValue'
    ELSE IF (BCVel%South=='NoSlip') THEN
      BCVec(uPosL)%South='ZeroValue'
      BCVec(uPosR)%South='ZeroValue'
      BCVec(vPosL)%South='ZeroValue'
      BCVec(vPosR)%South='ZeroValue'
      BCVec(wPosL)%South='ZeroValue'
      BCVec(wPosR)%South='ZeroValue'
    END IF

    IF (BCVel%North=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCVec(ic)%North='Function'
      END DO
    ELSE IF (BCVel%North=='MeanFlow') THEN ! Hinneburg
      DO ic=1,NumMet
        BCVec(ic)%North='MeanValue'
      END DO
      DO ic=NumMet+1,VectorComponentsM
        BCVec(ic)%North='Function'
      END DO
    ELSE IF (BCVel%North=='FreeSlip') THEN
      BCVec(vPosL)%North='ZeroValue'
      BCVec(vPosR)%North='ZeroValue'
    ELSE IF (BCVel%North=='NoSlip') THEN
      BCVec(uPosL)%North='ZeroValue'
      BCVec(uPosR)%North='ZeroValue'
      BCVec(vPosL)%North='ZeroValue'
      BCVec(vPosR)%North='ZeroValue'
      BCVec(wPosL)%North='ZeroValue'
      BCVec(wPosR)%North='ZeroValue'
    END IF

    IF (BCVel%Bottom=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCVec(ic)%Bottom='Function'
      END DO
    ELSE IF (BCVel%Bottom=='MeanFlow') THEN ! Hinneburg 
      BCVel%Bottom='OutFlow'
    ELSE IF (BCVel%Bottom=='FreeSlip') THEN
      BCVec(wPosL)%Bottom='ZeroValue'
      BCVec(wPosR)%Bottom='ZeroValue'
    ELSE IF (BCVel%Bottom=='NoSlip') THEN
      BCVec(uPosL)%Bottom='ZeroValue'
      BCVec(uPosR)%Bottom='ZeroValue'
      BCVec(vPosL)%Bottom='ZeroValue'
      BCVec(vPosR)%Bottom='ZeroValue'
      BCVec(wPosL)%Bottom='ZeroValue'
      BCVec(wPosR)%Bottom='ZeroValue'
    END IF

    IF (BCVel%Top=='InFlow') THEN
      DO ic=1,VectorComponentsM
        BCVec(ic)%Top='Function'
      END DO
    ELSE IF (BCVel%Top=='MeanFlow') THEN ! Hinneburg
      BCVel%Top='OutFlow' 
    ELSE IF (BCVel%Top=='FreeSlip') THEN
      BCVec(wPosL)%Top='ZeroValue'
      BCVec(wPosR)%Top='ZeroValue'
    ELSE IF (BCVel%Top=='NoSlip') THEN
      BCVec(uPosL)%Top='ZeroValue'
      BCVec(uPosR)%Top='ZeroValue'
      BCVec(vPosL)%Top='ZeroValue'
      BCVec(vPosR)%Top='ZeroValue'
      BCVec(wPosL)%Top='ZeroValue'
      BCVec(wPosR)%Top='ZeroValue'
    END IF
  END IF

! Find lines 'BCVec' (first appearance) and modify

  IF (uPosL>0) THEN
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=11) Line
      IF (INDEX(Line,'&ModelBCu')>0) THEN
        BACKSPACE(InputUnit)
        BCScal=BCVec(uPosL)
        READ(InputUnit,NML=ModelBCu)
        BCVec(uPosL)=BCScal
        BCVec(uPosR)=BCScal
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
        BCScal=BCVec(vPosL)
        READ(InputUnit,NML=ModelBCv)
        BCVec(vPosL)=BCScal
        BCVec(vPosR)=BCScal
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
        BCScal=BCVec(wPosL)
        READ(InputUnit,NML=ModelBCw)
        BCVec(wPosL)=BCScal
        BCVec(wPosR)=BCScal
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
        BCScal=BCVec(thPos)
        READ(InputUnit,NML=ModelBCth)
        BCVec(thPos)=BCScal
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
        BCScal=BCVec(tkePos)
        READ(InputUnit,NML=ModelBCtke)
        BCVec(tkePos)=BCScal
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
        BCScal=BCVec(disPos)
        READ(InputUnit,NML=ModelBCdis)
        BCVec(disPos)=BCScal
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
        BCScal=BCVec(RhoVPos)
        READ(InputUnit,NML=ModelBCRhoV)
        BCVec(RhoVPos)=BCScal
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
        BCScal=BCVec(RhoCPos)
        READ(InputUnit,NML=ModelBCRhoC)
        BCVec(RhoCPos)=BCScal
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
        BCScal=BCVec(RhoIPos)
        READ(InputUnit,NML=ModelBCRhoI)
        BCVec(RhoIPos)=BCScal
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
        BCScal=BCVec(RhoRPos)
        READ(InputUnit,NML=ModelBCRhoR)
        BCVec(RhoRPos)=BCScal
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
        BCScal=BCVec(RhoPos)
        READ(InputUnit,NML=ModelBCRho)
        BCVec(RhoPos)=BCScal
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
        BCScal=BCVec(tkeHPos)
        READ(InputUnit,NML=ModelBCtkeH)
        BCVec(tkeHPos)=BCScal
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
        BCScal=BCVec(tkeVPos)
        READ(InputUnit,NML=ModelBCtkeV)
        BCVec(tkeVPos)=BCScal
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
        BCScal=BCVec(LenPos)
        READ(InputUnit,NML=ModelBCLen)
        BCVec(LenPos)=BCScal
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
        BCScal=BCVec(EnPos)
        READ(InputUnit,NML=ModelBCEn)
        BCVec(EnPos)=BCScal
        EXIT
      END IF
    END DO
25  CONTINUE
  END IF

  IF (VectorComponentsM>NumMet) THEN ! Hinneburg
    REWIND(InputUnit)
    DO
      READ(InputUnit,*,END=26) Line
      IF (INDEX(Line,'&ModelBCgas')>0) THEN ! Hinneburg
        BACKSPACE(InputUnit)
        BCScal=BCVec(VectorComponentsM)
        READ(InputUnit,NML=ModelBCgas)
        DO ic=NumMet+1,VectorComponentsM
          BCVec(ic)=BCScal
        END DO
        EXIT
      END IF
    END DO
26  CONTINUE
  END IF

! Signing BCVel='MeanFlow' from BCVec
! (one 'MeanValue' --> 'MeanFlow') ! Hinneburg (for CALL MeanProfileCompute)

    IF (BCVel%West/='MeanFlow') THEN
      DO ic=1,VectorComponentsM
        IF (BCVec(ic)%West=='MeanValue') THEN
          BCVel%West='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%East/='MeanFlow') THEN
      DO ic=1,VectorComponentsM
        IF (BCVec(ic)%East=='MeanValue') THEN
          BCVel%East='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%South/='MeanFlow') THEN
      DO ic=1,VectorComponentsM
        IF (BCVec(ic)%South=='MeanValue') THEN
          BCVel%South='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%North/='MeanFlow') THEN
      DO ic=1,VectorComponentsM
        IF (BCVec(ic)%North=='MeanValue') THEN
          BCVel%North='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%Bottom/='MeanFlow') THEN
      DO ic=1,VectorComponentsM
        IF (BCVec(ic)%Bottom=='MeanValue') THEN
          BCVel%Bottom='MeanFlow'
          EXIT
        END IF
      END DO
    END IF
    IF (BCVel%Top/='MeanFlow') THEN
      DO ic=1,VectorComponentsM
        IF (BCVec(ic)%Top=='MeanValue') THEN
          BCVel%Top='MeanFlow'
          EXIT
        END IF
      END DO
    END IF

! Fixing of all BCVec=Period in case of BCVel=Period

  IF (BCVel%West=='Period') THEN
    DO ic=1,VectorComponentsM
      BCVec(ic)%West='Period'
    END DO
  END IF
  IF (BCVel%East=='Period') THEN
    DO ic=1,VectorComponentsM
      BCVec(ic)%East='Period'
    END DO
  END IF
  IF (BCVel%South=='Period') THEN
    DO ic=1,VectorComponentsM
      BCVec(ic)%South='Period'
    END DO
  END IF
  IF (BCVel%North=='Period') THEN
    DO ic=1,VectorComponentsM
      BCVec(ic)%North='Period'
    END DO
  END IF
  IF (BCVel%Bottom=='Period') THEN
    DO ic=1,VectorComponentsM
      BCVec(ic)%Bottom='Period'
    END DO
  END IF
  IF (BCVel%Top=='Period') THEN
    DO ic=1,VectorComponentsM
      BCVec(ic)%Top='Period'
    END DO
  END IF

! Derive BCP from (modified) BCVec

  IF (uPosL*vPosL*wPosL>0) THEN
    IF (BCVec(uPosL)%West  .NE.'ZeroGrad' .AND. BCVec(uPosL)%West  .NE.'MeanValue') THEN ! Hinneburg
      BCP%West  =-1
    END IF
    IF (BCVec(uPosL)%East  .NE.'ZeroGrad' .AND. BCVec(uPosL)%East  .NE.'MeanValue') THEN
      BCP%East  =-1
    END IF
    IF (BCVec(vPosL)%South .NE.'ZeroGrad' .AND. BCVec(vPosL)%South .NE.'MeanValue') THEN
      BCP%South =-1
    END IF
    IF (BCVec(vPosL)%North .NE.'ZeroGrad' .AND. BCVec(vPosL)%North .NE.'MeanValue') THEN
      BCP%North =-1
    END IF
    IF (BCVec(wPosL)%Bottom.NE.'ZeroGrad' .AND. BCVec(wPosL)%Bottom.NE.'MeanValue') THEN
      BCP%Bottom=-1
    END IF
    IF (BCVec(wPosL)%Top   .NE.'ZeroGrad' .AND. BCVec(wPosL)%Top   .NE.'MeanValue') THEN
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
  CLOSE(InputUnit)

  IF (MyId==0.AND.PrintNameLists) THEN
    WRITE(TermUnit,*) '         ','WEST      EAST      ','  ' &
                          ,'SOUTH     NORTH     ','  ' &
                          ,'BOTTOM    TOP       '
    WRITE(TermUnit,*)
    IF (  uPosL.GT. 0) THEN
      WRITE(TermUnit,*) '  uPosL  ',BCVec(  uPosL)%West  ,BCVec(  uPosL)%East  ,'  ' &
                            ,BCVec(  uPosL)%South ,BCVec(  uPosL)%North ,'  ' &
                            ,BCVec(  uPosL)%Bottom,BCVec(  uPosL)%Top
    END IF
    IF (  vPosL .GT. 0) THEN
      WRITE(TermUnit,*) '  vPosL   ',BCVec(  vPosL)%West  ,BCVec(  vPosL)%East  ,'  ' &
                            ,BCVec(  vPosL)%South ,BCVec(  vPosL)%North ,'  ' &
                            ,BCVec(  vPosL)%Bottom,BCVec(  vPosL)%Top
    END IF
    IF (  wPosL .GT. 0) THEN
      WRITE(TermUnit,*) '  wPosL   ',BCVec(  wPosL)%West  ,BCVec(  wPosL)%East  ,'  ' &
                            ,BCVec(  wPosL)%South ,BCVec(  wPosL)%North ,'  ' &
                            ,BCVec(  wPosL)%Bottom,BCVec(  wPosL)%Top
    END IF
    WRITE(TermUnit,*)
    IF ( thPos .GT. 0) THEN
      WRITE(TermUnit,*) ' thPos   ',BCVec( thPos)%West  ,BCVec( thPos)%East  ,'  ' &
                            ,BCVec( thPos)%South ,BCVec( thPos)%North ,'  ' &
                            ,BCVec( thPos)%Bottom,BCVec( thPos)%Top
    END IF
    IF (tkePos .GT. 0) THEN
      WRITE(TermUnit,*) 'tkePos   ',BCVec(tkePos)%West  ,BCVec(tkePos)%East  ,'  ' &
                            ,BCVec(tkePos)%South ,BCVec(tkePos)%North ,'  ' &
                            ,BCVec(tkePos)%Bottom,BCVec(tkePos)%Top
    END IF
    IF (disPos .GT. 0) THEN
      WRITE(TermUnit,*) 'disPos   ',BCVec(disPos)%West  ,BCVec(disPos)%East  ,'  ' &
                            ,BCVec(disPos)%South ,BCVec(disPos)%North ,'  ' &
                            ,BCVec(disPos)%Bottom,BCVec(disPos)%Top
    END IF
    IF (tkeHPos .GT. 0) THEN
      WRITE(TermUnit,*) 'tkeHPos  ',BCVec(tkeHPos)%West  ,BCVec(tkeHPos)%East  ,'  ' &
                            ,BCVec(tkeHPos)%South ,BCVec(tkeHPos)%North ,'  ' &
                            ,BCVec(tkeHPos)%Bottom,BCVec(tkeHPos)%Top
    END IF
    IF (tkeVPos .GT. 0) THEN
      WRITE(TermUnit,*) 'tkeVPos  ',BCVec(tkeVPos)%West  ,BCVec(tkeVPos)%East  ,'  ' &
                            ,BCVec(tkeVPos)%South ,BCVec(tkeVPos)%North ,'  ' &
                            ,BCVec(tkeVPos)%Bottom,BCVec(tkeVPos)%Top
    END IF
    IF (LenPos .GT. 0) THEN
      WRITE(TermUnit,*) 'LenPos   ',BCVec(LenPos)%West  ,BCVec(LenPos)%East  ,'  ' &
                            ,BCVec(LenPos)%South ,BCVec(LenPos)%North ,'  ' &
                            ,BCVec(LenPos)%Bottom,BCVec(LenPos)%Top
    END IF
    IF (EnPos .GT. 0) THEN
      WRITE(TermUnit,*) 'EnPos   ',BCVec(EnPos)%West  ,BCVec(EnPos)%East  ,'  ' &
                            ,BCVec(EnPos)%South ,BCVec(EnPos)%North ,'  ' &
                            ,BCVec(EnPos)%Bottom,BCVec(EnPos)%Top
    END IF
    IF ( RhoVPos .GT. 0) THEN
      WRITE(TermUnit,*) ' RhoVPos   ',BCVec( RhoVPos)%West  ,BCVec( RhoVPos)%East  ,'  ' &
                            ,BCVec( RhoVPos)%South ,BCVec( RhoVPos)%North ,'  ' &
                            ,BCVec( RhoVPos)%Bottom,BCVec( RhoVPos)%Top
    END IF
    IF ( RhoCPos .GT. 0) THEN
      WRITE(TermUnit,*) ' RhoCPos   ',BCVec( RhoCPos)%West  ,BCVec( RhoCPos)%East  ,'  ' &
                            ,BCVec( RhoCPos)%South ,BCVec( RhoCPos)%North ,'  ' &
                            ,BCVec( RhoCPos)%Bottom,BCVec( RhoCPos)%Top
    END IF
    IF ( RhoIPos .GT. 0) THEN
      WRITE(TermUnit,*) ' RhoIPos   ',BCVec( RhoIPos)%West  ,BCVec( RhoIPos)%East  ,'  ' &
                            ,BCVec( RhoIPos)%South ,BCVec( RhoIPos)%North ,'  ' &
                            ,BCVec( RhoIPos)%Bottom,BCVec( RhoIPos)%Top
    END IF
    IF ( RhoRPos .GT. 0) THEN
      WRITE(TermUnit,*) ' RhoRPos   ',BCVec( RhoRPos)%West  ,BCVec( RhoRPos)%East  ,'  ' &
                            ,BCVec( RhoRPos)%South ,BCVec( RhoRPos)%North ,'  ' &
                            ,BCVec( RhoRPos)%Bottom,BCVec( RhoRPos)%Top
    END IF
    IF (RhoPos .GT. 0) THEN
      WRITE(TermUnit,*) 'RhoPos   ',BCVec(RhoPos)%West  ,BCVec(RhoPos)%East  ,'  ' &
                            ,BCVec(RhoPos)%South ,BCVec(RhoPos)%North ,'  ' &
                            ,BCVec(RhoPos)%Bottom,BCVec(RhoPos)%Top
    END IF
    IF (VectorComponentsM .GT. NumMet) THEN ! Hinneburg
      WRITE(TermUnit,*)
      WRITE(TermUnit,*) 'gasPos   ',BCVec(VectorComponentsM)%West  ,BCVec(VectorComponentsM)%East  ,'  ' &
                                   ,BCVec(VectorComponentsM)%South ,BCVec(VectorComponentsM)%North ,'  ' &
                                   ,BCVec(VectorComponentsM)%Bottom,BCVec(VectorComponentsM)%Top
    END IF
    WRITE(TermUnit,*)
    WRITE(TermUnit,40) BCP%West,BCP%East,BCP%South,BCP%North,BCP%Bottom,BCP%Top
    WRITE(TermUnit,*)
  END IF
40 FORMAT(' P (boundaryCond_mod)        ',3(2(i2,8x),2x))

END SUBROUTINE InputModelBC

END MODULE BoundaryCond_Mod
