MODULE InputTool_Mod

  USE Kind_Mod
  USE Control_Mod

  IMPLICIT NONE

  INTEGER :: iS1,iS2

CONTAINS

SUBROUTINE OpenFile(FileName)

  CHARACTER(*) :: FileName

  iS1=0
  iS2=0
  OPEN(Unit=InputUnit,FILE=TRIM(FileName),STATUS='UNKNOWN')

END SUBROUTINE OpenFile

SUBROUTINE CloseFile

  iS1=0
  iS2=0
  CLOSE(Unit=InputUnit)

END SUBROUTINE CloseFile

SUBROUTINE ClearFile

  iS1=0
  iS2=0

END SUBROUTINE ClearFile

SUBROUTINE LineFile(Back,Start1,Start2,End &
                   ,Name1,Name2,Name3,Name4,Name5,Name6,Name7 &
                   ,R1,R2,R3,R4,R5,R6,R7,R8,R)

  LOGICAL :: Back
  CHARACTER(*), OPTIONAL :: Start1,Start2
  CHARACTER(*) :: End
  CHARACTER(*), OPTIONAL :: Name1,Name2,Name3,Name4,Name5,Name6,Name7
  REAL(RealKind), OPTIONAL :: R1,R2,R3,R4,R5,R6,R7,R8 
  REAL(RealKind), OPTIONAL :: R(:) 

  CHARACTER(300) :: Line
  INTEGER :: i,is
  INTEGER, PARAMETER :: LenWork=20
  REAL(RealKind) :: Work(LenWork)

  Back=.FALSE.
  IF (PRESENT(Start1)) THEN
    IF (iS1==0) THEN
      S1:DO 
        READ(InputUnit,'(a300)',END=1) Line
        IF (Line(1:1)=='#') THEN 
          CYCLE S1
        END IF
        iS1=INDEX(Line,TRIM(Start1))
        IF (iS1>0) THEN
          EXIT
        END IF 
      END DO S1
    END IF
  ELSE
    is1=1
  END IF
1 CONTINUE
  IF (PRESENT(Start2).AND.iS1>0) THEN
    IF (iS2==0) THEN
      S2:DO 
        READ(InputUnit,'(a300)',IOSTAT=is,END=2) Line
        IF (Line(1:1)=='#') THEN 
          CYCLE S2
        END IF
        iS2=INDEX(Line,TRIM(Start2))
        IF (iS2>0) THEN
          EXIT
        END IF 
      END DO S2
    END IF
  ELSE
    is2=1
  END IF
2 CONTINUE
  IF (iS1*iS2>0) THEN
    E:DO
      READ(InputUnit,'(a300)') Line
      IF (Line(1:1)=='#') THEN
        CYCLE E
      END IF
      IF (INDEX(Line,TRIM(End))>0) THEN
        Back=.TRUE.
      ELSE
        IF (PRESENT(Name1)) THEN
          READ(Line,*) Name1
          Line(INDEX(Line,TRIM(Name1)): &
               INDEX(Line,TRIM(Name1))+LEN(TRIM(Name1))-1)=' '
        END IF
        IF (PRESENT(Name2)) THEN
          READ(Line,*) Name2
          Line(INDEX(Line,TRIM(Name2)): &
               INDEX(Line,TRIM(Name2))+LEN(TRIM(Name2))-1)=' '
        END IF
        IF (PRESENT(Name3)) THEN
          READ(Line,*) Name3
          Line(INDEX(Line,TRIM(Name3)): &
               INDEX(Line,TRIM(Name3))+LEN(TRIM(Name3))-1)=' '
        END IF
        IF (PRESENT(Name4)) THEN
          READ(Line,*) Name4
          Line(INDEX(Line,TRIM(Name4)): &
               INDEX(Line,TRIM(Name4))+LEN(TRIM(Name4))-1)=' '
        END IF
        IF (PRESENT(Name5)) THEN
          READ(Line,*) Name5
          Line(INDEX(Line,TRIM(Name5)): &
               INDEX(Line,TRIM(Name5))+LEN(TRIM(Name5))-1)=' '
        END IF
        IF (PRESENT(Name6)) THEN
          READ(Line,*) Name6
          Line(INDEX(Line,TRIM(Name6)): &
               INDEX(Line,TRIM(Name6))+LEN(TRIM(Name6))-1)=' '
        END IF
        IF (PRESENT(Name7)) THEN
          READ(Line,*) Name7
          Line(INDEX(Line,TRIM(Name7)): &
               INDEX(Line,TRIM(Name7))+LEN(TRIM(Name7))-1)=' '
        END IF
        DO i=1,LenWork
          LINE=ADJUSTL(Line)
          IF (LEN(TRIM(Line))>0) THEN
            IF (Line(1:1)=='#') THEN
              EXIT
            ELSE
              READ(Line,*) Work(i)
              Line(1:INDEX(Line,' ')-1)=' '
            END IF
          ELSE 
            EXIT
          END IF
        END DO
      END IF
      EXIT
    END DO E
    
    IF (PRESENT(R1)) THEN
      R1=Work(1)
    END IF
    IF (PRESENT(R2)) THEN
      R2=Work(2)
    END IF
    IF (PRESENT(R3)) THEN
      R3=Work(3)
    END IF
    IF (PRESENT(R4)) THEN
      R4=Work(4)
    END IF
    IF (PRESENT(R5)) THEN
      R5=Work(5)
    END IF
    IF (PRESENT(R6)) THEN
      R6=Work(6)
    END IF
    IF (PRESENT(R7)) THEN
      R7=Work(7)
    END IF
    IF (PRESENT(R8)) THEN
      R8=Work(8)
    END IF
    IF (PRESENT(R)) THEN
      DO i=1,SIZE(R)
        R(i)=Work(i)
      END DO
    END IF
  ELSE
    Back=.TRUE.
  END IF
      
END SUBROUTINE LineFile

END MODULE InputTool_Mod

