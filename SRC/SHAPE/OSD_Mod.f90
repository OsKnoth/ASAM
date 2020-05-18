MODULE OSD_Mod

  USE Kind_Mod
  IMPLICIT NONE

  CHARACTER*15, PARAMETER :: GebVert='01101510151PGK4' 
  CHARACTER*6, PARAMETER :: NoGebVert1='011104' 
  INTEGER :: InputUnit=10
  INTEGER :: OutUnit=0 !20
  INTEGER :: OutUnitCheck1=0 !21
  INTEGER :: OutUnitCheck2=0 !22
  INTEGER :: OutUnitCheck3=0 !23
  INTEGER :: OutUnitCheck4=0 !24

  TYPE PointOSD_T
    REAL(RealKind) :: xP
    REAL(RealKind) :: yP
    REAL(RealKind) :: zP=0.0d0
    CHARACTER*120 :: ALKDESCode=''
    CHARACTER*100 :: SDLMCode=''
    TYPE(PointOSD_T), POINTER :: Next=>NULL()
  END TYPE PointOSD_T

  TYPE LineOSD_T
    TYPE(PointOSD_T) :: P1,P2
    INTEGER :: no_lp=0
    CHARACTER :: undef='n'
    CHARACTER*120 :: ALKDESCode=''
    CHARACTER*100 :: SDLMCode=''
    CHARACTER*20  :: ENUMCode=''
    TYPE(LineOSD_T), POINTER :: Next=>NULL()
  END TYPE LineOSD_T

  TYPE BuildOSD_T
    TYPE(LineOSD_T), POINTER :: ListOfLines
    INTEGER :: NumberOfLines=0
    TYPE(PointOSD_T), POINTER :: ListOfPoints
    TYPE(PointOSD_T) :: PLeft,PRight
    TYPE (BuildOSD_T), POINTER :: Next=>NULL()
    INTEGER :: anz_no_lp=0
    CHARACTER*7  :: OJ_id=''
    CHARACTER*5  :: STRASSE
    CHARACTER*3  :: HAUSNR
    CHARACTER*20 :: TXT
  END TYPE BuildOSD_T

  TYPE(BuildOSD_T), POINTER :: ListOfBuild
  INTEGER :: NumberOfBuild=0
  INTEGER :: anz_no_bp=0

  INTERFACE OPERATOR (.CROSS.)
    MODULE PROCEDURE CrossProduct
  END INTERFACE


CONTAINS 

FUNCTION CrossProduct(r1,r2)
  TYPE(PointOSD_T) :: CrossProduct
  TYPE(PointOSD_T), INTENT(IN) :: r1,r2

  CrossProduct%xP=r1%yP*r2%zP-r1%zP*r2%yP
  CrossProduct%yP=-r1%xP*r2%zP+r1%zP*r2%xP
  CrossProduct%zP=r1%xP*r2%yP-r1%yP*r2%xP
END FUNCTION CrossProduct

SUBROUTINE WritePoint(Point,Unit)
  TYPE(PointOSD_T) :: Point
  INTEGER :: Unit
  WRITE(Unit,*) Point%xP,Point%yP 
END SUBROUTINE WritePoint

SUBROUTINE WriteLine(nr_zeile,Line,Unit)
  INTEGER :: nr_zeile
  TYPE(LineOSD_T) :: Line
  INTEGER :: Unit
  WRITE(Unit,*) 'Line ',nr_zeile
  CALL WritePoint(Line%P1,Unit)
  CALL WritePoint(Line%P2,Unit)
END SUBROUTINE WriteLine

SUBROUTINE WriteCodeHeader(Unit)
  INTEGER :: Unit
  WRITE(Unit,'(6a80)') 'ALKDESCodeL','SDLMCodeL' &
                      ,'ALKDESCodeP1','SDLMCodeP1' &
                      ,'ALKDESCodeP2','SDLMCodeP2' 
END SUBROUTINE WriteCodeHeader

SUBROUTINE WriteCode(Line,Unit)
  TYPE(LineOSD_T) :: Line
  INTEGER :: Unit
    CHARACTER*80 :: ALKDESCode
    CHARACTER*80 :: SDLMCode
  WRITE(Unit,'(6a80)') Line%ALKDESCode,Line%SDLMCode &
                      ,Line%P1%ALKDESCode,Line%P1%SDLMCode &
                      ,Line%P2%ALKDESCode,Line%P2%SDLMCode 
END SUBROUTINE WriteCode


SUBROUTINE Read_Lpoints(InputFile)

  CHARACTER*80 :: InputFile
  INTEGER :: i

  CHARACTER(150) :: LineFile
  CHARACTER*4 :: NumFile
  CHARACTER*1 :: DummyChar1,DummyChar2
  TYPE(LineOSD_T) :: Line
  TYPE(LineOSD_T), POINTER :: NewLine
  TYPE(LineOSD_T), POINTER :: Current
  TYPE(PointOSD_T) :: Point
  TYPE(PointOSD_T), POINTER :: NewPoint
  TYPE(BuildOSD_T) :: Line_OJ
  TYPE(BuildOSD_T), POINTER :: NewBuild
  TYPE(BuildOSD_T), POINTER :: CurrBuild
  CHARACTER*20 :: puffer_txt

  OPEN(UNIT=InputUnit,FILE=TRIM(InputFile),STATUS='UNKNOWN')
  IF (OutUnitCheck1>0) OPEN(UNIT=OutUnitCheck1,FILE=TRIM(InputFile)//'.check1',STATUS='UNKNOWN')
  IF (OutUnitCheck2>0) OPEN(UNIT=OutUnitCheck2,FILE=TRIM(InputFile)//'.check2',STATUS='UNKNOWN')
  IF (OutUnitCheck3>0) OPEN(UNIT=OutUnitCheck3,FILE=TRIM(InputFile)//'.check3',STATUS='UNKNOWN')
  IF (OutUnitCheck4>0) OPEN(UNIT=OutUnitCheck4,FILE=TRIM(InputFile)//'.check4',STATUS='UNKNOWN')

  IF (OutUnit>0) THEN
    i=0
    DO
      READ(InputUnit,'(a150)',END=1) LineFile
      IF (INDEX(LineFile,'ETYP=OJ')>0) THEN
        EXIT
      END IF
    END DO 
    WRITE(NumFile,'(I4)') i+1000
    OPEN(UNIT=OutUnit,FILE=TRIM(InputFile)//'.'//NumFile,STATUS='UNKNOWN')
    WRITE(OutUnit,'(a150)') LineFile 
    DO
      READ(InputUnit,'(a150)',END=3) LineFile
      IF (INDEX(LineFile,'ETYP=OJ')>0) THEN
        i=i+1
        CLOSE(OutUnit)
        WRITE(NumFile,'(I4)') i+1000
        OPEN(UNIT=OutUnit,FILE=TRIM(InputFile)//'.'//NumFile,STATUS='UNKNOWN')
      END IF
      WRITE(OutUnit,'(a150)') LineFile 
    END DO 
3   CONTINUE
    CLOSE(UNIT=InputUnit)
  END IF 

  OPEN(UNIT=InputUnit,FILE=TRIM(InputFile),STATUS='UNKNOWN')
  SP5:DO
    READ(InputUnit,'(a150)',END=1) LineFile
    IF (INDEX(LineFile,'ETYP=OJ')>0) THEN
      Line_OJ%ListOfLines=>NULL()
      Line_OJ%NumberOfLines=0
      Line_OJ%anz_no_lp=0

      SP4:DO
        READ(InputUnit,'(a150)',END=1) LineFile
        IF (INDEX(LineFile,'OBJID')>0) THEN
           LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
           Line_OJ%OJ_id=''
           Line_OJ%OJ_id=LineFile(:INDEX(LineFile,"'")-1)
           Line_OJ%OJ_id=ADJUSTR(Line_OJ%OJ_id)
           Line_OJ%TXT=''
           Line_OJ%STRASSE=''
           Line_OJ%HAUSNR=''
           IF (OutUnitCheck1>0) WRITE(OutUnitCheck1,'(a7)') Line_OJ%OJ_id
           EXIT SP4 
        END IF
      END DO SP4 

      SP3:DO 
        READ(InputUnit,'(a150)',END=1) LineFile 
        IF (INDEX(LineFile,'ETYP=LI')>0.OR.INDEX(LineFile,'ETYP=BO')>0) THEN
        !IF (INDEX(LineFile,'ETYP=LI')>0) THEN
           LineFile(1:)=LineFile(INDEX(LineFile,"ENUM=")+5:)
           Line%ENUMCode=''
           Line%ENUMCode=LineFile(:INDEX(LineFile," ")-1)
           Line%ENUMCode=(Line%ENUMCode)
           DO 
             Line%no_lp=0
             READ(InputUnit,'(a150)',END=1) LineFile
             IF (INDEX(LineFile,'SDLMCODE')>0) THEN
               LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:) 
               Line%SDLMCode='' 
               Line%SDLMCode=LineFile(:INDEX(LineFile,"'")-1)
               Line%SDLMCode=ADJUSTR(Line%SDLMCode)
             ELSE IF (INDEX(LineFile,'ALKDES')>0) THEN
               LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:) 
               Line%ALKDESCode=''
               Line%ALKDESCode=LineFile(:INDEX(LineFile,"'")-1)
               Line%ALKDESCode=ADJUSTR(Line%ALKDESCode)
             ELSE IF (INDEX(LineFile,'ETYP')>0) THEN
               BACKSPACE InputUnit
               EXIT
             END IF
           END DO
         
           READ(InputUnit,'(a150)',END=1) LineFile
           IF (INDEX(LineFile,'ETYP=PG')>0) THEN
             !SP1:DO
             !  READ(InputUnit,'(a150)',END=1) LineFile
             !  IF (INDEX(LineFile,'ETYP=PG')>0) THEN
                 READ(InputUnit,*,END=1) DummyChar1,DummyChar2,Line%P1%xP 
                 READ(InputUnit,*,END=1) DummyChar1,DummyChar2,Line%P1%yP 
                 DO
                   READ(InputUnit,'(a150)',END=1) LineFile
                   IF (INDEX(LineFile,'SDLMCODE')>0) THEN
                     LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                     Line%P1%SDLMCode='' 
                     Line%P1%SDLMCode=LineFile(:INDEX(LineFile,"'")-1)
                     Line%P1%SDLMCode=ADJUSTR(Line%P1%SDLMCode)
                   ELSE IF (INDEX(LineFile,'ALKDES')>0) THEN
                     LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                     Line%P1%ALKDESCode=''
                     Line%P1%ALKDESCode=LineFile(:INDEX(LineFile,"'")-1)
                     Line%P1%ALKDESCode=ADJUSTR(Line%P1%ALKDESCode)
                   ELSE IF (INDEX(LineFile,'ETYP')>0) THEN
                     BACKSPACE InputUnit
                     !EXIT SP1
                      EXIT
                   END IF
                 END DO
               !END IF
             !END DO SP1
             SP2:DO
               READ(InputUnit,'(a150)',END=1) LineFile
               IF (INDEX(LineFile,'ETYP=PG')>0) THEN
                 READ(InputUnit,*,END=1) DummyChar1,DummyChar2,Line%P2%xP
                 READ(InputUnit,*,END=1) DummyChar1,DummyChar2,Line%P2%yP
                 DO
                   READ(InputUnit,'(a150)',END=2) LineFile
                   IF (INDEX(LineFile,'SDLMCODE')>0) THEN
                     LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                     Line%P2%SDLMCode='' 
                     Line%P2%SDLMCode=LineFile(:INDEX(LineFile,"'")-1)
                     Line%P2%SDLMCode=ADJUSTR(Line%P2%SDLMCode)
                   ELSE IF (INDEX(LineFile,'ALKDES')>0) THEN
                     LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                     Line%P2%ALKDESCode=''
                     Line%P2%ALKDESCode=LineFile(:INDEX(LineFile,"'")-1)
                     Line%P2%ALKDESCode=ADJUSTR(Line%P2%ALKDESCode)
                   ELSE IF (INDEX(LineFile,'ETYP')>0) THEN
                     BACKSPACE InputUnit
                     EXIT SP2
                   END IF
                 END DO
               END IF
             END DO SP2
2            CONTINUE
          !ELSE IF (INDEX(LineFile,'ETYP=LI')>0) THEN
           ELSE IF (INDEX(LineFile,'ETYP=LI')>0.OR.INDEX(LineFile,'ETYP=BO')>0) THEN
             Line%no_lp=1
             Line_OJ%anz_no_lp=Line_OJ%anz_no_lp+1
             BACKSPACE InputUnit
!            IF (OutUnitCheck1>0) WRITE(OutUnitCheck1,'(a120)') Line%ALKDESCode
           ELSE IF (INDEX(LineFile,'ETYP=TX')>0) THEN
             Line%no_lp=1
             Line_OJ%anz_no_lp=Line_OJ%anz_no_lp+1
!            IF (OutUnitCheck1>0) WRITE(OutUnitCheck1,'(a120)') Line%ALKDESCode
             DO
               READ(InputUnit,'(a150)',END=1) LineFile
               IF (INDEX(LineFile,'TXT')>0) THEN
                 LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                 Line_OJ%TXT=''
                 puffer_txt=LineFile(:INDEX(LineFile,"'")-1)
                 puffer_txt=ADJUSTR(puffer_txt)
               ELSE IF (INDEX(LineFile,'STRASSE')>0) THEN
                 LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                 Line_OJ%STRASSE=''
                 Line_OJ%STRASSE=LineFile(:INDEX(LineFile,"'")-1)
                 Line_OJ%STRASSE=ADJUSTR(Line_OJ%STRASSE)
                 Line_OJ%TXT=puffer_txt !!!da Reihenfolge unterschiedlich LI/OJ
               ELSE IF (INDEX(LineFile,'HAUSNR')>0) THEN
                 LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                 Line_OJ%HAUSNR=''
                 Line_OJ%HAUSNR=LineFile(:INDEX(LineFile,"'")-1)
                 Line_OJ%HAUSNR=ADJUSTR(Line_OJ%HAUSNR)
               ELSE IF (INDEX(LineFile,'ETYP')>0) THEN
                 BACKSPACE InputUnit
                 EXIT
               END IF
             END DO
           END IF    ! PG,LI,TXT

          !IF (INDEX(Line%P1%SDLMCode,NoGebVert1)*INDEX(Line%P2%SDLMCode,NoGebVert1)==0) THEN
           IF (INDEX(Line%ALKDESCode,'31K')==0) THEN
             Line_OJ%NumberOfLines=Line_OJ%NumberOfLines+1
             ALLOCATE(NewLine)
             IF (ASSOCIATED(Line_OJ%ListOfLines)) THEN
               Current%Next=>NewLine
               Current=>NewLine
             ELSE  
               Line_OJ%ListOfLines=>NewLine
               Current=>NewLine
             END IF 
             Current=Line
           END IF
 
          !ELSE IF (INDEX(LineFile,'ETYP=LI')>0) THEN
          ! fallen Linien ohne Punkte generell raus, Test wegen Ansicht
          !END IF
 
        ELSE IF (INDEX(LineFile,'ETYP=TX')>0) THEN
             DO
               READ(InputUnit,'(a150)',END=1) LineFile
               IF (INDEX(LineFile,'TXT')>0) THEN
                 LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                 Line_OJ%TXT=''
                 puffer_txt=LineFile(:INDEX(LineFile,"'")-1)
                 puffer_txt=ADJUSTR(puffer_txt)
               ELSE IF (INDEX(LineFile,'SRASSE')>0) THEN
                 LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                 Line_OJ%STRASSE=''
                 Line_OJ%STRASSE=LineFile(:INDEX(LineFile,"'")-1)
                 Line_OJ%STRASSE=ADJUSTR(Line_OJ%STRASSE)
                 Line_OJ%TXT=puffer_txt
               ELSE IF (INDEX(LineFile,'HAUSNR')>0) THEN
                 LineFile(1:)=LineFile(INDEX(LineFile,"'")+1:)
                 Line_OJ%HAUSNR=''
                 Line_OJ%HAUSNR=LineFile(:INDEX(LineFile,"'")-1)
                 Line_OJ%HAUSNR=ADJUSTR(Line_OJ%HAUSNR)
               ELSE IF (INDEX(LineFile,'ETYP')>0) THEN
                 BACKSPACE InputUnit
                 EXIT
               END IF
             END DO
        ELSE IF (INDEX(LineFile,'ETYP=OJ')>0) THEN
          BACKSPACE InputUnit
          EXIT SP3
        END IF  ! LI,TX,OJ
      END DO SP3   

      IF (Line_OJ%anz_no_lp>0) THEN
        anz_no_bp=anz_no_bp+1        
      END IF

     !IF (INDEX(LineFile,'ETYP=')>0) THEN !wenn bestimmte Objekte
        NumberOfBuild=NumberOfBuild+1
        ALLOCATE(NewBuild)
        IF (ASSOCIATED(ListOfBuild)) THEN
          CurrBuild%Next=>NewBuild
          CurrBuild=>NewBuild
        ELSE
          ListOfBuild=>NewBuild
          CurrBuild=>NewBuild
        END IF
        CurrBuild=Line_OJ 
     !END IF   ! INDEX 'OJ'     

    END IF  ! IF(...'ETYP=OJ'....)
  END DO SP5

1 CONTINUE
 
  IF (Line_OJ%anz_no_lp>0) THEN
    anz_no_bp=anz_no_bp+1
  END IF 
  NumberOfBuild=NumberOfBuild+1
  ALLOCATE(NewBuild)
  IF (ASSOCIATED(ListOfBuild)) THEN
    CurrBuild%Next=>NewBuild
    CurrBuild=>NewBuild
  ELSE
    ListOfBuild=>NewBuild
    CurrBuild=>NewBuild
  END IF
  CurrBuild=Line_OJ
  
  !.....................................................
  CALL SPointsEmptyLines
  CALL SortPointCorrectOrder
  CALL InitLeftRightPointBox
  
  CLOSE(UNIT=InputUnit)
  IF (OutUnit>0) CLOSE(Unit=OutUnit)
  IF (OutUnitCheck1>0) CLOSE(OutUnitCheck1)
  IF (OutUnitCheck2>0) Close(OutUnitCheck2)
  IF (OutUnitCheck3>0) Close(OutUnitCheck3)
  IF (OutUnitCheck4>0) CLOSE(OutUnitCheck4)
END SUBROUTINE Read_Lpoints


SUBROUTINE SPointsEmptyLines

  TYPE(LineOSD_T), POINTER :: CurrLine,SearchLine
  TYPE(BuildOSD_T), POINTER :: CurrBuild,SearchBuild
  INTEGER :: pos

  !1. OBJID nur LI ohne Pkt", "2. Search-ID OBJID und Line" 
  CurrBuild=>ListOfBuild
  SP5:DO
     IF (CurrBuild%anz_no_lp==0) THEN
       IF (ASSOCIATED(CurrBuild%Next)) THEN
         CurrBuild=>CurrBuild%Next
       ELSE
         EXIT SP5
       END IF
     ELSE  ! CurrBuild%anz_no_lp>0
       IF (OutUnitCheck2>0) WRITE(OutUnitCheck2,'(a7)') CurrBuild%OJ_id 
       CurrLine=>CurrBuild%ListOfLines
       SP6:DO
          IF (CurrLine%no_lp==0) THEN
            IF (ASSOCIATED(CurrLine%Next)) THEN
              CurrLine=>CurrLine%Next
            ELSE
              EXIT SP6
            END IF
          ELSE    ! CurrLine%no_lp==1
            IF (OutUnitCheck2>0) WRITE(OutUnitCheck2,'(a120)') CurrLine%ALKDESCode
            SearchBuild=>ListOfBuild
            SP7:DO
               SearchLine=>SearchBuild%ListOfLines
               SP8:DO
                 IF (SearchLine%no_lp==0) THEN
                    pos=INDEX(SearchLine%ALKDESCode,CurrLine%ALKDESCode)
                    IF (pos>0) THEN
                      IF (SearchLine%ENUMCode==CurrLine%ENUMCode) THEN 
                         IF (OutUnitCheck2>0) WRITE(OutUnitCheck2,'(a20,a7)') "  ....Parallele.....",SearchBuild%OJ_id
                         IF (OutUnitCheck2>0) WRITE(OutUnitCheck2,'(a120)') SearchLine%ALKDESCode
                         IF (OutUnitCheck3>0) WRITE(OutUnitCheck3,'(a120)') SearchLine%ALKDESCode
                         !Zuweisung 
                         CurrLine%P1=SearchLine%P1
                         CurrLine%P2=SearchLine%P2
                         CurrLine%no_lp=0
                         CurrBuild%anz_no_lp=CurrBuild%anz_no_lp-1
                         EXIT SP7
                      ELSE
                         IF (ASSOCIATED(SearchLine%Next)) THEN
                           SearchLine=>SearchLine%Next
                         ELSE
                           EXIT SP8
                         END IF
                      END IF  !(SearchLine%ENUMCode==CurrLine%ENUMCode)
                    ELSE
                      IF (ASSOCIATED(SearchLine%Next)) THEN
                        SearchLine=>SearchLine%Next
                      ELSE
                        EXIT SP8
                      END IF
                    END IF   !(pos>0)
                 ELSE
                   IF (ASSOCIATED(SearchLine%Next)) THEN
                      SearchLine=>SearchLine%Next
                   ELSE
                      EXIT SP8
                   END IF
                 END IF   !(SearchLine%no_lp==0)
               END DO SP8
                                
               IF (ASSOCIATED(SearchBuild%Next)) THEN
                 SearchBuild=>SearchBuild%Next
               ELSE
                 EXIT SP7
               END IF
            END DO SP7
            IF (CurrLine%no_lp==1) THEN
               CurrLine%undef='y'
               IF (OutUnitCheck2>0) WRITE(OutUnitCheck2,*) "             .......keine Parallele....."
            END IF 
 
            IF (ASSOCIATED(CurrLine%Next)) THEN
              CurrLine=>CurrLine%Next
            ELSE
              EXIT SP6
            END IF
          END IF   !CurrLine%no_lp
       END DO SP6

       IF (ASSOCIATED(CurrBuild%Next)) THEN
         CurrBuild=>CurrBuild%Next
       ELSE
         EXIT SP5
       END IF
     END IF  !CurrBuild%anz_no_lp
  END DO SP5
 
END SUBROUTINE SPointsEmptyLines  


SUBROUTINE SortPointCorrectOrder
  INTEGER :: nr_oj,nr_li 
  TYPE(PointOSD_T) :: P_zw
  TYPE(LineOSD_T), POINTER :: CurrLine,FirstLine
  TYPE(BuildOSD_T), POINTER :: CurrBuild

  CurrBuild=>ListOfBuild
  DO nr_oj=0,NumberOfBuild-1
    IF (ASSOCIATED(CurrBuild%ListOfLines)) THEN
      IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "NumberOBJ=", nr_oj
      IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~"
      IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) CurrBuild%OJ_id
      CurrLine=>CurrBuild%ListOfLines
      FirstLine=>CurrBuild%ListOfLines
      nr_li=0
      DO
        IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "NummerLine=",nr_li
        IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "P1   = ",CurrLine%P1%xP,CurrLine%P1%yP
        IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "P2   = ",CurrLine%P2%xP,CurrLine%P2%yP
        IF (ASSOCIATED(CurrLine%Next)) THEN 
          IF (CurrLine%P1%xP==CurrLine%Next%P1%xP.AND. &
              CurrLine%P1%yP==CurrLine%Next%P1%yP) THEN
            IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "Tausche Point "
            P_zw=CurrLine%P1
            CurrLine%P1=CurrLine%P2
            CurrLine%P2=P_zw
            !Gehe davon aus Code-Zuordnung stimmt          
          ELSE IF (CurrLine%P1%xP==CurrLine%Next%P2%xP.AND. &
                   CurrLine%P1%yP==CurrLine%Next%P2%yP) THEN
            IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "Tausche Point "
            P_zw=CurrLine%P1
            CurrLine%P1=CurrLine%P2
            CurrLine%P2=P_zw
          END IF
          IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "---"

          IF (ASSOCIATED(CurrLine%Next)) THEN
            Currline=>CurrLine%Next
          ELSE
            EXIT
          END IF
        ELSE
          IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "First", FirstLine%P1%xP,FirstLine%P1%yP
          IF (CurrLine%P1%xP==FirstLine%P1%xP.AND. &
              CurrLine%P1%yP==FirstLine%P1%yP) THEN
            IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "Tausche  Point "
            IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "---"

            P_zw=CurrLine%P1
            CurrLine%P1=CurrLine%P2
            CurrLine%P2=P_zw
            EXIT
          ELSE IF (CurrLine%P2%xP==FirstLine%P1%xP.AND. &
              CurrLine%P2%yP==FirstLine%P1%yP) THEN
              IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "---"
              EXIT
          ELSE
            IF (OutUnitCheck4>0) WRITE(OutUnitCheck4,*) "Fehler Points letzte Linie,kann nicht mit First verbunden werden"
            EXIT
          END IF
         
        END IF
        nr_li=nr_li+1
      END DO
    END IF

    IF (ASSOCIATED(CurrBuild%Next)) THEN
      CurrBuild=>CurrBuild%Next
    ELSE
      EXIT
    END IF
  END DO

END SUBROUTINE SortPointCorrectOrder



SUBROUTINE  InitLeftRightPointBox 
  INTEGER :: nr_z
  INTEGER :: i
  REAL(RealKind) :: Dist
  TYPE(LineOSD_T), POINTER :: Current
  TYPE(BuildOSD_T), POINTER :: CurrBuild

  CurrBuild=>ListOfBuild
  DO i=0,NumberOfBuild-1
    IF (ASSOCIATED(CurrBuild%ListOfLines)) THEN
      Current=>CurrBuild%ListOfLines
      CurrBuild%PLeft%xP=Current%P1%xP
      CurrBuild%PLeft%yP=Current%P1%yP
      CurrBuild%PRight%xP=Current%P1%xP
      CurrBuild%PRight%yP=Current%P1%yP
      nr_z=1
      DO
        CurrBuild%PLeft%xP=MIN(CurrBuild%PLeft%xP,Current%P1%xP)
        CurrBuild%PLeft%yP=MIN(CurrBuild%PLeft%yP,Current%P1%yP)
        CurrBuild%PLeft%xP=MIN(CurrBuild%PLeft%xP,Current%P2%xP)
        CurrBuild%PLeft%yP=MIN(CurrBuild%PLeft%yP,Current%P2%yP)
        CurrBuild%PRight%xP=MAX(CurrBuild%PRight%xP,Current%P1%xP)
        CurrBuild%PRight%yP=MAX(CurrBuild%PRight%yP,Current%P1%yP)
        CurrBuild%PRight%xP=MAX(CurrBuild%PRight%xP,Current%P2%xP)
        CurrBuild%PRight%yP=MAX(CurrBuild%PRight%yP,Current%P2%yP)
        IF (OutUnit>0) CALL WriteLine(nr_z,Current,OutUnit)
        IF (ASSOCIATED(Current%Next)) THEN
          Current=>Current%Next
        ELSE
          EXIT
        END IF
        nr_z=nr_z+1
      END DO
      Dist=CurrBuild%PRight%xP-CurrBuild%PLeft%xP
      CurrBuild%PRight%xP=CurrBuild%PRight%xP+0.1d0*Dist
      CurrBuild%PLeft%xP=CurrBuild%PLeft%xP-0.1d0*Dist
      Dist=CurrBuild%PRight%yP-CurrBuild%PLeft%yP
      CurrBuild%PRight%yP=CurrBuild%PRight%yP+0.1d0*Dist
      CurrBuild%PLeft%yP=CurrBuild%PLeft%yP-0.1d0*Dist
    END IF

    IF (ASSOCIATED(CurrBuild%Next)) THEN
      CurrBuild=>CurrBuild%Next
    ELSE
      EXIT
    END IF
  END DO
END SUBROUTINE InitLeftRightPointBox

END MODULE OSD_Mod

