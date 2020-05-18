MODULE DBase_Mod

  USE Shape_Mod, ONLY: MaxCol,ncols,ColName,ColValue &
                      ,ValuesNum,ValuesName,ValuesCol,ValuesIn,ValuesScale

  IMPLICIT NONE

  !INTEGER, PARAMETER :: MaxCol=100 !dh: see Shape_Mod
  !TYPE RecordDbf_T
  CHARACTER :: String(0:MaxCol)*40='' !dh see ReadRecordDbf
  !END   TYPE RecordDbf_T
  CHARACTER*1 :: Version
  INTEGER*1 :: I1Year,Month=4,Day=25
  CHARACTER*1 :: CI1Year,CMonth,CDay
  INTEGER*4 :: Year=2007
  INTEGER :: iRec=0
  INTEGER :: nRecs  ! Number of Rows
  !INTEGER :: ncols  ! Number of Columns !dh: see Shape_Mod
  INTEGER(kind=selected_int_kind(3)) :: lhead  ! =integer*2
  INTEGER(kind=selected_int_kind(3)) :: LenRec  ! =integer*2, Bytes per Row
  !CHARACTER :: ColName(MaxCol)*10 !dh: see Shape_Mod
  CHARACTER :: ColType(MaxCol)*1
  INTEGER :: ColWidth(0:MaxCol) !dh see ReadRecordDbf
  INTEGER :: ColDec(MaxCol)
  INTEGER :: ColOff(MaxCol)
  CHARACTER*1 :: CDummy
  INTEGER :: OutputUnitR=12,OutputUnitW=13


CONTAINS

SUBROUTINE OpenFileDbf(FileName)
  CHARACTER(*) :: FileName
  OPEN(UNIT=OutputUnitR,FILE=TRIM(FileName)//'.dbf',STATUS='UNKNOWN',ACCESS='STREAM')
  OPEN(UNIT=OutputUnitW,STATUS='UNKNOWN')
END SUBROUTINE OpenFileDbf

SUBROUTINE CloseFileDbf
  CLOSE(OutputUnitW)
END SUBROUTINE CloseFileDbf

SUBROUTINE WriteHeaderDbf

  INTEGER :: iCol
  INTEGER :: iOffset
  CHARACTER :: cWidth,cDec
  CHARACTER :: ca*4=''
  CHARACTER :: Flag*1='0'

!  0   	1 byte  	Valid dBASE for Windows table file, bits 0-2 indicate version number: 
!                       3 for dBASE Level 5, 4 for dBASE Level 7.
!                       Bit 3 and bit 7 indicate presence of a dBASE IV or dBASE for Windows memo file; 
!                       bits 4-6 indicate the presence of a dBASE IV SQL table; 
!                       bit 7 indicates the presence of any .DBT memo file (either a dBASE III PLUS type 
!                       or a dBASE IV or dBASE for Windows memo file).
  Version=CHAR(3)
  WRITE(OutputUnitW) Version 
!1-3 	3 bytes 	Date of last update; in YYMMDD format.  Each byte contains the number as a binary.  
!                       YY is added to a base of 1900 decimal to determine the actual year. 
!                       Therefore, YY has possible values from 0x00-0xFF, which allows for a range from 1900-2155.
  I1Year=Year-1900
  WRITE(OutputUnitW) CHAR(I1Year),CHAR(Month),CHAR(Day)
!4-7 	32-bit number 	Number of records in the table. (Least significant byte first.)
  WRITE(OutputUnitW) nRecs
  
!8-9 	16-bit number 	Number of bytes in the header. (Least significant byte first.)
  lhead=32+32*ncols+1
  WRITE(OutputUnitW) lHead
!10-11 	16-bit number 	Number of bytes in the record. (Least significant byte first.)
  LenRec=SUM(ColWidth)+1
  WRITE(OutputUnitW) LenRec
!12-13 	2 bytes 	Reserved; filled with zeros.
!14 	1 byte 	Flag indicating incomplete dBASE IV transaction.
!15 	1 byte 	dBASE IV encryption flag.
!16-27 	12 bytes 	Reserved for multi-user processing.
!28 	1 byte 	Production MDX flag; 0x01 if a production .MDX file exists for this table; 0x00 if no .MDX file exists.
!29 	1 byte 	Language driver ID.
!30-31 	2 bytes 	Reserved; filled with zeros.

  DO iCol = 1,nCols
    cWidth=CHAR(ColWidth(iCol))
    cDec=CHAR(ColDec(iCol))
    ColName(iCol)=ADJUSTR(ColName(iCol))
    WRITE(OutputUnitW,POS=32*iCol+1) ColName(iCol)//CHAR(0),ColType(iCol) &
                                   ,ca,cWidth,cDec
  END DO
  iOffset=32*nCols+35+iRec*LenRec
  WRITE(OutputUnitW,POS=32*nCols+33) CHAR(13)

END SUBROUTINE WriteHeaderDbf

SUBROUTINE ReadHeaderDbf(out_type)

  CHARACTER(10) :: out_type
  INTEGER :: iCol
  INTEGER :: iOffset
  CHARACTER :: cWidth,cDec
  CHARACTER :: ca*4=''
  CHARACTER :: Flag*1='0'

!  0   	1 byte  	Valid dBASE for Windows table file, bits 0-2 indicate version number: 
!                       3 for dBASE Level 5, 4 for dBASE Level 7.
!                       Bit 3 and bit 7 indicate presence of a dBASE IV or dBASE for Windows memo file; 
!                       bits 4-6 indicate the presence of a dBASE IV SQL table; 
!                       bit 7 indicates the presence of any .DBT memo file (either a dBASE III PLUS type 
!                       or a dBASE IV or dBASE for Windows memo file).
  Version=CHAR(3)
  READ(OutputUnitR) Version 
!1-3 	3 bytes 	Date of last update; in YYMMDD format.  Each byte contains the number as a binary.  
!                       YY is added to a base of 1900 decimal to determine the actual year. 
!                       Therefore, YY has possible values from 0x00-0xFF, which allows for a range from 1900-2155.
  I1Year=Year-1900
  READ(OutputUnitR) CI1Year,CMonth,CDay
  I1Year=ICHAR(CI1Year)
  Month=ICHAR(CMonth)
  Day=ICHAR(CDay)
!4-7 	32-bit number 	Number of records in the table. (Least significant byte first.)
  READ(OutputUnitR) nRecs
!8-9 	16-bit number 	Number of bytes in the header. (Least significant byte first.)
  lhead=32+32*ncols+1
  READ(OutputUnitR) lHead
  nCols=(lhead-33)/32
!10-11 	16-bit number 	Number of bytes in the record. (Least significant byte first.)
  LenRec=SUM(ColWidth)+1
  READ(OutputUnitR) LenRec
!12-13 	2 bytes 	Reserved; filled with zeros.
!14 	1 byte 	Flag indicating incomplete dBASE IV transaction.
!15 	1 byte 	dBASE IV encryption flag.
!16-27 	12 bytes 	Reserved for multi-user processing.
!28 	1 byte 	Production MDX flag; 0x01 if a production .MDX file exists for this table; 0x00 if no .MDX file exists.
!29 	1 byte 	Language driver ID.
!30-31 	2 bytes 	Reserved; filled with zeros.
  DO iCol = 1,nCols
!   cWidth=CHAR(ColWidth(iCol))
!   cDec=CHAR(ColDec(iCol))
!   ColName(iCol)=ADJUSTR(ColName(iCol))
    READ(OutputUnitR,POS=32*iCol+1) ColName(iCol),CDummy,ColType(iCol) &
                                   ,ca,cWidth,cDec
    ColWidth(iCol)=ICHAR(cWidth)
    ColDec(iCol)=ICHAR(cDec)
  END DO
  iOffset=32*nCols+35+iRec*LenRec
  READ(OutputUnitR,POS=32*nCols+33) CDummy

! Dummies for temporary write in ReadRecordDbf
  String(0)='   -1'
  ColWidth(0)=5

END SUBROUTINE ReadHeaderDbf


!SUBROUTINE WriteRecord(Record)
!
!  !TYPE(RecordDbf_T) :: Record
!
!  INTEGER :: iCol,iOffset,cW
!  CHARACTER :: Flag*1=''
!
!  iOffset=32*nCols+35+iRec*LenRec
!  WRITE(OutputUnitW,POS=iOffset) Flag
!  iOffset=iOffset+1
!  DO iCol=1,ncols
!    cW=ColWidth(iCol)
!    WRITE(OutputUnitW,Pos=iOffset) String(iCol)(1:cw)
!    iOffset=iOffset+ColWidth(iCol)
!  END DO
!  iRec=iRec+1
!
!END SUBROUTINE WriteRecordDbf

SUBROUTINE ReadRecordDbf(RecordEnd,RecordErr) !dh

  !TYPE(RecordDbf_T) :: Record
  LOGICAL :: RecordEnd,RecordErr
  INTEGER :: iCol,iOffset,cW
  CHARACTER :: Flag*1=''

! Sequence of READ(Dbf)-->WRITE(Num)-->READ(Num)

  RecordEnd=.FALSE.
  RecordErr=.FALSE.
  iOffset=32*nCols+35+iRec*LenRec
  READ(OutputUnitR,END=1,ERR=1,POS=iOffset) Flag
! iOffset=iOffset+1 OSSI
  DO iCol=1,ncols
    cW=ColWidth(iCol)
    READ(OutputUnitR,END=2,ERR=2,Pos=iOffset) String(iCol)(1:cw)
    iOffset=iOffset+ColWidth(iCol)
  END DO

  WRITE(OutputUnitW,*,ERR=2) (String(ValuesCol(iCol))(1:ColWidth(ValuesCol(iCol))),iCol=1,ValuesNum)

  BACKSPACE(OutputUnitW)
  READ(OutputUnitW,*,END=3,ERR=3) (ColValue(iCol),iCol=1,ValuesNum)
  DO iCol=1,ValuesNum
    IF (ValuesCol(iCol)>0) THEN !dh: else input values from ValuesIn (µg/ms)
      ValuesIn(iCol)=ColValue(iCol)*ValuesScale ! g/ma --> µg/ms
    END IF
  END DO

  iRec=iRec+1
  RETURN

1 RecordEnd=.TRUE.
  RecordErr=.TRUE.
  RETURN

2 RecordErr=.TRUE.
  iRec=iRec+1
  RETURN

3 RecordErr=.TRUE.
  iRec=iRec+1
  RETURN

END SUBROUTINE ReadRecordDbf

END MODULE DBase_Mod


