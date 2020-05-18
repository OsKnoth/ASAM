MODULE ActivityAimo_Mod

  USE Control_Mod
  USE Chemie_Mod
  USE Aerosol_Mod
  USE Microphysics_Mod


  REAL(RealKind),PARAMETER :: MM_H2O=18.01528E-3
  REAL(RealKind), PARAMETER :: MolH2O = 1.E0/MM_H2O     ! [Mol/l] Wasser

! Control
  INTEGER :: IowMode,mUni=1,rqMode=2

! ---integer variables
    INTEGER :: NoIonMR, NoA, NoC, NoMainGrpMR, NoSubGrpMR

! ---integer constants
    INTEGER :: IndH2O

 REAL(RealKind), ALLOCATABLE :: b1ki(:,:),        &
&                        b2ki(:,:),        &
&                        b3ki(:,:),        &
&                        bca1(:,:),        &
&                        bca2(:,:),        &
&                        bca3(:,:),        &
&                        cca1(:,:),        &
&                        cca2(:,:),        &
&                        cca3(:,:),        &
&                        rcc(:,:),         &                      
&                        qcc(:,:)

    
! ---real arrays: Data fields of MR systems
    REAL(RealKind), ALLOCATABLE :: bMR(:,:), cMR(:,:),   &   ! LIFAC coefficients (Organics, Ions)
&                           bcatan(:,:), ccatan(:,:)  ! LIFAC coefficients (Anions, Cations)
    REAL(RealKind), ALLOCATABLE :: MMk(:)                    ! Molar mass of Subgroups
    INTEGER, ALLOCATABLE :: IndMRic(:),       &       ! Transfo Cat: Species ==> Xcatan
&                           IndMRia(:),       &       !Transfo Ani: Species ==> Xcatan 
&                           IndBackMRic(:),   &       ! Back-Transfo Cat: Species <== Xcatan
&                           IndBackMRia (:),  &       ! Back-Transfo Ani: Species <== Xcatan
&                           IndBMRk(:),       &       ! Transfo Subgroups: Species ==> bMR, cMR
&                           IndMRk(:),        &
&                           IndMRi(:)                 ! Transfo Ions: Species ==> bMR, cMR

    CHARACTER(20), ALLOCATABLE :: IonName(:), NameCation(:), NameAnion(:)


!==================================================================
!===  UNIFAC-Modul: Approach of Ming/Russell (2002)
!==================================================================
!
!
!-----------------------------------------------------
!---  Additional UNIFAC variables: Ming/Russell (2002)
!-----------------------------------------------------
!
!---  integer
      INTEGER :: NoOrgs, NoIons, NoAddGroups
      INTEGER, ALLOCATABLE :: TypeMols(:),TypeOrg(:)
      INTEGER, ALLOCATABLE :: TypeAddGroups(:)
      INTEGER, ALLOCATABLE :: Group2Group(:,:)
      INTEGER, ALLOCATABLE :: OccurOrg(:,:),SubOccurOrg(:,:)

      INTEGER :: WaterIndUni
      INTEGER, ALLOCATABLE :: UnifacInd(:),UnifacBack(:)
      INTEGER, ALLOCATABLE :: UnifacIonInd(:),UnifacIonBack(:)
      INTEGER, ALLOCATABLE :: UnifacType(:)

!---  INTEGER constants
      INTEGER, PARAMETER:: IndSubWater     = 7          ! H2O subgroup No.
      INTEGER, PARAMETER:: MingIndSubWater = 7          ! H2O subgroup No.        
!
!---  real
      REAL(RealKind), ALLOCATABLE :: rk_Ion(:),qk_Ion(:)
      REAL(RealKind), ALLOCATABLE :: OW_IW(:,:), IW_OW(:,:)
      REAL(RealKind), ALLOCATABLE :: matrixA(:,:)

      REAL(RealKind), ALLOCATABLE :: f_star(:)
      REAL(RealKind), ALLOCATABLE :: f_inf(:)         ! molfrac infinity coefficient

      REAL(RealKind), ALLOCATABLE :: f_starIon(:)
      REAL(RealKind), ALLOCATABLE :: f_Ion(:,:)         ! molfrac infinity coefficient for Ions

      REAL(RealKind), ALLOCATABLE :: gamma_ref(:,:,:)   ! reference activity coefficients
!
!---  character
      CHARACTER(20), ALLOCATABLE :: NameMols(:),NameIon(:), NameOrg(:)
      CHARACTER(20), ALLOCATABLE :: NameAddGroups(:)

!-----------------------------------------------------
!---  Global program constants and variables
!-----------------------------------------------------
!
!---  strings
      CHARACTER(6) :: ProgramName = 'unifac'
    
!---  subtab names, initialized in block.inc
      CHARACTER(80) :: SubName,SubtabName

!---  control data values
      REAL(RealKind), PARAMETER :: ERROR_VALUE = -999.0D0, ERROR_BORDER = 1.D-10

!---  constants
      REAL(RealKind) :: KELVIN = 273.15D0, RGAS = 8.314D0

!---  rk and qk data of sub groups
      REAL(RealKind), ALLOCATABLE :: rkU(:), qk(:), mU(:,:)

      INTEGER, ALLOCATABLE :: SubGroups(:),    &  ! where subgroups start 
   &                          MainGroups(:)       ! main grps belonging to sub grp

!-----------------------------------------------------
!---  descriptors for calculation
      INTEGER :: NoGroups,NoSubGroups,NoMols
!     REAL(RealKind) :: Temperature

      INTEGER, ALLOCATABLE :: Occurrences_Ion(:,:), SubOccurrences_Ion(:,:)
      INTEGER, ALLOCATABLE :: Occurrences(:,:), SubOccurrences(:,:)
      REAL(RealKind), ALLOCATABLE :: MolFractions(:),MolFractions_Ion(:)

!---  arrays for sparse structure
      INTEGER, ALLOCATABLE :: SubcColPtr(:),  SubrRowPtr(:)
      INTEGER, ALLOCATABLE :: SubcRowInd(:),  SubrColInd(:)
      INTEGER, ALLOCATABLE :: SubcSparse(:),  SubrSparse(:)

!-----------------------------------------------------
!---  results
      REAL(RealKind) :: f_Water
      REAL(RealKind) :: Exzess, GExzess, SumRX,SumQX,SumLX   
      REAL(RealKind), ALLOCATABLE, PRIVATE ::  rU(:), q(:), l(:), s(:,:), nue(:),  &
   &   phi(:), xi(:), ln_gamma_c(:), x(:,:), TotalX(:),          &
   &   theta(:,:), Totaltheta(:), Theta2(:), GammaU(:,:),         &
   &   TotalGamma(:), ln_gamma_r(:), Results(:), qw(:),    &
   &   ln_gamma_c_ref(:),ln_gamma_r_ref(:), Results_ref(:)

REAL(RealKind), ALLOCATABLE ::  ln_gamma_r_First(:), ln_gamma_r_Second(:),  &
   &   Gamma_First(:,:),     Gamma_Second(:,:),   &
   &   TotalGamma_First(:), TotalGamma_Second(:)    


!-----------------------------------------------------
!---  file names and handles 
      INTEGER :: nonz_SubC, nonz_SubR
      INTEGER :: CurrentResultFile, DetFile
      INTEGER, PARAMETER :: FileRead = 0, FileWrite = 1, FileAppend=2
      INTEGER, PARAMETER :: CtrlFile = 70, MatrixFile = 71 

      CHARACTER(80) :: CtrlFileName,MatrixFileName

!-----------------------------------------------------
!---  Program control
!-----------------------------------------------------
!
!---  System Loop control
      INTEGER :: LoopSystem, MaxSystem


!---  dimensions
      INTEGER :: nc , na 
      INTEGER :: nac, nacc, naca, nc2, na2
      INTEGER :: zMax=2

!---  INTEGER variable arrays
      INTEGER, ALLOCATABLE :: nzc(:),nza(:)          ! charges
      INTEGER, ALLOCATABLE :: nuec(:,:),nuea(:,:)    ! ion number per salt molecule


!---  REAL(RealKind) parameter arrays
      REAL(RealKind), ALLOCATABLE :: b0(:,:),b1(:,:),a1(:,:),       &
!&                             b2(:,:),a2(:,:),c(:,:),        &
&                             b2(:,:),a2(:,:),               &
&                             psim(:,:,:),ThetaX(:,:),       &
&                             psix(:,:,:),ThetaM(:,:)         

!---  REAL(RealKind) variables 
      REAL(RealKind) :: stri, aw
      REAL(RealKind), ALLOCATABLE :: gc(:), ga(:)
      REAL(RealKind), ALLOCATABLE :: cma(:), cmc(:)

!--------------------------------------------------------
!---  arrays for activity coefficients
      INTEGER, ALLOCATABLE :: PitzInd(:)
      REAL(RealKind), ALLOCATABLE :: aw_Pitz(:,:)        ! water activity
      REAL(RealKind), ALLOCATABLE :: gamma_Pitz(:,:,:)   ! species

!--------------------------------------------------------
!---  CHARACTER  arrays with names
      CHARACTER(20), ALLOCATABLE :: mname(:), xname(:)
      CHARACTER(20), ALLOCATABLE :: mxname(:,:)
!--------------------------------------------------------
!---  arrays for activity coefficients
 
    REAL(RealKind), ALLOCATABLE :: aw_LR(:,:)          ! water activity
    REAL(RealKind), ALLOCATABLE :: gamma_LR(:,:,:)     ! LR molal activity coefficients   
    
!--------------------------------------------------------
!---  CHARACTER  arrays with names
    CHARACTER(20), ALLOCATABLE :: catname(:), anname(:)
   

    INTEGER, ALLOCATABLE :: LRInd(:)

!--------------------------------------------------------

CONTAINS


SUBROUTINE InitActivityAimo
  
!---  Define array for Pitzer and UNIFAC index arrays
  ALLOCATE (UnifacInd(nAqua))
  UnifacInd=0.d0    
  ALLOCATE (LRInd(nAqua))
  LRInd=0.0d0

  CALL InitUnifac
  CALL InitLR  
  CALL InitFinalMR

END SUBROUTINE InitActivityAimo

!
!--------------------------------------------------
!---     Integer-Function  IFIND                ---        
!--------------------------------------------------
   INTEGER FUNCTION ifindLR(nspc, string, namen)
!
!    #### Suchen des Indizes fuer Spezi  ####
!
      IMPLICIT NONE
!
      INTEGER  nspc, indspc, indspcLR
      character(20) :: string, namen(nspc)
!
      ifindLR = -1

      DO indspcLR = 1,nspc
         IF (string.eq.namen(indspcLR))  THEN
            ifindLR = indspcLR
            exit
         END IF
      END DO
!
!--------------------------------------------------
   END FUNCTION ifindLR
! ------------------------------------------------------------------------
   SUBROUTINE InitLR( )
!
! ------------------------------------------------------------------------
! The same system considered as per the A. Zeund et al (2008)
! -------------------------------------------------------------------------

      IMPLICIT NONE

!---  internal variables
      INTEGER :: ia, ic,  jt, indspcLR, na1LR, nc1LR, pos
      CHARACTER(20):: string
      
! -------------------------------------------------------------------------
! ---  Define set of Pitzer coefficients
! -------------------------------------------------------------------------
!
!---  Read Pitzer tables
      WRITE(*,801)
      Call init_actLR

! -------------------------------------------------------------------------
!---  Set index transformation 
! -------------------------------------------------------------------------
!
      na1LR = 0
      nc1LR = 0
!
      LRInd(:) = 0 
      DO jt=1,nAqua
         string = ADJUSTL(SpeciesNameAqua(jt))

!---     cations
         indspcLR = ifindLR(nc, string, catname)
         IF (indspcLR > 0) THEN
            nc1LR = nc1LR + 1
            pos = INDEX(string,' ') - 1
            IF ((string(pos:pos) /= 'p') .OR. (Charge(jt) <= 0))  THEN
               WRITE(*,*) 'InitLR ... Error: Wrong Cation Name =',string
               WRITE(*,*) '                    with given Charge =',Charge(jt)
               STOP       'InitLR ... Error: Wrong Cation Name !'
            END IF
            LRInd(jt) = indspcLR
         END IF

!---     anions
         indspcLR = ifindLR(na, string, anname)
         IF (indspcLR > 0) THEN
            na1LR = na1LR + 1
            pos = INDEX(string,' ') - 1
            IF ((string(pos:pos) /= 'm') .OR. (Charge(jt) >= 0))  THEN
               WRITE(*,*) 'InitLR ... Error: Wrong  Anion Name =',string
               WRITE(*,*) '                    with given Charge =',Charge(jt)
               STOP       'InitLR ... Error: Wrong Anion Name !'
            END IF
            LRInd(jt) = -indspcLR
         END IF
      END DO

! -------------------------------------------------------------------------
!---  Print  system considered for LongRange
! -------------------------------------------------------------------------
!
      WRITE(*,802) nc1LR,na1LR, '  No.','  LRInd','  Charge','     Species Name'
      DO jt=1,nAqua
         IF (LRInd(jt) <= 0)  CYCLE
         WRITE(*,803)  jt, LRInd(jt), Charge(jt), SpeciesNameAqua(jt) 
      END DO
      DO jt=1,nAqua
         IF (LRInd(jt) >= 0)  CYCLE
         WRITE(*,803)  jt, LRInd(jt), Charge(jt), SpeciesNameAqua(jt) 
      END DO
      WRITE(*,804)

! -------------------------------------------------------------------------
801   FORMAT(1x/1x,75('=')/ ' LongRange Initialization:' )
802   FORMAT(1x/1x,75('-')/           &
&            ' Considered Zuend System: Cations =',i3,'     Anions =',i3 // &
&            a5,a9,a9,a18 / 1x,40('-'))
!803   FORMAT(i4,i9,f9.1,7x,a18)
803   FORMAT(i4,i9,i9,7x,a18)
804   FORMAT(1x/1x,75('=')/1x)
!
! -------------------------------------------------------------------------
   END SUBROUTINE InitLR

! =========================================================================
! ===  Subroutines for Pitzer Initialization
! =========================================================================
!
!     Time-stamp:                                  <02/01/06 13:12:07 fm>
! ------------------------------------------------------------------------
   SUBROUTINE init_actLR
! ------------------------------------------------------------------------
!
!    LRMode - flag for the input data set;
!            1:  LIFAC
!            2: AIOMFAC
!
! -------------------------------------------------------------------------


!---  internal variables
      INTEGER :: ia, ic, ios
      INTEGER :: ir0=9, ir1=10, ir2=11, ir3=12

      CHARACTER(10)  ::  set
      CHARACTER(15) ::  Path
!
! -------------------------------------------------------------------------
! ---  Define set of Pitzer coefficients
! -------------------------------------------------------------------------
!
      write (6,'(/a/)') 'use General input data...'
      path = 'Data_LR/'
      set  = '.Long'     
!
!------------------------------------------------------------------
! --- Read Pitzer System: Dimensions, Names, Ions, Charges
!------------------------------------------------------------------
!
!      OPEN (ir0,FILE=TRIM(path)//'name_AIOMFAC'//'.dat',STATUS='old',iostat=ios)
       OPEN (ir0,FILE=TRIM(path)//'name'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACTLR...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'name'//'.dat'
         STOP  ' INIT_ACTLR...Error: Check LongRange =Name= File !!'
      END IF

!---  read and set dimensions
      READ(ir0,*)
      READ(ir0,*)  nc, na

      
!---  allocate arrays
      ALLOCATE ( nzc(nc),nza(na),nuec(nc,na),nuea(nc,na) )
      nza(:) = 0
      nzc(:) = 0
      nuea(:,:) = 0
      nuec(:,:) = 0
      ALLOCATE ( cmc(nc),cma(na) )                     ! Same will replace in LongRange interaction term also
      cmc(:) = 0.d0
      cma(:) = 0.d0
      

      ALLOCATE ( catname(nc),anname(na) )
!RW      ALLOCATE ( mxname(nc,na) )

!---  read names
      REWIND (ir0)
      CALL read_nameLR(ir0)
      CLOSE (ir0)
!
!------------------------------------------------------------------
! --- define maximum charge
      DO ic=1,nc
         zMax = MAX(zMax,nzc(ic))
      END DO

      DO ia=1,na
         zMax = MAX(zMax,ABS(nza(ia)))
      END DO

! --- determine the stoichiometric coefficients
      DO ic=1,nc
         DO ia=1,na
            nuec(ic,ia) = ABS(nza(ia))
            nuea(ic,ia) = nzc(ic)
         END DO
      END DO

!
! -------------------------------------------------------------------------
   END SUBROUTINE init_actLR

! =========================================================================

   SUBROUTINE read_nameLR (iread)
! ----------------------------------------------------------------------

!---  external variables
      INTEGER :: iread
!
!---  internal variables
      INTEGER :: ia, ic, nc1LR, na1LR
! ----------------------------------------------------------------------

      REWIND (iread)
      READ (iread,*)
      READ (iread,*) nc1LR,na1LR
      IF (nc1LR /= nc .OR. na1LR /= na) THEN
        print*,'dimensions of the parameter space does not coincident'  &
&              //' with the model dimenssions!'
        print*,'nc1LR: ',nc1LR,'  nc: ',nc
        print*,'na1LR: ',na1LR,'  na: ',na
        stop 'READ_nmae...'
      END IF

      READ (iread,*)
      READ (iread,*) (catname(ic),ic=1,nc)
      READ (iread,*) (nzc(ic),ic=1,nc)
      READ (iread,*)
      READ (iread,*) (anname(ia),ia=1,na)
      READ (iread,*) (nza(ia),ia=1,na)
!RW      READ (iread,*)
!RW      DO ic=1,nc
!RW        READ (iread,*) (mxname(ic,ia),ia=1,na)
!RW      END DO
!
      CLOSE (iread)

   END SUBROUTINE read_nameLR

!==============================================================
   SUBROUTINE UpperCase(String)
!
!-------------------------------------------
!---  Uppercases string
!-------------------------------------------
     INTEGER :: i
     CHARACTER*(*) :: String

     DO i=1 ,LEN_TRIM(String)
        IF (String(i:i) >= 'a' .AND. String(i:i) <= 'z')     &
  &             String(i:i) = CHAR (ICHAR(String(i:i)) - 32)
     END DO          
   END SUBROUTINE

!================================================================
!
!--------------------------------------------------
!---     Integer-Function  IFIND                ---        
!--------------------------------------------------
   INTEGER FUNCTION ifind(nspc, string, namen)
!
!    #### Suchen des Indizes fuer Spezi  ####
!
      IMPLICIT NONE
!
      INTEGER  nspc, indspc
      character(20) :: string, namen(nspc)
!
      ifind = -1

      DO indspc = 1,nspc
         IF (string.eq.namen(indspc))  THEN
            ifind = indspc
            exit
         END IF
      END DO
!
!--------------------------------------------------
   END FUNCTION ifind
!
!==============================================================
   LOGICAL FUNCTION IsAlphaNumeric(Char)
!
!-------------------------------------------
!       Tests for allphanumeric characters (0..9, a-z, A-Z).
!-------------------------------------------
        CHARACTER*(*) :: Char
!
     IsAlphaNumeric = ((char >= '0') .AND. (char <= '9'))       &
  &                 .OR.   ((char >= 'a') .AND. (char <= 'z'))  &
  &                 .OR.   ((char >= 'A') .AND. (char <= 'Z'))
   END FUNCTION IsAlphaNumeric
!
!==============================================================
   SUBROUTINE GetNextANWord(Ptr,String,Word)
!
!-------------------------------------------
!       Writes the next ALPHANUMERIC word in the string to the variable 
!       word and moves the pointer further.
!-------------------------------------------
!
     CHARACTER*(*) :: String,Word
     INTEGER :: i,Length,Ptr
!
     Length = LEN_TRIM(String)
!    --- search for start
     DO WHILE (.NOT. IsAlphaNumeric(String(Ptr:Ptr)) .AND. (Ptr <= Length))
        Ptr = Ptr + 1
     END DO

!    --- write word
     Word = ' '
     i    = 1
     DO WHILE (IsAlphaNumeric(String(Ptr:Ptr)) .AND. (Ptr <= Length))
        Word(i:i) = String(Ptr:Ptr)
        i   = i   + 1
        Ptr = Ptr + 1
     END DO          
   END SUBROUTINE  GetNextANWord
!
!==============================================================
   SUBROUTINE GetNextWord(Ptr,String,Word)
!
!-------------------------------------------
!       Writes the next word enclosed in blanks in the string to
!       the variable word and moves the pointer further.
!-------------------------------------------
!
     CHARACTER*(*) :: String,Word
     INTEGER :: i,Length,Ptr
!
     Length = LEN_TRIM(String)
!    --- search for start
     DO WHILE ((String(Ptr:Ptr) < '!' ) .AND. (Ptr <= Length))
        Ptr = Ptr + 1
     END DO

!    --- write word
     Word = ' '
     i    = 1
     DO WHILE ((String(Ptr:Ptr) >= '!') .AND. (Ptr <= Length))
        Word(i:i) = String(Ptr:Ptr)
        i   = i   + 1
        Ptr = Ptr + 1
     END DO          
   END SUBROUTINE  GetNextWord
!
!==============================================================
   SUBROUTINE PromptString(Prompt,Value)
!
!-------------------------------------------
!       User Interface routine
!       Prompts for a String. The first character has a special
!       meaning:
!       !               Execution of a shell command
!       Subroutine to all input routines.
!-------------------------------------------
!
     CHARACTER*(*) :: Value,Prompt
     INTEGER :: ParNo
     LOGICAL :: CommandEntered
!
     CommandEntered= .TRUE.
     DO WHILE (CommandEntered)
        Value = ' '
        WRITE (*,'(1X,A,$)') Prompt(1:LEN_TRIM(Prompt))
        READ  (*,'(A)',ERR=99,END=99)    Value
99      CONTINUE
!       -- check for special codes (UNIX only...)
        IF (Value(1:1) == '!') THEN
           CALL system(Value(2:))
           CommandEntered=.TRUE.
        ELSE
           CommandEntered=.FALSE.
        END IF
     END DO ! InputLoop
   END SUBROUTINE PromptString
!
!==============================================================
   SUBROUTINE InputString(Prompt,Value)
!
!-------------------------------------------
!       User interface routine
!       Prompts for String.
!-------------------------------------------
!
     ChARACTER*(*) :: Prompt,Value
     CHARACTER(80) :: FullPrompt
!
     FullPrompt = Prompt // '?:'
     CALL PromptString(FullPrompt,Value)
   END SUBROUTINE InputString
!
!==============================================================
   LOGICAL FUNCTION InputLogical(Prompt,Default)
!
!-------------------------------------------
!       User interface routine
!       Prompts for Boolean value. Returns default if an
!       empty string is entered.
!-------------------------------------------
!
     ChARACTER*(*) :: Prompt
     ChARACTER(80) :: FullPrompt,Input,HString
     LOGICAL :: Ok,Result,Default
!
     OK= .FALSE.
     HString = '? (Y/N) <'
     FullPrompt=Prompt
     IF (Default) THEN
        HString(10:) = 'Y>:'
     ELSE
        HString(10:) = 'N>:'
     END IF 
     FullPrompt(LEN_TRIM(FullPrompt)+1:) = HString
     DO WHILE (.NOT. OK)
        CALL PromptString(FullPrompt(1:LEN_TRIM(FullPrompt)),Input)
        IF (INDEX('yYjJ',Input(1:1)) /= 0) THEN
           Result= .True.            ! YES
           OK    = .True.
           ELSE
              IF (INDEX('nN',Input(1:1)) /= 0) THEN
                 Result= .False.        ! NO
                 OK    = .True.
              ELSE
                 IF (Input == ' ') THEN
                    Result= Default       ! Default
                    OK    = .True.
              END IF
           END IF     
        END IF
        IF (.NOT. OK) WRITE (*,'(1X,A)')  'Error, try again!'
     END DO ! InputLoop
     InputLogical = Result
   END FUNCTION InputLogical
!
!==============================================================
   INTEGER FUNCTION InputInteger(Prompt,Min,Max,Default)
!
!-------------------------------------------
!       User interface routine
!       Prompts for Integer until valid value entered (within min, max)
!       Returns default if empty string is entered.
!-------------------------------------------
    
     CHARACTER*(*) :: Prompt
     CHARACTER(80) :: FullPrompt,Input
     INTEGER :: Result,Min,Max,Default
     LOGICAL :: OutOfRange
!
     WRITE (FullPrompt,'(A,A,I5,A,I5,A,I5,A)')    &
  &    Prompt(1:LEN_TRIM(Prompt)),'? (',  Min,'..',Max,') <',Default,'> :'
     OutOfRange = .TRUE.   
     DO WHILE (OutOfRange)
        CALL PromptString(FullPrompt,Input)
        IF (LEN_TRIM(Input) == 0) THEN
           Result = Default
        ELSE
           READ (Input,*,ERR=1000) Result
        END IF
!           Check Range
        IF ((Result >= Min) .AND. (Result <= Max)) THEN
           OutOfRange = .FALSE.
           CYCLE          
        END IF
1000    OutOfRange = .TRUE.
        WRITE (*,'(1X,A)')  'Error: Number out of range or mispelt, try again!'
     END DO ! InputLoop
     InputInteger=Result
   END FUNCTION InputInteger
!
!==============================================================
   REAL(RealKind) FUNCTION InputReal(Prompt,Min,Max,Default)
!
!-------------------------------------------
!       User interface routine
!       Prompts for real until valid value entered (within min, max)
!       Returns default if empty string is entered.
!-------------------------------------------
!
     CHARACTER*(*) :: Prompt
     CHARACTER(80) :: FullPrompt,Input
     REAL(RealKind) :: Result,Min,Max,Default
     LOGICAL :: OutOfRange
!
     WRITE (FullPrompt,'(A,A,F12.5,A)')  Prompt(1:LEN_TRIM(Prompt)), ' <',Default,'> :'
     OutOfRange = .TRUE.   
     DO WHILE (OutOfRange)
        CALL PromptString(FullPrompt,Input)
        IF (LEN_TRIM(Input) == 0) THEN
           Result = Default
        ELSE
           READ (Input,*,ERR=1000) Result
        END IF
!            --- Check Range
        IF ((Result >= Min) .AND. (Result <= Max)) THEN
           OutOfRange = .FALSE.
           CYCLE         
        END IF
1000    OutOfRange = .TRUE.
        WRITE (*,'(1X,A)')  'Error: Number out of range or mispelt, try again!'
     END DO ! InputLoop
     InputReal=Result
   END FUNCTION InputReal
!
!==============================================================
   LOGICAL FUNCTION DoesFileExist(FileName)
!
!-------------------------------------------
!       File handling routine   
!       Returns True if the file given exists.
!-------------------------------------------
!         
     LOGICAL :: Exists
     CHARACTER*(*) :: FileName
!         
     INQUIRE (File=FileName,Exist=exists,Err=1000)
     GOTO 2000
1000 Exists=.FALSE.
2000 Continue
     DoesFileExist=Exists
   END FUNCTION DoesFileExist
!
!==============================================================
   LOGICAL FUNCTION FileOpen(UnitNo,FileName,FileMode)
!
!-------------------------------------------
!       File handling routine
!       Opens a file for sequential in/output as specified by mode 
!       (0 = READ, 1 = WRITE, 2 = Append)
!       Returns True on success.
!-------------------------------------------
!         
     INTEGER :: UnitNo,IoStat,FileMode
     LOGICAL :: ExistedBefore
     CHARACTER*(*)  :: FileName
     CHARACTER(10)  :: OpenStatus,OpenMode
     CHARACTER(500) :: Buffer
!         
! 
     IF (FileMode == 0) THEN      ! read 
        OpenStatus = 'OLD'
        OpenMode   = 'read'
     ELSE
        IF (FileMode == 1) THEN   ! write 
           OpenStatus='UNKNOWN'
           OpenMode   = 'write'
        ELSE ! append
           ExistedBefore = DoesFileExist(FileName)
           IF (ExistedBefore) THEN
              OpenStatus = 'OLD'
           ELSE
              OpenStatus='UNKNOWN'
           END IF
           OpenMode   = 'write'               
        END IF
     END IF

!     -- open!
     Open(UnitNo,FILE=FileName,ACCESS='SEQUENTIAL',   &
!#ifdef _NDP_
!  &          MODE=OpenMode,SHARE='DENYNONE',
!#endif
  &          STATUS=OpenStatus,ERR=1000,IOSTAT=IoStat)
     FileOpen = .TRUE.
     GOTO 1010
1000 FileOpen = .FALSE.
1010 CONTINUE

!    - if append: Forward file to end
     IF (FileMode == 2 .AND. ExistedBefore .AND. FileOpen) THEN
!       -- read till eof 
        DO WHILE  (.TRUE.)
           READ (UnitNo,'(A)',End=2000,Err=2000) Buffer
        END DO
2000    CONTINUE        
     END IF ! append

   END FUNCTION FileOpen
!
!==============================================================
   SUBROUTINE ErrorAndExit(ProgramName,ErrorMessage)
!
!-------------------------------------------
!       Displays error and exits program
!-------------------------------------------
!         
     CHARACTER*(*) :: ErrorMessage
     CHARACTER*(*) :: ProgramName
!         
     WRITE (*,'(/,1X,A,A,A,/,/,1X,A,A,/,/)')  ProgramName,': ', ErrorMessage,ProgramName,' quitting.'
     CALL EXIT(1)
   END SUBROUTINE 
!
!==============================================================
   LOGICAL FUNCTION StringToInteger(String,Number)
!         
!-------------------------------------------
!       Converts String into integer. Returns true on success
!-------------------------------------------
!         
     CHARACTER*(*) :: String
     INTEGER :: Number 
!         
     READ (String,'(I6)',Err=111) Number
     GO TO 112
111  StringToInteger = .FALSE.
     Number = 0
112  StringToInteger = .TRUE.        
   END FUNCTION StringToInteger

    Subroutine InitUnifac ( )
!
!================================================================
!===   Initialization  Subroutines
!================================================================

       INTEGER :: iCell, jt


!----------------------------------------------------------------
!--- open input files
      IF (rqmode ==1) THEN
!        
        MatrixFileName = 'Data_Unifac/UnifacFinal.uni'
        CtrlFileName   = 'Data_Unifac/Ming.ctr'
      ELSE IF (rqmode ==2) THEN
!       MatrixFileName = 'Data_Zuend/Zuend.uni'
!       CtrlFileName   = 'Data_Zuend/Zuend.ctr'
!       MatrixFileName = 'Data_Unifac/UnifacFinal.uni'
!       CtrlFileName   = 'Data_Unifac/MRJH_Full.ctr'
        MatrixFileName = 'Data_Unifac/UnifacFinal.uni'
        CtrlFileName   = 'Data_Unifac/NoRadical.ctr'
      END IF
       IF (.NOT. FileOpen(MatrixFile,MatrixFileName,FileRead)) THEN
          CALL ErrorAndExit(ProgramName,'Cannot open matrix file!')
       END IF
       IF (.NOT. FileOpen(CtrlFile,CtrlFileName,FileRead)) THEN
          CALL ErrorAndExit(ProgramName,'Cannot open matrix file!')
       END IF

!--- initialize group-dependend values
       CALL ReadParameters

!--- initialize chemical system
       IF (.NOT. ReadSystem( )) THEN
          CALL ErrorAndExit(ProgramName,'Cannot read CTRL File!')
       END IF
       CALL InitSystem ( )

!--- organize exchange between reaction scheme and UNIFAC
       CALL OrderSystem ( )

!----------------------------------------------------------------
!--- define rk and qk for ions (mo_UNIFAC: IndSubWater=6)
       IF (rqMode == 1)  THEN
          rkU(NoSubGroups-NoIons:NoSubGroups) = rkU(MingIndSubWater)
          qk(NoSubGroups-NoIons:NoSubGroups) = qk(MingIndSubWater)
       ELSE IF (rqMode == 2)  THEN
          rkU(NoSubGroups-NoIons:NoSubGroups) = rkU(NoSubGroups-NoIons:NoSubGroups)
          qk(NoSubGroups-NoIons:NoSubGroups) = qk(NoSubGroups-NoIons:NoSubGroups)
       END IF


    END Subroutine InitUnifac

!==========================================================================
!
   LOGICAL FUNCTION ReadSystem( )
!
! ------------------------------------------
! ---  Reads CTRL File
! ------------------------------------------

!---  internal variables
      INTEGER :: i, icol, ncol, maxcol, irow
      INTEGER :: Ptr,smd,EoL,Pos
      INTEGER :: No, ExtNo, Occ, MainGroup, SubGroup

      REAL(RealKind) ::  Sum

      CHARACTER(80)  :: Word
      CHARACTER(256) :: Buffer

      LOGICAL :: EndOfFile


!----------------------------------------------------------------
      EndOfFile  = .FALSE.
      ReadSystem = .FALSE.

      NoOrgs = 0
      NoIons = 0

!----------------------------------------------------------------
!--- read file 
      DO WHILE (.NOT. EndOfFile)
        READ (CtrlFile,'(A)',End=99,Err=99) Buffer

        !-- check for comment line
        IF (Buffer(1:1) == '#' .OR. LEN_TRIM(Buffer) == 0)  CYCLE

!--- eval header or smd file
        Ptr = 1
        CALL GetNextWord(Ptr,Buffer,Word)
        CALL UpCase(Word)

!----------------------------------------------------------------
!--- initialize organics
        IF (Word == 'ORGANICS') THEN

           !-- read size
           CALL GetNextWord(Ptr,Buffer,Word)
           READ (Word,*) NoOrgs

           IF (NoOrgs < 1)  CYCLE
           ALLOCATE (NameOrg(NoOrgs))
           ALLOCATE (TypeOrg(NoOrgs))
           ALLOCATE (OccurOrg(NoOrgs,NoGroups))
           ALLOCATE (SubOccurOrg(NoOrgs,NoSubGroups))
 
           NameOrg(:) = ' '
           TypeOrg(:) = 0

           OccurOrg(:,:)    = 0
           SubOccurOrg(:,:) = 0

           !-- read organics line
           No = 1
           DO WHILE (No <= NoOrgs)
              READ (CtrlFile,'(A)',End=99,Err=99) Buffer
              IF (Buffer(1:1) == '#' .OR. LEN_TRIM(Buffer) == 0)  CYCLE
              EoL = INDEX(Buffer,'#') - 1
              IF (EoL <= 0)  EoL = LEN_TRIM(Buffer) - 1

              Ptr = 1
              CALL GetNextWord(Ptr,Buffer,Word)         ! No.
              CALL GetNextWord(Ptr,Buffer,Word)         ! Short Name
              CALL UpCase(Word)
              READ (Word,*) NameOrg(No)
              CALL GetNextWord(Ptr,Buffer,Word)         ! Long Name
              CALL GetNextWord(Ptr,Buffer,Word)         ! Type
              READ (Word,*) TypeOrg(No)

              !-- read group information
              DO WHILE (Ptr < EoL)
                 CALL GetNextWord(Ptr,Buffer,Word)      ! Group infos
                 IF ( Word(1:1) /= '(' )  CYCLE
                 Pos = INDEX(Word,')') - 1
                 READ (Word(2:Pos),*)  ExtNo,Occ

                 MainGroup = (ExtNo / 10)
                 SubGroup  = ExtNo - MainGroup * 10

                 OccurOrg(No,MainGroup) = OccurOrg(No,MainGroup) + Occ
                 iCol = SubGroups(MainGroup) + SubGroup - 1
                 SubOccurOrg(No,iCol) = SubOccurOrg(No,iCol) + Occ
              END DO
              
              No = No + 1
           END DO

!----------------------------------------------------------------
!--- initialize ions
        ELSE IF (Word == 'IONS' .AND. IowMode >= 1) THEN

           !-- read size
           CALL GetNextWord(Ptr,Buffer,Word)
           READ (Word,*) NoIons

           IF (NoIons < 1)  CYCLE
           ALLOCATE (NameIon(NoIons))
           ALLOCATE (rk_Ion(NoIons))
           ALLOCATE (qk_Ion(NoIons))

           !-- read organics line
           No = 1
           DO WHILE (No <= NoIons)
              READ (CtrlFile,'(A)',End=99,Err=99) Buffer
              IF (Buffer(1:1) == '#' .OR. LEN_TRIM(Buffer) == 0)  CYCLE
              EoL = INDEX(Buffer,'#')
              IF (EoL <= 0)  EoL = LEN_TRIM(Buffer)

              Ptr = 1
              CALL GetNextWord(Ptr,Buffer,Word)         ! No.
              CALL GetNextWord(Ptr,Buffer,Word)         ! Short Name
              CALL UpCase(Word)
              READ (Word,*) NameIon(No)
              CALL GetNextWord(Ptr,Buffer,Word)         ! Long Name
              CALL GetNextWord(Ptr,Buffer,Word)         ! RK
              READ (Word,*) rk_Ion(No)
              CALL GetNextWord(Ptr,Buffer,Word)         ! QK
              READ (Word,*) qk_Ion(No)
              
              No = No + 1
           END DO

!----------------------------------------------------------------
!--- initialize additional organig functional groups
        ELSE IF (Word == 'GROUPS' .AND. IowMode >= 1) THEN

           !-- read size
           CALL GetNextWord(Ptr,Buffer,Word)
           READ (Word,*) NoAddGroups

           IF (NoAddGroups < 1)  CYCLE
           ALLOCATE (NameAddGroups(NoAddGroups))
           ALLOCATE (TypeAddGroups(NoAddGroups))
           ALLOCATE (Group2Group(NoAddGroups,2))

           !-- read additional groups 
           No = 1
           DO WHILE (No <= NoAddGroups)
              READ (CtrlFile,'(A)',End=99,Err=99) Buffer
              IF (Buffer(1:1) == '#' .OR. LEN_TRIM(Buffer) == 0)  CYCLE
              EoL = INDEX(Buffer,'#')
              IF (EoL <= 0)  EoL = LEN_TRIM(Buffer)

              Ptr = 1
              CALL GetNextWord(Ptr,Buffer,Word)         ! No.
              CALL GetNextWord(Ptr,Buffer,Word)         ! Short Name
              CALL UpCase(Word)
              READ (Word,*) NameAddGroups(No)
              CALL GetNextWord(Ptr,Buffer,Word)         ! TypeAddGroups
              READ (Word,*) TypeAddGroups(No)
              CALL GetNextWord(Ptr,Buffer,Word)         ! Group2Group
              READ (Word,*) Group2Group(No,1)           
              
              No = No + 1
           END DO

!----------------------------------------------------------------
!--- read interaction parameters IW-OW 
        ELSE IF (Word == 'IW-OW' .AND. IowMode >= 1) THEN

           !-- read size
           CALL GetNextWord(Ptr,Buffer,Word)
           READ (Word,*) MaxCol

           IF (MaxCol < 1)  CYCLE
           ALLOCATE (IW_OW(NoIons,NoAddGroups))

           !-- read additional groups 
           iCol = 0
           DO WHILE (iCol < NoAddGroups)
              nCol = MIN(MaxCol,NoAddGroups-iCol)
              iRow = 1

              DO WHILE (iRow <= NoIons)
                 READ (CtrlFile,'(A)',End=99,Err=99) Buffer
                 IF (Buffer(1:1) == '#' .OR. LEN_TRIM(Buffer) == 0)  CYCLE

                 Ptr = 1
                 CALL GetNextWord(Ptr,Buffer,Word)         ! Name (comments)s
                 Ptr = Ptr + 1
                 READ (Buffer(Ptr:),*) (IW_OW(iRow,iCol+i),i=1,nCol) 
              
                 iRow = iRow + 1
              END DO
              
              iCol = iCol + nCol
           END DO
           
           !-- LIFAC mode: Set IW-OW interactions to zero
           IF (mUni == 2) THEN
              IW_OW(:,:) = 0.E0
           END IF

!----------------------------------------------------------------
!--- read interaction parameters IW-OW 
        ELSE IF (Word == 'OW-IW' .AND. IowMode >= 1) THEN

           !-- read size
           CALL GetNextWord(Ptr,Buffer,Word)
           READ (Word,*) MaxCol

           IF (MaxCol < 1)  CYCLE
           ALLOCATE (OW_IW(NoAddGroups,NoIons))

           !-- read additional groups 
           iCol = 0
           DO WHILE (iCol < NoIons)
              nCol = MIN(MaxCol,NoIons-iCol)
              iRow = 1

              DO WHILE (iRow <= NoAddGroups)
                 READ (CtrlFile,'(A)',End=99,Err=99) Buffer
                 IF (Buffer(1:1) == '#' .OR. LEN_TRIM(Buffer) == 0)  CYCLE

                 Ptr = 1
                 CALL GetNextWord(Ptr,Buffer,Word)         ! Name (comments)s
                 Ptr = Ptr + 1
                 READ (Buffer(Ptr:),*) (OW_IW(iRow,iCol+i),i=1,nCol) 
              
                 iRow = iRow + 1
              END DO
              
              iCol = iCol + nCol
           END DO

           
           !-- LIFAC mode: Set IW-OW interactions to zero
           IF (mUni == 2) THEN
              OW_IW(:,:) = 0.E0
           END IF
        END IF  

        CYCLE
99      EndOfFile = .TRUE.
      END DO   ! 

      ReadSystem = .TRUE.
!
   END FUNCTION ReadSystem

!==========================================================================

   SUBROUTINE InitSystem( )
!
! ------------------------------------------
! ---  Initialize Chemical System for UNIFAC
! ------------------------------------------

      INTEGER :: iag, ia, ie, ns, NoGrp, NoSub
      INTEGER :: ig, igNew, igOld, No, NoNew, NoOld
      INTEGER :: NoGroupsOld, NoSubGOld, MainGrp
      INTEGER :: inonz

      INTEGER :: SubGroups0(NoGroups+1), MainGroups0(NoSubGroups)
      REAL(RealKind) :: rk0(NoSubGroups), qk0(NoSubGroups)

      CHARACTER(80)  :: Word
      CHARACTER(256) :: Buffer


!------------------------------------------
!--- Building of extended group arrays
      NoGroupsOld = NoGroups 
      NoSubGOld   = NoSubGroups 

      IF (NoIons > 0 .OR. NoAddGroups > 0) THEN

      !-- save old group arrays, rk and qk
         SubGroups0(:)  = SubGroups(:)
         MainGroups0(:) = MainGroups(:)

         rk0(:) = rkU(:)
         qk0(:) = qk(:)

         DEALLOCATE (SubGroups, MainGroups)
         DEALLOCATE (rkU, qk)

      !-- determine number of groups of the extended system
         NoGroups    = NoGroups + NoIons
         NoSubGroups = NoSubGroups + NoIons
         DO No=1,NoAddGroups
            IF (TypeAddGroups(No) <= 0) CYCLE
            MainGrp     = Group2Group(No,1)
            NoGroups    = NoGroups + 1
            NoSubGroups = NoSubGroups + SubGroups0(MainGrp+1)-SubGroups0(MainGrp)
         END DO
         Group2Group(:,2) = Group2Group(:,1)          

      !-- allocation of new group arrays rk and qk
         ALLOCATE (SubGroups(NoGroups+1), MainGroups(NoSubGroups))
         ALLOCATE (rkU(NoSubGroups), qk(NoSubGroups))

      !-- restore group arrays, rk and qk
         NoGrp = NoGroupsOld 
         NoSub = NoSubGOld

         SubGroups(1:NoGrp+1) = SubGroups0(:)
         MainGroups(1:NoSub)  = MainGroups0(:)

         rkU(1:NoSub) = rk0(:)
         qk(1:NoSub) = qk0(:)

      !-- add additional organics groups 
         DO No=1,NoAddGroups
            IF (TypeAddGroups(No) <= 0) CYCLE
            MainGrp = Group2Group(No,1)
            NoGrp   = NoGrp + 1

            ia = SubGroups0(MainGrp)
            ie = SubGroups0(MainGrp+1) - 1
            ns = SubGroups0(MainGrp+1)-SubGroups0(MainGrp)
            MainGroups(NoSub+1:NoSub+ns) = NoGrp
           
            rkU(NoSub+1:NoSub+ns) = rkU(ia:ie)
            qk(NoSub+1:NoSub+ns) = qk(ia:ie)

            NoSub = NoSub + ns

            SubGroups(NoGrp+1) = NoSub + 1
            Group2Group(No,2)  = NoGrp 
         END DO

      !-- add ions
         DO No=1,NoIons
            NoGrp = NoGrp + 1
            NoSub = NoSub + 1

            MainGroups(NoSub)  = NoGrp
            SubGroups(NoGrp+1) = NoSub 

            rkU(NoSub) = rk_Ion(No)
            qk(NoSub) = qk_Ion(No)
         END DO
      END IF

!------------------------------------------
!--- Building of the ínteraction matrix A
      ALLOCATE (MatrixA(NoGroups,NoGroups))
      MatrixA(:,:) = 0.D0

      !-- restore original UNIFAC matrix (Organics)
      MatrixA(1:NoGroupsOld,1:NoGroupsOld) = mU(:,:)

      !-- copy additional organics groups 
      DO No=1,NoAddGroups
         IF (TypeAddGroups(No) <= 0) CYCLE
         NoOld = Group2Group(No,1)
         NoNew = Group2Group(No,2)

         MatrixA(NoNew,:) = MatrixA(NoOld,:)
         MatrixA(:,NoNew) = MatrixA(:,NoOld)

         DO ig=1,NoAddGroups
            igOld = Group2Group(ig,1)
            igNew = Group2Group(ig,2)

            MatrixA(NoNew,igNew) = MatrixA(NoOld,igOld)
            MatrixA(igNew,NoNew) = MatrixA(igOld,NoOld)
         END DO
      END DO

      !-- copy ion-organics interactions into A
      ia = NoGroups - NoIons
      DO ig=1,NoAddGroups
         igNew = Group2Group(ig,2)

         DO No=1,NoIons
            NoNew = ia + No
            MatrixA(igNew,NoNew) = OW_IW(ig,No)
            MatrixA(NoNew,igNew) = IW_OW(No,ig)
         END DO
      END DO

!------------------------------------------
!--- Allocation and composition of species arrays
      NoMols = NoOrgs + NoIons
      IF (NoMols < 1 .OR. .NOT.AllocNoMols( ))   &
  &      CALL ErrorAndExit(ProgramName,'Invalid size!')

      !-- species names
      ia = NoOrgs + 1
      NameMols(1 :NoOrgs) = NameOrg(1:NoOrgs)
      NameMols(ia:NoMols) = NameIon(1:NoIons)

      !-- occurency of organics
      Occurrences   (1:NoOrgs,1:NoGroupsOld) = OccurOrg(:,:)
      SubOccurrences(1:NoOrgs,1:NoSubGOld)   = SubOccurOrg(:,:)
      DO No=1,NoOrgs
         IF (TypeOrg(No) <= 0) CYCLE
           
         DO iag=1,NoAddGroups
            IF (TypeOrg(No) /= TypeAddGroups(iag))  CYCLE
            DO ig=1,NoGroupsOld
               IF (Group2Group(iag,1) /= ig)  CYCLE
               igNew = Group2Group(iag,2)

               Occurrences(No,igNew) = Occurrences(No,ig)
               Occurrences(No,ig)    = 0

               SubOccurrences(No,SubGroups(igNew):SubGroups(igNew+1)-1) =   &
  &                           SubOccurrences(No,SubGroups(ig):SubGroups(ig+1)-1)
               SubOccurrences(No,SubGroups(ig):SubGroups(ig+1)-1) = 0
            END DO
         END DO
      END DO

      !-- occurency of ions
      DO No=1,NoIons
         igNew = NoGroups - NoIons + No
         Occurrences(NoOrgs+No,igNew) = 1
         
         igNew = NoSubGroups - NoIons + No
         SubOccurrences(NoOrgs+No,igNew) = 1
      END DO
!
!==========================================================================!
!---  Storage SubOccurrences in Sparse Form: Row- and Column-oriented form
!==========================================================================!
!
!--  column-oriented sparse structure 
   ALLOCATE(SubcColPtr(NoSubGroups+1))

   SubcColPtr(1)=1

!--  row pointer
   nonz_Subc = 0
   DO ig=1,NoSubGroups
      SubcColPtr(ig+1) = SubcColPtr(ig)
      DO No=1,NoMols
         IF (SubOccurrences(No,ig) > 0 ) THEN
            SubcColPtr(ig+1) = SubcColPtr(ig+1) + 1
            nonz_Subc = nonz_Subc + 1
         END IF
      END DO
   END DO

!--  column elements
   ALLOCATE(SubcRowInd(nonz_Subc))
   ALLOCATE(SubcSparse(nonz_Subc))
   inonz = 0
   DO ig=1,NoSubGroups
      DO No=1,NoMols
         IF (SubOccurrences(No,ig) > 0 ) THEN
            inonz = inonz+1
            SubcRowInd(inonz) = No
            SubcSparse(inonz) = SubOccurrences(No,ig)
         END IF
      END DO
   END DO

   IF (inonz /= nonz_Subc )    &
&     STOP  'UNIFAC..InitSystem: Error during packing of Occurences!'
!
!------------------------------------------------------------------------
!
!--  row-oriented sparse structure 
   ALLOCATE(SubrRowPtr(NoMols+1))

   SubrRowPtr(1)=1

!--  row pointer
   nonz_Subr = 0
   DO No=1,NoMols
      SubrRowPtr(No+1) = SubrRowPtr(No)
      DO ig=1,NoSubGroups
         IF (SubOccurrences(No,ig) > 0 ) THEN
            SubrRowPtr(No+1) = SubrRowPtr(No+1) + 1
            nonz_Subr = nonz_Subr + 1
         END IF
      END DO
   END DO

!--  column elements
   ALLOCATE(SubrColInd(nonz_Subr))
   ALLOCATE(SubrSparse(nonz_Subr))
   inonz = 0
   DO No=1,NoMols
      DO ig=1,NoSubGroups
         IF (SubOccurrences(No,ig) > 0 ) THEN
            inonz = inonz+1
            SubrColInd(inonz) = ig
            SubrSparse(inonz) = SubOccurrences(No,ig)
         END IF
      END DO
   END DO

   IF (inonz /= nonz_Subr )    &
&     STOP  'UNIFAC..InitSystem: Error during packing of SubOccurences!'
!
!------------------------------------------------------------------------
   END SUBROUTINE InitSystem

!==========================================================================

! ------------------------------------------------------------------------
   SUBROUTINE OrderSystem( )

! ------------------------------------------------------------------------
!--- organize communicartion between reaction scheme and UNIFAC
! -------------------------------------------------------------------------
!

!---  internal variables
      INTEGER :: i, ig, j, jt, No, indspc
      CHARACTER(20):: string

      INTEGER, PARAMETER :: Mx_Grp = 10
      INTEGER :: MGrps(Mx_Grp,3)


! -------------------------------------------------------------------------
! ---  Define set of Pitzer coefficients
! -------------------------------------------------------------------------
!
!---  Define index arrays 
      ALLOCATE (UnifacType(nAqua))
      ALLOCATE (UnifacBack(NoMols))

      UnifacType(:)  = 0
      WaterIndUni    = 0
      UnifacBack(:)  = 0
!
! -------------------------------------------------------------------------
!---  Set index transformation
! -------------------------------------------------------------------------
!
!---  water
      string = 'H2O'
      WaterIndUni = ifind(NoMols, string, NameMols)
      IF (WaterIndUni <= 0) THEN
         WRITE(*,*) 'InitUNIFAC ... Error: Water is not included in UNIFAC system!'
         STOP       'InitUNIFAC ... Error: Water is not included in UNIFAC system!'
      END IF

!---  forward transformation
      DO jt=1,nAqua
         string = ADJUSTL(SpeciesNameAqua(jt))
         CALL UpCase(string)

         indspc = ifind(NoMols, string, NameMols)
         IF (indspc > 0)  THEN
            UnifacInd(jt) = indspc
            IF (indspc <= NoOrgs)  UnifacType(jt) = TypeOrg(indspc)
         END IF
      END DO

!---  forward transformation
      DO No=1,NoMols
         DO jt=1,nAqua
            IF (UnifacInd(jt) == No)  THEN
               UnifacBack(No) = jt
               CYCLE
            END IF
         END DO
      END DO

! -------------------------------------------------------------------------
!---  Print UNIFAC system
! -------------------------------------------------------------------------
!
      WRITE(*,801) 
      WRITE(*,802) NoOrgs, NoGroups, NoIons, NoSubGroups,  &
&                  '  No.','  UnifacInd','    Type','     Species Name',  &
&                  '     Main Groups'

!---  water
      jt = 0
      DO ig=1,NoGroups
         IF (Occurrences(WaterIndUni,ig) <= 0 .OR. Mx_Grp <= jt)  CYCLE
         jt = jt + 1
         MGrps(jt,1) = ig
         MGrps(jt,2) = Occurrences(WaterIndUni,ig)
      END DO
      WRITE(*,803)  WaterIndUni, 'WATER             ',   &
&                   ((MGrps(j,i),i=1,2),j=1,jt)

!---  other species
      DO No=1,nAqua
         IF (UnifacInd(No) <= 0)  CYCLE
         jt = 0
         DO ig=1,NoGroups
            IF (Occurrences(UnifacInd(No),ig) <= 0 .OR. Mx_Grp <= jt)  CYCLE
            jt = jt + 1
            MGrps(jt,1) = ig
            MGrps(jt,2) = Occurrences(UnifacInd(No),ig)
            MGrps(jt,3) = Occurrences(UnifacInd(No),ig)
         END DO
         WRITE(*,804)  No, UnifacInd(No), UnifacType(No), SpeciesNameAqua(No),   &
&                      ((MGrps(j,i),i=1,2),j=1,jt)
      END DO
      WRITE(*,805)

! -------------------------------------------------------------------------
801   FORMAT(1x/1x,75('=')/ ' UNIFAC Initialization:' )
802   FORMAT(1x/1x,75('-')/           &
&            ' Considered Pitzer System: Organics =',i3,'     Groups    =',i3 /  &
&            '                           Ions     =',i3,'     Subgroups =',i3 // &
&            a5,a9,a9,a18,a20 / 1x,65('-'))
803   FORMAT(4x,i9,9x,7x,a18,4x,10(i5,' =>',i2))
804   FORMAT(i4,i9,i9,7x,a18,4x,10(i5,' =>',i2))
805   FORMAT(1x/1x,75('=')/1x)
!
! -------------------------------------------------------------------------
   END SUBROUTINE OrderSystem

!==========================================================================
!###########################################
!###   Subroutines
!###########################################
!
!==========================================================================

   SUBROUTINE NormalizeFractions

! -------------------------------------------
! ---   normalize fractions and copy to calc array
! -------------------------------------------
      REAL(RealKind) :: SumFractions
      INTEGER :: smd

      SumFractions = 0.D0
      DO smd=1,NoMols
         SumFractions = SumFractions + MolFractions(smd)
      END DO 
      IF (SumFractions < 0.99999)  MolFractions(NoMols) = 1.D0 - SumFractions

   END SUBROUTINE NormalizeFractions

!==========================================================================
   SUBROUTINE UpCase(String)
!
! -------------------------------------------
! ---    upcases a string
! -------------------------------------------
       CHARACTER*(*) :: String       
       INTEGER :: i

       DO i=1,LEN_TRIM(String)
          IF (String(i:i) >= 'a' .AND. String(i:i) <= 'z') THEN
             String(i:i) = char(ichar(String(i:i)) - 32)
          END IF      !          
       END DO      ! 
   END SUBROUTINE UpCase

!==========================================================================
   LOGICAL FUNCTION AllocNoMols( )
!
! ----------------------------------------------------
! ---  Allocation and initialization of species arrays
! ---------------------------------------------------

      INTEGER :: ios, iosum

!---------------------------------------------------
!---  allocation
      AllocNoMols = .TRUE.

      IF (LoopSystem == 0) THEN
         iosum = 0

         ALLOCATE (Occurrences(NoMols,NoGroups), STAT=ios)
         iosum = iosum + ios
         ALLOCATE (SubOccurrences(NoMols,NoSubGroups), STAT=ios)
         iosum = iosum + ios

         ALLOCATE (MolFractions (NoMols), STAT=ios)      !- real(8)
         iosum = iosum + ios

         ALLOCATE (NameMols(NoMols), STAT=ios)           !- species names
         iosum = iosum + ios

         !-- result arrays
         ALLOCATE (f_Inf(NoMols), f_Star(NoMols), STAT=ios)
         iosum = iosum + ios

         ALLOCATE (rU(NoMols), q(NoMols), l(NoMols), STAT=ios)
         iosum = iosum + ios
         ALLOCATE (s(NoSubGroups, NoMols), nue(NoSubGroups), STAT=ios)
         iosum = iosum + ios
         ALLOCATE (phi(NoMols), xi(NoMols), ln_gamma_c_ref(NoMols),ln_gamma_c(NoMols), STAT=ios)
         iosum = iosum + ios
         ALLOCATE (x(NoMols,NoSubGroups), TotalX(NoSubGroups), STAT=ios)
         iosum = iosum + ios
         ALLOCATE (theta(NoMols,NoSubGroups), Totaltheta(NoSubGroups), STAT=ios)
         iosum = iosum + ios
         ALLOCATE (Theta2(NoSubGroups), GammaU(NoMols,NoSubGroups), STAT=ios)
         iosum = iosum + ios
         ALLOCATE (Gamma_First(NoMols,NoSubGroups),Gamma_Second(NoMols,NoSubGroups), STAT=ios)
         iosum = iosum + ios
         ALLOCATE (TotalGamma(NoSubGroups), ln_gamma_r(NoMols), ln_gamma_r_ref(NoMols),STAT=ios)
!         iosum = iosum + ios
!         ALLOCATE (TotalGamma_First(NoSubGroups), TotalGamma_Second(NoSubGroups),ln_gamma_r_First(NoMols), ln_gamma_r_Second(NoMols),STAT=ios)
         iosum = iosum + ios        
         ALLOCATE (Results(NoMols), Results_ref(NoMols),STAT=ios)
         iosum = iosum + ios

         IF (iosum /= 0)  AllocNoMols = .FALSE.
      END IF

!---------------------------------------------------
!---  initialization
      SubOccurrences(:,:) = 0
      Occurrences(:,:)    = 0

      MolFractions(:)    = 0.D0

      !-- result arrays
      rU(:) = 0.D0 
      q(:) = 0.D0 
      l(:) = 0.D0

      s(:,:) = 0.D0 
      nue(:) = 0.D0
      phi(:) = 0.D0 
      xi(:)  = 0.D0 

      ln_gamma_c(:) = 0.D0
      ln_gamma_c_Ref(:) = 0.D0
      Totaltheta(:) = 0.D0

      x(:,:)     = 0.D0 
      TotalX(:)  = 0.D0
      theta(:,:) = 0.D0 
      Theta2(:)  = 0.D0 
      GammaU(:,:) = 0.D0

      TotalGamma(:) = 0.D0 
      ln_gamma_r(:) = 0.D0

      Results(:)  = 0.D0 
      Results_Ref(:)  = 0.D0 
!
! ---------------------------------------------------
   END FUNCTION AllocNoMols

!==========================================================================
!===   matrix: Reading and printing of transaction matrix 
!==========================================================================
!
  SUBROUTINE ReadMatrix(File)
!
!-------------------------------------------
!       writes matrix to file
!-------------------------------------------

!--- max columns per file line
    INTEGER, PARAMETER :: MAX_COLUMNS_PER_LINE = 10

    INTEGER :: File,ARows,AColumns,row,StartCol,EndCol,i,col
    REAL(RealKind) :: Buffer(MAX_COLUMNS_PER_LINE)
    CHARACTER(200) :: BufferString

!-------------------------------------------
!--- get header
    READ (File,*) AColumns,ARows

!--- allocation and initialization of matrix
    NoGroups = ARows
    ALLOCATE (SubGroups(NoGroups+1))
    ALLOCATE (mU(NoGroups,NoGroups))
    
    SubGroups(:) = 0
    mU(:,:)       = 0.E0

!-------------------------------------------
!--- rows 
    DO row = 1,ARows
      StartCol = 1
      DO WHILE (StartCol <= AColumns)
!       --- attempt to read  MAX_COLUMNS_PER_LINE vals, buffer in string
        READ (File,'(A)',End=88,Err=88) BufferString
88      CONTINUE
        READ (BufferString,*,End=99,Err=99) Buffer
99      CONTINUE
!       --- copy remaining values  
        EndCol = StartCol + MAX_COLUMNS_PER_LINE - 1
        IF (EndCol > AColumns) EndCol = AColumns
        i = 1
        DO col = StartCol,EndCol
           mU(row,col) = Buffer(i)
           i = i + 1
        END DO      ! copy 
        StartCol = StartCol + MAX_COLUMNS_PER_LINE
      END DO        ! read file lines         
    END DO          ! rows
!
  END SUBROUTINE  ReadMatrix

!======================================================================
  SUBROUTINE WriteMatrix(File,Matrix,ARows,AColumns)
!
!-------------------------------------------
!---  writes matrix to file
!-------------------------------------------

    INTEGER :: File,ARows,AColumns,row,col,no
    REAL(RealKind) :: Matrix(ARows,AColumns)

!--- max columns per file line
    INTEGER, PARAMETER :: MAX_COLUMNS_PER_LINE = 10

!-------------------------------------------
!--- header
    WRITE (File,'(2(I5,2X))') AColumns,ARows

!--- rows 
    DO row = 1,ARows
      no = 1
      DO col = 1,AColumns
        WRITE (File,'(F12.5,$)') Matrix(row,col) 
        no = no + 1
!    --- insert nl each maxcol or at end
        IF ((No > MAX_COLUMNS_PER_LINE) .OR.   &
 &                 (col == AColumns)) THEN
           WRITE (File,*)
           No = 1
        END IF      ! nl 
      END DO          
    END DO 
!
  END SUBROUTINE  WriteMatrix

!==========================================================================
!===   contribut: Reading of parameters
!==========================================================================
!
  SUBROUTINE ReadParameters()
!
! -------------------------------------------
! ---  copy the m-vectors to matrix m
! -------------------------------------------

     INTEGER :: No,Group,OldGroup
     LOGICAL :: EndOfFile
     CHARACTER(80) :: Buffer


! --- define temporary arrays for input
     INTEGER, PARAMETER :: MAX_SUB_GROUPS = 500 
     INTEGER :: MainGroups0(MAX_SUB_GROUPS)
     REAL(RealKind) :: rk0(MAX_SUB_GROUPS),qk0(MAX_SUB_GROUPS)

!-------------------------------------------
!--- open file
     IF (.NOT. FileOpen(MatrixFile,MatrixFileName,FileRead)) THEN
        CALL ErrorAndExit(ProgramName,'Cannot open matrix file!')
     END IF ! matrix file disok

!--- read name of tab & sub
     READ (MatrixFile,'(A)') SubtabName
     READ (MatrixFile,'(A)') SubName

!--- read matrix        
     CALL ReadMatrix(MatrixFile)

!--- read matrix        
!     CALL WriteMatrix(DetFile,m,NoGroups,NoGroups)

!-------------------------------------------
!--- read subgroups and set up index arrays
     EndOfFile   = .FALSE.
     NoSubGroups = 0
     OldGroup    = 0

     DO WHILE (.NOT. EndOfFile )
        READ (MatrixFile,'(A)',End=99,Err=99) Buffer
        IF (LEN_TRIM(Buffer) > 4) THEN
           NoSubGroups = NoSubGroups + 1
           IF (NoSubGroups > MAX_SUB_GROUPS)  &       
  &          STOP 'ReadParameters...Error: NoSubGroups to large! Increase MAX_SUB_GROUPS!'

           READ (Buffer,*) No,Group,rk0(NoSubGroups),qk0(NoSubGroups)
           !-- create index array group->subgroup
           IF (Group /= OldGroup) THEN
              SubGroups(group) = NoSubGroups
              OldGroup         = Group
           END IF
!          !-- create index array subgroup->group
           MainGroups0(NoSubGroups) = group
        END IF                                 ! string length OK
        CYCLE         
99      EndOfFile = .TRUE.
     END DO                                    ! .NOT. EndOfFile
     SubGroups(NoGroups+1) = NoSubGroups+1
     CLOSE (MatrixFile)

     WRITE (DetFile,'(1X,2(A,I3,1X))') 'Main groups: ',NoGroups,' Sub groups ',NoSubGroups

!-------------------------------------------
!--- allocation and setting of SubGroup arrays
     ALLOCATE (MainGroups(NoSubGroups))
     ALLOCATE (rkU(NoSubGroups),qk(NoSubGroups))

     MainGroups(:) = MainGroups0(1:NoSubGroups)

     rkU(:) = rk0(1:NoSubGroups)
     qk(:) = qk0(1:NoSubGroups)
!
!-------------------------------------------
  END SUBROUTINE ReadParameters
SUBROUTINE InitFinalMR

  CHARACTER(20) :: readaway 
  REAL(RealKind), ALLOCATABLE :: MolarWeight(:)
  INTEGER, ALLOCATABLE ::  Main(:), Sub(:)
  
  INTEGER :: IndSubMain(NoGroups)
  INTEGER :: i, j,jt

!   open(10,file='Data_Lifac/Extended_Final')
   OPEN(10,file='Data_Lifac/Extended_Test')            ! COOH Group values are set as Zero
!  open(10,file='Data_Lifac/InputMiddle_Final')

! --- Reading of Interaction Parameters
  j=0
  DO i=1,1000
    IF (i == 7) THEN
     READ(10,*,END=999) readaway, NoMainGrpMR
    ELSEIF (i == 8) THEN
     READ(10,*,END=999) readaway, NoSubGrpMR
    ELSEIF (i == 9) THEN
     READ(10,*,END=999) readaway, NoA
    ELSEIF (i == 10) THEN
     READ(10,*) readaway, NoC
       NoIonMR=NoA+NoC
 
       !-- allocation of arrays
    ALLOCATE  (b1ki(NoMainGrpMR,NoIonMR),   &  
&             b2ki(NoMainGrpMR,NoIonMR),    &
&             b3ki(NoMainGrpMR,NoIonMR),    &      
&             bca1(NoC,NoA),                &  
&             bca2(NoC,NoA),                &  
&             bca3(NoC,NoA),                &  
&             cca1(NoC,NoA),                &  
&             cca2(NoC,NoA),                &  
&             rcc(NoC,NoC),                 &  
&             qcc(NoC,NoA))                      
       
       ALLOCATE (IndBMRk(NoMainGrpMR))
       ALLOCATE (IndBackMRic(NoC))
       ALLOCATE (IndBackMRia(NoA))
       ALLOCATE (IonName(NoIonMR))
       ALLOCATE (NameCation(NoC))
       ALLOCATE (NameAnion(NoA))
       ALLOCATE (Main(NoSubGrpMR))
       ALLOCATE (Sub(NoSubGrpMR))
       ALLOCATE (MolarWeight(NoSubGrpMR))

     
       b1ki(:,:) = 0.E0
       b2ki(:,:) = 0.E0
       b3ki(:,:) = 0.E0
       bca1(:,:) = 0.E0
       bca2(:,:) = 0.E0
       bca3(:,:) = 0.E0
       cca1(:,:) = 0.E0
       cca2(:,:) = 0.E0
       rcc(:,:)  = 0.E0
       qcc(:,:)  = 0.E0
      
       
       IndBackMRic(:) = 0
       IndBackMRia(:)=0
       IndBMRk(:)=0
       IonName(:)=''
       NameCation(:)=''
       NameAnion (:)=''
       Main(:) = 0
       Sub (:) = 0
       MolarWeight(:) = 0.E0

    ELSEIF (i == 17) THEN
     READ(10,*,END=999) readaway, IndBMRk(1:NoMainGrpMR)
    ELSEIF ((i > 18) .and. (i < 43)) THEN
     j=j+1
     READ(10,*,END=999) IonName(j),b1ki(1:NoMainGrpMR,j)
    ELSEIF ((i > 49) .and. (i < 74)) THEN
     j=j+1  
     READ(10,*,END=999) readaway,b2ki(1:NoMainGrpMR,j)
    ELSEIF ((i > 80) .and. (i < 105)) THEN
     j=j+1  
     READ(10,*,END=999) readaway,b3ki(1:NoMainGrpMR,j)
    ELSEIF ((i == 109)) THEN
     READ(10,*,END=999) readaway, (NameAnion(j),j=1,NoA)
    ELSEIF ((i > 110) .and. (i < 122)) THEN
     j=j+1  
     READ(10,*,END=999) NameCation(j),bca1(j,1:NoA)
    ELSEIF ((i > 127) .and. (i < 138)) THEN
     j=j+1  
     READ(10,*,END=999) readaway,bca2(j,1:NoA)
    ELSEIF ((i > 143) .and. (i < 155)) THEN
     j=j+1  
     READ(10,*,END=999) readaway,bca3(j,1:NoA)
    ELSEIF ((i > 160) .and. (i < 172)) THEN
     j=j+1  
     READ(10,*,END=999) readaway,cca1(j,1:NoA)
    ELSEIF ((i > 177) .and. (i < 188)) THEN
     j=j+1  
     READ(10,*,END=999) readaway,cca2(j,1:NoA)
    ELSEIF ((i > 193) .and. (i < 205)) THEN
     j=j+1  
     READ(10,*,END=999) readaway,rcc(j,1:NoC) 
    ELSEIF ((i > 211) .and. (i < 223)) THEN
     j=j+1  
     READ(10,*,END=999) readaway,qcc(j,1:NoA)
    ELSEIF ((i > 228) .and. (i < 247)) THEN
     j=j+1
     READ(10,*,END=999) Main(j),Sub(j),MolarWeight(j)
    ELSE 
     READ(10,*,END=999)
     j=0
    ENDIF
  ENDDO
999  WRITE (*,*) 'InitMiddle: End of File reached! No. of Lines =', i 
  CLOSE (10)

! --- Allocation of index arrays
  ALLOCATE(IndMRk(NoSubGroups))  
  ALLOCATE(IndMRi(nAqua))
  ALLOCATE(IndMRia(nAqua))
  ALLOCATE(IndMRic(nAqua))
  ALLOCATE(MMk(NoSubGroups))
 
  IndMRk(:)=0
  IndMRi(:)=0
  IndMRia(:)=0
  IndMRic(:)=0  
  MMk(:) = 0

! --- Calculating of Indices for the different field matrices of coefficients  
! ------Index for Subgroups 
  IndSubMain(:)=0
  DO i=1,NoMainGrpMR
    IndSubMain(IndBMRk(i))=i
  ENDDO
  DO i=1,NoSubGrpMR
    IndMRk(Sub(i))=IndSubMain(Main(i))
  ENDDO
 
! ------Index for Ions
  DO i=1,nAqua
    DO j=1,NoIonMR
      IF (SpeciesNameAqua(i) /= IonName(j)) CYCLE
        IndMRi(i)=j
    ENDDO
  ENDDO

! ------Index for Cations
  DO i=1,nAqua
    DO j=1,NoC
      IF (SpeciesNameAqua(i) /= NameCation(j)) CYCLE
        IndBackMRic(j)=i
        IndMRic(i)=j       
    ENDDO         
  ENDDO

! ------Index for Anions
  DO i=1,nAqua
    DO j=1,NoA
      IF (SpeciesNameAqua(i) /= NameAnion(j)) CYCLE
        IndBackMRia(j)=i
        IndMRia(i)=j
    ENDDO
  ENDDO
    
! ------Index for MolarWeight
  DO i=1,NoSubGrpMR
     MMk(Sub(i))=MolarWeight(i)
  ENDDO
! ------Deallocate Arrays
  DEALLOCATE (MolarWeight)
  DEALLOCATE (IndBMRk, Main, Sub)
  DEALLOCATE (IonName, NameCation, NameAnion)

END SUBROUTINE InitFinalMR

SUBROUTINE NeuMR(ca,gamma_MR)
!
!=========================================================================
!===  Determine Middle-Range activity coefficients
!=========================================================================
!

!---  external variables
   REAL(RealKind) :: ca(:)      ! Molalities [mol/l]
   REAL(RealKind) :: gamma_MR(:) 

!---  internal variables
!---  for counting
   INTEGER :: jt, kt, jtc, jta, No, mol, k,i,j
            
!---  others
   REAL(RealKind) :: SumMol, aIonStr, MMs, AllOccurences

!*****************************************************
!----------- For AIOMFAC Interactions
!*****************************************************
   REAL(RealKind) :: BKI(NoMainGrpMR,NoIonMR),     & ! Binary interaction parameter of LIFAC (MainGroups, Ions)
&             BKIStr(NoMainGrpMR,NoIonMR),  & ! derivative of B with respect to ionic strength 
&             BBcatan(NoC,NoA),             & ! Binary interaction parameter of LIFAC (Cation, Anion)
&             BBcatanStr(NoC,NoA),          & ! derivative of Bcatan with respect to ionic strength
&             CCcatan(Noc,NoA),             & ! Interaction papramenter of LIFAC (Cation, Anion)
&             CCcatanStr(NoC,NoA)             ! derivative of CCcatan with respect to Ionic strength
!*****************************************************
!----------- For mod.LIFAC Interactions
!*****************************************************
   REAL(RealKind) :: BKI_Li(NoMainGrpMR,NoIonMR),     & ! Binary interaction parameter of LIFAC (MainGroups, Ions)
&             BKIStr_Li(NoMainGrpMR,NoIonMR),  & ! derivative of B with respect to ionic strength 
&             BBcatan_Li(NoC,NoA),             & ! Binary interaction parameter of LIFAC (Cation, Anion)
&             BBcatanStr_Li(NoC,NoA)             ! derivative of Bcatan with respect to ionic strength


   REAL(RealKind) :: SumIon1,Sumion2,SumIon3,SumIon4,SumIon5,SumIon6,SumIon7, SumIon8,SumIon16
   REAL(RealKind) :: SumIon1A,Sumion2A,SumIon3A,SumIon4A,SumIon5A,SumIon6A,SumIon7A ,SumIon8A
   REAL(RealKind) :: Charge1A,Charge2A,Charge1C, Charge2C,sum16ina,sumion16a,sum6ina
   REAL(RealKind) :: SumSolv1,SumSolv2,SumSolv3,SumSolv4, SumSolv5, SumSolv6,Sum4inA, Sum16in
   REAL(RealKind) :: sum6a,sum5a,sum7a,TotalCCA, Totalca	
   REAL(RealKind) :: LnGammaIon, LnGammaSpecie  
   REAL(RealKind) :: lnGammaMRk(NoSubGroups)
   REAL(RealKind) :: Sum2,sum3, sum4,sum5, Sum2in,Sum3in,Sum4j,Sum4in,Sum4k,Sum6in, Sum8in, Sum2inA
   REAL(RealKind) :: MolFracSubGrp(NoSubGroups)
   REAL(RealKind) :: MolFracMR(NoOrgs) ! solvent mol fractions 
    
!
! ------------------------------------------------------------------------
!

!    Temperature = Temp             ! temperature
!
!------------------------------------------------------------------------
!--- compute activity coefficients of ideal solution for normalizing
        MolFracMR  = 0.E0
        MolFracSubGrp= 0.E0
!--- compute mol fraction activity coefficients

     !--- transformation of molalities to mol fractions
        SumMol = MolH2O
        DO No=1,NoOrgs
          IF (UnifacBack(No) <= 0)  CYCLE
          SumMol = SumMol + ca(UnifacBack(No))
        END DO

        DO No=1,NoOrgs
          IF (UnifacBack(No) <= 0)  CYCLE
          MolFracMR(No) = ca(UnifacBack(No)) / SumMol 
        END DO
        MolFracMR(WaterIndUni) = MolH2O / SumMol

!???       CALL NormalizeFractions    ???
       
     !--- computation of mol fractions of subgroups
        DO kt=1,NoSubGroups
          DO mol=SubcColPtr(kt),SubcColPtr(kt+1)-1
            IF (SubcRowInd(mol) > NoOrgs) CYCLE
            AllOccurences=AllOccurences+SubcSparse(mol)*MolFracMR(SubcRowInd(mol))
          ENDDO
        ENDDO

     !--- H2O is included in subgroups 
        DO kt=1,NoSubGroups
          Sum2in=0.E0   
          DO mol=SubcColPtr(kt),SubcColPtr(kt+1)-1
            IF (SubcRowInd(mol) > NoOrgs) CYCLE
            Sum2in=Sum2in+SubcSparse(mol)*MolFracMR(SubcRowInd(mol))
          ENDDO
          MolFracSubGrp(kt)=Sum2in/AllOccurences       ! x' für Subgruppe k in Sum3 (12)
        ENDDO

!------------------------------------------------------------------------
!--- Ionic Strength
        aIonStr = 0.E0
        DO jt=1,SIZE(ca)
          IF (Charge(jt) == 0) CYCLE
          aIonStr = aIonStr + Charge(jt) * Charge(jt) *ca(jt)
        END DO
        aIonStr = 0.5E0 * aIonStr 
!
!--- Matrices for AIOMFAC interaction parameters
!
!************************************************************************
!*************** Write tha code according to the A. Zuend ***************
!***** Based on MolFraction scale, have to convert in to Molalities *****
!************************************************************************

!----Binary interaction coefficients between solvent main groups and ions
        BKI(:,:)        = 0.e0
        BBcatan(:,:)    = 0.e0
        CCcatan(:,:)    = 0.e0
        BKIStr(:,:)     = 0.e0     
        BBcatanStr(:,:) = 0.e0       
        CCcatanStr(:,:) = 0.e0
!-------
        DO k = 1, NoMainGrpMR-6
          DO i = 1, NoA+NoC
            BKI(k,i)    = b1ki(k,i) + b2ki(k,i)* exp(-1.E0 * b3ki(k,i)* sqrt(aIonStr))
            BKIStr(k,i) = (-1.E0 * b2ki(k,i)*b3ki(k,i)* exp(-1.E0 * b3ki(k,i)* sqrt(aIonStr)))/ (2.E0 *sqrt(aIonStr) )
          END DO
        END DO
!----Binary interaction coefficients between cations and anions     
        DO k= 1, NoC-6
          DO i=1,NoA-5
            BBcatan(k,i)    = bca1(k,i) + bca2(k,i)* exp(-1.E0 * bca3(k,i)* sqrt(aIonStr))
            BBcatanStr(k,i) = (-1.E0 * bca2(k,i)*bca3(k,i)* exp(-1.E0 * bca3(k,i)* sqrt(aIonStr)))/ (2.E0 *sqrt(aIonStr) )
          END DO
        END DO
!----Interaction coefficients between cation-anion pairs with respect to the total charge concentration
        DO k= 1, NoC-6
          DO i=1,NoA-5
            CCcatan(k,i)    =  cca1(k,i)* exp(-1.E0 * cca2(k,i)* sqrt(aIonStr))
            CCcatanStr(k,i) = -1.E0 * cca1(k,i)*cca2(k,i)* exp(-1.E0 * cca2(k,i)* sqrt(aIonStr))/ (2.E0 *sqrt(aIonStr) )
          END DO
        END DO
!******************************************************************************
!********************* Extended model *****************************************
!************** Calculating Matrices for LIFAC interaction parameters *********
!******************************************************************************
        BKI_Li(:,:)        = 0.e0
        BBcatan_Li(:,:)    = 0.e0
        BKIStr_Li(:,:)     = 0.e0     
        BBcatanStr_Li(:,:) = 0.e0       

       DO k= NoMainGrpMR-6,NoMainGrpMR
         DO i=1,NoA+NoC
           BKI_Li(k,i)    = b1ki(k,i) + b2ki(k,i)*exp(-1.2*sqrt(aIonStr + 0.13*aIonStr))
           BKIStr_Li(k,i) = b2ki(k,i)*exp(-1.2*sqrt(aIonStr) + 0.13*aIonStr)*(-0.6/sqrt(aIonStr) + 0.13)   
         ENDDO
       ENDDO

       DO k=NoC-6, NoC
         DO i=NoA-5, NoA
           BBcatan_Li(k,i)    = bca1(k,i) + bca2(k,i) *exp(-1.*sqrt(aIonStr) + 0.13*aIonStr)
           BBcatanStr_Li(k,i) = bca2(k,i)*exp(-1.*sqrt(aIonStr)+0.13*aIonStr)*(-0.5/sqrt(aIonStr)+0.13)
         ENDDO
       ENDDO
!========================================================================================================================================
!--- Molecular weight of solvent mixture
! To Change
!!!    MMs = MolFracSubGrp(WaterIndUni)*MM_H2O  ! in subgroups included
       MMs = 0.E0
       DO kt=1,NoSubGroups
         MMs=MMs + MolFracSubGrp(kt)*MMk(kt)
       ENDDO
!----------------------------
        Totalca = 0.E0
        DO jt=1,SIZE(ca)
          IF (IndMRi(jt) == 0) CYCLE !
          Totalca = Totalca +  (ca(jt)* ABS(Charge((jt))))
        END DO

!-----------------------------
!--- for solvent main group k*
        lnGammaMRk = 0.E0

!***********************************************************************************
!----SUMMATION OF ALL TERMS
!***********************************************************************************
!--- Set LnGammaMR of Subgroup
       DO k = 1, NoSubGroups-39
         IF (IndMRk(k) == 0) CYCLE
        
         SumSolv1 = 0.E0
         SumSolv2 = 0.E0
       
!Sum 1 in equation (3.17)
!------------------------------      
         
         DO jt = 1, SIZE(ca)
           IF (IndMRi(jt) == 0) CYCLE           
           SumSolv1 = SumSolv1 + BKI(IndMRk(k),IndMRi(jt))*ca(jt)
         END DO
         

!-----Sum 2 in equation (3.17)  
         DO kt=1,NoSubGroups
           IF (IndMRk(kt) == 0) CYCLE
           Sum2 = 0.E0
           DO jt=1,SIZE(ca)
             IF (IndMRi(jt) == 0) CYCLE
             Sum2 = Sum2 + (BKI(IndMRk(kt),IndMRi(jt))+aIonStr*BKIStr(IndMRk(kt),IndMRi(jt)))  &
&                              * ca(jt)
           ENDDO
           SumSolv2 = SumSolv2 + Sum2 * MolFracSubGrp(kt)
         ENDDO

!------- Summation terms for the calculation of solvent activity coefficients according to equation 17
!------------------------------       
        SumSolv3  = 0.E0 
        DO jtc=1,NoC-5
          IF (IndBackMRic(jtc) == 0)  CYCLE
          Sum3 = 0.E0
          DO jta=1,NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            Sum3 = Sum3 + (BBcatan(jtc,jta)+(aIonStr*BBcatanStr(jtc,jta)))   &
&                          *  ca(IndBackMRic(jtc))*ca(IndBackMRia(jta))
          ENDDO
          SumSolv3=SumSolv3+Sum3
        ENDDO

        DO jtc=1,NoC-5
          IF (IndBackMRic(jtc) == 0)  CYCLE
          Sum4 = 0.
          DO jta=1,NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            Sum4 = Sum4 + ((2.E0 * CCcatan(jtc,jta))+(aIonStr*CCcatanStr(jtc,jta)))   &
&                          *  ca(IndBackMRic(jtc))*ca(IndBackMRia(jta))
          ENDDO
          SumSolv4 = 0.E0
          DO jt=1,SIZE(ca)
            IF (IndMRi(jt) == 0) CYCLE
            SumSolv4 = SumSolv4+ Sum4 * ca(jt) * ABS(Charge(jt))
          END DO
        ENDDO
! ---------------------------
        SumSolv5 = 0.E0
        DO jtc = 1, NoC-5
          IF (IndBackMRic(jtc) == 0)  CYCLE
          DO jt = 1, NoC-8
          IF (IndBackMRic(jt) == 0)  CYCLE
            SumSolv5 = SumSolv5 + Rcc(jtc,jt)*ca(IndBackMRic(jtc))*ca(IndBackMRic(jt))
          END DO
        END DO
!----------------------------
        SumSolv6 = 0.E0
        DO jtc = 1, NoC-5 
          IF (IndBackMRic(jtc) == 0)  CYCLE
          DO jta = 1 , NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            SumSolv6 = SumSolv6 +2.E0* Qcc(jtc,jta)*ca(IndBackMRic(jtc))*ca(IndBackMRic(jtc))*ca(IndBackMRia(jta))
          END DO
        END DO


!--- Set LnGammaMR of Subgroup        
          lnGammaMRk(k)= SumSolv1 - (MMk(k)*(SumSolv2/MMs)) - MMk(k)*(SumSolv3 + SumSolv4+ SumSolv5 + SumSolv6)
       END DO
!===   H2O: Set MR water activity
       gamma_MR(iWater) = EXP(lnGammaMRk(IndSubWater))


!-----------------------------------------------------------------------------
!--- Calculating Gamma for Species k according to Equation 18
   
!       DO k=1,NoOrgs  ! because water is Number one in UNIFAC, Variables are from UNIFAC
!         IF (UnifacBack(k) <= 0)  CYCLE
!         lnGammaSpecie= 0.E0
!         DO i=SubrRowPtr(k),SubrRowPtr(k+1)-1
!           lnGammaSpecie= lnGammaSpecie+SubrSparse(i)*lnGammaMRk(SubrColInd(i))
!         ENDDO
!         gamma_MR(UnifacBack(k))= EXP(lnGammaSpecie)
!       ENDDO

!-- Set MR activity coefficient of water
       gamma_MR(iWater) = EXP(lnGammaMRk(IndSubWater)) 


!***********************************************************************************
!----SUMMATION OF ALL TERMS ACCORDING TO mod.LIFAC
!***********************************************************************************
!--- Set LnGammaMR of Subgroup
       DO k =  NoSubGroups+6, NoSubGroups
         IF (IndMRk(k) == 0) CYCLE
        
         SumSolv1 = 0.E0
         SumSolv2 = 0.E0
       
!Sum 1 in equation (3.17)
!------------------------------      
         DO jt = 1, SIZE(ca)
           IF (IndMRi(jt) == 0) CYCLE           
           SumSolv1 = SumSolv1 + BKI_Li(IndMRk(k),IndMRi(jt))*ca(jt)
         END DO
         

!-----Sum 2 in equation (3.17)  
         DO kt=NoSubGroups+6, NoSubGroups
           IF (IndMRk(kt) == 0) CYCLE
           Sum2 = 0.E0
           DO jt=1,SIZE(ca)
             IF (IndMRi(jt) == 0) CYCLE
             Sum2 = Sum2 + (BKI_Li(IndMRk(kt),IndMRi(jt))+aIonStr*BKIStr_Li(IndMRk(kt),IndMRi(jt)))  &
&                              * ca(jt)
           ENDDO
           SumSolv2 = SumSolv2 + Sum2 * MolFracSubGrp(kt)
         ENDDO

!------- Summation terms for the calculation of solvent activity coefficients according to equation 17
!------------------------------       
         SumSolv3  = 0.E0 
         DO jtc=NoC+6, NoC
           IF (IndBackMRic(jtc) == 0)  CYCLE
           Sum3 = 0.E0
           DO jta= NoA+5, NoA
             IF (IndBackMRia(jta) == 0)  CYCLE
             Sum3 = Sum3 + (BBcatan_Li(jtc,jta)+(aIonStr*BBcatanStr_Li(jtc,jta)))   &
&                          *  ca(IndBackMRic(jtc))*ca(IndBackMRia(jta))
           ENDDO
           SumSolv3=SumSolv3+Sum3
         ENDDO
!--- Set LnGammaMR of Subgroup        
          lnGammaMRk(k)= SumSolv1 - (MMk(k)*(SumSolv2/MMs)) - MMk(k)*(SumSolv3)
       END DO
!-----------------------------------------------------------------------------
!--- Calculating Gamma for Species k according to Equation 18
   
       DO k=1,NoOrgs  ! because water is Number one in UNIFAC, Variables are from UNIFAC
         IF (UnifacBack(k) <= 0)  CYCLE
         lnGammaSpecie= 0.E0
         DO i=SubrRowPtr(k),SubrRowPtr(k+1)-1
           lnGammaSpecie= lnGammaSpecie+SubrSparse(i)*lnGammaMRk(SubrColInd(i))
         ENDDO
         gamma_MR(UnifacBack(k))= EXP(lnGammaSpecie)
       ENDDO

!-- Set MR activity coefficient of water
       gamma_MR(iWater) = EXP(lnGammaMRk(IndSubWater)) 
       
!********************************************************************
!-------------------- FOR CATIONS ACCORDING TO AIOMFAC---------------
!********************************************************************

!--  For equation 19 
      DO j=1,SIZE(ca)
         IF (IndMRic(j) == 0) CYCLE 
     
         SumIon1=0.
         SumIon2=0.
         SumIon3=0.
!----------------------------      
! SumIon1 in (19 )
         DO kt=1,NoSubGroups-39
           IF (IndMRk(kt) == 0) CYCLE
           SumIon1=SumIon1+BKI(IndMRk(kt),IndMRic(j))*MolFracSubGrp(kt)
         ENDDO
!----------------------------      
! SumIon2 in (19)
         DO kt=1,NoSubGroups-39
           IF (IndMRk(kt) == 0) CYCLE
           Sum2in=0.
           DO jt=1,SIZE(ca)
             IF (IndMRi(jt) == 0) CYCLE
             Sum2in=Sum2in+BKIStr(IndMRk(kt),IndMRi(jt))*MolFracSubGrp(kt)*ca(jt)
           ENDDO
           SumIon2=SumIon2+Sum2in
         ENDDO
!----------------------------      
! SumIon3 in (19)          
         DO jta=1,NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            SumIon3=SumIon3+BBcatan(IndMRic(j),jta)*ca(IndBackMRia(jta))
         ENDDO
!----------------------------------------------------------------------------
!--- Equation 20 
        SumIon4  = 0.E0  
        DO jtc=1, NoC-5
          IF (IndBackMRic(jtc) == 0)  CYCLE
          Sum4in = 0.E0
          DO jta=1, NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            Sum4in = Sum4in + BBcatanStr(jtc,jta)*ca(IndBackMRic(jtc))*ca(IndBackMRia(jta))
          ENDDO
          SumIon4=SumIon4+Sum4in
        ENDDO
!----------------------------
! SumIon5 in (19)
!---------------------------- 
         SumIon5= 0.E0
         Sum5 = 0.E0          
         DO jta=1,1 , NoA-8
           IF (IndBackMRia(jta) == 0)  CYCLE         
           Sum5 = Sum5 + CCcatan(IndMRic(j),jta)*ca(IndBackMRia(jta))
         END DO
         DO jt=1,SIZE(ca)
           IF (IndMRi(jt) == 0) CYCLE
           SumIon5 =SumIon5+ Sum5 * ca(jt) * ABS(Charge(jt))
         END DO

!---------------------------- 
! SumIon6 in (19) Splitted
!---------------------------- 
        SumIon6 = 0.E0
        DO jtc = 1, NoC-5
          IF(IndBackMRic(jtc) == 0)  CYCLE 
          Sum6in = 0.E0
          DO jta=1, NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            Charge1A = ABS(Charge(IndBackMRia(jta)))*CCcatan(jtc,jta)
            Sum6in = Sum6in + (Charge1A *ca(IndBackMRic(jtc))*ca(IndBackMRia(jta)))
          END DO
          SumIon6 = SumIon6 + Sum6in
        END DO
!---------------------------- 
! SumIon6 in (19) Splitted
!---------------------------- 
        SumIon16 = 0.E0
        DO jtc = 1, NoC-5
          IF(IndBackMRic(jtc) == 0)  CYCLE 
          Sum16in = 0.E0
          DO jta=1,1 , NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            Charge2C = 0.5E0*CCcatanStr(jtc,jta)*Charge(IndBackMRic(jtc))**2.E0 
            Charge2C = Charge2C * Totalca 
            Sum16in = Sum16in + (Charge2C *ca(IndBackMRic(jtc))*ca(IndBackMRia(jta)))
          END DO
          SumIon16 = SumIon16 +Sum16in 
         END DO

!----------------------------
! SumIon7 in (19)
!---------------------------- 
         SumIon7= 0.E0
         DO jtc=1, NoC-5
           IF (IndBackMRic(jtc) == 0)  CYCLE
           SumIon7=SumIon7+Rcc(IndMRic(j),jtc)*ca(IndBackMRic(jtc))
         ENDDO    

!----------------------------
! SumIon8 in (19)
!---------------------------- 
         SumIon8= 0.E0

         DO jtc=1, NoC-5
           IF (IndBackMRic(jtc) == 0)  CYCLE
           Sum8in = 0.
           DO jta=1,NoA-8
             IF (IndBackMRia(jta) == 0)  CYCLE
             Sum8in = Sum8in + Qcc(jtc,jta)*ca(IndBackMRic(jtc))*ca(IndBackMRia(jta))
           ENDDO
           SumIon8=SumIon8 + Sum8in
         ENDDO

!------------------------------------------------------------------------
! Set LnGammaMR of Ions
      LnGammaIon = (SumIon1/MMs) + 0.5*((Charge(j)**2.E0)/MMs)*SumIon2 + SumIon3 + 0.5*(Charge(j)**2.E0)*SumIon4
      LnGammaIon=  LnGammaIon+SumIon5+SumIon6+SumIon16+SumIon7+SumIon8
      gamma_MR(j)=EXP(LnGammaIon)
!*******************************************************************************
      ENDDO
!++++++++++++++++++++++++++++++ End of Cation Loop  +++++++++++++++++++++++++++

!*****************************************************
!-------------------- FOR ANIONS ---------------------
!-----------------------------------------------------
      DO j=1,SIZE(ca)
          IF (IndMRia(j) == 0) CYCLE 

         SumIon1A=0.
         SumIon2A=0.
         SumIon3A=0.
!----------------------------      
! SumIon1 in (20 )
         DO kt=NoSubGroups+6, NoSubGroups
           IF (IndMRk(kt) == 0) CYCLE
           SumIon1A=SumIon1A+BKI(IndMRk(kt),IndMRia(j))*MolFracSubGrp(kt)
         ENDDO
!----------------------------      
! SumIon2A in (20)
         DO kt=NoSubGroups+6, NoSubGroups
           IF (IndMRk(kt) == 0) CYCLE
           Sum2inA=0.
           DO jt=1,SIZE(ca)
             IF (IndMRi(jt) == 0) CYCLE
             Sum2inA=Sum2inA+BKIStr(IndMRk(kt),IndMRi(jt))*MolFracSubGrp(kt)*ca(jt)
           ENDDO
           SumIon2A=SumIon2A+Sum2inA
         ENDDO
!----------------------------      
! SumIon3 in (20)
         DO jtc=1, NoC-5
           IF (IndBackMRic(jtc) == 0)  CYCLE
           SumIon3A=SumIon3A+BBcatan(jtc,IndMRia(j))*ca(IndBackMRic(jtc))
         ENDDO
!----------------------------------------------------------------------------
!--- Equation 20 
        SumIon4A  = 0.E0  
        DO jtc=1, NoC-5
          IF (IndBackMRic(jtc) == 0)  CYCLE
          Sum4inA = 0.E0
          DO jta=1, NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            Sum4inA = Sum4inA + BBcatanStr(jtc,jta)*ca(IndBackMRia(jta))*ca(IndBackMRic(jtc))
          ENDDO
          SumIon4A=SumIon4A+Sum4inA
        ENDDO
!----------------------------
! SumIon5 in (20)
!---------------------------- 
         SumIon5A= 0.E0  
         Sum5A= 0.E0
         DO jtc=1, NoC-5
           IF (IndBackMRic(jtc) == 0)  CYCLE         
           Sum5A = Sum5A + CCcatan(jtc,IndMRia(j))*ca(IndBackMRic(jtc))
         END DO
         DO jt=1,SIZE(ca)
           IF (IndMRi(jt) == 0) CYCLE
           SumIon5A =SumIon5A+ Sum5A * ca(jt) * ABS(Charge(jt))
         END DO

!---------------------------- 
! SumIon6 in (20) Splitted
!---------------------------- 
        SumIon6A = 0.E0
        DO jtc = 1, NoC-5
          IF(IndBackMRic(jtc) == 0)  CYCLE 
          Sum6inA = 0.E0
          DO jta=1, NoA-8
           IF (IndBackMRia(jta) == 0)  CYCLE
           Charge1A = ABS(Charge(IndBackMRia(jta)))*CCcatan(jtc,jta)
           Sum6inA = Sum6inA + (Charge1A *ca(IndBackMRic(jtc))*ca(IndBackMRia(jta)))
          END DO
          SumIon6A = SumIon6A + Sum6inA
         END DO
!---------------------------- 
! SumIon6 in (20) Splitted
!---------------------------- 
        SumIon16A = 0.E0
        DO jtc = 1, NoC-5
          IF(IndBackMRic(jtc) == 0)  CYCLE 
          Sum16inA = 0.E0
          DO jta=1, NoA-8
            IF (IndBackMRia(jta) == 0)  CYCLE
            Charge2A = 0.5E0*CCcatanStr(jtc,jta)*(Charge(IndBackMRia(jta))**2.E0 )
            Charge2A = Charge2A * Totalca 
            Sum16inA = Sum16inA + (Charge2A *ca(IndBackMRic(jtc))*ca(IndBackMRia(jta)))
          END DO
          SumIon16A = SumIon16A +Sum16inA
         END DO

!----------------------------
! SumIon8 in (19)
!---------------------------- 
         SumIon7A= 0.E0

         DO jtc=1, NoC-5
           IF (IndBackMRic(jtc) == 0)  CYCLE
           Sum7A= 0.
           DO jta=1, NoA-8
             IF (IndBackMRia(jta) == 0)  CYCLE
             Sum7A = Sum7A + Qcc(jtc,jta)*ca(IndBackMRic(jtc))*ca(IndBackMRia(jta))
           ENDDO
           SumIon7A=SumIon7A + Sum7A
         ENDDO

!------------------------------------------------------------------------
! Set LnGammaMR of Ions
         LnGammaIon = (SumIon1A/MMs) + 0.5*((Charge(j)**2.E0)/MMs)*SumIon2A + SumIon3A + 0.5*(Charge(j)**2.E0)*SumIon4A + SumIon5A
         LnGammaIon = LnGammaIon+SumIon6A+SumIon16A+SumIon7A
         gamma_MR(j)=EXP(LnGammaIon)

      ENDDO
!++++++++++++++++++++++++++++++ End of Anions Loop  +++++++++++++++++++++++++++
!====================================================================
!-------------------- EXTENSION PERFORMS HERE------------------------
!====================================================================

!********************************************************************
!-------------------- FOR CATIONS ACCORDING TO mod.LIFAC ------------
!********************************************************************

!--  For equation 19 
      DO j=1,SIZE(ca)
         IF (IndMRic(j) == 0) CYCLE 
     
         SumIon1=0.
         SumIon2=0.
         SumIon3=0.
!----------------------------      
! SumIon1 in (19 )
         DO kt=1,NoSubGroups
           IF (IndMRk(kt) == 0) CYCLE
           SumIon1=SumIon1+BKI_Li(IndMRk(kt),IndMRic(j))*MolFracSubGrp(kt)
         ENDDO
!----------------------------      
! SumIon2 in (19)
         DO kt=1,NoSubGroups
           IF (IndMRk(kt) == 0) CYCLE
           Sum2in=0.
           DO jt=1,SIZE(ca)
             IF (IndMRi(jt) == 0) CYCLE
             Sum2in=Sum2in+BKIStr_Li(IndMRk(kt),IndMRi(jt))*MolFracSubGrp(kt)*ca(jt)
           ENDDO
           SumIon2=SumIon2+Sum2in
         ENDDO
!----------------------------      
! SumIon3 in (19)          
         DO jta= NoA+5, NoA
            IF (IndBackMRia(jta) == 0)  CYCLE
            SumIon3=SumIon3+BBcatan_Li(IndMRic(j),jta)*ca(IndBackMRia(jta))
         ENDDO
!----------------------------------------------------------------------------
!--- Equation 20 
        SumIon4  = 0.E0  
        DO jtc= NoC+6, NoC
          IF (IndBackMRic(jtc) == 0)  CYCLE
          Sum4in = 0.E0
          DO jta= NoA+5, NoA
            IF (IndBackMRia(jta) == 0)  CYCLE
            Sum4in = Sum4in + BBcatanStr_Li(jtc,jta)*ca(IndBackMRic(jtc))*ca(IndBackMRia(jta))
          ENDDO
          SumIon4=SumIon4+Sum4in       
       END DO

!------------------------------------------------------------------------
! Set LnGammaMR of Ions
      LnGammaIon = (SumIon1/MMs) + 0.5*((Charge(j)**2.E0)/MMs)*SumIon2 + SumIon3 + 0.5*(Charge(j)**2.E0)*SumIon4     
      gamma_MR(j)=EXP(LnGammaIon)
!*******************************************************************************
      ENDDO
!++++++++++++++++++++++++++++++ End of Cation Loop  +++++++++++++++++++++++++++
!****************************************************************
!-------------------- FOR ANIONS according to mod.LIFAC----------
!****************************************************************
      DO j=1,SIZE(ca)
          IF (IndMRia(j) == 0) CYCLE 

         SumIon1A=0.
         SumIon2A=0.
         SumIon3A=0.
!----------------------------      
! SumIon1 in (20 )
         DO kt=1,NoSubGroups
           IF (IndMRk(kt) == 0) CYCLE
           SumIon1A=SumIon1A+BKI_Li(IndMRk(kt),IndMRia(j))*MolFracSubGrp(kt)
         ENDDO
!----------------------------      
! SumIon2A in (20)
         DO kt=1,NoSubGroups
           IF (IndMRk(kt) == 0) CYCLE
           Sum2inA=0.
           DO jt=1,SIZE(ca)
             IF (IndMRi(jt) == 0) CYCLE
             Sum2inA=Sum2inA+BKIStr_Li(IndMRk(kt),IndMRi(jt))*MolFracSubGrp(kt)*ca(jt)
           ENDDO
           SumIon2A=SumIon2A+Sum2inA
         ENDDO
!----------------------------      
! SumIon3 in (20)
         DO jtc=NoC+6, NoC
           IF (IndBackMRic(jtc) == 0)  CYCLE
           SumIon3A=SumIon3A+BBcatan_Li(jtc,IndMRia(j))*ca(IndBackMRic(jtc))
         ENDDO
!----------------------------------------------------------------------------
!--- Equation 20 
        SumIon4A  = 0.E0  
        DO jtc=NoC+6, NoC
          IF (IndBackMRic(jtc) == 0)  CYCLE
          Sum4inA = 0.E0
          DO jta=NoA+5, NoA
            IF (IndBackMRia(jta) == 0)  CYCLE
            Sum4inA = Sum4inA + BBcatanStr_Li(jtc,jta)*ca(IndBackMRia(jta))*ca(IndBackMRic(jtc))
          ENDDO
          SumIon4A=SumIon4A+Sum4inA
        ENDDO   
      
!------------------------------------------------------------------------
! Set LnGammaMR of Ions
         LnGammaIon = (SumIon1A/MMs) + 0.5*((Charge(j)**2.E0)/MMs)*SumIon2A + SumIon3A + 0.5*(Charge(j)**2.E0)*SumIon4A + SumIon5A        
         gamma_MR(j)=EXP(LnGammaIon)
      END DO
!++++++++++++++++++++++++++++++ End of Anions Loop  +++++++++++++++++++++++++++


!-----------------------------------------------------------------------------
 END SUBROUTINE NeuMR
SUBROUTINE SetLR(ca,Temp,Pres,gamma_LR)
!
!==========================================================
!===  Calculation of LongRange  Activitycoefficients 
!===   According to the AIMOFAC,A.Zuend et al(2008)
!==========================================================
!
!---  external variables
      REAL(RealKind) :: ca(:)             ! Molalities [mol/l]
      REAL(RealKind) :: Temp,Pres
      REAL(RealKind) :: gamma_LR(:)

      REAL(RealKind) :: Molal(SIZE(ca))

!---  internal variables
      INTEGER :: iSpc, ia, ie
      REAL(RealKind) :: Pre,tem,aw_DH,t

      REAL(RealKind) :: gam_DHLR(SIZE(ca))
      REAL(RealKind) :: SumMol
! ----------------------------------------------------------------------
!
         Tem = temp                     ! temperature [K]
         Pre = pres * 1013.e0           ! pressure    [at] ==>  [hPa]

! ----------------------------------------------------------------------
!--- determination of Debye-Hueckel term 
            gam_DHLR(:) = 1.e0        
            Molal(:)  = ca(:)     
            CALL Debye(Tem,Pre,Molal,gam_DHLR,aw_DH)

!--- set total coefficients 
            DO iSpc=1,SIZE(gamma_LR)
               gamma_LR(iSpc) = gamma_LR(iSpc) * gam_DHLR(iSpc)
            END DO

!---  set water activity 
            gamma_LR(iWater) = aw_DH  

     
END SUBROUTINE SetLR


! =========================================================================
! ===    Debye-Hueckel parameters
! =========================================================================

FUNCTION ADH(Temp)
  REAL(RealKind) :: ADH
  REAL(RealKind) :: Temp,p

  REAL(RealKind) ::  Pi,pp,b,c1,eps

  REAL(RealKind), PARAMETER :: u1=3.4279d2, u2=-5.0866d-3, u3=9.46900d-7  &
                              ,u4=-2.0525d0, u5=3.1159d3, u6=-1.8289d+2, u7=-8.0325d+3  &
                              ,u8=4.21142d6, u9=2.1417d0
  REAL(RealKind), PARAMETER :: t0 = 273.15d0
  REAL(RealKind), PARAMETER :: densw = 1.d0,            & ! water density [g/cm**3]
                               avo   = 6.02257d23,    & ! avogadro's number [molec./mol]
                               bol   = 1.3804d-16,    & ! bolzmann's constant [erg/K]
                               e     = 4.80298d-10      ! electronic charge [esu]
  Pi=4.d0*ATAN(1.d0)
  pp=p*1.d-3                                   ! [mbar] -> [bar]
  b=u7+u8/Temp+u9*Temp
  c1=u4+u5/(u6+Temp)
  eps=u1*EXP(u2*Temp+u3*Temp*Temp)
  eps=eps+c1*LOG((b+pp)/(b+1.d3))

  ADH=1.d0/3.d0*SQRT(2.0d0*pi*avo*densw*1.d-3)*(e*e/eps/bol/Temp)**(1.5d0)
END FUNCTION ADH

! =========================================================================
! FUNCTION ADH(Temp)
 
 ! APhiF calculates the Debye-Hueckel parameter aPhi
 
!     REAL(RealKind) :: ADH  
!     REAL(RealKind), PARAMETER :: densw = 1.E3               !995.010E0               ! Water density(kg/m**3)
!!    REAL(RealKind), PARAMETER :: epsw  = 80.1E0             !78.24787E0             ! static permitivity of water(C**2/J m)
!     REAL(RealKind), PARAMETER :: epsw  = 78.24787E0     
 
!     ADH = 1.327757E5 *SQRT(densw)/((epsw *Temp)**(3.E0/2.E0))
   
! END FUNCTION ADH
! =========================================================================
FUNCTION bDH(Temp)

! bDH calculates the Debye-Hueckel parameter b

  REAL(RealKind) :: Temp,p
  REAL(RealKind) :: bDH  
  
  REAL(RealKind), PARAMETER :: densw = 1.E3               !995.010E0               ! Water density(kg/m**3)
! REAL(RealKind), PARAMETER :: epsw  = 80.1E0             !78.24787E0             ! static permitivity of water(C**2/J m)
  REAL(RealKind), PARAMETER :: epsw  = 78.24787E0    

!   bDH = 6.359696E0 * SQRT(densw/(epsw*Temp))
  bDH = 1.2e0

END FUNCTION bDH

! =========================================================================
! ===  Subroutines for Debye-Hueckel term
! ===  Here water as a solvent
! =========================================================================
!


 SUBROUTINE Debye(Temp,p,ca,gamma1,aw_DH)
! ----------------------------------------------------------------------
!

!---  external variables
  REAL(RealKind) :: Temp,p,aw_DH, SumMolal, Osmatic,  Dummy
  REAL(RealKind) :: gamma1(:), ca(:)

!---  internal variables
  INTEGER :: iSpc
  REAL(RealKind) :: IonStr,SqrtIS
  REAL(RealKind) :: A, b                                         ! Dedye-Hueckel parameters
  REAL(RealKind) :: Term,Work, lngammak, SumM
  

!--------------------------------------
!---  Debye-Hueckel parameter
  A   = ADH(Temp)
  b   = bDH(Temp)

  
! ---  Ionic strength
  IonStr = IonicStrength(ca)
  SqrtIS = SQRT(IonStr)
!
!----------------------------------------------------------------------
!---  Debye-Hueckel term
     Work = 1.0E0+b*SqrtIS
     Term = -A *SqrtIS / Work           

!---  Set activity coefficients
  Gamma1(:) = 1.e0
  DO iSpc=1,SIZE(Gamma1)
     IF (Charge(iSpc) == 0.e0)  CYCLE
     Gamma1(iSpc) = (Charge(iSpc)*Charge(iSpc))*Term
     Gamma1(iSpc)= exp(Gamma1(iSpc))
  END DO
!---  Set water activity coefficient

  lngammak = (2.E0 * A *MM_H2O/(b**3.E0))*(Work-(1.E0/Work)-(2.E0*LOG(Work)))
  aw_DH    = EXP(lngammak)

END SUBROUTINE Debye


FUNCTION IonicStrength(yAqua)
!
!---  Computation of Dry Mass in a corresponding fraction
!
  REAL(RealKind) :: IonicStrength
  REAL(RealKind) :: yAqua(:)
  INTEGER :: jt

  IonicStrength=0.d0
  DO jt=1,nAqua
    IonicStrength=IonicStrength+Charge(jt)*Charge(jt)*yAqua(jt) 
  END DO
  IonicStrength=0.5d0*IonicStrength
END FUNCTION IonicStrength

! OSWALD ab hier alt



SUBROUTINE ComputeActivityAimo(cAqua,ActCoeff,TAbs)

  TYPE(Vec4_T), TARGET :: cAqua(0:),ActCoeff(0:),TAbs

  INTEGER :: i,iFrac,ix,iy,iz,is
! REAL(RealKind) :: xCat(nC),xAn(nA),xNeut(nN)
! REAL(RealKind) :: ActNeut(nN),ActCat(nC),ActAn(nA)

  REAL(RealKind) :: TempLoc,PresLoc
! REAL(RealKind) :: MolGes,MolWater,Fak,SumX
  REAL(RealKind) :: m(nAqua)
  REAL(RealKind) :: gamma_LR(nAqua)
  REAL(RealKind) :: gamma_Uni(nAqua)   ! molal activity coefficients
  REAL(RealKind) :: gamma_MR(nAqua)   ! molal activity coefficients
  REAL(RealKind) :: Mal2Mfr


! ActCoeff(:,:,:) =  gamma_LR(:,:,:) * gamma_Uni(:,:,:) * gamma_MR(:,:,:)
! DO jt=1,ntAqua
!   IF (ActIndex(jt) <= 0.E0) CYCLE
!   ActCoeff(ic,jt,ia:ie) = ActCoeff(ic,jt,ia:ie) / Mal2Mfr(ic,ia:ie)
! END DO

  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        TempLoc=TAbs%c(ix,iy,iz,1) 
        PresLoc=1.0d2
        DO iFrac=1,nFrac
          DO i=1,nAqua
            m(i)=cAqua(i)%c(ix,iy,iz,iFrac)/cAqua(iWater)%c(ix,iy,iz,iFrac)
          END DO   
          Mal2Mfr=MolH2O                          ! total mols
          DO i=iWater+1,nSoluble
            Mal2Mfr=Mal2Mfr+m(i)
          END DO  
          Mal2Mfr=Mal2Mfr/MolH2O          ! transformation factor
          gamma_LR=1.0d0
          gamma_Uni=1.0d0
          gamma_MR=1.0d0
          CALL SetLR(m,TempLoc,PresLoc,gamma_LR)
          CALL SetUnifac(m,TempLoc,gamma_UNI)
          CALL NeuMR(m,gamma_MR)
          DO jt=1,nAqua
            ActCoeff(jt)%c(ix,iy,iz,iFrac)=gamma_LR(jt)*gamma_Uni(jt)*gamma_MR(jt)
            IF (ActIndex(jt) <= 0.E0) CYCLE
            ActCoeff(jt)%c(ix,iy,iz,iFrac)=ActCoeff(jt)%c(ix,iy,iz,iFrac)/Mal2Mfr
          END DO  
!   ActCoeff(ic,jt,ia:ie) = ActCoeff(ic,jt,ia:ie) / Mal2Mfr(ic,ia:ie)
! END DO

!       CALL SetCoefficient_25(TempLoc)
!       DO iFrac=1,nFrac
!         IF (cAqua(iNC)%c(ix,iy,iz,iFrac)>Zero) THEN
!           IF (cAqua(iWater)%c(ix,iy,iz,iFrac)<Zero) THEN
!              WRITE(*,*) 'Negativ',cAqua(iWater)%c(ix,iy,iz,iFrac),iFrac,ix,iy,iz
!              STOP
!           END IF
!           SumX=Zero
!           xCat=Zero
!           xAn=Zero
!           xNeut=Zero
!           DO i=1,nAqua
!             IF (AimInd(i)>0) THEN
!               xCat(AimInd(i))=cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)
!               SumX=SumX+xCat(AimInd(i))
!             ELSE IF (AimInd(i)<0) THEN
!               xAn(ABS(AimInd(i)))=cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)
!               SumX=SumX+xAn(ABS(AimInd(i)))
!             END IF
!           END DO
!           xNeut(1)=cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(iWater)
!           SumX=SumX+xNeut(1)
!           xCat=xCat/SumX
!           xAn=xAn/SumX
!           xNeut=xNeut/SumX
!           CALL ExcessGibbs(xCAT,xAN,xNEUT,TempLoc &
!                           ,ActCat,ActAn,ActNeut &
!                           ,Excess)
!           MolWater=cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(iWater)
!           MolGes=0.0d0
!           DO i=iWater+1,nSoluble
!             MolGes=MolGes &
!                   +ABS(cAqua(i)%c(ix,iy,iz,iFrac))/MolMass(i)
!           END DO
!           Fak=One/(MolGes/MolWater+One)
!           DO i=1,nAqua
!             IF (AimInd(i)>0) THEN
!               ActCoeff(i)%c(ix,iy,iz,iFrac)=ActCat(AimInd(i))*Fak
!             ELSE IF (AimInd(i)<0) THEN
!               ActCoeff(i)%c(ix,iy,iz,iFrac)=ActAn(ABS(AimInd(i)))*Fak
!             END IF
!           END DO
!           ActCoeff(iWater)%c(ix,iy,iz,iFrac)=ActNeut(1)
!         END IF
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE ComputeActivityAimo

SUBROUTINE OutputActivityAimo(cAqua,ActCoeff,TAbs)

  TYPE(Vec4_T), TARGET :: cAqua(0:),ActCoeff(0:),TAbs

! INTEGER :: i,iFrac,ix,iy,iz,is
! REAL(RealKind) :: xCat(nC),xAn(nA),xNeut(nN)
! REAL(RealKind) :: ActNeut(nN),ActCat(nC),ActAn(nA)

! REAL(RealKind) :: TempLoc,Excess
! REAL(RealKind) :: Fak,SumX

! DO ix=ix0+1,ix1
!   DO iy=iy0+1,iy1
!     DO iz=iz0+1,iz1
!       TempLoc=TAbs%c(ix,iy,iz,1) ! OSSI Umrechnen
!       CALL SetCoefficient_25(TempLoc)
!       DO iFrac=1,nFrac
!         IF (cAqua(iNC)%c(ix,iy,iz,iFrac)>Zero) THEN
!           SumX=Zero
!           xCat=Zero
!           xAn=Zero
!           xNeut=Zero
!           DO i=1,nAqua
!             IF (AimInd(i)>0) THEN
!               xCat(AimInd(i))=cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)
!               SumX=SumX+xCat(AimInd(i))
!             ELSE IF (AimInd(i)<0) THEN
!               xAn(ABS(AimInd(i)))=cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)
!               SumX=SumX+xAn(ABS(AimInd(i)))
!             END IF
!           END DO
!           xNeut(1)=cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(iWater)
!           SumX=SumX+xNeut(1)
!           xCat=xCat/SumX
!           xAn=xAn/SumX
!           xNeut=xNeut/SumX
!           CALL ExcessGibbs(xCAT,xAN,xNEUT,TempLoc &
!                           ,ActCat,ActAn,ActNeut &
!                           ,Excess)
!           DO i=3,nAqua
!             IF (AimInd(i)>0) THEN
!               WRITE(*,*) SpeciesName(i),xCat(AimInd(i)),ActCat(AimInd(i))
!             ELSE IF (AimInd(i)<0) THEN
!               WRITE(*,*) SpeciesName(i),xAn(-AimInd(i)),ActAn(-AimInd(i))
!             END IF
!           END DO
!           WRITE(*,*) xNeut,ActNeut
!           DO i=3,nAqua
!             WRITE(*,*) SpeciesName(i),cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)/cAqua(iWater)%c(ix,iy,iz,iFrac), &
!                        ActCoeff(i)%c(ix,iy,iz,iFrac)
!           END DO
!         END IF
!       END DO
!     END DO
!   END DO
! END DO

END SUBROUTINE OutputActivityAimo

   SUBROUTINE SetUnifac(ca,Temp,gamma_Uni)

! =========================================================================
! ===  Determine UNIFAC activity coefficients
! =========================================================================
!
!---  external variables
      REAL(RealKind) :: Temp,ca(nAqua)
      REAL(RealKind) :: gamma_Uni(nAqua)

!---  internal variables
      INTEGER :: iSpc, jt, No

      REAL(RealKind) :: SumMol
!
! ------------------------------------------------------------------------
!
!

!--- compute activity coefficients of ideal solution for normalizing
            MolFractions(:) = 0.D0
            MolFractions(WaterIndUni) = 1.D0

            CALL Reference
!--- compute mole fraction activity coefficients
            CALL CombinatorialPart
            CALL ResidualPart
            CALL CalcResults
            f_Inf(:) = 1.E0
            f_Inf(:) = EXP(Results(:))
!RW         f_Inf(NoOrgs+1:NoMols) = EXP(Results(NoOrgs+1:NoMols))
!RW         f_Inf(:) = EXP(Results_Ref(:))      ! only temperature dependent
!------------------------------------------------------------------------
!--- compute mol fraction activity coefficients
      !--- transformation of molalities to mol fractions
            SumMol = MolH2O
            DO No=1,NoMols
               IF (UnifacBack(No) <= 0)  CYCLE
               SumMol = SumMol + ca(UnifacBack(No))
            END DO

            DO No=1,NoMols
               IF (UnifacBack(No) <= 0)  CYCLE
               MolFractions(No) = ca(UnifacBack(No)) / SumMol 
            END DO
            MolFractions(WaterIndUni) = MolH2O / SumMol
            CALL NormalizeFractions
       
      !--- compute mole fraction activity coefficients
            CALL CombinatorialPart           
            CALL ResidualPart
            CALL CalcResults              
!------------------------------------------------------------------------
!--- normalize activity coefficients
            DO jt=1,NoMols
               f_Star(:) = EXP(Results(:))/f_Inf(:)             ! Ions
            END DO

!--- water activity coefficients
!           aw_Uni(iFrac) = f_Star(IndSubWater)

!--- back-setting of activity coefficients
            DO jt=1,nAqua
              IF (UnifacInd(jt) <= 0)  CYCLE
              gamma_Uni(jt) = f_Star(UnifacInd(jt)) 
            END DO
            gamma_Uni(iWater) = f_Star(WaterIndUni)

   CONTAINS
!==========================================================================
!==========================================================================
!===  Calc.f90: calulation routines
!==========================================================================
!
  SUBROUTINE CombinatorialPart()
!
!-------------------------------------------
! calcs combinatorial part 
!-------------------------------------------
     INTEGER :: smd,group,j
     REAL(RealKind) :: Hr
     LOGICAL :: flag   
    
      
!
     DO smd=1,NoMols
       !-- calc sums
       rU(smd) = 0.0
       q(smd) = 0.0     
       DO group=SubrRowPtr(smd),SubrRowPtr(smd+1)-1
          rU(smd) = rU(smd) + SubrSparse(group) * rkU(SubrColInd(group))
          q(smd) = q(smd) + SubrSparse(group) * QK(SubrColInd(group))         
       END DO
       !-- l
       l(smd) = (5.E0 * (rU(smd) - q(smd))) - (rU(smd) - 1.E0) 
     END DO ! smd


     !-- add up SumRX,QX
     SumRX = 0.0
     SumQX = 0.0
     SumLX = 0.0
     DO smd=1,NoMols
       SumRX = SumRX + rU(smd) * MolFractions(smd)
       SumQX = SumQX + q(smd) * MolFractions(smd)
       SumLX = SumLX + l(smd) * MolFractions(smd)
     END DO

!-------------------------------------------------------------------
     !-- calc phi's and xi's 
     DO smd = 1, NoMols
       phi(smd) = rU(smd)  / SumRX
       xi (smd) = q(smd) / SumQX
     END DO
     
     
     DO smd=1,NoMols
       Hr = xi(smd) / phi(smd)
       ln_gamma_c(smd)      = log(phi(smd)) +  5.0*q(smd)*log(Hr) +  l(smd) - phi(smd)*SumLX            
     END DO    
!********************************************************************
     
   END SUBROUTINE CombinatorialPart 

!==========================================================================
!===  Calculation of Reference state
!==========================================================================
!
  SUBROUTINE Reference()
!
!-------------------------------------------
! calcs combinatorial part
!-------------------------------------------
     INTEGER :: smd,group,j
     REAL(RealKind) :: Hr
     MolFractions(:) = 0.D0
     MolFractions(WaterIndUni) = 1.D0

!--------------------------------------------------
! Compute the parameters for reference solution
!--------------------------------------------------
!      rw(IndSubWater) =  RK(IndSubWater)
!      qw(IndSubWater) =  QK(IndSubWater)

     rU(:) = 0.E0
     q(:) = 0.E0   

!---  r and q for H2O
     DO group=SubrRowPtr(WaterIndUni),SubrRowPtr(WaterIndUni+1)-1
        rU(WaterIndUni) = rU(WaterIndUni) + SubrSparse(group) * rkU(SubrColInd(group))
        q(WaterIndUni) = q(WaterIndUni) + SubrSparse(group) * QK(SubrColInd(group))         
     END DO

!---  r and q for ions
     DO smd=NoOrgs+1,NoMols
       !-- calc sums
       DO group=SubrRowPtr(smd),SubrRowPtr(smd+1)-1
          rU(smd) = rU(smd) + SubrSparse(group) * rkU(SubrColInd(group))
          q(smd) = q(smd) + SubrSparse(group) * QK(SubrColInd(group))         
       END DO
     END DO ! smd

!---  Calculate for combinatorical part
     DO smd=NoOrgs+1,NoMols
       ln_gamma_c_ref(smd)  = log((rU(smd)/rU(WaterIndUni)))+ 1.E0  - (rU(smd)/rU(WaterIndUni)) + 5.E0 *q(smd)*     &
&                        ((log(rU(WaterIndUni)*q(smd))/(rU(smd)*q(WaterIndUni)))-1.E0 + ((rU(smd)*q(WaterIndUni))/&
&(rU(WaterIndUni)*q(smd)))) 
     END DO

!---  Calculate for residual part   (Not needed due to the fact that ion-solvent interactions be zero!!))

!---  Set reference state
     Results_ref(:) = 0.E0
     DO smd=NoOrgs+1,NoMols
        Results_ref(smd) = ln_gamma_c_ref(smd) ! FOR IONS (Reference state calculation)
     END DO
     
!********************************************************************

 END SUBROUTINE Reference 
!===========================================
   SUBROUTINE ResidualPart()
!
!-------------------------------------------
!---  calcs residual parts
!-------------------------------------------
!
     INTEGER :: smd,i,j,group,group2,group3, SMDOccurrences(NoMols),  &
  &             SumSubOccurrences(NoSubGroups)
     REAL(RealKind) :: Denominator, Zaehler, Denominators(NoMols), Hr, Hr2,  &
  &             Hr3,Hr4,ExpAt

!-------------------------------------------
! --- count occurrences
     DO group = 1,NoSubGroups
       i = 0
       DO smd=SubcColPtr(group),SubcColPtr(group+1)-1
          i = i + SubcSparse(smd)
       END DO
       SumSubOccurrences(group) = i
     END DO  ! groups

!----------------------------------------------------------------------
! --- calc occurrences, denominator of X-expression
     Denominator = 0.0
     DO smd=1,NoMols
       SMDOccurrences(smd) = 0
       DO group=SubrRowPtr(smd),SubrRowPtr(smd+1)-1
          SMDOccurrences(smd) = SMDOccurrences(smd) + SubrSparse(group)               
       END DO ! groups
       Denominator = Denominator + SMDOccurrences(smd) * MolFractions(smd)
     END DO ! smd=1,NoMols

!----------------------------------------------------------------------
! --- calc X
     x(:,:) = 0.D0
     DO group = 1,NoSubGroups
     !--- upper
       Hr = 0.0
       DO smd=SubcColPtr(group),SubcColPtr(group+1)-1
         Hr  = Hr + MolFractions(SubcRowInd(smd)) * SubcSparse(smd)
     !--- convert to reaL
         Hr3 = SubcSparse(smd)
         Hr4 = SMDOccurrences(SubcRowInd(smd))
         x(SubcRowInd(smd),group) = Hr3 / Hr4
       END DO ! mols
       TotalX(group) = Hr / Denominator
       Theta2(group) = Hr
     END DO ! groups

     !-- calc theta, denominators
     DO smd=1,NoMols
       Denominators(smd) = 0.0
       DO group = 1,NoSubGroups
          Denominators(smd) = Denominators(smd) + qk(group) * x(smd,group)
       END DO ! group
     END DO ! mol

     !-- total denominator
     Denominator = 0.0
     DO group = 1,NoSubGroups
       Denominator = Denominator + qk(group) * TotalX(group)    
     END DO ! group

     !--- calc theta
     DO group = 1,NoSubGroups
       DO smd=1,NoMols
          theta(smd,group) = qk(group) * x(smd,group) / Denominators(smd)
       END DO ! mols
       TotalTheta(group) = qk(group) * TotalX(group) / Denominator
     END DO         ! groups

!----------------------------------------------------------------------
     DO group = 1,NoSubGroups  ! k-loop
     ! -- calc gamma k
       Hr  = 0.0
       Hr2 = 0.0
       TotalGamma(group) = 0
       IF (SumSubOccurrences(group) /= 0) THEN
         DO group2 = 1,NoSubGroups ! m-loop
           IF (SumSubOccurrences(group2) /= 0) THEN  
             !-- sum up den for 3rd term
             Denominator = 0.0
             DO group3 = 1,NoSubGroups ! n-loop
               IF (SumSubOccurrences(group3) /= 0) THEN
             !-- get interaction value of main group
                 ExpAt = GetInteraction(group3,group2,.FALSE.)
                 Denominator = Denominator + TotalTheta(group3) * ExpAT
               END IF ! (SumSubOccurrences(group3) /= 0)  
             END DO  ! group3

             !-- sum up  3rd term
             ExpAT = GetInteraction(group,group2,.TRUE.)
             Hr = Hr + TotalTheta(group2) * ExpAt / Denominator
             !-- sum up  2nd term
             ExpAT = GetInteraction(group2,group,.TRUE.)
             Hr2 = Hr2 + TotalTheta(group2) * ExpAt
           END IF    ! SumSubOccurrences(group2) /= 0  
         END DO ! group2

         !-- sum up
         TotalGamma(group) = qk(group) * (1.0 - log(Hr2) - Hr)
         nue(group) = 0.0
         DO smd=1,NoMols
           s(group,smd) = 0.0
           DO group2=SubrRowPtr(smd),SubrRowPtr(smd+1)-1
              s(group,smd) = s(group,smd) + SubrSparse(group2) * QK(SubrColInd(group2)) *  &
  &                          GetInteraction(SubrColInd(group2),group,.TRUE.)
           END DO      ! group2

           nue(group) = nue(group) + s(group,smd) * MolFractions(smd)
         END DO        ! smd

       END IF        ! occurrences there
     END DO  ! group

!----------------------------------------------------------------------
!---- calc smd Gammas
     GammaU(:,:) = 0
     DO smd=1,NoMols          
        DO group=SubrRowPtr(smd),SubrRowPtr(smd+1)-1
           !-- calc gamma k
           Hr  = 0.0
           Hr2 = 0.0
           DO group2=SubrRowPtr(smd),SubrRowPtr(smd+1)-1
              !-- sum up den for 3rd term
              Denominator = 0.0
              DO group3=SubrRowPtr(smd),SubrRowPtr(smd+1)-1
                 !-- get interaction value of main group
                 ExpAt = GetInteraction(SubrColInd(group3),SubrColInd(group2),.TRUE.)
                 Denominator = Denominator + Theta(smd,SubrColInd(group3)) * ExpAT
              END DO  ! group3
              !-- sum up  3rd term
              ExpAT = GetInteraction(SubrColInd(group),SubrColInd(group2),.TRUE.)
              Hr = Hr + Theta(smd,SubrColInd(group2)) * ExpAt / Denominator
              !-- sum up  2nd term
              ExpAT = GetInteraction(SubrColInd(group2),SubrColInd(group),.TRUE.)
              Hr2 = Hr2 + Theta(smd,SubrColInd(group2)) * ExpAt
           END DO ! group2

           !-- sum up
           IF (Hr2 <= 0) Hr = 2     ! interaction unknown error anyway,,
           GammaU(smd,SubrColInd(group)) = qk(SubrColInd(group)) * (1.0 - log(Hr2) - Hr)

        END DO ! group 
      
     END DO  ! mol

!----------------------------------------------------------------------
!---- calc residuals
     DO smd=1,NoMols
       Hr = 0.0
       DO group=SubrRowPtr(smd),SubrRowPtr(smd+1)-1
         Hr = Hr + SubrSparse(group) * (TotalGamma(SubrColInd(group)) - GammaU(smd,SubrColInd(group)))
       END DO ! 
       ln_gamma_r(smd) = Hr
       
     END DO ! mol
!      ln_gamma_r(:) =0.e0
   END SUBROUTINE ResidualPart
!
!===========================================
   REAL(RealKind) FUNCTION GetInteraction(i,j,Shutup)
!
!-------------------------------------------
!--- det Interaction value of subgroups
!-------------------------------------------
   INTEGER :: i,j,mi,mj,k1,k2
   LOGICAL :: Shutup
   REAL(RealKind) :: Value

!-- det main groups 
   mi = MainGroups(i)
   mj = MainGroups(j)
   Value = MatrixA(mi,mj)

!-- -999 matrix value?
   IF (IsErrorValue(Value)) THEN
      Value = 0.0
   END IF ! suppress error

   GetInteraction = exp(-Value / Temp)   

  END FUNCTION GetInteraction
!
!===========================================
  SUBROUTINE CalcResults
!
! -------------------------------------------
!---  calcs result values
! -------------------------------------------

    INTEGER :: smd, j
    REAL(RealKind) :: Hr        

    Exzess = 0.0
    DO smd=1,NoMols
        Results(smd) = ln_gamma_c(smd) + ln_gamma_r(smd) 
        Exzess = Exzess + MolFractions(smd) * Results(smd)
    END DO ! mols
    GExzess = Exzess * RGAS * Temp

  END SUBROUTINE CalcResults
!
!===========================================
  LOGICAL FUNCTION IsErrorValue(Value)
!
!----------------------------------------------------------
!--- Returns true if x = -999 
!-------------------------------------------
     REAL(RealKind) :: Value
     IF ( ABS(Value -  ERROR_VALUE) <  ERROR_BORDER) THEN
       IsErrorValue = .TRUE.
     ELSE
       IsErrorValue = .FALSE.
     END IF
  END FUNCTION IsErrorValue    ! LOGICAL FUNCTION IsErrorValue(x)
   END Subroutine SetUnifac
END MODULE ActivityAimo_Mod


