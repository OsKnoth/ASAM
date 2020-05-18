MODULE StackVector_Mod

  USE Kind_Mod
  USE DataType_Mod
  USE Variable_Mod

  IMPLICIT NONE

  TYPE StackVector_T
    TYPE(Vec4_T) :: c
    LOGICAL :: Temporary=.TRUE.
  END TYPE StackVector_T

  INTERFACE ASSIGNMENT(=)
    MODULE PROCEDURE Copy_StackVector,Copy_StackScalar
  END INTERFACE
  INTERFACE OPERATOR(+)
    MODULE PROCEDURE Add_StackVector
  END INTERFACE
  INTERFACE OPERATOR(-)
    MODULE PROCEDURE Sub_StackVector,Neg_StackVector
  END INTERFACE
  INTERFACE OPERATOR(*)
    MODULE PROCEDURE Mult_StackVector
  END INTERFACE
  INTERFACE OPERATOR(/)
    MODULE PROCEDURE Div_StackVector
  END INTERFACE
  INTERFACE EXP
    MODULE PROCEDURE Exp_StackVector
  END INTERFACE
  INTERFACE LOG
    MODULE PROCEDURE Log_StackVector
  END INTERFACE
  INTERFACE OPERATOR(**)
    MODULE PROCEDURE Pow_StackVector
  END INTERFACE
  INTERFACE Allocate
    MODULE PROCEDURE Allocate_StackVector
  END INTERFACE
  INTERFACE Deallocate
    MODULE PROCEDURE Deallocate_StackVector
  END INTERFACE
  TYPE(Vec4_T), POINTER, PRIVATE :: WorkLong(:)
  INTEGER :: iwLong=1
  TYPE(Vec4_T), POINTER, PRIVATE :: WorkShort(:)
  INTEGER :: iwShort=1

CONTAINS

SUBROUTINE Init_StackVector(Vec)

  TYPE(Vec4_T) :: Vec
  INTEGER :: i

  ALLOCATE(WorkLong(10))
  ALLOCATE(WorkShort(10))
  DO i=1,10
    ALLOCATE(WorkLong(i)%c(Size(Vec%cInt,1),Size(Vec%cInt,2),Size(Vec%cInt,3),Size(Vec%cInt,4)))
    ALLOCATE(WorkShort(i)%c(Size(Vec%cInt,1),Size(Vec%cInt,2),Size(Vec%cInt,3),1))
  END DO
  iwLong=1
  iwShort=1
    
END SUBROUTINE Init_StackVector

SUBROUTINE Close_StackVector

  INTEGER :: i

  DO i=1,10
    DEALLOCATE(WorkLong(i)%c)
    DEALLOCATE(WorkShort(i)%c)
  END DO
  DEALLOCATE(WorkLong)
  DEALLOCATE(WorkShort)

END SUBROUTINE Close_StackVector

SUBROUTINE Output_StackVector

  WRITE(*,*) 'iwLong',iwLong
  WRITE(*,*) 'iwShort',iwShort
  IF (iwLong>iwShort) THEN
    WRITE(*,*) WorkLong(iwLong-1)%c
  ELSE
    WRITE(*,*) WorkLong(iwLong-1)%c
  END IF

END SUBROUTINE Output_StackVector

SUBROUTINE Allocate_StackVector(Vec,Dim)

  TYPE(StackVector_T) :: Vec
  INTEGER :: Dim

  Vec%Temporary=.FALSE.
  ALLOCATE(Vec%c%c(1,1,1,Dim))

END SUBROUTINE Allocate_StackVector

SUBROUTINE Deallocate_StackVector(Vec)

  TYPE(StackVector_T) :: Vec

  DEALLOCATE(Vec%c%c)

END SUBROUTINE Deallocate_StackVector


FUNCTION Add_StackVector(Vec1,Vec2)

  TYPE(StackVector_T) :: Add_StackVector
  TYPE(StackVector_T), INTENT(IN) :: Vec1,Vec2

  INTEGER :: i,i1,i2,Len1,Len2,MaxLen,Shift1,Shift2

! Add=Vec1+Vec2

  Len1=SIZE(Vec1%c%c,4)
  Shift1=MIN(1,Len1-1)
  i1=1
  Len2=SIZE(Vec2%c%c,4)
  Shift2=MIN(1,Len2-1)
  i2=1
  
  IF (Vec1%Temporary) THEN
    iwShort=iwShort+Shift1-1
    iwLong=iwLong-Shift1
  END IF
  IF (Vec2%Temporary) THEN
    iwShort=iwShort+Shift2-1
    iwLong=iwLong-Shift2
  END IF
  Maxlen=MAX(Len1,Len2)
  IF (MaxLen==1) THEN
    Add_StackVector%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
  ELSE
    Add_StackVector%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
  END IF
  DO i=1,Maxlen
    Add_StackVector%c%c(:,:,:,i)=Vec1%c%c(:,:,:,i1)+Vec2%c%c(:,:,:,i2)
    i1=i1+Shift1
    i2=i2+Shift2
  END DO
      
END FUNCTION Add_StackVector
FUNCTION Sub_StackVector(Vec1,Vec2)

  TYPE(StackVector_T) :: Sub_StackVector
  TYPE(StackVector_T), INTENT(IN) :: Vec1,Vec2

  INTEGER :: i,i1,i2,Len1,Len2,MaxLen,Shift1,Shift2

! Sub=Vec1-Vec2

  Len1=SIZE(Vec1%c%c,4)
  Shift1=MIN(1,Len1-1)
  i1=1
  Len2=SIZE(Vec2%c%c,4)
  Shift2=MIN(1,Len2-1)
  i2=1

  IF (Vec1%Temporary) THEN
    iwShort=iwShort+Shift1-1
    iwLong=iwLong-Shift1
  END IF
  IF (Vec2%Temporary) THEN
    iwShort=iwShort+Shift2-1
    iwLong=iwLong-Shift2
  END IF
  Maxlen=MAX(Len1,Len2)
  IF (MaxLen==1) THEN
    Sub_StackVector%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
  ELSE
    Sub_StackVector%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
  END IF
  DO i=1,Maxlen
    Sub_StackVector%c%c(:,:,:,i)=Vec1%c%c(:,:,:,i1)-Vec2%c%c(:,:,:,i2)
    i1=i1+Shift1
    i2=i2+Shift2
  END DO

END FUNCTION Sub_StackVector

FUNCTION Mult_StackVector(Vec1,Vec2)

  TYPE(StackVector_T) :: Mult_StackVector
  TYPE(StackVector_T), INTENT(IN) :: Vec1,Vec2

  INTEGER :: i,i1,i2,Len1,Len2,MaxLen,Shift1,Shift2

! Mult=Vec1*Vec2

  Len1=SIZE(Vec1%c%c,4)
  Shift1=MIN(1,Len1-1)
  i1=1
  Len2=SIZE(Vec2%c%c,4)
  Shift2=MIN(1,Len2-1)
  i2=1

  IF (Vec1%Temporary) THEN
    iwShort=iwShort+Shift1-1
    iwLong=iwLong-Shift1
  END IF
  IF (Vec2%Temporary) THEN
    iwShort=iwShort+Shift2-1
    iwLong=iwLong-Shift2
  END IF
  Maxlen=MAX(Len1,Len2)
  IF (MaxLen==1) THEN
    Mult_StackVector%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
  ELSE
    Mult_StackVector%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
  END IF
  DO i=1,Maxlen
    Mult_StackVector%c%c(:,:,:,i)=Vec1%c%c(:,:,:,i1)*Vec2%c%c(:,:,:,i2)
    i1=i1+Shift1
    i2=i2+Shift2
  END DO

END FUNCTION Mult_StackVector

FUNCTION Div_StackVector(Vec1,Vec2)

  TYPE(StackVector_T) :: Div_StackVector
  TYPE(StackVector_T), INTENT(IN) :: Vec1,Vec2

  INTEGER :: i,i1,i2,Len1,Len2,MaxLen,Shift1,Shift2

! Div=Vec1/Vec2

  Len1=SIZE(Vec1%c%c,4)
  Shift1=MIN(1,Len1-1)
  i1=1
  Len2=SIZE(Vec2%c%c,4)
  Shift2=MIN(1,Len2-1)
  i2=1

  IF (Vec1%Temporary) THEN
    iwShort=iwShort+Shift1-1
    iwLong=iwLong-Shift1
  END IF
  IF (Vec2%Temporary) THEN
    iwShort=iwShort+Shift2-1
    iwLong=iwLong-Shift2
  END IF
  Maxlen=MAX(Len1,Len2)
  IF (MaxLen==1) THEN
    Div_StackVector%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
  ELSE
    Div_StackVector%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
  END IF
  DO i=1,Maxlen
    WHERE (Vec2%c%c(:,:,:,i2)/=Zero)
      Div_StackVector%c%c(:,:,:,i)=Vec1%c%c(:,:,:,i1)/Vec2%c%c(:,:,:,i2)
    ELSEWHERE
      Div_StackVector%c%c(:,:,:,i)=Zero
    END WHERE
    i1=i1+Shift1
    i2=i2+Shift2
  END DO
END FUNCTION Div_StackVector

FUNCTION Pow_StackVector(Vec1,Vec2)

  TYPE(StackVector_T) :: Pow_StackVector
  TYPE(StackVector_T), INTENT(IN) :: Vec1,Vec2

  INTEGER :: i,i1,i2,Len1,Len2,MaxLen,Shift1,Shift2

! Pow=Vec1/Vec2

  Len1=SIZE(Vec1%c%c,4)
  Shift1=MIN(1,Len1-1)
  i1=1
  Len2=SIZE(Vec2%c%c,4)
  Shift2=MIN(1,Len2-1)
  i2=1

  IF (Vec1%Temporary) THEN
    iwShort=iwShort+Shift1-1
    iwLong=iwLong-Shift1
  END IF
  IF (Vec2%Temporary) THEN
    iwShort=iwShort+Shift2-1
    iwLong=iwLong-Shift2
  END IF
  Maxlen=MAX(Len1,Len2)
  IF (MaxLen==1) THEN
    Pow_StackVector%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
  ELSE
    Pow_StackVector%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
  END IF
  DO i=1,Maxlen
    Pow_StackVector%c%c(:,:,:,i)=Vec1%c%c(:,:,:,i1)**Vec2%c%c(:,:,:,i2)
    i1=i1+Shift1
    i2=i2+Shift2
  END DO
END FUNCTION Pow_StackVector

FUNCTION Neg_StackVector(Vec1)

  TYPE(StackVector_T) :: Neg_StackVector
  TYPE(StackVector_T), INTENT(IN) :: Vec1

  REAL(RealKind) :: Scalar

! Neg=-Vec1

  IF (Vec1%Temporary) THEN
    IF (SIZE(Vec1%c%c,4)==1) THEN
      iwShort=iwShort-1
    ELSE
      iwLong=iwLong-1
    END IF
  END IF
  IF (SIZE(Vec1%c%c,4)==1) THEN
    Neg_StackVector%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
    Neg_StackVector%c%c=-Vec1%c%c
  ELSE
    Neg_StackVector%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
    Neg_StackVector%c%c=-Vec1%c%c
  END IF
END FUNCTION Neg_StackVector

FUNCTION Exp_StackVector(Vec1)

  TYPE(StackVector_T) :: Exp_StackVector
  TYPE(StackVector_T), INTENT(IN) :: Vec1

  REAL(RealKind) :: Scalar

! Exp=EXP(Vec1)

  IF (Vec1%Temporary) THEN
    IF (SIZE(Vec1%c%c,4)==1) THEN
      iwShort=iwShort-1
    ELSE
      iwLong=iwLong-1
    END IF
  END IF
  IF (SIZE(Vec1%c%c,4)==1) THEN
    Exp_StackVector%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
    Exp_StackVector%c%c=EXP(Vec1%c%c)
  ELSE
    Exp_StackVector%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
    Exp_StackVector%c%c=EXP(Vec1%c%c)
  END IF
END FUNCTION Exp_StackVector

FUNCTION Log_StackVector(Vec1)

  TYPE(StackVector_T) :: Log_StackVector
  TYPE(StackVector_T), INTENT(IN) :: Vec1

  REAL(RealKind) :: Scalar

! Log=LOG(Vec1)

  IF (Vec1%Temporary) THEN
    IF (SIZE(Vec1%c%c,4)==1) THEN
      iwShort=iwShort-1
    ELSE
      iwLong=iwLong-1
    END IF
  END IF
  IF (SIZE(Vec1%c%c,4)==1) THEN
    Log_StackVector%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
    Log_StackVector%c%c=LOG(Vec1%c%c)
  ELSE
    Log_StackVector%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
    Log_StackVector%c%c=LOG(Vec1%c%c)
  END IF
END FUNCTION Log_StackVector

SUBROUTINE Copy_StackScalar(Vec1,Scalar)
 
  TYPE(StackVector_T), INTENT(INOUT) :: Vec1
  TYPE(Vec4_T), INTENT(IN) :: Scalar

  IF (SIZE(Scalar%c,4)>1) THEN
    Vec1%c%c=>WorkLong(iwLong)%c
    iwLong=iwLong+1
  ELSE
    Vec1%c%c=>WorkShort(iwShort)%c
    iwShort=iwShort+1
  END IF
  Vec1%c%c=Scalar%c

END SUBROUTINE Copy_StackScalar 

SUBROUTINE Copy_StackVector(Vec1,Vec2)

  TYPE(StackVector_T), INTENT(INOUT) :: Vec1
  TYPE(StackVector_T), INTENT(IN) :: Vec2

  Vec1%c%c=>Vec2%c%c

END SUBROUTINE Copy_StackVector

END MODULE StackVector_Mod

MODULE Reaction_Mod
  USE Kind_Mod
  USE StackVector_Mod

  IMPLICIT NONE

  TYPE Reactant_T
    INTEGER :: species
    REAL(RealKind) :: Koeff
  END TYPE Reactant_T

  TYPE member_infix
    CHARACTER(len=5) :: op=''
    TYPE(Vec4_T) :: Value 
    CHARACTER(len=5) :: operand=''
    TYPE (member_infix), POINTER :: next=>NULL()
  END TYPE member_infix

  TYPE Reaction_T
    CHARACTER*15 :: ClassR='          '
    CHARACTER*15 :: TypeR='          '
    LOGICAL :: TypeE=.FALSE.
    CHARACTER(20), POINTER :: Factor=>NULL()
    INTEGER :: NumSpeciesLeft
    INTEGER :: NumSpeciesLeftAktiv
    INTEGER, POINTER :: SpeciesLeft(:)=>NULL()
    INTEGER, POINTER :: SpeciesLeftLoc(:)=>NULL() 
    INTEGER :: NumSpeciesRight
    INTEGER :: NumSpeciesRightAktiv
    INTEGER, POINTER :: SpeciesRight(:)=>NULL() 
    INTEGER, POINTER :: SpeciesRightLoc(:)=>NULL() 
    INTEGER :: NumSpecies
    TYPE(Reactant_T), POINTER :: Species(:)
    INTEGER :: NumConstants
    REAL(RealKind), POINTER :: Constants(:,:,:,:)=>NULL()
    INTEGER, POINTER :: struct(:) 
    TYPE (Reaction_T), POINTER :: next=>NULL()
    TYPE (member_infix), POINTER :: Infix=>NULL()
    REAL(RealKind) :: fac_exp
    REAL(RealKind) :: fac_A 
  END TYPE Reaction_T

! TYPE operator
!   CHARACTER(LEN=5) :: operator
!   INTEGER          :: rang
! END TYPE operator
! TYPE (operator), DIMENSION(10) :: stack_operator
! INTEGER p
  INTEGER, PARAMETER :: NumberOperanden=10
  CHARACTER(len=4), DIMENSION(10), PARAMETER :: name_operand= &
 (/'$c1 ','$c2 ','$c3 ','$c4 ','$c5 ','$c6 ','$c7 ','$c8 ','$c9 ','$c10'/)
! INTEGER, PARAMETER :: number_operanden=10
! CHARACTER(len=5), PARAMETER :: possible_operator='+-*^/'
! CHARACTER(len=2), PARAMETER :: bracket='()'
! CHARACTER(len=10), PARAMETER :: leer=' '
! LOGICAL, PARAMETER :: back=.TRUE.

  TYPE(StackVector_T), SAVE :: Stack_infix(20)
  TYPE(Vec4_T), POINTER :: Temperature
  TYPE(Vec4_T), POINTER :: cVecReac(:)
  INTEGER :: NumReactionTypes=14

CONTAINS

SUBROUTINE calculate(w,start_ausdruck,Constants)

!---------------------------------------------------------------
!---  Computation of Reaction Constants for
!---  Special Type Reactions 
!---------------------------------------------------------------

  TYPE(Vec4_T) :: w
  TYPE (member_infix), POINTER :: start_ausdruck
  REAL(RealKind) :: Constants(:,:,:,:)

  INTEGER :: i,pc
  TYPE (member_infix), POINTER :: current_infix

  CALL Replace(start_ausdruck,Constants)
    
  pc=0
  iwLong=1
  iwShort=1
  current_infix=>start_ausdruck
  DO while(associated(current_infix))
    SELECT CASE(current_infix%op)
    CASE('+')  
      stack_infix(pc-1)=stack_infix(pc-1)+stack_infix(pc)
      stack_infix(pc)%Temporary=.TRUE.
      stack_infix(pc-1)%Temporary=.TRUE.
      pc=pc-1
    CASE('-')  
      stack_infix(pc-1)=stack_infix(pc-1)-stack_infix(pc)
      stack_infix(pc)%Temporary=.TRUE.
      stack_infix(pc-1)%Temporary=.TRUE.
      pc=pc-1
    CASE('*')  
      stack_infix(pc-1)=stack_infix(pc-1)*stack_infix(pc)
      stack_infix(pc)%Temporary=.TRUE.
      stack_infix(pc-1)%Temporary=.TRUE.
      pc=pc-1
    CASE('/')  
      stack_infix(pc-1)=stack_infix(pc-1)/stack_infix(pc)
      stack_infix(pc)%Temporary=.TRUE.
      stack_infix(pc-1)%Temporary=.TRUE.
      pc=pc-1
    CASE('^')  
      stack_infix(pc-1)=stack_infix(pc-1)**stack_infix(pc)
      stack_infix(pc)%Temporary=.TRUE.
      stack_infix(pc-1)%Temporary=.TRUE.
      pc=pc-1
    CASE('--')  
      stack_infix(pc)=-stack_infix(pc)
      stack_infix(pc)%Temporary=.TRUE.
    CASE('exp')  
      stack_infix(pc)=EXP(stack_infix(pc))
      stack_infix(pc)%Temporary=.TRUE.
    CASE('log')  
      stack_infix(pc)=LOG(stack_infix(pc))
      stack_infix(pc)%Temporary=.TRUE.
    CASE DEFAULT
      pc=pc+1
      stack_infix(pc)=current_infix%value
      stack_infix(pc)%Temporary=.TRUE. !OSSI
    END SELECT 
    current_infix=>current_infix%next
  END DO
  IF (SIZE(w%cInt,4)==SIZE(stack_infix(1)%c%c,4)) THEN
    w%cInt=stack_infix(1)%c%c
  ELSE
    DO i=1,SIZE(w%cInt,4)
      w%cInt(:,:,:,i)=stack_infix(1)%c%c(:,:,:,1)
    END DO
  END IF
END SUBROUTINE calculate

SUBROUTINE Replace(start_ausdruck,Constants)

  TYPE (member_infix), POINTER :: start_ausdruck
  REAL(RealKind), TARGET :: Constants(:,:,:,:)

  INTEGER :: i
  TYPE (member_infix), POINTER :: current_infix

  current_infix=>start_ausdruck
  DO while(associated(current_infix))
    IF (current_infix%operand=='te') THEN 
      current_infix%value%c=>Temperature%cInt
    ELSE
      DO i=1,SIZE(SpeciesName)
        IF (current_infix%operand==SpeciesName(i)) THEN
          current_infix%value%c=>cVecReac(i)%cInt
          EXIT
        END IF
      END DO
      DO i=1,NumberOperanden
        IF (current_infix%operand==name_operand(i)) THEN
!         ALLOCATE(current_infix%value%c(1))
          current_infix%value%c=>Constants(:,:,:,i:i)
          EXIT
        END IF
      END DO
    END IF
    current_infix=>current_infix%next
  END DO
END SUBROUTINE replace

SUBROUTINE ReactionInput(Reaction,Unit)

  TYPE(Reaction_T) :: Reaction
  INTEGER :: Unit

  INTEGER :: i,l,is
  CHARACTER*200 :: string_aktuell
  TYPE (member_infix), POINTER :: current_infix
  TYPE operator
    CHARACTER(LEN=5) :: operator
    INTEGER          :: rang
  END TYPE operator
  TYPE (operator), DIMENSION(10) :: stack_operator
  INTEGER :: p
  LOGICAL :: Start
  CHARACTER(len=5), PARAMETER :: possible_operator='+-*^/'
  CHARACTER(len=2), PARAMETER :: bracket='()'
  CHARACTER(len=10), PARAMETER :: leer=' '
  LOGICAL, PARAMETER :: back=.TRUE.
  CHARACTER(90) :: chem_factor


  READ(Unit,*,END=100) Reaction%ClassR,Reaction%TypeR
  IF (INDEX(Reaction%TypeR,'_EXP')>0) THEN
    Reaction%TypeR=Reaction%TypeR(1:INDEX(Reaction%TypeR,'_EXP')-1)
    Reaction%TypeE=.TRUE.
  END IF  
  READ(Unit,*) Reaction%NumSpeciesLeft &
              ,Reaction%NumSpeciesRight &
              ,Reaction%NumSpeciesLeftAktiv &     
              ,Reaction%NumSpeciesRightAktiv     
  ALLOCATE(Reaction%SpeciesLeft(Reaction%NumSpeciesLeft))
  ALLOCATE(Reaction%SpeciesRight(Reaction%NumSpeciesRight))
  READ(Unit,*)  &                        
               (Reaction%SpeciesLeft(i),i=1,Reaction%NumSpeciesLeft) &
              ,(Reaction%SpeciesRight(i),i=1,Reaction%NumSpeciesRight) &
              ,Reaction%NumSpecies
  ALLOCATE(Reaction%Species(Reaction%NumSpecies))
  READ(Unit,*) (Reaction%Species(i),i=1,Reaction%NumSpecies)
  IF (Reaction%TypeE) THEN
    ALLOCATE(Reaction%SpeciesLeftLoc(Reaction%NumSpeciesLeft))
    ALLOCATE(Reaction%SpeciesRightLoc(Reaction%NumSpeciesRight))
    DO l=1,Reaction%NumSpeciesLeft
      is=Reaction%SpeciesLeft(l)
      DO i=1,Reaction%NumSpecies
        IF (is==Reaction%Species(i)%species) THEN
          Reaction%SpeciesLeftLoc(l)=i
          EXIT
        END IF  
      END DO  
    END DO  
    DO l=1,Reaction%NumSpeciesRight
      is=Reaction%SpeciesRight(l)
      DO i=1,Reaction%NumSpecies
        IF (is==Reaction%Species(i)%species) THEN
          Reaction%SpeciesRightLoc(l)=i
          EXIT
        END IF  
      END DO  
    END DO  
  END IF  
  IF (Reaction%TypeR.eq.'SPECIAL') THEN
    string_aktuell=' '
    READ(Unit,'(a200)') string_aktuell
  END IF
  READ(Unit,*) Reaction%NumConstants
!
  ALLOCATE(Reaction%Constants(1,1,1,Reaction%NumConstants))
  BACKSPACE Unit
  READ(Unit,*)  Reaction%NumConstants, &
               (Reaction%Constants(1,1,1,i),i=1, &
                Reaction%NumConstants)
  IF (Reaction%TypeR.eq.'SPECIAL') THEN
    CALL rem_blanks(string_aktuell)
    ALLOCATE(Current_infix)
    NULLIFY(Current_infix%next)
    Reaction%Infix=>Current_infix
    p=0
    Start=.TRUE.
    CALL delete_bracket(string_aktuell(1:len_trim(string_aktuell)))
  END IF
!--- Check if Factor is in Reaction
  READ(Unit,*,END=100) chem_factor
  IF (chem_factor(:1)=='F') THEN
    ALLOCATE(Reaction%Factor)
    BACKSPACE(Unit)
    READ(Unit,*) chem_factor,Reaction%Factor,Reaction%fac_exp,Reaction%fac_A
  ELSE 
   Backspace Unit 
   NULLIFY(Reaction%Factor)
 ENDIF
100 CONTINUE

CONTAINS

SUBROUTINE insert(operand,operator)
  CHARACTER(LEN=*) :: operand
  CHARACTER(LEN=*) :: operator

  IF (.NOT.Start) THEN
    ALLOCATE(current_infix%next)
    current_infix=>current_infix%next
    NULLIFY(current_infix%next)
  END IF
  Start=.FALSE.
  current_infix%operand=operand
  current_infix%op=operator
END SUBROUTINE insert

SUBROUTINE brack(string,brack_left,brack_right)
  CHARACTER(len=*) :: string
  INTEGER :: brack_left,brack_right

  INTEGER :: close
     
  brack_left=scan(string,'(')
  brack_right=brack_left
  close=1
  DO WHILE(close>0)
    brack_right=scan(string(brack_right+1:),'()')+brack_right
    IF (string(brack_right:brack_right).eq.'(') THEN
      close=close+1
    ELSE
      close=close-1
    END IF
  END DO 
END SUBROUTINE brack
      
SUBROUTINE rem_blanks(string)
  CHARACTER(len=*) :: string

  INTEGER :: i,j
  INTEGER :: len_str
  len_str=len(string)
  j=1
  DO i=1,len_str
    IF (string(i:i).ne.' ') THEN
      string(j:j)=string(i:i)
      j=j+1
    END IF
  END DO
  string(j:)=' '
END SUBROUTINE rem_blanks

FUNCTION rang(operator)
  INTEGER :: rang
  CHARACTER :: operator
       
  IF (scan(operator,'+-').eq.1) THEN
    rang=1
  ELSE IF (scan(operator,'*/').eq.1) THEN
    rang=2
  ELSE IF (scan(operator,'^').eq.1) THEN
    rang=3
  END IF
END function rang


SUBROUTINE analyze(string,ops,start,brack) 
  CHARACTER(len=*) :: string
  INTEGER :: ops,start,brack

  INTEGER :: pos_first_operator,current_rang
  INTEGER :: len_str
  CHARACTER :: first_operator

  first_operator=' '
  DO 
    pos_first_operator=scan(string,possible_operator)
    IF (pos_first_operator.gt.0) THEN
      first_operator=string(pos_first_operator:pos_first_operator)
      current_rang=rang(first_operator)
      IF (pos_first_operator.gt.1) THEN
        CALL insert(string(1:pos_first_operator-1),' ')
      END IF 
      string=string(pos_first_operator+1:)
    ELSE
      first_operator=' '
      current_rang=1
      IF (len_trim(string).gt.0) THEN
        CALL insert(string(1:len_trim(string)),leer)
      END IF
      string=' '
    END IF
    DO WHILE (ops.gt.0)
      IF (current_rang.le.stack_operator(p)%rang) THEN
        CALL insert(leer,stack_operator(p)%operator)
        ops=ops-1
        p=p-1
      ELSE
        EXIT
      END IF
    END DO
    IF (first_operator.ne.' ') THEN
      ops=ops+1
      p=p+1
      stack_operator(p)%operator=first_operator
      IF (pos_first_operator==1.and.first_operator=='-'.and.brack>0) THEN
        stack_operator(p)%operator='--'
      END IF
      stack_operator(p)%rang=current_rang
    END IF
    IF (len_trim(string).eq.0) THEN 
      EXIT
    END IF
  END DO
  start=start+1
END SUBROUTINE analyze

RECURSIVE SUBROUTINE delete_bracket(string)
  CHARACTER(len=*) :: string
  INTEGER :: brack_left,brack_right
  INTEGER :: len_str,current_rang,itemp
  INTEGER :: operator_before_brack
  INTEGER :: ops,start
  len_str=len_trim(string)

  IF (len_str>0) THEN
    ops=0
    start=0
    brack_left=1
    DO while (len_trim(string).gt.0)
      brack_left=scan(string,'(')
      IF (brack_left.gt.0) THEN
        operator_before_brack=scan(string(:brack_left),possible_operator,back)
        CALL analyze(string(:operator_before_brack),ops,start,brack_left)
        string=string(operator_before_brack+1:)
        brack_left=scan(string,'(')
        IF (brack_left.gt.1) THEN
          ops=ops+1
          p=p+1
          stack_operator(p)%operator=string(:brack_left-1)
          stack_operator(p)%rang=4
        END IF
        CALL brack(string,brack_left,brack_right)
        CALL delete_bracket(string(brack_left+1:brack_right-1))
        string=string(brack_right+1:)
      ELSE
        CALL analyze(string,ops,start,brack_left)
        string=' '
      END IF
    END DO
    CALL analyze(string,ops,start,brack_left)
  END IF
END SUBROUTINE delete_bracket
END SUBROUTINE ReactionInput

END MODULE Reaction_Mod


