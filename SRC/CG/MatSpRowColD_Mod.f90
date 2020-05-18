MODULE MatSpRowColD_Mod

  USE Kind_Mod

  IMPLICIT NONE

  TYPE SpRowColD
    INTEGER :: m,n
    INTEGER, POINTER :: RowPtr(:,:)=>NULL()
    INTEGER, POINTER :: ColInd(:)=>NULL()
    INTEGER, POINTER :: Permu(:)=>NULL()
    INTEGER, POINTER :: InvPer(:)=>NULL()
    INTEGER, POINTER :: Restr(:)=>NULL()
    INTEGER :: ep
    INTEGER :: last
    INTEGER :: len
  END TYPE SpRowColD

  INTEGER, PARAMETER, PRIVATE :: inilen=90

CONTAINS

SUBROUTINE SymbLU_SpRowColD(A)

  TYPE(SpRowColD) :: A

  INTEGER :: r(A%n),c(A%n),RowPiv(A%n)
  INTEGER :: i,j,l,jj,ip,ip1(1),iPiv
  REAL(RealKind) :: md 
  LOGICAL :: ins

  c=0
  DO i=1,A%n
    A%InvPer(i)=i
    A%Permu(i)=i
!   Compute initial Markowitz count
    r(i)=A%RowPtr(2,i)-A%RowPtr(1,i)+1
    DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
      c(A%ColInd(jj))=c(A%ColInd(jj))+1
    END DO
  END DO

  DO i=1,A%n
    ip=0
    md=1.e-40_RealKind
    DO j=i,A%n
      IF (((r(j)-1)*(c(j)-1)<=md).AND.(A%Restr(j)>=0)) THEN
        md=(r(j)-1)*(c(j)-1)
        ip=j
      END IF
    END DO
    IF (ip==0) THEN
      ip1(:)=MINLOC((r(i:A%n)-1)*(c(i:A%n)-1))+(i-1)
      ip=ip1(1)
    END IF
    CALL Swap(A%InvPer(i),A%InvPer(ip))
    CALL Swap(r(i),r(ip))
    CALL Swap(c(i),c(ip))
    CALL Swap(A%Restr(i),A%Restr(ip))
    CALL Swap(A%RowPtr(1,i),A%RowPtr(1,ip))
    CALL Swap(A%RowPtr(2,i),A%RowPtr(2,ip))
    A%Permu(A%InvPer(i))=i
    A%Permu(A%InvPer(ip))=ip
    IF (A%last==i) THEN
      A%last=ip
    ELSE IF (A%last==ip) THEN
      A%last=i
    END IF

!   Update
    iPiv=0
    DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
      IF (A%Permu(A%ColInd(jj))>i) THEN
        iPiv=iPiv+1 
        RowPiv(iPiv)=A%ColInd(jj)
        c(A%Permu(A%ColInd(jj)))=c(A%Permu(A%ColInd(jj)))-1 
      END IF
    END DO
    IF (iPiv>0) THEN
      DO j=i+1,A%n
        DO jj=A%RowPtr(1,j),A%RowPtr(2,j)
          IF (A%Permu(A%ColInd(jj))==i) THEN
            r(j)=r(j)-1 
            DO l=1,iPiv
              CALL Insert_SpRowColD(A,j,RowPiv(l),ins)
              IF (ins) THEN
                c(A%Permu(RowPiv(l)))=c(A%Permu(RowPiv(l)))+1
                r(j)=r(j)+1
              END IF
            END DO
            EXIT 
          END IF
        END DO
      END DO
    END IF
  END DO
  DO i=1,A%n
    DO jj=A%RowPtr(1,i),A%RowPtr(2,i)
      A%ColInd(jj)=A%Permu(A%ColInd(jj))
    END DO
    CALL Sort(A%ColInd(A%RowPtr(1,i):A%RowPtr(2,i)))
  END DO
CONTAINS

SUBROUTINE Swap(i,j)

  INTEGER :: i,j
  INTEGER :: iTemp
   
  iTemp=i
  i=j
  j=iTemp
END SUBROUTINE Swap

SUBROUTINE sort(vec)
  INTEGER :: vec(:)

  INTEGER :: i,itemp,j,n
  n=SIZE(Vec)
  DO i=1,n
    DO j=1,n-i
      IF (vec(j).gt.vec(j+1)) THEN
        itemp=vec(j)
        vec(j)=vec(j+1)
        vec(j+1)=itemp
      END IF
    END DO
  END DO
END SUBROUTINE sort
  
END SUBROUTINE SymbLU_SpRowColD
        
SUBROUTINE Insert_SpRowColD(A,iA,jA,ins)

  TYPE(SpRowColD) :: A
  INTEGER :: iA,jA
  LOGICAL :: ins

  INTEGER :: itemp,j,l
  INTEGER, ALLOCATABLE :: iWork(:)

! Test ob Element (ia,ja) bereits enthalten

  ins=.TRUE.
  DO j=A%RowPtr(1,iA),A%RowPtr(2,iA)
    IF (jA==A%ColInd(j)) THEN
      ins=.FALSE.
    END IF
  END DO

! Test auf freien Speicherplatz in der ia-ten
! Zeile von a

  IF (ins) THEN
    IF (A%ColInd(A%RowPtr(2,iA)+1)/=0) THEN
!     ja-te Zeile von a wird nach hinten
!     geschrieben
      itemp=A%ep
      DO l=A%RowPtr(1,iA),A%RowPtr(2,iA)
        A%ColInd(A%ep)=A%ColInd(l)
        A%ColInd(l)=0
        A%ep=A%ep+1
      END DO
      A%RowPtr(2,iA)=A%ep-1
      A%RowPtr(1,iA)=itemp
!     A%ep=A%ep+1
      A%last=iA
    ENDIF
    A%RowPtr(2,iA)=A%RowPtr(2,iA)+1
    A%ColInd(A%RowPtr(2,iA))=jA
    IF (iA==A%last) A%ep=A%ep+1
  END IF

  IF (A%ep.ge.A%len-A%n) THEN
    CALL gcmat_SpRowColD(A)
  END IF

  IF (A%ep.ge.A%len-A%n) THEN
!   Speicherplatz von A nicht ausreichend
    ALLOCATE(iWork(A%ep))
    iWork(1:A%ep)=A%ColInd(1:A%ep)
    DEALLOCATE(A%ColInd)
    A%len=2*A%len
    ALLOCATE(A%ColInd(A%len))
    A%ColInd(1:A%ep)=iWork(1:A%ep)
    DEALLOCATE(iWork)  
  END IF

END SUBROUTINE Insert_SpRowColD


SUBROUTINE gcmat_SpRowColD(A)
 
!   Externe Variable

  TYPE (SpRowColD) :: A
 
 
!   gcmat komprimiert eine zeilenorientierte, dynamische
!   Speicherstruktur einer schwachbesetzten Matrix a
!   der Ordnung n. Die Spaltenindizes der Nichtnull-
!   elemente der i-ten Zeile von a sind in A%ColInd(A%RowPtr(i,1)),
!   A%ColInd(A%RowPtr(1,i)+1)...,A%ColInd(A%RowPtr(2,i)) ent-
!   halten.
 
!    Beschreibung der Parameter
 
   
!    A%n      (i/o) integer
!                 Dimension of matrix a
 
!    A%RowPtr (i/o) integer(2,n)
!                 Pointerfeld zur Beschreibung von a. A%RowPtr muss
!                 durch das rufende Programm belegt werden.
 
!    A%ColInd (i/o) integer(len)
!                 Feld zur dynamischen Verwaltung der Spaltenindizes
!                 der Nichtnullelemente von a. 
 
!    A%ep     (i/o) integer
!                 Pointer der auf das erste freie Feld in a verweist,
!                 d.h A%ColInd(ep),...,A%ColInd(len) sind frei verfuegbar.
 
!    A%len    (i)   integer
!                 Gesamtlaenge des Feldes A%ColInd. 
 
!  Interne Variable
 
   INTEGER i,iz,j,k,l,pointr,rowlen,ep,len,n

 
   n=A%n
   ep=A%ep
   len=A%len
   pointr=1
   i=1
 
   DO 
      IF (i.ge.ep) EXIT
         IF (A%ColInd(i).ne.0) THEN

!           Ermittlung der aktuellen Zeile sowie deren Laenge
  
            DO l=1,n
               IF (A%RowPtr(1,l).le.i.and.i.le.A%RowPtr(2,l)) THEN
                  iz=l
               END IF
            END DO
            rowlen=A%RowPtr(2,iz)-A%RowPtr(1,iz)
 
!           Setzen der neuen Anfangs- und Endadresse der
!           aktuellen Zeile
 
            A%RowPtr(1,iz)=pointr
            A%RowPtr(2,iz)=pointr+rowlen
            DO j=pointr,pointr+rowlen
               A%ColInd(j)=A%ColInd(i)
               i=i+1
            END DO
            i=i-1
            pointr=A%RowPtr(2,iz)+1
         ENDIF
         i=i+1
    END DO

!   Belegung des freien Teils von A%ColInd mit 0
 
    ep=pointr
    DO i=1,n
       IF (A%RowPtr(1,i).gt.ep) THEN
          A%RowPtr(1,i)=ep
          A%RowPtr(2,i)=A%RowPtr(1,i)-1
          ep=ep+inilen
       END IF
    END DO
            

    DO i=pointr,len
       A%ColInd(i)=0
    END DO
    A%ep=ep
 
END SUBROUTINE gcmat_SpRowColD

SUBROUTINE SpNullify_SpRowColD(A)

!
! Nullify Pointers of A 
!
  TYPE(SpRowColD) :: A

  NULLIFY(A%RowPtr)
  NULLIFY(A%ColInd)
  NULLIFY(A%Permu)
  NULLIFY(A%InvPer)

END SUBROUTINE SpNullify_SpRowColD

SUBROUTINE SpDeallocate_SpRowColD(A)

!
! Deallocate A
!
  TYPE(SpRowColD) :: A

  A%m=0
  A%n=0
  IF (ASSOCIATED(A%RowPtr)) THEN
    DEALLOCATE(A%RowPtr)
  END IF
  IF (ASSOCIATED(A%ColInd)) THEN
    DEALLOCATE(A%ColInd)
  END IF
  IF (ASSOCIATED(A%Permu)) THEN
    DEALLOCATE(A%Permu)
  END IF
  IF (ASSOCIATED(A%InvPer)) THEN
    DEALLOCATE(A%InvPer)
  END IF

END SUBROUTINE SpDeallocate_SpRowColD

SUBROUTINE SpAllocate_SpRowColD(A)

!
! Allocate A
!
  TYPE(SpRowColD) :: A

  INTEGER :: i

  IF (.NOT.ASSOCIATED(A%Permu)) THEN
    ALLOCATE(A%Permu(A%n))
  END IF
  IF (.NOT.ASSOCIATED(A%InvPer)) THEN
    ALLOCATE(A%InvPer(A%n))
  END IF
  IF (.NOT.ASSOCIATED(A%Restr)) THEN
    ALLOCATE(A%Restr(A%n))
  END IF
  IF (.NOT.ASSOCIATED(A%RowPtr)) THEN
    ALLOCATE(A%RowPtr(2,A%n))
  END IF
  IF (.NOT.ASSOCIATED(A%ColInd)) THEN
    IF (A%len<A%n*inilen) THEN
      A%len=A%n*inilen+2*A%n
    END IF
    ALLOCATE(A%ColInd(A%len))
    A%ColInd=0
    DO i=1,A%n
      A%RowPtr(1,i)=(i-1)*inilen+1
      A%RowPtr(2,i)= A%RowPtr(1,i)-1
    END DO
    A%ep=A%RowPtr(2,A%n)+1
    A%last=A%n
  END IF

END SUBROUTINE SpAllocate_SpRowColD

END MODULE MatSpRowColD_Mod
  

