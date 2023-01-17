MODULE Argsort_Mod

USE Kindmod

INTERFACE merge_argsort

  MODULE PROCEDURE merge_argsort_int, merge_argsort_real
  
END INTERFACE


CONTAINS

SUBROUTINE merge_argsort_real(rnum,d)
  IMPLICIT NONE

  REAL(Realkind), INTENT(IN) :: rnum(:)
  INTEGER, INTENT(OUT) :: d(SIZE(rnum))
           
  INTEGER :: il(SIZE(rnum))

  INTEGER :: stepsize
  INTEGER :: i,j,n,left,k,ksize
            
  n = SIZE(rnum)
            
  DO i=1,n
    d(i)=i
  END DO
            
  IF ( n==1 ) RETURN
            
  stepsize = 1
  DO WHILE (stepsize<n)
    DO left=1,n-stepsize,stepsize*2
      i = left
      j = left+stepsize
      ksize = min(stepsize*2,n-left+1)
      k=1
                
      DO WHILE ( i<left+stepsize .AND. j<left+ksize )
        IF ( rnum(d(i))>rnum(d(j)) ) THEN
          il(k)=d(i)
          i=i+1
          k=k+1
        ELSE
          il(k)=d(j)
          j=j+1
          k=k+1
        ENDIF
      END DO
      
      IF ( i<left+stepsize ) THEN
        ! fill up remaining from left
        il(k:ksize) = d(i:left+stepsize-1)
      ELSE
        ! fill up remaining from right
        il(k:ksize) = d(j:left+ksize-1)
      ENDIF
      d(left:left+ksize-1) = il(1:ksize)
    END DO
    stepsize=stepsize*2
  END DO

  RETURN

END SUBROUTINE merge_argsort_real

SUBROUTINE merge_argsort_int(inum,d)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: inum(:)
  INTEGER, INTENT(OUT) :: d(SIZE(inum))

  INTEGER :: il(SIZE(inum))

  INTEGER :: stepsize
  INTEGER :: i,j,n,left,k,ksize

  n = SIZE(inum)

  DO i=1,n
    d(i)=i
  END DO

  IF ( n==1 ) RETURN

  stepsize = 1
  DO WHILE (stepsize<n)
    DO left=1,n-stepsize,stepsize*2
      i = left
      j = left+stepsize
      ksize = min(stepsize*2,n-left+1)
      k=1

      DO WHILE ( i<left+stepsize .AND. j<left+ksize )
        IF ( inum(d(i))>inum(d(j)) ) THEN
          il(k)=d(i)
          i=i+1
          k=k+1
        ELSE
          il(k)=d(j)
          j=j+1
          k=k+1
        ENDIF
      END DO

      IF ( i<left+stepsize ) THEN
        ! fill up remaining from left
        il(k:ksize) = d(i:left+stepsize-1)
      ELSE
        ! fill up remaining from right
        il(k:ksize) = d(j:left+ksize-1)
      ENDIF
      d(left:left+ksize-1) = il(1:ksize)
    END DO
    stepsize=stepsize*2
  END DO

  RETURN

END SUBROUTINE merge_argsort_int

END MODULE Argsort_Mod
