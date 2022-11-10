MODULE Sort_Mod

USE Kind_Mod

INTERFACE sortunique
  MODULE PROCEDURE rsortunique, isortunique
END INTERFACE sortunique

INTERFACE mergesort
  MODULE PROCEDURE rmergesort, imergesort
END INTERFACE mergesort

CONTAINS

SUBROUTINE rsortunique(a, unique)

  REAL(Realkind), INTENT(inout) :: a(:)
  REAL(Realkind), INTENT(inout), POINTER :: unique(:)
  REAL(Realkind) :: unique_tmp(SIZE(a, 1))
  INTEGER :: n, i, iunique

  IF (ASSOCIATED(unique)) THEN
    DEALLOCATE(unique)
  END IF

  n = SIZE(a, 1)

  IF (n .EQ. 1) THEN
    ALLOCATE(unique(n))
    unique(:) = a
    RETURN
  END IF


  CALL MergeSort(a)

  unique_tmp(1) = a(1)
  iunique = 1
  DO i = 2, n
    IF (a(i) .ne. unique_tmp(iunique)) THEN
      iunique = iunique + 1
      unique_tmp(iunique) = a(i) 
    END IF
  END DO

  ALLOCATE(unique(iunique))
  unique(1:iunique) = unique_tmp(1:iunique)

END SUBROUTINE  rsortunique


SUBROUTINE isortunique(a, unique)

  INTEGER, INTENT(inout) :: a(:)
  INTEGER, POINTER, INTENT(inout) :: unique(:)
  INTEGER :: unique_tmp(SIZE(a, 1))
  INTEGER :: n, i, iunique

  IF (ASSOCIATED(unique)) THEN
    DEALLOCATE(unique)
  END IF

  n = SIZE(a, 1)

  IF (n .EQ. 1) THEN
    ALLOCATE(unique(n))
    unique(:) = a
    RETURN
  END IF


  CALL MergeSort(a)

  unique_tmp(1) = a(1)
  iunique = 1
  DO i = 2, n
    IF (a(i) .ne. unique_tmp(iunique)) THEN
      iunique = iunique + 1
      unique_tmp(iunique) = a(i)
    END IF
  END DO

  ALLOCATE(unique(iunique))
  unique(1:iunique) = unique_tmp(1:iunique)

END SUBROUTINE isortunique


SUBROUTINE RMergeSort(a)
  REAL(Realkind), INTENT(inout) :: a(:)
  INTEGER :: n
  REAL(Realkind) :: t((SIZE(a, 1) + 1) / 2)

  n = SIZE(a, 1)

  CALL RMergeSort_rec(a, n, t)

END SUBROUTINE RMergeSort


SUBROUTINE IMergeSort(a)
  INTEGER, INTENT(inout) :: a(:)
  INTEGER :: n
  INTEGER :: t((SIZE(a, 1) + 1) / 2)

  n = SIZE(a, 1)

  CALL IMergeSort_rec(a, n, t)

END SUBROUTINE IMergeSort

RECURSIVE SUBROUTINE RMergeSort_rec(a, n, t)
  REAL(Realkind), INTENT(inout) :: a(:)
  INTEGER, INTENT(in) :: n
  REAL(Realkind), INTENT(out) :: t(:)
  REAL :: v
  INTEGER :: na, nb

  IF (n < 2) RETURN
  IF (n .EQ. 2) THEN
    IF (a(1) .GT. a(2)) THEN
      v = a(1)
      a(1) = a(2)
      a(2) = v
    END IF
    RETURN
  END IF

  na = (n + 1) / 2
  nb = n - na
  
  CALL RMergeSort_rec(a(1:na), na, t)
  CALL RMergeSort_rec(a(na + 1:n), nb, t)

  IF (a(na) .GT. a(na + 1)) THEN
    T(1:na) = a(1:na)
    CALL Rmerge(t, na, a(na + 1:n), nb, a, n)
  END IF

END SUBROUTINE RMergeSort_rec


RECURSIVE SUBROUTINE IMergeSort_rec(a, n, t)
  INTEGER, INTENT(inout) :: a(:)
  INTEGER, INTENT(in) :: n
  INTEGER, INTENT(out) :: t(:)
  INTEGER :: v
  INTEGER :: na, nb

  IF (n < 2) RETURN
  IF (n .EQ. 2) THEN
    IF (a(1) .GT. a(2)) THEN
      v = a(1)
      a(1) = a(2)
      a(2) = v
    END IF
    RETURN
  END IF

  na = (n + 1) / 2
  nb = n - na

  CALL IMergeSort_rec(a(1:na), na, t)
  CALL IMergeSort_rec(a(na + 1:n), nb, t)

  IF (a(na) .GT. a(na + 1)) THEN
    T(1:na) = a(1:na)
    CALL Imerge(t, na, a(na + 1:n), nb, a, n)
  END IF

END SUBROUTINE IMergeSort_rec


SUBROUTINE Rmerge(a, na, b, nb, c, nc)

  INTEGER, INTENT(in) :: na, nb, nc

  REAL(Realkind), INTENT(inout) :: a(:)
  REAL(Realkind), INTENT(in) :: b(:)
  REAL(Realkind), INTENT(inout) :: c(:)

  INTEGER :: i, j, k

  i = 1
  j = 1
  k = 1

  DO WHILE (i .LE. na .and. j .LE. nb) 
    IF (a(i) .LE. b(j)) THEN
      c(k) = a(i)
      i = i + 1
    ELSE
      c(k) = b(j)
      j = j + 1
    END IF
    k = k + 1
  END DO
  DO WHILE (i .LE. na) 
    c(k) = a(i)
    i = i + 1
    k = k + 1
  END DO

END SUBROUTINE Rmerge


SUBROUTINE Imerge(a, na, b, nb, c, nc)

  INTEGER, INTENT(in) :: na, nb, nc

  INTEGER, INTENT(inout) :: a(:)
  INTEGER, INTENT(in) :: b(:)
  INTEGER, INTENT(inout) :: c(:)

  INTEGER :: i, j, k

  i = 1
  j = 1
  k = 1

  DO WHILE (i .LE. na .and. j .LE. nb)
    IF (a(i) .LE. b(j)) THEN
      c(k) = a(i)
      i = i + 1
    ELSE
      c(k) = b(j)
      j = j + 1
    END IF
    k = k + 1
  END DO
  DO WHILE (i .LE. na)
    c(k) = a(i)
    i = i + 1
    k = k + 1
  END DO

END SUBROUTINE Imerge

END Module Sort_Mod
