MODULE List_Mod

USE Kind_Mod
USE Sort_Mod, ONLY : sortunique

IMPLICIT NONE

TYPE List

  TYPE(Element), POINTER :: first => NULL()
  INTEGER :: len = 0
  INTEGER :: len_flat = 0
  CONTAINS

  GENERIC :: append => &
             append_integer, &
             append_integer_array1d, &
             append_integer_array2d, &
             append_integer_array3d, & 
             append_real, &
             append_real_array1d, &
             append_real_array2d, &
             append_real_array3d, &
             append_list

  GENERIC :: insert => &
             insert_integer, &
             insert_integer_array1d, &
             insert_integer_array2d, &
             insert_integer_array3d, &
             insert_real, &
             insert_real_array1d, &
             insert_real_array2d, &
             insert_real_array3d, &
             insert_list

  GENERIC :: index => &
             iindex, &
             rindex 

  GENERIC :: toarray => &
             toiarray, &
             torarray

  GENERIC :: unique => &
             iunique, &
             runique

  PROCEDURE ::                 &
               append_integer, &
               append_integer_array1d, &
               append_integer_array2d, &
               append_integer_array3d, &
               append_real, &
               append_real_array1d, &
               append_real_array2d, &
               append_real_array3d, &
               append_list, &
               insert_integer, &
               insert_integer_array1d, &
               insert_integer_array2d, &
               insert_integer_array3d, &
               insert_real, &
               insert_real_array1d, &
               insert_real_array2d, &
               insert_real_array3d, &
               insert_list

  PROCEDURE :: new_sublist
  PROCEDURE :: get, count_deep

  PROCEDURE :: remove, delete, destroy, free
  PROCEDURE :: iindex, rindex

  PROCEDURE :: toiarray, icollect, torarray, rcollect
  PROCEDURE :: iunique, runique

END TYPE List



TYPE Element
  REAL(Realkind), ALLOCATABLE :: rvalue
  REAL(Realkind), POINTER :: rarray1d(:) => NULL()
  REAL(Realkind), POINTER :: rarray2d(:,:) => NULL()
  REAL(Realkind), POINTER :: rarray3d(:,:,:) => NULL()
  
  INTEGER, ALLOCATABLE :: ivalue
  INTEGER, POINTER :: iarray1d(:) => NULL()
  INTEGER, POINTER :: iarray2d(:,:) => NULL()
  INTEGER, POINTER :: iarray3d(:,:,:) => NULL()

  TYPE(List), POINTER :: sublist => NULL()
  TYPE(Element), POINTER :: next => NULL()

  CONTAINS 

  PROCEDURE :: nullify_pointers, remove_data

END TYPE Element

TYPE ListArray
   TYPE(List) :: lst
END TYPE ListArray

TYPE RArray
  REAL(Realkind), POINTER :: data(:) => NULL()
END TYPE RArray

TYPE IArray
  INTEGER, POINTER :: data(:) => NULL()
END TYPE IArray

TYPE IArray2d
  INTEGER, POINTER :: data(:,:) => NULL()
END TYPE IArray2d


CONTAINS


SUBROUTINE append_integer(self, i)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: i
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next)) 
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  ALLOCATE(current%ivalue)
  current%ivalue = i
  self%len = self%len + 1  
  
END SUBROUTINE append_integer


SUBROUTINE insert_integer(self, pos, i)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  INTEGER, INTENT(in) :: i
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first

  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  self%len = self%len + 1

  ALLOCATE(inserted)
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF
  
  ALLOCATE(inserted%ivalue)
  inserted%ivalue = i

END SUBROUTINE insert_integer


SUBROUTINE append_integer_array1d(self, i, copy)

  IMPLICIT NONE
  
  CLASS(LIST), TARGET :: self
  INTEGER, TARGET, INTENT(in) :: i(:)
  LOGICAL, OPTIONAL, INTENT(in) :: copy
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  self%len = self%len + 1

  IF (PRESENT(copy)) THEN
    IF (copy .EQV. .True.) THEN
      ALLOCATE(current%iarray1d(SIZE(i, 1)))
      current%iarray1d(:) = i
      RETURN
    END IF
  END IF

  current%iarray1d => i


END SUBROUTINE append_integer_array1d


SUBROUTINE insert_integer_array1d(self, pos, i)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  INTEGER, TARGET, INTENT(in) :: i(:)
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first
    
  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  ALLOCATE(inserted)
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF
  
  inserted%iarray1d => i
  self%len = self%len + 1   

END SUBROUTINE insert_integer_array1d


SUBROUTINE append_integer_array2d(self, i, copy)

  IMPLICIT NONE
  
  CLASS(LIST), TARGET :: self
  INTEGER, TARGET, INTENT(in) :: i(:, :)
  LOGICAL, OPTIONAL, INTENT(in) :: copy
  TYPE(Element), POINTER :: current


  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  IF (PRESENT(copy)) THEN
    IF (copy .EQV. .True.) THEN
      ALLOCATE(current%iarray2d(SIZE(i, 1), SIZE(i, 2)))
      current%iarray2d(:, :) = i
      RETURN
    END IF
  END IF

  current%iarray2d => i
  self%len = self%len + 1

END SUBROUTINE append_integer_array2d


SUBROUTINE insert_integer_array2d(self, pos, i)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  INTEGER, TARGET, INTENT(in) :: i(:, :)
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first

  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  ALLOCATE(inserted)
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF

  inserted%iarray2d => i
  self%len = self%len + 1

END SUBROUTINE insert_integer_array2d


SUBROUTINE append_integer_array3d(self, i)

  IMPLICIT NONE
  
  CLASS(LIST), TARGET :: self
  INTEGER, TARGET, INTENT(in) :: i(:, :, :)
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  current%iarray3d => i
  self%len = self%len + 1

END SUBROUTINE append_integer_array3d


SUBROUTINE insert_integer_array3d(self, pos, i)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  INTEGER, TARGET, INTENT(in) :: i(:, :, :)
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first

  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  ALLOCATE(inserted)
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF

  inserted%iarray3d => i
  self%len = self%len + 1

END SUBROUTINE insert_integer_array3d


SUBROUTINE append_real(self, r)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  REAL(Realkind), INTENT(in) :: r
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  ALLOCATE(current%rvalue)
  current%rvalue = r
  self%len = self%len + 1

END SUBROUTINE append_real


SUBROUTINE insert_real(self, pos, r)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  REAL(Realkind), INTENT(in) :: r
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first

  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  ALLOCATE(inserted)
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF

  ALLOCATE(inserted%rvalue)
  inserted%rvalue = r
  self%len = self%len + 1

END SUBROUTINE insert_real


SUBROUTINE append_real_array1d(self, r)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  REAL(Realkind), TARGET, INTENT(in) :: r(:)
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  current%rarray1d => r
  self%len = self%len + 1

END SUBROUTINE append_real_array1d


SUBROUTINE insert_real_array1d(self, pos, r)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  REAL(Realkind), TARGET, INTENT(in) :: r(:)
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first

  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  ALLOCATE(inserted)
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF

  inserted%rarray1d => r
  self%len = self%len + 1

END SUBROUTINE insert_real_array1d




SUBROUTINE append_real_array2d(self, r)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  REAL(Realkind), TARGET, INTENT(in) :: r(:,:)
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  current%rarray2d => r
  self%len = self%len + 1

END SUBROUTINE append_real_array2d


SUBROUTINE insert_real_array2d(self, pos, r)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  REAL(Realkind), TARGET, INTENT(in) :: r(:, :)
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first

  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  ALLOCATE(inserted)    
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF
  
  inserted%rarray2d => r
  self%len = self%len + 1

END SUBROUTINE insert_real_array2d


SUBROUTINE append_real_array3d(self, r)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  REAL(Realkind), TARGET, INTENT(in) :: r(:,:,:)
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  current%rarray3d => r
  self%len = self%len + 1

END SUBROUTINE append_real_array3d


SUBROUTINE insert_real_array3d(self, pos, r)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  REAL(Realkind), TARGET, INTENT(in) :: r(:,:,:)
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first

  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  ALLOCATE(inserted)    
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF
  
  inserted%rarray3d => r
  self%len = self%len + 1

END SUBROUTINE insert_real_array3d


SUBROUTINE append_list(self, lst)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  TYPE(List), TARGET, INTENT(in) :: lst
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  current%sublist => lst
  self%len = self%len + 1

END SUBROUTINE append_list


SUBROUTINE insert_list(self, pos, lst)

  IMPLICIT NONE

  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: pos
  TYPE(List), TARGET, INTENT(in) :: lst
  TYPE(Element), POINTER :: current, inserted
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) ALLOCATE(self%first)
  current => self%first

  DO j = 1, pos - 2
    IF (ASSOCIATED(current%next)) THEN
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  ALLOCATE(inserted)
  IF (ASSOCIATED(current%next)) THEN
    inserted%next => current%next
    current%next => inserted
  END IF

  inserted%sublist => lst
  self%len = self%len + 1

END SUBROUTINE insert_list

FUNCTION new_sublist(self)

  IMPLICIT NONE
  CLASS(LIST), TARGET :: self
  TYPE(List), POINTER :: new_sublist
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) THEN
    ALLOCATE(self%first)
    current => self%first
  ELSE 
    current => self%first
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
    ALLOCATE(current%next)
    current => current%next
  END IF

  ALLOCATE(current%sublist)

  new_sublist => current%sublist

  self%len = self%len + 1

END FUNCTION new_sublist


FUNCTION iindex(self, value) RESULT(i)
  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(in) :: value
  INTEGER :: i
  TYPE(Element), POINTER :: current
 
  i = 1

  IF (.NOT. ASSOCIATED(self%first)) THEN
    WRITE(*,*) 'Error: List is empty'
    STOP
  END IF
  current => self%first

  IF (ALLOCATED(current%ivalue)) THEN
    IF (value .EQ. current%ivalue) RETURN
  END IF
  
  DO WHILE(ASSOCIATED(current%next))
    current => current%next
    i = i + 1
    IF (ALLOCATED(current%ivalue)) THEN
      IF (value .EQ. current%ivalue) RETURN
    END IF
  END DO

  WRITE(*,*) "Index Error: List does not contain requested integer value"
  STOP

END FUNCTION iindex

FUNCTION rindex(self, value) RESULT(i)
  CLASS(LIST), TARGET :: self
  REAL(Realkind), INTENT(in) :: value
  INTEGER :: i
  TYPE(Element), POINTER :: current
 
  i = 1

  IF (.NOT. ASSOCIATED(self%first)) THEN
    WRITE(*,*) 'Error: List is empty'
    STOP
  END IF
  current => self%first

  IF (ALLOCATED(current%rvalue)) THEN
    IF (value .EQ. current%rvalue) RETURN
  END IF
  
  DO WHILE(ASSOCIATED(current%next))
    current => current%next
    i = i + 1
    IF (ALLOCATED(current%rvalue)) THEN
      IF (value .EQ. current%rvalue) RETURN
    END IF
  END DO

  WRITE(*,*) "Index Error: List does not contain requested real value"
  STOP

END FUNCTION rindex



SUBROUTINE remove(self, i)

  CLASS(LIST), TARGET :: self
  TYPE(Element), POINTER :: previous, current

  INTEGER, INTENT(in) :: i
  INTEGER :: j 

  IF (.NOT. ASSOCIATED(self%first)) THEN
    WRITE(*,*) 'Error: List is empty'
    STOP
  END IF
  current => self%first

  DO j = 1, i - 1
    IF (ASSOCIATED(current%next)) THEN
      previous => current
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO
  
  IF (ASSOCIATED(current%next)) previous%next => current%next

  CALL current%nullify_pointers()

  current => NULL()
  self%len = self%len - 1
END SUBROUTINE remove


SUBROUTINE delete(self, i)

  CLASS(LIST), TARGET :: self
  TYPE(Element), POINTER :: previous, current

  INTEGER, INTENT(in) :: i
  INTEGER :: j

  IF (.NOT. ASSOCIATED(self%first)) THEN
    WRITE(*,*) 'Error: List is empty'
    STOP
  END IF
  current => self%first

  DO j = 1, i - 1
    IF (ASSOCIATED(current%next)) THEN
      previous => current
      current => current%next
    ELSE
      WRITE(*,*) 'Error: Number of elements is smaller than position requested'
      STOP
    END IF
  END DO

  IF (ASSOCIATED(current%next)) previous%next => current%next

  CALL current%remove_data()

  current => NULL()
  self%len = self%len - 1

END SUBROUTINE delete


RECURSIVE SUBROUTINE destroy(self)
  CLASS(LIST), TARGET :: self
  TYPE(Element), POINTER :: current, next

  IF (.NOT. ASSOCIATED(self%first)) RETURN

  current => self%first
  CALL current%remove_data()

  DO WHILE (ASSOCIATED(current%next))

    next => current%next
    current%next => NULL()
    current => next
    CALL current%remove_data()
  END DO
  IF (ASSOCIATED(next)) next => NULL()
  current => NULL()
  self%first => NULL()
  self%len = 0
END SUBROUTINE destroy


RECURSIVE SUBROUTINE free(self)
  CLASS(LIST), TARGET :: self
  TYPE(Element), POINTER :: current, next

  IF (.NOT. ASSOCIATED(self%first)) RETURN

  current => self%first
  CALL current%nullify_pointers()

  DO WHILE (ASSOCIATED(current%next))
    
    next => current%next
    current%next => NULL()
    current => next
    CALL current%nullify_pointers()
  END DO
  IF (ASSOCIATED(next)) next => NULL()
  current => NULL()
  self%first => NULL()
  self%len = 0

END SUBROUTINE free


SUBROUTINE nullify_pointers(self)
  IMPLICIT NONE

  CLASS(Element) :: self

  IF (ALLOCATED(self%ivalue)) THEN
    DEALLOCATE(self%ivalue)
  ELSE  IF (ALLOCATED(self%rvalue)) THEN
    DEALLOCATE(self%rvalue)
  ELSE IF (ASSOCIATED(self%iarray1d)) THEN
    self%iarray1d => NULL()
  ELSE IF (ASSOCIATED(self%iarray2d)) THEN
    self%iarray2d => NULL()
  ELSE IF (ASSOCIATED(self%iarray3d)) THEN
    self%iarray3d => NULL()
  ELSE IF (ASSOCIATED(self%rarray1d)) THEN
    self%rarray1d => NULL()
  ELSE IF (ASSOCIATED(self%rarray2d)) THEN
    self%rarray2d => NULL()
  ELSE IF (ASSOCIATED(self%rarray3d)) THEN
    self%rarray3d => NULL()
  ELSE IF (ASSOCIATED(self%sublist)) THEN
    CALL self%sublist%free() 
    DEALLOCATE(self%sublist)
    self%sublist => NULL()
  END IF

END SUBROUTINE nullify_pointers

SUBROUTINE remove_data(self)
  IMPLICIT NONE

  CLASS(Element) :: self

  IF (ALLOCATED(self%ivalue)) THEN
    DEALLOCATE(self%ivalue)
  ELSE IF (ALLOCATED(self%rvalue)) THEN
    DEALLOCATE(self%rvalue)
  ELSE IF (ASSOCIATED(self%iarray1d)) THEN
    DEALLOCATE(self%iarray1d)
    self%iarray1d => NULL()
  ELSE IF (ASSOCIATED(self%iarray2d)) THEN
    DEALLOCATE(self%iarray2d)
    self%iarray2d => NULL()
  ELSE IF (ASSOCIATED(self%iarray3d)) THEN
    DEALLOCATE(self%iarray3d)
    self%iarray3d => NULL()
  ELSE IF (ASSOCIATED(self%rarray1d)) THEN
    DEALLOCATE(self%rarray1d)
    self%rarray1d => NULL()
  ELSE IF (ASSOCIATED(self%rarray2d)) THEN
    DEALLOCATE(self%rarray2d)
    self%rarray2d => NULL()
  ELSE IF (ASSOCIATED(self%rarray3d)) THEN
    DEALLOCATE(self%rarray3d)
    self%rarray3d => NULL()
  ELSE IF (ASSOCIATED(self%sublist)) THEN
    CALL self%sublist%destroy()
    DEALLOCATE(self%sublist)
    self%sublist => NULL()
  END IF

END SUBROUTINE remove_data


FUNCTION get(self, i) RESULT(current)
  CLASS(LIST), TARGET :: self
  TYPE(Element), POINTER :: current

  INTEGER, INTENT(in) :: i
  INTEGER :: j  

  IF (.NOT. ASSOCIATED(self%first)) THEN
    WRITE(*,*) 'Error: List is empty'
    STOP
  END IF
  current => self%first

  IF (i .EQ. -1) THEN
    DO WHILE(ASSOCIATED(current%next))
      current => current%next
    END DO
  ELSE
    DO j = 1, i - 1
      IF (ASSOCIATED(current%next)) THEN
        current => current%next
      ELSE
        WRITE(*,*) 'Error: Number of elements is smaller than position requested'
        STOP
      END IF
    END DO
  END IF

END FUNCTION get


RECURSIVE SUBROUTINE count_deep(self)
  CLASS(LIST), TARGET :: self
  TYPE(Element), POINTER :: current


  self%len_flat = 0

  IF (.NOT. ASSOCIATED(self%first)) RETURN

  current => self%first
  IF (ASSOCIATED(self%first%sublist)) THEN
    CALL self%first%sublist%count_deep()
    self%len_flat = self%len_flat + self%first%sublist%len_flat
  ELSE
    self%len_flat = self%len_flat + 1
  END IF
  DO WHILE(ASSOCIATED(current%next))
    current => current%next
    IF (ASSOCIATED(current%sublist)) THEN
      CALL current%sublist%count_deep()
      self%len_flat = self%len_flat + current%sublist%len_flat
    ELSE
      self%len_flat = self%len_flat + 1
    END IF
  END DO
END SUBROUTINE count_deep


SUBROUTINE iunique(self, unique)
  CLASS(LIST), TARGET :: self
  INTEGER, POINTER, INTENT(inout) :: unique(:)
  INTEGER, TARGET :: arr_tmp(self%len_flat)
  INTEGER, POINTER :: flattened(:)

  flattened => arr_tmp
  CALL self%toiarray(flattened)
  CALL sortunique(flattened, unique)

  flattened => NULL()

END SUBROUTINE iunique


SUBROUTINE toiarray(self, arr)
  CLASS(LIST), TARGET :: self
  INTEGER, POINTER, INTENT(inout) :: arr(:)
  INTEGER :: n

  IF (self%len_flat .EQ. 0) CALL self%count_deep()

  IF (.NOT. ASSOCIATED(arr)) ALLOCATE(arr(self%len_flat))
  n = 0
  CALL self%icollect(arr, n)
  
END SUBROUTINE toiarray


RECURSIVE SUBROUTINE icollect(self, array, n)
  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(inout) :: n
  INTEGER, POINTER, INTENT(inout)  :: array(:)
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) RETURN
  current => self%first
  IF (ASSOCIATED(current%sublist)) THEN
    CALL current%sublist%icollect(array, n)
  ELSE IF (ALLOCATED(current%ivalue)) THEN
    n = n + 1
    array(n) = current%ivalue
  END IF

  DO WHILE(ASSOCIATED(current%next))
    current => current%next
    IF (ASSOCIATED(current%sublist)) THEN
      CALL current%sublist%icollect(array, n)
    ELSE IF (ALLOCATED(current%ivalue)) THEN
      n = n + 1
      array(n) = current%ivalue
    END IF
  END DO

END SUBROUTINE icollect

SUBROUTINE  runique(self, unique)
  CLASS(LIST), TARGET :: self
  REAL(Realkind), POINTER, INTENT(inout) :: unique(:)
  REAL(Realkind), TARGET :: arr_tmp(self%len_flat)
  REAL(Realkind), POINTER :: flattened(:)

  flattened => arr_tmp
  CALL self%torarray(flattened)
  CALL sortunique(flattened, unique)

  flattened => NULL()

END SUBROUTINE runique


SUBROUTINE torarray(self, arr)
  CLASS(LIST), TARGET :: self
  REAL(Realkind), POINTER, INTENT(inout) :: arr(:)
  INTEGER :: n

  IF (self%len_flat .EQ. 0) CALL self%count_deep()

  IF (.NOT. ASSOCIATED(arr)) ALLOCATE(arr(self%len_flat))
  n = 0
  CALL self%rcollect(arr, n)

END SUBROUTINE torarray


RECURSIVE SUBROUTINE rcollect(self, array, n)
  CLASS(LIST), TARGET :: self
  INTEGER, INTENT(inout) :: n
  REAL(Realkind), POINTER, INTENT(inout)  :: array(:)
  TYPE(Element), POINTER :: current

  IF (.NOT. ASSOCIATED(self%first)) RETURN
  current => self%first
  IF (ASSOCIATED(current%sublist)) THEN
    CALL current%sublist%rcollect(array, n)
  ELSE IF (ALLOCATED(current%rvalue)) THEN
    n = n + 1
    array(n) = current%rvalue
  END IF

  DO WHILE(ASSOCIATED(current%next))
    current => current%next
    IF (ASSOCIATED(current%sublist)) THEN
      CALL current%sublist%rcollect(array, n)
    ELSE IF (ALLOCATED(current%rvalue)) THEN
      n = n + 1
      array(n) = current%rvalue
    END IF
  END DO

END SUBROUTINE rcollect


END MODULE List_Mod
