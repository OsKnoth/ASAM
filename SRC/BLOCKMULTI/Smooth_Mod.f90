MODULE Smooth_Mod

USE MatSpRowCol_Mod, ONLY: SpRowCol, SpDeallocate_SpRowCol, &
                           SpNullify_SpRowCol, SpMm_SpRowCol, &
                           SpAVec_SpRowCol, Resize_sprowcol_dataflds
USE Kind_Mod

USE Boundary_Mod, ONLY: BoundCommInterface, matmul_boundary_condition, &
                        matmul_rd_boundary_condition, matmul_bl_boundary_condition, &
                        init_two_sided_exchange, finish_two_sided_exchange

USE Block_Mod, ONLY: block


IMPLICIT NONE

TYPE :: SorSolver
  TYPE(SpRowCol), POINTER :: mat => NULL()
  TYPE(SpRowCol), POINTER :: mat_rd => NULL()
  TYPE(SpRowCol), POINTER :: mat_bl => NULL()
  TYPE(BoundCommInterface), POINTER :: comm_interf => NULL()
  REAL(Realkind), POINTER :: diaginv(:) => NULL()
  REAL(Realkind), POINTER :: diaginv_rd(:) => NULL()
  REAL(Realkind), POINTER :: diaginv_bl(:) => NULL()

  CONTAINS

  GENERIC :: init_from_sprowcol => &
             init_from_sprowcol_lexicographic, &
             init_from_sprowcol_redblack
  PROCEDURE ::                                   &
               init_from_sprowcol_lexicographic, &
               init_from_sprowcol_redblack
  PROCEDURE :: smooth
  PROCEDURE :: free
END TYPE SorSolver

CONTAINS

SUBROUTINE init_from_sprowcol_lexicographic(self, mat, omega, comm_interf)

  IMPLICIT NONE

  CLASS(SorSolver) self
  TYPE(SpRowCol), TARGET, INTENT(IN) :: mat
  REAL(Realkind), INTENT(IN) :: omega
  TYPE(BoundCommInterface), TARGET, OPTIONAL, INTENT(inout) :: comm_interf 
  INTEGER :: i

  self%mat => mat


  ALLOCATE(self%diaginv(mat%m))

  DO i = 1, mat%m
    IF (mat%Diagptr(i) == -1) THEN
      self%diaginv(i) = 0.0
    ELSE
      self%diaginv(i) = omega / mat%Val(mat%Diagptr(i))
    END IF
  END DO


  IF (PRESENT(comm_interf)) self%comm_interf => comm_interf

END SUBROUTINE init_from_sprowcol_lexicographic

SUBROUTINE init_from_sprowcol_redblack(self, mat,  Bl_crop, Rd_crop, omega, comm_interf)

  IMPLICIT NONE

  CLASS(SorSolver) self
  TYPE(SpRowCol), TARGET, INTENT(IN) :: mat, Bl_crop, Rd_crop
  REAL(Realkind), INTENT(IN) :: omega
  TYPE(BoundCommInterface), TARGET, OPTIONAL, INTENT(in) :: comm_interf 
  TYPE(SpRowCol), TARGET, SAVE :: mat_rd, mat_bl  
  INTEGER :: i

  CALL SpMm_SpRowCol(mat_rd, Rd_crop, mat)
  CALL SpMm_SpRowCol(mat_bl, Bl_crop, mat)

  self%mat_rd => mat_rd
  self%mat_bl => mat_bl

  ALLOCATE(self%diaginv_rd(mat_rd%m))
  ALLOCATE(self%diaginv_bl(mat_rd%m))

  self%diaginv_rd(:) = 0.0
  self%diaginv_bl(:) = 0.0

  DO i = 1, mat_rd%m
    IF (mat_rd%Diagptr(i) .GT. 0) THEN
      self%diaginv_rd(i) = omega / mat_rd%Val(mat_rd%Diagptr(i))
    END IF
  END DO

  DO i = 1, mat_bl%m
    IF (mat_bl%Diagptr(i) .GT. 0) THEN
      self%diaginv_bl(i) = omega / mat_bl%Val(mat_bl%Diagptr(i))
    END IF
  END DO

  IF (PRESENT(comm_interf)) self%comm_interf => comm_interf

END SUBROUTINE init_from_sprowcol_redblack


SUBROUTINE make_red_black(blocks, ncells_tmp, blockids_compghst, isblackijkzero, Bl, Rd)

  IMPLICIT NONE

  TYPE(Block), POINTER, INTENT(in) :: blocks(:)
  INTEGER, INTENT(in) :: ncells_tmp
  INTEGER, INTENT(in) :: blockids_compghst(:), isblackijkzero(:)

  TYPE(SpRowCol), INTENT(inout) :: Bl, Rd

  INTEGER :: nz, ny, nx, iblock, nblocks
  INTEGER :: i, j, k, ii
  INTEGER :: icolRd, icolBl
  INTEGER :: i_sign, j_sign, k_sign, ii_sign


  IF (.NOT. ASSOCIATED(blocks)) THEN
    nblocks = 0
  ELSE
    nblocks = SIZE(blocks, 1)
  END IF


  i_sign = 1
  j_sign = 1
  k_sign = -1

  icolRd = 1
  icolBl = 1

  Bl%n = ncells_tmp
  Bl%m = ncells_tmp

  Rd%n = ncells_tmp
  Rd%m = ncells_tmp

  Bl%RowPtr(1) = 1
  Rd%RowPtr(1) = 1

  ii = 1
  icolBl = 1
  icolRd = 1

  DO iblock = 1, nblocks

    nz = blocks(iblock)%fld_shape(1)
    ny = blocks(iblock)%fld_shape(2)
    nx = blocks(iblock)%fld_shape(3)

    i_sign = 1
    j_sign = 1

    IF (MOD(nz, 2) == 0) j_sign = -1
    IF (MOD(ny, 2) == 0) i_sign = -1

    ii_sign = -1 + 2 * isblackijkzero(blockids_compghst(iblock))

    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz
          IF (ii_sign > 0) THEN
            Bl%ColInd(icolBl) = ii
            icolBl = icolBl + 1
          ELSE
            Rd%ColInd(icolRd) = ii
            icolRd = icolRd + 1
          END IF
          ii_sign = ii_sign * k_sign
          ii = ii + 1
          Rd%RowPtr(ii) = icolRd
          Bl%RowPtr(ii) = icolBl
        END DO
        ii_sign = ii_sign * j_sign
!        IF (ny .GT. 1)  ii_sign = ii_sign * j_sign
      END DO
      ii_sign = ii_sign * i_sign
!      IF (nx .GT. 1) ii_sign = ii_sign * i_sign
    END DO
  END DO

END SUBROUTINE make_red_black


SUBROUTINE smooth(self, x, b, Numiter)

  IMPLICIT NONE

  CLASS(SorSolver) :: self
  REAL(Realkind), INTENT(IN) :: b(:)
  REAL(Realkind), INTENT(INOUT) :: x(:)
  INTEGER, INTENT(IN) :: NumIter
  REAL(Realkind) :: val
  REAL(Realkind) :: field_tmp(SIZE(x))
  INTEGER :: i, j, iter

  IF (ASSOCIATED(self%mat)) THEN
    DO iter = 1, NumIter
      field_tmp(:) = b
      IF (ASSOCIATED(self%comm_interf)) THEN
        CALL init_two_sided_exchange(self%comm_interf, x)
        CALL finish_two_sided_exchange(self%comm_interf)
        CALL matmul_boundary_condition(self%comm_interf, field_tmp, mode="sub")
      END IF
      DO i = 1, self%mat%m
        val = field_tmp(i)
        DO j = self%mat%RowPtr(i), self%mat%RowPtr(i + 1) - 1
            val = val - self%mat%Val(j) * x(self%mat%ColInd(j))
        END DO
        x(i) = x(i) + val * self%diaginv(i)
      END DO
    END DO

  ELSE IF (ASSOCIATED(self%mat_rd) .AND. ASSOCIATED(self%mat_bl)) THEN

   IF (ASSOCIATED(self%comm_interf)) THEN
      DO iter = 1, NumIter
        CALL init_two_sided_exchange(self%comm_interf, x)
        CALL SpAVec_SpRowCol(field_tmp, self%mat_rd, x)
        CALL finish_two_sided_exchange(self%comm_interf)
        CALL matmul_rd_boundary_condition(self%comm_interf, field_tmp)
        x(:) = x + (b - field_tmp) * self%diaginv_rd

        CALL init_two_sided_exchange(self%comm_interf, x)
        CALL SpAVec_SpRowCol(field_tmp, self%mat_bl, x)
        CALL finish_two_sided_exchange(self%comm_interf)
        CALL matmul_bl_boundary_condition(self%comm_interf, field_tmp)
        x(:) = x + (b - field_tmp) * self%diaginv_bl
      END DO
    ELSE
      DO iter = 1, NumIter
        CALL SpAVec_SpRowCol(field_tmp, self%mat_rd, x)
        x(:) = x + (b - field_tmp) * self%diaginv_rd
        CALL SpAVec_SpRowCol(field_tmp, self%mat_bl, x)
        x(:) = x + (b - field_tmp) * self%diaginv_bl
      END DO
    END IF
  END IF

END SUBROUTINE smooth


SUBROUTINE free(self)

  IMPLICIT NONE

  CLASS(SorSolver) :: self

  IF (ASSOCIATED(self%mat)) THEN
    self%mat => NULL()
  END IF
  IF (ASSOCIATED(self%mat_rd)) THEN
    self%mat_rd => NULL()
  END IF
  IF (ASSOCIATED(self%mat_bl)) THEN
    self%mat_bl => NULL()
  END IF
  IF (ASSOCIATED(self%comm_interf)) self%comm_interf => NULL()

END SUBROUTINE free

END MODULE Smooth_Mod
