MODULE Crop_Mod

USE Kind_Mod
USE Block_Mod, Only: Block
USE MatSpRowCol_Mod, ONLY: SpRowCol
USE Subdom_Mod, ONLY: SubDomain
USE Index_Mod, ONLY : IndexC3block, IndexU3block, Indexv3block, IndexW3block

CONTAINS


SUBROUTINE make_crop_operator_cell(subdom, crop_op, allocate_mat)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(IN) :: subdom
  TYPE(SpRowCol), INTENT(INOUT) :: crop_op
  INTEGER :: iblock
  INTEGER :: i, j, k
  INTEGER :: icell_crop, icell_full
  INTEGER :: ncrop
  LOGICAL, OPTIONAL, INTENT(in) :: allocate_mat
  LOGICAL :: allocate_mat_def

  allocate_mat_def = .TRUE.
  IF (PRESENT(allocate_mat)) THEN
    allocate_mat_def = allocate_mat
  END IF

  ncrop = subdom%ncells

  crop_op%m = ncrop
  crop_op%n = subdom%ncells_tmp

  IF (allocate_mat_def) THEN
    ALLOCATE(crop_op%RowPtr(ncrop + 1))
    ALLOCATE(crop_op%ColInd(ncrop))
    ALLOCATE(crop_op%Val(ncrop))
  END IF
  crop_op%RowPtr(1) = 1

  icell_crop = 1

  DO iblock = 1, subdom%nblocks
    IF (subdom%blockiscomp(iblock)) THEN
      DO i = 1, subdom%blocks(iblock)%fld_shape(3)
        DO j = 1, subdom%blocks(iblock)%fld_shape(2)
          DO k = 1, subdom%blocks(iblock)%fld_shape(1)
            icell_full = IndexC3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
            crop_op%ColInd(icell_crop) = icell_full
            crop_op%Val(icell_crop) = one
            crop_op%RowPtr(icell_crop + 1) = icell_crop + 1
            icell_crop = icell_crop + 1
          END DO
        END DO
      END DO
    END IF
  END DO

END SUBROUTINE make_crop_operator_cell


SUBROUTINE make_crop_operator_face(subdom, crop_op, allocate_mat)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(IN) :: subdom
  TYPE(SpRowCol), INTENT(INOUT) :: crop_op
  INTEGER :: iblock
  INTEGER :: i, j, k
  INTEGER :: iface_crop, iface_full
  INTEGER :: ncrop
  LOGICAL, OPTIONAL, INTENT(in) :: allocate_mat
  LOGICAL :: allocate_mat_def

  allocate_mat_def = .TRUE.
  IF (PRESENT(allocate_mat)) THEN
    allocate_mat_def = allocate_mat
  END IF

  ncrop = subdom%nfaces
  crop_op%m = ncrop
  crop_op%n = subdom%nfaces_tmp

  IF (allocate_mat_def) THEN
    ALLOCATE(crop_op%RowPtr(ncrop + 1))
    ALLOCATE(crop_op%ColInd(ncrop))
    ALLOCATE(crop_op%Val(ncrop))
  END IF  
  crop_op%RowPtr(1) = 1

  iface_crop = 1

  DO iblock = 1, subdom%nblocks
    IF (subdom%blockiscomp(iblock)) THEN
      DO i = 1, subdom%blocks(iblock)%fld_shape(3) + 1
        DO j = 1, subdom%blocks(iblock)%fld_shape(2)
          DO k = 1, subdom%blocks(iblock)%fld_shape(1)

            iface_full = IndexU3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
            crop_op%ColInd(iface_crop) = iface_full
            crop_op%Val(iface_crop) = one
            crop_op%RowPtr(iface_crop + 1) = iface_crop + 1
            iface_crop = iface_crop + 1

          END DO
        END DO
      END DO
      DO i = 1, subdom%blocks(iblock)%fld_shape(3)
        DO j = 1, subdom%blocks(iblock)%fld_shape(2) + 1
          DO k = 1, subdom%blocks(iblock)%fld_shape(1)

            iface_full = IndexV3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
            crop_op%ColInd(iface_crop) = iface_full
            crop_op%Val(iface_crop) = one
            crop_op%RowPtr(iface_crop + 1) = iface_crop + 1
            iface_crop = iface_crop + 1

          END DO
        END DO
      END DO
      DO i = 1, subdom%blocks(iblock)%fld_shape(3)
        DO j = 1, subdom%blocks(iblock)%fld_shape(2)
          DO k = 1, subdom%blocks(iblock)%fld_shape(1) + 1

            iface_full = IndexW3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
            crop_op%ColInd(iface_crop) = iface_full
            crop_op%Val(iface_crop) = one
            crop_op%RowPtr(iface_crop + 1) = iface_crop + 1
            iface_crop = iface_crop + 1

          END DO
        END DO
      END DO
    END IF
  END DO

END SUBROUTINE make_crop_operator_face


SUBROUTINE make_crop_operator_halo_T(subdom, crop_op)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(IN) :: subdom
  TYPE(SpRowCol), INTENT(INOUT) :: crop_op
  INTEGER :: iblock
  INTEGER :: i, j, k
  INTEGER :: irow, icell_crop, icell_full
  INTEGER :: icell


  crop_op%m = subdom%ncells_tmp
  crop_op%n = subdom%nghost
  ALLOCATE(crop_op%RowPtr(subdom%ncells_tmp + 1))
  crop_op%RowPtr(1) = 1
  ALLOCATE(crop_op%ColInd(subdom%nghost))
  ALLOCATE(crop_op%Val(subdom%nghost))

  irow = 1

  icell_crop = 1
  DO iblock = 1, subdom%nblocks
    DO i = 1, SIZE(subdom%blocks(iblock)%x)
      DO j = 1, SIZE(subdom%blocks(iblock)%y)
        DO k = 1, SIZE(subdom%blocks(iblock)%z)
          icell = IndexC3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
          IF (subdom%inds1d_type(icell) == 3) THEN
            crop_op%ColInd(icell_crop) = icell_crop
            crop_op%Val(icell_crop) = 1.0
            icell_crop = icell_crop + 1
          END IF
          irow = irow + 1
          crop_op%RowPtr(irow) = icell_crop
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE make_crop_operator_halo_T


SUBROUTINE make_crop_operator_boundary(subdom, crop_op)

  IMPLICIT NONE

  TYPE(SubDomain), INTENT(in) :: subdom
  TYPE(SpRowCol), INTENT(inout) :: crop_op
!  INTEGER, POINTER, INTENT(inout) :: bnd2block_inds
  INTEGER :: iblock
  INTEGER :: i, j, k

  INTEGER :: nz, ny, nx
  INTEGER :: icell_full, icell_crop


  crop_op%m = subdom%nbound
  crop_op%n = subdom%ncells_tmp
 ! ALLOCATE(crop_op%RowPtr(subdom%nbound + 1))
  crop_op%RowPtr(1) = 1
 ! ALLOCATE(crop_op%ColInd(subdom%nbound))
 ! ALLOCATE(crop_op%Val(subdom%nbound))
  crop_op%Val(:) = 1.0

  icell_crop = 1

  DO iblock = 1, subdom%nblocks

    nz = subdom%blocks(iblock)%fld_shape(1)
    ny = subdom%blocks(iblock)%fld_shape(2)
    nx = subdom%blocks(iblock)%fld_shape(3)

    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz
          icell_full = IndexC3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
          IF (subdom%inds1d_type(icell_full) == 2) THEN
            crop_op%ColInd(icell_crop) = icell_full
            icell_crop = icell_crop + 1
            crop_op%RowPtr(icell_crop) = icell_crop
          END IF
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE make_crop_operator_boundary


SUBROUTINE make_crop_operator_halo_face_T(subdom, crop_op)

  IMPLICIT NONE

  TYPE(SubDomain), INTENT(in) :: subdom
  TYPE(SpRowCol), INTENT(inout) :: crop_op
  INTEGER :: iblock
  INTEGER :: i, j, k
  INTEGER :: nz, ny, nx
  INTEGER :: iface_full, iface_crop, icol

!  INTEGER :: nbound

!  nbound = SUM(subdom%indsface1d_type)

  crop_op%m = subdom%nfaces_tmp
  crop_op%n = SIZE(crop_op%ColInd, 1)

!  ALLOCATE(crop_op%RowPtr(subdom%nfaces_tmp + 1))
!  ALLOCATE(crop_op%ColInd(nbound))
!  ALLOCATE(crop_op%Val(nbound))
  crop_op%Val(:) = 1.0
  crop_op%RowPtr(1) = 1

  iface_full = 1
  icol = 1

  DO iblock = 1, subdom%nblocks
    nz = subdom%blocks(iblock)%fld_shape(1)
    ny = subdom%blocks(iblock)%fld_shape(2)
    nx = subdom%blocks(iblock)%fld_shape(3)

    DO i = 1, nx + 1
      DO j = 1, ny
        DO k = 1, nz
          iface_crop = IndexU3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
          IF (subdom%indsface1d_type(iface_crop) == 1) THEN
            crop_op%ColInd(icol) = icol
            icol = icol + 1            
          END IF
          iface_full = iface_full + 1
          crop_op%RowPtr(iface_full) = icol
        END DO 
      END DO
    END DO

    DO i = 1, nx
      DO j = 1, ny + 1
        DO k = 1, nz
          iface_crop = IndexV3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
          IF (subdom%indsface1d_type(iface_crop) == 1) THEN
            crop_op%ColInd(icol) = icol
            icol = icol + 1
          END IF
          iface_full = iface_full + 1
          crop_op%RowPtr(iface_full) = icol
        END DO
      END DO
    END DO

    DO i = 1, nx 
      DO j = 1, ny
        DO k = 1, nz + 1
          iface_crop = IndexW3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
          IF (subdom%indsface1d_type(iface_crop) == 1) THEN
            crop_op%ColInd(icol) = icol
            icol = icol + 1
          END IF
          iface_full = iface_full + 1
          crop_op%RowPtr(iface_full) = icol
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE make_crop_operator_halo_face_T


SUBROUTINE make_crop_operator_boundary_face(subdom, crop_op)

  IMPLICIT NONE

  TYPE(SubDomain), INTENT(in) :: subdom
  TYPE(SpRowCol), INTENT(inout) :: crop_op
  INTEGER :: iblock
  INTEGER :: i, j, k
  INTEGER :: nz, ny, nx
  INTEGER :: iface_full, iface, icol

  crop_op%m = SIZE(crop_op%RowPtr, 1) - 1
  crop_op%n = subdom%nfaces_tmp

!  ALLOCATE(crop_op%RowPtr(nbound + 1))
!  ALLOCATE(crop_op%ColInd(nbound))
!  ALLOCATE(crop_op%Val(nbound))
  crop_op%Val(:) = 1.0
  crop_op%RowPtr(1) = 1
  icol = 1
  iface_full = 1

  DO iblock = 1, subdom%nblocks

    nz = subdom%blocks(iblock)%fld_shape(1)
    ny = subdom%blocks(iblock)%fld_shape(2)
    nx = subdom%blocks(iblock)%fld_shape(3)

    DO i = 1, nx + 1
      DO j = 1, ny
        DO k = 1, nz
          iface =  IndexU3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
          IF (subdom%indsface1d_type(iface) .EQ. 1) THEN
            crop_op%ColInd(icol) = iface
            icol = icol + 1
            iface_full = iface_full + 1 
            crop_op%RowPtr(iface_full) = icol
          END IF
        END DO
      END DO
    END DO

    DO i = 1, nx
      DO j = 1, ny + 1
        DO k = 1, nz
          iface =  IndexV3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
          IF (subdom%indsface1d_type(iface) .EQ. 1) THEN
            crop_op%ColInd(icol) = iface
            icol = icol + 1
            iface_full = iface_full + 1
            crop_op%RowPtr(iface_full) = icol
          END IF
        END DO
      END DO
    END DO

    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz + 1
          iface =  IndexW3block(i, j, k, subdom%blocks, iblock, check=.TRUE.)
          IF (subdom%indsface1d_type(iface) .EQ. 1) THEN
            crop_op%ColInd(icol) = iface
            icol = icol + 1
            iface_full = iface_full + 1
            crop_op%RowPtr(iface_full) = icol
          END IF
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE make_crop_operator_boundary_face


END MODULE Crop_Mod
