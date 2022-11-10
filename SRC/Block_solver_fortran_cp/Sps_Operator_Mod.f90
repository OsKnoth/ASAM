MODULE Sps_Operators_Mod

USE MatSpRowCol_Mod, ONLY: SpRowCol
USE Kind_Mod
USE Index_Mod, ONLY: vector_size, IndexC_IJK, &
                     IndexU_IJK, IndexV_IJK, IndexW_IJK, &
                     IndexC3, IndexU3, IndexV3, IndexW3

CONTAINS

SUBROUTINE sps_remove_halo_cell(ops_vol, nz, ny, nx, nhalo)
  !Constructs a sparse operator to remove nhalo
  !layers from a 3d cell-centred field

  IMPLICIT NONE
  
  TYPE(SpRowCol), INTENT(INOUT) :: ops_cell
  INTEGER, INTENT(IN) :: nz, ny, nx
  INTEGER, INTENT(IN) :: nhalo
  INTEGER :: nz_r, ny_r, nx_r
 
  INTEGER :: size_cell_halo, size_cell

  INTEGER :: ind_vec(3)
  INTEGER :: n, m

  nz_r = nz - 2 * nhalo
  ny_r = ny - 2 * nhalo
  nx_r = nx - 2 * nhalo

  size_cell_halo = nz * ny * nx
  size_cell = nz_r * ny_r * nx_r

  ops_cell%n = size_cell
  ops_cell%m = size_cell_halo
  ALLOCATE(ops_cell%Rowptr(size_cell + 1))
  ALLOCATE(ops_cell%Colind(size_cell))
  ALLOCATE(ops_cell%Val(size_cell))

  ops_cell%Rowptr(1) = 1
  DO n = 1, size_cell
    ind_vec =  IndexC_IJK(n, nx_r, ny_r, nz_r)
    ind_vec(1) = ind_vec(1) + nhalo
    ind_vec(2) = ind_vec(2) + nhalo
    ind_vec(3) = ind_vec(3) + nhalo
    m = IndexC3(ind_vec(1), ind_vec(2), ind_vec(3), nx, ny, nz)
    ops_cell%Colind(n) = m
    ops_cell%Val(n) = 1.0
    ops_cell%Rowptr(n + 1) = n + 1
  END DO

END SUBROUTINE sps_remove_halo_cell


SUBROUTINE sps_remove_halo_face(ops_face nz, ny, nx, nhalo)
  !Constructs a sparse operator to remove nhalo
  !layers from a 3d face-centred vector field

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(INOUT) :: ops_face
  INTEGER, INTENT(IN) :: nz, ny, nx
  INTEGER, INTENT(IN) :: nhalo
  INTEGER :: nz_r, ny_r, nx_r

  INTEGER :: size_face_halo, size_face
  INTEGER :: size_u, size_v, size_w

  INTEGER :: ind_vec(3)
  INTEGER :: n, m

  size_face_halo = vector_size(nx, ny, nz)
  size_face = vector_size(nx_r, ny_r, nz_r)
  size_u = (nx_r + 1) * ny_r  * nz_r
  size_v = nx_r * (ny_r + 1) * nz_r
  size_w = nx_r * ny_r * (nz_r + 1)

  ops_face%n = size_face
  ops_face%m = size_face_halo
  ALLOCATE(ops_face%Rowptr(size_face + 1))
  ALLOCATE(ops_face%Colind(size_face))
  ALLOCATE(ops_face%Val(size_face)) 

  ops_face%Rowptr(1) = 1
  DO n = 1, size_u
    ind_vec = IndexU3(n, nx_r, ny_r, nz_r)
    ind_vec(1) = ind_vec(1) + nhalo
    ind_vec(2) = ind_vec(2) + nhalo
    ind_vec(3) = ind_vec(3) + nhalo
    m = IndexU3(ind_vec(1), ind_vec(2), ind_vec(3), nx, ny, nz)
    ops_cell%Colind(n) = m
    ops_cell%Val(n) = 1.0
    ops_cell%Rowptr(n + 1) = n + 1
  END DO
  
  DO n = size_u + 1, size_u + size_v
    ind_vec = IndexV3(n, nx_r, ny_r, nz_r)
    ind_vec(1) = ind_vec(1) + nhalo
    ind_vec(2) = ind_vec(2) + nhalo
    ind_vec(3) = ind_vec(3) + nhalo
    m = IndexV3(ind_vec(1), ind_vec(2), ind_vec(3), nx, ny, nz)
    ops_cell%Colind(n) = m
    ops_cell%Val(n) = 1.0
    ops_cell%Rowptr(n + 1) = n + 1
  END DO

  DO n = size_u + size_v + 1, size_face
    ind_vec = IndexW3(n, nx_r, ny_r, nz_r)
    ind_vec(1) = ind_vec(1) + nhalo
    ind_vec(2) = ind_vec(2) + nhalo
    ind_vec(3) = ind_vec(3) + nhalo
    m = IndexW3(ind_vec(1), ind_vec(2), ind_vec(3), nx, ny, nz)
    ops_cell%Colind(n) = m
    ops_cell%Val(n) = 1.0
    ops_cell%Rowptr(n + 1) = n + 1
  END DO

END SUBROUTINE sps_remove_halo_face


SUBROUTINE sps_insert(ops_insert, field1d)
  ! Constructs a sparse operator to insert a 
  ! smaller cell-centred field into a larger cell-centred
  ! field

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(INOUT) :: ops_insert
  REAL(Realkind), INTENT(IN) :: field1d(:)

  INTEGER :: nnz, field_size

  INTEGER :: n, m

  nnz = 0
  field_size = SIZE(field1d)

  DO n = 1, field_size
    IF (field1d(n) .NOT. /= 0.0) nnz = nnz + 1 
  END DO

  ops_insert%n = field_size
  ops_insert%m = nnz
  ALLOCATE(ops_insert%Rowptr(field_size + 1))
  ALLOCATE(ops_insert%Colind(nnz))
  ALLOCATE(ops_insert%Val(nnz))

  m = 1
  ops_insert%Rowptr(1) = 1
  DO n = 1, field_size
    IF (field1d(n) .NOT. /= 0.0) THEN
       ops_insert%Colind(m) = n
       ops_insert%Val(m) = 1.0      
       m = m + 1
    END IF
    ops_insert%Rowptr(n + 1) = m
  END DO

END SUBROUTINE

SUBROUTINE make_eye(eye_op, n):

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(INOUT) :: eye_op
  INTEGER, INTENT(in) :: n
  INTEGER :: i

  eye_op%n = n
  eye_op%m = n

  ALLOCATE(eye_op%RowPtr(n + 1))
  ALLOCATE(eye_op%ColInd(n))
  ALLOCATE(eye_op%Val(n))

  eye_op%Val(:) = 1.0

  eye_op%RowPtr(1) = 1

  DO i = 1, n
    eye_op%ColInd(i) = i
    eye_op%RowPtr(i + 1) = i + 1
  END DO

END SUBROUTINE make_eye

END MODULE Sps_Operators_Mod
