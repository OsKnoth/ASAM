MODULE Laplace_sparse_Mod

USE MatSpRowCol_Mod, ONLY: SpRowCol, &
                           SpDeallocate_SpRowCol, SpNullify_SpRowCol, &
                           SpTrans_SpRowCol, SpMm_SpRowCol

USE Kind_Mod
USE Index_Mod, ONLY: vector_size, &
                     IndexU3, IndexV3, IndexW3

IMPLICIT NONE

CONTAINS

SUBROUTINE make_div3d_stencil(DIV3d_st, nz, ny, nx)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nz, ny, nx
  TYPE(SpRowCol), INTENT(INOUT) :: DIV3d_st

  INTEGER :: n, m
  INTEGER:: i, j, k
  INTEGER :: ii, irow

  m = nz * ny * nx
  n = vector_size(nx, ny, nz)

  DIV3d_st%m = m
  DIV3d_st%n = n

  ALLOCATE(DIV3d_st%RowPtr(m + 1))
  ALLOCATE(DIV3d_st%ColInd(m * 6))
  ALLOCATE(DIV3d_st%Val(m * 6))

  ii = 1
  irow = 1

  DIV3d_st%RowPtr(1) = 1
  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz 
        DIV3d_st%ColInd(ii) = IndexU3(i, j, k, nx, ny, nz)
        DIV3d_st%Val(ii) = -1.0
        ii = ii + 1
        DIV3d_st%ColInd(ii) = IndexU3(i + 1, j, k, nx, ny, nz) 
        DIV3d_st%Val(ii) = 1.0
        ii = ii + 1
        DIV3d_st%ColInd(ii) = IndexV3(i, j, k, nx, ny, nz) 
        DIV3d_st%Val(ii) = -1.0
        ii = ii + 1
        DIV3d_st%ColInd(ii) = IndexV3(i, j + 1, k, nx, ny, nz)
        DIV3d_st%Val(ii) = 1.0
        ii = ii + 1
        DIV3d_st%ColInd(ii) = IndexW3(i, j, k, nx, ny, nz)
        DIV3d_st%Val(ii) = -1.0
        ii = ii + 1
        DIV3d_st%ColInd(ii) = IndexW3(i, j, k + 1, nx, ny, nz)
        DIV3d_st%Val(ii) = 1.0
        ii = ii + 1
        DIV3d_st%RowPtr(irow + 1) = ii
        irow = irow + 1
      END DO
    END DO
  END DO

END SUBROUTINE 


SUBROUTINE make_grad_3d_stencil(DIV3d_st, GRAD3d_st)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(IN) :: DIV3d_st
  TYPE(SpRowCol), INTENT(INOUT) :: GRAD3d_st

  CALL SpTrans_SpRowCol(GRAD3d_st, DIV3d_st)
  GRAD3d_st%Val(:) = -GRAD3d_st%Val(:)

END SUBROUTINE 


SUBROUTINE make_div3d(DIV3d, DIV3d_st, ax, ay, az, v)

  IMPLICIT NONE
  
  TYPE(SpRowCol), INTENT(IN) :: DIV3d_st
  REAL(RealKind), DIMENSION(:, :, :), INTENT(IN) :: ax, ay, az, v
  TYPE(SpRowCol), INTENT(INOUT) :: DIV3d

  REAL(Realkind), POINTER :: a1d(:) 
  REAL(Realkind), POINTER :: vinv(:)

  TYPE(SpRowCol) :: A, Vi, DIV3dA

  ALLOCATE(vinv(SIZE(v)))

  vinv(:) = 1.0 / pack(v, .true.)

  CALL flatten_3dvec(a1d, ax, ay, az)
  CALL make_diag(A, a1d)
  CALL make_diag(Vi, vinv)

  DEALLOCATE(a1d)
  DEALLOCATE(vinv)

  CALL SpMm_SpRowCol(DIV3dA, DIV3d_st, A)
  CALL SpDeallocate_SpRowCol(A)  
  CALL SpMm_SpRowCol(DIV3d, Vi, DIV3dA)  
  CALL SpDeallocate_SpRowCol(Vi)
  CALL SpDeallocate_SpRowCol(DIV3dA)

END SUBROUTINE

SUBROUTINE make_grad3d(GRAD3d, GRAD3d_st, dx, dy, dz)

  IMPLICIT NONE
  
  TYPE(SpRowCol), INTENT(IN) :: GRAD3d_st
  REAL(Realkind), DIMENSION(:, :, :),  INTENT(IN) :: dx, dy, dz
  TYPE(SpRowCol), INTENT(INOUT) :: GRAD3d

  TYPE(SpRowCol) :: Dinv
  REAL(Realkind), POINTER :: d1d_inv(:)



  CALL flatten_3dvec(d1d_inv, dx, dy, dz)
  d1d_inv(:) = 1.0 / d1d_inv 
  CALL make_diag(Dinv, d1d_inv)
  DEALLOCATE(d1d_inv)
  CALL SpMm_SpRowCol(GRAD3d, Dinv, GRAD3d_st)
  CALL SpDeallocate_SpRowCol(Dinv)

END SUBROUTINE


SUBROUTINE make_lapl3d(LAPL3d, DIV3d, GRAD3d, ax, ay, az, v)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(INOUT) :: DIV3d, GRAD3d, LAPL3d
  REAL(Realkind), DIMENSION(:, :, :), INTENT(IN) :: ax, ay, az, v

  TYPE(SpRowCol) :: DIV3d_st, GRAD3d_st
  REAL(Realkind), POINTER :: dx(:, :, :), dy(:, :, :), dz(:, :, :)
  INTEGER :: shp(3) 
 
  shp(:) = SHAPE(v)

  CALL calc_dx(dx, dy, dz, ax, ay, az, v)
  CALL make_div3d_stencil(DIV3d_st, shp(1), shp(2), shp(3))
  CALL make_grad_3d_stencil(DIV3d_st, GRAD3d_st)  
  CALL make_div3d(DIV3d, DIV3d_st, ax, ay, az, v)
  CALL make_grad3d(GRAD3d, GRAD3d_st, dx, dy, dz)
  CALL SpMm_SpRowCol(LAPL3d, DIV3d_st, GRAD3d_st) 
  CALL SpDeallocate_SpRowCol(DIV3d_st)
  CALL SpDeallocate_SpRowCol(GRAD3d_st) 

  DEALLOCATE(dx)
  DEALLOCATE(dy)
  DEALLOCATE(dz)

END SUBROUTINE


SUBROUTINE flatten_3dvec(vec1d, vec_x, vec_y, vec_z)

  IMPLICIT NONE
  
  REAL(Realkind), INTENT(IN) :: vec_x(:, :, :), vec_y(:, :, :), vec_z(:, :, :)
  REAL(Realkind), POINTER, INTENT(INOUT) :: vec1d(:)
  
  INTEGER :: size_x, size_y, size_z

  size_x = SIZE(vec_x)
  size_y = SIZE(vec_y)
  size_z = SIZE(vec_z)

  ALLOCATE(vec1d(size_x + size_y + size_z))
  
  vec1d(1:size_x) = pack(vec_x, .true.)
  vec1d(size_x + 1:size_x + size_y) = pack(vec_y, .true.)
  vec1d(size_x + size_y + 1:size_x + size_y + size_z) = pack(vec_z, .true.)

END SUBROUTINE

  
SUBROUTINE make_diag(A, field1d)

  IMPLICIT NONE
  
  TYPE(SpRowCol), INTENT(INOUT) :: A
  REAL(Realkind), INTENT(IN) :: field1d(:)

  INTEGER :: n, ii

  n = SIZE(field1d)
  A%n = n
  A%m = n
  
  ALLOCATE(A%RowPtr(n + 1))
  ALLOCATE(A%ColInd(n))
  ALLOCATE(A%Val(n))

  A%RowPtr(1) = 1

  DO ii = 1, n
     A%ColInd(ii) = ii
     A%Val(ii) = field1d(ii)
     A%RowPtr(ii + 1) = ii + 1
  END DO

END SUBROUTINE


SUBROUTINE calc_dx(dx, dy, dz, ax, ay, az, v)

  IMPLICIT NONE

  REAL(Realkind), DIMENSION(:, :, :), INTENT(IN) :: ax, ay, az, v
  REAL(Realkind), POINTER, INTENT(INOUT) :: dx(:, :, :), dy(:, :, :), dz(:, :, :)

  INTEGER :: nz, ny, nx
  INTEGER :: shp(3)

  shp(:) = SHAPE(v)
  nz = shp(1)
  ny = shp(2)
  nx = shp(3)

  ALLOCATE(dx(nz, ny, nx + 1))
  ALLOCATE(dy(nz, ny + 1, nx))
  ALLOCATE(dz(nz + 1, ny, nx))

  dx(:, :, 2:nx) = 2.0 * ax(:, :, 2:nx) / (v(:, :, 1:nx - 1) + v(:, :, 2:nx))
  dx(:, :, 1) = ax(:, :, 1) / v(:, :, 1)       
  dx(:, :, nx + 1) = ax(:, :, nx + 1) / v(:, :, nx)

  dy(:, 2:ny, :) = 2.0 * ay(:, 2:ny, :) / (v(:, 1:ny - 1, :) + v(:, 2:ny, :))
  dy(:, 1, :) = ay(:, 1, :) / v(:, 1, :)
  dy(:, ny + 1, :) = ay(:, ny + 1, :) / v(:, ny, :)

  dz(2:nz, :, :) = 2.0 * az(2:nz, :, :) / (v(1:nz - 1, :, :) + v(2:nz, :, :))
  dz(1, :, :) = az(1, :, :) / v(1, :, :)
  dz(nz + 1, :, :) = az(nz + 1, :, :) / v(nz, :, :)

END SUBROUTINE


END MODULE Laplace_sparse_Mod
