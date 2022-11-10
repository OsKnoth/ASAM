MODULE Index_Mod

  USE Block_Mod, ONLY: Block
  USE Kind_Mod

  IMPLICIT NONE
  INTERFACE IndexC
    MODULE PROCEDURE IndexC3,IndexC4
  END INTERFACE


CONTAINS

FUNCTION IndexC3(ix,iy,iz,nx,ny,nz)

  INTEGER :: IndexC3
  INTEGER :: ix,iy,iz,nx,ny,nz

! SELECT CASE(GridType)
!  IndexC3=ix+((iy-1)+(iz-1)*ny)*nx 
 IndexC3=iz+((iy-1)+(ix-1)*ny)*nz 
! IndexC3=iz+((ix-1)+(iy-1)*nx)*nz 
! IndexC3=iy+((iz-1)+(ix-1)*nz)*ny 

END FUNCTION IndexC3
FUNCTION IndexC4(ip,ix,iy,iz,np,nx,ny,nz)

  INTEGER :: IndexC4
  INTEGER :: ip,ix,iy,iz,np,nx,ny,nz

!  IndexC4=ip+((ix-1)+((iy-1)+(iz-1)*ny)*nx)*np
 IndexC4=iz+((iy-1)+(ix-1)*ny)*nz 
! IndexC4=ip+((iz-1)+((ix-1)+(iy-1)*nx)*nz)*np 
! IndexC4=iy+((iz-1)+(ix-1)*nz)*ny 

END FUNCTION IndexC4

FUNCTION vector_size(nx, ny, nz)

  INTEGER :: vector_size
  INTEGER, INTENT(IN) :: nz, ny, nx

  vector_size = (nz + 1) * ny * nx + nz * (ny + 1) * nx + nz * ny * (nx + 1)

END FUNCTION vector_size

FUNCTION IndexU3(i, j, k, nx, ny, nz)

  INTEGER :: IndexU3
  INTEGER, INTENT(IN) :: i, j, k, nx, ny, nz

  IndexU3 = k + (j - 1) * nz + (i - 1) * nz * ny

END FUNCTION IndexU3


FUNCTION IndexV3(i, j, k, nx, ny, nz)

  INTEGER :: IndexV3
  INTEGER, INTENT(IN) :: i, j, k, nx, ny, nz

  IndexV3 = k + (j  - 1) * nz + (i - 1) * nz * (ny + 1) + &
            (nx + 1) * ny * nz 

END FUNCTION IndexV3

FUNCTION IndexW3(i, j, k, nx, ny, nz)

  INTEGER :: IndexW3
  INTEGER, INTENT(IN) :: i, j, k, nx, ny, nz

  IndexW3 = k + (j - 1) * (nz + 1) + (i - 1) * (nz + 1) * ny + &
            (nx + 1) * ny * nz + nx * (ny + 1) * nz 

END FUNCTION IndexW3


FUNCTION IndexC_IJK(n, nx, ny, nz)

  INTEGER :: IndexC_IJK(3)
  INTEGER, INTENT(IN) :: n, nx, ny, nz

  IndexC_IJK(1) = INT(n / (nz * ny))
  IndexC_IJK(2) = INT((n - IndexC_IJK(1) * nz * ny) / nz)
  IndexC_IJK(3) = n - IndexC_IJK(1) * nz * ny - IndexC_IJK(2) * nz

  IndexC_IJK(1) = IndexC_IJK(1) + 1
  IndexC_IJK(2) = IndexC_IJK(2) + 1
  IndexC_IJK(3) = IndexC_IJK(3) + 1  

END FUNCTION IndexC_IJK


FUNCTION IndexU_IJK(n, nx, ny, nz)

  INTEGER :: IndexU_IJK(3)
  INTEGER, INTENT(IN) :: n, nx, ny, nz

  IndexU_IJK(1) = INT((n - 1) / (nz * ny))
  IndexU_IJK(2) = INT(((n - 1) - IndexU_IJK(1) * nz * ny) / nz)
  IndexU_IJK(3) = n - IndexU_IJK(1) * nz * ny - IndexU_IJK(2) * nz

  IndexU_IJK(1) = IndexU_IJK(1) + 1
  IndexU_IJK(2) = IndexU_IJK(2) + 1

END FUNCTION IndexU_IJK


FUNCTION IndexV_IJK(n, nx, ny, nz)

  INTEGER :: IndexV_IJK(3)
  INTEGER, INTENT(IN) :: n, nx, ny, nz
  INTEGER :: m

  m = n - (nx + 1) * ny * nz - 1

  IndexV_IJK(1) = INT(m / (nz * (ny + 1)))
  IndexV_IJK(2) = INT((m - IndexV_IJK(1) * nz * (ny + 1)) / nz)
  IndexV_IJK(3) = m - IndexV_IJK(1) * nz * (ny + 1) - IndexV_IJK(2) * nz

  IndexV_IJK(1) = IndexV_IJK(1) + 1
  IndexV_IJK(2) = IndexV_IJK(2) + 1
  IndexV_IJK(3) = IndexV_IJK(3) + 1

END FUNCTION IndexV_IJK


FUNCTION IndexW_IJK(n, nx, ny, nz)

  INTEGER :: IndexW_IJK(3)
  INTEGER, INTENT(IN) :: n, nx, ny, nz
  INTEGER :: m

  m = n - (nx + 1) * ny * nz - nx * (ny + 1) * nz - 1

  IndexW_IJK(1) = INT(m / ((nz + 1) * ny + 1))
  IndexW_IJK(2) = INT((m - IndexW_IJK(1) * (nz + 1) * ny) / (nz + 1))
  IndexW_IJK(3) = m - IndexW_IJK(1) * (nz  + 1) * ny - IndexW_IJK(2) * (nz + 1)

  IndexW_IJK(1) = IndexW_IJK(1) + 1
  IndexW_IJK(2) = IndexW_IJK(2) + 1
  IndexW_IJK(3) = IndexW_IJK(3) + 1

END FUNCTION IndexW_IJK


FUNCTION IndexC3block(ix, iy, iz, blocks, iblock, check)

  INTEGER :: IndexC3block
  INTEGER, INTENT(in) :: ix, iy, iz, iblock
  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, OPTIONAL, INTENT(IN) :: check
  INTEGER :: n


  IF (PRESENT(check)) THEN
    IF (check .EQV. .TRUE.) THEN
      IF (ix .LT. 1 .OR. ix .GT. blocks(iblock)%fld_shape(3)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component x'
      END IF
      IF (iy .LT. 1 .OR. iy .GT. blocks(iblock)%fld_shape(2)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component y'
      END IF
      IF (iz .LT. 1 .OR. iz .GT. blocks(iblock)%fld_shape(1)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component z'
      END IF
    END IF
  END IF

  IndexC3block = IndexC3(ix, iy, iz, blocks(iblock)%fld_shape(3), blocks(iblock)%fld_shape(2), blocks(iblock)%fld_shape(1))

  DO n = 1, iblock - 1
    IndexC3block = IndexC3block + blocks(n)%ncells
  END DO

END FUNCTION IndexC3block


FUNCTION IndexU3block(ix, iy, iz, blocks, iblock, check)

  INTEGER :: IndexU3block
  INTEGER, INTENT(in) :: ix, iy, iz, iblock
  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, OPTIONAL, INTENT(IN) :: check
  INTEGER :: n


  IF (PRESENT(check)) THEN
    IF (check .EQV. .TRUE.) THEN
      IF (ix .LT. 1 .OR. ix .GT. blocks(iblock)%fld_shape(3) + 1) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component x'
      END IF
      IF (iy .LT. 1 .OR. iy .GT. blocks(iblock)%fld_shape(2)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component y'
      END IF
      IF (iz .LT. 1 .OR. iz .GT. blocks(iblock)%fld_shape(1)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component z'
      END IF
    END IF
  END IF

  IndexU3block = IndexU3(ix, iy, iz, blocks(iblock)%fld_shape(3), blocks(iblock)%fld_shape(2), blocks(iblock)%fld_shape(1))

  DO n = 1, iblock - 1
    IndexU3block = IndexU3block + blocks(n)%nfaces
  END DO

END FUNCTION IndexU3block

FUNCTION IndexV3block(ix, iy, iz, blocks, iblock, check)

  INTEGER :: IndexV3block
  INTEGER, INTENT(in) :: ix, iy, iz, iblock
  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, OPTIONAL, INTENT(IN) :: check
  INTEGER :: n


  IF (PRESENT(check)) THEN
    IF (check .EQV. .TRUE.) THEN
      IF (ix .LT. 1 .OR. ix .GT. blocks(iblock)%fld_shape(3)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component x'
      END IF
      IF (iy .LT. 1 .OR. iy .GT. blocks(iblock)%fld_shape(2) + 1) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component y'
      END IF
      IF (iz .LT. 1 .OR. iz .GT. blocks(iblock)%fld_shape(1)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component z'
      END IF
    END IF
  END IF

  IndexV3block = IndexV3(ix, iy, iz, blocks(iblock)%fld_shape(3), blocks(iblock)%fld_shape(2), blocks(iblock)%fld_shape(1))

  DO n = 1, iblock - 1
    IndexV3block = IndexV3block + blocks(n)%nfaces
  END DO

END FUNCTION IndexV3block


FUNCTION IndexW3block(ix, iy, iz, blocks, iblock, check)

  INTEGER :: IndexW3block
  INTEGER, INTENT(in) :: ix, iy, iz, iblock
  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, OPTIONAL, INTENT(IN) :: check
  INTEGER :: n


  IF (PRESENT(check)) THEN
    IF (check .EQV. .TRUE.) THEN
      IF (ix .LT. 1 .OR. ix .GT. blocks(iblock)%fld_shape(3)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component x'
      END IF
      IF (iy .LT. 1 .OR. iy .GT. blocks(iblock)%fld_shape(2)) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component y'
      END IF
      IF (iz .LT. 1 .OR. iz .GT. blocks(iblock)%fld_shape(1) + 1) THEN
        WRITE(*,*) 'ERROR IN IndexC3block component z'
      END IF
    END IF
  END IF

  IndexW3block = IndexW3(ix, iy, iz, blocks(iblock)%fld_shape(3), blocks(iblock)%fld_shape(2), blocks(iblock)%fld_shape(1))

  DO n = 1, iblock - 1
    IndexW3block = IndexW3block + blocks(n)%nfaces
  END DO

END FUNCTION IndexW3block


FUNCTION coord_ilbnd(coords2, lbnd_coord, nghost)
  IMPLICIT NONE

  INTEGER :: coord_ilbnd
  REAL(Realkind), INTENT(in) :: coords2(:)
  REAL(Realkind), INTENT(in) :: lbnd_coord
  INTEGER, INTENT(in) :: nghost
  INTEGER :: i

  DO i = 1, SIZE(coords2)
    IF (coords2(i) .GE. lbnd_coord) THEN
      coord_ilbnd = i
      EXIT
    END IF
  END DO

  !Add one extra ghost layer for communication
  coord_ilbnd = MAX(1, coord_ilbnd - nghost)
END FUNCTION coord_ilbnd


FUNCTION coord_iubnd(coords2, ubnd_coord, nghost)
  IMPLICIT NONE

  INTEGER :: coord_iubnd
  REAL(Realkind), INTENT(in) :: coords2(:)
  REAL(Realkind), INTENT(in) :: ubnd_coord
  INTEGER, INTENT(in) :: nghost
  INTEGER :: i

  DO i = SIZE(coords2), 1, -1
    IF (coords2(i) .LE. ubnd_coord) THEN
       coord_iubnd = i
       EXIT
    END IF
  END DO

  !Add one extra ghost layer for communication
 coord_iubnd = MIN(SIZE(coords2), coord_iubnd + nghost)

END FUNCTION coord_iubnd


INTEGER FUNCTION findloc_real(arr, val)

  IMPLICIT NONE
  REAL(Realkind), INTENT(in) :: arr(:)
  REAL(Realkind), INTENT(in) :: val

  INTEGER :: l_bnd, u_bnd, i

  l_bnd = LBOUND(arr, dim=1)
  u_bnd = UBOUND(arr, dim=1)  

  DO i = l_bnd, u_bnd
    IF (arr(i) .EQ. val) THEN
      findloc_real = i
      EXIT
    END IF
  END DO

END FUNCTION findloc_real

INTEGER FUNCTION findloc_integer(arr, val)

  IMPLICIT NONE
  INTEGER, INTENT(in) :: arr(:)
  INTEGER, INTENT(in) :: val

  INTEGER :: l_bnd, u_bnd, i

  l_bnd = LBOUND(arr, dim=1)
  u_bnd = UBOUND(arr, dim=1)

  DO i = l_bnd, u_bnd
    IF (arr(i) .EQ. val) THEN
      findloc_integer = i
      EXIT
    END IF
  END DO

END FUNCTION findloc_integer

FUNCTION between(arr, value) Result(i)

  REAL(Realkind), POINTER, INTENT(in) :: arr(:)
  REAL(Realkind), INTENT(in) :: value

  INTEGER :: i, lbnd, ubnd

  lbnd = LBOUND(arr, 1)
  ubnd = UBOUND(arr, 1) - 1
 
  DO i = lbnd, ubnd
    IF (arr(i) .LE. value .AND. arr(i + 1) .GE. value) EXIT
  END DO

END FUNCTION between

SUBROUTINE determine_overlap_1d(arr1, arr2, o12, o21)
  IMPLICIT NONE

  REAL(Realkind), INTENT(in) :: arr1(:), arr2(:)
  INTEGER, DIMENSION(2), INTENT(inout) :: o12, o21
  INTEGER :: n1, n2

  o12(:) = 0
  o21(:) = 0

  n1 = UBOUND(arr1, 1)
  n2 = UBOUND(arr2, 1)

  IF (arr1(1) .GE. arr2(n2) .OR. arr1(n1) .LE. arr2(1)) RETURN

  IF (arr1(1) .GE. arr2(1) .AND. arr1(1) .LT. arr2(n2)) THEN
    o21(1) = FINDLOC(arr2 .LE. arr1(1), .TRUE., BACK=.TRUE., DIM=1)
    o12(1) = 1
  ELSE
    o12(1) = FINDLOC(arr1, arr2(1), DIM=1)
    o21(1) = 1
  END IF

  IF (arr1(n1) .GT. arr2(1) .AND. arr1(n1) .LE. arr2(n2)) THEN
    o21(2) = MAX(FINDLOC(arr2, arr1(n1), DIM=1) - 1, 1)
    o12(2) = n1 - 1
  ELSE
    o12(2) = MAX(FINDLOC(arr1, arr2(n2), DIM=1) - 1, 1)
    o21(2) = n2 - 1
  END IF

END SUBROUTINE determine_overlap_1d


SUBROUTINE determine_overlap_1d_halo(arr1, arr2, o12, o21)
  IMPLICIT NONE

  REAL(Realkind), POINTER, INTENT(in) :: arr1(:), arr2(:)
  INTEGER, INTENT(inout) :: o12(2), o21(2)
  INTEGER :: n1, n2
  REAL(Realkind) :: minval, maxval

  o12(:) = 0.0
  o21(:) = 0.0

  n1 = UBOUND(arr1, 1)
  n2 = UBOUND(arr2, 1)

  IF (arr1(1) .GT. arr2(n2) .OR. arr1(n1) .LT. arr2(1)) RETURN

  maxval = MAX(arr1(1), arr2(1))
  o12(1) = MAX(1, FINDLOC(arr1, maxval, DIM=1) - 1)
  IF (arr2(1) .EQ. maxval) THEN
    o21(1) = 1
  ELSE
    o21(1) = FINDLOC(arr2 .LT. maxval, .TRUE., DIM=1, BACK=.True.)
  END IF

  minval = MIN(arr1(n1), arr2(n2))
  o12(2) = MIN(n1 - 1, FINDLOC(arr1, minval, DIM=1))
  IF (arr2(n2) .EQ. minval) THEN
    o21(2) = n2 - 1
  ELSE
    o21(2) = MAX(FINDLOC(arr2 .GT. minval, .True., DIM=1) - 1, 1)
  END IF


END SUBROUTINE determine_overlap_1d_halo


END MODULE Index_Mod
