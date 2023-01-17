MODULE Gradient_block_Mod

USE MatSpRowCol_Mod, ONLY: SpRowCol, &
                           SpDeallocate_SpRowCol, SpNullify_SpRowCol, &
                           SpTrans_SpRowCol, SpMm_SpRowCol, SpA_SpRowCol, &
                           resize_sprowcol_dataflds
USE Block_Mod, ONLY: Block
USE Kind_Mod
USE Index_Mod, ONLY: vector_size, &
                     IndexU3, IndexV3, IndexW3, IndexC3, &
                     IndexC3block, IndexU3block, IndexV3block, IndexW3block

USE Gradient_block_interpolation_Mod, ONLY: boundint_weights_yz, boundint_weights_xz, boundint_weights_xy, &
                                            find_indices_jk, find_indices_ik, find_indices_ij


IMPLICIT NONE

CONTAINS



SUBROUTINE make_grad3d(GRAD3d, blocks)
  !This routine creates the gradient operator for 
  !a 3d-field decomposed into unstructured blocks (which are
  !internally structured).
  !The blocks can have different horizontal refinement levels.
  !Input GRAD3d must have pre-allocated arrays

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: GRAD3d
  TYPE(Block), POINTER, INTENT(in) :: blocks(:)

  TYPE(SpRowCol) :: INT2d_bnds
  TYPE(SpRowCol) :: GRAD3d_a
  TYPE(SpRowCol) :: GRAD3d_b

  INTEGER, TARGET :: GRAD3d_a_RowPtr(GRAD3d%m + 1), &
                     GRAD3d_a_ColInd(SIZE(GRAD3d%ColInd)), &
                     GRAD3d_b_RowPtr(GRAD3d%m + 1), &
                     GRAD3d_b_ColInd(ubound_ncolids_int2d(blocks)), &
                     INT2d_bnds_RowPtr(GRAD3d%m + 1), &
                     INT2d_bnds_ColInd(ubound_ncolids_int2d(blocks))
                     

  REAL(Realkind), TARGET :: GRAD3d_a_Val(SIZE(GRAD3d%Val)), &
                            GRAD3d_b_Val(ubound_ncolids_int2d(blocks)), &
                            INT2d_bnds_Val(ubound_ncolids_int2d(blocks))

  INTEGER :: iblock, nblocks

  INTEGER :: k, nz
  INTEGER :: j, ny
  INTEGER :: i, nx

  INTEGER :: n, n2, n_tmp
  INTEGER :: irow, ist_cell
  INTEGER :: ucface_ncolids

  INTEGER :: ii, jj, kk, l

  INTEGER :: i_bnds(2), j_bnds(2), k_bnds(2)
  INTEGER :: bndblock_id

  REAL(Realkind), PARAMETER :: eps = 1e-40
  REAL :: area, area_fac, volume, gcoeff, gcoeff_sum

  REAL(Realkind), POINTER :: volseff(:,:,:), arseffx(:,:,:), arseffy(:,:,:), arseffz(:,:,:)
  REAL(Realkind), POINTER :: vols_bnd(:,:,:), volseff_bnd(:,:,:)
  REAL(Realkind), POINTER :: arsx_bnd(:,:,:), arseffx_bnd(:,:,:)
  REAL(Realkind), POINTER :: arsy_bnd(:,:,:), arseffy_bnd(:,:,:)
  REAL(Realkind), POINTER :: arsz_bnd(:,:,:), arseffz_bnd(:,:,:)

  REAL(Realkind), POINTER :: x(:), y(:), z(:)
  REAL(Realkind), POINTER :: x_bnd(:), y_bnd(:), z_bnd(:)
  REAL(Realkind), POINTER :: x2_bnd(:), y2_bnd(:), z2_bnd(:)

  REAL(Realkind) :: weights(3), point_dest(2), points_volsarea(3, 2)
  INTEGER :: points_loc(3, 4)

  INTEGER :: ipoint

  INTEGER ::  nz_bnd, ny_bnd, nx_bnd
  INTEGER :: uface_int2d_ncolids

  IF (.NOT. ASSOCIATED(blocks)) THEN
    nblocks = 0
  ELSE
    nblocks = SIZE(blocks)
  END IF  

  GRAD3d_a%m = GRAD3d%m
  GRAD3d_a%n = GRAD3d%n

  INT2d_bnds%m = GRAD3d%m
  INT2d_bnds%n = GRAD3d%m

!  ucface_ncolids = ubound_ncolids_blockgrad(blocks)
!  uface_int2d_ncolids= ubound_ncolids_int2d(blocks)

  GRAD3d_a%RowPtr => GRAD3d_a_RowPtr
  GRAD3d_a%ColInd => GRAD3d_a_ColInd
  GRAD3d_a%Val => GRAD3d_a_Val

  INT2d_bnds%RowPtr => INT2d_bnds_RowPtr
  INT2d_bnds%ColInd => INT2d_bnds_ColInd
  INT2d_bnds%Val => INT2d_bnds_Val

  GRAD3d_a%ColInd(:) = -1
  GRAD3d_a%RowPtr(:) = -1
  GRAD3d_a%Val(:) = -999.9
  GRAD3d_a%RowPtr(1) = 1

!  ALLOCATE(INT2d_bnds%RowPtr(INT2d_bnds%m + 1))
!  ALLOCATE(INT2d_bnds%ColInd(uface_int2d_ncolids))
!  ALLOCATE(INT2d_bnds%Val(uface_int2d_ncolids))

  INT2d_bnds%ColInd(:) = -1
  INT2d_bnds%RowPtr(:) = -1
  INT2d_bnds%Val(:) = -999.9
  INT2d_bnds%RowPtr(1) = 1

  GRAD3d_b%RowPtr => GRAD3d_b_RowPtr
  GRAD3d_b%ColInd => GRAD3d_b_ColInd
  GRAD3d_b%Val => GRAD3d_b_Val

  n = 1
  n2 = 1
  irow = 1

  ist_cell = 0

  DO iblock = 1, nblocks

    nz = blocks(iblock)%fld_shape(1)
    ny = blocks(iblock)%fld_shape(2)
    nx = blocks(iblock)%fld_shape(3)
 
    x => blocks(iblock)%x
    y => blocks(iblock)%y
    z => blocks(iblock)%z

    volseff => blocks(iblock)%volseff
    arseffx => blocks(iblock)%arseffx
    arseffy => blocks(iblock)%arseffy
    arseffz => blocks(iblock)%arseffz

    ! x-component
    DO i = 1, nx + 1

      IF (i == 1) THEN

        IF (ANY(blocks(iblock)%cface_w)) THEN

          DO j = 1, ny
            DO k = 1, nz

              CALL find_indices_jk(blocks, blocks(iblock)%cface_w, j, j_bnds, k, k_bnds, iblock, bndblock_id)

              IF (bndblock_id .GT. 99999999) THEN
                GRAD3d_a%RowPtr(irow + 1) = n
                INT2d_bnds%RowPtr(irow + 1) = n2
                irow = irow + 1              
                CYCLE
              END IF

              gcoeff_sum = 0.0

              vols_bnd => blocks(bndblock_id)%vols
              arsx_bnd => blocks(bndblock_id)%arsx
              volseff_bnd => blocks(bndblock_id)%volseff
              arseffx_bnd => blocks(bndblock_id)%arseffx

              y_bnd => blocks(bndblock_id)%y
              z_bnd => blocks(bndblock_id)%z
              y2_bnd => blocks(bndblock_id)%y2
              z2_bnd => blocks(bndblock_id)%z2
              nz_bnd = blocks(bndblock_id)%fld_shape(1)
              ny_bnd = blocks(bndblock_id)%fld_shape(2)
              nx_bnd = blocks(bndblock_id)%fld_shape(3)


              IF ((k_bnds(2) - k_bnds(1) > 0) .OR. (j_bnds(2) - j_bnds(1) > 0)) THEN
                DO jj = j_bnds(1), j_bnds(2)
                  DO kk = k_bnds(1), k_bnds(2)
                    gcoeff = 0.0
                    INT2d_bnds%Val(n2) = arseffx_bnd(kk, jj, nx_bnd  + 1) / (arseffx(k, j, i) + eps)
                    INT2d_bnds%ColInd(n2) = IndexU3block(nx_bnd + 1, jj, kk, blocks, bndblock_id, check=.TRUE.)
                    n2 = n2 + 1                       
                  END DO
                END DO
              ELSE
                IF (ABS(z_bnd(k_bnds(1)) - z(k)) + ABS(y_bnd(j_bnds(1)) - y(j)) > eps) THEN
                  CALL boundint_weights_yz(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                                           x, y, z, y_bnd, z_bnd, y2_bnd, z2_bnd, volseff_bnd, arseffx_bnd, &
                                           j, k,  j_bnds, k_bnds, -1, point_dest)
                  DO l = 1, ipoint
                    gcoeff = 2.0 * weights(l) * arseffx(k, j, i) * points_volsarea(l, 2) / &
                             (volseff(k, j, i) * points_volsarea(l, 2) + points_volsarea(l, 1) * arseffx(k, j, i) + eps)
                    gcoeff_sum = gcoeff_sum + gcoeff
                    GRAD3d_a%ColInd(n) = IndexC3block(points_loc(l, 2), points_loc(l, 3), points_loc(l, 4), &
                                                    blocks, points_loc(l, 1), check=.TRUE.)
                    GRAD3d_a%Val(n) = -gcoeff
                    n = n + 1
                  END DO
                ELSE
                  volume = 0.5 * (volseff(k, j, i) + volseff_bnd(k_bnds(1), j_bnds(1), nx_bnd))
                  gcoeff =  arseffx(k, j, i) / volume
                  gcoeff_sum = gcoeff
                  GRAD3d_a%ColInd(n) = IndexC3block(SIZE(blocks(bndblock_id)%x2) - 1, j_bnds(1), k_bnds(1), &
                                                  blocks, bndblock_id, check=.TRUE.)
                  GRAD3d_a%Val(n) = -gcoeff
                  n = n + 1
                END IF
              END IF

              GRAD3d_a%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
              GRAD3d_a%Val(n) = gcoeff_sum
              n = n + 1

              GRAD3d_a%RowPtr(irow + 1) = n
              INT2d_bnds%RowPtr(irow + 1) = n2
              irow = irow + 1
            END DO
          END DO
        ELSE
          DO j = 1, ny
            DO k = 1, nz
              GRAD3d_a%RowPtr(irow + 1) = n
              INT2d_bnds%RowPtr(irow + 1) = n2
              irow = irow + 1
            END DO
          END DO
        END IF

      ELSEIF (i == nx + 1) THEN

        IF (ANY(blocks(iblock)%cface_e)) THEN
          DO j = 1, ny
            DO k = 1, nz

              CALL find_indices_jk(blocks, blocks(iblock)%cface_e, j, j_bnds, k, k_bnds, iblock, bndblock_id)

              IF (bndblock_id .GT. 99999999) THEN
                GRAD3d_a%RowPtr(irow + 1) = n
                INT2d_bnds%RowPtr(irow + 1) = n2
                irow = irow + 1
                CYCLE
              END IF

              vols_bnd => blocks(bndblock_id)%vols
              arsx_bnd => blocks(bndblock_id)%arsx
              volseff_bnd => blocks(bndblock_id)%volseff
              arseffx_bnd => blocks(bndblock_id)%arseffx

              y_bnd => blocks(bndblock_id)%y
              z_bnd => blocks(bndblock_id)%z
              y2_bnd => blocks(bndblock_id)%y2
              z2_bnd => blocks(bndblock_id)%z2
              nz_bnd = blocks(bndblock_id)%fld_shape(1)
              ny_bnd = blocks(bndblock_id)%fld_shape(2)
              nx_bnd = blocks(bndblock_id)%fld_shape(3)

              GRAD3d_a%ColInd(n) = IndexC3block(i - 1, j, k, blocks, iblock, check=.TRUE.)

              n_tmp = n
              n = n + 1

              gcoeff_sum = 0.0              

              IF ((k_bnds(2) - k_bnds(1) > 0) .OR. (j_bnds(2) - j_bnds(1) > 0)) THEN
                DO jj = j_bnds(1), j_bnds(2)
                  DO kk = k_bnds(1), k_bnds(2)
                    gcoeff = 0.0
                    INT2d_bnds%Val(n2) = arseffx_bnd(kk, jj, 1) / (arseffx(k, j, i) + eps)
                    INT2d_bnds%ColInd(n2) = IndexU3block(1, jj, kk, blocks, bndblock_id, check=.TRUE.)
                    n2 = n2 + 1    
                  END DO
                END DO
              ELSE
                IF (ABS(z_bnd(k_bnds(1)) - z(k)) + ABS(y_bnd(j_bnds(1)) - y(j)) > eps) THEN

                  CALL boundint_weights_yz(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                                           x, y, z, y_bnd, z_bnd, y2_bnd, z2_bnd, volseff_bnd, arseffx_bnd, &
                                           j, k,  j_bnds, k_bnds, 1, point_dest)
                  DO l = 1, ipoint

                    gcoeff = 2.0 * weights(l) * arseffx(k, j, i) * points_volsarea(l, 2) / &
                             (volseff(k, j, i - 1) * points_volsarea(l, 2) + points_volsarea(l, 1) * arseffx(k, j, i) + eps)
                    gcoeff_sum = gcoeff_sum + gcoeff
                    GRAD3d_a%ColInd(n) = IndexC3block(points_loc(l, 2), points_loc(l, 3), points_loc(l, 4), &
                                                      blocks, points_loc(l, 1), check=.TRUE.)
                    GRAD3d_a%Val(n) = gcoeff
                    n = n + 1
                  END DO
                ELSE
                  volume = 0.5 * (volseff(k, j, i - 1) + volseff_bnd(k_bnds(1), j_bnds(1), 1))
                  gcoeff =  arseffx(k, j, i) / volume
                  gcoeff_sum = gcoeff
                  GRAD3d_a%ColInd(n) = IndexC3block(1, j_bnds(1), k_bnds(1), blocks, bndblock_id, check=.TRUE.)
                  GRAD3d_a%Val(n) = gcoeff
                  n = n + 1
                END IF
              END IF

              GRAD3d_a%Val(n_tmp) = -gcoeff_sum
              GRAD3d_a%RowPtr(irow + 1) = n
              INT2d_bnds%RowPtr(irow + 1) = n2
              irow = irow + 1
            END DO
          END DO
        ELSE
          DO j = 1, ny
            DO k = 1, nz
              GRAD3d_a%RowPtr(irow + 1) = n
              INT2d_bnds%RowPtr(irow + 1) = n2
              irow = irow + 1
            END DO
          END DO
        END IF

      ELSE

        DO j = 1, ny
          DO k = 1, nz
            area = blocks(iblock)%arseffx(k, j, i)
            volume = 0.5 * (blocks(iblock)%volseff(k, j, i - 1) + blocks(iblock)%volseff(k, j, i))
            gcoeff = area / volume
            GRAD3d_a%ColInd(n) = IndexC3block(i - 1, j, k, blocks, iblock, check=.TRUE.)
            GRAD3d_a%Val(n) = -gcoeff
            n = n + 1
            GRAD3d_a%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
            GRAD3d_a%Val(n) = gcoeff
            n = n + 1
            GRAD3d_a%RowPtr(irow + 1) = n
            INT2d_bnds%RowPtr(irow + 1) = n2
            irow = irow + 1
          END DO
        END DO
      END IF
    END DO


    ! y-component
    DO i = 1, nx
      DO j = 1, ny + 1
        
        IF (j == 1) THEN
          IF (ANY(blocks(iblock)%cface_s)) THEN
            DO k = 1, nz

              CALL find_indices_ik(blocks, blocks(iblock)%cface_s, i, i_bnds, k, k_bnds, iblock, bndblock_id)

              IF (bndblock_id .GT. 99999999) THEN
                GRAD3d_a%RowPtr(irow + 1) = n
                INT2d_bnds%RowPtr(irow + 1) = n2
                irow = irow + 1
                CYCLE
              END IF

              gcoeff_sum = 0.0

              vols_bnd => blocks(bndblock_id)%vols
              arsy_bnd => blocks(bndblock_id)%arsy
              volseff_bnd => blocks(bndblock_id)%volseff
              arseffy_bnd => blocks(bndblock_id)%arseffy

              x_bnd => blocks(bndblock_id)%x
              z_bnd => blocks(bndblock_id)%z
              x2_bnd => blocks(bndblock_id)%x2
              z2_bnd => blocks(bndblock_id)%z2
              nz_bnd = blocks(bndblock_id)%fld_shape(1)
              ny_bnd = blocks(bndblock_id)%fld_shape(2)
              nx_bnd = blocks(bndblock_id)%fld_shape(3)

              IF ((k_bnds(2) - k_bnds(1) > 0) .OR. (i_bnds(2) - i_bnds(1) > 0)) THEN
                DO ii = i_bnds(1), i_bnds(2)
                  DO kk = k_bnds(1), k_bnds(2)
                    gcoeff = 0.0
                    INT2d_bnds%Val(n2) = arseffy_bnd(kk, ny_bnd + 1, ii) / (arseffy(k, j, i) + eps)
                    INT2d_bnds%ColInd(n2) = IndexV3block(ii, ny_bnd + 1, kk, blocks, bndblock_id, check=.TRUE.)
                    n2 = n2 + 1
                  END DO
                END DO
              ELSE
                IF (ABS(z_bnd(k_bnds(1)) - z(k)) + ABS(x_bnd(i_bnds(1)) - x(i)) > eps) THEN
                  CALL boundint_weights_xz(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                                           x, y, z, x_bnd, z_bnd, x2_bnd, z2_bnd, volseff_bnd, arseffy_bnd, &
                                           i, k, i_bnds, k_bnds, -1, point_dest)
                  DO l = 1, ipoint
                    gcoeff = 2.0 * weights(l) * arseffy(k, j, i) * points_volsarea(l, 2) / &
                             (volseff(k, j, i) * points_volsarea(l, 2) + points_volsarea(l, 1) * arseffy(k, j, i) + eps)
                    gcoeff_sum = gcoeff_sum + gcoeff
                    GRAD3d_a%ColInd(n) = IndexC3block(points_loc(l, 2), points_loc(l, 3), points_loc(l, 4), &
                                                    blocks, points_loc(l, 1), check=.TRUE.)
                    GRAD3d_a%Val(n) = -gcoeff
                    n = n + 1
                  END DO
                ELSE
                  volume = 0.5 * (volseff(k, j, i) + volseff_bnd(k_bnds(1), ny_bnd, i_bnds(1)))
                  gcoeff =  arseffy(k, j, i) / volume
                  gcoeff_sum = gcoeff
                  GRAD3d_a%ColInd(n) = IndexC3block(i_bnds(1), SIZE(blocks(bndblock_id)%y2) - 1, k_bnds(1), &
                                                  blocks, bndblock_id, check=.TRUE.)
                  GRAD3d_a%Val(n) = -gcoeff
                  n = n + 1
                END IF
              END IF

              GRAD3d_a%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
              GRAD3d_a%Val(n) = gcoeff_sum
              n = n + 1

              GRAD3d_a%RowPtr(irow + 1) = n
              INT2d_bnds%RowPtr(irow + 1) = n2
              irow = irow + 1
            END DO

          ELSE
            DO k = 1, nz
              GRAD3d_a%RowPtr(irow + 1) = n
              INT2d_bnds%RowPtr(irow + 1) = n2
              irow = irow + 1
            END DO
          END IF
        ELSEIF (j == ny + 1) THEN
          IF (ANY(blocks(iblock)%cface_n)) THEN
            DO k = 1, nz

              CALL find_indices_ik(blocks, blocks(iblock)%cface_n, i, i_bnds, k, k_bnds, iblock, bndblock_id)

              IF (bndblock_id .GT. 99999999) THEN
                GRAD3d_a%RowPtr(irow + 1) = n
                INT2d_bnds%RowPtr(irow + 1) = n2
                irow = irow + 1
                CYCLE
              END IF

              vols_bnd => blocks(bndblock_id)%vols
              arsy_bnd => blocks(bndblock_id)%arsy
              volseff_bnd => blocks(bndblock_id)%volseff
              arseffy_bnd => blocks(bndblock_id)%arseffy

              x_bnd => blocks(bndblock_id)%x
              z_bnd => blocks(bndblock_id)%z
              x2_bnd => blocks(bndblock_id)%x2
              z2_bnd => blocks(bndblock_id)%z2
              nz_bnd = blocks(bndblock_id)%fld_shape(1)
              ny_bnd = blocks(bndblock_id)%fld_shape(2)
              nx_bnd = blocks(bndblock_id)%fld_shape(3)

              GRAD3d_a%ColInd(n) = IndexC3block(i, j - 1, k, blocks, iblock, check=.TRUE.)

              n_tmp = n
              n = n + 1

              gcoeff_sum = 0.0

              IF ((k_bnds(2) - k_bnds(1) > 0) .OR. (i_bnds(2) - i_bnds(1) > 0)) THEN

                DO ii = i_bnds(1), i_bnds(2)
                  DO kk = k_bnds(1), k_bnds(2)
                    gcoeff = 0.0
                    INT2d_bnds%Val(n2) = arseffy_bnd(kk, 1, ii) / (arseffy(k, j, i) + eps)
                    INT2d_bnds%ColInd(n2) = IndexV3block(ii, 1, kk, blocks, bndblock_id, check=.TRUE.)
                    n2 = n2 + 1
                  END DO
                END DO
              ELSE
                IF (ABS(z_bnd(k_bnds(1)) - z(k)) + ABS(x_bnd(i_bnds(1)) - x(i)) > eps) THEN

                  CALL boundint_weights_xz(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                                           x, y, z, x_bnd, z_bnd, x2_bnd, z2_bnd, volseff_bnd, arseffy_bnd, &
                                           i, k, i_bnds, k_bnds, 1, point_dest)
                  DO l = 1, ipoint

                    gcoeff = 2.0 * weights(l) * arseffy(k, j, i) * points_volsarea(l, 2) / &
                             (volseff(k, j - 1, i) * points_volsarea(l, 2) + points_volsarea(l, 1) * arseffy(k, j, i) + eps)
                    gcoeff_sum = gcoeff_sum + gcoeff
                    GRAD3d_a%ColInd(n) = IndexC3block(points_loc(l, 2), points_loc(l, 3), points_loc(l, 4), &
                                                      blocks, points_loc(l, 1), check=.TRUE.)
                    GRAD3d_a%Val(n) = gcoeff
                    n = n + 1

                  END DO
                ELSE
                  volume = 0.5 * (volseff(k, j - 1, i) + volseff_bnd(k_bnds(1), 1, i_bnds(1)))
                  gcoeff =  arseffy(k, j, i) / volume
                  gcoeff_sum = gcoeff
                  GRAD3d_a%ColInd(n) = IndexC3block(i_bnds(1), 1, k_bnds(1), blocks, bndblock_id, check=.TRUE.)
                  GRAD3d_a%Val(n) = gcoeff
                  n = n + 1
                END IF
                GRAD3d_a%Val(n_tmp) = -gcoeff_sum
  
              END IF
              GRAD3d_a%Val(n_tmp) = -gcoeff_sum
              GRAD3d_a%RowPtr(irow + 1) = n
              INT2d_bnds%RowPtr(irow + 1) = n2
              irow = irow + 1
            END DO
          ELSE
            DO k = 1, nz
              GRAD3d_a%RowPtr(irow + 1) = n
              INT2d_bnds%RowPtr(irow + 1) = n2
              irow = irow + 1
            END DO
          END IF
        ELSE
          DO k = 1, nz
            area = blocks(iblock)%arseffy(k, j, i)
            volume = 0.5 * (blocks(iblock)%volseff(k, j - 1, i) + blocks(iblock)%volseff(k, j, i))
            gcoeff = area / volume
            GRAD3d_a%ColInd(n) = IndexC3block(i, j - 1, k, blocks, iblock, check=.TRUE.)
            GRAD3d_a%Val(n) = -gcoeff
            n = n + 1
            GRAD3d_a%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
            GRAD3d_a%Val(n) = gcoeff
            n = n + 1
            GRAD3d_a%RowPtr(irow + 1) = n
            INT2d_bnds%RowPtr(irow + 1) = n2
            irow = irow + 1
          END DO
        END IF
      END DO
    END DO
    ! z-component 
    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz + 1

          IF (k == 1) THEN


            IF (ANY(blocks(iblock)%cface_b)) THEN
              
              CALL find_indices_ij(blocks, blocks(iblock)%cface_b, i, i_bnds, j, j_bnds, iblock, bndblock_id)

              IF (bndblock_id .GT. 99999999) THEN
                GRAD3d_a%RowPtr(irow + 1) = n
                INT2d_bnds%RowPtr(irow + 1) = n2
                irow = irow + 1
                CYCLE
              END IF

              gcoeff_sum = 0.0

              vols_bnd => blocks(bndblock_id)%vols
              arsz_bnd => blocks(bndblock_id)%arsz
              volseff_bnd => blocks(bndblock_id)%volseff
              arseffz_bnd => blocks(bndblock_id)%arseffz

              x_bnd => blocks(bndblock_id)%x
              y_bnd => blocks(bndblock_id)%y
              x2_bnd => blocks(bndblock_id)%x2
              y2_bnd => blocks(bndblock_id)%y2
              nz_bnd = blocks(bndblock_id)%fld_shape(1)
              ny_bnd = blocks(bndblock_id)%fld_shape(2)
              nx_bnd = blocks(bndblock_id)%fld_shape(3)

              IF ((j_bnds(2) - j_bnds(1) > 0) .OR. (i_bnds(2) - i_bnds(1) > 0)) THEN
                DO ii = i_bnds(1), i_bnds(2)
                  DO jj = j_bnds(1), j_bnds(2)
                    gcoeff = 0.0
                    INT2d_bnds%Val(n2) = arseffz_bnd(nz_bnd + 1, jj, ii) / (arseffz(k, j, i) + eps)
                    INT2d_bnds%ColInd(n2) = IndexW3block(ii, jj, nz_bnd + 1, blocks, bndblock_id, check=.TRUE.)
                    n2 = n2 + 1
                  END DO
                END DO
              ELSE

                IF (ABS(y_bnd(j_bnds(1)) - y(j)) + ABS(x_bnd(i_bnds(1)) - x(i)) > eps) THEN
                  CALL boundint_weights_xy(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                                           x, y, z, x_bnd, y_bnd, x2_bnd, y2_bnd, volseff_bnd, arseffz_bnd, &
                                           i, j, i_bnds, j_bnds, -1, point_dest)
                  DO l = 1, ipoint
                    gcoeff = 2.0 * weights(l) * arseffz(k, j, i) * points_volsarea(l, 2) / &
                             (volseff(k, j, i) * points_volsarea(l, 2) + points_volsarea(l, 1) * arseffz(k, j, i) + eps)
                    gcoeff_sum = gcoeff_sum + gcoeff
                    GRAD3d_a%ColInd(n) = IndexC3block(points_loc(l, 2), points_loc(l, 3), points_loc(l, 4), &
                                                    blocks, points_loc(l, 1), check=.TRUE.)
                    GRAD3d_a%Val(n) = -gcoeff
                    n = n + 1
                  END DO
                ELSE
                  volume = 0.5 * (volseff(k, j, i) + volseff_bnd(nz_bnd, j_bnds(1), i_bnds(1)))
                  gcoeff =  arseffz(k, j, i) / volume
                  gcoeff_sum = gcoeff
                  GRAD3d_a%ColInd(n) = IndexC3block(i_bnds(1), j_bnds(1), SIZE(blocks(bndblock_id)%z2) - 1, &
                                                  blocks, bndblock_id, check=.TRUE.)
                  GRAD3d_a%Val(n) = -gcoeff
                  n = n + 1
                END IF
              END IF

              GRAD3d_a%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
              GRAD3d_a%Val(n) = gcoeff_sum
              n = n + 1
            END IF

            GRAD3d_a%RowPtr(irow + 1) = n
            INT2d_bnds%RowPtr(irow + 1) = n2
            irow = irow + 1



          ELSE IF (k == nz + 1) THEN
            IF (ANY(blocks(iblock)%cface_t)) THEN
              CALL find_indices_ij(blocks, blocks(iblock)%cface_t, i, i_bnds, j, j_bnds, iblock, bndblock_id)

              IF (bndblock_id .GT. 99999999) THEN
                GRAD3d_a%RowPtr(irow + 1) = n
                INT2d_bnds%RowPtr(irow + 1) = n2
                irow = irow + 1
                CYCLE
              END IF

              gcoeff_sum = 0.0

              vols_bnd => blocks(bndblock_id)%vols
              arsz_bnd => blocks(bndblock_id)%arsz
              volseff_bnd => blocks(bndblock_id)%volseff
              arseffz_bnd => blocks(bndblock_id)%arseffz

              x_bnd => blocks(bndblock_id)%x
              y_bnd => blocks(bndblock_id)%y
              x2_bnd => blocks(bndblock_id)%x2
              y2_bnd => blocks(bndblock_id)%z
              nz_bnd = blocks(bndblock_id)%fld_shape(1)
              ny_bnd = blocks(bndblock_id)%fld_shape(2)
              nx_bnd = blocks(bndblock_id)%fld_shape(3)

              GRAD3d_a%ColInd(n) = IndexC3block(i, j, k - 1, blocks, iblock, check=.TRUE.)

              n_tmp = n
              n = n + 1

              gcoeff_sum = 0.0

              IF ((j_bnds(2) - j_bnds(1) > 0) .OR. (i_bnds(2) - i_bnds(1) > 0)) THEN
                DO ii = i_bnds(1), i_bnds(2)
                  DO jj = j_bnds(1), j_bnds(2)
                    gcoeff = 0.0
                    INT2d_bnds%Val(n2) = arseffz_bnd(1, jj, ii) / (arseffz(k, j, i) + eps)
                    INT2d_bnds%ColInd(n2) = IndexW3block(ii, jj, 1, blocks, bndblock_id, check=.TRUE.)
                    n2 = n2 + 1
                  END DO
                END DO
              ELSE
                IF (ABS(y_bnd(j_bnds(1)) - y(j)) + ABS(x_bnd(i_bnds(1)) - x(i)) > eps) THEN
                  CALL boundint_weights_xy(weights, points_loc, points_volsarea, ipoint, blocks, bndblock_id, &
                                           x, y, z, x_bnd, y_bnd, x2_bnd, y2_bnd, volseff_bnd, arseffz_bnd, &
                                           i, j, i_bnds, j_bnds, 1, point_dest)
                  DO l = 1, ipoint
                    gcoeff = 2.0 * weights(l) * arseffz(k, j, i) * points_volsarea(l, 2) / &
                             (volseff(k - 1, j, i) * points_volsarea(l, 2) + points_volsarea(l, 1) * arseffz(k, j, i) + eps)
                    gcoeff_sum = gcoeff_sum + gcoeff
                    GRAD3d_a%ColInd(n) = IndexC3block(points_loc(l, 2), points_loc(l, 3), points_loc(l, 4), &
                                                      blocks, points_loc(l, 1), check=.TRUE.)
                    GRAD3d_a%Val(n) = gcoeff
                    n = n + 1
                  END DO
                ELSE
                  volume = 0.5 * (volseff(k - 1, j, i) + volseff_bnd(1, j_bnds(1), i_bnds(1)))
                  gcoeff =  arseffz(k, j, i) / volume
                  gcoeff_sum = gcoeff
                  GRAD3d_a%ColInd(n) = IndexC3block(i_bnds(1), j_bnds(1), 1, blocks, bndblock_id, check=.TRUE.)
                  GRAD3d_a%Val(n) = gcoeff
                  n = n + 1
                END IF
              END IF

              GRAD3d_a%Val(n_tmp) = -gcoeff_sum
            END IF
            GRAD3d_a%RowPtr(irow + 1) = n
            INT2d_bnds%RowPtr(irow + 1) = n2
            irow = irow + 1

          ELSE
            area = blocks(iblock)%arseffz(k, j, i)
            volume = 0.5 * (blocks(iblock)%volseff(k - 1, j, i) + blocks(iblock)%volseff(k, j, i))
            gcoeff = area / volume
            GRAD3d_a%ColInd(n) = IndexC3block(i, j, k - 1, blocks, iblock, check=.TRUE.)
            GRAD3d_a%Val(n) = -gcoeff
            n = n + 1
            GRAD3d_a%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
            GRAD3d_a%Val(n) = gcoeff
            n = n + 1
            GRAD3d_a%RowPtr(irow + 1) = n
            INT2d_bnds%RowPtr(irow + 1) = n2
            irow = irow + 1
          END IF
        END DO
      END DO
    END DO
  END DO


!  CALL resize_sprowcol_dataflds(GRAD3d_a, n - 1)
!  CALL resize_sprowcol_dataflds(INT2d_bnds, n2 - 1)

  CALL SpMm_SpRowCol(Grad3d_b, INT2d_bnds, GRAD3d_a)
  CALL SpA_SpRowCol(Grad3d, Grad3d_a, Grad3d_b, allocate_mat=.FALSE.)

!  CALL SpDeallocate_SpRowCol(GRAD3d_a)
!  CALL SpDeallocate_SpRowCol(Grad3d_b)
!  CALL SpDeallocate_SpRowCol(INT2d_bnds)

  CALL SpNullify_SpRowCol(GRAD3d_a)
  CALL SpNullify_SpRowCol(GRAD3d_b)
  CALL SpNullify_SpRowCol(INT2d_bnds)

END SUBROUTINE make_grad3d



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



SUBROUTINE find_indices_j(blocks, bounding, j, j_bnds, myblock_id, bndblock_id)
  
  IMPLICIT NONE

  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  INTEGER, INTENT(in) :: j
  INTEGER, INTENT(in) :: myblock_id
  INTEGER, INTENT(inout) :: j_bnds(2)
  INTEGER, INTENT(inout) :: bndblock_id
  
  INTEGER :: iblock, nblocks
  INTEGER :: jj

  nblocks = SIZE(blocks)

  DO iblock = 1, nblocks
    IF (bounding(iblock) .EQV. .TRUE.) THEN
      IF (blocks(myblock_id)%y2(j + 1) .gt. blocks(iblock)%y2(1) .AND. &
          blocks(myblock_id)%y2(j) .lt. blocks(iblock)%y2(SIZE(blocks(iblock)%y2))) THEN
        j_bnds(1) = 1
        DO jj = 1, SIZE(blocks(iblock)%y2) - 1
          IF (blocks(iblock)%y2(jj) .gt. blocks(myblock_id)%y2(j)) EXIT         
          j_bnds(1) = jj
        END DO
        DO jj = j_bnds(1), SIZE(blocks(iblock)%y2) - 1
          j_bnds(2) = jj
          IF (blocks(iblock)%y2(jj + 1) .ge. blocks(myblock_id)%y2(j + 1)) EXIT
        END DO
        bndblock_id = iblock
        EXIT
      END IF
    END IF
  END DO

END SUBROUTINE find_indices_j



SUBROUTINE find_indices_i(blocks, bounding, i, i_bnds, myblock_id, bndblock_id)
  
  IMPLICIT NONE

  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  INTEGER, INTENT(in) :: i
  INTEGER, INTENT(in) :: myblock_id
  INTEGER, INTENT(inout) :: i_bnds(2)
  INTEGER, INTENT(inout) :: bndblock_id
  
  INTEGER :: iblock, nblocks
  INTEGER :: ii

  nblocks = SIZE(blocks)

  DO iblock = 1, nblocks
    IF (bounding(iblock) .EQV. .TRUE.) THEN 
      IF (blocks(myblock_id)%x2(i + 1) .gt. blocks(iblock)%x2(1) .AND. &
          blocks(myblock_id)%x2(i) .lt. blocks(iblock)%x2(SIZE(blocks(iblock)%x2))) THEN
        i_bnds(1) = 1
        DO ii = 1, SIZE(blocks(iblock)%x2) - 1
          IF (blocks(iblock)%x2(ii) .gt. blocks(myblock_id)%x2(i)) EXIT
          i_bnds(1) = ii
        END DO
        DO ii = i_bnds(1), SIZE(blocks(iblock)%x2) - 1
          i_bnds(2) = ii
          IF (blocks(iblock)%x2(ii + 1) .ge. blocks(myblock_id)%x2(i + 1)) EXIT
        END DO
        bndblock_id = iblock
        EXIT
      END IF
    END IF
  END DO

END SUBROUTINE find_indices_i


SUBROUTINE find_indices_k(blocks, bounding, i, i_bnds, myblock_id, bndblock_id)

  IMPLICIT NONE

  TYPE(block), POINTER, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  INTEGER, INTENT(in) :: i
  INTEGER, INTENT(in) :: myblock_id
  INTEGER, INTENT(inout) :: i_bnds(2)
  INTEGER, INTENT(inout) :: bndblock_id

  INTEGER :: iblock, nblocks
  INTEGER :: ii

  nblocks = SIZE(blocks)

  DO iblock = 1, nblocks
    IF (bounding(iblock) .EQV. .TRUE.) THEN
      IF (blocks(myblock_id)%z2(i + 1) .gt. blocks(iblock)%z2(1) .AND. &
          blocks(myblock_id)%z2(i) .lt. blocks(iblock)%z2(SIZE(blocks(iblock)%z2))) THEN
        i_bnds(1) = 1
        DO ii = 1, SIZE(blocks(iblock)%z2) - 1
          IF (blocks(iblock)%z2(ii) .gt. blocks(myblock_id)%z2(i)) EXIT
          i_bnds(1) = ii
        END DO
        DO ii = i_bnds(1), SIZE(blocks(iblock)%z2) - 1
          i_bnds(2) = ii
          IF (blocks(iblock)%z2(ii + 1) .ge. blocks(myblock_id)%z2(i + 1)) EXIT
        END DO
        bndblock_id = iblock
        EXIT
      END IF
    END IF
  END DO

END SUBROUTINE find_indices_k


PURE FUNCTION ubound_ncolids_blockgrad(blocks)

   TYPE(block), POINTER, INTENT(IN) :: blocks(:)
   INTEGER :: nblocks
   INTEGER :: ncolids
   INTEGER :: iblock
   INTEGER :: nz, ny, nx
   INTEGER :: ubound_ncolids_blockgrad

   IF (.NOT. ASSOCIATED(blocks)) THEN
     ubound_ncolids_blockgrad = 0
   ELSE

     nblocks = SIZE(blocks)

     ncolids = 0
   
     DO iblock = 1, nblocks
       nz = blocks(iblock)%fld_shape(1)
       ny = blocks(iblock)%fld_shape(2)
       nx = blocks(iblock)%fld_shape(3) 
       ncolids = ncolids + 2 * (nz + 2) * ny * nx + 2 * nz * (ny + 2) * nx + 2 * nz * ny * (nx + 2)
     
     END DO

     ubound_ncolids_blockgrad = ncolids
   END IF

END FUNCTION ubound_ncolids_blockgrad


PURE FUNCTION ubound_ncolids_int2d(blocks)

  TYPE(block), POINTER, INTENT(IN) :: blocks(:)
  INTEGER :: nblocks
  INTEGER :: ncolids
  INTEGER :: iblock
  INTEGER :: nz, ny, nx
  INTEGER :: ubound_ncolids_int2d


  IF (.NOT. ASSOCIATED(blocks)) THEN
    ubound_ncolids_int2d = 0
  ELSE

    nblocks = SIZE(blocks)

    ncolids = 0
  
    DO iblock = 1, nblocks
      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)
      ncolids = ncolids + 2 * nz * ny + 2 * nz * nx + 2 * ny * nx
 
    END DO

    ubound_ncolids_int2d = ncolids

  END IF    

END FUNCTION ubound_ncolids_int2d


END MODULE Gradient_block_Mod
