MODULE Block_Mod

USE Kind_Mod

USE List_Mod, ONLY: List

!A rectangular block
TYPE :: Block

  !Metadata from netcdf file
  INTEGER :: blockid    !Unique block identifier (blocks can be partitioned across different subdomains)
  INTEGER :: myrank
  INTEGER :: iglob_st, iglob_end !start and end i-index of block in global fine grid
  INTEGER :: jglob_st, jglob_end !start and end j-index of block in global fine grid
  INTEGER :: kglob_st, kglob_end !start and end k-index of block in global fine grid

  INTEGER :: fld_shape(3) !shape of 3D-field which may include ghost layers for communication (if the block is cut by a subdomain boundary).
  INTEGER :: ncells !Total number of computation cell-center nodes
  INTEGER :: nfaces !Total number of computation cell-face nodes
  INTEGER :: nghost
  INTEGER :: nbound
  INTEGER :: resx, resy, resz

  !Computation bounds
  INTEGER :: i_st, i_end !start and end i-index of vector x
  INTEGER :: j_st, j_end !start and end j-index of vector y
  INTEGER :: k_st, k_end !start and end k-index of vector z

  REAL(Realkind), POINTER :: x(:) => NULL()             !x-coordinate array of cell centres
  REAL(Realkind), POINTER :: y(:) => NULL()             !y-coordinate array of cell centres
  REAL(Realkind), POINTER :: z(:) => NULL()             !z-coordinate array of cell centres
  REAL(Realkind), POINTER :: x2(:) => NULL()            !x-coordinate array of x-normal cell-faces
  REAL(Realkind), POINTER :: y2(:) => NULL()            !y-coordinate array of y-normal cell-faces
  REAL(Realkind), POINTER :: z2(:) => NULL()            !z-coordinate array of z-normal cell-faces

  REAL(Realkind), POINTER :: volseff(:, :, :) => NULL() !effective volumes of cells
  REAL(Realkind), POINTER :: arseffx(:, :, :) => NULL() !effective areas of x-normal cell-faces
  REAL(Realkind), POINTER :: arseffy(:, :, :) => NULL() !effective areas of y-normal cell-faces
  REAL(Realkind), POINTER :: arseffz(:, :, :) => NULL() !effective areas of z-normal cell-faces
  REAL(Realkind), POINTER :: vols(:, :, :) => NULL()    !geometric volumes of cells
  REAL(Realkind), POINTER :: arsx(:, :, :) => NULL()    !geometric areas of x-normal cell-faces
  REAL(Realkind), POINTER :: arsy(:, :, :) => NULL()    !geometric areas of y-normal cell-faces
  REAL(Realkind), POINTER :: arsz(:, :, :) => NULL()    !geometric areas of z-normal cell-faces
  REAL(Realkind), POINTER :: scal(:, :, :) => NULL()    !volume scaling factor for residual scaling

  LOGICAL, POINTER :: cface_w(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one western face of this block 
  LOGICAL, POINTER :: cface_e(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one eastern face of this block 
  LOGICAL, POINTER :: cface_s(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one southern face of this block 
  LOGICAL, POINTER :: cface_n(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one northern face of this block 
  LOGICAL, POINTER :: cface_b(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one bottom face of this block 
  LOGICAL, POINTER :: cface_t(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one top face of this block 
  LOGICAL, POINTER :: cedge_w(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one western edge of this block 
  LOGICAL, POINTER :: cedge_e(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one eastern edge of this block 
  LOGICAL, POINTER :: cedge_s(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one southern edge of this block 
  LOGICAL, POINTER :: cedge_n(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one northern edge of this block 
  LOGICAL, POINTER :: cedge_b(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one bottom edge of this block 
  LOGICAL, POINTER :: cedge_t(:) => NULL()    !Boolean vector of size nblocks indicating which blocks share at least one top edge of this block 

  REAL(Realkind), POINTER :: volseff_tmp(:, :, :) => NULL() !effective volumes of cells
  REAL(Realkind), POINTER :: arseffx_tmp(:, :, :) => NULL() !effective areas of x-normal cell-faces
  REAL(Realkind), POINTER :: arseffy_tmp(:, :, :) => NULL() !effective areas of y-normal cell-faces
  REAL(Realkind), POINTER :: arseffz_tmp(:, :, :) => NULL() !effective areas of z-normal cell-faces
  REAL(Realkind), POINTER :: vols_tmp(:, :, :) => NULL()    !geometric volumes of cells
  REAL(Realkind), POINTER :: arsx_tmp(:, :, :) => NULL()    !geometric areas of x-normal cell-faces
  REAL(Realkind), POINTER :: arsy_tmp(:, :, :) => NULL()    !geometric areas of y-normal cell-faces
  REAL(Realkind), POINTER :: arsz_tmp(:, :, :) => NULL()    !geometric areas of z-normal cell-faces
  REAL(Realkind), POINTER :: scal_tmp(:, :, :) => NULL()    !volume scaling factor for residual scaling
  REAL(Realkind), POINTER :: x2_fine(:) => NULL()            !x-coordinate array of x-normal cell-faces
  REAL(Realkind), POINTER :: y2_fine(:) => NULL()            !y-coordinate array of y-normal cell-faces
  REAL(Realkind), POINTER :: z2_fine(:) => NULL()            !z-coordinate array of z-normal cell-faces


  INTEGER, POINTER :: global_neighbors(:) => NULL() !Block indices of all bounding blocks of the global domain


  CONTAINS

  PROCEDURE :: destroy

END TYPE Block


INTEGER, PARAMETER :: capacity = 2000
INTEGER, PARAMETER :: nfinemax = 20


TYPE :: Regular_grid
  INTEGER :: nx, ny, nz
  TYPE(gridding_cell), POINTER :: cells(:, :, :) => NULL()
  INTEGER, POINTER :: ix2(:) => NULL()
  INTEGER, POINTER :: jy2(:) => NULL()
  INTEGER, POINTER :: kz2(:) => NULL()
END TYPE Regular_grid

TYPE :: gridding_cell
  INTEGER :: contns(capacity) 
  INTEGER :: fgrid_inds(nfinemax, 3)
  INTEGER :: nfgrid_inds = 0
  INTEGER :: nconts = 0

END TYPE gridding_cell

CONTAINS

SUBROUTINE grid_blocks(finest_grid, block_ispartof, block_npartof, x2_glob, y2_glob, z2_glob, &
                       i_glob_st, i_glob_end, j_glob_st, j_glob_end, k_glob_st, k_glob_end)
  IMPLICIT NONE

  TYPE(Regular_grid), TARGET, INTENT(inout) :: finest_grid
  INTEGER, POINTER, INTENT(inout) :: block_ispartof(:, :, :)
  INTEGER, POINTER, INTENT(inout) :: block_npartof(:)


  REAL(Realkind), INTENT(in) :: x2_glob(:), &
                                y2_glob(:), &
                                z2_glob(:)
  INTEGER, INTENT(in) :: i_glob_st(:), i_glob_end(:), &
                         j_glob_st(:), j_glob_end(:), &
                         k_glob_st(:), k_glob_end(:)

  INTEGER :: nblocks, iblock

  INTEGER :: nx_glob, ny_glob, nz_glob
  REAL(Realkind) :: nx_block_avg, ny_block_avg, nz_block_avg

  TYPE(Regular_grid), ALLOCATABLE, TARGET :: grids_coarse(:)

  INTEGER :: ngrids_x, ngrids_y, ngrids_z
  INTEGER :: ngrids_coarse, igrid

  INTEGER :: i, j, k, ind, ii, jj, kk, icell
  REAL(Realkind) :: x2_first, x2_last, y2_first, y2_last, z2_first, z2_last

  TYPE(Regular_grid), POINTER :: grid_f, grid_c
  TYPE(gridding_cell), POINTER :: cell_f, cell_c

  nblocks = SIZE(i_glob_st)

  ALLOCATE(block_ispartof(nblocks, capacity, 3))
  ALLOCATE(block_npartof(nblocks))
  block_npartof(:) = 0

  nx_glob = SIZE(x2_glob) - 1
  ny_glob = SIZE(y2_glob) - 1
  nz_glob = SIZE(z2_glob) - 1

  nx_block_avg = REAL(SUM(i_glob_end - i_glob_st + 1)) / REAL(nblocks)
  ny_block_avg = REAL(SUM(j_glob_end - j_glob_st + 1)) / REAL(nblocks)
  nz_block_avg = REAL(SUM(k_glob_end - k_glob_st + 1)) / REAL(nblocks)


  finest_grid%nx = MAX(nx_glob / INT(nx_block_avg), 1)
  finest_grid%ny = MAX(ny_glob / INT(ny_block_avg), 1)
  finest_grid%nz = MAX(nz_glob / INT(nz_block_avg), 1)

  ngrids_x = MAX(INT(LOG(REAL(finest_grid%nx)) / LOG(2.0)), 1)
  ngrids_y = MAX(INT(LOG(REAL(finest_grid%ny)) / LOG(2.0)), 1)
  ngrids_z = MAX(INT(LOG(REAL(finest_grid%nz)) / LOG(2.0)), 1)

  ngrids_coarse = MAX(ngrids_x - 1, ngrids_y - 1 , ngrids_z - 1)
 
  ALLOCATE(grids_coarse(ngrids_coarse))

  ALLOCATE(finest_grid%cells(finest_grid%nz, finest_grid%ny, finest_grid%nx))
  ALLOCATE(finest_grid%ix2(finest_grid%nx + 1))
  ALLOCATE(finest_grid%jy2(finest_grid%ny + 1))
  ALLOCATE(finest_grid%kz2(finest_grid%nz + 1)) 

  DO i = 1, finest_grid%nx + 1
    finest_grid%ix2(i) = INT(REAL(nx_glob) / REAL(finest_grid%nx) * (i - 1)) + 1
  END DO 

  DO j = 1, finest_grid%ny + 1
    finest_grid%jy2(j) = INT(REAL(ny_glob) / REAL(finest_grid%ny) * (j - 1)) + 1
  END DO

  DO k = 1, finest_grid%nz + 1
    finest_grid%kz2(k) = INT(REAL(nz_glob) / REAL(finest_grid%nz) * (k - 1)) + 1
  END DO


  grid_f => finest_grid
  DO igrid = 1, ngrids_coarse

    grid_c => grids_coarse(igrid)

    grid_c%nx = MAX((grid_f%nx - 1) / 2 + 1, 1)
    grid_c%ny = MAX((grid_f%ny - 1) / 2 + 1, 1) 
    grid_c%nz = MAX((grid_f%nz - 1) / 2 + 1, 1)

    ALLOCATE(grid_c%cells(grid_c%nz, grid_c%ny, grid_c%nx))

    ALLOCATE(grid_c%ix2(grid_c%nx + 1))
    ALLOCATE(grid_c%jy2(grid_c%ny + 1))
    ALLOCATE(grid_c%kz2(grid_c%nz + 1))

    DO i = 1, grid_c%nx
      grid_c%ix2(i) = grid_f%ix2(2 * (i - 1) + 1)
    END DO
    grid_c%ix2(grid_c%nx + 1) = grid_f%ix2(grid_f%nx + 1)

    DO j = 1, grid_c%ny
      grid_c%jy2(j) = grid_f%jy2(2 * (j - 1) + 1)
    END DO
    grid_c%jy2(grid_c%ny + 1) = grid_f%jy2(grid_f%ny + 1)

    DO k = 1, grid_c%nz
      grid_c%kz2(k) = grid_f%kz2(2 * (k - 1) + 1)
    END DO
    grid_c%kz2(grid_c%nz + 1) = grid_f%kz2(grid_f%nz + 1)

    grid_f => grids_coarse(igrid)

  END DO

  IF (ngrids_coarse .GT. 0) THEN

    CALL connect_grids(finest_grid, grids_coarse(1))
    DO igrid = 1, ngrids_coarse - 1
      CALL connect_grids(grids_coarse(igrid), grids_coarse(igrid + 1))
    END DO

    grid_c => grids_coarse(ngrids_coarse)
  ELSE
    grid_c => finest_grid
  END IF

  DO i = 1, grid_c%nx 
    DO j = 1, grid_c%ny
      DO k = 1, grid_c%nz
        DO ind = 1, nblocks
          x2_first = x2_glob(i_glob_st(ind))
          x2_last = x2_glob(i_glob_end(ind) + 1)
          y2_first = y2_glob(j_glob_st(ind))
          y2_last = y2_glob(j_glob_end(ind) + 1)
          z2_first = z2_glob(k_glob_st(ind))
          z2_last = z2_glob(k_glob_end(ind) + 1)


          IF (.NOT. ANY((/(x2_first < x2_glob(grid_c%ix2(i)) .AND. x2_last < x2_glob(grid_c%ix2(i))), &
                        (x2_first > x2_glob(grid_c%ix2(i + 1)) .AND. x2_last > x2_glob(grid_c%ix2(i + 1))), &
                        (y2_first < y2_glob(grid_c%jy2(j)) .AND. y2_last < y2_glob(grid_c%jy2(j))), &
                        (y2_first > y2_glob(grid_c%jy2(j + 1)) .AND. y2_last > y2_glob(grid_c%jy2(j + 1))), &
                        (z2_first < z2_glob(grid_c%kz2(k)) .AND. z2_last < z2_glob(grid_c%kz2(k))), &
                        (z2_first > z2_glob(grid_c%kz2(k + 1)) .AND. z2_last > z2_glob(grid_c%kz2(k + 1))) &
                      /))) THEN
 
            cell_c => grid_c%cells(k, j, i)
            cell_c%nconts = cell_c%nconts + 1
            cell_c%contns(cell_c%nconts) = ind

            IF (ngrids_coarse .EQ. 0) THEN
              block_npartof(ind) = block_npartof(ind) + 1
              block_ispartof(ind, block_npartof(ind), 1) = k
              block_ispartof(ind, block_npartof(ind), 2) = j
              block_ispartof(ind, block_npartof(ind), 3) = i
            END IF
          END IF
        END DO 
      END DO
    END DO
  END DO

  

  IF (ngrids_coarse > 0) THEN
    DO igrid = ngrids_coarse - 1, 0, -1  
      IF (igrid == 0) THEN
        grid_f => finest_grid
        grid_c => grids_coarse(1)
      ELSE
        grid_f => grids_coarse(igrid)
        grid_c => grids_coarse(igrid + 1)
      END IF
      DO i = 1, grid_c%nx
        DO j = 1, grid_c%ny
          DO k = 1, grid_c%nz
            cell_c => grid_c%cells(k, j, i)
            IF (cell_c%nconts > 0) THEN
              DO ind = 1, cell_c%nconts

                x2_first = x2_glob(i_glob_st(cell_c%contns(ind)))
                x2_last = x2_glob(i_glob_end(cell_c%contns(ind)) + 1)
                y2_first = y2_glob(j_glob_st(cell_c%contns(ind)))
                y2_last = y2_glob(j_glob_end(cell_c%contns(ind)) + 1)
                z2_first = z2_glob(k_glob_st(cell_c%contns(ind)))
                z2_last = z2_glob(k_glob_end(cell_c%contns(ind)) + 1)

                DO icell = 1, cell_c%nfgrid_inds
                  ii = cell_c%fgrid_inds(icell, 3)
                  jj = cell_c%fgrid_inds(icell, 2)
                  kk = cell_c%fgrid_inds(icell, 1)

                  IF (.NOT. ANY((/(x2_first < x2_glob(grid_f%ix2(ii)) .AND. x2_last < x2_glob(grid_f%ix2(ii))), &
                                 (x2_first > x2_glob(grid_f%ix2(ii + 1)) .AND. x2_last > x2_glob(grid_f%ix2(ii + 1))), &
                                 (y2_first < y2_glob(grid_f%jy2(jj)) .AND. y2_last < y2_glob(grid_f%jy2(jj))), &
                                 (y2_first > y2_glob(grid_f%jy2(jj + 1)) .AND. y2_last > y2_glob(grid_f%jy2(jj + 1))), &
                                 (z2_first < z2_glob(grid_f%kz2(kk)) .AND. z2_last < z2_glob(grid_f%kz2(kk))), &
                                 (z2_first > z2_glob(grid_f%kz2(kk + 1)) .AND. z2_last > z2_glob(grid_f%kz2(kk + 1))) &
                               /))) THEN
                    
                    cell_f => grid_f%cells(kk, jj, ii)                  
                    cell_f%nconts = cell_f%nconts + 1
                    cell_f%contns(cell_f%nconts) = cell_c%contns(ind)

                    IF (igrid == 0) THEN
                      block_npartof(cell_c%contns(ind)) = block_npartof(cell_c%contns(ind)) + 1
                      block_ispartof(cell_c%contns(ind), block_npartof(cell_c%contns(ind)), 1) = kk
                      block_ispartof(cell_c%contns(ind), block_npartof(cell_c%contns(ind)), 2) = jj
                      block_ispartof(cell_c%contns(ind), block_npartof(cell_c%contns(ind)), 3) = ii
                    END IF
                  END IF

                END DO
              END DO 
           
            END IF
          END DO
        END DO
      END DO

    END DO
  END IF

  NULLIFY(cell_f)
  NULLIFY(cell_c)

  DO igrid = 1, ngrids_coarse
    grid_c => grids_coarse(igrid)
    DEALLOCATE(grid_c%cells)
    DEALLOCATE(grid_c%ix2)
    DEALLOCATE(grid_c%jy2)
    DEALLOCATE(grid_c%kz2)    
    grid_c%cells => NULL()
    grid_c%ix2 => NULL()
    grid_c%jy2 => NULL()
    grid_c%kz2 => NULL()
  END DO

  DEALLOCATE(grids_coarse)

!  DO i = 1, finest_grid%nx
!    DO j = 1, finest_grid%ny
!      DO k = 1, finest_grid%nz
!        WRITE(*,*) "CELL ", k, j, i, " CONTAINS BLOCKS" 
!        DO ind = 1, finest_grid%cells(k, j, i)%nconts

!          WRITE(*,*) finest_grid%cells(k, j, i)%contns(ind)

!        END DO
!      END DO
!    END DO
!  END DO
 
END SUBROUTINE grid_blocks


SUBROUTINE connect_grids(grid_f, grid_c)
  IMPLICIT NONE

  TYPE(Regular_grid), INTENT(in) :: grid_f
  TYPE(Regular_grid), INTENT(inout) :: grid_c

  INTEGER :: i, j, k, ii, jj, kk, ind_x, ind_y, ind_z
  REAL :: x_val, y_val, z_val
  INTEGER :: nfgrid_inds

  DO i = 1, grid_f%nx
    DO j = 1, grid_f%ny
      DO k = 1, grid_f%nz
        x_val = 0.5 * REAL(grid_f%ix2(i) + grid_f%ix2(i + 1))
        y_val = 0.5 * REAL(grid_f%jy2(j) + grid_f%jy2(j + 1))
        z_val = 0.5 * REAL(grid_f%kz2(k) + grid_f%kz2(k + 1))

        DO ii = 1, grid_c%nx
          IF (x_val .GE. grid_c%ix2(ii) .AND. x_val .LT. grid_c%ix2(ii + 1)) THEN
            ind_x = ii
            EXIT
          END IF
        END DO

        DO jj = 1, grid_c%ny
          IF (y_val .GE. grid_c%jy2(jj) .AND. y_val .LT. grid_c%jy2(jj + 1)) THEN
            ind_y = jj
          END IF
        END DO

        DO kk = 1, grid_c%nz
          IF (z_val .GE. grid_c%kz2(kk) .AND. z_val .LT. grid_c%kz2(kk + 1)) THEN
            ind_z = kk
          END IF
        END DO

        nfgrid_inds = grid_c%cells(ind_z, ind_y, ind_x)%nfgrid_inds
        nfgrid_inds = nfgrid_inds + 1
        grid_c%cells(ind_z, ind_y, ind_x)%fgrid_inds(nfgrid_inds, :) = (/k, j, i/)
        grid_c%cells(ind_z, ind_y, ind_x)%nfgrid_inds = nfgrid_inds

      END DO
    END DO
  END DO

END SUBROUTINE connect_grids

SUBROUTINE neighborblocks(neighbors, finest_grid, block_ispartof, block_npartof, thisblock, &
                          i_glob_st, i_glob_end, j_glob_st, j_glob_end, k_glob_st, k_glob_end, &
                          limit, cyclic_x, cyclic_y, cyclic_z) 
  IMPLICIT NONE

  TYPE(List), INTENT(inout) :: neighbors
  TYPE(Regular_Grid), INTENT(in) :: finest_grid
  INTEGER, POINTER, INTENT(in) :: block_ispartof(:,:,:)
  INTEGER, POINTER, INTENT(in) :: block_npartof(:)

  INTEGER, INTENT(in) :: thisblock
  INTEGER, INTENT(in) :: i_glob_st(:), i_glob_end(:), &
                         j_glob_st(:), j_glob_end(:), &
                         k_glob_st(:), k_glob_end(:)
  INTEGER, INTENT(in) :: limit
  LOGICAL, INTENT(in), OPTIONAL :: cyclic_x, cyclic_y, cyclic_z

  INTEGER :: nblocks_within
  INTEGER :: npartof
  INTEGER :: icell, iblock, ind
  INTEGER :: i, j, k
  LOGICAL :: cyclic_x_bnd, cyclic_y_bnd, cyclic_z_bnd
  TYPE(gridding_cell), POINTER :: cell


  IF (PRESENT(cyclic_x)) THEN
    cyclic_x_bnd = cyclic_x
  ELSE
    cyclic_x_bnd = .FALSE.
  END IF

  IF (PRESENT(cyclic_y)) THEN
    cyclic_y_bnd = cyclic_y
  ELSE
    cyclic_y_bnd = .FALSE.
  END IF

  IF (PRESENT(cyclic_z)) THEN
    cyclic_z_bnd = cyclic_z
  ELSE
    cyclic_z_bnd = .FALSE.
  END IF

  npartof = block_npartof(thisblock)

!  ineigh = 1

  DO icell = 1, npartof
    i = block_ispartof(thisblock, icell, 3)
    j = block_ispartof(thisblock, icell, 2)
    k = block_ispartof(thisblock, icell, 1)
    cell => finest_grid%cells(k, j, i)
    nblocks_within = cell%nconts

    DO ind = 1, nblocks_within
      iblock = cell%contns(ind)

      IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
        IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                        (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                        (j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                        (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                        (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                        (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock)), & 
                        (k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                        (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                        (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                        (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock)) &
                     /))) THEN 

          CALL neighbors%append(iblock)
        END IF
      END IF
    END DO 
  END DO

  ! for cyclic boundary conditions

  IF (cyclic_x_bnd) THEN
    IF (i_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1) THEN
          cell => finest_grid%cells(k, j, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end)) THEN
                IF (.NOT. ANY((/(j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                          (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                          (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                          (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock)), &
                          (k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                          (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                          (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                          (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock)) &
                       /))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx) THEN
          cell => finest_grid%cells(k, j, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                          (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                          (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                          (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock)), &
                          (k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                          (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                          (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                          (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock)) &
                       /))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF 
      END DO
    END IF
  END IF

  IF (cyclic_y_bnd) THEN
    IF (j_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (j .EQ. 1) THEN
          cell => finest_grid%cells(k, finest_grid%ny, i)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (iblock .LE. limit) THEN
              IF (.NOT. iblock .EQ. thisblock .AND. j_glob_end(iblock) .EQ. MAXVAL(j_glob_end)) THEN
                IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                          (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                          (i_glob_st(thisblock) .GT. i_glob_end(iblock) + 1) .AND. &
                          (i_glob_end(thisblock) + 1 .GT. i_glob_end(iblock)), &
                          (k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                          (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                          (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                          (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock)) &
                       /))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (j .EQ. finest_grid%ny) THEN
          cell => finest_grid%cells(k, 1, i)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (iblock .LE. limit) THEN
              IF (.NOT. iblock .EQ. thisblock .AND. j_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                          (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                          (i_glob_st(thisblock) .GT. i_glob_end(iblock) + 1) .AND. &
                          (i_glob_end(thisblock) + 1 .GT. i_glob_end(iblock)), &
                          (k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                          (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                          (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                          (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock)) &
                       /))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF
  END IF

  IF (cyclic_z_bnd) THEN
    IF (k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, j, i)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                          (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                          (i_glob_st(thisblock) .GT. i_glob_end(iblock) + 1) .AND. &
                          (i_glob_end(thisblock) + 1 .GT. i_glob_end(iblock)), &
                          (j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                          (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                          (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                          (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock)) &
                       /))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, j, i)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (k_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                          (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                          (i_glob_st(thisblock) .GT. i_glob_end(iblock) + 1) .AND. &
                          (i_glob_end(thisblock) + 1 .GT. i_glob_end(iblock)), &
                          (j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                          (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                          (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                          (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock)) &
                       /))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF
  END IF

  IF (cyclic_x_bnd .AND. cyclic_y_bnd) THEN
    IF (i_glob_st(thisblock) .EQ. 1 .AND. j_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1 .AND. j .EQ. 1) THEN
          cell => finest_grid%cells(k, finest_grid%ny, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end) .AND. j_glob_end(iblock) .EQ. MAXVAL(j_glob_end)) THEN
                IF (.NOT. ANY((/(k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                          (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                          (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                          (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (i_glob_st(thisblock) .EQ. 1 .AND. j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1 .AND. j .EQ. finest_grid%ny) THEN
          cell => finest_grid%cells(k, 1, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end) .AND. j_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                          (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                          (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                          (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end) .AND. j_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx .AND. j .EQ. 1) THEN
          cell => finest_grid%cells(k, finest_grid%ny, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1 .AND. j_glob_end(iblock) .EQ. MAXVAL(j_glob_end)) THEN
                IF (.NOT. ANY((/(k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                          (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                          (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                          (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end) .AND. j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx .AND. j .EQ. finest_grid%ny) THEN
          cell => finest_grid%cells(k, 1, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1 .AND. j_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(k_glob_st(thisblock) .LT. k_glob_st(iblock)) .AND. &
                          (k_glob_end(thisblock) + 1 .LT. k_glob_st(iblock)), &
                          (k_glob_st(thisblock) .GT. k_glob_end(iblock) + 1) .AND. &
                          (k_glob_end(thisblock) + 1 .GT. k_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF
  END IF

  IF (cyclic_x_bnd .AND. cyclic_z_bnd) THEN
    IF (i_glob_st(thisblock) .EQ. 1 .AND. k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1 .AND. k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, j, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end) .AND. k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                IF (.NOT. ANY((/(j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                          (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                          (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                          (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end) .AND. k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx .AND. k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, j, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1 .AND. k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                IF (.NOT. ANY((/(j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                          (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                          (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                          (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (i_glob_st(thisblock) .EQ. 1 .AND. k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1 .AND. k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, j, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end) .AND. k_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                          (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                          (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                          (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end) .AND. k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx .AND. k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, j, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1 .AND. k_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(j_glob_st(thisblock) .LT. j_glob_st(iblock)) .AND. &
                          (j_glob_end(thisblock) + 1 .LT. j_glob_st(iblock)), &
                          (j_glob_st(thisblock) .GT. j_glob_end(iblock) + 1) .AND. &
                          (j_glob_end(thisblock) + 1 .GT. j_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF
  END IF

  IF (cyclic_y_bnd .AND. cyclic_z_bnd) THEN
    IF (j_glob_st(thisblock) .EQ. 1 .AND. k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (j .EQ. 1 .AND. k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, finest_grid%ny, i)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (j_glob_end(iblock) .EQ. MAXVAL(j_glob_end) .AND. k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                          (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                          (i_glob_st(thisblock) .GT. i_glob_end(iblock) + 1) .AND. &
                          (i_glob_end(thisblock) + 1 .GT. i_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end) .AND. k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (j .EQ. finest_grid%ny .AND. k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, 1, i)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (j_glob_st(iblock) .EQ. 1 .AND. k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                          (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                          (i_glob_st(thisblock) .GT. i_glob_end(iblock) + 1) .AND. &
                          (i_glob_end(thisblock) + 1 .GT. i_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (j_glob_st(thisblock) .EQ. 1 .AND. k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (j .EQ. 1 .AND. k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, finest_grid%ny, i)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (j_glob_end(iblock) .EQ. MAXVAL(j_glob_end) .AND. k_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                          (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                          (i_glob_st(thisblock) .GT. i_glob_end(iblock) + 1) .AND. &
                          (i_glob_end(thisblock) + 1 .GT. i_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end) .AND. k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (j .EQ. finest_grid%ny .AND. k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, 1, i)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (j_glob_st(iblock) .EQ. 1 .AND. k_glob_st(iblock) .EQ. 1) THEN
                IF (.NOT. ANY((/(i_glob_st(thisblock) .LT. i_glob_st(iblock)) .AND. &
                          (i_glob_end(thisblock) + 1 .LT. i_glob_st(iblock)), &
                          (i_glob_st(thisblock) .GT. i_glob_end(iblock) + 1) .AND. &
                          (i_glob_end(thisblock) + 1 .GT. i_glob_end(iblock))/))) THEN
                  CALL neighbors%append(iblock)
                END IF
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF
  END IF

  IF (cyclic_x_bnd .AND. cyclic_y_bnd .AND. cyclic_z_bnd) THEN
    IF (i_glob_st(thisblock) .EQ. 1 .AND. j_glob_st(thisblock) .EQ. 1 .AND. k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1 .AND. j .EQ. 1 .AND. k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, finest_grid%ny, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end) .AND. &
                  j_glob_end(iblock) .EQ. MAXVAL(j_glob_end) .AND. &
                  k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                CALL neighbors%append(iblock)
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end) .AND. &
             j_glob_st(thisblock) .EQ. 1 .AND. k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx .AND. j .EQ. 1 .AND. k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, finest_grid%ny, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1 .AND. &
                  j_glob_end(iblock) .EQ. MAXVAL(j_glob_end) .AND. &
                  k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                CALL neighbors%append(iblock)
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end) .AND. &
             j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end) .AND. k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx .AND. j .EQ. finest_grid%ny .AND. k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, 1, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1 .AND. &
                  j_glob_st(iblock) .EQ. 1 .AND. &
                  k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                CALL neighbors%append(iblock)
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (i_glob_st(thisblock) .EQ. 1 .AND. &
             j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end) .AND. k_glob_st(thisblock) .EQ. 1) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1 .AND. j .EQ. finest_grid%ny .AND. k .EQ. 1) THEN
          cell => finest_grid%cells(finest_grid%nz, 1, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end) .AND. &
                  j_glob_st(iblock) .EQ. 1 .AND. &
                  k_glob_end(iblock) .EQ. MAXVAL(k_glob_end)) THEN
                CALL neighbors%append(iblock)
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (i_glob_st(thisblock) .EQ. 1 .AND. &
             j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end) .AND. k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof      
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1 .AND. j .EQ. finest_grid%ny .AND. k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, 1, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end) .AND. &
                  j_glob_st(iblock) .EQ. 1 .AND. &
                  k_glob_st(iblock) .EQ. 1) THEN
                CALL neighbors%append(iblock)
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (i_glob_st(thisblock) .EQ. 1 .AND. &
             j_glob_st(thisblock) .EQ. 1 .AND. k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. 1 .AND. j .EQ. 1 .AND. k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, finest_grid%ny, finest_grid%nx)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_end(iblock) .EQ. MAXVAL(i_glob_end) .AND. &
                  j_glob_end(iblock) .EQ. MAXVAL(j_glob_end) .AND. &
                  k_glob_st(iblock) .EQ. 1) THEN
                CALL neighbors%append(iblock)
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end) .AND. &
             j_glob_st(thisblock) .EQ. 1 .AND. k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx .AND. j .EQ. 1 .AND. k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, finest_grid%ny, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1 .AND. &
                  j_glob_end(iblock) .EQ. MAXVAL(j_glob_end) .AND. &
                  k_glob_st(iblock) .EQ. 1) THEN
                CALL neighbors%append(iblock)
              END IF
            END IF
          END DO
        END IF
      END DO

    ELSE IF (i_glob_end(thisblock) .EQ. MAXVAL(i_glob_end) .AND. &
             j_glob_end(thisblock) .EQ. MAXVAL(j_glob_end) .AND. k_glob_end(thisblock) .EQ. MAXVAL(k_glob_end)) THEN
      DO icell = 1, npartof
        i = block_ispartof(thisblock, icell, 3)
        j = block_ispartof(thisblock, icell, 2)
        k = block_ispartof(thisblock, icell, 1)
        IF (i .EQ. finest_grid%nx .AND. j .EQ. finest_grid%ny .AND. k .EQ. finest_grid%nz) THEN
          cell => finest_grid%cells(1, 1, 1)
          nblocks_within = cell%nconts
          DO ind = 1, nblocks_within
            iblock = cell%contns(ind)
            IF (.NOT. iblock .EQ. thisblock .AND. iblock .LE. limit) THEN
              IF (i_glob_st(iblock) .EQ. 1 .AND. &
                  j_glob_st(iblock) .EQ. 1 .AND. &
                  k_glob_st(iblock) .EQ. 1) THEN
                CALL neighbors%append(iblock)
              END IF
            END IF
          END DO
        END IF
      END DO
    END IF
  END IF

END SUBROUTINE neighborblocks


FUNCTION isin(array, val)

  IMPLICIT NONE

  INTEGER, INTENT(in) :: array(:)
  INTEGER, INTENT(in) :: val
  LOGICAL :: isin
  INTEGER :: i

  isin = .FALSE.
  DO i = 1, SIZE(array)
    IF (array(i) == val) THEN
      isin = .TRUE.
      EXIT
    END IF
  END DO

END FUNCTION isin


SUBROUTINE destroy(self) 
  CLASS(Block) :: self
   
  IF(ASSOCIATED(self%x)) THEN
    DEALLOCATE(self%x)
    self%x => NULL()
  END IF

  IF(ASSOCIATED(self%y)) THEN
    DEALLOCATE(self%y)
    self%y => NULL()
  END IF

  IF(ASSOCIATED(self%z)) THEN
    DEALLOCATE(self%z)
    self%z => NULL()
  END IF

  IF(ASSOCIATED(self%x2)) THEN
    DEALLOCATE(self%x2)
    self%x2 => NULL()
  END IF

  IF(ASSOCIATED(self%y2)) THEN
    DEALLOCATE(self%y2)
    self%y2 => NULL()
  END IF

  IF(ASSOCIATED(self%z2)) THEN
    DEALLOCATE(self%z2)
    self%z2 => NULL()
  END IF

  IF(ASSOCIATED(self%x)) THEN
    DEALLOCATE(self%x)
    self%x => NULL()
  END IF

  IF(ASSOCIATED(self%y)) THEN
    DEALLOCATE(self%y)
    self%y => NULL()
  END IF

  IF(ASSOCIATED(self%volseff)) THEN
    DEALLOCATE(self%volseff)
    self%volseff => NULL()
  END IF

  IF(ASSOCIATED(self%arseffx)) THEN
    DEALLOCATE(self%arseffx)
    self%arseffx => NULL()
  END IF

  IF(ASSOCIATED(self%arseffy)) THEN
    DEALLOCATE(self%arseffy)
    self%arseffy => NULL()
  END IF  

  IF(ASSOCIATED(self%arseffz)) THEN
    DEALLOCATE(self%arseffz)
    self%arseffz => NULL()
  END IF

  IF(ASSOCIATED(self%scal)) THEN
    DEALLOCATE(self%scal)
    self%scal => NULL()
  END IF


  IF(ASSOCIATED(self%vols)) THEN
    DEALLOCATE(self%vols)
    self%vols => NULL()
  END IF

  IF(ASSOCIATED(self%arsx)) THEN
    DEALLOCATE(self%arsx)
    self%arsx => NULL()
  END IF

  IF(ASSOCIATED(self%arsy)) THEN
    DEALLOCATE(self%arsy)
    self%arsy => NULL()
  END IF

  IF(ASSOCIATED(self%arsz)) THEN
    DEALLOCATE(self%arsz)
    self%arsz => NULL()
  END IF

  IF(ASSOCIATED(self%cface_w)) THEN
    DEALLOCATE(self%cface_w)
    self%cface_w => NULL()
  END IF

  IF(ASSOCIATED(self%cface_e)) THEN
    DEALLOCATE(self%cface_e)
    self%cface_e => NULL()
  END IF

  IF(ASSOCIATED(self%cface_s)) THEN
    DEALLOCATE(self%cface_s)
    self%cface_s => NULL()
  END IF

  IF(ASSOCIATED(self%cface_n)) THEN
    DEALLOCATE(self%cface_n)
    self%cface_n => NULL()
  END IF

  IF(ASSOCIATED(self%cface_b)) THEN
    DEALLOCATE(self%cface_b)
    self%cface_b => NULL()
  END IF

  IF(ASSOCIATED(self%cface_t)) THEN
    DEALLOCATE(self%cface_t)
    self%cface_t => NULL()
  END IF

  IF(ASSOCIATED(self%cedge_w)) THEN
    DEALLOCATE(self%cedge_w)
    self%cedge_w => NULL()
  END IF

  IF(ASSOCIATED(self%cedge_e)) THEN
    DEALLOCATE(self%cedge_e)
    self%cedge_e => NULL()
  END IF

  IF(ASSOCIATED(self%cedge_s)) THEN
    DEALLOCATE(self%cedge_s)
    self%cedge_s => NULL()
  END IF

  IF(ASSOCIATED(self%cedge_n)) THEN
    DEALLOCATE(self%cedge_n)
    self%cedge_n => NULL()
  END IF

  IF(ASSOCIATED(self%cedge_b)) THEN
    DEALLOCATE(self%cedge_b)
    self%cedge_b => NULL()
  END IF

  IF(ASSOCIATED(self%cedge_t)) THEN
    DEALLOCATE(self%cedge_t)
    self%cedge_t => NULL()
  END IF

  IF(ASSOCIATED(self%global_neighbors)) THEN
    DEALLOCATE(self%global_neighbors)
    self%global_neighbors => NULL()
  END IF

  IF(ASSOCIATED(self%volseff_tmp)) THEN
    DEALLOCATE(self%volseff_tmp)
    self%volseff_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%arseffx_tmp)) THEN
    DEALLOCATE(self%arseffx_tmp)
    self%arseffx_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%arseffy_tmp)) THEN
    DEALLOCATE(self%arseffy_tmp)
    self%arseffy_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%arseffz_tmp)) THEN
    DEALLOCATE(self%arseffz_tmp)
    self%arseffz_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%vols_tmp)) THEN
    DEALLOCATE(self%vols_tmp)
    self%vols_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%arsx_tmp)) THEN
    DEALLOCATE(self%arsx_tmp)
    self%arsx_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%arsy_tmp)) THEN
    DEALLOCATE(self%arsy_tmp)
    self%arsy_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%arsz_tmp)) THEN
    DEALLOCATE(self%arsz_tmp)
    self%arsz_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%scal_tmp)) THEN
    DEALLOCATE(self%scal_tmp)
    self%scal_tmp => NULL()
  END IF

  IF(ASSOCIATED(self%z2_fine)) THEN
    DEALLOCATE(self%z2_fine)
    self%z2_fine => NULL()
  END IF
    
  IF(ASSOCIATED(self%x2_fine)) THEN
    DEALLOCATE(self%x2_fine)
    self%x2_fine => NULL()
  END IF
    
  IF(ASSOCIATED(self%y2_fine)) THEN
    DEALLOCATE(self%y2_fine)
    self%y2_fine => NULL()
  END IF


END SUBROUTINE destroy



END MODULE Block_Mod
