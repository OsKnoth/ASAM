MODULE DomDecomp_Mod

USE Kind_Mod
USE netcdf
USE Index_Mod, ONLY: coord_ilbnd, coord_iubnd, findloc_real, &
                     IndexC3, IndexU3, IndexV3, IndexW3
USE Block_Mod, ONLY: Block, grid_blocks, neighborblocks, Regular_Grid
USE Subdom_Mod, ONLY: SubDomain
USE Boundary_Mod, ONLY: Bound, cyclic_x, cyclic_y, cyclic_z
USE MatSpRowCol_Mod, ONLY: SpRowCol
USE Part_Mod, ONLY: set_ranks_automatic
USE List_Mod, ONLY: List, Element

IMPLICIT NONE

INTEGER :: nblocks_total !Total number of blocks

!Edge coordinates of global fine grid
REAL(Realkind), ALLOCATABLE :: x2glob(:), y2glob(:), z2glob(:)

!Block boundary indices of global fine grid
INTEGER, ALLOCATABLE :: ix_first(:), ix_last(:), &
                        jy_first(:), jy_last(:), &
                        kz_first(:), kz_last(:)

INTEGER, ALLOCATABLE :: resx(:), resy(:), resz(:), btype(:)
INTEGER, ALLOCATABLE :: blockranks(:)

INTEGER, POINTER :: block_ispartof(:, :, :)
INTEGER, POINTER :: block_npartof(:)
TYPE(Regular_Grid) :: finest_grid


CONTAINS

SUBROUTINE read_blocks_from_netcdf(file_path, subdom)
  !Each process opens the netcdf-file with the block-grid descrioption and reads in
  !the block metadata.
  !If a given block overlaps with the subdomain, the overlapping section including
  ! a ghost layer is read in.
  !The blocks are "connected" by searching for common boundaries
  !The ghost cells of each block are defined, by setting the lower and upper bounds 
  !(i_st, i_end, j_st, j_end, k_st, k_end) of the computation fields.

  IMPLICIT NONE
  
  CHARACTER(LEN=128), INTENT(in) :: file_path ! path to file including file name
  TYPE(SubDomain), INTENT(inout) :: subdom    ! subdomain
  INTEGER :: nblocks_subdom                   ! Number of blocks part of the subdomain (including potential ghost blocks)

  INTEGER :: ncid, nc_err                                 !Some integers necessary to read in the netcdf file
  INTEGER :: var_ixf_id, var_jyf_id, var_kzf_id,       & 
             var_ixl_id, var_jyl_id, var_kzl_id,       &
             var_x2_id, var_y2_id, var_z2_id,          &
             var_farsx_id, var_farsy_id, var_farsz_id, &
             var_fvol_id, var_resx_id, var_resy_id,    &
             var_resz_id, var_btype_id, var_part_id

  INTEGER :: dimid_blocks                                
  INTEGER, DIMENSION(1) :: dimid_x, dimid_y, dimid_z
  INTEGER, DIMENSION(3) :: dimid_3d

  LOGICAL, ALLOCATABLE :: is_in_tmp(:), is_ghost(:)

  INTEGER :: n, m
  INTEGER :: i, j, k
  INTEGER :: nx, ny, nz
  INTEGER :: i1d
  INTEGER :: dimx_tmp_size, dimy_tmp_size, dimz_tmp_size
  INTEGER :: dimx_glob_size, dimy_glob_size, dimz_glob_size

  TYPE(List) :: neighbors
  TYPE(Element), POINTER :: elm => NulL()

  INTEGER, DIMENSION(1, 7) :: start_read
  INTEGER, DIMENSION(1) :: nread
  INTEGER :: i_st, i_end, j_st, j_end, k_st, k_end

  CHARACTER(LEN=128) :: var_str

  CHARACTER(LEN=8) :: numb2str
  REAL(Realkind), ALLOCATABLE :: fvol_tmp(:), fax_tmp(:), &
                                 fay_tmp(:), faz_tmp(:)


  nc_err = nf90_open(file_path, nf90_nowrite, ncid)
  nc_err = nf90_inq_dimid(ncid, 'nblocks', dimid_blocks)
  nc_err = nf90_inquire_dimension(ncid, dimid_blocks, len=nblocks_total)

!  CALL set_ranks_automatic(blockranks, nblocks_total)

!  nc_err = nf90_inq_varid(ncid, "partition", var_part_id)
! ALLOCATE(blockranks(nblocks_total))
!  nc_err = nf90_get_var(ncid, var_part_id, blockranks)
! blockranks(:) = 0

  ALLOCATE(is_in_tmp(nblocks_total))
  is_in_tmp(:) = .FALSE.

  nc_err = nf90_inq_varid(ncid, "x2glob", var_x2_id)
  nc_err = nf90_inq_varid(ncid, "y2glob", var_y2_id)
  nc_err = nf90_inq_varid(ncid, "z2glob", var_z2_id)
  nc_err = nf90_inquire_variable(ncid, var_x2_id, dimids=dimid_x)
  nc_err = nf90_inquire_variable(ncid, var_y2_id, dimids=dimid_y)
  nc_err = nf90_inquire_variable(ncid, var_z2_id, dimids=dimid_z)
  nc_err = nf90_inquire_dimension(ncid, dimid_x(1), len=dimx_glob_size)
  nc_err = nf90_inquire_dimension(ncid, dimid_y(1), len=dimy_glob_size)
  nc_err = nf90_inquire_dimension(ncid, dimid_z(1), len=dimz_glob_size)

  ALLOCATE(x2glob(dimx_glob_size))
  ALLOCATE(y2glob(dimy_glob_size))
  ALLOCATE(z2glob(dimz_glob_size))

  nc_err = nf90_get_var(ncid, var_x2_id, x2glob)
  nc_err = nf90_get_var(ncid, var_y2_id, y2glob)
  nc_err = nf90_get_var(ncid, var_z2_id, z2glob)  

  ALLOCATE(ix_first(nblocks_total))
  ALLOCATE(jy_first(nblocks_total))
  ALLOCATE(kz_first(nblocks_total))
  ALLOCATE(ix_last(nblocks_total))
  ALLOCATE(jy_last(nblocks_total))
  ALLOCATE(kz_last(nblocks_total))

  nc_err = nf90_inq_varid(ncid, 'ix_first', var_ixf_id)
  nc_err = nf90_inq_varid(ncid, 'jy_first', var_jyf_id)
  nc_err = nf90_inq_varid(ncid, 'kz_first', var_kzf_id)
  nc_err = nf90_inq_varid(ncid, 'ix_last', var_ixl_id)
  nc_err = nf90_inq_varid(ncid, 'jy_last', var_jyl_id)
  nc_err = nf90_inq_varid(ncid, 'kz_last', var_kzl_id)
  nc_err = nf90_get_var(ncid, var_ixf_id, ix_first)
  nc_err = nf90_get_var(ncid, var_jyf_id, jy_first)
  nc_err = nf90_get_var(ncid, var_kzf_id, kz_first)
  nc_err = nf90_get_var(ncid, var_ixl_id, ix_last)
  nc_err = nf90_get_var(ncid, var_jyl_id, jy_last)
  nc_err = nf90_get_var(ncid, var_kzl_id, kz_last)

  DO n = 1, nblocks_total
    ix_first(n) = ix_first(n) + 1
    jy_first(n) = jy_first(n) + 1
    kz_first(n) = kz_first(n) + 1
  END DO

  ALLOCATE(resx(nblocks_total))
  ALLOCATE(resy(nblocks_total))
  ALLOCATE(resz(nblocks_total))
  ALLOCATE(btype(nblocks_total))
  btype(:) = 1

  nc_err = nf90_inq_varid(ncid, 'resx_level', var_resx_id)
  nc_err = nf90_inq_varid(ncid, 'resy_level', var_resy_id)
  nc_err = nf90_inq_varid(ncid, 'resz_level', var_resz_id)
!  nc_err = nf90_inq_varid(ncid, 'btype', var_btype_id)

  nc_err = nf90_get_var(ncid, var_resx_id, resx)
  nc_err = nf90_get_var(ncid, var_resy_id, resy)
  nc_err = nf90_get_var(ncid, var_resz_id, resz)
!  nc_err = nf90_get_var(ncid, var_btype_id, btype)

  !Determine all ghost blocks of the subdomain block selection

  CALL  grid_blocks(finest_grid, block_ispartof, block_npartof, x2glob, y2glob, z2glob, &
                    ix_first, ix_last, jy_first, jy_last, kz_first, kz_last)


  m = 0
  DO n = 1, SIZE(blockranks)
    IF (blockranks(n) .EQ. subdom%myrank) THEN
      m = m + 1
    END IF
  END DO

  ALLOCATE(subdom%blockids_comp(m))
  
  m = 1 
  DO n = 1, SIZE(blockranks)
    IF (blockranks(n) .EQ. subdom%myrank) THEN
      subdom%blockids_comp(m) = n
      m = m + 1
    END IF
  END DO
 
  ALLOCATE(is_ghost(nblocks_total))
  is_ghost = .FALSE.
  DO n = 1, SIZE(subdom%blockids_comp)
    is_ghost(subdom%blockids_comp(n)) = .TRUE. 
!    neighbors => neighborblocks(finest_grid, block_ispartof, block_npartof, subdom%blockids_comp(n), &
!                     ix_first, ix_last, jy_first, jy_last, kz_first, kz_last, cyclic_x, cyclic_y, cyclic_z)
    CALL neighborblocks(neighbors, finest_grid, block_ispartof, block_npartof, subdom%blockids_comp(n), &
                        ix_first, ix_last, jy_first, jy_last, kz_first, kz_last, HUGE(n), cyclic_x, cyclic_y, cyclic_z)

    IF (ASSOCIATED(neighbors%first)) THEN
      elm => neighbors%first
      is_ghost(elm%ivalue) = .TRUE.
      DO WHILE(ASSOCIATED(elm%next))
        elm => elm%next
        is_ghost(elm%ivalue) = .TRUE.
      END DO
    END IF
    elm => NULL()
    CALL neighbors%destroy()
  END DO

  nblocks_subdom = COUNT(is_ghost)
 
  ALLOCATE(subdom%blockids_compghst(nblocks_subdom))
  m = 1
  DO n = 1, SIZE(is_ghost)
    IF (is_ghost(n) .EQV. .TRUE.) THEN
      subdom%blockids_compghst(m) = n
      m = m + 1
    END IF
  END DO

  ! Blocks part of the subdomain are not ghost blocks
  DO n = 1, SIZE(is_ghost)
    IF (isin(subdom%blockids_comp, n)) is_ghost(n) = .FALSE.
  END DO

  ! All blocks that are not ghost blocks are assigned True
  ALLOCATE(subdom%blockiscomp(nblocks_subdom))
  subdom%blockiscomp = .FALSE.
  DO n = 1, nblocks_subdom
    IF (is_ghost(subdom%blockids_compghst(n)) .EQV. .FALSE.) THEN
      subdom%blockiscomp(n) = .TRUE.
    END IF
  END DO 

  subdom%nblocks = nblocks_subdom
  ALLOCATE(subdom%blocks(nblocks_subdom))
  subdom%ncells_tmp = 0
  subdom%nfaces_tmp = 0
  subdom%ncells = 0
  subdom%nfaces = 0

  start_read(:, :) = 1
 
  var_str = 'x2b'
  nc_err = nf90_inq_varid(ncid, var_str, var_x2_id)
  var_str = 'y2b'
  nc_err = nf90_inq_varid(ncid, var_str, var_y2_id)
  var_str = 'z2b'
  nc_err = nf90_inq_varid(ncid, var_str, var_z2_id)

 ! Get the scaling fields
  var_str = 'fvolb'
  nc_err = nf90_inq_varid(ncid, var_str, var_fvol_id)
  var_str = 'faxb'
  nc_err = nf90_inq_varid(ncid, var_str, var_farsx_id)
  var_str = 'fayb'
  nc_err = nf90_inq_varid(ncid, var_str, var_farsy_id)
  var_str = 'fazb'
  nc_err = nf90_inq_varid(ncid, var_str, var_farsz_id)

  DO n = 1, nblocks_subdom

    ! Get the grid coordinates

    m = subdom%blockids_compghst(n)
    subdom%blocks(n)%myrank = blockranks(m)
    subdom%blocks(n)%resx = resx(m)
    subdom%blocks(n)%resy = resy(m)
    subdom%blocks(n)%resz = resz(m) 

    start_read(:, :) = 1
    DO i = 1, m - 1
      nx = (ix_last(i) - ix_first(i) + 1) / 2 ** resx(i)
      ny = (jy_last(i) - jy_first(i) + 1) / 2 ** resy(i)
      nz = (kz_last(i) - kz_first(i) + 1) / 2 ** resz(i)
      start_read(1, 1) = start_read(1, 1) + nx + 1
      start_read(1, 2) = start_read(1, 2) + ny + 1       
      start_read(1, 3) = start_read(1, 3) + nz + 1
      start_read(1, 4) = start_read(1, 4) + nx * ny * nz
      start_read(1, 5) = start_read(1, 5) + (nx + 1) * ny * nz
      start_read(1, 6) = start_read(1, 6) + nx * (ny + 1) * nz
      start_read(1, 7) = start_read(1, 7) + nx * ny * (nz + 1)
    END DO

    nx = (ix_last(m) - ix_first(m) + 1) / 2 ** resx(m)
    ny = (jy_last(m) - jy_first(m) + 1) / 2 ** resy(m)
    nz = (kz_last(m) - kz_first(m) + 1) / 2 ** resz(m)

    ALLOCATE(subdom%blocks(n)%x2(nx + 1))
    ALLOCATE(subdom%blocks(n)%y2(ny + 1))
    ALLOCATE(subdom%blocks(n)%z2(nz + 1))

    nread(1) = nx + 1
    nc_err = nf90_get_var(ncid, var_x2_id, subdom%blocks(n)%x2, start_read(:, 1), nread)
    nread(1) = ny + 1
    nc_err = nf90_get_var(ncid, var_y2_id, subdom%blocks(n)%y2, start_read(:, 2), nread)
    nread(1) = nz + 1
    nc_err = nf90_get_var(ncid, var_z2_id, subdom%blocks(n)%z2, start_read(:, 3), nread)

    subdom%blocks(n)%blockid = m

    subdom%blocks(n)%fld_shape(1) = nz    
    subdom%blocks(n)%fld_shape(2) = ny
    subdom%blocks(n)%fld_shape(3) = nx
    subdom%blocks(n)%ncells = nz * ny * nx
    subdom%blocks(n)%nfaces = (nz + 1) * ny * nx + nz * (ny + 1) * nx + nz * ny * (nx + 1)

    ALLOCATE(subdom%blocks(n)%x(nx))
    ALLOCATE(subdom%blocks(n)%y(ny))
    ALLOCATE(subdom%blocks(n)%z(nz))
    ALLOCATE(subdom%blocks(n)%arsx(nz, ny, nx + 1))
    ALLOCATE(subdom%blocks(n)%arseffx(nz, ny, nx + 1))
    ALLOCATE(subdom%blocks(n)%arsy(nz, ny + 1, nx))
    ALLOCATE(subdom%blocks(n)%arseffy(nz, ny + 1, nx))
    ALLOCATE(subdom%blocks(n)%arsz(nz + 1, ny, nx))
    ALLOCATE(subdom%blocks(n)%arseffz(nz + 1, ny, nx))
    ALLOCATE(subdom%blocks(n)%vols(nz, ny, nx))
    ALLOCATE(subdom%blocks(n)%volseff(nz, ny, nx))
    ALLOCATE(subdom%blocks(n)%scal(nz, ny, nx))

    ALLOCATE(fvol_tmp(nx * ny * nz))
    ALLOCATE(fax_tmp((nx + 1) * ny * nz))
    ALLOCATE(fay_tmp(nx * (ny + 1) * nz))
    ALLOCATE(faz_tmp(nx * ny * (nz + 1)))

    subdom%ncells_tmp = subdom%ncells_tmp + SIZE(subdom%blocks(n)%vols)
    subdom%nfaces_tmp = subdom%nfaces_tmp + SIZE(subdom%blocks(n)%arsx)
    subdom%nfaces_tmp = subdom%nfaces_tmp + SIZE(subdom%blocks(n)%arsy)
    subdom%nfaces_tmp = subdom%nfaces_tmp + SIZE(subdom%blocks(n)%arsz)

    IF (subdom%blockiscomp(n) .EQV. .TRUE.) THEN
      subdom%ncells = subdom%ncells + SIZE(subdom%blocks(n)%vols)
      subdom%nfaces = subdom%nfaces + SIZE(subdom%blocks(n)%arsx)
      subdom%nfaces = subdom%nfaces + SIZE(subdom%blocks(n)%arsy)
      subdom%nfaces = subdom%nfaces + SIZE(subdom%blocks(n)%arsz)
    END IF

    subdom%blocks(n)%ncells = SIZE(subdom%blocks(n)%vols)
    subdom%blocks(n)%nfaces = SIZE(subdom%blocks(n)%arsx)
    subdom%blocks(n)%nfaces = subdom%blocks(n)%nfaces + SIZE(subdom%blocks(n)%arsy)
    subdom%blocks(n)%nfaces = subdom%blocks(n)%nfaces + SIZE(subdom%blocks(n)%arsz)

    subdom%blocks(n)%x(:) = 0.5 * (subdom%blocks(n)%x2(1:nx) + subdom%blocks(n)%x2(2:nx + 1))
    subdom%blocks(n)%y(:) = 0.5 * (subdom%blocks(n)%y2(1:ny) + subdom%blocks(n)%y2(2:ny + 1))
    subdom%blocks(n)%z(:) = 0.5 * (subdom%blocks(n)%z2(1:nz) + subdom%blocks(n)%z2(2:nz + 1))
    
    subdom%blocks(n)%i_st = 1
    subdom%blocks(n)%i_end = nx
    subdom%blocks(n)%j_st = 1
    subdom%blocks(n)%j_end = ny
    subdom%blocks(n)%k_st = 1
    subdom%blocks(n)%k_end = nz

    nread(1) = nx * ny * nz
    nc_err = nf90_get_var(ncid, var_fvol_id, fvol_tmp, start_read(:, 4), nread)
    nread(1) = (nx + 1) * ny * nz
    nc_err = nf90_get_var(ncid, var_farsx_id, fax_tmp, start_read(:, 5), nread)
    nread = nx * (ny + 1) * nz
    nc_err = nf90_get_var(ncid, var_farsy_id, fay_tmp, start_read(:, 6), nread)
    nread = (nz + 1) * ny * nx
    nc_err = nf90_get_var(ncid, var_farsz_id, faz_tmp, start_read(:, 7), nread)

    i1d = 1
    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx
          subdom%blocks(n)%volseff(k, j, i) = MAX(fvol_tmp(i1d), 1e-20)
          i1d = i1d + 1
        END DO
      END DO
    END DO

    i1d = 1
    DO k = 1, nz
      DO j = 1, ny
        DO i = 1, nx + 1
          subdom%blocks(n)%arseffx(k, j, i) = fax_tmp(i1d)
          i1d = i1d + 1
        END DO
      END DO
    END DO

    i1d = 1
    DO k = 1, nz
      DO j = 1, ny + 1
        DO i = 1, nx
          subdom%blocks(n)%arseffy(k, j, i) = fay_tmp(i1d)
          i1d = i1d + 1
        END DO
      END DO
    END DO

    i1d = 1
    DO k = 1, nz + 1
      DO j = 1, ny
        DO i = 1, nx
          subdom%blocks(n)%arseffz(k, j, i) = faz_tmp(i1d)
          i1d = i1d + 1
        END DO
      END DO
    END DO
 
    subdom%blocks(n)%scal(:,:,:) = 1.0

    DEALLOCATE(fvol_tmp)
    DEALLOCATE(fax_tmp)
    DEALLOCATE(fay_tmp)
    DEALLOCATE(faz_tmp)

    ! Compute the geometric and effective areas and volumes

    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz
          subdom%blocks(n)%volseff(k, j, i) = subdom%blocks(n)%volseff(k, j, i) * &
                                              (subdom%blocks(n)%x2(i + 1) - subdom%blocks(n)%x2(i)) * &
                                              (subdom%blocks(n)%y2(j + 1) - subdom%blocks(n)%y2(j)) * &
                                              (subdom%blocks(n)%z2(k + 1) - subdom%blocks(n)%z2(k))
          subdom%blocks(n)%vols(k, j, i) = (subdom%blocks(n)%x2(i + 1) - subdom%blocks(n)%x2(i)) * &
                                           (subdom%blocks(n)%y2(j + 1) - subdom%blocks(n)%y2(j)) * &
                                           (subdom%blocks(n)%z2(k + 1) - subdom%blocks(n)%z2(k))
        END DO

        DO k = 1, nz + 1
          subdom%blocks(n)%arseffz(k, j, i) = subdom%blocks(n)%arseffz(k, j, i) * &
                                              (subdom%blocks(n)%x2(i + 1) - subdom%blocks(n)%x2(i)) * &
                                              (subdom%blocks(n)%y2(j + 1) - subdom%blocks(n)%y2(j))
          subdom%blocks(n)%arsz(k, j, i) = (subdom%blocks(n)%x2(i + 1) - subdom%blocks(n)%x2(i)) * &
                                          (subdom%blocks(n)%y2(j + 1) - subdom%blocks(n)%y2(j))
        END DO
      END DO
      DO j = 1, ny + 1
        DO k = 1, nz
          subdom%blocks(n)%arseffy(k, j, i) = subdom%blocks(n)%arseffy(k, j, i) * &
                                              (subdom%blocks(n)%x2(i + 1) - subdom%blocks(n)%x2(i)) * &
                                              (subdom%blocks(n)%z2(k + 1) - subdom%blocks(n)%z2(k)) 
          subdom%blocks(n)%arsy(k, j, i) = (subdom%blocks(n)%x2(i + 1) - subdom%blocks(n)%x2(i)) * &
                                           (subdom%blocks(n)%z2(k + 1) - subdom%blocks(n)%z2(k))
        END DO
      END DO
    END DO
    DO i = 1, nx + 1
      DO j = 1, ny
        DO k = 1, nz
          subdom%blocks(n)%arseffx(k, j, i) = subdom%blocks(n)%arseffx(k, j, i) * &
                                              (subdom%blocks(n)%y2(j + 1) - subdom%blocks(n)%y2(j)) * &
                                              (subdom%blocks(n)%z2(k + 1) - subdom%blocks(n)%z2(k))
          subdom%blocks(n)%arsx(k, j, i) = (subdom%blocks(n)%y2(j + 1) - subdom%blocks(n)%y2(j)) * &
                                           (subdom%blocks(n)%z2(k + 1) - subdom%blocks(n)%z2(k))
        END DO
      END DO
    END DO

  END DO
  DEALLOCATE(is_in_tmp)
  DEALLOCATE(is_ghost)
  nc_err = nf90_close(ncid) 

  CALL check_common_edges(subdom%blocks, cyclic_x, cyclic_y, cyclic_z)
  CALL check_common_faces(subdom%blocks, cyclic_x, cyclic_y, cyclic_z)


END SUBROUTINE read_blocks_from_netcdf


SUBROUTINE check_common_edges(blocks, cyclic_x, cyclic_y, cyclic_z)

  IMPLICIT NONE

  TYPE(Block), POINTER, INTENT(INOUT) :: blocks(:)
  LOGICAL, INTENT(in), OPTIONAL :: cyclic_x, cyclic_y, cyclic_z
  INTEGER :: nz, ny, nx
  INTEGER :: nblocks, iblock, jblock
  LOGICAL :: cyclic_x_def, cyclic_y_def, cyclic_z_def

  IF (PRESENT(cyclic_x)) THEN
    cyclic_x_def = cyclic_x
  ELSE
    cyclic_x_def = .FALSE.
  END IF

  IF (PRESENT(cyclic_y)) THEN
    cyclic_y_def = cyclic_y
  ELSE
    cyclic_y_def = .FALSE.
  END IF

  IF (PRESENT(cyclic_z)) THEN
    cyclic_z_def = cyclic_z
  ELSE
    cyclic_z_def = .FALSE.
  END IF

  IF(.NOT. ASSOCIATED(blocks)) RETURN

  nblocks = SIZE(blocks)

  DO iblock = 1, nblocks

    nz = blocks(iblock)%fld_shape(1)
    ny = blocks(iblock)%fld_shape(2)
    nx = blocks(iblock)%fld_shape(3)
    ALLOCATE(blocks(iblock)%cedge_w(nblocks))
    ALLOCATE(blocks(iblock)%cedge_e(nblocks))
    ALLOCATE(blocks(iblock)%cedge_s(nblocks))
    ALLOCATE(blocks(iblock)%cedge_n(nblocks))
    ALLOCATE(blocks(iblock)%cedge_b(nblocks))
    ALLOCATE(blocks(iblock)%cedge_t(nblocks))

    blocks(iblock)%cedge_w(:) = .FALSE.
    blocks(iblock)%cedge_e(:) = .FALSE.
    blocks(iblock)%cedge_s(:) = .FALSE.
    blocks(iblock)%cedge_n(:) = .FALSE.
    blocks(iblock)%cedge_b(:) = .FALSE.
    blocks(iblock)%cedge_t(:) = .FALSE.

    DO jblock = 1, nblocks
    
      IF (jblock == iblock) CYCLE
      !check western boundary for a touching surface
      IF (blocks(iblock)%x2(1) == blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) THEN
        IF ((((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
            ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))))) &
            .OR. &
            (((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
            ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
          blocks(iblock)%cedge_w(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cedge_w(jblock) = .TRUE.
        END IF
      !check eastern boundary for a touching surface
      END IF
      IF (blocks(iblock)%x2(nx + 1) == blocks(jblock)%x2(1)) THEN
        IF ((((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
            ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))))) &
            .OR. &
            (((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
            ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
          blocks(iblock)%cedge_e(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cedge_e(jblock) = .TRUE.
        END IF
      END IF
      !check southern boundary for a touching surface
      IF (blocks(iblock)%y2(1) == blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) THEN
        IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
            ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
            .OR. &
            (((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
            ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
          blocks(iblock)%cedge_s(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cedge_s(jblock) = .TRUE.
        END IF
      END IF
      !check northern boundary for a touching surface
      IF (blocks(iblock)%y2(ny + 1) == blocks(jblock)%y2(1)) THEN
        IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
            ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
            .OR. &
            (((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
            ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
          blocks(iblock)%cedge_n(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cedge_n(jblock) = .TRUE.
        END IF
      END IF
      !check bottom boundary for a touching surface
      IF (blocks(iblock)%z2(1) == blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) THEN
        IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
            ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
            .OR. &
            (((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
            ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
          blocks(iblock)%cedge_b(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cedge_b(jblock) = .TRUE.
        END IF
      !check top boundary for a touching surface
      END IF
      IF (blocks(iblock)%z2(nz + 1) == blocks(jblock)%z2(1)) THEN
        IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
            ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
            .OR. &
            (((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
            ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
          blocks(iblock)%cedge_t(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cedge_t(jblock) = .TRUE.
        END IF
      END IF
    END DO
  END DO

  ! check for cyclic boundaries
  IF (cyclic_x_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
!        IF (jblock == iblock) CYCLE
        !check western boundary
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
            blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1))) THEN
          IF ((((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))))) &
              .OR. &
              (((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cedge_w(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_w(jblock) = .TRUE.
          END IF
        END IF
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1)) THEN
          IF ((((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))))) &
              .OR. &
              (((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cedge_e(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_e(jblock) = .TRUE.
          END IF
        END IF
      END DO
    END DO
  END IF

  ! check for cyclic boundaries
  IF (cyclic_y_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
!        IF (jblock == iblock) CYCLE
        !check western boundary
        IF (blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
            blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1))) THEN
          IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
              .OR. &
              (((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cedge_s(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_s(jblock) = .TRUE.
          END IF
        END IF
        IF (blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(jblock)%y2(1) .EQ. y2glob(1)) THEN
          IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
              .OR. &
              (((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cedge_n(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_n(jblock) = .TRUE.
          END IF
        END IF
      END DO
    END DO
  END IF


  ! check for cyclic boundaries
  IF (cyclic_z_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
!        IF (jblock == iblock) CYCLE
        !check bottom boundary
        IF (blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
           blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
          IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
              .OR. &
              (((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
            blocks(iblock)%cedge_b(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_b(jblock) = .TRUE.
          END IF
        END IF
        IF (blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
                 blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
          IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
              .OR. &
              (((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
            blocks(iblock)%cedge_t(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_t(jblock) = .TRUE.
          END IF
        END IF
      END DO
    END DO
  END IF

  IF (cyclic_x_def .AND. cyclic_y_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
!        IF (jblock == iblock) CYCLE
        !check southwestern edge
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
           blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
           blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
           blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1))) THEN
          IF ((((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cedge_w(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_w(jblock) = .TRUE.
          END IF
        END IF
        !check northeastern edge
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(jblock)%y2(1) .EQ. y2glob(1)) THEN
          IF ((((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cedge_e(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_e(jblock) = .TRUE.
          END IF
        !check northwestern edge
        ELSE IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(jblock)%y2(1) .EQ. y2glob(1)) THEN
          IF ((((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cedge_n(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_n(jblock) = .TRUE.
          END IF
        END IF
        !check southeastern edge
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
                 blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1))) THEN
          IF ((((blocks(iblock)%z2(1) .lt. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .lt. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .gt. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cedge_s(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_s(jblock) = .TRUE.
          END IF
        END IF
      END DO
    END DO
  END IF

  IF (cyclic_x_def .AND. cyclic_z_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
!        IF (jblock == iblock) CYCLE
        !check bottomwestern edge
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
           blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
           blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
           blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
          IF ((((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
            blocks(iblock)%cedge_w(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_w(jblock) = .TRUE.
          END IF
        END IF
        !check bottomeastern edge
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
                 blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
          IF ((((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
            blocks(iblock)%cedge_e(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_e(jblock) = .TRUE.
          END IF
        END IF
        !check topwestern edge
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
           blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
           blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ.  z2glob(UBOUND(z2glob, 1)) .AND. &
           blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
          IF ((((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
            blocks(iblock)%cedge_w(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_w(jblock) = .TRUE.
          END IF
        END IF
        !check topeastern edge
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
                 blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
          IF ((((blocks(iblock)%y2(1) .lt. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .lt. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .gt. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
            blocks(iblock)%cedge_e(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_e(jblock) = .TRUE.
          END IF
        END IF
      END DO
    END DO
  END IF

  IF (cyclic_y_def .AND. cyclic_z_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
!        IF (jblock == iblock) CYCLE
        !check bottomsouthern edge
        IF (blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
           blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
           blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
           blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
          IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2)))))) THEN
            blocks(iblock)%cedge_s(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_s(jblock) = .TRUE.
          END IF
        END IF
        !check bottomnorthern edge
        IF (blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
           blocks(jblock)%y2(1) .EQ. y2glob(1) .AND. &
           blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
           blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
          IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2)))))) THEN
            blocks(iblock)%cedge_n(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_n(jblock) = .TRUE.
          END IF
        END IF
        !check topsouthern edge
        IF (blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
           blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
           blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
           blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
          IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2)))))) THEN
            blocks(iblock)%cedge_s(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_s(jblock) = .TRUE.
          END IF
        END IF
        !check topnorthern edge
        IF (blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
           blocks(jblock)%y2(1) .EQ. y2glob(1) .AND. &
           blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
           blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
          IF ((((blocks(iblock)%x2(1) .lt. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .lt. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .gt. blocks(jblock)%x2(SIZE(blocks(jblock)%x2)))))) THEN
            blocks(iblock)%cedge_n(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cedge_n(jblock) = .TRUE.
          END IF
        END IF
      END DO
    END DO
  END IF

  IF (cyclic_x_def .AND. cyclic_y_def .AND. cyclic_z_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
!        IF (jblock == iblock) CYCLE
        !bottom-southwestern corner
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
           blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
           blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
           blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
           blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
           blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
           blocks(iblock)%cedge_w(jblock) = .TRUE.
        END IF
        !bottom-southeastern corner
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
                 blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
                 blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
           blocks(iblock)%cedge_e(jblock) = .TRUE.
        END IF
        !bottom-northwestern corner
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(jblock)%y2(1) .EQ. y2glob(1) .AND. &
                 blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
                 blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
                 blocks(iblock)%cedge_w(jblock) = .TRUE.
        END IF
        !top-southwestern corner
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
                 blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
                 blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
                 blocks(iblock)%cedge_w(jblock) = .TRUE. 
        END IF
        !bottom-northeastern corner
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(jblock)%y2(1) .EQ. y2glob(1) .AND. &
                 blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
                 blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
           blocks(iblock)%cedge_e(jblock) = .TRUE.
        END IF
        !top-northwestern corner
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(jblock)%y2(1) .EQ. y2glob(1) .AND. &
                 blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
                 blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
                 blocks(iblock)%cedge_w(jblock) = .TRUE.
        !top-southeastern corner
        END IF
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
                 blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
                 blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
                 blocks(iblock)%cedge_e(jblock) = .TRUE.
        END IF
        !top-northeastern corner
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1) .AND. &
                 blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(jblock)%y2(1) .EQ. y2glob(1) .AND. &
                 blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
                 blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
                 blocks(iblock)%cedge_e(jblock) = .TRUE.
        END IF
      END DO
    END DO
  END IF

END SUBROUTINE check_common_edges


SUBROUTINE check_common_faces(blocks, cyclic_x, cyclic_y, cyclic_z)

  IMPLICIT NONE

  TYPE(Block), POINTER, INTENT(INOUT) :: blocks(:)
  LOGICAL, INTENT(in), OPTIONAL :: cyclic_x, cyclic_y, cyclic_z
  INTEGER :: nz, ny, nx
  INTEGER :: nblocks, iblock, jblock
  LOGICAL :: cyclic_x_def, cyclic_y_def, cyclic_z_def

  IF (PRESENT(cyclic_x)) THEN
    cyclic_x_def = cyclic_x
  ELSE
    cyclic_x_def = .FALSE.
  END IF

  IF (PRESENT(cyclic_y)) THEN
    cyclic_y_def = cyclic_y
  ELSE
    cyclic_y_def = .FALSE.
  END IF

  IF (PRESENT(cyclic_z)) THEN
    cyclic_z_def = cyclic_z
  ELSE
    cyclic_z_def = .FALSE.
  END IF

  IF(.NOT. ASSOCIATED(blocks)) RETURN

  nblocks = SIZE(blocks)

  DO iblock = 1, nblocks

    nz = blocks(iblock)%fld_shape(1)
    ny = blocks(iblock)%fld_shape(2)
    nx = blocks(iblock)%fld_shape(3)

    ALLOCATE(blocks(iblock)%cface_w(nblocks))
    ALLOCATE(blocks(iblock)%cface_e(nblocks))
    ALLOCATE(blocks(iblock)%cface_s(nblocks))
    ALLOCATE(blocks(iblock)%cface_n(nblocks))
    ALLOCATE(blocks(iblock)%cface_b(nblocks))
    ALLOCATE(blocks(iblock)%cface_t(nblocks))

    blocks(iblock)%cface_w(:) = .FALSE.
    blocks(iblock)%cface_e(:) = .FALSE.
    blocks(iblock)%cface_s(:) = .FALSE.
    blocks(iblock)%cface_n(:) = .FALSE.
    blocks(iblock)%cface_b(:) = .FALSE.
    blocks(iblock)%cface_t(:) = .FALSE.

    DO jblock = 1, nblocks

      IF (jblock == iblock) CYCLE
      !check western boundary for a touching surface
      IF (blocks(iblock)%x2(1) == blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) THEN
        IF ((((blocks(iblock)%y2(1) .le. blocks(jblock)%y2(1)) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .le. blocks(jblock)%y2(1))) .OR. &
            ((blocks(iblock)%y2(1) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))))) &
            .OR. &
            (((blocks(iblock)%z2(1) .le. blocks(jblock)%z2(1)) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .le. blocks(jblock)%z2(1))) .OR. &
            ((blocks(iblock)%z2(1) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
          blocks(iblock)%cface_w(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cface_w(jblock) = .TRUE.
        END IF
      !check eastern boundary for a touching surface
      ELSEIF (blocks(iblock)%x2(nx + 1) == blocks(jblock)%x2(1)) THEN
        IF ((((blocks(iblock)%y2(1) .le. blocks(jblock)%y2(1)) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .le. blocks(jblock)%y2(1))) .OR. &
            ((blocks(iblock)%y2(1) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))))) &
            .OR. &
            (((blocks(iblock)%z2(1) .le. blocks(jblock)%z2(1)) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .le. blocks(jblock)%z2(1))) .OR. &
            ((blocks(iblock)%z2(1) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
          blocks(iblock)%cface_e(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cface_e(jblock) = .TRUE.
        END IF
      !check southern boundary for a touching surface
      ELSEIF (blocks(iblock)%y2(1) == blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) THEN
        IF ((((blocks(iblock)%x2(1) .le. blocks(jblock)%x2(1)) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .le. blocks(jblock)%x2(1))) .OR. &
            ((blocks(iblock)%x2(1) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
            .OR. &
            (((blocks(iblock)%z2(1) .le. blocks(jblock)%z2(1)) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .le. blocks(jblock)%z2(1))) .OR. &
            ((blocks(iblock)%z2(1) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
          blocks(iblock)%cface_s(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cface_s(jblock) = .TRUE.
        END IF
      !check northern boundary for a touching surface
      ELSE IF (blocks(iblock)%y2(ny + 1) == blocks(jblock)%y2(1)) THEN
        IF ((((blocks(iblock)%x2(1) .le. blocks(jblock)%x2(1)) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .le. blocks(jblock)%x2(1))) .OR. &
            ((blocks(iblock)%x2(1) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
            .OR. &
            (((blocks(iblock)%z2(1) .le. blocks(jblock)%z2(1)) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .le. blocks(jblock)%z2(1))) .OR. &
            ((blocks(iblock)%z2(1) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
              (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
          blocks(iblock)%cface_n(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cface_n(jblock) = .TRUE.
        END IF
      !check bottom boundary for a touching surface
      ELSEIF (blocks(iblock)%z2(1) == blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) THEN
        IF ((((blocks(iblock)%x2(1) .le. blocks(jblock)%x2(1)) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .le. blocks(jblock)%x2(1))) .OR. &
            ((blocks(iblock)%x2(1) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
            .OR. &
            (((blocks(iblock)%y2(1) .le. blocks(jblock)%y2(1)) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .le. blocks(jblock)%y2(1))) .OR. &
            ((blocks(iblock)%y2(1) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
          blocks(iblock)%cface_b(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cface_b(jblock) = .TRUE.
        END IF
      !check top boundary for a touching surface
      ELSE IF (blocks(iblock)%z2(nz + 1) == blocks(jblock)%z2(1)) THEN
        IF ((((blocks(iblock)%x2(1) .le. blocks(jblock)%x2(1)) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .le. blocks(jblock)%x2(1))) .OR. &
            ((blocks(iblock)%x2(1) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
              (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
            .OR. &
            (((blocks(iblock)%y2(1) .le. blocks(jblock)%y2(1)) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .le. blocks(jblock)%y2(1))) .OR. &
            ((blocks(iblock)%y2(1) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
              (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
          blocks(iblock)%cface_t(jblock) = .FALSE.
        ELSE
          blocks(iblock)%cface_t(jblock) = .TRUE.
        END IF
      END IF
    END DO
  END DO

  ! check for cyclic boundaries in x-direction
  IF (cyclic_x_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
        !check western boundary
        IF (blocks(iblock)%x2(1) .EQ. x2glob(1) .AND. &
            blocks(jblock)%x2(blocks(jblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1))) THEN
          IF ((((blocks(iblock)%y2(1) .le. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .le. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))))) &
              .OR. &
              (((blocks(iblock)%z2(1) .le. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .le. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cface_w(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cface_w(jblock) = .TRUE.
          END IF
        END IF
        !check eastern boundary
        IF (blocks(iblock)%x2(blocks(iblock)%fld_shape(3) + 1) .EQ. x2glob(UBOUND(x2glob, 1)) .AND. &
                 blocks(jblock)%x2(1) .EQ. x2glob(1)) THEN
          IF ((((blocks(iblock)%y2(1) .le. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .le. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))))) &
              .OR. &
              (((blocks(iblock)%z2(1) .le. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .le. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cface_e(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cface_e(jblock) = .TRUE.
          END IF
        END IF 
      END DO 
    END DO
  END IF

  ! check for cyclic boundaries in y-direction
  IF (cyclic_y_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
        !check southern boundary
        IF (blocks(iblock)%y2(1) .EQ. y2glob(1) .AND. &
            blocks(jblock)%y2(blocks(jblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1))) THEN
          IF ((((blocks(iblock)%x2(1) .le. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .le. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
              .OR. &
              (((blocks(iblock)%z2(1) .le. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .le. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cface_s(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cface_s(jblock) = .TRUE.
          END IF
        END IF
        !check northern boundary
        IF (blocks(iblock)%y2(blocks(iblock)%fld_shape(2) + 1) .EQ. y2glob(UBOUND(y2glob, 1)) .AND. &
                 blocks(jblock)%y2(1) .EQ. y2glob(1)) THEN
          IF ((((blocks(iblock)%x2(1) .le. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .le. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
              .OR. &
              (((blocks(iblock)%z2(1) .le. blocks(jblock)%z2(1)) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .le. blocks(jblock)%z2(1))) .OR. &
              ((blocks(iblock)%z2(1) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2))) .AND. &
                (blocks(iblock)%z2(SIZE(blocks(iblock)%z2)) .ge. blocks(jblock)%z2(SIZE(blocks(jblock)%z2)))))) THEN
            blocks(iblock)%cface_n(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cface_n(jblock) = .TRUE.
          END IF
        END IF
      END DO 
    END DO
  END IF

  ! check for cyclic boundaries in z-direction
  IF (cyclic_z_def) THEN
    DO iblock = 1, nblocks

      nz = blocks(iblock)%fld_shape(1)
      ny = blocks(iblock)%fld_shape(2)
      nx = blocks(iblock)%fld_shape(3)

      DO jblock = 1, nblocks
        !check bottom boundary
        IF (blocks(iblock)%z2(1) .EQ. z2glob(1) .AND. &
            blocks(jblock)%z2(blocks(jblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1))) THEN
          IF ((((blocks(iblock)%x2(1) .le. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .le. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
              .OR. &
              (((blocks(iblock)%y2(1) .le. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .le. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
            blocks(iblock)%cface_b(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cface_b(jblock) = .TRUE.
          END IF
        END IF
        !check top boundary
        IF (blocks(iblock)%z2(blocks(iblock)%fld_shape(1) + 1) .EQ. z2glob(UBOUND(z2glob, 1)) .AND. &
                 blocks(jblock)%z2(1) .EQ. z2glob(1)) THEN
          IF ((((blocks(iblock)%x2(1) .le. blocks(jblock)%x2(1)) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .le. blocks(jblock)%x2(1))) .OR. &
              ((blocks(iblock)%x2(1) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))) .AND. &
                (blocks(iblock)%x2(SIZE(blocks(iblock)%x2)) .ge. blocks(jblock)%x2(SIZE(blocks(jblock)%x2))))) &
              .OR. &
              (((blocks(iblock)%y2(1) .le. blocks(jblock)%y2(1)) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .le. blocks(jblock)%y2(1))) .OR. &
              ((blocks(iblock)%y2(1) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2))) .AND. &
                (blocks(iblock)%y2(SIZE(blocks(iblock)%y2)) .ge. blocks(jblock)%y2(SIZE(blocks(jblock)%y2)))))) THEN
            blocks(iblock)%cface_t(jblock) = .FALSE.
          ELSE
            blocks(iblock)%cface_t(jblock) = .TRUE.
          END IF
        END IF
      END DO
    END DO
  END IF

END SUBROUTINE check_common_faces


SUBROUTINE vec1d_from_homog(subdom, u, v, w, x)
  IMPLICIT NONE
  TYPE(SubDomain), INTENT(in) :: subdom 
  REAL(Realkind), INTENT(in) :: u, v, w
  REAL(Realkind), INTENT(inout) :: x(:)
  REAL(Realkind), ALLOCATABLE :: u_fld(:, :, :), v_fld(:, :, :), w_fld(:, :, :)
  INTEGER :: iblock
  INTEGER :: nz, ny, nx
  INTEGER :: current_ind


  current_ind = 1
  DO iblock = 1, subdom%nblocks
    IF (subdom%blockiscomp(iblock) .EQV. .TRUE.) THEN      
      nx = subdom%blocks(iblock)%fld_shape(3)
      ny = subdom%blocks(iblock)%fld_shape(2)
      nz = subdom%blocks(iblock)%fld_shape(1)
      ALLOCATE(u_fld(nz, ny, nx + 1))
      ALLOCATE(v_fld(nz, ny + 1, nx))
      ALLOCATE(w_fld(nz + 1, ny, nx))

      u_fld(:, :, :) = u
      v_fld(:, :, :) = v
      w_fld(:, :, :) = w

      IF (iblock .EQ. 1) u_fld(1:1, 1:1, 5:10) = 1.0

      x(current_ind:current_ind + SIZE(u_fld) - 1) = PACK(u_fld, .true.)
      current_ind = current_ind + SIZE(u_fld)
      x(current_ind:current_ind + SIZE(v_fld) - 1) = PACK(v_fld, .true.)
      current_ind = current_ind + SIZE(v_fld)
      x(current_ind:current_ind + SIZE(w_fld) - 1) = PACK(w_fld, .true.)
      current_ind = current_ind + SIZE(w_fld)

      DEALLOCATE(u_fld)
      DEALLOCATE(v_fld)
      DEALLOCATE(w_fld)
    END IF
  END DO
  
END SUBROUTINE


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

END MODULE DomDecomp_Mod
