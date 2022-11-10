MODULE CoarseGrids_Mod

USE Subdom_Mod, ONLY: SubDomain, set_cell_type, set_face_type
USE DomDecomp_Mod, ONLY: x2glob, y2glob, z2glob, &
                         ix_first, ix_last, jy_first, jy_last, kz_first, kz_last, &
                         resx, resy, resz, btype, blockranks, &
                         check_common_faces, check_common_edges


USE MatSpRowCol_Mod, ONLY: SpRowCol, SpMm_SpRowCol, SpTrans_SpRowCol, &
                           SpDeallocate_SpRowCol, SpNullify_SpRowCol
USE Crop_Mod, ONLY: make_crop_operator_cell, make_crop_operator_face
USE Block_Mod, ONLY: Block, Regular_Grid, gridding_cell, grid_blocks, neighborblocks
USE List_Mod, ONLY : List, Element, ListArray, IArray, RArray, IArray2d, RArray
USE Sort_Mod, ONLY: sortunique
USE Transfer_Mod, ONLY: init_transfer
USE Laplace_block_Mod, ONLY: make_operators
USE Index_Mod, ONLY: findloc_real, between, determine_overlap_1d
USE Gradient_block_interpolation_Mod, ONLY: find_indices_jk, find_indices_ik, find_indices_ij

USE Kind_Mod
USE MPI_Mod

IMPLICIT NONE

TYPE(RArray), ALLOCATABLE :: x2glob_lev(:)
TYPE(RArray), ALLOCATABLE :: y2glob_lev(:)
TYPE(RArray), ALLOCATABLE :: z2glob_lev(:)

TYPE(IArray), ALLOCATABLE :: ix_st_lev(:)
TYPE(IArray), ALLOCATABLE :: ix_end_lev(:)
TYPE(IArray), ALLOCATABLE :: jy_st_lev(:)
TYPE(IArray), ALLOCATABLE :: jy_end_lev(:)
TYPE(IArray), ALLOCATABLE :: kz_st_lev(:)
TYPE(IArray), ALLOCATABLE :: kz_end_lev(:)

TYPE(IArray), ALLOCATABLE :: blockranks_lev(:)
TYPE(IArray), ALLOCATABLE :: blockglobids_lev(:)

TYPE(IArray), ALLOCATABLE :: resx_lev(:)
TYPE(IArray), ALLOCATABLE :: resy_lev(:)
TYPE(IArray), ALLOCATABLE :: resz_lev(:)

TYPE(IArray), ALLOCATABLE :: isblackijkzero_lev(:)

TYPE(IArray), ALLOCATABLE :: btype_lev(:)

TYPE(Subdomain), POINTER :: subdomain_lev
INTEGER :: ngrids

REAL(Realkind), PARAMETER :: w = 0.55

CONTAINS

SUBROUTINE define_coarse_grids()
  IMPLICIT NONE

  INTEGER :: nx_fine, ny_fine, nz_fine, dim_size_max
  INTEGER :: nblocks_finest, i

  nblocks_finest = SIZE(ix_first, 1)

  nx_fine = SIZE(x2glob, 1) - 1
  ny_fine = SIZE(y2glob, 1) - 1
  nz_fine = SIZE(z2glob, 1) - 1

  dim_size_max = MAXVAL((/nx_fine, ny_fine, nz_fine/))
 
  ngrids = INT(LOG(REAL(dim_size_max, kind=Realkind)) / LOG(2.0))

  ALLOCATE(x2glob_lev(ngrids))
  ALLOCATE(y2glob_lev(ngrids))
  ALLOCATE(z2glob_lev(ngrids))
  ALLOCATE(ix_st_lev(ngrids))
  ALLOCATE(ix_end_lev(ngrids))
  ALLOCATE(jy_st_lev(ngrids))
  ALLOCATE(jy_end_lev(ngrids))  
  ALLOCATE(kz_st_lev(ngrids))
  ALLOCATE(kz_end_lev(ngrids))
  ALLOCATE(blockranks_lev(ngrids))
  ALLOCATE(blockglobids_lev(ngrids))
  ALLOCATE(resx_lev(ngrids))
  ALLOCATE(resy_lev(ngrids))
  ALLOCATE(resz_lev(ngrids))
  ALLOCATE(btype_lev(ngrids))

  CALL coarsen_1d(x2glob, x2glob_lev, ix_first, ix_last, resx, btype)
  CALL coarsen_1d(y2glob, y2glob_lev, jy_first, jy_last, resy, btype)
  CALL coarsen_1d(z2glob, z2glob_lev, kz_first, kz_last, resz, btype)

END SUBROUTINE define_coarse_grids


SUBROUTINE coarsen_1d(x2, x2_glob_lst, ix_f, ix_l, res, bt)

  IMPLICIT NONE

  TYPE(RArray), INTENT(inout) :: x2_glob_lst(:)
  REAL(Realkind), INTENT(in) :: x2(:)
  INTEGER, INTENT(in) :: ix_f(:), ix_l(:)
  INTEGER, INTENT(in) :: res(:)
  INTEGER, INTENT(in) :: bt(:)

  INTEGER :: iblock, nblocks, iwork
  
  REAL(Realkind) :: coarsen_fac
  INTEGER :: i, j, n
  REAL(Realkind) :: work(SIZE(res, 1) * (SIZE(x2, 1) + 2))
  REAL(Realkind) :: line_prev, line_next
  INTEGER :: iprev, inext

  nblocks = SIZE(res, 1)

  work(1) = x2(1) 
  work(2) = x2(UBOUND(x2, 1))
  iwork = 2

  DO iblock = 1, nblocks
    IF (bt(iblock) .EQ. 0) THEN
      work(iwork + 1) = x2(ix_f(iblock))
      work(iwork + 2) = x2(ix_l(iblock) + 1)
      iwork = iwork + 2
    END IF
  END DO

  DO n = ngrids, 2, -1
    coarsen_fac = 2 ** (n - 1)
    DO iblock = 1, nblocks
      IF (res(iblock) .EQ. n - 1) THEN
        DO i = ix_f(iblock), ix_l(iblock), coarsen_fac
          work(iwork + 1) = x2(i)
          iwork = iwork + 1
        END DO
      END IF
    END DO
    CALL sortunique(work(1:iwork), x2_glob_lst(n)%data)
  END DO

  ALLOCATE(x2_glob_lst(1)%data(SIZE(x2, 1)))
  x2_glob_lst(1)%data = x2

  DO n = 2, ngrids
    iwork = SIZE(x2_glob_lst(n)%data, 1)
    work(1:iwork) = x2_glob_lst(n)%data
    line_prev = work(1)
    iprev = findloc_real(x2_glob_lst(n - 1)%data, line_prev)
    DO i = 2, SIZE(x2_glob_lst(n)%data, 1)
      line_next = x2_glob_lst(n)%data(i)
      inext = findloc_real(x2_glob_lst(n - 1)%data, line_next)
      DO j = iprev + 2, inext, 2
        work(iwork + 1) = x2_glob_lst(n - 1)%data(j)
        iwork = iwork + 1
     END DO
      iprev = inext
    END DO
    CALL sortunique(work(1:iwork), x2_glob_lst(n)%data)
  END DO
   
END SUBROUTINE coarsen_1d


SUBROUTINE partition_clevels()
  IMPLICIT NONE

  INTEGER :: i, nblocks_fine

  nblocks_fine = SIZE(resx, 1)

  ALLOCATE(resx_lev(1)%data(nblocks_fine))
  ALLOCATE(resy_lev(1)%data(nblocks_fine))
  ALLOCATE(resz_lev(1)%data(nblocks_fine))
  ALLOCATE(ix_st_lev(1)%data(nblocks_fine))
  ALLOCATE(ix_end_lev(1)%data(nblocks_fine))
  ALLOCATE(jy_st_lev(1)%data(nblocks_fine))
  ALLOCATE(jy_end_lev(1)%data(nblocks_fine))
  ALLOCATE(kz_st_lev(1)%data(nblocks_fine))
  ALLOCATE(kz_end_lev(1)%data(nblocks_fine))
  ALLOCATE(blockglobids_lev(1)%data(nblocks_fine))
  ALLOCATE(blockranks_lev(1)%data(nblocks_fine))
  ALLOCATE(btype_lev(1)%data(nblocks_fine))
  
  resx_lev(1)%data(:) = resx
  resy_lev(1)%data(:) = resy
  resz_lev(1)%data(:) = resz
 
  ix_st_lev(1)%data(:) = ix_first
  ix_end_lev(1)%data(:) = ix_last
  jy_st_lev(1)%data(:) = jy_first
  jy_end_lev(1)%data(:) = jy_last
  kz_st_lev(1)%data(:) = kz_first
  kz_end_lev(1)%data(:) = kz_last
  blockranks_lev(1)%data(:) = blockranks
  btype_lev(1)%data(:) = btype

  DO i = 1, nblocks_fine
    blockglobids_lev(1)%data(i) = i 
  END DO

  DO i = 2, ngrids
    CALL partition_clevel_similar_flevel(i)
  END DO

END SUBROUTINE partition_clevels


SUBROUTINE partition_clevel_similar_flevel(levelc)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: levelc
  INTEGER :: iblock, nblocks, n

  INTEGER, DIMENSION(:), POINTER :: blockranks_f, blockglobids_f, ix_st_f, ix_end_f, jy_st_f, jy_end_f, &
                                    kz_st_f, kz_end_f, resx_f, resy_f, resz_f, btype_f

  REAL(Realkind), DIMENSION(:), POINTER :: x2glob_f, y2glob_f, z2glob_f, x2glob_c, y2glob_c, z2glob_c
  REAL(Realkind) :: x2_st_f, x2_end_f, y2_st_f, y2_end_f, z2_st_f, z2_end_f

  INTEGER :: ix_stc, ix_endc, jy_stc, jy_endc, kz_stc, kz_endc

  INTEGER :: work(12, SIZE(resx_lev(levelc - 1)%data, 1))

  blockranks_f => blockranks_lev(levelc - 1)%data
  blockglobids_f => blockglobids_lev(levelc - 1)%data

  ix_st_f => ix_st_lev(levelc - 1)%data
  ix_end_f => ix_end_lev(levelc - 1)%data
  jy_st_f => jy_st_lev(levelc - 1)%data
  jy_end_f => jy_end_lev(levelc - 1)%data
  kz_st_f => kz_st_lev(levelc - 1)%data
  kz_end_f => kz_end_lev(levelc - 1)%data

  resx_f => resx_lev(levelc - 1)%data
  resy_f => resy_lev(levelc - 1)%data
  resz_f => resz_lev(levelc - 1)%data

  btype_f => btype_lev(levelc - 1)%data

  x2glob_f => x2glob_lev(levelc - 1)%data
  y2glob_f => y2glob_lev(levelc - 1)%data
  z2glob_f => z2glob_lev(levelc - 1)%data

  x2glob_c => x2glob_lev(levelc)%data
  y2glob_c => y2glob_lev(levelc)%data
  z2glob_c => z2glob_lev(levelc)%data

  nblocks = SIZE(ix_st_f, 1)

  n = 0

  DO iblock = 1, nblocks
    x2_st_f = x2glob_f(ix_st_f(iblock))
    x2_end_f = x2glob_f(ix_end_f(iblock) + 1)
    y2_st_f = y2glob_f(jy_st_f(iblock))
    y2_end_f = y2glob_f(jy_end_f(iblock) + 1)
    z2_st_f = z2glob_f(kz_st_f(iblock))
    z2_end_f = z2glob_f(kz_end_f(iblock) + 1)

    ix_stc = MINLOC(ABS(x2glob_c - x2_st_f), 1)
    ix_endc = MINLOC(ABS(x2glob_c - x2_end_f), 1) - 1
    jy_stc = MINLOC(ABS(y2glob_c - y2_st_f), 1)
    jy_endc = MINLOC(ABS(y2glob_c - y2_end_f), 1) - 1
    kz_stc = MINLOC(ABS(z2glob_c - z2_st_f), 1)
    kz_endc = MINLOC(ABS(z2glob_c - z2_end_f), 1) - 1


    IF (ANY((/ix_stc .GT. ix_endc, jy_stc .GT. jy_endc, kz_stc .GT. kz_endc/))) CYCLE

    n = n + 1

    work(1, n) = MAX(resx_f(iblock), levelc - 1)
    work(2, n) = MAX(resy_f(iblock), levelc - 1)
    work(3, n) = MAX(resz_f(iblock), levelc - 1) 
    work(4, n) = btype_f(iblock)
    work(5, n) = ix_stc
    work(6, n) = ix_endc
    work(7, n) = jy_stc
    work(8, n) = jy_endc
    work(9, n) = kz_stc
    work(10, n) = kz_endc
    work(11, n) = blockglobids_f(iblock)
    work(12, n) = blockranks_f(iblock)

  END DO
  
  ALLOCATE(resx_lev(levelc)%data(n))
  ALLOCATE(resy_lev(levelc)%data(n))
  ALLOCATE(resz_lev(levelc)%data(n))
  ALLOCATE(btype_lev(levelc)%data(n))
  ALLOCATE(ix_st_lev(levelc)%data(n))
  ALLOCATE(ix_end_lev(levelc)%data(n))
  ALLOCATE(jy_st_lev(levelc)%data(n))
  ALLOCATE(jy_end_lev(levelc)%data(n))
  ALLOCATE(kz_st_lev(levelc)%data(n))
  ALLOCATE(kz_end_lev(levelc)%data(n))
  ALLOCATE(blockglobids_lev(levelc)%data(n))
  ALLOCATE(blockranks_lev(levelc)%data(n))

  resx_lev(levelc)%data(:) = work(1, 1:n)
  resy_lev(levelc)%data(:) = work(2, 1:n)
  resz_lev(levelc)%data(:) = work(3, 1:n)
  btype_lev(levelc)%data(:) = work(4, 1:n)
  ix_st_lev(levelc)%data(:) = work(5, 1:n)
  ix_end_lev(levelc)%data(:) = work(6, 1:n)
  jy_st_lev(levelc)%data(:) = work(7, 1:n)
  jy_end_lev(levelc)%data(:) = work(8, 1:n)
  kz_st_lev(levelc)%data(:) = work(9, 1:n)
  kz_end_lev(levelc)%data(:) = work(10, 1:n)
  blockglobids_lev(levelc)%data(:) = work(11, 1:n)
  blockranks_lev(levelc)%data(:) = work(12, 1:n)


END SUBROUTINE partition_clevel_similar_flevel


SUBROUTINE intialize_coarse_subdomains(subdomain_finest)

  TYPE(Subdomain), TARGET, INTENT(inout) :: subdomain_finest
  TYPE(SubDomain), POINTER :: subdomain_fine
  TYPE(SubDomain), POINTER :: newsubdomain => NULL()
  INTEGER, DIMENSION(:), POINTER :: blockranks_f, blockranks, &
                                    ix_first, ix_last, jy_first, jy_last, kz_first, kz_last, &
                                    ix_first_f, ix_last_f, jy_first_f, jy_last_f, kz_first_f, kz_last_f, &
                                    btype, resx, resy, resz, resx_f, resy_f, resz_f, blockglobids, blockglobids_f
  REAL(Realkind), DIMENSION(:), POINTER :: x2glob, y2glob, z2glob, x2glob_f, y2glob_f, z2glob_f

  TYPE(SpRowCol) :: OP_crop_cell, OP_crop_face, OP_crop_cell_T, OP_crop_face_T, &
                    grad3d, div3d, lapl3d, intc2f, op_tmp

  INTEGER :: n, i, j
  INTEGER :: nblocks_total


  subdomain_lev => subdomain_finest
  subdomain_fine => subdomain_finest

  DO n = 2, ngrids

    ALLOCATE(subdomain_fine%next_coarse)
    newsubdomain => subdomain_fine%next_coarse

    newsubdomain%next_fine => subdomain_fine

    blockranks => blockranks_lev(n)%data
    blockranks_f => blockranks_lev(n - 1)%data

    nblocks_total = SIZE(blockranks, 1)

    ALLOCATE(newsubdomain%blockids_comp(COUNT(blockranks .EQ. subdomain_fine%myrank, 1)))
   
    j = 0
    DO i = 1, nblocks_total
      IF (blockranks(i) .EQ. subdomain_fine%myrank) THEN
        j = j + 1
        newsubdomain%blockids_comp(j) = i
      END IF
    END DO

    x2glob => x2glob_lev(n)%data
    y2glob => y2glob_lev(n)%data
    z2glob => z2glob_lev(n)%data
    x2glob_f => x2glob_lev(n - 1)%data
    y2glob_f => y2glob_lev(n - 1)%data
    z2glob_f => z2glob_lev(n - 1)%data
   
    ix_first => ix_st_lev(n)%data
    ix_last => ix_end_lev(n)%data
    jy_first => jy_st_lev(n)%data
    jy_last => jy_end_lev(n)%data
    kz_first => kz_st_lev(n)%data
    kz_last => kz_end_lev(n)%data

    ix_first_f => ix_st_lev(n - 1)%data
    ix_last_f => ix_end_lev(n - 1)%data
    jy_first_f => jy_st_lev(n - 1)%data
    jy_last_f => jy_end_lev(n - 1)%data
    kz_first_f => kz_st_lev(n - 1)%data
    kz_last_f => kz_end_lev(n - 1)%data  

    resx => resx_lev(n)%data
    resy => resy_lev(n)%data
    resz => resz_lev(n)%data

    btype => btype_lev(n)%data

    resx_f => resx_lev(n - 1)%data
    resy_f => resy_lev(n - 1)%data
    resz_f => resz_lev(n - 1)%data

    blockglobids => blockglobids_lev(n)%data
    blockglobids_f => blockglobids_lev(n - 1)%data

    CALL init_newblocks(newsubdomain, subdomain_fine, n, nblocks_total, &
                          x2glob, y2glob, z2glob, &
                          x2glob_f, y2glob_f, z2glob_f, ix_first, ix_last, jy_first, &
                          jy_last, kz_first, kz_last, ix_first_f, ix_last_f, jy_first_f, &
                          jy_last_f, kz_first_f, kz_last_f, resx, resy, resz, resx_f, &
                          resy_f, resz_f, blockglobids, blockglobids_f, blockranks, &
                          blockranks_f)



    CALL check_common_faces(newsubdomain%blocks)
    CALL check_common_edges(newsubdomain%blocks)

    CALL average_face_areas(newsubdomain%blocks)

    CALL newsubdomain%set_cell_type()
    CALL newsubdomain%set_face_type()

    CALL newsubdomain%init_boundcomm_interface(blockranks)
    CALL make_operators(newsubdomain, btype, isblackijkzero_lev(n)%data)

!    newsubdomain%mat => newsubdomain%lapl3d
    subdomain_fine => subdomain_fine%next_coarse

  END DO


END SUBROUTINE intialize_coarse_subdomains


SUBROUTINE init_newblocks(newsubdomain, subdomain_fine, n, nblocks_total, &
                          x2glob, y2glob, z2glob, &
                          x2glob_f, y2glob_f, z2glob_f, ix_first, ix_last, jy_first, &
                          jy_last, kz_first, kz_last, ix_first_f, ix_last_f, jy_first_f, &
                          jy_last_f, kz_first_f, kz_last_f, resx, resy, resz, resx_f, &
                          resy_f, resz_f, blockglobids, blockglobids_f, blockranks, &
                          blockranks_f)

  IMPLICIT NONE

  TYPE(SubDomain), POINTER, INTENT(inout) :: newsubdomain
  TYPE(SubDomain), POINTER, INTENT(in) :: subdomain_fine
  INTEGER, INTENT(in) :: n, nblocks_total

  REAL(Realkind), POINTER, DIMENSION(:), INTENT(in) :: x2glob, y2glob, z2glob, &
                                                       x2glob_f, y2glob_f, z2glob_f
  INTEGER, POINTER, DIMENSION(:), INTENT(in) :: ix_first, ix_last, jy_first, &
                                                jy_last, kz_first, kz_last, &
                                                ix_first_f, ix_last_f, jy_first_f, &
                                                jy_last_f, kz_first_f, kz_last_f, &
                                                resx, resy, resz, resx_f, resy_f, resz_f, &
                                                blockglobids, blockglobids_f, blockranks, &
                                                blockranks_f

  INTEGER :: m, i, j, k, kz, jy, ix, nz, ny, nx, iblock, jblock_fine, ind, ii, jj, kk

  INTEGER :: ix_first_c, ix_last_c, jy_first_c, jy_last_c, kz_first_c, kz_last_c, &
             ix_st_tmp, ix_end_tmp, jy_st_tmp, jy_end_tmp, kz_st_tmp, kz_end_tmp
  INTEGER :: o12(7), o21(7)
  INTEGER :: MPIErr, tag_send, tag_recv, data_size, buff_size, ibuff

  LOGICAL :: is_ghost(nblocks_total)
  LOGICAL :: overlaps

  INTEGER :: ix_first_tmp(nblocks_total + SIZE(subdomain_fine%blockids_comp, 1))
  INTEGER :: ix_last_tmp(nblocks_total + SIZE(subdomain_fine%blockids_comp, 1))
  INTEGER :: jy_first_tmp(nblocks_total + SIZE(subdomain_fine%blockids_comp, 1))
  INTEGER :: jy_last_tmp(nblocks_total + SIZE(subdomain_fine%blockids_comp, 1))
  INTEGER :: kz_first_tmp(nblocks_total + SIZE(subdomain_fine%blockids_comp, 1))
  INTEGER :: kz_last_tmp(nblocks_total + SIZE(subdomain_fine%blockids_comp, 1))

  INTEGER :: sendbufint(NumProcs), recvbufint(NumProcs), &
             reqs_send(NumProcs), reqs_recv(NumProcs), status(MPI_STATUS_SIZE)
  TYPE(IArray) :: jsblocks_from(NumProcs), jsblocks_to(NumProcs), &
                  reqs_send_arr(NumProcs), reqs_recv_arr(NumProcs)
  TYPE(IArray2d) :: sendbufint2d(NumProcs), blockbnds_to(NumProcs)
 
  TYPE(RArray) :: receivebuf1d(NumProcs)


  REAL(Realkind) :: x_st, x_end, y_st, y_end, z_st, z_end

  TYPE(Regular_Grid) :: gridding_grid
  INTEGER, POINTER :: block_ispartof(:,:,:) => NULL()
  INTEGER, POINTER :: block_npartof(:) => NULL()
  INTEGER, POINTER :: neighbor_ids_unique(:) => NULL()
  INTEGER, POINTER :: neighbors(:) => NULL()
  INTEGER, POINTER :: o12_ptr(:) => NULL()

  REAL(Realkind), POINTER :: x2_fine_tmp(:) => NULL(), &
                             y2_fine_tmp(:) => NULL(), &
                             z2_fine_tmp(:) => NULL()

  TYPE(List) :: neighbor_ids, neighbor_ids_transform, list_tmp
  TYPE(List), POINTER :: sublist => NULL(), sublist2 => NULL()
  TYPE(Element), POINTER :: elm => NULL()
  TYPE(gridding_cell), POINTER :: cell => NULL()

  TYPE(Block), POINTER :: neighbor_blocks_unique(:), &
                          thisblock

  TYPE(ListArray) :: jsblocks_from_la(NumProcs), &
                     receive_fineblockparts_la(NumProcs), &
                     blockbnds_from_la(NumProcs)

  ix_first_tmp(1:nblocks_total) = ix_first
  ix_last_tmp(1:nblocks_total) = ix_last
  jy_first_tmp(1:nblocks_total) = jy_first
  jy_last_tmp(1:nblocks_total) = jy_last
  kz_first_tmp(1:nblocks_total) = kz_first
  kz_last_tmp(1:nblocks_total) = kz_last

  j = 1

  DO i = 1, SIZE(subdomain_fine%blockids_compghst, 1)

    IF (subdomain_fine%blockiscomp(i) .EQV. .TRUE.) THEN
      x_st = subdomain_fine%blocks(i)%x2(1)
      x_end = subdomain_fine%blocks(i)%x2(UBOUND(subdomain_fine%blocks(i)%x2, 1))
      y_st = subdomain_fine%blocks(i)%y2(1)
      y_end = subdomain_fine%blocks(i)%y2(UBOUND(subdomain_fine%blocks(i)%y2, 1))
      z_st = subdomain_fine%blocks(i)%z2(1)
      z_end = subdomain_fine%blocks(i)%z2(UBOUND(subdomain_fine%blocks(i)%z2, 1))

      ix_first_tmp(nblocks_total + j) = MAX(FINDLOC(x2glob .GE. x_st, .TRUE., DIM=1) - 1, LBOUND(x2glob, 1))
      ix_last_tmp(nblocks_total + j) = MIN(FINDLOC(x2glob .LE. x_end, .TRUE., DIM=1, BACK=.TRUE.), UBOUND(x2glob, 1) - 1)
      jy_first_tmp(nblocks_total + j) = MAX(FINDLOC(y2glob .GE. y_st, .TRUE., DIM=1) - 1, LBOUND(y2glob, 1))
      jy_last_tmp(nblocks_total + j) = MIN(FINDLOC(y2glob .LE. y_end, .TRUE., DIM=1, BACK=.TRUE.), UBOUND(y2glob, 1) - 1)
      kz_first_tmp(nblocks_total + j) = MAX(FINDLOC(z2glob .GE. z_st, .TRUE., DIM=1) - 1, LBOUND(z2glob, 1))
      kz_last_tmp(nblocks_total + j) = MIN(FINDLOC(z2glob .LE. z_end, .TRUE., DIM=1, BACK=.TRUE.), UBOUND(z2glob, 1) - 1)
      j = j + 1
    END IF
  END DO  

  CALL  grid_blocks(gridding_grid, block_ispartof, block_npartof, x2glob, y2glob, z2glob, &
                     ix_first_tmp, ix_last_tmp, jy_first_tmp, jy_last_tmp, kz_first_tmp, kz_last_tmp)

  DO i = nblocks_total + 1, SIZE(ix_first_tmp, 1)
    sublist => neighbor_ids%new_sublist()
    DO j = 1, block_npartof(i)
      kz = block_ispartof(i, j, 1)
      jy = block_ispartof(i, j, 2)
      ix = block_ispartof(i, j, 3)
      cell => gridding_grid%cells(kz, jy, ix)
      DO k = 1, cell%nconts
        IF (cell%contns(k) .LE. SIZE(ix_first, 1)) CALL list_tmp%append(cell%contns(k))
      END DO
    END DO
    CALL list_tmp%count_deep()
    CALL list_tmp%unique(neighbor_ids_unique)
    CALL list_tmp%destroy()
    DO j = 1, SIZE(neighbor_ids_unique)
      CALL sublist%append(neighbor_ids_unique(j))
    END DO
  END DO

  CALL neighbor_ids%count_deep()
  CALL neighbor_ids%unique(neighbor_ids_unique)

  DO i = 1, neighbor_ids%len
    sublist => neighbor_ids_transform%new_sublist()
    elm => neighbor_ids%get(i)
    sublist2 => elm%sublist
    DO j = 1, sublist2%len
      elm => sublist2%get(j)
      m = FINDLOC(neighbor_ids_unique, elm%ivalue, DIM=1)
      CALL sublist%append(m)
    END DO
  END DO

  DO i = 1, neighbor_ids_transform%len
    elm => neighbor_ids_transform%get(i)
    sublist => elm%sublist
    DO j = 1, sublist%len
      elm => sublist%get(j)
    END DO
  END DO

  ALLOCATE(neighbor_blocks_unique(SIZE(neighbor_ids_unique)))

  DO i = 1, SIZE(neighbor_ids_unique)

    x_st = x2glob(ix_first(neighbor_ids_unique(i)))
    x_end = x2glob(ix_last(neighbor_ids_unique(i)) + 1)
    j = MIN(resx(neighbor_ids_unique(i)) + 1, ngrids)
    ix_first_c = FINDLOC(x2glob_lev(j)%data, x_st, DIM=1)
    ix_last_c = FINDLOC(x2glob_lev(j)%data, x_end, DIM=1)
    ALLOCATE(neighbor_blocks_unique(i)%x2(ix_last_c - ix_first_c + 1))
    neighbor_blocks_unique(i)%x2(:) = x2glob_lev(j)%data(ix_first_c:ix_last_c)

    y_st = y2glob(jy_first(neighbor_ids_unique(i)))
    y_end = y2glob(jy_last(neighbor_ids_unique(i)) + 1)
    j = MIN(resy(neighbor_ids_unique(i)) + 1, ngrids)
    jy_first_c = FINDLOC(y2glob_lev(j)%data, y_st, DIM=1)
    jy_last_c = FINDLOC(y2glob_lev(j)%data, y_end, DIM=1)
    ALLOCATE(neighbor_blocks_unique(i)%y2(jy_last_c - jy_first_c + 1))
    neighbor_blocks_unique(i)%y2(:) = y2glob_lev(j)%data(jy_first_c:jy_last_c)

    z_st = z2glob(kz_first(neighbor_ids_unique(i)))
    z_end = z2glob(kz_last(neighbor_ids_unique(i)) + 1)
    j = MIN(resz(neighbor_ids_unique(i)) + 1, ngrids)
    kz_first_c = FINDLOC(z2glob_lev(j)%data, z_st, DIM=1)
    kz_last_c = FINDLOC(z2glob_lev(j)%data, z_end, DIM=1)
    ALLOCATE(neighbor_blocks_unique(i)%z2(kz_last_c - kz_first_c + 1))
    neighbor_blocks_unique(i)%z2(:) = z2glob_lev(j)%data(kz_first_c:kz_last_c)

    ALLOCATE(neighbor_blocks_unique(i)%x(ix_last_c - ix_first_c))
    neighbor_blocks_unique(i)%x(:) = 0.5 * (neighbor_blocks_unique(i)%x2(1:UBOUND(neighbor_blocks_unique(i)%x2, 1) - 1) + &
                                            neighbor_blocks_unique(i)%x2(2:UBOUND(neighbor_blocks_unique(i)%x2, 1)))
    ALLOCATE(neighbor_blocks_unique(i)%y(jy_last_c - jy_first_c))
    neighbor_blocks_unique(i)%y(:) = 0.5 * (neighbor_blocks_unique(i)%y2(1:UBOUND(neighbor_blocks_unique(i)%y2, 1) - 1) + &
                                            neighbor_blocks_unique(i)%y2(2:UBOUND(neighbor_blocks_unique(i)%y2, 1)))
    ALLOCATE(neighbor_blocks_unique(i)%z(kz_last_c - kz_first_c))
    neighbor_blocks_unique(i)%z(:) = 0.5 * (neighbor_blocks_unique(i)%z2(1:UBOUND(neighbor_blocks_unique(i)%z2, 1) - 1) + &
                                            neighbor_blocks_unique(i)%z2(2:UBOUND(neighbor_blocks_unique(i)%z2, 1)))

    nx = SIZE(neighbor_blocks_unique(i)%x)
    ny = SIZE(neighbor_blocks_unique(i)%y)
    nz = SIZE(neighbor_blocks_unique(i)%z)

    neighbor_blocks_unique(i)%blockid = neighbor_ids_unique(i)
    neighbor_blocks_unique(i)%fld_shape(1) = nz
    neighbor_blocks_unique(i)%fld_shape(2) = ny
    neighbor_blocks_unique(i)%fld_shape(3) = nx
    neighbor_blocks_unique(i)%ncells = nz * ny * nx

  END DO

  CALL check_common_faces(neighbor_blocks_unique)


  newsubdomain%myrank = subdomain_fine%myrank
  IF (SIZE(newsubdomain%blockids_comp, 1) .GE. 1) THEN

    is_ghost(:) = .FALSE.
    DO i = 1, SIZE(newsubdomain%blockids_comp, 1)
      is_ghost(newsubdomain%blockids_comp(i)) = .TRUE.
      neighbors => neighborblocks(gridding_grid, block_ispartof, block_npartof, newsubdomain%blockids_comp(i), &
                                  ix_first_tmp, ix_last_tmp, jy_first_tmp, jy_last_tmp, kz_first_tmp, kz_last_tmp)

      DO j = 1, SIZE(neighbors)
        IF (neighbors(j) .LE. nblocks_total) is_ghost(neighbors(j)) = .TRUE.
      END DO

      DEALLOCATE(neighbors)
      neighbors => NULL()
    END DO

    newsubdomain%nblocks = COUNT(is_ghost, 1)
    ALLOCATE(newsubdomain%blockids_compghst(newsubdomain%nblocks))

    j = 1
    DO i = 1, nblocks_total
      IF (is_ghost(i) .EQV. .TRUE.) THEN
        newsubdomain%blockids_compghst(j) = i
        j = j + 1
      END IF
    END DO

    DO i = 1, SIZE(newsubdomain%blockids_comp, 1)
      is_ghost(newsubdomain%blockids_comp(i)) = .FALSE.
    END DO

    ALLOCATE(newsubdomain%blockiscomp(newsubdomain%nblocks))
    newsubdomain%blockiscomp = .FALSE.
    DO i = 1, SIZE(newsubdomain%blockids_compghst)
      j = newsubdomain%blockids_compghst(i)
      IF (is_ghost(j) .EQV. .FALSE.) newsubdomain%blockiscomp(i) = .TRUE.
    END DO
   
    ALLOCATE(newsubdomain%btype(newsubdomain%nblocks))
    DO i = 1, newsubdomain%nblocks
      newsubdomain%btype(i) = btype(blockglobids(newsubdomain%blockids_compghst(i)))
    END DO

    newsubdomain%ncells_tmp = 0
    newsubdomain%nfaces_tmp = 0
    newsubdomain%ncells = 0
    newsubdomain%nfaces = 0

    ALLOCATE(newsubdomain%blocks(newsubdomain%nblocks))

    DO i = 1, newsubdomain%nblocks

      iblock = newsubdomain%blockids_compghst(i)

      x_st = x2glob(ix_first(iblock))
      x_end = x2glob(ix_last(iblock) + 1)
      j = MIN(resx(iblock) + 1, ngrids)
      ix_first_c = FINDLOC(x2glob_lev(j)%data, x_st, DIM=1)
      ix_last_c = FINDLOC(x2glob_lev(j)%data, x_end, DIM=1)
      ALLOCATE(newsubdomain%blocks(i)%x2(ix_last_c - ix_first_c + 1))
      newsubdomain%blocks(i)%x2(:) = x2glob_lev(j)%data(ix_first_c:ix_last_c)

      y_st = y2glob(jy_first(iblock))
      y_end = y2glob(jy_last(iblock) + 1)
      j = MIN(resy(iblock) + 1, ngrids)
      jy_first_c = FINDLOC(y2glob_lev(j)%data, y_st, DIM=1)
      jy_last_c = FINDLOC(y2glob_lev(j)%data, y_end, DIM=1)
      ALLOCATE(newsubdomain%blocks(i)%y2(jy_last_c - jy_first_c + 1))
      newsubdomain%blocks(i)%y2(:) = y2glob_lev(j)%data(jy_first_c:jy_last_c)

      z_st = z2glob(kz_first(iblock))
      z_end = z2glob(kz_last(iblock) + 1)
      j = MIN(resz(iblock) + 1, ngrids)
      kz_first_c = FINDLOC(z2glob_lev(j)%data, z_st, DIM=1)
      kz_last_c = FINDLOC(z2glob_lev(j)%data, z_end, DIM=1)
      ALLOCATE(newsubdomain%blocks(i)%z2(kz_last_c - kz_first_c + 1))
      newsubdomain%blocks(i)%z2(:) = z2glob_lev(j)%data(kz_first_c:kz_last_c)
  
      ALLOCATE(newsubdomain%blocks(i)%x(ix_last_c - ix_first_c))
      newsubdomain%blocks(i)%x(:) = 0.5 * (newsubdomain%blocks(i)%x2(1:UBOUND(newsubdomain%blocks(i)%x2, 1) - 1) + &
                                           newsubdomain%blocks(i)%x2(2:UBOUND(newsubdomain%blocks(i)%x2, 1)))
      ALLOCATE(newsubdomain%blocks(i)%y(jy_last_c - jy_first_c))
      newsubdomain%blocks(i)%y(:) = 0.5 * (newsubdomain%blocks(i)%y2(1:UBOUND(newsubdomain%blocks(i)%y2, 1) - 1) + &
                                           newsubdomain%blocks(i)%y2(2:UBOUND(newsubdomain%blocks(i)%y2, 1)))
      ALLOCATE(newsubdomain%blocks(i)%z(kz_last_c - kz_first_c))
      newsubdomain%blocks(i)%z(:) = 0.5 * (newsubdomain%blocks(i)%z2(1:UBOUND(newsubdomain%blocks(i)%z2, 1) - 1) + &
                                           newsubdomain%blocks(i)%z2(2:UBOUND(newsubdomain%blocks(i)%z2, 1)))

      newsubdomain%blocks(i)%blockid = iblock

      nx = SIZE(newsubdomain%blocks(i)%x)
      ny = SIZE(newsubdomain%blocks(i)%y)
      nz = SIZE(newsubdomain%blocks(i)%z)

      newsubdomain%blocks(i)%fld_shape(1) = nz
      newsubdomain%blocks(i)%fld_shape(2) = ny
      newsubdomain%blocks(i)%fld_shape(3) = nx
      newsubdomain%blocks(i)%ncells = nz * ny * nx

    END DO

    DO i = 1, newsubdomain%nblocks
      iblock = newsubdomain%blockids_compghst(i)
      thisblock => newsubdomain%blocks(i)

      newsubdomain%blocks(i)%myrank = blockranks(iblock)

      nz = thisblock%fld_shape(1)
      ny = thisblock%fld_shape(2)
      nx = thisblock%fld_shape(3)

      newsubdomain%ncells_tmp = newsubdomain%ncells_tmp + nz * ny * nx
      newsubdomain%nfaces_tmp = newsubdomain%nfaces_tmp + (nz + 1) * ny * nx + nz * (ny + 1) * nx + nz * ny * (nx + 1)
      newsubdomain%blocks(i)%ncells = nz * ny * nx
      newsubdomain%blocks(i)%nfaces = (nz + 1) * ny * nx + nz * (ny + 1) * nx + nz * ny * (nx + 1)

      IF (newsubdomain%blockiscomp(i) .EQV. .TRUE.) THEN
        newsubdomain%ncells = newsubdomain%ncells + newsubdomain%blocks(i)%ncells
        newsubdomain%nfaces = newsubdomain%nfaces + newsubdomain%blocks(i)%nfaces
      END IF

      jblock_fine = FINDLOC(blockglobids_f, blockglobids(iblock), DIM=1)

      IF (ALL((/resx_f(jblock_fine) .EQ. resx(iblock), resy_f(jblock_fine) .EQ. resy(iblock), &
              resz_f(jblock_fine) .EQ. resz(iblock)/))) THEN

        j = FINDLOC(subdomain_fine%blockids_compghst, jblock_fine, DIM=1)

        thisblock%volseff_tmp => subdomain_fine%blocks(j)%volseff
        thisblock%arseffx_tmp => subdomain_fine%blocks(j)%arseffx      
        thisblock%arseffy_tmp => subdomain_fine%blocks(j)%arseffy
        thisblock%arseffz_tmp => subdomain_fine%blocks(j)%arseffz
        thisblock%vols_tmp => subdomain_fine%blocks(j)%vols
        thisblock%arsx_tmp => subdomain_fine%blocks(j)%arsx
        thisblock%arsy_tmp => subdomain_fine%blocks(j)%arsy
        thisblock%arsz_tmp => subdomain_fine%blocks(j)%arsz

        thisblock%x2 => subdomain_fine%blocks(j)%x2
        thisblock%y2 => subdomain_fine%blocks(j)%y2
        thisblock%z2 => subdomain_fine%blocks(j)%z2

        thisblock%x2_fine => subdomain_fine%blocks(j)%x2
        thisblock%y2_fine => subdomain_fine%blocks(j)%y2
        thisblock%z2_fine => subdomain_fine%blocks(j)%z2

      ELSE

        x_st = x2glob(ix_first(iblock))
        x_end = x2glob(ix_last(iblock) + 1)
        ix_st_tmp = FINDLOC(x2glob_f, x_st, DIM=1)
        ix_end_tmp = FINDLOC(x2glob_f, x_end, DIM=1)
        nx = ix_end_tmp - ix_st_tmp
        thisblock%x2_fine => x2glob_f(ix_st_tmp:ix_end_tmp)

        y_st = y2glob(jy_first(iblock))
        y_end = y2glob(jy_last(iblock) + 1)
        jy_st_tmp = FINDLOC(y2glob_f, y_st, DIM=1)
        jy_end_tmp = FINDLOC(y2glob_f, y_end, DIM=1)
        ny = jy_end_tmp - jy_st_tmp
        thisblock%y2_fine => y2glob_f(jy_st_tmp:jy_end_tmp)

        z_st = z2glob(kz_first(iblock))
        z_end = z2glob(kz_last(iblock) + 1)
        kz_st_tmp = FINDLOC(z2glob_f, z_st, DIM=1)
        kz_end_tmp = FINDLOC(z2glob_f, z_end, DIM=1)
        nz = kz_end_tmp - kz_st_tmp
        thisblock%z2_fine => z2glob_f(kz_st_tmp:kz_end_tmp)

        ALLOCATE(thisblock%volseff_tmp(nz, ny, nx))
        ALLOCATE(thisblock%arseffx_tmp(nz, ny, nx + 1))
        ALLOCATE(thisblock%arseffy_tmp(nz, ny + 1, nx))
        ALLOCATE(thisblock%arseffz_tmp(nz + 1, ny, nx))

        thisblock%volseff_tmp(:,:,:) = 0.0
        thisblock%arseffx_tmp(:,:,:) = 0.0
        thisblock%arseffy_tmp(:,:,:) = 0.0

        DO j = 1, SIZE(blockranks_f, 1)
          jblock_fine = blockranks_f(j)
          x_st = x2glob_f(ix_first_f(j))
          x_end = x2glob_f(ix_last_f(j) + 1)
          m = MIN(resx_f(j) + 1, ngrids)
          ix_st_tmp = FINDLOC(x2glob_lev(m)%data, x_st, DIM=1)
          ix_end_tmp = FINDLOC(x2glob_lev(m)%data, x_end, DIM=1)
          x2_fine_tmp => x2glob_lev(m)%data(ix_st_tmp:ix_end_tmp)
          CALL determine_overlap_1d(thisblock%x2_fine, x2_fine_tmp, o12(1:2), o21(1:2))

          y_st = y2glob_f(jy_first_f(j))
          y_end = y2glob_f(jy_last_f(j) + 1)
          m = MIN(resy_f(j) + 1, ngrids)
          jy_st_tmp = FINDLOC(y2glob_lev(m)%data, y_st, DIM=1)
          jy_end_tmp = FINDLOC(y2glob_lev(m)%data, y_end, DIM=1)
          y2_fine_tmp => y2glob_lev(m)%data(jy_st_tmp:jy_end_tmp)
          CALL determine_overlap_1d(thisblock%y2_fine, y2_fine_tmp, o12(3:4), o21(3:4))

          z_st = z2glob_f(kz_first_f(j))
          z_end = z2glob_f(kz_last_f(j) + 1)
          m = MIN(resz_f(j) + 1, ngrids)
          kz_st_tmp = FINDLOC(z2glob_lev(m)%data, z_st, DIM=1)
          kz_end_tmp = FINDLOC(z2glob_lev(m)%data, z_end, DIM=1)
          z2_fine_tmp => z2glob_lev(m)%data(kz_st_tmp:kz_end_tmp)
          CALL determine_overlap_1d(thisblock%z2_fine, z2_fine_tmp, o12(5:6), o21(5:6))

          o12(7) = i
          o21(7) = i

          overlaps = ALL(o12(1:6) .GE. 1)

          IF (overlaps .EQV. .TRUE.) THEN

            IF (ANY(j .EQ. subdomain_fine%blockids_compghst)) THEN
              m = FINDLOC(subdomain_fine%blockids_compghst, j, DIM=1)
              thisblock%volseff_tmp(o12(5):o12(6), o12(3):o12(4), o12(1):o12(2)) = &
                subdomain_fine%blocks(m)%volseff(o21(5):o21(6), o21(3):o21(4), o21(1):o21(2))
              thisblock%arseffx_tmp(o12(5):o12(6), o12(3):o12(4), o12(1):o12(2) + 1) = &
                subdomain_fine%blocks(m)%arseffx(o21(5):o21(6), o21(3):o21(4), o21(1):o21(2) + 1)
              thisblock%arseffy_tmp(o12(5):o12(6), o12(3):o12(4) + 1, o12(1):o12(2)) = &
                subdomain_fine%blocks(m)%arseffy(o21(5):o21(6), o21(3):o21(4) + 1, o21(1):o21(2))
              thisblock%arseffz_tmp(o12(5):o12(6) + 1, o12(3):o12(4), o12(1):o12(2)) = &
                subdomain_fine%blocks(m)%arseffz(o21(5):o21(6) + 1, o21(3):o21(4), o21(1):o21(2))

            ELSE
              m = blockranks_f(j) + 1
              CALL receive_fineblockparts_la(m)%lst%append(o12, copy=.True.)
              CALL jsblocks_from_la(m)%lst%append(j)
              CALL blockbnds_from_la(m)%lst%append(o21, copy=.True.)

            END IF
          END IF
        END DO
      END IF
    END DO
  ELSE
    newsubdomain%nblocks = 0
    newsubdomain%ncells = 0
    newsubdomain%nfaces = 0
    newsubdomain%ncells_tmp = 0
    newsubdomain%nfaces_tmp = 0
    subdomain_fine%ncoarse = 0
  END IF

  DO i = 1, NumProcs
    CALL jsblocks_from_la(i)%lst%toiarray(jsblocks_from(i)%data)
  END DO

  DO i = 1, NumProcs
    sendbufint(i) = SIZE(jsblocks_from(i)%data, 1)
  END DO

  CALL MPI_ALLTOALL(sendbufint, 1, MPI_INT, recvbufint, 1, MPI_INT, MPI_COMM_WORLD, MPIErr)

  DO i = 1, NumProcs
    tag_send = newsubdomain%myrank * NumProcs + i - 1
    tag_recv = (i - 1) * NumProcs + newsubdomain%myrank

    ALLOCATE(jsblocks_to(i)%data(recvbufint(i)))
    ALLOCATE(reqs_send_arr(i)%data(recvbufint(i)))
    ALLOCATE(reqs_recv_arr(i)%data(sendbufint(i)))

    CALL MPI_ISEND(jsblocks_from(i)%data, SIZE(jsblocks_from(i)%data, 1), MPI_INT, i - 1, &
                   tag_send, MPI_COMM_WORLD, reqs_send(i), MPIErr)

    CALL MPI_IRECV(jsblocks_to(i)%data, recvbufint(i), MPI_INT, i - 1, &
                   tag_recv, MPI_COMM_WORLD, reqs_recv(i), MPIErr)
  END DO

  DO i = 1, NumProcs
    CALL MPI_WAIT(reqs_send(i), status, MPIErr)
    CALL MPI_WAIT(reqs_recv(i), status, MPIErr)
  END DO

  DO i = 1, NumProcs
    tag_send = newsubdomain%myrank * NumProcs + i - 1
    tag_recv = (i - 1) * NumProcs + newsubdomain%myrank

    ALLOCATE(sendbufint2d(i)%data(7, SIZE(jsblocks_from(i)%data, 1)))
    ALLOCATE(blockbnds_to(i)%data(7, SIZE(jsblocks_to(i)%data, 1)))

    DO j = 1, blockbnds_from_la(i)%lst%len
      elm => blockbnds_from_la(i)%lst%get(j)
      sendbufint2d(i)%data(:, j) = elm%iarray1d
    END DO

    CALL MPI_ISEND(sendbufint2d(i)%data, SIZE(jsblocks_from(i)%data, 1) * 7, MPI_INT, i - 1, &
                   tag_send, MPI_COMM_WORLD, reqs_send(i), MPIErr)

    CALL MPI_IRECV(blockbnds_to(i)%data, recvbufint(i) * 7, MPI_INT, i - 1, &
                   tag_recv, MPI_COMM_WORLD, reqs_recv(i), MPIErr)

  END DO

  DO i = 1, NumProcs
    CALL MPI_WAIT(reqs_send(i), status, MPIErr)
    CALL MPI_WAIT(reqs_recv(i), status, MPIErr)
  END DO

  DO i  = 1, NumProcs
    buff_size = 0
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)
      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      buff_size = buff_size + (o12_ptr(6) - o12_ptr(5) + 2) * &
                  (o12_ptr(4) - o12_ptr(3) + 2) * (o12_ptr(2) - o12_ptr(1) + 2)
    END DO
    ALLOCATE(receivebuf1d(i)%data(buff_size))
  END DO


  DO i = 1, NumProcs
    tag_send = newsubdomain%myrank * NumProcs + i - 1
    tag_recv = (i - 1) * NumProcs + newsubdomain%myrank

    DO j = 1, SIZE(jsblocks_to(i)%data, 1)
      ind = FINDLOC(subdomain_fine%blockids_compghst, jsblocks_to(i)%data(j), DIM=1)
      o12_ptr => blockbnds_to(i)%data(:, j)
      data_size = (o12_ptr(6) - o12_ptr(5) + 1) * (o12_ptr(4) - o12_ptr(3) + 1) * (o12_ptr(2) - o12_ptr(1) + 1)

      tag_send = newsubdomain%myrank * NumProcs * SIZE(blockranks_f) * SIZE(blockranks) + &
                 (i - 1) * SIZE(blockranks_f) * SIZE(blockranks) + &
                 jsblocks_to(i)%data(j) * SIZE(blockranks) + &
                 o12_ptr(7)

      CALL MPI_ISEND(subdomain_fine%blocks(ind)%volseff(o12_ptr(5):o12_ptr(6), &
                                                        o12_ptr(3):o12_ptr(4), o12_ptr(1):o12_ptr(2)), &
                     data_size, MPI_RealKind, i - 1, tag_send, &
                     MPI_COMM_WORLD, reqs_send_arr(i)%data(j), MPIErr)
    END DO

    ibuff = 1
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)
      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      data_size = (o12_ptr(6) - o12_ptr(5) + 1) * (o12_ptr(4) - o12_ptr(3) + 1) * (o12_ptr(2) - o12_ptr(1) + 1)

      tag_recv = (i - 1) * NumProcs * SIZE(blockranks_f) * SIZE(blockranks) + &
                 newsubdomain%myrank * SIZE(blockranks_f) * SIZE(blockranks) + &
                 jsblocks_from(i)%data(j) * SIZE(blockranks) + &
                 o12_ptr(7)

      CALL MPI_IRECV(receivebuf1d(i)%data(ibuff:ibuff + data_size - 1), &
                     data_size, MPI_RealKind, i - 1, tag_recv, &
                     MPI_COMM_WORLD, reqs_recv_arr(i)%data(j), MPIErr)

      ibuff = ibuff + data_size
    END DO
  END DO

  DO i = 1, NumProcs
    DO j = 1, SIZE(reqs_send_arr(i)%data, 1)
      CALL MPI_WAIT(reqs_send_arr(i)%data(j), status, MPIErr)
    END DO
    DO j = 1, SIZE(reqs_recv_arr(i)%data, 1)
      CALL MPI_WAIT(reqs_recv_arr(i)%data(j), status, MPIErr)
    END DO
  END DO

  DO i  = 1, NumProcs
    ibuff = 1
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)
      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      DO ii = o12_ptr(1), o12_ptr(2)
        DO jj = o12_ptr(3), o12_ptr(4)
          DO kk = o12_ptr(5), o12_ptr(6)
            newsubdomain%blocks(o12_ptr(7))%volseff_tmp(kk, jj, ii) = receivebuf1d(i)%data(ibuff)
            ibuff = ibuff + 1     
          END DO
        END DO
      END DO
    END DO
  END DO

  DO i = 1, NumProcs

    DO j = 1, SIZE(jsblocks_to(i)%data, 1)
      ind = FINDLOC(subdomain_fine%blockids_compghst, jsblocks_to(i)%data(j), DIM=1)
      o12_ptr => blockbnds_to(i)%data(:, j)
      data_size = (o12_ptr(6) - o12_ptr(5) + 1) * (o12_ptr(4) - o12_ptr(3) + 1) * (o12_ptr(2) - o12_ptr(1) + 2)

      tag_send = newsubdomain%myrank * NumProcs * SIZE(blockranks_f) * SIZE(blockranks) + &
                 (i - 1) * SIZE(blockranks_f) * SIZE(blockranks) + &
                 jsblocks_to(i)%data(j) * SIZE(blockranks) + &
                 o12_ptr(7)

      CALL MPI_ISEND(subdomain_fine%blocks(ind)%arseffx(o12_ptr(5):o12_ptr(6), &
                                                        o12_ptr(3):o12_ptr(4), o12_ptr(1):o12_ptr(2) + 1), &
                     data_size, MPI_RealKind, i - 1, tag_send, &
                     MPI_COMM_WORLD, reqs_send_arr(i)%data(j), MPIErr)
    END DO

    ibuff = 1
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)
      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      data_size = (o12_ptr(6) - o12_ptr(5) + 1) * (o12_ptr(4) - o12_ptr(3) + 1) * (o12_ptr(2) - o12_ptr(1) + 2)

      tag_recv = (i - 1) * NumProcs * SIZE(blockranks_f) * SIZE(blockranks) + &
                 newsubdomain%myrank * SIZE(blockranks_f) * SIZE(blockranks) + &
                 jsblocks_from(i)%data(j) * SIZE(blockranks) + &
                 o12_ptr(7)

      CALL MPI_IRECV(receivebuf1d(i)%data(ibuff:ibuff + data_size - 1), &
                     data_size, MPI_RealKind, i - 1, tag_recv, &
                     MPI_COMM_WORLD, reqs_recv_arr(i)%data(j), MPIErr)

      ibuff = ibuff + data_size

    END DO
  END DO

  DO i = 1, NumProcs
    DO j = 1, SIZE(reqs_send_arr(i)%data, 1)
      CALL MPI_WAIT(reqs_send_arr(i)%data(j), status, MPIErr)
    END DO
    DO j = 1, SIZE(reqs_recv_arr(i)%data, 1)
      CALL MPI_WAIT(reqs_recv_arr(i)%data(j), status, MPIErr)
    END DO
  END DO

  DO i  = 1, NumProcs
    ibuff = 1
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)
      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      DO ii = o12_ptr(1), o12_ptr(2) + 1
        DO jj = o12_ptr(3), o12_ptr(4)
          DO kk = o12_ptr(5), o12_ptr(6)
            newsubdomain%blocks(o12_ptr(7))%arseffx_tmp(kk, jj, ii) = receivebuf1d(i)%data(ibuff)
            ibuff = ibuff + 1
          END DO
        END DO
      END DO
    END DO
  END DO 

  DO i = 1, NumProcs
    tag_send = newsubdomain%myrank * NumProcs + i - 1
    tag_recv = (i - 1) * NumProcs + newsubdomain%myrank
    DO j = 1, SIZE(jsblocks_to(i)%data, 1)
      ind = FINDLOC(subdomain_fine%blockids_compghst, jsblocks_to(i)%data(j), DIM=1)
      o12_ptr => blockbnds_to(i)%data(:, j)
      data_size = (o12_ptr(6) - o12_ptr(5) + 1) * (o12_ptr(4) - o12_ptr(3) + 2) * (o12_ptr(2) - o12_ptr(1) + 1)

      tag_send = newsubdomain%myrank * NumProcs * SIZE(blockranks_f) * SIZE(blockranks) + &
                 (i - 1) * SIZE(blockranks_f) * SIZE(blockranks) + &
                 jsblocks_to(i)%data(j) * SIZE(blockranks) + &
                 o12_ptr(7)

      CALL MPI_ISEND(subdomain_fine%blocks(ind)%arseffy(o12_ptr(5):o12_ptr(6), &
                                                        o12_ptr(3):o12_ptr(4) + 1, o12_ptr(1):o12_ptr(2)), &
                     data_size, MPI_RealKind, i - 1, tag_send, &
                     MPI_COMM_WORLD, reqs_send_arr(i)%data(j), MPIErr)
    END DO

    ibuff = 1
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)

      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      data_size = (o12_ptr(6) - o12_ptr(5) + 1) * (o12_ptr(4) - o12_ptr(3) + 2) * (o12_ptr(2) - o12_ptr(1) + 1)

      tag_recv = (i - 1) * NumProcs * SIZE(blockranks_f) * SIZE(blockranks) + &
                 newsubdomain%myrank * SIZE(blockranks_f) * SIZE(blockranks) + &
                 jsblocks_from(i)%data(j) * SIZE(blockranks) + &
                 o12_ptr(7)

      CALL MPI_IRECV(receivebuf1d(i)%data(ibuff:ibuff + data_size - 1), &
                     data_size, MPI_RealKind, i - 1, tag_recv, &
                     MPI_COMM_WORLD, reqs_recv_arr(i)%data(j), MPIErr)

      ibuff = ibuff + data_size

    END DO
  END DO

  DO i = 1, NumProcs
    DO j = 1, SIZE(reqs_send_arr(i)%data, 1)
      CALL MPI_WAIT(reqs_send_arr(i)%data(j), status, MPIErr)
    END DO
    DO j = 1, SIZE(reqs_recv_arr(i)%data, 1)
      CALL MPI_WAIT(reqs_recv_arr(i)%data(j), status, MPIErr)
    END DO
  END DO


  DO i  = 1, NumProcs
    ibuff = 1
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)
      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      DO ii = o12_ptr(1), o12_ptr(2)
        DO jj = o12_ptr(3), o12_ptr(4) + 1
          DO kk = o12_ptr(5), o12_ptr(6)
            newsubdomain%blocks(o12_ptr(7))%arseffy_tmp(kk, jj, ii) = receivebuf1d(i)%data(ibuff)
            ibuff = ibuff + 1
          END DO
        END DO
      END DO
    END DO
  END DO

  DO i = 1, NumProcs
    tag_send = newsubdomain%myrank * NumProcs + i - 1
    tag_recv = (i - 1) * NumProcs + newsubdomain%myrank
    DO j = 1, SIZE(jsblocks_to(i)%data, 1)
      ind = FINDLOC(subdomain_fine%blockids_compghst, jsblocks_to(i)%data(j), DIM=1)
      o12_ptr => blockbnds_to(i)%data(:, j)
      data_size = (o12_ptr(6) - o12_ptr(5) + 2) * (o12_ptr(4) - o12_ptr(3) + 1) * (o12_ptr(2) - o12_ptr(1) + 1)

      tag_send = newsubdomain%myrank * NumProcs * SIZE(blockranks_f) * SIZE(blockranks) + &
                 (i - 1) * SIZE(blockranks_f) * SIZE(blockranks) + &
                 jsblocks_to(i)%data(j) * SIZE(blockranks) + &
                 o12_ptr(7)

      CALL MPI_ISEND(subdomain_fine%blocks(ind)%arseffz(o12_ptr(5):o12_ptr(6) + 1, &
                                                        o12_ptr(3):o12_ptr(4), o12_ptr(1):o12_ptr(2)), &
                     data_size, MPI_RealKind, i - 1, tag_send, &
                     MPI_COMM_WORLD, reqs_send_arr(i)%data(j), MPIErr)
    END DO

    ibuff = 1
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)
      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      data_size = (o12_ptr(6) - o12_ptr(5) + 2) * (o12_ptr(4) - o12_ptr(3) + 1) * (o12_ptr(2) - o12_ptr(1) + 1)

      tag_recv = (i - 1) * NumProcs * SIZE(blockranks_f) * SIZE(blockranks) + &
                 newsubdomain%myrank * SIZE(blockranks_f) * SIZE(blockranks) + &
                 jsblocks_from(i)%data(j) * SIZE(blockranks) + &
                 o12_ptr(7)


      CALL MPI_IRECV(receivebuf1d(i)%data(ibuff:ibuff + data_size - 1), &
                     data_size, MPI_RealKind, i - 1, tag_recv, &
                     MPI_COMM_WORLD, reqs_recv_arr(i)%data(j), MPIErr)

      ibuff = ibuff + data_size

    END DO
  END DO

  DO i = 1, NumProcs
    DO j = 1, SIZE(reqs_send_arr(i)%data, 1)
      CALL MPI_WAIT(reqs_send_arr(i)%data(j), status, MPIErr)
    END DO
    DO j = 1, SIZE(reqs_recv_arr(i)%data, 1)
      CALL MPI_WAIT(reqs_recv_arr(i)%data(j), status, MPIErr)
    END DO
  END DO

  DO i  = 1, NumProcs
    ibuff = 1
    DO j = 1, SIZE(jsblocks_from(i)%data, 1)
      elm => receive_fineblockparts_la(i)%lst%get(j)
      o12_ptr => elm%iarray1d
      DO ii = o12_ptr(1), o12_ptr(2)
        DO jj = o12_ptr(3), o12_ptr(4)
          DO kk = o12_ptr(5), o12_ptr(6) + 1
            newsubdomain%blocks(o12_ptr(7))%arseffz_tmp(kk, jj, ii) = receivebuf1d(i)%data(ibuff)
            ibuff = ibuff + 1
          END DO
        END DO
      END DO
    END DO
  END DO 

  DO iblock = 1, newsubdomain%nblocks

    thisblock => newsubdomain%blocks(iblock)

    CALL coarsen_cfield(thisblock%volseff, thisblock%volseff_tmp, thisblock%x2_fine, thisblock%y2_fine, &
                        thisblock%z2_fine, thisblock%x2, thisblock%y2, thisblock%z2)

    CALL coarsen_ufield(thisblock%arseffx, thisblock%arseffx_tmp, thisblock%x2_fine, thisblock%y2_fine, &
                        thisblock%z2_fine, thisblock%x2, thisblock%y2, thisblock%z2)

    CALL coarsen_vfield(thisblock%arseffy, thisblock%arseffy_tmp, thisblock%x2_fine, thisblock%y2_fine, &
                        thisblock%z2_fine, thisblock%x2, thisblock%y2, thisblock%z2)

    CALL coarsen_wfield(thisblock%arseffz, thisblock%arseffz_tmp, thisblock%x2_fine, thisblock%y2_fine, &
                        thisblock%z2_fine, thisblock%x2, thisblock%y2, thisblock%z2)

    nz = thisblock%fld_shape(1)
    ny = thisblock%fld_shape(2)
    nx = thisblock%fld_shape(3)

    ALLOCATE(thisblock%vols(nz, ny, nx))
    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz
          thisblock%vols(k, j, i) = (thisblock%x2(i + 1) - thisblock%x2(i)) * &
                                    (thisblock%y2(j + 1) - thisblock%y2(j)) * &
                                    (thisblock%z2(k + 1) - thisblock%z2(k))
        END DO
      END DO
    END DO

    ALLOCATE(thisblock%arsx(nz, ny, nx + 1))
    DO j = 1, ny
      DO k = 1, nz
        thisblock%arsx(k, j, :) = (thisblock%y2(j + 1) - thisblock%y2(j)) * (thisblock%z2(k + 1) - thisblock%z2(k))
      END DO
    END DO

    ALLOCATE(thisblock%arsy(nz, ny + 1, nx))
    DO i = 1, nx
      DO k = 1, nz
        thisblock%arsy(k, :, i) = (thisblock%x2(i + 1) - thisblock%x2(i)) * (thisblock%z2(k + 1) - thisblock%z2(k))
      END DO
    END DO

    ALLOCATE(thisblock%arsz(nz + 1, ny, nx))
    DO i = 1, nx
      DO j = 1, ny
        thisblock%arsz(:, j, i) = (thisblock%x2(i + 1) - thisblock%x2(i)) * (thisblock%y2(j + 1) - thisblock%y2(j))
      END DO
    END DO

  END DO

  CALL init_transfer(newsubdomain, subdomain_fine, neighbor_ids, &
                     neighbor_ids_transform, neighbor_ids_unique, &
                     neighbor_blocks_unique, blockranks, n)

  subdomain_fine%ncoarse = subdomain_fine%Restr_samerank%m

  DO i = 1, SIZE(neighbor_blocks_unique)
    CALL neighbor_blocks_unique(i)%destroy() 
  END DO

  DEALLOCATE(neighbor_blocks_unique)
  neighbor_blocks_unique => NULL()
  DEALLOCATE(neighbor_ids_unique)
  neighbor_ids_unique => NULL()
  DEALLOCATE(block_ispartof)
  block_ispartof => NULL()
  DEALLOCATE(block_npartof)
  block_npartof => NULL()
  o12_ptr => NULL()
  x2_fine_tmp => NULL()
  y2_fine_tmp => NULL()
  z2_fine_tmp => NULL()


END SUBROUTINE init_newblocks



SUBROUTINE coarsen_cfield(field, field_tmp, x2_tmp, y2_tmp, z2_tmp, x2, y2, z2)

  IMPLICIT NONE
  REAL(Realkind), INTENT(in) :: field_tmp(:,:,:)
  REAL(Realkind), POINTER, INTENT(inout) :: field(:,:,:)
  REAL(Realkind), POINTER, INTENT(in), DIMENSION(:) :: x2_tmp, y2_tmp, z2_tmp, x2, y2, z2

  INTEGER :: nz, ny, nx, i, j, k, ii, jj, kk
  INTEGER :: nz_tmp, ny_tmp, nx_tmp

  nz = SIZE(z2, 1) - 1
  ny = SIZE(y2, 1) - 1
  nx = SIZE(x2, 1) - 1

  nz_tmp = SIZE(z2_tmp, 1) - 1
  ny_tmp = SIZE(y2_tmp, 1) - 1
  nx_tmp = SIZE(x2_tmp, 1) - 1

  ALLOCATE(field(nz, ny, nx))
  field(:,:,:) = 0.0

  DO i = 1, nx_tmp
    ii = between(x2, 0.5 * (x2_tmp(i) + x2_tmp(i + 1)))
    DO j = 1, ny_tmp
      jj = between(y2, 0.5 * (y2_tmp(j) + y2_tmp(j + 1)))
      DO k = 1, nz_tmp
        kk = between(z2, 0.5 * (z2_tmp(k) + z2_tmp(k + 1)))
        field(kk, jj, ii) = field(kk, jj, ii) + field_tmp(k, j, i)
      END DO
    END DO
  END DO

END SUBROUTINE coarsen_cfield


SUBROUTINE coarsen_ufield(field, field_tmp, x2_tmp, y2_tmp, z2_tmp, x2, y2, z2)

  IMPLICIT NONE
  REAL(Realkind), INTENT(in) :: field_tmp(:,:,:)
  REAL(Realkind), POINTER, INTENT(inout) :: field(:,:,:)
  REAL(Realkind), POINTER, INTENT(in), DIMENSION(:) :: x2_tmp, y2_tmp, z2_tmp, x2, y2, z2
  REAL(Realkind) :: weights(3)
  REAL(Realkind) :: fieldtmp(SIZE(z2, 1) - 1, SIZE(y2, 1) - 1, SIZE(x2_tmp, 1))

  INTEGER :: nz, ny, nx, i1, i2, i3, j, k, ii, jj, kk
  INTEGER :: nz_tmp, ny_tmp, nx_tmp

  weights(:) = (/(1.0 - w) / 2.0, w, (1.0 - w) / 2.0/)

  nz = SIZE(z2, 1) - 1
  ny = SIZE(y2, 1) - 1
  nx = SIZE(x2, 1) - 1

  nz_tmp = SIZE(z2_tmp, 1) - 1
  ny_tmp = SIZE(y2_tmp, 1) - 1
  nx_tmp = SIZE(x2_tmp, 1) - 1

  ALLOCATE(field(nz, ny, nx + 1))
  field(:,:,:) = 0.0
  fieldtmp(:,:,:) = 0.0

  DO ii = 1, nx_tmp + 1
    DO j = 1, ny_tmp
      jj = between(y2, 0.5 * (y2_tmp(j) + y2_tmp(j + 1)))
      DO k = 1, nz_tmp
        kk = between(z2, 0.5 * (z2_tmp(k) + z2_tmp(k + 1)))
        fieldtmp(kk, jj, ii) = fieldtmp(kk,jj,ii) + field_tmp(k, j, ii)
      END DO
    END DO
  END DO

  DO ii = 2, nx
    i1 = FINDLOC(x2_tmp, x2(ii - 1), DIM=1)
    i2 = FINDLOC(x2_tmp, x2(ii), DIM=1)
    i3 = FINDLOC(x2_tmp, x2(ii + 1), DIM=1)

    i1 = int((i1 + i2) / 2.0)
    i3 = int((i3 + i2) / 2.0)

    DO j = 1, ny
      DO k = 1, nz
        IF (fieldtmp(k, j, i2) .GT. 1e-20) THEN
          field(k, j, ii) = field(k, j, ii) + fieldtmp(k, j, i2) * weights(2)
          field(k, j, ii) = field(k, j, ii) + fieldtmp(k, j, i1) * weights(1)
          field(k, j, ii) = field(k, j, ii) + fieldtmp(k, j, i3) * weights(3)
        END IF
      END DO
    END DO
  END DO

  i2 = FINDLOC(x2_tmp, x2(1), DIM=1)
  i3 = FINDLOC(x2_tmp, x2(2), DIM=1)
  i3 = int((i2 + i3) / 2.0)

  DO j = 1, ny
    DO k = 1, nz
      IF (fieldtmp(k, j, i2) .GT. 1e-20) THEN
        field(k, j, 1) = field(k, j, 1) + fieldtmp(k, j, i2) * weights(2)
        field(k, j, 1) = field(k, j, 1) + fieldtmp(k, j, i3) * (weights(1) + weights(3))
      END IF 
    END DO
  END DO

  i1 = FINDLOC(x2_tmp, x2(UBOUND(x2, 1) - 1), DIM=1)
  i2 = FINDLOC(x2_tmp, x2(UBOUND(x2, 1)), DIM=1)
  i1 = int((i1 + i2) / 2.0)

  DO j = 1, ny
    DO k = 1, nz

      IF (fieldtmp(k, j, i2) .GT. 1e-20) THEN
        field(k, j, UBOUND(x2, 1)) = field(k, j, UBOUND(x2, 1)) + fieldtmp(k, j, i2) * weights(2)
        field(k, j, UBOUND(x2, 1)) = field(k, j, UBOUND(x2, 1)) + fieldtmp(k, j, i1) * (weights(1) + weights(3))
      END IF
    END DO
  END DO

END SUBROUTINE coarsen_ufield


SUBROUTINE coarsen_vfield(field, field_tmp, x2_tmp, y2_tmp, z2_tmp, x2, y2, z2)
  
  IMPLICIT NONE
  REAL(Realkind), INTENT(in) :: field_tmp(:,:,:)
  REAL(Realkind), POINTER, INTENT(inout) :: field(:,:,:)
  REAL(Realkind), POINTER, INTENT(in), DIMENSION(:) :: x2_tmp, y2_tmp, z2_tmp, x2, y2, z2 
  REAL(Realkind) :: fieldtmp(SIZE(z2, 1) - 1, SIZE(y2_tmp, 1), SIZE(x2, 1) - 1)
  REAL(Realkind) :: weights(3)
  
  INTEGER :: nz, ny, nx, i, j1, j2, j3, k, ii, jj, kk
  INTEGER :: nz_tmp, ny_tmp, nx_tmp
  
  weights(:) = (/(1.0 - w) / 2.0, w, (1.0 - w) / 2.0/)

  nz = SIZE(z2, 1) - 1
  ny = SIZE(y2, 1) - 1
  nx = SIZE(x2, 1) - 1 
  
  nz_tmp = SIZE(z2_tmp, 1) - 1
  ny_tmp = SIZE(y2_tmp, 1) - 1
  nx_tmp = SIZE(x2_tmp, 1) - 1
  
  ALLOCATE(field(nz, ny + 1, nx))
  field(:,:,:) = 0.0 
  fieldtmp(:,:,:) = 0.0

  DO jj = 1, ny_tmp + 1
    DO i = 1, nx_tmp
      ii = between(x2, 0.5 * (x2_tmp(i) + x2_tmp(i + 1)))
      DO k = 1, nz_tmp
        kk = between(z2, 0.5 * (z2_tmp(k) + z2_tmp(k + 1)))
        fieldtmp(kk, jj, ii) = fieldtmp(kk, jj, ii) + field_tmp(k, jj, i)
      END DO
    END DO
  END DO

  DO jj = 2, ny
    j1 = FINDLOC(y2_tmp, y2(jj - 1), DIM=1)
    j2 = FINDLOC(y2_tmp, y2(jj), DIM=1)
    j3 = FINDLOC(y2_tmp, y2(jj + 1), DIM=1)
    j1 = int((j1 + j2) / 2.0)
    j3 = int((j3 + j2) / 2.0)
    DO i = 1, nx
      DO k = 1, nz
        IF (fieldtmp(k, j2, i) .GT. 1e-20) THEN
          field(k, jj, i) = field(k, jj, i) + fieldtmp(k, j2, i) * weights(2)
          field(k, jj, i) = field(k, jj, i) + fieldtmp(k, j1, i) * weights(1)
          field(k, jj, i) = field(k, jj, i) + fieldtmp(k, j3, i) * weights(3)
        END IF
      END DO
    END DO
  END DO

  j2 = FINDLOC(y2_tmp, y2(1), DIM=1)
  j3 = FINDLOC(y2_tmp, y2(2), DIM=1)
  j3 = int((j2 + j3) / 2.0)

  DO i = 1, nx
    DO k = 1, nz
      IF (fieldtmp(k, j2, i) .GT. 1e-20) THEN
        field(k, 1, i) = field(k, 1, i) + fieldtmp(k, j2, i) * weights(2)
        field(k, 1, i) = field(k, 1, i) + fieldtmp(k, j3, i) * (weights(1) + weights(3))
      END IF
    END DO
  END DO

  j1 = FINDLOC(y2_tmp, y2(UBOUND(y2, 1) - 1), DIM=1)
  j2 = FINDLOC(y2_tmp, y2(UBOUND(y2, 1)), DIM=1)
  j1 = int((j1 + j2) / 2.0)

  DO i = 1, nx
    DO k = 1, nz
      IF (fieldtmp(k, j2, i) .GT. 1e-20) THEN
        field(k, UBOUND(y2, 1), i) = field(k, UBOUND(y2, 1), i) + fieldtmp(k, j2, i) * weights(2)
        field(k, UBOUND(y2, 1), i) = field(k, UBOUND(y2, 1), i) + fieldtmp(k, j1, i) * (weights(1) + weights(3))
      END IF
    END DO
  END DO

END SUBROUTINE coarsen_vfield


SUBROUTINE coarsen_wfield(field, field_tmp, x2_tmp, y2_tmp, z2_tmp, x2, y2, z2)
  
  IMPLICIT NONE
  REAL(Realkind), INTENT(in) :: field_tmp(:,:,:)
  REAL(Realkind), POINTER, INTENT(inout) :: field(:,:,:)
  REAL(Realkind), POINTER, INTENT(in), DIMENSION(:) :: x2_tmp, y2_tmp, z2_tmp, x2, y2, z2
  REAL(Realkind) :: fieldtmp(SIZE(z2_tmp, 1), SIZE(y2, 1) - 1, SIZE(x2, 1) - 1)
  REAL(Realkind) :: weights(3)
  
  INTEGER :: nz, ny, nx, i, j, k1, k2, k3, ii, jj, kk
  INTEGER :: nz_tmp, ny_tmp, nx_tmp
  
  weights(:) = (/(1.0 - w) / 2.0, w, (1.0 - w) / 2.0/)

  nz = SIZE(z2, 1) - 1
  ny = SIZE(y2, 1) - 1
  nx = SIZE(x2, 1) - 1
  
  nz_tmp = SIZE(z2_tmp, 1) - 1
  ny_tmp = SIZE(y2_tmp, 1) - 1
  nx_tmp = SIZE(x2_tmp, 1) - 1
  
  ALLOCATE(field(nz + 1, ny, nx))
  field(:,:,:) = 0.0
  fieldtmp(:,:,:) = 0.0
    
  DO kk = 1, nz_tmp + 1
    DO i = 1, nx_tmp
      ii = between(x2, 0.5 * (x2_tmp(i) + x2_tmp(i + 1)))
      DO j = 1, ny_tmp
        jj = between(y2, 0.5 * (y2_tmp(j) + y2_tmp(j + 1)))
        fieldtmp(kk, jj, ii) = fieldtmp(kk, jj, ii) + field_tmp(kk, j, i)
      END DO
    END DO
  END DO

  DO kk = 2, nz
    k2 = FINDLOC(z2_tmp, z2(kk), DIM=1)
    k1 = FINDLOC(z2_tmp, z2(kk - 1), DIM=1)
    k3 = FINDLOC(z2_tmp, z2(kk + 1), DIM=1)
    k1 = int((k1 + k2) / 2.0)
    k3 = int((k3 + k2) / 2.0)
    DO i = 1, nx
      DO j = 1, ny
        IF (fieldtmp(k2, j, i) .GT. 1e-20) THEN
          field(kk, j, i) = field(kk, j, i) + fieldtmp(k2, j, i) * weights(2)
          field(kk, j, i) = field(kk, j, i) + fieldtmp(k1, j, i) * weights(1)
          field(kk, j, i) = field(kk, j, i) + fieldtmp(k3, j, i) * weights(3)
        END IF
      END DO
    END DO
  END DO

  k2 = FINDLOC(z2_tmp, z2(1), DIM=1)
  k3 = FINDLOC(z2_tmp, z2(2), DIM=1)
  k3 = int((k3 + k2) / 2.0)

  DO i = 1, nx
    DO j = 1, ny
      IF (fieldtmp(k2, j, i) .GT. 1e-20) THEN
        field(1, j, i) = field(1, j, i) + fieldtmp(k2, j, i) * weights(2)
        field(1, j, i) = field(1, j, i) + fieldtmp(k3, j, i) * (weights(1) + weights(3))
      END IF
    END DO
  END DO

  k2 = FINDLOC(z2_tmp, z2(UBOUND(z2, 1)), DIM=1)
  k1 = FINDLOC(z2_tmp, z2(UBOUND(z2, 1) - 1), DIM=1)
  k1 = int((k2 + k1) / 2.0)

  DO i = 1, nx
    DO j = 1, ny
      IF (fieldtmp(k2, j, i) .GT. 1e-20) THEN
        field(UBOUND(z2, 1), j, i) = field(UBOUND(z2, 1), j, i) + fieldtmp(k2, j, i) * weights(2)
        field(UBOUND(z2, 1), j, i) = field(UBOUND(z2, 1), j, i) + fieldtmp(k1, j, i) * (weights(1) + weights(3))
      END IF
    END DO
  END DO


END SUBROUTINE coarsen_wfield

SUBROUTINE average_face_areas(blocks)
  IMPLICIT NONE
  
  TYPE(Block), POINTER, INTENT(inout) :: blocks(:)
  INTEGER :: iblock, nblocks
  INTEGER :: i, j, k, nx, ny, nz, bndblock_id
  INTEGER :: i_bnds(2), j_bnds(2), k_bnds(2)
  REAL(Realkind), POINTER :: areas(:, :, :), arsx_bnd(:, :, :)

  IF(.NOT. ASSOCIATED(blocks)) RETURN

  nblocks = SIZE(blocks)

  DO iblock = 1, nblocks
    nz = blocks(iblock)%fld_shape(1)
    ny = blocks(iblock)%fld_shape(2)
    nx = blocks(iblock)%fld_shape(3)
 
    IF (ANY(blocks(iblock)%cface_w)) THEN
      DO j = 1, ny
        DO k = 1, nz
          CALL find_indices_jk(blocks, blocks(iblock)%cface_w, j, j_bnds, k, k_bnds, iblock, bndblock_id)
          IF (bndblock_id .LE. 99999999) THEN
            IF ((k_bnds(2) - k_bnds(1) .EQ. 0) .AND. (j_bnds(2) - j_bnds(1) .EQ. 0)) THEN
              blocks(iblock)%arseffx(k, j, 1) = 0.5 * (blocks(iblock)%arseffx(k, j, 1) + &
                                                blocks(bndblock_id)%arseffx(k, j, UBOUND(blocks(bndblock_id)%arseffx, 3)))
              blocks(bndblock_id)%arseffx(k, j, UBOUND(blocks(bndblock_id)%arseffx, 3)) = blocks(iblock)%arseffx(k, j, 1)
            END IF
          END IF
        END DO
      END DO
    END IF

    IF (ANY(blocks(iblock)%cface_e)) THEN
      i = UBOUND(blocks(iblock)%arseffx, 3)
      DO j = 1, ny
        DO k = 1, nz
          CALL find_indices_jk(blocks, blocks(iblock)%cface_e, j, j_bnds, k, k_bnds, iblock, bndblock_id)
          IF (bndblock_id .LE. 99999999) THEN
            IF ((k_bnds(2) - k_bnds(1) .EQ. 0) .AND. (j_bnds(2) - j_bnds(1) .EQ. 0)) THEN
              blocks(iblock)%arseffx(k, j, i) = 0.5 * (blocks(iblock)%arseffx(k, j, i) + blocks(bndblock_id)%arseffx(k, j, 1))
              blocks(bndblock_id)%arseffx(k, j, 1) = blocks(iblock)%arseffx(k, j, i)
            END IF
          END IF
        END DO
      END DO
    END IF

    IF (ANY(blocks(iblock)%cface_s)) THEN
      DO i = 1, nx
        DO k = 1, nz
          CALL find_indices_ik(blocks, blocks(iblock)%cface_s, i, i_bnds, k, k_bnds, iblock, bndblock_id)
          IF (bndblock_id .LE. 99999999) THEN
            IF ((k_bnds(2) - k_bnds(1) .EQ. 0) .AND. (i_bnds(2) - i_bnds(1) .EQ. 0)) THEN
              blocks(iblock)%arseffy(k, 1, i) = 0.5 * (blocks(iblock)%arseffy(k, 1, i) + &
                             blocks(bndblock_id)%arseffy(k, UBOUND(blocks(bndblock_id)%arseffy, 2), i))
              blocks(bndblock_id)%arseffy(k, UBOUND(blocks(bndblock_id)%arseffy, 2), i) = blocks(iblock)%arseffy(k, 1, i)
            END IF
          END IF
        END DO
      END DO
    END IF

    IF (ANY(blocks(iblock)%cface_n)) THEN
      j = UBOUND(blocks(iblock)%arseffy, 2)
      DO i = 1, nx
        DO k = 1, nz
          CALL find_indices_ik(blocks, blocks(iblock)%cface_n, i, i_bnds, k, k_bnds, iblock, bndblock_id)
          IF (bndblock_id .LE. 99999999) THEN
            IF ((k_bnds(2) - k_bnds(1) .EQ. 0) .AND. (i_bnds(2) - i_bnds(1) .EQ. 0)) THEN
              blocks(iblock)%arseffy(k, j, i) = 0.5 * (blocks(iblock)%arseffy(k, j, i) + blocks(bndblock_id)%arseffy(k, 1, i))
              blocks(bndblock_id)%arseffy(k, 1, i) = blocks(iblock)%arseffy(k, j, i)
            END IF
          END IF
        END DO
      END DO
    END IF

    IF (ANY(blocks(iblock)%cface_b)) THEN
      DO i = 1, nx
        DO j = 1, ny
          CALL find_indices_ij(blocks, blocks(iblock)%cface_b, i, i_bnds, j, j_bnds, iblock, bndblock_id)
          IF (bndblock_id .LE. 99999999) THEN
            IF ((i_bnds(2) - i_bnds(1) .EQ. 0) .AND. (j_bnds(2) - j_bnds(1) .EQ. 0)) THEN
              blocks(iblock)%arseffz(1, j, i) = 0.5 * (blocks(iblock)%arseffz(1, j, i) + &
                                                blocks(bndblock_id)%arseffz(UBOUND(blocks(bndblock_id)%arseffz, 1), j, i))
              blocks(bndblock_id)%arseffz(UBOUND(blocks(bndblock_id)%arseffz, 1), j, i) = blocks(iblock)%arseffz(1, j, i) 
           END IF
          END IF
        END DO
      END DO
    END IF

    IF (ANY(blocks(iblock)%cface_t)) THEN
      k = UBOUND(blocks(iblock)%arseffz, 1)
      DO i = 1, nx
        DO j = 1, ny
          CALL find_indices_ij(blocks, blocks(iblock)%cface_t, i, i_bnds, j, j_bnds, iblock, bndblock_id)
          IF (bndblock_id .LE. 99999999) THEN
            IF ((i_bnds(2) - i_bnds(1) .EQ. 0) .AND. (j_bnds(2) - j_bnds(1) .EQ. 0)) THEN
              blocks(iblock)%arseffz(k, j, i) = 0.5 * (blocks(iblock)%arseffz(k, j, i) + blocks(bndblock_id)%arseffz(1, j, i))
              blocks(bndblock_id)%arseffz(1, j, i) = blocks(iblock)%arseffz(k, j, i)
            END IF
          END IF
        END DO
      END DO
    END IF
  END DO

END SUBROUTINE average_face_areas


SUBROUTINE determine_color_startcell()

  IMPLICIT NONE

  INTEGER :: igrid, iblock, nblocks
  INTEGER :: ix_st, jy_st, kz_st
  INTEGER :: modix2, modiy2, modiz2
  ALLOCATE(isblackijkzero_lev(ngrids))

  DO igrid = 1, ngrids
    nblocks = SIZE(resx_lev(igrid)%data, 1)
    ALLOCATE(isblackijkzero_lev(igrid)%data(nblocks))

    DO iblock = 1, nblocks
      ix_st = ix_st_lev(igrid)%data(iblock)
      jy_st = jy_st_lev(igrid)%data(iblock)
      kz_st = kz_st_lev(igrid)%data(iblock)
      modix2 = MOD(ix_st, 2)
      modiy2 = MOD(jy_st, 2)
      modiz2 = MOD(kz_st, 2)

      isblackijkzero_lev(igrid)%data(iblock) = (modix2 * modiy2 + (1 - modix2) * (1 - modiy2)) * modiz2 + &
                                           ((1 - modix2) * modiy2 + modix2 * (1 - modiy2)) * (1 - modiz2)

    END DO
  END DO 

END SUBROUTINE determine_color_startcell

  
END MODULE CoarseGrids_Mod
