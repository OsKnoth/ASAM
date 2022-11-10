MODULE Transfer_Mod

USE MatSpRowCol_Mod, ONLY: SpRowCol, &
                           SpDeallocate_SpRowCol, SpNullify_SpRowCol, &
                           SpTrans_SpRowCol, SpMm_SpRowCol, SpA_SpRowCol, &
                           SpAVec_SpRowCol, ATxpy_SpRowCol, SpAScalRow_SpRowCol, &
                           SpAScalCol_SpRowCol, Axpy_SpRowCol, sort


USE Index_Mod, ONLY: between, IndexC, IndexC3block, &
                     determine_overlap_1d, determine_overlap_1d_halo
USE Block_Mod, ONLY: Block
USE MPI_Mod
USE Boundary_Mod, ONLY: Bound, copy_data_on_bnd, copy_data_from_bnd
USE Subdom_Mod, ONLY: Subdomain, TransfCommInterface
USE List_Mod, ONLY: List, Element, ListArray, IArray, IArray2d
USE Kind_Mod

IMPLICIT NONE

CONTAINS


SUBROUTINE init_transfer(newsubdomain, subdomain_fine, neighbor_ids, &
                         neighbor_ids_transform, neighbor_ids_unique, &
                         neighbor_blocks_unique, blockranks, ilev)


  TYPE(SubDomain), INTENT(inout) :: newsubdomain
  TYPE(SubDomain), INTENT(inout) :: subdomain_fine
  TYPE(List), INTENT(in) :: neighbor_ids, neighbor_ids_transform
  INTEGER, INTENT(in) :: neighbor_ids_unique(:), &
                         blockranks(:) 
  TYPE(Block), INTENT(in), TARGET :: neighbor_blocks_unique(:)
  INTEGER, INTENT(in) :: ilev


  TYPE(TransfCommInterface), POINTER :: new_restrcomm_if  => NULL()

  TYPE(SpRowCol) :: Prol_fo, Prol_so, Restr, Cut_samerank, Cut_samerank_T, Cut_offrank, Cut_offrank_T

  INTEGER :: i, j, k, ip, n, m, ind, iblock, jblock, icell_st, bind, idata, nranks_comm, nnz

  TYPE(Element), POINTER :: elm => NULL()
  TYPE(List), POINTER :: possible_overlapblock_ids => NULL(), &
                         ind_trans => NULL()
  TYPE(Block), POINTER :: block_coarse => NULL(), &
                          block_fine => NULL()

  INTEGER :: o12(6, neighbor_ids%len_flat), &
             o21(6, neighbor_ids%len_flat)
  INTEGER :: nnzs(neighbor_ids%len_flat), &
             Rowptr_sizes(neighbor_ids%len_flat)
  LOGICAL :: overlaps(neighbor_ids%len_flat)

  TYPE(ListArray) :: jsblocks_to_la(NumProcs), &
                     blockbnds_to_la(NumProcs), &
                     send_transinds(NumProcs), &
                     overlap_blockids_transinds(NumProcs)

  INTEGER :: sendbufint(NumProcs), recvbufint(NumProcs), &
             reqs_send(NumProcs), reqs_recv(NumProcs), status(MPI_STATUS_SIZE)
  INTEGER :: blockrank, MPIErr, tag_send, tag_recv, data_size

  TYPE(IArray) :: jsblocks_from(NumProcs), jsblocks_to(NumProcs), &
                  reqs_send_arr(NumProcs), reqs_recv_arr(NumProcs), &
                  indices_unique_to(NumProcs), block_ptr_to(NumProcs), &
                  indices_unique_from(NumProcs), block_ptr_from(NumProcs)
  TYPE(IArray2d) :: blockbnds_from(NumProcs), blockbnds_to(NumProcs)

  INTEGER :: ncells_fine, ncells_coarse_z, ncells_coarse_zy, ncells_coarse_zyx, &
             prolz_nnz_ub, proly_nnz_ub, prolx_nnz_ub

  TYPE(Bound), POINTER :: current => NULL()
  TYPE(Block), POINTER :: thisblock => NULL()
  INTEGER, POINTER :: o21_ptr(:) => NULL()

  nnzs(:) = 0
  Rowptr_sizes(:) = 0

  ind = 1
  m = 1

  DO n = 1, subdomain_fine%nblocks

    IF (subdomain_fine%blockiscomp(n)) THEN

      block_fine => subdomain_fine%blocks(n)
      iblock = subdomain_fine%blockids_comp(ind)

      elm => neighbor_ids%get(ind)
      possible_overlapblock_ids => elm%sublist

      elm => neighbor_ids_transform%get(ind)
      ind_trans => elm%sublist

      ind = ind + 1

      DO bind = 1, possible_overlapblock_ids%len
        elm => possible_overlapblock_ids%get(bind)
        jblock = elm%ivalue

        elm => ind_trans%get(bind)
        block_coarse => neighbor_blocks_unique(elm%ivalue)

        CALL determine_overlap_1d_halo(block_fine%x2, block_coarse%x2, o12(1:2, m), o21(1:2, m))
        CALL determine_overlap_1d_halo(block_fine%y2, block_coarse%y2, o12(3:4, m), o21(3:4, m))
        CALL determine_overlap_1d_halo(block_fine%z2, block_coarse%z2, o12(5:6, m), o21(5:6, m))

!        CALL determine_overlap_1d(block_fine%x2, block_coarse%x2, o12(1:2, m), o21(1:2, m))
!        CALL determine_overlap_1d(block_fine%y2, block_coarse%y2, o12(3:4, m), o21(3:4, m))
!        CALL determine_overlap_1d(block_fine%z2, block_coarse%z2, o12(5:6, m), o21(5:6, m))

        IF (ALL(o12(:, m) .GT. 0)) THEN
          overlaps(m) = .TRUE.

          nnzs(m) = (o12(2, m) - o12(1, m) + 1) * (o12(4, m) - o12(3, m) + 1) * (o12(6, m) - o12(5, m) + 1)

          Rowptr_sizes(m) = block_fine%ncells + 1
        ELSE
          overlaps(m) = .FALSE.
        END IF

        IF (overlaps(m) .AND. .NOT. ANY(jblock .EQ. newsubdomain%blockids_comp)) THEN
          blockrank = blockranks(jblock)
          CALL overlap_blockids_transinds(blockrank + 1)%lst%append(elm%ivalue)
          CALL jsblocks_to_la(blockrank + 1)%lst%append(jblock)
          CALL blockbnds_to_la(blockrank + 1)%lst%append(o21(:, m))
          elm => ind_trans%get(bind)
          CALL send_transinds(blockrank + 1)%lst%append(elm%ivalue)
        END IF
        m = m  + 1
      END DO
    END IF
  END DO

  CALL make_prol(Prol_fo, subdomain_fine, neighbor_blocks_unique, neighbor_ids_unique, neighbor_ids, &
                   neighbor_ids_transform, overlaps, o12, o21, nnzs, Rowptr_sizes, order=1)

  CALL make_prol(Prol_so, subdomain_fine, neighbor_blocks_unique, neighbor_ids_unique, neighbor_ids, &
                   neighbor_ids_transform, overlaps, o12, o21, 8 * nnzs, Rowptr_sizes, order=2)

  CALL SpTrans_SpRowCol(Restr, Prol_fo)  

  nnz = 0

  DO ip = 1, NumProcs
    nnz = 0
    DO m = 1, jsblocks_to_la(ip)%lst%len
      elm => blockbnds_to_la(ip)%lst%get(m)
      o21_ptr => elm%iarray1d
      nnz = nnz + (o21_ptr(2) - o21_ptr(1) + 1) * (o21_ptr(4) - o21_ptr(3) + 1) * (o21_ptr(6) - o21_ptr(5) + 1)
    END DO
    CALL jsblocks_to_la(ip)%lst%count_deep()
    IF (jsblocks_to_la(ip)%lst%len_flat .GE. 1) THEN
      CALL jsblocks_to_la(ip)%lst%unique(jsblocks_to(ip)%data)
    ELSE
      ALLOCATE(jsblocks_to(ip)%data(0))
    END IF 
    CALL restr_offrank_indsto(jsblocks_to(ip)%data, neighbor_ids_unique, neighbor_blocks_unique, jsblocks_to_la(ip)%lst, &
                              blockbnds_to_la(ip)%lst, nnz, indices_unique_to(ip)%data, block_ptr_to(ip)%data)

  END DO

  CALL make_cut_samerank(Cut_samerank, Restr, newsubdomain, neighbor_blocks_unique, neighbor_ids_unique)
  CALL make_cut_offrank(Cut_offrank, Restr, jsblocks_to, indices_unique_to, block_ptr_to, &
                        neighbor_blocks_unique, neighbor_ids_unique)

  CALL SpTrans_SpRowCol(Cut_samerank_T, Cut_samerank)
  CALL SpTrans_SpRowCol(Cut_offrank_T, Cut_offrank)

  sendbufint(:) = 0
  DO i = 1, NumProcs
    sendbufint(i) = SIZE(jsblocks_to(i)%data, 1)
  END DO

  CALL MPI_ALLTOALL(sendbufint, 1, MPI_INT, recvbufint, 1, MPI_INT, MPI_COMM_WORLD, MPIErr)

  DO i = 1, NumProcs
    tag_send = newsubdomain%myrank * NumProcs + i - 1
    tag_recv = (i - 1) * NumProcs + newsubdomain%myrank

    ALLOCATE(jsblocks_from(i)%data(recvbufint(i)))
    ALLOCATE(reqs_send_arr(i)%data(recvbufint(i)))
    ALLOCATE(reqs_recv_arr(i)%data(sendbufint(i)))

    CALL MPI_ISEND(jsblocks_to(i)%data, SIZE(jsblocks_to(i)%data, 1), MPI_INT, i - 1, &
                   tag_send, MPI_COMM_WORLD, reqs_send(i), MPIErr)

    CALL MPI_IRECV(jsblocks_from(i)%data, recvbufint(i), MPI_INT, i - 1, &
                   tag_recv, MPI_COMM_WORLD, reqs_recv(i), MPIErr)
  END DO

  DO i = 1, NumProcs
    CALL MPI_WAIT(reqs_send(i), status, MPIErr)
    CALL MPI_WAIT(reqs_recv(i), status, MPIErr)
  END DO

  DO i = 1, NumProcs
    sendbufint(i) = SIZE(indices_unique_to(i)%data, 1)
  END DO

  CALL MPI_ALLTOALL(sendbufint, 1, MPI_INT, recvbufint, 1, MPI_INT, MPI_COMM_WORLD, MPIErr)

  DO i = 1, NumProcs
    tag_send = newsubdomain%myrank * NumProcs + i - 1
    tag_recv = (i - 1) * NumProcs + newsubdomain%myrank

    ALLOCATE(indices_unique_from(i)%data(recvbufint(i)))

    CALL MPI_ISEND(indices_unique_to(i)%data, SIZE(indices_unique_to(i)%data, 1), MPI_INT, i - 1, &
                   tag_send, MPI_COMM_WORLD, reqs_send(i), MPIErr)

    CALL MPI_IRECV(indices_unique_from(i)%data, recvbufint(i), MPI_INT, i - 1, &
                   tag_recv, MPI_COMM_WORLD, reqs_recv(i), MPIErr)
  END DO

  DO i = 1, NumProcs
    CALL MPI_WAIT(reqs_send(i), status, MPIErr)
    CALL MPI_WAIT(reqs_recv(i), status, MPIErr)
  END DO

  DO i = 1, NumProcs
    tag_send = newsubdomain%myrank * NumProcs + i - 1
    tag_recv = (i - 1) * NumProcs + newsubdomain%myrank
    
    ALLOCATE(block_ptr_from(i)%data(SIZE(jsblocks_from(i)%data, 1) + 1))

    CALL MPI_ISEND(block_ptr_to(i)%data, SIZE(block_ptr_to(i)%data, 1), MPI_INT, i - 1, &
                   tag_send, MPI_COMM_WORLD, reqs_send(i), MPIErr)
  
    CALL MPI_IRECV(block_ptr_from(i)%data, SIZE(block_ptr_from(i)%data, 1), MPI_INT, i - 1, &
                   tag_recv, MPI_COMM_WORLD, reqs_recv(i), MPIErr)
  END DO 
    
  DO i = 1, NumProcs
    CALL MPI_WAIT(reqs_send(i), status, MPIErr)
    CALL MPI_WAIT(reqs_recv(i), status, MPIErr)
  END DO

  ALLOCATE(subdomain_fine%tr_if)
  new_restrcomm_if => subdomain_fine%tr_if

  ALLOCATE(new_restrcomm_if%commbuff(Cut_offrank%m))

  nranks_comm = 0
  DO ip = 1, NumProcs
    IF (SIZE(jsblocks_to(ip)%data, 1) .GE. 1) THEN
      nranks_comm = nranks_comm + 1
    END IF
  END DO

  ALLOCATE(new_restrcomm_if%rank_ptr(nranks_comm + 1))
  new_restrcomm_if%rank_ptr(1) = 1

  i = 1
  DO ip = 1, NumProcs
    IF (SIZE(jsblocks_to(ip)%data, 1) .GE. 1) THEN
      IF (.NOT. ASSOCIATED(new_restrcomm_if%send_bnds)) THEN
        ALLOCATE(new_restrcomm_if%send_bnds)
        current => new_restrcomm_if%send_bnds
      ELSE
        ALLOCATE(current%next)
        current => current%next
      END IF
      current%rank_comm = ip - 1
      current%tag = (ip - 1) * NumProcs + subdomain_fine%myrank

      new_restrcomm_if%rank_ptr(i + 1) = new_restrcomm_if%rank_ptr(i) + SIZE(indices_unique_to(ip)%data, 1)

      current%fld_data => new_restrcomm_if%commbuff(new_restrcomm_if%rank_ptr(i):new_restrcomm_if%rank_ptr(i + 1) - 1)
      i = i + 1      
    END IF
  END DO

  DO ip = 1, NumProcs
    IF (SIZE(jsblocks_from(ip)%data, 1) .GE. 1) THEN
      IF (.NOT. ASSOCIATED(new_restrcomm_if%recv_bnds)) THEN
        ALLOCATE(new_restrcomm_if%recv_bnds)
        current => new_restrcomm_if%recv_bnds
      ELSE
        ALLOCATE(current%next)
        current => current%next
      END IF
      current%rank_comm = ip - 1
      current%tag = newsubdomain%myrank * NumProcs + (ip - 1)
      data_size = SIZE(indices_unique_from(ip)%data, 1)

      ALLOCATE(current%fld_inds(data_size))
      ALLOCATE(current%fld_data(data_size))
    
      idata = 1
      DO m = 1, SIZE(jsblocks_from(ip)%data, 1)        
        iblock = FINDLOC(newsubdomain%blockids_compghst, jsblocks_from(ip)%data(m), 1)
        thisblock => newsubdomain%blocks(iblock)

        icell_st = 0
        DO jblock = 1, iblock - 1
          IF (newsubdomain%blockiscomp(jblock) .EQV. .True.) THEN
            icell_st = icell_st + newsubdomain%blocks(jblock)%ncells
          END IF
        END DO

        DO i = block_ptr_from(ip)%data(m), block_ptr_from(ip)%data(m + 1) - 1
          current%fld_inds(idata) = indices_unique_from(ip)%data(i) + icell_st 
          idata = idata + 1
        END DO
      END DO
    END IF
  END DO

  CALL restr_scaling(Restr, Cut_samerank, Cut_offrank, subdomain_fine, newsubdomain)
  CALL SpMm_SpRowCol(subdomain_fine%Restr_samerank, Cut_samerank, Restr)
  CALL SpMm_SpRowCol(subdomain_fine%Prol_samerank, Prol_so, Cut_samerank_T)
  CALL SpMm_SpRowCol(new_restrcomm_if%Restr_offrank, Cut_offrank, Restr)
  CALL SpMm_SpRowCol(new_restrcomm_if%Prol_offrank, Prol_so, Cut_offrank_T)

END SUBROUTINE init_transfer


SUBROUTINE make_prol(Prol, subdomain_fine, blocks_coarse, neighbor_ids_unique, neighbor_ids, &
                       neighbor_ids_transform, overlaps, o12, o21, nnzs, Rowptr_sizes, order)
  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: Prol
  TYPE(SubDomain), INTENT(in) :: subdomain_fine
  TYPE(Block), POINTER, INTENT(in) :: blocks_coarse(:)
  INTEGER, INTENT(in) :: neighbor_ids_unique(:)
  TYPE(List), INTENT(in) :: neighbor_ids, neighbor_ids_transform
  LOGICAL, INTENT(in) :: overlaps(:)
  INTEGER, INTENT(in), DIMENSION(:,:) :: o21, o12
  INTEGER, INTENT(in), DIMENSION(:) :: nnzs, Rowptr_sizes
  INTEGER, INTENT(in) :: order

  INTEGER :: proltmp_colind(8, subdomain_fine%ncells)
  REAL(Realkind) :: proltmp_val(8, subdomain_fine%ncells)
  INTEGER :: colindptr(subdomain_fine%ncells)

  REAL(Realkind), TARGET :: Prolsub_val(MAXVAL(nnzs, 1))
  INTEGER, TARGET :: Prolsub_colind(MAXVAL(nnzs, 1))
  INTEGER, TARGET :: Prolsub_rowptr(subdomain_fine%ncells + 1)

  TYPE(SpRowCol) :: Prolsub

  INTEGER :: n, m, i, j, ii, iblock_fine, iblock_coarse, irow, ncol_off, icol
  INTEGER :: ncells_fine, ncells_coarse_z, ncells_coarse_zy, ncells_coarse_zyx, &
             prolz_nnz_ub, proly_nnz_ub, prolx_nnz_ub, prolxy_nnz_ub, prol_nnz_ub

  TYPE(Block), POINTER :: block_fine => NULL(), &
                          block_coarse => NULL()
  TYPE(Element), POINTER :: elm => NULL()
  TYPE(List), POINTER :: ind_trans, possible_overlapblock_ids

  Prolsub%Val => Prolsub_val
  Prolsub%colInd => Prolsub_colind

  m = 1
  iblock_fine = 1
  irow = 1
  colindptr(:) = 0

  Prol%n = 0
  DO n = 1, SIZE(blocks_coarse, 1)
    Prol%n = Prol%n + blocks_coarse(n)%ncells
  END DO

  DO n = 1, subdomain_fine%nblocks
    IF (subdomain_fine%blockiscomp(n)) THEN

      block_fine => subdomain_fine%blocks(n)

      elm => neighbor_ids%get(iblock_fine)
      possible_overlapblock_ids => elm%sublist

      elm => neighbor_ids_transform%get(iblock_fine)
      ind_trans => elm%sublist
      iblock_fine = iblock_fine + 1

      DO iblock_coarse = 1, possible_overlapblock_ids%len
        elm => ind_trans%get(iblock_coarse)
        block_coarse => blocks_coarse(elm%ivalue)

        IF (overlaps(m)) THEN

          ncells_fine = block_fine%ncells
          ncells_coarse_z = block_fine%fld_shape(3) * block_fine%fld_shape(2) * block_coarse%fld_shape(1)
          ncells_coarse_zy = block_fine%fld_shape(3) * block_coarse%fld_shape(2) * block_coarse%fld_shape(1)
          ncells_coarse_zyx = block_coarse%ncells

          prolz_nnz_ub = block_fine%fld_shape(3) * block_fine%fld_shape(2) * (o12(6, m) - o12(5, m) + 1)
          proly_nnz_ub = block_fine%fld_shape(3) * (o12(4, m) - o12(3, m) + 1) * block_coarse%fld_shape(1)
          prolx_nnz_ub = (o12(2, m) - o12(1, m) + 1) * block_coarse%fld_shape(2) * block_coarse%fld_shape(1)
          prolxy_nnz_ub = (o12(2, m) - o12(1, m) + 1) * (o12(4, m) - o12(3, m) + 1) * block_coarse%fld_shape(1)
          prol_nnz_ub = (o12(2, m) - o12(1, m) + 1) * (o12(4, m) - o12(3, m) + 1) * (o12(6, m) - o12(5, m) + 1)

          ncol_off = 0
          DO j = 1, elm%ivalue - 1
            ncol_off = ncol_off + blocks_coarse(j)%ncells
          END DO

          Prolsub%RowPtr => Prolsub_rowptr(1:Rowptr_sizes(m))

          IF (order .EQ. 1) THEN

            CALL make_prol_fo_part(Prolsub, block_fine, block_coarse, &
                                   prolz_nnz_ub, proly_nnz_ub, prolx_nnz_ub, prolxy_nnz_ub, prol_nnz_ub, &
                                   ncells_fine, ncells_coarse_z, ncells_coarse_zy, ncells_coarse_zyx)
          ELSE IF (order .EQ. 2) THEN
            CALL make_prol_so_part(Prolsub, block_fine, blocks_coarse, elm%ivalue, &
                                   2 * prolz_nnz_ub, 2 * proly_nnz_ub, 2 * prolx_nnz_ub, 4 * prolxy_nnz_ub, 8 * prol_nnz_ub, &
                                   ncells_fine, ncells_coarse_z, ncells_coarse_zy, ncells_coarse_zyx)
          END IF

          DO i = 1, Prolsub%m
            ii = i + irow - 1
            DO j = Prolsub%RowPtr(i), Prolsub%RowPtr(i + 1) - 1
              colindptr(ii) = colindptr(ii) + 1
              proltmp_colind(colindptr(ii), ii) = Prolsub%ColInd(j) + ncol_off
              proltmp_val(colindptr(ii), ii) = Prolsub%Val(j)
            END DO
          END DO
        END IF
        m = m + 1
      END DO
      irow = irow + block_fine%ncells
    END IF
  END DO

  Prol%m = irow - 1
  ALLOCATE(Prol%RowPtr(irow))
  ALLOCATE(Prol%ColInd(SUM(colindptr)))
  ALLOCATE(Prol%Val(SUM(colindptr)))

  icol = 1
  Prol%RowPtr(1) = icol
  DO i = 1, subdomain_fine%ncells
    CALL sort(proltmp_colind(1:colindptr(i), i), proltmp_val(1:colindptr(i), i))
    DO j = 1, colindptr(i)
      Prol%ColInd(icol) = proltmp_colind(j, i)
      Prol%Val(icol) = proltmp_val(j, i)
      icol = icol + 1
    END DO
    Prol%RowPtr(i + 1) = icol
  END DO

  CALL SpNullify_SpRowCol(Prolsub)

END SUBROUTINE make_prol


SUBROUTINE make_prol_fo_part(Prol, block_fine, block_coarse, &
                             prolz_nnz, proly_nnz, prolx_nnz, prolxy_nnz, prol_nnz, &
                             ncells_fine, ncells_coarse_z, ncells_coarse_zy, ncells_coarse_zyx)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: Prol
  TYPE(Block), INTENT(in) :: block_fine, block_coarse
  INTEGER, INTENT(in) :: prolz_nnz, proly_nnz, prolx_nnz, prolxy_nnz, prol_nnz, &
                         ncells_fine, ncells_coarse_z, ncells_coarse_zy, ncells_coarse_zyx

  TYPE(SpRowCol) :: Prolz, Proly, Prolx, Proltmp

  REAl(Realkind), TARGET :: prolz_val(prolz_nnz), &
                            proly_val(proly_nnz), &
                            prolx_val(prolx_nnz), &
                            proltmp_val(prolxy_nnz)

  INTEGER, TARGET :: prolz_colind(prolz_nnz), &
                     proly_colind(proly_nnz), &
                     prolx_colind(prolx_nnz), &
                     proltmp_colind(prolxy_nnz)

  INTEGER, TARGET :: prolz_indptr(ncells_fine + 1), &
                     proly_indptr(ncells_coarse_z + 1), &
                     prolx_indptr(ncells_coarse_zy + 1), &
                     proltmp_indptr(ncells_coarse_z + 1)

  REAL(Realkind), POINTER :: x(:), y(:), z(:)
  REAL(Realkind), POINTER :: x2_coarse(:), y2_coarse(:), z2_coarse(:)

  INTEGER :: nx, ny, nz, nx_coarse, ny_coarse, nz_coarse

  INTEGER :: m, iindptr, ind, c
  INTEGER :: k, i, j, k_ind, j_ind, i_ind

  Prolz%m = ncells_fine
  Prolz%n = ncells_coarse_z
  Proly%m = ncells_coarse_z
  Proly%n = ncells_coarse_zy
  Prolx%m = ncells_coarse_zy
  Prolx%n = ncells_coarse_zyx

  x => block_fine%x
  y => block_fine%y
  z => block_fine%z

  x2_coarse => block_coarse%x2
  y2_coarse => block_coarse%y2
  z2_coarse => block_coarse%z2

  nz = block_fine%fld_shape(1)
  ny = block_fine%fld_shape(2)
  nx = block_fine%fld_shape(3)

  nz_coarse = block_coarse%fld_shape(1)
  ny_coarse = block_coarse%fld_shape(2)
  nx_coarse = block_coarse%fld_shape(3)

  prolz_val(:) = 1.0
  iindptr = 1
  prolz_indptr(1) = 1
  c = 1


  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz
        IF (z(k) < z2_coarse(1).OR. z(k) > z2_coarse(UBOUND(z2_coarse, 1))) THEN
          iindptr = iindptr + 1
          prolz_indptr(iindptr) = c
        ELSE
          k_ind = between(z2_coarse, z(k))
          ind = IndexC(i, j, k_ind, nx, ny, nz_coarse)
          prolz_colind(c) = ind
          c = c + 1
          iindptr = iindptr + 1
          prolz_indptr(iindptr) = c
        END IF
      END DO
    END DO
  END DO

  Prolz%RowPtr => prolz_indptr(1:iindptr)
  Prolz%Val => prolz_val(1:Prolz%RowPtr(iindptr) - 1)
  Prolz%ColInd => prolz_colind(1:Prolz%RowPtr(iindptr) - 1)

  proly_val(:) = 1.0
  iindptr = 1
  proly_indptr(1) = 1
  c = 1

  DO i = 1, nx
    DO j = 1, ny
      IF (y(j) < y2_coarse(1).OR. y(j) > y2_coarse(UBOUND(y2_coarse, 1))) THEN
        DO k = 1, nz_coarse
          iindptr = iindptr + 1
          proly_indptr(iindptr) = c
        END DO
      ELSE
        j_ind = between(y2_coarse, y(j))
        DO k = 1, nz_coarse
          ind = IndexC(i, j_ind, k, nx, ny_coarse, nz_coarse)
          proly_colind(c) = ind
          c = c + 1
          iindptr = iindptr + 1
          proly_indptr(iindptr) = c
        END DO
      END IF
    END DO
  END DO

  Proly%RowPtr => proly_indptr(1:iindptr)
  Proly%Val => proly_val(1:c - 1)
  Proly%ColInd => proly_colind(1:c - 1)

  prolx_val(:) = 1.0
  iindptr = 1
  prolx_indptr(1) = 1
  c = 1

  DO i = 1, nx
    IF (x(i) < x2_coarse(1).OR. x(i) > x2_coarse(UBOUND(x2_coarse, 1))) THEN
      DO j = 1, ny_coarse
        DO k = 1, nz_coarse
          iindptr = iindptr + 1
          prolx_indptr(iindptr) = c
        END DO
      END DO
    ELSE
      i_ind = between(x2_coarse, x(i))
      DO j = 1, ny_coarse
        DO k = 1, nz_coarse
          ind = IndexC(i_ind, j, k, nx_coarse, ny_coarse, nz_coarse)
          prolx_colind(c) = ind
          c = c + 1
          iindptr = iindptr + 1
          prolx_indptr(iindptr) = c
        END DO
      END DO
    END IF
  END DO

  Prolx%RowPtr => prolx_indptr(1:iindptr)
  Prolx%Val => prolx_val(1:c - 1)
  Prolx%ColInd => prolx_colind(1:c - 1)

  Proltmp%RowPtr => proltmp_indptr
  Proltmp%ColInd => proltmp_colind
  Proltmp%Val => proltmp_val

  CALL SpMm_SpRowCol(Proltmp, Proly, Prolx, allocate_mat=.False., withdiagptr=.False.)
  CALL SpMm_SpRowCol(Prol, Prolz, Proltmp, allocate_mat=.False., withdiagptr=.False.)

  CALL SpNullify_SpRowCol(Prolz)
  CALL SpNullify_SpRowCol(Proly)
  CALL SpNullify_SpRowCol(Prolx)
  CALL SpNullify_SpRowCol(Proltmp)

END SUBROUTINE make_prol_fo_part


SUBROUTINE make_prol_so_part(Prol, block_fine, blocks_coarse, iblock, &
                             prolz_nnz, proly_nnz, prolx_nnz, prolxy_nnz, prol_nnz, &
                             ncells_fine, ncells_coarse_z, ncells_coarse_zy, ncells_coarse_zyx)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: Prol
  TYPE(Block), INTENT(in) :: block_fine
  TYPE(Block), TARGET, INTENT(in) :: blocks_coarse(:)
  INTEGER, INTENT(in) :: iblock

  INTEGER, INTENT(in) :: prolz_nnz, proly_nnz, prolx_nnz, prolxy_nnz, prol_nnz, &
                         ncells_fine, ncells_coarse_z, ncells_coarse_zy, ncells_coarse_zyx

  TYPE(SpRowCol) :: Prolz, Proly, Prolx, Proltmp, Proltmp2

  REAl(Realkind), TARGET :: prolz_val(prolz_nnz), &
                            proly_val(proly_nnz), &
                            prolx_val(prolx_nnz), &
                            proltmp_val(prolxy_nnz)

  INTEGER, TARGET :: prolz_colind(prolz_nnz), &
                     proly_colind(proly_nnz), &
                     prolx_colind(prolx_nnz), &
                     proltmp_colind(prolxy_nnz) 
                     
  INTEGER, TARGET :: prolz_indptr(ncells_fine + 1), &
                     proly_indptr(ncells_coarse_z + 1), &
                     prolx_indptr(ncells_coarse_zy + 1), &
                     proltmp_indptr(ncells_coarse_z + 1)

  REAL(Realkind), POINTER :: x(:), y(:), z(:), x2(:), y2(:), z2(:)
  REAL(Realkind), POINTER :: x_coarse(:), y_coarse(:), z_coarse(:), &
                             x2_coarse(:), y2_coarse(:), z2_coarse(:)
  REAL(Realkind), POINTER, DIMENSION(:) :: x_coarse_tmp, y_coarse_tmp, z_coarse_tmp
  REAL(Realkind) :: z_coarse_bnd, y_coarse_bnd, x_coarse_bnd
  REAL(Realkind) :: val1, val2

  TYPE(Block), POINTER :: block_coarse

  INTEGER :: nx, ny, nz, nx_coarse, ny_coarse, nz_coarse

  INTEGER :: i_bnds(2), j_bnds(2), k_bnds(2)  
  INTEGER :: m, iindptr, ind, c, bndblock_id
  INTEGER :: k, i, j, k_ind, j_ind, i_ind
  
  Prolz%m = ncells_fine
  Prolz%n = ncells_coarse_z
  Proly%m = ncells_coarse_z
  Proly%n = ncells_coarse_zy
  Prolx%m = ncells_coarse_zy
  Prolx%n = ncells_coarse_zyx
 
  block_coarse => blocks_coarse(iblock)
 
  x => block_fine%x
  y => block_fine%y
  z => block_fine%z

  x2 => block_fine%x2
  y2 => block_fine%y2
  z2 => block_fine%z2  

  x_coarse => block_coarse%x
  y_coarse => block_coarse%y
  z_coarse => block_coarse%z

  x2_coarse => block_coarse%x2
  y2_coarse => block_coarse%y2
  z2_coarse => block_coarse%z2
  
  nz = block_fine%fld_shape(1)
  ny = block_fine%fld_shape(2)
  nx = block_fine%fld_shape(3)
  
  nz_coarse = block_coarse%fld_shape(1)
  ny_coarse = block_coarse%fld_shape(2)
  nx_coarse = block_coarse%fld_shape(3)
  
  prolz_val(:) = 1.0
  iindptr = 1
  prolz_indptr(1) = 1
  c = 1

  DO i = 1, nx
    DO j = 1, ny
      DO k = 1, nz
        IF (z(k) .LT. z_coarse(1)) THEN
          IF (z2(k + 1) .LT. z2_coarse(1)) THEN
            iindptr = iindptr + 1
            prolz_indptr(iindptr) =  c
          ELSE IF (ANY(block_coarse%cface_b)) THEN
            CALL find_bndblock_ij(blocks_coarse, block_coarse%cface_b, x2(i), x2(i + 1), y2(j), y2(j + 1), iblock, bndblock_id)
            IF (bndblock_id .LT. HUGE(bndblock_id)) THEN
              z_coarse_bnd = blocks_coarse(bndblock_id)%z(UBOUND(blocks_coarse(bndblock_id)%z, 1))
              IF (z_coarse_bnd .NE. z(k)) THEN
                val1 = (z(k) - z_coarse_bnd) / (z_coarse(1) - z_coarse_bnd)
                ind = IndexC(i, j, 1, nx, ny, nz_coarse)
                prolz_colind(c) = ind
                prolz_val(c) = val1
                c = c + 1
              END IF
            END IF
            iindptr = iindptr + 1
            prolz_indptr(iindptr) =  c
          ELSE !Neumann boundary 
            val1 = 1.0
            ind = IndexC(i, j, 1, nx, ny, nz_coarse)
            prolz_colind(c) = ind
            prolz_val(c) = val1
            c = c + 1
            iindptr = iindptr + 1
            prolz_indptr(iindptr) =  c
          END IF
        ELSE IF (z(k) > z_coarse(UBOUND(z_coarse, 1))) THEN
          IF (z2(k) > z2_coarse(UBOUND(z2_coarse, 1))) THEN 
            iindptr = iindptr + 1
            prolz_indptr(iindptr) =  c
          ELSE IF (ANY(block_coarse%cface_t)) THEN
            CALL find_bndblock_ij(blocks_coarse, block_coarse%cface_t, x2(i), x2(i + 1), y2(j), y2(j + 1), iblock, bndblock_id)
            IF (bndblock_id .LT. HUGE(bndblock_id)) THEN
              z_coarse_bnd = blocks_coarse(bndblock_id)%z(1)
              IF (z_coarse_bnd .NE. z(k)) THEN
                val1 = (z_coarse_bnd - z(k)) / (z_coarse_bnd - z_coarse(UBOUND(z_coarse, 1)))
                ind = IndexC(i, j, UBOUND(z_coarse, 1), nx, ny, nz_coarse)
                prolz_colind(c) = ind
                prolz_val(c) = val1
                c = c + 1
              END IF
            END IF
            iindptr = iindptr + 1
            prolz_indptr(iindptr) =  c
          ELSE !Neumann boundary 
            val1 = 1.0
            ind = IndexC(i, j, UBOUND(z_coarse, 1), nx, ny, nz_coarse)
            prolz_colind(c) = ind
            prolz_val(c) = val1
            c = c + 1 
            iindptr = iindptr + 1
            prolz_indptr(iindptr) =  c
          END IF
        ELSE IF ((z(k) .EQ. z_coarse(1)) .AND. (z_coarse(1) .EQ. z_coarse(UBOUND(z_coarse, 1)))) THEN
          ind = IndexC(i, j, 1, nx, ny, nz_coarse)
          prolz_colind(c) = ind
          prolz_val(c) = 1.0
          c = c + 1 
          iindptr = iindptr + 1
          prolz_indptr(iindptr) =  c
        ELSE
          k_ind = between(z_coarse, z(k))
          IF (z_coarse(k_ind) .EQ. z(k)) THEN
            val1 = 1.0
            val2 = 0.0
          ELSE
            val1 = (z_coarse(k_ind + 1) - z(k)) / (z_coarse(k_ind + 1) - z_coarse(k_ind))
            val2 = 1.0 - val1
          END IF

          ind = IndexC(i, j, k_ind, nx, ny, nz_coarse)
          prolz_colind(c) = ind
          prolz_val(c) = val1
          c = c + 1 
          IF (val2 .NE. 0.0) THEN
            ind = IndexC(i, j, k_ind + 1, nx, ny, nz_coarse)
            prolz_colind(c) = ind
            prolz_val(c) = val2
            c = c + 1 
          END IF
          iindptr = iindptr + 1
          prolz_indptr(iindptr) =  c
        END IF
      END DO
    END DO
  END DO

  Prolz%RowPtr => prolz_indptr(1:iindptr)
  Prolz%Val => prolz_val(1:Prolz%RowPtr(iindptr) - 1)
  Prolz%ColInd => prolz_colind(1:Prolz%RowPtr(iindptr) - 1)

  proly_val(:) = 1.0
  iindptr = 1
  proly_indptr(1) = 1
  c = 1

  DO i = 1, nx
    DO j = 1, ny
      IF (y(j) .LT. y_coarse(1)) THEN
        IF (y2(j + 1) .LT. y2_coarse(1)) THEN
          DO k = 1, nz_coarse
            iindptr = iindptr + 1
            proly_indptr(iindptr) =  c
          END DO
        ELSE IF (ANY(block_coarse%cface_s)) THEN
          DO k = 1, nz_coarse
            CALL find_bndblock_ik(blocks_coarse, block_coarse%cface_s, x2(i), x2(i + 1), &
                                   z2_coarse(k), z2_coarse(k + 1), iblock, bndblock_id)
            IF (bndblock_id .LT. HUGE(bndblock_id)) THEN
              y_coarse_bnd = blocks_coarse(bndblock_id)%y(UBOUND(blocks_coarse(bndblock_id)%y, 1))
              IF (y_coarse_bnd .NE. y(j)) THEN
                val1 = (y(j) - y_coarse_bnd) / (y_coarse(1) - y_coarse_bnd)
                ind = IndexC(i, 1, k, nx, ny_coarse, nz_coarse)
                proly_colind(c) = ind
                Proly_val(c) = val1
                c = c + 1
              END IF
            END IF
            iindptr = iindptr + 1
            proly_indptr(iindptr) =  c
          END DO
        ELSE !Neumann boundary
          DO k = 1, nz_coarse
            val1 = 1.0
            ind = IndexC(i, 1, k, nx, ny_coarse, nz_coarse)
            proly_colind(c) = ind
            Proly_val(c) = val1
            c = c + 1
            iindptr = iindptr + 1
            proly_indptr(iindptr) =  c
          END DO
        END IF
      ELSE IF (y(j) .GT. y_coarse(UBOUND(y_coarse, 1))) THEN
        IF (y2(j) .GT. y2_coarse(UBOUND(y2_coarse, 1))) THEN
          DO k = 1, nz_coarse
            iindptr = iindptr + 1
            proly_indptr(iindptr) =  c
          END DO
        ELSE IF (ANY(block_coarse%cface_n)) THEN
          DO k = 1, nz_coarse
            CALL find_bndblock_ik(blocks_coarse, block_coarse%cface_n, x2(i), x2(i + 1), &
                                  z2_coarse(k), z2_coarse(k + 1), iblock, bndblock_id)
            IF (bndblock_id .LT. HUGE(bndblock_id)) THEN
              y_coarse_bnd = blocks_coarse(bndblock_id)%y(1)
              IF (y_coarse_bnd .NE. y(j)) THEN
                val1 = (y_coarse_bnd - y(j)) / (y_coarse_bnd - y_coarse(UBOUND(y_coarse, 1)))
                ind = IndexC(i, UBOUND(y_coarse, 1), k, nx, ny_coarse, nz_coarse)
                proly_colind(c) = ind
                Proly_val(c) = val1
                c = c + 1
              END IF
            END IF
            iindptr = iindptr + 1
            proly_indptr(iindptr) =  c
          END DO
        ELSE !Neumann boundary
          DO k = 1, nz_coarse
            val1 = 1.0
            ind = IndexC(i, UBOUND(y_coarse, 1), k, nx, ny_coarse, nz_coarse)
            proly_colind(c) = ind
            Proly_val(c) = val1
            c = c + 1
            iindptr = iindptr + 1
            proly_indptr(iindptr) =  c
          END DO
        END IF

      ELSE IF ((y(j) .EQ. y_coarse(1)) .AND. (y_coarse(1) .EQ. y_coarse(UBOUND(y_coarse, 1)))) THEN
        DO k = 1, nz_coarse
          ind = IndexC(i, 1, k, nx, ny_coarse, nz_coarse)
          proly_colind(c) = ind
          Proly_val(c) = val1
          c = c + 1
          iindptr = iindptr + 1
          proly_indptr(iindptr) =  c
        END DO
      ELSE
        j_ind = between(y_coarse, y(j))
        IF (y_coarse(j_ind) .EQ. y(j)) THEN
          val1 = 1.0
          val2 = 0.0
        ELSE
          val1 = (y_coarse(j_ind + 1) - y(j)) / (y_coarse(j_ind + 1) - y_coarse(j_ind))
          val2 = 1.0 - val1
        END IF
        DO k = 1, nz_coarse
          ind = IndexC(i, j_ind, k, nx, ny_coarse, nz_coarse)
          proly_colind(c) = ind
          Proly_val(c) = val1
          c = c + 1
          IF (val2 .NE. 0.0) THEN
            ind = IndexC(i, j_ind + 1, k, nx, ny_coarse, nz_coarse)
            proly_colind(c) = ind
            Proly_val(c) = val2
            c = c + 1
          END IF
          iindptr = iindptr + 1
          proly_indptr(iindptr) =  c
        END DO
      END IF
    END DO
  END DO

  Proly%RowPtr => proly_indptr(1:iindptr)
  Proly%Val => proly_val(1:c - 1)
  Proly%ColInd => proly_colind(1:c - 1)

  prolx_val(:) = 1.0
  iindptr = 1
  prolx_indptr(1) = 1
  c = 1

  DO i = 1, nx
    IF (x(i) .LT. x_coarse(1)) THEN
      IF (x2(i + 1) .LT. x2_coarse(1)) THEN
        DO j = 1, ny_coarse
          DO k = 1, nz_coarse
            iindptr = iindptr + 1
            Prolx_indptr(iindptr) =  c
          END DO
        END DO
      ELSE IF (ANY(block_coarse%cface_w)) THEN
        DO j = 1, ny_coarse
          DO k = 1, nz_coarse
            CALL find_bndblock_jk(blocks_coarse, block_coarse%cface_w, y2_coarse(j), y2_coarse(j + 1), &
                                  z2_coarse(k), z2_coarse(k + 1), iblock, bndblock_id)
            IF (bndblock_id .LT. HUGE(bndblock_id)) THEN
              x_coarse_bnd = blocks_coarse(bndblock_id)%x(UBOUND(blocks_coarse(bndblock_id)%x, 1))
              IF (x_coarse_bnd .NE. x(i)) THEN
                val1 = (x(i) - x_coarse_bnd) / (x_coarse(1) - x_coarse_bnd)
                ind = IndexC(1, j, k, nx_coarse, ny_coarse, nz_coarse)
                Prolx_colind(c) = ind
                Prolx_val(c) = val1
                c = c + 1
              END IF
            END IF
            iindptr = iindptr + 1
            Prolx_indptr(iindptr) =  c
          END DO
        END DO
      ELSE !Neumann boundary
        DO j = 1, ny_coarse
          DO k = 1, nz_coarse
            val1 = 1.0 
            ind = IndexC(1, j, k, nx_coarse, ny_coarse, nz_coarse)
            Prolx_colind(c) = ind
            Prolx_val(c) = val1
            c = c + 1
            iindptr = iindptr + 1
            Prolx_indptr(iindptr) =  c
          END DO
        END DO
      END IF
    ELSE IF (x(i) .GT. x_coarse(UBOUND(x_coarse, 1))) THEN
      IF (x2(i) .GT. x2_coarse(UBOUND(x2_coarse, 1))) THEN
        DO j = 1, ny_coarse
          DO k = 1, nz_coarse
            iindptr = iindptr + 1
            Prolx_indptr(iindptr) =  c
          END DO
        END DO
      ELSE IF (ANY(block_coarse%cface_e)) THEN
        DO j = 1, ny_coarse
          DO k = 1, nz_coarse
            CALL find_bndblock_jk(blocks_coarse, block_coarse%cface_e, y2_coarse(j), y2_coarse(j + 1), &
                                  z2_coarse(k), z2_coarse(k + 1), iblock, bndblock_id)
            IF (bndblock_id .LT. HUGE(bndblock_id)) THEN
              x_coarse_bnd = blocks_coarse(bndblock_id)%x(1)
              IF (x_coarse_bnd .NE. x(i)) THEN
                val1 = (x_coarse_bnd - x(i)) / (x_coarse_bnd - x_coarse(UBOUND(x_coarse, 1)))
                ind = IndexC(UBOUND(x_coarse, 1), j, k, nx_coarse, ny_coarse, nz_coarse)
                Prolx_colind(c) = ind
                Prolx_val(c) = val1
                c = c + 1
              END IF
            END IF
            iindptr = iindptr + 1
            Prolx_indptr(iindptr) =  c
          END DO
        END DO
      ELSE !Neumann boundary
        DO j = 1, ny_coarse
          DO k = 1, nz_coarse
            val1 = 1.0 
            ind = IndexC(UBOUND(x_coarse, 1), j, k, nx_coarse, ny_coarse, nz_coarse)
            Prolx_colind(c) = ind
            Prolx_val(c) = val1
            c = c + 1
            iindptr = iindptr + 1
            Prolx_indptr(iindptr) =  c
          END DO
        END DO
      END IF
    ELSE IF ((x(i) .EQ. x_coarse(1)) .AND. (x_coarse(1) .EQ. x_coarse(UBOUND(x_coarse, 1)))) THEN
      DO j = 1, ny_coarse
        DO k = 1, nz_coarse
          ind = IndexC(1, j, k, nx_coarse, ny_coarse, nz_coarse)
          Prolx_colind(c) = ind
          Prolx_val(c) = 1.0
          c = c + 1
          iindptr = iindptr + 1
          Prolx_indptr(iindptr) =  c
        END DO
      END DO
    ELSE
      i_ind = between(x_coarse, x(i))
      IF (x_coarse(i_ind) .EQ. x(i)) THEN
        val1 = 1.0 
        val2 = 0.0
      ELSE
        val1 = (x_coarse(i_ind + 1) - x(i)) / (x_coarse(i_ind + 1) - x_coarse(i_ind))
        val2 = 1.0 - val1
      END IF
      DO j = 1, ny_coarse
        DO k = 1, nz_coarse
          ind = IndexC(i_ind, j, k, nx_coarse, ny_coarse, nz_coarse)
          Prolx_colind(c) = ind
          Prolx_val(c) = val1
          c = c + 1
          IF (val2 .NE. 0.0) THEN
            ind = IndexC(i_ind + 1, j, k, nx_coarse, ny_coarse, nz_coarse)
            Prolx_colind(c) = ind
            Prolx_val(c) = val2
            c = c + 1
          END IF
          iindptr = iindptr + 1
          Prolx_indptr(iindptr) =  c
        END DO
      END DO
    END IF
  END DO         
  
  Prolx%RowPtr => prolx_indptr(1:iindptr)
  Prolx%Val => prolx_val(1:c - 1)
  Prolx%ColInd => prolx_colind(1:c - 1)

  Proltmp%RowPtr => proltmp_indptr
  Proltmp%ColInd => proltmp_colind
  Proltmp%Val => proltmp_val
  
  CALL SpMm_SpRowCol(Proltmp, Proly, Prolx, allocate_mat=.False., withdiagptr=.False.)
  CALL SpMm_SpRowCol(Prol, Prolz, Proltmp, allocate_mat=.False., withdiagptr=.False.)

  CALL SpNullify_SpRowCol(Prolz)
  CALL SpNullify_SpRowCol(Proly)
  CALL SpNullify_SpRowCol(Prolx)
  CALL SpNullify_SpRowCol(Proltmp)

END SUBROUTINE make_prol_so_part


SUBROUTINE make_cut_samerank(cut_samerank, Restr, newsubdomain, blockscoarse_tmp, blockids_tmp)
  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: cut_samerank
  TYPE(SpRowCol), INTENT(in) :: Restr

  TYPE(SubDomain), INTENT(in) :: newsubdomain
  TYPE(Block), INTENT(in) :: blockscoarse_tmp(:)
  INTEGER, INTENT(in) :: blockids_tmp(:)
  INTEGER :: nnz, icol, irow, iblock, jblock, icell_tmp, i, j
 
  nnz = 0
  icell_tmp = 1
  DO iblock = 1, SIZE(blockids_tmp)
    IF (ANY(blockids_tmp(iblock) .EQ. newsubdomain%blockids_comp)) THEN
      DO i = icell_tmp, icell_tmp + blockscoarse_tmp(iblock)%ncells - 1
        IF (Restr%RowPtr(i + 1) .GT. Restr%RowPtr(i)) THEN
          nnz = nnz + 1
        END IF
      END DO
    END IF
    icell_tmp = icell_tmp + blockscoarse_tmp(iblock)%ncells
  END DO
  cut_samerank%n = Restr%m

  ALLOCATE(cut_samerank%RowPtr(newsubdomain%ncells + 1))
  ALLOCATE(cut_samerank%ColInd(nnz))
  ALLOCATE(cut_samerank%Val(nnz))

  cut_samerank%m = newsubdomain%ncells

  icell_tmp = 1
  irow = 1
  icol = 1
  cut_samerank%RowPtr(1) = 1

  DO iblock = 1, newsubdomain%nblocks
    IF (newsubdomain%blockiscomp(iblock)) THEN
      IF (ANY(newsubdomain%blockids_compghst(iblock) .EQ. blockids_tmp)) THEN
        jblock = FINDLOC(blockids_tmp, newsubdomain%blockids_compghst(iblock), 1)
        icell_tmp = 1
        DO j = 1, jblock - 1
          icell_tmp = icell_tmp + blockscoarse_tmp(j)%ncells
        END DO
        DO i = icell_tmp, icell_tmp + blockscoarse_tmp(jblock)%ncells - 1
          IF (Restr%RowPtr(i + 1) .GT. Restr%RowPtr(i)) THEN
            cut_samerank%ColInd(icol) = i
            cut_samerank%Val(icol) = 1.0
            icol = icol + 1
          END IF
          cut_samerank%RowPtr(irow + 1) = icol
          irow = irow + 1
        END DO
      ELSE
        DO i = irow, irow + newsubdomain%blocks(iblock)%ncells - 1
          cut_samerank%RowPtr(i + 1) = icol
        END DO
        irow = irow + newsubdomain%blocks(iblock)%ncells
      END IF
    END IF
  END DO
  
END SUBROUTINE make_cut_samerank


SUBROUTINE make_cut_offrank(Cut_offrank, Restr, jsblocks_to, indices_unique_to, block_ptr_to, blocks, block_ids)
  IMPLICIT NONE 

  TYPE(SpRowCol), INTENT(inout) :: Cut_offrank
  TYPE(SpRowCol), INTENT(in) :: Restr
  TYPE(IArray), DIMENSION(:), INTENT(in) :: jsblocks_to, indices_unique_to, block_ptr_to
  TYPE(Block), POINTER, INTENT(in) :: blocks(:)
  INTEGER, INTENT(in) :: block_ids(:)

  INTEGER :: i, j, k, ip, NumProcs, m, nnz, iblock, jblock, iind_st, irow
  INTEGER, POINTER :: inds_block(:) => NULL()
  TYPE(Element), POINTER :: elm => NULL()
  TYPE(BLock), POINTER :: thisblock => NULL()
 

  NumProcs = SIZE(jsblocks_to, 1)

  nnz = 0

  DO ip = 1, NumProcs
    nnz = nnz + SIZE(indices_unique_to(ip)%data, 1)
  END DO

  Cut_offrank%m = nnz
  Cut_offrank%n = Restr%m
  ALLOCATE(Cut_offrank%RowPtr(nnz + 1))
  ALLOCATE(Cut_offrank%Val(nnz))
  ALLOCATE(Cut_offrank%ColInd(nnz))

  Cut_offrank%RowPtr(1) = 1
  irow = 1

  DO ip = 1, NumProcs
    DO m = 1, SIZE(jsblocks_to(ip)%data, 1)
      inds_block => indices_unique_to(ip)%data(block_ptr_to(ip)%data(m):block_ptr_to(ip)%data(m + 1) - 1)
      iblock = FINDLOC(block_ids, jsblocks_to(ip)%data(m), 1)
      thisblock => blocks(iblock)
      iind_st = 0
      DO jblock = 1, iblock - 1
        iind_st = iind_st + blocks(jblock)%ncells     
      END DO
      DO i = 1, SIZE(inds_block, 1)
        Cut_offrank%ColInd(irow) = inds_block(i) + iind_st
        Cut_offrank%Val(irow) = 1.0
        irow = irow + 1
        Cut_offrank%RowPtr(irow) = irow
      END DO
    END DO
  END DO

END SUBROUTINE make_cut_offrank


SUBROUTINE restr_scaling(Restr, Cut_samerank, Cut_offrank, subdomain_fine, newsubdomain)
  IMPLICIT NONE
  TYPE(SpRowCol), INTENT(in) :: Cut_samerank, Cut_offrank
  TYPE(SpRowCol), INTENT(inout) :: Restr     
  TYPE(SubDomain), INTENT(inout) :: newsubdomain
  TYPE(Subdomain), INTENT(in) :: subdomain_fine
  REAL(Realkind) :: vols_eff_fine1d(subdomain_fine%ncells)
  REAL(Realkind) :: vols_eff_coarse1d_offrank(newsubdomain%ncells)
  REAL(Realkind) :: vols_eff_coarse_tmp(Restr%m)
  INTEGER :: iblock, icell, i, j, k, MPIErr
  TYPE(Block), POINTER :: thisblock => NULL()
  TYPE(Bound), POINTER :: send_current => NULL(), &
                          recv_current => NULL()

  icell = 1
  DO iblock = 1, subdomain_fine%nblocks
    IF (subdomain_fine%blockiscomp(iblock)) THEN
      thisblock => subdomain_fine%blocks(iblock)
      DO i = 1, thisblock%fld_shape(3)
        DO j = 1, thisblock%fld_shape(2)
          DO k = 1, thisblock%fld_shape(1)
            vols_eff_fine1d(icell) = thisblock%volseff(k, j, i)
            icell = icell + 1
          END DO
        END DO
      END DO
    END IF
  END DO

  CALL SpAVec_SpRowCol(vols_eff_coarse_tmp, Restr, vols_eff_fine1d)  

  IF (ASSOCIATED(subdomain_fine%tr_if)) THEN
    vols_eff_coarse1d_offrank(:) = 0.0
    CALL SpAVec_SpRowCol(subdomain_fine%tr_if%commbuff, Cut_offrank, vols_eff_coarse_tmp)

    IF (ASSOCIATED(subdomain_fine%tr_if%send_bnds)) THEN
      send_current => subdomain_fine%tr_if%send_bnds
      CALL MPI_ISEND(send_current%fld_data, SIZE(send_current%fld_data), MPI_RealKind, &
                     send_current%rank_comm, send_current%tag, MPI_COMM_WORLD, send_current%req, MPIErr)

      DO WHILE(ASSOCIATED(send_current%next))
        send_current => send_current%next
        CALL MPI_ISEND(send_current%fld_data, SIZE(send_current%fld_data), MPI_RealKind, &
                       send_current%rank_comm, send_current%tag, MPI_COMM_WORLD, send_current%req, MPIErr)
      END DO
    END IF
   
    IF (ASSOCIATED(subdomain_fine%tr_if%recv_bnds)) THEN
      recv_current => subdomain_fine%tr_if%recv_bnds

      CALL MPI_IRECV(recv_current%fld_data, SIZE(recv_current%fld_data), MPI_RealKind, &
                     recv_current%rank_comm, recv_current%tag, MPI_COMM_WORLD, recv_current%req, MPIErr)

      DO WHILE(ASSOCIATED(recv_current%next))
        recv_current => recv_current%next

        CALL MPI_IRECV(recv_current%fld_data, SIZE(recv_current%fld_data), MPI_RealKind, &
                       recv_current%rank_comm, recv_current%tag, MPI_COMM_WORLD, recv_current%req, MPIErr)
      END DO
    END IF
    IF (ASSOCIATED(subdomain_fine%tr_if%send_bnds)) THEN
      send_current => subdomain_fine%tr_if%send_bnds
      CALL MPI_WAIT(send_current%req, status, MPIErr)
      DO WHILE(ASSOCIATED(send_current%next))
        send_current => send_current%next
        CALL MPI_WAIT(send_current%req, status, MPIErr)
      END DO
    END IF

    IF (ASSOCIATED(subdomain_fine%tr_if%recv_bnds)) THEN
      recv_current => subdomain_fine%tr_if%recv_bnds
      CALL MPI_WAIT(recv_current%req, status, MPIErr)
      CALL copy_data_from_bnd(recv_current, vols_eff_coarse1d_offrank, 'add')
      DO WHILE(ASSOCIATED(recv_current%next))
        recv_current => recv_current%next
        CALL MPI_WAIT(recv_current%req, status, MPIErr)
        CALL copy_data_from_bnd(recv_current, vols_eff_coarse1d_offrank, 'add')
      END DO
    END IF

    CALL ATxpy_SpRowCol(Cut_samerank, vols_eff_coarse1d_offrank, vols_eff_coarse_tmp)
    vols_eff_coarse1d_offrank(:) = 0.0
    CALL SpAVec_SpRowCol(vols_eff_coarse1d_offrank, Cut_samerank, vols_eff_coarse_tmp)
    
    IF (ASSOCIATED(subdomain_fine%tr_if%recv_bnds)) THEN
      send_current => subdomain_fine%tr_if%recv_bnds
      CALL copy_data_on_bnd(send_current, vols_eff_coarse1d_offrank, mode='rep')
      CALL MPI_ISEND(send_current%fld_data, SIZE(send_current%fld_data), MPI_RealKind, &
                     send_current%rank_comm, send_current%tag, MPI_COMM_WORLD, send_current%req, MPIErr)
      DO WHILE(ASSOCIATED(send_current%next))
        send_current => send_current%next
        CALL copy_data_on_bnd(send_current, vols_eff_coarse1d_offrank, mode='rep')
        CALL MPI_ISEND(send_current%fld_data, SIZE(send_current%fld_data), MPI_RealKind, &
                       send_current%rank_comm, send_current%tag, MPI_COMM_WORLD, send_current%req, MPIErr)
      END DO
    END IF
   
    IF (ASSOCIATED(subdomain_fine%tr_if%send_bnds)) THEN
      recv_current => subdomain_fine%tr_if%send_bnds
      CALL MPI_IRECV(recv_current%fld_data, SIZE(recv_current%fld_data), MPI_RealKind, &
                     recv_current%rank_comm, recv_current%tag, MPI_COMM_WORLD, recv_current%req, MPIErr)
      DO WHILE(ASSOCIATED(recv_current%next))
        recv_current => recv_current%next
        CALL MPI_IRECV(recv_current%fld_data, SIZE(recv_current%fld_data), MPI_RealKind, &
                       recv_current%rank_comm, recv_current%tag, MPI_COMM_WORLD, recv_current%req, MPIErr)
      END DO
    END IF
    IF (ASSOCIATED(subdomain_fine%tr_if%recv_bnds)) THEN
      send_current => subdomain_fine%tr_if%recv_bnds
      CALL MPI_WAIT(send_current%req, status, MPIErr)
      DO WHILE(ASSOCIATED(send_current%next))
        send_current => send_current%next
        CALL MPI_WAIT(send_current%req, status, MPIErr)
      END DO
    END IF
    IF (ASSOCIATED(subdomain_fine%tr_if%send_bnds)) THEN
      recv_current => subdomain_fine%tr_if%send_bnds
      CALL MPI_WAIT(recv_current%req, status, MPIErr)
      DO WHILE(ASSOCIATED(recv_current%next))
        recv_current => recv_current%next
        CALL MPI_WAIT(recv_current%req, status, MPIErr)
      END DO
    END IF

    DO i = 1, SIZE(Cut_offrank%ColInd)
      vols_eff_coarse_tmp(Cut_offrank%ColInd(i)) = 0.0
    END DO   

    CALL ATxpy_SpRowCol(Cut_offrank, subdomain_fine%tr_if%commbuff, vols_eff_coarse_tmp)
  END IF

  DO i = 1, SIZE(vols_eff_coarse_tmp)
    vols_eff_coarse_tmp(i) = 1.0 / vols_eff_coarse_tmp(i)
  END DO

  CALL SpAScalRow_SpRowCol(Restr, vols_eff_coarse_tmp)
  CALL SpAScalCol_SpRowCol(Restr, vols_eff_fine1d)
  
END SUBROUTINE restr_scaling


SUBROUTINE restr_offrank_indsto(js_blocks_to_unique, block_ids, blocks, jsblocks_to, blockbnds_to, nnz, indices, block_ptr)
  IMPLICIT NONE
  INTEGER, INTENT(in) :: js_blocks_to_unique(:)
  INTEGER, INTENT(in) :: block_ids(:)
  TYPE(Block), INTENT(in) :: blocks(:)
  TYPE(List), INTENT(in) :: jsblocks_to
  TYPE(List), INTENT(in) :: blockbnds_to
  INTEGER, INTENT(in) :: nnz
  INTEGER, POINTER, INTENT(inout) :: indices(:), block_ptr(:)

  INTEGER :: inds_unique_tmp(nnz)
  INTEGER :: iblock, iptr, index

  ALLOCATE(block_ptr(SIZE(js_blocks_to_unique, 1) + 1))

  iptr = 1
  block_ptr(1) = 1
  DO iblock = 1, SIZE(js_blocks_to_unique, 1)

    index = FINDLOC(block_ids, js_blocks_to_unique(iblock), 1)

    CALL restr_offrank_indsto_sb(js_blocks_to_unique(iblock), blocks(index), &
                                 jsblocks_to, blockbnds_to, inds_unique_tmp, iptr)
    block_ptr(iblock + 1) = iptr
  END DO

  ALLOCATE(indices(iptr - 1))
  indices = inds_unique_tmp(1:iptr - 1)

END SUBROUTINE restr_offrank_indsto


SUBROUTINE restr_offrank_indsto_sb(js_block, thisblock, jsblocks_to, blockbnds_to, inds_unique_tmp, iptr)
  IMPLICIT NONE

  INTEGER, INTENT(in) :: js_block
  TYPE(Block), INTENT(in) :: thisblock
  TYPE(List), INTENT(in) :: jsblocks_to
  TYPE(List), INTENT(in) :: blockbnds_to
  INTEGER, INTENT(inout) :: inds_unique_tmp(:)
  INTEGER, INTENT(inout) :: iptr
  INTEGER :: m, i, j, k
  TYPE(Element), POINTER :: elm => NULL()
  INTEGER, POINTER :: o21_ptr(:) => NULL()
  
  LOGICAL :: field3d(thisblock%fld_shape(1), thisblock%fld_shape(2), thisblock%fld_shape(3))

  field3d(:,:,:) = .FALSE.

  DO m = 1, jsblocks_to%len
    elm => jsblocks_to%get(m)
    IF (elm%ivalue .EQ. js_block) THEN
      elm => blockbnds_to%get(m)
      o21_ptr => elm%iarray1d
      DO i = o21_ptr(1), o21_ptr(2)
        DO j = o21_ptr(3), o21_ptr(4)
          DO k = o21_ptr(5), o21_ptr(6)
            field3d(k, j, i) = .True.
          END DO
        END DO
      END DO
    END IF
  END DO

  DO i = 1, thisblock%fld_shape(3)
    DO j = 1, thisblock%fld_shape(2)
      DO k = 1, thisblock%fld_shape(1)
        IF (field3d(k, j, i) .EQV. .True.) THEN
          inds_unique_tmp(iptr) = IndexC(i, j, k, thisblock%fld_shape(3), &
                                    thisblock%fld_shape(2), thisblock%fld_shape(1))
          iptr = iptr + 1
        END IF
      END DO
    END DO
  END DO

END SUBROUTINE restr_offrank_indsto_sb


SUBROUTINE find_bndblock_ij(blocks, bounding, x2_l, x2_u, y2_l, y2_u, block_id, bnd_block_id)
  IMPLICIT NONE
  
  TYPE(Block), TARGET, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  REAL(Realkind), INTENT(in) :: x2_l, x2_u, y2_l, y2_u
  INTEGER, INTENT(in) :: block_id
  INTEGER, INTENT(inout) :: bnd_block_id

  INTEGER :: iblock
  REAL(Realkind), POINTER :: x2_bnd(:), y2_bnd(:)

  bnd_block_id = HUGE(0)

  DO iblock = 1, SIZE(blocks, 1)
    
    IF (bounding(iblock) .EQV. .True.) THEN
     x2_bnd => blocks(block_id)%x2
     y2_bnd => blocks(block_id)%y2
 
     IF (ALL((/x2_u .GE. x2_bnd(1), x2_l .LE. x2_bnd(UBOUND(x2_bnd, 1)), &
               y2_u .GE. y2_bnd(1), y2_l .LE. y2_bnd(UBOUND(y2_bnd, 1))/))) THEN
       bnd_block_id = iblock
       RETURN 
     END IF
    END IF
  END DO

END SUBROUTINE find_bndblock_ij


SUBROUTINE find_bndblock_ik(blocks, bounding, x2_l, x2_u, z2_l, z2_u, block_id, bnd_block_id)
  IMPLICIT NONE

  TYPE(Block), TARGET, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  REAL(Realkind), INTENT(in) :: x2_l, x2_u, z2_l, z2_u
  INTEGER, INTENT(in) :: block_id
  INTEGER, INTENT(inout) :: bnd_block_id

  INTEGER :: iblock
  REAL(Realkind), POINTER :: x2_bnd(:), z2_bnd(:)

  bnd_block_id = HUGE(0)

  DO iblock = 1, SIZE(blocks, 1)

    IF (bounding(iblock) .EQV. .True.) THEN
     x2_bnd => blocks(block_id)%x2
     z2_bnd => blocks(block_id)%z2

     IF (ALL((/x2_u .GE. x2_bnd(1), x2_l .LE. x2_bnd(UBOUND(x2_bnd, 1)), &
               z2_u .GE. z2_bnd(1), z2_l .LE. z2_bnd(UBOUND(z2_bnd, 1))/))) THEN
       bnd_block_id = iblock
       RETURN
     END IF
    END IF
  END DO

END SUBROUTINE find_bndblock_ik


SUBROUTINE find_bndblock_jk(blocks, bounding, y2_l, y2_u, z2_l, z2_u, block_id, bnd_block_id)
  IMPLICIT NONE

  TYPE(Block), TARGET, INTENT(in) :: blocks(:)
  LOGICAL, INTENT(in) :: bounding(:)
  REAL(Realkind), INTENT(in) :: y2_l, y2_u, z2_l, z2_u
  INTEGER, INTENT(in) :: block_id
  INTEGER, INTENT(inout) :: bnd_block_id

  INTEGER :: iblock
  REAL(Realkind), POINTER :: y2_bnd(:), z2_bnd(:)

  bnd_block_id = HUGE(0)

  DO iblock = 1, SIZE(blocks, 1)

    IF (bounding(iblock) .EQV. .True.) THEN
     y2_bnd => blocks(block_id)%y2
     z2_bnd => blocks(block_id)%z2  

     IF (ALL((/y2_u .GE. y2_bnd(1), y2_l .LE. y2_bnd(UBOUND(y2_bnd, 1)), &
               z2_u .GE. z2_bnd(1), z2_l .LE. z2_bnd(UBOUND(z2_bnd, 1))/))) THEN
       bnd_block_id = iblock
       RETURN
     END IF
    END IF
  END DO

END SUBROUTINE find_bndblock_jk

END MODULE Transfer_Mod
