MODULE Laplace_block_Mod

USE MatSpRowCol_Mod, ONLY: SpRowCol, &
                           SpDeallocate_SpRowCol, SpNullify_SpRowCol, &
                           SpTrans_SpRowCol, SpMm_SpRowCol, resize_sprowcol_dataflds, &
                           SpMmBACT, SpMmBAC, SpMm_nnz, make_eye, SpMmBAC_tmp
USE Block_Mod, ONLY: Block
USE Crop_Mod, ONLY: make_crop_operator_cell, make_crop_operator_face, make_crop_operator_halo_T, &
                    make_crop_operator_boundary, make_crop_operator_halo_face_T, make_crop_operator_boundary_face
USE Kind_Mod
USE Index_Mod, ONLY: vector_size, &
                     IndexU3, IndexV3, IndexW3, IndexC3, &
                     IndexC3block, IndexU3block, IndexV3block, IndexW3block

USE Subdom_Mod, ONLY: Subdomain, BoundCommInterface

USE Gradient_block_Mod, ONLY: make_grad3d, ubound_ncolids_blockgrad

USE Gradient_block_interpolation_Mod, ONLY: find_indices_jk, find_indices_ik, find_indices_ij

USE Smooth_Mod, ONLY: make_red_black

IMPLICIT NONE

CONTAINS


SUBROUTINE make_operators_static(newsubdomain, btype, isblackijkzero)
  IMPLICIT NONE
 
  TYPE(SubDomain), TARGET, INTENT(inout) :: newsubdomain
  INTEGER, INTENT(in) :: btype(:), isblackijkzero(:)

  TYPE(Block), POINTER :: blocks(:) => NULL()
  TYPE(SpRowCol) :: Crop_domain_cell, Crop_domain_face
  TYPE(SpRowCol) :: grad3d_tmp, div3d_tmp, intc2f_tmp, Rd_tmp, Bl_tmp

  INTEGER, TARGET :: Crop_domain_cell_RowPtr(newsubdomain%ncells + 1), &
                     Crop_domain_cell_ColInd(newsubdomain%ncells)
  REAL(Realkind), TARGET :: Crop_domain_cell_Val(newsubdomain%ncells)

  INTEGER, TARGET :: Crop_domain_face_RowPtr(newsubdomain%nfaces + 1), &
                     Crop_domain_face_ColInd(newsubdomain%nfaces)
  REAL(Realkind), TARGET :: Crop_domain_face_Val(newsubdomain%nfaces)

  INTEGER, TARGET :: grad3d_tmp_RowPtr(newsubdomain%nfaces_tmp + 1), &
                     grad3d_tmp_ColInd(ubound_ncolids_blockgrad(newsubdomain%blocks)), &
                     div3d_tmp_RowPtr(newsubdomain%ncells_tmp + 1), &
                     div3d_tmp_ColInd(newsubdomain%ncells_tmp * 6), &
                     intc2f_tmp_RowPtr(newsubdomain%nfaces_tmp + 1), &
                     intc2f_tmp_ColInd(ubound_ncolids_blockgrad(newsubdomain%blocks)), &
                     rd_tmp_RowPtr(newsubdomain%ncells_tmp + 1), &
                     bl_tmp_RowPtr(newsubdomain%ncells_tmp + 1), &
                     rd_tmp_ColInd(newsubdomain%ncells_tmp), &
                     bl_tmp_ColInd(newsubdomain%ncells_tmp)

  REAL(Realkind), TARGET :: grad3d_tmp_Val(ubound_ncolids_blockgrad(newsubdomain%blocks)), &
                            div3d_tmp_Val(newsubdomain%ncells_tmp * 6), &
                            intc2f_tmp_Val(ubound_ncolids_blockgrad(newsubdomain%blocks)), &
                            rdbl_tmp_Val(newsubdomain%ncells_tmp)

  REAL(Realkind) :: cellisopen(newsubdomain%ncells)

  INTEGER :: nnz_lapl3d_tmp, n, i, j, k, icell

  blocks => newsubdomain%blocks

  Crop_domain_cell%RowPtr => Crop_domain_cell_RowPtr
  Crop_domain_cell%ColInd => Crop_domain_cell_ColInd
  Crop_domain_cell%Val => Crop_domain_cell_Val

  Crop_domain_face%RowPtr => Crop_domain_face_RowPtr
  Crop_domain_face%ColInd => Crop_domain_face_ColInd
  Crop_domain_face%Val => Crop_domain_face_Val             

  CALL make_crop_operator_cell(newsubdomain, Crop_domain_cell, allocate_mat=.FALSE.)
  CALL make_crop_operator_face(newsubdomain, Crop_domain_face, allocate_mat=.FALSE.) 

  grad3d_tmp%m = newsubdomain%nfaces_tmp
  grad3d_tmp%n = newsubdomain%ncells_tmp
  grad3d_tmp%RowPtr => grad3d_tmp_RowPtr
  grad3d_tmp%ColInd => grad3d_tmp_ColInd
  grad3d_tmp%Val => grad3d_tmp_Val
  CALL make_grad3d(grad3d_tmp, blocks)

  CALL SpMmBACT(newsubdomain%grad3d, Crop_domain_face, grad3d_tmp, Crop_domain_cell)

  div3d_tmp%m = newsubdomain%ncells_tmp
  div3d_tmp%n = newsubdomain%nfaces_tmp
  div3d_tmp%RowPtr => div3d_tmp_RowPtr
  div3d_tmp%ColInd => div3d_tmp_ColInd
  div3d_tmp%Val => div3d_tmp_Val
  CALL make_div3d(div3d_tmp, newsubdomain, btype)

  CALL SpMmBACT(newsubdomain%div3d, Crop_domain_cell, div3d_tmp, Crop_domain_face)

  nnz_lapl3d_tmp = SpMm_nnz(div3d_tmp, grad3d_tmp)
  CALL make_lapl3d(newsubdomain, div3d_tmp, grad3d_tmp, Crop_domain_cell, nnz_lapl3d_tmp)

  intc2f_tmp%m = newsubdomain%nfaces_tmp
  intc2f_tmp%n = newsubdomain%ncells_tmp
  intc2f_tmp%RowPtr => intc2f_tmp_RowPtr
  intc2f_tmp%Val => intc2f_tmp_Val
  intc2f_tmp%ColInd => intc2f_tmp_ColInd
  CALL make_intcell2face3d(intc2f_tmp, newsubdomain)
  
  CALL SpMmBACT(newsubdomain%intc2f, Crop_domain_face, intc2f_tmp, Crop_domain_cell)

  icell = 1
  DO n = 1, newsubdomain%nblocks
    IF (newsubdomain%blockiscomp(n)) THEN
      DO i = 1, newsubdomain%blocks(n)%fld_shape(3)
        DO j = 1, newsubdomain%blocks(n)%fld_shape(2)
          DO k = 1, newsubdomain%blocks(n)%fld_shape(1)
            IF (newsubdomain%blocks(n)%volseff(k, j, i) .LT. 1.0001e-20) THEN
              cellisopen(icell) = 0.0
            ELSE
              cellisopen(icell) = 1.0
            END IF
            icell = icell + 1
          END DO
        END DO
      END DO
    END IF
  END DO

  CALL make_eye(newsubdomain%eye, cellisopen)

  rdbl_tmp_Val(:) = 1.0
  Bl_tmp%RowPtr => bl_tmp_RowPtr
  Rd_tmp%RowPtr => rd_tmp_RowPtr
  Bl_tmp%ColInd => bl_tmp_ColInd
  Rd_tmp%ColInd => rd_tmp_ColInd
  Bl_tmp%Val => rdbl_tmp_Val
  Rd_tmp%Val => rdbl_tmp_Val

  CALL make_red_black(newsubdomain%blocks, newsubdomain%ncells_tmp, newsubdomain%blockids_compghst, isblackijkzero, Bl_tmp, Rd_tmp)

  CALL SpMmBACT(newsubdomain%Rd, Crop_domain_cell, Rd_tmp, Crop_domain_cell)
  CALL SpMmBACT(newsubdomain%Bl, Crop_domain_cell, Bl_tmp, Crop_domain_cell)

  IF (ASSOCIATED(newsubdomain%bc_if)) THEN
    CALL make_operators_commIf(newsubdomain, newsubdomain%bc_if, div3d_tmp, grad3d_tmp, intc2f_tmp, Bl_tmp, Rd_tmp)
  END IF

!  CALL SpDeallocate_SpRowCol(newsubdomain%div3d)
!  CALL SpNullify_SpRowCol(newsubdomain%div3d)
!  CALL SpDeallocate_SpRowCol(newsubdomain%grad3d)
!  CALL SpNullify_SpRowCol(newsubdomain%grad3d)
!  CALL SpDeallocate_SpRowCol(newsubdomain%lapl3d)
!  CALL SpNullify_SpRowCol(newsubdomain%lapl3d)

  CALL SpNullify_SpRowCol(Crop_domain_cell)
  CALL SpNullify_SpRowCol(Crop_domain_face)
  CALL SpNullify_SpRowCol(grad3d_tmp)
  CALL SpNullify_SpRowCol(div3d_tmp)
  CALL SpNullify_SpRowCol(intc2f_tmp)
  CALL SpNullify_SpRowCol(Bl_tmp)
  CALL SpNullify_SpRowCol(Rd_tmp)

END SUBROUTINE make_operators_static


SUBROUTINE make_operators_commIf(newsubdomain, comm_interface, div3d_tmp, grad3d_tmp, intc2f_tmp, Bl_tmp, Rd_tmp)
  IMPLICIT NONE
  TYPE(SubDomain), INTENT(in) :: newsubdomain
  TYPE(BoundCommInterface), INTENT(inout) :: comm_interface
  TYPE(SpRowCol), INTENT(in) :: div3d_tmp, grad3d_tmp, intc2f_tmp, Bl_tmp, Rd_tmp

  TYPE(SpRowCol) :: Crop_bound, Crop_halo_T, Crop_bound_face, Crop_halo_face_T, eye_tmp

  INTEGER, TARGET :: Crop_halo_T_RowPtr(newsubdomain%ncells_tmp + 1), &
                     Crop_halo_T_ColInd(newsubdomain%nghost), &
                     Crop_bound_RowPtr(newsubdomain%nbound + 1), &
                     Crop_bound_ColInd(newsubdomain%nbound), &
                     Crop_halo_face_T_RowPtr(newsubdomain%nfaces_tmp + 1), &
                     Crop_halo_face_T_ColInd(SUM(newsubdomain%indsface1d_type)), &
                     Crop_bound_face_RowPtr(SUM(newsubdomain%indsface1d_type) + 1), &
                     Crop_bound_face_ColInd(SUM(newsubdomain%indsface1d_type)), &
                     eye_tmp_RowPtr(div3d_tmp%m + 1), &
                     eye_tmp_ColInd(div3d_tmp%m)
                     

  REAL(Realkind), TARGET :: Crop_halo_T_Val(newsubdomain%nghost), &
                            Crop_bound_Val(newsubdomain%nbound), &
                            Crop_halo_face_T_Val(SUM(newsubdomain%indsface1d_type)), &
                            Crop_bound_face_Val(SUM(newsubdomain%indsface1d_type)), &
                            eye_tmp_Val(div3d_tmp%m)

  INTEGER :: nnz

  Crop_halo_T%RowPtr => Crop_halo_T_RowPtr
  Crop_halo_T%ColInd => Crop_halo_T_ColInd
  Crop_halo_T%Val => Crop_halo_T_Val

  Crop_bound%RowPtr => Crop_bound_RowPtr
  Crop_bound%ColInd => Crop_bound_ColInd
  Crop_bound%Val => Crop_bound_Val

  Crop_halo_face_T%RowPtr => Crop_halo_face_T_RowPtr
  Crop_halo_face_T%ColInd => Crop_halo_face_T_ColInd
  Crop_halo_face_T%Val => Crop_halo_face_T_Val

  Crop_bound_face%RowPtr => Crop_bound_face_RowPtr
  Crop_bound_face%ColInd => Crop_bound_face_ColInd
  Crop_bound_face%Val => Crop_bound_face_Val

  CALL make_crop_operator_halo_T(newsubdomain, Crop_halo_T)
  CALL make_crop_operator_boundary(newsubdomain, Crop_bound)
  CALL make_crop_operator_halo_face_T(newsubdomain, Crop_halo_face_T)
  CALL make_crop_operator_boundary_face(newsubdomain, Crop_bound_face)

  nnz = SpMm_nnz(div3d_tmp, Crop_halo_face_T)
  CALL SpMmBAC(comm_interface%halo%Div3d_halo, Crop_bound, div3d_tmp, Crop_halo_face_T, nnz)
  nnz = SpMm_nnz(grad3d_tmp, Crop_halo_T)
  CALL SpMmBAC(comm_interface%halo%Grad3d_halo, Crop_bound_face, grad3d_tmp, Crop_halo_T, nnz)
  nnz = SpMm_nnz(intc2f_tmp, Crop_halo_T)
  CALL SpMmBAC(comm_interface%halo%Intc2f_halo, Crop_bound_face, intc2f_tmp, Crop_halo_T, nnz)

  ALLOCATE(comm_interface%halo%Lapl3d_halo)
  CALL SpMm_SpRowCol(comm_interface%halo%Lapl3d_halo, comm_interface%halo%Div3d_halo, comm_interface%halo%Grad3d_halo)

!  comm_interface%halo%mat_halo => comm_interface%halo%Lapl3d_halo

  CALL SpMmBACT(comm_interface%halo%Rd_halo, Crop_bound, Rd_tmp, Crop_bound)
  CALL SpMmBACT(comm_interface%halo%Bl_halo, Crop_bound, Bl_tmp, Crop_bound)

!  CALL SpDeallocate_SpRowCol(comm_interface%halo%Eye)
!  CALL SpDeallocate_SpRowCol(comm_interface%halo%Div3d_halo)
!  CALL SpDeallocate_SpRowCol(comm_interface%halo%Grad3d_halo)
!  CALL SpDeallocate_SpRowCol(comm_interface%halo%Lapl3d_halo)
  
!  CALL SpNullify_SpRowCol(comm_interface%halo%Eye)
!  CALL SpNullify_SpRowCol(comm_interface%halo%Div3d_halo)
!  CALL SpNullify_SpRowCol(comm_interface%halo%Grad3d_halo)
!  CALL SpNullify_SpRowCol(comm_interface%halo%Lapl3d_halo)
  CALL SpNullify_SpRowCol(eye_tmp)
  CALL SpNullify_SpRowCol(Crop_halo_T)
  CALL SpNullify_SpRowCol(Crop_halo_face_T)
  CALL SpNullify_SpRowCol(Crop_bound)
  CALL SpNullify_SpRowCol(Crop_bound_face)

END SUBROUTINE make_operators_commIf


SUBROUTINE make_lapl3d(newsubdomain, div3d_tmp, grad3d_tmp, Crop_domain_cell, nnz_lapl_tmp)
  IMPLICIT NONE
  
  TYPE(SubDomain), INTENT(inout) :: newsubdomain
  TYPE(SpRowCol), INTENT(in) :: div3d_tmp, grad3d_tmp, Crop_domain_cell
  INTEGER, INTENT(in) :: nnz_lapl_tmp

  TYPE(SpRowCol) :: lapl3d_tmp

  INTEGER, TARGET :: lapl3d_tmp_RowPtr(div3d_tmp%m + 1), &
                     lapl3d_tmp_ColInd(nnz_lapl_tmp)
  REAL(Realkind), TARGET :: lapl3d_tmp_val(nnz_lapl_tmp)
  INTEGER :: nnz


  lapl3d_tmp%RowPtr => lapl3d_tmp_RowPtr
  lapl3d_tmp%ColInd => lapl3d_tmp_ColInd
  lapl3d_tmp%Val => lapl3d_tmp_val
  
  CALL SpMm_SpRowCol(lapl3d_tmp, div3d_tmp, grad3d_tmp)

  CALL SpMmBACT(newsubdomain%lapl3d, Crop_domain_cell, lapl3d_tmp, Crop_domain_cell)

  CALL SpNullify_SpRowCol(lapl3d_tmp)

END SUBROUTINE make_lapl3d


SUBROUTINE make_div3d_stencil(DIV3d_st, blocks)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(INOUT) :: DIV3d_st
  TYPE(Block), POINTER, INTENT(in) :: blocks(:)

  INTEGER :: iblock, nblocks
  INTEGER :: ist_face, irow

  INTEGER :: k, nz
  INTEGER :: j, ny
  INTEGER :: i, nx

  INTEGER :: ii

  nblocks = SIZE(blocks)

  DIV3d_st%m = 0
  DIV3d_st%n = 0

  DO iblock = 1, nblocks
    DIV3d_st%m = DIV3d_st%m + blocks(iblock)%ncells
    DIV3d_st%n = DIV3d_st%n + blocks(iblock)%nfaces
  END DO  

  ALLOCATE(DIV3d_st%RowPtr(DIV3d_st%m + 1))
  ALLOCATE(DIV3d_st%ColInd(DIV3d_st%m * 6))
  ALLOCATE(DIV3d_st%Val(DIV3d_st%m * 6))

  ii = 1
  irow = 1

  ist_face = 0

  DIV3d_st%RowPtr(1) = 1

  DO iblock = 1, nblocks

    nz = blocks(iblock)%fld_shape(1)
    ny = blocks(iblock)%fld_shape(2)
    nx = blocks(iblock)%fld_shape(3)

    DO i = 1, nx
      DO j = 1, ny
        DO k = 1, nz 

          DIV3d_st%ColInd(ii) = IndexU3(i, j, k, nx, ny, nz) + ist_face
          DIV3d_st%Val(ii) = -1.0
          ii = ii + 1
          DIV3d_st%ColInd(ii) = IndexU3(i + 1, j, k, nx, ny, nz) + ist_face
          DIV3d_st%Val(ii) = 1.0
          ii = ii + 1
          DIV3d_st%ColInd(ii) = IndexV3(i, j, k, nx, ny, nz) + ist_face
          DIV3d_st%Val(ii) = -1.0
          ii = ii + 1
          DIV3d_st%ColInd(ii) = IndexV3(i, j + 1, k, nx, ny, nz) + ist_face
          DIV3d_st%Val(ii) = 1.0
          ii = ii + 1
          DIV3d_st%ColInd(ii) = IndexW3(i, j, k, nx, ny, nz) + ist_face
          DIV3d_st%Val(ii) = -1.0
          ii = ii + 1
          DIV3d_st%ColInd(ii) = IndexW3(i, j, k + 1, nx, ny, nz) + ist_face
          DIV3d_st%Val(ii) = 1.0
          ii = ii + 1
          DIV3d_st%RowPtr(irow + 1) = ii
          irow = irow + 1

        END DO
      END DO
    END DO

    ist_face = ist_face + blocks(iblock)%nfaces

  END DO

END SUBROUTINE make_div3d_stencil


SUBROUTINE make_div3d(DIV3d, subdom, btype)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: DIV3d
  TYPE(Subdomain), INTENT(in) :: subdom
  INTEGER, INTENT(in) :: btype(:)

  INTEGER :: iblock, nblocks
  INTEGER :: ist_face, irow

  INTEGER :: k, nz
  INTEGER :: j, ny
  INTEGER :: i, nx
  INTEGER :: iu1, iu2, iv1, iv2, iw1, iw2, ic

  INTEGER :: ii

!  INTEGER :: icol

!  icol = 1

  nblocks = subdom%nblocks

!  DIV3d%m = 0
!  DIV3d%n = 0


!  DO iblock = 1, nblocks
!    IF (btype(iblock) .EQ. 1) THEN
!      icol = icol + subdom%blocks(iblock)%ncells * 6
!    END IF

!    DIV3d%m = DIV3d%m + subdom%blocks(iblock)%ncells
!    DIV3d%n = DIV3d%n + subdom%blocks(iblock)%nfaces
!  END DO

!  ALLOCATE(DIV3d%RowPtr(DIV3d%m + 1))
!  ALLOCATE(DIV3d%ColInd(icol))
!  ALLOCATE(DIV3d%Val(icol))

  ii = 1
  irow = 1

  ist_face = 0

  DIV3d%RowPtr(1) = 1

  DO iblock = 1, nblocks

    nz = subdom%blocks(iblock)%fld_shape(1)
    ny = subdom%blocks(iblock)%fld_shape(2)
    nx = subdom%blocks(iblock)%fld_shape(3)

    IF (btype(subdom%blockids_compghst(iblock)) .EQ. 1) THEN
      DO i = 1, nx
        DO j = 1, ny
          DO k = 1, nz

            iu1 = IndexU3(i, j, k, nx, ny, nz)
            iu2 = IndexU3(i + 1, j, k, nx, ny, nz)
            iv1 = IndexV3(i, j, k, nx, ny, nz)
            iv2 = IndexV3(i, j + 1, k, nx, ny, nz)
            iw1 = IndexW3(i, j, k, nx, ny, nz)
            iw2 = IndexW3(i, j, k + 1, nx, ny, nz)

            DIV3d%ColInd(ii) = iu1 + ist_face          
            DIV3d%Val(ii) = -subdom%blocks(iblock)%arseffx(k, j, i) / subdom%blocks(iblock)%volseff(k, j, i)
            ii = ii + 1
            DIV3d%ColInd(ii) = iu2 + ist_face
            DIV3d%Val(ii) = subdom%blocks(iblock)%arseffx(k, j, i + 1) / subdom%blocks(iblock)%volseff(k, j, i)
            ii = ii + 1
            DIV3d%ColInd(ii) = iv1 + ist_face
            DIV3d%Val(ii) = -subdom%blocks(iblock)%arseffy(k, j, i) / subdom%blocks(iblock)%volseff(k, j, i)
            ii = ii + 1
            DIV3d%ColInd(ii) = iv2 + ist_face
            DIV3d%Val(ii) = subdom%blocks(iblock)%arseffy(k, j + 1, i) / subdom%blocks(iblock)%volseff(k, j, i)
            ii = ii + 1
            DIV3d%ColInd(ii) = iw1 + ist_face
            DIV3d%Val(ii) = -subdom%blocks(iblock)%arseffz(k, j, i) / subdom%blocks(iblock)%volseff(k, j, i)
            ii = ii + 1
            DIV3d%ColInd(ii) = iw2 + ist_face
            DIV3d%Val(ii) = subdom%blocks(iblock)%arseffz(k + 1, j, i) / subdom%blocks(iblock)%volseff(k, j, i)          
            ii = ii + 1
            DIV3d%RowPtr(irow + 1) = ii
            irow = irow + 1

          END DO
        END DO
      END DO
    ELSE
      DO i = 1, nx
        DO j = 1, ny
          DO k = 1, nz
            DIV3d%RowPtr(irow + 1) = ii
            irow = irow + 1
          END DO
        END DO
      END DO
    END IF
    ist_face = ist_face + subdom%blocks(iblock)%nfaces
  END DO

END SUBROUTINE make_div3d


SUBROUTINE make_intcell2face3d(Intc2f, subdom)

  IMPLICIT NONE

  TYPE(SpRowCol), INTENT(inout) :: Intc2f
  TYPE(SubDomain), POINTER, INTENT(in) :: subdom
  TYPE(Block), POINTER :: blocks(:)

  INTEGER :: iblock, nblocks, bndblock_id

  INTEGER :: ist_cell, irow

  INTEGER :: ucface_ncolids

  INTEGER :: k, nz
  INTEGER :: j, ny
  INTEGER :: i, nx
  integer :: n, n_tmp

  INTEGER :: ii, jj, kk

  INTEGER :: i_bnds(2), j_bnds(2), k_bnds(2)

  REAL(Realkind), PARAMETER :: eps = 1e-20
  REAL(Realkind), POINTER :: arsx(:,:,:), arsy(:,:,:), arsz(:,:,:)
  REAL(Realkind), POINTER :: arseffx(:,:,:), arseffy(:,:,:), arseffz(:,:,:)
  REAL(Realkind), POINTER :: arsx_bnd(:,:,:), arsy_bnd(:,:,:), arsz_bnd(:,:,:)
  REAL(Realkind), POINTER :: arseffx_bnd(:,:,:), arseffy_bnd(:,:,:), arseffz_bnd(:,:,:)
  REAL(Realkind) :: area_fac, icoeff, icoeff_sum

  nblocks = subdom%nblocks
  blocks => subdom%blocks

  n = 1
  irow = 1

  ist_cell = 0

  Intc2f%RowPtr(1) = 1

  DO iblock = 1, nblocks
   
    nz = blocks(iblock)%fld_shape(1)
    ny = blocks(iblock)%fld_shape(2)
    nx = blocks(iblock)%fld_shape(3)
    
    arsx => blocks(iblock)%arsx
    arsy => blocks(iblock)%arsy
    arsz => blocks(iblock)%arsz

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
                Intc2f%RowPtr(irow + 1) = n
                irow = irow + 1
                CYCLE
              END IF

              arsx_bnd => blocks(bndblock_id)%arsx
              arseffx_bnd => blocks(bndblock_id)%arseffx
              icoeff_sum = 0.0

              DO jj = j_bnds(1), j_bnds(2)                
                DO kk = k_bnds(1), k_bnds(2)
                  area_fac = arsx_bnd(kk, jj, UBOUND(arsx_bnd, 3)) / (arsx(k, j, i) + eps)
                  IF (area_fac .GT. 1.0) THEN
                    icoeff =  0.5
                  ELSE
                    icoeff = 0.5 * arseffx_bnd(kk, jj, UBOUND(arseffx_bnd, 3)) / (arseffx(k, j, i) + eps)
                  END IF
                  icoeff_sum = icoeff_sum + icoeff
                  Intc2f%ColInd(n) = IndexC3block(SIZE(blocks(bndblock_id)%x2) - 1, jj, kk, blocks, bndblock_id, check=.TRUE.)
                  Intc2f%Val(n) = icoeff
                  n = n + 1
                END DO
              END DO
              Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
              Intc2f%Val(n) = icoeff_sum
              n = n + 1
              Intc2f%RowPtr(irow + 1) = n
              irow = irow + 1
            END DO
          END DO
        ELSE
          DO j = 1, ny
            DO k = 1, nz
              Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
              Intc2f%Val(n) = arseffx(k, j, i) / arsx(k, j, i)
              n = n + 1
              Intc2f%RowPtr(irow + 1) = n
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
                Intc2f%RowPtr(irow + 1) = n
                irow = irow + 1
                CYCLE
              END IF

              arsx_bnd => blocks(bndblock_id)%arsx
              arseffx_bnd => blocks(bndblock_id)%arseffx
              icoeff_sum = 0.0
              Intc2f%ColInd(n) = IndexC3block(i - 1, j, k, blocks, iblock, check=.TRUE.)
              n_tmp = n
              n = n + 1
              DO jj = j_bnds(1), j_bnds(2)
                DO kk = k_bnds(1), k_bnds(2)
                  area_fac = arsx_bnd(kk, jj, 1) / (arsx(k, j, i) + eps)
                  IF (area_fac .GT. 1.0) THEN
                    icoeff = 0.5
                  ELSE
                    icoeff = 0.5 * arseffx_bnd(kk, jj, 1) / (arseffx(k, j, i) + eps)
                  END IF
                  icoeff_sum = icoeff_sum + icoeff
                  Intc2f%ColInd(n) = IndexC3block(1, jj, kk, blocks, bndblock_id, check=.FALSE.)
                  Intc2f%Val(n) = icoeff
                  n = n + 1
                END DO
              END DO
              Intc2f%Val(n_tmp) = icoeff_sum
              Intc2f%RowPtr(irow + 1) = n
              irow = irow + 1
            END DO
          END DO
        ELSE
          DO j = 1, ny
            DO k = 1, nz
              Intc2f%ColInd(n) = IndexC3block(i - 1, j, k, blocks, iblock, check=.TRUE.)
              Intc2f%Val(n) = arseffx(k, j, i) / arsx(k, j, i)
              n = n + 1
              Intc2f%RowPtr(irow + 1) = n
              irow = irow + 1
            END DO
          END DO
        END IF
      ELSE
        DO j = 1, ny
          DO k = 1, nz
            Intc2f%ColInd(n) = IndexC3block(i - 1, j, k, blocks, iblock, check=.FALSE.)
            Intc2f%Val(n) = 0.5
            n = n + 1
            Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.FALSE.)
            Intc2f%Val(n) = 0.5
            n = n + 1
            Intc2f%RowPtr(irow + 1) = n
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
                Intc2f%RowPtr(irow + 1) = n
                irow = irow + 1
                CYCLE
              END IF

              arsy_bnd => blocks(bndblock_id)%arsy
              arseffy_bnd => blocks(bndblock_id)%arseffy
              icoeff_sum = 0.0

              DO ii = i_bnds(1), i_bnds(2)
                DO kk = k_bnds(1), k_bnds(2)
                  area_fac = arsy_bnd(kk, UBOUND(arsy_bnd, 2), ii) / (arsy(k, j, i) + eps)
                  IF (area_fac .GT. 1.0) THEN
                    icoeff = 0.5
                  ELSE
                    icoeff =  0.5 * arseffy_bnd(kk, UBOUND(arseffy_bnd, 2), ii) / (arseffy(k, j, i) + eps)
                  END IF
                  icoeff_sum = icoeff_sum + icoeff

                  Intc2f%ColInd(n) = IndexC3block(ii, SIZE(blocks(bndblock_id)%y2) - 1, kk, blocks, bndblock_id, check=.TRUE.)
                  Intc2f%Val(n) = icoeff
                  n = n + 1
                END DO
             END DO
              Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
              Intc2f%Val(n) = icoeff_sum
              n = n + 1
              Intc2f%RowPtr(irow + 1) = n
              irow = irow + 1
            END DO
          ELSE
            DO k = 1, nz
              Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
              Intc2f%Val(n) = arseffy(k, j, i) / arsy(k, j, i)
              n = n + 1
              Intc2f%RowPtr(irow + 1) = n
              irow = irow + 1
            END DO
          END IF
        ELSEIF (j == ny + 1) THEN
          IF (ANY(blocks(iblock)%cface_n)) THEN
            DO k = 1, nz
              CALL find_indices_ik(blocks, blocks(iblock)%cface_n, i, i_bnds, k, k_bnds, iblock, bndblock_id)
              IF (bndblock_id .GT. 99999999) THEN
                Intc2f%RowPtr(irow + 1) = n
                irow = irow + 1
                CYCLE
              END IF

              Intc2f%ColInd(n) = IndexC3block(i, j - 1, k, blocks, iblock, check=.TRUE.)
              arsy_bnd => blocks(bndblock_id)%arsy
              arseffy_bnd => blocks(bndblock_id)%arseffy
              icoeff_sum = 0.0
              n_tmp = n
              n = n + 1
            
              DO ii = i_bnds(1), i_bnds(2)
                DO kk = k_bnds(1), k_bnds(2)
                  area_fac = arsy_bnd(kk, 1, ii) / (arsy(k, j, i) + eps)
                  IF (area_fac .GT. 1.0) THEN
                    icoeff = 0.5
                  ELSE
                    icoeff = 0.5 * arseffy_bnd(kk, 1, ii) / (arseffy(k, j, i) + eps)
                  END IF
                  icoeff_sum = icoeff_sum + icoeff
                  Intc2f%ColInd(n) = IndexC3block(ii, 1, kk, blocks, bndblock_id, check=.TRUE.)
                  Intc2f%Val(n) = icoeff
                  n = n + 1
                END DO
              END DO
              Intc2f%Val(n_tmp) = icoeff_sum
              Intc2f%RowPtr(irow + 1) = n
              irow = irow + 1
            END DO
          ELSE
            DO k = 1, nz
              Intc2f%ColInd(n) = IndexC3block(i, j - 1, k, blocks, iblock, check=.TRUE.)
              Intc2f%Val(n) = arseffy(k, j, i) / arsy(k, j, i)
              n = n + 1
              Intc2f%RowPtr(irow + 1) = n
              irow = irow + 1
            END DO
          END IF
        ELSE
          DO k = 1, nz
            Intc2f%ColInd(n) = IndexC3block(i, j - 1, k, blocks, iblock, check=.TRUE.)
            Intc2f%Val(n) = 0.5
            n = n + 1
            Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock, check=.TRUE.)
            Intc2f%Val(n) = 0.5
            n = n + 1
            Intc2f%RowPtr(irow + 1) = n
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
                Intc2f%RowPtr(irow + 1) = n
                irow = irow + 1
                CYCLE
              END IF

              arsz_bnd => blocks(bndblock_id)%arsz
              arseffz_bnd => blocks(bndblock_id)%arseffz
              icoeff_sum = 0.0

              DO ii = i_bnds(1), i_bnds(2)
                DO jj = j_bnds(1), j_bnds(2)
                  area_fac = arsz_bnd(UBOUND(arsz_bnd, 1), jj, ii) / (arsz(k, j, i) + eps)
                  IF (area_fac .GT. 1.0) THEN
                    icoeff = 0.5
                  ELSE
                    icoeff = 0.5 * arseffz_bnd(UBOUND(arseffz_bnd, 1), jj, ii) / (arseffz(k, j, i) + eps)
                  END IF
                  Intc2f%ColInd(n) = IndexC3block(ii, jj, SIZE(blocks(bndblock_id)%z2) - 1, blocks, bndblock_id, check=.TRUE.)
                  Intc2f%Val(n) = icoeff
                  n = n + 1
                  icoeff_sum = icoeff_sum + icoeff
                END DO
              END DO
              Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock)
              Intc2f%Val(n) = icoeff_sum
              n = n + 1
            ELSE
              Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock)
              Intc2f%Val(n) = arseffz(k, j, i) / arsz(k, j, i)
              n = n + 1
            END IF
            Intc2f%RowPtr(irow + 1) = n
            irow = irow + 1
          ELSE IF (k == nz + 1) THEN   
            IF (ANY(blocks(iblock)%cface_t)) THEN
              CALL find_indices_ij(blocks, blocks(iblock)%cface_t, i, i_bnds, j, j_bnds, iblock, bndblock_id)
              IF (bndblock_id .GT. 99999999) THEN
                Intc2f%RowPtr(irow + 1) = n
                irow = irow + 1
                CYCLE
              END IF

              arsz_bnd => blocks(bndblock_id)%arsz
              arseffz_bnd => blocks(bndblock_id)%arseffz
              Intc2f%ColInd(n) = IndexC3block(i, j, k - 1, blocks, iblock)
              n_tmp = n 
              n = n + 1
              icoeff_sum = 0.0

              DO ii = i_bnds(1), i_bnds(2)
                DO jj = j_bnds(1), j_bnds(2)
                  area_fac = arsz_bnd(1, jj, ii) / (arsz(k, j, i) + eps)
                  IF (area_fac .GT. 1.0) THEN
                    icoeff = 0.5
                  ELSE
                    icoeff = 0.5 * arseffz_bnd(1, jj, ii) / (arseffz(k, j, i) + eps)
                  END IF
                  icoeff_sum = icoeff_sum + icoeff
                  Intc2f%ColInd(n) = IndexC3block(ii, jj, 1, blocks, bndblock_id, check=.TRUE.)
                  Intc2f%Val(n) = icoeff
                  n = n + 1
                END DO
              END DO 
              Intc2f%Val(n_tmp) = icoeff_sum
            ELSE
              Intc2f%ColInd(n) = IndexC3block(i, j, k - 1, blocks, iblock)
              Intc2f%Val(n) = arseffz(k, j, i) / arsz(k, j, i)
              n = n + 1
            END IF 
            Intc2f%RowPtr(irow + 1) = n
            irow = irow + 1
          ELSE
            Intc2f%ColInd(n) = IndexC3block(i, j, k - 1, blocks, iblock)
            Intc2f%Val(n) = 0.5
            n = n + 1
            Intc2f%ColInd(n) = IndexC3block(i, j, k, blocks, iblock)
            Intc2f%Val(n) = 0.5
            n = n + 1
            Intc2f%RowPtr(irow + 1) = n
            irow = irow + 1
          END IF
        END DO
      END DO
    END DO
  END DO

END SUBROUTINE make_intcell2face3d




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



END MODULE Laplace_block_Mod
