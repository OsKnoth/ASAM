MODULE Boundary_Mod

USE Kind_Mod
USE Block_Mod, ONLY: Block
USE Index_Mod, ONLY: IndexC3, IndexC3block, IndexU3block,IndexV3block, IndexW3block, &
                     findloc_real
USE MatSpRowCol_Mod, ONLY: SpRowCol, SpAVec_SpRowCol, SpMm_SpRowCol, SpDeallocate_SpRowCol, &
                           SpNullify_SpRowCol, SpTrans_SpRowCol, make_eye, Axpy_SpRowCol

USE Parallel_Mod

!Boundary structure for point to point communication of field parts
TYPE Bound
  INTEGER :: rank_comm ! process rank to communicate with
  INTEGER, POINTER :: fld_inds(:) => NULL() ! vector of field indices to communicate
  REAL(Realkind), POINTER :: fld_data(:) => NULL() ! buffer to send / receive
  INTEGER, POINTER :: blockfld_inds(:) => NULL()
  INTEGER, POINTER :: blockids(:) => NULL()
  INTEGER, POINTER :: blockptr(:) => NULL()
  INTEGER :: tag               ! tag
  TYPE(MPI_Request) :: req
  TYPE(Bound), POINTER :: next => NULL()
END TYPE Bound


!Halo structure that includes all ghost cells and gathers all received on halo_data
!The boundary condition is computed with the halo Laplace operator and stored on
!fld_data, which is finally added to the computation field
TYPE HaloType
  REAL(Realkind), POINTER :: hdata(:) => NULL()      !halo data array
  TYPE(SpRowCol), POINTER :: Lapl3d_halo => NULL()   !Laplace operator of shape (SIZE(fld_data) x SIZE(halo_data))
                                                     !used to compute the boundary condition from the halo data
  TYPE(SpRowCol), POINTER :: mat_halo => NULL()

  TYPE(SpRowCol) :: Grad3d_halo, Gradp3d_halo, Div3d_halo, Eye, Intc2f_halo, Rd_halo, Bl_halo
  TYPE(SpRowCol), POINTER :: OP_crop_bound => NULL()
  TYPE(SpRowCol), POINTER :: OP_crop_halo_T => NULL()
  TYPE(SpRowCol), POINTER :: mat_rd_halo => NULL()
  TYPE(SpRowCol), POINTER :: mat_bl_halo => NULL()

  REAL(Realkind), POINTER :: fdata(:) => NULL()   !computed boundary condition to add to the boundary layer of the inner computation field
  INTEGER, POINTER :: finds(:) => NULL()          !computation field indices of boundary layers where fld_data is added to
  REAL(Realkind), POINTER :: gdata(:) => NULL()   !equivalent data field for cell faces to compute gradient bc
  INTEGER, POINTER :: ginds(:) => NULL()          !equivalent field indices for cell faces to compute gradient bc

  INTEGER :: myrank
END TYPE HaloType

!Communication interface consisting of a halo structure and boundary objects
TYPE BoundCommInterface
  INTEGER :: myrank
  TYPE(HaloType) :: halo
  TYPE(Bound), POINTER :: send_bnds => NULL()
  TYPE(Bound), POINTER :: recv_bnds => NULL()
  INTEGER, POINTER :: block2halo_inds(:) => NULL()
END TYPE BoundCommInterface

LOGICAL :: cyclic_x, cyclic_y, cyclic_z


CONTAINS

!--------------------------------------------------------------
!Routines to exchange and compute the boundary condition 
!for smoothing

SUBROUTINE matmul_boundary_condition(comm_interf, ax, mode)

  IMPLICIT NONE
  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf
  REAL(Realkind), INTENT(inout) :: ax(:)
  CHARACTER(len=3), OPTIONAL, INTENT(in) :: mode


  CALL SpAVec_SpRowCol(comm_interf%halo%fdata, comm_interf%halo%mat_halo, comm_interf%halo%hdata)

  IF (PRESENT(mode)) THEN
    CALL copy_data_from_halo_cell(comm_interf%halo, ax, mode='sub')
  ELSE
    CALL copy_data_from_halo_cell(comm_interf%halo, ax, mode='add')
  END IF

END SUBROUTINE matmul_boundary_condition


SUBROUTINE laplace_boundary_condition(comm_interf, ax, mode)

  IMPLICIT NONE
  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf
  REAL(Realkind), INTENT(inout) :: ax(:)
  CHARACTER(len=3), OPTIONAL, INTENT(in) :: mode
  

  CALL SpAVec_SpRowCol(comm_interf%halo%fdata, comm_interf%halo%Lapl3d_halo, comm_interf%halo%hdata)
    
  IF (PRESENT(mode)) THEN
    CALL copy_data_from_halo_cell(comm_interf%halo, ax, mode='sub')
  ELSE
    CALL copy_data_from_halo_cell(comm_interf%halo, ax, mode='add')
  END IF

END SUBROUTINE laplace_boundary_condition


SUBROUTINE matmul_rd_boundary_condition(comm_interf, ax, mode)

  IMPLICIT NONE
  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf
  REAL(Realkind), INTENT(inout) :: ax(:)
  CHARACTER(len=3), OPTIONAL, INTENT(in) :: mode

  CALL SpAVec_SpRowCol(comm_interf%halo%fdata, comm_interf%halo%mat_rd_halo, comm_interf%halo%hdata)

  IF (PRESENT(mode)) THEN
    CALL copy_data_from_halo_cell(comm_interf%halo, ax, mode='sub')
  ELSE
    CALL copy_data_from_halo_cell(comm_interf%halo, ax, mode='add')
  END IF

END SUBROUTINE matmul_rd_boundary_condition

SUBROUTINE matmul_bl_boundary_condition(comm_interf, ax, mode)

  IMPLICIT NONE
  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf
  REAL(Realkind), INTENT(inout) :: ax(:)
  CHARACTER(len=3), OPTIONAL, INTENT(in) :: mode


  CALL SpAVec_SpRowCol(comm_interf%halo%fdata, comm_interf%halo%mat_bl_halo, comm_interf%halo%hdata)

  IF (PRESENT(mode)) THEN
    CALL copy_data_from_halo_cell(comm_interf%halo, ax, mode='sub')
  ELSE
    CALL copy_data_from_halo_cell(comm_interf%halo, ax, mode='add')
  END IF

END SUBROUTINE matmul_bl_boundary_condition


SUBROUTINE gradient_boundary_condition(comm_interf, g)

  IMPLICIT NONE
  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf
  REAL(Realkind), INTENT(inout) :: g(:)

  CALL SpAVec_SpRowCol(comm_interf%halo%gdata, comm_interf%halo%Grad3d_halo, comm_interf%halo%hdata)
  
  CALL copy_data_from_halo_face(comm_interf%halo, g, mode='add')

END SUBROUTINE gradient_boundary_condition


SUBROUTINE gradientp_boundary_condition(comm_interf, g)
  
  IMPLICIT NONE
  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf
  REAL(Realkind), INTENT(inout) :: g(:)

  CALL SpAVec_SpRowCol(comm_interf%halo%gdata, comm_interf%halo%Gradp3d_halo, comm_interf%halo%hdata)
  
  CALL copy_data_from_halo_face(comm_interf%halo, g, mode='add')

END SUBROUTINE gradientp_boundary_condition


SUBROUTINE intcell2face_boundary_condition(comm_interf, f)

  IMPLICIT NONE
  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf
  REAL(Realkind), INTENT(inout) :: f(:)

  CALL SpAVec_SpRowCol(comm_interf%halo%gdata, comm_interf%halo%Intc2f_halo, comm_interf%halo%hdata)
  
  CALL copy_data_from_halo_face(comm_interf%halo, f, mode='add')

END SUBROUTINE intcell2face_boundary_condition


SUBROUTINE init_two_sided_exchange(comm_interf, x)
  IMPLICIT NONE

  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf
  REAL(Realkind), INTENT(in) :: x(:)

  TYPE(Bound), POINTER :: send_current, recv_current

  send_current => comm_interf%send_bnds
  recv_current => comm_interf%recv_bnds

  !First step: Copy field data on bound buffer and exchange the data

  IF (ASSOCIATED(send_current)) THEN

    CALL copy_data_on_bnd(send_current, x, mode='rep')
    
    CALL MPI_ISEND(send_current%fld_data, SIZE(send_current%fld_data), MPI_RealKind, &
                   send_current%rank_comm, send_current%tag, MPI_COMM_WORLD, send_current%req, MPIErr)

    CALL MPI_IRECV(recv_current%fld_data, SIZE(recv_current%fld_data), MPI_RealKind, &
                   recv_current%rank_comm, recv_current%tag, MPI_COMM_WORLD, recv_current%req, MPIErr)
    
    DO WHILE (ASSOCIATED(send_current%next))

      send_current => send_current%next
      recv_current => recv_current%next
      CALL copy_data_on_bnd(send_current, x, mode='rep')
      
      CALL MPI_ISEND(send_current%fld_data, SIZE(send_current%fld_data), MPI_RealKind, &
                     send_current%rank_comm, send_current%tag, MPI_COMM_WORLD, send_current%req, MPIErr)

      CALL MPI_IRECV(recv_current%fld_data, SIZE(recv_current%fld_data), MPI_RealKind, &
                     recv_current%rank_comm, recv_current%tag, MPI_COMM_WORLD, recv_current%req, MPIErr)

    END DO

  END IF

END SUBROUTINE init_two_sided_exchange



SUBROUTINE finish_two_sided_exchange(comm_interf)

  IMPLICIT NONE

  TYPE(BoundCommInterface), INTENT(inout) :: comm_interf

  TYPE(Bound), POINTER :: send_current, recv_current
 
  send_current => comm_interf%send_bnds
  recv_current => comm_interf%recv_bnds

  IF (ASSOCIATED(send_current)) THEN  

    CALL MPI_WAIT(send_current%req, status, MPIErr)
    CALL MPI_WAIT(recv_current%req, status, MPIErr)

    CALL copy_data_from_bnd(recv_current, comm_interf%halo%hdata, mode='rep')
  
    DO WHILE (ASSOCIATED(send_current%next))
      
      send_current => send_current%next
      recv_current => recv_current%next
      CALL MPI_WAIT(send_current%req, status, MPIErr)
      CALL MPI_WAIT(recv_current%req, status, MPIErr)

      CALL copy_data_from_bnd(recv_current, comm_interf%halo%hdata, mode='rep')

    END DO
  END IF 


END SUBROUTINE finish_two_sided_exchange


SUBROUTINE copy_data_on_bnd(bnd, field1d, mode)
  IMPLICIT  NONE

  TYPE(Bound), POINTER, INTENT(inout) :: bnd
  REAL(Realkind), INTENT(in) :: field1d(:)
  CHARACTER(len=3), INTENT(in) :: mode
  INTEGER :: i

  SELECT CASE(mode)

    CASE('rep')   !replace data on bound

      DO i = 1, SIZE(bnd%fld_inds)        
        bnd%fld_data(i) = field1d(bnd%fld_inds(i))
      END DO

    CASE('add')    !add data to already existing data on bound

      DO i = 1, SIZE(bnd%fld_inds)
        bnd%fld_data(i) = bnd%fld_data(i) + field1d(bnd%fld_inds(i))
      END DO

    CASE('mul')    !multiply data with already existing data on bound

      DO i = 1, SIZE(bnd%fld_inds)
        bnd%fld_data(i) = bnd%fld_data(i) * field1d(bnd%fld_inds(i))
      END DO

    END SELECT

END SUBROUTINE copy_data_on_bnd

SUBROUTINE copy_data_from_bnd(bnd, field1d, mode)
  IMPLICIT  NONE

  TYPE(Bound), POINTER, INTENT(in) :: bnd
  REAL(Realkind), INTENT(inout) :: field1d(:)
  CHARACTER(len=3), INTENT(in) :: mode
  INTEGER :: i, j

  SELECT CASE(mode)

    CASE('rep')   !replace data 

      DO i = 1, SIZE(bnd%fld_inds)
        field1d(bnd%fld_inds(i)) = bnd%fld_data(i)
      END DO

    CASE('add')    !add data to already existing data 

      DO i = 1, SIZE(bnd%fld_inds)
        field1d(bnd%fld_inds(i)) = field1d(bnd%fld_inds(i)) + bnd%fld_data(i)
      END DO

    CASE('mul')    !multiply data with already existing data

      DO i = 1, SIZE(bnd%fld_inds)
        field1d(bnd%fld_inds(i)) = field1d(bnd%fld_inds(i)) * bnd%fld_data(i)
      END DO

  END SELECT


END SUBROUTINE copy_data_from_bnd


SUBROUTINE copy_data_from_halo_cell(halo, field1d, mode)
  IMPLICIT NONE

  TYPE(HaloType), INTENT(in) :: halo
  REAL(Realkind), INTENT(inout) :: field1d(:)
  CHARACTER(len=3), INTENT(in) :: mode

  INTEGER :: i

  SELECT CASE(mode)

    CASE('rep')

    DO i = 1, SIZE(halo%finds)
      field1d(halo%finds(i)) = halo%fdata(i)
    END DO

    CASE('add')

    DO i = 1, SIZE(halo%finds)
      field1d(halo%finds(i)) = field1d(halo%finds(i)) + halo%fdata(i)
    END DO

    CASE('sub')

    DO i = 1, SIZE(halo%finds)
      field1d(halo%finds(i)) = field1d(halo%finds(i)) - halo%fdata(i)
    END DO

    CASE('mul')

    DO i = 1, SIZE(halo%finds)
      field1d(halo%finds(i)) = field1d(halo%finds(i)) * halo%fdata(i)
    END DO

  END SELECT

END SUBROUTINE copy_data_from_halo_cell


SUBROUTINE copy_data_from_halo_face(halo, field1d, mode)
  IMPLICIT NONE

  TYPE(HaloType), INTENT(in) :: halo
  REAL(Realkind), INTENT(inout) :: field1d(:)
  CHARACTER(len=3), INTENT(in) :: mode

  INTEGER :: i

  SELECT CASE(mode)

    CASE('rep')

    DO i = 1, SIZE(halo%ginds)
      field1d(halo%ginds(i)) = halo%gdata(i)
    END DO

    CASE('add')

    DO i = 1, SIZE(halo%ginds)
      field1d(halo%ginds(i)) = field1d(halo%ginds(i)) + halo%gdata(i)
    END DO

    CASE('sub')

    DO i = 1, SIZE(halo%ginds)
      field1d(halo%ginds(i)) = field1d(halo%ginds(i)) - halo%gdata(i)
    END DO

    CASE('mul')

    DO i = 1, SIZE(halo%ginds)
      field1d(halo%ginds(i)) = field1d(halo%ginds(i)) * halo%gdata(i)
    END DO

  END SELECT

END SUBROUTINE copy_data_from_halo_face



END MODULE Boundary_Mod
