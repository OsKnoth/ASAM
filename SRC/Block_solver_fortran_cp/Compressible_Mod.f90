MODULE Compressible_Mod

USE Kind_Mod
USE MatSpRowCol_Mod, ONLY: SpRowCol, SpMm_SpRowCol, SpMmBAC, SpA_SpRowCol, SpAScalRow_SpRowCol
USE SubDom_Mod, ONLY: SubDomain

IMPLICIT NONE

CONTAINS

RECURSIVE SUBROUTINE setup_pressure_solver(subdom, c_sq_over_rhotheta, rhotheta, dt, create)
  IMPLICIT NONE
  
  TYPE(SubDomain), INTENT(inout) :: subdom
  REAL(Realkind), INTENT(in) :: c_sq_over_rhotheta(:), rhotheta(:)
  REAL(Realkind) :: dt
  LOGICAL, OPTIONAL, INTENT(in) :: create

  REAL(Realkind) :: face_val(subdom%nfaces)

  REAL(Realkind) :: c_sq_over_rhotheta_coarse(subdom%ncoarse), &
                    rhotheta_coarse(subdom%ncoarse)

  REAL(Realkind), TARGET :: laplp_val(SIZE(subdom%lapl3d%ColInd, 1))
  INTEGER, TARGET :: laplp_colind(SIZE(subdom%lapl3d%ColInd, 1))
  INTEGER, TARGET :: laplp_rowptr(subdom%lapl3d%m + 1)

  TYPE(SpRowCol) :: laplp

  LOGICAL :: create_def

  create_def = .True.

  IF (PRESENT(create)) create_def = create


  IF (create_def) THEN
    subdom%gradp3d%m = subdom%grad3d%m
    subdom%gradp3d%n = subdom%grad3d%n
    ALLOCATE(subdom%gradp3d%RowPtr(subdom%grad3d%m + 1))
    ALLOCATE(subdom%gradp3d%ColInd(SIZE(subdom%grad3d%ColInd, 1)))
    ALLOCATE(subdom%gradp3d%Val(SIZE(subdom%grad3d%Val, 1)))
  END IF

  subdom%gradp3d%RowPtr(:) = subdom%grad3d%RowPtr
  subdom%gradp3d%ColInd(:) = subdom%grad3d%ColInd
  subdom%gradp3d%Val(:) = subdom%grad3d%Val

  laplp%RowPtr => laplp_rowptr
  laplp%ColInd => laplp_colind
  laplp%Val => laplp_val

  CALL subdom%intc2f_mul(face_val, c_sq_over_rhotheta)

  CALL SpAScalRow_SpRowCol(subdom%gradp3d, face_val)
  CALL SpMm_SpRowCol(laplp, subdom%div3d, subdom%gradp3d, allocate_mat=.False., withdiagptr=.False.)
  CALL SpAScalRow_SpRowCol(laplp, -dt ** 2 * rhotheta)

  IF (create_def) ALLOCATE(subdom%mat)
  CALL SpA_SpRowCol(subdom%mat, subdom%eye, laplp, create_def, create_def)
!  CALL SpMm_SpRowCol(subdom%mat, subdom%div3d, subdom%grad3d, create_def, create_def)
 
  IF (ASSOCIATED(subdom%bc_if)) THEN
    CALL setup_pressure_solver_bc_if(subdom, face_val, rhotheta, dt, create_def)
  END IF

  IF (ASSOCIATED(subdom%next_coarse)) THEN
    CALL subdom%restrict(c_sq_over_rhotheta, c_sq_over_rhotheta_coarse)
    CALL subdom%restrict(rhotheta, rhotheta_coarse)
    CALL setup_pressure_solver(subdom%next_coarse, c_sq_over_rhotheta_coarse, rhotheta_coarse, dt, create_def)
  END IF

END SUBROUTINE setup_pressure_solver


SUBROUTINE setup_pressure_solver_bc_if(subdom, c_sq_over_rhotheta_face, rhotheta, dt, create)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(inout) :: subdom
  REAL(Realkind), INTENT(in) :: c_sq_over_rhotheta_face(:), rhotheta(:)
  REAL(Realkind), INTENT(in) :: dt
  LOGICAL, OPTIONAL, INTENT(in) :: create

  REAL(Realkind) :: c_sq_over_rhotheta_face_bound(SIZE(subdom%bc_if%halo%ginds, 1))
  REAL(Realkind) :: rhotheta_bound(SIZE(subdom%bc_if%halo%finds, 1)) 

  LOGICAL :: create_def

  INTEGER :: i

  IF (PRESENT(create)) create_def = create
  
  IF (create) THEN
    subdom%bc_if%halo%Gradp3d_halo%m = subdom%bc_if%halo%Grad3d_halo%m
    subdom%bc_if%halo%Gradp3d_halo%n = subdom%bc_if%halo%Grad3d_halo%n
    ALLOCATE(subdom%bc_if%halo%Gradp3d_halo%RowPtr(subdom%bc_if%halo%Grad3d_halo%m + 1))
    ALLOCATE(subdom%bc_if%halo%Gradp3d_halo%ColInd(SIZE(subdom%bc_if%halo%Grad3d_halo%ColInd, 1)))
    ALLOCATE(subdom%bc_if%halo%Gradp3d_halo%Val(SIZE(subdom%bc_if%halo%Grad3d_halo%Val, 1)))
  END IF

  subdom%bc_if%halo%Gradp3d_halo%RowPtr(:) = subdom%bc_if%halo%Grad3d_halo%RowPtr
  subdom%bc_if%halo%Gradp3d_halo%ColInd(:) = subdom%bc_if%halo%Grad3d_halo%ColInd
  subdom%bc_if%halo%Gradp3d_halo%Val(:) = subdom%bc_if%halo%Grad3d_halo%Val

  DO i = 1, SIZE(subdom%bc_if%halo%ginds, 1)
    c_sq_over_rhotheta_face_bound(i) = c_sq_over_rhotheta_face(subdom%bc_if%halo%ginds(i))
  END DO

  DO i = 1, SIZE(subdom%bc_if%halo%finds, 1)
    rhotheta_bound(i) = rhotheta(subdom%bc_if%halo%finds(i))
  END DO

  CALL SpAScalRow_SpRowCol(subdom%bc_if%halo%Gradp3d_halo, c_sq_over_rhotheta_face_bound)
  
  IF (create_def) ALLOCATE(subdom%bc_if%halo%mat_halo)
  CALL SpMm_SpRowCol(subdom%bc_if%halo%mat_halo, subdom%bc_if%halo%Div3d_halo, subdom%bc_if%halo%Gradp3d_halo, &
                       allocate_mat=create_def, restruct_mat=create_def, withdiagptr=.False.)
  CALL SpAScalRow_SpRowCol(subdom%bc_if%halo%mat_halo,  -dt ** 2 * rhotheta_bound)

END SUBROUTINE


END MODULE Compressible_Mod
