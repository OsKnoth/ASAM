MODULE Compressible_Mod

USE Kind_Mod
USE MatSpRowCol_Mod, ONLY: SpRowCol, SpMm_SpRowCol, SpMmBAC, SpA_SpRowCol, SpAScalRow_SpRowCol, SpAScalinvRow_SpRowCol
USE SubDom_Mod, ONLY: SubDomain
USE Constants_Mod, ONLY: Rd, kappa, p0
USE Parallel_mod

IMPLICIT NONE

REAL(Realkind), ALLOCATABLE :: theta_face(:)

CONTAINS


SUBROUTINE solveM(subdom, x_rhotheta, x_rho, x_rhovel, r_rhotheta, r_rho, r_rhovel, dt, gamm, NumIter, resmax_tol)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(inout), TARGET :: subdom
  REAL(Realkind), INTENT(inout) :: x_rhotheta(:), x_rho(:), x_rhovel(:)
  REAL(Realkind), INTENT(in) :: r_rhotheta(:), r_rho(:), r_rhovel(:)
  REAL(Realkind), INTENT(in) :: dt, gamm
  REAL(Realkind), OPTIONAL, INTENT(in) :: resmax_tol
  INTEGER, OPTIONAL, INTENT(in) :: NumIter
  REAL(Realkind) :: rhs(subdom%ncells)
  REAL(Realkind) :: face_val(subdom%nfaces)
  REAL(Realkind) :: resmax_tol_def
  INTEGER :: NumIter_def
  INTEGER :: i

  NumIter_def = 10       ! Upper bound of number of iterations
  resmax_tol_def = 1e-4  ! tolerance of maximum norm of residual

  IF (PRESENT(NumIter)) NumIter_def = NumIter
  IF (PRESENT(resmax_tol)) resmax_tol_def = resmax_tol

  DO i = 1, subdom%nfaces
    face_val(i) = -r_rhovel(i) * theta_face(i) * (dt * gamm)
  END DO

  CALL subdom%div_mul(rhs, face_val)

  DO i = 1, subdom%ncells
    rhs(i) = (rhs(i) + r_rhotheta(i))
  END DO

  x_rhotheta(:) = r_rhotheta 
  CALL subdom%solve_mg(x_rhotheta, rhs, NumIter_def, resmax_tol_def)
  
  CALL subdom%gradp_mul(x_rhovel, x_rhotheta)
  
  DO i = 1, subdom%nfaces
    x_rhovel(i) = r_rhovel(i) - dt * gamm * x_rhovel(i)
  END DO

  CALL subdom%div_mul(x_rho, x_rhovel)

  DO i = 1, subdom%ncells
     x_rho(i) = r_rho(i) - dt * gamm * x_rho(i)
  END DO

END SUBROUTINE solveM


SUBROUTINE solve_incompr(subdom, x_p, x_rhovel, r_rhovel, dt, gamm, NumIter, resmax_tol)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(inout), TARGET :: subdom
  REAL(Realkind), INTENT(inout) :: x_p(:)
  REAL(Realkind), INTENT(inout) :: x_rhovel(:)
  REAL(Realkind), INTENT(in) :: r_rhovel(:)
  REAL(Realkind), INTENT(in) :: dt, gamm
  REAL(Realkind), OPTIONAL, INTENT(in) :: resmax_tol
  INTEGER, OPTIONAL, INTENT(in) :: NumIter
  REAL(Realkind) :: rhs(subdom%ncells)
  REAL(Realkind) :: resmax_tol_def
  INTEGER :: NumIter_def
  INTEGER :: i

  NumIter_def = 10       ! Upper bound of number of iterations
  resmax_tol_def = 1e-4  ! tolerance of maximum norm of residual

  IF (PRESENT(NumIter)) NumIter_def = NumIter
  IF (PRESENT(resmax_tol)) resmax_tol_def = resmax_tol


  CALL subdom%div_mul(rhs, r_rhovel)

  DO i = 1, subdom%ncells
     rhs(i) = rhs(i) / (dt * gamm)
  END DO

  CALL subdom%solve_mg(x_p, rhs, NumIter_def, resmax_tol_def)
 
  CALL subdom%gradp_mul(x_rhovel, x_p)
 
  DO i = 1, subdom%nfaces
    x_rhovel(i) = r_rhovel(i) - dt * gamm * x_rhovel(i)
  END DO

END SUBROUTINE solve_incompr


SUBROUTINE solve_momentum(subdom, dveldt, x, b, NumIter, resmax_tol)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(inout), TARGET :: subdom
  REAL(Realkind), INTENT(inout) :: dveldt(:), x(:)
  REAL(Realkind), INTENT(in) :: b(:)
  REAL(Realkind), OPTIONAL, INTENT(in) :: resmax_tol

  TYPE(SubDomain), POINTER :: subdom_lev => NULL()
  INTEGER, OPTIONAL, INTENT(in) :: NumIter
  REAL(Realkind) :: resmax_tol_def

  INTEGER :: NumIter_def

  NumIter_def = 10       ! Upper bound of number of iterations
  resmax_tol_def = 1e-4  ! tolerance of maximum norm of residual

  IF (PRESENT(NumIter)) NumIter_def = NumIter
  IF (PRESENT(resmax_tol)) resmax_tol_def = resmax_tol

  CALL subdom%solve_mg(x, b, NumIter_def, resmax_tol_def)
!  CALL subdom%solve_mgbicgstab(rhotheta_new, rhs, NumIter_def, resmax_tol_def) 
  CALL subdom%gradp_mul(dveldt, -x)

END SUBROUTINE solve_momentum


RECURSIVE SUBROUTINE setup_pressure_solver_incompr(subdom)
  IMPLICIT NONE
  TYPE(SubDomain), INTENT(inout) :: subdom

  LOGICAL :: create

  create = .FALSE.

  IF (.NOT. ASSOCIATED(subdom%gradp3d%RowPtr)) THEN
    subdom%gradp3d%m = subdom%grad3d%m
    subdom%gradp3d%n = subdom%grad3d%n
    ALLOCATE(subdom%gradp3d%RowPtr(subdom%grad3d%m + 1))
    ALLOCATE(subdom%gradp3d%ColInd(SIZE(subdom%grad3d%ColInd, 1)))
    ALLOCATE(subdom%gradp3d%Val(SIZE(subdom%grad3d%Val, 1)))
    create = .TRUE.
  END IF

  subdom%gradp3d%RowPtr(:) = subdom%grad3d%RowPtr
  subdom%gradp3d%ColInd(:) = subdom%grad3d%ColInd
  subdom%gradp3d%Val(:) = subdom%grad3d%Val

  IF (create) ALLOCATE(subdom%mat)
  CALL SpMm_SpRowCol(subdom%mat, subdom%div3d, subdom%gradp3d, create, create)

  IF (ASSOCIATED(subdom%bc_if)) THEN
    CALL setup_pressure_incrompr_solver_bc_if(subdom, create)
  END IF

  IF (ASSOCIATED(subdom%next_coarse)) THEN
    CALL setup_pressure_solver_incompr(subdom%next_coarse)
  END IF

  CALL subdom%init_smoother()

END SUBROUTINE setup_pressure_solver_incompr


RECURSIVE SUBROUTINE setup_pressure_solver(subdom, dt, gamm, field1, rho)
  IMPLICIT NONE
  
  TYPE(SubDomain), INTENT(inout) :: subdom
  REAL(Realkind), INTENT(in) :: dt, gamm
  REAL(Realkind), INTENT(in) :: field1(:)
  REAL(Realkind), OPTIONAL, INTENT(in) :: rho(:)

  REAL(Realkind) :: field_cell_tmp(subdom%ncells)
  REAL(Realkind) :: field_face_tmp(subdom%nfaces)
  REAL(Realkind) :: c_sq_coarse(subdom%ncoarse)

  REAL(Realkind), TARGET :: laplp_val(SIZE(subdom%lapl3d%ColInd, 1))
  INTEGER, TARGET :: laplp_colind(SIZE(subdom%lapl3d%ColInd, 1))
  INTEGER, TARGET :: laplp_rowptr(subdom%lapl3d%m + 1)

  TYPE(SpRowCol) :: laplp
  REAL(Realkind) :: cons1, exp1, dt_sq
  REAL(Realkind), PARAMETER :: eps = 1e-20
  REAL(Realkind), PARAMETER :: one_real = 1.0
  INTEGER :: i
  LOGICAL :: create

  REAL(Realkind) :: lapl_diag(subdom%lapl3d%m)


  create = .FALSE.

  IF (.NOT. ASSOCIATED(subdom%gradp3d%RowPtr)) THEN
    subdom%gradp3d%m = subdom%grad3d%m
    subdom%gradp3d%n = subdom%grad3d%n
    ALLOCATE(subdom%gradp3d%RowPtr(subdom%grad3d%m + 1))
    ALLOCATE(subdom%gradp3d%ColInd(SIZE(subdom%grad3d%ColInd, 1)))
    ALLOCATE(subdom%gradp3d%Val(SIZE(subdom%grad3d%Val, 1)))
    create = .TRUE.
  END IF

  IF (.NOT. ALLOCATED(theta_face)) ALLOCATE(theta_face(subdom%nfaces))

  subdom%gradp3d%RowPtr(:) = subdom%grad3d%RowPtr
  subdom%gradp3d%ColInd(:) = subdom%grad3d%ColInd
  subdom%gradp3d%Val(:) = subdom%grad3d%Val

  laplp%RowPtr => laplp_rowptr
  laplp%ColInd => laplp_colind
  laplp%Val => laplp_val

  dt_sq = (dt * gamm) ** 2


  IF (PRESENT(rho)) THEN !field1 is rhotheta

    exp1 = one_real / (one_real - kappa)
    cons1 = Rd ** exp1 * p0 ** (one_real - exp1) / (one_real - kappa)

    DO i = 1, subdom%ncells  ! compute c_sq
      field_cell_tmp(i) = cons1 * field1(i) ** exp1 / (rho(i) + eps)
    END DO

    CALL subdom%intc2f_mul(field_face_tmp, field_cell_tmp)
    CALL SpAScalRow_SpRowCol(subdom%gradp3d, field_face_tmp)
    CALL SpMm_SpRowCol(laplp, subdom%div3d, subdom%gradp3d, create, withdiagptr=.True.)

    DO i = 1, laplp%RowPtr(laplp%m + 1) - 1
      laplp%Val(i) = -laplp%Val(i) * dt_sq
    END DO

    DO i = 1, subdom%eye%RowPtr(subdom%eye%m + 1) - 1
      subdom%eye%Val(i) = subdom%eye%Val(i)
    END DO

    IF (create) ALLOCATE(subdom%mat)

    DO i = 1, SIZE(lapl_diag, 1)
      lapl_diag(i) = laplp%Val(laplp%DiagPtr(i))
    END DO


    CALL SpA_SpRowCol(subdom%mat, subdom%eye, laplp, create, create)

    IF (ASSOCIATED(subdom%bc_if)) THEN
      CALL setup_pressure_solver_bc_if(subdom, field_face_tmp, dt_sq, create)
    END IF

    IF (ASSOCIATED(subdom%next_coarse)) THEN
      CALL subdom%restrict(field_cell_tmp, c_sq_coarse)
      CALL setup_pressure_solver(subdom%next_coarse, dt, gamm, c_sq_coarse)
    END IF 

    DO i = 1, subdom%ncells ! inverse scaling of gradient operator with theta
      field_cell_tmp(i) = field1(i) / (rho(i) + eps)
    END DO

    CALL subdom%intc2f_mul(theta_face, field_cell_tmp)
    CALL SpAScalinvRow_SpRowCol(subdom%gradp3d, theta_face)

    IF (ASSOCIATED(subdom%bc_if)) THEN   
      CALL scaleinv_gradient_bc_if(subdom, theta_face)
    END IF

  ELSE !field1 is already c_sq

    CALL subdom%intc2f_mul(field_face_tmp, field1)   
    CALL SpAScalRow_SpRowCol(subdom%gradp3d, field_face_tmp)
    CALL SpMm_SpRowCol(laplp, subdom%div3d, subdom%gradp3d, create, withdiagptr=.True.)

    DO i = 1, laplp%RowPtr(laplp%m + 1) - 1
      laplp%Val(i) = -laplp%Val(i) * dt_sq
    END DO

    DO i = 1, subdom%eye%RowPtr(subdom%eye%m + 1) - 1
      subdom%eye%Val(i) = subdom%eye%Val(i)
    END DO

    DO i = 1, SIZE(lapl_diag, 1)
      lapl_diag(i) = laplp%Val(laplp%DiagPtr(i))
    END DO

    IF (create) ALLOCATE(subdom%mat)
    CALL SpA_SpRowCol(subdom%mat, subdom%eye, laplp, create, create)
 
    IF (ASSOCIATED(subdom%bc_if)) THEN
      CALL setup_pressure_solver_bc_if(subdom, field_face_tmp, dt_sq, create)
    END IF

    IF (ASSOCIATED(subdom%next_coarse)) THEN
      CALL subdom%restrict(field1, c_sq_coarse)
      CALL setup_pressure_solver(subdom%next_coarse, dt, gamm, c_sq_coarse)
    END IF
  END IF

  CALL subdom%init_smoother()

END SUBROUTINE setup_pressure_solver


RECURSIVE SUBROUTINE update_pressure_solver_dt(subdom, dt_old, dt_new)
  TYPE(SubDomain), INTENT(inout) :: subdom
  REAL(Realkind), INTENT(in) :: dt_old, dt_new

  REAL(Realkind) :: scal
  REAL(Realkind), PARAMETER :: one = 1.0
  INTEGER :: i

  scal = (dt_new / dt_old) ** 2

  DO i = 1, subdom%mat%RowPtr(subdom%mat%m + 1) - 1
    subdom%mat%Val(i) = subdom%mat%Val(i) * scal
  END DO

  DO i = 1, subdom%mat%n
    subdom%mat%Val(subdom%mat%DiagPtr(i)) = subdom%mat%Val(subdom%mat%DiagPtr(i)) - (scal - one) * subdom%eye%Val(i)
  END DO

  IF (ASSOCIATED(subdom%bc_if)) THEN
    CALL update_dt_solver_bc_if(subdom, scal)
  END IF

  IF (ASSOCIATED(subdom%next_coarse)) THEN
    CALL update_pressure_solver_dt(subdom%next_coarse, dt_old, dt_new)
  END IF  

  CALL subdom%init_smoother()

END SUBROUTINE update_pressure_solver_dt

SUBROUTINE setup_pressure_incrompr_solver_bc_if(subdom, create)
  IMPLICIT NONE
  TYPE(SubDomain), INTENT(inout) :: subdom
  LOGICAL, OPTIONAL, INTENT(in) :: create
  LOGICAL :: create_def

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

  IF (create_def) ALLOCATE(subdom%bc_if%halo%mat_halo)
  CALL SpMm_SpRowCol(subdom%bc_if%halo%mat_halo, subdom%bc_if%halo%Div3d_halo, subdom%bc_if%halo%Gradp3d_halo, &
                       restruct_mat=create_def, withdiagptr=.False.)

END SUBROUTINE setup_pressure_incrompr_solver_bc_if

SUBROUTINE setup_pressure_solver_bc_if(subdom, c_sq_face, dt_sq, create)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(inout) :: subdom
  REAL(Realkind), INTENT(in) :: c_sq_face(:)
  REAL(Realkind), INTENT(in) :: dt_sq
  LOGICAL, OPTIONAL, INTENT(in) :: create

  REAL(Realkind) :: c_sq_face_bound(SIZE(subdom%bc_if%halo%ginds, 1))
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
    c_sq_face_bound(i) = c_sq_face(subdom%bc_if%halo%ginds(i))
  END DO

  CALL SpAScalRow_SpRowCol(subdom%bc_if%halo%Gradp3d_halo, c_sq_face_bound)
  
  IF (create_def) ALLOCATE(subdom%bc_if%halo%mat_halo)
  CALL SpMm_SpRowCol(subdom%bc_if%halo%mat_halo, subdom%bc_if%halo%Div3d_halo, subdom%bc_if%halo%Gradp3d_halo, &
                       restruct_mat=create_def, withdiagptr=.False.)
  DO i = 1, subdom%bc_if%halo%mat_halo%RowPtr(subdom%bc_if%halo%mat_halo%m + 1) - 1
    subdom%bc_if%halo%mat_halo%Val(i) = -subdom%bc_if%halo%mat_halo%Val(i) * dt_sq
  END DO
!  CALL SpMm_SpRowCol(subdom%bc_if%halo%mat_halo, subdom%bc_if%halo%Div3d_halo, subdom%bc_if%halo%Grad3d_halo, &
!                      restruct_mat=create_def, withdiagptr=.False.)

END SUBROUTINE


SUBROUTINE scaleinv_gradient_bc_if(subdom, scal_face)
  TYPE(SubDomain), INTENT(inout) :: subdom
  REAL(Realkind), INTENT(in) :: scal_face(:)
  REAL(Realkind) :: scal_face_bound(SIZE(subdom%bc_if%halo%ginds, 1))
  INTEGER :: i

  DO i = 1, SIZE(subdom%bc_if%halo%ginds, 1)
    scal_face_bound(i) = scal_face(subdom%bc_if%halo%ginds(i))
  END DO

  CALL SpAScalinvRow_SpRowCol(subdom%bc_if%halo%Gradp3d_halo, scal_face_bound)

END SUBROUTINE scaleinv_gradient_bc_if


SUBROUTINE update_dt_solver_bc_if(subdom, scal)
  IMPLICIT NONE

  TYPE(SubDomain), INTENT(inout) :: subdom
  REAL(Realkind), INTENT(in) :: scal
  INTEGER :: i

  DO i = 1, subdom%bc_if%halo%mat_halo%RowPtr(subdom%bc_if%halo%mat_halo%m + 1) - 1
    subdom%bc_if%halo%mat_halo%Val(i) = subdom%bc_if%halo%mat_halo%Val(i) * scal
  END DO

END SUBROUTINE update_dt_solver_bc_if

END MODULE Compressible_Mod
