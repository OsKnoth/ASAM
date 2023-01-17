PROGRAM test

USE Subdom_Mod, ONLY: SubDomain
USE DomDecomp_Mod, ONLY: read_blocks_from_netcdf, vec1d_from_homog, btype, blockranks, x2glob, y2glob, z2glob
USE Laplace_block_Mod, ONLY: make_operators_static
USE Compressible_Mod, ONLY: solve, setup_pressure_solver, update_pressure_solver_dt, setup_pressure_solver_incompr, solve_incompr
USE Kind_Mod
USE Parallel_Mod
USE Output_Mod, ONLY: output_serial
USE List_Mod, ONLY: RArray

USE CoarseGrids_Mod, ONLY: define_coarse_grids, partition_clevels, intialize_coarse_subdomains, &
                           subdomain_lev, isblackijkzero_lev, determine_color_startcell

USE Boundary_Mod, ONLY: cyclic_x, cyclic_y, cyclic_z

USE MatSpRowCol_Mod, ONLY: SpAVec_SpRowCol

USE Constants_Mod, ONLY: R, M_dry, kappa, p0

IMPLICIT NONE

CHARACTER(len=128), PARAMETER :: file_path = "CanyonHof.nc"
CHARACTER(len=3) :: file_nr 
TYPE(SubDomain), TARGET :: subdom
REAL(Realkind), ALLOCATABLE, TARGET :: rho_new(:), rhovel_new(:), rhotheta_new(:), x_coarse(:), x_fine(:)
REAL(Realkind), ALLOCATABLE :: rhovel(:), rho(:), rhotheta(:)

TYPE(RArray), ALLOCATABLE :: outfields(:)
INTEGER, ALLOCATABLE :: outtypes(:)
CHARACTER(LEN=20), ALLOCATABLE :: outnames(:)

TYPE(SubDomain), POINTER :: this_subdomain => NULL()

REAL(Realkind), PARAMETER :: urho0 = 0.0, &
                             vrho0 = 0.0, &
                             wrho0 = 0.0, &
                             rho0 = 1.2, &
                             theta0 = 293.0, &
                             dt = 1.0, &
                             gamm = 1.0

Real(Realkind) :: resmax_tol = 1e-12
INTEGER, PARAMETER :: NumIter = 100

INTEGER :: n, i, j, k, ngrids, ilev, icell

REAL(Realkind) :: wctime_st, wctime_end

! setting lateral boundary conditions to cyclic
cyclic_x = .TRUE.
cyclic_y = .TRUE.
cyclic_z = .FALSE.

CALL start_MPI(subdom%myrank)
CALL read_blocks_from_netcdf(file_path, subdom)

CALL define_coarse_grids()
CALL partition_clevels()
CALL determine_color_startcell()

CALL intialize_coarse_subdomains(subdom)
IF (subdom%myrank .EQ. 0) WRITE(*,*) "Blocks read"

this_subdomain => subdomain_lev
CALL this_subdomain%set_cell_type()
CALL this_subdomain%set_face_type()

!make the communication interface
CALL this_subdomain%init_boundcomm_interface(blockranks)

CALL make_operators_static(this_subdomain, btype, isblackijkzero_lev(1)%data)

ALLOCATE(rhovel(this_subdomain%nfaces))
ALLOCATE(rho(this_subdomain%ncells))
ALLOCATE(rhotheta(this_subdomain%ncells))
!ALLOCATE(x(this_subdomain%ncells))

CALL vec1d_from_homog(this_subdomain, urho0, vrho0, wrho0, rhovel)

icell = 1
DO n = 1, this_subdomain%nblocks
  IF (this_subdomain%blockiscomp(n)) THEN
    DO i = 1, this_subdomain%blocks(n)%fld_shape(3)
      DO j = 1, this_subdomain%blocks(n)%fld_shape(2)
        DO k = 1, this_subdomain%blocks(n)%fld_shape(1)
          IF (this_subdomain%blocks(n)%volseff(k, j, i) .LT. 1.001e-20) THEN
            rho(icell) = 0
            rhotheta(icell) = 0
!            x(icell) = this_subdomain%blocks(n)%x(i)
          ELSE  
            rho(icell) = rho0
            rhotheta(icell) = theta0 * rho0
!            x(icell) = this_subdomain%blocks(n)%x(i)
          END IF
          icell = icell + 1
        END DO
      END DO
    END DO
  END IF
END DO

ALLOCATE(rhotheta_new(this_subdomain%ncells))
ALLOCATE(rho_new(this_subdomain%ncells))
ALLOCATE(rhovel_new(this_subdomain%nfaces))

IF (this_subdomain%myrank .EQ. 0) WRITE(*,*) "Setup pressure matrix"
!CALL setup_pressure_solver(this_subdomain, dt, gamm, rhotheta, rho)
CALL setup_pressure_solver_incompr(this_subdomain)

IF (this_subdomain%myrank .EQ. 0) WRITE(*,*) "Solve"
!CALL solve(this_subdomain, rhotheta_new, rho_new, rhovel_new, rhotheta, rho, rhovel, dt, gamm, NumIter, resmax_tol)
rhotheta_new(:) = 0.0
CALL solve_incompr(subdom, rhotheta_new, rhovel_new, rhovel, dt, gamm, NumIter, resmax_tol)

!ALLOCATE(x_coarse(this_subdomain%next_coarse%ncells))
!CALL this_subdomain%restrict(rhotheta_new, x_coarse)
!this_subdomain => this_subdomain%next_coarse
!ALLOCATE(x_fine(this_subdomain%ncells))
!x_fine(:) = x_coarse
!DEALLOCATE(x_coarse)
!ALLOCATE(x_coarse(this_subdomain%next_coarse%ncells))
!CALL this_subdomain%restrict(x_fine, x_coarse)
!this_subdomain => this_subdomain%next_coarse
!DEALLOCATE(x_fine)
!ALLOCATE(x_fine(this_subdomain%ncells))
!x_fine(:) = x_coarse
!DEALLOCATE(x_coarse)
!ALLOCATE(x_coarse(this_subdomain%next_coarse%ncells))
!CALL this_subdomain%restrict(x_fine, x_coarse)
!this_subdomain => this_subdomain%next_coarse
!DEALLOCATE(x_fine)
!ALLOCATE(x_fine(this_subdomain%ncells))
!x_fine(:) = x_coarse
!DEALLOCATE(x_coarse)
!ALLOCATE(x_coarse(this_subdomain%next_coarse%ncells))
!CALL this_subdomain%restrict(x_fine, x_coarse)
!CALL this_subdomain%prolongate(x_coarse, x_fine, mode="rep")
!DEALLOCATE(x_coarse)
!ALLOCATE(x_coarse(this_subdomain%ncells))
!x_coarse(:) = x_fine
!this_subdomain => this_subdomain%next_fine
!DEALLOCATE(x_fine)
!ALLOCATE(x_fine(this_subdomain%ncells))
!CALL this_subdomain%prolongate(x_coarse, x_fine, mode="rep")
!DEALLOCATE(x_coarse)
!ALLOCATE(x_coarse(this_subdomain%ncells))
!x_coarse(:) = x_fine
!this_subdomain => this_subdomain%next_fine
!DEALLOCATE(x_fine)
!ALLOCATE(x_fine(this_subdomain%ncells))
!CALL this_subdomain%prolongate(x_coarse, x_fine, mode="rep")
!DEALLOCATE(x_coarse)
!ALLOCATE(x_coarse(this_subdomain%ncells))
!x_coarse(:) = x_fine
!this_subdomain => this_subdomain%next_fine
!DEALLOCATE(x_fine)
!ALLOCATE(x_fine(this_subdomain%ncells))
!CALL this_subdomain%prolongate(x_coarse, x_fine, mode="rep")





!DO i = 1, this_subdomain%blocks(1)%fld_shape(3)
!   rhotheta_new(i) = i
!END DO

!CALL this_subdomain%gradp_mul(rhovel_new, rhotheta_new)


!CALL update_pressure_solver_dt(this_subdomain, dt, 2.0 * dt)
!CALL solve(this_subdomain, rhotheta_new, rho_new, rhovel_new, rhotheta, rho, rhovel, 2.0 * dt, gamm, NumIter, resmax_tol)

WRITE(file_nr, '(I3.3)') this_subdomain%myrank

ALLOCATE(outfields(1))
ALLOCATE(outtypes(1))
ALLOCATE(outnames(1))

outtypes(:) = (/1/)!, 1, 3/)
outnames(:) = (/            &
                'rhotheta' & 
!                'rho     ', &
!                'rhovel  '  &
              /)

!outfields(1)%data => x_fine!rhotheta_new
outfields(1)%data => rhotheta_new
!outfields(2)%data => rho_new
!outfields(3)%data => rhovel_new

IF (ASSOCIATED(this_subdomain%blocks)) THEN
  CALL output_serial(this_subdomain, outfields, outtypes, outnames, "output"//trim(file_nr)//".nc")
END IF

CALL finalize_MPI()

END PROGRAM test
