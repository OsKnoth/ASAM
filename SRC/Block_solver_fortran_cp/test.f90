PROGRAM test

USE Subdom_Mod, ONLY: SubDomain
USE DomDecomp_Mod, ONLY: read_blocks_from_netcdf, vec1d_from_homog, btype, blockranks
USE Laplace_block_Mod, ONLY: make_operators
USE Compressible_Mod, ONLY: setup_pressure_solver
USE Kind_Mod
USE MPI_mod
USE Output_Mod, ONLY: output_serial
USE List_Mod, ONLY: RArray

USE CoarseGrids_Mod, ONLY: define_coarse_grids, partition_clevels, intialize_coarse_subdomains, &
                           subdomain_lev, isblackijkzero_lev, determine_color_startcell

USE MatSpRowCol_Mod, ONLY: SpAVec_SpRowCol

USE Part_Man_Mod, ONLY: set_ranks

IMPLICIT NONE

!CHARACTER(len=128), PARAMETER :: file_path = "Jena_5m_ffac_100blocks.nc"
!CHARACTER(len=128), PARAMETER :: file_path = "Haus.nc"
!CHARACTER(len=128), PARAMETER :: file_path = "CanyonHof.nc"
CHARACTER(len=128), PARAMETER :: file_path = "Eiba.nc"
CHARACTER(len=3) :: file_nr 

REAL(Realkind), PARAMETER :: theta = 280.0  !potential temperature
REAL(Realkind), PARAMETER :: c = 320.0      !speed of sound
REAL(Realkind), PARAMETER :: rho = 1.0      !air density
REAL(Realkind), PARAMETER :: dt = 0.1       !time step

REAL(Realkind), PARAMETER :: sor_param = 1.60
INTEGER, PARAMETER :: NumIter = 10

TYPE(SubDomain), TARGET :: subdom
REAL(Realkind), ALLOCATABLE, TARGET, DIMENSION(:) :: uvec1d, x, b, r, g
REAL(Realkind), ALLOCATABLE :: c_sq_over_rhotheta(:), &
                               rhotheta(:)

TYPE(RArray), ALLOCATABLE :: outfields(:)
INTEGER, ALLOCATABLE :: outtypes(:)
CHARACTER(LEN=20), ALLOCATABLE :: outnames(:)

TYPE(SubDomain), POINTER :: this_subdomain => NULL()

REAL(Realkind) :: u, v, w
INTEGER :: n, i, j, k, ngrids, ilev, icell

REAL(Realkind) :: wctime_st, wctime_end

u = 1.0
v = 2.0
w = 0.0

CALL start_MPI(subdom%myrank)
CALL set_ranks(blockranks)

CALL read_blocks_from_netcdf(file_path, subdom)

CALL define_coarse_grids()
CALL partition_clevels()
CALL determine_color_startcell()

CALL intialize_coarse_subdomains(subdom)
WRITE(*,*) "blocks read"

this_subdomain => subdomain_lev

CALL this_subdomain%set_cell_type()
CALL this_subdomain%set_face_type()

!make the communication interface
CALL this_subdomain%init_boundcomm_interface(blockranks)

CALL make_operators(this_subdomain, btype, isblackijkzero_lev(1)%data)

ALLOCATE(c_sq_over_rhotheta(this_subdomain%ncells))
ALLOCATE(rhotheta(this_subdomain%ncells))
c_sq_over_rhotheta(:) = c ** 2.0 / (rho * theta)
rhotheta(:) = rho * theta

CALL setup_pressure_solver(this_subdomain, c_sq_over_rhotheta, rhotheta, dt)

CALL this_subdomain%init_smoother(sor_param)

ngrids = 1
DO WHILE(ASSOCIATED(this_subdomain%next_coarse))
  this_subdomain => this_subdomain%next_coarse
  CALL this_subdomain%init_smoother(sor_param)
  ngrids = ngrids + 1
END DO

this_subdomain => subdomain_lev

ALLOCATE(x(this_subdomain%ncells))
ALLOCATE(b(this_subdomain%ncells))
ALLOCATE(r(this_subdomain%ncells))
ALLOCATE(g(this_subdomain%nfaces))

!Create solution, rhs and residual vector

ALLOCATE(uvec1d(this_subdomain%nfaces))

CALL vec1d_from_homog(this_subdomain, u, v, w, uvec1d)
CALL this_subdomain%div_mul(b, uvec1d) 

x(:) = 0.0

WRITE(*,*) "solve"

CALL CPU_TIME(wctime_st)

CALL this_subdomain%solve_mgbicgstab(x, b, NumIter)
!CALL this_subdomain%solve_mg(x, b, NumIter)

CALL CPU_TIME(wctime_end)
WRITE(*,*) wctime_end - wctime_st

WRITE(file_nr, '(I3.3)') this_subdomain%myrank

CALL this_subdomain%grad_mul(g, x)

ALLOCATE(outfields(2))
ALLOCATE(outtypes(2))
ALLOCATE(outnames(2))

outtypes(1:2) = (/1, 3/)
outnames(1:2) = (/'solution', 'gradient'/)

outfields(1)%data => x
outfields(2)%data => g

IF (ASSOCIATED(this_subdomain%blocks)) THEN
  CALL output_serial(this_subdomain, outfields, outtypes, outnames, "output"//trim(file_nr)//".nc")
END IF

CALL finalize_MPI()

END PROGRAM test
