MODULE SparseSolver_Mod

  USE Floor_Mod
  USE Parallel_Mod
  IMPLICIT NONE
  INCLUDE 'dmumps_struc.h'
  TYPE (DMUMPS_STRUC), POINTER :: ANeuMumps(:)
  TYPE (DMUMPS_STRUC), POINTER :: ADirMumps(:)

CONTAINS
SUBROUTINE InitMumps
! Init phase  
  WRITE(*,*) 'Init MUMPS'
  ALLOCATE(ANeuMumps(nbLoc))
  ALLOCATE(ADirMumps(nbLoc))

  DO ibLoc=1,nbLoc
    ANeuMumps(ibLoc)%COMM=MPI_COMM_SELF%mpi_val
    ANeuMumps(ibLoc)%SYM=0
    ANeuMumps(ibLoc)%PAR=1
    ANeuMumps(ibLoc)%JOB=-1
    ANeuMumps(ibLoc)%KEEP=0
    CALL DMUMPS(ANeuMumps(ibLoc))
    ADirMumps(ibLoc)%COMM=MPI_COMM_SELF%mpi_val
    ADirMumps(ibLoc)%SYM=0
    ADirMumps(ibLoc)%PAR=1
    ADirMumps(ibLoc)%JOB=-1
    ADirMumps(ibLoc)%KEEP=0
    CALL DMUMPS(ADirMumps(ibLoc))
  END DO

END SUBROUTINE InitMumps
END MODULE SparseSolver_Mod
