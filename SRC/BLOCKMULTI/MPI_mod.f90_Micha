MODULE MPI_Mod

USE Kind_Mod

include 'mpif.h'

INTEGER :: MPI_RealKind
INTEGER :: NumProcs
INTEGER :: MPIErr
LOGICAL :: flag
integer :: status(MPI_STATUS_SIZE)

CONTAINS


!-----------------------------------------------------------------------
!     MPI Initialization and finalization

SUBROUTINE start_MPI(myrank)

  INTEGER, INTENT(inout) :: myrank

      CALL MPI_INIT( MPIErr )
      IF (MPIErr .ne. MPI_SUCCESS) then
         WRITE(*,*) 'MPI Initialisierungsfehler '
         CALL MPI_FINALIZE(MPIErr)
         STOP
      END IF
!
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, myrank, MPIErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NumProcs, MPIErr )

      IF (RealKind==4) THEN
        MPI_RealKind=MPI_Real
      ELSE
        MPI_RealKind=MPI_Double_Precision
      END IF

END SUBROUTINE start_MPI

SUBROUTINE finalize_MPI()

  CALL MPI_FINALIZE(flag, MPIErr)

END SUBROUTINE finalize_MPI


END MODULE MPI_Mod
