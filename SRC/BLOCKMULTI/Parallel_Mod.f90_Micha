MODULE Parallel_Mod

   USE Kind_Mod
   USE Domain_Mod
   USE mpi_f08

   IMPLICIT NONE
   TYPE(mpi_datatype) :: MPI_RealKind

   TYPE block_procT
     INTEGER :: proc
     INTEGER :: ibLoc
   END TYPE block_procT
   TYPE CommGraphT
     INTEGER :: NumNgbr
     INTEGER, POINTER :: ProcId(:)
     INTEGER, POINTER :: BackId(:)
     INTEGER, POINTER :: Order(:)
     INTEGER, POINTER :: GetLen(:)
     INTEGER, POINTER :: PutLen(:)
   END TYPE CommGraphT
   TYPE kette
     INTEGER :: nlink
     INTEGER :: LenPack
     INTEGER :: LenUnPack
     TYPE(gliedT), POINTER :: glied(:)
   END TYPE kette

   INTEGER :: MyId
   INTEGER :: NumProcs
   INTEGER :: MPIErr
   LOGICAL :: flag
   INTEGER :: MaxBlkNgbr
   INTEGER :: nbGlob
   INTEGER :: nbLoc
   INTEGER, ALLOCATABLE :: LocGlob(:)
   TYPE (block_procT), ALLOCATABLE :: blMPI(:)
   INTEGER, POINTER :: tmpIb(:)               ! needed in ParMetis
   TYPE(CommGraphT), POINTER :: CoGr(:)
   INTEGER ::  MaxCommLevel
   TYPE (kette), POINTER :: chain(:)
   INTEGER :: NumberNeiProc
   INTEGER :: MaxPutBuf, MaxGetBuf
   

CONTAINS

SUBROUTINE start_MPI()

!-----------------------------------------------------------------------
!     MPI Initialization

      CALL MPI_INIT( MPIErr )
      IF (MPIErr .ne. MPI_SUCCESS) then
         WRITE(*,*) 'MPI Initialisierungsfehler '
         CALL MPI_FINALIZE(MPIErr)
         STOP
      END IF
!
      CALL MPI_COMM_RANK( MPI_COMM_WORLD, MyId, MPIErr )
      CALL MPI_COMM_SIZE( MPI_COMM_WORLD, NumProcs, MPIErr )

      IF (RealKind==4) THEN
        MPI_RealKind=MPI_Real
      ELSE
        MPI_RealKind=MPI_Double_Precision
      END IF
!
!---------------------------------------------
END SUBROUTINE start_MPI


SUBROUTINE finalize_MPI()

  CALL MPI_FINALIZE(flag, MPIErr)

END SUBROUTINE finalize_MPI

END MODULE Parallel_Mod
