MODULE Parallel_Mod
   
   USE mpi_f08
!   USE netcdf

   IMPLICIT NONE


   INTEGER :: ncid

   CHARACTER*128     :: FileNameMeanProfile = ""
   INTEGER :: MyId
   INTEGER :: NumProcs
   INTEGER :: MPIErr

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

!
!---------------------------------------------
END SUBROUTINE start_MPI

END MODULE Parallel_Mod

PROGRAM MainProgMPI_Test

  USE Parallel_Mod

  IMPLICIT NONE


  INTEGER, PARAMETER :: n=10
  CHARACTER(80) :: FileName='AAA'
  REAL(4) :: Val(n)
  INTEGER :: nReals
  INTEGER(KIND=MPI_OFFSET_KIND) :: Offset
  TYPE(MPI_File) :: fh
  TYPE(MPI_Status) :: status

  CALL start_MPI

  
  OffSet=(4*MyId*n)+1
  Val=MyId
  WRITE(*,*) 'Open mpi file'
  CALL MPI_FILE_OPEN(MPI_COMM_WORLD,FileName, &
                     MPI_MODE_CREATE+MPI_MODE_WRONLY, &
                     MPI_INFO_NULL,fh,MPIErr)

  nReals=n
  WRITE(*,*) 'write into the file'
  CALL MPI_FILE_WRITE_AT(fh,Offset,Val,nReals, &
                         MPI_REAL,Status,MPIErr)

  WRITE(*,*) 'Close the file' 
  CALL MPI_FILE_CLOSE(fh,MPIErr)

  WRITE(*,*) 'Finalize the writing process'
  CALL MPI_Finalize(MPIErr)

!  WRITE(*,*) 'netcdf_command'
!  CALL check(nf90_open(TRIM(FileNameMeanProfile)//'.nc',nf90_write,ncid))

END PROGRAM MainProgMPI_Test
