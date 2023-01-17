MODULE Part_Mod

USE Kind_Mod
USE Parallel_Mod

CONTAINS

SUBROUTINE set_ranks_automatic(subdom_blockids, nblocks)
  IMPLICIT NONE

  INTEGER, ALLOCATABLE, INTENT(inout) :: subdom_blockids(:)
  INTEGER, INTENT(in) :: nblocks
  INTEGER :: i

  ALLOCATE(subdom_blockids(nblocks))
  DO i = 1, nblocks
    subdom_blockids(i) = MOD(i, NumProcs)
  END DO

END SUBROUTINE set_ranks_automatic

END MODULE Part_Mod
