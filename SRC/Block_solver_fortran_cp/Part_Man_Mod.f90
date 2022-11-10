MODULE Part_Man_Mod

USE Kind_Mod

CONTAINS

SUBROUTINE set_ranks(subdom_blockids)
  !This function is only for testing the parallel
  !implementation with n cores. This has to be substituted 
  !by a partition algorithm.

  INTEGER, ALLOCATABLE, INTENT(inout) :: subdom_blockids(:)
  INTEGER :: i, n

  n = 100
  ALLOCATE(subdom_blockids(n))

! subdom_blockids(1) = 0
! subdom_blockids(2) = 0
! subdom_blockids(3) = 0
! subdom_blockids(4) = 0
! subdom_blockids(5) = 1
! subdom_blockids(6) = 0
! subdom_blockids(7) = 0
! subdom_blockids(8) = 0
! subdom_blockids(9) = 0
 
! n = 100

! ALLOCATE(subdom_blockids(n))

!  subdom_blockids(1) = 0  
!  subdom_blockids(2) = 0   
!  subdom_blockids(5) = 0 
!  subdom_blockids(6) = 0 
!  subdom_blockids(3) = 1  
!  subdom_blockids(4) = 1 
!  subdom_blockids(7) = 1 
!  subdom_blockids(8) = 1 
!  subdom_blockids(9) = 2
!  subdom_blockids(10) = 2
!  subdom_blockids(13) = 2 
!  subdom_blockids(14) = 2
!  subdom_blockids(11) = 3
!  subdom_blockids(12) = 3
!  subdom_blockids(15) = 3
!  subdom_blockids(16) = 3

  DO i = 1, n
    subdom_blockids(i) = MOD(i, 6)
  END DO


END SUBROUTINE set_ranks

END MODULE Part_Man_Mod
