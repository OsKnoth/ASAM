!************************************************************
!*                                                          *
!*    function cvmgmxd                                      *
!*                                                          *
!*    entspricht Gleichung (3.29)                           *
!*                                                          *
!************************************************************
DOUBLE PRECISION FUNCTION cvmgmxd(a,b,c)
  
  IMPLICIT NONE
  DOUBLE PRECISION :: a,b,c
  
  IF(c.LT.0.0D0) cvmgmxd = a
  IF(c.GT.0.0D0) cvmgmxd = b
  IF(c.EQ.0.0D0) cvmgmxd = 0.0D0
  
  RETURN
END FUNCTION cvmgmxd
