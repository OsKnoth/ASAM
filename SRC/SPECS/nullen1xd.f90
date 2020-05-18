SUBROUTINE nullen1xd(dw,dteta,dq,kz,ipmax)

IMPLICIT NONE

INTEGER :: kz,ipmax
DOUBLE PRECISION :: dw(kz,ipmax),dteta(kz,ipmax),dq(kz,ipmax)


dw = 0.0D0
dteta = 0.0D0
dq = 0.0D0


RETURN
END
