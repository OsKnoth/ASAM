SUBROUTINE nullen2xd(dnw,dqw,dqws,dqwa,kz,jmax,smax,ipmax)

IMPLICIT NONE

INTEGER :: kz,jmax,smax,ipmax
DOUBLE PRECISION :: dnw(kz,jmax,smax,ipmax),dqw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: dqws(kz,jmax,smax,ipmax),dqwa(kz,jmax,smax,ipmax)


dnw = 0.0D0
dqw = 0.0D0
dqws = 0.0D0
dqwa = 0.0D0


RETURN
END
