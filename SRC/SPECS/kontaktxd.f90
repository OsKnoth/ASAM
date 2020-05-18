SUBROUTINE kontaktxd(dateik,miv,itmax)

IMPLICIT NONE
INTEGER :: itmax
DOUBLE PRECISION :: miv(itmax)
CHARACTER :: dateik*100

INTEGER :: it


OPEN(10,FILE=dateik)
READ(10,*) (miv(it), it=1,itmax)
CLOSE(10)


RETURN
END
