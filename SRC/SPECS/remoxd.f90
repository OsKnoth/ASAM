SUBROUTINE remoxd(dateim,dateia)

IMPLICIT NONE

CHARACTER :: dateim*100,dateia*100

INTEGER :: i,k,anz
DOUBLE PRECISION, ALLOCATABLE :: druck(:),Temp(:),tau(:)


OPEN(10,FILE=dateim)
READ(10,*) anz
READ(10,*)
ALLOCATE(druck(anz),Temp(anz),tau(anz))

READ(10,*) (druck(i), i=1,anz)
DO k=1,6
  READ(10,*)
END DO
READ(10,*) (Temp(i), i=1,anz)
READ(10,*)
READ(10,*) (tau(i), i=1,anz)
CLOSE(10)

dateia = 'DATEN/profil.dat'
OPEN(20,FILE=dateia)
WRITE(20,*) anz
DO i=1,anz
  WRITE(20,100) druck(i),Temp(i),tau(i)
END DO
CLOSE(20)

100 FORMAT(3(E13.6,4X))


RETURN
END
