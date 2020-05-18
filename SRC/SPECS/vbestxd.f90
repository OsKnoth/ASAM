!*****************************************************************
!*                                                               *
!*    subroutine  vbestxd                                        *
!*                                                               *
!*****************************************************************
SUBROUTINE vbestxd(vt,kz,dz)

INCLUDE 'HEADER2'

USE data_parallel,      ONLY: my_cart_id

IMPLICIT NONE

INTEGER :: kz
DOUBLE PRECISION :: vt(kz,jmax),dz

INTEGER :: j,k
DOUBLE PRECISION :: rlinks,rrechts,mlinks,mrechts,z,m1
      
IF (my_cart_id == 0) THEN
  OPEN(11,FILE='RESULTATE/vt')
ENDIF
! Berechne Radien und Massen der einzelnen Klassen in Metern bzw. kg.
! Die Fallgeschwindigkeiten werden in m/s berechnet.
!PI = 2.0*asin(1.0)
!m1 = (1000.0D0*4.0D0*PI*(1.0D-9)**(3.0D0))/3.0D0
z = -dz
DO k=1,kz
  z = z + dz
  DO j=1,jmax
!    rlinks = 1.0D-9*2.0D0**((DBLE(FLOAT(j)) - 1.0D0)/3.0D0)
!    rrechts = 1.0D-9*2.0D0**(DBLE(FLOAT(j))/3.0D0)
!    zent(j) = (rlinks + rrechts)/2.0D0
!    mlinks = m1*2.0D0**(DBLE(FLOAT(j)) - 1.0D0)
!    mrechts = m1*2.0D0**(DBLE(FLOAT(j)))
!    mzent(j) = (mlinks + mrechts)/2.0D0
    vt(k,j) = 0.0D0
! Formulierung nach Tsias (1996).
!    IF(j.GE.28) vt(k,j) = -9.32D0*exp(40.5D-6*z)*(1.0D0              &
!       - exp(-(mzent(j)/2.9D-6)**0.382D0))
!
! Formulierung nach Silverman und Glass (1973).
! Diese Formeln benoetigen die Radien zent(j) in Mikrometern, deshalb
! der Faktor 1D6
!    IF(zent(j).LT.25.0D-6) vt(k,j) = -1.136D-4*                      &
!      (zent(j)*1.0D6)**2.0D0*EXP(1.72D-5*z)
!    IF(zent(j).GE.25.0D-6.AND.zent(j).LT.150.0D-6) vt(k,j) =         &
!      -1.88D0*EXP(2.56D-5*z)*(1.0D0 - EXP(-(zent(j)*1.0D6/           &
!      152.0D0)**1.819D0))
!    IF(zent(j).GE.150.0D-6) vt(k,j) = -9.58D0*EXP(3.54D-5*z)*        &
!      (1.0D0 - EXP(-(zent(j)*1.0D6/885.0D0)**1.147))
!
    IF(rmitte(j).LT.25.0D-6) vt(k,j) = -1.136D-4*(rmitte(j)*          &
      1.0D6)**2.0D0*EXP(1.72D-5*z)
    IF(rmitte(j).GE.25.0D-6.AND.rmitte(j).LT.150.0D-6) vt(k,j) =      &
      -1.88D0*EXP(2.56D-5*z)*(1.0D0 - EXP(-(rmitte(j)*1.0D6/          &
      152.0D0)**1.819D0))
    IF(rmitte(j).GE.150.0D-6) vt(k,j) = -9.58D0*EXP(3.54D-5*z)*       &
      (1.0D0 - EXP(-(rmitte(j)*1.0D6/885.0D0)**1.147))
!
IF (my_cart_id == 0) THEN
    WRITE(11,100) rmitte(j),z,vt(k,j)
ENDIF

  END DO

IF (my_cart_id == 0) THEN
  WRITE(11,*)               
ENDIF

END DO      
IF (my_cart_id == 0) THEN
  CLOSE(11)
ENDIF
100  FORMAT(D13.6,4X,D13.6,4X,D13.6)

           
RETURN
END
