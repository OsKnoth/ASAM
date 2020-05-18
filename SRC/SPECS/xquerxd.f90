SUBROUTINE xquerxd(mquer,rquer,squer,mfquer,rfquer,sfquer,nw,qw,qwa,nf,   &
                   qfw,qfa,kz)

INCLUDE 'HEADER3'

IMPLICIT NONE

INTEGER :: kz
DOUBLE PRECISION :: mquer(kz,jmax,smax,ipmax),rquer(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: squer(kz,jmax,smax,ipmax),mfquer(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: rfquer(kz,jmax,smax,ipmax),sfquer(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: nw(kz,jmax,smax,ipmax),qw(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qwa(kz,jmax,smax,ipmax),nf(kz,jmax,smax,ipmax)
DOUBLE PRECISION :: qfw(kz,jmax,smax,ipmax),qfa(kz,jmax,smax,ipmax)

INTEGER :: k,j,i,ip


DO ip=1,ipmax
  DO k=1,kz
    DO i=1,smax
      DO j=1,jmax
        IF(nw(k,j,i,ip).GT.small1) THEN
          mquer(k,j,i,ip) = qw(k,j,i,ip)/nw(k,j,i,ip)
          rquer(k,j,i,ip) = (mquer(k,j,i,ip)/fact)**qu1d3
          squer(k,j,i,ip) = qwa(k,j,i,ip)/nw(k,j,i,ip)
        ELSE
          mquer(k,j,i,ip) = 0.0D0
          rquer(k,j,i,ip) = 0.0D0
          squer(k,j,i,ip) = 0.0D0
        END IF
        IF(nf(k,j,i,ip).GT.small1) THEN
          mfquer(k,j,i,ip) = qfw(k,j,i,ip)/nf(k,j,i,ip)
          rfquer(k,j,i,ip) = (mfquer(k,j,i,ip)/facti)**qu1d3
          sfquer(k,j,i,ip) = qfa(k,j,i,ip)/nf(k,j,i,ip)
        ELSE
          mfquer(k,j,i,ip) = 0.0D0
          rfquer(k,j,i,ip) = 0.0D0
          sfquer(k,j,i,ip) = 0.0D0
        END IF
      END DO
    END DO
  END DO
END DO

!!$OPEN(10,FILE='RESULTATE/mftest66',POSITION='APPEND')
!!$OPEN(20,FILE='RESULTATE/mwtest66',POSITION='APPEND')
!!$OPEN(30,FILE='RESULTATE/mftest62',POSITION='APPEND')
!!$OPEN(40,FILE='RESULTATE/mwtest62',POSITION='APPEND')
!!$WRITE(10,100) (mfquer(k,jmax,1,1)/mgrenz(jmax+1), k=1,kz)
!!$WRITE(20,100) (mquer(k,jmax,1,1)/mgrenz(jmax+1), k=1,kz)
!!$WRITE(30,100) (mfquer(k,62,1,1)/mgrenz(63), k=1,kz)
!!$WRITE(40,100) (mquer(k,62,1,1)/mgrenz(63), k=1,kz)
!!$100 FORMAT(40(2X,E13.6))
!!$CLOSE(10)
!!$CLOSE(20)
!!$CLOSE(30)
!!$CLOSE(40)


RETURN
END

