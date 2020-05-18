PROGRAM Hyper

  IMPLICIT NONE
  
  INTEGER, PARAMETER :: RealKind=8

  INTEGER :: i
  REAL(RealKind) :: m,mL,mU
  REAL(RealKind) :: alpha,dalpha,func
  REAL(RealKind) :: a=1.0d0,b=1.0d0,da=.01
  REAL(RealKind) :: M0,M1

  
  
  mL=1.0d0
  mU=5.0d0

  WRITE(*,*) 'Eingabe m'
  READ(*,*) m

  dalpha=.01d0
  alpha=-100.0d0
  DO 
    func=f(alpha)
    WRITE(40,*) alpha,func,0.0d0
    alpha=alpha+dalpha
    IF (alpha>=-mU) THEN
      EXIT
    END IF
  END DO
  alpha=100.0d0
  DO 
    func=f(alpha)
    WRITE(41,*) alpha,func,0.0d0
    alpha=alpha-dalpha
    IF (alpha<=-mL) THEN
      EXIT
    END IF
  END DO

  IF (m<0.5d0*(mU+mL)) THEN
    alpha=-0.99d0*mL
  ELSE 
    alpha=-1.01d0*mU
  END IF
  DO  i=1,10
    WRITE(*,*) alpha,f(alpha)
    alpha=alpha-f(alpha)/fPrime(alpha)
  END DO

CONTAINS

FUNCTION f(alpha)
  REAL(RealKind) :: f,alpha
  f=(m+alpha)*LOG((alpha+mU)/(alpha+mL))-(mU-mL)
END FUNCTION f

FUNCTION Mom0(mL,mU,a,b)
  REAL(RealKind) :: Mom0
  REAL(RealKind) :: mL,mU,a,b
  Mom0=LOG((a+b*mU)/(a+b*mL))/b
END FUNCTION Mom0

FUNCTION Mom1(mL,mU,a,b)
  REAL(RealKind) :: Mom1
  REAL(RealKind) :: mL,mU,a,b
  Mom1=(mU-mL-a*LOG((a+b*mU)/(a+b*mL))/b)/b
END FUNCTION Mom1

FUNCTION fPrime(alpha)
  REAL(RealKind) :: fPrime,alpha
  fPrime=LOG((alpha+mU)/(alpha+mL)) &
        +(m+alpha)*(1.0d0/(alpha+mU)-1.0d0/(alpha+mL)) 
END FUNCTION fPrime

END PROGRAM Hyper
