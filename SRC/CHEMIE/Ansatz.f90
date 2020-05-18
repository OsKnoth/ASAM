SUBROUTINE AnsatzHyper(N,M,xL,xU,cL,cU)

  IMPLICIT NONE
  INTEGER, PARAMETER :: RealKind=8
  REAL(RealKind) :: N,M,xL,xU,cL,cU

  INTEGER :: Iter
  REAL(RealKind) :: xC,alpha,alphaOld,b
  REAL(RealKind), PARAMETER :: TolEps=1.d-12

  xC=M/N

  WRITE(*,*) xL,xC,xU
  IF (xC<0.5d0*(xU+xL)) THEN
    alpha=-0.99d0*xL
  ELSE 
    alpha=-1.01d0*xU
  END IF
  Iter=0
  DO
    alphaOld=alpha  
    alpha=alpha-f(alpha)/fPrime(alpha)
    WRITE(*,*) alpha,f(alpha)
    IF (ABS(alpha-alphaOld)<=TolEps*ABS(alpha)+TolEps) THEN
      EXIT
    END IF
    IF (Iter>20) THEN
      STOP 'Fehler in AnsatzHyper'
    END IF
  END DO
  b=LOG((alpha+xU)/(alpha+xL))/N
  cL=1.0d0/(b*(alpha+xL))
  cU=1.0d0/(b*(alpha+xU))

  DO iter=0,100
    xC=xL+iter*(xU-Xl)/100.0d0
    WRITE(40,*) xC,1.0d0/(b*(alpha+xC)),N/(xU-xL)
  END DO

  WRITE(*,*) 'N',N,Mom0(xL,xU,alpha,b)
  WRITE(*,*) 'M',M,Mom1(xL,xU,alpha,b)

CONTAINS

FUNCTION f(alpha)
  REAL(RealKind) :: f,alpha
  f=(xC+alpha)*LOG((alpha+xU)/(alpha+xL))-(xU-xL)
END FUNCTION f

FUNCTION fPrime(alpha)
  REAL(RealKind) :: fPrime,alpha
  fPrime=LOG((alpha+xU)/(alpha+xL)) &
        +(xC+alpha)*(1.0d0/(alpha+xU)-1.0d0/(alpha+xL)) 
END FUNCTION fPrime

FUNCTION Mom0(mL,mU,a,b)
  REAL(RealKind) :: Mom0
  REAL(RealKind) :: mL,mU,a,b
  Mom0=LOG((a+mU)/(a+mL))/b
END FUNCTION Mom0

FUNCTION Mom1(mL,mU,a,b)
  REAL(RealKind) :: Mom1
  REAL(RealKind) :: mL,mU,a,b
  Mom1=(mU-mL-a*LOG((a+mU)/(a+mL)))/b
END FUNCTION Mom1

END SUBROUTINE AnsatzHyper

PROGRAM Hyper

  IMPLICIT NONE

  INTEGER, PARAMETER :: RealKind=8

  INTEGER :: i
  REAL(RealKind) :: m,mL,mU
  REAL(RealKind) :: N,cL,cU

  mL=1.0d0
  mU=5.0d0

  WRITE(*,*) 'Eingabe m'
  READ(*,*) m

  N=1

  CALL AnsatzHyper(N,m,mL,mU,cL,cU)
  WRITE(*,*) 'cL',cL
  WRITE(*,*) 'cU',cU

END PROGRAM Hyper
