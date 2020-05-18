MODULE Ansatz_Mod

  USE DataType_Mod
  USE Aerosol_Mod

  IMPLICIT NONE

  REAL(RealKind), PRIVATE, PARAMETER :: TolEps=1.d-5
  INTEGER, PRIVATE, PARAMETER :: MaxIter=20

CONTAINS

SUBROUTINE fLogLog(f,fP,p,b)

  REAL(RealKind) :: f,fP
  REAL(RealKind) :: p,b
  
  REAL(RealKind) :: N,NP,Z,ZP,ff,ffP
  
  IF (b>Zero) THEN
    ff=p**(-(b+One))
    ffP=-LOG(p)*ff
    N=(p-ff)/(b+2.0d0)
    NP=(-ffp*(b+2.0d0)-(p-ff))/(b+2.0d0)**2
    Z=(One-ff)/(b+One)
    ZP=(-ffP*(b+One)-(One-ff))/(b+One)**2
  ELSE
    ff=p**(b+One)
    ffP=LOG(p)*ff
    IF (b==-One) THEN
      N=(ff*p-One)/(b+2.0d0)
      NP=(ffP*p*(b+2.0d0)-(ff*p-One))/(b+2.0d0)**2
      Z=LOG(p)
      ZP=Zero
    ELSE IF (b==-2.0d0) THEN
      N=LOG(p)
      NP=Zero
      Z=(ff-One)/(b+One)
      ZP=(ffP*(b+One)-(ff-One))/(b+One)**2
    ELSE
      N=(ff*p-One)/(b+2.0d0)
      NP=(ffP*p*(b+2.0d0)-(ff*p-One))/(b+2.0d0)**2
      Z=(ff-One)/(b+One)
      ZP=(ffP*(b+One)-(ff-One))/(b+One)**2
    END IF
  END IF
  f=N/Z
  fP=(NP*Z-N*ZP)/Z**2

END SUBROUTINE fLogLog

SUBROUTINE AnsatzLogLog(N,M,xL,xU,cL,cU)

  REAL(RealKind) :: N,M,xL,xU,cL,cU

  REAL(RealKind) :: A,b,p,bOld,f,fP,xC,xCC,xLL
  REAL(RealKind) :: ppFac,xM,ccc
  INTEGER :: Iter,j

  xCC=M/N
  xCC=MIN(MAX((1.0d0+.001d0)*xL,xCC),(1.0d0-.001d0)*xU)
  IF (xCC<xL.OR.xCC>xU) THEN
    WRITE(*,*) 'OSSI'
  END IF
  xLL=xL
  IF (xLL==Zero) THEN
    xLL=xU/100.0d0 
  END IF
  IF (xLL>Zero) THEN
    xC=xCC/xLL
    p=xU/xLL

    b=Zero
    Iter=0
    DO 
      Iter=Iter+1
      bOld=b
      CALL fLogLog(f,fP,p,bOld)
      b=bOld-(f-xC)/fP
      IF (ABS(b-bOld)<=TolEps*ABS(b)+TolEps) THEN
        EXIT
      END IF
      IF (Iter>5) THEN
        EXIT
      END IF
    END DO
    b=MAX(MIN(b,bMax),-bMax)
    IF (b/=-One) THEN
      A=N/((p**(b+One)-One)/(b+One))/xLL
      cL=MAX(Zero,A)
      cU=MAX(Zero,A*p**b)
    ELSE
      A=N/LOG(xU/xLL)/xLL
      cL=MAX(Zero,A)
      cU=MAX(Zero,A/p)
    END IF
  ELSE
    b=(2.0d0*xCC-xU)/(xU-xCC)
    cL=0.0d0
    cU=N*(b+1.0d0)/xU
  END IF
    
END SUBROUTINE AnsatzLogLog

SUBROUTINE AnsatzLogLogK(N,M,xL,xU,A,b)

  REAL(RealKind) :: N,M,xL,xU,A,b

  REAL(RealKind) :: p,bOld,f,fP,xC,xCC
  REAL(RealKind) :: ppFac,xM,ccc
  INTEGER :: Iter,j
  REAL(RealKind) :: N1,M1
  REAL(RealKind) :: N2,M2

  IF (N<=Zero) THEN
    A=Zero
    b=Zero
  ELSE
    xCC=M/N
    IF (xL>Zero) THEN
      xC=M/N/xL
      p=xU/xL

      b=Zero
      Iter=0
      DO 
        Iter=Iter+1
        bOld=b
        CALL fLogLog(f,fP,p,bOld)
        b=bOld-(f-xC)/fP
        IF (ABS(b-bOld)<=TolEps*ABS(b)+TolEps) THEN
          EXIT
        END IF
        IF (Iter>20) THEN
          EXIT
        END IF
      END DO
      b=MAX(MIN(b,bMax),-bMax)
      IF (b/=-One) THEN
        A=N/((p**(b+One)-One)/(b+One))/xL
      ELSE
        A=N/LOG(xU/xL)/xL
      END IF
    ELSE
      b=(2.0d0*xCC-xU)/(xU-xCC)
      A=N*(b+1.0d0)/xU**(b+1.0d0)
    END IF
  END IF
    
END SUBROUTINE AnsatzLogLogK

SUBROUTINE IntLog(N,M,A,b,xL,xU,xB)
  REAL(RealKind) :: N,M,A,b,xL,xU,xB
  
  IF (xB>Zero) THEN
    N=A/(b+1.0d0)*xB*((xU/xB)**(b+1.0d0)-(xL/xB)**(b+1.0d0))
    M=A/(b+2.0d0)*xB*xB*((xU/xB)**(b+2.0d0)-(xL/xB)**(b+2.0d0))
  ELSE
    N=A/(b+1.0d0)*(xU**(b+1.0d0)-xL**(b+1.0d0))
    M=A/(b+2.0d0)*(xU**(b+2.0d0)-xL**(b+2.0d0))
  END IF

END SUBROUTINE IntLog

SUBROUTINE LimiterCompute(c,c1)

  TYPE(Vec4_T) :: c(:)
  TYPE(Vec4_T) :: c1(:)

  INTEGER :: i,iFrac,iSpc
  INTEGER :: ix,iy,iz
  REAL(RealKind), POINTER :: NN(:)
  REAL(RealKind) :: mGes(nFrac),mGesOld(nFrac)
  REAL(RealKind) :: MassGes,MassGesNew

  IF (SIZE(c)>0) THEN
    mGes=Zero
    mGesOld=Zero
    DO iz=iz0+1,iz1
      DO iy=iy0+1,iy1
        DO ix=ix0+1,ix1
          DO iFrac=1,nFrac
            IF (c(iNC)%c(ix,iy,iz,iFrac)<=NuMin) THEN
              DO iSpc=1,Size(c)
                c(iSpc)%c(ix,iy,iz,iFrac)=Zero
                c1(iSpc)%c(ix,iy,iz,iFrac)=Zero
              END DO
            END IF
          END DO
          DO iSpc=3,SIZE(c)
            DO iFrac=1,nFrac
              IF (c(iSpc)%c(ix,iy,iz,iFrac)<Zero) THEN
                 c(iSpc)%c(ix,iy,iz,iFrac)=Zero
              END IF
            END DO
          END DO
        END DO
      END DO
    END DO
    IF (Limit) THEN
      DO iz=iz0+1,iz1
        DO iy=iy0+1,iy1
          DO ix=ix0+1,ix1
            DO iFrac=1,nFrac
              IF (c(iNC)%c(ix,iy,iz,iFrac)<=NuMin) THEN
                DO iSpc=1,Size(c)
                  c(iSpc)%c(ix,iy,iz,iFrac)=Zero
                  c1(iSpc)%c(ix,iy,iz,iFrac)=Zero
                END DO
              ELSE
                MassGes=Zero
                DO iSpc=3,Size(c)
                  c(iSpc)%c(ix,iy,iz,iFrac)=MAX(c(iSpc)%c(ix,iy,iz,iFrac),Zero)
                  MassGes=MassGes+c(iSpc)%c(ix,iy,iz,iFrac)
                END DO
                MassGes=MassGes/c(1)%c(ix,iy,iz,iFrac)
                IF (MassGes>Zero) THEN
                  mGesOld(iFrac)=MassGes
                  IF (MassGes<=m(iFrac).OR.MassGes>=m(iFrac+1)) THEN
                    MassGesNew=MIN(MAX(m(iFrac)*(1+MassEps),MassGes),m(iFrac+1)*(1-MassEps))
                    DO iSpc=3,Size(c)
                      c(iSpc)%c(ix,iy,iz,iFrac)=c(iSpc)%c(ix,iy,iz,iFrac)*MassGesNew/MassGes
                    END DO
                    mGesOld(iFrac)=MassGesNew
                  END IF
                ELSE
                  mGesOld(iFrac)=Zero
                  DO iSpc=1,Size(c)
                    c(iSpc)%c(ix,iy,iz,iFrac)=Zero
                  END DO
                END IF
              END IF
            END DO
          END DO
        END DO
      END DO
    END IF
  END IF

END SUBROUTINE LimiterCompute 


SUBROUTINE AnsatzHyper(N,M,xLL,xU,cL,cU)

  REAL(RealKind) :: N,M,xLL,xU,cL,cU

  INTEGER :: Iter
  REAL(RealKind) :: xC,alpha,alphaOld,b,xL

  xC=M/N
  xL=xLL
  IF (xLL==Zero) THEN
    xL=xU/100.0d0 
  END IF
  xC=MIN(MAX((1.0d0+.2d0)*xL,xC),(1.0d0-.2d0)*xU)

  IF (xC<0.5d0*(xU+xL)) THEN
    alpha=-0.999d0*xL+1.d-30
  ELSE 
    alpha=-1.001d0*xU
  END IF
  Iter=0
  DO
    alphaOld=alpha  
    alpha=alpha-f(alpha)/fPrime(alpha)
    IF (ABS(alpha-alphaOld)<=TolEps*ABS(alpha)+1.d-60) THEN
      EXIT
    END IF
    Iter=Iter+1
    IF (Iter>20) THEN
      IF (xC<0.5d0*(xU+xL)) THEN
        alpha=-0.999d0*xL+1.d-30
      ELSE
        alpha=-1.001d0*xU
      END IF
      Iter=0
      WRITE(*,*) 'xL',xL
      WRITE(*,*) 'xC',xC
      WRITE(*,*) 'xU',xU
      WRITE(*,*) 'alpha',alpha
      DO
        alphaOld=alpha
        WRITE(*,*) alpha,f(alpha)
        alpha=alpha-f(alpha)/fPrime(alpha)
        IF (ABS(alpha-alphaOld)<=TolEps*ABS(alpha)+1.d-60) THEN
          EXIT
        END IF
        Iter=Iter+1
        IF (Iter>20) THEN
          EXIT
          STOP 'Fehler in AnsatzHyper'
        END IF
      END DO
      EXIT
      STOP 'Fehler in AnsatzHyper'
    END IF
  END DO
  b=LOG((alpha+xU)/(alpha+xL))/N
  cL=1.0d0/(b*(alpha+xL))
  cU=1.0d0/(b*(alpha+xU))

CONTAINS

FUNCTION f(alpha)
  REAL(RealKind) :: f,alpha,temp,templ,tempu
  templ=alpha+xL
  tempu=alpha+xu
  temp=(alpha+xU)/(alpha+xL)
  f=(xC+alpha)*LOG((alpha+xU)/(alpha+xL))-(xU-xL)
END FUNCTION f

FUNCTION fPrime(alpha)
  REAL(RealKind) :: fPrime,alpha
  fPrime=LOG((alpha+xU)/(alpha+xL)) &
        +(xC+alpha)*(1.0d0/(alpha+xU)-1.0d0/(alpha+xL)) 
END FUNCTION fPrime

END SUBROUTINE AnsatzHyper

END MODULE Ansatz_Mod
