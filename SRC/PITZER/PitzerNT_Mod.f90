MODULE Pitzer_Mod

  IMPLICIT NONE

!---  dimensions
  INTEGER, PARAMETER :: RealKind=8
  INTEGER :: nc = 4, na = 5
  INTEGER :: nac, nacc, naca, nc2, na2
  INTEGER :: zMax=2

!---  INTEGER variable arrays
  INTEGER, ALLOCATABLE :: nzc(:),nza(:)          ! charges
  INTEGER, ALLOCATABLE :: nuec(:,:),nuea(:,:)    ! ion number per salt molecule

!--------------------------------------------------------
!---  REAL(RealKind) constants

      REAL(RealKind), PARAMETER :: cmw=1000.d0/18.015d0


!---  REAL(RealKind) parameter arrays
  REAL(RealKind), PRIVATE, ALLOCATABLE :: b0(:,:),b1(:,:),a1(:,:),        &
                                          b2(:,:),a2(:,:),c(:,:),        &
                                          psim(:,:,:),ThetaX(:,:),        &
                                          psix(:,:,:),ThetaM(:,:)         

!---  REAL(RealKind) variables 
  REAL(RealKind) :: stri,aW
  REAL(RealKind), ALLOCATABLE :: gmx(:,:), gc(:), ga(:)
  REAL(RealKind), ALLOCATABLE :: cma(:), cmc(:)

!--------------------------------------------------------
!---  CHARACTER  arrays with names
  CHARACTER(20), ALLOCATABLE :: mname(:), xname(:)
  CHARACTER(20), ALLOCATABLE :: mxname(:,:)

!--------------------------------------------------------
   INTEGER, ALLOCATABLE :: PitzInd(:)

!  HILFSGROESSEN

   REAL(RealKind), ALLOCATABLE :: Charge(:)
   CHARACTER*20, ALLOCATABLE :: SpeciesNameAqua(:)
   INTEGER :: nAqua

CONTAINS


SUBROUTINE ActivityCoeffCompute(Temp,p)

  REAL(RealKind) :: Temp,p

  INTEGER :: iA,i1A,i2A,iC,i1C,i2C
  INTEGER :: zi,zj
  REAL(RealKind) :: APhi,IonStr,SQRTIonStr,fPrime
  REAL(RealKind) :: Sum,SumZ,SumM
  REAL(RealKind) :: xbar,uBar,uPrimeBar
  REAL(RealKind) :: Temp1,Temp2
  REAL(RealKind) :: gCalc1,DgCalc1DIx
  REAL(RealKind) :: gCalc2,DgCalc2DIx
  REAL(RealKind) :: DSumDIx,DH1,DH2,DDHDIx
  REAL(RealKind) :: gaLoc(na),gcLoc(nc)
  REAL(RealKind) :: gaLoc1(na),gcLoc1(nc)
  REAL(RealKind) :: xij
  REAL(RealKind) :: JCalc(zMax,zMax),JCalcDx(zMax,zMax)
  REAL(RealKind) :: ThetaCalc(zMax,zMax),ThetaCalcDIx(zMax,zMax)
  REAL(RealKind) :: HOE,DHOEDix

  REAL(RealKind) :: t1

  APhi=APhiF(Temp,p)
  IonStr=IonicStrength()
  SQRTIonStr=SQRT(IonStr)
!     --------------------------------------
!    | (1) Charge * mole fraction sums that |
!    | we will need later.                  |
!     --------------------------------------
  SumZ=0.0d0
  SumM=0.0d0
  DO iC=1,nC
    SumZ=SumZ+cmc(iC)*nzc(iC)
    SumM=SumM+cmc(iC)
  END DO
  DO iA=1,nA
    SumZ=SumZ-cma(iA)*nza(iA)
    SumM=SumM+cma(iA)
  END DO

  gc=0.0d0
  ga=0.0d0
  

!
!   #(01)# Calculate Debye-Huckel (long-range) terms.
!          ** Act coeff contribution is +SUMDH ** 
!
  DH1=APhi*fF(IonStr)
  DDHDIx=APhi*fPrimeF(IonStr)
  t1=DH1-DDHDIx*IonStr

  gcLoc=0.0d0
  DH2=0.0d0
  DO iA=1,nA
    Sum=0.0d0
    DSumDIx=0.0d0
    DO iC=1,nC
      IF (a1(iC,iA)>0.0d0) THEN
        xBar=a1(iC,iA)*SQRTIonStr
        Temp1=EXP(-xBar)
        uBar=2.0d0*b1(iC,iA)*(1.0d0-(1.0d0+xBar)*Temp1)
        uPrimeBar=2.0d0*b1(iC,iA)*xBar*Temp1
        gCalc1=2.0d0*(uBar/(xBar*xBar)+b0(iC,iA))
        DgCalc1DIx=2.0d0*(uPrimeBar-2.0d0*uBar/xBar)*(0.5d0*a1(iC,iA)/SQRTIonStr)/(xBar*xBar)
      ELSE 
        gCalc1=2.0d0*(b1(iC,iA)+b0(iC,iA))
        DgCalc1DIx=0.0d0
      END IF
      IF (a2(iC,iA)>0.0d0) THEN
        xBar=a2(iC,iA)*SQRTIonStr
        Temp1=EXP(-xBar)
        uBar=2.0d0*b2(iC,iA)*(1.0d0-(1.0d0+xBar)*Temp1)
        uPrimeBar=2.0d0*b2(iC,iA)*xBar*Temp1
        gCalc2=2.0d0*uBar/(xBar*xBar)
        DgCalc2DIx=2.0d0*(uPrimeBar-2.0d0*uBar/xBar)*(0.5d0*a2(iC,iA)/SQRTIonStr)/(xBar*xBar)
      ELSE
        gCalc2=2.0d0*b2(iC,iA)
        DgCalc2DIx=0.0d0
      END IF
      Sum=Sum &
         +cmc(iC)*(gCalc1+gCalc2)
      gcLoc(iC)=gcLoc(iC) &
         +cma(iA)*(gCalc1+gCalc2)
      DSumDIx=DSumDIx &
         +cmc(iC)*(DgCalc1DIx+DgCalc2DIx)
    END DO
    gaLoc(iA)=Sum
    DH2=DH2+cma(iA)*Sum
    DDHDIx=DDHDIx+cma(iA)*DSumDIx
  END DO  

  aW=DH1-DH2-DDHDIx*IonStr
  DO iC=1,nC
    gc(iC)=gc(iC)-DDHDIx*(-0.5d0*nzc(iC)**2)+gcLoc(iC)
  END DO
  DO iA=1,nA
    ga(iA)=ga(iA)-DDHDIx*(-0.5d0*nza(iA)**2)+gaLoc(iA)
  END DO

  Temp2=0.0d0
  gcLoc=0.0d0
  DO iA=1,nA
    Temp1=0.0d0
    DO iC=1,nC
      Temp1=Temp1+c(iC,iA)*cmc(ic)
      gcLoc(iC)=gcLoc(iC)+c(iC,iA)*cma(iA) 
    END DO
    gaLoc(iA)=Temp1
    Temp2=Temp2+Temp1*cma(iA)
  END DO 
  aW=aW-2.0d0*SumZ*Temp2
  DO iC=1,nC
    gc(iC)=gc(iC)+nzc(ic)*Temp2+SumZ*gcLoc(iC)
  END DO
  DO iA=1,nA
    ga(iA)=ga(iA)-nza(iA)*Temp2+SumZ*gaLoc(iA)
  END DO
!     
!   #(02)# Calculate unsymmetrical mixing terms. 
!
  ThetaCalc=0.0d0
  ThetaCalcDIx=0.0d0
  IF (zMax>1) THEN
    DO zi=1,zMax
      DO zj=zi,zMax
        xij=6.0d0*zi*zj*APhi*SQRTIonStr
        CALL JFunc(xij,JCalc(zi,zj),JCalcDx(zi,zj))
      END DO
    END DO
    DO zi=1,zMax
      DO zj=zi+1,zMax
        xij=6.0d0*zi*zj*Aphi*SQRTIonStr
        ThetaCalc(zi,zj)=((zi*zj)/(4.0d0*IonStr))*(JCalc(zi,zj)-0.5d0*(JCalc(zi,zi)+JCalc(zj,zj)))
        ThetaCalc(zj,zi)=ThetaCalc(zi,zj)
        ThetaCalcDIx(zi,zj)= &
             -ThetaCalc(zi,zj)/IonStr &
             +((zi*zj)/(8.0d0*IonStr*IonStr)) &
             *(zi*zj*JCalcDx(zi,zj)-0.5d0*(zi*zi*JCalcDx(zi,zi)+zj*zj*JCalcDx(zj,zj))) &
             *6.0d0*APhi*SQRTIonStr
        ThetaCalcDIx(zj,zi)=ThetaCalcDIx(zi,zj)
      END DO
    END DO
  END IF
  ThetaCalc=0.0d0
  ThetaCalcDIx=0.0d0

  HOE=0.0d0
  DHOEDIx=0.0d0
  gaLoc=0.0d0
  DO iA=1,nA
    Sum=0.0d0
    DSumDIx=0.0d0
    zi=-nza(iA)
    DO i1A=iA+1,nA
      zj=-nza(i1A)
      Sum=Sum+cma(i1A)*(ThetaX(iA,i1A)+ThetaCalc(zi,zj))
      gaLoc(i1A)=gaLoc(i1A)+cma(iA)*(ThetaX(iA,i1A)+ThetaCalc(zi,zj))
      DSumDIx=DSumDIx+cma(i1A)*ThetaCalcDIx(zi,zj)
    END DO
    HOE=HOE+Sum*cma(iA)
    gaLoc(iA)=gaLoc(iA)+Sum
    DHOEDIx=DHOEDIx+DSumDIx*cma(iA)
  END DO 
  gcLoc=0.0d0
  DO iC=1,nC
    Sum=0.0d0
    DSumDIx=0.0d0
    zi=nzc(iC)
    DO i1C=iC+1,nC
      zj=nzc(i1C)
      Sum=Sum+cmc(i1C)*(ThetaM(iC,i1C)+ThetaCalc(zi,zj))
      gcLoc(i1C)=gcLoc(i1C)+cmc(iC)*(ThetaM(iC,i1C)+ThetaCalc(zi,zj))
      DSumDIx=DSumDIx+cmc(i1C)*ThetaCalcDIx(zi,zj)
    END DO
    HOE=HOE+Sum*cmc(iC)
    gcLoc(iC)=gcLoc(iC)+Sum
    DHOEDIx=DHOEDIx+DSumDIx*cmc(iC)
  END DO


  aW=aW-2.0d0*HOE-2.0d0*DHOEDIx*IonStr
  DO iC=1,nC
    gc(iC)=gc(iC)-DHOEDIx*(-nzc(iC)**2)+2.0d0*gcLoc(iC)
  END DO
  DO iA=1,nA
    ga(iA)=ga(iA)-DHOEDIx*(-nza(iA)**2)+2.0d0*gaLoc(iA)
  END DO
!
!   #(4)# Calculate summation for aac parameters: 
!
  Sum=0.0d0
  gcLoc=0.0d0
  gaLoc=0.0d0
  DO iC=1,nC
    Temp2=0.0d0
    gaLoc1=0.0d0
    DO i2A=1,nA
      Temp1=0.0d0
      DO i1A=i2A+1,nA
        Temp1=Temp1+cma(i1A)*PsiX(i1A,i2A,iC)
        gaLoc1(i1A)=gaLoc1(i1A)+cma(i2A)*PsiX(i1A,i2A,iC)
      END DO
      Temp2=Temp2+cma(i2A)*Temp1
      gaLoc1(i2A)=gaLoc1(i2A)+Temp1
    END DO
    Sum=Sum+cmc(iC)*Temp2
    gaLoc=gaLoc+gaLoc1*cmc(iC)
    gcLoc(iC)=Temp2
  END DO
!
  aW=aW-2.0d0*Sum
  DO iC=1,nC
    gc(iC)=gc(iC)+gcLoc(iC)
  END DO
  DO iA=1,nA
    ga(iA)=ga(iA)+gaLoc(iA)
  END DO

!
!   #(5)# Calculate summation for cca parameters:
!
  Sum=0.0d0
  gcLoc=0.0d0
  gaLoc=0.0d0
  DO iA=1,nA
    Temp2=0.0d0
    gcLoc1=0.0d0
    DO i2C=1,nC
      Temp1=0.0d0
      DO i1C=i2C+1,nC
        Temp1=Temp1+cmc(i1C)*PsiM(i1C,i2C,iA)
        gcLoc1(i1C)=gcLoc1(i1C)+cmc(i2C)*PsiM(i1C,i2C,iA)
      END DO
      Temp2=Temp2+cmc(i2C)*Temp1
      gcLoc1(i2C)=gcLoc1(i2C)+Temp1
    END DO
    Sum=Sum+cma(iA)*Temp2
    gcLoc=gcLoc+gcLoc1*cma(iA)
    gaLoc(iA)=Temp2
  END DO
!
  aW=aW-2.0d0*Sum
  DO iC=1,nC
    gc(iC)=gc(iC)+gcLoc(iC)
  END DO
  DO iA=1,nA
    ga(iA)=ga(iA)+gaLoc(iA)
  END DO


  aW=EXP((aW-SumM)/cmw)

  Temp2=0.0d0


CONTAINS

FUNCTION fF(IonStr)

  REAL(RealKind) :: fF
  REAL(RealKind) :: IonStr

  REAL(RealKind), PARAMETER :: b=1.2d0

  fF=-(4.0d0*IonStr/b)*LOG(1.0d0+b*SQRT(IonStr))

END FUNCTION fF

FUNCTION fPrimeF(IonStr)

  REAL(RealKind) :: fPrimeF
  REAL(RealKind) :: IonStr

  REAL(RealKind) :: Work
  REAL(RealKind), PARAMETER :: b=1.2d0

  Work=1.0d0+b*SQRT(IonStr)
  fPrimeF=-2.0d0*(SQRT(IonStr)/Work+2.0d0*LOG(Work)/b)

END FUNCTION fPrimeF

FUNCTION gF(x)
 
  REAL(RealKind) :: gF
  REAL(RealKind) :: x

  gF=2.0d0*(1.0d0-(1.0d0+x)*EXP(-x))/(x*x)

END FUNCTION gF

FUNCTION gPrimeF(x)
 
  REAL(RealKind) :: gPrimeF
  REAL(RealKind) :: x

  gPrimeF=-2.0d0*(1.0d0-(1.0d0+x+0.5d0*x*x)*EXP(-x))/(x*x)

END FUNCTION gPrimeF

SUBROUTINE  JFunc(x,J,Jprime)

  REAL(RealKind) :: x,J,Jprime

  INTEGER :: k
  REAL(RealKind) :: z,dzdx
  REAL(RealKind) :: bk(0:22),dk(0:22) 

  REAL(RealKind) :: ak1(0:20)=     &  
               (/1.925154014814667d0,-0.060076477753119d0,-0.029779077456514d0, &
              -0.007299499690937d0, 0.000388260636404d0, 0.000636874599598d0, &
               0.000036583601823d0,-0.000045036975204d0,-0.000004537895710d0, &
               0.000002937706971d0, 0.000000396566462d0,-0.000000202099617d0, &
              -0.000000025267769d0, 0.000000013522610d0, 0.000000001229405d0, &
              -0.000000000821969d0,-0.000000000050847d0, 0.000000000046333d0, &
               0.000000000001943d0,-0.000000000002563d0,-0.000000000010991d0/)
  REAL(RealKind) :: ak2(0:20)=     &  
               (/0.628023320520852d0, 0.462762985338493d0, 0.150044637187895d0, &
              -0.028796057604906d0,-0.036552745910311d0,-0.001668087945272d0, &
               0.006519840398744d0, 0.001130378079086d0,-0.000887171310131d0, &
              -0.000242107641309d0, 0.000087294451594d0, 0.000034682122751d0, &
              -0.000004583768938d0,-0.000003548684306d0,-0.000000250453880d0, &
               0.000000216991779d0, 0.000000080779570d0, 0.000000004558555d0, &
              -0.000000006944757d0,-0.000000002849257d0, 0.000000000237816d0/)
  
  bk=0.d0
  dk=0.0d0
  IF (x <=1.0d0) THEN
    z=4.d0*x**(0.2d0)-2.d0
    dzdx=0.8d0*x**(-0.8d0)
    DO k=20,0,-1
      bk(k)=z*bk(k+1)-bk(k+2)+ak1(k)
      dk(k)=bk(k+1)+z*dk(k+1)-dk(k+2)
    END DO
  ELSE
    z=40.d0/9.d0*x**(-0.1d0)-22.d0/9.d0
    dzdx=-4.d0/9.d0*x**(-1.1d0)
    DO k=20,0,-1
      bk(k)=z*bk(k+1)-bk(k+2)+ak2(k)
      dk(k)=bk(k+1)+z*dk(k+1)-dk(k+2)
    END DO
  END IF
  J=0.25d0*x-1.d0+0.5d0*(bk(0)-bk(2))
  JPrime=0.25d0+0.5d0*dzdx*(dk(0)-dk(2))

END SUBROUTINE  JFunc


FUNCTION APhiF(Temp,p)

! APhiF calculates the Debye-Hueckel parameter aPhi

  REAL(RealKind) :: APhiF
  REAL(RealKind) :: Temp,p

  REAL(RealKind) ::  Pi,pp,b,c1,eps

  REAL(RealKind), PARAMETER :: u1=3.4279d2, u2=-5.0866d-3, u3=9.46900d-7  &
                              ,u4=-2.0525d0, u5=3.1159d3, u6=-1.8289d+2, u7=-8.0325d+3  &
                              ,u8=4.21142d6, u9=2.1417d0
  REAL(RealKind), PARAMETER :: t0 = 273.15d0
  REAL(RealKind), PARAMETER :: densw = 1.d0,            & ! water density [g/cm**3]
                               avo   = 6.02257d23,    & ! avogadro's number [molec./mol]
                               bol   = 1.3804d-16,    & ! bolzmann's constant [erg/K]
                               e     = 4.80298d-10      ! electronic charge [esu]
  Pi=4.d0*ATAN(1.d0)
  pp=p*1.d-3                                   ! [mbar] -> [bar]
  b=u7+u8/Temp+u9*Temp
  c1=u4+u5/(u6+Temp)
  eps=u1*EXP(u2*Temp+u3*Temp*Temp)
  eps=eps+c1*LOG((b+pp)/(b+1.d3))

  APhiF=1.d0/3.d0*SQRT(2.0d0*pi*avo*densw*1.d-3)*(e*e/eps/bol/Temp)**(1.5d0)
END FUNCTION APhiF

FUNCTION IonicStrength( )

  REAL(RealKind) :: IonicStrength

! calculation of the total ionic strenght

  INTEGER :: ia, ic

  IonicStrength=0.0d0
  DO ia=1,na
    IonicStrength=IonicStrength+cma(ia)*nza(ia)*nza(ia)
  END DO
  DO ic=1,nc
    IonicStrength=IonicStrength+cmc(ic)*nzc(ic)*nzc(ic)
  END DO
  IonicStrength=0.5d0*IonicStrength

END FUNCTION IonicStrength

END SUBROUTINE ActivityCoeffCompute


! ----------------------------------------------------------------------
SUBROUTINE eteta (a,stri1,nzi,nzj,eThetaMn,etetaprmn)
! ----------------------------------------------------------------------
!
!     calculation of additional terms for unsymmetrical mixing:

!     considers the electrostatic forces between like signed ions but
!     unlike charges (mx & nx or mx & my)
!
!     ref: pitzer, k.s.: ion interaction approach: theory and data
!     correlation, in electrolyte solutions..., pp. 76-147, chap. 3
!
! ----------------------------------------------------------------------
      
!---  external variables
      INTEGER :: nzi,nzj
      REAL(RealKind) :: a,stri1,eThetaMn,etetaprmn

!---  internal variables
      INTEGER :: k, n
      REAL(RealKind) :: xj, zj, dzjdxj

      REAL(RealKind) :: ak(0:20,2),bk(0:22),dk(0:22),xx(3),xjxj(3),xjprxj(3)
      DATA ak /    &  ! ak(1)
&               1.925154014814667d0,-0.060076477753119d0,-0.029779077456514d0, &
&              -0.007299499690937d0, 0.000388260636404d0, 0.000636874599598d0, &
&               0.000036583601823d0,-0.000045036975204d0,-0.000004537895710d0, &
&               0.000002937706971d0, 0.000000396566462d0,-0.000000202099617d0, &
&              -0.000000025267769d0, 0.000000013522610d0, 0.000000001229405d0, &
&              -0.000000000821969d0,-0.000000000050847d0, 0.000000000046333d0, &
&               0.000000000001943d0,-0.000000000002563d0,-0.000000000010991d0, &
!                ! ak(2)
&               0.628023320520852d0, 0.462762985338493d0, 0.150044637187895d0, &
&              -0.028796057604906d0,-0.036552745910311d0,-0.001668087945272d0, &
&               0.006519840398744d0, 0.001130378079086d0,-0.000887171310131d0, &
&              -0.000242107641309d0, 0.000087294451594d0, 0.000034682122751d0, &
&              -0.000004583768938d0,-0.000003548684306d0,-0.000000250453880d0, &
&               0.000000216991779d0, 0.000000080779570d0, 0.000000004558555d0, &
&              -0.000000006944757d0,-0.000000002849257d0, 0.000000000237816d0/
      SAVE ak
! ----------------------------------------------------------------------
!
      bk = 0.d0
      dk = 0.d0
!
      xx(1) = 6.d0*nzi*nzi*a*sqrt(stri1)
      xx(2) = 6.d0*nzi*nzj*a*sqrt(stri1)
      xx(3) = 6.d0*nzj*nzj*a*sqrt(stri1)
!
!---  evaluate j-integrals
      DO n=1,3
        xj=xx(n)
        IF (xj > 1.e-30) THEN  ! xj>0.
          IF (xj <= 1.) THEN
            zj=4.d0*xj**(0.2d0)-2.d0
            dzjdxj=0.8d0*xj**(-0.8d0)
            DO k=20,0,-1
              bk(k)=zj*bk(k+1)-bk(k+2)+ak(k,1)
              dk(k)=bk(k+1)+zj*dk(k+1)-dk(k+2)
            END DO
          ELSE
            zj=40.d0/9.d0*xj**(-0.1d0)-22.d0/9.d0
            dzjdxj=-4.d0/9.d0*xj**(-1.1d0)
            DO k=20,0,-1
              bk(k)=zj*bk(k+1)-bk(k+2)+ak(k,2)
              dk(k)=bk(k+1)+zj*dk(k+1)-dk(k+2)             
            END DO
          END IF
          xjxj(n)=0.25d0*xj-1.d0+0.5d0*(bk(0)-bk(2))
          xjprxj(n)=0.25d0+0.5d0*dzjdxj*(dk(0)-dk(2))
        ELSE             
          xjxj(n)=0.d0
          xjprxj(n)=0.d0
        END IF     
      END DO                    ! n-loop
!     
      eThetaMn   = FLOAT(nzi*nzj)/4.d0/stri1*(xjxj(2)-0.5d0*(xjxj(1) + xjxj(3)))
      etetaprmn = -eThetaMn/stri1+float(nzi*nzj)/8.d0/stri1/stri1              &
&               * (xx(2)*xjprxj(2)-0.5d0*(xx(1)*xjprxj(1)+xx(3)*xjprxj(3)))
!     
! ----------------------------------------------------------------------
END SUBROUTINE eteta

SUBROUTINE InitPitzer(nPitz)
! ------------------------------------------------------------------------
!
!    nPitz - flag for the input data set;
!            1:  Judit
!            2:  harvie (1984)
!            3:  pitzer (1973)
!
! -------------------------------------------------------------------------

  INTEGER :: nPitz

!---  internal variables
  INTEGER :: ia, ic, ios, jt, indspc, na1, nc1, pos
  CHARACTER(20):: string

!---  external function
      INTEGER :: ifind
      EXTERNAL   ifind

! -------------------------------------------------------------------------
! ---  Define set of Pitzer coefficients
! -------------------------------------------------------------------------
!
!---  Define index array PitzInd, array for activity coefficients
      ALLOCATE (PitzInd(nAqua))

      PitzInd(:)      = 0
!
!---  Read Pitzer tables
      WRITE(*,801)
      call init_act(nPitz)

! -------------------------------------------------------------------------
!---  Set index transformation 
! -------------------------------------------------------------------------
!
      na1 = 0
      nc1 = 0
!
      PitzInd(:) = 0 
      DO jt=1,nc
         indspc=Position(mname(jt))
         IF (indspc > 0) THEN
            PitzInd(indspc)=jt
         END IF
      END DO
      DO jt=1,na
         indspc=Position(xname(jt))
         IF (indspc>0) THEN
            PitzInd(indspc)=-jt
         END IF
      END DO

! -------------------------------------------------------------------------
!---  Print Pitzer system
! -------------------------------------------------------------------------
!
      WRITE(*,802) nc1,na1, '  No.','  PitzInd','   Charge','     Species Name'
      DO jt=1,nAqua
         IF (PitzInd(jt) <= 0)  CYCLE
         WRITE(*,803)  jt, PitzInd(jt), Charge(jt), SpeciesNameAqua(jt) 
      END DO
      DO jt=1,nAqua
         IF (PitzInd(jt) >= 0)  CYCLE
         WRITE(*,803)  jt, PitzInd(jt), Charge(jt), SpeciesNameAqua(jt) 
      END DO
      WRITE(*,804)

! -------------------------------------------------------------------------
801   FORMAT(1x/1x,75('=')/ ' Pitzer Initialization:' )
802   FORMAT(1x/1x,75('-')/           &
&            ' Considered Pitzer System: Cations =',i3,'     Anions =',i3 // &
&            a5,a9,a9,a18 / 1x,40('-'))
803   FORMAT(i4,i9,f9.1,7x,a18)
804   FORMAT(1x/1x,75('=')/1x)
!
! -------------------------------------------------------------------------
END SUBROUTINE InitPitzer

! ------------------------------------------------------------------------
SUBROUTINE init_act (nPitz)
! ------------------------------------------------------------------------
!
!    nPitz - flag for the input data set;
!            1:  Judit
!            2:  harvie (1984)
!            3:  pitzer (1973)
!
! -------------------------------------------------------------------------
!---  exinternal variables
      INTEGER :: nPitz

!---  internal variables
      INTEGER :: ia, ic, ios
      INTEGER :: ir0=9, ir1=10, ir2=11, ir3=12

      CHARACTER(7)  ::  set
      CHARACTER(12) ::  Path
!
! -------------------------------------------------------------------------
! ---  Define set of Pitzer coefficients
! -------------------------------------------------------------------------
!
      IF (nPitz <= 0) RETURN
      IF (nPitz == 1) THEN
         write (6,'(/a/)') 'use standard pitzer input data...'
         path = 'Data_IfT/'
         set  = '.judit'
      ELSE IF (nPitz == 2) THEN
         write (6,'(/a/)') 'use harvie & standard pitzer input data...'
         path = 'Data_Harvie/'
         set  = '.harvie'
      ELSE IF (nPitz == 3) THEN
         write (6,'(/a/)') 'use standard pitzer input data...'
         path = 'Data_Pitzer/'
         set  = '.pitzer'
      END IF
!
!------------------------------------------------------------------
! --- Read Pitzer System: Dimensions, Names, Ions, Charges
!------------------------------------------------------------------
!
      OPEN (ir0,FILE=TRIM(path)//'name'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACT...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'name'//'.dat'
         STOP  ' INIT_ACT...Error: Check Pitzer =Name= File !!'
      END IF

!---  read and set dimensions
      READ(ir0,*)
      READ(ir0,*)  nc, na

      nac  = na*nc
      nacc = na * nac
      naca = nc * nac
      nc2  = nc * nc
      na2  = na * na

!---  allocate arrays
      ALLOCATE ( nzc(nc),nza(na),nuec(nc,na),nuea(nc,na) )
      nza(:) = 0
      nzc(:) = 0
      nuea(:,:) = 0
      nuec(:,:) = 0

      ALLOCATE ( b0(nc,na),c(nc,na) )
      ALLOCATE ( b1(nc,na),a1(nc,na) )
      ALLOCATE ( b2(nc,na),a2(nc,na) )
      ALLOCATE ( psim(nc,nc,na),ThetaX(na,na) )
      ALLOCATE ( psix(na,na,nc),ThetaM(nc,nc) )
      b0(:,:) = 0.d0
      a1(:,:) = 2.0d0
      b1(:,:) = 0.d0
      a2(:,:) = 2.0d0
      b2(:,:) = 0.d0
      c(:,:)  = 0.d0
      psim(:,:,:) = 0.d0
      psix(:,:,:) = 0.d0
      ThetaM(:,:)  = 0.d0
      ThetaX(:,:)  = 0.d0

      ALLOCATE ( gmx(nc,na),gc(nc),ga(na) )
      ALLOCATE ( cmc(nc),cma(na) )
      gc(:)  = 0.d0
      ga(:)  = 0.d0
      cma(:) = 0.d0
      cmc(:) = 0.d0
      gmx(:,:) = 0.d0

      ALLOCATE ( mname(nc),xname(na) )
      ALLOCATE ( mxname(nc,na) )

!---  read names
      REWIND (ir0)
      CALL read_name(ir0)
      CLOSE (ir0)
!
!------------------------------------------------------------------
! --- determine the stoichiometric coefficients
      DO ic=1,nc
         DO ia=1,na
            nuec(ic,ia) = ABS(nza(ia))
            nuea(ic,ia) = nzc(ic)
         END DO
      END DO
!
!------------------------------------------------------------------
! --- Read ion interaction characteristics: Pitzer coefficients
!------------------------------------------------------------------
!
! ---  B0, B1, C
      OPEN (ir1,FILE=TRIM(path)//'b01c'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACT...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'b01c'//'.dat'
         STOP  ' INIT_ACT...Error: Check Pitzer =B01C= File !!'
      END IF

!---  read values
      CALL read_b01c(ir1)
      CLOSE (ir1)
!
!------------------------------------------------------------------
! --- convert c's according to their charge
      DO ic=1,nc
         DO ia=1,na
            c(ic,ia) = 0.5d0*c(ic,ia)/sqrt(float(abs(nza(ia)*nzc(ic))))
         END DO
      END DO
!
!------------------------------------------------------------------
! ---  TetaM, TetaX
      OPEN (ir2,FILE=TRIM(path)//'teta'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACT...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'teta'//'.dat'
         STOP  ' INIT_ACT...Error: Check Pitzer =Teta= File !!'
      END IF

      CALL read_teta(ir2)
      CLOSE (ir2)
!
!------------------------------------------------------------------
! ---  PsiM, PsiX
      OPEN (ir3,FILE=TRIM(path)//'psinew'//'.dat',STATUS='old',iostat=ios)
      IF (ios /= 0) THEN
         PRINT *,' INIT_ACT...Error: IO-Stat = ',ios,' !'
         PRINT *,'     Check File: ',TRIM(path)//'psinew'//'.dat'
         STOP  ' INIT_ACT...Error: Check Pitzer =Psi= File !!'
      END IF

      CALL read_psi1(ir3)
      CLOSE (ir3)
!
! -------------------------------------------------------------------------
END SUBROUTINE init_act

! ----------------------------------------------------------------------
SUBROUTINE read_psi (iread)
! ----------------------------------------------------------------------
!--- external variables
      INTEGER :: iread
!
!--- internal variables
      INTEGER :: ia, ia1, ia2, ic, ic1, ic2
! ----------------------------------------------------------------------
!
      read (iread,*)
      DO ic1=1,nc
      DO ic2=1,nc
        read (iread,*) (psim(ic1,ic2,ia),ia=1,na)
      END DO
      END DO
      read (iread,*)
      read (iread,*)
      DO ia1=1,na
      DO ia2=1,na
        read (iread,*) (psix(ia1,ia2,ic),ic=1,nc)
      END DO
      END DO
!
      close (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_psi 
!
! ======================================================================
SUBROUTINE read_psi1 (iread)
! ----------------------------------------------------------------------
!--- external variables
      INTEGER :: iread
!
!--- internal variables
      INTEGER :: ia, ia1, ia2, ic, ic1, ic2
! ----------------------------------------------------------------------
!
      WRITE(*,801)  nc,mname
      WRITE(*,802)  na,xname
801   FORMAT(' Cations:',I3,/(11x,5a14))
802   FORMAT(' Anions :',I3,/(11x,5a14))
     
      READ (iread,*)
      DO ic1=1,nc-1
        DO ic2=ic1+1,nc
          READ (iread,*) (psim(ic1,ic2,ia),ia=1,na)
          DO ia=1,na
            psim(ic2,ic1,ia) = psim(ic1,ic2,ia)
          END DO 
        END DO
      END DO
      READ (iread,*)
      READ (iread,*)
      DO ia1=1,na-1
        DO ia2=ia1+1,na
          READ (iread,*) (psix(ia1,ia2,ic),ic=1,nc)
          DO ic=1,nc
            psix(ia2,ia1,ic) = psix(ia1,ia2,ic)
          END DO 
        END DO
      END DO
!
      CLOSE (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_psi1
! ----------------------------------------------------------------------
SUBROUTINE read_name (iread)
! ----------------------------------------------------------------------
!---  external variables
      INTEGER :: iread
!
!---  internal variables
      INTEGER :: ia, ic, nc1, na1
! ----------------------------------------------------------------------

      REWIND (iread)
      READ (iread,*)
      READ (iread,*) nc1,na1
      IF (nc1 /= nc .OR. na1 /= na) THEN
        print*,'dimensions of the parameter space does not coincident'  &
&              //' with the model dimenssions!'
        print*,'nc1: ',nc1,'  nc: ',nc
        print*,'na1: ',na1,'  na: ',na
        stop 'READ_nmae...'
      END IF

      READ (iread,*)
      READ (iread,*) (mname(ic),ic=1,nc)
      READ (iread,*) (nzc(ic),ic=1,nc)
      READ (iread,*)
      READ (iread,*) (xname(ia),ia=1,na)
      READ (iread,*) (nza(ia),ia=1,na)
      READ (iread,*)
      DO ic=1,nc
        READ (iread,*) (mxname(ic,ia),ia=1,na)
      END DO
!
      CLOSE (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_name
! ----------------------------------------------------------------------
SUBROUTINE read_b01c (iread)
! ----------------------------------------------------------------------
!---  external variables
      INTEGER :: iread
!
!---  internal variables
      INTEGER :: ia, ic

! ----------------------------------------------------------------------
      READ (iread,*)
      DO ic=1,nc
         READ (iread,*) (b0(ic,ia),ia=1,na)
      END DO
      READ (iread,*)
      READ (iread,*)
      DO ic=1,nc
         READ (iread,*) (b1(ic,ia),ia=1,na)
      END DO
      READ (iread,*)
      READ (iread,*)
      DO ic=1,nc
         READ (iread,*) (c(ic,ia),ia=1,na)
      END DO
!
      CLOSE (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_b01c 
! ----------------------------------------------------------------------
SUBROUTINE read_teta (iread)
! ----------------------------------------------------------------------
!---  external variables
      INTEGER :: iread

!---  internal variables
      INTEGER :: ia1, ic1, ia2, ic2
! ----------------------------------------------------------------------
!
      READ (iread,*)
      DO ic1=1,nc
        READ (iread,*) (ThetaM(ic1,ic2),ic2=1,nc)
      END DO
      READ (iread,*)
!
      READ (iread,*)
      DO ia1=1,na
        READ (iread,*) (ThetaX(ia1,ia2),ia2=1,na)
      END DO
!
      CLOSE (iread)
! ----------------------------------------------------------------------
END SUBROUTINE read_teta

FUNCTION Position(Species)

  INTEGER :: Position
  CHARACTER(*) :: Species

  INTEGER :: i

  Position=0
  DO i=1,nAqua
    IF (TRIM(Species)==TRIM(SpeciesNameAqua(i))) THEN
      Position=i
      EXIT
    END IF
  END DO

END FUNCTION

END MODULE Pitzer_Mod

PROGRAM TPitzer

  USE Pitzer_Mod
  
  IMPLICIT NONE

  INTEGER :: i
  REAL(RealKind) :: SumA,SumC
  REAL(RealKind) :: TempLoc,Pre
  REAL(RealKind), ALLOCATABLE :: xC(:)

  nAqua=17
  ALLOCATE(Charge(nAqua))
  ALLOCATE(SpeciesNameAqua(nAqua))
  SpeciesNameAqua(1)='NAp'
  Charge(1)=1
  SpeciesNameAqua(2)='NH4p'
  Charge(2)=1
  SpeciesNameAqua(3)='Hp'
  Charge(3)=1
  SpeciesNameAqua(4)='Kp'
  Charge(4)=1
  SpeciesNameAqua(5)='CApp'
  Charge(5)=2
  SpeciesNameAqua(6)='MGpp'
  Charge(6)=2
  SpeciesNameAqua(7)='FEpp'
  Charge(7)=2
  SpeciesNameAqua(8)='MNpp'
  Charge(8)=2
  SpeciesNameAqua(9)='CUpp'
  Charge(9)=2
  SpeciesNameAqua(10)='CLm'
  Charge(10)=-1
  SpeciesNameAqua(11)='SO4mm'
  Charge(11)=-2
  SpeciesNameAqua(12)='OHm'
  Charge(12)=-1
  SpeciesNameAqua(13)='NO3m'
  Charge(13)=-1
  SpeciesNameAqua(14)='HSO4m'
  Charge(14)=-1
  SpeciesNameAqua(15)='HCO3m'
  Charge(15)=-1
  SpeciesNameAqua(16)='CO3mm'
  Charge(16)=-2
  SpeciesNameAqua(17)='BRm'
  Charge(17)=-1
 
  CALL InitPitzer(1)
  ALLOCATE(xC(nAqua))
  SumA=0.0d0
  SumC=0.0d0
  DO i=1,nA
    xC(i)=i
    SumA=SumA+Charge(i)*xC(i)
  END DO
  DO i=1,nA
    xC(i)=xC(i)/SumA
  END DO
  DO i=nA+1,nAqua
    xC(i)=i
    SumC=SumC-Charge(i)*xC(i)
  END DO
  DO i=nA+1,nAqua
    xC(i)=xC(i)/SumC
  END DO
  SumA=0.0d0
  DO i=1,nAqua
    SumA=SumA+Charge(i)*xC(i)
  END DO

  DO i=1,nAqua
    IF (PitzInd(i)>0) THEN
      cmc(PitzInd(i))=xC(i)
    ELSE IF (PitzInd(i)<0) THEN
      cma(-PitzInd(i))=xC(i)
    END IF
  END DO
  WRITE(*,*) cma
  WRITE(*,*) cmc
  gc(:)=0.d0
  ga(:)=0.d0
  gmx(:,:)=0.d0
  TempLoc=289.15D0
  Pre=1000.0d0
  CALL ActivityCoeffCompute(TempLoc,Pre)
! CALL act_coeff (TempLoc,Pre)
! WRITE(*,*) 'gm',gm
! WRITE(*,*) 'gx',gx
! WRITE(*,*) 'aw',aw
  
END PROGRAM TPitzer
