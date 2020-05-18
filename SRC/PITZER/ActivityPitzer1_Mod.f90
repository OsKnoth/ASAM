MODULE ActivityPitzer1_Mod

  USE Control_Mod
  USE Chemie_Mod
  USE Aerosol_Mod
  USE Microphysics_Mod

  IMPLICIT NONE

!---  dimensions
  INTEGER :: nA,nC,nN
  INTEGER :: zMax=2


  INTEGER, ALLOCATABLE :: ChargeCat(:),ChargeAn(:)          ! charges
  CHARACTER*20, ALLOCATABLE :: NameCat(:),NameAn(:),NameNeut(:)

  REAL(RealKind), PARAMETER :: cmw=1000.d0/18.015d0

  REAL(RealKind), PRIVATE, ALLOCATABLE :: b0(:,:),b1(:,:),a1(:,:),        &
                                          b2(:,:),a2(:,:),c(:,:),        &
                                          psim(:,:,:),ThetaX(:,:),        &
                                          psix(:,:,:),ThetaM(:,:)         

!---  REAL(RealKind) variables 
  REAL(RealKind) :: stri,aW
  REAL(RealKind), ALLOCATABLE :: ActCat(:), ActAn(:)
  REAL(RealKind), ALLOCATABLE :: cAn(:), cCat(:)

!--------------------------------------------------------

!--------------------------------------------------------
  INTEGER, ALLOCATABLE :: PitzInd(:)

  TYPE Coefficient_T
    CHARACTER*20 :: Type
    INTEGER :: Ind(4)
    REAL(RealKind), POINTER :: Constants(:)
    TYPE(Reaction_T), POINTER :: Next=>NULL()
  END TYPE Coefficient_T
  TYPE ListCoefficient_T
    TYPE(Coefficient_T), POINTER :: Start=>NULL()
    TYPE(Coefficient_T), POINTER :: End=>NULL()
    INTEGER :: LenList=0
  END TYPE ListCoefficient_T

  TYPE(ListCoefficient_T), SAVE :: b0List


CONTAINS


SUBROUTINE ActivityCoeffCompute(cAn,cCat,cNeut, &
                                ActAn,ActCat,ActNeut, &
                                Temp,p)

  REAL(RealKind) :: cAn(:),cCat(:),cNeut(:)
  REAL(RealKind) :: ActAn(:),ActCat(:),ActNeut(:)
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
  REAL(RealKind) :: ActAnLoc(na),ActCatLoc(nc)
  REAL(RealKind) :: ActAnLoc1(na),ActCatLoc1(nc)
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
    SumZ=SumZ+cCat(iC)*ChargeCat(iC)
    SumM=SumM+cCat(iC)
  END DO
  DO iA=1,nA
    SumZ=SumZ-cAn(iA)*ChargeAn(iA)
    SumM=SumM+cAn(iA)
  END DO

  ActCat=0.0d0
  ActAn=0.0d0
  

!
!   #(01)# Calculate Debye-Huckel (long-range) terms.
!          ** Act coeff contribution is +SUMDH ** 
!
  DH1=APhi*fF(IonStr)
  DDHDIx=APhi*fPrimeF(IonStr)
  t1=DH1-DDHDIx*IonStr

  ActCatLoc=0.0d0
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
         +cCat(iC)*(gCalc1+gCalc2)
      ActCatLoc(iC)=ActCatLoc(iC) &
         +cAn(iA)*(gCalc1+gCalc2)
      DSumDIx=DSumDIx &
         +cCat(iC)*(DgCalc1DIx+DgCalc2DIx)
    END DO
    ActAnLoc(iA)=Sum
    DH2=DH2+cAn(iA)*Sum
    DDHDIx=DDHDIx+cAn(iA)*DSumDIx
  END DO  

  aW=DH1-DH2-DDHDIx*IonStr
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)-DDHDIx*(-0.5d0*ChargeCat(iC)**2)+ActCatLoc(iC)
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)-DDHDIx*(-0.5d0*ChargeAn(iA)**2)+ActAnLoc(iA)
  END DO

  Temp2=0.0d0
  ActCatLoc=0.0d0
  DO iA=1,nA
    Temp1=0.0d0
    DO iC=1,nC
      Temp1=Temp1+c(iC,iA)*cCat(ic)
      ActCatLoc(iC)=ActCatLoc(iC)+c(iC,iA)*cAn(iA) 
    END DO
    ActAnLoc(iA)=Temp1
    Temp2=Temp2+Temp1*cAn(iA)
  END DO 
  aW=aW-2.0d0*SumZ*Temp2
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ChargeCat(ic)*Temp2+SumZ*ActCatLoc(iC)
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)-ChargeAn(iA)*Temp2+SumZ*ActAnLoc(iA)
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

  HOE=0.0d0
  DHOEDIx=0.0d0
  ActAnLoc=0.0d0
  DO iA=1,nA
    Sum=0.0d0
    DSumDIx=0.0d0
    zi=-ChargeAn(iA)
    DO i1A=iA+1,nA
      zj=-ChargeAn(i1A)
      Sum=Sum+cAn(i1A)*(ThetaX(iA,i1A)+ThetaCalc(zi,zj))
      ActAnLoc(i1A)=ActAnLoc(i1A)+cAn(iA)*(ThetaX(iA,i1A)+ThetaCalc(zi,zj))
      DSumDIx=DSumDIx+cAn(i1A)*ThetaCalcDIx(zi,zj)
    END DO
    HOE=HOE+Sum*cAn(iA)
    ActAnLoc(iA)=ActAnLoc(iA)+Sum
    DHOEDIx=DHOEDIx+DSumDIx*cAn(iA)
  END DO 
  ActCatLoc=0.0d0
  DO iC=1,nC
    Sum=0.0d0
    DSumDIx=0.0d0
    zi=ChargeCat(iC)
    DO i1C=iC+1,nC
      zj=ChargeCat(i1C)
      Sum=Sum+cCat(i1C)*(ThetaM(iC,i1C)+ThetaCalc(zi,zj))
      ActCatLoc(i1C)=ActCatLoc(i1C)+cCat(iC)*(ThetaM(iC,i1C)+ThetaCalc(zi,zj))
      DSumDIx=DSumDIx+cCat(i1C)*ThetaCalcDIx(zi,zj)
    END DO
    HOE=HOE+Sum*cCat(iC)
    ActCatLoc(iC)=ActCatLoc(iC)+Sum
    DHOEDIx=DHOEDIx+DSumDIx*cCat(iC)
  END DO


  aW=aW-2.0d0*HOE-2.0d0*DHOEDIx*IonStr
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)-DHOEDIx*(-ChargeCat(iC)**2)+2.0d0*ActCatLoc(iC)
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)-DHOEDIx*(-ChargeAn(iA)**2)+2.0d0*ActAnLoc(iA)
  END DO
!
!   #(4)# Calculate summation for aac parameters: 
!
  Sum=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iC=1,nC
    Temp2=0.0d0
    ActAnLoc1=0.0d0
    DO i2A=1,nA
      Temp1=0.0d0
      DO i1A=i2A+1,nA
        Temp1=Temp1+cAn(i1A)*PsiX(i1A,i2A,iC)
        ActAnLoc1(i1A)=ActAnLoc1(i1A)+cAn(i2A)*PsiX(i1A,i2A,iC)
      END DO
      Temp2=Temp2+cAn(i2A)*Temp1
      ActAnLoc1(i2A)=ActAnLoc1(i2A)+Temp1
    END DO
    Sum=Sum+cCat(iC)*Temp2
    ActAnLoc=ActAnLoc+ActAnLoc1*cCat(iC)
    ActCatLoc(iC)=Temp2
  END DO
!
  aW=aW-2.0d0*Sum
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)
  END DO

!
!   #(5)# Calculate summation for cca parameters:
!
  Sum=0.0d0
  ActCatLoc=0.0d0
  ActAnLoc=0.0d0
  DO iA=1,nA
    Temp2=0.0d0
    ActCatLoc1=0.0d0
    DO i2C=1,nC
      Temp1=0.0d0
      DO i1C=i2C+1,nC
        Temp1=Temp1+cCat(i1C)*PsiM(i1C,i2C,iA)
        ActCatLoc1(i1C)=ActCatLoc1(i1C)+cCat(i2C)*PsiM(i1C,i2C,iA)
      END DO
      Temp2=Temp2+cCat(i2C)*Temp1
      ActCatLoc1(i2C)=ActCatLoc1(i2C)+Temp1
    END DO
    Sum=Sum+cAn(iA)*Temp2
    ActCatLoc=ActCatLoc+ActCatLoc1*cAn(iA)
    ActAnLoc(iA)=Temp2
  END DO
!
  aW=aW-2.0d0*Sum
  DO iC=1,nC
    ActCat(iC)=ActCat(iC)+ActCatLoc(iC)
    ActCat(iC)=EXP(ActCat(iC))
  END DO
  DO iA=1,nA
    ActAn(iA)=ActAn(iA)+ActAnLoc(iA)
    ActAn(iA)=EXP(ActAn(iA))
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

  REAL(RealKind) ::  Pi,b,c1,eps,pp

  REAL(RealKind), PARAMETER :: u1=3.4279d2, u2=-5.0866d-3, u3=9.46900d-7  &
                              ,u4=-2.0525d0, u5=3.1159d3, u6=-1.8289d+2, u7=-8.0325d+3  &
                              ,u8=4.21142d6, u9=2.1417d0
  REAL(RealKind), PARAMETER :: t0 = 273.15d0
  REAL(RealKind), PARAMETER :: densw = 1.d0,            & ! water density [g/cm**3]
                               avo   = 6.02257d23,    & ! avoActAndro's number [molec./mol]
                               bol   = 1.3804d-16,    & ! bolzmann's constant [erg/K]
                               e     = 4.80298d-10      ! electronic charge [esu]
  Pi=4.d0*ATAN(1.d0)
  b=u7+u8/Temp+u9*Temp
  c1=u4+u5/(u6+Temp)
  eps=u1*EXP(u2*Temp+u3*Temp*Temp)
  pp=p*1.-3   !mbar -> bar
  eps=eps+c1*LOG((b+pp)/(b+1.d3))

  APhiF=1.d0/3.d0*SQRT(2.0d0*pi*avo*densw*1.d-3)*(e*e/eps/bol/Temp)**(1.5d0)
END FUNCTION APhiF

FUNCTION IonicStrength( )

  REAL(RealKind) :: IonicStrength

! calculation of the total ionic strenght

  INTEGER :: ia, ic

  IonicStrength=0.0d0
  DO ia=1,na
    IonicStrength=IonicStrength+cAn(ia)*ChargeAn(ia)*ChargeAn(ia)
  END DO
  DO ic=1,nc
    IonicStrength=IonicStrength+cCat(ic)*ChargeCat(ic)*ChargeCat(ic)
  END DO
  IonicStrength=0.5d0*IonicStrength

END FUNCTION IonicStrength

END SUBROUTINE ActivityCoeffCompute

SUBROUTINE InitActivityPitzer1

  INTEGER :: InputUnit=10
  INTEGER :: ios
  INTEGER :: i,iPos,NumName,PosColon,LenLine,PosName
  INTEGER :: iA,iC 
  INTEGER :: PosAn(2),PosCat(2)
  CHARACTER(20) :: NameLoc(3)
  CHARACTER(1) :: Char
  CHARACTER(100) :: Line

! Set NameAn and NameCat
  nA=0
  nC=0
  nN=0
  DO i=3,nAqua
    IF (INDEX(SpeciesNameAqua(i),'p')>0) THEN
      nC=nC+1
    ELSE IF (INDEX(SpeciesNameAqua(i),'m')>0) THEN
      nA=nA+1
    ELSE IF (SpeciesNameAqua(i)(1:1)/='s') THEN
      nN=nN+1
    END IF
  END DO
  ALLOCATE(NameAn(nA))
  ALLOCATE(NameCat(nC))
  ALLOCATE(NameNeut(nN))
  ALLOCATE(ChargeAn(nA))
  ALLOCATE(ChargeCat(nC))
  nA=0
  nC=0
  nN=0
  DO i=3,nAqua
    IF (INDEX(SpeciesNameAqua(i),'p')>0) THEN
      nC=nC+1
      NameCat(nC)=SpeciesNameAqua(i)
      ChargeCat(nC)=Charge(i)
    ELSE IF (INDEX(SpeciesNameAqua(i),'m')>0) THEN
      nA=nA+1
      NameAn(nA)=SpeciesNameAqua(i)
      ChargeAn(nA)=Charge(i)
    ELSE IF (SpeciesNameAqua(i)(1:1)/='s') THEN
      NameNeut(nN)=SpeciesNameAqua(i)
      nN=nN+1
    END IF
  END DO
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
  OPEN (InputUnit,FILE='Pitzer.dat',STATUS='OLD',IOSTAT=ios)
  IF (ios/=0) THEN
    STOP ' File Pitzer.dat not available'
  END IF
  DO
    READ(InputUnit,'(a100)',END=1) Line
    PosColon=INDEX(Line,':')
    IF (PosColon>0) THEN
      Line(PosColon:PosColon)=''
      LenLine=LEN(TRIM(LINE))
      NumName=1
      NameLoc(:)=''
      iPos=1
      DO i=1,LenLine
        Char=Line(i:i)
        IF (Char/=' ') THEN
          NameLoc(NumName)(iPos:iPos)=Char
          iPos=iPos+1
        ELSE
          NumName=NumName+1 
          iPos=1
        END IF
      END DO
      iA=0
      iC=0
      DO i=1,NumName
        IF (INDEX(NameLoc(i),'m')>0) THEN
          iPos=PositionAn(NameLoc(i)) 
          IF (iPos>0) THEN
            iA=iA+1
            PosAn(iA)=iPos
          END IF 
        ELSE
          iPos=PositionCat(NameLoc(i)) 
          IF (iPos>0) THEN
            iC=iC+1
            PosCat(iC)=iPos
          END IF
        END IF
      END DO
    ELSE
      IF (iA+iC==NumName) THEN
        IF (iA==1.AND.iC==1) THEN
          PosName=INDEX(Line,'Beta0')
          IF (PosName>0) THEN
            READ(Line(PosName+5:),*) b0(PosCat(1),PosAn(1))
          END IF
          PosName=INDEX(Line,'Beta1')
          IF (PosName>0) THEN
            READ(Line(PosName+5:),*) b1(PosCat(1),PosAn(1))
          END IF
          PosName=INDEX(Line,'Beta2')
          IF (PosName>0) THEN
            READ(Line(PosName+5:),*) b2(PosCat(1),PosAn(1))
          END IF
          PosName=INDEX(Line,'CPhi')
          IF (PosName>0) THEN
            READ(Line(PosName+4:),*) c(PosCat(1),PosAn(1))
          END IF
        ELSE IF (iA==2.AND.iC==0) THEN
          PosName=INDEX(Line,'Theta')
          IF (PosName>0) THEN
            READ(Line(PosName+5:),*) ThetaX(PosAn(1),PosAn(2))
            ThetaX(PosAn(2),PosAn(1))=ThetaX(PosAn(1),PosAn(2)) 
          END IF
        ELSE IF (iA==0.AND.iC==2) THEN
          PosName=INDEX(Line,'Theta')
          IF (PosName>0) THEN
            READ(Line(PosName+5:),*) ThetaM(PosCat(1),PosCat(2))
            ThetaM(PosCat(2),PosCat(1))=ThetaM(PosCat(1),PosCat(2))
          END IF
        ELSE IF (iA==2.AND.iC==1) THEN
          PosName=INDEX(Line,'Psi')
          IF (PosName>0) THEN
            READ(Line(PosName+3:),*) PsiX(PosAn(1),PosAn(2),PosCat(1))
            PsiX(PosAn(2),PosAn(1),PosCat(1))=PsiX(PosAn(1),PosAn(2),PosCat(1))
          END IF
        ELSE IF (iA==1.AND.iC==2) THEN
          PosName=INDEX(Line,'Psi')
          IF (PosName>0) THEN
            READ(Line(PosName+3:),*) PsiM(PosCat(1),PosCat(2),PosAn(1))
            PsiM(PosCat(2),PosCat(1),PosAn(1))=PsiM(PosCat(1),PosCat(2),PosAn(1))
          END IF
        END IF
      END IF
    END IF
  END DO
1 CONTINUE

!--- convert c's according to their charge
  DO iC=1,nC
    DO iA=1,nA
      c(iC,iA) = 0.5d0*c(iC,iA)/SQRT(FLOAT(ABS(ChargeAn(iA)*ChargeCat(iC))))
    END DO
  END DO

  CALL OutputCoefficients

END SUBROUTINE InitActivityPitzer1

FUNCTION PositionAn(Name)
 
  INTEGER :: PositionAn
  CHARACTER(*) :: Name

  INTEGER :: i

  PositionAn=0
  DO i=1,nA
    IF (Name==NameAn(i))THEN
      PositionAn=i
      EXIT
    END IF
  END DO

END FUNCTION PositionAn

FUNCTION PositionCat(Name)

  INTEGER :: PositionCat
  CHARACTER(*) :: Name

  INTEGER :: i

  PositionCat=0
  DO i=1,nC
    IF (Name==NameCat(i))THEN
      PositionCat=i
      EXIT
    END IF
  END DO

END FUNCTION PositionCat

FUNCTION PositionNeut(Name)

  INTEGER :: PositionNeut
  CHARACTER(*) :: Name

  INTEGER :: i

  PositionNeut=0
  DO i=1,nN
    IF (Name==NameNeut(i))THEN
      PositionNeut=i
      EXIT
    END IF
  END DO

END FUNCTION PositionNeut


SUBROUTINE ComputeActivityPitzer1(cAqua,ActCoeff)

  TYPE(Vec4_T), TARGET :: cAqua(0:),ActCoeff(0:)

  INTEGER :: i,iFrac,ix,iy,iz,is,iA,iC
  REAL(RealKind) :: Pre,TempLoc
  REAL(RealKind) :: MM(nAqua)
  REAL(RealKind) :: cAn(nA),cCat(nC),cNeut(nN)
  REAL(RealKind) :: ActAn(nA),ActCat(nC),ActNeut(nN)


! =========================================================================
! ===  Determine of Pitzer1 coefficients
! =========================================================================
!
  ActCoeff=One
! -------------------------------------------------------------------------
! ---  Determine of Pitzer1 coefficients for all LWC fractions 
!
   DO ix=ix0+1,ix1
     DO iy=iy0+1,iy1
       DO iz=iz0+1,iz1
!        Pre=pres(iCell) * 1013.e0
         Pre=cAqua(prePos)%c(ix,iy,iz,1)*1.d-2
         TempLoc=cAqua(thPos)%c(ix,iy,iz,1) ! OSSI Umrechnen
         DO iFrac=1,nFrac
           IF (cAqua(iNC)%c(ix,iy,iz,iFrac)>Zero) THEN
             IF (cAqua(iWater)%c(ix,iy,iz,iFrac)<Zero) THEN
                WRITE(*,*) 'NeActAntiv',cAqua(iWater)%c(ix,iy,iz,iFrac),iFrac 
                STOP
             END IF
             iA=0
             iC=0
             DO i=3,nAqua
               IF (INDEX(SpeciesNameAqua(i),'p')>0) THEN
                 iC=iC+1
                 cCat(iC)=cAqua(i)%c(ix,iy,iz,iFrac)/cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(i)
               ELSE IF (INDEX(SpeciesNameAqua(i),'m')>0) THEN
                 iA=iA+1
                 cAn(iA)=cAqua(i)%c(ix,iy,iz,iFrac)/cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(i)
               END IF
             END DO
             CALL ActivityCoeffCompute(cAn,cCat,cNeut, &
                                ActAn,ActCat,ActNeut, &
                                TempLoc,Pre)
             DO is=1,nAqua
               MM(is)=cAqua(is)%c(ix,iy,iz,iFrac)/cAqua(iNC)%c(ix,iy,iz,iFrac)
             END DO
!            ActCoeff(3)%c(ix,iy,iz,iFrac)=aw/RaoultFak(MM)
             iA=0
             iC=0
             DO i=3,nAqua
               IF (INDEX(SpeciesNameAqua(i),'p')>0) THEN
                 iC=iC+1
!                ActCoeff(i)%c(ix,iy,iz,iFrac)=ActCat(iC)
               ELSE IF (INDEX(SpeciesNameAqua(i),'m')>0) THEN
                 iA=iA+1
!                ActCoeff(i)%c(ix,iy,iz,iFrac)=ActAn(iA)
               END IF
             END DO
           END IF 
         END DO
       END DO
     END DO
   END DO
END SUBROUTINE ComputeActivityPitzer1

SUBROUTINE OutputActivityPitzer1(cAqua,ActCoeff)

  TYPE(Vec4_T), TARGET :: cAqua(0:),ActCoeff(0:)

  INTEGER :: i,iFrac,ix,iy,iz,is,iA,iC
  REAL(RealKind) :: Pre,TempLoc
  REAL(RealKind) :: MM(nAqua)
  REAL(RealKind) :: cAn(nA),cCat(nC),cNeut(nN)
  REAL(RealKind) :: ActAn(nA),ActCat(nC),ActNeut(nN)

  WRITE(*,*) 'Output Activity 0'
  DO ix=ix0+1,ix1
    DO iy=iy0+1,iy1
      DO iz=iz0+1,iz1
        Pre=cAqua(prePos)%c(ix,iy,iz,1)*1.d-2
        TempLoc=cAqua(thPos)%c(ix,iy,iz,1) ! OSSI Umrechnen
        DO iFrac=1,nFrac
          IF (cAqua(iNC)%c(ix,iy,iz,iFrac)>Zero) THEN
            IF (cAqua(iWater)%c(ix,iy,iz,iFrac)<Zero) THEN
               WRITE(*,*) 'NeActAntiv',cAqua(iWater)%c(ix,iy,iz,iFrac),iFrac
               STOP
            END IF
            iA=0
            iC=0
            DO i=3,nAqua
              IF (INDEX(SpeciesNameAqua(i),'p')>0) THEN
                iC=iC+1
                cCat(iC)=cAqua(i)%c(ix,iy,iz,iFrac)/cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(i)
              ELSE IF (INDEX(SpeciesNameAqua(i),'m')>0) THEN
                iA=iA+1
                cAn(iA)=cAqua(i)%c(ix,iy,iz,iFrac)/cAqua(iWater)%c(ix,iy,iz,iFrac)/MolMass(i)
              END IF
            END DO
            ActCat(:)=0.e0
            ActAn(:)=0.e0
            WRITE(*,*) 'Output Activity 1'
            CALL ActivityCoeffCompute(cAn,cCat,cNeut, &
                                      ActAn,ActCat,ActNeut, &
                                      TempLoc,Pre)
            WRITE(*,*) 'Output Activity 2'
            DO is=1,nAqua
              MM(is)=cAqua(is)%c(ix,iy,iz,iFrac)/cAqua(iNC)%c(ix,iy,iz,iFrac)
            END DO
            ActCoeff(3)%c(ix,iy,iz,iFrac)=aw/RaoultFak(MM)
            iA=0
            iC=0
            DO i=3,nAqua
              IF (INDEX(SpeciesNameAqua(i),'p')>0) THEN
                iC=iC+1
                ActCoeff(i)%c(ix,iy,iz,iFrac)=ActCat(iC)
              ELSE IF (INDEX(SpeciesNameAqua(i),'m')>0) THEN
                iA=iA+1
                ActCoeff(i)%c(ix,iy,iz,iFrac)=ActAn(iA)
              END IF
            END DO
            WRITE(*,*) 'Output Activity 3'
            DO i=3,nAqua
              WRITE(*,*) SpeciesName(i),cAqua(i)%c(ix,iy,iz,iFrac)/MolMass(i)/cAqua(iWater)%c(ix,iy,iz,iFrac), &
                         ActCoeff(i)%c(ix,iy,iz,iFrac)
            END DO
            WRITE(*,*) 'Output Activity 4'

           END IF
         END DO
       END DO
     END DO
   END DO
END SUBROUTINE OutputActivityPitzer1

SUBROUTINE OutputCoefficients

  INTEGER :: i,j

  WRITE(56,*) 'b0'
  DO i=1,SIZE(b0,1)
    WRITE(56,'(10d17.5)') (b0(i,j),j=1,SIZE(b0,2))
  END DO
  WRITE(56,*) 'b1'
  DO i=1,SIZE(b1,1)
    WRITE(56,'(10d17.5)') (b1(i,j),j=1,SIZE(b1,2))
  END DO
  WRITE(56,*) 'cPhi'
  DO i=1,SIZE(c,1)
    WRITE(56,'(10d17.5)') (c(i,j),j=1,SIZE(c,2))
  END DO
  WRITE(56,*) 'ThetaM'
  DO i=1,SIZE(ThetaM,1)
    WRITE(56,'(10d17.5)') (ThetaM(i,j),j=1,SIZE(ThetaM,2))
  END DO
  WRITE(56,*) 'ThetaX'
  DO i=1,SIZE(ThetaX,1)
    WRITE(56,'(10d17.5)') (ThetaX(i,j),j=1,SIZE(ThetaX,2))
  END DO

END SUBROUTINE OutputCoefficients



END MODULE ActivityPitzer1_Mod

